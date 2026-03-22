//! Pipeline executor with Rayon-parallelized batch processing
//!
//! Supports both single-end and paired-end modes:
//! - Single-end: processes one FASTQ file, outputs one clean FASTQ
//! - Paired-end: processes R1/R2 together, outputs R1, R2, singletons, and optionally merged

use crate::config::{ProfileConfig, Thresholds};
use crate::modules::{
    AdapterTrimmer, AnalyticsSnapshot, ComplexityFilter, ContaminantScreener, HostFilter,
    LengthFilter, MergeConfig, MergeResult, ModuleReport, NFilter, PolyXTrimmer, QcModule,
    QualityTrimmer, ReadAnalytics, ReadMerger,
};
use crate::pipeline::record::AnnotatedRecord;
use crate::report::Passport;
use anyhow::Result;
use biometal::io::DataSource;
use biometal::{FastqStream, FastqWriter, PairedFastqStream};
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// Result of running the QC pipeline on a sample
pub struct PipelineResult {
    pub reads_input: u64,
    pub reads_passed: u64,
    pub reads_failed: u64,
    /// Pairs where both mates passed QC
    pub pairs_passed: u64,
    /// Reads where only one mate passed (singleton)
    pub singletons: u64,
    /// Pairs that were successfully merged
    pub pairs_merged: u64,
    pub module_reports: Vec<(String, ModuleReport)>,
    pub profile_name: String,
    pub thresholds: Thresholds,
    /// Comprehensive read analytics (per-position, distributions, duplication)
    pub qa_stats: Option<AnalyticsSnapshot>,
    /// Input file provenance
    pub provenance: Option<crate::report::passport::Provenance>,
}

impl PipelineResult {
    /// Overall survival rate (individual reads)
    pub fn survival_rate(&self) -> f64 {
        if self.reads_input == 0 {
            return 0.0;
        }
        self.reads_passed as f64 / self.reads_input as f64
    }

    /// Generate a QA passport from the pipeline results
    pub fn passport(&self) -> Passport {
        Passport::from_result(self)
    }
}

/// Context for paired-end block processing — avoids passing many arguments
struct PairedContext<'a> {
    modules: &'a [Box<dyn QcModule + Send + Sync>],
    merger: Option<&'a ReadMerger>,
    writer_r1: &'a std::sync::Mutex<FastqWriter>,
    writer_r2: &'a std::sync::Mutex<FastqWriter>,
    writer_singletons: &'a std::sync::Mutex<FastqWriter>,
    writer_ambiguous: &'a std::sync::Mutex<FastqWriter>,
    writer_merged: Option<&'a std::sync::Mutex<FastqWriter>>,
    reads_input: &'a AtomicU64,
    reads_passed: &'a AtomicU64,
    reads_failed: &'a AtomicU64,
    pairs_passed: &'a AtomicU64,
    singletons_count: &'a AtomicU64,
    pairs_merged: &'a AtomicU64,
    qa: &'a ReadAnalytics,
}

/// The main QC pipeline
pub struct Pipeline {
    config: ProfileConfig,
    threads: usize,
}

impl Pipeline {
    /// Create a new pipeline with the given profile and thread count
    pub fn new(config: ProfileConfig, threads: usize) -> Self {
        Self { config, threads }
    }

    /// Run single-end pipeline on input FASTQ files
    pub fn run(&self, input_files: &[PathBuf], output_dir: &Path) -> Result<PipelineResult> {
        std::fs::create_dir_all(output_dir)?;
        self.init_thread_pool();

        let modules = self.build_modules()?;
        let qa = ReadAnalytics::new();

        let mut total_input = 0u64;
        let mut total_passed = 0u64;
        let mut total_failed = 0u64;

        for input_path in input_files {
            let stem = input_path.file_stem().unwrap_or_default().to_string_lossy();
            // Strip .fastq/.fq from stem if double-extension (e.g., sample.fastq.gz → sample)
            let stem = stem
                .strip_suffix(".fastq")
                .or_else(|| stem.strip_suffix(".fq"))
                .unwrap_or(&stem);
            let output_path = output_dir.join(format!("clean_{stem}.fastq.gz"));

            let (input, passed, failed) =
                self.process_single_end(input_path, &output_path, &modules, &qa)?;

            total_input += input;
            total_passed += passed;
            total_failed += failed;
        }

        let module_reports = modules
            .iter()
            .map(|m| (m.name().to_string(), m.report()))
            .collect();

        let provenance = {
            use crate::report::passport::{InputFileInfo, Provenance};
            let timestamp = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs().to_string())
                .unwrap_or_default();
            let files = input_files
                .iter()
                .map(|p| InputFileInfo {
                    path: p.display().to_string(),
                    size_bytes: std::fs::metadata(p).map(|m| m.len()).unwrap_or(0),
                })
                .collect();
            Provenance {
                timestamp,
                input_files: files,
            }
        };

        Ok(PipelineResult {
            reads_input: total_input,
            reads_passed: total_passed,
            reads_failed: total_failed,
            pairs_passed: 0,
            singletons: 0,
            pairs_merged: 0,
            module_reports,
            profile_name: self.config.name.clone(),
            thresholds: self.config.thresholds.clone(),
            qa_stats: Some(qa.snapshot()),
            provenance: Some(provenance),
        })
    }

    /// Run paired-end pipeline on R1/R2 FASTQ files
    ///
    /// Outputs:
    /// - `{output_dir}/clean_R1.fastq` — R1 reads where both mates passed QC
    /// - `{output_dir}/clean_R2.fastq` — R2 reads where both mates passed QC
    /// - `{output_dir}/singletons.fastq` — reads where only one mate passed
    /// - `{output_dir}/merged.fastq` — merged overlapping pairs (if merge enabled)
    pub fn run_paired(
        &self,
        r1_path: &Path,
        r2_path: &Path,
        output_dir: &Path,
        merge: bool,
    ) -> Result<PipelineResult> {
        std::fs::create_dir_all(output_dir)?;
        self.init_thread_pool();

        let modules = self.build_modules()?;
        let merger = if merge {
            Some(ReadMerger::new(MergeConfig::default()))
        } else {
            None
        };

        // Open paired input stream
        let paired_stream = PairedFastqStream::from_paths(r1_path, r2_path)?;

        // Open output writers (gzip compressed)
        let writer_r1 =
            std::sync::Mutex::new(FastqWriter::create(output_dir.join("clean_R1.fastq.gz"))?);
        let writer_r2 =
            std::sync::Mutex::new(FastqWriter::create(output_dir.join("clean_R2.fastq.gz"))?);
        let writer_singletons =
            std::sync::Mutex::new(FastqWriter::create(output_dir.join("singletons.fastq.gz"))?);
        let writer_ambiguous =
            std::sync::Mutex::new(FastqWriter::create(output_dir.join("ambiguous.fastq.gz"))?);
        let writer_merged = if merge {
            Some(std::sync::Mutex::new(FastqWriter::create(
                output_dir.join("merged.fastq.gz"),
            )?))
        } else {
            None
        };

        // Counters and QA stats
        let reads_input = AtomicU64::new(0);
        let reads_passed = AtomicU64::new(0);
        let reads_failed = AtomicU64::new(0);
        let pairs_passed = AtomicU64::new(0);
        let singletons = AtomicU64::new(0);
        let pairs_merged = AtomicU64::new(0);
        let qa = ReadAnalytics::new();

        let ctx = PairedContext {
            modules: &modules,
            merger: merger.as_ref(),
            writer_r1: &writer_r1,
            writer_r2: &writer_r2,
            writer_singletons: &writer_singletons,
            writer_ambiguous: &writer_ambiguous,
            writer_merged: writer_merged.as_ref(),
            reads_input: &reads_input,
            reads_passed: &reads_passed,
            reads_failed: &reads_failed,
            pairs_passed: &pairs_passed,
            singletons_count: &singletons,
            pairs_merged: &pairs_merged,
            qa: &qa,
        };

        // Collect pairs into blocks for parallel processing
        let block_size = 5_000; // pairs per block
        let mut block: Vec<(biometal::FastqRecord, biometal::FastqRecord)> =
            Vec::with_capacity(block_size);
        let mut blocks_processed = 0u64;

        for pair_result in paired_stream {
            let (r1, r2) = pair_result?;
            block.push((r1, r2));

            if block.len() >= block_size {
                self.process_paired_block(&block, &ctx)?;
                block.clear();

                blocks_processed += 1;
                if blocks_processed % 100 == 0 {
                    eprint!("\r  {}", qa.progress_line());
                }
            }
        }

        // Process remaining pairs
        if !block.is_empty() {
            self.process_paired_block(&block, &ctx)?;
        }

        eprint!("\r  {}\r", " ".repeat(80));

        let module_reports = modules
            .iter()
            .map(|m| (m.name().to_string(), m.report()))
            .collect();

        let provenance = {
            use crate::report::passport::{InputFileInfo, Provenance};
            let timestamp = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs().to_string())
                .unwrap_or_default();
            Provenance {
                timestamp,
                input_files: vec![
                    InputFileInfo {
                        path: r1_path.display().to_string(),
                        size_bytes: std::fs::metadata(r1_path).map(|m| m.len()).unwrap_or(0),
                    },
                    InputFileInfo {
                        path: r2_path.display().to_string(),
                        size_bytes: std::fs::metadata(r2_path).map(|m| m.len()).unwrap_or(0),
                    },
                ],
            }
        };

        Ok(PipelineResult {
            reads_input: reads_input.load(Ordering::Relaxed),
            reads_passed: reads_passed.load(Ordering::Relaxed),
            reads_failed: reads_failed.load(Ordering::Relaxed),
            pairs_passed: pairs_passed.load(Ordering::Relaxed),
            singletons: singletons.load(Ordering::Relaxed),
            pairs_merged: pairs_merged.load(Ordering::Relaxed),
            module_reports,
            profile_name: self.config.name.clone(),
            thresholds: self.config.thresholds.clone(),
            qa_stats: Some(qa.snapshot()),
            provenance: Some(provenance),
        })
    }

    fn init_thread_pool(&self) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .ok();
    }

    /// Process a block of paired reads in parallel
    fn process_paired_block(
        &self,
        block: &[(biometal::FastqRecord, biometal::FastqRecord)],
        ctx: &PairedContext<'_>,
    ) -> Result<()> {
        // Record pre-QC stats
        for (r1, r2) in block {
            ctx.qa.record_input(&r1.sequence, &r1.quality);
            ctx.qa.record_input(&r2.sequence, &r2.quality);
        }

        // Process R1 and R2 in parallel — only modules needed (not writers)
        let modules = ctx.modules;
        let results: Vec<(AnnotatedRecord, AnnotatedRecord)> = block
            .par_iter()
            .map(|(r1, r2)| {
                let mut ann_r1 = AnnotatedRecord::new(r1.clone());
                let mut ann_r2 = AnnotatedRecord::new(r2.clone());

                for module in modules {
                    if !ann_r1.is_failed() {
                        module.process(&mut ann_r1);
                    }
                    if !ann_r2.is_failed() {
                        module.process(&mut ann_r2);
                    }
                }

                // Paired-end concordant flagging:
                // If one mate is contaminant or host, fail the other (same fragment)
                fn is_concordant_fail(ann: &AnnotatedRecord) -> bool {
                    matches!(
                        &ann.disposition,
                        crate::pipeline::Disposition::Fail(reason)
                            if reason.starts_with("contaminant_") || reason == "host"
                    )
                }
                if is_concordant_fail(&ann_r1) && !ann_r2.is_failed() {
                    ann_r2.fail("concordant_mate");
                }
                if is_concordant_fail(&ann_r2) && !ann_r1.is_failed() {
                    ann_r1.fail("concordant_mate");
                }

                // Propagate host_ambiguous flag to mate (same fragment)
                if ann_r1.is_flagged() && !ann_r2.is_flagged() && !ann_r2.is_failed() {
                    ann_r2.flag("host_ambiguous_mate");
                }
                if ann_r2.is_flagged() && !ann_r1.is_flagged() && !ann_r1.is_failed() {
                    ann_r1.flag("host_ambiguous_mate");
                }

                (ann_r1, ann_r2)
            })
            .collect();

        // Write results — serialized for I/O
        ctx.reads_input
            .fetch_add(results.len() as u64 * 2, Ordering::Relaxed);

        let mut w_r1 = ctx.writer_r1.lock().unwrap();
        let mut w_r2 = ctx.writer_r2.lock().unwrap();
        let mut w_sing = ctx.writer_singletons.lock().unwrap();

        for (ann_r1, ann_r2) in &results {
            // Record adapter metrics for both mates regardless of pass/fail
            for ann in [ann_r1, ann_r2] {
                if let Some(ref adapter_name) = ann.metrics.adapter_detected {
                    ctx.qa.record_adapter(adapter_name);
                }
                if ann.metrics.internal_adapter {
                    ctx.qa.record_internal_adapter();
                }
            }

            let r1_pass = !ann_r1.is_failed();
            let r2_pass = !ann_r2.is_failed();

            // Route flagged (ambiguous) reads to separate file
            let either_flagged = ann_r1.is_flagged() || ann_r2.is_flagged();

            match (r1_pass, r2_pass) {
                (true, true) if either_flagged => {
                    // Both passed but at least one flagged -- write to ambiguous
                    let mut w_amb = ctx.writer_ambiguous.lock().unwrap();
                    w_amb.write_record(&ann_r1.record)?;
                    w_amb.write_record(&ann_r2.record)?;
                    ctx.reads_passed.fetch_add(2, Ordering::Relaxed);
                    ctx.pairs_passed.fetch_add(1, Ordering::Relaxed);
                }
                (true, true) => {
                    if let Some(m) = ctx.merger {
                        match m.merge_pair(&ann_r1.record, &ann_r2.record) {
                            MergeResult::Merged(merged) => {
                                ctx.qa.record_insert_size(merged.sequence.len());
                                ctx.qa.record_passed(&merged.sequence, &merged.quality, 0);
                                if let Some(w_merged) = ctx.writer_merged {
                                    w_merged.lock().unwrap().write_record(&merged)?;
                                }
                                ctx.pairs_merged.fetch_add(1, Ordering::Relaxed);
                                ctx.reads_passed.fetch_add(2, Ordering::Relaxed);
                                ctx.pairs_passed.fetch_add(1, Ordering::Relaxed);
                            }
                            MergeResult::Unmerged(r1, r2) => {
                                ctx.qa.record_passed(
                                    &r1.sequence,
                                    &r1.quality,
                                    ann_r1.total_trimmed(),
                                );
                                ctx.qa.record_passed(
                                    &r2.sequence,
                                    &r2.quality,
                                    ann_r2.total_trimmed(),
                                );
                                w_r1.write_record(&r1)?;
                                w_r2.write_record(&r2)?;
                                ctx.reads_passed.fetch_add(2, Ordering::Relaxed);
                                ctx.pairs_passed.fetch_add(1, Ordering::Relaxed);
                            }
                        }
                    } else {
                        ctx.qa.record_passed(
                            &ann_r1.record.sequence,
                            &ann_r1.record.quality,
                            ann_r1.total_trimmed(),
                        );
                        ctx.qa.record_passed(
                            &ann_r2.record.sequence,
                            &ann_r2.record.quality,
                            ann_r2.total_trimmed(),
                        );
                        w_r1.write_record(&ann_r1.record)?;
                        w_r2.write_record(&ann_r2.record)?;
                        ctx.reads_passed.fetch_add(2, Ordering::Relaxed);
                        ctx.pairs_passed.fetch_add(1, Ordering::Relaxed);
                    }
                }
                (true, false) => {
                    ctx.qa.record_passed(
                        &ann_r1.record.sequence,
                        &ann_r1.record.quality,
                        ann_r1.total_trimmed(),
                    );
                    ctx.qa.record_failed();
                    w_sing.write_record(&ann_r1.record)?;
                    ctx.reads_passed.fetch_add(1, Ordering::Relaxed);
                    ctx.reads_failed.fetch_add(1, Ordering::Relaxed);
                    ctx.singletons_count.fetch_add(1, Ordering::Relaxed);
                }
                (false, true) => {
                    ctx.qa.record_failed();
                    ctx.qa.record_passed(
                        &ann_r2.record.sequence,
                        &ann_r2.record.quality,
                        ann_r2.total_trimmed(),
                    );
                    w_sing.write_record(&ann_r2.record)?;
                    ctx.reads_passed.fetch_add(1, Ordering::Relaxed);
                    ctx.reads_failed.fetch_add(1, Ordering::Relaxed);
                    ctx.singletons_count.fetch_add(1, Ordering::Relaxed);
                }
                (false, false) => {
                    ctx.qa.record_failed();
                    ctx.qa.record_failed();
                    ctx.reads_failed.fetch_add(2, Ordering::Relaxed);
                }
            }
        }

        Ok(())
    }

    /// Process a single-end FASTQ file
    fn process_single_end(
        &self,
        input: &Path,
        output: &Path,
        modules: &[Box<dyn QcModule + Send + Sync>],
        qa: &ReadAnalytics,
    ) -> Result<(u64, u64, u64)> {
        let source = DataSource::from_path(input.to_string_lossy().as_ref());
        let stream = FastqStream::new(source)?;

        let writer = FastqWriter::create(output)?;
        let writer = std::sync::Mutex::new(writer);

        // Ambiguous reads writer (host_ambiguous flagged reads)
        let ambiguous_path = output.with_file_name(
            output
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .replace("clean_", "ambiguous_"),
        );
        let ambiguous_writer = FastqWriter::create(&ambiguous_path)?;
        let ambiguous_writer = std::sync::Mutex::new(ambiguous_writer);

        let mut reads_input = 0u64;
        let mut reads_passed = 0u64;
        let mut reads_failed = 0u64;
        let mut blocks_processed = 0u64;

        let block_size = 10_000;
        let mut block = Vec::with_capacity(block_size);

        for record in stream {
            let record = record?;
            block.push(record);

            if block.len() >= block_size {
                let (input, passed, failed) =
                    self.process_single_block(&block, modules, &writer, &ambiguous_writer, qa)?;
                reads_input += input;
                reads_passed += passed;
                reads_failed += failed;
                block.clear();

                blocks_processed += 1;
                if blocks_processed % 100 == 0 {
                    eprint!("\r  {}", qa.progress_line());
                }
            }
        }

        if !block.is_empty() {
            let (input, passed, failed) =
                self.process_single_block(&block, modules, &writer, &ambiguous_writer, qa)?;
            reads_input += input;
            reads_passed += passed;
            reads_failed += failed;
        }

        // Clear progress line
        eprint!("\r  {}\r", " ".repeat(80));

        Ok((reads_input, reads_passed, reads_failed))
    }

    /// Process a block of single-end reads in parallel
    fn process_single_block(
        &self,
        block: &[biometal::FastqRecord],
        modules: &[Box<dyn QcModule + Send + Sync>],
        writer: &std::sync::Mutex<FastqWriter>,
        ambiguous_writer: &std::sync::Mutex<FastqWriter>,
        qa: &ReadAnalytics,
    ) -> Result<(u64, u64, u64)> {
        // Record pre-QC stats
        for record in block {
            qa.record_input(&record.sequence, &record.quality);
        }

        let results: Vec<AnnotatedRecord> = block
            .par_iter()
            .map(|record| {
                let mut annotated = AnnotatedRecord::new(record.clone());
                for module in modules {
                    if annotated.is_failed() {
                        break;
                    }
                    module.process(&mut annotated);
                }
                annotated
            })
            .collect();

        let input = results.len() as u64;
        let mut passed = 0u64;
        let mut failed = 0u64;

        let mut writer = writer.lock().unwrap();
        let mut amb_writer = ambiguous_writer.lock().unwrap();
        for annotated in &results {
            if annotated.is_failed() {
                failed += 1;
                qa.record_failed();
            } else if annotated.is_flagged() {
                // Flagged (ambiguous) reads go to separate file but count as passed
                passed += 1;
                qa.record_passed(
                    &annotated.record.sequence,
                    &annotated.record.quality,
                    annotated.total_trimmed(),
                );
                amb_writer.write_record(&annotated.record)?;
            } else {
                passed += 1;
                qa.record_passed(
                    &annotated.record.sequence,
                    &annotated.record.quality,
                    annotated.total_trimmed(),
                );
                if let Some(ref adapter_name) = annotated.metrics.adapter_detected {
                    qa.record_adapter(adapter_name);
                }
                if annotated.metrics.internal_adapter {
                    qa.record_internal_adapter();
                }
                writer.write_record(&annotated.record)?;
            }
        }

        Ok((input, passed, failed))
    }

    /// Build the ordered chain of QC modules based on profile configuration
    ///
    /// Module order is critical and intentional:
    /// 1. Adapter trim — remove non-biological sequence first
    /// 2. Poly-X trim — remove platform artifacts (poly-G has high Q, must precede quality trim)
    /// 3. N-filter — remove reads with excessive ambiguous bases
    /// 4. Quality trim — trim low-quality tails and filter by mean quality
    /// 5. Complexity filter — assess cleaned sequence entropy
    /// 6. Contaminant screening — rRNA, PhiX, vectors (k-mer index)
    /// 7. Host depletion — Super Bloom k-mer containment (if filter available)
    /// 8. Length filter — final catch for reads shortened by cumulative trimming
    ///
    /// Future phases:
    /// 9. Deduplication (Phase 4)
    fn build_modules(&self) -> Result<Vec<Box<dyn QcModule + Send + Sync>>> {
        let mut modules: Vec<Box<dyn QcModule + Send + Sync>> = Vec::new();
        let cfg = &self.config.modules;

        // 1. Adapter trimming — must be first (removes non-biological sequence)
        if cfg.adapter.enabled {
            modules.push(Box::new(AdapterTrimmer::new(&cfg.adapter)));
        }

        // 2. Poly-X trimming — before quality (poly-G has high Q scores on NovaSeq)
        if cfg.polyx.enabled {
            modules.push(Box::new(PolyXTrimmer::new(&cfg.polyx)));
        }

        // 3. N-filter — remove reads with excessive ambiguous bases
        modules.push(Box::new(NFilter::new(cfg.quality.max_n_fraction)));

        // 4. Quality trimming + mean quality filter
        if cfg.quality.enabled {
            modules.push(Box::new(QualityTrimmer::new(&cfg.quality)));
        }

        // 5. Complexity filter — assess cleaned sequence entropy
        if cfg.complexity.enabled {
            modules.push(Box::new(ComplexityFilter::new(&cfg.complexity)));
        }

        // 6. Contaminant screening — rRNA, PhiX, vectors (Phase 2)
        if cfg.contaminant.enabled {
            modules.push(Box::new(ContaminantScreener::new(&cfg.contaminant)));
        }

        // 7. Host depletion — Super Bloom k-mer containment (Phase 3)
        if cfg.host.enabled {
            match HostFilter::from_filter_file(&cfg.host) {
                Ok(filter) => modules.push(Box::new(filter)),
                Err(e) => {
                    log::warn!("Host depletion disabled: {}", e);
                    eprintln!("  Warning: Host filter not available: {}", e);
                }
            }
        }

        // 8. Final length filter — catch reads shortened by cumulative trimming
        if cfg.quality.enabled {
            modules.push(Box::new(LengthFilter::new(cfg.quality.min_length)));
        }

        Ok(modules)
    }
}
