//! Synthetic FASTQ corpus generator
//!
//! Generates realistic virome sequencing data with known ground truth
//! for benchmarking QC tools. Configurable composition mimics different
//! sample types (stool VLP, tissue, low-biomass, metagenomics).

use crate::corpus::labels::{ReadLabel, ReadLabels};
use crate::corpus::sequences::*;
use anyhow::Result;
use biometal::FastqRecord;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Random number generator (simple LCG for reproducibility without external dep)
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        // LCG parameters from Numerical Recipes
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn next_usize(&mut self, max: usize) -> usize {
        (self.next_u64() as usize) % max
    }

    fn next_base(&mut self) -> u8 {
        [b'A', b'C', b'G', b'T'][self.next_usize(4)]
    }

    /// Generate a random DNA sequence of given length
    fn random_sequence(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.next_base()).collect()
    }

    /// Generate quality scores with optional degradation
    fn quality_scores(&mut self, len: usize, mean_q: u8, tail_start: Option<usize>) -> Vec<u8> {
        let mut quals = Vec::with_capacity(len);
        for i in 0..len {
            let base_q = if let Some(ts) = tail_start {
                if i >= ts {
                    // Degrading quality tail
                    let decay = ((i - ts) as f64 * 2.0).min(mean_q as f64);
                    (mean_q as f64 - decay).max(2.0) as u8
                } else {
                    mean_q
                }
            } else {
                mean_q
            };

            // Add some noise (+/- 5)
            let noise = (self.next_usize(11) as i8) - 5;
            let q = (base_q as i8 + noise).clamp(2, 41) as u8;
            quals.push(q + 33); // Phred+33 encoding
        }
        quals
    }

    /// Generate approximately normal-distributed value using Box-Muller
    fn next_normal(&mut self, mean: f64, std_dev: f64) -> f64 {
        let u1 = self.next_f64().max(1e-10); // avoid log(0)
        let u2 = self.next_f64();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        mean + std_dev * z
    }
}

/// Composition of a synthetic sample (fractions must sum to ~1.0)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleComposition {
    /// Fraction of reads from viral genomes (true signal)
    pub viral: f64,
    /// Fraction of reads from host genome
    pub host: f64,
    /// Fraction of reads from rRNA
    pub rrna: f64,
    /// Fraction of reads from PhiX spike-in
    pub phix: f64,
    /// Fraction of low-complexity reads
    pub low_complexity: f64,
    /// Fraction of reads with random/novel sequence (background)
    pub background: f64,
}

impl SampleComposition {
    /// Stool VLP enrichment — mostly viral, some rRNA contamination
    pub fn stool_vlp() -> Self {
        Self {
            viral: 0.60,
            host: 0.05,
            rrna: 0.08,
            phix: 0.01,
            low_complexity: 0.02,
            background: 0.24,
        }
    }

    /// Tissue biopsy — mostly host
    pub fn tissue() -> Self {
        Self {
            viral: 0.02,
            host: 0.90,
            rrna: 0.03,
            phix: 0.01,
            low_complexity: 0.01,
            background: 0.03,
        }
    }

    /// Shotgun metagenomics — mixed
    pub fn metagenomics() -> Self {
        Self {
            viral: 0.15,
            host: 0.30,
            rrna: 0.10,
            phix: 0.01,
            low_complexity: 0.02,
            background: 0.42,
        }
    }

    /// Low biomass (plasma/serum) — very high host, trace viral
    pub fn low_biomass() -> Self {
        Self {
            viral: 0.005,
            host: 0.95,
            rrna: 0.02,
            phix: 0.01,
            low_complexity: 0.005,
            background: 0.01,
        }
    }
}

/// Configuration for corpus generation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorpusConfig {
    /// Number of reads to generate
    pub num_reads: usize,
    /// Read length (before artifact addition)
    pub read_length: usize,
    /// Sample composition
    pub composition: SampleComposition,
    /// Fraction of reads with 3' adapter contamination
    pub adapter_3prime_rate: f64,
    /// Fraction of reads with internal adapter contamination
    pub adapter_internal_rate: f64,
    /// Fraction of reads with random primer bias
    pub random_primer_rate: f64,
    /// Length of random primer to simulate
    pub random_primer_length: usize,
    /// Fraction of reads with quality degradation tail
    pub quality_tail_rate: f64,
    /// Fraction of reads with poly-G tail (NovaSeq artifact)
    pub poly_g_rate: f64,
    /// PCR duplicate rate
    pub pcr_dup_rate: f64,
    /// Mean insert size for paired-end generation
    pub insert_size_mean: f64,
    /// Standard deviation of insert size
    pub insert_size_std: f64,
    /// Random seed for reproducibility
    pub seed: u64,
}

impl Default for CorpusConfig {
    fn default() -> Self {
        Self {
            num_reads: 100_000,
            read_length: 150,
            composition: SampleComposition::stool_vlp(),
            adapter_3prime_rate: 0.08,
            adapter_internal_rate: 0.005,
            random_primer_rate: 0.0,
            random_primer_length: 8,
            quality_tail_rate: 0.15,
            poly_g_rate: 0.03,
            pcr_dup_rate: 0.05,
            insert_size_mean: 250.0,
            insert_size_std: 50.0,
            seed: 42,
        }
    }
}

impl CorpusConfig {
    /// Configuration mimicking stool VLP tagmentation prep
    pub fn stool_vlp() -> Self {
        Self {
            composition: SampleComposition::stool_vlp(),
            adapter_3prime_rate: 0.08,
            adapter_internal_rate: 0.003,
            poly_g_rate: 0.03,
            ..Default::default()
        }
    }

    /// Configuration mimicking tissue RNA-seq
    pub fn tissue() -> Self {
        Self {
            composition: SampleComposition::tissue(),
            adapter_3prime_rate: 0.05,
            adapter_internal_rate: 0.001,
            quality_tail_rate: 0.10,
            ..Default::default()
        }
    }

    /// Configuration mimicking low-biomass WGA
    pub fn low_biomass_wga() -> Self {
        Self {
            composition: SampleComposition::low_biomass(),
            adapter_3prime_rate: 0.12,
            adapter_internal_rate: 0.01,
            random_primer_rate: 0.80, // most reads have primer bias
            random_primer_length: 13,
            pcr_dup_rate: 0.30, // high duplication
            ..Default::default()
        }
    }
}

/// Generates synthetic FASTQ corpus with ground truth labels
pub struct CorpusGenerator {
    config: CorpusConfig,
    rng: Rng,
}

impl CorpusGenerator {
    pub fn new(config: CorpusConfig) -> Self {
        let rng = Rng::new(config.seed);
        Self { config, rng }
    }

    /// Generate the corpus and write to a FASTQ file
    pub fn generate_to_file(&mut self, path: &Path) -> Result<CorpusSummary> {
        let mut writer = biometal::FastqWriter::create(path)?;
        let mut summary = CorpusSummary::default();

        // Pre-generate some reads for PCR duplication
        let dup_source_count =
            (self.config.num_reads as f64 * self.config.pcr_dup_rate * 0.5) as usize;
        let mut dup_sources: Vec<(FastqRecord, ReadLabels)> = Vec::new();

        for i in 0..self.config.num_reads {
            let (record, labels) =
                if !dup_sources.is_empty() && self.rng.next_f64() < self.config.pcr_dup_rate {
                    // Generate a PCR duplicate
                    let source_idx = self.rng.next_usize(dup_sources.len());
                    let (source_record, source_labels) = &dup_sources[source_idx];
                    let mut labels = source_labels.clone();
                    labels.add(ReadLabel::PcrDuplicate(source_record.id.clone()));
                    summary.pcr_duplicates += 1;

                    let record = FastqRecord::new(
                        format!("read_{:06} {}", i, labels.to_comment()),
                        source_record.sequence.clone(),
                        source_record.quality.clone(),
                    );
                    (record, labels)
                } else {
                    // Generate a fresh read
                    let (mut seq, mut labels) = self.generate_base_read(i);

                    // Apply artifacts
                    self.maybe_add_adapter_3prime(&mut seq, &mut labels);
                    self.maybe_add_adapter_internal(&mut seq, &mut labels);
                    self.maybe_add_random_primer(&mut seq, &mut labels);
                    self.maybe_add_poly_g(&mut seq, &mut labels);

                    // Generate quality scores (with possible tail degradation)
                    let tail_start = if self.rng.next_f64() < self.config.quality_tail_rate {
                        let pos = self.config.read_length / 2
                            + self.rng.next_usize(self.config.read_length / 3);
                        labels.add(ReadLabel::QualityTail(pos));
                        summary.quality_tails += 1;
                        Some(pos.min(seq.len()))
                    } else {
                        None
                    };

                    let quality = self.rng.quality_scores(seq.len(), 35, tail_start);

                    // Save as potential dup source
                    if dup_sources.len() < dup_source_count {
                        let record =
                            FastqRecord::new(format!("read_{i:06}"), seq.clone(), quality.clone());
                        dup_sources.push((record, labels.clone()));
                    }

                    // Update summary
                    self.update_summary(&labels, &mut summary);

                    let record = FastqRecord::new(
                        format!("read_{:06} {}", i, labels.to_comment()),
                        seq,
                        quality,
                    );
                    (record, labels)
                };

            let _ = labels; // labels already encoded in record ID
            writer.write_record(&record)?;
        }

        summary.total_reads = self.config.num_reads;

        // Write summary as JSON alongside FASTQ
        let summary_path = path.with_extension("summary.json");
        let summary_json = serde_json::to_string_pretty(&summary)?;
        std::fs::write(summary_path, summary_json)?;

        Ok(summary)
    }

    /// Generate paired-end FASTQ corpus (R1 + R2 files)
    ///
    /// Generates a fragment of `insert_size` length, then derives R1 (forward)
    /// and R2 (reverse complement of fragment end). When insert < 2*read_length,
    /// reads overlap. When insert < read_length, adapter read-through occurs.
    pub fn generate_paired_to_files(
        &mut self,
        r1_path: &Path,
        r2_path: &Path,
    ) -> Result<CorpusSummary> {
        let mut writer_r1 = biometal::FastqWriter::create(r1_path)?;
        let mut writer_r2 = biometal::FastqWriter::create(r2_path)?;
        let mut summary = CorpusSummary::default();

        for i in 0..self.config.num_reads {
            // Generate fragment (insert) with biological content
            let insert_size = self
                .rng
                .next_normal(self.config.insert_size_mean, self.config.insert_size_std)
                .round()
                .max(50.0) as usize;

            // Generate base fragment sequence
            let (fragment, labels) = self.generate_base_fragment(i, insert_size);

            // Derive R1 (forward, first read_length bases of fragment)
            let read_len = self.config.read_length;
            let mut r1_seq: Vec<u8>;
            let mut r2_seq: Vec<u8>;
            let mut r1_labels = labels.clone();
            let mut r2_labels = labels.clone();

            if fragment.len() >= read_len {
                r1_seq = fragment[..read_len].to_vec();
            } else {
                // Insert shorter than read length — adapter read-through
                r1_seq = fragment.clone();
                // Append adapter sequence
                let adapter = &ADAPTERS[0]; // TruSeq R1
                let adapter_len = read_len - fragment.len();
                let adapter_append: Vec<u8> = adapter
                    .sequence
                    .iter()
                    .cycle()
                    .take(adapter_len)
                    .copied()
                    .collect();
                r1_seq.extend_from_slice(&adapter_append);
                r1_labels.add(ReadLabel::Adapter3Prime(
                    adapter.name.to_string(),
                    adapter_len,
                ));
                summary.adapter_3prime_reads += 1;
            }

            // Derive R2 (reverse complement of fragment end)
            let r2_fragment_start = fragment.len().saturating_sub(read_len);
            let r2_fragment = &fragment[r2_fragment_start..];
            let r2_rc = biometal::operations::reverse_complement(r2_fragment);
            if r2_rc.len() >= read_len {
                r2_seq = r2_rc[..read_len].to_vec();
            } else {
                r2_seq = r2_rc;
                // Adapter read-through on R2
                let adapter = &ADAPTERS[1]; // TruSeq R2
                let adapter_len = read_len - r2_seq.len();
                let adapter_append: Vec<u8> = adapter
                    .sequence
                    .iter()
                    .cycle()
                    .take(adapter_len)
                    .copied()
                    .collect();
                r2_seq.extend_from_slice(&adapter_append);
                r2_labels.add(ReadLabel::Adapter3Prime(
                    adapter.name.to_string(),
                    adapter_len,
                ));
            }

            // Apply additional artifacts to R1
            self.maybe_add_adapter_internal(&mut r1_seq, &mut r1_labels);
            self.maybe_add_random_primer(&mut r1_seq, &mut r1_labels);
            self.maybe_add_poly_g(&mut r1_seq, &mut r1_labels);

            // Apply additional artifacts to R2
            self.maybe_add_poly_g(&mut r2_seq, &mut r2_labels);

            // Add insert size to labels
            r1_labels.add(ReadLabel::Source(format!("insert_size={}", insert_size)));

            // Generate quality scores
            let r1_tail = if self.rng.next_f64() < self.config.quality_tail_rate {
                Some(read_len * 2 / 3 + self.rng.next_usize(read_len / 4))
            } else {
                None
            };
            let r2_tail = if self.rng.next_f64() < self.config.quality_tail_rate {
                Some(read_len * 2 / 3 + self.rng.next_usize(read_len / 4))
            } else {
                None
            };

            let r1_qual = self.rng.quality_scores(r1_seq.len(), 35, r1_tail);
            let r2_qual = self.rng.quality_scores(r2_seq.len(), 33, r2_tail); // R2 slightly lower quality

            self.update_summary(&labels, &mut summary);

            let r1_record = FastqRecord::new(
                format!("read_{:06}/1 {}", i, r1_labels.to_comment()),
                r1_seq,
                r1_qual,
            );
            let r2_record = FastqRecord::new(
                format!("read_{:06}/2 {}", i, r2_labels.to_comment()),
                r2_seq,
                r2_qual,
            );

            writer_r1.write_record(&r1_record)?;
            writer_r2.write_record(&r2_record)?;
        }

        summary.total_reads = self.config.num_reads * 2; // total reads across both files

        // Write summary
        let summary_path = r1_path.with_extension("summary.json");
        let summary_json = serde_json::to_string_pretty(&summary)?;
        std::fs::write(summary_path, summary_json)?;

        Ok(summary)
    }

    /// Generate a fragment (insert) of specified size with biological content
    fn generate_base_fragment(
        &mut self,
        read_idx: usize,
        insert_size: usize,
    ) -> (Vec<u8>, ReadLabels) {
        // Temporarily set read_length to insert_size to reuse generate_base_read
        let original_read_length = self.config.read_length;
        self.config.read_length = insert_size;
        let result = self.generate_base_read(read_idx);
        self.config.read_length = original_read_length;
        result
    }

    /// Generate the base biological sequence for a read (before artifacts)
    fn generate_base_read(&mut self, _read_idx: usize) -> (Vec<u8>, ReadLabels) {
        let mut labels = ReadLabels::new();
        let roll = self.rng.next_f64();
        let comp = &self.config.composition;

        let mut cumulative = 0.0;
        let seq;

        // Viral
        cumulative += comp.viral;
        if roll < cumulative {
            let vref = &VIRAL_REFERENCES[self.rng.next_usize(VIRAL_REFERENCES.len())];
            seq = self.sample_from_reference(vref.sequence);
            labels.add(ReadLabel::Source("viral".into()));
            labels.add(ReadLabel::Virus(vref.name.to_string()));
            let score = biometal::operations::complexity_score(&seq);
            labels.add(ReadLabel::Complexity(score));
            return (seq, labels);
        }

        // Host
        cumulative += comp.host;
        if roll < cumulative {
            seq = self.sample_from_reference(HUMAN_FRAGMENT);
            labels.add(ReadLabel::Source("host".into()));
            labels.add(ReadLabel::Host("human".into()));
            let score = biometal::operations::complexity_score(&seq);
            labels.add(ReadLabel::Complexity(score));
            return (seq, labels);
        }

        // rRNA
        cumulative += comp.rrna;
        if roll < cumulative {
            let (subunit, reference) = if self.rng.next_f64() < 0.6 {
                ("16S", RRNA_16S_FRAGMENT)
            } else {
                ("18S", RRNA_18S_FRAGMENT)
            };
            seq = self.sample_from_reference(reference);
            labels.add(ReadLabel::Source("rrna".into()));
            labels.add(ReadLabel::Rrna(subunit.into()));
            let score = biometal::operations::complexity_score(&seq);
            labels.add(ReadLabel::Complexity(score));
            return (seq, labels);
        }

        // PhiX
        cumulative += comp.phix;
        if roll < cumulative {
            seq = self.sample_from_reference(PHIX174_FRAGMENT);
            labels.add(ReadLabel::Source("phix".into()));
            labels.add(ReadLabel::PhiX);
            let score = biometal::operations::complexity_score(&seq);
            labels.add(ReadLabel::Complexity(score));
            return (seq, labels);
        }

        // Low complexity
        cumulative += comp.low_complexity;
        if roll < cumulative {
            let template = LOW_COMPLEXITY_SEQS[self.rng.next_usize(LOW_COMPLEXITY_SEQS.len())];
            // Repeat template to fill read length
            seq = template
                .iter()
                .cycle()
                .take(self.config.read_length)
                .copied()
                .collect();
            labels.add(ReadLabel::Source("low_complexity".into()));
            let score = biometal::operations::complexity_score(&seq);
            labels.add(ReadLabel::Complexity(score));
            return (seq, labels);
        }

        // Background (random sequence)
        seq = self.rng.random_sequence(self.config.read_length);
        labels.add(ReadLabel::Source("background".into()));
        let score = biometal::operations::complexity_score(&seq);
        labels.add(ReadLabel::Complexity(score));
        (seq, labels)
    }

    /// Sample a read-length fragment from a reference sequence
    fn sample_from_reference(&mut self, reference: &[u8]) -> Vec<u8> {
        let read_len = self.config.read_length;
        if reference.len() <= read_len {
            // Reference shorter than read — extend with random sequence
            let mut seq = reference.to_vec();
            while seq.len() < read_len {
                seq.push(self.rng.next_base());
            }
            return seq;
        }

        let start = self.rng.next_usize(reference.len() - read_len);
        let mut seq = reference[start..start + read_len].to_vec();

        // Add ~1% substitution errors
        for base in &mut seq {
            if self.rng.next_f64() < 0.01 {
                *base = self.rng.next_base();
            }
        }
        seq
    }

    /// Maybe append a 3' adapter to the read
    fn maybe_add_adapter_3prime(&mut self, seq: &mut Vec<u8>, labels: &mut ReadLabels) {
        if self.rng.next_f64() >= self.config.adapter_3prime_rate {
            return;
        }

        let adapter = &ADAPTERS[self.rng.next_usize(ADAPTERS.len())];
        // Simulate varying insert sizes: overlap = 10 to adapter.len()
        let overlap = 10
            + self
                .rng
                .next_usize(adapter.sequence.len().saturating_sub(10));
        let overlap = overlap.min(adapter.sequence.len());

        // Truncate read to make room, then append adapter
        let insert_end = seq.len().saturating_sub(overlap);
        seq.truncate(insert_end);
        seq.extend_from_slice(&adapter.sequence[..overlap]);

        labels.add(ReadLabel::Adapter3Prime(adapter.name.to_string(), overlap));
    }

    /// Maybe insert adapter sequence internally (chimeric read)
    fn maybe_add_adapter_internal(&mut self, seq: &mut [u8], labels: &mut ReadLabels) {
        if self.rng.next_f64() >= self.config.adapter_internal_rate {
            return;
        }

        let adapter = &ADAPTERS[self.rng.next_usize(ADAPTERS.len())];
        let adapter_fragment_len = 20.min(adapter.sequence.len());

        if seq.len() < 60 {
            return; // too short to insert internally
        }

        // Insert adapter fragment at a random internal position
        let insert_pos = 20 + self.rng.next_usize(seq.len() - 40);
        let end_pos = (insert_pos + adapter_fragment_len).min(seq.len());

        seq[insert_pos..end_pos].copy_from_slice(&adapter.sequence[..end_pos - insert_pos]);

        labels.add(ReadLabel::AdapterInternal(insert_pos));
    }

    /// Maybe add random primer bias to the 5' end
    fn maybe_add_random_primer(&mut self, seq: &mut [u8], labels: &mut ReadLabels) {
        if self.rng.next_f64() >= self.config.random_primer_rate {
            return;
        }

        let primer_len = self.config.random_primer_length.min(seq.len());
        // Replace first N bases with random sequence (simulating primer-derived composition)
        for base in seq[..primer_len].iter_mut() {
            *base = self.rng.next_base();
        }

        labels.add(ReadLabel::RandomPrimer(primer_len));
    }

    /// Maybe add poly-G tail (NovaSeq artifact)
    fn maybe_add_poly_g(&mut self, seq: &mut [u8], labels: &mut ReadLabels) {
        if self.rng.next_f64() >= self.config.poly_g_rate {
            return;
        }

        let poly_len = 10 + self.rng.next_usize(40);
        let poly_len = poly_len.min(seq.len() / 2);
        let start = seq.len() - poly_len;
        for base in seq[start..].iter_mut() {
            *base = b'G';
        }

        labels.add(ReadLabel::PolyX('G', poly_len));
    }

    fn update_summary(&self, labels: &ReadLabels, summary: &mut CorpusSummary) {
        match labels.source() {
            Some("viral") => summary.viral_reads += 1,
            Some("host") => summary.host_reads += 1,
            Some("rrna") => summary.rrna_reads += 1,
            Some("phix") => summary.phix_reads += 1,
            Some("low_complexity") => summary.low_complexity_reads += 1,
            Some("background") => summary.background_reads += 1,
            _ => {}
        }
        if labels.has_adapter_3prime() {
            summary.adapter_3prime_reads += 1;
        }
        if labels.has_internal_adapter() {
            summary.adapter_internal_reads += 1;
        }
        if labels.has_polyx() {
            summary.poly_g_reads += 1;
        }
    }
}

/// Summary statistics for a generated corpus
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CorpusSummary {
    pub total_reads: usize,
    pub viral_reads: usize,
    pub host_reads: usize,
    pub rrna_reads: usize,
    pub phix_reads: usize,
    pub low_complexity_reads: usize,
    pub background_reads: usize,
    pub adapter_3prime_reads: usize,
    pub adapter_internal_reads: usize,
    pub poly_g_reads: usize,
    pub quality_tails: usize,
    pub pcr_duplicates: usize,
}

impl CorpusSummary {
    pub fn print_report(&self) {
        println!("=== Corpus Summary ===");
        println!("Total reads:          {:>8}", self.total_reads);
        println!("--- Source composition ---");
        println!(
            "  Viral:              {:>8} ({:.1}%)",
            self.viral_reads,
            self.pct(self.viral_reads)
        );
        println!(
            "  Host:               {:>8} ({:.1}%)",
            self.host_reads,
            self.pct(self.host_reads)
        );
        println!(
            "  rRNA:               {:>8} ({:.1}%)",
            self.rrna_reads,
            self.pct(self.rrna_reads)
        );
        println!(
            "  PhiX:               {:>8} ({:.1}%)",
            self.phix_reads,
            self.pct(self.phix_reads)
        );
        println!(
            "  Low complexity:     {:>8} ({:.1}%)",
            self.low_complexity_reads,
            self.pct(self.low_complexity_reads)
        );
        println!(
            "  Background:         {:>8} ({:.1}%)",
            self.background_reads,
            self.pct(self.background_reads)
        );
        println!("--- Artifacts ---");
        println!(
            "  3' adapter:         {:>8} ({:.1}%)",
            self.adapter_3prime_reads,
            self.pct(self.adapter_3prime_reads)
        );
        println!(
            "  Internal adapter:   {:>8} ({:.1}%)",
            self.adapter_internal_reads,
            self.pct(self.adapter_internal_reads)
        );
        println!(
            "  Poly-G tail:        {:>8} ({:.1}%)",
            self.poly_g_reads,
            self.pct(self.poly_g_reads)
        );
        println!(
            "  Quality tail:       {:>8} ({:.1}%)",
            self.quality_tails,
            self.pct(self.quality_tails)
        );
        println!(
            "  PCR duplicates:     {:>8} ({:.1}%)",
            self.pcr_duplicates,
            self.pct(self.pcr_duplicates)
        );
    }

    fn pct(&self, count: usize) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            count as f64 / self.total_reads as f64 * 100.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_generate_small_corpus() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test_corpus.fastq");

        let config = CorpusConfig {
            num_reads: 1000,
            read_length: 150,
            composition: SampleComposition::stool_vlp(),
            seed: 42,
            ..Default::default()
        };

        let mut gen = CorpusGenerator::new(config);
        let summary = gen.generate_to_file(&path).unwrap();

        assert_eq!(summary.total_reads, 1000);
        assert!(summary.viral_reads > 0, "Should have viral reads");
        assert!(summary.host_reads > 0, "Should have host reads");
        assert!(
            summary.adapter_3prime_reads > 0,
            "Should have adapter reads"
        );

        // Verify file exists and is readable
        assert!(path.exists());
        let stream = biometal::FastqStream::from_path(path).unwrap();
        let count = stream.count();
        assert_eq!(count, 1000);
    }

    #[test]
    fn test_deterministic_generation() {
        let tmp = TempDir::new().unwrap();

        let config = CorpusConfig {
            num_reads: 100,
            seed: 12345,
            ..Default::default()
        };

        // Generate twice with same seed
        let path1 = tmp.path().join("corpus1.fastq");
        let path2 = tmp.path().join("corpus2.fastq");

        let mut gen1 = CorpusGenerator::new(config.clone());
        let mut gen2 = CorpusGenerator::new(config);

        let summary1 = gen1.generate_to_file(&path1).unwrap();
        let summary2 = gen2.generate_to_file(&path2).unwrap();

        assert_eq!(summary1.viral_reads, summary2.viral_reads);
        assert_eq!(summary1.host_reads, summary2.host_reads);
        assert_eq!(summary1.adapter_3prime_reads, summary2.adapter_3prime_reads);
    }

    #[test]
    fn test_composition_profiles() {
        let tmp = TempDir::new().unwrap();

        // Test tissue profile has high host fraction
        let config = CorpusConfig {
            num_reads: 10_000,
            composition: SampleComposition::tissue(),
            seed: 42,
            ..Default::default()
        };

        let path = tmp.path().join("tissue.fastq");
        let mut gen = CorpusGenerator::new(config);
        let summary = gen.generate_to_file(&path).unwrap();

        let host_fraction = summary.host_reads as f64 / summary.total_reads as f64;
        assert!(
            host_fraction > 0.80,
            "Tissue should be >80% host, got {:.1}%",
            host_fraction * 100.0
        );
    }

    #[test]
    fn test_labels_in_fastq_headers() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("labeled.fastq");

        let config = CorpusConfig {
            num_reads: 100,
            seed: 42,
            ..Default::default()
        };

        let mut gen = CorpusGenerator::new(config);
        gen.generate_to_file(&path).unwrap();

        // Read back and check labels
        let stream = biometal::FastqStream::from_path(path).unwrap();
        let mut found_labels = false;
        for record in stream {
            let record = record.unwrap();
            if record.id.contains("source=") {
                found_labels = true;
                let comment = record.id.split_once(' ').map(|(_, c)| c).unwrap_or("");
                let labels = ReadLabels::from_comment(comment);
                assert!(labels.source().is_some());
                break;
            }
        }
        assert!(found_labels, "Should find labeled reads in corpus");
    }

    #[test]
    fn test_generate_paired_end_corpus() {
        let tmp = TempDir::new().unwrap();
        let r1_path = tmp.path().join("test_R1.fastq");
        let r2_path = tmp.path().join("test_R2.fastq");

        let config = CorpusConfig {
            num_reads: 1000,
            read_length: 150,
            insert_size_mean: 250.0,
            insert_size_std: 50.0,
            seed: 42,
            ..Default::default()
        };

        let mut gen = CorpusGenerator::new(config);
        let summary = gen.generate_paired_to_files(&r1_path, &r2_path).unwrap();

        // Should have 2000 total reads (1000 pairs)
        assert_eq!(summary.total_reads, 2000);
        assert!(summary.viral_reads > 0);

        // Both files should exist and have the same number of records
        let r1_count = biometal::FastqStream::from_path(&r1_path).unwrap().count();
        let r2_count = biometal::FastqStream::from_path(&r2_path).unwrap().count();
        assert_eq!(r1_count, 1000);
        assert_eq!(r2_count, 1000);

        // Check that R1 and R2 have matching read names
        let r1_stream = biometal::FastqStream::from_path(&r1_path).unwrap();
        let r2_stream = biometal::FastqStream::from_path(&r2_path).unwrap();
        for (r1, r2) in r1_stream.zip(r2_stream) {
            let r1 = r1.unwrap();
            let r2 = r2.unwrap();
            let r1_name = r1.id.split('/').next().unwrap();
            let r2_name = r2.id.split('/').next().unwrap();
            assert_eq!(r1_name, r2_name, "R1 and R2 names should match");
            break; // just check the first pair
        }
    }

    #[test]
    fn test_paired_end_short_insert_has_adapters() {
        let tmp = TempDir::new().unwrap();
        let r1_path = tmp.path().join("short_insert_R1.fastq");
        let r2_path = tmp.path().join("short_insert_R2.fastq");

        // Very short inserts → guaranteed adapter read-through
        let config = CorpusConfig {
            num_reads: 500,
            read_length: 150,
            insert_size_mean: 100.0, // shorter than read length!
            insert_size_std: 10.0,
            adapter_3prime_rate: 0.0, // disable random adapters
            seed: 42,
            ..Default::default()
        };

        let mut gen = CorpusGenerator::new(config);
        let summary = gen.generate_paired_to_files(&r1_path, &r2_path).unwrap();

        // Most reads should have adapter contamination from read-through
        assert!(
            summary.adapter_3prime_reads > 200,
            "Short inserts should cause adapter read-through, got {} adapters",
            summary.adapter_3prime_reads
        );
    }
}
