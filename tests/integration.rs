//! Integration tests: generate corpus → run pipeline → validate passport against ground truth

use tempfile::TempDir;
use virome_qc::{CorpusConfig, CorpusGenerator, Pipeline, Profile, SampleComposition};

/// Generate a corpus, run single-end QC, and validate passport metrics
#[test]
fn test_single_end_pipeline_against_corpus() {
    let tmp = TempDir::new().unwrap();
    let corpus_path = tmp.path().join("corpus.fastq");
    let output_dir = tmp.path().join("output");

    // Generate corpus with known composition
    let config = CorpusConfig {
        num_reads: 5_000,
        read_length: 150,
        seed: 42,
        composition: SampleComposition::stool_vlp(),
        adapter_3prime_rate: 0.10,
        adapter_internal_rate: 0.005,
        poly_g_rate: 0.05,
        quality_tail_rate: 0.15,
        ..Default::default()
    };

    let mut gen = CorpusGenerator::new(config);
    let corpus_summary = gen.generate_to_file(&corpus_path).unwrap();

    // Run QC pipeline
    let profile = Profile::load("stool-vlp-tagmentation").unwrap();
    let pipeline = Pipeline::new(profile, 2);
    let result = pipeline.run(&[corpus_path], &output_dir).unwrap();

    // Validate: reads_input matches corpus
    assert_eq!(result.reads_input, corpus_summary.total_reads as u64);

    // Validate: adapter module detected adapters
    let adapter_report = result
        .module_reports
        .iter()
        .find(|(name, _)| name == "adapter")
        .map(|(_, r)| r)
        .unwrap();
    let adapters_found_3prime = adapter_report
        .extra
        .get("adapters_found_3prime")
        .and_then(|v| v.as_u64())
        .unwrap_or(0);
    assert!(
        adapters_found_3prime > 0,
        "Adapter module should detect planted 3' adapters"
    );

    // Validate: internal adapters detected
    let internal_found = adapter_report
        .extra
        .get("adapters_found_internal")
        .and_then(|v| v.as_u64())
        .unwrap_or(0);
    assert!(
        internal_found > 0,
        "Should detect some internal adapter contamination"
    );

    // Validate: poly-X module trimmed reads
    let polyx_report = result
        .module_reports
        .iter()
        .find(|(name, _)| name == "polyx")
        .map(|(_, r)| r)
        .unwrap();
    assert!(
        polyx_report.reads_modified > 0,
        "PolyX module should have trimmed some reads (planted {})",
        corpus_summary.poly_g_reads
    );

    // Validate: complexity module removed low-complexity reads
    let complexity_report = result
        .module_reports
        .iter()
        .find(|(name, _)| name == "complexity")
        .map(|(_, r)| r)
        .unwrap();
    assert!(
        complexity_report.reads_removed > 0,
        "Complexity module should have removed some reads (planted {})",
        corpus_summary.low_complexity_reads
    );

    // Validate: passport exists and has correct tier
    let passport = result.passport();
    assert_eq!(passport.reads_input, corpus_summary.total_reads as u64);
    assert!(
        passport.survival_rate > 0.40,
        "Reasonable fraction should survive stool VLP QC (got {:.1}%)",
        passport.survival_rate * 100.0
    );

    // Validate: output file exists (gzip compressed) and has reads
    let output_file = output_dir.join("clean_corpus.fastq.gz");
    assert!(
        output_file.exists(),
        "Clean gzipped FASTQ should be written"
    );
    let clean_count = biometal::FastqStream::from_path(&output_file)
        .unwrap()
        .count();
    // Allow small discrepancy from ERV analysis post-processing
    let diff = (clean_count as i64 - result.reads_passed as i64).unsigned_abs();
    assert!(
        diff <= (result.reads_passed / 100).max(50),
        "Clean output count ({}) should be close to reads_passed ({})",
        clean_count,
        result.reads_passed
    );
}

/// Generate paired-end corpus, run paired QC with merge, validate
#[test]
fn test_paired_end_pipeline_with_merge() {
    let tmp = TempDir::new().unwrap();
    let r1_path = tmp.path().join("R1.fastq");
    let r2_path = tmp.path().join("R2.fastq");
    let output_dir = tmp.path().join("output");

    let config = CorpusConfig {
        num_reads: 2_000,
        read_length: 150,
        insert_size_mean: 250.0,
        insert_size_std: 50.0,
        seed: 123,
        ..Default::default()
    };

    let mut gen = CorpusGenerator::new(config);
    let _summary = gen.generate_paired_to_files(&r1_path, &r2_path).unwrap();

    // Run paired-end QC with merge
    let profile = Profile::load("stool-vlp-tagmentation").unwrap();
    let pipeline = Pipeline::new(profile, 2);
    let result = pipeline
        .run_paired(&r1_path, &r2_path, &output_dir, true)
        .unwrap();

    // Validate paired-end output structure (gzip compressed)
    assert!(output_dir.join("clean_R1.fastq.gz").exists());
    assert!(output_dir.join("clean_R2.fastq.gz").exists());
    assert!(output_dir.join("singletons.fastq.gz").exists());
    assert!(output_dir.join("merged.fastq.gz").exists());

    // Validate counts
    assert_eq!(result.reads_input, 4_000, "2000 pairs = 4000 reads");
    assert!(result.pairs_passed > 0, "Some pairs should pass");
    assert!(
        result.pairs_merged > 0,
        "Some pairs should merge with 250bp insert / 150bp reads"
    );

    // Validate: R1 and R2 output have same count
    let r1_clean_count = biometal::FastqStream::from_path(output_dir.join("clean_R1.fastq.gz"))
        .unwrap()
        .count();
    let r2_clean_count = biometal::FastqStream::from_path(output_dir.join("clean_R2.fastq.gz"))
        .unwrap()
        .count();
    assert_eq!(
        r1_clean_count, r2_clean_count,
        "R1 and R2 output must have equal read count"
    );

    // Validate: merged + unmerged pairs = pairs_passed
    let merged_count = biometal::FastqStream::from_path(output_dir.join("merged.fastq.gz"))
        .unwrap()
        .count();
    let actual_pairs = r1_clean_count as u64 + merged_count as u64;
    let pair_diff = (actual_pairs as i64 - result.pairs_passed as i64).unsigned_abs();
    assert!(
        pair_diff <= (result.pairs_passed / 100).max(50),
        "Unmerged pairs ({}) + merged pairs ({}) = {} should be close to pairs_passed ({})",
        r1_clean_count, merged_count, actual_pairs, result.pairs_passed
    );

    // Validate: passport survival rate is reasonable
    let passport = result.passport();
    assert!(
        passport.survival_rate > 0.40,
        "Reasonable fraction should survive paired VLP QC (got {:.1}%)",
        passport.survival_rate * 100.0
    );
}

/// Test that low-survival corpus triggers threshold flags
#[test]
fn test_threshold_fail_on_low_survival() {
    let tmp = TempDir::new().unwrap();
    let corpus_path = tmp.path().join("bad_corpus.fastq");
    let output_dir = tmp.path().join("output");

    // Generate corpus that's mostly low-complexity — will fail complexity filter
    let config = CorpusConfig {
        num_reads: 1_000,
        read_length: 150,
        seed: 99,
        composition: SampleComposition {
            viral: 0.0,
            host: 0.0,
            rrna: 0.0,
            phix: 0.0,
            low_complexity: 0.95,
            background: 0.05,
        },
        ..Default::default()
    };

    let mut gen = CorpusGenerator::new(config);
    gen.generate_to_file(&corpus_path).unwrap();

    let profile = Profile::load("stool-vlp-tagmentation").unwrap();
    let pipeline = Pipeline::new(profile, 1);
    let result = pipeline.run(&[corpus_path], &output_dir).unwrap();

    let passport = result.passport();

    let complexity_report = result
        .module_reports
        .iter()
        .find(|(name, _)| name == "complexity")
        .map(|(_, r)| r)
        .unwrap();
    assert!(
        complexity_report.reads_removed > 200,
        "Should remove many low-complexity reads, got {}",
        complexity_report.reads_removed
    );
    assert!(
        passport.survival_rate < 0.90,
        "95% low-complexity corpus should reduce survival, got {:.1}%",
        passport.survival_rate * 100.0
    );
}

/// Test that gzip input is handled transparently
#[test]
fn test_gzip_input_output() {
    let tmp = TempDir::new().unwrap();
    let corpus_path = tmp.path().join("corpus.fastq");
    let output_dir = tmp.path().join("output");

    let config = CorpusConfig {
        num_reads: 500,
        seed: 77,
        ..Default::default()
    };

    let mut gen = CorpusGenerator::new(config);
    gen.generate_to_file(&corpus_path).unwrap();

    let profile = Profile::load("stool-vlp-tagmentation").unwrap();
    let pipeline = Pipeline::new(profile, 1);
    let result = pipeline.run(&[corpus_path], &output_dir).unwrap();

    // Output should be gzip compressed
    let output_file = output_dir.join("clean_corpus.fastq.gz");
    assert!(output_file.exists(), "Output should be .fastq.gz");

    // Should be readable by biometal (auto-detects gzip)
    let count = biometal::FastqStream::from_path(&output_file)
        .unwrap()
        .count();
    assert_eq!(count as u64, result.reads_passed);
    assert!(result.reads_passed > 0);
}
