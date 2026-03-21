//! Quality trimming and filtering module
//!
//! Uses biometal quality primitives for:
//! - Sliding window quality trimming (3' end)
//! - Leading low-quality base trimming (5' end)
//! - Mean quality filter (whole-read)
//! - Minimum length filter (post-trim)
//!
//! Supports bin-aware mode for NovaSeq/NextSeq quantized quality scores.
//! When quality_binned=true, sliding window uses fraction of low-bin bases
//! instead of mean quality, which is more robust to step-function Q-scores.

use crate::config::QualityConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use biometal::operations::{mean_quality, trim_quality_start};
use biometal::FastqRecord;
use std::sync::atomic::{AtomicU64, Ordering};

/// Quality trimming and filtering module
pub struct QualityTrimmer {
    window_size: usize,
    min_quality: u8,
    min_mean_quality: f64,
    min_length: usize,
    quality_binned: bool,
    stats: AtomicStats,
    failed_quality: AtomicU64,
    failed_length: AtomicU64,
}

impl QualityTrimmer {
    pub fn new(config: &QualityConfig) -> Self {
        Self {
            window_size: config.window_size,
            min_quality: config.min_quality,
            min_mean_quality: config.min_mean_quality,
            min_length: config.min_length,
            quality_binned: config.quality_binned,
            stats: AtomicStats::new(),
            failed_quality: AtomicU64::new(0),
            failed_length: AtomicU64::new(0),
        }
    }

    /// Bin-aware sliding window trim from 3' end
    ///
    /// Instead of checking mean quality (which behaves poorly with 3-4 discrete
    /// quality bins), counts the fraction of bases in the window that fall below
    /// the quality threshold. Trims when >50% of bases in the window are low-quality.
    ///
    /// This handles the NovaSeq/NextSeq step-function quality profile
    /// where Q37 is good, Q12 is marginal, and Q2 is bad. The trim point
    /// is where the bad bases start to dominate the window.
    fn trim_window_binaware(&self, record: &FastqRecord) -> FastqRecord {
        let quality = &record.quality;
        let seq_len = quality.len();

        if seq_len < self.window_size || self.window_size == 0 {
            return record.clone();
        }

        let phred_threshold = self.min_quality + 33; // Phred+33 encoding
        let max_low_fraction = 0.5; // Trim when >50% of window is below threshold

        // Scan from 3' end backward
        let mut trim_end = seq_len;

        for i in (self.window_size..=seq_len).rev() {
            let window_start = i - self.window_size;
            let window = &quality[window_start..i];

            let low_count = window.iter().filter(|&&q| q < phred_threshold).count();
            let low_fraction = low_count as f64 / self.window_size as f64;

            if low_fraction <= max_low_fraction {
                trim_end = i;
                break;
            }

            // If we've scanned all the way to the start and every window is bad
            if i == self.window_size {
                trim_end = 0;
            }
        }

        if trim_end == seq_len {
            return record.clone();
        }

        FastqRecord::new(
            record.id.clone(),
            record.sequence[..trim_end].to_vec(),
            record.quality[..trim_end].to_vec(),
        )
    }

    /// Standard mean-quality sliding window trim from 3' end (for non-binned data)
    fn trim_window_standard(&self, record: &FastqRecord) -> FastqRecord {
        match biometal::operations::trim_quality_window(record, self.min_quality, self.window_size)
        {
            Ok(trimmed) => trimmed,
            Err(_) => record.clone(),
        }
    }
}

impl QcModule for QualityTrimmer {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();
        let original_len = record.len();

        // Step 1: Trim low-quality 5' bases
        if let Ok(trimmed) = trim_quality_start(&record.record, self.min_quality) {
            let bases_trimmed = original_len - trimmed.sequence.len();
            if bases_trimmed > 0 {
                record.record = trimmed;
                record.metrics.bases_trimmed_5prime += bases_trimmed;
            }
        }

        // Step 2: Sliding window trim from 3' end
        // Use bin-aware mode for quantized quality scores (NovaSeq/NextSeq)
        if !record.is_empty() {
            let trimmed = if self.quality_binned {
                self.trim_window_binaware(&record.record)
            } else {
                self.trim_window_standard(&record.record)
            };

            let bases_trimmed = record.record.sequence.len() - trimmed.sequence.len();
            if bases_trimmed > 0 {
                record.record = trimmed;
                record.metrics.bases_trimmed_3prime += bases_trimmed;
            }
        }

        // Step 3: Mean quality filter
        if !record.is_empty() {
            let mean_q = mean_quality(&record.record.quality);
            record.metrics.mean_quality = Some(mean_q);

            if mean_q < self.min_mean_quality {
                record.fail("low_mean_quality");
                self.failed_quality.fetch_add(1, Ordering::Relaxed);
                self.stats.record_removed();
                return;
            }
        }

        // Step 4: Length filter
        if record.len() < self.min_length {
            record.fail("too_short");
            self.failed_length.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
            return;
        }

        // Track modifications
        let total_trimmed = original_len - record.len();
        if total_trimmed > 0 {
            self.stats.record_modified(total_trimmed);
        }
    }

    fn report(&self) -> ModuleReport {
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "failed_quality": self.failed_quality.load(Ordering::Relaxed),
                "failed_length": self.failed_length.load(Ordering::Relaxed),
                "quality_binned": self.quality_binned,
            }),
        )
    }

    fn name(&self) -> &str {
        "quality"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use biometal::FastqRecord;

    fn make_config() -> QualityConfig {
        QualityConfig {
            enabled: true,
            window_size: 4,
            min_quality: 20,
            min_mean_quality: 20.0,
            min_length: 20,
            quality_binned: false,
            max_n_fraction: 0.10,
        }
    }

    fn make_binned_config() -> QualityConfig {
        QualityConfig {
            quality_binned: true,
            ..make_config()
        }
    }

    fn make_record(seq: &[u8], qual: &[u8]) -> AnnotatedRecord {
        AnnotatedRecord::new(FastqRecord::new("test".into(), seq.to_vec(), qual.to_vec()))
    }

    #[test]
    fn test_quality_trim_3prime() {
        let trimmer = QualityTrimmer::new(&make_config());

        // Good quality followed by bad quality
        // Q40 = 'I' (73-33=40), Q2 = '#' (35-33=2)
        let seq = b"ATGCATGCATGCATGCATGCATGCATGC";
        let mut qual = vec![b'I'; 20]; // Q40
        qual.extend_from_slice(&[b'#'; 8]); // Q2

        let mut record = make_record(seq, &qual);
        trimmer.process(&mut record);

        assert!(record.len() < 28);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_quality_fail_too_short() {
        let trimmer = QualityTrimmer::new(&make_config());

        let seq = b"ATGCATGC";
        let qual = vec![b'I'; 8];
        let mut record = make_record(seq, &qual);
        trimmer.process(&mut record);

        assert!(record.is_failed());
    }

    #[test]
    fn test_quality_fail_low_mean() {
        let trimmer = QualityTrimmer::new(&make_config());

        let seq = vec![b'A'; 30];
        let qual = vec![b'#'; 30]; // Q2
        let mut record = make_record(&seq, &qual);
        trimmer.process(&mut record);

        assert!(record.is_failed());
    }

    #[test]
    fn test_quality_pass_good_read() {
        let trimmer = QualityTrimmer::new(&make_config());

        let seq = vec![b'A'; 100];
        let qual = vec![b'I'; 100]; // Q40
        let mut record = make_record(&seq, &qual);
        trimmer.process(&mut record);

        assert!(!record.is_failed());
        assert_eq!(record.len(), 100);
    }

    // --- Bin-aware trimming tests ---

    #[test]
    fn test_binaware_trim_novaseq_pattern() {
        // Simulate NovaSeq 3-bin pattern: Q37 then abrupt drop to Q2
        // Q37 = 'F' (70-33=37), Q2 = '#' (35-33=2)
        let trimmer = QualityTrimmer::new(&make_binned_config());

        let seq = vec![b'A'; 100];
        let mut qual = vec![b'F'; 80]; // Q37
        qual.extend_from_slice(&[b'#'; 20]); // Q2

        let mut record = make_record(&seq, &qual);
        trimmer.process(&mut record);

        // Should trim the Q2 tail
        assert!(
            record.len() <= 84,
            "Should trim Q2 tail, got len={}",
            record.len()
        );
        assert!(
            record.len() >= 76,
            "Should not over-trim, got len={}",
            record.len()
        );
        assert!(!record.is_failed());
    }

    #[test]
    fn test_binaware_mixed_bins() {
        // Simulate NovaSeq mixed quality: Q37, Q23, Q12 interleaved then Q2 tail
        // Q37='F', Q23='8'(56-33=23), Q12='/'(47-33=14 close enough), Q2='#'
        let trimmer = QualityTrimmer::new(&make_binned_config());

        let seq = vec![b'A'; 60];
        let mut qual = Vec::new();
        // Good region: mix of Q37 and Q23 (both above Q20 threshold)
        for i in 0..40 {
            qual.push(if i % 3 == 0 { b'8' } else { b'F' }); // Q23 and Q37
        }
        // Bad region: Q2
        qual.extend_from_slice(&[b'#'; 20]);

        let mut record = make_record(&seq, &qual);
        trimmer.process(&mut record);

        // Should trim the Q2 section but keep the mixed good region
        assert!(
            record.len() >= 36,
            "Should keep good region, got len={}",
            record.len()
        );
        assert!(
            record.len() <= 44,
            "Should trim bad tail, got len={}",
            record.len()
        );
    }

    #[test]
    fn test_binaware_all_high_quality() {
        let trimmer = QualityTrimmer::new(&make_binned_config());

        let seq = vec![b'A'; 100];
        let qual = vec![b'F'; 100]; // All Q37
        let mut record = make_record(&seq, &qual);
        trimmer.process(&mut record);

        assert!(!record.is_failed());
        assert_eq!(record.len(), 100); // no trimming
    }

    #[test]
    fn test_binaware_all_low_quality() {
        let trimmer = QualityTrimmer::new(&make_binned_config());

        let seq = vec![b'A'; 30];
        let qual = vec![b'#'; 30]; // All Q2
        let mut record = make_record(&seq, &qual);
        trimmer.process(&mut record);

        // Should trim everything then fail on mean quality or length
        assert!(record.is_failed());
    }
}
