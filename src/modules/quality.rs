//! Quality trimming and filtering module
//!
//! Uses biometal quality primitives for:
//! - Sliding window quality trimming (3' end)
//! - Leading low-quality base trimming (5' end)
//! - Mean quality filter (whole-read)
//! - Minimum length filter (post-trim)

use crate::config::QualityConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use biometal::operations::{mean_quality, trim_quality_start, trim_quality_window};
use std::sync::atomic::{AtomicU64, Ordering};

/// Quality trimming and filtering module
pub struct QualityTrimmer {
    window_size: usize,
    min_quality: u8,
    min_mean_quality: f64,
    min_length: usize,
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
            stats: AtomicStats::new(),
            failed_quality: AtomicU64::new(0),
            failed_length: AtomicU64::new(0),
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
        if !record.is_empty() {
            if let Ok(trimmed) =
                trim_quality_window(&record.record, self.min_quality, self.window_size)
            {
                let bases_trimmed = record.record.sequence.len() - trimmed.sequence.len();
                if bases_trimmed > 0 {
                    record.record = trimmed;
                    record.metrics.bases_trimmed_3prime += bases_trimmed;
                }
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

        // Should trim the low-quality tail
        assert!(record.len() < 28);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_quality_fail_too_short() {
        let trimmer = QualityTrimmer::new(&make_config());

        // Short read that will be below min_length after trimming
        let seq = b"ATGCATGC";
        let qual = vec![b'I'; 8];
        let mut record = make_record(seq, &qual);
        trimmer.process(&mut record);

        assert!(record.is_failed()); // too short (< 20bp)
    }

    #[test]
    fn test_quality_fail_low_mean() {
        let trimmer = QualityTrimmer::new(&make_config());

        // All low quality
        let seq = vec![b'A'; 30];
        let qual = vec![b'#'; 30]; // Q2, well below min_mean_quality of 20
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
        assert_eq!(record.len(), 100); // no trimming needed
    }
}
