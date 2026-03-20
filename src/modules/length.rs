//! Final length filter module
//!
//! Applied as the last step after all trimming modules (adapter, poly-X, quality).
//! Catches reads that were shortened below the useful threshold by the cumulative
//! effect of multiple trimming operations.
//!
//! This is separate from the quality module's internal length check because
//! poly-X trimming happens after quality trimming and can further shorten reads.

use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use std::sync::atomic::{AtomicU64, Ordering};

/// Final length filter module
pub struct LengthFilter {
    min_length: usize,
    stats: AtomicStats,
    failed_length: AtomicU64,
}

impl LengthFilter {
    pub fn new(min_length: usize) -> Self {
        Self {
            min_length,
            stats: AtomicStats::new(),
            failed_length: AtomicU64::new(0),
        }
    }
}

impl QcModule for LengthFilter {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        if record.len() < self.min_length {
            record.fail("too_short_final");
            self.failed_length.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "failed_length": self.failed_length.load(Ordering::Relaxed),
                "min_length_threshold": self.min_length,
            }),
        )
    }

    fn name(&self) -> &str {
        "length_filter"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use biometal::FastqRecord;

    fn make_record(seq: &[u8]) -> AnnotatedRecord {
        let qual = vec![b'I'; seq.len()];
        AnnotatedRecord::new(FastqRecord::new("test".into(), seq.to_vec(), qual))
    }

    #[test]
    fn test_pass_long_read() {
        let filter = LengthFilter::new(50);
        let mut record = make_record(&vec![b'A'; 100]);
        filter.process(&mut record);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_fail_short_read() {
        let filter = LengthFilter::new(50);
        let mut record = make_record(&vec![b'A'; 30]);
        filter.process(&mut record);
        assert!(record.is_failed());
    }

    #[test]
    fn test_exact_threshold() {
        let filter = LengthFilter::new(50);
        let mut record = make_record(&vec![b'A'; 50]);
        filter.process(&mut record);
        assert!(!record.is_failed()); // exactly at threshold = pass
    }
}
