//! N-base filtering module
//!
//! Removes reads with excessive ambiguous bases (N). High-N reads are
//! uninformative and can pollute k-mer indices and alignment results.
//! Some sequencers produce runs of N when clusters fail.

use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use std::sync::atomic::{AtomicU64, Ordering};

/// N-base filtering module
pub struct NFilter {
    max_n_fraction: f64,
    stats: AtomicStats,
    failed_n: AtomicU64,
    total_n_bases: AtomicU64,
}

impl NFilter {
    pub fn new(max_n_fraction: f64) -> Self {
        Self {
            max_n_fraction,
            stats: AtomicStats::new(),
            failed_n: AtomicU64::new(0),
            total_n_bases: AtomicU64::new(0),
        }
    }
}

impl QcModule for NFilter {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        if record.is_empty() {
            return;
        }

        let n_count = record
            .record
            .sequence
            .iter()
            .filter(|&&b| b == b'N' || b == b'n')
            .count();

        self.total_n_bases
            .fetch_add(n_count as u64, Ordering::Relaxed);

        let n_fraction = n_count as f64 / record.len() as f64;
        if n_fraction > self.max_n_fraction {
            record.fail("excessive_n_bases");
            self.failed_n.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "failed_n_content": self.failed_n.load(Ordering::Relaxed),
                "total_n_bases_observed": self.total_n_bases.load(Ordering::Relaxed),
                "max_n_fraction_threshold": self.max_n_fraction,
            }),
        )
    }

    fn name(&self) -> &str {
        "n_filter"
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
    fn test_pass_clean_read() {
        let filter = NFilter::new(0.10);
        let mut record = make_record(b"ATGCATGCATGCATGCATGCATGCATGC");
        filter.process(&mut record);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_fail_high_n_read() {
        let filter = NFilter::new(0.10);
        // 50% N's — well above threshold
        let mut record = make_record(b"NNNNNNNNNNNNNNNATGCATGCATGCATG");
        filter.process(&mut record);
        assert!(record.is_failed());
    }

    #[test]
    fn test_pass_borderline_n_read() {
        let filter = NFilter::new(0.10);
        // 2 N's in 30bp = 6.7% — below 10% threshold
        let mut record = make_record(b"ATGCATGCNATGCATGCATGCATGCNATG");
        filter.process(&mut record);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_all_n_read() {
        let filter = NFilter::new(0.10);
        let mut record = make_record(b"NNNNNNNNNNNNNNNNNNNN");
        filter.process(&mut record);
        assert!(record.is_failed());
    }
}
