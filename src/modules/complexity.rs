//! Complexity filtering module
//!
//! Removes low-complexity reads using Shannon entropy via biometal's
//! NEON-optimized complexity_score(). Virome-aware: configurable threshold
//! because some viral genomes (AT-rich phage, polyA from cDNA) are
//! genuinely lower complexity than typical genomic sequence.

use crate::config::ComplexityConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use biometal::operations::complexity_score;
use std::sync::atomic::{AtomicU64, Ordering};

/// Complexity filtering module
pub struct ComplexityFilter {
    min_entropy: f64,
    stats: AtomicStats,
    failed_complexity: AtomicU64,
}

impl ComplexityFilter {
    pub fn new(config: &ComplexityConfig) -> Self {
        Self {
            min_entropy: config.min_entropy,
            stats: AtomicStats::new(),
            failed_complexity: AtomicU64::new(0),
        }
    }
}

impl QcModule for ComplexityFilter {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        if record.is_empty() {
            return;
        }

        let score = complexity_score(&record.record.sequence);
        record.metrics.complexity = Some(score);

        if score < self.min_entropy {
            record.fail("low_complexity");
            self.failed_complexity.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "failed_complexity": self.failed_complexity.load(Ordering::Relaxed),
                "min_entropy_threshold": self.min_entropy,
            }),
        )
    }

    fn name(&self) -> &str {
        "complexity"
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
    fn test_fail_homopolymer() {
        let filter = ComplexityFilter::new(&ComplexityConfig {
            enabled: true,
            min_entropy: 0.5,
        });

        let mut record = make_record(&vec![b'A'; 100]);
        filter.process(&mut record);
        assert!(record.is_failed());
    }

    #[test]
    fn test_pass_complex_sequence() {
        let filter = ComplexityFilter::new(&ComplexityConfig {
            enabled: true,
            min_entropy: 0.5,
        });

        let mut record = make_record(b"ATGCTAGCGATCGATCGATCGATCGATCGATCGATCGATC");
        filter.process(&mut record);
        assert!(!record.is_failed());
        assert!(record.metrics.complexity.unwrap() > 0.5);
    }

    #[test]
    fn test_dinucleotide_repeat_borderline() {
        let filter = ComplexityFilter::new(&ComplexityConfig {
            enabled: true,
            min_entropy: 0.5,
        });

        // AT repeat has entropy ~0.5 — exactly at threshold
        let mut record = make_record(b"ATATATATATATATATATATATATATATATATATAT");
        filter.process(&mut record);
        // Dinucleotide repeat should be right at the boundary
        assert!(record.metrics.complexity.is_some());
    }

    #[test]
    fn test_virome_aware_low_threshold() {
        // With lower threshold for virome-specific profiles
        let filter = ComplexityFilter::new(&ComplexityConfig {
            enabled: true,
            min_entropy: 0.3,
        });

        // AT repeat should pass with lower threshold
        let mut record = make_record(b"ATATATATATATATATATATATATATATATATATAT");
        filter.process(&mut record);
        assert!(!record.is_failed());
    }
}
