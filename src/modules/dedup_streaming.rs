//! Streaming deduplication module -- integrates into the pipeline
//!
//! Uses a FxHashSet of read hashes to track seen sequences. First occurrence
//! passes, subsequent occurrences are failed as PCR duplicates.
//!
//! For paired-end: both mates are hashed together (R1 prefix + R2 prefix)
//! via the executor's paired processing. Each mate independently hashes
//! against the shared set, so this module works for single-end in the
//! streaming pipeline. Paired-end concordance is handled by the executor.
//!
//! Virome-aware: uses prefix-only hashing (skips first 5bp for trim tolerance)
//! so contained reads from differential trimming are correctly identified.
//!
//! Memory: ~16 bytes per unique read. 100M reads = ~1.6 GB.

use crate::config::DedupConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// Prefix hash skip and length (same as post-pipeline dedup)
const SKIP_BASES: usize = 5;
const PREFIX_LEN: usize = 50;

/// Streaming deduplication module
pub struct StreamingDedup {
    /// Thread-safe set of seen read hashes
    seen: Mutex<FxHashSet<u64>>,
    /// Standard stats
    stats: AtomicStats,
    /// Duplicate reads removed
    duplicates_removed: AtomicU64,
}

impl StreamingDedup {
    pub fn new(_config: &DedupConfig) -> Self {
        Self {
            seen: Mutex::new(FxHashSet::default()),
            stats: AtomicStats::new(),
            duplicates_removed: AtomicU64::new(0),
        }
    }
}

impl QcModule for StreamingDedup {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.len() < SKIP_BASES + 10 {
            return; // too short to hash meaningfully
        }

        let hash = hash_prefix(seq);

        // Try to insert -- if already present, it's a duplicate
        let is_new = {
            let mut set = self.seen.lock().unwrap();
            set.insert(hash)
        };

        if !is_new {
            record.fail("pcr_duplicate");
            self.duplicates_removed.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        let seen_count = self.seen.lock().unwrap().len() as u64;
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "duplicates_removed": self.duplicates_removed.load(Ordering::Relaxed),
                "unique_sequences": seen_count,
                "estimated_library_complexity": seen_count,
            }),
        )
    }

    fn name(&self) -> &str {
        "dedup"
    }
}

/// Hash read prefix (skip first 5 bases for trim tolerance)
fn hash_prefix(sequence: &[u8]) -> u64 {
    let start = SKIP_BASES.min(sequence.len());
    let end = (start + PREFIX_LEN).min(sequence.len());
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in &sequence[start..end] {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::DedupConfig;
    use biometal::FastqRecord;

    fn make_record(seq: &[u8]) -> AnnotatedRecord {
        let qual = vec![b'I'; seq.len()];
        AnnotatedRecord::new(FastqRecord::new("test".into(), seq.to_vec(), qual))
    }

    fn make_dedup() -> StreamingDedup {
        StreamingDedup::new(&DedupConfig {
            enabled: false,
            optical_distance: 2500,
            umi_aware: false,
        })
    }

    #[test]
    fn test_first_occurrence_passes() {
        let dedup = make_dedup();
        let seq: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let mut record = make_record(&seq);
        dedup.process(&mut record);
        assert!(!record.is_failed());
    }

    #[test]
    fn test_duplicate_fails() {
        let dedup = make_dedup();
        let seq: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();

        let mut first = make_record(&seq);
        dedup.process(&mut first);
        assert!(!first.is_failed());

        let mut second = make_record(&seq);
        dedup.process(&mut second);
        assert!(second.is_failed());
    }

    #[test]
    fn test_different_reads_pass() {
        let dedup = make_dedup();
        let seq1: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let seq2: Vec<u8> = (0..100).map(|i| [b'G', b'C', b'T', b'A'][i % 4]).collect();

        let mut r1 = make_record(&seq1);
        let mut r2 = make_record(&seq2);
        dedup.process(&mut r1);
        dedup.process(&mut r2);
        assert!(!r1.is_failed());
        assert!(!r2.is_failed());
    }

    #[test]
    fn test_contained_reads_detected() {
        let dedup = make_dedup();
        // Two reads with same prefix but different lengths (contained)
        let long_seq: Vec<u8> = (0..150).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let short_seq: Vec<u8> = (0..80).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();

        let mut long_read = make_record(&long_seq);
        dedup.process(&mut long_read);
        assert!(!long_read.is_failed());

        // Short read has same prefix -- should be duplicate
        let mut short_read = make_record(&short_seq);
        dedup.process(&mut short_read);
        assert!(
            short_read.is_failed(),
            "Contained read should be detected as duplicate"
        );
    }

    #[test]
    fn test_report_counts() {
        let dedup = make_dedup();
        let seq: Vec<u8> = (0..100).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();

        for _ in 0..5 {
            let mut r = make_record(&seq);
            dedup.process(&mut r);
        }

        let report = dedup.report();
        assert_eq!(report.reads_processed, 5);
        assert_eq!(report.reads_removed, 4); // 1 unique + 4 duplicates
    }
}
