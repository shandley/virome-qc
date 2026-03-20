//! Read merging module for overlapping paired-end reads
//!
//! Merges R1 and R2 reads when they overlap, producing a single longer
//! read with error-corrected overlap region. Quality scores in the
//! overlap are resolved by taking the higher quality base.
//!
//! Uses biometal primitives:
//! - `reverse_complement` for R2 orientation
//! - `hamming_distance` for fast mismatch counting in overlap candidates
//! - `mean_quality` for quality-aware consensus scoring
//!
//! This runs AFTER all QC modules on clean reads.

use biometal::operations::distance::hamming_distance;
use biometal::operations::reverse_complement;
use biometal::FastqRecord;
use std::sync::atomic::{AtomicU64, Ordering};

/// Result of attempting to merge a read pair
pub enum MergeResult {
    /// Reads were merged into a single record
    Merged(FastqRecord),
    /// Reads do not overlap sufficiently — return original pair
    Unmerged(FastqRecord, FastqRecord),
}

/// Configuration for read merging
pub struct MergeConfig {
    /// Minimum overlap length to consider merging
    pub min_overlap: usize,
    /// Maximum mismatch rate in the overlap region
    pub max_mismatch_rate: f64,
}

impl Default for MergeConfig {
    fn default() -> Self {
        Self {
            min_overlap: 15,
            max_mismatch_rate: 0.10,
        }
    }
}

/// Read merger with statistics tracking
pub struct ReadMerger {
    config: MergeConfig,
    pub pairs_attempted: AtomicU64,
    pub pairs_merged: AtomicU64,
    pub pairs_unmerged: AtomicU64,
}

impl ReadMerger {
    pub fn new(config: MergeConfig) -> Self {
        Self {
            config,
            pairs_attempted: AtomicU64::new(0),
            pairs_merged: AtomicU64::new(0),
            pairs_unmerged: AtomicU64::new(0),
        }
    }

    /// Attempt to merge a read pair
    ///
    /// R2 is reverse-complemented and aligned against R1 to find overlap.
    /// If sufficient overlap with acceptable mismatch rate is found,
    /// the reads are merged with quality-aware consensus in the overlap.
    pub fn merge_pair(&self, r1: &FastqRecord, r2: &FastqRecord) -> MergeResult {
        self.pairs_attempted.fetch_add(1, Ordering::Relaxed);

        // Reverse complement R2 to align with R1 orientation
        let r2_rc_seq = reverse_complement(&r2.sequence);
        let r2_rc_qual: Vec<u8> = r2.quality.iter().rev().copied().collect();

        // Find best suffix-prefix overlap using biometal's hamming_distance
        // R1: -------->
        // R2_RC:    -------->
        //        ^^^^ overlap
        if let Some((overlap_start, overlap_len)) = find_best_overlap(
            &r1.sequence,
            &r2_rc_seq,
            self.config.min_overlap,
            self.config.max_mismatch_rate,
        ) {
            // Build merged sequence with quality-aware consensus
            let merged_len = overlap_start + r2_rc_seq.len();
            let mut merged_seq = Vec::with_capacity(merged_len);
            let mut merged_qual = Vec::with_capacity(merged_len);

            // R1 prefix (before overlap)
            merged_seq.extend_from_slice(&r1.sequence[..overlap_start]);
            merged_qual.extend_from_slice(&r1.quality[..overlap_start]);

            // Overlap region — quality-aware consensus
            for i in 0..overlap_len {
                let r1_idx = overlap_start + i;
                let r2_idx = i;

                if r1_idx >= r1.sequence.len() || r2_idx >= r2_rc_seq.len() {
                    break;
                }

                let (base, qual) = consensus_base(
                    r1.sequence[r1_idx],
                    r1.quality[r1_idx],
                    r2_rc_seq[r2_idx],
                    r2_rc_qual[r2_idx],
                );
                merged_seq.push(base);
                merged_qual.push(qual);
            }

            // R2 suffix (beyond overlap)
            let r2_remaining_start = overlap_len;
            if r2_remaining_start < r2_rc_seq.len() {
                merged_seq.extend_from_slice(&r2_rc_seq[r2_remaining_start..]);
                merged_qual.extend_from_slice(&r2_rc_qual[r2_remaining_start..]);
            }

            let merged_id = format!("{} merged_overlap={overlap_len}", r1.id);
            self.pairs_merged.fetch_add(1, Ordering::Relaxed);

            MergeResult::Merged(FastqRecord::new(merged_id, merged_seq, merged_qual))
        } else {
            self.pairs_unmerged.fetch_add(1, Ordering::Relaxed);
            MergeResult::Unmerged(r1.clone(), r2.clone())
        }
    }

    pub fn merge_rate(&self) -> f64 {
        let attempted = self.pairs_attempted.load(Ordering::Relaxed);
        if attempted == 0 {
            return 0.0;
        }
        self.pairs_merged.load(Ordering::Relaxed) as f64 / attempted as f64
    }
}

/// Quality-aware consensus for a single position in the overlap region
///
/// When bases agree: keep the base, take the higher quality score.
/// When bases disagree: take the base with the higher quality score.
#[inline]
fn consensus_base(base1: u8, qual1: u8, base2: u8, qual2: u8) -> (u8, u8) {
    if base1 == base2 {
        // Agreement — boost confidence by taking max quality
        (base1, qual1.max(qual2))
    } else if qual1 >= qual2 {
        (base1, qual1)
    } else {
        (base2, qual2)
    }
}

/// Find the best suffix-prefix overlap between R1 and R2_RC
///
/// Uses biometal's `hamming_distance` for efficient mismatch counting
/// at each candidate overlap position. Iterates from longest to shortest
/// overlap and returns the first (longest) valid overlap found.
///
/// Returns (overlap_start_in_r1, overlap_length) if found.
fn find_best_overlap(
    r1: &[u8],
    r2_rc: &[u8],
    min_overlap: usize,
    max_mismatch_rate: f64,
) -> Option<(usize, usize)> {
    let r1_len = r1.len();
    let r2_len = r2_rc.len();
    let max_overlap = r1_len.min(r2_len);

    // Scan from longest overlap to shortest — return first valid match
    for overlap in (min_overlap..=max_overlap).rev() {
        let r1_start = r1_len - overlap;
        let r1_region = &r1[r1_start..];
        let r2_region = &r2_rc[..overlap];

        // Use biometal's hamming_distance for mismatch counting
        let mismatches = match hamming_distance(r1_region, r2_region) {
            Ok(d) => d,
            Err(_) => continue, // length mismatch, skip
        };
        let max_allowed = (overlap as f64 * max_mismatch_rate).ceil() as usize;

        if mismatches <= max_allowed {
            return Some((r1_start, overlap));
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(id: &str, seq: &[u8]) -> FastqRecord {
        FastqRecord::new(id.into(), seq.to_vec(), vec![b'I'; seq.len()])
    }

    #[test]
    fn test_merge_overlapping_pair() {
        // Non-repetitive fragment to avoid ambiguous overlap
        let fragment = b"ATGCTAGCGATCGTAGCTACGATCGTAG";
        assert_eq!(fragment.len(), 28);

        let r1 = make_record("r1", &fragment[..20]);
        // R2 = reverse complement of last 20bp of fragment
        let r2_seq = reverse_complement(&fragment[8..28]);
        let r2 = make_record("r2", &r2_seq);

        let merger = ReadMerger::new(MergeConfig {
            min_overlap: 10,
            max_mismatch_rate: 0.1,
        });

        match merger.merge_pair(&r1, &r2) {
            MergeResult::Merged(merged) => {
                assert_eq!(
                    merged.sequence.len(),
                    28,
                    "Merged length should be fragment length, got {}",
                    merged.sequence.len()
                );
                assert_eq!(&merged.sequence[..], &fragment[..]);
            }
            MergeResult::Unmerged(_, _) => panic!("Should have merged"),
        }
    }

    #[test]
    fn test_no_merge_when_no_overlap() {
        let r1 = make_record("r1", b"AAAAAAAAAAAAAAAAAAAAAAAAA");
        let r2 = make_record("r2", b"CCCCCCCCCCCCCCCCCCCCCCCCC");

        let merger = ReadMerger::new(MergeConfig::default());
        match merger.merge_pair(&r1, &r2) {
            MergeResult::Unmerged(_, _) => {} // expected
            MergeResult::Merged(_) => panic!("Should not have merged"),
        }
    }

    #[test]
    fn test_overlap_detection_with_hamming() {
        let r1 = b"ATGCATGCATGCATGC";
        let r2_rc = b"ATGCATGCTTTTTTTT"; // first 8bp overlap with r1 suffix

        let result = find_best_overlap(r1, r2_rc, 5, 0.1);
        assert!(result.is_some());
        let (start, len) = result.unwrap();
        assert_eq!(start, 8);
        assert_eq!(len, 8);
    }

    #[test]
    fn test_consensus_base_agreement() {
        // Same base, different quality — take max quality
        let (base, qual) = consensus_base(b'A', b'I', b'A', b'5');
        assert_eq!(base, b'A');
        assert_eq!(qual, b'I'); // higher quality wins
    }

    #[test]
    fn test_consensus_base_disagreement() {
        // Different bases — take the one with higher quality
        let (base, qual) = consensus_base(b'A', b'5', b'G', b'I');
        assert_eq!(base, b'G'); // G has higher quality
        assert_eq!(qual, b'I');
    }

    #[test]
    fn test_merge_statistics() {
        let merger = ReadMerger::new(MergeConfig::default());

        let r1 = make_record("r1", b"AAAAAAAAAAAAAAAAAAAAAAAAA");
        let r2 = make_record("r2", b"CCCCCCCCCCCCCCCCCCCCCCCCC");
        merger.merge_pair(&r1, &r2);

        assert_eq!(merger.pairs_attempted.load(Ordering::Relaxed), 1);
        assert_eq!(merger.pairs_unmerged.load(Ordering::Relaxed), 1);
        assert_eq!(merger.merge_rate(), 0.0);
    }
}
