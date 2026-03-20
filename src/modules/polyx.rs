//! Poly-X detection and trimming module
//!
//! Detects and trims homopolymer runs at the 3' end of reads.
//! Platform-aware: aggressive poly-G trimming for NextSeq/NovaSeq two-color chemistry
//! where G is called when signal drops (high-quality poly-G artifacts).

use crate::config::PolyXConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use std::sync::atomic::{AtomicU64, Ordering};

/// Poly-X trimming module
pub struct PolyXTrimmer {
    /// When true, use shorter threshold for poly-G (NovaSeq/NextSeq artifact)
    platform_aware: bool,
    min_length: usize,
    stats: AtomicStats,
    poly_g_trimmed: AtomicU64,
    poly_other_trimmed: AtomicU64,
    total_polyx_bases: AtomicU64,
}

impl PolyXTrimmer {
    pub fn new(config: &PolyXConfig) -> Self {
        Self {
            platform_aware: config.platform_aware,
            min_length: config.min_length,
            stats: AtomicStats::new(),
            poly_g_trimmed: AtomicU64::new(0),
            poly_other_trimmed: AtomicU64::new(0),
            total_polyx_bases: AtomicU64::new(0),
        }
    }
}

impl QcModule for PolyXTrimmer {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.is_empty() {
            return;
        }

        // Platform-aware: use shorter threshold for poly-G (NovaSeq/NextSeq two-color artifact)
        let effective_min = if self.platform_aware {
            // Poly-G gets a more aggressive (shorter) threshold since it's a known artifact
            // Other bases use the standard threshold
            self.min_length.saturating_sub(4).max(6)
        } else {
            self.min_length
        };

        // Scan from 3' end for homopolymer runs
        let (trim_len, base) = find_3prime_polyx(seq, effective_min);

        // If not platform-aware, apply standard threshold to all bases equally.
        // If platform-aware but non-G base, apply the standard (stricter) threshold.
        let trim_len = if self.platform_aware && base != b'G' && trim_len > 0 {
            let (standard_trim, _) = find_3prime_polyx(seq, self.min_length);
            standard_trim
        } else {
            trim_len
        };

        if trim_len > 0 {
            let new_len = seq.len() - trim_len;
            record.record.sequence.truncate(new_len);
            record.record.quality.truncate(new_len);
            record.metrics.bases_trimmed_3prime += trim_len;
            record.metrics.polyx_trimmed = trim_len;

            self.total_polyx_bases
                .fetch_add(trim_len as u64, Ordering::Relaxed);
            self.stats.record_modified(trim_len);

            if base == b'G' {
                self.poly_g_trimmed.fetch_add(1, Ordering::Relaxed);
            } else {
                self.poly_other_trimmed.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    fn report(&self) -> ModuleReport {
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "poly_g_trimmed": self.poly_g_trimmed.load(Ordering::Relaxed),
                "poly_other_trimmed": self.poly_other_trimmed.load(Ordering::Relaxed),
                "total_polyx_bases": self.total_polyx_bases.load(Ordering::Relaxed),
            }),
        )
    }

    fn name(&self) -> &str {
        "polyx"
    }
}

/// Find the longest homopolymer run at the 3' end of a sequence
///
/// Returns (number of bases to trim, the repeated base)
fn find_3prime_polyx(seq: &[u8], min_length: usize) -> (usize, u8) {
    if seq.is_empty() {
        return (0, 0);
    }

    let last_base = seq[seq.len() - 1];
    let mut run_len = 0;

    // Count consecutive identical bases from the 3' end
    for &base in seq.iter().rev() {
        if base.eq_ignore_ascii_case(&last_base) {
            run_len += 1;
        } else {
            break;
        }
    }

    if run_len >= min_length {
        (run_len, last_base)
    } else {
        (0, 0)
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
    fn test_find_poly_g_tail() {
        let seq = b"ATGCATGCATGCGGGGGGGGGGGGG";
        let (trim_len, base) = find_3prime_polyx(seq, 10);
        assert_eq!(trim_len, 13);
        assert_eq!(base, b'G');
    }

    #[test]
    fn test_find_poly_a_tail() {
        let seq = b"ATGCATGCAAAAAAAAAAAAAA";
        let (trim_len, base) = find_3prime_polyx(seq, 10);
        assert_eq!(trim_len, 14);
        assert_eq!(base, b'A');
    }

    #[test]
    fn test_no_polyx() {
        let seq = b"ATGCATGCATGCATGCATGC";
        let (trim_len, _) = find_3prime_polyx(seq, 10);
        assert_eq!(trim_len, 0);
    }

    #[test]
    fn test_short_run_below_threshold() {
        let seq = b"ATGCATGCGGGGG";
        let (trim_len, _) = find_3prime_polyx(seq, 10);
        assert_eq!(trim_len, 0); // 5 G's < threshold of 10
    }

    #[test]
    fn test_polyx_trimmer_module() {
        let config = PolyXConfig {
            enabled: true,
            platform_aware: true,
            min_length: 10,
        };

        let trimmer = PolyXTrimmer::new(&config);
        let mut record = make_record(b"ATGCATGCATGCATGCGGGGGGGGGGGGGGG");
        let original_len = record.len();
        trimmer.process(&mut record);

        assert!(record.len() < original_len);
        assert!(record.metrics.polyx_trimmed > 0);
        assert!(!record.is_failed());
    }
}
