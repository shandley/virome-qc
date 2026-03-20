//! Adapter trimming module — 3' overlap, internal k-mer scan, random primer head trim
//!
//! Handles three types of adapter contamination:
//! 1. **3' adapter read-through**: Standard case where insert < read length
//! 2. **Internal adapter**: Chimeric fragments with adapter sequence embedded
//! 3. **Random primer bias**: Compositional bias from random hexamer/octamer priming

use crate::config::AdapterConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicU64, Ordering};

/// K-mer size for internal adapter scanning
const INTERNAL_KMER_SIZE: usize = 15;

/// Known adapter sequences
struct AdapterSet {
    name: String,
    sequences: Vec<Vec<u8>>,
}

/// Adapter trimming module
pub struct AdapterTrimmer {
    adapter_sets: Vec<AdapterSet>,
    internal_scan: bool,
    /// Precomputed hash set of adapter k-mers + 1-mismatch neighborhood for O(1) lookup
    internal_kmer_hashes: FxHashSet<u64>,
    random_primer_trim: usize,
    min_overlap: usize,
    max_mismatch_rate: f64,
    stats: AtomicStats,
    adapters_found_3prime: AtomicU64,
    adapters_found_internal: AtomicU64,
    primers_trimmed: AtomicU64,
}

impl AdapterTrimmer {
    pub fn new(config: &AdapterConfig) -> Self {
        let adapter_sets: Vec<AdapterSet> = config
            .sequences
            .iter()
            .map(|name| load_adapter_set(name))
            .collect();

        // Precompute hash set of adapter k-mers + 1-mismatch variants for O(1) lookup
        let internal_kmer_hashes = if config.internal_scan {
            build_adapter_kmer_index(&adapter_sets, INTERNAL_KMER_SIZE)
        } else {
            FxHashSet::default()
        };

        Self {
            adapter_sets,
            internal_scan: config.internal_scan,
            internal_kmer_hashes,
            random_primer_trim: config.random_primer_trim,
            min_overlap: config.min_overlap,
            max_mismatch_rate: config.max_mismatch_rate,
            stats: AtomicStats::new(),
            adapters_found_3prime: AtomicU64::new(0),
            adapters_found_internal: AtomicU64::new(0),
            primers_trimmed: AtomicU64::new(0),
        }
    }

    /// Find 3' adapter by overlap alignment
    ///
    /// Scans the 3' end of the read for adapter sequence overlap.
    /// Returns the position where the adapter starts (trim point), or None.
    fn find_3prime_adapter(&self, sequence: &[u8]) -> Option<(usize, String)> {
        for adapter_set in &self.adapter_sets {
            for adapter_seq in &adapter_set.sequences {
                if let Some(pos) =
                    find_adapter_overlap(sequence, adapter_seq, self.min_overlap, self.max_mismatch_rate)
                {
                    return Some((pos, adapter_set.name.clone()));
                }
            }
        }
        None
    }

    /// Scan for internal adapter contamination using precomputed k-mer hash index
    ///
    /// Uses O(1) hash lookups per read k-mer instead of brute-force comparison.
    /// The adapter k-mer index (including 1-mismatch variants) is built once at startup.
    fn has_internal_adapter(&self, sequence: &[u8]) -> bool {
        if sequence.len() < 30 || self.internal_kmer_hashes.is_empty() {
            return false;
        }

        // Only scan the internal portion (skip first/last 15bp to avoid edge effects)
        let scan_start = INTERNAL_KMER_SIZE;
        let scan_end = sequence.len().saturating_sub(INTERNAL_KMER_SIZE);

        if scan_end <= scan_start + INTERNAL_KMER_SIZE {
            return false;
        }

        for kmer in sequence[scan_start..scan_end].windows(INTERNAL_KMER_SIZE) {
            let hash = hash_kmer(kmer);
            if self.internal_kmer_hashes.contains(&hash) {
                return true;
            }
        }
        false
    }
}

impl QcModule for AdapterTrimmer {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();
        let original_len = record.len();

        // Step 1: Random primer head trim (fixed 5' trim)
        if self.random_primer_trim > 0 && record.len() > self.random_primer_trim {
            let trim_amount = self.random_primer_trim;
            record.record.sequence = record.record.sequence[trim_amount..].to_vec();
            record.record.quality = record.record.quality[trim_amount..].to_vec();
            record.metrics.bases_trimmed_5prime += trim_amount;
            self.primers_trimmed.fetch_add(1, Ordering::Relaxed);
        }

        // Step 2: 3' adapter detection and trimming
        if let Some((trim_pos, adapter_name)) = self.find_3prime_adapter(&record.record.sequence) {
            let trimmed = record.record.sequence.len() - trim_pos;
            record.record.sequence.truncate(trim_pos);
            record.record.quality.truncate(trim_pos);
            record.metrics.bases_trimmed_3prime += trimmed;
            record.metrics.adapter_detected = Some(adapter_name);
            self.adapters_found_3prime.fetch_add(1, Ordering::Relaxed);
        }

        // Step 3: Internal adapter scan (flag, don't remove — the read is chimeric)
        if self.internal_scan && self.has_internal_adapter(&record.record.sequence) {
            record.metrics.internal_adapter = true;
            record.fail("internal_adapter_contamination");
            self.adapters_found_internal.fetch_add(1, Ordering::Relaxed);
        }

        // Track modification
        let bases_removed = original_len - record.len();
        if bases_removed > 0 {
            self.stats.record_modified(bases_removed);
        }
        if record.is_failed() {
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "adapters_found_3prime": self.adapters_found_3prime.load(Ordering::Relaxed),
                "adapters_found_internal": self.adapters_found_internal.load(Ordering::Relaxed),
                "primers_trimmed": self.primers_trimmed.load(Ordering::Relaxed),
            }),
        )
    }

    fn name(&self) -> &str {
        "adapter"
    }
}

/// Find adapter overlap at the 3' end of a read
///
/// Checks all possible overlap positions from the end of the read.
/// Returns the position in the read where the adapter starts.
fn find_adapter_overlap(
    read: &[u8],
    adapter: &[u8],
    min_overlap: usize,
    max_mismatch_rate: f64,
) -> Option<usize> {
    if read.is_empty() || adapter.is_empty() || min_overlap == 0 {
        return None;
    }

    let read_len = read.len();
    let adapter_len = adapter.len();

    // Check overlaps from longest to shortest
    let max_overlap = read_len.min(adapter_len);

    for overlap in (min_overlap..=max_overlap).rev() {
        let read_start = read_len - overlap;
        let read_region = &read[read_start..];
        let adapter_region = &adapter[..overlap];

        let mismatches = read_region
            .iter()
            .zip(adapter_region.iter())
            .filter(|(a, b)| !bases_match(**a, **b))
            .count();

        let mismatch_rate = mismatches as f64 / overlap as f64;
        if mismatch_rate <= max_mismatch_rate {
            return Some(read_start);
        }
    }

    None
}

/// Check if two bases match (case-insensitive, handles N)
#[inline]
fn bases_match(a: u8, b: u8) -> bool {
    if a == b'N' || b == b'N' {
        return true; // N matches anything
    }
    a.eq_ignore_ascii_case(&b)
}

/// Check if two k-mers match with at most `max_mismatches` differences
#[cfg(test)]
fn kmers_match(a: &[u8], b: &[u8], max_mismatches: usize) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mismatches = a
        .iter()
        .zip(b.iter())
        .filter(|(x, y)| !bases_match(**x, **y))
        .count();
    mismatches <= max_mismatches
}

/// Load adapter sequences by name
/// Hash a k-mer to u64 (FNV-1a)
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Build a hash set of all adapter k-mers plus 1-mismatch variants
///
/// For each k-mer in each adapter, generates all single-base substitution
/// variants and adds their hashes to the set. This allows O(1) lookup
/// with 1-mismatch tolerance during scanning.
fn build_adapter_kmer_index(adapter_sets: &[AdapterSet], k: usize) -> FxHashSet<u64> {
    let mut hashes = FxHashSet::default();
    let bases = [b'A', b'C', b'G', b'T'];

    for adapter_set in adapter_sets {
        for adapter_seq in &adapter_set.sequences {
            if adapter_seq.len() < k {
                continue;
            }
            for kmer in adapter_seq.windows(k) {
                // Add exact match
                hashes.insert(hash_kmer(kmer));

                // Add all 1-mismatch variants
                let mut variant = kmer.to_vec();
                for pos in 0..k {
                    let original = variant[pos];
                    for &base in &bases {
                        if base != original.to_ascii_uppercase() {
                            variant[pos] = base;
                            hashes.insert(hash_kmer(&variant));
                            variant[pos] = original;
                        }
                    }
                }
            }
        }
    }

    hashes
}

fn load_adapter_set(name: &str) -> AdapterSet {
    match name {
        "truseq" | "truseq_universal" => AdapterSet {
            name: "TruSeq".into(),
            sequences: vec![
                // TruSeq Universal Adapter
                b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
                // TruSeq Adapter, Index
                b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(),
            ],
        },
        "nextera" => AdapterSet {
            name: "Nextera".into(),
            sequences: vec![
                // Nextera Transposase Adapter A
                b"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG".to_vec(),
                // Nextera Transposase Adapter B
                b"GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG".to_vec(),
            ],
        },
        "nebnext" => AdapterSet {
            name: "NEBNext".into(),
            sequences: vec![
                // NEBNext Adapter
                b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
                // NEBNext Universal Adapter
                b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_vec(),
            ],
        },
        _ => {
            log::warn!("Unknown adapter set '{}', using empty set", name);
            AdapterSet {
                name: name.into(),
                sequences: vec![],
            }
        }
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
    fn test_find_adapter_overlap_exact() {
        let adapter = b"AGATCGGAAGAGC";
        // Read with adapter at 3' end
        let read = b"ATGCATGCATGCAGATCGGAAGAGC";
        let pos = find_adapter_overlap(read, adapter, 10, 0.1);
        assert_eq!(pos, Some(12)); // adapter starts at position 12
    }

    #[test]
    fn test_find_adapter_overlap_partial() {
        let adapter = b"AGATCGGAAGAGCACACGTCTG";
        // Read with partial adapter overlap (only first 13 bases of adapter)
        let read = b"ATGCATGCATGCAGATCGGAAGAGC";
        let pos = find_adapter_overlap(read, adapter, 10, 0.1);
        assert!(pos.is_some());
    }

    #[test]
    fn test_find_adapter_overlap_none() {
        let adapter = b"AGATCGGAAGAGC";
        let read = b"ATGCATGCATGCATGCATGC";
        let pos = find_adapter_overlap(read, adapter, 10, 0.1);
        assert_eq!(pos, None);
    }

    #[test]
    fn test_find_adapter_with_mismatches() {
        let adapter = b"AGATCGGAAGAGC";
        // Read with one mismatch in adapter region
        let read = b"ATGCATGCAGATCGNAAGAGC";
        let pos = find_adapter_overlap(read, adapter, 10, 0.15);
        assert!(pos.is_some());
    }

    #[test]
    fn test_adapter_trimmer_3prime() {
        let config = AdapterConfig {
            enabled: true,
            sequences: vec!["truseq".into()],
            internal_scan: false,
            random_primer_trim: 0,
            min_overlap: 10,
            max_mismatch_rate: 0.1,
        };

        let trimmer = AdapterTrimmer::new(&config);
        // Read with TruSeq adapter at position 20
        let mut record = make_record(b"ATGCATGCATGCATGCATGCAGATCGGAAGAGCACACGTCTG");
        trimmer.process(&mut record);

        assert!(record.len() < 42, "Read should be trimmed");
        assert!(!record.is_failed());
        assert!(record.metrics.adapter_detected.is_some());
    }

    #[test]
    fn test_random_primer_trim() {
        let config = AdapterConfig {
            enabled: true,
            sequences: vec![],
            internal_scan: false,
            random_primer_trim: 8,
            min_overlap: 10,
            max_mismatch_rate: 0.1,
        };

        let trimmer = AdapterTrimmer::new(&config);
        let mut record = make_record(b"NNNNNNNNATTGCATGCATGCATGC");
        let original_len = record.len();
        trimmer.process(&mut record);

        assert_eq!(record.len(), original_len - 8);
        assert_eq!(record.metrics.bases_trimmed_5prime, 8);
    }

    #[test]
    fn test_kmers_match_exact() {
        assert!(kmers_match(b"ATGCATGC", b"ATGCATGC", 0));
    }

    #[test]
    fn test_kmers_match_with_mismatch() {
        assert!(kmers_match(b"ATGCATGC", b"ATGNATGC", 1)); // N matches, 1 real mismatch (C->A at pos 4)
        assert!(!kmers_match(b"ATGCATGC", b"TTCCATGC", 1)); // 2 mismatches > max of 1
    }
}
