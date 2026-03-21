//! Contaminant screening module -- rRNA, PhiX, and vector detection/removal
//!
//! Uses k-mer hash indices to classify reads as contaminant or pass-through.
//! Architecture: unified index for fast screening, per-category indices for
//! reporting breakdown. Same FxHashSet pattern as internal adapter detection.
//!
//! Contaminant categories:
//! - rRNA (16S/23S/5S prokaryotic, 18S/28S/5.8S eukaryotic)
//! - PhiX174 spike-in control
//! - Cloning vectors (pUC19)

use super::contaminant_refs::*;
use crate::config::ContaminantConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use rustc_hash::FxHashSet;
use std::sync::atomic::{AtomicU64, Ordering};

/// K-mer size for contaminant screening
const CONTAMINANT_K: usize = 21;

/// Contaminant screening module
pub struct ContaminantScreener {
    /// Unified index containing all enabled contaminant k-mers (fast path)
    unified_index: FxHashSet<u64>,
    /// Per-category indices for classification reporting (slow path, only on hits)
    rrna_index: FxHashSet<u64>,
    phix_index: FxHashSet<u64>,
    vector_index: FxHashSet<u64>,
    /// Per-subunit rRNA indices for detailed breakdown
    rrna_subunit_indices: Vec<(String, FxHashSet<u64>)>,
    /// Minimum fraction of read k-mers matching to classify as contaminant
    min_kmer_fraction: f64,
    /// Standard stats
    stats: AtomicStats,
    /// Per-category removal counters
    rrna_removed: AtomicU64,
    phix_removed: AtomicU64,
    vector_removed: AtomicU64,
    /// Per-subunit rRNA counters (prokaryotic vs eukaryotic)
    rrna_prokaryotic: AtomicU64,
    rrna_eukaryotic: AtomicU64,
}

impl ContaminantScreener {
    /// Build the contaminant screener from config
    ///
    /// Constructs k-mer hash indices for all enabled contaminant categories.
    /// Index construction takes ~1ms for the embedded reference set.
    pub fn new(config: &ContaminantConfig) -> Self {
        let mut unified = FxHashSet::default();
        let mut rrna_idx = FxHashSet::default();
        let mut phix_idx = FxHashSet::default();
        let mut vector_idx = FxHashSet::default();

        let mut rrna_subunit_indices = Vec::new();

        if config.screen_rrna {
            for (name, seq) in RRNA_REFS {
                let mut subunit_idx = FxHashSet::default();
                add_sequence_kmers(seq, CONTAMINANT_K, &mut subunit_idx);
                rrna_idx.extend(subunit_idx.iter());
                rrna_subunit_indices.push((name.to_string(), subunit_idx));
            }
            unified.extend(rrna_idx.iter());
        }

        if config.screen_phix {
            add_sequence_kmers(PHIX174, CONTAMINANT_K, &mut phix_idx);
            unified.extend(phix_idx.iter());
        }

        if config.screen_vectors {
            for (_, seq) in VECTOR_REFS {
                add_sequence_kmers(seq, CONTAMINANT_K, &mut vector_idx);
            }
            unified.extend(vector_idx.iter());
        }

        log::debug!(
            "Contaminant index: {} total k-mers (rRNA: {}, PhiX: {}, vector: {})",
            unified.len(),
            rrna_idx.len(),
            phix_idx.len(),
            vector_idx.len()
        );

        Self {
            unified_index: unified,
            rrna_index: rrna_idx,
            phix_index: phix_idx,
            vector_index: vector_idx,
            rrna_subunit_indices,
            min_kmer_fraction: config.min_kmer_fraction,
            stats: AtomicStats::new(),
            rrna_removed: AtomicU64::new(0),
            phix_removed: AtomicU64::new(0),
            vector_removed: AtomicU64::new(0),
            rrna_prokaryotic: AtomicU64::new(0),
            rrna_eukaryotic: AtomicU64::new(0),
        }
    }

    /// Classify a contaminant read into a specific category
    ///
    /// Called only when the unified index has enough hits. Does a second pass
    /// against per-category indices to determine which contaminant type.
    /// For rRNA hits, also determines prokaryotic vs eukaryotic.
    fn classify_contaminant(&self, sequence: &[u8]) -> &'static str {
        if sequence.len() < CONTAMINANT_K {
            return "unknown_contaminant";
        }

        let mut rrna_hits = 0u32;
        let mut phix_hits = 0u32;
        let mut vector_hits = 0u32;

        for kmer in sequence.windows(CONTAMINANT_K) {
            let hash = hash_kmer(kmer);
            if self.rrna_index.contains(&hash) {
                rrna_hits += 1;
            }
            if self.phix_index.contains(&hash) {
                phix_hits += 1;
            }
            if self.vector_index.contains(&hash) {
                vector_hits += 1;
            }
        }

        if rrna_hits >= phix_hits && rrna_hits >= vector_hits {
            // Determine prokaryotic vs eukaryotic by checking subunit indices
            let mut prok_hits = 0u32;
            let mut euk_hits = 0u32;
            for (name, idx) in &self.rrna_subunit_indices {
                let hits: u32 = sequence
                    .windows(CONTAMINANT_K)
                    .filter(|kmer| idx.contains(&hash_kmer(kmer)))
                    .count() as u32;
                match name.as_str() {
                    "16S" | "23S" | "5S" => prok_hits += hits,
                    "18S" | "28S" | "5.8S" => euk_hits += hits,
                    _ => {}
                }
            }
            if prok_hits >= euk_hits {
                self.rrna_prokaryotic.fetch_add(1, Ordering::Relaxed);
            } else {
                self.rrna_eukaryotic.fetch_add(1, Ordering::Relaxed);
            }
            "contaminant_rrna"
        } else if phix_hits >= vector_hits {
            "contaminant_phix"
        } else {
            "contaminant_vector"
        }
    }
}

impl QcModule for ContaminantScreener {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.len() < CONTAMINANT_K || self.unified_index.is_empty() {
            return;
        }

        // Fast path: count k-mer hits against unified index
        let total_kmers = seq.len() - CONTAMINANT_K + 1;
        let mut hits = 0u32;

        for kmer in seq.windows(CONTAMINANT_K) {
            if self.unified_index.contains(&hash_kmer(kmer)) {
                hits += 1;
            }
        }

        let hit_fraction = hits as f64 / total_kmers as f64;

        if hit_fraction >= self.min_kmer_fraction {
            // Slow path: determine which category
            let category = self.classify_contaminant(seq);

            match category {
                "contaminant_rrna" => self.rrna_removed.fetch_add(1, Ordering::Relaxed),
                "contaminant_phix" => self.phix_removed.fetch_add(1, Ordering::Relaxed),
                "contaminant_vector" => self.vector_removed.fetch_add(1, Ordering::Relaxed),
                _ => 0,
            };

            record.fail(category);
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        let rrna = self.rrna_removed.load(Ordering::Relaxed);
        let phix = self.phix_removed.load(Ordering::Relaxed);
        let vector = self.vector_removed.load(Ordering::Relaxed);

        let rrna_prok = self.rrna_prokaryotic.load(Ordering::Relaxed);
        let rrna_euk = self.rrna_eukaryotic.load(Ordering::Relaxed);

        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "rrna_removed": rrna,
                "rrna_prokaryotic": rrna_prok,
                "rrna_eukaryotic": rrna_euk,
                "phix_removed": phix,
                "vector_removed": vector,
                "index_size_total": self.unified_index.len(),
                "index_size_rrna": self.rrna_index.len(),
                "index_size_phix": self.phix_index.len(),
                "index_size_vector": self.vector_index.len(),
            }),
        )
    }

    fn name(&self) -> &str {
        "contaminant"
    }
}

/// Hash a k-mer using FNV-1a (uppercase, same as adapter module)
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Add all k-mers from a sequence to a hash set (both orientations)
fn add_sequence_kmers(seq: &[u8], k: usize, set: &mut FxHashSet<u64>) {
    if seq.len() < k {
        return;
    }
    for kmer in seq.windows(k) {
        // Skip k-mers containing N
        if kmer.iter().any(|&b| b == b'N' || b == b'n') {
            continue;
        }
        // Add forward hash
        set.insert(hash_kmer(kmer));
        // Add reverse complement hash
        set.insert(hash_kmer_rc(kmer));
    }
}

/// Hash the reverse complement of a k-mer
#[inline]
fn hash_kmer_rc(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer.iter().rev() {
        let comp = match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        };
        hash ^= comp as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

#[cfg(test)]
mod tests {
    use super::*;
    use biometal::FastqRecord;

    fn make_config(rrna: bool, phix: bool, vectors: bool) -> ContaminantConfig {
        ContaminantConfig {
            enabled: true,
            screen_rrna: rrna,
            screen_phix: phix,
            screen_vectors: vectors,
            screen_kitome: false,
            min_kmer_fraction: 0.4,
        }
    }

    fn make_record(seq: &[u8]) -> AnnotatedRecord {
        let qual = vec![b'I'; seq.len()];
        AnnotatedRecord::new(FastqRecord::new("test".into(), seq.to_vec(), qual))
    }

    #[test]
    fn test_index_construction() {
        let screener = ContaminantScreener::new(&make_config(true, true, true));
        assert!(
            screener.unified_index.len() > 1000,
            "Should have substantial index"
        );
        assert!(screener.rrna_index.len() > 0);
        assert!(screener.phix_index.len() > 0);
        assert!(screener.vector_index.len() > 0);
    }

    #[test]
    fn test_detect_16s_rrna() {
        let screener = ContaminantScreener::new(&make_config(true, false, false));
        // Take a 150bp fragment from 16S
        let fragment = &ECOLI_16S[100..250];
        let mut record = make_record(fragment);
        screener.process(&mut record);
        assert!(record.is_failed(), "16S rRNA read should be detected");
    }

    #[test]
    fn test_detect_23s_rrna() {
        let screener = ContaminantScreener::new(&make_config(true, false, false));
        let fragment = &ECOLI_23S[100..250];
        let mut record = make_record(fragment);
        screener.process(&mut record);
        assert!(record.is_failed(), "23S rRNA read should be detected");
    }

    #[test]
    fn test_detect_18s_rrna() {
        let screener = ContaminantScreener::new(&make_config(true, false, false));
        let fragment = &HUMAN_18S[100..250];
        let mut record = make_record(fragment);
        screener.process(&mut record);
        assert!(record.is_failed(), "18S rRNA read should be detected");
    }

    #[test]
    fn test_detect_phix() {
        let screener = ContaminantScreener::new(&make_config(false, true, false));
        let fragment = &PHIX174[200..350];
        let mut record = make_record(fragment);
        screener.process(&mut record);
        assert!(record.is_failed(), "PhiX read should be detected");
    }

    #[test]
    fn test_detect_vector() {
        let screener = ContaminantScreener::new(&make_config(false, false, true));
        let fragment = &PUC19[100..250];
        let mut record = make_record(fragment);
        screener.process(&mut record);
        assert!(record.is_failed(), "pUC19 read should be detected");
    }

    #[test]
    fn test_random_read_passes() {
        let screener = ContaminantScreener::new(&make_config(true, true, true));
        // Random sequence should not match contaminants
        let random: Vec<u8> = (0..150)
            .map(|i| [b'A', b'C', b'G', b'T'][(i * 7 + 3) % 4])
            .collect();
        let mut record = make_record(&random);
        screener.process(&mut record);
        assert!(!record.is_failed(), "Random read should pass");
    }

    #[test]
    fn test_partial_match_below_threshold() {
        let screener = ContaminantScreener::new(&make_config(true, false, false));
        // Create a chimeric read: 50bp rRNA + 100bp random
        let mut chimera = ECOLI_16S[100..150].to_vec();
        let random: Vec<u8> = (0..100)
            .map(|i| [b'A', b'C', b'G', b'T'][(i * 7 + 3) % 4])
            .collect();
        chimera.extend_from_slice(&random);

        let mut record = make_record(&chimera);
        screener.process(&mut record);
        // 50/150 = 33% rRNA k-mers, below 40% threshold
        assert!(
            !record.is_failed(),
            "Partial match below threshold should pass"
        );
    }

    #[test]
    fn test_category_classification() {
        let screener = ContaminantScreener::new(&make_config(true, true, true));

        // rRNA read
        let mut rrna_record = make_record(&ECOLI_16S[100..250]);
        screener.process(&mut rrna_record);
        assert!(rrna_record.is_failed());

        // PhiX read
        let mut phix_record = make_record(&PHIX174[200..350]);
        screener.process(&mut phix_record);
        assert!(phix_record.is_failed());

        let report = screener.report();
        let rrna_count = report
            .extra
            .get("rrna_removed")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        let phix_count = report
            .extra
            .get("phix_removed")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        assert_eq!(rrna_count, 1, "Should count 1 rRNA removal");
        assert_eq!(phix_count, 1, "Should count 1 PhiX removal");
    }

    #[test]
    fn test_disabled_categories() {
        // Only PhiX enabled -- rRNA should pass
        let screener = ContaminantScreener::new(&make_config(false, true, false));
        let mut record = make_record(&ECOLI_16S[100..250]);
        screener.process(&mut record);
        assert!(
            !record.is_failed(),
            "rRNA should pass when screening disabled"
        );
    }
}
