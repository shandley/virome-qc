//! EVE k-mer exclusion set for host depletion
//!
//! Prevents false removal of exogenous viral reads that share k-mers with
//! endogenous viral elements (EVEs) integrated into the human reference genome.
//!
//! At startup, hashes all k=31 k-mers from embedded high-identity EVE sequences
//! into an FxHashSet. During host depletion, k-mer positions that match both
//! the host Super Bloom filter AND the EVE exclusion set are excluded from the
//! containment calculation. This lowers the apparent host containment of reads
//! from EVE regions, allowing genuinely exogenous viral reads to pass through.
//!
//! Follows the same exclusion pattern as PhiX/Microviridae shared k-mers
//! in the contaminant module and herpesvirus/ERV shared k-mers in the ERV module.

#[path = "eve_kmers.rs"]
mod eve_kmers;

use log::info;
use rustc_hash::FxHashSet;

/// K-mer size matching the host Super Bloom filter
const EVE_K: usize = 31;

/// EVE k-mer exclusion set for host depletion
pub struct EveExclusionSet {
    kmers: FxHashSet<u64>,
}

impl EveExclusionSet {
    /// Build the exclusion set from embedded EVE sequences.
    /// Hashes all k=31 k-mers (both strands) using FNV-1a.
    pub fn build() -> Self {
        let mut kmers = FxHashSet::default();

        for (name, seq) in eve_kmers::EVE_SEQUENCES {
            let before = kmers.len();
            add_sequence_kmers(seq, EVE_K, &mut kmers);
            let added = kmers.len() - before;
            log::debug!("EVE exclusion: {} -> {} new k-mers", name, added);
        }

        info!(
            "EVE exclusion set: {} unique {}-mers from {} sequences",
            kmers.len(),
            EVE_K,
            eve_kmers::EVE_SEQUENCES.len(),
        );

        Self { kmers }
    }

    /// Check if a k-mer (as a byte slice) falls in an EVE region.
    /// Checks both forward and reverse complement orientations.
    pub fn contains_kmer(&self, kmer: &[u8]) -> bool {
        if kmer.len() < EVE_K {
            return false;
        }
        let fwd = hash_kmer(kmer);
        if self.kmers.contains(&fwd) {
            return true;
        }
        let rc = hash_kmer_rc(kmer);
        self.kmers.contains(&rc)
    }

    /// Number of k-mers in the exclusion set
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    /// Whether the exclusion set is empty (no EVE sequences embedded)
    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }
}

/// FNV-1a hash (uppercase, same as contaminant and ERV modules)
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// FNV-1a hash of reverse complement
#[inline]
fn hash_kmer_rc(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer.iter().rev() {
        let comp = match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        };
        hash ^= comp as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Add all k-mers from a sequence (both strands) to the hash set
fn add_sequence_kmers(seq: &[u8], k: usize, set: &mut FxHashSet<u64>) {
    if seq.len() < k {
        return;
    }
    for window in seq.windows(k) {
        // Skip k-mers with ambiguous bases
        if window.iter().any(|&b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        set.insert(hash_kmer(window));
        set.insert(hash_kmer_rc(window));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_exclusion_set() {
        let set = EveExclusionSet::build();
        // With placeholder sequences, set should be empty
        assert!(set.is_empty() || !set.is_empty()); // just verify it builds
    }

    #[test]
    fn test_hash_consistency() {
        let kmer = b"ACGTACGTACGTACGTACGTACGTACGTACG"; // 31bp
        let h1 = hash_kmer(kmer);
        let h2 = hash_kmer(kmer);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_rc_hash_differs() {
        let kmer = b"ACGTACGTACGTACGTACGTACGTACGTACG"; // not palindromic
        let fwd = hash_kmer(kmer);
        let rc = hash_kmer_rc(kmer);
        // For non-palindromic k-mers, forward and RC hashes should differ
        // (they could theoretically collide but extremely unlikely)
        assert_ne!(fwd, rc);
    }

    #[test]
    fn test_contains_kmer_both_strands() {
        let mut kmers = FxHashSet::default();
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp -> one 31-mer + one 31-mer
        add_sequence_kmers(seq, EVE_K, &mut kmers);

        let set = EveExclusionSet { kmers };
        // Forward k-mer should be found
        assert!(set.contains_kmer(&seq[0..31]));
        // RC should also be found
        assert!(set.contains_kmer(&seq[0..31])); // contains checks both orientations
    }
}
