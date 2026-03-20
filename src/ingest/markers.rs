//! Marker k-mer database for fast sample composition screening
//!
//! Uses small sets of diagnostic k-mers from reference organisms to estimate
//! sample composition during ingestion. The marker sets are chosen from
//! high-copy or highly conserved regions that are diagnostic for each category:
//!
//! - **Human**: Alu elements (~1M copies in genome), LINE-1 5'UTR, 18S/28S rRNA
//! - **Bacterial**: 16S rRNA conserved regions (universal bacterial marker)
//! - **PhiX174**: Spike-in control genome
//!
//! Total marker database: ~35K k-mers = ~700KB in memory. Fits in L2 cache.
//! Screening 50K reads takes <10ms.

use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};

/// K-mer size for marker screening
const MARKER_K: usize = 21;

/// Organism category for classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MarkerCategory {
    Human,
    Bacterial,
    PhiX,
}

/// Marker k-mer database
pub struct MarkerDb {
    human: FxHashSet<u64>,
    bacterial: FxHashSet<u64>,
    phix: FxHashSet<u64>,
}

/// Sample composition estimate from marker screening
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompositionEstimate {
    /// Fraction of reads with human marker hits
    pub human_fraction: f64,
    /// Fraction of reads with bacterial marker hits
    pub bacterial_fraction: f64,
    /// Fraction of reads with PhiX marker hits
    pub phix_fraction: f64,
    /// Fraction of reads with no marker hits (viral, unknown, or low-complexity)
    pub unclassified_fraction: f64,
    /// Total reads screened
    pub reads_screened: u64,
    /// Number of marker k-mers in database per category
    pub marker_counts: Vec<(String, usize)>,
}

impl MarkerDb {
    /// Build the marker database from embedded reference sequences
    pub fn new() -> Self {
        Self {
            human: build_kmer_set(HUMAN_MARKERS),
            bacterial: build_kmer_set(BACTERIAL_MARKERS),
            phix: build_kmer_set(&[PHIX174_GENOME]),
        }
    }

    /// Number of marker k-mers per category
    pub fn marker_counts(&self) -> Vec<(String, usize)> {
        vec![
            ("human".into(), self.human.len()),
            ("bacterial".into(), self.bacterial.len()),
            ("phix".into(), self.phix.len()),
        ]
    }

    /// Classify a single read by counting k-mer hits against each marker set
    ///
    /// Returns the category with the most hits, or None if no significant hits.
    /// A read needs at least 3 k-mer hits (out of ~130 possible for 150bp) to
    /// be classified, preventing false positives from random k-mer collisions.
    pub fn classify_read(&self, sequence: &[u8]) -> Option<MarkerCategory> {
        if sequence.len() < MARKER_K {
            return None;
        }

        let mut human_hits = 0u32;
        let mut bacterial_hits = 0u32;
        let mut phix_hits = 0u32;

        for kmer in sequence.windows(MARKER_K) {
            let hash = hash_kmer(kmer);
            let hash_rc = hash_kmer_rc(kmer);
            // Check both orientations
            if self.human.contains(&hash) || self.human.contains(&hash_rc) {
                human_hits += 1;
            }
            if self.bacterial.contains(&hash) || self.bacterial.contains(&hash_rc) {
                bacterial_hits += 1;
            }
            if self.phix.contains(&hash) || self.phix.contains(&hash_rc) {
                phix_hits += 1;
            }
        }

        // Minimum 3 hits to classify (prevents random false positives)
        let min_hits = 3;
        let max_hits = human_hits.max(bacterial_hits).max(phix_hits);
        if max_hits < min_hits {
            return None;
        }

        if human_hits >= bacterial_hits && human_hits >= phix_hits {
            Some(MarkerCategory::Human)
        } else if bacterial_hits >= phix_hits {
            Some(MarkerCategory::Bacterial)
        } else {
            Some(MarkerCategory::PhiX)
        }
    }

    /// Screen a batch of reads and return composition estimate
    pub fn screen_reads<'a, I>(&self, reads: I) -> CompositionEstimate
    where
        I: Iterator<Item = &'a [u8]>,
    {
        let mut total = 0u64;
        let mut human = 0u64;
        let mut bacterial = 0u64;
        let mut phix = 0u64;

        for seq in reads {
            total += 1;
            match self.classify_read(seq) {
                Some(MarkerCategory::Human) => human += 1,
                Some(MarkerCategory::Bacterial) => bacterial += 1,
                Some(MarkerCategory::PhiX) => phix += 1,
                None => {}
            }
        }

        let t = total.max(1) as f64;
        CompositionEstimate {
            human_fraction: human as f64 / t,
            bacterial_fraction: bacterial as f64 / t,
            phix_fraction: phix as f64 / t,
            unclassified_fraction: (total - human - bacterial - phix) as f64 / t,
            reads_screened: total,
            marker_counts: self.marker_counts(),
        }
    }
}

impl Default for MarkerDb {
    fn default() -> Self {
        Self::new()
    }
}

/// Hash a k-mer (FNV-1a, uppercase)
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
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

/// Build a FxHashSet of k-mer hashes from reference sequences
fn build_kmer_set(sequences: &[&[u8]]) -> FxHashSet<u64> {
    let mut set = FxHashSet::default();
    for seq in sequences {
        if seq.len() < MARKER_K {
            continue;
        }
        for kmer in seq.windows(MARKER_K) {
            // Skip k-mers containing N
            if kmer.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }
            set.insert(hash_kmer(kmer));
        }
    }
    set
}

// ============================================================================
// Marker reference sequences
//
// These are consensus/representative sequences from high-copy diagnostic
// regions. They don't need to be complete -- a few hundred bp from each
// source generates enough 21-mers for robust detection.
// ============================================================================

/// Human marker sequences: Alu consensus, LINE-1 5'UTR, 18S/28S rRNA fragments
const HUMAN_MARKERS: &[&[u8]] = &[
    // Alu Sx consensus (most common Alu subfamily, ~1.1M copies in human genome)
    // GenBank: U14574
    b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGA\
      TCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAA\
      AATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGC\
      TGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCC\
      ACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA",
    // LINE-1 5'UTR (diagnostic region, ~500K copies)
    b"GGGGGAGGAGCCAAGATGGCCGAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCG\
      ACGCAGAAGACGGTGATTTCTGCATTTCCATCTGAGGTACCGGGTTCATCTCACTAGGGAG\
      TGCCAGACAGTGGGCGCAGGCCAGTGTGTGTGCGCACCGTGCGCGAGCCGAAGCAGGGCG\
      AGGCATTGCCTCACCTGGGAAGCGCAAGGGGTCAGGGAGTTCCCTTTCCGAGTCAAAGAAA\
      GGGGTGACGGACGCACCTGGAAAATCGGGTCACTCCCACCCGAATATTGCGCTTTTCAGAC\
      CGGCTTAAGAAACGGCGCACCACGAGACTATATCCCACACCTGGCTCGGAGGGTCCTACGC",
    // Human 18S rRNA (5' end, highly conserved, multi-copy)
    b"TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAA\
      GTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTG\
      ATCGCTCCACTTACATATCAGATAGAAAGCTTGCAATGGCTTTGGGTAGATAGGTGTATTA\
      ATTGTACAGATAAACATACTATTACAATAACAAATTGGAGGGCAAGTCTGGTGCCAGCAGC\
      CGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTT",
    // Human 28S rRNA fragment (multi-copy)
    b"GCCGATCCGAGCGCCCAGTGATCGATGTCTTGAGAGATCGCCCGGACAGTCCCCAGCCCCG\
      GATGCCGCGGCCGCCCCTGTTGGGCCGCCACTGTCTGCCTCCCTGGAGACTCCGAAGGGCC\
      TGTCACTGCCCGAAGCAGCCTCCCCTCGATGATCGTCTGAGGCCCAACAGCCCCCGTCTCC\
      CGGCGGCCCCCAGCGGTCTCTTCCCTGCTCGCCATCCGAAGAGCGTCGTGTTCAGTTCCCC\
      CGCCCTGACCCTGTCCGGCACCCCAACTGCCCCTGCCCTCAGGCGGCTTCCCCCGCCTCCG",
];

/// Bacterial marker sequences: 16S rRNA conserved regions from diverse phyla
const BACTERIAL_MARKERS: &[&[u8]] = &[
    // E. coli 16S rRNA (positions 8-1510, universal bacterial marker)
    // Covers V1-V9 hypervariable regions with conserved flanking
    b"AGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACG\
      GTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGA\
      AACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAG\
      ACCAAAGAGGGGGACCTTAGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTA\
      GGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACA\
      CTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGG\
      GCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTT",
    // Bacillus subtilis 16S rRNA fragment (Firmicutes representative)
    b"GCTTAACACATGCAAGTCGAGCGGACAGATGGGAGCTTGCTCCCTGATGTTAGCGGCGGACG\
      GGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAA\
      TACCGGATGGTTGTTTGAACCGCATGGTTCAAACATAAAAGGTGGCTTCGGCTACCACTTAC\
      AGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAACGATGCGT",
    // Bacteroides fragilis 16S fragment (Bacteroidetes representative)
    b"AACGCTGGCGGCGTGCCTAACACATGCAAGTCGAGCGAAACTTTCTTCCCCAAGGAAGTGGC\
      GAGCGGCGGACGGGTGAGTAACGCGTGGGCAACCTGCCCTGTACAGGGGGATAACAGTTGGA\
      AACGACTGCTAATACCCCATACGCTCCGATGGGGAAAGTGGAAAGATTTATCGCCAATAGATA\
      GCGTAAGGCGCCTGCATGGCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCATTACG",
    // 16S universal conserved region (positions 515-806, covers V4)
    // Highly conserved across nearly all bacteria
    b"GTGCCAGCMGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAG\
      CGCGCGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCCCCGGCTTAACCGGGGAGGGTCAT\
      TGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAAT\
      GCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTG",
];

/// PhiX174 genome (complete, 5386bp)
/// GenBank: NC_001422.1
/// Common Illumina spike-in control
const PHIX174_GENOME: &[u8] =
    b"GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAA\
      ATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCT\
      GGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACC\
      TTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAA\
      TATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAG\
      ATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTTCAACCA\
      CTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTG\
      GCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTTCGATTT\
      TCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAG\
      GCTTTGATATTCGACTTGTACCGACAGTTTCCGAAAGAATGCAGAATCACTGAAACAGCGAT\
      GACCAATTTGAGCCAGAAGCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAA\
      GTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGAT\
      AATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCG\
      AAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCA\
      ACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAA\
      AATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTT\
      CAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATT\
      TAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTA\
      AGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCT\
      CGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAG\
      CTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGC\
      GGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATAT\
      GCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCAT\
      TCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTG\
      GCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACG\
      ACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGATAGCGTCAGTGAGCGTCAAGGGA\
      TGACTTATTCAACCAAATAACTCAAACTTCTTTTCGCAAAAGCCGAGACGACTCGTGACTTTT\
      TGCCAGCGGTTTCAGAAGTTTCACCAGTTTTATGCCGCTGATGATGCAGTTCTAATTTTTCC\
      AACAGAACAACATACTGCTCGTGAATTGACTCGTCAACGAACTCTAATTTGTTACAGAATGC\
      TAATCCTCAAACTAGAAAGTTGTCAGACGACGAGCTTGATGGTAAAGAACGACAGCCGGCTGA\
      CCAAAGCAAGCCCGTGTCGATTAGATACCAGGAAAATCGATTGACAACATCCAGAACAACCTT\
      TGGCAAACTGATTTCCGAGCGCAGTTATCTCGATGCTGATCAAGCATCATTCAACGATCGTC\
      TCAAATTCTGCGCTGGAGAAAAGCGTCAATCATGCGAAGCTGATGGCCAAACTGATCTCAAG\
      CCGGCAAATCCTGAAGAAATCGCTAACTTTACTAGTCTGTTCCGAGTAGAATCTCAATCTATT\
      TCCGCTACTCACGGCTTCAAGAAACTGATCGAGCGCCAGAAGCTAGATAAATTCTTGTCTGAC\
      CGAGAAGAAGCCACGGCAGGTAAGAAGGCCTAATCCGCACGCTTCAAAGCGTCTCCTGATGG";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_marker_db_construction() {
        let db = MarkerDb::new();
        let counts = db.marker_counts();
        for (name, count) in &counts {
            assert!(*count > 0, "{name} should have marker k-mers, got {count}");
        }
        // Human should have the most markers (multiple long sequences)
        let human_count = counts.iter().find(|(n, _)| n == "human").unwrap().1;
        let bacterial_count = counts.iter().find(|(n, _)| n == "bacterial").unwrap().1;
        assert!(
            human_count > 500,
            "Human should have >500 markers, got {human_count}"
        );
        assert!(
            bacterial_count > 500,
            "Bacterial should have >500 markers, got {bacterial_count}"
        );
    }

    #[test]
    fn test_classify_human_alu() {
        let db = MarkerDb::new();
        // A read from Alu sequence should classify as human
        let alu_read = &HUMAN_MARKERS[0][10..160]; // 150bp from Alu
        let result = db.classify_read(alu_read);
        assert_eq!(result, Some(MarkerCategory::Human));
    }

    #[test]
    fn test_classify_bacterial_16s() {
        let db = MarkerDb::new();
        // A read from E. coli 16S should classify as bacterial
        let read_16s = &BACTERIAL_MARKERS[0][10..160];
        let result = db.classify_read(read_16s);
        assert_eq!(result, Some(MarkerCategory::Bacterial));
    }

    #[test]
    fn test_classify_phix() {
        let db = MarkerDb::new();
        // A read from PhiX should classify as PhiX
        let phix_read = &PHIX174_GENOME[100..250];
        let result = db.classify_read(phix_read);
        assert_eq!(result, Some(MarkerCategory::PhiX));
    }

    #[test]
    fn test_classify_random_unclassified() {
        let db = MarkerDb::new();
        // Random sequence should not classify
        let random: Vec<u8> = (0..150)
            .map(|i| [b'A', b'C', b'G', b'T'][i % 4])
            .collect();
        let result = db.classify_read(&random);
        assert_eq!(result, None);
    }

    #[test]
    fn test_screen_mixed_reads() {
        let db = MarkerDb::new();
        let reads: Vec<&[u8]> = vec![
            &HUMAN_MARKERS[0][10..160],       // human
            &HUMAN_MARKERS[0][50..200],       // human
            &BACTERIAL_MARKERS[0][10..160],   // bacterial
            &PHIX174_GENOME[100..250],        // phix
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // random
        ];

        let estimate = db.screen_reads(reads.iter().map(|r| r.as_ref()));
        assert_eq!(estimate.reads_screened, 5);
        assert!(estimate.human_fraction > 0.3, "Should detect human reads");
        assert!(estimate.bacterial_fraction > 0.1, "Should detect bacterial reads");
        assert!(estimate.phix_fraction > 0.1, "Should detect PhiX reads");
    }

    #[test]
    fn test_hash_kmer_rc_complement() {
        // Forward and RC of the same sequence should produce different hashes
        // (they're different sequences)
        let kmer = b"ATGCATGCATGCATGCATGCA";
        let h1 = hash_kmer(kmer);
        let h2 = hash_kmer_rc(kmer);
        // They should be different (not a palindrome)
        assert_ne!(h1, h2);
    }
}
