//! ERV three-signal classifier
//!
//! Classifies clusters of retroviral reads as endogenous or exogenous
//! using three independent signals:
//!   A. ORF integrity (intact = exogenous, frameshifts = endogenous)
//!   B. CpG depletion ratio (low = endogenous, high = exogenous)
//!   C. MinHash distance to ERV vs exogenous reference panels
//!
//! Phase 2 of the ERV analysis module.

use serde::{Deserialize, Serialize};

/// Classification result for a retroviral cluster
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ErvClassification {
    Endogenous,
    Ambiguous,
    Exogenous,
}

/// Detailed classification result with all three signals
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErvClassResult {
    /// Cluster identifier
    pub cluster_id: usize,
    /// Number of reads in this cluster
    pub read_count: usize,
    /// Best matching reference name
    pub best_match: String,

    // Signal A: ORF integrity
    /// Number of intact ORFs found (0-3 for gag/pol/env)
    pub intact_orfs: u8,
    /// ORF score (0.0 = all broken, 1.0 = all intact)
    pub orf_score: f64,

    // Signal B: CpG depletion
    /// CpG observed/expected ratio
    pub cpg_ratio: f64,
    /// CpG score (0.0 = depleted/endogenous, 1.0 = undepleted/exogenous)
    pub cpg_score: f64,

    // Signal C: MinHash distance (placeholder for Phase 2 full implementation)
    /// Distance to nearest ERV consensus
    pub dist_erv: f64,
    /// Name of nearest ERV consensus
    pub nearest_erv: String,
    /// Distance to nearest exogenous reference
    pub dist_exo: f64,
    /// Name of nearest exogenous reference
    pub nearest_exo: String,
    /// MinHash score (0.0 = closer to ERV, 1.0 = closer to exogenous)
    pub minhash_score: f64,

    // Combined
    /// Weighted combined score
    pub combined_score: f64,
    /// Final classification
    pub classification: ErvClassification,
    /// Read depth relative to genomic background
    pub depth_ratio: f64,
}

/// Weights for combining the three signals
/// CpG is dominant for short-read data (ORF and MinHash need assembled contigs).
/// For assembled contigs, ORF would be upweighted.
const WEIGHT_ORF: f64 = 0.15;
const WEIGHT_CPG: f64 = 0.60;
const WEIGHT_MINHASH: f64 = 0.25;

/// Thresholds for classification
const ENDOGENOUS_THRESHOLD: f64 = 0.3;
const EXOGENOUS_THRESHOLD: f64 = 0.7;

/// Minimum ORF lengths (amino acids) for retroviral proteins
const MIN_GAG_AA: usize = 400;
const MIN_POL_AA: usize = 700;
const MIN_ENV_AA: usize = 300;

/// Find the longest ORF in a nucleotide sequence (single frame)
fn longest_orf_in_frame(seq: &[u8], frame: usize) -> usize {
    if seq.len() < frame + 3 {
        return 0;
    }

    let mut max_orf = 0;
    let mut current_orf = 0;

    let mut i = frame;
    while i + 2 < seq.len() {
        let codon = [
            seq[i].to_ascii_uppercase(),
            seq[i + 1].to_ascii_uppercase(),
            seq[i + 2].to_ascii_uppercase(),
        ];

        // Check for stop codons (TAA, TAG, TGA)
        let is_stop = matches!(
            codon,
            [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
        );

        if is_stop {
            max_orf = max_orf.max(current_orf);
            current_orf = 0;
        } else {
            current_orf += 1;
        }
        i += 3;
    }
    max_orf = max_orf.max(current_orf);
    max_orf
}

/// Find the longest ORF across all six reading frames
pub fn longest_orf_six_frame(seq: &[u8]) -> usize {
    let mut max_orf = 0;

    // Forward strand: frames 0, 1, 2
    for frame in 0..3 {
        max_orf = max_orf.max(longest_orf_in_frame(seq, frame));
    }

    // Reverse complement: frames 0, 1, 2
    let rc: Vec<u8> = seq
        .iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect();

    for frame in 0..3 {
        max_orf = max_orf.max(longest_orf_in_frame(&rc, frame));
    }

    max_orf
}

/// Count intact retroviral ORFs in a sequence
/// Returns (count, orf_score)
pub fn count_intact_orfs(seq: &[u8]) -> (u8, f64) {
    let mut intact = 0u8;

    // Check all six frames for long ORFs
    let mut orfs = Vec::new();
    for frame in 0..3 {
        orfs.push(longest_orf_in_frame(seq, frame));
    }
    let rc: Vec<u8> = seq
        .iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect();
    for frame in 0..3 {
        orfs.push(longest_orf_in_frame(&rc, frame));
    }

    // Sort descending -- the three longest ORFs correspond to gag/pol/env
    orfs.sort_unstable_by(|a, b| b.cmp(a));

    // Check if top 3 ORFs meet minimum sizes
    if !orfs.is_empty() && orfs[0] >= MIN_POL_AA {
        intact += 1; // pol (longest)
    }
    if orfs.len() >= 2 && orfs[1] >= MIN_GAG_AA {
        intact += 1; // gag (second longest)
    }
    if orfs.len() >= 3 && orfs[2] >= MIN_ENV_AA {
        intact += 1; // env (third longest)
    }

    let score = match intact {
        3 => 1.0,
        2 => 0.7,
        1 => 0.4,
        _ => 0.0,
    };

    (intact, score)
}

/// Compute CpG score from CpG O/E ratio
pub fn cpg_score(cpg_ratio: f64) -> f64 {
    if cpg_ratio < 0.4 {
        0.0 // strongly endogenous
    } else if cpg_ratio > 0.8 {
        1.0 // likely exogenous
    } else {
        // Linear interpolation in 0.4-0.8 range
        (cpg_ratio - 0.4) / 0.4
    }
}

/// Input for cluster classification
pub struct ClusterInput<'a> {
    pub cluster_id: usize,
    pub consensus_seq: &'a [u8],
    pub read_count: usize,
    pub genomic_depth: f64,
    pub best_match_name: &'a str,
    pub dist_to_erv: f64,
    pub nearest_erv_name: &'a str,
    pub dist_to_exo: f64,
    pub nearest_exo_name: &'a str,
}

/// Classify a retroviral cluster using the three-signal approach
pub fn classify_cluster(input: &ClusterInput) -> ErvClassResult {
    let consensus_seq = input.consensus_seq;
    let read_count = input.read_count;
    // Signal A: ORF integrity
    let (intact_orfs, orf_score) = count_intact_orfs(consensus_seq);

    // Signal B: CpG depletion
    let cpg_ratio_val = super::erv::cpg_ratio(consensus_seq);
    let cpg_sc = cpg_score(cpg_ratio_val);

    // Signal C: MinHash distance
    let minhash_sc = if input.dist_to_erv + input.dist_to_exo > 0.0 {
        input.dist_to_exo / (input.dist_to_erv + input.dist_to_exo)
    } else {
        0.5
    };
    // Invert: closer to exogenous = higher score
    let minhash_sc = 1.0 - minhash_sc;

    // Combined score
    let combined = WEIGHT_ORF * orf_score + WEIGHT_CPG * cpg_sc + WEIGHT_MINHASH * minhash_sc;

    // Depth ratio
    let depth_ratio = if input.genomic_depth > 0.0 {
        read_count as f64 / input.genomic_depth
    } else {
        read_count as f64
    };

    // Classification
    let classification = if combined < ENDOGENOUS_THRESHOLD {
        ErvClassification::Endogenous
    } else if combined > EXOGENOUS_THRESHOLD {
        ErvClassification::Exogenous
    } else {
        ErvClassification::Ambiguous
    };

    ErvClassResult {
        cluster_id: input.cluster_id,
        read_count,
        best_match: input.best_match_name.to_string(),
        intact_orfs,
        orf_score,
        cpg_ratio: cpg_ratio_val,
        cpg_score: cpg_sc,
        dist_erv: input.dist_to_erv,
        nearest_erv: input.nearest_erv_name.to_string(),
        dist_exo: input.dist_to_exo,
        nearest_exo: input.nearest_exo_name.to_string(),
        minhash_score: minhash_sc,
        combined_score: combined,
        classification,
        depth_ratio,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_longest_orf_simple() {
        // A sequence with one long ORF (no stop codons)
        let seq = b"ATGAAACCCGGGTTTATGAAACCCGGGTTTATGAAACCCGGGTTT";
        let orf = longest_orf_in_frame(seq, 0);
        // 45 bases / 3 = 15 codons, no stops = 15 aa ORF
        assert_eq!(orf, 15);
    }

    #[test]
    fn test_longest_orf_with_stop() {
        // TAA at positions 9-11 splits into two ORFs
        let seq = b"ATGAAACCCTAAATGAAACCCGGGTTT";
        let orf = longest_orf_in_frame(seq, 0);
        // First ORF: 3 codons (ATG AAA CCC), second: 5 codons (ATG AAA CCC GGG TTT)
        assert_eq!(orf, 5);
    }

    #[test]
    fn test_count_intact_orfs_degraded() {
        // Short sequence with many stop codons -- no intact ORFs
        let seq = b"ATGTAAATGTAGATGTGA";
        let (intact, score) = count_intact_orfs(seq);
        assert_eq!(intact, 0);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_cpg_score_depleted() {
        assert_eq!(cpg_score(0.2), 0.0);
        assert_eq!(cpg_score(0.35), 0.0);
    }

    #[test]
    fn test_cpg_score_exogenous() {
        assert_eq!(cpg_score(0.9), 1.0);
        assert_eq!(cpg_score(1.0), 1.0);
    }

    #[test]
    fn test_cpg_score_intermediate() {
        let score = cpg_score(0.6);
        assert!(score > 0.4 && score < 0.6, "Expected ~0.5, got {score}");
    }

    #[test]
    fn test_classify_endogenous() {
        let result = classify_cluster(&ClusterInput {
            cluster_id: 1,
            consensus_seq: b"ATGTAAATGTAGATGTGAATGTAAATGTAGATGTGA",
            read_count: 45,
            genomic_depth: 30.0,
            best_match_name: "HERV-K_like",
            dist_to_erv: 0.12,
            nearest_erv_name: "HERV-K_con",
            dist_to_exo: 0.35,
            nearest_exo_name: "MMTV",
        });
        assert!(matches!(result.classification, ErvClassification::Endogenous));
    }

    #[test]
    fn test_classify_exogenous() {
        let result = classify_cluster(&ClusterInput {
            cluster_id: 2,
            consensus_seq: b"ACGTCGATCGCGATCGACGTACGTCGATCGCGATCG",
            read_count: 892,
            genomic_depth: 30.0,
            best_match_name: "HHV-6B",
            dist_to_erv: 0.35,
            nearest_erv_name: "ciHHV-6A",
            dist_to_exo: 0.01,
            nearest_exo_name: "HHV-6B",
        });
        assert!(
            result.combined_score > ENDOGENOUS_THRESHOLD,
            "Should not be classified as endogenous: score={}",
            result.combined_score
        );
    }
}
