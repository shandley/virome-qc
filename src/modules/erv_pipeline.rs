//! ERV analysis pipeline -- post-QC retroviral classification
//!
//! Connects the ERV screener (Phase 1) to the classifier (Phase 2):
//! 1. Collect reads flagged with retroviral signal
//! 2. Cluster by shared k-mer content
//! 3. Build consensus per cluster
//! 4. Classify each cluster using three-signal approach
//! 5. Output results for passport

use crate::modules::erv_classifier::{classify_cluster, ClusterInput, ErvClassResult};
use biometal::operations::sketching::{MinHashSketch, SketchConfig};
use rustc_hash::FxHashMap;
use std::path::Path;

/// Configuration for ERV analysis
pub struct ErvAnalysisConfig {
    /// MinHash sketch size for reference comparison
    pub sketch_size: usize,
    /// K-mer size for MinHash sketching
    pub sketch_k: usize,
    /// Minimum reads to form a cluster
    pub min_cluster_size: usize,
    /// Average genomic depth (for depth ratio calculation)
    pub genomic_depth: f64,
}

impl Default for ErvAnalysisConfig {
    fn default() -> Self {
        Self {
            sketch_size: 1000,
            sketch_k: 15,
            min_cluster_size: 3,
            genomic_depth: 30.0,
        }
    }
}

/// Pre-computed reference panel for MinHash comparison
pub struct ReferencePanel {
    pub sketches: Vec<(String, MinHashSketch)>,
}

impl ReferencePanel {
    /// Build reference panel from FASTA sequences
    pub fn from_sequences(sequences: &[(&str, &[u8])], config: &SketchConfig) -> Self {
        let sketches = sequences
            .iter()
            .filter_map(|(name, seq)| {
                // Filter to ACGT only (biometal's ntHash panics on ambiguous bases)
                let clean: Vec<u8> = seq
                    .iter()
                    .filter(|&&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
                    .copied()
                    .collect();
                if clean.len() < config.k * 2 {
                    return None; // Too short after cleaning
                }
                let mut sketch = MinHashSketch::new(config);
                sketch.add_sequence(&clean);
                Some((name.to_string(), sketch))
            })
            .collect();
        Self { sketches }
    }

    /// Find nearest reference and its distance
    pub fn nearest(&self, query: &MinHashSketch) -> (String, f64) {
        let mut best_name = "unknown".to_string();
        let mut best_dist = 1.0;

        for (name, ref_sketch) in &self.sketches {
            if ref_sketch.is_empty() || query.is_empty() {
                continue;
            }
            let jaccard = query.jaccard(ref_sketch);
            let dist = 1.0 - jaccard;
            if dist < best_dist {
                best_dist = dist;
                best_name = name.clone();
            }
        }

        (best_name, best_dist)
    }
}

/// A cluster of retroviral reads
pub struct ReadCluster {
    pub reads: Vec<Vec<u8>>,
    pub consensus: Vec<u8>,
}

impl ReadCluster {
    /// Build consensus by majority vote at each position
    fn build_consensus(reads: &[Vec<u8>]) -> Vec<u8> {
        if reads.is_empty() {
            return Vec::new();
        }

        let max_len = reads.iter().map(|r| r.len()).max().unwrap_or(0);
        let mut consensus = Vec::with_capacity(max_len);

        for pos in 0..max_len {
            let mut counts = [0u32; 5]; // A, C, G, T, N
            for read in reads {
                if pos < read.len() {
                    let idx = match read[pos].to_ascii_uppercase() {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => 4,
                    };
                    counts[idx] += 1;
                }
            }
            let best = counts
                .iter()
                .enumerate()
                .max_by_key(|(_, &c)| c)
                .map(|(i, _)| i)
                .unwrap_or(4);
            consensus.push([b'A', b'C', b'G', b'T', b'N'][best]);
        }

        consensus
    }
}

/// Simple k-mer clustering: group reads by shared k-mer content
///
/// Uses a union-find approach: reads sharing any k-mer are merged into
/// the same cluster. Fast but aggressive (may over-merge).
pub fn cluster_reads(reads: &[Vec<u8>], k: usize) -> Vec<ReadCluster> {
    if reads.is_empty() {
        return Vec::new();
    }

    // Build k-mer -> read indices mapping
    let mut kmer_to_reads: FxHashMap<u64, Vec<usize>> = FxHashMap::default();

    for (idx, read) in reads.iter().enumerate() {
        if read.len() < k {
            continue;
        }
        for kmer in read.windows(k) {
            if kmer.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }
            let hash = hash_kmer(kmer);
            kmer_to_reads.entry(hash).or_default().push(idx);
        }
    }

    // Union-find
    let n = reads.len();
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], i: usize) -> usize {
        if parent[i] != i {
            parent[i] = find(parent, parent[i]);
        }
        parent[i]
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }

    // Count shared k-mers between read pairs
    let mut pair_counts: FxHashMap<(usize, usize), u32> = FxHashMap::default();
    for indices in kmer_to_reads.values() {
        if indices.len() > 1 && indices.len() < 50 {
            // Skip very common k-mers (>50 reads = conserved domain)
            for i in 0..indices.len() {
                for j in (i + 1)..indices.len() {
                    let key = (indices[i].min(indices[j]), indices[i].max(indices[j]));
                    *pair_counts.entry(key).or_insert(0) += 1;
                }
            }
        }
    }

    // Only merge reads sharing at least MIN_SHARED_KMERS k-mers
    const MIN_SHARED_KMERS: u32 = 5;
    for ((a, b), count) in &pair_counts {
        if *count >= MIN_SHARED_KMERS {
            union(&mut parent, *a, *b);
        }
    }

    // Collect clusters
    let mut clusters: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..n {
        let root = find(&mut parent, i);
        clusters.entry(root).or_default().push(i);
    }

    // Build ReadCluster structs
    clusters
        .into_values()
        .map(|indices| {
            let cluster_reads: Vec<Vec<u8>> =
                indices.iter().map(|&i| reads[i].clone()).collect();
            let consensus = ReadCluster::build_consensus(&cluster_reads);
            ReadCluster {
                reads: cluster_reads,
                consensus,
            }
        })
        .collect()
}

/// Run the complete ERV analysis pipeline
pub fn analyze_erv_reads(
    retroviral_reads: &[Vec<u8>],
    erv_panel: &ReferencePanel,
    exo_panel: &ReferencePanel,
    config: &ErvAnalysisConfig,
) -> Vec<ErvClassResult> {
    if retroviral_reads.is_empty() {
        return Vec::new();
    }

    let sketch_config = SketchConfig::new()
        .with_k(config.sketch_k)
        .with_sketch_size(config.sketch_size);

    // Step 1: Cluster reads
    let clusters = cluster_reads(retroviral_reads, config.sketch_k);

    // Step 2: Classify each cluster
    let mut results = Vec::new();

    for (i, cluster) in clusters.iter().enumerate() {
        if cluster.reads.len() < config.min_cluster_size {
            continue;
        }

        // Compute MinHash sketch for cluster consensus (clean ACGT only)
        let clean_consensus: Vec<u8> = cluster
            .consensus
            .iter()
            .filter(|&&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .copied()
            .collect();

        if clean_consensus.len() < config.sketch_k * 2 {
            continue;
        }

        let mut query_sketch = MinHashSketch::new(&sketch_config);
        query_sketch.add_sequence(&clean_consensus);

        // Find nearest in both panels
        let (nearest_erv, dist_erv) = erv_panel.nearest(&query_sketch);
        let (nearest_exo, dist_exo) = exo_panel.nearest(&query_sketch);

        // Determine best match name (from whichever panel is closer)
        let best_match = if dist_erv < dist_exo {
            &nearest_erv
        } else {
            &nearest_exo
        };

        // Classify
        let result = classify_cluster(&ClusterInput {
            cluster_id: i,
            consensus_seq: &cluster.consensus,
            read_count: cluster.reads.len(),
            genomic_depth: config.genomic_depth,
            best_match_name: best_match,
            dist_to_erv: dist_erv,
            nearest_erv_name: &nearest_erv,
            dist_to_exo: dist_exo,
            nearest_exo_name: &nearest_exo,
        });

        results.push(result);
    }

    results
}

/// Load ERV reference panel from Dfam HERV consensus FASTA
pub fn load_erv_panel(fasta_path: &Path, sketch_config: &SketchConfig) -> ReferencePanel {
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();

    if let Ok(contents) = std::fs::read_to_string(fasta_path) {
        for line in contents.lines() {
            if line.starts_with('>') {
                if !current_seq.is_empty() {
                    sequences.push((current_name.clone(), current_seq.clone()));
                    current_seq.clear();
                }
                current_name = line.strip_prefix('>').unwrap_or(line).split_whitespace().next().unwrap_or("").to_string();
            } else {
                current_seq.extend_from_slice(line.trim().to_uppercase().as_bytes());
            }
        }
        if !current_seq.is_empty() {
            sequences.push((current_name, current_seq));
        }
    }

    let seq_refs: Vec<(&str, &[u8])> = sequences
        .iter()
        .map(|(name, seq)| (name.as_str(), seq.as_slice()))
        .collect();

    ReferencePanel::from_sequences(&seq_refs, sketch_config)
}

/// Build exogenous reference panel from embedded retroviral sequences
pub fn build_exo_panel(sketch_config: &SketchConfig) -> ReferencePanel {
    use super::erv::erv_sequences::RETROVIRAL_SEQUENCES;

    // Use simple numbered names (full names are in the source file comments)
    let sequences: Vec<(&str, &[u8])> = RETROVIRAL_SEQUENCES
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            // Generate a short name from the sequence
            let name: &str = match i {
                0..=18 => "Alpharetrovirus",
                19..=22 => "Betaretrovirus",
                23..=25 => "Deltaretrovirus",
                26 => "Epsilonretrovirus",
                27..=34 => "Gammaretrovirus",
                35..=43 => "Lentivirus",
                44 => "Spumavirus",
                45 => "HHV-6A",
                46 => "HHV-6B",
                _ => "Unknown",
            };
            (name, *seq)
        })
        .collect();

    ReferencePanel::from_sequences(&sequences, sketch_config)
}

/// FNV-1a hash for k-mer clustering
#[inline]
fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_identical_reads() {
        let reads = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),
        ];
        let clusters = cluster_reads(&reads, 15);
        // First two reads share k-mers, third is different
        assert!(clusters.len() >= 2, "Should have at least 2 clusters");
    }

    #[test]
    fn test_consensus_building() {
        let reads = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"ACGTACGT".to_vec(),
        ];
        let consensus = ReadCluster::build_consensus(&reads);
        assert_eq!(consensus, b"ACGTACGT");
    }

    #[test]
    fn test_reference_panel() {
        let config = SketchConfig::new().with_k(15).with_sketch_size(100);
        // Sequences must be longer than 2*k for sketching
        let seq1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq2 = b"TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAA";
        let sequences = vec![
            ("seq1", seq1.as_slice()),
            ("seq2", seq2.as_slice()),
        ];
        let panel = ReferencePanel::from_sequences(&sequences, &config);
        assert_eq!(panel.sketches.len(), 2);

        // Query matching seq1
        let mut query = MinHashSketch::new(&config);
        query.add_sequence(seq1);
        let (name, dist) = panel.nearest(&query);
        assert_eq!(name, "seq1");
        assert!(dist < 0.1, "Should be very close to seq1: {dist}");
    }

    #[test]
    fn test_build_exo_panel() {
        // Use smaller k to avoid issues with very short sequences
        let config = SketchConfig::new().with_k(15).with_sketch_size(100);
        let panel = build_exo_panel(&config);
        assert!(
            panel.sketches.len() > 40,
            "Should have 46+ sketches, got {}",
            panel.sketches.len()
        );
    }
}
