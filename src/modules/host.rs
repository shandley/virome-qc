//! Host depletion module -- Super Bloom k-mer containment screening
//!
//! Uses a pre-built Super Bloom filter of the host genome (e.g., T2T-CHM13)
//! to classify reads by k-mer containment fraction:
//!
//! - Host (>50% containment): remove
//! - Ambiguous (15-50%): flag for review
//! - Not host (<15%): keep
//!
//! The Super Bloom filter is built once via `virome-qc db build` and loaded
//! at pipeline start. ~4 GiB for the human genome, near-zero false positives
//! with the findere consistency scheme.

use crate::config::HostConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use superbloom::{FrozenSuperBloom, MinimizerMode, SuperBloom, SuperBloomConfig};

/// Host depletion module using Super Bloom k-mer containment
pub struct HostFilter {
    /// Frozen (read-only) Super Bloom filter
    filter: FrozenSuperBloom,
    /// Minimum containment fraction to classify as host
    host_threshold: f64,
    /// Minimum containment fraction to flag as ambiguous
    ambiguous_threshold: f64,
    /// Standard stats
    stats: AtomicStats,
    /// Host reads removed
    host_removed: AtomicU64,
    /// Ambiguous reads flagged
    ambiguous_flagged: AtomicU64,
    /// Containment distribution histogram: 10 bins (0-10%, 10-20%, ..., 90-100%)
    /// Only counts reads with containment > 0 (skips the vast majority at 0%)
    containment_hist: [AtomicU64; 10],
    /// Reads with zero containment
    zero_containment: AtomicU64,
}

impl HostFilter {
    /// Create a host filter from a pre-built Super Bloom filter file
    pub fn from_filter_file(config: &HostConfig) -> Result<Self, anyhow::Error> {
        let filter_path = resolve_filter_path(&config.reference)?;
        let sb = FrozenSuperBloom::load(&filter_path).map_err(|e| {
            anyhow::anyhow!(
                "Failed to load host filter from {}: {}",
                filter_path.display(),
                e
            )
        })?;

        Ok(Self {
            filter: sb,
            host_threshold: config.host_threshold,
            ambiguous_threshold: config.ambiguous_threshold,
            stats: AtomicStats::new(),
            host_removed: AtomicU64::new(0),
            ambiguous_flagged: AtomicU64::new(0),
            containment_hist: std::array::from_fn(|_| AtomicU64::new(0)),
            zero_containment: AtomicU64::new(0),
        })
    }

    /// Build a Super Bloom filter from a FASTA reference and save to disk
    ///
    /// Auto-sizes the filter based on reference size:
    /// - <100 Mbp: 128 MiB (exponent 30)
    /// - <500 Mbp: 512 MiB (exponent 32)
    /// - <2 Gbp: 2 GiB (exponent 34)
    /// - >=2 Gbp: 4 GiB (exponent 35)
    ///
    /// Target: ~11 bits per k-mer for good accuracy with 8 hash functions.
    pub fn build_filter(fasta_path: &Path, output_path: &Path) -> Result<(), anyhow::Error> {
        // Estimate reference size to auto-size filter
        let file_size = std::fs::metadata(fasta_path)?.len();
        // FASTA overhead ~2% (headers + newlines), rough estimate of sequence bp
        let est_bp = (file_size as f64 * 0.98) as u64;

        let exponent = if est_bp < 100_000_000 {
            30 // 128 MiB
        } else if est_bp < 500_000_000 {
            32 // 512 MiB
        } else if est_bp < 2_000_000_000 {
            34 // 2 GiB
        } else {
            35 // 4 GiB
        };

        let filter_size_mb = (1u64 << exponent) / 8 / 1024 / 1024;

        let config = SuperBloomConfig {
            k: 31,
            m: 21,
            s: 27,
            n_hashes: 8,
            bit_vector_size_exponent: exponent,
            block_size_exponent: 9,
            minimizer_mode: MinimizerMode::Simd,
        };

        eprintln!(
            "Building Super Bloom host filter ({} MiB, k=31, est. {} Mbp reference)...",
            filter_size_mb,
            est_bp / 1_000_000
        );
        let mut sb = SuperBloom::new(config)
            .map_err(|e| anyhow::anyhow!("Failed to create SuperBloom: {}", e))?;

        // Load and index FASTA sequences
        let contents = std::fs::read_to_string(fasta_path)?;
        let mut current_seq = Vec::new();
        let mut seq_count = 0;
        let mut total_kmers = 0u64;

        for line in contents.lines() {
            if line.starts_with('>') {
                if !current_seq.is_empty() {
                    match sb.add_sequence(&current_seq) {
                        Ok(added) => total_kmers += added,
                        Err(e) => eprintln!("  Warning: {}", e),
                    }
                    seq_count += 1;
                    if seq_count % 5 == 0 {
                        eprint!("\r  {} sequences, {} k-mers...", seq_count, total_kmers);
                    }
                    current_seq.clear();
                }
            } else {
                current_seq.extend_from_slice(line.trim().to_uppercase().as_bytes());
            }
        }
        if !current_seq.is_empty() {
            match sb.add_sequence(&current_seq) {
                Ok(added) => total_kmers += added,
                Err(e) => eprintln!("  Warning: {}", e),
            }
            seq_count += 1;
        }
        eprintln!(
            "\r  {} sequences, {} k-mers indexed",
            seq_count, total_kmers
        );

        // Save
        let frozen = sb.into_frozen();
        frozen
            .save(output_path)
            .map_err(|e| anyhow::anyhow!("Failed to save filter: {}", e))?;

        eprintln!("Filter saved to: {}", output_path.display());
        Ok(())
    }

    /// Compute k-mer containment fraction for a read
    fn containment(&self, sequence: &[u8]) -> f64 {
        let hits = self.filter.query_sequence(sequence);
        if hits.is_empty() {
            return 0.0;
        }
        let n_hits = hits.iter().filter(|&&h| h).count();
        n_hits as f64 / hits.len() as f64
    }
}

impl QcModule for HostFilter {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.len() < 31 {
            return; // too short for k=31
        }

        let containment = self.containment(seq);

        // Record containment distribution
        if containment > 0.0 {
            let bin = ((containment * 10.0) as usize).min(9);
            self.containment_hist[bin].fetch_add(1, Ordering::Relaxed);
        } else {
            self.zero_containment.fetch_add(1, Ordering::Relaxed);
        }

        if containment >= self.host_threshold {
            record.fail("host");
            self.host_removed.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        } else if containment >= self.ambiguous_threshold {
            record.flag("host_ambiguous");
            self.ambiguous_flagged.fetch_add(1, Ordering::Relaxed);
        }
    }

    fn report(&self) -> ModuleReport {
        let hist: Vec<u64> = self
            .containment_hist
            .iter()
            .map(|c| c.load(Ordering::Relaxed))
            .collect();

        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "host_removed": self.host_removed.load(Ordering::Relaxed),
                "ambiguous_flagged": self.ambiguous_flagged.load(Ordering::Relaxed),
                "host_threshold": self.host_threshold,
                "ambiguous_threshold": self.ambiguous_threshold,
                "containment_distribution": {
                    "bins": ["0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                             "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"],
                    "counts": hist,
                    "zero_containment": self.zero_containment.load(Ordering::Relaxed),
                },
            }),
        )
    }

    fn name(&self) -> &str {
        "host"
    }
}

/// Resolve a filter path from a reference name or direct path
fn resolve_filter_path(reference: &str) -> Result<std::path::PathBuf, anyhow::Error> {
    let path = std::path::PathBuf::from(reference);

    // Direct path to .sbf file
    if path.extension().map(|e| e == "sbf").unwrap_or(false) && path.exists() {
        return Ok(path);
    }

    // Check standard locations
    let home = std::env::var("HOME").unwrap_or_default();
    let mut standard_locations = vec![
        // Relative to current directory
        format!("{}.sbf", reference),
        // In a databases subdirectory
        format!("databases/host/{}.sbf", reference),
        // Common naming variants
        format!("benchmark_data/references/{}.sbf", reference),
        format!("benchmark_data/references/{}_t2t.sbf", reference),
        // Home directory
        format!("{}/.virome-qc/host/{}.sbf", home, reference),
        format!("{}/.virome-qc/host/{}_t2t.sbf", home, reference),
    ];

    // Also check VIROME_QC_DB environment variable
    if let Ok(db_dir) = std::env::var("VIROME_QC_DB") {
        standard_locations.push(format!("{}/{}.sbf", db_dir, reference));
        standard_locations.push(format!("{}/{}_t2t.sbf", db_dir, reference));
    }

    for loc in &standard_locations {
        let p = std::path::PathBuf::from(loc);
        if p.exists() {
            return Ok(p);
        }
    }

    anyhow::bail!(
        "Host filter not found for '{}'. Build it with: virome-qc db build --host {}",
        reference,
        reference
    )
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
    fn test_containment_calculation() {
        // Build a tiny filter with a known sequence
        let config = SuperBloomConfig {
            k: 15,
            m: 10,
            s: 12,
            n_hashes: 4,
            bit_vector_size_exponent: 20, // 1 MiB
            block_size_exponent: 9,
            minimizer_mode: MinimizerMode::Simd,
        };

        let mut sb = SuperBloom::new(config).unwrap();
        let ref_seq =
            b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
        sb.add_sequence(ref_seq).unwrap();
        let frozen = sb.into_frozen();

        let filter = HostFilter {
            filter: frozen,
            host_threshold: 0.50,
            ambiguous_threshold: 0.15,
            stats: AtomicStats::new(),
            host_removed: AtomicU64::new(0),
            ambiguous_flagged: AtomicU64::new(0),
            containment_hist: std::array::from_fn(|_| AtomicU64::new(0)),
            zero_containment: AtomicU64::new(0),
        };

        // A read from the reference should have high containment
        let containment = filter.containment(&ref_seq[10..60]);
        assert!(
            containment > 0.5,
            "Reference-derived read should have high containment, got {:.3}",
            containment
        );

        // A random read should have low containment
        let random: Vec<u8> = (0..50)
            .map(|i| [b'T', b'G', b'C', b'A'][(i * 7 + 3) % 4])
            .collect();
        let random_containment = filter.containment(&random);
        assert!(
            random_containment < 0.3,
            "Random read should have low containment, got {:.3}",
            random_containment
        );
    }

    #[test]
    fn test_host_classification() {
        let config = SuperBloomConfig {
            k: 15,
            m: 10,
            s: 12,
            n_hashes: 4,
            bit_vector_size_exponent: 20,
            block_size_exponent: 9,
            minimizer_mode: MinimizerMode::Simd,
        };

        let mut sb = SuperBloom::new(config).unwrap();
        let ref_seq = b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
        sb.add_sequence(ref_seq).unwrap();
        let frozen = sb.into_frozen();

        let filter = HostFilter {
            filter: frozen,
            host_threshold: 0.50,
            ambiguous_threshold: 0.15,
            stats: AtomicStats::new(),
            host_removed: AtomicU64::new(0),
            ambiguous_flagged: AtomicU64::new(0),
            containment_hist: std::array::from_fn(|_| AtomicU64::new(0)),
            zero_containment: AtomicU64::new(0),
        };

        // Host read
        let mut host_record = make_record(&ref_seq[5..55]);
        filter.process(&mut host_record);
        assert!(
            host_record.is_failed(),
            "Reference read should be classified as host"
        );

        // Random read
        let random: Vec<u8> = (0..50)
            .map(|i| [b'T', b'G', b'C', b'A'][(i * 7 + 3) % 4])
            .collect();
        let mut random_record = make_record(&random);
        filter.process(&mut random_record);
        assert!(!random_record.is_failed(), "Random read should pass");
    }
}
