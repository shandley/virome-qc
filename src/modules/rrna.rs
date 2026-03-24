//! rRNA screening module -- SILVA-based k-mer containment filter
//!
//! Uses a pre-built sorted hash vector of 21-mers from the SILVA rRNA
//! database (SSU + LSU). Reads with high k-mer containment against the
//! rRNA database are classified as rRNA and removed.
//!
//! Build the filter with: `virome-qc db --rrna`

use crate::config::RrnaConfig;
use crate::modules::{AtomicStats, ModuleReport, QcModule};
use crate::pipeline::AnnotatedRecord;
use anyhow::Result;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

/// K-mer size for rRNA screening (matches contaminant screener)
const RRNA_K: usize = 21;

/// Magic bytes for the rRNA filter file format
const RRNA_FILTER_MAGIC: &[u8; 8] = b"VRQCRRNA";

/// rRNA screening module using sorted hash vector
pub struct RrnaFilter {
    /// Sorted vector of FNV-1a hashes of all 21-mers from SILVA
    hashes: Vec<u64>,
    /// Minimum fraction of read k-mers matching to classify as rRNA
    min_kmer_fraction: f64,
    /// Stats
    stats: AtomicStats,
    rrna_removed: AtomicU64,
}

impl RrnaFilter {
    /// Load a pre-built rRNA filter from disk
    pub fn from_filter_file(config: &RrnaConfig) -> Result<Self> {
        let path = resolve_rrna_filter_path(&config.filter)?;
        let hashes = Self::load_hashes(&path)?;
        eprintln!(
            "  rRNA filter loaded: {} k-mers from {}",
            hashes.len(),
            path.display()
        );

        Ok(Self {
            hashes,
            min_kmer_fraction: config.min_kmer_fraction,
            stats: AtomicStats::new(),
            rrna_removed: AtomicU64::new(0),
        })
    }

    /// Load sorted hash vector from binary file
    fn load_hashes(path: &Path) -> Result<Vec<u64>> {
        let file = std::fs::File::open(path)?;
        let file_len = file.metadata()?.len();
        let mut reader = BufReader::with_capacity(8 * 1024 * 1024, file);

        // Read and verify magic bytes
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != RRNA_FILTER_MAGIC {
            anyhow::bail!(
                "Invalid rRNA filter file (wrong magic bytes). Rebuild with: virome-qc db --rrna"
            );
        }

        // Read hash count
        let mut count_bytes = [0u8; 8];
        reader.read_exact(&mut count_bytes)?;
        let count = u64::from_le_bytes(count_bytes) as usize;

        // Verify file size
        let expected_size = 8 + 8 + (count as u64 * 8);
        if file_len != expected_size {
            anyhow::bail!(
                "rRNA filter file is corrupt (expected {} bytes, got {})",
                expected_size,
                file_len
            );
        }

        // Read all hashes
        let mut hashes = vec![0u64; count];
        let hash_bytes = unsafe {
            std::slice::from_raw_parts_mut(hashes.as_mut_ptr() as *mut u8, count * 8)
        };
        reader.read_exact(hash_bytes)?;

        // Verify sorted (spot check)
        if count > 1 && hashes[0] > hashes[1] {
            anyhow::bail!("rRNA filter file is corrupt (not sorted)");
        }

        Ok(hashes)
    }

    /// Build rRNA filter from FASTA file(s)
    ///
    /// Extracts all 21-mers from both strands, hashes with FNV-1a,
    /// sorts and deduplicates, writes binary filter file.
    pub fn build_filter(fasta_paths: &[&Path], output_path: &Path) -> Result<()> {
        eprintln!("Building rRNA k-mer filter (k={})...", RRNA_K);

        let mut hashes: Vec<u64> = Vec::with_capacity(200_000_000);
        let mut seq_count = 0u64;

        for fasta_path in fasta_paths {
            let is_gzipped = fasta_path
                .extension()
                .map(|e| e == "gz")
                .unwrap_or(false);

            let file = std::fs::File::open(fasta_path)?;

            if is_gzipped {
                let decoder = flate2::read::GzDecoder::new(file);
                let reader = BufReader::with_capacity(4 * 1024 * 1024, decoder);
                seq_count += Self::process_fasta_reader(reader, &mut hashes)?;
            } else {
                let reader = BufReader::with_capacity(4 * 1024 * 1024, file);
                seq_count += Self::process_fasta_reader(reader, &mut hashes)?;
            }
        }

        eprintln!(
            "  {} sequences processed, {} raw k-mer hashes",
            seq_count,
            hashes.len()
        );

        // Sort and deduplicate
        eprint!("  Sorting and deduplicating...");
        hashes.sort_unstable();
        hashes.dedup();
        eprintln!(" {} unique hashes", hashes.len());

        // Write binary filter file
        eprint!("  Writing filter...");
        let file = std::fs::File::create(output_path)?;
        let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file);

        // Magic bytes
        writer.write_all(RRNA_FILTER_MAGIC)?;

        // Hash count
        writer.write_all(&(hashes.len() as u64).to_le_bytes())?;

        // Hash data
        let hash_bytes = unsafe {
            std::slice::from_raw_parts(hashes.as_ptr() as *const u8, hashes.len() * 8)
        };
        writer.write_all(hash_bytes)?;
        writer.flush()?;

        let file_size_mb = (8 + 8 + hashes.len() * 8) as f64 / 1e6;
        eprintln!(" {:.0} MB", file_size_mb);
        eprintln!("rRNA filter saved to: {}", output_path.display());

        Ok(())
    }

    /// Process a FASTA reader, extracting k-mer hashes
    fn process_fasta_reader<R: Read>(reader: R, hashes: &mut Vec<u64>) -> Result<u64> {
        use std::io::BufRead;
        let buf_reader = std::io::BufReader::new(reader);
        let mut current_seq = Vec::with_capacity(4000);
        let mut seq_count = 0u64;

        for line in buf_reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                if !current_seq.is_empty() {
                    Self::extract_kmers(&current_seq, hashes);
                    seq_count += 1;
                    if seq_count % 100_000 == 0 {
                        eprint!("\r  {} sequences, {} k-mers...", seq_count, hashes.len());
                    }
                    current_seq.clear();
                }
            } else {
                // Convert U->T (SILVA uses RNA alphabet), uppercase
                for &b in line.trim().as_bytes() {
                    let base = match b {
                        b'U' | b'u' => b'T',
                        _ => b.to_ascii_uppercase(),
                    };
                    current_seq.push(base);
                }
            }
        }
        // Last sequence
        if !current_seq.is_empty() {
            Self::extract_kmers(&current_seq, hashes);
            seq_count += 1;
        }

        Ok(seq_count)
    }

    /// Extract all k-mer hashes from a sequence (both strands)
    fn extract_kmers(seq: &[u8], hashes: &mut Vec<u64>) {
        if seq.len() < RRNA_K {
            return;
        }
        for kmer in seq.windows(RRNA_K) {
            if kmer.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }
            hashes.push(hash_kmer(kmer));
            hashes.push(hash_kmer_rc(kmer));
        }
    }

    /// Check if a hash is present in the filter
    #[inline]
    fn contains(&self, hash: u64) -> bool {
        self.hashes.binary_search(&hash).is_ok()
    }

    /// Compute k-mer containment of a read against the rRNA database
    fn containment(&self, sequence: &[u8]) -> f64 {
        if sequence.len() < RRNA_K || self.hashes.is_empty() {
            return 0.0;
        }
        let total = sequence.len() - RRNA_K + 1;
        let mut hits = 0u32;
        for kmer in sequence.windows(RRNA_K) {
            if self.contains(hash_kmer(kmer)) {
                hits += 1;
            }
        }
        hits as f64 / total as f64
    }
}

impl QcModule for RrnaFilter {
    fn process(&self, record: &mut AnnotatedRecord) {
        self.stats.record_processed();

        let seq = &record.record.sequence;
        if seq.len() < RRNA_K {
            return;
        }

        let containment = self.containment(seq);
        if containment >= self.min_kmer_fraction {
            record.fail("rrna");
            self.rrna_removed.fetch_add(1, Ordering::Relaxed);
            self.stats.record_removed();
        }
    }

    fn report(&self) -> ModuleReport {
        let removed = self.rrna_removed.load(Ordering::Relaxed);
        self.stats.to_report(
            self.name(),
            serde_json::json!({
                "rrna_removed": removed,
                "filter_size": self.hashes.len(),
            }),
        )
    }

    fn name(&self) -> &str {
        "rrna"
    }
}

/// FNV-1a hash (uppercase, same as contaminant screener)
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

/// Resolve rRNA filter file path, checking standard locations
fn resolve_rrna_filter_path(reference: &str) -> Result<std::path::PathBuf> {
    let path = std::path::PathBuf::from(reference);
    if path.exists() {
        return Ok(path);
    }

    // Check standard locations
    let mut candidates = vec![
        std::path::PathBuf::from("rrna.rrf"),
        std::path::PathBuf::from("databases/rrna/silva.rrf"),
        std::path::PathBuf::from("benchmark_data/references/rrna_silva.rrf"),
        std::env::var("HOME")
            .map(std::path::PathBuf::from)
            .unwrap_or_default()
            .join(".virome-qc/rrna/silva.rrf"),
    ];

    // Also check VIROME_QC_DB environment variable
    if let Ok(db_dir) = std::env::var("VIROME_QC_DB") {
        candidates.push(std::path::PathBuf::from(format!("{}/rrna_silva.rrf", db_dir)));
        candidates.push(std::path::PathBuf::from(format!("{}/silva.rrf", db_dir)));
    }

    for candidate in &candidates {
        if candidate.exists() {
            return Ok(candidate.clone());
        }
    }

    anyhow::bail!(
        "rRNA filter not found at '{}'. Build one with:\n  \
         virome-qc db --rrna\n\
         Or download SILVA and build manually:\n  \
         virome-qc db --rrna /path/to/SILVA_SSU_NR99.fasta.gz -o rrna.rrf",
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
    fn test_hash_deterministic() {
        let kmer = b"AGATCGGAAGAGCACACGTCT";
        let h1 = hash_kmer(kmer);
        let h2 = hash_kmer(kmer);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_hash_case_insensitive() {
        let h1 = hash_kmer(b"ACGTACGTACGTACGTACGTG");
        let h2 = hash_kmer(b"acgtacgtacgtacgtacgtg");
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_binary_search_containment() {
        // Create a small filter with known hashes
        let test_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut hashes = Vec::new();
        for kmer in test_seq.windows(RRNA_K) {
            hashes.push(hash_kmer(kmer));
            hashes.push(hash_kmer_rc(kmer));
        }
        hashes.sort_unstable();
        hashes.dedup();

        let filter = RrnaFilter {
            hashes,
            min_kmer_fraction: 0.3,
            stats: AtomicStats::new(),
            rrna_removed: AtomicU64::new(0),
        };

        // A read containing the same sequence should have high containment
        let containment = filter.containment(test_seq);
        assert!(
            containment > 0.9,
            "Same sequence should have high containment, got {}",
            containment
        );

        // Random sequence should have low containment
        let random: Vec<u8> = (0..150)
            .map(|i| [b'A', b'C', b'G', b'T'][(i * 7 + 3) % 4])
            .collect();
        let containment_rand = filter.containment(&random);
        assert!(
            containment_rand < 0.1,
            "Random sequence should have low containment, got {}",
            containment_rand
        );
    }

    #[test]
    fn test_filter_round_trip() {
        // Build a small filter, save, reload, verify
        let test_seqs = vec![
            b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"TGCATGCATGCATGCATGCATGCATGCA".to_vec(),
        ];
        let mut hashes = Vec::new();
        for seq in &test_seqs {
            RrnaFilter::extract_kmers(seq, &mut hashes);
        }
        hashes.sort_unstable();
        hashes.dedup();
        let original_len = hashes.len();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path();

        // Write
        {
            let file = std::fs::File::create(path).unwrap();
            let mut writer = BufWriter::new(file);
            writer.write_all(RRNA_FILTER_MAGIC).unwrap();
            writer
                .write_all(&(hashes.len() as u64).to_le_bytes())
                .unwrap();
            let hash_bytes = unsafe {
                std::slice::from_raw_parts(hashes.as_ptr() as *const u8, hashes.len() * 8)
            };
            writer.write_all(hash_bytes).unwrap();
        }

        // Reload
        let loaded = RrnaFilter::load_hashes(path).unwrap();
        assert_eq!(loaded.len(), original_len);
        assert_eq!(loaded, hashes);
    }
}
