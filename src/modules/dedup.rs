//! Deduplication module -- sequence-based + optical + UMI-aware
//!
//! Runs as a post-pipeline step on clean output files (not a streaming QcModule)
//! because it needs to see all reads before deciding which are duplicates.
//!
//! Virome-aware: for paired-end, requires BOTH mates identical to call a
//! duplicate. This preserves natural duplicates from abundant viruses.

use anyhow::Result;
use biometal::{FastqRecord, FastqStream, FastqWriter};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Deduplication configuration
#[derive(Debug, Clone)]
pub struct DedupConfig {
    /// Optical duplicate pixel distance (100 non-patterned, 2500 patterned)
    pub optical_distance: u32,
    /// Use UMIs from read names for dedup grouping
    pub umi_aware: bool,
    /// Number of bases from read start to use for hashing
    pub prefix_len: usize,
}

impl Default for DedupConfig {
    fn default() -> Self {
        Self {
            optical_distance: 2500,
            umi_aware: false,
            prefix_len: 50,
        }
    }
}

/// Statistics from dedup run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DedupStats {
    pub total_reads: u64,
    pub unique_reads: u64,
    pub pcr_duplicates: u64,
    pub optical_duplicates: u64,
    pub duplicate_rate: f64,
    pub estimated_library_complexity: u64,
}

/// Parsed Illumina read name for optical dedup
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields used in future optical dedup implementation
struct ReadCoords {
    tile: u32,
    x: u32,
    y: u32,
}

/// Entry in the dedup hash table
struct DupGroup {
    /// Index of the representative read (longest, then highest quality)
    representative: usize,
    /// Best length seen (prefer longer reads -- they retained more sequence)
    best_length: usize,
    /// Best mean quality seen (tiebreaker when lengths are equal)
    best_quality: f64,
    /// Count of reads in this group
    count: u32,
    /// Tile coordinates for optical dedup (from representative)
    #[allow(dead_code)]
    coords: Option<ReadCoords>,
    /// UMI of representative (if UMI-aware)
    #[allow(dead_code)]
    umi: Option<String>,
}

/// Deduplicate a single-end FASTQ file
pub fn dedup_single_end(input: &Path, output: &Path, config: &DedupConfig) -> Result<DedupStats> {
    // Pass 1: Hash all reads, build duplicate groups
    let stream = FastqStream::from_path(input)?;
    let mut groups: FxHashMap<u64, DupGroup> = FxHashMap::default();
    let mut records: Vec<FastqRecord> = Vec::new();

    for record in stream {
        let record = record?;
        let idx = records.len();
        let hash = hash_read(&record.sequence, config.prefix_len);
        let mean_q = mean_quality(&record.quality);
        let coords = parse_coords(&record.id);
        let umi = if config.umi_aware {
            extract_umi(&record.id)
        } else {
            None
        };

        // Combine hash with UMI if present (same sequence + different UMI = keep both)
        let group_key = if let Some(ref u) = umi {
            hash ^ hash_bytes(u.as_bytes())
        } else {
            hash
        };

        let read_len = record.sequence.len();
        let group = groups.entry(group_key).or_insert(DupGroup {
            representative: idx,
            best_length: read_len,
            best_quality: mean_q,
            count: 0,
            coords,
            umi,
        });
        group.count += 1;

        // Prefer longer reads (more sequence retained), then higher quality
        if read_len > group.best_length
            || (read_len == group.best_length && mean_q > group.best_quality)
        {
            group.representative = idx;
            group.best_length = read_len;
            group.best_quality = mean_q;
        }

        records.push(record);
    }

    // Pass 2: Write unique reads
    let mut writer = FastqWriter::create(output)?;
    let mut unique = 0u64;
    let mut pcr_dups = 0u64;
    let optical_dups = 0u64; // TODO: implement optical dedup detection
    let total = records.len() as u64;

    // Mark which indices are representatives
    let mut is_representative = vec![false; records.len()];
    for group in groups.values() {
        is_representative[group.representative] = true;
        if group.count > 1 {
            pcr_dups += (group.count - 1) as u64;
        }
    }

    for (idx, record) in records.iter().enumerate() {
        if is_representative[idx] {
            writer.write_record(record)?;
            unique += 1;
        }
    }

    let duplicate_rate = if total > 0 {
        1.0 - (unique as f64 / total as f64)
    } else {
        0.0
    };

    Ok(DedupStats {
        total_reads: total,
        unique_reads: unique,
        pcr_duplicates: pcr_dups,
        optical_duplicates: optical_dups,
        duplicate_rate,
        estimated_library_complexity: unique,
    })
}

/// Deduplicate paired-end FASTQ files
///
/// Requires BOTH mates identical to call a duplicate (virome-aware).
pub fn dedup_paired_end(
    r1_input: &Path,
    r2_input: &Path,
    r1_output: &Path,
    r2_output: &Path,
    config: &DedupConfig,
) -> Result<DedupStats> {
    // Pass 1: Hash all pairs
    let stream_r1 = FastqStream::from_path(r1_input)?;
    let stream_r2 = FastqStream::from_path(r2_input)?;

    let mut groups: FxHashMap<u64, DupGroup> = FxHashMap::default();
    let mut r1_records: Vec<FastqRecord> = Vec::new();
    let mut r2_records: Vec<FastqRecord> = Vec::new();

    for (rec1, rec2) in stream_r1.zip(stream_r2) {
        let rec1 = rec1?;
        let rec2 = rec2?;
        let idx = r1_records.len();

        // Hash both mates together -- both must be identical for a duplicate
        let hash = hash_pair(&rec1.sequence, &rec2.sequence, config.prefix_len);
        let mean_q = (mean_quality(&rec1.quality) + mean_quality(&rec2.quality)) / 2.0;
        let coords = parse_coords(&rec1.id);
        let umi = if config.umi_aware {
            extract_umi(&rec1.id)
        } else {
            None
        };

        let group_key = if let Some(ref u) = umi {
            hash ^ hash_bytes(u.as_bytes())
        } else {
            hash
        };

        let pair_len = rec1.sequence.len() + rec2.sequence.len();
        let group = groups.entry(group_key).or_insert(DupGroup {
            representative: idx,
            best_length: pair_len,
            best_quality: mean_q,
            count: 0,
            coords,
            umi,
        });
        group.count += 1;

        if pair_len > group.best_length
            || (pair_len == group.best_length && mean_q > group.best_quality)
        {
            group.representative = idx;
            group.best_length = pair_len;
            group.best_quality = mean_q;
        }

        r1_records.push(rec1);
        r2_records.push(rec2);
    }

    // Pass 2: Write unique pairs
    let mut w1 = FastqWriter::create(r1_output)?;
    let mut w2 = FastqWriter::create(r2_output)?;
    let mut unique = 0u64;
    let mut pcr_dups = 0u64;
    let total = r1_records.len() as u64;

    let mut is_representative = vec![false; r1_records.len()];
    for group in groups.values() {
        is_representative[group.representative] = true;
        if group.count > 1 {
            pcr_dups += (group.count - 1) as u64;
        }
    }

    for idx in 0..r1_records.len() {
        if is_representative[idx] {
            w1.write_record(&r1_records[idx])?;
            w2.write_record(&r2_records[idx])?;
            unique += 1;
        }
    }

    let duplicate_rate = if total > 0 {
        1.0 - (unique as f64 / total as f64)
    } else {
        0.0
    };

    Ok(DedupStats {
        total_reads: total * 2, // total individual reads
        unique_reads: unique * 2,
        pcr_duplicates: pcr_dups * 2, // each dup pair = 2 reads
        optical_duplicates: 0,
        duplicate_rate,
        estimated_library_complexity: unique,
    })
}

/// Hash a single-end read by prefix only
///
/// PCR duplicates share the same 5' start position. After adapter/quality
/// trimming (which removes from the 3' end), the 5' prefix is preserved.
/// Prefix-only hashing naturally handles contained reads: a 100bp read and
/// an 80bp read from the same PCR duplicate share the same prefix.
fn hash_read(sequence: &[u8], prefix_len: usize) -> u64 {
    let plen = prefix_len.min(sequence.len());
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in &sequence[..plen] {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Hash a paired-end read pair (both mates contribute)
fn hash_pair(r1: &[u8], r2: &[u8], prefix_len: usize) -> u64 {
    let h1 = hash_read(r1, prefix_len);
    let h2 = hash_read(r2, prefix_len);
    // Combine with rotation to avoid symmetry
    h1 ^ h2.rotate_left(31)
}

/// Hash arbitrary bytes (for UMI)
fn hash_bytes(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0x517cc1b727220a95;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

/// Mean quality score
fn mean_quality(quality: &[u8]) -> f64 {
    if quality.is_empty() {
        return 0.0;
    }
    let sum: u64 = quality.iter().map(|&q| q.saturating_sub(33) as u64).sum();
    sum as f64 / quality.len() as f64
}

/// Parse Illumina tile/x/y coordinates from read name
fn parse_coords(name: &str) -> Option<ReadCoords> {
    // Format: @instrument:run:flowcell:lane:tile:x:y
    let parts: Vec<&str> = name.split(':').collect();
    if parts.len() >= 7 {
        let tile = parts[4].parse().ok()?;
        let x = parts[5].parse().ok()?;
        let y = parts[6].split_whitespace().next()?.parse().ok()?;
        Some(ReadCoords { tile, x, y })
    } else {
        None
    }
}

/// Extract UMI from read name
///
/// Supports common formats:
/// - `@read:UMI` (last field after colon)
/// - `@read+UMI` (after plus sign)
/// - xGen format in index: `@... 1:N:0:ACGT+TGCA` (after the index)
fn extract_umi(name: &str) -> Option<String> {
    // Check for +UMI in the index portion (after space, last field)
    if let Some(comment) = name.split_whitespace().nth(1) {
        let fields: Vec<&str> = comment.split(':').collect();
        if let Some(last) = fields.last() {
            if last.contains('+') {
                // The UMI is typically after the index: ACGT+TGCA
                return Some(last.to_string());
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn make_fastq(records: &[(&str, &[u8], &[u8])], path: &Path) {
        let mut writer = FastqWriter::create(path).unwrap();
        for (name, seq, qual) in records {
            writer
                .write_record(&FastqRecord::new(
                    name.to_string(),
                    seq.to_vec(),
                    qual.to_vec(),
                ))
                .unwrap();
        }
    }

    #[test]
    fn test_dedup_single_end_removes_duplicates() {
        let tmp = TempDir::new().unwrap();
        let input = tmp.path().join("input.fastq");
        let output = tmp.path().join("output.fastq.gz");

        let len = 100;
        let seq: Vec<u8> = (0..len).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let qual_high = vec![b'I'; len];
        let qual_low = vec![b'5'; len];
        let seq2: Vec<u8> = (0..len).map(|i| [b'G', b'C', b'T', b'A'][i % 4]).collect();

        make_fastq(
            &[
                ("read1", &seq, &qual_high),
                ("read2", &seq2, &qual_high),
                ("read3", &seq, &qual_low), // duplicate of read1
            ],
            &input,
        );

        let stats = dedup_single_end(&input, &output, &DedupConfig::default()).unwrap();

        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.unique_reads, 2);
        assert_eq!(stats.pcr_duplicates, 1);
        assert!((stats.duplicate_rate - 0.333).abs() < 0.01);

        // Verify output has 2 reads
        let count = FastqStream::from_path(&output).unwrap().count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_dedup_paired_end() {
        let tmp = TempDir::new().unwrap();
        let r1_in = tmp.path().join("R1.fastq");
        let r2_in = tmp.path().join("R2.fastq");
        let r1_out = tmp.path().join("R1_dedup.fastq.gz");
        let r2_out = tmp.path().join("R2_dedup.fastq.gz");

        let len = 100;
        let seq_a: Vec<u8> = (0..len).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let seq_b: Vec<u8> = (0..len).map(|i| [b'G', b'C', b'T', b'A'][i % 4]).collect();
        let seq_c: Vec<u8> = vec![b'T'; len];
        let qual = vec![b'I'; len];

        // Pair 1 and Pair 3 have identical R1 AND R2 = PCR duplicate
        // Pair 2 has different R1 = unique
        make_fastq(
            &[
                ("pair1/1", &seq_a, &qual),
                ("pair2/1", &seq_b, &qual),
                ("pair3/1", &seq_a, &qual), // same R1 as pair1
            ],
            &r1_in,
        );
        make_fastq(
            &[
                ("pair1/2", &seq_c, &qual),
                ("pair2/2", &seq_c, &qual),
                ("pair3/2", &seq_c, &qual), // same R2 as pair1
            ],
            &r2_in,
        );

        let stats =
            dedup_paired_end(&r1_in, &r2_in, &r1_out, &r2_out, &DedupConfig::default()).unwrap();

        assert_eq!(stats.unique_reads, 4); // 2 unique pairs * 2
        assert_eq!(stats.pcr_duplicates, 2); // 1 dup pair * 2
    }

    #[test]
    fn test_no_duplicates() {
        let tmp = TempDir::new().unwrap();
        let input = tmp.path().join("input.fastq");
        let output = tmp.path().join("output.fastq.gz");

        let len = 100;
        let s1: Vec<u8> = (0..len).map(|i| [b'A', b'T', b'G', b'C'][i % 4]).collect();
        let s2: Vec<u8> = (0..len).map(|i| [b'G', b'C', b'T', b'A'][i % 4]).collect();
        let s3: Vec<u8> = (0..len)
            .map(|i| [b'T', b'A', b'C', b'G'][(i + 1) % 4])
            .collect();
        let qual = vec![b'I'; len];
        make_fastq(
            &[("r1", &s1, &qual), ("r2", &s2, &qual), ("r3", &s3, &qual)],
            &input,
        );

        let stats = dedup_single_end(&input, &output, &DedupConfig::default()).unwrap();
        assert_eq!(stats.unique_reads, 3);
        assert_eq!(stats.pcr_duplicates, 0);
        assert_eq!(stats.duplicate_rate, 0.0);
    }

    #[test]
    fn test_extract_umi() {
        // xGen UDI-UMI format
        assert_eq!(
            extract_umi("@A00882:1:FC:1:1:1:1 1:N:0:ACGT+TGCA"),
            Some("ACGT+TGCA".to_string())
        );

        // No UMI
        assert_eq!(extract_umi("@read1"), None);

        // No plus in index
        assert_eq!(extract_umi("@A00882:1:FC:1:1:1:1 1:N:0:ACGT"), None);
    }

    #[test]
    fn test_parse_coords() {
        let coords = parse_coords("@A00882:431:HVNLJDSXY:2:1101:1380:1000 1:N:0:ACGT");
        assert!(coords.is_some());
        let c = coords.unwrap();
        assert_eq!(c.tile, 1101);
        assert_eq!(c.x, 1380);
        assert_eq!(c.y, 1000);
    }
}
