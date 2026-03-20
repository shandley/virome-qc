//! Quick-scan ingestion: read first N records to characterize the data

use crate::ingest::platform::PlatformInfo;
use anyhow::Result;
use biometal::operations::{complexity_score, gc_content, mean_quality};
use biometal::FastqStream;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Number of reads to scan for quick characterization
const SCAN_SIZE: usize = 10_000;

/// Adapter core sequences for quick contamination check (first 12bp of each)
const ADAPTER_CORES: &[(&str, &[u8])] = &[
    ("TruSeq", b"AGATCGGAAGAG"),
    ("Nextera", b"CTGTCTCTTATACACATCT"),
    ("NEBNext", b"AGATCGGAAGAG"),
];

/// Result of ingestion analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IngestResult {
    /// Platform detection from FASTQ headers
    pub platform: Option<PlatformInfo>,
    /// Read characteristics
    pub reads: ReadInfo,
    /// Quick-scan statistics from first N reads
    pub quick_scan: QuickScanStats,
    /// Recommendations for QC configuration
    pub recommendations: Vec<Recommendation>,
    /// Warnings about potential issues
    pub warnings: Vec<String>,
}

/// Read-level characteristics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadInfo {
    /// Paired-end or single-end
    pub paired: bool,
    /// Read length (mode of observed lengths)
    pub read_length: usize,
    /// Whether read length is variable
    pub variable_length: bool,
    /// Min observed read length
    pub min_length: usize,
    /// Max observed read length
    pub max_length: usize,
    /// Quality encoding (33 = Phred+33, 64 = Phred+64)
    pub quality_offset: u8,
    /// Estimated total read count (from file size)
    pub estimated_read_count: Option<u64>,
}

/// Quick-scan statistics from first N reads
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuickScanStats {
    /// Number of reads scanned
    pub reads_scanned: usize,
    /// Mean quality across scanned reads
    pub mean_quality: f64,
    /// Mean GC content (0.0 - 1.0)
    pub mean_gc: f64,
    /// Fraction of bases that are N
    pub n_rate: f64,
    /// Estimated adapter contamination rate by type
    pub adapter_rates: HashMap<String, f64>,
    /// Most likely adapter type
    pub dominant_adapter: Option<String>,
    /// N-rate at position 0 (NextSeq cluster init indicator)
    pub n_rate_pos0: f64,
    /// Whether Q-score binning is detected (e.g., NovaSeq 3-bin)
    pub quality_binned: bool,
    /// Number of distinct quality values observed
    pub distinct_quality_values: usize,
    /// Complexity score distribution: 2nd percentile (data-driven threshold)
    pub complexity_p2: f64,
    /// Complexity score distribution: 5th percentile
    pub complexity_p5: f64,
    /// Complexity score distribution: median
    pub complexity_median: f64,
}

/// A recommendation for QC configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Recommendation {
    pub parameter: String,
    pub value: String,
    pub reason: String,
}

/// Run ingestion analysis on FASTQ file(s)
///
/// Scans the first 10,000 reads to detect platform, read characteristics,
/// and potential issues. Runs in <100ms for gzipped input.
pub fn ingest_fastq(r1_path: &Path, r2_path: Option<&Path>) -> Result<IngestResult> {
    let stream = FastqStream::from_path(r1_path)?;

    let paired = r2_path.is_some();
    let mut platform: Option<PlatformInfo> = None;

    // Accumulators
    let mut total_reads = 0usize;
    let mut total_bases = 0u64;
    let mut total_n_bases = 0u64;
    let mut total_gc = 0.0f64;
    let mut total_quality = 0.0f64;
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    let mut min_length = usize::MAX;
    let mut max_length = 0usize;
    let mut min_qual_byte = 255u8;

    // Adapter detection
    let mut adapter_hits: HashMap<String, usize> = HashMap::new();

    // Position 0 N-count
    let mut pos0_total = 0usize;
    let mut pos0_n = 0usize;

    // Quality value tracking (for binning detection)
    let mut quality_value_seen = [false; 94]; // Q0-Q93

    // Complexity score collection
    let mut complexity_scores: Vec<f64> = Vec::with_capacity(SCAN_SIZE);

    for record in stream {
        let record = record?;

        if total_reads == 0 {
            platform = PlatformInfo::from_header(&record.id);
        }

        let seq = &record.sequence;
        let qual = &record.quality;

        // Read length
        let len = seq.len();
        *length_counts.entry(len).or_insert(0) += 1;
        min_length = min_length.min(len);
        max_length = max_length.max(len);

        // Quality
        total_quality += mean_quality(qual);
        for &q in qual {
            min_qual_byte = min_qual_byte.min(q);
            let phred = q.saturating_sub(33) as usize;
            if phred < quality_value_seen.len() {
                quality_value_seen[phred] = true;
            }
        }

        // GC content
        total_gc += gc_content(seq);

        // Complexity score
        complexity_scores.push(complexity_score(seq));

        // N bases
        let n_count = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
        total_n_bases += n_count as u64;
        total_bases += len as u64;

        // Position 0 N-rate
        if !seq.is_empty() {
            pos0_total += 1;
            if seq[0] == b'N' || seq[0] == b'n' {
                pos0_n += 1;
            }
        }

        // Quick adapter scan (check last 20bp for adapter core)
        if len >= 20 {
            let tail = &seq[len.saturating_sub(20)..];
            for (name, core) in ADAPTER_CORES {
                let check_len = core.len().min(tail.len());
                let mismatches: usize = tail[..check_len]
                    .iter()
                    .zip(core[..check_len].iter())
                    .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
                    .count();
                if (mismatches as f64 / check_len as f64) < 0.15 {
                    *adapter_hits.entry(name.to_string()).or_insert(0) += 1;
                }
            }
        }

        total_reads += 1;
        if total_reads >= SCAN_SIZE {
            break;
        }
    }

    if total_reads == 0 {
        anyhow::bail!("No reads found in {}", r1_path.display());
    }

    // Compute statistics
    let mean_gc = total_gc / total_reads as f64;
    let mean_qual = total_quality / total_reads as f64;
    let n_rate = total_n_bases as f64 / total_bases as f64;
    let n_rate_pos0 = if pos0_total > 0 {
        pos0_n as f64 / pos0_total as f64
    } else {
        0.0
    };

    // Mode read length
    let read_length = length_counts
        .iter()
        .max_by_key(|(_, &count)| count)
        .map(|(&len, _)| len)
        .unwrap_or(0);
    let variable_length = length_counts.len() > 1;

    // Quality encoding: Phred+33 (modern standard) vs Phred+64 (Illumina 1.3-1.7, pre-2011)
    // If min quality byte < 59, definitely Phred+33 (can't be a valid Phred+64 quality)
    // Ambiguous range 59-64 exists but all modern data is Phred+33
    let quality_offset = if min_qual_byte < 59 { 33 } else if min_qual_byte >= 75 { 64 } else { 33 };

    // Quality binning detection
    let distinct_quality_values = quality_value_seen.iter().filter(|&&v| v).count();
    let quality_binned = distinct_quality_values <= 8; // NovaSeq uses 4 values (2,12,23,37)

    // Adapter rates
    let adapter_rates: HashMap<String, f64> = adapter_hits
        .iter()
        .map(|(name, &count)| (name.clone(), count as f64 / total_reads as f64))
        .collect();
    let dominant_adapter = adapter_rates
        .iter()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .filter(|(_, &rate)| rate > 0.005) // at least 0.5%
        .map(|(name, _)| name.clone());

    // Estimate total reads from file size
    let estimated_read_count = estimate_read_count(r1_path, total_reads, total_bases);

    // Compute complexity percentiles for data-driven threshold setting
    complexity_scores.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let complexity_p2 = percentile_sorted(&complexity_scores, 0.02);
    let complexity_p5 = percentile_sorted(&complexity_scores, 0.05);
    let complexity_median = percentile_sorted(&complexity_scores, 0.50);

    let quick_scan = QuickScanStats {
        reads_scanned: total_reads,
        mean_quality: mean_qual,
        mean_gc,
        n_rate,
        adapter_rates,
        dominant_adapter,
        n_rate_pos0,
        quality_binned,
        distinct_quality_values,
        complexity_p2,
        complexity_p5,
        complexity_median,
    };

    let reads = ReadInfo {
        paired,
        read_length,
        variable_length,
        min_length,
        max_length,
        quality_offset,
        estimated_read_count,
    };

    // Generate recommendations and warnings
    let mut recommendations = Vec::new();
    let mut warnings = Vec::new();

    if let Some(ref plat) = platform {
        // Poly-G
        if plat.has_poly_g_risk() {
            recommendations.push(Recommendation {
                parameter: "polyx.platform_aware".into(),
                value: "true".into(),
                reason: format!(
                    "{} uses 2-color chemistry: aggressive poly-G trimming recommended",
                    plat.model
                ),
            });
        }

        // Optical dup distance
        if plat.patterned_flowcell {
            recommendations.push(Recommendation {
                parameter: "dedup.optical_distance".into(),
                value: "2500".into(),
                reason: format!(
                    "{} has patterned flowcell: use 2500px optical dup distance",
                    plat.model
                ),
            });
        }

        // Q-score binning
        if let Some(bins) = plat.quality_bins {
            warnings.push(format!(
                "{} uses {}-bin Q-score quantization: quality values are approximate",
                plat.model, bins
            ));
        }
    }

    // Read length recommendations
    if read_length <= 100 {
        recommendations.push(Recommendation {
            parameter: "quality.min_length".into(),
            value: format!("{}", read_length * 2 / 5), // ~40% of read length
            reason: format!(
                "Short reads ({read_length}bp): lower min_length threshold recommended"
            ),
        });
        recommendations.push(Recommendation {
            parameter: "quality.window_size".into(),
            value: format!("{}", (read_length / 8).max(4)),
            reason: format!(
                "Short reads ({read_length}bp): smaller quality window recommended"
            ),
        });
    }

    // N-rate warnings
    if n_rate_pos0 > 0.01 {
        warnings.push(format!(
            "High N-rate at position 0 ({:.1}%): possible cluster initialization issue",
            n_rate_pos0 * 100.0
        ));
    }
    if n_rate > 0.005 {
        warnings.push(format!(
            "Elevated N-rate ({:.2}%): {:.0}% of bases are ambiguous",
            n_rate * 100.0,
            n_rate * 100.0
        ));
    }

    // Adapter warnings
    if let Some(ref adapter) = quick_scan.dominant_adapter {
        let rate = quick_scan.adapter_rates.get(adapter).copied().unwrap_or(0.0);
        if rate > 0.05 {
            warnings.push(format!(
                "High adapter contamination ({:.1}% {adapter}): consider verifying insert size distribution",
                rate * 100.0
            ));
        }
    }

    // GC content warnings
    if mean_gc < 0.30 {
        warnings.push(format!(
            "Low GC content ({:.1}%): possible AT-rich sample or GC bias",
            mean_gc * 100.0
        ));
    } else if mean_gc > 0.65 {
        warnings.push(format!(
            "High GC content ({:.1}%): possible GC-rich sample or GC bias",
            mean_gc * 100.0
        ));
    }

    // Quality encoding warning
    if quality_offset == 64 {
        warnings.push("Phred+64 quality encoding detected: old Illumina format".into());
    }

    Ok(IngestResult {
        platform,
        reads,
        quick_scan,
        recommendations,
        warnings,
    })
}

/// Get a percentile value from a sorted slice
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() as f64 * p) as usize).min(sorted.len() - 1);
    sorted[idx]
}

/// Estimate total read count from file size
fn estimate_read_count(path: &Path, scanned_reads: usize, scanned_bases: u64) -> Option<u64> {
    let file_size = std::fs::metadata(path).ok()?.len();
    if scanned_reads == 0 || scanned_bases == 0 {
        return None;
    }

    // Estimate: compressed bytes per base (from scan sample)
    // FASTQ overhead: ~4 lines per read, quality = sequence length
    // Rough: each read = 4 * read_length bytes uncompressed, ~0.4x compressed
    let avg_read_len = scanned_bases as f64 / scanned_reads as f64;
    let uncompressed_per_read = avg_read_len * 2.5; // seq + qual + headers + newlines
    let compression_ratio = 0.35; // typical gzip ratio for FASTQ
    let compressed_per_read = uncompressed_per_read * compression_ratio;

    Some((file_size as f64 / compressed_per_read) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ingest::platform::Chemistry;
    use tempfile::TempDir;

    fn write_test_fastq(path: &Path, instrument: &str, num_reads: usize, read_len: usize) {
        use std::io::Write;
        let mut f = std::fs::File::create(path).unwrap();
        for i in 0..num_reads {
            let seq: String = (0..read_len)
                .map(|j| ["A", "C", "G", "T"][(i + j) % 4])
                .collect();
            let qual = "I".repeat(read_len);
            writeln!(f, "@{instrument}:1:FC:1:1:{i}:1 1:N:0:ACGT").unwrap();
            writeln!(f, "{seq}").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{qual}").unwrap();
        }
    }

    #[test]
    fn test_ingest_novaseq_data() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test.fastq");
        write_test_fastq(&path, "A00882", 100, 150);

        let result = ingest_fastq(&path, None).unwrap();

        assert!(result.platform.is_some());
        let plat = result.platform.unwrap();
        assert_eq!(plat.model, "NovaSeq 6000");
        assert_eq!(plat.chemistry, Chemistry::TwoColor);
        assert!(plat.patterned_flowcell);

        assert_eq!(result.reads.read_length, 150);
        assert!(!result.reads.paired);
        assert_eq!(result.reads.quality_offset, 33);
    }

    #[test]
    fn test_ingest_miseq_data() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test.fastq");
        write_test_fastq(&path, "M01757", 100, 250);

        let result = ingest_fastq(&path, None).unwrap();

        let plat = result.platform.unwrap();
        assert_eq!(plat.model, "MiSeq");
        assert!(!plat.has_poly_g_risk());
        assert_eq!(result.reads.read_length, 250);
    }

    #[test]
    fn test_ingest_short_reads() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test.fastq");
        write_test_fastq(&path, "NB501138", 100, 75);

        let result = ingest_fastq(&path, None).unwrap();

        assert_eq!(result.reads.read_length, 75);
        // Should recommend lower min_length for short reads
        assert!(
            result
                .recommendations
                .iter()
                .any(|r| r.parameter == "quality.min_length"),
            "Should recommend min_length for short reads"
        );
    }

    #[test]
    fn test_ingest_paired() {
        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("R1.fastq");
        let r2 = tmp.path().join("R2.fastq");
        write_test_fastq(&r1, "A00882", 100, 150);
        write_test_fastq(&r2, "A00882", 100, 150);

        let result = ingest_fastq(&r1, Some(&r2)).unwrap();
        assert!(result.reads.paired);
    }
}
