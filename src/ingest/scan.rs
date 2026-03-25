//! Quick-scan ingestion: read first N records to characterize the data

use crate::ingest::platform::PlatformInfo;
use anyhow::Result;
use biometal::operations::{complexity_score, gc_content, mean_quality};
use biometal::FastqStream;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Number of reads to scan for quick characterization
const SCAN_SIZE: usize = 50_000;

/// Number of per-position quality buckets to track
const POSITION_BUCKETS: usize = 10;

/// Adapter families with distinguishing sequences (forward + RC for read-through)
/// Each family has a unique identifying sequence that doesn't overlap with others.
struct AdapterFamily {
    name: &'static str,
    /// Config name to use in profile (matches load_adapter_set keys)
    config_names: &'static [&'static str],
    /// Sequences to scan for (forward orientation, first N bp used as probe)
    probes: &'static [&'static [u8]],
}

const ADAPTER_FAMILIES: &[AdapterFamily] = &[
    AdapterFamily {
        name: "Nextera",
        config_names: &["nextera"],
        probes: &[
            // Nextera transposase -- unique sequence not shared with TruSeq
            b"TCGTCGGCAGCGTC",
            b"GTCTCGTGGGCTCGG",
            // RC forms (appear during read-through)
            b"CTGTCTCTTATACACATCT",
        ],
    },
    AdapterFamily {
        name: "TruSeq/NEBNext",
        // TruSeq and NEBNext share identical adapter sequences
        // Both config names included so either will match profile
        config_names: &["truseq", "truseq_universal", "nebnext"],
        probes: &[
            b"AGATCGGAAGAGCACACGTCTGAAC", // R1
            b"AGATCGGAAGAGCGTCGTGTAGGGA", // R2
        ],
    },
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
    /// Adapter detection rates by family name
    pub adapter_rates: HashMap<String, f64>,
    /// Most likely adapter family (highest hit rate above threshold)
    pub dominant_adapter: Option<String>,
    /// Config names to use for the detected adapter family
    pub detected_adapter_configs: Vec<String>,
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
    /// Per-position quality profile (10 evenly-spaced buckets across read length)
    /// Each value is mean Phred quality at that fractional position
    pub quality_profile: Vec<f64>,
    /// 3' adapter hit rate (read-through, separate from internal)
    pub adapter_3prime_rate: f64,
    /// Internal adapter chimera rate
    pub adapter_internal_rate: f64,
    /// R2 mean quality (None if single-end or R2 not scanned)
    pub r2_mean_quality: Option<f64>,
    /// R2 quality profile (same bucketing as R1)
    pub r2_quality_profile: Option<Vec<f64>>,
    /// Quality drop from R1 to R2 (negative = R2 is worse)
    pub r2_quality_delta: Option<f64>,
}

/// A recommendation for QC configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Recommendation {
    pub parameter: String,
    pub value: String,
    pub reason: String,
}

/// Accumulator for per-position quality across multiple reads
struct PositionQualityAccum {
    sums: Vec<f64>,
    counts: Vec<u64>,
    num_buckets: usize,
}

impl PositionQualityAccum {
    fn new(num_buckets: usize) -> Self {
        Self {
            sums: vec![0.0; num_buckets],
            counts: vec![0; num_buckets],
            num_buckets,
        }
    }

    fn add(&mut self, qual: &[u8]) {
        let len = qual.len();
        if len == 0 {
            return;
        }
        for bucket in 0..self.num_buckets {
            let pos = (bucket * (len - 1)) / (self.num_buckets - 1).max(1);
            let pos = pos.min(len - 1);
            let q = qual[pos].saturating_sub(33) as f64;
            self.sums[bucket] += q;
            self.counts[bucket] += 1;
        }
    }

    fn means(&self) -> Vec<f64> {
        self.sums
            .iter()
            .zip(self.counts.iter())
            .map(|(&s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
            .collect()
    }
}

/// Run ingestion analysis on FASTQ file(s)
///
/// Scans the first 50,000 reads to detect platform, read characteristics,
/// and potential issues. Runs in <1s for gzipped input.
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

    // Adapter detection -- separate 3' tail and internal counts
    let mut adapter_3prime_hits: HashMap<String, usize> = HashMap::new();
    let mut adapter_internal_hits: HashMap<String, usize> = HashMap::new();
    let mut adapter_any_hits: HashMap<String, usize> = HashMap::new();
    for family in ADAPTER_FAMILIES {
        adapter_3prime_hits.insert(family.name.to_string(), 0);
        adapter_internal_hits.insert(family.name.to_string(), 0);
        adapter_any_hits.insert(family.name.to_string(), 0);
    }

    // Position 0 N-count
    let mut pos0_total = 0usize;
    let mut pos0_n = 0usize;

    // Quality value tracking (for binning detection)
    let mut quality_value_seen = [false; 94]; // Q0-Q93

    // Complexity score collection
    let mut complexity_scores: Vec<f64> = Vec::with_capacity(SCAN_SIZE);

    // Per-position quality profile
    let mut pos_quality = PositionQualityAccum::new(POSITION_BUCKETS);

    let mut pe_hint_detected = false;

    for record in stream {
        let record = record?;

        if total_reads == 0 {
            platform = PlatformInfo::from_header(&record.id);

            // Check if read names suggest paired-end data when run as SE
            if !paired {
                let id = &record.id;
                if id.ends_with("/1") || id.ends_with("/2") || id.ends_with(" 1:") || id.ends_with(" 2:") {
                    pe_hint_detected = true;
                }
            }
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

        // Per-position quality
        pos_quality.add(qual);

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

        // Adapter family detection: scan for 3' tail and internal separately
        if len >= 20 {
            'family: for family in ADAPTER_FAMILIES {
                for probe in family.probes {
                    let probe_len = probe.len().min(15);
                    let probe_slice = &probe[..probe_len];

                    // Check 3' tail (last probe_len + 5 bases)
                    if len >= probe_len {
                        let tail_start = len.saturating_sub(probe_len + 5);
                        for window in seq[tail_start..].windows(probe_len) {
                            let mm: usize = window
                                .iter()
                                .zip(probe_slice.iter())
                                .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
                                .count();
                            if mm <= 2 {
                                *adapter_3prime_hits.get_mut(family.name).unwrap() += 1;
                                *adapter_any_hits.get_mut(family.name).unwrap() += 1;
                                continue 'family;
                            }
                        }
                    }

                    // Check internal (excluding last 20bp to avoid 3' overlap)
                    if len >= probe_len + 20 {
                        for window in seq[..len.saturating_sub(20)].windows(probe_len) {
                            let mm: usize = window
                                .iter()
                                .zip(probe_slice.iter())
                                .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
                                .count();
                            if mm <= 1 {
                                *adapter_internal_hits.get_mut(family.name).unwrap() += 1;
                                *adapter_any_hits.get_mut(family.name).unwrap() += 1;
                                continue 'family;
                            }
                        }
                    }
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

    // R2 scan (if paired-end)
    let r2_scan = if let Some(r2) = r2_path {
        scan_r2(r2, total_reads, total_quality / total_reads as f64)?
    } else {
        R2ScanResult {
            mean_quality: None,
            quality_profile: None,
            quality_delta: None,
        }
    };
    let r2_mean_quality = r2_scan.mean_quality;
    let r2_quality_profile = r2_scan.quality_profile;
    let r2_quality_delta = r2_scan.quality_delta;

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

    // Quality encoding
    let quality_offset = if min_qual_byte < 59 {
        33
    } else if min_qual_byte >= 75 {
        64
    } else {
        33
    };

    // Quality binning detection
    let distinct_quality_values = quality_value_seen.iter().filter(|&&v| v).count();
    let quality_binned = distinct_quality_values <= 8;

    // Adapter detection: determine dominant family and get config names
    let adapter_rates: HashMap<String, f64> = adapter_any_hits
        .iter()
        .filter(|(_, &count)| count > 0)
        .map(|(name, &count)| (name.clone(), count as f64 / total_reads as f64))
        .collect();

    let dominant_adapter = adapter_rates
        .iter()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .filter(|(_, &rate)| rate > 0.001)
        .map(|(name, _)| name.clone());

    let detected_adapter_configs = dominant_adapter
        .as_ref()
        .and_then(|name| {
            ADAPTER_FAMILIES
                .iter()
                .find(|f| f.name == name.as_str())
                .map(|f| f.config_names.iter().map(|s| s.to_string()).collect())
        })
        .unwrap_or_default();

    // Compute aggregate 3' and internal adapter rates across all families
    let total_3prime: usize = adapter_3prime_hits.values().sum();
    let total_internal: usize = adapter_internal_hits.values().sum();
    let adapter_3prime_rate = total_3prime as f64 / total_reads as f64;
    let adapter_internal_rate = total_internal as f64 / total_reads as f64;

    // Estimate total reads from file size
    let estimated_read_count = estimate_read_count(r1_path, total_reads, total_bases);

    // Compute complexity percentiles
    complexity_scores.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let complexity_p2 = percentile_sorted(&complexity_scores, 0.02);
    let complexity_p5 = percentile_sorted(&complexity_scores, 0.05);
    let complexity_median = percentile_sorted(&complexity_scores, 0.50);

    let quality_profile = pos_quality.means();

    let quick_scan = QuickScanStats {
        reads_scanned: total_reads,
        mean_quality: mean_qual,
        mean_gc,
        n_rate,
        adapter_rates,
        dominant_adapter,
        detected_adapter_configs,
        n_rate_pos0,
        quality_binned,
        distinct_quality_values,
        complexity_p2,
        complexity_p5,
        complexity_median,
        quality_profile,
        adapter_3prime_rate,
        adapter_internal_rate,
        r2_mean_quality,
        r2_quality_profile,
        r2_quality_delta,
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
            value: format!("{}", read_length * 2 / 5),
            reason: format!(
                "Short reads ({read_length}bp): lower min_length threshold recommended"
            ),
        });
        recommendations.push(Recommendation {
            parameter: "quality.window_size".into(),
            value: format!("{}", (read_length / 8).max(4)),
            reason: format!("Short reads ({read_length}bp): smaller quality window recommended"),
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
        let rate = quick_scan
            .adapter_rates
            .get(adapter)
            .copied()
            .unwrap_or(0.0);
        if rate > 0.05 {
            warnings.push(format!(
                "High adapter contamination ({:.1}% {adapter}): consider verifying insert size distribution",
                rate * 100.0
            ));
        }
    }

    // Internal adapter chimera warning (library prep signature)
    if adapter_internal_rate > 0.005 {
        warnings.push(format!(
            "Internal adapter chimeras detected ({:.2}%): typical for ligation-based library preps",
            adapter_internal_rate * 100.0
        ));
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

    // Quality profile shape warning: detect severe 3' degradation
    let qp = &quick_scan.quality_profile;
    if qp.len() >= 2 {
        let q_start = qp[0];
        let q_end = qp[qp.len() - 1];
        let q_drop = q_start - q_end;
        if q_drop > 10.0 {
            warnings.push(format!(
                "Severe 3' quality degradation: Q{:.0} at start, Q{:.0} at end (drop of {:.0})",
                q_start, q_end, q_drop
            ));
        }
    }

    // R2 quality warnings
    if let Some(delta) = quick_scan.r2_quality_delta {
        if delta < -5.0 {
            warnings.push(format!(
                "R2 quality is {:.1} Phred points lower than R1: expect higher singleton rate",
                -delta
            ));
        }
    }

    // PE data provided as SE warning
    if pe_hint_detected {
        warnings.push(
            "Read names contain /1 or /2 suffixes suggesting paired-end data, \
             but only one file was provided. Consider using -1 and -2 flags for \
             paired-end mode (enables concordant mate flagging and read merging)."
                .to_string(),
        );
    }

    Ok(IngestResult {
        platform,
        reads,
        quick_scan,
        recommendations,
        warnings,
    })
}

/// R2 scan results
struct R2ScanResult {
    mean_quality: Option<f64>,
    quality_profile: Option<Vec<f64>>,
    quality_delta: Option<f64>,
}

/// Quick scan of R2 file: compute mean quality and per-position profile
fn scan_r2(r2_path: &Path, max_reads: usize, r1_mean_quality: f64) -> Result<R2ScanResult> {
    let stream = FastqStream::from_path(r2_path)?;
    let mut total_quality = 0.0f64;
    let mut count = 0usize;
    let mut pos_quality = PositionQualityAccum::new(POSITION_BUCKETS);

    for record in stream {
        let record = record?;
        total_quality += mean_quality(&record.quality);
        pos_quality.add(&record.quality);
        count += 1;
        if count >= max_reads {
            break;
        }
    }

    if count == 0 {
        return Ok(R2ScanResult {
            mean_quality: None,
            quality_profile: None,
            quality_delta: None,
        });
    }

    let r2_mean = total_quality / count as f64;
    let delta = r2_mean - r1_mean_quality;

    Ok(R2ScanResult {
        mean_quality: Some(r2_mean),
        quality_profile: Some(pos_quality.means()),
        quality_delta: Some(delta),
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

    let avg_read_len = scanned_bases as f64 / scanned_reads as f64;
    let uncompressed_per_read = avg_read_len * 2.5;
    let compression_ratio = 0.35;
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

    fn write_test_fastq_degraded(
        path: &Path,
        instrument: &str,
        num_reads: usize,
        read_len: usize,
    ) {
        use std::io::Write;
        let mut f = std::fs::File::create(path).unwrap();
        for i in 0..num_reads {
            let seq: String = (0..read_len)
                .map(|j| ["A", "C", "G", "T"][(i + j) % 4])
                .collect();
            // Quality degrades from I (Q40) at start to 5 (Q20) at end
            let qual: String = (0..read_len)
                .map(|j| {
                    let frac = j as f64 / read_len as f64;
                    let q = 40.0 - (frac * 20.0);
                    (q as u8 + 33) as char
                })
                .collect();
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

        // Quality profile should have POSITION_BUCKETS entries
        assert_eq!(result.quick_scan.quality_profile.len(), POSITION_BUCKETS);
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
        assert!(result.quick_scan.r2_mean_quality.is_some());
        assert!(result.quick_scan.r2_quality_profile.is_some());
        // Same quality data -- delta should be ~0
        let delta = result.quick_scan.r2_quality_delta.unwrap();
        assert!(delta.abs() < 0.1, "Expected near-zero delta, got {delta}");
    }

    #[test]
    fn test_quality_degradation_detection() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test.fastq");
        write_test_fastq_degraded(&path, "M01757", 100, 250);

        let result = ingest_fastq(&path, None).unwrap();

        let qp = &result.quick_scan.quality_profile;
        // First bucket should be higher quality than last
        assert!(qp[0] > qp[qp.len() - 1], "Quality should degrade 3' to 5'");
        let drop = qp[0] - qp[qp.len() - 1];
        assert!(drop > 15.0, "Expected >15 Phred drop, got {drop:.1}");

        // Should generate a 3' degradation warning
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("3' quality degradation")),
            "Should warn about 3' degradation"
        );
    }

    #[test]
    fn test_r2_quality_degradation() {
        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("R1.fastq");
        let r2 = tmp.path().join("R2.fastq");
        write_test_fastq(&r1, "A00882", 100, 150); // high quality R1
        write_test_fastq_degraded(&r2, "A00882", 100, 150); // degraded R2

        let result = ingest_fastq(&r1, Some(&r2)).unwrap();

        let delta = result.quick_scan.r2_quality_delta.unwrap();
        assert!(delta < -5.0, "R2 should be >5 Phred worse than R1, got {delta:.1}");

        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("R2 quality")),
            "Should warn about R2 quality"
        );
    }

    #[test]
    fn test_adapter_rate_separation() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("test.fastq");
        write_test_fastq(&path, "A00882", 100, 150);

        let result = ingest_fastq(&path, None).unwrap();

        // With synthetic data, adapter rates should be 0
        assert_eq!(result.quick_scan.adapter_3prime_rate, 0.0);
        assert_eq!(result.quick_scan.adapter_internal_rate, 0.0);
    }
}
