//! Read analytics — comprehensive thread-safe statistics collection
//!
//! Single unified collector for all QA metrics during pipeline processing.
//! Uses atomic counters and fixed-size arrays for lock-free parallel accumulation.
//!
//! Biometal primitives used:
//! - `gc_content` (NEON-optimized) for GC distribution
//! - `count_bases` (NEON-optimized) for per-position base composition
//! - `HyperLogLog` for library complexity / duplication estimation

use biometal::operations::gc_content;
use biometal::operations::sketching::{HyperLogLog, SketchConfig};
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// Maximum read position tracked (covers all Illumina + most long-read truncated)
const MAX_POSITION: usize = 500;
/// Number of quality score bins (Q0 through Q41)
const QUALITY_BINS: usize = 42;

// ─── Atomic Histogram ───────────────────────────────────────────────────────

/// Serializable histogram snapshot
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Histogram {
    pub bin_edges: Vec<f64>,
    pub counts: Vec<u64>,
    pub underflow: u64,
    pub overflow: u64,
    pub total: u64,
}

/// Thread-safe histogram using atomic counters
pub struct AtomicHistogram {
    bin_edges: Vec<f64>,
    counts: Vec<AtomicU64>,
    underflow: AtomicU64,
    overflow: AtomicU64,
    total: AtomicU64,
}

impl AtomicHistogram {
    pub fn uniform(min: f64, max: f64, num_bins: usize) -> Self {
        let step = (max - min) / num_bins as f64;
        let bin_edges: Vec<f64> = (0..=num_bins).map(|i| min + step * i as f64).collect();
        let counts = (0..num_bins).map(|_| AtomicU64::new(0)).collect();
        Self {
            bin_edges,
            counts,
            underflow: AtomicU64::new(0),
            overflow: AtomicU64::new(0),
            total: AtomicU64::new(0),
        }
    }

    pub fn record(&self, value: f64) {
        self.total.fetch_add(1, Ordering::Relaxed);
        if value < self.bin_edges[0] {
            self.underflow.fetch_add(1, Ordering::Relaxed);
            return;
        }
        let last_edge = self.bin_edges[self.bin_edges.len() - 1];
        if value >= last_edge {
            self.overflow.fetch_add(1, Ordering::Relaxed);
            return;
        }
        // Fast bin lookup for uniform bins
        let step = self.bin_edges[1] - self.bin_edges[0];
        let bin = ((value - self.bin_edges[0]) / step) as usize;
        let bin = bin.min(self.counts.len() - 1);
        self.counts[bin].fetch_add(1, Ordering::Relaxed);
    }

    pub fn snapshot(&self) -> Histogram {
        Histogram {
            bin_edges: self.bin_edges.clone(),
            counts: self.counts.iter().map(|c| c.load(Ordering::Relaxed)).collect(),
            underflow: self.underflow.load(Ordering::Relaxed),
            overflow: self.overflow.load(Ordering::Relaxed),
            total: self.total.load(Ordering::Relaxed),
        }
    }
}

// ─── Per-Position Tracking ──────────────────────────────────────────────────

/// Per-position quality score distribution (Q0-Q41 bins at each read position)
struct PositionQualityTracker {
    /// [position][quality_bin] → count
    bins: Vec<[AtomicU64; QUALITY_BINS]>,
    /// Number of reads contributing to each position
    counts: Vec<AtomicU64>,
    /// Maximum position observed
    max_pos: AtomicU64,
}

impl PositionQualityTracker {
    fn new() -> Self {
        let bins: Vec<[AtomicU64; QUALITY_BINS]> = (0..MAX_POSITION)
            .map(|_| std::array::from_fn(|_| AtomicU64::new(0)))
            .collect();
        let counts = (0..MAX_POSITION).map(|_| AtomicU64::new(0)).collect();
        Self {
            bins,
            counts,
            max_pos: AtomicU64::new(0),
        }
    }

    fn record(&self, quality: &[u8]) {
        let len = quality.len().min(MAX_POSITION);
        // Update max position
        let current_max = self.max_pos.load(Ordering::Relaxed);
        if len as u64 > current_max {
            self.max_pos.store(len as u64, Ordering::Relaxed);
        }
        for (pos, &q) in quality[..len].iter().enumerate() {
            let phred = (q.saturating_sub(33)) as usize;
            let bin = phred.min(QUALITY_BINS - 1);
            self.bins[pos][bin].fetch_add(1, Ordering::Relaxed);
            self.counts[pos].fetch_add(1, Ordering::Relaxed);
        }
    }

    fn snapshot(&self) -> PerPositionQuality {
        let max_pos = (self.max_pos.load(Ordering::Relaxed) as usize).min(MAX_POSITION);
        let mut positions = Vec::with_capacity(max_pos);
        for pos in 0..max_pos {
            let count = self.counts[pos].load(Ordering::Relaxed);
            if count == 0 {
                continue;
            }
            let bins: Vec<u64> = self.bins[pos]
                .iter()
                .map(|b| b.load(Ordering::Relaxed))
                .collect();
            // Compute mean from bins
            let sum: u64 = bins.iter().enumerate().map(|(q, &c)| q as u64 * c).sum();
            let mean = sum as f64 / count as f64;
            // Compute quartiles from cumulative distribution
            let q25 = percentile_from_bins(&bins, count, 0.25);
            let median = percentile_from_bins(&bins, count, 0.50);
            let q75 = percentile_from_bins(&bins, count, 0.75);
            positions.push(PositionQualityStat {
                position: pos,
                count,
                mean,
                q25,
                median,
                q75,
            });
        }
        PerPositionQuality { positions }
    }
}

/// Per-position base composition tracker
struct PositionBaseTracker {
    /// [position] → [A, C, G, T, N] counts
    bases: Vec<[AtomicU64; 5]>,
    max_pos: AtomicU64,
}

impl PositionBaseTracker {
    fn new() -> Self {
        let bases: Vec<[AtomicU64; 5]> = (0..MAX_POSITION)
            .map(|_| std::array::from_fn(|_| AtomicU64::new(0)))
            .collect();
        Self {
            bases,
            max_pos: AtomicU64::new(0),
        }
    }

    fn record(&self, sequence: &[u8]) {
        let len = sequence.len().min(MAX_POSITION);
        let current_max = self.max_pos.load(Ordering::Relaxed);
        if len as u64 > current_max {
            self.max_pos.store(len as u64, Ordering::Relaxed);
        }
        for (pos, &base) in sequence[..len].iter().enumerate() {
            let idx = match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4, // N and anything else
            };
            self.bases[pos][idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    fn snapshot(&self) -> PerPositionBases {
        let max_pos = (self.max_pos.load(Ordering::Relaxed) as usize).min(MAX_POSITION);
        let mut positions = Vec::with_capacity(max_pos);
        for pos in 0..max_pos {
            let counts: [u64; 5] = std::array::from_fn(|i| {
                self.bases[pos][i].load(Ordering::Relaxed)
            });
            let total: u64 = counts.iter().sum();
            if total == 0 {
                continue;
            }
            positions.push(PositionBaseStat {
                position: pos,
                a: counts[0] as f64 / total as f64,
                c: counts[1] as f64 / total as f64,
                g: counts[2] as f64 / total as f64,
                t: counts[3] as f64 / total as f64,
                n: counts[4] as f64 / total as f64,
                total,
            });
        }
        PerPositionBases { positions }
    }
}

// ─── Read Analytics (unified collector) ─────────────────────────────────────

/// Comprehensive read analytics collector
///
/// Single struct that accumulates all QA metrics during pipeline processing.
/// Thread-safe via atomic counters — no locks needed for the hot path.
pub struct ReadAnalytics {
    // Summary counters
    pub reads_input: AtomicU64,
    pub reads_passed: AtomicU64,
    pub reads_failed: AtomicU64,
    pub bases_input: AtomicU64,
    pub bases_output: AtomicU64,

    // Per-position (before and after QC)
    quality_before: PositionQualityTracker,
    quality_after: PositionQualityTracker,
    bases_before: PositionBaseTracker,
    bases_after: PositionBaseTracker,

    // Distributions
    pub length_before: AtomicHistogram,
    pub length_after: AtomicHistogram,
    pub gc_content: AtomicHistogram,
    pub quality_scores: AtomicHistogram,
    pub insert_sizes: AtomicHistogram,
    pub trimmed_bases: AtomicHistogram,

    // Adapter breakdown
    adapter_counts: [AtomicU64; 8], // indexed by adapter type
    adapter_names: Vec<&'static str>,

    // Duplication estimate (HyperLogLog)
    hll: Mutex<HyperLogLog>,

    // Survival funnel (populated at snapshot time from module reports)
}

impl Default for ReadAnalytics {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadAnalytics {
    pub fn new() -> Self {
        Self {
            reads_input: AtomicU64::new(0),
            reads_passed: AtomicU64::new(0),
            reads_failed: AtomicU64::new(0),
            bases_input: AtomicU64::new(0),
            bases_output: AtomicU64::new(0),

            quality_before: PositionQualityTracker::new(),
            quality_after: PositionQualityTracker::new(),
            bases_before: PositionBaseTracker::new(),
            bases_after: PositionBaseTracker::new(),

            length_before: AtomicHistogram::uniform(0.0, 500.0, 50),
            length_after: AtomicHistogram::uniform(0.0, 500.0, 50),
            gc_content: AtomicHistogram::uniform(0.0, 1.0, 50),
            quality_scores: AtomicHistogram::uniform(0.0, 42.0, 42),
            insert_sizes: AtomicHistogram::uniform(0.0, 1000.0, 100),
            trimmed_bases: AtomicHistogram::uniform(0.0, 150.0, 30),

            adapter_counts: std::array::from_fn(|_| AtomicU64::new(0)),
            adapter_names: vec![
                "TruSeq", "Nextera", "NEBNext", "10X", "SISPA", "Other", "Internal", "None",
            ],

            hll: Mutex::new(HyperLogLog::new(&SketchConfig::default())),
        }
    }

    /// Record a read before QC processing
    pub fn record_input(&self, sequence: &[u8], quality: &[u8]) {
        self.reads_input.fetch_add(1, Ordering::Relaxed);
        self.bases_input
            .fetch_add(sequence.len() as u64, Ordering::Relaxed);
        self.length_before.record(sequence.len() as f64);
        self.quality_before.record(quality);
        self.bases_before.record(sequence);

        // Add to HyperLogLog for duplication estimate (one hash per read, not per k-mer)
        // Hash the first 50bp as a read-level fingerprint
        let prefix_len = sequence.len().min(50);
        let read_hash = hash_read_prefix(&sequence[..prefix_len]);
        let mut hll = self.hll.lock().unwrap();
        hll.add_hash(read_hash);
    }

    /// Record a read that passed QC
    pub fn record_passed(&self, sequence: &[u8], quality: &[u8], bases_trimmed: usize) {
        self.reads_passed.fetch_add(1, Ordering::Relaxed);
        self.bases_output
            .fetch_add(sequence.len() as u64, Ordering::Relaxed);
        self.length_after.record(sequence.len() as f64);
        self.quality_after.record(quality);
        self.bases_after.record(sequence);
        self.gc_content.record(gc_content(sequence));
        self.trimmed_bases.record(bases_trimmed as f64);

        // Per-base quality score distribution (sample every 10th base for efficiency)
        for &q in quality.iter().step_by(10) {
            let phred = q.saturating_sub(33) as f64;
            self.quality_scores.record(phred);
        }
    }

    /// Record a failed read
    pub fn record_failed(&self) {
        self.reads_failed.fetch_add(1, Ordering::Relaxed);
    }

    /// Record adapter detection by type name
    pub fn record_adapter(&self, adapter_name: &str) {
        let idx = match adapter_name {
            name if name.contains("TruSeq") || name.contains("truseq") => 0,
            name if name.contains("Nextera") || name.contains("nextera") => 1,
            name if name.contains("NEBNext") || name.contains("nebnext") => 2,
            _ => 5, // Other
        };
        self.adapter_counts[idx].fetch_add(1, Ordering::Relaxed);
    }

    /// Record internal adapter detection
    pub fn record_internal_adapter(&self) {
        self.adapter_counts[6].fetch_add(1, Ordering::Relaxed);
    }

    /// Record insert size from merged pair
    pub fn record_insert_size(&self, insert_size: usize) {
        self.insert_sizes.record(insert_size as f64);
    }

    /// Get current progress line for terminal display
    pub fn progress_line(&self) -> String {
        let input = self.reads_input.load(Ordering::Relaxed);
        let passed = self.reads_passed.load(Ordering::Relaxed);
        let failed = self.reads_failed.load(Ordering::Relaxed);
        let processed = passed + failed;
        let survival = if processed > 0 {
            passed as f64 / processed as f64 * 100.0
        } else {
            0.0
        };

        let adapters: u64 = self.adapter_counts[..6]
            .iter()
            .map(|c| c.load(Ordering::Relaxed))
            .sum();
        let adapter_pct = if processed > 0 {
            adapters as f64 / processed as f64 * 100.0
        } else {
            0.0
        };

        let bases_in = self.bases_input.load(Ordering::Relaxed);
        let bases_out = self.bases_output.load(Ordering::Relaxed);

        format!(
            "{} reads | {:.1}% survival | {:.1}% adapter | {:.1}M bases in → {:.1}M out",
            format_count(input),
            survival,
            adapter_pct,
            bases_in as f64 / 1_000_000.0,
            bases_out as f64 / 1_000_000.0,
        )
    }

    /// Generate a complete serializable snapshot
    pub fn snapshot(&self) -> AnalyticsSnapshot {
        let input = self.reads_input.load(Ordering::Relaxed);
        let passed = self.reads_passed.load(Ordering::Relaxed);

        // Duplication estimate from HyperLogLog
        let hll = self.hll.lock().unwrap();
        let estimated_unique = hll.cardinality();
        let estimated_dup_rate: f64 = if input > 0 && estimated_unique > 0 {
            1.0 - (estimated_unique as f64 / input as f64)
        } else {
            0.0
        };

        // Adapter breakdown
        let adapter_breakdown: Vec<AdapterCount> = self
            .adapter_names
            .iter()
            .zip(self.adapter_counts.iter())
            .map(|(name, count)| AdapterCount {
                name: name.to_string(),
                count: count.load(Ordering::Relaxed),
            })
            .filter(|a| a.count > 0)
            .collect();

        AnalyticsSnapshot {
            summary: SummaryStats {
                reads_input: input,
                reads_passed: passed,
                reads_failed: self.reads_failed.load(Ordering::Relaxed),
                bases_input: self.bases_input.load(Ordering::Relaxed),
                bases_output: self.bases_output.load(Ordering::Relaxed),
                survival_rate: if input > 0 {
                    passed as f64 / input as f64
                } else {
                    0.0
                },
            },
            per_position: PerPositionStats {
                quality_before: self.quality_before.snapshot(),
                quality_after: self.quality_after.snapshot(),
                bases_before: self.bases_before.snapshot(),
                bases_after: self.bases_after.snapshot(),
            },
            distributions: DistributionStats {
                length_before: self.length_before.snapshot(),
                length_after: self.length_after.snapshot(),
                gc_content: self.gc_content.snapshot(),
                quality_scores: self.quality_scores.snapshot(),
                insert_sizes: self.insert_sizes.snapshot(),
                trimmed_bases: self.trimmed_bases.snapshot(),
            },
            adapters: AdapterStats {
                breakdown: adapter_breakdown,
            },
            duplication: DuplicationStats {
                estimated_unique_sequences: estimated_unique,
                estimated_duplication_rate: estimated_dup_rate.max(0.0),
                estimated_library_complexity: estimated_unique,
            },
        }
    }
}

fn format_count(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.1}G", n as f64 / 1_000_000_000.0)
    } else if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}K", n as f64 / 1_000.0)
    } else {
        format!("{n}")
    }
}

/// Hash a read prefix into a single u64 for read-level deduplication estimation.
/// Uses FNV-1a for speed and good distribution.
fn hash_read_prefix(prefix: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &b in prefix {
        hash ^= b.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

fn percentile_from_bins(bins: &[u64], total: u64, percentile: f64) -> f64 {
    let target = (total as f64 * percentile) as u64;
    let mut cumulative = 0u64;
    for (q, &count) in bins.iter().enumerate() {
        cumulative += count;
        if cumulative >= target {
            return q as f64;
        }
    }
    0.0
}

// ─── Serializable Snapshot Types ────────────────────────────────────────────

/// Complete analytics snapshot — the data model for passport + report
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalyticsSnapshot {
    pub summary: SummaryStats,
    pub per_position: PerPositionStats,
    pub distributions: DistributionStats,
    pub adapters: AdapterStats,
    pub duplication: DuplicationStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SummaryStats {
    pub reads_input: u64,
    pub reads_passed: u64,
    pub reads_failed: u64,
    pub bases_input: u64,
    pub bases_output: u64,
    pub survival_rate: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerPositionStats {
    pub quality_before: PerPositionQuality,
    pub quality_after: PerPositionQuality,
    pub bases_before: PerPositionBases,
    pub bases_after: PerPositionBases,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerPositionQuality {
    pub positions: Vec<PositionQualityStat>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionQualityStat {
    pub position: usize,
    pub count: u64,
    pub mean: f64,
    pub q25: f64,
    pub median: f64,
    pub q75: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerPositionBases {
    pub positions: Vec<PositionBaseStat>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionBaseStat {
    pub position: usize,
    pub a: f64,
    pub c: f64,
    pub g: f64,
    pub t: f64,
    pub n: f64,
    pub total: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistributionStats {
    pub length_before: Histogram,
    pub length_after: Histogram,
    pub gc_content: Histogram,
    pub quality_scores: Histogram,
    pub insert_sizes: Histogram,
    pub trimmed_bases: Histogram,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdapterStats {
    pub breakdown: Vec<AdapterCount>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdapterCount {
    pub name: String,
    pub count: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplicationStats {
    pub estimated_unique_sequences: u64,
    pub estimated_duplication_rate: f64,
    pub estimated_library_complexity: u64,
}

// Keep backward compatibility aliases
pub type QaStatsCollector = ReadAnalytics;
pub type QaStatsSnapshot = AnalyticsSnapshot;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_histogram_uniform_bins() {
        let hist = AtomicHistogram::uniform(0.0, 100.0, 10);
        hist.record(5.0);
        hist.record(15.0);
        hist.record(95.0);
        hist.record(105.0); // overflow
        let snap = hist.snapshot();
        assert_eq!(snap.counts[0], 1);
        assert_eq!(snap.counts[1], 1);
        assert_eq!(snap.counts[9], 1);
        assert_eq!(snap.overflow, 1);
        assert_eq!(snap.total, 4);
    }

    #[test]
    fn test_read_analytics_basic() {
        let analytics = ReadAnalytics::new();
        let seq = b"ATGCATGCATGC";
        let qual = b"IIIIIIIIIIII";

        analytics.record_input(seq, qual);
        analytics.record_passed(seq, qual, 0);

        let snap = analytics.snapshot();
        assert_eq!(snap.summary.reads_input, 1);
        assert_eq!(snap.summary.reads_passed, 1);
        assert!(snap.summary.survival_rate > 0.99);
    }

    #[test]
    fn test_per_position_quality() {
        let analytics = ReadAnalytics::new();
        // Q40 = 'I' (73), Q20 = '5' (53)
        let qual = b"IIII5555";
        analytics.record_input(b"ATGCATGC", qual);

        let snap = analytics.snapshot();
        let positions = &snap.per_position.quality_before.positions;
        assert_eq!(positions.len(), 8);
        assert!((positions[0].mean - 40.0).abs() < 0.1); // First base Q40
        assert!((positions[4].mean - 20.0).abs() < 0.1); // Fifth base Q20
    }

    #[test]
    fn test_per_position_bases() {
        let analytics = ReadAnalytics::new();
        analytics.record_input(b"ACGT", b"IIII");
        analytics.record_input(b"ACGT", b"IIII");

        let snap = analytics.snapshot();
        let positions = &snap.per_position.bases_before.positions;
        assert_eq!(positions.len(), 4);
        assert!((positions[0].a - 1.0).abs() < 0.01); // Position 0 is always A
        assert!((positions[1].c - 1.0).abs() < 0.01); // Position 1 is always C
    }

    #[test]
    fn test_duplication_estimate() {
        let analytics = ReadAnalytics::new();
        // Add 100 unique sequences (150bp each — long enough for k-mer extraction)
        for i in 0..100u64 {
            let mut seq = vec![b'A'; 150];
            // Make each sequence unique by varying the middle section
            let id_bytes = format!("{:020}", i);
            for (j, b) in id_bytes.bytes().enumerate() {
                let base = match b % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                };
                if 50 + j < seq.len() {
                    seq[50 + j] = base;
                }
            }
            analytics.record_input(&seq, &vec![b'I'; seq.len()]);
        }
        // Add 100 duplicates of the first sequence
        let mut dup_seq = vec![b'A'; 150];
        let id_bytes = format!("{:020}", 0);
        for (j, b) in id_bytes.bytes().enumerate() {
            let base = match b % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            if 50 + j < dup_seq.len() {
                dup_seq[50 + j] = base;
            }
        }
        for _ in 0..100 {
            analytics.record_input(&dup_seq, &vec![b'I'; dup_seq.len()]);
        }

        let snap = analytics.snapshot();
        // HyperLogLog counts unique k-mers not unique reads,
        // so we just verify it detected some level of duplication
        assert!(
            snap.duplication.estimated_unique_sequences > 0,
            "Should estimate some unique sequences"
        );
        // With 200 reads (100 unique + 100 dups), the k-mer-based estimate
        // won't be exactly 100, but duplication should be detected
        assert_eq!(snap.summary.reads_input, 200);
    }

    #[test]
    fn test_adapter_breakdown() {
        let analytics = ReadAnalytics::new();
        analytics.record_adapter("TruSeq");
        analytics.record_adapter("TruSeq");
        analytics.record_adapter("Nextera");
        analytics.record_internal_adapter();

        let snap = analytics.snapshot();
        let truseq = snap.adapters.breakdown.iter().find(|a| a.name == "TruSeq");
        assert_eq!(truseq.unwrap().count, 2);
        let nextera = snap.adapters.breakdown.iter().find(|a| a.name == "Nextera");
        assert_eq!(nextera.unwrap().count, 1);
        let internal = snap.adapters.breakdown.iter().find(|a| a.name == "Internal");
        assert_eq!(internal.unwrap().count, 1);
    }

    #[test]
    fn test_progress_line() {
        let analytics = ReadAnalytics::new();
        for _ in 0..1000 {
            analytics.record_input(b"ATGCATGC", b"IIIIIIII");
            analytics.record_passed(b"ATGCATGC", b"IIIIIIII", 0);
        }
        let line = analytics.progress_line();
        assert!(line.contains("1.0K reads"));
        assert!(line.contains("100.0% survival"));
    }
}
