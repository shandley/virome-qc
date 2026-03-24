# rRNA Screening Design

## Problem

RNA virome sequencing requires cDNA conversion, which captures all RNA species in the sample. rRNA dominates total RNA (80-90% by mass). Even with rRNA depletion (Ribo-Zero, RNase H), 5-30% residual rRNA is typical. For virome QC, accurate rRNA detection and removal is essential because:

1. rRNA reads consume assembly resources and inflate contig counts
2. rRNA fraction is the primary quality metric for RNA virome library prep
3. Without accurate rRNA quantification, researchers cannot assess depletion efficiency

## Current state

The existing contaminant screener uses 9 embedded consensus rRNA sequences (~11K 21-mers) from model organisms. Benchmarking on Zhang 2025 undepleted mouse cecal metatranscriptome (~87% rRNA) showed only 0.17% detection -- the consensus sequences cover a tiny fraction of real microbial rRNA diversity.

## Design: Super Bloom rRNA filter

Use the same Super Bloom architecture proven for host depletion, but against the SILVA rRNA database instead of a host genome.

### Reference databases

| Source | Sequences | Description |
|---|---|---|
| SILVA 138.2 SSU NR99 | 510,495 | 16S/18S, clustered at 99% identity |
| SILVA 138.2 LSU NR99 | 95,279 | 23S/28S, clustered at 99% identity |
| 5S/5.8S | ~1,000 | From Rfam/existing embedded sequences |

Total: ~606K sequences, covering all domains of life.

### Filter construction

```bash
# Build from SILVA (one-time, ~5-10 min)
virome-qc db --rrna

# Or provide custom FASTA
virome-qc db --rrna /path/to/silva_ssu_nr99.fasta --output rrna.sbf
```

The `db` subcommand already supports `--host`. Adding `--rrna` follows the same pattern.

### Estimated filter size

- SILVA SSU NR99: 510K sequences x ~1500bp avg = ~765M bases
- SILVA LSU NR99: 95K sequences x ~3000bp avg = ~285M bases
- Total: ~1.05 Gbp of reference sequence
- At k=21: ~1.05B k-mers (both strands), many shared across sequences
- Estimated unique 21-mers: ~100-200M (high conservation within rRNA)
- Super Bloom filter at 100-200M k-mers: ~1-2 GB

### Integration with contaminant module

Two approaches:

**Option A: Replace embedded consensus k-mers when filter available**

```rust
if config.rrna_filter_path.is_some() {
    // Use Super Bloom filter for rRNA screening (high sensitivity)
    let filter = load_rrna_filter(path)?;
    // Query reads against filter
} else {
    // Fall back to embedded consensus k-mers (basic screening)
    // Current behavior
}
```

**Option B: Separate rRNA module**

Add a dedicated `RrnaFilter` module (like `HostFilter`) that loads the Super Bloom filter and runs after the contaminant screener. The contaminant screener keeps its embedded k-mers for PhiX and vector detection.

Option B is cleaner -- separation of concerns. The contaminant module handles PhiX and vectors (small, embedded databases). The rRNA module handles rRNA (large, external database). Same pattern as host depletion.

### Module ordering (updated)

```
1. Adapter trim
2. Poly-X trim
3. N-filter
4. Quality trim
5. Complexity filter
6. Streaming dedup
7. Contaminant screening (PhiX, vectors only)
8. rRNA screening (Super Bloom, optional)
9. Host depletion (Super Bloom, optional)
10. Length filter
```

rRNA screening runs before host depletion because rRNA reads can have partial host homology (eukaryotic rRNA is host-derived). Screening rRNA first prevents these from being double-counted as host.

### Classification threshold

Same approach as host depletion: k-mer containment fraction.

- rRNA classified: containment >= 0.3 (30% of read k-mers match rRNA database)
- Threshold is lower than host (50%) because rRNA reads from diverse organisms may share fewer k-mers with the NR99 database than host reads share with the reference genome

### Passport report additions

```json
{
  "rrna_screened": 1000000,
  "rrna_removed": 870000,
  "rrna_fraction": 0.87,
  "rrna_prokaryotic": 650000,
  "rrna_eukaryotic": 220000
}
```

The rRNA fraction becomes the primary quality metric for RNA virome samples, displayed prominently in the report.

### Profile changes

```yaml
modules:
  rrna:
    enabled: true
    filter: rrna.sbf    # or "auto" to check standard locations
    min_kmer_fraction: 0.3
```

Enable by default in RNA virome profiles. Disable in DNA virome profiles (where rRNA contamination is minimal and the current consensus approach is adequate).

## Validation plan

1. Run ribodetector on Zhang undepleted sample (1M reads) as ground truth
2. Build Super Bloom rRNA filter from SILVA SSU+LSU NR99
3. Run our filter on same sample, compare sensitivity and specificity
4. Compare speed: our filter vs ribodetector vs SortMeRNA
5. Test on Zhang depleted sample (15% residual rRNA) for partial-depletion accuracy
6. Test on Cook DNA virome (expected <0.1% rRNA) for false positive rate

### Success criteria

- Sensitivity: detect >90% of reads classified as rRNA by ribodetector
- Specificity: <0.1% false positive rate on known non-rRNA reads
- Speed: faster than ribodetector (target: <1 min for 1M reads)
- Filter size: <2 GB

## Experiment 1: Super Bloom rRNA filter validation (2026-03-22)

### Setup

- Built Super Bloom filter from SILVA 138.2 SSU NR99 (510K seqs) + LSU NR99 (95K seqs) = 606K sequences
- U->T conversion applied during loading (SILVA uses RNA alphabet)
- Filter: 2.1 GB, k=31, built in 38 seconds
- Test data: Zhang 2025 undepleted mouse cecal metatranscriptome (SRR33419012), 1M read pairs (2M reads total)

### Ground truth: ribodetector

ribodetector (v0.3.3, CPU mode, BiLSTM neural network) classified:
- rRNA: 968,844 / 1,000,000 = **96.9%**
- Non-rRNA: 31,156 / 1,000,000 = 3.1%

### Super Bloom results (k=31, SILVA NR99)

| Containment threshold | Flagged | % of reads | Sensitivity vs ribodetector |
|---|---|---|---|
| >0.30 (removed) | 1,197,601 | 59.9% | 61.8% |
| >0.15 (removed + ambiguous) | 1,423,878 | 71.2% | 73.5% |
| Remaining unflagged | 576,122 | 28.8% | -- |

### Analysis

**71.2% sensitivity at permissive threshold vs 96.9% ground truth = 25.7 percentage point gap.**

Root causes:

1. **k=31 is too long for divergent rRNA.** Super Bloom uses 31-mers designed for host genomes (3.1 Gbp, low diversity). rRNA sequences are short (1500bp for 16S, 3000bp for 23S) and vary substantially across taxa. A single nucleotide difference breaks a 31-mer match. At 95% sequence identity (common between genera), a 150bp read would share only ~30% of its 31-mers with the closest database entry.

2. **NR99 clustering leaves coverage gaps.** SILVA NR99 collapses sequences at 99% identity. Reads from organisms at 97-99% identity to the nearest NR99 representative (within-species variation) lose proportionally more 31-mer matches than they would with k=21.

3. **Super Bloom's findere consistency scheme is optimized for k>=31.** The data structure uses (k=31, m=21, s=27) where m is the minimizer size. Reducing k to 21 would require m<=21, which changes the false positive characteristics of the filter.

### Comparison with consensus k-mer approach

| Method | rRNA detected | Sensitivity | Speed | Size |
|---|---|---|---|---|
| Consensus k-mers (9 sequences, k=21) | 0.17% | 0.2% | Instant | 0 (embedded) |
| Super Bloom SILVA (606K sequences, k=31) | 71.2% | 73.5% | ~30s/2M reads | 2.1 GB |
| ribodetector (ML, BiLSTM) | 96.9% | 100% (reference) | ~2min/2M reads | ~200 MB model |

### Improvement vs current: 420x better sensitivity

The Super Bloom approach is **420x more sensitive** than the consensus k-mer approach (71.2% vs 0.17%). While it doesn't match ribodetector's 96.9%, it catches the vast majority of rRNA reads and is faster.

### Path to higher sensitivity

Options to close the remaining 25% gap:

1. **Shorter k-mer (k=21)**: Would require either modifying Super Bloom's parameters or using a different data structure (FxHashSet of 21-mers, ~800M entries, ~12 GB RAM). The FxHashSet approach is simpler but memory-intensive.

2. **Lower containment threshold (0.05-0.10)**: Risk of false positives on non-rRNA reads. Need to validate specificity on a known non-rRNA dataset.

3. **Use full SILVA Ref (2.2M SSU sequences) instead of NR99**: 4x more sequences would improve coverage of rare taxa. Filter would be ~4 GB.

4. **Hybrid approach**: Super Bloom for primary screening + ribodetector for borderline reads. But this adds a Python dependency.

5. **Accept 71% sensitivity for virome QC**: For virome samples with 5-30% rRNA, catching 71% of rRNA reduces it to 1.5-8.7%. Combined with the existing consensus k-mer screener, this is a substantial improvement over the current 0.17%.

### Recommendation

Implement Option B (separate rRNA module) with the SILVA Super Bloom filter at k=31 and threshold 0.15. For virome QC, 71% sensitivity is a massive improvement over 0.17% and adequate for:
- Quantifying rRNA depletion efficiency (5% vs 30% rRNA is clearly distinguishable)
- Removing the bulk of rRNA contamination before assembly
- Flagging samples with poor enrichment or depletion

For users needing >95% rRNA sensitivity (metatranscriptomics), recommend ribodetector or SortMeRNA as a preprocessing step. Document this in the tool's output.

### Next steps

1. **Build a k=21 sorted Vec<u64> rRNA filter from SILVA** -- this is the preferred approach over Super Bloom k=31. Estimated ~100M unique hashes, ~800 MB on disk. Binary search gives O(log n) = 27 lookups per k-mer query, fast enough for screening.

2. **Implement as a dedicated rRNA module** (not overloading the host module). New `RrnaFilter` struct similar to `HostFilter` but using a sorted hash vector instead of Super Bloom. This avoids the k=31 constraint and should achieve >90% sensitivity.

3. **Add `virome-qc db --rrna` command** to build the filter from SILVA SSU+LSU NR99. Auto-download from SILVA like `db --host` auto-downloads T2T-CHM13.

4. **Validate against ribodetector** on the Zhang undepleted sample. Target: >90% sensitivity with <0.1% false positives.

5. **Test on Cook DNA virome** for false positive rate on non-rRNA data.

## Design: k=21 sorted hash vector

The Super Bloom approach at k=31 achieves 71.2% sensitivity. The k=21 approach should achieve much higher sensitivity because:

- At 95% sequence identity, a 21-mer has 0.95^21 = 34% chance of exact match
- At 97% identity: 0.97^21 = 52% chance
- At 99% identity: 0.99^21 = 81% chance

Compared to k=31:
- At 95% identity: 0.95^31 = 20% chance
- At 97% identity: 0.97^31 = 39% chance
- At 99% identity: 0.99^31 = 73% chance

k=21 roughly doubles the per-k-mer match probability for divergent sequences.

### Data structure

```rust
pub struct RrnaFilter {
    /// Sorted vector of FNV-1a hashes of all 21-mers from SILVA
    hashes: Vec<u64>,
    /// Classification threshold (fraction of read k-mers matching)
    min_kmer_fraction: f64,
    /// Stats
    stats: AtomicStats,
    rrna_removed: AtomicU64,
}

impl RrnaFilter {
    fn contains(&self, hash: u64) -> bool {
        self.hashes.binary_search(&hash).is_ok()
    }
}
```

Build process:
1. Parse SILVA FASTA (gzipped), convert U->T
2. Extract all 21-mers from both strands
3. Hash each with FNV-1a (same hash as contaminant screener)
4. Sort and deduplicate
5. Write as binary file (sorted u64 array)

Query process:
1. For each read, extract all 21-mers
2. Hash each with FNV-1a
3. Binary search in sorted array
4. If fraction of hits exceeds threshold, classify as rRNA

### Validated performance (2026-03-22)

Built and tested on Zhang undepleted mouse cecal metatranscriptome (1M PE reads = 2M total):

**Build:**
- Input: SILVA 138.2 SSU NR99 (510K seqs) + LSU NR99 (95K seqs) = 606K sequences
- 2.0 billion raw 21-mer hashes extracted (both strands)
- 165,227,304 unique hashes after sort + dedup
- Filter file: **1.3 GB** on disk
- Build time: **8 minutes** (Rust, single-threaded sort)

**Query:**
- 2M reads processed in **32 seconds** (62.5K reads/sec)
- 1.3 GB filter loaded in <2 seconds

**Sensitivity: 98.3% vs ribodetector ground truth**

| Method | rRNA detected | Sensitivity | Speed | Filter size |
|---|---|---|---|---|
| Consensus k-mers (9 seqs, k=21) | 0.17% | 0.2% | instant | 0 (embedded) |
| Super Bloom SILVA (k=31) | 71.2% | 73.5% | 30s/2M reads | 2.1 GB |
| **Sorted hash SILVA (k=21)** | **95.2%** | **98.3%** | **32s/2M reads** | **1.3 GB** |
| ribodetector (ML, ground truth) | 96.9% | 100% | ~2min/2M reads | ~200 MB |

The k=21 approach achieves near-ribodetector sensitivity (98.3%) while being faster and requiring no Python/ML dependencies. The 1.7% gap (95.2% vs 96.9%) likely comes from highly divergent rRNA sequences not represented in SILVA NR99 at k=21 resolution.

### Implementation status

- `src/modules/rrna.rs`: Complete -- RrnaFilter QcModule with sorted hash vector
- `src/config/profiles.rs`: RrnaConfig added, `rrna` field on ModuleConfigs (Optional)
- `src/pipeline/executor.rs`: Integrated between contaminant screening and host depletion
- `src/main.rs`: `virome-qc db --rrna` command implemented
- Tests: 4 unit tests (hash determinism, case insensitivity, containment, round-trip)
- 114 tests total, all passing
