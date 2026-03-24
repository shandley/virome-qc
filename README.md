# virome-qc

High-performance, virome-specific quality control for sequencing data. Built for the Human Virome Program (HVP).

## Why virome-qc?

Standard QC tools (fastp, Trimmomatic, BBDuk) are designed for general-purpose metagenomics. Virome sequencing has unique challenges that these tools do not address:

- **Adapter contamination in reference databases** causes false-positive viral detections when even one read retains a partial adapter
- **Internal adapter contamination** from chimeric library fragments is not detected by standard 3' trimming
- **Poly-G artifacts** from NovaSeq/NextSeq two-color chemistry have high quality scores and evade Q-score filtering
- **Endogenous viral elements (EVEs)** in host genomes cause over-aggressive host read removal
- **Low-complexity viral sequences** (AT-rich phage, polyA from cDNA) are incorrectly filtered by generic complexity thresholds
- **Random primer bias** from SISPA/WGA preps requires targeted head-trimming

virome-qc addresses all of these with a single, opinionated pipeline that produces clean reads and a comprehensive analytics report.

## Features

- **11 QC modules** in biologically-motivated order: adapter trimming, poly-X removal, N-base filtering, quality trimming, complexity filtering, streaming deduplication, contaminant screening, rRNA screening (SILVA), host depletion, ERV analysis, length filtering
- **Internal adapter detection** using precomputed k-mer hash index with 1-mismatch tolerance
- **Platform-aware poly-G** trimming with lower thresholds for two-color chemistry artifacts
- **Contaminant screening** for PhiX spike-in (with Microviridae exclusion) and cloning vectors
- **rRNA screening** via SILVA-based sorted hash filter (165M k-mers, 78% sensitivity on diverse rRNA)
- **Host depletion** via Super Bloom k-mer containment with three-way classification (host/ambiguous/keep)
- **Paired-end support** with concordant mate handling, singleton output, and optional read merging
- **Comprehensive analytics** including per-position quality, base composition, GC distribution, library complexity estimation (HyperLogLog), adapter breakdown, and insert size distribution
- **Interactive HTML reports** using React + shadcn/ui + Recharts with Catppuccin theme and dark mode
- **Ingestion engine** that auto-detects platform, adapters, Q-score binning, per-position quality profile, and R2 quality degradation from the first 50K reads
- **Profile system** with presets for common virome sample types and data-driven parameter overrides
- **GC deviation detection** flags samples outside expected GC range for the profile
- **Gzip I/O** for both input and output
- **Parallel processing** via Rayon (160K+ reads/sec on Apple Silicon)
- **Synthetic test corpus generator** with ground truth labels for benchmarking
- **ERV analysis** classifies retroviral reads as endogenous (polymorphic HERV) or exogenous (active infection) using three-signal classifier (CpG depletion + MinHash + ORF integrity)
- **Validated on 27 datasets** (12 real public + 12 ViroForge reference + 3 ViroForge benchmarks) across 6 sequencing platforms, 5 library preps, and 5 sample types
- **Head-to-head with fastp**: virome-qc keeps more viral reads on clean data (+4.3%) and removes more contaminants on dirty data (-6.8%)

Built on [biometal](https://github.com/shandley/biometal) for NEON-optimized sequence operations.

## Installation

Requires Rust 1.75+ and local copies of [biometal](https://github.com/shandley/biometal) and [SuperBloom](https://github.com/shandley/SuperBloom).

```bash
# Clone all three repositories
git clone https://github.com/shandley/biometal.git
git clone https://github.com/shandley/SuperBloom.git
git clone https://github.com/shandley/virome-qc.git

# Update Cargo.toml paths to point to your local copies
cd virome-qc
# Edit biometal and superbloom paths in Cargo.toml

# Build
cargo build --release
```

The binary will be at `target/release/virome-qc`.

## Quickstart

**5 minutes from install to first report:**

```bash
# 1. Build reference databases (one-time, ~30 min)
#    Downloads SILVA rRNA + T2T-CHM13 human genome and builds filters
virome-qc db --setup

# 2. Run QC on your virome FASTQ files
virome-qc run \
  -p stool-vlp-tagmentation \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -i . --merge \
  -o results/ \
  -t 8

# 3. View the report
open results/report.html
```

That's it. The ingestion engine auto-detects your platform, adapters, and quality profile. All 11 QC modules run automatically. The output includes:

- `results/clean_R1.fastq.gz` / `clean_R2.fastq.gz` — QC-passed reads
- `results/merged.fastq.gz` — merged overlapping pairs
- `results/passport.json` — comprehensive QC analytics
- `results/report.html` — interactive HTML dashboard

Use `--report-only` to generate just the passport and report without writing clean FASTQ output (useful for QC assessment without disk overhead).

## Quick start

### Single-end

```bash
virome-qc run \
  -p stool-vlp-tagmentation \
  -i sample.fastq.gz \
  -o clean/ \
  -t 8
```

### Paired-end

```bash
virome-qc run \
  -p stool-vlp-tagmentation \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -i . \
  -o clean/ \
  --merge \
  -t 8
```

### Output

Single-end produces:

```
clean/
  clean_sample.fastq.gz    # QC-passed reads
  passport.json             # Analytics report
```

Paired-end with `--merge` produces:

```
clean/
  clean_R1.fastq.gz        # Paired reads (both mates passed, not merged)
  clean_R2.fastq.gz        # Paired reads (both mates passed, not merged)
  merged.fastq.gz          # Merged overlapping pairs
  singletons.fastq.gz      # Reads where one mate failed QC
  passport.json             # Analytics report
```

## Profiles

Profiles configure the QC pipeline for specific sample and library prep combinations. Each profile sets module parameters, adapter sequences, and quality thresholds appropriate for the expected data characteristics.

```bash
virome-qc profiles
```

Available profiles:

| Profile | Description |
|---------|-------------|
| `stool-vlp-tagmentation` | VLP-enriched stool virome (Nextera/tagmentation) |
| `tissue-truseq` | Tissue/biopsy (expects high host fraction) |
| `metagenomics-nextera` | Shotgun metagenomics |
| `low-biomass-wga` | Low-input samples with whole-genome amplification (MDA/SISPA) |

Custom profiles can be provided as YAML files. See `profiles/stool-vlp-tagmentation.yaml` for the format.

## QC pipeline

Modules run in this order:

1. **Adapter trimming** -- Removes 3' adapters by overlap detection (forward + RC for read-through), scans for internal adapter contamination via precomputed k-mer hash index, trims random primer bias from 5' end
2. **Poly-X trimming** -- Removes homopolymer runs at the 3' end. Platform-aware: aggressive poly-G threshold on 2-color chemistry (NextSeq/NovaSeq), classifies poly-X as likely artifact or likely genomic
3. **N-base filtering** -- Removes reads with excessive ambiguous bases. Threshold auto-set from ingestion N-rate.
4. **Quality trimming** -- Sliding window quality trim from 3' end, leading quality trim from 5' end, mean quality filter. Bin-aware mode for NovaSeq/NextSeq quantized Q-scores (auto-detected).
5. **Complexity filtering** -- Shannon entropy filter with data-driven threshold computed from the 2nd percentile of the sample's complexity distribution during ingestion.
6. **Streaming deduplication** -- Prefix hash deduplication using FxHashSet. Skips first 5bp (handles differential trimming) and uses 50bp prefix hash. Paired-end requires both mates to match. Disabled by default in single-end profiles.
7. **Contaminant screening** -- k-mer containment against PhiX174 (with Microviridae exclusion to prevent gut phage false positives) and cloning vectors (pUC19). Paired-end concordant flagging (if one mate is contaminant, the other is removed too).
8. **rRNA screening** -- SILVA-based sorted hash filter (165M FNV-1a hashes at k=21) for high-sensitivity rRNA detection across prokaryotic and eukaryotic species. Falls back to consensus k-mer screener when SILVA filter is not available. See [rRNA screening](#rrna-screening) below.
9. **Host depletion** -- Super Bloom k-mer containment screening against a pre-built host genome filter (T2T-CHM13). Three-way classification: host (>50% containment, removed), ambiguous (20-50%, written to separate file), not host (<20%, kept). Containment distribution reported in passport for transparency. Zero viral false positive rate at 50% threshold (verified with ground truth). Enabled by default.
10. **ERV analysis** -- Post-pipeline retroviral k-mer screen (101 references: 46 Retroviridae + 55 HERV consensus) with herpesvirus k-mer exclusion to prevent false positives. Classifies clusters as endogenous (polymorphic ERV, CpG-depleted) or exogenous (active infection, intact CpG) using three signals: ORF integrity, CpG depletion ratio, and MinHash distance. Reports per-locus classification in the passport. First virome QC tool with automated ERV vs exogenous retrovirus classification.
11. **Length filtering** -- Final safety net for reads shortened by cumulative trimming across all modules.

The ordering is intentional. Adapters must be removed first (non-biological sequence). Poly-X runs before quality trimming because NovaSeq poly-G artifacts have high quality scores and would survive Q-score-based trimming. Contaminant screening and host depletion run after all read cleaning to prevent adapter/quality artifacts from causing false classifications. Length filtering is last as a safety net.

## Analytics passport

Every run produces a `passport.json` containing:

- **Summary**: reads in/out, survival rate, bases in/out
- **Per-module reports**: reads processed, removed, modified, and bases removed by each module
- **Per-position quality**: mean, median, Q25, Q75 at each read position (before and after QC)
- **Per-position base composition**: A/C/G/T/N fractions at each position (before and after QC)
- **Distributions**: read length, GC content, quality scores, insert sizes, trimmed bases
- **Adapter breakdown**: counts by adapter type (TruSeq, Nextera, NEBNext, internal)
- **Duplication estimate**: library complexity from HyperLogLog unique k-mer counting
- **Quality tier**: PASS, WARN, or FAIL based on profile-specific thresholds
- **Flags**: specific warnings (low survival, high adapter contamination, etc.)

The passport is designed for both human review and programmatic consumption. Downstream tools can parse it to make decisions about sample inclusion or flag batch effects.

## HTML report

Every run generates an interactive `report.html` dashboard alongside the passport. Self-contained single file (React + shadcn/ui + Recharts, inlined) with Catppuccin theme and dark mode toggle. Includes:

- Summary stat cards: reads in/out, bases, survival, mean read length, mean quality, GC, library complexity, duplication, pairs/singletons/merge rate
- Quality flags with severity (PASS/WARN/FAIL)
- Survival funnel table with per-module breakdown
- Adapter and contaminant cards with prokaryotic/eukaryotic rRNA breakdown
- Per-position quality profiles with IQR band (before/after QC)
- Stacked base composition charts (before/after QC)
- Before/after read length overlay histogram
- GC content, quality score, and insert size distributions

Reports can also be generated from existing passport files:

```bash
virome-qc report -i passport.json -o report.html
```

Use `--report-only` with `virome-qc run` to generate only the passport and report without writing clean FASTQ output (useful for QC assessment without disk overhead).

## Host depletion

Host reads are removed using a Super Bloom k-mer containment filter. Build the filter once from any FASTA reference:

```bash
# Build filter (~5 min for human genome, 4 GiB output)
virome-qc db --host T2T-CHM13.fa -o human.sbf

# Use in a profile (set reference to the .sbf file path)
# host:
#   enabled: true
#   reference: /path/to/human.sbf
```

The filter uses the Super Bloom data structure (near-zero false positives via findere consistency scheme). Reads are classified by k-mer containment fraction: host (>50%, removed), ambiguous (15-50%, written to `ambiguous_*.fastq.gz`), not host (<15%, kept).

## Test corpus generator

virome-qc includes a synthetic FASTQ generator for benchmarking and development. Generated reads contain planted artifacts with ground truth labels in the FASTQ headers.

```bash
# Single-end corpus
virome-qc corpus -o test_corpus.fastq -s stool-vlp -n 100000

# Paired-end corpus
virome-qc corpus -o test_pe -s stool-vlp -n 100000 --paired
```

Sample types: `stool-vlp`, `tissue`, `low-biomass-wga`

Each read header encodes its ground truth:

```
@read_000042/1 source=viral;virus=crAssphage;adapter_3prime=truseq_r1:15;complexity=0.923
```

A summary JSON is written alongside the FASTQ with exact counts of each artifact type and source category.

## Performance

Benchmarked on Apple Silicon (M3 Max), all QC modules enabled including SILVA rRNA filter and T2T host depletion:

| Dataset | Reads | Wall time | Reads/sec |
|---------|-------|-----------|-----------|
| Zhang depleted (2M PE) | 2.0M | 43s | 47K/sec |
| Cook (2.5M PE) | 2.5M | 101s | 25K/sec |
| Shkoporov (12.6M PE) | 12.6M | 375s | 34K/sec |

Performance includes loading the SILVA rRNA filter (1.3 GB, 165M k-mers) and T2T-CHM13 host filter (4.1 GB). Without these reference databases, throughput is ~160K reads/sec.

### Comparison with fastp

| Metric | fastp 1.3.0 | virome-qc |
|--------|-------------|-----------|
| Speed | 3-10x faster | Slower (reference database loading) |
| Clean VLP virome (Shkoporov) | 93.3% survival | **97.6% survival** (+4.3% more viral reads kept) |
| Contaminated virome (Zhang) | 79.8% survival | **73.0% survival** (-6.8% more contaminants removed) |
| rRNA/host/PhiX detection | No | Yes |
| ERV classification | No | Yes |
| Data-driven thresholds | No | Yes (ingestion engine) |

On clean VLP data, fastp over-filters with fixed quality thresholds. On contaminated data, fastp passes through rRNA, host DNA, and PhiX that virome-qc removes. virome-qc's ingestion engine adapts thresholds to the data, producing more appropriate filtering for each sample.

## rRNA screening

virome-qc includes two levels of rRNA screening:

1. **Consensus rRNA screener** (built-in): k-mer matching against consensus rRNA sequences from 23 model organisms (~11,000 21-mers). Fast but limited sensitivity (~45% on diverse samples). Used as a fallback when the SILVA filter is not available.

2. **SILVA rRNA filter** (recommended): Sorted hash vector of 165M FNV-1a hashes at k=21, built from SILVA SSU+LSU NR99 database. Binary search for O(log n) per k-mer lookup. 78% sensitivity on diverse rRNA, validated against ViroForge ground truth. Build once with `virome-qc db --rrna`.

When the SILVA filter is available, it automatically replaces the consensus screener. The SILVA filter is enabled by default in all profiles and auto-discovers the filter file from standard locations.

**Validation**: On Zhang undepleted stool (95% rRNA), the SILVA filter correctly identified 95.2% of reads as rRNA. On ViroForge reference datasets with known 10.9% rRNA injection, virome-qc achieved 78% sensitivity. On VLP-enriched samples, rRNA is typically <1% and the filter catches the majority.

## Ingestion engine

Before QC runs, virome-qc scans the first 50,000 reads (R1 and R2 if paired) to auto-detect platform characteristics and configure the pipeline:

- **Platform detection**: Illumina (MiSeq, NextSeq, NovaSeq, HiSeq, iSeq), BGI/MGI (DNBSEQ-G400, DNBSEQ-T7, MGISEQ-2000), PacBio, ONT, Element AVITI. Handles ENA-reformatted and SRA-stripped headers.
- **Adapter auto-detection**: Identifies Nextera vs TruSeq/NEBNext from read content, with separate 3' tail and internal chimera rate tracking. Auto-enables internal adapter scanning when chimeras are detected.
- **Per-position quality profile**: Tracks quality decay across read positions to detect 3' degradation patterns (e.g., MiSeq 2x250).
- **R2 quality comparison**: Scans R2 separately and warns if quality is substantially lower than R1, indicating expected higher singleton rates.
- **GC validation**: Checks mean GC against the profile's expected range and warns if outside bounds (e.g., ocean virome at 38% GC vs gut virome at 43%).
- **Data-driven parameters**: Complexity threshold from sample distribution, quality window from read length, min read length from actual reads, poly-G behavior from detected chemistry, N-filter threshold from baseline N-rate.
- **Q-score binning detection**: Automatically enables bin-aware quality trimming for NovaSeq/NextSeq.

No user configuration needed -- the ingestion engine adjusts profile parameters based on the data.

## Validation

Validated on 27 datasets across 12 real public datasets (6 sequencing platforms, 5 library preps, 5 sample types) and 15 ViroForge synthetic datasets (12 reference atlas + 3 parameter sweep benchmarks). See [VALIDATION.md](VALIDATION.md) for complete results.

Key findings:
- **Zero viral false positive rate** in host depletion at 50% containment threshold (ground truth verified)
- **ERV-host correlation**: Pearson r=0.871 (p=0.001) — retroviral content tracks host background
- **Herpesvirus cross-reactivity**: identified and fixed via k-mer exclusion (99.3% FP reduction)
- **Giant virus assessment**: systematic scan of 40 virus families, only Orthoherpesviridae significant
- **fastp comparison**: virome-qc preserves more viral reads on clean data, removes more contaminants on dirty data
- **21 bugs found** through bidirectional ViroForge validation (11 virome-qc + 10 ViroForge)

## Development status

The tool is functional for production virome QC workflows. 133 tests, 24.8K lines Rust, clippy clean.

Planned additions:

- **Long-read support**: ONT and PacBio QC modules (adapter trimming, quality filtering for long-read chemistry)
- **Expected range reporting**: Embed ViroForge reference ranges into profiles, show user metrics vs expected range in reports
- **Population ERV database**: Known polymorphic HERV insertion sites from the HPRC pangenome
- **Coordinate-based EVE flagging**: BED annotation of known high-identity integration sites (ciHHV-6, recent HERV-K)

## License

MIT OR Apache-2.0
