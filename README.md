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

- **6 QC modules** in biologically-motivated order: adapter trimming, poly-X removal, N-base filtering, quality trimming, complexity filtering, length filtering
- **Internal adapter detection** using precomputed k-mer hash index with 1-mismatch tolerance
- **Platform-aware poly-G** trimming with lower thresholds for two-color chemistry artifacts
- **Paired-end support** with singleton handling and optional read merging
- **Comprehensive analytics** including per-position quality, base composition, GC distribution, library complexity estimation (HyperLogLog), adapter breakdown, and insert size distribution
- **Profile system** with presets for common virome sample types
- **Gzip I/O** for both input and output
- **Parallel processing** via Rayon (160K+ reads/sec on Apple Silicon)
- **Synthetic test corpus generator** with ground truth labels for benchmarking

Built on [biometal](https://github.com/shandley/biometal) for NEON-optimized sequence operations.

## Installation

Requires Rust 1.75+ and a local copy of [biometal](https://github.com/shandley/biometal).

```bash
git clone https://github.com/shandley/virome-qc.git
cd virome-qc

# Update Cargo.toml to point biometal path to your local copy
# Then build:
cargo build --release
```

The binary will be at `target/release/virome-qc`.

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
| `stool-vlp-tagmentation` | VLP-enriched stool with Nextera/tagmentation library prep |
| `stool-vlp-truseq` | VLP-enriched stool with TruSeq library prep |
| `tissue-truseq` | Tissue/biopsy with TruSeq (expects high host fraction) |
| `metagenomics-nextera` | Shotgun metagenomics with Nextera |
| `low-biomass-wga` | Low-input samples with whole-genome amplification (MDA/SISPA) |
| `rna-virome-truseq` | RNA virome with TruSeq |

Custom profiles can be provided as YAML files. See `profiles/stool-vlp-tagmentation.yaml` for the format.

## QC pipeline

Modules run in this order:

1. **Adapter trimming** -- Removes 3' adapters by overlap detection (forward + RC for read-through), scans for internal adapter contamination via precomputed k-mer hash index, trims random primer bias from 5' end
2. **Poly-X trimming** -- Removes homopolymer runs at the 3' end. Platform-aware: aggressive poly-G threshold on 2-color chemistry (NextSeq/NovaSeq), classifies poly-X as likely artifact or likely genomic
3. **N-base filtering** -- Removes reads with excessive ambiguous bases. Threshold auto-set from ingestion N-rate.
4. **Quality trimming** -- Sliding window quality trim from 3' end, leading quality trim from 5' end, mean quality filter. Bin-aware mode for NovaSeq/NextSeq quantized Q-scores (auto-detected).
5. **Complexity filtering** -- Shannon entropy filter with data-driven threshold computed from the 2nd percentile of the sample's complexity distribution during ingestion.
6. **Contaminant screening** -- k-mer containment against rRNA (16S/23S/5S prokaryotic, 18S/28S/5.8S eukaryotic), PhiX174, and cloning vectors (pUC19). Reports prokaryotic vs eukaryotic rRNA breakdown. Paired-end concordant flagging (if one mate is contaminant, the other is removed too).
7. **Host depletion** -- Super Bloom k-mer containment screening against a pre-built host genome filter. Three-way classification: host (>50% containment, removed), ambiguous (15-50%, written to separate file), not host (<15%, kept). Requires building a filter first with `virome-qc db`.
8. **Length filtering** -- Final safety net for reads shortened by cumulative trimming across all modules.

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

Every run also generates an interactive `report.html` dashboard alongside the passport. Self-contained single file with Chart.js charts, Catppuccin theme, and dark mode toggle. Includes per-position quality profiles, base composition, read length distributions, GC content, survival funnel, contaminant breakdown, and adapter analysis. Can also be generated from an existing passport:

```bash
virome-qc report -i passport.json -o report.html
```

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

Benchmarked on Apple Silicon (M3 Max) with 1M reads (150bp):

| Threads | Wall time | Reads/sec |
|---------|-----------|-----------|
| 1 | 8s | 125K/sec |
| 4 | 6.6s | 152K/sec |
| 8 | 6.2s | 161K/sec |

All QC modules enabled including internal adapter k-mer scan, per-position analytics, and HyperLogLog library complexity estimation.

## Ingestion engine

Before QC runs, virome-qc scans the first 50,000 reads to auto-detect platform characteristics and configure the pipeline:

- **Platform detection**: Illumina (MiSeq, NextSeq, NovaSeq, HiSeq, iSeq), BGI/MGI, PacBio, ONT, Element AVITI. Handles ENA-reformatted and SRA-stripped headers.
- **Adapter auto-detection**: Identifies Nextera vs TruSeq/NEBNext from read content. Overrides profile if mismatch detected.
- **Data-driven parameters**: Complexity threshold from sample distribution, quality window from read length, min read length from actual reads, poly-G behavior from detected chemistry, N-filter threshold from baseline N-rate.
- **Q-score binning detection**: Automatically enables bin-aware quality trimming for NovaSeq/NextSeq.

No user configuration needed -- the ingestion engine adjusts profile parameters based on the data.

## Development status

Phases 1-3 are complete. Planned additions:

- **Phase 4**: Duplicate detection (optical + PCR) with UMI support
- **Phase 3b**: SNP-level endogenous vs exogenous retrovirus discrimination on ambiguous host reads
- **EVE annotation**: Flagging reads that map to known endogenous viral element regions

## License

MIT OR Apache-2.0
