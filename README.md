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

1. **Adapter trimming** -- Removes 3' adapters by overlap detection, scans for internal adapter contamination via k-mer hash index, trims random primer bias from 5' end
2. **Poly-X trimming** -- Removes homopolymer runs at the 3' end. Platform-aware: uses a shorter detection threshold for poly-G on NovaSeq/NextSeq
3. **N-base filtering** -- Removes reads with more than 10% ambiguous bases
4. **Quality trimming** -- Sliding window quality trim from 3' end, leading quality trim from 5' end, mean quality filter
5. **Complexity filtering** -- Shannon entropy filter with virome-aware thresholds (lower than typical to retain AT-rich phage and cDNA-derived reads)
6. **Length filtering** -- Removes reads shortened below the minimum length by cumulative trimming

The ordering is intentional. Adapters must be removed first (non-biological sequence). Poly-X runs before quality trimming because NovaSeq poly-G artifacts have high quality scores and would survive Q-score-based trimming. Complexity assessment runs on fully cleaned sequence. Length filtering is last as a safety net.

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

## Development status

Phase 1 is complete. Planned additions:

- **Phase 2**: Contaminant screening (rRNA, PhiX, vectors, kitome) via unified k-mer index
- **Phase 3**: Host depletion with EVE-aware masking and viral read rescue
- **Phase 4**: Duplicate detection (optical + PCR) with UMI support
- **Phase 5**: Interactive HTML report generation from passport data

## License

MIT OR Apache-2.0
