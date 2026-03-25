# Publication Plan

## Target

Single comprehensive methods paper for Genome Biology, Genome Research, or Nature Methods.

## Title (tentative)

"Data-driven virome quality control: a platform for automated QC assessment, contaminant discrimination, and community benchmarking"

## Narrative arc

### 1. The problem

Virome sequencing spans tissue biopsies (high host) to ocean water (novel diversity), DNA viromes to RNA viromes (cDNA/rRNA contamination), VLP-enriched to bulk metagenomics. Each produces different contamination profiles, artifact rates, and quality expectations. Existing QC tools (fastp, BBDuk, Trimmomatic) are generic -- they don't distinguish lab contaminants from natural viral relatives (PhiX vs Microviridae), don't quantify rRNA contamination, and don't provide virome-aware quality assessment. No framework exists for defining what "good" virome data looks like across sample types.

### 2. The tool (virome-qc)

11-module Rust pipeline (24.8K lines, 133 tests) with ingestion engine that auto-tunes all processing parameters from the first 50K reads. Key innovations:
- PhiX/Microviridae unique k-mer discrimination (first quantitative solution)
- SILVA rRNA screening at 78% sensitivity via sorted hash filter (165M 21-mers)
- Super Bloom host depletion with three-way classification and containment distribution reporting
- ERV analysis: three-signal classifier (CpG + MinHash + ORF) for endogenous vs exogenous retrovirus discrimination (first in any QC tool), with herpesvirus k-mer exclusion
- Data-driven parameter setting (complexity, N-filter, quality, adapter overlap)
- Platform-aware poly-G, bin-aware quality trimming, R2 quality comparison
- Self-contained React HTML reports with interactive charts, ERV analysis panel, containment histograms
- Contamination summary separating biological (host + rRNA + contaminant) from technical (adapter + quality) removal

### 3. The methodology (ViroForge loop)

ViroForge synthetic ground truth enables quantitative optimization impossible with real data. Demonstrated closed loop:
- Synthetic data exposed PhiX cross-reactivity (93% FP rate on Microviridae)
- Exposed internal adapter double-counting (54K reads destroyed per run)
- Exposed truncated embedded reference (PhiX at 2,237bp not 5,386bp)
- Exposed adapter k-mer homology with viral genomes (1.72% baseline FP)
- Each bug fixed with measured improvement, validated on real data
- Tool testing simultaneously improved ViroForge (curation errors, synthetic sequences)

### 4. The validation

27 datasets: 12 real public, 12 ViroForge reference atlas, 3 ViroForge parameter sweep benchmarks:
- Zero viral false positive rate in host depletion (ground truth verified with synthetic corpus)
- ERV-host correlation: Pearson r=0.871 (p=0.001) across 11 real datasets
- Herpesvirus cross-reactivity: systematic scan of 40 virus families, only Orthoherpesviridae significant (194 shared k-mers, fixed with exclusion set)
- Giant virus assessment: Phycodnaviridae (12-14 shared k-mers) documented as marginal risk for environmental/stool viromes
- Santos KAPA vs Nextera: 580x rRNA difference from same sample — library prep quality detection
- TARA ocean: 0.18% host FP from giant virus eukaryotic gene homology — documented as known behavior
- Module ordering: rRNA before host creates funnel effect where host-origin rRNA reads are caught by SILVA filter first (97% true host sensitivity at read level, verified with containment diagnostic)

### 4b. fastp head-to-head

Direct comparison on 3 representative datasets (fastp 1.3.0, default parameters, 8 threads):
- **Clean VLP (Shkoporov)**: virome-qc 97.6% vs fastp 93.3% — fastp over-filters 704K reads as low quality
- **Contaminated (Zhang depleted)**: virome-qc 73.0% vs fastp 79.8% — fastp passes through 13.8% rRNA
- **MiSeq degraded (Cook)**: virome-qc 91.7% vs fastp 97.1% — virome-qc removes rRNA + host that fastp cannot detect
- fastp is 3-10x faster but misses all biological contaminants (rRNA, host, PhiX, ERVs)

### 5. The atlas (virome-qc-atlas)

Reference distributions derived from ViroForge synthetic viromes (20 collections) calibrated with real data. Auto-classification: "your sample most closely matches X, with deviations in Y and Z." Public repository of passport JSON files as community resource for defining virome data quality standards.

## Figures

1. Pipeline architecture with ingestion engine and 11-module flow
2. PhiX/Microviridae and herpesvirus/retrovirus k-mer exclusion analysis
3. Cross-dataset comparison (12 real + 12 ViroForge reference, heatmap or parallel coordinates)
4. ViroForge validation loop (21 bugs found, before/after metrics)
5. fastp comparison: survival rates across clean/contaminated/degraded datasets
6. ERV-host correlation scatter plot (r=0.871) with containment distributions
7. Reference atlas expected ranges by sample type (from ViroForge datasets)

## Supplementary

- Complete per-dataset validation results (VALIDATION.md content)
- ViroForge adapter sweep results (14 datasets, sensitivity/specificity)
- rRNA filter construction and ribodetector comparison
- Host depletion threshold analysis
- Module ordering justification

### 6. Future directions

**Clinical diagnostics mode**: The ERV endogenous/exogenous classifier is already a diagnostic-grade feature — it answers "is this retroviral signal from the patient's germline or an active infection?" A `--clinical` flag would tighten thresholds across all modules (stricter quality, mandatory dedup, conservative merging, expanded contaminant database) and add negative control comparison (`--blank`). Key additions: negative control subtraction, minimum evidence thresholds for pathogen calls, and confidence scoring for ERV classification. Note: rRNA as an extraction validation positive control is only applicable for non-sterile sites (stool, oral, respiratory) — sterile sites (blood, CSF, tissue) would not have rRNA.

**Probe-capture virome support**: Probe-capture/target-enrichment data requires fundamentally different QC metrics (on-target rate, enrichment fold, coverage uniformity) that are alignment-based rather than k-mer-based. This is a separate modality requiring panel BED files and a distinct reporting framework, best addressed as a v2 expansion or companion tool.

**Long-read QC**: ONT and PacBio virome data requires different adapter trimming, quality filtering, and error correction strategies.

## Authors

Scott Handley and collaborators

## Timeline

1. ~~Finalize module review and ViroForge optimization~~ COMPLETE
2. ~~Generate ViroForge reference datasets (12 sample types)~~ COMPLETE
3. ~~Run fastp comparison on 3 datasets~~ COMPLETE
4. ~~Embed ViroForge reference ranges into profiles~~ COMPLETE
5. ~~Final Tier 1 benchmark (13 datasets, HTCF)~~ COMPLETE
6. ~~Separate dedup from QC survival + conservative merge~~ COMPLETE
7. Generate publication figures
8. Write manuscript
9. Code/data deposition (GitHub, Zenodo for filters)

## Code and data availability

- virome-qc: https://github.com/shandley/virome-qc
- ViroForge: https://github.com/shandley/viroforge
- virome-qc-atlas: https://github.com/shandley/virome-qc-atlas (public passport repository)
- SILVA rRNA filter and T2T host filter: Zenodo
- All benchmark datasets: ENA/SRA accessions listed in VALIDATION.md
