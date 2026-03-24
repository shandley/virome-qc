# Benchmark Dataset Selection

Updated 2026-03-24. Organized into tiers for paper, extended validation, and HPC atlas generation.

## Design Principles

- Full runs (no `-X` truncation) — truncated data hides platform artifacts
- Paired-end where possible — PE is the standard for virome analysis
- One dataset per niche, chosen for being the best-characterized example
- Include ground-truth datasets (mocks, known composition)
- Cover every major virome sample type and platform
- Scalable: core set for the paper, extended set for HPC atlas generation

## Tier 1: Core Set (~15 datasets, for paper)

### Ground truth / Mocks

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 1 | Cook 15-phage mock | PRJEB56639 (ERR10359658) | MiSeq 2x250 | 2.5M PE | Known composition, adapter, quality, merge |
| 2 | 72-virus synthetic community | JoV 2024 | NovaSeq 2x150 | TBD | Most complex mock (21 families, all genome types) |

### Gut VLP (core use case)

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 3 | Shkoporov gut VLP | SRR9161520 | Illumina 2x151 | 12.6M PE | Clean gut VLP baseline |
| 4 | Santos KAPA | SRR8487022 | Illumina 2x151 | 27.4M PE | Clean VLP, PhiX spike-in detection |
| 5 | Santos Nextera | SRR8487034 | Illumina 2x151 | 26.8M PE | Library prep comparison (580x rRNA difference) |

### Clinical / High-host

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 6 | Buddle WGS (FULL) | ERR13480651 | NovaSeq SE 150 | ~10M | Host depletion, internal adapter chimeras |
| 7 | Buddle RNA (FULL) | ERR13480663 | NextSeq SE 150 | ~10M | Clinical RNA, host + adapter |

### rRNA / Depletion comparison

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 8 | Zhang undepleted | SRR33419012 | NovaSeq 2x150 | 2M PE | Extreme rRNA (95%), SILVA validation |
| 9 | Zhang depleted | SRR33419066 | NovaSeq 2x150 | 2M PE | Incomplete depletion, poly-G |

### Environmental

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 10 | TARA ocean | ERR599370 | Illumina 2x101 | 19.5M PE | Giant virus FPs, environmental baseline |
| 11 | Tisza wastewater | PRJNA966185 (1 sample) | NovaSeq 2x150 | TBD | Best-characterized wastewater virome |

### Non-Illumina platforms

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 12 | Chrisman DNBSEQ (FULL) | ERR9765742 | DNBSEQ 100 | ~100M | BGI/MGI platform |
| 13 | Cook ONT | PRJEB56639 | MinION | TBD | Long-read virome QC (same mock as #1) |

### Special cases

| # | Dataset | Accession | Platform | Reads | What it tests |
|---|---------|-----------|----------|-------|--------------|
| 14 | Negativeome blanks | PRJEB33578 | Illumina | TBD | Extraction controls, contaminant baseline |
| 15 | HIV+ gut virome | PRJNA1045584 (1 sample) | NextSeq 2x150 | TBD | Active retrovirus for ERV module |

## Tier 2: Extended Set (for HPC, +10 datasets)

| # | Dataset | Accession | Platform | What it tests |
|---|---------|-----------|----------|--------------|
| 16 | Cook PacBio | PRJEB56639 | PacBio Sequel | Third platform, same mock |
| 17 | Zolfo ONT gut virome | PRJEB47625 | ONT + Illumina | Paired long/short-read |
| 18 | Santos-Medellin soil VLP | PRJNA646773 | HiSeq 2x150 | Soil virome |
| 19 | Skin virome | PRJNA754140 | HiSeq 2x150 | Low biomass |
| 20 | Respiratory metatranscriptome | PRJNA671738 | Illumina | Nasal RNA virome |
| 21 | Fujimoto Japanese 4D (10 samples) | PRJNA862966 | NovaSeq 2x150 | Batch testing subset |
| 22 | Kleiner bacterial mock | ERR1877475 | NextSeq 2x76 | Non-virome control |
| 23 | Tweedy ciHHV-6 carriers | PRJNA412600 | Illumina | Known HHV-6 integration |
| 24 | LA wastewater metatranscriptome | PRJNA119800130 | Illumina | Deep RNA environmental |
| 25 | CRCbiome Norwegian (if accessible) | EGAS50000000170 | NovaSeq 2x151 | FIT sample extraction |

## Tier 3: HPC Atlas Generation (~5,000 samples)

Large-cohort processing to build empirical expected ranges from real data.

| Cohort | Accession | Samples | Purpose |
|--------|-----------|---------|---------|
| Fujimoto Japanese 4D | PRJNA862966 | 4,198 | Population baseline, batch effects, ERV rates |
| Tisza wastewater longitudinal | PRJNA966185 | 363 | Environmental temporal stability |
| Negativeome study (all blanks) | PRJEB33578 + others | 55 | Blank characterization |
| HIV+ mothers/infants | PRJNA1045584 | 306 | Clinical + retroviral baseline |
| Shkoporov longitudinal | PRJNA545408 | ~100 | Temporal QC metric stability |

**Total: ~5,000 samples → ~5,000 passport.json files**

This becomes the virome-qc-atlas: empirical expected ranges from thousands of real samples. Novel data includes population-level ERV rates (never reported across a virome cohort at scale) and batch effect quantification.

### HPC Requirements

- Storage: ~50 TB for raw FASTQ (download + delete after processing, stream if possible)
- Compute: virome-qc at ~34K reads/sec per thread, 6 GB RAM per process (SILVA + T2T filters)
- Output: ~5,000 passports at ~100 KB each = ~500 MB
- Estimated wall time: ~4 days at 8 parallel processes
- Bottleneck is SRA download bandwidth, not processing

## Changes from Previous Set

| Dataset | Action | Reason |
|---------|--------|--------|
| Buddle WGS/RNA | UPGRADE to full run | -X 5M truncation hides platform artifacts |
| Chrisman DNBSEQ | UPGRADE to full run | -X 10M truncation limits analysis |
| Zhang undepleted/depleted | KEEP at 1M spots | Intentional for manageable size |
| Kleiner | MOVE to Tier 2 | Bacterial mock, not a virome |
| Cook ONT/PacBio | ADD | Multi-platform comparison from same mock |
| Tisza wastewater | ADD | Best-characterized environmental virome |
| Negativeome blanks | ADD | Essential contamination baseline |
| HIV+ virome | ADD | Active retrovirus for ERV validation |
| 72-virus mock | ADD | Most complex synthetic community |

## Key Citations

- Cook et al. 2024. Microbial Genomics. PRJEB56639. Multi-platform phage mock.
- Santos-Medellin et al. 2021. ISME Journal. PRJNA646773. Soil VLP viromes.
- Shkoporov et al. 2019. Cell Host & Microbe. PRJNA545408. Longitudinal gut virome.
- Tisza et al. 2023. Nature Communications. PRJNA966185. Wastewater virome.
- Fujimoto et al. 2022. Nature Communications. PRJNA862966. Japanese 4D gut virome.
- Negativeome 2025. Nature Communications. PRJEB33578. Extraction controls.
- HIV+ virome 2024. PMC 11065063. PRJNA1045584. HIV+ mothers/infants.
- Tweedy et al. 2018. Scientific Reports. PRJNA412600. ciHHV-6 carriers.
- 72-virus synthetic community 2024. Journal of Virology.
- Zolfo et al. 2024. Microbial Genomics. PRJEB47625. Paired ONT/Illumina gut virome.
