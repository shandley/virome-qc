# Public datasets for EVE false-positive validation

Catalog of publicly available sequencing datasets and reference resources for future work
on endogenous-viral-element (EVE) false positives in virome-qc: quantifying how often host
EVEs (HERVs, endogenous bornavirus/parvovirus elements, ciHHV-6) get called as real
viruses, and validating an EVE-screening step across metagenomic DNA and RNAseq.

Provenance: assembled from a deep-research pass (six search angles, 26 sources fetched,
25 claims adversarially verified 25/25 confirmed) plus accessions verified directly
against NCBI/ENA/UCSC/Dfam earlier this session. Access status is the load-bearing field;
"controlled" means dbGaP/GDC application required. Before downloading any run, do a final
live accession + size check against SRA/ENA (byte-level status was not confirmed for every
entry here).

## Category 1: human gut virome metagenomes (the flagship input, fully open)

- Unified Human Gut Virome (UHGV) — integrates the four catalogs the field uses (MGV, GPD,
  CHVD, GVD) into 873,995 genomes / 168,536 species-level vOTUs. OPEN: JGI portal
  (uhgv.jgi.doe.gov/downloads) + Zenodo DOI 10.5281/zenodo.17402089 (v1.0, 2025, CC-BY);
  github.com/snayfach/UHGV. Best single starting point for the phage/EVE reference side.
- Gut Virome Database (GVD) — 2,697 VLP-enriched or bulk stool samples, 1,986 individuals,
  16 countries, collapsed to 33,242 viral populations. OPEN: iVirus/CyVerse DOI
  10.25739/12sq-k039; raw reads via SRA/iVirus/MG-RAST (Table S1). Gregory et al., Cell
  Host & Microbe 2020.
- Early-Life Gut Virome (ELGV) — 160,478 non-redundant DNA and RNA viral sequences from
  8,130 VLP-enriched or bulk gut metagenomes (first 3 years of life). OPEN. Nat Commun
  2024 (s41467-024-45793-z).
- Shkoporov et al. 2019 (PMID 31600503, Cell Host & Microbe) — longitudinal VLP-enriched
  adult gut virome; the repo already uses Shkoporov data. OPEN via SRA (confirm accession).
- Constituent catalogs (MGV, GPD, CHVD) exist standalone but their individual accessions
  were not independently re-verified here; UHGV supersedes them for most uses.

Recommended start: UHGV (reference/vOTUs) + GVD or Shkoporov raw reads (real VLP viromes to
screen against the EVE map).

## Category 2: human tissue / clinical mNGS (host-rich; partial, has gaps)

This category was the least covered by verified claims and is a genuine gap to fill. Leads
surfaced but not fully verified:
- mSystems 2023 (10.1128/msystems.00907-22) — "Elimination of foreign sequences in
  eukaryotic viral reference genomes"; documents HERV-H contaminating viral reference
  genomes. Directly on the EVE-in-reference problem; check for associated data.
- Clinical mNGS with public raw data: PMC10721142 and Nat Med 2024 (s41591-024-03275-1)
  were surfaced; access and EVE relevance need direct verification.
- TCGA tumor RNA/DNA is host-rich and EVE-abundant but CONTROLLED (GDC).

Open question to resolve: which fully-open tissue metagenomes / plasma-CSF mNGS with raw
reads exist for EVE testing. Worth a dedicated SRA/ENA query pass.

## Category 3: RNAseq viral detection (where HERV FPs are worst)

The problem is well quantified in the literature:
- HERV transcripts are pervasive in normal tissue: 13,889 locus-specific expressed
  hervRNAs across 42 GTEx body sites (Genome Biology 2022, s13059-022-02804-w); HML-2
  expressed >=1 TPM in every body site, 37 proviruses >1 TPM (PLOS Biology 2022,
  pbio.3001826).
- Mechanism + mitigation: <50% of TE-aligning fragments are uniquely assignable; naive
  best-count assignment misattributes ~12.1% of fragments to unexpressed HERV loci, vs
  <0.1% with Telescope's Bayesian EM (Bendall et al., PLoS Comput Biol 2018,
  pcbi.1006453). Telescope is the standard locus-resolution method to compare/beat.
- Documented real FPs in host-rich RNAseq: TCGA had HeLa-HPV18 in 131 non-cervical samples
  and XMRV in 96, traced to a QC reference RNA pool (BMC Genomics 2020, s12864-020-6483-6).
- GTEx raw reads are CONTROLLED (dbGaP phs000424, now v10.p1). Open routes around it:
  recount2/recount3 (uniformly processed, includes unmapped reads) and individual
  open-access human tissue RNAseq studies on SRA/ENA. GTEx summary expression matrices are
  open but not useful for read-level viral detection.

Recommended start: an open SRA tissue-RNAseq study (or recount3 unmapped reads) run through
viral detection, with Telescope as the HERV-resolution baseline; expect the highest EVE FP
fraction here.

## Category 4: ground-truth and benchmark datasets

- Cross-biome virus-ID benchmark (Genome Biology 2024, s13059-024-03236-4). OPEN:
  BioProjects PRJEB71789 (Antarctic seawater), PRJNA646773 (tomato soil), PRJNA389927
  (human gut) + Zenodo DOI 10.5281/zenodo.10886947 (assembled contigs + simulated
  positive/negative fragments, ~61,550 each). Ground truth by 0.22um size fractionation +
  DNase (contigs >=1500 bp). Host/prophage/EVE content is NOT removed from the microbial
  fraction, so a microbial-fraction contig called "viral" is a ready-made FP proxy.
- CAMI2 (Nat Methods 2022, s41592-022-01431-4). OPEN: computationally-defined truth from
  ~1,700 genomes + ~600 plasmids/viruses; viral detection underperformed bacterial.
- ViroForge (already used by the repo): synthetic reads with PLANTED EVEs + exogenous
  viruses and exact labels. The only source with true EVE ground truth. Best for measuring
  the EVE screen's TP/FP directly. No access barrier.
- Threshold method (transferable): a synthetic-virome study defines a per-molecule FP
  threshold at mean + 3 SD of background (JVI 2023, jvi.01300-23). Candidate calibration
  method, though derived for multiplex crosstalk not EVEs.

Recommended start: ViroForge (planted-EVE ground truth) for the controlled screen
validation; the cross-biome Zenodo set for an open real-data FP proxy.

## Category 5: endogenous-vs-exogenous discrimination (the FN side)

- Parrish et al. 2021 (PLOS Genetics, PMC8101998) — virus-derived variation in 3,332
  high-coverage WGS from 1000 Genomes + HGDP (LCLs). Documents heritable ciHHV-6, polymorphic
  HERV-K, and non-germline HIV-1/HTLV-1 integrations. The key OPEN EVE-integration reference
  catalog for this side.
- Verified ciHHV-6A carriers with OPEN unbiased WGS: NA18999 (Japanese, 35x, full-length
  endogenous HHV-6A) and HG00657 (Chinese CHS, chr22q), both in 1000 Genomes. Confirmed
  this session via PLOS Genetics pgen.1008915. Use the unbiased 30x WGS (NOT the enriched
  SRR6118055, which broke copy-number) for the copy-number signal.
- Active-infection counterparts for the exogenous side (HHV-6/HTLV/MLV/XMRV): still to be
  sourced as OPEN metagenomic/mNGS runs (open question). Note: HHV-6B has 218 SRA records
  under organism, but most are targeted, not metagenomic.

Recommended start: NA18999 or HG00657 30x 1000G WGS (open, unbiased, host-rich) for the
ciHHV-6 copy-number + host-confirmed-junction test.

## Reference resources (EVE annotation + viral panels)

- Dfam (dfam.org) — OPEN TE/HERV families + genome annotations; coordinate annotation API
  (used this session to name the chr1 HML-2 element and pull LTR5/LTR7/LTR17 accessions).
- RepeatMasker for T2T-CHM13 (hs1) — OPEN: hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/
  hs1.repeatMasker.out.gz (173 MB); used to build the genome-wide HERV set B.
- gEVE (genome-based endogenous viral element database, DB baw087) — OPEN: mammalian EVE
  ORFs (nucleotide + AA + loci) across 20 genomes. The route to source nrEVE (bornavirus/
  parvovirus) loci that RepeatMasker does not annotate.
- HERVd / RepBase — HERV catalogs; RepBase is license-restricted (commercial), Dfam is the
  open substitute.
- Exogenous viral panels: RVDB (github.com/ArifaKhanLab/RVDB, curated, handles endogenous),
  Virosaurus (clustered for clinical metagenomics), NCBI RefSeq viral. OPEN.

## What we already built this session (slots into the above)

- Set B (host EVE map): references/set_b_herv_rmsk.bed/.fa — 3,244 HERV loci (HML-2 +
  HERV-H + HERV-W), 10.2 Mb, from the T2T RepeatMasker source above. The EVE-screen reference.
- Set D anchors (exogenous): references/set_d_anchors.fa — HHV-6A NC_001664.4, HHV-6B
  NC_000898.1, MLV Moloney/Friend/typeC, MMTV, verified.
- Real K113: references/K113_NC_022518.1.fa.
- Scripts: build_set_b_rmsk.py (EVE map builder), cihhv6_real.py (discrimination engine on
  real reads), containment_test.py.

## Recommended overall starting set (2-3 fully open per goal)

- EVE-screen validation (controlled truth): ViroForge planted-EVE data.
- Gut virome real FP quantification: UHGV reference + GVD/Shkoporov raw reads.
- RNAseq HERV FP (highest signal): an open SRA tissue-RNAseq study or recount3, Telescope
  as baseline.
- Open real-data FP proxy: cross-biome benchmark (Zenodo 10.5281/zenodo.10886947).
- ciHHV-6 discrimination (FN side): NA18999 / HG00657 30x 1000G WGS.

## Gaps to close

1. Category 2: find fully-open tissue/clinical mNGS with raw reads (dedicated SRA/ENA pass).
2. Open active-infection runs for HHV-6/HTLV/MLV/XMRV (exogenous side of the pairs).
3. Confirm standalone MGV/GPD/CHVD accessions if UHGV integration is insufficient.
4. Live byte-level accession + size checks before any download.
