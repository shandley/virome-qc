# EVE screen: motivating evidence and module design

Host endogenous viral elements (EVEs), chiefly human endogenous retroviruses (HERVs), are
routinely mis-called as real viruses in virome and metagenomic analyses. This document
records the experiment that quantifies the problem across sample regimes (the motivating
evidence, intended to be publication-grade), then scopes an EVE screen as a real virome-qc
module.

Companion docs: EVE_DISCRIMINATION_REDESIGN.md (the broader subtract-then-discriminate
design and the endogenous-vs-exogenous / false-negative side), REFERENCE_CURATION.md (the
reference sets), EVE_VALIDATION_DATASETS.md (public datasets).

## Summary

The severity of EVE false positives is set by the host content of the sample, not by the
body site. A read-level EVE screen (k-mer containment against a T2T-derived HERV map) was
run across three regimes. On constructed ground truth it removed 100 percent of host-EVE
reads with zero real viruses lost. On a healthy gut virus-like-particle (VLP) metagenome
the EVE load was minor (0.003 percent of reads). On an active-IBD ileal-mucosa RNA-Seq
biopsy it was 0.559 percent, about 186 times higher, which extrapolates to roughly 106,000
host-HERV reads in a single sample, each of which a naive classifier assigns to a
retrovirus. The conclusion is that EVE false positives are a major, specificity-critical
problem in host-rich, disease-state, and RNAseq samples (the clinically important ones),
and are minor only in healthy VLP-enriched viromes.

## Motivating evidence (the study)

### Rationale

HERVs make up roughly 8 percent of the human genome and are actively transcribed. In
virome and metagenomic pipelines their reads are assigned to retroviruses, producing false
positives. Published HERV-detection false-discovery rates in RNAseq range from 8 to 55
percent (see EVE_VALIDATION_DATASETS.md). The question this study answers: how does the EVE
false-positive load vary across the sample types virome-qc actually processes, and does a
simple EVE screen catch it.

### Methods

EVE reference map (set B). All HERV-K (HML-2), HERV-H, and HERV-W loci in the T2T-CHM13
(hs1) assembly, extracted by their RepeatMasker coordinates
(hs1.repeatMasker.out.gz, UCSC) with scripts/build_set_b_rmsk.py: 3,244 loci, 10,200,594
bp, 0 percent N. Hashed to 4,772,615 canonical k=31 k-mers. This is 0.33 percent of the
3.1 Gb genome.

EVE screen. For each read, compute k-mer containment (fraction of the read's canonical
31-mers present in the set-B k-mer set); flag as EVE-derived if containment exceeds a
threshold (0.15; "strong" match at 0.50). scripts/eve_screen.py. A read flagged EVE would
be mis-assigned to a retrovirus by a naive classifier, because HERV sequence is retroviral.

Datasets and accessions (all verified against NCBI/ENA this session):

- Ground truth (constructed). ViroForge does not plant EVEs (its ground truth is viruses
  plus contamination), so EVE ground truth was constructed directly: host-EVE reads
  simulated from the set-B map, and real exogenous-virus reads simulated from HIV-1
  (NC_001802.1, 9,181 bp) and HTLV-1 (Primate T-lymphotropic virus 1, NC_000858.1, 9,028
  bp). Both classes are retroviral, so a naive retroviral classifier flags all of them;
  every host-EVE read is therefore a false positive. 5,000 reads per class, 150 bp.
- Healthy gut VLP. Shkoporov et al. 2019 (Cell Host and Microbe; BioProject PRJNA545408),
  run SRR9161520, VLP-enriched adult stool. First 1,000,000 reads.
- Active-IBD host-rich RNAseq. GSE57945 (the RISK cohort), treatment-naive pediatric IBD
  ileal-mucosa biopsy RNA-Seq, run SRR1783004, 51 bp single-end. First 1,000,000 reads.

Environment: HTCF, Python 3.11 (qa-genomics env). Reads streamed from ENA. Seed fixed (42).

### Results

| Sample | Regime | EVE match (>15%) | Strong (>50%) | Notes |
|---|---|---|---|---|
| Ground truth (host-EVE) | constructed | 5,000/5,000 = 100% | 100% | mean containment 100% |
| Ground truth (HIV/HTLV) | constructed | 0/5,000 = 0% | 0% | mean containment 0% |
| Shkoporov SRR9161520 | healthy gut VLP (DNA) | 29/1M = 0.003% | 3/1M | phage-dominated, host-depleted |
| IBD ileal mucosa SRR1783004 | active IBD biopsy (RNA) | 5,591/1M = 0.559% | 3,125/1M = 0.312% | mean containment 62% |

Derived figures:

- EVE screen on ground truth: recall 100 percent (all host-EVE flagged), specificity 100
  percent (no HIV/HTLV flagged). Naive retroviral false-positive rate 50 percent (the
  EVE reads among the calls); after the screen, 0.00 percent, with 0 real viruses lost.
- Host-rich vs healthy: 0.559 / 0.003 = ~186x more EVE-derived reads in the active-IBD
  biopsy than the healthy gut VLP.
- Absolute load: SRR1783004 is ~19 million reads; 0.559 percent is ~106,000 host-HERV
  reads, and the strong-match set (>50 percent containment, mean 62 percent) is ~59,000.

### Interpretation

- Regime, not body site, sets severity. The same gut compartment is benign when healthy
  and VLP-enriched (0.003 percent) and severe when inflamed and host-rich (0.559 percent).
  Active IBD floods the sample with human nucleic acid, moving a gut sample into the
  tissue-like, host-rich regime.
- The count is the metric, not the fraction. ~106,000 reads assigned to a retrovirus in
  one sample is an unambiguous, reportable false positive; for diagnostic metagenomics a
  far smaller count is a specificity failure.
- RNAseq amplifies the effect. The observed 0.559 percent exceeds the 0.33 percent
  genomic-fraction prediction because HERVs are transcribed; RNA-Seq over-represents
  expressed HERV loci beyond their genomic footprint. This is why the worst published
  HERV false-discovery figures come from transcriptomic data.
- The screen works. A simple k-mer containment screen separates host EVEs from real
  viruses cleanly, because host HERV reads saturate the EVE map (~100 percent containment)
  while HIV/HTLV, phylogenetically distant from all HERVs, register 0 percent.

### Limitations (to address before publication)

- The screen tested here is a read-level k-mer proxy, not the final competitive-identity
  module (below). It flags EVE-derived reads well but does not by itself distinguish a
  host EVE from a genuine exogenous retrovirus that is near-identical to one (the rare
  hard case; handled by the module's set C/D competitive step and coverage/junction).
- Ground truth was constructed by spiking, not generated by ViroForge, because ViroForge
  does not plant EVEs. A published version should add a planted-EVE simulator or use gEVE
  loci with known coordinates.
- One sample per regime is illustrative. A publication needs N samples per regime with
  distributions and statistics, ideally paired healthy-vs-active within cohorts (IBD).
- Set B covers HML-2, HERV-H, HERV-W only; total HERV is larger, and nrEVEs are excluded,
  so the true EVE false-positive load is an underestimate of what a full reference catches.
- Threshold-dependent; report a sweep. The transcript-vs-genomic mapping for RNAseq should
  be handled explicitly.

## Module scope: the EVE screen as a real virome-qc module

### Role and relationship to the redesign

The EVE screen is the false-positive side of the EVE_DISCRIMINATION_REDESIGN.md program:
at the point where reads or contigs would be called viral, flag those that are host-EVE
derived and suppress them from viral calls. The endogenous-vs-exogenous discrimination
engine (the false-negative side, for near-identical pairs like ciHHV-6 or mouse MLV) is a
separate stage; both draw on the same reference sets (B host EVE map, C endogenous
competitive, D exogenous panel).

### References

- Set B (built): HERV HML-2 + HERV-H + HERV-W loci from T2T. Extend to other HERV families
  (HERV-E, HERV-FRD, etc.) and to non-retroviral EVEs (EBLN, EPV) sourced from gEVE, since
  those are transcribed and matter for RNAseq.
- Set C / set D (partial): for the competitive form that also spares real retroviruses.

### Algorithm: alignment-based competitive screen

A k-mer / SuperBloom containment screen was evaluated and REJECTED for the production
screen (2026-07-16): (1) it is redundant with host depletion, because the EVE reference is
extracted from T2T, so its k-mers are already in the T2T host filter; (2) k-mer containment
needs ~95% identity to flag, so it has a blind spot on exactly the reads that matter -- the
EVE false positives that reach a classifier are the ones that SURVIVED k-mer host depletion
by being divergent, then got called viral by a homology (alignment/protein) classifier. A
k-mer screen uses the same k-mers that already failed to catch them. The k-mer screen
remains valid only as a fast, conservative lower-bound measure of raw host-EVE content (what
the regime study used), not as the authoritative FP screen.

The screen runs on viral CANDIDATES (contigs preferred; reads for RNAseq/unassembled), which
are few (post-classification), so alignment is affordable and is the sound choice.

Per candidate C, competitive alignment against EVE (host) and exogenous (set D) references:

1. Nucleotide: minimap2 C vs EVE-nt (full-HERV set_b_hervfull.fa + nrEVE) with -x asm20 for
   contigs (-x sr for reads); best hit -> (bitscore, %id, query-coverage, locus, family).
   minimap2 C vs set-D-nt -> best exogenous hit.
2. Protein (essential for RNAseq / deep divergence, and to match the protein-homology
   classifiers that drive the 8-55% RNAseq HERV FDR): diamond blastx C vs EVE-prot and vs
   set-D-prot; best bitscore each.
3. Competitive margin: delta = score(exogenous) - score(EVE), using bitscore (nt+protein
   combined; protein dominates for divergent hits).
4. Call:
   - strong EVE hit, no/weak exogenous hit, delta <= -d  -> EVE-DERIVED (suppress as FP)
   - strong exogenous hit, no/weak EVE hit, delta >= +d  -> GENUINE virus (keep)
   - both strong, |delta| < d (the near-identical pair)  -> AMBIGUOUS -> escalate to
     Stage 2 (coverage-breadth obs-vs-expected, copy-number-vs-host, host-confirmed
     junction) = the discrimination engine (ciHHV-6/HHV-6, mouse endo/exo MLV)
   - neither hit -> not an EVE; pass through
   Note: most HERVs have no exogenous counterpart, so the common case is "hits EVE only" ->
   EVE FP, unambiguous; HIV/HTLV are "hits set D only" -> genuine.

Thresholds (calibrate on ground truth + threshold sweep): min alignment length and an
identity floor to count a hit; diamond e-value <= 1e-5; competitive margin d starting near
delta-bitscore = 10 (the tale-of-caution EVE-like cutoff was <= -10 at contig level). Report
a sweep over d.

### Tooling and reference prep (per the redesign: established aligners, not k-mers)

- minimap2 (nucleotide) + DIAMOND blastx (protein) -- both established, matching the
  redesign's "alignment for load-bearing measurements" decision. Shell-out is fine; the
  screen is a post-classification analysis step, not the streaming hot path.
- References. All four screen references BUILT 2026-07-16 (provenance in
  references/protein_db_PROVENANCE.md; both discrimination tests passed):
  - protein/blastx: eve_prot.dmnd (657,297 peptides, internal HERV loci 6-frame), and
    set_d_prot.dmnd (262 NCBI-annotated CDS proteins, 14 accessions).
  - nucleotide/minimap2 (on HTCF /scratch/sahlab/shandley/eve_exclusion): eve_nt.{asm20,sr}.mmi
    (full set_b_hervfull.fa) and set_d_nt.{asm20,sr}.mmi (14-genome panel).
  Remaining: append translated nrEVE ORFs to eve_prot + nrEVE loci to eve_nt and rebuild
  those two in one pass when the detectEVE run finishes; then wire the screen into the tool.

### Outputs / passport

Per candidate: eve_hit (locus, family, %id, bitscore, qcov), exo_hit, competitive margin,
call (EVE_FP / genuine / ambiguous), reason. Passport reports the EVE-derived FP COUNT
(not just fraction) plus a per-family breakdown (HML-2 / HERV-H / HERV-W / other HERV /
nrEVE) and the ambiguous count. Regime-gated: run always, load-bearing for host-rich and
RNAseq profiles.

### Validation hook

The publication validation (below) uses THIS alignment screen -- not the k-mer proxy -- to
quantify FP and screen performance (ground-truth recall/specificity + the d threshold
sweep). The k-mer regime study stands as the motivating lower-bound evidence.

### Integration into virome-qc

- Placement: at or after the ERV-screener / viral-classification stage. Every read or
  contig that matches a viral reference is additionally screened against set B; matches are
  annotated EVE-derived and removed from the viral call set.
- Passport reporting: emit the EVE-derived read/contig count and fraction, a per-family
  breakdown (HML-2 / HERV-H / HERV-W / nrEVE), and a flag when the load is high. Report the
  count, not only the fraction, since the count is what drives a false call.
- Regime gating via profiles: run always (it is cheap), but treat it as load-bearing for
  host-rich and RNAseq profiles (tissue-truseq, RNA/metatranscriptome) where the study
  shows the load is 100x higher. For healthy gut VLP it is a low-cost safety net.

### Roadmap

1. Wire the validated v1 read-level screen into the pipeline with passport reporting.
2. Extend set B to all HERV families and add nrEVE (EBLN/EPV) for RNAseq.
3. Build v2 competitive contig-level screen (spares real retroviruses).
4. Publication validation: N samples per regime across healthy/disease and DNA/RNA, with a
   threshold sweep and comparison to a naive classifier and to Telescope-style locus
   resolution.

### Reproducibility (artifacts)

- Set B: scripts/build_set_b_rmsk.py; references/set_b_herv_rmsk.bed and .fa.
- Screen: scripts/eve_screen.py.
- Real K113 and set-D anchors: references/K113_NC_022518.1.fa, references/set_d_anchors.fa.
- Containment finding (host depletion already removes fixed HERV DNA reads):
  scripts/containment_test.py.
- Study results: benchmark_data/results/eve_screen/ (TSV, PROVENANCE.md, run logs).

## Validation run (publication design)

Turns the single-sample-per-regime illustration into distributions with statistics.

### Claims to test

- H1: EVE false-positive load rises with host content (regime): VLP virome < bulk stool
  metagenome < tissue biopsy.
- H2: within a cohort, active IBD > control (disease raises host content).
- H3: RNA > DNA at matched host content (HERV transcription amplifies the load).
- H4: the EVE screen removes host-EVE reads with high recall and specificity, driving a
  naive classifier's false-positive rate to near zero with no real viruses lost.
- H5 (mechanism): EVE load correlates with measured host fraction; H1-H3 all reduce to
  this. Host fraction is measured per sample (reads mapping to T2T-CHM13).

### Primary outcome

EVE load per sample = reads flagged by the screen (k-mer containment above the operating
threshold vs the FULL EVE reference), reported both as fraction and as count per million
reads (the count is what generates a false call). Secondary: EVE-derived fraction among
naive retroviral calls.

### Cohorts (chosen to isolate one factor each)

- Cohort 1, RISK / GSE57945 (OPEN, verified). Ileal-mucosa biopsy RNA-Seq; control, CD, UC.
  Isolates disease state at fixed tissue + protocol + RNA. Target N=25 per group.
- Cohort 2, HMP2 / IBDMDB (OPEN; ibdmdb.org and SRA; verify BioProject at build). Stool
  metagenome AND metatranscriptome from the SAME subjects (non-IBD, CD, UC; longitudinal).
  Isolates DNA-vs-RNA (paired, same sample) and disease state, in stool. Target N=25 per
  (nucleic-acid x disease) cell, drawn from active and remission/control timepoints.
- Cohort 3, healthy gut VLP DNA (low-host anchor): Shkoporov PRJNA545408 (verified) or GVD
  studies. Target N=25.

The regime axis (H1) spans Shkoporov VLP (lowest host) -> HMP2 stool MGX (mid) -> RISK
biopsy RNA (highest); pooled across all samples, EVE load is regressed on measured host
fraction (H5).

### Ground truth (H4, screen performance)

Extend the constructed planted-EVE benchmark: plant HERV + nrEVE reads + real-virus reads
(HIV/HTLV plus a broader viral panel) + microbial/phage background at controlled ratios,
sweeping the endo/exo identity gap to probe the hard boundary. Metrics: recall,
specificity, precision, AUC. Threshold sweep (containment 0.05-0.75) -> ROC -> operating
point.

### Reference upgrade (status)

Full HERV: DONE. build_set_b_rmsk.py gained a class-based mode (--herv-classes ERVK,ERV1,
ERVL); genome-wide on T2T it gives 99,035 loci / 183.7 Mb (~6% of the genome), vs 10.2 Mb
for HML-2/H/W. See benchmark_data/results/eve_screen_validation/. Screen note: ~184M
k-mers, so use a SuperBloom filter or a cached k-mer set, not a Python set.

nrEVE: NOT done, source problem. gEVE is defunct (domain dead); the repo's
eve_sequences_nonretro.rs is actually exogenous references (set-D material), not host EVEs.
Route: detectEVE/EEfinder on T2T to find EBLN/EPV loci, then extract by coordinate. Small
(~a handful of loci vs 99,035 HERV), so the HERV reference already covers the dominant load.

### Baselines and comparisons

- Naive retroviral classifier (reads matching a retroviral RefSeq panel) -> the
  false-positive inflation the screen removes.
- Telescope on the RNA-Seq samples -> position the EVE screen as a QC-stage filter,
  show concordance and that it flags what a naive viral classifier would call.
- With vs without the screen -> false-positive reduction, per-family breakdown.

### Statistics

- Per-group distributions (median + IQR, violin/box).
- H1: trend test across ordered regimes (Jonckheere-Terpstra) and regression on host
  fraction. H2: paired Wilcoxon within subject where longitudinal, else Mann-Whitney,
  with Cliff's delta and bootstrap CIs. H3: paired Wilcoxon DNA vs RNA (same sample).
  H5: Spearman on pooled samples, mixed model with subject random effect (longitudinal).
- Ground truth: recall/specificity/AUC with CIs. BH correction across the hypothesis family.

### Figures

1. EVE load per regime (headline). 2. EVE load vs measured host fraction (mechanism).
3. Disease state within RISK (control/CD/UC). 4. DNA vs RNA within HMP2 (paired).
5. Ground-truth ROC + threshold sweep. 6. False-positive reduction (naive vs post-screen).

### Compute plan (HTCF)

Per sample: stream a fixed 2M-read subset from ENA; run the EVE screen (k-mer, fast) and
host fraction (minimap2 to T2T or k-mer containment vs T2T). ~150-200 samples as a SLURM
array, a few hours total. Build the full EVE k-mer set once and cache it. Fixed seeds;
per-sample outputs to benchmark_data/results/eve_screen_validation/. Telescope, the
naive-classifier baseline, and the ground-truth sweep run as separate steps.

### Scope tiers (if effort-bounded)

- Core (sufficient for the main claim, ~100 samples): regime axis + host-fraction
  correlation (H1, H5) + RISK disease-state (H2) + ground-truth/threshold (H4).
- Full: add HMP2 DNA-vs-RNA (H3), Telescope comparison, nrEVE reference, and a second gut
  cohort for generalizability.

### Verify at build

HMP2/IBDMDB BioProject accession; RISK/GSE57945 sample-to-diagnosis mapping (GEO metadata);
Shkoporov PRJNA545408 run list; gEVE download for nrEVE; a retroviral RefSeq panel for the
naive baseline. Live accession + size checks before download, per the project discipline.
