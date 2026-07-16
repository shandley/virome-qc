# EVE redesign: findings and loose-ends tracker

Running log of the small but load-bearing things surfaced while building the EVE
screen / discrimination redesign, so they do not get lost. Cross-refs:
EVE_DISCRIMINATION_REDESIGN.md, EVE_SCREEN_MODULE.md, REFERENCE_CURATION.md,
EVE_VALIDATION_DATASETS.md.

## Set-D exogenous panel (accumulating)

Verified against NCBI this session, in references/set_d_anchors.fa:
- HHV-6A NC_001664.4 (159,378 bp), HHV-6B NC_000898.1 (162,114 bp) -- ciHHV-6 counterparts.
- MLV Moloney NC_001501.1, MLV Friend FB29 NC_001362.1, Murine type C NC_001702.1, MMTV
  NC_001503.1 -- mouse retroviruses.

Verified this session, used in the EVE-screen ground truth (HTCF hiv_htlv.fa):
- HIV-1 NC_001802.1 (9,181 bp), HTLV-1 NC_000858.1 (9,028 bp).

nrEVE-family exogenous counterparts, verified + fetched (references/set_d_nreve_exo.fa,
66,052 bp): Bornaviridae -- NC_030692.1 (Mammalian 1 orthobornavirus / Borna disease virus;
strain sequenced is BoDV-2, and there is NO standalone BoDV-1 RefSeq, so use this species
reference or add a specific BoDV-1 GenBank clone), VSBV-1 NC_030701.1; Parvoviridae --
parvovirus B19 NC_000883.2, AAV-2 NC_001401.2; Filoviridae -- Ebola Zaire NC_002549.1,
Marburg Musoke NC_001608.4 (no human EVE, detection-only). These supersede the mislabeled
sequences in src/modules/eve_sequences_nonretro.rs, which should be relocated to set D.

Still to source: XMRV (no RefSeq; GenBank VP62), HIV-1 group-M subtype spread + CRFs, HIV-2,
HTLV-2, HHV-7, human foamy/spumavirus, human bocavirus.

Protein DIAMOND DBs built 2026-07-16 (references/, provenance in protein_db_PROVENANCE.md):
- set_d_prot.dmnd: 262 NCBI-annotated CDS proteins from 14 accessions (adds HIV-1
  NC_001802.1 + HTLV-1/STLV-1 NC_000858.1 to the panel above); gag/pol/env + RdRP verified.
- eve_prot.dmnd: 657,297 peptides from 6-frame translation of the 21,302 internal HERV loci
  (ERV1/ERVK/ERVL) of the full-HERV reference; scripts/translate_orfs.py, >= 50 aa.
- Competitive test passed: HERV-K locus delta -402 (suppress); MLV Moloney delta +2885
  (keep) despite MLV hitting EVE-prot at 47.5% id -- the case that proves the competitive
  comparison is necessary. nrEVE ORFs to be appended + rebuilt when detectEVE finishes.

minimap2 nucleotide indices built 2026-07-16 on HTCF (/scratch/sahlab/shandley/eve_exclusion),
asm20 + sr presets: eve_nt.{asm20,sr}.mmi (full set_b_hervfull.fa) and set_d_nt.{asm20,sr}.mmi
(14-genome panel; set_d_nt.fa = anchors + nreve_exo + HIV-1 NC_001802.1 + HTLV-1 NC_000858.1).
Nucleotide test binary: HERV-K frag -> EVE-nt only; MLV frag -> set-D-nt only (no cross-align,
unlike the protein arm). All four screen references now exist; only nrEVE addition pending.

## Code cleanups / mislabelings to fix

- src/modules/eve_kmers.rs: the one embedded sequence is mislabeled
  "HERV-K_K113_AF164610.1". It is actually a chr1 HERVK-int (HML-2) locus (Dfam
  DF000000188); AF164610 is HERV-K102, and real K113 is NC_022518.1 (absent from CHM13).
  Regenerate with the fixed build_eve_exclusion.py or the RM-based builder.
- src/modules/eve_sequences_nonretro.rs: named "non-retroviral EVE" but holds exogenous
  references (set-D). Relabel and move to the set-D panel; it is NOT a host-EVE reference.
- Shipped 3-signal ERV classifier: DECOMMISSIONED 2026-07-16. It was wired live (ErvScreener
  every run -> run_erv_analysis -> report card) but its signals were inert (MinHash const 0.5,
  family always Retroviridae, polymorphism k-mers empty), so it emitted CpG-only endo/exo
  labels. Removed: run_erv_analysis (main.rs), erv_classifier.rs, erv_pipeline.rs,
  eve_polymorphisms.rs, passport.erv_analysis field, orphaned cpg_ratio, and the report-ui
  ErvAnalysisCard (bundle rebuilt). Retroviral read DETECTION (ErvScreener) retained. Build +
  clippy clean, 123 tests pass. Doc claims stripped from README/PUBLICATION_PLAN/VALIDATION/TODO.
- Cargo.toml: biometal and superbloom paths FIXED (2026-07-16) to
  /Users/shandley/Code/software/{biometal,SuperBloom}. Was /Users/scotthandley/...

## Incidental verified accessions / facts

- Real HERV-K113 = NC_022518.1 (RefSeq) = AY037928.1 (Turner 2001), 9,472 bp; polymorphic,
  absent from CHM13.
- The mislabeled chr1 element = HERVK-int (HML-2), Dfam DF000000188.
- NC_060925.1 = CHM13 chr1; full chr<->NC map from NCBI assembly report GCF_009914755.1.
- T2T RepeatMasker source: UCSC hs1.repeatMasker.out.gz (173 MB), cached on HTCF.
- Cohort accessions: RISK/GSE57945 = SRA study PRJNA248469 (CD 218 / UC 62 / Control 42);
  Shkoporov gut VLP = PRJNA545408 (run SRR9161520); HMP2 IBD stool MGX+MTX = PRJNA398089
  (verified); RISK biopsy host RNA-seq = PRJNA438663.
- gEVE database (geve.med.u-tokai.ac.jp) is DEFUNCT (domain does not resolve).
- After the Cargo.toml path fix, `cargo check` builds the whole tool cleanly on this
  machine (biometal + SuperBloom + virome-qc, incl. the uncommitted eve modules), 11s, only
  dead-code warnings. So the Rust side is unblocked (run the tool, wire in the screen).

## Design decisions / notes

- Host depletion (T2T k=31 containment) already removes fixed HERV DNA reads; even the
  polymorphic K113 gets ~92% containment via paralog pooling. So k-mer subtraction decoys
  are redundant with T2T; the EVE map's value is a classification-stage screen, not a decoy.
- EVE false-positive severity is set by REGIME (host content), not body site: healthy gut
  VLP 0.003%, active-IBD ileal biopsy RNA-Seq 0.559% (~186x).
- Copy-number-vs-host is valid only for unbiased (metagenomic/WGS) data, not targeted/
  enriched; the failed ciHHV-6 real-data test used an enriched library.
- Full HERV reference: ERVK + ERV1 + ERVL, MaLR excluded (degraded, not FP-prone).
- The production EVE screen is ALIGNMENT-based (minimap2 nt + DIAMOND blastx protein),
  competitive EVE-vs-set-D, run on viral candidates. A k-mer / SuperBloom containment screen
  was evaluated and REJECTED (2026-07-16): redundant with T2T host depletion (EVE k-mers are
  a subset of T2T) AND blind to the reads that matter -- the FPs that reach a classifier
  survived k-mer host depletion by being divergent, and a k-mer screen reuses the same
  k-mers. K-mer containment stays only as a fast lower-bound measure of raw host-EVE content
  (what the regime study used). Full spec in EVE_SCREEN_MODULE.md.
- ViroForge does not plant EVEs (ground truth is viruses + contamination), so EVE ground
  truth is constructed by spiking.
- No human endogenous filovirus exists; the Filoviridae EVE family is a non-entity for
  human samples.
- detectEVE on T2T-CHM13 COMPLETE (2026-07-16; output in benchmark_data/results/detecteve/
  validatEVEs.tsv + .fna, 144 loci). Curation: 74 loci are HERV hits, all LOW-confidence,
  retroviral -> redundant with set B, drop. 42 high-confidence, all non-retroviral. Of those
  the CREDIBLE human nrEVEs are ~21 bornavirus EBLN loci (Carbovirus + Orthobornavirus) and
  1 parvovirus EPV (AAV-2). The tail (plant Tymovirus/Marafivirus, fungal Fusarium flexivirus,
  fish carp cultervirus, influenza, enterovirus, animal herpesviruses) is cross-kingdom
  repeat/low-complexity artifact -- NOT genuine human nrEVE, exclude. So the fold-in set is
  the EBLN + EPV subset, not all 42. Needs a curation pass before adding to set B/C + rebuild
  of eve_prot/eve_nt.