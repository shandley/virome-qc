# Reference-set curation

Detailed scope for Step 1 of EVE_DISCRIMINATION_REDESIGN.md. Most of the method's
scientific correctness lives in what goes into these sets, so this is specified as
retrieval plus verification, not as a list of accessions. No accession in this document
or in the built sets is trusted from memory; every one is pulled from an authoritative
source and recorded with its provenance.

## The four sets and their distinct roles

The endo/exo problem is that endogenous and exogenous copies share sequence. The
reference sets resolve this by role, not by trying to make the sequences disjoint (they
are not). Four sets, three roles:

| Set | Role | Used to |
|---|---|---|
| A. Host pangenome | subtract | remove host reads including integrated EVE content |
| B. EVE decoy | subtract | remove EVE reads the pangenome misses (rare/private, divergent) |
| C. EVE competitive/consensus | discriminate | measure "closer to an EVE than to a circulating virus" |
| D. Exogenous panel | discriminate and rescue | measure "closer to a circulating virus", flag real infection |

Sets A and B are the subtraction reference. Sets C and D are the two sides of the
competitive-identity comparison in Stage 2. A sequence can legitimately appear in both a
subtraction set and a competitive set; the roles differ.

## The inclusion rule that ties it together

This is the scientific crux of curation. A subtraction set must never contain sequence
that closely matches a circulating exogenous virus, or real infections get subtracted
(false negatives). So placement of each EVE depends on whether it has a near-identical
circulating relative:

- Distant EVE, no close circulating relative (ancient HERVs, nrEVEs): put in the
  subtraction set (B) and the competitive set (C). Safe to subtract because no infection
  resembles it.
- Recent EVE with a near-identical circulating form (ciHHV-6 vs reactivating HHV-6;
  endogenous MLV vs exogenous MLV/XMRV): put in the competitive set (C) only, never in
  the subtraction set. Subtraction cannot separate these from a real infection, so the
  separation is deferred to Stage 2 coverage, copy-number, and junction signals.

Decision criterion for "near-identical", to be calibrated: maximum identity between the
EVE consensus and the nearest exogenous reference over the shared region above a
threshold (start near greater than 90 percent over greater than a few hundred bp). Every
EVE is scored against set D and routed to B+C or C-only by this rule, and the routing is
recorded in the manifest.

Note that the host pangenome (set A) already subtracts many polymorphic HERVs, because
they are real host sequence present in some individuals. So set B mainly needs to cover
EVEs the pangenome does not contain: rare or private insertions, and divergent ancient
elements a per-individual k-mer filter may miss.

## Set A: host pangenome (subtraction)

- Sources: T2T-CHM13v2.0 for the complete human reference, plus HPRC pangenome
  assemblies for polymorphic and structural content a single reference misses. For mouse
  samples, the mouse reference (GRCm39) plus available mouse pangenome assemblies.
- Retrieval: HPRC assemblies from the human-pangenomics distribution; T2T and GRCm39
  from NCBI. Record assembly accession and version for each.
- Packaging decision (open, from the redesign doc): full pangenome into the SuperBloom
  k-mer filter, or T2T plus a curated decoy FASTA. Measure filter size, build time, and
  false-positive reduction both ways on ground truth before choosing.
- Host-matched: the human and mouse host sets are separate and selected by the sample's
  profile. Do not mix.

## Set B: EVE decoy (subtraction side)

- Purpose: genomic EVE sequence (actual integrated copies), so reads from EVE loci are
  removed. This is the same idea as the T2T viral-homology map in EVE_MAP_DESIGN.md.
- Sources and tools:
  - Extract EVE-annotated intervals from T2T/HPRC using RepeatMasker HERV annotations
    and Dfam HMM models (HERV-K/HML-2, HERV-H, HERV-W, HERV-E, and others).
  - Run detectEVE or EEfinder on the host assemblies for a general, tool-defined EVE set
    rather than hand curation.
  - Cross-check coverage against the RepeatMasker HERV call set, per EVE_MAP_DESIGN.md
    (the built set should cover all RepeatMasker HERV calls).
- Apply the inclusion rule: only EVEs that pass the distant-relative test enter set B.
- Mask host-acquired genes and low-complexity or highly conserved domains that cause
  chance cross-mapping.

Prototype status (HERV-K/HML-2): `scripts/build_set_b_hervk.py` implements the
coordinate-based route and is validated genome-wide on HTCF. It pulls every HERVK
internal locus (Dfam family DF000000188) from the Dfam annotations API, translates UCSC
chrom names to the T2T RefSeq accessions, and extracts the genomic sequences with pysam.
Genome-wide result: 114 HML-2 loci, 417,306 bp, 0% N, sizes 531 bp to 9,553 bp, spread
across all 24 chromosomes (chr1/8/11 richest). Artifacts: references/set_b_hervk_genome.bed
and .fa (plus an eve_kmers.rs-format Rust file on HTCF). This replaces the alignment-based
`build_eve_exclusion.py` for fixed host EVEs, which cannot work here: a single HML-2
provirus query captures nothing at high identity because paralogs are only ~78-81%
identical.

LTR flanks added: the HML-2 LTR families were verified against Dfam and added to
`build_set_b_hervk.py` (default families now internal + all four LTRs):

    DF000000188  HERVK    (internal)
    DF000000540  LTR5
    DF000000556  LTR5A
    DF000000557  LTR5B
    DF000000558  LTR5_Hs

The script queries each family per window (Dfam takes one family per call and has no
genome-wide family endpoint), merges internal+LTR hits within --merge-distance into full
proviruses, and labels each locus provirus / solo_LTR / internal_only. On chr19 (test):
98 loci = 7 proviruses + 2 internal-only + 89 solo LTRs, so solo LTRs dominate. Scope
chosen: comprehensive (include solo LTRs) for the most complete host-HERV-K exclusion.

Genome-wide comprehensive result (HTCF job 43583825, ~2h53m, all 5 families): 1,755 raw
Dfam hits -> 1,182 merged loci, 1,645,380 bp, 0% N, ~3.29M k=31 k-mers. Structure: 109
proviruses (internal+LTR; max 12,367 bp), 1,059 solo LTRs (~1.08 Mb; typical ~1 kb), 14
internal-only. Richest on chr8/chr1/chr19. The 109 proviruses match the ~90-100 known
full-length HML-2 proviruses, a good sanity check. Artifacts:
references/set_b_hml2_genome.bed and .fa (eve_kmers.rs-format Rust on HTCF).

Extended to HERV-H and HERV-W (broader exclusion): the per-family Dfam API does not scale
to abundant families (HERV-H alone has ~5,800 internal copies + thousands of LTR7; broad
coverage would need tens of thousands of calls). Switched to the bulk CHM13 RepeatMasker
annotation (UCSC hs1.repeatMasker.out.gz, 173 MB) with a new extractor,
`scripts/build_set_b_rmsk.py`, which reads the .out once and pulls every targeted family.
Families matched by an exact repName allowlist (prefix-matching "LTR7" would wrongly grab
LTR78/LTR79, distinct ERV1 families):
    HML-2  (LTR/ERVK): HERVK-int, LTR5, LTR5A, LTR5B, LTR5_Hs
    HERV-H (LTR/ERV1): HERVH-int, HERVH48-int, LTR7, LTR7B, LTR7C, LTR7Y
    HERV-W (LTR/ERV1): HERV17-int, HERV17B-int, LTR17, LTR17b
Comprehensive result (13,110 RM rows -> 3,244 loci, 10,200,594 bp, 0% N, ~20.4M k-mers):
HML-2 1,184 loci / 1.78 Mb; HERV-H 1,412 loci / 7.05 Mb (1,241 proviruses); HERV-W 648
loci / 1.38 Mb. Provirus mean sizes 4.4-5.8 kb (correct for internal+LTR). Artifacts:
references/set_b_herv_rmsk.bed and .fa. Cross-validation: RM HML-2 = 1,184 loci vs the
Dfam-API run's 1,182 loci (within 0.2%), so the two independent methods agree.

The RM extractor supersedes the Dfam-API one for production (one consistent source, scales
to any HERV family); the Dfam-API version stays as the validated prototype. Total set B is
~10.2 Mb, ~0.33% of the genome, so the host-depletion impact of excluding these k-mers is
modest. The progress prints in build_set_b_hervk.py now flush (live per-chrom progress on
long runs).

## Set C: EVE competitive/consensus (identity side)

- Purpose: the EVE side of the competitive comparison. Needs representative consensus
  sequences plus the specific loci that matter for the hard cases.
- Contents:
  - HERV consensus models from Dfam (per HML/HERV family).
  - The polymorphic HERV-K(HML-2) loci named in eve_polymorphisms.rs (K113, K106, K116).
    Their sequences and coordinates come from the polymorphic-HERV-K literature already
    cited in that file (Wildschutte 2016, Thomas 2018). Retrieve the actual sequences
    from the accessions those papers deposit; do not reuse the labels currently in the
    code without verification (see the K113/K102 finding below).
  - ciHHV-6 integrated reference genomes for HHV-6A and HHV-6B (routed C-only by the
    inclusion rule, since exogenous HHV-6 is near-identical).
  - nrEVE elements: endogenous bornavirus-like elements (EBLN-1, EBLN-2) from the
    bornavirus endogenization literature; endogenous parvovirus elements (EPVs).
- Human endogenous filovirus: none exists. Humans have no endogenous filovirus, so a
  human filovirus EVE set is empty and the Filoviridae EVE family is not applicable to
  human samples. It is relevant only for non-human hosts (bats, rodents). Do not carry an
  empty human filovirus family as if it were populated.

## Set D: exogenous panel (identity side and rescue)

Purpose: the circulating-virus side of the Stage 2 competitive-identity comparison, and
the set a real infection is rescued or flagged against. Set D is a small, targeted
discrimination panel, not a general viral-detection database. The boundary matters: a
circulating virus belongs in set D when a viral read that survived host+EVE subtraction
could plausibly be either that virus or a host EVE. Pure detection targets with no human
endogenous counterpart (Ebola, Marburg) are the job of the tool's general viral
detection, not of set D's discrimination role.

Inclusion rule. Put a circulating virus in set D if it is either (a) the exogenous
counterpart of an EVE family the tool subtracts or screens, or (b) a human retrovirus.
Reason for (b): the ERV screener flags any read with retroviral k-mer signal, so a real
HIV or HTLV read arrives at Stage 2 flagged and must be comparable against real human
retroviruses; the HERVs themselves have no circulating exogenous form, so the exogenous
side of the retroviral comparison is the actual human retroviruses.

Confusable-pair map (the analytical core). For each set D member, whether competitive
identity alone resolves the endo/exo call or whether it must fall back to the Stage 2
coverage / copy-number / junction signals:

| Set D member | EVE counterpart | Gap | Resolution |
|---|---|---|---|
| HHV-6A, HHV-6B | ciHHV-6 (set C) | near-identical | identity FAILS; needs copy-number + telomeric junction |
| exogenous MLV, XMRV (mouse) | endogenous MLV (mouse, set B/C) | near-identical | identity FAILS; needs copy-number + junction + diversity |
| HIV-1, HIV-2, HTLV-1/2 | none in human | far | identity resolves trivially; the "must not lose" anchors |
| BoDV-1, VSBV-1 | EBLN (set C) | ~100 My | identity resolves trivially; include for real Borna detection |
| parvovirus B19, bocavirus | EPV (set C) | far | identity resolves trivially |

The two near-identical pairs (HHV-6, mouse MLV) are the only cases where set C and set D
cannot be separated by sequence; they are exactly the cases the coverage/junction signals
exist for. Everything else set D resolves by identity alone.

Panel contents, human (selected when the sample host is human):
- Retroviruses: HIV-1 across the major group M subtypes (A, B, C, D, F, G, H, plus the
  common circulating recombinants CRF01_AE and CRF02_AG; groups N/O/P optional), HIV-2
  (groups A and B), HTLV-1, HTLV-2. Human foamy/spumavirus optional (rare zoonotic).
- Roseolovirus (the ciHHV-6 hard case): HHV-6A (canonical strain U1102), HHV-6B
  (strain Z29).
- Bornavirus (EBLN counterpart): BoDV-1, VSBV-1.
- Parvovirus (EPV counterpart): parvovirus B19, human bocavirus.

Panel contents, mouse (separate panel, selected when host is mouse; the acute
near-identical case and the XMRV mouse-DNA-contamination scenario):
- Gammaretrovirus: exogenous MLV (ecotropic, xenotropic, polytropic, amphotropic), XMRV.
- Betaretrovirus: MMTV.

Sourcing and curation:
- Prefer a curated low-redundancy resource built for metagenomic specificity: RVDB
  (explicitly handles endogenous sequence) or Virosaurus (clustered at 90/98%). Fall back
  to NCBI RefSeq viral, one complete genome per species/subtype, retrieved by taxon via
  NCBI datasets or E-utilities.
- Verify every accession at build time against the authoritative record; do not hardcode
  from memory (the same discipline that caught the K113/K102 error). Record source,
  accession, version, taxid, and retrieval query in the provenance manifest.
- HIV needs subtype diversity: a competitive-identity comparison is only as sensitive as
  its nearest reference, and HIV-1 subtypes diverge enough (well over 10% in several
  genes) that a single reference such as HXB2 would miss non-B reads. Include the subtype
  spread, not one genome.
- Dereplicate at ~95-97% identity so no single taxon dominates the panel by redundant
  near-duplicates.

Size target: small and targeted, roughly 20-35 human sequences (mostly the HIV subtype
spread plus the pair anchors) and ~6-8 mouse. Set D is a scalpel, not a database.

Constraint: no set D sequence may appear in the host/EVE subtraction sets (A or B), or
real infections get subtracted (a direct false-negative source). Verify explicitly (see
QA below). For HHV-6 this is automatically safe because ciHHV-6 is polymorphic and absent
from CHM13, so the subtraction set contains no HHV-6, but verify anyway.

Pairing with set C: for the near-identical pairs, set C holds the endogenous side and set
D the exogenous side; the endo/exo call for those comes from Stage 2 coverage-breadth,
copy-number-vs-host, and junction signals, not from set D identity alone.

Sourced anchors (the two hard near-identical pairs, verified against NCBI RefSeq this
session; fetched lengths matched the records exactly). Files:
references/set_d_anchors.fa and references/set_d_anchors_manifest.tsv (6 records,
355,087 bp).

    HHV-6A   NC_001664.4  159,378 bp  isolate U1102   (roseolovirus, ciHHV-6 counterpart)
    HHV-6B   NC_000898.1  162,114 bp  strain Z29      (roseolovirus, ciHHV-6 counterpart)
    MLV Moloney      NC_001501.1  8,332 bp   (exogenous gammaretrovirus, mouse)
    MLV Friend FB29  NC_001362.1  8,323 bp   (exogenous gammaretrovirus, mouse)
    Murine type C    NC_001702.1  8,135 bp   (exogenous gammaretrovirus, mouse)
    MMTV             NC_001503.1  8,805 bp   (exogenous betaretrovirus, mouse)

Notes: MMTV RefSeq is 8,805 bp (a memory of ~9.9 kb was wrong; the accession is correct
and length verified). XMRV has no RefSeq (it is a lab-derived recombinant of two mouse
endogenous ERVs), so it must be sourced from GenBank (VP62 clone) and verified separately
if the contamination scenario needs it. Disjointness holds for these anchors: HHV-6 is
not in the HERV-only set B and is absent from CHM13 (ciHHV-6 polymorphic), and the mouse
retroviruses are not in the human host/EVE sets. The remaining set D members (HIV subtype
spread, HTLV-1/2, BoDV-1/VSBV-1, B19/bocavirus) are sourced the same way when the full
panel is built.

## Cross-set QA checks

Run these on the built sets before use:

- Disjointness: no accession and no near-identical sequence appears in both a subtraction
  set (A or B) and the exogenous panel (D). Quantify k-mer overlap between C and D; it
  should be high for the near-identical pairs and near zero for distant EVEs, which
  confirms the inclusion routing.
- Completeness: set B plus the pangenome covers all RepeatMasker HERV calls on T2T.
- Leakage test: map a known exogenous positive against the subtraction set. HIV must not
  be subtracted (zero loss). Confirm ciHHV-6/HHV-6 is absent from the subtraction set so
  a real HHV-6 infection survives to Stage 2.
- Negative test: non-viral human sequence must not hit set D.

## Provenance manifest

One row per sequence, so the entire resource is reproducible and citable and no
identifier is ever trusted from memory. Columns:

- set (A/B/C/D) and role
- organism, taxonomy id
- source database, accession, version
- region (coordinates, if extracted from an assembly)
- retrieval date, retrieval query or command
- near_identical_pair flag and the identity score that set it
- notes

Store as TSV plus JSON. This manifest is itself a publishable resource and is what makes
a validation reproducible.

## Finding from verification: the existing K113 label is wrong (resolved)

The one EVE sequence currently embedded in eve_kmers.rs is tagged
`NC_060925.1_154765238_154774358_HERV-K_K113_AF164610.1`. Both the locus label and the
accession are wrong. Verified against NCBI:

- NC_060925.1 is CHM13 chromosome 1 (248,387,328 bp). HERV-K113 is a chromosome 19
  provirus, so the embedded interval cannot be K113. Named from the T2T annotation
  (Dfam): it is a fixed HERV-K (HML-2) provirus internal region, Dfam family HERVK
  (accession DF000000188, class LTR/ERVK), at chr1:~154766199-154773441 on the minus
  strand (near-full internal model, bit scores 6912 and 2110). It coincides with the
  T2T-annotated lncRNA AL353807.5 (chr1:154765322-154776167). Correct decoy label:
  `HERVK-int (HML-2, Dfam DF000000188), T2T-CHM13 chr1:154766199-154773441, minus`.
- AF164610.1 is "HERV-K102, complete sequence" (9,178 bp), not K113.
- The correct K113 reference is NC_022518.1 (RefSeq) = AY037928.1 (GenBank, Turner et
  al. 2001), "Human endogenous retrovirus K113 complete genome," 9,472 bp, deposited as
  a standalone proviral sequence with no chromosomal coordinate because K113 is a
  polymorphic unfixed provirus.

Root cause confirmed empirically (HTCF, minimap2 asm20 of real K113 NC_022518.1 vs
T2T-CHM13v2.0): K113 aligns full-length to ~12 HML-2 loci genome-wide but all at only
71.5-81.5% identity (best 81.5% at NC_060930.1/chr6). A present allele would be ~99%, so
CHM13 does not carry K113; every hit is a diverged paralog. The chr1 locus that was
mislabeled matches K113 at 78.4%, confirming it is a paralog. This is the concrete
demonstration of why polymorphic HERV-K loci cannot be extracted from a linear reference
(set C sourcing rule above): pull K113 from NC_022518.1, not from T2T. Real K113 is now
saved at references/K113_NC_022518.1.fa.

Secondary finding with a design consequence: at the default 85% identity threshold, a
single HML-2 provirus query captures no host loci at all, because any one HML-2 provirus
is only ~78-81% identical to the others. So the exclusion set (set B) cannot be built by
aligning one provirus at high identity; it needs either a much lower threshold or, better,
coordinate-based extraction of all host HML-2 loci from Dfam/RepeatMasker annotations (the
approach EVE_MAP_DESIGN.md already favors). This is why the original run, using a query
that happened to match its own fixed genomic copy, extracted exactly one locus.

## Action items and open questions

- K113/K102 mismatch resolved (see finding above). The chr1 element is named (HERVK-int,
  Dfam DF000000188). `build_eve_exclusion.py` is fixed: it now fetches curated references
  by verified accession from NCBI (K113 -> NC_022518.1, tested: 9,472 bp), names each
  extracted region by its T2T locus instead of fusing the query label onto it (the cause
  of the original mislabel), and flags any curated reference whose best T2T hit falls
  below POLYMORPHIC_IDENTITY_FLAG as polymorphic/absent-from-CHM13, to be used in the
  competitive set rather than as a decoy. Remaining actions: (a) re-run the script against
  T2T (needs minimap2 + the 3.1 GB reference) to regenerate eve_kmers.rs, which still
  holds the old mislabeled data until then; the run will auto-report whether CHM13 carries
  the K113 allele; (b) add K106, K116, ciHHV-6A/B to CURATED_REFERENCES once their
  accessions are verified.
- Calibrate the near-identical inclusion threshold on ground truth (the identity cutoff
  that routes an EVE to B+C versus C-only).
- Decide the pangenome packaging (full versus T2T-plus-decoy), measured both ways.
- Confirm licensing for any resource used (RepBase is commercial; Dfam, HPRC, RVDB,
  Virosaurus, RefSeq are open), and record it in the manifest.
- Fix the biometal and superbloom paths in Cargo.toml before any build or test.
