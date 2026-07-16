# EVE/ERV discrimination redesign

Status: proposal. Supersedes the three-signal read-level classifier described in
ERV_ANALYSIS_DESIGN.md. Builds on the subtraction foundation in EVE_MAP_DESIGN.md.

## Why this document exists

The shipped ERV classifier (`erv_classifier.rs`, `erv_pipeline.rs`) makes an
endogenous/exogenous call from three signals computed on individual short reads or a
short unaligned consensus: ORF integrity, CpG depletion, and MinHash distance. A review
of the code and the literature found this is not the right design:

- Two of the three signals (ORF, CpG) need kilobases of contiguous sequence to be
  meaningful and are near-useless at read length.
- In the shipping configuration the MinHash signal is a constant 0.5 (the ERV and
  exogenous reference panels are the same object when no Dfam panel file is present),
  so classification collapses to a CpG-ratio threshold, silently.
- CpG depletion is retrovirus-specific. Circulating exogenous RNA viruses are often
  CpG-depleted for ZAP evasion (SARS-CoV-2, HCV, alphaviruses), so low CpG does not
  imply endogenous. For ciHHV-6 the composition signals invert: an inherited, intact
  integration reads as exogenous.
- Family routing, the polymorphism k-mer database, and most of the EVE exclusion set
  are scaffolded but run on empty or placeholder data.

The signals the field actually relies on are competitive reference identity, coverage
breadth and evenness, abundance/copy-number relative to host, within-sample diversity,
and host-virus junctions. This document scopes an implementation around those.

## Problem restatement

Two goals in tension:

1. Remove host-derived EVE reads so they do not surface as false-positive viruses in
   downstream classification.
2. Do not remove reads from a genuine exogenous infection that resembles an EVE.

The tension is real only where an endogenous copy and a circulating exogenous virus are
close relatives. For most cases they are far apart and the problem is easy:

- HIV has no endogenous relative in the human genome (no endogenous lentiviruses), so
  it is never at risk from EVE-aware host removal.
- Ancient nrEVEs (endogenous bornavirus, filovirus, parvovirus elements) are tens of
  millions of years old and have no close circulating relative, so a real infection
  maps to distant references and the EVE does not hide it.

The hard cases are recent integrations with a near-identical circulating form:

- ciHHV-6: a full-length HHV-6 genome inherited at a telomere in about 1% of people,
  greater than 99.9% identical to the reactivating exogenous virus.
- Cross-species or contaminant gammaretroviruses (mouse MLV, XMRV): endogenous and
  exogenous copies are nearly identical.

For the hard cases, sequence identity alone cannot separate endogenous from exogenous.
The discriminators are structural and quantitative: copy-number tracking the host
genome, the germline integration junction, and coverage evenness.

## Approach: subtract, then discriminate the residual

Two stages, matching the empirical finding in the tale-of-caution paper that host+EVE
subtraction alone is a blunt instrument (it removes real viral hits in EVE-similar
regions; they lost 213 viral contigs to over-subtraction).

Stage 1, subtract: remove reads that match the host genome including its endogenous
viral content, using a complete reference. This handles the false-positive goal for
everything that is sequence-separable, which is most EVEs.

Stage 2, discriminate the residual: for reads that survive Stage 1 and carry viral
signal, compute evidence along the axes below and report it. Do not render a verdict
where the evidence is ambiguous; flag it.

## Signals and the unit each operates on

Stating the unit is the discipline whose absence sank the original classifier. Read-
level bit-scores and ORF calls are noise; the published thresholds were calibrated on
contigs and proteins.

| Signal | Unit | What it separates | Metric |
|---|---|---|---|
| Competitive reference identity | contig (assemble first) | close vs distant exogenous relative | delta bit-score between exogenous-panel best hit and EVE/host-panel best hit; tale-of-caution uses <= -10 for EVE-like |
| Coverage breadth and evenness | per-reference, coverage | real full-genome presence vs conserved-domain cross-map | observed breadth vs expected breadth at the observed depth; evenness (for example Gini) of read distribution |
| Abundance / copy-number vs host | per-reference, coverage | germline integration (about 1 copy per cell) vs infection | viral mean depth / host mean depth |
| Host-virus junction | read (soft-clip) | integrated vs free/episomal, and where | split reads with one segment on virus and one on host; one split read is enough to flag |
| Within-sample diversity | per-reference, coverage | fixed host allele vs recent-ancestor infection | SNV density / allele-fraction spread across mapped reads |
| ORF integrity, CpG | contig only | supporting evidence for ancient degraded EVEs | keep as reported context, never as the deciding read-level signal; note both invert for ciHHV-6 |

## Gate discrimination by sample regime

Copy-number-vs-host and junctions are data-starved in host-depleted VLP viromes: there
is little host coverage to form the copy-number denominator and few reads spanning any
junction. This is not a flaw to engineer around, it is a reason to condition. The EVE
false-positive problem is mostly a tissue problem (VLP enrichment removes the host DNA
that carries EVEs), and tissue is exactly where host coverage is available. So the
expensive signals are feasible where they are needed and infeasible where they are not.

Use the existing profile system as the switch:

- `tissue-truseq` and similar high-host profiles: run the full discrimination stack
  including copy-number and junctions.
- `stool-vlp-tagmentation` and similar low-host profiles: run competitive identity and
  breadth only; report copy-number and junctions as not-assessable rather than
  computing an unstable value.

## Tooling: separate load-bearing measurements from primitives

The measurements that decide endogenous vs exogenous, and that a validation or an
independent publication would rest on, are: coordinate alignment, coverage breadth and
depth, split-read junctions, and competitive identity. These must rest on established,
benchmarked, citable tools. biometal has an in-house `mapping`, `pileup`, and `assembly`
module and using them would avoid a new dependency, but that is the wrong trade here:
biometal is unpublished and unbenchmarked, so building the load-bearing measurements on
it means validating the aligner and the method at the same time, with no way to tell
which one is responsible for a wrong call and no independent tool a reproducer or
reviewer already trusts. That confounds the exact thing this design is trying to
validate cleanly.

Use established tools for the load-bearing stage:

- Alignment, coverage breadth, depth, soft-clip junctions: minimap2 (short-read preset)
  or bwa-mem2, with coverage computed from the alignment. EVE_MAP_DESIGN.md already
  specifies minimap2 for the T2T homology map, so this is consistent. Pure-Rust
  integration is possible through minimap2 FFI bindings; confirm the binding crate and
  its API before committing, otherwise shell out to the binary in the post-pipeline
  discrimination stage (it is not on the streaming hot path).
- Assembly of flagged reads for the contig-level signals: an established assembler
  (metaSPAdes or MEGAHIT), not a first-pass in-house de Bruijn graph.
- Host subtraction: the benchmarked-best short-read host removers are bwa-mem2 or
  bowtie2 against T2T plus pangenome. SuperBloom k-mer containment can stay as a fast
  first pass, but validate its host-removal calls against bwa-mem2 on ground truth
  before trusting them for the false-positive claim.

Keep biometal where it is not carrying a novel quantitative claim: FASTQ I/O, trimming,
GC and complexity, the k-mer containment first-pass flag, and MinHash sketching for
initial panel assignment. These feed the streaming QC pipeline and a coarse first-pass
flag, and a wrong value there is caught by the established-tool stage that follows.

## Build order

Ordered by correctness and dependency, not by any release or manuscript need. Each step
is a prerequisite for the ones after it. Getting the analysis correct is the goal; a
clean validation of it (next section) could stand as an independent contribution.

### Step 0, decommission the current classifier

The shipped three-signal classifier produces calls its inputs do not support (MinHash
inert by construction, CpG doing the work, family routing and polymorphism k-mers on
empty data). Before building the replacement, stop the current one from emitting numbers
that read as verdicts: either remove it from the run path or reduce its output to a
first-pass flag with no endogenous/exogenous label. This is cleanup, not a gate.

### Step 1, curate the reference sets (foundational)

Most of the scientific correctness lives here: what is in the panels determines what can
be discriminated. Everything downstream depends on this being right. Detailed spec in
REFERENCE_CURATION.md.

- Pangenome host set: HPRC assemblies plus T2T, for subtraction.
- EVE/ERV decoy and consensus set: curated with detectEVE or EEfinder, including the
  polymorphic HERV-K loci (K113, K106, K116), ciHHV-6 reference genomes, and the nrEVE
  loci already scaffolded. Used both as decoys in subtraction and as one side of the
  competitive-identity comparison.
- Exogenous panel: curated, host-matched, small. For human: HIV-1/2, HTLV-1/2,
  reactivating HHV-6. For mouse: MLV, MMTV, XMRV. This is close to a closed set.

### Step 2, subtraction (Stage 1)

Remove host including EVE content with an established host remover (bwa-mem2 or bowtie2
against the Step 1 pangenome set), or keep SuperBloom as a fast first pass validated
against bwa-mem2. Do not over-subtract: route partially matching reads to a separate
ambiguous output for Stage 2 rather than deleting them, consistent with the current
three-way host classification.

### Step 3, read-level first-pass flag (Stage 2 entry)

Flag residual reads carrying viral signal with the biometal k-mer containment first pass
and MinHash panel assignment. Coarse and non-load-bearing; every flagged read is
re-examined by the established-tool stage that follows.

### Step 4, alignment-based discrimination (load-bearing)

- Map flagged reads to the exogenous panel and the EVE/host panel with minimap2 or
  bwa-mem2.
- Compute breadth, expected breadth at the observed depth, evenness, and depth from the
  alignment.
- Extract soft-clipped reads and test whether the clipped segment maps to host, to flag
  junctions (one split read is enough at low coverage).
- Gate copy-number and junctions by profile per the regime section.

### Step 5, assembly and competitive identity (load-bearing, contig level)

- Assemble flagged reads with an established assembler (metaSPAdes or MEGAHIT).
- Compute delta bit-score competitive identity at contig level against both panels.
- Only here do ORF integrity and CpG become legitimate supporting context, and only as
  supporting context, never as the deciding signal.

### Step 6, evidence-reporting contract

Emit the evidence schema below. Build alongside Steps 4 and 5 as the fields become
available. This is the product: report evidence, flag ambiguity, decide only what the
data supports.

### Out of scope, hand off

- Full breakpoint resolution and divergent-virus integration calling (ViFi, geNomad).
  The tool should surface junction evidence, not become an integration caller.
- Somatic integration clonality (HBV in HCC, HPV in cervical cancer). Different regime
  (subclonal allele fraction), usually a real finding to keep, better served by
  ViFi-style tools on assembled tumor contigs.
- Contig-level provirus/host boundary demarcation (geNomad plus CheckV already do this
  well). Emit contigs and reference bundle for the user to run these.

## Passport evidence schema

Per candidate reference hit, report evidence rather than a label:

- reference name, family, panel (exogenous or EVE/host)
- reads mapped, mean depth
- breadth observed, breadth expected at this depth, evenness score
- depth vs host (or not_assessable when host coverage is insufficient)
- competitor margin (delta bit-score vs the other panel; not_assessable if no assembly)
- junction reads (count; not_assessable in low-host regime)
- within-sample diversity (SNV density; optional)
- interpretation flag: consistent_with_infection, consistent_with_endogenous, or
  ambiguous_needs_orthogonal_data, with the reason string

The `ambiguous` flag is a first-class outcome, not a failure. The honest message for the
hard middle is "consistent with either, needs long reads or matched host WGS or a
targeted assay."

## Validation

This is the point of the work, and a clean version of it is what could stand alone.
Validate the method, not the tools under it, which is the reason the load-bearing stage
uses established aligners and assemblers.

Ground-truth datasets, spanning the easy and hard cases:

- Synthetic (ViroForge): plant known EVE reads and known exogenous reads at controlled
  ratios and coverage, with labels. Sweep the endo/exo identity gap from distant (HIV
  vs HERV) to near-identical (recent HERV-K, MLV) and sweep exogenous coverage from
  low-titer to high. This maps the confidence surface directly.
- ciHHV-6 real samples: the definitive hard positive. Expect about 1 copy per cell,
  whole-genome breadth, and a telomeric junction. The method should call endogenous
  where the shipped composition signals call exogenous.
- Mouse MLV / XMRV: the near-identical endo/exo case, ideally including the known
  mouse-DNA-contamination scenario that produced false XMRV calls.
- Exogenous-positive controls: HIV-positive plasma or spike-ins, where the answer is
  unambiguously exogenous and must not be subtracted.
- Negative controls: host-only tissue with no infection, where every viral-looking call
  must resolve to endogenous or be flagged, not called as an infection.

Metrics:

- False-positive virus rate: endogenous reads called as exogenous infection. The primary
  thing this whole design exists to reduce.
- False-negative infection rate: real exogenous reads lost to subtraction or called
  endogenous. The counter-goal, in tension with the above.
- Endo/exo call accuracy stratified by identity gap and by coverage, to show where the
  method is confident and where it correctly returns ambiguous.
- Calibration of the ambiguous flag: ambiguous calls should concentrate exactly in the
  low-identity-gap, low-coverage corner where no short-read method can decide.

Compare against the shipped three-signal classifier on the same data to quantify the
improvement, and against simple baselines (host subtraction alone, competitive identity
alone) to show each added signal earns its place.

## Prerequisites and open questions

- Build path: `Cargo.toml` points biometal and superbloom at `/Users/scotthandley/...`,
  but on this machine biometal is at `/Users/shandley/Code/software/biometal`. Fix the
  paths (or use a path override) before any build or test work.
- Confirm the minimap2 integration path: whether a maintained Rust FFI binding exists
  with a usable API, or whether the discrimination stage shells out to the binary. This
  decides how the load-bearing alignment stage is wired.
- Decide the pangenome packaging: full HPRC assembly set, or T2T plus a curated decoy
  FASTA. The decoy route is lighter and probably sufficient for the polymorphic-
  insertion gap; measure filter size, build time, and false-positive reduction both
  ways on ground truth.
- Recalibrate any threshold borrowed from contig-level literature (delta bit-score
  <= -10, ORF < 200 aa) before applying it, or assemble first so it applies as
  published.
- External dependency policy: the load-bearing stage adds minimap2 and an assembler as
  runtime dependencies for the post-pipeline discrimination step. Decide whether that is
  acceptable for the tool's distribution model (the streaming QC path stays pure Rust;
  the analysis stage already assumes external reference databases).

## References

- A tale of caution: how endogenous viral elements affect virus discovery in
  transcriptomic data. Virus Evolution 2024.
- Pangenome databases provide superior host removal and mycobacteria classification.
  GigaScience 2024.
- Use of the CHM13-T2T genome improves metagenomic analysis. mSystems 2025.
- Metapresence: species detection based on genome-wide distribution of mapping reads.
  mSystems 2024 (observed vs expected breadth, evenness).
- Virus-Clip; ViFi; VERSE (soft-clip and split-read integration detection).
- Inherited chromosomally integrated HHV-6 genomes are ancient, intact, telomeric.
  J Virol 2017.
- ZAP-driven CpG suppression in exogenous RNA viruses (SARS-CoV-2, HCV, alphaviruses).
- detectEVE; EEfinder (EVE curation); geNomad; CheckV (contig-level hand-off).
