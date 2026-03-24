# ViroForge + virome-qc Integration

## Overview

ViroForge (synthetic virome generator) and virome-qc (virome QC platform) form a closed-loop system for data-driven virome bioinformatics development. ViroForge provides controlled ground truth; virome-qc provides the QC pipeline under test. Together they enable quantitative optimization that neither tool can achieve alone.

## What we accomplished (2026-03-22)

### The validation loop in practice

1. **Generated** a ViroForge gut virome (collection 9, 134 genomes, NovaSeq, 5x, realistic contamination)
2. **Ran** virome-qc on the synthetic data
3. **Compared** virome-qc output against ViroForge ground truth
4. **Discovered** PhiX overcounting: 2,007 detected vs 138 actual (93% false positive rate)
5. **Diagnosed** root cause: k-mer cross-reactivity with natural Microviridae
6. **Implemented** fix: unique k-mer subtraction (334 shared hashes excluded)
7. **Identified** ViroForge curation issue: phiX174 assembly included as natural gut phage

This single cycle produced a concrete improvement to virome-qc (PhiX false positive elimination) and identified a curation error in ViroForge -- demonstrating bidirectional value.

### Findings from the first synthetic validation

| Module | Ground Truth | Detected | Issue Found | Fix |
|---|---|---|---|---|
| PhiX | 138 reads | 2,007 | Microviridae cross-reactivity | Unique k-mer subtraction |
| rRNA | 2,714 reads | 1,223 (45%) | Consensus k-mers insufficient | SILVA filter (98.3% sensitivity) |
| Host DNA | 345 reads | not tested | Host filter not loaded | Validated separately on Buddle WGS |
| Survival | 99.37% | 99.16% | Close match | N/A |
| GC | 28.8% | GC_DEVIATION flag | Correctly flagged AT-rich crAssphage | Working as designed |

## How ViroForge advances virome-qc

### 1. Module-level sensitivity/specificity optimization

For each QC module, generate ViroForge datasets that sweep the parameter of interest:

**Adapter detection**:
- Generate datasets with `--adapter-rate` at 0.01, 0.05, 0.10, 0.20, 0.50
- Vary `--adapter-type` (truseq, nextera)
- Measure: reads correctly trimmed / reads with adapter (sensitivity)
- Measure: reads incorrectly trimmed / reads without adapter (false positive rate)
- Find optimal `min_overlap` and `max_mismatch_rate` per adapter type

**rRNA screening**:
- Generate datasets with varying rRNA fractions (0.01, 0.05, 0.10, 0.30)
- Include diverse rRNA from ViroForge's bundled references (23 species)
- Compare consensus screener (11K k-mers) vs SILVA filter (165M k-mers)
- Determine optimal `min_kmer_fraction` threshold

**Host depletion**:
- Generate datasets with host DNA at 0.01, 0.05, 0.10, 0.30, 0.60, 0.90
- Use real T2T-CHM13 fragments (bundled in ViroForge)
- Sweep `host_threshold` and `ambiguous_threshold`
- Compute ROC curves for host classification

**Complexity filter**:
- Generate datasets with AT-rich phages (crAssphage 29% GC) and GC-rich phages
- Measure false positive rate on real viral sequences
- Optimize `min_entropy` threshold per GC range

**Deduplication**:
- Generate datasets with `--amplification rdab` (MDA) and `--amplification none`
- Known duplication rates from amplification modeling
- Validate dedup module sensitivity vs ground truth duplicate count

### 2. Contaminant discrimination validation

The PhiX/Microviridae finding generalizes to any lab contaminant with natural relatives:

- **pUC19 vs natural plasmids**: Generate dataset with both pUC19 contamination and natural mobilizable elements
- **Lambda vs Siphoviridae**: Test if Lambda screening cross-reacts with natural lambda-like phages
- **T7 vs Podoviridae**: Same pattern

For each, compute unique k-mer exclusion sets and validate precision/recall.

### 3. Regression testing

Generate a canonical ViroForge dataset with fixed seed. On every virome-qc code change:
1. Run virome-qc on the canonical dataset
2. Compare passport JSON against the baseline
3. Any unexpected change in read counts flags a regression

This can be integrated into CI/CD (GitHub Actions).

### 4. Profile validation

Generate datasets matching each virome-qc profile:
- Gut VLP (collection 9) -> stool-vlp-tagmentation profile
- Metagenomics (collection 9, `--no-vlp`) -> metagenomics-nextera profile
- RNA virome (collection 21/23) -> need RNA virome profile
- Low-biomass (collection 16, `--amplification mda`) -> low-biomass-wga profile

Verify that profile thresholds (survival, host fraction, rRNA fraction) produce sensible PASS/WARN/FAIL tiers.

### 5. Report showcase

Generate visually compelling datasets for publication figures:
- High adapter contamination (short inserts)
- Bimodal GC distribution (mixed community)
- High duplication from MDA
- Before/after host depletion
- rRNA fraction comparison (depleted vs undepleted)

## How virome-qc advances ViroForge

### 1. Collection curation validation

Running virome-qc on ViroForge collections exposes curation issues:
- **PhiX174 in gut collection**: virome-qc's contaminant screener flagged phiX174 as a lab contaminant, revealing that GCF_000819615.1 (a PhiX174 assembly) was incorrectly included as a natural gut phage
- **Adapter sequence accuracy**: virome-qc's adapter auto-detection validates that ViroForge's injected adapters match real Illumina sequences
- **Platform artifact fidelity**: virome-qc's poly-G detection validates that ISS's NovaSeq error model produces realistic 2-color artifacts

### 2. Contamination model validation

virome-qc's contaminant detection rates tell ViroForge whether its contamination model is realistic:
- If virome-qc detects more rRNA than ViroForge planted, the contamination abundance is underestimated
- If virome-qc detects less, either the reference sequences are wrong or the k-mer threshold needs tuning
- The 45% rRNA sensitivity finding (consensus k-mers) vs 98.3% (SILVA) tells ViroForge that its bundled rRNA references are detectable but the 23-sequence set doesn't cover all diversity

### 3. Read-level ground truth gaps

ViroForge currently has some gaps that virome-qc testing exposed:
- ISS read headers encode source genome but not artifact type (adapter vs clean)
- No per-read contamination labels (rRNA read vs host read vs viral read)
- Insert size is fixed (not drawn from a distribution)
- No optical duplicate modeling (only PCR duplicates)

virome-qc's module-level analysis reveals which of these gaps matter most for realistic benchmarking.

### 4. Platform-specific artifact validation

virome-qc's ingestion engine (platform detection, Q-score binning, per-position quality profile, poly-G detection) provides a checklist for ViroForge's platform simulation:
- Does ViroForge's NovaSeq model trigger virome-qc's poly-G detection? If not, the artifact isn't realistic enough.
- Does ViroForge's quality degradation produce the expected quality profile shape? virome-qc's per-position scan validates this.
- Does ViroForge's NextSeq model show N-base enrichment at position 0? This is a known NextSeq artifact that virome-qc detects.

## Joint publication potential

The bidirectional validation story has two natural publication angles:

### Application note: virome-qc

- 12 real datasets + ViroForge synthetic validation
- rRNA screening: 98.3% sensitivity with SILVA k=21 filter (validated against ribodetector)
- PhiX discrimination: unique k-mer subtraction eliminates Microviridae false positives
- Host depletion: Super Bloom T2T-CHM13 filter validated end-to-end
- Cross-platform analysis: 6 platforms, GC as diagnostic, library prep signatures

### Methods paper: ViroForge as benchmarking framework

- Demonstrate the closed-loop optimization methodology
- Show how synthetic data identified a real tool bug (PhiX overcounting)
- Show how tool testing identified a database curation error (phiX174 in gut collection)
- Propose the framework as a standard for virome tool development

## Completed work (2026-03-23)

### Bidirectional validation totals: 21 bugs found and fixed
- **11 virome-qc bugs**: PhiX/Microviridae FP, truncated PhiX, adapter double-counting, adapter FP from homology, contaminant threshold, complexity clamp, internal adapter auto-enable threshold, ERV flagging disposition, MinHash empty sketches, low ERV detection, over-clustering, herpesvirus ERV cross-reactivity
- **10 ViroForge bugs**: PhiX174 in gut collection, adapter stats mismatch, ISS no low-complexity, ISS platform models, no ERV injection, RNA template switching crash, viral fraction reporting, CLI missing flags, RNA abundance mismatch, null taxonomy, synthetic RNA contamination references

### Reference atlas: 12 ViroForge datasets
Generated and processed through virome-qc: gut (clean/realistic/heavy/bulk), oral, respiratory RNA, fecal RNA, marine, soil, wastewater, blood/plasma, IBD gut. Establishes expected QC metric ranges per sample type.

### fastp comparison: 3 datasets
virome-qc keeps more viral reads on clean data (+4.3%) and removes more contaminants on dirty data (-6.8%).

### Cross-reactivity analysis: 40 virus families scanned
Only Orthoherpesviridae has clinically significant ERV k-mer cross-reactivity. Fixed with 194-hash exclusion set. Poxviridae (marginal), Phycodnaviridae (low), and Schizomimiviridae (low) documented as known behaviors.

## Remaining next steps

1. **Embed reference ranges** in virome-qc profiles (from ViroForge atlas)
2. **Final benchmark regeneration** with latest code for all 11 real datasets
3. **CI integration**: Add ViroForge regression test to virome-qc GitHub Actions
6. **Build pUC19 exclusion set**: Same unique k-mer approach for vector screening
