# virome-qc Validation Log

Running record of real-data validation results. Each dataset tests different aspects of the QC pipeline and reveals areas for improvement.

## Datasets Tested

### 1. Cook et al. 2024 -- 15-phage mock (ERR10359658)

**Source**: Cook R et al. "The long and short of it: benchmarking viromics." Microbial Genomics 2024. ENA: PRJEB56639.

**Characteristics**: Nextera XT library prep, MiSeq 2x250, 15 known phages, GenomiPhi MDA amplification. 1.27M paired-end reads.

**What it tests**: Adapter detection (Nextera), quality trimming on degraded 3' ends, duplication from MDA.

**Results** (virome-qc v0.1.0, metagenomics-nextera profile):

| Metric | Value | Notes |
|--------|-------|-------|
| Survival rate | 91.9% | |
| 3' adapters found | 7,771 (0.31%) | Lower than expected for Nextera XT -- well-controlled insert sizes |
| Internal adapters | 129 (0.01%) | |
| Quality failures | 202K (8.0%) | Expected for MiSeq 2x250 3' degradation |
| - too short after trim | 172,700 | Dominant failure mode |
| - low mean quality | 29,401 | |
| N-base failures | 2,440 (0.1%) | Real failed MiSeq clusters |
| Poly-G trimmed | 10,210 (0.4%) | Real genomic poly-G, not platform artifact (MiSeq is 4-color) |
| Complexity failures | 15 (0.001%) | Negligible |
| Merge rate | 53.6% of passing pairs | Healthy for 250bp reads |
| Est. duplication | 55.9% | Consistent with MDA amplification |
| Pairs passed | 1,150,420 | |
| Singletons | 31,095 | 1.2% of input reads |

**Findings and fixes**:

1. **RC adapter sequences needed for read-through detection** (fixed in b53e3c0). Paired-end read-through produces the reverse complement of the opposite adapter. Adding RC forms tripled 3' adapter detection from 0.10% to 0.31%.

2. **RC sequences must not be used for internal adapter scanning** (fixed in b53e3c0). The shorter RC motifs match genomic sequence at high rates, causing 393x false positive inflation. Internal scanning uses forward-only sequences.

3. **MiSeq poly-G is genomic, not artifactual**. MiSeq uses 4-color chemistry, so poly-G tails are real sequence. The platform_aware feature should ideally detect the platform from FASTQ headers and adjust thresholds, but currently uses a blanket shorter threshold for poly-G regardless.

4. **Quality is the dominant filter for MiSeq 2x250**. The 3' quality degradation on MiSeq 250bp reads is severe (Q14 median at position 240). Most removed reads fail the length filter after aggressive quality trimming.

---

## Datasets Queued

### 2. Kleiner et al. 2017 -- 32-species mock with 5 phages (ERR1877475)

**Source**: Kleiner M et al. Nat Comms 2017; Ho et al. Microbiome 2023. ENA: PRJEB19901.

**Characteristics**: NEBNext Ultra II library prep, NextSeq 500 2x75bp, 32-species mock (23 bacteria, 3 archaea, 1 yeast, 5 phages). 18.3M paired-end reads.

**What it tests**: Poly-G artifacts on real 2-color chemistry, short-read handling, QC with mixed community.

**Results** (virome-qc v0.1.0, short-read-nebnext profile, min_length=40):

| Metric | Value | Notes |
|--------|-------|-------|
| Survival rate | 95.5% | |
| 3' adapters found | 763,419 (2.08%) | All NEBNext -- correctly identified |
| Internal adapters | 0 (0%) | Disabled for short reads (too few internal k-mers) |
| N-filter failures | 753,643 (2.06%) | Significant -- see finding #1 |
| Quality failures | 420,794 (1.15%) | |
| - too short after trim | 151,481 | |
| - low mean quality | 269,313 | |
| Complexity failures | 481,234 (1.31%) | Higher than expected -- see finding #3 |
| Poly-G trimmed | 37,698 (0.10%) | Real 2-color artifacts detected |
| Poly-other trimmed | 41,762 (0.11%) | Interesting -- more non-G than G |
| Merge rate | 54.6% of passing pairs | Good for short reads |
| Est. duplication | 15.9% | Reasonable for non-amplified library |
| Singletons | 47,767 (0.13%) | Low -- good R1/R2 concordance |

**Findings**:

1. **N-bases are a real problem on NextSeq**. 753K reads (2.06%) failed the N-filter -- 10x higher than MiSeq. Per-position base composition shows 2.1% N at position 0, confirming NextSeq cluster initialization issues. Our N-filter is catching real platform artifacts.

2. **Poly-G detection works on real 2-color data**, but the numbers are modest (37K reads, 0.10%). The poly-other count (41K) is actually higher -- suggesting many poly-A/T/C runs in the microbial community sequences, not just 2-color artifacts. May need to separate platform poly-G (artifact) from genomic homopolymers (real) in reporting.

3. **Complexity filter is aggressive for mixed communities**. 481K reads (1.31%) failed at entropy threshold 0.5. This is higher than the Cook phage-only data (0.001%). Mixed microbial communities include rRNA fragments and repetitive genomic regions that are biologically real but low-complexity. The 0.5 threshold may be too high for metagenomics. Consider whether this should be profile-dependent (lower for metagenomics, higher for VLP).

4. **Short-read profile was necessary**. The default min_length=90 would have killed nearly all 75bp reads. The short-read-nebnext profile with min_length=40 worked correctly. This confirms that profiles need read-length awareness.

5. **Per-position quality shows NextSeq Q-score binning**. Quality values cluster at Q32 and Q36 (two of the three NovaSeq/NextSeq quality bins). Median stays at 32-36 across all positions. This is the expected 2-color chemistry quantization.

6. **Adapter detection rate is 20x higher than Cook data** (2.08% vs 0.10%). NEBNext on NextSeq with 2x75 has more read-through than Nextera XT on MiSeq 2x250 -- shorter reads relative to insert size means more adapter overlap.

### 3. Buddle et al. 2024 -- ATCC virome mock (PRJEB74559, single sample)

**What it will test**: NovaSeq 6000 + NextSeq 2000 poly-G artifacts on human virome data. UMI adapter handling (xGen UDI-UMI adapters). High host background (human DNA/RNA). Multiple library prep methods (NEBNext, KAPA HyperPrep, Twist capture).

**Key questions**:
- Does NovaSeq poly-G detection and trimming work correctly?
- Can we handle UMI-containing adapters?
- What does the adapter breakdown look like across different library preps?
- How does host fraction (not yet removed in Phase 1) affect survival rate?

### 4. Wu/Piedade et al. 2024 -- Antarctic paired VLP/microbial (PRJEB71789, single sample)

**What it will test**: TruSeq Nano on NovaSeq 6000. Large-scale environmental virome. Paired viral and cellular fractions from same sample.

**Key questions**:
- Do QC metrics differ between VLP and cellular fractions from the same sample?
- Is TruSeq adapter detection working correctly?
- What does the GC distribution look like for a real environmental virome?

### 5. Lab data (Handley lab)

**What it will test**: Controlled samples where wet-lab details are fully known. Can validate against expected composition.

---

## Patterns to Watch Across Datasets

- Adapter detection rate by library prep (Nextera vs TruSeq vs NEBNext)
- Poly-G rate by platform (MiSeq vs NextSeq vs NovaSeq)
- Quality profile shape by platform and read length
- Duplication rate by amplification method (none vs PCR vs MDA)
- Complexity filter calibration across different viral genome GC ranges
- Insert size distributions from merged pairs
- Singleton rate as an indicator of asymmetric read quality (R1 vs R2)
