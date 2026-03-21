# Host Depletion Design Document

## Problem statement

Virome QC must remove host reads before downstream analysis. Traditional host depletion is binary (map to host, remove if mapped), but this is wrong for virome analysis because:

- **T2T-CHM13** resolves centromeric/pericentromeric regions loaded with HERVs and other viral integrations. More complete genome = more viral-like sequence in the reference = more true viral reads incorrectly removed.
- **Endogenous viral elements (EVEs)** span dozens of viral families: HERVs (~8% of genome), endogenous bornaviruses, parvoviruses, filoviruses, chromosomally integrated HHV-6 (~1% of people).
- **Conserved protein domains** shared between human and viral genes (kinases, polymerases, helicases) can cause real viral reads to map to host.
- **The question is not "does this read map to host" but "is this read more likely host or viral".**

## Design: classification not depletion

Three-way classification instead of binary removal:

| Classification | Criteria | Action |
|---|---|---|
| **Host** | Strong host signal, no viral signal | Remove |
| **Ambiguous** | Host signal AND viral signal (EVE/integration regions) | Flag for review |
| **Not host** | No significant host signal | Keep |

Ambiguous reads get flagged in the passport with mapping context (host region, viral family similarity), letting downstream analysis decide.

## Approach: dual-scoring

For each read:
1. Score against host reference (k-mer containment or alignment)
2. Score against viral reference (k-mer containment)
3. Classify based on relative scores

## Key decisions

### Reference genome
- GRCh38 vs T2T-CHM13 vs HPRC pan-genome
- T2T is more complete but contains more EVEs (higher false removal risk)
- Need to evaluate tradeoffs empirically

### EVE handling
- Catalog known EVE regions in host reference
- Reads mapping to EVE regions get dual-scored rather than automatically removed
- EVE catalog needs to be updatable as new integrations are discovered

### Index approach
Options evaluated:
- Full k-mer FxHashSet: ~48GB for human genome (too large)
- Minimizer subsampling (1/100): ~480MB (manageable)
- Bloom filter (1% FPR): ~3.6GB; (10% FPR): ~1.8GB
- Biometal FM-index + seed mapping: unknown memory footprint
- External aligner (minimap2): proven but adds dependency

### Multiple hosts
Config already supports human, mouse, rat, bat, etc. Host index needs to be organism-specific and downloadable.

---

## Experiment log

### Experiment 1: Biometal mapping vs minimap2

**Goal**: Validate that biometal's mapping primitives produce comparable results to minimap2 for the specific task of host read identification.

**Approach**: Take a subset of reads from a real dataset, map with both tools, compare classifications.

**Status**: Complete

**Setup**:
- Reference: Human chr22 (hg38, 50.8 Mbp)
- Reads: 10,000 reads from Buddle WGS (ERR13480651, NovaSeq 6000, CpG-depleted)
- minimap2 v2.30 with `-x sr` (short reads)
- biometal FmIndex with default sampling (sample_rate=32)

**Results**:

| Metric | minimap2 | biometal |
|---|---|---|
| Mapped reads | 1,416 (14.2%) | 1,182 (11.8%) |
| Index build time | n/a (on-the-fly) | 3.5s |
| Mapping time (10K reads) | 0.7s | 0.13s |
| Throughput | ~14K reads/sec | ~77K reads/sec |

**Analysis**:
- biometal maps 234 fewer reads (16.5% fewer than minimap2). These are likely lower-confidence alignments that minimap2 catches with more sensitive seeding.
- biometal is 5.5x faster than minimap2 for the mapping step (77K vs 14K reads/sec)
- FM-index build takes 3.5s for chr22 (50MB). For full human genome (~3GB), expect ~3-5 minutes build time.
- biometal MAPQ thresholds: 804 reads (8.0%) at MAPQ>=30. This is the high-confidence host subset.

**Conclusions**:
1. Biometal mapping is viable for host depletion -- slightly less sensitive than minimap2 but much faster
2. The sensitivity gap (234 reads / 2.4% of total) is acceptable for QC purposes where speed matters more than catching every last host read
3. For a dual-scoring approach (host score + viral score), the speed advantage of biometal is significant
4. FM-index memory for full human genome: ~5-7GB with default sampling, ~2-3GB with memory_optimized

### Detailed concordance analysis

| Category | Count | Notes |
|---|---|---|
| Both mapped | 1,004 | Core agreement |
| minimap2-only | 315 | biometal misses |
| biometal-only | 178 | minimap2 misses |

**minimap2-only reads (315) -- MAPQ distribution:**
- MAPQ 0: 64 (20%) -- no confidence
- MAPQ 1-9: 222 (70%) -- very low confidence
- MAPQ 10-29: 27 (9%) -- moderate
- MAPQ 30+: 2 (0.6%) -- only 2 high-confidence misses

**90% of the "sensitivity gap" is MAPQ < 10 reads.** These are multimapping, repetitive, or marginal alignments that would be filtered by any MAPQ threshold.

**biometal-only reads (178) -- MAPQ distribution:**
- MAPQ 0: 44 (25%)
- MAPQ 60: 134 (75%) -- high confidence

**biometal finds 134 high-confidence mappings that minimap2 misses entirely.** This is unexpected and needs investigation -- could be due to different seeding strategies or alignment heuristics.

**Concordance at MAPQ >= 10:**
- Both agree: 1,004
- minimap2-only: 29
- biometal-only: 134
- **At quality-filtered thresholds, biometal is actually MORE sensitive** (1,138 vs 1,033)

**Key insight**: The raw mapping rate comparison (14.2% vs 11.8%) is misleading. The apparent sensitivity gap is dominated by low-confidence (MAPQ < 10) alignments that minimap2 reports but biometal doesn't. At meaningful quality thresholds (MAPQ >= 10), biometal maps more reads than minimap2.

This is ideal for host depletion where we want confident classifications, not exhaustive low-quality mappings.

**Next steps**:
- Investigate the 134 biometal-only MAPQ 60 reads -- are these true host reads?
- Test on full human genome (not just chr22)
- Test k-mer containment approach as a fast pre-filter
- Design the dual-scoring (host vs viral) classification system

---

### Experiment 2: K-mer containment vs FM-index mapping

**Goal**: Evaluate whether k-mer containment screening (like our contaminant module) can replace full mapping for host identification, trading precision for speed and memory.

**Status**: Planned
