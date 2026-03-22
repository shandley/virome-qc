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

### Experiment 1b: Investigation of 134 biometal-only MAPQ 60 reads

**Finding**: These are NOT true host reads. They are short-insert adapter-contaminated reads.

- 115 of 178 (65%) contain TruSeq adapter sequence (`AGATCGGAAGAGCACACGTCTGAAC`)
- The genomic portion is only 30-80bp, rest is adapter read-through
- BLAST confirms: best alignments cover only 28-79bp of 150bp reads
- Biometal seeds on the short genomic fragment and calls it MAPQ 60
- Minimap2 default mode correctly rejects (insufficient aligned fraction)
- With sensitive settings (`-k 13 -w 5`), minimap2 finds 198 partial alignments

**Impact on pipeline design**: This is not a real problem because:
1. Adapter trimming (Module 1) runs before host depletion
2. After trimming, these reads become 30-80bp fragments
3. Most will fail the length filter (Module 7)
4. Those that survive are short genuine genomic fragments that should map correctly

**Biometal limitation to note**: MAPQ calculation doesn't account for the fraction of the read that aligns. A 40bp alignment out of 150bp should not get MAPQ 60. This could be addressed by adding a minimum aligned fraction threshold to the host depletion module.

### Revised concordance (excluding adapter-contaminated reads)

After removing the 115 adapter-contaminated reads from biometal-only:
- True biometal-only: ~63 reads (from 178 - 115)
- Of these, many are likely short partial alignments similar to the adapter reads

**At MAPQ >= 10, excluding adapter artifacts, both tools have comparable sensitivity for clean reads.** The pipeline ordering (adapter trim -> host mapping) naturally resolves the biometal false positive issue.

---

## Reference genome decision

**Choice: Unmasked T2T-CHM13**

Using T2T over GRCh38 because it resolves all centromeric/pericentromeric regions where many HERVs and repeats reside. More complete reference = better host depletion sensitivity. The increased EVE content risk is handled by the three-way classification, not by using an incomplete reference.

**Rejected: Masked reference approach (hecatomb-style)**

Hecatomb shredded RefSeq viral genomes and masked host regions where they mapped. Problems:
1. Viral database is incomplete -- unmasked novel EVEs still cause false removal
2. Masking creates alignment gaps -- reads spanning mask boundaries get poor scores
3. Static masking goes stale when RefSeq updates
4. Over-masking risk -- conserved domains (polymerases, kinases) shared between human genes and viruses would mask real coding regions
5. Conflates EVE regions (ancient integrations) with shared domains (convergent evolution)

The three-way classification approach handles ambiguity in the output (better) rather than modifying the input (fragile).

## Retrovirus handling: endogenous vs exogenous discrimination

### The problem

ERVs (~8% of human genome) are degraded viral fossils. Reads from ERV loci look viral but are host. However, some retroviruses (HHV-6, HTLV, HIV) can be actively infectious AND share sequence with endogenous elements. Binary classification fails here.

### SNP-level discrimination

**Endogenous reads** match the host reference perfectly at the ERV locus -- they ARE the reference at that position.

**Exogenous/infectious retroviral reads** show systematic divergence:
- Consistent SNPs across multiple reads (not random sequencing error)
- SNPs that preserve open reading frames (functional constraint)
- Divergence level >2-3% from the integrated copy (distinct lineage)
- RT error signature patterns

### Implementation plan (Phase 3b -- separate from host depletion)

This is a region-level analysis, not a per-read filter:

1. Collect all reads mapping to known EVE/HERV regions (from the "ambiguous" output)
2. Pileup and variant calling per EVE locus (biometal has these primitives)
3. Classify variant pattern:
   - No variants vs reference = endogenous (host reads from ERV locus)
   - Consistent SNPs preserving ORFs = possible exogenous retrovirus
   - Mix of perfect-match and variant reads = both present (co-mapping)
4. Report per-locus classification with evidence

This is genuinely novel -- no existing virome pipeline does SNP-level endogenous vs exogenous retrovirus discrimination.

**Separation rationale**: Keep the QC pipeline fast (per-read streaming). The retrovirus analysis runs optionally on the ambiguous reads as a post-hoc deep-dive.

---

## Phase 3 implementation plan

### Module: Host depletion

**Position in pipeline**: After contaminant screening (Module 6), before length filter (Module 7).

### Approach

1. Build FM-index from T2T-CHM13 reference (one-time, ~5GB index)
2. Map each read using biometal's ParallelMapper
3. Three-way classification:

| Read maps to host? | Maps to EVE region? | Viral k-mer signal? | Classification |
|---|---|---|---|
| No | -- | -- | KEEP |
| Yes, high MAPQ | No | -- | REMOVE (host) |
| Yes, high MAPQ | Yes | No | REMOVE (endogenous EVE) |
| Yes, high MAPQ | Yes | Yes | FLAG (ambiguous) |
| Yes, low MAPQ | -- | -- | FLAG (ambiguous) |
| Yes, any | -- | -- | Check aligned fraction >= 70% |

4. Minimum aligned fraction threshold (70%) to prevent the adapter-contamination false positives found in Experiment 1b

### Inputs required

- **T2T-CHM13 reference**: Downloaded via `virome-qc db build` or user-provided
- **EVE annotation BED file**: Coordinates of known EVE/HERV regions in T2T
- **Viral k-mer index**: Already built for contaminant screening module (rRNA/PhiX), extend with RefSeq viral k-mers for dual-scoring

### Outputs

- `clean_*.fastq.gz` -- reads classified as NOT host
- `ambiguous.fastq.gz` -- reads that map to host but with viral signal or low confidence
- Passport: host fraction, ambiguous fraction, EVE region hit distribution

### Experiment 3: Full T2T-CHM13 FM-index scalability test

**Result: FAILED -- OOM killed after 52 minutes**

- Reference: T2T-CHM13 v2.0 (3.117 Gbp)
- Process killed by OS (exit 137) after 52 min wall time
- Peak memory: ~5 GB observed, but suffix array construction requires ~25-30 GB intermediate
- System time: 2142s (35 min!) = severe memory thrashing / swapping

**Root cause**: Biometal's suffix array construction uses O(n * 8) bytes for 64-bit suffix array pointers. For 3.1 Gbp, this is ~25 GB for the SA alone, plus working memory. Exceeds typical workstation RAM.

**Works fine for**: Individual chromosomes (chr22 = 50 MB, builds in 3.3s, 5 GB peak)

### Revised approach options

1. **minimap2 as external dependency** -- proven, builds index in ~8 GB, most virome pipelines use this
2. **Chromosome-by-chromosome FM-index** -- build per-chromosome, map against each. Complex but pure biometal.
3. **K-mer containment screening** -- minimizer-subsampled host k-mer set (~500 MB). No suffix array needed. Less precise.
4. **Optimize biometal SA construction** -- use divsufsort/SA-IS (O(n) space). Biometal improvement.

### Performance targets (revised)

- Mapping approach TBD based on option selection
- Must work within 8-16 GB RAM on typical workstations

### Next steps

1. Test FM-index build on full T2T-CHM13 (memory and time validation)
2. Build EVE annotation BED from gEVE or similar database
3. Implement the host depletion module
4. Validate on Buddle data (known high-host content)
5. Test on Cook data (should remove near-zero reads -- no host)

---

### Experiment 2: K-mer containment vs FM-index mapping

**Goal**: Evaluate whether k-mer containment screening (like our contaminant module) can replace full mapping for host identification, trading precision for speed and memory.

**Status**: Planned
