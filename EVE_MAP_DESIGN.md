# T2T-CHM13 Viral Homology Map

## Concept

A comprehensive map of all regions in the T2T-CHM13 human reference genome that share sequence homology with known viral genomes. Built by aligning all RefSeq viral sequences against T2T-CHM13 using minimap2.

## Novelty

Previous EVE catalogs (gEVE, HERVd, pEVE) were built against GRCh38, which has gaps in centromeres, telomeres, and pericentromeric regions -- exactly where many HERVs and satellite-associated viral integrations reside. T2T-CHM13 resolves all these regions, providing the first complete viral homology landscape of a human genome.

Additionally, previous catalogs focused on specific EVE classes (HERVs, specific non-retroviral families). This map captures ALL homology between the host genome and the complete RefSeq viral database -- including ancient degraded elements, short fragments, and viral domains acquired through horizontal gene transfer.

## Method

### Step 1: Extract viral sequences from ViroForge database
- 14,423 RefSeq viral genomes (504 Mb total)
- Export as multi-FASTA with taxonomy annotations

### Step 2: Align to T2T-CHM13
```bash
minimap2 -x asm20 -N 1000 --secondary=yes \
  chm13v2.0.fa viral_genomes.fa \
  > viral_to_t2t.paf
```
- `-x asm20`: preset for divergent sequences (up to 20% divergence)
- `-N 1000`: report up to 1000 secondary alignments (EVEs may match multiple viral genomes)
- `--secondary=yes`: keep secondary alignments for multi-copy EVEs

### Step 3: Filter alignments
- Minimum alignment length: 100bp (short fragments may be chance matches)
- Minimum query coverage: 30% (partial homology is expected for degraded EVEs)
- Minimum identity: 50% (ancient EVEs may be highly diverged)
- These thresholds are deliberately permissive -- better to over-call and filter later

### Step 4: Generate BED file
For each alignment:
```
chr  start  end  viral_genome_id  identity  alignment_length  viral_family
```

### Step 5: Merge overlapping regions
Collapse overlapping intervals into EVE regions. Annotate each region with:
- Best-matching viral family
- Number of viral genomes with homology
- Mean and max identity
- Whether the region is a known HERV, endogenous bornavirus, etc.

### Step 6: Validate against RepeatMasker
Compare our EVE map against RepeatMasker's HERV annotations for T2T-CHM13 to assess completeness. Our map should contain all RepeatMasker HERV calls plus additional non-retroviral EVEs.

## Expected results

Based on known biology:
- ~450,000 HERV elements (~8% of genome, ~250 Mbp)
- Endogenous retroviruses: HERV-K (HML-2), HERV-H, HERV-W, HERV-E, etc.
- Non-retroviral EVEs: bornaviruses (~30 loci), parvoviruses (~4), filoviruses (~1)
- ciHHV-6 integration sites (telomeric)
- Possible novel EVE discoveries in T2T-resolved regions

## Output

Two files:
1. `t2t_eve_regions.bed` -- all viral homology regions with annotations
2. `t2t_eve_summary.json` -- statistics, family breakdown, comparison with RepeatMasker

## Use in virome-qc

During host depletion, reads mapping to EVE regions are treated differently:
- Instead of "host" classification, reads in EVE regions get "eve" classification
- The passport reports EVE region hits separately from host hits
- For approach 1 (k-mer exclusion): build exogenous-unique k-mers by subtracting EVE-region k-mers from viral k-mer sets
- For approach 2 (coordinate masking): reads within EVE coordinates flagged rather than removed

## Computational requirements

- minimap2 alignment: ~30-60 min on 8 cores
- Input: 504 Mb viral + 3.1 Gb host = 3.6 Gb total
- Output: PAF file (~100-500 Mb depending on number of hits)
- BED file: ~10-50 Mb
- Memory: ~10 Gb (minimap2 index for T2T)
