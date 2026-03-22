# Phase 4: Deduplication Design Document

## Problem

PCR duplicates inflate abundance estimates and waste downstream compute. But aggressive dedup destroys real signal from abundant viruses in virome data. Position-based dedup (picard/samtools) requires a reference genome, which doesn't exist for virome dark matter.

## Design: sequence-based + optical + UMI

### Virome-aware duplicate detection

PCR duplicates are exact copies (same sequence AND correlated quality scores from the same sequencing cluster). Natural duplicates from abundant organisms have the same sequence but independent quality scores.

**Key virome-aware feature**: For paired-end data, require BOTH mates to be identical to call a duplicate. Two independent fragments from an abundant phage might share R1 start position, but having identical R1 AND R2 (same insert boundaries) is exponentially less likely. This preserves abundant viral signal.

### Three dedup modes

1. **Sequence-based**: Hash read pairs by R1 prefix + R2 prefix (or full sequence for SE). Group identical hashes. Keep highest-quality representative.

2. **Optical**: Within duplicate groups, parse Illumina tile/x/y coordinates from read names. Reads within optical_distance pixels on same tile = optical duplicate (platform artifact, always remove).

3. **UMI-aware**: If UMIs present in read names, group by UMI first. Same UMI + same sequence = PCR duplicate. Different UMI + same sequence = natural duplicate (keep).

### Architecture

Dedup runs as a **post-pipeline step** (not a streaming QcModule) because it needs to see all reads before deciding which are duplicates. Operates on the clean output files.

```
Pipeline produces: clean_*.fastq.gz
    |
    v
Dedup step (two-pass):
  Pass 1: Hash all reads, build duplicate groups
  Pass 2: Mark duplicates, write deduplicated output
    |
    v
Output: dedup_*.fastq.gz + dedup stats in passport
```

### Memory

Hash table of read signatures. For 100M reads at 16 bytes/hash = 1.6 GB.
For paired-end: hash(R1_prefix_50bp + R2_prefix_50bp) = single u64 hash per pair.

### Read signature

```
Single-end: hash(first_50bp + last_50bp) -> u64
Paired-end: hash(R1_first_50bp + R2_first_50bp) -> u64
```

Using first+last captures both the start and end of the fragment, making natural duplicates less likely to collide.

### Optical duplicate detection

Parse Illumina read name: `@instrument:run:flowcell:lane:tile:x:y`
Optical duplicates: same tile, distance(x1,y1, x2,y2) < optical_distance
- Non-patterned flowcell: optical_distance = 100 pixels
- Patterned flowcell (NovaSeq/NextSeq): optical_distance = 2500 pixels

### UMI extraction

UMIs are typically in the read name after the last colon or in a specific field:
- `@read:UMI` or `@read+UMI`
- xGen UDI-UMI format: embedded in the index sequence

When UMIs are present: reads with the same UMI AND same sequence hash = PCR duplicate. Reads with different UMIs AND same sequence = natural duplicate (keep).

### Output

- Deduplicated FASTQ files
- Dedup statistics added to passport:
  - Total duplicates removed
  - PCR duplicate rate
  - Optical duplicate rate (if detectable)
  - Estimated library complexity
