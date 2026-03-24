# Lab Contaminant vs Natural Sequence Discrimination

## The Problem

Virome sequencing commonly encounters sequences that are both real lab contaminants AND natural members of viral communities. The canonical example is PhiX174:

- **As lab contaminant**: Illumina spike-in control, present at 0.1-1% in every sequencing run
- **As natural virus**: Microviridae phage, genuinely present in human gut viromes

Naive k-mer screening against the PhiX174 genome removes both -- discarding real viral diversity. ViroForge + virome-qc validation showed 1,869 false positives from natural gut Microviridae (coliphages) out of 2,007 "PhiX" detections (93% false positive rate).

Other lab contaminants with the same challenge:

| Contaminant | Lab origin | Natural relatives |
|---|---|---|
| PhiX174 | Illumina spike-in | Microviridae (gut, soil, marine) |
| pUC19 | Cloning vector | Naturally occurring mobilizable plasmids |
| Lambda phage | Cloning vector | Siphoviridae (gut, environmental) |
| T4/T7 phage | Expression systems | Myoviridae/Podoviridae |
| PhiX174 | Control DNA | Sinsheimervirus genus |

## Solution: Unique k-mer subtraction

For each lab contaminant, compute the set of k-mers unique to the lab strain that are NOT found in any natural relative. Only these unique k-mers are used for screening.

### PhiX174 implementation (completed 2026-03-22)

- PhiX174 genome: 5,386 bp, 10,732 21-mers (both strands)
- 20 other Microviridae in RefSeq: 154,784 21-mers
- Shared k-mers: 334 (3.1% of PhiX)
- PhiX-unique k-mers: 10,398 (96.9%)
- Genome coverage: 98.9% of PhiX positions covered by unique k-mers

Result: reads from the actual PhiX174 lab spike-in are correctly removed (96.9% k-mer sensitivity), while reads from natural Microviridae coliphages pass through unaffected.

### Extension to other contaminants

The same approach applies to any lab contaminant with natural relatives:

1. Identify the contaminant's taxonomic family
2. Collect all RefSeq genomes in that family
3. Compute k-mers shared between the contaminant and its relatives
4. Exclude shared k-mers from the contaminant index
5. Validate with ViroForge synthetic datasets containing both the contaminant and its relatives

### For pUC19/vectors

pUC19 derives from pBR322, which derives from natural E. coli plasmids. Some k-mers will match natural mobilizable elements. The same subtraction approach would prevent false positives on natural plasmid-containing bacteria in metagenomics data.

### For Lambda/T4/T7

These are common in E. coli expression work but also have natural relatives. Less critical for virome QC since they're not Illumina spike-ins, but relevant for labs that use these phages in molecular biology workflows adjacent to virome sequencing.

## Validation framework

ViroForge enables systematic validation of contaminant discrimination:

1. Generate a dataset WITH the lab contaminant at known abundance
2. Include the contaminant's natural relatives in the viral community
3. Run virome-qc and compute:
   - Sensitivity: fraction of true contaminant reads correctly removed
   - Specificity: fraction of natural relative reads correctly preserved
   - Precision: true contaminant reads / all reads flagged as contaminant
4. Sweep the containment threshold to find the optimal operating point

This is the first tool to quantitatively address the PhiX-vs-Microviridae problem in virome analysis.
