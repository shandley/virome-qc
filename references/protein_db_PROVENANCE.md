# EVE / set-D screen references (protein DBs + nucleotide indices)

Built 2026-07-16 for the alignment-based EVE screen. Spec: EVE_SCREEN_MODULE.md,
"Algorithm: alignment-based competitive screen". The protein (DIAMOND/blastx) arm and the
nucleotide (minimap2) arm are documented below.

## EVE-protein DB (host endogenous)  ->  eve_prot.dmnd

- Source: the 21,302 internal HERV loci (proviruses + internal-only; solo LTRs excluded,
  they encode nothing) from the validated full-HERV reference
  set_b_hervfull.fa (99,035 loci) on HTCF /scratch/sahlab/shandley/eve_exclusion.
  Families ERV1 (9,993) / ERVK (923) / ERVL (10,386). 88 Mb nucleotide.
- Method: 6-frame translation, split each frame on stop codons, keep peptides >= 50 aa
  (scripts/translate_orfs.py). Annotated CDS do not exist for degraded endogenous loci, so
  translation of the genomic sequence is the honest, reproducible source; the resulting
  fragmented ORFs are fine as blastx targets.
- Result: eve_prot.faa = 657,297 peptides, 48.8M aa letters. eve_prot.dmnd = 97 MB.
- diamond 2.2.4.
- TODO on detectEVE completion: append translated nrEVE ORFs (EBLN/EPV) and rebuild in one
  pass, so the protein DB matches the full nt reference.

## set-D-protein DB (exogenous panel)  ->  set_d_prot.dmnd

- Source: NCBI-annotated CDS proteins (efetch rettype=fasta_cds_aa) for 14 accessions --
  clean, correctly-bounded gag/pol/env and RdRP, not 6-frame guesses:
  HHV-6A NC_001664.4, HHV-6B NC_000898.1, MLV Moloney NC_001501.1, MLV Friend NC_001362.1,
  Murine type C NC_001702.1, MMTV NC_001503.1, Mammalian-1 orthobornavirus NC_030692.1,
  VSBV-1 NC_030701.1, parvovirus B19 NC_000883.2, AAV-2 NC_001401.2, Ebola NC_002549.1,
  Marburg NC_001608.4, HIV-1 NC_001802.1, HTLV-1/STLV-1 NC_000858.1.
- Verified present: MLV/MMTV/HIV/HTLV gag/pol/env, herpesvirus polymerase + envelope,
  bornavirus/filovirus L (RdRP).
- Result: set_d_prot.faa = 262 proteins, 127,785 aa letters. set_d_prot.dmnd = 183 KB.

## Functional discrimination test (diamond blastx, e<=1e-5, best hit)

| query | vs EVE-prot | vs set-D-prot | delta = exo-EVE | call |
|---|---|---|---|---|
| HERV-K internal locus (endogenous) | bit 674, 100% id | bit 272, 34.9% id (MMTV) | -402 | EVE-derived -> suppress |
| MLV Moloney genome (exogenous)      | bit 581, 47.5% id | bit 3466, 99.9% id (MLV pol) | +2885 | genuine -> keep |

The MLV case is the load-bearing one: a real exogenous gammaretrovirus DOES hit the EVE DB
(47.5% id to gammaretroviral HERVs), so "hits EVE -> FP" alone would wrongly suppress it.
The competitive margin spares it. This is why the screen compares EVE vs set-D rather than
thresholding a single EVE hit.

### Held-out recall test (the tests above are self-hits; this one is not)

K113 (NC_022518.1) is a real HERV-K113 provirus ABSENT from CHM13, so it is genuinely NOT
in the CHM13-derived eve_nt / eve_prot DBs -- an independent probe. Chopped into 18x 500 bp
fragments (2026-07-16):

| arm | K113 frags recovered | favoring EVE over set-D | mean identity |
|---|---|---|---|
| protein (blastx) | 18/18 (100%) | 18/18 | 97.7% aa |
| nucleotide (minimap2 -x sr) | 18/18 (100%) | 18/18 (0 hit set-D) | 94.0% nt |

100% recall of a polymorphic insertion absent from the reference, via paralog pooling
through the related HML-2 loci that ARE in the DB. This is the strong result: the screen is
robust to individual-specific / absent-from-reference HERV insertions.

Negative controls (2026-07-16): 20 random-DNA 500 bp fragments -> 0 hits to eve_prot and 0
to eve_nt; a human non-HERV genomic window (CHM13 NC_060925.1:3000000-3000499, verified 0
overlap with the HERV BED) -> 0 hits to eve_nt. No spurious competitive margin from the
657k-peptide DB at min-aa 50 / e<=1e-5.

## Nucleotide minimap2 indices (minimap2 2.30)

Built on HTCF, co-located with the reference FASTAs and the data at
/scratch/sahlab/shandley/eve_exclusion (that is where the screen runs; not copied local).

- EVE-nt from set_b_hervfull.fa (99,035 loci, 184 Mb -- full reference incl. solo LTRs, not
  just internal): eve_nt.asm20.mmi (537 MB, contig mode), eve_nt.sr.mmi (507 MB, read mode).
  8 s / 1.2 GB RAM each.
- set-D-nt from set_d_nt.fa (14 genomes = the 12-genome panel + HIV-1 NC_001802.1 + HTLV-1
  NC_000858.1 nucleotide, so set-D-nt matches set-D-prot): set_d_nt.asm20.mmi,
  set_d_nt.sr.mmi (~1.5 MB each).

Two presets per reference: -x asm20 for viral-candidate contigs, -x sr for reads
(RNAseq/unassembled). Note the nucleotide indices use the FULL EVE reference (solo LTRs
included) whereas the protein DB uses only internal loci -- correct, because solo LTRs have
nucleotide identity worth matching but no protein content.

### Functional nucleotide test (minimap2 -x sr, 500 bp fragments)

| query fragment | vs EVE-nt | vs set-D-nt |
|---|---|---|
| HERV-K interior | hit (465/485 match, mapq 60) | no hit |
| MLV Moloney interior | no hit | hit (493/493, mapq 60) |

Binary at the nucleotide level: MLV does NOT nucleotide-align to gammaretroviral HERVs
(contrast the 47.5% protein hit above). So the nucleotide arm is the high-specificity DNA
path; the protein arm supplies divergence sensitivity + competitive-margin handling.

## Rebuild

    # set-D proteins (re-fetch from NCBI):
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<accs>&rettype=fasta_cds_aa&retmode=text" -o set_d_prot.faa
    diamond makedb --in set_d_prot.faa -d set_d_prot
    # EVE proteins (from internal-loci nucleotide FASTA):
    python3 ../scripts/translate_orfs.py set_b_internal.fa -o eve_prot.faa --min-aa 50
    diamond makedb --in eve_prot.faa -d eve_prot
