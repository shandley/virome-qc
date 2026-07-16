#!/usr/bin/env python3
"""6-frame translate a nucleotide FASTA into peptide ORFs for a DIAMOND blastx target.

Splits each of the 6 frames on stop codons and keeps peptide runs >= --min-aa.
Used to build the EVE-protein DIAMOND database from endogenous HERV loci, where the
loci are degraded (stops, frameshifts) so annotated CDS do not exist. Standard genetic
code; ambiguous codons translate to X.
"""
import argparse
import sys

CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(seq: str) -> str:
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(CODON.get(seq[i:i + 3], "X"))
    return "".join(out)


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--min-aa", type=int, default=50)
    args = ap.parse_args()

    n_loci = n_pep = 0
    with open(args.out, "w") as out:
        for name, seq in read_fasta(args.fasta):
            seq = seq.upper()
            n_loci += 1
            rc = revcomp(seq)
            for strand, s in (("+", seq), ("-", rc)):
                for frame in range(3):
                    aa = translate(s[frame:])
                    for j, pep in enumerate(aa.split("*")):
                        if len(pep) >= args.min_aa:
                            n_pep += 1
                            out.write(f">{name}|{strand}{frame}|{j}\n{pep}\n")
            if n_loci % 5000 == 0:
                print(f"  {n_loci} loci, {n_pep} peptides", file=sys.stderr, flush=True)
    print(f"done: {n_loci} loci -> {n_pep} peptides (>= {args.min_aa} aa)", file=sys.stderr)


if __name__ == "__main__":
    main()
