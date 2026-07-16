#!/usr/bin/env python3
"""
EVE screen: quantify endogenous-viral-element (EVE) false positives among viral-candidate
reads, by k-mer containment against the host EVE map (set B, HERV loci from T2T).

Part A (ground truth): simulate host-EVE reads (drawn from the set-B HERV map) and real
exogenous-virus reads (HIV-1 + HTLV-1). A naive retroviral classifier flags BOTH as
"retrovirus" (both are retroviral k-mers), so every EVE read is a false positive. A good
EVE screen flags the EVE reads (removing the FPs) and spares the real viruses (HIV/HTLV
are phylogenetically distant from HERVs, so they should not match set B).

Part B (real gut virome): screen real reads against set B; the fraction matching is the
host-HERV load that a retroviral classifier would mis-call as virus.
"""
import argparse
import gzip
import random
import sys

K = 31
COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def canon(km):
    rc = km.translate(COMP)[::-1]
    return km if km < rc else rc


def load_seqs(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, buf = [], []
    with op(path, "rb") as f:
        for line in f:
            line = line.strip()
            if line.startswith(b">"):
                if buf:
                    seqs.append(b"".join(buf))
                    buf = []
            else:
                buf.append(line.upper())
        if buf:
            seqs.append(b"".join(buf))
    return seqs


def build_kmers(seqs):
    kmers = set()
    for s in seqs:
        for i in range(len(s) - K + 1):
            km = s[i:i + K]
            if b"N" in km:
                continue
            kmers.add(canon(km))
    return kmers


def containment(read, kmers):
    tot = hit = 0
    for i in range(len(read) - K + 1):
        km = read[i:i + K]
        if b"N" in km:
            continue
        tot += 1
        if canon(km) in kmers:
            hit += 1
    return hit / tot if tot else 0.0


def sim_reads(seqs, n, rlen, rng):
    seqs = [s for s in seqs if len(s) >= rlen]
    w = [len(s) for s in seqs]
    out = []
    for _ in range(n):
        s = rng.choices(seqs, weights=w)[0]
        p = rng.randint(0, len(s) - rlen)
        r = s[p:p + rlen]
        if rng.random() < 0.5:
            r = r.translate(COMP)[::-1]
        out.append(r)
    return out


def read_fastq(path, limit):
    op = gzip.open if path.endswith(".gz") else open
    n = 0
    with op(path, "rb") as f:
        while True:
            if not f.readline():
                break
            s = f.readline().strip().upper()
            f.readline()
            f.readline()
            yield s
            n += 1
            if limit and n >= limit:
                break


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--eve-ref", required=True, help="set B HERV map FASTA")
    p.add_argument("--virus-ref", help="real exogenous virus FASTA (HIV/HTLV) for Part A")
    p.add_argument("--reads", help="real reads FASTQ for Part B")
    p.add_argument("--threshold", type=float, default=0.15)
    p.add_argument("--n-sim", type=int, default=5000)
    p.add_argument("--rlen", type=int, default=150)
    p.add_argument("--limit", type=int, default=1_000_000)
    args = p.parse_args()
    rng = random.Random(42)
    th = args.threshold

    print("Building set-B (host EVE map) k-mer set...", flush=True)
    eve_seqs = load_seqs(args.eve_ref)
    eve_k = build_kmers(eve_seqs)
    print(f"  set B: {len(eve_seqs)} loci, {len(eve_k):,} canonical {K}-mers\n", flush=True)

    if args.virus_ref:
        print("=== Part A: ground truth (planted EVE vs real virus) ===", flush=True)
        vseqs = load_seqs(args.virus_ref)
        eve_reads = sim_reads(eve_seqs, args.n_sim, args.rlen, rng)
        vir_reads = sim_reads(vseqs, args.n_sim, args.rlen, rng)
        ec = [containment(r, eve_k) for r in eve_reads]
        vc = [containment(r, eve_k) for r in vir_reads]
        eve_flag = sum(1 for x in ec if x > th)
        vir_flag = sum(1 for x in vc if x > th)
        ne, nv = len(ec), len(vc)
        naive_fp = ne / (ne + nv)
        residual_fp_reads = ne - eve_flag          # EVE that slipped the screen
        kept_virus = nv - vir_flag                 # real virus retained
        fp_after = residual_fp_reads / (residual_fp_reads + kept_virus) \
            if (residual_fp_reads + kept_virus) else 0.0
        print(f"  EVE reads (n={ne}):   mean containment {sum(ec)/ne*100:.1f}%, "
              f"flagged {eve_flag}/{ne} ({eve_flag/ne*100:.1f}%)")
        print(f"  virus reads (n={nv}): mean containment {sum(vc)/nv*100:.1f}%, "
              f"flagged {vir_flag}/{nv} ({vir_flag/nv*100:.1f}%)")
        print(f"  naive retroviral-call FP rate (EVE among calls): {naive_fp*100:.1f}%")
        print(f"  after EVE screen -> FP rate: {fp_after*100:.2f}%  "
              f"(residual EVE FPs {residual_fp_reads}, real viruses lost {vir_flag})")
        print(f"  screen recall (EVE caught): {eve_flag/ne*100:.1f}%; "
              f"specificity (virus spared): {kept_virus/nv*100:.1f}%\n", flush=True)

    if args.reads:
        print("=== Part B: real gut virome EVE load ===", flush=True)
        n = 0
        match = 0
        strong = 0
        summ = 0.0
        for r in read_fastq(args.reads, args.limit):
            c = containment(r, eve_k)
            n += 1
            if c > th:
                match += 1
                summ += c
            if c > 0.5:
                strong += 1
        if n==0:
            print("  no reads (download failed?)"); return
        print(f"  reads screened: {n:,}")
        print(f"  EVE-map matches (>{th:.0%} containment): {match:,} "
              f"({match/n*100:.3f}%)")
        print(f"  strong EVE matches (>50%): {strong:,} ({strong/n*100:.3f}%)")
        if match:
            print(f"  mean containment of matches: {summ/match*100:.1f}%")
        print("  -> these host-HERV reads would be mis-called retroviral by a naive "
              "classifier; the EVE screen flags them.", flush=True)


if __name__ == "__main__":
    main()
