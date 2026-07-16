#!/usr/bin/env python3
"""
Stage 2 discrimination engine, vertical slice: ciHHV-6 vs active HHV-6.

The hardest endo/exo case. Both an inherited chromosomally-integrated HHV-6 (ciHHV-6) and
an active HHV-6 infection present the whole ~160 kb genome, so coverage BREADTH cannot
tell them apart. Two other signals can:

  copy-number vs host : ciHHV-6 sits at ~1 viral genome per cell, so HHV-6 depth tracks
                        host depth (ratio ~1); active infection runs much higher.
  host-virus junction : ciHHV-6 integrates at a telomere, so some reads split between
                        HHV-6 and a human telomeric-repeat (TTAGGG)n tail; an active,
                        episomal infection has no such junction.

This proves the engine computes those signals and calls the two scenarios correctly on
controlled simulated data, before it is run on real SRA. Simulation-first: prove the
logic, then validate on real data.

References: HHV-6B anchor NC_000898.1 (fetched), plus a human background chunk from T2T.
Alignment: minimap2 short-read preset (shell-out; defers the FFI decision). Analysis is
pure pysam over the SAM, so no samtools is needed.
"""

import argparse
import random
import subprocess
import sys
import urllib.request

TELO = "TTAGGG"
TELO_RC = "CCCTAA"


def fetch_ncbi(acc: str) -> str:
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=fasta&retmode=text")
    with urllib.request.urlopen(url, timeout=90) as r:
        data = r.read().decode()
    return "".join(l.strip() for l in data.splitlines()
                   if not l.startswith(">")).upper()


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(r: str, err: float, rng) -> str:
    if err <= 0:
        return r
    out = list(r)
    for i in range(len(out)):
        if rng.random() < err:
            out[i] = rng.choice("ACGT")
    return "".join(out)


def sim_from(seq: str, depth: float, rlen: int, rng, err, tag):
    """Uniform shotgun reads from a sequence at the given depth."""
    n = int(len(seq) * depth / rlen)
    reads = []
    span = max(1, len(seq) - rlen)
    for i in range(n):
        p = rng.randint(0, span)
        r = seq[p:p + rlen]
        if len(r) < rlen:
            continue
        if rng.random() < 0.5:
            r = rc(r)
        reads.append((f"{tag}{i}", mutate(r, err, rng)))
    return reads


def junction_reads(hhv6: str, n: int, rlen: int, rng):
    """Reads spanning an HHV-6 subterminal region into a human telomere tail. The
    HHV-6 part maps; the telomere tail soft-clips."""
    reads = []
    anchor = len(hhv6) - 2000  # near the 3' end but in mappable sequence
    hlen = 80
    tlen = rlen - hlen
    telo = (TELO * (tlen // len(TELO) + 1))[:tlen]
    for i in range(n):
        p = max(0, anchor + rng.randint(-200, 200))
        reads.append((f"junc{i}", hhv6[p:p + hlen] + telo))
    return reads


def write_fastq(reads, path):
    with open(path, "w") as f:
        for name, seq in reads:
            f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def write_fasta(records, path):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def is_telomeric(clip: str) -> bool:
    return clip.count(TELO) + clip.count(TELO_RC) >= 3


def analyze(sam_path, hhv6_name, host_name, hhv6_len, host_len):
    import pysam
    cov = bytearray()  # not used; use int array
    cov_hhv6 = [0] * hhv6_len
    host_aligned_bases = 0
    junctions = 0
    REF_CONSUME = {0, 2, 7, 8}  # M, D, =, X
    with pysam.AlignmentFile(sam_path, "r") as sam:
        for a in sam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            ref = a.reference_name
            span = sum(l for op, l in a.cigartuples if op in REF_CONSUME)
            if ref == hhv6_name:
                start = a.reference_start
                for p in range(start, min(start + span, hhv6_len)):
                    cov_hhv6[p] += 1
                # soft-clip telomere check (op 4 = S)
                ct = a.cigartuples
                q = a.query_sequence or ""
                if ct and ct[0][0] == 4 and ct[0][1] >= 20:
                    if is_telomeric(q[:ct[0][1]]):
                        junctions += 1
                elif ct and ct[-1][0] == 4 and ct[-1][1] >= 20:
                    if is_telomeric(q[-ct[-1][1]:]):
                        junctions += 1
            elif ref == host_name:
                host_aligned_bases += span
    covered = sum(1 for c in cov_hhv6 if c > 0)
    mean_hhv6 = sum(cov_hhv6) / hhv6_len
    breadth = covered / hhv6_len
    mean_host = host_aligned_bases / host_len if host_len else 0.0
    cn = mean_hhv6 / mean_host if mean_host > 0 else float("inf")
    return breadth, mean_hhv6, mean_host, cn, junctions


def call(breadth, cn, junctions):
    if breadth < 0.5:
        return "not_present / patchy (no whole-genome HHV-6)"
    if junctions >= 5 and cn <= 3:
        return "ENDOGENOUS (ciHHV-6): whole genome, ~1 copy/cell, telomeric junction"
    if cn >= 5 and junctions < 5:
        return "ACTIVE INFECTION: whole genome, high copy, no junction"
    return "AMBIGUOUS: needs orthogonal data"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--t2t", required=True)
    p.add_argument("--host-chrom", default="NC_060925.1")
    p.add_argument("--host-start", type=int, default=100_000_000)
    p.add_argument("--host-len", type=int, default=2_000_000)
    p.add_argument("--hhv6-acc", default="NC_000898.1")
    p.add_argument("--rlen", type=int, default=150)
    p.add_argument("--err", type=float, default=0.001)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--workdir", default=".")
    args = p.parse_args()

    import os
    os.makedirs(args.workdir, exist_ok=True)
    rng = random.Random(42)

    import pysam
    print("Loading references...", flush=True)
    ref = pysam.FastaFile(args.t2t)
    host = ref.fetch(args.host_chrom, args.host_start,
                     args.host_start + args.host_len).upper()
    ref.close()
    hhv6 = fetch_ncbi(args.hhv6_acc)
    host_name, hhv6_name = "host_chunk", args.hhv6_acc
    print(f"  host chunk {len(host):,} bp; HHV-6B {len(hhv6):,} bp", flush=True)

    comb = os.path.join(args.workdir, "combined_ref.fa")
    write_fasta([(host_name, host), (hhv6_name, hhv6)], comb)

    scenarios = {
        "ciHHV-6":  dict(host_d=10, hhv6_d=10,  junc=40),
        "active":   dict(host_d=10, hhv6_d=200, junc=0),
        "negative": dict(host_d=10, hhv6_d=0,   junc=0),
    }

    rows = []
    for name, s in scenarios.items():
        reads = sim_from(host, s["host_d"], args.rlen, rng, args.err, f"h_{name}_")
        if s["hhv6_d"] > 0:
            reads += sim_from(hhv6, s["hhv6_d"], args.rlen, rng, args.err, f"v_{name}_")
        if s["junc"] > 0:
            reads += junction_reads(hhv6, s["junc"], args.rlen, rng)
        rng.shuffle(reads)
        fq = os.path.join(args.workdir, f"reads_{name}.fq")
        sam = os.path.join(args.workdir, f"aln_{name}.sam")
        write_fastq(reads, fq)
        print(f"[{name}] {len(reads):,} reads -> minimap2", flush=True)
        with open(sam, "w") as out:
            subprocess.run(["minimap2", "-ax", "sr", "-t", str(args.threads),
                            comb, fq], stdout=out, stderr=subprocess.DEVNULL, check=True)
        breadth, md_v, md_h, cn, junc = analyze(
            sam, hhv6_name, host_name, len(hhv6), len(host))
        verdict = call(breadth, cn, junc)
        rows.append((name, breadth, md_v, md_h, cn, junc, verdict))

    print("\n=== ciHHV-6 vs active HHV-6: Stage 2 evidence ===")
    print(f"{'scenario':<10} {'breadth':>8} {'vDepth':>8} {'hDepth':>8} "
          f"{'CN(v/h)':>8} {'junc':>5}  call")
    for name, breadth, md_v, md_h, cn, junc, verdict in rows:
        cn_s = f"{cn:.2f}" if cn != float("inf") else "inf"
        print(f"{name:<10} {breadth*100:>7.1f}% {md_v:>8.1f} {md_h:>8.1f} "
              f"{cn_s:>8} {junc:>5}  {verdict}")
    print("\nNote: breadth is ~100% for BOTH ciHHV-6 and active (whole genome present),")
    print("so breadth cannot separate them; copy-number and junction do.")


if __name__ == "__main__":
    main()
