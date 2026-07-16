#!/usr/bin/env python3
"""
Run the ciHHV-6 vs active-HHV-6 discrimination engine on REAL reads.

Same three signals as the simulation slice (breadth, copy-number-vs-host, telomeric
junction), plus one addition forced by real data: observed-vs-expected breadth. Real WGS
of a ciHHV-6 carrier is modest coverage, so a present genome shows low RAW breadth (at
0.5x depth only ~39% of positions are covered). Comparing observed breadth to the Poisson
expectation 1 - exp(-depth) tells "uniformly present but low coverage" apart from
"patchy cross-mapping", which a fixed breadth cutoff cannot.

Maps reads to a combined reference (HHV-6 + a human chunk from T2T) with minimap2, so host
depth is the actual human coverage in the sample (works for WGS and for a pure-virus run,
where host depth ~0 -> copy-number -> infinity -> active).

Note on ciHHV-6 copy-number: an integrated genome is 1 copy in a diploid host, so viral
depth ~= 0.5x host depth (CN ~ 0.5), well below an active infection's CN >> 1.
"""

import argparse
import math
import os
import subprocess
import sys
import urllib.request

TELO, TELO_RC = "TTAGGG", "CCCTAA"


def fetch_ncbi(acc):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=fasta&retmode=text")
    with urllib.request.urlopen(url, timeout=90) as r:
        data = r.read().decode()
    return "".join(l.strip() for l in data.splitlines()
                   if not l.startswith(">")).upper()


def analyze(sam_path, hhv6_name, host_name, hhv6_len, host_len):
    import pysam
    cov = [0] * hhv6_len
    host_bases = 0
    junc = 0
    hhv6_reads = 0
    REF = {0, 2, 7, 8}
    with pysam.AlignmentFile(sam_path, "r") as sam:
        for a in sam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            span = sum(l for op, l in a.cigartuples if op in REF)
            if a.reference_name == hhv6_name:
                hhv6_reads += 1
                st = a.reference_start
                for p in range(st, min(st + span, hhv6_len)):
                    cov[p] += 1
                ct = a.cigartuples
                q = a.query_sequence or ""
                if ct and ct[0][0] == 4 and ct[0][1] >= 20:
                    c = q[:ct[0][1]]
                    if c.count(TELO) + c.count(TELO_RC) >= 3:
                        junc += 1
                elif ct and ct[-1][0] == 4 and ct[-1][1] >= 20:
                    c = q[-ct[-1][1]:]
                    if c.count(TELO) + c.count(TELO_RC) >= 3:
                        junc += 1
            elif a.reference_name == host_name:
                host_bases += span
    covered = sum(1 for c in cov if c > 0)
    mean_v = sum(cov) / hhv6_len
    breadth = covered / hhv6_len
    exp_breadth = 1 - math.exp(-mean_v) if mean_v > 0 else 0.0
    mean_h = host_bases / host_len if host_len else 0.0
    cn = mean_v / mean_h if mean_h > 0 else float("inf")
    return dict(hhv6_reads=hhv6_reads, breadth=breadth, exp_breadth=exp_breadth,
                mean_v=mean_v, mean_h=mean_h, cn=cn, junc=junc)


def call(m):
    if m["hhv6_reads"] < 30 or m["mean_v"] < 0.02:
        return "HHV-6 not detected (below coverage floor)"
    even = m["exp_breadth"] > 0 and m["breadth"] >= 0.6 * m["exp_breadth"]
    if not even:
        return "PATCHY coverage (cross-map / partial), not a whole genome present"
    cn = m["cn"]
    if cn <= 1.5 and m["junc"] >= 1:
        return "ENDOGENOUS (ciHHV-6): ~1 copy/cell (CN tracks host) + telomeric junction"
    if cn <= 1.5:
        return ("ENDOGENOUS-consistent: ~1 copy/cell (CN tracks host); "
                "junction not observed (coverage-limited)")
    if cn >= 3:
        return "ACTIVE INFECTION: high copy relative to host, no integration signature"
    return "AMBIGUOUS: needs orthogonal data"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reads", required=True, help="comma-separated FASTQ (1 or 2)")
    p.add_argument("--label", default="sample")
    p.add_argument("--hhv6-acc", default="NC_001664.4", help="HHV-6A 6A; NC_000898.1 6B")
    p.add_argument("--t2t", required=True)
    p.add_argument("--host-chrom", default="NC_060925.1")
    p.add_argument("--host-start", type=int, default=100_000_000)
    p.add_argument("--host-len", type=int, default=20_000_000)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--workdir", default=".")
    args = p.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    import pysam
    ref = pysam.FastaFile(args.t2t)
    host = ref.fetch(args.host_chrom, args.host_start,
                     args.host_start + args.host_len).upper()
    ref.close()
    hhv6 = fetch_ncbi(args.hhv6_acc)
    host_name, hhv6_name = "host_chunk", args.hhv6_acc

    comb = os.path.join(args.workdir, f"ref_{args.label}.fa")
    with open(comb, "w") as f:
        for nm, sq in [(host_name, host), (hhv6_name, hhv6)]:
            f.write(f">{nm}\n")
            for i in range(0, len(sq), 80):
                f.write(sq[i:i + 80] + "\n")

    sam = os.path.join(args.workdir, f"aln_{args.label}.sam")
    reads = args.reads.split(",")
    print(f"[{args.label}] mapping {len(reads)} FASTQ to HHV-6 ({args.hhv6_acc}) "
          f"+ {args.host_len//1_000_000}Mb host...", flush=True)
    with open(sam, "w") as out:
        subprocess.run(["minimap2", "-ax", "sr", "-t", str(args.threads), comb] + reads,
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    m = analyze(sam, hhv6_name, host_name, len(hhv6), len(host))
    verdict = call(m)

    print(f"\n=== {args.label}: ciHHV-6 vs active HHV-6 on REAL reads ===")
    print(f"HHV-6 reference : {args.hhv6_acc} ({len(hhv6):,} bp)")
    print(f"HHV-6 reads     : {m['hhv6_reads']:,}")
    print(f"HHV-6 depth     : {m['mean_v']:.3f}x")
    print(f"host depth      : {m['mean_h']:.3f}x")
    print(f"copy-number     : {m['cn']:.3f}  (viral/host; ciHHV-6~0.5, active>>1)")
    print(f"breadth obs     : {m['breadth']*100:.1f}%")
    print(f"breadth expected: {m['exp_breadth']*100:.1f}%  (Poisson at this depth)")
    print(f"junction reads  : {m['junc']}")
    print(f"CALL            : {verdict}")


if __name__ == "__main__":
    main()
