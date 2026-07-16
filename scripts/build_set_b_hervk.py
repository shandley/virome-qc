#!/usr/bin/env python3
"""
Prototype: coordinate-based EVE exclusion set (set B) for HERV-K (HML-2).

Instead of aligning a viral reference to T2T-CHM13 (which fails for polymorphic
proviruses like HERV-K113 that are absent from the assembly, and needs a very low
identity threshold to catch fixed paralogs), this pulls the coordinates of every
HERV-K/HML-2 element annotated in CHM13 directly from the Dfam annotation API and
extracts those genomic sequences from the T2T FASTA.

Families (HML-2, class LTR/ERVK), verified against Dfam:
    DF000000188  HERVK    (internal model)
    DF000000540  LTR5     (LTR)
    DF000000556  LTR5A    (LTR)
    DF000000557  LTR5B    (LTR)
    DF000000558  LTR5_Hs  (LTR, youngest)

Internal + LTR hits that merge within --merge-distance are joined into one locus, so a
full provirus (5' LTR + internal + 3' LTR) becomes a single record; a recombined solo LTR
stays its own record. Each locus is labeled by structure (provirus / solo_LTR /
internal_only) from the families that contributed to it.

Data source: Dfam annotations API (https://dfam.org/api/annotations). The API enforces a
per-request region-size limit (~1 Mb) and takes one family per call, so each chromosome
is tiled and each family is queried per window; the window halves on an over-size 400.

Coordinate note:
    Dfam returns UCSC chromosome names (chr1..chrY), 1-based, with seq_start > seq_end on
    the minus strand. The T2T FASTA uses RefSeq accessions (NC_060925.1..). This script
    translates chr -> RefSeq via the mapping below (NCBI assembly report for
    GCF_009914755.1) and extracts with pysam (0-based half-open).

Requirements: pysam, network access to Dfam, the T2T FASTA plus its .fai.

Usage:
    python build_set_b_hervk.py --t2t t2t_chm13.fna --chroms chr19 \
        --output-bed set_b_hervk.bed --output-fasta set_b_hervk.fa
    # genome-wide: --chroms all
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

# UCSC chromosome name -> RefSeq accession, T2T-CHM13v2.0.
# Verified against the NCBI assembly report for GCF_009914755.1.
CHR_TO_NC = {
    "chr1": "NC_060925.1", "chr2": "NC_060926.1", "chr3": "NC_060927.1",
    "chr4": "NC_060928.1", "chr5": "NC_060929.1", "chr6": "NC_060930.1",
    "chr7": "NC_060931.1", "chr8": "NC_060932.1", "chr9": "NC_060933.1",
    "chr10": "NC_060934.1", "chr11": "NC_060935.1", "chr12": "NC_060936.1",
    "chr13": "NC_060937.1", "chr14": "NC_060938.1", "chr15": "NC_060939.1",
    "chr16": "NC_060940.1", "chr17": "NC_060941.1", "chr18": "NC_060942.1",
    "chr19": "NC_060943.1", "chr20": "NC_060944.1", "chr21": "NC_060945.1",
    "chr22": "NC_060946.1", "chrX": "NC_060947.1", "chrY": "NC_060948.1",
}

DFAM_API = "https://dfam.org/api/annotations"

# HML-2 families (verified accessions). Internal first, then the LTR flanks.
FAMILY_NAMES = {
    "DF000000188": "HERVK-int",
    "DF000000540": "LTR5",
    "DF000000556": "LTR5A",
    "DF000000557": "LTR5B",
    "DF000000558": "LTR5_Hs",
}
INTERNAL_FAMILY = "DF000000188"
DEFAULT_FAMILIES = ",".join(FAMILY_NAMES)


def read_fai_lengths(fai_path: str) -> dict:
    lengths = {}
    with open(fai_path) as f:
        for line in f:
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)
    return lengths


def dfam_annotations(chrom: str, start: int, end: int, family: str,
                     timeout: int = 120, retries: int = 4) -> list:
    """Query Dfam for one window and one family. Raises HTTPError(400) so the
    caller can split the window; retries transient errors with backoff."""
    params = {"assembly": "hs1", "chrom": chrom, "start": start, "end": end,
              "family": family, "nrph": "true"}
    url = DFAM_API + "?" + urllib.parse.urlencode(params)
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return json.loads(resp.read().decode()).get("hits") or []
        except urllib.error.HTTPError as e:
            if e.code == 400:
                raise
            last_err = e
            time.sleep(1.0 * (attempt + 1))
        except Exception as e:  # noqa: BLE001 - transient network, retry
            last_err = e
            time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(f"Dfam query failed {chrom}:{start}-{end} {family}: {last_err}")


def collect_windows(chrom: str, chrom_len: int, families: list, window: int,
                    delay: float) -> list:
    """Tile a chromosome and collect hits for every family. Halves the window on
    an over-size 400. Returns (start, end, strand, bit_score, family), 1-based."""
    intervals = []
    pos = 0
    while pos < chrom_len:
        w = window
        end = min(pos + w, chrom_len)
        for family in families:
            placed = False
            while not placed:
                end = min(pos + w, chrom_len)
                try:
                    hits = dfam_annotations(chrom, pos, end, family)
                    placed = True
                except urllib.error.HTTPError as e:
                    if e.code == 400 and w > 50_000:
                        w //= 2
                        continue
                    raise
                time.sleep(delay)
            for h in hits:
                s, e2 = h["seq_start"], h["seq_end"]
                lo, hi = (s, e2) if s <= e2 else (e2, s)
                intervals.append((lo, hi, h.get("strand", "?"),
                                  float(h.get("bit_score", 0.0)), family))
        pos = end
    return intervals


def merge_intervals(intervals: list, merge_distance: int):
    """Merge overlapping/nearby intervals into loci, tracking contributing
    families so a locus can be labeled provirus / solo_LTR / internal_only."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    lo, hi, strand, bits, fam = intervals[0]
    cur = {"start": lo, "end": hi, "strands": {strand}, "best_bits": bits,
           "families": {fam}, "n_hits": 1}
    for lo, hi, strand, bits, fam in intervals[1:]:
        if lo <= cur["end"] + merge_distance:
            cur["end"] = max(cur["end"], hi)
            cur["strands"].add(strand)
            cur["best_bits"] = max(cur["best_bits"], bits)
            cur["families"].add(fam)
            cur["n_hits"] += 1
        else:
            merged.append(cur)
            cur = {"start": lo, "end": hi, "strands": {strand}, "best_bits": bits,
                   "families": {fam}, "n_hits": 1}
    merged.append(cur)
    return merged


def structure_label(families: set) -> str:
    has_int = INTERNAL_FAMILY in families
    has_ltr = any(f != INTERNAL_FAMILY for f in families)
    if has_int and has_ltr:
        return "provirus"
    if has_int:
        return "internal_only"
    return "solo_LTR"


def extract_and_write(loci_by_chrom: dict, nc_lengths: dict, t2t_ref: str,
                      out_bed: str, out_fasta: str, out_rust: str,
                      min_region_len: int):
    try:
        import pysam
    except ImportError:
        print("pysam required. Install: pip install pysam", file=sys.stderr)
        sys.exit(1)

    ref = pysam.FastaFile(t2t_ref)
    records, bed_rows = [], []
    counts = {"provirus": 0, "solo_LTR": 0, "internal_only": 0}
    idx = 0
    for chrom in sorted(loci_by_chrom, key=lambda c: list(CHR_TO_NC).index(c)):
        nc = CHR_TO_NC[chrom]
        for locus in loci_by_chrom[chrom]:
            if locus["end"] - locus["start"] + 1 < min_region_len:
                continue
            start0 = max(0, locus["start"] - 1)
            end0 = min(nc_lengths.get(nc, locus["end"]), locus["end"])
            seq = ref.fetch(nc, start0, end0)
            struct = structure_label(locus["families"])
            counts[struct] += 1
            fam_names = ",".join(sorted(FAMILY_NAMES.get(f, f)
                                        for f in locus["families"]))
            name = f"setb_{idx:04d}_{nc}_{start0}_{end0}_{struct}"
            records.append((name, seq.upper()))
            bed_rows.append(
                f"{nc}\t{start0}\t{end0}\t{name}\t{chrom}\t"
                f"{''.join(sorted(locus['strands']))}\t{struct}\t{fam_names}\t"
                f"{locus['n_hits']}\t{locus['best_bits']:.1f}\t{len(seq)}"
            )
            idx += 1
    ref.close()

    with open(out_bed, "w") as f:
        f.write("#nc_accession\tstart\tend\tname\tucsc_chrom\tstrand\t"
                "structure\tfamilies\tn_hits\tbest_bit_score\tlength\n")
        f.write("\n".join(bed_rows) + "\n")

    total_bp = sum(len(s) for _, s in records)
    with open(out_fasta, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    print(f"BED:   {out_bed} ({len(records)} loci)")
    print(f"FASTA: {out_fasta} ({len(records)} loci, {total_bp:,} bp)")
    if out_rust:
        _write_rust(records, out_rust, total_bp)
    return len(records), total_bp, counts


def _write_rust(records: list, out_rust: str, total_bp: int):
    with open(out_rust, "w") as f:
        f.write("//! HERV-K (HML-2) EVE exclusion set from T2T-CHM13, by Dfam "
                "coordinates\n//!\n")
        f.write("//! Generated by scripts/build_set_b_hervk.py (internal + LTR5 "
                "families).\n")
        f.write(f"//! {len(records)} loci, {total_bp:,} bp total.\n\n")
        f.write("pub const EVE_SEQUENCES: &[(&str, &[u8])] = &[\n")
        for name, seq in records:
            f.write(f'    ("{name}",\n     b"')
            for i in range(0, len(seq), 80):
                if i > 0:
                    f.write("\\\n       ")
                f.write(seq[i:i + 80])
            f.write('"),\n')
        f.write("];\n")
    print(f"Rust:  {out_rust} ({len(records)} loci, {total_bp:,} bp)")


def main():
    p = argparse.ArgumentParser(description="Coordinate-based HERV-K set B extraction")
    p.add_argument("--t2t", required=True, help="T2T-CHM13 FASTA (needs .fai)")
    p.add_argument("--chroms", default="chr19",
                   help="Comma-separated UCSC chrom names, or 'all'")
    p.add_argument("--families", default=DEFAULT_FAMILIES,
                   help="Comma-separated Dfam family accessions "
                        "(default: HERVK internal + LTR5/LTR5A/LTR5B/LTR5_Hs)")
    p.add_argument("--window", type=int, default=1_000_000)
    p.add_argument("--merge-distance", type=int, default=2000)
    p.add_argument("--min-region-len", type=int, default=500)
    p.add_argument("--delay", type=float, default=0.25)
    p.add_argument("--output-bed", default="set_b_hervk.bed")
    p.add_argument("--output-fasta", default="set_b_hervk.fa")
    p.add_argument("--output-rust", default=None)
    args = p.parse_args()

    nc_lengths = read_fai_lengths(args.t2t + ".fai")
    families = [f.strip() for f in args.families.split(",") if f.strip()]
    chroms = list(CHR_TO_NC) if args.chroms == "all" else \
        [c.strip() for c in args.chroms.split(",")]

    print(f"Families: {', '.join(FAMILY_NAMES.get(f, f) for f in families)}", flush=True)
    loci_by_chrom, grand_hits = {}, 0
    for chrom in chroms:
        if chrom not in CHR_TO_NC:
            print(f"skip unknown chrom {chrom}", file=sys.stderr)
            continue
        chrom_len = nc_lengths[CHR_TO_NC[chrom]]
        n_windows = -(-chrom_len // args.window)
        print(f"{chrom} ({chrom_len:,} bp, ~{n_windows} windows x {len(families)} fam)...",
              flush=True)
        raw = collect_windows(chrom, chrom_len, families, args.window, args.delay)
        loci = merge_intervals(raw, args.merge_distance)
        loci_by_chrom[chrom] = loci
        grand_hits += len(raw)
        print(f"  {len(raw)} raw hits -> {len(loci)} merged loci", flush=True)

    n_loci, total_bp, counts = extract_and_write(
        loci_by_chrom, nc_lengths, args.t2t,
        args.output_bed, args.output_fasta, args.output_rust, args.min_region_len)

    print("\n=== Summary ===")
    print(f"Chromosomes: {len(loci_by_chrom)}")
    print(f"Raw Dfam hits: {grand_hits}")
    print(f"Merged loci (>= {args.min_region_len} bp): {n_loci}")
    print(f"  provirus (internal+LTR): {counts['provirus']}")
    print(f"  solo LTR: {counts['solo_LTR']}")
    print(f"  internal only: {counts['internal_only']}")
    print(f"Total sequence: {total_bp:,} bp")
    print(f"Estimated k=31 k-mers: ~{total_bp * 2:,} (both strands)")


if __name__ == "__main__":
    main()
