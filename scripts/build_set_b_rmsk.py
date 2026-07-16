#!/usr/bin/env python3
"""
Coordinate-based EVE exclusion set (set B) from the CHM13 RepeatMasker annotation.

The Dfam-API version (build_set_b_hervk.py) works but does not scale past a few
families: it takes one family per request and enforces a ~1 Mb region limit, so broad
coverage (HERV-H alone has ~5,800 internal copies plus thousands of LTR7) would need
tens of thousands of calls. This version reads the bulk RepeatMasker .out for T2T-CHM13
once and extracts every targeted family in a single pass.

Source: hs1.repeatMasker.out.gz from UCSC
    https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz
    (RepeatMasker -s sensitive; uses UCSC chrom names chr1..chrY, 1-based coordinates.)

Families are matched by an EXACT repName allowlist, not a prefix: prefix-matching "LTR7"
would wrongly capture LTR78/LTR79, which are distinct ERV1 families, not HERV-H. Copy
counts below are from the CHM13 .out and are here so the allowlist is auditable.

Internal + LTR hits of the same HERV group that fall within --merge-distance are merged
into one locus (a full provirus), and each locus is labeled provirus / solo_LTR /
internal_only. Groups are merged independently so a HERV-H LTR never joins an HML-2
internal.

Requirements: pysam, the T2T FASTA plus its .fai, and hs1.repeatMasker.out.gz.
"""

import argparse
import gzip
import sys

# UCSC chromosome name -> RefSeq accession, T2T-CHM13v2.0 (NCBI assembly report
# GCF_009914755.1). chrM has no target EVEs and no RefSeq mapping here; it is skipped.
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

# HERV group -> (internal repNames, LTR repNames, expected class/family).
# Exact repNames, verified present in the CHM13 RepeatMasker .out.
TARGET_GROUPS = {
    "HML-2": {
        "internal": {"HERVK-int"},
        "ltr": {"LTR5", "LTR5A", "LTR5B", "LTR5_Hs"},
        "class": "LTR/ERVK",
    },
    "HERV-H": {
        "internal": {"HERVH-int", "HERVH48-int"},
        "ltr": {"LTR7", "LTR7B", "LTR7C", "LTR7Y"},
        "class": "LTR/ERV1",
    },
    "HERV-W": {
        "internal": {"HERV17-int", "HERV17B-int"},
        "ltr": {"LTR17", "LTR17b"},
        "class": "LTR/ERV1",
    },
}


def build_lookup(groups):
    """repName -> (group, is_internal)."""
    lut = {}
    for group, spec in groups.items():
        for name in spec["internal"]:
            lut[name] = (group, True)
        for name in spec["ltr"]:
            lut[name] = (group, False)
    return lut


def read_fai_lengths(fai_path: str) -> dict:
    lengths = {}
    with open(fai_path) as f:
        for line in f:
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)
    return lengths


def parse_rmsk(rmsk_path: str, lut: dict, class_set=None):
    """Yield (chrom, start, end, strand, repName, group, is_internal) for target
    families. RepeatMasker .out: 3 header lines, then whitespace-delimited columns
    5=query, 6=start, 7=end, 9=strand(+/C), 10=repName, 11=class/family. 1-based."""
    opener = gzip.open if rmsk_path.endswith(".gz") else open
    kept = 0
    class_mismatch = 0
    with opener(rmsk_path, "rt") as f:
        for i, line in enumerate(f):
            if i < 3 or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 11:
                continue
            repName = fields[9]
            cls = fields[10]
            if class_set is not None:
                # Class-based (full HERV): cls like "LTR/ERVK", "LTR/ERVL-MaLR".
                # Match the full post-slash token minus a trailing "?", so ERVL is
                # kept but ERVL-MaLR is excluded unless explicitly listed.
                if not cls.startswith("LTR/"):
                    continue
                token = cls.split("/", 1)[1].rstrip("?")
                if token not in class_set:
                    continue
                group = token
                is_internal = "-int" in repName
            else:
                hit = lut.get(repName)
                if hit is None:
                    continue
                group, is_internal = hit
                if cls != TARGET_GROUPS[group]["class"]:
                    class_mismatch += 1
                    continue
            chrom = fields[4]
            start, end = int(fields[5]), int(fields[6])
            strand = "-" if fields[8] == "C" else "+"
            kept += 1
            yield chrom, start, end, strand, repName, group, is_internal
    if class_mismatch:
        print(f"note: {class_mismatch} rows skipped on class mismatch", file=sys.stderr)
    print(f"RepeatMasker rows kept: {kept}")


def merge_group(intervals: list, merge_distance: int):
    """Merge intervals (same chrom+group) into loci, tracking families, strands,
    and whether internal/LTR contributed. Input: (start,end,strand,repName,is_internal)."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    s, e, strand, name, is_int = intervals[0]
    cur = {"start": s, "end": e, "strands": {strand}, "families": {name},
           "has_int": is_int, "has_ltr": not is_int, "n_hits": 1}
    for s, e, strand, name, is_int in intervals[1:]:
        if s <= cur["end"] + merge_distance:
            cur["end"] = max(cur["end"], e)
            cur["strands"].add(strand)
            cur["families"].add(name)
            cur["has_int"] |= is_int
            cur["has_ltr"] |= not is_int
            cur["n_hits"] += 1
        else:
            merged.append(cur)
            cur = {"start": s, "end": e, "strands": {strand}, "families": {name},
                   "has_int": is_int, "has_ltr": not is_int, "n_hits": 1}
    merged.append(cur)
    return merged


def structure(locus) -> str:
    if locus["has_int"] and locus["has_ltr"]:
        return "provirus"
    return "internal_only" if locus["has_int"] else "solo_LTR"


def main():
    p = argparse.ArgumentParser(description="RepeatMasker-based HERV set B extraction")
    p.add_argument("--rmsk", required=True, help="hs1.repeatMasker.out(.gz)")
    p.add_argument("--t2t", required=True, help="T2T-CHM13 FASTA (needs .fai)")
    p.add_argument("--herv-classes", default=None,
                   help="Class-based full-HERV mode: comma-separated RepeatMasker ERV "
                        "classes, e.g. 'ERVK,ERV1,ERVL' (excludes ERVL-MaLR unless listed). "
                        "Overrides --groups.")
    p.add_argument("--groups", default="HML-2,HERV-H,HERV-W",
                   help="Comma-separated HERV groups to include")
    p.add_argument("--merge-distance", type=int, default=2000)
    p.add_argument("--min-region-len", type=int, default=500)
    p.add_argument("--output-bed", default="set_b_rmsk.bed")
    p.add_argument("--output-fasta", default="set_b_rmsk.fa")
    p.add_argument("--output-rust", default=None)
    args = p.parse_args()

    class_set = None
    if args.herv_classes:
        class_set = {c.strip() for c in args.herv_classes.split(",") if c.strip()}
        groups = {c: {"class": f"LTR/{c}"} for c in sorted(class_set)}
        lut = {}
        print(f"Full-HERV class mode: {sorted(class_set)}", flush=True)
    else:
        wanted = [g.strip() for g in args.groups.split(",")]
        groups = {g: TARGET_GROUPS[g] for g in wanted if g in TARGET_GROUPS}
        lut = build_lookup(groups)
    print(f"Groups: {', '.join(groups)}", flush=True)

    # Collect target hits, bucketed by (chrom, group).
    buckets = {}
    for chrom, start, end, strand, name, group, is_int in parse_rmsk(
            args.rmsk, lut, class_set):
        if chrom not in CHR_TO_NC:
            continue
        buckets.setdefault((chrom, group), []).append(
            (start, end, strand, name, is_int))

    # Merge each bucket, then flatten to per-chrom loci.
    loci_by_chrom = {}
    for (chrom, group), intervals in buckets.items():
        for locus in merge_group(intervals, args.merge_distance):
            locus["group"] = group
            loci_by_chrom.setdefault(chrom, []).append(locus)

    try:
        import pysam
    except ImportError:
        print("pysam required. pip install pysam", file=sys.stderr)
        sys.exit(1)
    nc_lengths = read_fai_lengths(args.t2t + ".fai")
    ref = pysam.FastaFile(args.t2t)

    records, bed_rows = [], []
    from collections import defaultdict
    stats = defaultdict(lambda: {"provirus": 0, "solo_LTR": 0, "internal_only": 0,
                                 "bp": 0})
    idx = 0
    for chrom in sorted(loci_by_chrom, key=lambda c: list(CHR_TO_NC).index(c)):
        nc = CHR_TO_NC[chrom]
        for locus in sorted(loci_by_chrom[chrom], key=lambda x: x["start"]):
            if locus["end"] - locus["start"] + 1 < args.min_region_len:
                continue
            start0 = max(0, locus["start"] - 1)
            end0 = min(nc_lengths[nc], locus["end"])
            seq = ref.fetch(nc, start0, end0).upper()
            struct = structure(locus)
            grp = locus["group"]
            stats[grp][struct] += 1
            stats[grp]["bp"] += len(seq)
            name = f"setb_{idx:05d}_{grp}_{nc}_{start0}_{end0}_{struct}"
            fams = ",".join(sorted(locus["families"]))
            records.append((name, seq))
            bed_rows.append(
                f"{nc}\t{start0}\t{end0}\t{name}\t{chrom}\t{grp}\t{struct}\t{fams}\t"
                f"{''.join(sorted(locus['strands']))}\t{locus['n_hits']}\t{len(seq)}")
            idx += 1
    ref.close()

    with open(args.output_bed, "w") as f:
        f.write("#nc_accession\tstart\tend\tname\tucsc_chrom\tgroup\tstructure\t"
                "families\tstrand\tn_hits\tlength\n")
        f.write("\n".join(bed_rows) + "\n")
    total_bp = sum(len(s) for _, s in records)
    with open(args.output_fasta, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
    if args.output_rust:
        with open(args.output_rust, "w") as f:
            f.write("//! HERV EVE exclusion set from T2T-CHM13 RepeatMasker\n//!\n")
            f.write(f"//! Generated by scripts/build_set_b_rmsk.py. Groups: "
                    f"{', '.join(groups)}. {len(records)} loci, {total_bp:,} bp.\n\n")
            f.write("pub const EVE_SEQUENCES: &[(&str, &[u8])] = &[\n")
            for name, seq in records:
                f.write(f'    ("{name}",\n     b"')
                for i in range(0, len(seq), 80):
                    if i > 0:
                        f.write("\\\n       ")
                    f.write(seq[i:i + 80])
                f.write('"),\n')
            f.write("];\n")

    print(f"\nBED:   {args.output_bed} ({len(records)} loci)")
    print(f"FASTA: {args.output_fasta} ({len(records)} loci, {total_bp:,} bp)")
    print("\n=== Summary by group ===")
    for grp in groups:
        s = stats[grp]
        n = s["provirus"] + s["solo_LTR"] + s["internal_only"]
        print(f"{grp}: {n} loci, {s['bp']:,} bp "
              f"(provirus {s['provirus']}, solo_LTR {s['solo_LTR']}, "
              f"internal_only {s['internal_only']})")
    print(f"TOTAL: {len(records)} loci, {total_bp:,} bp, "
          f"~{total_bp * 2:,} k=31 k-mers (both strands)")


if __name__ == "__main__":
    main()
