#!/usr/bin/env python3
"""
Build EVE exclusion set for virome-qc host depletion.

Aligns curated high-identity viral references against T2T-CHM13 to identify
regions where viral sequence is integrated into the human genome. Extracts
the host-genome subsequences at these coordinates so that virome-qc can
exclude their k-mers from host containment calculations at runtime.

This prevents exogenous viral reads (e.g., HHV-6B in a patient with active
infection) from being falsely removed as "host" because they share k-mers
with the integrated copy (ciHHV-6 at chr17).

Reference sourcing:
    Curated references are fetched by accession from NCBI (nuccore) so their
    identity is authoritative and recorded, rather than pasted by hand. See
    CURATED_REFERENCES below. Additional references can be supplied with
    --viral-refs and are aligned alongside the fetched ones.

    Important: an extracted region is named by its T2T locus, not by the query
    that matched it. A query only records a "best_match" annotation with the
    alignment identity. This matters because a polymorphic provirus (e.g.
    HERV-K113) is often absent from the CHM13 haplotype, so its query hits a
    paralogous fixed element elsewhere in the genome. Naming the extracted
    locus after the query would mislabel that paralog (this is the bug that
    produced a chr1 HML-2 element mislabeled "HERV-K_K113"). When a curated
    reference's best T2T hit is below POLYMORPHIC_IDENTITY_FLAG, the script
    warns that the reference is likely polymorphic/absent from CHM13 and
    belongs in the competitive/consensus set (set C in REFERENCE_CURATION.md),
    sourced directly from its accession, not represented by a T2T decoy.

Requirements:
    - minimap2 in PATH
    - T2T-CHM13 reference FASTA (3.1 GB)
    - pysam (for FASTA extraction)
    - network access to NCBI E-utilities (for --fetch-curated)

Usage:
    python build_eve_exclusion.py \
        --t2t /path/to/chm13v2.0.fa \
        --output-bed references/t2t_eve_regions.bed \
        --output-fasta references/eve_exclusion_sequences.fa \
        --output-rust src/modules/eve_kmers.rs \
        --output-curated-fasta references/curated_viral_refs.fa \
        --email you@example.org \
        --min-identity 85 \
        --min-length 100

The --output-rust flag generates a Rust source file with embedded byte
sequences (same pattern as erv_sequences.rs) that virome-qc compiles into
the binary. At runtime, these sequences are hashed into an FxHashSet of
k=31 k-mers for the EVE exclusion set.
"""

import argparse
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict


# Curated viral references, sourced by verified NCBI accession.
# Only add an entry after the accession has been confirmed against the
# authoritative record (definition, organism, length). Do not paste
# accessions from memory. See REFERENCE_CURATION.md.
#
#   label -> accession.version
#
# HERV-K113: the real full-length polymorphic HERV-K(HML-2) provirus
# (Turner et al. 2001). NC_022518.1 is the RefSeq record, identical to the
# GenBank submission AY037928.1, "Human endogenous retrovirus K113 complete
# genome", 9,472 bp. It has no chromosomal coordinate because it is an
# unfixed (polymorphic) insertion, so it will not be present in CHM13 unless
# that haplotype happens to carry the allele.
#
# To add later, once each accession is verified: HERV-K106, HERV-K116,
# ciHHV-6A, ciHHV-6B. Leave them out until confirmed.
CURATED_REFERENCES = {
    "HERV-K113": "NC_022518.1",
}

# Below this alignment identity to the CHM13 reference, a curated reference is
# treated as likely polymorphic or absent from CHM13. Its best T2T hit is then
# a paralog, not the reference itself, so it should be used directly in the
# competitive/consensus set rather than as a T2T-extracted exclusion decoy.
POLYMORPHIC_IDENTITY_FLAG = 99.0


def fetch_ncbi_fasta(accession: str, email: str = None, api_key: str = None,
                     retries: int = 3) -> str:
    """Fetch a nucleotide sequence from NCBI by accession. Returns the
    uppercase sequence string (no header)."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": "virome-qc-build_eve_exclusion",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = base + "?" + urllib.parse.urlencode(params)

    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                data = resp.read().decode()
            lines = data.strip().splitlines()
            if not lines or not lines[0].startswith(">"):
                raise ValueError(
                    f"unexpected NCBI response for {accession}: {data[:200]!r}"
                )
            seq = "".join(line.strip() for line in lines[1:])
            if not seq:
                raise ValueError(f"empty sequence for {accession}")
            return seq.upper()
        except Exception as e:  # noqa: BLE001 - report and retry
            last_err = e
            if attempt < retries - 1:
                time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(f"failed to fetch {accession} from NCBI: {last_err}")


def build_curated_query_fasta(output_path: str, email: str = None,
                              api_key: str = None) -> list:
    """Fetch every curated reference and write a query FASTA. Returns a list of
    (query_name, accession, length)."""
    written = []
    # NCBI etiquette: <=3 requests/sec without an API key, <=10 with one.
    delay = 0.11 if api_key else 0.34
    with open(output_path, "w") as f:
        for label, accession in CURATED_REFERENCES.items():
            print(f"Fetching {label} from NCBI ({accession})...")
            seq = fetch_ncbi_fasta(accession, email, api_key)
            query_name = f"{label}_{accession}"
            f.write(f">{query_name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
            written.append((query_name, accession, len(seq)))
            time.sleep(delay)
    print(f"Curated references: {output_path} "
          f"({len(written)} sequences)")
    return written


def combine_query_fastas(curated_path: str, extra_path: str, combined_path: str):
    """Concatenate the curated query FASTA with an optional user-provided
    references FASTA into a single alignment query file."""
    with open(combined_path, "w") as out:
        for src in (curated_path, extra_path):
            if src and os.path.exists(src):
                with open(src) as f:
                    for line in f:
                        out.write(line)
                if not line.endswith("\n"):
                    out.write("\n")


def run_minimap2(t2t_ref: str, query_refs: str, output_paf: str, threads: int = 8):
    """Align viral references (query) to T2T-CHM13 (target)."""
    cmd = [
        "minimap2",
        "-x", "asm20",       # divergent sequence preset (up to 20% divergence)
        "-N", "1000",         # report up to 1000 secondary alignments
        "--secondary=yes",    # keep secondaries for multi-copy EVEs
        "-t", str(threads),
        t2t_ref,
        query_refs,
    ]
    print(f"Running: {' '.join(cmd)}")
    with open(output_paf, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"minimap2 error: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    print(f"Alignment complete: {output_paf}")


def parse_paf(paf_path: str, min_identity: float, min_length: int):
    """Parse PAF file and filter alignments."""
    regions = []
    with open(paf_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue

            query_name = fields[0]
            query_len = int(fields[1])
            target_name = fields[5]
            target_start = int(fields[7])
            target_end = int(fields[8])
            matches = int(fields[9])
            block_len = int(fields[10])

            if block_len < min_length:
                continue

            identity = matches / block_len * 100
            if identity < min_identity:
                continue

            regions.append({
                "chrom": target_name,
                "start": target_start,
                "end": target_end,
                "query": query_name,
                "identity": identity,
                "length": block_len,
                "query_len": query_len,
            })

    print(f"Filtered alignments: {len(regions)} "
          f"(>={min_identity}% identity, >={min_length}bp)")
    return regions


def merge_regions(regions: list, merge_distance: int = 1000):
    """Merge overlapping/adjacent regions on the same chromosome. The merged
    region keeps the query with the single best identity as its best_match,
    plus the full set of queries that overlapped it."""
    by_chrom = defaultdict(list)
    for r in regions:
        by_chrom[r["chrom"]].append(r)

    merged = []
    for chrom in sorted(by_chrom.keys()):
        intervals = sorted(by_chrom[chrom], key=lambda x: x["start"])
        current = intervals[0].copy()
        current["best_query"] = current["query"]
        current["best_identity"] = current["identity"]
        current["queries"] = {current["query"]}

        for r in intervals[1:]:
            if r["start"] <= current["end"] + merge_distance:
                current["end"] = max(current["end"], r["end"])
                current["queries"].add(r["query"])
                if r["identity"] > current["best_identity"]:
                    current["best_identity"] = r["identity"]
                    current["best_query"] = r["query"]
            else:
                merged.append(current)
                current = r.copy()
                current["best_query"] = current["query"]
                current["best_identity"] = current["identity"]
                current["queries"] = {current["query"]}
        merged.append(current)

    print(f"Merged regions: {len(merged)} (from {len(regions)} alignments)")
    return merged


def locus_name(index: int, region: dict) -> str:
    """Name an extracted region by its T2T locus only. The query that matched
    is recorded separately as an annotation, never fused into the name, so a
    paralogous hit is not mislabeled with the query's identity."""
    return f"eve_{index:04d}_{region['chrom']}_{region['start']}_{region['end']}"


def write_bed(regions: list, output_path: str):
    """Write BED file of EVE regions with best-match annotation."""
    with open(output_path, "w") as f:
        f.write("#chrom\tstart\tend\tname\tbest_match\tbest_identity\tlength\n")
        for i, r in enumerate(regions):
            f.write(
                f"{r['chrom']}\t{r['start']}\t{r['end']}\t{locus_name(i, r)}\t"
                f"{r['best_query']}\t{r['best_identity']:.1f}\t"
                f"{r['end'] - r['start']}\n"
            )
    print(f"BED file: {output_path} ({len(regions)} regions)")


def extract_sequences(regions: list, t2t_ref: str, output_fasta: str):
    """Extract T2T-CHM13 subsequences at EVE region coordinates. Header carries
    the locus name plus a best_match/identity annotation (not fused into the
    name)."""
    try:
        import pysam
    except ImportError:
        print("pysam required for FASTA extraction. Install: pip install pysam",
              file=sys.stderr)
        sys.exit(1)

    ref = pysam.FastaFile(t2t_ref)
    total_bp = 0

    with open(output_fasta, "w") as f:
        for i, r in enumerate(regions):
            seq = ref.fetch(r["chrom"], r["start"], r["end"])
            name = locus_name(i, r)
            f.write(
                f">{name} best_match={r['best_query']} "
                f"identity={r['best_identity']:.1f}\n{seq}\n"
            )
            total_bp += len(seq)

    ref.close()
    print(f"FASTA: {output_fasta} ({len(regions)} sequences, {total_bp:,} bp)")
    return total_bp


def write_rust_source(fasta_path: str, output_rs: str):
    """Generate Rust source with embedded EVE sequences. The Rust name is the
    locus name (first whitespace-delimited token of the FASTA header)."""
    sequences = []
    current_name = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences.append((current_name, "".join(current_seq)))
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_name:
            sequences.append((current_name, "".join(current_seq)))

    total_bp = sum(len(s) for _, s in sequences)

    with open(output_rs, "w") as f:
        f.write("//! Embedded high-identity EVE sequences from T2T-CHM13\n")
        f.write("//!\n")
        f.write("//! Generated by scripts/build_eve_exclusion.py\n")
        f.write(f"//! {len(sequences)} sequences, {total_bp:,} bp total\n")
        f.write("//!\n")
        f.write("//! These are regions of the human reference genome that share\n")
        f.write("//! high identity with curated viral references. Their k-mers are\n")
        f.write("//! excluded from host containment calculations to prevent false\n")
        f.write("//! removal of genuinely exogenous viral reads.\n")
        f.write("//!\n")
        f.write("//! Each entry is named by its T2T locus. The viral reference that\n")
        f.write("//! matched the locus, and the identity, are recorded in the\n")
        f.write("//! accompanying BED file, not in the name.\n\n")
        f.write("/// High-identity EVE sequences from T2T-CHM13.\n")
        f.write("/// Each entry is (name, sequence_bytes).\n")
        f.write("pub const EVE_SEQUENCES: &[(&str, &[u8])] = &[\n")

        for name, seq in sequences:
            f.write(f'    ("{name}",\n')
            f.write("     b\"")
            for i in range(0, len(seq), 80):
                chunk = seq[i:i + 80]
                if i > 0:
                    f.write("\\\n       ")
                f.write(chunk)
            f.write("\"),\n")

        f.write("];\n")

    print(f"Rust source: {output_rs} ({len(sequences)} sequences, {total_bp:,} bp)")


def best_raw_hits(paf_path: str) -> dict:
    """Best alignment identity per query from the UNFILTERED PAF. Used to report
    references whose hits fall below the extraction thresholds, e.g. a
    polymorphic provirus that is absent from the assembly and so aligns only to
    diverged paralogs."""
    best = {}
    with open(paf_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            target = fields[5]
            tstart, tend = int(fields[7]), int(fields[8])
            matches, block = int(fields[9]), int(fields[10])
            if block == 0:
                continue
            identity = matches / block * 100
            if query not in best or identity > best[query][0]:
                best[query] = (identity, target, tstart, tend)
    return best


def report_curated_coverage(curated: list, best_by_query: dict):
    """For each curated reference, report its best identity to CHM13 (from the
    unfiltered alignments) and flag references that are likely polymorphic or
    absent from the assembly."""
    if not curated:
        return
    print("\n=== Curated reference coverage in CHM13 ===")
    for query_name, accession, length in curated:
        hit = best_by_query.get(query_name)
        if hit is None:
            print(f"  {query_name}: no alignment to CHM13 at all. Absent from "
                  f"CHM13; use directly in the competitive set (set C).")
            continue
        identity, chrom, start, end = hit
        loc = f"{chrom}:{start}-{end}"
        if identity < POLYMORPHIC_IDENTITY_FLAG:
            print(f"  {query_name}: best T2T hit {identity:.1f}% at {loc}. "
                  f"Below {POLYMORPHIC_IDENTITY_FLAG}% -> likely polymorphic or "
                  f"absent from CHM13; the hit is a paralog, not this reference. "
                  f"Use {accession} directly in the competitive set (set C), not "
                  f"as a T2T decoy. See REFERENCE_CURATION.md.")
        else:
            print(f"  {query_name}: best T2T hit {identity:.1f}% at {loc} "
                  f"(present in CHM13).")


def main():
    parser = argparse.ArgumentParser(description="Build EVE exclusion set for virome-qc")
    parser.add_argument("--t2t", required=True, help="T2T-CHM13 reference FASTA")
    parser.add_argument("--viral-refs", default=None,
                        help="Optional extra viral references FASTA, aligned "
                             "alongside the curated (NCBI-fetched) references")
    parser.add_argument("--fetch-curated", action="store_true", default=True,
                        help="Fetch curated references from NCBI (default: on)")
    parser.add_argument("--no-fetch-curated", dest="fetch_curated",
                        action="store_false",
                        help="Skip NCBI fetch; use only --viral-refs")
    parser.add_argument("--email", default=None,
                        help="Contact email for NCBI E-utilities (recommended)")
    parser.add_argument("--api-key", default=None,
                        help="NCBI API key (raises the rate limit)")
    parser.add_argument("--output-bed", default="references/t2t_eve_regions.bed",
                        help="Output BED file")
    parser.add_argument("--output-fasta", default="references/eve_exclusion_sequences.fa",
                        help="Output FASTA of T2T subsequences")
    parser.add_argument("--output-rust", default="src/modules/eve_kmers.rs",
                        help="Output Rust source file")
    parser.add_argument("--output-curated-fasta",
                        default="references/curated_viral_refs.fa",
                        help="Where to save the fetched curated references "
                             "(reusable for the competitive/consensus set)")
    parser.add_argument("--min-identity", type=float, default=85,
                        help="Minimum alignment identity (%%)")
    parser.add_argument("--min-length", type=int, default=100,
                        help="Minimum alignment length (bp)")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_bed) or ".", exist_ok=True)

    # Step 0: Fetch curated references and assemble the alignment query set.
    curated = []
    if args.fetch_curated:
        os.makedirs(os.path.dirname(args.output_curated_fasta) or ".", exist_ok=True)
        curated = build_curated_query_fasta(
            args.output_curated_fasta, args.email, args.api_key)
    elif not args.viral_refs:
        print("Nothing to align: --no-fetch-curated set and no --viral-refs given.",
              file=sys.stderr)
        sys.exit(1)

    query_fasta = args.output_bed.replace(".bed", ".query.fa")
    combine_query_fastas(
        args.output_curated_fasta if args.fetch_curated else None,
        args.viral_refs,
        query_fasta,
    )

    # Step 1: Align
    paf_path = args.output_bed.replace(".bed", ".paf")
    run_minimap2(args.t2t, query_fasta, paf_path, args.threads)

    # Best raw (unfiltered) hit per query, for reporting polymorphic/absent refs.
    raw_best = best_raw_hits(paf_path)

    # Step 2: Filter and parse
    regions = parse_paf(paf_path, args.min_identity, args.min_length)
    if not regions:
        print(f"No regions passed the extraction thresholds "
              f"(>={args.min_identity}% identity, >={args.min_length}bp).")
        report_curated_coverage(curated, raw_best)
        print("\nNo exclusion set written; existing outputs left unchanged. A "
              "polymorphic reference with no high-identity host copy is an "
              "expected outcome, not an error: it belongs in the competitive "
              "set (already saved to the curated FASTA), not as a T2T decoy.")
        # Exit 0: this is a valid result, not a failure.
        sys.exit(0)

    # Step 3: Merge overlapping
    merged = merge_regions(regions)

    # Step 4: Write BED
    write_bed(merged, args.output_bed)

    # Step 5: Extract sequences
    os.makedirs(os.path.dirname(args.output_fasta) or ".", exist_ok=True)
    extract_sequences(merged, args.t2t, args.output_fasta)

    # Step 6: Generate Rust source
    os.makedirs(os.path.dirname(args.output_rust) or ".", exist_ok=True)
    write_rust_source(args.output_fasta, args.output_rust)

    # Step 7: Report curated-reference coverage and polymorphic flags
    report_curated_coverage(curated, raw_best)

    # Summary
    total_bp = sum(r["end"] - r["start"] for r in merged)
    print("\n=== Summary ===")
    print(f"EVE regions: {len(merged)}")
    print(f"Total bp: {total_bp:,}")
    print(f"Estimated k=31 k-mers: ~{total_bp * 2:,} (both strands)")
    print(f"Estimated FxHashSet memory: ~{total_bp * 2 * 8 / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
