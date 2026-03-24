#!/usr/bin/env python3
"""
Regenerate all benchmark reports with the current virome-qc code.

Downloads each dataset from SRA/ENA, runs virome-qc with --report-only,
generates the React HTML report, and cleans up raw FASTQ files.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

VIROME_QC_BIN = Path.home() / ".cargo" / "target" / "release" / "virome-qc"
VIROME_QC_DIR = Path(__file__).parent.parent
RESULTS_DIR = VIROME_QC_DIR / "benchmark_data" / "results"
REPORT_CSS = VIROME_QC_DIR / "report-ui" / "dist" / "report.css"
REPORT_JS = VIROME_QC_DIR / "report-ui" / "dist" / "report.js"

# Dataset definitions: (name, accession, paired, profile, fastq_dump_args)
DATASETS = [
    # Real datasets
    ("cook", "ERR10359658", True, "stool-vlp-tagmentation", []),
    ("kleiner", "ERR1877475", True, "profiles/short-read-nebnext.yaml", []),
    ("buddle_wgs", "ERR13480651", False, "profiles/short-read-nebnext.yaml", ["-X", "5000000"]),
    ("buddle_rna", "ERR13480663", False, "profiles/short-read-nebnext.yaml", ["-X", "5000000"]),
    ("zhang_undepleted", "SRR33419012", True, "profiles/short-read-nebnext.yaml", ["-X", "1000000"]),
    ("zhang_depleted", "SRR33419066", True, "profiles/short-read-nebnext.yaml", ["-X", "1000000"]),
    ("santos_kapa", "SRR8487022", True, "stool-vlp-tagmentation", []),
    ("santos_nextera", "SRR8487034", True, "stool-vlp-tagmentation", []),
    ("shkoporov_gut", "SRR9161520", True, "stool-vlp-tagmentation", []),
    ("tara_ocean", "ERR599370", True, "profiles/short-read-nebnext.yaml", []),
    ("chrisman_dnbseq", "ERR9765742", False, "profiles/short-read-nebnext.yaml", ["-X", "10000000"]),
]


def download_dataset(name, accession, paired, extra_args):
    """Download dataset from SRA using fastq-dump."""
    tmp_dir = VIROME_QC_DIR / "benchmark_data" / "tmp_download"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    r1 = tmp_dir / f"{accession}_1.fastq.gz"
    r2 = tmp_dir / f"{accession}_2.fastq.gz"

    if r1.exists():
        return (r1, r2 if r2.exists() and paired else None)

    cmd = ["fastq-dump", accession, "--split-files", "--gzip", "--clip"]
    cmd.extend(extra_args)

    print(f"  Downloading {accession}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd=str(tmp_dir))
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-200:]}")
        return (None, None)

    # Find output files
    r1_candidates = list(tmp_dir.glob(f"{accession}_1.fastq.gz"))
    r2_candidates = list(tmp_dir.glob(f"{accession}_2.fastq.gz"))

    r1_path = r1_candidates[0] if r1_candidates else None
    r2_path = r2_candidates[0] if r2_candidates and paired else None

    # If fastq-dump only produced one file (single-end or -X mode)
    if not r1_path:
        single = list(tmp_dir.glob(f"{accession}*.fastq.gz"))
        if single:
            r1_path = single[0]

    return (r1_path, r2_path)


def run_virome_qc(name, r1, r2, profile):
    """Run virome-qc with --report-only."""
    output_dir = RESULTS_DIR / name
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [str(VIROME_QC_BIN), "run", "-p", profile]

    if r2:
        cmd.extend(["-1", str(r1), "-2", str(r2), "-i", ".", "--merge"])
    else:
        cmd.extend(["-i", str(r1)])

    cmd.extend(["-o", str(output_dir), "--report-only", "-t", "8"])

    print(f"  Running virome-qc on {name}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-200:]}")
        return False

    # Print summary
    for line in result.stderr.split("\n"):
        if "QC complete" in line or "Platform" in line or "ERV analysis" in line:
            print(f"    {line.strip()}")

    return (output_dir / "passport.json").exists()


def generate_html_report(name):
    """Generate React HTML report from passport."""
    output_dir = RESULTS_DIR / name
    passport_path = output_dir / "passport.json"
    report_path = output_dir / "report.html"

    if not passport_path.exists():
        return False

    # virome-qc now generates reports automatically with --report-only
    # But let's verify it exists
    return report_path.exists()


def main():
    print("=" * 70)
    print("Regenerating all benchmark reports with current code")
    print("=" * 70)

    # Verify virome-qc is built
    if not VIROME_QC_BIN.exists():
        print("Building virome-qc...")
        subprocess.run(["cargo", "build", "--release"], cwd=str(VIROME_QC_DIR), check=True)

    successful = []
    failed = []

    for name, accession, paired, profile, extra_args in DATASETS:
        print(f"\n--- {name} ({accession}) ---")

        # Download
        r1, r2 = download_dataset(name, accession, paired, extra_args)
        if not r1:
            print(f"  SKIP: download failed")
            failed.append(name)
            continue

        # Run
        success = run_virome_qc(name, r1, r2, profile)
        if not success:
            print(f"  SKIP: virome-qc failed")
            failed.append(name)
            continue

        # Verify report
        if generate_html_report(name):
            successful.append(name)
            print(f"  OK: {name}")
        else:
            failed.append(name)
            print(f"  WARN: passport created but no report")

        # Clean up downloaded files for this dataset
        for f in (VIROME_QC_DIR / "benchmark_data" / "tmp_download").glob(f"{accession}*"):
            f.unlink()

    # Clean up tmp dir
    tmp_dir = VIROME_QC_DIR / "benchmark_data" / "tmp_download"
    if tmp_dir.exists():
        for f in tmp_dir.iterdir():
            f.unlink()
        tmp_dir.rmdir()

    print("\n" + "=" * 70)
    print(f"Results: {len(successful)} successful, {len(failed)} failed")
    print(f"Successful: {', '.join(successful)}")
    if failed:
        print(f"Failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
