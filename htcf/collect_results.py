#!/usr/bin/env python3
"""Collect and summarize benchmark results from HTCF runs."""

import json
import os
from pathlib import Path

RESULTS_DIR = Path("/scratch/sahlab/shandley/virome-qc-benchmark/results")

def main():
    if not RESULTS_DIR.exists():
        # Fall back to local path for testing
        results = Path("benchmark_data/results")
    else:
        results = RESULTS_DIR

    datasets = sorted(d for d in results.iterdir() if d.is_dir() and (d / "passport.json").exists())

    print(f"{'Dataset':30s} {'Reads':>10s} {'Surv':>7s} {'Host%':>8s} {'rRNA%':>8s} {'Contam%':>8s} {'ERV':>6s} {'GC':>6s} {'Tier':>5s}")
    print("-" * 95)

    for d in datasets:
        passport = json.loads((d / "passport.json").read_text())
        ri = passport.get("reads_input", 1)
        surv = passport.get("survival_rate", 0)

        host = sum(m.get("reads_removed", 0) for m in passport.get("modules", []) if m.get("name") == "host")
        rrna = sum(m.get("reads_removed", 0) for m in passport.get("modules", []) if m.get("name") == "rrna")
        contam = sum(m.get("reads_removed", 0) for m in passport.get("modules", []) if m.get("name") == "contaminant")
        erv = passport.get("erv_analysis", {}).get("retroviral_reads_flagged", 0)
        gc = passport.get("ingestion", {}).get("mean_gc", 0)
        tier = passport.get("quality_tier", "?")

        print(
            f"{d.name:30s} {ri:>10,} {surv*100:>6.1f}% "
            f"{host/ri*100:>7.3f}% {rrna/ri*100:>7.3f}% {contam/ri*100:>7.3f}% "
            f"{erv:>6,} {gc:>5.2f} {tier:>5s}"
        )

    print(f"\nTotal: {len(datasets)} datasets processed")

if __name__ == "__main__":
    main()
