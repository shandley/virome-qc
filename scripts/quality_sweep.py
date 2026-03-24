#!/usr/bin/env python3
"""
Quality trimming parameter sweep using ViroForge synthetic data.

Quality trimming doesn't have a traditional TP/FP framework because there's
no "ground truth bad quality" -- ISS models realistic quality degradation.
Instead, we measure the tradeoff between:
  - Viral read retention (higher is better)
  - Output quality (higher is better)
  - Bases retained (higher is better)

This produces a tradeoff curve, not a ROC curve.
"""

import gzip
import json
import os
import subprocess
from pathlib import Path

import numpy as np

VIROFORGE_DIR = Path(os.environ.get("VIROFORGE_DIR", "/Users/scotthandley/Code/tools/viroforge"))
VIROME_QC_DIR = Path(__file__).parent.parent
VIROME_QC_BIN = Path(os.environ.get(
    "VIROME_QC_BIN",
    Path.home() / ".cargo" / "target" / "release" / "virome-qc"
))
OUTPUT_DIR = VIROME_QC_DIR / "benchmark_data" / "quality_sweep"
VIROFORGE_PYTHON = VIROFORGE_DIR / ".venv" / "bin" / "python3"
GENERATE_SCRIPT = VIROFORGE_DIR / "scripts" / "generate_fastq_dataset.py"
DB_PATH = VIROFORGE_DIR / "viroforge" / "data" / "viral_genomes.db"

# Sweep parameters
MEAN_QUALITY_THRESHOLDS = [0, 15, 18, 20, 22, 25, 28, 30, 32]
WINDOW_SIZES = [4, 10, 15, 25]
MIN_QUALITY_VALUES = [10, 15, 20, 25]
PLATFORMS = [
    ("novaseq", 9),   # NovaSeq gut virome
    ("miseq", 9),     # MiSeq gut virome (different error model)
]
SEEDS = [42, 43]  # 2 replicates for speed
COVERAGE = 3


def generate_dataset(platform, seed):
    tag = f"{platform}_seed{seed}"
    out_dir = OUTPUT_DIR / tag
    fastq_dir = out_dir / "fastq"

    existing = list(fastq_dir.glob("*_R1.fastq.gz")) if fastq_dir.exists() else []
    if existing:
        return existing[0]

    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = str(VIROFORGE_DIR / ".venv" / "bin") + ":" + env.get("PATH", "")

    cmd = [
        str(VIROFORGE_PYTHON), str(GENERATE_SCRIPT),
        "--database", str(DB_PATH),
        "--collection-id", "9", "--platform", platform,
        "--coverage", str(COVERAGE), "--output", str(out_dir), "--seed", str(seed),
    ]

    print(f"  Generating {platform} seed={seed}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-200:]}")
        return None

    for fq in fastq_dir.glob("*.fastq"):
        subprocess.run(["gzip", str(fq)], check=True)
    for f in (out_dir / "fasta").glob("*"):
        f.unlink()

    r1_files = list(fastq_dir.glob("*_R1.fastq.gz"))
    return r1_files[0] if r1_files else None


def read_source_labels(r1_path):
    labels = {}
    with gzip.open(r1_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                parts = line.strip().split()
                read_id = parts[0][1:]
                source = "unknown"
                for part in parts[1:]:
                    if part.startswith("source="):
                        source = part.split("=")[1]
                labels[read_id] = source
    return labels


def run_virome_qc(r1_path, r2_path, output_dir, min_mean_quality, window_size, min_quality):
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = output_dir / "profile.yaml"

    # Only enable quality module
    with open(profile_path, "w") as f:
        f.write(f'name: "quality_sweep"\ndescription: "test"\nplatform: illumina\n')
        f.write('modules:\n')
        f.write('  adapter:\n    enabled: false\n    sequences: []\n    internal_scan: false\n')
        f.write('    random_primer_trim: 0\n    min_overlap: 8\n    max_mismatch_rate: 0.1\n')
        f.write(f'  quality:\n    enabled: true\n    window_size: {window_size}\n')
        f.write(f'    min_quality: {min_quality}\n    min_mean_quality: {min_mean_quality}\n')
        f.write('    min_length: 30\n')
        f.write('  polyx:\n    enabled: false\n    platform_aware: false\n    min_length: 8\n')
        f.write('  complexity:\n    enabled: false\n    min_entropy: 0.5\n')
        f.write('  contaminant:\n    enabled: false\n    screen_rrna: false\n    screen_phix: false\n')
        f.write('    screen_vectors: false\n    screen_kitome: false\n    min_kmer_fraction: 0.4\n')
        f.write('  host:\n    enabled: false\n    reference: human\n    eve_aware: false\n    rescue: false\n')
        f.write('  dedup:\n    enabled: false\n    optical_distance: 2500\n    umi_aware: false\n')
        f.write('  chimera:\n    enabled: false\n')
        f.write('thresholds:\n  min_survival_rate: 0.01\n  max_host_fraction: 0.99\n')
        f.write('  max_rrna_fraction: 0.99\n  max_duplicate_rate: 0.99\n')

    cmd = [
        str(VIROME_QC_BIN), "run", "-p", str(profile_path),
        "-1", str(r1_path), "-2", str(r2_path),
        "-i", ".", "-o", str(output_dir), "--report-only", "-t", "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    passport_path = output_dir / "passport.json"
    if passport_path.exists():
        with open(passport_path) as f:
            return json.load(f)
    return {}


def main():
    print("=" * 70)
    print("Quality Trimming Parameter Sweep")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for platform, collection_id in PLATFORMS:
        print(f"\n--- Platform: {platform} ---")

        for seed in SEEDS:
            r1_path = generate_dataset(platform, seed)
            if r1_path is None:
                continue
            r2_path = Path(str(r1_path).replace("_R1.", "_R2."))
            labels = read_source_labels(r1_path)
            n_viral = sum(1 for s in labels.values() if s == "viral")
            print(f"  Dataset: {len(labels):,} reads ({n_viral:,} viral)")

            # Sweep min_mean_quality with fixed window_size=4 and min_quality=15
            for mmq in MEAN_QUALITY_THRESHOLDS:
                tag = f"{platform}_s{seed}_mmq{mmq}_w4_q15"
                qc_dir = OUTPUT_DIR / tag

                passport = run_virome_qc(r1_path, r2_path, qc_dir, mmq, 4, 15)
                if not passport:
                    continue

                total = passport["reads_input"]
                passed = passport["reads_passed"]
                quality_mod = None
                for m in passport["modules"]:
                    if m["name"] == "quality":
                        quality_mod = m

                qa = passport.get("qa_stats", {})
                bases_out = qa.get("summary", {}).get("bases_output", 0)
                bases_in = qa.get("summary", {}).get("bases_input", 0)

                result = {
                    "platform": platform,
                    "seed": seed,
                    "min_mean_quality": mmq,
                    "window_size": 4,
                    "min_quality": 15,
                    "reads_input": total,
                    "reads_passed": passed,
                    "survival_rate": passed / max(total, 1),
                    "reads_removed": quality_mod["reads_removed"] if quality_mod else 0,
                    "bases_in": bases_in,
                    "bases_out": bases_out,
                    "base_retention": bases_out / max(bases_in, 1),
                    "viral_total": n_viral * 2,
                }
                all_results.append(result)

            # Also sweep window_size with fixed mmq=20
            for ws in WINDOW_SIZES:
                if ws == 4:
                    continue  # already covered above
                tag = f"{platform}_s{seed}_mmq20_w{ws}_q15"
                qc_dir = OUTPUT_DIR / tag

                passport = run_virome_qc(r1_path, r2_path, qc_dir, 20, ws, 15)
                if not passport:
                    continue

                total = passport["reads_input"]
                passed = passport["reads_passed"]
                quality_mod = None
                for m in passport["modules"]:
                    if m["name"] == "quality":
                        quality_mod = m

                qa = passport.get("qa_stats", {})
                bases_out = qa.get("summary", {}).get("bases_output", 0)
                bases_in = qa.get("summary", {}).get("bases_input", 0)

                result = {
                    "platform": platform,
                    "seed": seed,
                    "min_mean_quality": 20,
                    "window_size": ws,
                    "min_quality": 15,
                    "reads_input": total,
                    "reads_passed": passed,
                    "survival_rate": passed / max(total, 1),
                    "reads_removed": quality_mod["reads_removed"] if quality_mod else 0,
                    "bases_in": bases_in,
                    "bases_out": bases_out,
                    "base_retention": bases_out / max(bases_in, 1),
                    "viral_total": n_viral * 2,
                }
                all_results.append(result)

    # Save
    results_path = OUTPUT_DIR / "quality_sweep_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)

    # Summary
    print("\n" + "=" * 90)
    print("Quality tradeoff curve: min_mean_quality sweep (window=4, min_q=15)")
    print(f"{'Platform':>10s} {'MMQ':>5s} {'Removed':>10s} {'Survival':>10s} "
          f"{'Base Ret':>10s}")
    print("-" * 55)

    for platform in ["novaseq", "miseq"]:
        for mmq in MEAN_QUALITY_THRESHOLDS:
            subset = [r for r in all_results
                      if r["platform"] == platform
                      and r["min_mean_quality"] == mmq
                      and r["window_size"] == 4]
            if not subset:
                continue
            marker = " <-- default" if mmq == 20 else ""
            print(
                f"{platform:>10s} Q{mmq:>3d} "
                f"{np.mean([r['reads_removed'] for r in subset]):>10.0f} "
                f"{np.mean([r['survival_rate'] for r in subset])*100:>9.2f}% "
                f"{np.mean([r['base_retention'] for r in subset])*100:>9.2f}%"
                f"{marker}"
            )

    print("\n" + "=" * 90)
    print("Window size sweep (mmq=20, min_q=15)")
    print(f"{'Platform':>10s} {'Window':>7s} {'Removed':>10s} {'Survival':>10s} "
          f"{'Base Ret':>10s}")
    print("-" * 55)

    for platform in ["novaseq", "miseq"]:
        for ws in WINDOW_SIZES:
            subset = [r for r in all_results
                      if r["platform"] == platform
                      and r["min_mean_quality"] == 20
                      and r["window_size"] == ws]
            if not subset:
                continue
            print(
                f"{platform:>10s} {ws:>7d} "
                f"{np.mean([r['reads_removed'] for r in subset]):>10.0f} "
                f"{np.mean([r['survival_rate'] for r in subset])*100:>9.2f}% "
                f"{np.mean([r['base_retention'] for r in subset])*100:>9.2f}%"
            )


if __name__ == "__main__":
    main()
