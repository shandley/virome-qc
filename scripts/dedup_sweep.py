#!/usr/bin/env python3
"""
Dedup module validation using ViroForge PCR duplicate injection.

Tests streaming dedup accuracy at multiple duplicate rates using
ViroForge source labels (pcr_duplicate=true) as exact ground truth.
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
OUTPUT_DIR = VIROME_QC_DIR / "benchmark_data" / "dedup_sweep"
VIROFORGE_PYTHON = VIROFORGE_DIR / ".venv" / "bin" / "python3"
GENERATE_SCRIPT = VIROFORGE_DIR / "scripts" / "generate_fastq_dataset.py"
DB_PATH = VIROFORGE_DIR / "viroforge" / "data" / "viral_genomes.db"

DUPLICATE_RATES = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]
SEEDS = [42, 43]
COVERAGE = 3


def generate_dataset(dup_rate, seed):
    tag = f"dup{int(dup_rate*100):02d}_seed{seed}"
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
        "--collection-id", "9", "--platform", "novaseq",
        "--coverage", str(COVERAGE), "--output", str(out_dir), "--seed", str(seed),
    ]
    if dup_rate > 0:
        cmd.extend(["--duplicate-rate", str(dup_rate)])

    print(f"  Generating dup_rate={dup_rate:.0%} seed={seed}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600, env=env)
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-300:]}")
        return None

    for fq in fastq_dir.glob("*.fastq"):
        subprocess.run(["gzip", str(fq)], check=True)
    for f in (out_dir / "fasta").glob("*"):
        f.unlink()

    r1_files = list(fastq_dir.glob("*_R1.fastq.gz"))
    return r1_files[0] if r1_files else None


def read_labels(r1_path):
    """Read source and duplicate labels from headers."""
    labels = {}
    with gzip.open(r1_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                parts = line.strip().split()
                read_id = parts[0][1:]
                is_dup = any("pcr_duplicate=true" in p for p in parts)
                labels[read_id] = is_dup
    return labels


def run_virome_qc(r1_path, r2_path, output_dir):
    """Run virome-qc with only dedup enabled."""
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = output_dir / "profile.yaml"

    with open(profile_path, "w") as f:
        f.write('name: "dedup_sweep"\ndescription: "test"\nplatform: illumina\n')
        f.write('modules:\n')
        f.write('  adapter:\n    enabled: false\n    sequences: []\n    internal_scan: false\n')
        f.write('    random_primer_trim: 0\n    min_overlap: 8\n    max_mismatch_rate: 0.1\n')
        f.write('  quality:\n    enabled: false\n    window_size: 10\n    min_quality: 15\n')
        f.write('    min_mean_quality: 20.0\n    min_length: 30\n')
        f.write('  polyx:\n    enabled: false\n    platform_aware: false\n    min_length: 8\n')
        f.write('  complexity:\n    enabled: false\n    min_entropy: 0.5\n')
        f.write('  contaminant:\n    enabled: false\n    screen_rrna: false\n    screen_phix: false\n')
        f.write('    screen_vectors: false\n    screen_kitome: false\n    min_kmer_fraction: 0.25\n')
        f.write('  host:\n    enabled: false\n    reference: human\n    eve_aware: false\n    rescue: false\n')
        f.write('  dedup:\n    enabled: true\n    optical_distance: 2500\n    umi_aware: false\n')
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
    print("Dedup Module Validation Sweep")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for dup_rate in DUPLICATE_RATES:
        print(f"\n--- Duplicate rate: {dup_rate:.0%} ---")

        for seed in SEEDS:
            r1_path = generate_dataset(dup_rate, seed)
            if r1_path is None:
                continue
            r2_path = Path(str(r1_path).replace("_R1.", "_R2."))

            labels = read_labels(r1_path)
            n_dup = sum(1 for v in labels.values() if v)
            n_unique = sum(1 for v in labels.values() if not v)
            total_r1 = len(labels)
            print(f"  Dataset: {total_r1:,} R1 reads ({n_dup:,} duplicates, {n_unique:,} unique)")

            qc_dir = OUTPUT_DIR / f"dup{int(dup_rate*100):02d}_seed{seed}_qc"
            passport = run_virome_qc(r1_path, r2_path, qc_dir)
            if not passport:
                continue

            dedup_removed = 0
            for m in passport["modules"]:
                if m["name"] == "dedup":
                    dedup_removed = m["reads_removed"]

            total_pe = passport["reads_input"]
            gt_dup_pe = n_dup * 2
            gt_unique_pe = n_unique * 2

            # TP = duplicates correctly removed
            # FP = unique reads incorrectly removed
            tp = min(dedup_removed, gt_dup_pe)
            fp = max(0, dedup_removed - gt_dup_pe)
            fn = max(0, gt_dup_pe - tp)
            tn = gt_unique_pe - fp

            sens = tp / max(tp + fn, 1)
            spec = tn / max(tn + fp, 1)
            prec = tp / max(tp + fp, 1)
            fpr = fp / max(fp + tn, 1)
            f1 = 2 * prec * sens / max(prec + sens, 1e-10)

            result = {
                "duplicate_rate": dup_rate,
                "seed": seed,
                "total_reads": total_pe,
                "gt_duplicates": gt_dup_pe,
                "gt_unique": gt_unique_pe,
                "dedup_removed": dedup_removed,
                "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                "sensitivity": sens,
                "specificity": spec,
                "precision": prec,
                "fpr": fpr,
                "f1": f1,
            }
            all_results.append(result)

            print(f"    Removed={dedup_removed:,} Sens={sens:.3f} Spec={spec:.4f} "
                  f"FPR={fpr:.5f} Prec={prec:.3f} F1={f1:.3f}")

    results_path = OUTPUT_DIR / "dedup_sweep_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)

    print("\n" + "=" * 80)
    print(f"{'Dup Rate':>9s} {'GT Dups':>9s} {'Removed':>9s} "
          f"{'Sens':>8s} {'Spec':>8s} {'FPR':>10s} {'Prec':>8s} {'F1':>8s}")
    print("-" * 80)

    for dup_rate in DUPLICATE_RATES:
        subset = [r for r in all_results if r["duplicate_rate"] == dup_rate]
        if not subset:
            continue
        print(
            f"{dup_rate:>8.0%} "
            f"{np.mean([r['gt_duplicates'] for r in subset]):>9.0f} "
            f"{np.mean([r['dedup_removed'] for r in subset]):>9.0f} "
            f"{np.mean([r['sensitivity'] for r in subset]):>8.3f} "
            f"{np.mean([r['specificity'] for r in subset]):>8.4f} "
            f"{np.mean([r['fpr'] for r in subset]):>10.5f} "
            f"{np.mean([r['precision'] for r in subset]):>8.3f} "
            f"{np.mean([r['f1'] for r in subset]):>8.3f}"
        )


if __name__ == "__main__":
    main()
