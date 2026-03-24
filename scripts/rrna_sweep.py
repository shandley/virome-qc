#!/usr/bin/env python3
"""
rRNA SILVA filter threshold sweep using ViroForge synthetic data.

Sweeps min_kmer_fraction for the SILVA k=21 sorted hash filter.
Uses ViroForge source labels (source=rrna vs source=viral) as ground truth.
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
OUTPUT_DIR = VIROME_QC_DIR / "benchmark_data" / "rrna_sweep"
VIROFORGE_PYTHON = VIROFORGE_DIR / ".venv" / "bin" / "python3"
GENERATE_SCRIPT = VIROFORGE_DIR / "scripts" / "generate_fastq_dataset.py"
DB_PATH = VIROFORGE_DIR / "viroforge" / "data" / "viral_genomes.db"
RRNA_FILTER = VIROME_QC_DIR / "benchmark_data" / "references" / "rrna_silva.rrf"

THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50]
SEEDS = [42, 43, 44]
COVERAGE = 3


def generate_dataset(seed):
    out_dir = OUTPUT_DIR / f"gut_seed{seed}"
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
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
    if result.returncode != 0:
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


def run_virome_qc(r1_path, r2_path, output_dir, min_kmer_fraction):
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = output_dir / "profile.yaml"
    with open(profile_path, "w") as f:
        f.write(f'name: "rrna_sweep"\ndescription: "test"\nplatform: illumina\n')
        f.write('modules:\n')
        f.write('  adapter:\n    enabled: false\n    sequences: []\n    internal_scan: false\n')
        f.write('    random_primer_trim: 0\n    min_overlap: 8\n    max_mismatch_rate: 0.1\n')
        f.write('  quality:\n    enabled: false\n    window_size: 10\n    min_quality: 15\n')
        f.write('    min_mean_quality: 20.0\n    min_length: 30\n')
        f.write('  polyx:\n    enabled: false\n    platform_aware: false\n    min_length: 8\n')
        f.write('  complexity:\n    enabled: false\n    min_entropy: 0.5\n')
        # Disable consensus screener, use SILVA only
        f.write('  contaminant:\n    enabled: false\n    screen_rrna: false\n    screen_phix: false\n')
        f.write('    screen_vectors: false\n    screen_kitome: false\n    min_kmer_fraction: 0.4\n')
        f.write(f'  rrna:\n    enabled: true\n    filter: {RRNA_FILTER}\n')
        f.write(f'    min_kmer_fraction: {min_kmer_fraction}\n')
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
    print("rRNA SILVA Filter Threshold Sweep")
    print("=" * 70)

    if not RRNA_FILTER.exists():
        print(f"ERROR: rRNA filter not found at {RRNA_FILTER}")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for seed in SEEDS:
        r1_path = generate_dataset(seed)
        if r1_path is None:
            continue
        r2_path = Path(str(r1_path).replace("_R1.", "_R2."))
        labels = read_source_labels(r1_path)
        n_rrna = sum(1 for s in labels.values() if s == "rrna")
        n_viral = sum(1 for s in labels.values() if s == "viral")
        print(f"  Dataset seed={seed}: {len(labels):,} reads ({n_viral:,} viral, {n_rrna:,} rRNA)")

        for threshold in THRESHOLDS:
            qc_dir = OUTPUT_DIR / f"seed{seed}_thresh{int(threshold*100):02d}"
            passport = run_virome_qc(r1_path, r2_path, qc_dir, threshold)
            if not passport:
                continue

            total = passport["reads_input"]
            gt_rrna = n_rrna * 2
            gt_viral = n_viral * 2

            rrna_removed = 0
            for m in passport["modules"]:
                if m["name"] == "rrna":
                    rrna_removed = m["reads_removed"]

            tp = min(rrna_removed, gt_rrna)
            fp = max(0, rrna_removed - gt_rrna)
            fn = max(0, gt_rrna - tp)
            tn = gt_viral - fp

            sens = tp / max(tp + fn, 1)
            spec = tn / max(tn + fp, 1)
            fpr = fp / max(fp + tn, 1)
            prec = tp / max(tp + fp, 1)
            f1 = 2 * prec * sens / max(prec + sens, 1e-10)

            result = {
                "seed": seed, "min_kmer_fraction": threshold,
                "gt_rrna": gt_rrna, "gt_viral": gt_viral,
                "rrna_removed": rrna_removed,
                "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                "sensitivity": sens, "specificity": spec,
                "fpr": fpr, "precision": prec, "f1": f1,
            }
            all_results.append(result)

            print(f"    thresh={threshold:.2f}: Sens={sens:.3f} Spec={spec:.4f} "
                  f"FPR={fpr:.5f} Prec={prec:.3f} F1={f1:.3f}")

    results_path = OUTPUT_DIR / "rrna_sweep_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)

    print("\n" + "=" * 80)
    print(f"{'Threshold':>10s} {'Sens':>10s} {'Spec':>10s} {'FPR':>10s} {'Prec':>10s} {'F1':>10s}")
    print("-" * 80)

    for threshold in THRESHOLDS:
        subset = [r for r in all_results if r["min_kmer_fraction"] == threshold]
        if not subset:
            continue
        marker = " <-- current" if abs(threshold - 0.25) < 0.01 else ""
        print(
            f"{threshold:>10.2f} "
            f"{np.mean([r['sensitivity'] for r in subset]):>6.3f}+/-{np.std([r['sensitivity'] for r in subset]):.3f}"
            f"{np.mean([r['specificity'] for r in subset]):>7.4f}+/-{np.std([r['specificity'] for r in subset]):.4f}"
            f"{np.mean([r['fpr'] for r in subset]):>8.5f}+/-{np.std([r['fpr'] for r in subset]):.5f}"
            f"{np.mean([r['precision'] for r in subset]):>7.3f}+/-{np.std([r['precision'] for r in subset]):.3f}"
            f"{np.mean([r['f1'] for r in subset]):>7.3f}+/-{np.std([r['f1'] for r in subset]):.3f}"
            f"{marker}"
        )

    best_j = max(all_results, key=lambda r: r["sensitivity"] + r["specificity"] - 1)
    print(f"\nOptimal Youden's J: threshold={best_j['min_kmer_fraction']:.2f} "
          f"(J={best_j['sensitivity'] + best_j['specificity'] - 1:.4f})")
    best_f1 = max(all_results, key=lambda r: r["f1"])
    print(f"Optimal F1: threshold={best_f1['min_kmer_fraction']:.2f} "
          f"(F1={best_f1['f1']:.3f})")


if __name__ == "__main__":
    main()
