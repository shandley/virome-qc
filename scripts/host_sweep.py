#!/usr/bin/env python3
"""
Host depletion threshold sweep using ViroForge synthetic data.

Sweeps host_threshold and ambiguous_threshold using ViroForge source labels
(source=host_dna vs source=viral) as exact ground truth.

Produces ROC data for the three-way classification (host/ambiguous/keep).
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
OUTPUT_DIR = VIROME_QC_DIR / "benchmark_data" / "host_sweep"
VIROFORGE_PYTHON = VIROFORGE_DIR / ".venv" / "bin" / "python3"
GENERATE_SCRIPT = VIROFORGE_DIR / "scripts" / "generate_fastq_dataset.py"
DB_PATH = VIROFORGE_DIR / "viroforge" / "data" / "viral_genomes.db"
HOST_FILTER = VIROME_QC_DIR / "benchmark_data" / "references" / "human_t2t.sbf"

# Sweep parameters
HOST_THRESHOLDS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70]
AMBIG_THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
SEEDS = [42, 43, 44]
COVERAGE = 3


def generate_dataset(seed: int) -> Path:
    """Generate ViroForge gut virome with host contamination."""
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
        "--collection-id", "9",
        "--platform", "novaseq",
        "--coverage", str(COVERAGE),
        "--output", str(out_dir),
        "--seed", str(seed),
    ]

    print(f"  Generating seed={seed}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-300:]}")
        return None

    for fq in fastq_dir.glob("*.fastq"):
        subprocess.run(["gzip", str(fq)], check=True)
    fasta_dir = out_dir / "fasta"
    if fasta_dir.exists():
        for f in fasta_dir.glob("*"):
            f.unlink()

    r1_files = list(fastq_dir.glob("*_R1.fastq.gz"))
    return r1_files[0] if r1_files else None


def read_source_labels(r1_path: Path) -> dict:
    """Read source labels from FASTQ headers."""
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


def run_virome_qc(r1_path, r2_path, output_dir, host_threshold, ambig_threshold):
    """Run virome-qc with specific host thresholds."""
    output_dir.mkdir(parents=True, exist_ok=True)

    profile_path = output_dir / "profile.yaml"
    with open(profile_path, "w") as f:
        f.write(f'name: "host_sweep"\ndescription: "test"\nplatform: illumina\n')
        f.write('modules:\n')
        f.write('  adapter:\n    enabled: false\n    sequences: []\n    internal_scan: false\n')
        f.write('    random_primer_trim: 0\n    min_overlap: 8\n    max_mismatch_rate: 0.1\n')
        f.write('  quality:\n    enabled: false\n    window_size: 10\n    min_quality: 15\n')
        f.write('    min_mean_quality: 20.0\n    min_length: 30\n')
        f.write('  polyx:\n    enabled: false\n    platform_aware: false\n    min_length: 8\n')
        f.write('  complexity:\n    enabled: false\n    min_entropy: 0.5\n')
        f.write('  contaminant:\n    enabled: false\n    screen_rrna: false\n    screen_phix: false\n')
        f.write('    screen_vectors: false\n    screen_kitome: false\n    min_kmer_fraction: 0.4\n')
        f.write(f'  host:\n    enabled: true\n    reference: {HOST_FILTER}\n')
        f.write(f'    host_threshold: {host_threshold}\n    ambiguous_threshold: {ambig_threshold}\n')
        f.write('    eve_aware: false\n    rescue: false\n')
        f.write('  dedup:\n    enabled: false\n    optical_distance: 2500\n    umi_aware: false\n')
        f.write('  chimera:\n    enabled: false\n')
        f.write('thresholds:\n  min_survival_rate: 0.01\n  max_host_fraction: 0.99\n')
        f.write('  max_rrna_fraction: 0.99\n  max_duplicate_rate: 0.99\n')

    cmd = [
        str(VIROME_QC_BIN), "run",
        "-p", str(profile_path),
        "-1", str(r1_path), "-2", str(r2_path),
        "-i", ".", "-o", str(output_dir),
        "--report-only", "-t", "4",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    passport_path = output_dir / "passport.json"
    if passport_path.exists():
        with open(passport_path) as f:
            return json.load(f)
    return {}


def compute_metrics(passport, source_labels):
    """Compute host classification accuracy."""
    total = passport.get("reads_input", 0)

    gt_host = sum(1 for s in source_labels.values() if s == "host_dna") * 2
    gt_viral = sum(1 for s in source_labels.values() if s == "viral") * 2
    gt_other = total - gt_host - gt_viral  # rRNA, phix, reagent

    host_mod = None
    for m in passport.get("modules", []):
        if m["name"] == "host":
            host_mod = m
            break

    if not host_mod:
        return {}

    host_removed = host_mod["reads_removed"]
    ambig_flagged = host_mod["extra"].get("ambiguous_flagged", 0)

    # For host classification ROC:
    # TP = host reads correctly removed
    # FP = non-host reads incorrectly removed
    tp = min(host_removed, gt_host)
    fp = max(0, host_removed - gt_host)
    fn = max(0, gt_host - host_removed)
    tn = gt_viral + gt_other - fp

    sensitivity = tp / max(tp + fn, 1)
    specificity = tn / max(tn + fp, 1)
    precision = tp / max(tp + fp, 1)
    fpr = fp / max(fp + tn, 1)
    f1 = 2 * precision * sensitivity / max(precision + sensitivity, 1e-10)

    return {
        "gt_host": gt_host,
        "gt_viral": gt_viral,
        "host_removed": host_removed,
        "ambig_flagged": ambig_flagged,
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "fpr": fpr,
        "f1": f1,
    }


def main():
    print("=" * 70)
    print("Host Depletion Threshold Sweep")
    print("=" * 70)

    if not HOST_FILTER.exists():
        print(f"ERROR: Host filter not found at {HOST_FILTER}")
        print("Build with: virome-qc db --host benchmark_data/references/chm13v2.0.fa")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for seed in SEEDS:
        r1_path = generate_dataset(seed)
        if r1_path is None:
            continue

        r2_path = Path(str(r1_path).replace("_R1.", "_R2."))
        labels = read_source_labels(r1_path)
        n_host = sum(1 for s in labels.values() if s == "host_dna")
        n_viral = sum(1 for s in labels.values() if s == "viral")
        print(f"  Dataset seed={seed}: {len(labels):,} reads ({n_viral:,} viral, {n_host:,} host)")

        for host_thresh in HOST_THRESHOLDS:
            for ambig_thresh in AMBIG_THRESHOLDS:
                if ambig_thresh >= host_thresh:
                    continue  # ambiguous must be below host threshold

                tag = f"seed{seed}_h{int(host_thresh*100):02d}_a{int(ambig_thresh*100):02d}"
                qc_dir = OUTPUT_DIR / tag

                passport = run_virome_qc(r1_path, r2_path, qc_dir, host_thresh, ambig_thresh)
                if not passport:
                    continue

                metrics = compute_metrics(passport, labels)
                if not metrics:
                    continue

                result = {
                    "seed": seed,
                    "host_threshold": host_thresh,
                    "ambiguous_threshold": ambig_thresh,
                    **metrics,
                }
                all_results.append(result)

    # Save results
    results_path = OUTPUT_DIR / "host_sweep_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)

    # Summary: aggregate across seeds for each threshold pair
    print("\n" + "=" * 100)
    print("Host depletion ROC (mean across 3 replicates)")
    print(f"{'host_t':>7s} {'ambig_t':>8s} {'Removed':>8s} {'Ambig':>8s} "
          f"{'Sens':>8s} {'Spec':>8s} {'FPR':>10s} {'Prec':>8s} {'F1':>8s}")
    print("-" * 100)

    seen = set()
    for r in sorted(all_results, key=lambda x: (x["host_threshold"], x["ambiguous_threshold"])):
        key = (r["host_threshold"], r["ambiguous_threshold"])
        if key in seen:
            continue
        seen.add(key)

        subset = [x for x in all_results
                  if x["host_threshold"] == key[0] and x["ambiguous_threshold"] == key[1]]

        marker = " <-- current" if key == (0.50, 0.20) else ""
        print(
            f"{key[0]:>7.2f} {key[1]:>8.2f} "
            f"{np.mean([r['host_removed'] for r in subset]):>8.0f} "
            f"{np.mean([r['ambig_flagged'] for r in subset]):>8.0f} "
            f"{np.mean([r['sensitivity'] for r in subset]):>8.3f} "
            f"{np.mean([r['specificity'] for r in subset]):>8.4f} "
            f"{np.mean([r['fpr'] for r in subset]):>10.5f} "
            f"{np.mean([r['precision'] for r in subset]):>8.3f} "
            f"{np.mean([r['f1'] for r in subset]):>8.3f}"
            f"{marker}"
        )

    # Optimal by Youden's J
    print()
    best = max(all_results, key=lambda r: r["sensitivity"] + r["specificity"] - 1)
    print(f"Optimal by Youden's J: host={best['host_threshold']:.2f}, "
          f"ambig={best['ambiguous_threshold']:.2f} "
          f"(J={best['sensitivity'] + best['specificity'] - 1:.4f})")

    best_f1 = max(all_results, key=lambda r: r["f1"])
    print(f"Optimal by F1: host={best_f1['host_threshold']:.2f}, "
          f"ambig={best_f1['ambiguous_threshold']:.2f} "
          f"(F1={best_f1['f1']:.3f})")


if __name__ == "__main__":
    main()
