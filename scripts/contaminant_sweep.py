#!/usr/bin/env python3
"""
Contaminant screening (PhiX/rRNA) parameter sweep using ViroForge.

Sweeps min_kmer_fraction threshold and computes ROC curves using
ViroForge source labels as exact ground truth.
"""

import gzip
import json
import os
import subprocess
from pathlib import Path
from collections import defaultdict

import numpy as np

VIROFORGE_DIR = Path(os.environ.get("VIROFORGE_DIR", "/Users/scotthandley/Code/tools/viroforge"))
VIROME_QC_DIR = Path(__file__).parent.parent
VIROME_QC_BIN = Path(os.environ.get(
    "VIROME_QC_BIN",
    Path.home() / ".cargo" / "target" / "release" / "virome-qc"
))
OUTPUT_DIR = VIROME_QC_DIR / "benchmark_data" / "contaminant_sweep"
VIROFORGE_PYTHON = VIROFORGE_DIR / ".venv" / "bin" / "python3"
GENERATE_SCRIPT = VIROFORGE_DIR / "scripts" / "generate_fastq_dataset.py"
DB_PATH = VIROFORGE_DIR / "viroforge" / "data" / "viral_genomes.db"

THRESHOLDS = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60]
SEEDS = [42, 43, 44]
COVERAGE = 3


def generate_dataset(seed: int) -> Path:
    """Generate ViroForge gut virome dataset."""
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


def run_virome_qc(r1_path: Path, r2_path: Path, output_dir: Path,
                  min_kmer_fraction: float) -> dict:
    """Run virome-qc with specific contaminant threshold."""
    output_dir.mkdir(parents=True, exist_ok=True)

    profile_path = output_dir / "profile.yaml"
    with open(profile_path, "w") as f:
        f.write(f'name: "contaminant_sweep_{min_kmer_fraction}"\n')
        f.write('description: "Contaminant threshold sweep"\n')
        f.write('platform: illumina\n')
        f.write('modules:\n')
        f.write('  adapter:\n    enabled: false\n    sequences: []\n    internal_scan: false\n')
        f.write('    random_primer_trim: 0\n    min_overlap: 8\n    max_mismatch_rate: 0.1\n')
        f.write('  quality:\n    enabled: false\n    window_size: 10\n    min_quality: 15\n')
        f.write('    min_mean_quality: 20.0\n    min_length: 30\n')
        f.write('  polyx:\n    enabled: false\n    platform_aware: false\n    min_length: 8\n')
        f.write('  complexity:\n    enabled: false\n    min_entropy: 0.5\n')
        f.write(f'  contaminant:\n    enabled: true\n    screen_rrna: true\n    screen_phix: true\n')
        f.write(f'    screen_vectors: true\n    screen_kitome: false\n    min_kmer_fraction: {min_kmer_fraction}\n')
        f.write('  host:\n    enabled: false\n    reference: human\n    eve_aware: false\n    rescue: false\n')
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


def compute_roc_metrics(passport: dict, source_labels: dict) -> dict:
    """Compute TP/FP/TN/FN from passport and source labels."""
    total_input = passport.get("reads_input", 0)
    total_passed = passport.get("reads_passed", 0)

    # Ground truth counts (R1 only, PE = x2)
    gt_viral = sum(1 for s in source_labels.values() if s == "viral") * 2
    gt_rrna = sum(1 for s in source_labels.values() if s == "rrna") * 2
    gt_phix = sum(1 for s in source_labels.values() if s == "phix") * 2
    gt_host = sum(1 for s in source_labels.values() if s == "host_dna") * 2
    gt_reagent = sum(1 for s in source_labels.values() if s == "reagent_bacteria") * 2
    gt_contam = gt_rrna + gt_phix + gt_host + gt_reagent

    # Detected counts from passport
    contam_mod = None
    for m in passport.get("modules", []):
        if m["name"] == "contaminant":
            contam_mod = m
            break

    if not contam_mod:
        return {}

    detected_phix = contam_mod["extra"].get("phix_removed", 0)
    detected_rrna = contam_mod["extra"].get("rrna_removed", 0)
    detected_vector = contam_mod["extra"].get("vector_removed", 0)
    detected_total = contam_mod["reads_removed"]

    # For ROC: contaminant reads are "positive", viral reads are "negative"
    # TP = contaminant reads correctly removed
    # FP = viral reads incorrectly removed
    # FN = contaminant reads not removed
    # TN = viral reads correctly kept

    # We know: total removed = TP + FP
    # We know: gt_contam = TP + FN (but host/reagent aren't screened by contaminant module)
    # Screened contaminants = rrna + phix (not host, not reagent)
    gt_screened = gt_rrna + gt_phix

    # Estimate TP: min of detected and ground truth for each category
    tp_phix = min(detected_phix, gt_phix)
    tp_rrna = min(detected_rrna, gt_rrna)
    tp = tp_phix + tp_rrna

    fp = max(0, detected_total - tp)
    fn = max(0, gt_screened - tp)
    tn = gt_viral - fp

    sensitivity = tp / max(tp + fn, 1)
    specificity = tn / max(tn + fp, 1)
    precision = tp / max(tp + fp, 1)
    fpr = fp / max(fp + tn, 1)
    f1 = 2 * precision * sensitivity / max(precision + sensitivity, 1e-10)

    return {
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "fpr": fpr,
        "f1": f1,
        "gt_rrna": gt_rrna,
        "gt_phix": gt_phix,
        "gt_screened": gt_screened,
        "detected_rrna": detected_rrna,
        "detected_phix": detected_phix,
        "detected_total": detected_total,
        "phix_sensitivity": tp_phix / max(gt_phix, 1),
        "rrna_sensitivity": tp_rrna / max(gt_rrna, 1),
    }


def main():
    print("=" * 70)
    print("Contaminant Screening Parameter Sweep (min_kmer_fraction)")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for seed in SEEDS:
        r1_path = generate_dataset(seed)
        if r1_path is None:
            continue

        r2_path = Path(str(r1_path).replace("_R1.", "_R2."))
        labels = read_source_labels(r1_path)
        n_viral = sum(1 for s in labels.values() if s == "viral")
        n_rrna = sum(1 for s in labels.values() if s == "rrna")
        n_phix = sum(1 for s in labels.values() if s == "phix")
        print(f"  Dataset seed={seed}: {len(labels):,} reads "
              f"({n_viral:,} viral, {n_rrna:,} rRNA, {n_phix:,} PhiX)")

        for threshold in THRESHOLDS:
            qc_dir = OUTPUT_DIR / f"gut_seed{seed}_thresh{int(threshold*100):02d}"
            passport = run_virome_qc(r1_path, r2_path, qc_dir, threshold)
            if not passport:
                continue

            metrics = compute_roc_metrics(passport, labels)
            result = {
                "seed": seed,
                "min_kmer_fraction": threshold,
                **metrics,
            }
            all_results.append(result)

            print(
                f"    thresh={threshold:.2f}: "
                f"Sens={metrics['sensitivity']:.3f} "
                f"Spec={metrics['specificity']:.4f} "
                f"Prec={metrics['precision']:.3f} "
                f"F1={metrics['f1']:.3f} "
                f"FPR={metrics['fpr']:.5f} "
                f"PhiX={metrics['phix_sensitivity']:.3f} "
                f"rRNA={metrics['rrna_sensitivity']:.3f}"
            )

    # Save results
    results_path = OUTPUT_DIR / "contaminant_sweep_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)

    # ROC table (aggregated across seeds)
    print("\n" + "=" * 90)
    print("ROC Data (mean +/- SD across 3 replicates)")
    print(f"{'Threshold':>10s} {'Sens':>10s} {'Spec':>10s} {'FPR':>10s} "
          f"{'Prec':>10s} {'F1':>10s} {'PhiX':>10s} {'rRNA':>10s}")
    print("-" * 90)

    for threshold in THRESHOLDS:
        subset = [r for r in all_results if r["min_kmer_fraction"] == threshold]
        if not subset:
            continue
        sens = [r["sensitivity"] for r in subset]
        spec = [r["specificity"] for r in subset]
        fprs = [r["fpr"] for r in subset]
        prec = [r["precision"] for r in subset]
        f1s = [r["f1"] for r in subset]
        phix = [r["phix_sensitivity"] for r in subset]
        rrna = [r["rrna_sensitivity"] for r in subset]

        marker = " <-- current" if abs(threshold - 0.40) < 0.01 else ""
        print(
            f"{threshold:>10.2f} "
            f"{np.mean(sens):>6.3f}+/-{np.std(sens):.3f}"
            f"{np.mean(spec):>7.4f}+/-{np.std(spec):.4f}"
            f"{np.mean(fprs):>8.5f}+/-{np.std(fprs):.5f}"
            f"{np.mean(prec):>7.3f}+/-{np.std(prec):.3f}"
            f"{np.mean(f1s):>7.3f}+/-{np.std(f1s):.3f}"
            f"{np.mean(phix):>7.3f}+/-{np.std(phix):.3f}"
            f"{np.mean(rrna):>7.3f}+/-{np.std(rrna):.3f}"
            f"{marker}"
        )

    # Youden's J
    print("\nOptimal threshold (Youden's J = Sensitivity + Specificity - 1):")
    best_j = -1
    best_thresh = 0
    for threshold in THRESHOLDS:
        subset = [r for r in all_results if r["min_kmer_fraction"] == threshold]
        if not subset:
            continue
        j = np.mean([r["sensitivity"] + r["specificity"] - 1 for r in subset])
        if j > best_j:
            best_j = j
            best_thresh = threshold
    print(f"  Optimal: min_kmer_fraction={best_thresh:.2f} (J={best_j:.4f})")


if __name__ == "__main__":
    main()
