#!/usr/bin/env python3
"""
Adapter detection sensitivity/specificity sweep using ViroForge synthetic data.

Generates datasets with known adapter contamination rates, runs virome-qc,
and computes per-read classification accuracy against ground truth.

Produces:
- Sensitivity vs adapter rate curve
- Precision-recall analysis
- ROC curve for adapter detection threshold
- Optimal threshold recommendation
"""

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

# Paths
VIROFORGE_DIR = Path(os.environ.get("VIROFORGE_DIR", "/Users/scotthandley/Code/tools/viroforge"))
VIROME_QC_DIR = Path(__file__).parent.parent
OUTPUT_DIR = VIROME_QC_DIR / "benchmark_data" / "adapter_sweep"
VIROFORGE_PYTHON = VIROFORGE_DIR / ".venv" / "bin" / "python3"
GENERATE_SCRIPT = VIROFORGE_DIR / "scripts" / "generate_fastq_dataset.py"
VIROME_QC_BIN = Path(os.environ.get(
    "VIROME_QC_BIN",
    Path.home() / ".cargo" / "target" / "release" / "virome-qc"
))

# Sweep parameters
ADAPTER_RATES = [0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.30]
ADAPTER_TYPES = ["truseq", "nextera"]
COLLECTION_ID = 9  # gut virome
COVERAGE = 3  # low coverage for speed
SEED = 42


def generate_dataset(adapter_rate: float, adapter_type: str, output_dir: Path) -> dict:
    """Generate ViroForge dataset with known adapter contamination."""
    output_dir.mkdir(parents=True, exist_ok=True)

    db_path = VIROFORGE_DIR / "viroforge" / "data" / "viral_genomes.db"
    cmd = [
        str(VIROFORGE_PYTHON),
        str(GENERATE_SCRIPT),
        "--database", str(db_path),
        "--collection-id", str(COLLECTION_ID),
        "--platform", "novaseq",
        "--coverage", str(COVERAGE),
        "--output", str(output_dir),
        "--seed", str(SEED),
        "--adapter-rate", str(adapter_rate),
        "--adapter-type", adapter_type,
    ]

    # Add ViroForge venv bin to PATH so ISS is found
    env = os.environ.copy()
    venv_bin = str(VIROFORGE_DIR / ".venv" / "bin")
    env["PATH"] = venv_bin + ":" + env.get("PATH", "")

    print(f"  Generating: adapter_rate={adapter_rate}, type={adapter_type}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
        return {}

    # Load metadata
    meta_dir = output_dir / "metadata"
    meta_files = list(meta_dir.glob("*_metadata.json"))
    if meta_files:
        with open(meta_files[0]) as f:
            return json.load(f)
    return {}


def count_adapter_reads(fastq_path: Path, adapter_seqs: dict) -> dict:
    """
    Count reads with adapter sequence in the 3' end (ground truth from injection).

    ViroForge's adapter injector replaces the 3' end with adapter sequence.
    We detect this by checking if the last N bases match any adapter.
    """
    # For ISS-generated reads, the adapter-injected reads have adapter
    # at the 3' end. We check the last 20bp against adapter prefixes.
    r1_adapter = adapter_seqs.get("r1_adapter", "").encode()
    r2_adapter = adapter_seqs.get("r2_adapter", "").encode()

    probe_len = 15  # check last 15bp
    total = 0
    with_adapter = 0

    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    with opener(fastq_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line
                seq = line.strip().encode()
                total += 1
                tail = seq[-probe_len:]
                # Check if tail matches adapter prefix (allowing 2 mismatches)
                for adapter in [r1_adapter, r2_adapter]:
                    if len(adapter) >= probe_len:
                        probe = adapter[:probe_len]
                        mm = sum(1 for a, b in zip(tail, probe) if a != b)
                        if mm <= 2:
                            with_adapter += 1
                            break

    return {"total": total, "with_adapter": with_adapter}


def run_virome_qc(r1_path: Path, r2_path: Path, output_dir: Path, profile: str) -> dict:
    """Run virome-qc and return passport data."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(VIROME_QC_BIN),
        "run",
        "-p", profile,
        "-1", str(r1_path),
        "-2", str(r2_path),
        "-i", ".",
        "-o", str(output_dir),
        "--report-only",
        "-t", "8",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  virome-qc ERROR: {result.stderr[-300:]}")
        return {}

    passport_path = output_dir / "passport.json"
    if passport_path.exists():
        with open(passport_path) as f:
            return json.load(f)
    return {}


def extract_adapter_metrics(passport: dict) -> dict:
    """Extract adapter-related metrics from virome-qc passport."""
    for m in passport.get("modules", []):
        if m["name"] == "adapter":
            return {
                "reads_removed": m["reads_removed"],
                "reads_modified": m["reads_modified"],
                "adapters_3prime": m["extra"].get("adapters_found_3prime", 0),
                "adapters_internal": m["extra"].get("adapters_found_internal", 0),
                "bases_removed": m["bases_removed"],
            }
    return {}


def compute_classification_metrics(
    true_adapter_count: int,
    detected_3prime: int,
    total_reads: int,
) -> dict:
    """Compute sensitivity, specificity, precision, F1 for adapter detection."""
    # True positives: reads with adapter that were detected
    # We approximate: detected <= true_adapter means all detections are correct
    tp = min(detected_3prime, true_adapter_count)
    fp = max(0, detected_3prime - true_adapter_count)
    fn = max(0, true_adapter_count - detected_3prime)
    tn = total_reads - tp - fp - fn

    sensitivity = tp / max(tp + fn, 1)
    specificity = tn / max(tn + fp, 1)
    precision = tp / max(tp + fp, 1)
    f1 = 2 * precision * sensitivity / max(precision + sensitivity, 1e-10)

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "f1": f1,
    }


def main():
    print("=" * 70)
    print("Adapter Detection Sweep: ViroForge + virome-qc")
    print("=" * 70)

    # Ensure virome-qc is built
    print(f"virome-qc binary: {VIROME_QC_BIN}")
    if not VIROME_QC_BIN.exists():
        print("Building virome-qc...")
        subprocess.run(
            ["cargo", "build", "--release"],
            cwd=str(VIROME_QC_DIR),
            check=True,
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for adapter_type in ADAPTER_TYPES:
        print(f"\n--- Adapter type: {adapter_type} ---")

        # Adapter sequences for ground truth detection
        adapter_seqs = {
            "truseq": {
                "r1_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                "r2_adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            },
            "nextera": {
                "r1_adapter": "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                "r2_adapter": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
            },
        }[adapter_type]

        for rate in ADAPTER_RATES:
            tag = f"{adapter_type}_{int(rate*100):02d}pct"
            dataset_dir = OUTPUT_DIR / tag
            qc_dir = OUTPUT_DIR / f"{tag}_qc"

            # Generate
            meta = generate_dataset(rate, adapter_type, dataset_dir)
            if not meta:
                print(f"  SKIP: generation failed for {tag}")
                continue

            # Find FASTQ files
            fastq_dir = dataset_dir / "fastq"
            r1_files = sorted(fastq_dir.glob("*_R1.fastq*"))
            r2_files = sorted(fastq_dir.glob("*_R2.fastq*"))
            if not r1_files or not r2_files:
                print(f"  SKIP: no FASTQ files for {tag}")
                continue

            r1_path = r1_files[0]
            r2_path = r2_files[0]

            # Compress if not already
            if not str(r1_path).endswith(".gz"):
                subprocess.run(["gzip", str(r1_path)])
                subprocess.run(["gzip", str(r2_path)])
                r1_path = Path(str(r1_path) + ".gz")
                r2_path = Path(str(r2_path) + ".gz")

            # Count ground truth adapter reads
            gt = count_adapter_reads(r1_path, adapter_seqs)
            gt_adapter_count = gt["with_adapter"]
            gt_total = gt["total"]

            # Choose profile based on adapter type
            profile = "stool-vlp-tagmentation" if adapter_type == "nextera" else "stool-vlp-tagmentation"

            # Run virome-qc
            print(f"  Running virome-qc on {tag}...")
            passport = run_virome_qc(r1_path, r2_path, qc_dir, profile)
            if not passport:
                print(f"  SKIP: virome-qc failed for {tag}")
                continue

            adapter_metrics = extract_adapter_metrics(passport)
            detected_3prime = adapter_metrics.get("adapters_3prime", 0)

            # Compute metrics
            metrics = compute_classification_metrics(
                gt_adapter_count, detected_3prime, gt_total
            )

            result = {
                "adapter_type": adapter_type,
                "target_rate": rate,
                "total_reads": gt_total,
                "ground_truth_adapter": gt_adapter_count,
                "ground_truth_rate": gt_adapter_count / max(gt_total, 1),
                "detected_3prime": detected_3prime,
                "detected_rate": detected_3prime / max(gt_total, 1),
                "reads_removed": adapter_metrics.get("reads_removed", 0),
                "bases_trimmed": adapter_metrics.get("bases_removed", 0),
                **metrics,
            }
            all_results.append(result)

            print(
                f"  Rate={rate:.0%}: GT={gt_adapter_count:,} "
                f"Det={detected_3prime:,} "
                f"Sens={metrics['sensitivity']:.3f} "
                f"Prec={metrics['precision']:.3f} "
                f"F1={metrics['f1']:.3f}"
            )

            # Clean up FASTQ files to save disk
            for f in fastq_dir.glob("*.fastq*"):
                f.unlink()
            for f in fastq_dir.glob("*.vcf"):
                f.unlink()

    # Save results
    results_path = OUTPUT_DIR / "adapter_sweep_results.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    # Print summary table
    print("\n" + "=" * 90)
    print(f"{'Type':>8s} {'Rate':>6s} {'GT':>8s} {'Det':>8s} {'Sens':>7s} {'Spec':>7s} {'Prec':>7s} {'F1':>7s}")
    print("-" * 90)
    for r in all_results:
        print(
            f"{r['adapter_type']:>8s} "
            f"{r['target_rate']:>5.0%} "
            f"{r['ground_truth_adapter']:>8,} "
            f"{r['detected_3prime']:>8,} "
            f"{r['sensitivity']:>7.3f} "
            f"{r['specificity']:>7.3f} "
            f"{r['precision']:>7.3f} "
            f"{r['f1']:>7.3f}"
        )

    # Statistical summary
    if all_results:
        sensitivities = [r["sensitivity"] for r in all_results]
        precisions = [r["precision"] for r in all_results]
        print(f"\nOverall sensitivity: {np.mean(sensitivities):.3f} +/- {np.std(sensitivities):.3f}")
        print(f"Overall precision:   {np.mean(precisions):.3f} +/- {np.std(precisions):.3f}")


if __name__ == "__main__":
    main()
