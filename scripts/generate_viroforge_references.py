#!/usr/bin/env python3
"""
Generate ViroForge reference datasets for virome-qc expected range calibration.

Each dataset is generated with known contamination levels and processed through
virome-qc to establish expected QC metric ranges per sample type.
"""

import json
import subprocess
import sys
from pathlib import Path

VIROFORGE_DIR = Path("/Users/scotthandley/Code/tools/viroforge")
VIROME_QC_BIN = Path.home() / ".cargo" / "target" / "release" / "virome-qc"
OUTPUT_BASE = Path(__file__).parent.parent / "benchmark_data" / "viroforge_reference"

# Dataset definitions
# (name, collection_id, contamination_level, platform, coverage, vlp_protocol, molecule_type, virome_qc_profile, extra_args)
DATASETS = [
    # 1. Gut VLP clean
    ("01_gut_vlp_clean", 9, "clean", "novaseq", 10, "tangential_flow", "dna",
     "stool-vlp-tagmentation", []),
    # 2. Gut VLP realistic (already generated as 01_gut_vlp_realistic)
    ("02_gut_vlp_realistic", 9, "realistic", "novaseq", 10, "tangential_flow", "dna",
     "stool-vlp-tagmentation", []),
    # 3. Gut bulk (no VLP)
    ("03_gut_bulk_no_vlp", 9, "heavy", "novaseq", 10, None, "dna",
     "profiles/short-read-nebnext.yaml", []),
    # 4. IBD gut (disease + inflammation)
    ("04_gut_ibd", 18, "heavy", "novaseq", 10, "tangential_flow", "dna",
     "stool-vlp-tagmentation", []),
    # 5. Oral VLP
    ("05_oral_vlp", 10, "realistic", "novaseq", 10, "tangential_flow", "dna",
     "stool-vlp-tagmentation", []),
    # 6. Respiratory RNA virome
    ("06_respiratory_rna", 21, "realistic", "novaseq", 10, None, "rna",
     "profiles/short-read-nebnext.yaml", []),
    # 7. Fecal RNA virome
    ("07_fecal_rna", 23, "realistic", "novaseq", 10, None, "rna",
     "profiles/short-read-nebnext.yaml", []),
    # 8. Marine virome
    ("08_marine", 13, "clean", "miseq", 10, "tangential_flow", "dna",
     "profiles/short-read-nebnext.yaml", []),
    # 9. Soil virome
    ("09_soil", 14, "clean", "novaseq", 10, "tangential_flow", "dna",
     "profiles/short-read-nebnext.yaml", []),
    # 10. Wastewater virome
    ("10_wastewater", 17, "realistic", "novaseq", 10, "tangential_flow", "dna",
     "profiles/short-read-nebnext.yaml", []),
    # 11. Blood/plasma (high host)
    ("11_blood_plasma", 25, "heavy", "novaseq", 10, None, "dna",
     "profiles/short-read-nebnext.yaml", []),
    # 12. Gut VLP failed (worst case)
    ("12_gut_vlp_failed", 9, "failed", "novaseq", 10, "tangential_flow", "dna",
     "stool-vlp-tagmentation", []),
]


def generate_dataset(name: str, collection_id: int, contamination: str,
                     platform: str, coverage: int, vlp_protocol: str | None,
                     molecule_type: str, extra_args: list[str]) -> bool:
    """Generate a ViroForge dataset."""
    output_dir = OUTPUT_BASE / name

    if (output_dir / "metadata").exists():
        # Check if FASTQ files exist
        fastq_dir = output_dir / "fastq"
        if fastq_dir.exists() and any(fastq_dir.glob("*.fastq")):
            print(f"  SKIP: {name} (already generated)")
            return True

    cmd = [
        "uv", "run", "viroforge", "generate",
        "--collection-id", str(collection_id),
        "--platform", platform,
        "--coverage", str(coverage),
        "--seed", "42",
        "--output", str(output_dir),
        "-v",
    ]
    cmd.extend(extra_args)

    print(f"  Generating {name} (collection={collection_id}, contam={contamination}, vlp={vlp_protocol})...")
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=600,
        cwd=str(VIROFORGE_DIR),
        env={**__import__("os").environ,
             "VIROFORGE_CONTAMINATION": contamination,
             "VIROFORGE_VLP": vlp_protocol or "none",
             "VIROFORGE_MOLECULE": molecule_type},
    )

    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-500:]}")
        return False

    return True


def run_virome_qc(name: str, profile: str) -> bool:
    """Run virome-qc on a generated dataset."""
    dataset_dir = OUTPUT_BASE / name
    fastq_dir = dataset_dir / "fastq"
    qc_output = dataset_dir / "qc_output"

    if (qc_output / "passport.json").exists():
        print(f"  SKIP QC: {name} (already processed)")
        return True

    # Find FASTQ files
    r1_files = sorted(fastq_dir.glob("*_R1.fastq"))
    r2_files = sorted(fastq_dir.glob("*_R2.fastq"))

    if not r1_files:
        print(f"  ERROR: No R1 FASTQ found in {fastq_dir}")
        return False

    cmd = [str(VIROME_QC_BIN), "run", "-p", profile]
    if r2_files:
        cmd.extend(["-1", str(r1_files[0]), "-2", str(r2_files[0]), "-i", ".", "--merge"])
    else:
        cmd.extend(["-i", str(r1_files[0])])

    cmd.extend(["-o", str(qc_output), "--report-only", "-t", "8"])

    print(f"  Running virome-qc on {name}...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[-500:]}")
        return False

    # Print summary
    for line in result.stderr.split("\n"):
        if "QC complete" in line or "rRNA" in line or "ERV" in line:
            print(f"    {line.strip()}")

    return (qc_output / "passport.json").exists()


def summarize_results():
    """Print summary table of all processed datasets."""
    print("\n" + "=" * 100)
    print("REFERENCE DATASET SUMMARY")
    print("=" * 100)
    print(f"{'Name':30s} {'Reads':>10s} {'Surv':>7s} {'Host%':>8s} {'rRNA%':>8s} {'Contam%':>8s} {'ERV':>6s} {'GC':>6s}")
    print("-" * 100)

    for name, _, _, _, _, _, _, profile, _ in DATASETS:
        passport_path = OUTPUT_BASE / name / "qc_output" / "passport.json"
        if not passport_path.exists():
            print(f"{'  ' + name:30s} (not processed)")
            continue

        d = json.load(open(passport_path))
        ri = d.get("reads_input", 1)
        surv = d.get("survival_rate", 0)

        host = sum(m.get("reads_removed", 0) for m in d.get("modules", []) if m.get("name") == "host")
        rrna = sum(m.get("reads_removed", 0) for m in d.get("modules", []) if m.get("name") == "rrna")
        contam = sum(m.get("reads_removed", 0) for m in d.get("modules", []) if m.get("name") == "contaminant")
        erv = d.get("erv_analysis", {}).get("retroviral_reads_flagged", 0)
        gc = d.get("ingestion", {}).get("mean_gc", 0)

        print(
            f"{'  ' + name:30s} {ri:>10,} {surv*100:>6.1f}% "
            f"{host/ri*100:>7.3f}% {rrna/ri*100:>7.3f}% {contam/ri*100:>7.3f}% "
            f"{erv:>6,} {gc:>5.2f}"
        )


def main():
    print("=" * 70)
    print("Generating ViroForge reference datasets")
    print("=" * 70)

    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

    # For now, just summarize what needs to be done
    # ViroForge CLI doesn't support --contamination-level directly in all cases
    # We need to handle this per-dataset
    print("\nNote: ViroForge preset overrides may not support all contamination levels.")
    print("Datasets requiring custom contamination need manual ViroForge configuration.\n")

    successful = []
    failed = []

    for name, coll_id, contam, platform, coverage, vlp, molecule, profile, extra in DATASETS:
        print(f"\n--- {name} ---")

        # Generate
        ok = generate_dataset(name, coll_id, contam, platform, coverage, vlp, molecule, extra)
        if not ok:
            failed.append(name)
            continue

        # Run QC
        ok = run_virome_qc(name, profile)
        if ok:
            successful.append(name)
        else:
            failed.append(name)

    summarize_results()

    print(f"\nResults: {len(successful)} successful, {len(failed)} failed")
    if failed:
        print(f"Failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
