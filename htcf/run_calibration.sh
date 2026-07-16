#!/bin/bash
# Submit the full expected-range calibration pipeline to HTCF.
#
# Three stages:
#   Stage 1: Generate ViroForge reference datasets (4 jobs, ~30 min each)
#   Stage 2: Run virome-qc on each dataset (20 jobs, ~10 min each)
#   Stage 3: Derive ranges from passports (1 job, ~1 min, login node)
#
# Usage:
#   cd /scratch/sahlab/shandley/virome-qc-benchmark/src/virome-qc
#   bash htcf/run_calibration.sh

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
mkdir -p "$WORKDIR/logs"

echo "=== virome-qc Expected Range Calibration Pipeline ==="
echo ""

# Stage 1: Generate ViroForge references (4 profiles x 5 datasets)
echo "Stage 1: Submitting ViroForge reference generation (4 jobs)..."
JOB1=$(sbatch --array=1-4 --parsable htcf/generate_references.sh)
echo "  Job ID: $JOB1"

# Stage 2: Run virome-qc on generated datasets (depends on Stage 1)
echo "Stage 2: Submitting virome-qc processing (20 jobs, after Stage 1)..."
JOB2=$(sbatch --array=1-20 --dependency=afterok:${JOB1} --parsable htcf/run_qc_on_references.sh)
echo "  Job ID: $JOB2"

echo ""
echo "=== Submitted ==="
echo "Monitor: squeue -u shandley"
echo ""
echo "When complete, run on login node:"
echo "  bash htcf/derive_ranges.sh"
echo ""
echo "Expected timeline:"
echo "  Stage 1: ~30 min (ViroForge generation)"
echo "  Stage 2: ~30 min (virome-qc processing, needs SILVA + T2T filters)"
echo "  Stage 3: ~1 min (range derivation, run manually)"
