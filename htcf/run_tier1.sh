#!/bin/bash
# Submit Tier 1 benchmark jobs to SLURM
# Run from: /scratch/sahlab/shandley/virome-qc-benchmark/src/virome-qc/htcf/

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
mkdir -p "$WORKDIR/logs"

echo "Submitting Tier 1 benchmark jobs..."

# Submit all 15 Tier 1 datasets as array job
sbatch --array=1-15 htcf/run_benchmark.sh

echo "Submitted. Monitor with: squeue -u shandley"
echo "Check results: ls $WORKDIR/results/*/passport.json"
echo "Logs: ls $WORKDIR/logs/"
