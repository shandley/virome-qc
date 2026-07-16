#!/bin/bash
#SBATCH --job-name=vqc-ref
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/sahlab/shandley/virome-qc-benchmark/logs/vqc_ref_%A_%a.out
#SBATCH --error=/scratch/sahlab/shandley/virome-qc-benchmark/logs/vqc_ref_%A_%a.err
#
# Run virome-qc on ViroForge reference datasets to produce passports.
#
# Submit: sbatch --array=1-20 htcf/run_qc_on_references.sh
#   Array indices map to: 4 profiles x 5 datasets = 20 jobs
#
# Depends on: generate_references.sh having completed.

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
SRCDIR="$WORKDIR/src/virome-qc"
VIROME_QC="$SRCDIR/target/release/virome-qc"
DBDIR="$WORKDIR/databases"
REFDIR="$WORKDIR/viroforge_references"
RESULTS="$WORKDIR/viroforge_qc_results"

# Set up environment
source ~/.cargo/env 2>/dev/null || true
ENVBIN="/ref/sahlab/software/scott_conda/viroforge-bench/bin"
export PATH="$ENVBIN:$PATH"

mkdir -p "$RESULTS"

# Map array index to profile + dataset
# 1-5: stool_vlp_tagmentation, 6-10: tissue_truseq,
# 11-15: metagenomics_nextera, 16-20: low_biomass_wga
PROFILES=(
    "stool_vlp_tagmentation|stool-vlp-tagmentation"
    "tissue_truseq|tissue-truseq"
    "metagenomics_nextera|metagenomics-nextera"
    "low_biomass_wga|low-biomass-wga"
)

# Dataset names within each profile (must match batch config names)
DATASET_NAMES_0=("clean" "typical" "slightly_dirty" "high_duplication" "high_adapter")
DATASET_NAMES_1=("clean_tissue" "typical_tissue" "high_host_tissue" "degraded_tissue" "ffpe_like")
DATASET_NAMES_2=("clean_metag" "typical_metag" "marine_metag" "wastewater_metag" "dirty_metag")
DATASET_NAMES_3=("good_wga" "typical_wga" "high_dup_wga" "extreme_wga" "failed_wga")

IDX=${SLURM_ARRAY_TASK_ID:-${1:-1}}
PROFILE_IDX=$(( (IDX - 1) / 5 ))
DATASET_IDX=$(( (IDX - 1) % 5 ))

IFS='|' read -r PROFILE_DIR VQC_PROFILE <<< "${PROFILES[$PROFILE_IDX]}"

# Get dataset name from the appropriate array
eval "DATASET_NAME=\${DATASET_NAMES_${PROFILE_IDX}[$DATASET_IDX]}"

DATASET_DIR="$REFDIR/$PROFILE_DIR/$DATASET_NAME"
OUTDIR="$RESULTS/${PROFILE_DIR}/${DATASET_NAME}"

echo "=== Running virome-qc on ViroForge reference ==="
echo "Profile: $VQC_PROFILE (dir: $PROFILE_DIR)"
echo "Dataset: $DATASET_NAME"
echo "Input: $DATASET_DIR"
echo "Output: $OUTDIR"
echo "Time: $(date)"

# Find FASTQ files
R1=$(find "$DATASET_DIR" -name "*_R1.fastq*" -o -name "*_R1.fq*" | head -1)
R2=$(find "$DATASET_DIR" -name "*_R2.fastq*" -o -name "*_R2.fq*" | head -1)

if [ -z "$R1" ]; then
    echo "ERROR: No R1 FASTQ found in $DATASET_DIR"
    ls -lR "$DATASET_DIR" 2>/dev/null
    exit 1
fi

mkdir -p "$OUTDIR"

echo "R1: $R1"
echo "R2: ${R2:-none}"

# Run virome-qc
cd "$SRCDIR"
if [ -n "$R2" ] && [ -f "$R2" ]; then
    $VIROME_QC run -p "$VQC_PROFILE" -1 "$R1" -2 "$R2" -i . --merge \
        -o "$OUTDIR" --report-only -t 8 2>&1
else
    $VIROME_QC run -p "$VQC_PROFILE" -i "$R1" \
        -o "$OUTDIR" --report-only -t 8 2>&1
fi

# Verify
if [ -f "$OUTDIR/passport.json" ]; then
    echo "OK: passport.json created"
    python3 -c "
import json
d = json.load(open('$OUTDIR/passport.json'))
ri = d['reads_input']
surv = d['survival_rate']
host = sum(m.get('reads_removed',0) for m in d['modules'] if m['name']=='host')
rrna = sum(m.get('reads_removed',0) for m in d['modules'] if m['name']=='rrna')
dedup = sum(m.get('reads_removed',0) for m in d['modules'] if m['name']=='dedup')
print(f'  Reads: {ri:,}  Survival: {surv*100:.1f}%  Host: {host/ri*100:.3f}%  rRNA: {rrna/ri*100:.3f}%  Dedup: {dedup/ri*100:.1f}%')
" 2>/dev/null
else
    echo "FAILED: no passport.json"
    exit 1
fi

echo "=== Done: $PROFILE_DIR/$DATASET_NAME ($(date)) ==="
