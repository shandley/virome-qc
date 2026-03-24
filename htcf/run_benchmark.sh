#!/bin/bash
#SBATCH --job-name=vqc-bench
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/sahlab/shandley/virome-qc-benchmark/logs/vqc_%A_%a.out
#SBATCH --error=/scratch/sahlab/shandley/virome-qc-benchmark/logs/vqc_%A_%a.err
#
# virome-qc benchmark runner for HTCF
# Submit as array job: sbatch --array=1-15 run_benchmark.sh
#
# Or run individual sample: sbatch run_benchmark.sh <sample_index>

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
SRCDIR="$WORKDIR/src/virome-qc"
VIROME_QC="$SRCDIR/target/release/virome-qc"
DBDIR="$WORKDIR/databases"
RESULTS="$WORKDIR/results"
TMPDIR_BASE="$WORKDIR/tmp"

# Activate environments (don't use set -e here, conda activate can return non-zero)
source ~/.cargo/env 2>/dev/null || true
eval "$(conda shell.bash hook)" 2>/dev/null || source /ref/sahlab/software/miniforge3/bin/activate 2>/dev/null || true
conda activate virome-qc-bench 2>/dev/null || source /ref/sahlab/software/miniforge3/bin/activate virome-qc-bench 2>/dev/null || true
export PATH="/ref/sahlab/software/scott_conda/virome-qc-bench/bin:$PATH"

# Now enable strict mode
set -e

mkdir -p "$RESULTS" "$TMPDIR_BASE" logs

# Dataset definitions: index name accession paired profile extra_args
# Format: INDEX|NAME|ACCESSION|PAIRED|PROFILE|EXTRA_ARGS
DATASETS=(
    "1|cook_mock|ERR10359658|PE|stool-vlp-tagmentation|"
    "2|shkoporov_gut|SRR9161520|PE|stool-vlp-tagmentation|"
    "3|santos_kapa|SRR8487022|PE|stool-vlp-tagmentation|"
    "4|santos_nextera|SRR8487034|PE|stool-vlp-tagmentation|"
    "5|buddle_wgs|ERR13480651|SE|profiles/short-read-nebnext.yaml|"
    "6|buddle_rna|ERR13480663|SE|profiles/short-read-nebnext.yaml|"
    "7|zhang_undepleted|SRR33419012|PE|profiles/short-read-nebnext.yaml|"
    "8|zhang_depleted|SRR33419066|PE|profiles/short-read-nebnext.yaml|"
    "9|tara_ocean|ERR599370|PE|profiles/short-read-nebnext.yaml|"
    "10|chrisman_dnbseq|ERR9765742|SE|profiles/short-read-nebnext.yaml|"
    "11|tisza_wastewater|SRR24403399|PE|profiles/short-read-nebnext.yaml|"
    "12|cook_ont|ERR10359840|SE|profiles/short-read-nebnext.yaml|"
    "13|negativeome_blank|ERR3487542|PE|stool-vlp-tagmentation|"
    "14|hiv_gut_virome|SRR27430060|PE|profiles/short-read-nebnext.yaml|"
    "15|kleiner_mock|ERR1877475|PE|profiles/short-read-nebnext.yaml|"
)

# Get sample index (from array task or command line)
IDX=${SLURM_ARRAY_TASK_ID:-${1:-1}}
DATASET="${DATASETS[$((IDX-1))]}"

IFS='|' read -r _ NAME ACCESSION PAIRED PROFILE EXTRA <<< "$DATASET"

echo "=== Processing: $NAME ($ACCESSION) ==="
echo "Index: $IDX, Paired: $PAIRED, Profile: $PROFILE"
echo "Time: $(date)"

# Skip TBD datasets
if [ "$ACCESSION" = "TBD" ]; then
    echo "SKIP: accession TBD for $NAME"
    exit 0
fi

# Check if already processed
if [ -f "$RESULTS/$NAME/passport.json" ]; then
    echo "SKIP: $NAME already has passport.json"
    exit 0
fi

# Create temp directory for this sample
TMPDIR="$TMPDIR_BASE/$NAME"
mkdir -p "$TMPDIR"

# Download
echo "--- Downloading $ACCESSION ---"
cd "$TMPDIR"

if [ "$PAIRED" = "PE" ]; then
    fasterq-dump "$ACCESSION" --split-files --threads 4 --temp "$TMPDIR" 2>&1 | tail -5
    # Compress
    pigz -p 4 "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq" 2>/dev/null || \
    gzip "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
    R1="$TMPDIR/${ACCESSION}_1.fastq.gz"
    R2="$TMPDIR/${ACCESSION}_2.fastq.gz"
else
    fasterq-dump "$ACCESSION" --threads 4 --temp "$TMPDIR" 2>&1 | tail -5
    # fasterq-dump may produce _1 even for SE
    if [ -f "${ACCESSION}_1.fastq" ]; then
        pigz -p 4 "${ACCESSION}_1.fastq" 2>/dev/null || gzip "${ACCESSION}_1.fastq"
        R1="$TMPDIR/${ACCESSION}_1.fastq.gz"
    else
        pigz -p 4 "${ACCESSION}.fastq" 2>/dev/null || gzip "${ACCESSION}.fastq"
        R1="$TMPDIR/${ACCESSION}.fastq.gz"
    fi
    R2=""
fi

echo "Downloaded: $(ls -lh $TMPDIR/*.fastq.gz 2>/dev/null)"

# Run virome-qc
echo "--- Running virome-qc ---"
cd "$SRCDIR"
OUTDIR="$RESULTS/$NAME"
mkdir -p "$OUTDIR"

if [ "$PAIRED" = "PE" ] && [ -n "$R2" ] && [ -f "$R2" ]; then
    $VIROME_QC run -p "$PROFILE" -1 "$R1" -2 "$R2" -i . --merge \
        -o "$OUTDIR" --report-only -t 8 2>&1
else
    $VIROME_QC run -p "$PROFILE" -i "$R1" \
        -o "$OUTDIR" --report-only -t 8 2>&1
fi

# Verify
if [ -f "$OUTDIR/passport.json" ]; then
    echo "OK: passport.json created"
    # Print summary
    python3 -c "
import json
d = json.load(open('$OUTDIR/passport.json'))
ri = d['reads_input']
surv = d['survival_rate']
host = sum(m.get('reads_removed',0) for m in d['modules'] if m['name']=='host')
rrna = sum(m.get('reads_removed',0) for m in d['modules'] if m['name']=='rrna')
erv = d.get('erv_analysis',{}).get('retroviral_reads_flagged',0)
print(f'  Reads: {ri:,}  Survival: {surv*100:.1f}%  Host: {host/ri*100:.3f}%  rRNA: {rrna/ri*100:.3f}%  ERV: {erv}')
" 2>/dev/null
else
    echo "FAILED: no passport.json"
fi

# Clean up downloaded files
echo "--- Cleaning up ---"
rm -rf "$TMPDIR"

echo "=== Done: $NAME ($(date)) ==="
