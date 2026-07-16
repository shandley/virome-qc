#!/bin/bash
#SBATCH --job-name=vf-ref
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --output=/scratch/sahlab/shandley/virome-qc-benchmark/logs/vf_ref_%A_%a.out
#SBATCH --error=/scratch/sahlab/shandley/virome-qc-benchmark/logs/vf_ref_%A_%a.err
#
# Generate ViroForge reference datasets for virome-qc expected-range calibration.
#
# Submit: sbatch --array=1-4 htcf/generate_references.sh
#   1 = stool-vlp-tagmentation (5 datasets)
#   2 = tissue-truseq (5 datasets)
#   3 = metagenomics-nextera (5 datasets)
#   4 = low-biomass-wga (5 datasets)
#
# Each array task generates 5 datasets for one profile (~20-30 min each).

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
VFDIR="$WORKDIR/src/viroforge"
REFDIR="$WORKDIR/viroforge_references"
CONFIGDIR="$VFDIR/examples/batch_configs/virome_qc_profiles"
ENVBIN="/ref/sahlab/software/scott_conda/viroforge-bench/bin"
PYTHON="$ENVBIN/python3"

# Put conda env on PATH so ISS subprocess calls work
export PATH="$ENVBIN:$PATH"

mkdir -p "$REFDIR" "$WORKDIR/logs"

# Profile configs indexed by SLURM array task ID
PROFILES=(
    "stool_vlp_tagmentation"
    "tissue_truseq"
    "metagenomics_nextera"
    "low_biomass_wga"
)

IDX=${SLURM_ARRAY_TASK_ID:-${1:-1}}
PROFILE="${PROFILES[$((IDX-1))]}"
CONFIG="$CONFIGDIR/${PROFILE}.yaml"

echo "=== Generating ViroForge references: $PROFILE ==="
echo "Config: $CONFIG"
echo "Output: $REFDIR/$PROFILE/"
echo "Python: $PYTHON"
echo "ISS: $(which iss)"
echo "Time: $(date)"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config not found: $CONFIG"
    exit 1
fi

# Generate each dataset individually from the batch config
# (avoids the batch CLI's interactive prompts in SLURM)
cd "$VFDIR"
$PYTHON -c "
import yaml, sys, subprocess, os

config = yaml.safe_load(open('$CONFIG'))
datasets = config.get('datasets', [])
outbase = '$REFDIR/$PROFILE'
os.makedirs(outbase, exist_ok=True)

for i, ds in enumerate(datasets):
    name = ds['name']
    outdir = f'{outbase}/{name}'
    print(f'\\n=== Dataset {i+1}/{len(datasets)}: {name} ===')

    cmd = [
        '$PYTHON', 'scripts/generate_fastq_dataset.py',
        '--collection-id', str(ds['collection_id']),
        '--output', outdir,
        '--platform', ds.get('platform', 'novaseq'),
        '--coverage', '10',  # 10x is sufficient for calibration (~1M reads)
        '--seed', str(ds.get('seed', 42)),
    ]

    # VLP
    if ds.get('no_vlp'):
        cmd.append('--no-vlp')
    elif 'vlp_protocol' in ds:
        cmd.extend(['--vlp-protocol', ds['vlp_protocol']])

    # Contamination
    if 'contamination_level' in ds:
        cmd.extend(['--contamination-level', ds['contamination_level']])

    # Adapter / insert size
    if 'adapter_type' in ds:
        cmd.extend(['--adapter-type', ds['adapter_type']])
    if 'mean_insert_size' in ds:
        cmd.extend(['--mean-insert-size', str(ds['mean_insert_size'])])
    if 'insert_size_sd' in ds:
        cmd.extend(['--insert-size-sd', str(ds['insert_size_sd'])])
    if 'chimera_rate' in ds:
        cmd.extend(['--chimera-rate', str(ds['chimera_rate'])])

    # Duplicates
    if 'duplicate_rate' in ds:
        cmd.extend(['--duplicate-rate', str(ds['duplicate_rate'])])
    if 'duplicate_error_rate' in ds:
        cmd.extend(['--duplicate-error-rate', str(ds['duplicate_error_rate'])])
    if 'duplicate_max_copies' in ds:
        cmd.extend(['--duplicate-max-copies', str(ds['duplicate_max_copies'])])

    # Low complexity
    if 'low_complexity_rate' in ds:
        cmd.extend(['--low-complexity-rate', str(ds['low_complexity_rate'])])

    print(f'CMD: {\" \".join(cmd)}')
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        print(f'FAILED: {name}')
        sys.exit(1)
    print(f'OK: {name}')

print(f'\\n=== All {len(datasets)} datasets generated ===')
" 2>&1

echo ""
echo "--- Generation complete ---"
echo "Output: $REFDIR/$PROFILE/"
ls -la "$REFDIR/$PROFILE/" 2>/dev/null || echo "(no output yet)"
echo "Time: $(date)"
