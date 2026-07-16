#!/bin/bash
# ViroForge setup on HTCF for expected-range calibration
# Run on login.htcf.wustl.edu as shandley
set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
VFDIR="$WORKDIR/src/viroforge"

echo "=== ViroForge HTCF Setup ==="

# 1. Set up conda environment with InSilicoSeq
echo "--- Setting up conda environment ---"
source /ref/sahlab/software/miniforge3/bin/activate
if ! conda env list | grep -q virome-qc-bench; then
    mamba create -y -n virome-qc-bench -c bioconda -c conda-forge \
        sra-tools insilicoseq python=3.12
else
    # Ensure iss is installed in existing env
    conda activate virome-qc-bench
    mamba install -y -c bioconda -c conda-forge insilicoseq 2>/dev/null || true
fi
conda activate virome-qc-bench

# 2. Clone/update ViroForge
echo "--- Setting up ViroForge ---"
mkdir -p "$WORKDIR/src"
cd "$WORKDIR/src"

if [ ! -d "viroforge" ]; then
    git clone https://github.com/shandley/viroforge.git
fi
cd viroforge && git pull

# 3. Install ViroForge in the conda env
echo "--- Installing ViroForge ---"
pip install -e . 2>&1 | tail -5

# 4. Verify
echo "--- Verifying installation ---"
viroforge --version
iss --version 2>&1 | head -1
python3 -c "from viroforge.simulators.duplicates import add_pcr_duplicates; print('  duplicates module: OK')"
python3 -c "from viroforge.simulators.adapters import add_insert_size_adapters; print('  adapters module: OK')"
python3 -c "from viroforge.cli.summary import run_summary; print('  summary module: OK')"

# 5. Copy batch configs
echo "--- Batch configs ---"
ls -la examples/batch_configs/virome_qc_profiles/

echo ""
echo "=== ViroForge setup complete ==="
echo "ViroForge: $VFDIR"
echo "Conda env: virome-qc-bench"
echo ""
echo "Next: sbatch htcf/generate_references.sh"
