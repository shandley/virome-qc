#!/bin/bash
# HTCF Setup Script for virome-qc benchmarking
# Run on login.htcf.wustl.edu as shandley
set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
DBDIR="$WORKDIR/databases"
SRCDIR="$WORKDIR/src"

echo "=== virome-qc HTCF Setup ==="
mkdir -p "$WORKDIR" "$DBDIR" "$SRCDIR"

# 1. Set up Rust
echo "--- Setting up Rust ---"
source ~/.cargo/env 2>/dev/null || true
rustup default stable 2>/dev/null || {
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source ~/.cargo/env
}
rustc --version
cargo --version

# 2. Clone repositories
echo "--- Cloning repositories ---"
cd "$SRCDIR"

if [ ! -d "biometal" ]; then
    git clone https://github.com/shandley/biometal.git
fi

if [ ! -d "SuperBloom" ]; then
    git clone https://github.com/EtienneC-K/SuperBloom.git
fi

if [ ! -d "virome-qc" ]; then
    git clone https://github.com/shandley/virome-qc.git
fi

# Update all repos
cd "$SRCDIR/biometal" && git pull
cd "$SRCDIR/SuperBloom" && git pull
cd "$SRCDIR/virome-qc" && git pull

# 3. Patch Cargo.toml for cluster paths
echo "--- Patching Cargo.toml for cluster paths ---"
cd "$SRCDIR/virome-qc"
sed -i "s|path = \"/Users/scotthandley/Code/tools/biometal\"|path = \"$SRCDIR/biometal\"|" Cargo.toml
sed -i "s|path = \"/Users/scotthandley/Code/tools/SuperBloom\"|path = \"$SRCDIR/SuperBloom\"|" Cargo.toml

# 4. Build virome-qc
echo "--- Building virome-qc (release) ---"
cargo build --release 2>&1 | tail -5
VIROME_QC="$SRCDIR/virome-qc/target/release/virome-qc"
echo "Binary: $VIROME_QC"
$VIROME_QC --version

# 5. Set up conda environment with sra-tools
echo "--- Setting up conda environment ---"
source /ref/sahlab/software/miniforge3/bin/activate
if ! conda env list | grep -q virome-qc-bench; then
    mamba create -y -n virome-qc-bench -c bioconda -c conda-forge sra-tools
fi
conda activate virome-qc-bench
fastq-dump --version | head -1 || echo "fastq-dump not available"
fasterq-dump --version | head -1 || echo "fasterq-dump not available"

# 6. Download T2T-CHM13 and build host filter
echo "--- Building host filter ---"
if [ ! -f "$DBDIR/human_t2t.sbf" ]; then
    cd "$DBDIR"
    echo "Downloading T2T-CHM13..."
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
        -O t2t_chm13.fna.gz
    gunzip t2t_chm13.fna.gz
    echo "Building Super Bloom filter..."
    $VIROME_QC db --host t2t_chm13.fna -o "$DBDIR/human_t2t.sbf"
    rm t2t_chm13.fna
    echo "Host filter: $(ls -lh $DBDIR/human_t2t.sbf)"
else
    echo "Host filter exists: $(ls -lh $DBDIR/human_t2t.sbf)"
fi

# 7. Build SILVA rRNA filter
echo "--- Building rRNA filter ---"
if [ ! -f "$DBDIR/rrna_silva.rrf" ]; then
    cd "$DBDIR"
    echo "Downloading SILVA..."
    $VIROME_QC db --rrna -o "$DBDIR/rrna_silva.rrf"
    echo "rRNA filter: $(ls -lh $DBDIR/rrna_silva.rrf)"
else
    echo "rRNA filter exists: $(ls -lh $DBDIR/rrna_silva.rrf)"
fi

# 8. Create symlinks for filter auto-discovery
echo "--- Creating filter symlinks ---"
cd "$SRCDIR/virome-qc"
mkdir -p benchmark_data/references
ln -sf "$DBDIR/human_t2t.sbf" benchmark_data/references/human_t2t.sbf
ln -sf "$DBDIR/rrna_silva.rrf" benchmark_data/references/rrna_silva.rrf

echo ""
echo "=== Setup complete ==="
echo "virome-qc binary: $VIROME_QC"
echo "Host filter: $DBDIR/human_t2t.sbf"
echo "rRNA filter: $DBDIR/rrna_silva.rrf"
echo "Working directory: $WORKDIR"
