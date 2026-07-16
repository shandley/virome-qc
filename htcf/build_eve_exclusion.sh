#!/bin/bash
#SBATCH --job-name=eve-build
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/sahlab/shandley/virome-qc-benchmark/logs/eve_build_%j.out
#SBATCH --error=/scratch/sahlab/shandley/virome-qc-benchmark/logs/eve_build_%j.err
#
# Build EVE exclusion set for virome-qc host depletion.
#
# 1. Downloads high-identity viral references (ciHHV-6, HERV-K proviral)
# 2. Aligns them to T2T-CHM13 with minimap2
# 3. Extracts T2T subsequences at EVE coordinates
# 4. Generates Rust source with embedded sequences
#
# Run: sbatch htcf/build_eve_exclusion.sh

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
SRCDIR="$WORKDIR/src/virome-qc"
EVEDIR="$WORKDIR/eve_build"
DBDIR="$WORKDIR/databases"
ENVBIN="/ref/sahlab/software/scott_conda/viroforge-bench/bin"
MINIMAP2="/ref/sahlab/software/scott_conda/minimap2-2.24-h7132678_1/bin/minimap2"
SAMTOOLS="/ref/sahlab/software/scott_conda/hpp_gxd/bin/samtools"
export PATH="$ENVBIN:$(dirname $MINIMAP2):$(dirname $SAMTOOLS):$PATH"

mkdir -p "$EVEDIR" "$WORKDIR/logs"

echo "=== Building EVE Exclusion Set ($(date)) ==="

# Step 1: Get T2T-CHM13 reference
T2T="$EVEDIR/chm13v2.0.fa"
if [ ! -f "$T2T" ]; then
    echo "--- Downloading T2T-CHM13 ---"
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
        -O "$T2T.gz"
    gunzip "$T2T.gz"
    echo "T2T: $(ls -lh $T2T)"
else
    echo "T2T exists: $(ls -lh $T2T)"
fi

# Step 2: Curate high-identity viral references
# These are the viruses with known high-identity integrations in the human genome
# that can cause false host depletion of genuinely exogenous reads.
VIRAL_REFS="$EVEDIR/high_identity_eve_refs.fa"
echo "--- Downloading high-identity viral references ---"

"$ENVBIN/python3" - "$VIRAL_REFS" << 'PYEOF'
import sys
import urllib.request

outpath = sys.argv[1]

# Accessions for high-identity EVE source viruses
# These are the exogenous viruses whose reads could be falsely removed
# because the human genome contains integrated copies at >85% identity.
accessions = {
    # ciHHV-6: >99.9% identity integrations at chr17/chr18 telomeres (~0.5% of humans)
    "NC_001664.4": "HHV-6A",        # Human herpesvirus 6A (159 kb)
    "NC_000898.1": "HHV-6B",        # Human herpesvirus 6B (162 kb)
    # Recent HERV-K (HML-2): 85-95% identity, ~100 proviral loci
    # Use consensus/representative sequences from Dfam
    "AF164610.1": "HERV-K_K113",    # HERV-K113 provirus (chr19p12, ~30% of humans)
    "AY037928.1": "HERV-K_K106",    # HERV-K106 provirus (chr3p25.3, ~15%)
    # Known non-retroviral EVEs with high identity
    # Bornaviruses: EBLN-1 through EBLN-4 are endogenous but ancient (<70% identity)
    # - Not included here because they don't cause false host depletion at k=31
    # Parvoviruses/Filoviruses: similarly too diverged for k=31 cross-reactivity
}

print(f"Downloading {len(accessions)} viral reference sequences...")
sequences = []

for acc, name in accessions.items():
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={acc}&rettype=fasta&retmode=text"
    print(f"  {acc} ({name})...", end=" ", flush=True)
    try:
        response = urllib.request.urlopen(url)
        fasta = response.read().decode("utf-8")
        # Rename header to include our annotation
        lines = fasta.strip().split("\n")
        header = f">{name}_{acc} {lines[0][1:]}"
        seq = "".join(lines[1:])
        sequences.append((header, seq))
        print(f"{len(seq):,} bp")
    except Exception as e:
        print(f"FAILED: {e}")

with open(outpath, "w") as f:
    for header, seq in sequences:
        f.write(f"{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

total_bp = sum(len(s) for _, s in sequences)
print(f"\nWrote {len(sequences)} sequences ({total_bp:,} bp) to {outpath}")
PYEOF

echo "Viral refs: $(ls -lh $VIRAL_REFS)"

# Step 3: Verify tools
echo "minimap2: $(which minimap2) ($(minimap2 --version 2>&1))"
echo "samtools: $(which samtools)"
echo "pysam: $("$ENVBIN/python3" -c "import pysam; print('OK')" 2>&1)"

# Step 4: Index T2T if needed (for pysam FASTA extraction)
if [ ! -f "$T2T.fai" ]; then
    echo "--- Indexing T2T ---"
    samtools faidx "$T2T" 2>/dev/null || "$ENVBIN/python3" -c "import pysam; pysam.FastaFile('$T2T')"
fi

# Step 6: Run the build script
echo "--- Running build_eve_exclusion.py ---"
cd "$SRCDIR"
"$ENVBIN/python3" scripts/build_eve_exclusion.py \
    --t2t "$T2T" \
    --viral-refs "$VIRAL_REFS" \
    --output-bed "$EVEDIR/t2t_eve_regions.bed" \
    --output-fasta "$EVEDIR/eve_exclusion_sequences.fa" \
    --output-rust src/modules/eve_kmers.rs \
    --min-identity 85 \
    --min-length 100 \
    --threads 8 2>&1

# Step 7: Verify output
echo ""
echo "=== Results ==="
echo "BED: $(wc -l < $EVEDIR/t2t_eve_regions.bed) regions"
echo "FASTA: $(grep -c '^>' $EVEDIR/eve_exclusion_sequences.fa) sequences"
echo "Rust source: $(wc -l < src/modules/eve_kmers.rs) lines"
head -5 src/modules/eve_kmers.rs

echo ""
echo "=== Copy outputs ==="
# Copy BED and FASTA to references directory
mkdir -p references
cp "$EVEDIR/t2t_eve_regions.bed" references/
cp "$EVEDIR/eve_exclusion_sequences.fa" references/
echo "Copied to references/"

echo ""
echo "=== Done ($(date)) ==="
echo ""
echo "Next steps:"
echo "  1. Review src/modules/eve_kmers.rs (should have embedded EVE sequences)"
echo "  2. Rebuild virome-qc: cargo build --release"
echo "  3. Test with HHV-6B reads to verify EVE rescue"
