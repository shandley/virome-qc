#!/bin/bash
# Generate all ViroForge reference datasets and run virome-qc on each
set -e

VIROFORGE_DIR="/Users/scotthandley/Code/tools/viroforge"
VIROME_QC="/Users/scotthandley/.cargo/target/release/virome-qc"
OUTPUT_BASE="/Users/scotthandley/Desktop/hvp-bioinformatics/virome-qc/benchmark_data/viroforge_reference"
SCRIPT="$VIROFORGE_DIR/scripts/generate_fastq_dataset.py"

generate_and_qc() {
    local name="$1"
    local coll_id="$2"
    local contam="$3"
    local platform="$4"
    local vlp_args="$5"
    local molecule="$6"
    local qc_profile="$7"

    local outdir="$OUTPUT_BASE/$name"

    echo ""
    echo "=== $name ==="

    # Generate
    if [ -f "$outdir/qc_output/passport.json" ]; then
        echo "  SKIP: already complete"
        return 0
    fi

    echo "  Generating (collection=$coll_id, contam=$contam, platform=$platform)..."
    cd "$VIROFORGE_DIR"

    local gen_cmd="uv run python3 $SCRIPT --collection-id $coll_id --contamination-level $contam --platform $platform --coverage 10 --seed 42 --output $outdir"

    if [ "$vlp_args" = "none" ]; then
        gen_cmd="$gen_cmd --no-vlp"
    else
        gen_cmd="$gen_cmd --vlp-protocol $vlp_args"
    fi

    if [ "$molecule" = "rna" ]; then
        gen_cmd="$gen_cmd --molecule-type rna"
    fi

    eval $gen_cmd 2>&1 | grep -E "✓|ERROR|generation complete|Viral fraction|Contamination:"

    if [ $? -ne 0 ] && [ ! -d "$outdir/fastq" ]; then
        echo "  FAILED: generation error"
        return 1
    fi

    # Find FASTQ files
    local r1=$(ls "$outdir/fastq/"*_R1.fastq 2>/dev/null | head -1)
    local r2=$(ls "$outdir/fastq/"*_R2.fastq 2>/dev/null | head -1)

    if [ -z "$r1" ]; then
        echo "  FAILED: no R1 FASTQ found"
        return 1
    fi

    # Run virome-qc
    echo "  Running virome-qc..."
    cd /Users/scotthandley/Desktop/hvp-bioinformatics/virome-qc

    if [ -n "$r2" ]; then
        $VIROME_QC run -p "$qc_profile" -1 "$r1" -2 "$r2" -i . --merge -o "$outdir/qc_output" --report-only -t 8 2>&1 | grep -E "QC complete|ERV|Platform|Warning"
    else
        $VIROME_QC run -p "$qc_profile" -i "$r1" -o "$outdir/qc_output" --report-only -t 8 2>&1 | grep -E "QC complete|ERV|Platform|Warning"
    fi

    if [ -f "$outdir/qc_output/passport.json" ]; then
        echo "  OK"
    else
        echo "  FAILED: no passport generated"
        return 1
    fi

    # Clean up FASTQ files to save disk space (keep metadata + passport)
    # rm -f "$outdir/fastq/"*.fastq "$outdir/fasta/"*.fasta
}

echo "========================================"
echo "ViroForge Reference Dataset Generation"
echo "========================================"

mkdir -p "$OUTPUT_BASE"

# name, collection_id, contamination, platform, vlp_protocol, molecule, qc_profile
generate_and_qc "01_gut_vlp_clean"      9  clean     novaseq tangential_flow dna stool-vlp-tagmentation
generate_and_qc "02_gut_vlp_realistic"  9  realistic novaseq tangential_flow dna stool-vlp-tagmentation
generate_and_qc "03_gut_bulk_no_vlp"    9  heavy     novaseq none            dna profiles/short-read-nebnext.yaml
generate_and_qc "04_gut_ibd"            18 heavy     novaseq tangential_flow dna stool-vlp-tagmentation
generate_and_qc "05_oral_vlp"           10 realistic novaseq tangential_flow dna stool-vlp-tagmentation
generate_and_qc "06_respiratory_rna"    21 realistic novaseq none            rna profiles/short-read-nebnext.yaml
generate_and_qc "07_fecal_rna"          23 realistic novaseq none            rna profiles/short-read-nebnext.yaml
generate_and_qc "08_marine"             13 clean     miseq   tangential_flow dna profiles/short-read-nebnext.yaml
generate_and_qc "09_soil"               14 clean     novaseq tangential_flow dna profiles/short-read-nebnext.yaml
generate_and_qc "10_wastewater"         17 realistic novaseq tangential_flow dna profiles/short-read-nebnext.yaml
generate_and_qc "11_blood_plasma"       25 heavy     novaseq none            dna profiles/short-read-nebnext.yaml
generate_and_qc "12_gut_vlp_heavy"      9  heavy     novaseq tangential_flow dna stool-vlp-tagmentation

echo ""
echo "========================================"
echo "Generation complete. Running summary..."
echo "========================================"

# Summary
python3 -c "
import json, os
base = '$OUTPUT_BASE'
print(f'{\"Name\":30s} {\"Reads\":>10s} {\"Surv\":>7s} {\"Host%\":>8s} {\"rRNA%\":>8s} {\"ERV\":>6s}')
print('-' * 75)
for d in sorted(os.listdir(base)):
    p = os.path.join(base, d, 'qc_output', 'passport.json')
    if not os.path.isfile(p):
        print(f'{d:30s} (not processed)')
        continue
    data = json.load(open(p))
    ri = data.get('reads_input', 1)
    surv = data.get('survival_rate', 0)
    host = sum(m.get('reads_removed',0) for m in data.get('modules',[]) if m.get('name')=='host')
    rrna = sum(m.get('reads_removed',0) for m in data.get('modules',[]) if m.get('name')=='rrna')
    erv = data.get('erv_analysis',{}).get('retroviral_reads_flagged',0)
    print(f'{d:30s} {ri:>10,} {surv*100:>6.1f}% {host/ri*100:>7.3f}% {rrna/ri*100:>7.3f}% {erv:>6,}')
"
