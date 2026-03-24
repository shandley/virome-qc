#!/bin/bash
# Head-to-head comparison: virome-qc vs fastp on representative virome datasets
set -e

VIROME_QC="/Users/scotthandley/.cargo/target/release/virome-qc"
FASTP="conda run fastp"
TMPDIR="benchmark_data/fastp_comparison"
RESULTS="$TMPDIR/results.txt"

mkdir -p "$TMPDIR"
> "$RESULTS"

download_dataset() {
    local accession="$1"
    local name="$2"
    local extra_args="$3"

    local r1="$TMPDIR/${accession}_1.fastq.gz"
    if [ -f "$r1" ]; then
        echo "  $name already downloaded"
        return
    fi

    echo "  Downloading $name ($accession)..."
    fastq-dump "$accession" --split-files --gzip --clip $extra_args --outdir "$TMPDIR" 2>/dev/null
}

run_fastp() {
    local name="$1"
    local r1="$2"
    local r2="$3"
    local outdir="$TMPDIR/fastp_${name}"

    mkdir -p "$outdir"

    echo "  Running fastp on $name..."
    local start=$(date +%s%N)

    if [ -n "$r2" ] && [ -f "$r2" ]; then
        conda run fastp \
            -i "$r1" -I "$r2" \
            -o "$outdir/R1.fastq.gz" -O "$outdir/R2.fastq.gz" \
            --json "$outdir/fastp.json" --html "$outdir/fastp.html" \
            --thread 8 \
            --detect_adapter_for_pe \
            --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
            --length_required 50 \
            2>/dev/null
    else
        conda run fastp \
            -i "$r1" \
            -o "$outdir/R1.fastq.gz" \
            --json "$outdir/fastp.json" --html "$outdir/fastp.html" \
            --thread 8 \
            --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
            --length_required 50 \
            2>/dev/null
    fi

    local end=$(date +%s%N)
    local elapsed=$(( (end - start) / 1000000 ))

    # Parse fastp JSON
    python3 -c "
import json
d = json.load(open('$outdir/fastp.json'))
s = d['summary']
reads_in = s['before_filtering']['total_reads']
reads_out = s['after_filtering']['total_reads']
survival = reads_out / reads_in * 100
adapter = d.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0)
adapter_pct = adapter / reads_in * 100 if reads_in > 0 else 0
print(f'fastp,$name,{reads_in},{reads_out},{survival:.1f},{adapter_pct:.2f},{$elapsed}')
" >> "$RESULTS"

    echo "    fastp done (${elapsed}ms)"
}

run_virome_qc() {
    local name="$1"
    local r1="$2"
    local r2="$3"
    local profile="$4"
    local outdir="$TMPDIR/vqc_${name}"

    mkdir -p "$outdir"

    echo "  Running virome-qc on $name..."
    local start=$(date +%s%N)

    if [ -n "$r2" ] && [ -f "$r2" ]; then
        $VIROME_QC run -p "$profile" -1 "$r1" -2 "$r2" -i . --merge \
            -o "$outdir" --report-only -t 8 2>/dev/null
    else
        $VIROME_QC run -p "$profile" -i "$r1" \
            -o "$outdir" --report-only -t 8 2>/dev/null
    fi

    local end=$(date +%s%N)
    local elapsed=$(( (end - start) / 1000000 ))

    # Parse passport
    python3 -c "
import json
d = json.load(open('$outdir/passport.json'))
ri = d['reads_input']
rp = d['reads_passed']
survival = d['survival_rate'] * 100
adapter = 0
for m in d['modules']:
    if m['name'] == 'adapter':
        adapter = m.get('reads_modified', 0) + m.get('reads_removed', 0)
adapter_pct = adapter / ri * 100 if ri > 0 else 0
print(f'virome-qc,$name,{ri},{rp},{survival:.1f},{adapter_pct:.2f},{$elapsed}')
" >> "$RESULTS"

    echo "    virome-qc done (${elapsed}ms)"
}

echo "============================================="
echo "fastp vs virome-qc Head-to-Head Comparison"
echo "============================================="
echo ""

# Header
echo "tool,dataset,reads_in,reads_out,survival_pct,adapter_pct,time_ms" > "$RESULTS"

# Dataset 1: Cook (MiSeq 2x250, paired-end)
echo "--- Cook (ERR10359658, MiSeq 2x250, PE) ---"
download_dataset ERR10359658 cook ""
R1="$TMPDIR/ERR10359658_1.fastq.gz"
R2="$TMPDIR/ERR10359658_2.fastq.gz"
run_fastp cook "$R1" "$R2"
run_virome_qc cook "$R1" "$R2" stool-vlp-tagmentation

# Dataset 2: Zhang depleted (NovaSeq, 1M reads, PE)
echo "--- Zhang depleted (SRR33419066, NovaSeq, PE, 1M reads) ---"
download_dataset SRR33419066 zhang_depleted "-X 1000000"
R1="$TMPDIR/SRR33419066_1.fastq.gz"
R2="$TMPDIR/SRR33419066_2.fastq.gz"
run_fastp zhang_depleted "$R1" "$R2"
run_virome_qc zhang_depleted "$R1" "$R2" profiles/short-read-nebnext.yaml

# Dataset 3: Shkoporov (gut VLP, PE)
echo "--- Shkoporov (SRR9161520, gut VLP, PE) ---"
download_dataset SRR9161520 shkoporov ""
R1="$TMPDIR/SRR9161520_1.fastq.gz"
R2="$TMPDIR/SRR9161520_2.fastq.gz"
run_fastp shkoporov "$R1" "$R2"
run_virome_qc shkoporov "$R1" "$R2" stool-vlp-tagmentation

echo ""
echo "============================================="
echo "Results"
echo "============================================="
cat "$RESULTS"

echo ""
echo "============================================="
echo "Comparison Table"
echo "============================================="

python3 -c "
import csv, io

data = open('$RESULTS').read()
reader = csv.DictReader(io.StringIO(data))
rows = list(reader)

datasets = sorted(set(r['dataset'] for r in rows))

print(f'{\"\":15s} {\"\":>15s} | {\"fastp\":^30s} | {\"virome-qc\":^30s} | {\"Delta\":^15s}')
print(f'{\"Dataset\":15s} {\"Reads In\":>15s} | {\"Survival\":>10s} {\"Adapter\":>10s} {\"Time\":>8s} | {\"Survival\":>10s} {\"Adapter\":>10s} {\"Time\":>8s} | {\"Surv diff\":>15s}')
print('-' * 120)

for ds in datasets:
    fp = [r for r in rows if r['dataset'] == ds and r['tool'] == 'fastp']
    vq = [r for r in rows if r['dataset'] == ds and r['tool'] == 'virome-qc']
    if not fp or not vq:
        continue
    f, v = fp[0], vq[0]
    delta = float(v['survival_pct']) - float(f['survival_pct'])
    ft = int(f['time_ms']) / 1000
    vt = int(v['time_ms']) / 1000
    print(f'{ds:15s} {int(f[\"reads_in\"]):>15,} | {float(f[\"survival_pct\"]):>9.1f}% {float(f[\"adapter_pct\"]):>9.2f}% {ft:>7.1f}s | {float(v[\"survival_pct\"]):>9.1f}% {float(v[\"adapter_pct\"]):>9.2f}% {vt:>7.1f}s | {delta:>+14.1f}%')

print()
print('NOTE: virome-qc removes more reads because it screens for:')
print('  rRNA, host DNA, PhiX, vectors, ERVs, internal adapters, complexity')
print('  fastp only does: adapter trimming, quality filtering, length filtering')
"

# Clean up FASTQ files
echo ""
echo "Cleaning up downloaded FASTQs..."
rm -f "$TMPDIR"/*.fastq.gz "$TMPDIR"/fastp_*/*.fastq.gz "$TMPDIR"/vqc_*/*.fastq.gz
echo "Done."
