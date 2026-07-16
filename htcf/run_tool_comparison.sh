#!/bin/bash
#SBATCH --job-name=qc-compare
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/sahlab/shandley/virome-qc-benchmark/logs/qc_compare_%A_%a.out
#SBATCH --error=/scratch/sahlab/shandley/virome-qc-benchmark/logs/qc_compare_%A_%a.err
#
# Run fastp, BBDuk, Trimmomatic, and virome-qc on the same datasets for comparison.
#
# Submit: sbatch --array=1-3 htcf/run_tool_comparison.sh
#   1 = Cook mock (MiSeq 2x250, Nextera, clean VLP, MDA amplified)
#   2 = Shkoporov gut (NovaSeq 2x150, stool VLP, high duplication)
#   3 = Zhang depleted (NextSeq 2x75, stool, high rRNA + poly-G)

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
SRCDIR="$WORKDIR/src/virome-qc"
VIROME_QC="$SRCDIR/target/release/virome-qc"
ENVBIN="/ref/sahlab/software/scott_conda/virome-qc-bench/bin"
OUTDIR="$WORKDIR/tool_comparison"
TMPDIR_BASE="$WORKDIR/tmp"

export PATH="$ENVBIN:$PATH"
source ~/.cargo/env 2>/dev/null || true

# Adapter references
BBDUK_ADAPTERS="$ENVBIN/../opt/bbmap-39.01-1/resources/adapters.fa"
TRIMJAR="$ENVBIN/../share/trimmomatic-0.40-0/trimmomatic.jar"
TRIM_NEXTERA="$ENVBIN/../share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa"
TRIM_TRUSEQ="$ENVBIN/../share/trimmomatic-0.40-0/adapters/TruSeq3-PE-2.fa"

mkdir -p "$OUTDIR" "$TMPDIR_BASE" "$WORKDIR/logs"

# Dataset definitions
DATASETS=(
    "1|cook|ERR10359658|PE|stool-vlp-tagmentation|nextera"
    "2|shkoporov|SRR9161520|PE|stool-vlp-tagmentation|nextera"
    "3|zhang_depleted|SRR33419066|PE|profiles/short-read-nebnext.yaml|truseq"
)

IDX=${SLURM_ARRAY_TASK_ID:-${1:-1}}
DATASET="${DATASETS[$((IDX-1))]}"
IFS='|' read -r _ NAME ACCESSION PAIRED PROFILE ADAPTER_TYPE <<< "$DATASET"

echo "================================================================"
echo "Tool Comparison: $NAME ($ACCESSION)"
echo "Profile: $PROFILE | Adapters: $ADAPTER_TYPE"
echo "Time: $(date)"
echo "================================================================"

# Select Trimmomatic adapter file
if [ "$ADAPTER_TYPE" = "nextera" ]; then
    TRIM_ADAPTER="$TRIM_NEXTERA"
else
    TRIM_ADAPTER="$TRIM_TRUSEQ"
fi

# Download FASTQ
TMPDIR="$TMPDIR_BASE/$NAME"
mkdir -p "$TMPDIR"
cd "$TMPDIR"

echo ""
echo "=== Downloading $ACCESSION ==="
if [ ! -f "${ACCESSION}_1.fastq.gz" ]; then
    prefetch "$ACCESSION" --max-size 100G -O "$TMPDIR" 2>&1 | tail -3
    fasterq-dump "$TMPDIR/$ACCESSION/$ACCESSION.sra" --split-files --threads 4 --temp "$TMPDIR" --outdir "$TMPDIR" 2>&1 | tail -5
    if [ ! -f "${ACCESSION}_1.fastq" ]; then
        fasterq-dump "$ACCESSION" --split-files --threads 4 --temp "$TMPDIR" --outdir "$TMPDIR" 2>&1 | tail -5
    fi
    pigz -p 4 "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq" 2>/dev/null || \
    gzip "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
fi
R1="$TMPDIR/${ACCESSION}_1.fastq.gz"
R2="$TMPDIR/${ACCESSION}_2.fastq.gz"

READS_IN=$(zcat "$R1" | awk 'END{print NR/4}')
echo "Input reads: $READS_IN (per file)"

SAMPLE_OUT="$OUTDIR/$NAME"
mkdir -p "$SAMPLE_OUT"

# ============================================================
# 1. fastp (default parameters - adapter auto-detection)
# ============================================================
echo ""
echo "=== fastp ==="
FASTP_OUT="$SAMPLE_OUT/fastp"
mkdir -p "$FASTP_OUT"
START=$(date +%s%N)

"$ENVBIN/fastp" \
    -i "$R1" -I "$R2" \
    -o "$FASTP_OUT/R1.fastq.gz" -O "$FASTP_OUT/R2.fastq.gz" \
    -j "$FASTP_OUT/fastp.json" -h "$FASTP_OUT/fastp.html" \
    --thread 8 2>&1 | tail -5

END=$(date +%s%N)
FASTP_MS=$(( (END - START) / 1000000 ))

# Extract metrics from fastp JSON
python3 << PYEOF
import json
d = json.load(open("$FASTP_OUT/fastp.json"))
s = d["summary"]
filt = d["filtering_result"]
ri = s["before_filtering"]["total_reads"]
ro = s["after_filtering"]["total_reads"]
adapt = d.get("adapter_cutting", {})
adapt_r1 = adapt.get("adapter_trimmed_reads", 0)
# fastp reports adapter-trimmed reads (not pairs)
adapt_pct = adapt_r1 / max(ri, 1) * 100
surv = ro / max(ri, 1) * 100
print(f"fastp: {ri} in, {ro} out, {surv:.1f}% survival, {adapt_pct:.2f}% adapter, {$FASTP_MS}ms")
with open("$FASTP_OUT/metrics.json", "w") as f:
    json.dump({
        "tool": "fastp",
        "dataset": "$NAME",
        "reads_in": ri,
        "reads_out": ro,
        "survival_pct": round(surv, 2),
        "adapter_pct": round(adapt_pct, 2),
        "time_ms": $FASTP_MS,
    }, f, indent=2)
PYEOF

# ============================================================
# 2. BBDuk (adapter + quality + entropy filtering)
# ============================================================
echo ""
echo "=== BBDuk ==="
BBDUK_OUT="$SAMPLE_OUT/bbduk"
mkdir -p "$BBDUK_OUT"
START=$(date +%s%N)

"$ENVBIN/bbduk.sh" \
    in1="$R1" in2="$R2" \
    out1="$BBDUK_OUT/R1.fastq.gz" out2="$BBDUK_OUT/R2.fastq.gz" \
    ref="$BBDUK_ADAPTERS" \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=r trimq=15 \
    minlength=50 \
    entropy=0.5 \
    threads=8 -Xmx8g \
    stats="$BBDUK_OUT/stats.txt" 2>&1 | tee "$BBDUK_OUT/bbduk.log" | tail -10

END=$(date +%s%N)
BBDUK_MS=$(( (END - START) / 1000000 ))

# Extract metrics from BBDuk log
python3 << PYEOF
import re
log = open("$BBDUK_OUT/bbduk.log").read()

ri_match = re.search(r"Input:\s+(\d+) reads", log)
ro_match = re.search(r"Result:\s+(\d+) reads", log)
adapt_match = re.search(r"Total Removed:\s+(\d+) reads", log)

ri = int(ri_match.group(1)) if ri_match else $READS_IN * 2
ro = int(ro_match.group(1)) if ro_match else 0
removed = int(adapt_match.group(1)) if adapt_match else ri - ro
surv = ro / max(ri, 1) * 100

# BBDuk adapter stats
adapt_reads = 0
adapt_match2 = re.search(r"(\d+) reads .* adapters", log)
if adapt_match2:
    adapt_reads = int(adapt_match2.group(1))
adapt_pct = adapt_reads / max(ri, 1) * 100

print(f"BBDuk: {ri} in, {ro} out, {surv:.1f}% survival, {adapt_pct:.2f}% adapter, {$BBDUK_MS}ms")

import json
with open("$BBDUK_OUT/metrics.json", "w") as f:
    json.dump({
        "tool": "bbduk",
        "dataset": "$NAME",
        "reads_in": ri,
        "reads_out": ro,
        "survival_pct": round(surv, 2),
        "adapter_pct": round(adapt_pct, 2),
        "time_ms": $BBDUK_MS,
    }, f, indent=2)
PYEOF

# ============================================================
# 3. Trimmomatic (adapter + quality + length filtering)
# ============================================================
echo ""
echo "=== Trimmomatic ==="
TRIM_OUT="$SAMPLE_OUT/trimmomatic"
mkdir -p "$TRIM_OUT"
START=$(date +%s%N)

java -jar "$TRIMJAR" PE \
    -threads 8 -phred33 \
    -summary "$TRIM_OUT/summary.txt" \
    "$R1" "$R2" \
    "$TRIM_OUT/R1_paired.fastq.gz" "$TRIM_OUT/R1_unpaired.fastq.gz" \
    "$TRIM_OUT/R2_paired.fastq.gz" "$TRIM_OUT/R2_unpaired.fastq.gz" \
    ILLUMINACLIP:${TRIM_ADAPTER}:2:30:10:2:True \
    LEADING:3 TRAILING:15 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50 2>&1 | tee "$TRIM_OUT/trimmomatic.log" | tail -5

END=$(date +%s%N)
TRIM_MS=$(( (END - START) / 1000000 ))

# Extract metrics from Trimmomatic log
python3 << PYEOF
import re
log = open("$TRIM_OUT/trimmomatic.log").read()

# "Input Read Pairs: 1268684 Both Surviving: 1232508 (97.15%) ..."
m = re.search(r"Input Read Pairs: (\d+) Both Surviving: (\d+) .* Forward Only Surviving: (\d+)", log)
if m:
    pairs_in = int(m.group(1))
    pairs_surv = int(m.group(2))
    fwd_only = int(m.group(3))
    ri = pairs_in * 2
    ro = pairs_surv * 2 + fwd_only
else:
    ri = $READS_IN * 2
    ro = 0
    pairs_in = $READS_IN

surv = ro / max(ri, 1) * 100
print(f"Trimmomatic: {ri} in, {ro} out, {surv:.1f}% survival, {$TRIM_MS}ms")

import json
with open("$TRIM_OUT/metrics.json", "w") as f:
    json.dump({
        "tool": "trimmomatic",
        "dataset": "$NAME",
        "reads_in": ri,
        "reads_out": ro,
        "survival_pct": round(surv, 2),
        "adapter_pct": 0,  # Trimmomatic doesn't report adapter rate separately
        "time_ms": $TRIM_MS,
    }, f, indent=2)
PYEOF

# ============================================================
# 4. virome-qc
# ============================================================
echo ""
echo "=== virome-qc ==="
VQC_OUT="$SAMPLE_OUT/virome-qc"
mkdir -p "$VQC_OUT"
START=$(date +%s%N)

cd "$SRCDIR"
$VIROME_QC run -p "$PROFILE" -1 "$R1" -2 "$R2" -i . --merge \
    -o "$VQC_OUT" --report-only -t 8 2>&1 | tail -5

END=$(date +%s%N)
VQC_MS=$(( (END - START) / 1000000 ))

# Extract metrics from passport
python3 << PYEOF
import json
d = json.load(open("$VQC_OUT/passport.json"))
ri = d["reads_input"]
ro = d["reads_passed"]
surv = d["survival_rate"] * 100
qc_surv = d.get("qc_survival_rate", d["survival_rate"]) * 100

adapt = 0
for m in d["modules"]:
    if m["name"] == "adapter":
        adapt = m.get("extra", {}).get("adapters_found_3prime", 0) / max(ri, 1) * 100

dedup = sum(m["reads_removed"] for m in d["modules"] if m["name"] == "dedup") / max(ri, 1) * 100
host = sum(m["reads_removed"] for m in d["modules"] if m["name"] == "host") / max(ri, 1) * 100
rrna = sum(m["reads_removed"] for m in d["modules"] if m["name"] == "rrna") / max(ri, 1) * 100

print(f"virome-qc: {ri} in, {ro} out, {surv:.1f}% surv ({qc_surv:.1f}% QCSurv), {adapt:.2f}% adapt, {dedup:.1f}% dedup, {host:.2f}% host, {rrna:.2f}% rRNA, {$VQC_MS}ms")

with open("$VQC_OUT/metrics.json", "w") as f:
    json.dump({
        "tool": "virome-qc",
        "dataset": "$NAME",
        "reads_in": ri,
        "reads_out": ro,
        "survival_pct": round(surv, 2),
        "qc_survival_pct": round(qc_surv, 2),
        "adapter_pct": round(adapt, 2),
        "dedup_pct": round(dedup, 2),
        "host_pct": round(host, 2),
        "rrna_pct": round(rrna, 2),
        "time_ms": $VQC_MS,
    }, f, indent=2)
PYEOF

# ============================================================
# Summary
# ============================================================
echo ""
echo "================================================================"
echo "Summary: $NAME"
echo "================================================================"
for tool in fastp bbduk trimmomatic virome-qc; do
    MF="$SAMPLE_OUT/$tool/metrics.json"
    if [ -f "$MF" ]; then
        python3 -c "
import json
d = json.load(open('$MF'))
extra = ''
if 'qc_survival_pct' in d:
    extra = f' (QCSurv {d[\"qc_survival_pct\"]}%)'
if 'dedup_pct' in d and d['dedup_pct'] > 0:
    extra += f' dedup={d[\"dedup_pct\"]}%'
if 'host_pct' in d and d['host_pct'] > 0:
    extra += f' host={d[\"host_pct\"]}%'
if 'rrna_pct' in d and d['rrna_pct'] > 0:
    extra += f' rRNA={d[\"rrna_pct\"]}%'
print(f'  {d[\"tool\"]:<14} {d[\"reads_in\"]:>10} -> {d[\"reads_out\"]:>10}  {d[\"survival_pct\"]:>6.1f}%  adapt={d[\"adapter_pct\"]:.2f}%  {d[\"time_ms\"]/1000:.1f}s{extra}')
"
    fi
done

# Clean up downloaded FASTQ
echo ""
echo "=== Cleaning up ==="
rm -rf "$TMPDIR"

echo ""
echo "=== Done: $NAME ($(date)) ==="
