#!/bin/bash
#SBATCH --job-name=merge-cmp
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --output=logs/merge_cmp_%j.out

set -e
source /ref/sahlab/software/miniforge3/bin/activate virome-qc-bench
source ~/.cargo/env

WORKDIR=/scratch/sahlab/shandley/virome-qc-benchmark
VIROME_QC=$WORKDIR/src/virome-qc/target/release/virome-qc
DIR=$WORKDIR/merge_comparison
mkdir -p $DIR
cd $DIR

# Download Cook MiSeq 2x250
echo "=== Downloading Cook dataset ==="
if [ ! -f ERR10359658_1.fastq.gz ]; then
    fasterq-dump ERR10359658 --split-files --threads 4 --temp $DIR
    pigz -p 4 ERR10359658_1.fastq ERR10359658_2.fastq
fi

TOTAL_PAIRS=$(zcat ERR10359658_1.fastq.gz | awk 'END{print NR/4}')
echo "Total pairs: $TOTAL_PAIRS"

# Run BBMerge
echo ""
echo "=== BBMerge ==="
START=$(date +%s)
bbmerge.sh \
    in1=ERR10359658_1.fastq.gz in2=ERR10359658_2.fastq.gz \
    out=bbmerge_merged.fastq.gz \
    outu1=bbmerge_unmerged_R1.fastq.gz outu2=bbmerge_unmerged_R2.fastq.gz \
    ihist=bbmerge_insert_hist.txt \
    threads=8 -Xmx8g 2>&1
END=$(date +%s)
BB_TIME=$((END - START))

BB_MERGED=$(zcat bbmerge_merged.fastq.gz | awk 'END{print NR/4}')
BB_UNMERGED=$(zcat bbmerge_unmerged_R1.fastq.gz | awk 'END{print NR/4}')
BB_RATE=$(python3 -c "print(f'{$BB_MERGED / $TOTAL_PAIRS * 100:.1f}')")
echo ""
echo "BBMerge: merged=$BB_MERGED unmerged=$BB_UNMERGED rate=${BB_RATE}% time=${BB_TIME}s"

# Run virome-qc
echo ""
echo "=== virome-qc ==="
cd $WORKDIR/src/virome-qc
START=$(date +%s)
$VIROME_QC run \
    -p stool-vlp-tagmentation \
    -1 $DIR/ERR10359658_1.fastq.gz -2 $DIR/ERR10359658_2.fastq.gz \
    -i . --merge \
    -o $DIR/vqc_output \
    --report-only -t 8 2>&1
END=$(date +%s)
VQC_TIME=$((END - START))

# Parse passport
python3 << 'PYEOF'
import json
d = json.load(open("/scratch/sahlab/shandley/virome-qc-benchmark/merge_comparison/vqc_output/passport.json"))
ri = d["reads_input"]
pp = d.get("pairs_passed", 0)
merged = d.get("pairs_merged", 0)
singles = d.get("singletons", 0)
surv = d["survival_rate"]
print(f"virome-qc: reads_in={ri} pairs_passed={pp} merged={merged} singletons={singles} survival={surv*100:.1f}%")
print(f"virome-qc merge rate: {merged / max(pp, 1) * 100:.1f}% of passing pairs")
PYEOF

echo "virome-qc time: ${VQC_TIME}s"

# Compare merged read length distributions
echo ""
echo "=== Merged read length distributions ==="
echo "BBMerge (top 10 lengths):"
zcat bbmerge_merged.fastq.gz | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | sort -rn | head -10

echo ""
echo "virome-qc (top 10 lengths):"
if [ -f $DIR/vqc_output/merged.fastq.gz ]; then
    zcat $DIR/vqc_output/merged.fastq.gz | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | sort -rn | head -10
else
    echo "No merged output file"
fi

# Summary
echo ""
echo "=== SUMMARY ==="
echo "Dataset: Cook MiSeq 2x250, $TOTAL_PAIRS pairs"
echo "BBMerge:   $BB_MERGED merged (${BB_RATE}%), ${BB_TIME}s"
VQC_MERGED=$(python3 -c "import json; d=json.load(open('$DIR/vqc_output/passport.json')); print(d.get('pairs_merged',0))")
VQC_RATE=$(python3 -c "import json; d=json.load(open('$DIR/vqc_output/passport.json')); pp=d.get('pairs_passed',0); m=d.get('pairs_merged',0); print(f'{m/max(pp,1)*100:.1f}')")
echo "virome-qc: $VQC_MERGED merged (${VQC_RATE}% of passing pairs), ${VQC_TIME}s"

# Cleanup
rm -f $DIR/ERR10359658_*.fastq.gz
rm -f $DIR/bbmerge_*.fastq.gz $DIR/bbmerge_*.txt
rm -rf $DIR/vqc_output
