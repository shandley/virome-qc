#!/bin/bash
# Collect tool comparison results into a single CSV.
# Run on login node after run_tool_comparison.sh completes.

OUTDIR="/scratch/sahlab/shandley/virome-qc-benchmark/tool_comparison"

echo "=== Tool Comparison Results ==="
echo ""

python3 << 'PYEOF'
import json, glob, os

results = []
for mf in sorted(glob.glob("/scratch/sahlab/shandley/virome-qc-benchmark/tool_comparison/*/*/metrics.json")):
    d = json.load(open(mf))
    results.append(d)

if not results:
    print("No results found yet.")
    exit(0)

# Print table
print(f"{'Tool':<14} {'Dataset':<18} {'Reads In':>10} {'Reads Out':>10} {'Surv%':>7} {'QCSurv%':>8} {'Adapt%':>7} {'Dedup%':>7} {'Host%':>7} {'rRNA%':>7} {'Time(s)':>8}")
print("-" * 120)
for d in sorted(results, key=lambda x: (x["dataset"], x["tool"])):
    qcs = d.get("qc_survival_pct", "")
    dup = d.get("dedup_pct", "")
    host = d.get("host_pct", "")
    rrna = d.get("rrna_pct", "")
    print(f"{d['tool']:<14} {d['dataset']:<18} {d['reads_in']:>10} {d['reads_out']:>10} {d['survival_pct']:>7.1f} {str(qcs):>8} {d['adapter_pct']:>7.2f} {str(dup):>7} {str(host):>7} {str(rrna):>7} {d['time_ms']/1000:>8.1f}")

# Write CSV
csv_path = "/scratch/sahlab/shandley/virome-qc-benchmark/tool_comparison/results.csv"
with open(csv_path, "w") as f:
    f.write("tool,dataset,reads_in,reads_out,survival_pct,qc_survival_pct,adapter_pct,dedup_pct,host_pct,rrna_pct,time_ms\n")
    for d in sorted(results, key=lambda x: (x["dataset"], x["tool"])):
        f.write(f"{d['tool']},{d['dataset']},{d['reads_in']},{d['reads_out']},{d['survival_pct']},{d.get('qc_survival_pct','')},{d['adapter_pct']},{d.get('dedup_pct','')},{d.get('host_pct','')},{d.get('rrna_pct','')},{d['time_ms']}\n")

print(f"\nCSV: {csv_path}")
PYEOF
