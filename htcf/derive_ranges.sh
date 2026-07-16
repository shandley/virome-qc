#!/bin/bash
# Derive expected ranges from ViroForge QC results.
# Run after run_qc_on_references.sh completes (all 20 jobs).
#
# This produces YAML files that can be pasted into virome-qc profiles.
# Run on login node (lightweight, just reads JSON files).

set -e

WORKDIR="/scratch/sahlab/shandley/virome-qc-benchmark"
VFDIR="$WORKDIR/src/viroforge"
RESULTS="$WORKDIR/viroforge_qc_results"
RANGEDIR="$WORKDIR/derived_ranges"

# Set up environment
ENVBIN="/ref/sahlab/software/scott_conda/viroforge-bench/bin"
export PATH="$ENVBIN:$PATH"

mkdir -p "$RANGEDIR"

echo "=== Deriving expected ranges from ViroForge QC results ==="
echo ""

# Check what results are available
PROFILES=("stool_vlp_tagmentation" "tissue_truseq" "metagenomics_nextera" "low_biomass_wga")
VQC_NAMES=("stool-vlp-tagmentation" "tissue-truseq" "metagenomics-nextera" "low-biomass-wga")

for i in "${!PROFILES[@]}"; do
    PROFILE="${PROFILES[$i]}"
    VQC_NAME="${VQC_NAMES[$i]}"
    PROFILE_RESULTS="$RESULTS/$PROFILE"

    echo "--- $VQC_NAME ---"

    # Count passports
    N_PASSPORTS=$(find "$PROFILE_RESULTS" -name "passport.json" 2>/dev/null | wc -l)
    echo "  Passports found: $N_PASSPORTS"

    if [ "$N_PASSPORTS" -eq 0 ]; then
        echo "  SKIP: no results yet"
        continue
    fi

    # Extract metrics from passports
    python3 << PYEOF
import json, glob, os

passports = sorted(glob.glob("$PROFILE_RESULTS/*/passport.json"))
print(f"  Processing {len(passports)} passports...")

metrics = []
for p in passports:
    d = json.load(open(p))
    ri = max(d["reads_input"], 1)
    host = sum(m["reads_removed"] for m in d["modules"] if m["name"] == "host") / ri
    rrna = sum(m["reads_removed"] for m in d["modules"] if m["name"] == "rrna") / ri
    dedup = sum(m["reads_removed"] for m in d["modules"] if m["name"] == "dedup") / ri
    adapter_mod = [m for m in d["modules"] if m["name"] == "adapter"]
    adapter = 0
    if adapter_mod:
        am = adapter_mod[0]
        # Use adapters_found_3prime (correct), not reads_modified (inflated by primer trimming)
        adapter = am.get("extra", {}).get("adapters_found_3prime", 0) / ri

    # The profile `survival` range is QCSurv (qc_survival_rate: survival of unique
    # reads, dedup excluded from the denominator), NOT raw survival_rate. The tool
    # checks qc_survival_rate against this range, so derive it from the same metric.
    qcsurv = d.get("qc_survival_rate", d["survival_rate"])
    name = os.path.basename(os.path.dirname(p))
    print(f"    {name}: qcsurv={qcsurv:.4f} host={host:.4f} rrna={rrna:.4f} adapt={adapter:.4f} dedup={dedup:.4f}")

    metrics.append({
        "name": name,
        "survival": qcsurv,
        "host_fraction": host,
        "rrna_fraction": rrna,
        "adapter_rate": adapter,
        "duplication_rate": dedup,
    })

# Compute ranges (p5-p95, or min-max if <5 datasets)
import statistics
fields = ["survival", "host_fraction", "rrna_fraction", "adapter_rate", "duplication_rate"]
n = len(metrics)

print()
print(f"  # Expected ranges for $VQC_NAME")
print(f"  # Derived from {n} ViroForge reference datasets")
print(f"  expected_ranges:")
for field in fields:
    values = sorted([m[field] for m in metrics])
    if n >= 5:
        lo = values[max(0, int(n * 0.05))]
        hi = values[min(n - 1, int(n * 0.95))]
    else:
        lo = values[0]
        hi = values[-1]
    print(f"    {field}: [{lo:.4f}, {hi:.4f}]")

# Write YAML file
with open("$RANGEDIR/${VQC_NAME}.yaml", "w") as f:
    f.write(f"# Expected ranges for $VQC_NAME\n")
    f.write(f"# Derived from {n} ViroForge reference datasets\n")
    f.write(f"# Generated: $(date -I)\n")
    f.write(f"expected_ranges:\n")
    for field in fields:
        values = sorted([m[field] for m in metrics])
        if n >= 5:
            lo = values[max(0, int(n * 0.05))]
            hi = values[min(n - 1, int(n * 0.95))]
        else:
            lo = values[0]
            hi = values[-1]
        f.write(f"  {field}: [{lo:.4f}, {hi:.4f}]\n")

print(f"  Written: $RANGEDIR/${VQC_NAME}.yaml")
PYEOF

    echo ""
done

echo "=== Range derivation complete ==="
echo "Output: $RANGEDIR/"
ls -la "$RANGEDIR/"
echo ""
echo "Next steps:"
echo "  1. Review the YAML files in $RANGEDIR/"
echo "  2. Copy derived ranges into virome-qc profiles (profiles.rs + YAML files)"
echo "  3. Rebuild virome-qc and re-run reports"
