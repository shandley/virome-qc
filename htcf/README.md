# HTCF Benchmark Setup

## Quick Start

```bash
# 1. SSH to HTCF
ssh shandley@login.htcf.wustl.edu

# 2. Run setup (clones repos, builds virome-qc, downloads reference databases)
# This takes ~30 min (mostly building the T2T host filter)
cd /scratch/sahlab/shandley/virome-qc-benchmark/src/virome-qc
bash htcf/setup.sh

# 3. Submit Tier 1 benchmarks (11 datasets)
bash htcf/run_tier1.sh

# 4. Monitor
squeue -u shandley
ls /scratch/sahlab/shandley/virome-qc-benchmark/results/*/passport.json | wc -l

# 5. Collect results
python3 htcf/collect_results.py
```

## Directory Structure

```
/scratch/sahlab/shandley/virome-qc-benchmark/
├── src/
│   ├── virome-qc/          # virome-qc repo (with cluster Cargo.toml paths)
│   ├── biometal/            # biometal dependency
│   └── SuperBloom/          # SuperBloom dependency
├── databases/
│   ├── human_t2t.sbf       # T2T-CHM13 Super Bloom filter (~4.1 GB)
│   └── rrna_silva.rrf      # SILVA rRNA sorted hash filter (~1.3 GB)
├── results/
│   ├── cook_mock/
│   │   ├── passport.json
│   │   └── report.html
│   ├── shkoporov_gut/
│   │   └── ...
│   └── ...
├── tmp/                     # Per-sample FASTQ downloads (auto-cleaned)
└── logs/                    # SLURM job logs
```

## Resource Requirements

- **Per job:** 8 CPUs, 12 GB RAM, 4 hour time limit
- **Storage:** ~5-50 GB per sample during download (auto-cleaned after)
- **Persistent storage:** ~100 KB per passport + ~600 KB per report = ~10 MB total for Tier 1

## Tiers

- **Tier 1 (11 jobs):** Core benchmark datasets, ~2-4 hours total
- **Tier 2 (10 jobs):** Extended validation, submit after Tier 1 verification
- **Tier 3 (5000 jobs):** Atlas generation, submit as large array job
