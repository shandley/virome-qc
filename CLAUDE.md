# virome-qc

High-performance virome-specific QC platform (Rust). Streaming FASTQ QC pipeline that emits a
QA "passport" (JSON) and a self-contained HTML report. This file holds project-specific gotchas;
general Rust/workflow conventions live in the global config.

## Build and test

- Local path deps: `biometal` and `SuperBloom` at `/Users/shandley/Code/software/`. Run
  `cargo clippy` after changes, `cargo nextest run` for tests.
- Building on HTCF (x86): a transitive dep (`ensure_simd`) refuses to compile unless SIMD target
  features are enabled. Set `RUSTFLAGS="-C target-cpu=x86-64-v3"` or add a `.cargo/config.toml`
  with it. Apple Silicon does not need this (NEON is always present). HTCF has biometal/SuperBloom
  at `/ref/sahlab/software/{biometal,SuperBloom}` (versions match local); repoint the Cargo.toml
  paths there for an HTCF build.

## Filters (host + rRNA)

Auto-discovered from `~/.virome-qc/host/human.sbf` and `~/.virome-qc/rrna/silva.rrf` (or the
`VIROME_QC_DB` env var). On HTCF the filters live at
`/ref/sahlab/data/virome-qc-db/{host/human.sbf, rrna/silva.rrf}` and `~/.virome-qc` already
symlinks to them. The host `reference` name is `human`; the ambiguous containment threshold is
0.20 (host >= 0.50 remove, 0.20-0.50 ambiguous, < 0.20 keep).

## Current state (2026-07)

- The 3-signal ERV endo/exo classifier is DECOMMISSIONED. Retroviral read DETECTION
  (`ErvScreener`, a k-mer flag) remains. Do not describe the classifier as a working feature.
- Active EVE work is the discrimination redesign: an ALIGNMENT-based competitive screen
  (minimap2 + DIAMOND, EVE vs exogenous). References are built and validated (K113 100%
  held-out recall) but the screen is NOT yet wired into the tool. See EVE_SCREEN_MODULE.md,
  EVE_DISCRIMINATION_REDESIGN.md, REVIEW_FINDINGS.md (open Tier 1/2 items), TODO.md.
- Profile `expected_ranges` were re-derived 2026-07 from a fresh 20-dataset ViroForge calibration.

## Config

Profiles are defined in two places: `src/config/profiles.rs` (all profiles) and `profiles/*.yaml`
(only some). `expected_ranges` live in both, so update both when recalibrating.

## HTCF calibration (ViroForge)

- ViroForge is at `/scratch/sahlab/shandley/virome-qc-benchmark/src/viroforge`. Its only runtime
  data file is `viroforge/data/viral_genomes.db` (~500M, 14,423 genomes, all 20 collections).
  `iss` (InSilicoSeq) is in the `viroforge-bench` conda env.
- Scripts: `htcf/{run_calibration,generate_references,run_qc_on_references,derive_ranges}.sh`.
  `derive_ranges` reads `qc_survival_rate` (the metric the tool checks), not raw `survival_rate`.
- Gotcha: generation at coverage 10 takes ~1h+ per dataset; packing 5 datasets into one array job
  overruns the 6h SLURM wall. Use one dataset per array task.

## Repo notes

- `references/` holds large DIAMOND/minimap2 DBs; keep them out of commits. `benchmark_data/` is
  gitignored. Superseded design docs (old ERV classifier, EVE-map, publication plan) are in
  `archive/`.
