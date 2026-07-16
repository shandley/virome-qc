# virome-qc TODO

## Paper Submission
- [x] Final benchmark regeneration — Tier 1 complete on HTCF (13/13, job 38048929)
- [x] Evaluate and finalize benchmark dataset selection — BENCHMARK_DATASETS.md
- [x] Embed ViroForge reference ranges into profiles
- [x] Formal benchmark evaluation and weakness assessment
- [x] Separate dedup from QC survival rate (quality tier based on unique reads)
- [x] Conservative merge with ambiguity detection
- [x] ViroForge expected-range calibration — 20 datasets, 4 profiles, HTCF job 38149974/38149975
- [x] Update ExpectedRanges struct: added duplication_rate, gc_content, QCSurv for survival
- [x] Update profile ranges with ViroForge-derived values (all 4 profiles + 2 YAML files)
- [x] Fix WGA adapter_rate extraction bug (reads_modified inflated by random_primer_trim)
- [x] Reports show user metrics vs expected range with color-coding (UI update)
- [x] Chart clipping for low-coverage positions (quality profile + base composition)
- [ ] Generate publication-quality figures from existing data
- [ ] Formal sensitivity/specificity with bootstrap CIs from ViroForge ground truth
- [x] Add BBDuk/Trimmomatic to fastp comparison (HTCF job 38213984, shared-scope analysis)
- [ ] Write manuscript

## User Readiness
- [x] `virome-qc report` recomputes flags from passport data (single source of truth)
- [x] `virome-qc db --setup` auto-downloads SILVA + T2T and builds filters
- [x] Better error messages when filters are missing
- [x] Quickstart guide in README
- [x] Duplicate flag bug fixed
- [x] SE/PE ambiguity detection in ingestion engine
- [x] Dedup enabled by default in all profiles
- [x] Reports show expected ranges from profile (stat cards + host/rRNA cards)
- [x] CLI input validation (missing file, PE without -i, ingest -i flag)
- [x] Output directory overwrite protection (--force flag)
- [x] Profile descriptions in CLI (virome-qc profiles --show)
- [x] Document survival_rate vs qc_survival_rate in README

## Platform Expansion
- [ ] Long-read QC modules (ONT adapter trimming, PacBio CCS quality filtering)
- [ ] Probe-capture virome support (separate modality, v2 — requires panel BED, alignment, on-target metrics)
- [ ] Clinical diagnostics mode (negative control comparison, expanded contaminants, evidence thresholds)
- [x] Decommission the three-signal ERV endo/exo classifier [2026-07-16] (removed run_erv_analysis, erv_classifier, erv_pipeline, eve_polymorphisms, and the report card; retroviral read detection retained)
- [ ] Alignment-based EVE discrimination screen (subtract-then-discriminate) -- IN PROGRESS: DIAMOND (eve_prot/set_d_prot) + minimap2 (eve_nt/set_d_nt) reference DBs built and validated (K113 held-out recall + negative controls); tool integration pending. See EVE_SCREEN_MODULE.md, EVE_DISCRIMINATION_REDESIGN.md
- [ ] Fold detectEVE nrEVE loci into the EVE screen references, then wire the screen into the pipeline
- [ ] Non-retroviral EVE (nrEVE) detection -- the embedded Bornaviridae/Parvoviridae/Filoviridae sequences are exogenous references (set-D panel), not host EVEs; real nrEVE detection is outstanding
- [x] EVE-aware host-depletion rescue (eve_aware flag) -- exists but limited to the HERV-K K113 / ciHHV-6 loci; k-mer exclusion as an FP screen was rejected (redundant with T2T host depletion) in favor of the alignment screen
- [ ] Multi-sample comparison dashboard
- [ ] CI/CD regression testing with ViroForge canonical dataset

## Atlas / Community
- [ ] Publish ViroForge reference datasets + virome-qc passports as public resource
- [ ] Auto-classification: "your sample most closely matches X"
- [ ] Community-contributed sample type profiles
- [ ] Zenodo deposition for SILVA and T2T filter files

## HTCF Benchmarking
- [x] Tier 1: 13 core datasets complete (job 38048929)
- [x] Tier 1 evaluation complete
- [x] ViroForge calibration: 20 datasets generated + QC'd (jobs 38149974/38149975)
- [x] Expected ranges derived for all 4 profiles (see VALIDATION.md)
- [ ] Tier 2: 10 extended datasets
- [ ] Tier 3: ~5,000 samples for atlas generation
