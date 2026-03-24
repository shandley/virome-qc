# virome-qc TODO

## Paper Submission
- [x] Final benchmark regeneration with latest code — Tier 1 submitted on HTCF (job 38021183)
- [x] Evaluate and finalize benchmark dataset selection — BENCHMARK_DATASETS.md (3 tiers, 25+ datasets)
- [x] Embed ViroForge reference ranges into profiles — expected_ranges in all 4 builtins + 2 YAMLs
- [ ] Reports show user metrics vs expected range with color-coding — ranges in passport, UI needs update
- [ ] Generate publication-quality figures from existing data
- [ ] Formal sensitivity/specificity with bootstrap CIs from ViroForge ground truth
- [ ] Add BBDuk and/or Trimmomatic to the fastp comparison
- [ ] Write manuscript

## User Readiness
- [x] `virome-qc report` recomputes flags from passport data
- [x] `virome-qc db --setup` auto-downloads SILVA + T2T and builds filters
- [x] Better error messages when filters are missing (includes build instructions)
- [x] Quickstart guide in README
- [ ] Reports show expected ranges from profile (UI update needed)

## Platform Expansion
- [ ] Long-read QC modules (ONT adapter trimming, PacBio CCS quality filtering)
- [ ] Population ERV database from HPRC pangenome (known polymorphic HERV insertion sites)
- [ ] Coordinate-based EVE flagging at known high-identity integration sites (ciHHV-6, recent HERV-K)
- [ ] Multi-sample comparison dashboard (beyond basic batch report)
- [ ] CI/CD regression testing with ViroForge canonical dataset

## Atlas / Community
- [ ] Publish ViroForge reference datasets + virome-qc passports as public resource
- [ ] Auto-classification: "your sample most closely matches X, with deviations in Y and Z"
- [ ] Community-contributed sample type profiles
- [ ] Zenodo deposition for SILVA and T2T filter files

## HTCF Benchmarking
- [x] Tier 1: 15 core datasets submitted (job 38021183)
- [ ] Tier 1: Evaluate results
- [ ] Tier 2: 10 extended datasets
- [ ] Tier 3: ~5,000 samples for atlas generation
