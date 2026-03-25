# virome-qc TODO

## Paper Submission
- [x] Final benchmark regeneration with latest code — Tier 1 complete on HTCF (13/13 datasets)
- [x] Evaluate and finalize benchmark dataset selection — BENCHMARK_DATASETS.md
- [x] Embed ViroForge reference ranges into profiles
- [x] Formal benchmark evaluation and weakness assessment (documented)
- [ ] Reports show user metrics vs expected range with color-coding (UI update)
- [ ] Generate publication-quality figures from existing data
- [ ] Formal sensitivity/specificity with bootstrap CIs from ViroForge ground truth
- [ ] Add BBDuk and/or Trimmomatic to the fastp comparison
- [ ] Write manuscript

## User Readiness
- [x] `virome-qc report` recomputes flags from passport data (single source of truth)
- [x] `virome-qc db --setup` auto-downloads SILVA + T2T and builds filters
- [x] Better error messages when filters are missing
- [x] Quickstart guide in README
- [x] Duplicate flag bug fixed
- [x] SE/PE ambiguity detection in ingestion engine
- [ ] Reports show expected ranges from profile (UI update)

## Platform Expansion
- [ ] Long-read QC modules (ONT adapter trimming, PacBio CCS quality filtering)
- [ ] Probe-capture virome support (separate modality, v2 — requires panel BED, alignment, on-target metrics)
- [ ] Population ERV database from HPRC pangenome
- [ ] Coordinate-based EVE flagging (ciHHV-6, recent HERV-K)
- [ ] Multi-sample comparison dashboard
- [ ] CI/CD regression testing with ViroForge canonical dataset

## Atlas / Community
- [ ] Publish ViroForge reference datasets + virome-qc passports as public resource
- [ ] Auto-classification: "your sample most closely matches X"
- [ ] Community-contributed sample type profiles
- [ ] Zenodo deposition for SILVA and T2T filter files

## HTCF Benchmarking
- [x] Tier 1: 13 core datasets complete (Tisza wastewater removed — probe-capture not supported)
- [x] Tier 1 evaluation: duplicate flags bug, SE/PE detection, probe-capture decision
- [ ] Tier 2: 10 extended datasets
- [ ] Tier 3: ~5,000 samples for atlas generation
