# virome-qc TODO

## Paper Submission
- [ ] Final benchmark regeneration with latest code (all datasets, single code version)
- [ ] Evaluate and finalize benchmark dataset selection (see below)
- [ ] Embed ViroForge reference ranges into profiles
- [ ] Reports show user metrics vs expected range with color-coding
- [ ] Generate publication-quality figures from existing data
- [ ] Formal sensitivity/specificity with bootstrap CIs from ViroForge ground truth
- [ ] Add BBDuk and/or Trimmomatic to the fastp comparison
- [ ] Write manuscript

## User Readiness
- [ ] Auto-download mechanism for SILVA and T2T filters (`virome-qc db --setup`)
- [ ] Quickstart guide: install to first report in 5 minutes
- [ ] Better error messages when filters are missing
- [ ] `virome-qc report` command recomputes flags from passport data (not just re-renders HTML)

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
