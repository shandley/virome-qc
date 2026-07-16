# virome-qc review findings (2026-07-16)

Broad review requested by SH ("evaluate if there are other issues, big or small").
Three parallel audits (code wiring, doc consistency, load-bearing correctness) + empirical
validation of the new EVE screen. Nothing fixed yet -- this is the prioritized list to
choose from. Severity and, where noted, VERIFIED = I re-read the source and confirmed.

## Tier 0 -- publication integrity (fix before anything goes to a reviewer)

### A. [DONE 2026-07-16] The ERV endo/exo classifier is live and emits labels, but its "three signals" are one
RESOLVED: classifier decommissioned (code + docs). Removed run_erv_analysis, erv_classifier,
erv_pipeline, eve_polymorphisms, passport.erv_analysis, cpg_ratio, and the report-ui card
(bundle rebuilt); retroviral read detection retained; build/clippy clean, 123 tests pass.
Claims stripped from README, PUBLICATION_PLAN, VALIDATION (sec.17 + correlation), TODO.
Remaining in this tier: item B (the three correctness bugs) still open.
The classifier is NOT dead code: ErvScreener is pushed into every `run` (executor.rs:708),
gates run_erv_analysis (main.rs:750), and renders endo/exo in the report
(SampleReport.tsx:938). But MinHash is constant 0.5 (erv_classifier.rs:336), family is
hardcoded Retroviridae (erv_pipeline.rs:281), and polymorphism k-mers are empty
(eve_polymorphisms.rs:85). So CpG depletion silently does all the work -- and CpG is
directionally wrong for ciHHV-6 and confounded by ZAP in exogenous RNA viruses (per the
redesign). PUBLICATION_PLAN.md:82 calls it "diagnostic-grade"; README.md:36,170,283 and
VALIDATION.md sec.17 (662-690),:876 claim it as a validated first-in-class 3-signal
classifier. This is live-but-misleading output feeding manuscript claims.
Decommission cut points: main.rs:293 + executor.rs:708; relocate
eve_polymorphisms::match_cluster (runs inside run_erv_analysis, main.rs:869).

### B. [DONE 2026-07-16] Three correctness bugs that make the tool report WRONG QC numbers
FIXED + verified (build/clippy clean, 121 lib tests incl. 2 new dedup regression tests):
- qc_survival_rate: both copies (passport builder + recompute_flags) now count quality
  survivors as reads_passed + length-filter-only removals, so concordant-mate removals are
  correctly counted as failures instead of survivors.
- streaming dedup: added QcModule::process_pair (default = per-mate; StreamingDedup overrides
  to key on combined R1+R2 prefix); distinct fragments sharing one mate's start no longer
  collapse. Executor paired loop calls process_pair.
- HLL duplication: added hll_samples counter; every read offered to the sketch and the rate
  divides by reads actually added (robust to try_lock contention), not the full input.
CALIBRATION NOW STALE (open, publication-relevant): these fixes change reported survival and
duplication, and the pair-aware dedup lowers PE dedup_removed (fewer false duplicates -> higher
survival). So the profile expected_ranges and the VALIDATION.md survival/duplication tables are
now INCONSISTENT with the corrected code (i.e., wrong until re-derived). A rerun + re-derivation
is required before those numbers are used. Warning banner added at VALIDATION.md "Derived
Expected Ranges".

Original finding (all VERIFIED before fixing) below. These matter for publication because the
reported numbers feed the validation tables and the profile calibration ranges.

- **qc_survival_rate counts concordant-mate removals as survivors** (passport.rs:205-214).
  concordant_mate failures are applied in the executor (executor.rs:348,351), not by any
  module, so they are absent from qc_fail_non_dedup but remain in the denominator. 100%
  host-contaminated PE -> reports 50% survival (R1 fails "host" counted, R2 "concordant_mate"
  uncounted). Feeds LOW_QC_SURVIVAL flag -> can flip a fail into a pass. Worst in paired
  host-rich data (the clinically important regime). HIGH. VERIFIED.

- **Streaming dedup keys on a single mate, not the pair** (dedup_streaming.rs:48-69 +
  executor.rs:329-336). Design intent (doc lines 6-9, and dedup.rs) is "both mates identical
  to call a duplicate"; code hashes each mate's prefix alone against one shared set. Two
  distinct fragments from an abundant virus sharing an R1 prefix -> second dropped, then its
  R2 dropped via concordant logic (pcr_duplicate is in is_concordant_fail, executor.rs:344).
  Over-deduplicates exactly the natural viral duplicates the design meant to preserve.
  Enabled by default in all profiles. The doc comment is self-contradictory. HIGH. VERIFIED.

- **HLL duplication rate inflated: numerator subsampled, denominator not** (qa_stats.rs:
  308,404-406). Only every 8th read is added to the HLL; cardinality() returns the distinct
  count of what was added (~N/8), not rescaled; dup rate = 1 - unique/input divides by full
  input. A zero-duplicate library reports ~87.5% duplication. The comment (305-306) "HLL is
  accurate with subsampling" is the mental-model error. This is the analytics ESTIMATE;
  confirm whether it surfaces in the report and whether the profile duplication_rate ranges
  (VALIDATION.md ~0.26) came from this path or from the actual dedup module counts. HIGH
  (code); blast radius depends on surfacing. VERIFIED.

Consequence to check: if the survival / duplication numbers in VALIDATION.md and the profile
expected_ranges were produced by the buggy paths, they should be re-derived after the fix.

## Tier 1 -- reproducibility / consistency (a reviewer will notice)

- Correlation reported two ways for the same n=11 analysis: r=0.871 (p=0.001) in
  README.md:321 / PUB:43,69 vs r=0.922 (p=0.0004) in VALIDATION.md:865. Reconcile to one.
- Dataset count inconsistent: "27 / 12 real" (README:37, PUB:41) vs "28 / 13 real"
  (README:317). Reconcile.
- TODO.md:39-41 items marked [x] that are not done: "non-retroviral EVE detection" (the file
  holds exogenous refs, not host nrEVEs), "EVE k-mer exclusion" (approach rejected 2026-07-16),
  "polymorphism DB" (empty diagnostic_kmers scaffold).
- Profiles enable eve_aware/rescue by default (short-read-nebnext, stool-vlp-tagmentation
  :37-41) but the exclusion set is one K113 locus and the k-mer approach is superseded.
- Merge ambiguity uses a stale second_best (merge.rs:194-200): once set it never updates to a
  closer runner-up, so rates 0.10,0.03,0.04 merge a case that should be flagged ambiguous
  (potential chimera). MEDIUM. VERIFIED.
- qa_stats.summary.survival_rate disagrees with passport survival_rate for merged pairs
  (executor.rs:405 vs 410: reads_passed +=2 but record_passed once). Two conflicting survival
  numbers serialized into the same passport. MEDIUM/LOW.

## Tier 2 -- code health (low risk)

- Dead functions: erv_classifier.rs:200 longest_orf_six_frame (zero refs); :230, :111
  test-only.
- Two dedup implementations: dedup.rs (standalone command) vs dedup_streaming.rs (pipeline).
- Mislabeled data in the LIVE path (reached via #[path] submodules, so not orphaned):
  eve_kmers.rs (K113 vs the chr1 HERVK-int/K102 mislabel) read at eve_exclusion.rs:35;
  eve_sequences_nonretro.rs (exogenous set-D refs under a "non-retroviral EVE" name) read at
  erv.rs:149. eve_sequences_nonretro should move to the set-D panel.
- host.rs threshold is `>=` but the doc (host.rs:6-8) says ">50% / 15-50%"; a read at exactly
  0.50 containment is removed as host, not flagged ambiguous. Cosmetic. VERIFIED.

## Validated and worth preserving

- Retroviral read DETECTION (k-mer flagging) + herpesvirus exclusion: functional; only the
  endo/exo classification on top is the problem.
- New alignment EVE screen references: VALIDATED. K113 held-out recall 18/18 (100%) on both
  protein and nucleotide arms (a real absent-from-CHM13 probe, via HML-2 paralog pooling);
  0 false positives on random DNA and a real human non-HERV window. Provenance in
  references/protein_db_PROVENANCE.md.
- Tool builds clean; clippy clean (1 trivial warning).

## Cross-cutting

957 uncommitted insertions across 18 files + untracked modules/docs/references = two
entangled change-sets (dedup/survival-rate feature vs EVE redesign). Large unsaved surface;
wants splitting into two commits when ready. Not committed (SH has not asked).
