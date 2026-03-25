//! QA Passport — machine-readable per-sample QC report
//!
//! Uses profile-defined thresholds to determine PASS/WARN/FAIL for each metric.

use crate::config::Thresholds;
use crate::modules::{AnalyticsSnapshot, ModuleReport};
use crate::pipeline::PipelineResult;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::Path;

fn is_zero(val: &u64) -> bool {
    *val == 0
}

/// Quality tier for the sample
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "UPPERCASE")]
pub enum QualityTier {
    Pass,
    Warn,
    Fail,
}

/// A flag raised during QC
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcFlag {
    pub code: String,
    pub message: String,
    pub severity: QualityTier,
}

/// Input file provenance information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Provenance {
    pub timestamp: String,
    pub input_files: Vec<InputFileInfo>,
}

/// Information about an input file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputFileInfo {
    pub path: String,
    pub size_bytes: u64,
}

/// QA Passport — the primary structured output of virome-qc
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Passport {
    /// virome-qc version
    pub tool_version: String,
    /// Profile used
    pub profile: String,
    /// Total input reads
    pub reads_input: u64,
    /// Total reads passing QC
    pub reads_passed: u64,
    /// Overall survival rate
    pub survival_rate: f64,
    /// Pairs where both mates passed (paired-end only)
    #[serde(default, skip_serializing_if = "is_zero")]
    pub pairs_passed: u64,
    /// Singleton reads (mate failed QC, paired-end only)
    #[serde(default, skip_serializing_if = "is_zero")]
    pub singletons: u64,
    /// Pairs successfully merged (paired-end only)
    #[serde(default, skip_serializing_if = "is_zero")]
    pub pairs_merged: u64,
    /// Per-module reports
    pub modules: Vec<ModuleReport>,
    /// Quality flags raised
    pub flags: Vec<QcFlag>,
    /// Overall quality tier
    pub quality_tier: QualityTier,
    /// QA distribution statistics (read length, GC, insert size)
    /// Comprehensive read analytics (per-position, distributions, duplication)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qa_stats: Option<AnalyticsSnapshot>,
    /// Input file provenance for reproducibility
    #[serde(skip_serializing_if = "Option::is_none")]
    pub provenance: Option<Provenance>,
    /// Ingestion scan results (platform detection, sample characterization)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ingestion: Option<serde_json::Value>,
    /// Applied module parameters (after all overrides, with source annotations)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub parameters: Option<serde_json::Value>,
    /// ERV analysis results (post-pipeline retroviral classification)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub erv_analysis: Option<serde_json::Value>,
    /// Contamination summary: aggregates biological removal modules
    #[serde(skip_serializing_if = "Option::is_none")]
    pub contamination_summary: Option<serde_json::Value>,
}

impl Passport {
    /// Create a passport from pipeline results, evaluating against profile thresholds
    pub fn from_result(result: &PipelineResult) -> Self {
        let survival_rate = result.survival_rate();
        let _thresholds = &result.thresholds; // Used via parameters JSON in recompute_flags()

        let modules: Vec<ModuleReport> = result
            .module_reports
            .iter()
            .map(|(_, report)| report.clone())
            .collect();

        // Placeholder — recompute_flags() will set the real flags and tier
        let flags: Vec<QcFlag> = Vec::new();
        let quality_tier = if flags.iter().any(|f| f.severity == QualityTier::Fail) {
            QualityTier::Fail
        } else if flags.iter().any(|f| f.severity == QualityTier::Warn) {
            QualityTier::Warn
        } else {
            QualityTier::Pass
        };

        // Build ingestion section for passport
        let ingestion_json = result.ingestion.as_ref().map(|ingest| {
            serde_json::json!({
                "reads_scanned": ingest.quick_scan.reads_scanned,
                "platform": ingest.platform.as_ref().map(|p| serde_json::json!({
                    "model": p.model,
                    "chemistry": format!("{:?}", p.chemistry),
                    "patterned_flowcell": p.patterned_flowcell,
                })),
                "read_length": ingest.reads.read_length,
                "paired": ingest.reads.paired,
                "mean_quality": ingest.quick_scan.mean_quality,
                "mean_gc": ingest.quick_scan.mean_gc,
                "n_rate": ingest.quick_scan.n_rate,
                "complexity_p2": ingest.quick_scan.complexity_p2,
                "quality_binned": ingest.quick_scan.quality_binned,
                "dominant_adapter": ingest.quick_scan.dominant_adapter,
                "adapter_3prime_rate": ingest.quick_scan.adapter_3prime_rate,
                "adapter_internal_rate": ingest.quick_scan.adapter_internal_rate,
                "r2_quality_delta": ingest.quick_scan.r2_quality_delta,
                "quality_profile": ingest.quick_scan.quality_profile,
            })
        });

        // Build applied parameters section
        let parameters_json = result.applied_config.as_ref().map(|cfg| {
            serde_json::json!({
                "adapter": {
                    "min_overlap": cfg.modules.adapter.min_overlap,
                    "max_mismatch_rate": cfg.modules.adapter.max_mismatch_rate,
                    "internal_scan": cfg.modules.adapter.internal_scan,
                    "sequences": cfg.modules.adapter.sequences,
                },
                "quality": {
                    "window_size": cfg.modules.quality.window_size,
                    "min_quality": cfg.modules.quality.min_quality,
                    "min_mean_quality": cfg.modules.quality.min_mean_quality,
                    "min_length": cfg.modules.quality.min_length,
                    "quality_binned": cfg.modules.quality.quality_binned,
                },
                "complexity": {
                    "min_entropy": cfg.modules.complexity.min_entropy,
                },
                "polyx": {
                    "platform_aware": cfg.modules.polyx.platform_aware,
                    "min_length": cfg.modules.polyx.min_length,
                },
                "contaminant": {
                    "min_kmer_fraction": cfg.modules.contaminant.min_kmer_fraction,
                    "screen_rrna": cfg.modules.contaminant.screen_rrna,
                    "screen_phix": cfg.modules.contaminant.screen_phix,
                    "screen_vectors": cfg.modules.contaminant.screen_vectors,
                },
                "host": {
                    "enabled": cfg.modules.host.enabled,
                    "host_threshold": cfg.modules.host.host_threshold,
                    "ambiguous_threshold": cfg.modules.host.ambiguous_threshold,
                },
                "dedup": {
                    "enabled": cfg.modules.dedup.enabled,
                },
                "quality_assessment": {
                    "min_survival_rate": cfg.thresholds.min_survival_rate,
                    "expected_gc_range": cfg.thresholds.expected_gc_range,
                    "max_host_fraction": cfg.thresholds.max_host_fraction,
                    "max_rrna_fraction": cfg.thresholds.max_rrna_fraction,
                    "max_duplicate_rate": cfg.thresholds.max_duplicate_rate,
                    "expected_ranges": cfg.thresholds.expected_ranges,
                },
            })
        });

        // Compute contamination summary before moving modules
        let contamination_summary = {
            let biological_modules = ["contaminant", "rrna", "host"];
            let technical_modules = [
                "adapter", "polyx", "n_filter", "quality", "complexity", "dedup",
            ];

            let bio_removed: u64 = modules
                .iter()
                .filter(|m| biological_modules.contains(&m.name.as_str()))
                .map(|m| m.reads_removed)
                .sum();

            let tech_removed: u64 = modules
                .iter()
                .filter(|m| technical_modules.contains(&m.name.as_str()))
                .map(|m| m.reads_removed)
                .sum();

            let ri = result.reads_input as f64;
            Some(serde_json::json!({
                "biological_contamination_removed": bio_removed,
                "biological_contamination_fraction": if ri > 0.0 { bio_removed as f64 / ri } else { 0.0 },
                "technical_artifacts_removed": tech_removed,
                "technical_artifacts_fraction": if ri > 0.0 { tech_removed as f64 / ri } else { 0.0 },
                "total_removed": bio_removed + tech_removed,
                "total_removed_fraction": if ri > 0.0 { (bio_removed + tech_removed) as f64 / ri } else { 0.0 },
            }))
        };

        let mut passport = Self {
            tool_version: env!("CARGO_PKG_VERSION").to_string(),
            profile: result.profile_name.clone(),
            reads_input: result.reads_input,
            reads_passed: result.reads_passed,
            survival_rate,
            pairs_passed: result.pairs_passed,
            singletons: result.singletons,
            pairs_merged: result.pairs_merged,
            modules,
            flags,
            quality_tier,
            qa_stats: result.qa_stats.clone(),
            provenance: result.provenance.clone(),
            ingestion: ingestion_json,
            parameters: parameters_json,
            erv_analysis: None,
            contamination_summary,
        };

        // Single source of truth for all flags — recompute from module data
        passport.recompute_flags();
        passport
    }

    /// Legacy: Check per-module metrics against profile thresholds
    /// Superseded by recompute_flags() which uses passport's own embedded data.
    #[allow(dead_code)]
    fn check_module_thresholds(
        modules: &[ModuleReport],
        result: &PipelineResult,
        thresholds: &Thresholds,
        flags: &mut Vec<QcFlag>,
    ) {
        if result.reads_input == 0 {
            return;
        }

        for report in modules {
            let removal_rate = report.reads_removed as f64 / result.reads_input as f64;

            match report.name.as_str() {
                "complexity" => {
                    // High complexity removal suggests library issues
                    if removal_rate > 0.10 {
                        flags.push(QcFlag {
                            code: "HIGH_LOW_COMPLEXITY".into(),
                            message: format!(
                                "{:.1}% of reads failed complexity filter",
                                removal_rate * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                "adapter" => {
                    // Check internal adapter rate
                    // Nextera/tagmentation preps inherently produce ~1-2% chimeras,
                    // so use a higher threshold (2.5%) than ligation preps (1%)
                    if let Some(internal) = report.extra.get("adapters_found_internal") {
                        if let Some(count) = internal.as_u64() {
                            let internal_rate = count as f64 / result.reads_input as f64;
                            if internal_rate > 0.025 {
                                flags.push(QcFlag {
                                    code: "HIGH_INTERNAL_ADAPTER".into(),
                                    message: format!(
                                        "{:.2}% of reads have internal adapter contamination",
                                        internal_rate * 100.0
                                    ),
                                    severity: QualityTier::Warn,
                                });
                            }
                        }
                    }

                    // Check internal adapter contamination rate
                    if let Some(internal) = report.extra.get("adapters_found_internal") {
                        if let Some(count) = internal.as_u64() {
                            let internal_rate = count as f64 / result.reads_input as f64;
                            if internal_rate > 0.05 {
                                flags.push(QcFlag {
                                    code: "HIGH_INTERNAL_ADAPTER".into(),
                                    message: format!(
                                        "{:.1}% of reads contain internal adapter chimeras -- possible library quality issue",
                                        internal_rate * 100.0
                                    ),
                                    severity: if internal_rate > 0.20 {
                                        QualityTier::Fail
                                    } else {
                                        QualityTier::Warn
                                    },
                                });
                            }
                        }
                    }

                    // Check 3' adapter trimming rate -- high rates indicate short inserts
                    if let Some(trimmed) = report.extra.get("adapters_found_3prime") {
                        if let Some(count) = trimmed.as_u64() {
                            let trim_rate = count as f64 / result.reads_input as f64;
                            if trim_rate > 0.20 {
                                flags.push(QcFlag {
                                    code: "HIGH_ADAPTER_TRIMMING".into(),
                                    message: format!(
                                        "{:.1}% of reads had adapter trimmed (short inserts)",
                                        trim_rate * 100.0
                                    ),
                                    severity: QualityTier::Warn,
                                });
                            }
                        }
                    }
                }
                "polyx" => {
                    // High poly-X rate suggests platform artifact or fragmented library
                    let polyx_rate = report.reads_modified as f64 / result.reads_input as f64;
                    if polyx_rate > 0.05 {
                        flags.push(QcFlag {
                            code: "HIGH_POLYX".into(),
                            message: format!(
                                "{:.1}% of reads had poly-X tails trimmed",
                                polyx_rate * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                "contaminant" => {
                    // Check rRNA fraction against threshold
                    if let Some(rrna) = report.extra.get("rrna_removed") {
                        if let Some(rrna_count) = rrna.as_u64() {
                            let rrna_fraction = rrna_count as f64 / result.reads_input as f64;
                            if rrna_fraction > thresholds.max_rrna_fraction {
                                flags.push(QcFlag {
                                    code: "HIGH_RRNA".into(),
                                    message: format!(
                                        "{:.1}% rRNA contamination exceeds threshold ({:.1}%)",
                                        rrna_fraction * 100.0,
                                        thresholds.max_rrna_fraction * 100.0
                                    ),
                                    severity: QualityTier::Warn,
                                });
                            }
                        }
                    }
                    // Check PhiX fraction
                    if let Some(phix) = report.extra.get("phix_removed") {
                        if let Some(phix_count) = phix.as_u64() {
                            let phix_fraction = phix_count as f64 / result.reads_input as f64;
                            if phix_fraction > 0.05 {
                                flags.push(QcFlag {
                                    code: "HIGH_PHIX".into(),
                                    message: format!("{:.2}% PhiX detected", phix_fraction * 100.0),
                                    severity: QualityTier::Warn,
                                });
                            }
                        }
                    }
                }
                "rrna" => {
                    // Check SILVA rRNA filter results (separate from contaminant consensus rRNA)
                    let rrna_fraction = report.reads_removed as f64 / result.reads_input as f64;
                    if rrna_fraction > thresholds.max_rrna_fraction {
                        let context = if rrna_fraction > 0.50 {
                            " -- sample may lack rRNA depletion or VLP enrichment"
                        } else if rrna_fraction > 0.10 {
                            " -- consider whether rRNA depletion was incomplete or library prep captures small rRNA fragments"
                        } else {
                            ""
                        };
                        flags.push(QcFlag {
                            code: "HIGH_RRNA".into(),
                            message: format!(
                                "{:.1}% rRNA detected by SILVA filter (threshold: {:.1}%){}",
                                rrna_fraction * 100.0,
                                thresholds.max_rrna_fraction * 100.0,
                                context,
                            ),
                            severity: if rrna_fraction > 0.50 {
                                QualityTier::Fail
                            } else {
                                QualityTier::Warn
                            },
                        });
                    }
                }
                "host" => {
                    let host_fraction = report.reads_removed as f64 / result.reads_input as f64;
                    if host_fraction > thresholds.max_host_fraction {
                        flags.push(QcFlag {
                            code: "HIGH_HOST".into(),
                            message: format!(
                                "{:.1}% host reads exceeds threshold ({:.1}%) -- consider VLP enrichment or host depletion upstream",
                                host_fraction * 100.0,
                                thresholds.max_host_fraction * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                "dedup" => {
                    let dup_rate = report.reads_removed as f64 / result.reads_input as f64;
                    if dup_rate > thresholds.max_duplicate_rate {
                        flags.push(QcFlag {
                            code: "HIGH_DUPLICATION".into(),
                            message: format!(
                                "{:.1}% of reads removed as duplicates (threshold: {:.0}%)",
                                dup_rate * 100.0,
                                thresholds.max_duplicate_rate * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                _ => {}
            }
        }

        let _ = thresholds.max_host_fraction;
    }

    /// Recompute flags from the passport's own module data and embedded thresholds.
    /// Used by `virome-qc report` to ensure flags reflect the latest flag logic
    /// even for passports generated by older code versions.
    pub fn recompute_flags(&mut self) {
        let ri = self.reads_input as f64;
        if ri == 0.0 {
            return;
        }

        // Extract thresholds from embedded parameters
        let qa = self
            .parameters
            .as_ref()
            .and_then(|p| p.get("quality_assessment"));
        let min_survival = qa
            .and_then(|q| q.get("min_survival_rate"))
            .and_then(|v| v.as_f64())
            .unwrap_or(0.30);
        let max_rrna = qa
            .and_then(|q| q.get("max_rrna_fraction"))
            .and_then(|v| v.as_f64())
            .unwrap_or(0.05);
        let max_host = qa
            .and_then(|q| q.get("max_host_fraction"))
            .and_then(|v| v.as_f64())
            .unwrap_or(0.20);
        let max_dup = qa
            .and_then(|q| q.get("max_duplicate_rate"))
            .and_then(|v| v.as_f64())
            .unwrap_or(0.50);
        let gc_range: Option<(f64, f64)> = qa.and_then(|q| {
            let arr = q.get("expected_gc_range")?.as_array()?;
            Some((arr.first()?.as_f64()?, arr.get(1)?.as_f64()?))
        });

        let mut flags = Vec::new();

        // Survival
        if self.survival_rate < min_survival * 1.5 {
            let dominant = self
                .modules
                .iter()
                .max_by_key(|m| m.reads_removed)
                .map(|m| {
                    let pct = m.reads_removed as f64 / ri * 100.0;
                    format!(" -- largest contributor: {} ({:.1}%)", m.name, pct)
                })
                .unwrap_or_default();

            let is_fail = self.survival_rate < min_survival;
            flags.push(QcFlag {
                code: "LOW_SURVIVAL".into(),
                message: format!(
                    "{}{:.1}% of reads survived QC (threshold: {:.1}%){}",
                    if is_fail { "Only " } else { "" },
                    self.survival_rate * 100.0,
                    min_survival * 100.0,
                    dominant,
                ),
                severity: if is_fail {
                    QualityTier::Fail
                } else {
                    QualityTier::Warn
                },
            });
        }

        // GC deviation
        if let Some((gc_min, gc_max)) = gc_range {
            if let Some(ref ing) = self.ingestion {
                if let Some(gc) = ing.get("mean_gc").and_then(|v| v.as_f64()) {
                    if gc < gc_min || gc > gc_max {
                        flags.push(QcFlag {
                            code: "GC_DEVIATION".into(),
                            message: format!(
                                "Mean GC {:.1}% is outside expected range ({:.0}-{:.0}%) for this profile",
                                gc * 100.0,
                                gc_min * 100.0,
                                gc_max * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
            }
        }

        // Per-module flags
        for m in &self.modules {
            let removal_rate = m.reads_removed as f64 / ri;
            match m.name.as_str() {
                "adapter" => {
                    let internal = m
                        .extra
                        .get("adapters_found_internal")
                        .and_then(|v| v.as_u64())
                        .unwrap_or(0);
                    let internal_rate = internal as f64 / ri;
                    if internal_rate > 0.05 {
                        flags.push(QcFlag {
                            code: "HIGH_INTERNAL_ADAPTER".into(),
                            message: format!(
                                "{:.1}% of reads contain internal adapter chimeras -- possible library quality issue",
                                internal_rate * 100.0
                            ),
                            severity: if internal_rate > 0.20 {
                                QualityTier::Fail
                            } else {
                                QualityTier::Warn
                            },
                        });
                    }
                    let trimmed = m
                        .extra
                        .get("adapters_found_3prime")
                        .and_then(|v| v.as_u64())
                        .unwrap_or(0);
                    if trimmed as f64 / ri > 0.20 {
                        flags.push(QcFlag {
                            code: "HIGH_ADAPTER_TRIMMING".into(),
                            message: format!(
                                "{:.1}% of reads had adapter trimmed (short inserts)",
                                trimmed as f64 / ri * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                "contaminant" => {
                    if let Some(rrna_count) = m.extra.get("rrna_removed").and_then(|v| v.as_u64())
                    {
                        let frac = rrna_count as f64 / ri;
                        if frac > max_rrna {
                            flags.push(QcFlag {
                                code: "HIGH_RRNA".into(),
                                message: format!(
                                    "{:.1}% rRNA contamination exceeds threshold ({:.1}%)",
                                    frac * 100.0,
                                    max_rrna * 100.0
                                ),
                                severity: QualityTier::Warn,
                            });
                        }
                    }
                    if let Some(phix_count) = m.extra.get("phix_removed").and_then(|v| v.as_u64())
                    {
                        if phix_count as f64 / ri > 0.05 {
                            flags.push(QcFlag {
                                code: "HIGH_PHIX".into(),
                                message: format!(
                                    "{:.2}% PhiX detected",
                                    phix_count as f64 / ri * 100.0
                                ),
                                severity: QualityTier::Warn,
                            });
                        }
                    }
                }
                "rrna" => {
                    if removal_rate > max_rrna {
                        let context = if removal_rate > 0.50 {
                            " -- sample may lack rRNA depletion or VLP enrichment"
                        } else if removal_rate > 0.10 {
                            " -- consider whether rRNA depletion was incomplete or library prep captures small rRNA fragments"
                        } else {
                            ""
                        };
                        flags.push(QcFlag {
                            code: "HIGH_RRNA".into(),
                            message: format!(
                                "{:.1}% rRNA detected by SILVA filter (threshold: {:.1}%){}",
                                removal_rate * 100.0,
                                max_rrna * 100.0,
                                context,
                            ),
                            severity: if removal_rate > 0.50 {
                                QualityTier::Fail
                            } else {
                                QualityTier::Warn
                            },
                        });
                    }
                }
                "host" => {
                    if removal_rate > max_host {
                        flags.push(QcFlag {
                            code: "HIGH_HOST".into(),
                            message: format!(
                                "{:.1}% host reads exceeds threshold ({:.1}%) -- consider VLP enrichment or host depletion upstream",
                                removal_rate * 100.0,
                                max_host * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                "dedup" => {
                    if removal_rate > max_dup {
                        flags.push(QcFlag {
                            code: "HIGH_DUPLICATION".into(),
                            message: format!(
                                "{:.1}% of reads removed as duplicates (threshold: {:.0}%)",
                                removal_rate * 100.0,
                                max_dup * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                "complexity" => {
                    if removal_rate > 0.10 {
                        flags.push(QcFlag {
                            code: "HIGH_LOW_COMPLEXITY".into(),
                            message: format!(
                                "{:.1}% of reads failed complexity filter",
                                removal_rate * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                _ => {}
            }
        }

        // Update quality tier
        let has_fail = flags.iter().any(|f| f.severity == QualityTier::Fail);
        let has_warn = flags.iter().any(|f| f.severity == QualityTier::Warn);
        self.quality_tier = if has_fail {
            QualityTier::Fail
        } else if has_warn {
            QualityTier::Warn
        } else {
            QualityTier::Pass
        };
        self.flags = flags;
    }

    /// Write passport as JSON
    pub fn write_json(&self, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(path, json)?;
        Ok(())
    }

    /// Write passport as YAML
    pub fn write_yaml(&self, path: &Path) -> Result<()> {
        let yaml = serde_yaml::to_string(self)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }
}
