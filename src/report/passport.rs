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
}

impl Passport {
    /// Create a passport from pipeline results, evaluating against profile thresholds
    pub fn from_result(result: &PipelineResult) -> Self {
        let survival_rate = result.survival_rate();
        let thresholds = &result.thresholds;

        let modules: Vec<ModuleReport> = result
            .module_reports
            .iter()
            .map(|(_, report)| report.clone())
            .collect();

        let mut flags = Vec::new();

        // Check survival rate against profile threshold
        if survival_rate < thresholds.min_survival_rate {
            flags.push(QcFlag {
                code: "LOW_SURVIVAL".into(),
                message: format!(
                    "Only {:.1}% of reads survived QC (threshold: {:.1}%)",
                    survival_rate * 100.0,
                    thresholds.min_survival_rate * 100.0
                ),
                severity: QualityTier::Fail,
            });
        } else if survival_rate < thresholds.min_survival_rate * 1.5 {
            // Warn if within 1.5x of the fail threshold
            flags.push(QcFlag {
                code: "LOW_SURVIVAL".into(),
                message: format!(
                    "{:.1}% read survival rate is approaching threshold ({:.1}%)",
                    survival_rate * 100.0,
                    thresholds.min_survival_rate * 100.0
                ),
                severity: QualityTier::Warn,
            });
        }

        // Check per-module metrics against thresholds
        Self::check_module_thresholds(&modules, result, thresholds, &mut flags);

        // Determine overall tier from flags
        let quality_tier = if flags.iter().any(|f| f.severity == QualityTier::Fail) {
            QualityTier::Fail
        } else if flags.iter().any(|f| f.severity == QualityTier::Warn) {
            QualityTier::Warn
        } else {
            QualityTier::Pass
        };

        Self {
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
        }
    }

    /// Check per-module metrics against profile thresholds
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
                    // Check internal adapter rate specifically
                    if let Some(internal) = report.extra.get("adapters_found_internal") {
                        if let Some(count) = internal.as_u64() {
                            let internal_rate = count as f64 / result.reads_input as f64;
                            if internal_rate > 0.01 {
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
                }
                "polyx" => {
                    // High poly-X rate suggests platform issue (NovaSeq/NextSeq)
                    let polyx_rate = report.reads_modified as f64 / result.reads_input as f64;
                    if polyx_rate > 0.15 {
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
                "host" => {
                    let host_fraction = report.reads_removed as f64 / result.reads_input as f64;
                    if host_fraction > thresholds.max_host_fraction {
                        flags.push(QcFlag {
                            code: "HIGH_HOST".into(),
                            message: format!(
                                "{:.1}% host reads exceeds threshold ({:.1}%)",
                                host_fraction * 100.0,
                                thresholds.max_host_fraction * 100.0
                            ),
                            severity: QualityTier::Warn,
                        });
                    }
                }
                // Future: "dedup" -> check dup rate
                _ => {}
            }
        }

        let _ = thresholds.max_duplicate_rate;
        let _ = thresholds.max_host_fraction;
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
