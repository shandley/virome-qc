//! Apply ingestion results as platform-aware overrides to QC profile
//!
//! Adjusts profile parameters based on detected platform characteristics
//! without changing the user's sample-type selections (adapter sets, thresholds).

use crate::config::ProfileConfig;
use crate::ingest::platform::Chemistry;
use crate::ingest::scan::IngestResult;

/// Apply platform-detected overrides to a profile configuration
///
/// Modifies parameters that are platform-dependent (not sample-dependent):
/// - Poly-G trimming aggressiveness based on chemistry
/// - Optical duplicate distance based on flowcell type
/// - Quality window size based on read length
/// - Minimum read length based on actual read length
/// - Complexity threshold based on read characteristics
pub fn apply_overrides(mut profile: ProfileConfig, ingest: &IngestResult) -> ProfileConfig {
    if let Some(ref platform) = ingest.platform {
        // Poly-G: enable platform-aware mode for 2-color chemistry
        if platform.chemistry == Chemistry::TwoColor {
            profile.modules.polyx.platform_aware = true;
        } else if platform.chemistry == Chemistry::FourColor {
            // 4-color: no poly-G artifacts, use standard threshold
            profile.modules.polyx.platform_aware = false;
        }

        // Optical duplicate distance
        if platform.patterned_flowcell {
            profile.modules.dedup.optical_distance = 2500;
        } else {
            profile.modules.dedup.optical_distance = 100;
        }
    }

    // Q-score binning detection
    if ingest.quick_scan.quality_binned {
        profile.modules.quality.quality_binned = true;
    }

    // Read length adjustments
    let read_length = ingest.reads.read_length;
    if read_length > 0 {
        // Min length: ~40% of read length, at least 30
        let auto_min_length = (read_length * 2 / 5).max(30);
        if profile.modules.quality.min_length > read_length {
            // Profile min_length exceeds actual read length -- must override
            profile.modules.quality.min_length = auto_min_length;
        }

        // Quality window: ~1/8 of read length, at least 4, at most 25
        if read_length <= 100 {
            profile.modules.quality.window_size = (read_length / 8).clamp(4, 25);
        }
    }

    // Adapter overlap threshold: high adapter rate suggests short inserts,
    // lower the min_overlap to catch more short adapter overlaps
    if let Some(ref adapter) = ingest.quick_scan.dominant_adapter {
        let rate = ingest
            .quick_scan
            .adapter_rates
            .get(adapter)
            .copied()
            .unwrap_or(0.0);
        if rate > 0.10 {
            // >10% adapter contamination -- inserts are short, lower overlap threshold
            profile.modules.adapter.min_overlap = 6;
        } else if rate > 0.03 {
            profile.modules.adapter.min_overlap = 8;
        }
    }

    // Auto-detect adapter sequences from scan data.
    // If ingestion found a dominant adapter family, override the profile's adapter config.
    if !ingest.quick_scan.detected_adapter_configs.is_empty() {
        let detected = &ingest.quick_scan.detected_adapter_configs;
        // Only override if the detected adapters are different from profile
        let profile_has_detected = detected
            .iter()
            .any(|d| profile.modules.adapter.sequences.iter().any(|s| s == d));
        if !profile_has_detected {
            // Replace adapter sequences with detected set
            profile.modules.adapter.sequences = detected.clone();
        }
    }

    // Data-driven complexity threshold from ingestion scan.
    // Use the 2nd percentile of observed complexity scores as the threshold:
    // this removes only the extreme low-complexity tail while preserving
    // reads that are normal for this specific sample.
    let p2 = ingest.quick_scan.complexity_p2;
    if p2 > 0.0 {
        // Clamp to a reasonable range: at least 0.2 (catch true junk),
        // at most 0.6 (don't over-filter)
        let data_driven_threshold = p2.clamp(0.2, 0.6);
        profile.modules.complexity.min_entropy = data_driven_threshold;
    }

    // Data-driven N-filter threshold: set to 5x the baseline N-rate,
    // with a floor of 0.05 and ceiling of 0.20. This accommodates
    // platforms with higher baseline N-rates (NextSeq ~2%) without
    // removing reads that are normal for the platform.
    let n_rate = ingest.quick_scan.n_rate;
    if n_rate > 0.005 {
        let auto_threshold = (n_rate * 5.0).clamp(0.05, 0.20);
        profile.modules.quality.max_n_fraction = auto_threshold;
    }

    profile
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Profile;
    use crate::ingest::platform::PlatformInfo;
    use crate::ingest::scan::{QuickScanStats, ReadInfo};
    use std::collections::HashMap;

    fn make_ingest(model: &str, chemistry: Chemistry, read_length: usize, gc: f64) -> IngestResult {
        IngestResult {
            platform: Some(PlatformInfo {
                instrument_id: "TEST".into(),
                model: model.into(),
                chemistry,
                patterned_flowcell: chemistry == Chemistry::TwoColor,
                quality_bins: if chemistry == Chemistry::TwoColor {
                    Some(3)
                } else {
                    None
                },
                flowcell_id: None,
            }),
            reads: ReadInfo {
                paired: true,
                read_length,
                variable_length: false,
                min_length: read_length,
                max_length: read_length,
                quality_offset: 33,
                estimated_read_count: Some(1_000_000),
            },
            quick_scan: QuickScanStats {
                reads_scanned: 10_000,
                mean_quality: 35.0,
                mean_gc: gc,
                n_rate: 0.001,
                adapter_rates: HashMap::new(),
                dominant_adapter: None,
                n_rate_pos0: 0.01,
                quality_binned: false,
                distinct_quality_values: 30,
                detected_adapter_configs: vec![],
                complexity_p2: 0.45,
                complexity_p5: 0.55,
                complexity_median: 0.85,
            },
            recommendations: vec![],
            warnings: vec![],
        }
    }

    #[test]
    fn test_novaseq_enables_poly_g() {
        let profile = Profile::load("stool-vlp-tagmentation").unwrap();
        let ingest = make_ingest("NovaSeq 6000", Chemistry::TwoColor, 150, 0.48);
        let modified = apply_overrides(profile, &ingest);

        assert!(modified.modules.polyx.platform_aware);
        assert_eq!(modified.modules.dedup.optical_distance, 2500);
    }

    #[test]
    fn test_miseq_disables_poly_g_platform_aware() {
        let mut profile = Profile::load("stool-vlp-tagmentation").unwrap();
        profile.modules.polyx.platform_aware = true; // start with enabled
        let ingest = make_ingest("MiSeq", Chemistry::FourColor, 250, 0.47);
        let modified = apply_overrides(profile, &ingest);

        assert!(!modified.modules.polyx.platform_aware);
        assert_eq!(modified.modules.dedup.optical_distance, 100);
    }

    #[test]
    fn test_short_read_min_length_override() {
        let profile = Profile::load("stool-vlp-tagmentation").unwrap();
        // 75bp reads with profile min_length=90 -- must override
        let ingest = make_ingest("NextSeq 500", Chemistry::TwoColor, 75, 0.50);
        let modified = apply_overrides(profile, &ingest);

        assert!(
            modified.modules.quality.min_length <= 75,
            "min_length {} should be <= read_length 75",
            modified.modules.quality.min_length
        );
    }

    #[test]
    fn test_data_driven_complexity_threshold() {
        let profile = Profile::load("stool-vlp-tagmentation").unwrap();
        let ingest = make_ingest("NovaSeq 6000", Chemistry::TwoColor, 150, 0.48);
        let modified = apply_overrides(profile, &ingest);

        // Data-driven threshold is clamped p2 value (0.45 in mock -> clamped to [0.2, 0.6])
        assert_eq!(modified.modules.complexity.min_entropy, 0.45);
    }

    #[test]
    fn test_very_low_complexity_sample_clamped() {
        let profile = Profile::load("stool-vlp-tagmentation").unwrap();
        let mut ingest = make_ingest("NovaSeq 6000", Chemistry::TwoColor, 150, 0.48);
        // Simulate a sample where p2 is very low (lots of low-complexity reads)
        ingest.quick_scan.complexity_p2 = 0.10;
        let modified = apply_overrides(profile, &ingest);

        // Should be clamped to 0.2 minimum -- don't pass obvious junk
        assert_eq!(modified.modules.complexity.min_entropy, 0.2);
    }
}
