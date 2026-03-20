//! Platform detection from Illumina FASTQ headers
//!
//! Illumina instrument IDs encode the platform model:
//! ```text
//! @A00882:431:HVNLJDSXY:2:1101:1380:1000 1:N:0:ACGTACGT
//!  ^-- instrument ID
//! ```

use serde::{Deserialize, Serialize};

/// Sequencing chemistry type
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum Chemistry {
    /// 4-color SBS (MiSeq, HiSeq 2000/2500) -- no poly-G artifacts
    FourColor,
    /// 2-color SBS (NextSeq, NovaSeq) -- poly-G artifacts from no-signal = G
    TwoColor,
    /// Unknown or non-Illumina
    Unknown,
}

/// Detected platform information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlatformInfo {
    /// Raw instrument ID from FASTQ header
    pub instrument_id: String,
    /// Detected platform model
    pub model: String,
    /// Chemistry type (determines poly-G behavior)
    pub chemistry: Chemistry,
    /// Whether the flowcell uses patterned wells (affects optical duplicate distance)
    pub patterned_flowcell: bool,
    /// Expected Q-score binning (None = continuous, Some(3) = 3-bin NovaSeq)
    pub quality_bins: Option<u8>,
    /// Flowcell ID (if detected)
    pub flowcell_id: Option<String>,
}

impl PlatformInfo {
    /// Detect platform from a FASTQ header line
    ///
    /// Parses the instrument ID prefix to determine the platform model.
    /// Returns None if the header doesn't match Illumina format.
    pub fn from_header(header: &str) -> Option<Self> {
        let header = header.strip_prefix('@').unwrap_or(header);

        // Try two formats:
        // 1. Standard Illumina: @instrument:run:flowcell:lane:tile:x:y [comment]
        // 2. ENA-reformatted:   @ERR123456.1 instrument:run:flowcell:lane:tile:x:y/1
        //    (original header preserved in comment after space)

        // Split on space first to separate ID from comment
        let (primary, comment) = header.split_once(' ').unwrap_or((header, ""));

        // Try primary field first
        let primary_fields: Vec<&str> = primary.split(':').collect();
        let instrument_id = primary_fields[0];
        let (model, chemistry, patterned, q_bins) = detect_model(instrument_id);

        if chemistry != Chemistry::Unknown {
            // Primary field has a recognized instrument ID
            let flowcell_id = primary_fields.get(2).map(|s| s.to_string());
            return Some(PlatformInfo {
                instrument_id: instrument_id.to_string(),
                model,
                chemistry,
                patterned_flowcell: patterned,
                quality_bins: q_bins,
                flowcell_id,
            });
        }

        // Primary didn't match -- try the comment field (ENA-reformatted headers)
        if !comment.is_empty() {
            let comment_clean = comment.trim_end_matches("/1").trim_end_matches("/2");
            let comment_fields: Vec<&str> = comment_clean.split(':').collect();
            if comment_fields.len() >= 3 {
                let instrument_id = comment_fields[0];
                let (model, chemistry, patterned, q_bins) = detect_model(instrument_id);
                if chemistry != Chemistry::Unknown {
                    let flowcell_id = comment_fields.get(2).map(|s| s.to_string());
                    return Some(PlatformInfo {
                        instrument_id: instrument_id.to_string(),
                        model,
                        chemistry,
                        patterned_flowcell: patterned,
                        quality_bins: q_bins,
                        flowcell_id,
                    });
                }
            }
        }

        // Neither matched -- return with Unknown
        let flowcell_id = primary_fields.get(2).map(|s| s.to_string());
        Some(PlatformInfo {
            instrument_id: instrument_id.to_string(),
            model,
            chemistry,
            patterned_flowcell: patterned,
            quality_bins: q_bins,
            flowcell_id,
        })
    }

    /// Recommended optical duplicate pixel distance
    pub fn optical_dup_distance(&self) -> u32 {
        if self.patterned_flowcell {
            2500
        } else {
            100
        }
    }

    /// Whether aggressive poly-G trimming is warranted
    pub fn has_poly_g_risk(&self) -> bool {
        self.chemistry == Chemistry::TwoColor
    }
}

/// Detect platform model from instrument ID prefix
///
/// Returns (model_name, chemistry, patterned_flowcell, q_score_bins)
fn detect_model(instrument_id: &str) -> (String, Chemistry, bool, Option<u8>) {
    // Check prefix patterns (case-insensitive first 1-2 chars)
    let id_upper = instrument_id.to_uppercase();

    // NovaSeq X / X Plus (LH prefix)
    if id_upper.starts_with("LH") {
        return (
            "NovaSeq X".into(),
            Chemistry::TwoColor,
            true,
            Some(3),
        );
    }

    // NovaSeq 6000 (A prefix, typically A00xxx)
    if id_upper.starts_with("A0") {
        return (
            "NovaSeq 6000".into(),
            Chemistry::TwoColor,
            true,
            Some(3),
        );
    }

    // NextSeq 1000/2000 (VH prefix)
    if id_upper.starts_with("VH") {
        return (
            "NextSeq 2000".into(),
            Chemistry::TwoColor,
            true,
            Some(3),
        );
    }

    // NextSeq 500/550 (NB or NS prefix)
    if id_upper.starts_with("NB") || id_upper.starts_with("NS") {
        return (
            "NextSeq 500".into(),
            Chemistry::TwoColor,
            true,
            Some(3),
        );
    }

    // MiSeq (M prefix, typically M0xxxx)
    if id_upper.starts_with("M0") || id_upper.starts_with("M1") {
        return (
            "MiSeq".into(),
            Chemistry::FourColor,
            false,
            None,
        );
    }

    // HiSeq 3000/4000 (J prefix)
    if id_upper.starts_with('J') {
        return (
            "HiSeq 3000/4000".into(),
            Chemistry::FourColor,
            true,
            None,
        );
    }

    // HiSeq 2500 (D prefix in rapid mode, SN in high-output)
    if id_upper.starts_with('D') {
        return (
            "HiSeq 2500".into(),
            Chemistry::FourColor,
            false,
            None,
        );
    }

    // HiSeq 2000 (SN prefix)
    if id_upper.starts_with("SN") {
        return (
            "HiSeq 2000".into(),
            Chemistry::FourColor,
            false,
            None,
        );
    }

    // iSeq (FS prefix)
    if id_upper.starts_with("FS") {
        return (
            "iSeq 100".into(),
            Chemistry::TwoColor,
            true,
            Some(3),
        );
    }

    // Unknown
    (
        format!("Unknown ({})", instrument_id),
        Chemistry::Unknown,
        false,
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_novaseq_6000() {
        let header = "@A00882:431:HVNLJDSXY:2:1101:1380:1000 1:N:0:ACGTACGT";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "NovaSeq 6000");
        assert_eq!(info.chemistry, Chemistry::TwoColor);
        assert!(info.patterned_flowcell);
        assert!(info.has_poly_g_risk());
        assert_eq!(info.optical_dup_distance(), 2500);
        assert_eq!(info.flowcell_id, Some("HVNLJDSXY".into()));
    }

    #[test]
    fn test_miseq() {
        let header = "@M01757:7:000000000-AR185:1:1101:10000:19749/1";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "MiSeq");
        assert_eq!(info.chemistry, Chemistry::FourColor);
        assert!(!info.patterned_flowcell);
        assert!(!info.has_poly_g_risk());
        assert_eq!(info.optical_dup_distance(), 100);
    }

    #[test]
    fn test_nextseq_500() {
        let header = "@NB501138:21:HCTWLAFXX:1:11101:17138:1026/1";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "NextSeq 500");
        assert_eq!(info.chemistry, Chemistry::TwoColor);
        assert!(info.patterned_flowcell);
        assert!(info.has_poly_g_risk());
    }

    #[test]
    fn test_novaseq_x() {
        let header = "@LH00213:45:22FCNMLT3:1:1101:1000:1000 1:N:0:ACGT";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "NovaSeq X");
        assert_eq!(info.chemistry, Chemistry::TwoColor);
    }

    #[test]
    fn test_nextseq_2000() {
        let header = "@VH01234:10:AACJWN2M5:1:1101:1000:1000 1:N:0:ACGT";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "NextSeq 2000");
        assert_eq!(info.chemistry, Chemistry::TwoColor);
    }

    #[test]
    fn test_hiseq_2000() {
        let header = "@SN1234:100:ABCDEF:1:1101:1000:1000 1:N:0:ACGT";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "HiSeq 2000");
        assert_eq!(info.chemistry, Chemistry::FourColor);
    }

    #[test]
    fn test_unknown_platform() {
        let header = "@ZZZZ:1:FC:1:1:1:1";
        let info = PlatformInfo::from_header(header).unwrap();
        assert!(info.model.starts_with("Unknown"));
        assert_eq!(info.chemistry, Chemistry::Unknown);
    }

    #[test]
    fn test_ena_reformatted_miseq() {
        // ENA replaces instrument ID with run accession, preserves original in comment
        let header = "@ERR10359658.1 M01757:7:000000000-AR185:1:1101:10000:19749/1";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "MiSeq");
        assert_eq!(info.instrument_id, "M01757");
        assert_eq!(info.chemistry, Chemistry::FourColor);
    }

    #[test]
    fn test_ena_reformatted_nextseq() {
        let header = "@ERR1877475.1 NB501138:21:HCTWLAFXX:1:11101:17138:1026/1";
        let info = PlatformInfo::from_header(header).unwrap();
        assert_eq!(info.model, "NextSeq 500");
        assert_eq!(info.instrument_id, "NB501138");
        assert_eq!(info.chemistry, Chemistry::TwoColor);
        assert!(info.patterned_flowcell);
    }

    #[test]
    fn test_non_illumina_header() {
        // ONT header -- too few colon-separated fields
        let header = "@read_1234";
        let info = PlatformInfo::from_header(header);
        // Returns Some with Unknown model (not None) since we still parse what we can
        assert!(info.is_some());
        assert!(info.unwrap().model.starts_with("Unknown"));
    }
}
