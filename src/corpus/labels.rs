//! Ground truth labels for synthetic reads

use serde::{Deserialize, Serialize};
use std::fmt;

/// A single label on a read
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum ReadLabel {
    /// Source organism/category
    Source(String),
    /// Specific virus name (if viral)
    Virus(String),
    /// Position in reference genome
    RefPosition(usize),
    /// 3' adapter contamination: (adapter_name, overlap_length)
    Adapter3Prime(String, usize),
    /// Internal adapter contamination at position
    AdapterInternal(usize),
    /// Random primer (first N bases are primer-derived)
    RandomPrimer(usize),
    /// Quality tail starts degrading at this position
    QualityTail(usize),
    /// Poly-X tail: (base, length)
    PolyX(char, usize),
    /// Complexity score of the true sequence
    Complexity(f64),
    /// This read is a PCR duplicate of another read
    PcrDuplicate(String),
    /// rRNA subunit
    Rrna(String),
    /// Read is from PhiX spike-in
    PhiX,
    /// Read is from host genome
    Host(String),
}

/// Collection of labels for a single read
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReadLabels {
    pub labels: Vec<ReadLabel>,
}

impl ReadLabels {
    pub fn new() -> Self {
        Self { labels: Vec::new() }
    }

    pub fn add(&mut self, label: ReadLabel) {
        self.labels.push(label);
    }

    /// Check if this read has a specific label type
    pub fn has_adapter_3prime(&self) -> bool {
        self.labels
            .iter()
            .any(|l| matches!(l, ReadLabel::Adapter3Prime(_, _)))
    }

    pub fn has_internal_adapter(&self) -> bool {
        self.labels
            .iter()
            .any(|l| matches!(l, ReadLabel::AdapterInternal(_)))
    }

    pub fn has_polyx(&self) -> bool {
        self.labels
            .iter()
            .any(|l| matches!(l, ReadLabel::PolyX(_, _)))
    }

    pub fn is_low_complexity(&self) -> bool {
        self.labels.iter().any(|l| match l {
            ReadLabel::Complexity(score) => *score < 0.3,
            _ => false,
        })
    }

    pub fn is_host(&self) -> bool {
        self.labels.iter().any(|l| matches!(l, ReadLabel::Host(_)))
    }

    pub fn is_rrna(&self) -> bool {
        self.labels.iter().any(|l| matches!(l, ReadLabel::Rrna(_)))
    }

    pub fn is_phix(&self) -> bool {
        self.labels.iter().any(|l| matches!(l, ReadLabel::PhiX))
    }

    pub fn source(&self) -> Option<&str> {
        self.labels.iter().find_map(|l| match l {
            ReadLabel::Source(s) => Some(s.as_str()),
            _ => None,
        })
    }

    /// Encode labels into FASTQ comment string
    pub fn to_comment(&self) -> String {
        let parts: Vec<String> = self
            .labels
            .iter()
            .map(|l| match l {
                ReadLabel::Source(s) => format!("source={s}"),
                ReadLabel::Virus(v) => format!("virus={v}"),
                ReadLabel::RefPosition(p) => format!("pos={p}"),
                ReadLabel::Adapter3Prime(name, len) => format!("adapter_3prime={name}:{len}"),
                ReadLabel::AdapterInternal(pos) => format!("adapter_internal={pos}"),
                ReadLabel::RandomPrimer(len) => format!("random_primer={len}"),
                ReadLabel::QualityTail(pos) => format!("quality_tail={pos}"),
                ReadLabel::PolyX(base, len) => format!("polyx={base}:{len}"),
                ReadLabel::Complexity(score) => format!("complexity={score:.3}"),
                ReadLabel::PcrDuplicate(of) => format!("pcr_dup={of}"),
                ReadLabel::Rrna(subunit) => format!("rrna={subunit}"),
                ReadLabel::PhiX => "phix=true".to_string(),
                ReadLabel::Host(org) => format!("host={org}"),
            })
            .collect();
        parts.join(";")
    }

    /// Parse labels from a FASTQ comment string
    pub fn from_comment(comment: &str) -> Self {
        let mut labels = ReadLabels::new();
        for part in comment.split(';') {
            let Some((key, value)) = part.split_once('=') else {
                continue;
            };
            match key {
                "source" => labels.add(ReadLabel::Source(value.to_string())),
                "virus" => labels.add(ReadLabel::Virus(value.to_string())),
                "pos" => {
                    if let Ok(p) = value.parse() {
                        labels.add(ReadLabel::RefPosition(p));
                    }
                }
                "adapter_3prime" => {
                    if let Some((name, len)) = value.split_once(':') {
                        if let Ok(l) = len.parse() {
                            labels.add(ReadLabel::Adapter3Prime(name.to_string(), l));
                        }
                    }
                }
                "adapter_internal" => {
                    if let Ok(pos) = value.parse() {
                        labels.add(ReadLabel::AdapterInternal(pos));
                    }
                }
                "random_primer" => {
                    if let Ok(len) = value.parse() {
                        labels.add(ReadLabel::RandomPrimer(len));
                    }
                }
                "quality_tail" => {
                    if let Ok(pos) = value.parse() {
                        labels.add(ReadLabel::QualityTail(pos));
                    }
                }
                "polyx" => {
                    if let Some((base, len)) = value.split_once(':') {
                        if let (Some(b), Ok(l)) = (base.chars().next(), len.parse()) {
                            labels.add(ReadLabel::PolyX(b, l));
                        }
                    }
                }
                "complexity" => {
                    if let Ok(score) = value.parse() {
                        labels.add(ReadLabel::Complexity(score));
                    }
                }
                "pcr_dup" => labels.add(ReadLabel::PcrDuplicate(value.to_string())),
                "rrna" => labels.add(ReadLabel::Rrna(value.to_string())),
                "phix" => labels.add(ReadLabel::PhiX),
                "host" => labels.add(ReadLabel::Host(value.to_string())),
                _ => {} // ignore unknown labels
            }
        }
        labels
    }
}

impl fmt::Display for ReadLabels {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_comment())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip() {
        let mut labels = ReadLabels::new();
        labels.add(ReadLabel::Source("viral".into()));
        labels.add(ReadLabel::Virus("HIV-1".into()));
        labels.add(ReadLabel::Adapter3Prime("truseq".into(), 15));
        labels.add(ReadLabel::Complexity(0.95));

        let comment = labels.to_comment();
        let parsed = ReadLabels::from_comment(&comment);

        assert_eq!(parsed.labels.len(), 4);
        assert_eq!(parsed.source(), Some("viral"));
        assert!(parsed.has_adapter_3prime());
    }

    #[test]
    fn test_flags() {
        let mut labels = ReadLabels::new();
        labels.add(ReadLabel::Host("human".into()));
        labels.add(ReadLabel::Rrna("16S".into()));

        assert!(labels.is_host());
        assert!(labels.is_rrna());
        assert!(!labels.is_phix());
        assert!(!labels.has_adapter_3prime());
    }
}
