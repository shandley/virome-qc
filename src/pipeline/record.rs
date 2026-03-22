//! Annotated read record that carries QA metrics through the pipeline

use biometal::FastqRecord;
use serde::{Deserialize, Serialize};

/// Disposition of a read after processing
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Disposition {
    /// Read passed QC
    Pass,
    /// Read failed QC with reason
    Fail(String),
    /// Read flagged but not removed (for downstream consideration)
    Flag(String),
}

/// Per-read QA metrics accumulated through the pipeline
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReadMetrics {
    /// Bases trimmed from 5' end (adapter, primer, quality)
    pub bases_trimmed_5prime: usize,
    /// Bases trimmed from 3' end (adapter, quality, poly-X)
    pub bases_trimmed_3prime: usize,
    /// Adapter type detected (if any)
    pub adapter_detected: Option<String>,
    /// Internal adapter found
    pub internal_adapter: bool,
    /// Poly-X bases trimmed
    pub polyx_trimmed: usize,
    /// Shannon entropy score
    pub complexity: Option<f64>,
    /// Mean quality after trimming
    pub mean_quality: Option<f64>,
}

/// A FASTQ record annotated with QA tags that flows through the pipeline
#[derive(Debug, Clone)]
pub struct AnnotatedRecord {
    /// The FASTQ record (modified in-place by modules)
    pub record: FastqRecord,
    /// Current disposition
    pub disposition: Disposition,
    /// Accumulated QA metrics
    pub metrics: ReadMetrics,
    /// Original read length before any trimming
    pub original_length: usize,
}

impl AnnotatedRecord {
    /// Create a new annotated record from a raw FASTQ record
    pub fn new(record: FastqRecord) -> Self {
        let original_length = record.sequence.len();
        Self {
            record,
            disposition: Disposition::Pass,
            metrics: ReadMetrics::default(),
            original_length,
        }
    }

    /// Check if the read has been failed
    pub fn is_failed(&self) -> bool {
        matches!(self.disposition, Disposition::Fail(_))
    }

    /// Check if the read has been flagged (ambiguous)
    pub fn is_flagged(&self) -> bool {
        matches!(self.disposition, Disposition::Flag(_))
    }

    /// Current sequence length
    pub fn len(&self) -> usize {
        self.record.sequence.len()
    }

    /// Check if the record is empty
    pub fn is_empty(&self) -> bool {
        self.record.sequence.is_empty()
    }

    /// Mark this read as failed with a reason
    pub fn fail(&mut self, reason: impl Into<String>) {
        self.disposition = Disposition::Fail(reason.into());
    }

    /// Flag this read (but don't remove it)
    pub fn flag(&mut self, reason: impl Into<String>) {
        if matches!(self.disposition, Disposition::Pass) {
            self.disposition = Disposition::Flag(reason.into());
        }
    }

    /// Total bases trimmed from this read
    pub fn total_trimmed(&self) -> usize {
        self.metrics.bases_trimmed_5prime + self.metrics.bases_trimmed_3prime
    }
}
