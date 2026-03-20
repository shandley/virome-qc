//! Synthetic test corpus generator for virome QC benchmarking
//!
//! Generates FASTQ files with known, planted artifacts for measuring
//! sensitivity and specificity of each QC module. Each read is tagged
//! in the FASTQ header with ground truth labels.
//!
//! # Header format
//!
//! ```text
//! @read_00001 source=viral;virus=phiX174;pos=1234;adapter_3prime=truseq:10;quality_tail=20;complexity=0.95
//! ```
//!
//! Labels are semicolon-separated key=value pairs in the comment field.
//! This enables exact precision/recall computation for every QC module.

mod generator;
mod labels;
mod sequences;

pub use generator::{CorpusConfig, CorpusGenerator, SampleComposition};
pub use labels::{ReadLabel, ReadLabels};
pub use sequences::ADAPTERS;
