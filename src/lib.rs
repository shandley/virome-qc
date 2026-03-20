//! virome-qc: High-performance, virome-specific QC platform
//!
//! Built on biometal primitives for parallelized read processing with
//! dual QA/QC at every step. Designed for the Human Virome Program.

pub mod config;
pub mod corpus;
pub mod ingest;
pub mod modules;
pub mod pipeline;
pub mod report;

pub use config::{Profile, ProfileConfig};
pub use corpus::{CorpusConfig, CorpusGenerator, SampleComposition};
pub use ingest::{IngestResult, ingest_fastq};
pub use pipeline::{Pipeline, PipelineResult};
pub use report::Passport;
