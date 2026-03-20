//! Ingestion engine -- fast pre-scan of FASTQ data to detect platform,
//! read characteristics, and potential issues before QC processing.
//!
//! Reads the first N records (default 10,000) to determine:
//! - Sequencing platform from instrument ID in FASTQ headers
//! - Read length, quality encoding, paired-end layout
//! - Quick-scan statistics (GC, N-rate, adapter contamination, duplication)
//! - Platform-specific recommendations (poly-G, optical dup distance, Q-binning)

pub mod markers;
mod overrides;
mod platform;
mod scan;

pub use markers::{CompositionEstimate, MarkerDb};
pub use overrides::apply_overrides;
pub use platform::PlatformInfo;
pub use scan::{IngestResult, QuickScanStats, ReadInfo, Recommendation, ingest_fastq};
