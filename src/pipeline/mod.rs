//! Pipeline orchestration and parallel execution

mod executor;
mod record;

pub use executor::{Pipeline, PipelineResult};
pub use record::{AnnotatedRecord, Disposition, ReadMetrics};
