//! Configuration and profile system
//!
//! Profiles define which QC modules run, in what order, with what parameters.
//! Each profile targets a specific sample type + library prep combination.

mod profiles;

pub use profiles::{
    AdapterConfig, ChimeraConfig, ComplexityConfig, ContaminantConfig, DedupConfig,
    ExpectedRanges, HostConfig, ModuleConfigs, Platform, PolyXConfig, Profile, ProfileConfig,
    QualityConfig, RrnaConfig, Thresholds,
};
