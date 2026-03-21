//! Profile definitions and loading

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// Sequencing platform
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Platform {
    Illumina,
    Ont,
    Pacbio,
}

/// Top-level profile configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProfileConfig {
    pub name: String,
    pub description: String,
    pub platform: Platform,
    pub modules: ModuleConfigs,
    pub thresholds: Thresholds,
}

/// Module-level configurations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModuleConfigs {
    pub adapter: AdapterConfig,
    pub quality: QualityConfig,
    pub polyx: PolyXConfig,
    pub complexity: ComplexityConfig,
    pub contaminant: ContaminantConfig,
    pub host: HostConfig,
    pub dedup: DedupConfig,
    pub chimera: ChimeraConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdapterConfig {
    pub enabled: bool,
    /// Adapter sequence sets to check against
    pub sequences: Vec<String>,
    /// Scan for internal adapter contamination (not just 3')
    pub internal_scan: bool,
    /// Number of bases to trim from 5' end for random primer removal (0 = disabled)
    pub random_primer_trim: usize,
    /// Minimum overlap length for 3' adapter detection
    #[serde(default = "default_min_overlap")]
    pub min_overlap: usize,
    /// Maximum mismatch rate for adapter matching
    #[serde(default = "default_max_mismatch_rate")]
    pub max_mismatch_rate: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityConfig {
    pub enabled: bool,
    /// Sliding window size for quality trimming
    pub window_size: usize,
    /// Minimum quality score within window
    pub min_quality: u8,
    /// Minimum mean quality for entire read
    pub min_mean_quality: f64,
    /// Minimum read length after trimming
    pub min_length: usize,
    /// Quality scores are binned (NovaSeq/NextSeq). Auto-set by ingestion engine.
    /// When true, uses bin-aware trimming (fraction of low-bin bases) instead of
    /// mean quality in sliding window.
    #[serde(default)]
    pub quality_binned: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolyXConfig {
    pub enabled: bool,
    /// Auto-detect platform from FASTQ headers for poly-G behavior
    pub platform_aware: bool,
    /// Minimum homopolymer run length to trigger trimming
    #[serde(default = "default_min_polyx_len")]
    pub min_length: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexityConfig {
    pub enabled: bool,
    /// Minimum Shannon entropy (0.0 - 1.0)
    pub min_entropy: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContaminantConfig {
    pub enabled: bool,
    pub screen_rrna: bool,
    pub screen_phix: bool,
    pub screen_vectors: bool,
    pub screen_kitome: bool,
    /// Minimum fraction of read k-mers matching contaminant index to classify as contaminant
    #[serde(default = "default_contaminant_threshold")]
    pub min_kmer_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HostConfig {
    pub enabled: bool,
    /// Host reference name (e.g., "human", "mouse") or path to FASTA
    pub reference: String,
    /// Use EVE-masked reference to avoid scrubbing viral sequences
    pub eve_aware: bool,
    /// Rescue reads that hit EVE regions by re-checking against viral DB
    pub rescue: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DedupConfig {
    pub enabled: bool,
    /// Optical duplicate pixel distance (100 for non-patterned, 2500 for patterned flow cells)
    pub optical_distance: u32,
    /// Use UMIs from read names for deduplication
    pub umi_aware: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChimeraConfig {
    pub enabled: bool,
}

/// QA thresholds for PASS/WARN/FAIL determination
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Thresholds {
    /// Minimum fraction of reads surviving QC (below = FAIL)
    pub min_survival_rate: f64,
    /// Maximum host read fraction (above = WARN)
    pub max_host_fraction: f64,
    /// Maximum rRNA fraction (above = WARN)
    pub max_rrna_fraction: f64,
    /// Maximum duplicate rate (above = WARN)
    pub max_duplicate_rate: f64,
}

// Default value functions for serde
fn default_min_overlap() -> usize { 10 }
fn default_max_mismatch_rate() -> f64 { 0.1 }
fn default_min_polyx_len() -> usize { 10 }
fn default_contaminant_threshold() -> f64 { 0.4 }

/// Profile loader and manager
pub struct Profile;

impl Profile {
    /// Load a profile by name or file path
    pub fn load(name_or_path: &str) -> Result<ProfileConfig> {
        // First try as a file path
        let path = PathBuf::from(name_or_path);
        if path.exists() {
            return Self::load_from_file(&path);
        }

        // Try built-in profiles
        if let Some(config) = Self::builtin(name_or_path) {
            return Ok(config);
        }

        // Try profiles directory
        let profiles_dir = Self::profiles_dir();
        let profile_path = profiles_dir.join(format!("{name_or_path}.yaml"));
        if profile_path.exists() {
            return Self::load_from_file(&profile_path);
        }

        anyhow::bail!(
            "Profile '{}' not found. Use 'virome-qc profiles' to list available profiles.",
            name_or_path
        );
    }

    fn load_from_file(path: &Path) -> Result<ProfileConfig> {
        let contents = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read profile: {}", path.display()))?;
        serde_yaml::from_str(&contents)
            .with_context(|| format!("Failed to parse profile: {}", path.display()))
    }

    fn profiles_dir() -> PathBuf {
        // Look relative to executable, then fallback to known locations
        let exe_dir = std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|p| p.to_path_buf()));

        if let Some(dir) = exe_dir {
            let profiles = dir.join("profiles");
            if profiles.exists() {
                return profiles;
            }
        }

        PathBuf::from("profiles")
    }

    /// List all available profile names
    pub fn list_available() -> Vec<String> {
        let mut profiles: Vec<String> = Self::builtin_names()
            .into_iter()
            .map(String::from)
            .collect();

        // Also scan profiles directory
        if let Ok(entries) = std::fs::read_dir(Self::profiles_dir()) {
            for entry in entries.flatten() {
                if let Some(name) = entry.path().file_stem() {
                    let name = name.to_string_lossy().to_string();
                    if !profiles.contains(&name) {
                        profiles.push(name);
                    }
                }
            }
        }

        profiles.sort();
        profiles
    }

    fn builtin_names() -> Vec<&'static str> {
        vec![
            "stool-vlp-tagmentation",
            "stool-vlp-truseq",
            "tissue-truseq",
            "metagenomics-nextera",
            "low-biomass-wga",
            "rna-virome-truseq",
        ]
    }

    /// Get a built-in profile by name
    fn builtin(name: &str) -> Option<ProfileConfig> {
        match name {
            "stool-vlp-tagmentation" => Some(Self::stool_vlp_tagmentation()),
            "tissue-truseq" => Some(Self::tissue_truseq()),
            "metagenomics-nextera" => Some(Self::metagenomics_nextera()),
            "low-biomass-wga" => Some(Self::low_biomass_wga()),
            _ => None,
        }
    }

    fn stool_vlp_tagmentation() -> ProfileConfig {
        ProfileConfig {
            name: "Stool VLP - Tagmentation".into(),
            description: "VLP-enriched stool virome with Nextera/tagmentation library prep".into(),
            platform: Platform::Illumina,
            modules: ModuleConfigs {
                adapter: AdapterConfig {
                    enabled: true,
                    sequences: vec!["nextera".into(), "truseq_universal".into()],
                    internal_scan: true,
                    random_primer_trim: 0,
                    min_overlap: 10,
                    max_mismatch_rate: 0.1,
                },
                quality: QualityConfig {
                    enabled: true,
                    window_size: 25,
                    min_quality: 15,
                    min_mean_quality: 20.0,
                    min_length: 90,
                    quality_binned: false,
                },
                polyx: PolyXConfig {
                    enabled: true,
                    platform_aware: true,
                    min_length: 10,
                },
                complexity: ComplexityConfig {
                    enabled: true,
                    min_entropy: 0.5,
                },
                contaminant: ContaminantConfig {
                    enabled: true,
                    screen_rrna: true,
                    screen_phix: true,
                    screen_vectors: true,
                    screen_kitome: false,
                    min_kmer_fraction: 0.4,
                },
                host: HostConfig {
                    enabled: true,
                    reference: "human".into(),
                    eve_aware: true,
                    rescue: true,
                },
                dedup: DedupConfig {
                    enabled: true,
                    optical_distance: 2500,
                    umi_aware: false,
                },
                chimera: ChimeraConfig {
                    enabled: false,
                },
            },
            thresholds: Thresholds {
                min_survival_rate: 0.30,
                max_host_fraction: 0.20,
                max_rrna_fraction: 0.05,
                max_duplicate_rate: 0.50,
            },
        }
    }

    fn tissue_truseq() -> ProfileConfig {
        ProfileConfig {
            name: "Tissue - TruSeq".into(),
            description: "Tissue/biopsy virome with TruSeq library prep, high host expected".into(),
            platform: Platform::Illumina,
            modules: ModuleConfigs {
                adapter: AdapterConfig {
                    enabled: true,
                    sequences: vec!["truseq".into(), "truseq_universal".into()],
                    internal_scan: true,
                    random_primer_trim: 0,
                    min_overlap: 10,
                    max_mismatch_rate: 0.1,
                },
                quality: QualityConfig {
                    enabled: true,
                    window_size: 25,
                    min_quality: 15,
                    min_mean_quality: 20.0,
                    min_length: 90,
                    quality_binned: false,
                },
                polyx: PolyXConfig {
                    enabled: true,
                    platform_aware: true,
                    min_length: 10,
                },
                complexity: ComplexityConfig {
                    enabled: true,
                    min_entropy: 0.6,
                },
                contaminant: ContaminantConfig {
                    enabled: true,
                    screen_rrna: true,
                    screen_phix: true,
                    screen_vectors: true,
                    screen_kitome: false,
                    min_kmer_fraction: 0.4,
                },
                host: HostConfig {
                    enabled: true,
                    reference: "human".into(),
                    eve_aware: true,
                    rescue: true,
                },
                dedup: DedupConfig {
                    enabled: true,
                    optical_distance: 2500,
                    umi_aware: false,
                },
                chimera: ChimeraConfig {
                    enabled: false,
                },
            },
            thresholds: Thresholds {
                min_survival_rate: 0.01, // tissue may lose >99% to host
                max_host_fraction: 0.99,
                max_rrna_fraction: 0.10,
                max_duplicate_rate: 0.50,
            },
        }
    }

    fn metagenomics_nextera() -> ProfileConfig {
        ProfileConfig {
            name: "Metagenomics - Nextera".into(),
            description: "Shotgun metagenomics with Nextera/tagmentation library prep".into(),
            platform: Platform::Illumina,
            modules: ModuleConfigs {
                adapter: AdapterConfig {
                    enabled: true,
                    sequences: vec!["nextera".into()],
                    internal_scan: true,
                    random_primer_trim: 0,
                    min_overlap: 10,
                    max_mismatch_rate: 0.1,
                },
                quality: QualityConfig {
                    enabled: true,
                    window_size: 25,
                    min_quality: 15,
                    min_mean_quality: 20.0,
                    min_length: 90,
                    quality_binned: false,
                },
                polyx: PolyXConfig {
                    enabled: true,
                    platform_aware: true,
                    min_length: 10,
                },
                complexity: ComplexityConfig {
                    enabled: true,
                    min_entropy: 0.6,
                },
                contaminant: ContaminantConfig {
                    enabled: true,
                    screen_rrna: true,
                    screen_phix: true,
                    screen_vectors: true,
                    screen_kitome: false,
                    min_kmer_fraction: 0.4,
                },
                host: HostConfig {
                    enabled: true,
                    reference: "human".into(),
                    eve_aware: true,
                    rescue: true,
                },
                dedup: DedupConfig {
                    enabled: true,
                    optical_distance: 2500,
                    umi_aware: false,
                },
                chimera: ChimeraConfig {
                    enabled: false,
                },
            },
            thresholds: Thresholds {
                min_survival_rate: 0.10,
                max_host_fraction: 0.80,
                max_rrna_fraction: 0.15,
                max_duplicate_rate: 0.50,
            },
        }
    }

    fn low_biomass_wga() -> ProfileConfig {
        ProfileConfig {
            name: "Low Biomass - WGA".into(),
            description: "Low biomass samples with whole-genome amplification (MDA/SISPA)".into(),
            platform: Platform::Illumina,
            modules: ModuleConfigs {
                adapter: AdapterConfig {
                    enabled: true,
                    sequences: vec!["nextera".into(), "truseq".into()],
                    internal_scan: true,
                    random_primer_trim: 13, // random hexamer bias extends ~13bp
                    min_overlap: 10,
                    max_mismatch_rate: 0.1,
                },
                quality: QualityConfig {
                    enabled: true,
                    window_size: 25,
                    min_quality: 15,
                    min_mean_quality: 20.0,
                    min_length: 90,
                    quality_binned: false,
                },
                polyx: PolyXConfig {
                    enabled: true,
                    platform_aware: true,
                    min_length: 10,
                },
                complexity: ComplexityConfig {
                    enabled: true,
                    min_entropy: 0.4, // WGA can produce lower-complexity fragments
                },
                contaminant: ContaminantConfig {
                    enabled: true,
                    screen_rrna: true,
                    screen_phix: true,
                    screen_vectors: true,
                    screen_kitome: true, // important for low-biomass
                    min_kmer_fraction: 0.4,
                },
                host: HostConfig {
                    enabled: true,
                    reference: "human".into(),
                    eve_aware: true,
                    rescue: true,
                },
                dedup: DedupConfig {
                    enabled: true,
                    optical_distance: 2500,
                    umi_aware: false,
                },
                chimera: ChimeraConfig {
                    enabled: true, // WGA produces chimeras
                },
            },
            thresholds: Thresholds {
                min_survival_rate: 0.10,
                max_host_fraction: 0.90,
                max_rrna_fraction: 0.10,
                max_duplicate_rate: 0.70, // high duplication expected with WGA
            },
        }
    }
}
