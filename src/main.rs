use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "virome-qc")]
#[command(about = "High-performance virome-specific QC platform")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run QC pipeline on FASTQ files
    Run {
        /// QC profile to use (e.g., stool-vlp-tagmentation, tissue-truseq)
        #[arg(short, long)]
        profile: String,

        /// Input directory or single-end FASTQ file
        #[arg(short, long)]
        input: PathBuf,

        /// Output directory for clean reads and QA passport
        #[arg(short, long)]
        output: PathBuf,

        /// Forward reads (paired-end mode)
        #[arg(short = '1', long)]
        r1: Option<PathBuf>,

        /// Reverse reads (paired-end mode)
        #[arg(short = '2', long)]
        r2: Option<PathBuf>,

        /// Merge overlapping paired-end reads
        #[arg(long, default_value = "false")]
        merge: bool,

        /// Number of threads
        #[arg(short, long, default_value = "4")]
        threads: usize,
    },

    /// List available QC profiles
    Profiles,

    /// Generate synthetic test corpus with ground truth labels
    Corpus {
        /// Output FASTQ file path (for single-end) or prefix (for paired-end)
        #[arg(short, long)]
        output: PathBuf,

        /// Sample type preset
        #[arg(short, long, default_value = "stool-vlp")]
        sample_type: String,

        /// Number of reads (pairs for paired-end) to generate
        #[arg(short, long, default_value = "100000")]
        num_reads: usize,

        /// Read length
        #[arg(short, long, default_value = "150")]
        read_length: usize,

        /// Generate paired-end reads (R1 + R2 files)
        #[arg(long, default_value = "false")]
        paired: bool,

        /// Random seed for reproducibility
        #[arg(long, default_value = "42")]
        seed: u64,
    },

    /// Scan FASTQ files to detect platform, characteristics, and potential issues
    Ingest {
        /// Forward reads (or single-end file)
        #[arg(short = '1', long)]
        r1: PathBuf,

        /// Reverse reads (paired-end)
        #[arg(short = '2', long)]
        r2: Option<PathBuf>,

        /// Output results as JSON instead of human-readable text
        #[arg(long)]
        json: bool,
    },

    /// Build database files (host filter, etc.)
    Db {
        /// Build a host Super Bloom filter from a FASTA reference
        #[arg(long)]
        host: Option<PathBuf>,

        /// Output path for the filter file
        #[arg(short, long, default_value = "host.sbf")]
        output: PathBuf,
    },

    /// Generate report from existing passport files
    Report {
        /// Input directory containing passport JSON files
        #[arg(short, long)]
        input: PathBuf,

        /// Output HTML report path
        #[arg(short, long)]
        output: PathBuf,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Run {
            profile,
            input,
            output,
            r1,
            r2,
            merge,
            threads,
        } => {
            let mut profile_config = virome_qc::Profile::load(&profile)?;

            // Run ingestion scan and apply platform-aware overrides
            let primary_input = r1.as_ref().unwrap_or(&input);
            if let Ok(ingest) = virome_qc::ingest_fastq(primary_input, r2.as_deref()) {
                if let Some(ref plat) = ingest.platform {
                    eprintln!(
                        "Platform: {} | {}bp {} | ~{} reads",
                        plat.model,
                        ingest.reads.read_length,
                        if ingest.reads.paired { "PE" } else { "SE" },
                        ingest
                            .reads
                            .estimated_read_count
                            .map(|n| format!("{:.1}M", n as f64 / 1e6))
                            .unwrap_or_else(|| "?".into())
                    );
                }
                if let Some(ref adapter) = ingest.quick_scan.dominant_adapter {
                    let rate = ingest
                        .quick_scan
                        .adapter_rates
                        .get(adapter)
                        .copied()
                        .unwrap_or(0.0);
                    eprintln!("Adapters:  {} detected ({:.1}%)", adapter, rate * 100.0);
                }
                for w in &ingest.warnings {
                    eprintln!("  Warning: {w}");
                }
                // Apply platform-detected overrides to profile
                profile_config = virome_qc::apply_overrides(profile_config, &ingest);
            }

            let pipeline = virome_qc::Pipeline::new(profile_config, threads);

            let result = match (r1, r2) {
                (Some(r1_path), Some(r2_path)) => {
                    pipeline.run_paired(&r1_path, &r2_path, &output, merge)?
                }
                _ => {
                    let input_files = discover_fastq_files(&input)?;
                    pipeline.run(&input_files, &output)?
                }
            };

            // Write QA passport
            let passport = result.passport();
            passport.write_json(&output.join("passport.json"))?;
            virome_qc::report::generate_html_report(&passport, &output.join("report.html"))?;

            println!(
                "QC complete: {}/{} reads passed ({:.1}%)",
                result.reads_passed,
                result.reads_input,
                result.survival_rate() * 100.0
            );
            if result.pairs_passed > 0 {
                println!(
                    "  Pairs passed: {}  Singletons: {}  Merged: {}",
                    result.pairs_passed, result.singletons, result.pairs_merged
                );
            }
        }
        Commands::Profiles => {
            println!("Available QC profiles:");
            for name in virome_qc::Profile::list_available() {
                println!("  {name}");
            }
        }
        Commands::Corpus {
            output,
            sample_type,
            num_reads,
            read_length,
            paired,
            seed,
        } => {
            let config = match sample_type.as_str() {
                "stool-vlp" => virome_qc::CorpusConfig {
                    num_reads,
                    read_length,
                    seed,
                    ..virome_qc::CorpusConfig::stool_vlp()
                },
                "tissue" => virome_qc::CorpusConfig {
                    num_reads,
                    read_length,
                    seed,
                    ..virome_qc::CorpusConfig::tissue()
                },
                "low-biomass-wga" => virome_qc::CorpusConfig {
                    num_reads,
                    read_length,
                    seed,
                    ..virome_qc::CorpusConfig::low_biomass_wga()
                },
                other => {
                    anyhow::bail!(
                        "Unknown sample type '{}'. Available: stool-vlp, tissue, low-biomass-wga",
                        other
                    );
                }
            };

            let mut generator = virome_qc::CorpusGenerator::new(config);
            let summary = if paired {
                let r1_path = output.with_extension("R1.fastq.gz");
                let r2_path = output.with_extension("R2.fastq.gz");
                println!(
                    "Generating {} paired-end reads ({sample_type} profile)...",
                    num_reads
                );
                let s = generator.generate_paired_to_files(&r1_path, &r2_path)?;
                println!("\nWritten to:");
                println!("  R1: {}", r1_path.display());
                println!("  R2: {}", r2_path.display());
                s
            } else {
                println!("Generating {} reads ({sample_type} profile)...", num_reads);
                let s = generator.generate_to_file(&output)?;
                println!("\nWritten to: {}", output.display());
                s
            };
            summary.print_report();
        }
        Commands::Ingest { r1, r2, json } => {
            let result = virome_qc::ingest_fastq(&r1, r2.as_deref())?;

            if json {
                println!("{}", serde_json::to_string_pretty(&result)?);
            } else {
                // Human-readable output
                if let Some(ref plat) = result.platform {
                    println!("Platform:     {} ({})", plat.model, plat.instrument_id);
                    println!("Chemistry:    {:?}", plat.chemistry);
                    println!(
                        "Flowcell:     {}",
                        plat.flowcell_id.as_deref().unwrap_or("--")
                    );
                    println!("Patterned:    {}", plat.patterned_flowcell);
                } else {
                    println!("Platform:     Unknown (non-Illumina or unrecognized header)");
                }
                println!();

                let r = &result.reads;
                println!(
                    "Read length:  {}bp{}",
                    r.read_length,
                    if r.variable_length { " (variable)" } else { "" }
                );
                println!(
                    "Layout:       {}",
                    if r.paired { "paired-end" } else { "single-end" }
                );
                println!("Quality enc:  Phred+{}", r.quality_offset);
                if let Some(est) = r.estimated_read_count {
                    println!("Est. reads:   ~{:.1}M", est as f64 / 1_000_000.0);
                }
                println!();

                let qs = &result.quick_scan;
                println!("--- Quick Scan ({} reads) ---", qs.reads_scanned);
                println!("Mean quality: {:.1}", qs.mean_quality);
                println!("Mean GC:      {:.1}%", qs.mean_gc * 100.0);
                println!("N-rate:       {:.3}%", qs.n_rate * 100.0);
                println!("N-rate pos 0: {:.2}%", qs.n_rate_pos0 * 100.0);
                println!(
                    "Q-binned:     {}{}",
                    qs.quality_binned,
                    if qs.quality_binned {
                        format!(" ({} distinct values)", qs.distinct_quality_values)
                    } else {
                        String::new()
                    }
                );
                println!(
                    "Complexity:   median={:.3} p2={:.3} p5={:.3}",
                    qs.complexity_median, qs.complexity_p2, qs.complexity_p5
                );
                if let Some(ref adapter) = qs.dominant_adapter {
                    let rate = qs.adapter_rates.get(adapter).copied().unwrap_or(0.0);
                    println!("Adapters:     {adapter} ({:.1}%)", rate * 100.0);
                } else {
                    println!("Adapters:     none detected in scan");
                }
                println!();

                if !result.recommendations.is_empty() {
                    println!("--- Recommendations ---");
                    for rec in &result.recommendations {
                        println!("  {} = {} ({})", rec.parameter, rec.value, rec.reason);
                    }
                    println!();
                }

                if !result.warnings.is_empty() {
                    println!("--- Warnings ---");
                    for w in &result.warnings {
                        println!("  {w}");
                    }
                }
            }
        }
        Commands::Db { host, output } => {
            if let Some(fasta_path) = host {
                println!(
                    "Building host Super Bloom filter from: {}",
                    fasta_path.display()
                );
                let start = std::time::Instant::now();
                virome_qc::modules::host::HostFilter::build_filter(&fasta_path, &output)?;
                println!(
                    "Done in {:.1}s. Filter: {}",
                    start.elapsed().as_secs_f64(),
                    output.display()
                );
            } else {
                eprintln!("Specify --host <reference.fa> to build a host filter");
            }
        }
        Commands::Report { input, output } => {
            let passport_path = if input.is_dir() {
                input.join("passport.json")
            } else {
                input.clone()
            };
            let contents = std::fs::read_to_string(&passport_path)?;
            let passport: virome_qc::Passport = serde_json::from_str(&contents)?;
            virome_qc::report::generate_html_report(&passport, &output)?;
            println!("Report written to: {}", output.display());
        }
    }

    Ok(())
}

fn discover_fastq_files(dir: &std::path::Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    if dir.is_file() {
        files.push(dir.to_path_buf());
    } else if dir.is_dir() {
        for entry in std::fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if let Some(ext) = path.extension() {
                let ext = ext.to_string_lossy().to_lowercase();
                if ext == "fastq" || ext == "fq" || ext == "gz" {
                    files.push(path);
                }
            }
        }
        files.sort();
    }
    Ok(files)
}
