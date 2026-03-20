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
            let profile_config = virome_qc::Profile::load(&profile)?;
            let pipeline = virome_qc::Pipeline::new(profile_config, threads);

            let result = match (r1, r2) {
                (Some(r1_path), Some(r2_path)) => {
                    // Paired-end mode
                    println!("Running paired-end QC...");
                    pipeline.run_paired(&r1_path, &r2_path, &output, merge)?
                }
                _ => {
                    // Single-end mode
                    let input_files = discover_fastq_files(&input)?;
                    pipeline.run(&input_files, &output)?
                }
            };

            // Write QA passport
            let passport = result.passport();
            passport.write_json(&output.join("passport.json"))?;

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
                let r1_path = output.with_extension("R1.fastq");
                let r2_path = output.with_extension("R2.fastq");
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
        Commands::Report { input, output } => {
            eprintln!("Report generation not yet implemented");
            let _ = (input, output);
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
