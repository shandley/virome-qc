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

        /// Input directory or single-end FASTQ file (optional when -1/-2 are provided)
        #[arg(short, long)]
        input: Option<PathBuf>,

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

        /// Report only: generate passport and HTML report without writing clean FASTQ
        #[arg(long, default_value = "false")]
        report_only: bool,

        /// Overwrite existing results in output directory
        #[arg(long, default_value = "false")]
        force: bool,

        /// Number of threads
        #[arg(short, long, default_value = "4")]
        threads: usize,
    },

    /// List available QC profiles
    Profiles {
        /// Show full parameters for a specific profile
        #[arg(long)]
        show: Option<String>,
    },

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
        /// Input FASTQ file (single-end or R1)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Forward reads (or single-end file). Alias for --input.
        #[arg(short = '1', long)]
        r1: Option<PathBuf>,

        /// Reverse reads (paired-end)
        #[arg(short = '2', long)]
        r2: Option<PathBuf>,

        /// Output results as JSON instead of human-readable text
        #[arg(long)]
        json: bool,
    },

    /// Deduplicate clean FASTQ files (post-pipeline step)
    Dedup {
        /// Input FASTQ file (single-end) or R1 (paired-end)
        #[arg(short, long)]
        input: PathBuf,

        /// R2 input for paired-end dedup
        #[arg(short = '2', long)]
        r2: Option<PathBuf>,

        /// Output directory
        #[arg(short, long)]
        output: PathBuf,

        /// Optical duplicate pixel distance (100 non-patterned, 2500 patterned)
        #[arg(long, default_value = "2500")]
        optical_distance: u32,

        /// Enable UMI-aware dedup (extracts UMI from read names)
        #[arg(long)]
        umi: bool,
    },

    /// Generate multi-sample batch comparison report
    Batch {
        /// Directory containing sample result directories (each with passport.json)
        #[arg(short, long)]
        input: PathBuf,

        /// Output HTML report path
        #[arg(short, long, default_value = "batch_report.html")]
        output: PathBuf,
    },

    /// Build database files (host filter, etc.)
    Db {
        /// Build a host Super Bloom filter from a FASTA reference
        #[arg(long)]
        host: Option<PathBuf>,

        /// Build rRNA k-mer filter from SILVA FASTA file(s)
        /// Provide one or more gzipped FASTA files (SSU + LSU)
        #[arg(long, num_args = 0..)]
        rrna: Option<Vec<PathBuf>>,

        /// Output path for the filter file
        #[arg(short, long, default_value = "host.sbf")]
        output: PathBuf,

        /// Download and build all required databases (T2T human host + SILVA rRNA).
        /// Files are stored in ~/.virome-qc/ for auto-discovery.
        #[arg(long)]
        setup: bool,
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
            report_only,
            force,
            threads,
        } => {
            // Input validation: require either -i or -1/-2
            let primary_input = match (&r1, &input) {
                (Some(r1_path), _) => r1_path.clone(),
                (None, Some(input_path)) => input_path.clone(),
                (None, None) => {
                    anyhow::bail!(
                        "No input specified. Provide -i <file_or_dir> or -1 <R1> -2 <R2>."
                    );
                }
            };

            // Validate input files exist
            if !primary_input.exists() {
                anyhow::bail!(
                    "Input file not found: {}",
                    primary_input.display()
                );
            }
            if let Some(ref r2_path) = r2 {
                if !r2_path.exists() {
                    anyhow::bail!(
                        "R2 file not found: {}",
                        r2_path.display()
                    );
                }
            }

            // Output directory overwrite protection
            if output.join("passport.json").exists() && !force {
                anyhow::bail!(
                    "Output directory '{}' already contains results (passport.json found). \
                     Use --force to overwrite.",
                    output.display()
                );
            }

            let mut profile_config = virome_qc::Profile::load(&profile)?;

            // Run ingestion scan and apply platform-aware overrides
            let mut ingestion_result = None;
            if let Ok(ingest) = virome_qc::ingest_fastq(&primary_input, r2.as_deref()) {
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
                // Early GC validation against profile expected range
                if let Some((gc_min, gc_max)) = profile_config.thresholds.expected_gc_range {
                    let gc = ingest.quick_scan.mean_gc;
                    if gc < gc_min || gc > gc_max {
                        eprintln!(
                            "  Warning: Mean GC {:.1}% is outside expected range ({:.0}-{:.0}%) for this profile",
                            gc * 100.0,
                            gc_min * 100.0,
                            gc_max * 100.0
                        );
                    }
                }
                for w in &ingest.warnings {
                    eprintln!("  Warning: {w}");
                }
                // Apply platform-detected overrides to profile
                profile_config = virome_qc::apply_overrides(profile_config, &ingest);
                ingestion_result = Some(ingest);
            }

            let final_config = profile_config.clone();
            let pipeline = virome_qc::Pipeline::new(profile_config, threads);

            let mut result = match (r1, r2) {
                (Some(r1_path), Some(r2_path)) => {
                    pipeline.run_paired(&r1_path, &r2_path, &output, merge)?
                }
                _ => {
                    let input_files = discover_fastq_files(&primary_input)?;
                    pipeline.run(&input_files, &output)?
                }
            };

            // Attach ingestion results and applied config to pipeline result
            result.ingestion = ingestion_result;
            result.applied_config = Some(final_config);

            // Write QA passport and HTML report
            let passport = result.passport();

            passport.write_json(&output.join("passport.json"))?;
            virome_qc::report::generate_html_report(&passport, &output.join("report.html"))?;

            // In report-only mode, remove FASTQ output files (keep passport + report)
            if report_only {
                for entry in std::fs::read_dir(&output)? {
                    let entry = entry?;
                    let path = entry.path();
                    if let Some(ext) = path.extension() {
                        if ext == "gz" || ext == "fastq" || ext == "fq" {
                            std::fs::remove_file(&path)?;
                        }
                    }
                }
            }

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
        Commands::Profiles { show } => {
            if let Some(name) = show {
                let config = virome_qc::Profile::load(&name)?;
                println!("{}", serde_yaml::to_string(&config)?);
            } else {
                println!("Available QC profiles:");
                for (name, desc) in virome_qc::Profile::list_with_descriptions() {
                    println!("  {name:30} {desc}");
                }
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
        Commands::Ingest { input, r1, r2, json } => {
            let r1_path = r1.or(input).ok_or_else(|| {
                anyhow::anyhow!("No input specified. Provide -i <file> or -1 <R1>.")
            })?;
            if !r1_path.exists() {
                anyhow::bail!("Input file not found: {}", r1_path.display());
            }
            let result = virome_qc::ingest_fastq(&r1_path, r2.as_deref())?;

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
        Commands::Batch { input, output } => {
            virome_qc::report::batch::generate_batch_report(&input, &output)?;
            println!("Batch report written to: {}", output.display());
        }
        Commands::Dedup {
            input,
            r2,
            output,
            optical_distance,
            umi,
        } => {
            std::fs::create_dir_all(&output)?;
            let config = virome_qc::modules::dedup::DedupConfig {
                optical_distance,
                umi_aware: umi,
                prefix_len: 50,
            };

            let stats = if let Some(r2_path) = r2 {
                // Paired-end dedup
                let r1_out = output.join("dedup_R1.fastq.gz");
                let r2_out = output.join("dedup_R2.fastq.gz");
                println!("Deduplicating paired-end reads...");
                virome_qc::modules::dedup::dedup_paired_end(
                    &input, &r2_path, &r1_out, &r2_out, &config,
                )?
            } else {
                // Single-end dedup
                let se_out = output.join("dedup.fastq.gz");
                println!("Deduplicating single-end reads...");
                virome_qc::modules::dedup::dedup_single_end(&input, &se_out, &config)?
            };

            println!("\n=== Dedup Results ===");
            println!("Total reads:   {:>12}", stats.total_reads);
            println!("Unique reads:  {:>12}", stats.unique_reads);
            println!("PCR duplicates:{:>12}", stats.pcr_duplicates);
            println!("Duplicate rate:{:>11.1}%", stats.duplicate_rate * 100.0);
            println!(
                "Library complexity: ~{}",
                stats.estimated_library_complexity
            );

            // Save stats as JSON
            let stats_path = output.join("dedup_stats.json");
            let stats_json = serde_json::to_string_pretty(&stats)?;
            std::fs::write(&stats_path, stats_json)?;
            println!("\nStats written to: {}", stats_path.display());
        }
        Commands::Db {
            host,
            rrna,
            output,
            setup,
        } => {
            if setup {
                // Auto-setup: download T2T human genome + SILVA, build both filters
                let db_dir = std::env::var("HOME")
                    .map(PathBuf::from)
                    .unwrap_or_else(|_| PathBuf::from("."))
                    .join(".virome-qc");
                std::fs::create_dir_all(&db_dir)?;

                println!("=== virome-qc database setup ===");
                println!("Installing to: {}", db_dir.display());
                println!();

                // Build rRNA filter
                let rrna_path = db_dir.join("rrna").join("silva.rrf");
                if rrna_path.exists() {
                    println!("rRNA filter already exists: {}", rrna_path.display());
                } else {
                    println!("--- Building SILVA rRNA filter ---");
                    std::fs::create_dir_all(rrna_path.parent().unwrap())?;
                    // Trigger the same auto-download logic as --rrna with no args
                    let rrna_dir = rrna_path.parent().unwrap();
                    let ssu_url = "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz";
                    let lsu_url = "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva_trunc.fasta.gz";
                    let ssu_file = rrna_dir.join("SILVA_SSU.fasta.gz");
                    let lsu_file = rrna_dir.join("SILVA_LSU.fasta.gz");

                    for (url, dest) in [(ssu_url, &ssu_file), (lsu_url, &lsu_file)] {
                        if !dest.exists() {
                            println!("  Downloading {}...", dest.file_name().unwrap().to_string_lossy());
                            let status = std::process::Command::new("curl")
                                .args(["-L", "-o", &dest.to_string_lossy(), url])
                                .status()?;
                            if !status.success() {
                                anyhow::bail!("Failed to download SILVA");
                            }
                        }
                    }

                    virome_qc::modules::rrna::RrnaFilter::build_filter(
                        &[ssu_file.as_path(), lsu_file.as_path()],
                        &rrna_path,
                    )?;
                    println!("  rRNA filter: {}", rrna_path.display());
                }

                // Build host filter
                let host_path = db_dir.join("host").join("human_t2t.sbf");
                if host_path.exists() {
                    println!("Host filter already exists: {}", host_path.display());
                } else {
                    println!();
                    println!("--- Building T2T-CHM13 host filter ---");
                    std::fs::create_dir_all(host_path.parent().unwrap())?;
                    let host_dir = host_path.parent().unwrap();
                    let t2t_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz";
                    let t2t_gz = host_dir.join("t2t_chm13.fna.gz");
                    let t2t_fna = host_dir.join("t2t_chm13.fna");

                    if !t2t_fna.exists() {
                        println!("  Downloading T2T-CHM13 (~890 MB)...");
                        let status = std::process::Command::new("curl")
                            .args(["-L", "-o", &t2t_gz.to_string_lossy(), t2t_url])
                            .status()?;
                        if !status.success() {
                            anyhow::bail!("Failed to download T2T-CHM13");
                        }
                        println!("  Decompressing...");
                        let status = std::process::Command::new("gunzip")
                            .arg(&t2t_gz)
                            .status()?;
                        if !status.success() {
                            anyhow::bail!("Failed to decompress T2T-CHM13");
                        }
                    }

                    virome_qc::modules::host::HostFilter::build_filter(&t2t_fna, &host_path)?;
                    // Clean up uncompressed FASTA
                    let _ = std::fs::remove_file(&t2t_fna);
                    println!("  Host filter: {}", host_path.display());
                }

                println!();
                println!("=== Setup complete ===");
                println!("Filters installed to: {}", db_dir.display());
                println!("virome-qc will auto-discover these filters on future runs.");
                return Ok(());
            }
            if let Some(ref rrna_paths) = rrna {
                let rrna_output = if output.to_string_lossy() == "host.sbf" {
                    std::path::PathBuf::from("rrna.rrf")
                } else {
                    output.clone()
                };

                let downloaded_paths: Vec<PathBuf>;
                let paths: Vec<&std::path::Path> = if rrna_paths.is_empty() {
                    // Auto-download SILVA SSU + LSU NR99
                    let silva_dir = rrna_output.parent().unwrap_or(std::path::Path::new("."));
                    let ssu_url = "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz";
                    let lsu_url = "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva_trunc.fasta.gz";

                    let ssu_path = silva_dir.join("SILVA_138.2_SSU_NR99.fasta.gz");
                    let lsu_path = silva_dir.join("SILVA_138.2_LSU_NR99.fasta.gz");

                    for (url, path, name) in [
                        (ssu_url, &ssu_path, "SSU (16S/18S)"),
                        (lsu_url, &lsu_path, "LSU (23S/28S)"),
                    ] {
                        if !path.exists() {
                            eprintln!("Downloading SILVA 138.2 {}...", name);
                            let status = std::process::Command::new("curl")
                                .args(["-L", "-o", &path.to_string_lossy(), url])
                                .status()?;
                            if !status.success() {
                                anyhow::bail!("Failed to download SILVA {}", name);
                            }
                        } else {
                            eprintln!("Using cached {}", path.display());
                        }
                    }

                    downloaded_paths = vec![ssu_path, lsu_path];
                    downloaded_paths.iter().map(|p| p.as_path()).collect()
                } else {
                    rrna_paths.iter().map(|p| p.as_path()).collect()
                };

                let start = std::time::Instant::now();
                virome_qc::modules::rrna::RrnaFilter::build_filter(&paths, &rrna_output)?;
                println!(
                    "Done in {:.1}s. Filter: {}",
                    start.elapsed().as_secs_f64(),
                    rrna_output.display()
                );
            } else if let Some(ref host_arg) = host {
                let fasta_path = if host_arg.exists() {
                    // Direct file path
                    host_arg.clone()
                } else {
                    // Name-based: auto-download known references
                    let name = host_arg.to_string_lossy().to_lowercase();
                    let url = match name.as_str() {
                        "human" | "t2t" | "chm13" => {
                            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
                        }
                        other => {
                            anyhow::bail!(
                                "Unknown reference '{}'. Provide a FASTA file path, or use 'human' for T2T-CHM13.",
                                other
                            );
                        }
                    };
                    let download_path = output.with_extension("fa.gz");
                    println!("Downloading {} reference...", name);
                    let status = std::process::Command::new("curl")
                        .args(["-L", "-o", &download_path.to_string_lossy(), url])
                        .status()?;
                    if !status.success() {
                        anyhow::bail!("Download failed");
                    }
                    println!("Decompressing...");
                    let status = std::process::Command::new("gunzip")
                        .args(["-f", &download_path.to_string_lossy()])
                        .status()?;
                    if !status.success() {
                        anyhow::bail!("Decompression failed");
                    }
                    download_path.with_extension("") // remove .gz
                };

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
                eprintln!("Specify --host <reference.fa> or --rrna <silva.fasta.gz> to build a filter");
            }
        }
        Commands::Report { input, output } => {
            let passport_path = if input.is_dir() {
                input.join("passport.json")
            } else {
                input.clone()
            };
            let contents = std::fs::read_to_string(&passport_path)?;
            let mut passport: virome_qc::Passport = serde_json::from_str(&contents)?;

            // Recompute flags using latest logic (handles passports from older code)
            passport.recompute_flags();

            virome_qc::report::generate_html_report(&passport, &output)?;
            println!("Report written to: {}", output.display());
        }
    }

    Ok(())
}

fn discover_fastq_files(dir: &std::path::Path) -> Result<Vec<PathBuf>> {
    if !dir.exists() {
        anyhow::bail!("Input path does not exist: {}", dir.display());
    }
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
    if files.is_empty() {
        anyhow::bail!(
            "No FASTQ files found in: {}",
            dir.display()
        );
    }
    Ok(files)
}
