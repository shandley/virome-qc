//! Diagnostic: compute host containment distribution for input FASTQ
//!
//! Usage: host_containment_diag <fastq> <host.sbf>
//!
//! Outputs a histogram of k-mer containment fractions and lists reads
//! above the ambiguous threshold (0.20).

use std::path::PathBuf;
use superbloom::FrozenSuperBloom;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: host_containment_diag <fastq> <host.sbf>");
        std::process::exit(1);
    }

    let fastq_path = PathBuf::from(&args[1]);
    let filter_path = PathBuf::from(&args[2]);

    eprintln!("Loading host filter from {}...", filter_path.display());
    let filter = FrozenSuperBloom::load(&filter_path)
        .map_err(|e| anyhow::anyhow!("Failed to load filter: {}", e))?;
    eprintln!("Filter loaded.");

    // Histogram bins: 0-5%, 5-10%, 10-15%, ..., 95-100%
    let mut bins = [0u64; 20];
    let mut total_reads = 0u64;
    let mut above_ambiguous = Vec::new(); // (containment, name, seq_len)
    let mut above_host = Vec::new();

    let ambiguous_threshold = 0.20;
    let host_threshold = 0.50;

    let stream = biometal::FastqStream::from_path(&fastq_path)?;
    for record in stream.flatten() {
        total_reads += 1;
        let seq = &record.sequence;
        if seq.len() < 31 {
            continue;
        }

        let hits = filter.query_sequence(seq);
        if hits.is_empty() {
            bins[0] += 1;
            continue;
        }
        let n_hits = hits.iter().filter(|&&h| h).count();
        let containment = n_hits as f64 / hits.len() as f64;

        let bin = (containment * 20.0).min(19.0) as usize;
        bins[bin] += 1;

        if containment >= host_threshold {
            let name = record.id.clone();
            above_host.push((containment, name, seq.len()));
        } else if containment >= ambiguous_threshold {
            let name = record.id.clone();
            above_ambiguous.push((containment, name, seq.len()));
        }
    }

    println!("=== Host Containment Distribution ===");
    println!("Total reads: {}", total_reads);
    println!();
    println!("Containment Range  |  Count  |  Fraction");
    println!("-------------------|---------|----------");
    for (i, &count) in bins.iter().enumerate() {
        let lo = i as f64 * 5.0;
        let hi = lo + 5.0;
        let frac = count as f64 / total_reads as f64;
        let bar = "#".repeat((frac * 200.0) as usize);
        if count > 0 {
            println!(
                "{:5.0}% - {:5.0}%  | {:>7} | {:>8.4}%  {}",
                lo,
                hi,
                count,
                frac * 100.0,
                bar
            );
        }
    }

    println!();
    println!(
        "=== Reads above host threshold ({:.0}%) ===  Count: {}",
        host_threshold * 100.0,
        above_host.len()
    );
    for (c, name, len) in above_host.iter().take(20) {
        println!("  {:.1}% containment, {}bp: {}", c * 100.0, len, name);
    }
    if above_host.len() > 20 {
        println!("  ... and {} more", above_host.len() - 20);
    }

    println!();
    println!(
        "=== Reads in ambiguous zone ({:.0}%-{:.0}%) ===  Count: {}",
        ambiguous_threshold * 100.0,
        host_threshold * 100.0,
        above_ambiguous.len()
    );
    for (c, name, len) in above_ambiguous.iter().take(20) {
        println!("  {:.1}% containment, {}bp: {}", c * 100.0, len, name);
    }
    if above_ambiguous.len() > 20 {
        println!("  ... and {} more", above_ambiguous.len() - 20);
    }

    // Count by source category
    fn count_by_source(reads: &[(f64, String, usize)]) -> std::collections::HashMap<String, usize> {
        let mut counts = std::collections::HashMap::new();
        for (_, name, _) in reads {
            let source = name
                .split("source=")
                .nth(1)
                .and_then(|s| s.split(';').next())
                .unwrap_or("unknown")
                .to_string();
            *counts.entry(source).or_insert(0) += 1;
        }
        counts
    }

    println!();
    println!("=== Host zone (>{:.0}%) by source ===", host_threshold * 100.0);
    let host_by_source = count_by_source(&above_host);
    let mut host_sorted: Vec<_> = host_by_source.iter().collect();
    host_sorted.sort_by(|a, b| b.1.cmp(a.1));
    for (src, count) in &host_sorted {
        println!("  {:20}: {:>6}", src, count);
    }

    println!();
    println!("=== Ambiguous zone ({:.0}%-{:.0}%) by source ===", ambiguous_threshold * 100.0, host_threshold * 100.0);
    let ambig_by_source = count_by_source(&above_ambiguous);
    let mut ambig_sorted: Vec<_> = ambig_by_source.iter().collect();
    ambig_sorted.sort_by(|a, b| b.1.cmp(a.1));
    for (src, count) in &ambig_sorted {
        println!("  {:20}: {:>6}", src, count);
    }

    // True viral FP
    let viral_fp_host = host_by_source.get("viral").unwrap_or(&0)
        + host_by_source.get("background").unwrap_or(&0)
        + host_by_source.get("phage").unwrap_or(&0);
    let viral_fp_ambig = ambig_by_source.get("viral").unwrap_or(&0)
        + ambig_by_source.get("background").unwrap_or(&0)
        + ambig_by_source.get("phage").unwrap_or(&0);

    println!();
    println!("=== Summary ===");
    println!("Total reads: {}", total_reads);
    println!(
        "Host zone (>{:.0}%): {} reads ({:.4}%)",
        host_threshold * 100.0,
        above_host.len(),
        above_host.len() as f64 / total_reads as f64 * 100.0
    );
    println!(
        "  True viral FP in host zone: {} ({:.5}% of all reads)",
        viral_fp_host,
        viral_fp_host as f64 / total_reads as f64 * 100.0
    );
    println!(
        "Ambiguous zone ({:.0}%-{:.0}%): {} reads ({:.4}%)",
        ambiguous_threshold * 100.0,
        host_threshold * 100.0,
        above_ambiguous.len(),
        above_ambiguous.len() as f64 / total_reads as f64 * 100.0
    );
    println!(
        "  True viral FP in ambiguous zone: {} ({:.5}% of all reads)",
        viral_fp_ambig,
        viral_fp_ambig as f64 / total_reads as f64 * 100.0
    );

    Ok(())
}
