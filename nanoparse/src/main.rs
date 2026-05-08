//! nanoparse - Fast primer matching and BAM statistics
//!
//! A high-performance Rust CLI for nanopore amplicon analysis.

use anyhow::Result;
use clap::{Parser, Subcommand};
use nanoparse_core::{enrichment, matcher, pore_stats, stats, MatchMode};

#[derive(Parser)]
#[command(name = "nanoparse")]
#[command(author, version, about = "Fast primer matching and BAM statistics")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract read statistics from BAM
    Stats {
        /// Input BAM/FASTQ(.gz) file
        #[arg(short, long)]
        bam: String,

        /// Output file (JSON). Use - for stdout
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Number of threads
        #[arg(short, long, default_value = "8")]
        threads: usize,

        /// Minimum quality score filter
        #[arg(long, default_value = "0")]
        min_qs: f32,

        /// Minimum read length filter
        #[arg(long, default_value = "0")]
        min_len: usize,
    },

    /// Measure enrichment over BED regions and report pore statistics
    Enrichment {
        /// Input BAM file
        #[arg(short, long)]
        bam: String,

        /// BED file with target regions
        #[arg(short = 'd', long)]
        bed: String,

        /// Output file (JSON). Use - for stdout
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Number of threads for BAM decompression
        #[arg(short, long, default_value = "8")]
        threads: usize,

        /// Optional pore range from BAM cm:i tag, e.g. 1-2000
        #[arg(long)]
        cm_range: Option<String>,
    },

    /// Match reads to primers and identify amplicons
    Amplicons {
        /// Input BAM/FASTQ(.gz) file
        #[arg(short, long)]
        bam: String,

        /// Primers TSV file (name<TAB>sequence)
        #[arg(short, long)]
        primers: String,

        /// Output file (JSON). Use - for stdout
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Number of threads
        #[arg(short, long, default_value = "8")]
        threads: usize,

        /// Matching mode
        #[arg(short, long, value_enum, default_value = "semiglobal")]
        mode: MatchMode,

        /// Maximum edit distance for alignment
        #[arg(long, default_value = "3")]
        max_edit_dist: usize,

        /// Length of read ends to search for primers
        #[arg(long, default_value = "150")]
        end_length: usize,

        /// Tolerance in bp for fuzzy coordinate matching
        #[arg(long, default_value = "50")]
        primer_tolerance: i64,

        /// Minimum mean Q-score filter
        #[arg(long, default_value = "0")]
        min_qs: f32,

        /// Minimum read length filter
        #[arg(long, default_value = "0")]
        min_len: usize,

        /// Process at most this many reads (0 = no limit)
        #[arg(long, default_value = "0")]
        max_reads: usize,

        /// Keep only duplex reads (dx tag = 1)
        #[arg(long, default_value_t = false)]
        duplex_only: bool,

        /// Optional reference FASTA path (reserved for future use)
        #[arg(long)]
        reference: Option<String>,

        /// Optional GTF/GFF/BED path (reserved for future use)
        #[arg(long)]
        gtf: Option<String>,

        /// Print summary stats to stderr
        #[arg(long, default_value = "true")]
        summary: bool,
    },

    /// Calculate pore idle-time statistics
    PoreStats {
        /// Input BAM/FASTQ(.gz) file
        input: Option<String>,

        /// Optional sequencing summary TSV/TSV.GZ file
        #[arg(long)]
        sequencing_summary: Option<String>,

        /// Optional output file for full JSON results
        #[arg(short, long)]
        output: Option<String>,

        /// Number of threads for BAM decompression
        #[arg(short, long, default_value = "8")]
        threads: usize,

        /// Ignore idle times larger than this many seconds
        #[arg(long, default_value = "3600")]
        max_idle_s: f64,

        /// Threshold for counting an idle time as "long"
        #[arg(long, default_value = "60")]
        long_idle_s: f64,

        /// Sequencing speed in bases per second when duration must be estimated
        #[arg(long, default_value = "400")]
        speed_bps: f64,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Stats {
            bam,
            output,
            threads,
            min_qs,
            min_len,
        } => {
            stats::run_stats(&bam, &output, threads, min_qs, min_len)?;
        }
        Commands::Enrichment {
            bam,
            bed,
            output,
            threads,
            cm_range,
        } => {
            enrichment::run_enrichment(&bam, &bed, &output, threads, cm_range.as_deref())?;
        }
        Commands::Amplicons {
            bam,
            primers,
            output,
            threads,
            mode,
            max_edit_dist,
            end_length,
            primer_tolerance,
            min_qs,
            min_len,
            max_reads,
            duplex_only,
            reference,
            gtf,
            summary,
        } => {
            matcher::run_amplicons_to_output(
                &bam,
                &primers,
                &output,
                threads,
                mode,
                max_edit_dist,
                end_length,
                summary,
                primer_tolerance,
                min_qs,
                min_len,
                max_reads,
                duplex_only,
                reference.as_deref(),
                gtf.as_deref(),
            )?;
        }
        Commands::PoreStats {
            input,
            sequencing_summary,
            output,
            threads,
            max_idle_s,
            long_idle_s,
            speed_bps,
        } => {
            pore_stats::run_pore_stats(
                input.as_deref(),
                sequencing_summary.as_deref(),
                output.as_deref(),
                threads,
                max_idle_s,
                long_idle_s,
                speed_bps,
            )?;
        }
    }

    Ok(())
}
