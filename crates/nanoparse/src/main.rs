use anyhow::Result;
use clap::{Parser, Subcommand};
use nanoparse_core::{enrichment, matcher, pore_stats, stats, MatchMode};

#[derive(Parser)]
#[command(name = "nanoparse")]
#[command(author, version, about = "Fast primer matching and BAM statistics (Legacy wrapper for nanostream)")]
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
        #[arg(short, long, default_value = "-")]
        output: String,
        #[arg(short, long, default_value = "8")]
        threads: usize,
        #[arg(long, default_value = "0")]
        min_qs: f32,
        #[arg(long, default_value = "0")]
        min_len: usize,
    },
    /// Measure enrichment over BED regions
    Enrichment {
        #[arg(short, long)]
        bam: String,
        #[arg(short = 'd', long)]
        bed: String,
        #[arg(short, long, default_value = "-")]
        output: String,
        #[arg(short, long, default_value = "8")]
        threads: usize,
        #[arg(long)]
        cm_range: Option<String>,
    },
    /// Match reads to primers and identify amplicons
    Amplicons {
        #[arg(short, long)]
        bam: String,
        #[arg(short, long)]
        primers: String,
        #[arg(short, long, default_value = "-")]
        output: String,
        #[arg(short, long, default_value = "8")]
        threads: usize,
        #[arg(short, long, value_enum, default_value = "semiglobal")]
        mode: MatchMode,
        #[arg(long, default_value = "3")]
        max_edit_dist: usize,
        #[arg(long, default_value = "150")]
        end_length: usize,
        #[arg(long, default_value = "50")]
        primer_tolerance: i64,
        #[arg(long, default_value = "0")]
        min_qs: f32,
        #[arg(long, default_value = "0")]
        len: String,
        #[arg(long, default_value = "0")]
        max_reads: usize,
        #[arg(long, default_value_t = false)]
        duplex_only: bool,
        #[arg(long)]
        reference: Option<String>,
        #[arg(long)]
        gtf: Option<String>,
        #[arg(long, default_value = "true")]
        summary: bool,
        #[arg(long)]
        output_fastq: Option<String>,
        #[arg(long)]
        output_dimers: Option<String>,
        #[arg(long, default_value_t = false)]
        split_by_amplicon: bool,
        #[arg(long, default_value_t = false)]
        split_chimeras: bool,
    },

    /// Calculate pore idle-time statistics
    PoreStats {
        input: Option<String>,
        #[arg(long)]
        sequencing_summary: Option<String>,
        #[arg(short, long)]
        output: Option<String>,
        #[arg(short, long, default_value = "8")]
        threads: usize,
        #[arg(long, default_value = "3600")]
        max_idle_s: f64,
        #[arg(long, default_value = "60")]
        long_idle_s: f64,
        #[arg(long, default_value = "400")]
        speed_bps: f64,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Stats { bam, output, threads, min_qs, min_len } => {
            stats::run_stats(&bam, &output, threads, min_qs, min_len)?;
        }
        Commands::Enrichment { bam, bed, output, threads, cm_range } => {
            enrichment::run_enrichment(&bam, &bed, &output, threads, cm_range.as_deref())?;
        }
        Commands::Amplicons { bam, primers, output, threads, mode, max_edit_dist, end_length, primer_tolerance, min_qs, len, max_reads, duplex_only, reference, gtf, summary, output_fastq, output_dimers, split_by_amplicon, split_chimeras } => {
            let len_range = if len.is_empty() {
                (0, usize::MAX)
            } else {
                let parts: Vec<&str> = len
                    .split(|c| c == ',' || c == '-' || c == ' ')
                    .filter(|s| !s.is_empty())
                    .collect();
                match parts.as_slice() {
                    [min] => (min.parse::<usize>()?, usize::MAX),
                    [min, max] => (min.parse::<usize>()?, max.parse::<usize>()?),
                    _ => anyhow::bail!("Invalid --len value"),
                }
            };
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
                len_range,
                max_reads,
                duplex_only,
                reference.as_deref(),
                gtf.as_deref(),
                output_fastq.as_deref(),
                output_dimers.as_deref(),
                split_by_amplicon,
                split_chimeras,
            )?;
        }

        Commands::PoreStats { input, sequencing_summary, output, threads, max_idle_s, long_idle_s, speed_bps } => {
            pore_stats::run_pore_stats(input.as_deref(), sequencing_summary.as_deref(), output.as_deref(), threads, max_idle_s, long_idle_s, speed_bps)?;
        }
    }
    Ok(())
}
