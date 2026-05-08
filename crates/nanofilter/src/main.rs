use anyhow::Result;
use clap::{Parser, Subcommand};
use nanostream_lib::{cli, run_discover, run_extract, run_filter, run_make_pe, run_split, run_umi};

#[derive(Parser)]
#[command(
    name = "nanofilter",
    version,
    about = "Filtering and demultiplexing for FASTQ/BAM (Legacy wrapper for nanostream)"
)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Filter reads by QV, length, channel, or time
    Filter(cli::FilterArgs),
    /// Extract reads from a channel range
    Extract(cli::ExtractArgs),
    /// Split reads by barcode pairs
    Split(cli::SplitArgs),
    /// Make pseudo paired-end FASTQ reads
    MakePe(cli::MakePeArgs),
    /// Discover barcode-like sequences
    Discover(cli::DiscoverArgs),
    /// Detect and cluster UMIs
    Umi(cli::UmiArgs),
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();

    match args.command {
        Commands::Filter(args) => run_filter(args)?,
        Commands::Extract(args) => run_extract(args)?,
        Commands::Split(args) => run_split(args)?,
        Commands::MakePe(args) => run_make_pe(args)?,
        Commands::Discover(args) => run_discover(args)?,
        Commands::Umi(args) => run_umi(args)?,
    }

    Ok(())
}
