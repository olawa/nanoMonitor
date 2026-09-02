use anyhow::Result;
use clap::Parser;
use nanostream_lib::{cli, run_discover, run_extract, run_filter, run_make_pe, run_split, run_umi};
use nanoparse_core::{enrichment, matcher, pore_stats, stats};

fn main() -> Result<()> {
    env_logger::init();
    let cli = cli::Cli::parse();

    match cli.command {
        cli::Commands::Stats {
            input,
            output,
            threads,
            min_qs,
            min_len,
        } => {
            stats::run_stats(&input, &output, threads, min_qs, min_len)?;
        }
        cli::Commands::Enrichment {
            input,
            bed,
            output,
            threads,
            cm_range,
        } => {
            enrichment::run_enrichment(&input, &bed, &output, threads, cm_range.as_deref())?;
        }
        cli::Commands::Amplicons {
            input,
            primers,
            output,
            threads,
            mode,
            max_edit_dist,
            end_length,
            primer_tolerance,
            min_qs,
            len,
            max_reads,
            duplex_only,
            reference,
            gtf,
            summary,
            output_fastq,
            output_dimers,
            split_by_amplicon,
            split_chimeras,
        } => {
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
                &input,
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

        cli::Commands::PoreStats {
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
        cli::Commands::Filter(args) => run_filter(args)?,
        cli::Commands::Extract(args) => run_extract(args)?,
        cli::Commands::Split(args) => run_split(args)?,
        cli::Commands::MakePe(args) => run_make_pe(args)?,
        cli::Commands::Discover(args) => run_discover(args)?,
        cli::Commands::Umi(args) => run_umi(args)?,
        cli::Commands::Monitor => {
            anyhow::bail!(
                "The GUI is currently built as `nanomonitor`. Use `cargo run -p nanomonitor` while monitor embedding is wired into `nanostream`."
            );
        }
    }

    Ok(())
}
