use anyhow::{bail, Result};
use nanofilter_core::{
    bam,
    cluster::{ClusterFilterConfig, ClusterMode},
    consensus::{append_consensus_fastq_record, ConsensusBackend},
    fastq,
    report::{build_summary, write_cluster_stats_tsv, write_per_read_tsv, write_summary_tsv},
    umi::{run_umi_detection_fastq, AnchorIndex, IupacPattern, UmiConfig},
};
use nanoseq_core::filters::{parse_pore_range, parse_time_window};
use nanoseq_core::format::is_fastq_path;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub mod cli;

pub fn run_filter(args: cli::FilterArgs) -> Result<()> {
    let (min_len, max_len) = parse_len_args(&args.len)?;
    let channel_range = args
        .channel_range
        .as_deref()
        .map(parse_pore_range)
        .transpose()?;
    
    let time_window = match (args.time_start.as_deref(), args.time_end.as_deref()) {
        (Some(start), Some(end)) => {
            if chrono::DateTime::parse_from_rfc3339(start.trim()).is_ok()
                && chrono::DateTime::parse_from_rfc3339(end.trim()).is_ok()
            {
                Some(parse_time_window(start, end)?)
            } else {
                let (run_start, run_end) = if let Some(ref summary_path) = args.sequencing_summary {
                    eprintln!("Scanning sequencing summary file for run timing baseline...");
                    nanofilter_core::time_bounds::scan_sequencing_summary_bounds(summary_path)?
                } else {
                    eprintln!("Performing first-pass scan of input file for run timing baseline...");
                    let input_str = args.input.to_string_lossy().to_lowercase();
                    if input_str.ends_with(".bam") {
                        nanofilter_core::time_bounds::scan_bam_bounds(&args.input, args.threads)?
                    } else if is_fastq_path(&input_str) {
                        nanofilter_core::time_bounds::scan_fastq_bounds(&args.input)?
                    } else {
                        bail!("Unsupported input format for automatic time-filtering scan");
                    }
                };
                eprintln!("Discovered run baseline: start = {}, end = {}", run_start.to_rfc3339(), run_end.to_rfc3339());
                let resolved = nanoseq_core::filters::resolve_relative_window(start, end, run_start, run_end)?;
                eprintln!("Resolved relative time filter: start = {}, end = {}", resolved.start.to_rfc3339(), resolved.end.to_rfc3339());
                Some(resolved)
            }
        }
        (None, None) => None,
        _ => bail!("Use both --time-start and --time-end together"),
    };

    let settings = bam::FilterSettings {
        qv_threshold: args.qv,
        min_len,
        max_len,
        channel_range,
        time_window,
    };

    let output_bam = output_is_bam(&args.output, &args.output_format)?;
    run_filter_or_extract(
        &args.input,
        &args.output,
        &settings,
        output_bam,
        args.threads,
    )
}

pub fn run_extract(args: cli::ExtractArgs) -> Result<()> {
    let settings = bam::FilterSettings {
        qv_threshold: 0.0,
        min_len: 0,
        max_len: usize::MAX,
        channel_range: Some(parse_pore_range(&args.channel_range)?),
        time_window: None,
    };
    let output_bam = output_is_bam(&args.output, &args.output_format)?;
    run_filter_or_extract(
        &args.input,
        &args.output,
        &settings,
        output_bam,
        args.threads,
    )
}

fn run_filter_or_extract(
    input: &PathBuf,
    output: &PathBuf,
    settings: &bam::FilterSettings,
    output_bam: bool,
    threads: usize,
) -> Result<()> {
    let input_str = input.to_string_lossy().to_lowercase();
    if input_str.ends_with(".bam") {
        if output_bam {
            bam::filter_bam_with_settings_threads(input, output, settings, threads)?;
        } else {
            bam::bam_to_fastq_with_settings_threads(input, output, settings, threads)?;
        }
    } else if is_fastq_path(&input_str) {
        fastq::filter_fastq_with_settings(input, output, settings, output_bam)?;
    } else {
        bail!("Unsupported input format");
    }
    Ok(())
}

pub fn run_split(args: cli::SplitArgs) -> Result<()> {
    std::fs::create_dir_all(&args.output_dir)?;
    let input_str = args.input.to_string_lossy().to_lowercase();
    if input_str.ends_with(".bam") {
        bam::split_bam_by_barcodes(
            &args.input,
            &args.barcodes,
            &args.output_dir,
            args.mismatches,
            args.search_dist,
            args.threads,
            args.auto_discover,
            args.fast,
        )?;
    } else if is_fastq_path(&input_str) {
        fastq::split_fastq_by_barcodes(
            &args.input,
            &args.barcodes,
            &args.output_dir,
            args.mismatches,
            args.search_dist,
            args.threads,
            args.auto_discover,
            args.fast,
        )?;
    } else {
        bail!("Unsupported input format");
    }
    Ok(())
}

pub fn run_make_pe(args: cli::MakePeArgs) -> Result<()> {
    let input_str = args.input.to_string_lossy().to_lowercase();
    if !is_fastq_path(&input_str) {
        bail!("make-pe currently supports FASTQ/FASTQ.GZ input only");
    }
    let output_prefix = args
        .output
        .clone()
        .unwrap_or_else(|| default_make_pe_prefix(&args.input, args.insert, args.step));
    fastq::make_pe(
        &args.input,
        &output_prefix,
        args.len,
        args.insert,
        args.step,
    )
}

pub fn run_discover(args: cli::DiscoverArgs) -> Result<()> {
    let input_str = args.input.to_string_lossy().to_lowercase();
    if input_str.ends_with(".bam") {
        bam::discover_barcodes(
            &args.input,
            args.barcodes.as_deref(),
            args.sample_size,
            args.threads,
        )?;
    } else if is_fastq_path(&input_str) {
        fastq::discover_barcodes_fastq(
            &args.input,
            args.barcodes.as_deref(),
            args.sample_size,
            args.threads,
        )?;
    } else {
        bail!("Unsupported input format");
    }
    Ok(())
}

pub fn run_umi(args: cli::UmiArgs) -> Result<()> {
    std::fs::create_dir_all(&args.output_dir)?;
    let input_str = args.input.to_string_lossy().to_lowercase();
    if !is_fastq_path(&input_str) {
        bail!("UMI detection currently supports FASTQ/FASTQ.GZ input only");
    }

    let umi_config = UmiConfig {
        fwd_context: args.fwd_context.into_bytes(),
        rev_context: args.rev_context.into_bytes(),
        fwd_pattern: IupacPattern::new(&args.fwd_pattern),
        rev_pattern: IupacPattern::new(&args.rev_pattern),
        max_edit_dist: args.max_edit,
        window_len: args.window,
        min_umi_len: args.min_umi_len,
        max_umi_len: args.max_umi_len,
        normalize: args.normalize,
        min_read_len: args.min_read_len,
        max_read_len: if args.max_read_len == 0 {
            usize::MAX
        } else {
            args.max_read_len
        },
        min_mean_q: args.min_mean_q,
        amplicon_size: args.amplicon_size,
        size_tolerance: if args.size_tolerance == 0 && args.amplicon_size > 0 {
            args.amplicon_size / 5
        } else {
            args.size_tolerance
        },
        fwd_index: AnchorIndex::empty(),
        rev_index: AnchorIndex::empty(),
    }
    .build();

    let input_label = args
        .input
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| args.input.to_string_lossy().to_string());

    eprintln!(
        "[nanostream umi] Detecting UMIs in {} ({} threads)",
        input_label, args.threads
    );
    let all_records =
        run_umi_detection_fastq(&args.input, &input_label, &umi_config, args.threads)?;
    let umi_reads = all_records.iter().filter(|r| r.result.has_umi()).count();
    eprintln!(
        "  {} / {} reads had detectable UMIs",
        umi_reads,
        all_records.len()
    );

    let per_read_path = args.output_dir.join("detected_umis.tsv");
    write_per_read_tsv(&all_records, &per_read_path)?;
    eprintln!("  Per-read TSV -> {:?}", per_read_path);

    let cluster_mode = match args.cluster_mode.as_str() {
        "exact" => ClusterMode::Exact,
        "vsearch" => ClusterMode::Vsearch {
            identity: args.vsearch_identity,
        },
        _ => ClusterMode::Approximate {
            max_edit: args.cluster_edit,
        },
    };
    let all_families = nanofilter_core::cluster::cluster_by_umi(&all_records, &cluster_mode);
    eprintln!("  {} UMI families found", all_families.len());

    let cluster_config = ClusterFilterConfig {
        min_reads: args.min_reads,
        max_reads: args.max_reads,
        balance_strands: args.balance_strands,
    };
    let filtered = nanofilter_core::cluster::filter_and_downsample(all_families, &cluster_config);
    eprintln!(
        "  {} passing / {} failing families",
        filtered.passing.len(),
        filtered.failing.len()
    );

    let backend = parse_consensus_backend(
        &args.consensus_backend,
        args.medaka_model.as_deref(),
        args.vote_band,
    );
    let consensus_path = args
        .output_dir
        .join(format!("{}_{}_consensus.fastq", args.sample, args.amplicon));
    let mut consensus_writer = if matches!(backend, ConsensusBackend::None) {
        None
    } else {
        Some(BufWriter::new(File::create(&consensus_path)?))
    };
    let needs_family_fastq = matches!(backend, ConsensusBackend::Medaka { .. });
    for family in &filtered.passing {
        let fq_path = if needs_family_fastq {
            Some(nanofilter_core::cluster::write_family_reads(
                family,
                &args.output_dir,
            )?)
        } else {
            None
        };
        let cons = nanofilter_core::consensus::derive_consensus(
            family,
            fq_path.as_deref(),
            &backend,
            &args.output_dir,
        );
        if let (Some(seq), Some(writer)) = (cons.sequence.as_deref(), consensus_writer.as_mut()) {
            append_consensus_fastq_record(
                writer,
                cons.family_id,
                &cons.method,
                seq,
                cons.quality.as_deref(),
            )?;
        }
        if let Some(path) = fq_path {
            let _ = std::fs::remove_file(path);
        }
    }
    if let Some(writer) = consensus_writer.as_mut() {
        writer.flush()?;
        eprintln!("  Consensus FASTQ -> {:?}", consensus_path);
    }

    let stats = nanofilter_core::cluster::compute_family_stats(
        &filtered.passing,
        &filtered.failing,
        &cluster_config,
    );
    let stats_path = args
        .output_dir
        .join(format!("{}_cluster_stats.tsv", args.sample));
    write_cluster_stats_tsv(&stats, &stats_path)?;
    eprintln!("  Cluster stats TSV -> {:?}", stats_path);

    if let Some(summary_path) = args.summary {
        let summary = build_summary(
            &args.sample,
            &args.amplicon,
            &all_records,
            &filtered.passing,
            &filtered.failing,
            &cluster_config,
            &umi_config,
        );
        write_summary_tsv(&[summary], &summary_path)?;
        eprintln!("  Summary appended -> {:?}", summary_path);
    }

    Ok(())
}

fn output_is_bam(output: &PathBuf, output_format: &str) -> Result<bool> {
    let output_str = output.to_string_lossy().to_lowercase();
    match output_format {
        "bam" => Ok(true),
        "fastq" | "fq" => Ok(false),
        "auto" => Ok(output_str.ends_with(".bam")),
        _ => bail!(
            "Unknown output format '{}'. Use 'auto', 'fastq', or 'bam'",
            output_format
        ),
    }
}

fn parse_consensus_backend(
    name: &str,
    medaka_model: Option<&str>,
    vote_band: usize,
) -> ConsensusBackend {
    match name {
        "medoid" => ConsensusBackend::Medoid,
        "majority_vote" | "vote" => ConsensusBackend::MajorityVote { band: vote_band },
        "medaka" => ConsensusBackend::Medaka {
            model: medaka_model.map(|s| s.to_string()),
            min_length: 50,
            min_depth: 2,
        },
        "dorado" => ConsensusBackend::Dorado,
        _ => ConsensusBackend::None,
    }
}

fn parse_len_args(args: &[String]) -> Result<(usize, usize)> {
    if args.is_empty() {
        return Ok((0, usize::MAX));
    }

    let joined = args.join(" ");
    let parts: Vec<&str> = joined
        .split(|c| c == ',' || c == '-' || c == ' ')
        .filter(|s| !s.is_empty())
        .collect();

    match parts.as_slice() {
        [min] => Ok((min.parse()?, usize::MAX)),
        [min, max] => Ok((min.parse()?, max.parse()?)),
        _ => bail!("Invalid --len value"),
    }
}

fn default_make_pe_prefix(input: &PathBuf, insert: usize, step: usize) -> PathBuf {
    let path = input.as_path();
    let parent = path.parent().unwrap_or_else(|| std::path::Path::new(""));
    let stem = path.file_name().and_then(|n| n.to_str()).unwrap_or("reads");
    let base = stem
        .trim_end_matches(".fastq.gz")
        .trim_end_matches(".fq.gz")
        .trim_end_matches(".fastq")
        .trim_end_matches(".fq");
    parent.join(format!("{}_{}_{}", base, insert, step))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use std::io::Write;

    #[test]
    fn test_cli_relative_time_parsing() {
        let mut summary_file = NamedTempFile::new().unwrap();
        writeln!(
            summary_file,
            "read_id\tchannel\tstart_time\tduration\nread1\t101\t2026-03-17 12:00:00\t10.0\nread2\t102\t2026-03-17 12:30:00\t15.0"
        )
        .unwrap();

        // Construct mock FilterArgs
        let _args = cli::FilterArgs {
            input: PathBuf::from("mock.fastq"),
            output: PathBuf::from("mock_out.fastq"),
            qv: 7.0,
            len: vec![],
            output_format: "fastq".to_string(),
            channel_range: None,
            time_start: Some("10m".to_string()),
            time_end: Some("-5m".to_string()),
            threads: 4,
            sequencing_summary: Some(summary_file.path().to_path_buf()),
        };

        // Resolve time window using our summary scanning
        let (run_start, run_end) = nanofilter_core::time_bounds::scan_sequencing_summary_bounds(summary_file.path()).unwrap();
        let resolved = nanoseq_core::filters::resolve_relative_window("10m", "-5m", run_start, run_end).unwrap();

        let utc = chrono::FixedOffset::east_opt(0).unwrap();
        use chrono::TimeZone;
        assert_eq!(resolved.start, utc.with_ymd_and_hms(2026, 3, 17, 12, 10, 0).unwrap());
        assert_eq!(resolved.end, utc.with_ymd_and_hms(2026, 3, 17, 12, 25, 0).unwrap());
    }
}

