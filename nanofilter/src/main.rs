use anyhow::Result;
use clap::{Parser, Subcommand};
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
use std::io::BufWriter;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "nanofilter",
    version,
    about = "Filtering and demultiplexing for FASTQ/BAM"
)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Filter {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value_t = 20.0)]
        qv: f64,
        #[arg(long, num_args = 1.., value_delimiter = ' ')]
        len: Vec<String>,
        #[arg(long, default_value = "auto")]
        output_format: String,
        #[arg(long)]
        channel_range: Option<String>,
        #[arg(long)]
        time_start: Option<String>,
        #[arg(long)]
        time_end: Option<String>,
    },
    Extract {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        channel_range: String,
        #[arg(long, default_value = "auto")]
        output_format: String,
    },
    Split {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        barcodes: PathBuf,
        #[arg(short, long, default_value = ".")]
        output_dir: PathBuf,
        #[arg(short, long, default_value_t = 1)]
        mismatches: usize,
        #[arg(short, long, default_value_t = 1000)]
        search_dist: usize,
        #[arg(short, long, default_value_t = 1)]
        threads: usize,
        #[arg(long, default_value_t = 0)]
        auto_discover: usize,
        #[arg(long)]
        fast: bool,
    },
    MakePe {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long, default_value_t = 150)]
        len: usize,
        #[arg(long, default_value_t = 400)]
        insert: usize,
        #[arg(long, default_value_t = 50)]
        step: usize,
    },
    Discover {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        barcodes: Option<PathBuf>,
        #[arg(short, long, default_value_t = 10000)]
        sample_size: usize,
        #[arg(short, long, default_value_t = 1)]
        threads: usize,
    },
    Umi {
        /// Input FASTQ or FASTQ.GZ file.
        #[arg(short, long)]
        input: PathBuf,
        /// Directory for output TSVs and per-family FASTQs.
        #[arg(short, long, default_value = ".")]
        output_dir: PathBuf,
        /// Forward upstream context sequence (immediately before the forward UMI).
        #[arg(long, default_value = "GTATCGTGTAGAGACTGCGTAGG")]
        fwd_context: String,
        /// Reverse upstream context sequence (immediately before the reverse UMI).
        #[arg(long, default_value = "AGTGATCGAGTCAGTGCGAGTG")]
        rev_context: String,
        /// IUPAC pattern for the forward UMI (e.g. TTTVVVVTTVVVVTTVVVVTTVVVVTTT).
        #[arg(long, default_value = "TTTVVVVTTVVVVTTVVVVTTVVVVTTT")]
        fwd_pattern: String,
        /// IUPAC pattern for the reverse UMI (e.g. AAABBBBAABBBBAABBBBAABBBBAAA).
        #[arg(long, default_value = "AAABBBBAABBBBAABBBBAABBBBAAA")]
        rev_pattern: String,
        /// Maximum edit distance for both anchor and UMI approximate matching.
        #[arg(long, default_value_t = 4)]
        max_edit: usize,
        /// Number of bases inspected at each read end.
        #[arg(long, default_value_t = 250)]
        window: usize,
        /// Minimum combined UMI length (fwd + rev).
        #[arg(long, default_value_t = 40)]
        min_umi_len: usize,
        /// Maximum combined UMI length (fwd + rev).
        #[arg(long, default_value_t = 75)]
        max_umi_len: usize,
        /// Emit wildcard-position-only normalised UMI (strand-normalised).
        #[arg(long)]
        normalize: bool,
        /// Minimum read length filter (0 = disabled).
        #[arg(long, default_value_t = 0)]
        min_read_len: usize,
        /// Maximum read length filter (0 = disabled).
        #[arg(long, default_value_t = 0)]
        max_read_len: usize,
        /// Minimum mean phred quality filter (0.0 = disabled).
        #[arg(long, default_value_t = 0.0)]
        min_mean_q: f64,
        /// Minimum reads required for a UMI family to pass.
        #[arg(long, default_value_t = 4)]
        min_reads: usize,
        /// Maximum reads per family used for consensus (oversized families are downsampled).
        #[arg(long, default_value_t = 80)]
        max_reads: usize,
        /// Require roughly equal forward/reverse reads per family.
        #[arg(long)]
        balance_strands: bool,
        /// Clustering mode: exact | approximate | vsearch.
        #[arg(long, default_value = "approximate")]
        cluster_mode: String,
        /// Max edit distance for approximate clustering.
        #[arg(long, default_value_t = 3)]
        cluster_edit: usize,
        /// vsearch pairwise identity threshold (only used with cluster-mode=vsearch).
        #[arg(long, default_value_t = 0.85)]
        vsearch_identity: f64,
        /// Consensus backend: none | medoid | majority_vote | medaka | dorado.
        #[arg(long, default_value = "none")]
        consensus_backend: String,
        /// Medaka model name (only used when consensus-backend=medaka).
        #[arg(long)]
        medaka_model: Option<String>,
        /// Alignment band half-width for majority_vote backend.
        #[arg(long, default_value_t = 150)]
        vote_band: usize,
        /// Sample name used in output filenames and summary row.
        #[arg(long, default_value = "sample")]
        sample: String,
        /// Amplicon name used in summary row.
        #[arg(long, default_value = "amplicon")]
        amplicon: String,
        /// Number of parallel worker threads (0 = use all available cores).
        #[arg(short, long, default_value_t = 4)]
        threads: usize,
        /// Expected inter-anchor span in bp (fwd_ctx start → rev_ctx start).
        /// Spans this value ± size_tolerance are accepted; 0 = disabled.
        #[arg(long, default_value_t = 0)]
        amplicon_size: usize,
        /// Tolerance around amplicon_size in bp (default 20% of amplicon_size when non-zero).
        #[arg(long, default_value_t = 0)]
        size_tolerance: usize,
        /// Path to append/create summary TSV (optional; header written on first creation).
        #[arg(long)]
        summary: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let args = Args::parse();

    match args.command {
        Commands::Filter {
            input,
            output,
            qv,
            len,
            output_format,
            channel_range,
            time_start,
            time_end,
        } => {
            let (min_len, max_len) = parse_len_args(&len)?;
            let channel_range = channel_range.as_deref().map(parse_pore_range).transpose()?;
            let time_window = match (time_start.as_deref(), time_end.as_deref()) {
                (Some(start), Some(end)) => Some(parse_time_window(start, end)?),
                (None, None) => None,
                _ => anyhow::bail!("Use both --time-start and --time-end together"),
            };

            let settings = bam::FilterSettings {
                qv_threshold: qv,
                min_len,
                max_len,
                channel_range,
                time_window,
            };

            let input_str = input.to_string_lossy().to_lowercase();
            let output_str = output.to_string_lossy().to_lowercase();
            let output_bam = match output_format.as_str() {
                "bam" => true,
                "fastq" | "fq" => false,
                "auto" => output_str.ends_with(".bam"),
                _ => anyhow::bail!(
                    "Unknown output format '{}'. Use 'auto', 'fastq', or 'bam'",
                    output_format
                ),
            };

            if input_str.ends_with(".bam") {
                if output_bam {
                    bam::filter_bam_with_settings(&input, &output, &settings)?;
                } else {
                    bam::bam_to_fastq_with_settings(&input, &output, &settings)?;
                }
            } else if is_fastq_path(&input_str) {
                fastq::filter_fastq_with_settings(&input, &output, &settings, output_bam)?;
            } else {
                anyhow::bail!("Unsupported input format");
            }
        }

        Commands::Extract {
            input,
            output,
            channel_range,
            output_format,
        } => {
            let output_str = output.to_string_lossy().to_lowercase();
            let output_bam = match output_format.as_str() {
                "bam" => true,
                "fastq" | "fq" => false,
                "auto" => output_str.ends_with(".bam"),
                _ => anyhow::bail!(
                    "Unknown output format '{}'. Use 'auto', 'fastq', or 'bam'",
                    output_format
                ),
            };

            let settings = bam::FilterSettings {
                qv_threshold: 0.0,
                min_len: 0,
                max_len: usize::MAX,
                channel_range: Some(parse_pore_range(&channel_range)?),
                time_window: None,
            };

            let input_str = input.to_string_lossy().to_lowercase();
            if input_str.ends_with(".bam") {
                if output_bam {
                    bam::filter_bam_with_settings(&input, &output, &settings)?;
                } else {
                    bam::bam_to_fastq_with_settings(&input, &output, &settings)?;
                }
            } else if is_fastq_path(&input_str) {
                fastq::filter_fastq_with_settings(&input, &output, &settings, output_bam)?;
            } else {
                anyhow::bail!("Unsupported input format");
            }
        }

        Commands::Split {
            input,
            barcodes,
            output_dir,
            mismatches,
            search_dist,
            threads,
            auto_discover,
            fast,
        } => {
            std::fs::create_dir_all(&output_dir)?;
            let input_str = input.to_string_lossy().to_lowercase();
            if input_str.ends_with(".bam") {
                bam::split_bam_by_barcodes(
                    &input,
                    &barcodes,
                    &output_dir,
                    mismatches,
                    search_dist,
                    threads,
                    auto_discover,
                    fast,
                )?;
            } else if is_fastq_path(&input_str) {
                fastq::split_fastq_by_barcodes(
                    &input,
                    &barcodes,
                    &output_dir,
                    mismatches,
                    search_dist,
                    threads,
                    auto_discover,
                    fast,
                )?;
            } else {
                anyhow::bail!("Unsupported input format");
            }
        }

        Commands::MakePe { input, output, len, insert, step } => {
            let input_str = input.to_string_lossy().to_lowercase();
            if is_fastq_path(&input_str) {
                let output_prefix =
                    output.unwrap_or_else(|| default_make_pe_prefix(&input, insert, step));
                fastq::make_pe(&input, &output_prefix, len, insert, step)?;
            } else {
                anyhow::bail!("make-PE currently supports FASTQ/FASTQ.GZ input only");
            }
        }

        Commands::Discover { input, barcodes, sample_size, threads } => {
            let input_str = input.to_string_lossy().to_lowercase();
            if input_str.ends_with(".bam") {
                bam::discover_barcodes(&input, barcodes.as_deref(), sample_size, threads)?;
            } else if is_fastq_path(&input_str) {
                fastq::discover_barcodes_fastq(
                    &input,
                    barcodes.as_deref(),
                    sample_size,
                    threads,
                )?;
            } else {
                anyhow::bail!("Unsupported input format");
            }
        }

        Commands::Umi {
            input,
            output_dir,
            fwd_context,
            rev_context,
            fwd_pattern,
            rev_pattern,
            max_edit,
            window,
            min_umi_len,
            max_umi_len,
            normalize,
            min_read_len,
            max_read_len,
            min_mean_q,
            min_reads,
            max_reads,
            balance_strands,
            cluster_mode,
            cluster_edit,
            vsearch_identity,
            consensus_backend,
            medaka_model,
            vote_band,
            sample,
            amplicon,
            threads,
            amplicon_size,
            size_tolerance,
            summary,
        } => {
            std::fs::create_dir_all(&output_dir)?;

            // Build UmiConfig — call .build() to compile kmer anchor indices.
            let umi_config = UmiConfig {
                fwd_context: fwd_context.into_bytes(),
                rev_context: rev_context.into_bytes(),
                fwd_pattern: IupacPattern::new(&fwd_pattern),
                rev_pattern: IupacPattern::new(&rev_pattern),
                max_edit_dist: max_edit,
                window_len: window,
                min_umi_len,
                max_umi_len,
                normalize,
                min_read_len,
                max_read_len: if max_read_len == 0 { usize::MAX } else { max_read_len },
                min_mean_q,
                amplicon_size,
                // If size_tolerance is 0 and amplicon_size is set, default to 20%.
                size_tolerance: if size_tolerance == 0 && amplicon_size > 0 {
                    amplicon_size / 5
                } else {
                    size_tolerance
                },
                fwd_index: AnchorIndex::empty(),
                rev_index: AnchorIndex::empty(),
            }
            .build();

            let input_label = input
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| input.to_string_lossy().to_string());

            let input_str = input.to_string_lossy().to_lowercase();
            if !is_fastq_path(&input_str) {
                anyhow::bail!("UMI detection currently supports FASTQ/FASTQ.GZ input only");
            }

            eprintln!("[nanofilter umi] Detecting UMIs in {} ({} threads)", input_label, threads);
            let all_records =
                run_umi_detection_fastq(&input, &input_label, &umi_config, threads)?;
            let umi_reads = all_records.iter().filter(|r| r.result.has_umi()).count();
            eprintln!("  {} / {} reads had detectable UMIs", umi_reads, all_records.len());

            // Per-read TSV
            let per_read_path = output_dir.join("detected_umis.tsv");
            write_per_read_tsv(&all_records, &per_read_path)?;
            eprintln!("  Per-read TSV -> {:?}", per_read_path);

            // Cluster
            let cluster_mode_parsed = match cluster_mode.as_str() {
                "exact" => ClusterMode::Exact,
                "vsearch" => ClusterMode::Vsearch { identity: vsearch_identity },
                _ => ClusterMode::Approximate { max_edit: cluster_edit },
            };
            // Consume all_families into filter — no clone needed (P5 fix).
            let all_families =
                nanofilter_core::cluster::cluster_by_umi(&all_records, &cluster_mode_parsed);
            eprintln!("  {} UMI families found", all_families.len());

            let cluster_config = ClusterFilterConfig { min_reads, max_reads, balance_strands };
            let filtered = nanofilter_core::cluster::filter_and_downsample(
                all_families, // consumed — no clone
                &cluster_config,
            );
            eprintln!(
                "  {} passing / {} failing families",
                filtered.passing.len(),
                filtered.failing.len()
            );

            // Write per-family FASTQs + optional consensus
            let backend =
                parse_consensus_backend(&consensus_backend, medaka_model.as_deref(), vote_band);
            let consensus_path =
                output_dir.join(format!("{}_{}_consensus.fastq", sample, amplicon));
            let mut consensus_writer = if matches!(backend, ConsensusBackend::None) {
                None
            } else {
                Some(BufWriter::new(File::create(&consensus_path)?))
            };
            let needs_family_fastq = matches!(backend, ConsensusBackend::Medaka { .. });
            for family in &filtered.passing {
                let fq_path = if needs_family_fastq {
                    Some(nanofilter_core::cluster::write_family_reads(family, &output_dir)?)
                } else {
                    None
                };
                let cons = nanofilter_core::consensus::derive_consensus(
                    family,
                    fq_path.as_deref(),
                    &backend,
                    &output_dir,
                );
                if let (Some(seq), Some(writer)) =
                    (cons.sequence.as_deref(), consensus_writer.as_mut())
                {
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
                use std::io::Write as _;
                writer.flush()?;
                eprintln!("  Consensus FASTQ -> {:?}", consensus_path);
            }

            // Cluster stats TSV — no all_families clone needed (P5 fix).
            let stats = nanofilter_core::cluster::compute_family_stats(
                &filtered.passing,
                &filtered.failing,
                &cluster_config,
            );
            let stats_path = output_dir.join(format!("{}_cluster_stats.tsv", sample));
            write_cluster_stats_tsv(&stats, &stats_path)?;
            eprintln!("  Cluster stats TSV -> {:?}", stats_path);

            // Optional summary row
            if let Some(summary_path) = summary {
                let summary_rec = build_summary(
                    &sample,
                    &amplicon,
                    &all_records,
                    &filtered.passing,
                    &filtered.failing,
                    &cluster_config,
                    &umi_config,
                );
                write_summary_tsv(&[summary_rec], &summary_path)?;
                eprintln!("  Summary appended -> {:?}", summary_path);
            }
        }
    }

    Ok(())
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
        _ => anyhow::bail!("Invalid --len value"),
    }
}

fn default_make_pe_prefix(input: &PathBuf, insert: usize, step: usize) -> PathBuf {
    let path = input.as_path();
    let parent = path.parent().unwrap_or_else(|| std::path::Path::new(""));
    let stem = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("reads");
    let base = stem
        .trim_end_matches(".fastq.gz")
        .trim_end_matches(".fq.gz")
        .trim_end_matches(".fastq")
        .trim_end_matches(".fq");
    parent.join(format!("{}_{}_{}", base, insert, step))
}
