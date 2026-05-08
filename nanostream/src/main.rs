//! Unified CLI for nanoStream Rust tools.

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
use nanoparse_core::{enrichment, matcher, pore_stats, stats, MatchMode};
use nanoseq_core::filters::{parse_pore_range, parse_time_window};
use nanoseq_core::format::is_fastq_path;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "nanostream")]
#[command(
    author,
    version,
    about = "Unified Rust CLI for nanopore BAM/FASTQ workflows"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract read statistics from BAM/FASTQ
    Stats {
        /// Input BAM/FASTQ(.gz) file
        input: String,
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
    /// Measure enrichment over BED regions
    Enrichment {
        /// Input BAM file
        input: String,
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
        input: String,
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
        /// Optional reference FASTA path
        #[arg(long)]
        reference: Option<String>,
        /// Optional GTF/GFF/BED path
        #[arg(long)]
        gtf: Option<String>,
        /// Print summary stats to stderr
        #[arg(long, default_value = "true")]
        summary: bool,
    },
    /// Calculate pore idle-time statistics
    PoreStats {
        /// Optional input BAM/FASTQ(.gz) file
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
        /// Threshold for counting an idle time as long
        #[arg(long, default_value = "60")]
        long_idle_s: f64,
        /// Sequencing speed in bases per second when duration must be estimated
        #[arg(long, default_value = "400")]
        speed_bps: f64,
    },
    /// Filter reads by QV, length, channel, or time
    Filter(FilterArgs),
    /// Extract reads from a channel range
    Extract(ExtractArgs),
    /// Split reads by barcode pairs
    Split(SplitArgs),
    /// Make pseudo paired-end FASTQ reads
    MakePe(MakePeArgs),
    /// Discover barcode-like sequences
    Discover(DiscoverArgs),
    /// Detect and cluster UMIs
    Umi(UmiArgs),
    /// Launch the native GUI monitor
    Monitor,
}

#[derive(Parser)]
struct FilterArgs {
    /// Input BAM/FASTQ(.gz) file
    input: PathBuf,
    /// Output BAM/FASTQ(.gz) file
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
    /// Number of threads for threaded BAM IO and parallel FASTQ operations where supported
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}

#[derive(Parser)]
struct ExtractArgs {
    input: PathBuf,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long)]
    channel_range: String,
    #[arg(long, default_value = "auto")]
    output_format: String,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}

#[derive(Parser)]
struct SplitArgs {
    input: PathBuf,
    #[arg(short, long)]
    barcodes: PathBuf,
    #[arg(short, long, default_value = ".")]
    output_dir: PathBuf,
    #[arg(short, long, default_value_t = 1)]
    mismatches: usize,
    #[arg(short, long, default_value_t = 1000)]
    search_dist: usize,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value_t = 0)]
    auto_discover: usize,
    #[arg(long)]
    fast: bool,
}

#[derive(Parser)]
struct MakePeArgs {
    input: PathBuf,
    #[arg(short, long)]
    output: Option<PathBuf>,
    #[arg(long, default_value_t = 150)]
    len: usize,
    #[arg(long, default_value_t = 400)]
    insert: usize,
    #[arg(long, default_value_t = 50)]
    step: usize,
}

#[derive(Parser)]
struct DiscoverArgs {
    input: PathBuf,
    #[arg(short, long)]
    barcodes: Option<PathBuf>,
    #[arg(short, long, default_value_t = 10000)]
    sample_size: usize,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}

#[derive(Parser)]
struct UmiArgs {
    input: PathBuf,
    #[arg(short, long, default_value = ".")]
    output_dir: PathBuf,
    #[arg(long, default_value = "GTATCGTGTAGAGACTGCGTAGG")]
    fwd_context: String,
    #[arg(long, default_value = "AGTGATCGAGTCAGTGCGAGTG")]
    rev_context: String,
    #[arg(long, default_value = "TTTVVVVTTVVVVTTVVVVTTVVVVTTT")]
    fwd_pattern: String,
    #[arg(long, default_value = "AAABBBBAABBBBAABBBBAABBBBAAA")]
    rev_pattern: String,
    #[arg(long, default_value_t = 4)]
    max_edit: usize,
    #[arg(long, default_value_t = 250)]
    window: usize,
    #[arg(long, default_value_t = 40)]
    min_umi_len: usize,
    #[arg(long, default_value_t = 75)]
    max_umi_len: usize,
    #[arg(long)]
    normalize: bool,
    #[arg(long, default_value_t = 0)]
    min_read_len: usize,
    #[arg(long, default_value_t = 0)]
    max_read_len: usize,
    #[arg(long, default_value_t = 0.0)]
    min_mean_q: f64,
    #[arg(long, default_value_t = 4)]
    min_reads: usize,
    #[arg(long, default_value_t = 80)]
    max_reads: usize,
    #[arg(long)]
    balance_strands: bool,
    #[arg(long, default_value = "approximate")]
    cluster_mode: String,
    #[arg(long, default_value_t = 3)]
    cluster_edit: usize,
    #[arg(long, default_value_t = 0.85)]
    vsearch_identity: f64,
    #[arg(long, default_value = "none")]
    consensus_backend: String,
    #[arg(long)]
    medaka_model: Option<String>,
    #[arg(long, default_value_t = 150)]
    vote_band: usize,
    #[arg(long, default_value = "sample")]
    sample: String,
    #[arg(long, default_value = "amplicon")]
    amplicon: String,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value_t = 0)]
    amplicon_size: usize,
    #[arg(long, default_value_t = 0)]
    size_tolerance: usize,
    #[arg(long)]
    summary: Option<PathBuf>,
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Stats {
            input,
            output,
            threads,
            min_qs,
            min_len,
        } => {
            stats::run_stats(&input, &output, threads, min_qs, min_len)?;
        }
        Commands::Enrichment {
            input,
            bed,
            output,
            threads,
            cm_range,
        } => {
            enrichment::run_enrichment(&input, &bed, &output, threads, cm_range.as_deref())?;
        }
        Commands::Amplicons {
            input,
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
        Commands::Filter(args) => run_filter(args)?,
        Commands::Extract(args) => run_extract(args)?,
        Commands::Split(args) => run_split(args)?,
        Commands::MakePe(args) => run_make_pe(args)?,
        Commands::Discover(args) => run_discover(args)?,
        Commands::Umi(args) => run_umi(args)?,
        Commands::Monitor => {
            anyhow::bail!(
                "The GUI is currently built as `nanomonitor`. Use `cargo run -p nanomonitor` while monitor embedding is wired into `nanostream`."
            );
        }
    }

    Ok(())
}

fn run_filter(args: FilterArgs) -> Result<()> {
    let (min_len, max_len) = parse_len_args(&args.len)?;
    let channel_range = args
        .channel_range
        .as_deref()
        .map(parse_pore_range)
        .transpose()?;
    let time_window = match (args.time_start.as_deref(), args.time_end.as_deref()) {
        (Some(start), Some(end)) => Some(parse_time_window(start, end)?),
        (None, None) => None,
        _ => anyhow::bail!("Use both --time-start and --time-end together"),
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

fn run_extract(args: ExtractArgs) -> Result<()> {
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
        if threads > 1 {
            eprintln!("FASTQ filter currently uses threaded IO only where the backend supports it; split/UMI FASTQ paths use worker threads.");
        }
    } else {
        anyhow::bail!("Unsupported input format");
    }
    Ok(())
}

fn run_split(args: SplitArgs) -> Result<()> {
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
        anyhow::bail!("Unsupported input format");
    }
    Ok(())
}

fn run_make_pe(args: MakePeArgs) -> Result<()> {
    let input_str = args.input.to_string_lossy().to_lowercase();
    if !is_fastq_path(&input_str) {
        anyhow::bail!("make-pe currently supports FASTQ/FASTQ.GZ input only");
    }
    let output_prefix = args
        .output
        .unwrap_or_else(|| default_make_pe_prefix(&args.input, args.insert, args.step));
    fastq::make_pe(
        &args.input,
        &output_prefix,
        args.len,
        args.insert,
        args.step,
    )
}

fn run_discover(args: DiscoverArgs) -> Result<()> {
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
        anyhow::bail!("Unsupported input format");
    }
    Ok(())
}

fn run_umi(args: UmiArgs) -> Result<()> {
    std::fs::create_dir_all(&args.output_dir)?;
    let input_str = args.input.to_string_lossy().to_lowercase();
    if !is_fastq_path(&input_str) {
        anyhow::bail!("UMI detection currently supports FASTQ/FASTQ.GZ input only");
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
        use std::io::Write as _;
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
        _ => anyhow::bail!(
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
        _ => anyhow::bail!("Invalid --len value"),
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
