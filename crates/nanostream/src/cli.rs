use clap::{Parser, Subcommand};
use nanoparse_core::MatchMode;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "nanostream")]
#[command(
    author,
    version,
    about = "Unified Rust CLI for nanopore BAM/FASTQ workflows"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
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
        /// Read length filter (e.g. 1000 or 1000-3000)
        #[arg(long, default_value = "0")]
        len: String,
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
        /// Optional path to write valid matched FASTQ reads
        #[arg(long)]
        output_fastq: Option<String>,
        /// Optional path to write primer dimers/invalid pairs
        #[arg(long)]
        output_dimers: Option<String>,
        /// If set, split matched FASTQ reads by amplicon name
        #[arg(long, default_value_t = false)]
        split_by_amplicon: bool,
        /// If set, detect internal primers and split chimeric concatemers into individual sub-reads
        #[arg(long, default_value_t = false)]
        split_chimeras: bool,
        /// Maximum reads per amplicon to output, retaining the highest QV/quality reads (0 = no limit)
        #[arg(long, default_value_t = 0)]
        top_n: usize,
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
    /// BCR-ABL1 live mode: watch a MinKNOW run, filter, splice-align, keep per-barcode BAMs
    #[command(subcommand, name = "bcr-abl")]
    BcrAbl(BcrAblCommand),
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

#[derive(Parser, Clone)]
pub struct FilterArgs {
    /// Input BAM/FASTQ(.gz) file
    pub input: PathBuf,
    /// Output BAM/FASTQ(.gz) file
    #[arg(short, long)]
    pub output: PathBuf,
    #[arg(long, default_value_t = 20.0)]
    pub qv: f64,
    #[arg(long, num_args = 1.., value_delimiter = ' ')]
    pub len: Vec<String>,
    #[arg(long, default_value = "auto")]
    pub output_format: String,
    #[arg(long)]
    pub channel_range: Option<String>,
    #[arg(long)]
    pub time_start: Option<String>,
    #[arg(long)]
    pub time_end: Option<String>,
    /// Number of threads for threaded BAM IO and parallel FASTQ operations where supported
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    /// Path to ONT sequencing summary file for relative time-filtering baseline discovery
    #[arg(long)]
    pub sequencing_summary: Option<PathBuf>,
}

#[derive(Parser, Clone)]
pub struct ExtractArgs {
    pub input: PathBuf,
    #[arg(short, long)]
    pub output: PathBuf,
    #[arg(long)]
    pub channel_range: String,
    #[arg(long, default_value = "auto")]
    pub output_format: String,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
}

#[derive(Parser, Clone)]
pub struct SplitArgs {
    pub input: PathBuf,
    #[arg(short, long)]
    pub barcodes: PathBuf,
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,
    #[arg(short, long, default_value_t = 1)]
    pub mismatches: usize,
    #[arg(short, long, default_value_t = 1000)]
    pub search_dist: usize,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, default_value_t = 0)]
    pub auto_discover: usize,
    #[arg(long)]
    pub fast: bool,
}

#[derive(Parser, Clone)]
pub struct MakePeArgs {
    /// Input BAM/FASTQ(.gz) file
    pub input: PathBuf,
    /// Output prefix or path. Use '-' for stdout (interleaved)
    #[arg(short, long)]
    pub output: Option<PathBuf>,
    /// Length of synthetic PE reads (default: 150)
    #[arg(short = 'l', long, default_value_t = 150)]
    pub len: usize,
    /// Insert size (default: 2 * len, e.g. 300 for 150 bp reads)
    #[arg(short = 'i', long)]
    pub insert: Option<usize>,
    /// Step/shift size between consecutive fragments (default: insert, for 100% full coverage tiling)
    #[arg(short = 's', long, alias = "shift")]
    pub step: Option<usize>,
    /// Write R1 and R2 interleaved into a single FASTQ file or stdout
    #[arg(long, default_value_t = false)]
    pub interleaved: bool,
    /// Number of threads for BAM decompression
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Clone)]
pub struct DiscoverArgs {
    pub input: PathBuf,
    #[arg(short, long)]
    pub barcodes: Option<PathBuf>,
    #[arg(short, long, default_value_t = 10000)]
    pub sample_size: usize,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
}

#[derive(Parser, Clone)]
pub struct UmiArgs {
    pub input: PathBuf,
    #[arg(short, long, default_value = ".")]
    pub output_dir: PathBuf,
    #[arg(long, default_value = "GTATCGTGTAGAGACTGCGTAGG")]
    pub fwd_context: String,
    #[arg(long, default_value = "AGTGATCGAGTCAGTGCGAGTG")]
    pub rev_context: String,
    #[arg(long, default_value = "TTTVVVVTTVVVVTTVVVVTTVVVVTTT")]
    pub fwd_pattern: String,
    #[arg(long, default_value = "AAABBBBAABBBBAABBBBAABBBBAAA")]
    pub rev_pattern: String,
    #[arg(long, default_value_t = 4)]
    pub max_edit: usize,
    #[arg(long, default_value_t = 250)]
    pub window: usize,
    #[arg(long, default_value_t = 40)]
    pub min_umi_len: usize,
    #[arg(long, default_value_t = 75)]
    pub max_umi_len: usize,
    #[arg(long)]
    pub normalize: bool,
    #[arg(long, default_value_t = 0)]
    pub min_read_len: usize,
    #[arg(long, default_value_t = 0)]
    pub max_read_len: usize,
    #[arg(long, default_value_t = 0.0)]
    pub min_mean_q: f64,
    #[arg(long, default_value_t = 4)]
    pub min_reads: usize,
    #[arg(long, default_value_t = 80)]
    pub max_reads: usize,
    #[arg(long)]
    pub balance_strands: bool,
    #[arg(long, default_value = "approximate")]
    pub cluster_mode: String,
    #[arg(long, default_value_t = 3)]
    pub cluster_edit: usize,
    #[arg(long, default_value_t = 0.85)]
    pub vsearch_identity: f64,
    #[arg(long, default_value = "none")]
    pub consensus_backend: String,
    #[arg(long)]
    pub medaka_model: Option<String>,
    #[arg(long, default_value_t = 150)]
    pub vote_band: usize,
    #[arg(long, default_value = "sample")]
    pub sample: String,
    #[arg(long, default_value = "amplicon")]
    pub amplicon: String,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, default_value_t = 0)]
    pub amplicon_size: usize,
    #[arg(long, default_value_t = 0)]
    pub size_tolerance: usize,
    #[arg(long)]
    pub summary: Option<PathBuf>,
}

#[derive(Subcommand, Clone)]
pub enum BcrAblCommand {
    /// Watch a MinKNOW directory and keep sorted, indexed per-barcode BAMs up to date
    Watch(BcrWatchArgs),
    /// Fold every barcode's pending segments into its sorted, indexed BAM (engine need not be running)
    Snapshot(BcrSnapshotArgs),
    /// Run nagini's BCR-ABL scripts (fusion isoform, kinase-domain mutations, phasing) on the current snapshots
    Analyze(BcrAnalyzeArgs),
    /// Gather the final BAMs into the layout nagini reads (`barcodeNN.bam` + a barcode map); the engine must be stopped
    Handoff(BcrHandoffArgs),
    /// Start, check, follow or stop the final nagini run on a hand-off
    #[command(subcommand)]
    Nagini(BcrNaginiCommand),
    /// Compare the preliminary results with nagini's final results
    Compare(BcrCompareArgs),
    /// Feed a finished run into a directory the way MinKNOW would (for testing)
    Replay(BcrReplayArgs),
}

#[derive(Parser, Clone)]
pub struct BcrWatchArgs {
    /// MinKNOW run directory (or its fastq_pass directory)
    pub input: PathBuf,
    /// Output directory (state.json, provenance.json, per-barcode BAMs)
    #[arg(short, long)]
    pub out: PathBuf,
    /// Reference FASTA (indexed once with `-x splice`) or a prebuilt .mmi
    #[arg(short, long)]
    pub reference: PathBuf,
    /// Minimum mean read quality (dorado qs)
    #[arg(long, default_value_t = 20.0)]
    pub qv: f64,
    /// Read length range, MIN-MAX
    #[arg(long, default_value = "500-3000")]
    pub len: String,
    /// Only read fastq_pass (nagini includes fastq_fail by default)
    #[arg(long, default_value_t = false)]
    pub pass_only: bool,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    /// A FASTQ is read once it has not changed for this many seconds
    #[arg(long, default_value_t = 10)]
    pub stable_secs: u64,
    #[arg(long, default_value_t = 2)]
    pub poll_secs: u64,
    /// Minimum seconds between alignment cycles
    #[arg(long, default_value_t = 20)]
    pub cycle_secs: u64,
    /// Seconds between compactions into sorted, indexed BAMs (each new snapshot is analysed when --nagini is set)
    #[arg(long, default_value_t = 300)]
    pub snapshot_secs: u64,
    /// Process what is there, compact, and exit
    #[arg(long, default_value_t = false)]
    pub once: bool,
    #[arg(long, default_value = "minimap2")]
    pub minimap2: String,
    #[arg(long, default_value = "samtools")]
    pub samtools: String,
    /// Extra minimap2 arguments after the preset (space separated)
    #[arg(long, default_value = "-y -Y --MD")]
    pub minimap2_extra: String,
    /// nagini checkout: run its BCR-ABL scripts on every new snapshot (preliminary analysis)
    #[arg(long)]
    pub nagini: Option<PathBuf>,
    /// Python with pysam for the nagini scripts (default: nagini's .venv, else python3)
    #[arg(long)]
    pub python: Option<String>,
    /// Reference FASTA for the Branford caller's k-mer check (default: --reference when it is a FASTA)
    #[arg(long)]
    pub fasta: Option<PathBuf>,
    /// Reference build name for nagini's resources
    #[arg(long, default_value = "hg19")]
    pub build: String,
}

#[derive(Parser, Clone)]
pub struct BcrAnalyzeArgs {
    /// Engine output directory
    pub out: PathBuf,
    /// nagini checkout
    #[arg(long)]
    pub nagini: PathBuf,
    #[arg(long)]
    pub python: Option<String>,
    #[arg(long)]
    pub fasta: Option<PathBuf>,
    #[arg(long, default_value = "hg19")]
    pub build: String,
    /// Only this barcode
    #[arg(long)]
    pub barcode: Option<String>,
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
}

#[derive(Parser, Clone)]
pub struct BcrSnapshotArgs {
    /// Engine output directory
    pub out: PathBuf,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, default_value = "samtools")]
    pub samtools: String,
}

#[derive(Parser, Clone)]
pub struct BcrReplayArgs {
    /// Finished run directory to replay
    pub src: PathBuf,
    /// Directory to write into (point `bcr-abl watch` at it)
    pub dest: PathBuf,
    /// Seconds between files
    #[arg(long, default_value_t = 4.0)]
    pub interval: f64,
    /// Expose each file half-written for a moment
    #[arg(long, default_value_t = false)]
    pub partial: bool,
    /// Write final_summary_replay.txt at the end, as MinKNOW does
    #[arg(long, default_value_t = false)]
    pub finish: bool,
}

#[derive(Parser, Clone)]
pub struct BcrHandoffArgs {
    /// Engine output directory
    pub out: PathBuf,
    /// Barcode map CSV/TSV (barcode, sample, type, assay, library_type)
    #[arg(long)]
    pub barcode_map: Option<PathBuf>,
    /// Write a barcode map template with every barcode the engine saw, then stop
    #[arg(long)]
    pub init_map: Option<PathBuf>,
    /// Hand-off folder (default: <out>/handoff)
    #[arg(long)]
    pub handoff: Option<PathBuf>,
    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,
    #[arg(long, default_value = "samtools")]
    pub samtools: String,
}

#[derive(Subcommand, Clone)]
pub enum BcrNaginiCommand {
    /// Ask snakemake whether the workflow can be built here (`--dry-run`), without running it
    Check(BcrNaginiArgs),
    /// Start the run in the background
    Start(BcrNaginiArgs),
    /// Show progress and outcome
    Status {
        /// nagini analysis folder
        analysis: PathBuf,
    },
    /// Interrupt a running pipeline
    Stop {
        analysis: PathBuf,
    },
}

#[derive(Parser, Clone)]
pub struct BcrNaginiArgs {
    /// Engine output directory (its hand-off folder is used)
    pub out: PathBuf,
    #[arg(long)]
    pub nagini: PathBuf,
    /// nagini analysis folder (default: <out>/nagini_run)
    #[arg(long)]
    pub analysis: Option<PathBuf>,
    #[arg(long)]
    pub handoff: Option<PathBuf>,
    #[arg(long, default_value = "workstation")]
    pub site: String,
    #[arg(long, default_value = "hg19")]
    pub build: String,
    #[arg(long)]
    pub cores: Option<usize>,
    #[arg(long, default_value_t = false)]
    pub enable_rindels: bool,
}

#[derive(Parser, Clone)]
pub struct BcrCompareArgs {
    /// Engine output directory
    pub out: PathBuf,
    /// nagini analysis folder (default: <out>/nagini_run)
    #[arg(long)]
    pub analysis: Option<PathBuf>,
    #[arg(long)]
    pub handoff: Option<PathBuf>,
}
