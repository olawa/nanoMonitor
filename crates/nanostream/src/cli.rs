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
    pub input: PathBuf,
    #[arg(short, long)]
    pub output: Option<PathBuf>,
    #[arg(long, default_value_t = 150)]
    pub len: usize,
    #[arg(long, default_value_t = 400)]
    pub insert: usize,
    #[arg(long, default_value_t = 50)]
    pub step: usize,
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
