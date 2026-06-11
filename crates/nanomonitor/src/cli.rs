use clap::{Parser, ValueEnum};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliMode {
    Amplicon,
}

#[derive(Debug, Clone, Parser)]
#[command(name = "nanomonitor", version, about = "nanoMonitor GUI launcher")]
pub struct NanoMonitorCli {
    /// Mode to initialize in
    #[arg(long, value_enum, default_value = "amplicon")]
    pub mode: Option<CliMode>,

    /// Input data path (BAM/FASTQ file or directory)
    #[arg(long)]
    pub input: Option<String>,

    /// Monitor directory for BAM/FASTQ files
    #[arg(long = "monitor-dir")]
    pub monitor_dir: Option<String>,

    /// Reference FASTA path
    #[arg(long)]
    pub reference: Option<String>,

    /// GTF/GFF/BED path
    #[arg(long)]
    pub gtf: Option<String>,

    /// Primers TSV path
    #[arg(long)]
    pub primers: Option<String>,

    /// nanostream executable path
    #[arg(long = "nanostream-bin")]
    pub nanostream_bin: Option<String>,

    /// Start analysis immediately on launch
    #[arg(long)]
    pub start: bool,
}
