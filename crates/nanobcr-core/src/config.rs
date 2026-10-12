use crate::prelim::PrelimConfig;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Everything the engine needs. Defaults mirror nagini's `bcr_abl` assay
/// (`read_filter`: qs >= 20, length 500-3000, both fastq_pass and fastq_fail;
/// `minimap2 -a -x splice -y -Y --MD`), because the BAMs the engine produces are
/// handed to nagini as `aligned_bam`, which neither filters nor realigns.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineConfig {
    /// MinKNOW run directory (or `fastq_pass`); searched recursively for FASTQ.
    pub input_dir: PathBuf,
    /// Everything the engine writes lives here.
    pub out_dir: PathBuf,
    /// Reference FASTA (indexed with `-x splice` on first use) or a prebuilt `.mmi`.
    pub reference: PathBuf,
    pub minimap2: String,
    pub samtools: String,
    pub preset: String,
    /// Passed verbatim to minimap2 after the preset.
    pub minimap2_extra: Vec<String>,
    pub threads: usize,
    pub min_qv: f64,
    pub min_len: usize,
    pub max_len: usize,
    /// Include `fastq_fail` (nagini does by default).
    pub include_fail: bool,
    /// A FASTQ is read only once its mtime is at least this old.
    pub stable_secs: u64,
    pub poll_secs: u64,
    /// Minimum time between alignment cycles (each loads the index).
    pub min_cycle_secs: u64,
    /// Compact the per-barcode inbox into a sorted, indexed BAM this often.
    pub snapshot_secs: u64,
    pub max_reads_per_cycle: usize,
    /// Process what is there, compact, and exit instead of watching.
    pub once: bool,
    /// Run nagini's BCR-ABL scripts on every new snapshot (needs a nagini checkout).
    pub prelim: Option<PrelimConfig>,
}

impl EngineConfig {
    pub fn new(input_dir: PathBuf, out_dir: PathBuf, reference: PathBuf) -> Self {
        Self {
            input_dir,
            out_dir,
            reference,
            minimap2: "minimap2".into(),
            samtools: "samtools".into(),
            preset: "splice".into(),
            minimap2_extra: vec!["-y".into(), "-Y".into(), "--MD".into()],
            threads: 8,
            min_qv: 20.0,
            min_len: 500,
            max_len: 3000,
            include_fail: true,
            stable_secs: 10,
            poll_secs: 2,
            min_cycle_secs: 20,
            snapshot_secs: 600,
            max_reads_per_cycle: 200_000,
            once: false,
            prelim: None,
        }
    }
}
