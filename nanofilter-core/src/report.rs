//! TSV report writers for UMI pipeline output.
//!
//! Three output files:
//! 1. `detected_umis.tsv`        — one row per read, appended across runs.
//! 2. `*_cluster_stats.tsv`       — one row per UMI family.
//! 3. summary TSV                 — one row per sample/amplicon combination.

use anyhow::Result;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::cluster::{ClusterFilterConfig, FamilyStats, UmiFamily};
use crate::umi::{FullReadUmiRecord, UmiConfig};

// ---------------------------------------------------------------------------
// 1. Per-read TSV
// ---------------------------------------------------------------------------

const PER_READ_HEADER: &str = "input_file\tread_id\tstrand\t\
umi_fwd_edit_distance\tumi_rev_edit_distance\t\
umi_fwd_seq\tumi_rev_seq\tcombined_umi\tumi_normalised\tread_length\tinsert_size";

pub fn write_per_read_tsv(records: &[FullReadUmiRecord], path: &Path) -> Result<()> {
    let needs_header = !path.exists() || std::fs::metadata(path)?.len() == 0;
    let file = OpenOptions::new().create(true).append(true).open(path)?;
    let mut w = BufWriter::new(file);

    if needs_header {
        writeln!(w, "{}", PER_READ_HEADER)?;
    }

    for rec in records {
        let r = &rec.result;
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.input_file,
            r.read_id,
            r.strand.map(|s| s.as_str()).unwrap_or("unknown"),
            opt_usize(r.umi_fwd_edit_dist),
            opt_usize(r.umi_rev_edit_dist),
            opt_seq(r.umi_fwd_seq.as_deref()),
            opt_seq(r.umi_rev_seq.as_deref()),
            r.combined_umi.as_deref().unwrap_or("NA"),
            r.umi_normalised.as_deref().unwrap_or("NA"),
            r.read_length,
            opt_usize(r.insert_size),
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// 2. Per-family cluster stats TSV
// ---------------------------------------------------------------------------

const CLUSTER_STATS_HEADER: &str = "id_cluster\tconsensus_umi\tn_cluster_consensus\t\
min_reads_required\tmax_reads_allowed\tbalance_strands\t\
n_fwd\tn_rev\twritten_fwd\twritten_rev\tn\tskipped\twritten\tcluster_written";

pub fn write_cluster_stats_tsv(stats: &[FamilyStats], path: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "{}", CLUSTER_STATS_HEADER)?;

    for s in stats {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.id_cluster,
            s.consensus_umi,
            s.n_cluster_consensus,
            s.min_reads_required,
            s.max_reads_allowed,
            s.balance_strands,
            s.n_fwd,
            s.n_rev,
            s.written_fwd,
            s.written_rev,
            s.n,
            s.skipped,
            s.written,
            s.cluster_written,
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// 3. Summary table
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct SummaryRecord {
    pub sample: String,
    pub amplicon: String,
    pub umi_normalized: bool,
    pub detected_umi_reads: usize,
    pub unique_detected_umis: usize,
    pub total_clusters: usize,
    pub passing_clusters: usize,
    pub failing_clusters: usize,
    pub median_reads_per_cluster: f64,
    pub median_reads_per_passing_cluster: f64,
    pub max_cluster_reads: usize,
    pub reads_in_clusters: usize,
    pub reads_in_passing_clusters: usize,
    pub reads_in_failing_clusters: usize,
    pub reads_written_for_consensus: usize,
    pub reads_trimmed_from_large_clusters: usize,
    pub reads_marked_skipped: usize,
    pub reads_in_small_clusters: usize,
    pub oversized_clusters: usize,
    pub passing_cluster_pct: f64,
    pub reads_used_for_consensus_pct: f64,
    pub reads_in_small_clusters_pct: f64,
    pub mean_forward_fraction: f64,
    pub min_reads_required: usize,
    pub max_reads_allowed: usize,
}

const SUMMARY_HEADER: &str = "sample\tamplicon\tumi_normalized\t\
detected_umi_reads\tunique_detected_umis\ttotal_clusters\tpassing_clusters\t\
failing_clusters\tmedian_reads_per_cluster\tmedian_reads_per_passing_cluster\t\
max_cluster_reads\treads_in_clusters\treads_in_passing_clusters\t\
reads_in_failing_clusters\treads_written_for_consensus\t\
reads_trimmed_from_large_clusters\treads_marked_skipped\t\
reads_in_small_clusters\toversized_clusters\tpassing_cluster_pct\t\
reads_used_for_consensus_pct\treads_in_small_clusters_pct\t\
mean_forward_fraction\tmin_reads_required\tmax_reads_allowed";

/// Build one `SummaryRecord` from the filtered pipeline outputs.
/// Derives "all families" by chaining passing + failing — no separate list needed.
pub fn build_summary(
    sample: &str,
    amplicon: &str,
    all_records: &[FullReadUmiRecord],
    passing_families: &[UmiFamily],
    failing_families: &[UmiFamily],
    cluster_config: &ClusterFilterConfig,
    umi_config: &UmiConfig,
) -> SummaryRecord {
    let detected_umi_reads = all_records.iter().filter(|r| r.result.has_umi()).count();
    let unique_detected_umis = {
        let mut keys: Vec<_> = all_records.iter().filter_map(|r| r.cluster_key()).collect();
        keys.sort_unstable();
        keys.dedup();
        keys.len()
    };

    // Derive all_families by chaining (avoids the caller needing a clone).
    let all_iter = || passing_families.iter().chain(failing_families.iter());

    let total_clusters = passing_families.len() + failing_families.len();
    let passing_clusters = passing_families.len();
    let failing_clusters = failing_families.len();

    let all_sizes: Vec<usize> = all_iter().map(|f| f.total_reads()).collect();
    let pass_sizes: Vec<usize> = passing_families.iter().map(|f| f.total_reads()).collect();

    let median_reads_per_cluster = median(&all_sizes);
    let median_reads_per_passing_cluster = median(&pass_sizes);
    let max_cluster_reads = all_sizes.iter().copied().max().unwrap_or(0);

    let reads_in_clusters: usize = all_iter().map(|f| f.total_reads()).sum();
    let reads_in_passing_clusters: usize = passing_families.iter().map(|f| f.total_reads()).sum();
    let reads_in_failing_clusters: usize = failing_families.iter().map(|f| f.total_reads()).sum();

    let oversized_clusters = all_iter()
        .filter(|f| f.total_reads() > cluster_config.max_reads)
        .count();
    let reads_written_for_consensus: usize = passing_families
        .iter()
        .map(|f| f.total_reads().min(cluster_config.max_reads))
        .sum();
    let reads_trimmed_from_large_clusters =
        reads_in_passing_clusters.saturating_sub(reads_written_for_consensus);
    let reads_marked_skipped = reads_in_failing_clusters;
    let reads_in_small_clusters = reads_in_failing_clusters;

    let passing_cluster_pct = if total_clusters > 0 {
        passing_clusters as f64 / total_clusters as f64 * 100.0
    } else {
        0.0
    };
    let reads_used_for_consensus_pct = if reads_in_clusters > 0 {
        reads_written_for_consensus as f64 / reads_in_clusters as f64 * 100.0
    } else {
        0.0
    };
    let reads_in_small_clusters_pct = if reads_in_clusters > 0 {
        reads_in_small_clusters as f64 / reads_in_clusters as f64 * 100.0
    } else {
        0.0
    };

    let mean_forward_fraction = if reads_in_clusters > 0 {
        let total_fwd: usize = all_iter().map(|f| f.n_fwd).sum();
        total_fwd as f64 / reads_in_clusters as f64
    } else {
        0.0
    };

    SummaryRecord {
        sample: sample.to_string(),
        amplicon: amplicon.to_string(),
        umi_normalized: umi_config.normalize,
        detected_umi_reads,
        unique_detected_umis,
        total_clusters,
        passing_clusters,
        failing_clusters,
        median_reads_per_cluster,
        median_reads_per_passing_cluster,
        max_cluster_reads,
        reads_in_clusters,
        reads_in_passing_clusters,
        reads_in_failing_clusters,
        reads_written_for_consensus,
        reads_trimmed_from_large_clusters,
        reads_marked_skipped,
        reads_in_small_clusters,
        oversized_clusters,
        passing_cluster_pct,
        reads_used_for_consensus_pct,
        reads_in_small_clusters_pct,
        mean_forward_fraction,
        min_reads_required: cluster_config.min_reads,
        max_reads_allowed: cluster_config.max_reads,
    }
}

/// Write (or append) summary records.  Creates the file with a header if it
/// does not yet exist.
pub fn write_summary_tsv(records: &[SummaryRecord], path: &Path) -> Result<()> {
    let needs_header = !path.exists();
    let file = OpenOptions::new().create(true).append(true).open(path)?;
    let mut w = BufWriter::new(file);

    if needs_header {
        writeln!(w, "{}", SUMMARY_HEADER)?;
    }

    for s in records {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.4}\t{}\t{}",
            s.sample,
            s.amplicon,
            s.umi_normalized,
            s.detected_umi_reads,
            s.unique_detected_umis,
            s.total_clusters,
            s.passing_clusters,
            s.failing_clusters,
            s.median_reads_per_cluster,
            s.median_reads_per_passing_cluster,
            s.max_cluster_reads,
            s.reads_in_clusters,
            s.reads_in_passing_clusters,
            s.reads_in_failing_clusters,
            s.reads_written_for_consensus,
            s.reads_trimmed_from_large_clusters,
            s.reads_marked_skipped,
            s.reads_in_small_clusters,
            s.oversized_clusters,
            s.passing_cluster_pct,
            s.reads_used_for_consensus_pct,
            s.reads_in_small_clusters_pct,
            s.mean_forward_fraction,
            s.min_reads_required,
            s.max_reads_allowed,
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn opt_usize(v: Option<usize>) -> String {
    v.map(|n| n.to_string()).unwrap_or_else(|| "NA".to_string())
}

fn opt_seq(v: Option<&[u8]>) -> String {
    v.map(|s| String::from_utf8_lossy(s).to_string())
        .unwrap_or_else(|| "NA".to_string())
}

fn median(data: &[usize]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted = data.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    } else {
        sorted[mid] as f64
    }
}
