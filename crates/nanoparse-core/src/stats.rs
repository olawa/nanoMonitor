//! BAM statistics extraction module

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use nanoseq_core::format::{is_fastq_path, trim_line_ending};
use nanoseq_core::quality::mean_qv_from_fastq_ascii;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use serde::Serialize;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::output::write_output;
use crate::qv::qv_from_record;

#[derive(Debug, Serialize)]
pub struct ReadStats {
    pub id: String,
    pub len: usize,
    pub qs: f32,
    pub acc: Option<f32>,
    pub chrom: Option<String>,
    pub pos: Option<i64>,
    pub is_reverse: bool,
}

#[derive(Debug, Serialize)]
pub struct StatsResult {
    pub total_reads: usize,
    pub passed_reads: usize,
    pub length_stats: DistributionStats,
    pub quality_stats: DistributionStats,
    pub per_read: Vec<ReadStats>,
}

#[derive(Debug, Serialize)]
pub struct DistributionStats {
    pub mean: f64,
    pub median: f64,
    pub std: f64,
    pub min: f64,
    pub max: f64,
}

impl DistributionStats {
    pub fn from_values(values: &[f64]) -> Self {
        if values.is_empty() {
            return Self {
                mean: 0.0,
                median: 0.0,
                std: 0.0,
                min: 0.0,
                max: 0.0,
            };
        }

        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let sum: f64 = sorted.iter().sum();
        let mean = sum / sorted.len() as f64;
        let median = sorted[sorted.len() / 2];
        let min = sorted[0];
        let max = sorted[sorted.len() - 1];

        let variance: f64 =
            sorted.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / sorted.len() as f64;
        let std = variance.sqrt();

        Self {
            mean,
            median,
            std,
            min,
            max,
        }
    }
}

fn read_fastq_stats_from_reader<R: BufRead>(
    mut reader: R,
    min_qs: f32,
    min_len: usize,
) -> Result<StatsResult> {
    let mut total_reads = 0usize;
    let mut per_read = Vec::<ReadStats>::new();

    let mut h = String::new();
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();

    loop {
        h.clear();
        if reader.read_line(&mut h)? == 0 {
            break;
        }
        if !h.starts_with('@') {
            continue;
        }
        seq.clear();
        plus.clear();
        qual.clear();
        if reader.read_line(&mut seq)? == 0
            || reader.read_line(&mut plus)? == 0
            || reader.read_line(&mut qual)? == 0
        {
            break;
        }

        total_reads += 1;
        let seq_trim = trim_line_ending(&seq);
        let qual_trim = trim_line_ending(&qual);
        if seq_trim.len() != qual_trim.len() {
            continue;
        }

        let len = seq_trim.len();
        if len < min_len {
            continue;
        }

        let qs = mean_qv_from_fastq_ascii(qual_trim.as_bytes()) as f32;
        if qs < min_qs {
            continue;
        }

        let id = trim_line_ending(&h)[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        if id.is_empty() {
            continue;
        }

        per_read.push(ReadStats {
            id,
            len,
            qs,
            acc: None,
            chrom: None,
            pos: None,
            is_reverse: false,
        });
    }

    let lengths: Vec<f64> = per_read.iter().map(|r| r.len as f64).collect();
    let qualities: Vec<f64> = per_read.iter().map(|r| r.qs as f64).collect();

    Ok(StatsResult {
        total_reads,
        passed_reads: per_read.len(),
        length_stats: DistributionStats::from_values(&lengths),
        quality_stats: DistributionStats::from_values(&qualities),
        per_read,
    })
}

/// Calculate accuracy from CIGAR if available
fn calculate_accuracy(read: &bam::Record) -> Option<f32> {
    if read.is_unmapped() {
        return None;
    }

    let cigar = read.cigar();
    let mut matches = 0u64;
    let mut total = 0u64;

    for op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(n) | Equal(n) => {
                matches += *n as u64;
                total += *n as u64;
            }
            Diff(n) | Ins(n) | Del(n) => {
                total += *n as u64;
            }
            SoftClip(n) => {
                total += *n as u64;
            }
            _ => {}
        }
    }

    if total > 0 {
        Some((matches as f32 / total as f32) * 100.0)
    } else {
        None
    }
}

pub fn run_stats(
    bam_path: &str,
    output_path: &str,
    threads: usize,
    min_qs: f32,
    min_len: usize,
) -> Result<()> {
    if is_fastq_path(bam_path) {
        log::info!("Reading FASTQ: {}", bam_path);
        let file = File::open(bam_path).context("Failed to open FASTQ file")?;
        let result = if bam_path.to_ascii_lowercase().ends_with(".gz") {
            let gz = MultiGzDecoder::new(file);
            read_fastq_stats_from_reader(BufReader::new(gz), min_qs, min_len)?
        } else {
            read_fastq_stats_from_reader(BufReader::new(file), min_qs, min_len)?
        };
        write_output(&result, output_path)?;
        log::info!("Done. {} reads passed filters.", result.passed_reads);
        return Ok(());
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();

    log::info!("Reading BAM: {}", bam_path);

    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM file")?;
    bam.set_threads(threads).ok();

    // Collect all reads first (for parallel processing)
    let mut all_reads: Vec<bam::Record> = Vec::new();
    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        all_reads.push(record);
    }

    log::info!("Processing {} reads...", all_reads.len());

    // Process in parallel
    let read_stats: Vec<ReadStats> = all_reads
        .par_iter()
        .filter_map(|record| {
            let len = record.seq_len();
            let qs = qv_from_record(record);

            // Apply filters
            if qs < min_qs || len < min_len {
                return None;
            }

            let id = String::from_utf8_lossy(record.qname()).to_string();
            let acc = calculate_accuracy(record);
            let chrom = if record.is_unmapped() {
                None
            } else {
                Some(record.tid().to_string()) // TODO: resolve to actual name
            };
            let pos = if record.is_unmapped() {
                None
            } else {
                Some(record.pos())
            };

            Some(ReadStats {
                id,
                len,
                qs,
                acc,
                chrom,
                pos,
                is_reverse: record.is_reverse(),
            })
        })
        .collect();

    // Calculate distributions
    let lengths: Vec<f64> = read_stats.iter().map(|r| r.len as f64).collect();
    let qualities: Vec<f64> = read_stats.iter().map(|r| r.qs as f64).collect();

    let result = StatsResult {
        total_reads: all_reads.len(),
        passed_reads: read_stats.len(),
        length_stats: DistributionStats::from_values(&lengths),
        quality_stats: DistributionStats::from_values(&qualities),
        per_read: read_stats,
    };

    write_output(&result, output_path)?;

    log::info!("Done. {} reads passed filters.", result.passed_reads);
    Ok(())
}
