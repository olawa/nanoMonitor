//! Core amplicon matching logic
//!
//! Supports multiple matching modes:
//! - Semi-global alignment using triple_accel SIMD (default)
//! - Coordinate-based assignment for aligned BAMs (very fast)

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use nanoseq_core::format::{is_fastq_path, trim_line_ending};
use nanoseq_core::quality::mean_qv_from_fastq_ascii;
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use triple_accel::levenshtein::levenshtein_simd_k;

use crate::output::write_output;
use crate::primers::{load_primers, KmerIndex, Primer};
use crate::qv::qv_from_record;
use crate::MatchMode;

#[derive(Debug, Serialize)]
pub struct AmpliconResult {
    pub amplicons: HashMap<String, AmpliconStats>,
    pub chimera_count: usize,
    pub unmatched_count: usize,
    pub total_reads: usize,
    pub rescued_count: usize,
    pub distributions: ReadDistributions,
}

#[derive(Debug, Serialize)]
pub struct ReadDistributions {
    pub length_bins: Vec<DistributionBin>,
    pub qs_bins: Vec<DistributionBin>,
    pub accuracy_bins: Vec<DistributionBin>,
    pub length_median: f64,
    pub qs_mode: f64,
    pub accuracy_mode: f64,
}

#[derive(Debug, Serialize)]
pub struct DistributionBin {
    pub start: f64,
    pub end: f64,
    pub count: usize,
}

#[derive(Debug, Serialize, Default, Clone)]
pub struct AmpliconStats {
    pub count: usize,
    pub median_length: usize,
    pub std_length: f64,
    pub avg_qs: f32,
    pub chrom: Option<String>,
    pub start: Option<i64>,
    pub end: Option<i64>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub read_ids: Vec<String>,
    #[serde(skip)]
    lengths: Vec<usize>,
    #[serde(skip)]
    qualities: Vec<f32>,
}

impl AmpliconStats {
    pub fn finalize(&mut self) {
        if !self.lengths.is_empty() {
            let mut sorted = self.lengths.clone();
            sorted.sort();
            self.median_length = sorted[sorted.len() / 2];

            let mean: f64 = sorted.iter().map(|&x| x as f64).sum::<f64>() / sorted.len() as f64;
            let variance: f64 = sorted
                .iter()
                .map(|&x| (x as f64 - mean).powi(2))
                .sum::<f64>()
                / sorted.len() as f64;
            self.std_length = variance.sqrt();
        }

        if !self.qualities.is_empty() {
            self.avg_qs = self.qualities.iter().sum::<f32>() / self.qualities.len() as f32;
        }
    }
}

/// Match result for a single read
#[derive(Debug, Clone)]
struct MatchResult {
    // Optimization: Delay read_id allocation
    read_id_index: usize,
    amplicon_name: Option<String>,
    start_primer: Option<String>, // For debug stats
    end_primer: Option<String>,   // For debug stats
    is_chimera: bool,
    length: usize,
    quality: f32,
    chrom: Option<String>,
    start: Option<i64>,
    end: Option<i64>,
}

#[derive(Debug, Clone)]
struct FastqRead {
    id: String,
    header: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

#[derive(Debug, Clone)]
struct PendingUnassigned {
    id: String,
    res: MatchResult,
}

/// Semi-global alignment with early exit optimization using triple_accel
fn semi_global_align(pattern: &[u8], text: &[u8], max_dist: usize) -> Option<u32> {
    if pattern.is_empty() || text.is_empty() {
        return None;
    }

    let pattern_len = pattern.len();
    let mut best_dist: Option<u32> = None;

    let max_start = if text.len() > pattern_len {
        text.len() - pattern_len + max_dist
    } else {
        max_dist
    };

    for start in 0..=max_start.min(text.len().saturating_sub(1)) {
        let end = (start + pattern_len + max_dist).min(text.len());
        if end <= start {
            continue;
        }

        let text_slice = &text[start..end];

        if let Some(dist) = levenshtein_simd_k(pattern, text_slice, max_dist as u32) {
            if dist == 0 {
                return Some(0); // Perfect match - early exit
            }
            if best_dist.is_none() || dist < best_dist.unwrap() {
                best_dist = Some(dist);
            }
        }
    }

    best_dist
}

/// Match using sequence alignment and K-mer filtering
fn match_read_sequence(
    record: &bam::Record,
    read_index: usize,
    primers: &[Primer],
    kmer_index: Option<&KmerIndex>,
    end_length: usize,
    max_edit_dist: usize,
    early_exit: bool,
    detect_chimeras: bool,
    primer_tolerance: i64,
) -> MatchResult {
    // Optimization: Decode full sequence once using as_bytes()
    let full_seq_bytes = record.seq().as_bytes();
    let length = full_seq_bytes.len();
    let quality = qv_from_record(record);

    let tid = record.tid();
    let start = record.pos();
    let end = record.pos() + record.seq_len() as i64;

    // Keep chromosome id as numeric tid string to avoid header pointer issues.
    let chrom_name = if tid >= 0 {
        Some(tid.to_string())
    } else {
        None
    };
    let chrom = chrom_name.clone(); // Keep original clone for result

    let extract_len = std::cmp::min(length, end_length);
    let start_bytes = &full_seq_bytes[..extract_len];
    let end_start_idx = length.saturating_sub(extract_len);
    let end_bytes = &full_seq_bytes[end_start_idx..];

    let mut best_start_primer: Option<&str> = None;
    let mut best_start_dist: u32 = (max_edit_dist + 1) as u32;

    let mut best_end_primer: Option<&str> = None;
    let mut best_end_dist: u32 = (max_edit_dist + 1) as u32;

    // START Primer Matching
    // ---------------------

    // FAST PATH: Check Coordinates Coverage
    // If we have mapped coordinates and they match a primer region, we can skip expensive alignment.
    if primer_tolerance > 0 {
        if let Some(ref c) = chrom {
            let r_c_norm = c.replace("chr", "");

            // Check all primers for start/end match
            // We need to iterate primers.
            // Note: 'primers' slice contains individual primers, usually Fwd and Rev.
            // If a primer has coordinates, we check them.

            for primer in primers {
                if let (Some(p_c), Some(p_start), Some(p_end)) =
                    (&primer.chrom, primer.start, primer.end)
                {
                    let p_c_norm = p_c.replace("chr", "");
                    if r_c_norm == p_c_norm {
                        // Check Start
                        let d_start = (start - p_start).abs();
                        let d_end = (end - p_end).abs();

                        // If this primer matches the START of the read
                        if d_start <= primer_tolerance {
                            // This is a candidate for start primer
                            // If we assume perfect coordinate trust, we can set it and break.
                            best_start_primer = Some(&primer.name);
                            best_start_dist = 0; // "Perfect" match by coords
                        }

                        // If this primer matches the END of the read
                        if d_end <= primer_tolerance {
                            best_end_primer = Some(&primer.name);
                            best_end_dist = 0;
                        }
                    }
                }
            }

            // If we found both by coordinates, we can return early!
            // Note: This logic assumes individual primers have start/end.
            // If the primers file defines Amplicons as "PrimerA:chr:start-end"
            // Then usually PrimerA covers the whole region.
            // But if we defined Fwd and Rev primers separately with coord...
            // The current `primers.rs` loads each line as a Primer.
            // If the user provided "AmpliconName -> Sequence -> Region",
            // Then we have 1 Primer struct per line.
            // If that struct has start/end, and read matches BOTH start/end to this single struct,
            // then this single struct represents the whole amplicon.
            // This is the common case for the TSV format "AmpliconID Seq Region".

            if let (Some(s), Some(e)) = (best_start_primer, best_end_primer) {
                if s == e {
                    // Single amplicon entry covers both ends (Standard Amplicon definition)
                    return MatchResult {
                        read_id_index: read_index,
                        amplicon_name: Some(s.to_string()),
                        start_primer: Some(s.to_string()),
                        end_primer: Some(e.to_string()), // Same name
                        length: length,
                        quality,
                        chrom: chrom.clone(),
                        start: Some(start),
                        end: Some(end),
                        is_chimera: false,
                    };
                } else {
                    // Distinct primers matched (Forward and Reverse defined separately?)
                    // If so, combine them.
                    // But typically `nanoparse` expects to find pairs.
                    // If we found a pair by coords, proceed.
                    let name = if s < e {
                        format!("{}-{}", s, e)
                    } else {
                        format!("{}-{}", e, s)
                    };
                    return MatchResult {
                        read_id_index: read_index,
                        amplicon_name: Some(name),
                        start_primer: Some(s.to_string()),
                        end_primer: Some(e.to_string()),
                        length: length,
                        quality,
                        chrom: chrom.clone(),
                        start: Some(start),
                        end: Some(end),
                        is_chimera: false,
                    };
                }
            }
        }
    }

    // SEQUENCE MATCHING (Original Logic)
    if let Some(idx) = kmer_index {
        let start_str = unsafe { std::str::from_utf8_unchecked(start_bytes) };
        let candidates = idx.find_candidates(start_str, 0.1);

        // Fallback to all primers if candidates empty? No, rely on index.
        for &p_idx in &candidates {
            let primer = &primers[p_idx];
            let primer_bytes = primer.sequence.as_bytes();
            if let Some(dist) = semi_global_align(primer_bytes, start_bytes, max_edit_dist) {
                if dist < best_start_dist {
                    best_start_dist = dist;
                    best_start_primer = Some(&primer.name);
                    if early_exit && dist == 0 {
                        break;
                    }
                }
            }
        }
    } else {
        // Linear scan (all primers)
        for primer in primers {
            let primer_bytes = primer.sequence.as_bytes();
            if let Some(dist) = semi_global_align(primer_bytes, start_bytes, max_edit_dist) {
                if dist < best_start_dist {
                    best_start_dist = dist;
                    best_start_primer = Some(&primer.name);
                    if early_exit && dist == 0 {
                        break;
                    }
                }
            }
        }
    }

    // END Primer Matching (Reverse Complement)
    // ----------------------------------------
    if let Some(idx) = kmer_index {
        let end_str = unsafe { std::str::from_utf8_unchecked(end_bytes) };
        let candidates = idx.find_candidates(end_str, 0.1);

        for &p_idx in &candidates {
            let primer = &primers[p_idx];
            let primer_rc_bytes = primer.sequence_rc.as_bytes();
            if let Some(dist) = semi_global_align(primer_rc_bytes, end_bytes, max_edit_dist) {
                if dist < best_end_dist {
                    best_end_dist = dist;
                    best_end_primer = Some(&primer.name);
                }
            }
        }
    } else {
        for primer in primers {
            let primer_rc_bytes = primer.sequence_rc.as_bytes();
            if let Some(dist) = semi_global_align(primer_rc_bytes, end_bytes, max_edit_dist) {
                if dist < best_end_dist {
                    best_end_dist = dist;
                    best_end_primer = Some(&primer.name);
                }
            }
        }
    }

    // Chimera detection
    let mut is_chimera = false;
    if detect_chimeras && length > end_length * 2 {
        let middle_start = end_length;
        let middle_end = length.saturating_sub(end_length);
        if middle_end > middle_start {
            let middle_bytes = &full_seq_bytes[middle_start..middle_end];
            for primer in primers {
                let p_bytes = primer.sequence.as_bytes();
                let rc_bytes = primer.sequence_rc.as_bytes();
                let found = semi_global_align(p_bytes, middle_bytes, max_edit_dist).is_some()
                    || semi_global_align(rc_bytes, middle_bytes, max_edit_dist).is_some();
                if found {
                    is_chimera = true;
                    break;
                }
            }
        }
    }

    let mut amplicon_name = match (best_start_primer, best_end_primer) {
        (Some(p1), Some(p2)) => {
            let mut names = vec![p1, p2];
            names.sort();
            Some(format!("{}-{}", names[0], names[1]))
        }
        _ => None,
    };

    // Fallback: Fuzzy Coordinate Matching
    // If we failed to find a primer pair by sequence, check alignment coordinates
    if amplicon_name.is_none() && !record.is_unmapped() && primer_tolerance > 0 {
        if let Some(r_chrom) = &chrom_name {
            let mut best_coord_dist = primer_tolerance * 2 + 1; // Sum of distances
            let mut best_coord_primer_name = None;

            let r_c_norm = r_chrom.replace("chr", ""); // Normalize for loose matching

            for primer in primers {
                if let (Some(p_chrom), Some(p_start), Some(p_end)) =
                    (&primer.chrom, primer.start, primer.end)
                {
                    let p_c_norm = p_chrom.replace("chr", "");

                    if r_c_norm == p_c_norm {
                        let d_start = (start - p_start).abs();
                        let d_end = (end - p_end).abs();

                        if d_start <= primer_tolerance && d_end <= primer_tolerance {
                            let total_dist = d_start + d_end;
                            if total_dist < best_coord_dist {
                                best_coord_dist = total_dist;
                                best_coord_primer_name = Some(&primer.name);
                            }
                        }
                    }
                }
            }

            if let Some(name) = best_coord_primer_name {
                amplicon_name = Some(name.clone());
            }
        }
    }

    MatchResult {
        read_id_index: read_index,
        amplicon_name,
        start_primer: best_start_primer.map(|s| s.to_string()),
        end_primer: best_end_primer.map(|s| s.to_string()),
        is_chimera,
        length,
        quality,
        chrom,
        start: Some(start),
        end: Some(end),
    }
}

/// Match using mapping coordinates
fn match_read_coords(record: &bam::Record, index: usize) -> MatchResult {
    let length = record.seq_len();
    let quality = qv_from_record(record);

    if record.is_unmapped() {
        return MatchResult {
            read_id_index: index,
            amplicon_name: None,
            start_primer: None,
            end_primer: None,
            is_chimera: false,
            length,
            quality,
            chrom: None,
            start: None,
            end: None,
        };
    }

    let chrom = format!("{}", record.tid());
    let start = record.pos();
    let end = record.pos() + record.seq_len() as i64;

    // For now, use position-based naming as placeholder
    let amplicon_name = Some(format!("{}:{}-{}", chrom, start, end));

    MatchResult {
        read_id_index: index,
        amplicon_name,
        start_primer: None,
        end_primer: None,
        is_chimera: false,
        length,
        quality,
        chrom: Some(chrom),
        start: Some(start),
        end: Some(end),
    }
}

fn match_fastq_read(
    read: &FastqRead,
    read_index: usize,
    primers: &[Primer],
    kmer_index: Option<&KmerIndex>,
    end_length: usize,
    max_edit_dist: usize,
    early_exit: bool,
    detect_chimeras: bool,
) -> MatchResult {
    let full_seq_bytes = &read.seq;
    let length = full_seq_bytes.len();
    let quality = 0.0;

    let extract_len = std::cmp::min(length, end_length);
    let start_bytes = &full_seq_bytes[..extract_len];
    let end_start_idx = length.saturating_sub(extract_len);
    let end_bytes = &full_seq_bytes[end_start_idx..];

    let mut best_start_primer: Option<&str> = None;
    let mut best_start_dist: u32 = (max_edit_dist + 1) as u32;
    let mut best_end_primer: Option<&str> = None;
    let mut best_end_dist: u32 = (max_edit_dist + 1) as u32;

    if let Some(idx) = kmer_index {
        let start_str = unsafe { std::str::from_utf8_unchecked(start_bytes) };
        let candidates = idx.find_candidates(start_str, 0.1);
        for &p_idx in &candidates {
            let primer = &primers[p_idx];
            let primer_bytes = primer.sequence.as_bytes();
            if let Some(dist) = semi_global_align(primer_bytes, start_bytes, max_edit_dist) {
                if dist < best_start_dist {
                    best_start_dist = dist;
                    best_start_primer = Some(&primer.name);
                    if early_exit && dist == 0 {
                        break;
                    }
                }
            }
        }
    } else {
        for primer in primers {
            let primer_bytes = primer.sequence.as_bytes();
            if let Some(dist) = semi_global_align(primer_bytes, start_bytes, max_edit_dist) {
                if dist < best_start_dist {
                    best_start_dist = dist;
                    best_start_primer = Some(&primer.name);
                    if early_exit && dist == 0 {
                        break;
                    }
                }
            }
        }
    }

    if let Some(idx) = kmer_index {
        let end_str = unsafe { std::str::from_utf8_unchecked(end_bytes) };
        let candidates = idx.find_candidates(end_str, 0.1);
        for &p_idx in &candidates {
            let primer = &primers[p_idx];
            let primer_rc_bytes = primer.sequence_rc.as_bytes();
            if let Some(dist) = semi_global_align(primer_rc_bytes, end_bytes, max_edit_dist) {
                if dist < best_end_dist {
                    best_end_dist = dist;
                    best_end_primer = Some(&primer.name);
                }
            }
        }
    } else {
        for primer in primers {
            let primer_rc_bytes = primer.sequence_rc.as_bytes();
            if let Some(dist) = semi_global_align(primer_rc_bytes, end_bytes, max_edit_dist) {
                if dist < best_end_dist {
                    best_end_dist = dist;
                    best_end_primer = Some(&primer.name);
                }
            }
        }
    }

    let mut is_chimera = false;
    if detect_chimeras && length > end_length * 2 {
        let middle_start = end_length;
        let middle_end = length.saturating_sub(end_length);
        if middle_end > middle_start {
            let middle_bytes = &full_seq_bytes[middle_start..middle_end];
            for primer in primers {
                let p_bytes = primer.sequence.as_bytes();
                let rc_bytes = primer.sequence_rc.as_bytes();
                let found = semi_global_align(p_bytes, middle_bytes, max_edit_dist).is_some()
                    || semi_global_align(rc_bytes, middle_bytes, max_edit_dist).is_some();
                if found {
                    is_chimera = true;
                    break;
                }
            }
        }
    }

    let amplicon_name = match (best_start_primer, best_end_primer) {
        (Some(p1), Some(p2)) => {
            let mut names = vec![p1, p2];
            names.sort();
            Some(format!("{}-{}", names[0], names[1]))
        }
        _ => None,
    };

    MatchResult {
        read_id_index: read_index,
        amplicon_name,
        start_primer: best_start_primer.map(|s| s.to_string()),
        end_primer: best_end_primer.map(|s| s.to_string()),
        is_chimera,
        length,
        quality,
        chrom: None,
        start: None,
        end: None,
    }
}


fn phred_to_accuracy_pct(qs: f64) -> f64 {
    let p_err = 10f64.powf(-qs / 10.0);
    ((1.0 - p_err) * 100.0).clamp(0.0, 100.0)
}

fn build_histogram(values: &[f64], start: f64, end: f64, width: f64) -> Vec<DistributionBin> {
    if values.is_empty() || width <= 0.0 || end <= start {
        return Vec::new();
    }
    let n_bins = (((end - start) / width).ceil() as usize).max(1);
    let mut counts = vec![0usize; n_bins];

    for &v in values {
        if v < start || v > end {
            continue;
        }
        let mut idx = ((v - start) / width).floor() as usize;
        if idx >= n_bins {
            idx = n_bins - 1;
        }
        counts[idx] += 1;
    }

    counts
        .into_iter()
        .enumerate()
        .map(|(i, count)| DistributionBin {
            start: start + i as f64 * width,
            end: start + (i + 1) as f64 * width,
            count,
        })
        .collect()
}

fn median_from_sorted(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) * 0.5
    } else {
        values[mid]
    }
}

fn mode_from_bins(bins: &[DistributionBin]) -> f64 {
    if bins.is_empty() {
        return 0.0;
    }
    if let Some(bin) = bins.iter().max_by_key(|b| b.count) {
        (bin.start + bin.end) * 0.5
    } else {
        0.0
    }
}


fn build_distributions_from_values(lengths: &[f64], qs_vals: &[f64]) -> ReadDistributions {
    let mut lengths = lengths.to_vec();
    let qs_vals = qs_vals.to_vec();
    let mut acc_vals = Vec::with_capacity(qs_vals.len());
    for &qs in &qs_vals {
        acc_vals.push(phred_to_accuracy_pct(qs));
    }

    let max_len = lengths.iter().copied().fold(0.0_f64, f64::max).max(1000.0);
    let len_width = if max_len <= 3000.0 {
        50.0
    } else if max_len <= 12000.0 {
        100.0
    } else {
        250.0
    };
    let len_end = (max_len / len_width).ceil() * len_width;

    let length_bins = build_histogram(&lengths, 0.0, len_end, len_width);
    let qs_bins = build_histogram(&qs_vals, 0.0, 50.0, 1.0);
    let accuracy_bins = build_histogram(&acc_vals, 90.0, 100.0, 0.1);
    let qs_mode = mode_from_bins(&qs_bins);
    let accuracy_mode = mode_from_bins(&accuracy_bins);

    lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());

    ReadDistributions {
        length_bins,
        qs_bins,
        accuracy_bins,
        length_median: median_from_sorted(&lengths),
        qs_mode,
        accuracy_mode,
    }
}

fn is_fwd_rev_pair(p1: &str, p2: &str) -> bool {
    let p1_lower = p1.to_lowercase();
    let p2_lower = p2.to_lowercase();
    let p1_is_fwd = p1_lower.contains("fwd") || p1_lower.contains("forward");
    let p1_is_rev = p1_lower.contains("rev") || p1_lower.contains("reverse");
    let p2_is_fwd = p2_lower.contains("fwd") || p2_lower.contains("forward");
    let p2_is_rev = p2_lower.contains("rev") || p2_lower.contains("reverse");
    (p1_is_fwd && p2_is_rev) || (p1_is_rev && p2_is_fwd)
}

fn parse_qs_from_header(header: &str) -> Option<f32> {
    for part in header.split_whitespace().skip(1) {
        if let Some((key, value)) = part.split_once('=') {
            if key == "qs" {
                if let Ok(v) = value.parse::<f32>() {
                    return Some(v);
                }
            }
        }
    }
    None
}

fn get_split_path(base_path: &str, suffix: &str) -> String {
    let path = std::path::Path::new(base_path);
    let parent = path.parent().unwrap_or_else(|| std::path::Path::new(""));
    let file_name = path.file_name().and_then(|n| n.to_str()).unwrap_or("output");
    let (stem, ext) = if file_name.ends_with(".fastq.gz") {
        (&file_name[..file_name.len() - 9], ".fastq.gz")
    } else if file_name.ends_with(".fq.gz") {
        (&file_name[..file_name.len() - 6], ".fq.gz")
    } else if file_name.ends_with(".fastq") {
        (&file_name[..file_name.len() - 6], ".fastq")
    } else if file_name.ends_with(".fq") {
        (&file_name[..file_name.len() - 3], ".fq")
    } else {
        (file_name, "")
    };
    let new_name = format!("{}_{}{}", stem, suffix, ext);
    parent.join(new_name).to_string_lossy().to_string()
}

fn create_fastq_writer(path: &str) -> Result<Box<dyn std::io::Write>> {
    let file = File::create(path).with_context(|| format!("Failed to create output file: {}", path))?;
    if path.to_ascii_lowercase().ends_with(".gz") {
        Ok(Box::new(GzEncoder::new(
            std::io::BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Box::new(std::io::BufWriter::new(file)))
    }
}

pub fn run_amplicons(
    bam_path: &str,
    primers_path: &str,
    threads: usize,
    mode: MatchMode,
    max_edit_dist: usize,
    end_length: usize,
    print_summary: bool,
    primer_tolerance: i64,
    min_qs: f32,
    len_range: (usize, usize),
    max_reads: usize,
    duplex_only: bool,
    reference: Option<&str>,
    gtf: Option<&str>,
    output_fastq: Option<&str>,
    output_dimers: Option<&str>,
    split_by_amplicon: bool,
) -> Result<AmpliconResult> {
    run_amplicons_with_callback(
        bam_path,
        primers_path,
        threads,
        mode,
        max_edit_dist,
        end_length,
        print_summary,
        primer_tolerance,
        min_qs,
        len_range,
        max_reads,
        duplex_only,
        reference,
        gtf,
        output_fastq,
        output_dimers,
        split_by_amplicon,
        |_| {},
    )
}

pub fn run_amplicons_with_callback<F>(
    bam_path: &str,
    primers_path: &str,
    threads: usize,
    mode: MatchMode,
    max_edit_dist: usize,
    end_length: usize,
    print_summary: bool,
    primer_tolerance: i64,
    min_qs: f32,
    len_range: (usize, usize),
    max_reads: usize,
    duplex_only: bool,
    reference: Option<&str>,
    gtf: Option<&str>,
    output_fastq: Option<&str>,
    output_dimers: Option<&str>,
    split_by_amplicon: bool,
    mut progress_callback: F,
) -> Result<AmpliconResult>
where
    F: FnMut(&AmpliconResult),
{
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();

    let primers = load_primers(primers_path)?;

    // Build K-mer index for optimization
    let kmer_index = if matches!(mode, MatchMode::Semiglobal) {
        log::info!("Building k-mer index for {} primers...", primers.len());
        Some(KmerIndex::build(primers.clone(), 8))
    } else {
        None
    };

    if is_fastq_path(bam_path) {
        if duplex_only {
            return Err(anyhow::anyhow!(
                "duplex-only filtering is only supported for BAM input (dx tag)"
            ));
        }
        if matches!(mode, MatchMode::Coords) {
            log::warn!("Coordinate mode requested for FASTQ; falling back to semiglobal mode");
        }
        log::info!("Reading FASTQ: {}", bam_path);

        let file = File::open(bam_path).with_context(|| format!("Failed to open FASTQ file: {}", bam_path))?;
        let mut reader: Box<dyn BufRead> = if bam_path.to_ascii_lowercase().ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut main_fq_writer = if let Some(ref path) = output_fastq {
            if !split_by_amplicon {
                Some(create_fastq_writer(path)?)
            } else {
                None
            }
        } else {
            None
        };

        let mut dimers_fq_writer = if let Some(ref path) = output_dimers {
            Some(create_fastq_writer(path)?)
        } else {
            None
        };

        let mut writers: HashMap<String, Box<dyn std::io::Write>> = HashMap::new();

        let chunk_size = 500;
        let mut amplicons: HashMap<String, AmpliconStats> = HashMap::new();
        let mut primer_counts: HashMap<String, usize> = HashMap::new();
        let mut chimera_count = 0usize;
        let mut unmatched_count = 0usize;
        let mut dist_lengths = Vec::<f64>::new();
        let mut dist_qs = Vec::<f64>::new();
        let mut total_reads = 0;

        loop {
            let mut chunk = Vec::new();
            let mut line = String::new();
            let mut reached_eof = false;

            for _ in 0..chunk_size {
                if max_reads > 0 && total_reads >= max_reads {
                    break;
                }

                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    reached_eof = true;
                    break;
                }
                let h = line.clone();

                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    reached_eof = true;
                    break;
                }
                let seq = line.clone();

                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    reached_eof = true;
                    break;
                }

                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    reached_eof = true;
                    break;
                }
                let qual = line.clone();

                let seq_trim = trim_line_ending(&seq);
                let qual_trim = trim_line_ending(&qual);
                if seq_trim.len() != qual_trim.len() {
                    continue;
                }

                let len = seq_trim.len();
                if len < len_range.0 || len > len_range.1 {
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

                chunk.push(FastqRead {
                    id,
                    header: h,
                    seq: seq_trim.as_bytes().to_vec(),
                    qual: qual_trim.as_bytes().to_vec(),
                });
            }

            if chunk.is_empty() {
                break;
            }

            let chunk_results: Vec<MatchResult> = chunk
                .par_iter()
                .enumerate()
                .map(|(i, r)| {
                    match_fastq_read(
                        r,
                        i,
                        &primers,
                        kmer_index.as_ref(),
                        end_length,
                        max_edit_dist,
                        true,
                        false,
                    )
                })
                .collect();

            for mut res in chunk_results {
                total_reads += 1;
                if res.is_chimera {
                    chimera_count += 1;
                }
                if let Some(p) = &res.start_primer {
                    *primer_counts.entry(p.clone()).or_default() += 1;
                }
                if let Some(p) = &res.end_primer {
                    *primer_counts.entry(p.clone()).or_default() += 1;
                }

                if let Some(name) = &res.amplicon_name {
                    let r = &chunk[res.read_id_index];
                    let qs = if let Some(header_qs) = parse_qs_from_header(&r.header) {
                        header_qs
                    } else {
                        mean_qv_from_fastq_ascii(&r.qual) as f32
                    };

                    if qs < min_qs {
                        unmatched_count += 1;
                        continue;
                    }

                    res.quality = qs;
                    dist_lengths.push(res.length as f64);
                    dist_qs.push(res.quality as f64);

                    let is_valid = if let (Some(p1), Some(p2)) = (&res.start_primer, &res.end_primer) {
                        is_fwd_rev_pair(p1, p2)
                    } else {
                        false
                    };

                    if is_valid {
                        let stats = amplicons.entry(name.clone()).or_default();
                        stats.count += 1;
                        stats.lengths.push(res.length);
                        stats.qualities.push(res.quality);
                        stats.read_ids.push(r.id.clone());

                        if let Some(ref out_path) = output_fastq {
                            let writer = if split_by_amplicon {
                                let path = get_split_path(out_path, name);
                                writers.entry(name.clone()).or_insert_with(|| {
                                    create_fastq_writer(&path).expect("Failed to create split FASTQ writer")
                                })
                            } else {
                                main_fq_writer.as_mut().unwrap()
                            };
                            writer.write_all(b"@")?;
                            writer.write_all(r.id.as_bytes())?;
                            writer.write_all(b"\n")?;
                            writer.write_all(&r.seq)?;
                            writer.write_all(b"\n+\n")?;
                            writer.write_all(&r.qual)?;
                            writer.write_all(b"\n")?;
                        }
                    } else {
                        unmatched_count += 1;
                        if let Some(ref mut writer) = dimers_fq_writer {
                            writer.write_all(b"@")?;
                            writer.write_all(r.id.as_bytes())?;
                            writer.write_all(b"\n")?;
                            writer.write_all(&r.seq)?;
                            writer.write_all(b"\n+\n")?;
                            writer.write_all(&r.qual)?;
                            writer.write_all(b"\n")?;
                        }
                    }
                } else {
                    unmatched_count += 1;
                }
            }

            let mut current_amplicons = amplicons.clone();
            for stats in current_amplicons.values_mut() {
                stats.finalize();
            }
            let intermediate = AmpliconResult {
                amplicons: current_amplicons,
                chimera_count,
                unmatched_count,
                total_reads,
                rescued_count: 0,
                distributions: build_distributions_from_values(&dist_lengths, &dist_qs),
            };
            progress_callback(&intermediate);

            if reached_eof || (max_reads > 0 && total_reads >= max_reads) {
                break;
            }
        }

        for stats in amplicons.values_mut() {
            stats.finalize();
        }

        let result = AmpliconResult {
            amplicons,
            chimera_count,
            unmatched_count,
            total_reads,
            rescued_count: 0,
            distributions: build_distributions_from_values(&dist_lengths, &dist_qs),
        };

        if print_summary {
            eprintln!("\n=== nanoparse Summary ===");
            eprintln!("Total reads:    {:>8}", result.total_reads);
            let total_matched = result.total_reads.saturating_sub(result.unmatched_count);
            let denom = result.total_reads.max(1) as f64;
            eprintln!(
                "Matched:        {:>8} ({:.1}%)",
                total_matched,
                total_matched as f64 / denom * 100.0
            );
            eprintln!(
                "Unmatched:      {:>8} ({:.1}%)",
                result.unmatched_count,
                result.unmatched_count as f64 / denom * 100.0
            );
            eprintln!("Amplicon types: {:>8}", result.amplicons.len());
        }

        return Ok(result);
    }

    log::info!("Reading BAM: {}", bam_path);
    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM file")?;
    bam.set_threads(threads).ok();

    if let Some(reference) = reference {
        log::info!("Reference provided: {}", reference);
    }
    if let Some(gtf) = gtf {
        log::info!("Annotation provided: {}", gtf);
    }

    let mut total_seen = 0usize;
    let mut total_processed = 0usize;
    let mut amplicons: HashMap<String, AmpliconStats> = HashMap::new();
    let mut primer_counts: HashMap<String, usize> = HashMap::new();
    let mut chimera_count = 0usize;
    let mut unassigned_results: Vec<PendingUnassigned> = Vec::new();
    let mut amplicon_coords: HashMap<String, (String, Vec<i64>, Vec<i64>)> = HashMap::new();
    let mut dist_lengths = Vec::<f64>::new();
    let mut dist_qs = Vec::<f64>::new();

    for result in bam.records() {
        let record = result?;
        total_seen += 1;

        if max_reads > 0 && total_processed >= max_reads {
            break;
        }
        let seq_len = record.seq_len();
        if seq_len < len_range.0 || seq_len > len_range.1 {
            continue;
        }
        let qs = qv_from_record(&record);
        if qs < min_qs {
            continue;
        }
        if duplex_only {
            let is_duplex = match record.aux(b"dx") {
                Ok(Aux::I8(v)) => v == 1,
                Ok(Aux::U8(v)) => v == 1,
                Ok(Aux::I16(v)) => v == 1,
                Ok(Aux::U16(v)) => v == 1,
                Ok(Aux::I32(v)) => v == 1,
                Ok(Aux::U32(v)) => v == 1,
                _ => false,
            };
            if !is_duplex {
                continue;
            }
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let res = if matches!(mode, MatchMode::Coords) {
            match_read_coords(&record, 0)
        } else {
            match_read_sequence(
                &record,
                0,
                &primers,
                kmer_index.as_ref(),
                end_length,
                max_edit_dist,
                true,
                false,
                primer_tolerance,
            )
        };

        total_processed += 1;
        dist_lengths.push(res.length as f64);
        dist_qs.push(res.quality as f64);

        if res.is_chimera {
            chimera_count += 1;
        }
        if let Some(p) = &res.start_primer {
            *primer_counts.entry(p.clone()).or_default() += 1;
        }
        if let Some(p) = &res.end_primer {
            *primer_counts.entry(p.clone()).or_default() += 1;
        }

        if let Some(name) = &res.amplicon_name {
            let stats = amplicons.entry(name.clone()).or_default();
            stats.count += 1;
            stats.lengths.push(res.length);
            stats.qualities.push(res.quality);
            stats.read_ids.push(read_id);

            if let (Some(c), Some(s), Some(e)) = (&res.chrom, res.start, res.end) {
                if stats.chrom.is_none() {
                    stats.chrom = Some(c.clone());
                }
                let entry = amplicon_coords
                    .entry(name.clone())
                    .or_insert_with(|| (c.clone(), Vec::new(), Vec::new()));
                entry.1.push(s);
                entry.2.push(e);
            }
        } else {
            unassigned_results.push(PendingUnassigned { id: read_id, res });
        }

        if total_seen % 500 == 0 {
            let mut current_amplicons = amplicons.clone();
            for stats in current_amplicons.values_mut() {
                stats.finalize();
            }
            let intermediate = AmpliconResult {
                amplicons: current_amplicons,
                chimera_count,
                unmatched_count: unassigned_results.len(),
                total_reads: total_processed,
                rescued_count: 0,
                distributions: build_distributions_from_values(&dist_lengths, &dist_qs),
            };
            progress_callback(&intermediate);
        }
    }

    log::info!(
        "Processed {} reads after filters (seen={}, mode={}, tol={}, dist={}, end={}, min_qs={}, min_len={}, duplex_only={})...",
        total_processed,
        total_seen,
        if matches!(mode, MatchMode::Coords) { "coords" } else { "semiglobal" },
        primer_tolerance,
        max_edit_dist,
        end_length,
        len_range.0,
        len_range.1,
        duplex_only
    );

    // Compute Consensus Regions (Median)
    let mut inferred_regions: HashMap<String, (String, i64, i64)> = HashMap::new();
    for (name, (chrom, mut starts, mut ends)) in amplicon_coords {
        if starts.is_empty() {
            continue;
        }
        starts.sort();
        ends.sort();
        let med_start = starts[starts.len() / 2];
        let med_end = ends[ends.len() / 2];
        inferred_regions.insert(name, (chrom, med_start, med_end));
    }

    // Rescue Step (Pass 2)
    let mut rescued_count = 0;
    let mut unmatched_count = 0;

    // User requested "Start with setting primer tolerance to 200 for phase 2"
    // We strictly use 200 if tolerance is enabled, or larger if user specified larger.
    let rescue_tolerance = if primer_tolerance > 0 {
        std::cmp::max(primer_tolerance, 200)
    } else {
        0
    };

    if rescue_tolerance > 0 && !unassigned_results.is_empty() && !amplicons.is_empty() {
        for pending in unassigned_results {
            let res = pending.res;
            if let (Some(chrom), Some(start), Some(end)) = (&res.chrom, res.start, res.end) {
                let r_c_norm = chrom.replace("chr", "");

                let mut best_dist = rescue_tolerance * 2 + 1;
                let mut best_match = None;

                for (name, (t_chrom, t_start, t_end)) in &inferred_regions {
                    let t_c_norm = t_chrom.replace("chr", "");
                    if r_c_norm == t_c_norm {
                        let d_start = (start - t_start).abs();
                        let d_end = (end - t_end).abs();
                        if d_start <= rescue_tolerance && d_end <= rescue_tolerance {
                            let total = d_start + d_end;
                            if total < best_dist {
                                best_dist = total;
                                best_match = Some(name.clone());
                            }
                        }
                    }
                }

                if let Some(match_name) = best_match {
                    // Rescued!
                    rescued_count += 1;
                    let stats = amplicons.entry(match_name).or_default();
                    stats.count += 1;
                    stats.lengths.push(res.length);
                    stats.qualities.push(res.quality);
                    stats.read_ids.push(pending.id);
                } else {
                    unmatched_count += 1;
                }
            } else {
                unmatched_count += 1;
            }
        }
    } else {
        unmatched_count = unassigned_results.len();
    }


    if rescued_count > 0 {
        log::info!("Rescued {} reads using inferred coordinates", rescued_count);
    }

    // Finalize stats (compute medians, set coordinates)
    for (name, stats) in amplicons.iter_mut() {
        if let Some((chrom, start, end)) = inferred_regions.get(name) {
            stats.chrom = Some(chrom.clone());
            stats.start = Some(*start);
            stats.end = Some(*end);
        }
        stats.finalize();
    }

    let distributions = build_distributions_from_values(&dist_lengths, &dist_qs);

    let result = AmpliconResult {
        amplicons,
        chimera_count,
        unmatched_count,
        total_reads: total_processed,
        rescued_count,
        distributions,
    };

    if print_summary {
        eprintln!("\n=== nanoparse Summary ===");
        eprintln!("Total reads:    {:>8}", result.total_reads);

        let total_matched = result.total_reads - result.unmatched_count;
        let direct_matched = total_matched - rescued_count;

        let denom = result.total_reads.max(1) as f64;
        eprintln!(
            "Matched:        {:>8} ({:.1}%)",
            total_matched,
            total_matched as f64 / denom * 100.0
        );
        eprintln!(
            "  Direct:       {:>8} ({:.1}%)",
            direct_matched,
            direct_matched as f64 / denom * 100.0
        );
        eprintln!(
            "  Rescued:      {:>8} ({:.1}%)",
            rescued_count,
            rescued_count as f64 / denom * 100.0
        );

        eprintln!(
            "Unmatched:      {:>8} ({:.1}%)",
            result.unmatched_count,
            result.unmatched_count as f64 / denom * 100.0
        );
        eprintln!("Amplicon types: {:>8}", result.amplicons.len());

        // Primer counts table (include ALL primers)
        let mut sorted_primers: Vec<_> = primers
            .iter()
            .map(|p| (p.name.as_str(), *primer_counts.get(&p.name).unwrap_or(&0)))
            .collect();
        sorted_primers.sort_by(|a, b| b.1.cmp(&a.1));

        if !sorted_primers.is_empty() {
            eprintln!("\nPrimer Hit Counts (Any hit):");
            eprintln!("{:<40} | {:>8}", "Primer Name", "Hits");
            eprintln!("{}", "-".repeat(51));
            for (name, count) in sorted_primers {
                eprintln!("{:<40} | {:>8}", name, count);
            }
        }

        let mut sorted: Vec<_> = result.amplicons.iter().collect();
        sorted.sort_by(|a, b| b.1.count.cmp(&a.1.count));

        if !sorted.is_empty() {
            eprintln!("\nTop amplicons:");
            for (name, stats) in sorted.iter().take(10) {
                eprintln!(
                    "  {:40} {:>6} reads (median len: {}bp)",
                    name, stats.count, stats.median_length
                );
            }
        }
        eprintln!("=========================\n");
    }

    Ok(result)
}

pub fn run_amplicons_to_output(
    bam_path: &str,
    primers_path: &str,
    output_path: &str,
    threads: usize,
    mode: MatchMode,
    max_edit_dist: usize,
    end_length: usize,
    print_summary: bool,
    primer_tolerance: i64,
    min_qs: f32,
    len_range: (usize, usize),
    max_reads: usize,
    duplex_only: bool,
    reference: Option<&str>,
    gtf: Option<&str>,
    output_fastq: Option<&str>,
    output_dimers: Option<&str>,
    split_by_amplicon: bool,
) -> Result<()> {
    let result = run_amplicons(
        bam_path,
        primers_path,
        threads,
        mode,
        max_edit_dist,
        end_length,
        print_summary,
        primer_tolerance,
        min_qs,
        len_range,
        max_reads,
        duplex_only,
        reference,
        gtf,
        output_fastq,
        output_dimers,
        split_by_amplicon,
    )?;
    write_output(&result, output_path)?;
    Ok(())
}
