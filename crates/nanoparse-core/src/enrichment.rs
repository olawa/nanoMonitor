//! BED enrichment and pore statistics for read-until experiments.

use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use std::cmp::{max, min};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::output::write_output;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Interval {
    start: i64,
    end: i64,
}

#[derive(Debug, Clone, Copy)]
struct CmRange {
    min: i64,
    max: i64,
}

#[derive(Debug, Default, Serialize)]
pub struct PoreStats {
    pub pore: i64,
    pub reads: u64,
    pub aligned_bases: u64,
    pub on_target_bases: u64,
    pub off_target_bases: u64,
}

#[derive(Debug, Serialize)]
pub struct EnrichmentResult {
    pub bam_path: String,
    pub bed_path: String,
    pub cm_range: Option<String>,
    pub total_records: u64,
    pub mapped_primary_records: u64,
    pub filtered_records: u64,
    pub reads_without_cm: u64,
    pub aligned_bases: u64,
    pub on_target_bases: u64,
    pub off_target_bases: u64,
    pub target_size_bases: u64,
    pub off_target_size_bases: u64,
    pub mean_target_coverage: f64,
    pub mean_off_target_coverage: f64,
    pub enrichment_ratio: Option<f64>,
    pub pores_with_reads: usize,
    pub pore_stats: Vec<PoreStats>,
}

fn parse_cm_range(raw: &str) -> Result<CmRange> {
    let trimmed = raw.trim();
    let (start, end) = trimmed
        .split_once('-')
        .ok_or_else(|| anyhow!("Invalid cm range '{trimmed}', expected START-END"))?;
    let min = start
        .trim()
        .parse::<i64>()
        .with_context(|| format!("Invalid cm start value in '{trimmed}'"))?;
    let max = end
        .trim()
        .parse::<i64>()
        .with_context(|| format!("Invalid cm end value in '{trimmed}'"))?;
    if min > max {
        return Err(anyhow!("Invalid cm range '{trimmed}': start > end"));
    }
    Ok(CmRange { min, max })
}

fn parse_bed(path: &str) -> Result<BTreeMap<String, Vec<Interval>>> {
    let file = File::open(path).with_context(|| format!("Failed to open BED file: {path}"))?;
    let reader = BufReader::new(file);
    let mut regions: BTreeMap<String, Vec<Interval>> = BTreeMap::new();

    for (idx, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed reading BED line {}", idx + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 3 {
            return Err(anyhow!(
                "Invalid BED line {}: expected at least 3 columns",
                idx + 1
            ));
        }

        let chrom = fields[0].to_string();
        let start = fields[1]
            .parse::<i64>()
            .with_context(|| format!("Invalid start at BED line {}", idx + 1))?;
        let end = fields[2]
            .parse::<i64>()
            .with_context(|| format!("Invalid end at BED line {}", idx + 1))?;
        if end <= start {
            continue;
        }

        regions
            .entry(chrom)
            .or_default()
            .push(Interval { start, end });
    }

    Ok(regions)
}

fn merge_intervals(intervals: &mut Vec<Interval>, chrom_len: i64) -> Vec<Interval> {
    intervals.sort_by_key(|iv| (iv.start, iv.end));
    let mut merged = Vec::<Interval>::new();

    for iv in intervals.iter() {
        let clipped = Interval {
            start: max(0, min(iv.start, chrom_len)),
            end: max(0, min(iv.end, chrom_len)),
        };
        if clipped.end <= clipped.start {
            continue;
        }

        if let Some(last) = merged.last_mut() {
            if clipped.start <= last.end {
                last.end = max(last.end, clipped.end);
            } else {
                merged.push(clipped);
            }
        } else {
            merged.push(clipped);
        }
    }

    merged
}

fn query_aligned_blocks(record: &bam::Record) -> Vec<Interval> {
    use rust_htslib::bam::record::Cigar::*;

    let mut blocks = Vec::<Interval>::new();
    let mut ref_pos = record.pos();

    for op in record.cigar().iter() {
        match *op {
            Match(n) | Equal(n) | Diff(n) => {
                let len = i64::from(n);
                blocks.push(Interval {
                    start: ref_pos,
                    end: ref_pos + len,
                });
                ref_pos += len;
            }
            Del(n) | RefSkip(n) => {
                ref_pos += i64::from(n);
            }
            Ins(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
    }

    blocks
}

fn overlap_bases(blocks: &[Interval], targets: &[Interval]) -> u64 {
    if blocks.is_empty() || targets.is_empty() {
        return 0;
    }

    let mut target_idx = 0usize;
    let mut total = 0u64;

    for block in blocks {
        while target_idx < targets.len() && targets[target_idx].end <= block.start {
            target_idx += 1;
        }

        let mut idx = target_idx;
        while idx < targets.len() && targets[idx].start < block.end {
            let start = max(block.start, targets[idx].start);
            let end = min(block.end, targets[idx].end);
            if end > start {
                total += (end - start) as u64;
            }
            if targets[idx].end >= block.end {
                break;
            }
            idx += 1;
        }
    }

    total
}

fn aux_to_i64(aux: Aux<'_>) -> Option<i64> {
    match aux {
        Aux::I8(v) => Some(i64::from(v)),
        Aux::U8(v) => Some(i64::from(v)),
        Aux::I16(v) => Some(i64::from(v)),
        Aux::U16(v) => Some(i64::from(v)),
        Aux::I32(v) => Some(i64::from(v)),
        Aux::U32(v) => Some(i64::from(v)),
        _ => None,
    }
}

pub fn run_enrichment(
    bam_path: &str,
    bed_path: &str,
    output_path: &str,
    threads: usize,
    cm_range_raw: Option<&str>,
) -> Result<EnrichmentResult> {
    let cm_range = cm_range_raw.map(parse_cm_range).transpose()?;

    let bed_regions = parse_bed(bed_path)?;
    let mut bam = bam::Reader::from_path(bam_path).context("Failed to open BAM file")?;
    bam.set_threads(threads).ok();

    let header = bam.header().to_owned();
    let mut merged_targets_by_tid: BTreeMap<i32, Vec<Interval>> = BTreeMap::new();
    let mut target_size_bases = 0u64;
    let mut total_reference_bases = 0u64;

    for tid in 0..header.target_count() {
        let tid_u32 = tid;
        let chrom_name = String::from_utf8_lossy(header.tid2name(tid_u32)).to_string();
        let chrom_len = header.target_len(tid_u32).unwrap_or(0) as i64;
        total_reference_bases += chrom_len as u64;

        if let Some(mut regions) = bed_regions.get(&chrom_name).cloned() {
            let merged = merge_intervals(&mut regions, chrom_len);
            let merged_size: u64 = merged.iter().map(|iv| (iv.end - iv.start) as u64).sum();
            if !merged.is_empty() {
                merged_targets_by_tid.insert(tid_u32 as i32, merged);
                target_size_bases += merged_size;
            }
        }
    }

    let off_target_size_bases = total_reference_bases.saturating_sub(target_size_bases);
    let mut pore_stats = BTreeMap::<i64, PoreStats>::new();
    let mut total_records = 0u64;
    let mut mapped_primary_records = 0u64;
    let mut filtered_records = 0u64;
    let mut reads_without_cm = 0u64;
    let mut aligned_bases = 0u64;
    let mut on_target_bases = 0u64;

    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        total_records += 1;

        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        mapped_primary_records += 1;

        let pore = match record.aux(b"cm").ok().and_then(aux_to_i64) {
            Some(value) => value,
            None => {
                reads_without_cm += 1;
                if cm_range.is_some() {
                    filtered_records += 1;
                    continue;
                }
                -1
            }
        };

        if let Some(range) = cm_range {
            if pore < range.min || pore > range.max {
                filtered_records += 1;
                continue;
            }
        }

        let blocks = query_aligned_blocks(&record);
        let record_aligned_bases: u64 = blocks.iter().map(|iv| (iv.end - iv.start) as u64).sum();
        if record_aligned_bases == 0 {
            continue;
        }

        let record_on_target_bases = merged_targets_by_tid
            .get(&record.tid())
            .map(|targets| overlap_bases(&blocks, targets))
            .unwrap_or(0);
        let record_off_target_bases = record_aligned_bases.saturating_sub(record_on_target_bases);

        aligned_bases += record_aligned_bases;
        on_target_bases += record_on_target_bases;

        if pore >= 0 {
            let entry = pore_stats.entry(pore).or_insert_with(|| PoreStats {
                pore,
                ..PoreStats::default()
            });
            entry.reads += 1;
            entry.aligned_bases += record_aligned_bases;
            entry.on_target_bases += record_on_target_bases;
            entry.off_target_bases += record_off_target_bases;
        }
    }

    let off_target_bases = aligned_bases.saturating_sub(on_target_bases);
    let mean_target_coverage = if target_size_bases > 0 {
        on_target_bases as f64 / target_size_bases as f64
    } else {
        0.0
    };
    let mean_off_target_coverage = if off_target_size_bases > 0 {
        off_target_bases as f64 / off_target_size_bases as f64
    } else {
        0.0
    };
    let enrichment_ratio = if mean_off_target_coverage > 0.0 {
        Some(mean_target_coverage / mean_off_target_coverage)
    } else {
        None
    };

    let result = EnrichmentResult {
        bam_path: bam_path.to_string(),
        bed_path: bed_path.to_string(),
        cm_range: cm_range_raw.map(str::to_string),
        total_records,
        mapped_primary_records,
        filtered_records,
        reads_without_cm,
        aligned_bases,
        on_target_bases,
        off_target_bases,
        target_size_bases,
        off_target_size_bases,
        mean_target_coverage,
        mean_off_target_coverage,
        enrichment_ratio,
        pores_with_reads: pore_stats.len(),
        pore_stats: pore_stats.into_values().collect(),
    };

    write_output(&result, output_path)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::{merge_intervals, overlap_bases, parse_cm_range, Interval};

    #[test]
    fn parses_cm_range() {
        let range = parse_cm_range("1-2000").unwrap();
        assert_eq!(range.min, 1);
        assert_eq!(range.max, 2000);
    }

    #[test]
    fn merges_intervals_and_clips() {
        let mut intervals = vec![
            Interval { start: -5, end: 10 },
            Interval { start: 5, end: 20 },
            Interval { start: 40, end: 60 },
        ];
        let merged = merge_intervals(&mut intervals, 50);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start, 0);
        assert_eq!(merged[0].end, 20);
        assert_eq!(merged[1].start, 40);
        assert_eq!(merged[1].end, 50);
    }

    #[test]
    fn computes_overlap_across_blocks() {
        let blocks = vec![
            Interval { start: 10, end: 20 },
            Interval { start: 30, end: 40 },
        ];
        let targets = vec![Interval { start: 15, end: 35 }];
        assert_eq!(overlap_bases(&blocks, &targets), 10);
    }
}
