//! Pore idle-time statistics for BAM and FASTQ inputs.

use anyhow::{Context, Result};
use chrono::DateTime;
use flate2::read::MultiGzDecoder;
use nanoseq_core::format::{is_fastq_path, trim_line_ending};
use nanoseq_core::header::parse_nanopore_header;
use rust_htslib::bam::{self, Read, record::Aux};
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::output::write_output;

#[derive(Debug, Clone)]
struct ReadTiming {
    ch: i64,
    mx: i64,
    start_s: f64,
    end_s: f64,
}

#[derive(Debug, Clone)]
struct SequencingSummaryEntry {
    channel: Option<i64>,
    mux: Option<i64>,
    start_s: Option<f64>,
    duration_s: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct PoreStatsResult {
    pub total_reads: usize,
    pub reads_with_timing: usize,
    pub reads_with_explicit_duration: usize,
    pub reads_with_estimated_duration: usize,
    pub sequencing_summary_rows: usize,
    pub sequencing_summary_matches: usize,
    pub estimated_speed_bps: Option<f64>,
    pub groups_with_idles: usize,
    pub speed_bps: f64,
    pub long_idle_threshold_s: f64,
    pub max_idle_s: f64,
    pub global_mean_idle_s: f64,
    pub global_median_idle_s: f64,
    pub long_idle_count: usize,
    pub long_idle_pct: f64,
    pub channel_stats: Vec<ChannelIdleStats>,
}

#[derive(Debug, Serialize)]
pub struct ChannelIdleStats {
    pub ch: i64,
    pub mx: i64,
    pub idle_count: usize,
    pub mean_idle_s: f64,
    pub median_idle_s: f64,
    pub long_idles: usize,
}

#[derive(Debug, Serialize)]
pub struct PoreStatsSummary {
    pub total_reads: usize,
    pub reads_with_timing: usize,
    pub reads_with_explicit_duration: usize,
    pub reads_with_estimated_duration: usize,
    pub sequencing_summary_rows: usize,
    pub sequencing_summary_matches: usize,
    pub estimated_speed_bps: Option<f64>,
    pub groups_with_idles: usize,
    pub speed_bps: f64,
    pub long_idle_threshold_s: f64,
    pub max_idle_s: f64,
    pub global_mean_idle_s: f64,
    pub global_median_idle_s: f64,
    pub long_idle_count: usize,
    pub long_idle_pct: f64,
    pub top_channels: Vec<ChannelIdleStats>,
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

fn parse_iso_timestamp(value: &str) -> Option<f64> {
    DateTime::parse_from_rfc3339(value)
        .ok()
        .map(|dt| dt.timestamp_millis() as f64 / 1000.0)
}

fn aux_to_start_seconds(aux: Aux<'_>) -> Option<f64> {
    match aux {
        Aux::Float(v) => Some(v as f64),
        Aux::Double(v) => Some(v),
        Aux::I8(v) => Some(f64::from(v)),
        Aux::U8(v) => Some(f64::from(v)),
        Aux::I16(v) => Some(f64::from(v)),
        Aux::U16(v) => Some(f64::from(v)),
        Aux::I32(v) => Some(v as f64),
        Aux::U32(v) => Some(v as f64),
        Aux::String(v) => parse_iso_timestamp(v).or_else(|| v.parse::<f64>().ok()),
        _ => None,
    }
}

fn split_summary_fields(line: &str) -> Vec<&str> {
    if line.contains('\t') {
        line.split('\t').collect()
    } else {
        line.split_whitespace().collect()
    }
}

fn parse_summary_time(value: &str) -> Option<f64> {
    parse_iso_timestamp(value).or_else(|| value.parse::<f64>().ok())
}

fn load_sequencing_summary(
    path: &str,
) -> Result<(
    HashMap<String, SequencingSummaryEntry>,
    Vec<SequencingSummaryEntry>,
    usize,
    Option<f64>,
)> {
    let file = File::open(path).context("Failed to open sequencing summary file")?;
    let mut reader: Box<dyn BufRead> = if path.to_ascii_lowercase().ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut header_line = String::new();
    loop {
        header_line.clear();
        if reader.read_line(&mut header_line)? == 0 {
            return Ok((HashMap::new(), Vec::new(), 0, None));
        }
        if !header_line.trim().is_empty() {
            break;
        }
    }

    let headers = split_summary_fields(header_line.trim_end());
    let mut columns = HashMap::<&str, usize>::new();
    for (idx, col) in headers.iter().enumerate() {
        columns.insert(col.trim(), idx);
    }

    let read_id_idx = columns.get("read_id").copied();
    let parent_read_id_idx = columns.get("parent_read_id").copied();
    let channel_idx = columns.get("channel").copied();
    let mux_idx = columns.get("mux").copied();
    let start_idx = columns.get("start_time").copied();
    let duration_idx = columns.get("duration").copied();
    let template_duration_idx = columns.get("template_duration").copied();
    let seq_len_idx = columns.get("sequence_length_template").copied();

    let mut entries = HashMap::new();
    let mut rows = Vec::new();
    let mut row_count = 0usize;
    let mut estimated_bases = 0f64;
    let mut estimated_duration = 0f64;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        row_count += 1;
        let fields = split_summary_fields(trimmed);
        let get = |idx: Option<usize>| idx.and_then(|i| fields.get(i).copied()).map(str::trim);

        let read_id = get(read_id_idx)
            .filter(|s| !s.is_empty() && *s != "-")
            .map(|s| s.to_string())
            .or_else(|| {
                get(parent_read_id_idx)
                    .filter(|s| !s.is_empty() && *s != "-")
                    .map(|s| s.to_string())
            });

        let Some(read_id) = read_id else {
            continue;
        };

        let channel = get(channel_idx).and_then(|v| v.parse::<i64>().ok());
        let mux = get(mux_idx).and_then(|v| v.parse::<i64>().ok());
        let start_s = get(start_idx).and_then(parse_summary_time);
        let duration_s = get(duration_idx)
            .and_then(parse_summary_time)
            .or_else(|| get(template_duration_idx).and_then(parse_summary_time));
        let seq_len = get(seq_len_idx).and_then(|v| v.parse::<f64>().ok());

        if let (Some(len), Some(dur)) = (seq_len, duration_s) {
            if dur > 0.0 {
                estimated_bases += len;
                estimated_duration += dur;
            }
        }

        let entry = SequencingSummaryEntry {
            channel,
            mux,
            start_s,
            duration_s,
        };

        rows.push(entry.clone());
        entries.insert(read_id.clone(), entry.clone());
        if let Some(parent_read_id) = get(parent_read_id_idx)
            .filter(|s| !s.is_empty() && *s != "-")
            .map(|s| s.to_string())
        {
            entries.insert(parent_read_id, entry);
        }
    }

    let estimated_speed_bps = if estimated_bases > 0.0 && estimated_duration > 0.0 {
        Some(estimated_bases / estimated_duration)
    } else {
        None
    };

    Ok((entries, rows, row_count, estimated_speed_bps))
}

fn resolve_read_timing(
    read_id: &str,
    input_ch: Option<i64>,
    input_mx: Option<i64>,
    input_start_s: Option<f64>,
    input_end_s: Option<f64>,
    seq_len: f64,
    summary_index: Option<&HashMap<String, SequencingSummaryEntry>>,
    speed_bps: f64,
) -> Option<(ReadTiming, bool, bool)> {
    if let Some(summary_index) = summary_index {
        let summary = summary_index.get(read_id)?;
        let ch = summary.channel.or(input_ch)?;
        let mx = summary.mux.or(input_mx).unwrap_or(0);
        let start_s = summary.start_s?;
        let duration_s = summary.duration_s?;
        let end_s = start_s + duration_s;

        return Some((
            ReadTiming {
                ch,
                mx,
                start_s,
                end_s,
            },
            true,
            true,
        ));
    }

    let ch = input_ch?;
    let mx = input_mx.unwrap_or(0);
    let start_s = input_start_s?;
    let (end_s, explicit) = if let Some(end_s) = input_end_s {
        (end_s, true)
    } else {
        (start_s + seq_len / speed_bps, false)
    };

    Some((
        ReadTiming {
            ch,
            mx,
            start_s,
            end_s,
        },
        explicit,
        false,
    ))
}

fn read_summary_only_timings(
    summary_rows: &[SequencingSummaryEntry],
    total_rows: usize,
) -> (usize, usize, usize, usize, usize, Vec<ReadTiming>) {
    let mut reads_with_timing = 0usize;
    let mut timings = Vec::new();

    for summary in summary_rows {
        let (Some(ch), Some(start_s), Some(duration_s)) =
            (summary.channel, summary.start_s, summary.duration_s)
        else {
            continue;
        };

        timings.push(ReadTiming {
            ch,
            mx: summary.mux.unwrap_or(0),
            start_s,
            end_s: start_s + duration_s,
        });
        reads_with_timing += 1;
    }

    (
        total_rows,
        reads_with_timing,
        reads_with_timing,
        0,
        reads_with_timing,
        timings,
    )
}

fn median(sorted_values: &[f64]) -> f64 {
    if sorted_values.is_empty() {
        return 0.0;
    }
    let mid = sorted_values.len() / 2;
    if sorted_values.len() % 2 == 0 {
        (sorted_values[mid - 1] + sorted_values[mid]) / 2.0
    } else {
        sorted_values[mid]
    }
}

fn summarize_timings(
    timings: Vec<ReadTiming>,
    max_idle_s: f64,
    long_idle_s: f64,
    speed_bps: f64,
    reads_with_timing: usize,
    reads_with_explicit_duration: usize,
    reads_with_estimated_duration: usize,
    sequencing_summary_rows: usize,
    sequencing_summary_matches: usize,
    estimated_speed_bps: Option<f64>,
) -> PoreStatsResult {
    let total_reads = timings.len();
    let mut grouped: HashMap<(i64, i64), Vec<ReadTiming>> = HashMap::new();
    for timing in timings {
        grouped.entry((timing.ch, timing.mx)).or_default().push(timing);
    }

    let mut channel_stats = Vec::new();
    let mut all_idles = Vec::new();

    for ((ch, mx), reads) in grouped {
        let mut reads = reads;
        reads.sort_by(|a, b| {
            a.start_s
                .partial_cmp(&b.start_s)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let idles: Vec<f64> = reads
            .windows(2)
            .map(|window| (window[1].start_s - window[0].end_s).max(0.0))
            .filter(|idle| idle.is_finite() && *idle < max_idle_s)
            .collect();

        if idles.is_empty() {
            continue;
        }

        let mut sorted_idles = idles.clone();
        sorted_idles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idle_count = sorted_idles.len();
        let mean_idle_s = sorted_idles.iter().sum::<f64>() / idle_count as f64;
        let median_idle_s = median(&sorted_idles);
        let long_idles = sorted_idles
            .iter()
            .filter(|idle| **idle > long_idle_s)
            .count();

        all_idles.extend(sorted_idles.iter().copied());
        channel_stats.push(ChannelIdleStats {
            ch,
            mx,
            idle_count,
            mean_idle_s,
            median_idle_s,
            long_idles,
        });
    }

    channel_stats.sort_by(|a, b| {
        b.long_idles
            .cmp(&a.long_idles)
            .then_with(|| a.ch.cmp(&b.ch))
            .then_with(|| a.mx.cmp(&b.mx))
    });

    all_idles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let long_idle_count = all_idles.iter().filter(|idle| **idle > long_idle_s).count();
    let long_idle_pct = if all_idles.is_empty() {
        0.0
    } else {
        long_idle_count as f64 / all_idles.len() as f64 * 100.0
    };
    let global_mean_idle_s = if all_idles.is_empty() {
        0.0
    } else {
        all_idles.iter().sum::<f64>() / all_idles.len() as f64
    };
    let global_median_idle_s = median(&all_idles);

    PoreStatsResult {
        total_reads,
        reads_with_timing,
        reads_with_explicit_duration,
        reads_with_estimated_duration,
        sequencing_summary_rows,
        sequencing_summary_matches,
        estimated_speed_bps,
        groups_with_idles: channel_stats.len(),
        speed_bps,
        long_idle_threshold_s: long_idle_s,
        max_idle_s,
        global_mean_idle_s,
        global_median_idle_s,
        long_idle_count,
        long_idle_pct,
        channel_stats,
    }
}

fn format_summary(result: &PoreStatsResult) -> String {
    let mut lines = Vec::new();
    lines.push("Pore idle-time summary".to_string());
    lines.push(format!(
        "reads: total={}, with_timing={}, explicit_duration={}, estimated_duration={}",
        result.total_reads,
        result.reads_with_timing,
        result.reads_with_explicit_duration,
        result.reads_with_estimated_duration
    ));
    if result.sequencing_summary_rows > 0 {
        lines.push(format!(
            "sequencing summary: rows={}, matched_reads={}",
            result.sequencing_summary_rows, result.sequencing_summary_matches
        ));
    }
    if let Some(estimated_speed_bps) = result.estimated_speed_bps {
        lines.push(format!(
            "estimated speed from summary: {:.1} bases/s",
            estimated_speed_bps
        ));
    }
    lines.push(format!(
        "idle: mean={:.3}s, median={:.3}s, long>{:.0}s={} ({:.2}%)",
        result.global_mean_idle_s,
        result.global_median_idle_s,
        result.long_idle_threshold_s,
        result.long_idle_count,
        result.long_idle_pct
    ));
    lines.push(format!(
        "groups: {} channel/mux groups with idle measurements",
        result.groups_with_idles
    ));
    lines.push(format!(
        "estimation: speed={} bases/s, max_idle={:.0}s",
        result.speed_bps, result.max_idle_s
    ));
    if !result.channel_stats.is_empty() {
        lines.push("top channels by long idles:".to_string());
        for s in result.channel_stats.iter().take(10) {
            lines.push(format!(
                "  ch={} mx={} idle_count={} mean={:.3}s median={:.3}s long={}",
                s.ch, s.mx, s.idle_count, s.mean_idle_s, s.median_idle_s, s.long_idles
            ));
        }
        if result.channel_stats.len() > 10 {
            lines.push(format!(
                "  ... and {} more groups",
                result.channel_stats.len() - 10
            ));
        }
    }
    lines.join("\n")
}

fn read_fastq_timings_from_reader<R: BufRead>(
    mut reader: R,
    summary_index: Option<&HashMap<String, SequencingSummaryEntry>>,
    speed_bps: f64,
) -> Result<(usize, usize, usize, usize, usize, Vec<ReadTiming>)> {
    let mut total_reads = 0usize;
    let mut reads_with_timing = 0usize;
    let mut reads_with_explicit_duration = 0usize;
    let mut reads_with_estimated_duration = 0usize;
    let mut summary_matches = 0usize;
    let mut timings = Vec::new();

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
        let header = trim_line_ending(&h);
        let meta = parse_nanopore_header(header.as_bytes());
        let seq_len = trim_line_ending(&seq).len() as f64;
        if seq_len == 0.0 {
            continue;
        }

        let Some((timing, explicit, used_summary)) = resolve_read_timing(
            &meta.read_id,
            meta.channel,
            Some(0),
            meta.start_time.as_deref().and_then(parse_iso_timestamp),
            None,
            seq_len,
            summary_index,
            speed_bps,
        ) else {
            continue;
        };

        reads_with_timing += 1;
        if explicit {
            reads_with_explicit_duration += 1;
        } else {
            reads_with_estimated_duration += 1;
        }
        if used_summary {
            summary_matches += 1;
        }
        timings.push(timing);
    }

    Ok((
        total_reads,
        reads_with_timing,
        reads_with_explicit_duration,
        reads_with_estimated_duration,
        summary_matches,
        timings,
    ))
}

fn read_bam_timings(
    path: &str,
    threads: usize,
    summary_index: Option<&HashMap<String, SequencingSummaryEntry>>,
    speed_bps: f64,
) -> Result<(usize, usize, usize, usize, usize, Vec<ReadTiming>)> {
    let mut bam = bam::Reader::from_path(path).context("Failed to open BAM file")?;
    bam.set_threads(threads).ok();

    let mut total_reads = 0usize;
    let mut reads_with_timing = 0usize;
    let mut reads_with_explicit_duration = 0usize;
    let mut reads_with_estimated_duration = 0usize;
    let mut summary_matches = 0usize;
    let mut timings = Vec::new();

    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        total_reads += 1;

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let ch = record.aux(b"ch").ok().and_then(aux_to_i64);
        let start_s = record
            .aux(b"st_ts")
            .ok()
            .and_then(aux_to_start_seconds)
            .or_else(|| record.aux(b"st").ok().and_then(aux_to_start_seconds));
        let mx = record.aux(b"mx").ok().and_then(aux_to_i64);
        let explicit_end = record
            .aux(b"et")
            .ok()
            .and_then(aux_to_start_seconds)
            .filter(|end_s| start_s.map(|s| *end_s >= s).unwrap_or(true));
        let explicit_duration = record
            .aux(b"du")
            .ok()
            .and_then(aux_to_start_seconds)
            .filter(|duration_s| *duration_s >= 0.0);
        let seq_len = record.seq_len() as f64;
        if seq_len == 0.0 {
            continue;
        }

        let input_end_s = explicit_end.or_else(|| explicit_duration.zip(start_s).map(|(d, s)| s + d));
        let input_start_s = start_s;
        let Some((timing, explicit, used_summary)) = resolve_read_timing(
            &read_id,
            ch,
            mx,
            input_start_s,
            input_end_s,
            seq_len,
            summary_index,
            speed_bps,
        ) else {
            continue;
        };

        reads_with_timing += 1;
        if explicit {
            reads_with_explicit_duration += 1;
        } else {
            reads_with_estimated_duration += 1;
        }
        if used_summary {
            summary_matches += 1;
        }
        timings.push(timing);
    }

    Ok((
        total_reads,
        reads_with_timing,
        reads_with_explicit_duration,
        reads_with_estimated_duration,
        summary_matches,
        timings,
    ))
}

pub fn run_pore_stats(
    input_path: Option<&str>,
    sequencing_summary_path: Option<&str>,
    output_path: Option<&str>,
    threads: usize,
    max_idle_s: f64,
    long_idle_s: f64,
    speed_bps: f64,
) -> Result<()> {
    let (sequencing_summary_index, sequencing_summary_rows_data, sequencing_summary_rows, estimated_speed_bps) =
        if let Some(summary_path) = sequencing_summary_path {
            log::info!("Reading sequencing summary: {}", summary_path);
            load_sequencing_summary(summary_path)?
        } else {
            (HashMap::new(), Vec::new(), 0, None)
        };
    let summary_ref = sequencing_summary_path.map(|_| &sequencing_summary_index);

    let (total_reads, reads_with_timing, reads_with_explicit_duration, reads_with_estimated_duration, sequencing_summary_matches, timings) =
        match input_path {
            Some(input_path) if is_fastq_path(input_path) => {
                log::info!("Reading FASTQ: {}", input_path);
                let file = File::open(input_path).context("Failed to open FASTQ file")?;
                if input_path.to_ascii_lowercase().ends_with(".gz") {
                    read_fastq_timings_from_reader(
                        BufReader::new(MultiGzDecoder::new(file)),
                        summary_ref,
                        speed_bps,
                    )?
                } else {
                    read_fastq_timings_from_reader(BufReader::new(file), summary_ref, speed_bps)?
                }
            }
            Some(input_path) => {
                log::info!("Reading BAM: {}", input_path);
                read_bam_timings(input_path, threads, summary_ref, speed_bps)?
            }
            None => {
                if summary_ref.is_none() {
                    anyhow::bail!("Provide either an input BAM/FASTQ file or --sequencing-summary");
                }
                read_summary_only_timings(&sequencing_summary_rows_data, sequencing_summary_rows)
            }
        };

    let mut result = summarize_timings(
        timings,
        max_idle_s,
        long_idle_s,
        speed_bps,
        reads_with_timing,
        reads_with_explicit_duration,
        reads_with_estimated_duration,
        sequencing_summary_rows,
        sequencing_summary_matches,
        estimated_speed_bps,
    );
    result.total_reads = total_reads;

    if let Some(output_path) = output_path {
        write_output(&result, output_path)?;
    }

    println!("{}", format_summary(&result));
    log::info!(
        "Done. {} / {} reads had timing metadata.",
        result.reads_with_timing,
        result.total_reads
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{ReadTiming, format_summary, median, summarize_timings};

    #[test]
    fn calculates_idle_stats_per_channel_and_mux() {
        let result = summarize_timings(
            vec![
                ReadTiming {
                    ch: 100,
                    mx: 1,
                    start_s: 10.0,
                    end_s: 12.0,
                },
                ReadTiming {
                    ch: 100,
                    mx: 1,
                    start_s: 25.0,
                    end_s: 28.0,
                },
                ReadTiming {
                    ch: 100,
                    mx: 1,
                    start_s: 40.0,
                    end_s: 41.0,
                },
                ReadTiming {
                    ch: 100,
                    mx: 2,
                    start_s: 5.0,
                    end_s: 10.0,
                },
                ReadTiming {
                    ch: 100,
                    mx: 2,
                    start_s: 80.0,
                    end_s: 81.0,
                },
            ],
            3600.0,
            60.0,
            400.0,
            5,
            2,
            3,
            0,
            0,
            None,
        );

        assert_eq!(result.reads_with_timing, 5);
        assert_eq!(result.reads_with_explicit_duration, 2);
        assert_eq!(result.reads_with_estimated_duration, 3);
        assert_eq!(result.groups_with_idles, 2);
        assert_eq!(result.long_idle_count, 1);
        assert!((result.global_mean_idle_s - 31.6666666667).abs() < 1e-6);
        assert!((result.global_median_idle_s - 13.0).abs() < 1e-6);
        assert_eq!(result.channel_stats[0].mx, 2);
        assert_eq!(result.channel_stats[0].long_idles, 1);
    }

    #[test]
    fn clamps_negative_idle_to_zero() {
        let result = summarize_timings(
            vec![
                ReadTiming {
                    ch: 7,
                    mx: 1,
                    start_s: 10.0,
                    end_s: 20.0,
                },
                ReadTiming {
                    ch: 7,
                    mx: 1,
                    start_s: 15.0,
                    end_s: 18.0,
                },
            ],
            3600.0,
            60.0,
            400.0,
            2,
            0,
            2,
            0,
            0,
            None,
        );

        assert_eq!(result.channel_stats[0].idle_count, 1);
        assert_eq!(result.channel_stats[0].mean_idle_s, 0.0);
    }

    #[test]
    fn median_handles_even_and_odd_counts() {
        assert_eq!(median(&[]), 0.0);
        assert_eq!(median(&[5.0]), 5.0);
        assert_eq!(median(&[1.0, 3.0, 9.0]), 3.0);
        assert_eq!(median(&[1.0, 3.0, 5.0, 7.0]), 4.0);
    }

    #[test]
    fn summary_mentions_idle_statistics() {
        let result = summarize_timings(
            vec![
                ReadTiming {
                    ch: 1,
                    mx: 1,
                    start_s: 10.0,
                    end_s: 11.0,
                },
                ReadTiming {
                    ch: 1,
                    mx: 1,
                    start_s: 20.0,
                    end_s: 21.0,
                },
            ],
            3600.0,
            60.0,
            400.0,
            2,
            0,
            2,
            0,
            0,
            None,
        );

        let summary = format_summary(&result);
        assert!(summary.contains("Pore idle-time summary"));
        assert!(summary.contains("mean="));
        assert!(summary.contains("estimated_duration"));
    }

    #[test]
    fn sequencing_summary_entries_mark_explicit_duration() {
        let mut summary_index = std::collections::HashMap::new();
        summary_index.insert(
            "read1".to_string(),
            super::SequencingSummaryEntry {
                channel: Some(1537),
                mux: Some(2),
                start_s: Some(963.2436),
                duration_s: Some(3.753),
            },
        );

        let resolved = super::resolve_read_timing(
            "read1",
            None,
            None,
            None,
            None,
            1500.0,
            Some(&summary_index),
            400.0,
        )
        .unwrap();

        assert!(resolved.1);
        assert!(resolved.2);
        assert_eq!(resolved.0.ch, 1537);
        assert_eq!(resolved.0.mx, 2);
        assert!((resolved.0.end_s - 966.9966).abs() < 1e-4);
    }
}
