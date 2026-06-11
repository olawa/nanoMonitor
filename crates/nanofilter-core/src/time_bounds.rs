use anyhow::{anyhow, Context, Result};
use chrono::{DateTime, FixedOffset};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use needletail::parse_fastx_file;
use nanoseq_core::header::parse_nanopore_header;
use noodles::bam;
use crate::bam::get_start_time_from_record;

fn parse_summary_timestamp(s: &str) -> Option<DateTime<FixedOffset>> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    // 1. Try RFC3339 (e.g. 2026-03-17T12:00:00Z)
    if let Ok(dt) = DateTime::parse_from_rfc3339(s) {
        return Some(dt);
    }
    // 2. Try YYYY-MM-DD HH:MM:SS (e.g., 2026-03-17 12:00:00)
    if let Ok(naive) = chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%d %H:%M:%S") {
        let offset = FixedOffset::east_opt(0).unwrap();
        return Some(DateTime::from_naive_utc_and_offset(naive, offset));
    }
    // 3. Try YYYY-MM-DD HH:MM:SS.3f (e.g., 2026-03-17 12:00:00.123)
    if let Ok(naive) = chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%d %H:%M:%S%.3f") {
        let offset = FixedOffset::east_opt(0).unwrap();
        return Some(DateTime::from_naive_utc_and_offset(naive, offset));
    }
    // 4. Try YYYY-MM-DD HH:MM:SS.6f
    if let Ok(naive) = chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%d %H:%M:%S%.6f") {
        let offset = FixedOffset::east_opt(0).unwrap();
        return Some(DateTime::from_naive_utc_and_offset(naive, offset));
    }
    // 5. Try YYYY-MM-DDTHH:MM:SS
    if let Ok(naive) = chrono::NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S") {
        let offset = FixedOffset::east_opt(0).unwrap();
        return Some(DateTime::from_naive_utc_and_offset(naive, offset));
    }
    // 6. Try parsing as a float (seconds since epoch)
    if let Ok(sec) = s.parse::<f64>() {
        let millis = (sec * 1000.0).round() as i64;
        if let Some(utc) = DateTime::from_timestamp_millis(millis) {
            let offset = FixedOffset::east_opt(0).unwrap();
            return Some(utc.with_timezone(&offset));
        }
    }
    None
}

pub fn scan_sequencing_summary_bounds(
    path: &Path,
) -> Result<(DateTime<FixedOffset>, DateTime<FixedOffset>)> {
    let file = File::open(path).context("Failed to open sequencing summary file")?;
    let path_str = path.to_string_lossy().to_lowercase();
    let mut reader: Box<dyn BufRead> = if path_str.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut header_line = String::new();
    loop {
        header_line.clear();
        if reader.read_line(&mut header_line)? == 0 {
            return Err(anyhow!("Empty sequencing summary file"));
        }
        if !header_line.trim().is_empty() {
            break;
        }
    }

    let headers: Vec<&str> = if header_line.contains('\t') {
        header_line.trim_end().split('\t').collect()
    } else {
        header_line.trim_end().split_whitespace().collect()
    };

    let start_time_idx = headers
        .iter()
        .position(|&h| h.trim() == "start_time")
        .ok_or_else(|| anyhow!("sequencing summary is missing 'start_time' column"))?;

    let mut min_time: Option<DateTime<FixedOffset>> = None;
    let mut max_time: Option<DateTime<FixedOffset>> = None;

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = if line.contains('\t') {
            trimmed.split('\t').collect()
        } else {
            trimmed.split_whitespace().collect()
        };

        if let Some(field) = fields.get(start_time_idx) {
            if let Some(dt) = parse_summary_timestamp(field) {
                min_time = Some(min_time.map_or(dt, |m| std::cmp::min(m, dt)));
                max_time = Some(max_time.map_or(dt, |m| std::cmp::max(m, dt)));
            }
        }
        line.clear();
    }

    match (min_time, max_time) {
        (Some(min), Some(max)) => Ok((min, max)),
        _ => Err(anyhow!("No valid start_time found in sequencing summary file")),
    }
}

pub fn scan_fastq_bounds(path: &Path) -> Result<(DateTime<FixedOffset>, DateTime<FixedOffset>)> {
    let mut reader = parse_fastx_file(path).map_err(|e| anyhow!("Failed to parse FASTQ: {e}"))?;
    let mut min_time: Option<DateTime<FixedOffset>> = None;
    let mut max_time: Option<DateTime<FixedOffset>> = None;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| anyhow!("Failed to read FASTQ record: {e}"))?;
        let meta = parse_nanopore_header(record.id());
        if let Some(ref st_str) = meta.start_time {
            if let Some(dt) = parse_summary_timestamp(st_str) {
                min_time = Some(min_time.map_or(dt, |m| std::cmp::min(m, dt)));
                max_time = Some(max_time.map_or(dt, |m| std::cmp::max(m, dt)));
            }
        }
    }

    match (min_time, max_time) {
        (Some(min), Some(max)) => Ok((min, max)),
        _ => Err(anyhow!("No valid start_time found in FASTQ headers")),
    }
}

pub fn scan_bam_bounds(
    path: &Path,
    threads: usize,
) -> Result<(DateTime<FixedOffset>, DateTime<FixedOffset>)> {
    let mut reader = bam::io::Reader::from(noodles::bgzf::MultithreadedReader::with_worker_count(
        std::num::NonZeroUsize::new(std::cmp::max(1, threads)).unwrap(),
        File::open(path)?,
    ));
    let _header = reader.read_header()?;

    let mut min_time: Option<DateTime<FixedOffset>> = None;
    let mut max_time: Option<DateTime<FixedOffset>> = None;
    let mut record = bam::Record::default();

    while reader.read_record(&mut record)? != 0 {
        if let Some(dt) = get_start_time_from_record(&record) {
            min_time = Some(min_time.map_or(dt, |m| std::cmp::min(m, dt)));
            max_time = Some(max_time.map_or(dt, |m| std::cmp::max(m, dt)));
        }
    }

    match (min_time, max_time) {
        (Some(min), Some(max)) => Ok((min, max)),
        _ => Err(anyhow!("No valid start_time found in BAM records")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use std::io::Write;
    use chrono::TimeZone;

    #[test]
    fn test_parse_summary_timestamp() {
        let utc = FixedOffset::east_opt(0).unwrap();
        
        // RFC3339
        let t1 = parse_summary_timestamp("2026-03-17T12:00:00Z").unwrap();
        assert_eq!(t1, utc.with_ymd_and_hms(2026, 3, 17, 12, 0, 0).unwrap());

        // YYYY-MM-DD HH:MM:SS
        let t2 = parse_summary_timestamp("2026-03-17 12:00:00").unwrap();
        assert_eq!(t2, utc.with_ymd_and_hms(2026, 3, 17, 12, 0, 0).unwrap());

        // Seconds since epoch as float
        let t3 = parse_summary_timestamp("1773748800.0").unwrap(); // 2026-03-17T12:00:00Z
        assert_eq!(t3, utc.with_ymd_and_hms(2026, 3, 17, 12, 0, 0).unwrap());
    }

    #[test]
    fn test_scan_sequencing_summary_bounds() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "read_id\tchannel\tstart_time\tduration\nread1\t101\t2026-03-17 12:00:00\t10.0\nread2\t102\t2026-03-17 12:30:00\t15.0"
        )
        .unwrap();

        let (min, max) = scan_sequencing_summary_bounds(file.path()).unwrap();
        let utc = FixedOffset::east_opt(0).unwrap();
        assert_eq!(min, utc.with_ymd_and_hms(2026, 3, 17, 12, 0, 0).unwrap());
        assert_eq!(max, utc.with_ymd_and_hms(2026, 3, 17, 12, 30, 0).unwrap());
    }
}
