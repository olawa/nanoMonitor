use anyhow::Result;
use chrono::{DateTime, FixedOffset};
use crossbeam_channel::{bounded, Receiver, Sender};
use nanoseq_core::filters::{PoreRange, TimeWindow};
use nanoseq_core::header::NanoporeMetadata;
use nanoseq_core::quality::{mean_qv_from_phred_scores, select_qv};
use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::QualityScores as _;
use rayon::prelude::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use crate::barcode;

#[derive(Debug, Clone)]
pub struct FilterSettings {
    pub qv_threshold: f64,
    pub min_len: usize,
    pub max_len: usize,
    pub channel_range: Option<PoreRange>,
    pub time_window: Option<TimeWindow>,
}

fn create_writer(path: &Path, header: &sam::Header) -> Result<bam::io::Writer<bgzf::Writer<File>>> {
    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(header)?;
    Ok(writer)
}

pub fn discover_barcodes(
    input_path: &Path,
    barcodes_path: Option<&Path>,
    sample_size: usize,
    threads: usize,
) -> Result<()> {
    let mut reader = bam::io::Reader::from(noodles::bgzf::MultithreadedReader::with_worker_count(
        std::num::NonZeroUsize::new(std::cmp::max(1, threads)).unwrap(),
        File::open(input_path)?,
    ));
    let _header = reader.read_header()?;

    let mut pocket_counts: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut user_barcode_matches: HashMap<String, usize> = HashMap::new();

    let user_pairs = if let Some(p) = barcodes_path {
        Some(barcode::parse_barcodes(p)?)
    } else {
        None
    };

    let anchors: [(&str, &[u8], i32, usize); 3] = [
        ("P5-Anchor", b"AATGATACGGCGACCACCGAGATCTACAC", 29, 10),
        ("P7-RC-Anchor", b"ATCTCGTATGCCGTCTTCTGCTTG", -10, 10),
        (
            "Read1-Adapter",
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            33,
            10,
        ),
    ];

    let mut total = 0;
    let mut anchor_hits = 0;

    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        total += 1;
        if total > sample_size {
            break;
        }

        let seq_buf = record.sequence();
        let seq: Vec<u8> = seq_buf.iter().collect();

        for (_name, anchor_seq, offset, len) in &anchors {
            let mut found_pos = None;
            let mut is_rc = false;

            if let Some(pos) = seq.windows(anchor_seq.len()).position(|w| w == *anchor_seq) {
                found_pos = Some(pos);
            } else {
                let rc = barcode::reverse_complement(*anchor_seq);
                if let Some(pos) = seq.windows(rc.len()).position(|w| w == rc) {
                    found_pos = Some(pos);
                    is_rc = true;
                }
            }

            if let Some(pos) = found_pos {
                anchor_hits += 1;
                let (p_start, p_end) = if is_rc {
                    if *offset >= 0 {
                        let start = pos.saturating_sub(*len);
                        (start, pos)
                    } else {
                        let start = pos + anchor_seq.len();
                        (start, start + *len)
                    }
                } else if *offset >= 0 {
                    (pos + *offset as usize, pos + *offset as usize + *len)
                } else {
                    let start = pos.saturating_sub(-*offset as usize);
                    (start, pos)
                };

                if p_end <= seq.len() {
                    let mut pocket = seq[p_start..p_end].to_vec();
                    if is_rc {
                        pocket = barcode::reverse_complement(&pocket);
                    }
                    if pocket.len() == 10 {
                        *pocket_counts.entry(pocket.clone()).or_insert(0) += 1;
                        if let Some(pairs) = &user_pairs {
                            for pair in pairs {
                                if barcode::edit_distance(&pocket, &pair.bc1, 1) <= 1 {
                                    *user_barcode_matches
                                        .entry(format!("{}_BC1", pair.sample))
                                        .or_insert(0) += 1;
                                }
                                if barcode::edit_distance(&pocket, &pair.bc2, 1) <= 1 {
                                    *user_barcode_matches
                                        .entry(format!("{}_BC2", pair.sample))
                                        .or_insert(0) += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Discovery Results (sampled {} reads):", total);
    eprintln!("  Anchor hits: {}", anchor_hits);
    if !user_barcode_matches.is_empty() {
        eprintln!("\nUser Barcode Matches in identified pockets:");
        let mut sorted: Vec<_> = user_barcode_matches.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));
        for (label, count) in sorted {
            eprintln!("  {}: {}", label, count);
        }
    }
    let mut sorted_pockets: Vec<_> = pocket_counts.into_iter().collect();
    sorted_pockets.sort_by(|a, b| b.1.cmp(&a.1));
    eprintln!("\nTop 10-mers found in anchor-flanking pockets:");
    for (seq, count) in sorted_pockets.iter().take(20) {
        eprintln!("  {}: {}", String::from_utf8_lossy(seq), count);
    }
    Ok(())
}

pub fn split_bam_by_barcodes(
    input_path: &Path,
    barcodes_path: &Path,
    output_dir: &Path,
    max_mismatches: usize,
    search_dist: usize,
    threads: usize,
    _auto_discover: usize,
    fast: bool,
) -> Result<()> {
    let pool_threads = std::cmp::max(1, threads);
    rayon::ThreadPoolBuilder::new()
        .num_threads(pool_threads)
        .build_global()
        .ok();

    let barcode_pairs = Arc::new(barcode::parse_barcodes(barcodes_path)?);
    let matcher = if fast {
        eprintln!("Fast mode enabled: precomputing Hamming distance variants...");
        Some(Arc::new(barcode::BarcodeMatcher::new(
            (*barcode_pairs).clone(),
            max_mismatches,
        )))
    } else {
        None
    };

    let (out_tx, out_rx): (
        Sender<(bam::Record, Option<String>)>,
        Receiver<(bam::Record, Option<String>)>,
    ) = bounded(2000);

    let writer_out_dir = output_dir.to_path_buf();
    let mut temp_reader = bam::io::Reader::new(File::open(input_path)?);
    let header = temp_reader.read_header()?;
    let writer_header = header.clone();

    let writer_handle = std::thread::spawn(move || -> Result<HashMap<String, usize>> {
        let mut writers = HashMap::new();
        let mut counts = HashMap::new();
        let mut unassigned_writer =
            create_writer(&writer_out_dir.join("unassigned.bam"), &writer_header)?;

        let mut total = 0;
        let mut last_log = std::time::Instant::now();

        while let Ok((record, matched_sample)) = out_rx.recv() {
            total += 1;

            if let Some(sample) = matched_sample {
                *counts.entry(sample.clone()).or_insert(0) += 1;
                let writer = writers.entry(sample.clone()).or_insert_with(|| {
                    let path = writer_out_dir.join(format!("{}.bam", sample));
                    create_writer(&path, &writer_header).expect("Failed to create writer")
                });
                writer.write_record(&writer_header, &record)?;
            } else {
                unassigned_writer.write_record(&writer_header, &record)?;
            }

            if last_log.elapsed().as_secs() >= 5 {
                eprintln!("Processed {} reads...", total);
                last_log = std::time::Instant::now();
            }
        }
        Ok(counts)
    });

    let (batch_tx, batch_rx) = bounded::<Vec<bam::Record>>(pool_threads * 4);
    let reader_path = input_path.to_path_buf();
    std::thread::spawn(move || {
        let file = File::open(reader_path).unwrap();
        let reader_threads = std::cmp::min(pool_threads / 2, 4).max(1);
        let mut reader = bam::io::Reader::from(bgzf::MultithreadedReader::with_worker_count(
            std::num::NonZeroUsize::new(reader_threads).unwrap(),
            file,
        ));
        let _ = reader.read_header().unwrap();

        let mut batch = Vec::with_capacity(1000);
        let mut record = bam::Record::default();

        while reader.read_record(&mut record).unwrap() != 0 {
            batch.push(record.clone());
            record = bam::Record::default();

            if batch.len() >= 1000 {
                let to_send = std::mem::replace(&mut batch, Vec::with_capacity(1000));
                if batch_tx.send(to_send).is_err() {
                    return;
                }
            }
        }
        if !batch.is_empty() {
            let _ = batch_tx.send(batch);
        }
    });

    let matcher_clone = matcher.clone();
    batch_rx
        .into_iter()
        .par_bridge()
        .for_each_with(out_tx, |tx, batch| {
            for record in batch {
                let matched = match_one_record(
                    &record,
                    &barcode_pairs,
                    matcher_clone.as_deref(),
                    search_dist,
                    max_mismatches,
                );
                let _ = tx.send((record, matched));
            }
        });

    let counts = writer_handle.join().unwrap()?;
    eprintln!("Done. Results:");
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));
    for (s, c) in sorted {
        eprintln!("  {}: {}", s, c);
    }
    Ok(())
}

thread_local! {
    static SCRATCH: RefCell<Vec<u8>> = RefCell::new(Vec::with_capacity(2048));
}

fn match_one_record(
    record: &bam::Record,
    pairs: &[barcode::BarcodePair],
    matcher: Option<&barcode::BarcodeMatcher>,
    search_dist: usize,
    max_mismatches: usize,
) -> Option<String> {
    let seq_buf = record.sequence();
    let seq_len = seq_buf.len();
    if seq_len < 20 {
        return None;
    }

    if let Some(m) = matcher {
        return SCRATCH.with(|scratch| {
            let mut s = scratch.borrow_mut();
            s.clear();

            let s_len = std::cmp::min(search_dist, seq_len);
            for b in seq_buf.iter().take(s_len) {
                s.push(u8::from(b));
            }

            let has_middle_gap = seq_len > search_dist * 2;
            let start_of_suffix = if has_middle_gap {
                let current_len = s.len();
                let e_len = search_dist;
                for b in seq_buf.iter().skip(seq_len - e_len) {
                    s.push(u8::from(b));
                }
                current_len
            } else {
                for b in seq_buf.iter().skip(s_len) {
                    s.push(u8::from(b));
                }
                s_len
            };

            match_fast_with_anchors(&s, start_of_suffix, m)
        });
    }

    let s_len = std::cmp::min(search_dist, seq_len);
    let start_seq: Vec<u8> = seq_buf.iter().take(s_len).map(u8::from).collect();
    let e_start = seq_len.saturating_sub(search_dist);
    let end_seq: Vec<u8> = seq_buf.iter().skip(e_start).map(u8::from).collect();

    for pair in pairs {
        if barcode::match_regions(&start_seq, &end_seq, pair, max_mismatches) {
            return Some(pair.sample.clone());
        }
    }
    None
}

fn match_fast_with_anchors(
    seq: &[u8],
    suffix_start: usize,
    matcher: &barcode::BarcodeMatcher,
) -> Option<String> {
    barcode::match_fast_with_anchors(seq, suffix_start, matcher)
}

pub fn identify_top_unassigned_pairs(
    input_path: &Path,
    existing_pairs: &[barcode::BarcodePair],
    n: usize,
    sample_size: usize,
    threads: usize,
) -> Result<Vec<barcode::BarcodePair>> {
    let mut reader = bam::io::Reader::from(noodles::bgzf::MultithreadedReader::with_worker_count(
        std::num::NonZeroUsize::new(std::cmp::max(1, threads)).unwrap(),
        File::open(input_path)?,
    ));
    let _header = reader.read_header()?;
    let mut pair_counts: HashMap<(Vec<u8>, Vec<u8>), usize> = HashMap::new();
    let anchors: [(&str, &[u8], i32, usize); 2] = [
        ("BC2-Anchor", b"AATGATACGGCGACCACCGAGATCTACAC", 29, 30),
        ("BC1-Anchor", b"ATCTCGTATGCCGTCTTCTGCTTG", -10, 30),
    ];

    let mut total = 0;
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        total += 1;
        if total > sample_size {
            break;
        }

        let seq_buf = record.sequence();
        let full_seq: Vec<u8> = seq_buf.iter().collect();
        let mut bc2_extracted = None;
        let mut bc1_extracted = None;

        for (name, anchor_seq, offset, len) in &anchors {
            let mut found_pos = None;
            let mut is_rc = false;
            if let Some(pos) = full_seq
                .windows(anchor_seq.len())
                .position(|w| w == *anchor_seq)
            {
                found_pos = Some(pos);
            } else {
                let rc = barcode::reverse_complement(*anchor_seq);
                if let Some(pos) = full_seq.windows(rc.len()).position(|w| w == rc) {
                    found_pos = Some(pos);
                    is_rc = true;
                }
            }

            if let Some(pos) = found_pos {
                let (p_start, p_end) = if is_rc {
                    if *offset >= 0 {
                        let start = pos.saturating_sub(*len);
                        (start, pos)
                    } else {
                        let start = pos + anchor_seq.len();
                        (start, start + *len)
                    }
                } else if *offset >= 0 {
                    (pos + *offset as usize, pos + *offset as usize + *len)
                } else {
                    let start = pos.saturating_sub(-*offset as usize);
                    (start, pos)
                };

                if p_end <= full_seq.len() {
                    let mut pocket = full_seq[p_start..p_end].to_vec();
                    if is_rc {
                        pocket = barcode::reverse_complement(&pocket);
                    }
                    if pocket.len() >= 10 {
                        if *name == "BC2-Anchor" {
                            bc2_extracted = Some(pocket);
                        } else {
                            bc1_extracted = Some(pocket);
                        }
                    }
                }
            }
        }
        if let (Some(b2), Some(b1)) = (bc2_extracted, bc1_extracted) {
            *pair_counts.entry((b1, b2)).or_insert(0) += 1;
        }
    }

    let mut sorted: Vec<_> = pair_counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));
    let mut found = Vec::new();
    for ((bc1, bc2), count) in sorted {
        let mut exists = false;
        for existing in existing_pairs {
            if (barcode::edit_distance(&bc1, &existing.bc1, 1) <= 1
                && barcode::edit_distance(&bc2, &existing.bc2, 1) <= 1)
                || (barcode::edit_distance(&bc1, &existing.bc2, 1) <= 1
                    && barcode::edit_distance(&bc2, &existing.bc1, 1) <= 1)
            {
                exists = true;
                break;
            }
        }
        if !exists && count > 5 {
            let bc1_s = String::from_utf8_lossy(&bc1).to_string();
            let bc2_s = String::from_utf8_lossy(&bc2).to_string();
            found.push(barcode::BarcodePair::new(
                format!("Auto_{}+{}", bc1_s, bc2_s),
                bc1_s,
                bc2_s,
            ));
            if found.len() >= n {
                break;
            }
        }
    }
    Ok(found)
}

pub fn filter_bam(
    input_path: &Path,
    output_path: &Path,
    qv_threshold: f64,
    min_len: usize,
    max_len: usize,
) -> Result<()> {
    let settings = FilterSettings {
        qv_threshold,
        min_len,
        max_len,
        channel_range: None,
        time_window: None,
    };
    filter_bam_with_settings(input_path, output_path, &settings)
}

pub fn filter_bam_with_settings(
    input_path: &Path,
    output_path: &Path,
    settings: &FilterSettings,
) -> Result<()> {
    filter_bam_with_settings_threads(input_path, output_path, settings, 1)
}

pub fn filter_bam_with_settings_threads(
    input_path: &Path,
    output_path: &Path,
    settings: &FilterSettings,
    threads: usize,
) -> Result<()> {
    let worker_count = std::num::NonZeroUsize::new(std::cmp::max(1, threads)).unwrap();
    let mut reader = bam::io::Reader::from(bgzf::MultithreadedReader::with_worker_count(
        worker_count,
        File::open(input_path)?,
    ));
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::from(bgzf::MultithreadedWriter::with_worker_count(
        worker_count,
        File::create(output_path)?,
    ));
    writer.write_header(&header)?;

    let mut kept = 0;
    let mut total = 0;
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        total += 1;
        let seq_len = record.sequence().len();
        if seq_len < settings.min_len || seq_len > settings.max_len {
            continue;
        }
        if !matches_channel_range(get_channel_from_record(&record), settings.channel_range) {
            continue;
        }
        if !matches_time_window_from_record(&record, settings.time_window.as_ref()) {
            continue;
        }
        let avg_qv = get_quality_score_lowlevel(&record)?;
        if avg_qv >= settings.qv_threshold {
            writer.write_record(&header, &record)?;
            kept += 1;
        }
    }
    eprintln!("BAM: Kept {} / {} reads", kept, total);
    Ok(())
}

pub fn bam_to_fastq(
    input_path: &Path,
    output_path: &Path,
    qv_threshold: f64,
    min_len: usize,
    max_len: usize,
) -> Result<()> {
    let settings = FilterSettings {
        qv_threshold,
        min_len,
        max_len,
        channel_range: None,
        time_window: None,
    };
    bam_to_fastq_with_settings(input_path, output_path, &settings)
}

pub fn bam_to_fastq_with_settings(
    input_path: &Path,
    output_path: &Path,
    settings: &FilterSettings,
) -> Result<()> {
    bam_to_fastq_with_settings_threads(input_path, output_path, settings, 1)
}

pub fn bam_to_fastq_with_settings_threads(
    input_path: &Path,
    output_path: &Path,
    settings: &FilterSettings,
    threads: usize,
) -> Result<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::BufWriter;

    let worker_count = std::num::NonZeroUsize::new(std::cmp::max(1, threads)).unwrap();
    let mut reader = bam::io::Reader::from(bgzf::MultithreadedReader::with_worker_count(
        worker_count,
        File::open(input_path)?,
    ));
    let _header = reader.read_header()?;

    let output_file = File::create(output_path)?;
    let mut writer: Box<dyn std::io::Write> = if output_path.to_string_lossy().ends_with(".gz") {
        Box::new(GzEncoder::new(
            BufWriter::new(output_file),
            Compression::default(),
        ))
    } else {
        Box::new(BufWriter::new(output_file))
    };

    let mut kept = 0;
    let mut total = 0;
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        total += 1;
        let seq_len = record.sequence().len();
        if seq_len < settings.min_len || seq_len > settings.max_len {
            continue;
        }
        if !matches_channel_range(get_channel_from_record(&record), settings.channel_range) {
            continue;
        }
        if !matches_time_window_from_record(&record, settings.time_window.as_ref()) {
            continue;
        }
        let avg_qv = get_quality_score_lowlevel(&record)?;
        if avg_qv >= settings.qv_threshold {
            use std::io::Write as _;
            writer.write_all(b"@")?;
            writer.write_all(
                record
                    .name()
                    .ok_or_else(|| anyhow::anyhow!("Missing read name"))?,
            )?;
            writer.write_all(b"\n")?;
            let seq = record.sequence();
            for b in seq.iter() {
                writer.write_all(&[b])?;
            }
            writer.write_all(b"\n+\n")?;
            let qual = record.quality_scores();
            for b in qual.iter() {
                writer.write_all(&[b + 33])?;
            }
            writer.write_all(b"\n")?;
            kept += 1;
        }
    }
    eprintln!("BAM->FASTQ: Kept {} / {} reads", kept, total);
    Ok(())
}

fn matches_channel_range(channel: Option<i64>, range: Option<PoreRange>) -> bool {
    match range {
        Some(range) => channel.map(|c| range.contains(c)).unwrap_or(false),
        None => true,
    }
}

fn get_channel_from_record(record: &bam::Record) -> Option<i64> {
    let tag = sam::alignment::record::data::field::Tag::new(b'c', b'm');
    match record.data().get(&tag) {
        Some(Ok(Value::Int8(v))) => Some(i64::from(v)),
        Some(Ok(Value::UInt8(v))) => Some(i64::from(v)),
        Some(Ok(Value::Int16(v))) => Some(i64::from(v)),
        Some(Ok(Value::UInt16(v))) => Some(i64::from(v)),
        Some(Ok(Value::Int32(v))) => Some(i64::from(v)),
        Some(Ok(Value::UInt32(v))) => Some(i64::from(v)),
        _ => None,
    }
}

pub fn get_start_time_from_record(record: &bam::Record) -> Option<DateTime<FixedOffset>> {
    let tag = sam::alignment::record::data::field::Tag::new(b's', b't');
    match record.data().get(&tag) {
        Some(Ok(value)) => parse_start_time_value(&value),
        _ => None,
    }
}

fn parse_start_time_value(value: &Value<'_>) -> Option<DateTime<FixedOffset>> {
    match value {
        Value::String(v) => {
            let raw = std::str::from_utf8(v.as_ref()).ok()?;
            DateTime::parse_from_rfc3339(raw).ok()
        }
        Value::Int8(v) => unix_seconds_to_datetime(*v as f64),
        Value::UInt8(v) => unix_seconds_to_datetime(*v as f64),
        Value::Int16(v) => unix_seconds_to_datetime(*v as f64),
        Value::UInt16(v) => unix_seconds_to_datetime(*v as f64),
        Value::Int32(v) => unix_seconds_to_datetime(*v as f64),
        Value::UInt32(v) => unix_seconds_to_datetime(*v as f64),
        Value::Float(v) => unix_seconds_to_datetime(*v as f64),
        _ => None,
    }
}

fn unix_seconds_to_datetime(seconds: f64) -> Option<DateTime<FixedOffset>> {
    let millis = (seconds * 1000.0).round() as i64;
    let utc = DateTime::from_timestamp_millis(millis)?;
    let offset = FixedOffset::east_opt(0)?;
    Some(utc.with_timezone(&offset))
}

fn matches_time_window_from_record(record: &bam::Record, window: Option<&TimeWindow>) -> bool {
    match window {
        Some(window) => get_start_time_from_record(record)
            .map(|t| window.contains(&t))
            .unwrap_or(false),
        None => true,
    }
}

pub fn get_quality_score_lowlevel(record: &bam::Record) -> Result<f64> {
    let tag = sam::alignment::record::data::field::Tag::new(b'Q', b'S');
    let explicit_qs = match record.data().get(&tag) {
        Some(Ok(Value::Float(qv))) => Some(qv as f64),
        Some(Ok(Value::Int32(qv))) => Some(qv as f64),
        _ => None,
    };

    let qual = record.quality_scores();
    if qual.is_empty() {
        return Ok(explicit_qs.unwrap_or(0.0));
    }

    let phred_scores: Vec<u8> = qual.iter().collect();
    Ok(select_qv(explicit_qs, &phred_scores))
}

pub fn create_unaligned_header(
    metadata: &NanoporeMetadata,
    filter_settings: &FilterSettings,
) -> Result<sam::Header> {
    let run_id = metadata.run_id.as_deref().unwrap_or("unknown");
    let header_text = format!(
        "@HD\tVN:1.6\tSO:unknown\n@RG\tID:{}\tPL:ONT\tSM:SAMPLE\tDS:run_id={}\n@PG\tID:nanofilter\tPN:nanofilter\tVN:0.1.0\tCL:qv_threshold={};min_len={};max_len={}",
        run_id, run_id, filter_settings.qv_threshold, filter_settings.min_len, filter_settings.max_len
    );
    header_text
        .parse()
        .map_err(|e| anyhow::anyhow!("Failed to parse header: {}", e))
}

pub struct UnalignedBamWriter {
    header: sam::Header,
    writer: bam::io::Writer<noodles::bgzf::MultithreadedWriter<File>>,
}

impl UnalignedBamWriter {
    pub fn new(
        output_path: &Path,
        metadata: &NanoporeMetadata,
        filter_settings: &FilterSettings,
    ) -> Result<Self> {
        let header = create_unaligned_header(metadata, filter_settings)?;
        let mut writer =
            bam::io::Writer::from(noodles::bgzf::MultithreadedWriter::with_worker_count(
                std::num::NonZeroUsize::new(1).unwrap(),
                File::create(output_path)?,
            ));
        writer.write_header(&header)?;
        Ok(Self { header, writer })
    }

    pub fn write_record(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        use sam::alignment::record_buf::RecordBuf;

        let mut record = RecordBuf::default();
        *record.name_mut() = Some(name.to_vec().into());
        *record.sequence_mut() = seq.to_vec().into();
        *record.quality_scores_mut() = qual.to_vec().into();

        let tag = sam::alignment::record::data::field::Tag::new(b'Q', b'S');
        record.data_mut().insert(
            tag,
            noodles::sam::alignment::record_buf::data::field::Value::Float(
                mean_qv_from_phred_scores(
                    &qual
                        .iter()
                        .map(|q| q.saturating_sub(33))
                        .collect::<Vec<_>>(),
                ) as f32,
            ),
        );

        sam::alignment::io::Write::write_alignment_record(&mut self.writer, &self.header, &record)
            .map_err(|e| e.into())
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_start_time_value, unix_seconds_to_datetime};
    use chrono::DateTime;
    use noodles::sam::alignment::record::data::field::Value;

    #[test]
    fn parses_rfc3339_start_time_value() {
        let value = Value::String("2026-03-17T10:30:00+00:00".into());
        assert!(parse_start_time_value(&value).is_some());
    }

    #[test]
    fn parses_numeric_start_time_as_unix_seconds() {
        let dt = unix_seconds_to_datetime(1_710_633_600.0).unwrap();
        let expected = DateTime::parse_from_rfc3339("2024-03-17T00:00:00+00:00").unwrap();
        assert_eq!(dt, expected);
    }
}
