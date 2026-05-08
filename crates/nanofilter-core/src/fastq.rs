use anyhow::Result;
use chrono::DateTime;
use flate2::write::GzEncoder;
use flate2::Compression;
use nanoseq_core::filters::{PoreRange, TimeWindow};
use nanoseq_core::header::{parse_nanopore_header, NanoporeMetadata};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use crate::bam::{FilterSettings, UnalignedBamWriter};
use crate::barcode;
use crate::filter::calculate_phred_avg;

struct OwnedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

pub fn filter_fastq(
    input_path: &Path,
    output_path: &Path,
    qv_threshold: f64,
    min_len: usize,
    max_len: usize,
    output_bam: bool,
) -> Result<()> {
    let settings = FilterSettings {
        qv_threshold,
        min_len,
        max_len,
        channel_range: None,
        time_window: None,
    };
    filter_fastq_with_settings(input_path, output_path, &settings, output_bam)
}

pub fn filter_fastq_with_settings(
    input_path: &Path,
    output_path: &Path,
    settings: &FilterSettings,
    output_bam: bool,
) -> Result<()> {
    let mut reader = parse_fastx_file(input_path)?;

    let mut kept = 0;
    let mut total = 0;

    if output_bam {
        let mut bam_writer: Option<UnalignedBamWriter> = None;
        while let Some(record) = reader.next() {
            let record = record?;
            total += 1;
            let seq = record.seq();
            let len = seq.len();
            let metadata = parse_nanopore_header(record.id());

            if bam_writer.is_none() {
                bam_writer = Some(UnalignedBamWriter::new(output_path, &metadata, settings)?);
            }

            if len >= settings.min_len
                && len <= settings.max_len
                && matches_fastq_filters(
                    &metadata,
                    settings.channel_range,
                    settings.time_window.as_ref(),
                )
            {
                let qual = record
                    .qual()
                    .ok_or_else(|| anyhow::anyhow!("No quality scores"))?;
                let avg_qv = calculate_phred_avg(qual);
                if avg_qv >= settings.qv_threshold {
                    let id = record.id();
                    let read_name = id.split(|&b| b == b' ').next().unwrap_or(id);
                    bam_writer
                        .as_mut()
                        .unwrap()
                        .write_record(read_name, &seq, qual)?;
                    kept += 1;
                }
            }
        }
        eprintln!("FASTQ→BAM: Kept {} / {} reads", kept, total);
    } else {
        let output_file = File::create(output_path)?;
        let mut writer: Box<dyn Write> = if output_path.to_string_lossy().ends_with(".gz") {
            Box::new(GzEncoder::new(
                BufWriter::new(output_file),
                Compression::default(),
            ))
        } else {
            Box::new(BufWriter::new(output_file))
        };

        while let Some(record) = reader.next() {
            let record = record?;
            total += 1;
            let seq = record.seq();
            let len = seq.len();
            let metadata = parse_nanopore_header(record.id());
            if len >= settings.min_len
                && len <= settings.max_len
                && matches_fastq_filters(
                    &metadata,
                    settings.channel_range,
                    settings.time_window.as_ref(),
                )
            {
                let qual = record
                    .qual()
                    .ok_or_else(|| anyhow::anyhow!("No quality scores"))?;
                let avg_qv = calculate_phred_avg(qual);
                if avg_qv >= settings.qv_threshold {
                    record.write(&mut writer, None)?;
                    kept += 1;
                }
            }
        }
        eprintln!("FASTQ: Kept {} / {} reads", kept, total);
    }
    Ok(())
}

pub fn make_pe(
    input_path: &Path,
    output_prefix: &Path,
    len: usize,
    insert: usize,
    step: usize,
) -> Result<()> {
    if len == 0 {
        anyhow::bail!("--len must be greater than 0");
    }
    if step == 0 {
        anyhow::bail!("--step must be greater than 0");
    }
    if insert < len {
        anyhow::bail!("--insert must be at least as large as --len");
    }

    let (r1_path, r2_path) = make_pe_output_paths(output_prefix)?;
    let mut r1_writer = create_fastq_writer(&r1_path)?;
    let mut r2_writer = create_fastq_writer(&r2_path)?;
    let mut reader = parse_fastx_file(input_path)?;

    let mut total_reads = 0usize;
    let mut total_pairs = 0usize;
    let window = insert;
    let r2_offset = insert - len;

    while let Some(record) = reader.next() {
        let record = record?;
        total_reads += 1;

        let seq = record.seq();
        let qual = record
            .qual()
            .ok_or_else(|| anyhow::anyhow!("FASTQ record is missing quality scores"))?;

        if seq.len() < window || qual.len() < window {
            continue;
        }

        let read_id = record
            .id()
            .split(|&b| b == b' ')
            .next()
            .unwrap_or(record.id());
        let mut pair_idx = 0usize;

        let mut start = 0usize;
        while start + window <= seq.len() && start + window <= qual.len() {
            let r1_seq = &seq[start..start + len];
            let r1_qual = &qual[start..start + len];
            let r2_src = &seq[start + r2_offset..start + r2_offset + len];
            let r2_qual_src = &qual[start + r2_offset..start + r2_offset + len];
            let r2_seq = barcode::reverse_complement(r2_src);
            let mut r2_qual = r2_qual_src.to_vec();
            r2_qual.reverse();

            let pair_tag = format!("{}_{}", String::from_utf8_lossy(read_id), pair_idx);
            write_fastq_record(&mut r1_writer, pair_tag.as_bytes(), r1_seq, r1_qual, b"/1")?;
            write_fastq_record(
                &mut r2_writer,
                pair_tag.as_bytes(),
                &r2_seq,
                &r2_qual,
                b"/2",
            )?;

            total_pairs += 1;
            pair_idx += 1;
            start += step;
        }
    }

    eprintln!(
        "make-PE: processed {} reads and wrote {} read pairs",
        total_reads, total_pairs
    );
    Ok(())
}

fn matches_fastq_filters(
    metadata: &NanoporeMetadata,
    channel_range: Option<PoreRange>,
    time_window: Option<&TimeWindow>,
) -> bool {
    let channel_ok = match channel_range {
        Some(range) => metadata.channel.map(|c| range.contains(c)).unwrap_or(false),
        None => true,
    };
    if !channel_ok {
        return false;
    }

    match time_window {
        Some(window) => metadata
            .start_time
            .as_deref()
            .and_then(|t| DateTime::parse_from_rfc3339(t).ok())
            .map(|t| window.contains(&t))
            .unwrap_or(false),
        None => true,
    }
}

pub fn split_fastq_by_barcodes(
    input_path: &Path,
    barcodes_path: &Path,
    output_dir: &Path,
    max_mismatches: usize,
    search_dist: usize,
    threads: usize,
    auto_discover: usize,
    fast: bool,
) -> Result<()> {
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }

    let mut barcode_pairs = barcode::parse_barcodes(barcodes_path)?;

    if auto_discover > 0 {
        eprintln!(
            "Auto-discovery enabled: identifying top {} unassigned barcode pairs...",
            auto_discover
        );
        let discovered =
            identify_top_unassigned_pairs_fastq(input_path, &barcode_pairs, auto_discover, 10000)?;
        if !discovered.is_empty() {
            eprintln!(
                "  Added {} discovered pairs to splitting.",
                discovered.len()
            );
            barcode_pairs.extend(discovered);
        }
    }

    let matcher = if fast {
        eprintln!("Fast mode enabled: precomputing Hamming distance variants...");
        Some(Arc::new(barcode::BarcodeMatcher::new(
            barcode_pairs.clone(),
            max_mismatches,
        )))
    } else {
        None
    };

    let (in_tx, in_rx) = crossbeam_channel::bounded::<Vec<OwnedRecord>>(threads * 2);
    let (out_tx, out_rx) =
        crossbeam_channel::bounded::<Vec<(OwnedRecord, Option<String>)>>(threads * 2);

    let input_path_buf = input_path.to_path_buf();
    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = parse_fastx_file(&input_path_buf)?;
        let mut batch = Vec::with_capacity(2000);
        while let Some(record) = reader.next() {
            let record = record?;
            batch.push(OwnedRecord {
                id: record.id().to_vec(),
                seq: record.seq().to_vec(),
                qual: record.qual().map(|q| q.to_vec()),
            });
            if batch.len() >= 2000 {
                if in_tx.send(batch).is_err() {
                    return Ok(());
                }
                batch = Vec::with_capacity(2000);
            }
        }
        if !batch.is_empty() {
            let _ = in_tx.send(batch);
        }
        Ok(())
    });

    let worker_barcode_pairs = barcode_pairs.clone();
    let matcher_clone = matcher.clone();
    let worker_handle = std::thread::spawn(move || {
        in_rx.into_iter().for_each(|batch| {
            let processed: Vec<(OwnedRecord, Option<String>)> = batch
                .into_par_iter()
                .map(|record| {
                    let seq_len = record.seq.len();
                    if seq_len < 20 {
                        return (record, None);
                    }

                    if let Some(m) = &matcher_clone {
                        let suffix_idx = if seq_len > search_dist * 2 {
                            search_dist
                        } else {
                            seq_len / 2
                        };
                        if let Some(sample) =
                            barcode::match_fast_with_anchors(&record.seq, suffix_idx, m)
                        {
                            return (record, Some(sample));
                        }
                    } else {
                        let s_len = std::cmp::min(search_dist, seq_len);
                        let e_len = std::cmp::min(search_dist, seq_len);
                        let start_seq = &record.seq[..s_len];
                        let end_seq = &record.seq[seq_len - e_len..];

                        for pair in &worker_barcode_pairs {
                            if barcode::match_regions(start_seq, end_seq, pair, max_mismatches) {
                                return (record, Some(pair.sample.clone()));
                            }
                        }
                    }
                    (record, None)
                })
                .collect();
            let _ = out_tx.send(processed);
        });
    });

    let writer_out_dir = output_dir.to_path_buf();
    let is_gz = input_path.to_string_lossy().ends_with(".gz")
        || input_path.to_string_lossy().ends_with(".fq.gz");
    let writer_handle = std::thread::spawn(move || -> Result<HashMap<String, usize>> {
        let mut writers: HashMap<String, Box<dyn Write>> = HashMap::new();
        let mut unassigned_writer = create_fastq_writer(&writer_out_dir.join(if is_gz {
            "unassigned.fastq.gz"
        } else {
            "unassigned.fastq"
        }))?;
        let mut counts = HashMap::new();
        let mut total = 0;

        for batch in out_rx {
            for (record, matched_sample) in batch {
                total += 1;
                let writer = if let Some(sample) = matched_sample {
                    *counts.entry(sample.clone()).or_insert(0) += 1;
                    if let Some(w) = writers.get_mut(&sample) {
                        w
                    } else {
                        let path = writer_out_dir.join(format!(
                            "{}.fastq{}",
                            sample,
                            if is_gz { ".gz" } else { "" }
                        ));
                        writers.insert(sample.clone(), create_fastq_writer(&path)?);
                        writers.get_mut(&sample).unwrap()
                    }
                } else {
                    &mut unassigned_writer
                };

                writer.write_all(b"@")?;
                writer.write_all(&record.id)?;
                writer.write_all(b"\n")?;
                writer.write_all(&record.seq)?;
                writer.write_all(b"\n+\n")?;
                if let Some(q) = &record.qual {
                    writer.write_all(q)?;
                }
                writer.write_all(b"\n")?;

                if total % 100000 == 0 {
                    eprint!("\rProcessed {} reads...", total);
                }
            }
        }
        eprintln!("\rProcessed {} reads. Done.", total);
        Ok(counts)
    });

    reader_handle.join().unwrap()?;
    worker_handle.join().unwrap();
    let counts = writer_handle.join().unwrap()?;

    eprintln!("FASTQ Split Results:");
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));
    for (sample, count) in sorted {
        eprintln!("  {}: {}", sample, count);
    }
    Ok(())
}

fn create_fastq_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path)?;
    if path.to_string_lossy().ends_with(".gz") {
        Ok(Box::new(GzEncoder::new(
            BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn make_pe_output_paths(output_prefix: &Path) -> Result<(std::path::PathBuf, std::path::PathBuf)> {
    let prefix = output_prefix.to_string_lossy();
    let (base, suffix) = if let Some(stripped) = prefix.strip_suffix(".fastq.gz") {
        (stripped.to_string(), ".fastq.gz")
    } else if let Some(stripped) = prefix.strip_suffix(".fq.gz") {
        (stripped.to_string(), ".fastq.gz")
    } else if let Some(stripped) = prefix.strip_suffix(".fastq") {
        (stripped.to_string(), ".fastq")
    } else if let Some(stripped) = prefix.strip_suffix(".fq") {
        (stripped.to_string(), ".fastq")
    } else {
        (prefix.to_string(), ".fastq")
    };

    Ok((
        std::path::PathBuf::from(format!("{}_R1{}", base, suffix)),
        std::path::PathBuf::from(format!("{}_R2{}", base, suffix)),
    ))
}

fn write_fastq_record(
    writer: &mut Box<dyn Write>,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    suffix: &[u8],
) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(suffix)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn discover_barcodes_fastq(
    input_path: &Path,
    barcodes_path: Option<&Path>,
    sample_size: usize,
    _threads: usize,
) -> Result<()> {
    let mut reader = parse_fastx_file(input_path)?;

    let user_pairs = if let Some(p) = barcodes_path {
        Some(barcode::parse_barcodes(p)?)
    } else {
        None
    };

    let mut pocket_counts: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut user_barcode_matches: HashMap<String, usize> = HashMap::new();

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

    while let Some(record) = reader.next() {
        let record = record?;
        total += 1;
        if total > sample_size {
            break;
        }

        let seq = record.seq();
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

    eprintln!("FASTQ Discovery Results (sampled {} reads):", total);
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

fn identify_top_unassigned_pairs_fastq(
    input_path: &Path,
    existing_pairs: &[barcode::BarcodePair],
    n: usize,
    sample_size: usize,
) -> Result<Vec<barcode::BarcodePair>> {
    let mut reader = parse_fastx_file(input_path)?;
    let mut pair_counts: HashMap<(Vec<u8>, Vec<u8>), usize> = HashMap::new();
    let anchors: [(&str, &[u8], i32, usize); 2] = [
        ("BC2-Anchor", b"AATGATACGGCGACCACCGAGATCTACAC", 29, 10),
        ("BC1-Anchor", b"ATCTCGTATGCCGTCTTCTGCTTG", -10, 10),
    ];

    let mut total = 0;
    while let Some(record) = reader.next() {
        let record = record?;
        total += 1;
        if total > sample_size {
            break;
        }

        let full_seq = record.seq();
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
                    if pocket.len() == 10 {
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
