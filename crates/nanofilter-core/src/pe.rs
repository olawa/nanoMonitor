use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::record::QualityScores as _;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::barcode;
use nanoseq_core::quality::phred_to_fastq_ascii;

#[derive(Debug, Clone)]
pub struct MakePeConfig {
    pub len: usize,
    pub insert: usize,
    pub step: usize,
    pub interleaved: bool,
    pub threads: usize,
}

impl MakePeConfig {
    pub fn new(
        len: usize,
        insert: Option<usize>,
        step: Option<usize>,
        interleaved: bool,
        threads: usize,
    ) -> Result<Self> {
        if len == 0 {
            bail!("--len must be greater than 0");
        }
        let insert = insert.unwrap_or(2 * len);
        if insert < len {
            bail!("--insert ({insert}) must be at least as large as --len ({len})");
        }
        let step = step.unwrap_or(insert);
        if step == 0 {
            bail!("--step / --shift must be greater than 0");
        }
        let threads = std::cmp::max(1, threads);
        Ok(Self {
            len,
            insert,
            step,
            interleaved,
            threads,
        })
    }
}

pub enum PeOutputTarget<'a> {
    Prefix(&'a Path),
    SingleFile(&'a Path),
    Stdout,
}

pub struct PeWriter {
    r1_writer: Box<dyn Write>,
    r2_writer: Option<Box<dyn Write>>,
    interleaved: bool,
}

impl PeWriter {
    pub fn new(output_target: PeOutputTarget, interleaved: bool) -> Result<Self> {
        match output_target {
            PeOutputTarget::Stdout => {
                let stdout = io::stdout();
                let writer = Box::new(BufWriter::new(stdout));
                Ok(Self {
                    r1_writer: writer,
                    r2_writer: None,
                    interleaved: true,
                })
            }
            PeOutputTarget::SingleFile(path) => {
                let writer = create_fastq_writer(path)?;
                Ok(Self {
                    r1_writer: writer,
                    r2_writer: None,
                    interleaved: true,
                })
            }
            PeOutputTarget::Prefix(prefix) => {
                if interleaved {
                    let path = resolve_interleaved_path(prefix);
                    let writer = create_fastq_writer(&path)?;
                    Ok(Self {
                        r1_writer: writer,
                        r2_writer: None,
                        interleaved: true,
                    })
                } else {
                    let (r1_path, r2_path) = make_pe_output_paths(prefix)?;
                    let r1_writer = create_fastq_writer(&r1_path)?;
                    let r2_writer = create_fastq_writer(&r2_path)?;
                    Ok(Self {
                        r1_writer,
                        r2_writer: Some(r2_writer),
                        interleaved: false,
                    })
                }
            }
        }
    }

    pub fn write_pair(
        &mut self,
        read_id: &[u8],
        frag_idx: usize,
        r1_seq: &[u8],
        r1_qual: &[u8],
        r2_seq: &[u8],
        r2_qual: &[u8],
    ) -> Result<()> {
        let clean_id = read_id
            .split(|&b| b == b' ' || b == b'\t')
            .next()
            .unwrap_or(read_id);
        let clean_id = clean_id.strip_prefix(b"@").unwrap_or(clean_id);
        let frag_tag = format!("{}_f{}", String::from_utf8_lossy(clean_id), frag_idx);
        let frag_bytes = frag_tag.as_bytes();

        if self.interleaved || self.r2_writer.is_none() {
            write_fastq_record(&mut self.r1_writer, frag_bytes, r1_seq, r1_qual, b"/1")?;
            write_fastq_record(&mut self.r1_writer, frag_bytes, r2_seq, r2_qual, b"/2")?;
        } else {
            write_fastq_record(&mut self.r1_writer, frag_bytes, r1_seq, r1_qual, b"/1")?;
            if let Some(r2) = self.r2_writer.as_mut() {
                write_fastq_record(r2, frag_bytes, r2_seq, r2_qual, b"/2")?;
            }
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.r1_writer.flush()?;
        if let Some(r2) = self.r2_writer.as_mut() {
            r2.flush()?;
        }
        Ok(())
    }
}

pub fn write_fastq_record(
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

pub struct PeSummary {
    pub total_reads: usize,
    pub total_pairs: usize,
    pub skipped_short: usize,
}

pub fn process_read_into_pe(
    read_id: &[u8],
    seq: &[u8],
    qual: &[u8],
    config: &MakePeConfig,
    writer: &mut PeWriter,
) -> Result<(usize, bool)> {
    let window = config.insert;
    let len = config.len;
    let step = config.step;
    let r2_offset = config.insert - len;

    if seq.len() < window || qual.len() < window {
        return Ok((0, false));
    }

    let mut start = 0usize;
    let mut frag_idx = 0usize;
    let mut pairs_written = 0usize;

    while start + window <= seq.len() && start + window <= qual.len() {
        let r1_seq = &seq[start..start + len];
        let r1_qual = &qual[start..start + len];

        let r2_raw_seq = &seq[start + r2_offset..start + r2_offset + len];
        let r2_raw_qual = &qual[start + r2_offset..start + r2_offset + len];

        let r2_seq = barcode::reverse_complement(r2_raw_seq);
        let mut r2_qual = r2_raw_qual.to_vec();
        r2_qual.reverse();

        writer.write_pair(read_id, frag_idx, r1_seq, r1_qual, &r2_seq, &r2_qual)?;
        pairs_written += 1;
        frag_idx += 1;
        start += step;
    }

    Ok((pairs_written, true))
}

pub fn make_pe(
    input_path: &Path,
    output_target: PeOutputTarget,
    config: &MakePeConfig,
) -> Result<PeSummary> {
    let mut writer = PeWriter::new(output_target, config.interleaved)?;
    let input_str = input_path.to_string_lossy().to_lowercase();
    let summary = if input_str.ends_with(".bam") {
        run_make_pe_bam(input_path, config, &mut writer)?
    } else {
        run_make_pe_fastq(input_path, config, &mut writer)?
    };
    Ok(summary)
}

fn run_make_pe_bam(
    input_path: &Path,
    config: &MakePeConfig,
    writer: &mut PeWriter,
) -> Result<PeSummary> {
    let worker_count = std::num::NonZeroUsize::new(std::cmp::max(1, config.threads)).unwrap();
    let mut reader = bam::io::Reader::from(bgzf::MultithreadedReader::with_worker_count(
        worker_count,
        File::open(input_path)
            .with_context(|| format!("Failed to open BAM file {:?}", input_path))?,
    ));
    let _header = reader.read_header()?;

    let mut total_reads = 0usize;
    let mut total_pairs = 0usize;
    let mut skipped_short = 0usize;

    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        total_reads += 1;
        let read_name: &[u8] = record.name().map(|n| n.as_ref()).unwrap_or(b"read");
        let seq: Vec<u8> = record.sequence().iter().collect();
        let qual: Vec<u8> = record
            .quality_scores()
            .iter()
            .map(phred_to_fastq_ascii)
            .collect();

        let (pairs, processed) = process_read_into_pe(read_name, &seq, &qual, config, writer)?;
        if processed {
            total_pairs += pairs;
        } else {
            skipped_short += 1;
        }
    }
    writer.flush()?;
    Ok(PeSummary {
        total_reads,
        total_pairs,
        skipped_short,
    })
}

fn run_make_pe_fastq(
    input_path: &Path,
    config: &MakePeConfig,
    writer: &mut PeWriter,
) -> Result<PeSummary> {
    let mut reader = needletail::parse_fastx_file(input_path)
        .with_context(|| format!("Failed to open FASTQ file {:?}", input_path))?;

    let mut total_reads = 0usize;
    let mut total_pairs = 0usize;
    let mut skipped_short = 0usize;

    while let Some(record) = reader.next() {
        let record = record?;
        total_reads += 1;

        let seq = record.seq();
        let qual = record
            .qual()
            .ok_or_else(|| anyhow::anyhow!("FASTQ record is missing quality scores"))?;

        let id = record.id();
        let (pairs, processed) = process_read_into_pe(id, &seq, qual, config, writer)?;
        if processed {
            total_pairs += pairs;
        } else {
            skipped_short += 1;
        }
    }
    writer.flush()?;
    Ok(PeSummary {
        total_reads,
        total_pairs,
        skipped_short,
    })
}

pub fn create_fastq_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create output file {:?}", path))?;
    if path.to_string_lossy().ends_with(".gz") {
        Ok(Box::new(GzEncoder::new(
            BufWriter::new(file),
            Compression::default(),
        )))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

pub fn make_pe_output_paths(output_prefix: &Path) -> Result<(PathBuf, PathBuf)> {
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
        PathBuf::from(format!("{}_R1{}", base, suffix)),
        PathBuf::from(format!("{}_R2{}", base, suffix)),
    ))
}

pub fn resolve_interleaved_path(output_prefix: &Path) -> PathBuf {
    let prefix = output_prefix.to_string_lossy();
    if prefix.ends_with(".fastq")
        || prefix.ends_with(".fastq.gz")
        || prefix.ends_with(".fq")
        || prefix.ends_with(".fq.gz")
    {
        output_prefix.to_path_buf()
    } else {
        PathBuf::from(format!("{}.fastq", prefix))
    }
}

pub fn default_make_pe_prefix(input: &Path, insert: usize, step: usize) -> PathBuf {
    let parent = input.parent().unwrap_or_else(|| Path::new(""));
    let stem = input.file_name().and_then(|n| n.to_str()).unwrap_or("reads");
    let base = stem
        .trim_end_matches(".fastq.gz")
        .trim_end_matches(".fq.gz")
        .trim_end_matches(".fastq")
        .trim_end_matches(".fq")
        .trim_end_matches(".bam");
    parent.join(format!("{}_{}_{}", base, insert, step))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_make_pe_default_contiguous() {
        let mut fq = NamedTempFile::new().unwrap();
        // 600 bp read
        let read_seq = "A".repeat(150) + &"C".repeat(150) + &"G".repeat(150) + &"T".repeat(150);
        let read_qual = "I".repeat(600);
        writeln!(fq, "@read1\n{}\n+\n{}", read_seq, read_qual).unwrap();

        let out_dir = tempfile::tempdir().unwrap();
        let prefix = out_dir.path().join("pe_out");

        let config = MakePeConfig::new(150, None, None, false, 1).unwrap();
        assert_eq!(config.insert, 300);
        assert_eq!(config.step, 300);

        let summary = make_pe(fq.path(), PeOutputTarget::Prefix(&prefix), &config).unwrap();
        assert_eq!(summary.total_reads, 1);
        assert_eq!(summary.total_pairs, 2);
        assert_eq!(summary.skipped_short, 0);

        let (r1_path, r2_path) = make_pe_output_paths(&prefix).unwrap();
        let mut r1_content = String::new();
        File::open(r1_path).unwrap().read_to_string(&mut r1_content).unwrap();
        let mut r2_content = String::new();
        File::open(r2_path).unwrap().read_to_string(&mut r2_content).unwrap();

        // Check R1
        assert!(r1_content.contains("@read1_f0/1\n"));
        assert!(r1_content.contains(&"A".repeat(150)));
        assert!(r1_content.contains("@read1_f1/1\n"));
        assert!(r1_content.contains(&"G".repeat(150)));

        // Check R2 (revcomp of C is G, revcomp of T is A)
        assert!(r2_content.contains("@read1_f0/2\n"));
        assert!(r2_content.contains(&"G".repeat(150))); // revcomp of C..C
        // Qualities survive FASTQ -> BAM -> pseudo-PE unchanged ('I' = Q40).
        assert_eq!(r1_content.lines().nth(3).unwrap(), "I".repeat(150));
        assert_eq!(r2_content.lines().nth(3).unwrap(), "I".repeat(150));
        assert!(r2_content.contains("@read1_f1/2\n"));
        assert!(r2_content.contains(&"A".repeat(150))); // revcomp of T..T
    }

    #[test]
    fn test_make_pe_sliding_window_1bp() {
        let mut fq = NamedTempFile::new().unwrap();
        // 302 bp read -> with insert 300 and step 1, produces 3 pairs (start=0, 1, 2)
        let read_seq = "A".repeat(302);
        let read_qual = "I".repeat(302);
        writeln!(fq, "@read_tile\n{}\n+\n{}", read_seq, read_qual).unwrap();

        let out_dir = tempfile::tempdir().unwrap();
        let prefix = out_dir.path().join("pe_tile");

        let config = MakePeConfig::new(150, Some(300), Some(1), false, 1).unwrap();
        let summary = make_pe(fq.path(), PeOutputTarget::Prefix(&prefix), &config).unwrap();
        assert_eq!(summary.total_pairs, 3);
    }

    #[test]
    fn test_make_pe_interleaved() {
        let mut fq = NamedTempFile::new().unwrap();
        let read_seq = "A".repeat(150) + &"C".repeat(150);
        let read_qual = "I".repeat(300);
        writeln!(fq, "@read_interleaved\n{}\n+\n{}", read_seq, read_qual).unwrap();

        let out_dir = tempfile::tempdir().unwrap();
        let out_file = out_dir.path().join("interleaved.fastq");

        let config = MakePeConfig::new(150, None, None, true, 1).unwrap();
        let summary = make_pe(fq.path(), PeOutputTarget::SingleFile(&out_file), &config).unwrap();
        assert_eq!(summary.total_pairs, 1);

        let mut content = String::new();
        File::open(out_file).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.contains("@read_interleaved_f0/1\n"));
        assert!(content.contains("@read_interleaved_f0/2\n"));
    }

    #[test]
    fn test_make_pe_bam_input() {
        use crate::bam::{FilterSettings, UnalignedBamWriter};
        use nanoseq_core::header::NanoporeMetadata;

        let out_dir = tempfile::tempdir().unwrap();
        let bam_path = out_dir.path().join("test.bam");

        let metadata = NanoporeMetadata {
            read_id: "bam_read_1".to_string(),
            run_id: Some("test_run".to_string()),
            ..Default::default()
        };
        let settings = FilterSettings {
            qv_threshold: 0.0,
            min_len: 0,
            max_len: usize::MAX,
            channel_range: None,
            time_window: None,
        };
        let mut writer = UnalignedBamWriter::new(&bam_path, &metadata, &settings).unwrap();

        // 300 bp read
        let read_seq = "A".repeat(150) + &"C".repeat(150);
        let read_qual = vec![b'I'; 300]; // Phred ASCII 'I' = 73 -> Phred 40
        writer.write_record(b"bam_read_1", read_seq.as_bytes(), &read_qual).unwrap();
        drop(writer);

        let prefix = out_dir.path().join("bam_pe_out");
        let config = MakePeConfig::new(150, None, None, false, 1).unwrap();
        let summary = make_pe(&bam_path, PeOutputTarget::Prefix(&prefix), &config).unwrap();
        assert_eq!(summary.total_reads, 1);
        assert_eq!(summary.total_pairs, 1);

        let (r1_path, r2_path) = make_pe_output_paths(&prefix).unwrap();
        let mut r1_content = String::new();
        File::open(r1_path).unwrap().read_to_string(&mut r1_content).unwrap();
        let mut r2_content = String::new();
        File::open(r2_path).unwrap().read_to_string(&mut r2_content).unwrap();

        assert!(r1_content.contains("@bam_read_1_f0/1\n"));
        assert!(r1_content.contains(&"A".repeat(150)));
        assert!(r2_content.contains("@bam_read_1_f0/2\n"));
        assert!(r2_content.contains(&"G".repeat(150))); // revcomp of C..C
        // Qualities survive FASTQ -> BAM -> pseudo-PE unchanged ('I' = Q40).
        assert_eq!(r1_content.lines().nth(3).unwrap(), "I".repeat(150));
        assert_eq!(r2_content.lines().nth(3).unwrap(), "I".repeat(150));
    }
}
