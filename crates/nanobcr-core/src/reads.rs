//! Reading one FASTQ file and applying the same filter as `nanostream filter`.

use crate::watcher::barcode_from_path;
use anyhow::{anyhow, Context, Result};
use nanofilter_core::bam::FilterSettings;
use nanofilter_core::fastq::fastq_record_passes;
use nanoseq_core::header::parse_nanopore_header;
use needletail::parse_fastx_file;
use std::collections::BTreeMap;
use std::path::Path;

/// A read that passed the filter, kept exactly as MinKNOW wrote it (the whole
/// header line, so minimap2 `-y` carries the tags into the BAM).
#[derive(Debug, Clone)]
pub struct Read {
    pub barcode: String,
    /// Header line without the leading `@`.
    pub header: Vec<u8>,
    /// First whitespace-delimited token of the header.
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadCounts {
    pub seen: u64,
    pub passed: u64,
    pub failed_len: u64,
    pub failed_qv: u64,
    pub bases_passed: u64,
}

impl ReadCounts {
    pub fn add(&mut self, other: &ReadCounts) {
        self.seen += other.seen;
        self.passed += other.passed;
        self.failed_len += other.failed_len;
        self.failed_qv += other.failed_qv;
        self.bases_passed += other.bases_passed;
    }
}

pub fn filter_settings(min_qv: f64, min_len: usize, max_len: usize) -> FilterSettings {
    FilterSettings {
        qv_threshold: min_qv,
        min_len,
        max_len,
        channel_range: None,
        time_window: None,
    }
}

pub struct LoadedFile {
    pub reads: Vec<Read>,
    pub counts: BTreeMap<String, ReadCounts>,
}

/// Read the whole file. Any parse error (a truncated gzip, a record cut in half)
/// fails the call and nothing is returned, so a file that is still being written
/// is retried later instead of being half processed.
pub fn load_fastq(path: &Path, settings: &FilterSettings) -> Result<LoadedFile> {
    let path_barcode = barcode_from_path(path);
    let mut reader =
        parse_fastx_file(path).with_context(|| format!("opening {}", path.display()))?;
    let mut reads = Vec::new();
    let mut counts: BTreeMap<String, ReadCounts> = BTreeMap::new();

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("reading {}", path.display()))?;
        let header = record.id();
        let qual = record
            .qual()
            .ok_or_else(|| anyhow!("{}: record without qualities", path.display()))?;
        let seq = record.seq();
        if seq.len() != qual.len() {
            return Err(anyhow!(
                "{}: sequence and quality lengths differ in {}",
                path.display(),
                String::from_utf8_lossy(header)
            ));
        }

        let metadata = parse_nanopore_header(header);
        let barcode = path_barcode
            .clone()
            .or_else(|| metadata.barcode.clone())
            .unwrap_or_else(|| "all".to_string());
        let c = counts.entry(barcode.clone()).or_default();
        c.seen += 1;

        if fastq_record_passes(settings, &metadata, seq.len(), qual) {
            c.passed += 1;
            c.bases_passed += seq.len() as u64;
            let name = header
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or(header)
                .to_vec();
            reads.push(Read {
                barcode,
                header: header.to_vec(),
                name,
                seq: seq.to_vec(),
                qual: qual.to_vec(),
            });
        } else if seq.len() < settings.min_len || seq.len() > settings.max_len {
            c.failed_len += 1;
        } else {
            c.failed_qv += 1;
        }
    }
    Ok(LoadedFile { reads, counts })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn fastq_record(name: &str, comment: &str, len: usize, q: u8) -> String {
        format!(
            "@{name} {comment}\n{}\n+\n{}\n",
            &"ACGT".repeat(len / 4 + 1)[..len],
            (q as char).to_string().repeat(len)
        )
    }

    #[test]
    fn filters_like_nanostream_and_keeps_the_header() {
        let dir = tempfile::tempdir().unwrap();
        let bc = dir.path().join("fastq_pass/barcode03");
        std::fs::create_dir_all(&bc).unwrap();
        let path = bc.join("x.fastq");
        let mut f = std::fs::File::create(&path).unwrap();
        // qs tag decides quality, not the quality string: Q40 string but qs 12.
        write!(f, "{}", fastq_record("good", "qs:f:18.5 ch:i:4 SM:Z:barcode03", 1000, b'I')).unwrap();
        write!(f, "{}", fastq_record("lowqs", "qs:f:12.0 ch:i:5", 1000, b'I')).unwrap();
        write!(f, "{}", fastq_record("short", "qs:f:25.0", 100, b'I')).unwrap();
        write!(f, "{}", fastq_record("long", "qs:f:25.0", 4000, b'I')).unwrap();
        // No qs tag: falls back to the quality string (Q40 passes).
        write!(f, "{}", fastq_record("noqs", "", 800, b'I')).unwrap();
        drop(f);

        let loaded = load_fastq(&path, &filter_settings(20.0, 500, 3000)).unwrap();
        // 'good' has qs 18.5 < 20 -> fails; only the untagged Q40 read passes.
        let names: Vec<_> = loaded.reads.iter().map(|r| String::from_utf8_lossy(&r.name).to_string()).collect();
        assert_eq!(names, vec!["noqs"]);
        let c = &loaded.counts["barcode03"];
        assert_eq!((c.seen, c.passed, c.failed_len, c.failed_qv), (5, 1, 2, 2));
        assert_eq!(c.bases_passed, 800);

        let loaded = load_fastq(&path, &filter_settings(18.0, 500, 3000)).unwrap();
        let r = &loaded.reads[0];
        assert_eq!(r.header, b"good qs:f:18.5 ch:i:4 SM:Z:barcode03");
        assert_eq!(r.name, b"good");
        assert_eq!(r.barcode, "barcode03");
    }

    #[test]
    fn truncated_gzip_is_an_error_not_a_partial_load() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.fastq.gz");
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        for i in 0..200 {
            write!(gz, "{}", fastq_record(&format!("r{i}"), "qs:f:25.0", 1000, b'I')).unwrap();
        }
        let bytes = gz.finish().unwrap();
        std::fs::write(&path, &bytes[..bytes.len() / 2]).unwrap();
        assert!(load_fastq(&path, &filter_settings(10.0, 0, 5000)).is_err());

        std::fs::write(&path, &bytes).unwrap();
        assert_eq!(load_fastq(&path, &filter_settings(10.0, 0, 5000)).unwrap().reads.len(), 200);
    }
}
