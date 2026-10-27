//! Read selection from an ONT `sequencing_summary.txt`.
//!
//! MinKNOW writes one row per read with its barcode, the score of that barcode
//! at the front and at the rear of the read, and the basecaller's mean quality.
//! The FASTQ headers carry only the barcode name, and MinKNOW's BAM only the
//! quality, so the summary is the one place that says whether a read has its
//! barcode at both ends. Selecting on it makes `--barcode-both-ends` possible
//! for a run that was demultiplexed without it, without demultiplexing again.

use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// What the summary says about one read.
#[derive(Debug, Clone, PartialEq)]
pub struct SummaryRead {
    pub mean_qscore: Option<f64>,
    pub barcode: Option<String>,
    pub front_score: Option<f64>,
    pub rear_score: Option<f64>,
}

#[derive(Debug, Default)]
pub struct SequencingSummary {
    reads: HashMap<String, SummaryRead>,
}

/// `barcode13`, `bc13`, `NB13`, `13` -> `barcode13`; anything else as-is, lowercased.
pub fn normalize_barcode(value: &str) -> String {
    let lower = value.trim().to_ascii_lowercase();
    let digits = ["barcode", "bc", "nb", "rb", ""].iter().find_map(|prefix| {
        lower
            .strip_prefix(prefix)
            .filter(|rest| !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()))
    });
    match digits.and_then(|d| d.parse::<u32>().ok()) {
        Some(n) => format!("barcode{n:02}"),
        None => lower,
    }
}

fn parse_score(value: &str) -> Option<f64> {
    value.trim().parse::<f64>().ok().filter(|v| v.is_finite())
}

impl SequencingSummary {
    /// Load the reads of `path` (plain or gzipped). With `only_barcode`, keep
    /// just the reads MinKNOW assigned to that barcode: one barcode's reads are
    /// all a per-barcode filter needs, and a full run can hold tens of millions.
    pub fn load(path: &Path, only_barcode: Option<&str>) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open sequencing summary {}", path.display()))?;
        let reader: Box<dyn BufRead> = if path.to_string_lossy().to_lowercase().ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        Self::from_reader(reader, only_barcode)
    }

    pub fn from_reader(mut reader: impl BufRead, only_barcode: Option<&str>) -> Result<Self> {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            return Err(anyhow!("Empty sequencing summary"));
        }
        let headers: Vec<String> = line.trim_end().split('\t').map(|h| h.trim().to_string()).collect();
        let column = |name: &str| headers.iter().position(|h| h == name);
        let read_id = column("read_id").ok_or_else(|| anyhow!("sequencing summary has no read_id column"))?;
        let qscore = column("mean_qscore_template");
        let barcode = column("barcode_arrangement");
        let front = column("barcode_front_score");
        let rear = column("barcode_rear_score");
        let wanted = only_barcode.map(normalize_barcode);
        if wanted.is_some() && barcode.is_none() {
            return Err(anyhow!("sequencing summary has no barcode_arrangement column"));
        }

        let mut reads = HashMap::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let fields: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
            let Some(id) = fields.get(read_id).map(|s| s.trim()).filter(|s| !s.is_empty()) else {
                continue;
            };
            let read_barcode = barcode
                .and_then(|i| fields.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty() && *s != "-")
                .map(normalize_barcode);
            if let Some(wanted) = &wanted {
                if read_barcode.as_deref() != Some(wanted.as_str()) {
                    continue;
                }
            }
            reads.insert(
                id.to_string(),
                SummaryRead {
                    mean_qscore: qscore.and_then(|i| fields.get(i)).and_then(|s| parse_score(s)),
                    barcode: read_barcode,
                    front_score: front.and_then(|i| fields.get(i)).and_then(|s| parse_score(s)),
                    rear_score: rear.and_then(|i| fields.get(i)).and_then(|s| parse_score(s)),
                },
            );
        }
        Ok(Self { reads })
    }

    pub fn get(&self, read_id: &str) -> Option<&SummaryRead> {
        self.reads.get(read_id)
    }

    pub fn len(&self) -> usize {
        self.reads.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }
}

/// Why a read was kept or dropped by a [`ReadSelection`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Verdict {
    /// Kept. `mean_qscore` is the summary's, which overrides the read's own.
    Keep { mean_qscore: Option<f64> },
    /// The summary does not list the read; with `--both-ends` its ends cannot
    /// be checked, so it is dropped.
    NotInSummary,
    /// MinKNOW assigned the read to another barcode than `--barcode`.
    OtherBarcode,
    /// The barcode was not found at both ends with at least the minimum score.
    NotBothEnds,
}

/// Which reads of a summary to keep, and which quality to judge them by.
#[derive(Debug)]
pub struct ReadSelection {
    pub summary: SequencingSummary,
    /// `Some(score)`: require the assigned barcode at the front and at the rear,
    /// each scoring at least `score`. The rear score is that of the barcode the
    /// read was assigned, so a read whose rear holds another barcode fails too.
    pub both_ends_min_score: Option<f64>,
    /// Keep only reads MinKNOW assigned to this barcode.
    pub barcode: Option<String>,
}

impl ReadSelection {
    pub fn judge(&self, read_id: &str) -> Verdict {
        let Some(read) = self.summary.get(read_id) else {
            return if self.both_ends_min_score.is_some() || self.barcode.is_some() {
                Verdict::NotInSummary
            } else {
                Verdict::Keep { mean_qscore: None }
            };
        };
        if let Some(wanted) = &self.barcode {
            if read.barcode.as_deref() != Some(normalize_barcode(wanted).as_str()) {
                return Verdict::OtherBarcode;
            }
        }
        if let Some(min) = self.both_ends_min_score {
            let both = matches!((read.front_score, read.rear_score), (Some(f), Some(r)) if f >= min && r >= min);
            if !both {
                return Verdict::NotBothEnds;
            }
        }
        Verdict::Keep { mean_qscore: read.mean_qscore }
    }
}

/// Counts of what a selection dropped, for the filter's report.
#[derive(Debug, Default, Clone, Copy)]
pub struct SelectionStats {
    pub not_in_summary: u64,
    pub other_barcode: u64,
    pub not_both_ends: u64,
}

impl SelectionStats {
    pub fn record(&mut self, verdict: Verdict) {
        match verdict {
            Verdict::Keep { .. } => {}
            Verdict::NotInSummary => self.not_in_summary += 1,
            Verdict::OtherBarcode => self.other_barcode += 1,
            Verdict::NotBothEnds => self.not_both_ends += 1,
        }
    }

    pub fn report(&self) -> String {
        format!(
            "dropped by the sequencing summary: {} not both ends, {} other barcode, {} not in summary",
            self.not_both_ends, self.other_barcode, self.not_in_summary
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SUMMARY: &str = "read_id\tmean_qscore_template\tbarcode_arrangement\tbarcode_front_score\tbarcode_rear_score\n\
        both\t18.5\tbarcode13\t97.3\t88.0\n\
        front_only\t21.0\tbarcode13\t100.0\t46.2\n\
        rear_only\t19.0\tbarcode13\t40.0\t95.0\n\
        other\t20.0\tbarcode15\t99.0\t99.0\n\
        unclassified\t12.0\tunclassified\t-\t-\n";

    fn selection(both_ends: Option<f64>, barcode: Option<&str>) -> ReadSelection {
        ReadSelection {
            summary: SequencingSummary::from_reader(SUMMARY.as_bytes(), None).unwrap(),
            both_ends_min_score: both_ends,
            barcode: barcode.map(str::to_string),
        }
    }

    #[test]
    fn both_ends_keeps_only_reads_scored_at_front_and_rear() {
        let s = selection(Some(60.0), None);
        assert_eq!(s.judge("both"), Verdict::Keep { mean_qscore: Some(18.5) });
        assert_eq!(s.judge("front_only"), Verdict::NotBothEnds);
        assert_eq!(s.judge("rear_only"), Verdict::NotBothEnds);
        assert_eq!(s.judge("unclassified"), Verdict::NotBothEnds);
        assert_eq!(s.judge("missing"), Verdict::NotInSummary);
    }

    #[test]
    fn barcode_must_match_the_assigned_one() {
        let s = selection(None, Some("bc13"));
        assert_eq!(s.judge("front_only"), Verdict::Keep { mean_qscore: Some(21.0) });
        assert_eq!(s.judge("other"), Verdict::OtherBarcode);
    }

    #[test]
    fn without_both_ends_or_barcode_it_only_supplies_the_quality() {
        let s = selection(None, None);
        assert_eq!(s.judge("front_only"), Verdict::Keep { mean_qscore: Some(21.0) });
        assert_eq!(s.judge("missing"), Verdict::Keep { mean_qscore: None });
    }

    #[test]
    fn loading_one_barcode_keeps_only_its_reads() {
        let summary = SequencingSummary::from_reader(SUMMARY.as_bytes(), Some("NB13")).unwrap();
        assert_eq!(summary.len(), 3);
        assert!(summary.get("other").is_none());
    }

    #[test]
    fn normalizes_barcode_spellings() {
        assert_eq!(normalize_barcode("barcode05"), "barcode05");
        assert_eq!(normalize_barcode("BC5"), "barcode05");
        assert_eq!(normalize_barcode("NB13"), "barcode13");
        assert_eq!(normalize_barcode("13"), "barcode13");
        assert_eq!(normalize_barcode("unclassified"), "unclassified");
    }
}
