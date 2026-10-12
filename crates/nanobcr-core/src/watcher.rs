//! Finds FASTQ files MinKNOW has finished writing.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Candidate {
    pub path: PathBuf,
    /// Relative to the watched directory; the ledger key.
    pub rel_path: String,
    pub size: u64,
}

pub fn is_fastq(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    n.ends_with(".fastq") || n.ends_with(".fq") || n.ends_with(".fastq.gz") || n.ends_with(".fq.gz")
}

/// `barcode01`, `barcode24`, `unclassified`: the directory MinKNOW put a read in.
pub fn barcode_from_path(path: &Path) -> Option<String> {
    for comp in path.components() {
        let name = comp.as_os_str().to_string_lossy();
        let lower = name.to_ascii_lowercase();
        if lower == "unclassified" {
            return Some(lower);
        }
        if let Some(rest) = lower.strip_prefix("barcode") {
            if !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()) {
                return Some(lower);
            }
        }
    }
    None
}

pub struct Watcher {
    root: PathBuf,
    include_fail: bool,
    stable: Duration,
    last_seen: HashMap<PathBuf, (u64, SystemTime)>,
    waiting: usize,
}

impl Watcher {
    pub fn new(root: &Path, include_fail: bool, stable_secs: u64) -> Self {
        Self {
            root: root.to_path_buf(),
            include_fail,
            stable: Duration::from_secs(stable_secs),
            last_seen: HashMap::new(),
            waiting: 0,
        }
    }

    /// Candidates seen by the last scan that were not yet ready (too young, still
    /// changing, or empty). Lets a one-shot run know whether anything is left.
    pub fn waiting(&self) -> usize {
        self.waiting
    }

    /// Files that look finished: non-empty, unchanged since the previous scan and
    /// not modified for `stable_secs`. A writer that pauses longer than that can
    /// still be caught mid-file; the reader treats a truncated file as not ready.
    pub fn scan(&mut self) -> Vec<Candidate> {
        let mut found = Vec::new();
        self.walk(&self.root.clone(), &mut found);
        found.sort_by(|a, b| a.0.cmp(&b.0));

        let now = SystemTime::now();
        let mut ready = Vec::new();
        let mut waiting = 0;
        let mut seen_now = HashMap::new();
        for (path, size, mtime) in found {
            let age = now.duration_since(mtime).unwrap_or(Duration::ZERO);
            let unchanged = self
                .last_seen
                .get(&path)
                .map(|&(s, m)| s == size && m == mtime)
                .unwrap_or(true);
            if size > 0 && unchanged && age >= self.stable {
                let rel = path
                    .strip_prefix(&self.root)
                    .unwrap_or(&path)
                    .to_string_lossy()
                    .to_string();
                ready.push(Candidate {
                    path: path.clone(),
                    rel_path: rel,
                    size,
                });
            } else {
                waiting += 1;
            }
            seen_now.insert(path, (size, mtime));
        }
        self.last_seen = seen_now;
        self.waiting = waiting;
        ready
    }

    fn walk(&self, dir: &Path, out: &mut Vec<(PathBuf, u64, SystemTime)>) {
        let Ok(entries) = fs::read_dir(dir) else {
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            let name = entry.file_name().to_string_lossy().to_string();
            if name.starts_with('.') {
                continue;
            }
            let Ok(meta) = entry.metadata() else { continue };
            if meta.is_dir() {
                if !self.include_fail && name == "fastq_fail" {
                    continue;
                }
                self.walk(&path, out);
            } else if meta.is_file() && is_fastq(&name) {
                if let Ok(mtime) = meta.modified() {
                    out.push((path, meta.len(), mtime));
                }
            }
        }
    }
}

/// MinKNOW writes `final_summary_*.txt` next to `fastq_pass` when the run ends.
pub fn run_finished(input_dir: &Path) -> bool {
    let has_summary = |dir: &Path| {
        fs::read_dir(dir)
            .map(|rd| {
                rd.flatten().any(|e| {
                    let n = e.file_name().to_string_lossy().to_string();
                    n.starts_with("final_summary") && n.ends_with(".txt")
                })
            })
            .unwrap_or(false)
    };
    has_summary(input_dir) || input_dir.parent().map(has_summary).unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn barcode_from_directory_name() {
        let p = Path::new("/run/fastq_pass/barcode07/PBE_pass_barcode07_0.fastq.gz");
        assert_eq!(barcode_from_path(p).as_deref(), Some("barcode07"));
        assert_eq!(
            barcode_from_path(Path::new("fastq_fail/unclassified/x.fastq")).as_deref(),
            Some("unclassified")
        );
        assert_eq!(barcode_from_path(Path::new("fastq_pass/x.fastq")), None);
        assert_eq!(barcode_from_path(Path::new("barcodes_x/y.fastq")), None);
    }

    #[test]
    fn recent_and_changing_files_are_not_ready() {
        let dir = tempfile::tempdir().unwrap();
        let bc = dir.path().join("fastq_pass/barcode01");
        fs::create_dir_all(&bc).unwrap();
        fs::write(bc.join("a.fastq"), "@r\nA\n+\nI\n").unwrap();
        fs::write(bc.join(".hidden.fastq"), "x").unwrap();
        fs::write(bc.join("notes.txt"), "x").unwrap();

        // Just written: too young for a 3600 s stability window.
        let mut w = Watcher::new(dir.path(), true, 3600);
        assert!(w.scan().is_empty());

        // No stability requirement: ready, with a path relative to the root.
        let mut w = Watcher::new(dir.path(), true, 0);
        let ready = w.scan();
        assert_eq!(ready.len(), 1);
        assert_eq!(ready[0].rel_path, "fastq_pass/barcode01/a.fastq");

        // A file that grew between scans is held back for one more round.
        fs::write(bc.join("a.fastq"), "@r\nA\n+\nI\n@s\nA\n+\nI\n").unwrap();
        assert!(w.scan().is_empty());
        assert_eq!(w.scan().len(), 1);
    }

    #[test]
    fn fail_directory_can_be_excluded() {
        let dir = tempfile::tempdir().unwrap();
        for d in ["fastq_pass/barcode01", "fastq_fail/barcode01"] {
            fs::create_dir_all(dir.path().join(d)).unwrap();
            fs::write(dir.path().join(d).join("a.fastq"), "@r\nA\n+\nI\n").unwrap();
        }
        let mut w = Watcher::new(dir.path(), false, 0);
        let rels: Vec<_> = w.scan().into_iter().map(|c| c.rel_path).collect();
        assert_eq!(rels, vec!["fastq_pass/barcode01/a.fastq".to_string()]);
    }
}
