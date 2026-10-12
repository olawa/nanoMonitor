//! Append-only record of which input files have been processed, so a restarted
//! engine continues instead of re-aligning (and double-counting) everything.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct LedgerEntry {
    /// Path relative to the input directory.
    pub path: String,
    pub size: u64,
    pub reads: u64,
    pub passed: u64,
    pub ts_unix: u64,
}

pub struct Ledger {
    path: PathBuf,
    done: HashMap<String, LedgerEntry>,
}

impl Ledger {
    pub fn open(path: &Path) -> Result<Self> {
        let mut done = HashMap::new();
        if path.exists() {
            let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
            for (n, line) in BufReader::new(file).lines().enumerate() {
                let line = line?;
                if line.trim().is_empty() {
                    continue;
                }
                match serde_json::from_str::<LedgerEntry>(&line) {
                    Ok(entry) => {
                        done.insert(entry.path.clone(), entry);
                    }
                    // A torn last line from a crash: ignore it, the file is simply redone.
                    Err(e) => log::warn!("ledger line {} unreadable ({e}); ignoring", n + 1),
                }
            }
        }
        // After a crash the last line may lack its newline; terminate it so the next
        // record does not get glued onto the torn one.
        if let Ok(bytes) = std::fs::read(path) {
            if bytes.last().map(|&b| b != b'\n').unwrap_or(false) {
                OpenOptions::new().append(true).open(path)?.write_all(b"\n")?;
            }
        }
        Ok(Self {
            path: path.to_path_buf(),
            done,
        })
    }

    pub fn is_done(&self, rel_path: &str) -> bool {
        self.done.contains_key(rel_path)
    }

    pub fn get(&self, rel_path: &str) -> Option<&LedgerEntry> {
        self.done.get(rel_path)
    }

    pub fn len(&self) -> usize {
        self.done.len()
    }

    pub fn is_empty(&self) -> bool {
        self.done.is_empty()
    }

    pub fn record(&mut self, entry: LedgerEntry) -> Result<()> {
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        writeln!(file, "{}", serde_json::to_string(&entry)?)?;
        file.flush()?;
        self.done.insert(entry.path.clone(), entry);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(path: &str) -> LedgerEntry {
        LedgerEntry {
            path: path.into(),
            size: 10,
            reads: 4,
            passed: 3,
            ts_unix: 1,
        }
    }

    #[test]
    fn survives_a_restart_and_a_torn_line() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ledger.jsonl");
        let mut l = Ledger::open(&p).unwrap();
        l.record(entry("barcode01/a.fastq.gz")).unwrap();
        l.record(entry("barcode01/b.fastq.gz")).unwrap();
        // Simulate a crash mid-write.
        std::fs::OpenOptions::new()
            .append(true)
            .open(&p)
            .unwrap()
            .write_all(b"{\"path\":\"barcode01/c.fa")
            .unwrap();

        let mut l = Ledger::open(&p).unwrap();
        l.record(entry("barcode01/d.fastq.gz")).unwrap();
        let l = Ledger::open(&p).unwrap();
        assert!(l.is_done("barcode01/d.fastq.gz"));
        assert!(l.is_done("barcode01/a.fastq.gz"));
        assert!(l.is_done("barcode01/b.fastq.gz"));
        assert!(!l.is_done("barcode01/c.fastq.gz"));
        assert_eq!(l.len(), 3);
    }
}
