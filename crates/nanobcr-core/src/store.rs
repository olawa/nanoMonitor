//! Per-barcode BAM store.
//!
//! `inbox/` receives unsorted segments from alignment cycles. `compact` folds them
//! into a sorted, indexed `compacted.<generation>.bam`. The manifest names the
//! current BAM and the segments it already contains, and is replaced atomically
//! *before* anything is deleted, so a crash at any point loses nothing and
//! double-counts nothing.

use crate::align::run_checked;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct Manifest {
    pub generation: u64,
    pub bam: Option<String>,
    /// Segments already folded into `bam` (deleted, or about to be).
    pub consumed: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Compacted {
    pub bam: PathBuf,
    pub generation: u64,
    pub segments_folded: usize,
}

pub struct Store {
    dir: PathBuf,
}

pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let tmp = path.with_extension("tmp");
    fs::write(&tmp, bytes)?;
    fs::rename(&tmp, path)?;
    Ok(())
}

impl Store {
    pub fn new(dir: &Path) -> Self {
        Self { dir: dir.to_path_buf() }
    }

    pub fn inbox(&self) -> PathBuf {
        self.dir.join("inbox")
    }

    fn manifest_path(&self) -> PathBuf {
        self.dir.join("manifest.json")
    }

    pub fn manifest(&self) -> Manifest {
        fs::read(self.manifest_path())
            .ok()
            .and_then(|b| serde_json::from_slice(&b).ok())
            .unwrap_or_default()
    }

    pub fn current_bam(&self) -> Option<PathBuf> {
        self.manifest().bam.map(|n| self.dir.join(n))
    }

    /// Segments waiting to be compacted, oldest first.
    pub fn pending(&self) -> Vec<PathBuf> {
        let consumed = self.manifest().consumed;
        let mut v: Vec<PathBuf> = fs::read_dir(self.inbox())
            .map(|rd| {
                rd.flatten()
                    .map(|e| e.path())
                    .filter(|p| {
                        let name = p.file_name().unwrap().to_string_lossy().to_string();
                        name.ends_with(".bam") && !name.starts_with('.') && !consumed.contains(&name)
                    })
                    .collect()
            })
            .unwrap_or_default();
        v.sort();
        v
    }

    /// Clean up after a crash: delete segments the manifest says are already
    /// folded in, half-written temporaries, and compacted files that are not current.
    pub fn recover(&self) -> Result<()> {
        let m = self.manifest();
        if let Ok(rd) = fs::read_dir(self.inbox()) {
            for e in rd.flatten() {
                let name = e.file_name().to_string_lossy().to_string();
                if name.starts_with('.') || m.consumed.contains(&name) {
                    let _ = fs::remove_file(e.path());
                }
            }
        }
        if let Ok(rd) = fs::read_dir(&self.dir) {
            for e in rd.flatten() {
                let name = e.file_name().to_string_lossy().to_string();
                let is_compacted = name.starts_with("compacted.");
                let is_current = m.bam.as_deref().map(|b| name == b || name == format!("{b}.bai")).unwrap_or(false);
                if name.starts_with('.') || (is_compacted && !is_current) {
                    let _ = fs::remove_file(e.path());
                }
            }
        }
        Ok(())
    }

    /// Fold all pending segments into the compacted BAM. `None` if there are none.
    pub fn compact(&self, samtools: &str, threads: usize) -> Result<Option<Compacted>> {
        let segments = self.pending();
        if segments.is_empty() {
            return Ok(None);
        }
        let t = threads.max(1).to_string();
        let work = self.dir.join(".work");
        let _ = fs::remove_dir_all(&work);
        fs::create_dir_all(&work)?;

        let to_s = |p: &Path| p.to_string_lossy().to_string();
        let sorted_new = work.join("new.sorted.bam");
        if segments.len() == 1 {
            run_checked(
                samtools,
                &["sort".into(), "-@".into(), t.clone(), "-o".into(), to_s(&sorted_new), to_s(&segments[0])],
            )?;
        } else {
            let catted = work.join("new.cat.bam");
            let mut args = vec!["cat".to_string(), "-o".into(), to_s(&catted)];
            args.extend(segments.iter().map(|s| to_s(s)));
            run_checked(samtools, &args)?;
            run_checked(
                samtools,
                &["sort".into(), "-@".into(), t.clone(), "-o".into(), to_s(&sorted_new), to_s(&catted)],
            )?;
            let _ = fs::remove_file(&catted);
        }

        let old = self.manifest();
        let generation = old.generation + 1;
        let name = format!("compacted.{generation:06}.bam");
        let dest = self.dir.join(&name);
        let tmp_dest = self.dir.join(format!(".{name}"));
        match old.bam.as_ref().map(|b| self.dir.join(b)) {
            Some(prev) if prev.exists() => run_checked(
                samtools,
                &["merge".into(), "-f".into(), "-@".into(), t.clone(), to_s(&tmp_dest), to_s(&prev), to_s(&sorted_new)],
            )?,
            _ => fs::copy(&sorted_new, &tmp_dest).map(|_| ())?,
        }
        fs::rename(&tmp_dest, &dest)?;
        run_checked(samtools, &["index".into(), "-@".into(), t, to_s(&dest)])?;
        let _ = fs::remove_dir_all(&work);

        let consumed: Vec<String> = segments
            .iter()
            .map(|s| s.file_name().unwrap().to_string_lossy().to_string())
            .collect();
        let manifest = Manifest {
            generation,
            bam: Some(name),
            consumed: consumed.clone(),
        };
        // Commit point: from here the new BAM is the truth.
        write_atomic(&self.manifest_path(), &serde_json::to_vec_pretty(&manifest)?)?;

        for s in &segments {
            let _ = fs::remove_file(s);
        }
        if let Some(prev) = old.bam {
            let _ = fs::remove_file(self.dir.join(&prev));
            let _ = fs::remove_file(self.dir.join(format!("{prev}.bai")));
        }
        Ok(Some(Compacted {
            bam: dest,
            generation,
            segments_folded: segments.len(),
        }))
    }
}
