//! Feeds a finished run into a directory the way MinKNOW would, so the engine
//! (and later the GUI) can be exercised without a sequencer.

use crate::watcher::is_fastq;
use anyhow::{Context, Result};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, SystemTime};

#[derive(Debug, Clone)]
pub struct ReplayConfig {
    pub src: PathBuf,
    pub dest: PathBuf,
    /// Seconds between files.
    pub interval_secs: f64,
    /// Write each file in two halves with a pause between, so it is observably
    /// "being written" (a truncated gzip) for a moment.
    pub partial: bool,
    pub partial_pause_secs: f64,
    /// Write `final_summary_replay.txt` in `dest` when done, as MinKNOW does.
    pub finish_marker: bool,
}

impl ReplayConfig {
    pub fn new(src: PathBuf, dest: PathBuf) -> Self {
        Self {
            src,
            dest,
            interval_secs: 4.0,
            partial: false,
            partial_pause_secs: 1.5,
            finish_marker: false,
        }
    }
}

fn collect(dir: &Path, out: &mut Vec<(PathBuf, SystemTime)>) {
    let Ok(rd) = fs::read_dir(dir) else { return };
    for e in rd.flatten() {
        let p = e.path();
        let Ok(meta) = e.metadata() else { continue };
        if meta.is_dir() {
            collect(&p, out);
        } else if is_fastq(&e.file_name().to_string_lossy()) {
            out.push((p, meta.modified().unwrap_or(SystemTime::UNIX_EPOCH)));
        }
    }
}

fn sleep_unless_stopped(secs: f64, stop: &AtomicBool) {
    let mut left = Duration::from_secs_f64(secs.max(0.0));
    while left > Duration::ZERO && !stop.load(Ordering::Relaxed) {
        let step = left.min(Duration::from_millis(200));
        std::thread::sleep(step);
        left -= step;
    }
}

/// Copy every FASTQ under `src` into `dest`, in modification-time order and with
/// the same relative paths. Returns the number of files written.
pub fn replay(cfg: &ReplayConfig, stop: &AtomicBool) -> Result<usize> {
    let mut files = Vec::new();
    collect(&cfg.src, &mut files);
    files.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
    fs::create_dir_all(&cfg.dest)?;

    let mut written = 0;
    for (i, (path, _)) in files.iter().enumerate() {
        if stop.load(Ordering::Relaxed) {
            break;
        }
        if i > 0 {
            sleep_unless_stopped(cfg.interval_secs, stop);
        }
        let rel = path.strip_prefix(&cfg.src).unwrap_or(path);
        let dest = cfg.dest.join(rel);
        fs::create_dir_all(dest.parent().unwrap())?;
        let bytes = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        // Written under a hidden name first when not partial, so the watcher never sees
        // a half-copied file; partial mode deliberately exposes one.
        if cfg.partial {
            let mut f = fs::File::create(&dest)?;
            f.write_all(&bytes[..bytes.len() / 2])?;
            f.flush()?;
            sleep_unless_stopped(cfg.partial_pause_secs, stop);
            f.write_all(&bytes[bytes.len() / 2..])?;
            f.flush()?;
        } else {
            let tmp = dest.with_file_name(format!(".{}", dest.file_name().unwrap().to_string_lossy()));
            fs::write(&tmp, &bytes)?;
            fs::rename(&tmp, &dest)?;
        }
        written += 1;
        log::info!("replayed {} ({written}/{})", rel.display(), files.len());
    }
    if cfg.finish_marker && !stop.load(Ordering::Relaxed) {
        fs::write(cfg.dest.join("final_summary_replay.txt"), "replay finished\n")?;
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copies_in_order_with_relative_paths_and_marker() {
        let src = tempfile::tempdir().unwrap();
        let dest = tempfile::tempdir().unwrap();
        for (rel, body) in [("fastq_pass/barcode01/a.fastq", "A"), ("fastq_pass/barcode02/b.fastq", "B"), ("other.txt", "x")] {
            let p = src.path().join(rel);
            fs::create_dir_all(p.parent().unwrap()).unwrap();
            fs::write(p, body).unwrap();
        }
        let mut cfg = ReplayConfig::new(src.path().into(), dest.path().join("run"));
        cfg.interval_secs = 0.0;
        cfg.partial = true;
        cfg.partial_pause_secs = 0.0;
        cfg.finish_marker = true;
        assert_eq!(replay(&cfg, &AtomicBool::new(false)).unwrap(), 2);
        assert_eq!(fs::read(dest.path().join("run/fastq_pass/barcode01/a.fastq")).unwrap(), b"A");
        assert_eq!(fs::read(dest.path().join("run/fastq_pass/barcode02/b.fastq")).unwrap(), b"B");
        assert!(!dest.path().join("run/other.txt").exists());
        assert!(dest.path().join("run/final_summary_replay.txt").exists());
    }
}
