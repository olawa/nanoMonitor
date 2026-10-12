//! Checks worth showing before a run starts, so a missing tool, an unusable
//! reference or too little memory is found in seconds rather than at the first
//! alignment cycle.

use crate::align::tool_version;
use crate::config::EngineConfig;
use crate::watcher::is_fastq;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;
use std::process::Command;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CheckStatus {
    Pass,
    Warn,
    Fail,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Check {
    pub name: String,
    pub status: CheckStatus,
    pub detail: String,
}

fn check(name: &str, status: CheckStatus, detail: impl Into<String>) -> Check {
    Check { name: name.into(), status, detail: detail.into() }
}

/// Number of FASTQ files below `dir`, stopping early at `cap` (a MinKNOW run can
/// hold thousands).
pub fn count_fastq(dir: &Path, cap: usize) -> usize {
    let mut n = 0;
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        let Ok(rd) = fs::read_dir(&d) else { continue };
        for e in rd.flatten() {
            let Ok(meta) = e.metadata() else { continue };
            if meta.is_dir() {
                stack.push(e.path());
            } else if is_fastq(&e.file_name().to_string_lossy()) {
                n += 1;
                if n >= cap {
                    return n;
                }
            }
        }
    }
    n
}

/// Free space in GB on the filesystem holding `path`, via `df -Pk`.
pub fn free_disk_gb(path: &Path) -> Option<f64> {
    let out = Command::new("df").args(["-Pk"]).arg(path).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    let kb: f64 = text.lines().nth(1)?.split_whitespace().nth(3)?.parse().ok()?;
    Some(kb / 1024.0 / 1024.0)
}

/// Installed RAM in GB (macOS `sysctl`, Linux `/proc/meminfo`).
pub fn total_ram_gb() -> Option<f64> {
    if let Ok(out) = Command::new("sysctl").args(["-n", "hw.memsize"]).output() {
        if let Ok(b) = String::from_utf8_lossy(&out.stdout).trim().parse::<f64>() {
            return Some(b / 1e9);
        }
    }
    let text = fs::read_to_string("/proc/meminfo").ok()?;
    let kb: f64 = text.lines().find(|l| l.starts_with("MemTotal:"))?.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb / 1e6)
}

/// Peak memory of one alignment cycle, from what was measured on hg19: a 11.6 GB
/// splice index peaked at 14-18 GB, and building it from the 3.1 GB FASTA at 23 GB.
pub fn estimated_cycle_ram_gb(reference: &Path) -> f64 {
    let size_gb = fs::metadata(reference).map(|m| m.len() as f64 / 1e9).unwrap_or(0.0);
    if reference.extension().map(|e| e == "mmi").unwrap_or(false) {
        size_gb * 1.6
    } else {
        // The index is built once (peak ~7x the FASTA) and then loaded (~4x).
        size_gb * 7.0
    }
}

pub fn preflight(cfg: &EngineConfig) -> Vec<Check> {
    use CheckStatus::*;
    let mut v = Vec::new();

    if cfg.input_dir.is_dir() {
        let n = count_fastq(&cfg.input_dir, 100_000);
        if n > 0 {
            v.push(check("MinKNOW directory", Pass, format!("{} FASTQ files found", n)));
        } else {
            v.push(check(
                "MinKNOW directory",
                Warn,
                "no FASTQ files yet; the engine will wait for them (is this the run folder or fastq_pass?)",
            ));
        }
    } else {
        v.push(check("MinKNOW directory", Fail, format!("{} does not exist", cfg.input_dir.display())));
    }

    for (label, prog) in [("minimap2", &cfg.minimap2), ("samtools", &cfg.samtools)] {
        match tool_version(prog) {
            Ok(ver) if !ver.is_empty() => v.push(check(label, Pass, ver)),
            Ok(_) => v.push(check(label, Warn, format!("`{prog}` ran but printed no version"))),
            Err(e) => v.push(check(label, Fail, format!("{e:#}"))),
        }
    }

    let is_mmi = cfg.reference.extension().map(|e| e == "mmi").unwrap_or(false);
    if !cfg.reference.is_file() {
        v.push(check("Reference", Fail, format!("{} does not exist", cfg.reference.display())));
    } else if is_mmi {
        v.push(check("Reference", Pass, "prebuilt minimap2 index (it must have been built with -x splice)"));
    } else {
        let mut fai = cfg.reference.clone().into_os_string();
        fai.push(".fai");
        let have_fai = Path::new(&fai).is_file();
        let idx = cfg.out_dir.join("index");
        let cached = fs::read_dir(&idx)
            .map(|rd| rd.flatten().any(|e| e.file_name().to_string_lossy().ends_with(".mmi")))
            .unwrap_or(false);
        let mut detail = if cached {
            "FASTA; splice index already built".to_string()
        } else {
            "FASTA; the splice index is built once on first start (about a minute, several GB of RAM)".to_string()
        };
        if !have_fai {
            detail.push_str(" - no .fai next to it, so the provenance cannot fingerprint the reference");
        }
        v.push(check("Reference", if have_fai { Pass } else { Warn }, detail));
    }

    if let Some(total) = total_ram_gb() {
        let need = estimated_cycle_ram_gb(&cfg.reference);
        if need > 0.0 && total < need {
            v.push(check("Memory", Fail, format!("{total:.0} GB installed; each alignment cycle needs about {need:.0} GB")));
        } else if need > 0.0 && total < need * 1.5 {
            v.push(check("Memory", Warn, format!("{total:.0} GB installed; each alignment cycle needs about {need:.0} GB, little room for MinKNOW")));
        } else {
            v.push(check("Memory", Pass, format!("{total:.0} GB installed; a cycle needs about {need:.0} GB")));
        }
    }

    match fs::create_dir_all(&cfg.out_dir).and_then(|_| {
        let probe = cfg.out_dir.join(".write_test");
        fs::write(&probe, b"x")?;
        fs::remove_file(probe)
    }) {
        Ok(()) => {
            let free = free_disk_gb(&cfg.out_dir);
            match free {
                Some(g) if g < 20.0 => v.push(check("Output directory", Warn, format!("writable, but only {g:.0} GB free"))),
                Some(g) => v.push(check("Output directory", Pass, format!("writable, {g:.0} GB free"))),
                None => v.push(check("Output directory", Pass, "writable")),
            }
        }
        Err(e) => v.push(check("Output directory", Fail, format!("{}: {e}", cfg.out_dir.display()))),
    }
    if cfg.out_dir.join("ledger.jsonl").is_file() {
        let n = fs::read_to_string(cfg.out_dir.join("ledger.jsonl")).map(|t| t.lines().count()).unwrap_or(0);
        v.push(check("Previous run", Warn, format!("this output directory already holds {n} processed files; the run will resume, not start over")));
    }

    match &cfg.prelim {
        Some(p) => match p.check() {
            Ok(()) => v.push(check("nagini scripts", Pass, format!("{} (Python: {})", p.nagini_dir.display(), p.python))),
            Err(e) => v.push(check("nagini scripts", Fail, format!("{e:#}"))),
        },
        None => v.push(check("nagini scripts", Warn, "not configured: alignments are kept, but there will be no preliminary analysis")),
    }
    v
}

pub fn worst(checks: &[Check]) -> CheckStatus {
    if checks.iter().any(|c| c.status == CheckStatus::Fail) {
        CheckStatus::Fail
    } else if checks.iter().any(|c| c.status == CheckStatus::Warn) {
        CheckStatus::Warn
    } else {
        CheckStatus::Pass
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn status_of(checks: &[Check], name: &str) -> CheckStatus {
        checks.iter().find(|c| c.name == name).unwrap().status
    }

    #[test]
    fn reports_missing_input_reference_and_tools() {
        let tmp = tempfile::tempdir().unwrap();
        let mut cfg = EngineConfig::new(tmp.path().join("nope"), tmp.path().join("out"), tmp.path().join("ref.fa"));
        cfg.minimap2 = "definitely-not-installed-minimap2".into();
        let c = preflight(&cfg);
        assert_eq!(status_of(&c, "MinKNOW directory"), CheckStatus::Fail);
        assert_eq!(status_of(&c, "Reference"), CheckStatus::Fail);
        assert_eq!(status_of(&c, "minimap2"), CheckStatus::Fail);
        assert_eq!(status_of(&c, "Output directory"), CheckStatus::Pass);
        assert_eq!(status_of(&c, "nagini scripts"), CheckStatus::Warn);
        assert_eq!(worst(&c), CheckStatus::Fail);
    }

    #[test]
    fn counts_fastq_and_flags_a_resumed_output_directory() {
        let tmp = tempfile::tempdir().unwrap();
        let run = tmp.path().join("run/fastq_pass/barcode01");
        fs::create_dir_all(&run).unwrap();
        fs::write(run.join("a.fastq.gz"), "x").unwrap();
        fs::write(run.join("b.fastq"), "x").unwrap();
        fs::write(run.join("c.txt"), "x").unwrap();
        assert_eq!(count_fastq(&tmp.path().join("run"), 100), 2);
        assert_eq!(count_fastq(&tmp.path().join("run"), 1), 1);

        let fa = tmp.path().join("ref.fa");
        fs::write(&fa, ">c\nACGT\n").unwrap();
        let cfg = EngineConfig::new(tmp.path().join("run"), tmp.path().join("out"), fa);
        fs::create_dir_all(&cfg.out_dir).unwrap();
        fs::write(cfg.out_dir.join("ledger.jsonl"), "{}\n{}\n").unwrap();
        let c = preflight(&cfg);
        assert_eq!(status_of(&c, "MinKNOW directory"), CheckStatus::Pass);
        assert_eq!(status_of(&c, "Reference"), CheckStatus::Warn, "no .fai");
        assert!(c.iter().any(|x| x.name == "Previous run" && x.detail.contains("2 processed")));
    }

    #[test]
    fn memory_estimate_follows_the_measured_hg19_numbers() {
        let tmp = tempfile::tempdir().unwrap();
        let f = tmp.path().join("x.mmi");
        fs::write(&f, vec![0u8; 1_000_000]).unwrap();
        assert!((estimated_cycle_ram_gb(&f) - 0.0016).abs() < 1e-6);
    }
}
