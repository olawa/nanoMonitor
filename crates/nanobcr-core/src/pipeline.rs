//! Running nagini on the hand-off and following it.
//!
//! The run is `start_nagini.sh --input-type aligned_bam ...`, started in its own
//! process group with output in `nagini_run.log`. A small wrapper records the exit
//! code in `.nanomonitor_exit`, so a window (or a later session) that did not start
//! the run can still tell success from failure, and so a vanished process is
//! distinguishable from a finished one.

use crate::state::now_unix;
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    pub nagini_dir: PathBuf,
    pub analysis_dir: PathBuf,
    pub handoff_dir: PathBuf,
    pub barcode_map: PathBuf,
    /// nagini site: `workstation`, `local` or `marvin`.
    pub site: String,
    pub build: String,
    pub assay: String,
    pub library_type: String,
    pub cores: Option<usize>,
    pub enable_rindels: bool,
    pub sequence_id: Option<String>,
    /// Extra launcher arguments, passed through untouched.
    pub extra: Vec<String>,
}

impl PipelineConfig {
    pub fn new(nagini_dir: PathBuf, analysis_dir: PathBuf, handoff_dir: PathBuf, barcode_map: PathBuf) -> Self {
        Self {
            nagini_dir,
            analysis_dir,
            handoff_dir,
            barcode_map,
            site: "workstation".into(),
            build: "hg19".into(),
            assay: "bcr_abl".into(),
            library_type: "cdna".into(),
            cores: None,
            enable_rindels: false,
            sequence_id: None,
            extra: Vec::new(),
        }
    }

    pub fn launcher(&self) -> PathBuf {
        self.nagini_dir.join("pipeline_start_scripts/start_nagini.sh")
    }

    /// `bash start_nagini.sh ...`; the first element is the program.
    pub fn command(&self, dry_run: bool) -> Vec<String> {
        let s = |p: &Path| p.to_string_lossy().to_string();
        let mut c = vec![
            "bash".to_string(),
            s(&self.launcher()),
            "--inbox-path".into(),
            s(&self.handoff_dir),
            "--analysis-dir".into(),
            s(&self.analysis_dir),
            "--site".into(),
            self.site.clone(),
            "--input-type".into(),
            "aligned_bam".into(),
            "--assay".into(),
            self.assay.clone(),
            "--library-type".into(),
            self.library_type.clone(),
            "--reference-build".into(),
            self.build.clone(),
            "--barcode-map".into(),
            s(&self.barcode_map),
        ];
        if let Some(id) = &self.sequence_id {
            c.extend(["--sequence-id".into(), id.clone()]);
        }
        if let Some(n) = self.cores {
            c.extend(["--cores".into(), n.to_string()]);
        }
        if self.enable_rindels {
            c.push("--enable-rindels".into());
        }
        c.extend(self.extra.iter().cloned());
        if dry_run {
            c.push("--dry-run".into());
        }
        c
    }

    pub fn check(&self) -> Result<()> {
        if !self.launcher().is_file() {
            bail!("{} not found: is {} a nagini checkout?", self.launcher().display(), self.nagini_dir.display());
        }
        if !self.barcode_map.is_file() {
            bail!("barcode map {} does not exist", self.barcode_map.display());
        }
        if !self.handoff_dir.is_dir() {
            bail!("hand-off folder {} does not exist: prepare the hand-off first", self.handoff_dir.display());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunRecord {
    pub pid: u32,
    pub started_unix: u64,
    pub command: Vec<String>,
    pub analysis_dir: PathBuf,
    pub handoff_dir: PathBuf,
    pub barcode_map: PathBuf,
}

pub fn run_record_path(analysis_dir: &Path) -> PathBuf {
    analysis_dir.join("nanomonitor_run.json")
}
pub fn exit_file(analysis_dir: &Path) -> PathBuf {
    analysis_dir.join(".nanomonitor_exit")
}
pub fn log_file(analysis_dir: &Path) -> PathBuf {
    analysis_dir.join("nagini_run.log")
}

/// Output of `start_nagini.sh ... --dry-run`: whether snakemake can build the
/// workflow on this machine, without running anything.
#[derive(Debug, Clone)]
pub struct DryRun {
    pub ok: bool,
    pub tail: String,
}

pub fn dry_run(cfg: &PipelineConfig) -> Result<DryRun> {
    cfg.check()?;
    fs::create_dir_all(&cfg.analysis_dir)?;
    let c = cfg.command(true);
    let out = Command::new(&c[0]).args(&c[1..]).stdin(Stdio::null()).output().with_context(|| format!("starting `{}`", c[0]))?;
    let mut text = String::from_utf8_lossy(&out.stdout).to_string();
    text.push_str(&String::from_utf8_lossy(&out.stderr));
    let tail: Vec<&str> = text.lines().rev().take(25).collect();
    Ok(DryRun { ok: out.status.success(), tail: tail.into_iter().rev().collect::<Vec<_>>().join("\n") })
}

pub fn pid_alive(pid: u32) -> bool {
    Command::new("kill")
        .args(["-0", &pid.to_string()])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Start the pipeline in the background. Refuses if one is already running there.
pub fn start(cfg: &PipelineConfig) -> Result<RunRecord> {
    use std::os::unix::process::CommandExt;
    cfg.check()?;
    if let PipelineState::Running { .. } = status(&cfg.analysis_dir).state {
        bail!("a nagini run is already in progress in {}", cfg.analysis_dir.display());
    }
    fs::create_dir_all(&cfg.analysis_dir)?;
    let _ = fs::remove_file(exit_file(&cfg.analysis_dir));

    let command = cfg.command(false);
    let log = fs::OpenOptions::new().create(true).append(true).open(log_file(&cfg.analysis_dir))?;
    // `sh -c '"$@"; echo $? > "$0"' EXITFILE cmd...` records the exit code.
    let child = Command::new("sh")
        .arg("-c")
        .arg("\"$@\"; echo $? > \"$0\"")
        .arg(exit_file(&cfg.analysis_dir))
        .args(&command)
        .current_dir(&cfg.analysis_dir)
        .stdin(Stdio::null())
        .stdout(log.try_clone()?)
        .stderr(log)
        .process_group(0)
        .spawn()
        .context("starting the nagini launcher")?;
    let rec = RunRecord {
        pid: child.id(),
        started_unix: now_unix(),
        command,
        analysis_dir: cfg.analysis_dir.clone(),
        handoff_dir: cfg.handoff_dir.clone(),
        barcode_map: cfg.barcode_map.clone(),
    };
    crate::store::write_atomic(&run_record_path(&cfg.analysis_dir), &serde_json::to_vec_pretty(&rec)?)?;
    // The run outlives this call (and the window that made it); a reaper thread
    // only keeps it from lingering as a zombie while this process is alive.
    let mut child = child;
    std::thread::spawn(move || {
        let _ = child.wait();
    });
    Ok(rec)
}

/// Interrupt the run (SIGINT to the whole process group, as Ctrl-C would).
pub fn stop(analysis_dir: &Path) -> Result<()> {
    let rec = read_record(analysis_dir).context("no nagini run recorded here")?;
    let ok = Command::new("kill").args(["-INT", &format!("-{}", rec.pid)]).status()?.success();
    if !ok {
        bail!("could not signal process group {} (already finished?)", rec.pid);
    }
    Ok(())
}

pub fn read_record(analysis_dir: &Path) -> Option<RunRecord> {
    serde_json::from_slice(&fs::read(run_record_path(analysis_dir)).ok()?).ok()
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Progress {
    pub done: u32,
    pub total: u32,
    pub percent: u32,
    /// The last few rules snakemake announced (what is running or just ran).
    pub recent_rules: Vec<String>,
    pub failed_rules: Vec<String>,
}

/// Parse snakemake's log: `N of M steps (P%) done`, `rule X:` announcements and
/// `Error in rule X`. Works for the Snakemake 7 and 9 wording.
pub fn parse_progress(log: &str) -> Progress {
    let mut p = Progress::default();
    for line in log.lines() {
        let l = line.trim();
        if let Some(rest) = l.split_once(" steps (").and_then(|(head, tail)| {
            let (a, b) = head.split_once(" of ")?;
            let pct = tail.split_once("%)")?.0;
            Some((a.trim().parse::<u32>().ok()?, b.trim().parse::<u32>().ok()?, pct.trim().parse::<u32>().ok()?))
        }) {
            (p.done, p.total, p.percent) = rest;
        } else if let Some(name) = l.strip_prefix("rule ").and_then(|r| r.strip_suffix(':')) {
            if !name.contains(' ') {
                p.recent_rules.retain(|r| r != name);
                p.recent_rules.push(name.to_string());
                if p.recent_rules.len() > 4 {
                    p.recent_rules.remove(0);
                }
            }
        } else if let Some(name) = l.strip_prefix("Error in rule ").and_then(|r| r.strip_suffix(':')) {
            if !p.failed_rules.contains(&name.to_string()) {
                p.failed_rules.push(name.to_string());
            }
        }
    }
    p
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PipelineState {
    NotStarted,
    Running { pid: u32 },
    Succeeded,
    /// Exited non-zero, or exited 0 without producing the summary workbook.
    Failed { exit_code: Option<i32>, reason: String },
    /// No exit code was recorded and the process is gone (killed, machine rebooted).
    Lost,
}

#[derive(Debug, Clone)]
pub struct PipelineStatus {
    pub state: PipelineState,
    pub progress: Progress,
    pub started_unix: Option<u64>,
    pub log_tail: Vec<String>,
}

/// The newest snakemake log written since the run started.
fn newest_snakemake_log(analysis_dir: &Path, since_unix: u64) -> Option<String> {
    let dir = analysis_dir.join(".snakemake/log");
    let mut best: Option<(std::time::SystemTime, PathBuf)> = None;
    for e in fs::read_dir(dir).ok()?.flatten() {
        let m = e.metadata().ok()?.modified().ok()?;
        let unix = m.duration_since(std::time::UNIX_EPOCH).ok()?.as_secs();
        if unix + 5 >= since_unix && best.as_ref().map(|(t, _)| m > *t).unwrap_or(true) {
            best = Some((m, e.path()));
        }
    }
    fs::read_to_string(best?.1).ok()
}

pub fn summary_workbook(analysis_dir: &Path) -> Option<PathBuf> {
    let name = analysis_dir.file_name()?.to_string_lossy().to_string();
    [
        analysis_dir.join(format!("results/reports/{name}_summary_report.xlsx")),
        analysis_dir.join("results/reports/nagini_summary_report.xlsx"),
        analysis_dir.join("reports/nagini_summary_report.xlsx"),
    ]
    .into_iter()
    .find(|p| p.is_file())
}

pub fn status(analysis_dir: &Path) -> PipelineStatus {
    let Some(rec) = read_record(analysis_dir) else {
        return PipelineStatus { state: PipelineState::NotStarted, progress: Progress::default(), started_unix: None, log_tail: Vec::new() };
    };
    let log_text = fs::read_to_string(log_file(analysis_dir)).unwrap_or_default();
    let log_tail: Vec<String> = {
        let lines: Vec<&str> = log_text.lines().collect();
        lines[lines.len().saturating_sub(200)..].iter().map(|s| s.to_string()).collect()
    };
    let sm_log = newest_snakemake_log(analysis_dir, rec.started_unix).unwrap_or_default();
    let progress = parse_progress(&format!("{sm_log}\n{log_text}"));

    let state = match fs::read_to_string(exit_file(analysis_dir)).ok().and_then(|t| t.trim().parse::<i32>().ok()) {
        Some(0) if summary_workbook(analysis_dir).is_some() => PipelineState::Succeeded,
        Some(0) => PipelineState::Failed {
            exit_code: Some(0),
            reason: "the launcher exited normally but no summary workbook was produced".into(),
        },
        Some(code) => PipelineState::Failed {
            exit_code: Some(code),
            reason: if progress.failed_rules.is_empty() {
                format!("exit code {code}")
            } else {
                format!("exit code {code}; failed: {}", progress.failed_rules.join(", "))
            },
        },
        None if pid_alive(rec.pid) => PipelineState::Running { pid: rec.pid },
        None => PipelineState::Lost,
    };
    PipelineStatus { state, progress, started_unix: Some(rec.started_unix), log_tail }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    #[test]
    fn builds_the_launcher_command_for_an_aligned_bam_run() {
        let mut c = PipelineConfig::new("/n".into(), "/a/run1".into(), "/h".into(), "/h/map.csv".into());
        c.cores = Some(32);
        c.enable_rindels = true;
        c.sequence_id = Some("run1".into());
        c.extra = vec!["--force-inputs".into()];
        let cmd = c.command(false).join(" ");
        assert!(cmd.starts_with("bash /n/pipeline_start_scripts/start_nagini.sh --inbox-path /h --analysis-dir /a/run1 --site workstation"), "{cmd}");
        for want in ["--input-type aligned_bam", "--assay bcr_abl", "--library-type cdna", "--reference-build hg19", "--barcode-map /h/map.csv", "--cores 32", "--enable-rindels", "--sequence-id run1", "--force-inputs"] {
            assert!(cmd.contains(want), "{want} missing in {cmd}");
        }
        assert!(!cmd.contains("--dry-run"));
        assert!(c.command(true).last().unwrap() == "--dry-run");
    }

    #[test]
    fn parses_progress_and_failures_from_snakemake_wording() {
        let log = "Building DAG\nrule a:\n    input: x\n[t]\nFinished job 1.\n1 of 8 steps (12%) done\nrule b:\nrule a:\n\
                   Error in rule bcr_abl_fusion_counts:\n3 of 8 steps (38%) done\nExiting because a job execution failed.";
        let p = parse_progress(log);
        assert_eq!((p.done, p.total, p.percent), (3, 8, 38));
        assert_eq!(p.recent_rules, ["b", "a"]);
        assert_eq!(p.failed_rules, ["bcr_abl_fusion_counts"]);
        assert_eq!(parse_progress("nothing here"), Progress::default());
    }

    /// A launcher that behaves like a pipeline: prints progress, then either
    /// writes the summary workbook or fails.
    fn fake_nagini(root: &Path, outcome: &str) -> PathBuf {
        let dir = root.join("pipeline_start_scripts");
        fs::create_dir_all(&dir).unwrap();
        let script = dir.join("start_nagini.sh");
        let body = format!(
            "#!/bin/bash\nwhile [[ $# -gt 0 ]]; do case \"$1\" in --analysis-dir) A=\"$2\";; --dry-run) D=1;; esac; shift; done\n\
             if [[ -n \"$D\" ]]; then echo 'dry run ok'; exit 0; fi\n\
             mkdir -p \"$A/.snakemake/log\" \"$A/results/reports\"\n\
             echo 'rule one:' | tee -a \"$A/.snakemake/log/2026-01-01T000000.snakemake.log\"\n\
             echo '1 of 2 steps (50%) done' | tee -a \"$A/.snakemake/log/2026-01-01T000000.snakemake.log\"\n\
             case '{outcome}' in\n  ok) echo '2 of 2 steps (100%) done'; touch \"$A/results/reports/nagini_summary_report.xlsx\";;\n  \
             fail) echo 'Error in rule two:'; exit 3;;\n  hollow) echo 'done';;\nesac\n"
        );
        fs::write(&script, body).unwrap();
        fs::set_permissions(&script, fs::Permissions::from_mode(0o755)).unwrap();
        root.to_path_buf()
    }

    fn wait_finished(a: &Path) -> PipelineStatus {
        for _ in 0..200 {
            let s = status(a);
            if !matches!(s.state, PipelineState::Running { .. }) {
                return s;
            }
            std::thread::sleep(std::time::Duration::from_millis(50));
        }
        panic!("run did not finish");
    }

    fn setup(outcome: &str) -> (tempfile::TempDir, PipelineConfig) {
        let d = tempfile::tempdir().unwrap();
        let nagini = fake_nagini(&d.path().join("nagini"), outcome);
        let handoff = d.path().join("handoff");
        fs::create_dir_all(&handoff).unwrap();
        let map = handoff.join("barcode_map.csv");
        fs::write(&map, "barcode,sample\nbarcode01,A\n").unwrap();
        let cfg = PipelineConfig::new(nagini, d.path().join("analysis"), handoff, map);
        (d, cfg)
    }

    #[test]
    fn a_successful_run_is_followed_to_the_summary_workbook() {
        let (_d, cfg) = setup("ok");
        assert!(dry_run(&cfg).unwrap().ok);
        assert_eq!(status(&cfg.analysis_dir).state, PipelineState::NotStarted);
        let rec = start(&cfg).unwrap();
        assert!(rec.command.contains(&"aligned_bam".to_string()));
        let s = wait_finished(&cfg.analysis_dir);
        assert_eq!(s.state, PipelineState::Succeeded);
        assert_eq!((s.progress.done, s.progress.total), (2, 2));
        assert!(s.log_tail.iter().any(|l| l.contains("100%")));
    }

    #[test]
    fn a_failing_run_reports_the_exit_code_and_rule() {
        let (_d, cfg) = setup("fail");
        start(&cfg).unwrap();
        match wait_finished(&cfg.analysis_dir).state {
            PipelineState::Failed { exit_code, reason } => {
                assert_eq!(exit_code, Some(3));
                assert!(reason.contains("two"), "{reason}");
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn exiting_zero_without_the_workbook_is_not_success() {
        let (_d, cfg) = setup("hollow");
        start(&cfg).unwrap();
        assert!(matches!(wait_finished(&cfg.analysis_dir).state, PipelineState::Failed { exit_code: Some(0), .. }));
    }

    #[test]
    fn a_vanished_process_is_lost_and_a_second_start_is_refused_while_running() {
        let (_d, cfg) = setup("ok");
        fs::create_dir_all(&cfg.analysis_dir).unwrap();
        let rec = RunRecord {
            pid: 4_000_000,
            started_unix: now_unix(),
            command: vec![],
            analysis_dir: cfg.analysis_dir.clone(),
            handoff_dir: cfg.handoff_dir.clone(),
            barcode_map: cfg.barcode_map.clone(),
        };
        fs::write(run_record_path(&cfg.analysis_dir), serde_json::to_vec(&rec).unwrap()).unwrap();
        assert_eq!(status(&cfg.analysis_dir).state, PipelineState::Lost);

        let rec = RunRecord { pid: std::process::id(), ..rec };
        fs::write(run_record_path(&cfg.analysis_dir), serde_json::to_vec(&rec).unwrap()).unwrap();
        assert!(matches!(status(&cfg.analysis_dir).state, PipelineState::Running { .. }));
        assert!(start(&cfg).unwrap_err().to_string().contains("already in progress"));
    }

    #[test]
    fn missing_launcher_or_map_is_reported_before_anything_runs() {
        let (_d, mut cfg) = setup("ok");
        cfg.barcode_map = "/no/map.csv".into();
        assert!(cfg.check().unwrap_err().to_string().contains("barcode map"));
        cfg.nagini_dir = "/no/nagini".into();
        assert!(cfg.check().unwrap_err().to_string().contains("nagini checkout"));
    }
}
