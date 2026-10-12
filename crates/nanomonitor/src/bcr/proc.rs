//! Starting, stopping and attaching to the `nanostream bcr-abl watch` engine.
//!
//! The engine is a separate process that outlives the GUI: closing the window
//! does not interrupt alignment. The GUI only ever starts it, asks it to stop (by
//! creating `STOP` in its output directory) and reads what it writes.

use anyhow::{bail, Context, Result};
use nanobcr_core::prelim::PrelimConfig;
use nanobcr_core::EngineConfig;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};

#[derive(Debug, Clone)]
pub struct BcrSetup {
    pub input_dir: String,
    pub out_dir: String,
    pub reference: String,
    /// FASTA for the Branford k-mer check when `reference` is a `.mmi`.
    pub fasta: String,
    pub nagini_dir: String,
    pub python: String,
    pub build: String,
    pub qv: f64,
    pub min_len: usize,
    pub max_len: usize,
    pub threads: usize,
    pub pass_only: bool,
    pub snapshot_secs: u64,
    pub nanostream_bin: String,
}

impl Default for BcrSetup {
    fn default() -> Self {
        let defaults = EngineConfig::new(PathBuf::new(), PathBuf::new(), PathBuf::new());
        Self {
            input_dir: String::new(),
            out_dir: String::new(),
            reference: String::new(),
            fasta: String::new(),
            nagini_dir: String::new(),
            python: String::new(),
            build: "hg19".into(),
            qv: defaults.min_qv,
            min_len: defaults.min_len,
            max_len: defaults.max_len,
            threads: defaults.threads,
            pass_only: !defaults.include_fail,
            snapshot_secs: 300,
            nanostream_bin: String::new(),
        }
    }
}

fn nonempty(s: &str) -> Option<&str> {
    let t = s.trim();
    (!t.is_empty()).then_some(t)
}

impl BcrSetup {
    /// The configuration `nanostream bcr-abl watch` will run with, for preflight.
    pub fn engine_config(&self) -> EngineConfig {
        let mut cfg = EngineConfig::new(
            PathBuf::from(self.input_dir.trim()),
            PathBuf::from(self.out_dir.trim()),
            PathBuf::from(self.reference.trim()),
        );
        cfg.min_qv = self.qv;
        cfg.min_len = self.min_len;
        cfg.max_len = self.max_len;
        cfg.threads = self.threads;
        cfg.include_fail = !self.pass_only;
        cfg.snapshot_secs = self.snapshot_secs;
        if let Some(n) = nonempty(&self.nagini_dir) {
            let mut p = PrelimConfig::new(PathBuf::from(n));
            if let Some(py) = nonempty(&self.python) {
                p.python = py.to_string();
            }
            p.reference_build = self.build.clone();
            p.reference_fasta = nonempty(&self.fasta).map(PathBuf::from).or_else(|| {
                let r = PathBuf::from(self.reference.trim());
                (r.extension().map(|e| e != "mmi").unwrap_or(true)).then_some(r)
            });
            cfg.prelim = Some(p);
        }
        cfg
    }

    /// Arguments after the program name.
    pub fn args(&self) -> Vec<String> {
        let mut a: Vec<String> = vec![
            "bcr-abl".into(),
            "watch".into(),
            self.input_dir.trim().into(),
            "--out".into(),
            self.out_dir.trim().into(),
            "--reference".into(),
            self.reference.trim().into(),
            "--qv".into(),
            self.qv.to_string(),
            "--len".into(),
            format!("{}-{}", self.min_len, self.max_len),
            "--threads".into(),
            self.threads.to_string(),
            "--snapshot-secs".into(),
            self.snapshot_secs.to_string(),
        ];
        if self.pass_only {
            a.push("--pass-only".into());
        }
        if let Some(n) = nonempty(&self.nagini_dir) {
            a.extend(["--nagini".into(), n.into(), "--build".into(), self.build.clone()]);
            if let Some(py) = nonempty(&self.python) {
                a.extend(["--python".into(), py.into()]);
            }
            if let Some(f) = nonempty(&self.fasta) {
                a.extend(["--fasta".into(), f.into()]);
            }
        }
        a
    }
}

/// The configured binary, else a `nanostream` next to this executable (a cargo
/// build puts both in the same directory), else whatever is on PATH.
pub fn resolve_nanostream(configured: &str) -> String {
    if let Some(c) = nonempty(configured) {
        return c.to_string();
    }
    if let Some(dir) = std::env::current_exe().ok().and_then(|e| e.parent().map(Path::to_path_buf)) {
        let sibling = dir.join("nanostream");
        if sibling.is_file() {
            return sibling.to_string_lossy().into();
        }
    }
    "nanostream".into()
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

pub fn log_path(out_dir: &Path) -> PathBuf {
    out_dir.join("engine.log")
}

/// Last `max` lines of the engine log.
pub fn tail_log(out_dir: &Path, max: usize) -> Vec<String> {
    let text = fs::read_to_string(log_path(out_dir)).unwrap_or_default();
    let lines: Vec<&str> = text.lines().collect();
    lines[lines.len().saturating_sub(max)..].iter().map(|s| s.to_string()).collect()
}

#[derive(Default)]
pub struct EngineProc {
    child: Option<Child>,
}

impl EngineProc {
    /// Start the engine, logging to `<out>/engine.log`.
    pub fn start(&mut self, setup: &BcrSetup) -> Result<()> {
        if self.running() {
            bail!("the engine started from this window is still running");
        }
        let out = PathBuf::from(setup.out_dir.trim());
        fs::create_dir_all(&out).with_context(|| format!("creating {}", out.display()))?;
        let log = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(log_path(&out))?;
        let bin = resolve_nanostream(&setup.nanostream_bin);
        let child = Command::new(&bin)
            .args(setup.args())
            .env("RUST_LOG", "info")
            .stdin(Stdio::null())
            .stdout(log.try_clone()?)
            .stderr(log)
            .spawn()
            .with_context(|| format!("starting `{bin}` (set the nanostream path in the setup panel)"))?;
        self.child = Some(child);
        Ok(())
    }

    /// True while the child we started is alive.
    pub fn running(&mut self) -> bool {
        match self.child.as_mut() {
            Some(c) => match c.try_wait() {
                Ok(None) => true,
                _ => {
                    self.child = None;
                    false
                }
            },
            None => false,
        }
    }

    /// Exit status of a child that has just ended, if we started it.
    pub fn exit_message(&mut self) -> Option<String> {
        let c = self.child.as_mut()?;
        match c.try_wait() {
            Ok(Some(status)) => {
                self.child = None;
                Some(format!("engine exited: {status}"))
            }
            _ => None,
        }
    }

    /// Ask for a clean stop: the engine finishes its cycle, compacts and exits.
    pub fn request_stop(out_dir: &Path) -> Result<()> {
        fs::write(out_dir.join("STOP"), b"stop\n").context("writing STOP")
    }

    /// Last resort for a hung engine this window started.
    pub fn kill(&mut self) -> bool {
        let had_child = self.child.is_some();
        if let Some(c) = self.child.as_mut() {
            let _ = c.kill();
            let _ = c.wait();
        }
        self.child = None;
        had_child
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    fn setup() -> BcrSetup {
        BcrSetup {
            input_dir: " /run/minknow ".into(),
            out_dir: "/out".into(),
            reference: "/ref/hg19.fa".into(),
            nagini_dir: "/nagini".into(),
            python: "/nagini/.venv/bin/python".into(),
            ..Default::default()
        }
    }

    #[test]
    fn builds_the_watch_command_with_nagini_defaults() {
        let a = setup().args();
        assert_eq!(&a[..4], ["bcr-abl", "watch", "/run/minknow", "--out"]);
        let joined = a.join(" ");
        assert!(joined.contains("--qv 20") && joined.contains("--len 500-3000"), "{joined}");
        assert!(joined.contains("--nagini /nagini --build hg19 --python /nagini/.venv/bin/python"), "{joined}");
        assert!(!joined.contains("--pass-only"));

        let mut s = setup();
        s.pass_only = true;
        s.nagini_dir.clear();
        let joined = s.args().join(" ");
        assert!(joined.contains("--pass-only") && !joined.contains("--nagini"), "{joined}");
    }

    #[test]
    fn engine_config_uses_the_fasta_for_branford_unless_the_reference_is_an_index() {
        assert_eq!(
            setup().engine_config().prelim.unwrap().reference_fasta,
            Some(PathBuf::from("/ref/hg19.fa"))
        );
        let mut s = setup();
        s.reference = "/ref/hg19.splice.mmi".into();
        assert_eq!(s.engine_config().prelim.unwrap().reference_fasta, None);
        s.fasta = "/ref/hg19.fa".into();
        assert_eq!(s.engine_config().prelim.unwrap().reference_fasta, Some(PathBuf::from("/ref/hg19.fa")));
    }

    #[test]
    fn pid_alive_knows_itself_and_a_dead_pid() {
        assert!(pid_alive(std::process::id()));
        assert!(!pid_alive(4_000_000));
    }

    #[test]
    fn starts_a_process_with_those_args_and_captures_its_log() {
        let tmp = tempfile::tempdir().unwrap();
        let fake = tmp.path().join("fake-nanostream");
        fs::write(&fake, "#!/bin/sh\necho \"args: $@\"\necho \"rust_log=$RUST_LOG\" >&2\n").unwrap();
        fs::set_permissions(&fake, fs::Permissions::from_mode(0o755)).unwrap();

        let mut s = setup();
        s.out_dir = tmp.path().join("out").to_string_lossy().into();
        s.nanostream_bin = fake.to_string_lossy().into();
        let mut p = EngineProc::default();
        p.start(&s).unwrap();
        for _ in 0..100 {
            if !p.running() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(20));
        }
        let log = tail_log(&tmp.path().join("out"), 10).join("\n");
        assert!(log.contains("args: bcr-abl watch /run/minknow --out"), "{log}");
        assert!(log.contains("rust_log=info"), "{log}");

        EngineProc::request_stop(&tmp.path().join("out")).unwrap();
        assert!(tmp.path().join("out/STOP").exists());
    }

    #[test]
    fn a_missing_binary_is_a_clear_error() {
        let tmp = tempfile::tempdir().unwrap();
        let mut s = setup();
        s.out_dir = tmp.path().join("out").to_string_lossy().into();
        s.nanostream_bin = "/no/such/nanostream".into();
        let err = EngineProc::default().start(&s).unwrap_err();
        assert!(format!("{err:#}").contains("nanostream path"), "{err:#}");
    }
}
