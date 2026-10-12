//! One alignment cycle: pipe filtered reads through minimap2, split the SAM by
//! barcode in-process (counting as we go) and write one BAM segment per barcode.

use crate::config::EngineConfig;
use crate::reads::Read;
use anyhow::{anyhow, bail, Context, Result};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::{BufRead, BufReader, Read as _, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};
use std::thread;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignCounts {
    /// SAM records, including supplementary and secondary.
    pub records: u64,
    pub primary: u64,
    pub primary_mapped: u64,
    pub primary_unmapped: u64,
    pub supplementary: u64,
    pub secondary: u64,
}

impl AlignCounts {
    pub fn add(&mut self, o: &AlignCounts) {
        self.records += o.records;
        self.primary += o.primary;
        self.primary_mapped += o.primary_mapped;
        self.primary_unmapped += o.primary_unmapped;
        self.supplementary += o.supplementary;
        self.secondary += o.secondary;
    }
}

#[derive(Debug)]
pub struct AlignResult {
    pub counts: BTreeMap<String, AlignCounts>,
    /// (barcode, finished segment BAM)
    pub segments: Vec<(String, PathBuf)>,
}

/// Run a tool to completion and return an error carrying its stderr on failure.
pub fn run_checked(program: &str, args: &[String]) -> Result<()> {
    let out = Command::new(program)
        .args(args)
        .stdin(Stdio::null())
        .output()
        .with_context(|| format!("starting `{program}`"))?;
    if !out.status.success() {
        bail!(
            "`{program} {}` failed ({}): {}",
            args.join(" "),
            out.status,
            tail(&String::from_utf8_lossy(&out.stderr))
        );
    }
    Ok(())
}

fn tail(s: &str) -> String {
    let s = s.trim();
    if s.len() > 1500 {
        format!("...{}", &s[s.len() - 1500..])
    } else {
        s.to_string()
    }
}

pub fn tool_version(program: &str) -> Result<String> {
    let out = Command::new(program)
        .arg("--version")
        .output()
        .with_context(|| format!("`{program}` not found on PATH; install it or pass its path"))?;
    let text = String::from_utf8_lossy(&out.stdout);
    Ok(text.lines().next().unwrap_or("").trim().to_string())
}

/// The minimap2 index to align against. A `.mmi` is used as given. A FASTA is
/// indexed once with the same preset (a prebuilt index silently keeps the `-k/-w`
/// it was built with, so it must be built with `-x splice` to align exactly like
/// aligning against the FASTA) and cached under `out_dir/index`.
pub fn ensure_index(cfg: &EngineConfig) -> Result<PathBuf> {
    if cfg.reference.extension().map(|e| e == "mmi").unwrap_or(false) {
        if !cfg.reference.exists() {
            bail!("index {} does not exist", cfg.reference.display());
        }
        return Ok(cfg.reference.clone());
    }
    let fasta_mtime = fs::metadata(&cfg.reference)
        .with_context(|| format!("reference {}", cfg.reference.display()))?
        .modified()?;
    let dir = cfg.out_dir.join("index");
    fs::create_dir_all(&dir)?;
    let stem = cfg
        .reference
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "reference".into());
    let mmi = dir.join(format!("{stem}.{}.mmi", cfg.preset));
    if let Ok(meta) = fs::metadata(&mmi) {
        if meta.modified()? >= fasta_mtime {
            return Ok(mmi);
        }
    }
    log::info!("building minimap2 -x {} index {} (once)", cfg.preset, mmi.display());
    let tmp = dir.join(format!(".{stem}.{}.mmi.tmp", cfg.preset));
    run_checked(
        &cfg.minimap2,
        &[
            "-x".into(),
            cfg.preset.clone(),
            "-t".into(),
            cfg.threads.to_string(),
            "-d".into(),
            tmp.to_string_lossy().into(),
            cfg.reference.to_string_lossy().into(),
        ],
    )?;
    fs::rename(&tmp, &mmi)?;
    Ok(mmi)
}

pub fn minimap2_args(cfg: &EngineConfig, index: &Path) -> Vec<String> {
    let mut args = vec![
        "-a".to_string(),
        "-x".into(),
        cfg.preset.clone(),
        "-t".into(),
        cfg.threads.to_string(),
    ];
    args.extend(cfg.minimap2_extra.iter().cloned());
    args.push(index.to_string_lossy().into());
    args.push("-".into());
    args
}

struct SegmentWriter {
    child: Child,
    stdin: Option<ChildStdin>,
    tmp: PathBuf,
    dest: PathBuf,
}

/// Align `reads` and write one BAM segment per barcode into `inbox_for(barcode)`,
/// named `{tag}.bam`. Nothing is left behind on failure.
pub fn align_reads(
    cfg: &EngineConfig,
    index: &Path,
    reads: &[Read],
    inbox_for: &dyn Fn(&str) -> PathBuf,
    tag: &str,
) -> Result<AlignResult> {
    let mut result = AlignResult {
        counts: BTreeMap::new(),
        segments: Vec::new(),
    };
    if reads.is_empty() {
        return Ok(result);
    }

    let barcode_of: HashMap<&[u8], &str> = reads
        .iter()
        .map(|r| (r.name.as_slice(), r.barcode.as_str()))
        .collect();

    let mut mm = Command::new(&cfg.minimap2)
        .args(minimap2_args(cfg, index))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .with_context(|| format!("starting `{}`", cfg.minimap2))?;
    let mm_stdin = mm.stdin.take().unwrap();
    let mm_stdout = mm.stdout.take().unwrap();
    let mut mm_stderr = mm.stderr.take().unwrap();

    let mut writers: HashMap<String, SegmentWriter> = HashMap::new();
    let mut header: Vec<u8> = Vec::new();

    let outcome: Result<()> = thread::scope(|scope| {
        let feeder = scope.spawn(move || -> std::io::Result<()> {
            let mut w = std::io::BufWriter::with_capacity(1 << 20, mm_stdin);
            for r in reads {
                w.write_all(b"@")?;
                w.write_all(&r.header)?;
                w.write_all(b"\n")?;
                w.write_all(&r.seq)?;
                w.write_all(b"\n+\n")?;
                w.write_all(&r.qual)?;
                w.write_all(b"\n")?;
            }
            w.flush()
        });
        let stderr_reader = scope.spawn(move || {
            let mut s = String::new();
            let _ = mm_stderr.read_to_string(&mut s);
            s
        });

        let mut rd = BufReader::with_capacity(1 << 20, mm_stdout);
        let mut line = Vec::new();
        let read_loop: Result<()> = (|| {
            loop {
                line.clear();
                if rd.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                if line.first() == Some(&b'@') {
                    header.extend_from_slice(&line);
                    continue;
                }
                let mut fields = line.splitn(3, |&b| b == b'\t');
                let qname = fields.next().unwrap_or(&[]);
                let flag: u16 = fields
                    .next()
                    .and_then(|f| std::str::from_utf8(f).ok())
                    .and_then(|f| f.parse().ok())
                    .ok_or_else(|| anyhow!("malformed SAM line from minimap2"))?;
                let barcode = *barcode_of.get(qname).ok_or_else(|| {
                    anyhow!(
                        "minimap2 returned unknown read {}",
                        String::from_utf8_lossy(qname)
                    )
                })?;

                let c = result.counts.entry(barcode.to_string()).or_default();
                c.records += 1;
                if flag & 0x800 != 0 {
                    c.supplementary += 1;
                } else if flag & 0x100 != 0 {
                    c.secondary += 1;
                } else {
                    c.primary += 1;
                    if flag & 0x4 != 0 {
                        c.primary_unmapped += 1;
                    } else {
                        c.primary_mapped += 1;
                    }
                }

                if !writers.contains_key(barcode) {
                    let inbox = inbox_for(barcode);
                    fs::create_dir_all(&inbox)?;
                    let tmp = inbox.join(format!(".{tag}.tmp"));
                    let dest = inbox.join(format!("{tag}.bam"));
                    let mut child = Command::new(&cfg.samtools)
                        .args(["view", "-b", "-o"])
                        .arg(&tmp)
                        .arg("-")
                        .stdin(Stdio::piped())
                        .stdout(Stdio::null())
                        .stderr(Stdio::piped())
                        .spawn()
                        .with_context(|| format!("starting `{}`", cfg.samtools))?;
                    let mut stdin = child.stdin.take().unwrap();
                    stdin.write_all(&header)?;
                    writers.insert(
                        barcode.to_string(),
                        SegmentWriter {
                            child,
                            stdin: Some(stdin),
                            tmp,
                            dest,
                        },
                    );
                }
                writers
                    .get_mut(barcode)
                    .unwrap()
                    .stdin
                    .as_mut()
                    .unwrap()
                    .write_all(&line)?;
            }
            Ok(())
        })();

        // Make sure minimap2 is not left blocked on a full pipe if we bailed early.
        drop(rd);
        let feed = feeder.join().expect("feeder thread panicked");
        let mm_err = stderr_reader.join().unwrap_or_default();
        let status = mm.wait()?;
        if !status.success() {
            bail!("minimap2 failed ({status}): {}", tail(&mm_err));
        }
        read_loop?;
        feed.context("feeding reads to minimap2")?;
        Ok(())
    });

    // Close and reap every samtools writer, then publish or discard.
    let mut failure = outcome.err();
    for (barcode, mut w) in writers.drain() {
        drop(w.stdin.take());
        let ok = match w.child.wait_with_output() {
            Ok(out) if out.status.success() => true,
            Ok(out) => {
                failure.get_or_insert_with(|| {
                    anyhow!(
                        "samtools view failed for {barcode}: {}",
                        tail(&String::from_utf8_lossy(&out.stderr))
                    )
                });
                false
            }
            Err(e) => {
                failure.get_or_insert_with(|| anyhow!("samtools view: {e}"));
                false
            }
        };
        if ok && failure.is_none() {
            fs::rename(&w.tmp, &w.dest)?;
            result.segments.push((barcode, w.dest));
        } else {
            let _ = fs::remove_file(&w.tmp);
        }
    }
    if let Some(e) = failure {
        // A failed cycle must not leave half its segments behind: the files are not
        // ledgered and will be redone, which would double-count the finished ones.
        for (_, seg) in &result.segments {
            let _ = fs::remove_file(seg);
        }
        return Err(e);
    }
    Ok(result)
}
