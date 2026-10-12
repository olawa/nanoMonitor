//! Preliminary BCR-ABL1 analysis on a snapshot BAM.
//!
//! The analysis is nagini's own: the same three scripts, called with the same
//! arguments as its Snakemake rules (`bcr_abl.smk`), reading nagini's resources
//! and thresholds. Nothing is re-implemented here, so a preliminary result and the
//! final pipeline result can only differ because the data differ, not the method.
//! `bcr_abl_exon_coverage.py` is deliberately left out: it loops over every exon
//! base in Python and is too slow to repeat every few minutes.

use crate::store::write_atomic;
use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread::JoinHandle;
use std::time::Instant;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrelimConfig {
    /// nagini repository root (contains `workflow/` and `config/`).
    pub nagini_dir: PathBuf,
    /// Python with pysam. Defaults to nagini's `.venv`, else `python3`.
    pub python: String,
    pub reference_build: String,
    /// Reference FASTA for the Branford caller's k-mer context check (optional).
    pub reference_fasta: Option<PathBuf>,
    pub threads: usize,
    /// Overrides for nagini's `bcr_abl:` settings.
    pub overrides: BTreeMap<String, String>,
}

impl PrelimConfig {
    pub fn new(nagini_dir: PathBuf) -> Self {
        let venv = nagini_dir.join(".venv/bin/python");
        let python = if venv.exists() {
            venv.to_string_lossy().into()
        } else {
            "python3".to_string()
        };
        Self {
            nagini_dir,
            python,
            reference_build: "hg19".into(),
            reference_fasta: None,
            threads: 4,
            overrides: BTreeMap::new(),
        }
    }

    fn script(&self, name: &str) -> PathBuf {
        self.nagini_dir.join("workflow/scripts").join(name)
    }

    fn resource(&self, name: &str) -> PathBuf {
        self.nagini_dir.join("workflow/resources/bcr_abl").join(name)
    }

    fn branford_bed(&self) -> PathBuf {
        let name = if self.reference_build.contains("19") {
            "branford_positions_hg19.bed"
        } else {
            "branford_positions_hg38.bed"
        };
        self.resource(name)
    }

    /// Fail early, with the fix, if nagini or its Python environment is unusable.
    pub fn check(&self) -> Result<()> {
        for s in ["bcr_abl_fusion_counter.py", "bcr_abl_branford_caller.py", "bcr_abl_phaser.py"] {
            if !self.script(s).is_file() {
                bail!("{} not found: is --nagini {} a nagini checkout?", self.script(s).display(), self.nagini_dir.display());
            }
        }
        for f in ["bcr_abl_features.json"] {
            if !self.resource(f).is_file() {
                bail!("{} not found", self.resource(f).display());
            }
        }
        if !self.branford_bed().is_file() {
            bail!("{} not found", self.branford_bed().display());
        }
        let out = Command::new(&self.python)
            .args(["-c", "import pysam"])
            .output()
            .with_context(|| format!("starting `{}`", self.python))?;
        if !out.status.success() {
            bail!(
                "`{}` cannot import pysam ({}); point --python at nagini's environment",
                self.python,
                String::from_utf8_lossy(&out.stderr).trim()
            );
        }
        Ok(())
    }

    /// nagini's `bcr_abl:` settings: its config.yaml, then our overrides.
    pub fn settings(&self) -> BTreeMap<String, String> {
        let mut m = read_bcr_abl_section(&self.nagini_dir.join("config/config.yaml"));
        m.extend(self.overrides.clone());
        m
    }
}

/// The flat `bcr_abl:` map from nagini's config.yaml. Deliberately not a YAML
/// parser: the section is `key: value` lines under an unindented `bcr_abl:`.
pub fn read_bcr_abl_section(path: &Path) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    let Ok(text) = fs::read_to_string(path) else {
        return out;
    };
    let mut inside = false;
    for line in text.lines() {
        if line.starts_with("bcr_abl:") {
            inside = true;
            continue;
        }
        if !inside {
            continue;
        }
        if !line.starts_with(' ') && !line.trim().is_empty() && !line.trim_start().starts_with('#') {
            break;
        }
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if let Some((k, v)) = t.split_once(':') {
            let v = v.split('#').next().unwrap_or("").trim().trim_matches('"').to_string();
            if !v.is_empty() {
                out.insert(k.trim().to_string(), v);
            }
        }
    }
    out
}

// ---------------------------------------------------------------- results

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct IsoformRow {
    pub isoform: String,
    pub protein_form: String,
    pub read_count: u64,
    pub fraction_percent: f64,
    pub status: String,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct FusionSummary {
    pub total_fusion_reads: u64,
    pub primary_isoform: Option<String>,
    pub isoforms: Vec<IsoformRow>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct MutationRow {
    pub mutation: String,
    pub aa_change: String,
    pub chrom: String,
    pub pos: u64,
    pub depth: u64,
    pub ad_alt: u64,
    pub vaf_percent: f64,
    pub status: String,
    pub tki: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct BranfordSummary {
    pub sites_total: u64,
    pub sites_at_min_depth: u64,
    pub min_depth_setting: u64,
    pub min_site_depth: u64,
    pub median_site_depth: u64,
    pub mutations: Vec<MutationRow>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct PhasePair {
    pub mutation_1: String,
    pub mutation_2: String,
    pub phase_state: String,
    pub cis_fraction_percent: f64,
    pub total_spanning_reads: u64,
}

/// Is there enough data to read the result? A heuristic for the UI, with its
/// inputs recorded so it can be judged: a PRIMARY fusion call and nearly every
/// kinase-domain site covered at nagini's minimum calling depth.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Sufficiency {
    pub fusion_called: bool,
    pub site_coverage_fraction: f64,
    pub required_site_coverage_fraction: f64,
    pub enough: bool,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PrelimSummary {
    pub schema: u32,
    pub barcode: String,
    pub generation: u64,
    pub reads_passed: u64,
    pub written_unix: u64,
    pub duration_secs: f64,
    /// `ok`, `partial` (some stage failed) or `failed`.
    pub status: String,
    pub errors: Vec<String>,
    pub fusion: Option<FusionSummary>,
    pub branford: Option<BranfordSummary>,
    pub phasing: Option<Vec<PhasePair>>,
    pub sufficiency: Option<Sufficiency>,
    /// Where the full nagini script outputs are.
    pub dir: String,
}

fn unix_now() -> u64 {
    crate::state::now_unix()
}

/// Rows of a nagini TSV: `#` comment lines skipped, first other line is the header.
pub fn read_tsv(path: &Path) -> Result<Vec<BTreeMap<String, String>>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut header: Option<Vec<&str>> = None;
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match &header {
            None => header = Some(fields),
            Some(h) => rows.push(
                h.iter()
                    .zip(fields.iter())
                    .map(|(k, v)| (k.to_string(), v.to_string()))
                    .collect(),
            ),
        }
    }
    Ok(rows)
}

fn num<T: std::str::FromStr + Default>(row: &BTreeMap<String, String>, key: &str) -> T {
    row.get(key).and_then(|v| v.parse().ok()).unwrap_or_default()
}

fn text(row: &BTreeMap<String, String>, key: &str) -> String {
    row.get(key).cloned().unwrap_or_default()
}

pub fn parse_fusion(json_path: &Path) -> Result<FusionSummary> {
    let v: serde_json::Value = serde_json::from_slice(&fs::read(json_path)?)?;
    let isoforms = v["isoforms"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|r| IsoformRow {
                    isoform: r["fusion_isoform"].as_str().unwrap_or("").into(),
                    protein_form: r["protein_form"].as_str().unwrap_or("").into(),
                    read_count: r["read_count"].as_u64().unwrap_or(0),
                    fraction_percent: r["fraction_percent"].as_f64().unwrap_or(0.0),
                    status: r["status"].as_str().unwrap_or("").into(),
                })
                .collect()
        })
        .unwrap_or_default();
    Ok(FusionSummary {
        total_fusion_reads: v["total_fusion_reads"].as_u64().unwrap_or(0),
        primary_isoform: v["primary_isoform"].as_str().filter(|s| !s.is_empty()).map(String::from),
        isoforms,
    })
}

const TKIS: [&str; 6] = ["imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib", "asciminib"];

pub fn parse_branford(all_sites: &Path, mutations: &Path, min_depth: u64) -> Result<BranfordSummary> {
    // One row per (site, alt); depth is a property of the site.
    let mut depth_by_site: BTreeMap<(String, u64), u64> = BTreeMap::new();
    for r in read_tsv(all_sites)? {
        depth_by_site.insert((text(&r, "chrom"), num(&r, "pos")), num(&r, "depth"));
    }
    let mut depths: Vec<u64> = depth_by_site.values().copied().collect();
    depths.sort_unstable();
    let muts = read_tsv(mutations)?
        .into_iter()
        .map(|r| MutationRow {
            mutation: text(&r, "mutation"),
            aa_change: text(&r, "aa_change"),
            chrom: text(&r, "chrom"),
            pos: num(&r, "pos"),
            depth: num(&r, "depth"),
            ad_alt: num(&r, "ad_alt"),
            vaf_percent: num(&r, "vaf_percent"),
            status: text(&r, "status"),
            tki: TKIS.iter().map(|t| (t.to_string(), text(&r, t))).collect(),
        })
        .collect();
    Ok(BranfordSummary {
        sites_total: depths.len() as u64,
        sites_at_min_depth: depths.iter().filter(|&&d| d >= min_depth).count() as u64,
        min_depth_setting: min_depth,
        min_site_depth: depths.first().copied().unwrap_or(0),
        median_site_depth: depths.get(depths.len() / 2).copied().unwrap_or(0),
        mutations: muts,
    })
}

pub fn parse_phasing(json_path: &Path) -> Result<Vec<PhasePair>> {
    let v: serde_json::Value = serde_json::from_slice(&fs::read(json_path)?)?;
    Ok(v["phased_pairs"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|p| PhasePair {
                    mutation_1: p["mutation_1"].as_str().unwrap_or("").into(),
                    mutation_2: p["mutation_2"].as_str().unwrap_or("").into(),
                    phase_state: p["phase_state"].as_str().unwrap_or("").into(),
                    cis_fraction_percent: p["cis_fraction_percent"].as_f64().unwrap_or(0.0),
                    total_spanning_reads: p["total_spanning_reads"].as_u64().unwrap_or(0),
                })
                .collect()
        })
        .unwrap_or_default())
}

pub const REQUIRED_SITE_COVERAGE: f64 = 0.95;

pub fn sufficiency(fusion: Option<&FusionSummary>, branford: Option<&BranfordSummary>) -> Sufficiency {
    let fusion_called = fusion.map(|f| f.primary_isoform.is_some()).unwrap_or(false);
    let site_coverage_fraction = branford
        .filter(|b| b.sites_total > 0)
        .map(|b| b.sites_at_min_depth as f64 / b.sites_total as f64)
        .unwrap_or(0.0);
    Sufficiency {
        fusion_called,
        site_coverage_fraction,
        required_site_coverage_fraction: REQUIRED_SITE_COVERAGE,
        enough: fusion_called && site_coverage_fraction >= REQUIRED_SITE_COVERAGE,
    }
}

// ---------------------------------------------------------------- running

fn run_script(cfg: &PrelimConfig, script: &str, args: &[String], log: &Path) -> Result<()> {
    let logfile = fs::File::create(log)?;
    let status = Command::new(&cfg.python)
        .arg(cfg.script(script))
        .args(args)
        .stdin(Stdio::null())
        .stdout(logfile.try_clone()?)
        .stderr(logfile)
        .status()
        .with_context(|| format!("starting `{}`", cfg.python))?;
    if !status.success() {
        let tail: String = fs::read_to_string(log)
            .unwrap_or_default()
            .lines()
            .rev()
            .take(6)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
            .join(" | ");
        bail!("{script} failed ({status}): {tail}");
    }
    Ok(())
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into()
}

/// Run the three scripts on `bam` (coordinate-sorted, indexed) and write
/// `summary.json` into `out_dir`. Stages that fail are reported, not fatal.
pub fn run_prelim(
    cfg: &PrelimConfig,
    barcode: &str,
    bam: &Path,
    generation: u64,
    reads_passed: u64,
    out_dir: &Path,
) -> Result<PrelimSummary> {
    let started = Instant::now();
    fs::create_dir_all(out_dir)?;
    let set = cfg.settings();
    let get = |k: &str, d: &str| set.get(k).cloned().unwrap_or_else(|| d.to_string());

    // The scripts read a BAM that compaction may replace while they run, so work
    // from a private (hard-linked, else copied) copy of this generation.
    let input_dir = out_dir.join("input");
    let _ = fs::remove_dir_all(&input_dir);
    fs::create_dir_all(&input_dir)?;
    let my_bam = input_dir.join(format!("{barcode}.gen{generation}.bam"));
    let bai = PathBuf::from(format!("{}.bai", bam.display()));
    for (from, to) in [(bam.to_path_buf(), my_bam.clone()), (bai, PathBuf::from(format!("{}.bai", my_bam.display())))] {
        if fs::hard_link(&from, &to).is_err() {
            fs::copy(&from, &to).with_context(|| format!("copying {}", from.display()))?;
        }
    }

    let threads = cfg.threads.max(1).to_string();
    let features = cfg.resource("bcr_abl_features.json");
    let mut errors = Vec::new();

    // 1. fusion isoforms
    let fusion_tsv = out_dir.join("fusion_summary.tsv");
    let fusion_json = out_dir.join("fusion_summary.json");
    let r = run_script(
        cfg,
        "bcr_abl_fusion_counter.py",
        &[
            "--bam".into(), s(&my_bam), "--sample".into(), barcode.into(),
            "--genome-build".into(), cfg.reference_build.clone(),
            "--features-json".into(), s(&features),
            "--output".into(), s(&fusion_tsv), "--json-output".into(), s(&fusion_json),
            "--arriba-output".into(), s(&out_dir.join("arriba_fusions.tsv")),
            "--tolerance".into(), get("tolerance", "150"),
            "--min-mapq".into(), get("min_mapq", "1"),
            "--threads".into(), threads.clone(),
        ],
        &out_dir.join("fusion.log"),
    );
    let fusion = match r.and_then(|_| parse_fusion(&fusion_json)) {
        Ok(f) => Some(f),
        Err(e) => {
            errors.push(format!("fusion: {e:#}"));
            None
        }
    };

    // 2. kinase-domain mutations
    let all_sites = out_dir.join("branford_all_sites.tsv");
    let muts = out_dir.join("branford_mutations.tsv");
    let mut args: Vec<String> = vec![
        "--bam".into(), s(&my_bam), "--bed".into(), s(&cfg.branford_bed()),
        "--sample".into(), barcode.into(),
        "--output-all".into(), s(&all_sites), "--output-mutations".into(), s(&muts),
        "--output-vcf".into(), s(&out_dir.join("branford_mutations.vcf")),
        "--features-json".into(), s(&features),
    ];
    if let Some(fa) = &cfg.reference_fasta {
        args.extend(["--reference-fasta".into(), s(fa)]);
    }
    for (flag, key, default) in [
        ("--min-af", "min_af", "0.01"), ("--subclonal-af", "subclonal_af", "0.002"),
        ("--min-depth", "min_depth", "50"), ("--min-alt-reads", "min_alt_reads", "5"),
        ("--max-depth", "max_depth", "1000000"), ("--min-base-quality", "min_base_quality", "12"),
        ("--min-snr", "min_snr", "3.0"), ("--min-minor-strand-ratio", "min_minor_strand_ratio", "0.05"),
        ("--min-mean-baseq", "min_mean_baseq", "14.0"), ("--min-kmer-rate", "min_kmer_rate", "0.60"),
        ("--max-sor", "max_sor", "3.0"), ("--kmer-window", "kmer_window", "3"),
    ] {
        args.extend([flag.into(), get(key, default)]);
    }
    if let Some(kv) = set.get("known_variants") {
        args.extend(["--known-variants".into(), kv.clone()]);
    }
    args.extend(["--threads".into(), threads.clone()]);
    let min_depth: u64 = get("min_depth", "50").parse().unwrap_or(50);
    let r = run_script(cfg, "bcr_abl_branford_caller.py", &args, &out_dir.join("branford.log"));
    let branford = match r.and_then(|_| parse_branford(&all_sites, &muts, min_depth)) {
        Ok(b) => Some(b),
        Err(e) => {
            errors.push(format!("branford: {e:#}"));
            None
        }
    };

    // 3. phasing of detected mutations (needs the mutations TSV)
    let phasing = if branford.is_some() {
        let ph_tsv = out_dir.join("phasing.tsv");
        let ph_json = out_dir.join("phasing.json");
        let r = run_script(
            cfg,
            "bcr_abl_phaser.py",
            &[
                "--bam".into(), s(&my_bam), "--mutations-tsv".into(), s(&muts),
                "--sample".into(), barcode.into(),
                "--output-tsv".into(), s(&ph_tsv), "--output-json".into(), s(&ph_json),
                "--min-spanning-reads".into(), get("min_spanning_reads", "10"),
                "--min-mapq".into(), get("min_mapq", "1"),
                "--threads".into(), threads,
            ],
            &out_dir.join("phasing.log"),
        );
        match r.and_then(|_| parse_phasing(&ph_json)) {
            Ok(p) => Some(p),
            Err(e) => {
                errors.push(format!("phasing: {e:#}"));
                None
            }
        }
    } else {
        None
    };

    let _ = fs::remove_dir_all(&input_dir);
    let status = if errors.is_empty() {
        "ok"
    } else if fusion.is_some() || branford.is_some() {
        "partial"
    } else {
        "failed"
    };
    let suff = sufficiency(fusion.as_ref(), branford.as_ref());
    let summary = PrelimSummary {
        schema: 1,
        barcode: barcode.into(),
        generation,
        reads_passed,
        written_unix: unix_now(),
        duration_secs: started.elapsed().as_secs_f64(),
        status: status.into(),
        errors,
        fusion,
        branford,
        phasing,
        sufficiency: Some(suff),
        dir: s(out_dir),
    };
    write_atomic(&out_dir.join("summary.json"), &serde_json::to_vec_pretty(&summary)?)?;

    // One line per snapshot, so the GUI can plot how the call firms up with depth.
    let history = serde_json::json!({
        "unix": summary.written_unix,
        "generation": generation,
        "reads_passed": reads_passed,
        "status": summary.status,
        "primary_isoform": summary.fusion.as_ref().and_then(|f| f.primary_isoform.clone()),
        "fusion_reads": summary.fusion.as_ref().map(|f| f.total_fusion_reads),
        "mutations": summary.branford.as_ref().map(|b| b.mutations.iter()
            .map(|m| serde_json::json!({"mutation": m.mutation, "vaf_percent": m.vaf_percent, "depth": m.depth, "status": m.status}))
            .collect::<Vec<_>>()),
        "enough": summary.sufficiency.as_ref().map(|x| x.enough),
    });
    let mut h = fs::OpenOptions::new().create(true).append(true).open(out_dir.join("history.jsonl"))?;
    writeln!(h, "{history}")?;
    Ok(summary)
}

/// Analyse the current snapshot of every barcode under an engine output directory
/// (the engine need not be running). `only` restricts it to one barcode.
pub fn analyze_out_dir(cfg: &PrelimConfig, out_dir: &Path, only: Option<&str>) -> Result<Vec<PrelimSummary>> {
    cfg.check()?;
    let state: Option<crate::state::State> = fs::read(out_dir.join("state.json"))
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok());
    let mut out = Vec::new();
    for entry in fs::read_dir(out_dir)?.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if only.map(|o| o != name).unwrap_or(false) || !entry.path().join("inbox").is_dir() {
            continue;
        }
        let store = crate::store::Store::new(&entry.path());
        let Some(bam) = store.current_bam() else {
            log::warn!("{name}: no compacted BAM yet (run `bcr-abl snapshot` first)");
            continue;
        };
        let reads = state.as_ref().and_then(|s| s.barcodes.get(&name)).map(|b| b.reads_passed).unwrap_or(0);
        out.push(run_prelim(cfg, &name, &bam, store.manifest().generation, reads, &entry.path().join("prelim"))?);
    }
    Ok(out)
}

// ---------------------------------------------------------------- background worker

pub struct PrelimJob {
    pub barcode: String,
    pub bam: PathBuf,
    pub generation: u64,
    pub reads_passed: u64,
}

#[derive(Debug, Clone)]
pub struct PrelimDone {
    pub barcode: String,
    pub generation: u64,
    pub status: String,
    pub duration_secs: f64,
    pub error: Option<String>,
}

/// Runs preliminary analyses one at a time off the engine's thread, so a slow
/// Python run never delays an alignment cycle. If several snapshots of one
/// barcode queue up, only the newest is analysed.
pub struct PrelimRunner {
    tx: Option<Sender<PrelimJob>>,
    pub done: Receiver<PrelimDone>,
    handle: Option<JoinHandle<()>>,
}

impl PrelimRunner {
    pub fn spawn(cfg: PrelimConfig, out_root: PathBuf) -> Self {
        let (tx, rx) = mpsc::channel::<PrelimJob>();
        let (done_tx, done) = mpsc::channel::<PrelimDone>();
        let handle = std::thread::spawn(move || {
            while let Ok(first) = rx.recv() {
                let mut latest: BTreeMap<String, PrelimJob> = BTreeMap::new();
                latest.insert(first.barcode.clone(), first);
                while let Ok(j) = rx.try_recv() {
                    latest.insert(j.barcode.clone(), j);
                }
                for (barcode, job) in latest {
                    let dir = out_root.join(&barcode).join("prelim");
                    let res = run_prelim(&cfg, &barcode, &job.bam, job.generation, job.reads_passed, &dir);
                    let msg = match res {
                        Ok(s) => PrelimDone {
                            barcode: barcode.clone(),
                            generation: job.generation,
                            status: s.status.clone(),
                            duration_secs: s.duration_secs,
                            error: s.errors.first().cloned(),
                        },
                        Err(e) => PrelimDone {
                            barcode: barcode.clone(),
                            generation: job.generation,
                            status: "failed".into(),
                            duration_secs: 0.0,
                            error: Some(format!("{e:#}")),
                        },
                    };
                    let _ = done_tx.send(msg);
                }
            }
        });
        Self { tx: Some(tx), done, handle: Some(handle) }
    }

    pub fn submit(&self, job: PrelimJob) -> Result<()> {
        self.tx
            .as_ref()
            .ok_or_else(|| anyhow!("prelim runner already shut down"))?
            .send(job)
            .map_err(|_| anyhow!("prelim worker has stopped"))
    }

    /// Finish queued work and stop the worker.
    pub fn shutdown(&mut self) {
        self.tx.take();
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

impl Drop for PrelimRunner {
    fn drop(&mut self) {
        self.shutdown();
    }
}

/// A stand-in for a nagini checkout: the three scripts write canned results and
/// record their argv, so the plumbing can be tested without pysam or real data.
#[cfg(test)]
pub(crate) mod stub {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    pub fn nagini(root: &Path, fail_fusion: bool) -> PrelimConfig {
        let scripts = root.join("workflow/scripts");
        let res = root.join("workflow/resources/bcr_abl");
        fs::create_dir_all(&scripts).unwrap();
        fs::create_dir_all(&res).unwrap();
        fs::create_dir_all(root.join("config")).unwrap();
        fs::write(res.join("bcr_abl_features.json"), "{}").unwrap();
        fs::write(res.join("branford_positions_hg19.bed"), "").unwrap();
        fs::write(root.join("config/config.yaml"), "bcr_abl:\n  min_af: 0.02\n  min_depth: 100\n").unwrap();
        let common = "import sys,json\nargs=sys.argv[1:]\nopen(sys.argv[0]+'.argv','w').write(' '.join(args))\ndef a(k): return args[args.index(k)+1]\n";
        fs::write(
            scripts.join("bcr_abl_fusion_counter.py"),
            format!("{common}{}\njson.dump({{'total_fusion_reads':30,'primary_isoform':'e14a2','isoforms':[{{'fusion_isoform':'e14a2','protein_form':'p210','read_count':30,'fraction_percent':100.0,'status':'PRIMARY'}}]}},open(a('--json-output'),'w'))\nopen(a('--output'),'w').write('x')\n",
                if fail_fusion { "sys.stderr.write('boom\\n'); sys.exit(3)" } else { "" }),
        ).unwrap();
        fs::write(
            scripts.join("bcr_abl_branford_caller.py"),
            format!("{common}h='sample\\tchrom\\tpos\\tmutation\\taa_change\\tdepth\\tad_alt\\tvaf_percent\\tstatus\\timatinib\\n'\nopen(a('--output-all'),'w').write('# c\\nsample\\tchrom\\tpos\\tdepth\\nx\\tchr9\\t1\\t500\\nx\\tchr9\\t2\\t60\\n')\nopen(a('--output-mutations'),'w').write('#c\\n'+h+'x\\tchr9\\t1\\tT315I\\tp.T315I\\t500\\t20\\t4.0\\tDETECTED\\tResistant\\n')\n"),
        ).unwrap();
        fs::write(
            scripts.join("bcr_abl_phaser.py"),
            format!("{common}json.dump({{'phased_pairs':[]}},open(a('--output-json'),'w'))\nopen(a('--output-tsv'),'w').write('x')\n"),
        ).unwrap();
        // `python` that accepts the `import pysam` probe and otherwise runs python3.
        let py = root.join("fakepy");
        fs::write(&py, "#!/bin/sh\nif [ \"$1\" = \"-c\" ]; then exit 0; fi\nexec python3 \"$@\"\n").unwrap();
        fs::set_permissions(&py, fs::Permissions::from_mode(0o755)).unwrap();
        let mut cfg = PrelimConfig::new(root.to_path_buf());
        cfg.python = py.to_string_lossy().into();
        cfg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_flat_bcr_abl_section_of_nagini_config() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("config.yaml");
        fs::write(
            &p,
            "other:\n  min_af: 0.5\n\nbcr_abl:\n  tolerance: 150\n  min_af: 0.01 # one percent\n  # comment\n  known_variants: \"/x/y.bed\"\n\namplicon:\n  min_af: 9\n",
        )
        .unwrap();
        let m = read_bcr_abl_section(&p);
        assert_eq!(m.get("tolerance").map(String::as_str), Some("150"));
        assert_eq!(m.get("min_af").map(String::as_str), Some("0.01"));
        assert_eq!(m.get("known_variants").map(String::as_str), Some("/x/y.bed"));
        assert_eq!(m.len(), 3);
    }

    #[test]
    fn parses_branford_and_decides_sufficiency() {
        let dir = tempfile::tempdir().unwrap();
        let all = dir.path().join("all.tsv");
        let mut f = fs::File::create(&all).unwrap();
        writeln!(f, "# sample=x\nsample\tchrom\tpos\tdepth").unwrap();
        // Two alts at one site share a depth; one shallow site.
        for (pos, depth) in [(100, 500), (100, 500), (200, 400), (300, 10)] {
            writeln!(f, "x\tchr9\t{pos}\t{depth}").unwrap();
        }
        let muts = dir.path().join("m.tsv");
        let mut f = fs::File::create(&muts).unwrap();
        writeln!(f, "sample\tchrom\tpos\tmutation\taa_change\tdepth\tad_alt\tvaf_percent\tstatus\timatinib").unwrap();
        writeln!(f, "x\tchr9\t100\tT315I\tp.Thr315Ile\t500\t12\t2.4\tDETECTED\tResistant").unwrap();
        drop(f);

        let b = parse_branford(&all, &muts, 50).unwrap();
        assert_eq!((b.sites_total, b.sites_at_min_depth, b.min_site_depth), (3, 2, 10));
        assert_eq!(b.mutations[0].mutation, "T315I");
        assert_eq!(b.mutations[0].tki["imatinib"], "Resistant");

        let fusion = FusionSummary { primary_isoform: Some("e14a2".into()), ..Default::default() };
        let s = sufficiency(Some(&fusion), Some(&b));
        assert!(s.fusion_called && !s.enough, "2/3 sites covered is below 95 %");
        assert!(!sufficiency(None, Some(&b)).enough);
    }

    fn fake_bam(dir: &Path) -> PathBuf {
        let bam = dir.join("compacted.000001.bam");
        fs::write(&bam, "bam").unwrap();
        fs::write(dir.join("compacted.000001.bam.bai"), "bai").unwrap();
        bam
    }

    #[test]
    fn runs_the_three_scripts_with_nagini_settings_and_writes_the_summary() {
        let tmp = tempfile::tempdir().unwrap();
        let cfg = stub::nagini(&tmp.path().join("nagini"), false);
        let bam = fake_bam(tmp.path());
        let out = tmp.path().join("prelim");

        let s = run_prelim(&cfg, "barcode07", &bam, 1, 321, &out).unwrap();
        assert_eq!(s.status, "ok");
        assert_eq!(s.fusion.as_ref().unwrap().primary_isoform.as_deref(), Some("e14a2"));
        let b = s.branford.as_ref().unwrap();
        assert_eq!((b.sites_total, b.sites_at_min_depth, b.min_depth_setting), (2, 1, 100));
        assert_eq!(b.mutations[0].mutation, "T315I");
        assert!(s.phasing.as_ref().unwrap().is_empty());
        assert!(out.join("summary.json").exists());
        assert!(!out.join("input").exists(), "private BAM copy is removed");

        // The scripts got nagini's own thresholds and resources, not our defaults.
        let argv = fs::read_to_string(tmp.path().join("nagini/workflow/scripts/bcr_abl_branford_caller.py.argv")).unwrap();
        assert!(argv.contains("--min-af 0.02") && argv.contains("--min-depth 100"), "{argv}");
        assert!(argv.contains("branford_positions_hg19.bed") && argv.contains("--sample barcode07"), "{argv}");
        let fusion_argv = fs::read_to_string(tmp.path().join("nagini/workflow/scripts/bcr_abl_fusion_counter.py.argv")).unwrap();
        assert!(fusion_argv.contains("--genome-build hg19") && fusion_argv.contains("--tolerance 150"), "{fusion_argv}");

        let history = fs::read_to_string(out.join("history.jsonl")).unwrap();
        let line: serde_json::Value = serde_json::from_str(history.lines().next().unwrap()).unwrap();
        assert_eq!(line["reads_passed"], 321);
        assert_eq!(line["primary_isoform"], "e14a2");
        assert_eq!(line["mutations"][0]["mutation"], "T315I");
    }

    #[test]
    fn a_failing_stage_is_reported_and_the_others_still_run() {
        let tmp = tempfile::tempdir().unwrap();
        let cfg = stub::nagini(&tmp.path().join("nagini"), true);
        let bam = fake_bam(tmp.path());
        let s = run_prelim(&cfg, "barcode07", &bam, 2, 5, &tmp.path().join("prelim")).unwrap();
        assert_eq!(s.status, "partial");
        assert!(s.fusion.is_none());
        assert!(s.errors[0].starts_with("fusion:") && s.errors[0].contains("boom"), "{:?}", s.errors);
        assert_eq!(s.branford.as_ref().unwrap().mutations.len(), 1);
        assert!(!s.sufficiency.as_ref().unwrap().enough);
    }

    #[test]
    fn check_reports_a_missing_nagini_checkout() {
        let cfg = PrelimConfig::new(PathBuf::from("/nonexistent/nagini"));
        assert!(cfg.check().unwrap_err().to_string().contains("nagini checkout"));
    }
}
