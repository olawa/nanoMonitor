use crate::align::{align_reads, ensure_index, tool_version, AlignCounts};
use crate::config::EngineConfig;
use crate::ledger::{Ledger, LedgerEntry};
use crate::prelim::{PrelimJob, PrelimRunner};
use crate::reads::{filter_settings, load_fastq, LoadedFile, Read, ReadCounts};
use crate::state::{now_unix, write_provenance, State};
use crate::store::Store;
use crate::watcher::{run_finished, Watcher};
use anyhow::{bail, Context, Result};
use nanofilter_core::bam::FilterSettings;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

/// Files that fail to parse (still being written) are retried with a growing delay.
const RETRY_BASE_SECS: u64 = 5;
const RETRY_MAX_SECS: u64 = 60;
/// In `--once` mode an unreadable file is given up on after this many attempts.
const ONCE_MAX_ATTEMPTS: u32 = 5;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct StepReport {
    pub files: usize,
    pub reads_seen: u64,
    pub reads_passed: u64,
    pub aligned: bool,
    pub compacted: bool,
    /// Candidates not yet ready, or being retried.
    pub waiting: usize,
}

pub struct Engine {
    pub cfg: EngineConfig,
    index: PathBuf,
    settings: FilterSettings,
    ledger: Ledger,
    watcher: Watcher,
    stores: BTreeMap<String, Store>,
    pub state: State,
    retries: HashMap<String, (u32, Instant)>,
    last_cycle: Option<Instant>,
    last_compact: Instant,
    cycle_counter: u64,
    prelim: Option<PrelimRunner>,
}

impl Engine {
    pub fn new(mut cfg: EngineConfig) -> Result<Self> {
        if !cfg.input_dir.is_dir() {
            bail!("input directory {} does not exist", cfg.input_dir.display());
        }
        std::fs::create_dir_all(&cfg.out_dir)?;
        // Absolute paths, because state.json and provenance.json are read by other
        // processes (the GUI, nagini) that do not share this working directory.
        cfg.input_dir = std::fs::canonicalize(&cfg.input_dir)?;
        cfg.out_dir = std::fs::canonicalize(&cfg.out_dir)?;
        cfg.reference = std::fs::canonicalize(&cfg.reference)
            .with_context(|| format!("reference {}", cfg.reference.display()))?;
        tool_version(&cfg.minimap2)?;
        tool_version(&cfg.samtools)?;
        let index = ensure_index(&cfg)?;

        let ledger = Ledger::open(&cfg.out_dir.join("ledger.jsonl"))?;
        let state = State::resume(&cfg, &cfg.out_dir.join("state.json"));
        write_provenance(&cfg, &index, &cfg.out_dir.join("provenance.json"))?;

        let mut stores = BTreeMap::new();
        for barcode in state.barcodes.keys() {
            let s = Store::new(&cfg.out_dir.join(barcode));
            s.recover()?;
            stores.insert(barcode.clone(), s);
        }
        // Stores that exist on disk but not in the (possibly lost) state file.
        if let Ok(rd) = std::fs::read_dir(&cfg.out_dir) {
            for e in rd.flatten() {
                let name = e.file_name().to_string_lossy().to_string();
                if e.path().join("inbox").is_dir() && !stores.contains_key(&name) {
                    let s = Store::new(&e.path());
                    s.recover()?;
                    stores.insert(name, s);
                }
            }
        }

        let prelim = match &cfg.prelim {
            Some(p) => {
                p.check()?;
                Some(PrelimRunner::spawn(p.clone(), cfg.out_dir.clone()))
            }
            None => None,
        };
        let settings = filter_settings(cfg.min_qv, cfg.min_len, cfg.max_len);
        let watcher = Watcher::new(&cfg.input_dir, cfg.include_fail, cfg.stable_secs);
        Ok(Self {
            index,
            settings,
            ledger,
            watcher,
            stores,
            state,
            retries: HashMap::new(),
            last_cycle: None,
            last_compact: Instant::now(),
            cycle_counter: 0,
            prelim,
            cfg,
        })
    }

    fn state_path(&self) -> PathBuf {
        self.cfg.out_dir.join("state.json")
    }

    fn store_for(&mut self, barcode: &str) -> &Store {
        let dir = self.cfg.out_dir.join(barcode);
        self.stores
            .entry(barcode.to_string())
            .or_insert_with(|| Store::new(&dir))
    }

    pub fn save_state(&mut self) -> Result<()> {
        let path = self.state_path();
        self.state.write(&path)
    }

    /// One scan and, if there is new data and the cycle interval allows, one
    /// filter -> align -> segment cycle (plus a compaction when due).
    pub fn step(&mut self, force_cycle: bool) -> Result<StepReport> {
        let mut report = StepReport::default();
        let ready = self.watcher.scan();
        let mut waiting = self.watcher.waiting();

        let now = Instant::now();
        let todo: Vec<_> = ready
            .into_iter()
            .filter(|c| !self.ledger.is_done(&c.rel_path))
            .filter(|c| match self.retries.get(&c.rel_path) {
                // Given up on in a one-shot run: no longer counts as outstanding.
                Some((attempts, _)) if self.cfg.once && *attempts >= ONCE_MAX_ATTEMPTS => false,
                Some((_, next)) if *next > now => {
                    waiting += 1;
                    false
                }
                _ => true,
            })
            .collect();

        let cycle_due = force_cycle
            || self
                .last_cycle
                .map(|t| t.elapsed() >= Duration::from_secs(self.cfg.min_cycle_secs))
                .unwrap_or(true);

        if !todo.is_empty() && cycle_due {
            let mut batch: Vec<(String, u64, LoadedFile)> = Vec::new();
            let mut total_reads = 0usize;
            for cand in todo {
                if total_reads >= self.cfg.max_reads_per_cycle {
                    waiting += 1;
                    continue;
                }
                match load_fastq(&cand.path, &self.settings) {
                    Ok(loaded) => {
                        self.retries.remove(&cand.rel_path);
                        total_reads += loaded.reads.len();
                        batch.push((cand.rel_path, cand.size, loaded));
                    }
                    Err(e) => {
                        let attempts = self.retries.get(&cand.rel_path).map(|r| r.0).unwrap_or(0) + 1;
                        let delay = (RETRY_BASE_SECS * attempts as u64).min(RETRY_MAX_SECS);
                        log::warn!(
                            "{} not readable yet (attempt {attempts}): {e:#}",
                            cand.rel_path
                        );
                        self.retries.insert(
                            cand.rel_path.clone(),
                            (attempts, Instant::now() + Duration::from_secs(delay)),
                        );
                        if !(self.cfg.once && attempts >= ONCE_MAX_ATTEMPTS) {
                            waiting += 1;
                        }
                    }
                }
            }

            if !batch.is_empty() {
                self.run_cycle(batch, &mut report)?;
            }
        } else if !todo.is_empty() {
            waiting += todo.len();
        }

        let compact_due = self.last_compact.elapsed() >= Duration::from_secs(self.cfg.snapshot_secs);
        if compact_due {
            report.compacted = self.compact_all()?;
        }
        report.waiting = waiting;
        self.collect_prelim();
        self.refresh_pending();
        self.state.status = if report.aligned { "running" } else { "idle" }.into();
        self.save_state()?;
        Ok(report)
    }

    fn run_cycle(&mut self, mut batch: Vec<(String, u64, LoadedFile)>, report: &mut StepReport) -> Result<()> {
        let mut all_reads: Vec<Read> = Vec::new();
        for (_, _, loaded) in batch.iter_mut() {
            all_reads.append(&mut loaded.reads);
        }

        let started = Instant::now();
        self.cycle_counter += 1;
        let tag = format!("seg_{}_{:04}", now_unix(), self.cycle_counter);
        let out_dir = self.cfg.out_dir.clone();
        let result = align_reads(
            &self.cfg,
            &self.index,
            &all_reads,
            &|barcode| out_dir.join(barcode).join("inbox"),
            &tag,
        )
        .context("alignment cycle")?;
        self.last_cycle = Some(Instant::now());

        // Only now, with the segments safely on disk, do the files count as done.
        let now = now_unix();
        let mut read_totals: BTreeMap<String, ReadCounts> = BTreeMap::new();
        for (rel, size, loaded) in &batch {
            let seen: u64 = loaded.counts.values().map(|c| c.seen).sum();
            let passed: u64 = loaded.counts.values().map(|c| c.passed).sum();
            self.ledger.record(LedgerEntry {
                path: rel.clone(),
                size: *size,
                reads: seen,
                passed,
                ts_unix: now,
            })?;
            report.reads_seen += seen;
            report.reads_passed += passed;
            for (barcode, c) in &loaded.counts {
                read_totals.entry(barcode.clone()).or_default().add(c);
                let b = self.state.barcode(barcode);
                b.files += 1;
                b.last_file = Some(rel.clone());
            }
        }
        for (barcode, c) in &read_totals {
            let b = self.state.barcode(barcode);
            b.add_reads(c);
            b.last_update_unix = now;
        }
        let mut align_totals: BTreeMap<String, AlignCounts> = BTreeMap::new();
        for (barcode, c) in &result.counts {
            align_totals.entry(barcode.clone()).or_default().add(c);
        }
        for (barcode, c) in &align_totals {
            self.state.barcode(barcode).add_align(c);
        }
        let passed_in_cycle: u64 = report.reads_passed;
        let primary_in_cycle: u64 = result.counts.values().map(|c| c.primary).sum();
        if primary_in_cycle != passed_in_cycle {
            log::warn!(
                "cycle {}: {passed_in_cycle} reads passed the filter but minimap2 returned {primary_in_cycle} primary records",
                self.cycle_counter
            );
        }
        for (barcode, _) in &result.segments {
            self.store_for(barcode);
        }

        self.state.files_done = self.ledger.len() as u64;
        self.state.cycles += 1;
        log::info!(
            "cycle {}: {} files, {} reads, {} passed the filter, {} primary records, aligned in {:.1}s",
            self.cycle_counter,
            batch.len(),
            report.reads_seen,
            report.reads_passed,
            primary_in_cycle,
            started.elapsed().as_secs_f64()
        );
        report.files = batch.len();
        report.aligned = true;
        Ok(())
    }

    /// Fold every barcode's inbox into its sorted, indexed BAM.
    pub fn compact_all(&mut self) -> Result<bool> {
        let mut any = false;
        let threads = self.cfg.threads;
        let samtools = self.cfg.samtools.clone();
        let names: Vec<String> = self.stores.keys().cloned().collect();
        for name in names {
            let store = &self.stores[&name];
            if let Some(c) = store.compact(&samtools, threads)? {
                any = true;
                let b = self.state.barcode(&name);
                b.compacted_bam = Some(c.bam.to_string_lossy().into());
                b.compacted_generation = c.generation;
                b.compacted_unix = Some(now_unix());
                log::info!("{name}: compacted {} segments into {}", c.segments_folded, c.bam.display());
                if let Some(runner) = &self.prelim {
                    let reads_passed = self.state.barcode(&name).reads_passed;
                    let job = PrelimJob {
                        barcode: name.clone(),
                        bam: c.bam.clone(),
                        generation: c.generation,
                        reads_passed,
                    };
                    let b = self.state.barcode(&name);
                    b.prelim_status = Some("queued".into());
                    b.prelim_summary = Some(
                        self.cfg.out_dir.join(&name).join("prelim/summary.json").to_string_lossy().into(),
                    );
                    if let Err(e) = runner.submit(job) {
                        log::warn!("{name}: preliminary analysis not queued: {e:#}");
                    }
                }
            }
        }
        self.last_compact = Instant::now();
        self.refresh_pending();
        Ok(any)
    }

    /// Record finished preliminary analyses in the state.
    fn collect_prelim(&mut self) {
        let mut finished = Vec::new();
        if let Some(r) = &self.prelim {
            while let Ok(d) = r.done.try_recv() {
                finished.push(d);
            }
        }
        for d in finished {
            let b = self.state.barcode(&d.barcode);
            b.prelim_generation = Some(d.generation);
            b.prelim_status = Some(d.status);
            b.prelim_unix = Some(now_unix());
            b.prelim_error = d.error;
        }
    }

    fn refresh_pending(&mut self) {
        let pending: Vec<(String, usize)> = self
            .stores
            .iter()
            .map(|(n, s)| (n.clone(), s.pending().len()))
            .collect();
        for (n, p) in pending {
            self.state.barcode(&n).segments_pending = p;
        }
    }

    /// Run until stopped (or, with `once`, until the directory is drained).
    pub fn run(&mut self, stop: &AtomicBool) -> Result<()> {
        let stop_file = self.cfg.out_dir.join("STOP");
        let _ = std::fs::remove_file(&stop_file);
        self.save_state()?;
        loop {
            if stop.load(Ordering::Relaxed) || stop_file.exists() {
                log::info!("stop requested");
                break;
            }
            let report = match self.step(self.cfg.once) {
                Ok(r) => r,
                Err(e) => {
                    log::error!("{e:#}");
                    self.state.status = "error".into();
                    self.state.last_error = Some(format!("{e:#}"));
                    self.save_state()?;
                    std::thread::sleep(Duration::from_secs(10));
                    continue;
                }
            };
            if report.aligned {
                self.state.last_error = None;
            }
            if self.cfg.once && report.waiting == 0 && !report.aligned {
                break;
            }
            if !self.cfg.once && report.waiting == 0 && !report.aligned && run_finished(&self.cfg.input_dir) {
                log::info!("MinKNOW wrote its final summary and nothing is left to process");
                break;
            }
            std::thread::sleep(Duration::from_secs(self.cfg.poll_secs.max(1)));
        }
        self.compact_all()?;
        // Let the last preliminary analysis finish: it is the result of the run.
        if let Some(r) = self.prelim.as_mut() {
            r.shutdown();
        }
        self.collect_prelim();
        self.state.status = "finished".into();
        self.save_state()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::collections::BTreeSet;
    use std::io::Write;
    use std::path::Path;
    use std::process::Command;

    fn tools_available() -> bool {
        let ok = |p: &str| Command::new(p).arg("--version").output().map(|o| o.status.success()).unwrap_or(false);
        let have = ok("minimap2") && ok("samtools");
        if !have {
            eprintln!("skipping: minimap2/samtools not on PATH");
        }
        have
    }

    /// Small deterministic RNG so the tests need no extra crate.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn random_dna(rng: &mut Rng, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[rng.below(4)]).collect()
    }

    fn mutate(rng: &mut Rng, s: &[u8], per_mille: usize) -> Vec<u8> {
        s.iter()
            .map(|&b| if rng.below(1000) < per_mille { b"ACGT"[rng.below(4)] } else { b })
            .collect()
    }

    struct Fixture {
        _tmp: tempfile::TempDir,
        root: PathBuf,
        reference: PathBuf,
        genome: Vec<Vec<u8>>,
        rng: Rng,
        counter: usize,
    }

    impl Fixture {
        fn new() -> Self {
            let tmp = tempfile::tempdir().unwrap();
            let root = tmp.path().to_path_buf();
            let mut rng = Rng(0x9E3779B97F4A7C15);
            let genome = vec![random_dna(&mut rng, 40_000), random_dna(&mut rng, 30_000)];
            let reference = root.join("ref.fa");
            let mut f = std::fs::File::create(&reference).unwrap();
            for (i, g) in genome.iter().enumerate() {
                writeln!(f, ">chr{}", i + 1).unwrap();
                f.write_all(g).unwrap();
                f.write_all(b"\n").unwrap();
            }
            Fixture { _tmp: tmp, root, reference, genome, rng, counter: 0 }
        }

        /// A FASTQ chunk the way MinKNOW writes it. Returns the names expected to pass
        /// the default filter (qs >= 20, 500-3000 bp).
        fn write_chunk(&mut self, rel: &str, n: usize, gz: bool) -> Vec<String> {
            let mut text = Vec::new();
            let mut passing = Vec::new();
            for _ in 0..n {
                self.counter += 1;
                let name = format!("read-{:05}", self.counter);
                let contig = self.rng.below(2);
                let len = 300 + self.rng.below(3500); // some too short, some too long
                let start = self.rng.below(self.genome[contig].len() - len);
                let seq = mutate(&mut self.rng, &self.genome[contig][start..start + len], 30);
                let qs = if self.rng.below(5) == 0 { 12.0 } else { 18.0 + self.rng.below(100) as f32 / 10.0 };
                let barcode = Path::new(rel).components().find_map(|c| {
                    let s = c.as_os_str().to_string_lossy().to_string();
                    s.starts_with("barcode").then_some(s)
                }).unwrap_or_default();
                writeln!(
                    text,
                    "@{name} qs:f:{qs:.4} ch:i:{} st:Z:2026-08-28T10:15:56.552323+00:00 SM:Z:{barcode}",
                    1 + self.rng.below(500)
                ).unwrap();
                text.extend_from_slice(&seq);
                text.extend_from_slice(b"\n+\n");
                text.extend(std::iter::repeat_n(b'5', len));
                text.push(b'\n');
                if qs >= 20.0 && (500..=3000).contains(&len) {
                    passing.push(name);
                }
            }
            let path = self.root.join("run").join(rel);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            if gz {
                let mut e = GzEncoder::new(Vec::new(), Compression::fast());
                e.write_all(&text).unwrap();
                std::fs::write(&path, e.finish().unwrap()).unwrap();
            } else {
                std::fs::write(&path, &text).unwrap();
            }
            passing
        }

        fn config(&self) -> EngineConfig {
            let mut cfg = EngineConfig::new(self.root.join("run"), self.root.join("out"), self.reference.clone());
            cfg.threads = 2;
            cfg.stable_secs = 0;
            cfg.poll_secs = 1;
            cfg.min_cycle_secs = 0;
            cfg.snapshot_secs = 0;
            cfg.once = true;
            cfg
        }
    }

    fn bam_names(bam: &Path) -> Vec<String> {
        let out = Command::new("samtools").args(["view", "-F", "0x900"]).arg(bam).output().unwrap();
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    }

    fn bam_records(bam: &Path) -> BTreeSet<String> {
        let out = Command::new("samtools").args(["view"]).arg(bam).output().unwrap();
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .map(|l| l.split('\t').take(6).collect::<Vec<_>>().join("\t"))
            .collect()
    }

    #[test]
    fn end_to_end_filter_align_compact_resume() {
        if !tools_available() {
            return;
        }
        let mut fx = Fixture::new();
        let mut expect: BTreeMap<&str, Vec<String>> = BTreeMap::new();
        expect.entry("barcode01").or_default().extend(fx.write_chunk("fastq_pass/barcode01/a_0.fastq.gz", 60, true));
        expect.entry("barcode01").or_default().extend(fx.write_chunk("fastq_pass/barcode01/a_1.fastq", 40, false));
        expect.entry("barcode01").or_default().extend(fx.write_chunk("fastq_fail/barcode01/f_0.fastq.gz", 20, true));
        expect.entry("barcode02").or_default().extend(fx.write_chunk("fastq_pass/barcode02/b_0.fastq.gz", 50, true));

        let cfg = fx.config();
        let mut engine = Engine::new(cfg.clone()).unwrap();
        engine.run(&AtomicBool::new(false)).unwrap();

        for (barcode, names) in &expect {
            let b = &engine.state.barcodes[*barcode];
            assert_eq!(b.reads_passed as usize, names.len(), "{barcode} passed");
            assert_eq!(b.primary as usize, names.len(), "{barcode} primary records");
            assert_eq!(b.segments_pending, 0);
            let bam = PathBuf::from(b.compacted_bam.clone().expect("compacted bam"));
            assert!(bam.with_extension("bam.bai").exists(), "index next to {}", bam.display());
            let got: BTreeSet<_> = bam_names(&bam).into_iter().collect();
            let want: BTreeSet<_> = names.iter().cloned().collect();
            assert_eq!(got, want, "{barcode} read set");
        }
        assert_eq!(engine.state.files_done, 4);
        assert_eq!(engine.state.status, "finished");
        assert!(cfg.out_dir.join("provenance.json").exists());

        // Restart: nothing is redone, the counters carry on.
        let before = engine.state.barcodes["barcode01"].reads_passed;
        drop(engine);
        let mut engine = Engine::new(cfg.clone()).unwrap();
        let r = engine.step(true).unwrap();
        assert_eq!((r.files, r.aligned), (0, false));
        assert_eq!(engine.state.barcodes["barcode01"].reads_passed, before);

        // A new chunk is the only thing processed, and is merged into the existing BAM.
        let more = fx.write_chunk("fastq_pass/barcode01/a_2.fastq.gz", 30, true);
        let r = engine.step(true).unwrap();
        assert_eq!(r.files, 1);
        engine.compact_all().unwrap();
        let b = &engine.state.barcodes["barcode01"];
        assert_eq!(b.reads_passed, before + more.len() as u64);
        let bam = PathBuf::from(b.compacted_bam.clone().unwrap());
        let mut want: BTreeSet<String> = expect["barcode01"].iter().cloned().collect();
        want.extend(more);
        assert_eq!(bam_names(&bam).into_iter().collect::<BTreeSet<_>>(), want);
        assert_eq!(b.compacted_generation, 2);
    }

    #[test]
    fn alignments_equal_one_shot_minimap2_on_the_fasta() {
        if !tools_available() {
            return;
        }
        let mut fx = Fixture::new();
        let mut names = fx.write_chunk("fastq_pass/barcode01/a_0.fastq", 80, false);
        names.extend(fx.write_chunk("fastq_pass/barcode01/a_1.fastq", 80, false));

        let mut engine = Engine::new(fx.config()).unwrap();
        engine.run(&AtomicBool::new(false)).unwrap();
        let engine_bam = PathBuf::from(engine.state.barcodes["barcode01"].compacted_bam.clone().unwrap());

        // The reference: nagini's way -- filter the whole FASTQ with the same code,
        // then one `minimap2 -a -x splice -y -Y --MD` against the FASTA.
        let all = fx.root.join("all.fastq");
        let mut cat = Vec::new();
        for f in ["a_0.fastq", "a_1.fastq"] {
            cat.extend(std::fs::read(fx.root.join("run/fastq_pass/barcode01").join(f)).unwrap());
        }
        std::fs::write(&all, cat).unwrap();
        let filtered = fx.root.join("filtered.fastq");
        nanofilter_core::fastq::filter_fastq(&all, &filtered, 20.0, 500, 3000, false).unwrap();
        let sam = Command::new("minimap2")
            .args(["-a", "-x", "splice", "-t", "2", "-y", "-Y", "--MD"])
            .arg(&fx.reference)
            .arg(&filtered)
            .output()
            .unwrap();
        assert!(sam.status.success());
        let oneshot_sam = fx.root.join("oneshot.sam");
        std::fs::write(&oneshot_sam, &sam.stdout).unwrap();

        let oneshot = bam_records(&oneshot_sam);
        let ours = bam_records(&engine_bam);
        assert!(!oneshot.is_empty());
        assert_eq!(ours, oneshot, "name/flag/contig/pos/mapq/cigar of every record");
        assert_eq!(bam_names(&engine_bam).len(), names.len());
    }

    #[test]
    fn a_file_still_being_written_is_retried_not_lost() {
        if !tools_available() {
            return;
        }
        let mut fx = Fixture::new();
        let names = fx.write_chunk("fastq_pass/barcode01/a_0.fastq.gz", 60, true);
        let path = fx.root.join("run/fastq_pass/barcode01/a_0.fastq.gz");
        let full = std::fs::read(&path).unwrap();
        std::fs::write(&path, &full[..full.len() / 2]).unwrap(); // writer is mid-file

        let mut cfg = fx.config();
        cfg.once = false;
        let mut engine = Engine::new(cfg).unwrap();
        let r = engine.step(true).unwrap();
        assert_eq!((r.files, r.waiting), (0, 1));
        assert_eq!(engine.state.files_done, 0);

        std::fs::write(&path, &full).unwrap(); // writer finishes
        engine.retries.clear(); // skip the backoff delay
        // The size changed since the last scan, so it is held back one more round.
        let r = engine.step(true).unwrap();
        assert_eq!((r.files, r.waiting), (0, 1));
        let r = engine.step(true).unwrap();
        assert_eq!(r.files, 1);
        engine.compact_all().unwrap();
        let bam = PathBuf::from(engine.state.barcodes["barcode01"].compacted_bam.clone().unwrap());
        assert_eq!(bam_names(&bam).len(), names.len());
    }

    #[test]
    fn a_crash_between_manifest_and_cleanup_does_not_double_count() {
        if !tools_available() {
            return;
        }
        let mut fx = Fixture::new();
        let names = fx.write_chunk("fastq_pass/barcode01/a_0.fastq.gz", 50, true);
        let mut engine = Engine::new(fx.config()).unwrap();
        engine.run(&AtomicBool::new(false)).unwrap();
        let inbox = fx.root.join("out/barcode01/inbox");

        // Simulate the crash: the segment the manifest already consumed is still there.
        let consumed = Store::new(&fx.root.join("out/barcode01")).manifest().consumed;
        assert_eq!(consumed.len(), 1);
        std::fs::copy(
            fx.root.join("out/barcode01").join(Store::new(&fx.root.join("out/barcode01")).manifest().bam.unwrap()),
            inbox.join(&consumed[0]),
        )
        .unwrap();

        let mut cfg = fx.config();
        cfg.once = true;
        let mut engine = Engine::new(cfg).unwrap();
        engine.compact_all().unwrap();
        let bam = PathBuf::from(engine.state.barcodes["barcode01"].compacted_bam.clone().unwrap());
        assert_eq!(bam_names(&bam).len(), names.len());
    }

    #[test]
    fn each_snapshot_gets_a_preliminary_analysis_and_the_state_says_so() {
        if !tools_available() {
            return;
        }
        let mut fx = Fixture::new();
        fx.write_chunk("fastq_pass/barcode01/a_0.fastq.gz", 60, true);
        fx.write_chunk("fastq_pass/barcode02/b_0.fastq.gz", 60, true);
        let mut cfg = fx.config();
        cfg.prelim = Some(crate::prelim::stub::nagini(&fx.root.join("nagini"), false));

        let mut engine = Engine::new(cfg.clone()).unwrap();
        engine.run(&AtomicBool::new(false)).unwrap();

        for barcode in ["barcode01", "barcode02"] {
            let b = &engine.state.barcodes[barcode];
            assert_eq!(b.prelim_status.as_deref(), Some("ok"), "{barcode}");
            assert_eq!(b.prelim_generation, Some(1));
            let summary: serde_json::Value =
                serde_json::from_slice(&std::fs::read(b.prelim_summary.clone().unwrap()).unwrap()).unwrap();
            assert_eq!(summary["barcode"], barcode);
            assert_eq!(summary["reads_passed"].as_u64().unwrap(), b.reads_passed);
            assert_eq!(summary["fusion"]["primary_isoform"], "e14a2");
        }
        // A second generation appends to the history instead of replacing it.
        drop(engine);
        fx.write_chunk("fastq_pass/barcode01/a_1.fastq.gz", 30, true);
        let mut engine = Engine::new(cfg).unwrap();
        engine.run(&AtomicBool::new(false)).unwrap();
        let history = std::fs::read_to_string(fx.root.join("out/barcode01/prelim/history.jsonl")).unwrap();
        assert_eq!(history.lines().count(), 2);
        assert_eq!(engine.state.barcodes["barcode01"].prelim_generation, Some(2));
    }
}
