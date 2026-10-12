//! The "Final analysis" tab: hand the alignments to nagini, follow its run, and
//! show its results next to the preliminary ones.

use super::view::{fusion_section, mutation_section, phasing_section, section_title, AMBER, GREEN, GREY, RED};
use super::data::{format_age, format_count};
use eframe::egui::{self, RichText};
use nanobcr_core::handoff::{self, Handoff, MapRow};
use nanobcr_core::pipeline::{self, PipelineConfig, PipelineState, PipelineStatus};
use nanobcr_core::prelim::read_bcr_abl_section;
use nanobcr_core::results::{self, Comparison, Difference, FinalSample};
use nanobcr_core::state::now_unix;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver};
use std::time::{Duration, Instant};

const POLL_EVERY: Duration = Duration::from_millis(1500);

pub struct FinalEnv<'a> {
    pub out_dir: Option<&'a Path>,
    pub nagini_dir: Option<&'a Path>,
    pub engine_alive: bool,
    pub barcodes: &'a [String],
    pub threads: usize,
    pub build: &'a str,
}

enum TaskOutcome {
    Handoff(Result<Handoff, String>),
    DryRun(Result<pipeline::DryRun, String>),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum TaskKind {
    Handoff,
    DryRun,
}

pub struct FinalRun {
    map_path: String,
    handoff_override: String,
    analysis_override: String,
    site: String,
    cores: usize,
    rindels: bool,
    map: Option<(String, Result<Vec<MapRow>, String>)>,
    handoff: Option<Handoff>,
    pipeline: Option<PipelineStatus>,
    dry: Option<(bool, String)>,
    task: Option<(TaskKind, Receiver<TaskOutcome>)>,
    message: Option<(bool, String)>,
    results: Vec<FinalSample>,
    comparisons: BTreeMap<String, Comparison>,
    results_for_run: Option<u64>,
    selected: Option<String>,
    textures: HashMap<PathBuf, egui::TextureHandle>,
    last_poll: Option<Instant>,
}

impl Default for FinalRun {
    fn default() -> Self {
        Self {
            map_path: String::new(),
            handoff_override: String::new(),
            analysis_override: String::new(),
            site: "workstation".into(),
            cores: 0,
            rindels: false,
            map: None,
            handoff: None,
            pipeline: None,
            dry: None,
            task: None,
            message: None,
            results: Vec::new(),
            comparisons: BTreeMap::new(),
            results_for_run: None,
            selected: None,
            textures: HashMap::new(),
            last_poll: None,
        }
    }
}

/// Open a file or folder with the system's default application.
pub fn open_path(p: &Path) {
    let opener = if cfg!(target_os = "macos") { "open" } else { "xdg-open" };
    let _ = std::process::Command::new(opener).arg(p).spawn();
}

impl FinalRun {
    #[cfg(test)]
    pub fn results_len(&self) -> usize {
        self.results.len()
    }

    #[cfg(test)]
    pub fn comparison_identical(&self, barcode: &str) -> Option<bool> {
        self.comparisons.get(barcode).map(|c| c.identical())
    }

    fn handoff_dir(&self, env: &FinalEnv) -> Option<PathBuf> {
        let t = self.handoff_override.trim();
        if !t.is_empty() {
            return Some(PathBuf::from(t));
        }
        env.out_dir.map(|o| o.join("handoff"))
    }

    fn analysis_dir(&self, env: &FinalEnv) -> Option<PathBuf> {
        let t = self.analysis_override.trim();
        if !t.is_empty() {
            return Some(PathBuf::from(t));
        }
        env.out_dir.map(|o| o.join("nagini_run"))
    }

    fn pipeline_config(&self, env: &FinalEnv) -> Option<PipelineConfig> {
        let (nagini, handoff, analysis) = (env.nagini_dir?, self.handoff_dir(env)?, self.analysis_dir(env)?);
        let map = {
            // nagini reads the map the hand-off wrote (only the barcodes handed over).
            let used = handoff.join("barcode_map.csv");
            if used.is_file() { used } else { PathBuf::from(self.map_path.trim()) }
        };
        let mut c = PipelineConfig::new(nagini.to_path_buf(), analysis.clone(), handoff, map);
        c.site = self.site.clone();
        c.build = env.build.to_string();
        c.cores = (self.cores > 0).then_some(self.cores);
        c.enable_rindels = self.rindels;
        c.sequence_id = analysis.file_name().map(|n| n.to_string_lossy().to_string());
        Some(c)
    }

    fn min_depth(&self, env: &FinalEnv) -> u64 {
        env.nagini_dir
            .map(|n| read_bcr_abl_section(&n.join("config/config.yaml")))
            .and_then(|m| m.get("min_depth").and_then(|v| v.parse().ok()))
            .unwrap_or(50)
    }

    /// Refresh what lives on disk; cheap, so it runs every 1.5 s.
    pub fn poll(&mut self, env: &FinalEnv) {
        if let Some((kind, rx)) = &self.task {
            if let Ok(outcome) = rx.try_recv() {
                let kind = *kind;
                self.task = None;
                match (kind, outcome) {
                    (_, TaskOutcome::Handoff(Ok(h))) => {
                        self.message = Some((true, format!("Hand-off ready: {} sample(s)", h.barcodes.len())));
                        self.handoff = Some(h);
                    }
                    (_, TaskOutcome::Handoff(Err(e))) => self.message = Some((false, e)),
                    (_, TaskOutcome::DryRun(Ok(d))) => {
                        self.message = Some((d.ok, if d.ok { "Dry run succeeded".into() } else { "Dry run failed: see the output below".into() }));
                        self.dry = Some((d.ok, d.tail));
                    }
                    (_, TaskOutcome::DryRun(Err(e))) => self.message = Some((false, e)),
                }
            }
        }
        if self.last_poll.map(|t| t.elapsed() < POLL_EVERY).unwrap_or(false) {
            return;
        }
        self.last_poll = Some(Instant::now());

        if !self.map_path.trim().is_empty() && self.map.as_ref().map(|(p, _)| p != &self.map_path).unwrap_or(true) {
            let r = handoff::read_map(Path::new(self.map_path.trim())).map_err(|e| format!("{e:#}"));
            self.map = Some((self.map_path.clone(), r));
        }
        if let Some(dir) = self.handoff_dir(env) {
            if let Some(h) = handoff::load(&dir) {
                self.handoff = Some(h);
            }
        }
        let Some(analysis) = self.analysis_dir(env) else { return };
        let st = pipeline::status(&analysis);
        if st.state == PipelineState::Succeeded && st.started_unix != self.results_for_run {
            if let (Some(h), Some(out)) = (self.handoff.clone(), env.out_dir) {
                self.results = results::load_all(&analysis, &h, self.min_depth(env));
                self.comparisons.clear();
                for r in &self.results {
                    let snapshot_gen = h.barcodes.iter().find(|b| b.barcode == r.barcode).map(|b| b.snapshot_generation).unwrap_or(0);
                    if let Ok(Some(p)) = results::load_prelim(out, &r.barcode) {
                        self.comparisons.insert(r.barcode.clone(), results::compare(&p, r, snapshot_gen));
                    }
                }
                self.textures.clear();
                self.results_for_run = st.started_unix;
                if self.selected.is_none() {
                    self.selected = self.results.first().map(|r| r.barcode.clone());
                }
            }
        }
        self.pipeline = Some(st);
    }

    pub fn show(&mut self, ui: &mut egui::Ui, env: &FinalEnv) {
        self.poll(env);
        egui::ScrollArea::vertical().show(ui, |ui| {
            if let Some((ok, msg)) = &self.message {
                ui.label(RichText::new(msg).color(if *ok { GREEN } else { RED }));
            }
            if self.task.is_some() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Working…");
                });
            }
            self.step_map(ui, env);
            ui.add_space(10.0);
            self.step_handoff(ui, env);
            ui.add_space(10.0);
            self.step_run(ui, env);
            ui.add_space(10.0);
            self.step_results(ui, env);
        });
    }

    // ---------------------------------------------------------------- 1

    fn step_map(&mut self, ui: &mut egui::Ui, env: &FinalEnv) {
        section_title(ui, "1 · Barcode map");
        ui.label("Which sample each barcode is. nagini needs this to name the results (CSV: barcode, sample, type T/N/R, assay, library_type).");
        ui.horizontal(|ui| {
            ui.add(egui::TextEdit::singleline(&mut self.map_path).desired_width(420.0).hint_text("barcode_map.csv"));
            if ui.button("…").clicked() {
                if let Some(p) = rfd::FileDialog::new().pick_file() {
                    self.map_path = p.to_string_lossy().into();
                }
            }
            let can = env.out_dir.is_some() && !env.barcodes.is_empty();
            if ui.add_enabled(can, egui::Button::new("Create template from this run")).clicked() {
                let path = env.out_dir.unwrap().join("barcode_map.csv");
                if path.exists() {
                    self.message = Some((false, format!("{} already exists: not overwritten", path.display())));
                } else {
                    match handoff::write_map(&path, &handoff::template_for(env.barcodes)) {
                        Ok(()) => {
                            self.map_path = path.to_string_lossy().into();
                            self.message = Some((true, "Template written: put the real sample names in the `sample` column, then save".into()));
                            open_path(&path);
                        }
                        Err(e) => self.message = Some((false, format!("{e:#}"))),
                    }
                }
            }
        });
        match self.map.as_ref().map(|(_, r)| r) {
            Some(Ok(rows)) => {
                ui.label(RichText::new(format!("{} sample(s): {}", rows.len(), rows.iter().map(|r| format!("{}={}", r.barcode, r.sample)).collect::<Vec<_>>().join(", "))).color(GREEN));
                let missing: Vec<_> = env.barcodes.iter().filter(|b| *b != "unclassified" && !rows.iter().any(|r| &r.barcode == *b)).cloned().collect();
                if !missing.is_empty() {
                    ui.label(RichText::new(format!("Not in the map, so not handed to nagini: {}", missing.join(", "))).color(AMBER));
                }
            }
            Some(Err(e)) => {
                ui.label(RichText::new(e).color(RED));
            }
            None => {}
        }
    }

    // ---------------------------------------------------------------- 2

    fn step_handoff(&mut self, ui: &mut egui::Ui, env: &FinalEnv) {
        section_title(ui, "2 · Hand over the alignments");
        ui.label("Folds the last reads into each barcode's sorted BAM and gathers them for nagini. The engine must be stopped first (a clean stop does the final compaction).");
        let map_ok = matches!(self.map.as_ref().map(|(_, r)| r), Some(Ok(_)));
        let can = map_ok && !env.engine_alive && env.out_dir.is_some() && self.task.is_none();
        ui.horizontal(|ui| {
            let b = ui.add_enabled(can, egui::Button::new("Prepare hand-off"));
            if b.clicked() {
                let (out, dir, map, threads) = (env.out_dir.unwrap().to_path_buf(), self.handoff_dir(env).unwrap(), PathBuf::from(self.map_path.trim()), env.threads);
                let (tx, rx) = mpsc::channel();
                std::thread::spawn(move || {
                    let _ = tx.send(TaskOutcome::Handoff(handoff::prepare(&out, &dir, &map, "samtools", threads).map_err(|e| format!("{e:#}"))));
                });
                self.task = Some((TaskKind::Handoff, rx));
                self.message = None;
            }
            if env.engine_alive {
                ui.label(RichText::new("The engine is still running: stop it first.").color(AMBER));
            } else if !map_ok {
                ui.label(RichText::new("Needs a valid barcode map.").color(GREY));
            }
        });
        let Some(h) = &self.handoff else { return };
        ui.label(RichText::new(format!("Hand-off folder: {}", h.dir.display())).small().color(GREY));
        egui::Grid::new("bcr_handoff").striped(true).spacing([14.0, 4.0]).show(ui, |ui| {
            for c in ["Barcode", "Sample", "Type", "Snapshot", "Primary alignments", "Reads passed (engine)"] {
                ui.label(RichText::new(c).strong());
            }
            ui.end_row();
            for b in &h.barcodes {
                ui.label(&b.barcode);
                ui.label(&b.sample);
                ui.label(&b.sample_type);
                ui.label(format!("gen {}", b.snapshot_generation));
                ui.label(format_count(b.primary_records));
                ui.label(format_count(b.reads_passed));
                ui.end_row();
            }
        });
        for (b, why) in &h.skipped {
            ui.label(RichText::new(format!("Skipped {b}: {why}")).color(AMBER));
        }
        for w in &h.warnings {
            ui.label(RichText::new(w).color(AMBER));
        }
        if let Ok(text) = std::fs::read_to_string(h.dir.join("provenance.json")) {
            if let Ok(v) = serde_json::from_str::<serde_json::Value>(&text) {
                ui.label(
                    RichText::new(format!(
                        "Filtered qs ≥ {} · {}–{} bp · {} · -x {} {}  (nagini does not filter or realign these BAMs; its own read filter does not apply)",
                        v["read_filter"]["min_qv"], v["read_filter"]["min_len"], v["read_filter"]["max_len"], v["minimap2"].as_str().unwrap_or("minimap2"),
                        v["alignment"]["preset"].as_str().unwrap_or("splice"),
                        v["alignment"]["extra"].as_array().map(|a| a.iter().filter_map(|x| x.as_str()).collect::<Vec<_>>().join(" ")).unwrap_or_default()
                    ))
                    .small()
                    .color(GREY),
                );
            }
        }
    }

    // ---------------------------------------------------------------- 3

    fn step_run(&mut self, ui: &mut egui::Ui, env: &FinalEnv) {
        section_title(ui, "3 · Run nagini");
        let Some(cfg) = self.pipeline_config(env) else {
            ui.label(RichText::new("Set the output folder and the nagini checkout in the setup panel.").color(GREY));
            return;
        };
        ui.horizontal(|ui| {
            ui.label("Site");
            egui::ComboBox::from_id_salt("bcr_site").selected_text(self.site.clone()).show_ui(ui, |ui| {
                for s in ["workstation", "local", "marvin"] {
                    ui.selectable_value(&mut self.site, s.to_string(), s);
                }
            });
            ui.label("Cores (0 = profile default)");
            ui.add(egui::DragValue::new(&mut self.cores).range(0..=256));
            ui.checkbox(&mut self.rindels, "Rindels calls (--enable-rindels)");
        });
        ui.collapsing("Command", |ui| {
            ui.label(RichText::new(cfg.command(false).join(" ")).monospace().small());
            ui.label(RichText::new(format!("Analysis folder: {}", cfg.analysis_dir.display())).small().color(GREY));
        });

        let status = self.pipeline.clone();
        let running = matches!(status.as_ref().map(|s| &s.state), Some(PipelineState::Running { .. }));
        let have_handoff = self.handoff.as_ref().map(|h| !h.barcodes.is_empty()).unwrap_or(false);
        ui.horizontal(|ui| {
            if ui.add_enabled(have_handoff && !running && self.task.is_none(), egui::Button::new("Check (dry run)")).on_hover_text("Asks snakemake whether the workflow can be built here, without running it").clicked() {
                let c = cfg.clone();
                let (tx, rx) = mpsc::channel();
                std::thread::spawn(move || {
                    let _ = tx.send(TaskOutcome::DryRun(pipeline::dry_run(&c).map_err(|e| format!("{e:#}"))));
                });
                self.task = Some((TaskKind::DryRun, rx));
            }
            if ui.add_enabled(have_handoff && !running && self.task.is_none(), egui::Button::new("Start nagini")).clicked() {
                match pipeline::start(&cfg) {
                    Ok(rec) => self.message = Some((true, format!("nagini started (pid {})", rec.pid))),
                    Err(e) => self.message = Some((false, format!("{e:#}"))),
                }
                self.last_poll = None;
            }
            if ui.add_enabled(running, egui::Button::new("Stop run")).clicked() {
                match pipeline::stop(&cfg.analysis_dir) {
                    Ok(()) => self.message = Some((true, "Interrupt sent; snakemake finishes the running jobs".into())),
                    Err(e) => self.message = Some((false, format!("{e:#}"))),
                }
            }
            if !have_handoff {
                ui.label(RichText::new("Prepare the hand-off first.").color(GREY));
            }
        });
        if let Some((ok, text)) = &self.dry {
            ui.collapsing(RichText::new(if *ok { "Dry run output (ok)" } else { "Dry run output (FAILED)" }).color(if *ok { GREEN } else { RED }), |ui| {
                ui.label(RichText::new(text).monospace().small());
            });
        }

        let Some(st) = status else { return };
        let (text, color) = match &st.state {
            PipelineState::NotStarted => ("Not started".to_string(), GREY),
            PipelineState::Running { pid } => (format!("Running (pid {pid}) · started {} ago", st.started_unix.map(|s| format_age(now_unix().saturating_sub(s))).unwrap_or_default()), GREEN),
            PipelineState::Succeeded => ("Finished: the summary workbook was produced".to_string(), GREEN),
            PipelineState::Failed { reason, .. } => (format!("FAILED: {reason}"), RED),
            PipelineState::Lost => ("The nagini process is gone without recording an exit code (killed, or the machine restarted)".to_string(), RED),
        };
        ui.label(RichText::new(text).color(color).strong());
        if st.progress.total > 0 {
            ui.add(egui::ProgressBar::new(st.progress.done as f32 / st.progress.total as f32).text(format!("{} of {} steps ({}%)", st.progress.done, st.progress.total, st.progress.percent)));
            if !st.progress.recent_rules.is_empty() {
                ui.label(RichText::new(format!("Recent: {}", st.progress.recent_rules.join(" → "))).small().color(GREY));
            }
        }
        if !st.progress.failed_rules.is_empty() {
            ui.label(RichText::new(format!("Failed rules: {}", st.progress.failed_rules.join(", "))).color(RED));
        }
        if !st.log_tail.is_empty() {
            ui.collapsing("nagini log", |ui| {
                egui::ScrollArea::vertical().max_height(220.0).stick_to_bottom(true).show(ui, |ui| {
                    for l in st.log_tail.iter().rev().take(80).rev() {
                        ui.label(RichText::new(l).monospace().small());
                    }
                });
            });
        }
    }

    // ---------------------------------------------------------------- 4

    fn texture(&mut self, ctx: &egui::Context, path: &Path) -> Option<egui::TextureHandle> {
        if let Some(t) = self.textures.get(path) {
            return Some(t.clone());
        }
        let img = image::open(path).ok()?.to_rgba8();
        let size = [img.width() as usize, img.height() as usize];
        let tex = ctx.load_texture(path.to_string_lossy(), egui::ColorImage::from_rgba_unmultiplied(size, img.as_raw()), Default::default());
        self.textures.insert(path.to_path_buf(), tex.clone());
        Some(tex)
    }

    fn step_results(&mut self, ui: &mut egui::Ui, env: &FinalEnv) {
        section_title(ui, "4 · nagini's results");
        if self.results.is_empty() {
            ui.label(RichText::new("Shown here when the run has finished.").color(GREY));
            return;
        }
        let analysis = self.analysis_dir(env);
        ui.horizontal(|ui| {
            if let Some(w) = analysis.as_deref().and_then(results::summary_workbook) {
                if ui.button("Open summary workbook").clicked() {
                    open_path(&w);
                }
            }
            if let Some(m) = analysis.as_deref().and_then(results::multiqc_report) {
                if ui.button("Open MultiQC").clicked() {
                    open_path(&m);
                }
            }
            if let Some(d) = self.results.first().and_then(|r| r.files.reports_dir.clone()) {
                if ui.button("Open reports folder").clicked() {
                    open_path(&d);
                }
            }
        });
        ui.horizontal(|ui| {
            ui.label("Sample");
            for r in &self.results {
                ui.selectable_value(&mut self.selected, Some(r.barcode.clone()), format!("{} ({})", r.sample, r.barcode));
            }
        });
        let Some(r) = self.results.iter().find(|r| Some(&r.barcode) == self.selected.as_ref()).cloned() else { return };
        ui.add_space(6.0);
        ui.heading(format!("{} · {}", r.sample, r.sample_type));
        for p in &r.problems {
            ui.label(RichText::new(p).color(AMBER));
        }
        if let Some(w) = &r.files.sample_workbook {
            if ui.button("Open this sample's workbook").clicked() {
                open_path(w);
            }
        }
        self.comparison_section(ui, &r.barcode);
        ui.add_space(8.0);
        fusion_section(ui, r.fusion.as_ref(), false);
        if let Some(png) = &r.files.fusion_plot_png {
            if let Some(tex) = self.texture(ui.ctx(), png) {
                let max_w = ui.available_width().min(900.0);
                let scale = (max_w / tex.size_vec2().x).min(1.0);
                ui.add(egui::Image::new((tex.id(), tex.size_vec2() * scale)));
            }
        }
        ui.add_space(8.0);
        mutation_section(ui, r.branford.as_ref());
        ui.add_space(8.0);
        phasing_section(ui, r.phasing.as_ref());
        ui.add_space(8.0);
        section_title(ui, "Exon coverage");
        if r.exons.is_empty() {
            ui.label("Not available.");
        } else {
            egui::Grid::new("bcr_exons").striped(true).spacing([14.0, 4.0]).show(ui, |ui| {
                for c in ["Exon", "Mean depth", "≥100×", "Spanning reads"] {
                    ui.label(RichText::new(c).strong());
                }
                ui.end_row();
                for e in r.exons.iter().filter(|e| e.mean_depth > 0.0) {
                    ui.label(&e.exon);
                    ui.label(format!("{:.0}×", e.mean_depth));
                    ui.label(format!("{:.0}%", e.pct_gte_100x));
                    ui.label(format_count(e.spanning_reads));
                    ui.end_row();
                }
            });
        }
    }

    fn comparison_section(&self, ui: &mut egui::Ui, barcode: &str) {
        section_title(ui, "Preliminary vs final");
        let Some(c) = self.comparisons.get(barcode) else {
            ui.label(RichText::new("No preliminary result was recorded for this sample.").color(GREY));
            return;
        };
        if c.identical() {
            ui.label(RichText::new(format!("Identical: the final result matches the preliminary one (snapshot gen {} vs hand-off gen {}).", c.preliminary_generation, c.handoff_generation)).color(GREEN));
            return;
        }
        if c.preliminary_generation < c.handoff_generation {
            ui.label(RichText::new(format!("The preliminary result is from an older snapshot (gen {} of {}), so differences are expected.", c.preliminary_generation, c.handoff_generation)).color(AMBER));
        } else {
            ui.label(RichText::new("The preliminary result was computed on the same snapshot: any difference means nagini used different data or settings.").color(RED));
        }
        for d in &c.differences {
            let text = match d {
                Difference::Isoform { preliminary, final_call } => format!("Isoform: {} → {}", preliminary.as_deref().unwrap_or("none"), final_call.as_deref().unwrap_or("none")),
                Difference::MutationOnlyInPreliminary { mutation, vaf_percent } => format!("{mutation} ({vaf_percent:.2}%) was called preliminarily but not in the final result"),
                Difference::MutationOnlyInFinal { mutation, vaf_percent } => format!("{mutation} ({vaf_percent:.2}%) appears only in the final result"),
                Difference::VafShift { mutation, preliminary, final_vaf } => format!("{mutation}: VAF {preliminary:.2}% → {final_vaf:.2}%"),
                Difference::PhasingChanged { pair, preliminary, final_state } => format!("Phasing {pair}: {} → {}", preliminary.as_deref().unwrap_or("none"), final_state.as_deref().unwrap_or("none")),
            };
            ui.label(RichText::new(format!("• {text}")).color(AMBER));
        }
    }
}
