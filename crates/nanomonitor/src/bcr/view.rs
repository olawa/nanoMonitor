use super::data::{self, BcrData, Liveness};
use super::final_run::{FinalEnv, FinalRun};
use super::proc::{self, BcrSetup, EngineProc};
use eframe::egui::{self, Color32, RichText};
use egui_plot::{Legend, Line, Plot, PlotPoints};
use nanobcr_core::preflight::{self, Check, CheckStatus};
use nanobcr_core::prelim::{BranfordSummary, FusionSummary, PhasePair, PrelimSummary};
use nanobcr_core::state::{now_unix, BarcodeState};
use std::path::PathBuf;
use std::time::{Duration, Instant};

const RELOAD_EVERY: Duration = Duration::from_millis(1500);

pub(super) const GREEN: Color32 = Color32::from_rgb(80, 170, 100);
pub(super) const AMBER: Color32 = Color32::from_rgb(230, 170, 40);
pub(super) const RED: Color32 = Color32::from_rgb(220, 70, 70);
pub(super) const GREY: Color32 = Color32::from_rgb(140, 140, 150);

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tab {
    Samples,
    Detail,
    Final,
    Log,
}

pub struct BcrApp {
    setup: BcrSetup,
    data: BcrData,
    last_load: Option<Instant>,
    engine: EngineProc,
    checks: Vec<Check>,
    selected: Option<String>,
    tab: Tab,
    message: Option<(bool, String)>,
    log_tail: Vec<String>,
    final_run: FinalRun,
}

pub struct BcrStartup {
    pub setup: BcrSetup,
}

impl BcrApp {
    pub fn new(startup: BcrStartup) -> Self {
        let mut app = Self {
            setup: startup.setup,
            data: BcrData::default(),
            last_load: None,
            engine: EngineProc::default(),
            checks: Vec::new(),
            selected: None,
            tab: Tab::Samples,
            message: None,
            log_tail: Vec::new(),
            final_run: FinalRun::default(),
        };
        app.reload();
        app
    }

    fn out_dir(&self) -> Option<PathBuf> {
        let t = self.setup.out_dir.trim();
        (!t.is_empty()).then(|| PathBuf::from(t))
    }

    fn reload(&mut self) {
        self.last_load = Some(Instant::now());
        match self.out_dir() {
            Some(dir) => {
                self.data = data::load(&dir);
                self.log_tail = proc::tail_log(&dir, 400);
            }
            None => {
                self.data = BcrData::default();
                self.log_tail.clear();
            }
        }
        if self.selected.is_none() || !self.known_barcodes().contains(self.selected.as_ref().unwrap()) {
            self.selected = self.known_barcodes().into_iter().next();
        }
    }

    fn known_barcodes(&self) -> Vec<String> {
        let mut v: Vec<String> = self
            .data
            .state
            .as_ref()
            .map(|s| s.barcodes.keys().cloned().collect())
            .unwrap_or_default();
        for k in self.data.prelim.keys() {
            if !v.contains(k) {
                v.push(k.clone());
            }
        }
        v.sort();
        v
    }

    fn liveness(&self) -> Liveness {
        data::liveness(self.data.state.as_ref(), now_unix())
    }

    fn engine_alive(&mut self) -> bool {
        if self.engine.running() {
            return true;
        }
        matches!(self.liveness(), Liveness::Running { .. })
            && self.data.state.as_ref().map(|s| proc::pid_alive(s.pid)).unwrap_or(false)
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        if self.last_load.map(|t| t.elapsed() >= RELOAD_EVERY).unwrap_or(true) {
            self.reload();
        }
        if let Some(msg) = self.engine.exit_message() {
            self.message = Some((false, msg));
        }
        ctx.request_repaint_after(Duration::from_secs(1));

        egui::TopBottomPanel::top("bcr_status").show(ctx, |ui| self.status_bar(ui));
        egui::SidePanel::left("bcr_setup")
            .default_width(360.0)
            .resizable(true)
            .show(ctx, |ui| egui::ScrollArea::vertical().show(ui, |ui| self.setup_panel(ui)));
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.tab, Tab::Samples, "Samples");
                ui.selectable_value(&mut self.tab, Tab::Detail, "Sample detail");
                ui.selectable_value(&mut self.tab, Tab::Final, "Final analysis (nagini)");
                ui.selectable_value(&mut self.tab, Tab::Log, "Engine log");
            });
            ui.separator();
            match self.tab {
                Tab::Samples => self.samples_tab(ui),
                Tab::Detail => self.detail_tab(ui),
                Tab::Final => self.final_tab(ui),
                Tab::Log => self.log_tab(ui),
            }
        });
    }

    // ------------------------------------------------------------ status

    fn status_bar(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.heading("nanoMonitor · BCR-ABL1");
            ui.separator();
            let live = self.liveness();
            let (text, color) = match &live {
                Liveness::NotStarted => ("Engine not started".to_string(), GREY),
                Liveness::Running { heartbeat_age } => (format!("Engine running · updated {} ago", data::format_age(*heartbeat_age)), GREEN),
                Liveness::Stale { heartbeat_age } => (format!("Engine NOT RESPONDING · last update {} ago", data::format_age(*heartbeat_age)), RED),
                Liveness::Finished => ("Run finished".to_string(), GREEN),
                Liveness::Error(e) => (format!("Engine error: {e}"), RED),
            };
            ui.label(RichText::new(text).color(color).strong());
            if let Some(s) = &self.data.state {
                ui.separator();
                ui.label(format!(
                    "{} files · {} alignment cycles",
                    data::format_count(s.files_done),
                    data::format_count(s.cycles)
                ));
            }
            if let Some((ok, msg)) = &self.message {
                ui.separator();
                ui.label(RichText::new(msg).color(if *ok { GREEN } else { RED }));
            }
        });
        for p in &self.data.problems {
            ui.label(RichText::new(p).color(AMBER));
        }
        ui.add_space(4.0);
    }

    // ------------------------------------------------------------ setup

    fn path_row(ui: &mut egui::Ui, label: &str, value: &mut String, folder: bool) {
        ui.label(label);
        ui.horizontal(|ui| {
            ui.add(egui::TextEdit::singleline(value).desired_width(250.0));
            if ui.button("…").clicked() {
                let dialog = rfd::FileDialog::new();
                let picked = if folder { dialog.pick_folder() } else { dialog.pick_file() };
                if let Some(p) = picked {
                    *value = p.to_string_lossy().into();
                }
            }
        });
    }

    fn setup_panel(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.heading("Run setup");
        Self::path_row(ui, "MinKNOW run folder (or fastq_pass)", &mut self.setup.input_dir, true);
        Self::path_row(ui, "Output folder", &mut self.setup.out_dir, true);
        Self::path_row(ui, "Reference (hg19 FASTA, or a -x splice .mmi)", &mut self.setup.reference, false);
        Self::path_row(ui, "nagini checkout (for the preliminary analysis)", &mut self.setup.nagini_dir, true);
        egui::CollapsingHeader::new("Advanced").show(ui, |ui| {
            Self::path_row(ui, "Python with pysam (default: nagini's .venv)", &mut self.setup.python, false);
            Self::path_row(ui, "Reference FASTA (when the reference is a .mmi)", &mut self.setup.fasta, false);
            Self::path_row(ui, "nanostream binary (default: next to this app, then PATH)", &mut self.setup.nanostream_bin, false);
            ui.horizontal(|ui| {
                ui.label("Reference build");
                ui.text_edit_singleline(&mut self.setup.build);
            });
            ui.horizontal(|ui| {
                ui.label("Min qs");
                ui.add(egui::DragValue::new(&mut self.setup.qv).speed(0.1).range(0.0..=60.0));
                ui.label("Length");
                ui.add(egui::DragValue::new(&mut self.setup.min_len).range(0..=100_000));
                ui.label("–");
                ui.add(egui::DragValue::new(&mut self.setup.max_len).range(0..=1_000_000));
            });
            ui.horizontal(|ui| {
                ui.label("Threads");
                ui.add(egui::DragValue::new(&mut self.setup.threads).range(1..=128));
                ui.label("Snapshot every (s)");
                ui.add(egui::DragValue::new(&mut self.setup.snapshot_secs).range(30..=7200));
            });
            ui.checkbox(&mut self.setup.pass_only, "fastq_pass only (nagini includes fastq_fail)");
            ui.label(
                RichText::new("The final nagini run takes its filter from its own config: keep these equal to it.")
                    .small()
                    .color(GREY),
            );
        });

        ui.add_space(8.0);
        ui.horizontal(|ui| {
            if ui.button("Check setup").clicked() {
                self.checks = preflight::preflight(&self.setup.engine_config());
            }
            let blocked = preflight::worst(&self.checks) == CheckStatus::Fail && !self.checks.is_empty();
            let alive = self.engine_alive();
            if ui.add_enabled(!alive && !blocked, egui::Button::new("Start engine")).clicked() {
                self.start_engine();
            }
            if ui.add_enabled(alive, egui::Button::new("Stop (clean)")).clicked() {
                self.stop_engine();
            }
            if ui.add_enabled(self.engine.running(), egui::Button::new("Force stop")).on_hover_text("Kills the engine this window started. Use only if a clean stop does not finish.").clicked() {
                self.engine.kill();
                self.message = Some((false, "Engine killed; the next start resumes from the ledger".into()));
            }
        });
        for c in &self.checks {
            let (mark, color) = match c.status {
                CheckStatus::Pass => ("✔", GREEN),
                CheckStatus::Warn => ("▲", AMBER),
                CheckStatus::Fail => ("✖", RED),
            };
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new(mark).color(color));
                ui.label(RichText::new(&c.name).strong());
                ui.label(&c.detail);
            });
        }
        ui.add_space(6.0);
        ui.label(
            RichText::new("The engine keeps running if you close this window. Reopen it and point at the same output folder to watch it again.")
                .small()
                .color(GREY),
        );
    }

    fn start_engine(&mut self) {
        // Always re-check just before starting: the state of the machine may have changed.
        self.checks = preflight::preflight(&self.setup.engine_config());
        if preflight::worst(&self.checks) == CheckStatus::Fail {
            self.message = Some((false, "Not started: fix the failed checks first".into()));
            return;
        }
        if self.engine_alive() {
            self.message = Some((false, "An engine is already running on this output folder".into()));
            return;
        }
        match self.engine.start(&self.setup) {
            Ok(()) => self.message = Some((true, "Engine started".into())),
            Err(e) => self.message = Some((false, format!("{e:#}"))),
        }
    }

    fn stop_engine(&mut self) {
        let Some(dir) = self.out_dir() else { return };
        match EngineProc::request_stop(&dir) {
            Ok(()) => self.message = Some((true, "Stop requested: the engine finishes its cycle, compacts and exits".into())),
            Err(e) => self.message = Some((false, format!("{e:#}"))),
        }
    }

    // ------------------------------------------------------------ samples

    fn samples_tab(&mut self, ui: &mut egui::Ui) {
        let Some(state) = self.data.state.clone() else {
            ui.add_space(20.0);
            ui.label("No engine output yet. Fill in the run setup, check it and start the engine.");
            return;
        };
        if state.barcodes.is_empty() {
            ui.label("Waiting for the first FASTQ files…");
            return;
        }
        let now = now_unix();
        egui::ScrollArea::both().show(ui, |ui| {
            egui::Grid::new("bcr_samples").striped(true).num_columns(12).spacing([14.0, 6.0]).show(ui, |ui| {
                for h in ["Barcode", "Files", "Reads", "Passed filter", "Mapped", "Supplementary", "Last data", "Snapshot", "Analysis", "Fusion", "Kinase mutations", "Enough?"] {
                    ui.label(RichText::new(h).strong());
                }
                ui.end_row();
                for (barcode, b) in &state.barcodes {
                    let summary = self.data.prelim.get(barcode);
                    if ui.selectable_label(self.selected.as_deref() == Some(barcode), barcode).clicked() {
                        self.selected = Some(barcode.clone());
                        self.tab = Tab::Detail;
                    }
                    ui.label(data::format_count(b.files));
                    ui.label(data::format_count(b.reads_seen));
                    let pct = if b.reads_seen > 0 { 100.0 * b.reads_passed as f64 / b.reads_seen as f64 } else { 0.0 };
                    ui.label(format!("{} ({pct:.0}%)", data::format_count(b.reads_passed)));
                    ui.label(data::format_count(b.primary_mapped));
                    ui.label(data::format_count(b.supplementary));
                    let age = now.saturating_sub(b.last_update_unix);
                    ui.label(if b.last_update_unix == 0 { "-".into() } else { format!("{} ago", data::format_age(age)) });
                    ui.label(if b.compacted_generation == 0 { "-".into() } else { format!("gen {}", b.compacted_generation) });
                    prelim_status_label(ui, b, summary);
                    fusion_label(ui, summary);
                    mutation_label(ui, summary);
                    enough_label(ui, summary);
                    ui.end_row();
                }
            });
        });
    }

    // ------------------------------------------------------------ detail

    fn detail_tab(&mut self, ui: &mut egui::Ui) {
        let barcodes = self.known_barcodes();
        if barcodes.is_empty() {
            ui.label("No samples yet.");
            return;
        }
        ui.horizontal(|ui| {
            ui.label("Sample");
            egui::ComboBox::from_id_salt("bcr_sample_pick")
                .selected_text(self.selected.clone().unwrap_or_default())
                .show_ui(ui, |ui| {
                    for b in &barcodes {
                        ui.selectable_value(&mut self.selected, Some(b.clone()), b);
                    }
                });
        });
        let Some(barcode) = self.selected.clone() else { return };
        let state = self.data.state.as_ref().and_then(|s| s.barcodes.get(&barcode)).cloned();
        let Some(summary) = self.data.prelim.get(&barcode).cloned() else {
            ui.add_space(12.0);
            match state.as_ref().and_then(|b| b.prelim_status.clone()) {
                Some(s) => ui.label(format!("Preliminary analysis: {s}. The first one runs after the first snapshot.")),
                None => ui.label("No preliminary analysis yet. It needs the nagini checkout in the setup and runs after each snapshot."),
            };
            return;
        };
        let history = self.data.history.get(&barcode).cloned().unwrap_or_default();

        egui::ScrollArea::vertical().show(ui, |ui| {
            summary_header(ui, &summary, state.as_ref());
            ui.add_space(8.0);
            fusion_section(ui, summary.fusion.as_ref(), data::isoform_changed(&history));
            ui.add_space(8.0);
            mutation_section(ui, summary.branford.as_ref());
            ui.add_space(8.0);
            phasing_section(ui, summary.phasing.as_ref());
            ui.add_space(8.0);
            coverage_section(ui, &summary);
            ui.add_space(8.0);
            history_plots(ui, &history);
        });
    }

    fn final_tab(&mut self, ui: &mut egui::Ui) {
        let alive = self.engine_alive();
        let barcodes = self.known_barcodes();
        let out = self.out_dir();
        let nagini = {
            let t = self.setup.nagini_dir.trim();
            (!t.is_empty()).then(|| PathBuf::from(t))
        };
        let env = FinalEnv {
            out_dir: out.as_deref(),
            nagini_dir: nagini.as_deref(),
            engine_alive: alive,
            barcodes: &barcodes,
            threads: self.setup.threads,
            build: &self.setup.build,
        };
        self.final_run.show(ui, &env);
    }

    // ------------------------------------------------------------ log

    fn log_tab(&mut self, ui: &mut egui::Ui) {
        if self.log_tail.is_empty() {
            ui.label("The engine log (engine.log in the output folder) is empty.");
            return;
        }
        egui::ScrollArea::vertical().stick_to_bottom(true).show(ui, |ui| {
            for l in &self.log_tail {
                let color = if l.contains("ERROR") {
                    RED
                } else if l.contains("WARN") {
                    AMBER
                } else {
                    ui.visuals().text_color()
                };
                ui.label(RichText::new(l).monospace().small().color(color));
            }
        });
    }
}

// ---------------------------------------------------------------- pieces

fn prelim_status_label(ui: &mut egui::Ui, b: &BarcodeState, summary: Option<&PrelimSummary>) {
    let text = match (b.prelim_status.as_deref(), summary) {
        (Some("queued"), _) => RichText::new("queued").color(AMBER),
        (Some("ok"), Some(s)) => RichText::new(format!("ok · gen {}", s.generation)).color(GREEN),
        (Some("partial"), _) => RichText::new("partial").color(AMBER),
        (Some("failed"), _) => RichText::new("failed").color(RED),
        (_, Some(s)) => RichText::new(format!("gen {}", s.generation)).color(GREY),
        _ => RichText::new("-").color(GREY),
    };
    ui.label(text);
}

fn fusion_label(ui: &mut egui::Ui, s: Option<&PrelimSummary>) {
    match s.and_then(|s| s.fusion.as_ref()) {
        Some(f) => match &f.primary_isoform {
            Some(iso) => ui.label(RichText::new(format!("{iso} ({} reads)", data::format_count(f.total_fusion_reads))).strong()),
            None => ui.label(RichText::new("none called").color(GREY)),
        },
        None => ui.label("-"),
    };
}

fn mutation_label(ui: &mut egui::Ui, s: Option<&PrelimSummary>) {
    let Some(b) = s.and_then(|s| s.branford.as_ref()) else {
        ui.label("-");
        return;
    };
    if b.mutations.is_empty() {
        ui.label(RichText::new("none detected").color(GREEN));
    } else {
        let names: Vec<String> = b.mutations.iter().map(|m| format!("{} {:.1}%", m.mutation, m.vaf_percent)).collect();
        ui.label(RichText::new(names.join(", ")).color(RED));
    }
}

fn enough_label(ui: &mut egui::Ui, s: Option<&PrelimSummary>) {
    match s.and_then(|s| s.sufficiency.as_ref()) {
        Some(x) if x.enough => ui.label(RichText::new("yes").color(GREEN)),
        Some(x) => ui.label(RichText::new(format!("not yet ({:.0}% of sites)", 100.0 * x.site_coverage_fraction)).color(AMBER)),
        None => ui.label("-"),
    };
}

fn summary_header(ui: &mut egui::Ui, s: &PrelimSummary, b: Option<&BarcodeState>) {
    ui.horizontal_wrapped(|ui| {
        ui.heading(&s.barcode);
        ui.label(format!(
            "snapshot generation {} · {} reads passed the filter · analysed {} ago in {:.1} s",
            s.generation,
            data::format_count(s.reads_passed),
            data::format_age(now_unix().saturating_sub(s.written_unix)),
            s.duration_secs
        ));
        let behind = b.map(|b| b.compacted_generation > s.generation).unwrap_or(false);
        if behind {
            ui.label(RichText::new("· a newer snapshot is being analysed").color(AMBER));
        }
    });
    ui.label(
        RichText::new("Preliminary: computed on the reads sequenced so far, with the same scripts the final nagini run uses.")
            .small()
            .color(GREY),
    );
    for e in &s.errors {
        ui.label(RichText::new(e).color(RED));
    }
}

pub(super) fn section_title(ui: &mut egui::Ui, title: &str) {
    ui.label(RichText::new(title).strong().size(16.0));
}

pub(super) fn fusion_section(ui: &mut egui::Ui, f: Option<&FusionSummary>, changed: bool) {
    section_title(ui, "Fusion isoform");
    let Some(f) = f else {
        ui.label("Not available.");
        return;
    };
    ui.horizontal_wrapped(|ui| {
        match &f.primary_isoform {
            Some(iso) => ui.label(RichText::new(format!("Primary: {iso}")).strong().size(18.0)),
            None => ui.label(RichText::new("No isoform has reached the calling threshold yet").color(AMBER)),
        };
        ui.label(format!("{} fusion-spanning reads", data::format_count(f.total_fusion_reads)));
        if changed {
            ui.label(RichText::new("The leading isoform changed between snapshots: the call is not stable yet.").color(AMBER));
        }
    });
    egui::Grid::new("bcr_isoforms").striped(true).spacing([14.0, 4.0]).show(ui, |ui| {
        for h in ["Isoform", "Protein", "Reads", "Fraction", "", "Status"] {
            ui.label(RichText::new(h).strong());
        }
        ui.end_row();
        for r in f.isoforms.iter().filter(|r| r.read_count > 0 || r.status == "PRIMARY") {
            ui.label(&r.isoform);
            ui.label(&r.protein_form);
            ui.label(data::format_count(r.read_count));
            ui.label(format!("{:.1}%", r.fraction_percent));
            ui.add(egui::ProgressBar::new((r.fraction_percent / 100.0) as f32).desired_width(120.0));
            let color = match r.status.as_str() {
                "PRIMARY" => GREEN,
                "SECONDARY" => AMBER,
                _ => GREY,
            };
            ui.label(RichText::new(&r.status).color(color));
            ui.end_row();
        }
    });
}

fn tki_color(class: &str) -> Color32 {
    match class {
        "Resistant" => RED,
        "Intermediate" => AMBER,
        "Sensitive" => GREEN,
        _ => GREY,
    }
}

pub(super) fn mutation_section(ui: &mut egui::Ui, b: Option<&BranfordSummary>) {
    section_title(ui, "ABL1 kinase-domain mutations");
    let Some(b) = b else {
        ui.label("Not available.");
        return;
    };
    if b.mutations.is_empty() {
        ui.label(RichText::new("No mutation above nagini's calling thresholds.").color(GREEN));
        return;
    }
    egui::Grid::new("bcr_mutations").striped(true).spacing([12.0, 4.0]).show(ui, |ui| {
        for h in ["Mutation", "Change", "Position", "Depth", "Alt reads", "VAF", "Status", "Imat", "Dasa", "Nilo", "Bosu", "Pona", "Asci"] {
            ui.label(RichText::new(h).strong());
        }
        ui.end_row();
        for m in &b.mutations {
            ui.label(RichText::new(&m.mutation).strong());
            ui.label(&m.aa_change);
            ui.label(format!("{}:{}", m.chrom, m.pos));
            ui.label(data::format_count(m.depth));
            ui.label(data::format_count(m.ad_alt));
            ui.label(format!("{:.2}%", m.vaf_percent));
            ui.label(RichText::new(&m.status).color(if m.status == "DETECTED" { RED } else { AMBER }));
            for tki in ["imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib", "asciminib"] {
                let class = m.tki.get(tki).map(String::as_str).unwrap_or("");
                ui.label(RichText::new(match class { "Resistant" => "R", "Intermediate" => "I", "Sensitive" => "S", _ => "?" }).color(tki_color(class)).strong())
                    .on_hover_text(format!("{tki}: {}", if class.is_empty() { "unknown" } else { class }));
            }
            ui.end_row();
        }
    });
}

pub(super) fn phasing_section(ui: &mut egui::Ui, p: Option<&Vec<PhasePair>>) {
    section_title(ui, "Compound mutations (single-molecule phasing)");
    match p {
        None => {
            ui.label("Not available.");
        }
        Some(p) if p.is_empty() => {
            ui.label("No pair of detected mutations to phase.");
        }
        Some(p) => {
            egui::Grid::new("bcr_phasing").striped(true).spacing([14.0, 4.0]).show(ui, |ui| {
                for h in ["Mutation 1", "Mutation 2", "State", "In cis", "Spanning reads"] {
                    ui.label(RichText::new(h).strong());
                }
                ui.end_row();
                for x in p {
                    ui.label(&x.mutation_1);
                    ui.label(&x.mutation_2);
                    let cis = x.phase_state.starts_with("CIS");
                    ui.label(RichText::new(&x.phase_state).color(if cis { RED } else { AMBER }));
                    ui.label(format!("{:.1}%", x.cis_fraction_percent));
                    ui.label(data::format_count(x.total_spanning_reads));
                    ui.end_row();
                }
            });
        }
    }
}

fn coverage_section(ui: &mut egui::Ui, s: &PrelimSummary) {
    section_title(ui, "Depth over the kinase-domain sites");
    let Some(b) = &s.branford else {
        ui.label("Not available.");
        return;
    };
    ui.label(format!(
        "{} of {} sites at ≥ {}× · lowest {}× · median {}×",
        b.sites_at_min_depth,
        b.sites_total,
        b.min_depth_setting,
        data::format_count(b.min_site_depth),
        data::format_count(b.median_site_depth)
    ));
    if let Some(x) = &s.sufficiency {
        let (text, color) = if x.enough {
            ("Enough data for a reading: a fusion is called and ≥ 95% of sites are covered.".to_string(), GREEN)
        } else {
            (
                format!(
                    "Not enough yet: fusion {} · {:.0}% of sites covered (needs {:.0}%).",
                    if x.fusion_called { "called" } else { "not called" },
                    100.0 * x.site_coverage_fraction,
                    100.0 * x.required_site_coverage_fraction
                ),
                AMBER,
            )
        };
        ui.label(RichText::new(text).color(color));
    }
}

fn history_plots(ui: &mut egui::Ui, history: &[data::HistoryPoint]) {
    section_title(ui, "How the result firms up with more reads");
    if history.len() < 2 {
        ui.label("Needs at least two snapshots.");
        return;
    }
    ui.label("Mutation VAF (%) against reads that passed the filter");
    Plot::new("bcr_vaf_history")
        .height(200.0)
        .legend(Legend::default())
        .x_axis_label("reads passed")
        .y_axis_label("VAF %")
        .show(ui, |plot_ui| {
            for (name, pts) in data::mutation_series(history) {
                plot_ui.line(Line::new(PlotPoints::from(pts)).name(name));
            }
        });
    ui.add_space(6.0);
    ui.label("Fusion-spanning reads against reads that passed the filter");
    Plot::new("bcr_fusion_history")
        .height(160.0)
        .x_axis_label("reads passed")
        .y_axis_label("fusion reads")
        .show(ui, |plot_ui| {
            plot_ui.line(Line::new(PlotPoints::from(data::fusion_series(history))).name("fusion reads"));
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanobcr_core::EngineConfig;
    use nanobcr_core::state::State;
    use std::fs;
    use std::path::Path;

    /// Lay out a finished engine directory with a realistic preliminary result.
    fn fixture(dir: &Path) {
        let cfg = EngineConfig::new(dir.into(), dir.into(), PathBuf::from("/ref.fa"));
        let mut st = State::new(&cfg);
        st.status = "running".into();
        let b = st.barcode("barcode24");
        b.files = 8;
        b.reads_seen = 31_110;
        b.reads_passed = 2_742;
        b.primary_mapped = 2_742;
        b.supplementary = 2_274;
        b.compacted_generation = 3;
        b.last_update_unix = now_unix();
        b.prelim_status = Some("ok".into());
        st.write(&dir.join("state.json")).unwrap();

        let p = dir.join("barcode24/prelim");
        fs::create_dir_all(&p).unwrap();
        let summary = serde_json::json!({
            "schema": 1, "barcode": "barcode24", "generation": 3, "reads_passed": 2742,
            "written_unix": now_unix(), "duration_secs": 0.7, "status": "ok", "errors": [],
            "fusion": {"total_fusion_reads": 1607, "primary_isoform": "e1a2", "isoforms": [
                {"isoform":"e1a2","protein_form":"p190","read_count":1566,"fraction_percent":97.45,"status":"PRIMARY"},
                {"isoform":"e13a2","protein_form":"p210","read_count":0,"fraction_percent":0.0,"status":"NOT_DETECTED"}]},
            "branford": {"sites_total": 73, "sites_at_min_depth": 73, "min_depth_setting": 50, "min_site_depth": 1446, "median_site_depth": 1572,
                "mutations": [{"mutation":"T315I","aa_change":"p.Thr315Ile","chrom":"chr9","pos":133748283,"depth":1618,"ad_alt":38,"vaf_percent":2.349,"status":"DETECTED",
                    "tki":{"imatinib":"Resistant","dasatinib":"Resistant","nilotinib":"Resistant","bosutinib":"Resistant","ponatinib":"Sensitive","asciminib":"Sensitive"}}]},
            "phasing": [{"mutation_1":"Y253H","mutation_2":"T315I","phase_state":"CIS (Compound)","cis_fraction_percent":90.2,"total_spanning_reads":1585}],
            "sufficiency": {"fusion_called": true, "site_coverage_fraction": 1.0, "required_site_coverage_fraction": 0.95, "enough": true},
            "dir": "x"
        });
        fs::write(p.join("summary.json"), summary.to_string()).unwrap();
        let hist = [
            r#"{"unix":1,"generation":1,"reads_passed":329,"status":"ok","primary_isoform":"e1a2","fusion_reads":220,"mutations":[],"enough":true}"#,
            r#"{"unix":2,"generation":2,"reads_passed":1869,"status":"ok","primary_isoform":"e1a2","fusion_reads":1598,"mutations":[{"mutation":"T315I","vaf_percent":2.36,"depth":900,"status":"DETECTED"}],"enough":true}"#,
            r#"{"unix":3,"generation":3,"reads_passed":2742,"status":"ok","primary_isoform":"e1a2","fusion_reads":1607,"mutations":[{"mutation":"T315I","vaf_percent":2.35,"depth":1618,"status":"DETECTED"}],"enough":true}"#,
        ];
        fs::write(p.join("history.jsonl"), hist.join("\n")).unwrap();
        fs::write(dir.join("engine.log"), "[INFO] cycle 1\n[WARN] x\n[ERROR] y\n").unwrap();
    }

    fn run_frames(app: &mut BcrApp, n: usize) {
        let ctx = egui::Context::default();
        for _ in 0..n {
            let _ = ctx.run(egui::RawInput::default(), |ctx| app.show(ctx));
        }
    }

    #[test]
    fn every_tab_renders_without_panicking_on_real_looking_data() {
        let tmp = tempfile::tempdir().unwrap();
        fixture(tmp.path());
        let setup = BcrSetup { out_dir: tmp.path().to_string_lossy().into(), ..Default::default() };
        let mut app = BcrApp::new(BcrStartup { setup });
        assert_eq!(app.selected.as_deref(), Some("barcode24"));
        for tab in [Tab::Samples, Tab::Detail, Tab::Final, Tab::Log] {
            app.tab = tab;
            run_frames(&mut app, 3);
        }
    }

    #[test]
    fn empty_and_missing_output_directories_render_too() {
        let tmp = tempfile::tempdir().unwrap();
        for out in [String::new(), tmp.path().to_string_lossy().into(), "/no/such/dir".to_string()] {
            let mut app = BcrApp::new(BcrStartup { setup: BcrSetup { out_dir: out, ..Default::default() } });
            for tab in [Tab::Samples, Tab::Detail, Tab::Final, Tab::Log] {
                app.tab = tab;
                run_frames(&mut app, 2);
            }
        }
    }

    #[test]
    fn a_broken_summary_is_shown_as_a_problem_not_a_crash() {
        let tmp = tempfile::tempdir().unwrap();
        fixture(tmp.path());
        fs::write(tmp.path().join("barcode24/prelim/summary.json"), "{").unwrap();
        let mut app = BcrApp::new(BcrStartup { setup: BcrSetup { out_dir: tmp.path().to_string_lossy().into(), ..Default::default() } });
        app.tab = Tab::Detail;
        run_frames(&mut app, 2);
        assert!(!app.data.problems.is_empty());
    }

    /// A finished nagini run next to the engine output, with the same T315I/e1a2 result.
    fn finished_run(dir: &Path, same: bool) {
        use nanobcr_core::handoff::{Handoff, HandoffBarcode};
        let hand = dir.join("handoff");
        fs::create_dir_all(&hand).unwrap();
        let h = Handoff {
            dir: hand.clone(),
            barcode_map: hand.join("barcode_map.csv"),
            written_unix: now_unix(),
            barcodes: vec![HandoffBarcode {
                barcode: "barcode24".into(),
                sample: "PAT".into(),
                sample_type: "T".into(),
                bam: hand.join("barcode24.bam"),
                snapshot_generation: 3,
                primary_records: 2742,
                reads_passed: 2742,
            }],
            skipped: vec![("barcode99".into(), "no sample in the barcode map".into())],
            warnings: vec!["barcode09 (GHOST) is in the barcode map but the engine has no reads for it".into()],
        };
        fs::write(hand.join("handoff.json"), serde_json::to_vec(&h).unwrap()).unwrap();
        fs::write(hand.join("provenance.json"), r#"{"minimap2":"2.31","alignment":{"preset":"splice","extra":["-y","-Y","--MD"]},"read_filter":{"min_qv":20.0,"min_len":500,"max_len":3000}}"#).unwrap();

        let a = dir.join("nagini_run");
        let rep = a.join("results/reports/bcr_abl");
        fs::create_dir_all(&rep).unwrap();
        let iso = if same { "e1a2" } else { "e14a2" };
        fs::write(rep.join("PAT_T.fusion_summary.json"), format!(r#"{{"total_fusion_reads":1607,"primary_isoform":"{iso}","isoforms":[{{"fusion_isoform":"{iso}","protein_form":"p190","read_count":1566,"fraction_percent":97.4,"status":"PRIMARY"}}]}}"#)).unwrap();
        fs::write(rep.join("PAT_T.branford_all_sites.tsv"), "sample\tchrom\tpos\tdepth\nx\tchr9\t1\t900\n").unwrap();
        fs::write(rep.join("PAT_T.branford_mutations.tsv"), "sample\tchrom\tpos\tmutation\taa_change\tdepth\tad_alt\tvaf_percent\tstatus\timatinib\nx\tchr9\t133748283\tT315I\tp.Thr315Ile\t1618\t38\t2.349\tDETECTED\tResistant\n").unwrap();
        fs::write(rep.join("PAT_T.phasing.json"), r#"{"phased_pairs":[{"mutation_1":"Y253H","mutation_2":"T315I","phase_state":"CIS (Compound)","cis_fraction_percent":90.2,"total_spanning_reads":1585}]}"#).unwrap();
        fs::write(rep.join("PAT_T.exon_coverage.tsv"), "sample\texon\tmean_depth\tpct_gte_100x\tspanning_reads\nx\tABL1_e2\t1626.1\t100.0\t1625\n").unwrap();
        // A real (tiny) PNG so the texture path is exercised.
        image::RgbaImage::from_pixel(8, 4, image::Rgba([10, 120, 200, 255])).save(rep.join("PAT_T.fusion_plot.png")).unwrap();
        fs::write(a.join("results/reports/nagini_summary_report.xlsx"), "x").unwrap();
        let rec = nanobcr_core::pipeline::RunRecord {
            pid: 1, started_unix: now_unix() - 600, command: vec![], analysis_dir: a.clone(), handoff_dir: hand.clone(), barcode_map: hand.join("barcode_map.csv"),
        };
        fs::write(nanobcr_core::pipeline::run_record_path(&a), serde_json::to_vec(&rec).unwrap()).unwrap();
        fs::write(nanobcr_core::pipeline::exit_file(&a), "0\n").unwrap();
        fs::create_dir_all(a.join(".snakemake/log")).unwrap();
        fs::write(a.join(".snakemake/log/x.snakemake.log"), "rule one:\n8 of 8 steps (100%) done\n").unwrap();
    }

    #[test]
    fn the_final_tab_shows_nagini_results_and_an_identical_comparison() {
        let tmp = tempfile::tempdir().unwrap();
        fixture(tmp.path());
        finished_run(tmp.path(), true);
        let setup = BcrSetup { out_dir: tmp.path().to_string_lossy().into(), nagini_dir: "/some/nagini".into(), ..Default::default() };
        let mut app = BcrApp::new(BcrStartup { setup });
        app.tab = Tab::Final;
        run_frames(&mut app, 3);
        let f = &app.final_run;
        assert_eq!(f.results_len(), 1);
        assert_eq!(f.comparison_identical("barcode24"), Some(true));
    }

    #[test]
    fn a_difference_between_preliminary_and_final_is_surfaced() {
        let tmp = tempfile::tempdir().unwrap();
        fixture(tmp.path());
        finished_run(tmp.path(), false);
        let setup = BcrSetup { out_dir: tmp.path().to_string_lossy().into(), nagini_dir: "/some/nagini".into(), ..Default::default() };
        let mut app = BcrApp::new(BcrStartup { setup });
        app.tab = Tab::Final;
        run_frames(&mut app, 3);
        assert_eq!(app.final_run.comparison_identical("barcode24"), Some(false));
    }
}
