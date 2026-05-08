use crate::model::{
    AnalysisMode, DashboardData, FilterConfig, HistogramBin, MainTab, ResultRow, RunConfig,
    RunSource,
};
use crate::nanoparse_cli::NanoparseConfig;
use crate::remote::{MonitorEvent, MonitorRequest, RemoteConfig, RemoteStatus};
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke, Vec2};
use egui_plot::{Bar, BarChart, Line, Plot, PlotPoints, Points};
use nanofilter_core::{bam as filter_bam, fastq as filter_fastq};
use nanoparse_core::{MatchMode, matcher, matcher::AmpliconResult};
use nanoseq_core::filters::parse_pore_range;
use rfd::FileDialog;
use std::collections::{HashSet, VecDeque};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;
use std::time::Duration;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RunState {
    Idle,
    Running,
}

enum WorkerMessage {
    Log(String),
    Completed {
        input_path: String,
        result: Result<DashboardData, String>,
    },
}

enum MonitorMessage {
    Discovered(String),
    Log(String),
}

enum FileOpMessage {
    Completed(Result<String, String>),
}

#[derive(Debug, Clone, Default)]
pub struct AppStartupConfig {
    pub mode: Option<AnalysisMode>,
    pub input_path: Option<String>,
    pub monitor_dir: Option<String>,
    pub reference_path: Option<String>,
    pub gtf_path: Option<String>,
    pub primers_path: Option<String>,
    pub nanoparse_bin: Option<String>,
    pub run_on_start: bool,
}

pub struct NanoMonitorApp {
    mode: AnalysisMode,
    tab: MainTab,
    filters: FilterConfig,
    run: RunConfig,
    remote: RemoteConfig,
    nanoparse: NanoparseConfig,
    reference_path: String,
    gtf_path: String,
    data: DashboardData,
    log_lines: Vec<String>,
    run_state: RunState,
    worker_rx: Option<Receiver<WorkerMessage>>,
    monitor_rx: Option<Receiver<MonitorMessage>>,
    monitor_stop_tx: Option<Sender<()>>,
    monitor_active: bool,
    queued_files: HashSet<String>,
    processed_files: HashSet<String>,
    failed_files: HashSet<String>,
    pending_files: VecDeque<String>,
    current_input: Option<String>,
    last_error: Option<String>,
    file_op_rx: Option<Receiver<FileOpMessage>>,
    file_op_running: bool,
    file_op_output_path: String,
    file_op_max_len: u32,
    file_op_channel_range: String,
    barcode_file_path: String,
    barcode_output_dir: String,
}

impl NanoMonitorApp {
    fn is_supported_data_file(path: &Path) -> bool {
        let s = path.to_string_lossy().to_ascii_lowercase();
        s.ends_with(".bam")
            || s.ends_with(".fastq")
            || s.ends_with(".fq")
            || s.ends_with(".fastq.gz")
            || s.ends_with(".fq.gz")
    }

    fn collect_data_files(dir: &Path, out: &mut Vec<PathBuf>) {
        if let Ok(entries) = fs::read_dir(dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    Self::collect_data_files(&path, out);
                } else if Self::is_supported_data_file(&path) {
                    out.push(path);
                }
            }
        }
    }

    fn resolve_run_input_path(&mut self) -> Result<String, String> {
        match self.run.source {
            RunSource::SingleFile => {
                let p = PathBuf::from(self.run.input_path.trim());
                if !p.exists() {
                    return Err("Input path does not exist".into());
                }
                if !p.is_file() {
                    return Err("Input path is not a file; use Monitor Dir for folders".into());
                }
                if !Self::is_supported_data_file(&p) {
                    return Err("Unsupported input file. Expected BAM/FASTQ(.gz)".into());
                }
                Ok(p.to_string_lossy().to_string())
            }
            RunSource::MonitorDirectory => {
                let dir = PathBuf::from(self.run.monitor_dir.trim());
                if !dir.exists() || !dir.is_dir() {
                    return Err("Monitor directory is missing or invalid".into());
                }
                let mut files = Vec::new();
                Self::collect_data_files(&dir, &mut files);
                if files.is_empty() {
                    return Err("No BAM/FASTQ files found in monitor directory".into());
                }
                files.sort_by(|a, b| {
                    let ma = fs::metadata(a).and_then(|m| m.modified()).ok();
                    let mb = fs::metadata(b).and_then(|m| m.modified()).ok();
                    mb.cmp(&ma)
                });
                let selected = files[0].to_string_lossy().to_string();
                self.log_lines.push(format!(
                    "Monitor dir: found {} supported files, using latest: {}",
                    files.len(),
                    selected
                ));
                Ok(selected)
            }
        }
    }

    fn list_monitor_files_sorted(&self) -> Result<Vec<String>, String> {
        let dir = PathBuf::from(self.run.monitor_dir.trim());
        if !dir.exists() || !dir.is_dir() {
            return Err("Monitor directory is missing or invalid".into());
        }
        let mut files = Vec::new();
        Self::collect_data_files(&dir, &mut files);
        files.sort_by(|a, b| {
            let ma = fs::metadata(a).and_then(|m| m.modified()).ok();
            let mb = fs::metadata(b).and_then(|m| m.modified()).ok();
            ma.cmp(&mb)
        });
        Ok(files
            .into_iter()
            .map(|p| p.to_string_lossy().to_string())
            .collect())
    }

    fn validate_run_prerequisites(&mut self) -> Result<(), String> {
        if self.mode != AnalysisMode::Amplicon {
            return Err("Only Amplicon mode is wired to nanoparse for now".into());
        }
        if !self.filters.use_nanoparse {
            return Err("Enable 'Use Rust (nanoparse)' to run analysis".into());
        }
        if self.nanoparse.primers_path.trim().is_empty() {
            return Err("Primers path is required".into());
        }
        let primer_path = PathBuf::from(self.nanoparse.primers_path.trim());
        if !primer_path.exists() {
            return Err(format!(
                "Primers path does not exist: {}",
                self.nanoparse.primers_path
            ));
        }
        if !self.reference_path.trim().is_empty() {
            let reference = PathBuf::from(self.reference_path.trim());
            if !reference.exists() {
                return Err(format!(
                    "Reference path does not exist: {}",
                    self.reference_path
                ));
            }
        }
        if !self.gtf_path.trim().is_empty() {
            let gtf = PathBuf::from(self.gtf_path.trim());
            if !gtf.exists() {
                return Err(format!("GTF/BED path does not exist: {}", self.gtf_path));
            }
        }
        Ok(())
    }

    fn find_in_path(executable: &str) -> Option<PathBuf> {
        let path_var = std::env::var_os("PATH")?;
        for dir in std::env::split_paths(&path_var) {
            let candidate = dir.join(executable);
            if candidate.is_file() {
                return Some(candidate);
            }
            #[cfg(windows)]
            {
                let candidate_exe = dir.join(format!("{}.exe", executable));
                if candidate_exe.is_file() {
                    return Some(candidate_exe);
                }
            }
        }
        None
    }

    fn resolve_nanoparse_executable(&self) -> Result<String, String> {
        let configured = self.nanoparse.executable.trim();
        if configured.is_empty() {
            return Err("nanoparse binary path is empty".into());
        }

        let configured_path = PathBuf::from(configured);
        if configured_path.is_absolute() || configured.contains('/') || configured.contains('\\') {
            if configured_path.is_file() {
                return Ok(configured_path.to_string_lossy().to_string());
            }
            return Err(format!("nanoparse binary not found: {}", configured));
        }

        if let Some(found) = Self::find_in_path(configured) {
            return Ok(found.to_string_lossy().to_string());
        }

        let exe_name = {
            #[cfg(windows)]
            {
                "nanoparse.exe"
            }
            #[cfg(not(windows))]
            {
                "nanoparse"
            }
        };

        let cwd = std::env::current_dir().ok();
        let mut candidates = Vec::new();
        if let Some(cwd) = cwd {
            candidates.push(cwd.join("target").join("debug").join(exe_name));
            candidates.push(cwd.join("target").join("release").join(exe_name));
            candidates.push(
                cwd.join("nanoparse")
                    .join("target")
                    .join("debug")
                    .join(exe_name),
            );
            candidates.push(
                cwd.join("nanoparse")
                    .join("target")
                    .join("release")
                    .join(exe_name),
            );
        }

        if let Ok(current_exe) = std::env::current_exe() {
            if let Some(bin_dir) = current_exe.parent() {
                candidates.push(bin_dir.join(exe_name));
                if let Some(target_dir) = bin_dir.parent() {
                    candidates.push(target_dir.join("debug").join(exe_name));
                    candidates.push(target_dir.join("release").join(exe_name));
                }
            }
        }

        for candidate in candidates {
            if candidate.is_file() {
                return Ok(candidate.to_string_lossy().to_string());
            }
        }

        Err(format!(
            "Could not locate 'nanoparse'. Set nanoparse binary explicitly or build it (`cargo build -p nanoparse`)."
        ))
    }

    fn enqueue_file(&mut self, path: String) {
        if self.processed_files.contains(&path)
            || self.failed_files.contains(&path)
            || self.queued_files.contains(&path)
            || self.current_input.as_deref() == Some(path.as_str())
        {
            return;
        }
        self.queued_files.insert(path.clone());
        self.pending_files.push_back(path.clone());
        self.log_lines.push(format!("Queued file: {}", path));
    }

    fn maybe_start_next_analysis(&mut self) {
        if self.run_state == RunState::Running {
            return;
        }
        if let Some(next) = self.pending_files.pop_front() {
            self.queued_files.remove(&next);
            self.start_analysis_for_file(next);
        }
    }

    fn start_directory_watcher(&mut self, initial_known: Vec<String>) -> Result<(), String> {
        if self.monitor_active {
            return Ok(());
        }
        let dir = PathBuf::from(self.run.monitor_dir.trim());
        if !dir.exists() || !dir.is_dir() {
            return Err("Monitor directory is missing or invalid".into());
        }

        let (tx, rx) = mpsc::channel::<MonitorMessage>();
        let (stop_tx, stop_rx) = mpsc::channel::<()>();
        let dir_for_thread = dir.clone();
        let mut known: HashSet<String> = initial_known.into_iter().collect();

        thread::spawn(move || {
            loop {
                if stop_rx.try_recv().is_ok() {
                    break;
                }

                let mut files = Vec::new();
                NanoMonitorApp::collect_data_files(&dir_for_thread, &mut files);
                for p in files {
                    let s = p.to_string_lossy().to_string();
                    if known.insert(s.clone()) {
                        let _ = tx.send(MonitorMessage::Discovered(s));
                    }
                }
                thread::sleep(Duration::from_secs(2));
            }
            let _ = tx.send(MonitorMessage::Log("Directory watcher stopped".into()));
        });

        self.monitor_rx = Some(rx);
        self.monitor_stop_tx = Some(stop_tx);
        self.monitor_active = true;
        self.log_lines.push(format!(
            "Monitoring directory for new files: {}",
            dir.to_string_lossy()
        ));
        Ok(())
    }

    fn stop_directory_watcher(&mut self) {
        if let Some(stop_tx) = self.monitor_stop_tx.take() {
            let _ = stop_tx.send(());
        }
        self.monitor_rx = None;
        self.monitor_active = false;
    }

    fn nanoparse_supports_file(path: &str) -> bool {
        let s = path.to_ascii_lowercase();
        s.ends_with(".bam")
            || s.ends_with(".fastq")
            || s.ends_with(".fq")
            || s.ends_with(".fastq.gz")
            || s.ends_with(".fq.gz")
    }

    fn start_analysis_for_file(&mut self, input_path: String) {
        if !Self::nanoparse_supports_file(&input_path) {
            let msg = format!(
                "Skipping unsupported input for nanoparse (expected BAM/FASTQ): {}",
                input_path
            );
            self.last_error = Some(msg.clone());
            self.log_lines.push(msg);
            self.failed_files.insert(input_path);
            return;
        }

        let (tx, rx) = mpsc::channel::<WorkerMessage>();
        let command_line = self
            .nanoparse
            .build_amplicon_command(
                &input_path,
                &self.filters,
                Some(self.reference_path.as_str()),
                Some(self.gtf_path.as_str()),
            )
            .as_shell_line();
        self.log_lines
            .push(format!("Running in-process analysis for {}", input_path));
        self.log_lines
            .push(format!("CLI preview> {}", command_line));
        self.worker_rx = Some(rx);
        self.run_state = RunState::Running;
        self.current_input = Some(input_path.clone());
        let primers_path = self.nanoparse.primers_path.clone();
        let threads = self.nanoparse.threads;
        let primer_tolerance = self.nanoparse.primer_tolerance;
        let min_qs = self.filters.min_qs;
        let min_len = self.filters.min_len as usize;
        let max_reads = self.filters.max_reads as usize;
        let duplex_only = self.filters.duplex_only;
        let reference =
            (!self.reference_path.trim().is_empty()).then(|| self.reference_path.clone());
        let gtf = (!self.gtf_path.trim().is_empty()).then(|| self.gtf_path.clone());

        thread::spawn(move || {
            let _ = tx.send(WorkerMessage::Log(format!("CLI preview> {}", command_line)));
            let result = match matcher::run_amplicons(
                &input_path,
                &primers_path,
                threads,
                MatchMode::Semiglobal,
                3,
                150,
                true,
                primer_tolerance,
                min_qs,
                min_len,
                max_reads,
                duplex_only,
                reference.as_deref(),
                gtf.as_deref(),
            ) {
                Ok(parsed) => Ok(build_dashboard_from_nanoparse(parsed)),
                Err(e) => Err(format!("nanoparse core failed: {}", e)),
            };

            let _ = tx.send(WorkerMessage::Completed { input_path, result });
        });
    }

    fn pick_input_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("BAM/FASTQ", &["bam", "fastq", "fq", "gz"])
            .pick_file()
        {
            self.run.source = RunSource::SingleFile;
            self.run.input_path = path.to_string_lossy().to_string();
        }
    }

    fn pick_monitor_dir(&mut self) {
        if let Some(path) = FileDialog::new().pick_folder() {
            self.run.source = RunSource::MonitorDirectory;
            self.run.monitor_dir = path.to_string_lossy().to_string();
        }
    }

    fn pick_primers_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("TSV/TXT", &["tsv", "txt"])
            .pick_file()
        {
            self.nanoparse.primers_path = path.to_string_lossy().to_string();
        }
    }

    fn pick_filter_output_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .set_file_name("filtered_output.bam")
            .save_file()
        {
            self.file_op_output_path = path.to_string_lossy().to_string();
        }
    }

    fn pick_barcodes_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("TXT", &["txt", "tsv"])
            .pick_file()
        {
            self.barcode_file_path = path.to_string_lossy().to_string();
        }
    }

    fn pick_barcode_output_dir(&mut self) {
        if let Some(path) = FileDialog::new().pick_folder() {
            self.barcode_output_dir = path.to_string_lossy().to_string();
        }
    }

    fn selected_file_for_operations(&self) -> Result<String, String> {
        match self.run.source {
            RunSource::SingleFile => {
                let p = PathBuf::from(self.run.input_path.trim());
                if !p.exists() || !p.is_file() {
                    return Err("Select a valid input file first".into());
                }
                Ok(p.to_string_lossy().to_string())
            }
            RunSource::MonitorDirectory => self
                .current_input
                .clone()
                .ok_or_else(|| "File operations currently require an active file".into()),
        }
    }

    fn start_filter_export(&mut self) {
        if self.file_op_running {
            return;
        }
        let input = match self.selected_file_for_operations() {
            Ok(v) => v,
            Err(msg) => {
                self.last_error = Some(msg.clone());
                self.log_lines.push(msg);
                return;
            }
        };
        if self.file_op_output_path.trim().is_empty() {
            let msg = "Output path is required for filter/export".to_string();
            self.last_error = Some(msg.clone());
            self.log_lines.push(msg);
            return;
        }

        let settings = match self.build_file_filter_settings() {
            Ok(v) => v,
            Err(msg) => {
                self.last_error = Some(msg.clone());
                self.log_lines.push(msg);
                return;
            }
        };

        let output = self.file_op_output_path.clone();
        let (tx, rx) = mpsc::channel::<FileOpMessage>();
        self.file_op_rx = Some(rx);
        self.file_op_running = true;
        self.log_lines
            .push(format!("Filtering/exporting {}", input));

        thread::spawn(move || {
            let input_lower = input.to_ascii_lowercase();
            let output_lower = output.to_ascii_lowercase();
            let output_bam = output_lower.ends_with(".bam");
            let result = if input_lower.ends_with(".bam") {
                if output_bam {
                    filter_bam::filter_bam_with_settings(
                        Path::new(&input),
                        Path::new(&output),
                        &settings,
                    )
                    .map(|_| format!("Wrote filtered BAM to {}", output))
                } else {
                    filter_bam::bam_to_fastq_with_settings(
                        Path::new(&input),
                        Path::new(&output),
                        &settings,
                    )
                    .map(|_| format!("Wrote filtered FASTQ to {}", output))
                }
            } else if Self::nanoparse_supports_file(&input) {
                filter_fastq::filter_fastq_with_settings(
                    Path::new(&input),
                    Path::new(&output),
                    &settings,
                    output_bam,
                )
                .map(|_| format!("Wrote filtered output to {}", output))
            } else {
                Err(anyhow::anyhow!("Unsupported input format"))
            };

            let _ = tx.send(FileOpMessage::Completed(
                result.map_err(|e| format!("Filter/export failed: {}", e)),
            ));
        });
    }

    fn start_barcode_split(&mut self) {
        if self.file_op_running {
            return;
        }
        let input = match self.selected_file_for_operations() {
            Ok(v) => v,
            Err(msg) => {
                self.last_error = Some(msg.clone());
                self.log_lines.push(msg);
                return;
            }
        };
        if self.barcode_file_path.trim().is_empty() {
            let msg = "Barcode file is required".to_string();
            self.last_error = Some(msg.clone());
            self.log_lines.push(msg);
            return;
        }
        if self.barcode_output_dir.trim().is_empty() {
            let msg = "Barcode output directory is required".to_string();
            self.last_error = Some(msg.clone());
            self.log_lines.push(msg);
            return;
        }

        let barcode_file = self.barcode_file_path.clone();
        let output_dir = self.barcode_output_dir.clone();
        let (tx, rx) = mpsc::channel::<FileOpMessage>();
        self.file_op_rx = Some(rx);
        self.file_op_running = true;
        self.log_lines
            .push(format!("Splitting {} by barcodes", input));

        thread::spawn(move || {
            let input_lower = input.to_ascii_lowercase();
            let result = if input_lower.ends_with(".bam") {
                filter_bam::split_bam_by_barcodes(
                    Path::new(&input),
                    Path::new(&barcode_file),
                    Path::new(&output_dir),
                    1,
                    1000,
                    1,
                    0,
                    false,
                )
                .map(|_| format!("Barcode split complete in {}", output_dir))
            } else if Self::nanoparse_supports_file(&input) {
                filter_fastq::split_fastq_by_barcodes(
                    Path::new(&input),
                    Path::new(&barcode_file),
                    Path::new(&output_dir),
                    1,
                    1000,
                    1,
                    0,
                    false,
                )
                .map(|_| format!("Barcode split complete in {}", output_dir))
            } else {
                Err(anyhow::anyhow!("Unsupported input format"))
            };

            let _ = tx.send(FileOpMessage::Completed(
                result.map_err(|e| format!("Barcode split failed: {}", e)),
            ));
        });
    }

    fn build_file_filter_settings(&self) -> Result<filter_bam::FilterSettings, String> {
        let channel_range = if self.file_op_channel_range.trim().is_empty() {
            None
        } else {
            Some(parse_pore_range(self.file_op_channel_range.trim()).map_err(|e| e.to_string())?)
        };
        Ok(filter_bam::FilterSettings {
            qv_threshold: self.filters.min_qs as f64,
            min_len: self.filters.min_len as usize,
            max_len: if self.file_op_max_len == 0 {
                usize::MAX
            } else {
                self.file_op_max_len as usize
            },
            channel_range,
            time_window: None,
        })
    }

    fn pick_gtf_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("GTF/GFF/BED", &["gtf", "gff", "gff3", "bed", "gz"])
            .pick_file()
        {
            self.gtf_path = path.to_string_lossy().to_string();
        }
    }

    fn pick_reference_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("FASTA", &["fa", "fasta", "fna", "fa.gz", "fasta.gz"])
            .pick_file()
        {
            self.reference_path = path.to_string_lossy().to_string();
        }
    }

    pub fn new(cc: &eframe::CreationContext<'_>, startup: AppStartupConfig) -> Self {
        cc.egui_ctx.set_visuals(egui::Visuals::light());
        let mut app = Self {
            mode: startup.mode.unwrap_or(AnalysisMode::Amplicon),
            tab: MainTab::Results,
            filters: FilterConfig::default(),
            run: RunConfig::default(),
            remote: RemoteConfig::default(),
            nanoparse: NanoparseConfig::default(),
            reference_path: String::new(),
            gtf_path: String::new(),
            data: DashboardData::empty(),
            log_lines: vec![
                "nanoMonitor initialized".into(),
                "Ready for local and remote analysis".into(),
                "Press Start Monitor to run nanoparse and refresh distributions".into(),
            ],
            run_state: RunState::Idle,
            worker_rx: None,
            monitor_rx: None,
            monitor_stop_tx: None,
            monitor_active: false,
            queued_files: HashSet::new(),
            processed_files: HashSet::new(),
            failed_files: HashSet::new(),
            pending_files: VecDeque::new(),
            current_input: None,
            last_error: None,
            file_op_rx: None,
            file_op_running: false,
            file_op_output_path: String::new(),
            file_op_max_len: 0,
            file_op_channel_range: String::new(),
            barcode_file_path: String::new(),
            barcode_output_dir: String::new(),
        };

        if let Some(bin) = startup.nanoparse_bin {
            app.nanoparse.executable = bin;
        }
        if let Some(primers) = startup.primers_path {
            app.nanoparse.primers_path = primers;
        }
        if let Some(reference) = startup.reference_path {
            app.reference_path = reference;
        }
        if let Some(gtf) = startup.gtf_path {
            app.gtf_path = gtf;
        }
        if let Some(input) = startup.input_path {
            let p = PathBuf::from(&input);
            if p.is_dir() {
                app.run.source = RunSource::MonitorDirectory;
                app.run.monitor_dir = input;
            } else {
                app.run.source = RunSource::SingleFile;
                app.run.input_path = input;
            }
        }
        if let Some(dir) = startup.monitor_dir {
            app.run.source = RunSource::MonitorDirectory;
            app.run.monitor_dir = dir;
        }
        if startup.run_on_start {
            app.start_monitor();
        }
        app
    }

    fn start_monitor(&mut self) {
        if self.run_state == RunState::Running || self.monitor_active {
            return;
        }
        self.last_error = None;
        self.data = DashboardData::empty();
        self.pending_files.clear();
        self.queued_files.clear();
        self.processed_files.clear();
        self.failed_files.clear();
        if let Err(msg) = self.validate_run_prerequisites() {
            self.last_error = Some(msg.clone());
            self.log_lines.push(msg);
            return;
        }

        match self.run.source {
            RunSource::SingleFile => match self.resolve_run_input_path() {
                Ok(path) => self.start_analysis_for_file(path),
                Err(msg) => {
                    self.last_error = Some(msg.clone());
                    self.log_lines.push(msg);
                }
            },
            RunSource::MonitorDirectory => {
                let initial_files = match self.list_monitor_files_sorted() {
                    Ok(files) => files,
                    Err(msg) => {
                        self.last_error = Some(msg.clone());
                        self.log_lines.push(msg);
                        return;
                    }
                };
                if let Err(msg) = self.start_directory_watcher(initial_files.clone()) {
                    self.last_error = Some(msg.clone());
                    self.log_lines.push(msg);
                    return;
                }
                for file in initial_files {
                    self.enqueue_file(file);
                }
                self.maybe_start_next_analysis();
            }
        }
    }

    fn stop_monitor(&mut self) {
        self.stop_directory_watcher();
        self.pending_files.clear();
        self.queued_files.clear();
        self.log_lines.push("Monitoring stopped".into());
        if self.run_state == RunState::Running {
            self.log_lines.push(
                "Analysis is running in-process; stop will prevent new runs after current file completes"
                    .into(),
            );
        }
    }

    fn poll_worker(&mut self, ctx: &egui::Context) {
        let mut file_op_msgs = Vec::new();
        if let Some(rx) = &self.file_op_rx {
            while let Ok(msg) = rx.try_recv() {
                file_op_msgs.push(msg);
            }
        }
        for msg in file_op_msgs {
            match msg {
                FileOpMessage::Completed(result) => {
                    self.file_op_running = false;
                    match result {
                        Ok(line) => self.log_lines.push(line),
                        Err(err) => {
                            self.last_error = Some(err.clone());
                            self.log_lines.push(err);
                        }
                    }
                }
            }
        }
        if !self.file_op_running {
            self.file_op_rx = None;
        }

        let mut monitor_msgs = Vec::new();
        if let Some(rx) = &self.monitor_rx {
            while let Ok(msg) = rx.try_recv() {
                monitor_msgs.push(msg);
            }
        }
        for msg in monitor_msgs {
            match msg {
                MonitorMessage::Discovered(path) => self.enqueue_file(path),
                MonitorMessage::Log(msg) => self.log_lines.push(msg),
            }
        }

        let mut completed = false;
        let mut completed_path: Option<String> = None;
        let mut completed_ok = false;
        if let Some(rx) = &self.worker_rx {
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    WorkerMessage::Log(line) => {
                        for chunk in line.lines() {
                            self.log_lines.push(chunk.to_string());
                        }
                    }
                    WorkerMessage::Completed { input_path, result } => {
                        completed = true;
                        completed_path = Some(input_path.clone());
                        self.run_state = RunState::Idle;
                        match result {
                            Ok(data) => {
                                self.data = data;
                                completed_ok = true;
                                self.last_error = None;
                                self.log_lines
                                    .push(format!("Analysis complete: {}", input_path));
                                if self.run.auto_scan_variants {
                                    self.log_lines.push(format!(
                                        "Auto-variant hook queued for {} (implementation pending)",
                                        input_path
                                    ));
                                }
                            }
                            Err(err) => {
                                self.last_error = Some(err.clone());
                                self.log_lines
                                    .push(format!("Analysis failed for {}: {}", input_path, err));
                            }
                        }
                    }
                }
            }
        }
        if completed {
            self.worker_rx = None;
            self.current_input = None;
            if let Some(path) = completed_path {
                if completed_ok {
                    self.processed_files.insert(path);
                } else {
                    self.failed_files.insert(path);
                }
            }
            self.maybe_start_next_analysis();
        } else if self.run_state == RunState::Running {
            ctx.request_repaint();
        } else if self.monitor_active && !self.pending_files.is_empty() {
            self.maybe_start_next_analysis();
            ctx.request_repaint();
        }
    }

    fn draw_top_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.heading("nanoMonitor");
            ui.label(
                RichText::new("Rust / egui prototype")
                    .color(Color32::from_rgb(70, 70, 70))
                    .small(),
            );
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let run_label = if self.monitor_active {
                    "Run: Monitoring"
                } else {
                    match self.run_state {
                        RunState::Idle => "Run: Idle",
                        RunState::Running => "Run: Running",
                    }
                };
                ui.colored_label(
                    if self.run_state == RunState::Running || self.monitor_active {
                        Color32::from_rgb(210, 140, 20)
                    } else {
                        Color32::from_rgb(70, 70, 70)
                    },
                    run_label,
                );
                ui.separator();
                let (color, label) = match self.remote.status {
                    RemoteStatus::Connected => {
                        (Color32::from_rgb(34, 139, 34), "Remote: Connected")
                    }
                    RemoteStatus::Connecting => {
                        (Color32::from_rgb(210, 140, 20), "Remote: Connecting")
                    }
                    RemoteStatus::Disconnected => {
                        (Color32::from_rgb(160, 40, 40), "Remote: Offline")
                    }
                };
                ui.colored_label(color, label);
            });
        });
        ui.separator();
    }

    fn draw_sidebar(&mut self, ui: &mut egui::Ui) {
        ui.spacing_mut().item_spacing = Vec2::new(8.0, 8.0);

        ui.group(|ui| {
            ui.label(RichText::new("Mode").strong());
            egui::ComboBox::from_id_salt("mode_combo")
                .selected_text(self.mode.label())
                .show_ui(ui, |ui| {
                    for mode in AnalysisMode::ALL {
                        ui.selectable_value(&mut self.mode, mode, mode.label());
                    }
                });
        });

        ui.group(|ui| {
            ui.label(RichText::new("Resources").strong());
            if ui.button("Load Primers").clicked() {
                self.pick_primers_file();
            }
            if ui.button("Load GTF/BED").clicked() {
                self.pick_gtf_file();
            }
            if ui.button("Load Reference FASTA").clicked() {
                self.pick_reference_file();
            }
            ui.horizontal(|ui| {
                ui.label("nanoparse binary:");
                ui.text_edit_singleline(&mut self.nanoparse.executable);
            });
            ui.horizontal(|ui| {
                ui.label("primers file:");
                ui.text_edit_singleline(&mut self.nanoparse.primers_path);
                if ui.button("...").clicked() {
                    self.pick_primers_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("reference:");
                ui.text_edit_singleline(&mut self.reference_path);
                if ui.button("...").clicked() {
                    self.pick_reference_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("gtf/bed:");
                ui.text_edit_singleline(&mut self.gtf_path);
                if ui.button("...").clicked() {
                    self.pick_gtf_file();
                }
            });
        });

        ui.group(|ui| {
            ui.label(RichText::new("Run Control").strong());
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.run.source, RunSource::SingleFile, "Single File");
                ui.selectable_value(
                    &mut self.run.source,
                    RunSource::MonitorDirectory,
                    "Monitor Dir",
                );
            });
            ui.label("Input path");
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut self.run.input_path);
                if ui.button("Browse").clicked() {
                    self.pick_input_file();
                }
            });
            ui.label("Monitor directory");
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut self.run.monitor_dir);
                if ui.button("Browse").clicked() {
                    self.pick_monitor_dir();
                }
            });
            ui.checkbox(&mut self.run.auto_scan_variants, "Auto-scan variants");
            ui.label(format!(
                "Queue: {} | Processed: {} | Failed: {}",
                self.pending_files.len(),
                self.processed_files.len(),
                self.failed_files.len()
            ));
            if let Some(err) = &self.last_error {
                ui.colored_label(
                    Color32::from_rgb(170, 35, 35),
                    format!("Last error: {}", err),
                );
            }
            if let Some(current) = &self.current_input {
                ui.label(format!("Current: {}", current));
            }
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(
                        self.run_state == RunState::Idle && !self.monitor_active,
                        egui::Button::new("Start Monitor"),
                    )
                    .clicked()
                {
                    self.start_monitor();
                }
                if ui
                    .add_enabled(
                        self.run_state == RunState::Running || self.monitor_active,
                        egui::Button::new("Stop"),
                    )
                    .clicked()
                {
                    self.stop_monitor();
                }
            });
        });

        ui.group(|ui| {
            ui.label(RichText::new("File Ops").strong());
            ui.label("Filter/export output");
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut self.file_op_output_path);
                if ui.button("Browse").clicked() {
                    self.pick_filter_output_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("Channel range");
                ui.text_edit_singleline(&mut self.file_op_channel_range);
            });
            ui.horizontal(|ui| {
                ui.label("Max len");
                ui.add(egui::DragValue::new(&mut self.file_op_max_len).speed(10.0));
                ui.label("0 = no max");
            });
            if ui
                .add_enabled(!self.file_op_running, egui::Button::new("Filter / Extract"))
                .clicked()
            {
                self.start_filter_export();
            }

            ui.separator();
            ui.label("Barcode file");
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut self.barcode_file_path);
                if ui.button("Browse").clicked() {
                    self.pick_barcodes_file();
                }
            });
            ui.label("Barcode output dir");
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut self.barcode_output_dir);
                if ui.button("Browse").clicked() {
                    self.pick_barcode_output_dir();
                }
            });
            if ui
                .add_enabled(!self.file_op_running, egui::Button::new("Split Barcodes"))
                .clicked()
            {
                self.start_barcode_split();
            }
        });

        ui.group(|ui| {
            ui.label(RichText::new("Remote (Secondary)").strong());
            ui.checkbox(&mut self.remote.enabled, "Enable remote server");
            ui.label("Endpoint");
            ui.text_edit_singleline(&mut self.remote.endpoint);
            ui.label("Token");
            ui.add(egui::TextEdit::singleline(&mut self.remote.auth_token).password(true));
            ui.horizontal(|ui| {
                if ui.button("Connect").clicked() {
                    self.remote.status = RemoteStatus::Connecting;
                    self.log_lines
                        .push(format!("Connecting to {}", self.remote.endpoint));
                    if let Ok(msg) = serde_json::to_string(&MonitorRequest::Ping) {
                        self.log_lines.push(format!("Protocol example -> {}", msg));
                    }
                    self.remote.status = RemoteStatus::Connected;
                    let ack = MonitorEvent::Pong;
                    if let Ok(msg) = serde_json::to_string(&ack) {
                        self.log_lines
                            .push(format!("Remote event example <- {}", msg));
                    }
                }
                if ui.button("Disconnect").clicked() {
                    self.remote.status = RemoteStatus::Disconnected;
                    self.log_lines.push("Disconnected".into());
                }
            });
            ui.label(format!("Status: {}", self.remote.status.label()));
        });

        ui.group(|ui| {
            ui.label(RichText::new("Quick Actions").strong());
            ui.horizontal(|ui| {
                let _ = ui.button("Snap");
                let _ = ui.button("Variants");
                let _ = ui.button("Matrix");
            });
        });
    }

    fn draw_filter_strip(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Interactive Filters").strong());
                ui.add(egui::DragValue::new(&mut self.filters.min_qs).speed(0.1).prefix("Min QS "));
                ui.add(egui::DragValue::new(&mut self.filters.min_len).speed(10.0).prefix("Min Len "));
                ui.add(egui::DragValue::new(&mut self.filters.max_reads).speed(100.0).prefix("Max Reads "));
                ui.label(RichText::new("(0 = all reads)").small());
                ui.checkbox(&mut self.filters.duplex_only, "Duplex only");
                ui.checkbox(&mut self.filters.use_nanoparse, "Use Rust (nanoparse)");
                if ui.button("Recalculate").clicked() {
                    self.log_lines
                        .push("Recalculate requested with current filters".into());
                }
                if ui.button("Auto-Variant").clicked() {
                    self.run.auto_scan_variants = true;
                    self.log_lines.push(
                        "Auto-Variant enabled. Variant pipeline hook will run after successful analysis."
                            .into(),
                    );
                }
                if ui.button("Run Duplex Discovery").clicked() {
                    self.log_lines.push(
                        "Duplex discovery action requested (backend hook pending)".into(),
                    );
                }
            });
        });
    }

    fn draw_result_table(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.tab, MainTab::Results, "Results");
            ui.selectable_value(&mut self.tab, MainTab::Log, "Log");
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                if ui.button("Build nanoparse command").clicked() {
                    match self.resolve_run_input_path() {
                        Ok(path) => {
                            let mut cmd = self.nanoparse.build_amplicon_command(
                                &path,
                                &self.filters,
                                Some(self.reference_path.as_str()),
                                Some(self.gtf_path.as_str()),
                            );
                            match self.resolve_nanoparse_executable() {
                                Ok(bin) => {
                                    cmd.program = bin;
                                    self.log_lines.push(format!("CLI> {}", cmd.as_shell_line()));
                                }
                                Err(msg) => {
                                    self.last_error = Some(msg.clone());
                                    self.log_lines
                                        .push(format!("Cannot build command: {}", msg));
                                }
                            }
                        }
                        Err(msg) => self
                            .log_lines
                            .push(format!("Cannot build command: {}", msg)),
                    }
                }
            });
        });
        ui.separator();

        match self.tab {
            MainTab::Results => {
                egui::ScrollArea::vertical()
                    .max_height(340.0)
                    .show(ui, |ui| {
                        let stroke = Stroke::new(1.0, Color32::from_gray(210));
                        egui::Frame::default().stroke(stroke).show(ui, |ui| {
                            ui.horizontal(|ui| {
                                ui.add_sized(
                                    [28.0, 18.0],
                                    egui::Label::new(RichText::new("#").strong()),
                                );
                                ui.add_sized(
                                    [560.0, 18.0],
                                    egui::Label::new(RichText::new("Amplicon Name").strong()),
                                );
                                ui.add_sized(
                                    [80.0, 18.0],
                                    egui::Label::new(RichText::new("Count").strong()),
                                );
                                ui.add_sized(
                                    [90.0, 18.0],
                                    egui::Label::new(RichText::new("Med Len").strong()),
                                );
                                ui.add_sized(
                                    [80.0, 18.0],
                                    egui::Label::new(RichText::new("SD Len").strong()),
                                );
                                ui.add_sized(
                                    [80.0, 18.0],
                                    egui::Label::new(RichText::new("Avg QS").strong()),
                                );
                                ui.add_sized(
                                    [80.0, 18.0],
                                    egui::Label::new(RichText::new("Vars").strong()),
                                );
                            });
                            ui.separator();

                            if self.data.rows.is_empty() {
                                ui.label(
                                    RichText::new(
                                        "No results yet. Select input and press Start Monitor.",
                                    )
                                    .italics(),
                                );
                            } else {
                                for (idx, row) in self.data.rows.iter().enumerate() {
                                    let selected = self.data.selected_row == Some(idx);
                                    ui.horizontal(|ui| {
                                        if ui
                                            .selectable_label(selected, format!("{}", idx + 1))
                                            .clicked()
                                        {
                                            self.data.selected_row = Some(idx);
                                        }
                                        let label =
                                            ui.selectable_label(selected, &row.amplicon_name);
                                        if label.clicked() {
                                            self.data.selected_row = Some(idx);
                                        }
                                        ui.add_sized(
                                            [80.0, 18.0],
                                            egui::Label::new(format!("{}", row.count)),
                                        );
                                        ui.add_sized(
                                            [90.0, 18.0],
                                            egui::Label::new(format!("{}", row.median_length)),
                                        );
                                        ui.add_sized(
                                            [80.0, 18.0],
                                            egui::Label::new(format!("{:.1}", row.sd_length)),
                                        );
                                        ui.add_sized(
                                            [80.0, 18.0],
                                            egui::Label::new(format!("{:.1}", row.avg_qs)),
                                        );
                                        ui.add_sized(
                                            [80.0, 18.0],
                                            egui::Label::new(format!("{}", row.variants)),
                                        );
                                    });
                                }
                            }
                        });
                    });
            }
            MainTab::Log => {
                egui::ScrollArea::vertical()
                    .max_height(340.0)
                    .show(ui, |ui| {
                        for line in self.log_lines.iter().rev().take(120) {
                            ui.monospace(line);
                        }
                    });
            }
        }

        let (sel_amplicons, sel_reads) = self.data.selected_counts();
        ui.separator();
        ui.horizontal(|ui| {
            ui.label(format!(
                "Selected: {} amplicons, {} reads",
                sel_amplicons, sel_reads
            ));
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let _ = ui.add_enabled(false, egui::Button::new("Export Selected"));
            });
        });
    }

    fn draw_bottom_plots(&mut self, ui: &mut egui::Ui) {
        ui.columns(3, |columns| {
            columns[0].group(|ui| {
                ui.label(RichText::new("Accuracy Density").strong());
                if self.data.accuracy_bins.is_empty() {
                    ui.label(RichText::new("No data").italics());
                } else {
                    ui.label(
                        RichText::new(format!("Mode: {:.2}", self.data.accuracy_mode))
                            .color(Color32::from_rgb(180, 40, 40)),
                    );
                }
                let line = Line::new(PlotPoints::from_iter(density_points(
                    &self.data.accuracy_bins,
                )));
                Plot::new("acc_density")
                    .allow_drag(true)
                    .allow_zoom(true)
                    .allow_scroll(true)
                    .height(220.0)
                    .show(ui, |plot_ui| {
                        if !self.data.accuracy_bins.is_empty() {
                            plot_ui.line(line.color(Color32::from_rgb(76, 175, 80)));
                        }
                    });
            });

            columns[1].group(|ui| {
                ui.label(RichText::new("Q-Score Density").strong());
                if self.data.qs_bins.is_empty() {
                    ui.label(RichText::new("No data").italics());
                } else {
                    ui.label(
                        RichText::new(format!("Mode: {:.2}", self.data.qs_mode))
                            .color(Color32::from_rgb(180, 40, 40)),
                    );
                }
                let line = Line::new(PlotPoints::from_iter(density_points(&self.data.qs_bins)));
                Plot::new("qs_density")
                    .allow_drag(true)
                    .allow_zoom(true)
                    .allow_scroll(true)
                    .height(220.0)
                    .show(ui, |plot_ui| {
                        if !self.data.qs_bins.is_empty() {
                            plot_ui.line(line.color(Color32::from_rgb(30, 136, 229)));
                        }
                    });
            });

            columns[2].group(|ui| {
                if self.mode == AnalysisMode::Wgs {
                    ui.label(RichText::new("CNV (WGS)").strong());
                    let pts = self
                        .data
                        .cnv_bins
                        .iter()
                        .map(|p| [p.position_mb, p.log2_ratio])
                        .collect::<Vec<_>>();
                    Plot::new("cnv_plot")
                        .height(220.0)
                        .allow_drag(true)
                        .allow_zoom(true)
                        .allow_scroll(true)
                        .show(ui, |plot_ui| {
                            if pts.is_empty() {
                                // No data loaded yet for WGS.
                            } else {
                                plot_ui.points(
                                    Points::new(pts)
                                        .radius(1.5)
                                        .color(Color32::from_rgb(136, 84, 208)),
                                );
                            }
                        });
                } else {
                    ui.label(RichText::new("Length Histogram").strong());
                    if self.data.length_bins.is_empty() {
                        ui.label(RichText::new("No data").italics());
                    } else {
                        ui.label(
                            RichText::new(format!("Median: {:.0}", self.data.length_median))
                                .color(Color32::from_rgb(230, 140, 0)),
                        );
                    }
                    let bars = self
                        .data
                        .length_bins
                        .iter()
                        .map(|b| {
                            let center = (b.start + b.end) * 0.5;
                            Bar::new(center, b.count as f64).width((b.end - b.start).max(1.0))
                        })
                        .collect::<Vec<_>>();
                    Plot::new("len_hist")
                        .allow_drag(true)
                        .allow_zoom(true)
                        .allow_scroll(true)
                        .height(220.0)
                        .show(ui, |plot_ui| {
                            if !self.data.length_bins.is_empty() {
                                plot_ui.bar_chart(
                                    BarChart::new(bars).color(Color32::from_rgb(126, 87, 194)),
                                );
                            }
                        });
                }
            });
        });
    }
}

impl eframe::App for NanoMonitorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_worker(ctx);
        if self.monitor_active {
            ctx.request_repaint_after(Duration::from_millis(300));
        }

        egui::TopBottomPanel::top("top_panel")
            .resizable(false)
            .show(ctx, |ui| self.draw_top_bar(ui));

        egui::SidePanel::left("left_sidebar")
            .resizable(true)
            .default_width(340.0)
            .show(ctx, |ui| self.draw_sidebar(ui));

        egui::CentralPanel::default().show(ctx, |ui| {
            self.draw_filter_strip(ui);
            ui.add_space(8.0);
            self.draw_result_table(ui);
            ui.add_space(8.0);
            self.draw_bottom_plots(ui);
            ui.add_space(6.0);
            ui.label(format!(
                "Total: {} | Filtered: {}",
                format_count(self.data.total_reads),
                format_count(self.data.filtered_reads)
            ));
        });
    }
}

impl Drop for NanoMonitorApp {
    fn drop(&mut self) {
        self.stop_directory_watcher();
    }
}

fn build_dashboard_from_nanoparse(output: AmpliconResult) -> DashboardData {
    let mut rows = output
        .amplicons
        .into_iter()
        .map(|(name, s)| ResultRow {
            amplicon_name: name,
            count: s.count as u32,
            median_length: s.median_length as u32,
            sd_length: s.std_length as f32,
            avg_qs: s.avg_qs,
            variants: 0,
        })
        .collect::<Vec<_>>();
    rows.sort_by(|a, b| b.count.cmp(&a.count));

    let d = output.distributions;
    let (length_bins, qs_bins, accuracy_bins, length_median, qs_mode, accuracy_mode) = (
        map_bins(d.length_bins),
        map_bins(d.qs_bins),
        map_bins(d.accuracy_bins),
        d.length_median,
        d.qs_mode,
        d.accuracy_mode,
    );

    let filtered_reads = output.total_reads.saturating_sub(output.unmatched_count);

    DashboardData {
        rows,
        selected_row: None,
        total_reads: output.total_reads as u64,
        filtered_reads: filtered_reads as u64,
        cnv_bins: Vec::new(),
        length_bins,
        qs_bins,
        accuracy_bins,
        length_median,
        qs_mode,
        accuracy_mode,
    }
}

fn map_bins(bins: Vec<matcher::DistributionBin>) -> Vec<HistogramBin> {
    bins.into_iter()
        .map(|b| HistogramBin {
            start: b.start,
            end: b.end,
            count: b.count as u64,
        })
        .collect()
}

fn density_points(bins: &[HistogramBin]) -> Vec<[f64; 2]> {
    let total: f64 = bins.iter().map(|b| b.count as f64).sum();
    if total <= 0.0 {
        return Vec::new();
    }
    bins.iter()
        .filter(|b| b.count > 0)
        .map(|b| {
            let center = (b.start + b.end) * 0.5;
            let width = (b.end - b.start).max(1e-9);
            let density = (b.count as f64) / total / width;
            [center, density]
        })
        .collect()
}

fn format_count(value: u64) -> String {
    let text = value.to_string();
    let mut out = String::with_capacity(text.len() + text.len() / 3);
    for (i, ch) in text.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}
