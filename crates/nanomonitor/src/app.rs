use crate::integrations::{self, ClinicalMutation};
use crate::model::{
    AnalysisMode, DashboardData, FilterConfig, HistogramBin, MainTab, ResultRow, RunConfig,
    RunSource, ThemeSkin, VariantRow,
};
use crate::nanostream_cli::NanostreamConfig;
use crate::remote::{self, MonitorEvent, MonitorRequest, RemoteConfig, RemoteStatus};
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke, Vec2};
use egui_plot::{Bar, BarChart, Line, Plot, PlotPoints};
use nanofilter_core::{bam as filter_bam, fastq as filter_fastq};
use nanoparse_core::{MatchMode, matcher};
use nanoseq_core::filters::parse_pore_range;
use rfd::FileDialog;
use std::collections::{HashMap, HashSet, VecDeque};
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
    Progress {
        input_path: String,
        processed_reads: usize,
        current_data: DashboardData,
    },
    Completed {
        input_path: String,
        result: Result<DashboardData, String>,
    },
    SnapCompleted {
        path: String,
    },
    VariantsCompleted {
        variants: Vec<VariantRow>,
    },
    Error(String),
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
    pub nanostream_bin: Option<String>,
    pub run_on_start: bool,
}

pub struct NanoMonitorApp {
    mode: AnalysisMode,
    tab: MainTab,
    filters: FilterConfig,
    run: RunConfig,
    remote: RemoteConfig,
    pub nanostream: NanostreamConfig,
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

    // Premium Aesthetics & Integrations
    active_theme: ThemeSkin,
    theme_applied: bool,
    clinical_db: HashMap<(String, u64), ClinicalMutation>,
    tool_rs_qc: String,
    tool_rindels: String,
    remote_req_tx: Option<Sender<MonitorRequest>>,
    remote_stop_tx: Option<Sender<()>>,
    remote_rx: Option<Receiver<MonitorEvent>>,
    
    // Background pipelines state
    snap_running: bool,
    variant_calling_running: bool,
    active_center_sub_tab: usize, // 0 = Histograms, 1 = Alignment Snapshot, 2 = Called Variants

    accuracy_plot_hovered: bool,
    qs_plot_hovered: bool,
    len_plot_hovered: bool,
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
        if !self.filters.use_nanostream {
            return Err("Enable 'Use Rust (nanostream)' to run analysis".into());
        }
        if self.nanostream.primers_path.trim().is_empty() {
            return Err("Primers path is required".into());
        }
        let primer_path = PathBuf::from(self.nanostream.primers_path.trim());
        if !primer_path.exists() {
            return Err(format!(
                "Primers path does not exist: {}",
                self.nanostream.primers_path
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

    fn resolve_nanostream_executable(&self) -> Result<String, String> {
        let configured = self.nanostream.executable.trim();
        if configured.is_empty() {
            return Err("nanostream binary path is empty".into());
        }

        let configured_path = PathBuf::from(configured);
        if configured_path.is_absolute() || configured.contains('/') || configured.contains('\\') {
            if configured_path.is_file() {
                return Ok(configured_path.to_string_lossy().to_string());
            }
            return Err(format!("nanostream binary not found: {}", configured));
        }

        if let Some(found) = Self::find_in_path(configured) {
            return Ok(found.to_string_lossy().to_string());
        }

        let exe_name = {
            #[cfg(windows)]
            {
                "nanostream.exe"
            }
            #[cfg(not(windows))]
            {
                "nanostream"
            }
        };

        let cwd = std::env::current_dir().ok();
        let mut candidates = Vec::new();
        if let Some(cwd) = cwd {
            candidates.push(cwd.join("target").join("debug").join(exe_name));
            candidates.push(cwd.join("target").join("release").join(exe_name));
            candidates.push(
                cwd.join("nanostream")
                    .join("target")
                    .join("debug")
                    .join(exe_name),
            );
            candidates.push(
                cwd.join("nanostream")
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
            "Could not locate 'nanostream'. Set nanostream binary explicitly or build it (`cargo build -p nanostream`)."
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
                Self::collect_data_files(&dir_for_thread, &mut files);
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

    fn nanostream_supports_file(path: &str) -> bool {
        let s = path.to_ascii_lowercase();
        s.ends_with(".bam")
            || s.ends_with(".fastq")
            || s.ends_with(".fq")
            || s.ends_with(".fastq.gz")
            || s.ends_with(".fq.gz")
    }

    fn start_analysis_for_file(&mut self, input_path: String) {
        if !Self::nanostream_supports_file(&input_path) {
            let msg = format!(
                "Skipping unsupported input for nanostream (expected BAM/FASTQ): {}",
                input_path
            );
            self.last_error = Some(msg.clone());
            self.log_lines.push(msg);
            self.failed_files.insert(input_path);
            return;
        }

        let (tx, rx) = mpsc::channel::<WorkerMessage>();
        let command_line = self
            .nanostream
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
        let primers_path = self.nanostream.primers_path.clone();
        let threads = self.nanostream.threads;
        let primer_tolerance = self.nanostream.primer_tolerance;
        let min_qs = self.filters.min_qs;
        let min_len = self.filters.min_len as usize;
        let max_reads = self.filters.max_reads as usize;
        let duplex_only = self.filters.duplex_only;
        let reference =
            (!self.reference_path.trim().is_empty()).then(|| self.reference_path.clone());
        let gtf = (!self.gtf_path.trim().is_empty()).then(|| self.gtf_path.clone());

        let tx_progress = tx.clone();
        let input_path_clone = input_path.clone();

        thread::spawn(move || {
            let _ = tx.send(WorkerMessage::Log(format!("CLI preview> {}", command_line)));
            let result = match matcher::run_amplicons_with_callback(
                &input_path,
                &primers_path,
                threads,
                MatchMode::Semiglobal,
                3,
                150,
                true,
                primer_tolerance,
                min_qs,
                (min_len, usize::MAX),
                max_reads,
                duplex_only,
                reference.as_deref(),
                gtf.as_deref(),
                None,
                None,
                false,
                false,
                |intermediate| {

                    let new_data = build_dashboard_from_nanostream(intermediate);
                    let _ = tx_progress.send(WorkerMessage::Progress {
                        input_path: input_path_clone.clone(),
                        processed_reads: intermediate.total_reads,
                        current_data: new_data,
                    });
                }
            ) {
                Ok(parsed) => Ok(build_dashboard_from_nanostream(&parsed)),
                Err(e) => Err(format!("nanostream core failed: {}", e)),
            };

            let _ = tx.send(WorkerMessage::Completed { input_path, result });
        });
    }

    fn trigger_region_snapshot(&mut self) {
        let row_idx = match self.data.selected_row {
            Some(i) => i,
            None => {
                self.log_lines.push("Select an amplicon row in the table first".into());
                return;
            }
        };
        let region = self.data.rows[row_idx].amplicon_name.clone();
        let input_bam = match self.selected_file_for_operations() {
            Ok(path) => path,
            Err(e) => {
                self.log_lines.push(format!("Cannot generate snapshot: {}", e));
                return;
            }
        };

        if !input_bam.to_ascii_lowercase().ends_with(".bam") {
            self.log_lines.push("Snapshot generation requires a loaded BAM file".into());
            return;
        }

        self.snap_running = true;
        self.log_lines.push(format!("Queued rs-qc snapshot for amplicon: {}", region));
        
        // We can spawn a quick thread sending results back to worker channel if mapped, or log it.
        // Let's create a dynamic channel.
        let worker_tx = self.create_worker_sender_if_none();

        let rs_qc_bin = self.tool_rs_qc.clone();
        let gtf = self.gtf_path.clone();
        let reference = self.reference_path.clone();
        
        // Output folder snapshot setup
        let output_folder = Path::new("snapshots");
        let _ = fs::create_dir_all(output_folder);
        let clean_name = region.replace(':', "_").replace('-', "_");
        let output_png = output_folder.join(format!("{}_snapshot.png", clean_name));
        let output_png_str = output_png.to_string_lossy().to_string();

        thread::spawn(move || {
            match integrations::run_rs_qc_snap(&rs_qc_bin, &input_bam, &region, &gtf, &reference, &output_png_str) {
                Ok(_) => {
                    let _ = worker_tx.send(WorkerMessage::SnapCompleted { path: output_png_str });
                }
                Err(e) => {
                    let _ = worker_tx.send(WorkerMessage::Error(format!("rs-qc Snap failed: {}", e)));
                }
            }
        });
    }

    fn trigger_variant_calling(&mut self) {
        let row_idx = match self.data.selected_row {
            Some(i) => i,
            None => {
                self.log_lines.push("Select an amplicon row in the table first".into());
                return;
            }
        };
        let region = self.data.rows[row_idx].amplicon_name.clone();
        let input_bam = match self.selected_file_for_operations() {
            Ok(path) => path,
            Err(e) => {
                self.log_lines.push(format!("Cannot call variants: {}", e));
                return;
            }
        };

        if !input_bam.to_ascii_lowercase().ends_with(".bam") {
            self.log_lines.push("Variant calling requires a mapped BAM file".into());
            return;
        }

        self.variant_calling_running = true;
        self.log_lines.push(format!("Queued rindels variant calling for: {}", region));

        let worker_tx = self.create_worker_sender_if_none();
        
        let rindels_bin = self.tool_rindels.clone();
        let reference = self.reference_path.clone();
        let clinical_db = self.clinical_db.clone();

        thread::spawn(move || {
            // Write temporary BED file containing the region if parsed, or format a quick bed file
            let bed_path = "temp_amplicon.bed";
            // Simple VCF output setup
            let output_vcf = "temp_variants.vcf";
            
            // Format Region to BED format chrom \t start \t end
            let mut parts = region.split(':');
            let bed_content = if let Some(chrom) = parts.next() {
                if let Some(range) = parts.next() {
                    let mut coords = range.split('-');
                    let start = coords.next().unwrap_or("0");
                    let end = coords.next().unwrap_or("0");
                    format!("{}\t{}\t{}\n", chrom, start, end)
                } else {
                    format!("{}\t0\t500000000\n", chrom)
                }
            } else {
                format!("{}\t0\t500000000\n", region)
            };
            let _ = fs::write(bed_path, bed_content);

            match integrations::run_rindels(&rindels_bin, &input_bam, &reference, bed_path, output_vcf) {
                Ok(_) => {
                    // Parse VCF
                    match integrations::parse_vcf(output_vcf, &clinical_db) {
                        Ok(vars) => {
                            let _ = worker_tx.send(WorkerMessage::VariantsCompleted { variants: vars });
                        }
                        Err(e) => {
                            let _ = worker_tx.send(WorkerMessage::Error(format!("VCF Parsing failed: {}", e)));
                        }
                    }
                }
                Err(e) => {
                    let _ = worker_tx.send(WorkerMessage::Error(format!("rindels calling failed: {}", e)));
                }
            }
            
            // Cleanup temp files
            let _ = fs::remove_file(bed_path);
            let _ = fs::remove_file(output_vcf);
        });
    }

    fn create_worker_sender_if_none(&mut self) -> Sender<WorkerMessage> {
        let (tx, rx) = mpsc::channel::<WorkerMessage>();
        // Add rx to worker channel checks
        self.worker_rx = Some(rx);
        tx
    }

    fn pick_input_file(&mut self) {
        if let Some(path) = FileDialog::new()
            .add_filter("BAM/FASTQ", &["bam", "fastq", "fq", "gz"])
            .pick_file()
        {
            self.run.source = RunSource::SingleFile;
            self.run.input_path = path.to_string_lossy().to_string();
            self.data = DashboardData::empty();
            self.processed_files.clear();
            self.failed_files.clear();
            self.pending_files.clear();
            self.queued_files.clear();
            self.current_input = None;
            self.last_error = None;
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
            self.nanostream.primers_path = path.to_string_lossy().to_string();
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
            } else if Self::nanostream_supports_file(&input) {
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
            } else if Self::nanostream_supports_file(&input) {
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
        // Bio-Teal (Dark) is the default stunning aesthetic
        let default_theme = ThemeSkin::BioTeal;

        let mut app = Self {
            mode: startup.mode.unwrap_or(AnalysisMode::Amplicon),
            tab: MainTab::Results,
            filters: FilterConfig::default(),
            run: RunConfig::default(),
            remote: RemoteConfig::default(),
            nanostream: NanostreamConfig::default(),
            reference_path: String::new(),
            gtf_path: String::new(),
            data: DashboardData::empty(),
            log_lines: vec![
                "nanoMonitor PCR Suite initialized".into(),
                "Ready for local and remote amplicon analysis".into(),
                "Select data inputs and press Start Monitor.".into(),
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

            // Aesthetic themes and integrations defaults
            active_theme: default_theme,
            theme_applied: false,
            clinical_db: integrations::load_clinical_mutations(),
            tool_rs_qc: String::new(),
            tool_rindels: String::new(),
            remote_req_tx: None,
            remote_stop_tx: None,
            remote_rx: None,
            snap_running: false,
            variant_calling_running: false,
            active_center_sub_tab: 0,
            accuracy_plot_hovered: false,
            qs_plot_hovered: false,
            len_plot_hovered: false,
        };

        if let Some(bin) = startup.nanostream_bin {
            app.nanostream.executable = bin;
        }
        if let Some(primers) = startup.primers_path {
            app.nanostream.primers_path = primers;
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
        
        apply_theme_visuals(&cc.egui_ctx, default_theme);
        app.theme_applied = true;
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

        // Poll LAN Remote channel
        if let Some(rx) = &self.remote_rx {
            while let Ok(event) = rx.try_recv() {
                match event {
                    MonitorEvent::Pong => {
                        self.remote.status = RemoteStatus::Connected;
                        self.log_lines.push("Remote LAN node connected!".into());
                    }
                    MonitorEvent::Progress { reads_processed, percent } => {
                        self.data.total_reads = reads_processed;
                        self.log_lines.push(format!("Remote progress: {} reads ({:.1}%)", reads_processed, percent));
                    }
                    MonitorEvent::ResultSummary { total_reads, filtered_reads } => {
                        self.data.total_reads = total_reads;
                        self.data.filtered_reads = filtered_reads;
                        self.log_lines.push(format!("Remote Summary: Total={}, Filtered={}", total_reads, filtered_reads));
                    }
                    MonitorEvent::Error { message } => {
                        if message.contains("Connection failed") || message.contains("closed") {
                            self.remote.status = RemoteStatus::Disconnected;
                        }
                        self.last_error = Some(message.clone());
                        self.log_lines.push(format!("Remote error: {}", message));
                    }
                }
            }
        }

        let mut completed = false;
        let mut completed_path: Option<String> = None;
        let mut completed_ok = false;
        
        let mut worker_msgs = Vec::new();
        if let Some(rx) = &self.worker_rx {
            while let Ok(msg) = rx.try_recv() {
                worker_msgs.push(msg);
            }
        }
        for msg in worker_msgs {
            match msg {
                WorkerMessage::Log(line) => {
                    for chunk in line.lines() {
                        self.log_lines.push(chunk.to_string());
                    }
                }
                WorkerMessage::Progress { input_path, processed_reads: _, current_data } => {
                    self.data = current_data;
                    self.data.accumulated_files = 1;
                    self.current_input = Some(input_path);
                }
                WorkerMessage::Completed { input_path, result } => {
                    completed = true;
                    completed_path = Some(input_path.clone());
                    self.run_state = RunState::Idle;
                    match result {
                        Ok(new_data) => {
                            completed_ok = true;
                            self.last_error = None;
                            self.log_lines.push(format!("Analysis complete: {}", input_path));
                            
                            // Dynamic stats aggregation for live watch mode
                            if self.run.source == RunSource::MonitorDirectory && self.data.accumulated_files > 0 {
                                self.log_lines.push("Merging statistics with existing run data".into());
                                self.data.merge(new_data);
                            } else {
                                self.data = new_data;
                                self.data.accumulated_files = 1;
                            }

                            if self.run.auto_scan_variants {
                                self.log_lines.push("Auto-variant calling sequence queued".into());
                                self.trigger_variant_calling();
                            }
                        }
                        Err(err) => {
                            self.last_error = Some(err.clone());
                            self.log_lines.push(format!("Analysis failed for {}: {}", input_path, err));
                        }
                    }
                }
                WorkerMessage::SnapCompleted { path } => {
                    self.snap_running = false;
                    self.data.snapshot_img_path = Some(path);
                    self.log_lines.push("rs-qc Alignment snapshot generation complete".into());
                    self.active_center_sub_tab = 0; // Switch tab to alignment snapshot view
                }
                WorkerMessage::VariantsCompleted { variants } => {
                    self.variant_calling_running = false;
                    self.data.variants = variants;
                    self.log_lines.push("rindels Variant calling pipeline complete".into());
                    
                    // Update result row variants count
                    if let Some(idx) = self.data.selected_row {
                        if let Some(row) = self.data.rows.get_mut(idx) {
                            row.variants = self.data.variants.len() as u32;
                        }
                    }
                    self.active_center_sub_tab = 1; // Switch tab to variants view
                }
                WorkerMessage::Error(err) => {
                    self.snap_running = false;
                    self.variant_calling_running = false;
                    self.last_error = Some(err.clone());
                    self.log_lines.push(format!("Background pipeline failed: {}", err));
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
        } else if self.run_state == RunState::Running || self.snap_running || self.variant_calling_running {
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
                RichText::new("PCR Amplicon Suite")
                    .color(Color32::from_rgb(120, 120, 125))
                    .small(),
            );
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let run_label = if self.monitor_active {
                    "Watcher: Monitoring Run"
                } else {
                    match self.run_state {
                        RunState::Idle => "Watcher: Idle",
                        RunState::Running => "Watcher: Running Analysis",
                    }
                };
                ui.colored_label(
                    if self.run_state == RunState::Running || self.monitor_active {
                        Color32::from_rgb(45, 212, 191) // Biotech Teal highlight
                    } else {
                        Color32::from_rgb(120, 120, 125)
                    },
                    run_label,
                );
                ui.separator();
                let (color, label) = match self.remote.status {
                    RemoteStatus::Connected => {
                        (Color32::from_rgb(45, 212, 191), "Remote Node: Connected")
                    }
                    RemoteStatus::Connecting => {
                        (Color32::from_rgb(234, 179, 8), "Remote Node: Connecting...")
                    }
                    RemoteStatus::Disconnected => {
                        (Color32::from_rgb(239, 68, 68), "Remote Node: Offline")
                    }
                };
                ui.colored_label(color, label);
            });
        });
        ui.separator();
    }

    fn draw_sidebar(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.spacing_mut().item_spacing = Vec2::new(8.0, 8.0);

        // Visual Theme Skin selector
        ui.group(|ui| {
            ui.label(RichText::new("Aesthetics & Theme Skin").strong());
            let current_skin_label = self.active_theme.label();
            
            let mut changed = false;
            egui::ComboBox::from_id_salt("theme_selector_combobox")
                .selected_text(current_skin_label)
                .show_ui(ui, |ui| {
                    for theme in ThemeSkin::ALL {
                        if ui.selectable_value(&mut self.active_theme, theme, theme.label()).clicked() {
                            changed = true;
                        }
                    }
                });

            if changed {
                apply_theme_visuals(ctx, self.active_theme);
                self.log_lines.push(format!("Interface theme set to: {}", self.active_theme.label()));
            }
        });

        // Loaded inputs & annotation resources
        ui.group(|ui| {
            ui.label(RichText::new("Annotation Resources").strong());
            
            ui.horizontal(|ui| {
                ui.label("Primers TSV:");
                ui.text_edit_singleline(&mut self.nanostream.primers_path);
                if ui.button("...").clicked() {
                    self.pick_primers_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("Ref FASTA:");
                ui.text_edit_singleline(&mut self.reference_path);
                if ui.button("...").clicked() {
                    self.pick_reference_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("GTF / BED:");
                ui.text_edit_singleline(&mut self.gtf_path);
                if ui.button("...").clicked() {
                    self.pick_gtf_file();
                }
            });
            ui.label(format!("Clinical Mutations Loaded: {} entries", self.clinical_db.len()));
        });

        // Pipeline Tool Executables Paths configuration
        ui.collapsing("Pipeline Tool Paths", |ui| {
            ui.label("Local directories override PATH lookup:");
            ui.horizontal(|ui| {
                ui.label("nanostream:");
                ui.text_edit_singleline(&mut self.nanostream.executable);
            });
            ui.horizontal(|ui| {
                ui.label("rs-qc:");
                ui.text_edit_singleline(&mut self.tool_rs_qc);
            });
            ui.horizontal(|ui| {
                ui.label("rindels:");
                ui.text_edit_singleline(&mut self.tool_rindels);
            });
        });

        // Sequence run watcher controls
        ui.group(|ui| {
            ui.label(RichText::new("Amplicon Run Control").strong());
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.run.source, RunSource::SingleFile, "Single File");
                ui.selectable_value(
                    &mut self.run.source,
                    RunSource::MonitorDirectory,
                    "Monitor Dir",
                );
            });
            
            match self.run.source {
                RunSource::SingleFile => {
                    ui.label("Sequencing Input file (.bam / .fastq)");
                    ui.horizontal(|ui| {
                        let prev_path = self.run.input_path.clone();
                        if ui.text_edit_singleline(&mut self.run.input_path).changed() {
                            if self.run.input_path != prev_path {
                                self.data = DashboardData::empty();
                                self.processed_files.clear();
                                self.failed_files.clear();
                                self.pending_files.clear();
                                self.queued_files.clear();
                                self.current_input = None;
                                self.last_error = None;
                            }
                        }
                        if ui.button("Browse").clicked() {
                            self.pick_input_file();
                        }
                    });
                }
                RunSource::MonitorDirectory => {
                    ui.label("Live monitoring directory");
                    ui.horizontal(|ui| {
                        ui.text_edit_singleline(&mut self.run.monitor_dir);
                        if ui.button("Browse").clicked() {
                            self.pick_monitor_dir();
                        }
                    });
                    ui.label(format!("Accumulated Files: {}", self.data.accumulated_files));
                }
            }

            ui.checkbox(&mut self.run.auto_scan_variants, "Auto-variant calling hook");
            ui.label(format!(
                "Queue: {} | Processed: {} | Failed: {}",
                self.pending_files.len(),
                self.processed_files.len(),
                self.failed_files.len()
            ));
            if let Some(err) = &self.last_error {
                ui.colored_label(
                    Color32::from_rgb(239, 68, 68),
                    format!("Last failure: {}", err),
                );
            }
            if let Some(current) = &self.current_input {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Processing: {}", current));
                });
                if self.data.total_reads > 0 {
                    ui.label(format!("Analyzed so far: {} reads", self.data.total_reads));
                }
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

        // Remote LAN streaming configuration
        ui.group(|ui| {
            ui.label(RichText::new("LAN Remote Streaming").strong());
            ui.checkbox(&mut self.remote.enabled, "Enable LAN client");
            ui.label("Node Endpoint");
            ui.text_edit_singleline(&mut self.remote.endpoint);
            ui.label("Secret Handshake Token");
            ui.add(egui::TextEdit::singleline(&mut self.remote.auth_token).password(true));
            
            ui.horizontal(|ui| {
                if ui.button("Connect").clicked() {
                    self.remote.status = RemoteStatus::Connecting;
                    self.log_lines.push(format!("Connecting to LAN node: {}", self.remote.endpoint));
                    
                    let (tx, rx) = mpsc::channel::<MonitorEvent>();
                    self.remote_rx = Some(rx);
                    let (req_tx, stop_tx) = remote::spawn_remote_client(
                        self.remote.endpoint.clone(),
                        self.remote.auth_token.clone(),
                        tx,
                    );
                    self.remote_req_tx = Some(req_tx);
                    self.remote_stop_tx = Some(stop_tx);
                }
                if ui.button("Disconnect").clicked() {
                    if let Some(stop_tx) = self.remote_stop_tx.take() {
                        let _ = stop_tx.send(());
                    }
                    self.remote.status = RemoteStatus::Disconnected;
                    self.remote_rx = None;
                    self.remote_req_tx = None;
                    self.log_lines.push("LAN connection closed".into());
                }
            });
            ui.label(format!("Status: {}", self.remote.status.label()));
        });
    }

    fn draw_filter_strip(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Interactive Quality Filters").strong());
                ui.add(egui::DragValue::new(&mut self.filters.min_qs).speed(0.1).prefix("Min QS "));
                ui.add(egui::DragValue::new(&mut self.filters.min_len).speed(10.0).prefix("Min Len "));
                ui.add(egui::DragValue::new(&mut self.filters.max_reads).speed(100.0).prefix("Max Reads "));
                ui.label(RichText::new("(0 = all reads)").small());
                ui.checkbox(&mut self.filters.duplex_only, "Duplex only");
                ui.checkbox(&mut self.filters.use_nanostream, "Use Rust core (nanostream)");
                if ui.button("Recalculate").clicked() {
                    self.log_lines.push("Recalculate requested with updated quality thresholds".into());
                }
            });
        });
    }

    fn draw_result_table(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.tab, MainTab::Results, "Dashboard Table");
            ui.selectable_value(&mut self.tab, MainTab::Log, "Diagnostic logs");
        });
        ui.separator();

        match self.tab {
            MainTab::Results => {
                egui::ScrollArea::vertical()
                    .max_height(260.0)
                    .show(ui, |ui| {
                        let text_color = ui.visuals().widgets.noninteractive.fg_stroke.color;
                        let stroke = Stroke::new(1.0, ui.visuals().widgets.noninteractive.bg_stroke.color);
                        
                        egui::Frame::default().stroke(stroke).show(ui, |ui| {
                            ui.horizontal(|ui| {
                                ui.add_sized(
                                    [40.0, 18.0],
                                    egui::Label::new(RichText::new("#").strong()),
                                );
                                ui.separator();
                                ui.add_sized(
                                    [400.0, 18.0],
                                    egui::Label::new(RichText::new("PCR Target Amplicon Name").strong()),
                                );
                                ui.separator();
                                ui.add_sized(
                                    [100.0, 18.0],
                                    egui::Label::new(RichText::new("Read Count").strong()),
                                );
                                ui.separator();
                                ui.add_sized(
                                    [110.0, 18.0],
                                    egui::Label::new(RichText::new("Median Length").strong()),
                                );
                                ui.separator();
                                ui.add_sized(
                                    [100.0, 18.0],
                                    egui::Label::new(RichText::new("SD Length").strong()),
                                );
                                ui.separator();
                                ui.add_sized(
                                    [100.0, 18.0],
                                    egui::Label::new(RichText::new("Avg Q-Score").strong()),
                                );
                                ui.separator();
                                ui.add_sized(
                                    [100.0, 18.0],
                                    egui::Label::new(RichText::new("Called Variants").strong()),
                                );
                            });
                            ui.separator();

                            if self.data.rows.is_empty() {
                                ui.label(
                                    RichText::new(
                                        "No amplicon matches loaded. Start monitor/load inputs to process BAM/FASTQ sequencing files.",
                                    )
                                    .italics()
                                    .color(Color32::from_rgb(120, 120, 120)),
                                );
                            } else {
                                for (idx, row) in self.data.rows.iter().enumerate() {
                                    let selected = self.data.selected_row == Some(idx);
                                    
                                    ui.horizontal(|ui| {
                                        let num_res = ui.add_sized(
                                            [40.0, 18.0],
                                            egui::SelectableLabel::new(selected, format!("{}", idx + 1)),
                                        );
                                        if num_res.clicked() {
                                            self.data.selected_row = Some(idx);
                                            self.data.snapshot_img_path = None;
                                        }
                                        ui.separator();
                                        
                                        let name_res = ui.add_sized(
                                            [400.0, 18.0],
                                            egui::SelectableLabel::new(selected, &row.amplicon_name),
                                        );
                                        if name_res.clicked() {
                                            self.data.selected_row = Some(idx);
                                            self.data.snapshot_img_path = None;
                                        }
                                        ui.separator();
                                        
                                        ui.add_sized(
                                            [100.0, 18.0],
                                            egui::Label::new(format!("{}", row.count)),
                                        );
                                        ui.separator();
                                        ui.add_sized(
                                            [110.0, 18.0],
                                            egui::Label::new(format!("{} bp", row.median_length)),
                                        );
                                        ui.separator();
                                        ui.add_sized(
                                            [100.0, 18.0],
                                            egui::Label::new(format!("{:.1}", row.sd_length)),
                                        );
                                        ui.separator();
                                        ui.add_sized(
                                            [100.0, 18.0],
                                            egui::Label::new(format!("{:.1}", row.avg_qs)),
                                        );
                                        ui.separator();
                                        ui.add_sized(
                                            [100.0, 18.0],
                                            egui::Label::new(if row.variants > 0 {
                                                RichText::new(format!("{} mutations", row.variants))
                                                    .color(Color32::from_rgb(239, 68, 68))
                                                    .strong()
                                            } else {
                                                RichText::new("0").color(text_color)
                                            }),
                                        );
                                    });
                                    ui.separator();
                                }
                            }
                        });
                    });
            }
            MainTab::Log => {
                egui::ScrollArea::vertical()
                    .max_height(260.0)
                    .show(ui, |ui| {
                        for line in self.log_lines.iter().rev().take(120) {
                            ui.monospace(line);
                        }
                    });
            }
        }

        // Row details & manual hooks triggers
        let (sel_amplicons, sel_reads) = self.data.selected_counts();
        ui.separator();
        ui.horizontal(|ui| {
            ui.label(format!(
                "Active selection: {} amplicon(s), {} mapped reads",
                sel_amplicons, sel_reads
            ));
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let has_selection = self.data.selected_row.is_some();
                
                // Alignment Snap hook button
                let btn_snap = ui.add_enabled(
                    has_selection && !self.snap_running,
                    egui::Button::new("Generate Region Alignment Snapshot"),
                );
                if btn_snap.clicked() {
                    self.trigger_region_snapshot();
                }
                
                // Variants hook button
                let btn_vars = ui.add_enabled(
                    has_selection && !self.variant_calling_running,
                    egui::Button::new("Call Regional Mutations (rindels)"),
                );
                if btn_vars.clicked() {
                    self.trigger_variant_calling();
                }
            });
        });
    }

    fn draw_center_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.active_center_sub_tab, 0, "Alignment Snapshot (rs-qc)");
            ui.selectable_value(&mut self.active_center_sub_tab, 1, "Called Variants Table (rindels)");
        });
        ui.separator();

        match self.active_center_sub_tab {
            0 => {
                ui.group(|ui| {
                    ui.label(RichText::new("Genomic Alignment Snap browser representation").strong());
                    if self.snap_running {
                        ui.horizontal(|ui| {
                            ui.spinner();
                            ui.label("Rendering local alignment region PNG using rs-qc snap...");
                        });
                    } else if let Some(img_path) = &self.data.snapshot_img_path {
                        ui.label(format!("Alignment File location: {}", img_path));
                        
                        // Load image using local eframe file provider
                        egui::ScrollArea::both().show(ui, |ui| {
                            let img = egui::Image::from_uri(format!("file://{}", img_path))
                                .fit_to_original_size(1.0)
                                .max_width(ui.available_width());
                            ui.add(img);
                        });
                    } else {
                        ui.label(
                            RichText::new("No alignment snap loaded. Select an amplicon above and click 'Generate Region Alignment Snapshot'.").italics()
                        );
                    }
                });
            }
            1 => {
                ui.group(|ui| {
                    ui.label(RichText::new("rindels Local Assembly-based Variant Calls").strong());
                    if self.variant_calling_running {
                        ui.horizontal(|ui| {
                            ui.spinner();
                            ui.label("Calling local insertions, deletions and single nucleotide polymorphisms...");
                        });
                    } else if self.data.variants.is_empty() {
                        ui.label(
                            RichText::new("No variants called for this region. Select an amplicon above and press 'Call Regional Mutations'.").italics()
                        );
                    } else {
                        egui::ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                            let stroke = Stroke::new(1.0, ui.visuals().widgets.noninteractive.bg_stroke.color);
                            egui::Frame::default().stroke(stroke).show(ui, |ui| {
                                ui.horizontal(|ui| {
                                    ui.add_sized([100.0, 18.0], egui::Label::new(RichText::new("Contig").strong()));
                                    ui.add_sized([110.0, 18.0], egui::Label::new(RichText::new("Position").strong()));
                                    ui.add_sized([90.0, 18.0], egui::Label::new(RichText::new("Reference").strong()));
                                    ui.add_sized([90.0, 18.0], egui::Label::new(RichText::new("Mutation").strong()));
                                    ui.add_sized([90.0, 18.0], egui::Label::new(RichText::new("Depth").strong()));
                                    ui.add_sized([90.0, 18.0], egui::Label::new(RichText::new("VAF (%)").strong()));
                                    ui.add_sized([150.0, 18.0], egui::Label::new(RichText::new("ClinVar Annotation").strong()));
                                    ui.add_sized([220.0, 18.0], egui::Label::new(RichText::new("Classification Verdict").strong()));
                                    ui.allocate_space(ui.available_size());
                                });
                                ui.separator();

                                for var in &self.data.variants {
                                    ui.horizontal(|ui| {
                                        ui.add_sized([100.0, 18.0], egui::Label::new(&var.chrom));
                                        ui.add_sized([110.0, 18.0], egui::Label::new(format!("{}", var.position)));
                                        ui.add_sized([90.0, 18.0], egui::Label::new(&var.ref_allele));
                                        ui.add_sized([90.0, 18.0], egui::Label::new(&var.alt_allele));
                                        ui.add_sized([90.0, 18.0], egui::Label::new(format!("{}", var.depth)));
                                        ui.add_sized([90.0, 18.0], egui::Label::new(format!("{:.2}%", var.vaf)));
                                        
                                        // Highlights ClinVar matching
                                        let clinvar_color = if var.clinvar != "Unknown / VUS" {
                                            Color32::from_rgb(239, 68, 68) // Bright red clinical annotation
                                        } else {
                                            ui.visuals().widgets.noninteractive.fg_stroke.color
                                        };
                                        ui.add_sized([150.0, 18.0], egui::Label::new(
                                            RichText::new(&var.clinvar).color(clinvar_color).strong()
                                        ));

                                        // Highlights clinical hot-spots
                                        let verdict_color = if var.verdict.contains("★ HOTSPOT") {
                                            Color32::from_rgb(234, 179, 8) // Golden yellow highlights for clinical mutations!
                                        } else {
                                            ui.visuals().widgets.noninteractive.fg_stroke.color
                                        };
                                        ui.add_sized([220.0, 18.0], egui::Label::new(
                                            RichText::new(&var.verdict).color(verdict_color).strong()
                                        ));
                                        ui.allocate_space(ui.available_size());
                                    });
                                }
                            });
                        });
                    }
                });
            }
            _ => {}
        }
    }

    fn draw_bottom_plots(&mut self, ui: &mut egui::Ui) {
        let (acc_color, qs_color, len_color) = theme_colors(self.active_theme);

        let mut accuracy_hover = None;
        let mut qs_hover = None;
        let mut len_hover = None;

        let mut accuracy_hovered = false;
        let mut qs_hovered = false;
        let mut len_hovered = false;

        let mut clicked_accuracy = None;
        let mut clicked_qs = None;
        let mut clicked_len = None;

        let prev_acc_hovered = self.accuracy_plot_hovered;
        let prev_qs_hovered = self.qs_plot_hovered;
        let prev_len_hovered = self.len_plot_hovered;

        ui.columns(3, |columns| {
            columns[0].group(|ui| {
                ui.label(RichText::new("Sequence Alignment Accuracy Density").strong());
                if self.data.accuracy_bins.is_empty() {
                    ui.label(RichText::new("No sequence run data loaded.").italics());
                } else {
                    ui.label(
                        RichText::new(format!("Accuracy Mode: {:.2}%", self.data.accuracy_mode * 100.0))
                            .color(acc_color)
                            .strong(),
                    );
                }
                
                let points = density_points(&self.data.accuracy_bins);
                let line = Line::new(PlotPoints::from_iter(points));
                let target_acc = Self::phred_to_accuracy_pct(self.filters.min_qs as f64);

                let plot_response = Plot::new("accuracy_distribution_density")
                    .allow_drag(false)
                    .allow_zoom(false)
                    .allow_scroll(false)
                    .height(200.0)
                    .show(ui, |plot_ui| {
                        if !self.data.accuracy_bins.is_empty() {
                            plot_ui.line(line.color(acc_color).width(2.0));
                            
                            // Draw red vertical line for current threshold
                            plot_ui.vline(
                                egui_plot::VLine::new(target_acc)
                                    .color(Color32::from_rgb(239, 68, 68))
                                    .style(egui_plot::LineStyle::Dashed { length: 4.0 })
                            );
                        }

                        if prev_acc_hovered {
                            if let Some(coord) = plot_ui.pointer_coordinate() {
                                accuracy_hover = Some(coord);
                                // Draw hover crossbar
                                plot_ui.vline(
                                    egui_plot::VLine::new(coord.x)
                                        .color(Color32::GRAY.linear_multiply(0.4))
                                );
                            }
                        }
                    });

                accuracy_hovered = plot_response.response.hovered();

                if plot_response.response.clicked() {
                    if let Some(coord) = accuracy_hover {
                        clicked_accuracy = Some(coord.x);
                    }
                }
            });

            columns[1].group(|ui| {
                ui.label(RichText::new("Base Quality Q-Score Density").strong());
                if self.data.qs_bins.is_empty() {
                    ui.label(RichText::new("No sequence run data loaded.").italics());
                } else {
                    ui.label(
                        RichText::new(format!("Phred Score Mode: Q{:.1}", self.data.qs_mode))
                            .color(qs_color)
                            .strong(),
                    );
                }
                
                let points = density_points(&self.data.qs_bins);
                let line = Line::new(PlotPoints::from_iter(points));
                let target_qs = self.filters.min_qs as f64;

                let plot_response = Plot::new("qs_distribution_density")
                    .allow_drag(false)
                    .allow_zoom(false)
                    .allow_scroll(false)
                    .height(200.0)
                    .show(ui, |plot_ui| {
                        if !self.data.qs_bins.is_empty() {
                            plot_ui.line(line.color(qs_color).width(2.0));
                            
                            // Draw red vertical line for current threshold
                            plot_ui.vline(
                                egui_plot::VLine::new(target_qs)
                                    .color(Color32::from_rgb(239, 68, 68))
                                    .style(egui_plot::LineStyle::Dashed { length: 4.0 })
                            );
                        }

                        if prev_qs_hovered {
                            if let Some(coord) = plot_ui.pointer_coordinate() {
                                qs_hover = Some(coord);
                                // Draw hover crossbar
                                plot_ui.vline(
                                    egui_plot::VLine::new(coord.x)
                                        .color(Color32::GRAY.linear_multiply(0.4))
                                );
                            }
                        }
                    });

                qs_hovered = plot_response.response.hovered();

                if plot_response.response.clicked() {
                    if let Some(coord) = qs_hover {
                        clicked_qs = Some(coord.x);
                    }
                }
            });

            columns[2].group(|ui| {
                ui.label(RichText::new("Fragment Length Distribution Histogram").strong());
                if self.data.length_bins.is_empty() {
                    ui.label(RichText::new("No sequence run data loaded.").italics());
                } else {
                    ui.label(
                        RichText::new(format!("Median read length: {:.0} bp", self.data.length_median))
                            .color(len_color)
                            .strong(),
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
                let target_len = self.filters.min_len as f64;

                let plot_response = Plot::new("fragment_length_hist")
                    .allow_drag(false)
                    .allow_zoom(false)
                    .allow_scroll(false)
                    .height(200.0)
                    .show(ui, |plot_ui| {
                        if !self.data.length_bins.is_empty() {
                            plot_ui.bar_chart(
                                BarChart::new(bars).color(len_color),
                            );
                            
                            // Draw red vertical line for current threshold
                            plot_ui.vline(
                                egui_plot::VLine::new(target_len)
                                    .color(Color32::from_rgb(239, 68, 68))
                                    .style(egui_plot::LineStyle::Dashed { length: 4.0 })
                            );
                        }

                        if prev_len_hovered {
                            if let Some(coord) = plot_ui.pointer_coordinate() {
                                len_hover = Some(coord);
                                // Draw hover crossbar
                                plot_ui.vline(
                                    egui_plot::VLine::new(coord.x)
                                        .color(Color32::GRAY.linear_multiply(0.4))
                                );
                            }
                        }
                    });

                len_hovered = plot_response.response.hovered();

                if plot_response.response.clicked() {
                    if let Some(coord) = len_hover {
                        clicked_len = Some(coord.x);
                    }
                }
            });
        });

        // Save hovered states for the next frame
        self.accuracy_plot_hovered = accuracy_hovered;
        self.qs_plot_hovered = qs_hovered;
        self.len_plot_hovered = len_hovered;

        // Tooltips on hover
        if accuracy_hovered {
            if let Some(coord) = accuracy_hover {
                let q_equiv = Self::accuracy_pct_to_phred(coord.x);
                egui::show_tooltip_at_pointer(ui.ctx(), ui.layer_id(), egui::Id::new("acc_plot_tooltip"), |ui| {
                    ui.label(format!("Accuracy: {:.2}% (equiv. to Q{:.1})", coord.x, q_equiv));
                    ui.label("Click to set Q-Score filter threshold.");
                });
            }
        }
        if qs_hovered {
            if let Some(coord) = qs_hover {
                let acc_equiv = Self::phred_to_accuracy_pct(coord.x);
                egui::show_tooltip_at_pointer(ui.ctx(), ui.layer_id(), egui::Id::new("qs_plot_tooltip"), |ui| {
                    ui.label(format!("Phred Q-Score: Q{:.1} (equiv. to {:.2}% accuracy)", coord.x, acc_equiv));
                    ui.label("Click to set Q-Score filter threshold.");
                });
            }
        }
        if len_hovered {
            if let Some(coord) = len_hover {
                egui::show_tooltip_at_pointer(ui.ctx(), ui.layer_id(), egui::Id::new("len_plot_tooltip"), |ui| {
                    ui.label(format!("Read Length: {:.0} bp", coord.x));
                    ui.label("Click to set minimum length filter.");
                });
            }
        }

        // Apply setting updates from clicks
        if let Some(x) = clicked_accuracy {
            let q = Self::accuracy_pct_to_phred(x);
            self.filters.min_qs = q.clamp(0.0, 40.0) as f32;
        }
        if let Some(x) = clicked_qs {
            self.filters.min_qs = x.clamp(0.0, 40.0) as f32;
        }
        if let Some(x) = clicked_len {
            self.filters.min_len = x.clamp(0.0, 100000.0) as u32;
        }
    }

fn phred_to_accuracy_pct(qs: f64) -> f64 {
    let p_err = 10f64.powf(-qs / 10.0);
    ((1.0 - p_err) * 100.0).clamp(0.0, 100.0)
}

fn accuracy_pct_to_phred(acc: f64) -> f64 {
    let p_err = 1.0 - acc / 100.0;
    let p_err = p_err.clamp(1e-10, 1.0);
    -10.0 * p_err.log10()
}


    fn draw_file_ops_panel(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            ui.label(RichText::new("Post-Run Amplicon Extraction & barcode Demultiplexing").strong());
            
            ui.label("Export Mapped & Filtered Output Reads:");
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut self.file_op_output_path);
                if ui.button("Browse").clicked() {
                    self.pick_filter_output_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("Channel Filter range:");
                ui.text_edit_singleline(&mut self.file_op_channel_range);
                ui.label(RichText::new("(e.g. 1-128)").small());
            });
            ui.horizontal(|ui| {
                ui.label("Maximum Read Length limit:");
                ui.add(egui::DragValue::new(&mut self.file_op_max_len).speed(10.0));
                ui.label("bp (0 = none)");
            });
            
            if ui
                .add_enabled(!self.file_op_running, egui::Button::new("Filter & Extract Reads"))
                .clicked()
            {
                self.start_filter_export();
            }

            ui.separator();
            ui.label("Demultiplex by Barcodes List:");
            ui.horizontal(|ui| {
                ui.label("Barcodes file:");
                ui.text_edit_singleline(&mut self.barcode_file_path);
                if ui.button("...").clicked() {
                    self.pick_barcodes_file();
                }
            });
            ui.horizontal(|ui| {
                ui.label("Output Directory:");
                ui.text_edit_singleline(&mut self.barcode_output_dir);
                if ui.button("...").clicked() {
                    self.pick_barcode_output_dir();
                }
            });
            
            if ui
                .add_enabled(!self.file_op_running, egui::Button::new("Split Mapped Barcodes"))
                .clicked()
            {
                self.start_barcode_split();
            }
        });
    }
}

impl eframe::App for NanoMonitorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_worker(ctx);
        if self.monitor_active {
            ctx.request_repaint_after(Duration::from_millis(300));
        }

        // Renders visual theme skin colors on context repaint
        if !self.theme_applied {
            apply_theme_visuals(ctx, self.active_theme);
            self.theme_applied = true;
        }

        egui::TopBottomPanel::top("header_strip_panel")
            .resizable(false)
            .show(ctx, |ui| self.draw_top_bar(ui));

        egui::SidePanel::left("input_controls_sidebar")
            .resizable(true)
            .default_width(360.0)
            .show(ctx, |ui| self.draw_sidebar(ui, ctx));

        egui::TopBottomPanel::bottom("fixed_bottom_plots")
            .resizable(true)
            .default_height(280.0)
            .show(ctx, |ui| {
                ui.add_space(4.0);
                self.draw_bottom_plots(ui);
                ui.add_space(4.0);
                
                ui.separator();
                ui.horizontal(|ui| {
                    ui.label(format!(
                        "Total Reads Loaded: {} | Quality Filter Mapped Reads: {}",
                        format_count(self.data.total_reads),
                        format_count(self.data.filtered_reads)
                    ));
                    
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui.button("Extract Demux Options").clicked() {
                            self.active_center_sub_tab = 0;
                            self.tab = MainTab::Log;
                        }
                    });
                });
                ui.add_space(4.0);
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            self.draw_filter_strip(ui);
            ui.add_space(8.0);
            
            self.draw_result_table(ui);
            ui.add_space(8.0);
            
            self.draw_center_panel(ui);
        });
    }
}

impl Drop for NanoMonitorApp {
    fn drop(&mut self) {
        self.stop_directory_watcher();
        if let Some(stop_tx) = self.remote_stop_tx.take() {
            let _ = stop_tx.send(());
        }
    }
}

fn build_dashboard_from_nanostream(output: &matcher::AmpliconResult) -> DashboardData {
    let mut rows = output
        .amplicons
        .iter()
        .map(|(name, s)| ResultRow {
            amplicon_name: name.clone(),
            count: s.count as u32,
            median_length: s.median_length as u32,
            sd_length: s.std_length as f32,
            avg_qs: s.avg_qs,
            variants: 0,
        })
        .collect::<Vec<_>>();
    rows.sort_by(|a, b| b.count.cmp(&a.count));

    let d = &output.distributions;
    let (length_bins, qs_bins, accuracy_bins, length_median, qs_mode, accuracy_mode) = (
        map_bins(&d.length_bins),
        map_bins(&d.qs_bins),
        map_bins(&d.accuracy_bins),
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
        length_bins,
        qs_bins,
        accuracy_bins,
        length_median,
        qs_mode,
        accuracy_mode,
        accumulated_files: 1,
        variants: Vec::new(),
        snapshot_img_path: None,
    }
}

fn map_bins(bins: &[matcher::DistributionBin]) -> Vec<HistogramBin> {
    bins.iter()
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

/// Applies visuals colors and layouts dynamically for customizable UI theme skins
fn apply_theme_visuals(ctx: &egui::Context, theme: ThemeSkin) {
    let mut visuals = match theme {
        ThemeSkin::Solarized | ThemeSkin::ClassicLight => egui::Visuals::light(),
        _ => egui::Visuals::dark(),
    };

    match theme {
        ThemeSkin::BioTeal => {
            // Obsidian Charcoal + Biotech Teal
            visuals.widgets.noninteractive.bg_fill = Color32::from_rgb(18, 18, 20);
            visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(33, 33, 38));
            visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(205, 205, 210));
            
            visuals.widgets.inactive.bg_fill = Color32::from_rgb(26, 26, 32);
            visuals.widgets.hovered.bg_fill = Color32::from_rgb(38, 38, 46);
            visuals.widgets.hovered.fg_stroke = Stroke::new(1.2, Color32::from_rgb(45, 212, 191));
            
            visuals.widgets.active.bg_fill = Color32::from_rgb(13, 148, 136);
            visuals.widgets.active.fg_stroke = Stroke::new(1.0, Color32::WHITE);
            visuals.selection.bg_fill = Color32::from_rgb(13, 148, 136);
        }
        ThemeSkin::Matrix => {
            // Pure Console Black + Terminal Lime Green
            visuals.widgets.noninteractive.bg_fill = Color32::BLACK;
            visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.5, Color32::from_rgb(0, 180, 0));
            visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(0, 255, 0));
            
            visuals.widgets.inactive.bg_fill = Color32::from_rgb(12, 12, 12);
            visuals.widgets.hovered.bg_fill = Color32::from_rgb(25, 25, 25);
            visuals.widgets.hovered.fg_stroke = Stroke::new(1.5, Color32::from_rgb(50, 255, 50));
            
            visuals.widgets.active.bg_fill = Color32::from_rgb(0, 150, 0);
            visuals.widgets.active.fg_stroke = Stroke::new(1.0, Color32::BLACK);
            visuals.selection.bg_fill = Color32::from_rgb(0, 100, 0);
        }
        ThemeSkin::Cyberpunk => {
            // Deep purple space + glowing pink neon
            visuals.widgets.noninteractive.bg_fill = Color32::from_rgb(16, 12, 30);
            visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(244, 63, 94));
            visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(254, 240, 138));
            
            visuals.widgets.inactive.bg_fill = Color32::from_rgb(24, 18, 44);
            visuals.widgets.hovered.bg_fill = Color32::from_rgb(38, 28, 68);
            visuals.widgets.hovered.fg_stroke = Stroke::new(1.2, Color32::from_rgb(236, 72, 153));
            
            visuals.widgets.active.bg_fill = Color32::from_rgb(244, 63, 94);
            visuals.widgets.active.fg_stroke = Stroke::new(1.0, Color32::WHITE);
            visuals.selection.bg_fill = Color32::from_rgb(244, 63, 94);
        }
        ThemeSkin::Solarized => {
            // Solarized Light Cream + amber details
            visuals.widgets.noninteractive.bg_fill = Color32::from_rgb(253, 246, 227);
            visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(238, 232, 213));
            visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(88, 110, 117));
            
            visuals.widgets.inactive.bg_fill = Color32::from_rgb(238, 232, 213);
            visuals.widgets.hovered.bg_fill = Color32::from_rgb(220, 214, 195);
            visuals.widgets.hovered.fg_stroke = Stroke::new(1.0, Color32::from_rgb(38, 139, 210));
            
            visuals.widgets.active.bg_fill = Color32::from_rgb(181, 137, 0);
            visuals.widgets.active.fg_stroke = Stroke::new(1.0, Color32::WHITE);
            visuals.selection.bg_fill = Color32::from_rgb(181, 137, 0);
        }
        ThemeSkin::ClassicLight => {
            // Clean Gray + Royal Slate Blue
            visuals.widgets.noninteractive.bg_fill = Color32::from_rgb(245, 245, 248);
            visuals.widgets.noninteractive.bg_stroke = Stroke::new(1.0, Color32::from_rgb(218, 218, 224));
            visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, Color32::from_rgb(45, 45, 48));
            
            visuals.widgets.inactive.bg_fill = Color32::from_rgb(232, 232, 238);
            visuals.widgets.hovered.bg_fill = Color32::from_rgb(220, 220, 228);
            visuals.widgets.hovered.fg_stroke = Stroke::new(1.0, Color32::from_rgb(37, 99, 235));
            
            visuals.widgets.active.bg_fill = Color32::from_rgb(37, 99, 235);
            visuals.widgets.active.fg_stroke = Stroke::new(1.0, Color32::WHITE);
            visuals.selection.bg_fill = Color32::from_rgb(37, 99, 235);
        }
    }

    ctx.set_visuals(visuals);
}

/// Dynamic plot theme color resolver
fn theme_colors(theme: ThemeSkin) -> (Color32, Color32, Color32) {
    match theme {
        ThemeSkin::BioTeal => (
            Color32::from_rgb(45, 212, 191), // Teal
            Color32::from_rgb(59, 130, 246), // Blue
            Color32::from_rgb(139, 92, 246), // Purple
        ),
        ThemeSkin::Matrix => (
            Color32::from_rgb(0, 255, 0),     // Green
            Color32::from_rgb(0, 180, 0),     // Medium Green
            Color32::from_rgb(50, 255, 50),   // Light Green
        ),
        ThemeSkin::Cyberpunk => (
            Color32::from_rgb(244, 63, 94),   // Pink
            Color32::from_rgb(6, 182, 212),   // Cyan
            Color32::from_rgb(250, 204, 21),  // Yellow
        ),
        ThemeSkin::Solarized => (
            Color32::from_rgb(42, 161, 152),  // Cyan
            Color32::from_rgb(38, 139, 210),  // Blue
            Color32::from_rgb(203, 75, 22),   // Orange
        ),
        ThemeSkin::ClassicLight => (
            Color32::from_rgb(76, 175, 80),   // Green
            Color32::from_rgb(33, 150, 243),  // Blue
            Color32::from_rgb(156, 39, 176),  // Purple
        ),
    }
}
