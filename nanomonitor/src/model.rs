#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisMode {
    Amplicon,
    RnaSeq,
    Wgs,
}

impl AnalysisMode {
    pub const ALL: [AnalysisMode; 3] = [
        AnalysisMode::Amplicon,
        AnalysisMode::RnaSeq,
        AnalysisMode::Wgs,
    ];

    pub fn label(self) -> &'static str {
        match self {
            AnalysisMode::Amplicon => "Amplicon",
            AnalysisMode::RnaSeq => "RNA-Seq",
            AnalysisMode::Wgs => "WGS + CNV",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MainTab {
    Results,
    Log,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunSource {
    SingleFile,
    MonitorDirectory,
}

#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub min_qs: f32,
    pub min_len: u32,
    pub max_reads: u32,
    pub duplex_only: bool,
    pub use_nanoparse: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_qs: 10.0,
            min_len: 300,
            max_reads: 0,
            duplex_only: false,
            use_nanoparse: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub source: RunSource,
    pub input_path: String,
    pub monitor_dir: String,
    pub auto_scan_variants: bool,
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            source: RunSource::SingleFile,
            input_path: String::new(),
            monitor_dir: String::new(),
            auto_scan_variants: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ResultRow {
    pub amplicon_name: String,
    pub count: u32,
    pub median_length: u32,
    pub sd_length: f32,
    pub avg_qs: f32,
    pub variants: u32,
}

#[derive(Debug, Clone, Copy)]
pub struct CnvBin {
    pub position_mb: f64,
    pub log2_ratio: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct HistogramBin {
    pub start: f64,
    pub end: f64,
    pub count: u64,
}

#[derive(Debug, Clone)]
pub struct DashboardData {
    pub rows: Vec<ResultRow>,
    pub selected_row: Option<usize>,
    pub total_reads: u64,
    pub filtered_reads: u64,
    pub cnv_bins: Vec<CnvBin>,
    pub length_bins: Vec<HistogramBin>,
    pub qs_bins: Vec<HistogramBin>,
    pub accuracy_bins: Vec<HistogramBin>,
    pub length_median: f64,
    pub qs_mode: f64,
    pub accuracy_mode: f64,
}

impl DashboardData {
    pub fn empty() -> Self {
        Self {
            rows: Vec::new(),
            selected_row: None,
            total_reads: 0,
            filtered_reads: 0,
            cnv_bins: Vec::new(),
            length_bins: Vec::new(),
            qs_bins: Vec::new(),
            accuracy_bins: Vec::new(),
            length_median: 0.0,
            qs_mode: 0.0,
            accuracy_mode: 0.0,
        }
    }

    pub fn selected_counts(&self) -> (usize, u64) {
        if let Some(i) = self.selected_row {
            if let Some(row) = self.rows.get(i) {
                return (1, row.count as u64);
            }
        }
        (0, 0)
    }
}
