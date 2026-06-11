use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisMode {
    Amplicon,
}

impl AnalysisMode {
    pub const ALL: [AnalysisMode; 1] = [
        AnalysisMode::Amplicon,
    ];

    pub fn label(self) -> &'static str {
        match self {
            AnalysisMode::Amplicon => "PCR Amplicon",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThemeSkin {
    BioTeal,
    Matrix,
    Cyberpunk,
    Solarized,
    ClassicLight,
}

impl ThemeSkin {
    pub const ALL: [ThemeSkin; 5] = [
        ThemeSkin::BioTeal,
        ThemeSkin::Matrix,
        ThemeSkin::Cyberpunk,
        ThemeSkin::Solarized,
        ThemeSkin::ClassicLight,
    ];

    pub fn label(self) -> &'static str {
        match self {
            ThemeSkin::BioTeal => "Bio-Teal (Dark)",
            ThemeSkin::Matrix => "Matrix Terminal",
            ThemeSkin::Cyberpunk => "Cyberpunk Neon",
            ThemeSkin::Solarized => "Solarized Amber",
            ThemeSkin::ClassicLight => "Classic Light",
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
    pub use_nanostream: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_qs: 10.0,
            min_len: 300,
            max_reads: 0,
            duplex_only: false,
            use_nanostream: true,
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct HistogramBin {
    pub start: f64,
    pub end: f64,
    pub count: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantRow {
    pub chrom: String,
    pub position: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub depth: u32,
    pub vaf: f32,
    pub clinvar: String,
    pub verdict: String,
}

#[derive(Debug, Clone)]
pub struct DashboardData {
    pub rows: Vec<ResultRow>,
    pub selected_row: Option<usize>,
    pub total_reads: u64,
    pub filtered_reads: u64,
    pub length_bins: Vec<HistogramBin>,
    pub qs_bins: Vec<HistogramBin>,
    pub accuracy_bins: Vec<HistogramBin>,
    pub length_median: f64,
    pub qs_mode: f64,
    pub accuracy_mode: f64,
    pub accumulated_files: usize,
    pub variants: Vec<VariantRow>,
    pub snapshot_img_path: Option<String>,
}

impl DashboardData {
    pub fn empty() -> Self {
        Self {
            rows: Vec::new(),
            selected_row: None,
            total_reads: 0,
            filtered_reads: 0,
            length_bins: Vec::new(),
            qs_bins: Vec::new(),
            accuracy_bins: Vec::new(),
            length_median: 0.0,
            qs_mode: 0.0,
            accuracy_mode: 0.0,
            accumulated_files: 0,
            variants: Vec::new(),
            snapshot_img_path: None,
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

    pub fn merge(&mut self, other: DashboardData) {
        self.accumulated_files += 1;
        self.total_reads += other.total_reads;
        self.filtered_reads += other.filtered_reads;

        // Merge amplicon rows
        for other_row in other.rows {
            if let Some(existing) = self.rows.iter_mut().find(|r| r.amplicon_name == other_row.amplicon_name) {
                let total_cnt = existing.count + other_row.count;
                if total_cnt > 0 {
                    existing.avg_qs = (existing.avg_qs * (existing.count as f32) + other_row.avg_qs * (other_row.count as f32)) / (total_cnt as f32);
                    existing.median_length = ((existing.median_length as f32 * existing.count as f32 + other_row.median_length as f32 * other_row.count as f32) / total_cnt as f32) as u32;
                    existing.sd_length = (existing.sd_length * (existing.count as f32) + other_row.sd_length * (other_row.count as f32)) / (total_cnt as f32);
                }
                existing.count = total_cnt;
                existing.variants += other_row.variants;
            } else {
                self.rows.push(other_row);
            }
        }
        self.rows.sort_by(|a, b| b.count.cmp(&a.count));

        // Merge histograms
        fn merge_hists(existing: &mut Vec<HistogramBin>, other_bins: Vec<HistogramBin>) {
            if existing.is_empty() {
                *existing = other_bins;
                return;
            }
            for ob in other_bins {
                if let Some(eb) = existing.iter_mut().find(|b| (b.start - ob.start).abs() < 1e-5 && (b.end - ob.end).abs() < 1e-5) {
                    eb.count += ob.count;
                } else {
                    existing.push(ob);
                }
            }
            existing.sort_by(|a, b| a.start.partial_cmp(&b.start).unwrap_or(std::cmp::Ordering::Equal));
        }

        merge_hists(&mut self.length_bins, other.length_bins);
        merge_hists(&mut self.qs_bins, other.qs_bins);
        merge_hists(&mut self.accuracy_bins, other.accuracy_bins);

        // Recalculate summary metrics: length median, Q-score mode, accuracy mode
        if !self.length_bins.is_empty() {
            let total_len_count: u64 = self.length_bins.iter().map(|b| b.count).sum();
            let half = total_len_count / 2;
            let mut acc = 0;
            for b in &self.length_bins {
                acc += b.count;
                if acc >= half {
                    self.length_median = (b.start + b.end) * 0.5;
                    break;
                }
            }
        }

        fn find_mode(bins: &[HistogramBin]) -> f64 {
            let mut max_cnt = 0;
            let mut mode_val = 0.0;
            for b in bins {
                if b.count > max_cnt {
                    max_cnt = b.count;
                    mode_val = (b.start + b.end) * 0.5;
                }
            }
            mode_val
        }

        self.qs_mode = find_mode(&self.qs_bins);
        self.accuracy_mode = find_mode(&self.accuracy_bins);
    }
}
