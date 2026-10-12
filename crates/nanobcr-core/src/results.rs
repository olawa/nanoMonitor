//! nagini's final BCR-ABL1 results, and how they differ from the preliminary ones.

use crate::handoff::Handoff;
use crate::prelim::{
    parse_branford, parse_fusion, parse_phasing, read_tsv, BranfordSummary, FusionSummary, PhasePair, PrelimSummary,
};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct ExonRow {
    pub exon: String,
    pub mean_depth: f64,
    pub pct_gte_100x: f64,
    pub spanning_reads: u64,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ResultFiles {
    pub fusion_plot_png: Option<PathBuf>,
    pub fusion_plot_pdf: Option<PathBuf>,
    pub sample_workbook: Option<PathBuf>,
    pub reports_dir: Option<PathBuf>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct FinalSample {
    pub barcode: String,
    pub sample: String,
    pub sample_type: String,
    pub fusion: Option<FusionSummary>,
    pub branford: Option<BranfordSummary>,
    pub phasing: Option<Vec<PhasePair>>,
    pub exons: Vec<ExonRow>,
    pub files: ResultFiles,
    /// What could not be read, so a missing piece is never silently an empty one.
    pub problems: Vec<String>,
}

/// `results/reports/bcr_abl` (delivered), else `reports/bcr_abl` (working copy).
pub fn reports_dir(analysis_dir: &Path) -> Option<PathBuf> {
    ["results/reports/bcr_abl", "reports/bcr_abl"]
        .iter()
        .map(|p| analysis_dir.join(p))
        .find(|p| p.is_dir())
}

fn existing(p: PathBuf) -> Option<PathBuf> {
    p.is_file().then_some(p)
}

pub fn load_sample(analysis_dir: &Path, barcode: &str, sample: &str, sample_type: &str, min_depth: u64) -> FinalSample {
    let mut r = FinalSample {
        barcode: barcode.into(),
        sample: sample.into(),
        sample_type: sample_type.into(),
        ..Default::default()
    };
    let Some(dir) = reports_dir(analysis_dir) else {
        r.problems.push("no results/reports/bcr_abl folder yet".into());
        return r;
    };
    let key = format!("{sample}_{sample_type}");
    let f = |ext: &str| dir.join(format!("{key}.{ext}"));
    r.files.reports_dir = Some(dir.clone());
    r.files.fusion_plot_png = existing(f("fusion_plot.png"));
    r.files.fusion_plot_pdf = existing(f("fusion_plot.pdf"));
    r.files.sample_workbook = [
        dir.parent().map(|p| p.join(format!("{key}.nagini_report.xlsx"))),
        Some(analysis_dir.join(format!("reports/{key}.nagini_report.xlsx"))),
    ]
    .into_iter()
    .flatten()
    .find(|p| p.is_file());

    match parse_fusion(&f("fusion_summary.json")) {
        Ok(x) => r.fusion = Some(x),
        Err(e) => r.problems.push(format!("fusion_summary.json: {e:#}")),
    }
    match parse_branford(&f("branford_all_sites.tsv"), &f("branford_mutations.tsv"), min_depth) {
        Ok(x) => r.branford = Some(x),
        Err(e) => r.problems.push(format!("branford: {e:#}")),
    }
    match parse_phasing(&f("phasing.json")) {
        Ok(x) => r.phasing = Some(x),
        Err(e) => r.problems.push(format!("phasing.json: {e:#}")),
    }
    match read_tsv(&f("exon_coverage.tsv")) {
        Ok(rows) => {
            r.exons = rows
                .iter()
                .map(|x| ExonRow {
                    exon: x.get("exon").cloned().unwrap_or_default(),
                    mean_depth: x.get("mean_depth").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                    pct_gte_100x: x.get("pct_gte_100x").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                    spanning_reads: x.get("spanning_reads").and_then(|v| v.parse().ok()).unwrap_or(0),
                })
                .collect();
        }
        Err(e) => r.problems.push(format!("exon_coverage.tsv: {e:#}")),
    }
    r
}

/// Load every sample of a hand-off.
pub fn load_all(analysis_dir: &Path, handoff: &Handoff, min_depth: u64) -> Vec<FinalSample> {
    handoff
        .barcodes
        .iter()
        .map(|b| load_sample(analysis_dir, &b.barcode, &b.sample, &b.sample_type, min_depth))
        .collect()
}

pub fn summary_workbook(analysis_dir: &Path) -> Option<PathBuf> {
    crate::pipeline::summary_workbook(analysis_dir)
}

pub fn multiqc_report(analysis_dir: &Path) -> Option<PathBuf> {
    existing(analysis_dir.join("results/qc/multiqc/multiqc_nanopore.html"))
        .or_else(|| existing(analysis_dir.join("qc/multiqc/multiqc_nanopore.html")))
}

// ---------------------------------------------------------------- comparison

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Difference {
    Isoform { preliminary: Option<String>, final_call: Option<String> },
    MutationOnlyInPreliminary { mutation: String, vaf_percent: f64 },
    MutationOnlyInFinal { mutation: String, vaf_percent: f64 },
    VafShift { mutation: String, preliminary: f64, final_vaf: f64 },
    PhasingChanged { pair: String, preliminary: Option<String>, final_state: Option<String> },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Comparison {
    pub barcode: String,
    pub sample: String,
    /// Snapshot generation the preliminary result was computed on, against the one handed over.
    pub preliminary_generation: u64,
    pub handoff_generation: u64,
    pub differences: Vec<Difference>,
}

impl Comparison {
    pub fn identical(&self) -> bool {
        self.differences.is_empty()
    }
}

/// A VAF difference smaller than this (percentage points) is not worth reporting:
/// the scripts are the same and the reads should be too, so any shift means the
/// final run saw different data or settings.
pub const VAF_TOLERANCE: f64 = 0.05;

pub fn compare(prelim: &PrelimSummary, fin: &FinalSample, handoff_generation: u64) -> Comparison {
    let mut d = Vec::new();
    let p_iso = prelim.fusion.as_ref().and_then(|f| f.primary_isoform.clone());
    let f_iso = fin.fusion.as_ref().and_then(|f| f.primary_isoform.clone());
    if p_iso != f_iso {
        d.push(Difference::Isoform { preliminary: p_iso, final_call: f_iso });
    }

    let none = Vec::new();
    let pm = prelim.branford.as_ref().map(|b| &b.mutations).unwrap_or(&none);
    let fm = fin.branford.as_ref().map(|b| &b.mutations).unwrap_or(&none);
    for m in pm {
        match fm.iter().find(|x| x.mutation == m.mutation) {
            None => d.push(Difference::MutationOnlyInPreliminary { mutation: m.mutation.clone(), vaf_percent: m.vaf_percent }),
            Some(x) if (x.vaf_percent - m.vaf_percent).abs() > VAF_TOLERANCE => d.push(Difference::VafShift {
                mutation: m.mutation.clone(),
                preliminary: m.vaf_percent,
                final_vaf: x.vaf_percent,
            }),
            Some(_) => {}
        }
    }
    for x in fm {
        if !pm.iter().any(|m| m.mutation == x.mutation) {
            d.push(Difference::MutationOnlyInFinal { mutation: x.mutation.clone(), vaf_percent: x.vaf_percent });
        }
    }

    let pp = prelim.phasing.as_ref().unwrap_or(&Vec::new()).clone();
    let fp = fin.phasing.clone().unwrap_or_default();
    let key = |p: &PhasePair| format!("{}+{}", p.mutation_1, p.mutation_2);
    let mut keys: Vec<String> = pp.iter().chain(fp.iter()).map(key).collect();
    keys.sort();
    keys.dedup();
    for k in keys {
        let a = pp.iter().find(|p| key(p) == k).map(|p| p.phase_state.clone());
        let b = fp.iter().find(|p| key(p) == k).map(|p| p.phase_state.clone());
        if a != b {
            d.push(Difference::PhasingChanged { pair: k, preliminary: a, final_state: b });
        }
    }
    Comparison {
        barcode: fin.barcode.clone(),
        sample: fin.sample.clone(),
        preliminary_generation: prelim.generation,
        handoff_generation,
        differences: d,
    }
}

/// The preliminary summary the engine wrote for `barcode`, if any.
pub fn load_prelim(out_dir: &Path, barcode: &str) -> Result<Option<PrelimSummary>> {
    let p = out_dir.join(barcode).join("prelim/summary.json");
    if !p.is_file() {
        return Ok(None);
    }
    Ok(Some(serde_json::from_slice(&fs::read(p)?)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelim::{IsoformRow, MutationRow};

    fn mutation(name: &str, vaf: f64) -> MutationRow {
        MutationRow { mutation: name.into(), vaf_percent: vaf, ..Default::default() }
    }

    fn prelim(iso: &str, muts: Vec<MutationRow>, phase: &str) -> PrelimSummary {
        PrelimSummary {
            generation: 3,
            fusion: Some(FusionSummary {
                total_fusion_reads: 10,
                primary_isoform: Some(iso.into()),
                isoforms: vec![IsoformRow::default()],
            }),
            branford: Some(BranfordSummary { mutations: muts, ..Default::default() }),
            phasing: Some(vec![PhasePair { mutation_1: "Y253H".into(), mutation_2: "T315I".into(), phase_state: phase.into(), ..Default::default() }]),
            ..Default::default()
        }
    }

    fn fin(iso: &str, muts: Vec<MutationRow>, phase: &str) -> FinalSample {
        let p = prelim(iso, muts, phase);
        FinalSample { barcode: "barcode01".into(), sample: "A".into(), fusion: p.fusion, branford: p.branford, phasing: p.phasing, ..Default::default() }
    }

    #[test]
    fn identical_results_compare_equal_even_with_negligible_vaf_noise() {
        let p = prelim("e14a2", vec![mutation("T315I", 2.349)], "CIS (Compound)");
        let f = fin("e14a2", vec![mutation("T315I", 2.351)], "CIS (Compound)");
        assert!(compare(&p, &f, 3).identical());
    }

    #[test]
    fn every_kind_of_difference_is_reported() {
        let p = prelim("e13a2", vec![mutation("T315I", 2.0), mutation("Y253H", 1.0)], "CIS (Compound)");
        let f = fin("e14a2", vec![mutation("T315I", 3.0), mutation("E255K", 1.2)], "TRANS (Polyclonal)");
        let c = compare(&p, &f, 4);
        assert_eq!((c.preliminary_generation, c.handoff_generation), (3, 4));
        assert!(c.differences.contains(&Difference::Isoform { preliminary: Some("e13a2".into()), final_call: Some("e14a2".into()) }));
        assert!(c.differences.contains(&Difference::VafShift { mutation: "T315I".into(), preliminary: 2.0, final_vaf: 3.0 }));
        assert!(c.differences.contains(&Difference::MutationOnlyInPreliminary { mutation: "Y253H".into(), vaf_percent: 1.0 }));
        assert!(c.differences.contains(&Difference::MutationOnlyInFinal { mutation: "E255K".into(), vaf_percent: 1.2 }));
        assert!(c.differences.iter().any(|d| matches!(d, Difference::PhasingChanged { .. })));
    }

    /// Lay out what nagini writes and read it back.
    #[test]
    fn loads_a_nagini_reports_folder() {
        let d = tempfile::tempdir().unwrap();
        let rep = d.path().join("results/reports/bcr_abl");
        fs::create_dir_all(&rep).unwrap();
        fs::write(
            rep.join("PAT_T.fusion_summary.json"),
            r#"{"total_fusion_reads":30,"primary_isoform":"e14a2","isoforms":[{"fusion_isoform":"e14a2","protein_form":"p210","read_count":30,"fraction_percent":100.0,"status":"PRIMARY"}]}"#,
        ).unwrap();
        fs::write(rep.join("PAT_T.branford_all_sites.tsv"), "# c\nsample\tchrom\tpos\tdepth\nx\tchr9\t1\t900\n").unwrap();
        fs::write(rep.join("PAT_T.branford_mutations.tsv"), "sample\tchrom\tpos\tmutation\taa_change\tdepth\tad_alt\tvaf_percent\tstatus\nx\tchr9\t1\tT315I\tp.T315I\t900\t20\t2.2\tDETECTED\n").unwrap();
        fs::write(rep.join("PAT_T.phasing.json"), r#"{"phased_pairs":[]}"#).unwrap();
        fs::write(rep.join("PAT_T.exon_coverage.tsv"), "sample\texon\tmean_depth\tpct_gte_100x\tspanning_reads\nx\tABL1_e2\t812.5\t100.0\t77\n").unwrap();
        fs::write(rep.join("PAT_T.fusion_plot.png"), "png").unwrap();
        fs::write(d.path().join("results/reports/PAT_T.nagini_report.xlsx"), "x").unwrap();

        let r = load_sample(d.path(), "barcode01", "PAT", "T", 50);
        assert!(r.problems.is_empty(), "{:?}", r.problems);
        assert_eq!(r.fusion.unwrap().primary_isoform.as_deref(), Some("e14a2"));
        assert_eq!(r.branford.unwrap().mutations[0].mutation, "T315I");
        assert_eq!(r.exons[0], ExonRow { exon: "ABL1_e2".into(), mean_depth: 812.5, pct_gte_100x: 100.0, spanning_reads: 77 });
        assert!(r.files.fusion_plot_png.is_some() && r.files.sample_workbook.is_some());
        assert!(r.files.fusion_plot_pdf.is_none());
    }

    #[test]
    fn missing_outputs_are_problems_not_empty_results() {
        let d = tempfile::tempdir().unwrap();
        assert!(load_sample(d.path(), "b", "S", "T", 50).problems[0].contains("no results/reports"));
        fs::create_dir_all(d.path().join("reports/bcr_abl")).unwrap();
        let r = load_sample(d.path(), "b", "S", "T", 50);
        assert!(r.fusion.is_none() && r.problems.len() == 4, "{:?}", r.problems);
    }
}
