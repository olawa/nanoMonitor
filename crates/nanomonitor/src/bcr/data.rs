//! What the BCR-ABL view shows, read from the engine's output directory.
//!
//! The GUI never talks to the engine: it reads `state.json`, each barcode's
//! `prelim/summary.json` and `history.jsonl`. All of those are replaced
//! atomically by their writers, so a read sees a complete old or new file.

use nanobcr_core::prelim::PrelimSummary;
use nanobcr_core::state::State;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

/// The engine updates `state.json` at least every poll; an alignment cycle or a
/// compaction can take a while, so only a long silence means it is gone.
pub const STALE_AFTER_SECS: u64 = 120;

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct HistoryMutation {
    pub mutation: String,
    pub vaf_percent: f64,
    pub depth: u64,
    pub status: String,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct HistoryPoint {
    pub unix: u64,
    pub generation: u64,
    pub reads_passed: u64,
    pub status: String,
    pub primary_isoform: Option<String>,
    pub fusion_reads: Option<u64>,
    pub mutations: Option<Vec<HistoryMutation>>,
    pub enough: Option<bool>,
}

#[derive(Debug, Clone, Default)]
pub struct BcrData {
    pub state: Option<State>,
    pub prelim: BTreeMap<String, PrelimSummary>,
    pub history: BTreeMap<String, Vec<HistoryPoint>>,
    /// Files that exist but could not be read (surfaced, never silently ignored).
    pub problems: Vec<String>,
}

pub fn load(out_dir: &Path) -> BcrData {
    let mut data = BcrData::default();
    // A missing state.json just means the engine has not started.
    if let Ok(bytes) = fs::read(out_dir.join("state.json")) {
        match serde_json::from_slice::<State>(&bytes) {
            Ok(s) => data.state = Some(s),
            Err(e) => data.problems.push(format!("state.json unreadable: {e}")),
        }
    }
    let Ok(rd) = fs::read_dir(out_dir) else {
        return data;
    };
    for entry in rd.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        let prelim_dir = entry.path().join("prelim");
        if !prelim_dir.is_dir() {
            continue;
        }
        if let Ok(bytes) = fs::read(prelim_dir.join("summary.json")) {
            match serde_json::from_slice::<PrelimSummary>(&bytes) {
                Ok(s) => {
                    data.prelim.insert(name.clone(), s);
                }
                Err(e) => data.problems.push(format!("{name}/prelim/summary.json unreadable: {e}")),
            }
        }
        if let Ok(text) = fs::read_to_string(prelim_dir.join("history.jsonl")) {
            let points: Vec<HistoryPoint> = text
                .lines()
                .filter_map(|l| serde_json::from_str(l).ok())
                .collect();
            if !points.is_empty() {
                data.history.insert(name, points);
            }
        }
    }
    data
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Liveness {
    NotStarted,
    Running { heartbeat_age: u64 },
    /// State says running but nothing has been written for a long time.
    Stale { heartbeat_age: u64 },
    Finished,
    Error(String),
}

pub fn liveness(state: Option<&State>, now: u64) -> Liveness {
    let Some(s) = state else {
        return Liveness::NotStarted;
    };
    let age = now.saturating_sub(s.heartbeat_unix);
    match s.status.as_str() {
        "finished" => Liveness::Finished,
        _ if age > STALE_AFTER_SECS => Liveness::Stale { heartbeat_age: age },
        "error" => Liveness::Error(s.last_error.clone().unwrap_or_default()),
        _ => Liveness::Running { heartbeat_age: age },
    }
}

pub fn format_age(secs: u64) -> String {
    match secs {
        0..=59 => format!("{secs} s"),
        60..=3599 => format!("{} min", secs / 60),
        3600..=86_399 => format!("{:.1} h", secs as f64 / 3600.0),
        _ => format!("{:.1} d", secs as f64 / 86_400.0),
    }
}

pub fn format_count(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out.chars().rev().collect()
}

/// Mutation VAF over snapshots, one series per mutation: x = reads that passed the
/// filter at that snapshot, y = VAF %.
pub fn mutation_series(history: &[HistoryPoint]) -> BTreeMap<String, Vec<[f64; 2]>> {
    let mut out: BTreeMap<String, Vec<[f64; 2]>> = BTreeMap::new();
    for p in history {
        for m in p.mutations.iter().flatten() {
            out.entry(m.mutation.clone())
                .or_default()
                .push([p.reads_passed as f64, m.vaf_percent]);
        }
    }
    out
}

/// Fusion-spanning reads over snapshots (x = reads that passed the filter).
pub fn fusion_series(history: &[HistoryPoint]) -> Vec<[f64; 2]> {
    history
        .iter()
        .filter_map(|p| p.fusion_reads.map(|f| [p.reads_passed as f64, f as f64]))
        .collect()
}

/// Has the leading isoform changed between snapshots? Worth flagging: it means the
/// call is not yet stable.
pub fn isoform_changed(history: &[HistoryPoint]) -> bool {
    let calls: Vec<_> = history.iter().filter_map(|p| p.primary_isoform.as_deref()).collect();
    calls.windows(2).any(|w| w[0] != w[1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use nanobcr_core::EngineConfig;
    use std::path::PathBuf;

    fn state(status: &str, heartbeat: u64) -> State {
        let cfg = EngineConfig::new(PathBuf::from("/in"), PathBuf::from("/out"), PathBuf::from("/ref.fa"));
        let mut s = State::new(&cfg);
        s.status = status.into();
        s.heartbeat_unix = heartbeat;
        s
    }

    #[test]
    fn liveness_distinguishes_running_stale_finished_and_error() {
        assert_eq!(liveness(None, 100), Liveness::NotStarted);
        assert_eq!(liveness(Some(&state("running", 95)), 100), Liveness::Running { heartbeat_age: 5 });
        assert_eq!(liveness(Some(&state("idle", 100)), 600), Liveness::Stale { heartbeat_age: 500 });
        assert_eq!(liveness(Some(&state("finished", 0)), 10_000), Liveness::Finished);
        let mut s = state("error", 99);
        s.last_error = Some("boom".into());
        assert_eq!(liveness(Some(&s), 100), Liveness::Error("boom".into()));
    }

    #[test]
    fn formats_ages_and_counts() {
        assert_eq!(format_age(5), "5 s");
        assert_eq!(format_age(150), "2 min");
        assert_eq!(format_age(7200), "2.0 h");
        assert_eq!(format_age(172_800), "2.0 d");
        assert_eq!(format_count(1234567), "1,234,567");
        assert_eq!(format_count(12), "12");
    }

    #[test]
    fn loads_the_engine_directory_and_builds_series() {
        let dir = tempfile::tempdir().unwrap();
        let cfg = EngineConfig::new(dir.path().into(), dir.path().into(), PathBuf::from("/ref.fa"));
        let mut st = State::new(&cfg);
        st.barcode("barcode01").reads_passed = 7;
        st.write(&dir.path().join("state.json")).unwrap();

        let prelim = dir.path().join("barcode01/prelim");
        fs::create_dir_all(&prelim).unwrap();
        let hist = [
            r#"{"unix":1,"generation":1,"reads_passed":100,"status":"ok","primary_isoform":"e13a2","fusion_reads":40,"mutations":[],"enough":false}"#,
            r#"{"unix":2,"generation":2,"reads_passed":300,"status":"ok","primary_isoform":"e14a2","fusion_reads":120,"mutations":[{"mutation":"T315I","vaf_percent":2.5,"depth":900,"status":"DETECTED"}],"enough":true}"#,
            "not json at all",
        ];
        fs::write(prelim.join("history.jsonl"), hist.join("\n")).unwrap();
        fs::write(prelim.join("summary.json"), "{ broken").unwrap();

        let d = load(dir.path());
        assert_eq!(d.state.unwrap().barcodes["barcode01"].reads_passed, 7);
        let h = &d.history["barcode01"];
        assert_eq!(h.len(), 2, "the garbage line is skipped");
        assert_eq!(mutation_series(h)["T315I"], vec![[300.0, 2.5]]);
        assert_eq!(fusion_series(h), vec![[100.0, 40.0], [300.0, 120.0]]);
        assert!(isoform_changed(h));
        assert!(d.problems.iter().any(|p| p.contains("summary.json")), "a broken summary is reported");
    }

    #[test]
    fn an_empty_directory_loads_as_not_started() {
        let dir = tempfile::tempdir().unwrap();
        let d = load(dir.path());
        assert!(d.state.is_none() && d.problems.is_empty());
        assert_eq!(liveness(d.state.as_ref(), 5), Liveness::NotStarted);
    }
}
