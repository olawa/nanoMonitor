//! Handing the engine's alignments to nagini.
//!
//! nagini takes `--input-type aligned_bam` as is: it does not realign or filter,
//! it only stages (read group, sort, index). The BAMs therefore have to be exactly
//! what nagini's own FASTQ path would have produced, which the engine guarantees
//! (same filter predicate, same minimap2 invocation); this module only gathers the
//! final snapshots into the layout nagini's input preparation understands:
//! `<handoff>/barcodeNN.bam` (+`.bai`) and a barcode map CSV.

use crate::state::{now_unix, State};
use crate::store::{compact_out_dir, write_atomic};
use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// A heartbeat older than this means the engine is gone (same rule as the GUI).
pub const ENGINE_STALE_SECS: u64 = 120;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MapRow {
    pub barcode: String,
    pub sample: String,
    /// `T` (tumour), `N` (normal) or `R` — nagini's sample type.
    pub sample_type: String,
    pub assay: String,
    pub library_type: String,
}

fn valid_barcode(b: &str) -> bool {
    b.strip_prefix("barcode")
        .map(|d| !d.is_empty() && d.bytes().all(|c| c.is_ascii_digit()))
        .unwrap_or(false)
}

/// Read a nagini barcode map (CSV or TSV; `barcode`, `sample`, optional `type`,
/// `assay`, `library_type`). Barcodes may be written `barcode01`, `bc01` or `01`.
pub fn read_map(path: &Path) -> Result<Vec<MapRow>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading barcode map {}", path.display()))?;
    let delim = if text.lines().next().map(|l| l.contains('\t')).unwrap_or(false) { '\t' } else { ',' };
    let mut lines = text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let header: Vec<String> = lines
        .next()
        .ok_or_else(|| anyhow!("{} is empty", path.display()))?
        .split(delim)
        .map(|h| h.trim().trim_matches('"').to_ascii_lowercase())
        .collect();
    let col = |name: &str| header.iter().position(|h| h == name);
    let (Some(bc), Some(sm)) = (col("barcode"), col("sample")) else {
        bail!("{}: needs `barcode` and `sample` columns (found: {})", path.display(), header.join(", "));
    };
    let (ty, assay, lib) = (col("type"), col("assay"), col("library_type"));
    let mut rows: Vec<MapRow> = Vec::new();
    for (n, line) in lines.enumerate() {
        let f: Vec<&str> = line.split(delim).map(|s| s.trim().trim_matches('"')).collect();
        let get = |i: Option<usize>, d: &str| i.and_then(|i| f.get(i)).filter(|s| !s.is_empty()).map(|s| s.to_string()).unwrap_or_else(|| d.to_string());
        let raw = get(Some(bc), "").to_ascii_lowercase();
        let barcode = if valid_barcode(&raw) {
            raw
        } else if let Some(d) = raw.strip_prefix("bc").filter(|d| !d.is_empty() && d.bytes().all(|c| c.is_ascii_digit())) {
            format!("barcode{d}")
        } else if !raw.is_empty() && raw.bytes().all(|c| c.is_ascii_digit()) {
            format!("barcode{raw}")
        } else {
            bail!("{}: row {}: `{}` is not a barcode (expected barcode01, bc01 or 01)", path.display(), n + 2, raw);
        };
        let sample = get(Some(sm), "");
        if sample.is_empty() {
            bail!("{}: row {}: {barcode} has no sample name", path.display(), n + 2);
        }
        let row = MapRow {
            barcode,
            sample,
            sample_type: get(ty, "T"),
            assay: get(assay, "bcr_abl"),
            library_type: get(lib, "cdna"),
        };
        if rows.iter().any(|r| r.barcode == row.barcode) {
            bail!("{}: {} appears twice", path.display(), row.barcode);
        }
        if rows.iter().any(|r| r.sample == row.sample && r.sample_type == row.sample_type) {
            bail!("{}: sample {} ({}) appears twice; nagini needs unique sample names", path.display(), row.sample, row.sample_type);
        }
        rows.push(row);
    }
    if rows.is_empty() {
        bail!("{} has a header but no rows", path.display());
    }
    Ok(rows)
}

pub fn write_map(path: &Path, rows: &[MapRow]) -> Result<()> {
    let mut s = String::from("barcode,sample,type,assay,library_type\n");
    for r in rows {
        s.push_str(&format!("{},{},{},{},{}\n", r.barcode, r.sample, r.sample_type, r.assay, r.library_type));
    }
    write_atomic(path, s.as_bytes())
}

/// A starting point for the user to edit: every barcode the engine saw, with a
/// placeholder sample name (the barcode itself).
pub fn template_for(barcodes: &[String]) -> Vec<MapRow> {
    barcodes
        .iter()
        .filter(|b| valid_barcode(b))
        .map(|b| MapRow {
            barcode: b.clone(),
            sample: b.clone(),
            sample_type: "T".into(),
            assay: "bcr_abl".into(),
            library_type: "cdna".into(),
        })
        .collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HandoffBarcode {
    pub barcode: String,
    pub sample: String,
    pub sample_type: String,
    pub bam: PathBuf,
    pub snapshot_generation: u64,
    /// Primary alignments in the BAM (`samtools view -c -F 0x900`).
    pub primary_records: u64,
    /// Reads the engine counted as passing the filter.
    pub reads_passed: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Handoff {
    pub dir: PathBuf,
    pub barcode_map: PathBuf,
    pub written_unix: u64,
    pub barcodes: Vec<HandoffBarcode>,
    /// (barcode, why) for barcodes the engine has data for but that are not handed over.
    pub skipped: Vec<(String, String)>,
    pub warnings: Vec<String>,
}

/// True while an engine is (or may still be) working in `out_dir`.
pub fn engine_active(out_dir: &Path) -> bool {
    let Ok(bytes) = fs::read(out_dir.join("state.json")) else { return false };
    let Ok(state) = serde_json::from_slice::<State>(&bytes) else { return false };
    state.status != "finished" && now_unix().saturating_sub(state.heartbeat_unix) <= ENGINE_STALE_SECS
}

fn samtools_count(samtools: &str, bam: &Path, flags: &[&str]) -> Result<u64> {
    let out = Command::new(samtools).arg("view").arg("-c").args(flags).arg(bam).output()?;
    if !out.status.success() {
        bail!("samtools view -c failed on {}: {}", bam.display(), String::from_utf8_lossy(&out.stderr).trim());
    }
    String::from_utf8_lossy(&out.stdout).trim().parse().map_err(|e| anyhow!("unexpected samtools output: {e}"))
}

fn link_or_copy(from: &Path, to: &Path) -> Result<()> {
    let _ = fs::remove_file(to);
    if fs::hard_link(from, to).is_err() {
        fs::copy(from, to).with_context(|| format!("copying {}", from.display()))?;
    }
    Ok(())
}

/// Fold the last segments into each barcode's snapshot and gather the BAMs that
/// have a sample in `map_path` into `handoff_dir`. The engine must be stopped (a
/// clean stop compacts too); a running engine would still be adding reads.
pub fn prepare(out_dir: &Path, handoff_dir: &Path, map_path: &Path, samtools: &str, threads: usize) -> Result<Handoff> {
    if engine_active(out_dir) {
        bail!("the engine is still running in {}: stop it first (a clean stop does the final compaction)", out_dir.display());
    }
    let map = read_map(map_path)?;
    let by_barcode: BTreeMap<&str, &MapRow> = map.iter().map(|r| (r.barcode.as_str(), r)).collect();
    let state: Option<State> = fs::read(out_dir.join("state.json")).ok().and_then(|b| serde_json::from_slice(&b).ok());

    let bams = compact_out_dir(out_dir, samtools, threads).context("final compaction")?;
    if bams.is_empty() {
        bail!("no alignments in {}: has the engine processed any reads?", out_dir.display());
    }
    fs::create_dir_all(handoff_dir)?;
    // Absolute paths: handoff.json is read later from other places.
    let handoff_dir = &fs::canonicalize(handoff_dir)?;
    let map_path = &fs::canonicalize(map_path)?;
    // Drop stale BAMs from an earlier hand-off: nagini globs the directory.
    for e in fs::read_dir(handoff_dir)?.flatten() {
        let n = e.file_name().to_string_lossy().to_string();
        if n.ends_with(".bam") || n.ends_with(".bam.bai") {
            let _ = fs::remove_file(e.path());
        }
    }

    let mut h = Handoff {
        dir: handoff_dir.to_path_buf(),
        barcode_map: handoff_dir.join("barcode_map.csv"),
        written_unix: now_unix(),
        barcodes: Vec::new(),
        skipped: Vec::new(),
        warnings: Vec::new(),
    };
    for (barcode, bam) in &bams {
        let Some(row) = by_barcode.get(barcode.as_str()) else {
            h.skipped.push((barcode.clone(), "no sample in the barcode map".into()));
            continue;
        };
        let ok = Command::new(samtools).arg("quickcheck").arg(bam).status().map(|s| s.success()).unwrap_or(false);
        if !ok {
            h.skipped.push((barcode.clone(), "samtools quickcheck failed on its BAM".into()));
            continue;
        }
        let primary = samtools_count(samtools, bam, &["-F", "0x900"])?;
        let passed = state.as_ref().and_then(|s| s.barcodes.get(barcode)).map(|b| b.reads_passed).unwrap_or(0);
        if passed != 0 && primary != passed {
            h.warnings.push(format!(
                "{barcode}: the BAM has {primary} primary alignments but the engine counted {passed} reads passing the filter"
            ));
        }
        let dest = handoff_dir.join(format!("{barcode}.bam"));
        link_or_copy(bam, &dest)?;
        link_or_copy(&PathBuf::from(format!("{}.bai", bam.display())), &PathBuf::from(format!("{}.bai", dest.display())))?;
        let generation = crate::store::Store::new(&out_dir.join(barcode)).manifest().generation;
        h.barcodes.push(HandoffBarcode {
            barcode: barcode.clone(),
            sample: row.sample.clone(),
            sample_type: row.sample_type.clone(),
            bam: dest,
            snapshot_generation: generation,
            primary_records: primary,
            reads_passed: passed,
        });
    }
    for r in &map {
        if !bams.contains_key(&r.barcode) {
            h.warnings.push(format!("{} ({}) is in the barcode map but the engine has no reads for it", r.barcode, r.sample));
        }
    }
    if h.barcodes.is_empty() {
        bail!("none of the barcodes with data ({}) is in {}", bams.keys().cloned().collect::<Vec<_>>().join(", "), map_path.display());
    }
    // The map nagini reads: only the barcodes actually handed over.
    let used: Vec<MapRow> = map.into_iter().filter(|r| h.barcodes.iter().any(|b| b.barcode == r.barcode)).collect();
    write_map(&h.barcode_map, &used)?;
    if out_dir.join("provenance.json").is_file() {
        fs::copy(out_dir.join("provenance.json"), handoff_dir.join("provenance.json"))?;
    }
    write_atomic(&handoff_dir.join("handoff.json"), &serde_json::to_vec_pretty(&h)?)?;
    Ok(h)
}

pub fn load(handoff_dir: &Path) -> Option<Handoff> {
    serde_json::from_slice(&fs::read(handoff_dir.join("handoff.json")).ok()?).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
        let p = dir.join(name);
        fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn reads_nagini_style_maps_in_csv_and_tsv_with_loose_barcode_names() {
        let d = tempfile::tempdir().unwrap();
        let csv = write(d.path(), "m.csv", "barcode,sample,type\nbarcode01,PATIENT_A,T\nbc02,PATIENT_B,\n03,PATIENT_C,N\n");
        let rows = read_map(&csv).unwrap();
        assert_eq!(rows.iter().map(|r| r.barcode.as_str()).collect::<Vec<_>>(), ["barcode01", "barcode02", "barcode03"]);
        assert_eq!(rows[1].sample_type, "T", "type defaults to T");
        assert_eq!((rows[0].assay.as_str(), rows[0].library_type.as_str()), ("bcr_abl", "cdna"));
        assert_eq!(rows[2].sample_type, "N");

        let tsv = write(d.path(), "m.tsv", "Barcode\tSample\nbarcode07\tX\n");
        assert_eq!(read_map(&tsv).unwrap()[0].sample, "X");
    }

    #[test]
    fn rejects_maps_nagini_would_reject() {
        let d = tempfile::tempdir().unwrap();
        for (body, needle) in [
            ("barcode,sample\nbarcode01,A\nbarcode01,B\n", "appears twice"),
            ("barcode,sample\nbarcode01,A\nbarcode02,A\n", "unique sample"),
            ("barcode,sample\nfoo,A\n", "not a barcode"),
            ("barcode,sample\nbarcode01,\n", "no sample"),
            ("code,name\nx,y\n", "needs `barcode` and `sample`"),
            ("barcode,sample\n", "no rows"),
        ] {
            let p = write(d.path(), "bad.csv", body);
            let e = read_map(&p).unwrap_err().to_string();
            assert!(e.contains(needle), "{body:?} -> {e}");
        }
    }

    #[test]
    fn template_round_trips_and_skips_unclassified() {
        let rows = template_for(&["barcode02".into(), "unclassified".into(), "barcode01".into()]);
        assert_eq!(rows.len(), 2);
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("t.csv");
        write_map(&p, &rows).unwrap();
        assert_eq!(read_map(&p).unwrap(), rows);
    }

    #[test]
    fn a_fresh_or_missing_state_decides_whether_the_engine_is_active() {
        let d = tempfile::tempdir().unwrap();
        assert!(!engine_active(d.path()));
        let cfg = crate::EngineConfig::new(d.path().into(), d.path().into(), PathBuf::from("/r.fa"));
        let mut s = State::new(&cfg);
        s.status = "idle".into();
        s.write(&d.path().join("state.json")).unwrap();
        assert!(engine_active(d.path()));
        s.status = "finished".into();
        s.write(&d.path().join("state.json")).unwrap();
        assert!(!engine_active(d.path()));
        s.status = "running".into();
        s.write(&d.path().join("state.json")).unwrap();
        // Age the heartbeat past the stale limit.
        let mut v: serde_json::Value = serde_json::from_slice(&fs::read(d.path().join("state.json")).unwrap()).unwrap();
        v["heartbeat_unix"] = (now_unix() - ENGINE_STALE_SECS - 5).into();
        fs::write(d.path().join("state.json"), v.to_string()).unwrap();
        assert!(!engine_active(d.path()));
    }
}
