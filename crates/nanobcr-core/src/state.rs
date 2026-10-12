//! `state.json` (what the GUI reads) and `provenance.json` (what a later nagini
//! run can be checked against).

use crate::align::{tool_version, AlignCounts};
use crate::config::EngineConfig;
use crate::reads::ReadCounts;
use crate::store::write_atomic;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

pub const STATE_SCHEMA: u32 = 1;

pub fn now_unix() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BarcodeState {
    pub files: u64,
    pub reads_seen: u64,
    pub reads_passed: u64,
    pub reads_failed_len: u64,
    pub reads_failed_qv: u64,
    pub bases_passed: u64,
    pub records: u64,
    pub primary: u64,
    pub primary_mapped: u64,
    pub primary_unmapped: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub segments_pending: usize,
    /// Latest sorted, indexed BAM (absolute path), once a compaction has run.
    pub compacted_bam: Option<String>,
    pub compacted_generation: u64,
    pub compacted_unix: Option<u64>,
    pub last_file: Option<String>,
    pub last_update_unix: u64,
}

impl BarcodeState {
    pub fn add_reads(&mut self, c: &ReadCounts) {
        self.reads_seen += c.seen;
        self.reads_passed += c.passed;
        self.reads_failed_len += c.failed_len;
        self.reads_failed_qv += c.failed_qv;
        self.bases_passed += c.bases_passed;
    }

    pub fn add_align(&mut self, c: &AlignCounts) {
        self.records += c.records;
        self.primary += c.primary;
        self.primary_mapped += c.primary_mapped;
        self.primary_unmapped += c.primary_unmapped;
        self.supplementary += c.supplementary;
        self.secondary += c.secondary;
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct State {
    pub schema: u32,
    /// `starting`, `running`, `idle`, `finished` or `error`.
    pub status: String,
    pub pid: u32,
    pub started_unix: u64,
    pub heartbeat_unix: u64,
    pub last_error: Option<String>,
    pub input_dir: String,
    pub out_dir: String,
    pub files_done: u64,
    pub cycles: u64,
    pub barcodes: BTreeMap<String, BarcodeState>,
}

impl State {
    pub fn new(cfg: &EngineConfig) -> Self {
        let now = now_unix();
        Self {
            schema: STATE_SCHEMA,
            status: "starting".into(),
            pid: std::process::id(),
            started_unix: now,
            heartbeat_unix: now,
            last_error: None,
            input_dir: cfg.input_dir.to_string_lossy().into(),
            out_dir: cfg.out_dir.to_string_lossy().into(),
            files_done: 0,
            cycles: 0,
            barcodes: BTreeMap::new(),
        }
    }

    /// Resume the counters of an earlier run in the same output directory.
    pub fn resume(cfg: &EngineConfig, path: &Path) -> Self {
        let fresh = Self::new(cfg);
        let Ok(bytes) = std::fs::read(path) else {
            return fresh;
        };
        match serde_json::from_slice::<State>(&bytes) {
            Ok(prev) if prev.schema == STATE_SCHEMA => Self {
                files_done: prev.files_done,
                cycles: prev.cycles,
                barcodes: prev.barcodes,
                ..fresh
            },
            _ => fresh,
        }
    }

    pub fn barcode(&mut self, name: &str) -> &mut BarcodeState {
        self.barcodes.entry(name.to_string()).or_default()
    }

    pub fn write(&mut self, path: &Path) -> Result<()> {
        self.heartbeat_unix = now_unix();
        self.pid = std::process::id();
        write_atomic(path, &serde_json::to_vec_pretty(self)?)
    }
}

/// FNV-1a: a stable, dependency-free fingerprint (not a security hash) of the
/// reference's `.fai`, enough to tell two builds of a genome apart.
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

pub fn write_provenance(cfg: &EngineConfig, index: &Path, path: &Path) -> Result<()> {
    let fai = {
        let mut p = cfg.reference.clone().into_os_string();
        p.push(".fai");
        std::fs::read(std::path::PathBuf::from(p)).ok()
    };
    let doc = serde_json::json!({
        "schema": 1,
        "written_unix": now_unix(),
        "nanobcr_core": env!("CARGO_PKG_VERSION"),
        "minimap2": tool_version(&cfg.minimap2).unwrap_or_default(),
        "samtools": tool_version(&cfg.samtools).unwrap_or_default(),
        "alignment": {
            "preset": cfg.preset,
            "extra": cfg.minimap2_extra,
            "reference": cfg.reference,
            "index": index,
            "reference_fai_fnv1a": fai.as_deref().map(|b| format!("{:016x}", fnv1a(b))),
        },
        "read_filter": {
            "min_qv": cfg.min_qv,
            "min_len": cfg.min_len,
            "max_len": cfg.max_len,
            "include_fastq_fail": cfg.include_fail,
            "qv_source": "dorado qs header tag, else recomputed from the quality string",
        },
        "input_dir": cfg.input_dir,
    });
    write_atomic(path, &serde_json::to_vec_pretty(&doc)?)
}
