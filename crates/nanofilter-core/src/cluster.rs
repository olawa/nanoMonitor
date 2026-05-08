//! UMI family clustering, strand balance, and downsampling.
//!
//! Reads with a detected UMI key are grouped into families using either
//! exact key matching or a greedy star-clustering approximation that tolerates
//! sequencing errors in the UMI string.

use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::umi::{FullReadUmiRecord, Strand};

// ---------------------------------------------------------------------------
// Clustering mode
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub enum ClusterMode {
    /// Keys must match exactly.
    Exact,
    /// Greedy star-clustering: two UMIs are in the same family if their edit
    /// distance is `<= max_edit`.  Conceptually equivalent to vsearch
    /// `--cluster_fast` at ~0.85 identity on short UMI strings.
    Approximate { max_edit: usize },
    /// Shell out to vsearch for clustering (requires vsearch on PATH).
    Vsearch { identity: f64 },
}

// ---------------------------------------------------------------------------
// Core data structures
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct FamilyRead {
    pub read_id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
    pub strand: Option<Strand>,
}

#[derive(Debug, Clone)]
pub struct UmiFamily {
    /// 0-based family index.
    pub id: usize,
    /// Representative UMI string (first cluster centre).
    pub consensus_umi: String,
    pub reads: Vec<FamilyRead>,
    pub n_fwd: usize,
    pub n_rev: usize,
}

impl UmiFamily {
    pub fn new(id: usize, consensus_umi: String) -> Self {
        Self {
            id,
            consensus_umi,
            reads: Vec::new(),
            n_fwd: 0,
            n_rev: 0,
        }
    }

    pub fn push(&mut self, fr: FamilyRead) {
        match fr.strand {
            Some(Strand::Fwd) => self.n_fwd += 1,
            Some(Strand::Rev) => self.n_rev += 1,
            None => {}
        }
        self.reads.push(fr);
    }

    pub fn total_reads(&self) -> usize {
        self.reads.len()
    }
}

// ---------------------------------------------------------------------------
// Clustering
// ---------------------------------------------------------------------------

/// Group `records` (those that have a cluster key) into UMI families.
pub fn cluster_by_umi(records: &[FullReadUmiRecord], mode: &ClusterMode) -> Vec<UmiFamily> {
    // Collect only reads that have a UMI key
    let keyed: Vec<(&FullReadUmiRecord, &str)> = records
        .iter()
        .filter_map(|r| r.cluster_key().map(|k| (r, k)))
        .collect();

    match mode {
        ClusterMode::Exact => cluster_exact(&keyed),
        ClusterMode::Approximate { max_edit } => cluster_approximate(&keyed, *max_edit),
        ClusterMode::Vsearch { identity } => cluster_vsearch(&keyed, *identity),
    }
}

fn make_family_read(rec: &FullReadUmiRecord) -> FamilyRead {
    FamilyRead {
        read_id: rec.read_id().to_string(),
        seq: rec.seq.clone(),
        qual: rec.qual.clone(),
        strand: rec.strand(),
    }
}

fn cluster_exact(keyed: &[(&FullReadUmiRecord, &str)]) -> Vec<UmiFamily> {
    use std::collections::HashMap;
    let mut map: HashMap<&str, usize> = HashMap::new();
    let mut families: Vec<UmiFamily> = Vec::new();

    for &(rec, key) in keyed {
        let idx = map.entry(key).or_insert_with(|| {
            let id = families.len();
            families.push(UmiFamily::new(id, key.to_string()));
            id
        });
        families[*idx].push(make_family_read(rec));
    }
    families
}

/// Greedy O(n²) star-clustering.  For small UMI strings (≤ ~75 bases) this is
/// fast enough: we use the IUPAC-aware edit_distance from `barcode.rs` via a
/// simple loop, but with a short-circuit on the first good match.
fn cluster_approximate(keyed: &[(&FullReadUmiRecord, &str)], max_edit: usize) -> Vec<UmiFamily> {
    use crate::barcode::edit_distance;

    let mut families: Vec<UmiFamily> = Vec::new();
    // centres: Vec<Vec<u8>>
    let mut centres: Vec<Vec<u8>> = Vec::new();

    for &(rec, key) in keyed {
        let key_bytes = key.as_bytes();
        let mut assigned = None;

        for (fi, centre) in centres.iter().enumerate() {
            let d = edit_distance(key_bytes, centre, max_edit);
            if d <= max_edit {
                assigned = Some(fi);
                break;
            }
        }

        match assigned {
            Some(fi) => families[fi].push(make_family_read(rec)),
            None => {
                let id = families.len();
                let mut fam = UmiFamily::new(id, key.to_string());
                fam.push(make_family_read(rec));
                families.push(fam);
                centres.push(key_bytes.to_vec());
            }
        }
    }
    families
}

/// Shell out to vsearch for clustering; writes a temporary FASTA of UMI
/// sequences, runs vsearch `--cluster_fast`, reads back the UC file and
/// groups reads accordingly.  Requires vsearch on PATH.
fn cluster_vsearch(keyed: &[(&FullReadUmiRecord, &str)], identity: f64) -> Vec<UmiFamily> {
    use std::collections::HashMap;
    use std::process::Command;
    use std::time::{SystemTime, UNIX_EPOCH};

    // Write tmp FASTA of UMI strings
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or_default();
    let tmp_dir = std::env::temp_dir();
    let tmp_fasta = tmp_dir.join(format!(
        "nanofilter_umi_cluster_{}_{}.fasta",
        std::process::id(),
        unique
    ));
    let tmp_uc = tmp_dir.join(format!(
        "nanofilter_umi_cluster_{}_{}.uc",
        std::process::id(),
        unique
    ));
    {
        let mut f = BufWriter::new(File::create(&tmp_fasta).expect("Failed to create tmp FASTA"));
        for (i, &(_rec, key)) in keyed.iter().enumerate() {
            writeln!(f, ">umi_{}\n{}", i, key).unwrap();
        }
    }

    let status = Command::new("vsearch")
        .args([
            "--cluster_fast",
            tmp_fasta.to_str().unwrap(),
            "--id",
            &identity.to_string(),
            "--uc",
            tmp_uc.to_str().unwrap(),
            "--quiet",
        ])
        .status();

    if status.map(|s| !s.success()).unwrap_or(true) {
        eprintln!(
            "vsearch clustering failed or vsearch not on PATH; falling back to approximate mode"
        );
        let _ = std::fs::remove_file(&tmp_fasta);
        let _ = std::fs::remove_file(&tmp_uc);
        return cluster_approximate(keyed, 3);
    }

    // Parse UC file: fields 0=type (H=hit,S=seed), 1=cluster_idx, 8=read_label
    let uc_content = std::fs::read_to_string(&tmp_uc).unwrap_or_default();
    let mut cluster_idx_to_fam: HashMap<usize, usize> = HashMap::new();
    let mut families: Vec<UmiFamily> = Vec::new();

    for line in uc_content.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let rec_type = fields[0];
        if rec_type != "H" && rec_type != "S" {
            continue;
        }
        let cluster_num: usize = fields[1].parse().unwrap_or(0);
        let label = fields[8]; // "umi_N"
        let read_idx: usize = label
            .strip_prefix("umi_")
            .and_then(|n| n.parse().ok())
            .unwrap_or(0);

        if read_idx >= keyed.len() {
            continue;
        }
        let (rec, key) = keyed[read_idx];

        let fam_idx = *cluster_idx_to_fam.entry(cluster_num).or_insert_with(|| {
            let id = families.len();
            families.push(UmiFamily::new(id, key.to_string()));
            id
        });
        families[fam_idx].push(make_family_read(rec));
    }

    let _ = std::fs::remove_file(&tmp_fasta);
    let _ = std::fs::remove_file(&tmp_uc);

    families
}

// ---------------------------------------------------------------------------
// Filter and downsample
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct ClusterFilterConfig {
    pub min_reads: usize,
    pub max_reads: usize,
    pub balance_strands: bool,
}

impl Default for ClusterFilterConfig {
    fn default() -> Self {
        Self {
            min_reads: 4,
            max_reads: 80,
            balance_strands: false,
        }
    }
}

pub struct FilteredFamilies {
    pub passing: Vec<UmiFamily>,
    pub failing: Vec<UmiFamily>,
}

/// Apply `min_reads` threshold and downsample oversized families.
pub fn filter_and_downsample(
    families: Vec<UmiFamily>,
    config: &ClusterFilterConfig,
) -> FilteredFamilies {
    let mut passing = Vec::new();
    let mut failing = Vec::new();

    for mut fam in families {
        if fam.total_reads() < config.min_reads {
            failing.push(fam);
            continue;
        }

        if fam.total_reads() > config.max_reads {
            fam = downsample(fam, config.max_reads, config.balance_strands);
        }
        passing.push(fam);
    }
    FilteredFamilies { passing, failing }
}

fn downsample(mut fam: UmiFamily, max_reads: usize, balance: bool) -> UmiFamily {
    if !balance {
        // Simple truncation (reads already in arrival order)
        fam.reads.truncate(max_reads);
    } else {
        let cap_per_strand = max_reads / 2;
        let mut fwd_seen = 0;
        let mut rev_seen = 0;
        fam.reads.retain(|r| match r.strand {
            Some(Strand::Fwd) => {
                if fwd_seen < cap_per_strand {
                    fwd_seen += 1;
                    true
                } else {
                    false
                }
            }
            Some(Strand::Rev) => {
                if rev_seen < cap_per_strand {
                    rev_seen += 1;
                    true
                } else {
                    false
                }
            }
            None => fwd_seen + rev_seen < max_reads,
        });
    }
    // Recount after downsample
    fam.n_fwd = fam
        .reads
        .iter()
        .filter(|r| r.strand == Some(Strand::Fwd))
        .count();
    fam.n_rev = fam
        .reads
        .iter()
        .filter(|r| r.strand == Some(Strand::Rev))
        .count();
    fam
}

// ---------------------------------------------------------------------------
// Write per-family FASTQ
// ---------------------------------------------------------------------------

/// Write all reads of `family` to `<output_dir>/family_<id>.fastq`.
/// Returns the path of the created file.
pub fn write_family_reads(family: &UmiFamily, output_dir: &Path) -> Result<PathBuf> {
    let path = output_dir.join(format!("family_{:06}.fastq", family.id));
    let mut writer = BufWriter::new(File::create(&path)?);

    for r in &family.reads {
        writer.write_all(b"@")?;
        writer.write_all(r.read_id.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&r.seq)?;
        writer.write_all(b"\n+\n")?;
        match &r.qual {
            Some(q) => writer.write_all(q)?,
            None => writer.write_all(&vec![b'I'; r.seq.len()])?,
        }
        writer.write_all(b"\n")?;
    }
    Ok(path)
}

// ---------------------------------------------------------------------------
// Per-family stats
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct FamilyStats {
    pub id_cluster: usize,
    pub consensus_umi: String,
    pub n_cluster_consensus: usize, // total reads in family before downsampling
    pub min_reads_required: usize,
    pub max_reads_allowed: usize,
    pub balance_strands: bool,
    pub n_fwd: usize,
    pub n_rev: usize,
    pub written_fwd: usize,
    pub written_rev: usize,
    /// Total reads in family after downsampling.
    pub n: usize,
    /// Reads skipped/trimmed.
    pub skipped: usize,
    /// Reads written (= n after downsampling, same as n unless trimmed).
    pub written: usize,
    /// 1 if family passed min_reads threshold, 0 otherwise.
    pub cluster_written: u8,
}

/// Build per-family stats rows from the filtered pipeline outputs.
///
/// Derives "all families" by chaining `passing_families` and `failing_families`,
/// eliminating the need for a separate pre-filter copy of the family list.
pub fn compute_family_stats(
    passing_families: &[UmiFamily],
    failing_families: &[UmiFamily],
    config: &ClusterFilterConfig,
) -> Vec<FamilyStats> {
    use std::collections::HashSet;
    let passing_ids: HashSet<usize> = passing_families.iter().map(|f| f.id).collect();

    // Lookup for post-downsampled counts.
    let post_lookup: std::collections::HashMap<usize, &UmiFamily> =
        passing_families.iter().map(|f| (f.id, f)).collect();

    let mut stats = Vec::new();

    // All families = passing ∪ failing; no separate Vec needed.
    for orig in passing_families.iter().chain(failing_families.iter()) {
        let cluster_written = if passing_ids.contains(&orig.id) {
            1u8
        } else {
            0u8
        };
        let (written_n, written_fwd, written_rev) = if let Some(post) = post_lookup.get(&orig.id) {
            (post.reads.len(), post.n_fwd, post.n_rev)
        } else {
            (0, 0, 0)
        };
        let skipped = orig.reads.len().saturating_sub(written_n);

        stats.push(FamilyStats {
            id_cluster: orig.id,
            consensus_umi: orig.consensus_umi.clone(),
            n_cluster_consensus: orig.reads.len(),
            min_reads_required: config.min_reads,
            max_reads_allowed: config.max_reads,
            balance_strands: config.balance_strands,
            n_fwd: orig.n_fwd,
            n_rev: orig.n_rev,
            written_fwd,
            written_rev,
            n: written_n,
            skipped,
            written: written_n,
            cluster_written,
        });
    }

    stats.sort_by_key(|s| s.id_cluster);
    stats
}
