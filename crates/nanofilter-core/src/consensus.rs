//! Consensus derivation for UMI families.
//!
//! Backends:
//! - `None`        – no consensus; per-family FASTQs only.
//! - `Medoid`      – pick the read whose UMI has minimum total edit distance
//!                   to all others (zero-dependency, fast, ~O(n²·L_umi)).
//! - `MajorityVote`– banded Needleman-Wunsch alignment of all reads to the
//!                   medoid anchor, then column-wise plurality vote.
//!                   Pure Rust, no external tools (~O(n·L·band)).
//! - `Medaka`      – shell out to `medaka smolecule`.
//! - `Dorado`      – shell out to `dorado polish` (future; stubbed).

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::cluster::UmiFamily;

// ---------------------------------------------------------------------------
// Backend enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub enum ConsensusBackend {
    None,
    Medoid,
    MajorityVote {
        /// Half-width of the alignment band (default 150 works well for ONT
        /// amplicons with <5% error rate).
        band: usize,
    },
    Medaka {
        model: Option<String>,
        /// `medaka smolecule` `--length` arg (min read length accepted).
        min_length: usize,
        /// `medaka smolecule` `--depth` arg.
        min_depth: usize,
    },
    Dorado,
}

impl Default for ConsensusBackend {
    fn default() -> Self {
        ConsensusBackend::None
    }
}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct ConsensusResult {
    pub family_id: usize,
    /// The consensus (or representative) sequence, if derivation succeeded.
    pub sequence: Option<Vec<u8>>,
    /// Optional per-base quality string aligned to `sequence`.
    pub quality: Option<Vec<u8>>,
    /// Human-readable description of the method used.
    pub method: String,
    pub success: bool,
    pub error: Option<String>,
    /// Path to temporary per-family FASTQ written for this consensus call.
    pub family_fastq: Option<PathBuf>,
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Derive consensus for `family` using the requested `backend`.
/// `family_fastq_path` is the on-disk FASTQ of family reads (already written
/// by `cluster::write_family_reads`).
pub fn derive_consensus(
    family: &UmiFamily,
    family_fastq_path: Option<&Path>,
    backend: &ConsensusBackend,
    output_dir: &Path,
) -> ConsensusResult {
    match backend {
        ConsensusBackend::None => ConsensusResult {
            family_id: family.id,
            sequence: None,
            quality: None,
            method: "none".to_string(),
            success: true,
            error: None,
            family_fastq: family_fastq_path.map(|p| p.to_path_buf()),
        },

        ConsensusBackend::Medoid => medoid_consensus(family),

        ConsensusBackend::MajorityVote { band } => majority_vote_consensus(family, *band),

        ConsensusBackend::Medaka {
            model,
            min_length,
            min_depth,
        } => medaka_consensus(
            family,
            family_fastq_path,
            model.as_deref(),
            *min_length,
            *min_depth,
            output_dir,
        ),

        ConsensusBackend::Dorado => ConsensusResult {
            family_id: family.id,
            sequence: None,
            quality: None,
            method: "dorado (not yet implemented)".to_string(),
            success: false,
            error: Some("Dorado backend is not yet implemented".to_string()),
            family_fastq: family_fastq_path.map(|p| p.to_path_buf()),
        },
    }
}

/// Append one consensus record to a FASTQ file.
///
/// Consensus backends in this crate do not currently compute per-base quality,
/// so the output uses synthetic `I` qualities to remain FASTQ-compatible.
pub fn append_consensus_fastq_record<W: Write>(
    writer: &mut W,
    family_id: usize,
    method: &str,
    sequence: &[u8],
    quality: Option<&[u8]>,
) -> anyhow::Result<()> {
    let qual = quality.filter(|q| q.len() == sequence.len());
    writeln!(
        writer,
        "@family_{:06} method={} length={}",
        family_id,
        method,
        sequence.len()
    )?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n+\n")?;
    if let Some(q) = qual {
        writer.write_all(q)?;
    } else {
        writer.write_all(&vec![b'I'; sequence.len()])?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Medoid backend
// ---------------------------------------------------------------------------

/// Pick the read whose sequence is most central to all others.
/// Uses the UMI strings (≤75 chars) for pairwise distance to keep it cheap.
/// Falls back to the longest read if the UMI isn't available.
fn medoid_consensus(family: &UmiFamily) -> ConsensusResult {
    let reads = &family.reads;
    if reads.is_empty() {
        return ConsensusResult {
            family_id: family.id,
            sequence: None,
            quality: None,
            method: "medoid".to_string(),
            success: false,
            error: Some("Family has no reads".to_string()),
            family_fastq: None,
        };
    }
    if reads.len() == 1 {
        return ConsensusResult {
            family_id: family.id,
            sequence: Some(reads[0].seq.clone()),
            quality: reads[0].qual.clone(),
            method: "medoid (singleton)".to_string(),
            success: true,
            error: None,
            family_fastq: None,
        };
    }

    // Compute pairwise distances using seq slices capped at 75 chars
    // (avoids the 64-char limit in the bit-parallel edit_distance by using
    // our own simple DP here).
    let sigs: Vec<&[u8]> = reads
        .iter()
        .map(|r| {
            let end = r.seq.len().min(75);
            &r.seq[..end]
        })
        .collect();

    let n = sigs.len();
    let mut total_dists = vec![0usize; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let d = simple_edit_dist_capped(sigs[i], sigs[j], usize::MAX);
            total_dists[i] += d;
            total_dists[j] += d;
        }
    }

    let best = total_dists
        .iter()
        .enumerate()
        .min_by_key(|(_, &d)| d)
        .map(|(i, _)| i)
        .unwrap_or(0);

    ConsensusResult {
        family_id: family.id,
        sequence: Some(reads[best].seq.clone()),
        quality: reads[best].qual.clone(),
        method: format!("medoid (n={})", n),
        success: true,
        error: None,
        family_fastq: None,
    }
}

// ---------------------------------------------------------------------------
// MajorityVote backend
// ---------------------------------------------------------------------------

fn majority_vote_consensus(family: &UmiFamily, band: usize) -> ConsensusResult {
    let reads = &family.reads;
    if reads.is_empty() {
        return ConsensusResult {
            family_id: family.id,
            sequence: None,
            quality: None,
            method: "majority_vote".to_string(),
            success: false,
            error: Some("Family has no reads".to_string()),
            family_fastq: None,
        };
    }
    if reads.len() == 1 {
        return ConsensusResult {
            family_id: family.id,
            sequence: Some(reads[0].seq.clone()),
            quality: reads[0].qual.clone(),
            method: "majority_vote (singleton)".to_string(),
            success: true,
            error: None,
            family_fastq: None,
        };
    }

    // Use the longest read as anchor reference
    let anchor_idx = reads
        .iter()
        .enumerate()
        .max_by_key(|(_, r)| r.seq.len())
        .map(|(i, _)| i)
        .unwrap_or(0);
    let anchor = &reads[anchor_idx].seq;
    let ref_len = anchor.len();

    // Column base-count table: counts[col][base_idx]
    // base_idx: A=0, C=1, G=2, T=3, gap=4
    let mut counts: Vec<[u32; 5]> = vec![[0; 5]; ref_len];

    // Contribute anchor itself
    for (col, &b) in anchor.iter().enumerate() {
        counts[col][base_to_idx(b)] += 1;
    }

    // Align each other read to anchor via banded NW
    for (i, read) in reads.iter().enumerate() {
        if i == anchor_idx {
            continue;
        }
        let aligned = banded_nw_align(anchor, &read.seq, band);
        // `aligned` maps anchor column → read base (or gap=b'-')
        for (col, b) in aligned.into_iter().enumerate() {
            if col < ref_len {
                if b == b'-' {
                    counts[col][4] += 1;
                } else {
                    counts[col][base_to_idx(b)] += 1;
                }
            }
        }
    }

    // Emit plurality base at each column (gaps only win if a majority),
    // and derive a simple confidence score from the winning fraction.
    const BASES: [u8; 5] = [b'A', b'C', b'G', b'T', b'-'];
    let mut consensus: Vec<u8> = Vec::with_capacity(ref_len);
    let mut qualities: Vec<u8> = Vec::with_capacity(ref_len);
    let total = reads.len();
    for col in &counts {
        let (best_idx, &best_count) = col.iter().enumerate().max_by_key(|(_, &v)| v).unwrap();
        let b = BASES[best_idx];
        if b != b'-' {
            consensus.push(b);
            qualities.push(phred_from_support(best_count as usize, total));
        }
    }

    ConsensusResult {
        family_id: family.id,
        sequence: Some(consensus),
        quality: Some(qualities),
        method: format!("majority_vote (n={}, band={})", reads.len(), band),
        success: true,
        error: None,
        family_fastq: None,
    }
}

#[inline]
fn base_to_idx(b: u8) -> usize {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum TraceMove {
    Diag,
    Up,
    Left,
}

#[derive(Debug)]
struct TraceRow {
    j_lo: usize,
    moves: Vec<TraceMove>,
}

impl TraceRow {
    fn j_hi(&self) -> usize {
        self.j_lo + self.moves.len().saturating_sub(1)
    }
}

/// Banded Needleman-Wunsch alignment of `query` against `reference`.
/// Returns a Vec of length `reference.len()` where each element is the
/// aligned query base at that reference column (or b'-' for a gap).
fn banded_nw_align(reference: &[u8], query: &[u8], band: usize) -> Vec<u8> {
    let r = reference.len();
    let q = query.len();

    if r == 0 {
        return Vec::new();
    }

    // If the band cannot even cover the length difference, fall back to the
    // original full-matrix alignment to preserve correctness for outliers.
    if r.abs_diff(q) > band {
        return full_nw_align(reference, query);
    }

    const GAP: i16 = -2;
    const MATCH: i16 = 1;
    const MISMATCH: i16 = -1;
    const NEG_INF: i16 = i16::MIN / 2;

    let mut prev = vec![NEG_INF; q + 1];
    let mut curr = vec![NEG_INF; q + 1];
    let mut trace: Vec<TraceRow> = Vec::with_capacity(r);

    // Row 0: pure insertions in the query.
    prev[0] = 0;
    for j in 1..=q.min(band) {
        prev[j] = prev[j - 1] + GAP;
    }

    for i in 1..=r {
        for cell in &mut curr {
            *cell = NEG_INF;
        }

        let j_lo = i.saturating_sub(band);
        let j_hi = q.min(i + band);
        if j_lo > j_hi {
            return full_nw_align(reference, query);
        }

        let mut moves = vec![TraceMove::Diag; j_hi - j_lo + 1];

        if j_lo == 0 {
            if i <= band {
                curr[0] = (i as i16) * GAP;
                moves[0] = TraceMove::Up;
            }
        }

        let start_j = j_lo.max(1);
        for j in start_j..=j_hi {
            let mut best = NEG_INF;
            let mut best_move = TraceMove::Diag;

            let diag_score = prev[j - 1];
            if diag_score != NEG_INF {
                let score = diag_score
                    + if reference[i - 1].to_ascii_uppercase() == query[j - 1].to_ascii_uppercase()
                    {
                        MATCH
                    } else {
                        MISMATCH
                    };
                if score > best {
                    best = score;
                    best_move = TraceMove::Diag;
                }
            }

            let up_score = prev[j];
            if up_score != NEG_INF {
                let score = up_score + GAP;
                if score > best {
                    best = score;
                    best_move = TraceMove::Up;
                }
            }

            let left_score = curr[j - 1];
            if left_score != NEG_INF {
                let score = left_score + GAP;
                if score > best {
                    best = score;
                    best_move = TraceMove::Left;
                }
            }

            curr[j] = best;
            moves[j - j_lo] = best_move;
        }

        trace.push(TraceRow { j_lo, moves });
        std::mem::swap(&mut prev, &mut curr);
    }

    if prev[q] == NEG_INF {
        return full_nw_align(reference, query);
    }

    // Traceback over the banded directions.
    let mut aligned = vec![b'-'; r];
    let mut i = r;
    let mut j = q;

    while i > 0 {
        if j == 0 {
            aligned[i - 1] = b'-';
            i -= 1;
            continue;
        }

        let row = &trace[i - 1];
        if j < row.j_lo || j > row.j_hi() {
            return full_nw_align(reference, query);
        }

        match row.moves[j - row.j_lo] {
            TraceMove::Diag => {
                aligned[i - 1] = query[j - 1];
                i -= 1;
                j -= 1;
            }
            TraceMove::Up => {
                aligned[i - 1] = b'-';
                i -= 1;
            }
            TraceMove::Left => {
                j -= 1;
            }
        }
    }

    aligned
}

/// Original full-matrix alignment retained as a correctness fallback.
fn full_nw_align(reference: &[u8], query: &[u8]) -> Vec<u8> {
    let r = reference.len();
    let q = query.len();

    // Full score matrix: (r+1) x (q+1).
    // Use i16 for memory efficiency; gap = -2, match = 1, mismatch = -1
    const GAP: i16 = -2;
    const MATCH: i16 = 1;
    const MISMATCH: i16 = -1;
    const NEG_INF: i16 = i16::MIN / 2;

    let mut dp = vec![vec![NEG_INF; q + 1]; r + 1];

    // Initialise
    for i in 0..=r {
        dp[i][0] = -(i as i16) * (-GAP);
    }
    for j in 0..=q {
        dp[0][j] = -(j as i16) * (-GAP);
    }

    for i in 1..=r {
        for j in 1..=q {
            let diag_score = dp[i - 1][j - 1];
            let cost = if reference[i - 1].to_ascii_uppercase() == query[j - 1].to_ascii_uppercase()
            {
                MATCH
            } else {
                MISMATCH
            };
            let from_diag = if diag_score == NEG_INF {
                NEG_INF
            } else {
                diag_score + cost
            };
            let from_up = if dp[i - 1][j] == NEG_INF {
                NEG_INF
            } else {
                dp[i - 1][j] + GAP
            };
            let from_left = if dp[i][j - 1] == NEG_INF {
                NEG_INF
            } else {
                dp[i][j - 1] + GAP
            };
            dp[i][j] = from_diag.max(from_up).max(from_left);
        }
    }

    // Traceback
    let mut aligned = vec![b'-'; r];
    let (mut i, mut j) = (r, q);

    while i > 0 && j > 0 {
        let score = dp[i][j];
        let up = if i > 0 && dp[i - 1][j] != NEG_INF {
            dp[i - 1][j] + GAP
        } else {
            NEG_INF
        };
        let left = if j > 0 && dp[i][j - 1] != NEG_INF {
            dp[i][j - 1] + GAP
        } else {
            NEG_INF
        };
        let diag_cost =
            if reference[i - 1].to_ascii_uppercase() == query[j - 1].to_ascii_uppercase() {
                MATCH
            } else {
                MISMATCH
            };
        let diag = if i > 0 && j > 0 && dp[i - 1][j - 1] != NEG_INF {
            dp[i - 1][j - 1] + diag_cost
        } else {
            NEG_INF
        };

        if score == diag {
            aligned[i - 1] = query[j - 1];
            i -= 1;
            j -= 1;
        } else if score == up {
            aligned[i - 1] = b'-'; // gap in query
            i -= 1;
        } else if score == left {
            j -= 1; // insertion in query (skip, doesn't map to ref column)
        } else {
            // Fallback: treat as diagonal
            aligned[i - 1] = query[j - 1];
            i -= 1;
            j -= 1;
        }
    }

    aligned
}

// ---------------------------------------------------------------------------
// Medaka backend
// ---------------------------------------------------------------------------

fn medaka_consensus(
    family: &UmiFamily,
    family_fastq_path: Option<&Path>,
    model: Option<&str>,
    min_length: usize,
    min_depth: usize,
    _output_dir: &Path,
) -> ConsensusResult {
    let fastq = match family_fastq_path {
        Some(p) => p,
        None => {
            return ConsensusResult {
                family_id: family.id,
                sequence: None,
                quality: None,
                method: "medaka".to_string(),
                success: false,
                error: Some("No family FASTQ path provided to medaka backend".to_string()),
                family_fastq: None,
            };
        }
    };

    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or_default();
    let consensus_out = std::env::temp_dir().join(format!(
        "nanofilter_medaka_{}_{}_{}",
        std::process::id(),
        family.id,
        unique
    ));
    if let Err(e) = std::fs::create_dir_all(&consensus_out) {
        return ConsensusResult {
            family_id: family.id,
            sequence: None,
            quality: None,
            method: "medaka".to_string(),
            success: false,
            error: Some(format!("Failed to create medaka temp dir: {}", e)),
            family_fastq: Some(fastq.to_path_buf()),
        };
    }

    let mut cmd = Command::new("medaka");
    cmd.args(["smolecule", "--threads", "1"])
        .args(["--length", &min_length.to_string()])
        .args(["--depth", &min_depth.to_string()])
        .args(["--method", "spoa"]);
    if let Some(m) = model {
        cmd.args(["--model", m]);
    }
    cmd.arg(fastq.to_str().unwrap_or(""))
        .arg(consensus_out.join("consensus.fastq").to_str().unwrap_or(""));

    let result = match cmd.output() {
        Ok(out) if out.status.success() => {
            // Read back consensus sequence (first record in output fastq)
            let seq =
                read_first_fastq_seq(&consensus_out.join("consensus.fastq")).unwrap_or_default();
            ConsensusResult {
                family_id: family.id,
                sequence: if seq.is_empty() { None } else { Some(seq) },
                quality: None,
                method: format!(
                    "medaka smolecule (model={}, depth={})",
                    model.unwrap_or("default"),
                    min_depth
                ),
                success: true,
                error: None,
                family_fastq: Some(fastq.to_path_buf()),
            }
        }
        Ok(out) => {
            let stderr = String::from_utf8_lossy(&out.stderr).to_string();
            ConsensusResult {
                family_id: family.id,
                sequence: None,
                quality: None,
                method: "medaka".to_string(),
                success: false,
                error: Some(format!("medaka exited with error: {}", stderr)),
                family_fastq: Some(fastq.to_path_buf()),
            }
        }
        Err(e) => ConsensusResult {
            family_id: family.id,
            sequence: None,
            quality: None,
            method: "medaka".to_string(),
            success: false,
            error: Some(format!("Failed to launch medaka: {}", e)),
            family_fastq: Some(fastq.to_path_buf()),
        },
    };

    let _ = std::fs::remove_dir_all(&consensus_out);
    result
}

fn read_first_fastq_seq(path: &Path) -> Option<Vec<u8>> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(path).ok()?;
    let mut lines = BufReader::new(f).lines();
    lines.next(); // @header
    let seq_line = lines.next()?.ok()?;
    Some(seq_line.into_bytes())
}

fn phred_from_support(best_count: usize, total_count: usize) -> u8 {
    if total_count == 0 {
        return 0;
    }
    if best_count >= total_count {
        return 60;
    }
    let support = best_count as f64 / total_count as f64;
    let error_prob = (1.0 - support).max(1e-6);
    let phred = (-10.0 * error_prob.log10()).round();
    phred.clamp(0.0, 60.0) as u8
}

// ---------------------------------------------------------------------------
// Simple unbounded edit distance (used by medoid for >64 char seqs)
// ---------------------------------------------------------------------------

pub(crate) fn simple_edit_dist_capped(a: &[u8], b: &[u8], max: usize) -> usize {
    let n = a.len();
    let m = b.len();
    if n.abs_diff(m) > max {
        return usize::MAX;
    }
    let mut prev = (0..=m).collect::<Vec<_>>();
    let mut curr = vec![0usize; m + 1];
    for i in 1..=n {
        curr[0] = i;
        let mut row_min = curr[0];
        for j in 1..=m {
            let cost = if a[i - 1].to_ascii_uppercase() == b[j - 1].to_ascii_uppercase() {
                0
            } else {
                1
            };
            curr[j] = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
            row_min = row_min.min(curr[j]);
        }
        if row_min > max {
            return usize::MAX;
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    let d = prev[m];
    if d <= max {
        d
    } else {
        usize::MAX
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn append_consensus_fastq_record_writes_fastq() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let path = tmp.path().join("consensus.fastq");
        let file = std::fs::File::create(&path).expect("create file");
        let mut writer = std::io::BufWriter::new(file);

        append_consensus_fastq_record(&mut writer, 7, "majority_vote", b"ACGT", None)
            .expect("write consensus");
        writer.flush().expect("flush");

        let content = std::fs::read_to_string(&path).expect("read file");
        assert!(content.contains("@family_000007 method=majority_vote length=4"));
        assert!(content.contains("\nACGT\n+\nIIII\n"));
    }

    #[test]
    fn banded_nw_align_preserves_internal_gap() {
        let aligned = banded_nw_align(b"ACGTAC", b"ACTAC", 2);
        assert_eq!(aligned, b"AC-TAC");
    }

    #[test]
    fn majority_vote_emits_quality_scores() {
        use crate::cluster::{FamilyRead, UmiFamily};

        let mut family = UmiFamily::new(1, "ACGT".to_string());
        family.push(FamilyRead {
            read_id: "r1".to_string(),
            seq: b"ACGT".to_vec(),
            qual: Some(b"IIII".to_vec()),
            strand: None,
        });
        family.push(FamilyRead {
            read_id: "r2".to_string(),
            seq: b"ACGT".to_vec(),
            qual: Some(b"JJJJ".to_vec()),
            strand: None,
        });
        family.push(FamilyRead {
            read_id: "r3".to_string(),
            seq: b"ACGT".to_vec(),
            qual: Some(b"KKKK".to_vec()),
            strand: None,
        });

        let result = majority_vote_consensus(&family, 8);
        assert_eq!(result.sequence.as_deref(), Some(b"ACGT".as_ref()));
        let qual = result.quality.as_deref().expect("quality");
        assert_eq!(qual.len(), 4);
        assert!(qual.iter().all(|&q| q >= 20));
    }
}
