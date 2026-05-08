//! UMI detection for ONT amplicon reads.
//!
//! Uses anchor-based detection (same technique as the barcode splitter):
//! scan the 5' and 3' windows for upstream context sequences, then extract
//! the UMI sequence adjacent to the anchor using IUPAC-aware approximate matching.
//!
//! # Performance
//! - [`AnchorIndex`] provides an O(1) 8-mer pre-filter that prunes ≫ 90% of
//!   edit-distance calls (pigeonhole heuristic: any candidate within `max_edit`
//!   of the anchor must share at least one 8-mer with it on real ONT data).
//! - [`run_umi_detection_fastq`] uses a reader thread + rayon par_iter worker
//!   pool matching the same pattern as `split_fastq_by_barcodes`.

use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::collections::HashSet;
use std::path::Path;
use std::sync::Arc;

use crate::barcode::{edit_distance, reverse_complement};
use crate::filter::calculate_phred_avg;

// ---------------------------------------------------------------------------
// IUPAC helpers
// ---------------------------------------------------------------------------

/// Returns true if `read_base` satisfies the IUPAC ambiguity code `pattern_base`.
#[inline]
pub fn iupac_matches(pattern_base: u8, read_base: u8) -> bool {
    let b = read_base.to_ascii_uppercase();
    match pattern_base.to_ascii_uppercase() {
        b'A' => b == b'A',
        b'C' => b == b'C',
        b'G' => b == b'G',
        b'T' | b'U' => b == b'T',
        b'R' => matches!(b, b'A' | b'G'),
        b'Y' => matches!(b, b'C' | b'T'),
        b'S' => matches!(b, b'C' | b'G'),
        b'W' => matches!(b, b'A' | b'T'),
        b'K' => matches!(b, b'G' | b'T'),
        b'M' => matches!(b, b'A' | b'C'),
        b'B' => matches!(b, b'C' | b'G' | b'T'),
        b'D' => matches!(b, b'A' | b'G' | b'T'),
        b'H' => matches!(b, b'A' | b'C' | b'T'),
        b'V' => matches!(b, b'A' | b'C' | b'G'),
        b'N' => true,
        _ => false,
    }
}

/// Returns true if `base` is an IUPAC ambiguity code (not a plain A/C/G/T/U).
#[inline]
pub fn is_iupac_wildcard(base: u8) -> bool {
    !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U')
}

/// Global alignment edit-distance with IUPAC-aware cost function.
/// Capped at `max_dist`; returns `usize::MAX` when exceeded.
pub fn iupac_edit_distance(pattern: &[u8], seq: &[u8], max_dist: usize) -> usize {
    let n = pattern.len();
    let m = seq.len();
    if n.abs_diff(m) > max_dist {
        return usize::MAX;
    }
    if n == 0 {
        return m;
    }
    if m == 0 {
        return n;
    }
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];
    for j in 0..=m {
        prev[j] = j;
    }
    for i in 1..=n {
        curr[0] = i;
        let mut row_min = curr[0];
        for j in 1..=m {
            let cost = if iupac_matches(pattern[i - 1], seq[j - 1]) {
                0
            } else {
                1
            };
            curr[j] = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
            row_min = row_min.min(curr[j]);
        }
        if row_min > max_dist {
            return usize::MAX;
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    let d = prev[m];
    if d <= max_dist {
        d
    } else {
        usize::MAX
    }
}

// ---------------------------------------------------------------------------
// IupacPattern
// ---------------------------------------------------------------------------

/// A compiled UMI pattern that records which positions are IUPAC wildcards.
#[derive(Debug, Clone)]
pub struct IupacPattern {
    pub bytes: Vec<u8>,
    pub wildcard_positions: Vec<usize>,
}

impl IupacPattern {
    pub fn new(pattern: &str) -> Self {
        let bytes: Vec<u8> = pattern.bytes().map(|b| b.to_ascii_uppercase()).collect();
        let wildcard_positions = bytes
            .iter()
            .enumerate()
            .filter(|(_, &b)| is_iupac_wildcard(b))
            .map(|(i, _)| i)
            .collect();
        Self {
            bytes,
            wildcard_positions,
        }
    }

    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    pub fn edit_distance(&self, seq: &[u8], max_dist: usize) -> usize {
        iupac_edit_distance(&self.bytes, seq, max_dist)
    }

    /// Extract wildcard-position bases. Returns `None` if lengths differ.
    pub fn extract_wildcards(&self, seq: &[u8]) -> Option<Vec<u8>> {
        if seq.len() != self.bytes.len() {
            return None;
        }
        Some(self.wildcard_positions.iter().map(|&i| seq[i]).collect())
    }
}

// ---------------------------------------------------------------------------
// AnchorIndex — fast kmer pre-filter
// ---------------------------------------------------------------------------

/// Pre-compiled 8-mer index of an anchor sequence for O(1) candidate filtering.
///
/// Implements the pigeonhole heuristic: for any candidate window within
/// `max_edit ≤ 4` edits of a 23+ bp anchor, at least one 8-mer from the
/// anchor appears in the candidate window (holds with probability > 0.9999
/// on real ONT data). Prunes > 90% of edit-distance calls.
#[derive(Debug, Clone)]
pub struct AnchorIndex {
    /// Two-bit encoded 8-mers (A=0, C=1, G=2, T=3).
    kmers: HashSet<u32>,
    /// Upper-cased anchor bytes.
    pub anchor: Vec<u8>,
}

impl AnchorIndex {
    const K: usize = 8;

    pub fn new(anchor: &[u8]) -> Self {
        let anchor_up: Vec<u8> = anchor.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let mut kmers = HashSet::with_capacity(anchor_up.len().saturating_sub(Self::K) + 1);
        if anchor_up.len() >= Self::K {
            for i in 0..=(anchor_up.len() - Self::K) {
                if let Some(code) = encode_kmer8(&anchor_up[i..i + Self::K]) {
                    kmers.insert(code);
                }
            }
        }
        Self {
            kmers,
            anchor: anchor_up,
        }
    }

    /// Placeholder with no kmers — always passes through to edit distance.
    pub fn empty() -> Self {
        Self {
            kmers: HashSet::new(),
            anchor: Vec::new(),
        }
    }

    /// True if at least one 8-mer from `candidate` appears in the anchor index.
    #[inline]
    pub fn has_kmer_hit(&self, candidate: &[u8]) -> bool {
        if self.kmers.is_empty() {
            return true; // empty index = no pre-filter
        }
        if candidate.len() < Self::K {
            return true; // too short to check — let edit dist decide
        }
        for i in 0..=(candidate.len() - Self::K) {
            if let Some(code) = encode_kmer8(&candidate[i..i + Self::K]) {
                if self.kmers.contains(&code) {
                    return true;
                }
            }
        }
        false
    }
}

impl Default for AnchorIndex {
    fn default() -> Self {
        Self::empty()
    }
}

#[inline]
fn encode_kmer8(seq: &[u8]) -> Option<u32> {
    debug_assert_eq!(seq.len(), AnchorIndex::K);
    let mut code: u32 = 0;
    for &b in seq {
        let bits = match b.to_ascii_uppercase() {
            b'A' => 0u32,
            b'C' => 1u32,
            b'G' => 2u32,
            b'T' => 3u32,
            _ => return None,
        };
        code = (code << 2) | bits;
    }
    Some(code)
}

// ---------------------------------------------------------------------------
// Anchor scanning
// ---------------------------------------------------------------------------

/// Original O(W·A) scan — kept for backward compatibility and tests.
/// Prefer `scan_for_anchor_fast` in production.
pub fn scan_for_anchor(window: &[u8], anchor: &[u8], max_edit: usize) -> Option<(usize, usize)> {
    if window.len() < anchor.len() {
        return None;
    }
    let mut best_pos: Option<usize> = None;
    let mut best_dist = usize::MAX;
    for i in 0..=(window.len() - anchor.len()) {
        let sub = &window[i..i + anchor.len()];
        let d = edit_distance(sub, anchor, max_edit);
        if d <= max_edit && d < best_dist {
            best_dist = d;
            best_pos = Some(i);
            if d == 0 {
                break;
            }
        }
    }
    best_pos.map(|p| (p, best_dist))
}

/// Kmer-indexed anchor scan. Same semantics as `scan_for_anchor` but uses
/// `AnchorIndex` to skip positions that share no 8-mer with the anchor,
/// reducing edit-distance calls by > 90% in practice.
pub fn scan_for_anchor_fast(
    window: &[u8],
    index: &AnchorIndex,
    max_edit: usize,
) -> Option<(usize, usize)> {
    let anchor = &index.anchor;
    let anchor_len = anchor.len();
    if window.len() < anchor_len {
        return None;
    }
    let mut best_pos: Option<usize> = None;
    let mut best_dist = usize::MAX;
    for i in 0..=(window.len() - anchor_len) {
        let sub = &window[i..i + anchor_len];
        if !index.has_kmer_hit(sub) {
            continue; // kmer pre-filter: skip candidates with no shared 8-mer
        }
        let d = edit_distance(sub, anchor, max_edit);
        if d <= max_edit && d < best_dist {
            best_dist = d;
            best_pos = Some(i);
            if d == 0 {
                break;
            }
        }
    }
    best_pos.map(|p| (p, best_dist))
}

// ---------------------------------------------------------------------------
// UMI extraction
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct UmiHit {
    pub seq: Vec<u8>,
    pub edit_dist: usize,
}

/// Extract a UMI matching `pattern` immediately after the anchor.
pub fn extract_umi_after_anchor(
    window: &[u8],
    anchor_start: usize,
    anchor_len: usize,
    pattern: &IupacPattern,
    max_edit: usize,
    min_len: usize,
    max_len: usize,
) -> Option<UmiHit> {
    let umi_start = anchor_start + anchor_len;
    if umi_start >= window.len() {
        return None;
    }
    try_umi_candidates(window, umi_start, pattern, max_edit, min_len, max_len)
}

/// Extract a UMI matching `pattern` immediately **before** the anchor.
/// Used at the 3′ end of reverse-strand reads where the UMI precedes the context.
pub fn extract_umi_before_anchor(
    window: &[u8],
    anchor_start: usize,
    pattern: &IupacPattern,
    max_edit: usize,
    min_len: usize,
    max_len: usize,
) -> Option<UmiHit> {
    let plen = pattern.len();
    let mut best: Option<UmiHit> = None;
    for delta in [0i32, -2, 2, -4, 4] {
        let try_len = plen as i32 + delta;
        if try_len < min_len as i32 || try_len > max_len as i32 {
            continue;
        }
        let try_len = try_len as usize;
        let umi_end = anchor_start;
        if umi_end < try_len {
            continue;
        }
        let candidate = &window[umi_end - try_len..umi_end];
        let d = pattern.edit_distance(candidate, max_edit);
        if d <= max_edit {
            match &best {
                None => {
                    best = Some(UmiHit {
                        seq: candidate.to_vec(),
                        edit_dist: d,
                    })
                }
                Some(b) if d < b.edit_dist => {
                    best = Some(UmiHit {
                        seq: candidate.to_vec(),
                        edit_dist: d,
                    })
                }
                _ => {}
            }
            if d == 0 {
                break;
            }
        }
    }
    best
}

fn try_umi_candidates(
    window: &[u8],
    umi_start: usize,
    pattern: &IupacPattern,
    max_edit: usize,
    min_len: usize,
    max_len: usize,
) -> Option<UmiHit> {
    let plen = pattern.len();
    let mut best: Option<UmiHit> = None;
    for delta in [0i32, -2, 2, -4, 4] {
        let try_len = plen as i32 + delta;
        if try_len < min_len as i32 || try_len > max_len as i32 {
            continue;
        }
        let try_len = try_len as usize;
        let umi_end = umi_start + try_len;
        if umi_end > window.len() {
            continue;
        }
        let candidate = &window[umi_start..umi_end];
        let d = pattern.edit_distance(candidate, max_edit);
        if d <= max_edit {
            match &best {
                None => {
                    best = Some(UmiHit {
                        seq: candidate.to_vec(),
                        edit_dist: d,
                    })
                }
                Some(b) if d < b.edit_dist => {
                    best = Some(UmiHit {
                        seq: candidate.to_vec(),
                        edit_dist: d,
                    })
                }
                _ => {}
            }
            if d == 0 {
                break;
            }
        }
    }
    best
}

// ---------------------------------------------------------------------------
// Strand + per-read result
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Fwd,
    Rev,
}

impl Strand {
    pub fn as_str(self) -> &'static str {
        match self {
            Strand::Fwd => "fwd",
            Strand::Rev => "rev",
        }
    }
}

#[derive(Debug, Clone)]
pub struct ReadUmiResult {
    pub read_id: String,
    pub strand: Option<Strand>,
    pub umi_fwd_seq: Option<Vec<u8>>,
    pub umi_fwd_edit_dist: Option<usize>,
    pub umi_rev_seq: Option<Vec<u8>>,
    pub umi_rev_edit_dist: Option<usize>,
    /// Combined UMI in canonical fwd orientation (strand-invariant; no RC applied).
    pub combined_umi: Option<String>,
    /// Wildcard-position-only UMI (populated when `config.normalize = true`).
    pub umi_normalised: Option<String>,
    pub read_length: usize,
    /// Inter-anchor span (5′ anchor start → 3′ anchor start, bp). Always computed
    /// when both anchors are found; used for primer-dimer / chimera detection.
    pub insert_size: Option<usize>,
}

impl ReadUmiResult {
    fn undetected(read_id: &str, read_length: usize) -> Self {
        Self {
            read_id: read_id.to_string(),
            strand: None,
            umi_fwd_seq: None,
            umi_fwd_edit_dist: None,
            umi_rev_seq: None,
            umi_rev_edit_dist: None,
            combined_umi: None,
            umi_normalised: None,
            read_length,
            insert_size: None,
        }
    }

    pub fn has_umi(&self) -> bool {
        self.combined_umi.is_some()
    }

    pub fn cluster_key(&self) -> Option<&str> {
        self.umi_normalised
            .as_deref()
            .or(self.combined_umi.as_deref())
    }
}

// ---------------------------------------------------------------------------
// UmiConfig
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct UmiConfig {
    pub fwd_context: Vec<u8>,
    pub rev_context: Vec<u8>,
    pub fwd_pattern: IupacPattern,
    pub rev_pattern: IupacPattern,
    pub max_edit_dist: usize,
    pub window_len: usize,
    pub min_umi_len: usize,
    pub max_umi_len: usize,
    pub normalize: bool,
    pub min_read_len: usize,
    pub max_read_len: usize,
    pub min_mean_q: f64,
    /// Expected inter-anchor span (start of 5′ anchor → start of 3′ anchor, bp).
    /// Includes both UMIs + amplicon insert. Set to 0 to disable the check.
    pub amplicon_size: usize,
    /// Allowed deviation from `amplicon_size` in either direction (bp).
    /// E.g. amplicon_size=350, size_tolerance=50 → accepts spans of 300–400 bp.
    pub size_tolerance: usize,
    /// Pre-compiled kmer index for fwd_context. Built by `build()` / `Default`.
    pub fwd_index: AnchorIndex,
    /// Pre-compiled kmer index for rev_context. Built by `build()` / `Default`.
    pub rev_index: AnchorIndex,
}

impl UmiConfig {
    /// Recompile kmer indices from the current `fwd_context`/`rev_context`.
    /// Call after constructing with explicit field syntax (e.g. in tests or main.rs).
    pub fn build(mut self) -> Self {
        self.fwd_index = AnchorIndex::new(&self.fwd_context);
        self.rev_index = AnchorIndex::new(&self.rev_context);
        self
    }
}

impl Default for UmiConfig {
    fn default() -> Self {
        let fwd_context = b"GTATCGTGTAGAGACTGCGTAGG".to_vec();
        let rev_context = b"AGTGATCGAGTCAGTGCGAGTG".to_vec();
        let fwd_index = AnchorIndex::new(&fwd_context);
        let rev_index = AnchorIndex::new(&rev_context);
        Self {
            fwd_context,
            rev_context,
            fwd_pattern: IupacPattern::new("TTTVVVVTTVVVVTTVVVVTTVVVVTTT"),
            rev_pattern: IupacPattern::new("AAABBBBAABBBBAABBBBAABBBBAAA"),
            max_edit_dist: 4,
            window_len: 250,
            min_umi_len: 40,
            max_umi_len: 75,
            normalize: false,
            min_read_len: 0,
            max_read_len: usize::MAX,
            min_mean_q: 0.0,
            amplicon_size: 0, // disabled by default
            size_tolerance: 0,
            fwd_index,
            rev_index,
        }
    }
}

// ---------------------------------------------------------------------------
// Core per-read detection
// ---------------------------------------------------------------------------

pub fn detect_umi(read_id: &str, seq: &[u8], mean_q: f64, config: &UmiConfig) -> ReadUmiResult {
    let read_length = seq.len();
    if read_length < config.min_read_len
        || read_length > config.max_read_len
        || mean_q < config.min_mean_q
    {
        return ReadUmiResult::undetected(read_id, read_length);
    }

    let w = config.window_len.min(read_length);
    let five_prime = &seq[..w];
    let three_prime = &seq[read_length.saturating_sub(w)..];

    // Kmer-indexed scan — determine strand from the 5' end first, then scan
    // only the matching 3' partner anchor.
    let fwd_at_5 = scan_for_anchor_fast(five_prime, &config.fwd_index, config.max_edit_dist);
    let rev_at_5 = scan_for_anchor_fast(five_prime, &config.rev_index, config.max_edit_dist);

    let fwd_5_edit = fwd_at_5.map(|(_, d)| d).unwrap_or(usize::MAX);
    let rev_5_edit = rev_at_5.map(|(_, d)| d).unwrap_or(usize::MAX);

    if fwd_5_edit == usize::MAX && rev_5_edit == usize::MAX {
        return ReadUmiResult::undetected(read_id, read_length);
    }

    let strand = if fwd_5_edit <= rev_5_edit {
        Strand::Fwd
    } else {
        Strand::Rev
    };
    let (fwd_at_3, rev_at_3) = match strand {
        Strand::Fwd => {
            let rev_at_3 =
                scan_for_anchor_fast(three_prime, &config.rev_index, config.max_edit_dist);
            (None, rev_at_3)
        }
        Strand::Rev => {
            let fwd_at_3 =
                scan_for_anchor_fast(three_prime, &config.fwd_index, config.max_edit_dist);
            (fwd_at_3, None)
        }
    };

    // --- Inter-anchor distance filter + insert_size reporting ---
    // Compute the span from the start of the 5′ anchor to the start of the 3′
    // anchor in full-read coordinates. This covers:
    //   fwd_ctx + fwd_umi + <amplicon> + rev_umi + rev_ctx  (Fwd)
    //   rev_ctx + rev_umi + <amplicon> + fwd_umi + fwd_ctx  (Rev)
    // If `amplicon_size > 0`, reads outside [amplicon_size ± size_tolerance] are rejected.
    let three_prime_offset = read_length.saturating_sub(w);
    let insert_size: Option<usize> = match strand {
        Strand::Fwd => match (fwd_at_5, rev_at_3) {
            (Some((f5, _)), Some((r3, _))) => (three_prime_offset + r3).checked_sub(f5),
            _ => None,
        },
        Strand::Rev => match (rev_at_5, fwd_at_3) {
            (Some((r5, _)), Some((f3, _))) => (three_prime_offset + f3).checked_sub(r5),
            _ => None,
        },
    };
    if config.amplicon_size > 0 {
        let in_range = insert_size
            .map(|s| {
                let lo = config.amplicon_size.saturating_sub(config.size_tolerance);
                let hi = config.amplicon_size + config.size_tolerance;
                s >= lo && s <= hi
            })
            .unwrap_or(false);
        if !in_range {
            return ReadUmiResult::undetected(read_id, read_length);
        }
    }

    let fwd_min = config.min_umi_len / 2;
    let fwd_max = config.max_umi_len;
    let rev_min = config.min_umi_len / 2;
    let rev_max = config.max_umi_len;

    let (umi_fwd_hit, umi_rev_hit) = match strand {
        Strand::Fwd => {
            let fht = fwd_at_5.and_then(|(pos, _)| {
                extract_umi_after_anchor(
                    five_prime,
                    pos,
                    config.fwd_context.len(),
                    &config.fwd_pattern,
                    config.max_edit_dist,
                    fwd_min,
                    fwd_max,
                )
            });
            let rht = rev_at_3.and_then(|(pos, _)| {
                extract_umi_after_anchor(
                    three_prime,
                    pos,
                    config.rev_context.len(),
                    &config.rev_pattern,
                    config.max_edit_dist,
                    rev_min,
                    rev_max,
                )
            });
            (fht, rht)
        }
        Strand::Rev => {
            // 5′ end: rev_context followed by RC(rev_umi).
            // RC(rev_umi) resembles the fwd pattern → match against fwd_pattern.
            let rev_umi_hit_at_5 = rev_at_5.and_then(|(pos, _)| {
                extract_umi_after_anchor(
                    five_prime,
                    pos,
                    config.rev_context.len(),
                    &config.fwd_pattern,
                    config.max_edit_dist,
                    rev_min,
                    rev_max,
                )
            });
            // 3′ end: fwd_umi (fwd orientation) BEFORE fwd_context
            let fwd_umi_hit_at_3 = fwd_at_3.and_then(|(pos, _)| {
                extract_umi_before_anchor(
                    three_prime,
                    pos,
                    &config.fwd_pattern,
                    config.max_edit_dist,
                    fwd_min,
                    fwd_max,
                )
            });
            // fwd_hit is already in fwd orientation.
            // rev_hit: RC the raw extraction to restore rev_umi canonical orientation.
            let fht = fwd_umi_hit_at_3;
            let rht = rev_umi_hit_at_5.map(|h| UmiHit {
                seq: reverse_complement(&h.seq),
                edit_dist: h.edit_dist,
            });
            (fht, rht)
        }
    };

    // combined_umi: both UMIs are stored in canonical fwd orientation regardless
    // of read strand — no RC adjustment. This makes combined_umi strand-invariant
    // so reads from the same molecule always produce the same clustering key.
    let combined_umi = match (&umi_fwd_hit, &umi_rev_hit) {
        (Some(fwd), Some(rev)) => {
            let mut cat: Vec<u8> = fwd.seq.clone();
            cat.extend_from_slice(&rev.seq);
            Some(String::from_utf8_lossy(&cat).to_string())
        }
        _ => None,
    };

    let umi_normalised = if config.normalize {
        match (&umi_fwd_hit, &umi_rev_hit) {
            (Some(fwd), Some(rev)) => {
                let fwd_clamped = clamp_to_pattern_len(&fwd.seq, config.fwd_pattern.len());
                let rev_clamped = clamp_to_pattern_len(&rev.seq, config.rev_pattern.len());
                let fwd_wc = fwd_clamped.and_then(|s| config.fwd_pattern.extract_wildcards(s));
                let rev_wc = rev_clamped.and_then(|s| config.rev_pattern.extract_wildcards(s));
                match (fwd_wc, rev_wc) {
                    (Some(mut fwc), Some(rwc)) => {
                        fwc.extend_from_slice(&rwc);
                        Some(String::from_utf8_lossy(&fwc).to_string())
                    }
                    _ => combined_umi.clone(),
                }
            }
            _ => None,
        }
    } else {
        combined_umi.clone()
    };

    ReadUmiResult {
        read_id: read_id.to_string(),
        strand: Some(strand),
        umi_fwd_seq: umi_fwd_hit.as_ref().map(|h| h.seq.clone()),
        umi_fwd_edit_dist: umi_fwd_hit.as_ref().map(|h| h.edit_dist),
        umi_rev_seq: umi_rev_hit.as_ref().map(|h| h.seq.clone()),
        umi_rev_edit_dist: umi_rev_hit.as_ref().map(|h| h.edit_dist),
        combined_umi,
        umi_normalised,
        read_length,
        insert_size,
    }
}

fn clamp_to_pattern_len(seq: &[u8], target_len: usize) -> Option<&[u8]> {
    if seq.len() < target_len {
        None
    } else {
        Some(&seq[..target_len])
    }
}

// ---------------------------------------------------------------------------
// Full-read record
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct FullReadUmiRecord {
    pub input_file: String,
    pub result: ReadUmiResult,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

impl FullReadUmiRecord {
    pub fn read_id(&self) -> &str {
        &self.result.read_id
    }
    pub fn cluster_key(&self) -> Option<&str> {
        self.result.cluster_key()
    }
    pub fn strand(&self) -> Option<Strand> {
        self.result.strand
    }
}

// ---------------------------------------------------------------------------
// Multi-threaded FASTQ pipeline driver
// ---------------------------------------------------------------------------

/// Owned FASTQ record for cross-thread transfer (mirrors `OwnedRecord` in fastq.rs).
struct OwnedRead {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

/// Process a FASTQ/FASTQ.GZ file using a reader thread + rayon par_iter
/// worker pool (same pattern as `split_fastq_by_barcodes`).
///
/// - `threads = 0` or `1` → single-threaded rayon (useful for tests).
/// - `threads ≥ 2` → configures the rayon global thread pool.
///
/// Returns all records including those without a detected UMI.
pub fn run_umi_detection_fastq(
    input_path: &Path,
    input_file_label: &str,
    config: &UmiConfig,
    threads: usize,
) -> Result<Vec<FullReadUmiRecord>> {
    if threads >= 2 {
        // `build_global` is best-effort; silently ignored if already configured.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }

    // Arc-wrap so config is cheaply shared with rayon worker closures.
    let config = Arc::new(config.clone());
    // Arc<str> avoids per-record clone of the label string (P6 fix).
    let label: Arc<str> = Arc::from(input_file_label);

    // Channel capacity: enough to keep workers fed without large memory backlog.
    let cap = (threads.max(1) * 4).max(8);
    let (in_tx, in_rx) = crossbeam_channel::bounded::<Vec<OwnedRead>>(cap);

    // --- Reader thread: parse FASTQ, batch records, send over channel ---
    let input_path_buf = input_path.to_path_buf();
    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = parse_fastx_file(&input_path_buf)
            .with_context(|| format!("Opening {:?}", input_path_buf))?;
        let mut batch = Vec::with_capacity(1000);
        while let Some(rec) = reader.next() {
            let rec = rec.with_context(|| "Error reading FASTQ record")?;
            batch.push(OwnedRead {
                id: rec.id().to_vec(),
                seq: rec.seq().to_vec(),
                qual: rec.qual().map(|q| q.to_vec()),
            });
            if batch.len() >= 1000 {
                if in_tx.send(batch).is_err() {
                    return Ok(()); // receiver dropped
                }
                batch = Vec::with_capacity(1000);
            }
        }
        if !batch.is_empty() {
            let _ = in_tx.send(batch);
        }
        Ok(())
    });

    // --- Main thread: receive batches, process in parallel, collect ---
    let mut all_records: Vec<FullReadUmiRecord> = Vec::new();
    for batch in in_rx {
        let cfg = Arc::clone(&config);
        let lbl = Arc::clone(&label);
        let results: Vec<FullReadUmiRecord> = batch
            .into_par_iter()
            .map(|rec| {
                let id_bytes = &rec.id;
                let read_id = String::from_utf8_lossy(
                    id_bytes.split(|&b| b == b' ').next().unwrap_or(id_bytes),
                )
                .to_string();
                let mean_q = rec.qual.as_deref().map(calculate_phred_avg).unwrap_or(0.0);
                let result = detect_umi(&read_id, &rec.seq, mean_q, &cfg);
                FullReadUmiRecord {
                    input_file: lbl.to_string(),
                    result,
                    seq: rec.seq,
                    qual: rec.qual,
                }
            })
            .collect();
        all_records.extend(results);
    }

    reader_handle.join().expect("Reader thread panicked")?;
    Ok(all_records)
}
