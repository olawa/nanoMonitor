use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

const RC_TABLE: [u8; 256] = {
    let mut table = [b'N'; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = match i as u8 {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        };
        i += 1;
    }
    table
};

#[derive(Debug, Clone)]
pub struct BarcodePair {
    pub sample: String,
    pub bc1: Vec<u8>,
    pub bc2: Vec<u8>,
    pub bc1_rc: Vec<u8>,
    pub bc2_rc: Vec<u8>,
}

impl BarcodePair {
    pub fn new(sample: String, bc1: String, bc2: String) -> Self {
        let bc1_bytes = bc1.into_bytes();
        let bc2_bytes = bc2.into_bytes();
        let bc1_rc = reverse_complement(&bc1_bytes);
        let bc2_rc = reverse_complement(&bc2_bytes);
        Self {
            sample,
            bc1: bc1_bytes,
            bc2: bc2_bytes,
            bc1_rc,
            bc2_rc,
        }
    }
}

pub fn parse_barcodes(path: &Path) -> Result<Vec<BarcodePair>> {
    let file = File::open(path).with_context(|| format!("Failed to open {:?}", path))?;
    let mut lines = BufReader::new(file).lines();
    let mut barcodes = Vec::new();

    while let Some(line) = lines.next() {
        let sample = line?.trim().to_string();
        if sample.is_empty() {
            continue;
        }
        let bc1 = lines.next().context("Missing BC1")??.trim().to_uppercase();
        let bc2 = lines.next().context("Missing BC2")??.trim().to_uppercase();
        barcodes.push(BarcodePair::new(sample, bc1, bc2));
    }
    Ok(barcodes)
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| RC_TABLE[b as usize]).collect()
}

pub fn edit_distance(a: &[u8], b: &[u8], max_dist: usize) -> usize {
    let n = a.len();
    let m = b.len();

    if n.abs_diff(m) > max_dist {
        return usize::MAX;
    }
    if n == 0 {
        return m;
    }
    if m == 0 {
        return n;
    }
    if m > 64 {
        return edit_distance_generic(a, b, max_dist);
    }

    let mut prev_row = [0u8; 65];
    let mut curr_row = [0u8; 65];

    for j in 0..=m {
        prev_row[j] = j as u8;
    }

    for i in 1..=n {
        curr_row[0] = i as u8;
        let mut min_in_row = i as u8;

        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            let val = (prev_row[j] + 1)
                .min(curr_row[j - 1] + 1)
                .min(prev_row[j - 1] + cost);

            curr_row[j] = val;
            if val < min_in_row {
                min_in_row = val;
            }
        }

        if min_in_row as usize > max_dist {
            return usize::MAX;
        }

        prev_row[..=m].copy_from_slice(&curr_row[..=m]);
    }

    let final_dist = prev_row[m] as usize;
    if final_dist <= max_dist {
        final_dist
    } else {
        usize::MAX
    }
}

fn edit_distance_generic(a: &[u8], b: &[u8], max_dist: usize) -> usize {
    let n = a.len();
    let m = b.len();

    if n.abs_diff(m) > max_dist {
        return usize::MAX;
    }
    if n == 0 {
        return m;
    }
    if m == 0 {
        return n;
    }

    let mut prev = (0..=m).collect::<Vec<_>>();
    let mut curr = vec![0usize; m + 1];

    for i in 1..=n {
        curr[0] = i;
        let mut min_in_row = i;
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            let val = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
            curr[j] = val;
            min_in_row = min_in_row.min(val);
        }
        if min_in_row > max_dist {
            return usize::MAX;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let final_dist = prev[m];
    if final_dist <= max_dist {
        final_dist
    } else {
        usize::MAX
    }
}

pub fn find_barcode_fuzzy(seq: &[u8], barcode: &[u8], max_dist: usize) -> bool {
    let blen = barcode.len();
    let slen = seq.len();

    if slen < blen.saturating_sub(max_dist) {
        return false;
    }

    for i in 0..=(slen.saturating_sub(blen)) {
        let sub = &seq[i..i + blen];
        if edit_distance(sub, barcode, max_dist) <= max_dist {
            return true;
        }
    }
    false
}

pub fn match_regions(start: &[u8], end: &[u8], pair: &BarcodePair, max: usize) -> bool {
    if find_barcode_fuzzy(start, &pair.bc1, max) {
        if find_barcode_fuzzy(end, &pair.bc2_rc, max) {
            return true;
        }
        if find_barcode_fuzzy(end, &pair.bc1, max) {
            return true;
        }
    }
    if find_barcode_fuzzy(start, &pair.bc2, max) {
        if find_barcode_fuzzy(end, &pair.bc1_rc, max) {
            return true;
        }
        if find_barcode_fuzzy(end, &pair.bc1, max) {
            return true;
        }
    }
    if find_barcode_fuzzy(start, &pair.bc1_rc, max) {
        if find_barcode_fuzzy(end, &pair.bc2_rc, max) {
            return true;
        }
    }
    false
}

pub fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    if a.len() != b.len() {
        return usize::MAX;
    }
    a.iter().zip(b.iter()).filter(|&(x, y)| x != y).count()
}

pub fn generate_variants(seq: &[u8], max_mismatches: usize) -> Vec<Vec<u8>> {
    let mut variants = Vec::new();
    let mut current = seq.to_vec();

    fn recurse(
        pos: usize,
        mismatches: usize,
        max: usize,
        current: &mut [u8],
        variants: &mut Vec<Vec<u8>>,
    ) {
        if pos == current.len() {
            variants.push(current.to_vec());
            return;
        }

        let original = current[pos];
        recurse(pos + 1, mismatches, max, current, variants);

        if mismatches < max {
            for &base in &[b'A', b'C', b'G', b'T'] {
                if base != original {
                    current[pos] = base;
                    recurse(pos + 1, mismatches + 1, max, current, variants);
                }
            }
            current[pos] = original;
        }
    }

    recurse(0, 0, max_mismatches, &mut current, &mut variants);
    variants
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarcodeId {
    BC1,
    BC2,
    BC1RC,
    BC2RC,
}

#[derive(Debug, Clone)]
pub struct MatchEntry {
    pub sample: String,
    pub barcode_id: BarcodeId,
    pub original: Vec<u8>,
}

pub struct BarcodeMatcher {
    pub map: HashMap<Vec<u8>, Vec<MatchEntry>>,
    pub max_mismatches: usize,
    pub barcodes: Vec<BarcodePair>,
}

impl BarcodeMatcher {
    pub fn new(barcodes: Vec<BarcodePair>, max_mismatches: usize) -> Self {
        let mut map: HashMap<Vec<u8>, Vec<MatchEntry>> = HashMap::new();

        for pair in &barcodes {
            let sequences = [
                (&pair.bc1, BarcodeId::BC1),
                (&pair.bc2, BarcodeId::BC2),
                (&pair.bc1_rc, BarcodeId::BC1RC),
                (&pair.bc2_rc, BarcodeId::BC2RC),
            ];

            for (seq, id) in sequences {
                let variants = generate_variants(seq, max_mismatches);
                for var in variants {
                    map.entry(var).or_default().push(MatchEntry {
                        sample: pair.sample.clone(),
                        barcode_id: id,
                        original: seq.to_vec(),
                    });
                }
            }
        }

        Self {
            map,
            max_mismatches,
            barcodes,
        }
    }

    pub fn match_sample(&self, start_region: &[u8], end_region: &[u8]) -> Option<String> {
        if start_region.len() < 10 || end_region.len() < 10 {
            return None;
        }

        let mut best_sample = None;
        let mut min_total_dist = usize::MAX;

        for i in 0..=(start_region.len().saturating_sub(10)) {
            let s_cand = &start_region[i..i + 10];
            if let Some(s_entries) = self.map.get(s_cand) {
                for j in 0..=(end_region.len().saturating_sub(10)) {
                    let e_cand = &end_region[j..j + 10];
                    if let Some(e_entries) = self.map.get(e_cand) {
                        for s_entry in s_entries {
                            for e_entry in e_entries {
                                if s_entry.sample == e_entry.sample
                                    && is_valid_pair(s_entry.barcode_id, e_entry.barcode_id)
                                {
                                    let total = hamming_distance(s_cand, &s_entry.original)
                                        + hamming_distance(e_cand, &e_entry.original);
                                    if total < min_total_dist {
                                        min_total_dist = total;
                                        best_sample = Some(s_entry.sample.clone());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        best_sample
    }
}

pub fn match_fast_with_anchors(
    seq: &[u8],
    suffix_start: usize,
    matcher: &BarcodeMatcher,
) -> Option<String> {
    const BC2_FWD: &[u8] = b"ATCTACAC";
    const BC1_FWD: &[u8] = b"ATCTCGTA";
    const BC2_REV: &[u8] = b"GTGTAGAT";
    const BC1_REV: &[u8] = b"TACGAGAT";

    let prefix = &seq[..suffix_start.min(seq.len())];
    let suffix = &seq[suffix_start.min(seq.len())..];

    fn get_pocket<'a>(target_seq: &'a [u8], anchor: &[u8], is_after: bool) -> Option<&'a [u8]> {
        target_seq
            .windows(anchor.len())
            .position(|w| w == anchor)
            .and_then(|pos| {
                let (start, end) = if is_after {
                    (
                        pos + anchor.len().saturating_sub(2),
                        pos + anchor.len() + 12,
                    )
                } else {
                    (pos.saturating_sub(12), pos + 2)
                };
                if start < end && end <= target_seq.len() {
                    Some(&target_seq[start..end])
                } else {
                    None
                }
            })
    }

    let bc1_raw_pocket = get_pocket(prefix, BC1_FWD, false)
        .or_else(|| get_pocket(prefix, BC1_REV, true))
        .or_else(|| get_pocket(suffix, BC1_FWD, false))
        .or_else(|| get_pocket(suffix, BC1_REV, true));

    let bc2_raw_pocket = get_pocket(prefix, BC2_FWD, true)
        .or_else(|| get_pocket(prefix, BC2_REV, false))
        .or_else(|| get_pocket(suffix, BC2_FWD, true))
        .or_else(|| get_pocket(suffix, BC2_REV, false));

    if let (Some(b1), Some(b2)) = (bc1_raw_pocket, bc2_raw_pocket) {
        return matcher.match_sample(b1, b2);
    }
    None
}

fn is_valid_pair(a: BarcodeId, b: BarcodeId) -> bool {
    use BarcodeId::*;
    matches!(
        (a, b),
        (BC1, BC2RC)
            | (BC2RC, BC1)
            | (BC1, BC2)
            | (BC2, BC1)
            | (BC2, BC1RC)
            | (BC1RC, BC2)
            | (BC1RC, BC2RC)
            | (BC2RC, BC1RC)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_variants() {
        let variants = generate_variants(b"ATCG", 1);
        assert_eq!(variants.len(), 13);
        assert!(variants.contains(&b"ATCG".to_vec()));
        assert!(variants.contains(&b"GTCG".to_vec()));
        assert!(variants.contains(&b"ATCA".to_vec()));
    }

    #[test]
    fn test_barcode_matcher() {
        let barcodes = vec![
            BarcodePair::new(
                "Sample1".to_string(),
                "ATGCATGCAT".to_string(),
                "GCGCTATAAG".to_string(),
            ),
            BarcodePair::new(
                "Sample2".to_string(),
                "CCGGTTAAAA".to_string(),
                "TTAACCGGGG".to_string(),
            ),
        ];
        let matcher = BarcodeMatcher::new(barcodes, 1);
        let bc1 = b"ATGCATGCAT";
        let bc2_rc = b"CTTATAGCGC";
        assert_eq!(
            matcher.match_sample(bc1, bc2_rc).as_deref(),
            Some("Sample1")
        );
    }
}
