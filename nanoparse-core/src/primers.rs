//! Primer loading and k-mer indexing

use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct Primer {
    pub name: String,
    pub sequence: String,
    pub sequence_rc: String,
    pub chrom: Option<String>,
    pub start: Option<i64>,
    pub end: Option<i64>,
}

impl Primer {
    pub fn new(name: &str, sequence: &str, region: Option<&str>) -> Self {
        let seq_upper = sequence.to_uppercase();
        let seq_rc = reverse_complement(&seq_upper);

        let mut chrom = None;
        let mut start = None;
        let mut end = None;

        if let Some(r) = region {
            // Parse region: chrom:start-end  or  chrom:start-end(Gene)
            // Example: chr17:43000-43150
            let r_clean = r.split('(').next().unwrap_or(r); // Remove (Gene) suffix
            let parts: Vec<&str> = r_clean.split(':').collect();
            if parts.len() == 2 {
                chrom = Some(parts[0].to_string());
                let coords: Vec<&str> = parts[1].split('-').collect();
                if coords.len() == 2 {
                    if let (Ok(s), Ok(e)) = (coords[0].parse::<i64>(), coords[1].parse::<i64>()) {
                        start = Some(s);
                        end = Some(e);
                    }
                }
            }
        }

        Self {
            name: name.to_string(),
            sequence: seq_upper,
            sequence_rc: seq_rc,
            chrom,
            start,
            end,
        }
    }
}

/// Reverse complement a DNA sequence
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N',
        })
        .collect()
}

/// K-mer index for fast primer lookup
#[derive(Debug)]
pub struct KmerIndex {
    pub kmer_size: usize,
    /// Maps k-mer -> set of primer indices
    pub index: HashMap<String, HashSet<usize>>,
    pub primers: Vec<Primer>,
}

impl KmerIndex {
    /// Build k-mer index from primers
    pub fn build(primers: Vec<Primer>, kmer_size: usize) -> Self {
        let mut index: HashMap<String, HashSet<usize>> = HashMap::new();

        for (idx, primer) in primers.iter().enumerate() {
            // Index both forward and reverse complement
            for seq in [&primer.sequence, &primer.sequence_rc] {
                if seq.len() >= kmer_size {
                    for i in 0..=(seq.len() - kmer_size) {
                        let kmer = &seq[i..i + kmer_size];
                        index.entry(kmer.to_string()).or_default().insert(idx);
                    }
                }
            }
        }

        Self {
            kmer_size,
            index,
            primers,
        }
    }

    /// Find candidate primers for a sequence based on k-mer hits
    /// Returns primer indices sorted by hit count (descending)
    ///
    /// Note: Optimized for small number of primers (<256) and reasonable k-mer size
    pub fn find_candidates(&self, seq: &str, min_hit_ratio: f32) -> Vec<usize> {
        if seq.len() < self.kmer_size {
            return vec![];
        }

        let seq_upper = seq.to_uppercase();
        let num_primers = self.primers.len();

        // Optimize: Use small fixed-size buffer or Vec for counts to avoid HashMap overhead
        // Assuming num_primers is relatively small (e.g. 50-1000)
        let mut hit_counts = vec![0u8; num_primers];
        let mut total_hits = 0;

        for i in 0..=(seq_upper.len() - self.kmer_size) {
            let kmer = &seq_upper[i..i + self.kmer_size];
            if let Some(primer_indices) = self.index.get(kmer) {
                for &idx in primer_indices {
                    if idx < num_primers {
                        hit_counts[idx] = hit_counts[idx].saturating_add(1);
                        total_hits += 1;
                    }
                }
            }
        }

        if total_hits == 0 {
            return vec![];
        }

        let mut candidates = Vec::new();
        for (idx, &hits) in hit_counts.iter().enumerate() {
            if hits > 0 {
                let primer_len = self.primers[idx].sequence.len();
                let max_possible = primer_len.saturating_sub(self.kmer_size) + 1;
                let ratio = hits as f32 / max_possible as f32;

                if ratio >= min_hit_ratio {
                    candidates.push((idx, hits));
                }
            }
        }

        // Sort by hit count descending
        candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        candidates.into_iter().map(|(idx, _)| idx).collect()
    }

    pub fn get_primer(&self, idx: usize) -> Option<&Primer> {
        self.primers.get(idx)
    }
}

/// Load primers from TSV file (name<TAB>sequence[<TAB>Region])
pub fn load_primers(path: &str) -> Result<Vec<Primer>> {
    let file = File::open(path).context("Failed to open primers file")?;
    let reader = BufReader::new(file);
    let mut primers = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Split on any whitespace (tabs or spaces)
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let name = parts[0];
            let sequence = parts[1];
            let region = if parts.len() >= 3 {
                Some(parts[2])
            } else {
                None
            };

            if !name.is_empty()
                && !sequence.is_empty()
                && sequence.chars().all(|c| "ACGTNacgtn".contains(c))
            {
                primers.push(Primer::new(name, sequence, region));
            }
        }
    }

    log::info!("Loaded {} primers from {}", primers.len(), path);
    Ok(primers)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
    }

    #[test]
    fn test_kmer_index() {
        let primers = vec![
            Primer::new("P1", "ATCGATCG", None),
            Primer::new("P2", "GGCCGGCC", None),
        ];
        let index = KmerIndex::build(primers, 4);

        let candidates = index.find_candidates("ATCGATCG", 0.5);
        assert!(!candidates.is_empty());
        assert_eq!(candidates[0], 0); // P1 should be first
    }
}
