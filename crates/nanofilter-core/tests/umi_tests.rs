//! Integration tests for the UMI pipeline.
//!
//! Run with: `cargo test -p nanofilter-core`

use nanofilter_core::{
    cluster::{
        cluster_by_umi, compute_family_stats, filter_and_downsample, ClusterFilterConfig,
        ClusterMode,
    },
    report::{build_summary, write_summary_tsv},
    umi::{
        detect_umi, is_iupac_wildcard, iupac_matches, scan_for_anchor, IupacPattern, Strand,
        UmiConfig,
    },
};

// ---------------------------------------------------------------------------
// Helper: build a minimal FullReadUmiRecord with a known UMI key
// ---------------------------------------------------------------------------
use nanofilter_core::umi::FullReadUmiRecord;

fn make_record(
    read_id: &str,
    combined_umi: Option<&str>,
    strand: Option<Strand>,
) -> FullReadUmiRecord {
    use nanofilter_core::umi::ReadUmiResult;
    FullReadUmiRecord {
        input_file: "test.fastq".to_string(),
        result: ReadUmiResult {
            read_id: read_id.to_string(),
            strand,
            umi_fwd_seq: combined_umi.map(|u| u.as_bytes()[..u.len() / 2].to_vec()),
            umi_fwd_edit_dist: Some(0),
            umi_rev_seq: combined_umi.map(|u| u.as_bytes()[u.len() / 2..].to_vec()),
            umi_rev_edit_dist: Some(0),
            combined_umi: combined_umi.map(|s| s.to_string()),
            umi_normalised: combined_umi.map(|s| s.to_string()),
            read_length: 3000,
            insert_size: None,
        },
        seq: b"ACGT".repeat(750),
        qual: Some(b"I".repeat(3000).to_vec()),
    }
}

// ---------------------------------------------------------------------------
// 1. IUPAC wildcard detection
// ---------------------------------------------------------------------------
#[test]
fn test_iupac_wildcard_detection() {
    assert!(!is_iupac_wildcard(b'A'));
    assert!(!is_iupac_wildcard(b'C'));
    assert!(!is_iupac_wildcard(b'G'));
    assert!(!is_iupac_wildcard(b'T'));
    assert!(is_iupac_wildcard(b'V')); // A, C, G
    assert!(is_iupac_wildcard(b'B')); // C, G, T
    assert!(is_iupac_wildcard(b'N'));
}

#[test]
fn test_iupac_matches() {
    // V = A, C, G (not T)
    assert!(iupac_matches(b'V', b'A'));
    assert!(iupac_matches(b'V', b'C'));
    assert!(iupac_matches(b'V', b'G'));
    assert!(!iupac_matches(b'V', b'T'));
    // B = C, G, T (not A)
    assert!(iupac_matches(b'B', b'C'));
    assert!(iupac_matches(b'B', b'G'));
    assert!(iupac_matches(b'B', b'T'));
    assert!(!iupac_matches(b'B', b'A'));
    // N = any
    assert!(iupac_matches(b'N', b'A'));
    assert!(iupac_matches(b'N', b'T'));
    // Fixed
    assert!(iupac_matches(b'A', b'A'));
    assert!(!iupac_matches(b'A', b'C'));
}

// ---------------------------------------------------------------------------
// 2. IupacPattern wildcard position extraction
// ---------------------------------------------------------------------------
#[test]
fn test_iupac_wildcard_extraction() {
    let pat = IupacPattern::new("TTTVVVVTTT"); // 10 chars; positions 3,4,5,6 are V
    assert_eq!(pat.wildcard_positions, vec![3, 4, 5, 6]);

    // seq must be exactly 10 chars to match the 10-char pattern
    let seq = b"TTTACGTTT"; // only 9 – should return None
    let wc = pat.extract_wildcards(seq);
    assert!(
        wc.is_none(),
        "extract_wildcards should return None when lengths differ"
    );

    // Correct: exactly 10 chars
    let seq10 = b"TTTACGTTTN"; // 10 chars (trailing N is in position 9 = T position)
                               // Pattern pos 9 = 'T' (fixed); 'N' won't match 'T', but extract_wildcards only
                               // reads at wildcard positions (3-6) regardless – so we still get the V bases
    let wc10 = pat.extract_wildcards(seq10);
    assert!(
        wc10.is_some(),
        "extract_wildcards should succeed for exact-length (10-char) seq"
    );
    assert_eq!(wc10.unwrap(), b"ACGT");
}

// ---------------------------------------------------------------------------
// 3. Anchor scanning
// ---------------------------------------------------------------------------
#[test]
fn test_scan_for_anchor_exact() {
    let anchor = b"GTATCGTGTAG";
    let window = b"NNNNNNNGTATCGTGTAGNNNNNNNN";
    let result = scan_for_anchor(window, anchor, 0);
    assert!(result.is_some());
    let (pos, dist) = result.unwrap();
    assert_eq!(pos, 7);
    assert_eq!(dist, 0);
}

#[test]
fn test_scan_for_anchor_approximate() {
    let anchor = b"GTATCGTGTAG";
    // One mismatch: position 4 changed A→T
    let window = b"NNNNNNNGTATCTTGTAGNNNNNNNN";
    let result = scan_for_anchor(window, anchor, 2);
    assert!(result.is_some());
    let (_, dist) = result.unwrap();
    assert!(dist <= 2, "edit distance should be ≤ 2");
}

#[test]
fn test_scan_for_anchor_not_found() {
    let anchor = b"GTATCGTGTAG";
    let window = b"AAAAAAAAAAAAAAAA";
    let result = scan_for_anchor(window, anchor, 1);
    assert!(result.is_none());
}

// ---------------------------------------------------------------------------
// 4. Per-read UMI extraction on a synthetic forward-strand read
// ---------------------------------------------------------------------------
fn build_synthetic_read_fwd() -> Vec<u8> {
    // Structure: [250 random] [fwd_context] [fwd_umi] [amplicon] [rev_context] [rev_umi] [250 random]
    let prefix = b"A".repeat(50);
    let fwd_ctx = b"GTATCGTGTAGAGACTGCGTAGG";
    // fwd_umi is 29 chars matching TTTVVVVTTVVVVTTVVVVTTVVVVTTT
    let fwd_umi = b"TTTACGCTTACGCTTACGCTTACGCTTT";
    let amplicon = b"G".repeat(2000);
    let rev_ctx = b"AGTGATCGAGTCAGTGCGAGTG";
    // rev_umi is 28 chars matching AAABBBBAABBBBAABBBBAABBBBAAA
    let rev_umi = b"AAACGTTAACGTTAACGTTAACGTTAAA";
    let suffix = b"T".repeat(50);

    let mut read = Vec::new();
    read.extend_from_slice(&prefix);
    read.extend_from_slice(fwd_ctx);
    read.extend_from_slice(fwd_umi);
    read.extend_from_slice(&amplicon);
    read.extend_from_slice(rev_ctx);
    read.extend_from_slice(rev_umi);
    read.extend_from_slice(&suffix);
    read
}

#[test]
fn test_umi_extraction_fwd() {
    let seq = build_synthetic_read_fwd();
    let config = UmiConfig::default();
    let result = detect_umi("read_fwd_001", &seq, 30.0, &config);

    assert!(
        result.combined_umi.is_some(),
        "Forward read should have a combined UMI; got strand={:?}",
        result.strand
    );
    assert_eq!(result.strand, Some(Strand::Fwd));
    assert!(
        result.umi_fwd_seq.is_some(),
        "Forward UMI should be extracted"
    );
    assert!(
        result.umi_rev_seq.is_some(),
        "Reverse UMI should be extracted"
    );

    let fwd_edit = result.umi_fwd_edit_dist.unwrap_or(usize::MAX);
    let rev_edit = result.umi_rev_edit_dist.unwrap_or(usize::MAX);
    assert!(
        fwd_edit <= 4,
        "Forward UMI edit dist should be ≤ 4, got {}",
        fwd_edit
    );
    assert!(
        rev_edit <= 4,
        "Reverse UMI edit dist should be ≤ 4, got {}",
        rev_edit
    );
}

// ---------------------------------------------------------------------------
// 5. Strand normalization: same molecule, both strands → same normalised UMI
// ---------------------------------------------------------------------------

#[test]
fn test_strand_normalization() {
    let fwd_seq = build_synthetic_read_fwd();
    let mut config = UmiConfig::default();
    config.normalize = true;

    let fwd_result = detect_umi("read_fwd", &fwd_seq, 30.0, &config);
    assert!(
        fwd_result.umi_normalised.is_some(),
        "Forward read normalised UMI should be present"
    );

    // Build a reverse-strand read of the SAME physical molecule.
    // On the rev strand, the 5' end has:
    //   RC(rev_umi) after rev_ctx (because reading opposite strand)
    // And the 3' end has:
    //   fwd_umi (same literal) before fwd_ctx
    let rev_ctx = b"AGTGATCGAGTCAGTGCGAGTG";
    // fwd-strand rev_umi =  AAACGTTAACGTTAACGTTAACGTTAAA
    // RC of that          =  TTTAACGTTAACGTTAACGTTAACGTTT
    let rev_umi_rc: Vec<u8> =
        nanofilter_core::barcode::reverse_complement(b"AAACGTTAACGTTAACGTTAACGTTAAA");
    let amplicon_r = b"C".repeat(2000);
    let fwd_umi = b"TTTACGCTTACGCTTACGCTTACGCTTT";
    let fwd_ctx = b"GTATCGTGTAGAGACTGCGTAGG";
    // 150-bp suffix ensures fwd_ctx + fwd_umi fit comfortably inside the 250-bp 3' window
    let suffix = b"T".repeat(150);
    let mut rev_seq = Vec::new();
    rev_seq.extend_from_slice(&b"A".repeat(50));
    rev_seq.extend_from_slice(rev_ctx);
    rev_seq.extend_from_slice(&rev_umi_rc);
    rev_seq.extend_from_slice(&amplicon_r);
    rev_seq.extend_from_slice(fwd_umi);
    rev_seq.extend_from_slice(fwd_ctx);
    rev_seq.extend_from_slice(&suffix);

    let rev_result = detect_umi("read_rev", &rev_seq, 30.0, &config);
    assert!(
        rev_result.umi_normalised.is_some(),
        "Reverse read normalised UMI should be present; strand={:?}",
        rev_result.strand
    );

    let fwd_norm = fwd_result.umi_normalised.as_deref().unwrap_or("");
    let rev_norm = rev_result.umi_normalised.as_deref().unwrap_or("");
    let d = nanofilter_core::umi::iupac_edit_distance(fwd_norm.as_bytes(), rev_norm.as_bytes(), 8);
    assert!(
        d <= 4,
        "Normalised UMIs of both strands should match within edit dist 4; got '{}' vs '{}' (d={})",
        fwd_norm,
        rev_norm,
        d
    );
}

// ---------------------------------------------------------------------------
// 6. Exact clustering
// ---------------------------------------------------------------------------
#[test]
fn test_cluster_exact() {
    let records: Vec<FullReadUmiRecord> = (0..10)
        .map(|i| {
            make_record(
                &format!("r{}", i),
                Some("AAACCCGGGTTTNNN"),
                Some(Strand::Fwd),
            )
        })
        .chain((10..13).map(|i| {
            make_record(
                &format!("r{}", i),
                Some("TTTGGGCCCAAA___"),
                Some(Strand::Rev),
            )
        }))
        .collect();

    let families = cluster_by_umi(&records, &ClusterMode::Exact);
    assert_eq!(
        families.len(),
        2,
        "Exact clustering should produce 2 families"
    );
    let sizes: std::collections::HashSet<usize> =
        families.iter().map(|f| f.total_reads()).collect();
    assert!(sizes.contains(&10) && sizes.contains(&3));
}

// ---------------------------------------------------------------------------
// 7. Approximate clustering
// ---------------------------------------------------------------------------
#[test]
fn test_cluster_approximate() {
    // UMIs that differ by ≤ 3 bases should merge
    let umis = [
        "AAACCCGGGTTT", // seed
        "AAACCCGGGTTG", // 1 edit from seed
        "AAACCCGGGTCC", // 2 edits from seed
        "TTTGGGCCCAAA", // different family
        "TTTGGGCCCAAB", // 1 edit from second seed
    ];
    let records: Vec<FullReadUmiRecord> = umis
        .iter()
        .enumerate()
        .map(|(i, &u)| make_record(&format!("r{}", i), Some(u), Some(Strand::Fwd)))
        .collect();

    let families = cluster_by_umi(&records, &ClusterMode::Approximate { max_edit: 3 });
    assert_eq!(
        families.len(),
        2,
        "Approximate clustering should produce 2 families"
    );
}

// ---------------------------------------------------------------------------
// 8. Below minimum reads → failing family
// ---------------------------------------------------------------------------
#[test]
fn test_below_min_reads() {
    let records: Vec<FullReadUmiRecord> = (0..3)
        .map(|i| make_record(&format!("r{}", i), Some("AAACCCGGGTTT"), Some(Strand::Fwd)))
        .collect();

    let families = cluster_by_umi(&records, &ClusterMode::Exact);
    assert_eq!(families.len(), 1);

    let config = ClusterFilterConfig {
        min_reads: 4,
        max_reads: 80,
        balance_strands: false,
    };
    let filtered = filter_and_downsample(families, &config);

    assert_eq!(
        filtered.passing.len(),
        0,
        "Family with 3 reads should fail min_reads=4"
    );
    assert_eq!(filtered.failing.len(), 1);
}

// ---------------------------------------------------------------------------
// 9. Oversized cluster → downsampled to max_reads
// ---------------------------------------------------------------------------
#[test]
fn test_oversized_cluster() {
    let records: Vec<FullReadUmiRecord> = (0..100)
        .map(|i| make_record(&format!("r{}", i), Some("AAACCCGGGTTT"), Some(Strand::Fwd)))
        .collect();

    let families = cluster_by_umi(&records, &ClusterMode::Exact);
    let config = ClusterFilterConfig {
        min_reads: 4,
        max_reads: 80,
        balance_strands: false,
    };
    let filtered = filter_and_downsample(families, &config);

    assert_eq!(filtered.passing.len(), 1);
    assert_eq!(
        filtered.passing[0].total_reads(),
        80,
        "Oversized family should be downsampled to max_reads=80"
    );
}

// ---------------------------------------------------------------------------
// 10. Strand balance
// ---------------------------------------------------------------------------
#[test]
fn test_strand_balance() {
    // 60 fwd + 40 rev = 100 reads; balance_strands=true, max_reads=80 → 40 + 40 = 80
    let records: Vec<FullReadUmiRecord> = (0..60)
        .map(|i| make_record(&format!("f{}", i), Some("AAACCCGGGTTT"), Some(Strand::Fwd)))
        .chain(
            (0..40)
                .map(|i| make_record(&format!("r{}", i), Some("AAACCCGGGTTT"), Some(Strand::Rev))),
        )
        .collect();

    let families = cluster_by_umi(&records, &ClusterMode::Exact);
    let config = ClusterFilterConfig {
        min_reads: 4,
        max_reads: 80,
        balance_strands: true,
    };
    let filtered = filter_and_downsample(families, &config);

    assert_eq!(filtered.passing.len(), 1);
    let fam = &filtered.passing[0];
    assert!(
        fam.n_fwd <= 40,
        "Balanced fwd should be ≤ 40, got {}",
        fam.n_fwd
    );
    assert!(
        fam.n_rev <= 40,
        "Balanced rev should be ≤ 40, got {}",
        fam.n_rev
    );
    assert_eq!(fam.total_reads(), fam.n_fwd + fam.n_rev);
}

// ---------------------------------------------------------------------------
// 11. Family stats generation round-trip
// ---------------------------------------------------------------------------
#[test]
fn test_family_stats_generation() {
    let records: Vec<FullReadUmiRecord> = (0..10)
        .map(|i| make_record(&format!("r{}", i), Some("AAACCCGGGTTT"), Some(Strand::Fwd)))
        .collect();
    let all_families = cluster_by_umi(&records, &ClusterMode::Exact);
    let config = ClusterFilterConfig {
        min_reads: 4,
        max_reads: 80,
        balance_strands: false,
    };
    let filtered = filter_and_downsample(all_families.clone(), &config);

    let stats = compute_family_stats(&filtered.passing, &filtered.failing, &config);

    assert_eq!(stats.len(), 1);
    let s = &stats[0];
    assert_eq!(s.n_cluster_consensus, 10);
    assert_eq!(s.written, 10);
    assert_eq!(s.cluster_written, 1);
    assert_eq!(s.skipped, 0);
    assert_eq!(s.min_reads_required, 4);
    assert_eq!(s.max_reads_allowed, 80);
}

// ---------------------------------------------------------------------------
// 12. Summary table — multi-sample / multi-amplicon
// ---------------------------------------------------------------------------
#[test]
fn test_summary_table_multi_sample() {
    // SummaryRecord is returned directly from build_summary, no import needed.

    let tmp = tempfile::tempdir().expect("tempdir");
    let summary_path = tmp.path().join("summary.tsv");

    let umi_config = UmiConfig::default();
    let cluster_config = ClusterFilterConfig::default();

    for (sample, amplicon) in &[("S1", "AMP1"), ("S1", "AMP2"), ("S2", "AMP1")] {
        let records: Vec<FullReadUmiRecord> = (0..20)
            .map(|i| make_record(&format!("r{}", i), Some("AAACCCGGGTTT"), Some(Strand::Fwd)))
            .collect();
        let all_fam = cluster_by_umi(&records, &ClusterMode::Exact);
        let filtered = filter_and_downsample(all_fam.clone(), &cluster_config);

        let rec = build_summary(
            sample,
            amplicon,
            &records,
            &filtered.passing,
            &filtered.failing,
            &cluster_config,
            &umi_config,
        );
        write_summary_tsv(&[rec], &summary_path).expect("write summary");
    }

    let content = std::fs::read_to_string(&summary_path).expect("read summary");
    let lines: Vec<&str> = content.lines().collect();
    // 1 header + 3 data rows
    assert_eq!(lines.len(), 4, "Summary should have header + 3 data rows");
    assert!(lines[0].starts_with("sample\tamplicon"));
    assert!(lines[1].contains("S1\tAMP1"));
    assert!(lines[2].contains("S1\tAMP2"));
    assert!(lines[3].contains("S2\tAMP1"));
}
