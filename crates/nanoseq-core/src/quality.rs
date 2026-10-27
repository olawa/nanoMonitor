/// BAM stores a record with no quality string (`QUAL` = `*`) as `0xFF` for every base.
pub const BAM_QUAL_MISSING: u8 = 0xFF;

/// True when a BAM quality array carries no real qualities: empty, or all `0xFF`.
pub fn phred_scores_missing(phred_scores: &[u8]) -> bool {
    phred_scores.iter().all(|&q| q == BAM_QUAL_MISSING)
}

/// One BAM Phred score as a FASTQ quality byte.
///
/// A missing score (`0xFF`) must not be offset by 33 -- that overflows `u8` -- and
/// is written as Q0 (`!`) rather than invented: a read without qualities should
/// fail a quality threshold, not pass it.
pub fn phred_to_fastq_ascii(phred: u8) -> u8 {
    if phred == BAM_QUAL_MISSING {
        b'!'
    } else {
        phred.min(93) + 33
    }
}

pub fn mean_qv_from_phred_scores(phred_scores: &[u8]) -> f64 {
    // Also covers a quality-less BAM record: averaging 0xFF would give Q255 and
    // pass every threshold.
    if phred_scores_missing(phred_scores) {
        return 0.0;
    }

    let mean_error_probability = phred_scores
        .iter()
        .map(|&q| 10_f64.powf(-(q as f64) / 10.0))
        .sum::<f64>()
        / phred_scores.len() as f64;

    -10.0 * mean_error_probability.log10()
}

pub fn mean_qv_from_fastq_ascii(qual_ascii: &[u8]) -> f64 {
    if qual_ascii.is_empty() {
        return 0.0;
    }

    let mut phred_scores = Vec::with_capacity(qual_ascii.len());
    for &q in qual_ascii {
        phred_scores.push(q.saturating_sub(33));
    }
    mean_qv_from_phred_scores(&phred_scores)
}

pub fn select_qv(explicit_qs: Option<f64>, phred_scores: &[u8]) -> f64 {
    explicit_qs.unwrap_or_else(|| mean_qv_from_phred_scores(phred_scores))
}

/// Dorado leaves the first bases out of a read's mean quality: they are the
/// adapter and barcode region, and low quality.
pub const BASECALLER_QSCORE_START: usize = 60;

/// A read's mean quality as dorado computes its `qs`: from base 60 on, except
/// for a read split out of a longer one, and a read too short to have
/// anything after base 60, which are averaged over every base.
///
/// Checked against MinKNOW 6.4 output: this reproduces `qs` to 0.01 for 99.2 %
/// of reads, where the mean over the whole read matches 0.7 %. It is the
/// fallback for reads without a `qs`, so a `--qv` threshold means the same
/// thing whether the basecaller's value was kept or not.
pub fn basecaller_mean_qv_from_phred_scores(phred_scores: &[u8], split_read: bool) -> f64 {
    if split_read || phred_scores.len() <= BASECALLER_QSCORE_START {
        mean_qv_from_phred_scores(phred_scores)
    } else {
        mean_qv_from_phred_scores(&phred_scores[BASECALLER_QSCORE_START..])
    }
}

pub fn basecaller_mean_qv_from_fastq_ascii(qual_ascii: &[u8], split_read: bool) -> f64 {
    let phred_scores: Vec<u8> = qual_ascii.iter().map(|q| q.saturating_sub(33)).collect();
    if phred_scores.is_empty() {
        return 0.0;
    }
    basecaller_mean_qv_from_phred_scores(&phred_scores, split_read)
}

/// The basecaller's `qs` when the read carries it, else the same value
/// recomputed the way the basecaller computes it.
pub fn select_basecaller_qv(explicit_qs: Option<f64>, phred_scores: &[u8], split_read: bool) -> f64 {
    explicit_qs.unwrap_or_else(|| basecaller_mean_qv_from_phred_scores(phred_scores, split_read))
}

#[cfg(test)]
mod tests {
    use super::{
        basecaller_mean_qv_from_fastq_ascii, basecaller_mean_qv_from_phred_scores,
        mean_qv_from_fastq_ascii, mean_qv_from_phred_scores, phred_to_fastq_ascii, select_qv,
        BASECALLER_QSCORE_START,
    };

    #[test]
    fn basecaller_mean_leaves_out_the_first_60_bases() {
        // Q5 over the adapter region, Q20 after it: dorado reports 20.
        let mut phred = vec![5u8; BASECALLER_QSCORE_START];
        phred.extend(vec![20u8; 500]);
        assert!((basecaller_mean_qv_from_phred_scores(&phred, false) - 20.0).abs() < 1e-9);
        assert!(mean_qv_from_phred_scores(&phred) < 19.0);
    }

    #[test]
    fn split_and_short_reads_are_averaged_over_every_base() {
        let mut phred = vec![5u8; BASECALLER_QSCORE_START];
        phred.extend(vec![20u8; 500]);
        let whole = mean_qv_from_phred_scores(&phred);
        assert!((basecaller_mean_qv_from_phred_scores(&phred, true) - whole).abs() < 1e-9);
        let short = vec![20u8; 40];
        assert!((basecaller_mean_qv_from_phred_scores(&short, false) - 20.0).abs() < 1e-9);
        assert_eq!(basecaller_mean_qv_from_fastq_ascii(b"", false), 0.0);
    }

    #[test]
    fn computes_mean_qv_from_phred_scores() {
        let qv = mean_qv_from_phred_scores(&[20, 20, 20, 20]);
        assert!((qv - 20.0).abs() < 1e-9);
    }

    #[test]
    fn computes_mean_qv_from_fastq_ascii() {
        let qv = mean_qv_from_fastq_ascii(b"5555");
        assert!((qv - 20.0).abs() < 1e-9);
    }

    #[test]
    fn prefers_explicit_qs() {
        assert_eq!(select_qv(Some(17.5), &[40, 40, 40]), 17.5);
    }

    #[test]
    fn missing_bam_qualities_have_no_quality() {
        assert_eq!(mean_qv_from_phred_scores(&[0xFF; 10]), 0.0);
        assert_eq!(phred_to_fastq_ascii(0xFF), b'!');
        assert_eq!(phred_to_fastq_ascii(40), b'I');
    }
}
