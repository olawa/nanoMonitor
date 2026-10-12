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

#[cfg(test)]
mod tests {
    use super::{
        mean_qv_from_fastq_ascii, mean_qv_from_phred_scores, phred_to_fastq_ascii, select_qv,
    };

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
