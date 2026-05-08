pub fn mean_qv_from_phred_scores(phred_scores: &[u8]) -> f64 {
    if phred_scores.is_empty() {
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
    use super::{mean_qv_from_fastq_ascii, mean_qv_from_phred_scores, select_qv};

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
}
