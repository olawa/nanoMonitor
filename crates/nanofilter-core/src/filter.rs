use nanoseq_core::header::NanoporeMetadata;
use nanoseq_core::quality::mean_qv_from_fastq_ascii;

pub fn calculate_phred_avg(qual_string: &[u8]) -> f64 {
    mean_qv_from_fastq_ascii(qual_string)
}

/// The quality a `--qv` threshold is compared against.
///
/// Dorado writes its own mean read quality into the `qs` header tag, and that
/// is the number MinKNOW, the sequencing report and `--min-qscore` all use, so
/// it must win when present: recomputing from the quality string gives a
/// slightly different value and would silently disagree with the basecaller.
/// Reads without the tag (older guppy output, trimmed headers) fall back to
/// the computed mean, which is the correct expectation-of-error average and
/// not the arithmetic mean of the Phred scores.
pub fn read_mean_qv(metadata: &NanoporeMetadata, qual_string: &[u8]) -> f64 {
    metadata
        .mean_qscore
        .unwrap_or_else(|| mean_qv_from_fastq_ascii(qual_string))
}

#[cfg(test)]
mod tests {
    use super::{calculate_phred_avg, read_mean_qv};
    use nanoseq_core::header::parse_nanopore_header;

    #[test]
    fn prefers_the_qs_tag_over_recomputing() {
        let meta = parse_nanopore_header(b"@read1 qs:f:18.9348 ch:i:566");
        // All-Q40 quality string, so a recomputed value would be 40.
        let qv = read_mean_qv(&meta, b"IIIIIIIIII");
        assert!((qv - 18.9348).abs() < 1e-9, "got {qv}");
    }

    #[test]
    fn falls_back_to_the_quality_string_without_a_tag() {
        let meta = parse_nanopore_header(b"@read1 runid=abc");
        let qual = b"5555";
        assert!((read_mean_qv(&meta, qual) - calculate_phred_avg(qual)).abs() < 1e-9);
    }
}
