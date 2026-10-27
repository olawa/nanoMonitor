use nanoseq_core::header::NanoporeMetadata;
use nanoseq_core::quality::{basecaller_mean_qv_from_fastq_ascii, mean_qv_from_fastq_ascii};

pub fn calculate_phred_avg(qual_string: &[u8]) -> f64 {
    mean_qv_from_fastq_ascii(qual_string)
}

/// The quality a `--qv` threshold is compared against.
///
/// Dorado writes its own mean read quality into the `qs` header tag, and that
/// is the number MinKNOW, the sequencing report and `--min-qscore` all use, so
/// it must win when present: recomputing from the quality string gives a
/// slightly different value and would silently disagree with the basecaller.
/// Reads without the tag (MinKNOW's FASTQ, older guppy output) fall back to
/// the value recomputed the way dorado computes it -- an expectation-of-error
/// average from base 60 on -- so the threshold means the same either way.
/// Averaging every base instead passed 1 250 of the 4 084 reads MinKNOW's own
/// qs puts at Q20 or above in one MinKNOW FASTQ.
pub fn read_mean_qv(metadata: &NanoporeMetadata, qual_string: &[u8]) -> f64 {
    metadata.mean_qscore.unwrap_or_else(|| {
        basecaller_mean_qv_from_fastq_ascii(qual_string, metadata.is_split_read())
    })
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
        // Short enough that the adapter rule does not apply.
        let meta = parse_nanopore_header(b"@read1 runid=abc");
        let qual = b"5555";
        assert!((read_mean_qv(&meta, qual) - calculate_phred_avg(qual)).abs() < 1e-9);
    }

    #[test]
    fn fallback_leaves_out_the_adapter_region_like_dorado() {
        // Q5 for 60 bases, then Q20: dorado's qs would be 20.
        let mut qual = vec![b'&'; 60];
        qual.extend(vec![b'5'; 300]);
        let meta = parse_nanopore_header(b"@read1 runid=abc parent_read_id=read1");
        assert!((read_mean_qv(&meta, &qual) - 20.0).abs() < 1e-9);
        let split = parse_nanopore_header(b"@child runid=abc parent_read_id=parent");
        assert!(read_mean_qv(&split, &qual) < 19.0);
    }
}
