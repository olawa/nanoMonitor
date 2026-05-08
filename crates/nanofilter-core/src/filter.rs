use nanoseq_core::quality::mean_qv_from_fastq_ascii;

pub fn calculate_phred_avg(qual_string: &[u8]) -> f64 {
    mean_qv_from_fastq_ascii(qual_string)
}
