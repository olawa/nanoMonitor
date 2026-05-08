pub fn is_fastq_path(path: &str) -> bool {
    let p = path.to_ascii_lowercase();
    p.ends_with(".fastq") || p.ends_with(".fq") || p.ends_with(".fastq.gz") || p.ends_with(".fq.gz")
}

pub fn trim_line_ending(s: &str) -> &str {
    s.trim_end_matches(['\n', '\r'])
}

#[cfg(test)]
mod tests {
    use super::{is_fastq_path, trim_line_ending};

    #[test]
    fn detects_fastq_paths() {
        assert!(is_fastq_path("a.fastq"));
        assert!(is_fastq_path("a.fq.gz"));
        assert!(!is_fastq_path("a.bam"));
    }

    #[test]
    fn trims_line_endings() {
        assert_eq!(trim_line_ending("abc\r\n"), "abc");
        assert_eq!(trim_line_ending("abc\n"), "abc");
        assert_eq!(trim_line_ending("abc"), "abc");
    }
}
