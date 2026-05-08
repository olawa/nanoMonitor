use crate::model::FilterConfig;

#[derive(Debug, Clone)]
pub struct NanostreamConfig {
    pub executable: String,
    pub primers_path: String,
    pub threads: usize,
    pub primer_tolerance: i64,
}

impl Default for NanostreamConfig {
    fn default() -> Self {
        Self {
            executable: "nanostream".into(),
            primers_path: "primers.tsv".into(),
            threads: 8,
            primer_tolerance: 50,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CommandSpec {
    pub program: String,
    pub args: Vec<String>,
}

impl CommandSpec {
    pub fn as_shell_line(&self) -> String {
        let mut out = self.program.clone();
        for arg in &self.args {
            if arg.contains(' ') {
                out.push(' ');
                out.push('"');
                out.push_str(arg);
                out.push('"');
            } else {
                out.push(' ');
                out.push_str(arg);
            }
        }
        out
    }
}

impl NanostreamConfig {
    pub fn build_amplicon_command(
        &self,
        bam_path: &str,
        filters: &FilterConfig,
        reference: Option<&str>,
        gtf: Option<&str>,
    ) -> CommandSpec {
        let mut args = vec![
            "amplicons".into(),
            bam_path.into(), // Positional in nanostream
            "--primers".into(),
            self.primers_path.clone(),
            "--threads".into(),
            self.threads.to_string(),
            "--primer-tolerance".into(),
            self.primer_tolerance.to_string(),
            "--min-qs".into(),
            format!("{}", filters.min_qs),
            "--min-len".into(),
            format!("{}", filters.min_len),
            "--max-reads".into(),
            format!("{}", filters.max_reads),
            "--output".into(),
            "-".into(),
            "--end-length".into(),
            "150".into(),
            "--max-edit-dist".into(),
            "3".into(),
        ];
        if filters.duplex_only {
            args.push("--duplex-only".into());
        }
        if let Some(reference) = reference {
            if !reference.trim().is_empty() {
                args.push("--reference".into());
                args.push(reference.to_string());
            }
        }
        if let Some(gtf) = gtf {
            if !gtf.trim().is_empty() {
                args.push("--gtf".into());
                args.push(gtf.to_string());
            }
        }

        CommandSpec {
            program: self.executable.clone(),
            args,
        }
    }
}
