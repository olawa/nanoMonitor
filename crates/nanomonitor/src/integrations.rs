use crate::model::VariantRow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct ClinicalMutation {
    pub chrom: String,
    pub start: u64,
    pub stop: u64,
    pub bas_variant: String,
    pub aa_variant: String,
    pub exon: String,
}

/// Loads clinical_mutations.tsv from the workspace root or standard relative paths.
pub fn load_clinical_mutations() -> HashMap<(String, u64), ClinicalMutation> {
    let mut db = HashMap::new();
    
    // Look in current working directory and common locations
    let paths = vec![
        PathBuf::from("clinical_mutations.tsv"),
        PathBuf::from("../../clinical_mutations.tsv"),
        PathBuf::from("crates/nanomonitor/clinical_mutations.tsv"),
    ];
    
    let mut resolved_path = None;
    for p in paths {
        if p.exists() && p.is_file() {
            resolved_path = Some(p);
            break;
        }
    }
    
    let path = match resolved_path {
        Some(p) => p,
        None => return db, // Return empty if not found
    };
    
    if let Ok(file) = File::open(path) {
        let reader = BufReader::new(file);
        // Header: Chr	Start	Stop	Bas-variant	AA-variant	Exon
        for (idx, line) in reader.lines().flatten().enumerate() {
            if idx == 0 || line.trim().is_empty() {
                continue; // Skip header and empty lines
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let chrom = parts[0].trim().to_string();
                let start: u64 = parts[1].trim().parse().unwrap_or(0);
                let stop: u64 = parts[2].trim().parse().unwrap_or(0);
                let bas_variant = parts[3].trim().to_string();
                let aa_variant = parts[4].trim().to_string();
                let exon = parts.get(5).map(|s| s.trim().to_string()).unwrap_or_default();
                
                db.insert((chrom.clone(), start), ClinicalMutation {
                    chrom,
                    start,
                    stop,
                    bas_variant,
                    aa_variant,
                    exon,
                });
            }
        }
    }
    db
}

/// Helper to discover bin path or fall back to system PATH
fn resolve_binary_path(user_configured: &str, default_executable: &str, sibling_path: &str) -> String {
    let configured = user_configured.trim();
    if !configured.is_empty() {
        if Path::new(configured).exists() {
            return configured.to_string();
        }
    }

    // Try sibling workspace release directory
    let cwd = std::env::current_dir().unwrap_or_default();
    let sibling_abs = cwd.join(sibling_path);
    if sibling_abs.exists() && sibling_abs.is_file() {
        return sibling_abs.to_string_lossy().to_string();
    }

    // Check system path fallback
    default_executable.to_string()
}

/// Runs `rs-qc snap` to render a regional alignment screenshot.
pub fn run_rs_qc_snap(
    rs_qc_bin: &str,
    bam_path: &str,
    region: &str,
    gtf_path: &str,
    reference_path: &str,
    output_png: &str,
) -> Result<String, String> {
    // Resolve binary path
    let bin = resolve_binary_path(
        rs_qc_bin,
        "rs-qc",
        "../rs-qc/target/release/rs-qc",
    );

    let mut cmd = Command::new(&bin);
    cmd.arg("snap")
       .arg("-i")
       .arg(bam_path)
       .arg("-r")
       .arg(region)
       .arg("-o")
       .arg(output_png);

    if !gtf_path.trim().is_empty() {
        cmd.arg("-a").arg(gtf_path);
    }
    if !reference_path.trim().is_empty() {
        cmd.arg("--reference").arg(reference_path);
    }

    let output = cmd.output().map_err(|e| format!("Failed to execute rs-qc snap: {}", e))?;
    if !output.status.success() {
        return Err(format!(
            "rs-qc snap failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    Ok(format!("Generated snapshot at: {}", output_png))
}

/// Runs `rindels` variant calling.
pub fn run_rindels(
    rindels_bin: &str,
    bam_path: &str,
    reference_path: &str,
    bed_path: &str,
    output_vcf: &str,
) -> Result<String, String> {
    // Resolve binary path
    let bin = resolve_binary_path(
        rindels_bin,
        "rindels",
        "../rindels/target/release/rindels",
    );

    let mut cmd = Command::new(&bin);
    cmd.arg("-b")
       .arg(bam_path)
       .arg("-o")
       .arg(output_vcf)
       .arg("--single-end"); // Nanopore reads are single-end in this pipeline context

    if !reference_path.trim().is_empty() {
        cmd.arg("-r").arg(reference_path);
    }
    if !bed_path.trim().is_empty() {
        cmd.arg("--bed").arg(bed_path);
    }

    let output = cmd.output().map_err(|e| format!("Failed to execute rindels: {}", e))?;
    if !output.status.success() {
        return Err(format!(
            "rindels failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    Ok(format!("Variant calling complete. VCF saved to: {}", output_vcf))
}

/// Parses output VCF and structures it into VariantRow representations with clinical annotation lookup.
pub fn parse_vcf(
    vcf_path: &str,
    clinical_db: &HashMap<(String, u64), ClinicalMutation>,
) -> Result<Vec<VariantRow>, String> {
    let path = Path::new(vcf_path);
    if !path.exists() || !path.is_file() {
        return Err(format!("VCF file not found: {}", vcf_path));
    }

    let file = File::open(path).map_err(|e| format!("Failed to open VCF: {}", e))?;
    let reader = BufReader::new(file);
    let mut variants = Vec::new();

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Failed to read VCF line: {}", e))?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue; // Skip headers & whitespace
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 8 {
            continue;
        }

        let chrom = parts[0].trim().to_string();
        let position: u64 = parts[1].trim().parse().unwrap_or(0);
        let ref_allele = parts[3].trim().to_string();
        let alt_allele = parts[4].trim().to_string();
        let info = parts[7].trim();

        // Extract DP (depth) and AF/VAF (allele frequency) from the INFO column
        let mut depth = 0;
        let mut af = 0.0;

        for item in info.split(';') {
            if item.starts_with("DP=") {
                if let Some(val_str) = item.strip_prefix("DP=") {
                    depth = val_str.parse().unwrap_or(0);
                }
            } else if item.starts_with("AF=") {
                if let Some(val_str) = item.strip_prefix("AF=") {
                    af = val_str.parse().unwrap_or(0.0);
                }
            }
        }

        let vaf = af * 100.0; // Convert fraction to percentage

        // Perform Clinical Annotations lookup
        let mut clinvar = "Unknown / VUS".to_string();
        let mut verdict = if vaf >= 10.0 {
            "High Confidence Somatic".to_string()
        } else if vaf >= 1.0 {
            "Somatic Candidate".to_string()
        } else {
            "Sub-clonal / Noise".to_string()
        };

        // Match against clinical mutation hotspots database
        if let Some(clinical) = clinical_db.get(&(chrom.clone(), position)) {
            clinvar = clinical.aa_variant.clone();
            let mut clinical_verdict = if !clinical.exon.is_empty() {
                clinical.exon.clone()
            } else {
                "Known Clinical Hotspot".to_string()
            };
            // Append visual marker
            clinical_verdict = format!("★ HOTSPOT ({})", clinical_verdict);
            verdict = clinical_verdict;
        }

        variants.push(VariantRow {
            chrom,
            position,
            ref_allele,
            alt_allele,
            depth,
            vaf,
            clinvar,
            verdict,
        });
    }

    // Sort variants by position
    variants.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| a.position.cmp(&b.position)));

    Ok(variants)
}
