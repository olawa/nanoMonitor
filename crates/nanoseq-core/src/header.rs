use std::collections::HashMap;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct NanoporeMetadata {
    pub read_id: String,
    pub run_id: Option<String>,
    pub flow_cell_id: Option<String>,
    pub barcode: Option<String>,
    pub basecall_model: Option<String>,
    pub basecall_gpu: Option<String>,
    pub channel: Option<i64>,
    pub start_time: Option<String>,
    /// The basecaller's own mean read quality, from the `qs` tag. Preferred
    /// over recomputing from the quality string: it is what dorado reports and
    /// what a `--qv` threshold is normally understood to mean.
    pub mean_qscore: Option<f64>,
    /// The read this one was split from (`parent_read_id=` / `pi:Z:`). Equal to
    /// `read_id`, or absent, for a read that was not split.
    pub parent_read_id: Option<String>,
}

impl NanoporeMetadata {
    /// Whether dorado split this read out of a longer one. Dorado computes the
    /// mean quality of a split read over all of it, not from base 60.
    pub fn is_split_read(&self) -> bool {
        self.parent_read_id
            .as_deref()
            .map(|parent| !parent.is_empty() && parent != self.read_id)
            .unwrap_or(false)
    }

    pub fn read_group_id(&self) -> String {
        match (&self.run_id, &self.barcode) {
            (Some(run), Some(bc)) => format!("{}_{}", run, bc),
            (Some(run), None) => run.clone(),
            _ => self.read_id.clone(),
        }
    }
}

/// The last `_`-separated field of a dorado read group when it is a barcode or
/// alias label (`barcode01`, `unclassified`, `RB12`, `NB03`, ...), not model text.
fn looks_like_barcode_suffix(field: &str) -> bool {
    let lower = field.to_ascii_lowercase();
    lower == "unclassified"
        || lower == "no_barcode"
        || ["barcode", "rb", "nb", "bc"].iter().any(|p| {
            lower
                .strip_prefix(p)
                .map(|rest| !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()))
                .unwrap_or(false)
        })
}

pub fn parse_nanopore_header(header: &[u8]) -> NanoporeMetadata {
    let header_str = String::from_utf8_lossy(header);
    let parts: Vec<&str> = header_str.split_whitespace().collect();

    let read_id = parts
        .first()
        .map(|s| s.trim_start_matches('@').to_string())
        .unwrap_or_default();

    // Two header dialects are in the wild and they do not overlap:
    //
    //   guppy / older dorado   runid=... ch=566 start_time=... barcode=...
    //   MinKNOW / dorado now   qs:f:18.93 ch:i:566 st:Z:... al:Z:barcode01
    //
    // Parsing only the first meant every metadata-driven feature degraded
    // silently on current MinKNOW output: no qs, so --qv recomputed instead of
    // using the basecaller's own score, and no channel or start time, so
    // --channel-range and the time window matched nothing.
    let mut tags: HashMap<&str, &str> = HashMap::new();
    for part in parts.iter().skip(1) {
        let mut fields = part.splitn(3, ':');
        if let (Some(tag), Some(kind), Some(value)) = (fields.next(), fields.next(), fields.next())
        {
            if matches!(kind, "i" | "f" | "Z" | "A" | "H" | "B") {
                tags.insert(tag, value);
                continue;
            }
        }
        if let Some((key, value)) = part.split_once('=') {
            tags.insert(key, value);
        }
    }

    let first = |keys: &[&str]| -> Option<String> {
        keys.iter()
            .find_map(|key| tags.get(key))
            .map(|value| value.to_string())
    };

    // RG is "<acquisition id>_<model>_<barcode>" in dorado's output, so it
    // carries the run id and model when the standalone tags are absent.
    let read_group = tags.get("RG").copied().unwrap_or_default();
    let run_id_from_rg = read_group.split('_').next().filter(|s| !s.is_empty());
    // The trailing barcode field is optional, so only strip it when it looks
    // like one; otherwise the last model field (e.g. `sup@v5.2.0`) is lost.
    let model_from_rg = read_group
        .split_once('_')
        .map(|(_, rest)| rest)
        .map(|rest| match rest.rsplit_once('_') {
            Some((model, tail)) if looks_like_barcode_suffix(tail) => model,
            _ => rest,
        })
        .filter(|s| !s.is_empty());

    NanoporeMetadata {
        read_id,
        run_id: first(&["runid"]).or_else(|| run_id_from_rg.map(|s| s.to_string())),
        flow_cell_id: first(&["flow_cell_id", "PU"]),
        barcode: first(&["barcode", "BC", "al", "SM"]),
        basecall_model: first(&["basecall_model_version_id"])
            .or_else(|| model_from_rg.map(|s| s.to_string())),
        basecall_gpu: first(&["basecall_gpu"]),
        channel: tags.get("ch").and_then(|s| s.parse::<i64>().ok()),
        start_time: first(&["start_time", "st"]),
        mean_qscore: tags.get("qs").and_then(|s| s.parse::<f64>().ok()),
        parent_read_id: first(&["parent_read_id", "pi"]),
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_nanopore_header, NanoporeMetadata};

    #[test]
    fn parses_nanopore_header() {
        let header = b"@57de3efd-e78b-4492-9387-e891d23d604b runid=9a421bdc flow_cell_id=FBE72298 barcode=barcode01 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.2.0 basecall_gpu=RTX4090 ch=1766 start_time=2026-03-17T12:00:00Z";
        let meta = parse_nanopore_header(header);

        assert_eq!(meta.read_id, "57de3efd-e78b-4492-9387-e891d23d604b");
        assert_eq!(meta.run_id.as_deref(), Some("9a421bdc"));
        assert_eq!(meta.flow_cell_id.as_deref(), Some("FBE72298"));
        assert_eq!(meta.barcode.as_deref(), Some("barcode01"));
        assert_eq!(
            meta.basecall_model.as_deref(),
            Some("dna_r10.4.1_e8.2_400bps_sup@v5.2.0")
        );
        assert_eq!(meta.basecall_gpu.as_deref(), Some("RTX4090"));
        assert_eq!(meta.channel, Some(1766));
        assert_eq!(meta.start_time.as_deref(), Some("2026-03-17T12:00:00Z"));
    }

    #[test]
    fn parses_minknow_sam_tag_header() {
        // Verbatim from a MinKNOW run (bcr_abl_exp7_9, PBE82562).
        let header = b"@afb174c8-3078-4db2-bc17-dea9a9b09f55 qs:f:18.9348 mx:i:4 ch:i:566 rn:i:375 st:Z:2026-08-28T10:15:56.552323+00:00 ns:i:36251 ts:i:10 sm:f:93.6924 sd:f:23.506744 sv:Z:pa du:f:7.2502 dx:i:0 RG:Z:ef769e0a-6c2d-4f1f-9a0c-3b8e5d2a7c41_dna_r10.4.1_e8.2_400bps_sup@v5.2.0_barcode01 DT:Z:2026-08-28T10:15:56.552323+00:00 PU:Z:PBE82562 LB:Z:bcr_abl_exp7_9 SM:Z:barcode01 al:Z:barcode01";
        let meta = parse_nanopore_header(header);

        assert_eq!(meta.read_id, "afb174c8-3078-4db2-bc17-dea9a9b09f55");
        assert_eq!(meta.mean_qscore, Some(18.9348));
        assert_eq!(meta.channel, Some(566));
        assert_eq!(
            meta.start_time.as_deref(),
            Some("2026-08-28T10:15:56.552323+00:00")
        );
        assert_eq!(meta.barcode.as_deref(), Some("barcode01"));
        assert_eq!(meta.flow_cell_id.as_deref(), Some("PBE82562"));
        assert_eq!(
            meta.run_id.as_deref(),
            Some("ef769e0a-6c2d-4f1f-9a0c-3b8e5d2a7c41")
        );
        assert_eq!(
            meta.basecall_model.as_deref(),
            Some("dna_r10.4.1_e8.2_400bps_sup@v5.2.0")
        );
    }

    #[test]
    fn read_group_without_barcode_keeps_the_whole_model() {
        let header = b"@r1 RG:Z:ef769e0a-6c2d-4f1f-9a0c-3b8e5d2a7c41_dna_r10.4.1_e8.2_400bps_sup@v5.2.0";
        let meta = parse_nanopore_header(header);
        assert_eq!(
            meta.basecall_model.as_deref(),
            Some("dna_r10.4.1_e8.2_400bps_sup@v5.2.0")
        );
        assert_eq!(
            meta.run_id.as_deref(),
            Some("ef769e0a-6c2d-4f1f-9a0c-3b8e5d2a7c41")
        );
    }

    #[test]
    fn legacy_header_has_no_qscore_tag() {
        let header = b"@read1 runid=9a421bdc ch=12 start_time=2026-03-17T12:00:00Z";
        let meta = parse_nanopore_header(header);
        assert_eq!(meta.mean_qscore, None);
        assert_eq!(meta.channel, Some(12));
        // The timestamp's own colons must not be mistaken for a SAM tag.
        assert_eq!(meta.start_time.as_deref(), Some("2026-03-17T12:00:00Z"));
    }

    #[test]
    fn recognises_a_split_read_in_either_header_dialect() {
        let minknow = parse_nanopore_header(b"@child runid=r parent_read_id=parent barcode=barcode01");
        assert!(minknow.is_split_read());
        let dorado = parse_nanopore_header(b"@child qs:f:18.0 pi:Z:parent");
        assert!(dorado.is_split_read());
        let unsplit = parse_nanopore_header(b"@read1 runid=r parent_read_id=read1");
        assert!(!unsplit.is_split_read());
        assert!(!parse_nanopore_header(b"@read1 qs:f:18.0").is_split_read());
    }

    #[test]
    fn builds_read_group_id() {
        let meta = NanoporeMetadata {
            read_id: "read1".to_string(),
            run_id: Some("run123".to_string()),
            barcode: Some("bc01".to_string()),
            ..Default::default()
        };
        assert_eq!(meta.read_group_id(), "run123_bc01");
    }
}
