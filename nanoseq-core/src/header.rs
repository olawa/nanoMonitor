use std::collections::HashMap;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct NanoporeMetadata {
    pub read_id: String,
    pub run_id: Option<String>,
    pub flow_cell_id: Option<String>,
    pub barcode: Option<String>,
    pub basecall_model: Option<String>,
    pub basecall_gpu: Option<String>,
    pub channel: Option<i64>,
    pub start_time: Option<String>,
}

impl NanoporeMetadata {
    pub fn read_group_id(&self) -> String {
        match (&self.run_id, &self.barcode) {
            (Some(run), Some(bc)) => format!("{}_{}", run, bc),
            (Some(run), None) => run.clone(),
            _ => self.read_id.clone(),
        }
    }
}

pub fn parse_nanopore_header(header: &[u8]) -> NanoporeMetadata {
    let header_str = String::from_utf8_lossy(header);
    let parts: Vec<&str> = header_str.split_whitespace().collect();

    let read_id = parts
        .first()
        .map(|s| s.trim_start_matches('@').to_string())
        .unwrap_or_default();

    let mut tags: HashMap<&str, &str> = HashMap::new();
    for part in parts.iter().skip(1) {
        if let Some((key, value)) = part.split_once('=') {
            tags.insert(key, value);
        }
    }

    NanoporeMetadata {
        read_id,
        run_id: tags.get("runid").map(|s| s.to_string()),
        flow_cell_id: tags.get("flow_cell_id").map(|s| s.to_string()),
        barcode: tags.get("barcode").map(|s| s.to_string()),
        basecall_model: tags.get("basecall_model_version_id").map(|s| s.to_string()),
        basecall_gpu: tags.get("basecall_gpu").map(|s| s.to_string()),
        channel: tags.get("ch").and_then(|s| s.parse::<i64>().ok()),
        start_time: tags.get("start_time").map(|s| s.to_string()),
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
