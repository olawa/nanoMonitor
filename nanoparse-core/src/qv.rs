use nanoseq_core::quality::select_qv;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;

fn aux_to_f64(aux: Aux<'_>) -> Option<f64> {
    match aux {
        Aux::Float(v) => Some(v as f64),
        Aux::Double(v) => Some(v),
        Aux::I8(v) => Some(f64::from(v)),
        Aux::U8(v) => Some(f64::from(v)),
        Aux::I16(v) => Some(f64::from(v)),
        Aux::U16(v) => Some(f64::from(v)),
        Aux::I32(v) => Some(v as f64),
        Aux::U32(v) => Some(v as f64),
        _ => None,
    }
}

pub fn qv_from_aux_and_phred_scores(qs_aux: Option<Aux<'_>>, phred_scores: &[u8]) -> f32 {
    let explicit_qs = qs_aux.and_then(aux_to_f64);
    select_qv(explicit_qs, phred_scores) as f32
}

pub fn qv_from_record(record: &bam::Record) -> f32 {
    let qs_aux = record.aux(b"QS").ok().or_else(|| record.aux(b"qs").ok());
    qv_from_aux_and_phred_scores(qs_aux, record.qual())
}
