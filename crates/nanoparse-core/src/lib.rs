use clap::ValueEnum;

pub mod enrichment;
pub mod matcher;
pub mod output;
pub mod pore_stats;
pub mod primers;
pub mod qv;
pub mod stats;

#[derive(Clone, Copy, ValueEnum, Debug, PartialEq)]
pub enum MatchMode {
    /// Semi-global alignment using triple_accel SIMD
    Semiglobal,
    /// Use mapping coordinates with learned cache (very fast)
    Coords,
}
