//! Live BCR-ABL1 engine.
//!
//! Watches a MinKNOW output directory, filters FASTQ exactly like
//! `nanostream filter`, aligns with `minimap2 -x splice` against the same hg19
//! reference nagini uses, and keeps one sorted, indexed BAM per barcode that can
//! be analysed while the run is going and later handed to nagini as
//! `--input-type aligned_bam`.

pub mod align;
pub mod config;
pub mod engine;
pub mod ledger;
pub mod reads;
pub mod replay;
pub mod state;
pub mod store;
pub mod watcher;

pub use config::EngineConfig;
pub use engine::{Engine, StepReport};
