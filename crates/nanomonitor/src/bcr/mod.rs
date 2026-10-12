//! BCR-ABL1 live mode: run setup, live per-barcode monitoring and the
//! preliminary analysis, all read from the engine's output directory.

mod final_run;
pub mod data;
pub mod proc;
mod view;

pub use view::{BcrApp, BcrStartup};
