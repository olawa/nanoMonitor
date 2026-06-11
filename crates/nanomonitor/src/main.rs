mod app;
mod cli;
mod integrations;
mod model;
mod nanostream_cli;
mod remote;

use app::{AppStartupConfig, NanoMonitorApp};
use clap::Parser;
use cli::{CliMode, NanoMonitorCli};
use model::AnalysisMode;

fn main() -> Result<(), eframe::Error> {
    let args = NanoMonitorCli::parse();
    let startup = AppStartupConfig {
        mode: args.mode.map(|m| match m {
            CliMode::Amplicon => AnalysisMode::Amplicon,
        }),
        input_path: args.input,
        monitor_dir: args.monitor_dir,
        reference_path: args.reference,
        gtf_path: args.gtf,
        primers_path: args.primers,
        nanostream_bin: args.nanostream_bin,
        run_on_start: args.start,
    };

    let native_options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1600.0, 980.0])
            .with_title("nanoMonitor - PCR Amplicon Analysis"),
        ..Default::default()
    };

    eframe::run_native(
        "nanoMonitor",
        native_options,
        Box::new(move |cc| Ok(Box::new(NanoMonitorApp::new(cc, startup.clone())))),
    )
}
