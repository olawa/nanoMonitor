mod app;
mod bcr;
mod cli;
mod integrations;
mod model;
mod nanostream_cli;
mod remote;

use app::{AppStartupConfig, NanoMonitorApp};
use clap::Parser;
use cli::{CliMode, NanoMonitorCli};
use model::AnalysisMode;

fn run_bcr(args: NanoMonitorCli) -> Result<(), eframe::Error> {
    let setup = bcr::proc::BcrSetup {
        input_dir: args.monitor_dir.or(args.input).unwrap_or_default(),
        out_dir: args.bcr_out.unwrap_or_default(),
        reference: args.reference.unwrap_or_default(),
        nagini_dir: args.nagini.unwrap_or_default(),
        nanostream_bin: args.nanostream_bin.unwrap_or_default(),
        ..Default::default()
    };
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1500.0, 950.0])
            .with_title("nanoMonitor - BCR-ABL1"),
        ..Default::default()
    };
    eframe::run_native(
        "nanoMonitor BCR-ABL1",
        options,
        Box::new(move |_cc| {
            struct Wrapper(bcr::BcrApp);
            impl eframe::App for Wrapper {
                fn update(&mut self, ctx: &eframe::egui::Context, _: &mut eframe::Frame) {
                    self.0.show(ctx);
                }
            }
            Ok(Box::new(Wrapper(bcr::BcrApp::new(bcr::BcrStartup { setup }))))
        }),
    )
}

fn main() -> Result<(), eframe::Error> {
    let args = NanoMonitorCli::parse();
    if matches!(args.mode, Some(CliMode::BcrAbl)) {
        return run_bcr(args);
    }
    let startup = AppStartupConfig {
        mode: args.mode.map(|m| match m {
            CliMode::Amplicon | CliMode::BcrAbl => AnalysisMode::Amplicon,
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
