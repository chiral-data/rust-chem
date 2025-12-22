mod app;
mod dataset;
mod fingerprint_view;
mod molecule_view;
mod search;

use app::ChemFpDemoApp;

fn main() -> eframe::Result {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_title("ChemFP Demo - Molecular Fingerprint Search"),
        ..Default::default()
    };

    eframe::run_native(
        "ChemFP Demo",
        options,
        Box::new(|cc| Ok(Box::new(ChemFpDemoApp::new(cc)))),
    )
}
