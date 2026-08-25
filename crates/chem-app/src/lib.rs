mod app;
mod dataset;
mod fingerprint_view;
mod molecule_view;
mod state;
mod structure_view;
mod svg;
mod task;
mod views;
mod windows;

pub use app::WorkbenchApp;

/// Id of the `<canvas>` element eframe renders into on web builds (see
/// crates/chem-app/index.html).
pub const CANVAS_ID: &str = "the_canvas_id";
