mod app;
mod dataset;
mod fingerprint_view;
mod molecule_view;
mod search;
mod state;
mod structure_view;
mod task;
mod views;
mod windows;

pub use app::WorkbenchApp;

/// Id of the `<canvas>` element eframe renders into on web builds (see
/// chem_app/index.html).
pub const CANVAS_ID: &str = "the_canvas_id";
