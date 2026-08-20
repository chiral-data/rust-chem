//! The app's views, each owning its own UI state.
//!
//! A view is a struct holding what only it cares about — a text box's contents,
//! a slider's value, which row is expanded — plus a `ui` method that draws it
//! against a `&mut AppState` for everything shared. Nothing here reaches into
//! another view.
//!
//! The module names are the windows v0.5.0 is heading for
//! ([tracking issue #112]), not the panels that exist today: the four files
//! below still hold today's panel boundaries, and #104/#105/#106 move content
//! between them as the windows take shape. `operations.rs` in particular does
//! not yet hold the fingerprint controls — those are still drawn by the data
//! panel that has always drawn them, reading parameters the operations view
//! owns, until #105 moves them.
//!
//! [tracking issue #112]: https://github.com/chiral-data/rust-chem/issues/112

pub mod data_sources;
pub mod detail;
pub mod operations;
pub mod visualization;

pub use data_sources::DataSourcesView;
pub use detail::DetailView;
pub use operations::OperationsView;
pub use visualization::VisualizationView;
