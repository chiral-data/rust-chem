//! The app's views, each owning its own UI state.
//!
//! A view is a struct holding what only it cares about — a text box's contents,
//! a slider's value, which row is expanded — plus a `ui` method that draws it
//! against a `&mut AppState` for everything shared. Nothing here reaches into
//! another view.
//!
//! One module per window, plus `detail` for the molecule windows, of which
//! there may be several at once.
//!
//! Nothing here reaches into another view. That was not free: moving the query's
//! structure into the Inspector needed `AppState` to remember the SMILES it
//! parsed (#106), because the alternative was reading the text box that belongs
//! to `operations`. Where two views need the same fact, the fact belongs in
//! [`crate::state::AppState`].

pub mod datasets;
pub mod detail;
pub mod inspector;
pub mod operations;
pub mod settings;

pub use datasets::DatasetsView;
pub use detail::DetailView;
pub use inspector::InspectorView;
pub use operations::OperationsView;
pub use settings::SettingsView;
