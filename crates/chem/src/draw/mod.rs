//! 2D molecule depiction, with no GUI attached.
//!
//! Two layers:
//!
//! - [`structure`] describes a molecule as a list of drawable primitives.
//! - [`svg`] is one backend for that description, writing an SVG document.
//!
//! A GUI backend lives in the application rather than here, since painting is
//! the part that needs a toolkit. Anything else — a CLI, a server, a test
//! asserting on a depiction — can depend on this crate alone.

pub mod structure;
pub mod svg;

/// The geometry and colour types this module's signatures use.
///
/// They belong to `emath` and `ecolor`, which egui publishes separately and
/// re-exports — so a GUI consumer gets egui's own types back with no conversion,
/// while a command-line one links no toolkit at all.
///
/// Re-exported because a caller must be able to *name* the types a public
/// function takes. [`svg::structure_to_svg`] needs a [`Vec2`] for its viewport,
/// and without this a consumer of `chem` would have to add `emath` as a direct
/// dependency and match its version by hand.
pub use ecolor::Color32;
pub use emath::{Align2, Pos2, Rect, Vec2};

pub use structure::{
    AtomVisualization, MeasureText, ShowCarbons, StructureOptions, StructureShape, StructureTheme,
    describe_structure,
};
pub use svg::{structure_to_svg, suggested_filename};
