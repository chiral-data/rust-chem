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

pub use structure::{
    AtomVisualization, MeasureText, ShowCarbons, StructureOptions, StructureShape, StructureTheme,
    describe_structure,
};
pub use svg::{structure_to_svg, suggested_filename};
