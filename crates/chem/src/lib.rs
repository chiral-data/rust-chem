//! Cheminformatics in Rust.
//!
//! Transitional: this re-exports the separate library crates so every consumer
//! can move to the final `chem::core::` spelling while the workspace still
//! compiles. The next commit replaces these re-exports with the real modules
//! and deletes the six crates.

pub use chemcore as core;
pub use chemdraw as draw;
pub use chemfp as fp;
pub use chemgpu as gpu;
pub use chemio as io;
pub use chemsearch as search;
