//! Cheminformatics in Rust: molecules, SMILES and SDF, fingerprints, similarity
//! search, and 2D depiction.
//!
//! ```no_run
//! use chem::io::smiles::parse_smiles;
//! use chem::fp::morgan::MorganFingerprint;
//!
//! let mol = parse_smiles("c1ccccc1O")?;
//! let fp = MorganFingerprint::get_fingerprint_as_bitvec(
//!     &mol, 2, 2048, None, None, false, true, false,
//! )?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Features
//!
//! | feature | brings in |
//! | --- | --- |
//! | *(default)* | [`core`], [`io`], [`fp`], [`draw`], [`search`] — all CPU, no heavy dependencies |
//! | `gpu` | [`gpu`], plus the GPU paths inside [`search`]. Pulls `wgpu`. |
//! | `cli` | the `chem` binary. Implies `gpu`. |
//!
//! `gpu` is off by default because `wgpu` is a large tree and irrelevant to
//! anyone parsing SMILES or fingerprinting on a CPU. Without it, `search` still
//! offers its whole API — it simply never finds a device, and says so through
//! [`search::FingerprintSearch::gpu_init_error`].
//!
//! # A note for contributors
//!
//! This crate has a root module named `core`, which collides with Rust's
//! built-in `core` crate. Inside this crate, write `::core::` with the leading
//! colons to reach the built-in. Nothing here needs it today — everything goes
//! through `std` — and the collision is a compile error rather than a silent
//! misresolution, so it announces itself.

pub mod core;
pub mod draw;
pub mod fp;
pub mod io;
pub mod search;

#[cfg(feature = "gpu")]
pub mod gpu;

pub mod prelude {
    pub use crate::core::prelude::*;
}
