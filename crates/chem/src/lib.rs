#![doc = include_str!("../README.md")]
//!
//! # A note for contributors
//!
//! This crate has a root module named `core`, which collides with Rust's
//! built-in `core` crate. Inside this crate, write `::core::` with the leading
//! colons to reach the built-in. Nothing here needs it today — everything goes
//! through `std` — and the collision is a compile error rather than a silent
//! misresolution, so it announces itself.

#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod core;
pub mod draw;
pub mod fp;
pub mod io;
pub mod search;

#[cfg(feature = "gpu")]
#[cfg_attr(docsrs, doc(cfg(feature = "gpu")))]
pub mod gpu;

pub mod prelude {
    pub use crate::core::prelude::*;
}
