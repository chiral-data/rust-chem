pub mod buffers;
pub mod context;
pub mod error;
pub mod morgan;
pub mod pipeline;
pub mod tanimoto;

pub use buffers::BufferManager;
pub use context::GpuContext;
pub use error::GpuError;
pub use morgan::GpuMorganFingerprint;
pub use pipeline::ComputePipeline;
pub use tanimoto::{GpuTanimoto, GpuTargetSet};

// `init_logging()` used to live here, calling `env_logger::builder().init()`.
// That panics if a logger is already installed, which makes it a landmine in a
// library — the caller cannot know whether something else got there first.
// Installing a logger is the application's job; its only caller was an example,
// which now does it itself.
