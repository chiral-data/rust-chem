pub mod buffers;
pub mod context;
pub mod error;
pub mod morgan;
pub mod pipeline;
pub mod tanimoto;

pub use buffers::BufferManager;
pub use context::GpuContext;
pub use error::GpuError;
pub use morgan::{GpuMorganFingerprint, MoleculeEncoder};
pub use pipeline::ComputePipeline;
pub use tanimoto::GpuTanimoto;

/// Initialize logging for GPU operations.
pub fn init_logging() {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();
}
