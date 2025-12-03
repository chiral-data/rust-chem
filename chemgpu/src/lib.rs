pub mod buffers;
pub mod context;
pub mod error;
pub mod pipeline;

pub use buffers::BufferManager;
pub use context::GpuContext;
pub use error::GpuError;
pub use pipeline::ComputePipeline;

/// Initialize logging for GPU operations.
pub fn init_logging() {
    env_logger::builder()
        .filter_level(log::LevelFilter::Info)
        .init();
}
