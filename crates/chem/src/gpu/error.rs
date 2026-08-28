use thiserror::Error;

#[derive(Error, Debug)]
/// Adding a variant to a public enum is a breaking change unless callers
/// are told not to match it exhaustively. This is the attribute that says
/// so, and it has to be present from the first published version: adding it
/// later invalidates every exhaustive match written against the earlier one.
#[non_exhaustive]
pub enum GpuError {
    #[error("No suitable GPU adapter found")]
    NoAdapter,

    #[error("Failed to request device: {0}")]
    DeviceRequest(String),

    #[error("GPU operation failed: {0}")]
    OperationFailed(String),

    #[error("Buffer error: {0}")]
    BufferError(String),

    #[error("Shader compilation failed: {0}")]
    ShaderError(String),

    #[error("Invalid buffer size: expected {expected}, got {actual}")]
    InvalidBufferSize { expected: usize, actual: usize },
}
