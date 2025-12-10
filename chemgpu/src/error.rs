use thiserror::Error;

#[derive(Error, Debug)]
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
