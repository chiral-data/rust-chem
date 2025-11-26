use thiserror::Error;

#[derive(Error, Debug)]
pub enum FingerprintError {
    #[error("Invalid molecule")]
    InvalidMolecule,

    #[error("Invalid radius: {0}")]
    InvalidRadius(u32),

    #[error("Invalid bit size: {0}")]
    InvalidBitSize(usize),

    #[error("Fingerprint size mismatch: {0} vs {1}")]
    SizeMismatch(usize, usize),
}
