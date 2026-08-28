use thiserror::Error;

#[derive(Error, Debug)]
/// Adding a variant to a public enum is a breaking change unless callers
/// are told not to match it exhaustively. This is the attribute that says
/// so, and it has to be present from the first published version: adding it
/// later invalidates every exhaustive match written against the earlier one.
#[non_exhaustive]
pub enum FingerprintError {
    #[error("Invalid molecule")]
    InvalidMolecule,

    #[error("Invalid radius: {0}")]
    InvalidRadius(u32),

    #[error("Invalid bit size: {0}")]
    InvalidBitSize(usize),

    #[error("Fingerprint size mismatch: {0} vs {1}")]
    SizeMismatch(usize, usize),

    #[error("Molecule has no atoms")]
    EmptyMolecule,
}
