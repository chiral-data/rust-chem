use thiserror::Error;

#[derive(Error, Debug)]
pub enum SmilesError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    #[error("Unclosed ring: {0}")]
    UnclosedRing(u32),

    #[error("Invalid ring number: {0}")]
    InvalidRing(u32),

    #[error("Mismatched branches")]
    MismatchedBranches,
}
