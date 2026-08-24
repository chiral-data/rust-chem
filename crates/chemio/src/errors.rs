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

    /// The input tokenized but described no atoms.
    ///
    /// Bond and branch characters are legal tokens on their own, so a string
    /// made only of them used to parse "successfully" into a molecule with
    /// nothing in it. `$$$$` is the case that matters — it is the SDF record
    /// terminator, and an SDF read as SMILES produced one atomless molecule per
    /// record instead of failing.
    #[error("No atoms in SMILES")]
    NoAtoms,
}

#[derive(Error, Debug)]
pub enum SdfError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    #[error("Invalid bond line: {0}")]
    InvalidBondLine(String),

    #[error("Missing counts line")]
    MissingCountsLine,
}
