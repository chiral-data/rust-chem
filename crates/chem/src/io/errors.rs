use thiserror::Error;

#[derive(Error, Debug)]
/// Adding a variant to a public enum is a breaking change unless callers
/// are told not to match it exhaustively. This is the attribute that says
/// so, and it has to be present from the first published version: adding it
/// later invalidates every exhaustive match written against the earlier one.
#[non_exhaustive]
pub enum SmilesError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    /// Every ring left open, sorted.
    ///
    /// All of them rather than one: `c1ccnc2ccccc12` leaves rings 1 and 2 open,
    /// and naming one tells the reader least about a molecule with two problems.
    /// It used to name an arbitrary `HashMap` key, so the same input reported a
    /// different number between runs (#153).
    #[error("Unclosed rings: {}", .0.iter().map(u32::to_string).collect::<Vec<_>>().join(", "))]
    UnclosedRings(Vec<u32>),

    #[error("Invalid ring number: {0}")]
    InvalidRing(u32),

    #[error("Mismatched branches")]
    MismatchedBranches,

    /// A bond symbol with nothing on one side of it — `=CC` or `CC=`.
    ///
    /// These used to be dropped silently, so both parsed as ethane: a real
    /// molecule, quietly different from the one written. A rejection is worse
    /// than a correct parse and far better than a plausible wrong answer (#155).
    #[error("A bond has no atom to attach to")]
    DanglingBond,

    /// The input tokenized but described no atoms.
    ///
    /// Bond and branch characters are legal tokens on their own, so a string
    /// made only of them used to parse "successfully" into a molecule with
    /// nothing in it. `$$$$` is the case that matters — it is the SDF record
    /// terminator, and an SDF read as SMILES produced one atomless molecule per
    /// record instead of failing.
    #[error("No atoms in SMILES")]
    NoAtoms,

    /// A `.` with no component on one side of it — `CC.` or `.CC`.
    ///
    /// Named rather than folded into [`SmilesError::DanglingBond`]: a dot is
    /// not a bond, and "a bond has no atom to attach to" sends the reader
    /// hunting for a `=` that was never there.
    #[error("A dot has no component on one side of it")]
    DanglingDot,
}

#[derive(Error, Debug)]
#[non_exhaustive]
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
