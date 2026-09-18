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

    /// A bond symbol immediately followed by another one, with no atom
    /// between them — `C##C`, `C###C`, `C==C`.
    ///
    /// These used to collapse silently into whichever symbol arrived last, so
    /// `C##C` and `C###C` both parsed as ethyne. Silent acceptance of
    /// malformed input is the failure mode hardest to notice downstream: a
    /// typo in a generated file becomes a plausible molecule rather than a
    /// reported skip (#190).
    #[error("Two bond symbols in a row with no atom between them")]
    RepeatedBondSymbol,

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

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum XyzError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    #[error("Declared {expected} atoms but found {got}")]
    AtomCountMismatch { expected: usize, got: usize },
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PdbError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    /// No `ATOM`/`HETATM` record was ever recognized, so zero atoms were
    /// read.
    ///
    /// Every other record (`HEADER`, `TITLE`, `SEQRES`, arbitrary garbage...)
    /// falls into the parser's catch-all arm and is silently ignored, so any
    /// text at all used to "parse" into an empty, valid molecule (#268).
    #[error("No atoms in PDB")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum MmcifError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom_site row: {0}")]
    InvalidAtomRow(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    /// No `_atom_site.*` `loop_` was ever seen at all, so zero atoms were
    /// read.
    ///
    /// A `loop_` genuinely tagged `_atom_site.*` with zero data rows is legal
    /// mmCIF (an intentionally empty structure) and is not this — only the
    /// absence of any `_atom_site.*` loop makes this unreadable garbage
    /// rather than a real, empty structure (#268).
    #[error("No atoms in mmCIF")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum CifCoreError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom_site row: {0}")]
    InvalidAtomRow(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    /// No `_atom_site_label` `loop_` was ever seen at all, so zero atoms were
    /// read -- the same distinction #268 draws for mmCIF between "genuinely
    /// no atom-site loop" and "an atom-site loop stating zero rows."
    #[error("No atoms in CIF core")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PsfError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomRow(String),

    /// PSF states no element directly -- only a name/type label (ambiguous:
    /// `CA` is alpha-carbon in every protein PSF, not calcium) and a mass,
    /// which this crate infers the element from instead. This is what a
    /// mass matching no real element within tolerance reports (#321).
    #[error("No element matches mass {0}")]
    InvalidElement(f64),

    /// No `!NATOM` section was ever seen at all.
    #[error("No atoms in PSF")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum TopError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    /// GROMACS states no element directly -- only a force-field-specific
    /// type string (`opls_135`, ambiguous by design) -- and this reader does
    /// not resolve `[ atomtypes ]` defaults, so mass must be stated inline.
    /// This is what a mass matching no real element within tolerance
    /// reports (#323), the same fallback PSF's own `element_from_mass`
    /// provides.
    #[error("No element matches mass {0}")]
    InvalidElement(f64),

    /// No `[ moleculetype ]` was ever seen at all.
    #[error("No atoms in TOP")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum LammpsError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    /// The `Atoms` section carries no recognized `# style` comment, no
    /// explicit `LammpsReadOptions::atom_style` was given, and this many
    /// columns matches more than one supported atom style (#324) --
    /// `charge` and `molecular`/`bond`/`angle` share a column count both
    /// with and without the optional image-flag triplet. A guess here is
    /// how a charge column becomes a molecule id.
    #[error(
        "{0} columns in the Atoms section is ambiguous between atom styles; specify LammpsReadOptions::atom_style"
    )]
    AmbiguousAtomStyle(usize),

    /// A recognized `# style` comment names a real LAMMPS atom style this
    /// reader does not model (`sphere`, `ellipsoid`, `electron`, ...).
    #[error("Unsupported LAMMPS atom style: {0}")]
    UnsupportedAtomStyle(String),

    /// No `Atoms` section was ever seen, or it stated zero atoms.
    #[error("No atoms in LAMMPS data file")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PrmtopError {
    #[error("Parse error: {0}")]
    ParseError(String),

    /// PRMTOP's `ATOMIC_NUMBER` section is optional (AmberTools 12+); when
    /// absent, or when it names nothing in range, the element is inferred
    /// from mass instead, the same fallback PSF's own `element_from_mass`
    /// provides (#322). This is what neither route resolving reports.
    #[error("No element matches atomic number {atomic_number:?} or mass {mass}")]
    InvalidElement {
        atomic_number: Option<i64>,
        mass: f64,
    },

    /// No `%FLAG ATOM_NAME` section was ever seen at all.
    #[error("No atoms in PRMTOP")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum Mol2Error {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    #[error("Invalid bond line: {0}")]
    InvalidBondLine(String),

    #[error("Unsupported SYBYL atom type: {0}")]
    UnsupportedAtomType(String),

    /// No `@<TRIPOS>ATOM` section was ever seen, so zero atoms were read.
    ///
    /// Any text without recognized `@<TRIPOS>` section headers used to
    /// "parse" into an empty, valid molecule instead of being rejected (#268).
    #[error("No atoms in Mol2")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum PdbqtError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    #[error("Invalid AutoDock atom type: {0}")]
    InvalidAtomType(String),

    #[error("Invalid torsion tree line: {0}")]
    InvalidTorsionTree(String),

    /// No `ATOM`/`HETATM` record was ever recognized, so zero atoms were
    /// read.
    ///
    /// Every other record (`ROOT`, `ENDROOT`, `TORSDOF`, arbitrary
    /// garbage...) falls into the parser's catch-all arm and is silently
    /// ignored, so any text at all used to "parse" into an empty, valid
    /// molecule (#268).
    #[error("No atoms in PDBQT")]
    NoAtoms,
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum GroError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom line: {0}")]
    InvalidAtomLine(String),

    #[error("Invalid element for atom name: {0}")]
    InvalidElement(String),

    #[error("Declared {expected} atoms but the file ended early")]
    AtomCountMismatch { expected: usize },
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum CommonchemError {
    #[error("Invalid JSON: {0}")]
    Json(#[from] serde_json::Error),

    #[error("No commonchem or rdkitjson version header")]
    MissingHeader,

    #[error("Unsupported {key} version {version}")]
    UnsupportedVersion { key: String, version: u32 },

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atomic number: {0}")]
    InvalidElement(u8),

    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(u32),

    #[error("Invalid atom stereo: {0}")]
    InvalidAtomStereo(String),

    #[error("Invalid bond stereo: {0}")]
    InvalidBondStereo(String),

    #[error("Bond references atom {atom}, but the molecule has {num_atoms}")]
    BondIndexOutOfRange { atom: usize, num_atoms: usize },

    #[error("{what} references index {index}, which does not exist")]
    ExtensionIndexOutOfRange { what: &'static str, index: usize },

    #[error("A dim-{dim} conformer needs {expected} coordinates, got {got}")]
    ConformerLengthMismatch {
        dim: u8,
        expected: usize,
        got: usize,
    },

    #[error("Unsupported conformer dimensionality: {0}")]
    UnsupportedConformerDim(u8),
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum CmlError {
    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid atom element: {0}")]
    InvalidAtomElement(String),

    #[error("Invalid bond element: {0}")]
    InvalidBondElement(String),

    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    #[error("Bond references unknown atom id: {0}")]
    UnknownAtomReference(String),
}

/// A record failed to read while streaming through a [`crate::io::supplier::Supplier`].
///
/// Unlike [`crate::io::reader::Skipped`] (used by the one-shot `read()`,
/// where a bad record is data to report and move on from), a `Supplier`
/// surfaces this as its iterator's `Err` — a genuinely new failure mode
/// streaming introduces that one-shot reading never had: the underlying
/// `Read` itself can fail (a broken pipe, a permissions error mid-file),
/// not just a malformed record.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum ReadError {
    #[error("I/O error at record {position}: {source}")]
    Io {
        position: usize,
        #[source]
        source: std::io::Error,
    },

    #[error("record {position} failed to parse: {message}")]
    Parse { position: usize, message: String },
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum BcifError {
    #[error("Parse error: {0}")]
    ParseError(String),

    /// Malformed MessagePack itself -- a truncated buffer, a length prefix
    /// that runs past the end of the input, and so on.
    #[error("Invalid MessagePack: {0}")]
    InvalidMessagePack(String),

    /// A MessagePack type tag this crate does not decode: ext types,
    /// timestamps, str64/array64/map64. None of these appear in a real
    /// BinaryCIF file (#319) -- a hard error here means a producer this
    /// crate has never seen, not a silently wrong read.
    #[error("Unsupported MessagePack tag: 0x{0:02x}")]
    UnsupportedTag(u8),

    /// An encoding step's `"kind"` string is not one of the seven BinaryCIF
    /// defines.
    #[error("Unknown BinaryCIF encoding kind: {0}")]
    UnknownEncoding(String),

    #[error(transparent)]
    Mmcif(#[from] MmcifError),
}

/// What `crate::io::xdr`'s `XdrReader`/`XdrWriter` return internally --
/// shared by every format built on top of XDR framing (TRR, XTC), so a
/// truncated-input failure while parsing one doesn't surface as another
/// format's own error type (#325).
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum XdrError {
    #[error("Parse error: {0}")]
    ParseError(String),
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum TrrError {
    #[error("Parse error: {0}")]
    ParseError(String),

    /// GROMACS's XDR trajectory magic number is a fixed `1993` -- anything
    /// else at byte 0 is not a TRR file.
    #[error("Not a TRR file: expected magic number 1993, got {0}")]
    InvalidMagicNumber(i32),

    /// A header's array size resolved to neither 4 nor 8 bytes per real
    /// number, the two precisions the classic `xdrfile` format supports.
    #[error("Unsupported TRR precision")]
    UnsupportedPrecision,

    #[error(transparent)]
    Xdr(#[from] XdrError),

    #[error(transparent)]
    Trajectory(#[from] crate::core::trajectory::TrajectoryError),
}

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum XtcError {
    #[error("Parse error: {0}")]
    ParseError(String),

    /// GROMACS's XDR compressed-trajectory magic number is a fixed `1995`
    /// -- anything else at byte 0 is not an XTC file.
    #[error("Not an XTC file: expected magic number 1995, got {0}")]
    InvalidMagicNumber(i32),

    #[error(transparent)]
    Xdr(#[from] XdrError),

    #[error(transparent)]
    Trajectory(#[from] crate::core::trajectory::TrajectoryError),
}
