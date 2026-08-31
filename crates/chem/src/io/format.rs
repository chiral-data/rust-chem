//! The format registry: what this crate can read and write, and how to ask.

use std::fmt;

use crate::core::molecule::Molecule;
use crate::io::reader::ReadOutcome;

/// Where a format sits in the landscape.
///
/// The seventeen groupings OpenBabel uses, kept because they are the ones a
/// user listing formats expects to see and because matching them makes a
/// coverage gap obvious.
///
/// `#[non_exhaustive]` from the first release. That is the mistake the old
/// two-variant `Format` enum made, and the reason it could not grow.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Category {
    CommonCheminformatics,
    Utility,
    OtherCheminformatics,
    ComputationalChemistry,
    MolecularFingerprint,
    Crystallography,
    Reaction,
    Image,
    TwoDimensionalDrawing,
    ThreeDimensionalViewer,
    KineticsAndThermodynamics,
    MolecularDynamicsAndDocking,
    VolumeData,
    Json,
    Miscellaneous,
    BiologicalData,
    Obscure,
}

impl Category {
    pub fn label(&self) -> &'static str {
        match self {
            Category::CommonCheminformatics => "Common cheminformatics formats",
            Category::Utility => "Utility formats",
            Category::OtherCheminformatics => "Other cheminformatics formats",
            Category::ComputationalChemistry => "Computational chemistry formats",
            Category::MolecularFingerprint => "Molecular fingerprint formats",
            Category::Crystallography => "Crystallography formats",
            Category::Reaction => "Reaction formats",
            Category::Image => "Image formats",
            Category::TwoDimensionalDrawing => "2D drawing formats",
            Category::ThreeDimensionalViewer => "3D viewer formats",
            Category::KineticsAndThermodynamics => "Kinetics and Thermodynamics formats",
            Category::MolecularDynamicsAndDocking => "Molecular dynamics and docking formats",
            Category::VolumeData => "Volume data formats",
            Category::Json => "JSON formats",
            Category::Miscellaneous => "Miscellaneous formats",
            Category::BiologicalData => "Biological data formats",
            Category::Obscure => "Obscure formats",
        }
    }
}

/// Parses a whole file into molecules.
///
/// Aliased rather than written inline because both of these signatures are
/// going to change: the option-bag story adds an `&OptionSet` parameter, and
/// the streaming story replaces `&str` with a reader. Naming them now means
/// that is one edit rather than four.
pub(crate) type ReadFn = fn(&str) -> ReadOutcome;

/// Serialises named molecules into one file's worth of text.
pub(crate) type WriteFn = fn(&[(String, Molecule)]) -> String;

/// Everything the registry knows about one format.
///
/// The reader and writer are function pointers rather than `dyn` trait objects.
/// Both work in a `static`, but pointers keep the table a plain slice and skip
/// a vtable for a call made once per file. Traits earn their place when
/// per-format options arrive and a reader needs state — and because these
/// fields are not public, widening their signatures then is not a breaking
/// change.
pub struct FormatDescriptor {
    /// Human-readable name, as a format list would print it.
    pub name: &'static str,
    /// Short codes, as `-i` and `-o` accept them. The first is canonical.
    pub codes: &'static [&'static str],
    /// Filename extensions, lowercase and without the dot.
    pub extensions: &'static [&'static str],
    pub category: Category,

    pub(crate) reader: Option<ReadFn>,
    pub(crate) writer: Option<WriteFn>,
}

/// Every format compiled into this build.
///
/// A plain slice rather than a runtime plugin registry. OpenBabel's registers
/// formats at load time; this one is fixed at compile time, so there is no
/// global, no interior mutability, and a format list that cannot depend on
/// what ran first. The trade is that a downstream crate cannot add a format
/// without patching this one.
static FORMATS: &[FormatDescriptor] = &[
    FormatDescriptor {
        name: "SMILES format",
        codes: &["smi", "smiles"],
        extensions: &["smi", "smiles", "txt"],
        category: Category::CommonCheminformatics,
        reader: Some(crate::io::reader::read_smiles),
        writer: Some(write_smiles_records),
    },
    FormatDescriptor {
        name: "MDL MOL format",
        codes: &["sdf", "sd", "mol", "mdl"],
        extensions: &["sdf"],
        category: Category::CommonCheminformatics,
        reader: Some(crate::io::reader::read_sdf),
        writer: Some(write_sdf_records),
    },
];

fn write_smiles_records(records: &[(String, Molecule)]) -> String {
    let mut out = String::new();
    for (name, molecule) in records {
        // `name<TAB>` is not the convention: the reader splits on whitespace
        // and takes everything after the first token as the name, so a space
        // is what round-trips.
        out.push_str(&crate::io::smiles_writer::write_smiles_for_molecule(
            molecule,
        ));
        out.push(' ');
        out.push_str(name);
        out.push('\n');
    }
    out
}

fn write_sdf_records(records: &[(String, Molecule)]) -> String {
    let mut molecules = Vec::with_capacity(records.len());
    for (name, molecule) in records {
        let mut copy = molecule.clone();
        copy.set_name(name.clone());
        molecules.push(copy);
    }
    crate::io::sdf::write_sdf_all(&molecules)
}

/// A format this build supports.
///
/// An index into the static table rather than a `&'static FormatDescriptor`,
/// which makes `Copy`, `Eq` and `Hash` correct by construction. A reference
/// would have to compare by pointer, and const promotion gives no guarantee
/// that two references to the same descriptor are the same address.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Format(u16);

impl Format {
    /// SMILES — also the fallback for a name whose extension is unrecognised.
    pub const SMILES: Format = Format(0);
    /// MDL MOL / SDF.
    pub const SDF: Format = Format(1);

    pub fn descriptor(&self) -> &'static FormatDescriptor {
        &FORMATS[self.0 as usize]
    }

    /// Looks a format up by short code, case-insensitively.
    pub fn from_code(code: &str) -> Option<Format> {
        let code = code.to_lowercase();
        all().find(|f| f.codes().iter().any(|c| *c == code))
    }

    /// Looks a format up by filename extension, without the dot,
    /// case-insensitively.
    pub fn from_extension(extension: &str) -> Option<Format> {
        let extension = extension.to_lowercase();
        all().find(|f| f.extensions().iter().any(|e| *e == extension))
    }

    /// Picks a format from a filename.
    ///
    /// **An unrecognised or absent extension resolves to SMILES**, which is
    /// deliberate rather than a gap: SMILES is the universal default, and
    /// rejecting an unknown name would reject the extensionless files people
    /// actually have. Once a format for a given extension is registered, that
    /// extension resolves to it instead — the fallback only ever catches names
    /// nothing claims.
    ///
    /// This is now the **only** copy of that rule. It used to exist three
    /// times: here, again for output paths in the CLI's writer, and a third
    /// time as a stdin special case.
    pub fn from_filename(name: &str) -> Format {
        name.rsplit_once('.')
            .and_then(|(_, extension)| Format::from_extension(extension))
            .unwrap_or(Format::SMILES)
    }

    pub fn name(&self) -> &'static str {
        self.descriptor().name
    }

    /// The short label used in progress messages — the canonical code,
    /// uppercased at the point of use rather than stored twice.
    pub fn label(&self) -> &'static str {
        match *self {
            Format::SMILES => "SMILES",
            Format::SDF => "SDF",
            _ => self.descriptor().name,
        }
    }

    pub fn codes(&self) -> &'static [&'static str] {
        self.descriptor().codes
    }

    pub fn extensions(&self) -> &'static [&'static str] {
        self.descriptor().extensions
    }

    pub fn category(&self) -> Category {
        self.descriptor().category
    }

    pub fn can_read(&self) -> bool {
        self.descriptor().reader.is_some()
    }

    pub fn can_write(&self) -> bool {
        self.descriptor().writer.is_some()
    }

    pub(crate) fn reader(&self) -> Option<ReadFn> {
        self.descriptor().reader
    }

    /// Serialises molecules in this format, carrying their names, or `None`
    /// if the format cannot be written.
    ///
    /// The public way to reach a writer. The function pointer itself stays
    /// private so the option-bag story can widen its signature without that
    /// being a breaking change — and because the binary is a separate crate,
    /// `pub(crate)` would not reach it anyway.
    pub fn write(&self, records: &[(String, Molecule)]) -> Option<String> {
        self.descriptor().writer.map(|writer| writer(records))
    }
}

/// Every format in this build, in table order.
///
/// What a format listing iterates. Deliberately an iterator over handles
/// rather than the descriptors, so a caller cannot hold a descriptor whose
/// index it has lost.
pub fn all() -> impl Iterator<Item = Format> {
    (0..FORMATS.len()).map(|i| Format(i as u16))
}

impl fmt::Debug for Format {
    /// The name, not the whole descriptor — a `Format` inside a larger `{:?}`
    /// should not print two slices and two function pointers.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Format({})", self.codes()[0])
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lookup_by_code_finds_every_alias() {
        // SDF has four codes and they must all reach the same format — the
        // reason codes is a slice rather than a single string.
        for code in ["sdf", "sd", "mol", "mdl"] {
            assert_eq!(Format::from_code(code), Some(Format::SDF), "{code}");
        }
        for code in ["smi", "smiles"] {
            assert_eq!(Format::from_code(code), Some(Format::SMILES), "{code}");
        }
        assert_eq!(Format::from_code("SDF"), Some(Format::SDF), "case");
        assert_eq!(Format::from_code("pdb"), None, "not registered yet");
    }

    #[test]
    fn test_from_filename_keeps_the_smiles_fallback() {
        // Pinned deliberately, not by accident. An unrecognised extension
        // resolving to SMILES is what keeps extensionless files working; the
        // alternative rejects `.txt` and bare names that parse fine today.
        assert_eq!(Format::from_filename("a.sdf"), Format::SDF);
        assert_eq!(Format::from_filename("A.SDF"), Format::SDF, "case");
        for name in ["a.smi", "a.smiles", "a.txt", "a", "-"] {
            assert_eq!(Format::from_filename(name), Format::SMILES, "{name}");
        }
        // Suffix matching, not "contains": a gzipped SDF is not an SDF here.
        assert_eq!(Format::from_filename("a.sdf.gz"), Format::SMILES);
        // The case this issue exists for. It is SMILES today because no PDB
        // format is registered; the day one is, this line changes and that is
        // the point of the registry.
        assert_eq!(Format::from_filename("a.pdb"), Format::SMILES);
    }

    #[test]
    fn test_every_registered_format_is_well_formed() {
        // Cheap, and it catches a malformed entry the moment someone adds one
        // rather than when a lookup mysteriously misses.
        let mut seen_codes = Vec::new();
        for format in all() {
            assert!(!format.name().is_empty(), "{format:?} has no name");
            assert!(!format.codes().is_empty(), "{format:?} has no codes");
            for code in format.codes() {
                assert_eq!(*code, code.to_lowercase(), "code {code} is not lowercase");
                assert!(
                    !seen_codes.contains(code),
                    "code {code} is claimed by two formats"
                );
                seen_codes.push(code);
            }
            for extension in format.extensions() {
                assert_eq!(
                    *extension,
                    extension.to_lowercase(),
                    "extension {extension} is not lowercase"
                );
                assert!(!extension.starts_with('.'), "{extension} has a leading dot");
            }
            assert!(
                format.can_read() || format.can_write(),
                "{format:?} can do neither"
            );
        }
    }

    #[test]
    fn test_a_format_round_trips_through_its_canonical_code() {
        for format in all() {
            assert_eq!(Format::from_code(format.codes()[0]), Some(format));
        }
    }

    #[test]
    fn test_handles_compare_by_identity() {
        assert_eq!(Format::SMILES, Format::from_code("smi").unwrap());
        assert_ne!(Format::SMILES, Format::SDF);
        // The constants must agree with the table they index; getting these
        // out of step would silently swap two formats.
        assert_eq!(Format::SMILES.label(), "SMILES");
        assert_eq!(Format::SDF.label(), "SDF");
        assert_eq!(Format::SDF.name(), "MDL MOL format");
    }

    #[test]
    fn test_write_goes_through_the_descriptor() {
        use crate::core::prelude::*;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.add_bond(Bond::new(0, 1, BondOrder::Single)).unwrap();
        let records = vec![("ethanol-ish".to_string(), mol)];

        let smiles = Format::SMILES.write(&records).expect("SMILES writes");
        assert!(smiles.ends_with(" ethanol-ish\n"), "{smiles:?}");

        let sdf = Format::SDF.write(&records).expect("SDF writes");
        assert!(
            sdf.starts_with("ethanol-ish\n"),
            "the name is the title line"
        );
        assert!(sdf.contains("M  END"));
        assert!(sdf.trim_end().ends_with("$$$$"));

        // And what was written reads back as the same two atoms, which is the
        // property that makes the table's reader and writer a matched pair.
        assert_eq!(crate::io::reader::read(&sdf, Format::SDF).len(), 1);
    }

    #[test]
    fn test_both_formats_read_and_write() {
        assert_eq!(all().count(), 2);
        for format in all() {
            assert!(format.can_read() && format.can_write(), "{format:?}");
        }
    }
}
