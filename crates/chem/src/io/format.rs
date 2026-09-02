//! The format registry: what this crate can read and write, and how to ask.

use std::fmt;
use std::ops::{BitAnd, BitOr, Not};

use crate::core::atom::Chirality;
use crate::core::bond::{BondOrder, BondStereo};
use crate::core::molecule::Molecule;
use crate::io::reader::ReadOutcome;

/// What a format can carry across a conversion.
///
/// A flag set rather than a list of booleans because the interesting operation
/// is set difference: what a molecule holds, minus what the target can keep, is
/// exactly what a conversion is about to throw away.
///
/// # This describes our writer, not the specification
///
/// A descriptor declares what this crate round-trips **today**, which is often
/// less than the file format permits. MDL molfiles have columns for charges and
/// isotopes; [`crate::io::sdf::write_sdf`] documents that it does not write
/// them, so [`Format::SDF`] does not claim them. Declaring the spec's
/// capability instead would make the drop report lie about the data it exists
/// to protect — reporting nothing lost while the writer quietly dropped it,
/// which is the exact failure mode this type was added to prevent.
///
/// A format story widens its own mask when it implements an attribute, and
/// `test_declared_masks_match_what_actually_survives` refuses to let the two
/// drift apart.
///
/// Hand-rolled rather than `bitflags` deliberately: the crate has no such
/// dependency, the lean build's tree is a pinned property, and the surface
/// needed here is four operators and an iterator.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Carries(u32);

impl Carries {
    /// Atoms, elements and the bonds between them — the one thing every
    /// molecular format holds, and the reason an empty mask is a table typo.
    pub const TOPOLOGY: Carries = Carries(1 << 0);
    pub const COORDS_2D: Carries = Carries(1 << 1);
    pub const COORDS_3D: Carries = Carries(1 << 2);
    pub const FORMAL_CHARGE: Carries = Carries(1 << 3);
    pub const PARTIAL_CHARGE: Carries = Carries(1 << 4);
    pub const ISOTOPE: Carries = Carries(1 << 5);
    pub const STEREO_ATOM: Carries = Carries(1 << 6);
    pub const STEREO_BOND: Carries = Carries(1 << 7);
    pub const AROMATICITY: Carries = Carries(1 << 8);
    pub const RESIDUES: Carries = Carries(1 << 9);
    pub const B_FACTOR: Carries = Carries(1 << 10);
    pub const OCCUPANCY: Carries = Carries(1 << 11);
    pub const UNIT_CELL: Carries = Carries(1 << 12);
    pub const PROPERTIES: Carries = Carries(1 << 13);
    /// Query features — atom and bond expressions rather than concrete atoms.
    ///
    /// Declared so the numbering is stable, and deliberately never set by
    /// [`held`]: the data model has nothing to detect until query molecules
    /// arrive in the reaction and 2D-drawing wave.
    pub const QUERY: Carries = Carries(1 << 14);

    /// Every flag above, in the order the report prints them.
    const ALL: &'static [(Carries, &'static str)] = &[
        (Carries::TOPOLOGY, "topology"),
        (Carries::COORDS_2D, "coords_2d"),
        (Carries::COORDS_3D, "coords_3d"),
        (Carries::FORMAL_CHARGE, "formal_charge"),
        (Carries::PARTIAL_CHARGE, "partial_charge"),
        (Carries::ISOTOPE, "isotope"),
        (Carries::STEREO_ATOM, "stereo_atom"),
        (Carries::STEREO_BOND, "stereo_bond"),
        (Carries::AROMATICITY, "aromaticity"),
        (Carries::RESIDUES, "residues"),
        (Carries::B_FACTOR, "b_factor"),
        (Carries::OCCUPANCY, "occupancy"),
        (Carries::UNIT_CELL, "unit_cell"),
        (Carries::PROPERTIES, "properties"),
        (Carries::QUERY, "query"),
    ];

    pub const fn empty() -> Carries {
        Carries(0)
    }

    pub const fn is_empty(&self) -> bool {
        self.0 == 0
    }

    /// True when every flag in `other` is also set here.
    ///
    /// Vacuously true for an empty `other`, which is what lets a command that
    /// produced nothing in particular accept any format.
    pub const fn contains(&self, other: Carries) -> bool {
        self.0 & other.0 == other.0
    }

    /// What is set here and not in `other` — a conversion's losses.
    pub const fn difference(&self, other: Carries) -> Carries {
        Carries(self.0 & !other.0)
    }

    pub const fn intersection(&self, other: Carries) -> Carries {
        Carries(self.0 & other.0)
    }

    /// Union, usable in a `const`.
    ///
    /// `BitOr` is not a const trait, and the descriptor table is a `static`, so
    /// the masks there are built with this rather than `|`.
    pub const fn or(self, other: Carries) -> Carries {
        Carries(self.0 | other.0)
    }

    /// The flags set here, lowest bit first, as the report names them.
    pub fn names(&self) -> impl Iterator<Item = &'static str> + '_ {
        Carries::ALL
            .iter()
            .filter(move |(flag, _)| self.contains(*flag))
            .map(|(_, name)| *name)
    }
}

impl BitOr for Carries {
    type Output = Carries;
    fn bitor(self, rhs: Carries) -> Carries {
        Carries(self.0 | rhs.0)
    }
}

impl BitAnd for Carries {
    type Output = Carries;
    fn bitand(self, rhs: Carries) -> Carries {
        Carries(self.0 & rhs.0)
    }
}

impl Not for Carries {
    type Output = Carries;
    fn not(self) -> Carries {
        Carries(!self.0)
    }
}

impl fmt::Debug for Carries {
    /// The names, not the bits — a failing assertion should say `coords_3d`
    /// rather than leave the reader decoding a hex mask.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            return write!(f, "Carries()");
        }
        write!(f, "Carries({})", self.names().collect::<Vec<_>>().join("|"))
    }
}

/// What this molecule actually holds.
///
/// A free function rather than a method on [`Molecule`] so `core` does not gain
/// a dependency on `io` for a fact that only the registry cares about.
///
/// Every predicate here already existed; nothing about a molecule had to change
/// for this. [`crate::core::site::AtomSite`] keeps its fields as separate
/// `Option`s, so a b-factor is detectable without assuming that a molecule with
/// sites has all of them.
pub fn held(molecule: &Molecule) -> Carries {
    let mut carries = Carries::empty();

    if molecule.num_atoms() > 0 {
        carries = carries | Carries::TOPOLOGY;
    }
    if molecule.has_coords() {
        carries = carries | Carries::COORDS_2D;
    }
    if molecule.has_coords3() {
        carries = carries | Carries::COORDS_3D;
    }
    if molecule.has_topology() {
        carries = carries | Carries::RESIDUES;
    }
    if molecule.has_cell() {
        carries = carries | Carries::UNIT_CELL;
    }
    if !molecule.properties().is_empty() {
        carries = carries | Carries::PROPERTIES;
    }

    for atom in molecule.atoms() {
        if atom.formal_charge() != 0 {
            carries = carries | Carries::FORMAL_CHARGE;
        }
        if atom.isotope().is_some() {
            carries = carries | Carries::ISOTOPE;
        }
        if atom.is_aromatic() {
            carries = carries | Carries::AROMATICITY;
        }
        if atom.chirality() != Chirality::None {
            carries = carries | Carries::STEREO_ATOM;
        }
    }

    for bond in molecule.bonds() {
        if bond.stereo() != BondStereo::None {
            carries = carries | Carries::STEREO_BOND;
        }
        // Three channels carry the same chemical fact, and they do not travel
        // together: an SDF round trip loses the atom flag and the bond flag but
        // keeps `BondOrder::Aromatic`, because type 4 is how a molfile says it.
        // Reading only the flags would report aromaticity as lost while it sat
        // in the file.
        if bond.is_aromatic() || bond.order() == BondOrder::Aromatic {
            carries = carries | Carries::AROMATICITY;
        }
    }

    if let Some(sites) = molecule.sites() {
        for site in sites {
            if site.partial_charge.is_some() {
                carries = carries | Carries::PARTIAL_CHARGE;
            }
            if site.occupancy.is_some() {
                carries = carries | Carries::OCCUPANCY;
            }
            if site.b_factor.is_some() {
                carries = carries | Carries::B_FACTOR;
            }
        }
    }

    carries
}

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
    /// What survives a write in this format, as this crate implements it
    /// today — see [`Carries`] on why that is not the same as what the
    /// specification allows.
    pub carries: Carries,

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
        // Bracket atoms carry charge, isotope and the chirality marker;
        // lowercase carries aromaticity; `/` and `\` carry double-bond
        // geometry. Stereo arrived with #191 — this mask claimed no stereo
        // until then, and `test_declared_masks_match_what_actually_survives`
        // is what refused to let it stay that way.
        carries: Carries::TOPOLOGY
            .or(Carries::AROMATICITY)
            .or(Carries::FORMAL_CHARGE)
            .or(Carries::ISOTOPE)
            .or(Carries::STEREO_ATOM)
            .or(Carries::STEREO_BOND),
        reader: Some(crate::io::reader::read_smiles),
        writer: Some(write_smiles_records),
    },
    FormatDescriptor {
        name: "MDL MOL format",
        codes: &["sdf", "sd", "mol", "mdl"],
        extensions: &["sdf"],
        category: Category::CommonCheminformatics,
        // An atom block holds one set of positions, and the program line's
        // dimensional code says which. Charges, isotopes, chirality and bond
        // stereo are all listed under "What is not written" in
        // `sdf::write_sdf`; aromaticity survives as bond type 4. The data
        // fields the reader parses are not written back, so no PROPERTIES.
        // Widened in #197, when the V2000 writer stopped emitting the smallest
        // record that would parse. Charges and isotopes ride on `M  CHG` and
        // `M  ISO`, chirality on the atom parity column and a wedge bond,
        // data fields on the block after `M  END`.
        //
        // Every flag here was added only after
        // `test_declared_masks_match_what_actually_survives` demanded it. The
        // mask describes what a round trip *does*, not what the specification
        // permits, so it is written last and by observation.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_2D)
            .or(Carries::COORDS_3D)
            .or(Carries::AROMATICITY)
            .or(Carries::FORMAL_CHARGE)
            .or(Carries::ISOTOPE)
            .or(Carries::STEREO_ATOM)
            .or(Carries::PROPERTIES),
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

    /// What survives a write in this format — see [`Carries`].
    pub fn carries(&self) -> Carries {
        self.descriptor().carries
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

    /// One minimal molecule per attribute, each holding that attribute and as
    /// little else as possible.
    ///
    /// Per attribute rather than one maximal molecule, because attributes
    /// interfere: an SDF atom block holds one set of positions, so a molecule
    /// carrying both a layout and a conformer loses the layout — which says
    /// nothing about whether SDF can carry 2D coordinates, only that it cannot
    /// carry both at once. Isolating them makes a failure name the attribute
    /// whose claim is wrong.
    fn one_per_attribute() -> Vec<(Carries, Molecule)> {
        use crate::core::cell::UnitCell;
        use crate::core::geometry::{Point2, Point3};
        use crate::core::residue::{Chain, Residue};
        use crate::core::site::AtomSite;
        use crate::io::smiles::parse_smiles;

        let ethane = || parse_smiles("CC").expect("valid SMILES");

        let with_2d = {
            let mut m = ethane();
            m.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
                .expect("one per atom");
            m
        };
        let with_3d = {
            let mut m = ethane();
            m.set_coords3(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.5, 0.0, 0.9)])
                .expect("one per atom");
            m
        };
        let with_sites = |mutate: fn(&mut AtomSite)| {
            let mut m = ethane();
            let mut site = AtomSite::empty();
            mutate(&mut site);
            m.set_sites(vec![site, AtomSite::empty()])
                .expect("one per atom");
            m
        };
        let with_cell = {
            let mut m = ethane();
            m.set_cell(UnitCell::cubic(10.0)).expect("valid cell");
            m
        };
        let with_properties = {
            let mut m = ethane();
            m.set_property("comment".to_string(), "kept".to_string());
            m
        };
        let with_residues = {
            let mut m = ethane();
            m.set_topology(
                vec![Chain {
                    id: "A".to_string(),
                    label_id: None,
                    residues: 0..1,
                }],
                vec![Residue {
                    name: "LIG".to_string(),
                    sequence: 1,
                    insertion_code: None,
                    label_seq: None,
                    chain_ix: 0,
                    is_hetero: false,
                    atoms: 0..2,
                }],
            )
            .expect("valid topology");
            m
        };
        // Parsed rather than hand-built, since #191. The hand-built versions
        // were quietly nonsense: the bond fixture put `BondStereo::E` on
        // ethane's single C-C bond, which no format can express because the
        // configuration of a single bond means nothing. It passed only while
        // no writer emitted stereo at all, and the moment one did it reported
        // SMILES as losing something that was never expressible.
        let with_stereo_atom = parse_smiles("N[C@@H](C)C(=O)O").expect("valid SMILES");
        let with_stereo_bond = parse_smiles("F/C=C/F").expect("valid SMILES");

        vec![
            (Carries::TOPOLOGY, ethane()),
            (Carries::COORDS_2D, with_2d),
            (Carries::COORDS_3D, with_3d),
            (
                Carries::FORMAL_CHARGE,
                parse_smiles("C[NH3+]").expect("valid SMILES"),
            ),
            (
                Carries::ISOTOPE,
                parse_smiles("[13CH4]").expect("valid SMILES"),
            ),
            (
                Carries::AROMATICITY,
                parse_smiles("c1ccccc1").expect("valid SMILES"),
            ),
            (Carries::STEREO_ATOM, with_stereo_atom),
            (Carries::STEREO_BOND, with_stereo_bond),
            (
                Carries::PARTIAL_CHARGE,
                with_sites(|s| s.partial_charge = Some(-0.35)),
            ),
            (Carries::OCCUPANCY, with_sites(|s| s.occupancy = Some(0.5))),
            (Carries::B_FACTOR, with_sites(|s| s.b_factor = Some(42.0))),
            (Carries::UNIT_CELL, with_cell),
            (Carries::PROPERTIES, with_properties),
            (Carries::RESIDUES, with_residues),
        ]
    }

    #[test]
    fn test_every_flag_has_a_fixture_that_actually_holds_it() {
        // Guards the mask test from the other side. A fixture that failed to
        // set its attribute would let every format pass by never being asked.
        for (flag, molecule) in one_per_attribute() {
            assert!(
                held(&molecule).contains(flag),
                "the {flag:?} fixture does not hold it; held {:?}",
                held(&molecule)
            );
        }
        // Nothing in the data model expresses a query yet, so there is
        // deliberately no fixture and `held` must never set it.
        assert!(!held(&one_per_attribute()[0].1).contains(Carries::QUERY));
    }

    #[test]
    fn test_declared_masks_match_what_actually_survives() {
        // The masks were derived by reading the writers, which is exactly the
        // kind of claim that rots. So they are checked against reality: a mask
        // that overstates makes the drop report lie about the very data it
        // exists to protect.
        for format in all() {
            if !format.can_write() || !format.can_read() {
                continue;
            }
            for (flag, molecule) in one_per_attribute() {
                let records = vec![("probe".to_string(), molecule)];
                let text = format.write(&records).expect("can_write said so");
                let outcome = crate::io::reader::read(&text, format);
                let Some(back) = outcome.records.first() else {
                    panic!("{format:?} wrote nothing readable for {flag:?}");
                };
                let survived = held(&back.molecule).contains(flag);
                let claimed = format.carries().contains(flag);

                assert_eq!(
                    claimed,
                    survived,
                    "{} claims {flag:?}={claimed} but a round trip gives {survived}",
                    format.name()
                );
            }
        }
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
            // An empty mask means the entry was added without one, which would
            // silently report every molecule as losing everything.
            assert!(
                format.carries().contains(Carries::TOPOLOGY),
                "{format:?} declares no topology, so its mask is missing"
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
