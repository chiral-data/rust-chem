//! The format registry: what this crate can read and write, and how to ask.

use std::fmt;
use std::io::{BufRead, Write};
use std::ops::{BitAnd, BitOr, Not};

use crate::core::atom::Chirality;
use crate::core::bond::{BondOrder, BondStereo};
use crate::core::molecule::Molecule;
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::ReadOutcome;
use crate::io::supplier::{Supplier, Writer};

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
    /// Enhanced stereochemistry groups
    /// ([`crate::core::stereo_group::StereoGroup`]) — a relative, rather
    /// than absolute, configuration asserted across a set of stereocentres.
    /// CXSMILES is the only format that carries these today (#221).
    pub const STEREO_GROUP: Carries = Carries(1 << 15);

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
        (Carries::STEREO_GROUP, "stereo_group"),
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
    if !molecule.stereo_groups().is_empty() {
        carries = carries | Carries::STEREO_GROUP;
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
pub(crate) type ReadFn = fn(&str, &ReadOptions) -> ReadOutcome;

/// Serialises named molecules into one file's worth of text.
pub(crate) type WriteFn = fn(&[(String, Molecule)], &WriteOptions) -> String;

/// Builds a streaming [`Supplier`] over a boxed reader (#213).
pub(crate) type SupplierCtor = fn(Box<dyn BufRead>, &ReadOptions) -> Box<dyn Supplier>;

/// Builds a streaming [`Writer`] over a boxed writer (#213).
pub(crate) type WriterCtor = fn(Box<dyn Write>, &WriteOptions) -> Box<dyn Writer>;

/// Everything the registry knows about one format.
///
/// The reader and writer are function pointers rather than `dyn` trait
/// objects. Both work in a `static`, but pointers keep the table a plain
/// slice and skip a vtable for a call made once per file. `supplier`/
/// `writer_stream` build a `Box<dyn Supplier>`/`Box<dyn Writer>` instead —
/// those really are trait objects, since a streaming reader keeps state
/// across records — and because none of these four fields are public,
/// widening any of their signatures is not a breaking change.
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
    pub(crate) supplier: Option<SupplierCtor>,
    pub(crate) writer_stream: Option<WriterCtor>,
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
        reader: Some(crate::io::reader::read_smiles_with_options),
        writer: Some(write_smiles_records),
        supplier: Some(smiles_supplier),
        writer_stream: Some(smiles_writer_stream),
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
            // Earned in #198. Nothing is written for it — V2000 has no field
            // for double-bond geometry — so it survives only because the layout
            // now draws cis and trans differently and `perceive_bond_stereo`
            // reads the difference back.
            .or(Carries::STEREO_BOND)
            .or(Carries::PROPERTIES),
        reader: Some(crate::io::reader::read_sdf_with_options),
        writer: Some(write_sdf_records),
        supplier: Some(sdf_supplier),
        writer_stream: Some(sdf_writer_stream),
    },
    FormatDescriptor {
        name: "CXSMILES",
        codes: &["cxsmiles"],
        extensions: &["cxsmiles"],
        // One of the formats RDKit has that OpenBabel lacks (#173's Goal
        // section) — deliberately not lumped in with plain SMILES's
        // category, since it is a different toolkit's extension of it.
        category: Category::OtherCheminformatics,
        // A strict superset of plain SMILES's own mask (#221 is layered on
        // top of it, never modifies it), plus the one thing CXSMILES adds
        // that plain SMILES cannot express: enhanced stereo groups. Every
        // other CXSMILES feature (coordinates, atom labels, atomProp,
        // aromatic-bond markers, radical/valence flags) is out of scope for
        // this crate today and so is not claimed here.
        //
        // Appended after SDF rather than placed next to plain SMILES:
        // `Format::SMILES`/`Format::SDF` are hardcoded indices into this
        // slice (`Format(0)`/`Format(1)`), so a new entry has to go at the
        // end, not wherever it reads best, or it silently renumbers them.
        carries: Carries::TOPOLOGY
            .or(Carries::AROMATICITY)
            .or(Carries::FORMAL_CHARGE)
            .or(Carries::ISOTOPE)
            .or(Carries::STEREO_ATOM)
            .or(Carries::STEREO_BOND)
            .or(Carries::STEREO_GROUP),
        reader: Some(crate::io::reader::read_cxsmiles_with_options),
        writer: Some(write_cxsmiles_records),
        supplier: Some(cxsmiles_supplier),
        writer_stream: Some(cxsmiles_writer_stream),
    },
    FormatDescriptor {
        name: "XYZ",
        codes: &["xyz"],
        extensions: &["xyz"],
        category: Category::ComputationalChemistry,
        // The smallest mask registered so far: nothing else about an atom —
        // charge, isotope, stereo, aromaticity, residue membership — has a
        // column in this format at all. Deliberately no UNIT_CELL: extended
        // XYZ's `Lattice=` names Cartesian basis vectors in whatever
        // orientation the file's writer used, not necessarily this crate's
        // fixed orientation convention (`core/cell.rs`), and reading it
        // correctly would mean rotating every atom's coordinate to match —
        // not attempted here (#222), see `io/xyz.rs`'s module doc.
        carries: Carries::TOPOLOGY.or(Carries::COORDS_3D),
        reader: Some(crate::io::reader::read_xyz_with_options),
        writer: Some(write_xyz_records),
        supplier: Some(xyz_supplier),
        writer_stream: Some(xyz_writer_stream),
    },
    FormatDescriptor {
        name: "PDB",
        codes: &["pdb", "ent"],
        extensions: &["pdb", "ent"],
        category: Category::BiologicalData,
        // No FORMAL_CHARGE/ISOTOPE/STEREO_ATOM/STEREO_BOND/AROMATICITY/
        // PARTIAL_CHARGE/PROPERTIES: none of those have a real, reliably
        // populated column in this format as this crate reads it. No
        // synthesised bonds beyond explicit CONECT records either — real
        // files leave standard polymer backbone bonding to residue-template
        // inference this crate does not implement (#223), see
        // `io/pdb.rs`'s module doc for why.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::RESIDUES)
            .or(Carries::OCCUPANCY)
            .or(Carries::B_FACTOR)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_pdb_with_options),
        writer: Some(write_pdb_records),
        supplier: Some(pdb_supplier),
        writer_stream: Some(pdb_writer_stream),
    },
    FormatDescriptor {
        name: "mmCIF",
        codes: &["mmcif", "cif"],
        extensions: &["cif", "mmcif"],
        category: Category::BiologicalData,
        // Same mask as PDB, and for the same reasons -- see `io/pdb.rs`'s
        // module doc. Stricter on bonds: not even CONECT's rough
        // equivalent (`_struct_conn`) is read here (#224), so no bonds are
        // claimed or produced at all.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::RESIDUES)
            .or(Carries::OCCUPANCY)
            .or(Carries::B_FACTOR)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_mmcif_with_options),
        writer: Some(write_mmcif_records),
        supplier: Some(mmcif_supplier),
        writer_stream: Some(mmcif_writer_stream),
    },
    FormatDescriptor {
        name: "Mol2",
        codes: &["mol2"],
        extensions: &["mol2"],
        // Docking-prep, not general comp-chem, per #173's own framing.
        category: Category::MolecularDynamicsAndDocking,
        // No FORMAL_CHARGE/ISOTOPE/STEREO_ATOM/STEREO_BOND/PROPERTIES --
        // none of those have a column in this format. Mirrors SDF's own
        // COORDS_2D|COORDS_3D pair: Mol2 has no dimensionality header
        // either, so `io/mol2.rs` uses the same all-zero-z heuristic
        // `parse_sdf` already does. The first format to claim
        // PARTIAL_CHARGE -- `AtomSite::partial_charge`'s own doc comment
        // already named Mol2 as the format this was modelled for (#225).
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_2D)
            .or(Carries::COORDS_3D)
            .or(Carries::PARTIAL_CHARGE)
            .or(Carries::AROMATICITY)
            .or(Carries::RESIDUES)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_mol2_with_options),
        writer: Some(write_mol2_records),
        supplier: Some(mol2_supplier),
        writer_stream: Some(mol2_writer_stream),
    },
    FormatDescriptor {
        name: "PDBQT",
        codes: &["pdbqt"],
        extensions: &["pdbqt"],
        category: Category::MolecularDynamicsAndDocking,
        // PDB's own mask plus PARTIAL_CHARGE (the new AutoDock charge
        // column) and AROMATICITY (the `A` atom type). No UNIT_CELL --
        // real ligand PDBQT carries no cell; AutoGrid's cell lives in a
        // separate grid-parameter file, not the ligand PDBQT itself.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::RESIDUES)
            .or(Carries::OCCUPANCY)
            .or(Carries::B_FACTOR)
            .or(Carries::PARTIAL_CHARGE)
            .or(Carries::AROMATICITY),
        reader: Some(crate::io::reader::read_pdbqt_with_options),
        writer: Some(write_pdbqt_records),
        supplier: Some(pdbqt_supplier),
        writer_stream: Some(pdbqt_writer_stream),
    },
];

fn smiles_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::SmilesSupplier::new(reader, options))
}

fn cxsmiles_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::CxSmilesSupplier::new(reader, options))
}

fn sdf_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::SdfSupplier::new(reader, options))
}

fn xyz_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::XyzSupplier::new(reader, options))
}

fn pdb_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::PdbSupplier::new(reader, options))
}

fn mmcif_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::MmcifSupplier::new(reader, options))
}

fn mol2_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::Mol2Supplier::new(reader, options))
}

fn pdbqt_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::PdbqtSupplier::new(reader, options))
}

fn smiles_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::SmilesWriter::new(writer, options))
}

fn cxsmiles_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::CxSmilesWriter::new(writer, options))
}

fn sdf_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::SdfWriter::new(writer, options))
}

fn xyz_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::XyzWriter::new(writer, options))
}

fn pdb_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::PdbWriter::new(writer, options))
}

fn mmcif_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::MmcifWriter::new(writer, options))
}

fn mol2_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::Mol2Writer::new(writer, options))
}

fn pdbqt_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::PdbqtWriter::new(writer, options))
}

fn write_smiles_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // SMILES has no write options today.
    let mut out = String::new();
    for (name, molecule) in records {
        // `name<TAB>` is not the convention: the reader splits on whitespace
        // and takes everything after the first token as the name, so a space
        // is what round-trips.
        //
        // Canonical (#220): the same molecule, built with atoms in a
        // different order, always writes the same string.
        out.push_str(&crate::io::smiles_writer::write_smiles_for_molecule_canonical(molecule));
        out.push(' ');
        out.push_str(name);
        out.push('\n');
    }
    out
}

fn write_cxsmiles_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // CXSMILES has no write options of its own today, beyond what plain
    // SMILES already has none of.
    let mut out = String::new();
    for (name, molecule) in records {
        out.push_str(&crate::io::cxsmiles::write_cxsmiles(molecule));
        out.push(' ');
        out.push_str(name);
        out.push('\n');
    }
    out
}

fn write_sdf_records(records: &[(String, Molecule)], options: &WriteOptions) -> String {
    let mut molecules = Vec::with_capacity(records.len());
    for (name, molecule) in records {
        let mut copy = molecule.clone();
        copy.set_name(name.clone());
        molecules.push(copy);
    }
    crate::io::sdf::write_sdf_all_with_options(&molecules, &options.sdf)
}

fn write_xyz_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // XYZ has no write options today.
    let mut out = String::new();
    for (name, molecule) in records {
        let mut copy = molecule.clone();
        copy.set_name(name.clone());
        out.push_str(&crate::io::xyz::write_xyz(&copy));
    }
    out
}

fn write_pdb_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // PDB has no write options today, and no per-record name to thread
    // through -- see `PdbWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::pdb::write_pdb(molecule));
    }
    out
}

fn write_mmcif_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // mmCIF has no write options today, and no per-record name to thread
    // through -- see `MmcifWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::mmcif::write_mmcif(molecule));
    }
    out
}

fn write_mol2_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // Mol2 has no write options today.
    let mut out = String::new();
    for (name, molecule) in records {
        let mut copy = molecule.clone();
        copy.set_name(name.clone());
        out.push_str(&crate::io::mol2::write_mol2(&copy));
    }
    out
}

fn write_pdbqt_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // PDBQT has no write options today, and no per-record name to thread
    // through -- see `PdbqtWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::pdbqt::write_pdbqt(molecule));
    }
    out
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
    /// CXSMILES (#221) — enhanced stereo groups only, see
    /// [`crate::io::cxsmiles`].
    pub const CXSMILES: Format = Format(2);
    /// XYZ (#222) — coordinates only, no unit cell, see [`crate::io::xyz`].
    pub const XYZ: Format = Format(3);
    /// PDB (#223) — no synthesised bonds beyond `CONECT`, see
    /// [`crate::io::pdb`].
    pub const PDB: Format = Format(4);
    /// mmCIF (#224) — no bonds at all, see [`crate::io::mmcif`].
    pub const MMCIF: Format = Format(5);
    /// Mol2 (#225) — a curated SYBYL type subset, see [`crate::io::mol2`].
    pub const MOL2: Format = Format(6);
    /// PDBQT (#226) — a new rotatable-bond/fragment-decomposition
    /// algorithm on write, see [`crate::io::pdbqt`].
    pub const PDBQT: Format = Format(7);

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

    /// Serialises molecules in this format with default options, carrying
    /// their names, or `None` if the format cannot be written.
    ///
    /// The public way to reach a writer. The function pointer itself stays
    /// private — the binary is a separate crate, so `pub(crate)` would not
    /// reach it anyway.
    pub fn write(&self, records: &[(String, Molecule)]) -> Option<String> {
        self.write_with_options(records, &WriteOptions::default())
    }

    /// [`Self::write`], with explicit per-format options (#212) — e.g. which
    /// molfile dialect to write for SDF.
    pub fn write_with_options(
        &self,
        records: &[(String, Molecule)],
        options: &WriteOptions,
    ) -> Option<String> {
        self.descriptor()
            .writer
            .map(|writer| writer(records, options))
    }

    /// Streams molecules from `reader` one at a time, rather than
    /// materializing the whole file first (#213), or `None` if the format
    /// cannot be read.
    pub fn supplier(
        &self,
        reader: impl BufRead + 'static,
        options: &ReadOptions,
    ) -> Option<Box<dyn Supplier>> {
        self.descriptor()
            .supplier
            .map(|ctor| ctor(Box::new(reader), options))
    }

    /// Streams molecules to `writer` one at a time, rather than
    /// materializing one `String` for the whole output first (#213), or
    /// `None` if the format cannot be written.
    pub fn writer_stream(
        &self,
        writer: impl Write + 'static,
        options: &WriteOptions,
    ) -> Option<Box<dyn Writer>> {
        self.descriptor()
            .writer_stream
            .map(|ctor| ctor(Box::new(writer), options))
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
        // #223: this used to assert `None`, "not registered yet" -- the day
        // one is, per that comment's own point, is exactly this one.
        for code in ["pdb", "ent"] {
            assert_eq!(Format::from_code(code), Some(Format::PDB), "{code}");
        }
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
        // The case this test used to pin: `.pdb` fell back to SMILES because
        // no PDB format was registered. #223 registered one, so this changed
        // -- exactly the point the old comment here was making.
        assert_eq!(Format::from_filename("a.pdb"), Format::PDB);
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
        use crate::core::stereo_group::{StereoGroup, StereoGroupKind};
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
        let with_stereo_group = {
            let mut m = parse_smiles("N[C@@H](C)C(=O)O").expect("valid SMILES");
            m.set_stereo_groups(vec![StereoGroup::new(
                StereoGroupKind::And,
                Some(1),
                vec![1],
            )])
            .expect("valid stereo groups");
            m
        };

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
            (Carries::STEREO_GROUP, with_stereo_group),
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
    fn test_write_with_options_reaches_the_sdf_writer() {
        // #212: the registry's WriteFn now threads WriteOptions through to
        // the actual writer -- default options must produce the same output
        // as the no-options convenience.
        use crate::core::prelude::*;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        let records = vec![("m".to_string(), mol)];

        assert_eq!(
            Format::SDF.write(&records),
            Format::SDF.write_with_options(&records, &WriteOptions::default())
        );
    }

    #[test]
    fn test_every_registered_format_reads_and_writes() {
        // Named for what it checks rather than a count that was only ever
        // true while there happened to be exactly two: every format the
        // registry has grown since (#221's CXSMILES included) has to keep
        // satisfying this, not just the first two.
        assert_eq!(all().count(), 8);
        for format in all() {
            assert!(format.can_read() && format.can_write(), "{format:?}");
        }
    }
}
