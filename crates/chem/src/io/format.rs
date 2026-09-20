//! The format registry: what this crate can read and write, and how to ask.

use std::fmt;
use std::io::{BufRead, Write};
use std::ops::{BitAnd, BitOr, Not};

use crate::core::atom::Chirality;
use crate::core::bond::{BondOrder, BondStereo};
use crate::core::mesh::Mesh;
use crate::core::molecule::Molecule;
use crate::core::table::Table;
use crate::core::trajectory::{Trajectory, TrajectoryError};
use crate::core::volume::VolumeGrid;
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
    /// Atoms and their elements — the one thing every molecular format holds,
    /// and the reason an empty mask is a table typo.
    ///
    /// Bonds are [`Carries::BONDS`] and not this flag. They used to be: this
    /// doc comment read "atoms, elements and the bonds between them" while four
    /// formats claiming it read back none at all.
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
    /// The bonds between atoms, separate from [`Carries::TOPOLOGY`] because
    /// four registered formats carry atoms and no bonds at all — XYZ and GRO
    /// have no bond block, mmCIF's `_struct_conn` is unread (#224), and PDBQT
    /// gives only `BRANCH` pivots (#226).
    ///
    /// Split out in #257. Before it, a conversion into any of the four turned
    /// benzene into six unbonded carbons and the drop report named only the
    /// coordinates: [`held`] set `TOPOLOGY` on the atom count alone, so the
    /// claim was never testable. Appended rather than inserted so the existing
    /// bit numbering stays put; the report order is a separate table.
    pub const BONDS: Carries = Carries(1 << 16);
    /// A force-field atom type name (`"CT"`, `"OW"`), from
    /// [`crate::core::force_field::ForceFieldAtom::atom_type`]. PSF and TOP
    /// state one per atom; PRMTOP does too, in `AMBER_ATOM_TYPE`.
    pub const ATOM_TYPE: Carries = Carries(1 << 17);
    /// A force-field atomic mass, from
    /// [`crate::core::force_field::ForceFieldAtom::mass`]. All three of PSF,
    /// TOP and PRMTOP state one per atom.
    pub const MASS: Carries = Carries(1 << 18);
    /// Three-atom angle terms, from
    /// [`crate::core::force_field::ForceFieldTopology::angles`].
    pub const ANGLES: Carries = Carries(1 << 19);
    /// Four-atom proper-torsion terms, from
    /// [`crate::core::force_field::ForceFieldTopology::dihedrals`].
    pub const DIHEDRALS: Carries = Carries(1 << 20);
    /// Four-atom improper (out-of-plane) torsion terms, from
    /// [`crate::core::force_field::ForceFieldTopology::impropers`].
    pub const IMPROPERS: Carries = Carries(1 << 21);
    /// Nonbonded-exclusion atom pairs, from
    /// [`crate::core::force_field::ForceFieldTopology::exclusions`].
    pub const EXCLUSIONS: Carries = Carries(1 << 22);
    /// Per-atom velocities, from [`crate::core::trajectory::Frame::velocities`].
    /// TRR carries these; XTC never does.
    pub const VELOCITIES: Carries = Carries(1 << 23);
    /// Per-atom forces, from [`crate::core::trajectory::Frame::forces`].
    pub const FORCES: Carries = Carries(1 << 24);
    /// A frame's simulation time, from
    /// [`crate::core::trajectory::Frame::time`]. Named `FRAME_TIME` rather
    /// than `TIME` to stay unambiguous next to
    /// [`Category::KineticsAndThermodynamics`], a category nothing in this
    /// type otherwise names.
    pub const FRAME_TIME: Carries = Carries(1 << 25);
    /// A grid's sampled scalar values, from
    /// [`crate::core::volume::VolumeGrid::values`] — the one thing every
    /// volumetric format holds, the same role [`Carries::TOPOLOGY`] plays
    /// for a molecule.
    pub const SAMPLES: Carries = Carries(1 << 26);
    /// A mesh's vertex positions, from [`crate::core::mesh::Mesh::vertices`]
    /// — the one thing every mesh format holds; see [`Carries::FACES`] for
    /// why the faces are a separate flag.
    pub const VERTICES: Carries = Carries(1 << 27);
    /// A mesh's faces, separate from [`Carries::VERTICES`] the same reason
    /// [`Carries::BONDS`] is separate from [`Carries::TOPOLOGY`]: a format
    /// could in principle carry vertex positions with no face list at all
    /// (a point cloud).
    pub const FACES: Carries = Carries(1 << 28);
    /// A table's columns, from [`crate::core::table::Table::columns`] — the
    /// one thing every tabular format holds.
    pub const COLUMNS: Carries = Carries(1 << 29);
    /// Hydrogen-bond donor pairs, from
    /// [`crate::core::force_field::ForceFieldTopology::donors`]. PSF is the
    /// first format to state these (#321).
    pub const DONORS: Carries = Carries(1 << 30);
    /// Hydrogen-bond acceptor pairs, from
    /// [`crate::core::force_field::ForceFieldTopology::acceptors`]. PSF is
    /// the first format to state these (#321). The last bit this `u32` has
    /// room for — a 33rd flag needs a wider representation, not this one's
    /// problem.
    pub const ACCEPTORS: Carries = Carries(1 << 31);

    /// Every flag above, in the order the report prints them.
    const ALL: &'static [(Carries, &'static str)] = &[
        (Carries::TOPOLOGY, "topology"),
        (Carries::BONDS, "bonds"),
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
        // Force-field attributes (#315), grouped together here even though
        // their bits were appended at the end of the numeric list above.
        (Carries::ATOM_TYPE, "atom_type"),
        (Carries::MASS, "mass"),
        (Carries::ANGLES, "angles"),
        (Carries::DIHEDRALS, "dihedrals"),
        (Carries::IMPROPERS, "impropers"),
        (Carries::EXCLUSIONS, "exclusions"),
        (Carries::DONORS, "donors"),
        (Carries::ACCEPTORS, "acceptors"),
        // Non-molecule kinds (#316): a trajectory's frames, a volume grid's
        // samples, a mesh's surface, a table's columns.
        (Carries::VELOCITIES, "velocities"),
        (Carries::FORCES, "forces"),
        (Carries::FRAME_TIME, "frame_time"),
        (Carries::SAMPLES, "samples"),
        (Carries::VERTICES, "vertices"),
        (Carries::FACES, "faces"),
        (Carries::COLUMNS, "columns"),
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
    if molecule.num_bonds() > 0 {
        carries = carries | Carries::BONDS;
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
        // Three channels carry the same chemical fact. Since #261 a reader
        // reconciles them at the boundary, so a molecule this crate produced
        // sets all three -- but a mask is a claim about a *format*, and CML's
        // only aromatic channel is the bond order (`order="A"`). Reading just
        // the flags would report aromaticity as lost while it sat in the file,
        // which is what this OR is for. (It used to cite SDF, whose reader has
        // perceived since #197.)
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

    if let Some(force_field) = molecule.force_field() {
        if let Some(atoms) = &force_field.atoms {
            for atom in atoms {
                if atom.atom_type.is_some() {
                    carries = carries | Carries::ATOM_TYPE;
                }
                if atom.mass.is_some() {
                    carries = carries | Carries::MASS;
                }
                // Ors into the same flag the `sites` loop above already
                // sets -- a Mol2 charge and an Amber charge are different
                // provenances for the same fact, and the mask only claims
                // that a partial charge survives, not which table it lives
                // in.
                if atom.partial_charge.is_some() {
                    carries = carries | Carries::PARTIAL_CHARGE;
                }
            }
        }
        if !force_field.angles.is_empty() {
            carries = carries | Carries::ANGLES;
        }
        if !force_field.dihedrals.is_empty() {
            carries = carries | Carries::DIHEDRALS;
        }
        if !force_field.impropers.is_empty() {
            carries = carries | Carries::IMPROPERS;
        }
        if !force_field.exclusions.is_empty() {
            carries = carries | Carries::EXCLUSIONS;
        }
        if !force_field.donors.is_empty() {
            carries = carries | Carries::DONORS;
        }
        if !force_field.acceptors.is_empty() {
            carries = carries | Carries::ACCEPTORS;
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
    /// A triangulated surface, no chemistry at all -- OBJ (#335), PLY
    /// (#336).
    MeshData,
    /// Named columns and rows, no structure implied -- CSV (#337).
    TabularData,
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
            Category::MeshData => "Mesh data formats",
            Category::TabularData => "Tabular data formats",
            Category::Json => "JSON formats",
            Category::Miscellaneous => "Miscellaneous formats",
            Category::BiologicalData => "Biological data formats",
            Category::Obscure => "Obscure formats",
        }
    }
}

/// Whether a format's canonical bytes are UTF-8 text or an arbitrary binary
/// layout.
///
/// Introduced by #309 so the registry can widen to binary formats (XTC, DCD,
/// CCP4, ...) without repointing the eleven existing text `reader`/`writer`
/// function pointers, which stay exactly as typed and shaped as they are
/// today.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    Text,
    Binary,
}

/// What a format's records *are*.
///
/// `test_every_registered_format_is_well_formed`'s `Carries::TOPOLOGY`
/// assertion used to be unconditional -- a fair typo-catcher while every
/// registered format was a molecule format. #310 is the first format-shaped
/// story where that stops being true: a density map, a mesh and a table have
/// no atoms to declare, so the invariant has to know which formats are
/// making a claim about atoms at all before it can enforce one.
///
/// `#[non_exhaustive]`, like [`Category`]: a variant added later (this enum
/// names all five kinds the v0.9.0 milestone needs, and as of #314 every one
/// has a container type -- but no format is registered as anything but
/// `Molecules` yet) forces every exhaustive match inside this crate to be
/// revisited rather than silently compiling with a wrong assumption.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// One conformer, one topology -- what every format registered before
    /// #310 already is.
    Molecules,
    /// One topology, many conformers -- a trajectory
    /// ([`crate::core::trajectory::Trajectory`], #311). No format is
    /// registered as this yet; the first will be #329 (LAMMPS trajectory) or
    /// #326 (TRR).
    Frames,
    /// A scalar field on a grid -- a density map
    /// ([`crate::core::volume::VolumeGrid`], #312). No format is registered
    /// as this yet; the first will be #331 (CUBE).
    Volume,
    /// Vertices, normals, faces -- no chemistry at all
    /// ([`crate::core::mesh::Mesh`], #313). No format is registered as this
    /// yet; the first will be #335 (OBJ) or #336 (PLY).
    Mesh,
    /// Typed columns, no structure implied ([`crate::core::table::Table`],
    /// #314). No format is registered as this yet; the first will be #337
    /// (CSV).
    Table,
}

/// Parses a whole file into molecules.
pub(crate) type ReadFn = fn(&str, &ReadOptions) -> ReadOutcome;

/// Serialises named molecules into one file's worth of text.
pub(crate) type WriteFn = fn(&[(String, Molecule)], &WriteOptions) -> String;

/// Parses a whole file into molecules, from raw bytes rather than decoded
/// text — what a binary format's reader is shaped like (#309).
pub(crate) type ByteReadFn = fn(&[u8], &ReadOptions) -> ReadOutcome;

/// Serialises named molecules into one file's worth of bytes — what a binary
/// format's writer is shaped like (#309).
pub(crate) type ByteWriteFn = fn(&[(String, Molecule)], &WriteOptions) -> Vec<u8>;

/// Serialises a whole trajectory into one file's worth of bytes (#326) —
/// singular, unlike [`ByteWriteFn`]'s `&[(String, Molecule)]` list, since one
/// trajectory file holds exactly one trajectory, not a list of named
/// records. A new, parallel field rather than a widened `WriteFn`/
/// `ByteWriteFn`: those are fundamentally shaped for `Molecule`, and
/// generalising them to be payload-polymorphic is separate, larger work
/// this format does not need (see [`crate::io::trr`]'s module doc).
pub(crate) type ByteWriteTrajectoryFn = fn(&mut Trajectory, &WriteOptions) -> Vec<u8>;

/// Serialises a whole [`VolumeGrid`] into one file's worth of bytes (#331) —
/// the same "new, parallel field" reasoning [`ByteWriteTrajectoryFn`] already
/// documents, since a grid is neither a `Molecule` list nor a `Trajectory`.
/// Byte-returning like that field too, so a later binary volumetric format
/// (CCP4/DSN6, #332-#333) reuses this same field rather than needing a
/// second one just because CUBE happened to be text.
pub(crate) type ByteWriteVolumeFn = fn(&VolumeGrid, &WriteOptions) -> Vec<u8>;

/// Serialises a whole [`Mesh`] into one file's worth of bytes (#335) — the
/// same "new, parallel field" reasoning [`ByteWriteVolumeFn`] already
/// documents, since a mesh is neither a `Molecule` list nor a `VolumeGrid`.
pub(crate) type ByteWriteMeshFn = fn(&Mesh, &WriteOptions) -> Vec<u8>;

/// Serialises a whole [`Table`] into one file's worth of bytes (#337) — the
/// same "new, parallel field" reasoning [`ByteWriteMeshFn`] already
/// documents, since a table is neither a `Molecule` list nor any of the
/// other three payload shapes.
pub(crate) type ByteWriteTableFn = fn(&Table, &WriteOptions) -> Vec<u8>;

/// A byte pattern identifying a format's content, independent of any
/// filename: the exact bytes expected starting at `offset` (#317).
///
/// CCP4's `MAP ` sits at byte 208, not byte 0 — the offset is part of the
/// signature, not always zero the way gzip's `\x1f\x8b` is.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Signature {
    pub offset: usize,
    pub bytes: &'static [u8],
}

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
    /// Whether this format's canonical bytes are UTF-8 text or binary.
    pub encoding: Encoding,
    /// What this format's records *are* -- see [`Kind`].
    pub kind: Kind,
    /// Content signatures identifying this format independent of its
    /// filename — checked by `crate::io::open::open_supplier` only when
    /// nothing claims the file's extension, before falling back to SMILES
    /// (#317). Empty for every format registered today; a binary format
    /// story (#319, #325-#337) is what will populate this.
    pub(crate) magic: &'static [Signature],

    pub(crate) reader: Option<ReadFn>,
    pub(crate) writer: Option<WriteFn>,
    /// Set only for a format whose canonical reader takes raw bytes (#309) —
    /// `None` for every text format, including all eleven registered today.
    pub(crate) reader_bytes: Option<ByteReadFn>,
    /// Set only for a format whose canonical writer produces raw bytes
    /// (#309) — `None` for every text format, including all eleven
    /// registered today.
    pub(crate) writer_bytes: Option<ByteWriteFn>,
    pub(crate) supplier: Option<SupplierCtor>,
    pub(crate) writer_stream: Option<WriterCtor>,
    /// Set only for a format whose writer takes a whole [`Trajectory`]
    /// rather than a `&[(String, Molecule)]` list (#326) — `None` for every
    /// `Kind::Molecules` format, including all seventeen registered before
    /// TRR.
    pub(crate) writer_trajectory: Option<ByteWriteTrajectoryFn>,
    /// Set only for a format whose writer takes a whole [`VolumeGrid`]
    /// (#331) — `None` for every format registered before CUBE.
    pub(crate) writer_volume: Option<ByteWriteVolumeFn>,
    /// Set only for a format whose writer takes a whole [`Mesh`] (#335) —
    /// `None` for every format registered before OBJ.
    pub(crate) writer_mesh: Option<ByteWriteMeshFn>,
    /// Set only for a format whose writer takes a whole [`Table`] (#337) —
    /// `None` for every format registered before CSV. Unlike
    /// `writer_volume`/`writer_mesh`, this can be populated *alongside* an
    /// ordinary `writer`/`writer_bytes` on the same descriptor: CSV writes
    /// both a `Table` (its declared kind) and, as a disclosed opt-in, a
    /// molecule list with a `smiles` column.
    pub(crate) writer_table: Option<ByteWriteTableFn>,
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
            .or(Carries::BONDS)
            .or(Carries::AROMATICITY)
            .or(Carries::FORMAL_CHARGE)
            .or(Carries::ISOTOPE)
            .or(Carries::STEREO_ATOM)
            .or(Carries::STEREO_BOND),
        reader: Some(crate::io::reader::read_smiles_with_options),
        writer: Some(write_smiles_records),
        supplier: Some(smiles_supplier),
        writer_stream: Some(smiles_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "MDL MOL format",
        codes: &["sdf", "sd", "mol", "mdl"],
        // `mol` is also a code (`-imol` already worked), but was missing
        // here -- a `.mol` file, the more common name for a single-molecule
        // molfile than `.sdf`, silently read as SMILES with every line
        // skipped (#318).
        extensions: &["sdf", "mol"],
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
            .or(Carries::BONDS)
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
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
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
            .or(Carries::BONDS)
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
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
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
        // not attempted here (#222), see `io/xyz.rs`'s module doc. No BONDS:
        // the format has no bond block, so this reads back atoms and nothing
        // joining them.
        carries: Carries::TOPOLOGY.or(Carries::COORDS_3D),
        reader: Some(crate::io::reader::read_xyz_with_options),
        writer: Some(write_xyz_records),
        supplier: Some(xyz_supplier),
        writer_stream: Some(xyz_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
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
            .or(Carries::BONDS)
            .or(Carries::COORDS_3D)
            .or(Carries::RESIDUES)
            .or(Carries::OCCUPANCY)
            .or(Carries::B_FACTOR)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_pdb_with_options),
        writer: Some(write_pdb_records),
        supplier: Some(pdb_supplier),
        writer_stream: Some(pdb_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "mmCIF",
        codes: &["mmcif", "cif"],
        extensions: &["cif", "mmcif"],
        category: Category::BiologicalData,
        // Same mask as PDB, and for the same reasons -- see `io/pdb.rs`'s
        // module doc. Stricter on bonds: not even CONECT's rough
        // equivalent (`_struct_conn`) is read here (#224), so no bonds are
        // claimed or produced at all -- which is what the absent BONDS says
        // now that #257 has split it out of TOPOLOGY.
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
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
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
            .or(Carries::BONDS)
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
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "PDBQT",
        codes: &["pdbqt"],
        extensions: &["pdbqt"],
        category: Category::MolecularDynamicsAndDocking,
        // PDB's own mask plus PARTIAL_CHARGE (the new AutoDock charge
        // column) and AROMATICITY (the `A` atom type). No UNIT_CELL --
        // real ligand PDBQT carries no cell; AutoGrid's cell lives in a
        // separate grid-parameter file, not the ligand PDBQT itself. No BONDS
        // either: reading recovers only the BRANCH pivots, not a bond graph.
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
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "GRO",
        codes: &["gro"],
        extensions: &["gro"],
        category: Category::MolecularDynamicsAndDocking,
        // No BONDS/OCCUPANCY/B_FACTOR/PARTIAL_CHARGE/FORMAL_CHARGE/ISOTOPE/
        // STEREO_ATOM/STEREO_BOND/AROMATICITY/PROPERTIES -- none of those
        // have a column in this format at all; GROMACS keeps connectivity in
        // the topology file, not here.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::RESIDUES)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_gro_with_options),
        writer: Some(write_gro_records),
        supplier: Some(gro_supplier),
        writer_stream: Some(gro_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "CML",
        codes: &["cml"],
        extensions: &["cml"],
        category: Category::CommonCheminformatics,
        // No RESIDUES/UNIT_CELL/PARTIAL_CHARGE/OCCUPANCY/B_FACTOR/
        // PROPERTIES/STEREO_ATOM/STEREO_BOND -- this crate only reads and
        // writes the core <molecule> element (#228), see `io/cml.rs`'s
        // module doc for the full scope cut. AROMATICITY survives through
        // bond order alone (`order="A"`), the same channel SDF's type-4
        // bonds use -- no atom-level aromaticity flag in this format.
        carries: Carries::TOPOLOGY
            .or(Carries::BONDS)
            .or(Carries::COORDS_2D)
            .or(Carries::COORDS_3D)
            .or(Carries::FORMAL_CHARGE)
            .or(Carries::ISOTOPE)
            .or(Carries::AROMATICITY),
        reader: Some(crate::io::reader::read_cml_with_options),
        writer: Some(write_cml_records),
        supplier: Some(cml_supplier),
        writer_stream: Some(cml_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "commonchem JSON",
        // `json` is claimed as a bare code and as the extension because this
        // is the only JSON format in the crate today. #230 (utility formats)
        // may want ChemDoodle JSON, which is a different schema on the same
        // extension -- that story reassigns these, it does not add a second
        // claimant.
        codes: &["commonchem", "cjson", "json"],
        extensions: &["json"],
        category: Category::Json,
        // No PARTIAL_CHARGE/RESIDUES/B_FACTOR/OCCUPANCY/UNIT_CELL -- the
        // schema has no field for any of them. No STEREO_GROUP: commonchem
        // does have `stereoGroups`, but mapping it onto #221's `StereoGroup`
        // is its own story. AROMATICITY rides on the `rdkitRepresentation`
        // extension rather than on a bond order, the only channel this format
        // has for it -- see `io/commonchem.rs`'s module doc.
        carries: Carries::TOPOLOGY
            .or(Carries::BONDS)
            .or(Carries::COORDS_2D)
            .or(Carries::COORDS_3D)
            .or(Carries::FORMAL_CHARGE)
            .or(Carries::ISOTOPE)
            .or(Carries::STEREO_ATOM)
            .or(Carries::STEREO_BOND)
            .or(Carries::AROMATICITY)
            .or(Carries::PROPERTIES),
        reader: Some(crate::io::reader::read_commonchem_with_options),
        writer: Some(write_commonchem_records),
        supplier: Some(commonchem_supplier),
        writer_stream: Some(commonchem_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "BinaryCIF",
        codes: &["bcif", "binarycif"],
        extensions: &["bcif"],
        category: Category::BiologicalData,
        // Same mask as mmCIF, and for the same reasons -- see `io/mmcif.rs`'s
        // module doc. Both formats go through the same
        // `cif_model::build_molecule`/`build_rows` (#319), so what survives
        // a round trip is identical by construction.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::RESIDUES)
            .or(Carries::OCCUPANCY)
            .or(Carries::B_FACTOR)
            .or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(bcif_supplier),
        writer_stream: Some(bcif_writer_stream),
        encoding: Encoding::Binary,
        kind: Kind::Molecules,
        // MessagePack has no fixed leading bytes -- its first byte encodes
        // the top-level value's own type and length, so there is nothing
        // to sniff. Resolved by extension only, same as every format
        // registered before #317 added `magic` at all.
        magic: &[],
        reader_bytes: Some(crate::io::bcif::read_bcif_with_options),
        writer_bytes: Some(crate::io::bcif::write_bcif_records),
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "CIF core",
        // Not "cif" -- mmCIF's descriptor already claims that code
        // (`codes: &["mmcif", "cif"]` above), and codes must be unique
        // (`test_every_registered_format_is_well_formed`). Extensions are
        // not required to be unique, which is exactly what lets this share
        // `.cif` with mmCIF at all -- content, not the name, decides which
        // one a `.cif` file resolves to (#320, see `cif_core::
        // is_small_molecule_cif` and its two call sites in `io::open` and
        // `bin::chem::stream`).
        codes: &["cif-core"],
        extensions: &["cif"],
        category: Category::Crystallography,
        // No BONDS, no RESIDUES: this dictionary has no bond loop and no
        // chain/residue notion at all (no `label_asym_id`/`auth_seq_id`
        // machinery) -- see `io/cif_core.rs`'s module doc.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::OCCUPANCY)
            .or(Carries::B_FACTOR)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_cif_core_with_options),
        writer: Some(write_cif_core_records),
        supplier: Some(cif_core_supplier),
        writer_stream: Some(cif_core_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        // ASCII `data_` text, same as mmCIF -- nothing fixed-offset to
        // sniff; the `.cif` dispatch special case does its own scan
        // instead of going through this mechanism, see the module doc.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "PSF",
        codes: &["psf"],
        extensions: &["psf"],
        category: Category::MolecularDynamicsAndDocking,
        // Deliberately no COORDS_2D/COORDS_3D -- PSF states no coordinates
        // at all, the same "topology-only, zero geometry" shape SMILES
        // already establishes as legitimate. The first real user of
        // `ForceFieldTopology` (#315), and of its DONORS/ACCEPTORS flags
        // (#321) -- see `io/psf.rs`'s module doc.
        carries: Carries::TOPOLOGY
            .or(Carries::BONDS)
            .or(Carries::RESIDUES)
            .or(Carries::ATOM_TYPE)
            .or(Carries::MASS)
            .or(Carries::PARTIAL_CHARGE)
            .or(Carries::ANGLES)
            .or(Carries::DIHEDRALS)
            .or(Carries::IMPROPERS)
            .or(Carries::EXCLUSIONS)
            .or(Carries::DONORS)
            .or(Carries::ACCEPTORS),
        reader: Some(crate::io::reader::read_psf_with_options),
        writer: Some(write_psf_records),
        supplier: Some(psf_supplier),
        writer_stream: Some(psf_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        // ASCII `PSF` text -- nothing fixed-offset to sniff, and no
        // extension ambiguity to resolve (`.psf` is claimed by nothing
        // else), so `magic` buys nothing here.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "PRMTOP",
        codes: &["prmtop"],
        // `"top"` included deliberately (#323): AMBER topologies are
        // sometimes saved under `.top` too, alongside the native
        // `prmtop`/`parm7`. `Format::from_filename_checked` resolves a bare
        // `.top` to whichever format is registered first -- this one, so
        // no behavior change for anyone already relying on it -- and
        // `io::open::resolve_format`/`bin/chem/stream.rs`'s own copy
        // override to `Format::TOP` by content when `top::is_gromacs_top`
        // says so, the same shape `.cif`'s mmCIF/CIF-core disambiguation
        // already uses.
        extensions: &["prmtop", "parm7", "top"],
        category: Category::MolecularDynamicsAndDocking,
        // Topology only (#322) -- PRMTOP states substantially more of the
        // force field than PSF does (force constants, equilibrium values,
        // Lennard-Jones coefficients), but those have no home in
        // `ForceFieldTopology`, which already frames them as out of scope
        // for this milestone, and `Carries` is completely full at 32/32
        // bits. Same topology mask PSF uses, minus DONORS/ACCEPTORS --
        // PRMTOP has no donor/acceptor section at all.
        carries: Carries::TOPOLOGY
            .or(Carries::BONDS)
            .or(Carries::RESIDUES)
            .or(Carries::ATOM_TYPE)
            .or(Carries::MASS)
            .or(Carries::PARTIAL_CHARGE)
            .or(Carries::ANGLES)
            .or(Carries::DIHEDRALS)
            .or(Carries::IMPROPERS)
            .or(Carries::EXCLUSIONS),
        reader: Some(crate::io::reader::read_prmtop_with_options),
        writer: Some(write_prmtop_records),
        supplier: Some(prmtop_supplier),
        writer_stream: Some(prmtop_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        // ASCII `%VERSION` text -- nothing fixed-offset to sniff, and no
        // extension ambiguity to resolve (`.prmtop`/`.parm7` are claimed by
        // nothing else), so `magic` buys nothing here.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "TOP",
        codes: &["top"],
        // Shares `.top` with PRMTOP -- see PRMTOP's own `extensions`
        // comment for the disambiguation story (#323).
        extensions: &["top"],
        category: Category::MolecularDynamicsAndDocking,
        // Topology only, same reasoning as PRMTOP: force constants,
        // equilibrium values and the `[ *types ]` tables that carry them
        // have no home in `ForceFieldTopology`, and `Carries` is
        // completely full at 32/32 bits. Same mask PSF/PRMTOP use, minus
        // DONORS/ACCEPTORS -- GROMACS has no such section either.
        carries: Carries::TOPOLOGY
            .or(Carries::BONDS)
            .or(Carries::RESIDUES)
            .or(Carries::ATOM_TYPE)
            .or(Carries::MASS)
            .or(Carries::PARTIAL_CHARGE)
            .or(Carries::ANGLES)
            .or(Carries::DIHEDRALS)
            .or(Carries::IMPROPERS)
            .or(Carries::EXCLUSIONS),
        reader: Some(crate::io::reader::read_top_with_options),
        writer: Some(write_top_records),
        supplier: Some(top_supplier),
        writer_stream: Some(top_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        // ASCII `[ section ]` text -- nothing fixed-offset to sniff, and
        // the one extension ambiguity that exists (`.top`, with PRMTOP) is
        // resolved by `is_gromacs_top` at the `io::open`/CLI entry points,
        // not through this mechanism.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "LAMMPS Data",
        codes: &["lammps"],
        // No reliable extension convention really exists for this format
        // (confirmed by research -- real files are often named
        // `data.<system>` with no suffix at all); these are a reasonable,
        // non-colliding default, not a claim of authority (#324).
        extensions: &["lammps", "lmp", "data"],
        category: Category::MolecularDynamicsAndDocking,
        // The first force-field-topology format that also states
        // coordinates -- see `io/lammps.rs`'s module doc. Topology only,
        // same reasoning as PSF/PRMTOP/TOP. No RESIDUES (a molecule-id
        // column exists but is a bare integer, no name to build a residue
        // from) and no EXCLUSIONS (LAMMPS has no such section at all --
        // exclusions are an input-script `special_bonds` setting, not
        // file data).
        carries: Carries::TOPOLOGY
            .or(Carries::BONDS)
            .or(Carries::ATOM_TYPE)
            .or(Carries::MASS)
            .or(Carries::PARTIAL_CHARGE)
            .or(Carries::ANGLES)
            .or(Carries::DIHEDRALS)
            .or(Carries::IMPROPERS)
            .or(Carries::COORDS_3D)
            .or(Carries::UNIT_CELL),
        reader: Some(crate::io::reader::read_lammps_data_with_options),
        writer: Some(write_lammps_data_records),
        supplier: Some(lammps_data_supplier),
        writer_stream: Some(lammps_data_writer_stream),
        encoding: Encoding::Text,
        kind: Kind::Molecules,
        // ASCII text with an arbitrary first line -- nothing fixed-offset
        // to sniff, and no extension ambiguity to resolve (nothing else
        // claims `lammps`/`lmp`/`data`).
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "TRR",
        codes: &["trr"],
        extensions: &["trr"],
        category: Category::MolecularDynamicsAndDocking,
        // The first `Kind::Frames` format -- topology (a bare atom count,
        // see `io/trr.rs`'s module doc), positions, and whatever a given
        // frame happened to also state: velocities, forces, simulation
        // time, a box. No RESIDUES/ATOM_TYPE/MASS/BONDS -- TRR states none
        // of a molecule's identity, only its trajectory.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::VELOCITIES)
            .or(Carries::FORCES)
            .or(Carries::FRAME_TIME)
            .or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(trr_supplier),
        writer_stream: None,
        encoding: Encoding::Binary,
        kind: Kind::Frames,
        // GROMACS's own fixed magic number, big-endian -- the first real,
        // populated signature in this crate (#317 built the mechanism;
        // #319/#325-#337 are what populate it for a real format).
        magic: &[Signature {
            offset: 0,
            bytes: &[0x00, 0x00, 0x07, 0xC9],
        }],
        reader_bytes: Some(crate::io::trr::read_trr_bytes),
        writer_bytes: None,
        writer_trajectory: Some(crate::io::trr::write_trr_bytes),
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "XTC",
        codes: &["xtc"],
        extensions: &["xtc"],
        category: Category::MolecularDynamicsAndDocking,
        // The second `Kind::Frames` format -- same shared topology
        // (`Element::UNKNOWN`) as TRR, but no VELOCITIES/FORCES at all:
        // XTC's compressed coordinate block only ever carries positions.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::FRAME_TIME)
            .or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(xtc_supplier),
        writer_stream: None,
        encoding: Encoding::Binary,
        kind: Kind::Frames,
        // GROMACS's own fixed magic number, big-endian -- distinct from
        // TRR's `1993` by exactly 2 (#325).
        magic: &[Signature {
            offset: 0,
            bytes: &[0x00, 0x00, 0x07, 0xCB],
        }],
        reader_bytes: Some(crate::io::xtc::read_xtc_bytes),
        writer_bytes: None,
        writer_trajectory: Some(crate::io::xtc::write_xtc_bytes),
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "DCD",
        codes: &["dcd"],
        extensions: &["dcd"],
        category: Category::MolecularDynamicsAndDocking,
        // The third `Kind::Frames` format -- same shared topology
        // (`Element::UNKNOWN`) as TRR/XTC, no VELOCITIES/FORCES: DCD never
        // carries either.
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::FRAME_TIME)
            .or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(dcd_supplier),
        writer_stream: None,
        encoding: Encoding::Binary,
        kind: Kind::Frames,
        // "CORD" at byte 4, not byte 0 -- the leading four bytes are a
        // Fortran record-length marker whose own byte order depends on
        // the file's endianness, but "CORD" itself is plain ASCII and
        // endianness-independent (#327).
        magic: &[Signature {
            offset: 4,
            bytes: b"CORD",
        }],
        reader_bytes: Some(crate::io::dcd::read_dcd_bytes),
        writer_bytes: None,
        writer_trajectory: Some(crate::io::dcd::write_dcd_bytes),
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "NCTRAJ",
        codes: &["nctraj"],
        extensions: &["nc"],
        category: Category::MolecularDynamicsAndDocking,
        // The fourth `Kind::Frames` format -- same shared topology
        // (`Element::UNKNOWN`) as TRR/XTC/DCD. Unlike XTC/DCD, this format
        // genuinely can carry velocities (like TRR); unlike TRR, it has no
        // per-frame time/step at all -- `time`/`forces` are both out of
        // scope for this story (#328).
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::VELOCITIES)
            .or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(nctraj_supplier),
        writer_stream: None,
        encoding: Encoding::Binary,
        kind: Kind::Frames,
        // Every NetCDF-3 file starts with "CDF" then a version byte (`1`
        // classic, `2` 64-bit offset) -- this only proves "a NetCDF-3
        // file", not specifically an Amber one; the `Conventions ==
        // "AMBER"` check happens one layer up, in `io::nctraj` itself.
        magic: &[
            Signature {
                offset: 0,
                bytes: b"CDF\x01",
            },
            Signature {
                offset: 0,
                bytes: b"CDF\x02",
            },
        ],
        reader_bytes: Some(crate::io::nctraj::read_nctraj_bytes),
        writer_bytes: None,
        writer_trajectory: Some(crate::io::nctraj::write_nctraj_bytes),
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "LAMMPS Trajectory",
        codes: &["lammpstrj"],
        extensions: &["lammpstrj", "dump"],
        category: Category::MolecularDynamicsAndDocking,
        // The fifth `Kind::Frames` format, and the first that's text: same
        // shared topology (`Element::UNKNOWN`) as TRR/XTC/DCD/NCTRAJ. The
        // `ATOMS` column list can name `vx/vy/vz`/`fx/fy/fz`, so the mask
        // states what the format is *capable* of carrying, the same
        // always-on posture TRR's own VELOCITIES|FORCES already
        // established -- not every file states them. No FRAME_TIME: a
        // dump's TIMESTEP is a step count, not a time (#329).
        carries: Carries::TOPOLOGY
            .or(Carries::COORDS_3D)
            .or(Carries::VELOCITIES)
            .or(Carries::FORCES)
            .or(Carries::UNIT_CELL),
        // Plain text, so this uses `reader` (the `&str` entry point) rather
        // than `reader_bytes` -- keeps `Format::encoding` and the presence
        // of a byte reader in agreement, the same invariant every
        // `Kind::Molecules` text format already satisfies. Writing a
        // trajectory has only one registry field at all
        // (`writer_trajectory`, always byte-returning), so the writer
        // stays there regardless of encoding.
        reader: Some(crate::io::lammpstrj::read_lammpstrj),
        writer: None,
        supplier: Some(lammpstrj_supplier),
        writer_stream: None,
        encoding: Encoding::Text,
        kind: Kind::Frames,
        // Text, no magic bytes -- resolved by extension only, the same
        // posture LAMMPS Data already takes.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: Some(crate::io::lammpstrj::write_lammpstrj_bytes),
        writer_volume: None,
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "CUBE",
        codes: &["cube"],
        extensions: &["cube", "cub"],
        category: Category::VolumeData,
        // The first `Kind::Volume` format, and the first to carry
        // `TOPOLOGY`/`COORDS_3D` alongside `SAMPLES` -- CUBE genuinely
        // states both a grid and the atoms that produced it
        // (`VolumeGrid::atoms`, #331). No UNIT_CELL: CUBE states a
        // sampling box, not a separate crystallographic cell.
        carries: Carries::SAMPLES
            .or(Carries::TOPOLOGY)
            .or(Carries::COORDS_3D),
        reader: Some(crate::io::cube::read_cube_with_options),
        writer: None,
        supplier: Some(cube_supplier),
        writer_stream: None,
        encoding: Encoding::Text,
        kind: Kind::Volume,
        // Text, no magic bytes -- two free-text comment lines up front,
        // resolved by extension only, the same posture every other text
        // format in this crate already takes.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: Some(crate::io::cube::write_cube_bytes),
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "CCP4/MRC",
        codes: &["ccp4", "mrc"],
        extensions: &["ccp4", "mrc", "map"],
        category: Category::VolumeData,
        // The second `Kind::Volume` format, and the first to carry a real
        // crystallographic UNIT_CELL alongside SAMPLES -- CCP4 states one
        // directly, unlike CUBE's atoms-shaped TOPOLOGY/COORDS_3D. No
        // atoms at all here.
        carries: Carries::SAMPLES.or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(ccp4_supplier),
        writer_stream: None,
        encoding: Encoding::Binary,
        kind: Kind::Volume,
        // "MAP " at byte 208, confirmed directly against a real file this
        // story built with `gemmi` -- proves "this is a CCP4/MRC-shaped
        // file", not that it's well-formed past that (see io::ccp4).
        magic: &[Signature {
            offset: 208,
            bytes: b"MAP ",
        }],
        reader_bytes: Some(crate::io::ccp4::read_ccp4_bytes),
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: Some(crate::io::ccp4::write_ccp4_bytes),
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "DX",
        codes: &["dx"],
        extensions: &["dx"],
        category: Category::VolumeData,
        // The third `Kind::Volume` format, and the simplest -- no atoms,
        // no axis permutation, no endianness. Deliberately no UNIT_CELL:
        // DX states no crystallographic cell at all, so a CCP4-to-DX
        // round trip loses it, and this absence is exactly how that loss
        // is declared (#333).
        carries: Carries::SAMPLES,
        // Plain text, so this uses `reader` (the `&str` entry point)
        // rather than `reader_bytes` -- keeps `Format::encoding` and the
        // presence of a byte reader in agreement, the same posture CUBE
        // and LAMMPS Trajectory already established.
        reader: Some(crate::io::dx::read_dx_with_options),
        writer: None,
        supplier: Some(dx_supplier),
        writer_stream: None,
        encoding: Encoding::Text,
        kind: Kind::Volume,
        // Text, no magic bytes -- resolved by extension only, the same
        // posture CUBE already takes.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: Some(crate::io::dx::write_dx_bytes),
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "DSN6",
        codes: &["dsn6"],
        extensions: &["dsn6", "omap"],
        category: Category::VolumeData,
        // The fourth and last `Kind::Volume` format -- always carries a
        // real crystallographic cell (like CCP4, unlike DX), but with no
        // axis permutation at all (a fixed c=X, r=Y, s=Z mapping) and no
        // atoms (#334).
        carries: Carries::SAMPLES.or(Carries::UNIT_CELL),
        reader: None,
        writer: None,
        supplier: Some(dsn6_supplier),
        writer_stream: None,
        encoding: Encoding::Binary,
        kind: Kind::Volume,
        // No fixed byte signature at all -- the closest thing DSN6 has is a
        // documented constant word, checked in io::dsn6's own reader, not a
        // byte pattern this table can express. Resolved by extension only.
        magic: &[],
        reader_bytes: Some(crate::io::dsn6::read_dsn6_bytes),
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: Some(crate::io::dsn6::write_dsn6_bytes),
        writer_mesh: None,
        writer_table: None,
    },
    FormatDescriptor {
        name: "OBJ",
        codes: &["obj"],
        extensions: &["obj"],
        category: Category::MeshData,
        // The first `Kind::Mesh` format, and the first registered format
        // with no chemistry at all. Always carries vertex positions and a
        // face list -- OBJ's own `f` lines are what a mesh format is for
        // (#335).
        carries: Carries::VERTICES.or(Carries::FACES),
        // Plain text, so this uses `reader` (the `&str` entry point)
        // rather than `reader_bytes` -- the same posture CUBE/DX already
        // established.
        reader: Some(crate::io::obj::read_obj_with_options),
        writer: None,
        supplier: Some(obj_supplier),
        writer_stream: None,
        encoding: Encoding::Text,
        kind: Kind::Mesh,
        // Text, no magic bytes -- resolved by extension only, the same
        // posture every other text format already takes.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: Some(crate::io::obj::write_obj_bytes),
        writer_table: None,
    },
    FormatDescriptor {
        name: "PLY",
        codes: &["ply"],
        extensions: &["ply"],
        category: Category::MeshData,
        // The second and last `Kind::Mesh` format -- same defining mask as
        // OBJ (#336).
        carries: Carries::VERTICES.or(Carries::FACES),
        reader: None,
        writer: None,
        supplier: Some(ply_supplier),
        writer_stream: None,
        // Binary, since two of PLY's three wire encodings (declared in its
        // own header, sniffed by io::ply's reader) are raw bytes -- the
        // ASCII variant is still read/written correctly through the same
        // byte-based entry points, the same posture CCP4/DSN6 already take.
        encoding: Encoding::Binary,
        kind: Kind::Mesh,
        // "ply\n" is byte-identical across all three encodings -- a real,
        // reliable fixed signature, unlike OBJ/DX's extension-only
        // resolution.
        magic: &[Signature {
            offset: 0,
            bytes: b"ply\n",
        }],
        reader_bytes: Some(crate::io::ply::read_ply_bytes),
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: Some(crate::io::ply::write_ply_bytes),
        writer_table: None,
    },
    FormatDescriptor {
        name: "CSV",
        codes: &["csv"],
        extensions: &["csv"],
        category: Category::TabularData,
        // The declared kind is Table -- the honest, always-safe default. A
        // structure column is a disclosed, explicit read-option opt-in
        // (io::csv::CsvReadOptions), the same "declared kind stays fixed,
        // an option flips the produced Payload" mechanism #330 already
        // proved safe for XYZ/PDB/PDBQT/GRO's own opt-in Frames reading
        // (#337).
        carries: Carries::COLUMNS,
        reader: Some(crate::io::csv::read_csv_with_options),
        // The opt-in molecule-list write shape -- a `smiles` column plus
        // one column per distinct molecule property. Populated *alongside*
        // `writer_table` below, a first for this registry: every other
        // format with a writer_volume/writer_mesh leaves `writer` `None`.
        writer: Some(crate::io::csv::write_csv_records),
        supplier: Some(csv_supplier),
        writer_stream: None,
        encoding: Encoding::Text,
        kind: Kind::Table,
        // Text, no magic bytes -- resolved by extension only, the same
        // posture every other text format already takes.
        magic: &[],
        reader_bytes: None,
        writer_bytes: None,
        writer_trajectory: None,
        writer_volume: None,
        writer_mesh: None,
        // The declared/primary write shape.
        writer_table: Some(crate::io::csv::write_csv_table_bytes),
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

fn psf_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::PsfSupplier::new(reader, options))
}

fn prmtop_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::PrmtopSupplier::new(reader, options))
}

fn lammps_data_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::LammpsDataSupplier::new(
        reader, options,
    ))
}

fn top_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::TopSupplier::new(reader, options))
}

fn cif_core_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::CifCoreSupplier::new(reader, options))
}

fn mol2_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::Mol2Supplier::new(reader, options))
}

fn pdbqt_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::PdbqtSupplier::new(reader, options))
}

fn gro_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::GroSupplier::new(reader, options))
}

fn commonchem_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::CommonchemSupplier::new(
        reader, options,
    ))
}

fn cml_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::supplier::CmlSupplier::new(reader, options))
}

fn bcif_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::bcif::BcifSupplier::new(reader, options))
}

fn trr_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::trr::TrrSupplier::new(reader, options))
}

fn xtc_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::xtc::XtcSupplier::new(reader, options))
}

fn dcd_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::dcd::DcdSupplier::new(reader, options))
}

fn nctraj_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::nctraj::NctrajSupplier::new(reader, options))
}

fn lammpstrj_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::lammpstrj::LammpstrjSupplier::new(
        reader, options,
    ))
}

fn csv_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::csv::CsvSupplier::new(reader, options))
}

fn cube_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::cube::CubeSupplier::new(reader, options))
}

fn ccp4_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::ccp4::Ccp4Supplier::new(reader, options))
}

fn dx_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::dx::DxSupplier::new(reader, options))
}

fn obj_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::obj::ObjSupplier::new(reader, options))
}

fn ply_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::ply::PlySupplier::new(reader, options))
}

fn dsn6_supplier(reader: Box<dyn BufRead>, options: &ReadOptions) -> Box<dyn Supplier> {
    Box::new(crate::io::dsn6::Dsn6Supplier::new(reader, options))
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

fn psf_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::PsfWriter::new(writer, options))
}

fn prmtop_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::PrmtopWriter::new(writer, options))
}

fn lammps_data_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::LammpsDataWriter::new(writer, options))
}

fn top_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::TopWriter::new(writer, options))
}

fn cif_core_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::CifCoreWriter::new(writer, options))
}

fn mol2_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::Mol2Writer::new(writer, options))
}

fn pdbqt_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::PdbqtWriter::new(writer, options))
}

fn gro_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::GroWriter::new(writer, options))
}

fn commonchem_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::CommonchemWriter::new(writer, options))
}

fn cml_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::supplier::CmlWriter::new(writer, options))
}

fn bcif_writer_stream(writer: Box<dyn Write>, options: &WriteOptions) -> Box<dyn Writer> {
    Box::new(crate::io::bcif::BcifWriter::new(writer, options))
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
    //
    // Framed rather than concatenated: `ENDMDL` is the record boundary both
    // readers split on, and without it N structures read back as one merged
    // molecule (#267). One structure is returned unframed, as it always was.
    let structures: Vec<String> = records
        .iter()
        .map(|(_, molecule)| crate::io::pdb::write_pdb(molecule))
        .collect();
    crate::io::pdb::frame_models(&structures)
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

fn write_psf_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // No write options today, and no per-record name to thread through --
    // see `PsfWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::psf::write_psf(molecule));
    }
    out
}

fn write_prmtop_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // No write options today, and no per-record name to thread through --
    // see `PrmtopWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::prmtop::write_prmtop(molecule));
    }
    out
}

fn write_lammps_data_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // No write options today, and no per-record name to thread through --
    // see `LammpsDataWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::lammps::write_lammps_data(molecule));
    }
    out
}

fn write_top_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // Unlike PSF/PRMTOP (looping and concatenating), this one takes all
    // records at once, the same reason `write_commonchem_records` does --
    // see `TopWriter`'s own doc comment: the `[ system ]`/`[ molecules ]`
    // footer can only be written once every record has arrived.
    crate::io::top::write_top(records)
}

fn write_cif_core_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // No write options today, and no per-record name to thread through --
    // see `CifCoreWriter`'s own doc comment.
    let mut out = String::new();
    for (_, molecule) in records {
        out.push_str(&crate::io::cif_core::write_cif_core(molecule));
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
    //
    // Framed like PDB's, which is where PDBQT takes this convention from --
    // and where AutoDock Vina's own multi-pose output takes it too (#267).
    let ligands: Vec<String> = records
        .iter()
        .map(|(_, molecule)| crate::io::pdbqt::write_pdbqt(molecule))
        .collect();
    crate::io::pdb::frame_models(&ligands)
}

fn write_gro_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // GRO has no write options today.
    let mut out = String::new();
    for (name, molecule) in records {
        let mut copy = molecule.clone();
        copy.set_name(name.clone());
        out.push_str(&crate::io::gro::write_gro(&copy));
    }
    out
}

fn write_cml_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // CML has no write options today. Always wrapped in a single `<cml>`
    // root, regardless of record count: several sibling `<molecule>`
    // elements with no enclosing root is not valid XML (a strict parser
    // stops after the first) -- confirmed against OpenBabel's own CML
    // reader during this story's verification, and matches what RDKit's
    // own CML writer does even for a single molecule.
    let mut out = String::from("<cml xmlns=\"http://www.xml-cml.org/schema\">\n");
    for (name, molecule) in records {
        let mut copy = molecule.clone();
        copy.set_name(name.clone());
        out.push_str(&crate::io::cml::write_cml(&copy));
    }
    out.push_str("</cml>\n");
    out
}

fn write_commonchem_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    // commonchem has no write options today. Unlike every other writer here
    // this one takes all the records at once rather than looping and
    // concatenating: the output is a single JSON document with one
    // `molecules` array, and two documents back to back are not a document.
    crate::io::commonchem::write_commonchem(records)
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
    /// GRO (#227) — the one format needing a real nm/Å unit conversion,
    /// see [`crate::io::gro`].
    pub const GRO: Format = Format(8);
    /// CML (#228) — the `<molecule>` element only, see [`crate::io::cml`].
    /// The first format needing an XML parser; deliberately not gated
    /// behind a cargo feature — see `Cargo.toml`'s `roxmltree` dependency
    /// comment and #184.
    pub const CML: Format = Format(9);
    /// commonchem JSON (#229) — the first structured-interchange format, and
    /// the first where aromaticity survives a round trip in both directions,
    /// see [`crate::io::commonchem`]. Not gated behind a cargo feature, on
    /// the ruling #184 closed with — see `Cargo.toml`'s `serde_json`
    /// dependency comment.
    pub const COMMONCHEM: Format = Format(10);
    /// BinaryCIF (#319) — MessagePack-encoded mmCIF, the same
    /// `cif_model::build_molecule`/`build_rows` mmCIF text goes through,
    /// see [`crate::io::bcif`]. The first format to use `Encoding::Binary`,
    /// `reader_bytes`/`writer_bytes` (#309's plumbing, unused until now),
    /// and a hand-rolled codec ([`crate::io::msgpack`]) rather than a
    /// dependency — the crate's existing style for a small, bounded
    /// surface, the same reasoning `Carries` is hand-rolled instead of
    /// depending on `bitflags`.
    pub const BCIF: Format = Format(11);
    /// CIF core (#320) — the small-molecule crystallography dictionary,
    /// old-style flat tags and fractional coordinates rather than mmCIF's
    /// dot-namespace and Cartesian ones, see [`crate::io::cif_core`]. Reads
    /// and writes only the asymmetric unit a file states, no symmetry
    /// expansion. Shares the `.cif` extension with mmCIF — the first
    /// format registered that does — disambiguated by content, not name,
    /// in [`crate::io::open`] and the CLI's own input path.
    pub const CIF_CORE: Format = Format(12);
    /// PSF (#321) — CHARMM/NAMD's topology format, see [`crate::io::psf`].
    /// The first real consumer of
    /// [`crate::core::force_field::ForceFieldTopology`] (#315, built for
    /// exactly this but never wired to a format until now). No coordinates
    /// at all — a PSF is paired with a DCD or a PDB for geometry.
    pub const PSF: Format = Format(13);
    /// PRMTOP (#322) -- AMBER's topology format, also named `.parm7`, see
    /// [`crate::io::prmtop`]. Topology only: PRMTOP states substantially
    /// more of the force field than PSF does, but the actual parameters
    /// (force constants, equilibrium values, Lennard-Jones coefficients)
    /// have no home in `ForceFieldTopology` and are parsed-and-discarded.
    pub const PRMTOP: Format = Format(14);
    /// GROMACS topology (#323), see [`crate::io::top`]. Topology only, same
    /// scope as PSF/PRMTOP. Shares the `.top` extension with PRMTOP,
    /// disambiguated by content (`top::is_gromacs_top`) at the `io::open`/
    /// CLI entry points, the same way `.cif` disambiguates mmCIF from CIF
    /// core. One `Molecule` per `[ moleculetype ]` block, unlike PSF/PRMTOP
    /// which hold exactly one topology each.
    pub const TOP: Format = Format(15);
    /// LAMMPS data (#324), see [`crate::io::lammps`]. The first
    /// force-field-topology format that also states coordinates -- PSF/
    /// PRMTOP/TOP never touch `coords3`, but LAMMPS combines topology and
    /// geometry in one file. Every atom reads as
    /// [`crate::core::atom::Element::UNKNOWN`] -- this format states a
    /// numeric type and mass, never a name, and inventing an element from
    /// mass would fail outright for coarse-grained/reduced-unit systems.
    pub const LAMMPS_DATA: Format = Format(16);
    /// TRR (#326), see [`crate::io::trr`]. The first format registered as
    /// [`Kind::Frames`] rather than [`Kind::Molecules`] -- a shared topology
    /// (every atom [`crate::core::atom::Element::UNKNOWN`], TRR states no
    /// chemical identity at all) plus lazily-read frames, each carrying
    /// positions and whatever else it happened to state: velocities,
    /// forces, a box. The first format with a real, populated `magic`
    /// signature, and the first with a trajectory writer (this format
    /// descriptor's own `writer_trajectory` field) rather than a
    /// `&[(String, Molecule)]`-shaped one.
    pub const TRR: Format = Format(17);
    /// XTC (#325), see [`crate::io::xtc`]. GROMACS's compressed trajectory
    /// format -- the same shared-topology `Kind::Frames` shape TRR
    /// established, wrapped around a genuine lossy coordinate compressor
    /// (a full, faithful port of GROMACS's own encoder heuristic, not a
    /// simplified always-absolute variant). No velocities or forces at
    /// all -- only positions, time, step and a box.
    pub const XTC: Format = Format(18);
    /// DCD (#327), see [`crate::io::dcd`]. CHARMM/NAMD's binary trajectory,
    /// and the oldest format in the milestone -- Fortran unformatted
    /// records, a per-file endianness this crate detects rather than
    /// assumes, and a header dialect flag governing both the timestep's
    /// storage type and whether a unit cell can be present at all. No
    /// velocities or forces, same as XTC.
    pub const DCD: Format = Format(19);
    /// NCTRAJ (#328), see [`crate::io::nctraj`]. Amber's NetCDF trajectory
    /// -- the Amber convention layered on a hand-rolled, generic NetCDF-3
    /// classic container ([`crate::io::netcdf3`]), the only route
    /// available since the real `netcdf` crate is a `-sys` crate this
    /// repo's own CI purity gate refuses. Coordinates are already Å, no
    /// unit conversion. Can carry velocities, like TRR; unlike TRR, no
    /// per-frame time or step at all.
    pub const NCTRAJ: Format = Format(20);
    /// LAMMPS Trajectory (#329), see [`crate::io::lammpstrj`]. The dump
    /// format's `ITEM:`-delimited sections, one set per frame -- the first
    /// [`Kind::Frames`] format that's text rather than binary. The column
    /// list is read fresh every frame, never cached, and a coordinate
    /// convention (unscaled/scaled/unwrapped) is resolved the same way.
    /// Can carry velocities and forces, like TRR; no per-frame time, only
    /// a step count, same reasoning DCD/XTC/NCTRAJ already settled for
    /// whichever of the two their own format states.
    pub const LAMMPS_TRAJECTORY: Format = Format(21);
    /// CUBE (#331), see [`crate::io::cube`]. Gaussian's volumetric format --
    /// the first ever registered as [`Kind::Volume`], and the first to
    /// carry both a grid (`Carries::SAMPLES`) and the atoms that produced
    /// it (`TOPOLOGY`/`COORDS_3D`, via [`crate::core::volume::VolumeGrid::atoms`]).
    pub const CUBE: Format = Format(22);
    /// CCP4/MRC (#332), see [`crate::io::ccp4`]. The standard electron-
    /// density/cryo-EM map format -- a real crystallographic `UnitCell`
    /// this time (unlike CUBE), and no atoms. Endianness comes from a
    /// direct machine-stamp byte-pattern match (no DCD-style guessing);
    /// `MAPC`/`MAPR`/`MAPS` permute which canonical axis each file
    /// dimension is, the second real exercise of
    /// [`crate::core::volume::VolumeGrid::from_source_order`].
    pub const CCP4: Format = Format(23);
    /// DX (#333), see [`crate::io::dx`]. OpenDX's grid format -- the
    /// simplest of the four volumetric formats: no atoms, no cell, no axis
    /// permutation. Confirms [`crate::core::volume::VolumeGrid::axes`]
    /// genuinely stores three full vectors, not three scalars -- DX's own
    /// deltas need not be axis-aligned.
    pub const DX: Format = Format(24);
    /// DSN6 (#334), see [`crate::io::dsn6`]. Frodo/O's bricked,
    /// byte-quantized electron-density format -- the last of the four
    /// volumetric formats, and the one with no live oracle available
    /// anywhere in this environment to verify against. A real
    /// crystallographic `UnitCell` (like CCP4), but a fixed axis mapping
    /// (no `MAPC`/`MAPR`/`MAPS`-style permutation) and density values
    /// quantized to a single byte per voxel via a `prod`/`plus` linear
    /// scale.
    pub const DSN6: Format = Format(25);
    /// OBJ (#335), see [`crate::io::obj`]. Wavefront's plain-text mesh
    /// format, first story of Phase 5 and the first registered format with
    /// no chemistry at all -- [`Kind::Mesh`], not `Kind::Molecules`. Faces
    /// are fan-triangulated into [`crate::core::mesh::Mesh`]'s
    /// triangles-only shape, and a vertex referenced with two different
    /// normals across faces (a hard edge) is split into two output
    /// vertices, since that type's normals are per-vertex, not per-face-
    /// corner like OBJ's own `v/vt/vn` indexing.
    pub const OBJ: Format = Format(26);
    /// PLY (#336), see [`crate::io::ply`]. The Stanford polygon format --
    /// second and last mesh story, and the first format whose header is a
    /// genuine schema (`element`/`property` lines) rather than a fixed
    /// field layout, the same class of work as PRMTOP's `%FORMAT` lines.
    /// One descriptor covers all three of its wire encodings (`ascii`,
    /// `binary_little_endian`, `binary_big_endian`), sniffed internally
    /// from the header's own `format` line -- this crate's first "ASCII or
    /// binary, both handled by one reader" format, unlike CCP4's
    /// little/big-endian-only sniff.
    pub const PLY: Format = Format(27);
    /// CSV (#337), see [`crate::io::csv`]. The tabular format, last story of
    /// Phase 5 -- declared [`Kind::Table`], with molecule production as an
    /// explicit, disclosed read option (`structure_column`) rather than a
    /// column-name heuristic, mirroring the exact "declared kind stays
    /// fixed, an option flips the produced payload" mechanism #330 already
    /// proved safe for XYZ/PDB/PDBQT/GRO's own opt-in `Frames` reading.
    /// Parses via [`crate::core::table::Table::from_csv`] (#314) directly,
    /// not a second RFC4180 parser. The first format to populate both
    /// `writer` (an opt-in molecule-list shape) and `writer_table` (its
    /// declared/primary shape) on the same descriptor.
    pub const CSV: Format = Format(28);

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
        Self::from_filename_checked(name).unwrap_or(Format::SMILES)
    }

    /// Like [`Self::from_filename`], but `None` rather than the SMILES
    /// default when nothing claims the extension — what a caller that wants
    /// to try something else first (content-sniffing, in
    /// `crate::io::open::open_supplier`, #317) needs instead of the default
    /// already baked in.
    pub fn from_filename_checked(name: &str) -> Option<Format> {
        name.rsplit_once('.')
            .and_then(|(_, extension)| Format::from_extension(extension))
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

    /// Whether this format's canonical bytes are UTF-8 text or binary.
    pub fn encoding(&self) -> Encoding {
        self.descriptor().encoding
    }

    /// What this format's records are — see [`Kind`].
    pub fn kind(&self) -> Kind {
        self.descriptor().kind
    }

    pub fn can_read(&self) -> bool {
        let d = self.descriptor();
        d.reader.is_some() || d.reader_bytes.is_some()
    }

    pub fn can_write(&self) -> bool {
        let d = self.descriptor();
        d.writer.is_some()
            || d.writer_bytes.is_some()
            || d.writer_trajectory.is_some()
            || d.writer_volume.is_some()
            || d.writer_mesh.is_some()
            || d.writer_table.is_some()
    }

    /// Parses a whole file into molecules, from raw bytes (#309) — the
    /// canonical read path every format goes through, text or binary, or
    /// `None` if the format cannot be read.
    ///
    /// A binary format's `reader_bytes` is called directly. A text format
    /// has none, so the bytes are decoded as UTF-8 first; invalid UTF-8
    /// becomes a `Skipped` entry rather than a panic, matching the rest of
    /// this crate's "reading a file cannot fail as a whole" contract.
    pub fn read_bytes(&self, bytes: &[u8]) -> Option<ReadOutcome> {
        self.read_bytes_with_options(bytes, &ReadOptions::default())
    }

    /// [`Self::read_bytes`], with explicit per-format options.
    pub fn read_bytes_with_options(
        &self,
        bytes: &[u8],
        options: &ReadOptions,
    ) -> Option<ReadOutcome> {
        let d = self.descriptor();
        if let Some(reader_bytes) = d.reader_bytes {
            return Some(reader_bytes(bytes, options));
        }
        d.reader.map(|reader| match std::str::from_utf8(bytes) {
            Ok(text) => reader(text, options),
            Err(e) => ReadOutcome {
                records: Vec::new(),
                skipped: vec![crate::io::reader::Skipped {
                    position: 1,
                    input: String::new(),
                    error: format!("{} is not valid UTF-8: {e}", self.name()),
                }],
            },
        })
    }

    /// Serialises molecules in this format with default options, carrying
    /// their names, or `None` if the format cannot be written as text.
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
        // `None` for a binary format (BinaryCIF, #319, the first one) rather
        // than forcing its bytes through `from_utf8` -- they are not text,
        // and were never going to decode as any. `write_bytes`/
        // `write_bytes_with_options` is the canonical path for a caller that
        // wants a binary format's own bytes.
        if self.encoding() == Encoding::Binary {
            return None;
        }
        self.write_bytes_with_options(records, options)
            .map(|bytes| String::from_utf8(bytes).expect("text writer produced valid UTF-8"))
    }

    /// Serialises molecules in this format into raw bytes (#309) — the
    /// canonical write path every format goes through, text or binary, or
    /// `None` if the format cannot be written.
    pub fn write_bytes(&self, records: &[(String, Molecule)]) -> Option<Vec<u8>> {
        self.write_bytes_with_options(records, &WriteOptions::default())
    }

    /// [`Self::write_bytes`], with explicit per-format options.
    pub fn write_bytes_with_options(
        &self,
        records: &[(String, Molecule)],
        options: &WriteOptions,
    ) -> Option<Vec<u8>> {
        let d = self.descriptor();
        if let Some(writer_bytes) = d.writer_bytes {
            return Some(writer_bytes(records, options));
        }
        d.writer.map(|writer| writer(records, options).into_bytes())
    }

    /// Serialises a whole trajectory into raw bytes (#326), or `None` if
    /// this format has no trajectory writer — every `Kind::Molecules`
    /// format, and any `Kind::Frames` format that has not implemented one.
    pub fn write_trajectory_bytes(&self, trajectory: &mut Trajectory) -> Option<Vec<u8>> {
        self.write_trajectory_bytes_with_options(trajectory, &WriteOptions::default())
    }

    /// [`Self::write_trajectory_bytes`], with explicit per-format options.
    pub fn write_trajectory_bytes_with_options(
        &self,
        trajectory: &mut Trajectory,
        options: &WriteOptions,
    ) -> Option<Vec<u8>> {
        let writer_trajectory = self.descriptor().writer_trajectory?;
        Some(writer_trajectory(trajectory, options))
    }

    /// Serialises a whole [`VolumeGrid`] into raw bytes (#331), or `None` if
    /// this format has no volume writer — every `Kind::Molecules`/
    /// `Kind::Frames` format, and any `Kind::Volume` format that has not
    /// implemented one.
    pub fn write_volume_bytes(&self, grid: &VolumeGrid) -> Option<Vec<u8>> {
        self.write_volume_bytes_with_options(grid, &WriteOptions::default())
    }

    /// [`Self::write_volume_bytes`], with explicit per-format options.
    pub fn write_volume_bytes_with_options(
        &self,
        grid: &VolumeGrid,
        options: &WriteOptions,
    ) -> Option<Vec<u8>> {
        let writer_volume = self.descriptor().writer_volume?;
        Some(writer_volume(grid, options))
    }

    /// Serialises a whole [`Mesh`] into raw bytes (#335), or `None` if this
    /// format has no mesh writer — every `Kind::Molecules`/`Kind::Frames`/
    /// `Kind::Volume` format, and any `Kind::Mesh` format that has not
    /// implemented one.
    pub fn write_mesh_bytes(&self, mesh: &Mesh) -> Option<Vec<u8>> {
        self.write_mesh_bytes_with_options(mesh, &WriteOptions::default())
    }

    /// [`Self::write_mesh_bytes`], with explicit per-format options.
    pub fn write_mesh_bytes_with_options(
        &self,
        mesh: &Mesh,
        options: &WriteOptions,
    ) -> Option<Vec<u8>> {
        let writer_mesh = self.descriptor().writer_mesh?;
        Some(writer_mesh(mesh, options))
    }

    /// Serialises a whole [`Table`] into raw bytes (#337), or `None` if this
    /// format has no table writer — every format other than CSV today.
    pub fn write_table_bytes(&self, table: &Table) -> Option<Vec<u8>> {
        self.write_table_bytes_with_options(table, &WriteOptions::default())
    }

    /// [`Self::write_table_bytes`], with explicit per-format options.
    pub fn write_table_bytes_with_options(
        &self,
        table: &Table,
        options: &WriteOptions,
    ) -> Option<Vec<u8>> {
        let writer_table = self.descriptor().writer_table?;
        Some(writer_table(table, options))
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

/// Attributes a writer manufactures when the input has none.
///
/// A structural format's atom line has fixed columns, so writing one means
/// putting *something* in them. The value is invented, and a reader cannot
/// tell it from a measured one -- `chem convert x.smi --to pdb` produces a
/// B-factor of 0.00 for every atom, which is the same class of defect #173
/// records against OpenBabel, arriving from the other direction.
///
/// Keyed by target and attribute rather than by pair, because that is how it
/// behaves: every source lacking the attribute gains it, and no source
/// carrying it is affected. `test_every_format_pair_carries_the_intersection_of_its_masks`
/// asserts this is exactly the set that appears.
static SUPPLIED: &[(Format, Carries, &str)] = &[
    // Coordinates used to be here -- all seven of them, five claiming a
    // conformer and two a drawing, and every one of them zeros a writer had put
    // in a column it could not leave empty. They came out with #270: the readers
    // no longer mistake a format's placeholder for a measurement, so there is
    // nothing supplied left to record.
    // Residue identity. Every atom line names one; absent input writes a
    // placeholder, `LIG` for a docking ligand and `UNK` elsewhere.
    (
        Format::PDB,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    (
        Format::MMCIF,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    (
        // Same reason as MMCIF, verbatim: both go through
        // `cif_model::build_rows` (#319), so what one supplies, the other
        // does too, by construction.
        Format::BCIF,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    (
        Format::MOL2,
        Carries::RESIDUES,
        "every atom line names a substructure; absent input writes UNK",
    ),
    (
        Format::GRO,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    (
        Format::PDBQT,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes LIG",
    ),
    (
        Format::PSF,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    (
        // `write_prmtop` (#322) always states a residue label, the same
        // reason PSF does.
        Format::PRMTOP,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    (
        // `write_top` (#323) always states a residue name, the same
        // reason PSF/PRMTOP do.
        Format::TOP,
        Carries::RESIDUES,
        "every atom line names a residue; absent input writes UNK",
    ),
    // PSF's atom-type and mass columns are required fields in every real
    // file, the same "must state something" shape as the columns below --
    // absent input falls back to the element symbol and its standard
    // atomic weight.
    (
        Format::PSF,
        Carries::ATOM_TYPE,
        "fixed column; absent input writes the element symbol",
    ),
    (
        Format::PSF,
        Carries::MASS,
        "fixed column; absent input writes the standard atomic weight",
    ),
    // PRMTOP's AMBER_ATOM_TYPE/MASS sections are equally required fields
    // in every real file, the same shape as PSF's own columns above.
    (
        Format::PRMTOP,
        Carries::ATOM_TYPE,
        "fixed column; absent input writes the element symbol",
    ),
    (
        Format::PRMTOP,
        Carries::MASS,
        "fixed column; absent input writes the standard atomic weight",
    ),
    // GROMACS TOP's `type`/`mass` columns are equally required fields in
    // every real file, the same shape as PSF/PRMTOP's own columns above.
    (
        Format::TOP,
        Carries::ATOM_TYPE,
        "fixed column; absent input writes the element symbol",
    ),
    (
        Format::TOP,
        Carries::MASS,
        "fixed column; absent input writes the standard atomic weight",
    ),
    // The PDB family's fixed occupancy and temperature-factor columns.
    (
        Format::PDB,
        Carries::OCCUPANCY,
        "fixed column; absent input writes 1.00",
    ),
    (
        Format::MMCIF,
        Carries::OCCUPANCY,
        "fixed column; absent input writes 1.00",
    ),
    (
        // Same reason as MMCIF -- shares `cif_model::build_rows` (#319).
        Format::BCIF,
        Carries::OCCUPANCY,
        "fixed column; absent input writes 1.00",
    ),
    (
        Format::PDBQT,
        Carries::OCCUPANCY,
        "fixed column; absent input writes 1.00",
    ),
    (
        // `write_cif_core` (#320) defaults occupancy to 1.00 the same way
        // mmCIF's writer does, for the same reason.
        Format::CIF_CORE,
        Carries::OCCUPANCY,
        "fixed column; absent input writes 1.00",
    ),
    (
        Format::PDB,
        Carries::B_FACTOR,
        "fixed column; absent input writes 0.00",
    ),
    (
        Format::MMCIF,
        Carries::B_FACTOR,
        "fixed column; absent input writes 0.00",
    ),
    (
        // Same reason as MMCIF -- shares `cif_model::build_rows` (#319).
        Format::BCIF,
        Carries::B_FACTOR,
        "fixed column; absent input writes 0.00",
    ),
    (
        Format::PDBQT,
        Carries::B_FACTOR,
        "fixed column; absent input writes 0.00",
    ),
    (
        // `write_cif_core` (#320) defaults B-factor to 0.00 the same way
        // mmCIF's writer does, for the same reason.
        Format::CIF_CORE,
        Carries::B_FACTOR,
        "fixed column; absent input writes 0.00",
    ),
    // A per-atom charge column, written whether or not one was computed.
    (
        Format::MOL2,
        Carries::PARTIAL_CHARGE,
        "per-atom charge column; absent input writes 0.000",
    ),
    (
        Format::PDBQT,
        Carries::PARTIAL_CHARGE,
        "per-atom charge column; absent input writes 0.000",
    ),
    (
        // `write_psf` (#321) always states a charge, the same reason Mol2
        // and PDBQT do.
        Format::PSF,
        Carries::PARTIAL_CHARGE,
        "per-atom charge column; absent input writes 0.000",
    ),
    (
        // `write_prmtop` (#322) always states a charge, the same reason
        // PSF does.
        Format::PRMTOP,
        Carries::PARTIAL_CHARGE,
        "per-atom charge column; absent input writes 0.000",
    ),
    (
        // `write_top` (#323) always states a charge, the same reason
        // PSF/PRMTOP do.
        Format::TOP,
        Carries::PARTIAL_CHARGE,
        "per-atom charge column; absent input writes 0.000",
    ),
    (
        // `write_lammps_data` (#324) always states a position, even for a
        // source with no coordinates at all -- the origin, the same
        // "something has to go here" reasoning `cif_core`'s own
        // cell-less-molecule writer already uses.
        Format::LAMMPS_DATA,
        Carries::COORDS_3D,
        "every atom line states a position; absent input writes 0.0 0.0 0.0",
    ),
    (
        // Every atom's LAMMPS type/mass is synthesized on write, whether
        // or not the source stated one -- the same shape PSF/PRMTOP/TOP's
        // own atom-type/mass columns already have.
        Format::LAMMPS_DATA,
        Carries::ATOM_TYPE,
        "a type is always synthesized; absent input still gets one",
    ),
    (
        Format::LAMMPS_DATA,
        Carries::MASS,
        "fixed column; absent input writes the standard atomic weight",
    ),
    (
        // Box bounds are mandatory in every real LAMMPS data file; a
        // source with no periodic cell at all still gets a placeholder
        // box (a bounding box of its coordinates, padded to be
        // non-degenerate -- see `write_lammps_data`'s own doc comment).
        Format::LAMMPS_DATA,
        Carries::UNIT_CELL,
        "box bounds are mandatory; absent input gets a placeholder box",
    ),
];

/// Attributes a conversion delivers that neither mask claims.
///
/// Pair-specific, unlike [`SUPPLIED`]: what a writer manufactures depends on
/// what actually arrived, not only on the target's columns.
static PAIR_EXTRAS: &[(Format, Format, Carries, &str)] = &[
    // A molfile records a stereo field for every bond, and a double bond
    // claiming nothing comes back `BondStereo::Unspecified` rather than
    // `None` -- an explicit "not stated", which is #198's fix rather than a
    // regression of it, and which `held` counts as bond stereo.
    //
    // Pair-specific because it needs a *double* bond to arrive, which is not a
    // `Carries` flag: PDB carries bonds too, but writes only CONECT singles, so
    // `pdb -> sdf` gains nothing here. Mol2 and CML are the two bond-carrying
    // formats that state a bond order without claiming STEREO_BOND.
    (
        Format::MOL2,
        Format::SDF,
        Carries::STEREO_BOND,
        "a molfile states, explicitly, that an arriving double bond claims nothing",
    ),
    (
        Format::CML,
        Format::SDF,
        Carries::STEREO_BOND,
        "a molfile states, explicitly, that an arriving double bond claims nothing",
    ),
];

/// Attributes both masks claim but a conversion loses anyway, with the issue
/// that owns each.
///
/// Pinned rather than deleted, the same discipline as the corpus's
/// `known-gap-*` entries -- closing one fails the test that pins it, which is
/// the prompt to remove the line.
static PAIR_LOSSES: &[(Format, Format, Carries, &str)] = &[
    // Not a defect: PDBQT reads no bonds at all, so a bond-based target has
    // nothing for atom aromaticity to ride on. Attribute interference, exactly
    // what `one_per_attribute`'s doc warns about -- here between AROMATICITY
    // and BONDS.
    (
        Format::PDBQT,
        Format::SDF,
        Carries::AROMATICITY,
        "no bonds survive PDBQT for aromaticity to ride on",
    ),
    (
        Format::PDBQT,
        Format::CML,
        Carries::AROMATICITY,
        "no bonds survive PDBQT for aromaticity to ride on",
    ),
    // Four rows came out with #261. Three were CML, whose only aromatic
    // channel is the bond order, and one was Mol2 -- whose row blamed the
    // missing bond flag when the real cause was an unstated hydrogen count
    // reaching `kekulize` (#281). Readers now reconcile all three channels at
    // the boundary, so neither format loses what both masks claim.
];

/// Conversions that lose *atoms*, each naming the issue that owns it.
///
/// A pair, not a format: what breaks here is a property of the conversion, so
/// no per-format mask can express it -- and `held` sets `TOPOLOGY` on the atom
/// count alone, so one atom of six satisfies every claim a mask makes.
static PAIR_GAPS: &[(Format, Format, &str)] = &[
    // Empty from #259 to #324. PDBQT's writer kept only the largest connected
    // component, so a source carrying no bonds arrived as N one-atom fragments
    // and left as one atom -- four pairs, `pdbqt -> pdbqt` among them, which is
    // why the A -> A diagonal could not see it. The writer now writes every
    // atom, so there is no gap to pin.
    //
    // The mechanism stays: this is the only way to express a loss no mask can,
    // and it took a milestone to notice the first one.
    //
    // LAMMPS data (#324) states no element at all, only a numeric type --
    // every atom reads back as `Element::UNKNOWN`, whose symbol is the empty
    // string. Every format below needs a real atomic symbol to write an atom
    // at all, so a LAMMPS-sourced molecule fails to round-trip through any of
    // them. Not a bug on either side: the format genuinely never stated what
    // these targets need to recover.
    (
        Format::LAMMPS_DATA,
        Format::SMILES,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::SDF,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::CXSMILES,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::XYZ,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::PDB,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::MMCIF,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::MOL2,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::PDBQT,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::GRO,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::CML,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::BCIF,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::CIF_CORE,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::PSF,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
    (
        Format::LAMMPS_DATA,
        Format::TOP,
        "LAMMPS states no element, only a numeric type; Element::UNKNOWN has no atomic symbol to write",
    ),
];

/// Everything `target`'s writer manufactures when the input has none.
pub fn supplied(target: Format) -> Carries {
    SUPPLIED
        .iter()
        .filter(|(format, _, _)| *format == target)
        .fold(Carries::empty(), |acc, (_, flag, _)| acc | *flag)
}

/// Attributes `source` to `target` delivers that neither mask claims.
pub fn pair_extra(source: Format, target: Format) -> Carries {
    PAIR_EXTRAS
        .iter()
        .filter(|(from, to, _, _)| *from == source && *to == target)
        .fold(Carries::empty(), |acc, (_, _, flag, _)| acc | *flag)
}

/// Why `source` to `target` loses what [`pair_loss`] reports.
pub fn pair_loss_reason(source: Format, target: Format) -> Option<&'static str> {
    PAIR_LOSSES
        .iter()
        .find(|(from, to, _, _)| *from == source && *to == target)
        .map(|(_, _, _, why)| *why)
}

/// Attributes `source` to `target` loses despite both masks claiming them.
pub fn pair_loss(source: Format, target: Format) -> Carries {
    PAIR_LOSSES
        .iter()
        .filter(|(from, to, _, _)| *from == source && *to == target)
        .fold(Carries::empty(), |acc, (_, _, flag, _)| acc | *flag)
}

/// Why converting `source` to `target` loses atoms, or `None` when it does not.
///
/// The drop report answers "can the target hold this?" one hop at a time. A
/// conversion is two, and a handful of pairs lose something no mask mentions.
pub fn pair_gap(source: Format, target: Format) -> Option<&'static str> {
    PAIR_GAPS
        .iter()
        .find(|(from, to, _)| *from == source && *to == target)
        .map(|(_, _, why)| *why)
}

/// What a conversion from `source` to `target` keeps of what a molecule holds.
///
/// Subtract this from [`held`] and the remainder is the conversion's losses:
///
/// ```
/// # use chem::io::format::{self, Format, held};
/// # let molecule = chem::io::smiles::parse_smiles("c1ccccc1").unwrap();
/// let dropped = held(&molecule).difference(format::kept(Format::CML, Format::SMILES));
/// ```
///
/// `target.carries()` alone is not the answer, which is the whole of #276: a
/// conversion is a read and a write, and #257 pinned six pairs that lose an
/// attribute *both* masks claim. CML to SMILES hands back cyclohexane while a
/// target-only report says nothing (#261).
///
/// Distinct from [`fidelity`], which folds in what the target *manufactures*
/// and so answers what the output will contain. This answers what the input
/// keeps -- a drop report needs the second, and using the first would under-
/// report exactly the pairs this exists for.
///
/// One function because there were two: the workbench had this formula and the
/// command line had `target.carries()`, and the two disagreeing is what #276
/// was filed about.
pub fn kept(source: Format, target: Format) -> Carries {
    target.carries().difference(pair_loss(source, target))
}

/// What a conversion from `source` to `target` actually delivers.
///
/// The two masks intersected, less what the pair is known to lose, plus what
/// the target's writer supplies of its own. This is one cell of the fidelity
/// matrix, and `test_every_format_pair_carries_the_intersection_of_its_masks`
/// is what keeps it true.
pub fn fidelity(source: Format, target: Format) -> Carries {
    let both = source.carries() & target.carries();
    let manufactured = supplied(target) & target.carries();
    (both | manufactured | pair_extra(source, target)) & !pair_loss(source, target)
}

/// A short, human phrase for what a `Kind`'s records *are* — used only in
/// [`kinds_compatible`]'s own refusal message.
fn kind_description(kind: Kind) -> &'static str {
    match kind {
        Kind::Molecules => "a set of molecules",
        Kind::Frames => "a trajectory",
        Kind::Volume => "a volumetric grid",
        Kind::Mesh => "a mesh",
        Kind::Table => "a table",
    }
}

/// Whether `source` can feed `target` at all — the "refuse before reading
/// any bytes" gate `chem convert` needs now that the registry has five
/// `Kind`s (#338).
///
/// `Ok(())` for a same-`Kind` pair, and exactly one deliberately enumerated
/// exception: a `Kind::Volume` source that also states real atoms
/// (`Carries::TOPOLOGY`) can feed a `Kind::Molecules` target — CUBE's own
/// dual nature (#331), and the *only* one in the whole registry, confirmed
/// directly against every other `Kind::Volume`/`Mesh`/`Table` format's own
/// mask before writing this. Never the reverse direction: no format can
/// synthesize real grid samples from a bare molecule.
///
/// Every other cross-`Kind` pair is a genuine mismatch, refused by name
/// rather than left to degrade into a per-record "not a molecule" skip and
/// a generic empty-output exit — which is what every one of these pairs did
/// before this function existed.
pub fn kinds_compatible(source: Format, target: Format) -> Result<(), String> {
    if source.kind() == target.kind() {
        return Ok(());
    }
    if source.kind() == Kind::Volume
        && source.carries().contains(Carries::TOPOLOGY)
        && target.kind() == Kind::Molecules
    {
        return Ok(());
    }
    Err(format!(
        "{} holds {}, {} holds {} -- there is no conversion between them",
        source.name(),
        kind_description(source.kind()),
        target.name(),
        kind_description(target.kind()),
    ))
}

/// What a trajectory actually holds, the same "inspect the real instance,
/// not the format's declared mask" discipline [`held`] already follows for
/// a molecule (#338).
///
/// `TOPOLOGY`/`COORDS_3D` unconditionally — a trajectory with no positions
/// is not a trajectory. `VELOCITIES`/`FORCES`/`FRAME_TIME`/`UNIT_CELL` only
/// if frame 0 actually states them, since [`crate::core::trajectory::Frame`]
/// keeps every one of those fields independently `Option` (TRR carries
/// velocities and forces, XTC never does) and this crate's own binary
/// trajectory formats never vary that per frame within one file.
///
/// # Errors
/// Whatever reading frame 0 itself can fail with.
pub fn held_from_trajectory(trajectory: &mut Trajectory) -> Result<Carries, TrajectoryError> {
    let mut carries = Carries::TOPOLOGY.or(Carries::COORDS_3D);
    if trajectory.frame_count() > 0 {
        let frame = trajectory.frame(0)?;
        if frame.velocities.is_some() {
            carries = carries.or(Carries::VELOCITIES);
        }
        if frame.forces.is_some() {
            carries = carries.or(Carries::FORCES);
        }
        if frame.time.is_some() {
            carries = carries.or(Carries::FRAME_TIME);
        }
        if frame.cell.is_some() {
            carries = carries.or(Carries::UNIT_CELL);
        }
    }
    Ok(carries)
}

/// What a volume grid actually holds, the same discipline [`held`]/
/// [`held_from_trajectory`] already follow (#338). `SAMPLES` unconditionally
/// — a grid with no values is not a grid — plus `UNIT_CELL`/`TOPOLOGY.or(
/// COORDS_3D)` only if this specific grid actually states a cell/atoms
/// ([`VolumeGrid::cell`]/[`VolumeGrid::atoms`] — CUBE states atoms,
/// CCP4/DSN6 state a cell but no atoms, DX states neither).
pub fn held_from_volume(grid: &VolumeGrid) -> Carries {
    let mut carries = Carries::SAMPLES;
    if grid.cell().is_some() {
        carries = carries.or(Carries::UNIT_CELL);
    }
    if grid.atoms().is_some() {
        carries = carries.or(Carries::TOPOLOGY).or(Carries::COORDS_3D);
    }
    carries
}

/// Whether `bytes` contains `pattern` starting at `offset` — never panics on
/// a buffer shorter than `offset + pattern.len()`, since `slice::get` on an
/// out-of-range range answers `None` rather than indexing (#317).
fn has_magic_at(bytes: &[u8], offset: usize, pattern: &[u8]) -> bool {
    bytes.get(offset..offset + pattern.len()) == Some(pattern)
}

/// The dispatch [`sniff`] runs, factored out so it can be proven against a
/// hand-built candidate list rather than the real registry (#317) — a
/// `Format` is always a real, registered handle (see its own doc comment),
/// so there is no way to hand this a fake one the way a bare
/// [`FormatDescriptor`] could stand in for an unregistered `Kind` in #310's
/// own tests. Pairing a *real* `Format` with a *made-up* signature list
/// sidesteps that: the matching logic under test is identical either way.
fn find_signature(
    bytes: &[u8],
    candidates: impl Iterator<Item = (Format, &'static [Signature])>,
) -> Option<Format> {
    candidates
        .filter(|(_, sigs)| {
            sigs.iter()
                .any(|sig| has_magic_at(bytes, sig.offset, sig.bytes))
        })
        .map(|(format, _)| format)
        .next()
}

/// Matches `bytes` against every registered format's content signature, or
/// `None` if nothing recognizes it (#317).
///
/// The last resort `crate::io::open::open_supplier` reaches for before
/// `Format::from_filename`'s SMILES default, and only when the file's own
/// extension claimed nothing. Every registered format's `magic` is empty
/// today — #309 built the byte-reading path a binary format needs, but none
/// has registered a signature yet (#319, #325-#337) — so this always
/// answers `None`, the same "the mechanism lands now, a later story
/// populates it" shape every Phase-1 story since #311 has followed.
pub(crate) fn sniff(bytes: &[u8]) -> Option<Format> {
    find_signature(bytes, all().map(|f| (f, f.descriptor().magic)))
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
        // Another case that used to fall back: `.mol` files register `mol`
        // as a code (`-imol` already worked) but not as an extension, so
        // every line of a real molfile silently skipped. #318 added it.
        assert_eq!(Format::from_filename("a.mol"), Format::SDF);
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
        use crate::core::force_field::{ForceFieldAtom, ForceFieldTopology};
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
        let with_force_field_atom = |mutate: fn(&mut ForceFieldAtom)| {
            let mut m = ethane();
            let mut atom = ForceFieldAtom::empty();
            mutate(&mut atom);
            m.set_force_field(ForceFieldTopology {
                atoms: Some(vec![atom, ForceFieldAtom::empty()]),
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };
        // Three and four atoms respectively -- ethane's two aren't enough
        // for an angle or a dihedral/improper term. Propane and butane's
        // real, linear connectivity is exactly the shape an angle and a
        // proper torsion need; isobutane's branch point is the shape an
        // improper needs (a central atom plus three substituents).
        let with_angle = {
            let mut m = parse_smiles("CCC").expect("valid SMILES");
            m.set_force_field(ForceFieldTopology {
                angles: vec![[0, 1, 2]],
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };
        let with_dihedral = {
            let mut m = parse_smiles("CCCC").expect("valid SMILES");
            m.set_force_field(ForceFieldTopology {
                dihedrals: vec![[0, 1, 2, 3]],
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };
        let with_improper = {
            let mut m = parse_smiles("CC(C)C").expect("valid SMILES");
            m.set_force_field(ForceFieldTopology {
                impropers: vec![[0, 1, 2, 3]],
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };
        let with_exclusion = {
            let mut m = ethane();
            m.set_force_field(ForceFieldTopology {
                exclusions: vec![[0, 1]],
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };
        let with_donor = {
            let mut m = ethane();
            m.set_force_field(ForceFieldTopology {
                donors: vec![[0, 1]],
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };
        let with_acceptor = {
            let mut m = ethane();
            m.set_force_field(ForceFieldTopology {
                acceptors: vec![[0, 1]],
                ..ForceFieldTopology::default()
            })
            .expect("valid force field");
            m
        };

        vec![
            (Carries::TOPOLOGY, ethane()),
            // Ethane holds TOPOLOGY and BONDS both, which is fine: each row
            // asserts only its own flag. Isolating bonds from atoms is not
            // possible anyway -- a bond needs two atoms to join.
            (Carries::BONDS, ethane()),
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
            (
                Carries::ATOM_TYPE,
                with_force_field_atom(|a| a.atom_type = Some("CT".to_string())),
            ),
            (
                Carries::MASS,
                with_force_field_atom(|a| a.mass = Some(12.011)),
            ),
            (Carries::ANGLES, with_angle),
            (Carries::DIHEDRALS, with_dihedral),
            (Carries::IMPROPERS, with_improper),
            (Carries::EXCLUSIONS, with_exclusion),
            (Carries::DONORS, with_donor),
            (Carries::ACCEPTORS, with_acceptor),
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
    fn test_partial_charge_is_detected_from_a_force_field_alone() {
        // `PARTIAL_CHARGE` already has a fixture via `AtomSite` above; this
        // is the other source (#315) -- proving the `OR` in `held` actually
        // fires, not just that it compiles, by using a molecule with no
        // site data at all.
        use crate::core::force_field::{ForceFieldAtom, ForceFieldTopology};
        use crate::io::smiles::parse_smiles;

        let mut m = parse_smiles("CC").expect("valid SMILES");
        assert!(m.sites().is_none());
        m.set_force_field(ForceFieldTopology {
            atoms: Some(vec![
                ForceFieldAtom {
                    partial_charge: Some(-0.1),
                    ..ForceFieldAtom::default()
                },
                ForceFieldAtom::empty(),
            ]),
            ..ForceFieldTopology::default()
        })
        .expect("valid force field");

        assert!(held(&m).contains(Carries::PARTIAL_CHARGE));
    }

    #[test]
    fn test_declared_masks_match_what_actually_survives() {
        // The masks were derived by reading the writers, which is exactly the
        // kind of claim that rots. So they are checked against reality: a mask
        // that overstates makes the drop report lie about the very data it
        // exists to protect.
        for format in all() {
            // `one_per_attribute`'s fixtures are all `Molecule`s -- a
            // `Kind::Frames` format (TRR, #326) has no Molecule-shaped
            // writer to probe here; it is proven correct on its own terms
            // in `crate::io::trr`'s own test module instead.
            if format.kind() != Kind::Molecules || !format.can_write() || !format.can_read() {
                continue;
            }
            for (flag, molecule) in one_per_attribute() {
                let records = vec![("probe".to_string(), molecule)];
                let bytes = format.write_bytes(&records).expect("can_write said so");
                let outcome = format.read_bytes(&bytes).expect("can_read said so");
                let Some(back) = outcome.records.first() else {
                    panic!("{format:?} wrote nothing readable for {flag:?}");
                };
                let survived = held(back.molecule().expect("fixture format is Kind::Molecules"))
                    .contains(flag);
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

    /// One conversion, exactly as `chem convert` performs it: read the source
    /// file, write the target, read it back.
    ///
    /// Goes through `write_bytes`/`read_bytes` rather than the text-only
    /// `write`/`reader::read` -- the canonical path every format goes
    /// through regardless of encoding (#309), and the only one BinaryCIF
    /// (#319) can use at all, since its bytes are not valid UTF-8.
    ///
    /// Returns `None` when a format writes something it cannot read back,
    /// which is a failure rather than a loss and is reported as one.
    fn convert(source: Format, target: Format, molecule: &Molecule) -> Option<Molecule> {
        let records = vec![("probe".to_string(), molecule.clone())];
        let as_source = source.write_bytes(&records)?;
        let read_source = source.read_bytes(&as_source)?;
        let intermediate = read_source
            .records
            .first()?
            .molecule()
            .expect("fixture format is Kind::Molecules");

        let records = vec![("probe".to_string(), intermediate.clone())];
        let as_target = target.write_bytes(&records)?;
        let read_target = target.read_bytes(&as_target)?;
        Some(
            read_target
                .records
                .first()?
                .molecule()
                .expect("fixture format is Kind::Molecules")
                .clone(),
        )
    }

    /// Every ordered `(source, target)` pair sharing a `Kind`, restricted to
    /// `Kind::Molecules` -- the only kind [`one_per_attribute`]'s fixtures
    /// (and this whole fidelity matrix) are shaped for. TRR (#326) is the
    /// first format registered outside it; generalising this matrix to work
    /// across kinds is #338/#339's job, not this story's.
    fn pairs_within_kind() -> impl Iterator<Item = (Format, Format)> {
        all()
            .filter(|f| f.kind() == Kind::Molecules)
            .flat_map(|source| {
                all()
                    .filter(move |target| target.kind() == source.kind())
                    .map(move |target| (source, target))
            })
    }

    #[test]
    fn test_every_format_pair_delivers_what_the_matrix_says() {
        // The masks are per-format; a conversion is a pair. What survives
        // A -> B is what both masks claim, less what the pair is known to lose
        // and plus what B's writer supplies of its own -- and until #257
        // nothing checked it, because
        // `test_declared_masks_match_what_actually_survives` only walks the
        // diagonal.
        let fixtures = one_per_attribute();
        let pairs: Vec<(Format, Format)> = pairs_within_kind().collect();
        // 17 `Kind::Molecules` formats, 17 x 17 -- TRR (#326) is the first
        // format `pairs_within_kind` excludes, being the first registered
        // outside that kind, so this count is unchanged by its arrival.
        assert_eq!(pairs.len(), 289);
        for (source, target) in pairs {
            let predicted_mask = fidelity(source, target);
            for (flag, molecule) in &fixtures {
                let back = match convert(source, target, molecule) {
                    Some(back) => back,
                    // A pair pinned in `PAIR_GAPS` may fail to produce
                    // anything readable at all, not just lose an
                    // attribute -- LAMMPS's `Element::UNKNOWN` has no
                    // atomic symbol for a target that requires one. Any
                    // other pair failing here is a real bug, not a known
                    // gap, so it still panics.
                    None => {
                        assert!(
                            pair_gap(source, target).is_some(),
                            "{} -> {} wrote nothing readable and is not pinned in PAIR_GAPS",
                            source.name(),
                            target.name()
                        );
                        continue;
                    }
                };

                let predicted = predicted_mask.contains(*flag);
                let survived = held(&back).contains(*flag);

                assert_eq!(
                    predicted,
                    survived,
                    "{} -> {}: the matrix says {flag:?}={predicted} but the conversion gives {survived}",
                    source.name(),
                    target.name()
                );
            }
        }
    }

    #[test]
    fn test_a_converted_molecule_can_still_be_drawn() {
        // The defect #270 was reported for, end to end: a SMILES molecule
        // converted to any format and read back must still produce a usable
        // layout. It did not, and every mask assertion in this file passed
        // anyway -- `held` reports COORDS_2D whether the layout distinguishes
        // the atoms or stacks them all on one point.
        //
        // Counting *distinct* positions is what makes this fail. `has_coords`
        // was true the whole time.
        let outcome = crate::io::reader::read("CCO ethanol\n", Format::SMILES);
        let records: Vec<(String, Molecule)> = outcome
            .records
            .iter()
            .map(|r| {
                (
                    r.name.clone(),
                    r.molecule()
                        .expect("fixture format is Kind::Molecules")
                        .clone(),
                )
            })
            .collect();

        for format in all().filter(|f| f.kind() == Kind::Molecules) {
            let bytes = format.write_bytes(&records).expect("every format writes");
            let back = format
                .read_bytes(&bytes)
                .expect("every format reads its own bytes");
            let mut molecule = back.records[0]
                .molecule()
                .expect("fixture format is Kind::Molecules")
                .clone();
            crate::core::layout::ensure_coords(&mut molecule);

            let coords = molecule.coords().expect("a layout, computed or read");
            let mut seen: Vec<String> = coords
                .iter()
                .map(|p| format!("{:.3},{:.3}", p.x, p.y))
                .collect();
            seen.sort();
            seen.dedup();

            assert_eq!(
                seen.len(),
                molecule.num_atoms(),
                "{} leaves atoms on top of each other, so it cannot be drawn",
                format.name()
            );
        }
    }

    #[test]
    fn test_no_matrix_exception_is_stale() {
        // Both tables are claims about behaviour, and a claim nothing exercises
        // rots. Each row must be reachable: a supplied attribute the writer
        // does not actually manufacture, or a loss that stopped happening,
        // would otherwise sit in the table describing a crate that moved on.
        let fixtures = one_per_attribute();

        for (target, flag, why) in SUPPLIED {
            let fires = all().any(|source| {
                !source.carries().contains(*flag)
                    && fixtures.iter().any(|(_, molecule)| {
                        convert(source, *target, molecule)
                            .is_some_and(|back| held(&back).contains(*flag))
                    })
            });
            assert!(
                fires,
                "{} is pinned as supplying {flag:?} ({why:?}) but never does -- delete the line",
                target.name()
            );
        }

        for (source, target, flag, why) in PAIR_EXTRAS {
            assert!(
                !(source.carries().contains(*flag) && target.carries().contains(*flag)),
                "{} -> {} is pinned as an extra {flag:?} ({why:?}), but both masks already claim it",
                source.name(),
                target.name()
            );
        }

        for (source, target, flag, why) in PAIR_LOSSES {
            assert!(
                source.carries().contains(*flag) && target.carries().contains(*flag),
                "{} -> {} is pinned as losing {flag:?} ({why:?}), but a mask never claimed it",
                source.name(),
                target.name()
            );
        }
    }

    #[test]
    fn test_every_reader_leaves_the_aromaticity_channels_agreeing() {
        // The invariant, as a rule rather than a CML test. Aromaticity travels
        // on three channels and readers used to set different subsets, so a
        // file that stated it correctly came back as a different compound --
        // benzene through CML wrote out as cyclohexane, because the SMILES
        // writer reads only the atom flag (#261).
        //
        // Nothing here suspects a particular format, which is the point: this
        // would have caught it without anyone thinking to look at CML.
        let outcome = crate::io::reader::read("c1ccccc1 benzene\n", Format::SMILES);
        let records: Vec<(String, Molecule)> = outcome
            .records
            .iter()
            .map(|r| {
                (
                    r.name.clone(),
                    r.molecule()
                        .expect("fixture format is Kind::Molecules")
                        .clone(),
                )
            })
            .collect();

        for format in all().filter(|f| f.kind() == Kind::Molecules) {
            let bytes = format.write_bytes(&records).expect("every format writes");
            let back = format
                .read_bytes(&bytes)
                .expect("every format reads its own bytes");
            let molecule = back.records[0]
                .molecule()
                .expect("fixture format is Kind::Molecules");

            let aromatic_atoms = molecule.atoms().iter().filter(|a| a.is_aromatic()).count();
            let flagged = molecule.bonds().iter().filter(|b| b.is_aromatic()).count();
            let ordered = molecule
                .bonds()
                .iter()
                .filter(|b| b.order() == BondOrder::Aromatic)
                .count();

            assert_eq!(
                flagged,
                ordered,
                "{}: {flagged} bonds flagged aromatic but {ordered} carry the aromatic order",
                format.name()
            );

            // An aromatic bond's atoms must be aromatic too. The reverse does
            // not hold: PDBQT carries no bonds, so its atom flags survive with
            // nothing to ride on -- correct, and why this is one-directional.
            if flagged > 0 {
                assert!(
                    aromatic_atoms > 0,
                    "{}: aromatic bonds but no aromatic atom",
                    format.name()
                );
            }
        }
    }

    /// The formats that can carry a *stated zero* on an atom with no bonds.
    ///
    /// Measured, not assumed. Each has somewhere to say it and a writer that
    /// does: SMILES' `[C]`, a molfile's valence field, commonchem's defaulted
    /// `impHs`. So a count there is read back rather than invented.
    ///
    /// CML is deliberately absent even though it reads `hydrogenCount`: our
    /// writer omits the attribute when the count is zero, so a stated zero does
    /// not survive the round trip. Invisible in SMILES output today, since a
    /// bracketed atom with no `H` is how both `None` and `Some(0)` are written
    /// -- a separate gap from this one.
    ///
    /// Every other format must leave a bondless atom at `None`. Implying a
    /// count there hands it the free-atom valence, because there are no bond
    /// orders to subtract: `[C]` comes back as methane and `[Na+].[Cl-]` as
    /// sodium metal and hydrogen chloride, which is what #282 did to Mol2 and
    /// #291 undid.
    const KEEPS_A_STATED_ZERO: &[Format] = &[
        Format::SMILES,
        Format::CXSMILES,
        Format::SDF,
        Format::COMMONCHEM,
    ];

    #[test]
    fn test_a_count_is_implied_only_where_the_bonds_are_the_whole_story() {
        // Two opposite failures, one rule, and the matrix is blind to both:
        // `held` never reads hydrogens and there is no `Carries` bit for them,
        // so each was green for a whole milestone.
        //
        // Under-filling: PDB stated no count, so every atom kept `None`, the
        // SMILES writer compared 0 against the implied count and bracketed the
        // lot -- ethanol from a PDB read back as `[C][C][O]` (#285).
        //
        // Over-filling: #282 filled every Mol2 atom, and a bondless one has no
        // bond orders to sum, so it took its free-atom valence -- `[C]` came
        // back as methane (#291).
        //
        // Not asserted here, and deliberately: PDBQT and mmCIF are absent from
        // both halves because neither reads back a bonded atom at all today.
        // If either gains one it will land in the first assertion, which is the
        // point -- but note their bond information is partial by construction
        // (PDBQT states only rotatable pivots, `_struct_conn` only specific
        // links), so implying a count from it would invent hydrogens rather
        // than recover them. The fix for whoever gets there is a stated count,
        // not this rule.
        let bonded = crate::io::reader::read("CCO ethanol\n", Format::SMILES);
        let bonded: Vec<(String, Molecule)> = bonded
            .records
            .iter()
            .map(|r| {
                (
                    r.name.clone(),
                    r.molecule()
                        .expect("fixture format is Kind::Molecules")
                        .clone(),
                )
            })
            .collect();
        // Two atoms, no bond between them, each stating it has no hydrogens.
        let lone = crate::io::reader::read("[C].[Cl] lone\n", Format::SMILES);
        let lone: Vec<(String, Molecule)> = lone
            .records
            .iter()
            .map(|r| {
                (
                    r.name.clone(),
                    r.molecule()
                        .expect("fixture format is Kind::Molecules")
                        .clone(),
                )
            })
            .collect();

        for format in all().filter(|f| f.kind() == Kind::Molecules && f.can_read() && f.can_write())
        {
            let bytes = format.write_bytes(&bonded).expect("every format writes");
            let back = format
                .read_bytes(&bytes)
                .expect("every format reads its own bytes");
            let molecule = back.records[0]
                .molecule()
                .expect("fixture format is Kind::Molecules");

            let mut bonded_atoms = 0;
            for index in 0..molecule.num_atoms() {
                if molecule.neighbors(index).is_empty() {
                    continue;
                }
                bonded_atoms += 1;
                assert!(
                    molecule.atom(index).hydrogens().is_some(),
                    "{}: a bonded atom came back with no hydrogen count, so the \
                     writer will bracket it",
                    format.name()
                );
            }
            assert_eq!(
                bonded_atoms > 0,
                format.carries().contains(Carries::BONDS),
                "{}: bonds read back disagree with the mask",
                format.name()
            );

            let bytes = format.write_bytes(&lone).expect("every format writes");
            let back = format
                .read_bytes(&bytes)
                .expect("every format reads its own bytes");
            let molecule = back.records[0]
                .molecule()
                .expect("fixture format is Kind::Molecules");

            for index in 0..molecule.num_atoms() {
                assert!(
                    molecule.neighbors(index).is_empty(),
                    "{}: the bondless fixture came back with a bond",
                    format.name()
                );
                assert_eq!(
                    molecule.atom(index).hydrogens().is_some(),
                    KEEPS_A_STATED_ZERO.contains(&format),
                    "{}: a bondless atom's count was invented rather than read",
                    format.name()
                );
            }
        }
    }

    #[test]
    fn test_every_format_writes_as_many_records_as_it_was_given() {
        // The dimension the matrix was blind to. Every fixture in
        // `one_per_attribute` is a single molecule and `convert` reads back
        // `.records.first()`, so a writer that merged N records into one --
        // which PDB and PDBQT both did, for want of `MODEL`/`ENDMDL` framing --
        // satisfied every assertion in this module (#267).
        //
        // Written as the general rule rather than a PDB test, because nobody
        // thought to check PDB specifically for a whole milestone.
        let outcome = crate::io::reader::read("CCO a\nc1ccccc1 b\nCCN c\n", Format::SMILES);
        let records: Vec<(String, Molecule)> = outcome
            .records
            .iter()
            .map(|r| {
                (
                    r.name.clone(),
                    r.molecule()
                        .expect("fixture format is Kind::Molecules")
                        .clone(),
                )
            })
            .collect();
        assert_eq!(records.len(), 3, "the fixture must be multi-record");

        for format in all().filter(|f| f.kind() == Kind::Molecules) {
            let bytes = format.write_bytes(&records).expect("every format writes");
            let back = format
                .read_bytes(&bytes)
                .expect("every format reads its own bytes");
            assert_eq!(
                back.records.len(),
                records.len(),
                "{} wrote {} records and read back {}",
                format.name(),
                records.len(),
                back.records.len()
            );
        }
    }

    #[test]
    fn test_no_format_pair_silently_loses_atoms() {
        // The assertion `Carries` cannot make. `held` sets TOPOLOGY on the
        // atom count alone, so a conversion keeping one atom of six claims
        // everything the mask promised -- which is how PDBQT's
        // largest-component rule stayed invisible through a whole milestone,
        // on its own diagonal included.
        let fixtures = one_per_attribute();
        let mut found: Vec<(Format, Format)> = Vec::new();

        for (source, target) in pairs_within_kind() {
            let lost = fixtures.iter().any(|(_, molecule)| {
                convert(source, target, molecule)
                    .is_none_or(|back| back.num_atoms() != molecule.num_atoms())
            });
            if lost {
                found.push((source, target));
            }
        }

        for (source, target) in &found {
            assert!(
                pair_gap(*source, *target).is_some(),
                "{} -> {} loses atoms and is not pinned in PAIR_GAPS",
                source.name(),
                target.name()
            );
        }
        for (source, target, why) in PAIR_GAPS {
            assert!(
                found.contains(&(*source, *target)),
                "{} -> {} is pinned as {why:?} but no longer loses atoms -- delete the line",
                source.name(),
                target.name()
            );
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
            // An empty mask means the entry was added without a defining
            // flag, which would silently report every record as losing
            // everything. Every kind gets exactly one mandatory "you exist"
            // flag (#316) -- TOPOLOGY for a molecule, SAMPLES for a grid,
            // VERTICES for a mesh, COLUMNS for a table -- mirroring how
            // BONDS/FACES are separate, optional facts once the mandatory
            // one is satisfied.
            match format.kind() {
                Kind::Molecules | Kind::Frames => {
                    assert!(
                        format.carries().contains(Carries::TOPOLOGY),
                        "{format:?} declares no topology, so its mask is missing"
                    );
                }
                Kind::Volume => {
                    assert!(
                        format.carries().contains(Carries::SAMPLES),
                        "{format:?} declares no samples, so its mask is missing"
                    );
                }
                Kind::Mesh => {
                    assert!(
                        format.carries().contains(Carries::VERTICES),
                        "{format:?} declares no vertices, so its mask is missing"
                    );
                }
                Kind::Table => {
                    assert!(
                        format.carries().contains(Carries::COLUMNS),
                        "{format:?} declares no columns, so its mask is missing"
                    );
                }
            }
        }
    }

    #[test]
    fn test_the_topology_invariant_is_kind_aware_not_disabled() {
        // `test_every_registered_format_is_well_formed`'s relaxation above,
        // exercised against real `FormatDescriptor` values that never enter
        // `FORMATS` -- none of `Kind::Volume`/`Mesh`/`Table` has a real
        // container type yet (#311-#314), so this is the only way to prove
        // the relaxation works for the kinds it targets without waiting on
        // them. Same condition as the real test, applied directly to a
        // descriptor's own fields rather than duplicated in a helper.
        fn bare(kind: Kind) -> FormatDescriptor {
            FormatDescriptor {
                name: "probe",
                codes: &["probe"],
                extensions: &["probe"],
                category: Category::Miscellaneous,
                carries: Carries::empty(),
                encoding: Encoding::Text,
                kind,
                magic: &[],
                reader: None,
                writer: None,
                reader_bytes: None,
                writer_bytes: None,
                supplier: None,
                writer_stream: None,
                writer_trajectory: None,
                writer_volume: None,
                writer_mesh: None,
                writer_table: None,
            }
        }

        let requires_topology =
            |d: &FormatDescriptor| matches!(d.kind, Kind::Molecules | Kind::Frames);

        // A volume descriptor with nothing declared: the relaxed check
        // accepts it -- an empty mask is meaningless for this kind.
        assert!(!requires_topology(&bare(Kind::Volume)));
        assert!(!requires_topology(&bare(Kind::Mesh)));
        assert!(!requires_topology(&bare(Kind::Table)));
        // A molecule/frames descriptor with nothing declared: still flagged
        // -- the relaxation did not turn the check off for the kinds that
        // need it.
        assert!(requires_topology(&bare(Kind::Molecules)));
        assert!(requires_topology(&bare(Kind::Frames)));
    }

    #[test]
    fn test_every_kind_requires_its_own_defining_flag() {
        // #316's half of the same proof, for the three kinds that gained a
        // real container type (#312-#314) but no registered format yet --
        // SAMPLES/VERTICES/COLUMNS are each kind's own "you exist" flag, the
        // role TOPOLOGY plays for a molecule. The predicate here is written
        // out independently of the real match in
        // `test_every_registered_format_is_well_formed` rather than sharing
        // it, the same discipline `requires_topology` above already follows
        // -- a test that calls the function it is meant to catch bugs in
        // proves nothing.
        fn bare(kind: Kind, carries: Carries) -> FormatDescriptor {
            FormatDescriptor {
                name: "probe",
                codes: &["probe"],
                extensions: &["probe"],
                category: Category::Miscellaneous,
                carries,
                encoding: Encoding::Text,
                kind,
                magic: &[],
                reader: None,
                writer: None,
                reader_bytes: None,
                writer_bytes: None,
                supplier: None,
                writer_stream: None,
                writer_trajectory: None,
                writer_volume: None,
                writer_mesh: None,
                writer_table: None,
            }
        }

        fn passes(d: &FormatDescriptor) -> bool {
            match d.kind {
                Kind::Molecules | Kind::Frames => d.carries.contains(Carries::TOPOLOGY),
                Kind::Volume => d.carries.contains(Carries::SAMPLES),
                Kind::Mesh => d.carries.contains(Carries::VERTICES),
                Kind::Table => d.carries.contains(Carries::COLUMNS),
            }
        }

        // Each kind's own defining flag satisfies it ...
        assert!(passes(&bare(Kind::Volume, Carries::SAMPLES)));
        assert!(passes(&bare(Kind::Mesh, Carries::VERTICES)));
        assert!(passes(&bare(Kind::Table, Carries::COLUMNS)));

        // ... an empty mask does not -- the loophole the issue opens with.
        assert!(!passes(&bare(Kind::Volume, Carries::empty())));
        assert!(!passes(&bare(Kind::Mesh, Carries::empty())));
        assert!(!passes(&bare(Kind::Table, Carries::empty())));

        // ... and declaring a different kind's defining flag does not
        // satisfy this one -- the flags are not interchangeable.
        assert!(!passes(&bare(Kind::Volume, Carries::VERTICES)));
        assert!(!passes(&bare(Kind::Mesh, Carries::COLUMNS)));
        assert!(!passes(&bare(Kind::Table, Carries::SAMPLES)));
    }

    #[test]
    fn test_has_magic_at_matches_and_never_panics_on_short_input() {
        assert!(has_magic_at(b"\x1f\x8bxxxx", 0, b"\x1f\x8b"));
        assert!(!has_magic_at(b"nope", 0, b"\x1f\x8b"));
        // CCP4's `MAP ` sits at byte 208, not byte 0 -- a signature is not
        // always anchored at the start of the file.
        let mut buf = vec![0u8; 212];
        buf[208..212].copy_from_slice(b"MAP ");
        assert!(has_magic_at(&buf, 208, b"MAP "));
        // Too short to hold the pattern at that offset at all -- must not
        // panic, only answer false.
        assert!(!has_magic_at(b"MA", 208, b"MAP "));
        assert!(!has_magic_at(b"", 0, b"\x1f\x8b"));
    }

    #[test]
    fn test_find_signature_matches_the_right_candidate() {
        // `Format` is always a real, registered handle (see its own doc
        // comment), so there is no way to hand this a fake one -- pairing a
        // real `Format` with a made-up signature list proves the same
        // matching/dispatch logic `sniff` uses without waiting on a real
        // binary format to register one (#317).
        let candidates: Vec<(Format, &'static [Signature])> = vec![
            (
                Format::SMILES,
                &[Signature {
                    offset: 0,
                    bytes: b"AA",
                }],
            ),
            (
                Format::SDF,
                &[Signature {
                    offset: 2,
                    bytes: b"BB",
                }],
            ),
        ];

        assert_eq!(
            find_signature(b"AAxx", candidates.iter().copied()),
            Some(Format::SMILES)
        );
        assert_eq!(
            find_signature(b"xxBB", candidates.iter().copied()),
            Some(Format::SDF)
        );
        assert_eq!(find_signature(b"zzzz", candidates.iter().copied()), None);
    }

    #[test]
    fn test_sniff_resolves_every_binary_formats_real_signature_and_nothing_else_has_one_yet() {
        // TRR (#326) was the first format to populate `magic` for real --
        // #309 built the byte-reading path a binary format needs, and
        // #317 built this matching mechanism -- XTC (#325) is the second,
        // DCD (#327) the third, NCTRAJ (#328) the fourth, CCP4 (#332) the
        // fifth -- the first outside `Kind::Frames` -- and PLY (#336) the
        // sixth -- the first outside `Kind::Volume` too, and the first
        // `Kind::Mesh` format with a real fixed signature at all (OBJ has
        // none). Every other format leaves `magic` empty. Pinned explicitly
        // rather than trusted silently, the same discipline #316's
        // `pairs_within_kind` count assertion follows.
        for format in all() {
            if format == Format::TRR
                || format == Format::XTC
                || format == Format::DCD
                || format == Format::NCTRAJ
                || format == Format::CCP4
                || format == Format::PLY
            {
                continue;
            }
            assert!(
                format.descriptor().magic.is_empty(),
                "{format:?} already has a signature"
            );
        }
        assert!(sniff(b"\x1f\x8b\x08\x00").is_none());
        assert!(sniff(b"CORD").is_none());
        assert!(sniff(b"").is_none());

        // GROMACS's own fixed magic numbers, big-endian, 2 apart.
        assert_eq!(sniff(b"\x00\x00\x07\xc9REST"), Some(Format::TRR));
        assert_eq!(sniff(b"\x00\x00\x07\xcbREST"), Some(Format::XTC));
        // DCD's "CORD" sits at byte 4, after the leading Fortran record
        // marker (whatever it is) -- not byte 0.
        assert_eq!(sniff(b"\x54\x00\x00\x00CORD"), Some(Format::DCD));
        // Every NetCDF-3 file, classic or 64-bit offset -- this alone
        // doesn't prove it's an Amber trajectory, just a NetCDF-3 file.
        assert_eq!(sniff(b"CDF\x01REST"), Some(Format::NCTRAJ));
        assert_eq!(sniff(b"CDF\x02REST"), Some(Format::NCTRAJ));
        // CCP4/MRC's "MAP " sits at byte 208, not byte 0.
        let mut ccp4_like = vec![0u8; 212];
        ccp4_like[208..212].copy_from_slice(b"MAP ");
        assert_eq!(sniff(&ccp4_like), Some(Format::CCP4));
        // PLY's "ply\n" is byte-identical across all three of its wire
        // encodings, at byte 0.
        assert_eq!(sniff(b"ply\nformat ascii 1.0\n"), Some(Format::PLY));
    }

    #[test]
    fn test_every_text_format_stays_text_encoded() {
        // #309 added the byte-level path; #319 (BinaryCIF) and #326 (TRR)
        // are the formats that actually use it -- skipped here by their
        // `Encoding`, not by name, since every *other* format must still be
        // plain text with no binary reader/writer wired in.
        for format in all() {
            if format.encoding() == Encoding::Binary {
                continue;
            }
            let d = format.descriptor();
            assert_eq!(
                d.encoding,
                Encoding::Text,
                "{format:?} is not registered as text"
            );
            assert!(
                d.reader_bytes.is_none(),
                "{format:?} has a byte reader already"
            );
            assert!(
                d.writer_bytes.is_none(),
                "{format:?} has a byte writer already"
            );
        }
    }

    #[test]
    fn test_read_bytes_agrees_with_read() {
        // Proves `read`/`read_with_options` genuinely delegate to
        // `read_bytes_with_options` (#309) rather than sitting next to it as
        // dead code: the same input through either path must yield the same
        // records.
        let sdf = "ethanol-ish\n  -ish-\n\nM  END\n$$$$\n";
        let via_read = crate::io::reader::read(sdf, Format::SDF);
        let via_bytes = Format::SDF
            .read_bytes(sdf.as_bytes())
            .expect("SDF can be read");
        assert_eq!(via_read.records.len(), via_bytes.records.len());
        assert_eq!(via_read.skipped.len(), via_bytes.skipped.len());
    }

    #[test]
    fn test_read_bytes_reports_invalid_utf8_as_skipped_not_a_panic() {
        let invalid = [b'C', 0xff, 0xfe];
        let outcome = Format::SMILES
            .read_bytes(&invalid)
            .expect("SMILES can be read");
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(outcome.skipped[0].error.contains("UTF-8"));
    }

    #[test]
    fn test_write_bytes_agrees_with_write() {
        use crate::core::prelude::*;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        let records = vec![("m".to_string(), mol)];

        let via_write = Format::SDF.write(&records).expect("SDF writes");
        let via_bytes = Format::SDF.write_bytes(&records).expect("SDF writes bytes");
        assert_eq!(via_write.into_bytes(), via_bytes);
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
        assert_eq!(all().count(), 29);
        for format in all() {
            assert!(format.can_read() && format.can_write(), "{format:?}");
        }
    }

    #[test]
    fn test_kinds_compatible_allows_a_same_kind_pair() {
        assert!(kinds_compatible(Format::TRR, Format::XTC).is_ok());
        assert!(kinds_compatible(Format::CCP4, Format::DX).is_ok());
        assert!(kinds_compatible(Format::OBJ, Format::PLY).is_ok());
        assert!(kinds_compatible(Format::CSV, Format::CSV).is_ok());
    }

    #[test]
    fn test_kinds_compatible_allows_cubes_own_dual_nature_one_direction_only() {
        // CUBE genuinely carries atoms (Carries::TOPOLOGY) alongside its
        // grid -- the one deliberately enumerated exception.
        assert!(kinds_compatible(Format::CUBE, Format::PDB).is_ok());
        // Never the reverse: a molecule has no grid samples to offer.
        let err = kinds_compatible(Format::PDB, Format::CUBE).unwrap_err();
        assert!(err.contains("PDB"), "{err}");
        assert!(err.contains("CUBE"), "{err}");
    }

    #[test]
    fn test_kinds_compatible_refuses_a_genuine_cross_kind_pair_naming_both() {
        let err = kinds_compatible(Format::CSV, Format::DCD).unwrap_err();
        assert!(err.contains("CSV"), "{err}");
        assert!(err.contains("DCD"), "{err}");
        assert!(err.contains("table"), "{err}");
        assert!(err.contains("trajectory"), "{err}");
    }

    #[test]
    fn test_kinds_compatible_refuses_a_volume_format_with_no_atoms_into_molecules() {
        // CCP4 is Kind::Volume but declares no Carries::TOPOLOGY -- unlike
        // CUBE, it never has real atoms to offer, so it gets no exception.
        let err = kinds_compatible(Format::CCP4, Format::PDB).unwrap_err();
        assert!(err.contains("CCP4"), "{err}");
        assert!(err.contains("PDB"), "{err}");
    }

    #[test]
    fn test_held_from_volume_reports_only_what_this_specific_grid_states() {
        use crate::core::atom::{Atom, Element};
        use crate::core::cell::UnitCell;
        use crate::core::geometry::Point3;
        use crate::core::molecule::Molecule;
        use crate::core::volume::VolumeGrid;

        // No cell, no atoms -- DX's own shape.
        let bare = VolumeGrid::new(
            [1, 1, 1],
            Point3::ORIGIN,
            [
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            vec![0.0],
            None,
        )
        .unwrap();
        assert_eq!(held_from_volume(&bare), Carries::SAMPLES);

        // A cell but no atoms -- CCP4/DSN6's own shape.
        let with_cell = VolumeGrid::new(
            [1, 1, 1],
            Point3::ORIGIN,
            [
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            vec![0.0],
            Some(UnitCell::cubic(10.0)),
        )
        .unwrap();
        assert_eq!(
            held_from_volume(&with_cell),
            Carries::SAMPLES.or(Carries::UNIT_CELL)
        );

        // Atoms but no cell -- CUBE's own shape.
        let mut with_atoms = bare.clone();
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        with_atoms.set_atoms(mol);
        assert_eq!(
            held_from_volume(&with_atoms),
            Carries::SAMPLES
                .or(Carries::TOPOLOGY)
                .or(Carries::COORDS_3D)
        );
    }

    #[test]
    fn test_held_from_trajectory_reports_only_what_frame_zero_states() {
        use crate::core::atom::{Atom, Element};
        use crate::core::geometry::Point3;
        use crate::core::molecule::Molecule;
        use crate::core::trajectory::{Frame, FrameSource};

        struct OneFrame(Frame);
        impl FrameSource for OneFrame {
            fn frame_count(&self) -> usize {
                1
            }
            fn num_atoms(&self) -> usize {
                self.0.num_atoms()
            }
            fn frame(&mut self, _index: usize) -> std::io::Result<Frame> {
                Ok(self.0.clone())
            }
        }

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));

        // XTC's own shape: positions only.
        let bare_frame = Frame {
            positions: vec![Point3::ORIGIN],
            velocities: None,
            forces: None,
            time: None,
            step: None,
            cell: None,
        };
        let mut bare = Trajectory::new(mol.clone(), Box::new(OneFrame(bare_frame))).unwrap();
        assert_eq!(
            held_from_trajectory(&mut bare).unwrap(),
            Carries::TOPOLOGY.or(Carries::COORDS_3D)
        );

        // TRR's own shape: velocities and forces too.
        let full_frame = Frame {
            positions: vec![Point3::ORIGIN],
            velocities: Some(vec![Point3::ORIGIN]),
            forces: Some(vec![Point3::ORIGIN]),
            time: Some(0.0),
            step: None,
            cell: None,
        };
        let mut full = Trajectory::new(mol, Box::new(OneFrame(full_frame))).unwrap();
        assert_eq!(
            held_from_trajectory(&mut full).unwrap(),
            Carries::TOPOLOGY
                .or(Carries::COORDS_3D)
                .or(Carries::VELOCITIES)
                .or(Carries::FORCES)
                .or(Carries::FRAME_TIME)
        );
    }
}
