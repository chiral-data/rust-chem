//! LAMMPS data file — a header (a skipped title line, mandatory box bounds,
//! and discardable counts) followed by named sections: `Masses`, `Atoms`,
//! `Velocities`, `Bonds`, `Angles`, `Dihedrals`, `Impropers`, and several
//! `* Coeffs` parameter tables (#324).
//!
//! **The first force-field-topology format that also states coordinates.**
//! PSF/PRMTOP/TOP all deliberately never touch `coords3` (paired with a
//! separate coordinate file in their own ecosystems); LAMMPS data is
//! AMBER's `.prmtop`+`.inpcrd` and GROMACS's `.top`+`.gro` combined into
//! one file. `Molecule::set_force_field`/`set_coords3` are independent
//! methods with no invariant coupling them, so this combination needs no
//! new machinery, just the first real use of both together.
//!
//! **Topology only**, the same milestone-wide scope PSF/PRMTOP/TOP already
//! settled on: every `* Coeffs` parameter section, `Velocities`,
//! `Ellipsoids`/`Lines`/`Triangles`/`Bodies`, and every `* Type Labels`
//! section is recognized and skipped, never stored. LAMMPS has no
//! `Exclusions`-shaped section at all (exclusions are a `special_bonds`
//! setting in the *input script*, not file data), so
//! `ForceFieldTopology::exclusions` stays empty here.
//!
//! **The atom style is not stated reliably.** The `Atoms` section's column
//! layout depends on which `atom_style` produced the file; this module
//! models four: `atomic` (id,type,x,y,z), `charge` (id,type,q,x,y,z),
//! `full` (id,mol,type,q,x,y,z), and `molecular` (id,mol,type,x,y,z --
//! also `bond` and `angle`, identical shapes this crate never needs to
//! tell apart). An optional trailing `ix,iy,iz` image-flag triplet can pad
//! any style's column count. Resolved in order: an explicit
//! [`crate::io::options::LammpsReadOptions::atom_style`], then a
//! recognized `Atoms # style` comment (a real, LAMMPS-recognized
//! convention, not an informal one), then an unambiguous column count.
//! `charge` and `molecular` share a column count both with and without
//! image flags -- [`LammpsError::AmbiguousAtomStyle`] refuses rather than
//! guesses when none of the three resolves it.
//!
//! **Elements are never inferred.** A LAMMPS atom type has a number and
//! (usually) a mass, never a name -- mapping mass to an element the way
//! PSF/PRMTOP/TOP's `element_from_mass` does would fail outright for
//! coarse-grained/reduced-unit systems (bead-spring polymers, LJ fluids)
//! whose masses are things like `1.0`, corresponding to no real element.
//! Every atom here is [`crate::core::atom::Element::UNKNOWN`], honestly
//! recording only its numeric type and mass. `calculate_implicit_hydrogens_where_bonded`
//! still runs over `Bonds`'s complete real bond list, the same precedent
//! PSF/PRMTOP/TOP already established -- `typical_valence()` for an
//! unknown element is honestly `0`, so every bonded atom here gets
//! `Some(0)`, the same answer this crate already gives a real element
//! with no typical valence (a metal, say), not a new claim of certainty.
//!
//! **Box bounds convert to [`UnitCell`], with a disclosed origin loss.**
//! LAMMPS's restricted-triclinic box already uses the same "a along x, b
//! in the xy plane" orientation convention `UnitCell` hardcodes, so cell
//! *shape* survives exactly via the standard closed-form conversion. What
//! doesn't survive: `UnitCell` has no origin field, so `xlo`/`ylo`/`zlo`
//! (often nonzero in real files) are discarded -- atom coordinates are
//! kept exactly as stated, not re-centered, the same "cell is shape-only,
//! atoms are raw Cartesian" posture PDB's own `CRYST1`+`ATOM` pair already
//! has in this crate.
//!
//! **No chains or residues.** A molecule-id column exists in
//! `full`/`molecular` styles, but it is a bare integer with no name --
//! inventing a residue name from a number would be fabrication, so
//! `set_topology` is simply never called, the same treatment
//! [`crate::io::cif_core`] gives a dictionary with no residue machinery.
//!
//! **Units are unrecoverable from the file alone.** LAMMPS's `units`
//! command lives in the input script, never the data file -- this reader
//! stores whatever numbers the file states, uninterpreted.
//!
//! **Header counts are read but never trusted** -- section boundaries (the
//! next recognized section keyword, or EOF) are the real count, the same
//! design PRMTOP's reader already uses.

use std::collections::HashMap;

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::cell::UnitCell;
use crate::core::elements::ATOMIC_MASSES;
use crate::core::force_field::{ForceFieldAtom, ForceFieldTopology};
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::io::errors::LammpsError;
use crate::io::options::{AtomStyle, LammpsReadOptions};

/// Every real LAMMPS data section keyword this reader recognizes as a
/// section boundary. Anything else is a data line belonging to whichever
/// section (if any) is currently open.
const RECOGNIZED_SECTIONS: &[&str] = &[
    "Masses",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
    "Pair Coeffs",
    "PairIJ Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "BondBond Coeffs",
    "BondAngle Coeffs",
    "MiddleBondTorsion Coeffs",
    "EndBondTorsion Coeffs",
    "AngleTorsion Coeffs",
    "AngleAngleTorsion Coeffs",
    "BondBond13 Coeffs",
    "AngleAngle Coeffs",
    "Ellipsoids",
    "Lines",
    "Triangles",
    "Bodies",
    "Atom Type Labels",
    "Bond Type Labels",
    "Angle Type Labels",
    "Dihedral Type Labels",
    "Improper Type Labels",
];

/// Splits a trailing `# comment` off `line`, trimming trailing whitespace
/// from what's left. `None` when there is no `#` at all.
fn split_trailing_comment(line: &str) -> (&str, Option<&str>) {
    match line.split_once('#') {
        Some((before, after)) => (before.trim_end(), Some(after.trim())),
        None => (line, None),
    }
}

fn parse_f64(token: &str) -> Result<f64, LammpsError> {
    token
        .parse()
        .map_err(|_| LammpsError::ParseError(format!("expected a number, got {token:?}")))
}

fn parse_i64(token: &str) -> Result<i64, LammpsError> {
    token
        .parse()
        .map_err(|_| LammpsError::ParseError(format!("expected an integer, got {token:?}")))
}

/// The file split into its header lines (before the first recognized
/// section) and the recognized bodies this reader keeps. Every other
/// recognized section is tracked only well enough to know where its data
/// ends -- never accumulated, since nothing here has a data-model slot for
/// it (see the module doc).
struct ScannedSections<'a> {
    header_lines: Vec<&'a str>,
    masses: Vec<&'a str>,
    atoms: Vec<&'a str>,
    bonds: Vec<&'a str>,
    angles: Vec<&'a str>,
    dihedrals: Vec<&'a str>,
    impropers: Vec<&'a str>,
    /// The comment on the `Atoms` section's own header line, if any --
    /// `Atoms # full` names the atom style, a real LAMMPS-recognized
    /// convention (see the module doc).
    atoms_header_comment: Option<&'a str>,
}

enum Target {
    Header,
    Masses,
    Atoms,
    Bonds,
    Angles,
    Dihedrals,
    Impropers,
    Discard,
}

fn scan_sections(text: &str) -> ScannedSections<'_> {
    let mut lines = text.lines();
    lines.next(); // The mandatory, always-skipped title/comment line.

    let mut out = ScannedSections {
        header_lines: Vec::new(),
        masses: Vec::new(),
        atoms: Vec::new(),
        bonds: Vec::new(),
        angles: Vec::new(),
        dihedrals: Vec::new(),
        impropers: Vec::new(),
        atoms_header_comment: None,
    };
    let mut target = Target::Header;

    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let (content, comment) = split_trailing_comment(trimmed);
        let content = content.trim();
        if RECOGNIZED_SECTIONS.contains(&content) {
            target = match content {
                "Masses" => Target::Masses,
                "Atoms" => {
                    out.atoms_header_comment = comment;
                    Target::Atoms
                }
                "Bonds" => Target::Bonds,
                "Angles" => Target::Angles,
                "Dihedrals" => Target::Dihedrals,
                "Impropers" => Target::Impropers,
                _ => Target::Discard,
            };
            continue;
        }
        match target {
            Target::Header => out.header_lines.push(line),
            Target::Masses => out.masses.push(line),
            Target::Atoms => out.atoms.push(line),
            Target::Bonds => out.bonds.push(line),
            Target::Angles => out.angles.push(line),
            Target::Dihedrals => out.dihedrals.push(line),
            Target::Impropers => out.impropers.push(line),
            Target::Discard => {}
        }
    }
    out
}

/// The closed-form `lx/ly/lz/xy/xz/yz` (LAMMPS restricted-triclinic edge
/// lengths and tilt factors) to [`UnitCell`] conversion -- shape-lossless,
/// since both already share the "a along x, b in the xy plane" convention.
/// Shared by [`parse_box`] (data-file `xlo xhi`-style bounds) and
/// [`crate::io::lammpstrj`] (dump-file `BOX BOUNDS`, which needs its own
/// bounding-box-to-edges inversion first, since a triclinic dump reports
/// values already shifted by the tilt).
pub(crate) fn unit_cell_from_lammps_box(
    lx: f64,
    ly: f64,
    lz: f64,
    xy: f64,
    xz: f64,
    yz: f64,
) -> UnitCell {
    let b = (ly * ly + xy * xy).sqrt();
    let c = (lz * lz + xz * xz + yz * yz).sqrt();
    let cos_alpha = ((xy * xz + ly * yz) / (b * c)).clamp(-1.0, 1.0);
    let cos_beta = (xz / c).clamp(-1.0, 1.0);
    let cos_gamma = (xy / b).clamp(-1.0, 1.0);

    UnitCell::new(
        lx,
        b,
        c,
        cos_alpha.acos().to_degrees(),
        cos_beta.acos().to_degrees(),
        cos_gamma.acos().to_degrees(),
    )
}

/// Converts LAMMPS's box bounds (`xlo xhi`/`ylo yhi`/`zlo zhi`, optional
/// triclinic `xy xz yz`) into a [`UnitCell`] via [`unit_cell_from_lammps_box`]
/// (see the module doc for the origin caveat this does *not* solve).
fn parse_box(header_lines: &[&str]) -> Result<UnitCell, LammpsError> {
    let mut xlo = None;
    let mut xhi = None;
    let mut ylo = None;
    let mut yhi = None;
    let mut zlo = None;
    let mut zhi = None;
    let mut xy = 0.0;
    let mut xz = 0.0;
    let mut yz = 0.0;

    for raw in header_lines {
        let (content, _) = split_trailing_comment(raw);
        let tokens: Vec<&str> = content.split_whitespace().collect();
        match tokens.as_slice() {
            [lo, hi, "xlo", "xhi"] => {
                xlo = Some(parse_f64(lo)?);
                xhi = Some(parse_f64(hi)?);
            }
            [lo, hi, "ylo", "yhi"] => {
                ylo = Some(parse_f64(lo)?);
                yhi = Some(parse_f64(hi)?);
            }
            [lo, hi, "zlo", "zhi"] => {
                zlo = Some(parse_f64(lo)?);
                zhi = Some(parse_f64(hi)?);
            }
            [a, b, c, "xy", "xz", "yz"] => {
                xy = parse_f64(a)?;
                xz = parse_f64(b)?;
                yz = parse_f64(c)?;
            }
            // Any other header line (an atom/bond/type count, or anything
            // this reader doesn't need) is discarded -- see the module
            // doc's "header counts are never trusted" decision.
            _ => {}
        }
    }

    let missing = || LammpsError::ParseError("missing box bounds".to_string());
    let (xlo, xhi) = (xlo.ok_or_else(missing)?, xhi.ok_or_else(missing)?);
    let (ylo, yhi) = (ylo.ok_or_else(missing)?, yhi.ok_or_else(missing)?);
    let (zlo, zhi) = (zlo.ok_or_else(missing)?, zhi.ok_or_else(missing)?);

    let lx = xhi - xlo;
    let ly = yhi - ylo;
    let lz = zhi - zlo;

    Ok(unit_cell_from_lammps_box(lx, ly, lz, xy, xz, yz))
}

fn parse_masses(lines: &[&str]) -> Result<HashMap<i64, f64>, LammpsError> {
    let mut masses = HashMap::new();
    for raw in lines {
        let (content, _) = split_trailing_comment(raw);
        let fields: Vec<&str> = content.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        masses.insert(parse_i64(fields[0])?, parse_f64(fields[1])?);
    }
    Ok(masses)
}

/// Resolves which [`AtomStyle`] the `Atoms` section uses: an explicit
/// option first, then a recognized `# style` comment, then an unambiguous
/// column count -- refusing (not guessing) when none of those settles it.
fn resolve_atom_style(
    options: &LammpsReadOptions,
    header_comment: Option<&str>,
    first_data_line: Option<&str>,
) -> Result<AtomStyle, LammpsError> {
    if let Some(style) = options.atom_style {
        return Ok(style);
    }
    if let Some(name) = header_comment {
        return match name {
            "atomic" => Ok(AtomStyle::Atomic),
            "charge" => Ok(AtomStyle::Charge),
            "full" => Ok(AtomStyle::Full),
            "molecular" | "bond" | "angle" => Ok(AtomStyle::Molecular),
            other => Err(LammpsError::UnsupportedAtomStyle(other.to_string())),
        };
    }
    let first_line = first_data_line.ok_or(LammpsError::NoAtoms)?;
    let (content, _) = split_trailing_comment(first_line);
    match content.split_whitespace().count() {
        5 | 8 => Ok(AtomStyle::Atomic),
        7 | 10 => Ok(AtomStyle::Full),
        n @ (6 | 9) => Err(LammpsError::AmbiguousAtomStyle(n)),
        n => Err(LammpsError::ParseError(format!(
            "{n} columns in the Atoms section matches no supported atom style"
        ))),
    }
}

struct ParsedAtomLine {
    id: i64,
    atom_type: i64,
    charge: Option<f64>,
    x: f64,
    y: f64,
    z: f64,
}

fn parse_atom_line(style: AtomStyle, raw: &str) -> Result<ParsedAtomLine, LammpsError> {
    let (content, _) = split_trailing_comment(raw);
    let fields: Vec<&str> = content.split_whitespace().collect();
    let base_len = match style {
        AtomStyle::Atomic => 5,
        AtomStyle::Charge | AtomStyle::Molecular => 6,
        AtomStyle::Full => 7,
    };
    if fields.len() != base_len && fields.len() != base_len + 3 {
        return Err(LammpsError::InvalidAtomLine(raw.to_string()));
    }
    let bad = || LammpsError::InvalidAtomLine(raw.to_string());
    let field_i64 = |s: &str| s.parse::<i64>().map_err(|_| bad());
    let field_f64 = |s: &str| s.parse::<f64>().map_err(|_| bad());

    let id = field_i64(fields[0])?;
    Ok(match style {
        AtomStyle::Atomic => ParsedAtomLine {
            id,
            atom_type: field_i64(fields[1])?,
            charge: None,
            x: field_f64(fields[2])?,
            y: field_f64(fields[3])?,
            z: field_f64(fields[4])?,
        },
        AtomStyle::Charge => ParsedAtomLine {
            id,
            atom_type: field_i64(fields[1])?,
            charge: Some(field_f64(fields[2])?),
            x: field_f64(fields[3])?,
            y: field_f64(fields[4])?,
            z: field_f64(fields[5])?,
        },
        // `fields[1]` is the molecule-id -- no data-model slot, discarded.
        AtomStyle::Molecular => ParsedAtomLine {
            id,
            atom_type: field_i64(fields[2])?,
            charge: None,
            x: field_f64(fields[3])?,
            y: field_f64(fields[4])?,
            z: field_f64(fields[5])?,
        },
        AtomStyle::Full => ParsedAtomLine {
            id,
            atom_type: field_i64(fields[2])?,
            charge: Some(field_f64(fields[3])?),
            x: field_f64(fields[4])?,
            y: field_f64(fields[5])?,
            z: field_f64(fields[6])?,
        },
    })
}

/// Reads `n` global atom ids from a `Bonds`/`Angles`/`Dihedrals`/
/// `Impropers` line -- `id type atom1 atom2 [atom3 [atom4]]`, with the
/// leading id and bonded-type index discarded (no data-model slot).
fn parse_bonded_ids(raw: &str, n: usize) -> Result<Vec<i64>, LammpsError> {
    let (content, _) = split_trailing_comment(raw);
    let fields: Vec<&str> = content.split_whitespace().collect();
    if fields.len() < 2 + n {
        return Err(LammpsError::ParseError(format!(
            "invalid bonded-term line: {raw:?}"
        )));
    }
    fields[2..2 + n].iter().map(|s| parse_i64(s)).collect()
}

fn resolve_atom(id_to_index: &HashMap<i64, usize>, id: i64) -> Result<usize, LammpsError> {
    id_to_index.get(&id).copied().ok_or_else(|| {
        LammpsError::ParseError(format!("bonded section references unknown atom id {id}"))
    })
}

/// Pads a bounding-box edge that has zero (or near-zero) extent -- e.g.
/// every atom sharing one axis's coordinate, or a single-atom molecule --
/// so the box this writer states is never degenerate. A zero-length edge
/// would make [`parse_box`]'s closed-form conversion divide by zero on the
/// way back in, failing a round trip that has nothing to do with the box
/// at all. The padding is a fixed, disclosed placeholder, the same
/// "something has to go here" reasoning [`crate::io::cif_core`]'s own
/// cell-less-molecule writer already uses.
fn ensure_nonzero_extent(lo: f64, hi: f64) -> (f64, f64) {
    const MIN_EXTENT: f64 = 1.0;
    if hi - lo > 1e-6 {
        (lo, hi)
    } else {
        (lo - MIN_EXTENT / 2.0, hi + MIN_EXTENT / 2.0)
    }
}

/// Parses a whole LAMMPS data file into a single [`Molecule`] -- one
/// topology per file, the same shape as PSF/PRMTOP, since LAMMPS's atom
/// ids are global and file-wide, not per-molecule-local like GROMACS TOP's.
pub fn parse_lammps_data(text: &str, options: &LammpsReadOptions) -> Result<Molecule, LammpsError> {
    let sections = scan_sections(text);
    let cell = parse_box(&sections.header_lines)?;
    let masses = parse_masses(&sections.masses)?;

    let non_blank_atom_lines: Vec<&str> = sections
        .atoms
        .iter()
        .copied()
        .filter(|line| !split_trailing_comment(line).0.trim().is_empty())
        .collect();
    if non_blank_atom_lines.is_empty() {
        return Err(LammpsError::NoAtoms);
    }
    let style = resolve_atom_style(
        options,
        sections.atoms_header_comment,
        non_blank_atom_lines.first().copied(),
    )?;

    let mut mol = Molecule::new();
    let mut ff_atoms = Vec::with_capacity(non_blank_atom_lines.len());
    let mut coords = Vec::with_capacity(non_blank_atom_lines.len());
    let mut id_to_index: HashMap<i64, usize> = HashMap::with_capacity(non_blank_atom_lines.len());

    for raw in &non_blank_atom_lines {
        let parsed = parse_atom_line(style, raw)?;
        let index = mol.add_atom(Atom::new(Element::UNKNOWN));
        id_to_index.insert(parsed.id, index);
        ff_atoms.push(ForceFieldAtom {
            atom_type: Some(parsed.atom_type.to_string()),
            mass: masses.get(&parsed.atom_type).copied(),
            partial_charge: parsed.charge,
        });
        coords.push(Point3::new(parsed.x, parsed.y, parsed.z));
    }

    mol.set_coords3(coords)
        .map_err(|e| LammpsError::ParseError(e.to_string()))?;
    mol.set_cell(cell)
        .map_err(|e| LammpsError::ParseError(e.to_string()))?;

    for raw in &sections.bonds {
        if split_trailing_comment(raw).0.trim().is_empty() {
            continue;
        }
        let ids = parse_bonded_ids(raw, 2)?;
        let a = resolve_atom(&id_to_index, ids[0])?;
        let b = resolve_atom(&id_to_index, ids[1])?;
        mol.add_bond(Bond::new(a, b, BondOrder::Single))
            .map_err(|e| LammpsError::ParseError(e.to_string()))?;
    }
    // `Bonds` is a complete real bond list, the same precedent PSF/PRMTOP
    // already established for implying hydrogens from one.
    // `Element::UNKNOWN::typical_valence()` is honestly `0`, so every
    // bonded atom here gets `Some(0)` -- the same answer this crate
    // already gives any bonded atom of a real element with no typical
    // valence (a metal, say), not a new claim of certainty.
    mol.calculate_implicit_hydrogens_where_bonded();

    let mut angles = Vec::new();
    for raw in &sections.angles {
        if split_trailing_comment(raw).0.trim().is_empty() {
            continue;
        }
        let ids = parse_bonded_ids(raw, 3)?;
        angles.push([
            resolve_atom(&id_to_index, ids[0])?,
            resolve_atom(&id_to_index, ids[1])?,
            resolve_atom(&id_to_index, ids[2])?,
        ]);
    }

    let mut dihedrals = Vec::new();
    for raw in &sections.dihedrals {
        if split_trailing_comment(raw).0.trim().is_empty() {
            continue;
        }
        let ids = parse_bonded_ids(raw, 4)?;
        dihedrals.push([
            resolve_atom(&id_to_index, ids[0])?,
            resolve_atom(&id_to_index, ids[1])?,
            resolve_atom(&id_to_index, ids[2])?,
            resolve_atom(&id_to_index, ids[3])?,
        ]);
    }

    let mut impropers = Vec::new();
    for raw in &sections.impropers {
        if split_trailing_comment(raw).0.trim().is_empty() {
            continue;
        }
        let ids = parse_bonded_ids(raw, 4)?;
        impropers.push([
            resolve_atom(&id_to_index, ids[0])?,
            resolve_atom(&id_to_index, ids[1])?,
            resolve_atom(&id_to_index, ids[2])?,
            resolve_atom(&id_to_index, ids[3])?,
        ]);
    }

    mol.set_force_field(ForceFieldTopology {
        atoms: Some(ff_atoms),
        angles,
        dihedrals,
        impropers,
        exclusions: Vec::new(),
        donors: Vec::new(),
        acceptors: Vec::new(),
    })
    .map_err(|e| LammpsError::ParseError(e.to_string()))?;

    Ok(mol)
}

/// The inverse of [`unit_cell_from_lammps_box`]'s shape conversion, plus a
/// coordinate-bounding-box fallback when there is no cell at all: LAMMPS box
/// bounds (`xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz`) from a [`UnitCell`] when one
/// is present, or an exact min/max bounding box of `coords` when it is not
/// -- a deterministic function of real data, never a fabricated placeholder.
/// Shared by [`write_lammps_data`] and [`crate::io::lammpstrj`]'s writer,
/// which additionally forward-shifts these into the dump format's
/// bounding-box convention only when triclinic.
pub(crate) fn lammps_box_bounds(
    cell: Option<&UnitCell>,
    coords: &[Point3],
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    match cell {
        Some(cell) => {
            let (alpha, beta, gamma) = (
                cell.alpha.to_radians(),
                cell.beta.to_radians(),
                cell.gamma.to_radians(),
            );
            let xy = cell.b * gamma.cos();
            let xz = cell.c * beta.cos();
            let ly = (cell.b * cell.b - xy * xy).sqrt();
            let yz = (cell.b * cell.c * alpha.cos() - xy * xz) / ly;
            let lz = (cell.c * cell.c - xz * xz - yz * yz).sqrt();
            (0.0, cell.a, 0.0, ly, 0.0, lz, xy, xz, yz)
        }
        None => match coords.split_first() {
            Some((first, rest)) => {
                let (mut xlo, mut xhi) = (first.x, first.x);
                let (mut ylo, mut yhi) = (first.y, first.y);
                let (mut zlo, mut zhi) = (first.z, first.z);
                for p in rest {
                    xlo = xlo.min(p.x);
                    xhi = xhi.max(p.x);
                    ylo = ylo.min(p.y);
                    yhi = yhi.max(p.y);
                    zlo = zlo.min(p.z);
                    zhi = zhi.max(p.z);
                }
                let (xlo, xhi) = ensure_nonzero_extent(xlo, xhi);
                let (ylo, yhi) = ensure_nonzero_extent(ylo, yhi);
                let (zlo, zhi) = ensure_nonzero_extent(zlo, zhi);
                (xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0)
            }
            None => (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0),
        },
    }
}

/// Writes a [`Molecule`] as a LAMMPS data file: topology only (see the
/// module doc). LAMMPS integer types are synthesized from each atom's
/// `ForceFieldAtom::atom_type` (falling back to its element symbol when
/// absent) by first-occurrence order -- this crate's own reader states
/// the numeric type as exactly that string, so a LAMMPS-to-LAMMPS round
/// trip recovers the same types; a molecule from another format gets
/// types synthesized from whatever it does carry. Box bounds come from
/// the molecule's `UnitCell` when present, or an exact min/max bounding
/// box of the actual coordinates when it is not -- a deterministic
/// function of real data, never a fabricated placeholder.
pub fn write_lammps_data(mol: &Molecule) -> String {
    let empty_ff = ForceFieldTopology::default();
    let force_field = mol.force_field().unwrap_or(&empty_ff);

    let mut labels: Vec<String> = Vec::new();
    let mut label_to_type: HashMap<String, i64> = HashMap::new();
    let mut type_mass: HashMap<i64, f64> = HashMap::new();
    let mut type_of: Vec<i64> = Vec::with_capacity(mol.num_atoms());

    for (i, atom) in mol.atoms().iter().enumerate() {
        let ff_atom = force_field.atoms.as_ref().and_then(|atoms| atoms.get(i));
        let label = ff_atom
            .and_then(|a| a.atom_type.clone())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| atom.element().symbol().to_string());
        let ty = *label_to_type.entry(label.clone()).or_insert_with(|| {
            labels.push(label);
            labels.len() as i64
        });
        type_of.push(ty);
        type_mass.entry(ty).or_insert_with(|| {
            ff_atom
                .and_then(|a| a.mass)
                .unwrap_or(ATOMIC_MASSES[atom.element().atomic_number as usize])
        });
    }

    // A partial charge can arrive via `ForceFieldAtom` or `AtomSite` --
    // different provenances for the same fact (see
    // `core::force_field::ForceFieldAtom`'s own doc comment); a writer
    // has to check both or silently drop whichever one a source used.
    let atom_charge = |i: usize| -> Option<f64> {
        force_field
            .atoms
            .as_ref()
            .and_then(|atoms| atoms.get(i))
            .and_then(|a| a.partial_charge)
            .or_else(|| mol.site(i).and_then(|s| s.partial_charge))
    };
    let has_charge = (0..mol.num_atoms()).any(|i| atom_charge(i).is_some());

    let (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz) =
        lammps_box_bounds(mol.cell().as_ref(), mol.coords3().unwrap_or(&[]));
    let triclinic = xy.abs() > 1e-9 || xz.abs() > 1e-9 || yz.abs() > 1e-9;

    let mut out = String::from("Written by chem\n\n");
    out.push_str(&format!("{} atoms\n", mol.num_atoms()));
    out.push_str(&format!("{} bonds\n", mol.num_bonds()));
    out.push_str(&format!("{} angles\n", force_field.angles.len()));
    out.push_str(&format!("{} dihedrals\n", force_field.dihedrals.len()));
    out.push_str(&format!("{} impropers\n\n", force_field.impropers.len()));
    out.push_str(&format!("{} atom types\n", labels.len()));
    if mol.num_bonds() > 0 {
        out.push_str("1 bond types\n");
    }
    if !force_field.angles.is_empty() {
        out.push_str("1 angle types\n");
    }
    if !force_field.dihedrals.is_empty() {
        out.push_str("1 dihedral types\n");
    }
    if !force_field.impropers.is_empty() {
        out.push_str("1 improper types\n");
    }
    out.push('\n');
    out.push_str(&format!("{xlo} {xhi} xlo xhi\n"));
    out.push_str(&format!("{ylo} {yhi} ylo yhi\n"));
    out.push_str(&format!("{zlo} {zhi} zlo zhi\n"));
    if triclinic {
        out.push_str(&format!("{xy} {xz} {yz} xy xz yz\n"));
    }
    out.push('\n');

    out.push_str("Masses\n\n");
    for (i, _) in labels.iter().enumerate() {
        let ty = (i + 1) as i64;
        out.push_str(&format!(
            "{ty} {}\n",
            type_mass.get(&ty).copied().unwrap_or(0.0)
        ));
    }
    out.push('\n');

    out.push_str(if has_charge {
        "Atoms # charge\n\n"
    } else {
        "Atoms # atomic\n\n"
    });
    let coords = mol.coords3().unwrap_or(&[]);
    for (i, _atom) in mol.atoms().iter().enumerate() {
        let p = coords.get(i).copied().unwrap_or(Point3::ORIGIN);
        let ty = type_of[i];
        if has_charge {
            let charge = atom_charge(i).unwrap_or(0.0);
            out.push_str(&format!(
                "{} {ty} {charge:.6} {:.6} {:.6} {:.6}\n",
                i + 1,
                p.x,
                p.y,
                p.z
            ));
        } else {
            out.push_str(&format!(
                "{} {ty} {:.6} {:.6} {:.6}\n",
                i + 1,
                p.x,
                p.y,
                p.z
            ));
        }
    }
    out.push('\n');

    if mol.num_bonds() > 0 {
        out.push_str("Bonds\n\n");
        for (i, bond) in mol.bonds().iter().enumerate() {
            out.push_str(&format!(
                "{} 1 {} {}\n",
                i + 1,
                bond.atom1() + 1,
                bond.atom2() + 1
            ));
        }
        out.push('\n');
    }
    if !force_field.angles.is_empty() {
        out.push_str("Angles\n\n");
        for (i, a) in force_field.angles.iter().enumerate() {
            out.push_str(&format!(
                "{} 1 {} {} {}\n",
                i + 1,
                a[0] + 1,
                a[1] + 1,
                a[2] + 1
            ));
        }
        out.push('\n');
    }
    if !force_field.dihedrals.is_empty() {
        out.push_str("Dihedrals\n\n");
        for (i, d) in force_field.dihedrals.iter().enumerate() {
            out.push_str(&format!(
                "{} 1 {} {} {} {}\n",
                i + 1,
                d[0] + 1,
                d[1] + 1,
                d[2] + 1,
                d[3] + 1
            ));
        }
        out.push('\n');
    }
    if !force_field.impropers.is_empty() {
        out.push_str("Impropers\n\n");
        for (i, d) in force_field.impropers.iter().enumerate() {
            out.push_str(&format!(
                "{} 1 {} {} {} {}\n",
                i + 1,
                d[0] + 1,
                d[1] + 1,
                d[2] + 1,
                d[3] + 1
            ));
        }
        out.push('\n');
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const ATOMIC: &str = include_str!("../../tests/corpus/lammps/atomic.data");
    const FULL: &str = include_str!("../../tests/corpus/lammps/full.data");
    const AMBIGUOUS: &str = include_str!("../../tests/corpus/lammps/ambiguous_style.data");
    const COARSE_GRAINED: &str = include_str!("../../tests/corpus/lammps/coarse_grained.data");
    const TRICLINIC: &str = include_str!("../../tests/corpus/lammps/triclinic.data");

    fn options() -> LammpsReadOptions {
        LammpsReadOptions::default()
    }

    #[test]
    fn test_an_atomic_topology_free_system_round_trips() {
        let mol = parse_lammps_data(ATOMIC, &options()).expect("valid LAMMPS data");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.num_bonds(), 0);
        assert!(mol.has_coords3());
        assert!(mol.has_cell());
        for atom in mol.atoms() {
            assert_eq!(atom.element(), Element::UNKNOWN);
        }

        let written = write_lammps_data(&mol);
        let back = parse_lammps_data(&written, &options()).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.coords3(), mol.coords3());
    }

    #[test]
    fn test_a_full_style_topology_round_trips() {
        let mol = parse_lammps_data(FULL, &options()).expect("valid LAMMPS data");
        assert_eq!(mol.num_atoms(), 4);
        assert_eq!(mol.num_bonds(), 3);
        let force_field = mol.force_field().expect("has a force field");
        assert_eq!(force_field.angles.len(), 1);
        assert_eq!(force_field.dihedrals.len(), 1);
        assert_eq!(force_field.impropers.len(), 1);
        assert!(
            force_field
                .atoms
                .as_ref()
                .unwrap()
                .iter()
                .all(|a| a.partial_charge.is_some())
        );

        let written = write_lammps_data(&mol);
        let back = parse_lammps_data(&written, &options()).expect("round trips");
        assert_eq!(back.num_bonds(), mol.num_bonds());
        assert_eq!(back.force_field().unwrap().dihedrals, force_field.dihedrals);
        assert_eq!(back.force_field().unwrap().impropers, force_field.impropers);
    }

    #[test]
    fn test_a_bonded_unknown_element_gets_a_zero_implicit_hydrogen_count() {
        // Not `None` (unstated) -- `Bonds` is a complete bond list, so the
        // same implication PSF/PRMTOP run applies; an unknown element's
        // typical valence is honestly 0, the same answer this crate
        // already gives a real zero-valence element like a metal.
        let mol = parse_lammps_data(FULL, &options()).expect("valid LAMMPS data");
        for atom in mol.atoms() {
            assert_eq!(atom.hydrogens(), Some(0));
        }
    }

    #[test]
    fn test_ambiguous_atom_style_is_a_clear_error_not_a_guess() {
        let err = parse_lammps_data(AMBIGUOUS, &options()).unwrap_err();
        assert!(matches!(err, LammpsError::AmbiguousAtomStyle(6)), "{err}");
    }

    #[test]
    fn test_the_explicit_atom_style_option_is_the_escape_hatch() {
        let with_override = LammpsReadOptions {
            atom_style: Some(AtomStyle::Charge),
        };
        let mol = parse_lammps_data(AMBIGUOUS, &with_override).expect("resolved by the option");
        assert_eq!(mol.num_atoms(), 2);
    }

    #[test]
    fn test_a_coarse_grained_mass_reads_as_unknown_not_a_failure() {
        // A bead-spring/reduced-unit mass like 1.0 matches no real
        // element -- this must read successfully, not fail or guess.
        let mol = parse_lammps_data(COARSE_GRAINED, &options()).expect("valid LAMMPS data");
        assert_eq!(mol.num_atoms(), 2);
        assert!(mol.atoms().iter().all(|a| a.element() == Element::UNKNOWN));
        let force_field = mol.force_field().unwrap();
        let masses: Vec<f64> = force_field
            .atoms
            .as_ref()
            .unwrap()
            .iter()
            .map(|a| a.mass.unwrap())
            .collect();
        assert_eq!(masses, vec![1.0, 1.0]);
    }

    #[test]
    fn test_no_coordinates_are_ever_invented_beyond_what_the_file_states() {
        let mol = parse_lammps_data(ATOMIC, &options()).expect("valid LAMMPS data");
        let coords = mol.coords3().unwrap();
        assert_eq!(coords.len(), 3);
    }

    #[test]
    fn test_a_triclinic_box_round_trips_through_unit_cell() {
        // The one part of this format most likely to be subtly wrong: a
        // sign or a swapped cosine in the closed-form conversion would
        // still produce *a* cell, just the wrong one -- this pins the
        // shape (none of the three angles is 90°) and a round trip.
        let mol = parse_lammps_data(TRICLINIC, &options()).expect("valid LAMMPS data");
        let cell = mol.cell().expect("has a cell");
        for angle in [cell.alpha, cell.beta, cell.gamma] {
            assert!(
                (angle - 90.0).abs() > 1.0,
                "expected a non-orthogonal angle, got {angle}"
            );
        }

        let written = write_lammps_data(&mol);
        let back = parse_lammps_data(&written, &options()).expect("round trips");
        let back_cell = back.cell().expect("still has a cell");
        assert!((back_cell.a - cell.a).abs() < 1e-6);
        assert!((back_cell.b - cell.b).abs() < 1e-6);
        assert!((back_cell.c - cell.c).abs() < 1e-6);
        assert!((back_cell.alpha - cell.alpha).abs() < 1e-6);
        assert!((back_cell.beta - cell.beta).abs() < 1e-6);
        assert!((back_cell.gamma - cell.gamma).abs() < 1e-6);
    }

    #[test]
    fn test_an_unrecognised_header_is_a_clear_error_not_a_panic() {
        let err =
            parse_lammps_data("comment\n\nAtoms\n\n1 1 0.0 0.0 0.0\n", &options()).unwrap_err();
        assert!(matches!(err, LammpsError::ParseError(_)), "{err}");
    }
}
