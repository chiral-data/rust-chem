//! PRMTOP — AMBER's topology file, also named `.parm7`: a flat sequence of
//! `%FLAG <NAME>` sections, each preceded by a `%FORMAT(<fortran-format>)`
//! line stating exactly how to parse the data that follows (#322).
//!
//! **Topology only.** PRMTOP states substantially more of the force field
//! than PSF does -- not just topology but the actual force-field
//! *parameters* (force constants, equilibrium values, Lennard-Jones
//! coefficients). Those have no home in
//! [`crate::core::force_field::ForceFieldTopology`], which already frames
//! them as "out of scope for this milestone", and `Carries` is completely
//! full at 32/32 bits -- widening it is a separate change this story does
//! not take on. Every parameter section here is read far enough to skip
//! past correctly and never stored, the same treatment PSF gives `!NGRP`.
//!
//! **The `%FORMAT` line is the parser, not a hint.** Every section states
//! its own Fortran format (`nAw` strings, `nIw` integers, `nEw.d` floats) --
//! `n` values per line, each exactly `w` columns wide, with no guaranteed
//! whitespace between adjacent values (two 4-character names can butt up
//! against each other). This reader slices by the declared width rather
//! than whitespace-tokenizing, the opposite of PSF's approach, and reads a
//! section's data until the next `%FLAG` line rather than trusting a
//! separately-stated count -- the boundary *is* the count.
//!
//! **Charges are pre-scaled.** AMBER stores partial charge × this
//! module's own `AMBER_CHARGE_SCALE` (18.2223, the literal constant real
//! tools use -- not a fuller-precision recomputation of `sqrt(k_e)`).
//! Converted at the boundary in both directions.
//!
//! **`EXCLUDED_ATOMS_LIST`/`NUMBER_EXCLUDED_ATOMS` is a different shape
//! than PSF's `!NNB`.** Not a cumulative pointer array -- a per-atom
//! *count* array consumed sequentially against a flat partner list, where
//! an atom with no real exclusions is still given one placeholder entry:
//! count `1` paired with a literal `0` (a non-existent atom) in the list.
//! See this module's own `decode_exclusions`/`encode_exclusions`.
//!
//! **Bond/angle/dihedral atom fields are `3×(0-based atom index)`, not a
//! serial** -- an AMBER-internal optimization (these were originally
//! coordinate-array offsets). Dihedral tuples carry two sign tricks: a
//! negative 3rd atom index means "skip this torsion's 1-4 nonbonded
//! interaction" (no data-model slot -- discarded), and a negative 4th atom
//! index means "this is actually an improper", which *does* have a home:
//! `ForceFieldTopology` already separates `dihedrals`/`impropers`. See
//! this module's own `decode_dihedrals`.
//!
//! **No chain concept exists in PRMTOP** -- only `RESIDUE_LABEL`/
//! `RESIDUE_POINTER`. Every atom gets the same placeholder chain id when
//! grouped via `cif_model`'s own `group_into_chains_and_residues`, the same
//! algorithm PSF reuses, honestly reflecting that the format states no
//! more than one implicit chain.
//!
//! Bonds go through [`Molecule`]'s ordinary bond list, not
//! `ForceFieldTopology` -- the same precedent PSF and PDB's `CONECT`
//! already established. `BONDS_INC_HYDROGEN`/`BONDS_WITHOUT_HYDROGEN`
//! (and the angle/dihedral equivalents) are two arrays forming one
//! complete list -- an AMBER-internal split, not a semantic one --
//! concatenated on read and re-split by actual hydrogen content on write.

use std::collections::HashMap;

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::elements::ATOMIC_MASSES;
use crate::core::force_field::{ForceFieldAtom, ForceFieldTopology};
use crate::core::molecule::Molecule;
use crate::core::site::AtomSite;
use crate::io::cif_model::{ResidueKey, group_into_chains_and_residues};
use crate::io::errors::PrmtopError;

/// Amber's own literal charge-scaling constant (ParmEd's own
/// `AMBER_ELECTROSTATIC`) -- not a fuller-precision recomputation of
/// `sqrt(k_e)`, so this crate's round trip matches what real AMBER tools
/// use bit-for-bit.
const AMBER_CHARGE_SCALE: f64 = 18.2223;

/// PRMTOP states no chain/segment concept at all -- every atom gets this
/// same placeholder when grouped into chains/residues.
const PLACEHOLDER_CHAIN_ID: &str = "";

/// Parses a `%FORMAT(...)` line into the field width it declares -- the
/// only part of the Fortran format spec this reader needs. Section
/// boundaries (the next `%FLAG` line) already say how many fields there
/// are; interpreting a token as an integer, float or string happens
/// afterward, keyed by the section's own name.
fn parse_fortran_format(line: &str) -> Result<usize, PrmtopError> {
    let inner = line
        .trim()
        .strip_prefix("%FORMAT(")
        .and_then(|s| s.strip_suffix(')'))
        .ok_or_else(|| PrmtopError::ParseError(format!("malformed %FORMAT line: {line:?}")))?;
    let digit_end = inner
        .find(|c: char| !c.is_ascii_digit())
        .filter(|&n| n > 0)
        .ok_or_else(|| {
            PrmtopError::ParseError(format!("missing repeat count in %FORMAT: {line:?}"))
        })?;
    let mut chars = inner[digit_end..].chars();
    chars.next().ok_or_else(|| {
        PrmtopError::ParseError(format!("missing format kind in %FORMAT: {line:?}"))
    })?;
    let width_and_decimals = chars.as_str();
    let width_str = width_and_decimals
        .split('.')
        .next()
        .unwrap_or(width_and_decimals);
    width_str
        .parse()
        .map_err(|_| PrmtopError::ParseError(format!("bad field width in %FORMAT: {line:?}")))
}

/// Slices `line` into `width`-character chunks, trimmed. PRMTOP's fields
/// have no guaranteed separator between them -- this is the whole reason
/// the reader cannot whitespace-tokenize the way PSF does.
fn split_fixed_width(line: &str, width: usize) -> Vec<String> {
    if width == 0 {
        return Vec::new();
    }
    let chars: Vec<char> = line.chars().collect();
    let mut out = Vec::new();
    let mut pos = 0;
    while pos < chars.len() {
        let end = (pos + width).min(chars.len());
        out.push(
            chars[pos..end]
                .iter()
                .collect::<String>()
                .trim()
                .to_string(),
        );
        pos = end;
    }
    out
}

/// Scans a PRMTOP file into its `%FLAG` sections, keyed by name, as their
/// raw fixed-width tokens (already sliced by the section's own `%FORMAT`
/// and trimmed). A section's data runs until the next `%FLAG` line or
/// end of file -- there is no separate stated count to trust or distrust.
fn scan_sections(text: &str) -> Result<HashMap<String, Vec<String>>, PrmtopError> {
    let mut lines = text.lines().peekable();

    let header = lines
        .next()
        .ok_or_else(|| PrmtopError::ParseError("empty file".to_string()))?;
    if !header.trim_start().starts_with("%VERSION") {
        return Err(PrmtopError::ParseError(format!(
            "not a PRMTOP file: {header:?}"
        )));
    }

    let mut sections = HashMap::new();
    while let Some(line) = lines.next() {
        let trimmed = line.trim_start();
        if !trimmed.starts_with("%FLAG") {
            continue;
        }
        let name = trimmed.trim_start_matches("%FLAG").trim().to_string();

        while lines
            .peek()
            .is_some_and(|l| l.trim_start().starts_with("%COMMENT"))
        {
            lines.next();
        }

        let format_line = lines
            .next()
            .ok_or_else(|| PrmtopError::ParseError(format!("{name} has no %FORMAT line")))?;
        if !format_line.trim_start().starts_with("%FORMAT") {
            return Err(PrmtopError::ParseError(format!(
                "expected %FORMAT after %FLAG {name}, found: {format_line:?}"
            )));
        }
        let width = parse_fortran_format(format_line)?;

        let mut tokens = Vec::new();
        while let Some(&data_line) = lines.peek() {
            if data_line.trim_start().starts_with('%') {
                break;
            }
            lines.next();
            tokens.extend(split_fixed_width(data_line, width));
        }
        sections.insert(name, tokens);
    }
    Ok(sections)
}

fn ints(sections: &HashMap<String, Vec<String>>, name: &str) -> Result<Vec<i64>, PrmtopError> {
    match sections.get(name) {
        Some(tokens) => tokens
            .iter()
            .map(|t| {
                t.parse().map_err(|_| {
                    PrmtopError::ParseError(format!("{name}: expected an integer, got {t:?}"))
                })
            })
            .collect(),
        None => Ok(Vec::new()),
    }
}

fn floats(sections: &HashMap<String, Vec<String>>, name: &str) -> Result<Vec<f64>, PrmtopError> {
    match sections.get(name) {
        Some(tokens) => tokens
            .iter()
            .map(|t| {
                t.parse().map_err(|_| {
                    PrmtopError::ParseError(format!("{name}: expected a number, got {t:?}"))
                })
            })
            .collect(),
        None => Ok(Vec::new()),
    }
}

fn strings<'a>(sections: &'a HashMap<String, Vec<String>>, name: &str) -> &'a [String] {
    sections.get(name).map(Vec::as_slice).unwrap_or(&[])
}

/// Infers an element from `ATOMIC_NUMBER` when present and in range,
/// falling back to PSF's own mass-based inference otherwise -- older
/// PRMTOP files predate the `ATOMIC_NUMBER` section (AmberTools 12+), and
/// every real reader (ParmEd included) falls back to mass the same way.
fn element_from_atomic_number_or_mass(atomic_number: Option<i64>, mass: f64) -> Option<Element> {
    atomic_number
        .and_then(|n| u8::try_from(n).ok())
        .and_then(Element::new)
        .or_else(|| crate::io::psf::element_from_mass(mass))
}

/// Converts PRMTOP's 1-based atom serial (as used by `RESIDUE_POINTER` and
/// `EXCLUDED_ATOMS_LIST`) into a 0-based index, checked against
/// `num_atoms` rather than trusted.
fn atom_index(serial: i64, num_atoms: usize) -> Result<usize, PrmtopError> {
    if serial < 1 || serial as usize > num_atoms {
        return Err(PrmtopError::ParseError(format!(
            "atom serial {serial} is out of range for {num_atoms} atoms"
        )));
    }
    Ok(serial as usize - 1)
}

/// Converts a bonded-term's stored coordinate-array offset (`3 ×` a
/// 0-based atom index) into that index, checked against `num_atoms`.
/// AMBER's own optimization: these integers were originally offsets into
/// a flat x/y/z coordinate array, not atom serials. The sign a caller may
/// have already stripped (dihedral tuples use it for the improper/skip-1-4
/// tricks) -- this only ever sees the magnitude.
fn coord_index(raw: i64, num_atoms: usize) -> Result<usize, PrmtopError> {
    let magnitude = raw.unsigned_abs();
    if !magnitude.is_multiple_of(3) {
        return Err(PrmtopError::ParseError(format!(
            "bonded-term index {raw} is not a multiple of 3"
        )));
    }
    let idx = (magnitude / 3) as usize;
    if idx >= num_atoms {
        return Err(PrmtopError::ParseError(format!(
            "bonded-term index {raw} is out of range for {num_atoms} atoms"
        )));
    }
    Ok(idx)
}

/// Reads a flat, concatenated bonds/angles section (each entry is `N`
/// coordinate-offset atom fields followed by one type-index field this
/// story does not store) into `N`-atom terms.
fn decode_bonded_terms<const N: usize>(
    raw: &[i64],
    num_atoms: usize,
) -> Result<Vec<[usize; N]>, PrmtopError> {
    if !raw.len().is_multiple_of(N + 1) {
        return Err(PrmtopError::ParseError(format!(
            "a bonded-term section's length is not a multiple of {}",
            N + 1
        )));
    }
    let mut terms = Vec::with_capacity(raw.len() / (N + 1));
    for chunk in raw.chunks_exact(N + 1) {
        let mut atoms = [0usize; N];
        for (slot, &value) in atoms.iter_mut().zip(chunk) {
            *slot = coord_index(value, num_atoms)?;
        }
        terms.push(atoms);
    }
    Ok(terms)
}

/// Proper torsions, then impropers -- both `[usize; 4]` term lists, split
/// by [`decode_dihedrals`] on the stored 4th field's sign.
type Dihedrals = (Vec<[usize; 4]>, Vec<[usize; 4]>);

/// Reads a flat, concatenated dihedrals section into proper/improper
/// 4-atom terms, splitting on the 4th field's sign. The 3rd field's own
/// sign ("skip this torsion's 1-4 nonbonded interaction") has no
/// data-model slot and is discarded once read.
fn decode_dihedrals(raw: &[i64], num_atoms: usize) -> Result<Dihedrals, PrmtopError> {
    if !raw.len().is_multiple_of(5) {
        return Err(PrmtopError::ParseError(
            "a dihedral section's length is not a multiple of 5".to_string(),
        ));
    }
    let mut dihedrals = Vec::new();
    let mut impropers = Vec::new();
    for chunk in raw.chunks_exact(5) {
        let is_improper = chunk[3] < 0;
        let atoms = [
            coord_index(chunk[0], num_atoms)?,
            coord_index(chunk[1], num_atoms)?,
            coord_index(chunk[2], num_atoms)?,
            coord_index(chunk[3], num_atoms)?,
        ];
        if is_improper {
            impropers.push(atoms);
        } else {
            dihedrals.push(atoms);
        }
    }
    Ok((dihedrals, impropers))
}

/// Reconstructs `ForceFieldTopology::exclusions` from PRMTOP's own shape:
/// `NUMBER_EXCLUDED_ATOMS[i]` says how many of the next entries in
/// `EXCLUDED_ATOMS_LIST` belong to atom `i`, consumed sequentially (not a
/// cumulative pointer array like PSF's `!NNB`). An atom with no real
/// exclusions is still given one placeholder entry: count `1` paired with
/// a literal `0` (a non-existent atom), which this skips rather than
/// treating as a real exclusion against atom zero.
fn decode_exclusions(
    counts: &[i64],
    partners: &[i64],
    num_atoms: usize,
) -> Result<Vec<[usize; 2]>, PrmtopError> {
    if counts.len() != num_atoms {
        return Err(PrmtopError::ParseError(format!(
            "NUMBER_EXCLUDED_ATOMS has {} entries for {num_atoms} atoms",
            counts.len()
        )));
    }
    let mut exclusions = Vec::new();
    let mut pos = 0usize;
    for (atom_ix, &count) in counts.iter().enumerate() {
        let count = usize::try_from(count)
            .map_err(|_| PrmtopError::ParseError(format!("negative exclusion count {count}")))?;
        let end = pos + count;
        if end > partners.len() {
            return Err(PrmtopError::ParseError(
                "EXCLUDED_ATOMS_LIST is shorter than NUMBER_EXCLUDED_ATOMS promises".to_string(),
            ));
        }
        for &partner in &partners[pos..end] {
            if partner != 0 {
                exclusions.push([atom_ix, atom_index(partner, num_atoms)?]);
            }
        }
        pos = end;
    }
    Ok(exclusions)
}

/// The inverse of [`decode_exclusions`]: builds `NUMBER_EXCLUDED_ATOMS`/
/// `EXCLUDED_ATOMS_LIST`, keyed by the smaller of each pair (the format's
/// own convention -- an exclusion is stated only once, against the atom
/// with the smaller index), using the `1`/`0` sentinel for an atom with no
/// real exclusions rather than a literal `0` count.
fn encode_exclusions(exclusions: &[[usize; 2]], num_atoms: usize) -> (Vec<i64>, Vec<i64>) {
    let mut partners: Vec<Vec<usize>> = vec![Vec::new(); num_atoms];
    for &[a, b] in exclusions {
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        if lo < num_atoms {
            partners[lo].push(hi);
        }
    }
    for list in &mut partners {
        list.sort_unstable();
    }

    let mut counts = Vec::with_capacity(num_atoms);
    let mut flat = Vec::new();
    for list in &partners {
        if list.is_empty() {
            counts.push(1);
            flat.push(0);
        } else {
            counts.push(list.len() as i64);
            flat.extend(list.iter().map(|&p| p as i64 + 1));
        }
    }
    (counts, flat)
}

/// Parses a whole PRMTOP file into a [`Molecule`].
pub fn parse_prmtop(text: &str) -> Result<Molecule, PrmtopError> {
    let sections = scan_sections(text)?;

    let atom_names = sections.get("ATOM_NAME").ok_or(PrmtopError::NoAtoms)?;
    let natom = atom_names.len();
    if natom == 0 {
        return Err(PrmtopError::NoAtoms);
    }

    let charges = floats(&sections, "CHARGE")?;
    let masses = floats(&sections, "MASS")?;
    let atomic_numbers = ints(&sections, "ATOMIC_NUMBER")?;
    let amber_types = strings(&sections, "AMBER_ATOM_TYPE");
    let residue_labels = strings(&sections, "RESIDUE_LABEL");
    let residue_pointers = ints(&sections, "RESIDUE_POINTER")?;

    if charges.len() != natom || masses.len() != natom {
        return Err(PrmtopError::ParseError(format!(
            "ATOM_NAME/CHARGE/MASS lengths disagree: {natom}/{}/{}",
            charges.len(),
            masses.len()
        )));
    }

    // Residue boundaries: RESIDUE_POINTER gives each residue's first atom
    // (1-based); the next residue's pointer (or NATOM for the last) gives
    // the span's end.
    let mut residue_of_atom = vec![0usize; natom];
    for (residue_ix, &start) in residue_pointers.iter().enumerate() {
        let start = atom_index(start, natom)?;
        let end = residue_pointers
            .get(residue_ix + 1)
            .map(|&next| atom_index(next, natom))
            .transpose()?
            .unwrap_or(natom);
        if end < start {
            return Err(PrmtopError::ParseError(format!(
                "RESIDUE_POINTER is not ascending at residue {residue_ix}"
            )));
        }
        for slot in residue_of_atom.iter_mut().skip(start).take(end - start) {
            *slot = residue_ix;
        }
    }

    let mut mol = Molecule::new();
    let mut sites = Vec::with_capacity(natom);
    let mut ff_atoms = Vec::with_capacity(natom);

    for i in 0..natom {
        let mass = masses[i];
        let atomic_number = atomic_numbers.get(i).copied();
        let element = element_from_atomic_number_or_mass(atomic_number, mass).ok_or(
            PrmtopError::InvalidElement {
                atomic_number,
                mass,
            },
        )?;
        mol.add_atom(Atom::new(element));

        sites.push(AtomSite {
            name: Some(atom_names[i].clone()),
            ..AtomSite::empty()
        });
        ff_atoms.push(ForceFieldAtom {
            atom_type: amber_types.get(i).cloned(),
            mass: Some(mass),
            partial_charge: Some(charges[i] / AMBER_CHARGE_SCALE),
        });
    }

    mol.set_sites(sites)
        .map_err(|e| PrmtopError::ParseError(e.to_string()))?;

    let keys: Vec<ResidueKey> = (0..natom)
        .map(|i| {
            let residue_ix = residue_of_atom[i];
            ResidueKey {
                chain_id: PLACEHOLDER_CHAIN_ID.to_string(),
                name: residue_labels
                    .get(residue_ix)
                    .cloned()
                    .unwrap_or_else(|| "UNK".to_string()),
                sequence: residue_ix as i32 + 1,
                insertion_code: None,
                is_hetero: false,
            }
        })
        .collect();
    let (chains, residues) = group_into_chains_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| PrmtopError::ParseError(e.to_string()))?;

    let mut bond_raw = ints(&sections, "BONDS_INC_HYDROGEN")?;
    bond_raw.extend(ints(&sections, "BONDS_WITHOUT_HYDROGEN")?);
    for [a, b] in decode_bonded_terms::<2>(&bond_raw, natom)? {
        mol.add_bond(Bond::new(a, b, BondOrder::Single))
            .map_err(|e| PrmtopError::ParseError(e.to_string()))?;
    }
    // BONDS_INC_HYDROGEN + BONDS_WITHOUT_HYDROGEN together are PRMTOP's
    // complete bond list -- an AMBER-internal optimisation split, not a
    // semantic one -- so implying hydrogens from it is legitimate, the
    // same precedent PSF's `!NBOND` reader already established (#321).
    mol.calculate_implicit_hydrogens_where_bonded();

    let mut angle_raw = ints(&sections, "ANGLES_INC_HYDROGEN")?;
    angle_raw.extend(ints(&sections, "ANGLES_WITHOUT_HYDROGEN")?);
    let angles = decode_bonded_terms::<3>(&angle_raw, natom)?;

    let mut dihedral_raw = ints(&sections, "DIHEDRALS_INC_HYDROGEN")?;
    dihedral_raw.extend(ints(&sections, "DIHEDRALS_WITHOUT_HYDROGEN")?);
    let (dihedrals, impropers) = decode_dihedrals(&dihedral_raw, natom)?;

    let counts = ints(&sections, "NUMBER_EXCLUDED_ATOMS")?;
    let partners = ints(&sections, "EXCLUDED_ATOMS_LIST")?;
    let exclusions = decode_exclusions(&counts, &partners, natom)?;

    mol.set_force_field(ForceFieldTopology {
        atoms: Some(ff_atoms),
        angles,
        dihedrals,
        impropers,
        exclusions,
        donors: Vec::new(),
        acceptors: Vec::new(),
    })
    .map_err(|e| PrmtopError::ParseError(e.to_string()))?;

    Ok(mol)
}

fn push_int_section(out: &mut String, flag: &str, values: &[i64]) {
    out.push_str(&format!("%FLAG {flag}\n%FORMAT(10I8)\n"));
    for chunk in values.chunks(10) {
        let line: String = chunk.iter().map(|v| format!("{v:>8}")).collect();
        out.push_str(&line);
        out.push('\n');
    }
}

fn push_float_section(out: &mut String, flag: &str, values: &[f64]) {
    out.push_str(&format!("%FLAG {flag}\n%FORMAT(5E16.8)\n"));
    for chunk in values.chunks(5) {
        let line: String = chunk
            .iter()
            .map(|&v| format!("{:>16}", format_amber_float(v)))
            .collect();
        out.push_str(&line);
        out.push('\n');
    }
}

fn push_string_section(out: &mut String, flag: &str, values: &[String], width: usize) {
    out.push_str(&format!("%FLAG {flag}\n%FORMAT(20a{width})\n"));
    for chunk in values.chunks(20) {
        let line: String = chunk
            .iter()
            .map(|v| {
                let truncated: String = v.chars().take(width).collect();
                format!("{truncated:<width$}")
            })
            .collect();
        out.push_str(&line);
        out.push('\n');
    }
}

/// Formats a float the way AMBER's own `%FORMAT(5E16.8)` fields look:
/// `[sign]D.DDDDDDDDE±DD`, one nonzero digit before the decimal point,
/// eight after, a two-digit signed exponent -- not merely a same-width
/// decimal, since a real reader (this one included) trusts the column
/// width, not the notation.
fn format_amber_float(value: f64) -> String {
    if value == 0.0 {
        return "0.00000000E+00".to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let magnitude = value.abs();
    let mut exponent = magnitude.log10().floor() as i32;
    let mut mantissa = magnitude / 10f64.powi(exponent);
    if mantissa >= 10.0 {
        mantissa /= 10.0;
        exponent += 1;
    } else if mantissa < 1.0 {
        mantissa *= 10.0;
        exponent -= 1;
    }
    let exp_sign = if exponent < 0 { "-" } else { "+" };
    format!("{sign}{mantissa:.8}E{exp_sign}{:02}", exponent.abs())
}

/// Swaps a would-be improper's atoms so atom index `0`, if present, never
/// sits in the 4th (sign-bearing) slot -- AMBER's own constraint, since
/// `-0` is indistinguishable from `0` and a reader could not tell the term
/// was ever meant to be an improper. Real files satisfy this by ordering
/// their torsions so the first atom in the topology is always listed
/// first or second; this reorders on the way out instead of requiring the
/// caller to have done so.
fn improper_atoms_for_write(atoms: [usize; 4]) -> [usize; 4] {
    let mut atoms = atoms;
    if atoms[3] == 0 {
        atoms.swap(0, 3);
    }
    atoms
}

/// Writes a [`Molecule`] as a PRMTOP file: topology only (see the module
/// doc) -- no force-field parameter section is ever emitted, so `POINTERS`
/// states zero for every parameter-table size it names.
pub fn write_prmtop(mol: &Molecule) -> String {
    let num_atoms = mol.num_atoms();
    let empty_ff = ForceFieldTopology::default();
    let force_field = mol.force_field().unwrap_or(&empty_ff);

    let mut out = String::from("%VERSION  VERSION_STAMP = V0001.000  DATE = 00/00/00  00:00:00\n");
    out.push_str("%FLAG TITLE\n%FORMAT(20a4)\n");

    let is_hydrogen = |ix: usize| mol.atoms()[ix].element().atomic_number == 1;

    let mut bonds_h = Vec::new();
    let mut bonds_no_h = Vec::new();
    for bond in mol.bonds() {
        let (a, b) = (bond.atom1(), bond.atom2());
        let entry = [a as i64 * 3, b as i64 * 3, 0];
        if is_hydrogen(a) || is_hydrogen(b) {
            bonds_h.push(entry);
        } else {
            bonds_no_h.push(entry);
        }
    }

    let mut angles_h = Vec::new();
    let mut angles_no_h = Vec::new();
    for &[a, b, c] in &force_field.angles {
        let entry = [a as i64 * 3, b as i64 * 3, c as i64 * 3, 0];
        if is_hydrogen(a) || is_hydrogen(b) || is_hydrogen(c) {
            angles_h.push(entry);
        } else {
            angles_no_h.push(entry);
        }
    }

    let mut dihedrals_h = Vec::new();
    let mut dihedrals_no_h = Vec::new();
    for &[a, b, c, d] in &force_field.dihedrals {
        let entry = [a as i64 * 3, b as i64 * 3, c as i64 * 3, d as i64 * 3, 0];
        if [a, b, c, d].iter().any(|&ix| is_hydrogen(ix)) {
            dihedrals_h.push(entry);
        } else {
            dihedrals_no_h.push(entry);
        }
    }
    for &atoms in &force_field.impropers {
        let [a, b, c, d] = improper_atoms_for_write(atoms);
        let entry = [a as i64 * 3, b as i64 * 3, c as i64 * 3, -(d as i64 * 3), 0];
        if [a, b, c, d].iter().any(|&ix| is_hydrogen(ix)) {
            dihedrals_h.push(entry);
        } else {
            dihedrals_no_h.push(entry);
        }
    }

    let (exclusion_counts, exclusion_partners) =
        encode_exclusions(&force_field.exclusions, num_atoms);

    let residues = mol.residues();
    let (residue_pointers, residue_labels): (Vec<i64>, Vec<String>) = if residues.is_empty() {
        (vec![1], vec!["UNK".to_string()])
    } else {
        (
            residues.iter().map(|r| r.atoms.start as i64 + 1).collect(),
            residues.iter().map(|r| r.name.clone()).collect(),
        )
    };
    let max_residue_atoms = if residues.is_empty() {
        num_atoms
    } else {
        residues.iter().map(|r| r.num_atoms()).max().unwrap_or(0)
    };

    let pointers: Vec<i64> = vec![
        num_atoms as i64,                // 1 NATOM
        0,                               // 2 NTYPES
        bonds_h.len() as i64,            // 3 NBONH
        bonds_no_h.len() as i64,         // 4 MBONA
        angles_h.len() as i64,           // 5 NTHETH
        angles_no_h.len() as i64,        // 6 MTHETA
        dihedrals_h.len() as i64,        // 7 NPHIH
        dihedrals_no_h.len() as i64,     // 8 MPHIA
        0,                               // 9 NHPARM
        0,                               // 10 NPARM
        exclusion_partners.len() as i64, // 11 NNB
        residue_pointers.len() as i64,   // 12 NRES
        bonds_no_h.len() as i64,         // 13 NBONA
        angles_no_h.len() as i64,        // 14 NTHETA
        dihedrals_no_h.len() as i64,     // 15 NPHIA
        0,
        0,
        0,
        0,
        0, // 16-20 NUMBND,NUMANG,NPTRA,NATYP,NPHB
        0,
        0,
        0,
        0,
        0,
        0,
        0,                        // 21-27 IFPERT,NBPER,NGPER,NDPER,MBPER,MGPER,MDPER
        0,                        // 28 IFBOX
        max_residue_atoms as i64, // 29 NMXRS
        0,                        // 30 IFCAP
        0,                        // 31 NUMEXTRA
    ];
    push_int_section(&mut out, "POINTERS", &pointers);

    let names: Vec<String> = (0..num_atoms)
        .map(|i| {
            mol.site(i)
                .and_then(|s| s.name.as_deref())
                .filter(|s| !s.is_empty())
                .unwrap_or(mol.atoms()[i].element().symbol())
                .to_string()
        })
        .collect();
    push_string_section(&mut out, "ATOM_NAME", &names, 4);

    let charges: Vec<f64> = (0..num_atoms)
        .map(|i| {
            let ff_atom = force_field.atoms.as_ref().and_then(|atoms| atoms.get(i));
            ff_atom.and_then(|a| a.partial_charge).unwrap_or(0.0) * AMBER_CHARGE_SCALE
        })
        .collect();
    push_float_section(&mut out, "CHARGE", &charges);

    let atomic_numbers: Vec<i64> = mol
        .atoms()
        .iter()
        .map(|atom| atom.element().atomic_number as i64)
        .collect();
    push_int_section(&mut out, "ATOMIC_NUMBER", &atomic_numbers);

    let masses: Vec<f64> = (0..num_atoms)
        .map(|i| {
            let ff_atom = force_field.atoms.as_ref().and_then(|atoms| atoms.get(i));
            ff_atom
                .and_then(|a| a.mass)
                .unwrap_or(ATOMIC_MASSES[mol.atoms()[i].element().atomic_number as usize])
        })
        .collect();
    push_float_section(&mut out, "MASS", &masses);

    let amber_types: Vec<String> = (0..num_atoms)
        .map(|i| {
            let ff_atom = force_field.atoms.as_ref().and_then(|atoms| atoms.get(i));
            ff_atom
                .and_then(|a| a.atom_type.as_deref())
                .filter(|s| !s.is_empty())
                .unwrap_or(mol.atoms()[i].element().symbol())
                .to_string()
        })
        .collect();
    push_string_section(&mut out, "AMBER_ATOM_TYPE", &amber_types, 4);

    push_string_section(&mut out, "RESIDUE_LABEL", &residue_labels, 4);
    push_int_section(&mut out, "RESIDUE_POINTER", &residue_pointers);

    push_int_section(
        &mut out,
        "BONDS_INC_HYDROGEN",
        &bonds_h.into_iter().flatten().collect::<Vec<_>>(),
    );
    push_int_section(
        &mut out,
        "BONDS_WITHOUT_HYDROGEN",
        &bonds_no_h.into_iter().flatten().collect::<Vec<_>>(),
    );
    push_int_section(
        &mut out,
        "ANGLES_INC_HYDROGEN",
        &angles_h.into_iter().flatten().collect::<Vec<_>>(),
    );
    push_int_section(
        &mut out,
        "ANGLES_WITHOUT_HYDROGEN",
        &angles_no_h.into_iter().flatten().collect::<Vec<_>>(),
    );
    push_int_section(
        &mut out,
        "DIHEDRALS_INC_HYDROGEN",
        &dihedrals_h.into_iter().flatten().collect::<Vec<_>>(),
    );
    push_int_section(
        &mut out,
        "DIHEDRALS_WITHOUT_HYDROGEN",
        &dihedrals_no_h.into_iter().flatten().collect::<Vec<_>>(),
    );

    push_int_section(&mut out, "NUMBER_EXCLUDED_ATOMS", &exclusion_counts);
    push_int_section(&mut out, "EXCLUDED_ATOMS_LIST", &exclusion_partners);

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = include_str!("../../tests/corpus/prmtop/water.prmtop");
    const EXCLUSIONS: &str = include_str!("../../tests/corpus/prmtop/exclusions.prmtop");
    const NO_ATOMIC_NUMBER: &str =
        include_str!("../../tests/corpus/prmtop/no_atomic_number.prmtop");

    #[test]
    fn test_parse_fortran_format_reads_repeat_kind_and_width() {
        assert_eq!(parse_fortran_format("%FORMAT(10I8)").unwrap(), 8);
        assert_eq!(parse_fortran_format("%FORMAT(5E16.8)").unwrap(), 16);
        assert_eq!(parse_fortran_format("%FORMAT(20a4)").unwrap(), 4);
    }

    #[test]
    fn test_a_water_topology_round_trips() {
        let mol = parse_prmtop(WATER).expect("valid PRMTOP");
        assert_eq!(mol.num_atoms(), 3);
        assert!(!mol.has_coords3(), "PRMTOP states no coordinates");
        assert_eq!(mol.num_bonds(), 2);
        assert_eq!(mol.chains().len(), 1);
        assert_eq!(mol.residues().len(), 1);
        assert_eq!(mol.residues()[0].name, "WAT");

        let force_field = mol.force_field().expect("has a force field");
        assert_eq!(force_field.angles.len(), 1);
        assert_eq!(
            force_field.atoms.as_ref().unwrap()[0].atom_type.as_deref(),
            Some("OW")
        );

        let written = write_prmtop(&mol);
        let back = parse_prmtop(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.num_bonds(), mol.num_bonds());
        assert_eq!(
            back.force_field().unwrap().angles,
            mol.force_field().unwrap().angles
        );
    }

    #[test]
    fn test_charges_are_unscaled_and_rescaled_through_the_amber_constant() {
        let mol = parse_prmtop(WATER).expect("valid PRMTOP");
        let charge = mol.force_field().unwrap().atoms.as_ref().unwrap()[0]
            .partial_charge
            .unwrap();
        // The fixture states -1.64200000E+00 for the oxygen -- confirm it
        // came back divided by AMBER_CHARGE_SCALE, not left in Amber units.
        assert!((charge - (-1.642 / AMBER_CHARGE_SCALE)).abs() < 1e-6);
    }

    #[test]
    fn test_the_exclusion_sentinel_survives_a_round_trip() {
        // The one part of this format most likely to be subtly wrong: an
        // atom with no real exclusions is stated as count 1 / partner 0,
        // not count 0 -- reading a literal atom 0 as a real exclusion
        // would silently invent a bond-like relationship to a nonexistent
        // atom.
        let mol = parse_prmtop(EXCLUSIONS).expect("valid PRMTOP");
        let force_field = mol.force_field().expect("has a force field");
        assert!(
            force_field.exclusions.len() >= 3,
            "fixture should exercise more than one atom's worth of exclusions"
        );

        let written = write_prmtop(&mol);
        let back = parse_prmtop(&written).expect("round trips");

        let mut original: Vec<[usize; 2]> = force_field.exclusions.clone();
        let mut round_tripped: Vec<[usize; 2]> = back.force_field().unwrap().exclusions.clone();
        for pair in original.iter_mut().chain(round_tripped.iter_mut()) {
            pair.sort_unstable();
        }
        original.sort_unstable();
        round_tripped.sort_unstable();
        assert_eq!(round_tripped, original);
    }

    #[test]
    fn test_no_coordinates_are_ever_invented() {
        let mol = parse_prmtop(WATER).expect("valid PRMTOP");
        assert!(!mol.has_coords3());
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_element_falls_back_to_mass_when_atomic_number_is_absent() {
        let mol = parse_prmtop(NO_ATOMIC_NUMBER).expect("valid PRMTOP");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.atoms()[0].element().symbol(), "O");
        assert_eq!(mol.atoms()[1].element().symbol(), "H");
    }

    #[test]
    fn test_an_unrecognised_header_is_a_clear_error_not_a_panic() {
        let err = parse_prmtop("NOT PRMTOP AT ALL\n").unwrap_err();
        assert!(matches!(err, PrmtopError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_mass_matching_no_element_is_a_clear_error() {
        let text = "%VERSION  VERSION_STAMP = V0001.000  DATE = 00/00/00  00:00:00\n\
                    %FLAG ATOM_NAME\n%FORMAT(20a4)\nX   \n\
                    %FLAG CHARGE\n%FORMAT(5E16.8)\n  0.00000000E+00\n\
                    %FLAG MASS\n%FORMAT(5E16.8)\n  9.99000000E+02\n\
                    %FLAG RESIDUE_LABEL\n%FORMAT(20a4)\nRES \n\
                    %FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n       1\n";
        let err = parse_prmtop(text).unwrap_err();
        assert!(matches!(err, PrmtopError::InvalidElement { .. }), "{err}");
    }

    #[test]
    fn test_writing_an_improper_never_leaves_atom_zero_in_the_sign_bearing_slot() {
        // AMBER's own constraint: atom index 0 can never sit in a
        // dihedral's 4th (sign-bearing) slot, since `-0` is
        // indistinguishable from `0`. The writer reorders rather than
        // silently mis-encoding such a term as a proper torsion.
        let reordered = improper_atoms_for_write([2, 1, 3, 0]);
        assert_ne!(reordered[3], 0, "atom 0 must not remain in the 4th slot");
        assert_eq!(reordered[0], 0, "the swapped atom lands in the 1st slot");
        assert_eq!(
            reordered
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>(),
            [0, 1, 2, 3].into_iter().collect(),
            "the same four atoms are still named, just reordered"
        );
    }
}
