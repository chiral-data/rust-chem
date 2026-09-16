//! PSF — CHARMM/NAMD's protein structure file: atom types, masses, partial
//! charges, and the explicit bond/angle/dihedral/improper/exclusion/donor/
//! acceptor lists. No coordinates at all (#321) — a PSF is paired with a
//! DCD or a PDB for geometry, never read alone for it.
//!
//! The first real consumer of [`crate::core::force_field::ForceFieldTopology`]
//! (built in #315, never wired to a format until now).
//!
//! **Three on-disk shapes, one extension, unreliable header flags.** The
//! header line is `PSF` plus optional flags (`EXT`, `CHEQ`, `XPLOR`, `CMAP`,
//! `DRUDE`) — but real files lie by omission (a real fetched fixture uses
//! string atom types with a bare `PSF` header, no `XPLOR` flag at all).
//! Every real parser copes by never trusting field *width*: this one
//! whitespace-tokenizes rather than replicating Fortran fixed-column
//! slicing, and reads the atom-type field as a string unconditionally
//! (`ForceFieldAtom::atom_type` is already `Option<String>` — an
//! XPLOR-style integer type code just becomes a numeric-looking string,
//! which is fine).
//!
//! **PSF states no element** — only a name label (`CA`, `OH2`, ambiguous:
//! alpha-carbon in every protein PSF, not calcium) and a force-field type
//! (`CT`, `OT`, equally ambiguous). This module's own `element_from_mass`
//! infers it from the stated mass instead, against the crate's own existing
//! [`crate::core::elements::ATOMIC_MASSES`] table — unambiguous, since PSF
//! always states mass.
//!
//! **`!NNB` is two concatenated, differently-shaped arrays.** After the
//! `!NNB` count and that many exclusion-partner integers (`INB14`), the
//! file continues — with no header line of its own — into exactly `NATOM`
//! more integers: a per-atom cumulative pointer array (`IBLO14`) that
//! carves `INB14` into per-atom exclusion lists. See this module's own
//! `decode_exclusions`/`encode_exclusions`.
//!
//! Bonds go through [`Molecule`]'s ordinary bond list, not
//! `ForceFieldTopology` — bonds are not a force-field-specific concept.
//! `!NBOND` pairs become [`BondOrder::Single`], the same "connectivity
//! stated, order not" precedent [`crate::io::pdb`]'s `CONECT` already
//! established. Per-atom segid/resid/resname group into
//! [`crate::core::residue::Chain`]/[`crate::core::residue::Residue`] via
//! `cif_model`'s own `group_into_chains_and_residues`, the same algorithm
//! mmCIF/CIF-core use, keyed by `(segid, resid, resname)` instead of their
//! own `auth_*` fields.
//!
//! `!NGRP` (charge groups) and any `CMAP`/`DRUDE` sections are skipped by
//! count, never stored — no data-model slot exists for them, the same
//! treatment mmCIF gives `_struct_conn`. The writer emits a minimal
//! placeholder `!NGRP` record for compatibility with other tools reading
//! this crate's own output.
//!
//! A PSF file is always exactly one topology — no multi-model convention
//! exists anywhere in the format.

use std::collections::VecDeque;

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::elements::ATOMIC_MASSES;
use crate::core::force_field::{ForceFieldAtom, ForceFieldTopology};
use crate::core::molecule::Molecule;
use crate::core::site::AtomSite;
use crate::io::cif_model::{ResidueKey, group_into_chains_and_residues};
use crate::io::errors::PsfError;

/// The largest difference (Ångström-scale atomic mass units) between a
/// stated mass and a table entry that still counts as a match. Loose enough
/// to tolerate a force field's own rounding, tight enough that no two
/// elements this crate covers are ever confused for each other.
const MASS_TOLERANCE: f64 = 0.3;

/// Infers an element from a stated atomic mass -- PSF's substitute for an
/// element column it does not have. `None` when nothing in
/// [`ATOMIC_MASSES`] is close enough, rather than guessing.
pub(crate) fn element_from_mass(mass: f64) -> Option<Element> {
    ATOMIC_MASSES
        .iter()
        .enumerate()
        .skip(1) // index 0 is the placeholder "no element" entry
        .map(|(i, &table_mass)| (i, (mass - table_mass).abs()))
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .filter(|(_, diff)| *diff <= MASS_TOLERANCE)
        .and_then(|(i, _)| Element::new(i as u8))
}

/// A flat, line-boundary-agnostic reader over a PSF file's sections.
///
/// Two reading modes, because PSF genuinely mixes them: `!NATOM` is one
/// atom per physical line, never wrapped -- but every bonded/exclusion
/// section is a flat stream of integers that wraps across an arbitrary
/// number of entries per line (confirmed against a real fixture, whose
/// last line of a section is often short). `line()` serves the first,
/// `int()` the second.
struct Sections<'a> {
    lines: std::str::Lines<'a>,
    buffer: VecDeque<&'a str>,
}

impl<'a> Sections<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            lines: text.lines(),
            buffer: VecDeque::new(),
        }
    }

    /// Skips blank lines and returns the next non-blank one, whole.
    fn next_line(&mut self) -> Option<&'a str> {
        loop {
            let line = self.lines.next()?;
            if !line.trim().is_empty() {
                return Some(line);
            }
        }
    }

    /// Reads a section header: skips blank lines, expects the next
    /// non-blank line to contain `marker`, and returns the leading integer
    /// count -- the first whitespace token, whatever else the line says
    /// (`!NGRP`'s header carries a second count this crate does not need,
    /// which this simply never looks at).
    fn header(&mut self, marker: &str) -> Result<usize, PsfError> {
        let line = self
            .next_line()
            .ok_or_else(|| PsfError::ParseError(format!("missing {marker} section")))?;
        if !line.contains(marker) {
            return Err(PsfError::ParseError(format!(
                "expected {marker}, found: {line:?}"
            )));
        }
        line.split_whitespace()
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| PsfError::ParseError(format!("no count before {marker}")))
    }

    /// One raw atom line, for `!NATOM`.
    fn line(&mut self) -> Result<&'a str, PsfError> {
        self.lines
            .next()
            .ok_or_else(|| PsfError::ParseError("truncated !NATOM section".to_string()))
    }

    /// The next integer from the flat, line-wrapping token stream.
    fn int(&mut self) -> Result<i64, PsfError> {
        loop {
            if let Some(token) = self.buffer.pop_front() {
                return token.parse().map_err(|_| {
                    PsfError::ParseError(format!("expected an integer, got {token:?}"))
                });
            }
            let line = self
                .lines
                .next()
                .ok_or_else(|| PsfError::ParseError("unexpected end of section".to_string()))?;
            self.buffer.extend(line.split_whitespace());
        }
    }
}

/// Converts a PSF file's 1-based atom serial into a 0-based index, checked
/// against `num_atoms` rather than trusted.
fn atom_index(serial: i64, num_atoms: usize) -> Result<usize, PsfError> {
    if serial < 1 || serial as usize > num_atoms {
        return Err(PsfError::ParseError(format!(
            "atom serial {serial} is out of range for {num_atoms} atoms"
        )));
    }
    Ok(serial as usize - 1)
}

fn read_terms<const N: usize>(
    sections: &mut Sections,
    count: usize,
    num_atoms: usize,
) -> Result<Vec<[usize; N]>, PsfError> {
    let mut terms = Vec::with_capacity(count);
    for _ in 0..count {
        let mut term = [0usize; N];
        for slot in &mut term {
            *slot = atom_index(sections.int()?, num_atoms)?;
        }
        terms.push(term);
    }
    Ok(terms)
}

/// Reconstructs `ForceFieldTopology::exclusions` from `!NNB`'s two
/// concatenated arrays: `INB14` (the exclusion count and that many partner
/// serials) followed immediately, with no header of its own, by `IBLO14`
/// (`num_atoms` cumulative pointers). Atom `i`'s exclusion partners are the
/// slice of `INB14` between `IBLO14[i-1]` (or `0`) and `IBLO14[i]`.
fn decode_exclusions(
    sections: &mut Sections,
    num_atoms: usize,
) -> Result<Vec<[usize; 2]>, PsfError> {
    let nnb = sections.header("!NNB")?;
    let mut inb14 = Vec::with_capacity(nnb);
    for _ in 0..nnb {
        inb14.push(sections.int()?);
    }
    let mut iblo14 = Vec::with_capacity(num_atoms);
    for _ in 0..num_atoms {
        iblo14.push(sections.int()?);
    }

    let mut exclusions = Vec::new();
    let mut start = 0usize;
    for (atom_ix, &cumulative) in iblo14.iter().enumerate() {
        let end = usize::try_from(cumulative)
            .map_err(|_| PsfError::ParseError(format!("negative IBLO14 pointer {cumulative}")))?;
        if end < start || end > inb14.len() {
            return Err(PsfError::ParseError(format!(
                "IBLO14 pointer {end} is out of range for {} INB14 entries",
                inb14.len()
            )));
        }
        for &partner in &inb14[start..end] {
            exclusions.push([atom_ix, atom_index(partner, num_atoms)?]);
        }
        start = end;
    }
    Ok(exclusions)
}

/// The inverse of [`decode_exclusions`]: builds `INB14`/`IBLO14` from
/// `exclusions`, keyed by the smaller of each pair (CHARMM's own
/// convention, and the one that makes the cumulative-pointer structure
/// well-formed regardless of which order a pair arrived in).
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

    let mut inb14 = Vec::new();
    let mut iblo14 = Vec::with_capacity(num_atoms);
    let mut cumulative = 0i64;
    for list in &partners {
        for &partner in list {
            inb14.push(partner as i64 + 1);
        }
        cumulative += list.len() as i64;
        iblo14.push(cumulative);
    }
    (inb14, iblo14)
}

/// Parses a whole PSF file into a [`Molecule`].
pub fn parse_psf(text: &str) -> Result<Molecule, PsfError> {
    let mut sections = Sections::new(text);

    let header = sections
        .next_line()
        .ok_or_else(|| PsfError::ParseError("empty file".to_string()))?;
    if !header.trim_start().starts_with("PSF") {
        return Err(PsfError::ParseError(format!("not a PSF file: {header:?}")));
    }

    let ntitle = sections.header("!NTITLE")?;
    for _ in 0..ntitle {
        sections.line()?;
    }

    let natom = sections.header("!NATOM")?;
    let mut mol = Molecule::new();
    let mut sites = Vec::with_capacity(natom);
    let mut ff_atoms = Vec::with_capacity(natom);
    let mut keys = Vec::with_capacity(natom);

    for _ in 0..natom {
        let line = sections.line()?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 9 {
            return Err(PsfError::InvalidAtomRow(line.to_string()));
        }
        let segid = fields[1];
        let resid: i32 = fields[2]
            .parse()
            .map_err(|_| PsfError::InvalidAtomRow(line.to_string()))?;
        let resname = fields[3];
        let atomname = fields[4];
        let atomtype = fields[5];
        let charge: f64 = fields[6]
            .parse()
            .map_err(|_| PsfError::InvalidAtomRow(line.to_string()))?;
        let mass: f64 = fields[7]
            .parse()
            .map_err(|_| PsfError::InvalidAtomRow(line.to_string()))?;

        let element = element_from_mass(mass).ok_or(PsfError::InvalidElement(mass))?;
        mol.add_atom(Atom::new(element));
        sites.push(AtomSite {
            name: Some(atomname.to_string()),
            ..AtomSite::empty()
        });
        ff_atoms.push(ForceFieldAtom {
            atom_type: Some(atomtype.to_string()),
            mass: Some(mass),
            partial_charge: Some(charge),
        });
        keys.push(ResidueKey {
            chain_id: segid.to_string(),
            name: resname.to_string(),
            sequence: resid,
            insertion_code: None,
            is_hetero: false,
        });
    }

    mol.set_sites(sites)
        .map_err(|e| PsfError::ParseError(e.to_string()))?;
    let (chains, residues) = group_into_chains_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| PsfError::ParseError(e.to_string()))?;

    let nbond = sections.header("!NBOND")?;
    for _ in 0..nbond {
        let a = atom_index(sections.int()?, natom)?;
        let b = atom_index(sections.int()?, natom)?;
        mol.add_bond(Bond::new(a, b, BondOrder::Single))
            .map_err(|e| PsfError::ParseError(e.to_string()))?;
    }
    // `!NBOND` is PSF's complete, real bond list -- not PDBQT's
    // rotatable-pivots-only info or mmCIF's absent one -- so implying
    // hydrogens from it is legitimate, the same precedent PDB's own
    // `CONECT` reader already established (#285).
    mol.calculate_implicit_hydrogens_where_bonded();

    let ntheta = sections.header("!NTHETA")?;
    let angles = read_terms::<3>(&mut sections, ntheta, natom)?;

    let nphi = sections.header("!NPHI")?;
    let dihedrals = read_terms::<4>(&mut sections, nphi, natom)?;

    let nimphi = sections.header("!NIMPHI")?;
    let impropers = read_terms::<4>(&mut sections, nimphi, natom)?;

    let ndon = sections.header("!NDON")?;
    let donors = read_terms::<2>(&mut sections, ndon, natom)?;

    let nacc = sections.header("!NACC")?;
    let acceptors = read_terms::<2>(&mut sections, nacc, natom)?;

    let exclusions = decode_exclusions(&mut sections, natom)?;

    mol.set_force_field(ForceFieldTopology {
        atoms: Some(ff_atoms),
        angles,
        dihedrals,
        impropers,
        exclusions,
        donors,
        acceptors,
    })
    .map_err(|e| PsfError::ParseError(e.to_string()))?;

    if mol.num_atoms() == 0 {
        return Err(PsfError::NoAtoms);
    }

    Ok(mol)
}

fn push_terms(out: &mut String, header: &str, terms: &[Vec<i64>]) {
    let count: usize = terms.len();
    out.push_str(&format!("{count:>8} {header}\n"));
    let flat: Vec<i64> = terms.iter().flatten().copied().collect();
    for chunk in flat.chunks(8) {
        let row: Vec<String> = chunk.iter().map(|n| format!("{n:>8}")).collect();
        out.push_str(&row.join(""));
        out.push('\n');
    }
    out.push('\n');
}

fn one_based_terms<const N: usize>(terms: &[[usize; N]]) -> Vec<Vec<i64>> {
    terms
        .iter()
        .map(|term| term.iter().map(|&ix| ix as i64 + 1).collect())
        .collect()
}

/// Writes a [`Molecule`] as a PSF file.
///
/// A bare `PSF` header, whitespace-separated (not fixed-width) columns --
/// matching what real, modern-generated files already look like, and
/// consistent with every real parser whitespace-tokenizing regardless of
/// declared flags (see the module doc).
pub fn write_psf(mol: &Molecule) -> String {
    let num_atoms = mol.num_atoms();
    let empty_ff = ForceFieldTopology::default();
    let force_field = mol.force_field().unwrap_or(&empty_ff);

    let mut out = String::from("PSF\n\n");
    out.push_str("       1 !NTITLE\n REMARKS written by chem\n\n");

    out.push_str(&format!("{num_atoms:>8} !NATOM\n"));
    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let residue = mol.residue_of(i);
        let chain = mol.chain_of(i);
        let ff_atom = force_field.atoms.as_ref().and_then(|atoms| atoms.get(i));

        let segid = chain
            .map(|c| c.id.as_str())
            .filter(|s| !s.is_empty())
            .unwrap_or("SEG");
        let resid = residue.map(|r| r.sequence).unwrap_or(1);
        let resname = residue
            .map(|r| r.name.as_str())
            .filter(|s| !s.is_empty())
            .unwrap_or("UNK");
        let atomname = site
            .and_then(|s| s.name.as_deref())
            .filter(|s| !s.is_empty())
            .unwrap_or(atom.element().symbol());
        let atomtype = ff_atom
            .and_then(|a| a.atom_type.as_deref())
            .filter(|s| !s.is_empty())
            .unwrap_or(atom.element().symbol());
        let charge = ff_atom.and_then(|a| a.partial_charge).unwrap_or(0.0);
        let mass = ff_atom
            .and_then(|a| a.mass)
            .unwrap_or(ATOMIC_MASSES[atom.element().atomic_number as usize]);

        out.push_str(&format!(
            "{:>8} {segid:<4} {resid:<4} {resname:<4} {atomname:<4} {atomtype:<4} {charge:>10.6} {mass:>10.4} {:>8}\n",
            i + 1,
            0,
        ));
    }
    out.push('\n');

    let bond_count = mol.bonds().len();
    out.push_str(&format!("{bond_count:>8} !NBOND: bonds\n"));
    let bond_pairs: Vec<i64> = mol
        .bonds()
        .iter()
        .flat_map(|b| [b.atom1() as i64 + 1, b.atom2() as i64 + 1])
        .collect();
    for chunk in bond_pairs.chunks(8) {
        let row: Vec<String> = chunk.iter().map(|n| format!("{n:>8}")).collect();
        out.push_str(&row.join(""));
        out.push('\n');
    }
    out.push('\n');

    push_terms(
        &mut out,
        "!NTHETA: angles",
        &one_based_terms(&force_field.angles),
    );
    push_terms(
        &mut out,
        "!NPHI: dihedrals",
        &one_based_terms(&force_field.dihedrals),
    );
    push_terms(
        &mut out,
        "!NIMPHI: impropers",
        &one_based_terms(&force_field.impropers),
    );
    push_terms(
        &mut out,
        "!NDON: donors",
        &one_based_terms(&force_field.donors),
    );
    push_terms(
        &mut out,
        "!NACC: acceptors",
        &one_based_terms(&force_field.acceptors),
    );

    let (inb14, iblo14) = encode_exclusions(&force_field.exclusions, num_atoms);
    out.push_str(&format!("{:>8} !NNB\n", inb14.len()));
    for chunk in inb14.chunks(8) {
        let row: Vec<String> = chunk.iter().map(|n| format!("{n:>8}")).collect();
        out.push_str(&row.join(""));
        out.push('\n');
    }
    out.push('\n');
    for chunk in iblo14.chunks(8) {
        let row: Vec<String> = chunk.iter().map(|n| format!("{n:>8}")).collect();
        out.push_str(&row.join(""));
        out.push('\n');
    }
    out.push('\n');

    // A minimal placeholder -- charge groups have no data-model slot (see
    // the module doc), but a real `!NGRP` record's presence is what some
    // other tools expect structurally.
    out.push_str("       1       0 !NGRP\n       0       0       0\n");

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = include_str!("../../tests/corpus/psf/water.psf");
    const EXCLUSIONS: &str = include_str!("../../tests/corpus/psf/exclusions.psf");

    #[test]
    fn test_element_from_mass_matches_common_elements() {
        assert_eq!(element_from_mass(1.008).map(|e| e.symbol()), Some("H"));
        assert_eq!(element_from_mass(12.011).map(|e| e.symbol()), Some("C"));
        assert_eq!(element_from_mass(15.9994).map(|e| e.symbol()), Some("O"));
        assert_eq!(element_from_mass(14.007).map(|e| e.symbol()), Some("N"));
    }

    #[test]
    fn test_element_from_mass_refuses_a_mass_matching_nothing() {
        assert_eq!(element_from_mass(999.0), None);
    }

    #[test]
    fn test_a_water_topology_round_trips() {
        let mol = parse_psf(WATER).expect("valid PSF");
        assert_eq!(mol.num_atoms(), 3);
        assert!(!mol.has_coords3(), "PSF states no coordinates");
        assert_eq!(mol.num_bonds(), 2);
        assert_eq!(mol.chains().len(), 1);
        assert_eq!(mol.residues().len(), 1);
        assert_eq!(mol.residues()[0].name, "TIP3");

        let force_field = mol.force_field().expect("has a force field");
        assert_eq!(force_field.angles.len(), 1);
        assert_eq!(
            force_field.atoms.as_ref().unwrap()[0].atom_type.as_deref(),
            Some("OT")
        );

        let written = write_psf(&mol);
        let back = parse_psf(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.num_bonds(), mol.num_bonds());
        assert_eq!(
            back.force_field().unwrap().angles,
            mol.force_field().unwrap().angles
        );
    }

    #[test]
    fn test_the_nnb_reconstruction_survives_a_round_trip() {
        // The one part of this format most likely to be subtly wrong: an
        // off-by-one in IBLO14's cumulative pointer, or reading it before
        // INB14, both parse "successfully" while getting every exclusion
        // pair after the first atom wrong.
        let mol = parse_psf(EXCLUSIONS).expect("valid PSF");
        let force_field = mol.force_field().expect("has a force field");
        assert!(
            force_field.exclusions.len() >= 3,
            "fixture should exercise more than one atom's worth of exclusions"
        );

        let written = write_psf(&mol);
        let back = parse_psf(&written).expect("round trips");

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
    fn test_donors_and_acceptors_survive_a_round_trip() {
        let mol = parse_psf(WATER).expect("valid PSF");
        let force_field = mol.force_field().unwrap();
        assert!(!force_field.donors.is_empty());
        assert!(!force_field.acceptors.is_empty());

        let written = write_psf(&mol);
        let back = parse_psf(&written).expect("round trips");
        assert_eq!(back.force_field().unwrap().donors, force_field.donors);
        assert_eq!(back.force_field().unwrap().acceptors, force_field.acceptors);
    }

    #[test]
    fn test_no_coordinates_are_ever_invented() {
        let mol = parse_psf(WATER).expect("valid PSF");
        assert!(!mol.has_coords3());
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_an_unrecognised_header_is_a_clear_error_not_a_panic() {
        let err = parse_psf("NOT PSF AT ALL\n").unwrap_err();
        assert!(matches!(err, PsfError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_mass_matching_no_element_is_a_clear_error() {
        let text = "PSF\n\n       1 !NTITLE\n REMARKS\n\n       1 !NATOM\n\
                    1 A 1 RES NM X 0.0 999.0 0\n\n       0 !NBOND: bonds\n\n\n\
                    0 !NTHETA: angles\n\n\n0 !NPHI: dihedrals\n\n\n\
                    0 !NIMPHI: impropers\n\n\n0 !NDON: donors\n\n\n\
                    0 !NACC: acceptors\n\n\n0 !NNB\n\n0\n";
        let err = parse_psf(text).unwrap_err();
        assert!(matches!(err, PsfError::InvalidElement(_)), "{err}");
    }
}
