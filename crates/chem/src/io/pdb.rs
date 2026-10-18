//! PDB — fixed-column text, `ATOM`/`HETATM` per atom, `CRYST1` for the unit
//! cell, `CONECT` for explicit bonds (#223).
//!
//! **By default no bonds are synthesised beyond what `CONECT` states.** Real
//! PDB files leave standard protein/nucleic-acid bonding to residue
//! templates rather than stating it. With
//! [`PdbReadOptions::residue_templates`] on (#403), standard residues are
//! bonded from [`crate::core::residue_templates`] after `CONECT`, then the
//! `SSBOND` and `LINK` records are bonded too. A `LINK` to a metal, or across
//! a symmetry operator, is skipped: it is coordination or a contact with a
//! crystal mate, not a bond inside this model. Under `MODEL` splitting these
//! header records reach the first model only, so in
//! [`MultiFrameMode::Molecules`](crate::io::options::MultiFrameMode::Molecules)
//! later models get templates but no links.
//!
//! [`PdbReadOptions::perceive_bonds`] runs last, after `CONECT` and the
//! templates: it bonds by distance only the pairs still unbonded (see below).
//!
//! Hydrogen counts follow from that. PDB names no count of its own, so an
//! atom arrives with none, and
//! [`Molecule::calculate_implicit_hydrogens`] cannot be used to fill it in:
//! it sums the bond orders around an atom, and with atoms left bondless it
//! would hand each one its *free-atom* count — a backbone
//! carbon reading as methane. What is called instead is
//! [`Molecule::calculate_implicit_hydrogens_where_bonded`], which fills only
//! the atoms some bond describes and leaves the rest honestly blank (#285).
//! Before it, a fully-`CONECT`ed benzene wrote back as
//! `[C]1[C][C][C][C][C]1` — six carbons stated to have no hydrogens, which is
//! not the molecule the file described.
//!
//! With [`PdbReadOptions::perceive_bonds`] set, atom pairs `CONECT` left
//! unbonded are bonded by distance ([`crate::core::perceive::perceive_bonds`])
//! before that fill, so the perceived bonds count toward it too.
//!
//! `CRYST1`'s orientation convention already matches this crate's fixed
//! [`crate::core::cell::UnitCell`] (its own doc comment: "the near-universal
//! crystallographic setting... which is what PDB, CIF and essentially every
//! toolkit assume") — unlike extended XYZ's `Lattice=` (#222, arbitrary
//! per-writer orientation, deliberately not read), so `CRYST1` is read
//! directly with no rotation step and no correctness risk.
//!
//! **A multi-structure write frames every record in `MODEL`/`ENDMDL`**
//! ([`frame_models`], used by `write_pdb_records` and
//! [`crate::io::supplier::PdbWriter`]) — `ENDMDL` is the boundary the reader
//! splits on, and structures concatenated without it read back as one
//! molecule with every atom merged (#267). One structure is written
//! unframed, as a real single deposition is. OpenBabel frames the same way,
//! and gemmi reads three models out of the result.
//!
//! A file may hold more than one structure via `MODEL`/`ENDMDL` (an NMR
//! ensemble); splitting those frames is the reader's/supplier's job
//! ([`crate::io::reader::read_pdb_with_options`],
//! [`crate::io::supplier::PdbSupplier`]), the same division SDF's `$$$$`
//! and XYZ's atom-count boundary already have. This module parses and
//! writes exactly one already-isolated frame.

use std::collections::HashMap;

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::cell::{SpaceGroup, UnitCell};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::{Point3, is_placeholder_3d};
use crate::core::molecule::Molecule;
use crate::core::perceive::perceive_bonds;
use crate::core::residue::{Chain, Residue};
use crate::core::residue_templates::{add_residue_template_bonds, bond_conformer_pairs};
use crate::core::site::AtomSite;
use crate::io::errors::PdbError;
use crate::io::options::PdbReadOptions;

/// 1-based, inclusive column range, matching the spec's own numbering.
/// `None` for a line too short to have this field, or a wholly blank one —
/// real files routinely omit trailing optional columns.
fn column(line: &str, start_1based: usize, end_1based: usize) -> Option<&str> {
    let bytes = line.as_bytes();
    if start_1based == 0 || start_1based > bytes.len() {
        return None;
    }
    let start = start_1based - 1;
    let end = end_1based.min(bytes.len());
    line.get(start..end).filter(|s| !s.trim().is_empty())
}

fn parse_f64_column(line: &str, start: usize, end: usize) -> Result<f64, PdbError> {
    column(line, start, end)
        .map(str::trim)
        .ok_or_else(|| PdbError::InvalidAtomLine(line.to_string()))?
        .parse()
        .map_err(|_| PdbError::InvalidAtomLine(line.to_string()))
}

/// Looks a symbol up case-insensitively against [`ELEMENT_SYMBOLS`]'
/// canonical mixed-case spelling (`"Ca"`, `"Fe"`) — real PDB files write
/// their element/atom-name columns in upper case throughout.
fn element_from_symbol(sym: &str) -> Option<Element> {
    let sym = sym.trim();
    let mut chars = sym.chars();
    let normalised = match (chars.next(), chars.next()) {
        (Some(a), Some(b)) if chars.next().is_none() => {
            format!("{}{}", a.to_ascii_uppercase(), b.to_ascii_lowercase())
        }
        (Some(a), None) => a.to_ascii_uppercase().to_string(),
        _ => return None,
    };
    ELEMENT_SYMBOLS
        .iter()
        .position(|&s| s == normalised)
        .and_then(|n| Element::new(n as u8))
}

/// Resolves an atom's element. Prefers the dedicated column (77-78, added
/// to the spec in a later revision and absent from many older files); when
/// absent, falls back to the *positional* convention real parsers use on
/// the atom-name field itself — column 13 (the name field's first
/// character) blank means a single-letter element right in column 14
/// (`" CA "` → carbon), non-blank means a two-letter element occupying
/// columns 13-14 (`"CA  "` → calcium). A documented column-position rule,
/// not a content-based guess — the same "count columns, don't guess"
/// discipline #202 already established for SDF's atom-count field.
fn resolve_element(element_col: Option<&str>, name_field: &str) -> Result<Element, PdbError> {
    if let Some(sym) = element_col {
        return element_from_symbol(sym)
            .ok_or_else(|| PdbError::InvalidElement(sym.trim().to_string()));
    }

    let first_col_blank = name_field.chars().next().is_none_or(|c| c == ' ');
    let mut chars = name_field.chars();
    let candidate = if first_col_blank {
        chars.next(); // the blank column itself
        chars.next().map(|c| c.to_string())
    } else {
        match (chars.next(), chars.next()) {
            (Some(a), Some(b)) if b != ' ' => Some(format!("{a}{b}")),
            (Some(a), _) => Some(a.to_string()),
            _ => None,
        }
    };

    candidate
        .as_deref()
        .and_then(element_from_symbol)
        .ok_or_else(|| PdbError::InvalidElement(name_field.trim().to_string()))
}

/// The fields that define a residue's identity, one per atom in file order
/// — collected during the main pass, then grouped into [`Residue`]/[`Chain`]
/// in a second pass. Keeping the two passes separate is what makes the
/// grouping logic simple: a residue is just a maximal run of atoms sharing
/// one of these.
#[derive(Clone, PartialEq, Eq)]
struct ResidueKey {
    chain_id: String,
    name: String,
    sequence: i32,
    insertion_code: Option<char>,
    is_hetero: bool,
}

fn group_into_chains_and_residues(keys: &[ResidueKey]) -> (Vec<Chain>, Vec<Residue>) {
    let mut chains = Vec::new();
    let mut residues = Vec::new();
    let mut current_chain_id: Option<&str> = None;
    let mut chain_start = 0;

    let mut i = 0;
    while i < keys.len() {
        let key = &keys[i];
        let start = i;
        while i < keys.len() && keys[i] == *key {
            i += 1;
        }

        if current_chain_id != Some(key.chain_id.as_str()) {
            if let Some(id) = current_chain_id {
                chains.push(Chain {
                    id: id.to_string(),
                    label_id: None,
                    residues: chain_start..residues.len(),
                });
            }
            current_chain_id = Some(&key.chain_id);
            chain_start = residues.len();
        }

        residues.push(Residue {
            name: key.name.clone(),
            sequence: key.sequence,
            insertion_code: key.insertion_code,
            label_seq: None,
            chain_ix: chains.len(),
            is_hetero: key.is_hetero,
            atoms: start..i,
        });
    }
    if let Some(id) = current_chain_id {
        chains.push(Chain {
            id: id.to_string(),
            label_id: None,
            residues: chain_start..residues.len(),
        });
    }

    (chains, residues)
}

/// One end of an `SSBOND` or `LINK` record.
struct LinkEnd {
    atom: String,
    alt_loc: Option<char>,
    residue: String,
    chain_id: String,
    sequence: i32,
    insertion_code: Option<char>,
}

/// A `SSBOND`/`LINK` pair, or `None` for one across a symmetry operator,
/// whose partner is a crystal mate rather than an atom in this file.
fn parse_link(line: &str, is_ssbond: bool) -> Option<(LinkEnd, LinkEnd)> {
    let field = |a, b| column(line, a, b).map(str::trim);
    let (sym1, sym2) = (field(60, 65), field(67, 72));
    if sym1.is_some() && sym2.is_some() && sym1 != sym2 {
        return None;
    }
    let end = |name: Option<&str>, alt: usize, res: usize, chain: usize, seq: usize| {
        Some(LinkEnd {
            atom: name?.to_string(),
            alt_loc: column(line, alt, alt).and_then(|s| s.chars().next()),
            residue: field(res, res + 2)?.to_string(),
            chain_id: field(chain, chain).unwrap_or("").to_string(),
            sequence: field(seq, seq + 3)?.parse().ok()?,
            insertion_code: column(line, seq + 4, seq + 4).and_then(|s| s.chars().next()),
        })
    };
    if is_ssbond {
        // SSBOND has no atom or alt-loc columns: the partners are always SG.
        Some((
            end(Some("SG"), 0, 12, 16, 18)?,
            end(Some("SG"), 0, 26, 30, 32)?,
        ))
    } else {
        Some((
            end(field(13, 16), 17, 18, 22, 23)?,
            end(field(43, 46), 47, 48, 52, 53)?,
        ))
    }
}

/// The atoms one link end names, every alternate location of it.
fn resolve_link_end(mol: &Molecule, end: &LinkEnd) -> Vec<usize> {
    let chains = mol.chains();
    let sites = mol.sites().unwrap_or(&[]);
    mol.residues()
        .iter()
        .filter(|r| {
            r.name == end.residue
                && r.sequence == end.sequence
                && r.insertion_code == end.insertion_code
                && chains[r.chain_ix].id == end.chain_id
        })
        .flat_map(|r| r.atoms.clone())
        .filter(|&i| {
            let site = &sites[i];
            site.name.as_deref() == Some(end.atom.as_str())
                && (end.alt_loc.is_none() || site.alt_loc.is_none() || site.alt_loc == end.alt_loc)
        })
        .collect()
}

/// Metals coordinate rather than bond covalently; PDB's `LINK` does not say
/// which it is, so the element decides.
fn is_metal(atomic_number: u8) -> bool {
    !matches!(
        atomic_number,
        1 | 2 | 5..=10 | 14..=18 | 32..=36 | 51..=54 | 85 | 86
    )
}

/// Parses one PDB structure (a single `MODEL`, or a whole file with none).
pub fn parse_pdb(text: &str) -> Result<Molecule, PdbError> {
    parse_pdb_with_options(text, &PdbReadOptions::default())
}

/// [`parse_pdb`] with PDB's own read options.
pub fn parse_pdb_with_options(text: &str, options: &PdbReadOptions) -> Result<Molecule, PdbError> {
    let mut mol = Molecule::new();
    let mut links: Vec<(LinkEnd, LinkEnd)> = Vec::new();
    let mut sites = Vec::new();
    let mut coords = Vec::new();
    let mut keys = Vec::new();
    let mut serial_to_index: HashMap<i64, usize> = HashMap::new();
    let mut conect_pairs: Vec<(usize, usize)> = Vec::new();
    let mut cell: Option<UnitCell> = None;
    let mut space_group: Option<SpaceGroup> = None;

    for line in text.lines() {
        if line.len() < 6 {
            continue;
        }
        match line[0..6].trim() {
            "ATOM" | "HETATM" => {
                let is_hetero = line[0..6].trim() == "HETATM";
                let serial: i64 = column(line, 7, 11)
                    .map(str::trim)
                    .ok_or_else(|| PdbError::InvalidAtomLine(line.to_string()))?
                    .parse()
                    .map_err(|_| PdbError::InvalidAtomLine(line.to_string()))?;
                let name_field = line.get(12..16.min(line.len())).unwrap_or("");
                let name = name_field.trim().to_string();
                let alt_loc = column(line, 17, 17).and_then(|s| s.chars().next());
                let res_name = column(line, 18, 20).unwrap_or("").trim().to_string();
                let chain_id = column(line, 22, 22).unwrap_or("").trim().to_string();
                let res_seq: i32 = column(line, 23, 26)
                    .map(str::trim)
                    .ok_or_else(|| PdbError::InvalidAtomLine(line.to_string()))?
                    .parse()
                    .map_err(|_| PdbError::InvalidAtomLine(line.to_string()))?;
                let insertion_code = column(line, 27, 27).and_then(|s| s.chars().next());
                let x = parse_f64_column(line, 31, 38)?;
                let y = parse_f64_column(line, 39, 46)?;
                let z = parse_f64_column(line, 47, 54)?;
                let occupancy = column(line, 55, 60).and_then(|s| s.trim().parse().ok());
                let b_factor = column(line, 61, 66).and_then(|s| s.trim().parse().ok());
                let element_col = column(line, 77, 78);

                let element = resolve_element(element_col, name_field)?;
                let atom_idx = mol.add_atom(Atom::new(element));
                serial_to_index.insert(serial, atom_idx);
                coords.push(Point3::new(x, y, z));
                sites.push(AtomSite {
                    name: if name.is_empty() { None } else { Some(name) },
                    alt_loc,
                    partial_charge: None,
                    occupancy,
                    b_factor,
                    radius: None,
                });
                keys.push(ResidueKey {
                    chain_id,
                    name: res_name,
                    sequence: res_seq,
                    insertion_code,
                    is_hetero,
                });
            }
            "CRYST1" => {
                let a = parse_f64_column(line, 7, 15)?;
                let b = parse_f64_column(line, 16, 24)?;
                let c = parse_f64_column(line, 25, 33)?;
                let alpha = parse_f64_column(line, 34, 40)?;
                let beta = parse_f64_column(line, 41, 47)?;
                let gamma = parse_f64_column(line, 48, 54)?;
                cell = Some(UnitCell::new(a, b, c, alpha, beta, gamma));
                if let Some(sym) = column(line, 56, 66) {
                    space_group = Some(SpaceGroup::from_symbol(sym.trim()));
                }
            }
            "SSBOND" | "LINK" => {
                links.extend(parse_link(line, line.starts_with("SSBOND")));
            }
            "CONECT" => {
                let src_serial: i64 = column(line, 7, 11)
                    .map(str::trim)
                    .ok_or_else(|| PdbError::ParseError(line.to_string()))?
                    .parse()
                    .map_err(|_| PdbError::ParseError(line.to_string()))?;
                let src_idx = *serial_to_index.get(&src_serial).ok_or_else(|| {
                    PdbError::ParseError(format!(
                        "CONECT references unknown atom serial {src_serial}"
                    ))
                })?;
                for (start, end) in [(12, 16), (17, 21), (22, 26), (27, 31)] {
                    let Some(field) = column(line, start, end) else {
                        continue;
                    };
                    let bonded_serial: i64 = field.trim().parse().map_err(|_| {
                        PdbError::ParseError(format!("invalid CONECT serial {field:?}"))
                    })?;
                    let bonded_idx = *serial_to_index.get(&bonded_serial).ok_or_else(|| {
                        PdbError::ParseError(format!(
                            "CONECT references unknown atom serial {bonded_serial}"
                        ))
                    })?;
                    conect_pairs.push((src_idx.min(bonded_idx), src_idx.max(bonded_idx)));
                }
            }
            _ => {
                // HEADER, TITLE, SEQRES, HELIX, SHEET, REMARK, MODEL,
                // ENDMDL, TER, END and everything else: out of scope for
                // this story, and MODEL/ENDMDL framing is the caller's job
                // (see the module doc).
            }
        }
    }

    // Zeros in these columns are what a format with no room to say "unknown"
    // writes for a molecule that has no conformer, so believing them back is
    // how a converted molecule ended up undrawable (#270).
    if !is_placeholder_3d(&coords) {
        mol.set_coords3(coords)
            .map_err(|e| PdbError::ParseError(e.to_string()))?;
    }
    mol.set_sites(sites)
        .map_err(|e| PdbError::ParseError(e.to_string()))?;

    let (chains, residues) = group_into_chains_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| PdbError::ParseError(e.to_string()))?;

    if let Some(cell) = cell {
        mol.set_cell(cell)
            .map_err(|e| PdbError::ParseError(e.to_string()))?;
    }
    if let Some(sg) = space_group {
        mol.set_space_group(sg);
    }

    conect_pairs.sort_unstable();
    conect_pairs.dedup();
    for (a, b) in conect_pairs {
        mol.add_bond(Bond::new(a, b, BondOrder::Single))
            .map_err(|e| PdbError::ParseError(e.to_string()))?;
    }

    if options.residue_templates {
        add_residue_template_bonds(&mut mol);
        for (first, second) in &links {
            let left = resolve_link_end(&mol, first);
            let right = resolve_link_end(&mol, second);
            let metal = left
                .iter()
                .chain(&right)
                .any(|&i| is_metal(mol.atom(i).atomic_number()));
            if !metal {
                bond_conformer_pairs(&mut mol, &left, &right, BondOrder::Single);
            }
        }
    }

    if mol.num_atoms() == 0 {
        return Err(PdbError::NoAtoms);
    }

    if options.perceive_bonds {
        perceive_bonds(&mut mol);
    }

    // Only where CONECT said what an atom is bonded to. The unguarded form
    // cannot be used here -- with most atoms bondless by construction (see the
    // module doc) it would give each one its free-atom count -- but refusing to
    // fill any was its own wrong answer: every atom kept `hydrogens: None`, so
    // the SMILES writer compared 0 against the implied count and bracketed the
    // lot, and ethanol from a PDB read back as `[C][C][O]` (#285).
    //
    // A bondless atom still keeps `None`. That is not the same claim as
    // `Some(0)`: the file stated no connectivity for it, so its hydrogens are
    // unknown rather than absent.
    mol.calculate_implicit_hydrogens_where_bonded();

    Ok(mol)
}

/// `name` padded into the 4-column atom-name field. A best-effort
/// approximation of the convention most tools use (a single-letter element
/// gets a leading blank in column 13), not exact for every historical
/// naming quirk -- but this crate's own reader never depends on it, since
/// this writer always also emits the dedicated element column (77-78).
fn format_atom_name(name: &str, element_symbol: &str) -> String {
    if element_symbol.len() == 1 && name.len() < 4 {
        format!(" {name:<3}")
    } else {
        format!("{name:<4}")
    }
}

/// One structure inside `MODEL n`/`ENDMDL`.
///
/// What [`frame_models`] does per entry, exposed for the streaming writers,
/// which frame as they go rather than over a slice they can see all of.
pub fn frame_model(serial: usize, body: &str) -> String {
    let mut out = format!("MODEL     {serial:>4}\n");
    // `END` terminates a *file*; inside a model `ENDMDL` does, so the one the
    // single-structure writer appended comes off.
    for line in body.lines() {
        if line.trim_end() != "END" {
            out.push_str(line);
            out.push('\n');
        }
    }
    out.push_str("ENDMDL\n");
    out
}

/// Wraps each of several structures in `MODEL`/`ENDMDL`, or returns the one
/// structure untouched.
///
/// Both readers here split records on `ENDMDL`, and both writers used to
/// concatenate without it -- so `chem convert three.smi --to pdb` wrote three
/// structures that read back as one molecule with all their atoms merged
/// (#267). OpenBabel frames the same way and gemmi reads three models out of
/// the result, so the convention is borrowed rather than invented.
///
/// **One structure is left exactly as it was**, unframed: a real
/// single-structure PDB carries no `MODEL`, and every existing fixture and test
/// would otherwise move. CML wraps unconditionally (#228) for a reason that
/// does not apply here -- sibling `<molecule>` elements with no root are not
/// valid XML, whereas an unframed PDB structure is perfectly valid.
///
/// **No trailing `END` after the last `ENDMDL`.** Measured: a bare `END` there
/// reads back as an extra, empty record, the artefact visible in OpenBabel's
/// own multi-model output.
pub fn frame_models(structures: &[String]) -> String {
    if structures.len() < 2 {
        return structures.first().cloned().unwrap_or_default();
    }

    let mut out = String::new();
    for (i, structure) in structures.iter().enumerate() {
        out.push_str(&frame_model(i + 1, structure));
    }
    out
}

/// Writes one PDB structure.
pub fn write_pdb(mol: &Molecule) -> String {
    let mut out = String::new();

    if let Some(cell) = mol.cell() {
        let sg = mol
            .space_group()
            .and_then(|g| g.symbol.as_deref())
            .unwrap_or("P 1");
        out.push_str(&format!(
            "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} {:<11}\n",
            cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma, sg
        ));
    }

    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let residue = mol.residue_of(i);
        let chain = mol.chain_of(i);
        let p = mol.coord3(i).unwrap_or(Point3::new(0.0, 0.0, 0.0));

        let is_hetero = residue.map(|r| r.is_hetero).unwrap_or(false);
        let record: &str = if is_hetero { "HETATM" } else { "ATOM" };
        let res_name = residue.map(|r| r.name.as_str()).unwrap_or("UNK");
        let res_seq = residue.map(|r| r.sequence).unwrap_or(1);
        let icode = residue.and_then(|r| r.insertion_code).unwrap_or(' ');
        let chain_id = chain.and_then(|c| c.id.chars().next()).unwrap_or(' ');
        let alt_loc = site.and_then(|s| s.alt_loc).unwrap_or(' ');
        let element = atom.element().symbol();
        let name = site
            .and_then(|s| s.name.as_deref())
            .unwrap_or(element)
            .to_string();
        let occupancy = site.and_then(|s| s.occupancy).unwrap_or(1.0);
        let b_factor = site.and_then(|s| s.b_factor).unwrap_or(0.0);

        out.push_str(&format!(
            "{record:<6}{serial:>5} {name4}{alt_loc}{res_name:<3} {chain_id}{res_seq:>4}{icode}   {x:>8.3}{y:>8.3}{z:>8.3}{occupancy:>6.2}{b_factor:>6.2}          {element:>2}\n",
            serial = i + 1,
            name4 = format_atom_name(&name, element),
            x = p.x,
            y = p.y,
            z = p.z,
        ));
    }
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); mol.num_atoms()];
    for bond in mol.bonds() {
        neighbors[bond.atom1()].push(bond.atom2());
        neighbors[bond.atom2()].push(bond.atom1());
    }
    for (atom_idx, mut list) in neighbors.into_iter().enumerate() {
        if list.is_empty() {
            continue;
        }
        list.sort_unstable();
        for chunk in list.chunks(4) {
            let mut line = format!("CONECT{:>5}", atom_idx + 1);
            for &n in chunk {
                line.push_str(&format!("{:>5}", n + 1));
            }
            out.push_str(&line);
            out.push('\n');
        }
    }

    // Last, because `END` terminates the file: a reader that honours it never
    // sees connectivity written after it, and every PDB this crate wrote until
    // #263 put `CONECT` there -- OpenBabel read ethanol back as `C.C.O`. Our
    // own parser ignores `END` outright, so a round trip through it could not
    // have caught this.
    out.push_str("END\n");

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_PDB: &str = "\
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  H1  HOH A   1       0.759   0.000   0.504  1.00  0.00           H
HETATM    3  H2  HOH A   1       0.759   0.000  -0.504  1.00  0.00           H
END
";

    #[test]
    fn test_a_single_structure_round_trips() {
        let mol = parse_pdb(WATER_PDB).expect("valid PDB");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.residues().len(), 1);
        assert!(mol.residues()[0].is_hetero);
        assert_eq!(mol.residues()[0].name, "HOH");
        assert_eq!(mol.chains().len(), 1);
        assert_eq!(mol.chains()[0].id, "A");
        assert_eq!(mol.site(0).unwrap().occupancy, Some(1.0));

        let written = write_pdb(&mol);
        let back = parse_pdb(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.residues()[0].name, "HOH");
        assert_eq!(back.chains()[0].id, "A");
    }

    #[test]
    fn test_cryst1_round_trips_through_unit_cell_and_space_group() {
        let text = format!(
            "CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 21 21 21\n{WATER_PDB}"
        );
        let mol = parse_pdb(&text).expect("valid PDB");
        let cell = mol.cell().expect("has a cell");
        assert_eq!((cell.a, cell.b, cell.c), (10.0, 20.0, 30.0));
        assert_eq!(
            mol.space_group().and_then(|g| g.symbol.as_deref()),
            Some("P 21 21 21")
        );

        let written = write_pdb(&mol);
        assert!(written.starts_with("CRYST1"), "{written}");
        let back = parse_pdb(&written).expect("round trips");
        assert_eq!(back.cell(), mol.cell());
    }

    #[test]
    fn test_conect_bonds_land_on_the_right_atoms_and_dedup() {
        // Both directions of the O-H1 bond, deliberately -- real files
        // often list both, and this must not create the bond twice.
        let text = format!("{WATER_PDB}CONECT    1    2\nCONECT    2    1\n");
        let mol = parse_pdb(&text).expect("valid PDB");
        assert_eq!(mol.num_bonds(), 1);
        let bond = &mol.bonds()[0];
        assert_eq!((bond.atom1(), bond.atom2()), (0, 1));
    }

    #[test]
    fn test_conect_bonds_survive_a_round_trip() {
        let text = format!("{WATER_PDB}CONECT    1    2\nCONECT    1    3\n");
        let mol = parse_pdb(&text).expect("valid PDB");
        assert_eq!(mol.num_bonds(), 2);

        let written = write_pdb(&mol);
        assert!(written.contains("CONECT"), "{written}");
        let back = parse_pdb(&written).expect("round trips");
        assert_eq!(back.num_bonds(), 2);
    }

    /// Distinct per-atom occupancies and B-factors, none of them the value a
    /// missing site would produce.
    ///
    /// [`WATER_PDB`] cannot serve here: its B-factors are all `0.00`, which is
    /// exactly what `write_pdb` emits for an atom with no site at all, so a
    /// writer that dropped the column entirely would round-trip through it
    /// perfectly.
    const SITES_PDB: &str = "\
ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 42.50           N
ATOM      2  CA  ALA A   1      12.560  13.207  10.000  0.80 37.25           C
ATOM      3  C   ALA A   1      13.100  14.600  10.000  0.55 55.00           C
END
";

    #[test]
    fn test_end_is_the_last_record_so_conect_is_not_stranded_after_it() {
        // #263. `END` terminates the file, so connectivity written after it is
        // invisible to any reader that honours it -- OpenBabel read ethanol
        // back as three disconnected atoms. Asserted on the *text*, because a
        // round trip through our own parser is precisely what hid this: it
        // ignores `END` along with every other record it does not model.
        let text = format!("{WATER_PDB}CONECT    1    2\nCONECT    1    3\n");
        let mol = parse_pdb(&text).expect("valid PDB");
        let written = write_pdb(&mol);

        let records: Vec<&str> = written.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(
            records.last(),
            Some(&"END"),
            "END must terminate the file: {written}"
        );

        let end = records.iter().position(|l| *l == "END").expect("an END");
        let last_conect = records
            .iter()
            .rposition(|l| l.starts_with("CONECT"))
            .expect("this molecule has bonds");
        assert!(
            last_conect < end,
            "CONECT after END is unreachable to a reader: {written}"
        );
    }

    #[test]
    fn test_one_structure_is_written_unframed() {
        // A real single-deposition PDB carries no `MODEL`, and framing one
        // would move every fixture in the crate. CML wraps unconditionally
        // (#228) because sibling elements with no root are not valid XML; an
        // unframed PDB structure is perfectly valid, so the reason does not
        // carry over (#267).
        let mol = parse_pdb(WATER_PDB).expect("valid PDB");
        let one = write_pdb(&mol);
        assert_eq!(frame_models(std::slice::from_ref(&one)), one);
        assert!(!one.contains("MODEL"), "{one}");
    }

    #[test]
    fn test_several_structures_are_framed_and_read_back_as_several() {
        // The defect: `ENDMDL` is the boundary the reader splits on, so
        // concatenated structures came back as one merged molecule (#267).
        let mol = parse_pdb(WATER_PDB).expect("valid PDB");
        let framed = frame_models(&[write_pdb(&mol), write_pdb(&mol)]);

        assert_eq!(framed.matches("MODEL     ").count(), 2, "{framed}");
        assert_eq!(framed.matches("ENDMDL").count(), 2, "{framed}");

        let back = crate::io::reader::read(&framed, crate::io::format::Format::PDB);
        assert_eq!(back.records.len(), 2, "{framed}");
        assert_eq!(
            back.records[0].molecule().unwrap().num_atoms(),
            mol.num_atoms()
        );
    }

    #[test]
    fn test_a_framed_structure_ends_with_endmdl_not_end() {
        // `END` terminates a *file*. Left inside a model it is both wrong and
        // visible: a bare `END` after the last `ENDMDL` reads back as an extra,
        // empty record -- measured, and the artefact OpenBabel's own multi-model
        // output shows.
        let mol = parse_pdb(WATER_PDB).expect("valid PDB");
        assert!(
            write_pdb(&mol).contains("END\n"),
            "the single form keeps it"
        );

        let framed = frame_models(&[write_pdb(&mol), write_pdb(&mol)]);
        let records: Vec<&str> = framed.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(records.last(), Some(&"ENDMDL"), "{framed}");
        assert!(
            !records.iter().any(|l| l.trim_end() == "END"),
            "an END inside a model: {framed}"
        );
    }

    #[test]
    fn test_occupancy_and_b_factor_keep_their_values_not_just_their_flags() {
        // `Carries` is presence-based, so `held` reports B_FACTOR whether the
        // value came back as 42.50 or as a manufactured 0.00 -- which is
        // exactly the defect #173 records against OpenBabel, whose PDB writer
        // zeroes this column while preserving the occupancy beside it. Only a
        // value comparison can tell those apart (#258).
        let mol = parse_pdb(SITES_PDB).expect("valid PDB");
        let expected = [(1.00, 42.50), (0.80, 37.25), (0.55, 55.00)];

        let read_back = |mol: &Molecule| -> Vec<(f64, f64)> {
            (0..mol.num_atoms())
                .map(|i| {
                    let site = mol.site(i).expect("every atom line carries both columns");
                    (
                        site.occupancy.expect("an occupancy"),
                        site.b_factor.expect("a b-factor"),
                    )
                })
                .collect()
        };

        // Twice, following mmcif.rs's #260 regression: a column that survives
        // one pass can still be misaligned enough to move on the next.
        let mut current = mol;
        for pass in 1..=2 {
            let measured = read_back(&current);
            for (i, (occupancy, b_factor)) in expected.iter().enumerate() {
                assert!(
                    (measured[i].0 - occupancy).abs() < 1e-9,
                    "pass {pass}, atom {i}: occupancy {} is not {occupancy}",
                    measured[i].0
                );
                assert!(
                    (measured[i].1 - b_factor).abs() < 1e-9,
                    "pass {pass}, atom {i}: b-factor {} is not {b_factor}",
                    measured[i].1
                );
            }
            current = parse_pdb(&write_pdb(&current)).expect("our own output reads back");
        }
    }

    #[test]
    fn test_an_atom_with_no_site_is_written_as_full_occupancy_and_no_b_factor() {
        // The other half: these columns are not optional, so writing an atom
        // that never had them manufactures 1.00 and 0.00 -- values a reader
        // cannot tell from measured ones. Pinned because #257's `SUPPLIED`
        // table records exactly this, and the two descriptions have to agree.
        let mol = crate::io::smiles::parse_smiles("CC").expect("valid SMILES");
        assert!(mol.sites().is_none(), "no site data to begin with");

        let back = parse_pdb(&write_pdb(&mol)).expect("our own output reads back");
        for i in 0..back.num_atoms() {
            let site = back.site(i).expect("the columns are always written");
            assert_eq!(site.occupancy, Some(1.0));
            assert_eq!(site.b_factor, Some(0.0));
        }
    }

    #[test]
    fn test_element_falls_back_to_column_position_when_the_element_column_is_absent() {
        // " CA " (blank column 13): single-letter element, carbon.
        let ca_carbon = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n";
        let mol = parse_pdb(ca_carbon).expect("valid PDB");
        assert_eq!(mol.atoms()[0].element().symbol(), "C");

        // "CA  " (non-blank column 13): two-letter element, calcium.
        let ca_calcium = "HETATM    1 CA    CA A   1       0.000   0.000   0.000  1.00  0.00\n";
        let mol = parse_pdb(ca_calcium).expect("valid PDB");
        assert_eq!(mol.atoms()[0].element().symbol(), "Ca");
    }

    #[test]
    fn test_a_multi_model_file_reads_as_multiple_structures_when_pre_split() {
        // Frame-splitting on MODEL/ENDMDL is the reader's/supplier's job
        // (see the module doc) -- this only confirms parse_pdb handles an
        // already-isolated single MODEL's body correctly.
        let one_model = WATER_PDB;
        let mol = parse_pdb(one_model).expect("valid PDB");
        assert_eq!(mol.num_atoms(), 3);
    }

    #[test]
    fn test_a_short_line_missing_trailing_columns_does_not_panic() {
        // No element column, no charge column -- common in older files.
        let text = "HETATM    1  O   HOH A   1       0.000   0.000   0.000\n";
        let mol = parse_pdb(text).expect("valid PDB");
        assert_eq!(mol.num_atoms(), 1);
    }

    #[test]
    fn test_an_out_of_range_conect_serial_is_a_clear_error_not_a_panic() {
        let text = format!("{WATER_PDB}CONECT    1   99\n");
        let err = parse_pdb(&text).unwrap_err();
        assert!(matches!(err, PdbError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_all_zero_coordinates_are_not_read_as_a_conformer() {
        // #270: the coordinate columns are mandatory, so a molecule with no
        // conformer is written as zeros.
        let text = "\
ATOM      1  C   UNK     1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  O   UNK     1       0.000   0.000   0.000  1.00  0.00           O
END
";
        let mol = parse_pdb(text).expect("valid PDB");
        assert_eq!(mol.num_atoms(), 2);
        assert!(
            !mol.has_coords3(),
            "all-zero is a placeholder, not a conformer"
        );
        // The site data beside it is real and must survive.
        assert_eq!(mol.site(0).unwrap().occupancy, Some(1.0));
    }

    #[test]
    fn test_garbage_text_reports_no_atoms_instead_of_an_empty_molecule() {
        for input in [
            "",
            "HEADER  SOME PROTEIN\nTITLE   NOTHING HERE\n",
            "not a pdb file at all\n",
        ] {
            assert!(
                matches!(parse_pdb(input), Err(PdbError::NoAtoms)),
                "{input:?} should report NoAtoms"
            );
        }
    }

    #[test]
    fn test_a_bonded_atom_gets_the_count_conect_implies() {
        // The reported symptom (#285). Ethanol through PDB read back as
        // `[C][C][O]` -- the right topology and the wrong molecule, since `[C]`
        // states a carbon with no hydrogens.
        //
        // The answer is what the file says, not what it once meant: RDKit reads
        // our benzene PDB as `C1CCCCC1` too. `CONECT` carries no bond order, so
        // the aromaticity went when the PDB was written and cyclohexane is the
        // honest reading; recovering benzene would need ring perception, which
        // this reader deliberately does not do.
        let molecule = crate::io::smiles::parse_smiles("CCO").expect("valid SMILES");
        let text = write_pdb(&molecule);
        let back = parse_pdb(&text).expect("reads back");

        let counts: Vec<Option<u8>> = (0..back.num_atoms())
            .map(|i| back.atom(i).hydrogens())
            .collect();
        assert_eq!(counts, vec![Some(3), Some(2), Some(1)]);
        assert_eq!(
            crate::io::smiles_writer::write_smiles_for_molecule_canonical(&back),
            "CCO"
        );
    }

    #[test]
    fn test_an_atom_conect_never_mentioned_keeps_an_unknown_count() {
        // The other half, and the reason the unguarded
        // `calculate_implicit_hydrogens` cannot be used here: with no bond
        // orders to subtract it hands an atom its free-atom valence, so every
        // bondless atom would claim a full complement of hydrogens -- a
        // backbone carbon reading as methane.
        //
        // `WATER_PDB` has explicit H atoms and no CONECT at all, so all three
        // atoms are bondless and all three stay unknown.
        let molecule = parse_pdb(WATER_PDB).expect("reads");
        assert_eq!(molecule.num_atoms(), 3);
        assert_eq!(molecule.num_bonds(), 0);
        for i in 0..molecule.num_atoms() {
            assert_eq!(
                molecule.atom(i).hydrogens(),
                None,
                "atom {i} was given a count from no bonds at all"
            );
        }
    }

    #[test]
    fn test_a_partly_connected_file_fills_only_the_atoms_it_described() {
        // A real file rather than a fixture: eight atoms and one CONECT pair,
        // which is what an unguarded fill would turn into a bag of small
        // alkanes. Six atoms stay unknown; the pair becomes a methanol
        // fragment, which is exactly what the file states about them.
        let text = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/corpus/pdb/dipeptide-with-ligand.pdb"
        ))
        .expect("the corpus file is there");
        let molecule = parse_pdb(&text).expect("reads");

        let stated = (0..molecule.num_atoms())
            .filter(|&i| molecule.atom(i).hydrogens().is_some())
            .count();
        assert_eq!(molecule.num_atoms(), 8);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(stated, 2, "only the CONECT pair can imply a count");
        assert_eq!(
            crate::io::smiles_writer::write_smiles_for_molecule_canonical(&molecule),
            "[C].[C].[C].[N].[N].[O].CO"
        );
    }

    #[test]
    fn test_an_explicit_hydrogen_atom_is_not_counted_twice() {
        // Hydrogens in a PDB are ordinary atoms, so the count must come from
        // the bonds they form rather than being added on top: this carbon's
        // four bonds leave it needing none.
        let text = "\
HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  H1  LIG A   1       0.629   0.629   0.629  1.00  0.00           H
HETATM    3  H2  LIG A   1      -0.629  -0.629   0.629  1.00  0.00           H
HETATM    4  H3  LIG A   1      -0.629   0.629  -0.629  1.00  0.00           H
HETATM    5  H4  LIG A   1       0.629  -0.629  -0.629  1.00  0.00           H
CONECT    1    2    3    4    5
END
";
        let molecule = parse_pdb(text).expect("reads");
        assert_eq!(molecule.atom(0).hydrogens(), Some(0));
        assert_eq!(molecule.num_atoms(), 5);
    }

    /// An 80-column record with each `(column, text)` placed at its 1-based
    /// column, so a test states the spec's positions rather than counting
    /// blanks.
    fn record(fields: &[(usize, &str)]) -> String {
        let mut line = vec![b' '; 80];
        for &(col, text) in fields {
            line[col - 1..col - 1 + text.len()].copy_from_slice(text.as_bytes());
        }
        String::from_utf8(line).unwrap().trim_end().to_string()
    }

    fn atom(
        serial: usize,
        het: bool,
        name: &str,
        res: &str,
        seq: &str,
        x: f64,
        el: &str,
    ) -> String {
        let coords = format!("{x:>8.3}{:>8.3}{:>8.3}", 0.0, 0.0);
        record(&[
            (1, if het { "HETATM" } else { "ATOM" }),
            (7, &format!("{serial:>5}")),
            (14, name),
            (18, res),
            (22, "A"),
            (23, &format!("{seq:>4}")),
            (31, &coords),
            (77, &format!("{el:>2}")),
        ])
    }

    #[allow(clippy::too_many_arguments)]
    fn link(
        name1: &str,
        res1: &str,
        seq1: &str,
        name2: &str,
        res2: &str,
        seq2: &str,
        sym2: &str,
    ) -> String {
        record(&[
            (1, "LINK"),
            (14, name1),
            (18, res1),
            (22, "A"),
            (23, &format!("{seq1:>4}")),
            (44, name2),
            (48, res2),
            (52, "A"),
            (53, &format!("{seq2:>4}")),
            (60, "  1555"),
            (67, sym2),
        ])
    }

    fn linked_structure() -> String {
        [
            record(&[
                (1, "SSBOND"),
                (12, "CYS"),
                (16, "A"),
                (18, "   3"),
                (26, "CYS"),
                (30, "A"),
                (32, "  40"),
                (60, "  1555"),
                (67, "  1555"),
            ]),
            link("NE2", "HIS", "50", "ZN", "ZN", "101", "  1555"),
            link("ND2", "ASN", "60", "C1", "NAG", "102", "  1555"),
            link("O4", "NAG", "102", "C1", "NAG", "103", "  2555"),
            atom(1, false, "SG", "CYS", "3", 0.0, "S"),
            atom(2, false, "SG", "CYS", "40", 2.0, "S"),
            atom(3, false, "NE2", "HIS", "50", 5.0, "N"),
            atom(4, false, "ND2", "ASN", "60", 9.0, "N"),
            atom(5, true, "ZN", "ZN", "101", 7.0, "ZN"),
            atom(6, true, "C1", "NAG", "102", 10.4, "C"),
            atom(7, true, "O4", "NAG", "102", 12.0, "O"),
            atom(8, true, "C1", "NAG", "103", 13.0, "C"),
            "END".to_string(),
        ]
        .join("\n")
    }

    #[test]
    fn test_links_bond_only_covalent_partners_inside_the_model() {
        let options = PdbReadOptions {
            residue_templates: true,
            ..PdbReadOptions::default()
        };
        let mol = parse_pdb_with_options(&linked_structure(), &options).expect("reads");
        let bonded = |a: usize, b: usize| mol.graph().has_edge(a, b);
        assert!(bonded(0, 1), "SSBOND");
        assert!(bonded(3, 5), "LINK between non-metals");
        assert!(!bonded(2, 4), "a LINK to a metal is coordination");
        assert!(!bonded(6, 7), "a LINK across a symmetry operator");
        assert_eq!(mol.num_bonds(), 2);
    }

    #[test]
    fn test_links_and_templates_are_off_by_default() {
        let text = linked_structure();
        let mol = parse_pdb(&text).expect("reads");
        assert_eq!(mol.num_bonds(), 0);
        let with_default = parse_pdb_with_options(&text, &PdbReadOptions::default()).unwrap();
        assert_eq!(write_pdb(&mol), write_pdb(&with_default));
    }

    #[test]
    fn test_perception_keeps_a_conect_bond_once_and_fills_hydrogens_from_both() {
        // O-H1 is stated; O-H2 is only perceived. A stated 2.9 Å O...O contact
        // is far past any covalent cutoff and must survive as stated.
        let text = "\
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  H1  HOH A   1       0.759   0.000   0.504  1.00  0.00           H
HETATM    3  H2  HOH A   1       0.759   0.000  -0.504  1.00  0.00           H
HETATM    4  O   HOH A   2       2.900   0.000   0.000  1.00  0.00           O
CONECT    1    2    4
END
";
        let options = PdbReadOptions {
            perceive_bonds: true,
            ..Default::default()
        };
        let mol = parse_pdb_with_options(text, &options).unwrap();
        let mut pairs: Vec<_> = mol.bonds().iter().map(|b| b.atoms()).collect();
        pairs.sort_unstable();
        assert_eq!(pairs, vec![(0, 1), (0, 2), (0, 3)]);
        assert_eq!(mol.atom(0).hydrogens(), Some(0));

        assert_eq!(parse_pdb(text).unwrap().num_bonds(), 2, "off by default");
    }
}
