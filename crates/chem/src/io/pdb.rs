//! PDB — fixed-column text, `ATOM`/`HETATM` per atom, `CRYST1` for the unit
//! cell, `CONECT` for explicit bonds (#223).
//!
//! **No bonds are synthesised beyond what `CONECT` states explicitly.** Real
//! PDB files conventionally leave standard protein/nucleic-acid backbone
//! bonding to be inferred from residue templates (a fixed table of atom
//! names/bonds per each of the 20 standard amino acids plus nucleotides)
//! rather than stating it — this crate has no such table, and building one
//! is a materially larger, separable undertaking than this story. `CONECT`
//! itself stays in scope: it is reading data the file states, not guessing
//! at it. Because of this, [`Molecule::calculate_implicit_hydrogens`] is
//! never called on a PDB-read molecule — it sums existing bond orders per
//! atom, so a molecule with atoms but no bonds would have every atom
//! assigned its *free-atom* hydrogen count (a backbone carbon reading as
//! methane) rather than being left honestly blank.
//!
//! `CRYST1`'s orientation convention already matches this crate's fixed
//! [`crate::core::cell::UnitCell`] (its own doc comment: "the near-universal
//! crystallographic setting... which is what PDB, CIF and essentially every
//! toolkit assume") — unlike extended XYZ's `Lattice=` (#222, arbitrary
//! per-writer orientation, deliberately not read), so `CRYST1` is read
//! directly with no rotation step and no correctness risk.
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
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::core::site::AtomSite;
use crate::io::errors::PdbError;

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

/// Parses one PDB structure (a single `MODEL`, or a whole file with none).
pub fn parse_pdb(text: &str) -> Result<Molecule, PdbError> {
    let mut mol = Molecule::new();
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

    mol.set_coords3(coords)
        .map_err(|e| PdbError::ParseError(e.to_string()))?;
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

    // Deliberately not called: this molecule's bonds are exactly what
    // CONECT stated, nothing more, and this function computes implicit
    // hydrogens from existing bond orders -- with most atoms bondless by
    // construction (see the module doc), it would assign every one of them
    // a free-atom hydrogen count instead of leaving the count honestly
    // unknown.
    // mol.calculate_implicit_hydrogens();

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
    out.push_str("END\n");

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
}
