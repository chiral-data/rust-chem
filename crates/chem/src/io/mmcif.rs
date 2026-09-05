//! mmCIF — tag-value text (`_category.item value`), with `loop_` blocks for
//! repeated records. The modern PDB successor (#224).
//!
//! Data-model-wise this is #223 (PDB) again: [`AtomSite`], [`Chain`]/
//! [`Residue`] (already carrying mmCIF's own dual `auth_*`/`label_*`
//! numbering — built for exactly this), [`UnitCell`]/[`SpaceGroup`] all
//! already exist and need no changes. What's new is the *grammar*, not the
//! model.
//!
//! **No bonds at all** — stricter than PDB's `CONECT`-only cut. mmCIF's
//! rough analogue, `_struct_conn`, is a separate loop keyed by
//! `label_asym_id`/`label_seq_id`/`label_atom_id` triples rather than a
//! simple serial number, and rarer in ordinary deposited structures than
//! `CONECT`. Deferred entirely, stated here rather than silently dropped.
//!
//! **Grammar subset**: bare (whitespace-delimited) tokens, single- and
//! double-quoted tokens. Multi-line semicolon (`;...;`) text fields are
//! **not** parsed — real files never need one inside `_atom_site`/`_cell`/
//! `_symmetry` data, only for free-text descriptions this module doesn't
//! read, so the limitation costs nothing in practice.
//!
//! **Multiple `data_` blocks** in one file each become a separate record —
//! splitting is the reader's/supplier's job
//! ([`crate::io::reader::read_mmcif_with_options`],
//! [`crate::io::supplier::MmcifSupplier`]), the same division every prior
//! format has. **A varying `pdbx_PDB_model_num` within one block's
//! `_atom_site` loop** (mmCIF's actual NMR-ensemble convention — unlike
//! PDB, models share one loop rather than living in separate blocks) is
//! *not* split into multiple records: only the first model's atoms are
//! read, later ones skipped rather than merged in. Stated plainly rather
//! than silently combining distinct models into one incoherent structure;
//! most deposited structures are single-model, and NMR ensembles
//! specifically are the case this doesn't handle.

use crate::core::atom::{Atom, Element};
use crate::core::cell::{SpaceGroup, UnitCell};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::core::site::AtomSite;
use crate::io::errors::MmcifError;

/// Splits one line into its whitespace-delimited or quoted tokens. A quote
/// is a delimiter only where it starts a token (preceded by whitespace or
/// nothing) and ends one (followed by whitespace or nothing) — the
/// standard CIF quoting rule, so `O5'` (a common atom name containing an
/// apostrophe) is not mistaken for the start of a quoted string.
fn tokenize_line(line: &str) -> Vec<String> {
    let chars: Vec<char> = line.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        while i < chars.len() && chars[i].is_whitespace() {
            i += 1;
        }
        if i >= chars.len() {
            break;
        }
        if chars[i] == '\'' || chars[i] == '"' {
            let quote = chars[i];
            i += 1;
            let start = i;
            while i < chars.len()
                && !(chars[i] == quote && chars.get(i + 1).is_none_or(|c| c.is_whitespace()))
            {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
            i += 1; // the closing quote
        } else {
            let start = i;
            while i < chars.len() && !chars[i].is_whitespace() {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
        }
    }
    tokens
}

/// `.` and `?` both mean "not provided" for an optional mmCIF value.
fn cif_value(row: &[String], idx: Option<usize>) -> Option<&str> {
    let s = row.get(idx?)?.as_str();
    (s != "." && s != "?").then_some(s)
}

fn find_tag(tags: &[String], full: &str) -> Option<usize> {
    tags.iter().position(|t| t == full)
}

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

/// Parses one mmCIF `data_` block (or a whole file with only one).
pub fn parse_mmcif(text: &str) -> Result<Molecule, MmcifError> {
    let mut mol = Molecule::new();
    let mut singles: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut sites = Vec::new();
    let mut coords = Vec::new();
    let mut keys = Vec::new();

    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        let trimmed = lines[i].trim();
        let lower = trimmed.to_ascii_lowercase();
        if trimmed.is_empty() || trimmed.starts_with('#') || lower.starts_with("data_") {
            i += 1;
            continue;
        }

        if lower == "loop_" {
            i += 1;
            let mut tags = Vec::new();
            while i < lines.len() && lines[i].trim_start().starts_with('_') {
                tags.push(lines[i].trim().to_string());
                i += 1;
            }
            let mut rows: Vec<Vec<String>> = Vec::new();
            while i < lines.len() {
                let row_trimmed = lines[i].trim();
                let row_lower = row_trimmed.to_ascii_lowercase();
                if row_trimmed.is_empty()
                    || row_trimmed.starts_with('#')
                    || row_trimmed.starts_with('_')
                    || row_lower == "loop_"
                    || row_lower.starts_with("data_")
                {
                    break;
                }
                rows.push(tokenize_line(lines[i]));
                i += 1;
            }

            if tags.first().is_some_and(|t| t.starts_with("_atom_site.")) {
                read_atom_site_loop(&tags, &rows, &mut mol, &mut sites, &mut coords, &mut keys)?;
            }
            // Every other loop (e.g. `_struct_conn`) is out of scope --
            // see the module doc.
            continue;
        }

        if trimmed.starts_with('_') {
            let tokens = tokenize_line(lines[i]);
            if tokens.len() >= 2 {
                singles.insert(tokens[0].clone(), tokens[1].clone());
            }
        }
        i += 1;
    }

    mol.set_coords3(coords)
        .map_err(|e| MmcifError::ParseError(e.to_string()))?;
    mol.set_sites(sites)
        .map_err(|e| MmcifError::ParseError(e.to_string()))?;

    let (chains, residues) = group_into_chains_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| MmcifError::ParseError(e.to_string()))?;

    if let (Some(a), Some(b), Some(c), Some(alpha), Some(beta), Some(gamma)) = (
        singles.get("_cell.length_a").and_then(|s| s.parse().ok()),
        singles.get("_cell.length_b").and_then(|s| s.parse().ok()),
        singles.get("_cell.length_c").and_then(|s| s.parse().ok()),
        singles
            .get("_cell.angle_alpha")
            .and_then(|s| s.parse().ok()),
        singles.get("_cell.angle_beta").and_then(|s| s.parse().ok()),
        singles
            .get("_cell.angle_gamma")
            .and_then(|s| s.parse().ok()),
    ) {
        mol.set_cell(UnitCell::new(a, b, c, alpha, beta, gamma))
            .map_err(|e| MmcifError::ParseError(e.to_string()))?;
    }
    if let Some(symbol) = singles.get("_symmetry.space_group_name_H-M") {
        mol.set_space_group(SpaceGroup::from_symbol(symbol.as_str()));
    }

    Ok(mol)
}

#[allow(clippy::too_many_arguments)]
fn read_atom_site_loop(
    tags: &[String],
    rows: &[Vec<String>],
    mol: &mut Molecule,
    sites: &mut Vec<AtomSite>,
    coords: &mut Vec<Point3>,
    keys: &mut Vec<ResidueKey>,
) -> Result<(), MmcifError> {
    let tag = |suffix: &str| find_tag(tags, &format!("_atom_site.{suffix}"));

    let group_pdb = tag("group_PDB");
    let type_symbol = tag("type_symbol").ok_or_else(|| {
        MmcifError::ParseError("_atom_site loop has no type_symbol column".to_string())
    })?;
    let label_atom_id = tag("label_atom_id");
    let label_alt_id = tag("label_alt_id");
    let label_comp_id = tag("label_comp_id");
    let label_asym_id = tag("label_asym_id");
    let auth_seq_id = tag("auth_seq_id").or(tag("label_seq_id"));
    let auth_asym_id = tag("auth_asym_id").or(label_asym_id);
    let auth_comp_id = tag("auth_comp_id").or(label_comp_id);
    let ins_code = tag("pdbx_PDB_ins_code");
    let cartn_x = tag("Cartn_x").ok_or_else(|| {
        MmcifError::ParseError("_atom_site loop has no Cartn_x column".to_string())
    })?;
    let cartn_y = tag("Cartn_y").ok_or_else(|| {
        MmcifError::ParseError("_atom_site loop has no Cartn_y column".to_string())
    })?;
    let cartn_z = tag("Cartn_z").ok_or_else(|| {
        MmcifError::ParseError("_atom_site loop has no Cartn_z column".to_string())
    })?;
    let occupancy = tag("occupancy");
    let b_iso = tag("B_iso_or_equiv");
    let model_num = tag("pdbx_PDB_model_num");

    // A varying pdbx_PDB_model_num is mmCIF's NMR-ensemble convention --
    // only the first model's rows are read, see the module doc.
    let first_model = model_num.and_then(|idx| rows.first().and_then(|r| r.get(idx).cloned()));

    for row in rows {
        if let (Some(idx), Some(first)) = (model_num, &first_model)
            && row.get(idx).is_some_and(|m| m != first)
        {
            continue;
        }

        let is_hetero = group_pdb
            .and_then(|idx| row.get(idx))
            .is_some_and(|g| g.eq_ignore_ascii_case("HETATM"));
        let symbol = row.get(type_symbol).ok_or_else(|| {
            MmcifError::InvalidAtomRow("row has no type_symbol value".to_string())
        })?;
        let element = element_from_symbol(symbol)
            .ok_or_else(|| MmcifError::InvalidElement(symbol.to_string()))?;

        let x: f64 = row
            .get(cartn_x)
            .ok_or_else(|| MmcifError::InvalidAtomRow(row.join(" ")))?
            .parse()
            .map_err(|_| MmcifError::InvalidAtomRow(row.join(" ")))?;
        let y: f64 = row
            .get(cartn_y)
            .ok_or_else(|| MmcifError::InvalidAtomRow(row.join(" ")))?
            .parse()
            .map_err(|_| MmcifError::InvalidAtomRow(row.join(" ")))?;
        let z: f64 = row
            .get(cartn_z)
            .ok_or_else(|| MmcifError::InvalidAtomRow(row.join(" ")))?
            .parse()
            .map_err(|_| MmcifError::InvalidAtomRow(row.join(" ")))?;

        let atom_idx = mol.add_atom(Atom::new(element));
        coords.push(Point3::new(x, y, z));
        sites.push(AtomSite {
            name: cif_value(row, label_atom_id).map(str::to_string),
            alt_loc: cif_value(row, label_alt_id).and_then(|s| s.chars().next()),
            partial_charge: None,
            occupancy: cif_value(row, occupancy).and_then(|s| s.parse().ok()),
            b_factor: cif_value(row, b_iso).and_then(|s| s.parse().ok()),
            radius: None,
        });
        keys.push(ResidueKey {
            chain_id: cif_value(row, auth_asym_id).unwrap_or("").to_string(),
            name: cif_value(row, auth_comp_id).unwrap_or("UNK").to_string(),
            sequence: cif_value(row, auth_seq_id)
                .and_then(|s| s.parse().ok())
                .unwrap_or(1),
            insertion_code: cif_value(row, ins_code).and_then(|s| s.chars().next()),
            is_hetero,
        });
        let _ = atom_idx;
    }

    Ok(())
}

/// Writes one mmCIF `data_` block.
pub fn write_mmcif(mol: &Molecule) -> String {
    let mut out = String::from("data_chem\n");

    if let Some(cell) = mol.cell() {
        out.push_str(&format!(
            "_cell.length_a {:.3}\n_cell.length_b {:.3}\n_cell.length_c {:.3}\n\
             _cell.angle_alpha {:.2}\n_cell.angle_beta {:.2}\n_cell.angle_gamma {:.2}\n",
            cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma
        ));
    }
    if let Some(sg) = mol.space_group().and_then(|g| g.symbol.as_deref()) {
        out.push_str(&format!("_symmetry.space_group_name_H-M '{sg}'\n"));
    }

    out.push_str("loop_\n");
    for tag in [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "auth_seq_id",
        "auth_asym_id",
        "auth_comp_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_PDB_model_num",
    ] {
        out.push_str(&format!("_atom_site.{tag}\n"));
    }

    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let residue = mol.residue_of(i);
        let chain = mol.chain_of(i);
        let p = mol.coord3(i).unwrap_or(Point3::new(0.0, 0.0, 0.0));

        let is_hetero = residue.map(|r| r.is_hetero).unwrap_or(false);
        let group = if is_hetero { "HETATM" } else { "ATOM" };
        let res_name = residue.map(|r| r.name.as_str()).unwrap_or("UNK");
        let res_seq = residue.map(|r| r.sequence).unwrap_or(1);
        let icode = residue
            .and_then(|r| r.insertion_code)
            .map(String::from)
            .unwrap_or_else(|| "?".to_string());
        let chain_id = chain.map(|c| c.id.as_str()).unwrap_or(".");
        let alt_loc = site
            .and_then(|s| s.alt_loc)
            .map(String::from)
            .unwrap_or_else(|| ".".to_string());
        let name = site
            .and_then(|s| s.name.as_deref())
            .unwrap_or(atom.element().symbol());
        let occupancy = site.and_then(|s| s.occupancy).unwrap_or(1.0);
        let b_factor = site.and_then(|s| s.b_factor).unwrap_or(0.0);

        out.push_str(&format!(
            "{group} {serial} {element} {name} {alt_loc} {res_name} {chain_id} \
             {res_seq} {chain_id} {res_name} {icode} {x:.3} {y:.3} {z:.3} \
             {occupancy:.2} {b_factor:.2} 1\n",
            serial = i + 1,
            element = atom.element().symbol(),
            x = p.x,
            y = p.y,
            z = p.z,
        ));
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_MMCIF: &str = "\
data_water
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.auth_comp_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
HETATM 1 O O . HOH A 1 A HOH . 0.000 0.000 0.000 1.00 20.00 1
HETATM 2 H H1 . HOH A 1 A HOH . 0.759 0.000 0.504 1.00 20.00 1
HETATM 3 H H2 . HOH A 1 A HOH . 0.759 0.000 -0.504 1.00 20.00 1
";

    #[test]
    fn test_a_single_block_round_trips() {
        let mol = parse_mmcif(WATER_MMCIF).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.residues().len(), 1);
        assert!(mol.residues()[0].is_hetero);
        assert_eq!(mol.residues()[0].name, "HOH");
        assert_eq!(mol.chains()[0].id, "A");
        assert_eq!(mol.site(0).unwrap().occupancy, Some(1.0));

        let written = write_mmcif(&mol);
        let back = parse_mmcif(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.residues()[0].name, "HOH");
        assert_eq!(back.chains()[0].id, "A");
    }

    #[test]
    fn test_cell_and_space_group_round_trip() {
        let text = format!(
            "data_x\n_cell.length_a 10.000\n_cell.length_b 20.000\n_cell.length_c 30.000\n\
             _cell.angle_alpha 90.00\n_cell.angle_beta 90.00\n_cell.angle_gamma 90.00\n\
             _symmetry.space_group_name_H-M 'P 21 21 21'\n{}",
            &WATER_MMCIF["data_water\n".len()..]
        );
        let mol = parse_mmcif(&text).expect("valid mmCIF");
        let cell = mol.cell().expect("has a cell");
        assert_eq!((cell.a, cell.b, cell.c), (10.0, 20.0, 30.0));
        assert_eq!(
            mol.space_group().and_then(|g| g.symbol.as_deref()),
            Some("P 21 21 21")
        );

        let written = write_mmcif(&mol);
        let back = parse_mmcif(&written).expect("round trips");
        assert_eq!(back.cell(), mol.cell());
    }

    #[test]
    fn test_a_quoted_value_containing_a_space_parses_correctly() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_comp_id\n\
                    _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
                    C 'LIG AND' 0.0 0.0 0.0\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.residues()[0].name, "LIG AND");
    }

    #[test]
    fn test_question_mark_and_dot_both_mean_absent() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n\
                    _atom_site.B_iso_or_equiv\nC 0.0 0.0 0.0 ? .\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.site(0).unwrap().occupancy, None);
        assert_eq!(mol.site(0).unwrap().b_factor, None);
    }

    #[test]
    fn test_multiple_data_blocks_are_this_modules_callers_job() {
        // parse_mmcif itself only ever sees one already-isolated block --
        // confirm it does not choke on a bare `data_` line at the top,
        // which is all that distinguishes "isolated" from "not".
        let mol = parse_mmcif(WATER_MMCIF).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 3);
    }

    #[test]
    fn test_only_the_first_model_is_read_when_the_model_number_varies() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n\
                    C 0.0 0.0 0.0 1\nC 1.0 1.0 1.0 2\n";
        let mol = parse_mmcif(text).expect("valid mmCIF");
        assert_eq!(mol.num_atoms(), 1);
        assert_eq!(mol.coord3(0), Some(Point3::new(0.0, 0.0, 0.0)));
    }

    #[test]
    fn test_a_malformed_row_is_a_clear_error_not_a_panic() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\nC not-a-number 0.0 0.0\n";
        let err = parse_mmcif(text).unwrap_err();
        assert!(matches!(err, MmcifError::InvalidAtomRow(_)), "{err}");
    }

    #[test]
    fn test_an_unrecognised_element_is_a_clear_error() {
        let text = "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n\
                    _atom_site.Cartn_y\n_atom_site.Cartn_z\nXx 0.0 0.0 0.0\n";
        let err = parse_mmcif(text).unwrap_err();
        assert!(matches!(err, MmcifError::InvalidElement(_)), "{err}");
    }
}
