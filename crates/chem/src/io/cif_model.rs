//! The CIF data model, shared between mmCIF text ([`crate::io::mmcif`]) and
//! BinaryCIF ([`crate::io::bcif`], #319).
//!
//! Both formats describe the same thing: a `singles` map of one-off
//! `_category.item value` pairs (`_cell.length_a`, `_symmetry.space_group_
//! name_H-M`, ...) and, for the one loop this crate reads, an `_atom_site`
//! table as parallel `tags`/`rows`. mmCIF text tokenizes lines into that
//! shape; BinaryCIF decodes MessagePack categories/columns into the same
//! shape (see `bcif::decode_bcif`). Everything below this line — turning
//! that shape into a [`Molecule`] and back — used to live inside
//! `mmcif.rs`'s text-specific parser and writer; extracted here so a second
//! format can reuse it instead of re-deriving atom-site/chain/cell handling
//! from scratch.

use std::collections::HashMap;

use crate::core::atom::{Atom, Element};
use crate::core::cell::{SpaceGroup, UnitCell};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::{Point3, is_placeholder_3d};
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::core::site::AtomSite;
use crate::io::errors::MmcifError;

/// An `_atom_site` table: one CIF tag per column, one row per atom -- the
/// shape both mmCIF text and BinaryCIF decode into before
/// [`build_molecule`] ever runs.
pub(crate) type AtomSiteRows = (Vec<String>, Vec<Vec<String>>);

/// `.` and `?` both mean "not provided" for an optional CIF value.
pub(crate) fn cif_value(row: &[String], idx: Option<usize>) -> Option<&str> {
    let s = row.get(idx?)?.as_str();
    (s != "." && s != "?").then_some(s)
}

pub(crate) fn find_tag(tags: &[String], full: &str) -> Option<usize> {
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
    // only the first model's rows are read, see `mmcif`'s module doc.
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

        mol.add_atom(Atom::new(element));
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
    }

    Ok(())
}

/// Builds a [`Molecule`] from a CIF-shaped `singles` map and an optional
/// `_atom_site` table, independent of which format the shape came from.
///
/// `atom_site` being `Some` (even with zero rows) means the category/loop
/// was present at all — mirrors mmCIF text's own `saw_atom_site_loop`
/// distinction between "a legitimately empty structure" and "no atoms were
/// ever declared" (#268).
pub(crate) fn build_molecule(
    singles: &HashMap<String, String>,
    atom_site: Option<(&[String], &[Vec<String>])>,
) -> Result<Molecule, MmcifError> {
    let mut mol = Molecule::new();
    let mut sites = Vec::new();
    let mut coords = Vec::new();
    let mut keys = Vec::new();

    if let Some((tags, rows)) = atom_site {
        read_atom_site_loop(tags, rows, &mut mol, &mut sites, &mut coords, &mut keys)?;
    }

    // Zeros in these columns are what a format with no room to say "unknown"
    // writes for a molecule that has no conformer, so believing them back is
    // how a converted molecule ended up undrawable (#270).
    if !is_placeholder_3d(&coords) {
        mol.set_coords3(coords)
            .map_err(|e| MmcifError::ParseError(e.to_string()))?;
    }
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

    if mol.num_atoms() == 0 && atom_site.is_none() {
        return Err(MmcifError::NoAtoms);
    }

    Ok(mol)
}

/// A CIF value, or `fallback` when there is nothing to write.
///
/// CIF has no empty token: a value is a word, `.` (inapplicable) or `?`
/// (unknown). An empty string reaches the file as *nothing*, so the row ends
/// up with fewer values than its `loop_` header declares and every column
/// after it shifts -- which is how a second mmCIF write used to move the
/// B-factor into the model-number column (#260).
fn cif_token<'a>(value: Option<&'a str>, fallback: &'a str) -> &'a str {
    match value {
        Some(text) if !text.is_empty() => text,
        _ => fallback,
    }
}

/// How a numeric value in [`build_rows`] is turned into the CIF string
/// bridge — mmCIF text and BinaryCIF must not share one.
///
/// mmCIF text has always truncated coordinates/cell lengths to three decimal
/// places and angles/occupancy/B-factor to two — that's `Text`. BinaryCIF's
/// own encoding chain is exact, so nothing upstream of it should truncate
/// first (#319); `Full` uses `f64`'s `Display`, which round-trips exactly
/// (the shortest decimal that reads back to the same bits).
pub(crate) enum RowPrecision {
    Text,
    Full,
}

impl RowPrecision {
    fn length(&self, v: f64) -> String {
        match self {
            RowPrecision::Text => format!("{v:.3}"),
            RowPrecision::Full => v.to_string(),
        }
    }

    fn angle_or_score(&self, v: f64) -> String {
        match self {
            RowPrecision::Text => format!("{v:.2}"),
            RowPrecision::Full => v.to_string(),
        }
    }
}

/// The `_atom_site` tags this crate writes, in column order — the inverse of
/// [`build_molecule`]'s `read_atom_site_loop`, and the same set #223 (PDB)
/// established mmCIF should carry.
const ATOM_SITE_COLUMNS: &[&str] = &[
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
];

/// Builds the CIF-shaped `singles`/`_atom_site` rows for `mol` — the inverse
/// of [`build_molecule`]. Infallible, like the text writer this replaces:
/// every `Molecule` this crate can hold has a row for every atom.
pub(crate) fn build_rows(
    mol: &Molecule,
    precision: RowPrecision,
) -> (HashMap<String, String>, AtomSiteRows) {
    let mut singles = HashMap::new();
    if let Some(cell) = mol.cell() {
        singles.insert("_cell.length_a".to_string(), precision.length(cell.a));
        singles.insert("_cell.length_b".to_string(), precision.length(cell.b));
        singles.insert("_cell.length_c".to_string(), precision.length(cell.c));
        singles.insert(
            "_cell.angle_alpha".to_string(),
            precision.angle_or_score(cell.alpha),
        );
        singles.insert(
            "_cell.angle_beta".to_string(),
            precision.angle_or_score(cell.beta),
        );
        singles.insert(
            "_cell.angle_gamma".to_string(),
            precision.angle_or_score(cell.gamma),
        );
    }
    if let Some(sg) = mol.space_group().and_then(|g| g.symbol.as_deref()) {
        singles.insert("_symmetry.space_group_name_H-M".to_string(), sg.to_string());
    }

    let tags: Vec<String> = ATOM_SITE_COLUMNS
        .iter()
        .map(|t| format!("_atom_site.{t}"))
        .collect();

    let mut rows = Vec::with_capacity(mol.num_atoms());
    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let residue = mol.residue_of(i);
        let chain = mol.chain_of(i);
        let p = mol.coord3(i).unwrap_or(Point3::new(0.0, 0.0, 0.0));

        let is_hetero = residue.map(|r| r.is_hetero).unwrap_or(false);
        let group = if is_hetero { "HETATM" } else { "ATOM" };
        let res_name = cif_token(residue.map(|r| r.name.as_str()), "UNK");
        let res_seq = residue.map(|r| r.sequence).unwrap_or(1);
        let icode = residue
            .and_then(|r| r.insertion_code)
            .map(String::from)
            .unwrap_or_else(|| "?".to_string());
        let chain_id = cif_token(chain.map(|c| c.id.as_str()), ".");
        let alt_loc = site
            .and_then(|s| s.alt_loc)
            .map(String::from)
            .unwrap_or_else(|| ".".to_string());
        let name = site
            .and_then(|s| s.name.as_deref())
            .unwrap_or(atom.element().symbol());
        let occupancy = site.and_then(|s| s.occupancy).unwrap_or(1.0);
        let b_factor = site.and_then(|s| s.b_factor).unwrap_or(0.0);

        rows.push(vec![
            group.to_string(),
            (i + 1).to_string(),
            atom.element().symbol().to_string(),
            name.to_string(),
            alt_loc,
            res_name.to_string(),
            chain_id.to_string(),
            res_seq.to_string(),
            chain_id.to_string(),
            res_name.to_string(),
            icode,
            precision.length(p.x),
            precision.length(p.y),
            precision.length(p.z),
            precision.angle_or_score(occupancy),
            precision.angle_or_score(b_factor),
            "1".to_string(),
        ]);
    }

    (singles, (tags, rows))
}
