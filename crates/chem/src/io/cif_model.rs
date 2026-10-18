//! The CIF data model, shared between mmCIF text ([`crate::io::mmcif`]) and
//! BinaryCIF ([`crate::io::bcif`], #319).
//!
//! Both formats describe the same thing: a `singles` map of one-off
//! `_category.item value` pairs (`_cell.length_a`, `_symmetry.space_group_
//! name_H-M`, ...), an `_atom_site` table as parallel `tags`/`rows`, and
//! every other loop in the block, by category, in the same shape. mmCIF text tokenizes lines into that
//! shape; BinaryCIF decodes MessagePack categories/columns into the same
//! shape (see `bcif::decode_bcif`). Everything below this line — turning
//! that shape into a [`Molecule`] and back — used to live inside
//! `mmcif.rs`'s text-specific parser and writer; extracted here so a second
//! format can reuse it instead of re-deriving atom-site/chain/cell handling
//! from scratch.

use std::collections::HashMap;

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::cell::{SpaceGroup, UnitCell};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::{Point3, is_placeholder_3d};
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::core::site::AtomSite;
use crate::io::errors::MmcifError;
use crate::io::options::MmcifReadOptions;

/// An `_atom_site` table: one CIF tag per column, one row per atom -- the
/// shape both mmCIF text and BinaryCIF decode into before
/// [`build_molecule`] ever runs.
pub(crate) type AtomSiteRows = (Vec<String>, Vec<Vec<String>>);

/// Every other loop in the block (`_struct_conn`, `_chem_comp_bond`, ...),
/// keyed by its lower-cased category name, in the same tags/rows shape.
pub(crate) type CifLoops = HashMap<String, AtomSiteRows>;

/// The category a CIF tag belongs to, lower-cased: `_struct_conn.id` →
/// `_struct_conn`.
pub(crate) fn category_of(tag: &str) -> String {
    tag.split('.').next().unwrap_or(tag).to_ascii_lowercase()
}

/// `.` and `?` both mean "not provided" for an optional CIF value.
pub(crate) fn cif_value(row: &[String], idx: Option<usize>) -> Option<&str> {
    let s = row.get(idx?)?.as_str();
    (s != "." && s != "?").then_some(s)
}

pub(crate) fn find_tag(tags: &[String], full: &str) -> Option<usize> {
    tags.iter().position(|t| t == full)
}

/// `pub(crate)`: shared with [`crate::io::cif_core`] (#320), which parses
/// the same `_atom_site_type_symbol`-shaped element symbols.
pub(crate) fn element_from_symbol(sym: &str) -> Option<Element> {
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

/// `pub(crate)`: shared with [`crate::io::psf`] (#321), which groups atoms
/// into chains/residues by its own segid/resid/resname triple rather than
/// mmCIF's `auth_*` fields, but through this same grouping algorithm.
#[derive(Clone, PartialEq, Eq)]
pub(crate) struct ResidueKey {
    pub(crate) chain_id: String,
    pub(crate) name: String,
    pub(crate) sequence: i32,
    pub(crate) insertion_code: Option<char>,
    pub(crate) is_hetero: bool,
}

pub(crate) fn group_into_chains_and_residues(keys: &[ResidueKey]) -> (Vec<Chain>, Vec<Residue>) {
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

/// One later model of an NMR ensemble (#404): its `pdbx_PDB_model_num`, and
/// its positions in the first model's atom order, or why it was rejected.
pub(crate) struct LaterModel {
    pub(crate) number: String,
    pub(crate) positions: Result<Vec<Point3>, String>,
}

/// A block's first model as a [`Molecule`], plus every later model.
pub(crate) struct CifModels {
    pub(crate) molecule: Molecule,
    pub(crate) later: Vec<LaterModel>,
}

/// [`build_molecule`] over the first `pdbx_PDB_model_num` only, with every
/// later model checked against it (#404). Topology and bonds come from the
/// first model alone, so a later model contributes nothing but positions.
pub(crate) fn build_models(
    singles: &HashMap<String, String>,
    atom_site: Option<(&[String], &[Vec<String>])>,
    loops: &CifLoops,
    options: &MmcifReadOptions,
) -> Result<CifModels, MmcifError> {
    let Some((tags, rows)) = atom_site else {
        let molecule = build_molecule(singles, None, loops, options)?;
        return Ok(CifModels {
            molecule,
            later: Vec::new(),
        });
    };

    let groups = group_rows_by_model(tags, rows);
    if groups.len() <= 1 {
        let molecule = build_molecule(singles, Some((tags, rows)), loops, options)?;
        return Ok(CifModels {
            molecule,
            later: Vec::new(),
        });
    }

    let first: Vec<Vec<String>> = groups[0].1.iter().map(|&i| rows[i].clone()).collect();
    let molecule = build_molecule(singles, Some((tags, &first)), loops, options)?;
    let later = groups[1..]
        .iter()
        .map(|(number, indices)| {
            let model: Vec<&Vec<String>> = indices.iter().map(|&i| &rows[i]).collect();
            LaterModel {
                number: number.clone(),
                positions: later_model_positions(tags, &first, &model, &groups[0].0),
            }
        })
        .collect();
    Ok(CifModels { molecule, later })
}

/// Row indices per model number, models in order of first appearance; one
/// group when the column is absent.
fn group_rows_by_model(tags: &[String], rows: &[Vec<String>]) -> Vec<(String, Vec<usize>)> {
    let Some(col) = find_tag(tags, "_atom_site.pdbx_PDB_model_num") else {
        return vec![(String::new(), (0..rows.len()).collect())];
    };
    let mut groups: Vec<(String, Vec<usize>)> = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let number = row.get(col).map(String::as_str).unwrap_or("");
        match groups.iter_mut().find(|(n, _)| n == number) {
            Some((_, indices)) => indices.push(i),
            None => groups.push((number.to_string(), vec![i])),
        }
    }
    groups
}

/// A later model's positions, provided every atom matches the first model's
/// on element, atom name, residue name, chain, sequence number and alt_loc.
/// Stricter than PDB's count-only check because mmCIF rows state identity.
fn later_model_positions(
    tags: &[String],
    first: &[Vec<String>],
    model: &[&Vec<String>],
    first_number: &str,
) -> Result<Vec<Point3>, String> {
    if model.len() != first.len() {
        return Err(format!(
            "{} atoms where model {first_number} has {}",
            model.len(),
            first.len()
        ));
    }

    let tag = |suffix: &str| find_tag(tags, &format!("_atom_site.{suffix}"));
    let identity_cols = [
        tag("type_symbol"),
        tag("label_atom_id"),
        tag("auth_comp_id").or(tag("label_comp_id")),
        tag("auth_asym_id").or(tag("label_asym_id")),
        tag("auth_seq_id").or(tag("label_seq_id")),
        tag("label_alt_id"),
    ];
    let identity =
        |row: &[String]| identity_cols.map(|idx| cif_value(row, idx).unwrap_or("").to_string());
    let xyz = [tag("Cartn_x"), tag("Cartn_y"), tag("Cartn_z")];

    let mut positions = Vec::with_capacity(model.len());
    for (i, (expected, row)) in first.iter().zip(model).enumerate() {
        let (want, got) = (identity(expected), identity(row));
        if want != got {
            return Err(format!(
                "atom {} is {} where model {first_number} has {}",
                i + 1,
                got.join(" ").trim(),
                want.join(" ").trim()
            ));
        }
        let [x, y, z] = xyz.map(|idx| idx.and_then(|c| row.get(c)?.parse::<f64>().ok()));
        let (Some(x), Some(y), Some(z)) = (x, y, z) else {
            return Err(format!("atom {} has no readable coordinates", i + 1));
        };
        positions.push(Point3::new(x, y, z));
    }
    Ok(positions)
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

    for row in rows {
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

/// A `_chem_comp_bond.value_order`. The aromatic flag is ignored: the PDB
/// flags Kekulé orders too, and keeping them avoids half-aromatic rings.
fn chem_comp_bond_order(value: Option<&str>) -> BondOrder {
    match value.map(str::to_ascii_lowercase).as_deref() {
        Some("doub") => BondOrder::Double,
        Some("trip") => BondOrder::Triple,
        Some("quad") => BondOrder::Quadruple,
        Some("arom") => BondOrder::Aromatic,
        _ => BondOrder::Single,
    }
}

/// Bonds every residue instance whose name has `_chem_comp_bond` rows,
/// matched by atom name within that instance only (#406). A name the
/// instance lacks is skipped; an already-bonded pair keeps its order.
fn apply_chem_comp_bonds(
    mol: &mut Molecule,
    (tags, rows): &AtomSiteRows,
) -> Result<(), MmcifError> {
    let tag = |suffix: &str| find_tag(tags, &format!("_chem_comp_bond.{suffix}"));
    let (comp_id, atom_1, atom_2) = (tag("comp_id"), tag("atom_id_1"), tag("atom_id_2"));
    let value_order = tag("value_order");

    let mut by_comp: HashMap<&str, Vec<(&str, &str, BondOrder)>> = HashMap::new();
    for row in rows {
        if let (Some(comp), Some(a), Some(b)) = (
            cif_value(row, comp_id),
            cif_value(row, atom_1),
            cif_value(row, atom_2),
        ) {
            let order = chem_comp_bond_order(cif_value(row, value_order));
            by_comp.entry(comp).or_default().push((a, b, order));
        }
    }

    let mut pairs = Vec::new();
    for residue in mol.residues() {
        let Some(entries) = by_comp.get(residue.name.as_str()) else {
            continue;
        };
        // A name maps to several atoms when it has alternate locations.
        let mut by_name: HashMap<&str, Vec<usize>> = HashMap::new();
        for i in residue.atoms.clone() {
            if let Some(name) = mol.site(i).and_then(|s| s.name.as_deref()) {
                by_name.entry(name).or_default().push(i);
            }
        }
        for &(a, b, order) in entries {
            let (Some(xs), Some(ys)) = (by_name.get(a), by_name.get(b)) else {
                continue;
            };
            for &x in xs {
                for &y in ys {
                    if x != y && alt_locs_compatible(mol, x, y) {
                        pairs.push((x, y, order));
                    }
                }
            }
        }
    }

    for (x, y, order) in pairs {
        if !mol.graph().has_edge(x, y) {
            mol.add_bond(Bond::new(x, y, order))
                .map_err(|e| MmcifError::ParseError(e.to_string()))?;
        }
    }
    Ok(())
}

/// Two atoms can bond unless they sit in different alternate locations.
fn alt_locs_compatible(mol: &Molecule, x: usize, y: usize) -> bool {
    let alt = |i| mol.site(i).and_then(|s| s.alt_loc);
    match (alt(x), alt(y)) {
        (Some(p), Some(q)) => p == q,
        _ => true,
    }
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
    loops: &CifLoops,
    options: &MmcifReadOptions,
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

    // The bond pipeline: each step only adds pairs not yet bonded.
    if options.chem_comp_bonds
        && let Some(table) = loops.get("_chem_comp_bond")
    {
        apply_chem_comp_bonds(&mut mol, table)?;
    }
    // A no-op on a bondless read, so the default output is unchanged.
    mol.calculate_implicit_hydrogens_where_bonded();

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::mmcif::parse_mmcif_with_options;

    const ON: MmcifReadOptions = MmcifReadOptions {
        chem_comp_bonds: true,
        multi_frame: crate::io::options::MultiFrameMode::Molecules,
    };

    /// An mmCIF block with the given `_chem_comp_bond` rows and
    /// `(name, alt_id, comp, seq)` atoms, all in chain A.
    fn block(bond_rows: &[&str], atoms: &[(&str, &str, &str, i32)]) -> String {
        let mut text = String::from("data_test\n");
        if !bond_rows.is_empty() {
            text.push_str(
                "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
                 _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n",
            );
            for row in bond_rows {
                text.push_str(row);
                text.push('\n');
            }
        }
        text.push_str(
            "loop_\n_atom_site.group_PDB\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
             _atom_site.label_alt_id\n_atom_site.auth_comp_id\n_atom_site.auth_seq_id\n\
             _atom_site.auth_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        );
        for (i, (name, alt, comp, seq)) in atoms.iter().enumerate() {
            let element = &name[..1];
            text.push_str(&format!(
                "HETATM {element} {name} {alt} {comp} {seq} A {i}.0 1.0 1.0\n"
            ));
        }
        text
    }

    fn bonds(mol: &Molecule) -> Vec<(usize, usize, BondOrder)> {
        let mut out: Vec<_> = mol
            .bonds()
            .iter()
            .map(|b| {
                let (x, y) = b.atoms();
                (x.min(y), x.max(y), b.order())
            })
            .collect();
        out.sort_by_key(|&(x, y, _)| (x, y));
        out
    }

    #[test]
    fn test_a_doub_row_gives_a_double_bond() {
        let text = block(
            &["LIG C1 O1 doub"],
            &[("C1", ".", "LIG", 1), ("O1", ".", "LIG", 1)],
        );
        let mol = parse_mmcif_with_options(&text, &ON).unwrap();
        assert_eq!(bonds(&mol), vec![(0, 1, BondOrder::Double)]);
        assert_eq!(mol.atom(0).hydrogens(), Some(2));
    }

    #[test]
    fn test_the_option_off_reads_no_bonds() {
        let text = block(
            &["LIG C1 O1 doub"],
            &[("C1", ".", "LIG", 1), ("O1", ".", "LIG", 1)],
        );
        let mol = parse_mmcif_with_options(&text, &MmcifReadOptions::default()).unwrap();
        assert_eq!(mol.num_bonds(), 0);
        assert_eq!(mol.atom(0).hydrogens(), None);
    }

    #[test]
    fn test_a_row_naming_a_missing_hydrogen_is_skipped() {
        let text = block(
            &["LIG C1 O1 sing", "LIG O1 HO1 sing"],
            &[("C1", ".", "LIG", 1), ("O1", ".", "LIG", 1)],
        );
        let mol = parse_mmcif_with_options(&text, &ON).unwrap();
        assert_eq!(bonds(&mol), vec![(0, 1, BondOrder::Single)]);
    }

    #[test]
    fn test_a_comp_with_no_rows_gets_no_bonds() {
        let text = block(
            &["LIG C1 O1 sing"],
            &[("C1", ".", "ETH", 1), ("O1", ".", "ETH", 1)],
        );
        let mol = parse_mmcif_with_options(&text, &ON).unwrap();
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_two_instances_are_bonded_separately() {
        let text = block(
            &["LIG C1 O1 sing"],
            &[
                ("C1", ".", "LIG", 1),
                ("O1", ".", "LIG", 1),
                ("C1", ".", "LIG", 2),
                ("O1", ".", "LIG", 2),
            ],
        );
        let mol = parse_mmcif_with_options(&text, &ON).unwrap();
        assert_eq!(
            bonds(&mol),
            vec![(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)]
        );
    }

    #[test]
    fn test_alternate_locations_bond_only_within_their_own_location() {
        let text = block(
            &["LIG C1 O1 sing"],
            &[
                ("C1", ".", "LIG", 1),
                ("O1", "A", "LIG", 1),
                ("O1", "B", "LIG", 1),
                ("C1", "A", "LIG", 2),
                ("O1", "B", "LIG", 2),
            ],
        );
        let mol = parse_mmcif_with_options(&text, &ON).unwrap();
        assert_eq!(
            bonds(&mol),
            vec![(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)]
        );
    }

    #[test]
    fn test_value_orders_map_case_insensitively_and_unknown_is_single() {
        assert_eq!(chem_comp_bond_order(Some("DOUB")), BondOrder::Double);
        assert_eq!(chem_comp_bond_order(Some("trip")), BondOrder::Triple);
        assert_eq!(chem_comp_bond_order(Some("quad")), BondOrder::Quadruple);
        assert_eq!(chem_comp_bond_order(Some("arom")), BondOrder::Aromatic);
        assert_eq!(chem_comp_bond_order(Some("delo")), BondOrder::Single);
        assert_eq!(chem_comp_bond_order(None), BondOrder::Single);
    }

    #[test]
    fn test_a_repeated_pair_keeps_the_first_order() {
        let text = block(
            &["LIG C1 O1 doub", "LIG O1 C1 sing"],
            &[("C1", ".", "LIG", 1), ("O1", ".", "LIG", 1)],
        );
        let mol = parse_mmcif_with_options(&text, &ON).unwrap();
        assert_eq!(bonds(&mol), vec![(0, 1, BondOrder::Double)]);
    }
}
