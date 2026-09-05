//! Mol2 (Tripos) — text, `@<TRIPOS>`-delimited sections. Heavily used in
//! docking-prep pipelines (#225).
//!
//! **`am` (amide) bonds have nothing to map to.** [`BondOrder`] has exactly
//! five variants (`Single`/`Double`/`Triple`/`Quadruple`/`Aromatic`) and no
//! amide concept — reading `am` gives [`BondOrder::Single`] (an amide's
//! formal bond order), and this crate never *writes* `am`, always a
//! concrete order. Stated as a limitation, not silently invented otherwise.
//!
//! **No aromaticity perception on read, and no Kekulé conversion on
//! write.** SYBYL atom types state aromaticity directly per atom (`C.ar`),
//! and Mol2's bond records have a native `ar` type — unlike SDF's V2000,
//! which has neither and needs [`crate::io::aromaticity::detect_aromaticity`]/
//! [`crate::io::aromaticity::kekulize`] to bridge the gap. Running either
//! here would be redundant at best, and could contradict a file's own
//! aromatic-bond assertion at worst. An atom's own aromaticity and its
//! bonds' `ar` type are trusted independently, exactly as written — no
//! cross-validation, matching this crate's "read only what's explicit"
//! convention rather than reconciling a real file's internal
//! inconsistencies.
//!
//! **SYBYL atom types are a curated common subset, not the full ~70-type
//! spec.** Ordinary organic chemistry (C/N/O/S/P with standard
//! hybridization suffixes), halogens, common metals, and solvent-model
//! variants (`O.spc`, `H.t3p`) folding into their plain equivalent are all
//! covered. True dummy/pseudo-atom types (`Du`, `Any`, `Hal`, `Het`, `LP`)
//! are a query-molecule concept this crate has no model for at all — the
//! same reasoning CXSMARTS was deferred under (#221) — and are a clear
//! error naming the type, not a guess. Writing back a type this crate
//! narrowed on read (`N.pl3`/`N.am` → `SP2` → `N.2`, `O.co2` → `SP2`/`SP3`
//! → `O.2`/`O.3`) does not reproduce the original spelling — the same
//! "canonical, not byte-identical" honesty #220's SMILES writer already
//! established.
//!
//! **`@<TRIPOS>SUBSTRUCTURE` is parsed** for chain identity and a
//! hetero-like flag: `is_hetero` is `true` iff a substructure's
//! `subst_type` is present and not (case-insensitively) the literal
//! `"RESIDUE"` — a stated convention, not discovered data, matching how
//! docking-prep tools typically label ligands/cofactors `"GROUP"` and
//! standard polymer residues `"RESIDUE"`. Absent the section (or a
//! substructure's own `subst_type`/`chain` fields), residues fall into one
//! synthetic chain with `is_hetero: false`.
//!
//! A file may hold several `@<TRIPOS>MOLECULE` blocks back to back (a
//! ligand library, this format's normal batch shape) — splitting those is
//! the reader's/supplier's job
//! ([`crate::io::reader::read_mol2_with_options`],
//! [`crate::io::supplier::Mol2Supplier`]), the same division every prior
//! format has. This module parses and writes exactly one already-isolated
//! block.

use std::collections::HashMap;

use crate::core::atom::{Atom, Element, Hybridization};
use crate::core::bond::{Bond, BondOrder};
use crate::core::cell::{SpaceGroup, UnitCell};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::core::site::AtomSite;
use crate::io::errors::Mol2Error;

/// One `@<TRIPOS>SECTION` header and the lines making up its body.
struct Section<'a> {
    name: &'a str,
    lines: &'a [&'a str],
}

fn is_section_header(line: &str) -> bool {
    line.trim_start().starts_with("@<TRIPOS>")
}

fn section_name(line: &str) -> &str {
    line.trim_start()
        .strip_prefix("@<TRIPOS>")
        .unwrap_or("")
        .trim()
}

/// Splits one `@<TRIPOS>MOLECULE` block into its sections. More robust
/// than trusting the `MOLECULE` record's declared counts (real files
/// sometimes have them wrong) — a section's body simply runs until the
/// next `@<TRIPOS>` line or the block's end.
fn split_sections<'a>(lines: &'a [&'a str]) -> Vec<Section<'a>> {
    let mut sections = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if !is_section_header(lines[i]) {
            i += 1;
            continue;
        }
        let name_line = lines[i];
        i += 1;
        let start = i;
        while i < lines.len() && !is_section_header(lines[i]) {
            i += 1;
        }
        sections.push(Section {
            name: section_name(name_line),
            lines: &lines[start..i],
        });
    }
    sections
}

fn tokens(line: &str) -> Vec<&str> {
    line.split_whitespace().collect()
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

/// Resolves a SYBYL atom type (`"C.ar"`, `"N.pl3"`, `"Cl"`, ...) into an
/// element, hybridization and aromaticity — the curated subset described
/// in the module doc. Elements needing a suffix in the real spec but
/// found without one here still resolve (defaulting to `Unknown`
/// hybridization) rather than erroring, since a curated subset should fail
/// only on what it genuinely cannot represent (dummy/pseudo types), not on
/// every spelling variant of what it can.
fn resolve_sybyl_type(sybyl: &str) -> Result<(Element, Hybridization, bool), Mol2Error> {
    let (elem_part, suffix) = match sybyl.split_once('.') {
        Some((e, s)) => (e, Some(s)),
        None => (sybyl, None),
    };

    match elem_part {
        "Du" | "Any" | "Hal" | "Het" | "LP" => {
            return Err(Mol2Error::UnsupportedAtomType(sybyl.to_string()));
        }
        _ => {}
    }

    let element = element_from_symbol(elem_part)
        .ok_or_else(|| Mol2Error::UnsupportedAtomType(sybyl.to_string()))?;

    let (hybridization, is_aromatic) = match suffix {
        Some("ar") => (Hybridization::SP2, true),
        Some("1") => (Hybridization::SP, false),
        Some("2") | Some("pl3") | Some("am") | Some("cat") | Some("co2") => {
            (Hybridization::SP2, false)
        }
        Some("3") | Some("O") | Some("O2") | Some("spc") | Some("t3p") => {
            (Hybridization::SP3, false)
        }
        Some("4") => (Hybridization::SP3, false),
        _ => (Hybridization::Unknown, false),
    };

    Ok((element, hybridization, is_aromatic))
}

/// The reverse of [`resolve_sybyl_type`]. Only these five elements take a
/// hybridization suffix in the curated subset; every other element
/// (halogens, `H`, metals) writes as its bare symbol.
fn sybyl_type_for(atom: &Atom) -> String {
    let symbol = atom.element().symbol();
    let takes_suffix = matches!(symbol, "C" | "N" | "O" | "S" | "P");
    if !takes_suffix {
        return symbol.to_string();
    }
    if atom.is_aromatic() {
        return format!("{symbol}.ar");
    }
    let suffix = match atom.hybridization() {
        Hybridization::SP => "1",
        Hybridization::SP2 => "2",
        // SP3D/SP3D2/Unknown all fall back to the most common suffix,
        // since these five elements always need one in the real spec.
        Hybridization::SP3
        | Hybridization::SP3D
        | Hybridization::SP3D2
        | Hybridization::Unknown => "3",
    };
    format!("{symbol}.{suffix}")
}

fn bond_order_from_type(bond_type: &str) -> Option<BondOrder> {
    match bond_type {
        "1" => Some(BondOrder::Single),
        "2" => Some(BondOrder::Double),
        "3" => Some(BondOrder::Triple),
        "ar" => Some(BondOrder::Aromatic),
        // Amide: no equivalent, see the module doc -- formal bond order.
        "am" => Some(BondOrder::Single),
        // Dummy/unknown, rare in practice: same fallback.
        "du" | "un" => Some(BondOrder::Single),
        // "Not connected" states there is no bond at all.
        "nc" => None,
        _ => Some(BondOrder::Single),
    }
}

fn bond_type_for(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "1",
        BondOrder::Double => "2",
        BondOrder::Triple => "3",
        BondOrder::Aromatic => "ar",
        // No Mol2 equivalent; written as the closest concrete order.
        BondOrder::Quadruple => "1",
    }
}

struct SubstMeta {
    chain: String,
    is_hetero: bool,
}

#[derive(Clone, PartialEq, Eq)]
struct ResidueKey {
    chain_id: String,
    name: String,
    sequence: i32,
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
            insertion_code: None,
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

/// Parses one `@<TRIPOS>MOLECULE` block.
pub fn parse_mol2(text: &str) -> Result<Molecule, Mol2Error> {
    let lines: Vec<&str> = text.lines().collect();
    let sections = split_sections(&lines);

    let mut mol = Molecule::new();
    let mut coords = Vec::new();
    let mut sites = Vec::new();
    let mut keys = Vec::new();
    let mut serial_to_index: HashMap<i64, usize> = HashMap::new();
    let mut subst_meta: HashMap<i64, SubstMeta> = HashMap::new();

    for section in &sections {
        if section.name == "MOLECULE"
            && let Some(name_line) = section.lines.first()
            && !name_line.trim().is_empty()
        {
            mol.set_name(name_line.trim().to_string());
        }
    }

    for section in &sections {
        if section.name != "SUBSTRUCTURE" {
            continue;
        }
        for line in section.lines {
            let fields = tokens(line);
            if fields.len() < 2 {
                continue;
            }
            let Ok(subst_id) = fields[0].parse::<i64>() else {
                continue;
            };
            let subst_type = fields.get(3).copied();
            let is_hetero = subst_type.is_some_and(|t| !t.eq_ignore_ascii_case("RESIDUE"));
            let chain = fields.get(5).map(|s| s.to_string()).unwrap_or_default();
            subst_meta.insert(subst_id, SubstMeta { chain, is_hetero });
        }
    }

    for section in &sections {
        if section.name != "ATOM" {
            continue;
        }
        for line in section.lines {
            let fields = tokens(line);
            if fields.len() < 6 {
                continue;
            }
            let atom_id: i64 = fields[0]
                .parse()
                .map_err(|_| Mol2Error::InvalidAtomLine(line.to_string()))?;
            let name = fields[1];
            let x: f64 = fields[2]
                .parse()
                .map_err(|_| Mol2Error::InvalidAtomLine(line.to_string()))?;
            let y: f64 = fields[3]
                .parse()
                .map_err(|_| Mol2Error::InvalidAtomLine(line.to_string()))?;
            let z: f64 = fields[4]
                .parse()
                .map_err(|_| Mol2Error::InvalidAtomLine(line.to_string()))?;
            let sybyl = fields[5];
            let subst_id: i64 = fields.get(6).and_then(|s| s.parse().ok()).unwrap_or(1);
            let subst_name = fields.get(7).copied().unwrap_or("UNK");
            let charge: Option<f64> = fields.get(8).and_then(|s| s.parse().ok());

            let (element, hybridization, is_aromatic) = resolve_sybyl_type(sybyl)?;
            let atom = Atom::new(element)
                .with_hybridization(hybridization)
                .with_aromatic(is_aromatic);
            let atom_idx = mol.add_atom(atom);
            serial_to_index.insert(atom_id, atom_idx);
            coords.push(Point3::new(x, y, z));
            sites.push(AtomSite {
                name: Some(name.to_string()),
                alt_loc: None,
                partial_charge: charge,
                occupancy: None,
                b_factor: None,
                radius: None,
            });

            let meta = subst_meta.get(&subst_id);
            keys.push(ResidueKey {
                chain_id: meta.map(|m| m.chain.clone()).unwrap_or_default(),
                name: subst_name.to_string(),
                sequence: subst_id as i32,
                is_hetero: meta.map(|m| m.is_hetero).unwrap_or(false),
            });
        }
    }

    // The 2D-vs-3D distinction is exactly the one `parse_sdf` already
    // makes: an all-zero z is a flat drawing, anything else is geometry.
    // Mol2 has no dimensionality header of its own either.
    if coords.iter().all(|p: &Point3| p.z == 0.0) {
        let flat: Vec<_> = coords.iter().map(|p| p.to_2d()).collect();
        mol.set_coords(flat)
            .map_err(|e| Mol2Error::ParseError(e.to_string()))?;
    } else {
        mol.set_coords3(coords)
            .map_err(|e| Mol2Error::ParseError(e.to_string()))?;
    }
    mol.set_sites(sites)
        .map_err(|e| Mol2Error::ParseError(e.to_string()))?;

    let (chains, residues) = group_into_chains_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| Mol2Error::ParseError(e.to_string()))?;

    for section in &sections {
        if section.name != "BOND" {
            continue;
        }
        for line in section.lines {
            let fields = tokens(line);
            if fields.len() < 4 {
                continue;
            }
            let origin: i64 = fields[1]
                .parse()
                .map_err(|_| Mol2Error::InvalidBondLine(line.to_string()))?;
            let target: i64 = fields[2]
                .parse()
                .map_err(|_| Mol2Error::InvalidBondLine(line.to_string()))?;
            let Some(order) = bond_order_from_type(fields[3]) else {
                continue;
            };
            let origin_idx = *serial_to_index
                .get(&origin)
                .ok_or_else(|| Mol2Error::InvalidBondLine(line.to_string()))?;
            let target_idx = *serial_to_index
                .get(&target)
                .ok_or_else(|| Mol2Error::InvalidBondLine(line.to_string()))?;
            mol.add_bond(Bond::new(origin_idx, target_idx, order))
                .map_err(|e| Mol2Error::ParseError(e.to_string()))?;
        }
    }

    for section in &sections {
        if section.name != "CRYSIN" {
            continue;
        }
        if let Some(line) = section.lines.first() {
            let fields = tokens(line);
            if fields.len() >= 6 {
                let parsed: Result<Vec<f64>, _> = fields[0..6].iter().map(|s| s.parse()).collect();
                if let Ok(values) = parsed {
                    let cell = UnitCell::new(
                        values[0], values[1], values[2], values[3], values[4], values[5],
                    );
                    mol.set_cell(cell)
                        .map_err(|e| Mol2Error::ParseError(e.to_string()))?;
                }
                if let Some(sg) = fields.get(6).and_then(|s| s.parse::<u16>().ok()) {
                    mol.set_space_group(SpaceGroup::from_number(sg));
                }
            }
        }
    }

    Ok(mol)
}

/// Writes one `@<TRIPOS>MOLECULE` block.
pub fn write_mol2(mol: &Molecule) -> String {
    let mut out = String::from("@<TRIPOS>MOLECULE\n");
    out.push_str(mol.name().unwrap_or("*"));
    out.push('\n');

    let num_bonds = mol.num_bonds();
    out.push_str(&format!("{} {} 0 0 0\n", mol.num_atoms(), num_bonds));
    out.push_str("SMALL\n");
    let has_charges = mol
        .sites()
        .is_some_and(|s| s.iter().any(|a| a.partial_charge.is_some()));
    out.push_str(if has_charges {
        "USER_CHARGES\n"
    } else {
        "NO_CHARGES\n"
    });

    if let Some(cell) = mol.cell() {
        out.push_str("@<TRIPOS>CRYSIN\n");
        let sg = mol.space_group().and_then(|g| g.number).unwrap_or(1);
        out.push_str(&format!(
            "{:.4} {:.4} {:.4} {:.4} {:.4} {:.4} {} 1\n",
            cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma, sg
        ));
    }

    out.push_str("@<TRIPOS>ATOM\n");
    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let residue = mol.residue_of(i);
        let name = site
            .and_then(|s| s.name.as_deref())
            .unwrap_or(atom.element().symbol());
        let p3 = mol.coord3(i);
        let p2 = mol.coord(i);
        let (x, y, z) = match (p3, p2) {
            (Some(p), _) => (p.x, p.y, p.z),
            (None, Some(p)) => (p.x, p.y, 0.0),
            (None, None) => (0.0, 0.0, 0.0),
        };
        let charge = site.and_then(|s| s.partial_charge).unwrap_or(0.0);
        let subst_id = residue.map(|r| r.sequence).unwrap_or(1);
        let subst_name = residue.map(|r| r.name.as_str()).unwrap_or("UNK");

        out.push_str(&format!(
            "{:>7} {:<8} {:>9.4} {:>9.4} {:>9.4} {:<8} {:>4} {:<8} {:>9.4}\n",
            i + 1,
            name,
            x,
            y,
            z,
            sybyl_type_for(atom),
            subst_id,
            subst_name,
            charge,
        ));
    }

    out.push_str("@<TRIPOS>BOND\n");
    for (i, bond) in mol.bonds().iter().enumerate() {
        out.push_str(&format!(
            "{:>6} {:>6} {:>6} {}\n",
            i + 1,
            bond.atom1() + 1,
            bond.atom2() + 1,
            bond_type_for(bond.order()),
        ));
    }

    out.push_str("@<TRIPOS>SUBSTRUCTURE\n");
    for (i, residue) in mol.residues().iter().enumerate() {
        let chain = mol
            .chains()
            .get(residue.chain_ix)
            .map(|c| c.id.as_str())
            .unwrap_or("");
        let root_atom = residue.atoms.start + 1;
        let subst_type = if residue.is_hetero {
            "GROUP"
        } else {
            "RESIDUE"
        };
        out.push_str(&format!(
            "{:>6} {:<8} {:>6} {} 1 {}\n",
            i + 1,
            residue.name,
            root_atom,
            subst_type,
            chain,
        ));
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const BENZENE_MOL2: &str = "\
@<TRIPOS>MOLECULE
benzene
6 6 1 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 1.396 0.000 C.ar 1 BNZ 0.0
2 C2 1.209 0.698 0.000 C.ar 1 BNZ 0.0
3 C3 1.209 -0.698 0.000 C.ar 1 BNZ 0.0
4 C4 0.000 -1.396 0.000 C.ar 1 BNZ 0.0
5 C5 -1.209 -0.698 0.000 C.ar 1 BNZ 0.0
6 C6 -1.209 0.698 0.000 C.ar 1 BNZ 0.0
@<TRIPOS>BOND
1 1 2 ar
2 2 3 ar
3 3 4 ar
4 4 5 ar
5 5 6 ar
6 6 1 ar
@<TRIPOS>SUBSTRUCTURE
1 BNZ 1 RESIDUE 1 A
";

    #[test]
    fn test_a_single_molecule_round_trips() {
        let mol = parse_mol2(BENZENE_MOL2).expect("valid Mol2");
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.num_bonds(), 6);
        assert!(mol.atoms()[0].is_aromatic());
        assert_eq!(mol.bonds()[0].order(), BondOrder::Aromatic);
        assert_eq!(mol.chains()[0].id, "A");
        assert!(!mol.residues()[0].is_hetero);

        let written = write_mol2(&mol);
        let back = parse_mol2(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.num_bonds(), mol.num_bonds());
        assert!(back.atoms()[0].is_aromatic());
    }

    #[test]
    fn test_partial_charges_round_trip() {
        let text = "\
@<TRIPOS>MOLECULE
charged
2 1 1 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 0.000 0.000 C.3 1 LIG -0.1500
2 O1 1.500 0.000 0.000 O.3 1 LIG -0.4000
@<TRIPOS>BOND
1 1 2 1
";
        let mol = parse_mol2(text).expect("valid Mol2");
        assert_eq!(mol.site(0).unwrap().partial_charge, Some(-0.15));
        assert_eq!(mol.site(1).unwrap().partial_charge, Some(-0.4));

        let written = write_mol2(&mol);
        let back = parse_mol2(&written).expect("round trips");
        assert_eq!(back.site(0).unwrap().partial_charge, Some(-0.15));
    }

    #[test]
    fn test_an_amide_bond_reads_as_a_plain_single_bond() {
        let text = "\
@<TRIPOS>MOLECULE
amide
2 1 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 0.000 0.000 C.2 1 LIG 0.0
2 N1 1.500 0.000 0.000 N.am 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 am
";
        let mol = parse_mol2(text).expect("valid Mol2");
        assert_eq!(mol.bonds()[0].order(), BondOrder::Single);
    }

    #[test]
    fn test_a_not_connected_bond_adds_no_edge() {
        let text = "\
@<TRIPOS>MOLECULE
nc
2 1 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 0.000 0.000 C.3 1 LIG 0.0
2 C2 1.500 0.000 0.000 C.3 1 LIG 0.0
@<TRIPOS>BOND
1 1 2 nc
";
        let mol = parse_mol2(text).expect("valid Mol2");
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_substructure_group_type_is_hetero() {
        let text = "\
@<TRIPOS>MOLECULE
lig
1 0 1 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 0.000 0.000 C.3 1 LIG 0.0
@<TRIPOS>SUBSTRUCTURE
1 LIG 1 GROUP 1 B
";
        let mol = parse_mol2(text).expect("valid Mol2");
        assert!(mol.residues()[0].is_hetero);
        assert_eq!(mol.chains()[0].id, "B");
    }

    #[test]
    fn test_no_substructure_section_falls_back_to_one_synthetic_chain() {
        let text = "\
@<TRIPOS>MOLECULE
lig
1 0 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0.000 0.000 0.000 C.3 1 LIG 0.0
";
        let mol = parse_mol2(text).expect("valid Mol2");
        assert_eq!(mol.chains().len(), 1);
        assert!(!mol.residues()[0].is_hetero);
    }

    #[test]
    fn test_crysin_round_trips_through_unit_cell_and_space_group() {
        let text = "\
@<TRIPOS>MOLECULE
cell-test
1 0 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>CRYSIN
10.0000 20.0000 30.0000 90.0000 90.0000 90.0000 19 1
@<TRIPOS>ATOM
1 C1 0.000 0.000 1.000 C.3 1 LIG 0.0
";
        let mol = parse_mol2(text).expect("valid Mol2");
        let cell = mol.cell().expect("has a cell");
        assert_eq!((cell.a, cell.b, cell.c), (10.0, 20.0, 30.0));
        assert_eq!(mol.space_group().and_then(|g| g.number), Some(19));

        let written = write_mol2(&mol);
        assert!(written.contains("@<TRIPOS>CRYSIN"), "{written}");
        let back = parse_mol2(&written).expect("round trips");
        assert_eq!(back.cell(), mol.cell());
    }

    #[test]
    fn test_an_unrecognised_dummy_atom_type_is_a_clear_error() {
        let text = "\
@<TRIPOS>MOLECULE
dummy
1 0 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 X1 0.000 0.000 0.000 Du 1 LIG 0.0
";
        let err = parse_mol2(text).unwrap_err();
        assert!(matches!(err, Mol2Error::UnsupportedAtomType(_)), "{err}");
    }

    #[test]
    fn test_a_malformed_coordinate_is_a_clear_error_not_a_panic() {
        let text = "\
@<TRIPOS>MOLECULE
bad
1 0 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 not-a-number 0.000 0.000 C.3 1 LIG 0.0
";
        let err = parse_mol2(text).unwrap_err();
        assert!(matches!(err, Mol2Error::InvalidAtomLine(_)), "{err}");
    }
}
