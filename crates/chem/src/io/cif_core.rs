//! CIF core — the small-molecule crystallography dictionary over the same
//! CIF grammar [`crate::io::mmcif`] already tokenizes, but a different
//! vocabulary entirely (#320).
//!
//! Old-style flat tags (`_cell_length_a`, `_atom_site_fract_x`) rather than
//! mmCIF's dot-namespaced ones (`_cell.length_a`, `_atom_site.Cartn_x`) —
//! that difference is also what tells the two dictionaries apart in a
//! `.cif` file with no other way to name which one it is; see
//! [`is_small_molecule_cif`], consulted by [`crate::io::open`] and by the
//! CLI's own `bin/chem/stream.rs::read_input` before either commits to
//! mmCIF.
//!
//! Three real differences from mmCIF:
//!
//! **Coordinates are fractional**, converted through [`UnitCell::to_cartesian`]
//! on read and back through [`UnitCell::to_fractional`] on write. A
//! cell-less molecule has nowhere to convert to or from, so its raw
//! Cartesian coordinates are written into (and read back from) the
//! `_atom_site_fract_*` columns unconverted -- a disclosed escape hatch
//! (`WriteFn` cannot fail, so *something* has to go in those columns),
//! mirroring [`crate::io::gro`]'s own all-zero-box-for-no-cell precedent.
//!
//! **Symmetry is operator strings, and the file states only the asymmetric
//! unit.** This module reads and writes exactly that -- no expansion into
//! the full unit cell. A `_symmetry_equiv_pos_as_xyz`/
//! `_space_group_symop_operation_xyz` loop, if present, is parsed only far
//! enough to be discarded, the same treatment mmCIF gives `_struct_conn`.
//! `SpaceGroup`'s number/symbol are still read -- no operator strings
//! needed for that.
//!
//! **Numbers carry an inline estimated standard deviation**, `3.8402(1)`,
//! on nearly every value in a real file -- this module's own `strip_esd`
//! trims it before every `.parse::<f64>()`.
//!
//! **`_atom_site_U_iso_or_equiv` is converted to a B-factor** (`B = 8π²U`)
//! and stored in [`AtomSite::b_factor`] directly, rather than a dedicated
//! field -- it is common, real data, and the conversion is exact. The cost:
//! a read-then-write cycle re-derives U from B rather than preserving the
//! original bit-for-bit, the same kind of disclosed precision boundary
//! `cif_model::RowPrecision` documents elsewhere.
//!
//! No chains or residues: this dictionary has no `label_asym_id`/
//! `auth_seq_id` machinery, so [`Molecule::set_topology`] is simply never
//! called here, the same as [`crate::io::xyz`].

use std::collections::HashMap;

use crate::core::atom::Atom;
use crate::core::cell::{SpaceGroup, UnitCell};
use crate::core::geometry::{Point3, is_placeholder_3d};
use crate::core::molecule::Molecule;
use crate::core::site::AtomSite;
use crate::io::cif_model::{cif_value, element_from_symbol, find_tag};
use crate::io::errors::CifCoreError;
use crate::io::mmcif::tokenize_line;

/// `8π²` — converts an isotropic mean-square displacement `U` (Å²) into the
/// equivalent B-factor (Å²).
const U_TO_B: f64 = 8.0 * std::f64::consts::PI * std::f64::consts::PI;

/// Trims a trailing estimated-standard-deviation group -- `3.8402(1)` reads
/// as `3.8402`. A CIF numeric token's only legal parenthetical content is a
/// trailing esd digit group, so splitting on the first `(` is exactly as
/// robust as a purpose-built parser here.
fn strip_esd(raw: &str) -> &str {
    match raw.split_once('(') {
        Some((value, _)) => value,
        None => raw,
    }
}

fn parse_f64(raw: &str) -> Option<f64> {
    strip_esd(raw).trim().parse().ok()
}

/// Whether `text` looks like a small-molecule CIF-core file rather than
/// mmCIF -- scans for the first atom-site tag of either dictionary's shape
/// and lets whichever appears first decide. `false` (mmCIF, the existing
/// default) when neither appears, which keeps this an override only when
/// there's real evidence, never a guess.
///
/// `pub`, not `pub(crate)`: called from both `io::open::resolve_format`
/// (inside this crate) and `bin/chem/stream.rs::read_input` (the CLI
/// binary, a separate crate target that only sees this library's public
/// API).
pub fn is_small_molecule_cif(text: &str) -> bool {
    for line in text.lines() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("_atom_site.") {
            return false;
        }
        if trimmed.starts_with("_atom_site_") {
            return true;
        }
    }
    false
}

fn read_cell(singles: &HashMap<String, String>) -> Option<UnitCell> {
    let get = |key: &str| singles.get(key).and_then(|s| parse_f64(s));
    Some(UnitCell::new(
        get("_cell_length_a")?,
        get("_cell_length_b")?,
        get("_cell_length_c")?,
        get("_cell_angle_alpha")?,
        get("_cell_angle_beta")?,
        get("_cell_angle_gamma")?,
    ))
}

fn read_space_group(singles: &HashMap<String, String>) -> Option<SpaceGroup> {
    let symbol = singles
        .get("_symmetry_space_group_name_H-M")
        .or_else(|| singles.get("_space_group_name_H-M"))
        .cloned();
    let number = singles
        .get("_symmetry_Int_Tables_number")
        .or_else(|| singles.get("_space_group_IT_number"))
        .and_then(|s| s.trim().parse::<u16>().ok());
    (symbol.is_some() || number.is_some()).then_some(SpaceGroup { number, symbol })
}

fn build_molecule(
    singles: &HashMap<String, String>,
    atom_site: Option<(&[String], &[Vec<String>])>,
) -> Result<Molecule, CifCoreError> {
    let mut mol = Molecule::new();
    let mut sites = Vec::new();
    let mut coords = Vec::new();

    let cell = read_cell(singles);

    if let Some((tags, rows)) = atom_site {
        let type_symbol = find_tag(tags, "_atom_site_type_symbol").ok_or_else(|| {
            CifCoreError::ParseError("_atom_site loop has no type_symbol column".to_string())
        })?;
        let label = find_tag(tags, "_atom_site_label");
        let fract_x = find_tag(tags, "_atom_site_fract_x").ok_or_else(|| {
            CifCoreError::ParseError("_atom_site loop has no fract_x column".to_string())
        })?;
        let fract_y = find_tag(tags, "_atom_site_fract_y").ok_or_else(|| {
            CifCoreError::ParseError("_atom_site loop has no fract_y column".to_string())
        })?;
        let fract_z = find_tag(tags, "_atom_site_fract_z").ok_or_else(|| {
            CifCoreError::ParseError("_atom_site loop has no fract_z column".to_string())
        })?;
        let occupancy = find_tag(tags, "_atom_site_occupancy");
        let b_iso = find_tag(tags, "_atom_site_B_iso_or_equiv");
        let u_iso = find_tag(tags, "_atom_site_U_iso_or_equiv");

        for row in rows {
            let symbol = row.get(type_symbol).ok_or_else(|| {
                CifCoreError::InvalidAtomRow("row has no type_symbol value".to_string())
            })?;
            let element = element_from_symbol(symbol)
                .ok_or_else(|| CifCoreError::InvalidElement(symbol.to_string()))?;

            let x = row
                .get(fract_x)
                .and_then(|s| parse_f64(s))
                .ok_or_else(|| CifCoreError::InvalidAtomRow(row.join(" ")))?;
            let y = row
                .get(fract_y)
                .and_then(|s| parse_f64(s))
                .ok_or_else(|| CifCoreError::InvalidAtomRow(row.join(" ")))?;
            let z = row
                .get(fract_z)
                .and_then(|s| parse_f64(s))
                .ok_or_else(|| CifCoreError::InvalidAtomRow(row.join(" ")))?;
            let frac = Point3::new(x, y, z);
            let cart = cell.map(|c| c.to_cartesian(frac)).unwrap_or(frac);

            mol.add_atom(Atom::new(element));
            coords.push(cart);
            sites.push(AtomSite {
                name: cif_value(row, label).map(str::to_string),
                occupancy: cif_value(row, occupancy).and_then(parse_f64),
                b_factor: cif_value(row, b_iso).and_then(parse_f64).or_else(|| {
                    cif_value(row, u_iso)
                        .and_then(parse_f64)
                        .map(|u| u * U_TO_B)
                }),
                ..AtomSite::empty()
            });
        }
    }

    // Zeros here are what a format with no room to say "unknown" writes
    // for a molecule with no conformer, the same trap #270 found in mmCIF.
    if !is_placeholder_3d(&coords) {
        mol.set_coords3(coords)
            .map_err(|e| CifCoreError::ParseError(e.to_string()))?;
    }
    mol.set_sites(sites)
        .map_err(|e| CifCoreError::ParseError(e.to_string()))?;

    if let Some(cell) = cell {
        mol.set_cell(cell)
            .map_err(|e| CifCoreError::ParseError(e.to_string()))?;
    }
    if let Some(sg) = read_space_group(singles) {
        mol.set_space_group(sg);
    }

    if mol.num_atoms() == 0 && atom_site.is_none() {
        return Err(CifCoreError::NoAtoms);
    }

    Ok(mol)
}

/// Parses one CIF-core `data_` block (or a whole file with only one).
pub fn parse_cif_core(text: &str) -> Result<Molecule, CifCoreError> {
    let mut singles: HashMap<String, String> = HashMap::new();
    let mut atom_site: Option<(Vec<String>, Vec<Vec<String>>)> = None;

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

            // Any tag, not just the first -- a real file's atom-site loop
            // does not reliably lead with `_atom_site_label`, unlike
            // mmCIF's own `_atom_site.`-prefixed convention this mirrors.
            if tags.iter().any(|t| t.starts_with("_atom_site_")) {
                atom_site = Some((tags, rows));
            }
            // Every other loop -- the symmetry-operator one included -- is
            // out of scope, see the module doc.
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

    build_molecule(
        &singles,
        atom_site
            .as_ref()
            .map(|(t, r)| (t.as_slice(), r.as_slice())),
    )
}

/// Writes one CIF-core `data_` block.
pub fn write_cif_core(mol: &Molecule) -> String {
    let mut out = String::from("data_chem\n");

    if let Some(cell) = mol.cell() {
        out.push_str(&format!(
            "_cell_length_a {:.4}\n_cell_length_b {:.4}\n_cell_length_c {:.4}\n\
             _cell_angle_alpha {:.2}\n_cell_angle_beta {:.2}\n_cell_angle_gamma {:.2}\n",
            cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma
        ));
    }
    if let Some(sg) = mol.space_group() {
        if let Some(symbol) = &sg.symbol {
            out.push_str(&format!("_symmetry_space_group_name_H-M '{symbol}'\n"));
        }
        if let Some(number) = sg.number {
            out.push_str(&format!("_symmetry_Int_Tables_number {number}\n"));
        }
    }

    out.push_str("loop_\n");
    for tag in [
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
        "_atom_site_B_iso_or_equiv",
    ] {
        out.push_str(tag);
        out.push('\n');
    }

    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let cart = mol.coord3(i).unwrap_or(Point3::ORIGIN);
        let frac = mol.cell().map(|c| c.to_fractional(cart)).unwrap_or(cart);

        let name = site
            .and_then(|s| s.name.as_deref())
            .unwrap_or(atom.element().symbol());
        let occupancy = site.and_then(|s| s.occupancy).unwrap_or(1.0);
        let b_factor = site.and_then(|s| s.b_factor).unwrap_or(0.0);

        out.push_str(&format!(
            // B-factor gets more precision than occupancy: it is what
            // `_atom_site_U_iso_or_equiv` gets converted into on read
            // (`U_TO_B`), and CIF (unlike PDB's fixed-width columns) has no
            // reason to round that conversion's output back down again --
            // doing so at only two decimals is what used to turn a stated
            // U=0.0140 into 0.0141 after a round trip.
            "{name} {element} {x:.5} {y:.5} {z:.5} {occupancy:.2} {b_factor:.6}\n",
            element = atom.element().symbol(),
            x = frac.x,
            y = frac.y,
            z = frac.z,
        ));
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // Corpus fixtures rather than duplicated inline text, so there is one
    // copy of this content to keep correct -- see `tests/corpus/cif_core/`.
    const QUARTZ: &str = include_str!("../../tests/corpus/cif_core/quartz.cif");
    const QUARTZ_WITH_SYMOPS: &str =
        include_str!("../../tests/corpus/cif_core/symmetry-operators.cif");

    #[test]
    fn test_a_cell_and_esd_bearing_structure_round_trips() {
        let mol = parse_cif_core(QUARTZ).expect("valid CIF core");
        assert_eq!(mol.num_atoms(), 2);
        let cell = mol.cell().expect("has a cell");
        assert!(
            (cell.a - 4.9134).abs() < 1e-9,
            "esd not stripped: {}",
            cell.a
        );
        assert_eq!(
            mol.space_group().and_then(|g| g.symbol.as_deref()),
            Some("P 32 2 1")
        );
        assert_eq!(mol.space_group().and_then(|g| g.number), Some(154));

        let written = write_cif_core(&mol);
        let back = parse_cif_core(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.cell(), mol.cell());
    }

    #[test]
    fn test_u_iso_converts_to_the_equivalent_b_factor() {
        let mol = parse_cif_core(QUARTZ).expect("valid CIF core");
        let expected = 0.0090 * U_TO_B;
        let got = mol
            .site(0)
            .and_then(|s| s.b_factor)
            .expect("has a b-factor");
        assert!((got - expected).abs() < 1e-6, "{got} vs {expected}");
    }

    #[test]
    fn test_a_symmetry_operator_loop_is_parsed_and_ignored_not_expanded() {
        // The decision this story settled on: the asymmetric unit only,
        // exactly as stated -- one atom in, one atom out, not three.
        let mol = parse_cif_core(QUARTZ_WITH_SYMOPS).expect("valid CIF core");
        assert_eq!(mol.num_atoms(), 1);
    }

    #[test]
    fn test_a_cell_less_molecule_round_trips_through_raw_cartesian() {
        let text = "data_x\nloop_\n_atom_site_label\n_atom_site_type_symbol\n\
                    _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                    C1 C 1.5 2.5 3.5\n";
        let mol = parse_cif_core(text).expect("valid CIF core");
        assert!(mol.cell().is_none());
        assert_eq!(mol.coord3(0), Some(Point3::new(1.5, 2.5, 3.5)));

        let written = write_cif_core(&mol);
        let back = parse_cif_core(&written).expect("round trips");
        assert_eq!(back.coord3(0), mol.coord3(0));
    }

    #[test]
    fn test_an_unstripped_esd_is_a_clear_error_not_a_panic() {
        // strip_esd is what makes 3.8402(1) parse at all -- confirm the raw
        // token really would fail on its own, so the strip is proven
        // necessary rather than merely present.
        assert!("3.8402(1)".parse::<f64>().is_err());
        assert_eq!(strip_esd("3.8402(1)"), "3.8402");
    }

    #[test]
    fn test_is_small_molecule_cif_distinguishes_the_two_dictionaries() {
        assert!(is_small_molecule_cif(QUARTZ));
        assert!(!is_small_molecule_cif(
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\nC 0.0\n"
        ));
        assert!(!is_small_molecule_cif("data_x\n_cell_length_a 1.0\n"));
    }

    #[test]
    fn test_garbage_text_reports_no_atoms_instead_of_an_empty_molecule() {
        for input in [
            "",
            "data_x\n_cell_length_a 10.000\n",
            "not a cif file at all\n",
        ] {
            assert!(
                matches!(parse_cif_core(input), Err(CifCoreError::NoAtoms)),
                "{input:?} should report NoAtoms"
            );
        }
    }

    #[test]
    fn test_an_atom_site_loop_with_zero_rows_is_a_legitimately_empty_structure() {
        let text = "data_x\nloop_\n_atom_site_type_symbol\n_atom_site_fract_x\n\
                    _atom_site_fract_y\n_atom_site_fract_z\n";
        let mol = parse_cif_core(text).expect("a present-but-empty atom_site loop is legal");
        assert_eq!(mol.num_atoms(), 0);
    }

    #[test]
    fn test_an_unrecognised_element_is_a_clear_error() {
        let text = "data_x\nloop_\n_atom_site_type_symbol\n_atom_site_fract_x\n\
                    _atom_site_fract_y\n_atom_site_fract_z\nXx 0.0 0.0 0.0\n";
        let err = parse_cif_core(text).unwrap_err();
        assert!(matches!(err, CifCoreError::InvalidElement(_)), "{err}");
    }

    #[test]
    fn test_a_malformed_row_is_a_clear_error_not_a_panic() {
        let text = "data_x\nloop_\n_atom_site_type_symbol\n_atom_site_fract_x\n\
                    _atom_site_fract_y\n_atom_site_fract_z\nC not-a-number 0.0 0.0\n";
        let err = parse_cif_core(text).unwrap_err();
        assert!(matches!(err, CifCoreError::InvalidAtomRow(_)), "{err}");
    }
}
