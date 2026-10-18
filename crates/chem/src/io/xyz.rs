//! XYZ — the simplest format in scope (#222): an atom count, a free-text
//! comment line, then one `<element> <x> <y> <z>` line per atom.
//!
//! Extended XYZ's `key="value"` comment line (`Lattice=`, `Properties=`,
//! per-atom velocity/force/charge columns) is accepted without erroring —
//! this format's comment line is free text either way, and extra columns
//! past the fourth are simply ignored — but none of it is specially parsed.
//! In particular **no unit cell is read or written**: `Lattice=` names three
//! Cartesian basis vectors in whatever orientation the file's writer used,
//! not necessarily this crate's fixed "a along x, b in the xy-plane"
//! convention ([`crate::core::cell::UnitCell`], #178). Reading it correctly
//! would mean rotating every atom's coordinate into that canonical frame
//! too, or [`crate::core::molecule::Molecule::cell`] and
//! [`crate::core::molecule::Molecule::coords3`] would silently disagree for
//! any file not already in that orientation — a correctness risk this story
//! deliberately does not take on, rather than getting it quietly wrong.
//!
//! A file may hold more than one frame back to back — an MD trajectory's
//! normal shape. Each frame's own atom count is a self-delimiting record
//! boundary, the same role SDF's `$$$$` terminator plays; splitting frames
//! is the reader's/supplier's job ([`crate::io::reader::read_xyz_with_options`],
//! [`crate::io::supplier::XyzSupplier`]), not this module's — this file
//! parses and writes exactly one frame, mirroring [`crate::io::sdf::parse_sdf`].
//!
//! XYZ states no bonds, so a molecule arrives with none unless
//! [`XyzReadOptions::perceive_bonds`] is set, which bonds atoms by distance
//! ([`crate::core::perceive::perceive_bonds`]), all as single bonds.

use crate::core::geometry::is_placeholder_3d;
use crate::core::molecule::Molecule;
use crate::core::perceive::perceive_bonds;
use crate::io::errors::XyzError;
use crate::io::options::XyzReadOptions;

/// Parses one XYZ frame: the count line, the comment line, then that many
/// atom lines.
pub fn parse_xyz(text: &str) -> Result<Molecule, XyzError> {
    let mut lines = text.lines();

    let count_line = lines
        .next()
        .ok_or_else(|| XyzError::ParseError("empty frame".to_string()))?;
    let count: usize = count_line
        .trim()
        .parse()
        .map_err(|_| XyzError::ParseError(format!("invalid atom count: {count_line:?}")))?;

    // The comment line is free text, always -- extended XYZ's `key="value"`
    // form included. Nothing here inspects it beyond taking it as the name.
    let comment = lines
        .next()
        .ok_or_else(|| XyzError::ParseError("missing comment line".to_string()))?
        .trim();

    let mut mol = Molecule::new();
    if !comment.is_empty() {
        mol.set_name(comment.to_string());
    }

    let mut coords = Vec::with_capacity(count);
    for _ in 0..count {
        let line = lines.next().ok_or(XyzError::AtomCountMismatch {
            expected: count,
            got: coords.len(),
        })?;
        let (atom, point) = parse_atom_line(line)?;
        mol.add_atom(atom);
        coords.push(point);
    }

    if lines.next().is_some() {
        // A frame this module was handed should be exactly `2 + count`
        // lines -- the caller (reader/supplier) is what finds that
        // boundary, so trailing content here means the split was wrong,
        // not that this file is malformed. Fatal either way.
        return Err(XyzError::ParseError(
            "more lines than the declared atom count accounts for".to_string(),
        ));
    }

    // Zeros in these columns are what a format with no room to say "unknown"
    // writes for a molecule that has no conformer, so believing them back is
    // how a converted molecule ended up undrawable (#270).
    if !is_placeholder_3d(&coords) {
        mol.set_coords3(coords)
            .map_err(|e| XyzError::ParseError(e.to_string()))?;
    }
    Ok(mol)
}

/// [`parse_xyz`] with XYZ's own read options.
pub fn parse_xyz_with_options(text: &str, options: &XyzReadOptions) -> Result<Molecule, XyzError> {
    let mut mol = parse_xyz(text)?;
    if options.perceive_bonds {
        perceive_bonds(&mut mol);
        mol.calculate_implicit_hydrogens_where_bonded();
    }
    Ok(mol)
}

/// `<element> <x> <y> <z> [ignored extra columns]`.
fn parse_atom_line(
    line: &str,
) -> Result<(crate::core::atom::Atom, crate::core::geometry::Point3), XyzError> {
    use crate::core::atom::{Atom, Element};
    use crate::core::elements::ELEMENT_SYMBOLS;
    use crate::core::geometry::Point3;

    let mut fields = line.split_whitespace();
    let symbol = fields
        .next()
        .ok_or_else(|| XyzError::InvalidAtomLine(line.to_string()))?;

    let atomic_number = ELEMENT_SYMBOLS
        .iter()
        .position(|&s| s == symbol)
        .ok_or_else(|| XyzError::InvalidElement(symbol.to_string()))?;
    let element = Element::new(atomic_number as u8)
        .ok_or_else(|| XyzError::InvalidElement(symbol.to_string()))?;

    let mut coord = || -> Result<f64, XyzError> {
        fields
            .next()
            .ok_or_else(|| XyzError::InvalidAtomLine(line.to_string()))?
            .parse()
            .map_err(|_| XyzError::InvalidAtomLine(line.to_string()))
    };
    let x = coord()?;
    let y = coord()?;
    let z = coord()?;
    // Everything after z (velocities, forces, per-atom charges in extended
    // XYZ) is deliberately not read -- see the module doc.

    Ok((Atom::new(element), Point3::new(x, y, z)))
}

/// Writes one XYZ frame.
///
/// A molecule with no 3D conformer is still written with zeros, matching
/// [`crate::io::sdf::write_sdf`]'s own "an ordinary molecule with no
/// coordinates is still written with zeros, as the caller expects." Unlike
/// SDF, though, XYZ has no all-zero-means-2D reading to fall back on — a
/// zero-filled write reads back claiming a (degenerate) conformer rather
/// than no geometry at all. Stated here rather than silently inherited from
/// the precedent it's borrowed from.
pub fn write_xyz(mol: &Molecule) -> String {
    let coords = mol.coords3();
    let mut out = format!("{}\n{}\n", mol.num_atoms(), mol.name().unwrap_or(""));
    for (i, atom) in mol.atoms().iter().enumerate() {
        let p = coords
            .and_then(|c| c.get(i))
            .copied()
            .unwrap_or(crate::core::geometry::Point3::new(0.0, 0.0, 0.0));
        out.push_str(&format!(
            "{} {:.6} {:.6} {:.6}\n",
            atom.element().symbol(),
            p.x,
            p.y,
            p.z
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::geometry::Point3;

    const WATER: &str = "3\nwater\nO 0.000000 0.000000 0.000000\nH 0.758602 0.000000 0.504284\nH 0.758602 0.000000 -0.504284\n";

    #[test]
    fn test_a_single_frame_round_trips() {
        let mol = parse_xyz(WATER).expect("valid XYZ");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.name(), Some("water"));
        assert_eq!(
            mol.coords3().expect("has coords")[0],
            Point3::new(0.0, 0.0, 0.0)
        );

        let written = write_xyz(&mol);
        let back = parse_xyz(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.coords3(), mol.coords3());
    }

    #[test]
    fn test_an_extended_xyz_comment_line_does_not_error() {
        let text =
            "1\nLattice=\"1 0 0 0 1 0 0 0 1\" Properties=species:S:1:pos:R:3\nC 0.0 0.0 0.0\n";
        let mol = parse_xyz(text).expect("an extended-style comment is still just a comment");
        assert_eq!(mol.num_atoms(), 1);
        // Deliberately not read into a unit cell -- see the module doc.
        assert!(mol.cell().is_none());
    }

    #[test]
    fn test_extra_per_atom_columns_are_ignored_not_misread() {
        let text = "1\nwith velocity\nC 0.0 0.0 0.0 1.5 0.0 0.0\n";
        let mol = parse_xyz(text).expect("valid XYZ");
        assert_eq!(mol.coords3().unwrap()[0], Point3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_an_unrecognised_element_is_a_clear_error() {
        let text = "1\n\nXx 0.0 0.0 0.0\n";
        let err = parse_xyz(text).unwrap_err();
        assert!(matches!(err, XyzError::InvalidElement(_)), "{err}");
    }

    #[test]
    fn test_too_few_atom_lines_is_a_clear_error_not_a_panic() {
        let text = "2\nonly one atom follows\nC 0.0 0.0 0.0\n";
        let err = parse_xyz(text).unwrap_err();
        assert!(matches!(err, XyzError::AtomCountMismatch { .. }), "{err}");
    }

    #[test]
    fn test_a_non_numeric_count_line_is_a_clear_error() {
        let err = parse_xyz("not a number\n\n").unwrap_err();
        assert!(matches!(err, XyzError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_malformed_coordinate_is_a_clear_error() {
        let text = "1\n\nC 0.0 0.0 not-a-number\n";
        let err = parse_xyz(text).unwrap_err();
        assert!(matches!(err, XyzError::InvalidAtomLine(_)), "{err}");
    }

    #[test]
    fn test_writing_a_molecule_with_no_conformer_falls_back_to_zeros() {
        // Matches write_sdf's own precedent, stated in this module's docs
        // rather than silently copied.
        let mut mol = Molecule::new();
        mol.add_atom(crate::core::atom::Atom::new(
            crate::core::atom::Element::carbon(),
        ));
        let written = write_xyz(&mol);
        assert!(
            written.contains("C 0.000000 0.000000 0.000000"),
            "{written}"
        );
    }

    #[test]
    fn test_all_zero_coordinates_are_not_read_as_a_conformer() {
        // A format with mandatory coordinate columns writes zeros for a
        // molecule that has none; believing them back leaves every atom at the
        // origin (#270).
        let text = "3\nethanol\nC 0.0 0.0 0.0\nC 0.0 0.0 0.0\nO 0.0 0.0 0.0\n";
        let mol = parse_xyz(text).expect("valid XYZ");
        assert_eq!(mol.num_atoms(), 3);
        assert!(
            !mol.has_coords3(),
            "all-zero is a placeholder, not a conformer"
        );
    }

    #[test]
    fn test_real_coordinates_are_still_read() {
        let text = "2\nreal\nC 0.0 0.0 0.0\nO 1.4 0.0 0.0\n";
        let mol = parse_xyz(text).expect("valid XYZ");
        assert!(mol.has_coords3());
    }

    /// PubChem CID 2519's 3D conformer.
    const CAFFEINE: &str = "24\ncaffeine\nO 0.4700 2.5688 0.0006\nO -3.1271 -0.4436 -0.0003\nN -0.9686 -1.3125 0.0000\nN 2.2182 0.1412 -0.0003\nN -1.3477 1.0797 -0.0001\nN 1.4119 -1.9372 0.0002\nC 0.8579 0.2592 -0.0008\nC 0.3897 -1.0264 -0.0004\nC 0.0307 1.4220 -0.0006\nC -1.9061 -0.2495 -0.0004\nC 2.5032 -1.1998 0.0003\nC -1.4276 -2.6960 0.0008\nC 3.1926 1.2061 0.0003\nC -2.2969 2.1881 0.0007\nH 3.5163 -1.5787 0.0008\nH -1.0451 -3.1973 -0.8937\nH -2.5186 -2.7596 0.0011\nH -1.0447 -3.1963 0.8957\nH 4.1992 0.7801 0.0002\nH 3.0468 1.8092 -0.8992\nH 3.0466 1.8083 0.9004\nH -1.8087 3.1651 -0.0003\nH -2.9322 2.1027 0.8881\nH -2.9346 2.1021 -0.8849\n";

    #[test]
    fn test_bonds_are_not_perceived_by_default() {
        let mol = parse_xyz_with_options(CAFFEINE, &XyzReadOptions::default()).unwrap();
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_perceived_caffeine_matches_its_stated_connectivity() {
        let options = XyzReadOptions {
            perceive_bonds: true,
            ..Default::default()
        };
        let mol = parse_xyz_with_options(CAFFEINE, &options).unwrap();
        let mut pairs: Vec<(usize, usize)> = mol
            .bonds()
            .iter()
            .map(|b| {
                let (a, c) = b.atoms();
                (a.min(c), a.max(c))
            })
            .collect();
        pairs.sort_unstable();
        // PubChem's own bond block for the same record, orders dropped.
        let expected = vec![
            (0, 8),
            (1, 9),
            (2, 7),
            (2, 9),
            (2, 11),
            (3, 6),
            (3, 10),
            (3, 12),
            (4, 8),
            (4, 9),
            (4, 13),
            (5, 7),
            (5, 10),
            (6, 7),
            (6, 8),
            (10, 14),
            (11, 15),
            (11, 16),
            (11, 17),
            (12, 18),
            (12, 19),
            (12, 20),
            (13, 21),
            (13, 22),
            (13, 23),
        ];
        assert_eq!(pairs, expected);
        // Bonded atoms get a hydrogen count, so the write isn't all brackets.
        assert!(mol.atoms().iter().all(|a| a.hydrogens().is_some()));
    }
}
