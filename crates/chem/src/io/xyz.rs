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

use crate::core::molecule::Molecule;
use crate::io::errors::XyzError;

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

    mol.set_coords3(coords)
        .map_err(|e| XyzError::ParseError(e.to_string()))?;
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
}
