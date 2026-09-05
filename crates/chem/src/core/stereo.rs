//! Reading double-bond geometry out of a drawing, and back into
//! [`BondStereo`].
//!
//! A coordinate format does not state a configuration; it states positions, and
//! the configuration follows from them. So a reader that takes such a file
//! literally hands back a molecule that is chemically cis while claiming
//! nothing — the same gap [`crate::io::aromaticity::detect_aromaticity`] fills
//! for aromatic rings, and filled for the same reason.

use crate::core::geometry::Point2;
use crate::core::prelude::*;

/// Sets [`BondStereo`] on every double bond whose drawn geometry determines
/// one.
///
/// Call it after coordinates are in place. Bonds whose geometry says nothing
/// are left alone rather than guessed at: a double bond with substituents on
/// only one end has no configuration to read, and asserting `E` there would
/// invent stereochemistry the file never claimed.
///
/// A molecule with no 2D layout is left entirely untouched.
pub fn perceive_bond_stereo(molecule: &mut Molecule) {
    let Some(points) = molecule.coords() else {
        return;
    };
    let points = points.to_vec();

    let mut updates: Vec<(usize, BondStereo)> = Vec::new();
    for bond_idx in 0..molecule.num_bonds() {
        if molecule.bond(bond_idx).order() != BondOrder::Double {
            continue;
        }
        // The file said "either" explicitly. Its coordinates put the
        // substituents somewhere, as any drawing must, but that placement is
        // not a claim and must not be read as one.
        if molecule.bond(bond_idx).stereo() == BondStereo::Unspecified {
            continue;
        }
        let (left, right) = (
            molecule.bond(bond_idx).atom1(),
            molecule.bond(bond_idx).atom2(),
        );
        if let Some(stereo) = geometry_of(molecule, &points, left, right) {
            updates.push((bond_idx, stereo));
        }
    }

    for (bond_idx, stereo) in updates {
        let bond = molecule.bond(bond_idx);
        let replacement = Bond::new(bond.atom1(), bond.atom2(), bond.order())
            .with_aromatic(bond.is_aromatic())
            .with_stereo(stereo);
        *molecule.bond_mut(bond_idx) = replacement;
    }
}

/// Which side of the double-bond axis each end's reference substituent sits,
/// and therefore the configuration — or `None` when the drawing does not say.
///
/// The reference substituent is the lowest-indexed neighbour that is not the
/// other end of the double bond, matching the choice `write_sdf` makes for
/// wedges so the two cannot disagree about which atom they are describing.
pub(crate) fn geometry_of(
    molecule: &Molecule,
    points: &[Point2],
    left: usize,
    right: usize,
) -> Option<BondStereo> {
    let left_ref = reference_substituent(molecule, left, right)?;
    let right_ref = reference_substituent(molecule, right, left)?;

    let axis = points[right] - points[left];
    let left_side = side_of(axis, points[left_ref] - points[left]);
    let right_side = side_of(axis, points[right_ref] - points[right]);

    match (left_side, right_side) {
        // Same side of the axis is cis; opposite sides is trans. Anything
        // collinear leaves the drawing unable to say, which is not the same as
        // saying trans.
        (Some(a), Some(b)) if a == b => Some(BondStereo::Z),
        (Some(_), Some(_)) => Some(BondStereo::E),
        _ => None,
    }
}

/// The neighbour whose position speaks for this end of the double bond.
fn reference_substituent(molecule: &Molecule, atom: usize, partner: usize) -> Option<usize> {
    molecule
        .neighbors(atom)
        .iter()
        .map(|n| n.atom_idx)
        .filter(|&n| n != partner)
        .min()
}

/// Which side of `axis` a substituent falls, by the sign of the cross product.
///
/// The same test `write_sdf`'s wedge placement uses. `None` when the two are
/// collinear, where there is no side to be on.
fn side_of(axis: Point2, to_substituent: Point2) -> Option<bool> {
    let cross = axis.x * to_substituent.y - axis.y * to_substituent.x;
    if cross.abs() < 1e-9 {
        None
    } else {
        Some(cross > 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::layout::layout;
    use crate::io::smiles::parse_smiles;

    /// The configuration a laid-out molecule's drawing actually shows.
    fn drawn(smiles: &str) -> BondStereo {
        let mut mol = parse_smiles(smiles).expect(smiles);
        layout(&mut mol);
        let points = mol.coords().expect("laid out").to_vec();
        let (idx, bond) = mol
            .bonds()
            .iter()
            .enumerate()
            .find(|(_, b)| b.order() == BondOrder::Double)
            .expect("a double bond");
        let _ = idx;
        geometry_of(&mol, &points, bond.atom1(), bond.atom2()).unwrap_or(BondStereo::None)
    }

    #[test]
    fn test_the_layout_draws_cis_and_trans_differently() {
        // The bug: every double bond came out drawn trans, so a cis molecule
        // was written as a valid file asserting the wrong isomer — and drawn
        // as the wrong isomer on screen (#198).
        assert_eq!(drawn("F/C=C/F"), BondStereo::E);
        assert_eq!(drawn("F/C=C\\F"), BondStereo::Z);
    }

    #[test]
    fn test_the_reversed_spelling_draws_the_same_isomer() {
        // `C(/F)=C/F` is cis despite its two forward slashes.
        assert_eq!(drawn("C(/F)=C/F"), BondStereo::Z);
        assert_eq!(drawn("C(\\F)=C/F"), BondStereo::E);
    }

    #[test]
    fn test_a_longer_chain_keeps_its_geometry() {
        assert_eq!(drawn("CC/C=C/CC"), BondStereo::E);
        assert_eq!(drawn("CC/C=C\\CC"), BondStereo::Z);
    }

    #[test]
    fn test_perception_recovers_what_the_layout_drew() {
        for (smiles, expected) in [
            ("F/C=C/F", BondStereo::E),
            ("F/C=C\\F", BondStereo::Z),
            ("CC/C=C\\CC", BondStereo::Z),
        ] {
            let mut mol = parse_smiles(smiles).expect(smiles);
            layout(&mut mol);
            // Clear what the parser knew, so only the drawing can answer.
            for i in 0..mol.num_bonds() {
                let b = mol.bond(i);
                *mol.bond_mut(i) = Bond::new(b.atom1(), b.atom2(), b.order());
            }
            perceive_bond_stereo(&mut mol);
            let found = mol
                .bonds()
                .iter()
                .find(|b| b.order() == BondOrder::Double)
                .expect("a double bond")
                .stereo();
            assert_eq!(found, expected, "{smiles}");
        }
    }

    #[test]
    fn test_an_undetermined_double_bond_is_left_alone() {
        // A double bond with a substituent on one end only has no
        // configuration to read. Whatever the layout happened to draw, saying
        // `E` here would invent stereochemistry the molecule never claimed.
        let mut mol = parse_smiles("CC=C").expect("valid SMILES");
        layout(&mut mol);
        perceive_bond_stereo(&mut mol);
        let stereo = mol
            .bonds()
            .iter()
            .find(|b| b.order() == BondOrder::Double)
            .expect("a double bond")
            .stereo();
        assert_eq!(stereo, BondStereo::None);
    }

    #[test]
    fn test_an_undetermined_double_bond_is_marked_either_in_the_file() {
        // Any 2D drawing puts the substituents on *some* side, so writing
        // coordinates alone would let a reader — ours included — derive a
        // configuration the molecule never claimed. Field 3 says "either",
        // which is the file declining to answer rather than answering.
        let mol = parse_smiles("FC=CF").expect("valid SMILES");
        let text = crate::io::sdf::write_sdf(&mol);
        assert!(
            text.lines().any(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                f.len() == 7 && f[2] == "2" && f[3] == "3"
            }),
            "no `either` marker written:\n{text}"
        );

        // And the round trip must not invent one.
        let back = crate::io::sdf::parse_sdf(&text).expect("round trips");
        let stereo = back
            .bonds()
            .iter()
            .find(|b| b.order() == BondOrder::Double)
            .expect("a double bond")
            .stereo();
        assert_ne!(stereo, BondStereo::E);
        assert_ne!(stereo, BondStereo::Z);
    }

    #[test]
    fn test_a_stated_configuration_is_not_marked_either() {
        let mol = parse_smiles("F/C=C\\F").expect("valid SMILES");
        let text = crate::io::sdf::write_sdf(&mol);
        assert!(
            !text.lines().any(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                f.len() == 7 && f[2] == "2" && f[3] == "3"
            }),
            "a stated configuration was written as `either`:\n{text}"
        );
        let back = crate::io::sdf::parse_sdf(&text).expect("round trips");
        assert_eq!(
            back.bonds()
                .iter()
                .find(|b| b.order() == BondOrder::Double)
                .expect("a double bond")
                .stereo(),
            BondStereo::Z
        );
    }

    #[test]
    fn test_a_molecule_with_no_drawing_is_untouched() {
        let mut mol = parse_smiles("F/C=C\\F").expect("valid SMILES");
        assert!(mol.coords().is_none());
        let before: Vec<BondStereo> = mol.bonds().iter().map(|b| b.stereo()).collect();
        perceive_bond_stereo(&mut mol);
        let after: Vec<BondStereo> = mol.bonds().iter().map(|b| b.stereo()).collect();
        assert_eq!(before, after);
    }

    #[test]
    fn test_a_ring_double_bond_still_closes_its_ring() {
        // A double bond inside a ring has its geometry forced by the ring, so
        // the corrective pass must leave it alone — reflecting one side would
        // tear the ring open.
        let mut mol = parse_smiles("C1=CCCCC1").expect("valid SMILES");
        layout(&mut mol);
        let points = mol.coords().expect("laid out");
        let closing = mol
            .bonds()
            .iter()
            .map(|b| (points[b.atom1()] - points[b.atom2()]).length())
            .fold(0.0_f64, f64::max);
        assert!(closing < 3.0, "the ring was pulled apart: {closing}");
    }
}
