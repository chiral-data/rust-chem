//! Bonds from interatomic distances, for structures whose file states none
//! (#405).

use std::collections::HashMap;

use crate::core::bond::{Bond, BondOrder};
use crate::core::molecule::Molecule;

/// Added to r₁ + r₂ when deciding whether two atoms are bonded; OpenBabel's
/// and RDKit's value.
const TOLERANCE: f64 = 0.4;

/// Closer than this is two overlapping atoms, a modelling artefact, not a bond.
const MIN_DISTANCE: f64 = 0.4;

/// Adds a single bond between every pair of atoms closer than the sum of their
/// covalent radii plus 0.4 Å, returning how many it added.
///
/// Only adds: a pair already bonded is left with its stated order, and nothing
/// is added at all without [`Molecule::coords3`]. Skipped pairs are those with
/// an element lacking a covalent radius ([`crate::core::atom::Element::UNKNOWN`]
/// included), H–H, two different non-blank `alt_loc`s, and anything under
/// 0.4 Å. A hydrogen bonds at most once, to its nearest partner, and not at all
/// if it already has a bond.
///
/// Every added bond is [`BondOrder::Single`]: geometry gives connectivity, not
/// order, so a perceived benzene ring is six single bonds and not aromatic.
///
/// Runs in roughly linear time by binning atoms on a uniform grid.
pub fn perceive_bonds(mol: &mut Molecule) -> usize {
    let Some(coords) = mol.coords3() else {
        return 0;
    };
    let n = mol.num_atoms();
    let radii: Vec<Option<f64>> = mol
        .atoms()
        .iter()
        .zip(coords)
        .map(|(atom, p)| {
            let finite = p.x.is_finite() && p.y.is_finite() && p.z.is_finite();
            atom.element().covalent_radius().filter(|_| finite)
        })
        .collect();
    let Some(max_radius) = radii.iter().flatten().copied().reduce(f64::max) else {
        return 0;
    };
    let is_h: Vec<bool> = mol
        .atoms()
        .iter()
        .map(|a| a.element().atomic_number == 1)
        .collect();
    let alt: Vec<Option<char>> = (0..n)
        .map(|i| {
            mol.site(i)
                .and_then(|s| s.alt_loc)
                .filter(|c| !c.is_whitespace())
        })
        .collect();

    // No bond can be longer than the cell edge, so every partner lies in the
    // 27 cells around an atom.
    let cell = 2.0 * max_radius + TOLERANCE;
    let key = |i: usize| {
        let p = coords[i];
        (
            (p.x / cell).floor() as i32,
            (p.y / cell).floor() as i32,
            (p.z / cell).floor() as i32,
        )
    };
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    for i in (0..n).filter(|&i| radii[i].is_some()) {
        grid.entry(key(i)).or_default().push(i);
    }

    let mut candidates: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        let Some(ri) = radii[i] else { continue };
        let (cx, cy, cz) = key(i);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbour = (
                        cx.saturating_add(dx),
                        cy.saturating_add(dy),
                        cz.saturating_add(dz),
                    );
                    let Some(bin) = grid.get(&neighbour) else {
                        continue;
                    };
                    for &j in bin.iter().filter(|&&j| j > i) {
                        if is_h[i] && is_h[j] {
                            continue;
                        }
                        if let (Some(a), Some(b)) = (alt[i], alt[j])
                            && a != b
                        {
                            continue;
                        }
                        let rj = radii[j].expect("only atoms with a radius are binned");
                        let d = coords[i].distance(coords[j]);
                        if d >= MIN_DISTANCE && d < ri + rj + TOLERANCE {
                            candidates.push((i, j, d));
                        }
                    }
                }
            }
        }
    }

    // Each unbonded hydrogen keeps only its nearest candidate, ties to the
    // lower index so the result doesn't depend on grid order.
    let mut nearest: Vec<Option<(f64, usize)>> = vec![None; n];
    for &(i, j, d) in &candidates {
        for (h, other) in [(i, j), (j, i)] {
            if is_h[h] && nearest[h].is_none_or(|best| (d, other) < best) {
                nearest[h] = Some((d, other));
            }
        }
    }
    let keeps = |h: usize, other: usize| {
        !is_h[h] || (mol.degree(h) == 0 && nearest[h].is_some_and(|(_, o)| o == other))
    };

    let mut accepted: Vec<(usize, usize)> = candidates
        .into_iter()
        .filter(|&(i, j, _)| keeps(i, j) && keeps(j, i) && !mol.graph().has_edge(i, j))
        .map(|(i, j, _)| (i, j))
        .collect();
    accepted.sort_unstable();

    for &(i, j) in &accepted {
        mol.add_bond(Bond::new(i, j, BondOrder::Single))
            .expect("indices come from this molecule");
    }
    accepted.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{Atom, Element};
    use crate::core::geometry::Point3;
    use crate::core::site::AtomSite;

    fn molecule(atoms: &[(u8, [f64; 3])]) -> Molecule {
        let mut mol = Molecule::new();
        for &(z, _) in atoms {
            mol.add_atom(Atom::new(Element::new(z).unwrap()));
        }
        mol.set_coords3(
            atoms
                .iter()
                .map(|&(_, [x, y, z])| Point3::new(x, y, z))
                .collect(),
        )
        .unwrap();
        mol
    }

    fn pairs(mol: &Molecule) -> Vec<(usize, usize)> {
        let mut out: Vec<_> = mol
            .bonds()
            .iter()
            .map(|b| {
                let (a, c) = b.atoms();
                (a.min(c), a.max(c))
            })
            .collect();
        out.sort_unstable();
        out
    }

    fn water(offset: f64) -> Vec<(u8, [f64; 3])> {
        vec![
            (8, [offset, 0.0, 0.0]),
            (1, [offset + 0.757, 0.586, 0.0]),
            (1, [offset - 0.757, 0.586, 0.0]),
        ]
    }

    #[test]
    fn test_water_gets_two_bonds() {
        let mut mol = molecule(&water(0.0));
        assert_eq!(perceive_bonds(&mut mol), 2);
        assert_eq!(pairs(&mol), vec![(0, 1), (0, 2)]);
        assert!(mol.bonds().iter().all(|b| b.order() == BondOrder::Single));
    }

    #[test]
    fn test_ethane_gets_seven_bonds_and_no_h_h() {
        let mut mol = molecule(&[
            (6, [0.0, 0.0, 0.765]),
            (6, [0.0, 0.0, -0.765]),
            (1, [1.018, 0.0, 1.162]),
            (1, [-0.509, 0.882, 1.162]),
            (1, [-0.509, -0.882, 1.162]),
            (1, [-1.018, 0.0, -1.162]),
            (1, [0.509, -0.882, -1.162]),
            (1, [0.509, 0.882, -1.162]),
        ]);
        assert_eq!(perceive_bonds(&mut mol), 7);
        assert!(mol.bonds().iter().all(|b| {
            let (a, c) = b.atoms();
            mol.atom(a).element().atomic_number != 1 || mol.atom(c).element().atomic_number != 1
        }));
    }

    #[test]
    fn test_two_waters_three_angstroms_apart_stay_separate() {
        let mut atoms = water(0.0);
        atoms.extend(water(3.0));
        let mut mol = molecule(&atoms);
        assert_eq!(perceive_bonds(&mut mol), 4);
        assert_eq!(pairs(&mol), vec![(0, 1), (0, 2), (3, 4), (3, 5)]);
    }

    #[test]
    fn test_a_hydrogen_bonds_only_to_its_nearest_partner() {
        // H sits within bonding range of both oxygens, nearer the first.
        let mut mol = molecule(&[
            (8, [0.0, 0.0, 0.0]),
            (1, [1.0, 0.0, 0.0]),
            (8, [2.2, 0.0, 0.0]),
        ]);
        perceive_bonds(&mut mol);
        assert_eq!(pairs(&mol), vec![(0, 1)]);
    }

    #[test]
    fn test_differing_alt_locs_are_not_bonded() {
        let mut mol = molecule(&[
            (6, [0.0, 0.0, 0.0]),
            (6, [1.5, 0.0, 0.0]),
            (6, [3.0, 0.0, 0.0]),
        ]);
        let site = |alt: Option<char>| AtomSite {
            alt_loc: alt,
            ..AtomSite::default()
        };
        mol.set_sites(vec![site(Some('A')), site(Some('B')), site(Some(' '))])
            .unwrap();
        perceive_bonds(&mut mol);
        assert_eq!(pairs(&mol), vec![(1, 2)]);
    }

    #[test]
    fn test_an_element_with_no_radius_is_left_alone() {
        let mut mol = molecule(&[(6, [0.0, 0.0, 0.0]), (6, [1.5, 0.0, 0.0])]);
        mol.atoms_mut()[1] = Atom::new(Element::UNKNOWN);
        assert_eq!(perceive_bonds(&mut mol), 0);
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_overlapping_atoms_are_not_bonded() {
        let mut mol = molecule(&[(6, [0.0, 0.0, 0.0]), (6, [0.2, 0.0, 0.0])]);
        assert_eq!(perceive_bonds(&mut mol), 0);
    }

    #[test]
    fn test_a_stated_double_bond_keeps_its_order_and_is_not_duplicated() {
        let mut mol = molecule(&[(6, [0.0, 0.0, 0.0]), (8, [1.21, 0.0, 0.0])]);
        mol.add_bond(Bond::new(1, 0, BondOrder::Double)).unwrap();
        assert_eq!(perceive_bonds(&mut mol), 0);
        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bond(0).order(), BondOrder::Double);
    }

    #[test]
    fn test_a_hydrogen_already_bonded_gains_nothing() {
        let mut mol = molecule(&[
            (8, [0.0, 0.0, 0.0]),
            (1, [1.0, 0.0, 0.0]),
            (8, [2.2, 0.0, 0.0]),
        ]);
        mol.add_bond(Bond::new(1, 2, BondOrder::Single)).unwrap();
        assert_eq!(perceive_bonds(&mut mol), 0);
    }

    #[test]
    fn test_no_coordinates_adds_nothing() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::carbon()));
        assert_eq!(perceive_bonds(&mut mol), 0);
    }

    /// Deterministic xorshift, so the cross-check needs no `rand`.
    fn random_points(n: usize, extent: f64, mut seed: u64) -> Vec<[f64; 3]> {
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n)
            .map(|_| {
                [
                    next() * extent - extent / 2.0,
                    next() * extent,
                    next() * extent,
                ]
            })
            .collect()
    }

    #[test]
    fn test_binning_matches_an_all_pairs_reference() {
        // Heavy atoms only, so the hydrogen cap doesn't hide a missed pair.
        let elements = [6u8, 7, 8, 16, 26, 55];
        let atoms: Vec<(u8, [f64; 3])> = random_points(400, 14.0, 0x9e37_79b9_7f4a_7c15)
            .into_iter()
            .enumerate()
            .map(|(i, p)| (elements[i % elements.len()], p))
            .collect();
        let mut mol = molecule(&atoms);
        perceive_bonds(&mut mol);

        let mut expected = Vec::new();
        for i in 0..atoms.len() {
            for j in i + 1..atoms.len() {
                let r = |k: usize| Element::new(atoms[k].0).unwrap().covalent_radius().unwrap();
                let [a, b] = [atoms[i].1, atoms[j].1].map(|[x, y, z]| Point3::new(x, y, z));
                let d = a.distance(b);
                if d >= MIN_DISTANCE && d < r(i) + r(j) + TOLERANCE {
                    expected.push((i, j));
                }
            }
        }
        assert!(expected.len() > 50, "the sample should exercise many pairs");
        assert_eq!(pairs(&mol), expected);
    }

    #[test]
    #[ignore = "timing; run with --release --ignored"]
    fn test_a_100k_atom_grid_is_fast() {
        // A simple-cubic lattice of carbons 1.5 Å apart: the six axial
        // neighbours bond, the 2.12 Å diagonals don't. Only perception is timed.
        let side = 47; // 47³ ≈ 103k atoms
        let atoms: Vec<(u8, [f64; 3])> = (0..side * side * side)
            .map(|k| {
                let (x, y, z) = (k % side, (k / side) % side, k / (side * side));
                (6, [x as f64 * 1.5, y as f64 * 1.5, z as f64 * 1.5])
            })
            .collect();
        let mut mol = molecule(&atoms);
        let start = std::time::Instant::now();
        let added = perceive_bonds(&mut mol);
        let elapsed = start.elapsed();
        eprintln!("{} atoms, {added} bonds, {elapsed:?}", atoms.len());
        assert_eq!(added, 3 * side * side * (side - 1));
        assert!(elapsed.as_secs_f64() < 1.0, "{elapsed:?}");
    }
}
