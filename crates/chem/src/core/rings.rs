//! Ring perception — finds the smallest set of smallest rings (SSSR).
//!
//! Structure layout places ring systems as polygons and hangs acyclic
//! substituents off them, so it needs to know which atoms and bonds form rings
//! before it can position anything.

use crate::core::molecule::Molecule;
use std::collections::{HashMap, HashSet, VecDeque};

/// A ring: its atoms in cycle order, and the bonds connecting them.
///
/// `atoms[i]` is bonded to `atoms[i + 1]`, and `atoms[len - 1]` back to
/// `atoms[0]`, so the sequence can be walked directly when placing a polygon.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
}

impl Ring {
    pub fn atoms(&self) -> &[usize] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[usize] {
        &self.bonds
    }

    /// Number of atoms in the ring — 6 for benzene.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    pub fn contains_atom(&self, atom_idx: usize) -> bool {
        self.atoms.contains(&atom_idx)
    }

    pub fn contains_bond(&self, bond_idx: usize) -> bool {
        self.bonds.contains(&bond_idx)
    }
}

/// How many independent rings a molecule has — its circuit rank,
/// `bonds - atoms + connected components`.
///
/// This is the exact size of the SSSR, so it doubles as the target when
/// selecting rings.
pub fn ring_count(molecule: &Molecule) -> usize {
    if molecule.num_atoms() == 0 {
        return 0;
    }
    let components = molecule.graph().connected_components().len();
    (molecule.num_bonds() + components).saturating_sub(molecule.num_atoms())
}

/// Finds the smallest set of smallest rings.
///
/// For each bond, the shortest path between its two atoms *without* using that
/// bond closes the smallest ring through it. Collecting those candidates,
/// taking them smallest-first, and keeping only the ones that are linearly
/// independent of those already chosen yields a smallest ring basis of exactly
/// [`ring_count`] rings.
///
/// Independence is checked over GF(2) in bond space: a ring is a set of bonds,
/// two rings combine by symmetric difference, and a ring that's the sum of
/// rings already chosen adds nothing. Without that check, a fused system like
/// naphthalene would also yield its 10-membered perimeter — a real cycle, but
/// just the sum of the two 6-rings, and laying it out as an extra ring would be
/// wrong.
pub fn find_sssr(molecule: &Molecule) -> Vec<Ring> {
    let target = ring_count(molecule);
    if target == 0 {
        return Vec::new();
    }

    let mut candidates: Vec<Ring> = Vec::new();
    let mut seen: HashSet<Vec<usize>> = HashSet::new();

    for bond_idx in 0..molecule.num_bonds() {
        let (a, b) = molecule.bond(bond_idx).atoms();
        let Some((atoms, mut bonds)) = shortest_path_excluding(molecule, a, b, bond_idx) else {
            // Removing the bond disconnects its atoms, so it's a bridge and
            // lies on no ring.
            continue;
        };
        bonds.push(bond_idx);

        // Key on the bond set so rotations and reversals of the same ring
        // collapse to one candidate.
        let mut key = bonds.clone();
        key.sort_unstable();
        if seen.insert(key) {
            candidates.push(Ring { atoms, bonds });
        }
    }

    candidates.sort_by_key(|r| r.len());

    let words = molecule.num_bonds().div_ceil(64).max(1);
    let mut basis: HashMap<usize, Vec<u64>> = HashMap::new();
    let mut sssr = Vec::with_capacity(target);

    for ring in candidates {
        if sssr.len() == target {
            break;
        }
        let vector = bond_vector(&ring, words);
        if insert_if_independent(&mut basis, vector) {
            sssr.push(ring);
        }
    }

    sssr
}

/// Finds the SSSR and records ring membership on the molecule's bonds.
///
/// `Bond::in_ring` exists but nothing previously set it, so it was always
/// false; this gives it a real value.
pub fn perceive_rings(molecule: &mut Molecule) -> Vec<Ring> {
    let rings = find_sssr(molecule);

    let mut in_ring = vec![false; molecule.num_bonds()];
    for ring in &rings {
        for &bond_idx in &ring.bonds {
            in_ring[bond_idx] = true;
        }
    }
    for (bond_idx, flag) in in_ring.into_iter().enumerate() {
        molecule.bond_mut(bond_idx).set_in_ring(flag);
    }

    rings
}

/// Shortest path from `from` to `to` without traversing `excluded_bond`,
/// as (atoms in order, bonds in order). `None` if they're disconnected without
/// it.
///
/// Bonds are recorded as traversed rather than looked up afterwards by atom
/// pair, so a molecule with two bonds between the same atoms resolves to the
/// one actually used.
fn shortest_path_excluding(
    molecule: &Molecule,
    from: usize,
    to: usize,
    excluded_bond: usize,
) -> Option<(Vec<usize>, Vec<usize>)> {
    let n = molecule.num_atoms();
    let mut visited = vec![false; n];
    let mut prev_atom = vec![usize::MAX; n];
    let mut prev_bond = vec![usize::MAX; n];
    let mut queue = VecDeque::new();

    visited[from] = true;
    queue.push_back(from);

    while let Some(current) = queue.pop_front() {
        if current == to {
            break;
        }
        for neighbor in molecule.neighbors(current) {
            if neighbor.bond_idx == excluded_bond || visited[neighbor.atom_idx] {
                continue;
            }
            visited[neighbor.atom_idx] = true;
            prev_atom[neighbor.atom_idx] = current;
            prev_bond[neighbor.atom_idx] = neighbor.bond_idx;
            queue.push_back(neighbor.atom_idx);
        }
    }

    if !visited[to] {
        return None;
    }

    let mut atoms = vec![to];
    let mut bonds = Vec::new();
    let mut current = to;
    while current != from {
        bonds.push(prev_bond[current]);
        current = prev_atom[current];
        atoms.push(current);
    }
    atoms.reverse();
    bonds.reverse();

    Some((atoms, bonds))
}

/// A ring as a bit per bond, for GF(2) arithmetic.
fn bond_vector(ring: &Ring, words: usize) -> Vec<u64> {
    let mut v = vec![0u64; words];
    for &bond_idx in &ring.bonds {
        v[bond_idx / 64] |= 1u64 << (bond_idx % 64);
    }
    v
}

fn highest_set_bit(v: &[u64]) -> Option<usize> {
    for (word_idx, word) in v.iter().enumerate().rev() {
        if *word != 0 {
            return Some(word_idx * 64 + (63 - word.leading_zeros() as usize));
        }
    }
    None
}

/// Reduces `vector` against the basis; if anything remains it's independent of
/// what's already there, so it joins the basis keyed by its leading bit.
fn insert_if_independent(basis: &mut HashMap<usize, Vec<u64>>, mut vector: Vec<u64>) -> bool {
    while let Some(lead) = highest_set_bit(&vector) {
        match basis.get(&lead) {
            Some(row) => {
                for (v, r) in vector.iter_mut().zip(row.iter()) {
                    *v ^= *r;
                }
            }
            None => {
                basis.insert(lead, vector);
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::prelude::*;

    /// Builds a molecule from a list of carbon-carbon bonds.
    fn carbon_skeleton(num_atoms: usize, bonds: &[(usize, usize)]) -> Molecule {
        let mut mol = Molecule::new();
        for _ in 0..num_atoms {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        for &(a, b) in bonds {
            mol.add_bond(Bond::new(a, b, BondOrder::Single)).unwrap();
        }
        mol
    }

    fn cyclohexane() -> Molecule {
        carbon_skeleton(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)])
    }

    /// Two six-membered rings sharing the 0-1 bond.
    fn naphthalene() -> Molecule {
        carbon_skeleton(
            10,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 0),
                (1, 6),
                (6, 7),
                (7, 8),
                (8, 9),
                (9, 0),
            ],
        )
    }

    #[test]
    fn test_acyclic_has_no_rings() {
        let propane = carbon_skeleton(3, &[(0, 1), (1, 2)]);
        assert_eq!(ring_count(&propane), 0);
        assert!(find_sssr(&propane).is_empty());
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new();
        assert_eq!(ring_count(&mol), 0);
        assert!(find_sssr(&mol).is_empty());
    }

    #[test]
    fn test_single_ring() {
        let mol = cyclohexane();
        assert_eq!(ring_count(&mol), 1);

        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
        assert_eq!(rings[0].bonds().len(), 6);
    }

    #[test]
    fn test_ring_atoms_are_in_cycle_order() {
        let rings = find_sssr(&cyclohexane());
        let atoms = rings[0].atoms();
        let mol = cyclohexane();

        // Consecutive atoms — including the wrap-around — must be bonded, which
        // is what lets a layout walk the ring as a polygon.
        for i in 0..atoms.len() {
            let a = atoms[i];
            let b = atoms[(i + 1) % atoms.len()];
            assert!(
                mol.graph().get_bond(a, b).is_some(),
                "atoms {a} and {b} are adjacent in the ring but not bonded"
            );
        }
    }

    #[test]
    fn test_ring_size_not_capped_at_seven() {
        // The previous DFS-based finder capped rings at 7 atoms, so
        // macrocycles were invisible to it.
        let macrocycle = carbon_skeleton(
            12,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 8),
                (8, 9),
                (9, 10),
                (10, 11),
                (11, 0),
            ],
        );
        let rings = find_sssr(&macrocycle);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 12);
    }

    #[test]
    fn test_fused_rings_exclude_the_perimeter() {
        let mol = naphthalene();
        // 11 bonds - 10 atoms + 1 component = 2.
        assert_eq!(ring_count(&mol), 2);

        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 2);
        // Both are six-membered: the 10-membered perimeter is a real cycle but
        // the sum of these two, so it's excluded as linearly dependent.
        for ring in &rings {
            assert_eq!(ring.len(), 6, "expected two 6-rings, got {}", ring.len());
        }
    }

    #[test]
    fn test_spiro_rings() {
        // Two rings sharing a single atom (atom 0) rather than a bond.
        let mol = carbon_skeleton(
            9,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 0),
                (0, 5),
                (5, 6),
                (6, 7),
                (7, 8),
                (8, 0),
            ],
        );
        assert_eq!(ring_count(&mol), 2);
        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 2);
        for ring in &rings {
            assert_eq!(ring.len(), 5);
        }
    }

    #[test]
    fn test_ring_with_substituent() {
        // Cyclohexane with a methyl: the bridging bond is on no ring.
        let mut mol = cyclohexane();
        let methyl = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(0, methyl, BondOrder::Single))
            .unwrap();

        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
        assert!(!rings[0].contains_atom(methyl));
    }

    #[test]
    fn test_disconnected_components() {
        // Two separate cyclohexanes: components count in the circuit rank, so
        // this is 2 rings rather than 1.
        let mut mol = cyclohexane();
        let base = mol.num_atoms();
        for _ in 0..6 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        for i in 0..6 {
            mol.add_bond(Bond::new(base + i, base + (i + 1) % 6, BondOrder::Single))
                .unwrap();
        }

        assert_eq!(ring_count(&mol), 2);
        assert_eq!(find_sssr(&mol).len(), 2);
    }

    #[test]
    fn test_perceive_rings_sets_in_ring_flags() {
        let mut mol = cyclohexane();
        let methyl = mol.add_atom(Atom::new(Element::carbon()));
        let bridge = mol
            .add_bond(Bond::new(0, methyl, BondOrder::Single))
            .unwrap();

        // Nothing sets this before perception runs.
        assert!(mol.bonds().iter().all(|b| !b.in_ring()));

        let rings = perceive_rings(&mut mol);
        assert_eq!(rings.len(), 1);

        for bond_idx in 0..mol.num_bonds() {
            let expected = bond_idx != bridge;
            assert_eq!(
                mol.bond(bond_idx).in_ring(),
                expected,
                "bond {bond_idx} in_ring flag"
            );
        }
    }

    #[test]
    fn test_contains_bond() {
        let mol = cyclohexane();
        let rings = find_sssr(&mol);
        for bond_idx in 0..mol.num_bonds() {
            assert!(rings[0].contains_bond(bond_idx));
        }
    }

    #[test]
    fn test_bridged_bicyclic() {
        // Norbornane skeleton: bridgeheads 0 and 3 joined by three chains, of
        // two, two and one intermediate atoms.
        let mol = carbon_skeleton(
            7,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (0, 4),
                (4, 5),
                (5, 3),
                (0, 6),
                (6, 3),
            ],
        );
        // 8 bonds - 7 atoms + 1 = 2.
        assert_eq!(ring_count(&mol), 2);

        let rings = find_sssr(&mol);
        assert_eq!(rings.len(), 2);

        // Pairing the chains gives three cycles, of 6, 5 and 5 atoms. The two
        // 5-rings are independent of each other — their symmetric difference
        // is the 6-ring — so the basis is both of them, and the 6-ring is left
        // out as their sum. This is the textbook SSSR for norbornane, and the
        // reason the smallest-first ordering matters: taking the 6-ring first
        // would give an equally valid ring basis but not the *smallest* one.
        let mut sizes: Vec<usize> = rings.iter().map(|r| r.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![5, 5]);
    }
}
