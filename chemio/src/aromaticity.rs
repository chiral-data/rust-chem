use chemcore::{bond::BondOrder, molecule::Molecule};

const MIN_RING_SIZE: usize = 3;
const MAX_RING_SIZE: usize = 7;

/// Detects aromatic rings in the given molecule and marks
/// both the atoms and bonds belonging to aromatic rings.
///
/// # Parameters
/// - `mol`: A mutable reference to the molecule in which aromaticity should be detected.
///
/// # Behavior
/// - Finds all rings in the molecule.
/// - Evaluates each ring with a simplified aromaticity heuristic.
/// - Marks aromatic rings by setting aromatic flags on the atoms and bonds.
///
/// # Notes
/// - This uses a lightweight ring-search and aromaticity heuristic.
/// - It is NOT a full, industry-grade aromaticity model.
pub fn detect_aromaticity(mol: &mut Molecule) {
    let rings = find_rings(mol);

    for ring in rings {
        if is_aromatic_ring(mol, &ring) {
            mark_ring_aromatic(mol, &ring);
        }
    }
}

/// Finds all unique rings in the molecule.
///
/// # Parameters
/// - `mol`: The molecule to search for rings.
///
/// # Returns
/// A `Vec<Vec<usize>>` where each inner vector is a sequence of atom indices
/// representing a detected ring.
///
/// # Notes
/// - Calls `find_rings_from` starting at every atom.
/// - Rings are deduplicated and sorted to remove duplicates (e.g. rotated versions).
fn find_rings(mol: &Molecule) -> Vec<Vec<usize>> {
    let mut rings = Vec::new();

    for start in 0..mol.num_atoms() {
        let paths = find_rings_from(mol, start, MIN_RING_SIZE, MAX_RING_SIZE);
        rings.extend(paths);
    }

    rings.sort();
    rings.dedup();
    rings
}

/// Attempts to detect all rings reachable from a given start atom.
///
/// # Parameters
/// - `mol`: Molecule being analyzed.
/// - `start`: Atom index at which DFS ring detection begins.
/// - `min_size`: Minimum allowed ring size.
/// - `max_size`: Maximum allowed ring size.
///
/// # Returns
/// A vector of atom-index sequences, each representing a ring found starting from `start`.
///
/// # Behavior
/// - Performs DFS via `dfs_rings`.
/// - Filters out rings that do not satisfy the size constraints.
fn find_rings_from(
    mol: &Molecule,
    start: usize,
    min_size: usize,
    max_size: usize,
) -> Vec<Vec<usize>> {
    let mut rings = Vec::new();
    let mut path = vec![start];
    let mut visited = vec![false; mol.num_atoms()];
    visited[start] = true;

    dfs_rings(
        mol,
        start,
        start,
        &mut path,
        &mut visited,
        &mut rings,
        max_size,
    );

    // Filter out rings that are too small or too large
    rings
        .into_iter()
        .filter(|r| r.len() > min_size && r.len() <= max_size)
        .collect()
}

/// Depth-first search for cyclic paths (rings).
///
/// # Parameters
/// - `mol`: Molecule graph to traverse.
/// - `current`: Current atom index in DFS traversal.
/// - `start`: Atom index where the DFS began (potential ring closure point).
/// - `path`: Mutable list of visited atoms in the current DFS path.
/// - `visited`: Boolean flag vector preventing revisiting atoms in the same DFS branch.
/// - `rings`: Accumulator for detected rings.
/// - `max_size`: Maximum allowed ring size.
///
/// # Behavior
/// - Explores neighbors recursively.
/// - If a neighbor equals `start`, and path length ≥ 3, a ring is detected.
/// - Avoids infinite recursion by marking atoms as visited only during the path.
fn dfs_rings(
    mol: &Molecule,
    current: usize,
    start: usize,
    path: &mut Vec<usize>,
    visited: &mut [bool],
    rings: &mut Vec<Vec<usize>>,
    max_size: usize,
) {
    if path.len() > max_size {
        return;
    }

    for neighbour in mol.neighbors(current) {
        let next = neighbour.atom_idx;

        // If current_atom == atom we started with |> We found a Ring
        if next == start && path.len() >= 3 {
            rings.push(path.clone());
            continue;
        }

        if !visited[next] && path.len() <= max_size {
            visited[next] = true;
            path.push(next);
            dfs_rings(mol, next, start, path, visited, rings, max_size);
            path.pop();
            visited[next] = false;
        }
    }
}

/// Determines whether a given ring satisfies a simplified aromaticity model.
///
/// # Parameters
/// - `mol`: The molecule containing the ring.
/// - `ring`: Slice of atom indices representing the ring.
///
/// # Returns
/// `true` if the ring is considered aromatic; `false` otherwise.
///
/// # Aromaticity Model Used
/// - Ring must be 5 or 6 members long.
/// - Atoms must be C, N, O, S, P (atomic numbers: 5,6,7,8,15,16).
/// - π-electron count must satisfy Hückel's rule: `4n + 2` (6 or 10 used here).
///
/// # Notes
/// This is a **heuristic** and not a fully correct aromaticity method.
fn is_aromatic_ring(mol: &Molecule, ring: &[usize]) -> bool {
    // INFO:Simple heuristic: 5 or 6 membered rings with alternating bonds
    // TODO:Is this correct?? Verify
    if ring.len() != 5 && ring.len() != 6 {
        return false;
    }

    // Check if all atoms are C, N, O, S
    for &atom_idx in ring {
        let atomic_num = mol.atom(atom_idx).atomic_number();
        // FIX:Hardcoded - Generic approach later
        if ![5, 6, 7, 8, 15, 16].contains(&atomic_num) {
            return false;
        }
    }

    // Count pi electrons (simplified Hückel rule)
    let mut pi_electrons = 0;
    for &atom_idx in ring {
        let atom = mol.atom(atom_idx);
        let atomic_num = atom.atomic_number();

        match atomic_num {
            // A ring carbon only has a p-orbital to contribute to a
            // conjugated pi system if it's actually part of an in-ring
            // double (or already-aromatic) bond -- a fully saturated
            // carbon (e.g. every carbon in cyclohexane) has none, and its
            // presence breaks conjugation around the whole ring, so the
            // ring can't be aromatic regardless of what the other atoms
            // are. Previously this contributed 1 pi-electron
            // unconditionally, so any 6-membered all-carbon ring summed to
            // exactly 6 and passed the Hückel check whether or not it had
            // any unsaturation at all.
            6 => {
                if !has_in_ring_pi_bond(mol, atom_idx, ring) {
                    return false;
                }
                pi_electrons += 1;
            }
            7 => {
                // N: check if it has lone pair or is sp2
                if mol.degree(atom_idx) == 2 {
                    pi_electrons += 2; // Pyridine-like
                } else {
                    pi_electrons += 1; // Pyrrole-like
                }
            }
            // O, S with lone pair
            8 | 16 if mol.degree(atom_idx) == 2 => {
                pi_electrons += 2;
            }
            _ => {}
        }
    }

    // Hückel's rule: 4n + 2 pi electrons
    pi_electrons == 6 || pi_electrons == 10
}

/// Whether `atom_idx` has a double or aromatic bond to another atom that's
/// also in `ring` -- i.e. an in-ring pi bond, not just any double bond
/// (an exocyclic one, like a ring carbon's C=O in a cyclohexanone, doesn't
/// make the ring's own bonds conjugated).
fn has_in_ring_pi_bond(mol: &Molecule, atom_idx: usize, ring: &[usize]) -> bool {
    mol.neighbors(atom_idx).iter().any(|neighbor| {
        ring.contains(&neighbor.atom_idx)
            && matches!(
                mol.bond(neighbor.bond_idx).order(),
                BondOrder::Double | BondOrder::Aromatic
            )
    })
}

/// Marks the atoms and bonds in a ring as aromatic.
///
/// # Parameters
/// - `mol`: Mutable reference to the molecule.
/// - `ring`: List of atom indices representing a ring.
///
/// # Behavior
/// - Sets aromatic flags on each atom in the ring.
/// - Sets aromatic flags on each bond in the ring.
/// - Bond order is explicitly set to `BondOrder::Aromatic`.
fn mark_ring_aromatic(mol: &mut Molecule, ring: &[usize]) {
    // Mark atoms -> aromatic
    for &atom_idx in ring {
        mol.atom_mut(atom_idx).set_aromatic(true);
    }

    // Mark bonds -> aromatic
    for i in 0..ring.len() {
        let next = (i + 1) % ring.len(); // index of new Atom
        if let Some(bond_idx) = mol.graph().get_bond(ring[i], ring[next]) {
            mol.bond_mut(bond_idx).set_aromatic(true);
            mol.bond_mut(bond_idx).set_order(BondOrder::Aromatic);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_benzene_aromaticity() {
        let mut mol = crate::smiles::parse_smiles("C1=CC=CC=C1").unwrap();
        detect_aromaticity(&mut mol);

        for i in 0..6 {
            assert!(mol.atom(i).is_aromatic());
        }
    }

    #[test]
    fn test_cyclopentane_not_aromatic() {
        let mut mol = crate::smiles::parse_smiles("C1CCCC1").unwrap();
        detect_aromaticity(&mut mol);

        assert!(!mol.atom(0).is_aromatic());
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        // Fully saturated 6-membered all-carbon ring, no double bonds at
        // all -- previously flagged aromatic because every carbon
        // contributed a pi-electron unconditionally, landing on exactly 6
        // (the same count real benzene gets) regardless of saturation.
        let mut mol = crate::smiles::parse_smiles("C1CCCCC1").unwrap();
        detect_aromaticity(&mut mol);

        for i in 0..6 {
            assert!(!mol.atom(i).is_aromatic(), "atom {} wrongly aromatic", i);
        }
    }

    #[test]
    fn test_cyclohexanone_not_aromatic() {
        // A saturated ring carbon with an *exocyclic* double bond (the
        // ketone's C=O) shouldn't count as an in-ring pi bond -- the ring's
        // own bonds are still all single.
        let mut mol = crate::smiles::parse_smiles("O=C1CCCCC1").unwrap();
        detect_aromaticity(&mut mol);

        for i in 0..mol.num_atoms() {
            assert!(!mol.atom(i).is_aromatic(), "atom {} wrongly aromatic", i);
        }
    }

    #[test]
    fn test_complex_polycyclic_aromaticity() {
        // Naphthalene: Two fused benzene rings
        // Two aromatic 6-membered rings sharing edge 4-5
        // Ring 1: [0, 1, 2, 3, 4, 5, 9, 8, 7, 6] (outer perimeter - NOT a simple ring)
        // Ring 2: [0, 1, 2, 3, 4, 5] (right benzene)
        // Ring 3: [5, 6, 7, 8, 9, 0] (left benzene) - shares edge with Ring 2

        let mut mol = crate::smiles::parse_smiles("c1ccc2ccccc2c1").unwrap();

        // Verify molecule structure
        assert_eq!(mol.num_atoms(), 10, "Naphthalene should have 10 carbons");

        // Detect aromaticity
        detect_aromaticity(&mut mol);

        // After detection - ALL atoms should be aromatic
        // INFO: All Atoms are alreay Aromatic because of `parse_smiles` fn.
        for i in 0..10 {
            assert!(
                mol.atom(i).is_aromatic(),
                "Atom {} in naphthalene should be aromatic (fused aromatic rings)",
                i
            );
        }

        // Verify bonds are marked aromatic
        // Check some key bonds in the structure
        let aromatic_bond_count = (0..mol.num_bonds())
            .filter(|&bond_idx| mol.bond(bond_idx).is_aromatic())
            .count();

        assert!(
            aromatic_bond_count >= 10,
            "Expected at least 10 aromatic bonds in naphthalene, found {}",
            aromatic_bond_count
        );

        // Verify bond orders are set to Aromatic
        for bond_idx in 0..mol.num_bonds() {
            if mol.bond(bond_idx).is_aromatic() {
                assert_eq!(
                    mol.bond(bond_idx).order(),
                    BondOrder::Aromatic,
                    "Aromatic bond {} should have BondOrder::Aromatic",
                    bond_idx
                );
            }
        }

        // Test ring detection specifically
        let rings = find_rings(&mol);

        // Should find at least 2 rings (the two fused benzenes)
        // Might find more depending on the algorithm (e.g., outer perimeter)
        assert!(
            rings.len() >= 2,
            "Should detect at least 2 rings in naphthalene, found {}",
            rings.len()
        );

        // Check that we found 6-membered rings
        let six_membered_rings: Vec<_> = rings.iter().filter(|ring| ring.len() == 6).collect();

        assert!(
            six_membered_rings.len() >= 2,
            "Should find at least 2 six-membered rings in naphthalene, found {}",
            six_membered_rings.len()
        );

        // Verify each 6-membered ring is aromatic
        for ring in six_membered_rings {
            assert!(
                is_aromatic_ring(&mol, ring),
                "Six-membered ring {:?} should be aromatic",
                ring
            );

            // Verify all atoms in the ring are marked aromatic
            for &atom_idx in ring.iter() {
                assert!(
                    mol.atom(atom_idx).is_aromatic(),
                    "Atom {} in aromatic ring should be marked aromatic",
                    atom_idx
                );
            }
        }
    }
}
