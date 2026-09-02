use crate::core::{bond::BondOrder, molecule::Molecule, rings::find_sssr};

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
    // Computed once: marking a ring changes bond *orders*, not topology, so the
    // ring set stays valid across passes.
    let rings = find_sssr(mol);

    // Repeated until nothing new is found, because marking one ring can enable
    // another. `mark_ring_aromatic` sets its bonds to `BondOrder::Aromatic`,
    // and `has_in_ring_pi_bond` accepts `Aromatic` — so a fused neighbour that
    // failed earlier in a pass may succeed after its partner is marked.
    //
    // A single pass made detection depend on the order `find_sssr` happened to
    // return rings, and therefore on how the input SMILES was written: Kekulé
    // naphthalene needed two passes while anthracene succeeded in one. Two
    // spellings of the same compound could disagree.
    //
    // Terminates because each pass either marks at least one previously
    // unmarked ring or stops, so it runs at most `rings.len() + 1` times.
    loop {
        let mut marked_any = false;
        for ring in &rings {
            if ring_bonds_all_aromatic(mol, ring.atoms()) {
                continue;
            }
            if is_aromatic_ring(mol, ring.atoms()) {
                mark_ring_aromatic(mol, ring.atoms());
                marked_any = true;
            }
        }
        if !marked_any {
            break;
        }
    }
}

/// Whether this ring has already been marked.
///
/// Tests the ring's own bonds rather than its atoms: in a fused system an atom
/// can be aromatic by virtue of a neighbouring ring while this ring's bonds are
/// untouched, and treating that as "already done" is what would stop the loop
/// one ring short.
fn ring_bonds_all_aromatic(mol: &Molecule, ring: &[usize]) -> bool {
    (0..ring.len()).all(|i| {
        let next = (i + 1) % ring.len();
        mol.graph()
            .get_bond(ring[i], ring[next])
            .is_some_and(|idx| mol.bond(idx).order() == BondOrder::Aromatic)
    })
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
            // The discriminator is an in-ring pi bond, not the degree.
            //
            // A pyridine-type nitrogen is part of a ring double bond and
            // contributes one electron, exactly as a ring carbon does. A
            // pyrrole-type nitrogen has no ring double bond and donates its
            // lone pair, contributing two.
            //
            // This used to key off `degree == 2`, with the two cases labelled
            // the wrong way round. Both types have heavy-atom degree 2, so
            // degree cannot tell them apart: pyridine's five carbons plus a
            // nitrogen credited 2 summed to 7 and failed Huckel, while pyrrole
            // came out right by coincidence. Imidazole has one nitrogen of each
            // type and needs them counted differently, which no degree-based
            // rule can do.
            7 => {
                if has_in_ring_pi_bond(mol, atom_idx, ring) {
                    pi_electrons += 1;
                } else {
                    pi_electrons += 2;
                }
            }
            // Same rule for oxygen and sulfur, for the same reason. Furan and
            // thiophene are unaffected — their heteroatom has no ring double
            // bond and still contributes two — but a ring oxygen that does
            // carry one (a pyrylium, say) now contributes one rather than being
            // miscounted.
            8 | 16 => {
                if has_in_ring_pi_bond(mol, atom_idx, ring) {
                    pi_electrons += 1;
                } else {
                    pi_electrons += 2;
                }
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

/// Assigns a Kekulé form to a molecule's aromatic bonds, without changing it.
///
/// Returns one [`BondOrder`] per bond, or `None` when no valid assignment
/// exists. The molecule's own aromatic perception is correct and stays; only a
/// *file* needs the alternating form, so this hands back orders rather than
/// mutating.
///
/// # Why a writer needs this
///
/// Molfile bond type 4 means "aromatic", which is a **query** bond type: it
/// belongs in a substructure search, not in a structure record. Benzene
/// survives it because a reader can kekulise `c1ccccc1` unambiguously. Pyrrole
/// does not — without knowing the nitrogen carries a hydrogen there is no valid
/// Kekulé form to recover, so RDKit rejects the record outright (#194).
///
/// # How
///
/// An aromatic atom needs exactly one double bond unless it donates a lone pair
/// to the ring instead. That distinction is the whole difficulty, and it is
/// decided by valence rather than by element: count what the atom already has,
/// against what its element and charge allow. Pyrrole's N-H reaches its valence
/// of three with two ring bonds and a hydrogen, so it needs nothing more;
/// pyridine's N has only two, so it needs a double bond. Furan's O is full at
/// two.
///
/// What remains is a perfect matching over the aromatic subgraph restricted to
/// the atoms that need one, solved by backtracking. The molecules this crate
/// sees have a handful of aromatic rings each, so the search is small.
///
/// Returning `None` rather than guessing is deliberate. A wrong Kekulé form is
/// wrong chemistry in a file that reads perfectly; a type-4 bond is a
/// portability gap that is already recorded and understood.
pub fn kekulize(mol: &Molecule) -> Option<Vec<BondOrder>> {
    let mut orders: Vec<BondOrder> = (0..mol.num_bonds()).map(|i| mol.bond(i).order()).collect();

    let aromatic: Vec<usize> = (0..mol.num_bonds())
        .filter(|&i| mol.bond(i).order() == BondOrder::Aromatic)
        .collect();
    if aromatic.is_empty() {
        return Some(orders);
    }

    // Does this atom still need a double bond, or is it already satisfied?
    let needs_double = |atom_idx: usize| -> bool {
        let atom = mol.atom(atom_idx);
        let typical = atom.element().typical_valence();
        if typical == 0 {
            return false;
        }
        // An aromatic bond counts as a single bond here: the question is what
        // the atom would still be missing if every aromatic bond were single.
        let used: f64 = mol
            .neighbors(atom_idx)
            .iter()
            .map(|n| match mol.bond(n.bond_idx).order() {
                BondOrder::Aromatic => 1.0,
                other => other.value(),
            })
            .sum();
        let charge = i32::from(atom.formal_charge());
        let allowed = i32::from(typical) + charge;
        let have = used as i32 + i32::from(atom.total_hydrogens());
        allowed - have == 1
    };

    let mut wanted: Vec<usize> = Vec::new();
    for &bond_idx in &aromatic {
        for atom in [mol.bond(bond_idx).atom1(), mol.bond(bond_idx).atom2()] {
            if !wanted.contains(&atom) && needs_double(atom) {
                wanted.push(atom);
            }
        }
    }

    // An odd number of atoms needing a partner cannot be matched at all, and
    // saying so here saves the search from proving it the slow way.
    if !wanted.len().is_multiple_of(2) {
        return None;
    }

    let mut matched: Vec<bool> = vec![false; mol.num_atoms()];
    if match_atoms(mol, &aromatic, &wanted, &mut matched, &mut orders) {
        // Everything left aromatic in a matched system becomes single.
        for &bond_idx in &aromatic {
            if orders[bond_idx] == BondOrder::Aromatic {
                orders[bond_idx] = BondOrder::Single;
            }
        }
        Some(orders)
    } else {
        None
    }
}

/// Backtracking perfect matching over the aromatic bonds.
///
/// Takes the first unmatched atom that needs a partner and tries each aromatic
/// bond to another such atom, undoing the choice when the rest cannot be
/// completed. Depth is half the number of atoms in the aromatic system.
fn match_atoms(
    mol: &Molecule,
    aromatic: &[usize],
    wanted: &[usize],
    matched: &mut [bool],
    orders: &mut [BondOrder],
) -> bool {
    let Some(&atom) = wanted.iter().find(|&&a| !matched[a]) else {
        return true; // everyone has a partner
    };

    for &bond_idx in aromatic {
        if orders[bond_idx] != BondOrder::Aromatic {
            continue;
        }
        let bond = mol.bond(bond_idx);
        let (a, b) = (bond.atom1(), bond.atom2());
        let partner = if a == atom {
            b
        } else if b == atom {
            a
        } else {
            continue;
        };
        if matched[partner] || !wanted.contains(&partner) {
            continue;
        }

        matched[atom] = true;
        matched[partner] = true;
        orders[bond_idx] = BondOrder::Double;

        if match_atoms(mol, aromatic, wanted, matched, orders) {
            return true;
        }

        matched[atom] = false;
        matched[partner] = false;
        orders[bond_idx] = BondOrder::Aromatic;
    }

    false
}

#[cfg(test)]
mod tests {

    use crate::io::smiles::parse_smiles;

    fn kekule_orders(smiles: &str) -> Vec<BondOrder> {
        let mol = parse_smiles(smiles).expect(smiles);
        kekulize(&mol).unwrap_or_else(|| panic!("{smiles} has no Kekule form"))
    }

    #[test]
    fn test_every_aromatic_system_in_the_corpus_kekulises() {
        // `None` here means the needs-a-double-bond rule is wrong, which is a
        // cheaper diagnostic than reading the matching by hand.
        for smiles in [
            "c1ccccc1",
            "c1ccncc1",
            "c1ccoc1",
            "c1ccsc1",
            "c1cc[nH]c1",
            "c1cnc[nH]1",
            "C1CCc2ccccc2C1",
            "[O-][N+](=O)c1ccccc1",
        ] {
            let mol = parse_smiles(smiles).expect(smiles);
            assert!(kekulize(&mol).is_some(), "{smiles} has no Kekule form");
        }
    }

    #[test]
    fn test_no_aromatic_bond_survives_kekulisation() {
        for smiles in ["c1ccccc1", "c1cc[nH]c1", "c1cnc[nH]1"] {
            assert!(
                !kekule_orders(smiles).contains(&BondOrder::Aromatic),
                "{smiles} kept a query bond type"
            );
        }
    }

    #[test]
    fn test_benzene_alternates_three_and_three() {
        let orders = kekule_orders("c1ccccc1");
        assert_eq!(
            orders.iter().filter(|o| **o == BondOrder::Double).count(),
            3
        );
        assert_eq!(
            orders.iter().filter(|o| **o == BondOrder::Single).count(),
            3
        );
    }

    #[test]
    fn test_a_pyrrole_nitrogen_takes_no_double_bond() {
        // The case that makes this more than counting: the N-H donates its lone
        // pair, so it is already satisfied and must stay single-bonded. Treating
        // it like pyridine's N is what leaves pyrrole unwritable.
        let mol = parse_smiles("c1cc[nH]c1").expect("valid SMILES");
        let orders = kekulize(&mol).expect("has a Kekule form");
        let nitrogen = (0..mol.num_atoms())
            .find(|&i| mol.atom(i).element().symbol() == "N")
            .expect("a nitrogen");
        for neighbour in mol.neighbors(nitrogen) {
            assert_eq!(
                orders[neighbour.bond_idx],
                BondOrder::Single,
                "the pyrrole nitrogen took a double bond"
            );
        }
    }

    #[test]
    fn test_a_molecule_with_no_aromatic_bonds_is_returned_unchanged() {
        let mol = parse_smiles("CCO").expect("valid SMILES");
        let orders = kekulize(&mol).expect("nothing to do");
        assert!(orders.iter().all(|o| *o == BondOrder::Single));
    }
    use super::*;

    #[test]
    fn test_benzene_aromaticity() {
        let mut mol = crate::io::smiles::parse_smiles("C1=CC=CC=C1").unwrap();
        detect_aromaticity(&mut mol);

        for i in 0..6 {
            assert!(mol.atom(i).is_aromatic());
        }
    }

    #[test]
    fn test_cyclopentane_not_aromatic() {
        let mut mol = crate::io::smiles::parse_smiles("C1CCCC1").unwrap();
        detect_aromaticity(&mut mol);

        assert!(!mol.atom(0).is_aromatic());
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        // Fully saturated 6-membered all-carbon ring, no double bonds at
        // all -- previously flagged aromatic because every carbon
        // contributed a pi-electron unconditionally, landing on exactly 6
        // (the same count real benzene gets) regardless of saturation.
        let mut mol = crate::io::smiles::parse_smiles("C1CCCCC1").unwrap();
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
        let mut mol = crate::io::smiles::parse_smiles("O=C1CCCCC1").unwrap();
        detect_aromaticity(&mut mol);

        for i in 0..mol.num_atoms() {
            assert!(!mol.atom(i).is_aromatic(), "atom {} wrongly aromatic", i);
        }
    }

    /// Every case here starts from **Kekulé** SMILES — explicit alternating
    /// double bonds, uppercase atoms — because that is the only input that
    /// reaches this code. Lowercase aromatic SMILES is already flagged aromatic
    /// by the parser, so a test written that way passes whatever
    /// `detect_aromaticity` does.
    ///
    /// The test this replaces used `c1ccc2ccccc2c1` and said so in its own
    /// comment: "All Atoms are already Aromatic because of `parse_smiles`". It
    /// asserted a property the parser had already established, which is why
    /// #160 survived having a polycyclic test at all.
    fn aromatic_count(smiles: &str) -> usize {
        let mut mol = crate::io::smiles::parse_smiles(smiles).expect(smiles);
        assert_eq!(
            mol.atoms().iter().filter(|a| a.is_aromatic()).count(),
            0,
            "{smiles} must start non-aromatic or this test proves nothing"
        );
        detect_aromaticity(&mut mol);
        mol.atoms().iter().filter(|a| a.is_aromatic()).count()
    }

    #[test]
    fn test_pyridine_type_nitrogen_is_aromatic() {
        // The inverted electron count made all of these return zero: five ring
        // carbons contributing 1 each plus a nitrogen credited 2 summed to 7,
        // which fails Hückel.
        assert_eq!(aromatic_count("C1=CC=NC=C1"), 6, "pyridine");
        assert_eq!(aromatic_count("C1=CN=CN=C1"), 6, "pyrimidine");
        assert_eq!(aromatic_count("C1=CN=CC=N1"), 6, "pyrazine");
    }

    #[test]
    fn test_pyrrole_type_nitrogen_still_donates_two() {
        // Right before the fix as well, but by coincidence — the old rule keyed
        // off degree, and a pyrrole nitrogen also has heavy-atom degree 2. This
        // pins it so the new rule is not credited for luck.
        assert_eq!(aromatic_count("C1=CNC=C1"), 5, "pyrrole");
    }

    #[test]
    fn test_a_ring_with_one_nitrogen_of_each_type() {
        // Imidazole is the case no degree-based rule can get right: one
        // pyridine-type nitrogen contributing 1, one pyrrole-type contributing
        // 2, and three carbons — six in total.
        assert_eq!(aromatic_count("C1=CN=CN1"), 5, "imidazole");
    }

    #[test]
    fn test_oxygen_and_sulfur_heteroaromatics() {
        assert_eq!(aromatic_count("C1=COC=C1"), 5, "furan");
        assert_eq!(aromatic_count("C1=CSC=C1"), 5, "thiophene");
    }

    #[test]
    fn test_fused_rings_are_all_found_regardless_of_ring_order() {
        // A single pass marked naphthalene's rings one at a time, so it found
        // 6 of 10 — while anthracene happened to come out complete. Whether
        // detection worked depended on how the SMILES was written.
        assert_eq!(aromatic_count("C1=CC2=CC=CC=C2C=C1"), 10, "naphthalene");
        assert_eq!(aromatic_count("C1=CC2=CC=CC=C2N=C1"), 10, "quinoline");
        assert_eq!(
            aromatic_count("C1=CC=C2C=C3C=CC=CC3=CC2=C1"),
            14,
            "anthracene"
        );
    }

    #[test]
    fn test_substituents_do_not_disturb_the_ring() {
        assert_eq!(aromatic_count("CC1=CC=CC=C1"), 6, "toluene");
        assert_eq!(aromatic_count("OC1=CC=CC=C1"), 6, "phenol");
    }

    #[test]
    fn test_saturated_and_partly_saturated_rings_stay_non_aromatic() {
        // The other direction, and the one a looser electron count breaks: #65
        // was this function calling cyclohexane aromatic. Heteroatom rings are
        // included because the new rule credits a lone-pair heteroatom 2
        // electrons, and 2 plus nothing must not reach 6.
        for smiles in [
            "C1CCCCC1",   // cyclohexane
            "C1=CCCCC1",  // cyclohexene
            "C1CCCC1",    // cyclopentane
            "O=C1CCCCC1", // cyclohexanone
            "C1CCNCC1",   // piperidine
            "C1CCOC1",    // tetrahydrofuran
            "C1COCCO1",   // dioxane
        ] {
            assert_eq!(aromatic_count(smiles), 0, "{smiles} must not be aromatic");
        }
    }

    #[test]
    fn test_detection_is_idempotent() {
        // The fixed-point loop must settle rather than oscillate, and a second
        // call must not widen the result.
        let mut mol = crate::io::smiles::parse_smiles("C1=CC2=CC=CC=C2N=C1").expect("quinoline");
        detect_aromaticity(&mut mol);
        let first = mol.atoms().iter().filter(|a| a.is_aromatic()).count();
        detect_aromaticity(&mut mol);
        let second = mol.atoms().iter().filter(|a| a.is_aromatic()).count();
        assert_eq!(first, 10);
        assert_eq!(first, second);
    }
}
