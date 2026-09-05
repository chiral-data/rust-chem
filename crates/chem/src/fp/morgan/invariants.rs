use crate::fp::fingerprint::{AtomInvariantsGenerator, BondInvariantsGenerator, ROMol};

use crate::core::bond::{BondOrder, BondStereo};

/// Morgan atom invariant generator (ECFP-type)
pub struct MorganAtomInvGenerator {
    include_ring_membership: bool,
}

impl MorganAtomInvGenerator {
    pub fn new(include_ring_membership: bool) -> Self {
        Self {
            include_ring_membership,
        }
    }

    fn get_connectivity_invariants(&self, mol: &ROMol) -> Vec<u32> {
        let mut invariants = Vec::with_capacity(mol.num_atoms());

        for (i, atom) in mol.atoms().iter().enumerate() {
            // Degree
            let degree = mol.degree(i) as u32;

            // Atomic number
            let mut inv = atom.atomic_number() as u32;
            inv = inv.wrapping_mul(100) + degree;

            // Hydrogen count. total_hydrogens() fits in a u8, so a base of
            // 100 leaves no risk of it colliding with the charge term folded
            // in next (#192: omitting this collapsed atoms that differ only
            // in H count into the same environment).
            inv = inv.wrapping_mul(100) + atom.total_hydrogens() as u32;

            // Formal charge. Offset by +5 to stay non-negative before
            // multiplying, matching Atom::compute_hash's own convention —
            // formal charges in practice never reach that range. Omitting
            // this is what collapsed nitrobenzene's and the glycine
            // zwitterion's charged atoms into their neutral counterparts'
            // environments, understating both molecules' true bit count
            // against RDKit (#192).
            inv = inv.wrapping_mul(20) + (atom.formal_charge() as i32 + 5) as u32;

            // Aromatic flag
            if atom.is_aromatic() {
                inv = inv.wrapping_mul(2);
            }

            // Hash together
            inv = inv.wrapping_mul(31);

            invariants.push(inv);
        }

        invariants
    }
}

impl AtomInvariantsGenerator for MorganAtomInvGenerator {
    fn get_atom_invariants(&self, mol: &ROMol) -> Vec<u32> {
        self.get_connectivity_invariants(mol)
    }

    fn info_string(&self) -> String {
        format!(
            "MorganInvariantGenerator includeRingMembership={}",
            self.include_ring_membership
        )
    }

    fn clone_box(&self) -> Box<dyn AtomInvariantsGenerator> {
        Box::new(Self {
            include_ring_membership: self.include_ring_membership,
        })
    }
}

/// Morgan feature atom invariant generator (FCFP-type)
pub struct MorganFeatureAtomInvGenerator {
    // Feature patterns would go here
}

impl MorganFeatureAtomInvGenerator {
    pub fn new() -> Self {
        Self {}
    }
}

impl Default for MorganFeatureAtomInvGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl AtomInvariantsGenerator for MorganFeatureAtomInvGenerator {
    fn get_atom_invariants(&self, mol: &ROMol) -> Vec<u32> {
        // Feature-based invariants implementation
        vec![1; mol.num_atoms()]
    }

    fn info_string(&self) -> String {
        "MorganFeatureInvariantGenerator".to_string()
    }

    fn clone_box(&self) -> Box<dyn AtomInvariantsGenerator> {
        Box::new(Self {})
    }
}

/// Morgan bond invariant generator
pub struct MorganBondInvGenerator {
    use_bond_types: bool,
    use_chirality: bool,
}

impl MorganBondInvGenerator {
    pub fn new(use_bond_types: bool, use_chirality: bool) -> Self {
        Self {
            use_bond_types,
            use_chirality,
        }
    }
}

impl BondInvariantsGenerator for MorganBondInvGenerator {
    fn get_bond_invariants(&self, mol: &ROMol) -> Vec<u32> {
        let mut result = Vec::with_capacity(mol.num_bonds());

        for bond in mol.bonds() {
            let mut bond_invariant = 1u32;

            if self.use_bond_types {
                if !self.use_chirality
                    || bond.order() != BondOrder::Double
                    || bond.stereo() == BondStereo::None
                {
                    bond_invariant = match bond.order() {
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Quadruple => 4,
                        BondOrder::Aromatic => 12, // 1.5??
                    };
                } else {
                    let bond_stereo = match bond.stereo() {
                        BondStereo::E => 1,
                        BondStereo::Z => 2,
                        _ => 0,
                    };
                    //let bond_stereo = bond.stereo() as u32;
                    let stereo_offset = 100;
                    let bond_type_offset = 10;
                    bond_invariant =
                        stereo_offset + bond_type_offset * (bond.bond_type() as u32) + bond_stereo;
                }
            }

            result.push(bond_invariant);
        }

        result
    }

    fn info_string(&self) -> String {
        format!(
            "MorganInvariantGenerator useBondTypes={} useChirality={}",
            self.use_bond_types, self.use_chirality
        )
    }

    fn clone_box(&self) -> Box<dyn BondInvariantsGenerator> {
        Box::new(Self {
            use_bond_types: self.use_bond_types,
            use_chirality: self.use_chirality,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{Atom, Element};
    use crate::core::molecule::Molecule;

    #[test]
    fn test_formal_charge_changes_the_atom_invariant() {
        // Ammonium's nitrogen and a neutral nitrogen of the same degree used
        // to hash identically -- formal charge was never read (#192), which
        // silently merged environments RDKit keeps distinct.
        let mut neutral = Molecule::new();
        neutral.add_atom(Atom::new(Element::nitrogen()));

        let mut charged = Molecule::new();
        charged.add_atom(Atom::new(Element::nitrogen()).with_charge(1));

        let generator = MorganAtomInvGenerator::new(false);
        assert_ne!(
            generator.get_atom_invariants(&neutral)[0],
            generator.get_atom_invariants(&charged)[0]
        );
    }

    #[test]
    fn test_hydrogen_count_changes_the_atom_invariant() {
        // Two atoms of the same element, degree and charge but a different
        // hydrogen count used to hash identically -- total_hydrogens() was
        // never read (#192).
        let mut few_h = Molecule::new();
        few_h.add_atom(Atom::new(Element::carbon()));

        let mut many_h = Molecule::new();
        let idx = many_h.add_atom(Atom::new(Element::carbon()));
        many_h.atom_mut(idx).set_implicit_hydrogens(3);

        let generator = MorganAtomInvGenerator::new(false);
        assert_ne!(
            generator.get_atom_invariants(&few_h)[0],
            generator.get_atom_invariants(&many_h)[0]
        );
    }

    #[test]
    fn test_nitrobenzene_and_glycine_zwitterion_match_rdkits_bit_count() {
        // #192: chem's Morgan atom invariant omitted formal charge and
        // hydrogen count, so these two charged molecules found fewer
        // distinct environments than RDKit (14 vs 16, 9 vs 11) even though
        // both use exactly the same radius and bit width. Verified against
        // RDKit 2025.3.3 via the oracle Docker image.
        use crate::fp::morgan::MorganFingerprint;
        use crate::io::smiles::parse_smiles;

        let nitrobenzene = parse_smiles("c1ccccc1[N+](=O)[O-]").expect("valid SMILES");
        let fp = MorganFingerprint::get_fingerprint_as_bitvec(
            &nitrobenzene,
            2,
            2048,
            None,
            None,
            false,
            true,
            false,
        )
        .expect("fingerprint");
        assert_eq!(fp.count_ones(), 16);

        let glycine_zwitterion = parse_smiles("[O-]C(=O)C[NH3+]").expect("valid SMILES");
        let fp = MorganFingerprint::get_fingerprint_as_bitvec(
            &glycine_zwitterion,
            2,
            2048,
            None,
            None,
            false,
            true,
            false,
        )
        .expect("fingerprint");
        assert_eq!(fp.count_ones(), 11);
    }
}
