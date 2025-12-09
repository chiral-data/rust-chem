use crate::fingerprint::{AtomInvariantsGenerator, BondInvariantsGenerator, ROMol};

use chemcore::bond::{BondOrder, BondStereo};

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
            let mut inv = 0u32;

            // Atomic number
            inv = atom.atomic_number() as u32;

            // Degree
            let degree = mol.degree(i) as u32;
            inv = inv.wrapping_mul(100) + degree;

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
                let bond_type_value = match bond.order() {
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Quadruple => 4,
                    BondOrder::Aromatic => 12, // RDKit uses 12!
                };

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
