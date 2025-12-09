mod arguments;
mod atom_env;
mod env_generator;
mod invariants;

pub use arguments::MorganArguments;
pub use atom_env::MorganAtomEnv;
pub use env_generator::MorganEnvGenerator;
pub use invariants::{
    MorganAtomInvGenerator, MorganBondInvGenerator, MorganFeatureAtomInvGenerator,
};

use crate::errors::FingerprintError;
use crate::fingerprint::{FingerprintFuncArguments, FingerprintGenerator, ROMol};
use bitvec::prelude::*;
use std::collections::HashMap;

/// Morgan fingerprint main interface
pub struct MorganFingerprint;

impl MorganFingerprint {
    /// Get Morgan fingerprint as BitVec
    pub fn get_fingerprint_as_bitvec(
        mol: &ROMol,
        radius: u32,
        n_bits: u32,
        invariants: Option<&[u32]>,
        from_atoms: Option<&[u32]>,
        use_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
    ) -> Result<BitVec, FingerprintError> {
        let generator = Self::get_morgan_generator(
            radius,
            false,
            use_chirality,
            use_bond_types,
            only_nonzero_invariants,
            false,
            n_bits,
        );

        let args = FingerprintFuncArguments {
            custom_atom_invariants: invariants,
            from_atoms,
            additional_output: None,
        };

        generator.get_fingerprint(mol, args)
    }

    /// Get sparse count fingerprint
    pub fn get_sparse_count_fingerprint(
        mol: &ROMol,
        radius: u32,
        invariants: Option<&[u32]>,
        from_atoms: Option<&[u32]>,
        use_chirality: bool,
        use_bond_types: bool,
        use_counts: bool,
        only_nonzero_invariants: bool,
    ) -> Result<HashMap<u32, u32>, FingerprintError> {
        // Implementation for sparse count fingerprint
        Ok(HashMap::new())
    }

    fn get_morgan_generator(
        radius: u32,
        count_simulation: bool,
        include_chirality: bool,
        use_bond_types: bool,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        fp_size: u32,
    ) -> FingerprintGenerator<u32> {
        let env_generator = Box::new(MorganEnvGenerator::new());
        let arguments = Box::new(MorganArguments::new(
            radius,
            count_simulation,
            include_chirality,
            only_nonzero_invariants,
            vec![1, 2, 4, 8],
            fp_size,
            include_redundant_environments,
            use_bond_types,
        ));
        let atom_inv_gen = Box::new(MorganAtomInvGenerator::new(true));
        let bond_inv_gen = Box::new(MorganBondInvGenerator::new(
            use_bond_types,
            include_chirality,
        ));

        FingerprintGenerator::new(env_generator, arguments, atom_inv_gen, bond_inv_gen)
    }
}
