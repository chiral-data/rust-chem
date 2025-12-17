use super::AdditionalOutputData;
use crate::errors::FingerprintError;
use bitvec::prelude::*;
use chemcore::molecule::Molecule;

pub type ROMol = Molecule;

/// Arguments for fingerprint generation
pub trait FingerprintArguments: Send + Sync {
    fn info_string(&self) -> String;
    fn fp_size(&self) -> u32;
    fn count_simulation(&self) -> bool;
    fn as_any(&self) -> &dyn std::any::Any;
}

/// Atom environment trait
pub trait AtomEnvironment<OutputType>: Send + Sync {
    fn get_bit_id(
        &self,
        arguments: &dyn FingerprintArguments,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
        additional_output: Option<&mut AdditionalOutputData>,
        hash_results: bool,
        fp_size: u64,
    ) -> OutputType;

    fn update_additional_output(&self, output: &mut AdditionalOutputData, bit_id: u64);
}

/// Atom invariants generator trait
pub trait AtomInvariantsGenerator: Send + Sync {
    fn get_atom_invariants(&self, mol: &ROMol) -> Vec<u32>;
    fn info_string(&self) -> String;
    fn clone_box(&self) -> Box<dyn AtomInvariantsGenerator>;
}

/// Bond invariants generator trait
pub trait BondInvariantsGenerator: Send + Sync {
    fn get_bond_invariants(&self, mol: &ROMol) -> Vec<u32>;
    fn info_string(&self) -> String;
    fn clone_box(&self) -> Box<dyn BondInvariantsGenerator>;
}

/// Atom environment generator trait
pub trait AtomEnvironmentGenerator<OutputType>: Send + Sync {
    fn get_environments(
        &self,
        mol: &ROMol,
        arguments: &dyn FingerprintArguments,
        from_atoms: Option<&[u32]>,
        ignore_atoms: Option<&[u32]>,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
    ) -> Result<Vec<Box<dyn AtomEnvironment<OutputType>>>, FingerprintError>;

    fn get_result_size(&self) -> OutputType;
    fn info_string(&self) -> String;
}

/// Function arguments for fingerprint generation
pub struct FingerprintFuncArguments<'a> {
    pub custom_atom_invariants: Option<&'a [u32]>,
    pub from_atoms: Option<&'a [u32]>,
    pub additional_output: Option<&'a mut AdditionalOutputData>,
}

/// Additional output interface
pub trait AdditionalOutput {
    fn allocate_bit_info_map(&mut self);
    fn allocate_atom_counts(&mut self, num_atoms: usize);
    fn allocate_atom_to_bits(&mut self, num_atoms: usize);
}

/// Main fingerprint generator
pub struct FingerprintGenerator<OutputType> {
    env_generator: Box<dyn AtomEnvironmentGenerator<OutputType>>,
    arguments: Box<dyn FingerprintArguments>,
    atom_inv_generator: Box<dyn AtomInvariantsGenerator>,
    bond_inv_generator: Box<dyn BondInvariantsGenerator>,
}

impl<OutputType: Copy + std::hash::Hash + Eq + Into<u64>> FingerprintGenerator<OutputType> {
    pub fn new(
        env_generator: Box<dyn AtomEnvironmentGenerator<OutputType>>,
        arguments: Box<dyn FingerprintArguments>,
        atom_inv_generator: Box<dyn AtomInvariantsGenerator>,
        bond_inv_generator: Box<dyn BondInvariantsGenerator>,
    ) -> Self {
        Self {
            env_generator,
            arguments,
            atom_inv_generator,
            bond_inv_generator,
        }
    }

    pub fn get_fingerprint(
        &self,
        mol: &ROMol,
        args: FingerprintFuncArguments,
    ) -> Result<BitVec, FingerprintError> {
        let atom_invariants = if let Some(inv) = args.custom_atom_invariants {
            inv.to_vec()
        } else {
            self.atom_inv_generator.get_atom_invariants(mol)
        };

        let bond_invariants = self.bond_inv_generator.get_bond_invariants(mol);

        let environments = self.env_generator.get_environments(
            mol,
            self.arguments.as_ref(),
            args.from_atoms,
            None,
            &atom_invariants,
            &bond_invariants,
        )?;

        let fp_size = self.arguments.fp_size() as usize;
        let mut result = bitvec![0; fp_size];

        for env in environments {
            let bit_id = env.get_bit_id(
                self.arguments.as_ref(),
                &atom_invariants,
                &bond_invariants,
                None,
                true,
                fp_size as u64,
            );

            // Hash to fit in fp_size (simplified)
            let bit_pos = self.hash_to_size(bit_id, fp_size);
            result.set(bit_pos, true);
        }

        Ok(result)
    }

    fn hash_to_size(&self, value: OutputType, fp_size: usize) -> usize
    where
        OutputType: Into<u64>,
    {
        // Simplified hashing - in production use proper hash function
        //0 % fp_size
        let hash = value.into();
        (hash % fp_size as u64) as usize
    }
}
