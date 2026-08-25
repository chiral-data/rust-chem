mod generator;

pub use generator::{
    AdditionalOutput, AtomEnvironment, AtomEnvironmentGenerator, AtomInvariantsGenerator,
    BondInvariantsGenerator, FingerprintArguments, FingerprintFuncArguments, FingerprintGenerator,
};

pub use crate::core::molecule::Molecule as ROMol;

use std::collections::HashMap;

/// Bit info map type: bitId -> Vec<(atomId, radius)>
pub type BitInfoMap = HashMap<u32, Vec<(u32, u32)>>;

/// Additional output structure for fingerprint generation
#[derive(Default)]
pub struct AdditionalOutputData {
    pub bit_info_map: Option<BitInfoMap>,
    pub atom_counts: Option<Vec<u32>>,
    pub atom_to_bits: Option<Vec<Vec<u64>>>,
    pub atoms_per_bit: Option<HashMap<u64, Vec<Vec<u32>>>>,
}

impl AdditionalOutputData {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn allocate_bit_info_map(&mut self) {
        self.bit_info_map = Some(HashMap::new());
    }

    pub fn allocate_atom_counts(&mut self, num_atoms: usize) {
        self.atom_counts = Some(vec![0; num_atoms]);
    }

    pub fn allocate_atom_to_bits(&mut self, num_atoms: usize) {
        self.atom_to_bits = Some(vec![Vec::new(); num_atoms]);
    }

    pub fn allocate_atoms_per_bit(&mut self) {
        self.atoms_per_bit = Some(HashMap::new());
    }
}
