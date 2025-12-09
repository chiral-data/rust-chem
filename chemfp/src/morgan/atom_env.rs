use crate::fingerprint::{AdditionalOutputData, AtomEnvironment, FingerprintArguments, ROMol};

/// Morgan atom environment
pub struct MorganAtomEnv {
    code: u32,
    atom_id: u32,
    layer: u32,
    num_atoms: usize,
}

impl MorganAtomEnv {
    pub fn new(code: u32, atom_id: u32, layer: u32, mol: &ROMol) -> Self {
        Self {
            code,
            atom_id,
            layer,
            num_atoms: mol.num_atoms(),
        }
    }
}

impl AtomEnvironment<u32> for MorganAtomEnv {
    fn get_bit_id(
        &self,
        arguments: &dyn FingerprintArguments,
        _atom_invariants: &[u32],
        _bond_invariants: &[u32],
        _additional_output: Option<&mut AdditionalOutputData>,
        _hash_results: bool,
        _fp_size: u64,
    ) -> u32 {
        let fp_size = if _fp_size == 0 {
            arguments.fp_size() as u64
        } else {
            _fp_size
        };

        let bit_pos = (self.code as u64 % fp_size) as u32;

        eprintln!(
            "DEBUG: code={} fp_size={} → bit={}",
            self.code, fp_size, bit_pos
        );

        bit_pos
        //eprintln!(
        //    "DEBUG: MorganAtomEnv code={} atom={} layer={}",
        //    self.code, self.atom_id, self.layer
        //);
        ////self.code
        //(self.code as u64 % _fp_size) as u32
    }

    fn update_additional_output(&self, output: &mut AdditionalOutputData, bit_id: u64) {
        if let Some(ref mut bit_info) = output.bit_info_map {
            bit_info
                .entry(bit_id as u32)
                .or_insert_with(Vec::new)
                .push((self.atom_id, self.layer));
        }

        if let Some(ref mut counts) = output.atom_counts {
            counts[self.atom_id as usize] += 1;
        }

        if let Some(ref mut atom_to_bits) = output.atom_to_bits {
            atom_to_bits[self.atom_id as usize].push(bit_id);
        }
    }
}
