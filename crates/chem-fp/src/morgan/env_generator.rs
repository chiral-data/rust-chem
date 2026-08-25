use crate::errors::FingerprintError;
use crate::fingerprint::{AtomEnvironment, AtomEnvironmentGenerator, FingerprintArguments, ROMol};
use crate::hash;
use crate::morgan::{MorganArguments, MorganAtomEnv};
use bitvec::prelude::*;
use chem_core::atom::Chirality;
use fxhash::FxHashSet;

/// Morgan environment generator
pub struct MorganEnvGenerator;

impl Default for MorganEnvGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl MorganEnvGenerator {
    pub fn new() -> Self {
        Self
    }
}

impl AtomEnvironmentGenerator<u32> for MorganEnvGenerator {
    fn get_environments(
        &self,
        mol: &ROMol,
        arguments: &dyn FingerprintArguments,
        from_atoms: Option<&[u32]>,
        _ignore_atoms: Option<&[u32]>,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
    ) -> Result<Vec<Box<dyn AtomEnvironment<u32>>>, FingerprintError> {
        let morgan_args = arguments
            .as_any()
            .downcast_ref::<MorganArguments>()
            .ok_or(FingerprintError::InvalidMolecule)?;

        let n_atoms = mol.num_atoms();
        let max_results = (morgan_args.radius + 1) as usize * n_atoms;

        let mut result: Vec<Box<dyn AtomEnvironment<u32>>> = Vec::with_capacity(max_results);

        let mut current_invariants = atom_invariants.to_vec();
        let mut next_layer_invariants = vec![0u32; n_atoms];

        let include_atoms = if let Some(from) = from_atoms {
            let mut bs = bitvec![0; n_atoms];
            for &idx in from {
                bs.set(idx as usize, true);
            }
            bs
        } else {
            bitvec![1; n_atoms]
        };

        let mut neighborhoods: FxHashSet<BitVec> = FxHashSet::default();
        neighborhoods.reserve(max_results);

        let mut atom_neighborhoods = vec![bitvec![0; mol.num_bonds()]; n_atoms];
        let mut round_atom_neighborhoods = atom_neighborhoods.clone();
        let mut dead_atoms = bitvec![0; n_atoms];

        // Add round 0 invariants
        for i in 0..n_atoms {
            if include_atoms[i]
                && (!morgan_args.only_nonzero_invariants || current_invariants[i] != 0)
            {
                result.push(Box::new(MorganAtomEnv::new(
                    current_invariants[i],
                    i as u32,
                    0,
                    //mol,
                )));
            }
        }

        // Iterate through layers
        for layer in 0..morgan_args.radius {
            let mut all_neighborhoods_this_round: Vec<(BitVec, u32, u32)> = Vec::new();

            for atom_idx in 0..n_atoms {
                if dead_atoms[atom_idx] {
                    continue;
                }

                let atom = &mol.atoms()[atom_idx];
                if mol.degree(atom_idx) == 0 {
                    dead_atoms.set(atom_idx, true);
                    continue;
                }

                let mut neighborhood_invariants: Vec<(i32, u32)> = Vec::new();

                // Get neighbors through bonds
                for neighbor in mol.neighbors(atom_idx) {
                    let other_idx = neighbor.atom_idx;

                    round_atom_neighborhoods[atom_idx].set(neighbor.bond_idx, true);
                    round_atom_neighborhoods[atom_idx] |= &atom_neighborhoods[other_idx];

                    let bt = bond_invariants[neighbor.bond_idx] as i32;
                    neighborhood_invariants.push((bt, current_invariants[other_idx]));
                }

                // Sort neighbors
                neighborhood_invariants.sort_unstable();

                // Calculate new invariant
                let mut invar = layer;
                let current_inv = current_invariants[atom_idx];

                hash::hash_combine(&mut invar, &current_inv);

                for (bt, neighbor_inv) in &neighborhood_invariants {
                    hash::hash_pair(&mut invar, *bt, *neighbor_inv);
                }

                // Handle chirality if needed
                if morgan_args.include_chirality
                    && !matches!(atom.chirality(), Chirality::None | Chirality::Unspecified)
                //&& atom.chiral_tag != crate::fingerprint::ChiralTag::Unspecified
                {
                    let chiral_code = match atom.chirality() {
                        Chirality::Clockwise => 3u32,
                        Chirality::CounterClockwise => 2u32,
                        _ => 1u32,
                    };
                    hash::hash_combine(&mut invar, &chiral_code);
                    //if let Some(ref cip) = atom.cip_code {
                    //    let chiral_code = match cip.as_str() {
                    //        "R" => 3,
                    //        "S" => 2,
                    //        _ => 1,
                    //    };
                    //    Self::hash_combine(&mut invar, chiral_code);
                    //}
                }

                next_layer_invariants[atom_idx] = invar;
                all_neighborhoods_this_round.push((
                    round_atom_neighborhoods[atom_idx].clone(),
                    invar,
                    atom_idx as u32,
                ));
            }

            // Sort neighborhoods
            all_neighborhoods_this_round.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

            // Add unique environments
            for (neighborhood, invar, atom_idx) in all_neighborhoods_this_round {
                if morgan_args.include_redundant_environments
                    || !neighborhoods.contains(&neighborhood)
                {
                    if (!morgan_args.only_nonzero_invariants
                        || atom_invariants[atom_idx as usize] != 0)
                        && include_atoms[atom_idx as usize]
                    {
                        result.push(Box::new(MorganAtomEnv::new(
                            invar,
                            atom_idx,
                            layer + 1,
                            //mol,
                        )));
                        neighborhoods.insert(neighborhood);
                    }
                } else {
                    dead_atoms.set(atom_idx as usize, true);
                }
            }

            // Swap invariants for next iteration
            std::mem::swap(&mut current_invariants, &mut next_layer_invariants);
            next_layer_invariants.fill(0);
            atom_neighborhoods = round_atom_neighborhoods.clone();
        }

        Ok(result)
    }

    fn get_result_size(&self) -> u32 {
        u32::MAX
    }

    fn info_string(&self) -> String {
        "MorganEnvironmentGenerator".to_string()
    }
}
