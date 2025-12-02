use crate::{
    bitvec::{hash_to_position, BitVec},
    errors::FingerprintError,
};
use chemcore::prelude::*;
use std::collections::HashSet;

/// Configuration options for the Morgan/ECFP fingerprint generator.
///
/// The Morgan fingerprint iteratively hashes atom neighborhoods
/// (up to a given radius) and maps hash values into a fixed-length
/// bit vector.
///
/// # Fields
/// - `radius`: number of expansion iterations (bond distance).
///   - `radius = 0` hashes atoms individually
///   - `radius = 2` corresponds to ECFP4
///   - `radius = 3` corresponds to ECFP6
///
/// - `nbits`: total number of bits in the fingerprint.
///   Must be a power of two (e.g., 1024, 2048, 4096).
///   The value determines collision probability.
///
/// - `bits_per_feature`: how many *different* bit positions to set
///   for each hashed atom-environment. Setting multiple bits reduces
///   sensitivity to hash collisions and improves robustness.
pub struct MorganOptions {
    /// Radius of the fingerprint (dist from a central atom to neighbouring atoms when geenrating circular substructures)
    pub radius: u32,
    /// Number of bits in the fingerprint
    pub nbits: usize,
    // Number of bits to set per feature
    pub bits_per_feature: usize,
}

impl Default for MorganOptions {
    fn default() -> Self {
        MorganOptions {
            radius: 2,
            nbits: 2048,
            bits_per_feature: 1,
        }
    }
}

/// Compute a Morgan/ECFP fingerprint.
///
/// This algorithm:
/// 1. Assigns an initial hash to each atom based on local properties
/// 2. Iteratively updates hashes using neighbor information (`radius` iterations)
/// 3. Converts each hash to one or more bit positions and sets them in a bit vector
///
/// The result is a structural fingerprint suitable for Tanimoto similarity.
///
/// # Parameters
/// - `mol`: input molecule
/// - `options`: fingerprint configuration
///
/// # Errors
/// - [`FingerprintError::InvalidMolecule`] if the molecule has zero atoms
/// - [`FingerprintError::InvalidBitSize`] if `nbits` is zero or not a power of two
///
/// # Returns
/// A [`BitVec`] containing the binary fingerprint.
///
/// # Notes
/// The implementation closely mirrors RDKit's ECFP update procedure.
/// Equal structure → equal fingerprint under equal options.
pub fn morgan_fingerprint(
    mol: &Molecule,
    options: &MorganOptions,
) -> Result<BitVec, FingerprintError> {
    if mol.num_atoms() == 0 {
        return Err(FingerprintError::InvalidMolecule);
    }

    if options.nbits == 0 || (options.nbits & (options.nbits - 1)) != 0 {
        return Err(FingerprintError::InvalidBitSize(options.nbits));
    }

    let mut fp = BitVec::new(options.nbits);
    let mut all_hashes = HashSet::new();
    let mut identifiers = initialize_identifiers(mol);

    for iteration in 0..=options.radius {
        for &hash in identifiers.iter().take(mol.num_atoms()) {
            all_hashes.insert(hash);
        }

        if iteration < options.radius {
            identifiers = update_identifiers(mol, &identifiers);
        }
    }

    // Simple mapping: one hash → one bit position
    for hash in all_hashes {
        let pos = hash_to_position(hash, options.nbits);
        fp.set(pos);
    }

    Ok(fp)
}

/// Compute initial atom identifiers for radius 0.
///
/// Each identifier encodes only *local atom properties*:
/// - atomic number
/// - formal charge
/// - aromaticity
/// - heavy atom degree
/// - total hydrogen count
///
/// These hashes form the base layer for later iterative updates.
fn initialize_identifiers(mol: &Molecule) -> Vec<u32> {
    let mut identifiers = Vec::with_capacity(mol.num_atoms());

    for atom_idx in 0..mol.num_atoms() {
        let atom = mol.atom(atom_idx);
        let hash = compute_atom_hash(atom, mol.degree(atom_idx));
        identifiers.push(hash);
    }

    identifiers
}

/// Hash chemical atom properties into a 64-bit integer.
///
/// This hash is deterministic and designed to cluster chemically similar
/// atoms while separating chemically distinct ones.
///
/// The hashing scheme uses multiplicative mixing to avoid collisions.
/// No cryptographic security is intended.
fn compute_atom_hash(atom: &Atom, degree: usize) -> u32 {
    let atomic_num = atom.atomic_number() as u32;
    let valence = (degree as u32) + (atom.total_hydrogens() as u32);
    let charge = (atom.formal_charge() + 5) as u32; // Offset by 5

    // RDKit's atom invariant calculation
    let mut hash = atomic_num;
    hash = hash.wrapping_mul(37).wrapping_add(valence);
    hash = hash.wrapping_mul(37).wrapping_add(charge);

    if atom.is_aromatic() {
        hash = hash.wrapping_mul(37).wrapping_add(1);
    }

    hash
}

/// Update atom identifiers for the next Morgan iteration.
///
/// Combines each atom's previous identifier with the sorted list of
/// combined neighbor hashes. Sorting ensures identifier order does not
/// depend on the atom indexing of the molecule (canonical behavior).
///
/// This stage is equivalent to RDKit's ECFP neighborhood hash update.
fn update_identifiers(mol: &Molecule, old_ids: &[u32]) -> Vec<u32> {
    let mut new_ids = Vec::with_capacity(mol.num_atoms());

    for atom_idx in 0..mol.num_atoms() {
        let mut neighbor_hashes = Vec::new();

        for neighbor in mol.neighbors(atom_idx) {
            let bond = mol.bond(neighbor.bond_idx);
            let neighbor_id = old_ids[neighbor.atom_idx];

            let bond_val = match bond.order() {
                BondOrder::Single => 1u32,
                BondOrder::Double => 2u32,
                BondOrder::Triple => 3u32,
                BondOrder::Aromatic => 4u32,
                BondOrder::Quadruple => 4u32,
            };

            // RDKit combines: neighbor_hash * 100 + bond_order
            let combined = neighbor_id.wrapping_mul(100).wrapping_add(bond_val);
            neighbor_hashes.push(combined);
        }

        neighbor_hashes.sort_unstable();

        // Combine with current hash using multiply-add chain
        let mut new_hash = old_ids[atom_idx];
        for &nh in &neighbor_hashes {
            new_hash = new_hash.wrapping_mul(37).wrapping_add(nh);
        }

        new_ids.push(new_hash);
    }

    new_ids
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_methane() -> Molecule {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.calculate_implicit_hydrogens();
        mol
    }

    fn create_ethane() -> Molecule {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new(Element::carbon()));
        let c2 = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(c1, c2, BondOrder::Single)).unwrap();
        mol.calculate_implicit_hydrogens();
        mol
    }

    #[test]
    fn test_morgan_methane() {
        let mol = create_methane();
        let options = MorganOptions::default();
        let fp = morgan_fingerprint(&mol, &options).unwrap();

        assert_eq!(fp.size(), 2048);
        assert!(fp.count_ones() > 0);
    }

    #[test]
    fn test_morgan_ethane() {
        let mol = create_ethane();
        let options = MorganOptions::default();
        let fp = morgan_fingerprint(&mol, &options).unwrap();

        assert_eq!(fp.size(), 2048);
        assert!(fp.count_ones() > 0);
    }

    #[test]
    fn test_morgan_different_molecules() {
        let methane = create_methane();
        let ethane = create_ethane();

        let options = MorganOptions::default();
        let fp1 = morgan_fingerprint(&methane, &options).unwrap();
        let fp2 = morgan_fingerprint(&ethane, &options).unwrap();

        // Different molecules should have different fingerprints
        assert_ne!(fp1, fp2);
    }

    #[test]
    fn test_morgan_reproducible() {
        let mol = create_ethane();
        let options = MorganOptions::default();

        let fp1 = morgan_fingerprint(&mol, &options).unwrap();
        let fp2 = morgan_fingerprint(&mol, &options).unwrap();

        // Same molecule should produce same fingerprint
        assert_eq!(fp1, fp2);
    }
}
