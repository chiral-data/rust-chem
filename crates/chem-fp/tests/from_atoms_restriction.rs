use chem_core::atom::{Atom, Element};
use chem_core::bond::{Bond, BondOrder};
use chem_core::molecule::Molecule;
use chem_fp::morgan::MorganFingerprint;

/// Regression test for issue #13: `include_atoms[i] && !only_nonzero_invariants
/// || current_invariants[i] != 0` parsed as `(include_atoms[i] && !only_nonzero)
/// || (invariant != 0)` instead of `include_atoms[i] && (!only_nonzero || invariant
/// != 0)`, so the `from_atoms` restriction was silently ignored for round-0
/// environments whenever an atom's invariant was nonzero (i.e. always, for real
/// atoms).
#[test]
fn test_from_atoms_restricts_round_zero_environments() {
    // Two atoms with different elements -> different round-0 invariants, so
    // (with a large enough fp_size) they set different bits.
    let mut mol = Molecule::new();
    let carbon = mol.add_atom(Atom::new(Element::carbon()));
    let oxygen = mol.add_atom(Atom::new(Element::oxygen()));
    mol.add_bond(Bond::new(carbon, oxygen, BondOrder::Single))
        .unwrap();
    mol.calculate_implicit_hydrogens();

    let fp_size = 2048;

    let fp_all = MorganFingerprint::get_fingerprint_as_bitvec(
        &mol, 0, fp_size, None, None, false, true, false,
    )
    .unwrap();

    let fp_carbon_only = MorganFingerprint::get_fingerprint_as_bitvec(
        &mol,
        0,
        fp_size,
        None,
        Some(&[carbon as u32]),
        false,
        true,
        false,
    )
    .unwrap();

    let bits_all: Vec<usize> = fp_all.iter_ones().collect();
    let bits_carbon_only: Vec<usize> = fp_carbon_only.iter_ones().collect();

    // Different elements must produce different round-0 bits, or this test
    // can't distinguish "restricted" from "unrestricted".
    assert_eq!(
        bits_all.len(),
        2,
        "expected carbon and oxygen to set distinct round-0 bits, got {:?}",
        bits_all
    );

    assert_eq!(
        bits_carbon_only.len(),
        1,
        "from_atoms=[carbon] must restrict to a single round-0 bit, got {:?} \
         (full set was {:?}) — from_atoms gate is not restricting environments",
        bits_carbon_only,
        bits_all
    );
    assert!(
        bits_all.contains(&bits_carbon_only[0]),
        "the single restricted bit should be one of the two full-set bits"
    );
}
