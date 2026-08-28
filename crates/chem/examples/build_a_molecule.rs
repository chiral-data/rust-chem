//! Building a molecule by hand, without parsing anything.
//!
//! ```sh
//! cargo run --example build_a_molecule
//! ```
//!
//! Most of the time you want [`chem::io::smiles::parse_smiles`]. This is the
//! layer underneath it — useful when molecules arrive from somewhere that is
//! not a file, and the fastest way to see what a `Molecule` actually holds.

use chem::core::prelude::*;

fn main() -> Result<(), MoleculeError> {
    // Ethanol, CCO. Atoms first: `add_atom` returns the index you refer to
    // them by afterwards.
    let mut mol = Molecule::new();
    let c1 = mol.add_atom(Atom::new(Element::carbon()));
    let c2 = mol.add_atom(Atom::new(Element::carbon()));
    let o = mol.add_atom(Atom::new(Element::oxygen()));

    // Then bonds, by index. `add_bond` validates that both atoms exist, which
    // is why it returns a Result.
    mol.add_bond(Bond::new(c1, c2, BondOrder::Single))?;
    mol.add_bond(Bond::new(c2, o, BondOrder::Single))?;

    // Hydrogens are implicit: they are counted from each atom's valence rather
    // than stored, so a molecule holds only its heavy atoms.
    mol.calculate_implicit_hydrogens();

    println!("{} heavy atoms, {} bonds", mol.num_atoms(), mol.num_bonds());
    for i in 0..mol.num_atoms() {
        let atom = mol.atom(i);
        println!(
            "  {i}: {:<2} {} implicit hydrogens, degree {}",
            ELEMENT_SYMBOLS[atom.atomic_number() as usize],
            atom.implicit_hydrogens(),
            mol.degree(i),
        );
    }

    // 2D coordinates are not intrinsic to a molecule — they are computed, and
    // only when something needs to draw it.
    println!("\nhas coordinates: {}", mol.has_coords());
    ensure_coords(&mut mol);
    println!("after ensure_coords: {}", mol.has_coords());
    if let Some(points) = mol.coords() {
        for (i, p) in points.iter().enumerate() {
            println!("  {i}: ({:.3}, {:.3})", p.x, p.y);
        }
    }

    Ok(())
}
