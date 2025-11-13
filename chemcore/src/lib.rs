pub mod atom;
pub mod bond;
pub mod graph;
pub mod molecule;

mod elements;

pub mod prelude {
    pub use crate::atom::{Atom, Chirality, Element, Hybridization};
    pub use crate::bond::{Bond, BondOrder, BondStereo, BondType};
    pub use crate::graph::{MoleculeGraph, Neighbor};
    pub use crate::molecule::{Molecule, MoleculeError};
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn test_basic_workflow() {
        let mut mol = Molecule::new();
        let c = mol.add_atom(Atom::new(Element::carbon()));
        mol.atom_mut(c).set_implicit_hydrogens(4);
        assert_eq!(mol.num_atoms(), 1);
        assert_eq!(mol.formula(), "CH4");
    }

    #[test]
    fn test_ethane() {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new(Element::carbon()));
        let c2 = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(c1, c2, BondOrder::Single)).unwrap();
        mol.calculate_implicit_hydrogens();
        assert_eq!(mol.formula(), "C2H6");
        assert_eq!(mol.num_bonds(), 1);
    }
}
