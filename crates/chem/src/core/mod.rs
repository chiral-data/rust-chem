pub mod atom;
pub mod bond;
pub mod cell;
pub mod geometry;
pub mod graph;
pub mod layout;
pub mod molecule;
pub mod residue;
pub mod rings;
pub mod site;
pub mod stereo;

pub mod elements;

pub mod prelude {
    pub use crate::core::atom::{Atom, Chirality, Element, Hybridization};
    pub use crate::core::bond::{Bond, BondOrder, BondStereo, BondType};
    pub use crate::core::cell::{SpaceGroup, UnitCell};
    pub use crate::core::elements::{ATOMIC_MASSES, ELEMENT_NAMES, ELEMENT_SYMBOLS};
    pub use crate::core::geometry::{BoundingBox, Point2, Point3};
    pub use crate::core::graph::{MoleculeGraph, Neighbor};
    pub use crate::core::layout::{ensure_coords, layout};
    pub use crate::core::molecule::{Molecule, MoleculeError};
    pub use crate::core::residue::{Chain, Residue};
    pub use crate::core::rings::{Ring, find_sssr, perceive_rings, ring_count};
    pub use crate::core::site::AtomSite;
    pub use crate::core::stereo::perceive_bond_stereo;
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
