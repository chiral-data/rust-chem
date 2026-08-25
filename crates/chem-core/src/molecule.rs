use crate::atom::Atom;
use crate::bond::Bond;
use crate::elements::ATOMIC_MASSES;
use crate::geometry::Point2;
use crate::graph::MoleculeGraph;

use std::{collections::HashMap, fmt};

use thiserror::Error;

/// Errors that can occur during molecule operations.
#[derive(Error, Debug)]
pub enum MoleculeError {
    #[error("Invalid atom index: {0}")]
    InvalidAtomIndex(usize),

    #[error("Invalid bond index: {0}")]
    InvalidBondIndex(usize),

    #[error("Bond references invalid atom indices")]
    InvalidBondAtoms,

    #[error("Molecule is empty")]
    EmptyMolecule,

    #[error("Invalid valence for atom {0}")]
    InvalidValence(usize),

    #[error("Expected {expected} coordinates (one per atom), got {got}")]
    CoordinateCountMismatch { expected: usize, got: usize },
}

/// Represents a complete molecule with atoms, bonds, and connectivity.
#[derive(Debug, Clone)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    graph: MoleculeGraph,
    name: Option<String>,
    properties: HashMap<String, String>,
    /// 2D coordinates for depiction, one per atom, indexed in parallel with
    /// `atoms`. `None` when the molecule has no layout — SMILES carries no
    /// coordinates, so this is the common case until a layout pass runs or
    /// the molecule came from a format (SDF) that supplies them.
    ///
    /// Stored per-molecule rather than on [`Atom`] for two reasons: `Atom`
    /// derives `Eq`, which floats can't satisfy, and coordinates are
    /// all-or-nothing in practice — there's no meaningful state where only
    /// some atoms have positions.
    coords: Option<Vec<Point2>>,
}

impl Molecule {
    pub fn new() -> Self {
        Molecule {
            atoms: Vec::new(),
            bonds: Vec::new(),
            graph: MoleculeGraph::new(0),
            name: None,
            properties: HashMap::new(),
            coords: None,
        }
    }

    pub fn with_capacity(num_atoms: usize, num_bonds: usize) -> Self {
        Molecule {
            atoms: Vec::with_capacity(num_atoms),
            bonds: Vec::with_capacity(num_bonds),
            graph: MoleculeGraph::new(0),
            name: None,
            properties: HashMap::new(),
            coords: None,
        }
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);

        // The new atom has no position, so any existing coordinate set is now
        // one short and no longer indexable in parallel with `atoms`. Drop it
        // rather than leave the two out of sync — the layout has to be redone
        // to place the new atom anyway.
        self.coords = None;

        let mut new_graph = MoleculeGraph::new(self.atoms.len());
        for (bond_idx, bond) in self.bonds.iter().enumerate() {
            new_graph.add_edge(bond.atom1(), bond.atom2(), bond_idx);
        }
        self.graph = new_graph;

        idx
    }

    pub fn add_bond(&mut self, bond: Bond) -> Result<usize, MoleculeError> {
        let (atom1, atom2) = bond.atoms();

        if atom1 >= self.atoms.len() || atom2 >= self.atoms.len() {
            return Err(MoleculeError::InvalidBondAtoms);
        }

        let bond_idx = self.bonds.len();
        self.graph.add_edge(atom1, atom2, bond_idx);
        self.bonds.push(bond);

        Ok(bond_idx)
    }

    pub fn atom(&self, idx: usize) -> &Atom {
        &self.atoms[idx]
    }

    pub fn atom_mut(&mut self, idx: usize) -> &mut Atom {
        &mut self.atoms[idx]
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn atoms_mut(&mut self) -> &mut [Atom] {
        &mut self.atoms
    }

    pub fn bond(&self, idx: usize) -> &Bond {
        &self.bonds[idx]
    }

    pub fn bond_mut(&mut self, idx: usize) -> &mut Bond {
        &mut self.bonds[idx]
    }

    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    pub fn graph(&self) -> &MoleculeGraph {
        &self.graph
    }

    pub fn degree(&self, atom_idx: usize) -> usize {
        self.graph.degree(atom_idx)
    }

    pub fn neighbors(&self, atom_idx: usize) -> &[crate::graph::Neighbor] {
        self.graph.neighbors(atom_idx)
    }

    /// This molecule's 2D coordinates, one per atom, or `None` if it has no
    /// layout. See the [`Molecule::coords`] field docs for why coordinates
    /// are all-or-nothing.
    pub fn coords(&self) -> Option<&[Point2]> {
        self.coords.as_deref()
    }

    /// The coordinate of a single atom, or `None` if this molecule has no
    /// layout or `atom_idx` is out of range.
    pub fn coord(&self, atom_idx: usize) -> Option<Point2> {
        self.coords.as_ref()?.get(atom_idx).copied()
    }

    pub fn has_coords(&self) -> bool {
        self.coords.is_some()
    }

    /// Sets this molecule's 2D coordinates.
    ///
    /// # Errors
    /// [`MoleculeError::CoordinateCountMismatch`] if `coords.len()` isn't
    /// exactly one per atom — coordinates are indexed in parallel with
    /// `atoms`, so a mismatched length would silently misattribute positions.
    pub fn set_coords(&mut self, coords: Vec<Point2>) -> Result<(), MoleculeError> {
        if coords.len() != self.atoms.len() {
            return Err(MoleculeError::CoordinateCountMismatch {
                expected: self.atoms.len(),
                got: coords.len(),
            });
        }
        self.coords = Some(coords);
        Ok(())
    }

    /// Discards any coordinates, e.g. to force a fresh layout.
    pub fn clear_coords(&mut self) {
        self.coords = None;
    }

    /// Mutable access to the coordinates, for a layout pass refining positions
    /// in place. `None` if this molecule has no layout yet — use
    /// [`Self::set_coords`] to establish one first.
    pub fn coords_mut(&mut self) -> Option<&mut [Point2]> {
        self.coords.as_deref_mut()
    }

    pub fn set_name(&mut self, name: String) {
        self.name = Some(name);
    }

    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    pub fn set_property(&mut self, key: String, value: String) {
        self.properties.insert(key, value);
    }

    pub fn property(&self, key: &str) -> Option<&str> {
        self.properties.get(key).map(|s| s.as_str())
    }

    pub fn properties(&self) -> &std::collections::HashMap<String, String> {
        &self.properties
    }

    pub fn calculate_implicit_hydrogens(&mut self) {
        for atom_idx in 0..self.atoms.len() {
            let atom = &self.atoms[atom_idx];

            // Skip if explicit H count is already set
            if atom.explicit_hydrogens() > 0 {
                continue;
            }

            let element = atom.element();
            let typical_valence = element.typical_valence();

            if typical_valence == 0 {
                continue;
            }

            let mut explicit_valence = 0.0_f64;
            for neighbor in self.graph.neighbors(atom_idx) {
                let bond = &self.bonds[neighbor.bond_idx];
                explicit_valence += bond.order().value();
            }

            let charge = atom.formal_charge();
            let adjusted_valence = (typical_valence as i16 - charge as i16) as u8;

            if (explicit_valence.round() as u8) < adjusted_valence {
                let implicit_h = adjusted_valence - explicit_valence.round() as u8;
                self.atoms[atom_idx].set_implicit_hydrogens(implicit_h);
            }
        }
    }

    pub fn formula(&self) -> String {
        let mut counts = std::collections::HashMap::new();

        for atom in &self.atoms {
            let symbol = atom.element().symbol();
            *counts.entry(symbol).or_insert(0) += 1;

            let h_count = atom.total_hydrogens();
            if h_count > 0 {
                *counts.entry("H").or_insert(0) += h_count as usize;
            }
        }

        let mut formula = String::new();

        if let Some(&count) = counts.get("C") {
            formula.push('C');
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }

        if let Some(&count) = counts.get("H") {
            formula.push('H');
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }

        let mut elements: Vec<_> = counts
            .iter()
            .filter(|&(k, _)| *k != "C" && *k != "H")
            .collect();
        elements.sort_by_key(|(k, _)| *k);

        for (element, &count) in elements {
            formula.push_str(element);
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }

        formula
    }

    pub fn molecular_weight(&self) -> f64 {
        let mut weight = 0.0;

        for atom in &self.atoms {
            weight += ATOMIC_MASSES[atom.atomic_number() as usize];
            weight += atom.total_hydrogens() as f64 * ATOMIC_MASSES[1];
        }

        weight
    }

    pub fn validate(&self) -> Result<(), MoleculeError> {
        if self.atoms.is_empty() {
            return Err(MoleculeError::EmptyMolecule);
        }

        for (idx, bond) in self.bonds.iter().enumerate() {
            let (a1, a2) = bond.atoms();
            if a1 >= self.atoms.len() || a2 >= self.atoms.len() {
                return Err(MoleculeError::InvalidBondIndex(idx));
            }
        }

        Ok(())
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(name) = &self.name {
            writeln!(f, "Molecule: {}", name)?;
        }
        writeln!(f, "Formula: {}", self.formula())?;
        writeln!(
            f,
            "Atoms: {}, Bonds: {}",
            self.num_atoms(),
            self.num_bonds()
        )?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Element;
    use crate::bond::BondOrder;

    #[test]
    fn test_molecule_creation() {
        let mol = Molecule::new();
        assert_eq!(mol.num_atoms(), 0);
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_add_atoms() {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new(Element::carbon()));
        let c2 = mol.add_atom(Atom::new(Element::carbon()));
        assert_eq!(c1, 0);
        assert_eq!(c2, 1);
        assert_eq!(mol.num_atoms(), 2);
    }

    #[test]
    fn test_add_bonds() {
        let mut mol = Molecule::new();
        let a1 = mol.add_atom(Atom::new(Element::carbon()));
        let a2 = mol.add_atom(Atom::new(Element::carbon()));
        let bond_idx = mol.add_bond(Bond::new(a1, a2, BondOrder::Single)).unwrap();
        assert_eq!(bond_idx, 0);
        assert_eq!(mol.num_bonds(), 1);
    }

    #[test]
    fn test_invalid_bond() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        let result = mol.add_bond(Bond::new(0, 5, BondOrder::Single));
        assert!(result.is_err());
    }

    #[test]
    fn test_formula() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.atoms_mut()[0].set_implicit_hydrogens(4);
        assert_eq!(mol.formula(), "CH4");
    }

    #[test]
    fn test_molecular_weight() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.atoms_mut()[0].set_implicit_hydrogens(2);
        let weight = mol.molecular_weight();
        assert!((weight - 18.016).abs() < 0.1);
    }

    #[test]
    fn test_neighbors() {
        let mut mol = Molecule::new();
        let a1 = mol.add_atom(Atom::new(Element::carbon()));
        let a2 = mol.add_atom(Atom::new(Element::carbon()));
        let a3 = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(a1, a2, BondOrder::Single)).unwrap();
        mol.add_bond(Bond::new(a1, a3, BondOrder::Single)).unwrap();
        assert_eq!(mol.degree(a1), 2);
        assert_eq!(mol.degree(a2), 1);
    }

    fn two_atom_molecule() -> Molecule {
        let mut mol = Molecule::new();
        let a1 = mol.add_atom(Atom::new(Element::carbon()));
        let a2 = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(a1, a2, BondOrder::Single)).unwrap();
        mol
    }

    #[test]
    fn test_no_coords_by_default() {
        let mol = two_atom_molecule();
        assert!(!mol.has_coords());
        assert!(mol.coords().is_none());
        assert!(mol.coord(0).is_none());
    }

    #[test]
    fn test_set_and_read_coords() {
        let mut mol = two_atom_molecule();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();

        assert!(mol.has_coords());
        assert_eq!(mol.coords().unwrap().len(), 2);
        assert_eq!(mol.coord(0), Some(Point2::new(0.0, 0.0)));
        assert_eq!(mol.coord(1), Some(Point2::new(1.5, 0.0)));
        // Out of range, even though the molecule does have a layout.
        assert!(mol.coord(2).is_none());
    }

    #[test]
    fn test_set_coords_rejects_wrong_count() {
        let mut mol = two_atom_molecule();

        let too_few = mol.set_coords(vec![Point2::new(0.0, 0.0)]);
        assert!(matches!(
            too_few,
            Err(MoleculeError::CoordinateCountMismatch {
                expected: 2,
                got: 1
            })
        ));

        let too_many = mol.set_coords(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
        ]);
        assert!(matches!(
            too_many,
            Err(MoleculeError::CoordinateCountMismatch {
                expected: 2,
                got: 3
            })
        ));

        // A rejected set leaves the molecule without coordinates rather than
        // half-applying them.
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_adding_atom_invalidates_coords() {
        let mut mol = two_atom_molecule();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();
        assert!(mol.has_coords());

        // The new atom has no position, so the coordinate set can no longer be
        // indexed in parallel with `atoms` -- it's dropped rather than left
        // one short.
        mol.add_atom(Atom::new(Element::oxygen()));
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_adding_bond_keeps_coords() {
        let mut mol = two_atom_molecule();
        let a3 = mol.add_atom(Atom::new(Element::carbon()));
        mol.set_coords(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.5, 0.0),
            Point2::new(3.0, 0.0),
        ])
        .unwrap();

        // A new bond doesn't change the atom count, so coordinates stay
        // structurally valid (the geometry may be less ideal, but it's still
        // one position per atom).
        mol.add_bond(Bond::new(0, a3, BondOrder::Single)).unwrap();
        assert!(mol.has_coords());
        assert_eq!(mol.coord(2), Some(Point2::new(3.0, 0.0)));
    }

    #[test]
    fn test_clear_coords() {
        let mut mol = two_atom_molecule();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();
        mol.clear_coords();
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_coords_mut() {
        let mut mol = two_atom_molecule();
        assert!(mol.coords_mut().is_none());

        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();
        for p in mol.coords_mut().unwrap() {
            p.x += 10.0;
        }
        assert_eq!(mol.coord(0), Some(Point2::new(10.0, 0.0)));
        assert_eq!(mol.coord(1), Some(Point2::new(11.5, 0.0)));
    }

    #[test]
    fn test_coords_survive_clone() {
        let mut mol = two_atom_molecule();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();

        let cloned = mol.clone();
        assert_eq!(cloned.coord(1), Some(Point2::new(1.5, 0.0)));
    }
}
