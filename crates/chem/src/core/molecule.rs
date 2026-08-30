use crate::core::atom::Atom;
use crate::core::bond::Bond;
use crate::core::elements::ATOMIC_MASSES;
use crate::core::geometry::{Point2, Point3};
use crate::core::graph::MoleculeGraph;
use crate::core::residue::{Chain, Residue};
use crate::core::site::AtomSite;

use std::{collections::HashMap, fmt};

use thiserror::Error;

/// Errors that can occur during molecule operations.
#[derive(Error, Debug)]
/// Adding a variant to a public enum is a breaking change unless callers
/// are told not to match it exhaustively. This is the attribute that says
/// so, and it has to be present from the first published version: adding it
/// later invalidates every exhaustive match written against the earlier one.
#[non_exhaustive]
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

    #[error("Expected {expected} atom sites (one per atom), got {got}")]
    SiteCountMismatch { expected: usize, got: usize },

    /// A single variant rather than one per failure mode: topology has five
    /// distinct ways to be wrong and structured variants for each would be
    /// permanent API bought for very little. The message names the offending
    /// index and the rule it broke.
    #[error("Invalid residue topology: {0}")]
    InvalidTopology(String),
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
    /// 3D coordinates — a conformer, one per atom, indexed in parallel with
    /// `atoms`. `None` when the molecule has no geometry, which is the common
    /// case: SMILES carries none, and a flat SDF carries a drawing rather than
    /// a conformer.
    ///
    /// Kept alongside `coords` rather than replacing it because the two are
    /// different artefacts and neither is derivable from the other. A layout
    /// is computed for drawing; a conformer is physical. A molecule read from
    /// a 3D file and then laid out for the structure view legitimately has
    /// both, and depiction wants the first while a conversion to XYZ wants the
    /// second.
    coords3: Option<Vec<Point3>>,
    /// Per-atom file data — names, partial charges, occupancies, temperature
    /// factors — one per atom, indexed in parallel with `atoms`. `None` for
    /// every molecule that did not come from a format carrying such columns,
    /// which is most of them.
    ///
    /// A third side table rather than fields on [`Atom`] for the reason the
    /// coordinates give: `Atom` derives `Eq`, and these are floats. It is also
    /// the honest shape — an element is a fact about a species, while a
    /// B-factor is a fact about one observation of it, and two files of the
    /// same molecule will disagree about the second while agreeing on the
    /// first.
    sites: Option<Vec<AtomSite>>,
    /// Chains, in file order. Empty for every molecule that did not come from a
    /// format organised by chain, which is most of them.
    ///
    /// A plain `Vec` rather than an `Option<Vec>` like the three tables above:
    /// chains are not indexed in parallel with `atoms`, so "empty" is already
    /// the natural absent state and an `Option` would add noise for nothing.
    chains: Vec<Chain>,
    /// Residues, in file order, each owning a contiguous range of atoms.
    ///
    /// Ranges rather than a per-atom residue index because PDB and mmCIF
    /// already order their records that way and it costs nothing per atom. The
    /// price is that a residue's atoms must be contiguous and the ranges must
    /// ascend — the latter is what makes [`Molecule::residue_of`]'s binary
    /// search correct, so both are enforced by [`Molecule::set_topology`]
    /// rather than trusted.
    residues: Vec<Residue>,
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
            coords3: None,
            sites: None,
            chains: Vec::new(),
            residues: Vec::new(),
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
            coords3: None,
            sites: None,
            chains: Vec::new(),
            residues: Vec::new(),
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
        // to place the new atom anyway. The conformer goes for the same
        // reason, and there is no way to invent a position for the new atom.
        // The site table is dropped on the same grounds: it is indexed in
        // parallel with `atoms` too, and a file supplied it for a set of atoms
        // this is no longer.
        self.coords = None;
        self.coords3 = None;
        self.sites = None;

        // Topology goes too. Appending does not strictly invalidate the
        // existing ranges — indices 0..n keep their meaning — but it produces
        // an atom belonging to no residue, which a PDB write would silently
        // drop. One rule for everything is easier to reason about than an
        // exception here: what a file said about a set of atoms stops applying
        // when that is no longer the set.
        self.chains.clear();
        self.residues.clear();

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

    pub fn neighbors(&self, atom_idx: usize) -> &[crate::core::graph::Neighbor] {
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

    /// This molecule's 3D coordinates, one per atom, or `None` if it carries
    /// no conformer. Independent of [`Self::coords`]: a molecule may have
    /// either, both or neither.
    pub fn coords3(&self) -> Option<&[Point3]> {
        self.coords3.as_deref()
    }

    /// The 3D coordinate of a single atom, or `None` if this molecule has no
    /// conformer or `atom_idx` is out of range.
    pub fn coord3(&self, atom_idx: usize) -> Option<Point3> {
        self.coords3.as_ref()?.get(atom_idx).copied()
    }

    pub fn has_coords3(&self) -> bool {
        self.coords3.is_some()
    }

    /// Sets this molecule's 3D coordinates.
    ///
    /// # Errors
    /// [`MoleculeError::CoordinateCountMismatch`] if `coords.len()` isn't
    /// exactly one per atom — the same all-or-nothing rule
    /// [`Self::set_coords`] enforces, and for the same reason.
    pub fn set_coords3(&mut self, coords: Vec<Point3>) -> Result<(), MoleculeError> {
        if coords.len() != self.atoms.len() {
            return Err(MoleculeError::CoordinateCountMismatch {
                expected: self.atoms.len(),
                got: coords.len(),
            });
        }
        self.coords3 = Some(coords);
        Ok(())
    }

    /// Discards any conformer, leaving a 2D layout untouched.
    pub fn clear_coords3(&mut self) {
        self.coords3 = None;
    }

    /// Mutable access to the conformer, for a pass refining positions in
    /// place. `None` if this molecule has no conformer yet — use
    /// [`Self::set_coords3`] to establish one first.
    pub fn coords3_mut(&mut self) -> Option<&mut [Point3]> {
        self.coords3.as_deref_mut()
    }

    /// This molecule's per-atom file data, one entry per atom, or `None` if it
    /// carries none. Independent of both coordinate sets.
    pub fn sites(&self) -> Option<&[AtomSite]> {
        self.sites.as_deref()
    }

    /// The site record for a single atom, or `None` if this molecule has no
    /// site data or `atom_idx` is out of range.
    ///
    /// Returns a reference rather than a copy, unlike [`Self::coord3`]:
    /// `AtomSite` holds a `String` and is not `Copy`.
    pub fn site(&self, atom_idx: usize) -> Option<&AtomSite> {
        self.sites.as_ref()?.get(atom_idx)
    }

    pub fn has_sites(&self) -> bool {
        self.sites.is_some()
    }

    /// Sets this molecule's per-atom file data.
    ///
    /// # Errors
    /// [`MoleculeError::SiteCountMismatch`] if `sites.len()` isn't exactly one
    /// per atom — the same all-or-nothing rule [`Self::set_coords`] enforces,
    /// and for the same reason: the table is indexed in parallel with `atoms`,
    /// so a mismatched length would silently misattribute a charge or a
    /// B-factor to the wrong atom.
    pub fn set_sites(&mut self, sites: Vec<AtomSite>) -> Result<(), MoleculeError> {
        if sites.len() != self.atoms.len() {
            return Err(MoleculeError::SiteCountMismatch {
                expected: self.atoms.len(),
                got: sites.len(),
            });
        }
        self.sites = Some(sites);
        Ok(())
    }

    /// Discards any per-atom file data, leaving both coordinate sets untouched.
    pub fn clear_sites(&mut self) {
        self.sites = None;
    }

    /// Mutable access to the site table, for a pass filling in a column the
    /// reader could not. `None` if this molecule has no site data yet — use
    /// [`Self::set_sites`] to establish it first.
    pub fn sites_mut(&mut self) -> Option<&mut [AtomSite]> {
        self.sites.as_deref_mut()
    }

    /// This molecule's chains, in file order. Empty if it carries no topology.
    pub fn chains(&self) -> &[Chain] {
        &self.chains
    }

    /// This molecule's residues, in file order, each owning a contiguous and
    /// ascending range of atoms.
    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }

    pub fn has_topology(&self) -> bool {
        !self.residues.is_empty()
    }

    /// Sets this molecule's chain and residue topology.
    ///
    /// The ranges are the only thing tying residues to atoms, so they are
    /// validated rather than trusted — a bad range would silently attribute an
    /// atom to the wrong residue, which is exactly the class of error this
    /// structure exists to prevent.
    ///
    /// Ranges need not cover every atom: a ligand appended without residue
    /// information is a legal molecule, and [`Self::residue_of`] answers `None`
    /// for it.
    ///
    /// # Errors
    /// [`MoleculeError::InvalidTopology`], naming the offending index, if an
    /// atom range runs past the last atom, if residue ranges overlap or do not
    /// ascend, if a chain's residue range is out of bounds or does not ascend,
    /// or if a residue's `chain_ix` disagrees with the chain that claims it.
    pub fn set_topology(
        &mut self,
        chains: Vec<Chain>,
        residues: Vec<Residue>,
    ) -> Result<(), MoleculeError> {
        let invalid = |msg: String| Err(MoleculeError::InvalidTopology(msg));

        // Residue atom ranges: in bounds, ascending, non-overlapping. Ascending
        // is not cosmetic — `residue_of` binary-searches these.
        let mut previous_end = 0usize;
        for (ix, residue) in residues.iter().enumerate() {
            if residue.atoms.start > residue.atoms.end {
                return invalid(format!(
                    "residue {ix} has a backwards atom range {}..{}",
                    residue.atoms.start, residue.atoms.end
                ));
            }
            if residue.atoms.end > self.atoms.len() {
                return invalid(format!(
                    "residue {ix} covers atoms {}..{} but the molecule has {}",
                    residue.atoms.start,
                    residue.atoms.end,
                    self.atoms.len()
                ));
            }
            if ix > 0 && residue.atoms.start < previous_end {
                return invalid(format!(
                    "residue {ix} starts at atom {} but residue {} already ends at {previous_end}; \
                     residue atom ranges must ascend and must not overlap",
                    residue.atoms.start,
                    ix - 1
                ));
            }
            previous_end = residue.atoms.end;

            if residue.chain_ix >= chains.len() {
                return invalid(format!(
                    "residue {ix} names chain {} but there are {} chains",
                    residue.chain_ix,
                    chains.len()
                ));
            }
        }

        // Chain residue ranges: same three rules, one level up.
        let mut previous_end = 0usize;
        for (ix, chain) in chains.iter().enumerate() {
            if chain.residues.start > chain.residues.end {
                return invalid(format!(
                    "chain {ix} has a backwards residue range {}..{}",
                    chain.residues.start, chain.residues.end
                ));
            }
            if chain.residues.end > residues.len() {
                return invalid(format!(
                    "chain {ix} covers residues {}..{} but there are {}",
                    chain.residues.start,
                    chain.residues.end,
                    residues.len()
                ));
            }
            if ix > 0 && chain.residues.start < previous_end {
                return invalid(format!(
                    "chain {ix} starts at residue {} but chain {} already ends at {previous_end}; \
                     chain residue ranges must ascend and must not overlap",
                    chain.residues.start,
                    ix - 1
                ));
            }
            previous_end = chain.residues.end;

            // The two directions have to agree, or a writer grouping by chain
            // and a reader following `chain_ix` would disagree about the same
            // structure.
            for residue_ix in chain.residues.clone() {
                if residues[residue_ix].chain_ix != ix {
                    return invalid(format!(
                        "chain {ix} claims residue {residue_ix}, but that residue names chain {}",
                        residues[residue_ix].chain_ix
                    ));
                }
            }
        }

        self.chains = chains;
        self.residues = residues;
        Ok(())
    }

    /// Discards chain and residue topology, leaving every per-atom table alone.
    pub fn clear_topology(&mut self) {
        self.chains.clear();
        self.residues.clear();
    }

    /// The residue owning `atom_idx`, or `None` if this molecule has no
    /// topology or the atom belongs to no residue.
    ///
    /// `O(log n)` — a binary search over the ascending ranges, which is the
    /// cost of storing ranges instead of a residue index on every atom.
    pub fn residue_of(&self, atom_idx: usize) -> Option<&Residue> {
        let found = self
            .residues
            .partition_point(|residue| residue.atoms.end <= atom_idx);
        self.residues
            .get(found)
            .filter(|residue| residue.contains_atom(atom_idx))
    }

    /// The chain owning `atom_idx`, or `None` on the same terms as
    /// [`Self::residue_of`].
    pub fn chain_of(&self, atom_idx: usize) -> Option<&Chain> {
        let residue = self.residue_of(atom_idx)?;
        self.chains.get(residue.chain_ix)
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
    use crate::core::atom::Element;
    use crate::core::bond::BondOrder;
    use std::ops::Range;

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
    fn test_conformer_is_independent_of_the_layout() {
        // Both, at once, holding different values. A molecule read from a 3D
        // file and then laid out for drawing is exactly this state, and
        // neither set may overwrite the other.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::carbon()));

        assert!(!mol.has_coords3());
        assert!(mol.coords3().is_none());

        mol.set_coords3(vec![
            Point3::new(0.0, 0.0, -1.0),
            Point3::new(1.0, 0.0, 1.0),
        ])
        .unwrap();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();

        assert_eq!(mol.coord3(1), Some(Point3::new(1.0, 0.0, 1.0)));
        assert_eq!(mol.coord(1), Some(Point2::new(1.5, 0.0)));

        // Clearing one leaves the other standing.
        mol.clear_coords3();
        assert!(!mol.has_coords3());
        assert!(mol.has_coords());
    }

    #[test]
    fn test_set_coords3_rejects_wrong_count() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::carbon()));

        let too_few = mol.set_coords3(vec![Point3::ORIGIN]);
        assert!(matches!(
            too_few,
            Err(MoleculeError::CoordinateCountMismatch {
                expected: 2,
                got: 1
            })
        ));
        assert!(!mol.has_coords3());
    }

    #[test]
    fn test_sites_are_independent_of_both_coordinate_sets() {
        // The state a PDB reader produces: geometry from the file, site data
        // from the same file, and a layout computed later for drawing. All
        // three are indexed in parallel with `atoms` and none may disturb the
        // others.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));

        assert!(!mol.has_sites());
        assert!(mol.sites().is_none());
        assert!(mol.site(0).is_none());

        mol.set_coords3(vec![
            Point3::new(0.0, 0.0, -1.0),
            Point3::new(1.0, 0.0, 1.0),
        ])
        .unwrap();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.5, 0.0)])
            .unwrap();
        mol.set_sites(vec![
            AtomSite {
                name: Some("CA".to_string()),
                occupancy: Some(1.0),
                b_factor: Some(23.45),
                ..AtomSite::default()
            },
            AtomSite {
                name: Some("OD1".to_string()),
                partial_charge: Some(-0.4157),
                ..AtomSite::default()
            },
        ])
        .unwrap();

        assert_eq!(mol.site(0).unwrap().name.as_deref(), Some("CA"));
        assert_eq!(mol.site(1).unwrap().partial_charge, Some(-0.4157));
        assert_eq!(mol.coord3(1), Some(Point3::new(1.0, 0.0, 1.0)));
        assert_eq!(mol.coord(1), Some(Point2::new(1.5, 0.0)));

        // Clearing one leaves the other two standing.
        mol.clear_sites();
        assert!(!mol.has_sites());
        assert!(mol.has_coords());
        assert!(mol.has_coords3());
    }

    #[test]
    fn test_a_b_factor_survives_untouched() {
        // Pinned deliberately. OpenBabel zeroes this column on every PDB write,
        // and a predicted structure reuses it for per-atom confidence — so the
        // value has to come back exactly, and being different from that
        // behaviour is the point rather than an accident.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.set_sites(vec![AtomSite {
            b_factor: Some(87.31),
            ..AtomSite::default()
        }])
        .unwrap();

        let cloned = mol.clone();
        assert_eq!(cloned.site(0).unwrap().b_factor, Some(87.31));

        if let Some(sites) = mol.sites_mut() {
            sites[0].occupancy = Some(0.5);
        }
        assert_eq!(mol.site(0).unwrap().b_factor, Some(87.31));
        assert_eq!(mol.site(0).unwrap().occupancy, Some(0.5));
    }

    #[test]
    fn test_set_sites_rejects_wrong_count() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::carbon()));

        let too_few = mol.set_sites(vec![AtomSite::empty()]);
        assert!(matches!(
            too_few,
            Err(MoleculeError::SiteCountMismatch {
                expected: 2,
                got: 1
            })
        ));
        assert!(!mol.has_sites());

        let too_many = mol.set_sites(vec![
            AtomSite::empty(),
            AtomSite::empty(),
            AtomSite::empty(),
        ]);
        assert!(matches!(
            too_many,
            Err(MoleculeError::SiteCountMismatch {
                expected: 2,
                got: 3
            })
        ));
        assert!(!mol.has_sites());
    }

    // ---- residue and chain topology ------------------------------------

    /// A two-chain structure: chain A with LYS 58 (3 atoms) and GLY 59 (2),
    /// chain B with its own LYS 58 (2), then a trailing water in no residue.
    fn two_chain_molecule() -> Molecule {
        let mut mol = Molecule::new();
        for _ in 0..8 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        mol
    }

    fn residue(name: &str, sequence: i32, chain_ix: usize, atoms: Range<usize>) -> Residue {
        Residue {
            name: name.to_string(),
            sequence,
            insertion_code: None,
            label_seq: None,
            chain_ix,
            is_hetero: false,
            atoms,
        }
    }

    fn chain(id: &str, residues: Range<usize>) -> Chain {
        Chain {
            id: id.to_string(),
            label_id: None,
            residues,
        }
    }

    fn two_chain_topology() -> (Vec<Chain>, Vec<Residue>) {
        (
            vec![chain("A", 0..2), chain("B", 2..3)],
            vec![
                residue("LYS", 58, 0, 0..3),
                residue("GLY", 59, 0, 3..5),
                residue("LYS", 58, 1, 5..7),
            ],
        )
    }

    #[test]
    fn test_the_same_residue_number_in_two_chains_stays_distinguishable() {
        // The case a flat atom list cannot represent, and the reason topology
        // exists: LYS 58 of chain A and LYS 58 of chain B are different
        // residues, and a structure that conflates them is wrong in a way no
        // atom-level check would catch.
        let mut mol = two_chain_molecule();
        let (chains, residues) = two_chain_topology();
        mol.set_topology(chains, residues).unwrap();

        assert!(mol.has_topology());
        assert_eq!(mol.chains().len(), 2);
        assert_eq!(mol.residues().len(), 3);

        let a58 = mol.residue_of(0).unwrap();
        let b58 = mol.residue_of(5).unwrap();
        assert_eq!(a58.to_string(), "LYS 58");
        assert_eq!(b58.to_string(), "LYS 58");
        assert_ne!(a58.chain_ix, b58.chain_ix);
        assert_eq!(mol.chain_of(0).unwrap().id, "A");
        assert_eq!(mol.chain_of(5).unwrap().id, "B");
    }

    #[test]
    fn test_residue_of_resolves_every_atom_and_admits_none() {
        let mut mol = two_chain_molecule();
        let (chains, residues) = two_chain_topology();
        mol.set_topology(chains, residues).unwrap();

        for atom in 0..3 {
            assert_eq!(mol.residue_of(atom).unwrap().name, "LYS");
        }
        for atom in 3..5 {
            assert_eq!(mol.residue_of(atom).unwrap().name, "GLY");
        }
        assert_eq!(mol.residue_of(6).unwrap().sequence, 58);

        // Atoms 7 is past every range — a ligand appended without residue
        // information is a legal molecule, not an error.
        assert!(mol.residue_of(7).is_none());
        assert!(mol.chain_of(7).is_none());
        assert!(mol.residue_of(999).is_none());
    }

    #[test]
    fn test_insertion_codes_make_three_residues_from_one_number() {
        // 58, 58A and 58B are distinct residues. Without the code, a structure
        // with an insertion collapses to a single residue and the extra atoms
        // are silently reattributed.
        let mut mol = Molecule::new();
        for _ in 0..3 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        let residues = vec![
            residue("LYS", 58, 0, 0..1),
            Residue {
                insertion_code: Some('A'),
                ..residue("SER", 58, 0, 1..2)
            },
            Residue {
                insertion_code: Some('B'),
                ..residue("THR", 58, 0, 2..3)
            },
        ];
        mol.set_topology(vec![chain("A", 0..3)], residues).unwrap();

        let rendered: Vec<String> = mol.residues().iter().map(|r| r.to_string()).collect();
        assert_eq!(rendered, vec!["LYS 58", "SER 58A", "THR 58B"]);
        assert_ne!(mol.residues()[0], mol.residues()[1]);
    }

    #[test]
    fn test_both_numbering_schemes_are_kept() {
        // mmCIF's label_* and auth_* frequently disagree and neither is
        // derivable from the other, so storing one renumbers the structure on
        // write. A water carries label_seq: None, which is what mmCIF itself
        // writes for a non-polymer.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));

        let residues = vec![
            Residue {
                label_seq: Some(12),
                ..residue("LYS", 58, 0, 0..1)
            },
            Residue {
                is_hetero: true,
                label_seq: None,
                ..residue("HOH", 401, 0, 1..2)
            },
        ];
        mol.set_topology(
            vec![Chain {
                label_id: Some("C".to_string()),
                ..chain("A", 0..2)
            }],
            residues,
        )
        .unwrap();

        let lys = &mol.residues()[0];
        assert_eq!((lys.sequence, lys.label_seq), (58, Some(12)));
        let water = &mol.residues()[1];
        assert_eq!((water.sequence, water.label_seq), (401, None));
        assert!(water.is_hetero && !lys.is_hetero);
        assert_eq!(mol.chains()[0].id, "A");
        assert_eq!(mol.chains()[0].label_id.as_deref(), Some("C"));
    }

    #[test]
    fn test_set_topology_rejects_every_broken_shape() {
        // The ranges are the only thing tying residues to atoms, so each of
        // these would silently attribute an atom to the wrong residue.
        let base = two_chain_molecule;

        // 1. an atom range past the last atom
        let mut mol = base();
        let err = mol
            .set_topology(vec![chain("A", 0..1)], vec![residue("LYS", 1, 0, 0..99)])
            .unwrap_err();
        assert!(format!("{err}").contains("but the molecule has 8"), "{err}");
        assert!(!mol.has_topology());

        // 2. overlapping residue ranges
        let mut mol = base();
        let err = mol
            .set_topology(
                vec![chain("A", 0..2)],
                vec![residue("LYS", 1, 0, 0..4), residue("GLY", 2, 0, 2..6)],
            )
            .unwrap_err();
        assert!(format!("{err}").contains("must ascend"), "{err}");
        assert!(!mol.has_topology());

        // 3. a residue naming a chain that does not exist
        let mut mol = base();
        let err = mol
            .set_topology(vec![chain("A", 0..1)], vec![residue("LYS", 1, 7, 0..3)])
            .unwrap_err();
        assert!(format!("{err}").contains("there are 1 chains"), "{err}");
        assert!(!mol.has_topology());

        // 4. a chain claiming more residues than exist
        let mut mol = base();
        let err = mol
            .set_topology(vec![chain("A", 0..5)], vec![residue("LYS", 1, 0, 0..3)])
            .unwrap_err();
        assert!(format!("{err}").contains("but there are 1"), "{err}");
        assert!(!mol.has_topology());

        // 5. the two directions disagreeing about who owns a residue
        let mut mol = base();
        let err = mol
            .set_topology(
                vec![chain("A", 0..1), chain("B", 1..2)],
                vec![residue("LYS", 1, 0, 0..3), residue("GLY", 2, 0, 3..5)],
            )
            .unwrap_err();
        assert!(format!("{err}").contains("names chain 0"), "{err}");
        assert!(!mol.has_topology());
    }

    /// Many residues, arranged to include every awkward shape: a gap before
    /// the first, single-atom residues, gaps between residues, several chains,
    /// and a tail of atoms past the last residue.
    fn gappy_structure() -> Molecule {
        const RESIDUES: usize = 200;
        const CHAINS: usize = 5;
        const PER_CHAIN: usize = RESIDUES / CHAINS;

        let mut residues = Vec::with_capacity(RESIDUES);
        let mut cursor = 2; // leading gap: atoms 0 and 1 belong to nothing
        for ix in 0..RESIDUES {
            // Alternating sizes, so single-atom residues are covered.
            let size = if ix % 2 == 0 { 1 } else { 3 };
            residues.push(residue(
                "RES",
                ix as i32,
                ix / PER_CHAIN,
                cursor..cursor + size,
            ));
            cursor += size;
            // A gap after every fifth residue.
            if ix % 5 == 4 {
                cursor += 2;
            }
        }
        let chains: Vec<Chain> = (0..CHAINS)
            .map(|c| {
                chain(
                    "ABCDE".get(c..c + 1).unwrap(),
                    c * PER_CHAIN..(c + 1) * PER_CHAIN,
                )
            })
            .collect();

        let mut mol = Molecule::new();
        for _ in 0..cursor + 3 {
            // trailing atoms past the last residue
            mol.add_atom(Atom::new(Element::carbon()));
        }
        mol.set_topology(chains, residues).unwrap();
        mol
    }

    #[test]
    fn test_residue_of_agrees_with_a_linear_scan_at_every_atom() {
        // `residue_of` is a binary search over the ascending ranges, and a
        // wrong predicate is an off-by-one that shows only at a boundary.
        // Rather than guess which boundaries matter, check every atom against
        // the obvious slow answer over a deliberately awkward layout. Size is
        // not the point here — coverage of the boundaries is.
        let mol = gappy_structure();
        assert_eq!(mol.residues().len(), 200);
        assert_eq!(mol.chains().len(), 5);

        let mut hits = 0;
        let mut misses = 0;
        for atom in 0..mol.num_atoms() {
            let expected = mol.residues().iter().find(|r| r.atoms.contains(&atom));
            assert_eq!(mol.residue_of(atom), expected, "residue_of({atom})");

            match expected {
                Some(res) => {
                    hits += 1;
                    assert_eq!(
                        mol.chain_of(atom),
                        mol.chains().get(res.chain_ix),
                        "chain_of({atom})"
                    );
                }
                None => {
                    misses += 1;
                    assert!(mol.chain_of(atom).is_none(), "chain_of({atom})");
                }
            }
        }

        // The fixture has to actually contain both kinds of atom, or this
        // whole test could pass vacuously.
        assert!(
            hits > 0 && misses > 0,
            "{hits} in residues, {misses} in gaps"
        );
        assert!(mol.residue_of(mol.num_atoms()).is_none());
        assert!(mol.residue_of(usize::MAX).is_none());
    }

    #[test]
    fn test_an_atom_in_a_gap_answers_none_rather_than_its_neighbour() {
        // The `filter` after the binary search is what guarantees this: if the
        // search lands on the wrong residue, membership rejects it and the
        // answer is a miss rather than a confidently wrong residue. Anyone
        // optimising that filter away would turn misses into wrong answers,
        // so the property is pinned rather than left implicit.
        let mut mol = Molecule::new();
        for _ in 0..6 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        // Atoms 2 and 3 are in neither residue.
        mol.set_topology(
            vec![chain("A", 0..2)],
            vec![residue("LYS", 1, 0, 0..2), residue("GLY", 2, 0, 4..6)],
        )
        .unwrap();

        assert_eq!(mol.residue_of(1).unwrap().name, "LYS");
        assert!(mol.residue_of(2).is_none());
        assert!(mol.residue_of(3).is_none());
        assert_eq!(mol.residue_of(4).unwrap().name, "GLY");
    }

    #[test]
    fn test_clear_topology_leaves_the_per_atom_tables_alone() {
        let mut mol = two_chain_molecule();
        let (chains, residues) = two_chain_topology();
        mol.set_topology(chains, residues).unwrap();
        mol.set_coords(vec![Point2::ORIGIN; 8]).unwrap();

        mol.clear_topology();
        assert!(!mol.has_topology());
        assert!(mol.residues().is_empty() && mol.chains().is_empty());
        assert!(mol.has_coords());
    }

    #[test]
    fn test_adding_an_atom_drops_everything_a_file_supplied() {
        // The conformer goes for the same reason the layout does: it would be
        // one short, and there is no position to invent for the new atom. The
        // site table goes with them, and so does the topology — `add_atom` is
        // the only method that can change the atom count, so it is the only
        // place this invariant has to be maintained, and missing one table
        // here is how a charge silently ends up on the wrong atom or an atom
        // vanishes from a PDB write by belonging to no residue.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.set_coords3(vec![Point3::ORIGIN]).unwrap();
        mol.set_coords(vec![Point2::ORIGIN]).unwrap();
        mol.set_sites(vec![AtomSite {
            b_factor: Some(12.0),
            ..AtomSite::default()
        }])
        .unwrap();

        mol.set_topology(
            vec![Chain {
                id: "A".to_string(),
                label_id: None,
                residues: 0..1,
            }],
            vec![Residue {
                name: "LIG".to_string(),
                sequence: 1,
                insertion_code: None,
                label_seq: None,
                chain_ix: 0,
                is_hetero: true,
                atoms: 0..1,
            }],
        )
        .unwrap();

        mol.add_atom(Atom::new(Element::oxygen()));

        assert!(!mol.has_coords3());
        assert!(!mol.has_coords());
        assert!(!mol.has_sites());
        assert!(!mol.has_topology());
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
