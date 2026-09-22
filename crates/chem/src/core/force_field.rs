//! Per-atom force-field facts, and the angle/dihedral/improper/exclusion
//! lists that make a simulation reproducible (#315).
//!
//! PSF, CHARMM/Gromacs `TOP` and Amber `PRMTOP` all state these; `Molecule`
//! has bonds and nothing else on that list. Held as a side table on
//! [`crate::core::molecule::Molecule`], the way
//! [`crate::core::site::AtomSite`] is: `Atom` derives `Eq`, and masses and
//! charges are floats.

/// Per-atom force-field facts: a type name, a mass, a partial charge.
///
/// Every field independently optional, the same shape and the same reason
/// as [`crate::core::site::AtomSite`]: PSF, TOP and PRMTOP each state a
/// different subset.
///
/// `partial_charge` genuinely duplicates a field `AtomSite` already has —
/// deliberately. A Mol2 charge and an Amber charge are different
/// provenances; [`crate::io::format::held`] reports
/// [`crate::io::format::Carries::PARTIAL_CHARGE`] if *either* source has
/// one, so keeping them distinct loses nothing.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ForceFieldAtom {
    pub atom_type: Option<String>,
    pub mass: Option<f64>,
    pub partial_charge: Option<f64>,
}

impl ForceFieldAtom {
    /// An entry with nothing recorded — the same as [`Default`], named for
    /// the places where `ForceFieldAtom::empty()` reads better than
    /// `default()`.
    pub fn empty() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.atom_type.is_none() && self.mass.is_none() && self.partial_charge.is_none()
    }
}

/// A force field's topology: which atoms participate in which angle,
/// dihedral, improper and exclusion terms, plus their per-atom facts.
///
/// **Not the force-field parameters.** A term here names the atoms it
/// relates, nothing more — no force constant, no equilibrium geometry. Those
/// live in a separate parameter file (CHARMM's `.prm`, for one) that is out
/// of scope for this milestone. "Topology" in the name is the connectivity,
/// the same sense a PSF (protein *structure* file) uses it in.
///
/// Built plain and unvalidated — validating a term's atom indices needs the
/// owning [`crate::core::molecule::Molecule`]'s own atom count, which this
/// type does not have. [`crate::core::molecule::Molecule::set_force_field`]
/// is where that happens, the same division [`crate::core::residue::Chain`]/
/// [`crate::core::residue::Residue`] have with
/// [`crate::core::molecule::Molecule::set_topology`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ForceFieldTopology {
    /// One entry per atom, parallel to `Molecule::atoms`, when present.
    pub atoms: Option<Vec<ForceFieldAtom>>,
    /// Three-atom terms: which atoms are angled together.
    pub angles: Vec<[usize; 3]>,
    /// Four-atom terms: proper torsions.
    pub dihedrals: Vec<[usize; 4]>,
    /// Four-atom terms: out-of-plane (improper) torsions.
    pub impropers: Vec<[usize; 4]>,
    /// Atom pairs excluded from nonbonded interactions.
    pub exclusions: Vec<[usize; 2]>,
    /// Hydrogen-bond donor pairs: `[donor, hydrogen]` (#321).
    pub donors: Vec<[usize; 2]>,
    /// Hydrogen-bond acceptor pairs: `[acceptor, antecedent]` (#321).
    pub acceptors: Vec<[usize; 2]>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_new_atom_records_nothing() {
        let atom = ForceFieldAtom::default();
        assert!(atom.is_empty());
        assert_eq!(atom, ForceFieldAtom::empty());
        assert_eq!(atom.atom_type, None);
        assert_eq!(atom.mass, None);
    }

    #[test]
    fn test_fields_are_independently_optional() {
        // The point of three `Option`s rather than one around the struct: a
        // PSF supplies a type and a mass, an Amber PRMTOP supplies a charge
        // scaled by its own 18.2223 factor, and neither should force the
        // other to invent a value.
        let from_psf = ForceFieldAtom {
            atom_type: Some("CT".to_string()),
            mass: Some(12.011),
            ..ForceFieldAtom::default()
        };
        assert!(!from_psf.is_empty());
        assert_eq!(from_psf.partial_charge, None);
    }

    #[test]
    fn test_a_new_topology_states_no_terms() {
        let topology = ForceFieldTopology::default();
        assert!(topology.atoms.is_none());
        assert!(topology.angles.is_empty());
        assert!(topology.dihedrals.is_empty());
        assert!(topology.impropers.is_empty());
        assert!(topology.exclusions.is_empty());
        assert!(topology.donors.is_empty());
        assert!(topology.acceptors.is_empty());
    }
}
