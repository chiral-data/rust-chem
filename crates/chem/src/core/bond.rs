use std::fmt;

/// The order of a chemical bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
}

impl BondOrder {
    pub const fn value(&self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Quadruple => 4.0,
            BondOrder::Aromatic => 1.5,
        }
    }

    pub const fn is_rotatable(&self) -> bool {
        matches!(self, BondOrder::Single)
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BondOrder::Single => write!(f, "-"),
            BondOrder::Double => write!(f, "="),
            BondOrder::Triple => write!(f, "≡"),
            BondOrder::Quadruple => write!(f, "≣"),
            BondOrder::Aromatic => write!(f, ":"),
        }
    }
}

/// Stereochemistry of a bond (E/Z isomerism).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    None,
    E,
    Z,
    Unspecified,
}

/// The type of a bond for more specific classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondType {
    Covalent,
    Ionic,
    Hydrogen,
    Dative,
    Any,
}

/// Represents a chemical bond between two atoms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    atom1: usize,
    atom2: usize,
    order: BondOrder,
    is_aromatic: bool,
    stereo: BondStereo,
    bond_type: BondType,
    in_ring: bool,
}

impl Bond {
    pub const fn new(atom1: usize, atom2: usize, order: BondOrder) -> Self {
        Bond {
            atom1,
            atom2,
            order,
            is_aromatic: false,
            stereo: BondStereo::None,
            bond_type: BondType::Covalent,
            in_ring: false,
        }
    }

    pub const fn atom1(&self) -> usize {
        self.atom1
    }

    pub const fn atom2(&self) -> usize {
        self.atom2
    }

    pub const fn atoms(&self) -> (usize, usize) {
        (self.atom1, self.atom2)
    }

    pub const fn order(&self) -> BondOrder {
        self.order
    }

    pub fn set_order(&mut self, order: BondOrder) {
        self.order = order;
    }

    pub const fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    pub const fn with_aromatic(mut self, aromatic: bool) -> Self {
        self.is_aromatic = aromatic;
        self
    }

    pub fn set_aromatic(&mut self, aromatic: bool) {
        self.is_aromatic = aromatic;
    }

    pub const fn stereo(&self) -> BondStereo {
        self.stereo
    }

    pub const fn with_stereo(mut self, stereo: BondStereo) -> Self {
        self.stereo = stereo;
        self
    }

    pub const fn bond_type(&self) -> BondType {
        self.bond_type
    }

    pub const fn with_type(mut self, bond_type: BondType) -> Self {
        self.bond_type = bond_type;
        self
    }

    pub const fn in_ring(&self) -> bool {
        self.in_ring
    }

    pub fn set_in_ring(&mut self, in_ring: bool) {
        self.in_ring = in_ring;
    }

    pub const fn connects_to(&self, atom_idx: usize) -> bool {
        self.atom1 == atom_idx || self.atom2 == atom_idx
    }

    pub const fn other_atom(&self, atom_idx: usize) -> Option<usize> {
        if self.atom1 == atom_idx {
            Some(self.atom2)
        } else if self.atom2 == atom_idx {
            Some(self.atom1)
        } else {
            None
        }
    }

    pub const fn is_rotatable(&self) -> bool {
        self.order.is_rotatable() && !self.in_ring
    }

    pub fn compute_hash(&self) -> u64 {
        let mut hash = (self.order.value() * 10.0) as u64;
        hash = hash.wrapping_mul(31).wrapping_add(self.is_aromatic as u64);
        hash
    }
}

impl fmt::Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}{}", self.atom1, self.atom2, self.order)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_order_value() {
        assert_eq!(BondOrder::Single.value(), 1.0);
        assert_eq!(BondOrder::Double.value(), 2.0);
        assert_eq!(BondOrder::Triple.value(), 3.0);
        assert_eq!(BondOrder::Aromatic.value(), 1.5);
    }

    #[test]
    fn test_bond_creation() {
        let bond = Bond::new(0, 1, BondOrder::Single);
        assert_eq!(bond.atom1(), 0);
        assert_eq!(bond.atom2(), 1);
        assert_eq!(bond.order(), BondOrder::Single);
        assert!(!bond.is_aromatic());
    }

    #[test]
    fn test_bond_connectivity() {
        let bond = Bond::new(5, 10, BondOrder::Double);
        assert!(bond.connects_to(5));
        assert!(bond.connects_to(10));
        assert!(!bond.connects_to(3));
        assert_eq!(bond.other_atom(5), Some(10));
        assert_eq!(bond.other_atom(10), Some(5));
        assert_eq!(bond.other_atom(3), None);
    }

    #[test]
    fn test_bond_rotatable() {
        let single = Bond::new(0, 1, BondOrder::Single);
        assert!(single.is_rotatable());
        let double = Bond::new(0, 1, BondOrder::Double);
        assert!(!double.is_rotatable());
        let mut ring_bond = Bond::new(0, 1, BondOrder::Single);
        ring_bond.set_in_ring(true);
        assert!(!ring_bond.is_rotatable());
    }

    #[test]
    fn test_bond_builder() {
        let bond = Bond::new(0, 1, BondOrder::Single)
            .with_aromatic(true)
            .with_stereo(BondStereo::E);
        assert!(bond.is_aromatic());
        assert_eq!(bond.stereo(), BondStereo::E);
    }
}
