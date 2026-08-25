use std::fmt;

use crate::core::elements::{ELEMENT_NAMES, ELEMENT_SYMBOLS};

/// Represents a chemical element.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    pub atomic_number: u8,
}

impl Element {
    pub const fn new(atomic_number: u8) -> Option<Self> {
        if atomic_number == 0 || atomic_number > 118 {
            None
        } else {
            Some(Element { atomic_number })
        }
    }

    pub const fn symbol(&self) -> &'static str {
        ELEMENT_SYMBOLS[self.atomic_number as usize]
    }

    pub const fn name(&self) -> &'static str {
        ELEMENT_NAMES[self.atomic_number as usize]
    }

    pub const fn typical_valence(&self) -> u8 {
        match self.atomic_number {
            1 => 1,  // H
            5 => 3,  // B
            6 => 4,  // C
            7 => 3,  // N
            8 => 2,  // O
            9 => 1,  // F
            15 => 3, // P
            16 => 2, // S
            17 => 1, // Cl
            35 => 1, // Br
            53 => 1, // I
            _ => 0,
        }
    }

    pub const fn hydrogen() -> Self {
        Element { atomic_number: 1 }
    }
    pub const fn carbon() -> Self {
        Element { atomic_number: 6 }
    }
    pub const fn nitrogen() -> Self {
        Element { atomic_number: 7 }
    }
    pub const fn oxygen() -> Self {
        Element { atomic_number: 8 }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

/// Hybridization state of an atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Hybridization {
    SP,
    SP2,
    SP3,
    SP3D,
    SP3D2,
    Unknown,
}

/// Chirality/stereochemistry of an atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Chirality {
    None,
    Clockwise,
    CounterClockwise,
    Unspecified,
}

/// Represents an atom in a molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    element: Element,
    formal_charge: i8,
    isotope: Option<u16>,
    is_aromatic: bool,
    hybridization: Hybridization,
    chirality: Chirality,
    implicit_hydrogens: u8,
    explicit_hydrogens: u8,
}

impl Atom {
    pub const fn new(element: Element) -> Self {
        Atom {
            element,
            formal_charge: 0,
            isotope: None,
            is_aromatic: false,
            hybridization: Hybridization::Unknown,
            chirality: Chirality::None,
            implicit_hydrogens: 0,
            explicit_hydrogens: 0,
        }
    }

    pub const fn element(&self) -> Element {
        self.element
    }

    pub const fn atomic_number(&self) -> u8 {
        self.element.atomic_number
    }

    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    pub const fn with_charge(mut self, charge: i8) -> Self {
        self.formal_charge = charge;
        self
    }

    pub const fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    pub const fn with_isotope(mut self, mass: u16) -> Self {
        self.isotope = Some(mass);
        self
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

    pub const fn hybridization(&self) -> Hybridization {
        self.hybridization
    }

    pub const fn with_hybridization(mut self, hyb: Hybridization) -> Self {
        self.hybridization = hyb;
        self
    }

    pub const fn chirality(&self) -> Chirality {
        self.chirality
    }

    pub const fn with_chirality(mut self, chir: Chirality) -> Self {
        self.chirality = chir;
        self
    }

    pub const fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
    }

    pub fn set_implicit_hydrogens(&mut self, count: u8) {
        self.implicit_hydrogens = count;
    }

    pub const fn explicit_hydrogens(&self) -> u8 {
        self.explicit_hydrogens
    }

    pub fn set_explicit_hydrogens(&mut self, count: u8) {
        self.explicit_hydrogens = count;
    }

    pub const fn total_hydrogens(&self) -> u8 {
        self.implicit_hydrogens + self.explicit_hydrogens
    }

    pub fn compute_hash(&self) -> u64 {
        let mut hash = self.atomic_number() as u64;
        hash = hash
            .wrapping_mul(31)
            .wrapping_add((self.formal_charge + 5) as u64);
        hash = hash.wrapping_mul(31).wrapping_add(self.is_aromatic as u64);
        hash
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.element.symbol())?;
        if self.formal_charge != 0 {
            write!(f, "{:+}", self.formal_charge)?;
        }
        if let Some(iso) = self.isotope {
            write!(f, "[{}]", iso)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_element_creation() {
        let carbon = Element::new(6).unwrap();
        assert_eq!(carbon.atomic_number, 6);
        assert_eq!(carbon.symbol(), "C");
        assert_eq!(carbon.name(), "Carbon");
        assert!(Element::new(0).is_none());
        assert!(Element::new(119).is_none());
    }

    #[test]
    fn test_element_valence() {
        assert_eq!(Element::carbon().typical_valence(), 4);
        assert_eq!(Element::nitrogen().typical_valence(), 3);
        assert_eq!(Element::oxygen().typical_valence(), 2);
    }

    #[test]
    fn test_atom_creation() {
        let atom = Atom::new(Element::carbon());
        assert_eq!(atom.atomic_number(), 6);
        assert_eq!(atom.formal_charge(), 0);
        assert!(!atom.is_aromatic());
    }

    #[test]
    fn test_atom_builder() {
        let atom = Atom::new(Element::nitrogen())
            .with_charge(1)
            .with_aromatic(true)
            .with_isotope(15);
        assert_eq!(atom.formal_charge(), 1);
        assert!(atom.is_aromatic());
        assert_eq!(atom.isotope(), Some(15));
    }

    #[test]
    fn test_atom_hash() {
        let c1 = Atom::new(Element::carbon());
        let c2 = Atom::new(Element::carbon()).with_charge(1);
        let n = Atom::new(Element::nitrogen());
        assert_eq!(
            c1.compute_hash(),
            Atom::new(Element::carbon()).compute_hash()
        );
        assert_ne!(c1.compute_hash(), c2.compute_hash());
        assert_ne!(c1.compute_hash(), n.compute_hash());
    }
}
