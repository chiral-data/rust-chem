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

    /// The valence this element supports while carrying `charge`.
    ///
    /// Deliberately not `typical_valence() + charge`. Which way a charge moves
    /// the valence depends on where the element sits in its period, and two
    /// call sites disagreeing about that was #240 -- so the rule lives here,
    /// once, rather than being re-derived wherever a hydrogen count is filled
    /// in. Verified against RDKit 2025.3.3 for every element
    /// [`Self::typical_valence`] covers.
    ///
    /// Zero for an element with no typical valence, and for a charge that
    /// would drive the valence below it: a chloride has no bonds left to give.
    pub const fn valence_for_charge(&self, charge: i8) -> u8 {
        let base = self.typical_valence() as i16;
        if base == 0 {
            return 0;
        }
        let charge = charge as i16;
        let adjusted = match self.atomic_number {
            // A proton has no electrons and a hydride no bonds, so both carry
            // none; carbon loses a bond to an ion of either sign, `[CH3-]` and
            // `[CH3+]` being three apiece.
            1 | 6 => base - charge.abs(),
            // Boron is electron-deficient, so an extra electron buys it another
            // bond rather than costing one: `[BH4-]` is four.
            5 => base - charge,
            // N, O, P, S and the halogens, the majority and the ones #240
            // reported: `[NH4+]` is four, `[O-]` is one.
            _ => base + charge,
        };
        if adjusted < 0 { 0 } else { adjusted as u8 }
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
    /// Hydrogens attached to this atom but not present as their own graph
    /// nodes -- and whether the source said anything about them at all.
    ///
    /// `None` means the input was silent, so a reader may fill it in by
    /// valence ([`crate::core::molecule::Molecule::calculate_implicit_hydrogens`]).
    /// `Some(0)` means the input said *none*, and nothing may add any. A bare
    /// SMILES `C` is the first; `[C]` and `[CH0]` are the second. Before #244
    /// both were a plain `0` and `[C]` parsed as methane.
    hydrogens: Option<u8>,
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
            hydrogens: None,
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

    /// The hydrogen count as the source stated it, or `None` if it did not.
    ///
    /// Use this to tell "said none" from "said nothing"; use
    /// [`Self::total_hydrogens`] when only the number matters, which is most
    /// callers.
    pub const fn hydrogens(&self) -> Option<u8> {
        self.hydrogens
    }

    pub fn set_hydrogens(&mut self, count: u8) {
        self.hydrogens = Some(count);
    }

    pub const fn with_hydrogens(mut self, count: u8) -> Self {
        self.hydrogens = Some(count);
        self
    }

    /// Forget the count, so a reader may fill it in again.
    pub fn clear_hydrogens(&mut self) {
        self.hydrogens = None;
    }

    /// How many hydrogens this atom carries, treating "unstated" as none.
    ///
    /// The number every consumer that just wants a count should ask for --
    /// formula and mass, the Morgan invariants, AutoDock typing,
    /// kekulisation, and the SMILES writer's bracket decision.
    pub const fn total_hydrogens(&self) -> u8 {
        match self.hydrogens {
            Some(count) => count,
            None => 0,
        }
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
    fn test_hydrogens_distinguishes_saying_none_from_saying_nothing() {
        // The distinction #244 exists for. Both answer 0 to
        // `total_hydrogens`, and only one of them invites a reader to fill it.
        let unsaid = Atom::new(Element::carbon());
        assert_eq!(unsaid.hydrogens(), None);
        assert_eq!(unsaid.total_hydrogens(), 0);

        let mut stated_none = Atom::new(Element::carbon());
        stated_none.set_hydrogens(0);
        assert_eq!(stated_none.hydrogens(), Some(0));
        assert_eq!(stated_none.total_hydrogens(), 0);

        assert_ne!(unsaid, stated_none, "the two states are not the same atom");
    }

    #[test]
    fn test_hydrogens_can_be_set_cleared_and_built() {
        let mut atom = Atom::new(Element::nitrogen()).with_hydrogens(3);
        assert_eq!(atom.hydrogens(), Some(3));
        assert_eq!(atom.total_hydrogens(), 3);

        atom.set_hydrogens(1);
        assert_eq!(atom.hydrogens(), Some(1));

        atom.clear_hydrogens();
        assert_eq!(atom.hydrogens(), None, "back to fillable");
        assert_eq!(atom.total_hydrogens(), 0);
    }

    #[test]
    fn test_valence_for_charge_matches_the_oracle() {
        // The table as RDKit 2025.3.3 reports it: total implicit hydrogens on
        // a bare atom at each charge, for every element `typical_valence`
        // covers. Written out rather than computed, because a test that
        // re-derives the implementation proves only that it equals itself.
        //
        //                     Z    -1   0  +1
        let table: &[(u8, [u8; 3])] = &[
            (1, [0, 1, 0]),  // H  -- a proton has none, a hydride none
            (5, [4, 3, 2]),  // B  -- electron-deficient, gains with -1
            (6, [3, 4, 3]),  // C  -- loses either way
            (7, [2, 3, 4]),  // N
            (8, [1, 2, 3]),  // O
            (9, [0, 1, 2]),  // F
            (15, [2, 3, 4]), // P
            (16, [1, 2, 3]), // S
            (17, [0, 1, 2]), // Cl
            (35, [0, 1, 2]), // Br
            (53, [0, 1, 2]), // I
        ];

        for (z, expected) in table {
            let element = Element::new(*z).expect("real element");
            for (i, charge) in [-1i8, 0, 1].iter().enumerate() {
                assert_eq!(
                    element.valence_for_charge(*charge),
                    expected[i],
                    "{} at charge {charge}",
                    element.symbol()
                );
            }
        }
    }

    #[test]
    fn test_valence_for_charge_is_not_a_single_sign_flip() {
        // The two rows a naive `typical_valence() + charge` gets wrong, and
        // the reason #240 was not a one-character fix. Boron is the one the
        // *original* code happened to get right.
        assert_eq!(Element::new(5).unwrap().valence_for_charge(-1), 4, "[BH4-]");
        assert_eq!(Element::carbon().valence_for_charge(1), 3, "[CH3+]");
        assert_eq!(Element::carbon().valence_for_charge(-1), 3, "[CH3-]");
    }

    #[test]
    fn test_valence_for_charge_floors_at_zero() {
        // A chloride has no bonds left to give, and a doubly charged one is
        // not a reason to underflow.
        assert_eq!(Element::new(17).unwrap().valence_for_charge(-1), 0);
        assert_eq!(Element::new(17).unwrap().valence_for_charge(-3), 0);
        assert_eq!(Element::oxygen().valence_for_charge(-5), 0);
    }

    #[test]
    fn test_valence_for_charge_is_zero_where_there_is_no_valence() {
        // Sodium, magnesium and every metal: `typical_valence` declines to
        // guess, and a charge must not talk it into one.
        for z in [11u8, 12, 26, 79] {
            let element = Element::new(z).expect("real element");
            assert_eq!(element.typical_valence(), 0, "{}", element.symbol());
            for charge in [-2i8, -1, 0, 1, 2] {
                assert_eq!(
                    element.valence_for_charge(charge),
                    0,
                    "{}",
                    element.symbol()
                );
            }
        }
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
