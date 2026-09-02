//! Writes a SMILES string for a molecule — the reverse direction of
//! [`crate::io::smiles::parse_smiles`].
//!
//! The traversal algorithm (two passes: find ring-closure digits, then
//! build the string) is from *Tutorials in Chemoinformatics*, pages 422-425,
//! by Alexandre Varnek.
//!
//! # Examples
//!
//! ```rust
//! use chem::io::smiles::parse_smiles;
//! use chem::io::smiles_writer::write_smiles_for_molecule;
//!
//! let mol = parse_smiles("Oc1ccccc1").unwrap();
//! assert_eq!(write_smiles_for_molecule(&mol), "Oc1ccccc1");
//! ```
//!
//! or write from raw atom/bond data directly, without a [`crate::core::molecule::Molecule`]:
//!
//! ```rust
//! use chem::io::smiles_writer::write_smiles;
//!
//! let atom_symbols: Vec<String> = vec!["c", "c", "c", "c", "c", "c", "O"].into_iter().map(|s| s.to_string()).collect();
//! let atom_neighbours: Vec<Vec<usize>> = vec![vec![1, 5], vec![2, 0], vec![1, 3], vec![2, 4], vec![3, 5], vec![0, 4, 6], vec![5]];
//! let atom_rankings: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6];
//! let bond_symbols: std::collections::HashMap<String, String> = [("0,1", ""), ("0,5", ""), ("1,0", ""), ("1,2", ""), ("2,1", ""), ("2,3", ""), ("3,2", ""), ("3,4", ""), ("4,3", ""), ("4,5", ""), ("5,4", ""), ("5,0", ""), ("5,6", "="), ("6,5", "=")].iter().map(|s| (s.0.to_string(), s.1.to_string())).collect();
//! assert_eq!(write_smiles(atom_symbols, atom_neighbours, atom_rankings, bond_symbols), "c1ccccc1=O".to_string());
//! ```

use crate::core::atom::Chirality;
use crate::core::bond::{BondOrder, BondStereo};
use crate::core::molecule::Molecule;

/// Manager of ring digit
struct DigitHeap {
    digits: Vec<usize>,
}

impl DigitHeap {
    fn init() -> Self {
        Self { digits: vec![] }
    }

    fn find(&mut self) -> usize {
        let mut digit: usize = 1;
        while digit < 100 {
            if !self.digits.contains(&digit) {
                break;
            }
            digit += 1;
        }

        self.digits.push(digit);
        digit
    }

    fn remove(&mut self, digit: &usize) {
        let index = self.digits.iter().position(|x| *x == *digit).unwrap();
        self.digits.remove(index);
    }
}

/// Intermediate Data
type OpenAtomDigitTuple = (usize, usize);
type OpeningClosures = std::collections::HashMap<usize, Vec<usize>>;
type ClosingClosures = std::collections::HashMap<usize, Vec<OpenAtomDigitTuple>>;

struct DataPool {
    pub ancestors: Vec<usize>,
    pub visited: Vec<usize>,
    pub opening_closures: OpeningClosures,
    pub closing_closures: ClosingClosures,
    pub dh: DigitHeap,
}

impl DataPool {
    fn init() -> Self {
        Self {
            ancestors: vec![],
            visited: vec![],
            opening_closures: OpeningClosures::new(),
            closing_closures: ClosingClosures::new(),
            dh: DigitHeap::init(),
        }
    }
}

/// Implement this trait to write SMILES for your own molecule type via
/// [`write_smiles_for_mol`]. Already implemented for
/// [`crate::core::molecule::Molecule`] below.
pub trait TraitMoleculeForSMILES {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize>;
    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String;
    fn get_atom_symbol(&self, atom: &usize) -> String;

    /// The atom's token, given the order this traversal will write its
    /// neighbours in.
    ///
    /// Only stereocentres care. [`Chirality`] is stored against sorted
    /// neighbour order — see `normalise_chirality` in the parser — so a writer
    /// has to convert into whatever order it is about to emit, and that order
    /// is a property of the traversal rather than of the molecule.
    ///
    /// Defaults to ignoring it, so an implementation with no stereochemistry
    /// need not care.
    fn get_atom_symbol_written(&self, atom: &usize, written: &[usize]) -> String {
        let _ = written;
        self.get_atom_symbol(atom)
    }
    fn get_atom_ranking(&self, atom: &usize, rankings: &[usize]) -> usize;
    fn count_of_atoms(&self) -> usize;
}

struct MoleculeForSmiles {
    atom_symbols: Vec<String>,
    atom_neighbours: Vec<Vec<usize>>,
    bond_symbols: std::collections::HashMap<String, String>,
}

impl TraitMoleculeForSMILES for MoleculeForSmiles {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
        self.atom_neighbours[*atom].clone()
    }

    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
        match self.bond_symbols.get(&format!("{},{}", *atom_1, *atom_2)) {
            Some(symbol) => symbol.clone(),
            None => panic!("No such bond"),
        }
    }

    fn get_atom_symbol(&self, atom: &usize) -> String {
        self.atom_symbols[*atom].clone()
    }

    fn get_atom_ranking(&self, atom: &usize, rankings: &[usize]) -> usize {
        rankings[*atom]
    }

    fn count_of_atoms(&self) -> usize {
        self.atom_symbols.len()
    }
}

// Maps a bond's order to the symbol written between its two atoms. Single
// and Aromatic both write as "" (implicit) rather than an explicit "-" or
// ":" — parse_smiles fills in exactly this default itself (Aromatic if the
// atom being connected to is aromatic, Single otherwise), so leaving it
// implicit is what makes the output re-parse back to the same bond order.
/// The implicit hydrogen's place in a written neighbour order.
///
/// Matches the parser's `IMPLICIT_H`: `usize::MAX`, so it sorts last, which is
/// where every index-based convention puts it.
const IMPLICIT_H_SLOT: usize = usize::MAX;

/// The atom's token, with any chirality marker converted out of storage order
/// into the order this traversal writes.
///
/// [`Chirality`] means "looking from the lowest-indexed neighbour, the rest in
/// increasing index order". A SMILES marker means "looking from the first
/// neighbour *as written*". The two agree when the written order is an even
/// permutation of the sorted one and disagree when it is odd, which is exactly
/// what the parser corrected for on the way in.
///
/// Without this, writing was self-consistent and meaningless: the stored value
/// went out unchanged and came back unchanged, so a round trip through this
/// crate looked perfect while `[C@@H](N)(C)C(=O)O` and `N[C@@H](C)C(=O)O`
/// claimed different configurations for the same molecule.
fn atom_symbol_in_order(mol: &Molecule, atom_idx: usize, written: &[usize]) -> String {
    let chirality = mol.atom(atom_idx).chirality();
    if chirality != Chirality::Clockwise && chirality != Chirality::CounterClockwise {
        return atom_symbol(mol, atom_idx);
    }

    // Only the slots that actually name something: the hydrogen placeholder is
    // present in the written order only when the atom has one.
    let mut order: Vec<usize> = written.to_vec();
    if mol.atom(atom_idx).total_hydrogens() == 0 {
        order.retain(|slot| *slot != IMPLICIT_H_SLOT);
    }
    let mut sorted = order.clone();
    sorted.sort_unstable();

    let effective = if crate::io::smiles::permutation_is_odd(&order, &sorted) {
        match chirality {
            Chirality::Clockwise => Chirality::CounterClockwise,
            _ => Chirality::Clockwise,
        }
    } else {
        chirality
    };

    atom_symbol_with(mol, atom_idx, effective)
}

fn bond_order_symbol(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single | BondOrder::Aromatic => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Quadruple => "$",
    }
}

// Writes an atom's SMILES token. Plain element symbol (lowercase if
// aromatic) for the common case; bracket notation only when charge or
// isotope need representing, since those can't be expressed in the
// organic-subset shorthand.
fn atom_symbol(mol: &Molecule, atom_idx: usize) -> String {
    atom_symbol_with(mol, atom_idx, mol.atom(atom_idx).chirality())
}

/// [`atom_symbol`], with the chirality marker stated explicitly.
///
/// Split out so a caller that has converted the configuration into its own
/// neighbour order can say so, without cloning the molecule to carry one field.
fn atom_symbol_with(mol: &Molecule, atom_idx: usize, chirality: Chirality) -> String {
    let atom = mol.atom(atom_idx);
    let symbol = atom.element().symbol();
    let symbol = if atom.is_aromatic() {
        symbol.to_lowercase()
    } else {
        symbol.to_string()
    };

    let needs_brackets = atom.formal_charge() != 0
        || atom.isotope().is_some()
        || atom.explicit_hydrogens() > 0
        || chirality != Chirality::None;
    if !needs_brackets {
        return symbol;
    }

    let mut bracket = String::from("[");
    if let Some(isotope) = atom.isotope() {
        bracket += &isotope.to_string();
    }
    bracket += &symbol;
    // Between the symbol and the hydrogen count, which is the order the spec
    // fixes and the order `parse_bracket_atom` reads. Kept adjacent in both
    // files for that reason.
    //
    // This re-emits the marker as it was read. It does *not* recompute parity
    // for a traversal that visits this atom's neighbours in a different order
    // than the input wrote them, so a round trip preserves the marker rather
    // than provably the configuration. Fixing that needs a canonical
    // neighbour order, which belongs to the canonical-writer story.
    bracket += match chirality {
        Chirality::CounterClockwise => "@",
        Chirality::Clockwise => "@@",
        Chirality::None | Chirality::Unspecified => "",
    };
    match atom.explicit_hydrogens() {
        0 => {}
        1 => bracket += "H",
        h => bracket += &format!("H{}", h),
    }
    match atom.formal_charge() {
        0 => {}
        1 => bracket += "+",
        -1 => bracket += "-",
        c if c > 0 => bracket += &format!("+{}", c),
        c => bracket += &format!("-{}", -c),
    }
    bracket += "]";
    bracket
}

/// The reference substituent whose bond carries the marker for `anchor`.
///
/// One per end of the double bond, chosen deterministically so both ends agree
/// without passing state between them: the lowest-indexed neighbour that is not
/// the double-bond partner.
fn reference_substituent(mol: &Molecule, anchor: usize, partner: usize) -> Option<usize> {
    mol.neighbors(anchor)
        .iter()
        .map(|n| n.atom_idx)
        .filter(|&n| n != partner)
        .min()
}

/// `/` or `\` for a single bond written `first` → `second`, if it is the bond
/// carrying a neighbouring double bond's configuration.
///
/// The inverse of `resolve_directional_bonds` in the parser, and it has to
/// agree with it exactly or a round trip flips cis and trans. Both use the same
/// rule: a substituent's *side* of the double-bond axis is the marker's lean
/// when the double-bond atom is written first, and the negated lean when it is
/// written second. Same sides is Z, opposite sides is E.
///
/// Here the sides are chosen and the leans derived. The first end's substituent
/// is arbitrarily placed below the axis; the second end's follows from the
/// configuration. Then each lean is whatever writes that side, given which way
/// this particular traversal happens to be emitting the bond.
///
/// Returns `None` when the configuration cannot be written on this bond, which
/// is not an error: a molecule whose reference substituent is reached by a ring
/// closure rather than by this edge simply loses the marker. That degrades the
/// output, unlike a dropped component, which would corrupt it.
fn directional_marker(mol: &Molecule, first: usize, second: usize) -> Option<char> {
    for bond in mol.bonds() {
        if bond.order() != BondOrder::Double {
            continue;
        }
        let stereo = bond.stereo();
        if stereo != BondStereo::E && stereo != BondStereo::Z {
            continue;
        }
        let (left, right) = (bond.atom1(), bond.atom2());

        // Which end of the double bond does this single bond hang off, and
        // which atom is the substituent?
        let (anchor, substituent, partner, is_second_end) = if first == left || second == left {
            let substituent = if first == left { second } else { first };
            (left, substituent, right, false)
        } else if first == right || second == right {
            let substituent = if first == right { second } else { first };
            (right, substituent, left, true)
        } else {
            continue;
        };

        // The substituent must not be the other end of the double bond itself.
        if substituent == partner {
            continue;
        }
        if reference_substituent(mol, anchor, partner) != Some(substituent) {
            continue;
        }

        // The first end's substituent is placed below the axis by fiat — the
        // choice is arbitrary, only the relationship matters. The second end's
        // is then above for E (opposite sides) and below for Z (same side).
        let side: i8 = if is_second_end && stereo == BondStereo::E {
            1
        } else {
            -1
        };

        // Convert the side back into a lean, using this traversal's direction.
        let lean = if anchor == first { side } else { -side };
        return Some(if lean > 0 { '/' } else { '\\' });
    }
    None
}

impl TraitMoleculeForSMILES for Molecule {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
        self.neighbors(*atom).iter().map(|n| n.atom_idx).collect()
    }

    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
        match self.graph().get_bond(*atom_1, *atom_2) {
            Some(bond_idx) => {
                let bond = self.bond(bond_idx);
                // A single bond next to a stereo double bond carries the
                // configuration, so it is asked about first: `/` and `\` take
                // the slot a `-` would have, and `bond_order_symbol` writes
                // nothing for a single bond anyway.
                if bond.order() == BondOrder::Single
                    && let Some(marker) = directional_marker(self, *atom_1, *atom_2)
                {
                    return marker.to_string();
                }
                bond_order_symbol(bond.order()).to_string()
            }
            None => panic!("No such bond"),
        }
    }

    fn get_atom_symbol(&self, atom: &usize) -> String {
        atom_symbol(self, *atom)
    }

    fn get_atom_symbol_written(&self, atom: &usize, written: &[usize]) -> String {
        atom_symbol_in_order(self, *atom, written)
    }

    fn get_atom_ranking(&self, atom: &usize, rankings: &[usize]) -> usize {
        rankings[*atom]
    }

    fn count_of_atoms(&self) -> usize {
        self.num_atoms()
    }
}

#[inline(always)]
fn get_neighbours_excluding_parent<T: TraitMoleculeForSMILES>(
    mol: &T,
    atom_current: usize,
    atom_parent_opt: Option<usize>,
) -> Vec<usize> {
    mol.get_neighbours_of_atom(&atom_current)
        .into_iter()
        .filter(|&idx| match atom_parent_opt {
            Some(atom_parent) => idx != atom_parent,
            None => true,
        })
        .collect()
}

/// A recursive function for detecting opening closures by the 1st traversing
fn get_closures_for_atom<T: TraitMoleculeForSMILES>(
    mol: &T,
    rankings: &[usize],
    atom_current: usize,
    atom_parent_opt: Option<usize>,
    dp: &mut DataPool,
) {
    dp.ancestors.push(atom_current);
    dp.visited.push(atom_current);

    let mut nbors: Vec<usize> = get_neighbours_excluding_parent(mol, atom_current, atom_parent_opt);
    nbors.sort_by_key(|idx| mol.get_atom_ranking(idx, rankings));

    for nb in nbors.iter() {
        if dp.ancestors.contains(nb) {
            dp.opening_closures
                .entry(*nb)
                .or_default()
                .push(atom_current);
        } else if !dp.visited.contains(nb) {
            get_closures_for_atom(mol, rankings, *nb, Some(atom_current), dp);
        }
    }

    let index = dp
        .ancestors
        .iter()
        .position(|x| *x == atom_current)
        .unwrap();
    dp.ancestors.remove(index);
}

/// A recursive function for building smiles by the 2nd traversing
fn build_smiles_for_atom<T: TraitMoleculeForSMILES>(
    mol: &T,
    rankings: &[usize],
    atom_current: usize,
    atom_parent_opt: Option<usize>,
    dp: &mut DataPool,
) -> String {
    dp.visited.push(atom_current);
    let mut seq: String = String::from("");

    if let Some(atom_parent) = atom_parent_opt {
        seq += &mol.get_bond_symbol(&atom_parent, &atom_current);
    }

    // The atom's own token is emitted last, once the order of everything after
    // it is known: a chirality marker is defined against the order this
    // traversal is about to write the neighbours in, and that is not settled
    // until the closures and branches below have been laid out.
    //
    // `written` accumulates that order. The atom hung off comes first, then the
    // bracket's hydrogen, then ring-closure digits, then branches.
    let mut written: Vec<usize> = Vec::new();
    if let Some(atom_parent) = atom_parent_opt {
        written.push(atom_parent);
    }
    written.push(IMPLICIT_H_SLOT);

    let mut after_symbol = String::new();

    if let Some(oadts) = dp.closing_closures.get_mut(&atom_current) {
        oadts.sort_by_key(|oadt| oadt.1); // close multiple rings, start from smaller digits
        for oadt in oadts.iter() {
            after_symbol += &mol.get_bond_symbol(&atom_current, &oadt.0);
            if oadt.1 > 9 {
                after_symbol += "%";
            }
            after_symbol += &oadt.1.to_string();
            written.push(oadt.0);
            dp.dh.remove(&oadt.1);
        }
    }

    if let Some(ocs) = dp.opening_closures.get(&atom_current) {
        for oc in ocs.iter() {
            let digit = dp.dh.find();
            if digit > 9 {
                after_symbol += "%";
            }
            after_symbol += &digit.to_string();
            written.push(*oc);
            let oadts = dp.closing_closures.entry(*oc).or_default();
            oadts.push((atom_current, digit));
        }
    }

    let mut nbors: Vec<usize> = get_neighbours_excluding_parent(mol, atom_current, atom_parent_opt);
    nbors.sort_by_key(|idx| mol.get_atom_ranking(idx, rankings));

    let mut branches: Vec<String> = vec![];
    for n in nbors.iter() {
        if !dp.visited.contains(n) {
            written.push(*n);
            branches.push(build_smiles_for_atom(
                mol,
                rankings,
                *n,
                Some(atom_current),
                dp,
            ));
        }
    }

    seq += &mol.get_atom_symbol_written(&atom_current, &written);
    seq += &after_symbol;

    if branches.len() > 1 {
        for branch in branches[..(branches.len() - 1)].iter() {
            seq += &format!("({})", branch);
        }
    }

    if let Some(last) = branches.last() {
        seq += last;
    }

    seq
}

/// SMILES writer with trait Molecule
///
/// Traverses once per connected component and joins them with `.`. A single
/// traversal from the lowest-ranked atom was enough while nothing could parse a
/// dot, but it meant a two-component molecule wrote back as its first component
/// alone — `[Na+].[Cl-]` became `[Na+]`, a different molecule returned as a
/// success. #191 made that reachable by teaching the parser `.`, so the two
/// changes ship together.
///
/// Components are found by exhaustion rather than by
/// [`crate::core::graph::MoleculeGraph::connected_components`]: this function
/// is generic over [`TraitMoleculeForSMILES`], which exposes neighbours but no
/// graph, so anything not reached by the previous traversal starts the next.
pub fn write_smiles_for_mol<T: TraitMoleculeForSMILES>(mol: &T, rankings: &[usize]) -> String {
    // Ranking order decides which component comes first and where each
    // traversal starts, exactly as it did for the single-component case.
    let mut atom_indexes: Vec<usize> = (0..mol.count_of_atoms()).collect();
    atom_indexes.sort_by_key(|idx| mol.get_atom_ranking(idx, rankings));

    let mut written: Vec<usize> = Vec::new();
    let mut components: Vec<String> = Vec::new();

    for &start in &atom_indexes {
        if written.contains(&start) {
            continue;
        }

        // A fresh pool per component. Ring digits are freed as closures close,
        // so sharing one would work — but a component that somehow leaked a
        // digit would silently renumber the next one's rings, and a salt is not
        // where anyone would look for that.
        let mut dp = DataPool::init();
        get_closures_for_atom(mol, rankings, start, None, &mut dp);

        // What the closure pass reached *is* the component.
        let reached = dp.visited.clone();
        dp.visited.clear();

        components.push(build_smiles_for_atom(mol, rankings, start, None, &mut dp));
        written.extend(reached);
    }

    components.join(".")
}

/// SMILES writer with raw data
pub fn write_smiles(
    atom_symbols: Vec<String>,
    atom_neighbours: Vec<Vec<usize>>,
    atom_rankings: Vec<usize>,
    bond_symbols: std::collections::HashMap<String, String>,
) -> String {
    let mol = MoleculeForSmiles {
        atom_symbols,
        atom_neighbours,
        bond_symbols,
    };
    write_smiles_for_mol(&mol, &atom_rankings)
}

/// Writes a SMILES string for a [`crate::core::molecule::Molecule`], using
/// atom index as the traversal ranking. Not canonical — the same molecule
/// built in a different atom order can write a different (but equally
/// valid) SMILES string.
pub fn write_smiles_for_molecule(mol: &Molecule) -> String {
    let rankings: Vec<usize> = (0..mol.num_atoms()).collect();
    write_smiles_for_mol(mol, &rankings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::smiles::parse_smiles;

    /// The configuration of the one double bond, after a full round trip.
    ///
    /// Asserting on the *string* would pin a spelling rather than a molecule:
    /// `C(/F)=C/F` legitimately comes back as `C(\F)=C\F`, which is the same
    /// isomer written from a different starting atom. What has to survive is
    /// the E or the Z.
    fn stereo_after_round_trip(smiles: &str) -> BondStereo {
        let written = write_smiles_for_molecule(&parse_smiles(smiles).expect(smiles));
        parse_smiles(&written)
            .unwrap_or_else(|e| panic!("{smiles} wrote {written}, which will not parse: {e}"))
            .bonds()
            .iter()
            .find(|b| b.order() == BondOrder::Double)
            .expect("a double bond")
            .stereo()
    }

    #[test]
    fn test_every_component_is_written() {
        // The silent-corruption case. A single traversal wrote the component
        // containing the lowest-ranked atom and dropped the rest, so a salt
        // came back as one of its ions — a different molecule, returned as a
        // success (#191).
        let mol = parse_smiles("[Na+].[Cl-]").expect("valid SMILES");
        let written = write_smiles_for_molecule(&mol);
        assert_eq!(written, "[Na+].[Cl-]");

        let back = parse_smiles(&written).expect("round trips");
        assert_eq!(back.num_atoms(), 2);
        assert_eq!(back.graph().connected_components().len(), 2);
    }

    #[test]
    fn test_components_are_written_in_ranking_order() {
        let mol = parse_smiles("CC.O").expect("valid SMILES");
        assert_eq!(write_smiles_for_molecule(&mol), "CC.O");
    }

    #[test]
    fn test_a_ring_closed_across_a_dot_writes_as_one_component() {
        // `C1.C1` is bonded, so there is nothing to separate and no dot to
        // write.
        let mol = parse_smiles("C1.C1").expect("valid SMILES");
        assert_eq!(write_smiles_for_molecule(&mol), "CC");
    }

    #[test]
    fn test_chirality_round_trips() {
        for smiles in ["N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O"] {
            let mol = parse_smiles(smiles).expect(smiles);
            let written = write_smiles_for_molecule(&mol);
            assert_eq!(written, smiles, "chirality marker lost or flipped");
        }
    }

    #[test]
    fn test_e_and_z_survive_a_round_trip() {
        assert_eq!(stereo_after_round_trip("F/C=C/F"), BondStereo::E);
        assert_eq!(stereo_after_round_trip("F/C=C\\F"), BondStereo::Z);
    }

    #[test]
    fn test_e_and_z_survive_from_the_reversed_spelling_too() {
        // The writer and `resolve_directional_bonds` have to agree about the
        // sign, and the reversed form is where disagreement shows: this
        // molecule is written back with both markers flipped, which is only
        // correct if both sides use the same rule.
        assert_eq!(stereo_after_round_trip("C(/F)=C/F"), BondStereo::Z);
        assert_eq!(stereo_after_round_trip("C(\\F)=C/F"), BondStereo::E);
    }

    #[test]
    fn test_writing_is_idempotent_for_stereo() {
        // A second pass must not flip anything. The first round trip may
        // respell the molecule; the spelling it lands on has to be stable.
        for smiles in ["F/C=C/F", "F/C=C\\F", "C(/F)=C/F"] {
            let once = write_smiles_for_molecule(&parse_smiles(smiles).expect(smiles));
            let twice = write_smiles_for_molecule(&parse_smiles(&once).expect(&once));
            assert_eq!(once, twice, "{smiles} is not stable under rewriting");
        }
    }

    #[test]
    fn test_a_double_bond_without_stereo_writes_no_markers() {
        // The common case must not gain slashes.
        let mol = parse_smiles("FC=CF").expect("valid SMILES");
        let written = write_smiles_for_molecule(&mol);
        assert!(!written.contains('/'), "{written}");
        assert!(!written.contains('\\'), "{written}");
    }

    #[test]
    fn test_write_smiles() {
        let atom_symbols: Vec<String> = vec!["c", "c", "c", "c", "c", "c", "O"]
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        let atom_neighbours: Vec<Vec<usize>> = vec![
            vec![1, 5],
            vec![2, 0],
            vec![1, 3],
            vec![2, 4],
            vec![3, 5],
            vec![0, 4, 6],
            vec![5],
        ];
        let atom_rankings: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6];
        let bond_symbols: std::collections::HashMap<String, String> = [
            ("0,1", ""),
            ("0,5", ""),
            ("1,0", ""),
            ("1,2", ""),
            ("2,1", ""),
            ("2,3", ""),
            ("3,2", ""),
            ("3,4", ""),
            ("4,3", ""),
            ("4,5", ""),
            ("5,4", ""),
            ("5,0", ""),
            ("5,6", "="),
            ("6,5", "="),
        ]
        .iter()
        .map(|s| (s.0.to_string(), s.1.to_string()))
        .collect();

        assert_eq!(
            write_smiles(atom_symbols, atom_neighbours, atom_rankings, bond_symbols),
            "c1ccccc1=O".to_string()
        );
    }

    // Round-trips: parse a SMILES string, write it back out, and expect the
    // exact same string. Not every valid SMILES round-trips exactly this
    // way (e.g. stereochemistry markers, or which ring-closure digit gets
    // used when an atom closes more than one ring) -- these are known,
    // pre-existing limitations of this traversal algorithm, not new bugs
    // from this port. Cases exercising those are commented out, same as
    // they were in the original pre-port file.
    #[test]
    fn test_round_trip_smiles() {
        let test_data: Vec<&str> = vec![
            "Oc1ccccc1",
            // Re-enabled by #154: the parser used to mis-tokenize back-to-back
            // bare ring-closure digits ("12") as one ring number 12 rather than
            // two closures. Both this and the fused-ring case below round-trip
            // exactly now, writer output included.
            "Oc1cccc2ccccc12",
            "CC(C)(CCCOc1cc(Cl)c(OCCCC(C)(C)C(=O)O)cc1Cl)C(=O)O", // 4631
            "OCCCCCNCc1c2ccccc2c(CNCCCCCO)c2ccccc12",
            "CC1(C)c2ccc([nH]2)C2(C)CCCCNC(=O)c3cccc(n3)C(=O)NCCCCC(C)(c3ccc1[nH]3)c1ccc([nH]1)C(C)(C)c1ccc2[nH]1", // 4971
            // "C1CC1N1CN2c3nonc3N3CN(C4CC4)CN4c5nonc5N(C1)C2C34" is omitted
            // for the same reason (ends in "C34").
            "BrC1CCC(Br)C(Br)CCC(Br)C(Br)CCC1Br", // 377203
            "OC(=O)c1cc2Cc3cc(Cc4cc(Cc5cc(Cc(c2)c1)cc(c5)C(O)=O)cc(c4)C(O)=O)cc(c3)C(O)=O", // graph reduction demo
            "C1C2CC3CC1CC(C2)C3", // example from nauty, https://pallini.di.uniroma1.it/Introduction.html
            "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 168272 -- exercises bracket-atom charge notation
            "O=P1(O)OC2C3OP(=O)(O)OP(=O)(O)OC3C3OP(=O)(O)OP(=O)(O)OC3C2OP(=O)(O)O1", // 208361
        ];

        for smiles in test_data.iter() {
            let mol = parse_smiles(smiles).unwrap_or_else(|e| {
                panic!("failed to parse {smiles}: {e}");
            });
            let written = write_smiles_for_molecule(&mol);
            assert_eq!(&written, smiles, "round-trip mismatch for {smiles}");
        }
    }
}
