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

use crate::core::bond::BondOrder;
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
    let atom = mol.atom(atom_idx);
    let symbol = atom.element().symbol();
    let symbol = if atom.is_aromatic() {
        symbol.to_lowercase()
    } else {
        symbol.to_string()
    };

    let needs_brackets =
        atom.formal_charge() != 0 || atom.isotope().is_some() || atom.explicit_hydrogens() > 0;
    if !needs_brackets {
        return symbol;
    }

    let mut bracket = String::from("[");
    if let Some(isotope) = atom.isotope() {
        bracket += &isotope.to_string();
    }
    bracket += &symbol;
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

impl TraitMoleculeForSMILES for Molecule {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
        self.neighbors(*atom).iter().map(|n| n.atom_idx).collect()
    }

    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
        match self.graph().get_bond(*atom_1, *atom_2) {
            Some(bond_idx) => bond_order_symbol(self.bond(bond_idx).order()).to_string(),
            None => panic!("No such bond"),
        }
    }

    fn get_atom_symbol(&self, atom: &usize) -> String {
        atom_symbol(self, *atom)
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

    seq += &mol.get_atom_symbol(&atom_current);

    if let Some(oadts) = dp.closing_closures.get_mut(&atom_current) {
        oadts.sort_by_key(|oadt| oadt.1); // close multiple rings, start from smaller digits
        for oadt in oadts.iter() {
            seq += &mol.get_bond_symbol(&atom_current, &oadt.0);
            if oadt.1 > 9 {
                seq += "%";
            }
            seq += &oadt.1.to_string();
            dp.dh.remove(&oadt.1);
        }
    }

    if let Some(ocs) = dp.opening_closures.get(&atom_current) {
        for oc in ocs.iter() {
            let digit = dp.dh.find();
            if digit > 9 {
                seq += "%";
            }
            seq += &digit.to_string();
            let oadts = dp.closing_closures.entry(*oc).or_default();
            oadts.push((atom_current, digit));
        }
    }

    let mut nbors: Vec<usize> = get_neighbours_excluding_parent(mol, atom_current, atom_parent_opt);
    nbors.sort_by_key(|idx| mol.get_atom_ranking(idx, rankings));

    let mut branches: Vec<String> = vec![];
    for n in nbors.iter() {
        if !dp.visited.contains(n) {
            branches.push(build_smiles_for_atom(
                mol,
                rankings,
                *n,
                Some(atom_current),
                dp,
            ));
        }
    }

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
pub fn write_smiles_for_mol<T: TraitMoleculeForSMILES>(mol: &T, rankings: &[usize]) -> String {
    let mut dp = DataPool::init();

    // find the atom with minimum ranking to start
    let mut atom_indexes: Vec<usize> = (0..mol.count_of_atoms()).collect();
    atom_indexes.sort_by_key(|idx| mol.get_atom_ranking(idx, rankings));

    get_closures_for_atom(mol, rankings, atom_indexes[0], None, &mut dp);

    dp.visited.clear();
    build_smiles_for_atom(mol, rankings, atom_indexes[0], None, &mut dp)
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
