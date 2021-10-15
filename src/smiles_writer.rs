// Copyright 2021 Chiral Ltd.
// Licensed under the Apache license 2.0 (http://www.apache.org/licenses/) 
// This file may not be copied, modified, or distributed
// except according to those terms.

//! 
//! A two-pass traversing algorithm is from the Book [Tutorials in Chemoinformatics](https://www.google.co.jp/books/edition/Tutorials_in_Chemoinformatics/L8toswEACAAJ?hl=en), Page 422 - 425, by Alexandre Varnek
//!
//! # Examples
//! 
//! build SMILES from atom symbols, bond symbols, topological connectivities and atom rankings
//! ```rust
//! let atom_symbols: Vec<String> = vec!["c", "c", "c", "c", "c", "c", "O"].into_iter().map(|s| s.to_string()).collect();
//! let atom_neighbours: Vec<Vec<usize>> = vec![vec![1, 5], vec![2, 0], vec![1, 3], vec![2, 4], vec![3, 5], vec![0, 4, 6], vec![5]];
//! let atom_rankings: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6];
//! let bond_symbols: std::collections::HashMap<String, String> = [("0,1", ""), ("0,5", ""), ("1,0", ""), ("1,2", ""), ("2,1", ""), ("2,3", ""), ("3,2", ""), ("3,4", ""), ("4,3", ""), ("4,5", ""), ("5,4", ""), ("5,0", ""), ("5,6", "="), ("6,5", "=")].iter().map(|s| (s.0.to_string(), s.1.to_string())).collect();
//! assert_eq!(write_smiles(atom_symbols, atom_neighbours, atom_rankings, bond_symbols), "c1ccccc1=O".to_string());
//! ``` 
//!
//! or implement the trait [TraitMoleculeForSMILES](trait.TraitMoleculeForSMILES.html) for your own Molecule Type, e.g. a Molecule Type from crate [purr](https://github.com/rapodaca/purr)
//! ```rust
//! struct Molecule {
//!     pub atoms: Vec<purr::graph::Atom>,
//!     pub bond_table: std::collections::HashMap<String, String>
//! }
//! 
//! impl Molecule {
//!     pub fn from_smiles(smiles: &str) -> Self {
//!         let mut builder = purr::graph::Builder::new();
//! 
//!         match purr::read::read(smiles, &mut builder, None) {
//!             Ok(_) => {
//!                 let mut atoms = builder.build().expect("atoms");
//!                 let mut bond_table = std::collections::HashMap::new();
//!                 for atom_idx in 0..(atoms.len()) {
//!                     for bond in atoms[atom_idx].bonds.iter_mut() {
//!                         bond_table.insert(format!("{},{}", atom_idx, bond.tid), bond.kind.to_string());
//!                     }
//!                 }
//! 
//!                 Self { atoms, bond_table }
//!             },
//!             Err(e) => {
//!                 panic!("error smiles parsing: {:?}", e);
//!             }
//!         }
//!     }
//! }
//! 
//! impl TraitMoleculeForSMILES for Molecule {
//!     fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
//!         self.atoms[*atom].bonds.iter()
//!             .map(|b| b.tid)
//!             .collect()
//!     }
//! 
//!     fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
//!         match self.bond_table.get(&format!("{},{}", atom_1, atom_2)) {
//!             Some(bond_string) => bond_string.clone(),
//!             None => String::from("No such bond")
//!         }
//!     }
//! 
//!     fn get_atom_symbol(&self, atom: &usize) -> String {
//!         self.atoms[*atom].kind.to_string()
//!     }
//! 
//!     fn get_atom_ranking(&self, atom: &usize) -> usize {
//!         atom.clone()
//!     }
//! 
//!     fn count_of_atoms(&self) -> usize {
//!         self.atoms.len()
//!     }
//! }
//! 
//! let smiles = String::from("Oc1ccccc1")
//! let mol = Molecule::from_smiles(&smiles);
//! assert_eq!(write_smiles_for_mol(&mol), smiles);
//! ```

/// Manager of ring digit
struct DigitHeap {
    digits: Vec<usize>
}

impl DigitHeap {
    fn init() -> Self {
        Self { digits: vec![] }
    }

    fn find(&mut self) -> usize {
        let mut digit: usize = 1;
        while digit < 100 {
            if !self.digits.contains(&digit) { break; }
            digit += 1;
        }

        self.digits.push(digit);
        digit
    }

    fn remove(&mut self, digit: &usize)  {
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
    pub dh: DigitHeap
}

impl DataPool {
    fn init() -> Self {
        Self {
            ancestors: vec![],
            visited: vec![],
            opening_closures: OpeningClosures::new(),
            closing_closures: ClosingClosures::new(),
            dh: DigitHeap::init() 
        }
    }
}

/// implement trait TraitMoleculeForSMILES to use [`write_smiles_for_mol`](fn@write_smiles_for_mol)
pub trait TraitMoleculeForSMILES {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize>;
    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String;
    fn get_atom_symbol(&self, atom: &usize) -> String;
    fn get_atom_ranking(&self, atom: &usize) -> usize;
    fn count_of_atoms(&self) -> usize;
}

struct MoleculeForSmiles {
    atom_symbols: Vec<String>,
    atom_neighbours: Vec<Vec<usize>>,
    atom_rankings: Vec<usize>,
    bond_symbols: std::collections::HashMap<String, String>,
}

impl TraitMoleculeForSMILES for MoleculeForSmiles {
    fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
        self.atom_neighbours[*atom].clone()
    }

    fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
        match self.bond_symbols.get(&format!("{},{}", *atom_1, *atom_2)) {
            Some(symbol) => symbol.clone(),
            None => panic!("No such bond")
        }
    }

    fn get_atom_symbol(&self, atom: &usize) -> String {
        self.atom_symbols[*atom].clone()
    }

    fn get_atom_ranking(&self, atom: &usize) -> usize {
        self.atom_rankings[*atom]
    }

    fn count_of_atoms(&self) -> usize {
        self.atom_symbols.len()
    }
}

#[inline(always)]
fn get_neighbours_excluding_parent<T: TraitMoleculeForSMILES>(
    mol: &T,
    atom_current: usize, 
    atom_parent_opt: Option<usize>, 
) -> Vec<usize> {
    mol.get_neighbours_of_atom(&atom_current).clone().into_iter()
        .filter(|&idx| match atom_parent_opt {
            Some(atom_parent) => idx != atom_parent,
            None => true
        })
        .collect()
}

/// A recursive function for detecting opening closures by the 1st traversing
fn get_closures_for_atom<T: TraitMoleculeForSMILES>(
    mol: &T,
    atom_current: usize, 
    atom_parent_opt: Option<usize>, 
    dp: &mut DataPool
) {
    dp.ancestors.push(atom_current);
    dp.visited.push(atom_current);
    
    let mut nbors: Vec<usize> = get_neighbours_excluding_parent(mol, atom_current, atom_parent_opt); 
    nbors.sort_by_key(|idx| mol.get_atom_ranking(idx));

    for nb in nbors.iter() {
        if dp.ancestors.contains(nb) { 
            dp.opening_closures.entry(*nb).or_insert(vec![]).push(atom_current); 
        }
        else {
            if !dp.visited.contains(nb) { 
                get_closures_for_atom(mol, *nb, Some(atom_current), dp); 
            }
        }
    }

    let index = dp.ancestors.iter().position(|x| *x == atom_current).unwrap();
    dp.ancestors.remove(index);
}

/// A recursive function for building smiles by the 2nd traversing
fn build_smiles_for_atom<T: TraitMoleculeForSMILES>(
    mol: &T,
    atom_current: usize,
    atom_parent_opt: Option<usize>,
    dp: &mut DataPool
) -> String {
    dp.visited.push(atom_current);
    let mut seq: String = String::from("");

    match atom_parent_opt {
        Some(atom_parent) => { seq += &mol.get_bond_symbol(&atom_parent, &atom_current); },

        None => {}
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
            let oadts = dp.closing_closures.entry(*oc).or_insert(vec![]);
            oadts.push((atom_current, digit));
        }
    }

    let mut nbors: Vec<usize> = get_neighbours_excluding_parent(mol, atom_current, atom_parent_opt); 
    nbors.sort_by_key(|idx| mol.get_atom_ranking(idx));

    let mut branches: Vec<String> = vec![];
    for n in nbors.iter() {
        if !dp.visited.contains(&n) {
            branches.push(build_smiles_for_atom(mol, *n, Some(atom_current), dp));
        }
    }

    if branches.len() > 1 {
        for branch in branches[..(branches.len()-1)].iter() {
            seq += &format!("({})", branch);
        }
    }

    if branches.len() > 0 {
        seq += &branches[branches.len()-1];
    }

    seq
}

/// SMILES writer with trait Molecule 
pub fn write_smiles_for_mol<T: TraitMoleculeForSMILES>(
    mol: &T
) -> String {
    let mut dp = DataPool::init();

    // find the atom with minimum ranking to start
    let mut atom_indexes: Vec<usize> = (0..mol.count_of_atoms()).collect();
    atom_indexes.sort_by_key(|idx| mol.get_atom_ranking(idx));

    get_closures_for_atom(mol, atom_indexes[0], None, &mut dp);

    dp.visited.clear();
    build_smiles_for_atom(mol, atom_indexes[0], None, &mut dp)
}

/// SMILES writer with raw data
pub fn write_smiles(
    atom_symbols: Vec<String>,
    atom_neighbours: Vec<Vec<usize>>,
    atom_rankings: Vec<usize>,
    bond_symbols: std::collections::HashMap<String, String>,
) -> String {
    let mol = MoleculeForSmiles { atom_symbols, atom_neighbours, bond_symbols, atom_rankings };
    write_smiles_for_mol(&mol)
}


#[cfg(test)]
mod tests {
    use super::*;

    struct Molecule {
        pub atoms: Vec<purr::graph::Atom>,
        pub bond_table: std::collections::HashMap<String, String>
    }
    
    impl Molecule {
        pub fn from_smiles(smiles: &str) -> Self {
            let mut builder = purr::graph::Builder::new();
    
            match purr::read::read(smiles, &mut builder, None) {
                Ok(_) => {
                    let mut atoms = builder.build().expect("atoms");
                    let mut bond_table = std::collections::HashMap::new();
                    for atom_idx in 0..(atoms.len()) {
                        for bond in atoms[atom_idx].bonds.iter_mut() {
                            bond_table.insert(format!("{},{}", atom_idx, bond.tid), bond.kind.to_string());
                        }
                    }
    
                    Self { atoms, bond_table }
                },
                Err(e) => {
                    panic!("error smiles parsing: {:?}", e);
                }
            }
        }
    }

    impl TraitMoleculeForSMILES for Molecule {
        fn get_neighbours_of_atom(&self, atom: &usize) -> Vec<usize> {
            self.atoms[*atom].bonds.iter()
                .map(|b| b.tid)
                .collect()
        }

        fn get_bond_symbol(&self, atom_1: &usize, atom_2: &usize) -> String {
            match self.bond_table.get(&format!("{},{}", atom_1, atom_2)) {
                Some(bond_string) => bond_string.clone(),
                None => String::from("No such bond")
            }
        }

        fn get_atom_symbol(&self, atom: &usize) -> String {
            self.atoms[*atom].kind.to_string()
        }

        fn get_atom_ranking(&self, atom: &usize) -> usize {
            atom.clone()
        }

        fn count_of_atoms(&self) -> usize {
            self.atoms.len()
        }
    }

    #[test]
    fn test_write_smiles_for_mol() {
        type InputType1 = String;
        let test_data: Vec<InputType1> = vec![
            //
            //  Commented failure cases are due to the following two reason
            //      - stereochemistry C@H vs C@@H
            //      - the order of ring digits when an atom closes multiple rings, e.g. Oc1cccc2ccccc12 vs Oc1cccc2ccccc21 
            //
            "Oc1ccccc1",
            "Oc1cccc2ccccc12",
            // "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)-c1cccc(c1C(=O)O)-c1ccc(cc1)C[n+]1ccn(c1)Cc1ccc(cc1)C3=O", // chembl 15,
            "CC(C)(CCCOc1cc(Cl)c(OCCCC(C)(C)C(=O)O)cc1Cl)C(=O)O", // 4631
            // "C[N+](C)(CCCCCC[N+](C)(C)CCCN1C(=O)C2C3c4ccccc4C(c4ccccc43)C2C1=O)CCCN1C(=O)c2ccccc2C1=O", // 6053 
            // "N[C@@H](Cc1cnc(C23CC4CC(CC(C4)C2)C3)[nH]1)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1cnc(C23CC4CC(CC(C4)C2)C3)[nH]1)C(=O)NCc1ccccc1", // 7844 
            "OCCCCCNCc1c2ccccc2c(CNCCCCCO)c2ccccc12", // 23218 
            // "NC[C@@H]1O[C@H](O[C@@H]2[C@@H](CSCCNC(=S)NCCCCn3c(=O)c4ccc5c6ccc7c(=O)n(CCCCNC(=S)NCCSC[C@H]8O[C@@H](O[C@@H]9[C@@H](O)[C@H](N)C[C@H](N)[C@H]9O[C@H]9O[C@H](CN)[C@@H](O)[C@H](O)[C@H]9N)[C@H](O)[C@@H]8O[C@H]8O[C@@H](CN)[C@@H](O)[C@H](O)[C@H]8N)c(=O)c8ccc(c9ccc(c3=O)c4c59)c6c78)O[C@@H](O[C@@H]3[C@@H](O)[C@H](N)C[C@H](N)[C@H]3O[C@H]3O[C@H](CN)[C@@H](O)[C@H](O)[C@H]3N)[C@@H]2O)[C@H](N)[C@@H](O)[C@@H]1O", // 52881 
            "CC1(C)c2ccc([nH]2)C2(C)CCCCNC(=O)c3cccc(n3)C(=O)NCCCCC(C)(c3ccc1[nH]3)c1ccc([nH]1)C(C)(C)c1ccc2[nH]1", // 4971 
            // "O=C1NNC(=O)c2ccccc2SSc2ccccc2C(=O)NNC(=O)c2ccccc2SSc2ccccc21", // 140635
            "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 168272
            "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 171007
            "C1CC1N1CN2c3nonc3N3CN(C4CC4)CN4c5nonc5N(C1)C2C34", // 199821
            "O=P1(O)OC2C3OP(=O)(O)OP(=O)(O)OC3C3OP(=O)(O)OP(=O)(O)OC3C2OP(=O)(O)O1", // 208361
            "CC[n+]1ccc(-c2cc[n+](Cc3cc(C[n+]4ccc(-c5cc[n+](CC)cc5)cc4)cc(C[n+]4ccc(-c5cc[n+](Cc6cc(C[n+]7ccc(-c8cc[n+](Cc9cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)c9)cc8)cc7)cc(-[n+]7ccc(-c8cc[n+](-c9cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)c9)cc8)cc7)c6)cc5)cc4)c3)cc2)cc1", // 826428
            "CC[n+]1ccc(-c2cc[n+](Cc3cc(C[n+]4ccc(-c5cc[n+](CC)cc5)cc4)cc(C[n+]4ccc(-c5cc[n+](Cc6cc(C[n+]7ccc(-c8cc[n+](Cc9cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)c9)cc8)cc7)cc(-[n+]7ccc(-c8cc[n+](-c9cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)c9)cc8)cc7)c6)cc5)cc4)c3)cc2)cc1", // 1246825
            "BrC1CCC(Br)C(Br)CCC(Br)C(Br)CCC1Br", // 377203
            // "C[N+]1(C)CC23c4c5c6c7c8c4c4c2c2c9c%10c%11c%12c%13c9c9c%14c%15c%16c%17c%18c%19c(c8c%17c4c%16c29)C7C2c4c-%19c7c8c9c(c%14c%13c%13c9c9c8c4c4c2c6c2c5c(c=%11c5c2c4c9c5c%12%13)C%103C1)C%15C%187", // CHEMBL415840 
            // "CCC[C@H]1CC[C@H]([C@H]2CC[C@H](OC(=O)[C@H]3[C@@H](c4ccc(O)cc4)[C@H](C(=O)O[C@H]4CC[C@H]([C@H]5CC[C@H](CCC)CC5)CC4)[C@@H]3c3ccc(O)cc3)CC2)CC1", // CHEMBL2348759
            // "OC(c1ccccc1)C1(c2ccccc2)C23c4c5c6c7c8c9c(c%10c%11c2c2c4c4c%12c5c5c6c6c8c8c%13c9c9c%10c%10c%11c%11c2c2c4c4c%12c%12c5c5c6c8c6c8c%13c9c9c%10c%10c%11c2c2c4c4c%12c5c6c5c8c9c%10c2c45)C731", // 408840 
            // "O=C(CCCc1ccc(C2(c3ccccc3)C34c5c6c7c8c9c%10c(c%11c%12c3c3c5c5c%13c6c6c7c7c9c9c%14c%10c%10c%11c%11c%12c%12c3c3c5c5c%13c%13c6c6c7c9c7c9c%14c%10c%10c%11c%11c%12c3c3c5c5c%13c6c7c6c9c%10c%11c3c56)C824)cc1)NC(CO)(CO)CO", // 267348 
            // r#"C[C@H](CC[C@@H]([C@@H]([C@H](C)C[C@H](C(=C)/C(=C/CO)/C)O)O)OS(=O)(=O)[O-])[C@H]([C@@H](C)[C@H]1[C@@H]([C@@H]([C@H]2[C@H](O1)[C@@H](C[C@]3([C@H](O2)C[C@H]4[C@H](O3)C[C@]5([C@H](O4)[C@H]([C@H]6[C@H](O5)C[C@H]([C@H](O6)[C@@H]([C@H](C[C@H]7[C@@H]([C@@H]([C@H]8[C@H](O7)C[C@H]9[C@H](O8)C[C@H]1[C@H](O9)[C@H]([C@@H]2[C@@H](O1)[C@@H]([C@H]([C@@H](O2)[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@@H](C[C@H](C[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@H]1[C@@H](C[C@]2([C@H](O1)[C@@H]([C@]1([C@H](O2)C[C@]2([C@H](O1)CC[C@]1([C@H](O2)C[C@]2([C@H](O1)C[C@H]1[C@H](O2)CC[C@H](O1)[C@]1([C@@H](C[C@H]2[C@](O1)(C[C@H]1[C@](O2)(CC[C@]2([C@H](O1)C[C@H]1[C@](O2)(C[C@H]2[C@H](O1)C/C=C\[C@H]1[C@H](O2)C[C@H]2[C@](O1)(C[C@]1([C@H](O2)C[C@H]2[C@](O1)(CC[C@H](O2)[C@H]([C@@H](C[C@@H](C)[C@@H](C)CC=C)O)O)C)C)C)C)C)C)C)O)C)C)C)C)C)O)C)O)O)O)O)O)O)O)O)O)O)O)O)O)OS(=O)(=O)[O-])O)O)O)O)C)C)O)O)O)O"#, // Maitotoxin
            "OC(=O)c1cc2Cc3cc(Cc4cc(Cc5cc(Cc(c2)c1)cc(c5)C(O)=O)cc(c4)C(O)=O)cc(c3)C(O)=O", // graph reduction demo
            "C1C2CC3CC1CC(C2)C3", // example from nauty, https://pallini.di.uniroma1.it/Introduction.html
       ].into_iter().map(|s| s.to_string()).collect();

       for td in test_data.iter() {
            let smiles = td.clone();
            let mol = Molecule::from_smiles(&smiles);
            assert_eq!(write_smiles_for_mol(&mol), smiles);
       }
    }

    #[test]
    fn test_write_smiles() {
        type InputType1 = Vec<String>;
        type InputType2 = Vec<Vec<usize>>;
        type InputType3 = Vec<usize>;
        type InputType4 = std::collections::HashMap<String, String>;
        type OutputType = String;

        let test_data: Vec<(InputType1, InputType2, InputType3, InputType4, OutputType)> = vec![
            (
                vec!["c", "c", "c", "c", "c", "c", "O"].into_iter().map(|s| s.to_string()).collect(),
                vec![vec![1, 5], vec![2, 0], vec![1, 3], vec![2, 4], vec![3, 5], vec![0, 4, 6], vec![5]],
                vec![0, 1, 2, 3, 4, 5, 6],
                [("0,1", ""), ("0,5", ""), ("1,0", ""), ("1,2", ""), ("2,1", ""), ("2,3", ""), ("3,2", ""), ("3,4", ""), ("4,3", ""), ("4,5", ""), ("5,4", ""), ("5,0", ""), ("5,6", "="), ("6,5", "=")].iter().map(|s| (s.0.to_string(), s.1.to_string())).collect(),
                "c1ccccc1=O".to_string()
            )
        ];


       for td in test_data.iter() {
            let (atom_symbols, atom_neighbours, atom_rankings, bond_symbols, smiles) = td.clone();
            assert_eq!(write_smiles(atom_symbols, atom_neighbours, atom_rankings, bond_symbols), smiles);
       }
    }

    #[test]
    fn name() {
    }
}
