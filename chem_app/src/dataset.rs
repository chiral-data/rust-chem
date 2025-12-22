use chemcore::molecule::Molecule;
use chemio::smiles::parse_smiles;
use std::path::Path;

#[derive(Clone)]
pub struct MoleculeDataset {
    pub molecules: Vec<Molecule>,
    pub smiles: Vec<String>,
    pub names: Vec<String>,
}

impl MoleculeDataset {
    pub fn new() -> Self {
        Self {
            molecules: Vec::new(),
            smiles: Vec::new(),
            names: Vec::new(),
        }
    }

    pub fn load_from_smiles_file<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let mut dataset = Self::new();

        for (line_num, line) in content.lines().enumerate() {
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }

            let smiles_str = parts[0];
            let name = if parts.len() > 1 {
                parts[1..].join(" ")
            } else {
                format!("Molecule_{}", line_num + 1)
            };

            match parse_smiles(smiles_str) {
                Ok(mol) => {
                    dataset.molecules.push(mol);
                    dataset.smiles.push(smiles_str.to_string());
                    dataset.names.push(name);
                }
                Err(e) => {
                    log::warn!(
                        "Failed to parse SMILES on line {}: {} - Error: {}",
                        line_num + 1,
                        smiles_str,
                        e
                    );
                }
            }
        }

        log::info!(
            "Loaded {} molecules from SMILES file",
            dataset.molecules.len()
        );
        Ok(dataset)
    }

    pub fn len(&self) -> usize {
        self.molecules.len()
    }

    pub fn is_empty(&self) -> bool {
        self.molecules.is_empty()
    }

    pub fn example_dataset() -> anyhow::Result<Self> {
        let examples = vec![
            ("C", "Methane"),
            ("CC", "Ethane"),
            ("CCC", "Propane"),
            ("C=C", "Ethene"),
            ("C#C", "Ethyne"),
            ("c1ccccc1", "Benzene"),
            ("CC(C)C", "Isobutane"),
            ("CCO", "Ethanol"),
            ("CC(=O)O", "Acetic_acid"),
            ("c1ccccc1O", "Phenol"),
            ("c1ccc(O)cc1", "Phenol_alt"),
            ("c1ccccc1N", "Aniline"),
            ("CC(=O)c1ccccc1", "Acetophenone"),
            ("c1ccc2ccccc2c1", "Naphthalene"),
            ("CC(C)(C)C", "Neopentane"),
        ];

        let mut dataset = Self::new();
        for (smiles, name) in examples {
            match parse_smiles(smiles) {
                Ok(mol) => {
                    dataset.molecules.push(mol);
                    dataset.smiles.push(smiles.to_string());
                    dataset.names.push(name.to_string());
                }
                Err(e) => {
                    log::warn!("Failed to parse example {}: {}", name, e);
                }
            }
        }

        Ok(dataset)
    }
}

impl Default for MoleculeDataset {
    fn default() -> Self {
        Self::new()
    }
}
