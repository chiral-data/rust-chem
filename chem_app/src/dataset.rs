use chemcore::molecule::Molecule;
use chemio::smiles::parse_smiles;

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

    /// Parses a SMILES dataset from already-loaded file content. Takes content
    /// instead of a path so callers can supply bytes from any source (native
    /// filesystem, browser file picker, etc.) without this function needing
    /// to know how they were obtained.
    pub fn load_from_smiles_str(content: &str) -> anyhow::Result<Self> {
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

pub struct LoadedFile {
    pub name: String,
    pub dataset: MoleculeDataset,
}

/// Every dataset loaded this session, and which one is currently active.
/// Loading a file or the example set adds an entry rather than replacing
/// the previous one, so the Data window can switch back to something
/// loaded earlier instead of losing it.
pub struct LoadedFiles {
    entries: Vec<LoadedFile>,
    active: usize,
}

impl LoadedFiles {
    pub fn new(initial_name: String, initial: MoleculeDataset) -> Self {
        Self {
            entries: vec![LoadedFile {
                name: initial_name,
                dataset: initial,
            }],
            active: 0,
        }
    }

    /// Adds a new entry and makes it active, or, if an entry with this name
    /// already exists (e.g. reloading the same file), replaces its dataset
    /// in place and activates that instead of appending a duplicate.
    pub fn add_and_activate(&mut self, name: String, dataset: MoleculeDataset) {
        if let Some(idx) = self.entries.iter().position(|e| e.name == name) {
            self.entries[idx].dataset = dataset;
            self.active = idx;
        } else {
            self.entries.push(LoadedFile { name, dataset });
            self.active = self.entries.len() - 1;
        }
    }

    pub fn activate(&mut self, index: usize) {
        if index < self.entries.len() {
            self.active = index;
        }
    }

    pub fn active_index(&self) -> usize {
        self.active
    }

    pub fn active_dataset(&self) -> &MoleculeDataset {
        &self.entries[self.active].dataset
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|e| e.name.as_str())
    }
}
