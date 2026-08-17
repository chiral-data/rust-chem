use chemcore::molecule::Molecule;
use chemio::sdf::parse_sdf;
use chemio::smiles::parse_smiles;

/// Which file format a loaded dataset came from. Doesn't affect anything
/// downstream (fingerprinting/search work the same either way) — it's
/// tracked per `LoadedFile` purely so the Data window can label entries and
/// status messages accurately.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DatasetFormat {
    Smiles,
    Sdf,
}

impl DatasetFormat {
    /// Picks a format from a filename's extension. Defaults to SMILES for
    /// anything not recognized as SDF, matching the loader's prior behavior
    /// (`.smi`/`.smiles`/`.txt`/no extension all went through the SMILES path).
    pub fn from_filename(name: &str) -> Self {
        if name.to_lowercase().ends_with(".sdf") {
            DatasetFormat::Sdf
        } else {
            DatasetFormat::Smiles
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            DatasetFormat::Smiles => "SMILES",
            DatasetFormat::Sdf => "SDF",
        }
    }
}

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

    /// Parses a multi-molecule SDF dataset from already-loaded file content.
    /// `chemio::sdf::parse_sdf` only parses a single `$$$$`-terminated
    /// record, so this splits the file into records first (real-world SDF
    /// files are almost always a batch of these), mirroring how
    /// `load_from_smiles_str` splits its input into one molecule per line.
    pub fn load_from_sdf_str(content: &str) -> anyhow::Result<Self> {
        let mut dataset = Self::new();
        let mut record_lines: Vec<&str> = Vec::new();
        let mut record_num = 0;

        for line in content.lines() {
            record_lines.push(line);
            if line.trim() == "$$$$" {
                record_num += 1;
                dataset.push_sdf_record(&record_lines, record_num);
                record_lines.clear();
            }
        }
        // A trailing record with no `$$$$` terminator (e.g. a single-molecule
        // file) still has real content sitting in `record_lines`.
        if record_lines.iter().any(|line| !line.trim().is_empty()) {
            record_num += 1;
            dataset.push_sdf_record(&record_lines, record_num);
        }

        log::info!("Loaded {} molecules from SDF file", dataset.molecules.len());
        Ok(dataset)
    }

    fn push_sdf_record(&mut self, lines: &[&str], record_num: usize) {
        let record = lines.join("\n");
        match parse_sdf(&record) {
            Ok(mol) => {
                let name = mol
                    .name()
                    .map(str::to_owned)
                    .unwrap_or_else(|| format!("Molecule_{}", record_num));
                self.molecules.push(mol);
                // SDF stores 2D/3D coordinates and connectivity, not a SMILES
                // string, so there's nothing to put here — a placeholder
                // makes that clear rather than showing a blank code snippet.
                self.smiles.push("(SDF)".to_string());
                self.names.push(name);
            }
            Err(e) => {
                log::warn!("Failed to parse SDF record {}: {}", record_num, e);
            }
        }
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
    pub format: DatasetFormat,
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
    pub fn new(initial_name: String, initial: MoleculeDataset, format: DatasetFormat) -> Self {
        Self {
            entries: vec![LoadedFile {
                name: initial_name,
                dataset: initial,
                format,
            }],
            active: 0,
        }
    }

    /// Adds a new entry and makes it active, or, if an entry with this name
    /// already exists (e.g. reloading the same file), replaces its dataset
    /// in place and activates that instead of appending a duplicate.
    pub fn add_and_activate(
        &mut self,
        name: String,
        dataset: MoleculeDataset,
        format: DatasetFormat,
    ) {
        if let Some(idx) = self.entries.iter().position(|e| e.name == name) {
            self.entries[idx].dataset = dataset;
            self.entries[idx].format = format;
            self.active = idx;
        } else {
            self.entries.push(LoadedFile {
                name,
                dataset,
                format,
            });
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

    pub fn active_dataset_mut(&mut self) -> &mut MoleculeDataset {
        &mut self.entries[self.active].dataset
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|e| e.name.as_str())
    }
}
