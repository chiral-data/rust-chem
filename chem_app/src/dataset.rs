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

    /// Every loaded dataset, in load order.
    ///
    /// [`LoadedFiles::names`] was enough while the list showed only names.
    /// Showing what each entry *is* — its format, how many molecules it holds —
    /// needs the entries themselves, and both facts are already here.
    pub fn entries(&self) -> &[LoadedFile] {
        &self.entries
    }

    /// Whether anything can be removed.
    ///
    /// False at one entry. [`LoadedFiles::active_dataset`] indexes
    /// `entries[active]` and every view calls it without checking, so an empty
    /// list would panic — and nothing can produce one anyway, since the app
    /// loads the examples at startup. Making emptiness representable would mean
    /// changing every caller for a state that does not occur.
    pub fn can_remove(&self) -> bool {
        self.entries.len() > 1
    }

    /// Removes a loaded dataset, returning whether the *active* one changed.
    ///
    /// The return value is what the caller needs: fingerprints, results, open
    /// detail windows and row indices all belong to whichever dataset is active,
    /// so they have to be discarded when a different one takes over — and left
    /// alone when it doesn't.
    ///
    /// Removing an entry *before* the active one is the case worth care. Every
    /// entry after it shifts down a slot, so leaving `active` where it was would
    /// silently point it at a different dataset — the same index, a different
    /// molecule set, and nothing invalidated because as far as the caller knows
    /// nothing moved.
    pub fn remove(&mut self, index: usize) -> bool {
        if index >= self.entries.len() || !self.can_remove() {
            return false;
        }

        self.entries.remove(index);

        if index < self.active {
            // Same dataset, new index.
            self.active -= 1;
            false
        } else if index > self.active {
            // Untouched.
            false
        } else {
            // The active entry itself went. Whatever slid into its slot takes
            // over, or the new last entry if it was at the end.
            self.active = self.active.min(self.entries.len() - 1);
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dataset_of(smiles: &[&str]) -> MoleculeDataset {
        let content = smiles.join("\n");
        MoleculeDataset::load_from_smiles_str(&content).expect("valid SMILES")
    }

    #[test]
    fn test_entries_report_name_format_and_size() {
        let mut files = LoadedFiles::new(
            "first.smi".to_string(),
            dataset_of(&["C", "CC"]),
            DatasetFormat::Smiles,
        );
        files.add_and_activate(
            "second.sdf".to_string(),
            dataset_of(&["c1ccccc1"]),
            DatasetFormat::Sdf,
        );

        // The list shows all three of these per entry; before `entries` it
        // could only show the name.
        let described: Vec<(&str, &str, usize)> = files
            .entries()
            .iter()
            .map(|e| (e.name.as_str(), e.format.label(), e.dataset.len()))
            .collect();

        assert_eq!(
            described,
            vec![("first.smi", "SMILES", 2), ("second.sdf", "SDF", 1)]
        );
    }

    #[test]
    fn test_entries_are_in_load_order_not_activation_order() {
        let mut files =
            LoadedFiles::new("a".to_string(), dataset_of(&["C"]), DatasetFormat::Smiles);
        files.add_and_activate("b".to_string(), dataset_of(&["CC"]), DatasetFormat::Smiles);
        files.activate(0);

        // Switching back must not reorder the list under the user.
        let names: Vec<&str> = files.entries().iter().map(|e| e.name.as_str()).collect();
        assert_eq!(names, vec!["a", "b"]);
        assert_eq!(files.active_index(), 0);
    }

    /// Three entries whose molecule counts identify them: 1, 2 and 3.
    fn three_files() -> LoadedFiles {
        let mut files =
            LoadedFiles::new("a".to_string(), dataset_of(&["C"]), DatasetFormat::Smiles);
        files.add_and_activate(
            "b".to_string(),
            dataset_of(&["C", "CC"]),
            DatasetFormat::Smiles,
        );
        files.add_and_activate(
            "c".to_string(),
            dataset_of(&["C", "CC", "CCC"]),
            DatasetFormat::Smiles,
        );
        files
    }

    #[test]
    fn test_removing_before_the_active_entry_keeps_the_same_dataset_active() {
        let mut files = three_files();
        files.activate(2); // "c", three molecules

        let active_changed = files.remove(0);

        // This is the quiet bug in the feature: every entry after the removed
        // one shifts down a slot, so an unadjusted index would now point at "b"
        // while the caller was told nothing changed.
        assert!(!active_changed);
        assert_eq!(files.active_index(), 1);
        assert_eq!(files.names().nth(files.active_index()), Some("c"));
        assert_eq!(files.active_dataset().len(), 3);
    }

    #[test]
    fn test_removing_after_the_active_entry_changes_nothing() {
        let mut files = three_files();
        files.activate(0);

        let active_changed = files.remove(2);

        assert!(!active_changed);
        assert_eq!(files.active_index(), 0);
        assert_eq!(files.active_dataset().len(), 1);
    }

    #[test]
    fn test_removing_the_active_entry_promotes_the_one_after_it() {
        let mut files = three_files();
        files.activate(1); // "b"

        let active_changed = files.remove(1);

        // "c" slid into the slot, so it takes over — and the caller is told,
        // because everything derived from "b" is now stale.
        assert!(active_changed);
        assert_eq!(files.names().nth(files.active_index()), Some("c"));
        assert_eq!(files.active_dataset().len(), 3);
    }

    #[test]
    fn test_removing_the_active_last_entry_falls_back_to_the_new_last() {
        let mut files = three_files();
        files.activate(2); // "c", at the end

        let active_changed = files.remove(2);

        // Nothing slid into the slot, so the index has to come back rather than
        // point past the end.
        assert!(active_changed);
        assert_eq!(files.active_index(), 1);
        assert_eq!(files.names().nth(1), Some("b"));
    }

    #[test]
    fn test_the_last_entry_cannot_be_removed() {
        let mut files = LoadedFiles::new(
            "only".to_string(),
            dataset_of(&["C"]),
            DatasetFormat::Smiles,
        );

        assert!(!files.can_remove());
        assert!(!files.remove(0));
        // `active_dataset` indexes `entries[active]`, and every view calls it
        // without checking, so an empty list would panic rather than show
        // nothing.
        assert_eq!(files.entries().len(), 1);
        assert_eq!(files.active_dataset().len(), 1);
    }

    #[test]
    fn test_removing_out_of_range_is_ignored() {
        let mut files = three_files();
        assert!(!files.remove(99));
        assert_eq!(files.entries().len(), 3);
    }

    #[test]
    fn test_removing_down_to_one_then_stopping() {
        let mut files = three_files();
        files.activate(0);

        assert!(!files.remove(2));
        assert!(!files.remove(1));
        assert!(!files.can_remove());
        assert!(!files.remove(0));

        assert_eq!(files.entries().len(), 1);
        assert_eq!(files.names().next(), Some("a"));
    }

    #[test]
    fn test_reloading_a_name_replaces_it_rather_than_appending() {
        let mut files =
            LoadedFiles::new("a".to_string(), dataset_of(&["C"]), DatasetFormat::Smiles);
        files.add_and_activate(
            "a".to_string(),
            dataset_of(&["C", "CC", "CCC"]),
            DatasetFormat::Smiles,
        );

        assert_eq!(files.entries().len(), 1);
        assert_eq!(files.entries()[0].dataset.len(), 3);
    }
}
