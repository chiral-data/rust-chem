use chem::core::molecule::Molecule;
use chem::io::reader::{self, ReadOutcome};

/// Which file format a loaded dataset came from.
///
/// An alias rather than a type of its own: the app used to keep a parallel enum
/// beside `chem::io`'s parsers, and two enums meaning the same thing is one more
/// place for the app and the CLI to disagree about what `.txt` means.
pub type DatasetFormat = reader::Format;

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

    /// Builds a dataset from what [`chem::io::reader`] read.
    ///
    /// The parallel-vector shape is what the views consume; the skipped records
    /// are the caller's to report, which is why they are not swallowed here.
    pub fn from_outcome(outcome: &ReadOutcome) -> Self {
        let mut dataset = Self::new();
        for record in &outcome.records {
            dataset.molecules.push(record.molecule.clone());
            // SDF carries coordinates and connectivity rather than a SMILES
            // string. The placeholder is a display decision, so it is made
            // here rather than in the library that read the file.
            dataset
                .smiles
                .push(record.smiles.clone().unwrap_or_else(|| "(SDF)".to_owned()));
            dataset.names.push(record.name.clone());
        }
        dataset
    }

    pub fn len(&self) -> usize {
        self.molecules.len()
    }

    pub fn is_empty(&self) -> bool {
        self.molecules.is_empty()
    }

    /// The built-in sampler, as a SMILES file rather than a third parse loop.
    ///
    /// Written as the text a user could have supplied so it goes through
    /// exactly the path a loaded file does — anything that breaks the examples
    /// breaks real files too, and is therefore visible immediately.
    pub fn example_dataset() -> Self {
        const EXAMPLES: &str = "\
C Methane
CC Ethane
CCC Propane
C=C Ethene
C#C Ethyne
c1ccccc1 Benzene
CC(C)C Isobutane
CCO Ethanol
CC(=O)O Acetic_acid
c1ccccc1O Phenol
c1ccc(O)cc1 Phenol_alt
c1ccccc1N Aniline
CC(=O)c1ccccc1 Acetophenone
c1ccc2ccccc2c1 Naphthalene
CC(C)(C)C Neopentane
";
        let outcome = reader::read_smiles(EXAMPLES);
        for skipped in &outcome.skipped {
            log::warn!(
                "Built-in example on line {} failed to parse: {}",
                skipped.position,
                skipped.error
            );
        }
        Self::from_outcome(&outcome)
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
        let outcome = reader::read_smiles(&content);
        assert!(outcome.skipped.is_empty(), "fixture should parse cleanly");
        MoleculeDataset::from_outcome(&outcome)
    }

    #[test]
    fn test_entries_report_name_format_and_size() {
        let mut files = LoadedFiles::new(
            "first.smi".to_string(),
            dataset_of(&["C", "CC"]),
            DatasetFormat::SMILES,
        );
        files.add_and_activate(
            "second.sdf".to_string(),
            dataset_of(&["c1ccccc1"]),
            DatasetFormat::SDF,
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
            LoadedFiles::new("a".to_string(), dataset_of(&["C"]), DatasetFormat::SMILES);
        files.add_and_activate("b".to_string(), dataset_of(&["CC"]), DatasetFormat::SMILES);
        files.activate(0);

        // Switching back must not reorder the list under the user.
        let names: Vec<&str> = files.entries().iter().map(|e| e.name.as_str()).collect();
        assert_eq!(names, vec!["a", "b"]);
        assert_eq!(files.active_index(), 0);
    }

    /// Three entries whose molecule counts identify them: 1, 2 and 3.
    fn three_files() -> LoadedFiles {
        let mut files =
            LoadedFiles::new("a".to_string(), dataset_of(&["C"]), DatasetFormat::SMILES);
        files.add_and_activate(
            "b".to_string(),
            dataset_of(&["C", "CC"]),
            DatasetFormat::SMILES,
        );
        files.add_and_activate(
            "c".to_string(),
            dataset_of(&["C", "CC", "CCC"]),
            DatasetFormat::SMILES,
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
            DatasetFormat::SMILES,
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
            LoadedFiles::new("a".to_string(), dataset_of(&["C"]), DatasetFormat::SMILES);
        files.add_and_activate(
            "a".to_string(),
            dataset_of(&["C", "CC", "CCC"]),
            DatasetFormat::SMILES,
        );

        assert_eq!(files.entries().len(), 1);
        assert_eq!(files.entries()[0].dataset.len(), 3);
    }
}
