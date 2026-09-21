use chem::core::molecule::Molecule;
use chem::core::trajectory::Trajectory;
use chem::io::format::Carries;
use chem::io::reader::{self, Payload, ReadOutcome, Record};
use chem::io::smiles_writer::write_smiles_for_molecule_canonical;

/// Which file format a loaded dataset came from.
///
/// An alias rather than a type of its own: the app used to keep a parallel enum
/// beside `chem::io`'s parsers, and two enums meaning the same thing is one more
/// place for the app and the CLI to disagree about what `.txt` means.
pub type DatasetFormat = reader::Format;

pub struct MoleculeDataset {
    pub molecules: Vec<Molecule>,
    pub smiles: Vec<String>,
    pub names: Vec<String>,
    /// Which entries of `smiles` this app wrote rather than read.
    ///
    /// [`chem::io::reader::Record::smiles`] is `None` for every format that is
    /// not itself a SMILES string, and its contract says so deliberately: a
    /// library should not invent text for a value it does not have. The app
    /// does invent it -- so the app also says which ones, and the views mark
    /// them.
    pub generated: Vec<bool>,
    /// What the loaded file held, when it wasn't molecules at all (#342).
    ///
    /// `None` for every `Kind::Molecules` format -- the four fields above
    /// keep meaning exactly what they always have. `Some` means
    /// `molecules`/`smiles`/`names`/`generated` are all empty not because
    /// the file was empty, but because it held a trajectory, a volume, a
    /// mesh or a table instead -- a fact a view must check before it
    /// reports "0 molecules" for a 50MB DCD.
    pub non_molecule: Option<NonMoleculeRecord>,
}

/// What a loaded record was, when it wasn't a molecule at all (#342).
///
/// Deliberately not a rendering: #307's own milestone-wide rule is that
/// this wave makes a viewer's formats *readable*, not drawn. Every field
/// here is a fact a view can print, never a picture -- the same "a dash,
/// not a placeholder" honesty #283 already established for a molecule row
/// with no bond model, generalised to a record with no molecule at all.
pub enum NonMoleculeRecord {
    /// One shared topology, many frames -- not `MoleculeDataset`-shaped at
    /// all (a frame carries positions/velocities riding on one topology,
    /// never its own bonds or elements), so this gets its own summary
    /// rather than a row, or ten thousand of them.
    Trajectory {
        atom_count: usize,
        frame_count: usize,
        has_cell: bool,
        /// Kept live, not just summarised, so a view's frame slider can
        /// seek and show a specific frame's own facts on demand. Boxed:
        /// `Trajectory` makes this variant far larger than the others
        /// (`clippy::large_enum_variant`), and every other variant would
        /// pay for that space regardless of which one is active.
        trajectory: Box<Trajectory>,
        /// The frame last sought to.
        selected_frame: usize,
        /// `describe_frame`'s text for `selected_frame`, cached rather than
        /// recomputed on every repaint -- `Trajectory::frame` can mean real
        /// I/O, and egui redraws far more often than the slider moves.
        /// [`MoleculeDataset::seek_trajectory_frame`] is the only thing
        /// that updates it.
        current_frame_summary: String,
    },
    Volume {
        dims: [usize; 3],
        has_cell: bool,
        has_atoms: bool,
    },
    Mesh {
        vertex_count: usize,
        face_count: usize,
    },
    Table {
        row_count: usize,
        columns: Vec<String>,
    },
}

/// Whether a SMILES written from what this format's reader produced is
/// chemistry rather than a plausible-looking fragment list.
///
/// The registry answers it: a format needs bonds *and* an aromatic model, so
/// `bonds | aromaticity` selects exactly SMILES, CXSMILES, SDF, Mol2, CML and
/// commonchem. Measured, benzene and aspirin written to each format and read
/// back:
///
/// | source | bonds read back | written back as |
/// |---|---|---|
/// | mol2, cml, sdf, commonchem | all | `c1ccccc1`, `CC(=O)Oc1ccccc1C(=O)O` |
/// | pdb, full `CONECT` | all | `[C]1[C][C][C][C][C]1` -- adjacency, no order |
/// | pdbqt (aspirin) | 3 of 13 | `[C].[c].[c]...[C][c].[C]O[c]` |
/// | xyz, mmcif, gro | none | `[C].[C].[C].[C].[C].[C]` |
///
/// So the question is not whether *this molecule* has bonds, which #283
/// suggested and which is wrong in both directions: a real PDB
/// (`corpus/pdb/dipeptide-with-ligand.pdb`, 8 atoms and one `CONECT` pair)
/// passes it and writes `[C].[C].[C].[N].[N].[O].[C][O]`, while `[Na+].[Cl-]`
/// from an SDF fails it and would have its correct answer suppressed.
///
/// `bonds` alone is not enough because PDB claims it and states no bond order;
/// `aromaticity` alone is not enough because PDBQT claims it and carries no
/// bonds. If a third caller ever needs this it belongs in the registry, and
/// the honest form there is splitting `BOND_ORDER` out of `BONDS` the way #257
/// split `BONDS` out of `TOPOLOGY`.
/// Above this many atoms the column keeps the placeholder instead of writing a
/// SMILES.
///
/// Two reasons, and either would do. **Cost**: canonical ranking individualises
/// tied atoms in a loop, which is quadratic — measured on a symmetric ring,
/// native release, one record: 0.02s at 500 atoms, 0.07s at 1000, 0.30s at
/// 2000, 1.20s at 4000. Loading is synchronous, and on the web build that is
/// the only thread there is, so a docking receptor would freeze the tab rather
/// than load slowly. Record *count* is not the problem: 2000 drug-like records
/// together cost 0.02s.
///
/// **Legibility**: a five-hundred-character string in a table cell is not
/// information. 500 atoms covers ligands and peptides, which is what the
/// column is for.
const GENERATE_UP_TO_ATOMS: usize = 500;

pub(crate) fn plural(n: usize, word: &str) -> String {
    if n == 1 {
        word.to_string()
    } else {
        format!("{word}s")
    }
}

fn states_a_bond_model(format: DatasetFormat) -> bool {
    format
        .carries()
        .contains(Carries::BONDS.or(Carries::AROMATICITY))
}

/// The facts a trajectory's frame slider shows for whichever frame it last
/// sought to (#342) -- step/time where the format states them, the frame's
/// own cell, and a position bounding box, cheap to compute and the concrete
/// evidence a seek actually moved rather than a redraw of the same frame.
/// Never a picture, per #307's own "readable, not drawn" rule for this wave.
fn describe_frame(frame: &chem::core::trajectory::Frame) -> String {
    let mut parts = Vec::new();
    if let Some(step) = frame.step {
        parts.push(format!("step {step}"));
    }
    if let Some(time) = frame.time {
        parts.push(format!("time {time:.3}"));
    }
    if let Some(cell) = &frame.cell {
        parts.push(format!("cell {:.2}x{:.2}x{:.2} Å", cell.a, cell.b, cell.c));
    }
    if let Some((lo, hi)) = frame.positions.split_first().map(|(first, rest)| {
        rest.iter().fold((*first, *first), |(lo, hi), p| {
            (
                chem::core::geometry::Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
                chem::core::geometry::Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
            )
        })
    }) {
        parts.push(format!(
            "bounding box ({:.2}, {:.2}, {:.2}) to ({:.2}, {:.2}, {:.2})",
            lo.x, lo.y, lo.z, hi.x, hi.y, hi.z
        ));
    }
    if parts.is_empty() {
        "no per-frame facts stated".to_string()
    } else {
        parts.join(", ")
    }
}

impl MoleculeDataset {
    pub fn new() -> Self {
        Self {
            molecules: Vec::new(),
            smiles: Vec::new(),
            names: Vec::new(),
            generated: Vec::new(),
            non_molecule: None,
        }
    }

    /// Builds a dataset from what [`chem::io::reader`] read.
    ///
    /// The parallel-vector shape is what the views consume; the skipped records
    /// are the caller's to report, which is why they are not swallowed here.
    ///
    /// Takes `outcome` by value, not `&ReadOutcome`: a `Kind::Frames` record's
    /// [`Trajectory`] has to be *moved* out of its `Payload` to stay seekable
    /// afterwards (#342) -- the borrow-returning [`Record::trajectory`] can't
    /// give that. A caller that still needs `outcome.skipped` afterwards reads
    /// it before calling this, not after.
    pub fn from_outcome(outcome: ReadOutcome, format: DatasetFormat) -> Self {
        let mut dataset = Self::new();
        for record in outcome.records {
            let Record {
                payload,
                name,
                smiles,
            } = record;
            let molecule = match payload {
                Payload::Molecule(m) => m,
                // #310 widened `Record` to hold a payload that is not
                // necessarily a molecule. `MoleculeDataset` itself stays
                // molecule-shaped -- these four kinds get their own summary
                // instead of a row (#342), the same "readable, not drawn"
                // rule #307 states for the whole milestone.
                Payload::Frames(mut trajectory) => {
                    let atom_count = trajectory.num_atoms();
                    let frame_count = trajectory.frame_count();
                    let frame0 = trajectory.frame(0).ok();
                    let has_cell = frame0.as_ref().is_some_and(|f| f.cell.is_some());
                    let current_frame_summary = frame0
                        .as_ref()
                        .map(describe_frame)
                        .unwrap_or_else(|| "frame 0 could not be read".to_string());
                    dataset.non_molecule = Some(NonMoleculeRecord::Trajectory {
                        atom_count,
                        frame_count,
                        has_cell,
                        trajectory: Box::new(trajectory),
                        selected_frame: 0,
                        current_frame_summary,
                    });
                    continue;
                }
                Payload::Volume(grid) => {
                    dataset.non_molecule = Some(NonMoleculeRecord::Volume {
                        dims: grid.dims(),
                        has_cell: grid.cell().is_some(),
                        has_atoms: grid.atoms().is_some(),
                    });
                    continue;
                }
                Payload::Mesh(mesh) => {
                    dataset.non_molecule = Some(NonMoleculeRecord::Mesh {
                        vertex_count: mesh.num_vertices(),
                        face_count: mesh.num_faces(),
                    });
                    continue;
                }
                Payload::Table(table) => {
                    dataset.non_molecule = Some(NonMoleculeRecord::Table {
                        row_count: table.num_rows(),
                        columns: table.columns().iter().map(|c| c.name.clone()).collect(),
                    });
                    continue;
                }
                // `Payload` is `#[non_exhaustive]` from outside `chem` --
                // a future fifth kind lands here until this match learns it,
                // dropped the same honest way an unparseable record already is.
                _ => {
                    log::warn!(
                        "Dropped a record of an unrecognised kind from a {} dataset",
                        format.label()
                    );
                    continue;
                }
            };
            // Most formats carry coordinates and connectivity rather than a
            // SMILES string, so the column needs something to say. Both the
            // generated string and the placeholder are display decisions, made
            // here rather than in the library that read the file.
            //
            // Where the format states a bond model, writing one is a structure
            // instead of a label (#283) -- and for Mol2 that is only true since
            // #282, which is why this was not done with #266. Where it does
            // not, the placeholder names the format it came from; that used to
            // be the literal `(SDF)`, true while SDF was the only structure
            // format registered and a lie for the nine v0.8.0 added -- open a
            // PDB and every row claimed to be an SDF record (#266).
            let write_one =
                states_a_bond_model(format) && molecule.num_atoms() <= GENERATE_UP_TO_ATOMS;
            let written = smiles.clone().or_else(|| {
                // An atomless molecule writes an empty string, and four readers
                // accept any text as exactly that (#268), so an empty result
                // falls through to the placeholder -- an empty cell reads as a
                // rendering fault.
                write_one
                    .then(|| write_smiles_for_molecule_canonical(&molecule))
                    .filter(|smiles| !smiles.is_empty())
            });
            dataset
                .generated
                .push(smiles.is_none() && written.is_some());
            dataset
                .smiles
                .push(written.unwrap_or_else(|| format!("({})", format.label())));
            dataset.names.push(name);
            dataset.molecules.push(molecule);
        }
        dataset
    }

    /// What this dataset actually holds, in a form the Files list and the
    /// load status line both quote directly (#342).
    ///
    /// `"0 molecules"` for a 10,000-frame trajectory used to be this app's
    /// only answer, because nothing checked `non_molecule` first -- reported
    /// truthfully, indistinguishable from opening an actually-empty file.
    pub fn describe(&self) -> String {
        match &self.non_molecule {
            None => format!(
                "{} {}",
                self.molecules.len(),
                plural(self.molecules.len(), "molecule")
            ),
            Some(NonMoleculeRecord::Trajectory {
                atom_count,
                frame_count,
                has_cell,
                ..
            }) => format!(
                "1 trajectory ({atom_count} atoms, {frame_count} frames{})",
                if *has_cell { ", has a cell" } else { "" }
            ),
            Some(NonMoleculeRecord::Volume {
                dims,
                has_cell,
                has_atoms,
            }) => format!(
                "a {}x{}x{} grid{}{}",
                dims[0],
                dims[1],
                dims[2],
                if *has_cell { ", has a cell" } else { "" },
                if *has_atoms { ", has atoms" } else { "" }
            ),
            Some(NonMoleculeRecord::Mesh {
                vertex_count,
                face_count,
            }) => format!("a mesh ({vertex_count} vertices, {face_count} faces)"),
            Some(NonMoleculeRecord::Table { row_count, columns }) => {
                format!("a table ({row_count} rows, {} columns)", columns.len())
            }
        }
    }

    /// Seeks the active trajectory to `index`, refreshing its cached
    /// per-frame summary. A no-op for every other kind (#342).
    ///
    /// `Trajectory::frame` needs `&mut self`, which a render pass reading
    /// `non_molecule` immutably to draw the slider cannot also hold -- the
    /// view calls this afterwards, once it knows the slider actually moved.
    pub fn seek_trajectory_frame(&mut self, index: usize) {
        let Some(NonMoleculeRecord::Trajectory {
            trajectory,
            selected_frame,
            current_frame_summary,
            ..
        }) = &mut self.non_molecule
        else {
            return;
        };
        *selected_frame = index;
        *current_frame_summary = match trajectory.frame(index) {
            Ok(frame) => describe_frame(&frame),
            Err(e) => format!("could not read frame {index}: {e}"),
        };
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
        Self::from_outcome(outcome, DatasetFormat::SMILES)
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

    /// The format the active dataset was read from.
    ///
    /// What a conversion's losses depend on as much as the target does: a
    /// read and a write each lose their own things, and some pairs lose
    /// something neither mask predicts (#275).
    pub fn active_format(&self) -> DatasetFormat {
        self.entries()[self.active_index()].format
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
        MoleculeDataset::from_outcome(outcome, DatasetFormat::SMILES)
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

    /// Writes a molecule in `format` and reads it back the way a loaded file
    /// is read, so the column under test is the one a user would see.
    ///
    /// Whether the *reader* stated a SMILES is the only honest way to check
    /// which strings the app wrote itself -- read from the record before
    /// `from_outcome` takes `outcome` by value, since it moves a non-molecule
    /// payload's contents out (#342) and `ReadOutcome`/`Record` are not
    /// `Clone`.
    fn round_trip(format: DatasetFormat, name: &str, smiles: &str) -> (bool, MoleculeDataset) {
        let molecule = chem::io::smiles::parse_smiles(smiles).expect("valid SMILES");
        let text = format
            .write(&[(name.to_string(), molecule)])
            .unwrap_or_else(|| panic!("{} writes", format.label()));
        let outcome = reader::read(&text, format);
        let stated = outcome.records[0].smiles.is_some();
        let dataset = MoleculeDataset::from_outcome(outcome, format);
        (stated, dataset)
    }

    #[test]
    fn test_every_format_fills_the_column_and_names_itself_only_when_it_must() {
        // The general guard, in the spirit of #257's mask assertions: one loop
        // over the registry rather than a test per format, so a new format
        // is covered the day it registers.
        //
        // BinaryCIF (#319) is excluded: its canonical bytes are not text, so
        // `DatasetFormat::write` correctly answers `None` for it (same as
        // `chem::io::format::Format::write`) rather than the `round_trip`
        // helper's text-only path applying.
        //
        // `Kind::Frames` formats (TRR/XTC/DCD/NCTRAJ/LAMMPS Trajectory,
        // #325-#329) are excluded too, for the same reason regardless of
        // encoding: their writer lives in `writer_trajectory`
        // (`&mut Trajectory -> Vec<u8>`), not `writer`/`writer_bytes`
        // (`&[(String, Molecule)] -> ...`), so `Format::write` -- and this
        // test's `round_trip` helper, which only ever builds a `Molecule`
        // -- correctly has nothing to call for them. LAMMPS Trajectory is
        // the first of these that's also `Encoding::Text`, so the encoding
        // filter alone stopped being enough to exclude every non-Molecules
        // format once it registered. The app now shows one of these as a
        // `NonMoleculeRecord` summary instead of a molecule row (#342) --
        // irrelevant here regardless, since this loop only exercises the
        // `Molecule`-shaped writer/reader path `round_trip` calls.
        for format in chem::io::format::all()
            .filter(|f| f.can_read() && f.can_write())
            .filter(|f| f.encoding() == chem::io::format::Encoding::Text)
            .filter(|f| f.kind() == chem::io::format::Kind::Molecules)
        {
            let (stated, dataset) = round_trip(format, "benzene", "c1ccccc1");
            let cell = &dataset.smiles[0];
            let placeholder = *cell == format!("({})", format.label());

            assert!(
                !cell.is_empty(),
                "{} leaves the column blank",
                format.label()
            );
            assert_eq!(
                placeholder,
                !states_a_bond_model(format),
                "{} shows {cell:?}",
                format.label()
            );
            assert_eq!(
                dataset.generated[0],
                !stated && !placeholder,
                "{} marks {cell:?} wrongly",
                format.label()
            );
        }
    }

    #[test]
    fn test_a_bond_model_is_written_as_the_molecule_not_the_format_name() {
        // The reported symptom (#283). Mol2 is only correct here since #282 --
        // before it, the same round trip produced `[c]1[c][c][c][c][c]1` (#281)
        // and this column would have displayed something wrong.
        for format in [DatasetFormat::MOL2, DatasetFormat::CML, DatasetFormat::SDF] {
            let (_, dataset) = round_trip(format, "benzene", "c1ccccc1");
            assert_eq!(dataset.smiles[0], "c1ccccc1", "{}", format.label());
            assert!(dataset.generated[0], "{}", format.label());
        }
    }

    #[test]
    fn test_a_format_stating_no_bond_model_keeps_its_name_even_holding_bonds() {
        // PDB's `CONECT` is adjacency without bond order, so a SMILES written
        // from it is the right topology and the wrong molecule. This is the
        // case that makes the predicate a format property: the molecule here
        // *has* every bond, and `num_bonds() > 0` -- what #283 proposed --
        // would have written `[C]1[C][C][C][C][C]1` into a column headed
        // SMILES.
        let (_, dataset) = round_trip(DatasetFormat::PDB, "benzene", "c1ccccc1");
        assert!(
            dataset.molecules[0].num_bonds() > 0,
            "the fixture must hold bonds or it passes for the wrong reason"
        );
        assert_eq!(dataset.smiles[0], "(PDB)");
        assert!(!dataset.generated[0]);
    }

    #[test]
    fn test_a_molecule_with_no_bonds_is_still_written_when_its_format_stated_so() {
        // The other direction of the same point. A salt read from an SDF has
        // zero bonds and the file said so, which is information -- unlike zero
        // bonds from XYZ, which is the format's silence. `num_bonds() > 0`
        // would suppress a correct answer.
        let (_, dataset) = round_trip(DatasetFormat::SDF, "salt", "[Na+].[Cl-]");
        assert_eq!(dataset.molecules[0].num_bonds(), 0);
        assert_eq!(dataset.smiles[0], "[Na+].[Cl-]");
        assert!(dataset.generated[0]);
    }

    #[test]
    fn test_a_record_with_no_atoms_is_skipped_rather_than_shown_as_an_empty_row() {
        // Before #268, the four lenient structure readers accepted arbitrary
        // text as one atomless molecule, and the canonical writer's empty
        // string for it fell back to a `(Mol2)`-style placeholder cell --
        // indistinguishable from a real, if unusual, zero-atom read. Now the
        // reader reports it as skipped instead, so `from_outcome` never sees
        // a record to build a row from at all.
        let outcome = reader::read("not a molecule\n", DatasetFormat::MOL2);
        assert_eq!(outcome.skipped.len(), 1);

        let dataset = MoleculeDataset::from_outcome(outcome, DatasetFormat::MOL2);
        assert_eq!(dataset.len(), 0);
    }

    #[test]
    fn test_a_molecule_too_large_to_write_keeps_the_placeholder() {
        // The ceiling is a cost decision (canonical ranking is quadratic, and
        // loading blocks the only thread the web build has) that also happens
        // to be a legibility one. A ring one atom over it is the cheapest
        // fixture that crosses the line.
        let ring = format!("C1{}1", "C".repeat(GENERATE_UP_TO_ATOMS));
        let (_, dataset) = round_trip(DatasetFormat::MOL2, "big", &ring);

        assert!(dataset.molecules[0].num_atoms() > GENERATE_UP_TO_ATOMS);
        assert_eq!(dataset.smiles[0], "(Mol2)");
        assert!(!dataset.generated[0]);
    }
}
