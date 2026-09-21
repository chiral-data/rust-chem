//! The loaded datasets, and the active one's molecules.
//!
//! Only what is about data. The fingerprint controls, aromaticity and the
//! structure display options all used to be wedged in between the file list and
//! the table; #105 and #106 took them to the windows that own them.
//!
//! The table is virtualised: it draws the rows on screen and no others, so a
//! dataset of ten thousand molecules scrolls rather than being silently cut to
//! the first twenty.

use crate::dataset::NonMoleculeRecord;
use crate::molecule_view::GENERATED_NOTE;
use crate::state::AppState;
use crate::structure_view::StructureView;
use chem::draw::structure::{ShowCarbons, StructureOptions};
use egui::{RichText, Vec2};
use egui_extras::{Column, TableBuilder};

/// Height of a table row. Fixed rather than measured, because virtualisation
/// needs to know which rows fall in the viewport before drawing any of them.
const ROW_HEIGHT: f32 = 22.0;

/// Row height when thumbnails are on, giving each structure a legible cell.
const ROW_HEIGHT_WITH_THUMBNAILS: f32 = 48.0;

#[derive(Default)]
pub struct DatasetsView;

impl DatasetsView {
    pub fn ui(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        // No outer scroll area: the table owns the scrolling, and nesting one
        // inside another makes the wheel ambiguous. Collapse Files if the
        // window is too short for both.
        self.files_section(ui, state);
        ui.separator();
        self.table(ui, state);
    }

    fn files_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let active = state.loaded_files.active_index();
        let summary = format!(
            "{} loaded \u{b7} {} active",
            state.loaded_files.entries().len(),
            state.loaded_files.active_dataset().describe()
        );

        egui::CollapsingHeader::new(RichText::new("Files").strong())
            .default_open(true)
            .show(ui, |ui| {
                ui.label(RichText::new(&summary).small().weak());

                ui.horizontal(|ui| {
                    if ui.button("📂 Load File").clicked() {
                        state.load_dataset_from_file();
                    }
                    if ui.button("📋 Load Examples").clicked() {
                        state.load_example_dataset();
                    }
                });

                ui.label(&state.dataset_status);

                // Name, format and what it holds per entry. `LoadedFiles`
                // always knew the latter two and the list only ever showed
                // the name, so telling two SMILES files apart meant
                // switching to them. `describe()` (#342) is what keeps a
                // trajectory or a grid from reading as "0 molecules" here.
                let described: Vec<(String, &'static str, String)> = state
                    .loaded_files
                    .entries()
                    .iter()
                    .map(|e| (e.name.clone(), e.format.label(), e.dataset.describe()))
                    .collect();

                let can_remove = state.loaded_files.can_remove();
                let mut clicked = None;
                let mut removed = None;
                for (i, (name, format, description)) in described.iter().enumerate() {
                    ui.horizontal(|ui| {
                        // Disabled at one entry rather than hidden, so the
                        // control doesn't appear and vanish as files are loaded.
                        let remove = ui
                            .add_enabled_ui(can_remove, remove_button)
                            .inner
                            .on_hover_text("Remove this dataset")
                            .on_disabled_hover_text("The last dataset can't be removed");
                        if remove.clicked() {
                            removed = Some(i);
                        }

                        let label = format!("{}  ({}, {})", name, format, description);
                        if ui.selectable_label(i == active, label).clicked() {
                            clicked = Some(i);
                        }
                    });
                }

                // Removal first: acting on a click into a list that has just
                // changed length would activate the wrong entry.
                if let Some(i) = removed {
                    state.remove_loaded_file(i);
                } else if let Some(i) = clicked {
                    state.activate_loaded_file(i);
                }
            });
    }

    fn table(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let dataset = state.loaded_files.active_dataset();
        if dataset.is_empty() {
            // A non-`Kind::Molecules` file also has zero `molecules` -- but
            // not for the reason an empty file does (#342). Checked first,
            // so a 10,000-frame trajectory never again reads as "No
            // molecules", indistinguishable from opening a blank file.
            if dataset.non_molecule.is_some() {
                self.non_molecule_summary(ui, state);
            } else {
                ui.label(RichText::new("No molecules — load a file or the examples.").weak());
            }
            return;
        }

        let thumbnails = state.display.show_thumbnails;
        let fingerprinted = state.dataset_fingerprints.len();
        let total = dataset.len();
        let row_height = if thumbnails {
            ROW_HEIGHT_WITH_THUMBNAILS
        } else {
            ROW_HEIGHT
        };

        // A thumbnail is read as a shape rather than for its labels: hydrogens
        // are dropped and the structure fits its cell rather than using a
        // shared bond length. Carbons stay on Default rather than None, since
        // None would leave a lone carbon with neither a label nor a bond and
        // methane's cell would draw empty.
        let thumbnail_options = StructureOptions {
            padding: 2.0,
            show_carbons: ShowCarbons::Default,
            explicit_hydrogens: false,
            scale: 0.0,
            bond_length: 0.0,
            ..state.display.structure
        };

        ui.label(RichText::new(format!("{} molecules", total)).small().weak());

        let mut clicked_row = None;

        let mut table = TableBuilder::new(ui)
            .striped(true)
            .resizable(true)
            .vscroll(true)
            .auto_shrink([false, false])
            .cell_layout(egui::Layout::left_to_right(egui::Align::Center));

        if thumbnails {
            table = table.column(Column::exact(56.0));
        }
        table = table
            .column(Column::initial(120.0).at_least(60.0)) // Name
            .column(Column::initial(180.0).at_least(80.0)) // SMILES
            .column(Column::initial(90.0).at_least(50.0)) // Formula
            .column(Column::initial(70.0).at_least(45.0)) // MW
            .column(Column::initial(80.0).at_least(45.0)) // Fingerprint
            .column(Column::remainder().at_least(60.0)); // Aromatic

        table
            .header(20.0, |mut header| {
                if thumbnails {
                    header.col(|ui| {
                        ui.label(RichText::new("Structure").strong());
                    });
                }
                for title in ["Name", "SMILES", "Formula", "MW", "Fingerprint", "Aromatic"] {
                    header.col(|ui| {
                        ui.label(RichText::new(title).strong());
                    });
                }
            })
            .body(|body| {
                // Only the rows in the viewport are built, so the cost of a row
                // is paid for the screenful being looked at rather than for the
                // dataset.
                body.rows(row_height, total, |mut row| {
                    let i = row.index();
                    let mol = &dataset.molecules[i];

                    if thumbnails {
                        row.col(|ui| {
                            // Laid out on demand and shared with the detail
                            // windows and the result rows, so all three draw
                            // the same picture. This used to be a dash and an
                            // instruction to run 2D Coordinates, while the
                            // detail window drew the same molecule fine (#273).
                            if let Some(drawable) = state.drawable(i) {
                                ui.add(
                                    StructureView::new(&drawable, Vec2::new(52.0, 44.0))
                                        .with_options(thumbnail_options),
                                );
                            }
                        });
                    }

                    row.col(|ui| {
                        // Lit while this molecule's detail window is open — with
                        // several open at once, the highlight is what says which.
                        if ui
                            .selectable_label(state.is_detail_open(i), &dataset.names[i])
                            .clicked()
                        {
                            clicked_row = Some(i);
                        }
                    });
                    row.col(|ui| {
                        let cell = ui.label(RichText::new(&dataset.smiles[i]).code().small());
                        // The column is 180px and a label clips rather than
                        // wraps, so hovering is also how a long string is read
                        // in full -- but the reason it is here is that this
                        // string may be one the app wrote rather than one the
                        // file stated, and the two are otherwise identical
                        // (#283).
                        if dataset.generated[i] {
                            cell.on_hover_text(GENERATED_NOTE);
                        }
                    });
                    row.col(|ui| {
                        ui.label(mol.formula());
                    });
                    row.col(|ui| {
                        ui.label(format!("{:.2}", mol.molecular_weight()));
                    });
                    row.col(|ui| {
                        ui.label(if i < fingerprinted { "Yes" } else { "No" });
                    });
                    row.col(|ui| {
                        let aromatic = mol.atoms().iter().any(|atom| atom.is_aromatic());
                        ui.label(if aromatic { "Yes" } else { "No" });
                    });
                });
            });

        if let Some(i) = clicked_row {
            state.toggle_detail(i);
        }
    }

    /// The active dataset's `NonMoleculeRecord`, as text -- never a picture,
    /// per #307's own "readable, not drawn" rule for this whole wave (#342).
    ///
    /// Reads what it needs into owned locals before drawing, ending the
    /// borrow of `state.loaded_files` before the slider's seek (which needs
    /// `&mut`) is applied -- an immutable read and a `Trajectory::frame`
    /// seek cannot both be live at once.
    fn non_molecule_summary(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let dataset = state.loaded_files.active_dataset();
        let Some(non_molecule) = dataset.non_molecule.as_ref() else {
            return;
        };

        match non_molecule {
            NonMoleculeRecord::Trajectory {
                atom_count,
                frame_count,
                has_cell,
                selected_frame,
                current_frame_summary,
                ..
            } => {
                let (atom_count, frame_count, has_cell) = (*atom_count, *frame_count, *has_cell);
                let mut frame = *selected_frame;
                let current_frame_summary = current_frame_summary.clone();

                ui.label(format!(
                    "Trajectory: {atom_count} atoms, {frame_count} frames{}",
                    if has_cell { ", has a cell" } else { "" }
                ));
                let response = ui.add(
                    egui::Slider::new(&mut frame, 0..=frame_count.saturating_sub(1)).text("frame"),
                );
                ui.label(RichText::new(&current_frame_summary).small().weak());

                if response.changed() {
                    state
                        .loaded_files
                        .active_dataset_mut()
                        .seek_trajectory_frame(frame);
                }
            }
            NonMoleculeRecord::Volume {
                dims,
                has_cell,
                has_atoms,
            } => {
                ui.label(format!(
                    "Volume grid: {}x{}x{}{}{}",
                    dims[0],
                    dims[1],
                    dims[2],
                    if *has_cell { ", has a cell" } else { "" },
                    if *has_atoms { ", has atoms" } else { "" }
                ));
            }
            NonMoleculeRecord::Mesh {
                vertex_count,
                face_count,
            } => {
                ui.label(format!("Mesh: {vertex_count} vertices, {face_count} faces"));
            }
            NonMoleculeRecord::Table { row_count, columns } => {
                ui.label(format!(
                    "Table: {row_count} rows \u{b7} columns: {}",
                    columns.join(", ")
                ));
            }
        }
    }
}

/// A close control, painted rather than lettered.
///
/// egui paints its own window close button from two line segments instead of
/// using a character, and the reason shows up the moment you try: the
/// multiplication-x glyphs aren't in the bundled font subset, so `✕` renders as
/// a missing-glyph box. Painting it takes the stroke from the interaction
/// visuals, so it also follows the theme and greys itself when disabled without
/// being told to.
fn remove_button(ui: &mut egui::Ui) -> egui::Response {
    let size = egui::Vec2::splat(ui.spacing().icon_width);
    let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click());

    if ui.is_rect_visible(rect) {
        let visuals = ui.style().interact(&response);
        let stroke = visuals.fg_stroke;
        let cross = rect.shrink(rect.width() * 0.28);
        let painter = ui.painter();
        painter.line_segment([cross.left_top(), cross.right_bottom()], stroke);
        painter.line_segment([cross.right_top(), cross.left_bottom()], stroke);
    }

    response
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dataset::DatasetFormat;
    use crate::state::AppState;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;

    /// A harness over the Files window, with whatever `prepare` puts in it.
    ///
    /// Drives the view rather than the whole app: `WorkbenchApp::new` wants an
    /// `eframe::CreationContext` and probes for a GPU, while a view wants a
    /// `Ui` and nothing else -- so this needs neither a window nor an adapter,
    /// which is what lets it run in CI (#303).
    ///
    /// State is mutated through `AppState`'s own methods rather than by
    /// simulating the input that would call them. The plumbing already has
    /// tests in `state.rs`; what has none is whether the values reach the
    /// screen, and that is all these assert. (`poll_pending_work` reads the
    /// context `AppState` stored at construction, which is detached here, so
    /// input injected into the harness would not reach it anyway.)
    fn files_window(prepare: impl FnOnce(&mut AppState)) -> Harness<'static, AppState> {
        let mut state = AppState::cpu_only();
        prepare(&mut state);
        let mut view = DatasetsView;
        let mut harness = Harness::new_ui_state(move |ui, state| view.ui(ui, state), state);
        // Not load-bearing today -- measured, the default already renders all
        // 15 examples including the last -- but pinned so that a change to the
        // harness default cannot clip a row and turn a real assertion into a
        // spurious failure. `body.rows` builds only what is visible.
        harness.set_size(egui::vec2(1400.0, 1200.0));
        harness.run();
        harness
    }

    /// Writes a molecule in `format` and loads it the way a file load does.
    fn load(state: &mut AppState, name: &str, format: DatasetFormat, smiles: &str) {
        let molecule = chem::io::smiles::parse_smiles(smiles).expect("valid SMILES");
        let text = format
            .write(&[(name.to_string(), molecule)])
            .unwrap_or_else(|| panic!("{} writes", format.label()));
        state.apply_loaded_file_bytes(name.to_string(), text.into_bytes());
    }

    #[test]
    fn test_the_examples_reach_the_table() {
        // The floor: if this fails, nothing below it means anything. Fifteen
        // built-in molecules, and the assertion is on the *last* of them --
        // the row a clipped viewport would lose.
        let harness = files_window(|_| {});

        assert!(harness.query_by_label("Methane").is_some());
        assert!(
            harness.query_by_label("Neopentane").is_some(),
            "the last example row did not render"
        );
        for header in ["Name", "SMILES", "Formula", "MW"] {
            assert!(
                harness.query_by_label(header).is_some(),
                "the {header} header did not render"
            );
        }
    }

    #[test]
    fn test_a_dropped_file_becomes_rows_on_screen() {
        // `apply_dropped_files` has tests; that the table then shows the
        // molecules does not (#296).
        let harness = files_window(|state| {
            state.apply_dropped_files(vec![(
                "dropped.smi".to_string(),
                b"CCO ethanol\nCC(=O)O acetic\n".to_vec(),
            )]);
        });

        assert!(harness.query_by_label("ethanol").is_some());
        assert!(harness.query_by_label("acetic").is_some());
        // The list entry is `name  (FORMAT, n molecules)`, so this pins the
        // format and the count #104 added alongside the name.
        assert!(
            harness
                .query_by_label_contains("dropped.smi  (SMILES, 2 molecules)")
                .is_some(),
            "the Files list did not describe the dropped file"
        );
    }

    #[test]
    fn test_the_smiles_column_shows_the_molecule_or_the_format() {
        // #283's rule, on screen: a format stating a bond model gets a written
        // SMILES, one that does not keeps its own name. Both strings are
        // pinned in `dataset.rs`; that the column displays them is not.
        let harness = files_window(|state| {
            load(state, "rings.mol2", DatasetFormat::MOL2, "c1ccccc1");
        });
        assert!(
            harness.query_by_label("c1ccccc1").is_some(),
            "a Mol2 should show the molecule it holds"
        );

        let harness = files_window(|state| {
            load(state, "ethanol.pdb", DatasetFormat::PDB, "CCO");
        });
        assert!(
            harness.query_by_label("(PDB)").is_some(),
            "a PDB has no aromatic model, so the column keeps the format name"
        );
    }

    #[test]
    fn test_the_formula_column_counts_the_hydrogens_a_pdb_implies() {
        // The value #294 moved. Before it, a PDB-read ethanol had no hydrogen
        // count and this cell read C2O; the fix is only visible to a user
        // here and in MW.
        let harness = files_window(|state| {
            load(state, "ethanol.pdb", DatasetFormat::PDB, "CCO");
        });

        assert!(
            harness.query_by_label("C2H6O").is_some(),
            "the Formula column should count the hydrogens CONECT implies"
        );
        assert!(
            harness.query_by_label("C2O").is_none(),
            "C2O is the pre-#294 answer and must not come back"
        );
    }

    #[test]
    fn test_a_refused_drop_names_the_file_in_the_status_line() {
        // The status line is the only place this can appear: a refused file
        // leaves no Files entry to notice afterwards (#296).
        let harness = files_window(|state| {
            state.apply_dropped_files(vec![
                ("good.smi".to_string(), b"CCO ethanol\n".to_vec()),
                ("logo.png".to_string(), vec![0x89, b'P', b'N', b'G']),
            ]);
        });

        assert!(
            harness
                .query_by_label_contains("logo.png: not a format this build reads")
                .is_some(),
            "the refusal did not reach the status line"
        );
        assert!(
            harness.query_by_label("ethanol").is_some(),
            "the file beside it should still have loaded"
        );
    }

    /// A two-frame LAMMPS dump, plain text and dependency-free -- the
    /// easiest of the twelve non-`Kind::Molecules` formats to hand-author.
    const TWO_FRAME_TRAJECTORY: &str = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
ITEM: TIMESTEP
1
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
1 1 5.0 0.0 0.0
";

    #[test]
    fn test_a_trajectory_reports_frames_not_zero_molecules() {
        // Before #342, `MoleculeDataset::from_outcome` dropped this record
        // silently and the app reported "Loaded 0 molecules from
        // 'traj.lammpstrj' (LAMMPS Trajectory)" -- indistinguishable from
        // opening an actually-empty file.
        let harness = files_window(|state| {
            state.apply_loaded_file_bytes(
                "traj.lammpstrj".to_string(),
                TWO_FRAME_TRAJECTORY.as_bytes().to_vec(),
            );
        });

        assert!(
            harness
                .query_all_by_label_contains("1 trajectory")
                .next()
                .is_some(),
            "a loaded trajectory must say what it is"
        );
        assert!(
            harness
                .query_all_by_label_contains("2 frames")
                .next()
                .is_some(),
            "the frame count must reach the screen"
        );
        assert!(
            harness
                .query_all_by_label_contains("0 molecules")
                .next()
                .is_none(),
            "a real trajectory must never read like an empty file"
        );
    }

    #[test]
    fn test_the_frame_slider_shows_a_different_frame_on_seek() {
        // The slider itself can't be dragged through this harness
        // (`files_window`'s own doc: injected input never reaches a context
        // detached at construction) -- so this drives the same method the
        // slider's `changed()` branch calls, and asserts the effect reaches
        // the screen, the same "state through AppState's own methods"
        // convention every other test here already uses.
        let mut state = AppState::cpu_only();
        state.apply_loaded_file_bytes(
            "two.lammpstrj".to_string(),
            TWO_FRAME_TRAJECTORY.as_bytes().to_vec(),
        );

        let mut view = DatasetsView;
        let mut harness = Harness::new_ui_state(move |ui, state| view.ui(ui, state), state);
        harness.set_size(egui::vec2(1400.0, 1200.0));
        harness.run();
        assert!(
            harness
                .query_by_label_contains("(0.00, 0.00, 0.00) to (0.00, 0.00, 0.00)")
                .is_some(),
            "frame 0 should show the single atom at the origin"
        );

        harness
            .state_mut()
            .loaded_files
            .active_dataset_mut()
            .seek_trajectory_frame(1);
        harness.run();
        assert!(
            harness
                .query_by_label_contains("(5.00, 0.00, 0.00) to (5.00, 0.00, 0.00)")
                .is_some(),
            "seeking to frame 1 should show that frame's own atom position, not frame 0's"
        );
    }

    #[test]
    fn test_a_volume_reports_its_grid_not_zero_molecules() {
        let cube = "\
Test CUBE
comment
    1    0.000000    0.000000    0.000000
    1    1.000000    0.000000    0.000000
    1    0.000000    1.000000    0.000000
    1    0.000000    0.000000    1.000000
    6    6.000000    0.500000    0.500000    0.500000
0.0
";
        let harness = files_window(|state| {
            state.apply_loaded_file_bytes("density.cube".to_string(), cube.as_bytes().to_vec());
        });

        assert!(
            harness
                .query_all_by_label_contains("1x1x1 grid")
                .next()
                .is_some(),
            "a loaded volume must describe its own grid"
        );
        assert!(
            harness
                .query_all_by_label_contains("0 molecules")
                .next()
                .is_none()
        );
    }

    #[test]
    fn test_a_mesh_reports_vertices_and_faces_not_zero_molecules() {
        let obj = "\
v 0 0 0
v 1 0 0
v 0 1 0
f 1 2 3
";
        let harness = files_window(|state| {
            state.apply_loaded_file_bytes("shape.obj".to_string(), obj.as_bytes().to_vec());
        });

        assert!(
            harness
                .query_all_by_label_contains("3 vertices")
                .next()
                .is_some(),
            "a loaded mesh must describe its own geometry"
        );
        assert!(
            harness
                .query_all_by_label_contains("0 molecules")
                .next()
                .is_none()
        );
    }

    #[test]
    fn test_a_table_reports_rows_and_columns_not_zero_molecules() {
        // No structure column, so #337's own rule reads this as a table
        // rather than a molecule set.
        let csv = "a,b\n1,2\n3,4\n";
        let harness = files_window(|state| {
            state.apply_loaded_file_bytes("data.csv".to_string(), csv.as_bytes().to_vec());
        });

        assert!(
            harness
                .query_all_by_label_contains("2 rows")
                .next()
                .is_some(),
            "a loaded table must describe its own shape"
        );
        assert!(
            harness
                .query_all_by_label_contains("0 molecules")
                .next()
                .is_none()
        );
    }
}
