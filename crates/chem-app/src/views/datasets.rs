//! The loaded datasets, and the active one's molecules.
//!
//! Only what is about data. The fingerprint controls, aromaticity and the
//! structure display options all used to be wedged in between the file list and
//! the table; #105 and #106 took them to the windows that own them.
//!
//! The table is virtualised: it draws the rows on screen and no others, so a
//! dataset of ten thousand molecules scrolls rather than being silently cut to
//! the first twenty.

use crate::state::AppState;
use crate::structure_view::StructureView;
use chem_draw::structure::{ShowCarbons, StructureOptions};
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
            "{} loaded \u{b7} {} molecules active",
            state.loaded_files.entries().len(),
            state.loaded_files.active_dataset().len()
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

                // Name, format and size per entry. `LoadedFiles` always knew
                // the latter two and the list only ever showed the name, so
                // telling two SMILES files apart meant switching to them.
                let described: Vec<(String, &'static str, usize)> = state
                    .loaded_files
                    .entries()
                    .iter()
                    .map(|e| (e.name.clone(), e.format.label(), e.dataset.len()))
                    .collect();

                let can_remove = state.loaded_files.can_remove();
                let mut clicked = None;
                let mut removed = None;
                for (i, (name, format, len)) in described.iter().enumerate() {
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

                        let label = format!("{}  ({}, {} molecules)", name, format, len);
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
            ui.label(RichText::new("No molecules — load a file or the examples.").weak());
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
                            // Drawn only where coordinates already exist —
                            // generating them needs the dataset mutably, which
                            // the table cannot have while reading it. Run 2D
                            // Coordinates in Operations to fill them in; SDF
                            // files arrive with their own.
                            if mol.has_coords() {
                                ui.add(
                                    StructureView::new(mol, Vec2::new(52.0, 44.0))
                                        .with_options(thumbnail_options),
                                );
                            } else {
                                ui.label(RichText::new("\u{2014}").weak());
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
                        ui.label(RichText::new(&dataset.smiles[i]).code().small());
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
