//! Where you look at your data.
//!
//! Owns the display options *and* the things they govern, which is the point:
//! before this the options that decide how every structure in the app is drawn
//! lived in the data panel, under the dataset table, while the structures
//! themselves were in two other panels and a floating window.
//!
//! Three sections — how to draw, the query, and the ranking a search produced.
//! The per-molecule structures stay in detail windows, which #107 turns into one
//! window per molecule; drawing them here as well would only duplicate that.

use crate::fingerprint_view::{fingerprint_compact, fingerprint_full};
use crate::molecule_view::{molecule_compact, show_atom_list, show_bond_list, show_molecule_info};
use crate::state::AppState;
use crate::structure_view::{structure_option_controls, structure_panel_with_options};
use egui::{Color32, RichText};

#[derive(Default)]
pub struct VisualizationView {
    /// Which result row is expanded, and the results epoch it refers to.
    ///
    /// Indexes into `AppState::search_results`, which is replaced wholesale
    /// every time a search finishes — so the epoch is what tells us the index
    /// we're holding belongs to results that no longer exist.
    selected_result: Option<usize>,
    results_epoch: u64,
}

impl VisualizationView {
    pub fn ui(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        self.sync(state);

        // Both axes, no auto-shrink: a fingerprint grid or a wide result row
        // should scroll inside the window rather than stretch it over its
        // neighbours.
        egui::ScrollArea::both()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                display_section(ui, state);
                query_section(ui, state);
                self.results_section(ui, state);
            });
    }

    fn results_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let count = state.search_results.len();
        let summary = if count == 0 {
            "\u{2014}".to_string()
        } else {
            format!("{} hits", count)
        };

        collapsing(ui, "Results", true, &summary, |ui| {
            if state.search_results.is_empty() {
                ui.label(
                    RichText::new("Nothing searched yet — run Similarity Search in Operations.")
                        .weak(),
                );
                return;
            }
            self.result_rows(ui, state);
        });
    }

    fn result_rows(&mut self, ui: &mut egui::Ui, state: &AppState) {
        let active_dataset = state.loaded_files.active_dataset();
        let mut toggled = None;

        for (rank, result) in state.search_results.iter().enumerate() {
            let idx = result.index;
            let Some(mol) = active_dataset.molecules.get(idx) else {
                continue;
            };
            let smiles = &active_dataset.smiles[idx];
            let name = &active_dataset.names[idx];

            let is_selected = self.selected_result == Some(rank);

            let response = ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("#{}", rank + 1)).strong().size(16.0));
                    ui.separator();

                    ui.vertical(|ui| {
                        molecule_compact(ui, smiles, name, &mol.formula());
                        ui.label(format!(
                            "Similarity: {:.3} ({:.1}%)",
                            result.similarity,
                            result.similarity * 100.0
                        ));
                    });

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .button(if is_selected { "▲ Hide" } else { "▼ Show" })
                            .clicked()
                        {
                            toggled = Some(rank);
                        }
                    });
                });
            });

            if response.response.hovered() {
                ui.painter().rect_stroke(
                    response.response.rect,
                    2.0,
                    egui::Stroke::new(2.0_f32, Color32::from_rgb(70, 130, 180)),
                    egui::epaint::StrokeKind::Outside,
                );
            }

            if is_selected {
                ui.indent(format!("result_{}", rank), |ui| {
                    ui.horizontal_top(|ui| {
                        ui.vertical(|ui| {
                            show_molecule_info(ui, mol, smiles, name);
                            show_atom_list(ui, mol);
                            show_bond_list(ui, mol);
                        });

                        ui.separator();

                        if let Some(fp) = state.dataset_fingerprints.get(idx) {
                            ui.vertical(|ui| {
                                fingerprint_compact(ui, fp, "Molecule Fingerprint");

                                if let Some(query_fp) = &state.query_fingerprint {
                                    ui.separator();
                                    fingerprint_compact(ui, query_fp, "Query Fingerprint");
                                }
                            });
                        }
                    });
                });
            }

            ui.separator();
        }

        if let Some(rank) = toggled {
            self.selected_result = if self.selected_result == Some(rank) {
                None
            } else {
                Some(rank)
            };
        }
    }

    /// Drops an expansion that refers to results which have since been
    /// replaced.
    fn sync(&mut self, state: &AppState) {
        if self.results_epoch != state.results_epoch() {
            self.results_epoch = state.results_epoch();
            self.selected_result = None;
        }
    }
}

/// How structures are drawn, everywhere.
///
/// Closed by default: it is configuration you set once, not something to read.
/// The one display choice *not* here is the dataset table's thumbnail toggle,
/// which stays beside the table it governs.
fn display_section(ui: &mut egui::Ui, state: &mut AppState) {
    collapsing(ui, "Display", false, "all structures", |ui| {
        structure_option_controls(ui, &mut state.display.structure);
    });
}

/// The parsed query: its structure, its details, and its fingerprint.
///
/// Drawn from `AppState::query_source` rather than the query box's live text —
/// that belongs to Operations, and differs from this the moment you type.
fn query_section(ui: &mut egui::Ui, state: &mut AppState) {
    let summary = if state.query_molecule.is_some() {
        state.query_source.clone()
    } else {
        "\u{2014}".to_string()
    };

    collapsing(ui, "Query", true, &summary, |ui| {
        let Some(mol) = state.query_molecule.clone() else {
            ui.label(RichText::new("Nothing parsed yet — enter a SMILES in Operations.").weak());
            return;
        };

        structure_panel_with_options(ui, &mol, 200.0, state.display.structure);
        show_molecule_info(ui, &mol, &state.query_source, "Query");

        if let Some(fp) = &state.query_fingerprint {
            fingerprint_full(ui, fp);
        }
    });
}

/// A section header carrying a one-line summary, so a collapsed section still
/// says what is inside it.
///
/// Same shape as the Operations window's sections, deliberately — but not the
/// same function: those summarise an `OperationOutcome`, these summarise
/// whatever the section happens to hold, and forcing one helper to do both would
/// mean inventing an outcome for things that never ran.
fn collapsing(
    ui: &mut egui::Ui,
    title: &str,
    default_open: bool,
    summary: &str,
    add_contents: impl FnOnce(&mut egui::Ui),
) {
    egui::CollapsingHeader::new(RichText::new(title).strong())
        .default_open(default_open)
        .show(ui, |ui| {
            ui.label(RichText::new(summary).small().weak());
            add_contents(ui);
        });
}
