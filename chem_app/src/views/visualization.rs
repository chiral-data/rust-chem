//! Ranked search results, and the per-hit detail an expanded row shows.
//!
//! #106 turns this into the window that owns the display settings and every
//! structure view alongside the results; for now it is the results panel, moved
//! off `app.rs` with its own expansion state.

use crate::fingerprint_view::fingerprint_compact;
use crate::molecule_view::{molecule_compact, show_atom_list, show_bond_list, show_molecule_info};
use crate::state::AppState;
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

        ui.heading("Search Results");
        ui.separator();

        if state.search_results.is_empty() {
            ui.centered_and_justified(|ui| {
                ui.label(RichText::new("No search results yet").size(20.0).weak());
            });
            return;
        }

        let active_dataset = state.loaded_files.active_dataset();
        let mut toggled = None;

        egui::ScrollArea::vertical().show(ui, |ui| {
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
        });

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
