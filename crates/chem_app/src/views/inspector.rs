//! Where you look at your data.
//!
//! The query and the ranking a search produced.
//!
//! Result rows draw the molecule they found. Similarity is a visual question —
//! two molecules can score 0.9 for reasons plain in a drawing and invisible in a
//! SMILES string, which is how a fingerprinting bug once went unnoticed while
//! phenol scored 0.27 against itself.
//!
//! Held the display options too, briefly. They govern three surfaces and only
//! one of them is this window, so #121 moved them to Settings; what is left here
//! is the data itself. The per-molecule structures live in detail windows, one
//! per molecule, rather than being duplicated here.

use crate::fingerprint_view::{fingerprint_comparison, fingerprint_full};
use crate::molecule_view::{molecule_compact, show_molecule_info};
use crate::state::AppState;
use crate::structure_view::{StructureView, structure_panel_with_options};
use crate::svg::save_svg;
use chemdraw::structure::{ShowCarbons, StructureOptions, StructureTheme};
use chemdraw::svg::{structure_to_svg, suggested_filename};
use egui::{Color32, RichText, Vec2};

#[derive(Default)]
pub struct InspectorView {
    /// Which result row is expanded, and the results epoch it refers to.
    ///
    /// Indexes into `AppState::search_results`, which is replaced wholesale
    /// every time a search finishes — so the epoch is what tells us the index
    /// we're holding belongs to results that no longer exist.
    selected_result: Option<usize>,
    results_epoch: u64,
}

impl InspectorView {
    pub fn ui(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        self.sync(state);

        // Both axes, no auto-shrink: a fingerprint grid or a wide result row
        // should scroll inside the window rather than stretch it over its
        // neighbours.
        egui::ScrollArea::both()
            .auto_shrink([false, false])
            .show(ui, |ui| {
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
        // Larger than the dataset table's 64x48, which is read as a shape at a
        // glance. A result is being compared against the query, so it wants
        // enough room for the difference between two similar molecules to be
        // visible — but not so much that the next result is pushed off screen.
        let structure_size = Vec2::new(96.0, 72.0);
        // Same reasoning as the table's override: at this size a structure is
        // read as a shape, so hydrogens are dropped and it fits its box rather
        // than sharing a bond length with anything else. Carbons stay on
        // `Default` so a lone carbon isn't an empty box.
        let options = StructureOptions {
            padding: 2.0,
            show_carbons: ShowCarbons::Default,
            explicit_hydrogens: false,
            scale: 0.0,
            bond_length: 0.0,
            ..state.display.structure
        };
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

                    // Drawn only where coordinates already exist. Generating
                    // them needs the dataset mutably, which this cannot have
                    // while reading it — the same constraint the dataset table
                    // works under. Run 2D Coordinates in Operations, or load an
                    // SDF, which brings its own.
                    if mol.has_coords() {
                        ui.add(StructureView::new(mol, structure_size).with_options(options));
                    } else {
                        ui.add_sized(
                            structure_size,
                            egui::Label::new(RichText::new("\u{2014}").weak()),
                        );
                    }
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
                            .button(if is_selected {
                                "\u{25b2} Hide"
                            } else {
                                "\u{25bc} Why?"
                            })
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
                // The row now says what the molecule *is*, so the expansion
                // answers the other question: why it scored what it did. That
                // is the fingerprints, since Tanimoto is computed from exactly
                // these bits and nothing else. The atom and bond lists that
                // used to be here are in a detail window, which shows them
                // better and is a click away in the Datasets table.
                ui.indent(format!("result_{}", rank), |ui| {
                    if let Some(fingerprint) = state.dataset_fingerprints.get(idx) {
                        fingerprint_comparison(ui, fingerprint, state.query_fingerprint.as_ref());
                    } else {
                        ui.label(
                            RichText::new("No fingerprint for this molecule.")
                                .small()
                                .weak(),
                        );
                    }
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

        if ui
            .button("Export SVG")
            .on_hover_text("Save the query structure as an SVG file")
            .clicked()
        {
            let svg = structure_to_svg(
                &mol,
                egui::vec2(360.0, 300.0),
                &state.display.structure,
                &StructureTheme::light(),
            );
            save_svg(&suggested_filename(&state.query_source), &svg);
        }

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
