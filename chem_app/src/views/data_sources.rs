//! Loaded files and the active dataset's table of molecules.
//!
//! Still draws the fingerprint controls and the structure-display section that
//! the left panel has always drawn — #105 and #106 move those to the views that
//! will own them. The parameters those controls edit already live elsewhere
//! ([`OperationsView::fp_params`], [`AppState::display`]), so the move is a
//! matter of relocating the widgets, not of untangling state.
//!
//! [`OperationsView::fp_params`]: crate::views::OperationsView::fp_params

use crate::state::{AppState, FingerprintParams};
use crate::structure_view::{
    ShowCarbons, StructureOptions, StructureView, structure_display_section,
};
use chemcore::layout::ensure_coords;
use egui::{RichText, Vec2};

/// Rows of the dataset table drawn at once. Datasets run to tens of thousands
/// of molecules, so the table shows a window rather than all of them — and it
/// bounds how many structures a thumbnail pass has to lay out.
const MAX_TABLE_ROWS: usize = 20;

#[derive(Default)]
pub struct DataSourcesView;

impl DataSourcesView {
    pub fn ui(
        &mut self,
        ui: &mut egui::Ui,
        state: &mut AppState,
        fp_params: &mut FingerprintParams,
    ) {
        // Both axes, and no auto-shrink. `ScrollArea::vertical` leaves the
        // horizontal axis disabled, and egui sizes a disabled axis to its
        // content — so the seven-column table widened the window until it
        // covered its neighbour instead of scrolling. One region for the whole
        // window, too, rather than a scrolling table inside it: as a resizable
        // window this can be short, and the controls above the table would
        // otherwise be unreachable.
        egui::ScrollArea::both()
            .auto_shrink([false, false])
            .show(ui, |ui| self.contents(ui, state, fp_params));
    }

    fn contents(
        &mut self,
        ui: &mut egui::Ui,
        state: &mut AppState,
        fp_params: &mut FingerprintParams,
    ) {
        // No heading: the window's title bar names it now.
        ui.horizontal(|ui| {
            if ui.button("📂 Load File").clicked() {
                state.load_dataset_from_file();
            }
            if ui.button("📋 Load Examples").clicked() {
                state.load_example_dataset();
            }
        });

        ui.label(&state.dataset_status);
        ui.separator();

        ui.heading("Loaded Files");
        let active_index = state.loaded_files.active_index();
        let names: Vec<String> = state.loaded_files.names().map(str::to_owned).collect();
        let mut clicked_index = None;
        for (i, name) in names.iter().enumerate() {
            if ui.selectable_label(i == active_index, name).clicked() {
                clicked_index = Some(i);
            }
        }
        if let Some(i) = clicked_index {
            state.activate_loaded_file(i);
        }
        ui.separator();

        ui.heading("Fingerprint Settings");
        ui.add(egui::Slider::new(&mut fp_params.radius, 0..=5).text("Radius"));
        ui.add(
            egui::Slider::new(&mut fp_params.size, 512..=4096)
                .text("Size")
                .logarithmic(true),
        );

        if ui.button("⚡ Compute Fingerprints").clicked() {
            state.precompute_dataset_fingerprints(*fp_params);
        }
        if ui.button("🔬 Detect Aromaticity").clicked() {
            state.detect_aromaticity_for_dataset();
        }

        ui.separator();

        if state.loaded_files.active_dataset().is_empty() {
            return;
        }

        ui.label(format!(
            "Molecules in dataset: {}",
            state.loaded_files.active_dataset().len()
        ));
        ui.label(format!(
            "Fingerprints computed: {}",
            state.dataset_fingerprints.len()
        ));
        structure_display_section(
            ui,
            &mut state.display.structure,
            &mut state.display.show_thumbnails,
        );

        ui.separator();

        // Laid out before anything is drawn, so the mutable pass over the
        // dataset finishes before the draw code borrows it immutably.
        prepare_thumbnails(state);
        self.table(ui, state);
    }

    fn table(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let num_fingerprints = state.dataset_fingerprints.len();
        let selected_row = state.selected_row;
        let active_dataset = state.loaded_files.active_dataset();
        let shown = active_dataset.len().min(MAX_TABLE_ROWS);

        // A thumbnail is ~64x48, so it's read as a shape rather than for its
        // labels: hydrogens are dropped and the structure fits its cell rather
        // than using a shared bond length. Fixed-length sizing keeps a column
        // comparable (#78) but lets a large molecule overflow, and at this size
        // containment matters more.
        //
        // Carbons stay on Default rather than None: None would leave a lone
        // carbon with neither a label nor a bond, so methane's cell would draw
        // empty.
        let thumbnail_options = StructureOptions {
            padding: 2.0,
            show_carbons: ShowCarbons::Default,
            explicit_hydrogens: false,
            scale: 0.0,
            bond_length: 0.0,
            ..state.display.structure
        };
        let show_thumbnails = state.display.show_thumbnails;
        let mut clicked_row = None;

        egui::Grid::new("dataset_table")
            .num_columns(if show_thumbnails { 7 } else { 6 })
            .spacing([8.0, 4.0])
            .striped(true)
            .show(ui, |ui| {
                if show_thumbnails {
                    ui.label(RichText::new("Structure").strong());
                }
                ui.label(RichText::new("Name").strong());
                ui.label(RichText::new("SMILES").strong());
                ui.label(RichText::new("Formula").strong());
                ui.label(RichText::new("MW").strong());
                ui.label(RichText::new("Fingerprint").strong());
                ui.label(RichText::new("Aromatic").strong());
                ui.end_row();

                for i in 0..shown {
                    let mol = &active_dataset.molecules[i];
                    let is_selected = selected_row == Some(i);

                    if show_thumbnails {
                        ui.add(
                            StructureView::new(mol, Vec2::new(64.0, 48.0))
                                .with_options(thumbnail_options),
                        );
                    }
                    if ui
                        .selectable_label(is_selected, &active_dataset.names[i])
                        .clicked()
                    {
                        clicked_row = Some(i);
                    }
                    ui.label(RichText::new(&active_dataset.smiles[i]).code().small());
                    ui.label(mol.formula());
                    ui.label(format!("{:.2}", mol.molecular_weight()));
                    ui.label(if i < num_fingerprints { "Yes" } else { "No" });
                    let is_aromatic = mol.atoms().iter().any(|atom| atom.is_aromatic());
                    ui.label(if is_aromatic { "Yes" } else { "No" });
                    ui.end_row();
                }
            });

        if active_dataset.len() > shown {
            ui.label(format!("... and {} more", active_dataset.len() - shown));
        }

        if let Some(i) = clicked_row {
            state.selected_row = if selected_row == Some(i) {
                None
            } else {
                Some(i)
            };
        }
    }
}

/// Lays out the molecules the dataset table is about to draw.
///
/// Only the visible window is touched — a dataset can hold tens of thousands of
/// molecules, and laying all of them out to show twenty would stall the frame.
/// `ensure_coords` returns immediately for anything already positioned, so this
/// costs nothing after the first frame and leaves SDF-supplied coordinates
/// alone.
fn prepare_thumbnails(state: &mut AppState) {
    if !state.display.show_thumbnails {
        return;
    }
    let dataset = state.loaded_files.active_dataset_mut();
    for mol in dataset.molecules.iter_mut().take(MAX_TABLE_ROWS) {
        ensure_coords(mol);
    }
}
