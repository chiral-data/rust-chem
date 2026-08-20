//! Query entry and similarity search.
//!
//! Owns the operation *parameters* — the SMILES being typed, the fingerprint
//! radius and width, how many hits to return — while [`AppState`] owns what
//! running an operation produces. #105 gathers the rest of the operations here;
//! for now the fingerprint controls are still drawn by the data panel, reading
//! [`OperationsView::fp_params`].

use crate::state::{AppState, FingerprintParams, format_elapsed_ms};
use crate::structure_view::structure_panel_with_options;
use egui::Color32;
use web_time::{Duration, Instant};

use crate::fingerprint_view::fingerprint_full;
use crate::molecule_view::show_molecule_info;

/// How long the query SMILES box must sit idle before we run fingerprint
/// generation, so a blocking GPU dispatch doesn't fire on every keystroke.
const QUERY_DEBOUNCE: Duration = Duration::from_millis(300);

pub struct OperationsView {
    query_smiles: String,
    /// When the query box was last edited, or `None` if the edit has already
    /// been acted on. Drives the debounce in [`OperationsView::tick`].
    query_dirty_since: Option<Instant>,
    pub fp_params: FingerprintParams,
    top_k: usize,
}

impl Default for OperationsView {
    fn default() -> Self {
        Self {
            query_smiles: String::from("c1ccccc1"),
            query_dirty_since: None,
            fp_params: FingerprintParams::default(),
            top_k: 10,
        }
    }
}

impl OperationsView {
    /// Fires the debounced query parse once the box has been idle long enough.
    ///
    /// Called every frame from the frame loop rather than from `ui` below,
    /// because a debounce that only runs while its view is being drawn would
    /// silently stop working the moment the view can be closed (#103).
    pub fn tick(&mut self, ctx: &egui::Context, state: &mut AppState) {
        let Some(dirty_since) = self.query_dirty_since else {
            return;
        };
        let elapsed = dirty_since.elapsed();
        if elapsed >= QUERY_DEBOUNCE {
            self.query_dirty_since = None;
            state.parse_query(&self.query_smiles, self.fp_params);
        } else {
            // Wake up exactly when the debounce window closes instead of
            // repainting every frame while the query box sits idle.
            ctx.request_repaint_after(QUERY_DEBOUNCE - elapsed);
        }
    }

    pub fn ui(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        ui.heading("Query Molecule");
        ui.separator();

        ui.horizontal(|ui| {
            ui.label("SMILES:");
            let response = ui.text_edit_singleline(&mut self.query_smiles);

            if response.changed() {
                self.query_dirty_since = Some(Instant::now());
            }
            if ui.button("Parse").clicked() {
                self.query_dirty_since = None;
                state.parse_query(&self.query_smiles, self.fp_params);
            }
        });

        if let Some(error) = &state.query_error {
            ui.colored_label(Color32::RED, error);
        }

        if let Some(mol) = state.query_molecule.clone() {
            ui.separator();
            // Above the text details: seeing the structure is how you tell at a
            // glance whether the SMILES you typed is the molecule you meant.
            structure_panel_with_options(ui, &mol, 180.0, state.display.structure);
            show_molecule_info(ui, &mol, &self.query_smiles, "Query");
        }

        if let Some(fp) = &state.query_fingerprint {
            ui.separator();
            fingerprint_full(ui, fp);

            if let Some(time) = state.last_fp_gen_time {
                ui.label(format!("Generated in {}", format_elapsed_ms(time)));
            }
        }

        ui.separator();

        ui.horizontal(|ui| {
            ui.label("Top K:");
            ui.add(egui::Slider::new(&mut self.top_k, 1..=50));
        });

        if ui
            .add_enabled(state.can_search(), egui::Button::new("🔍 Search"))
            .clicked()
        {
            state.run_search(self.top_k);
        }

        if let Some(time) = state.last_search_time {
            ui.label(format!("Search completed in {}", format_elapsed_ms(time)));
        }
    }
}
