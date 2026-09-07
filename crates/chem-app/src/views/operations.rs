//! The operations: what you run against the active dataset.
//!
//! One section per operation — its parameters, a button to run it, and what
//! happened last time. #99 found that these four share no signature worth a
//! trait; what they do share is being the things a user *does* here, and
//! reporting a duration and a backend when they're done. That's a place in the
//! UI, not a type.
//!
//! Owns the parameters. [`AppState`] owns what running an operation produces.

use crate::dataset::DatasetFormat;
use crate::save::save_text;
use crate::state::{AppState, FingerprintParams, OperationOutcome};
use chem::io::format;
use egui::{Color32, RichText};
use web_time::{Duration, Instant};

/// How long the query SMILES box must sit idle before we run fingerprint
/// generation, so a blocking GPU dispatch doesn't fire on every keystroke.
const QUERY_DEBOUNCE: Duration = Duration::from_millis(300);

pub struct OperationsView {
    pub(crate) fp_params: FingerprintParams,
    query_smiles: String,
    /// When the query box was last edited, or `None` if the edit has already
    /// been acted on. Drives the debounce in [`OperationsView::tick`].
    query_dirty_since: Option<Instant>,
    pub(crate) top_k: usize,
    /// What Convert writes. Not persisted: `Format` is an index into a table
    /// this build fixes at compile time, so a saved one would be a number whose
    /// meaning could move under it.
    convert_target: DatasetFormat,
}

impl Default for OperationsView {
    fn default() -> Self {
        Self {
            fp_params: FingerprintParams::default(),
            query_smiles: String::from("c1ccccc1"),
            query_dirty_since: None,
            top_k: 10,
            // Not SMILES: the default should say something useful before the
            // user touches it, and most datasets arrive as SMILES.
            convert_target: DatasetFormat::SDF,
        }
    }
}

impl OperationsView {
    /// Fires the debounced query parse once the box has been idle long enough.
    ///
    /// Called every frame from the frame loop rather than from `ui` below,
    /// because a debounce that only ran while its window was open would stop
    /// working the moment that window was closed.
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
        // Both axes, no auto-shrink: a wide fingerprint grid should scroll
        // inside the window rather than stretch it over its neighbours.
        egui::ScrollArea::both()
            .auto_shrink([false, false])
            .show(ui, |ui| self.contents(ui, state));
    }

    fn contents(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        self.backend_section(ui, state);
        ui.separator();
        self.fingerprints_section(ui, state);
        self.aromaticity_section(ui, state);
        self.coordinates_section(ui, state);
        self.search_section(ui, state);
        self.convert_section(ui, state);
        self.export_section(ui, state);
    }

    fn convert_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let outcome = outcome_line(&state.convert);
        let target = self.convert_target;
        let source = state.loaded_files.active_format();
        // Read before the section takes `state` mutably, like every other
        // section's outcome line.
        let losses = state.conversion_losses(target);

        section(ui, "Convert", false, outcome, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Write as:");
                egui::ComboBox::from_id_salt("convert_target")
                    .selected_text(target.label())
                    .show_ui(ui, |ui| {
                        // Writable only, the mirror of the load dialog's
                        // readable-only rule (#266): offering a target that
                        // cannot be written just fails later.
                        for format in format::all().filter(|f| f.can_write()) {
                            ui.selectable_value(&mut self.convert_target, format, format.label());
                        }
                    });
            });

            // Before converting, not after: knowing what it costs is the point,
            // and this is the panel's most important thing -- it used to render
            // small and grey, which is how you hide something.
            //
            // Amber rather than the failure red: this app already uses
            // (220, 120, 50) for "GPU unavailable" and "not comparable", which
            // is the right register. A conversion dropping what the format
            // cannot hold is correct behaviour, not an error.
            if losses.is_empty() {
                ui.label(
                    RichText::new(format!(
                        "{} keeps everything this dataset holds.",
                        target.label()
                    ))
                    .small()
                    .weak(),
                );
            } else {
                ui.label(
                    RichText::new("Will lose")
                        .strong()
                        .color(Color32::from_rgb(220, 120, 50)),
                );
                for (attribute, count) in &losses {
                    ui.label(
                        RichText::new(format!("    {attribute} — {count} molecule(s)"))
                            .color(Color32::from_rgb(220, 120, 50)),
                    );
                }
                // The registry knows why, for the pairs where both formats
                // claim the attribute and the conversion loses it anyway.
                if let Some(reason) = format::pair_loss_reason(source, target) {
                    ui.label(RichText::new(format!("    {reason}")).small().weak());
                }
            }

            if ui
                .button("⟳ Convert")
                .on_hover_text("Adds the result as a new dataset, leaving this one loaded")
                .clicked()
            {
                state.convert_dataset(target);
            }
        });
    }

    fn export_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let format = state.loaded_files.active_format();
        // Not an `OperationOutcome`: neither save arm can report whether the
        // file landed -- the browser's download is fire-and-forget -- so a
        // header claiming success would be true on one platform only.
        section(ui, "Export", false, (String::new(), false), |ui| {
            ui.label(
                RichText::new(format!(
                    "Writes the active dataset as {}. Convert first to change format.",
                    format.label()
                ))
                .small()
                .weak(),
            );
            if ui.button("💾 Export…").clicked()
                && let Some((name, text)) = state.export_active_dataset()
            {
                save_text(
                    &name,
                    &text,
                    (format.label(), format.extensions()),
                    "text/plain",
                );
            }
        });
    }

    fn backend_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        ui.horizontal(|ui| {
            ui.label(RichText::new("Backend").strong());

            let using_gpu = state.search_engine.is_using_gpu();
            if ui.radio(using_gpu, "🚀 GPU").clicked() {
                // `force_gpu` reports false when there is no GPU context to
                // switch to, which makes this a first (or renewed) attempt.
                if !state.search_engine.force_gpu() {
                    state.retry_gpu();
                }
            }
            if ui.radio(!using_gpu, "💻 CPU").clicked() {
                state.search_engine.force_cpu();
            }
        });

        // The reason as readable text, not a tooltip you have to know to hover
        // for. This is what the menu bar can't carry.
        if let Some(err) = state.search_engine.gpu_init_error().map(str::to_owned) {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    RichText::new(format!("GPU unavailable: {}", err))
                        .small()
                        .color(Color32::from_rgb(220, 120, 50)),
                );
                if ui.small_button("Retry").clicked() {
                    state.retry_gpu();
                }
            });
        }
    }

    fn fingerprints_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let outcome = outcome_line(&state.fingerprints);
        section(ui, "Fingerprints", true, outcome, |ui| {
            ui.add(egui::Slider::new(&mut self.fp_params.radius, 0..=5).text("Radius"));
            ui.add(
                egui::Slider::new(&mut self.fp_params.size, 512..=4096)
                    .text("Size")
                    .logarithmic(true),
            );
            if ui.button("⚡ Compute Fingerprints").clicked() {
                state.precompute_dataset_fingerprints(self.fp_params);
            }
        });
    }

    fn aromaticity_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        // Collapsed by default: no parameters, so the header line is the whole
        // story once it has run.
        let outcome = outcome_line(&state.aromaticity);
        section(ui, "Aromaticity", false, outcome, |ui| {
            if ui.button("🔬 Detect Aromaticity").clicked() {
                state.detect_aromaticity_for_dataset();
            }
        });
    }

    fn coordinates_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let outcome = outcome_line(&state.coordinates);
        section(ui, "2D Coordinates", false, outcome, |ui| {
            if ui.button("📐 Generate Coordinates").clicked() {
                state.generate_coordinates_for_dataset();
            }
            ui.label(
                RichText::new("Molecules with coordinates from a file keep them.")
                    .small()
                    .weak(),
            );
        });
    }

    fn search_section(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        let outcome = outcome_line(&state.search);
        section(ui, "Similarity Search", true, outcome, |ui| {
            // Controls only. The query's structure and fingerprint are drawn
            // in the Inspector window, which owns everything the
            // display options govern — so this section stays two rows however
            // large the molecule is.
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

            // With the input, since it is about what was typed.
            if let Some(error) = &state.query_error {
                ui.colored_label(Color32::RED, error);
            }

            let can_search = state.can_search();
            ui.horizontal(|ui| {
                ui.label("Top K:");
                ui.add(egui::Slider::new(&mut self.top_k, 1..=50));
                if ui
                    .add_enabled(can_search, egui::Button::new("🔍 Search"))
                    .clicked()
                {
                    state.run_search(self.top_k);
                }
            });

            // Says which prerequisite is missing rather than presenting a greyed
            // button and leaving you to work it out. Search needs two things,
            // and one of them is another operation in this same window.
            if !can_search {
                let missing = if state.query_fingerprint.is_none() {
                    "Needs a parsed query molecule"
                } else {
                    "Needs dataset fingerprints — run Fingerprints above"
                };
                ui.label(RichText::new(missing).small().weak());
            }
        });
    }
}

/// One operation: a collapsing header carrying its last outcome, and its
/// controls inside.
///
/// The outcome sits in the header so a collapsed section still reports itself —
/// which is the point of collapsing the ones that have no parameters.
fn section(
    ui: &mut egui::Ui,
    title: &str,
    default_open: bool,
    outcome: (String, bool),
    add_contents: impl FnOnce(&mut egui::Ui),
) {
    // Taken as a rendered line rather than a `&OperationOutcome`: the borrow of
    // `AppState` has to end before the closure, which needs it mutably to run
    // the operation.
    let (text, failed) = outcome;
    let summary = RichText::new(text).small();
    let summary = if failed {
        summary.color(Color32::from_rgb(220, 80, 50))
    } else {
        summary.weak()
    };

    egui::CollapsingHeader::new(RichText::new(title).strong())
        .default_open(default_open)
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| ui.label(summary));
            add_contents(ui);
        });
}

/// The header line for an operation's last outcome, taken before the section
/// borrows `AppState` mutably.
fn outcome_line(outcome: &OperationOutcome) -> (String, bool) {
    (outcome.summary(), outcome.failed())
}
