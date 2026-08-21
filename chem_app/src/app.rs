//! The app struct and its frame loop.
//!
//! Everything shared lives in [`AppState`]; everything a single view cares about
//! lives in that view. What's left here is the shell: the top bar, the order the
//! views are drawn in, and the per-frame bookkeeping that isn't any view's job.

use crate::state::AppState;
use crate::views::{DataSourcesView, DetailView, OperationsView, VisualizationView};
use egui::{Color32, RichText};

/// How many frames the FPS readout averages over.
const FPS_WINDOW: usize = 60;

/// The views the shell draws, in the order it draws them.
///
/// A plain struct for now; #103 replaces it with a registry that can open and
/// close them and knows how to iterate them.
#[derive(Default)]
struct Views {
    data_sources: DataSourcesView,
    operations: OperationsView,
    visualization: VisualizationView,
    detail: DetailView,
}

pub struct ChemFpDemoApp {
    state: AppState,
    views: Views,
    fps_counter: f64,
    frame_times: Vec<f64>,
}

impl ChemFpDemoApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            state: AppState::new(),
            views: Views::default(),
            fps_counter: 0.0,
            frame_times: Vec::with_capacity(FPS_WINDOW),
        }
    }

    fn track_frame_rate(&mut self, ctx: &egui::Context) {
        let frame_time = ctx.input(|i| i.stable_dt as f64);
        self.frame_times.push(frame_time);
        if self.frame_times.len() > FPS_WINDOW {
            self.frame_times.remove(0);
        }
        let avg_frame_time = self.frame_times.iter().sum::<f64>() / self.frame_times.len() as f64;
        self.fps_counter = if avg_frame_time > 0.0 {
            1.0 / avg_frame_time
        } else {
            0.0
        };
    }

    fn top_panel(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("🧪 ChemFP Demo - Molecular Fingerprint Search");

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(format!("FPS: {:.0}", self.fps_counter));
                    ui.separator();

                    let using_gpu = self.state.search_engine.is_using_gpu();
                    let gpu_error = self.state.search_engine.gpu_init_error().map(str::to_owned);

                    // GPU chip: green when active or available-but-unselected,
                    // red when a real init attempt failed. Clicking switches
                    // to it if a GPU context already exists, or kicks off a
                    // (re)init attempt otherwise.
                    let gpu_response = if let Some(err) = &gpu_error {
                        ui.selectable_label(
                            using_gpu,
                            RichText::new("⚠ GPU").color(Color32::from_rgb(220, 50, 50)),
                        )
                        .on_hover_text(format!("GPU unavailable: {}", err))
                    } else {
                        ui.selectable_label(
                            using_gpu,
                            RichText::new("🚀 GPU").color(Color32::from_rgb(50, 200, 50)),
                        )
                    };
                    if gpu_response.clicked() && !self.state.search_engine.force_gpu() {
                        self.state.retry_gpu();
                    }

                    // CPU chip: always available, never styled as an alert —
                    // CPU is a working, legitimate mode, whether chosen
                    // deliberately or by GPU being unavailable.
                    if ui
                        .selectable_label(
                            !using_gpu,
                            RichText::new("💻 CPU").color(Color32::from_rgb(50, 200, 50)),
                        )
                        .clicked()
                    {
                        self.state.search_engine.force_cpu();
                    }
                });
            });
        });
    }
}

impl eframe::App for ChemFpDemoApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Results of anything that finished since the last frame are applied
        // before a view can read them, so no view sees half of one.
        self.state.poll_pending_work();
        self.track_frame_rate(ctx);
        // Deliberately not inside the operations view's `ui`: a debounce that
        // only ran while its view was drawn would stop working as soon as views
        // can be closed (#103).
        self.views.operations.tick(ctx, &mut self.state);

        self.top_panel(ctx);

        let state = &mut self.state;
        let views = &mut self.views;

        egui::SidePanel::left("dataset_panel")
            .default_width(300.0)
            .show(ctx, |ui| {
                views
                    .data_sources
                    .ui(ui, state, &mut views.operations.fp_params)
            });

        egui::SidePanel::right("query_panel")
            .default_width(400.0)
            .show(ctx, |ui| views.operations.ui(ui, state));

        egui::CentralPanel::default().show(ctx, |ui| views.visualization.ui(ui, state));

        views.detail.show(ctx, state);
    }
}
