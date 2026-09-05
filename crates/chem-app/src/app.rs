//! The app struct and its frame loop.
//!
//! Everything shared lives in [`AppState`]; everything a single view cares about
//! lives in that view; which windows exist and whether they are open lives in
//! [`WindowRegistry`]. What's left here is the shell — the menu bar, the empty
//! workspace the windows float over, and the per-frame bookkeeping that is no
//! view's job.

use crate::state::AppState;
use crate::views::{DatasetsView, DetailView, InspectorView, OperationsView, SettingsView};
use crate::windows::WindowRegistry;
use egui::{Color32, RichText};

/// How many frames the FPS readout averages over.
const FPS_WINDOW: usize = 60;

/// Storage key for [`Persisted`]. Changing it discards saved state rather than
/// misreading it.
const PERSISTED_KEY: &str = "chem_workbench_v1";

/// What survives a restart.
///
/// Deliberately not here: loaded datasets, computed fingerprints, search
/// results, the query text, and which molecule windows were open. Those are
/// data rather than preferences, and restoring a fingerprint cache raises
/// questions a layout does not — one computed under a different radius is worse
/// than none.
///
/// Window geometry is also absent, because it is egui's: `Memory::areas` is
/// persisted by egui itself, so there is no second copy here to fall out of step
/// with where a window actually is.
#[derive(Default, serde::Serialize, serde::Deserialize)]
#[serde(default)]
struct Persisted {
    open_windows: crate::windows::OpenWindows,
    display: crate::state::DisplaySettings,
    fingerprint: crate::state::FingerprintParams,
    top_k: usize,
}

/// The views, each owning its own UI state.
///
/// Held separately from [`WindowRegistry`], which owns the *shells* those views
/// are drawn into: a view doesn't know or care whether it's in a window.
#[derive(Default)]
struct Views {
    datasets: DatasetsView,
    operations: OperationsView,
    inspector: InspectorView,
    settings: SettingsView,
    detail: DetailView,
}

pub struct WorkbenchApp {
    state: AppState,
    views: Views,
    windows: WindowRegistry,
    fps_counter: f64,
    frame_times: Vec<f64>,
}

impl WorkbenchApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self {
            // The frame loop itself, so work that finishes off-frame can wake
            // it rather than waiting for the user's next input (#186).
            state: AppState::new(&cc.egui_ctx),
            views: Views::default(),
            windows: WindowRegistry::default(),
            fps_counter: 0.0,
            frame_times: Vec::with_capacity(FPS_WINDOW),
        };

        // A stored value from an older build must not wedge the app, so a failed
        // decode is simply nothing: `eframe::get_value` returns None and the
        // defaults already in place stand.
        if let Some(saved) = cc
            .storage
            .and_then(|storage| eframe::get_value::<Persisted>(storage, PERSISTED_KEY))
        {
            app.windows.set_open_windows(saved.open_windows);
            app.state.display = saved.display;
            app.views.operations.fp_params = saved.fingerprint;
            // Zero would mean a search that returns nothing, which a saved
            // file from a future schema could otherwise impose.
            app.views.operations.top_k = saved.top_k.max(1);
        }

        app
    }

    /// Restores the default arrangement.
    ///
    /// Two halves, because the two halves are owned by different places: the
    /// open flags are ours, and the geometry is egui's. Without the second, a
    /// window dragged mostly off-screen would come back exactly where it was,
    /// which is the situation this exists for.
    fn reset_layout(&mut self, ctx: &egui::Context) {
        self.windows.reset_open_windows();
        ctx.memory_mut(|memory| memory.reset_areas());
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

    fn menu_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                // On web there is no OS title bar, so the menu bar is the only
                // place the app names itself.
                ui.label(RichText::new("\u{1f9ea} Chem Workbench").strong());
                ui.separator();

                ui.menu_button("File", |ui| {
                    if ui.button("📂 Load File…").clicked() {
                        self.state.load_dataset_from_file();
                        ui.close();
                    }
                    if ui.button("📋 Load Examples").clicked() {
                        self.state.load_example_dataset();
                        ui.close();
                    }
                });

                // Generated from the registry rather than listing windows here,
                // so adding a window doesn't mean editing the menu.
                ui.menu_button("View", |ui| {
                    for entry in self.windows.entries_mut() {
                        let title = entry.title();
                        ui.checkbox(&mut entry.open, title);
                    }

                    ui.separator();
                    if ui.button("Reset layout").clicked() {
                        self.reset_layout(ui.ctx());
                        ui.close();
                    }

                    // Molecule detail windows are not in the registry — they
                    // aren't toggleable singletons — so this is how they are
                    // represented among the window controls at all. Rows can be
                    // clicked much faster than windows can be closed.
                    ui.separator();
                    let open = self.state.open_details().len();
                    if ui
                        .add_enabled(
                            open > 0,
                            egui::Button::new(format!("Close all molecule windows ({})", open)),
                        )
                        .clicked()
                    {
                        self.state.close_all_details();
                        ui.close();
                    }
                });

                // Beside the menus rather than out at the right-hand edge,
                // where it was findable only by someone who already knew to
                // look for it.
                //
                // Plain text, like File and View. The app does put emoji on
                // labels — but on action buttons inside windows and on the
                // status chips, not on items in this bar.
                //
                // A toggle for the Settings window, not a menu of its own: a
                // menu would overlay the canvas and close the moment a control
                // was touched, and the whole point of a window is adjusting a
                // setting while watching a structure change. It is a
                // `selectable_label` rather than a button so it still shows
                // whether the window is open.
                let settings_open = self.windows.settings.open;
                if ui
                    .selectable_label(settings_open, "Settings")
                    .on_hover_text("Theme, and how structures are drawn")
                    .clicked()
                {
                    self.windows.settings.open = !settings_open;
                }

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(format!("FPS: {:.0}", self.fps_counter));
                    ui.separator();
                    // Set by crates/chem-app/e2e.sh, absent from an ordinary build.
                    // Without it there is no way to tell from the app which
                    // commit you are looking at, which is how a stale bundle
                    // went unnoticed for eleven minutes of commits during
                    // v0.5.0.
                    if let Some(build) = option_env!("CHEM_BUILD_ID") {
                        ui.label(RichText::new(build).small().weak())
                            .on_hover_text("Build under test");
                        ui.separator();
                    }
                    self.backend_chips(ui);
                });
            });
        });
    }

    /// GPU/CPU status and quick toggle, right-aligned in the menu bar where it
    /// has always been. These carry behaviour, not just status: clicking
    /// switches backend, or kicks off an init attempt if there is no context to
    /// switch to.
    fn backend_chips(&mut self, ui: &mut egui::Ui) {
        let using_gpu = self.state.search_engine.is_using_gpu();
        let gpu_error = self.state.search_engine.gpu_init_error().map(str::to_owned);

        // GPU chip: green when active or available-but-unselected, red when a
        // real init attempt failed.
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

        // CPU chip: always available, never styled as an alert — CPU is a
        // working, legitimate mode, whether chosen deliberately or by GPU being
        // unavailable.
        if ui
            .selectable_label(
                !using_gpu,
                RichText::new("💻 CPU").color(Color32::from_rgb(50, 200, 50)),
            )
            .clicked()
        {
            self.state.search_engine.force_cpu();
        }
    }

    /// The workspace the windows float over. Deliberately empty: no content
    /// lives here, which is what makes a window free to be any size rather than
    /// squeezed into a panel's fixed width.
    fn workspace(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            // The one exception, and it isn't content: a bare canvas gives no
            // clue that the View menu is where windows come back from, so
            // closing the last one would otherwise look like a broken app.
            if self.windows.all_content_closed() {
                ui.centered_and_justified(|ui| {
                    ui.label(
                        RichText::new("No windows open — reopen them from the View menu")
                            .size(16.0)
                            .weak(),
                    );
                });
            }
        });
    }

    fn show_windows(&mut self, ctx: &egui::Context) {
        // Destructured so each window's contents can borrow just the view it
        // needs rather than the whole app.
        let Self {
            state,
            views,
            windows,
            ..
        } = self;
        let Views {
            datasets,
            operations,
            inspector,
            settings,
            detail,
        } = views;

        windows.datasets.show(ctx, |ui| datasets.ui(ui, state));
        windows.operations.show(ctx, |ui| operations.ui(ui, state));
        windows.inspector.show(ctx, |ui| inspector.ui(ui, state));
        windows.settings.show(ctx, |ui| settings.ui(ui, state));

        // Not in the registry and not in the View menu: these aren't toggleable
        // singletons. They open by clicking a table row, close with their own
        // button, and there are between zero and MAX_OPEN_DETAILS of them.
        // "Close all molecule windows" in the View menu is how they are
        // represented among the window controls at all.
        detail.show(ctx, state);
    }
}

impl eframe::App for WorkbenchApp {
    /// Called by eframe periodically and on exit.
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(
            storage,
            PERSISTED_KEY,
            &Persisted {
                open_windows: self.windows.open_windows(),
                display: self.state.display,
                fingerprint: self.views.operations.fp_params,
                top_k: self.views.operations.top_k,
            },
        );
    }

    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Results of anything that finished since the last frame are applied
        // before a view can read them, so no view sees half of one.
        self.state.poll_pending_work();
        // Applied every frame rather than on change: it is one comparison, and
        // it means the preference is the single source of truth even if egui's
        // own theme is altered from elsewhere.
        ctx.set_theme(self.state.display.theme);
        self.track_frame_rate(ctx);
        // Deliberately not inside the operations view's `ui`: a debounce that
        // only ran while its view was drawn would stop the moment that window
        // could be closed — which, as of this story, it can be.
        self.views.operations.tick(ctx, &mut self.state);

        self.menu_bar(ctx);
        // After the menu bar has taken its strip, so the windows tile the space
        // actually left for them.
        self.windows.ensure_layout(ctx.available_rect());
        self.workspace(ctx);
        self.show_windows(ctx);
    }
}
