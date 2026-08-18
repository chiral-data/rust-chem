use crate::dataset::{DatasetFormat, LoadedFiles, MoleculeDataset};
use crate::fingerprint_view::{fingerprint_compact, fingerprint_full};
use crate::molecule_view::{molecule_compact, show_atom_list, show_bond_list, show_molecule_info};
use crate::search::{FingerprintSearch, SearchResult};
use crate::structure_view::structure_panel;
use bitvec::prelude::BitVec;
use chemcore::layout::ensure_coords;
use chemcore::molecule::Molecule;
use chemio::aromaticity::detect_aromaticity;
use chemio::smiles::parse_smiles;
use egui::{Color32, RichText};
// std::time::Instant panics at runtime on wasm32-unknown-unknown (no clock
// source there); web-time is a drop-in replacement that re-exports std on
// native and uses performance.now() on web.
use web_time::{Duration, Instant};

#[cfg(target_arch = "wasm32")]
use chemgpu::{GpuMorganFingerprint, GpuTanimoto};
#[cfg(target_arch = "wasm32")]
use std::{cell::RefCell, rc::Rc};

/// How long the query SMILES box must sit idle before we run fingerprint
/// generation, so a blocking GPU dispatch doesn't fire on every keystroke.
const QUERY_DEBOUNCE: Duration = Duration::from_millis(300);

// A pending async task's slot: `None` while in flight, filled with its
// result (and how long it took, in ms) once `spawn_local`'s future resolves.
#[cfg(target_arch = "wasm32")]
type PendingSlot<T> = Rc<RefCell<Option<(anyhow::Result<T>, f64)>>>;
#[cfg(target_arch = "wasm32")]
type PendingGpuInit = Rc<RefCell<Option<Result<(GpuMorganFingerprint, GpuTanimoto), String>>>>;
#[cfg(target_arch = "wasm32")]
type PendingFileLoad = Rc<RefCell<Option<(String, Vec<u8>)>>>;

/// Formats a duration for display, switching to microseconds below 1ms so
/// fast operations (a single small-molecule fingerprint, say) don't just
/// show as "0.00ms".
fn format_elapsed_ms(ms: f64) -> String {
    if ms < 1.0 {
        format!("{:.1}\u{b5}s", ms * 1000.0)
    } else {
        format!("{:.2}ms", ms)
    }
}

pub struct ChemFpDemoApp {
    loaded_files: LoadedFiles,
    dataset_fingerprints: Vec<BitVec>,
    dataset_status: String,
    search_engine: FingerprintSearch,
    query_smiles: String,
    query_molecule: Option<Molecule>,
    query_fingerprint: Option<BitVec>,
    query_error: Option<String>,
    search_results: Vec<SearchResult>,
    selected_result: Option<usize>,
    // Which row in the Data window's molecule table is expanded to show its
    // detail view. Indexes into the active dataset, so it's reset alongside
    // dataset_fingerprints/search_results whenever the active dataset changes.
    selected_dataset_row: Option<usize>,
    // The detail window rebuilds its contents every frame, and a molecule from
    // SMILES has to be laid out before it can be drawn, so the laid-out copy is
    // cached against the row it came from rather than recomputed each frame.
    detail_molecule: Option<(usize, Molecule)>,
    fp_radius: u32,
    fp_size: u32,
    top_k: usize,
    last_search_time: Option<f64>,
    last_fp_gen_time: Option<f64>,
    fps_counter: f64,
    frame_times: Vec<f64>,
    query_dirty_since: Option<Instant>,
    // Browsers have no blocking main-thread file picker, so the wasm file
    // load has to happen on a spawned future and hand its result back here
    // to be picked up by the next `update()` poll.
    #[cfg(target_arch = "wasm32")]
    pending_file_load: PendingFileLoad,
    // GPU init can't happen inside `FingerprintSearch::new()` on wasm32 (see
    // its doc comment), so it's kicked off here instead and polled the same
    // way. Outer Option = has the attempt resolved yet; inner Option = did
    // it succeed.
    #[cfg(target_arch = "wasm32")]
    pending_gpu_init: PendingGpuInit,
    #[cfg(target_arch = "wasm32")]
    pending_dataset_fingerprints: PendingSlot<Vec<BitVec>>,
    #[cfg(target_arch = "wasm32")]
    pending_query_fingerprint: PendingSlot<BitVec>,
    #[cfg(target_arch = "wasm32")]
    pending_search_results: PendingSlot<Vec<SearchResult>>,
}

impl ChemFpDemoApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let dataset = MoleculeDataset::example_dataset().unwrap_or_default();
        let dataset_status = format!("Loaded {} example molecules", dataset.len());
        let loaded_files = LoadedFiles::new("Examples".to_string(), dataset, DatasetFormat::Smiles);

        #[cfg(target_arch = "wasm32")]
        let pending_gpu_init = Rc::new(RefCell::new(None));
        #[cfg(target_arch = "wasm32")]
        {
            let slot = pending_gpu_init.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = FingerprintSearch::try_init_gpu_async().await;
                *slot.borrow_mut() = Some(result);
            });
        }

        Self {
            loaded_files,
            dataset_fingerprints: Vec::new(),
            dataset_status,
            search_engine: FingerprintSearch::new(),
            query_smiles: String::from("c1ccccc1"),
            query_molecule: None,
            query_fingerprint: None,
            query_error: None,
            search_results: Vec::new(),
            selected_result: None,
            selected_dataset_row: None,
            detail_molecule: None,
            fp_radius: 2,
            fp_size: 2048,
            top_k: 10,
            last_search_time: None,
            last_fp_gen_time: None,
            fps_counter: 0.0,
            frame_times: Vec::with_capacity(60),
            query_dirty_since: None,
            #[cfg(target_arch = "wasm32")]
            pending_file_load: Rc::new(RefCell::new(None)),
            #[cfg(target_arch = "wasm32")]
            pending_gpu_init,
            #[cfg(target_arch = "wasm32")]
            pending_dataset_fingerprints: Rc::new(RefCell::new(None)),
            #[cfg(target_arch = "wasm32")]
            pending_query_fingerprint: Rc::new(RefCell::new(None)),
            #[cfg(target_arch = "wasm32")]
            pending_search_results: Rc::new(RefCell::new(None)),
        }
    }

    // Goes through AsyncFileDialog and reads content as bytes (rather than
    // FileDialog + a path) so file loading works on both native and web,
    // where there's no filesystem to read a path from.
    #[cfg(not(target_arch = "wasm32"))]
    fn load_dataset_from_file(&mut self) {
        // pollster::block_on keeps the dialog's existing blocking-call UX;
        // this is safe on native since blocking the calling thread doesn't
        // stop other threads from driving the future forward.
        let picked = pollster::block_on(async {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("SMILES", &["smi", "smiles", "txt"])
                .add_filter("SDF", &["sdf"])
                .pick_file()
                .await?;
            let name = file.file_name();
            Some((name, file.read().await))
        });

        if let Some((name, bytes)) = picked {
            self.apply_loaded_file_bytes(name, bytes);
        }
    }

    // Browsers have only one JS thread, and it also drives the file picker's
    // Promise machinery, so blocking on it would deadlock the tab. Spawn the
    // dialog as a non-blocking task instead and hand its result to
    // `pending_file_load`, polled from `update()` on the next frame.
    #[cfg(target_arch = "wasm32")]
    fn load_dataset_from_file(&mut self) {
        let slot = self.pending_file_load.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("SMILES", &["smi", "smiles", "txt"])
                .add_filter("SDF", &["sdf"])
                .pick_file()
                .await;
            if let Some(file) = file {
                let name = file.file_name();
                let bytes = file.read().await;
                *slot.borrow_mut() = Some((name, bytes));
            }
        });
    }

    fn apply_loaded_file_bytes(&mut self, name: String, bytes: Vec<u8>) {
        let content = match String::from_utf8(bytes) {
            Ok(content) => content,
            Err(e) => {
                self.dataset_status = "Failed to load file: not valid UTF-8".to_string();
                log::error!("Dataset load failed: {}", e);
                return;
            }
        };

        let format = DatasetFormat::from_filename(&name);
        let result = match format {
            DatasetFormat::Sdf => MoleculeDataset::load_from_sdf_str(&content),
            DatasetFormat::Smiles => MoleculeDataset::load_from_smiles_str(&content),
        };

        match result {
            Ok(dataset) => {
                self.dataset_status = format!(
                    "Loaded {} molecules from '{}' ({})",
                    dataset.len(),
                    name,
                    format.label()
                );
                self.loaded_files.add_and_activate(name, dataset, format);
                self.dataset_fingerprints.clear();
                self.search_engine.invalidate_target_dataset();
                self.search_results.clear();
                self.selected_result = None;
                self.selected_dataset_row = None;
                self.detail_molecule = None;
                log::info!("Dataset loaded successfully");
            }
            Err(e) => {
                self.dataset_status = format!("Failed to load file: {}", e);
                log::error!("Dataset load failed: {}", e);
            }
        }
    }

    fn load_example_dataset(&mut self) {
        match MoleculeDataset::example_dataset() {
            Ok(dataset) => {
                self.dataset_status = format!("Loaded {} example molecules", dataset.len());
                self.loaded_files.add_and_activate(
                    "Examples".to_string(),
                    dataset,
                    DatasetFormat::Smiles,
                );
                self.dataset_fingerprints.clear();
                self.search_engine.invalidate_target_dataset();
                self.search_results.clear();
                self.selected_result = None;
                self.selected_dataset_row = None;
                self.detail_molecule = None;
            }
            Err(e) => {
                self.dataset_status = format!("Failed to load examples: {}", e);
            }
        }
    }

    /// Switches the active dataset to an already-loaded entry, e.g. the user
    /// clicking a name in the Data window's loaded-files list. Runs the same
    /// reset as freshly loading a file, since fingerprints/search results
    /// belong to whichever dataset was active when they were computed.
    fn activate_loaded_file(&mut self, index: usize) {
        self.loaded_files.activate(index);
        self.dataset_status = format!(
            "Switched to '{}' ({} molecules)",
            self.loaded_files.names().nth(index).unwrap_or_default(),
            self.loaded_files.active_dataset().len()
        );
        self.dataset_fingerprints.clear();
        self.search_engine.invalidate_target_dataset();
        self.search_results.clear();
        self.selected_result = None;
        self.selected_dataset_row = None;
        self.detail_molecule = None;
    }

    fn precompute_dataset_fingerprints(&mut self) {
        if self.loaded_files.active_dataset().is_empty() {
            self.dataset_status = "No dataset loaded".to_string();
            return;
        }
        self.precompute_dataset_fingerprints_dispatch();
    }

    // CPU-only, no GPU implementation exists or is needed for this — it's a
    // simple ring search, nowhere near the cost of fingerprint generation.
    // Wired directly rather than through any new operation abstraction, per
    // this milestone's decision to defer that until there are more
    // operations to generalize from (see the v0.3.0 tracking issue).
    fn detect_aromaticity_for_dataset(&mut self) {
        let dataset = self.loaded_files.active_dataset_mut();
        if dataset.is_empty() {
            self.dataset_status = "No dataset loaded".to_string();
            return;
        }
        for mol in dataset.molecules.iter_mut() {
            detect_aromaticity(mol);
        }
        self.dataset_status = format!("Detected aromaticity for {} molecules", dataset.len());
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn precompute_dataset_fingerprints_dispatch(&mut self) {
        let start = Instant::now();
        let result = self.search_engine.generate_fingerprints_batch(
            &self.loaded_files.active_dataset().molecules,
            self.fp_radius,
            self.fp_size,
        );
        if let Ok(fps) = &result
            && let Err(e) = self.search_engine.set_target_dataset(fps)
        {
            log::warn!("Failed to upload dataset to GPU: {}", e);
        }
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        self.apply_dataset_fingerprints_result(result, elapsed_ms);
    }

    // Clones a snapshot of search_engine (cheap — see FingerprintSearch's
    // doc comment) to move into the spawned future instead of borrowing
    // `self` across the await boundary. The snapshot's own GPU-target-cache
    // upload (if any) is discarded once the task completes — `search_async`
    // re-uploads lazily as needed regardless, so this only costs one extra
    // upload on the next search rather than any actual bug.
    #[cfg(target_arch = "wasm32")]
    fn precompute_dataset_fingerprints_dispatch(&mut self) {
        let search_snapshot = self.search_engine.clone();
        let molecules = self.loaded_files.active_dataset().molecules.clone();
        let radius = self.fp_radius;
        let fp_size = self.fp_size;
        let slot = self.pending_dataset_fingerprints.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let start = Instant::now();
            let result = search_snapshot
                .generate_fingerprints_batch_async(&molecules, radius, fp_size)
                .await;
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            *slot.borrow_mut() = Some((result, elapsed_ms));
        });
    }

    fn apply_dataset_fingerprints_result(
        &mut self,
        result: anyhow::Result<Vec<BitVec>>,
        elapsed_ms: f64,
    ) {
        match result {
            Ok(fps) => {
                self.dataset_fingerprints = fps;
                self.dataset_status = format!(
                    "Computed {} fingerprints in {} ({} mode)",
                    self.dataset_fingerprints.len(),
                    format_elapsed_ms(elapsed_ms),
                    if self.search_engine.is_using_gpu() {
                        "GPU"
                    } else {
                        "CPU"
                    }
                );
                log::info!("Fingerprints computed in {:.2}ms", elapsed_ms);
            }
            Err(e) => {
                self.dataset_status = format!("Failed to compute fingerprints: {}", e);
                log::error!("Fingerprint computation failed: {}", e);
            }
        }
    }

    fn parse_query(&mut self) {
        let smiles = self.query_smiles.trim();
        if smiles.is_empty() {
            self.query_error = Some("SMILES string is empty".to_string());
            self.query_molecule = None;
            self.query_fingerprint = None;
            return;
        }

        match parse_smiles(smiles) {
            Ok(mol) => {
                self.query_molecule = Some(mol.clone());
                self.query_error = None;
                self.generate_query_fingerprint(mol);
            }
            Err(e) => {
                self.query_error = Some(format!("Invalid SMILES: {}", e));
                self.query_molecule = None;
                self.query_fingerprint = None;
            }
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn generate_query_fingerprint(&mut self, mol: Molecule) {
        let start = Instant::now();
        let result = self
            .search_engine
            .generate_fingerprint(&mol, self.fp_radius, self.fp_size);
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        self.apply_query_fingerprint_result(result, elapsed_ms);
    }

    #[cfg(target_arch = "wasm32")]
    fn generate_query_fingerprint(&mut self, mol: Molecule) {
        let search_snapshot = self.search_engine.clone();
        let radius = self.fp_radius;
        let fp_size = self.fp_size;
        let slot = self.pending_query_fingerprint.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let start = Instant::now();
            let result = search_snapshot
                .generate_fingerprint_async(&mol, radius, fp_size)
                .await;
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            *slot.borrow_mut() = Some((result, elapsed_ms));
        });
    }

    fn apply_query_fingerprint_result(&mut self, result: anyhow::Result<BitVec>, elapsed_ms: f64) {
        match result {
            Ok(fp) => {
                self.query_fingerprint = Some(fp);
                self.last_fp_gen_time = Some(elapsed_ms);
            }
            Err(e) => {
                self.query_error = Some(format!("Fingerprint generation failed: {}", e));
                self.query_fingerprint = None;
            }
        }
    }

    fn run_search(&mut self) {
        if self.dataset_fingerprints.is_empty() {
            self.query_error = Some("Please compute dataset fingerprints first".to_string());
            return;
        }

        if let Some(query_fp) = self.query_fingerprint.clone() {
            self.run_search_dispatch(query_fp);
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn run_search_dispatch(&mut self, query_fp: BitVec) {
        let start = Instant::now();
        let result = self
            .search_engine
            .search(&query_fp, &self.dataset_fingerprints, self.top_k);
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        self.apply_search_result(result, elapsed_ms);
    }

    #[cfg(target_arch = "wasm32")]
    fn run_search_dispatch(&mut self, query_fp: BitVec) {
        let mut search_snapshot = self.search_engine.clone();
        let target_fps = self.dataset_fingerprints.clone();
        let top_k = self.top_k;
        let slot = self.pending_search_results.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let start = Instant::now();
            let result = search_snapshot
                .search_async(&query_fp, &target_fps, top_k)
                .await;
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            *slot.borrow_mut() = Some((result, elapsed_ms));
        });
    }

    fn apply_search_result(&mut self, result: anyhow::Result<Vec<SearchResult>>, elapsed_ms: f64) {
        match result {
            Ok(results) => {
                self.search_results = results;
                self.last_search_time = Some(elapsed_ms);
                self.selected_result = None;
            }
            Err(e) => {
                self.query_error = Some(format!("Search failed: {}", e));
                self.search_results.clear();
            }
        }
    }

    // Kicks off a (re)attempt at GPU init, e.g. from the top bar's "Retry
    // GPU" button. Native does this synchronously (matching new()'s own
    // startup behavior); wasm32 can't block the browser's single JS thread,
    // so it spawns the same async task new() kicks off at startup and polls
    // pending_gpu_init the same way.
    #[cfg(not(target_arch = "wasm32"))]
    fn retry_gpu(&mut self) {
        if let Err(e) = self.search_engine.retry_gpu_init() {
            log::warn!("GPU retry failed: {}", e);
        } else {
            self.dataset_status = "GPU acceleration is now active".to_string();
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn retry_gpu(&mut self) {
        let slot = self.pending_gpu_init.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let result = FingerprintSearch::try_init_gpu_async().await;
            *slot.borrow_mut() = Some(result);
        });
    }

    fn top_panel(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("🧪 ChemFP Demo - Molecular Fingerprint Search");

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(format!("FPS: {:.0}", self.fps_counter));
                    ui.separator();

                    let using_gpu = self.search_engine.is_using_gpu();
                    let gpu_error = self.search_engine.gpu_init_error().map(str::to_owned);

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
                    if gpu_response.clicked() && !self.search_engine.force_gpu() {
                        self.retry_gpu();
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
                        self.search_engine.force_cpu();
                    }
                });
            });
        });
    }

    fn dataset_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("dataset_panel")
            .default_width(300.0)
            .show(ctx, |ui| {
                ui.heading("Dataset");
                ui.separator();

                ui.horizontal(|ui| {
                    if ui.button("📂 Load File").clicked() {
                        self.load_dataset_from_file();
                    }
                    if ui.button("📋 Load Examples").clicked() {
                        self.load_example_dataset();
                    }
                });

                ui.label(&self.dataset_status);
                ui.separator();

                ui.heading("Loaded Files");
                let active_index = self.loaded_files.active_index();
                let names: Vec<String> = self.loaded_files.names().map(str::to_owned).collect();
                let mut clicked_index = None;
                for (i, name) in names.iter().enumerate() {
                    if ui.selectable_label(i == active_index, name).clicked() {
                        clicked_index = Some(i);
                    }
                }
                if let Some(i) = clicked_index {
                    self.activate_loaded_file(i);
                }
                ui.separator();

                ui.heading("Fingerprint Settings");
                ui.add(egui::Slider::new(&mut self.fp_radius, 0..=5).text("Radius"));
                ui.add(
                    egui::Slider::new(&mut self.fp_size, 512..=4096)
                        .text("Size")
                        .logarithmic(true),
                );

                if ui.button("⚡ Compute Fingerprints").clicked() {
                    self.precompute_dataset_fingerprints();
                }
                if ui.button("🔬 Detect Aromaticity").clicked() {
                    self.detect_aromaticity_for_dataset();
                }

                ui.separator();

                let active_dataset = self.loaded_files.active_dataset();
                if !active_dataset.is_empty() {
                    ui.label(format!("Molecules in dataset: {}", active_dataset.len()));
                    ui.label(format!(
                        "Fingerprints computed: {}",
                        self.dataset_fingerprints.len()
                    ));

                    ui.separator();

                    let num_fingerprints = self.dataset_fingerprints.len();
                    let shown = active_dataset.len().min(20);
                    let mut clicked_row = None;

                    egui::ScrollArea::vertical().show(ui, |ui| {
                        egui::Grid::new("dataset_table")
                            .num_columns(6)
                            .spacing([8.0, 4.0])
                            .striped(true)
                            .show(ui, |ui| {
                                ui.label(RichText::new("Name").strong());
                                ui.label(RichText::new("SMILES").strong());
                                ui.label(RichText::new("Formula").strong());
                                ui.label(RichText::new("MW").strong());
                                ui.label(RichText::new("Fingerprint").strong());
                                ui.label(RichText::new("Aromatic").strong());
                                ui.end_row();

                                for i in 0..shown {
                                    let mol = &active_dataset.molecules[i];
                                    let is_selected = self.selected_dataset_row == Some(i);

                                    if ui
                                        .selectable_label(is_selected, &active_dataset.names[i])
                                        .clicked()
                                    {
                                        clicked_row = Some(i);
                                    }
                                    ui.label(
                                        RichText::new(&active_dataset.smiles[i]).code().small(),
                                    );
                                    ui.label(mol.formula());
                                    ui.label(format!("{:.2}", mol.molecular_weight()));
                                    ui.label(if i < num_fingerprints { "Yes" } else { "No" });
                                    let is_aromatic =
                                        mol.atoms().iter().any(|atom| atom.is_aromatic());
                                    ui.label(if is_aromatic { "Yes" } else { "No" });
                                    ui.end_row();
                                }
                            });

                        if active_dataset.len() > shown {
                            ui.label(format!("... and {} more", active_dataset.len() - shown));
                        }
                    });

                    if let Some(i) = clicked_row {
                        self.selected_dataset_row = if self.selected_dataset_row == Some(i) {
                            None
                        } else {
                            Some(i)
                        };
                    }
                }
            });
    }

    // A floating window rather than an inline expand below the table: the
    // Data window's side panel has a fixed height, and the detail view
    // (info grid + atom list + bond list) can be taller than that on a
    // short viewport with no way to scroll to the rest of it. A window is
    // independently movable/resizable and scrolls on its own if needed.
    fn molecule_detail_window(&mut self, ctx: &egui::Context) {
        let Some(i) = self.selected_dataset_row else {
            return;
        };
        let active_dataset = self.loaded_files.active_dataset();
        let Some(source) = active_dataset.molecules.get(i) else {
            return;
        };
        let name = active_dataset.names[i].clone();
        let smiles = active_dataset.smiles[i].clone();

        // Molecules parsed from SMILES carry no coordinates, so one is
        // generated here; SDF-sourced molecules keep the layout their file
        // supplied. Done once per selected row rather than every frame, since
        // laying out a large molecule isn't free.
        let needs_layout = !matches!(&self.detail_molecule, Some((cached, _)) if *cached == i);
        if needs_layout {
            let mut prepared = source.clone();
            ensure_coords(&mut prepared);
            self.detail_molecule = Some((i, prepared));
        }
        let Some((_, mol)) = &self.detail_molecule else {
            return;
        };
        let mol = mol.clone();

        // Caps the window itself to the available viewport height (minus
        // some margin) so it can never grow past the screen on a short
        // window/monitor — content beyond that scrolls inside instead of
        // the window just growing off-screen with no way to reach it.
        let max_height = (ctx.content_rect().height() - 80.0).max(200.0);

        let mut open = true;
        egui::Window::new(format!("Molecule: {}", name))
            .id(egui::Id::new("molecule_detail_window"))
            .open(&mut open)
            .resizable(true)
            .collapsible(true)
            .default_height(500.0)
            .max_height(max_height)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    structure_panel(ui, &mol, 220.0);
                    show_molecule_info(ui, &mol, &smiles, &name);
                    show_atom_list(ui, &mol);
                    show_bond_list(ui, &mol);
                });
            });

        if !open {
            self.selected_dataset_row = None;
            self.detail_molecule = None;
        }
    }

    fn query_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::right("query_panel")
            .default_width(400.0)
            .show(ctx, |ui| {
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
                        self.parse_query();
                    }
                });

                if let Some(ref error) = self.query_error {
                    ui.colored_label(Color32::RED, error);
                }

                if let Some(ref mol) = self.query_molecule {
                    ui.separator();
                    show_molecule_info(ui, mol, &self.query_smiles, "Query");
                }

                if let Some(ref fp) = self.query_fingerprint {
                    ui.separator();
                    fingerprint_full(ui, fp);

                    if let Some(time) = self.last_fp_gen_time {
                        ui.label(format!("Generated in {}", format_elapsed_ms(time)));
                    }
                }

                ui.separator();

                ui.horizontal(|ui| {
                    ui.label("Top K:");
                    ui.add(egui::Slider::new(&mut self.top_k, 1..=50));
                });

                if ui
                    .add_enabled(
                        self.query_fingerprint.is_some() && !self.dataset_fingerprints.is_empty(),
                        egui::Button::new("🔍 Search"),
                    )
                    .clicked()
                {
                    self.run_search();
                }

                if let Some(time) = self.last_search_time {
                    ui.label(format!("Search completed in {}", format_elapsed_ms(time)));
                }
            });
    }

    fn results_panel(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Search Results");
            ui.separator();

            if self.search_results.is_empty() {
                ui.centered_and_justified(|ui| {
                    ui.label(RichText::new("No search results yet").size(20.0).weak());
                });
            } else {
                let active_dataset = self.loaded_files.active_dataset();
                egui::ScrollArea::vertical().show(ui, |ui| {
                    for (rank, result) in self.search_results.iter().enumerate() {
                        let idx = result.index;
                        let mol = &active_dataset.molecules[idx];
                        let smiles = &active_dataset.smiles[idx];
                        let name = &active_dataset.names[idx];

                        let is_selected = self.selected_result == Some(rank);

                        let response = ui.group(|ui| {
                            ui.horizontal(|ui| {
                                ui.label(
                                    RichText::new(format!("#{}", rank + 1)).strong().size(16.0),
                                );
                                ui.separator();

                                ui.vertical(|ui| {
                                    molecule_compact(ui, smiles, name, &mol.formula());
                                    ui.label(format!(
                                        "Similarity: {:.3} ({:.1}%)",
                                        result.similarity,
                                        result.similarity * 100.0
                                    ));
                                });

                                ui.with_layout(
                                    egui::Layout::right_to_left(egui::Align::Center),
                                    |ui| {
                                        if ui
                                            .button(if is_selected {
                                                "▲ Hide"
                                            } else {
                                                "▼ Show"
                                            })
                                            .clicked()
                                        {
                                            self.selected_result =
                                                if is_selected { None } else { Some(rank) };
                                        }
                                    },
                                );
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

                                    if let Some(fp) = self.dataset_fingerprints.get(idx) {
                                        ui.vertical(|ui| {
                                            fingerprint_compact(ui, fp, "Molecule Fingerprint");

                                            if let Some(ref query_fp) = self.query_fingerprint {
                                                ui.separator();
                                                fingerprint_compact(
                                                    ui,
                                                    query_fp,
                                                    "Query Fingerprint",
                                                );
                                            }
                                        });
                                    }
                                });
                            });
                        }

                        ui.separator();
                    }
                });
            }
        });
    }
}

impl eframe::App for ChemFpDemoApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        #[cfg(target_arch = "wasm32")]
        {
            let loaded = self.pending_file_load.borrow_mut().take();
            if let Some((name, bytes)) = loaded {
                self.apply_loaded_file_bytes(name, bytes);
            }

            let gpu_init = self.pending_gpu_init.borrow_mut().take();
            match gpu_init {
                Some(Ok((morgan, tanimoto))) => {
                    self.search_engine.install_gpu(morgan, tanimoto);
                    self.dataset_status = "GPU acceleration is now active".to_string();
                }
                Some(Err(e)) => {
                    self.search_engine.record_gpu_init_failure(e);
                }
                None => {} // still pending
            }

            let fingerprints = self.pending_dataset_fingerprints.borrow_mut().take();
            if let Some((result, elapsed_ms)) = fingerprints {
                self.apply_dataset_fingerprints_result(result, elapsed_ms);
            }

            let query_fp = self.pending_query_fingerprint.borrow_mut().take();
            if let Some((result, elapsed_ms)) = query_fp {
                self.apply_query_fingerprint_result(result, elapsed_ms);
            }

            let search_results = self.pending_search_results.borrow_mut().take();
            if let Some((result, elapsed_ms)) = search_results {
                self.apply_search_result(result, elapsed_ms);
            }
        }

        let frame_time = ctx.input(|i| i.stable_dt as f64);
        self.frame_times.push(frame_time);
        if self.frame_times.len() > 60 {
            self.frame_times.remove(0);
        }
        let avg_frame_time = self.frame_times.iter().sum::<f64>() / self.frame_times.len() as f64;
        self.fps_counter = if avg_frame_time > 0.0 {
            1.0 / avg_frame_time
        } else {
            0.0
        };

        if let Some(dirty_since) = self.query_dirty_since {
            let elapsed = dirty_since.elapsed();
            if elapsed >= QUERY_DEBOUNCE {
                self.query_dirty_since = None;
                self.parse_query();
            } else {
                // Wake up exactly when the debounce window closes instead of
                // repainting every frame while the query box sits idle.
                ctx.request_repaint_after(QUERY_DEBOUNCE - elapsed);
            }
        }

        self.top_panel(ctx);
        self.dataset_panel(ctx);
        self.query_panel(ctx);
        self.results_panel(ctx);
        self.molecule_detail_window(ctx);
    }
}
