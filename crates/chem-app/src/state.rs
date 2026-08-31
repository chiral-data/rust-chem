//! Data shared by every view, and the operations that produce it.
//!
//! The app used to be one struct: 30 fields, and every panel a method taking
//! `&mut self`. That works for panels drawn in a fixed order and nothing else —
//! a view can't own its own state if all state lives on the app, and two views
//! of the same kind can't exist at all.
//!
//! What's here is the half that really is shared: the loaded datasets, what has
//! been computed from them, and the operations that do the computing. What's
//! *not* here is any view's own state — a text box's contents, a slider's
//! value, which row is expanded. Those live with the view that owns them, in
//! [`crate::views`], and are handed to these operations as arguments.

use crate::dataset::{DatasetFormat, LoadedFiles, MoleculeDataset};
use crate::task::Task;
use bitvec::prelude::BitVec;
use chem::core::layout::ensure_coords;
use chem::core::molecule::Molecule;
use chem::draw::structure::StructureOptions;
use chem::io::aromaticity::detect_aromaticity;
use chem::io::smiles::parse_smiles;
use chem::search::{FingerprintSearch, SearchResult};

#[cfg(target_arch = "wasm32")]
use chem::gpu::{GpuMorganFingerprint, GpuTanimoto};
#[cfg(target_arch = "wasm32")]
use std::{cell::RefCell, rc::Rc};

#[cfg(target_arch = "wasm32")]
type PendingGpuInit = Rc<RefCell<Option<Result<(GpuMorganFingerprint, GpuTanimoto), String>>>>;
#[cfg(target_arch = "wasm32")]
type PendingFileLoad = Rc<RefCell<Option<(String, Vec<u8>)>>>;

/// Formats a duration for display, switching to microseconds below 1ms so fast
/// operations (a single small-molecule fingerprint, say) don't just show as
/// "0.00ms".
pub fn format_elapsed_ms(ms: f64) -> String {
    if ms < 1.0 {
        format!("{:.1}\u{b5}s", ms * 1000.0)
    } else {
        format!("{:.2}ms", ms)
    }
}

/// What happened the last time one operation ran.
///
/// The app used to carry a single `dataset_status` string that loading, both
/// fingerprint paths, aromaticity and GPU init all wrote to, and one label read.
/// With an operation per section that's wrong twice over: running aromaticity
/// wiped the fingerprint result off the screen, and an operation that has moved
/// to another window was still reporting into the one it left.
#[derive(Clone, Debug, Default)]
pub struct OperationOutcome {
    /// What happened, in words. Empty until the operation has run at all.
    message: String,
    elapsed_ms: Option<f64>,
    /// Which backend ran it, where there is a choice. `None` for the CPU-only
    /// operations, whose backend isn't news.
    used_gpu: Option<bool>,
    failed: bool,
}

impl OperationOutcome {
    fn ok(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            ..Default::default()
        }
    }

    fn failure(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            failed: true,
            ..Default::default()
        }
    }

    fn timed(mut self, elapsed_ms: f64) -> Self {
        self.elapsed_ms = Some(elapsed_ms);
        self
    }

    fn on_gpu(mut self, used_gpu: bool) -> Self {
        self.used_gpu = Some(used_gpu);
        self
    }

    pub fn has_run(&self) -> bool {
        !self.message.is_empty()
    }

    pub fn failed(&self) -> bool {
        self.failed
    }

    /// One line for a collapsed section's header: what happened, how long it
    /// took, and which backend ran it, to whatever extent those are known.
    pub fn summary(&self) -> String {
        if !self.has_run() {
            return "\u{2014}".to_string(); // em dash: hasn't run
        }
        let mut summary = self.message.clone();
        if let Some(ms) = self.elapsed_ms {
            summary.push_str(&format!(" \u{b7} {}", format_elapsed_ms(ms)));
        }
        if let Some(used_gpu) = self.used_gpu {
            summary.push_str(if used_gpu {
                " \u{b7} GPU"
            } else {
                " \u{b7} CPU"
            });
        }
        summary
    }
}

/// How many molecule detail windows may be open at once.
///
/// Table rows can be clicked far faster than windows can be closed, so this is
/// bounded. Opening one past the cap closes the oldest rather than refusing the
/// click: a click that silently does nothing reads as a bug, and the window you
/// opened first is the one you are least likely to still be reading.
pub const MAX_OPEN_DETAILS: usize = 8;

/// Radius and bit-width for Morgan fingerprint generation.
///
/// Owned by whichever view offers the controls, and passed to the operations
/// that need it — a dataset's fingerprints and a query's have to be generated
/// under the same parameters to be comparable, so there is one set of them
/// rather than one per operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct FingerprintParams {
    pub radius: u32,
    pub size: u32,
}

impl Default for FingerprintParams {
    fn default() -> Self {
        Self {
            radius: 2,
            size: 2048,
        }
    }
}

/// How structures are drawn, anywhere in the app.
///
/// Shared rather than owned by the view that edits them: the dataset table's
/// thumbnails, the query structure and every detail window all draw with these,
/// so the controls being in one place doesn't make the values that view's
/// private business.
#[derive(Clone, Copy, Debug, Default, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct DisplaySettings {
    /// Light, dark, or follow the system.
    ///
    /// egui's own type, so `Context::set_theme` takes it directly and the
    /// picker is egui's own widget. Held here rather than read back out of
    /// egui because #108 persists this struct, and a preference has to be
    /// somewhere it can be saved from.
    pub theme: egui::ThemePreference,
    pub structure: StructureOptions,
    /// Whether the dataset table draws a structure per row. Off by default:
    /// it's a per-row render, and most of the time the table is being scanned
    /// for names and numbers rather than shapes.
    pub show_thumbnails: bool,
}

pub struct AppState {
    pub loaded_files: LoadedFiles,
    pub dataset_fingerprints: Vec<BitVec>,
    /// Dataset-level messages only — loaded, switched, failed to load. What the
    /// operations report goes in their own outcome, below.
    pub dataset_status: String,
    pub search_engine: FingerprintSearch,
    pub search_results: Vec<SearchResult>,

    /// The parsed query, its fingerprint, and why parsing failed if it did.
    ///
    /// Outputs, not inputs: the SMILES text being typed belongs to the view
    /// with the text box, but the molecule it parses to is drawn in one place
    /// and searched with in another.
    pub query_molecule: Option<Molecule>,
    pub query_fingerprint: Option<BitVec>,
    pub query_error: Option<String>,
    /// The SMILES that produced [`AppState::query_molecule`].
    ///
    /// Not the same thing as the text in the query box, which belongs to the
    /// view that owns the box and differs from this the moment you start
    /// typing. This is part of the *output*: whatever is drawing the parsed
    /// molecule needs the string it came from, and shouldn't have to reach into
    /// another view for it.
    pub query_source: String,

    /// Last outcome of each operation, reported by its own section.
    pub fingerprints: OperationOutcome,
    pub aromaticity: OperationOutcome,
    pub coordinates: OperationOutcome,
    pub query: OperationOutcome,
    pub search: OperationOutcome,

    pub display: DisplaySettings,

    /// Rows of the active dataset whose detail window is open, in the order
    /// they were opened.
    ///
    /// Shared because two views need it: the table highlights these rows, and
    /// the detail windows draw them. Order matters — it is what makes the cap
    /// close the *oldest* window rather than an arbitrary one.
    open_details: Vec<usize>,

    /// Bumped whenever the active dataset is replaced or swapped.
    ///
    /// Row indices, fingerprints and results all belong to whichever dataset
    /// was active when they were made. Rather than every dataset-changing path
    /// reaching into every view to clear what it holds — which is what the old
    /// single struct did, in a six-line block duplicated three times — the
    /// paths bump this, and a view holding indices notices its own state is
    /// stale. That also means a view which isn't being drawn this frame still
    /// finds out.
    dataset_epoch: u64,
    /// Bumped whenever new search results land, for the same reason.
    results_epoch: u64,

    // GPU-capable work, which is async on wasm32. Each is started in one place
    // and collected in `update()`; see `task::Task`.
    dataset_fingerprint_task: Task<Vec<BitVec>>,
    query_fingerprint_task: Task<BitVec>,
    search_task: Task<Vec<SearchResult>>,

    // Browsers have no blocking main-thread file picker, so the wasm file load
    // has to happen on a spawned future and hand its result back here to be
    // picked up by the next `update()` poll.
    #[cfg(target_arch = "wasm32")]
    pending_file_load: PendingFileLoad,
    // GPU init can't happen inside `FingerprintSearch::new()` on wasm32 (see
    // its doc comment), so it's kicked off here instead and polled the same
    // way. Outer Option = has the attempt resolved yet; inner Option = did it
    // succeed.
    #[cfg(target_arch = "wasm32")]
    pending_gpu_init: PendingGpuInit,
}

impl AppState {
    pub fn new() -> Self {
        Self::with_engine(FingerprintSearch::new())
    }

    /// A CPU-only instance, for tests that exercise state transitions and have
    /// no business depending on whether the machine running them has a GPU —
    /// or on GPU init running at all, which does not take kindly to being
    /// driven from parallel tests (#19). Mirrors
    /// [`FingerprintSearch::new_cpu_only`], which exists for the same reason.
    #[cfg(test)]
    fn cpu_only() -> Self {
        Self::with_engine(FingerprintSearch::new_cpu_only())
    }

    fn with_engine(search_engine: FingerprintSearch) -> Self {
        let dataset = MoleculeDataset::example_dataset();
        let dataset_status = format!("Loaded {} example molecules", dataset.len());
        let loaded_files = LoadedFiles::new("Examples".to_string(), dataset, DatasetFormat::SMILES);

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
            search_engine,
            search_results: Vec::new(),
            query_molecule: None,
            query_fingerprint: None,
            query_error: None,
            query_source: String::new(),
            fingerprints: OperationOutcome::default(),
            aromaticity: OperationOutcome::default(),
            coordinates: OperationOutcome::default(),
            query: OperationOutcome::default(),
            search: OperationOutcome::default(),
            display: DisplaySettings::default(),
            open_details: Vec::new(),
            dataset_epoch: 0,
            results_epoch: 0,
            dataset_fingerprint_task: Task::new(),
            query_fingerprint_task: Task::new(),
            search_task: Task::new(),
            #[cfg(target_arch = "wasm32")]
            pending_file_load: Rc::new(RefCell::new(None)),
            #[cfg(target_arch = "wasm32")]
            pending_gpu_init,
        }
    }

    /// Rows whose detail window is open, oldest first.
    pub fn open_details(&self) -> &[usize] {
        &self.open_details
    }

    pub fn is_detail_open(&self, row: usize) -> bool {
        self.open_details.contains(&row)
    }

    /// Opens a row's detail window, or closes it if it is already open.
    ///
    /// Toggling rather than raising an existing window: the table row is lit
    /// while its window is open, so a second click on a lit row reads as
    /// "put that away".
    pub fn toggle_detail(&mut self, row: usize) {
        if let Some(pos) = self.open_details.iter().position(|&r| r == row) {
            self.open_details.remove(pos);
        } else {
            if self.open_details.len() >= MAX_OPEN_DETAILS {
                self.open_details.remove(0);
            }
            self.open_details.push(row);
        }
    }

    pub fn close_detail(&mut self, row: usize) {
        self.open_details.retain(|&r| r != row);
    }

    pub fn close_all_details(&mut self) {
        self.open_details.clear();
    }

    /// Version of the active dataset. See [`AppState::dataset_epoch`] field docs.
    pub fn dataset_epoch(&self) -> u64 {
        self.dataset_epoch
    }

    /// Version of [`AppState::search_results`].
    pub fn results_epoch(&self) -> u64 {
        self.results_epoch
    }

    /// Drops everything derived from the dataset that just went away, and tells
    /// the views holding indices into it to do the same.
    fn invalidate_active_dataset(&mut self) {
        self.dataset_fingerprints.clear();
        self.search_engine.invalidate_target_dataset();
        self.search_results.clear();
        // Keyed by row index, which means nothing against a different dataset.
        self.open_details.clear();
        // These describe a dataset that is no longer active — "2048
        // fingerprints · GPU" against a dataset that has been swapped out is
        // worse than saying nothing. The query outcome survives: it is about
        // the query, not the dataset.
        self.fingerprints = OperationOutcome::default();
        self.aromaticity = OperationOutcome::default();
        self.coordinates = OperationOutcome::default();
        self.search = OperationOutcome::default();
        self.dataset_epoch += 1;
        self.results_epoch += 1;
    }

    // Goes through AsyncFileDialog and reads content as bytes (rather than
    // FileDialog + a path) so file loading works on both native and web, where
    // there's no filesystem to read a path from.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_dataset_from_file(&mut self) {
        // pollster::block_on keeps the dialog's existing blocking-call UX; this
        // is safe on native since blocking the calling thread doesn't stop
        // other threads from driving the future forward.
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
    pub fn load_dataset_from_file(&mut self) {
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

    pub fn apply_loaded_file_bytes(&mut self, name: String, bytes: Vec<u8>) {
        let content = match String::from_utf8(bytes) {
            Ok(content) => content,
            Err(e) => {
                self.dataset_status = "Failed to load file: not valid UTF-8".to_string();
                log::error!("Dataset load failed: {}", e);
                return;
            }
        };

        let format = DatasetFormat::from_filename(&name);
        let outcome = chem::io::reader::read(&content, format);
        let dataset = MoleculeDataset::from_outcome(&outcome);

        // Records that failed used to be logged and never surfaced, so a file
        // that half-loaded looked like a file that fully loaded. Reading now
        // reports them, so the status line can too.
        self.dataset_status = if outcome.skipped.is_empty() {
            format!(
                "Loaded {} molecules from '{}' ({})",
                dataset.len(),
                name,
                format.label()
            )
        } else {
            format!(
                "Loaded {} molecules from '{}' ({}) — {} skipped",
                dataset.len(),
                name,
                format.label(),
                outcome.skipped.len()
            )
        };
        for skipped in &outcome.skipped {
            log::warn!(
                "Skipped record {} in '{}': {}",
                skipped.position,
                name,
                skipped.error
            );
        }

        self.loaded_files.add_and_activate(name, dataset, format);
        self.invalidate_active_dataset();
    }

    pub fn load_example_dataset(&mut self) {
        let dataset = MoleculeDataset::example_dataset();
        self.dataset_status = format!("Loaded {} example molecules", dataset.len());
        self.loaded_files
            .add_and_activate("Examples".to_string(), dataset, DatasetFormat::SMILES);
        self.invalidate_active_dataset();
    }

    /// Switches the active dataset to an already-loaded entry, e.g. the user
    /// clicking a name in the loaded-files list. Runs the same invalidation as
    /// freshly loading a file, since fingerprints and search results belong to
    /// whichever dataset was active when they were computed.
    pub fn activate_loaded_file(&mut self, index: usize) {
        self.loaded_files.activate(index);
        self.dataset_status = format!(
            "Switched to '{}' ({} molecules)",
            self.loaded_files.names().nth(index).unwrap_or_default(),
            self.loaded_files.active_dataset().len()
        );
        self.invalidate_active_dataset();
    }

    /// Removes a loaded dataset, discarding what belonged to it if it was the
    /// active one.
    pub fn remove_loaded_file(&mut self, index: usize) {
        let Some(name) = self
            .loaded_files
            .entries()
            .get(index)
            .map(|entry| entry.name.clone())
        else {
            return;
        };

        if self.loaded_files.remove(index) {
            // A different dataset is active now, so everything derived from the
            // old one goes — the same reset that switching files performs.
            self.invalidate_active_dataset();
            self.dataset_status = format!(
                "Removed '{}' — now showing '{}' ({} molecules)",
                name,
                self.loaded_files
                    .names()
                    .nth(self.loaded_files.active_index())
                    .unwrap_or_default(),
                self.loaded_files.active_dataset().len()
            );
        } else {
            self.dataset_status = format!("Removed '{}'", name);
        }
    }

    pub fn precompute_dataset_fingerprints(&mut self, params: FingerprintParams) {
        if self.loaded_files.active_dataset().is_empty() {
            self.fingerprints = OperationOutcome::failure("No dataset loaded");
            return;
        }

        // Clones a snapshot of the search engine to hand to the task rather
        // than borrowing `self`, which a spawned future can't do. The clone is
        // cheap — the GPU handles behind it are `Arc`-backed — and its own GPU
        // target cache is discarded when the task ends, costing one extra
        // upload on the next search rather than any wrong answer.
        let engine = self.search_engine.clone();
        let molecules = self.loaded_files.active_dataset().molecules.clone();

        self.dataset_fingerprint_task.start(async move {
            engine
                .generate_fingerprints_batch_async(&molecules, params.radius, params.size)
                .await
        });
    }

    // CPU-only, no GPU implementation exists or is needed for this — it's a
    // simple ring search, nowhere near the cost of fingerprint generation.
    // Wired directly rather than through any operation abstraction: #99 found
    // that the four operations' signatures share nothing worth a trait.
    pub fn detect_aromaticity_for_dataset(&mut self) {
        let dataset = self.loaded_files.active_dataset_mut();
        if dataset.is_empty() {
            self.aromaticity = OperationOutcome::failure("No dataset loaded");
            return;
        }
        for mol in dataset.molecules.iter_mut() {
            detect_aromaticity(mol);
        }
        let aromatic = dataset
            .molecules
            .iter()
            .filter(|mol| mol.atoms().iter().any(|atom| atom.is_aromatic()))
            .count();
        self.aromaticity =
            OperationOutcome::ok(format!("{} of {} aromatic", aromatic, dataset.len()));
    }

    /// Generates 2D coordinates across the dataset.
    ///
    /// Was never an action: coordinates appeared as a side effect of selecting a
    /// row or turning on thumbnails, which made it the one operation of the four
    /// with no way to run it. Reports how many it generated against how many it
    /// left alone, since coordinates a file supplied are deliberately kept
    /// (#88) and that is otherwise invisible.
    pub fn generate_coordinates_for_dataset(&mut self) {
        let dataset = self.loaded_files.active_dataset_mut();
        if dataset.is_empty() {
            self.coordinates = OperationOutcome::failure("No dataset loaded");
            return;
        }

        let mut generated = 0;
        let mut kept = 0;
        for mol in dataset.molecules.iter_mut() {
            if mol.has_coords() {
                kept += 1;
            } else if ensure_coords(mol) {
                generated += 1;
            }
        }
        self.coordinates = OperationOutcome::ok(if kept > 0 {
            format!("{} generated, {} kept from file", generated, kept)
        } else {
            format!("{} generated", generated)
        });
    }

    pub fn parse_query(&mut self, smiles: &str, params: FingerprintParams) {
        let smiles = smiles.trim();
        if smiles.is_empty() {
            self.query_error = Some("SMILES string is empty".to_string());
            self.query = OperationOutcome::failure("No query");
            self.query_molecule = None;
            self.query_fingerprint = None;
            self.query_source.clear();
            return;
        }

        match parse_smiles(smiles) {
            Ok(mut mol) => {
                // SMILES has no geometry. Laying out here rather than at draw
                // time means it happens once per parse (already debounced)
                // instead of every frame.
                ensure_coords(&mut mol);
                self.query_molecule = Some(mol.clone());
                self.query_source = smiles.to_string();
                self.query_error = None;
                self.generate_query_fingerprint(mol, params);
            }
            Err(e) => {
                self.query_error = Some(format!("Invalid SMILES: {}", e));
                self.query = OperationOutcome::failure("Invalid SMILES");
                self.query_molecule = None;
                self.query_fingerprint = None;
                self.query_source.clear();
            }
        }
    }

    fn generate_query_fingerprint(&mut self, mol: Molecule, params: FingerprintParams) {
        let engine = self.search_engine.clone();

        self.query_fingerprint_task.start(async move {
            engine
                .generate_fingerprint_async(&mol, params.radius, params.size)
                .await
        });
    }

    pub fn run_search(&mut self, top_k: usize) {
        if self.dataset_fingerprints.is_empty() {
            self.search = OperationOutcome::failure("No dataset fingerprints");
            return;
        }

        let Some(query_fp) = self.query_fingerprint.clone() else {
            return;
        };

        let mut engine = self.search_engine.clone();
        let target_fps = self.dataset_fingerprints.clone();

        self.search_task
            .start(async move { engine.search_async(&query_fp, &target_fps, top_k).await });
    }

    /// True once a search has everything it needs: a query fingerprint to look
    /// for, and dataset fingerprints to look through.
    pub fn can_search(&self) -> bool {
        self.query_fingerprint.is_some() && !self.dataset_fingerprints.is_empty()
    }

    // Kicks off a (re)attempt at GPU init, e.g. from the top bar's GPU chip.
    // Native does this synchronously (matching new()'s own startup behavior);
    // wasm32 can't block the browser's single JS thread, so it spawns the same
    // async task new() kicks off at startup and polls pending_gpu_init the same
    // way.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn retry_gpu(&mut self) {
        // No status message: which backend is live is shown by the Operations
        // window's backend section and the menu bar chips, so it doesn't need
        // to be announced into a dataset's status line.
        if let Err(e) = self.search_engine.retry_gpu_init() {
            log::warn!("GPU retry failed: {}", e);
        }
    }

    #[cfg(target_arch = "wasm32")]
    pub fn retry_gpu(&mut self) {
        let slot = self.pending_gpu_init.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let result = FingerprintSearch::try_init_gpu_async().await;
            *slot.borrow_mut() = Some(result);
        });
    }

    /// Collects anything that finished since the last frame.
    ///
    /// Called once per frame, before any view is drawn, so a view never sees a
    /// half-applied result.
    pub fn poll_pending_work(&mut self) {
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
                }
                Some(Err(e)) => {
                    self.search_engine.record_gpu_init_failure(e);
                }
                None => {} // still pending
            }
        }

        // Collected on both platforms alike: the task ran on a spawned future
        // (wasm32) or inline (native), but either way its result is applied
        // here rather than at the call site.
        if let Some((result, elapsed_ms)) = self.dataset_fingerprint_task.poll() {
            self.apply_dataset_fingerprints_result(result, elapsed_ms);
        }
        if let Some((result, elapsed_ms)) = self.query_fingerprint_task.poll() {
            self.apply_query_fingerprint_result(result, elapsed_ms);
        }
        if let Some((result, elapsed_ms)) = self.search_task.poll() {
            self.apply_search_result(result, elapsed_ms);
        }
    }

    fn apply_dataset_fingerprints_result(
        &mut self,
        result: anyhow::Result<Vec<BitVec>>,
        elapsed_ms: f64,
    ) {
        match result {
            Ok(fps) => {
                self.dataset_fingerprints = fps;
                self.fingerprints = OperationOutcome::ok(format!(
                    "{} fingerprints",
                    self.dataset_fingerprints.len()
                ))
                .timed(elapsed_ms)
                .on_gpu(self.search_engine.is_using_gpu());
                log::info!("Fingerprints computed in {:.2}ms", elapsed_ms);
            }
            Err(e) => {
                self.fingerprints = OperationOutcome::failure(format!("Failed: {}", e));
                log::error!("Fingerprint computation failed: {}", e);
            }
        }
    }

    fn apply_query_fingerprint_result(&mut self, result: anyhow::Result<BitVec>, elapsed_ms: f64) {
        match result {
            Ok(fp) => {
                self.query_fingerprint = Some(fp);
                self.query = OperationOutcome::ok("Parsed")
                    .timed(elapsed_ms)
                    .on_gpu(self.search_engine.is_using_gpu());
            }
            Err(e) => {
                self.query_error = Some(format!("Fingerprint generation failed: {}", e));
                self.query = OperationOutcome::failure("Fingerprint failed");
                self.query_fingerprint = None;
            }
        }
    }

    fn apply_search_result(&mut self, result: anyhow::Result<Vec<SearchResult>>, elapsed_ms: f64) {
        match result {
            Ok(results) => {
                self.search = OperationOutcome::ok(format!("{} hits", results.len()))
                    .timed(elapsed_ms)
                    .on_gpu(self.search_engine.is_using_gpu());
                self.search_results = results;
            }
            Err(e) => {
                self.query_error = Some(format!("Search failed: {}", e));
                self.search = OperationOutcome::failure(format!("Failed: {}", e));
                self.search_results.clear();
            }
        }
        // Either way the old results are gone, so a view holding an index into
        // them has to let go of it.
        self.results_epoch += 1;
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Stands in for real derived data — the tests here care about *whether*
    /// it survives a dataset change, not what's in it.
    fn with_derived_data(state: &mut AppState) {
        state.dataset_fingerprints = vec![BitVec::new()];
        state.search_results = vec![SearchResult {
            index: 0,
            similarity: 1.0,
        }];
        state.toggle_detail(0);
    }

    #[test]
    fn test_a_row_opens_and_a_second_click_closes_it() {
        let mut state = AppState::cpu_only();

        state.toggle_detail(3);
        assert!(state.is_detail_open(3));
        assert_eq!(state.open_details(), [3]);

        state.toggle_detail(3);
        assert!(!state.is_detail_open(3));
        assert!(state.open_details().is_empty());
    }

    #[test]
    fn test_several_rows_stay_open_together() {
        let mut state = AppState::cpu_only();
        for row in [4, 1, 7] {
            state.toggle_detail(row);
        }

        // Comparing two molecules is the point; opening a second must not
        // evict the first, which is what one shared slot used to do.
        assert_eq!(state.open_details(), [4, 1, 7]);
    }

    #[test]
    fn test_closing_one_leaves_the_others() {
        let mut state = AppState::cpu_only();
        for row in [4, 1, 7] {
            state.toggle_detail(row);
        }

        state.close_detail(1);

        assert_eq!(state.open_details(), [4, 7]);
    }

    #[test]
    fn test_past_the_cap_the_oldest_window_closes() {
        let mut state = AppState::cpu_only();
        for row in 0..MAX_OPEN_DETAILS {
            state.toggle_detail(row);
        }
        assert_eq!(state.open_details().len(), MAX_OPEN_DETAILS);

        state.toggle_detail(99);

        // The click is never refused; the oldest gives way.
        assert_eq!(state.open_details().len(), MAX_OPEN_DETAILS);
        assert!(!state.is_detail_open(0), "the oldest should have closed");
        assert!(state.is_detail_open(99));
        assert!(state.is_detail_open(1), "the second-oldest should remain");
    }

    #[test]
    fn test_reopening_an_open_row_does_not_reorder_the_rest() {
        let mut state = AppState::cpu_only();
        for row in [2, 5] {
            state.toggle_detail(row);
        }

        // Closing and reopening moves it to the back of the queue, which is
        // what decides who the cap evicts next.
        state.toggle_detail(2);
        state.toggle_detail(2);

        assert_eq!(state.open_details(), [5, 2]);
    }

    #[test]
    fn test_close_all_clears_every_window() {
        let mut state = AppState::cpu_only();
        for row in 0..5 {
            state.toggle_detail(row);
        }

        state.close_all_details();

        assert!(state.open_details().is_empty());
    }

    #[test]
    fn test_starts_on_the_example_dataset_with_nothing_derived() {
        let state = AppState::cpu_only();

        assert!(!state.loaded_files.active_dataset().is_empty());
        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
    }

    #[test]
    fn test_loading_a_dataset_drops_what_the_old_one_derived() {
        let mut state = AppState::cpu_only();
        with_derived_data(&mut state);
        let epoch = state.dataset_epoch();

        state.load_example_dataset();

        // Fingerprints and results belong to whichever dataset was active when
        // they were computed, and a row index means nothing against a different
        // dataset.
        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
        // Views hold indices too, and find out the same way.
        assert!(state.dataset_epoch() > epoch);
    }

    #[test]
    fn test_switching_between_loaded_files_invalidates_the_same_way() {
        let mut state = AppState::cpu_only();
        state.load_example_dataset(); // a second entry to switch between
        with_derived_data(&mut state);
        let epoch = state.dataset_epoch();

        state.activate_loaded_file(0);

        assert_eq!(state.loaded_files.active_index(), 0);
        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
        assert!(state.dataset_epoch() > epoch);
    }

    #[test]
    fn test_a_failed_load_leaves_the_active_dataset_alone() {
        let mut state = AppState::cpu_only();
        with_derived_data(&mut state);
        let epoch = state.dataset_epoch();
        let files_before = state.loaded_files.names().count();

        // Not valid UTF-8, so it never reaches a parser.
        state.apply_loaded_file_bytes("broken.smi".to_string(), vec![0xff, 0xfe]);

        assert!(state.dataset_status.contains("Failed to load"));
        // Nothing was replaced, so nothing derived from it is stale.
        assert_eq!(state.loaded_files.names().count(), files_before);
        assert_eq!(state.dataset_epoch(), epoch);
        assert_eq!(state.open_details(), [0]);
        assert!(!state.dataset_fingerprints.is_empty());
    }

    /// A second entry under a different name.
    ///
    /// `load_example_dataset` cannot be used for this: it always loads under the
    /// name "Examples", and `add_and_activate` replaces a same-named entry in
    /// place rather than appending, so the list would still hold one.
    fn add_second_file(state: &mut AppState) {
        state.apply_loaded_file_bytes("second.smi".to_string(), b"C\nCC\n".to_vec());
        assert_eq!(state.loaded_files.entries().len(), 2, "need two entries");
    }

    #[test]
    fn test_removing_the_active_file_discards_what_belonged_to_it() {
        let mut state = AppState::cpu_only();
        add_second_file(&mut state); // now active
        with_derived_data(&mut state);
        let epoch = state.dataset_epoch();

        state.remove_loaded_file(state.loaded_files.active_index());

        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
        assert!(state.dataset_epoch() > epoch);
        assert!(state.dataset_status.contains("Removed"));
    }

    #[test]
    fn test_removing_another_file_leaves_the_active_one_alone() {
        let mut state = AppState::cpu_only();
        add_second_file(&mut state);
        // The second entry stays active; the first is removed.
        with_derived_data(&mut state);
        let epoch = state.dataset_epoch();

        state.remove_loaded_file(0);

        // The dataset on screen didn't change, so its fingerprints, results and
        // open windows are still about the right molecules. Discarding them here
        // would be the easy over-correction.
        assert!(!state.dataset_fingerprints.is_empty());
        assert!(!state.search_results.is_empty());
        assert_eq!(state.open_details(), [0]);
        assert_eq!(state.dataset_epoch(), epoch);
    }

    #[test]
    fn test_removing_the_only_file_does_nothing() {
        let mut state = AppState::cpu_only();
        with_derived_data(&mut state);
        let epoch = state.dataset_epoch();

        state.remove_loaded_file(0);

        assert_eq!(state.loaded_files.entries().len(), 1);
        assert_eq!(state.dataset_epoch(), epoch);
        assert!(!state.dataset_fingerprints.is_empty());
    }

    #[test]
    fn test_switching_dataset_clears_what_the_operations_reported() {
        let mut state = AppState::cpu_only();
        state.detect_aromaticity_for_dataset();
        state.generate_coordinates_for_dataset();
        assert!(state.aromaticity.has_run());
        assert!(state.coordinates.has_run());

        state.load_example_dataset();

        // "15 of 15 aromatic" against a dataset that has been swapped out is
        // worse than saying nothing.
        assert!(!state.aromaticity.has_run());
        assert!(!state.coordinates.has_run());
        assert!(!state.fingerprints.has_run());
        assert!(!state.search.has_run());
    }

    #[test]
    fn test_generating_coordinates_reports_generated_against_kept() {
        let mut state = AppState::cpu_only();
        let total = state.loaded_files.active_dataset().len();

        // The example set is SMILES, so nothing arrives with coordinates.
        state.generate_coordinates_for_dataset();
        let first = state.coordinates.summary();
        assert!(first.contains(&format!("{} generated", total)), "{}", first);
        assert!(!first.contains("kept"), "{}", first);

        // Second time everything already has them, which is the same branch
        // that keeps coordinates an SDF supplied (#88).
        state.generate_coordinates_for_dataset();
        let second = state.coordinates.summary();
        assert!(second.contains("0 generated"), "{}", second);
        assert!(
            second.contains(&format!("{} kept from file", total)),
            "{}",
            second
        );
    }

    #[test]
    fn test_an_operation_that_has_not_run_says_so() {
        let state = AppState::cpu_only();
        assert!(!state.fingerprints.has_run());
        // An em dash, not an empty header line.
        assert_eq!(state.fingerprints.summary(), "\u{2014}");
    }

    #[test]
    fn test_an_operation_on_an_empty_dataset_fails_in_its_own_section() {
        let mut state = AppState::cpu_only();
        state.loaded_files.active_dataset_mut().molecules.clear();
        state.loaded_files.active_dataset_mut().smiles.clear();
        state.loaded_files.active_dataset_mut().names.clear();

        state.detect_aromaticity_for_dataset();

        assert!(state.aromaticity.failed());
        // And not into the dataset's status line, which belongs to another
        // window now.
        assert!(!state.dataset_status.contains("No dataset"));
    }

    #[test]
    fn test_a_parsed_query_remembers_the_smiles_it_came_from() {
        let mut state = AppState::cpu_only();
        state.parse_query("c1ccccc1", FingerprintParams::default());

        // Whatever draws the parsed molecule needs the string it came from, and
        // shouldn't have to reach into the view that owns the text box for it.
        assert!(state.query_molecule.is_some());
        assert_eq!(state.query_source, "c1ccccc1");
    }

    #[test]
    fn test_a_failed_parse_leaves_no_stale_query_source() {
        let mut state = AppState::cpu_only();
        state.parse_query("c1ccccc1", FingerprintParams::default());
        assert_eq!(state.query_source, "c1ccccc1");

        state.parse_query("not a molecule", FingerprintParams::default());

        // Otherwise the Query section would label an absent structure with the
        // last SMILES that happened to work.
        assert!(state.query_molecule.is_none());
        assert!(state.query_source.is_empty());
        assert!(state.query.failed());
    }

    #[test]
    fn test_the_parsed_smiles_is_trimmed_not_the_raw_box_contents() {
        let mut state = AppState::cpu_only();
        state.parse_query("  c1ccccc1  ", FingerprintParams::default());
        assert_eq!(state.query_source, "c1ccccc1");
    }

    #[test]
    fn test_search_needs_both_a_query_and_a_fingerprinted_dataset() {
        let mut state = AppState::cpu_only();
        assert!(!state.can_search());

        state.query_fingerprint = Some(BitVec::new());
        assert!(!state.can_search(), "nothing to search through yet");

        state.dataset_fingerprints = vec![BitVec::new()];
        assert!(state.can_search());
    }

    #[test]
    fn test_new_results_invalidate_a_view_holding_an_index_into_the_old_ones() {
        let mut state = AppState::cpu_only();
        let epoch = state.results_epoch();

        state.apply_search_result(
            Ok(vec![SearchResult {
                index: 0,
                similarity: 0.5,
            }]),
            1.0,
        );

        assert_eq!(state.search_results.len(), 1);
        assert!(state.search.summary().contains("1 hits"));
        assert!(state.search.summary().contains("1.00ms"));
        assert!(!state.search.failed());
        assert!(state.results_epoch() > epoch);
    }

    #[test]
    fn test_a_failed_search_also_invalidates_the_results() {
        let mut state = AppState::cpu_only();
        state.apply_search_result(
            Ok(vec![SearchResult {
                index: 0,
                similarity: 0.5,
            }]),
            1.0,
        );
        let epoch = state.results_epoch();

        state.apply_search_result(Err(anyhow::anyhow!("boom")), 1.0);

        // The old results are gone either way, so an index into them is stale
        // whether the new search succeeded or not.
        assert!(state.search_results.is_empty());
        assert!(state.query_error.is_some());
        assert!(state.search.failed());
        assert!(state.results_epoch() > epoch);
    }
}
