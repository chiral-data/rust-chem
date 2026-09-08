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
use chem::io::format;
use chem::io::smiles::parse_smiles;
use chem::search::{FingerprintSearch, SearchResult};
use std::cell::RefCell;
use std::collections::HashMap;

#[cfg(target_arch = "wasm32")]
use chem::gpu::{GpuMorganFingerprint, GpuTanimoto};
#[cfg(target_arch = "wasm32")]
use std::rc::Rc;

#[cfg(target_arch = "wasm32")]
type PendingGpuInit = Rc<RefCell<Option<Result<(GpuMorganFingerprint, GpuTanimoto), String>>>>;
#[cfg(target_arch = "wasm32")]
type PendingFileLoad = Rc<RefCell<Vec<(String, Vec<u8>)>>>;

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
    pub convert: OperationOutcome,
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

    /// Bumped whenever new search results land.
    ///
    /// Results, and the indices into them, belong to whichever search produced
    /// them. Rather than every path reaching into every view to clear what it
    /// holds — which is what the old single struct did, in a six-line block
    /// duplicated three times — the paths bump this, and a view holding an
    /// index notices its own state is stale. That also means a view which isn't
    /// being drawn this frame still finds out.
    ///
    /// There was a `dataset_epoch` beside this until #273, for views holding
    /// row indices. Its only reader was the detail window's own layout cache,
    /// and that cache is now [`AppState::layouts`] — shared, and cleared
    /// directly by the paths that invalidate it.
    results_epoch: u64,

    /// Laid-out copies of the active dataset's molecules, keyed by row.
    ///
    /// A view cannot lay a molecule out itself: generating coordinates needs
    /// the dataset mutably and a table is reading it. So the table and the
    /// result rows used to show a dash for a molecule the detail window drew
    /// perfectly well, and each ran its own cache -- which meant running
    /// 2D Coordinates could leave two views showing two different, both valid,
    /// layouts of one molecule (#273).
    ///
    /// Beside the dataset rather than in it, so `has_coords()` keeps meaning
    /// "the file carried a layout" (#270) and this answers the separate
    /// question of what it looks like. `RefCell` because callers hold `&self`
    /// through a row closure; every borrow is taken and dropped inside one
    /// method, never held across drawing.
    layouts: RefCell<HashMap<usize, Molecule>>,

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

    /// The frame loop, so work that finishes off-frame can wake it.
    ///
    /// Everything above that lands in a slot is collected by
    /// `poll_pending_work`, which only runs inside `update()` — and eframe is
    /// reactive, so it only calls `update()` when it paints. Without a repaint
    /// request the result waits for the user to move the mouse (#186). Cloning
    /// is cheap: `egui::Context` is a handle, not the state itself.
    repaint: egui::Context,
}

/// What separates a converted dataset's name from where it came from.
///
/// `add_and_activate` replaces a same-named entry *in place*, so naming a
/// conversion `aromatics.sdf` would silently destroy a loaded file called that.
/// Carrying the provenance instead cannot collide with a filename -- and does
/// collide with itself, so converting twice to one target replaces rather than
/// piling up entries (#275).
const DERIVED_SEPARATOR: &str = " \u{2192} ";

/// What a dataset converted from `source_name` to `target` is called.
fn derived_dataset_name(source_name: &str, target: DatasetFormat) -> String {
    // From the original rather than the chain, so converting A to B to C reads
    // as "A -> C" instead of accumulating every step it went through.
    let origin = source_name
        .split(DERIVED_SEPARATOR)
        .next()
        .unwrap_or(source_name);
    format!("{origin}{DERIVED_SEPARATOR}{}", target.label())
}

/// The file dialog's filters: everything readable first, then one per format.
///
/// Generated from the registry rather than written out. The list used to be two
/// `add_filter` calls naming SMILES and SDF, duplicated across the two `cfg`
/// arms of [`AppState::load_dataset_from_file`] -- so the nine formats v0.8.0
/// added could not be picked at all, and the two copies were free to drift
/// (#266).
///
/// Readable formats only: offering a file the reader would refuse is worse than
/// not offering it. The combined entry comes first because that is the one rfd
/// preselects, and someone opening a `.pdb` should not have to know which entry
/// claims it.
fn molecule_file_filters() -> Vec<(&'static str, Vec<&'static str>)> {
    let readable = || format::all().filter(|f| f.can_read());

    let mut every: Vec<&'static str> = Vec::new();
    for extension in readable().flat_map(|f| f.extensions().iter().copied()) {
        // Two formats may claim one extension -- the registry pins codes as
        // unique but not extensions -- and a repeat in this list would show up
        // in the dialog.
        if !every.contains(&extension) {
            every.push(extension);
        }
    }

    let mut filters = vec![("Molecule files", every)];
    filters.extend(readable().map(|f| (f.label(), f.extensions().to_vec())));
    filters
}

/// A file dialog that offers every format this build can read.
fn molecule_file_dialog() -> rfd::AsyncFileDialog {
    let mut dialog = rfd::AsyncFileDialog::new();
    for (name, extensions) in molecule_file_filters() {
        dialog = dialog.add_filter(name, &extensions);
    }
    dialog
}

/// The name and bytes of a dropped file, whichever half the backend filled in.
///
/// The two backends describe a drop differently and neither fills in the other's
/// fields. `egui-winit` sets `path` and leaves `name` **empty**; eframe's web
/// backend sets `name` and `bytes` and has no path to give. So on native the
/// name has to come from the path — using `DroppedFile::name` there would hand
/// `DatasetFormat::from_filename` an empty string, which resolves to SMILES, and
/// a dropped `.pdb` would parse as SMILES, skip every line, and look like an
/// empty file rather than a bug.
fn dropped_file_contents(file: &egui::DroppedFile) -> Option<(String, Vec<u8>)> {
    if let Some(bytes) = &file.bytes {
        return Some((file.name.clone(), bytes.to_vec()));
    }

    #[cfg(not(target_arch = "wasm32"))]
    if let Some(path) = &file.path {
        let name = path.file_name()?.to_string_lossy().into_owned();
        match std::fs::read(path) {
            Ok(bytes) => return Some((name, bytes)),
            Err(e) => {
                log::error!("Could not read dropped file {name}: {e}");
                return None;
            }
        }
    }

    None
}

/// One line for a whole batch, naming what did not load.
///
/// The counts are what the Files list already shows per entry, so the summary
/// leads with them and then spends its length on what the list *cannot* show: a
/// refused file has no entry, a skipped record has no entry, and a replaced name
/// looks identical to a file that simply loaded.
fn summarise_load(outcomes: &[FileLoad]) -> String {
    // Deduplicated by name, keeping the last: a file that replaced an entry of
    // the same name did not add one, and counting both would report "2 files, 3
    // molecules" for a Files list holding one entry of two.
    let mut surviving: Vec<(&str, usize)> = Vec::new();
    for outcome in outcomes {
        if let FileLoad::Loaded {
            name, molecules, ..
        } = outcome
        {
            match surviving.iter_mut().find(|(seen, _)| *seen == name) {
                Some(entry) => entry.1 = *molecules,
                None => surviving.push((name, *molecules)),
            }
        }
    }
    let molecules: usize = surviving.iter().map(|(_, n)| n).sum();

    let mut summary = match surviving.len() {
        0 => "Loaded nothing".to_string(),
        1 => format!("Loaded {molecules} {}", plural(molecules, "molecule")),
        files => format!(
            "Loaded {files} files, {molecules} {}",
            plural(molecules, "molecule")
        ),
    };

    for outcome in outcomes {
        match outcome {
            FileLoad::Loaded {
                name,
                skipped,
                replaced,
                ..
            } => {
                if *skipped > 0 {
                    summary.push_str(&format!(" \u{b7} {name}: {skipped} skipped"));
                }
                if *replaced {
                    summary.push_str(&format!(" \u{b7} {name}: replaced"));
                }
            }
            FileLoad::Refused { name, reason } => {
                summary.push_str(&format!(" \u{b7} {name}: {reason}"));
            }
        }
    }
    summary
}

fn plural(n: usize, word: &str) -> String {
    if n == 1 {
        word.to_string()
    } else {
        format!("{word}s")
    }
}

/// What became of one file in a load.
///
/// [`AppState::apply_loaded_file_bytes`] used to return `()`, which was enough
/// while a load was one file: it wrote `dataset_status` and the caller had
/// nothing to decide. A batch has to say what happened to *each* file, because
/// the Files list cannot — a refused file leaves no entry in it at all, so
/// mixed with a good one it would vanish without this (#296).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FileLoad {
    Loaded {
        name: String,
        /// Where it landed, so a batch can activate the first one it loaded.
        index: usize,
        molecules: usize,
        skipped: usize,
        /// An entry of this name already existed and was replaced in place.
        /// Silent until now, and much easier to hit when several files arrive
        /// at once: two directories can each hold a `d.smi`.
        replaced: bool,
    },
    Refused {
        name: String,
        reason: &'static str,
    },
}

impl AppState {
    pub fn new(ctx: &egui::Context) -> Self {
        Self::with_engine(ctx, FingerprintSearch::new())
    }

    /// A CPU-only instance, for tests that exercise state transitions and have
    /// no business depending on whether the machine running them has a GPU —
    /// or on GPU init running at all, which does not take kindly to being
    /// driven from parallel tests (#19). Mirrors
    /// [`FingerprintSearch::new_cpu_only`], which exists for the same reason.
    ///
    /// Its context is detached from any window, which is what a test wants: a
    /// repaint request can be observed but never has to be serviced.
    #[cfg(test)]
    pub(crate) fn cpu_only() -> Self {
        Self::with_engine(&egui::Context::default(), FingerprintSearch::new_cpu_only())
    }

    fn with_engine(ctx: &egui::Context, search_engine: FingerprintSearch) -> Self {
        let dataset = MoleculeDataset::example_dataset();
        let dataset_status = format!("Loaded {} example molecules", dataset.len());
        let loaded_files = LoadedFiles::new("Examples".to_string(), dataset, DatasetFormat::SMILES);

        #[cfg(target_arch = "wasm32")]
        let pending_gpu_init = Rc::new(RefCell::new(None));
        #[cfg(target_arch = "wasm32")]
        {
            let slot = pending_gpu_init.clone();
            // `self` does not exist yet, so this one takes the parameter.
            let ctx = ctx.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = FingerprintSearch::try_init_gpu_async().await;
                *slot.borrow_mut() = Some(result);
                ctx.request_repaint();
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
            convert: OperationOutcome::default(),
            query: OperationOutcome::default(),
            search: OperationOutcome::default(),
            display: DisplaySettings::default(),
            open_details: Vec::new(),
            results_epoch: 0,
            layouts: RefCell::new(HashMap::new()),
            dataset_fingerprint_task: Task::new(),
            query_fingerprint_task: Task::new(),
            search_task: Task::new(),
            #[cfg(target_arch = "wasm32")]
            pending_file_load: Rc::new(RefCell::new(Vec::new())),
            #[cfg(target_arch = "wasm32")]
            pending_gpu_init,
            repaint: ctx.clone(),
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

    /// The molecule at `row`, laid out so it can be drawn, or `None` when there
    /// is no such row.
    ///
    /// Laid out once per dataset and kept, so every view that draws this
    /// molecule draws the same picture -- `layout` is not deterministic across
    /// runs, so two views computing their own would disagree (#273).
    ///
    /// A molecule whose file supplied a layout is returned untouched:
    /// `ensure_coords` computes one only where there is none.
    ///
    /// Owned rather than borrowed, so no `RefCell` guard is held while the
    /// caller paints. It costs no more than the caller is about to spend --
    /// `StructureView` rebuilds the depiction every frame regardless.
    pub fn drawable(&self, row: usize) -> Option<Molecule> {
        if let Some(cached) = self.layouts.borrow().get(&row) {
            return Some(cached.clone());
        }

        let source = self
            .loaded_files
            .active_dataset()
            .molecules
            .get(row)?
            .clone();
        let mut prepared = source;
        ensure_coords(&mut prepared);
        self.layouts.borrow_mut().insert(row, prepared.clone());
        Some(prepared)
    }

    /// The molecules changed under their cached layouts, so drop them.
    ///
    /// For an operation that mutates the active dataset *in place* --
    /// coordinates, aromaticity -- as distinct from one that replaces it, which
    /// goes through [`AppState::invalidate_active_dataset`].
    ///
    /// Aromaticity is the one observable today: it changes the molecule, so a
    /// clone cached before perception ran keeps drawing Kekulé bonds. The
    /// coordinates call is the same rule applied consistently rather than a
    /// bug being fixed -- see the note there.
    fn dataset_molecules_changed(&mut self) {
        self.layouts.get_mut().clear();
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
        self.convert = OperationOutcome::default();
        self.search = OperationOutcome::default();
        // Keyed by row index too, and a row means nothing against a different
        // dataset.
        self.layouts.get_mut().clear();
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
            let files = molecule_file_dialog().pick_files().await?;
            let mut picked = Vec::with_capacity(files.len());
            for file in files {
                let name = file.file_name();
                picked.push((name, file.read().await));
            }
            Some(picked)
        });

        if let Some(picked) = picked {
            self.apply_loaded_files(picked);
        }
    }

    // Browsers have only one JS thread, and it also drives the file picker's
    // Promise machinery, so blocking on it would deadlock the tab. Spawn the
    // dialog as a non-blocking task instead and hand its result to
    // `pending_file_load`, polled from `update()` on the next frame.
    #[cfg(target_arch = "wasm32")]
    pub fn load_dataset_from_file(&mut self) {
        let slot = self.pending_file_load.clone();
        let ctx = self.repaint.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let files = molecule_file_dialog().pick_files().await;
            if let Some(files) = files {
                let mut picked = Vec::with_capacity(files.len());
                for file in files {
                    let name = file.file_name();
                    picked.push((name, file.read().await));
                }
                // Extends rather than replaces: the slot used to hold one file,
                // so a second pick before the next frame silently dropped the
                // first (#296).
                slot.borrow_mut().extend(picked);
                // Closing the picker is not itself an input event the canvas
                // sees, so without this the file stays unloaded until the user
                // moves the mouse (#186).
                ctx.request_repaint();
            }
        });
    }

    /// Whether any format this build can read claims this file's extension.
    ///
    /// Only drops ask. A file chosen through the dialog came from a list built
    /// from the registry (#266), so an unrecognised extension there is someone
    /// overriding the filter on purpose and keeps the documented behaviour of
    /// being read as SMILES. A dropped file passed no filter at all, so this is
    /// the only place the answer can come from.
    fn is_readable_extension(name: &str) -> bool {
        let Some(extension) = std::path::Path::new(name).extension() else {
            return false;
        };
        let Some(extension) = extension.to_str() else {
            return false;
        };
        let extension = extension.to_ascii_lowercase();
        format::all()
            .filter(|f| f.can_read())
            .flat_map(|f| f.extensions().iter().copied())
            .any(|known| known.eq_ignore_ascii_case(&extension))
    }

    pub fn apply_loaded_file_bytes(&mut self, name: String, bytes: Vec<u8>) -> FileLoad {
        let content = match String::from_utf8(bytes) {
            Ok(content) => content,
            Err(e) => {
                self.dataset_status = "Failed to load file: not valid UTF-8".to_string();
                log::error!("Dataset load failed: {}", e);
                return FileLoad::Refused {
                    name,
                    reason: "not valid UTF-8",
                };
            }
        };

        let format = DatasetFormat::from_filename(&name);
        let outcome = chem::io::reader::read(&content, format);
        let dataset = MoleculeDataset::from_outcome(&outcome, format);

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

        let molecules = dataset.len();
        let skipped = outcome.skipped.len();
        let replaced = self.loaded_files.names().any(|existing| existing == name);
        self.loaded_files
            .add_and_activate(name.clone(), dataset, format);
        self.invalidate_active_dataset();
        FileLoad::Loaded {
            name,
            index: self.loaded_files.active_index(),
            molecules,
            skipped,
            replaced,
        }
    }

    /// Loads several files as one action, from the dialog or from a drop.
    ///
    /// Two things it does that a loop over
    /// [`AppState::apply_loaded_file_bytes`] would not.
    ///
    /// It writes **one** status for the batch. `dataset_status` is a single
    /// string rendered by a single label, so per-file messages overwrite each
    /// other and the last file wins — which would silently swallow a refusal,
    /// the one outcome that leaves no Files entry to notice afterwards.
    ///
    /// And it activates the **first** file that loaded rather than the last.
    /// `add_and_activate` activates whatever it just added, so without this the
    /// batch would leave you looking at the file you happened to select last.
    /// Going back through [`AppState::activate_loaded_file`] rather than
    /// `LoadedFiles::activate` matters: it runs the same invalidation a click in
    /// the list does, and fingerprints belong to whichever dataset was active
    /// when they were computed.
    pub fn apply_loaded_files(&mut self, files: Vec<(String, Vec<u8>)>) {
        let outcomes: Vec<FileLoad> = files
            .into_iter()
            .map(|(name, bytes)| self.apply_loaded_file_bytes(name, bytes))
            .collect();
        self.finish_batch(outcomes);
    }

    /// Refuses a dropped file whose extension no readable format claims, and
    /// otherwise loads it.
    ///
    /// The refusal is the whole difference from the dialog path. Refusals go
    /// into the same list as the loads, in the order the files arrived, so the
    /// summary reads as one sentence about one gesture -- building a second
    /// summary for them and appending it produced `Loaded 2 files, 3 molecules
    /// \u{b7} Loaded nothing \u{b7} logo.png: ...`, which is two answers to one
    /// question.
    pub fn apply_dropped_files(&mut self, files: Vec<(String, Vec<u8>)>) {
        let outcomes: Vec<FileLoad> = files
            .into_iter()
            .map(|(name, bytes)| {
                if Self::is_readable_extension(&name) {
                    self.apply_loaded_file_bytes(name, bytes)
                } else {
                    FileLoad::Refused {
                        name,
                        reason: "not a format this build reads",
                    }
                }
            })
            .collect();
        self.finish_batch(outcomes);
    }

    /// Activates the first file that loaded, then says what became of the batch.
    ///
    /// The first rather than the last because `add_and_activate` activates
    /// whatever it just added, so a batch would otherwise leave you looking at
    /// whichever file happened to be selected last. Going back through
    /// [`AppState::activate_loaded_file`] rather than `LoadedFiles::activate`
    /// runs the same invalidation a click in the list does -- and it writes its
    /// own status, which is why the summary is written after it rather than
    /// before.
    fn finish_batch(&mut self, outcomes: Vec<FileLoad>) {
        if outcomes.is_empty() {
            return;
        }

        if let Some(FileLoad::Loaded { index, .. }) = outcomes
            .iter()
            .find(|outcome| matches!(outcome, FileLoad::Loaded { .. }))
        {
            self.activate_loaded_file(*index);
        }

        self.dataset_status = summarise_load(&outcomes);
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

        self.dataset_fingerprint_task
            .start(&self.repaint, async move {
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
        let summary = format!("{} of {} aromatic", aromatic, dataset.len());
        // Perception changes the picture, so anything laid out before it ran is
        // now drawing the wrong bonds.
        self.dataset_molecules_changed();
        self.aromaticity = OperationOutcome::ok(summary);
    }

    /// What converting the active dataset to `target` would discard: each
    /// attribute, and how many molecules lose it.
    ///
    /// The command line asks one question -- can the target hold this? -- and
    /// that is not the whole of it. A conversion is a read and a write, and
    /// #257 pinned six pairs that lose an attribute *both* masks claim: CML to
    /// SMILES drops aromaticity and hands back cyclohexane, which a target-only
    /// report cannot see (#261). `chem convert` missed them too until #276
    /// pointed both at one function.
    ///
    /// Deliberately not `format::fidelity`, which folds in what the *target*
    /// manufactures. That answers what the output will contain; this answers
    /// what the input loses.
    ///
    /// First-seen order with a count, matching the CLI's own tracker so the two
    /// can be read side by side -- and, since #276, from the same formula:
    /// `format::kept` is what both call, so they cannot drift again.
    pub fn conversion_losses(&self, target: DatasetFormat) -> Vec<(&'static str, usize)> {
        let source = self.loaded_files.active_format();
        let kept = format::kept(source, target);

        let mut losses: Vec<(&'static str, usize)> = Vec::new();
        for molecule in &self.loaded_files.active_dataset().molecules {
            for attribute in format::held(molecule).difference(kept).names() {
                match losses.iter_mut().find(|(name, _)| *name == attribute) {
                    Some((_, count)) => *count += 1,
                    None => losses.push((attribute, 1)),
                }
            }
        }
        losses
    }

    /// Converts the active dataset and adds the result as a dataset of its own.
    ///
    /// Deliberately a round trip -- written, then read back -- so what the user
    /// looks at is what the conversion actually produced rather than the
    /// molecules it started from. Converting CML to SMILES puts `C1CCCCC1` in
    /// the table where the original had `c1ccccc1`: the drop report predicts
    /// that loss, and this makes it visible (#275).
    ///
    /// The original stays in the file list, so the two can be compared.
    pub fn convert_dataset(&mut self, target: DatasetFormat) -> bool {
        let dataset = self.loaded_files.active_dataset();
        if dataset.is_empty() {
            self.convert = OperationOutcome::failure("No dataset loaded");
            return false;
        }

        let records: Vec<(String, Molecule)> = dataset
            .names
            .iter()
            .cloned()
            .zip(dataset.molecules.iter().cloned())
            .collect();

        let Some(text) = target.write(&records) else {
            // Every registered format writes today, so this is a format that
            // grew a reader and no writer -- the picker filters on `can_write`,
            // making this the belt to that braces.
            self.convert =
                OperationOutcome::failure(format!("{} cannot be written", target.label()));
            return false;
        };

        let losses = self.conversion_losses(target);
        let source_name = self.loaded_files.entries()[self.loaded_files.active_index()]
            .name
            .clone();
        let derived = derived_dataset_name(&source_name, target);

        let outcome = chem::io::reader::read(&text, target);
        let converted = MoleculeDataset::from_outcome(&outcome, target);
        let written = converted.len();

        self.dataset_status = if outcome.skipped.is_empty() {
            format!("Converted {written} molecules to {}", target.label())
        } else {
            // Our own output failing to read back is worth saying out loud
            // rather than logging, the same as it is for a loaded file.
            format!(
                "Converted {written} molecules to {} — {} did not read back",
                target.label(),
                outcome.skipped.len()
            )
        };
        for skipped in &outcome.skipped {
            log::warn!(
                "Conversion to {} produced record {} that did not read back: {}",
                target.label(),
                skipped.position,
                skipped.error
            );
        }

        self.loaded_files
            .add_and_activate(derived, converted, target);
        // Before the outcome, not after: this resets `self.convert` along with
        // everything else the old dataset derived, so writing the summary first
        // would leave the section's header blank after a conversion that worked.
        self.invalidate_active_dataset();

        let summary = if losses.is_empty() {
            format!("{written} converted to {}", target.label())
        } else {
            let named: Vec<String> = losses
                .iter()
                .map(|(attribute, count)| format!("{attribute} ({count})"))
                .collect();
            format!(
                "{written} converted to {} — lost {}",
                target.label(),
                named.join(", ")
            )
        };
        self.convert = OperationOutcome::ok(summary);
        true
    }

    /// The active dataset written in its own format, and a filename for it.
    ///
    /// Separate from converting: this writes whatever dataset is active, so it
    /// works for one that was merely loaded. Changing format is what Convert is
    /// for.
    ///
    /// Returns the text rather than saving it -- the dialog belongs to the
    /// view, as it does for an SVG export, and on web there is no dialog at all,
    /// just a download the browser takes over.
    pub fn export_active_dataset(&self) -> Option<(String, String)> {
        let format = self.loaded_files.active_format();
        let dataset = self.loaded_files.active_dataset();
        if dataset.is_empty() {
            return None;
        }

        let records: Vec<(String, Molecule)> = dataset
            .names
            .iter()
            .cloned()
            .zip(dataset.molecules.iter().cloned())
            .collect();
        let text = format.write(&records)?;
        Some((self.suggested_output_name(format), text))
    }

    /// A filename for the active dataset: its own stem with the format's
    /// extension, and the provenance suffix a converted dataset carries
    /// removed.
    ///
    /// Sanitised for the same reason `draw::svg::suggested_filename` is -- a
    /// name comes from a file's own records and can hold anything, and a slash
    /// would quietly redirect where the file lands.
    fn suggested_output_name(&self, format: DatasetFormat) -> String {
        let entry = &self.loaded_files.entries()[self.loaded_files.active_index()];
        let displayed = entry.name.split(DERIVED_SEPARATOR).next().unwrap_or("");
        let stem: String = displayed
            .rsplit_once('.')
            .map_or(displayed, |(stem, _)| stem)
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        let stem = stem.trim_matches('_');
        let stem = if stem.is_empty() { "molecules" } else { stem };
        format!("{stem}.{}", format.extensions().first().unwrap_or(&"txt"))
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
        // The dataset now holds layouts of its own. Not observable today --
        // `layout` is deterministic, so the cached layout and the one this just
        // computed are identical -- but the rule is that an operation mutating
        // the dataset in place drops the cache, and #110's layout refinement
        // would make the difference real.
        self.dataset_molecules_changed();
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

        self.query_fingerprint_task
            .start(&self.repaint, async move {
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

        self.search_task.start(&self.repaint, async move {
            engine.search_async(&query_fp, &target_fps, top_k).await
        });
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
        let ctx = self.repaint.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let result = FingerprintSearch::try_init_gpu_async().await;
            *slot.borrow_mut() = Some(result);
            ctx.request_repaint();
        });
    }

    /// Collects anything that finished since the last frame.
    ///
    /// Called once per frame, before any view is drawn, so a view never sees a
    /// half-applied result.
    pub fn poll_pending_work(&mut self) {
        #[cfg(target_arch = "wasm32")]
        {
            // Drained whole, so the Files list never shows half a batch.
            let loaded: Vec<(String, Vec<u8>)> =
                self.pending_file_load.borrow_mut().drain(..).collect();
            self.apply_loaded_files(loaded);

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

        // Dropped files are raw input, read through the stored context rather
        // than a parameter -- it is the same one eframe passes to `update`.
        // egui clears them each frame, so this is the one chance to take them.
        let dropped: Vec<egui::DroppedFile> = self.repaint.input(|i| i.raw.dropped_files.clone());
        if !dropped.is_empty() {
            let files: Vec<(String, Vec<u8>)> =
                dropped.iter().filter_map(dropped_file_contents).collect();
            self.apply_dropped_files(files);
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

        state.load_example_dataset();

        // Fingerprints and results belong to whichever dataset was active when
        // they were computed, and a row index means nothing against a different
        // dataset.
        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
        // Views hold indices too, and find out the same way.
    }

    #[test]
    fn test_switching_between_loaded_files_invalidates_the_same_way() {
        let mut state = AppState::cpu_only();
        state.load_example_dataset(); // a second entry to switch between
        with_derived_data(&mut state);

        state.activate_loaded_file(0);

        assert_eq!(state.loaded_files.active_index(), 0);
        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
    }

    #[test]
    fn test_a_failed_load_leaves_the_active_dataset_alone() {
        let mut state = AppState::cpu_only();
        with_derived_data(&mut state);
        let files_before = state.loaded_files.names().count();

        // Not valid UTF-8, so it never reaches a parser.
        state.apply_loaded_file_bytes("broken.smi".to_string(), vec![0xff, 0xfe]);

        assert!(state.dataset_status.contains("Failed to load"));
        // Nothing was replaced, so nothing derived from it is stale.
        assert_eq!(state.loaded_files.names().count(), files_before);
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

        state.remove_loaded_file(state.loaded_files.active_index());

        assert!(state.dataset_fingerprints.is_empty());
        assert!(state.search_results.is_empty());
        assert!(state.open_details().is_empty());
        assert!(state.dataset_status.contains("Removed"));
    }

    #[test]
    fn test_removing_another_file_leaves_the_active_one_alone() {
        let mut state = AppState::cpu_only();
        add_second_file(&mut state);
        // The second entry stays active; the first is removed.
        with_derived_data(&mut state);

        state.remove_loaded_file(0);

        // The dataset on screen didn't change, so its fingerprints, results and
        // open windows are still about the right molecules. Discarding them here
        // would be the easy over-correction.
        assert!(!state.dataset_fingerprints.is_empty());
        assert!(!state.search_results.is_empty());
        assert_eq!(state.open_details(), [0]);
    }

    #[test]
    fn test_removing_the_only_file_does_nothing() {
        let mut state = AppState::cpu_only();
        with_derived_data(&mut state);

        state.remove_loaded_file(0);

        assert_eq!(state.loaded_files.entries().len(), 1);
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

    /// The rounded positions of a drawable molecule, for comparing one layout
    /// against another. `layout` is not deterministic across runs, so equality
    /// here means "the same cached layout", not merely "both laid out".
    fn positions(state: &AppState, row: usize) -> Vec<String> {
        state
            .drawable(row)
            .expect("a row")
            .coords()
            .expect("laid out")
            .iter()
            .map(|p| format!("{:.4},{:.4}", p.x, p.y))
            .collect()
    }

    #[test]
    fn test_the_loss_the_command_line_cannot_see() {
        // The story. PDBQT and SDF both claim aromaticity, so a report asking
        // only `held(m).difference(target.carries())` reports nothing -- the
        // loss is structural, since PDBQT carries no bonds for atom
        // aromaticity to ride on. Only the `pair_loss` term catches it.
        //
        // This used to pin `cml -> smi`, which really did hand back
        // cyclohexane; #261 fixed that at the reader boundary, so the pair that
        // demonstrates the mechanism is now one where the loss is correct
        // behaviour rather than a defect.
        let mut state = AppState::cpu_only();
        let pdbqt = DatasetFormat::PDBQT
            .write(&[(
                "benzene".to_string(),
                parse_smiles("c1ccccc1").expect("valid SMILES"),
            )])
            .expect("PDBQT writes");
        state.apply_loaded_file_bytes("rings.pdbqt".to_string(), pdbqt.into_bytes());
        assert_eq!(state.loaded_files.active_format(), DatasetFormat::PDBQT);

        let losses = state.conversion_losses(DatasetFormat::SDF);
        assert!(
            losses.iter().any(|(name, _)| *name == "aromaticity"),
            "the pair loss is not reported: {losses:?}"
        );
    }

    #[test]
    fn test_an_ordinary_loss_is_still_reported() {
        // The other term must still work: XYZ has no bond block at all, which
        // is what `Carries::BONDS` was added to be able to say (#257).
        let state = AppState::cpu_only();
        let losses = state.conversion_losses(DatasetFormat::XYZ);
        assert!(
            losses.iter().any(|(name, _)| *name == "bonds"),
            "{losses:?}"
        );
    }

    #[test]
    fn test_a_conversion_that_keeps_everything_reports_nothing() {
        // So the report is not merely always non-empty. CXSMILES carries
        // everything SMILES does and more.
        let state = AppState::cpu_only();
        assert_eq!(state.conversion_losses(DatasetFormat::CXSMILES), Vec::new());
    }

    #[test]
    fn test_the_count_is_molecules_not_attributes() {
        let mut state = AppState::cpu_only();
        state.apply_loaded_file_bytes(
            "two.smi".to_string(),
            b"CCO ethanol\nCCN ethylamine\n".to_vec(),
        );

        let losses = state.conversion_losses(DatasetFormat::XYZ);
        let bonds = losses.iter().find(|(name, _)| *name == "bonds");
        assert_eq!(bonds, Some(&("bonds", 2)), "{losses:?}");
    }

    /// A dataset read from CML, whose benzene is aromatic on the way in.
    fn cml_benzene(state: &mut AppState) {
        let cml = DatasetFormat::CML
            .write(&[(
                "benzene".to_string(),
                parse_smiles("c1ccccc1").expect("valid SMILES"),
            )])
            .expect("CML writes");
        state.apply_loaded_file_bytes("rings.cml".to_string(), cml.into_bytes());
    }

    #[test]
    fn test_the_conversion_is_visible_in_what_it_produces() {
        // A report predicting a loss is weaker than a dataset that shows it, so
        // this asserts on the molecules rather than the flags (#275).
        //
        // It used to convert CML to SMILES and assert the *loss* was visible --
        // `C1CCCCC1` where the original had `c1ccccc1`. Since #261 that round
        // trip is lossless, so what is asserted is the other half of the same
        // property: the converted dataset is what the conversion actually
        // produced, and here that means the aromaticity survived.
        let mut state = AppState::cpu_only();
        cml_benzene(&mut state);
        assert!(
            format::held(&state.loaded_files.active_dataset().molecules[0])
                .contains(format::Carries::AROMATICITY),
            "the fixture must start aromatic, or this tests nothing"
        );

        assert!(state.convert_dataset(DatasetFormat::SMILES));

        let converted = state.loaded_files.active_dataset();
        assert_eq!(converted.len(), 1);
        assert!(
            format::held(&converted.molecules[0]).contains(format::Carries::AROMATICITY),
            "aromaticity was lost: CML states it in the bond order and the \
             SMILES writer reads the atom flag, which is what #261 reconciled"
        );
        // The column the user reads, and the whole of #261 in one assertion.
        assert_eq!(converted.smiles[0], "c1ccccc1");
    }

    #[test]
    fn test_the_original_dataset_survives_the_conversion() {
        // The point is comparing them, so the one converted from has to still
        // be there and switchable.
        let mut state = AppState::cpu_only();
        cml_benzene(&mut state);
        let before = state.loaded_files.entries().len();

        assert!(state.convert_dataset(DatasetFormat::SMILES));

        assert_eq!(state.loaded_files.entries().len(), before + 1);
        assert!(
            state.loaded_files.names().any(|n| n == "rings.cml"),
            "the source dataset went away"
        );
        assert_eq!(state.loaded_files.active_format(), DatasetFormat::SMILES);
    }

    #[test]
    fn test_converting_twice_replaces_rather_than_piling_up() {
        let mut state = AppState::cpu_only();
        cml_benzene(&mut state);

        assert!(state.convert_dataset(DatasetFormat::SMILES));
        let after_one = state.loaded_files.entries().len();

        // Back to the CML entry: index 0 is the example dataset every session
        // starts on, and converting *that* would add a differently-named entry.
        let source = state
            .loaded_files
            .names()
            .position(|n| n == "rings.cml")
            .expect("the cml entry");
        state.activate_loaded_file(source);
        assert!(state.convert_dataset(DatasetFormat::SMILES));

        assert_eq!(state.loaded_files.entries().len(), after_one);
    }

    #[test]
    fn test_converting_does_not_clobber_a_loaded_file_of_that_name() {
        // `add_and_activate` replaces a same-named entry *in place*, so naming
        // the result `two.sdf` would silently destroy a file the user loaded
        // under that name. The derived name carries its provenance instead.
        let mut state = AppState::cpu_only();
        state.apply_loaded_file_bytes("two.smi".to_string(), b"CCO a\nCCN b\n".to_vec());
        let sdf = DatasetFormat::SDF
            .write(&[("mine".to_string(), parse_smiles("C").expect("valid"))])
            .expect("SDF writes");
        state.apply_loaded_file_bytes("two.sdf".to_string(), sdf.into_bytes());
        let smi = state
            .loaded_files
            .names()
            .position(|n| n == "two.smi")
            .expect("the smi entry");
        state.activate_loaded_file(smi);

        assert!(state.convert_dataset(DatasetFormat::SDF));

        let theirs = state
            .loaded_files
            .entries()
            .iter()
            .find(|e| e.name == "two.sdf")
            .expect("the loaded file is still there");
        assert_eq!(theirs.dataset.len(), 1, "their file was overwritten");
    }

    #[test]
    fn test_the_outcome_survives_the_invalidation_that_precedes_it() {
        // Adding a dataset invalidates what the old one derived, and that
        // resets `convert` along with the rest. Writing the summary first would
        // leave the section's header blank after a conversion that worked.
        let mut state = AppState::cpu_only();
        assert!(state.convert_dataset(DatasetFormat::XYZ));

        assert!(state.convert.has_run(), "the summary was wiped");
        assert!(!state.convert.failed());
        assert!(
            state.convert.summary().contains("bonds"),
            "{}",
            state.convert.summary()
        );
    }

    #[test]
    fn test_exporting_writes_the_active_dataset_named_for_its_format() {
        let mut state = AppState::cpu_only();
        state.apply_loaded_file_bytes("two.smi".to_string(), b"CCO a\nCCN b\n".to_vec());

        assert!(state.convert_dataset(DatasetFormat::SDF));
        let (name, text) = state.export_active_dataset().expect("SDF writes");

        // From the format, not from the display name the conversion gave it.
        assert_eq!(name, "two.sdf");
        assert_eq!(text.matches("$$$$").count(), 2, "both records written");
    }

    #[test]
    fn test_a_filename_from_a_hostile_dataset_name_is_still_a_filename() {
        // Names come from a file's own records; a slash would quietly redirect
        // where the file lands.
        let mut state = AppState::cpu_only();
        state.apply_loaded_file_bytes("../../etc/passwd.smi".to_string(), b"C methane\n".to_vec());

        let (name, _) = state.export_active_dataset().expect("writes");
        assert!(!name.contains('/'), "{name}");
        assert!(name.ends_with(".smi"), "{name}");
    }

    #[test]
    fn test_a_molecule_with_no_layout_is_still_drawable() {
        // The bug: the table and the result rows showed a dash for a molecule
        // the detail window drew perfectly well (#273). The example dataset is
        // SMILES, so nothing in it carries a layout.
        let state = AppState::cpu_only();
        assert!(
            !state.loaded_files.active_dataset().molecules[0].has_coords(),
            "the fixture must start without a layout, or this tests nothing"
        );

        let drawable = state.drawable(0).expect("row 0 exists");
        assert!(drawable.has_coords());

        // Beside the dataset, not in it: `has_coords()` keeps meaning "the file
        // carried a layout" (#270), so 2D Coordinates still has work to report.
        assert!(!state.loaded_files.active_dataset().molecules[0].has_coords());
    }

    #[test]
    fn test_the_layout_distinguishes_the_atoms() {
        // Presence is not enough -- #270 shipped a layout that satisfied
        // `has_coords()` with every atom stacked on one point, undrawable.
        let state = AppState::cpu_only();
        let drawn = positions(&state, 0);
        let mut distinct = drawn.clone();
        distinct.sort();
        distinct.dedup();
        assert_eq!(distinct.len(), drawn.len(), "atoms share a position");
    }

    #[test]
    fn test_every_view_asking_for_a_row_gets_the_same_layout() {
        // The subject of #273. Two views laying out independently would each
        // get a valid but different picture of one molecule.
        let state = AppState::cpu_only();
        assert_eq!(positions(&state, 0), positions(&state, 0));
    }

    #[test]
    fn test_a_new_dataset_drops_the_cached_layouts() {
        // A row index means nothing against a different dataset. Without the
        // clear, row 0 would keep answering with the previous dataset's
        // molecule.
        let mut state = AppState::cpu_only();
        let _ = state.drawable(0);

        state.apply_loaded_file_bytes(
            "one.smi".to_string(),
            b"c1ccccc1 benzene
"
            .to_vec(),
        );

        let drawable = state.drawable(0).expect("row 0 of the new dataset");
        assert_eq!(
            drawable.num_atoms(),
            state.loaded_files.active_dataset().molecules[0].num_atoms(),
            "a stale layout from the previous dataset"
        );
    }

    #[test]
    fn test_detecting_aromaticity_drops_the_cached_layouts() {
        // Perception changes the picture, so a layout cached before it ran
        // draws the wrong bonds. Nothing else pins this, which is exactly why
        // it is the invalidation that would be forgotten.
        let mut state = AppState::cpu_only();
        state.apply_loaded_file_bytes(
            "ring.smi".to_string(),
            b"C1=CC=CC=C1 benzene
"
            .to_vec(),
        );
        let drawn = state.drawable(0).expect("a row");
        assert!(!drawn.atoms().iter().any(|a| a.is_aromatic()));

        state.detect_aromaticity_for_dataset();

        let redrawn = state.drawable(0).expect("a row");
        assert!(
            redrawn.atoms().iter().any(|a| a.is_aromatic()),
            "the cache is still handing out the pre-perception molecule"
        );
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
        // Replaced rather than cleared vector by vector: `new()` is the one
        // constructor, so a parallel vector added later cannot be left behind
        // here still holding rows.
        *state.loaded_files.active_dataset_mut() = MoleculeDataset::new();

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

    #[test]
    fn test_the_file_dialog_offers_every_readable_format() {
        // Asserted against the registry rather than a count, so registering a
        // twelfth format cannot silently go un-offered -- which is exactly what
        // happened to the nine v0.8.0 added (#266).
        let filters = molecule_file_filters();
        let (combined, every) = &filters[0];
        assert_eq!(*combined, "Molecule files");

        for expected in format::all().filter(|f| f.can_read()) {
            assert!(
                filters.iter().any(|(name, _)| *name == expected.label()),
                "{} has no filter entry",
                expected.label()
            );
            for extension in expected.extensions() {
                assert!(
                    every.contains(extension),
                    "{extension} missing from the combined filter"
                );
            }
        }
    }

    #[test]
    fn test_the_dialog_offers_nothing_it_cannot_read() {
        // Offering a file the reader would refuse is worse than not offering
        // it: the picker would accept it and the load would fail afterwards.
        let filters = molecule_file_filters();
        for format in format::all().filter(|f| !f.can_read()) {
            assert!(
                !filters.iter().any(|(name, _)| *name == format.label()),
                "{} cannot be read but is offered",
                format.label()
            );
        }
    }

    #[test]
    fn test_the_combined_filter_lists_each_extension_once() {
        // The registry pins codes as unique but not extensions, so two formats
        // may claim one and a repeat would reach the dialog.
        let filters = molecule_file_filters();
        let mut seen = filters[0].1.clone();
        let before = seen.len();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(
            before,
            seen.len(),
            "the combined filter repeats an extension"
        );
    }

    #[test]
    fn test_a_loaded_mol2_shows_the_molecule_where_it_used_to_show_its_format() {
        // Through the whole load path rather than `from_outcome`, because that
        // is where a user meets it (#283). Mol2 and CML are the two formats
        // the issue was filed about.
        let mut state = AppState::cpu_only();
        let mol2 = DatasetFormat::MOL2
            .write(&[(
                "benzene".to_string(),
                parse_smiles("c1ccccc1").expect("valid SMILES"),
            )])
            .expect("Mol2 writes");
        state.apply_loaded_file_bytes("rings.mol2".to_string(), mol2.into_bytes());

        let dataset = state.loaded_files.active_dataset();
        assert_eq!(dataset.smiles[0], "c1ccccc1");
        assert!(
            dataset.generated[0],
            "the app wrote this string and the views say so"
        );

        cml_benzene(&mut state);
        let dataset = state.loaded_files.active_dataset();
        assert_eq!(dataset.smiles[0], "c1ccccc1");
        assert!(dataset.generated[0]);
    }

    #[test]
    fn test_a_structure_file_loads_and_is_labelled_by_its_own_format() {
        // The `(SDF)` placeholder was applied to every molecule that arrived
        // without a SMILES string, so a PDB's rows claimed to be SDF records
        // (#266). PDB is the right fixture precisely because it has no SMILES.
        let mut state = AppState::cpu_only();
        let pdb = "\
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O
ATOM      2  H1  HOH A   1       0.759   0.000   0.504  1.00 20.00           H
CONECT    1    2
END
";
        state.apply_loaded_file_bytes("water.pdb".to_string(), pdb.as_bytes().to_vec());

        let dataset = state.loaded_files.active_dataset();
        assert_eq!(dataset.len(), 1);
        // Still the placeholder after #283, and deliberately: the fixture has
        // a bond, but PDB's `CONECT` is adjacency with no bond order, so a
        // SMILES written from it would be the right topology and the wrong
        // molecule. Asserting the bond keeps this from passing for the wrong
        // reason -- under the `num_bonds() > 0` predicate #283 proposed, this
        // row would read `[O][H]`.
        assert!(dataset.molecules[0].num_bonds() > 0);
        assert_eq!(dataset.smiles[0], "(PDB)");
        assert!(!dataset.generated[0]);
        assert!(
            !dataset.smiles[0].contains("SDF"),
            "a PDB must not describe itself as an SDF"
        );
        assert!(
            state.dataset_status.contains("PDB"),
            "{}",
            state.dataset_status
        );
    }
    /// Two atoms and a CONECT, so it is a real PDB rather than something the
    /// reader would reject for having no atoms (#292).
    const WATER_PDB: &str = "\
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O
ATOM      2  H1  HOH A   1       0.759   0.000   0.504  1.00 20.00           H
CONECT    1    2
END
";

    #[test]
    fn test_a_dropped_file_is_named_by_its_path_when_that_is_all_the_backend_gave() {
        // The native half, and the one that fails silently. `egui-winit` sets
        // `path` and leaves `name` empty, so reading `name` here would hand
        // `from_filename` an empty string -- which resolves to SMILES, and a
        // dropped PDB would then skip every line and look like an empty file
        // rather than a bug.
        let dir = std::env::temp_dir().join("chem-app-drop-native");
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("water.pdb");
        std::fs::write(&path, WATER_PDB).expect("fixture");

        let dropped = egui::DroppedFile {
            path: Some(path.clone()),
            ..Default::default()
        };
        let (name, bytes) = dropped_file_contents(&dropped).expect("reads from the path");

        assert_eq!(name, "water.pdb");
        assert_eq!(DatasetFormat::from_filename(&name), DatasetFormat::PDB);
        assert_eq!(bytes, WATER_PDB.as_bytes());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_a_dropped_file_uses_the_bytes_when_the_backend_gave_those_instead() {
        // The web half: eframe reads the file itself and hands over name and
        // bytes with no path, so nothing may touch the filesystem here.
        let dropped = egui::DroppedFile {
            name: "water.pdb".to_string(),
            bytes: Some(WATER_PDB.as_bytes().to_vec().into()),
            ..Default::default()
        };
        let (name, bytes) = dropped_file_contents(&dropped).expect("uses the bytes");

        assert_eq!(name, "water.pdb");
        assert_eq!(bytes, WATER_PDB.as_bytes());
    }

    #[test]
    fn test_a_batch_activates_the_first_file_it_loaded() {
        // `add_and_activate` activates whatever it just added, so without the
        // batch stepping back the user would be left looking at whichever file
        // happened to be selected last.
        let mut state = AppState::cpu_only();
        let before = state.loaded_files.entries().len();

        state.apply_loaded_files(vec![
            (
                "first.smi".to_string(),
                b"CCO
"
                .to_vec(),
            ),
            (
                "second.smi".to_string(),
                b"CC
CCC
"
                .to_vec(),
            ),
            (
                "third.smi".to_string(),
                b"C
"
                .to_vec(),
            ),
        ]);

        let names: Vec<&str> = state.loaded_files.names().collect();
        assert_eq!(names.len(), before + 3);
        assert!(names.ends_with(&["first.smi", "second.smi", "third.smi"]));
        assert_eq!(
            state.loaded_files.entries()[state.loaded_files.active_index()].name,
            "first.smi"
        );
        assert!(
            state.dataset_status.contains("Loaded 3 files"),
            "{}",
            state.dataset_status
        );
    }

    #[test]
    fn test_a_refused_file_is_named_and_the_ones_beside_it_still_load() {
        // A drop passes no filter, so this is the only place a user is told
        // that a file is not one this build reads. The status is the only
        // place it can be said: a refused file leaves no Files entry.
        let mut state = AppState::cpu_only();
        let before = state.loaded_files.entries().len();

        state.apply_dropped_files(vec![
            ("logo.png".to_string(), vec![0x89, b'P', b'N', b'G']),
            (
                "good.smi".to_string(),
                b"CCO
"
                .to_vec(),
            ),
        ]);

        assert_eq!(state.loaded_files.entries().len(), before + 1);
        assert!(
            state.dataset_status.contains("logo.png")
                && state
                    .dataset_status
                    .contains("not a format this build reads"),
            "{}",
            state.dataset_status
        );
        assert!(
            state.dataset_status.contains("Loaded"),
            "{}",
            state.dataset_status
        );
    }

    #[test]
    fn test_a_file_that_is_not_utf8_is_still_reported_beside_one_that_loaded() {
        // The regression the checklist rests on: "a binary file is refused with
        // 'not valid UTF-8' and the current dataset is left alone". One status
        // line means a naive loop would let the good file overwrite that, and
        // the refusal leaves no entry to notice afterwards. `.smi` so it gets
        // past the extension check and fails where it is meant to.
        let mut state = AppState::cpu_only();

        state.apply_dropped_files(vec![
            ("broken.smi".to_string(), vec![0xff, 0xfe]),
            (
                "good.smi".to_string(),
                b"CCO
"
                .to_vec(),
            ),
        ]);

        assert!(
            state.dataset_status.contains("broken.smi")
                && state.dataset_status.contains("not valid UTF-8"),
            "{}",
            state.dataset_status
        );
    }

    #[test]
    fn test_two_files_of_one_name_collapse_and_the_status_says_so() {
        // `add_and_activate` replaces a same-named entry in place, which is
        // right for reloading a file and surprising when two directories each
        // hold a `d.smi` -- much easier to hit now that a drop can carry both.
        let mut state = AppState::cpu_only();
        let before = state.loaded_files.entries().len();

        state.apply_loaded_files(vec![
            (
                "d.smi".to_string(),
                b"CCO
"
                .to_vec(),
            ),
            (
                "d.smi".to_string(),
                b"CC
CCC
"
                .to_vec(),
            ),
        ]);

        assert_eq!(state.loaded_files.entries().len(), before + 1);
        assert!(
            state.dataset_status.contains("d.smi: replaced"),
            "{}",
            state.dataset_status
        );
    }

    #[test]
    fn test_only_extensions_a_readable_format_claims_are_accepted_from_a_drop() {
        for name in ["a.smi", "a.pdb", "a.mol2", "a.CIF", "a.json"] {
            assert!(
                AppState::is_readable_extension(name),
                "{name} should be accepted"
            );
        }
        for name in ["logo.png", "notes", "report.docx", "archive.tar.gz"] {
            assert!(
                !AppState::is_readable_extension(name),
                "{name} should be refused"
            );
        }
    }
    #[test]
    fn test_the_status_says_what_became_of_every_file_in_the_batch() {
        // Exact strings, not `contains`. The first version of this feature
        // built a second summary for the refusals and appended it, producing
        // `Loaded 2 files, 3 molecules \u{b7} Loaded nothing \u{b7} logo.png: ...`
        // -- two answers to one question -- and a `contains("logo.png")`
        // assertion passed the whole way through.
        /// A batch to drop, and the one line it should produce.
        type Case = (Vec<(String, Vec<u8>)>, &'static str);

        let cases: Vec<Case> = vec![
            (
                vec![("a.smi".into(), b"CCO\n".to_vec())],
                "Loaded 1 molecule",
            ),
            (
                vec![
                    ("a.smi".into(), b"CCO\n".to_vec()),
                    ("b.smi".into(), b"CC\nC\n".to_vec()),
                ],
                "Loaded 2 files, 3 molecules",
            ),
            (
                vec![
                    ("a.smi".into(), b"CCO\n".to_vec()),
                    ("logo.png".into(), vec![0x89]),
                ],
                "Loaded 1 molecule \u{b7} logo.png: not a format this build reads",
            ),
            (
                vec![("logo.png".into(), vec![0x89])],
                "Loaded nothing \u{b7} logo.png: not a format this build reads",
            ),
            (
                vec![
                    ("a.smi".into(), b"CCO\n".to_vec()),
                    ("broken.smi".into(), vec![0xff, 0xfe]),
                ],
                "Loaded 1 molecule \u{b7} broken.smi: not valid UTF-8",
            ),
            (
                // Two entries went in and one came out, so the count is the
                // survivor's rather than the sum -- "2 files, 3 molecules"
                // would describe a Files list that does not exist.
                vec![
                    ("d.smi".into(), b"CCO\n".to_vec()),
                    ("d.smi".into(), b"CC\nC\n".to_vec()),
                ],
                "Loaded 2 molecules \u{b7} d.smi: replaced",
            ),
            (
                vec![("s.smi".into(), b"CCO\nnot a molecule!!\n".to_vec())],
                "Loaded 1 molecule \u{b7} s.smi: 1 skipped",
            ),
        ];

        for (files, expected) in cases {
            let mut state = AppState::cpu_only();
            state.apply_dropped_files(files);
            assert_eq!(state.dataset_status, expected);
        }
    }
}
