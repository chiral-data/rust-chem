use crate::dataset::MoleculeDataset;
use crate::fingerprint_view::{fingerprint_compact, fingerprint_full};
use crate::molecule_view::{molecule_compact, show_atom_list, show_bond_list, show_molecule_info};
use crate::search::{FingerprintSearch, SearchResult};
use bitvec::prelude::BitVec;
use chemcore::molecule::Molecule;
use chemio::smiles::parse_smiles;
use egui::{Color32, RichText};
use std::time::{Duration, Instant};

/// How long the query SMILES box must sit idle before we run fingerprint
/// generation, so a blocking GPU dispatch doesn't fire on every keystroke.
const QUERY_DEBOUNCE: Duration = Duration::from_millis(300);

pub struct ChemFpDemoApp {
    dataset: MoleculeDataset,
    dataset_fingerprints: Vec<BitVec>,
    dataset_status: String,
    search_engine: FingerprintSearch,
    query_smiles: String,
    query_molecule: Option<Molecule>,
    query_fingerprint: Option<BitVec>,
    query_error: Option<String>,
    search_results: Vec<SearchResult>,
    selected_result: Option<usize>,
    fp_radius: u32,
    fp_size: u32,
    top_k: usize,
    last_search_time: Option<f64>,
    last_fp_gen_time: Option<f64>,
    fps_counter: f64,
    frame_times: Vec<f64>,
    query_dirty_since: Option<Instant>,
}

impl ChemFpDemoApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let dataset = MoleculeDataset::example_dataset().unwrap_or_default();
        let dataset_status = format!("Loaded {} example molecules", dataset.len());

        Self {
            dataset,
            dataset_fingerprints: Vec::new(),
            dataset_status,
            search_engine: FingerprintSearch::new(),
            query_smiles: String::from("c1ccccc1"),
            query_molecule: None,
            query_fingerprint: None,
            query_error: None,
            search_results: Vec::new(),
            selected_result: None,
            fp_radius: 2,
            fp_size: 2048,
            top_k: 10,
            last_search_time: None,
            last_fp_gen_time: None,
            fps_counter: 0.0,
            frame_times: Vec::with_capacity(60),
            query_dirty_since: None,
        }
    }

    fn load_dataset_from_file(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("SMILES", &["smi", "smiles", "txt"])
            .pick_file()
        {
            match MoleculeDataset::load_from_smiles_file(&path) {
                Ok(dataset) => {
                    self.dataset_status = format!("Loaded {} molecules from file", dataset.len());
                    self.dataset = dataset;
                    self.dataset_fingerprints.clear();
                    self.search_engine.invalidate_target_dataset();
                    self.search_results.clear();
                    self.selected_result = None;
                    log::info!("Dataset loaded successfully");
                }
                Err(e) => {
                    self.dataset_status = format!("Failed to load file: {}", e);
                    log::error!("Dataset load failed: {}", e);
                }
            }
        }
    }

    fn load_example_dataset(&mut self) {
        match MoleculeDataset::example_dataset() {
            Ok(dataset) => {
                self.dataset_status = format!("Loaded {} example molecules", dataset.len());
                self.dataset = dataset;
                self.dataset_fingerprints.clear();
                self.search_engine.invalidate_target_dataset();
                self.search_results.clear();
                self.selected_result = None;
            }
            Err(e) => {
                self.dataset_status = format!("Failed to load examples: {}", e);
            }
        }
    }

    fn precompute_dataset_fingerprints(&mut self) {
        if self.dataset.is_empty() {
            self.dataset_status = "No dataset loaded".to_string();
            return;
        }

        let start = Instant::now();

        match self.search_engine.generate_fingerprints_batch(
            &self.dataset.molecules,
            self.fp_radius,
            self.fp_size,
        ) {
            Ok(fps) => {
                self.dataset_fingerprints = fps;
                if let Err(e) = self
                    .search_engine
                    .set_target_dataset(&self.dataset_fingerprints)
                {
                    log::warn!("Failed to upload dataset to GPU: {}", e);
                }
                let elapsed = start.elapsed().as_secs_f64();
                self.dataset_status = format!(
                    "Computed {} fingerprints in {:.2}ms ({} mode)",
                    self.dataset_fingerprints.len(),
                    elapsed * 1000.0,
                    if self.search_engine.is_using_gpu() {
                        "GPU"
                    } else {
                        "CPU"
                    }
                );
                log::info!("Fingerprints computed in {:.2}ms", elapsed * 1000.0);
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

                let start = Instant::now();
                match self
                    .search_engine
                    .generate_fingerprint(&mol, self.fp_radius, self.fp_size)
                {
                    Ok(fp) => {
                        self.query_fingerprint = Some(fp);
                        self.last_fp_gen_time = Some(start.elapsed().as_secs_f64() * 1000.0);
                    }
                    Err(e) => {
                        self.query_error = Some(format!("Fingerprint generation failed: {}", e));
                        self.query_fingerprint = None;
                    }
                }
            }
            Err(e) => {
                self.query_error = Some(format!("Invalid SMILES: {}", e));
                self.query_molecule = None;
                self.query_fingerprint = None;
            }
        }
    }

    fn run_search(&mut self) {
        if self.dataset_fingerprints.is_empty() {
            self.query_error = Some("Please compute dataset fingerprints first".to_string());
            return;
        }

        if let Some(ref query_fp) = self.query_fingerprint {
            let start = Instant::now();

            match self
                .search_engine
                .search(query_fp, &self.dataset_fingerprints, self.top_k)
            {
                Ok(results) => {
                    self.search_results = results;
                    self.last_search_time = Some(start.elapsed().as_secs_f64() * 1000.0);
                    self.selected_result = None;
                }
                Err(e) => {
                    self.query_error = Some(format!("Search failed: {}", e));
                    self.search_results.clear();
                }
            }
        }
    }

    fn top_panel(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("🧪 ChemFP Demo - Molecular Fingerprint Search");

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(format!("FPS: {:.0}", self.fps_counter));
                    ui.separator();

                    let mode_text = if self.search_engine.is_using_gpu() {
                        RichText::new("🚀 GPU").color(Color32::from_rgb(50, 200, 50))
                    } else {
                        RichText::new("💻 CPU").color(Color32::from_rgb(200, 200, 50))
                    };
                    ui.label(mode_text);
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

                ui.separator();

                if !self.dataset.is_empty() {
                    ui.label(format!("Molecules in dataset: {}", self.dataset.len()));
                    ui.label(format!(
                        "Fingerprints computed: {}",
                        self.dataset_fingerprints.len()
                    ));

                    ui.separator();

                    egui::ScrollArea::vertical().show(ui, |ui| {
                        for (i, (smiles, name)) in self
                            .dataset
                            .smiles
                            .iter()
                            .zip(self.dataset.names.iter())
                            .enumerate()
                            .take(20)
                        {
                            ui.horizontal(|ui| {
                                ui.label(format!("{}.", i + 1));
                                ui.label(RichText::new(smiles).code().small());
                                ui.label(RichText::new(name).small());
                            });
                        }
                        if self.dataset.len() > 20 {
                            ui.label(format!("... and {} more", self.dataset.len() - 20));
                        }
                    });
                }
            });
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
                        ui.label(format!("Generated in {:.2}ms", time));
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
                    ui.label(format!("Search completed in {:.2}ms", time));
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
                egui::ScrollArea::vertical().show(ui, |ui| {
                    for (rank, result) in self.search_results.iter().enumerate() {
                        let idx = result.index;
                        let mol = &self.dataset.molecules[idx];
                        let smiles = &self.dataset.smiles[idx];
                        let name = &self.dataset.names[idx];

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
                                egui::Stroke::new(2.0, Color32::from_rgb(70, 130, 180)),
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

                                    if let Some(ref fp) = self.dataset_fingerprints.get(idx) {
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
    }
}
