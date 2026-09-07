//! Floating windows showing molecules in full, one per molecule.
//!
//! A window rather than an inline expand below the table: the detail view (info
//! grid + atom list + bond list) can be taller than the window holding the
//! table, with no way to scroll to the rest of it. A window is independently
//! movable and resizable, and scrolls on its own.
//!
//! Several at once, because comparing two molecules is what a similarity search
//! exists to set up, and one window at a time meant clicking back and forth and
//! holding the first structure in your head.
//!
//! Not in the [`crate::windows::WindowRegistry`]: these aren't toggleable
//! singletons with a place in the View menu. They open by clicking a table row,
//! close with their own button, and there are between zero and
//! [`crate::state::MAX_OPEN_DETAILS`] of them.

use crate::molecule_view::{show_atom_list, show_bond_list, show_molecule_info};
use crate::state::AppState;
use crate::structure_view::structure_panel_with_options;
use crate::svg::save_svg;
use chem::draw::structure::StructureTheme;
use chem::draw::svg::{structure_to_svg, suggested_filename};

/// Size of an exported SVG, in points.
///
/// Larger than the on-screen panel: an export is destined for a document or a
/// slide rather than a 220px window, and the geometry is proportional so a
/// bigger viewport simply gives the labels more room.
const EXPORT_SIZE: (f32, f32) = (360.0, 300.0);

/// Offset between successive windows, so the second one to open is visible
/// rather than exactly beneath the first.
const CASCADE_STEP: f32 = 28.0;

/// Where the first window opens, as a fraction of the workspace. Toward the
/// middle, since the default layout leaves canvas there.
const CASCADE_ORIGIN: (f32, f32) = (0.28, 0.16);

#[derive(Default)]
pub struct DetailView {}

impl DetailView {
    pub fn show(&mut self, ctx: &egui::Context, state: &mut AppState) {
        // Cloned so the windows can borrow `state` mutably to close themselves.
        let rows: Vec<usize> = state.open_details().to_vec();
        if rows.is_empty() {
            return;
        }

        let options = state.display.structure;
        let workspace = ctx.available_rect();
        // Caps a window to the viewport height, so it can never grow past the
        // screen on a short monitor — content beyond that scrolls inside
        // instead of the window growing off-screen with no way to reach it.
        let max_height = (workspace.height() - 80.0).max(200.0);

        let mut to_close = Vec::new();

        for (n, &row) in rows.iter().enumerate() {
            let dataset = state.loaded_files.active_dataset();
            let name = dataset.names.get(row).cloned();
            let smiles = dataset.smiles.get(row).cloned();

            // Laid out once per dataset and shared with the table and the
            // result rows, so all three draw the same picture (#273). `None` is
            // the row having gone away under us; forget it rather than drawing
            // an empty window.
            let (Some(mol), Some(name), Some(smiles)) = (state.drawable(row), name, smiles) else {
                to_close.push(row);
                continue;
            };

            let offset = CASCADE_STEP * n as f32;
            let default_pos = workspace.min
                + egui::vec2(
                    workspace.width() * CASCADE_ORIGIN.0 + offset,
                    workspace.height() * CASCADE_ORIGIN.1 + offset,
                );

            let mut open = true;
            egui::Window::new(format!("Molecule: {}", name))
                // Keyed by row, so each molecule keeps its own position rather
                // than every window sharing one.
                .id(egui::Id::new(("molecule_detail", row)))
                .open(&mut open)
                .resizable(true)
                .collapsible(true)
                .constrain(true)
                .default_pos(default_pos)
                .default_height(440.0)
                .max_height(max_height)
                .show(ctx, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        // Beside the structure it exports, rather than in a menu
                        // that would have to guess which of up to
                        // MAX_OPEN_DETAILS molecules was meant.
                        if ui
                            .button("Export SVG")
                            .on_hover_text("Save this structure as an SVG file")
                            .clicked()
                        {
                            // The theme the file gets is decided here rather
                            // than read from the live one: an SVG bound for a
                            // light document should not carry a dark palette
                            // because that is what the app happened to be
                            // showing.
                            let svg = structure_to_svg(
                                &mol,
                                egui::vec2(EXPORT_SIZE.0, EXPORT_SIZE.1),
                                &options,
                                &StructureTheme::light(),
                            );
                            save_svg(&suggested_filename(&name), &svg);
                        }

                        structure_panel_with_options(ui, &mol, 220.0, options);
                        show_molecule_info(ui, &mol, &smiles, &name);
                        show_atom_list(ui, &mol);
                        show_bond_list(ui, &mol);
                    });
                });

            if !open {
                to_close.push(row);
            }
        }

        for row in to_close {
            state.close_detail(row);
        }
    }
}
