//! A floating window showing one molecule in full.
//!
//! A window rather than an inline expand below the table: the detail view (info
//! grid + atom list + bond list) can be taller than the panel holding the
//! table, with no way to scroll to the rest of it. A window is independently
//! movable and resizable, and scrolls on its own.
//!
//! One molecule at a time, for now — the egui id below is fixed, so a second
//! selection replaces the first. #107 makes this a collection.

use crate::molecule_view::{show_atom_list, show_bond_list, show_molecule_info};
use crate::state::AppState;
use crate::structure_view::structure_panel_with_options;
use chemcore::layout::ensure_coords;
use chemcore::molecule::Molecule;

#[derive(Default)]
pub struct DetailView {
    /// The laid-out copy of the selected molecule, and the row it came from.
    ///
    /// The window rebuilds its contents every frame and a molecule from SMILES
    /// has to be laid out before it can be drawn, so the laid-out copy is
    /// cached against its row rather than recomputed each frame.
    molecule: Option<(usize, Molecule)>,
    dataset_epoch: u64,
}

impl DetailView {
    pub fn show(&mut self, ctx: &egui::Context, state: &mut AppState) {
        self.sync(state);

        let Some(i) = state.selected_row else {
            return;
        };
        let active_dataset = state.loaded_files.active_dataset();
        let Some(source) = active_dataset.molecules.get(i) else {
            return;
        };
        let name = active_dataset.names[i].clone();
        let smiles = active_dataset.smiles[i].clone();

        // Molecules parsed from SMILES carry no coordinates, so one is
        // generated here; SDF-sourced molecules keep the layout their file
        // supplied. Done once per selected row rather than every frame, since
        // laying out a large molecule isn't free.
        let needs_layout = !matches!(&self.molecule, Some((cached, _)) if *cached == i);
        if needs_layout {
            let mut prepared = source.clone();
            ensure_coords(&mut prepared);
            self.molecule = Some((i, prepared));
        }
        let Some((_, mol)) = &self.molecule else {
            return;
        };
        let mol = mol.clone();
        let options = state.display.structure;

        // Caps the window itself to the available viewport height (minus some
        // margin) so it can never grow past the screen on a short
        // window/monitor — content beyond that scrolls inside instead of the
        // window just growing off-screen with no way to reach it.
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
                    structure_panel_with_options(ui, &mol, 220.0, options);
                    show_molecule_info(ui, &mol, &smiles, &name);
                    show_atom_list(ui, &mol);
                    show_bond_list(ui, &mol);
                });
            });

        if !open {
            state.selected_row = None;
            self.molecule = None;
        }
    }

    /// Drops a layout cached against a row of a dataset that is no longer
    /// active. `AppState` clears the selection itself; this is the copy of the
    /// molecule that went with it.
    fn sync(&mut self, state: &AppState) {
        if self.dataset_epoch != state.dataset_epoch() {
            self.dataset_epoch = state.dataset_epoch();
            self.molecule = None;
        }
    }
}
