//! Preferences: the things that govern the whole app rather than one window.
//!
//! The structure options used to live in the Inspector, on the principle that
//! options belong with what they govern (#106). They govern *three* surfaces
//! though — the dataset table's thumbnails, the query structure, and every
//! detail window — and only one of those is the Inspector, so hosting them
//! there was a milder version of the problem #106 fixed. Being closed by
//! default in a window whose job is showing data was the hint.
//!
//! A window rather than a menu or a modal, because the decisive constraint is
//! live preview: you turn hydrogens on to *see* what it does, so this has to be
//! visible at the same time as a structure. A modal blocks interaction and a
//! menu closes on it.
//!
//! Closed by default, so the default layout stays the three content windows.

use crate::state::AppState;
use crate::structure_view::structure_option_controls;
use egui::RichText;

#[derive(Default)]
pub struct SettingsView;

impl SettingsView {
    pub fn ui(&mut self, ui: &mut egui::Ui, state: &mut AppState) {
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                ui.label(RichText::new("Theme").strong());
                // egui's own picker: System, Dark, Light, with the detected
                // system theme on hover.
                state.display.theme.radio_buttons(ui);

                ui.add_space(8.0);
                ui.separator();

                ui.label(RichText::new("Structures").strong());
                ui.label(
                    RichText::new("Applies everywhere a structure is drawn.")
                        .small()
                        .weak(),
                );
                structure_option_controls(ui, &mut state.display.structure);

                ui.add_space(4.0);
                // Back alongside the structure options. #106 split this off to
                // sit beside the table it governs, because a toggle in the
                // Inspector affecting the Datasets table was worse. In a
                // preferences window neither is "somewhere else".
                ui.checkbox(
                    &mut state.display.show_thumbnails,
                    "Show structures in the dataset table",
                );
            });
    }
}
