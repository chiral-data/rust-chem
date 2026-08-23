//! The egui side of 2D structure depiction.
//!
//! Describing a structure is [`chemdraw`]; this is the backend that paints the
//! description, plus the pieces that need egui itself — measuring labels with
//! the app's own font, following the app's light/dark setting, and the option
//! controls.

use chemcore::molecule::Molecule;
use chemdraw::structure::{
    AtomVisualization, ShowCarbons, StructureOptions, StructureShape, StructureTheme,
    describe_structure,
};
use egui::{Color32, FontId, Response, Sense, Shape, Stroke, Ui, Vec2, Widget};

/// A widget drawing a molecule's 2D structure.
///
/// Requires the molecule to carry coordinates — SDF supplies them, SMILES
/// doesn't. Molecules without a layout render a placeholder rather than
/// nothing, so the reason for an empty panel is visible.
pub struct StructureView<'a> {
    molecule: &'a Molecule,
    options: StructureOptions,
    desired_size: Vec2,
}

impl<'a> StructureView<'a> {
    pub fn new(molecule: &'a Molecule, desired_size: Vec2) -> Self {
        Self {
            molecule,
            options: StructureOptions::default(),
            desired_size,
        }
    }

    pub fn with_options(mut self, options: StructureOptions) -> Self {
        self.options = options;
        self
    }

    // A `with_theme` builder for pinning a fixed palette (rather than
    // following the app's light/dark setting) belongs here, but nothing needs
    // it until there's an export path; adding it now would just be dead code.
}

/// A structure palette following egui's current visuals, so structures track
/// the app's theme rather than needing a setting of their own.
///
/// A free function rather than `StructureTheme::from_visuals`: the type is
/// `chemdraw`'s now, and an inherent method would have to live with it — which
/// would put `egui` back in the crate the extraction took it out of.
pub fn theme_from_visuals(visuals: &egui::Visuals) -> StructureTheme {
    if visuals.dark_mode {
        StructureTheme::dark()
    } else {
        StructureTheme::light()
    }
}

/// Draws a description with egui.
fn paint_structure(painter: &egui::Painter, shapes: &[StructureShape]) {
    for shape in shapes {
        match shape {
            StructureShape::Line {
                from,
                to,
                width,
                color,
            } => {
                painter.line_segment([*from, *to], Stroke::new(*width, *color));
            }
            StructureShape::DashedLine {
                from,
                to,
                width,
                color,
                dash,
            } => {
                painter.extend(Shape::dashed_line(
                    &[*from, *to],
                    Stroke::new(*width, *color),
                    *dash,
                    *dash,
                ));
            }
            StructureShape::Disc {
                center,
                radius,
                color,
            } => {
                painter.circle_filled(*center, *radius, *color);
            }
            StructureShape::Text {
                pos,
                align,
                text,
                size,
                color,
            } => {
                painter.text(*pos, *align, text, FontId::proportional(*size), *color);
            }
        }
    }
}

impl Widget for StructureView<'_> {
    fn ui(self, ui: &mut Ui) -> Response {
        let (rect, response) = ui.allocate_exact_size(self.desired_size, Sense::hover());

        if !ui.is_rect_visible(rect) {
            return response;
        }

        let theme = theme_from_visuals(ui.visuals());
        let weak_color = ui.visuals().weak_text_color();

        // Scoped so the measuring borrow ends before the painting one begins.
        let shapes = {
            let measure = |text: &str, size: f32| {
                ui.painter()
                    .layout_no_wrap(
                        text.to_owned(),
                        FontId::proportional(size),
                        Color32::PLACEHOLDER,
                    )
                    .size()
            };
            describe_structure(
                self.molecule,
                rect,
                &self.options,
                &theme,
                weak_color,
                &measure,
            )
        };

        paint_structure(ui.painter(), &shapes);
        response
    }
}

/// Controls for the atom-display options.
///
/// These apply to every structure the app draws — table thumbnails, the query,
/// detail windows — so they belong in one findable place rather than beside any
/// one of them. #106 puts that place in the Inspector window.
///
/// `show_thumbnails` used to ride along here. It doesn't any more: it governs
/// one table, in another window, and a toggle whose effect appears somewhere
/// else is worse than one extra checkbox.
pub fn structure_option_controls(ui: &mut Ui, options: &mut StructureOptions) {
    ui.horizontal_wrapped(|ui| {
        ui.label("Carbons:");
        egui::ComboBox::from_id_salt("show_carbons")
            .selected_text(options.show_carbons.label())
            .show_ui(ui, |ui| {
                for mode in ShowCarbons::ALL {
                    ui.selectable_value(&mut options.show_carbons, mode, mode.label());
                }
            });

        ui.separator();
        ui.label("Atoms:");
        egui::ComboBox::from_id_salt("atom_visualization")
            .selected_text(options.atom_visualization.label())
            .show_ui(ui, |ui| {
                for mode in AtomVisualization::ALL {
                    ui.selectable_value(&mut options.atom_visualization, mode, mode.label());
                }
            });

        ui.separator();
        ui.checkbox(&mut options.explicit_hydrogens, "Hydrogens");
    });
}

/// Draws a molecule's structure in a bordered group with a heading, fitted to
/// the available width.
pub fn structure_panel_with_options(
    ui: &mut Ui,
    molecule: &Molecule,
    height: f32,
    options: StructureOptions,
) {
    ui.group(|ui| {
        ui.label(egui::RichText::new("Structure").strong());
        let width = ui.available_width();
        ui.add(StructureView::new(molecule, Vec2::new(width, height)).with_options(options));
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_theme_follows_visuals() {
        assert_eq!(
            theme_from_visuals(&egui::Visuals::dark()).foreground,
            StructureTheme::dark().foreground
        );
        assert_eq!(
            theme_from_visuals(&egui::Visuals::light()).foreground,
            StructureTheme::light().foreground
        );
    }
}
