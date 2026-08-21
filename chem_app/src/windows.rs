//! Which windows exist, whether they are open, and where they start.
//!
//! The shell used to draw four regions unconditionally, so "which windows
//! exist" was implicit in the order of five calls at the end of `update()`.
//! Holding it as data instead is what lets the View menu be generated rather
//! than hand-written, and gives #108 one place to save and restore.
//!
//! Deliberately *not* a `dyn Window` trait. The views' `ui` signatures aren't
//! uniform yet — Data Sources still draws the fingerprint controls whose
//! parameters Operations owns, until #105 moves those widgets — so a trait
//! would have to carry that argument to every view to accommodate one caller.
//! #99 turned down an abstraction for the same reason. What every window really
//! does share is its *shell*, so that is what this holds; the calls that draw
//! their contents stay concrete and typed in `app.rs`.

use egui::{Pos2, Vec2, pos2, vec2};

/// One window's shell: identity, title, whether it is open, and where it opens
/// the first time.
pub struct WindowEntry {
    /// Stable across renames and reorderings. It keys egui's own memory for the
    /// window's position and size, so changing it moves a user's window; #108
    /// will persist layout against it too.
    id: &'static str,
    title: &'static str,
    pub open: bool,
    default_pos: Pos2,
    default_size: Vec2,
}

impl WindowEntry {
    pub fn title(&self) -> &'static str {
        self.title
    }

    /// Draws the window if it is open, and does nothing if it is not.
    pub fn show(&mut self, ctx: &egui::Context, add_contents: impl FnOnce(&mut egui::Ui)) {
        if !self.open {
            return;
        }
        // Destructured because `.open()` needs `&mut self.open` while the title
        // and id are read from the same struct.
        let Self {
            id,
            title,
            open,
            default_pos,
            default_size,
        } = self;

        egui::Window::new(*title)
            .id(egui::Id::new(*id))
            .open(open)
            .default_pos(*default_pos)
            .default_size(*default_size)
            .resizable(true)
            .collapsible(true)
            // Keeps a window reachable when the viewport is smaller than the
            // default layout assumes — a narrow browser canvas, or a laptop
            // screen rather than the 1400x900 the native build asks for.
            .constrain(true)
            .show(ctx, add_contents);
    }
}

pub struct WindowRegistry {
    pub data_sources: WindowEntry,
    pub operations: WindowEntry,
    pub visualization: WindowEntry,
}

impl Default for WindowRegistry {
    fn default() -> Self {
        // Laid out as the milestone's sketch has it: sources and operations
        // side by side along the top, visualization across the bottom. All
        // three start open — a first frame showing an empty canvas would be a
        // poor introduction to an app that has just gained three windows.
        //
        // Positions assume the native build's 1400x900 viewport and are only
        // *defaults*: egui remembers where a window was dragged, and
        // `constrain` keeps them on-screen on a smaller canvas.
        Self {
            data_sources: WindowEntry {
                id: "window_data_sources",
                title: "Data Sources",
                open: true,
                default_pos: pos2(16.0, 44.0),
                default_size: vec2(540.0, 400.0),
            },
            operations: WindowEntry {
                id: "window_operations",
                title: "Operations",
                open: true,
                default_pos: pos2(572.0, 44.0),
                default_size: vec2(420.0, 400.0),
            },
            visualization: WindowEntry {
                id: "window_visualization",
                title: "Data Visualization",
                open: true,
                default_pos: pos2(16.0, 464.0),
                default_size: vec2(976.0, 400.0),
            },
        }
    }
}

impl WindowRegistry {
    /// Every window, in menu order. What the View menu is built from, and what
    /// #108 will iterate to save open state and geometry.
    pub fn entries_mut(&mut self) -> [&mut WindowEntry; 3] {
        [
            &mut self.data_sources,
            &mut self.operations,
            &mut self.visualization,
        ]
    }

    /// True when the workspace has nothing on it. The canvas draws a way back
    /// in that case, since a bare canvas offers no clue that the View menu is
    /// where windows come from.
    pub fn all_closed(&self) -> bool {
        !self.data_sources.open && !self.operations.open && !self.visualization.open
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_every_window_starts_open() {
        let mut registry = WindowRegistry::default();
        assert!(!registry.all_closed());
        assert!(registry.entries_mut().iter().all(|e| e.open));
    }

    #[test]
    fn test_closing_every_window_is_noticed() {
        let mut registry = WindowRegistry::default();
        for entry in registry.entries_mut() {
            entry.open = false;
        }
        // Otherwise the canvas would sit blank with nothing pointing at the
        // View menu.
        assert!(registry.all_closed());
    }

    #[test]
    fn test_window_ids_are_unique() {
        let mut registry = WindowRegistry::default();
        let ids: Vec<&str> = registry.entries_mut().iter().map(|e| e.id).collect();
        let mut deduped = ids.clone();
        deduped.sort_unstable();
        deduped.dedup();
        // A shared id would make two windows share one position in egui's
        // memory, and would collide again in #108's saved layout.
        assert_eq!(ids.len(), deduped.len(), "duplicate window id: {:?}", ids);
    }
}
