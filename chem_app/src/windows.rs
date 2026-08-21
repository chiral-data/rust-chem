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

use egui::{Rect, Vec2, pos2, vec2};

/// One window's shell: identity, title, whether it is open, and where it opens
/// the first time.
pub struct WindowEntry {
    /// Stable across renames and reorderings. It keys egui's own memory for the
    /// window's position and size, so changing it moves a user's window; #108
    /// will persist layout against it too.
    id: &'static str,
    title: &'static str,
    pub open: bool,
    /// Where the window opens the first time. Assigned from the real viewport
    /// by [`WindowRegistry::ensure_layout`] rather than hard-coded, so the
    /// default layout tiles whatever canvas it is given; the value here is only
    /// the fallback for a viewport too small to tile.
    default_rect: Rect,
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
            default_rect,
        } = self;

        egui::Window::new(*title)
            .id(egui::Id::new(*id))
            .open(open)
            .default_rect(*default_rect)
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
    /// Whether [`WindowRegistry::ensure_layout`] has run. The default rects
    /// depend on the viewport, which egui doesn't know on the very first frame.
    laid_out: bool,
}

impl Default for WindowRegistry {
    fn default() -> Self {
        // All three start open — a first frame showing an empty canvas would be
        // a poor introduction to an app that has just gained three windows.
        //
        // These rects are fallbacks only, for a viewport too small to tile.
        // `ensure_layout` replaces them with a tiling of the real workspace
        // before any window is shown for the first time.
        Self {
            data_sources: WindowEntry {
                id: "window_data_sources",
                title: "Data Sources",
                open: true,
                default_rect: Rect::from_min_size(pos2(16.0, 44.0), vec2(540.0, 400.0)),
            },
            operations: WindowEntry {
                id: "window_operations",
                title: "Operations",
                open: true,
                default_rect: Rect::from_min_size(pos2(572.0, 44.0), vec2(420.0, 400.0)),
            },
            visualization: WindowEntry {
                id: "window_visualization",
                title: "Data Visualization",
                open: true,
                default_rect: Rect::from_min_size(pos2(16.0, 464.0), vec2(976.0, 400.0)),
            },
            laid_out: false,
        }
    }
}

/// Gap between tiled windows, and between a window and the workspace edge.
const PAD: f32 = 8.0;

/// A workspace smaller than this in either direction isn't worth tiling; the
/// fallback rects and `constrain` handle it instead.
const MIN_TILEABLE: f32 = 240.0;

impl WindowRegistry {
    /// Tiles the three windows across the workspace, once, on the first frame
    /// that reports a usable size.
    ///
    /// The defaults used to be hard-coded against the native build's 1400x900
    /// viewport, which meant they were wrong everywhere else — a 1366x768 laptop
    /// or a browser canvas put Data Visualization partly off the bottom, and
    /// `constrain` then pulled it back on top of its neighbours. Deriving the
    /// layout from the real workspace tiles any size without overlap.
    pub fn ensure_layout(&mut self, workspace: Rect) {
        if self.laid_out || workspace.width() < MIN_TILEABLE || workspace.height() < MIN_TILEABLE {
            return;
        }
        self.laid_out = true;

        // Two windows side by side along the top, one across the bottom: three
        // gaps horizontally (edge, middle, edge) and three vertically.
        let inner = workspace.size() - Vec2::splat(PAD * 3.0);
        let top_h = inner.y * 0.55;
        let left_w = inner.x * 0.56;
        let origin = workspace.min + Vec2::splat(PAD);

        self.data_sources.default_rect = Rect::from_min_size(origin, vec2(left_w, top_h));
        self.operations.default_rect = Rect::from_min_size(
            origin + vec2(left_w + PAD, 0.0),
            vec2(inner.x - left_w, top_h),
        );
        self.visualization.default_rect = Rect::from_min_size(
            origin + vec2(0.0, top_h + PAD),
            vec2(inner.x + PAD, inner.y - top_h),
        );
    }

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

    /// Common viewport sizes: the native build's request, a 1366x768 laptop, a
    /// modest browser canvas, and something small enough to be awkward.
    const VIEWPORTS: [(f32, f32); 4] = [
        (1400.0, 874.0),
        (1366.0, 742.0),
        (1024.0, 640.0),
        (800.0, 500.0),
    ];

    fn tiled(width: f32, height: f32) -> (Rect, WindowRegistry) {
        // y offset stands in for the menu bar above the workspace.
        let workspace = Rect::from_min_size(pos2(0.0, 26.0), vec2(width, height));
        let mut registry = WindowRegistry::default();
        registry.ensure_layout(workspace);
        (workspace, registry)
    }

    #[test]
    fn test_default_layout_never_overlaps() {
        for (w, h) in VIEWPORTS {
            let (_, mut registry) = tiled(w, h);
            let rects: Vec<Rect> = registry
                .entries_mut()
                .iter()
                .map(|e| e.default_rect)
                .collect();

            for i in 0..rects.len() {
                for j in (i + 1)..rects.len() {
                    assert!(
                        !rects[i].intersects(rects[j]),
                        "windows {} and {} overlap at {}x{}: {:?} vs {:?}",
                        i,
                        j,
                        w,
                        h,
                        rects[i],
                        rects[j]
                    );
                }
            }
        }
    }

    #[test]
    fn test_default_layout_stays_inside_the_workspace() {
        for (w, h) in VIEWPORTS {
            let (workspace, mut registry) = tiled(w, h);
            for entry in registry.entries_mut() {
                // A window placed outside gets dragged back by `constrain`,
                // which is how the hard-coded layout ended up stacking windows
                // on a viewport shorter than it assumed.
                assert!(
                    workspace.contains_rect(entry.default_rect),
                    "{} escapes a {}x{} workspace: {:?}",
                    entry.title,
                    w,
                    h,
                    entry.default_rect
                );
            }
        }
    }

    #[test]
    fn test_a_workspace_too_small_to_tile_is_left_alone() {
        let mut registry = WindowRegistry::default();
        let before = registry.data_sources.default_rect;
        // egui reports a degenerate rect before it knows the viewport; laying
        // out against that would place every window at the origin.
        registry.ensure_layout(Rect::from_min_size(pos2(0.0, 0.0), vec2(0.0, 0.0)));
        assert_eq!(registry.data_sources.default_rect, before);
        assert!(!registry.laid_out, "should retry on a later frame");
    }

    #[test]
    fn test_layout_is_computed_once() {
        let (_, mut registry) = tiled(1400.0, 874.0);
        let first = registry.data_sources.default_rect;
        // Re-tiling every frame would yank a window the user had dragged back
        // to where it started.
        registry.ensure_layout(Rect::from_min_size(pos2(0.0, 26.0), vec2(600.0, 400.0)));
        assert_eq!(registry.data_sources.default_rect, first);
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
