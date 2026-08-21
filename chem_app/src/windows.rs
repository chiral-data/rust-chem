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

use egui::{Rect, pos2, vec2};

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
                default_rect: Rect::from_min_size(pos2(24.0, 48.0), vec2(400.0, 300.0)),
            },
            operations: WindowEntry {
                id: "window_operations",
                title: "Operations",
                open: true,
                default_rect: Rect::from_min_size(pos2(448.0, 48.0), vec2(340.0, 300.0)),
            },
            visualization: WindowEntry {
                id: "window_visualization",
                title: "Data Visualization",
                open: true,
                default_rect: Rect::from_min_size(pos2(120.0, 380.0), vec2(560.0, 260.0)),
            },
            laid_out: false,
        }
    }
}

/// A workspace smaller than this in either direction isn't worth laying out; the
/// fallback rects and `constrain` handle it instead.
const MIN_TILEABLE: f32 = 240.0;

// The default layout, as fractions of the workspace. The windows deliberately
// don't fill it: free canvas between and around them is what makes them read as
// floating windows rather than panels that happen to have title bars, and it
// leaves somewhere to drop a molecule detail window (#107) without it landing
// on top of something.
//
// Whatever these leave unused becomes margin, split evenly, so the cluster sits
// in the middle of any canvas rather than hugging a corner.
const DATA_SOURCES_W: f32 = 0.34;
const OPERATIONS_W: f32 = 0.28;
const VISUALIZATION_W: f32 = 0.46;
const COL_GAP: f32 = 0.06;
const TOP_ROW_H: f32 = 0.40;
const BOTTOM_ROW_H: f32 = 0.30;
const ROW_GAP: f32 = 0.08;

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

        let (w, h) = (workspace.width(), workspace.height());
        // Fractions in, absolute positions out.
        let at = |fx: f32, fy: f32| workspace.min + vec2(w * fx, h * fy);
        let size = |fw: f32, fh: f32| vec2(w * fw, h * fh);

        // Data Sources and Operations side by side, Data Visualization centred
        // beneath the pair. Unused fractions become even margins.
        let top_row_w = DATA_SOURCES_W + COL_GAP + OPERATIONS_W;
        let x0 = (1.0 - top_row_w) / 2.0;
        let y0 = (1.0 - (TOP_ROW_H + ROW_GAP + BOTTOM_ROW_H)) / 2.0;

        self.data_sources.default_rect =
            Rect::from_min_size(at(x0, y0), size(DATA_SOURCES_W, TOP_ROW_H));
        self.operations.default_rect = Rect::from_min_size(
            at(x0 + DATA_SOURCES_W + COL_GAP, y0),
            size(OPERATIONS_W, TOP_ROW_H),
        );
        self.visualization.default_rect = Rect::from_min_size(
            at(
                x0 + (top_row_w - VISUALIZATION_W) / 2.0,
                y0 + TOP_ROW_H + ROW_GAP,
            ),
            size(VISUALIZATION_W, BOTTOM_ROW_H),
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

    /// Clearance between two disjoint rects, along whichever axis separates
    /// them. Negative when they overlap.
    fn separation(a: Rect, b: Rect) -> f32 {
        let dx = (b.min.x - a.max.x).max(a.min.x - b.max.x);
        let dy = (b.min.y - a.max.y).max(a.min.y - b.max.y);
        dx.max(dy)
    }

    #[test]
    fn test_windows_leave_visible_space_between_them() {
        for (w, h) in VIEWPORTS {
            let (_, mut registry) = tiled(w, h);
            let rects: Vec<Rect> = registry
                .entries_mut()
                .iter()
                .map(|e| e.default_rect)
                .collect();
            // Free canvas is the point of the arrangement: windows that merely
            // fail to overlap still read as panels butted together.
            let min_gap = w.min(h) * 0.04;

            for i in 0..rects.len() {
                for j in (i + 1)..rects.len() {
                    let gap = separation(rects[i], rects[j]);
                    assert!(
                        gap >= min_gap,
                        "windows {} and {} sit {:.0}px apart at {}x{}, want >= {:.0}px",
                        i,
                        j,
                        gap,
                        w,
                        h,
                        min_gap
                    );
                }
            }
        }
    }

    #[test]
    fn test_windows_leave_a_margin_around_the_workspace() {
        for (w, h) in VIEWPORTS {
            let (workspace, mut registry) = tiled(w, h);
            let inset = workspace.shrink(w.min(h) * 0.04);
            for entry in registry.entries_mut() {
                assert!(
                    inset.contains_rect(entry.default_rect),
                    "{} reaches the workspace edge at {}x{}: {:?}",
                    entry.title,
                    w,
                    h,
                    entry.default_rect
                );
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
