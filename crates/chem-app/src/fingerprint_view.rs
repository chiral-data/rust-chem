use bitvec::prelude::BitVec;
use chem_fp::tanimoto::tanimoto_similarity;
use egui::{Color32, Pos2, Rect, Response, Sense, Ui, Vec2, Widget};

/// Bits per row, everywhere a fingerprint is drawn.
///
/// One number rather than one per call site: two grids of the same bits are only
/// comparable if a bit lands in the same place in both, and they were 64 wide in
/// one view and 32 in another — the same 2048 bits as a landscape grid and a
/// portrait one, with nothing in common to read across.
const GRID_WIDTH: usize = 64;

/// Cell edge, in points, for a grid shown on its own.
const CELL_SIZE: f32 = 8.0;

/// Cell edge for the comparison grid inside a result row.
const CELL_SIZE_COMPACT: f32 = 5.0;

/// A bit set in both fingerprints — the numerator of the score.
const BOTH: Color32 = Color32::from_rgb(70, 130, 180);
/// Set in the molecule but not the query.
const MOLECULE_ONLY: Color32 = Color32::from_rgb(214, 148, 48);
/// Set in the query but not the molecule.
const QUERY_ONLY: Color32 = Color32::from_rgb(142, 110, 190);

/// Background for a bit set in neither.
///
/// Taken from the theme rather than fixed: this was a hardcoded near-white,
/// which was invisible against a light background and glaring against a dark
/// one. It only became wrong when #121 made the theme a choice.
fn unset_color(ui: &Ui) -> Color32 {
    ui.visuals().faint_bg_color
}

pub struct FingerprintView<'a> {
    fingerprint: &'a BitVec,
    cell_size: f32,
}

impl<'a> FingerprintView<'a> {
    pub fn new(fingerprint: &'a BitVec) -> Self {
        Self {
            fingerprint,
            cell_size: CELL_SIZE,
        }
    }

    pub fn compact(mut self) -> Self {
        self.cell_size = CELL_SIZE_COMPACT;
        self
    }

    pub fn stats(&self) -> FingerprintStats {
        let total_bits = self.fingerprint.len();
        let set_bits = self.fingerprint.count_ones();
        let density = if total_bits > 0 {
            (set_bits as f64 / total_bits as f64) * 100.0
        } else {
            0.0
        };

        FingerprintStats {
            total_bits,
            set_bits,
            density,
        }
    }
}

pub struct FingerprintStats {
    pub total_bits: usize,
    pub set_bits: usize,
    pub density: f64,
}

impl<'a> Widget for FingerprintView<'a> {
    fn ui(self, ui: &mut Ui) -> Response {
        let total_bits = self.fingerprint.len();
        let grid_height = total_bits.div_ceil(GRID_WIDTH);

        let cell_size = self.cell_size;
        let spacing = (cell_size / 8.0).max(0.5);
        let total_width = (cell_size + spacing) * GRID_WIDTH as f32;
        let total_height = (cell_size + spacing) * grid_height as f32;

        let (rect, response) =
            ui.allocate_exact_size(Vec2::new(total_width, total_height), Sense::hover());

        let unset = unset_color(ui);
        if ui.is_rect_visible(rect) {
            let painter = ui.painter();

            for i in 0..total_bits {
                let row = i / GRID_WIDTH;
                let col = i % GRID_WIDTH;

                let x = rect.min.x + col as f32 * (cell_size + spacing);
                let y = rect.min.y + row as f32 * (cell_size + spacing);

                let cell_rect =
                    Rect::from_min_size(Pos2::new(x, y), Vec2::new(cell_size, cell_size));

                let color = if *self.fingerprint.get(i).unwrap() {
                    BOTH
                } else {
                    unset
                };

                painter.rect_filled(cell_rect, 0.0, color);
            }
        }

        response
    }
}

/// Both fingerprints in one grid, each bit coloured by which of them has it set.
///
/// This *is* the Tanimoto calculation drawn. Every bit is in one of four states,
/// and three of them are the sum: set in both is the numerator, set in either
/// alone joins it to make the denominator, set in neither is not counted. Two
/// stacked grids made you do that intersection by eye across two rows; one grid
/// has already done it.
struct FingerprintDiffView<'a> {
    molecule: &'a BitVec,
    query: &'a BitVec,
}

impl Widget for FingerprintDiffView<'_> {
    fn ui(self, ui: &mut Ui) -> Response {
        let total_bits = self.molecule.len().min(self.query.len());
        let grid_height = total_bits.div_ceil(GRID_WIDTH);

        let cell_size = CELL_SIZE_COMPACT;
        let spacing = 1.0;
        let (rect, response) = ui.allocate_exact_size(
            Vec2::new(
                (cell_size + spacing) * GRID_WIDTH as f32,
                (cell_size + spacing) * grid_height as f32,
            ),
            Sense::hover(),
        );

        let unset = unset_color(ui);
        if ui.is_rect_visible(rect) {
            let painter = ui.painter();
            for i in 0..total_bits {
                let in_molecule = self.molecule[i];
                let in_query = self.query[i];
                let color = match (in_molecule, in_query) {
                    (true, true) => BOTH,
                    (true, false) => MOLECULE_ONLY,
                    (false, true) => QUERY_ONLY,
                    (false, false) => unset,
                };

                let x = rect.min.x + (i % GRID_WIDTH) as f32 * (cell_size + spacing);
                let y = rect.min.y + (i / GRID_WIDTH) as f32 * (cell_size + spacing);
                painter.rect_filled(
                    Rect::from_min_size(Pos2::new(x, y), Vec2::splat(cell_size)),
                    0.0,
                    color,
                );
            }
        }

        response
    }
}

/// A colour swatch and its count, for the grid's legend.
fn legend_entry(ui: &mut Ui, color: Color32, label: &str, count: usize) {
    ui.label(egui::RichText::new("\u{25a0}").color(color));
    ui.label(egui::RichText::new(format!("{label} {count}")).small());
}

/// Why a result scored what it did: the arithmetic, and both fingerprints in one
/// colour-coded grid.
///
/// One grid rather than two. Side by side they did not fit a window's width, so
/// one was cut off; stacked they fit but left the reader intersecting two rows by
/// eye. Colouring a single grid by which fingerprint holds each bit does that
/// intersection for them, and is the same three numbers the score is made of.
pub fn fingerprint_comparison(ui: &mut Ui, molecule: &BitVec, query: Option<&BitVec>) {
    ui.group(|ui| {
        let Some(query) = query else {
            ui.label(egui::RichText::new("Molecule fingerprint").strong());
            ui.add(FingerprintView::new(molecule).compact());
            ui.label(
                egui::RichText::new("No query fingerprint to compare against.")
                    .small()
                    .weak(),
            );
            return;
        };

        if molecule.len() != query.len() {
            // Changing the fingerprint size after computing the dataset leaves
            // the two incomparable. Saying so beats drawing two grids of
            // different heights and letting them look comparable.
            ui.label(
                egui::RichText::new(format!(
                    "Not comparable: this molecule has {} bits, the query {}. \
                     Recompute fingerprints.",
                    molecule.len(),
                    query.len()
                ))
                .small()
                .color(Color32::from_rgb(220, 120, 50)),
            );
            return;
        }

        // Tanimoto is the intersection over the union, so these numbers *are*
        // the score, which is what the expansion is being asked for.
        //
        // The ratio comes from the same function that ranked the results rather
        // than being recomputed here. Reimplementing it would let the
        // explanation drift from the score it is explaining — silently, and in
        // the one place a user goes to be told why.
        let intersection = (molecule.clone() & query.clone()).count_ones();
        let union = (molecule.clone() | query.clone()).count_ones();
        let tanimoto = tanimoto_similarity(molecule, query).unwrap_or(0.0);

        let molecule_only = (molecule.clone() & !query.clone()).count_ones();
        let query_only = (query.clone() & !molecule.clone()).count_ones();

        ui.horizontal_wrapped(|ui| {
            legend_entry(ui, BOTH, "both", intersection);
            ui.separator();
            legend_entry(ui, MOLECULE_ONLY, "this only", molecule_only);
            ui.separator();
            legend_entry(ui, QUERY_ONLY, "query only", query_only);
        });
        ui.label(
            egui::RichText::new(format!(
                "{intersection} shared / {union} in either = {tanimoto:.3}"
            ))
            .strong(),
        );

        ui.add_space(2.0);
        ui.add(FingerprintDiffView { molecule, query });
    });
}

pub fn fingerprint_full(ui: &mut Ui, fingerprint: &BitVec) {
    ui.group(|ui| {
        ui.label(
            egui::RichText::new("Fingerprint Visualization")
                .strong()
                .size(16.0),
        );

        let stats = FingerprintView::new(fingerprint).stats();

        ui.horizontal(|ui| {
            ui.label(format!("Total bits: {}", stats.total_bits));
            ui.separator();
            ui.label(format!("Bits set: {}", stats.set_bits));
            ui.separator();
            ui.label(format!("Density: {:.2}%", stats.density));
        });

        ui.separator();

        egui::ScrollArea::vertical()
            .max_height(300.0)
            .show(ui, |ui| {
                ui.add(FingerprintView::new(fingerprint));
            });
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fingerprint(bits: &[u8]) -> BitVec {
        bits.iter().map(|&b| b != 0).collect()
    }

    #[test]
    fn test_the_counts_shown_reproduce_the_score() {
        // The panel shows "intersection/union = score". If the similarity
        // measure ever changed — to Dice, say — the counts would still be
        // right and the caption would quietly become a lie, so this pins the
        // relationship the caption claims.
        let a = fingerprint(&[1, 1, 0, 1, 0, 0, 1, 0]);
        let b = fingerprint(&[1, 0, 0, 1, 1, 0, 1, 1]);

        let intersection = (a.clone() & b.clone()).count_ones();
        let union = (a.clone() | b.clone()).count_ones();
        let score = tanimoto_similarity(&a, &b).unwrap();

        assert_eq!(intersection, 3);
        assert_eq!(union, 6);
        assert!((score - intersection as f64 / union as f64).abs() < 1e-12);
    }

    #[test]
    fn test_identical_fingerprints_score_one() {
        // Two spellings of the same molecule must land here. They once scored
        // 0.27 (#95), and this panel is where that would now be visible.
        let a = fingerprint(&[1, 0, 1, 1, 0, 0, 1, 0]);
        assert!((tanimoto_similarity(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_every_grid_shares_one_width() {
        // The views used to take grid width as an argument, so the same 2048
        // bits were a landscape grid in one place and a portrait one in another,
        // with no bit in the same position in both. There is one width now.
        assert_eq!(128_usize.div_ceil(GRID_WIDTH), 2);
        assert_eq!(2048_usize.div_ceil(GRID_WIDTH), 32);
    }

    #[test]
    fn test_the_four_bit_states_partition_the_grid() {
        // The three coloured states are exactly the union, and the fourth is
        // everything the score ignores. If that stopped holding, the legend
        // would be counting something other than what the grid draws.
        let a = fingerprint(&[1, 1, 0, 1, 0, 0, 1, 0]);
        let b = fingerprint(&[1, 0, 0, 1, 1, 0, 1, 1]);

        let both = (a.clone() & b.clone()).count_ones();
        let a_only = (a.clone() & !b.clone()).count_ones();
        let b_only = (b.clone() & !a.clone()).count_ones();
        let neither = a.len() - (both + a_only + b_only);

        assert_eq!(both, 3);
        assert_eq!(a_only, 1);
        assert_eq!(b_only, 2);
        assert_eq!(both + a_only + b_only, (a.clone() | b.clone()).count_ones());
        assert_eq!(both + a_only + b_only + neither, a.len());
    }
}
