use bitvec::prelude::BitVec;
use chemfp::tanimoto::tanimoto_similarity;
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

/// Cell edge for grids compared against each other inside a result row, where
/// two have to fit one above the other.
const CELL_SIZE_COMPACT: f32 = 4.0;

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
                    Color32::from_rgb(70, 130, 180)
                } else {
                    Color32::from_rgb(240, 240, 240)
                };

                painter.rect_filled(cell_rect, 0.0, color);
            }
        }

        response
    }
}

/// Two fingerprints, one above the other, and the arithmetic behind their score.
///
/// Stacked rather than side by side: two grids are 64 cells wide, and a pair of
/// them will not fit a window's width, so one was being cut off. Stacking also
/// puts bit *n* at the same x in both, which is what makes a difference readable
/// by running an eye down a column.
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

        ui.horizontal_wrapped(|ui| {
            ui.label(egui::RichText::new("Shared bits").strong());
            ui.label(format!("{intersection}"));
            ui.separator();
            ui.label(egui::RichText::new("Either").strong());
            ui.label(format!("{union}"));
            ui.separator();
            ui.label(
                egui::RichText::new(format!("{intersection}/{union} = {tanimoto:.3}")).strong(),
            );
        });

        ui.add_space(2.0);
        ui.label(egui::RichText::new("This molecule").small());
        ui.add(FingerprintView::new(molecule).compact());
        ui.label(egui::RichText::new("Query").small());
        ui.add(FingerprintView::new(query).compact());
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
    fn test_a_grid_is_the_same_shape_whatever_its_cell_size() {
        // The two views differed in grid width, so the same bits were a
        // landscape grid in one place and a portrait one in another, with no bit
        // in the same position in both. Only the cell size varies now.
        let fp = fingerprint(&[1; 128]);
        let full = FingerprintView::new(&fp);
        let compact = FingerprintView::new(&fp).compact();

        assert_eq!(full.stats().total_bits, compact.stats().total_bits);
        assert_eq!(128_usize.div_ceil(GRID_WIDTH), 2);
    }
}
