use bitvec::prelude::BitVec;
use egui::{Color32, Pos2, Rect, Response, Sense, Ui, Vec2, Widget};

pub struct FingerprintView<'a> {
    fingerprint: &'a BitVec,
    grid_width: usize,
}

impl<'a> FingerprintView<'a> {
    pub fn new(fingerprint: &'a BitVec, grid_width: usize) -> Self {
        Self {
            fingerprint,
            grid_width,
        }
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
        let grid_height = total_bits.div_ceil(self.grid_width);

        let cell_size = 8.0;
        let spacing = 1.0;
        let total_width = (cell_size + spacing) * self.grid_width as f32;
        let total_height = (cell_size + spacing) * grid_height as f32;

        let (rect, response) =
            ui.allocate_exact_size(Vec2::new(total_width, total_height), Sense::hover());

        if ui.is_rect_visible(rect) {
            let painter = ui.painter();

            for i in 0..total_bits {
                let row = i / self.grid_width;
                let col = i % self.grid_width;

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

pub fn fingerprint_compact(ui: &mut Ui, fingerprint: &BitVec, label: &str) {
    ui.group(|ui| {
        ui.label(egui::RichText::new(label).strong());

        let stats = FingerprintView::new(fingerprint, 32).stats();

        ui.horizontal(|ui| {
            ui.label(format!("Size: {} bits", stats.total_bits));
            ui.separator();
            ui.label(format!("Set: {}", stats.set_bits));
            ui.separator();
            ui.label(format!("Density: {:.2}%", stats.density));
        });

        ui.add(FingerprintView::new(fingerprint, 32));
    });
}

pub fn fingerprint_full(ui: &mut Ui, fingerprint: &BitVec) {
    ui.group(|ui| {
        ui.label(
            egui::RichText::new("Fingerprint Visualization")
                .strong()
                .size(16.0),
        );

        let stats = FingerprintView::new(fingerprint, 64).stats();

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
                ui.add(FingerprintView::new(fingerprint, 64));
            });
    });
}
