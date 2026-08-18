//! 2D structure depiction — draws a molecule's atoms and bonds to an egui
//! painter, given coordinates from `chemcore::Molecule`.

use chemcore::bond::BondOrder;
use chemcore::geometry::{BoundingBox, Point2};
use chemcore::molecule::Molecule;
use chemcore::rings::find_sssr;
use egui::{Align2, FontId, Pos2, Rect, Response, Sense, Shape, Stroke, Ui, Vec2, Widget};

/// Rendering options for [`StructureView`].
///
/// Proportional options are expressed relative to the on-screen bond length
/// rather than in absolute pixels. smilesDrawer can use absolute values
/// (`fontSizeLarge: 11`, `bondSpacing: 0.17 * 30`) because it also pins
/// `bondLength` to 30; here the structure is scaled to fit whatever rect it's
/// given, so the on-screen bond length varies and anything fixed in pixels
/// would look wrong at one end of the range — oversized in a table thumbnail,
/// hairline in a large detail view.
#[derive(Debug, Clone, Copy)]
pub struct StructureOptions {
    /// Pixels per coordinate unit. `0.0` means derive one instead — see
    /// [`StructureOptions::bond_length`]. smilesDrawer's `scale` defaults to
    /// `0.0` with the same meaning.
    pub scale: f32,
    /// Target on-screen bond length in pixels, used when [`Self::scale`] is
    /// `0.0` and this is greater than `0.0`. Drawing every molecule at the
    /// same bond length keeps a row of structures visually consistent, which
    /// fit-to-rect can't do — it sizes each molecule to its own box, so a
    /// two-atom molecule and a thirty-atom one come out the same width and
    /// read as though they were the same size.
    ///
    /// `0.0` falls back to fitting the structure to its rect.
    pub bond_length: f32,
    /// Margin between the structure's bounding box and the edge of its rect.
    pub padding: f32,
    /// Gap between the parallel lines of a double/triple/aromatic bond, as a
    /// fraction of bond length. smilesDrawer's `bondSpacing` is the same 0.17.
    pub bond_spacing_ratio: f32,
    /// Bond stroke width in pixels (smilesDrawer's `bondThickness`).
    pub bond_thickness: f32,
    /// Length of the secondary line of a ring bond, as a fraction of the full
    /// bond (smilesDrawer's `shortBondLength`). Trimming it stops adjacent ring
    /// bonds from running into each other at the vertices.
    ///
    /// Applies only where there's an interior to draw into: a double or
    /// aromatic bond *off* a ring puts its two lines symmetrically about the
    /// bond axis, so neither is the inner one to shorten.
    pub short_bond_length: f32,
    /// Atom label height as a fraction of bond length.
    pub font_size_ratio: f32,
    /// Clamp for the derived font size, so labels stay legible on very small
    /// renders and don't dominate very large ones.
    pub font_size_range: (f32, f32),
    /// Extra clearance between an atom label and the bonds meeting it, as a
    /// fraction of font size — keeps strokes from touching the glyphs.
    pub label_margin_ratio: f32,
}

impl Default for StructureOptions {
    fn default() -> Self {
        Self {
            scale: 0.0,
            bond_length: 0.0,
            padding: 10.0,
            bond_spacing_ratio: 0.17,
            bond_thickness: 1.0,
            short_bond_length: 0.8,
            font_size_ratio: 0.42,
            font_size_range: (7.0, 22.0),
            label_margin_ratio: 0.18,
        }
    }
}

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
}

/// Maps molecule coordinates onto screen positions.
struct Transform {
    scale: f32,
    /// Coordinate-space point that maps to the center of the target rect.
    src_center: Point2,
    dst_center: Pos2,
}

impl Transform {
    /// Uniform scale fitting `bbox` inside `rect` less `padding`, centered.
    ///
    /// Scale is uniform on both axes deliberately: scaling x and y
    /// independently would distort bond angles, which misrepresents the
    /// chemistry rather than merely looking stretched.
    fn fit(bbox: BoundingBox, rect: Rect, padding: f32) -> Self {
        let available_w = (rect.width() - 2.0 * padding).max(1.0);
        let available_h = (rect.height() - 2.0 * padding).max(1.0);

        let w = bbox.width() as f32;
        let h = bbox.height() as f32;

        // A molecule can legitimately have zero extent on an axis — ethane laid
        // out horizontally has no height, a single atom has neither — so an
        // axis with no extent places no constraint on the scale instead of
        // dividing by zero.
        let scale_x = if w > f32::EPSILON {
            available_w / w
        } else {
            f32::INFINITY
        };
        let scale_y = if h > f32::EPSILON {
            available_h / h
        } else {
            f32::INFINITY
        };

        let scale = match scale_x.min(scale_y) {
            // Both axes degenerate (a single atom): nothing to fit, so any
            // finite scale does. Pick one that leaves the glyph a sane size.
            s if s.is_infinite() => 20.0,
            s => s,
        };

        Self {
            scale,
            src_center: bbox.center(),
            dst_center: rect.center(),
        }
    }

    /// A fixed number of pixels per coordinate unit, centered on `rect`.
    fn fixed(bbox: BoundingBox, rect: Rect, scale: f32) -> Self {
        Self {
            scale,
            src_center: bbox.center(),
            dst_center: rect.center(),
        }
    }

    /// Picks a transform per the sizing options: an explicit `scale` wins, then
    /// a target `bond_length`, otherwise fit-to-rect.
    ///
    /// A target bond length can't be honoured for a molecule with no bonds to
    /// measure, or whose bonds are all zero-length; those fall back to fitting.
    fn for_options(
        molecule: &Molecule,
        bbox: BoundingBox,
        rect: Rect,
        options: &StructureOptions,
    ) -> Self {
        if options.scale > 0.0 {
            return Self::fixed(bbox, rect, options.scale);
        }

        if options.bond_length > 0.0
            && let Some(src_len) = mean_source_bond_length(molecule)
            && src_len > f64::EPSILON
        {
            return Self::fixed(bbox, rect, options.bond_length / src_len as f32);
        }

        Self::fit(bbox, rect, options.padding)
    }

    /// Molecule coordinate to screen position.
    ///
    /// The y axis is negated: chemistry coordinates put +y upward, screen
    /// coordinates put +y downward, so without this the structure renders
    /// vertically mirrored.
    fn apply(&self, p: Point2) -> Pos2 {
        let dx = (p.x - self.src_center.x) as f32 * self.scale;
        let dy = (p.y - self.src_center.y) as f32 * self.scale;
        Pos2::new(self.dst_center.x + dx, self.dst_center.y - dy)
    }
}

/// Whether this atom gets a drawn label.
///
/// Carbons are implicit vertices — a benzene ring is six lines, not six "C"
/// glyphs — so only heteroatoms are labeled. Making this configurable
/// (smilesDrawer's `showCarbons` modes, explicit hydrogens) is #80's job; this
/// is the baseline chemical convention it will override.
fn is_labeled(molecule: &Molecule, atom_idx: usize) -> bool {
    molecule.atom(atom_idx).atomic_number() != 6
}

fn atom_label(molecule: &Molecule, atom_idx: usize) -> String {
    molecule.atom(atom_idx).element().symbol().to_string()
}

/// Distance from the center of an axis-aligned box with the given half-extents
/// to its edge, along the unit direction `dir`.
///
/// Used to pull a bond's endpoint back to the edge of the atom label it meets,
/// so strokes stop at the glyph instead of running under it. Exact for a box,
/// which is a closer fit to laid-out text than a circle would be.
fn box_edge_distance(dir: Vec2, half_w: f32, half_h: f32) -> f32 {
    let tx = if dir.x.abs() > 1e-6 {
        half_w / dir.x.abs()
    } else {
        f32::INFINITY
    };
    let ty = if dir.y.abs() > 1e-6 {
        half_h / dir.y.abs()
    } else {
        f32::INFINITY
    };
    tx.min(ty)
}

/// For each bond, the screen-space centre of the ring it belongs to, or `None`
/// if it isn't in one.
///
/// This is what tells a multiple bond which side of itself is the ring
/// interior, so its second line can be drawn there rather than on an arbitrary
/// side. A bond shared between fused rings gets the centre of the *smallest*
/// ring containing it, matching how such bonds are conventionally drawn.
fn ring_bond_centers(molecule: &Molecule, transform: &Transform) -> Vec<Option<Pos2>> {
    let mut centers: Vec<Option<Pos2>> = vec![None; molecule.num_bonds()];
    let mut chosen_ring_size: Vec<usize> = vec![usize::MAX; molecule.num_bonds()];

    for ring in find_sssr(molecule) {
        let positions: Vec<Pos2> = ring
            .atoms()
            .iter()
            .filter_map(|&a| molecule.coord(a))
            .map(|p| transform.apply(p))
            .collect();
        if positions.is_empty() {
            continue;
        }

        let sum = positions
            .iter()
            .fold(Vec2::ZERO, |acc, p| acc + p.to_vec2());
        let center = (sum / positions.len() as f32).to_pos2();

        for &bond_idx in ring.bonds() {
            if ring.len() < chosen_ring_size[bond_idx] {
                chosen_ring_size[bond_idx] = ring.len();
                centers[bond_idx] = Some(center);
            }
        }
    }

    centers
}

/// Mean distance between bonded atoms in the molecule's own coordinate space,
/// or `None` if there are no bonds to measure.
///
/// Deriving a scale from a target on-screen bond length needs the source
/// length to divide by, before any transform exists.
fn mean_source_bond_length(molecule: &Molecule) -> Option<f64> {
    let mut total = 0.0;
    let mut count = 0usize;

    for bond in molecule.bonds() {
        let (a, b) = bond.atoms();
        if let (Some(pa), Some(pb)) = (molecule.coord(a), molecule.coord(b)) {
            total += pa.distance(pb);
            count += 1;
        }
    }

    if count == 0 {
        None
    } else {
        Some(total / count as f64)
    }
}

/// Mean distance between bonded atoms on screen, or `None` if there are no
/// bonds to measure. Proportional options are derived from this.
fn mean_bond_length(molecule: &Molecule, transform: &Transform) -> Option<f32> {
    let mut total = 0.0;
    let mut count = 0usize;

    for bond in molecule.bonds() {
        let (a, b) = bond.atoms();
        if let (Some(pa), Some(pb)) = (molecule.coord(a), molecule.coord(b)) {
            total += transform.apply(pa).distance(transform.apply(pb));
            count += 1;
        }
    }

    if count == 0 {
        None
    } else {
        Some(total / count as f32)
    }
}

impl Widget for StructureView<'_> {
    fn ui(self, ui: &mut Ui) -> Response {
        let (rect, response) = ui.allocate_exact_size(self.desired_size, Sense::hover());

        if !ui.is_rect_visible(rect) {
            return response;
        }

        let color = ui.visuals().text_color();
        let weak_color = ui.visuals().weak_text_color();

        let Some(coords) = self.molecule.coords() else {
            // No layout to draw. Say so rather than leaving a blank panel that
            // looks like a rendering failure.
            ui.painter().text(
                rect.center(),
                Align2::CENTER_CENTER,
                "No 2D coordinates",
                FontId::proportional(12.0),
                weak_color,
            );
            return response;
        };

        let Some(bbox) = BoundingBox::from_points(coords) else {
            // Coordinates present but no atoms — nothing to draw.
            return response;
        };

        let transform = Transform::for_options(self.molecule, bbox, rect, &self.options);

        // Derive the proportional options from the bond length this molecule
        // actually ended up with on screen. A single-atom molecule has no bonds
        // to measure, so fall back to the scale itself.
        let bond_len = mean_bond_length(self.molecule, &transform).unwrap_or(transform.scale);
        let font_size = (bond_len * self.options.font_size_ratio).clamp(
            self.options.font_size_range.0,
            self.options.font_size_range.1,
        );
        let bond_spacing = bond_len * self.options.bond_spacing_ratio;
        let font_id = FontId::proportional(font_size);
        let stroke = Stroke::new(self.options.bond_thickness, color);

        // Lay out the labels first: their sizes determine how far to pull the
        // bond endpoints back, which has to be known before any bond is drawn.
        let label_margin = font_size * self.options.label_margin_ratio;
        let label_half_extents: Vec<Option<Vec2>> = (0..self.molecule.num_atoms())
            .map(|atom_idx| {
                if !is_labeled(self.molecule, atom_idx) {
                    return None;
                }
                let galley = ui.painter().layout_no_wrap(
                    atom_label(self.molecule, atom_idx),
                    font_id.clone(),
                    color,
                );
                Some(galley.size() / 2.0 + Vec2::splat(label_margin))
            })
            .collect();

        let painter = ui.painter();

        // Which way is "into the ring" for each ring bond, so the second line
        // of a multiple bond can be drawn on that side as convention requires.
        let ring_centers = ring_bond_centers(self.molecule, &transform);

        // Bonds under labels, so the labels read cleanly on top.
        for (bond_idx, bond) in self.molecule.bonds().iter().enumerate() {
            let (a, b) = bond.atoms();
            let (Some(pa), Some(pb)) = (self.molecule.coord(a), self.molecule.coord(b)) else {
                continue;
            };

            let mut start = transform.apply(pa);
            let mut end = transform.apply(pb);

            let Some(dir) = (end - start).normalized_or_zero_check() else {
                // Superimposed atoms have no bond direction to draw along.
                // A 3D SDF flattened to 2D can produce these.
                continue;
            };

            // Pull each end back to the edge of its label, if it has one,
            // without letting the two insets cross past each other on a short
            // bond with large labels.
            let full_len = start.distance(end);
            let inset_a = label_half_extents[a]
                .map(|h| box_edge_distance(dir, h.x, h.y))
                .unwrap_or(0.0);
            let inset_b = label_half_extents[b]
                .map(|h| box_edge_distance(dir, h.x, h.y))
                .unwrap_or(0.0);
            if inset_a + inset_b >= full_len {
                continue;
            }
            start += dir * inset_a;
            end -= dir * inset_b;

            let mut offset = Vec2::new(-dir.y, dir.x) * bond_spacing;

            // On a ring bond, flip the offset if it points out of the ring, so
            // the secondary line always lands in the interior.
            let in_ring = match ring_centers.get(bond_idx).copied().flatten() {
                Some(center) => {
                    let midpoint = start + (end - start) * 0.5;
                    if offset.dot(center - midpoint) < 0.0 {
                        offset = -offset;
                    }
                    true
                }
                None => false,
            };

            // Trim applied to a secondary line drawn inside a ring, so adjacent
            // ring bonds don't run into each other at the vertices.
            let trim = (1.0 - self.options.short_bond_length).clamp(0.0, 1.0) / 2.0;

            match bond.order() {
                BondOrder::Single => {
                    painter.line_segment([start, end], stroke);
                }
                BondOrder::Double if in_ring => {
                    // Ring convention: one line on the ring perimeter, the
                    // second inside it and shortened.
                    painter.line_segment([start, end], stroke);
                    let inner_a = start + offset;
                    let inner_b = end + offset;
                    let shrink = (inner_b - inner_a) * trim;
                    painter.line_segment([inner_a + shrink, inner_b - shrink], stroke);
                }
                BondOrder::Double => {
                    // Off a ring there's no interior to favour, so both lines
                    // sit symmetrically about the bond axis.
                    painter.line_segment([start + offset * 0.5, end + offset * 0.5], stroke);
                    painter.line_segment([start - offset * 0.5, end - offset * 0.5], stroke);
                }
                BondOrder::Triple => {
                    painter.line_segment([start, end], stroke);
                    painter.line_segment([start + offset, end + offset], stroke);
                    painter.line_segment([start - offset, end - offset], stroke);
                }
                BondOrder::Aromatic => {
                    // Solid line on the ring perimeter, dashed one inside it.
                    // (The other common convention draws a single circle inside
                    // the ring instead; this keeps the per-bond form.)
                    let (solid_a, solid_b, inner_a, inner_b) = if in_ring {
                        (start, end, start + offset, end + offset)
                    } else {
                        // Aromatic bonds do occur outside a perceived ring
                        // (a lone aromatic-flagged bond, a ring the SSSR basis
                        // didn't pick); fall back to a symmetric pair.
                        (
                            start + offset * 0.5,
                            end + offset * 0.5,
                            start - offset * 0.5,
                            end - offset * 0.5,
                        )
                    };

                    painter.line_segment([solid_a, solid_b], stroke);

                    // Trimmed equally at both ends, so adjacent aromatic bonds
                    // don't collide at ring vertices.
                    let shrink = (inner_b - inner_a) * trim;
                    let dash = (bond_spacing * 1.5).max(2.0);
                    painter.extend(Shape::dashed_line(
                        &[inner_a + shrink, inner_b - shrink],
                        stroke,
                        dash,
                        dash,
                    ));
                }
                BondOrder::Quadruple => {
                    painter.line_segment([start + offset, end + offset], stroke);
                    painter.line_segment([start, end], stroke);
                    painter.line_segment([start - offset, end - offset], stroke);
                    painter.line_segment([start + offset * 2.0, end + offset * 2.0], stroke);
                }
            }
        }

        for (atom_idx, half_extent) in label_half_extents.iter().enumerate() {
            if half_extent.is_none() {
                continue;
            }
            let Some(p) = self.molecule.coord(atom_idx) else {
                continue;
            };
            painter.text(
                transform.apply(p),
                Align2::CENTER_CENTER,
                atom_label(self.molecule, atom_idx),
                font_id.clone(),
                color,
            );
        }

        response
    }
}

/// `Vec2::normalized` returns a NaN vector for a zero-length input rather than
/// signalling it, so a zero-length bond would silently produce NaN endpoints
/// and a shape egui can't tessellate. This makes that case explicit.
trait NormalizedOrZero {
    fn normalized_or_zero_check(self) -> Option<Vec2>;
}

impl NormalizedOrZero for Vec2 {
    fn normalized_or_zero_check(self) -> Option<Vec2> {
        let len = self.length();
        if len > f32::EPSILON {
            Some(self / len)
        } else {
            None
        }
    }
}

/// Draws a molecule's structure in a bordered group with a heading, fitted to
/// the available width.
pub fn structure_panel(ui: &mut Ui, molecule: &Molecule, height: f32) {
    structure_panel_with_options(ui, molecule, height, StructureOptions::default());
}

/// [`structure_panel`] with explicit rendering options.
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
    use egui::pos2;

    fn rect(w: f32, h: f32) -> Rect {
        Rect::from_min_size(pos2(0.0, 0.0), Vec2::new(w, h))
    }

    #[test]
    fn test_fit_centers_structure() {
        let bbox = BoundingBox {
            min: Point2::new(-1.0, -1.0),
            max: Point2::new(1.0, 1.0),
        };
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);

        // The bbox center maps to the rect center.
        let center = t.apply(Point2::new(0.0, 0.0));
        assert!((center.x - 50.0).abs() < 1e-4);
        assert!((center.y - 50.0).abs() < 1e-4);
    }

    #[test]
    fn test_fit_respects_padding() {
        let bbox = BoundingBox {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(1.0, 1.0),
        };
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);

        // 1 coordinate unit spans the 80px available after 10px padding a side.
        assert!((t.scale - 80.0).abs() < 1e-4);
    }

    #[test]
    fn test_fit_uses_uniform_scale() {
        // A wide, short molecule in a square rect: the constrained axis (x)
        // sets the scale for both, preserving bond angles.
        let bbox = BoundingBox {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(4.0, 1.0),
        };
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);
        assert!((t.scale - 20.0).abs() < 1e-4);
    }

    #[test]
    fn test_fit_handles_zero_height() {
        // Ethane laid out horizontally: no vertical extent at all.
        let bbox = BoundingBox {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(2.0, 0.0),
        };
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);
        assert!(t.scale.is_finite());
        assert!((t.scale - 40.0).abs() < 1e-4);
    }

    #[test]
    fn test_fit_handles_single_point() {
        // A lone atom constrains neither axis.
        let bbox = BoundingBox {
            min: Point2::new(3.0, 3.0),
            max: Point2::new(3.0, 3.0),
        };
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);
        assert!(t.scale.is_finite() && t.scale > 0.0);
        let p = t.apply(Point2::new(3.0, 3.0));
        assert!((p.x - 50.0).abs() < 1e-4);
        assert!((p.y - 50.0).abs() < 1e-4);
    }

    #[test]
    fn test_apply_flips_y() {
        let bbox = BoundingBox {
            min: Point2::new(0.0, 0.0),
            max: Point2::new(1.0, 1.0),
        };
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);

        // Chemistry +y is up; screen +y is down. A point above the center must
        // land above it on screen, i.e. at a smaller y.
        let above = t.apply(Point2::new(0.5, 1.0));
        let center = t.apply(Point2::new(0.5, 0.5));
        assert!(
            above.y < center.y,
            "higher y coordinate should render higher on screen"
        );
    }

    #[test]
    fn test_box_edge_distance() {
        // Straight along +x: the horizontal half-extent is the distance.
        let d = box_edge_distance(Vec2::new(1.0, 0.0), 5.0, 3.0);
        assert!((d - 5.0).abs() < 1e-4);

        // Straight down: the vertical half-extent.
        let d = box_edge_distance(Vec2::new(0.0, 1.0), 5.0, 3.0);
        assert!((d - 3.0).abs() < 1e-4);

        // Diagonal: whichever edge is reached first (here the shorter y).
        let inv = 1.0 / 2.0f32.sqrt();
        let d = box_edge_distance(Vec2::new(inv, inv), 5.0, 3.0);
        assert!((d - 3.0 * 2.0f32.sqrt()).abs() < 1e-3);
    }

    #[test]
    fn test_normalized_or_zero_check() {
        assert!(Vec2::new(0.0, 0.0).normalized_or_zero_check().is_none());
        let n = Vec2::new(3.0, 4.0).normalized_or_zero_check().unwrap();
        assert!((n.length() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_is_labeled_skips_carbon() {
        use chemcore::prelude::*;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));

        assert!(!is_labeled(&mol, 0), "carbon is an implicit vertex");
        assert!(is_labeled(&mol, 1), "heteroatoms are labeled");
        assert_eq!(atom_label(&mol, 1), "O");
    }

    #[test]
    fn test_mean_bond_length() {
        use chemcore::prelude::*;

        let mut mol = Molecule::new();
        let a = mol.add_atom(Atom::new(Element::carbon()));
        let b = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(a, b, BondOrder::Single)).unwrap();
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)])
            .unwrap();

        let bbox = BoundingBox::from_points(mol.coords().unwrap()).unwrap();
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);
        // One coordinate unit maps to 80px, and the bond spans exactly one.
        let len = mean_bond_length(&mol, &t).unwrap();
        assert!((len - 80.0).abs() < 1e-3);
    }

    /// Two carbons a fixed distance apart, for exercising the sizing modes.
    fn linear_molecule(coord_bond_length: f64) -> Molecule {
        use chemcore::prelude::*;

        let mut mol = Molecule::new();
        let a = mol.add_atom(Atom::new(Element::carbon()));
        let b = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(a, b, BondOrder::Single)).unwrap();
        mol.set_coords(vec![
            Point2::new(0.0, 0.0),
            Point2::new(coord_bond_length, 0.0),
        ])
        .unwrap();
        mol
    }

    fn bbox_of(mol: &Molecule) -> BoundingBox {
        BoundingBox::from_points(mol.coords().unwrap()).unwrap()
    }

    #[test]
    fn test_explicit_scale_overrides_fit() {
        let mol = linear_molecule(2.0);
        let options = StructureOptions {
            scale: 12.0,
            ..Default::default()
        };
        let t = Transform::for_options(&mol, bbox_of(&mol), rect(100.0, 100.0), &options);
        // Honoured verbatim, not recomputed to fill the rect.
        assert!((t.scale - 12.0).abs() < 1e-4);
    }

    #[test]
    fn test_bond_length_derives_scale() {
        // Source bonds are 2.0 units; asking for 30px on screen means 15px per
        // unit, regardless of how big the rect is.
        let mol = linear_molecule(2.0);
        let options = StructureOptions {
            bond_length: 30.0,
            ..Default::default()
        };
        let t = Transform::for_options(&mol, bbox_of(&mol), rect(500.0, 500.0), &options);
        assert!((t.scale - 15.0).abs() < 1e-4);

        let rendered = mean_bond_length(&mol, &t).unwrap();
        assert!((rendered - 30.0).abs() < 1e-3);
    }

    #[test]
    fn test_bond_length_is_independent_of_rect() {
        // The point of a target bond length: molecules render at a consistent
        // size across differently-shaped rects, which fit-to-rect can't do.
        let mol = linear_molecule(1.5);
        let options = StructureOptions {
            bond_length: 24.0,
            ..Default::default()
        };
        let small = Transform::for_options(&mol, bbox_of(&mol), rect(80.0, 60.0), &options);
        let large = Transform::for_options(&mol, bbox_of(&mol), rect(900.0, 700.0), &options);
        assert!((small.scale - large.scale).abs() < 1e-4);
    }

    #[test]
    fn test_explicit_scale_beats_bond_length() {
        let mol = linear_molecule(2.0);
        let options = StructureOptions {
            scale: 7.0,
            bond_length: 30.0,
            ..Default::default()
        };
        let t = Transform::for_options(&mol, bbox_of(&mol), rect(100.0, 100.0), &options);
        assert!((t.scale - 7.0).abs() < 1e-4);
    }

    #[test]
    fn test_bond_length_falls_back_to_fit_without_bonds() {
        use chemcore::prelude::*;

        // A lone atom has no bond to measure, so a target bond length can't be
        // honoured — it fits instead of dividing by zero.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.set_coords(vec![Point2::new(0.0, 0.0)]).unwrap();

        let options = StructureOptions {
            bond_length: 30.0,
            ..Default::default()
        };
        let t = Transform::for_options(&mol, bbox_of(&mol), rect(100.0, 100.0), &options);
        assert!(t.scale.is_finite() && t.scale > 0.0);
    }

    #[test]
    fn test_default_options_fit_to_rect() {
        let mol = linear_molecule(1.0);
        let options = StructureOptions::default();
        let t = Transform::for_options(&mol, bbox_of(&mol), rect(100.0, 100.0), &options);
        // Same as fit(): one unit spans the 80px left after padding.
        assert!((t.scale - 80.0).abs() < 1e-4);
    }

    fn benzene_with_coords() -> Molecule {
        use chemcore::prelude::*;

        let mut mol = Molecule::new();
        for _ in 0..6 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        for i in 0..6 {
            mol.add_bond(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic))
                .unwrap();
        }
        chemcore::layout::layout(&mut mol);
        mol
    }

    #[test]
    fn test_ring_bonds_get_a_center() {
        let mol = benzene_with_coords();
        let t = Transform::fit(bbox_of(&mol), rect(200.0, 200.0), 10.0);
        let centers = ring_bond_centers(&mol, &t);

        assert_eq!(centers.len(), mol.num_bonds());
        for (bond_idx, center) in centers.iter().enumerate() {
            assert!(
                center.is_some(),
                "ring bond {bond_idx} should have a centre"
            );
        }
    }

    #[test]
    fn test_ring_center_is_inside_the_ring() {
        let mol = benzene_with_coords();
        let t = Transform::fit(bbox_of(&mol), rect(200.0, 200.0), 10.0);
        let centers = ring_bond_centers(&mol, &t);
        let center = centers[0].unwrap();

        // The ring centre must be nearer every vertex than the vertices are to
        // the far side of the ring — i.e. genuinely interior, not off to a side.
        let ring_radius = (0..6)
            .map(|i| t.apply(mol.coord(i).unwrap()).distance(center))
            .fold(0.0f32, f32::max);
        for i in 0..6 {
            let p = t.apply(mol.coord(i).unwrap());
            assert!((p.distance(center) - ring_radius).abs() < 1.0);
        }
    }

    #[test]
    fn test_offset_flips_toward_the_ring_interior() {
        // The whole point of the fix: whichever way the perpendicular happens
        // to face, the secondary line ends up on the interior side.
        let mol = benzene_with_coords();
        let t = Transform::fit(bbox_of(&mol), rect(200.0, 200.0), 10.0);
        let centers = ring_bond_centers(&mol, &t);

        for (bond_idx, bond) in mol.bonds().iter().enumerate() {
            let (a, b) = bond.atoms();
            let start = t.apply(mol.coord(a).unwrap());
            let end = t.apply(mol.coord(b).unwrap());
            let dir = (end - start).normalized_or_zero_check().unwrap();
            let mut offset = Vec2::new(-dir.y, dir.x) * 5.0;

            let center = centers[bond_idx].unwrap();
            let midpoint = start + (end - start) * 0.5;
            if offset.dot(center - midpoint) < 0.0 {
                offset = -offset;
            }

            // After the flip, moving along the offset must get closer to the
            // ring centre than the bond midpoint is.
            assert!(
                (midpoint + offset).distance(center) < midpoint.distance(center),
                "bond {bond_idx} offset points away from the ring interior"
            );
        }
    }

    #[test]
    fn test_non_ring_bonds_have_no_center() {
        use chemcore::prelude::*;

        // Ethane: one acyclic bond, so there's no interior to favour.
        let mut mol = Molecule::new();
        let a = mol.add_atom(Atom::new(Element::carbon()));
        let b = mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(a, b, BondOrder::Double)).unwrap();
        chemcore::layout::layout(&mut mol);

        let t = Transform::fit(bbox_of(&mol), rect(200.0, 200.0), 10.0);
        let centers = ring_bond_centers(&mol, &t);
        assert_eq!(centers, vec![None]);
    }

    #[test]
    fn test_fused_bond_uses_smallest_ring() {
        use chemcore::prelude::*;

        // Naphthalene: the shared bond belongs to both 6-rings, and either
        // centre is valid; what matters is that it gets one at all rather than
        // being treated as acyclic.
        let mut mol = Molecule::new();
        for _ in 0..10 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        for &(a, b) in &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ] {
            mol.add_bond(Bond::new(a, b, BondOrder::Aromatic)).unwrap();
        }
        chemcore::layout::layout(&mut mol);

        let t = Transform::fit(bbox_of(&mol), rect(300.0, 300.0), 10.0);
        let centers = ring_bond_centers(&mol, &t);
        assert!(centers.iter().all(|c| c.is_some()));
    }

    #[test]
    fn test_mean_source_bond_length() {
        let mol = linear_molecule(2.5);
        assert!((mean_source_bond_length(&mol).unwrap() - 2.5).abs() < 1e-9);
    }

    #[test]
    fn test_mean_source_bond_length_no_bonds() {
        use chemcore::prelude::*;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.set_coords(vec![Point2::new(0.0, 0.0)]).unwrap();
        assert!(mean_source_bond_length(&mol).is_none());
    }

    #[test]
    fn test_mean_bond_length_no_bonds() {
        use chemcore::prelude::*;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.set_coords(vec![Point2::new(0.0, 0.0)]).unwrap();

        let bbox = BoundingBox::from_points(mol.coords().unwrap()).unwrap();
        let t = Transform::fit(bbox, rect(100.0, 100.0), 10.0);
        assert!(mean_bond_length(&mol, &t).is_none());
    }
}
