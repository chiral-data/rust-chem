//! 2D coordinate generation for molecules that carry none.
//!
//! SMILES is a connectivity string with no geometry, so anything parsed from
//! it needs coordinates invented before it can be drawn. SDF supplies its own
//! and doesn't need this.
//!
//! The approach is the conventional one: place ring systems as regular
//! polygons, then hang acyclic substituents off them at tetrahedral-ish
//! angles, then push apart anything that landed on top of something else.
//!
//! # Limitations
//!
//! There's no force-directed refinement pass (the Kamada-Kawai step
//! smilesDrawer exposes through its `kk*` options). Ordinary drug-like
//! molecules — rings, fused rings, chains, substituents — come out fine.
//! Heavily bridged, spiro or cage systems can still overlap, because placing
//! each ring as its own regular polygon can't satisfy several rings sharing
//! atoms in three dimensions. That's a known gap, not a bug.

use crate::geometry::Point2;
use crate::molecule::Molecule;
use crate::rings::{Ring, find_sssr};
use std::collections::{HashSet, VecDeque};
use std::f64::consts::PI;

/// Target distance between bonded atoms, in coordinate units. Chosen to sit in
/// the same range as the C-C distances SDF files carry, so laid-out and
/// file-supplied molecules render at comparable scales.
pub const BOND_LENGTH: f64 = 1.5;

/// Gap between disconnected components, as a multiple of bond length.
const COMPONENT_GAP: f64 = 2.0;

/// Atoms closer than this fraction of a bond length are treated as overlapping.
/// smilesDrawer's equivalent knob is `overlapSensitivity`, defaulting to 0.42.
const OVERLAP_THRESHOLD: f64 = 0.42;

/// Number of overlap-resolution passes (smilesDrawer's
/// `overlapResolutionIterations`).
const OVERLAP_ITERATIONS: usize = 12;

/// Generates 2D coordinates, replacing any the molecule already has.
///
/// Returns `false` for an empty molecule, which has nothing to lay out.
pub fn layout(molecule: &mut Molecule) -> bool {
    if molecule.num_atoms() == 0 {
        return false;
    }

    let rings = find_sssr(molecule);
    let mut state = Layout::new(molecule, rings);
    state.run();

    let coords = state.finish();
    molecule
        .set_coords(coords)
        .expect("layout produces one coordinate per atom");
    true
}

/// Generates coordinates only if the molecule doesn't already have some, so a
/// file-supplied layout isn't discarded in favour of a generated one.
///
/// Returns `true` if the molecule has coordinates afterwards.
pub fn ensure_coords(molecule: &mut Molecule) -> bool {
    if molecule.has_coords() {
        return true;
    }
    layout(molecule)
}

struct Layout<'a> {
    molecule: &'a Molecule,
    rings: Vec<Ring>,
    coords: Vec<Option<Point2>>,
    /// Alternates the turn direction along a chain, so it zig-zags rather than
    /// spiralling consistently one way.
    flip: Vec<bool>,
}

impl<'a> Layout<'a> {
    fn new(molecule: &'a Molecule, rings: Vec<Ring>) -> Self {
        Self {
            molecule,
            rings,
            coords: vec![None; molecule.num_atoms()],
            flip: vec![false; molecule.num_atoms()],
        }
    }

    fn run(&mut self) {
        // Each connected component is laid out independently, then offset so
        // they don't share space.
        let components = self.molecule.graph().connected_components();
        let mut x_offset = 0.0;

        for component in components {
            self.layout_component(&component);

            // Shift this component clear of everything already placed.
            if x_offset > 0.0 {
                for &atom in &component {
                    if let Some(p) = self.coords[atom].as_mut() {
                        p.x += x_offset;
                    }
                }
            }
            let width = self.component_width(&component);
            x_offset += width + COMPONENT_GAP * BOND_LENGTH;
        }

        self.resolve_overlaps();
    }

    fn component_width(&self, component: &[usize]) -> f64 {
        let xs: Vec<f64> = component
            .iter()
            .filter_map(|&a| self.coords[a].map(|p| p.x))
            .collect();
        match (
            xs.iter().cloned().fold(f64::INFINITY, f64::min),
            xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        ) {
            (min, max) if min.is_finite() && max.is_finite() => max - min,
            _ => 0.0,
        }
    }

    fn layout_component(&mut self, component: &[usize]) {
        // Rings first: they're the rigid part, and everything acyclic hangs
        // off them. Doing it the other way round would mean re-deriving ring
        // geometry from wherever the chain happened to end up.
        let ring_indices: Vec<usize> = (0..self.rings.len())
            .filter(|&i| {
                self.rings[i]
                    .atoms()
                    .first()
                    .is_some_and(|a| component.contains(a))
            })
            .collect();

        self.place_rings(&ring_indices);

        // Seed acyclic-only components, which have no ring to start from.
        if !component.iter().any(|&a| self.coords[a].is_some())
            && let Some(&first) = component.first()
        {
            self.coords[first] = Some(Point2::ORIGIN);
        }

        self.place_remaining(component);
    }

    /// Places ring systems as regular polygons, fusing each subsequent ring
    /// onto whichever ring is already positioned.
    fn place_rings(&mut self, ring_indices: &[usize]) {
        let mut pending: VecDeque<usize> = ring_indices.iter().copied().collect();
        let mut stalled = 0;

        while let Some(ring_idx) = pending.pop_front() {
            let ring = self.rings[ring_idx].clone();
            let placed: Vec<usize> = ring
                .atoms()
                .iter()
                .copied()
                .filter(|&a| self.coords[a].is_some())
                .collect();

            if placed.is_empty() {
                // Nothing to attach to yet. If every remaining ring is in this
                // state we've run out of anchors, so seed one at the origin.
                if stalled > pending.len() {
                    self.place_first_ring(&ring);
                    stalled = 0;
                } else {
                    pending.push_back(ring_idx);
                    stalled += 1;
                }
                continue;
            }

            stalled = 0;
            self.place_fused_ring(&ring, &placed);
        }
    }

    /// Places a ring with no already-positioned atoms: a regular polygon
    /// centered on the origin.
    fn place_first_ring(&mut self, ring: &Ring) {
        let n = ring.len();
        let radius = circumradius(n);
        let step = 2.0 * PI / n as f64;

        for (i, &atom) in ring.atoms().iter().enumerate() {
            let angle = step * i as f64;
            self.coords[atom] = Some(Point2::new(radius * angle.cos(), radius * angle.sin()));
        }
    }

    /// Places a ring that shares atoms with one already positioned.
    fn place_fused_ring(&mut self, ring: &Ring, placed: &[usize]) {
        let n = ring.len();

        // Fused across a bond: the shared edge fixes the new ring's position
        // entirely — its center sits on the far side of that edge.
        if placed.len() >= 2
            && let Some((u, v)) = self.adjacent_pair_in_ring(ring, placed)
        {
            let pu = self.coords[u].expect("placed");
            let pv = self.coords[v].expect("placed");

            let mid = pu.midpoint(pv);
            let Some(edge_dir) = (pv - pu).normalized() else {
                return;
            };
            let perp = edge_dir.perpendicular();
            let apothem = apothem(n);

            // Two candidate centers; take the one pointing away from the
            // atoms already placed, so the new ring doesn't land on top of the
            // existing one.
            let c1 = mid + perp * apothem;
            let c2 = mid - perp * apothem;
            let reference = self.centroid_of_placed();
            let center = if c1.distance(reference) >= c2.distance(reference) {
                c1
            } else {
                c2
            };

            self.place_ring_around(ring, center, u, v);
            return;
        }

        // Spiro, or a single shared atom: pivot the polygon about it, aiming
        // away from what's already placed.
        let anchor = placed[0];
        let anchor_pos = self.coords[anchor].expect("placed");
        let reference = self.centroid_of_placed();
        let away = (anchor_pos - reference)
            .normalized()
            .unwrap_or(Point2::new(1.0, 0.0));

        let center = anchor_pos + away * circumradius(n);
        let neighbor = self.ring_neighbor(ring, anchor);
        self.place_ring_around(ring, center, anchor, neighbor);
    }

    /// Positions every unplaced atom of `ring` on a circle about `center`,
    /// walking the ring's cycle order starting from `from` toward `toward`.
    fn place_ring_around(&mut self, ring: &Ring, center: Point2, from: usize, toward: usize) {
        let n = ring.len();
        let radius = circumradius(n);
        let ordered = order_ring_from(ring, from, toward);

        let start_angle = (self.coords[from].expect("placed") - center).angle();

        // Direction of travel is whichever way puts the second atom where it
        // already is, so a fused ring joins up with its neighbour instead of
        // mirroring away from it.
        let sign = if self.coords[toward].is_some() {
            let want = (self.coords[toward].expect("placed") - center).angle();
            let delta = normalize_angle(want - start_angle);
            if delta >= 0.0 { 1.0 } else { -1.0 }
        } else {
            1.0
        };

        let step = 2.0 * PI / n as f64;
        for (i, &atom) in ordered.iter().enumerate() {
            if self.coords[atom].is_some() {
                continue;
            }
            let angle = start_angle + sign * step * i as f64;
            self.coords[atom] =
                Some(center + Point2::new(radius * angle.cos(), radius * angle.sin()));
        }
    }

    /// Finds two atoms of `ring` that are both placed and adjacent within it —
    /// a shared fusion edge.
    fn adjacent_pair_in_ring(&self, ring: &Ring, placed: &[usize]) -> Option<(usize, usize)> {
        let atoms = ring.atoms();
        let placed_set: HashSet<usize> = placed.iter().copied().collect();
        for i in 0..atoms.len() {
            let a = atoms[i];
            let b = atoms[(i + 1) % atoms.len()];
            if placed_set.contains(&a) && placed_set.contains(&b) {
                return Some((a, b));
            }
        }
        None
    }

    fn ring_neighbor(&self, ring: &Ring, atom: usize) -> usize {
        let atoms = ring.atoms();
        let pos = atoms.iter().position(|&a| a == atom).unwrap_or(0);
        atoms[(pos + 1) % atoms.len()]
    }

    fn centroid_of_placed(&self) -> Point2 {
        let placed: Vec<Point2> = self.coords.iter().flatten().copied().collect();
        if placed.is_empty() {
            return Point2::ORIGIN;
        }
        let sum = placed.iter().fold(Point2::ORIGIN, |acc, &p| acc + p);
        sum / placed.len() as f64
    }

    /// Places every atom not covered by ring placement, walking outward from
    /// what's already positioned.
    fn place_remaining(&mut self, component: &[usize]) {
        let mut queue: VecDeque<usize> = component
            .iter()
            .copied()
            .filter(|&a| self.coords[a].is_some())
            .collect();

        while let Some(current) = queue.pop_front() {
            let unplaced: Vec<usize> = self
                .molecule
                .neighbors(current)
                .iter()
                .map(|n| n.atom_idx)
                .filter(|&a| self.coords[a].is_none())
                .collect();

            if unplaced.is_empty() {
                continue;
            }

            let directions = self.substituent_directions(current, unplaced.len());
            for (&atom, dir) in unplaced.iter().zip(directions) {
                self.coords[atom] = Some(self.coords[current].expect("placed") + dir * BOND_LENGTH);
                self.flip[atom] = !self.flip[current];
                queue.push_back(atom);
            }
        }
    }

    /// Directions to fan `count` new substituents off `atom`.
    fn substituent_directions(&self, atom: usize, count: usize) -> Vec<Point2> {
        let pos = self.coords[atom].expect("placed");
        let placed_dirs: Vec<Point2> = self
            .molecule
            .neighbors(atom)
            .iter()
            .filter_map(|n| self.coords[n.atom_idx])
            .filter_map(|p| (p - pos).normalized())
            .collect();

        // Point away from everything already attached.
        let away = if placed_dirs.is_empty() {
            Point2::new(1.0, 0.0)
        } else {
            let sum = placed_dirs.iter().fold(Point2::ORIGIN, |acc, &d| acc + d);
            (sum * -1.0)
                .normalized()
                // Substituents that cancel out (a straight chain, a symmetric
                // ring atom) leave no "away" direction; go perpendicular.
                .unwrap_or_else(|| placed_dirs[0].perpendicular())
        };

        let base = away.angle();

        // With exactly one existing bond, aiming straight down `away` would
        // draw a 180° bond angle. Real chains are ~120°, so turn off the axis —
        // alternating each step, which is what makes a chain zig-zag.
        let base = if placed_dirs.len() == 1 {
            let turn = PI / 3.0;
            if self.flip[atom] {
                base + turn
            } else {
                base - turn
            }
        } else {
            base
        };

        if count == 1 {
            return vec![Point2::new(base.cos(), base.sin())];
        }

        // Several substituents on one atom: spread them evenly about `base`.
        let spread = 2.0 * PI / 3.0;
        let start = base - spread / 2.0;
        let step = spread / (count - 1) as f64;
        (0..count)
            .map(|i| {
                let a = start + step * i as f64;
                Point2::new(a.cos(), a.sin())
            })
            .collect()
    }

    /// Nudges apart atoms that ended up on top of each other.
    ///
    /// Placing each ring as its own regular polygon can't account for how
    /// several fused rings interact, so some structures land atoms in the same
    /// spot. This doesn't fix the underlying geometry — that's what a
    /// force-directed pass would do — but it stops atoms from being exactly
    /// coincident, which would otherwise render as a bond of zero length.
    fn resolve_overlaps(&mut self) {
        let threshold = BOND_LENGTH * OVERLAP_THRESHOLD;
        let n = self.coords.len();

        for _ in 0..OVERLAP_ITERATIONS {
            let mut moved = false;
            for i in 0..n {
                for j in (i + 1)..n {
                    let (Some(pi), Some(pj)) = (self.coords[i], self.coords[j]) else {
                        continue;
                    };
                    // Bonded atoms are meant to be one bond apart; only
                    // non-bonded pairs count as overlapping.
                    if self.molecule.graph().has_edge(i, j) {
                        continue;
                    }

                    let dist = pi.distance(pj);
                    if dist >= threshold {
                        continue;
                    }

                    let push = match (pj - pi).normalized() {
                        Some(dir) => dir,
                        // Exactly coincident: no direction to separate along,
                        // so pick one deterministically rather than randomly,
                        // which would make layouts unreproducible.
                        None => Point2::new(1.0, 0.0),
                    };
                    let shift = push * ((threshold - dist) / 2.0);
                    self.coords[i] = Some(pi - shift);
                    self.coords[j] = Some(pj + shift);
                    moved = true;
                }
            }
            if !moved {
                break;
            }
        }
    }

    /// Any atom that somehow escaped placement gets the origin, so the result
    /// is always one coordinate per atom as `set_coords` requires.
    fn finish(self) -> Vec<Point2> {
        self.coords
            .into_iter()
            .map(|c| c.unwrap_or(Point2::ORIGIN))
            .collect()
    }
}

/// Distance from a regular polygon's center to a vertex, for unit bond length.
fn circumradius(n: usize) -> f64 {
    if n < 3 {
        return BOND_LENGTH / 2.0;
    }
    BOND_LENGTH / (2.0 * (PI / n as f64).sin())
}

/// Distance from a regular polygon's center to an edge midpoint.
fn apothem(n: usize) -> f64 {
    if n < 3 {
        return 0.0;
    }
    BOND_LENGTH / (2.0 * (PI / n as f64).tan())
}

/// The ring's atoms rotated to start at `from`, and reversed if needed so the
/// next atom is `toward`.
fn order_ring_from(ring: &Ring, from: usize, toward: usize) -> Vec<usize> {
    let atoms = ring.atoms();
    let n = atoms.len();
    let start = atoms.iter().position(|&a| a == from).unwrap_or(0);

    let forward: Vec<usize> = (0..n).map(|i| atoms[(start + i) % n]).collect();
    if n > 1 && forward[1] == toward {
        return forward;
    }

    let mut backward = vec![atoms[start]];
    for i in 1..n {
        backward.push(atoms[(start + n - i) % n]);
    }
    backward
}

/// Wraps an angle into `(-π, π]`.
fn normalize_angle(mut angle: f64) -> f64 {
    while angle > PI {
        angle -= 2.0 * PI;
    }
    while angle <= -PI {
        angle += 2.0 * PI;
    }
    angle
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn carbon_skeleton(num_atoms: usize, bonds: &[(usize, usize)]) -> Molecule {
        let mut mol = Molecule::new();
        for _ in 0..num_atoms {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        for &(a, b) in bonds {
            mol.add_bond(Bond::new(a, b, BondOrder::Single)).unwrap();
        }
        mol
    }

    /// Every bond's length, for checking geometry came out sane.
    fn bond_lengths(mol: &Molecule) -> Vec<f64> {
        mol.bonds()
            .iter()
            .map(|b| {
                let (u, v) = b.atoms();
                mol.coord(u).unwrap().distance(mol.coord(v).unwrap())
            })
            .collect()
    }

    fn min_nonbonded_distance(mol: &Molecule) -> f64 {
        let mut min = f64::INFINITY;
        for i in 0..mol.num_atoms() {
            for j in (i + 1)..mol.num_atoms() {
                if mol.graph().has_edge(i, j) {
                    continue;
                }
                min = min.min(mol.coord(i).unwrap().distance(mol.coord(j).unwrap()));
            }
        }
        min
    }

    #[test]
    fn test_empty_molecule() {
        let mut mol = Molecule::new();
        assert!(!layout(&mut mol));
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_single_atom() {
        let mut mol = carbon_skeleton(1, &[]);
        assert!(layout(&mut mol));
        assert_eq!(mol.coords().unwrap().len(), 1);
    }

    #[test]
    fn test_every_atom_gets_a_coordinate() {
        let mut mol = carbon_skeleton(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert!(layout(&mut mol));
        assert_eq!(mol.coords().unwrap().len(), mol.num_atoms());
    }

    #[test]
    fn test_chain_bond_lengths_are_uniform() {
        let mut mol = carbon_skeleton(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]);
        layout(&mut mol);

        for len in bond_lengths(&mol) {
            assert!(
                (len - BOND_LENGTH).abs() < 1e-6,
                "bond length {len} should be {BOND_LENGTH}"
            );
        }
    }

    #[test]
    fn test_chain_zigzags_rather_than_folding_back() {
        let mut mol = carbon_skeleton(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        layout(&mut mol);

        // A chain drawn straight down the "away" direction would give 180°
        // angles; folding back on itself would give small ones. Real chains
        // sit near 120°.
        for i in 1..4 {
            let prev = mol.coord(i - 1).unwrap();
            let cur = mol.coord(i).unwrap();
            let next = mol.coord(i + 1).unwrap();
            let a = (prev - cur).normalized().unwrap();
            let b = (next - cur).normalized().unwrap();
            let angle = (a.x * b.x + a.y * b.y).clamp(-1.0, 1.0).acos().to_degrees();
            assert!(
                (100.0..=140.0).contains(&angle),
                "chain angle {angle} not near 120°"
            );
        }
    }

    #[test]
    fn test_benzene_is_a_regular_hexagon() {
        let mut mol = carbon_skeleton(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]);
        layout(&mut mol);

        for len in bond_lengths(&mol) {
            assert!((len - BOND_LENGTH).abs() < 1e-6, "ring bond {len}");
        }

        // Every vertex the same distance from the centroid is what makes it
        // regular.
        let coords = mol.coords().unwrap();
        let centroid = coords.iter().fold(Point2::ORIGIN, |a, &p| a + p) / 6.0;
        let radii: Vec<f64> = coords.iter().map(|p| p.distance(centroid)).collect();
        for r in &radii {
            assert!((r - radii[0]).abs() < 1e-6, "radius {r} vs {}", radii[0]);
        }
    }

    #[test]
    fn test_ring_atoms_are_distinct() {
        let mut mol = carbon_skeleton(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]);
        layout(&mut mol);
        assert!(min_nonbonded_distance(&mol) > BOND_LENGTH * 0.5);
    }

    #[test]
    fn test_ring_with_substituent() {
        // Toluene skeleton: the methyl must land outside the ring, not in it.
        let mut mol = carbon_skeleton(7, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6)]);
        layout(&mut mol);

        let ring_centroid = (0..6)
            .map(|i| mol.coord(i).unwrap())
            .fold(Point2::ORIGIN, |a, p| a + p)
            / 6.0;
        let methyl = mol.coord(6).unwrap();
        let attach = mol.coord(0).unwrap();

        assert!(
            methyl.distance(ring_centroid) > attach.distance(ring_centroid),
            "substituent should point away from the ring centre"
        );
    }

    #[test]
    fn test_fused_rings_share_an_edge() {
        // Naphthalene: rings 0-1-2-3-4-5 and 1-6-7-8-9-0 share the 0-1 bond.
        let mut mol = carbon_skeleton(
            10,
            &[
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
            ],
        );
        layout(&mut mol);

        // Fusing onto the shared edge should keep every bond near the target
        // length rather than stretching the second ring to reach.
        for len in bond_lengths(&mol) {
            assert!(
                (len - BOND_LENGTH).abs() < 0.3,
                "fused-ring bond length {len} strays from {BOND_LENGTH}"
            );
        }
    }

    #[test]
    fn test_fused_rings_do_not_superimpose() {
        let mut mol = carbon_skeleton(
            10,
            &[
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
            ],
        );
        layout(&mut mol);

        // The second ring must land on the far side of the shared edge; if the
        // centre choice were wrong it would fold back onto the first.
        assert!(
            min_nonbonded_distance(&mol) > BOND_LENGTH * 0.4,
            "fused rings overlap"
        );
    }

    #[test]
    fn test_disconnected_components_are_separated() {
        // Two ethanes with no bond between them.
        let mut mol = carbon_skeleton(4, &[(0, 1), (2, 3)]);
        layout(&mut mol);

        let a = mol.coord(0).unwrap().midpoint(mol.coord(1).unwrap());
        let b = mol.coord(2).unwrap().midpoint(mol.coord(3).unwrap());
        assert!(
            a.distance(b) > BOND_LENGTH,
            "components should not sit on top of each other"
        );
    }

    #[test]
    fn test_layout_is_deterministic() {
        // Same input, same output — otherwise a molecule would redraw
        // differently every time it was opened.
        let build =
            || carbon_skeleton(7, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6)]);
        let mut a = build();
        let mut b = build();
        layout(&mut a);
        layout(&mut b);
        assert_eq!(a.coords().unwrap(), b.coords().unwrap());
    }

    #[test]
    fn test_ensure_coords_preserves_existing() {
        let mut mol = carbon_skeleton(2, &[(0, 1)]);
        let supplied = vec![Point2::new(10.0, 10.0), Point2::new(11.0, 10.0)];
        mol.set_coords(supplied.clone()).unwrap();

        assert!(ensure_coords(&mut mol));
        assert_eq!(
            mol.coords().unwrap(),
            supplied.as_slice(),
            "an existing layout (e.g. from SDF) must not be replaced"
        );
    }

    #[test]
    fn test_ensure_coords_generates_when_missing() {
        let mut mol = carbon_skeleton(3, &[(0, 1), (1, 2)]);
        assert!(!mol.has_coords());
        assert!(ensure_coords(&mut mol));
        assert!(mol.has_coords());
    }

    #[test]
    fn test_polygon_geometry() {
        // A regular hexagon of unit-bond edges: circumradius equals the edge
        // length, apothem is edge * sqrt(3) / 2.
        assert!((circumradius(6) - BOND_LENGTH).abs() < 1e-9);
        assert!((apothem(6) - BOND_LENGTH * 3f64.sqrt() / 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_order_ring_from_reverses_when_needed() {
        let ring_mol = carbon_skeleton(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        let rings = find_sssr(&ring_mol);
        let ring = &rings[0];
        let atoms = ring.atoms().to_vec();

        let forward = order_ring_from(ring, atoms[0], atoms[1]);
        assert_eq!(forward[0], atoms[0]);
        assert_eq!(forward[1], atoms[1]);

        // Asking to travel the other way reverses the walk.
        let backward = order_ring_from(ring, atoms[0], atoms[3]);
        assert_eq!(backward[0], atoms[0]);
        assert_eq!(backward[1], atoms[3]);
    }

    #[test]
    fn test_normalize_angle() {
        assert!((normalize_angle(3.0 * PI) - PI).abs() < 1e-9);
        assert!((normalize_angle(-3.0 * PI) - PI).abs() < 1e-9);
        assert!((normalize_angle(0.5) - 0.5).abs() < 1e-9);
    }
}
