//! The unit cell and space group of a periodic structure.

use std::fmt;

use crate::core::geometry::Point3;

/// A crystallographic unit cell: three edge lengths and the three angles
/// between them.
///
/// Lengths are Ångström and angles are **degrees**, which is what every format
/// in scope writes. Trigonometry happens in radians internally and the
/// conversion is confined to this module — mixing the two is a bug that
/// survives a cubic test, because every angle is 90° there and the error
/// cancels.
///
/// # The orientation is a convention, and getting it wrong is silent
///
/// A lattice can be placed in Cartesian space in more than one way, and the
/// choices differ by a rotation. This uses the near-universal crystallographic
/// setting — **a along x, b in the xy plane, c taking up the rest** — which is
/// what PDB, CIF and essentially every toolkit assume.
///
/// Choose differently and every coordinate is wrong by a rotation: bond
/// lengths still check out, the volume still checks out, nothing throws, and
/// the structure is quietly in the wrong frame. That is why the tests pin the
/// convention directly rather than only round-tripping.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitCell {
    /// Edge lengths, Ångström.
    pub a: f64,
    pub b: f64,
    pub c: f64,
    /// Angle between **b** and **c**, degrees.
    pub alpha: f64,
    /// Angle between **a** and **c**, degrees.
    pub beta: f64,
    /// Angle between **a** and **b**, degrees.
    pub gamma: f64,
}

impl UnitCell {
    pub const fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        UnitCell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        }
    }

    /// A cube of edge `a`. Convenient, and a reminder that a cubic cell is a
    /// **bad** test fixture for anything about orientation: its off-diagonal
    /// terms are zero, so several wrong implementations pass.
    pub const fn cubic(a: f64) -> Self {
        Self::new(a, a, a, 90.0, 90.0, 90.0)
    }

    /// Whether this describes a real lattice.
    ///
    /// # Errors
    /// A message naming the offending parameter and its value. Three families
    /// of failure: a non-positive or non-finite length, an angle outside the
    /// open interval `(0, 180)`, and — the one that is not obvious — an angle
    /// *triple* that no lattice can have even though each angle is individually
    /// legal, such as 170°/170°/170°.
    pub fn validate(&self) -> Result<(), String> {
        for (name, value) in [("a", self.a), ("b", self.b), ("c", self.c)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(format!(
                    "edge {name} must be finite and positive, got {value}"
                ));
            }
        }
        for (name, value) in [
            ("alpha", self.alpha),
            ("beta", self.beta),
            ("gamma", self.gamma),
        ] {
            if !value.is_finite() || value <= 0.0 || value >= 180.0 {
                return Err(format!(
                    "angle {name} must be finite and between 0 and 180 degrees, got {value}"
                ));
            }
        }

        // Each angle can be legal on its own while the three together describe
        // no cell at all. The radicand below is what detects that, and it is
        // the whole reason validation cannot be done parameter by parameter.
        // Spelled out rather than `!(radicand > 0.0)`: clippy objects to the
        // negated comparison, and the tempting rewrite `radicand <= 0.0` would
        // quietly stop rejecting NaN, since every comparison against NaN is
        // false. The angles above are already checked finite, so NaN should be
        // unreachable — but this is the check that would have to notice.
        let radicand = self.volume_radicand();
        if !radicand.is_finite() || radicand <= 0.0 {
            return Err(format!(
                "angles {}, {} and {} describe no real cell (the volume term is {radicand})",
                self.alpha, self.beta, self.gamma
            ));
        }
        Ok(())
    }

    /// `1 - cos²α - cos²β - cos²γ + 2·cosα·cosβ·cosγ`, the term under the root
    /// in the cell volume. Positive for every real lattice.
    fn volume_radicand(&self) -> f64 {
        let (ca, cb, cg) = self.cosines();
        1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    }

    fn cosines(&self) -> (f64, f64, f64) {
        (
            self.alpha.to_radians().cos(),
            self.beta.to_radians().cos(),
            self.gamma.to_radians().cos(),
        )
    }

    /// Cell volume in cubic Ångström. `NaN` for a cell that fails
    /// [`Self::validate`].
    pub fn volume(&self) -> f64 {
        self.a * self.b * self.c * self.volume_radicand().sqrt()
    }

    /// The fractional-to-Cartesian matrix, in the setting described on the
    /// type: `m[i][j]` is axis `j`'s contribution to Cartesian component `i`.
    ///
    /// Upper-triangular by construction, which is exactly what "a along x, b in
    /// the xy plane" means — and what lets [`Self::to_fractional`] invert it by
    /// back-substitution rather than by a general matrix inverse.
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        let (ca, cb, cg) = self.cosines();
        let sg = self.gamma.to_radians().sin();
        let v = self.volume_radicand().sqrt();

        [
            [self.a, self.b * cg, self.c * cb],
            [0.0, self.b * sg, self.c * (ca - cb * cg) / sg],
            [0.0, 0.0, self.c * v / sg],
        ]
    }

    /// The **a**, **b** and **c** basis vectors in Cartesian space — the
    /// columns of [`Self::matrix`].
    pub fn basis(&self) -> [Point3; 3] {
        let m = self.matrix();
        [
            Point3::new(m[0][0], m[1][0], m[2][0]),
            Point3::new(m[0][1], m[1][1], m[2][1]),
            Point3::new(m[0][2], m[1][2], m[2][2]),
        ]
    }

    /// Fractional coordinates to Cartesian.
    pub fn to_cartesian(&self, frac: Point3) -> Point3 {
        let m = self.matrix();
        Point3::new(
            m[0][0] * frac.x + m[0][1] * frac.y + m[0][2] * frac.z,
            m[1][1] * frac.y + m[1][2] * frac.z,
            m[2][2] * frac.z,
        )
    }

    /// Cartesian coordinates to fractional.
    ///
    /// Back-substitution up the triangular matrix rather than a general
    /// inverse: exact, cheap, and it cannot quietly return nonsense for a cell
    /// that passed [`Self::validate`].
    pub fn to_fractional(&self, cart: Point3) -> Point3 {
        let m = self.matrix();
        let z = cart.z / m[2][2];
        let y = (cart.y - m[1][2] * z) / m[1][1];
        let x = (cart.x - m[0][1] * y - m[0][2] * z) / m[0][0];
        Point3::new(x, y, z)
    }
}

impl fmt::Display for UnitCell {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.4} {:.4} {:.4}  {:.3} {:.3} {:.3}",
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
    }
}

/// A space group, as the source file stated it.
///
/// **Deliberately not validated.** A CIF often gives both a number and a
/// Hermann-Mauguin symbol, and real files disagree with themselves often
/// enough that this is a case rather than a curiosity. Round-tripping what was
/// supplied is most of what correctness means here, and a converter that
/// refuses a self-inconsistent file is worse than one that reads it faithfully.
/// Checking the two against each other would also mean carrying a table of all
/// 230 groups for very little gain.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SpaceGroup {
    /// International Tables number, 1..=230 in a well-formed file.
    pub number: Option<u16>,
    /// Hermann-Mauguin symbol, as written — `P 21 21 21`, `P-1`.
    pub symbol: Option<String>,
}

impl SpaceGroup {
    pub fn from_number(number: u16) -> Self {
        SpaceGroup {
            number: Some(number),
            symbol: None,
        }
    }

    pub fn from_symbol(symbol: impl Into<String>) -> Self {
        SpaceGroup {
            number: None,
            symbol: Some(symbol.into()),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.number.is_none() && self.symbol.is_none()
    }
}

impl fmt::Display for SpaceGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (&self.symbol, self.number) {
            (Some(symbol), Some(number)) => write!(f, "{symbol} (No. {number})"),
            (Some(symbol), None) => write!(f, "{symbol}"),
            (None, Some(number)) => write!(f, "No. {number}"),
            (None, None) => write!(f, "(unknown)"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// All three angles different and none of them 90° — the only shape that
    /// exercises every off-diagonal term of the matrix.
    fn triclinic() -> UnitCell {
        UnitCell::new(9.42, 10.15, 11.32, 88.7, 79.4, 64.2)
    }

    fn angle_between(u: Point3, v: Point3) -> f64 {
        (u.dot(v) / (u.length() * v.length()))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }

    #[test]
    fn test_the_basis_reproduces_the_cell_it_came_from() {
        // The strongest check available, and the one that catches a wrong
        // orientation convention: build the basis vectors, then measure their
        // lengths and the angles between them. Those must be the six numbers
        // that went in. This tests the matrix against the *definition* of a
        // unit cell rather than against the formula used to build it.
        let cell = triclinic();
        let [va, vb, vc] = cell.basis();

        assert!((va.length() - cell.a).abs() < 1e-10, "a");
        assert!((vb.length() - cell.b).abs() < 1e-10, "b");
        assert!((vc.length() - cell.c).abs() < 1e-10, "c");

        assert!((angle_between(vb, vc) - cell.alpha).abs() < 1e-10, "alpha");
        assert!((angle_between(va, vc) - cell.beta).abs() < 1e-10, "beta");
        assert!((angle_between(va, vb) - cell.gamma).abs() < 1e-10, "gamma");
    }

    #[test]
    fn test_the_orientation_convention_is_pinned() {
        // "a along x, b in the xy plane" is not a derivable fact — it is the
        // choice this type makes, and every other toolkit makes. Asserting it
        // directly is what stops someone "simplifying" the matrix into a
        // different, equally self-consistent frame.
        let cell = triclinic();

        let along_a = cell.to_cartesian(Point3::new(1.0, 0.0, 0.0));
        assert!((along_a.x - cell.a).abs() < 1e-12);
        assert_eq!((along_a.y, along_a.z), (0.0, 0.0));

        let along_b = cell.to_cartesian(Point3::new(0.0, 1.0, 0.0));
        let expected_x = cell.b * cell.gamma.to_radians().cos();
        let expected_y = cell.b * cell.gamma.to_radians().sin();
        assert!((along_b.x - expected_x).abs() < 1e-12);
        assert!((along_b.y - expected_y).abs() < 1e-12);
        assert_eq!(along_b.z, 0.0, "b must lie in the xy plane");

        // c is the only axis with a z component.
        assert!(cell.to_cartesian(Point3::new(0.0, 0.0, 1.0)).z > 0.0);
    }

    #[test]
    fn test_fractional_and_cartesian_round_trip() {
        let cell = triclinic();
        for frac in [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.5, 0.5),
            Point3::new(0.137, 0.892, 0.041),
            Point3::new(-0.25, 1.75, 0.5),
        ] {
            let back = cell.to_fractional(cell.to_cartesian(frac));
            assert!((back.x - frac.x).abs() < 1e-12, "{frac} -> {back}");
            assert!((back.y - frac.y).abs() < 1e-12, "{frac} -> {back}");
            assert!((back.z - frac.z).abs() < 1e-12, "{frac} -> {back}");
        }
    }

    #[test]
    fn test_volume_agrees_with_the_triple_product() {
        // Two independent routes: the closed-form angle expression, and the
        // scalar triple product of the basis vectors. They must agree, and a
        // cubic cell must simply be a cubed.
        let cell = triclinic();
        let [va, vb, vc] = cell.basis();
        let triple = va.dot(vb.cross(vc)).abs();
        assert!(
            (cell.volume() - triple).abs() < 1e-9,
            "{} vs {triple}",
            cell.volume()
        );

        assert!((UnitCell::cubic(4.0).volume() - 64.0).abs() < 1e-12);
        let ortho = UnitCell::new(2.0, 3.0, 5.0, 90.0, 90.0, 90.0);
        assert!((ortho.volume() - 30.0).abs() < 1e-12);
    }

    #[test]
    fn test_validate_rejects_every_degenerate_cell() {
        let ok = triclinic();
        assert!(ok.validate().is_ok());

        let cases = [
            (
                "zero length",
                UnitCell::new(0.0, 5.0, 5.0, 90.0, 90.0, 90.0),
            ),
            (
                "negative length",
                UnitCell::new(-1.0, 5.0, 5.0, 90.0, 90.0, 90.0),
            ),
            ("angle 0", UnitCell::new(5.0, 5.0, 5.0, 0.0, 90.0, 90.0)),
            ("angle 180", UnitCell::new(5.0, 5.0, 5.0, 180.0, 90.0, 90.0)),
            (
                "nan length",
                UnitCell::new(f64::NAN, 5.0, 5.0, 90.0, 90.0, 90.0),
            ),
            (
                "nan angle",
                UnitCell::new(5.0, 5.0, 5.0, f64::NAN, 90.0, 90.0),
            ),
            // Each angle is legal on its own; together they describe nothing.
            (
                "impossible triple",
                UnitCell::new(5.0, 5.0, 5.0, 170.0, 170.0, 170.0),
            ),
        ];
        for (label, cell) in cases {
            assert!(cell.validate().is_err(), "{label} should be rejected");
        }
    }

    #[test]
    fn test_a_space_group_is_stored_exactly_as_given() {
        // Number 19 is P 21 21 21; this pairs it with P-1 on purpose. Real
        // files carry contradictions like this and a converter that refuses
        // them is less useful than one that round-trips them faithfully.
        let contradictory = SpaceGroup {
            number: Some(19),
            symbol: Some("P-1".to_string()),
        };
        assert_eq!(contradictory.number, Some(19));
        assert_eq!(contradictory.symbol.as_deref(), Some("P-1"));
        assert_eq!(contradictory.to_string(), "P-1 (No. 19)");

        assert_eq!(SpaceGroup::from_number(19).to_string(), "No. 19");
        assert_eq!(
            SpaceGroup::from_symbol("P 21 21 21").to_string(),
            "P 21 21 21"
        );
        assert_eq!(SpaceGroup::default().to_string(), "(unknown)");
        assert!(SpaceGroup::default().is_empty());
    }
}
