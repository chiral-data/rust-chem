//! 2D geometry primitives for structure coordinates and depiction.

use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

/// A point (or vector) in 2D space.
///
/// Used for atom coordinates — either read from a file that carries them
/// (SDF) or produced by a layout algorithm for molecules that don't (SMILES).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2 {
    pub x: f64,
    pub y: f64,
}

impl Point2 {
    pub const ORIGIN: Point2 = Point2 { x: 0.0, y: 0.0 };

    pub const fn new(x: f64, y: f64) -> Self {
        Point2 { x, y }
    }

    /// Straight-line distance to `other`.
    pub fn distance(&self, other: Point2) -> f64 {
        self.distance_squared(other).sqrt()
    }

    /// Squared distance to `other`. Cheaper than [`Self::distance`] when only
    /// comparing distances against each other or against a squared threshold,
    /// since it skips the square root.
    pub fn distance_squared(&self, other: Point2) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// The point halfway between this one and `other`.
    pub fn midpoint(&self, other: Point2) -> Point2 {
        Point2::new((self.x + other.x) / 2.0, (self.y + other.y) / 2.0)
    }

    /// Length of this point treated as a vector from the origin.
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// This vector scaled to unit length, or `None` if it has zero length
    /// (which has no meaningful direction to normalize toward).
    pub fn normalized(&self) -> Option<Point2> {
        let len = self.length();
        if len == 0.0 {
            None
        } else {
            Some(Point2::new(self.x / len, self.y / len))
        }
    }

    /// This vector rotated 90° counter-clockwise.
    ///
    /// Offsetting along a bond's perpendicular is how the parallel lines of a
    /// double or triple bond get placed either side of the bond axis.
    pub fn perpendicular(&self) -> Point2 {
        Point2::new(-self.y, self.x)
    }

    /// Angle from the positive x-axis, in radians, in `(-π, π]`.
    pub fn angle(&self) -> f64 {
        self.y.atan2(self.x)
    }

    /// This point rotated about the origin by `radians` counter-clockwise.
    pub fn rotated(&self, radians: f64) -> Point2 {
        let (sin, cos) = radians.sin_cos();
        Point2::new(self.x * cos - self.y * sin, self.x * sin + self.y * cos)
    }
}

impl Add for Point2 {
    type Output = Point2;
    fn add(self, rhs: Point2) -> Point2 {
        Point2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Sub for Point2 {
    type Output = Point2;
    fn sub(self, rhs: Point2) -> Point2 {
        Point2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f64> for Point2 {
    type Output = Point2;
    fn mul(self, rhs: f64) -> Point2 {
        Point2::new(self.x * rhs, self.y * rhs)
    }
}

impl Div<f64> for Point2 {
    type Output = Point2;
    fn div(self, rhs: f64) -> Point2 {
        Point2::new(self.x / rhs, self.y / rhs)
    }
}

impl fmt::Display for Point2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.4}, {:.4})", self.x, self.y)
    }
}

/// An axis-aligned bounding box over a set of points.
///
/// Fitting a structure into an on-screen rect needs its extent first, so the
/// renderer can derive a uniform scale and center it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub min: Point2,
    pub max: Point2,
}

impl BoundingBox {
    /// The bounding box enclosing `points`, or `None` if `points` is empty
    /// (an empty set has no extent).
    pub fn from_points<'a, I>(points: I) -> Option<Self>
    where
        I: IntoIterator<Item = &'a Point2>,
    {
        let mut iter = points.into_iter();
        let first = iter.next()?;
        let mut bbox = BoundingBox {
            min: *first,
            max: *first,
        };

        for p in iter {
            bbox.min.x = bbox.min.x.min(p.x);
            bbox.min.y = bbox.min.y.min(p.y);
            bbox.max.x = bbox.max.x.max(p.x);
            bbox.max.y = bbox.max.y.max(p.y);
        }

        Some(bbox)
    }

    pub fn width(&self) -> f64 {
        self.max.x - self.min.x
    }

    pub fn height(&self) -> f64 {
        self.max.y - self.min.y
    }

    pub fn center(&self) -> Point2 {
        self.min.midpoint(self.max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(3.0, 4.0);
        assert_eq!(a.distance(b), 5.0);
        assert_eq!(a.distance_squared(b), 25.0);
    }

    #[test]
    fn test_midpoint() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(4.0, 2.0);
        assert_eq!(a.midpoint(b), Point2::new(2.0, 1.0));
    }

    #[test]
    fn test_arithmetic() {
        let a = Point2::new(1.0, 2.0);
        let b = Point2::new(3.0, 5.0);
        assert_eq!(a + b, Point2::new(4.0, 7.0));
        assert_eq!(b - a, Point2::new(2.0, 3.0));
        assert_eq!(a * 3.0, Point2::new(3.0, 6.0));
        assert_eq!(b / 2.0, Point2::new(1.5, 2.5));
    }

    #[test]
    fn test_normalized() {
        let n = Point2::new(3.0, 4.0).normalized().unwrap();
        assert!((n.length() - 1.0).abs() < 1e-12);
        // Zero-length has no direction to normalize toward.
        assert!(Point2::ORIGIN.normalized().is_none());
    }

    #[test]
    fn test_perpendicular() {
        // 90° CCW: +x becomes +y.
        assert_eq!(Point2::new(1.0, 0.0).perpendicular(), Point2::new(0.0, 1.0));
        // Perpendicular vectors have a zero dot product.
        let v = Point2::new(2.5, -1.5);
        let p = v.perpendicular();
        assert!((v.x * p.x + v.y * p.y).abs() < 1e-12);
    }

    #[test]
    fn test_rotated() {
        let rotated = Point2::new(1.0, 0.0).rotated(std::f64::consts::FRAC_PI_2);
        assert!(rotated.distance(Point2::new(0.0, 1.0)) < 1e-12);
    }

    #[test]
    fn test_angle() {
        assert!((Point2::new(1.0, 0.0).angle() - 0.0).abs() < 1e-12);
        assert!((Point2::new(0.0, 1.0).angle() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn test_bounding_box() {
        let points = vec![
            Point2::new(1.0, 2.0),
            Point2::new(-3.0, 5.0),
            Point2::new(4.0, -1.0),
        ];
        let bbox = BoundingBox::from_points(&points).unwrap();
        assert_eq!(bbox.min, Point2::new(-3.0, -1.0));
        assert_eq!(bbox.max, Point2::new(4.0, 5.0));
        assert_eq!(bbox.width(), 7.0);
        assert_eq!(bbox.height(), 6.0);
        assert_eq!(bbox.center(), Point2::new(0.5, 2.0));
    }

    #[test]
    fn test_bounding_box_empty() {
        let points: Vec<Point2> = Vec::new();
        assert!(BoundingBox::from_points(&points).is_none());
    }

    #[test]
    fn test_bounding_box_single_point() {
        let points = vec![Point2::new(2.0, 3.0)];
        let bbox = BoundingBox::from_points(&points).unwrap();
        assert_eq!(bbox.min, bbox.max);
        assert_eq!(bbox.width(), 0.0);
        assert_eq!(bbox.height(), 0.0);
    }
}
