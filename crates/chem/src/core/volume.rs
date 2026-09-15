//! A scalar sampled on a regular 3D grid (#312).
//!
//! CCP4/MRC, DSN6, DX and CUBE all carry the same thing — a density (or other
//! scalar field) sampled on a grid — and [`crate::core::cell::UnitCell`]
//! describes a crystallographic cell, not a sampling. [`VolumeGrid`] is that
//! type: sample counts, an origin, three axis vectors and the values, kept
//! separate from any crystallographic cell a format may also state.

use thiserror::Error;

use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;

/// One of the three spatial axes a grid is sampled along.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    /// This axis's position in [`VolumeGrid`]'s own canonical order
    /// (`X` fastest-varying, then `Y`, then `Z`).
    pub fn index(&self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }
}

/// The recomputed min/max/mean/RMS of a [`VolumeGrid`]'s values.
///
/// A map header's own stated statistics are frequently stale in files
/// written by other tools, so this is never read from one — only computed,
/// by [`VolumeGrid::statistics`], from the values actually held.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeStatistics {
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub rms: f64,
}

/// Errors constructing a [`VolumeGrid`].
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum VolumeGridError {
    #[error("expected {expected} values ({dims:?}), got {got}")]
    ValueCountMismatch {
        dims: [usize; 3],
        expected: usize,
        got: usize,
    },

    /// The three axis vectors are coplanar (or zero) — a scalar triple
    /// product of (near) zero, the same volume-radicand-style check
    /// [`UnitCell::validate`] applies to a cell's angles, applied here to
    /// vectors instead.
    #[error("the three axis vectors are degenerate (span zero volume)")]
    DegenerateAxes,
}

/// A scalar sampled on a regular 3D grid.
///
/// # Canonical order
///
/// `dims`/`values` are always stored **X fastest, then Y, then Z** —
/// `values[ix + dims[0] * (iy + dims[1] * iz)]` — regardless of what order a
/// source file used. CCP4's `MAPC`/`MAPR`/`MAPS` words are the reason this
/// matters: files that permute which physical axis is columns/rows/sections
/// are common, and passing that permutation through as this type's own order
/// would put every consumer on the hook for un-permuting it themselves.
/// [`VolumeGrid::from_source_order`] is how a file's own order gets
/// canonicalized on the way in, recording what it was in
/// [`VolumeGrid::source_axis_order`].
///
/// # The grid and the cell are different things
///
/// A CCP4 map states a unit cell *and* a sampling over some sub-box of it;
/// conflating them puts the density in the wrong position for any non-cubic
/// cell. So `axes` — the Cartesian displacement of *one* grid step along
/// each axis — is independent of [`VolumeGrid::cell`], which is only ever
/// the crystallographic cell a format happened to also state, or `None` for
/// one (like DX) that states no cell at all.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeGrid {
    dims: [usize; 3],
    /// The grid's origin in Cartesian space.
    ///
    /// **Decided here, once, for whichever format needs it**: CCP4's
    /// `NXSTART`/`NYSTART`/`NZSTART` (a start index, relative to the cell)
    /// and its `ORIGIN` record (a Cartesian position) disagree in real
    /// files. A reader should prefer an explicit `ORIGIN` record when the
    /// format states one, falling back to `NXSTART * step` only when it
    /// does not — `ORIGIN` names an unambiguous Cartesian position, while
    /// `NXSTART` is stated relative to a cell that may itself be in
    /// question.
    origin: Point3,
    axes: [Point3; 3],
    values: Vec<f64>,
    cell: Option<UnitCell>,
    source_axis_order: Option<[Axis; 3]>,
}

impl VolumeGrid {
    /// # Errors
    /// [`VolumeGridError::ValueCountMismatch`] if `values.len()` isn't
    /// exactly `dims[0] * dims[1] * dims[2]`, and
    /// [`VolumeGridError::DegenerateAxes`] if `axes` spans no volume.
    pub fn new(
        dims: [usize; 3],
        origin: Point3,
        axes: [Point3; 3],
        values: Vec<f64>,
        cell: Option<UnitCell>,
    ) -> Result<Self, VolumeGridError> {
        Self::validate(dims, axes, values.len())?;
        Ok(Self {
            dims,
            origin,
            axes,
            values,
            cell,
            source_axis_order: None,
        })
    }

    /// Builds a grid from data stored in `source_order` — the file's own
    /// column/row/section axis assignment (CCP4's `MAPC`/`MAPR`/`MAPS`, or
    /// equivalent) — permuting `source_values` into this type's canonical
    /// order and recording `source_order` in [`Self::source_axis_order`].
    ///
    /// `source_order[k]` names which canonical axis the file's `k`-th
    /// stored dimension (`source_dims[k]`, `source_axes[k]`) corresponds
    /// to; `source_values` is stored fastest-varying in dimension 0, as
    /// `source_dims` itself is X-fastest in [`Self::new`].
    ///
    /// ```
    /// use chem::core::prelude::*;
    ///
    /// // A 2x1x1 grid whose file order is (Z, X, Y): the file's fastest
    /// // dimension (length 2) is actually the canonical Z axis.
    /// let grid = VolumeGrid::from_source_order(
    ///     [2, 1, 1],
    ///     [Axis::Z, Axis::X, Axis::Y],
    ///     Point3::ORIGIN,
    ///     [Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
    ///     vec![10.0, 20.0],
    ///     None,
    /// ).unwrap();
    ///
    /// // Canonical dims are (1, 1, 2): the length-2 axis landed at Z.
    /// assert_eq!(grid.dims(), [1, 1, 2]);
    /// assert_eq!(grid.value(0, 0, 0), 10.0);
    /// assert_eq!(grid.value(0, 0, 1), 20.0);
    /// ```
    ///
    /// # Errors
    /// Same as [`Self::new`], checked against the already-permuted shape.
    pub fn from_source_order(
        source_dims: [usize; 3],
        source_order: [Axis; 3],
        origin: Point3,
        source_axes: [Point3; 3],
        source_values: Vec<f64>,
        cell: Option<UnitCell>,
    ) -> Result<Self, VolumeGridError> {
        let expected = source_dims[0] * source_dims[1] * source_dims[2];
        if source_values.len() != expected {
            return Err(VolumeGridError::ValueCountMismatch {
                dims: source_dims,
                expected,
                got: source_values.len(),
            });
        }

        let mut dims = [0usize; 3];
        let mut axes = [Point3::ORIGIN; 3];
        for k in 0..3 {
            let canonical = source_order[k].index();
            dims[canonical] = source_dims[k];
            axes[canonical] = source_axes[k];
        }

        let mut values = vec![0.0; source_values.len()];
        for i2 in 0..source_dims[2] {
            for i1 in 0..source_dims[1] {
                for i0 in 0..source_dims[0] {
                    let file_index = i0 + source_dims[0] * (i1 + source_dims[1] * i2);
                    let file_coords = [i0, i1, i2];
                    let mut canonical_coords = [0usize; 3];
                    for k in 0..3 {
                        canonical_coords[source_order[k].index()] = file_coords[k];
                    }
                    let canonical_index = canonical_coords[0]
                        + dims[0] * (canonical_coords[1] + dims[1] * canonical_coords[2]);
                    values[canonical_index] = source_values[file_index];
                }
            }
        }

        Self::validate(dims, axes, values.len())?;
        Ok(Self {
            dims,
            origin,
            axes,
            values,
            cell,
            source_axis_order: Some(source_order),
        })
    }

    fn validate(
        dims: [usize; 3],
        axes: [Point3; 3],
        values_len: usize,
    ) -> Result<(), VolumeGridError> {
        let expected = dims[0] * dims[1] * dims[2];
        if values_len != expected {
            return Err(VolumeGridError::ValueCountMismatch {
                dims,
                expected,
                got: values_len,
            });
        }
        let triple = axes[0].dot(axes[1].cross(axes[2]));
        if !triple.is_finite() || triple.abs() < 1e-12 {
            return Err(VolumeGridError::DegenerateAxes);
        }
        Ok(())
    }

    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    pub fn origin(&self) -> Point3 {
        self.origin
    }

    pub fn axes(&self) -> [Point3; 3] {
        self.axes
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// The crystallographic cell, if the format also states one — separate
    /// from [`Self::axes`], see the type's own doc comment.
    pub fn cell(&self) -> Option<UnitCell> {
        self.cell
    }

    /// The file's own column/row/section axis assignment, if this grid was
    /// built from one that permuted axes — see [`Self::from_source_order`].
    pub fn source_axis_order(&self) -> Option<[Axis; 3]> {
        self.source_axis_order
    }

    /// The linear index into [`Self::values`] for grid coordinate
    /// `(ix, iy, iz)`, in canonical (X-fastest) order.
    ///
    /// # Panics
    /// If any coordinate is out of range for [`Self::dims`].
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        assert!(ix < self.dims[0] && iy < self.dims[1] && iz < self.dims[2]);
        ix + self.dims[0] * (iy + self.dims[1] * iz)
    }

    /// The value at grid coordinate `(ix, iy, iz)`.
    ///
    /// # Panics
    /// If any coordinate is out of range for [`Self::dims`].
    pub fn value(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.values[self.index(ix, iy, iz)]
    }

    /// The Cartesian position of grid coordinate `(ix, iy, iz)`.
    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> Point3 {
        self.origin + self.axes[0] * ix as f64 + self.axes[1] * iy as f64 + self.axes[2] * iz as f64
    }

    /// Recomputes min/max/mean/RMS from [`Self::values`] — never trusted
    /// from a file header, which is frequently stale. A future writer calls
    /// this to populate its own header correctly.
    ///
    /// # Panics
    /// If [`Self::values`] is empty (a grid with no samples has no
    /// statistics to report) — unreachable through [`Self::new`]/
    /// [`Self::from_source_order`], both of which reject a zero-length
    /// `dims` product only when it also disagrees with `values.len()`; an
    /// explicitly empty grid (all `dims` zero, `values` empty) is
    /// well-formed by those constructors and simply has nothing to
    /// summarize.
    pub fn statistics(&self) -> VolumeStatistics {
        let n = self.values.len() as f64;
        let min = self.values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = self
            .values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let mean = self.values.iter().sum::<f64>() / n;
        let rms = (self.values.iter().map(|v| v * v).sum::<f64>() / n).sqrt();
        VolumeStatistics {
            min,
            max,
            mean,
            rms,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_axes() -> [Point3; 3] {
        [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ]
    }

    #[test]
    fn test_a_well_formed_grid_indexes_and_positions_correctly() {
        let values = vec![
            0.0, 1.0, // z=0: (x=0,y=0), (x=1,y=0)
            2.0, 3.0, // z=0: (x=0,y=1), (x=1,y=1)
            4.0, 5.0, // z=1: (x=0,y=0), (x=1,y=0)
            6.0, 7.0, // z=1: (x=0,y=1), (x=1,y=1)
        ];
        let grid = VolumeGrid::new([2, 2, 2], Point3::ORIGIN, unit_axes(), values, None).unwrap();

        assert_eq!(grid.value(0, 0, 0), 0.0);
        assert_eq!(grid.value(1, 0, 0), 1.0);
        assert_eq!(grid.value(0, 1, 0), 2.0);
        assert_eq!(grid.value(1, 1, 1), 7.0);
        assert_eq!(grid.position(1, 1, 1), Point3::new(1.0, 1.0, 1.0));
        assert!(grid.cell().is_none());
        assert!(grid.source_axis_order().is_none());
    }

    #[test]
    fn test_value_count_mismatch_is_refused() {
        match VolumeGrid::new([2, 2, 2], Point3::ORIGIN, unit_axes(), vec![0.0; 7], None) {
            Err(VolumeGridError::ValueCountMismatch {
                dims,
                expected,
                got,
            }) => {
                assert_eq!(dims, [2, 2, 2]);
                assert_eq!(expected, 8);
                assert_eq!(got, 7);
            }
            other => panic!("expected ValueCountMismatch, got {other:?}"),
        }
    }

    #[test]
    fn test_degenerate_axes_are_refused() {
        // All three axes in the xy plane: zero volume.
        let coplanar = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        match VolumeGrid::new([1, 1, 1], Point3::ORIGIN, coplanar, vec![0.0], None) {
            Err(VolumeGridError::DegenerateAxes) => {}
            other => panic!("expected DegenerateAxes, got {other:?}"),
        }
    }

    #[test]
    fn test_from_source_order_permutes_into_canonical_order() {
        // File order is (Z, X, Y): the file's fastest-varying dimension
        // (length 3) is canonical Z; its next (length 2) is canonical X;
        // its slowest (length 1) is canonical Y. Values below are written
        // in file order (dimension 0 fastest), by hand.
        let source_dims = [3, 2, 1];
        let source_order = [Axis::Z, Axis::X, Axis::Y];
        let source_axes = [
            Point3::new(0.0, 0.0, 1.0), // steps along file dim 0 -> canonical Z
            Point3::new(1.0, 0.0, 0.0), // steps along file dim 1 -> canonical X
            Point3::new(0.0, 1.0, 0.0), // steps along file dim 2 -> canonical Y
        ];
        // file_index = f0 + 3*(f1 + 2*f2), f2 always 0 here.
        let source_values = vec![
            0.0, 1.0, 2.0, // f1=0: f0=0,1,2
            3.0, 4.0, 5.0, // f1=1: f0=0,1,2
        ];

        let grid = VolumeGrid::from_source_order(
            source_dims,
            source_order,
            Point3::ORIGIN,
            source_axes,
            source_values,
            None,
        )
        .expect("well-formed permuted fixture");

        // Canonical dims: X gets source_dims[1]=2, Y gets source_dims[2]=1,
        // Z gets source_dims[0]=3.
        assert_eq!(grid.dims(), [2, 1, 3]);
        assert_eq!(grid.source_axis_order(), Some(source_order));
        assert_eq!(
            grid.axes(),
            [source_axes[1], source_axes[2], source_axes[0]]
        );

        // Each file coordinate (f0, f1, f2=0) maps to canonical (x=f1, y=f2, z=f0).
        // file (f0=0,f1=0) = 0.0 -> canonical (x=0,y=0,z=0)
        assert_eq!(grid.value(0, 0, 0), 0.0);
        // file (f0=1,f1=0) = 1.0 -> canonical (x=0,y=0,z=1)
        assert_eq!(grid.value(0, 0, 1), 1.0);
        // file (f0=2,f1=0) = 2.0 -> canonical (x=0,y=0,z=2)
        assert_eq!(grid.value(0, 0, 2), 2.0);
        // file (f0=0,f1=1) = 3.0 -> canonical (x=1,y=0,z=0)
        assert_eq!(grid.value(1, 0, 0), 3.0);
        // file (f0=1,f1=1) = 4.0 -> canonical (x=1,y=0,z=1)
        assert_eq!(grid.value(1, 0, 1), 4.0);
        // file (f0=2,f1=1) = 5.0 -> canonical (x=1,y=0,z=2)
        assert_eq!(grid.value(1, 0, 2), 5.0);
    }

    #[test]
    fn test_statistics_are_computed_not_trusted() {
        let values = vec![1.0, 2.0, 3.0, 4.0];
        let grid = VolumeGrid::new([4, 1, 1], Point3::ORIGIN, unit_axes(), values, None).unwrap();
        let stats = grid.statistics();
        assert_eq!(stats.min, 1.0);
        assert_eq!(stats.max, 4.0);
        assert_eq!(stats.mean, 2.5);
        assert!((stats.rms - (30.0f64 / 4.0).sqrt()).abs() < 1e-12);
    }
}
