//! A triangulated surface: vertices, optional normals, optional colour, and
//! faces indexing into them (#313).
//!
//! OBJ and PLY both carry this, and nothing in `core` could hold it before.

/// Re-exported the same way [`crate::draw`] does, so a caller building a
/// [`Mesh`] doesn't need a second import path for the type its own
/// per-vertex colour uses.
pub use ecolor::Color32;
use thiserror::Error;

use crate::core::geometry::Point3;

/// Errors constructing a [`Mesh`].
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum MeshError {
    #[error("expected {expected} normals (one per vertex), got {got}")]
    NormalCountMismatch { expected: usize, got: usize },

    #[error("expected {expected} colors (one per vertex), got {got}")]
    ColorCountMismatch { expected: usize, got: usize },

    #[error("face {face} references vertex {index}, but there are only {vertex_count}")]
    FaceIndexOutOfRange {
        face: usize,
        index: usize,
        vertex_count: usize,
    },

    #[error("a polygon needs at least 3 vertices, got {vertex_count}")]
    DegeneratePolygon { vertex_count: usize },
}

/// A triangulated surface.
///
/// # Faces are always triangles
///
/// OBJ and PLY both permit arbitrary polygons, but this type stores only
/// triangles (`faces: Vec<[usize; 3]>`). A polygon with more than three
/// vertices is fan-triangulated on the way in
/// ([`Mesh::from_polygons`]) — vertex 0 paired with every consecutive edge —
/// and **a round trip through this crate cannot reproduce a source file's
/// original polygon structure**: a quad becomes two triangles and stays two
/// triangles. Chosen over keeping arbitrary polygons because every actual
/// consumer (a renderer, `chem convert`) wants triangles uniformly, and
/// nothing in this crate has an N-gon path to justify the alternative.
///
/// # Normals and colour are per-vertex
///
/// One entry per vertex, parallel to [`Mesh::vertices`], not OBJ's separate
/// per-face-vertex index streams (`f v/vt/vn`, which can in principle give
/// one position two different normals across faces — a hard edge).
/// Reconciling that down to this type's one-normal-per-vertex model is
/// OBJ's own story's problem (#335), not this container's.
///
/// # No topology
///
/// A mesh has no atoms, so it carries no [`crate::io::format::Carries::TOPOLOGY`]
/// — the case #310's kind-aware invariant exists for.
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    vertices: Vec<Point3>,
    normals: Option<Vec<Point3>>,
    colors: Option<Vec<Color32>>,
    faces: Vec<[usize; 3]>,
}

impl Mesh {
    /// # Errors
    /// [`MeshError::NormalCountMismatch`]/[`MeshError::ColorCountMismatch`]
    /// if a present optional array doesn't have one entry per vertex, and
    /// [`MeshError::FaceIndexOutOfRange`] if any face references a vertex
    /// index that doesn't exist.
    pub fn new(
        vertices: Vec<Point3>,
        normals: Option<Vec<Point3>>,
        colors: Option<Vec<Color32>>,
        faces: Vec<[usize; 3]>,
    ) -> Result<Self, MeshError> {
        Self::validate(&vertices, &normals, &colors, &faces)?;
        Ok(Self {
            vertices,
            normals,
            colors,
            faces,
        })
    }

    /// Builds a mesh from polygonal faces (each at least 3 vertices),
    /// fan-triangulating every one — see the type's own doc comment for the
    /// cost this accepts.
    ///
    /// ```
    /// use chem::core::prelude::*;
    ///
    /// // A unit-square quad: 0-1-2-3 in order around the perimeter.
    /// let vertices = vec![
    ///     Point3::new(0.0, 0.0, 0.0),
    ///     Point3::new(1.0, 0.0, 0.0),
    ///     Point3::new(1.0, 1.0, 0.0),
    ///     Point3::new(0.0, 1.0, 0.0),
    /// ];
    /// let mesh = Mesh::from_polygons(vertices, None, None, vec![vec![0, 1, 2, 3]]).unwrap();
    ///
    /// // Fan-triangulated from vertex 0: two triangles, not one quad.
    /// assert_eq!(mesh.faces(), &[[0, 1, 2], [0, 2, 3]]);
    /// ```
    ///
    /// # Errors
    /// [`MeshError::DegeneratePolygon`] for a polygon with fewer than 3
    /// vertices, plus everything [`Self::new`] checks against the
    /// already-triangulated result.
    pub fn from_polygons(
        vertices: Vec<Point3>,
        normals: Option<Vec<Point3>>,
        colors: Option<Vec<Color32>>,
        polygons: Vec<Vec<usize>>,
    ) -> Result<Self, MeshError> {
        let mut faces = Vec::new();
        for polygon in &polygons {
            if polygon.len() < 3 {
                return Err(MeshError::DegeneratePolygon {
                    vertex_count: polygon.len(),
                });
            }
            for i in 1..polygon.len() - 1 {
                faces.push([polygon[0], polygon[i], polygon[i + 1]]);
            }
        }
        Self::new(vertices, normals, colors, faces)
    }

    fn validate(
        vertices: &[Point3],
        normals: &Option<Vec<Point3>>,
        colors: &Option<Vec<Color32>>,
        faces: &[[usize; 3]],
    ) -> Result<(), MeshError> {
        if let Some(normals) = normals
            && normals.len() != vertices.len()
        {
            return Err(MeshError::NormalCountMismatch {
                expected: vertices.len(),
                got: normals.len(),
            });
        }
        if let Some(colors) = colors
            && colors.len() != vertices.len()
        {
            return Err(MeshError::ColorCountMismatch {
                expected: vertices.len(),
                got: colors.len(),
            });
        }
        for (face_index, face) in faces.iter().enumerate() {
            for &vertex_index in face {
                if vertex_index >= vertices.len() {
                    return Err(MeshError::FaceIndexOutOfRange {
                        face: face_index,
                        index: vertex_index,
                        vertex_count: vertices.len(),
                    });
                }
            }
        }
        Ok(())
    }

    pub fn vertices(&self) -> &[Point3] {
        &self.vertices
    }

    pub fn normals(&self) -> Option<&[Point3]> {
        self.normals.as_deref()
    }

    pub fn colors(&self) -> Option<&[Color32]> {
        self.colors.as_deref()
    }

    pub fn faces(&self) -> &[[usize; 3]] {
        &self.faces
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn num_faces(&self) -> usize {
        self.faces.len()
    }

    /// The three actual vertex positions for `face_index`.
    ///
    /// # Panics
    /// If `face_index >= self.num_faces()`.
    pub fn face_vertices(&self, face_index: usize) -> [Point3; 3] {
        let [a, b, c] = self.faces[face_index];
        [self.vertices[a], self.vertices[b], self.vertices[c]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_triangle() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ]
    }

    fn unit_quad() -> Vec<Point3> {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ]
    }

    #[test]
    fn test_a_well_formed_triangle_mesh_constructs_and_reports_its_face() {
        let vertices = unit_triangle();
        let mesh = Mesh::new(vertices.clone(), None, None, vec![[0, 1, 2]]).unwrap();
        assert_eq!(mesh.num_vertices(), 3);
        assert_eq!(mesh.num_faces(), 1);
        assert_eq!(
            mesh.face_vertices(0),
            [vertices[0], vertices[1], vertices[2]]
        );
        assert!(mesh.normals().is_none());
        assert!(mesh.colors().is_none());
    }

    #[test]
    fn test_normal_count_mismatch_is_refused() {
        let vertices = unit_triangle();
        match Mesh::new(
            vertices,
            Some(vec![Point3::new(0.0, 0.0, 1.0)]),
            None,
            vec![[0, 1, 2]],
        ) {
            Err(MeshError::NormalCountMismatch { expected, got }) => {
                assert_eq!(expected, 3);
                assert_eq!(got, 1);
            }
            other => panic!("expected NormalCountMismatch, got {other:?}"),
        }
    }

    #[test]
    fn test_color_count_mismatch_is_refused() {
        let vertices = unit_triangle();
        match Mesh::new(
            vertices,
            None,
            Some(vec![
                Color32::from_rgb(255, 0, 0),
                Color32::from_rgb(0, 255, 0),
            ]),
            vec![[0, 1, 2]],
        ) {
            Err(MeshError::ColorCountMismatch { expected, got }) => {
                assert_eq!(expected, 3);
                assert_eq!(got, 2);
            }
            other => panic!("expected ColorCountMismatch, got {other:?}"),
        }
    }

    #[test]
    fn test_face_index_out_of_range_is_refused() {
        let vertices = unit_triangle();
        match Mesh::new(vertices, None, None, vec![[0, 1, 3]]) {
            Err(MeshError::FaceIndexOutOfRange {
                face,
                index,
                vertex_count,
            }) => {
                assert_eq!(face, 0);
                assert_eq!(index, 3);
                assert_eq!(vertex_count, 3);
            }
            other => panic!("expected FaceIndexOutOfRange, got {other:?}"),
        }
    }

    #[test]
    fn test_from_polygons_triangulates_a_quad() {
        let mesh = Mesh::from_polygons(unit_quad(), None, None, vec![vec![0, 1, 2, 3]]).unwrap();
        assert_eq!(mesh.faces(), &[[0, 1, 2], [0, 2, 3]]);
    }

    #[test]
    fn test_from_polygons_passes_a_triangle_through_unchanged() {
        let mesh = Mesh::from_polygons(unit_triangle(), None, None, vec![vec![0, 1, 2]]).unwrap();
        assert_eq!(mesh.faces(), &[[0, 1, 2]]);
    }

    #[test]
    fn test_from_polygons_refuses_a_degenerate_polygon() {
        match Mesh::from_polygons(unit_triangle(), None, None, vec![vec![0, 1]]) {
            Err(MeshError::DegeneratePolygon { vertex_count }) => assert_eq!(vertex_count, 2),
            other => panic!("expected DegeneratePolygon, got {other:?}"),
        }
    }
}
