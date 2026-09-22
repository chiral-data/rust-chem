//! Wavefront OBJ (#335) -- plain text, `v`/`vn`/`f` lines, and the first
//! registered format with no chemistry in it at all
//! ([`crate::core::mesh::Mesh`], #313).
//!
//! **Verified against two real, independent tools, not reconstructed from
//! memory**: `trimesh` and `meshio`, both run directly against hand-built
//! fixtures.
//!
//! **Negative indices** -- a face may reference a vertex by a negative
//! offset counting backward from the most recently declared one, rather
//! than the usual 1-based forward index. Confirmed via `trimesh`: `f -4 -3
//! -2` with 5 vertices already declared resolves to 0-based `(1,2,3)`,
//! exactly `count_so_far + negative_index` per corner (checked
//! independently by hand). **`meshio` has a confirmed bug here** -- it does
//! a blind `index - 1` with no negative check at all, producing silently
//! wrong negative array indices that Python/NumPy would then misinterpret
//! as counting from the end. This module's own index-resolution helper
//! implements the semantics `trimesh` gets right, not `meshio`'s.
//!
//! **Per-corner normals vs. [`Mesh`]'s per-vertex ones**: OBJ indexes a
//! face's position (`v`) and normal (`vn`) independently per corner, so the
//! same vertex position can carry two different normals across faces (a
//! hard edge) -- something [`Mesh`] cannot represent, by its own design (see
//! its doc comment, which names this exact story as where the reconciling
//! happens). Confirmed via `trimesh`: when a 3-vertex/2-normal fixture had
//! one position referenced with both normals across two faces, the real
//! tool's output mesh had 4 vertices, not 3 -- it duplicated the position
//! once per distinct `(position, normal)` pair actually used by a face.
//! This module does the same, deduping by `(v_index, vn_index)` while
//! walking faces. `meshio` cannot represent this at all and errors on the
//! identical fixture -- a real, disclosed limitation of that tool, not
//! evidence against the design.
//!
//! **Two disclosed consequences follow directly, not bugs**: a file with no
//! `vn` lines at all (the common case) round-trips every declared vertex
//! exactly, referenced by a face or not; a file that does have `vn` lines
//! drops a vertex referenced by no face at all, since splitting is defined
//! only over what faces actually reference -- the same behavior `trimesh`
//! itself has. A file where some face corners state a normal and others
//! don't is refused outright rather than guessed at.
//!
//! **Writing** mirrors `trimesh`'s own confirmed OBJ output shape: `v x y
//! z` per vertex, `vn x y z` per vertex when normals exist, and `f a b c`
//! (or `f a//a b//b c//c` when normals exist -- index reused for vertex and
//! normal, since [`Mesh`] normals are per-vertex) with 1-based indices,
//! never negative.
//!
//! **Deliberately unread**: `vt` (texture coordinates -- parsed only far
//! enough to keep the `v/vt/vn` token stream correctly indexed, since a
//! negative `vt` offset must still resolve correctly even though the value
//! itself is discarded; [`Mesh`] has no texture-coordinate field),
//! `mtllib`/`usemtl` (a second file a `&[u8]` reader can't follow and a
//! wasm build can't open), and `o`/`g`/`s` (grouping/smoothing, no
//! [`Mesh`] slot for either) -- all recognized and silently skipped, never
//! erroring.

use std::collections::HashMap;

use crate::core::geometry::Point3;
use crate::core::mesh::Mesh;
use crate::io::errors::{ObjError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

/// Resolves a 1-based (or negative, counting backward) OBJ index against
/// how many of that kind (`v`/`vt`/`vn`) have been declared so far at this
/// point in the file, to a 0-based index -- see the module doc for why this
/// must happen during parsing, not after: a negative index's meaning
/// depends on file position, not the eventual total.
fn resolve_index(raw: i64, count_so_far: usize) -> Result<usize, ObjError> {
    if raw > 0 {
        let idx = (raw - 1) as usize;
        if idx >= count_so_far {
            return Err(ObjError::ParseError(format!(
                "index {raw} is out of range: only {count_so_far} declared so far"
            )));
        }
        Ok(idx)
    } else if raw < 0 {
        let back = (-raw) as usize;
        if back > count_so_far {
            return Err(ObjError::ParseError(format!(
                "negative index {raw} underflows: only {count_so_far} declared so far"
            )));
        }
        Ok(count_so_far - back)
    } else {
        Err(ObjError::ParseError(
            "index 0 is not valid in OBJ (indices are 1-based)".to_string(),
        ))
    }
}

fn parse_f64(tokens: &[&str], index: usize, line_kind: &str) -> Result<f64, ObjError> {
    let tok = tokens.get(index).ok_or_else(|| {
        ObjError::ParseError(format!(
            "`{line_kind}` line is missing a value at position {index}"
        ))
    })?;
    tok.parse::<f64>().map_err(|e| {
        ObjError::ParseError(format!("invalid number {tok:?} in `{line_kind}` line: {e}"))
    })
}

/// One face corner: the resolved vertex index, and the resolved normal
/// index if this corner stated one.
type Corner = (usize, Option<usize>);

/// Parses a whole OBJ file into a [`Mesh`].
pub(crate) fn parse_obj(text: &str) -> Result<Mesh, ObjError> {
    let mut v_count = 0usize;
    let mut vn_count = 0usize;
    let mut vt_count = 0usize;
    let mut raw_vertices: Vec<Point3> = Vec::new();
    let mut raw_normals: Vec<Point3> = Vec::new();
    let mut raw_faces: Vec<Vec<Corner>> = Vec::new();
    let mut any_corner_has_normal = false;
    let mut any_corner_missing_normal = false;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let Some(keyword) = tokens.next() else {
            continue;
        };
        match keyword {
            "v" => {
                let rest: Vec<&str> = tokens.collect();
                let x = parse_f64(&rest, 0, "v")?;
                let y = parse_f64(&rest, 1, "v")?;
                let z = parse_f64(&rest, 2, "v")?;
                raw_vertices.push(Point3::new(x, y, z));
                v_count += 1;
            }
            "vn" => {
                let rest: Vec<&str> = tokens.collect();
                let x = parse_f64(&rest, 0, "vn")?;
                let y = parse_f64(&rest, 1, "vn")?;
                let z = parse_f64(&rest, 2, "vn")?;
                raw_normals.push(Point3::new(x, y, z));
                vn_count += 1;
            }
            "vt" => {
                vt_count += 1;
            }
            "f" => {
                let mut corners = Vec::new();
                for tok in tokens {
                    let parts: Vec<&str> = tok.split('/').collect();
                    let v_raw: i64 = parts[0].parse().map_err(|_| {
                        ObjError::ParseError(format!("invalid face vertex index {:?}", parts[0]))
                    })?;
                    let v_idx = resolve_index(v_raw, v_count)?;

                    if let Some(vt_raw) = parts.get(1).filter(|s| !s.is_empty()) {
                        let vt_raw: i64 = vt_raw.parse().map_err(|_| {
                            ObjError::ParseError(format!("invalid face texture index {vt_raw:?}"))
                        })?;
                        resolve_index(vt_raw, vt_count)?; // validated, then discarded
                    }

                    let vn_idx = match parts.get(2).filter(|s| !s.is_empty()) {
                        Some(vn_raw) => {
                            let vn_raw: i64 = vn_raw.parse().map_err(|_| {
                                ObjError::ParseError(format!(
                                    "invalid face normal index {vn_raw:?}"
                                ))
                            })?;
                            any_corner_has_normal = true;
                            Some(resolve_index(vn_raw, vn_count)?)
                        }
                        None => {
                            any_corner_missing_normal = true;
                            None
                        }
                    };

                    corners.push((v_idx, vn_idx));
                }
                raw_faces.push(corners);
            }
            // o/g/s (grouping, smoothing) and mtllib/usemtl (materials) --
            // deliberately unread, see the module doc.
            _ => {}
        }
    }

    if any_corner_has_normal && any_corner_missing_normal {
        return Err(ObjError::ParseError(
            "some face corners state a normal and others don't in the same file".to_string(),
        ));
    }

    if any_corner_has_normal {
        let mut seen: HashMap<(usize, usize), usize> = HashMap::new();
        let mut vertices = Vec::new();
        let mut normals = Vec::new();
        let mut polygons = Vec::new();
        for corners in &raw_faces {
            let mut polygon = Vec::with_capacity(corners.len());
            for &(v_idx, vn_idx) in corners {
                let vn_idx = vn_idx.expect("checked above: every corner has a normal");
                let key = (v_idx, vn_idx);
                let new_idx = match seen.get(&key) {
                    Some(&idx) => idx,
                    None => {
                        let idx = vertices.len();
                        vertices.push(raw_vertices[v_idx]);
                        normals.push(raw_normals[vn_idx]);
                        seen.insert(key, idx);
                        idx
                    }
                };
                polygon.push(new_idx);
            }
            polygons.push(polygon);
        }
        Ok(Mesh::from_polygons(
            vertices,
            Some(normals),
            None,
            polygons,
        )?)
    } else {
        let polygons: Vec<Vec<usize>> = raw_faces
            .iter()
            .map(|corners| corners.iter().map(|&(v, _)| v).collect())
            .collect();
        Ok(Mesh::from_polygons(raw_vertices, None, None, polygons)?)
    }
}

/// Writes a [`Mesh`] as an OBJ file -- see the module doc for why this
/// exact shape (confirmed against `trimesh`'s own writer output).
pub(crate) fn write_obj(mesh: &Mesh) -> String {
    let mut out = String::new();
    out.push_str("# Written by chem\n");
    for v in mesh.vertices() {
        out.push_str(&format!("v {} {} {}\n", v.x, v.y, v.z));
    }
    if let Some(normals) = mesh.normals() {
        for n in normals {
            out.push_str(&format!("vn {} {} {}\n", n.x, n.y, n.z));
        }
    }
    let has_normals = mesh.normals().is_some();
    for face in mesh.faces() {
        let [a, b, c] = [face[0] + 1, face[1] + 1, face[2] + 1];
        if has_normals {
            out.push_str(&format!("f {a}//{a} {b}//{b} {c}//{c}\n"));
        } else {
            out.push_str(&format!("f {a} {b} {c}\n"));
        }
    }
    out
}

/// [`crate::io::format::ReadFn`] for OBJ.
pub(crate) fn read_obj_with_options(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_obj(text) {
        Ok(mesh) => out.records.push(Record {
            payload: Payload::Mesh(mesh),
            name: "Molecule_1".to_string(),
            smiles: None,
        }),
        Err(e) => out.skipped.push(Skipped {
            position: 1,
            input: String::new(),
            error: e.to_string(),
        }),
    }
    out
}

/// [`crate::io::format::ByteWriteMeshFn`] for OBJ.
pub(crate) fn write_obj_bytes(mesh: &Mesh, _options: &WriteOptions) -> Vec<u8> {
    write_obj(mesh).into_bytes()
}

/// Buffers the whole input and parses it once -- an OBJ file is always
/// exactly one mesh, the same "one record" shape
/// [`crate::io::dx::DxSupplier`] already has.
pub struct ObjSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl ObjSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match parse_obj(&text) {
                Ok(mesh) => vec![Ok(Record {
                    payload: Payload::Mesh(mesh),
                    name: "Molecule_1".to_string(),
                    smiles: None,
                })],
                Err(e) => vec![Err(ReadError::Parse {
                    position: 1,
                    message: e.to_string(),
                })],
            },
        };
        Self {
            records: records.into_iter(),
        }
    }
}

impl Iterator for ObjSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_quad_face_is_fan_triangulated_the_same_way_mesh_does() {
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 1.0 1.0 0.0
v 0.0 1.0 0.0
f 1 2 3 4
";
        let mesh = parse_obj(text).expect("valid OBJ");
        assert_eq!(mesh.faces(), &[[0, 1, 2], [0, 2, 3]]);
    }

    #[test]
    fn test_negative_indices_resolve_by_counting_backward_from_the_most_recent_vertex() {
        // Confirmed against `trimesh` directly during this story's
        // research: `f -4 -3 -2` with 5 vertices declared resolves to
        // 0-based (1,2,3).
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 1.0 1.0 0.0
v 0.0 1.0 0.0
v 0.0 0.0 1.0
f 1 2 3
f -4 -3 -2
";
        let mesh = parse_obj(text).expect("valid OBJ");
        // No `vn` anywhere -- every declared vertex round-trips, referenced
        // or not (vertex 4, the apex, is never referenced by any face).
        assert_eq!(mesh.num_vertices(), 5);
        assert_eq!(mesh.faces(), &[[0, 1, 2], [1, 2, 3]]);
    }

    #[test]
    fn test_per_corner_normals_that_differ_across_faces_split_the_vertex() {
        // Confirmed against `trimesh` directly: this exact fixture (one
        // position referenced by two different vn indices across two
        // faces) produces 4 output vertices, not 3.
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 1.0 1.0 0.0
vn 0.0 0.0 1.0
vn 1.0 0.0 0.0
f 1//1 2//1 3//1
f 1//2 2//1 3//1
";
        let mesh = parse_obj(text).expect("valid OBJ");
        assert_eq!(mesh.num_vertices(), 4);
        let normals = mesh.normals().expect("normals survive");
        assert_eq!(normals.len(), 4);
        // Position 0 (v=1) appears twice, once per distinct normal it was
        // actually paired with by a face.
        assert_eq!(mesh.vertices()[0], mesh.vertices()[3]);
        assert_eq!(normals[0], Point3::new(0.0, 0.0, 1.0));
        assert_eq!(normals[3], Point3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh.faces(), &[[0, 1, 2], [3, 1, 2]]);
    }

    #[test]
    fn test_comments_blank_lines_and_unrecognized_keywords_are_ignored() {
        let text = "\
# a comment
mtllib scene.mtl
o my_object

v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.0 1.0 0.0
g group1
s 1
usemtl material1
f 1 2 3
";
        let mesh = parse_obj(text).expect("valid OBJ");
        assert_eq!(mesh.num_vertices(), 3);
        assert_eq!(mesh.faces(), &[[0, 1, 2]]);
    }

    #[test]
    fn test_inconsistent_per_corner_normals_are_refused() {
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.0 1.0 0.0
v 1.0 1.0 0.0
vn 0.0 0.0 1.0
f 1//1 2//1 3//1
f 1 2 4
";
        let err = parse_obj(text).unwrap_err();
        assert!(matches!(err, ObjError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_an_out_of_range_index_is_refused() {
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.0 1.0 0.0
f 1 2 5
";
        let err = parse_obj(text).unwrap_err();
        assert!(matches!(err, ObjError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_zero_index_is_refused() {
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.0 1.0 0.0
f 1 2 0
";
        let err = parse_obj(text).unwrap_err();
        assert!(matches!(err, ObjError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_negative_index_that_underflows_is_refused() {
        let text = "\
v 0.0 0.0 0.0
v 1.0 0.0 0.0
f 1 -10 2
";
        let err = parse_obj(text).unwrap_err();
        assert!(matches!(err, ObjError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_writer_round_trip_with_normals() {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let normals = vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let mesh = Mesh::new(
            vertices.clone(),
            Some(normals.clone()),
            None,
            vec![[0, 1, 2], [0, 2, 3]],
        )
        .expect("valid mesh");

        let text = write_obj(&mesh);
        assert!(text.starts_with("# Written by chem\n"));
        let back = parse_obj(&text).expect("valid OBJ");

        assert_eq!(back.vertices(), &vertices[..]);
        assert_eq!(back.normals(), Some(&normals[..]));
        assert_eq!(back.faces(), mesh.faces());
    }

    #[test]
    fn test_a_writer_round_trip_without_normals() {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let mesh = Mesh::new(vertices.clone(), None, None, vec![[0, 1, 2]]).expect("valid mesh");

        let text = write_obj(&mesh);
        assert!(!text.contains("vn "));
        let back = parse_obj(&text).expect("valid OBJ");

        assert_eq!(back.vertices(), &vertices[..]);
        assert!(back.normals().is_none());
        assert_eq!(back.faces(), mesh.faces());
    }
}
