//! The Stanford PLY format (#336) -- second and last mesh story of Phase 5.
//! Unlike OBJ (#335), PLY has no fixed field layout at all: the header
//! *is* the schema (`element vertex N` / `property float x` / `property
//! list uchar int vertex_indices` / ...), the same class of work as
//! PRMTOP's `%FORMAT` lines. It also has three wire encodings declared in
//! that same header (`ascii`, `binary_little_endian`, `binary_big_endian`),
//! all sharing one [`crate::io::format::FormatDescriptor`].
//!
//! **Verified against two real, independent tools, not reconstructed from
//! memory** -- `plyfile` (the purpose-built Python library for exactly this
//! format) and `trimesh`, both run directly against generated and
//! hand-built fixtures.
//!
//! **Header shape and byte-exact binary layout**, confirmed by generating
//! real ASCII/little-endian/big-endian fixtures with `plyfile` and
//! hex-dumping them directly: header lines are `\n`-terminated ASCII ending
//! in a literal `end_header\n`; properties are stored in *declared order*
//! with no padding, each at its declared type's fixed width (`float` = 4
//! bytes, `double` = 8 bytes, confirmed both ways); a list property is one
//! `count_type` value followed by exactly that many `value_type` values,
//! also confirmed byte-for-byte.
//!
//! **Face vertex indices are plain 0-based integers** -- no OBJ-style
//! 1-based/negative-offset convention exists in PLY at all, confirmed
//! directly against both tools' round trips.
//!
//! **Real files tolerate property-name variation and extra properties**: a
//! hand-built fixture using `vertex_index` (singular, an older real
//! convention) instead of `vertex_indices`, plus an unrecognized extra
//! vertex property (`confidence`), was read correctly by both `plyfile` and
//! `trimesh`. This module's own reader does the same: the face element's
//! index list is found *structurally* (the first `list`-kind property
//! declared in the `face` element), never by matching a fixed name, and any
//! other declared property -- on any element -- is consumed for exactly its
//! declared width/token count and then discarded, so binary alignment never
//! drifts regardless of what a file additionally states.
//!
//! **Double-precision coordinates** (`property double x`) confirmed
//! directly via a generated binary fixture -- this module's type system is
//! generic over property width and never hardcodes `float`.
//!
//! **Writing** always emits `float` (32-bit) coordinates/normals and
//! `uchar` colours -- the common real-world convention every reference tool
//! checked here defaults to -- and a `uchar`-count/`int`-value
//! `vertex_indices` list (always exactly 3 per face, since
//! [`crate::core::mesh::Mesh`] faces are fixed-size triangles). Binary
//! little-endian is the write default; ASCII is available through
//! [`crate::io::options::PlyWriteOptions`]. This crate never writes
//! big-endian PLY -- nothing asks for it, and the reader's big-endian
//! support exists solely to read a real file that states it.
//!
//! **Deliberately unread**: any vertex property beyond `x`/`y`/`z`,
//! `nx`/`ny`/`nz`, `red`/`green`/`blue`/`alpha`, and any element other than
//! `vertex`/`face` (e.g. a real file's `edge` element) -- all consumed to
//! stay aligned, never surfaced, the same disclosed-cut posture OBJ (#335)
//! already took with `vt`.

use crate::core::geometry::Point3;
use crate::core::mesh::{Color32, Mesh};
use crate::io::errors::{PlyError, ReadError};
use crate::io::options::{PlyWriteEncoding, PlyWriteOptions, ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

/// One of PLY's eight scalar property types, generic over width -- see the
/// module doc's double-precision finding on why this can't be hardcoded to
/// `float`/`int`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PlyType {
    Int8,
    UInt8,
    Int16,
    UInt16,
    Int32,
    UInt32,
    Float32,
    Float64,
}

impl PlyType {
    fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "char" | "int8" => PlyType::Int8,
            "uchar" | "uint8" => PlyType::UInt8,
            "short" | "int16" => PlyType::Int16,
            "ushort" | "uint16" => PlyType::UInt16,
            "int" | "int32" => PlyType::Int32,
            "uint" | "uint32" => PlyType::UInt32,
            "float" | "float32" => PlyType::Float32,
            "double" | "float64" => PlyType::Float64,
            _ => return None,
        })
    }

    fn size_bytes(self) -> usize {
        match self {
            PlyType::Int8 | PlyType::UInt8 => 1,
            PlyType::Int16 | PlyType::UInt16 => 2,
            PlyType::Int32 | PlyType::UInt32 | PlyType::Float32 => 4,
            PlyType::Float64 => 8,
        }
    }

    fn parse_ascii(self, tok: &str) -> Result<f64, PlyError> {
        tok.parse::<f64>()
            .map_err(|e| PlyError::ParseError(format!("invalid number {tok:?}: {e}")))
    }

    /// `bytes` must be exactly [`Self::size_bytes`] long.
    fn read_binary(self, bytes: &[u8], big_endian: bool) -> f64 {
        match self {
            PlyType::Int8 => bytes[0] as i8 as f64,
            PlyType::UInt8 => bytes[0] as f64,
            PlyType::Int16 => {
                let b: [u8; 2] = bytes.try_into().unwrap();
                (if big_endian {
                    i16::from_be_bytes(b)
                } else {
                    i16::from_le_bytes(b)
                }) as f64
            }
            PlyType::UInt16 => {
                let b: [u8; 2] = bytes.try_into().unwrap();
                (if big_endian {
                    u16::from_be_bytes(b)
                } else {
                    u16::from_le_bytes(b)
                }) as f64
            }
            PlyType::Int32 => {
                let b: [u8; 4] = bytes.try_into().unwrap();
                (if big_endian {
                    i32::from_be_bytes(b)
                } else {
                    i32::from_le_bytes(b)
                }) as f64
            }
            PlyType::UInt32 => {
                let b: [u8; 4] = bytes.try_into().unwrap();
                (if big_endian {
                    u32::from_be_bytes(b)
                } else {
                    u32::from_le_bytes(b)
                }) as f64
            }
            PlyType::Float32 => {
                let b: [u8; 4] = bytes.try_into().unwrap();
                (if big_endian {
                    f32::from_be_bytes(b)
                } else {
                    f32::from_le_bytes(b)
                }) as f64
            }
            PlyType::Float64 => {
                let b: [u8; 8] = bytes.try_into().unwrap();
                if big_endian {
                    f64::from_be_bytes(b)
                } else {
                    f64::from_le_bytes(b)
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
enum PropertyDecl {
    Scalar {
        name: String,
        ty: PlyType,
    },
    // No `name` field: the face element's index list is found structurally
    // (the first `List`-kind property declared), never by matching a fixed
    // name -- see the module doc's real-tool finding on `vertex_index` vs
    // `vertex_indices`. The name token is still consumed while parsing the
    // header line, just not retained.
    List {
        count_ty: PlyType,
        value_ty: PlyType,
    },
}

#[derive(Debug, Clone)]
struct ElementDecl {
    name: String,
    count: usize,
    properties: Vec<PropertyDecl>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WireEncoding {
    Ascii,
    BinaryLe,
    BinaryBe,
}

/// Reads header lines only as UTF-8 -- never the whole buffer, since a
/// binary payload may not be valid UTF-8 at all. Returns the parsed
/// encoding, element/property schema, and the byte offset where the header
/// ends and the payload begins.
fn parse_header(bytes: &[u8]) -> Result<(WireEncoding, Vec<ElementDecl>, usize), PlyError> {
    let mut offset = 0usize;
    let mut first_line = true;
    let mut encoding: Option<WireEncoding> = None;
    let mut elements: Vec<ElementDecl> = Vec::new();

    loop {
        let newline_pos = bytes[offset..]
            .iter()
            .position(|&b| b == b'\n')
            .ok_or_else(|| {
                PlyError::ParseError("header truncated: missing end_header".to_string())
            })?;
        let line_bytes = &bytes[offset..offset + newline_pos];
        let line = std::str::from_utf8(line_bytes)
            .map_err(|_| PlyError::ParseError("header line is not valid UTF-8".to_string()))?
            .trim_end_matches('\r')
            .trim();
        offset += newline_pos + 1;

        if first_line {
            first_line = false;
            if line != "ply" {
                return Err(PlyError::ParseError(format!(
                    "expected \"ply\" as the first line, got {line:?}"
                )));
            }
            continue;
        }

        let mut tokens = line.split_whitespace();
        let Some(keyword) = tokens.next() else {
            continue;
        };

        match keyword {
            "format" => {
                let enc = tokens.next().ok_or_else(|| {
                    PlyError::ParseError("format line is missing an encoding".to_string())
                })?;
                encoding = Some(match enc {
                    "ascii" => WireEncoding::Ascii,
                    "binary_little_endian" => WireEncoding::BinaryLe,
                    "binary_big_endian" => WireEncoding::BinaryBe,
                    other => {
                        return Err(PlyError::ParseError(format!(
                            "unrecognized PLY encoding {other:?}"
                        )));
                    }
                });
            }
            "comment" | "obj_info" => {}
            "element" => {
                let name = tokens.next().ok_or_else(|| {
                    PlyError::ParseError("element line is missing a name".to_string())
                })?;
                let count: usize = tokens
                    .next()
                    .ok_or_else(|| {
                        PlyError::ParseError("element line is missing a count".to_string())
                    })?
                    .parse()
                    .map_err(|_| {
                        PlyError::ParseError("element count is not a valid integer".to_string())
                    })?;
                elements.push(ElementDecl {
                    name: name.to_string(),
                    count,
                    properties: Vec::new(),
                });
            }
            "property" => {
                let element = elements.last_mut().ok_or_else(|| {
                    PlyError::ParseError("property line appears before any element".to_string())
                })?;
                let first = tokens.next().ok_or_else(|| {
                    PlyError::ParseError("property line is missing a type".to_string())
                })?;
                if first == "list" {
                    let count_ty_name = tokens.next().ok_or_else(|| {
                        PlyError::ParseError("list property is missing a count type".to_string())
                    })?;
                    let value_ty_name = tokens.next().ok_or_else(|| {
                        PlyError::ParseError("list property is missing a value type".to_string())
                    })?;
                    // Consumed (a list property's name is part of the
                    // line's own grammar) but not retained -- see
                    // `PropertyDecl::List`'s own doc comment.
                    tokens.next().ok_or_else(|| {
                        PlyError::ParseError("list property is missing a name".to_string())
                    })?;
                    let count_ty = PlyType::from_name(count_ty_name).ok_or_else(|| {
                        PlyError::ParseError(format!(
                            "unrecognized property type {count_ty_name:?}"
                        ))
                    })?;
                    let value_ty = PlyType::from_name(value_ty_name).ok_or_else(|| {
                        PlyError::ParseError(format!(
                            "unrecognized property type {value_ty_name:?}"
                        ))
                    })?;
                    element
                        .properties
                        .push(PropertyDecl::List { count_ty, value_ty });
                } else {
                    let name = tokens.next().ok_or_else(|| {
                        PlyError::ParseError("property line is missing a name".to_string())
                    })?;
                    let ty = PlyType::from_name(first).ok_or_else(|| {
                        PlyError::ParseError(format!("unrecognized property type {first:?}"))
                    })?;
                    element.properties.push(PropertyDecl::Scalar {
                        name: name.to_string(),
                        ty,
                    });
                }
            }
            "end_header" => {
                let encoding = encoding
                    .ok_or_else(|| PlyError::ParseError("header has no format line".to_string()))?;
                return Ok((encoding, elements, offset));
            }
            other => {
                return Err(PlyError::ParseError(format!(
                    "unrecognized header keyword {other:?}"
                )));
            }
        }
    }
}

/// One already-read property value -- a plain number, or a list's values.
enum PropValue {
    Scalar(f64),
    List(Vec<f64>),
}

/// Reads property values uniformly whether the payload is ASCII tokens or
/// binary bytes -- the one place this module's ASCII/binary duality is
/// resolved, so the element-walking logic above it never branches on
/// encoding at all.
enum Source<'a> {
    Ascii(std::str::SplitWhitespace<'a>),
    Binary {
        bytes: &'a [u8],
        cursor: usize,
        big_endian: bool,
    },
}

impl Source<'_> {
    fn read_scalar(&mut self, ty: PlyType) -> Result<f64, PlyError> {
        match self {
            Source::Ascii(tokens) => {
                let tok = tokens.next().ok_or_else(|| {
                    PlyError::ParseError("unexpected end of PLY data".to_string())
                })?;
                ty.parse_ascii(tok)
            }
            Source::Binary {
                bytes,
                cursor,
                big_endian,
            } => {
                let size = ty.size_bytes();
                let slice = bytes
                    .get(*cursor..*cursor + size)
                    .ok_or_else(|| PlyError::ParseError("PLY data truncated".to_string()))?;
                let v = ty.read_binary(slice, *big_endian);
                *cursor += size;
                Ok(v)
            }
        }
    }

    fn read_property(&mut self, decl: &PropertyDecl) -> Result<PropValue, PlyError> {
        match decl {
            PropertyDecl::Scalar { ty, .. } => Ok(PropValue::Scalar(self.read_scalar(*ty)?)),
            PropertyDecl::List {
                count_ty, value_ty, ..
            } => {
                let count = self.read_scalar(*count_ty)?.max(0.0) as usize;
                let mut values = Vec::with_capacity(count);
                for _ in 0..count {
                    values.push(self.read_scalar(*value_ty)?);
                }
                Ok(PropValue::List(values))
            }
        }
    }
}

fn scalar_at(values: &[PropValue], index: usize) -> f64 {
    match &values[index] {
        PropValue::Scalar(v) => *v,
        PropValue::List(_) => unreachable!("recognized scalar properties are never declared list"),
    }
}

/// Parses a whole PLY file (any of its three wire encodings) into a
/// [`Mesh`].
pub(crate) fn parse_ply(bytes: &[u8]) -> Result<Mesh, PlyError> {
    let (encoding, elements, data_offset) = parse_header(bytes)?;

    let mut source = match encoding {
        WireEncoding::Ascii => {
            let text = std::str::from_utf8(&bytes[data_offset..])
                .map_err(|_| PlyError::ParseError("PLY payload is not valid UTF-8".to_string()))?;
            Source::Ascii(text.split_whitespace())
        }
        WireEncoding::BinaryLe => Source::Binary {
            bytes,
            cursor: data_offset,
            big_endian: false,
        },
        WireEncoding::BinaryBe => Source::Binary {
            bytes,
            cursor: data_offset,
            big_endian: true,
        },
    };

    let mut vertices: Vec<Point3> = Vec::new();
    let mut normals: Vec<Point3> = Vec::new();
    let mut any_normals = false;
    let mut colors: Vec<Color32> = Vec::new();
    let mut any_colors = false;
    let mut polygons: Vec<Vec<usize>> = Vec::new();

    for element in &elements {
        let is_vertex = element.name == "vertex";
        let is_face = element.name == "face";

        let idx = |target: &str| {
            element
                .properties
                .iter()
                .position(|p| matches!(p, PropertyDecl::Scalar { name, .. } if name == target))
        };
        let (ix, iy, iz) = (idx("x"), idx("y"), idx("z"));
        let (inx, iny, inz) = (idx("nx"), idx("ny"), idx("nz"));
        let (ir, ig, ib, ia) = (idx("red"), idx("green"), idx("blue"), idx("alpha"));
        let has_normals = inx.is_some() && iny.is_some() && inz.is_some();
        let partial_normals = (inx.is_some() || iny.is_some() || inz.is_some()) && !has_normals;
        let has_colors = ir.is_some() && ig.is_some() && ib.is_some();
        let partial_colors = (ir.is_some() || ig.is_some() || ib.is_some()) && !has_colors;

        if is_vertex {
            if ix.is_none() || iy.is_none() || iz.is_none() {
                return Err(PlyError::ParseError(
                    "vertex element is missing x/y/z".to_string(),
                ));
            }
            if partial_normals {
                return Err(PlyError::ParseError(
                    "vertex element has some but not all of nx/ny/nz".to_string(),
                ));
            }
            if partial_colors {
                return Err(PlyError::ParseError(
                    "vertex element has some but not all of red/green/blue".to_string(),
                ));
            }
        }

        let list_index = if is_face {
            let pos = element
                .properties
                .iter()
                .position(|p| matches!(p, PropertyDecl::List { .. }));
            if pos.is_none() {
                return Err(PlyError::ParseError(
                    "face element has no list property to use as vertex indices".to_string(),
                ));
            }
            pos
        } else {
            None
        };

        for _ in 0..element.count {
            let mut values: Vec<PropValue> = Vec::with_capacity(element.properties.len());
            for prop in &element.properties {
                values.push(source.read_property(prop)?);
            }

            if is_vertex {
                vertices.push(Point3::new(
                    scalar_at(&values, ix.unwrap()),
                    scalar_at(&values, iy.unwrap()),
                    scalar_at(&values, iz.unwrap()),
                ));
                if has_normals {
                    normals.push(Point3::new(
                        scalar_at(&values, inx.unwrap()),
                        scalar_at(&values, iny.unwrap()),
                        scalar_at(&values, inz.unwrap()),
                    ));
                    any_normals = true;
                }
                if has_colors {
                    let r = scalar_at(&values, ir.unwrap()) as u8;
                    let g = scalar_at(&values, ig.unwrap()) as u8;
                    let b = scalar_at(&values, ib.unwrap()) as u8;
                    let color = match ia {
                        Some(ia) => {
                            let a = scalar_at(&values, ia) as u8;
                            Color32::from_rgba_unmultiplied(r, g, b, a)
                        }
                        None => Color32::from_rgb(r, g, b),
                    };
                    colors.push(color);
                    any_colors = true;
                }
            } else if is_face && let PropValue::List(list_values) = &values[list_index.unwrap()] {
                let mut polygon = Vec::with_capacity(list_values.len());
                for &v in list_values {
                    if v < 0.0 || v.fract() != 0.0 {
                        return Err(PlyError::ParseError(format!(
                            "face vertex index {v} is not a non-negative integer"
                        )));
                    }
                    polygon.push(v as usize);
                }
                polygons.push(polygon);
            }
            // Any other element (e.g. a real file's `edge`) is fully
            // consumed above and discarded here -- see the module doc.
        }
    }

    let normals = any_normals.then_some(normals);
    let colors = any_colors.then_some(colors);
    Ok(Mesh::from_polygons(vertices, normals, colors, polygons)?)
}

/// Writes a [`Mesh`] as a PLY file in the requested wire encoding -- see the
/// module doc for the exact property set and why writing never produces
/// big-endian.
pub(crate) fn write_ply(mesh: &Mesh, options: &PlyWriteOptions) -> Vec<u8> {
    let ascii = matches!(options.encoding, PlyWriteEncoding::Ascii);
    let normals = mesh.normals();
    let colors = mesh.colors();

    let mut header = String::new();
    header.push_str("ply\n");
    header.push_str(if ascii {
        "format ascii 1.0\n"
    } else {
        "format binary_little_endian 1.0\n"
    });
    header.push_str("comment Written by chem\n");
    header.push_str(&format!("element vertex {}\n", mesh.num_vertices()));
    header.push_str("property float x\nproperty float y\nproperty float z\n");
    if normals.is_some() {
        header.push_str("property float nx\nproperty float ny\nproperty float nz\n");
    }
    if colors.is_some() {
        header.push_str("property uchar red\nproperty uchar green\nproperty uchar blue\n");
    }
    header.push_str(&format!("element face {}\n", mesh.num_faces()));
    header.push_str("property list uchar int vertex_indices\n");
    header.push_str("end_header\n");

    let mut out = header.into_bytes();

    if ascii {
        for (i, v) in mesh.vertices().iter().enumerate() {
            let mut line = format!("{} {} {}", v.x, v.y, v.z);
            if let Some(normals) = normals {
                let n = normals[i];
                line.push_str(&format!(" {} {} {}", n.x, n.y, n.z));
            }
            if let Some(colors) = colors {
                let c = colors[i];
                line.push_str(&format!(" {} {} {}", c.r(), c.g(), c.b()));
            }
            line.push('\n');
            out.extend_from_slice(line.as_bytes());
        }
        for face in mesh.faces() {
            out.extend_from_slice(format!("3 {} {} {}\n", face[0], face[1], face[2]).as_bytes());
        }
    } else {
        for (i, v) in mesh.vertices().iter().enumerate() {
            out.extend_from_slice(&(v.x as f32).to_le_bytes());
            out.extend_from_slice(&(v.y as f32).to_le_bytes());
            out.extend_from_slice(&(v.z as f32).to_le_bytes());
            if let Some(normals) = normals {
                let n = normals[i];
                out.extend_from_slice(&(n.x as f32).to_le_bytes());
                out.extend_from_slice(&(n.y as f32).to_le_bytes());
                out.extend_from_slice(&(n.z as f32).to_le_bytes());
            }
            if let Some(colors) = colors {
                let c = colors[i];
                out.push(c.r());
                out.push(c.g());
                out.push(c.b());
            }
        }
        for face in mesh.faces() {
            out.push(3u8);
            for &idx in face {
                out.extend_from_slice(&(idx as i32).to_le_bytes());
            }
        }
    }

    out
}

/// [`crate::io::format::ByteReadFn`] for PLY.
pub(crate) fn read_ply_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_ply(bytes) {
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

/// [`crate::io::format::ByteWriteMeshFn`] for PLY -- narrows the shared
/// [`WriteOptions`] down to [`PlyWriteOptions`] before calling
/// [`write_ply`], the same shim shape `format.rs`'s own `write_sdf_records`
/// already uses for SDF.
pub(crate) fn write_ply_bytes(mesh: &Mesh, options: &WriteOptions) -> Vec<u8> {
    write_ply(mesh, &options.ply)
}

/// Buffers the whole input and parses it once -- a PLY file is always
/// exactly one mesh, the same "one record" shape
/// [`crate::io::obj::ObjSupplier`] already has.
pub struct PlySupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl PlySupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match std::io::Read::read_to_end(&mut reader, &mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match parse_ply(&bytes) {
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

impl Iterator for PlySupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_mesh_with_normals_and_colors() -> Mesh {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let normals = vec![Point3::new(0.0, 0.0, 1.0); 4];
        let colors = vec![
            Color32::from_rgb(255, 0, 0),
            Color32::from_rgb(0, 255, 0),
            Color32::from_rgb(0, 0, 255),
            Color32::from_rgb(255, 255, 0),
        ];
        Mesh::new(
            vertices,
            Some(normals),
            Some(colors),
            vec![[0, 1, 2], [0, 2, 3]],
        )
        .expect("valid mesh")
    }

    #[test]
    fn test_an_ascii_round_trip_is_exact() {
        let mesh = sample_mesh_with_normals_and_colors();
        let bytes = write_ply(
            &mesh,
            &PlyWriteOptions {
                encoding: PlyWriteEncoding::Ascii,
            },
        );
        assert!(bytes.starts_with(b"ply\nformat ascii 1.0\n"));
        let back = parse_ply(&bytes).expect("valid PLY");

        // f64's Display/FromStr round-trips exactly, so ASCII needs no
        // tolerance despite the header nominally saying "float".
        assert_eq!(back.vertices(), mesh.vertices());
        assert_eq!(back.normals(), mesh.normals());
        assert_eq!(back.colors(), mesh.colors());
        assert_eq!(back.faces(), mesh.faces());
    }

    #[test]
    fn test_a_binary_little_endian_round_trip_is_the_write_default() {
        let mesh = sample_mesh_with_normals_and_colors();
        let bytes = write_ply(&mesh, &PlyWriteOptions::default());
        assert!(bytes.starts_with(b"ply\nformat binary_little_endian 1.0\n"));
        let back = parse_ply(&bytes).expect("valid PLY");

        // Coordinates round-trip through f32 on the binary path; every
        // fixture value here is exactly representable in f32, so equality
        // holds without a tolerance.
        assert_eq!(back.vertices(), mesh.vertices());
        assert_eq!(back.normals(), mesh.normals());
        assert_eq!(back.colors(), mesh.colors());
        assert_eq!(back.faces(), mesh.faces());
    }

    #[test]
    fn test_a_real_plyfile_produced_big_endian_fixture_reads_correctly() {
        // Genuine bytes produced by Python's `plyfile` library (the
        // purpose-built reference for this format) during this story's
        // research -- not this crate's own writer, which never emits
        // big-endian. 3 vertices with normals+colors, 1 triangle face.
        #[rustfmt::skip]
        let bytes: &[u8] = &[
            0x70, 0x6c, 0x79, 0x0a, 0x66, 0x6f, 0x72, 0x6d, 0x61, 0x74, 0x20, 0x62, 0x69, 0x6e, 0x61, 0x72,
            0x79, 0x5f, 0x62, 0x69, 0x67, 0x5f, 0x65, 0x6e, 0x64, 0x69, 0x61, 0x6e, 0x20, 0x31, 0x2e, 0x30,
            0x0a, 0x65, 0x6c, 0x65, 0x6d, 0x65, 0x6e, 0x74, 0x20, 0x76, 0x65, 0x72, 0x74, 0x65, 0x78, 0x20,
            0x33, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x66, 0x6c, 0x6f, 0x61, 0x74,
            0x20, 0x78, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x66, 0x6c, 0x6f, 0x61,
            0x74, 0x20, 0x79, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x66, 0x6c, 0x6f,
            0x61, 0x74, 0x20, 0x7a, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x66, 0x6c,
            0x6f, 0x61, 0x74, 0x20, 0x6e, 0x78, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20,
            0x66, 0x6c, 0x6f, 0x61, 0x74, 0x20, 0x6e, 0x79, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74,
            0x79, 0x20, 0x66, 0x6c, 0x6f, 0x61, 0x74, 0x20, 0x6e, 0x7a, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65,
            0x72, 0x74, 0x79, 0x20, 0x75, 0x63, 0x68, 0x61, 0x72, 0x20, 0x72, 0x65, 0x64, 0x0a, 0x70, 0x72,
            0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x75, 0x63, 0x68, 0x61, 0x72, 0x20, 0x67, 0x72, 0x65,
            0x65, 0x6e, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x75, 0x63, 0x68, 0x61,
            0x72, 0x20, 0x62, 0x6c, 0x75, 0x65, 0x0a, 0x65, 0x6c, 0x65, 0x6d, 0x65, 0x6e, 0x74, 0x20, 0x66,
            0x61, 0x63, 0x65, 0x20, 0x31, 0x0a, 0x70, 0x72, 0x6f, 0x70, 0x65, 0x72, 0x74, 0x79, 0x20, 0x6c,
            0x69, 0x73, 0x74, 0x20, 0x75, 0x63, 0x68, 0x61, 0x72, 0x20, 0x69, 0x6e, 0x74, 0x20, 0x76, 0x65,
            0x72, 0x74, 0x65, 0x78, 0x5f, 0x69, 0x6e, 0x64, 0x69, 0x63, 0x65, 0x73, 0x0a, 0x65, 0x6e, 0x64,
            0x5f, 0x68, 0x65, 0x61, 0x64, 0x65, 0x72, 0x0a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0x80, 0x00, 0x00,
            0xff, 0x00, 0x00, 0x3f, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0x80, 0x00, 0x00, 0x00, 0xff, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x3f, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x3f, 0x80, 0x00, 0x00, 0x00, 0x00, 0xff, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x01, 0x00, 0x00, 0x00, 0x02,
        ];

        let mesh = parse_ply(bytes).expect("valid PLY");
        assert_eq!(mesh.num_vertices(), 3);
        assert_eq!(mesh.vertices()[0], Point3::new(0.0, 0.0, 0.0));
        assert_eq!(mesh.vertices()[1], Point3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh.vertices()[2], Point3::new(0.0, 1.0, 0.0));
        let normals = mesh.normals().expect("normals present");
        assert_eq!(normals[0], Point3::new(0.0, 0.0, 1.0));
        let colors = mesh.colors().expect("colors present");
        assert_eq!(colors[0], Color32::from_rgb(255, 0, 0));
        assert_eq!(colors[1], Color32::from_rgb(0, 255, 0));
        assert_eq!(colors[2], Color32::from_rgb(0, 0, 255));
        assert_eq!(mesh.faces(), &[[0, 1, 2]]);
    }

    #[test]
    fn test_an_alternate_face_property_name_and_an_extra_vertex_property_are_tolerated() {
        // Confirmed against both `plyfile` and `trimesh` directly: real
        // readers accept `vertex_index` (singular) as the face's index
        // list, and silently skip an extra declared property they don't
        // recognize (here, `confidence`) without losing alignment.
        let text = "\
ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
property float confidence
element face 1
property list uchar int vertex_index
end_header
0.0 0.0 0.0 0.9
1.0 0.0 0.0 0.8
0.0 1.0 0.0 0.7
3 0 1 2
";
        let mesh = parse_ply(text.as_bytes()).expect("valid PLY");
        assert_eq!(mesh.num_vertices(), 3);
        assert_eq!(mesh.vertices()[1], Point3::new(1.0, 0.0, 0.0));
        assert_eq!(mesh.faces(), &[[0, 1, 2]]);
    }

    #[test]
    fn test_an_unrecognized_element_between_vertex_and_face_is_skipped_without_losing_alignment() {
        let text = "\
ply
format ascii 1.0
element vertex 2
property float x
property float y
property float z
element edge 1
property int vertex1
property int vertex2
element face 1
property list uchar int vertex_indices
end_header
0.0 0.0 0.0
1.0 0.0 0.0
0 1
3 0 1 0
";
        let mesh = parse_ply(text.as_bytes()).expect("valid PLY");
        assert_eq!(mesh.num_vertices(), 2);
        assert_eq!(mesh.faces(), &[[0, 1, 0]]);
    }

    #[test]
    fn test_double_precision_coordinates_round_trip_exactly() {
        let mut header = String::from(
            "ply\nformat binary_little_endian 1.0\nelement vertex 1\nproperty double x\nproperty double y\nproperty double z\nelement face 0\nproperty list uchar int vertex_indices\nend_header\n",
        );
        let mut bytes = std::mem::take(&mut header).into_bytes();
        bytes.extend_from_slice(&0.123_456_789_012_345_f64.to_le_bytes());
        bytes.extend_from_slice(&1.0_f64.to_le_bytes());
        bytes.extend_from_slice(&2.0_f64.to_le_bytes());

        let mesh = parse_ply(&bytes).expect("valid PLY");
        assert_eq!(mesh.num_vertices(), 1);
        assert_eq!(
            mesh.vertices()[0],
            Point3::new(0.123_456_789_012_345, 1.0, 2.0)
        );
    }

    #[test]
    fn test_a_face_less_file_is_a_valid_point_cloud() {
        let text = "\
ply
format ascii 1.0
element vertex 2
property float x
property float y
property float z
end_header
0.0 0.0 0.0
1.0 1.0 1.0
";
        let mesh = parse_ply(text.as_bytes()).expect("valid PLY");
        assert_eq!(mesh.num_vertices(), 2);
        assert_eq!(mesh.num_faces(), 0);
    }

    #[test]
    fn test_a_vertex_element_missing_x_is_refused() {
        let text = "\
ply
format ascii 1.0
element vertex 1
property float y
property float z
end_header
0.0 0.0
";
        let err = parse_ply(text.as_bytes()).unwrap_err();
        assert!(matches!(err, PlyError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_partial_normal_triplet_is_refused() {
        let text = "\
ply
format ascii 1.0
element vertex 1
property float x
property float y
property float z
property float nx
property float ny
end_header
0.0 0.0 0.0 0.0 0.0
";
        let err = parse_ply(text.as_bytes()).unwrap_err();
        assert!(matches!(err, PlyError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_face_element_with_no_list_property_is_refused() {
        let text = "\
ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
element face 1
property int material_id
end_header
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
5
";
        let err = parse_ply(text.as_bytes()).unwrap_err();
        assert!(matches!(err, PlyError::ParseError(_)), "{err}");
    }
}
