//! OpenDX's grid format (#333) -- the simplest of the four volumetric
//! formats, and a check that [`VolumeGrid`] (#312) isn't over-fitted to
//! CCP4's own header shape: no atoms, no cell, no axis permutation, no
//! endianness, no data-mode zoo.
//!
//! Confirmed against two real, independent implementations, read directly:
//! `gridData.OpenDX` (the `GridDataFormats` package, purpose-built for this
//! format) and `pdb2pqr`'s own simpler `read_dx`. Both agree on the header
//! shape: `object 1 class gridpositions counts NX NY NZ`, `origin OX OY
//! OZ`, three separate `delta` lines, `object 3 class array ... items N
//! data follows`, then `N` whitespace-separated values whose line
//! boundaries carry no meaning at all -- both real readers consume tokens
//! until they have exactly `N`, never assuming a fixed count per line.
//!
//! **Values are Z-fastest** -- confirmed beyond the issue's own claim by
//! `gridData`'s own writer comment: *"Data is written in C array order: In
//! grid\[x,y,z\] the axis z is fastest varying, then y, then finally x, i.e.
//! z is the innermost loop."* The same convention CUBE already established
//! and the same [`VolumeGrid::from_source_order`] call shape (`source_dims`
//! and `source_axes` both fed in reverse to match).
//!
//! **`delta` is a full 3-vector, never collapsed to a diagonal scalar** --
//! confirmed in both real readers, which store the complete vector on every
//! `delta` line even though most real files happen to be axis-aligned.
//! Non-axis-aligned deltas are legal; a reader that stored three scalars
//! instead of three vectors could not represent such a file.
//!
//! **No unit ambiguity** (unlike CUBE's Bohr/Å trap) -- DX states no unit
//! convention at all; every real producer (APBS) writes Å directly, and
//! every value here is treated as already Å.
//!
//! **No axis permutation** (unlike CCP4's `MAPC`/`MAPR`/`MAPS`) -- `delta`
//! line *k* always corresponds positionally to grid dimension *k*; nothing
//! in the format lets a file state a different assignment.
//!
//! **Dispatched by class-name token, not object-ID number** --
//! `gridpositions`/`gridconnections`/`array`/`field`, matching `gridData`'s
//! own real parser (object IDs are serial numbers, never relied on for
//! meaning) rather than a simpler reference's "ID happens to be 1"
//! shortcut. `gridconnections` and `field` carry no data this crate's
//! model needs and are discarded, the same "parsed and has no data-model
//! slot" treatment `io/lammps.rs`'s own discarded sections already get.
//!
//! **No `Carries::UNIT_CELL`** -- this format's own descriptor simply never
//! declares it, since [`VolumeGrid::cell`] is never referenced by either
//! direction here. That absence *is* the disclosed "a round trip through
//! DX loses the cell" the issue asks for -- no separate warning mechanism
//! needed.
//!
//! **Writing** mirrors `gridData`'s own real, VMD/Chimera/PyMOL-compatible
//! structure: the redundant `object 2 class gridconnections` statement,
//! three values per line (a real `VMD` reader requirement, per its own
//! comment), and the trailing `field`/`component` objects a strict reader
//! may expect.

use crate::core::geometry::Point3;
use crate::core::volume::{Axis, VolumeGrid};
use crate::io::errors::{DxError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

fn field<'a>(tokens: &[&'a str], index: usize) -> Result<&'a str, DxError> {
    tokens
        .get(index)
        .copied()
        .ok_or_else(|| DxError::ParseError(format!("expected a token at position {index}")))
}

fn field_f64(tokens: &[&str], index: usize) -> Result<f64, DxError> {
    let tok = field(tokens, index)?;
    tok.parse()
        .map_err(|_| DxError::ParseError(format!("invalid number: {tok:?}")))
}

fn field_usize(tokens: &[&str], index: usize) -> Result<usize, DxError> {
    let tok = field(tokens, index)?;
    tok.parse()
        .map_err(|_| DxError::ParseError(format!("invalid integer: {tok:?}")))
}

fn parse_f64_tok(tok: &str) -> Result<f64, DxError> {
    tok.parse()
        .map_err(|_| DxError::ParseError(format!("invalid number: {tok:?}")))
}

/// Parses a whole OpenDX file into a [`VolumeGrid`]. `cell` is always
/// `None` -- DX states no crystallographic cell at all.
pub(crate) fn parse_dx(text: &str) -> Result<VolumeGrid, DxError> {
    let mut dims: Option<[usize; 3]> = None;
    let mut origin: Option<[f64; 3]> = None;
    let mut deltas: Vec<[f64; 3]> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut declared_items: Option<usize> = None;

    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        match tokens[0] {
            "object" => {
                let class_pos = tokens.iter().position(|&t| t == "class").ok_or_else(|| {
                    DxError::ParseError("'object' line has no 'class' keyword".to_string())
                })?;
                let class_name = field(&tokens, class_pos + 1)?;
                match class_name {
                    "gridpositions" => {
                        let counts_pos =
                            tokens.iter().position(|&t| t == "counts").ok_or_else(|| {
                                DxError::ParseError(
                                    "'gridpositions' line has no 'counts'".to_string(),
                                )
                            })?;
                        dims = Some([
                            field_usize(&tokens, counts_pos + 1)?,
                            field_usize(&tokens, counts_pos + 2)?,
                            field_usize(&tokens, counts_pos + 3)?,
                        ]);
                    }
                    "array" => {
                        let items_pos =
                            tokens.iter().position(|&t| t == "items").ok_or_else(|| {
                                DxError::ParseError("'array' line has no 'items'".to_string())
                            })?;
                        let n = field_usize(&tokens, items_pos + 1)?;
                        declared_items = Some(n);
                        // Line boundaries carry no meaning -- consume
                        // tokens from as many subsequent lines as needed,
                        // the same "don't assume a fixed count per row"
                        // discipline both real reference readers use.
                        'values: while values.len() < n {
                            let value_line = lines.next().ok_or_else(|| {
                                DxError::ParseError(
                                    "data ends before the stated item count".to_string(),
                                )
                            })?;
                            for tok in value_line.split_whitespace() {
                                values.push(parse_f64_tok(tok)?);
                                if values.len() >= n {
                                    break 'values;
                                }
                            }
                        }
                    }
                    _ => {} // gridconnections, field: no data-model slot
                }
            }
            "origin" => {
                origin = Some([
                    field_f64(&tokens, 1)?,
                    field_f64(&tokens, 2)?,
                    field_f64(&tokens, 3)?,
                ]);
            }
            "delta" => {
                deltas.push([
                    field_f64(&tokens, 1)?,
                    field_f64(&tokens, 2)?,
                    field_f64(&tokens, 3)?,
                ]);
            }
            _ => {} // attribute/component/anything else: discarded
        }
    }

    let dims =
        dims.ok_or_else(|| DxError::ParseError("missing gridpositions counts".to_string()))?;
    let origin = origin.ok_or_else(|| DxError::ParseError("missing origin".to_string()))?;
    if deltas.len() != 3 {
        return Err(DxError::ParseError(format!(
            "expected 3 delta vectors, got {}",
            deltas.len()
        )));
    }
    let declared_items =
        declared_items.ok_or_else(|| DxError::ParseError("missing array items".to_string()))?;
    let expected = dims[0] * dims[1] * dims[2];
    if declared_items != expected || values.len() != expected {
        return Err(DxError::ParseError(format!(
            "expected {expected} values ({dims:?}), the header declared {declared_items}, got {}",
            values.len()
        )));
    }

    let origin_point = Point3::new(origin[0], origin[1], origin[2]);
    let delta_points: Vec<Point3> = deltas
        .iter()
        .map(|d| Point3::new(d[0], d[1], d[2]))
        .collect();

    let grid = VolumeGrid::from_source_order(
        [dims[2], dims[1], dims[0]],
        [Axis::Z, Axis::Y, Axis::X],
        origin_point,
        [delta_points[2], delta_points[1], delta_points[0]],
        values,
        None,
    )?;
    Ok(grid)
}

/// Writes a [`VolumeGrid`] as an OpenDX file -- the real,
/// VMD/Chimera/PyMOL-compatible structure `gridData`'s own writer produces
/// (see the module doc). `grid.cell()` is never consulted: DX cannot state
/// one.
pub(crate) fn write_dx(grid: &VolumeGrid) -> String {
    let dims = grid.dims();
    let origin = grid.origin();
    let axes = grid.axes();
    let values = grid.values();
    let n = values.len();

    let mut out = String::new();
    out.push_str("# OpenDX density file written by chem\n");
    out.push_str(&format!(
        "object 1 class gridpositions counts {} {} {}\n",
        dims[0], dims[1], dims[2]
    ));
    out.push_str(&format!(
        "origin {:.6} {:.6} {:.6}\n",
        origin.x, origin.y, origin.z
    ));
    for axis in axes {
        out.push_str(&format!(
            "delta {:.6} {:.6} {:.6}\n",
            axis.x, axis.y, axis.z
        ));
    }
    out.push_str(&format!(
        "object 2 class gridconnections counts {} {} {}\n",
        dims[0], dims[1], dims[2]
    ));
    out.push_str(&format!(
        "object 3 class array type \"double\" rank 0 items {n} data follows\n"
    ));

    let mut col = 0;
    for ix in 0..dims[0] {
        for iy in 0..dims[1] {
            for iz in 0..dims[2] {
                out.push_str(&format!("{:.6e}", grid.value(ix, iy, iz)));
                col += 1;
                if col % 3 == 0 {
                    out.push('\n');
                } else {
                    out.push(' ');
                }
            }
        }
    }
    if col % 3 != 0 {
        out.push('\n');
    }

    out.push_str("attribute \"dep\" string \"positions\"\n");
    out.push_str("object \"chem density\" class field\n");
    out.push_str("component \"positions\" value 1\n");
    out.push_str("component \"connections\" value 2\n");
    out.push_str("component \"data\" value 3\n");

    out
}

/// [`crate::io::format::ReadFn`] for DX.
pub(crate) fn read_dx_with_options(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_dx(text) {
        Ok(grid) => out.records.push(Record {
            payload: Payload::Volume(grid),
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

/// [`crate::io::format::ByteWriteVolumeFn`] for DX.
pub(crate) fn write_dx_bytes(grid: &VolumeGrid, _options: &WriteOptions) -> Vec<u8> {
    write_dx(grid).into_bytes()
}

/// Buffers the whole input and parses it once -- a DX file is always
/// exactly one grid, the same "one record" shape
/// [`crate::io::cube::CubeSupplier`]/[`crate::io::ccp4::Ccp4Supplier`]
/// already have.
pub struct DxSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl DxSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match parse_dx(&text) {
                Ok(grid) => vec![Ok(Record {
                    payload: Payload::Volume(grid),
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

impl Iterator for DxSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::cell::UnitCell;

    #[test]
    fn test_a_round_trip_with_genuinely_non_axis_aligned_deltas() {
        // The issue's own named central risk: a reader that collapsed
        // deltas to three scalars could not represent this grid at all.
        let axes = [
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(2.0, 10.0, 0.0),
            Point3::new(1.0, 1.0, 8.0),
        ];
        // Sanity: this fixture is genuinely skewed, not accidentally
        // diagonal -- otherwise the test would pass even with the bug it
        // exists to catch.
        assert!(axes[1].x != 0.0 && axes[2].x != 0.0 && axes[2].y != 0.0);

        let grid = VolumeGrid::new(
            [2, 2, 2],
            Point3::new(1.0, 2.0, 3.0),
            axes,
            (0..8).map(|i| i as f64).collect(),
            None,
        )
        .expect("valid grid");

        let text = write_dx(&grid);
        let back = parse_dx(&text).expect("valid DX");

        assert_eq!(back.dims(), grid.dims());
        let back_axes = back.axes();
        for i in 0..3 {
            assert!((back_axes[i].x - axes[i].x).abs() < 1e-4, "{i}");
            assert!((back_axes[i].y - axes[i].y).abs() < 1e-4, "{i}");
            assert!((back_axes[i].z - axes[i].z).abs() < 1e-4, "{i}");
        }
        for ix in 0..2 {
            for iy in 0..2 {
                for iz in 0..2 {
                    assert!(
                        (back.value(ix, iy, iz) - grid.value(ix, iy, iz)).abs() < 1e-6,
                        "{ix},{iy},{iz}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_values_are_read_in_the_correct_z_fastest_order_across_irregular_lines() {
        // Hand-derived, independent of this module's own code, using the
        // same z-fastest convention CUBE already established: file order
        // is X slowest, Z fastest, so values 1..6 land at
        // grid.value(0,0,0)=1, (1,0,0)=4, (0,0,1)=2, (1,0,1)=5,
        // (0,0,2)=3, (1,0,2)=6. Deliberately split across lines with an
        // irregular count (2, then 1, then 3) -- line boundaries must
        // carry no meaning, matching both real reference readers.
        let text = "\
# comment
object 1 class gridpositions counts 2 1 3
origin 0.0 0.0 0.0
delta 1.0 0.0 0.0
delta 0.0 1.0 0.0
delta 0.0 0.0 1.0
object 2 class gridconnections counts 2 1 3
object 3 class array type \"double\" rank 0 items 6 data follows
1.0 2.0
3.0
4.0 5.0 6.0
attribute \"dep\" string \"positions\"
object \"test\" class field
component \"positions\" value 1
component \"connections\" value 2
component \"data\" value 3
";
        let grid = parse_dx(text).expect("valid DX");
        assert_eq!(grid.dims(), [2, 1, 3]);
        assert_eq!(grid.value(0, 0, 0), 1.0);
        assert_eq!(grid.value(1, 0, 0), 4.0);
        assert_eq!(grid.value(0, 0, 1), 2.0);
        assert_eq!(grid.value(1, 0, 1), 5.0);
        assert_eq!(grid.value(0, 0, 2), 3.0);
        assert_eq!(grid.value(1, 0, 2), 6.0);
    }

    #[test]
    fn test_a_cell_does_not_survive_a_dx_round_trip() {
        // The mask's honesty made concrete: DX cannot state a cell at
        // all, so this loss must be real, not just declared.
        let grid = VolumeGrid::new(
            [1, 1, 1],
            Point3::ORIGIN,
            [
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            vec![1.0],
            Some(UnitCell::cubic(10.0)),
        )
        .expect("valid grid");
        assert!(grid.cell().is_some());

        let back = parse_dx(&write_dx(&grid)).expect("valid DX");
        assert!(back.cell().is_none());
    }
}
