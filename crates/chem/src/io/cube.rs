//! Gaussian's CUBE format (#331) — a scalar sampled on a regular 3D grid,
//! plus (optionally) the atoms that produced it. The first format ever
//! registered as [`crate::io::format::Kind::Volume`], and the first to
//! populate [`VolumeGrid::atoms`] — reusing the same bonds-free `Molecule`
//! shape XYZ/GRO already use, rather than a second, parallel atoms-only
//! type nothing else understands.
//!
//! Confirmed directly against two real, independent implementations, not
//! reconstructed from prose: `ase.io.cube` and `pdb2pqr`'s own
//! `io.py::write_cube` (both read locally from their installed packages).
//!
//! **Layout**: two free-text comment lines (skipped, never parsed for
//! metadata); `natoms origin_x origin_y origin_z [num_val]`; three axis
//! rows, each `n_i x_i y_i z_i`; `natoms` atom lines
//! (`atomic_number charge x y z`, `charge` discarded — no data-model slot,
//! and typically `0.0` in real files); then every value, whitespace-
//! separated, **written with the third axis (Z) varying fastest** —
//! `ase.io.cube`'s own `"OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"`
//! convention. This is exactly the shape [`VolumeGrid::from_source_order`]
//! (#312) exists to canonicalize: `source_dims = [nz, ny, nx]` (the file's
//! fastest dimension first), `source_order = [Axis::Z, Axis::Y, Axis::X]`,
//! fed the raw value stream completely unrearranged.
//!
//! **Two independent sign conventions, both confirmed against real code**:
//! - `natoms < 0` means an extra orbital-index record precedes the values,
//!   and the values themselves carry an extra per-point dimension — a
//!   shape [`VolumeGrid`] does not model (one scalar per point, not N).
//!   Refused as [`crate::io::errors::CubeError::MultipleValuesUnsupported`]
//!   rather than half-supported, the same reasoning as a positive `natoms`
//!   with a stated `num_val != 1` on the header line (an independent way
//!   the same file shape states "more than one value per point").
//! - A negative axis-row sample count means that axis's step vector — and,
//!   in every real writer, the shared origin — is already Å, not Bohr.
//!   Confirmed concretely in `pdb2pqr`'s real `write_cube`, which
//!   deliberately writes a negated count when its source data needs no
//!   Bohr conversion. `ase.io.cube`'s own vendored *reader* does not
//!   actually implement reading this back correctly (it feeds the raw
//!   signed value straight into a `reshape` call, which would fail) — a
//!   real gap in that reference, not something replicated here. This
//!   reader treats the sign as one whole-file flag (any axis negative ⇒
//!   the whole header is Å), matching how every real writer emits it.
//!
//! **Writing** always uses Bohr and positive counts — the least-ambiguous
//! convention, the same posture every prior format's writer in this crate
//! already took (NCTRAJ always writing `"angstrom"`, LAMMPS Trajectory
//! always writing unscaled coordinates).

use crate::core::atom::{Atom, Element};
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::units::BOHR_TO_ANGSTROM;
use crate::core::volume::{Axis, VolumeGrid};
use crate::io::errors::{CubeError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

fn next_tokens<'a>(lines: &mut std::str::Lines<'a>) -> Result<Vec<&'a str>, CubeError> {
    let line = lines
        .next()
        .ok_or_else(|| CubeError::ParseError("unexpected end of file".to_string()))?;
    Ok(line.split_whitespace().collect())
}

fn parse_f64(tok: &str) -> Result<f64, CubeError> {
    tok.parse()
        .map_err(|_| CubeError::ParseError(format!("invalid number: {tok:?}")))
}

fn parse_i64(tok: &str) -> Result<i64, CubeError> {
    tok.parse()
        .map_err(|_| CubeError::ParseError(format!("invalid integer: {tok:?}")))
}

/// Parses a whole CUBE file into a [`VolumeGrid`], with
/// [`VolumeGrid::atoms`] always set (even to an empty [`Molecule`] for
/// `natoms == 0`) — CUBE always states this block, unlike CCP4/DSN6/DX,
/// which have no such concept at all.
pub(crate) fn parse_cube(text: &str) -> Result<VolumeGrid, CubeError> {
    let mut lines = text.lines();
    lines
        .next()
        .ok_or_else(|| CubeError::ParseError("empty file".to_string()))?;
    lines
        .next()
        .ok_or_else(|| CubeError::ParseError("empty file".to_string()))?;

    let header = next_tokens(&mut lines)?;
    if header.len() < 4 {
        return Err(CubeError::ParseError(
            "malformed atom-count/origin line".to_string(),
        ));
    }
    let raw_natoms = parse_i64(header[0])?;
    if raw_natoms < 0 {
        return Err(CubeError::MultipleValuesUnsupported(
            "a negative atom count indicates orbital data".to_string(),
        ));
    }
    if header.len() == 5 {
        let num_val = parse_i64(header[4])?;
        if num_val != 1 {
            return Err(CubeError::MultipleValuesUnsupported(format!(
                "the header states {num_val} values per grid point"
            )));
        }
    }
    let natoms = raw_natoms as usize;
    let origin_raw = [
        parse_f64(header[1])?,
        parse_f64(header[2])?,
        parse_f64(header[3])?,
    ];

    let mut counts = [0i64; 3];
    let mut axis_raw = [[0.0f64; 3]; 3];
    for (i, slot) in axis_raw.iter_mut().enumerate() {
        let row = next_tokens(&mut lines)?;
        if row.len() < 4 {
            return Err(CubeError::ParseError("malformed axis row".to_string()));
        }
        counts[i] = parse_i64(row[0])?;
        *slot = [parse_f64(row[1])?, parse_f64(row[2])?, parse_f64(row[3])?];
    }

    let is_angstrom = counts.iter().any(|&c| c < 0);
    let scale = if is_angstrom { 1.0 } else { BOHR_TO_ANGSTROM };
    let dims_file = [
        counts[0].unsigned_abs() as usize,
        counts[1].unsigned_abs() as usize,
        counts[2].unsigned_abs() as usize,
    ];
    let origin = Point3::new(origin_raw[0], origin_raw[1], origin_raw[2]) * scale;
    let axis_vectors: [Point3; 3] = std::array::from_fn(|i| {
        Point3::new(axis_raw[i][0], axis_raw[i][1], axis_raw[i][2]) * scale
    });

    let mut molecule = Molecule::new();
    let mut coords = Vec::with_capacity(natoms);
    for _ in 0..natoms {
        let row = next_tokens(&mut lines)?;
        if row.len() < 5 {
            return Err(CubeError::ParseError("malformed atom line".to_string()));
        }
        let atomic_number = parse_i64(row[0])?;
        let element = u8::try_from(atomic_number)
            .ok()
            .and_then(Element::new)
            .unwrap_or(Element::UNKNOWN);
        molecule.add_atom(Atom::new(element));
        let x = parse_f64(row[2])? * scale;
        let y = parse_f64(row[3])? * scale;
        let z = parse_f64(row[4])? * scale;
        coords.push(Point3::new(x, y, z));
    }
    molecule.set_coords3(coords)?;

    let remaining: String = lines.collect::<Vec<_>>().join(" ");
    let values: Vec<f64> = remaining
        .split_whitespace()
        .map(parse_f64)
        .collect::<Result<_, _>>()?;

    let expected = dims_file[0] * dims_file[1] * dims_file[2];
    if values.len() != expected {
        return Err(CubeError::ParseError(format!(
            "expected {expected} values ({dims_file:?}), got {}",
            values.len()
        )));
    }

    let mut grid = VolumeGrid::from_source_order(
        [dims_file[2], dims_file[1], dims_file[0]],
        [Axis::Z, Axis::Y, Axis::X],
        origin,
        [axis_vectors[2], axis_vectors[1], axis_vectors[0]],
        values,
        None,
    )?;
    grid.set_atoms(molecule);
    Ok(grid)
}

/// Writes a [`VolumeGrid`] as a CUBE file — always Bohr, always positive
/// axis counts (see the module doc).
pub(crate) fn write_cube(grid: &VolumeGrid) -> String {
    let dims = grid.dims();
    let origin = grid.origin() / BOHR_TO_ANGSTROM;
    let axes = grid.axes();
    let atoms = grid.atoms();
    let natoms = atoms.map(Molecule::num_atoms).unwrap_or(0);

    let mut out = String::new();
    out.push_str("Written by chem\n");
    out.push_str("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
    out.push_str(&format!(
        "{natoms:5}{:12.6}{:12.6}{:12.6}\n",
        origin.x, origin.y, origin.z
    ));
    for i in 0..3 {
        let a = axes[i] / BOHR_TO_ANGSTROM;
        out.push_str(&format!(
            "{:5}{:12.6}{:12.6}{:12.6}\n",
            dims[i], a.x, a.y, a.z
        ));
    }

    if let Some(mol) = atoms {
        let coords = mol.coords3().unwrap_or(&[]);
        for (i, atom) in mol.atoms().iter().enumerate() {
            let p = coords.get(i).copied().unwrap_or(Point3::ORIGIN) / BOHR_TO_ANGSTROM;
            out.push_str(&format!(
                "{:5}{:12.6}{:12.6}{:12.6}{:12.6}\n",
                atom.element().atomic_number,
                0.0,
                p.x,
                p.y,
                p.z
            ));
        }
    }

    let mut col = 0;
    for ix in 0..dims[0] {
        for iy in 0..dims[1] {
            for iz in 0..dims[2] {
                out.push_str(&format!("{:13.5e}", grid.value(ix, iy, iz)));
                col += 1;
                if col % 6 == 0 {
                    out.push('\n');
                } else {
                    out.push(' ');
                }
            }
        }
    }
    if col % 6 != 0 {
        out.push('\n');
    }

    out
}

/// [`crate::io::format::ReadFn`] for CUBE.
pub(crate) fn read_cube_with_options(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_cube(text) {
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

/// [`crate::io::format::ByteWriteVolumeFn`] for CUBE.
pub(crate) fn write_cube_bytes(grid: &VolumeGrid, _options: &WriteOptions) -> Vec<u8> {
    write_cube(grid).into_bytes()
}

/// Buffers the whole input and parses it once — a CUBE file is always
/// exactly one grid, the same "one record" shape
/// [`crate::io::dcd::DcdSupplier`]/[`crate::io::nctraj::NctrajSupplier`]
/// already have for their own always-one-record files.
pub struct CubeSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl CubeSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match parse_cube(&text) {
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

impl Iterator for CubeSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_grid() -> VolumeGrid {
        let mut grid = VolumeGrid::new(
            [1, 1, 2],
            Point3::new(0.0, 0.0, 0.0),
            [
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 0.5),
            ],
            vec![10.0, 20.0],
            None,
        )
        .expect("valid grid");

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::new(6).unwrap())); // C
        mol.add_atom(Atom::new(Element::new(8).unwrap())); // O
        mol.set_coords3(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.5, 0.0, 0.0)])
            .expect("valid coords");
        grid.set_atoms(mol);
        grid
    }

    #[test]
    fn test_an_orthogonal_round_trip_with_atoms() {
        let grid = sample_grid();
        let text = write_cube(&grid);
        let back = parse_cube(&text).expect("valid CUBE");

        assert_eq!(back.dims(), grid.dims());
        for iz in 0..2 {
            assert!(
                (back.value(0, 0, iz) - grid.value(0, 0, iz)).abs() < 1e-3,
                "{} vs {}",
                back.value(0, 0, iz),
                grid.value(0, 0, iz)
            );
        }

        let atoms = back.atoms().expect("atoms survive");
        assert_eq!(atoms.num_atoms(), 2);
        assert_eq!(atoms.atoms()[0].element(), Element::new(6).unwrap());
        assert_eq!(atoms.atoms()[1].element(), Element::new(8).unwrap());
        let coords = atoms.coords3().expect("coords survive");
        assert!((coords[0].x - 0.0).abs() < 1e-3);
        assert!((coords[1].x - 1.5).abs() < 1e-3, "{}", coords[1].x);
    }

    #[test]
    fn test_a_negative_axis_count_means_angstrom_not_bohr() {
        // Hand-computed: a negative axis count means the whole header
        // (origin, axis vectors, atom positions) is already Angstrom, no
        // Bohr conversion -- confirmed against `pdb2pqr`'s real
        // `write_cube`, which deliberately emits this exact convention.
        let text = "\
comment one
comment two
1   1.0 2.0 3.0
-2  0.5 0.0 0.0
3   0.0 0.5 0.0
4   0.0 0.0 0.5
6  0.0  1.0 2.0 3.0
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
";
        let grid = parse_cube(text).expect("valid CUBE");
        assert!((grid.origin().x - 1.0).abs() < 1e-6, "{}", grid.origin().x);
        assert!((grid.origin().y - 2.0).abs() < 1e-6);
        assert!((grid.origin().z - 3.0).abs() < 1e-6);

        let axes = grid.axes();
        assert!((axes[0].x - 0.5).abs() < 1e-6, "{}", axes[0].x);
        assert!((axes[1].y - 0.5).abs() < 1e-6);
        assert!((axes[2].z - 0.5).abs() < 1e-6);

        let atoms = grid.atoms().expect("atoms present");
        let coords = atoms.coords3().expect("coords present");
        assert!((coords[0].x - 1.0).abs() < 1e-6);
        assert!((coords[0].y - 2.0).abs() < 1e-6);
        assert!((coords[0].z - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_values_are_read_in_the_correct_axis_order() {
        // File order is X outer, Y middle, Z inner; nx=2, ny=1, nz=3
        // exercises a genuine permutation (unlike a grid where two of the
        // three dims are 1, which would round-trip correctly even with the
        // permutation logic broken).
        let text = "\
comment one
comment two
0   0.0 0.0 0.0
2   1.0 0.0 0.0
1   0.0 1.0 0.0
3   0.0 0.0 1.0
1.0 2.0 3.0 4.0 5.0 6.0
";
        let grid = parse_cube(text).expect("valid CUBE");
        assert_eq!(grid.dims(), [2, 1, 3]);
        assert_eq!(grid.value(0, 0, 0), 1.0);
        assert_eq!(grid.value(1, 0, 0), 4.0);
        assert_eq!(grid.value(0, 0, 1), 2.0);
        assert_eq!(grid.value(1, 0, 1), 5.0);
        assert_eq!(grid.value(0, 0, 2), 3.0);
        assert_eq!(grid.value(1, 0, 2), 6.0);
    }

    #[test]
    fn test_a_zero_atom_grid_still_gets_an_empty_molecule() {
        let text = "\
comment one
comment two
0   0.0 0.0 0.0
1   1.0 0.0 0.0
1   0.0 1.0 0.0
1   0.0 0.0 1.0
0.0
";
        let grid = parse_cube(text).expect("valid CUBE");
        let atoms = grid.atoms().expect("CUBE always states an atom block");
        assert_eq!(atoms.num_atoms(), 0);
    }

    #[test]
    fn test_a_negative_atom_count_is_refused() {
        let text = "\
comment one
comment two
-1   0.0 0.0 0.0
1   1.0 0.0 0.0
1   0.0 1.0 0.0
1   0.0 0.0 1.0
";
        let err = parse_cube(text).unwrap_err();
        assert!(
            matches!(err, CubeError::MultipleValuesUnsupported(_)),
            "{err}"
        );
    }

    #[test]
    fn test_a_stated_num_val_other_than_one_is_refused() {
        let text = "\
comment one
comment two
0   0.0 0.0 0.0 2
1   1.0 0.0 0.0
1   0.0 1.0 0.0
1   0.0 0.0 1.0
";
        let err = parse_cube(text).unwrap_err();
        assert!(
            matches!(err, CubeError::MultipleValuesUnsupported(_)),
            "{err}"
        );
    }
}
