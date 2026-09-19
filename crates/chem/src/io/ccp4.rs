//! CCP4/MRC (#332) -- the standard electron-density and cryo-EM map
//! format, and the one this crate takes most seriously: every real
//! structural-biology tool reads it.
//!
//! Every load-bearing detail here was confirmed empirically against
//! `gemmi` (the issue's own named oracle), not reconstructed from memory:
//! a synthetic non-cubic-cell map was built with `gemmi.Ccp4Map`, every
//! header word read back via its own `header_i32`/`header_float`, and the
//! raw post-header bytes inspected directly.
//!
//! **Layout** (1-based "word" numbering the spec itself uses; byte offset
//! = `(word-1)*4`): 1-3 `NX,NY,NZ`; 4 `MODE`; 5-7 `NXSTART,NYSTART,NZSTART`;
//! 8-10 `MX,MY,MZ` (the *whole* unit cell's own grid sampling, used only to
//! derive spacing -- not necessarily equal to `NX,NY,NZ`); 11-13 cell
//! `a,b,c`; 14-16 cell `alpha,beta,gamma`; 17-19 `MAPC,MAPR,MAPS`; 20-22
//! `AMIN,AMAX,AMEAN` (never trusted, see below); 23 `ISPG`; 24 `NSYMBT`
//! (bytes of symmetry records between the 1024-byte header and the density
//! data -- confirmed directly: a file with `NSYMBT=80` had its density
//! start exactly at byte `1024+80`); 25-49 skew/extra (never read or
//! written); 50-52 `ORIGIN`; 53 `"MAP "` (byte 208); 54 the machine stamp;
//! 55 `RMS`; 56 `NLABL`; 57-256 ten 80-char labels (parsed and discarded).
//!
//! **On-disk value order**: confirmed directly (four distinct values set
//! at known grid coordinates, raw post-header bytes inspected) -- columns
//! (`NX`) vary fastest, exactly `VolumeGrid`'s own canonical X-fastest
//! storage. `MAPC`/`MAPR`/`MAPS` (1/2/3 = crystallographic a/b/c) only
//! permute *which* canonical axis each file dimension is --
//! [`VolumeGrid::from_source_order`]'s exact contract, fed unreversed this
//! time (unlike CUBE, whose file order was the other way around).
//!
//! **Endianness**: the machine-stamp word (byte 212) is a direct byte
//! pattern, not a numeral needing a decode-both-ways guess the way DCD's
//! leading record length does -- little-endian IEEE is `44 41 00 00`
//! (confirmed directly), big-endian IEEE's well-documented complementary
//! pattern is `11 11 00 00`. Neither matching is refused, not guessed.
//!
//! **`RMS` is standard deviation around the mean, not `sqrt(mean(x^2))`.**
//! A genuine, non-obvious finding: `VolumeGrid::statistics`'s own `rms`
//! field is RMS around zero, a real and differently-defined quantity
//! confirmed by cross-checking `gemmi`'s own written `RMS` header word
//! against a fixture with a hand-computable mean. This module computes
//! CCP4's own real formula locally on write rather than reusing that field.
//!
//! **`ORIGIN` vs `NXSTART`, already decided -- before this story existed.**
//! `VolumeGrid::origin`'s own doc comment (#312) already states the
//! resolution: prefer `ORIGIN` when the format states one, falling back to
//! `NXSTART * step` only when it does not. "States one" here means any of
//! `ORIGIN`'s three components is nonzero -- the standard real-tool
//! heuristic, since an all-zero `ORIGIN` is indistinguishable from a
//! pre-`ORIGIN`-era file that never populated it. Applied identically on
//! write: this crate always writes a real `ORIGIN` and zero
//! `NXSTART`/`NYSTART`/`NZSTART`.
//!
//! **Modes** 0 (signed byte), 1 (16-bit int), 2 (32-bit float) and 6
//! (unsigned 16-bit) are all genuinely read, not just mode 2 -- the issue's
//! own text is explicit that the others are "where a reader quietly
//! produces noise." 3/4 (complex) are refused by name. Writing always uses
//! mode 2, the same least-ambiguous-convention posture every prior format's
//! writer in this crate already took.
//!
//! **Writing never permutes** -- always `MAPC=1,MAPR=2,MAPS=3`, so
//! `NX,NY,NZ` is `VolumeGrid::dims()` directly and the value stream needs
//! no reordering at all, unlike CUBE's reversed write loop. A grid with no
//! stated `VolumeGrid::cell` gets a synthesized orthogonal "cell equals
//! box" one (`a = |axis_x| * dims[0]`, etc., 90 degree angles) -- the same
//! disclosed no-cell fallback `io::lammps`'s own data-file writer already
//! established.

use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;
use crate::core::volume::{Axis, VolumeGrid};
use crate::io::errors::{Ccp4Error, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

const HEADER_LEN: usize = 1024;
const LE_MACHINE_STAMP: [u8; 4] = [0x44, 0x41, 0x00, 0x00];
const BE_MACHINE_STAMP: [u8; 4] = [0x11, 0x11, 0x00, 0x00];

fn word_offset(word: usize) -> usize {
    (word - 1) * 4
}

fn read_i32(bytes: &[u8], word: usize, big_endian: bool) -> Result<i32, Ccp4Error> {
    let off = word_offset(word);
    let raw: [u8; 4] = bytes
        .get(off..off + 4)
        .ok_or_else(|| Ccp4Error::ParseError("header truncated".to_string()))?
        .try_into()
        .unwrap();
    Ok(if big_endian {
        i32::from_be_bytes(raw)
    } else {
        i32::from_le_bytes(raw)
    })
}

fn read_f32(bytes: &[u8], word: usize, big_endian: bool) -> Result<f32, Ccp4Error> {
    let off = word_offset(word);
    let raw: [u8; 4] = bytes
        .get(off..off + 4)
        .ok_or_else(|| Ccp4Error::ParseError("header truncated".to_string()))?
        .try_into()
        .unwrap();
    Ok(if big_endian {
        f32::from_be_bytes(raw)
    } else {
        f32::from_le_bytes(raw)
    })
}

fn axis_from_map_code(code: i32) -> Result<Axis, Ccp4Error> {
    match code {
        1 => Ok(Axis::X),
        2 => Ok(Axis::Y),
        3 => Ok(Axis::Z),
        other => Err(Ccp4Error::ParseError(format!(
            "invalid MAPC/MAPR/MAPS axis code: {other}"
        ))),
    }
}

/// Parses a whole CCP4/MRC file into a [`VolumeGrid`].
pub(crate) fn parse_ccp4(bytes: &[u8]) -> Result<VolumeGrid, Ccp4Error> {
    if bytes.len() < HEADER_LEN {
        return Err(Ccp4Error::ParseError(
            "file shorter than the 1024-byte header".to_string(),
        ));
    }
    if &bytes[208..212] != b"MAP " {
        return Err(Ccp4Error::InvalidMagicNumber);
    }
    let stamp: [u8; 4] = bytes[212..216].try_into().unwrap();
    let big_endian = if stamp == LE_MACHINE_STAMP {
        false
    } else if stamp == BE_MACHINE_STAMP {
        true
    } else {
        return Err(Ccp4Error::UnknownMachineStamp);
    };

    let nx = read_i32(bytes, 1, big_endian)? as usize;
    let ny = read_i32(bytes, 2, big_endian)? as usize;
    let nz = read_i32(bytes, 3, big_endian)? as usize;
    let mode = read_i32(bytes, 4, big_endian)?;
    let nxstart = read_i32(bytes, 5, big_endian)? as f64;
    let nystart = read_i32(bytes, 6, big_endian)? as f64;
    let nzstart = read_i32(bytes, 7, big_endian)? as f64;
    let mx = read_i32(bytes, 8, big_endian)?;
    let my = read_i32(bytes, 9, big_endian)?;
    let mz = read_i32(bytes, 10, big_endian)?;
    if mx <= 0 || my <= 0 || mz <= 0 {
        return Err(Ccp4Error::ParseError(
            "MX/MY/MZ must be positive".to_string(),
        ));
    }

    let cell = UnitCell::new(
        read_f32(bytes, 11, big_endian)? as f64,
        read_f32(bytes, 12, big_endian)? as f64,
        read_f32(bytes, 13, big_endian)? as f64,
        read_f32(bytes, 14, big_endian)? as f64,
        read_f32(bytes, 15, big_endian)? as f64,
        read_f32(bytes, 16, big_endian)? as f64,
    );

    let mapc = axis_from_map_code(read_i32(bytes, 17, big_endian)?)?;
    let mapr = axis_from_map_code(read_i32(bytes, 18, big_endian)?)?;
    let maps = axis_from_map_code(read_i32(bytes, 19, big_endian)?)?;

    let nsymbt = read_i32(bytes, 24, big_endian)?.max(0) as usize;

    let origin_raw = [
        read_f32(bytes, 50, big_endian)? as f64,
        read_f32(bytes, 51, big_endian)? as f64,
        read_f32(bytes, 52, big_endian)? as f64,
    ];

    // MX/MY/MZ and the cell are always in terms of the crystallographic
    // a/b/c axes, independent of which file dimension (column/row/section)
    // MAPC/MAPR/MAPS says represents which -- so the *column* axis's own
    // step vector is whichever of step_a/step_b/step_c its own MAPC value
    // names, not positionally the first one. Getting this wrong silently
    // produces a plausible-looking, transposed grid on any permuted file,
    // exactly the "expensive mistake" this story exists to avoid -- caught
    // during this story's own fixture-verification, not assumed.
    let basis = cell.basis();
    let step_for = |axis: Axis| -> Point3 {
        match axis {
            Axis::X => basis[0] / mx as f64,
            Axis::Y => basis[1] / my as f64,
            Axis::Z => basis[2] / mz as f64,
        }
    };
    let column_step = step_for(mapc);
    let row_step = step_for(mapr);
    let section_step = step_for(maps);

    let origin = if origin_raw.iter().any(|&v| v != 0.0) {
        Point3::new(origin_raw[0], origin_raw[1], origin_raw[2])
    } else {
        column_step * nxstart + row_step * nystart + section_step * nzstart
    };

    let data_start = HEADER_LEN + nsymbt;
    let expected_count = nx * ny * nz;

    let values: Vec<f64> = match mode {
        0 => {
            let raw = bytes
                .get(data_start..data_start + expected_count)
                .ok_or_else(|| Ccp4Error::ParseError("density data truncated".to_string()))?;
            raw.iter().map(|&b| b as i8 as f64).collect()
        }
        1 => {
            let byte_len = expected_count * 2;
            let raw = bytes
                .get(data_start..data_start + byte_len)
                .ok_or_else(|| Ccp4Error::ParseError("density data truncated".to_string()))?;
            (0..expected_count)
                .map(|i| {
                    let b: [u8; 2] = raw[i * 2..i * 2 + 2].try_into().unwrap();
                    let v = if big_endian {
                        i16::from_be_bytes(b)
                    } else {
                        i16::from_le_bytes(b)
                    };
                    v as f64
                })
                .collect()
        }
        2 => {
            let byte_len = expected_count * 4;
            let raw = bytes
                .get(data_start..data_start + byte_len)
                .ok_or_else(|| Ccp4Error::ParseError("density data truncated".to_string()))?;
            (0..expected_count)
                .map(|i| {
                    let b: [u8; 4] = raw[i * 4..i * 4 + 4].try_into().unwrap();
                    let v = if big_endian {
                        f32::from_be_bytes(b)
                    } else {
                        f32::from_le_bytes(b)
                    };
                    v as f64
                })
                .collect()
        }
        6 => {
            let byte_len = expected_count * 2;
            let raw = bytes
                .get(data_start..data_start + byte_len)
                .ok_or_else(|| Ccp4Error::ParseError("density data truncated".to_string()))?;
            (0..expected_count)
                .map(|i| {
                    let b: [u8; 2] = raw[i * 2..i * 2 + 2].try_into().unwrap();
                    let v = if big_endian {
                        u16::from_be_bytes(b)
                    } else {
                        u16::from_le_bytes(b)
                    };
                    v as f64
                })
                .collect()
        }
        other => return Err(Ccp4Error::UnsupportedMode(other)),
    };

    let grid = VolumeGrid::from_source_order(
        [nx, ny, nz],
        [mapc, mapr, maps],
        origin,
        [column_step, row_step, section_step],
        values,
        Some(cell),
    )?;
    Ok(grid)
}

/// Standard deviation around the mean -- CCP4's own real `RMS` definition,
/// confirmed against `gemmi`'s own output (see the module doc). Different
/// from [`VolumeGrid::statistics`]'s own, differently-defined `rms` field.
fn ccp4_rms(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / n;
    let mean_sq = values.iter().map(|v| v * v).sum::<f64>() / n;
    (mean_sq - mean * mean).max(0.0).sqrt()
}

fn put_i32(out: &mut [u8], word: usize, v: i32) {
    let off = word_offset(word);
    out[off..off + 4].copy_from_slice(&v.to_le_bytes());
}

fn put_f32(out: &mut [u8], word: usize, v: f32) {
    let off = word_offset(word);
    out[off..off + 4].copy_from_slice(&v.to_le_bytes());
}

/// Writes a [`VolumeGrid`] as a CCP4/MRC file -- always mode 2, always
/// unpermuted, little-endian (see the module doc).
pub(crate) fn write_ccp4(grid: &VolumeGrid) -> Vec<u8> {
    let dims = grid.dims();
    let axes = grid.axes();
    let origin = grid.origin();
    let values = grid.values();

    let cell = grid.cell().unwrap_or_else(|| {
        UnitCell::new(
            axes[0].length() * dims[0] as f64,
            axes[1].length() * dims[1] as f64,
            axes[2].length() * dims[2] as f64,
            90.0,
            90.0,
            90.0,
        )
    });

    let basis = cell.basis();
    let sampling = |basis_len: f64, axis_len: f64| -> i32 {
        if axis_len <= 0.0 {
            1
        } else {
            (basis_len / axis_len).round().max(1.0) as i32
        }
    };
    let mx = sampling(basis[0].length(), axes[0].length());
    let my = sampling(basis[1].length(), axes[1].length());
    let mz = sampling(basis[2].length(), axes[2].length());

    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    };
    let rms = ccp4_rms(values);

    let mut out = vec![0u8; HEADER_LEN];
    put_i32(&mut out, 1, dims[0] as i32);
    put_i32(&mut out, 2, dims[1] as i32);
    put_i32(&mut out, 3, dims[2] as i32);
    put_i32(&mut out, 4, 2); // MODE: float32
    put_i32(&mut out, 5, 0); // NXSTART
    put_i32(&mut out, 6, 0); // NYSTART
    put_i32(&mut out, 7, 0); // NZSTART
    put_i32(&mut out, 8, mx);
    put_i32(&mut out, 9, my);
    put_i32(&mut out, 10, mz);
    put_f32(&mut out, 11, cell.a as f32);
    put_f32(&mut out, 12, cell.b as f32);
    put_f32(&mut out, 13, cell.c as f32);
    put_f32(&mut out, 14, cell.alpha as f32);
    put_f32(&mut out, 15, cell.beta as f32);
    put_f32(&mut out, 16, cell.gamma as f32);
    put_i32(&mut out, 17, 1); // MAPC: X
    put_i32(&mut out, 18, 2); // MAPR: Y
    put_i32(&mut out, 19, 3); // MAPS: Z
    put_f32(&mut out, 20, min as f32);
    put_f32(&mut out, 21, max as f32);
    put_f32(&mut out, 22, mean as f32);
    put_i32(&mut out, 23, 1); // ISPG: P1
    put_i32(&mut out, 24, 0); // NSYMBT: no symmetry records written
    put_f32(&mut out, 50, origin.x as f32);
    put_f32(&mut out, 51, origin.y as f32);
    put_f32(&mut out, 52, origin.z as f32);
    out[208..212].copy_from_slice(b"MAP ");
    out[212..216].copy_from_slice(&LE_MACHINE_STAMP);
    put_f32(&mut out, 55, rms as f32);
    put_i32(&mut out, 56, 1); // NLABL

    let label = b"Written by chem";
    let label_start = word_offset(57);
    out[label_start..label_start + label.len()].copy_from_slice(label);

    for &v in values {
        out.extend_from_slice(&(v as f32).to_le_bytes());
    }

    out
}

/// [`crate::io::format::ByteReadFn`] for CCP4/MRC.
pub(crate) fn read_ccp4_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_ccp4(bytes) {
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

/// [`crate::io::format::ByteWriteVolumeFn`] for CCP4/MRC.
pub(crate) fn write_ccp4_bytes(grid: &VolumeGrid, _options: &WriteOptions) -> Vec<u8> {
    write_ccp4(grid)
}

/// Buffers the whole input and parses it once -- a CCP4/MRC file is always
/// exactly one grid, the same "one record" shape
/// [`crate::io::cube::CubeSupplier`] already has.
pub struct Ccp4Supplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl Ccp4Supplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match std::io::Read::read_to_end(&mut reader, &mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match parse_ccp4(&bytes) {
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

impl Iterator for Ccp4Supplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a full 1024-byte CCP4/MRC header, byte-exact, independent of
    /// this module's own writer -- so a fixture-based test exercises the
    /// *reader* against ground truth, not the reader against the writer's
    /// own (possibly equally wrong) assumptions.
    #[allow(clippy::too_many_arguments)]
    fn build_header(
        nx: i32,
        ny: i32,
        nz: i32,
        mode: i32,
        nxstart: i32,
        nystart: i32,
        nzstart: i32,
        mx: i32,
        my: i32,
        mz: i32,
        cell: UnitCell,
        mapc: i32,
        mapr: i32,
        maps: i32,
        origin: [f32; 3],
        big_endian: bool,
    ) -> Vec<u8> {
        let mut out = vec![0u8; HEADER_LEN];
        let put_i32 = |out: &mut [u8], word: usize, v: i32| {
            let off = word_offset(word);
            let bytes = if big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            };
            out[off..off + 4].copy_from_slice(&bytes);
        };
        let put_f32 = |out: &mut [u8], word: usize, v: f32| {
            let off = word_offset(word);
            let bytes = if big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            };
            out[off..off + 4].copy_from_slice(&bytes);
        };

        put_i32(&mut out, 1, nx);
        put_i32(&mut out, 2, ny);
        put_i32(&mut out, 3, nz);
        put_i32(&mut out, 4, mode);
        put_i32(&mut out, 5, nxstart);
        put_i32(&mut out, 6, nystart);
        put_i32(&mut out, 7, nzstart);
        put_i32(&mut out, 8, mx);
        put_i32(&mut out, 9, my);
        put_i32(&mut out, 10, mz);
        put_f32(&mut out, 11, cell.a as f32);
        put_f32(&mut out, 12, cell.b as f32);
        put_f32(&mut out, 13, cell.c as f32);
        put_f32(&mut out, 14, cell.alpha as f32);
        put_f32(&mut out, 15, cell.beta as f32);
        put_f32(&mut out, 16, cell.gamma as f32);
        put_i32(&mut out, 17, mapc);
        put_i32(&mut out, 18, mapr);
        put_i32(&mut out, 19, maps);
        put_i32(&mut out, 24, 0); // NSYMBT
        put_f32(&mut out, 50, origin[0]);
        put_f32(&mut out, 51, origin[1]);
        put_f32(&mut out, 52, origin[2]);
        out[208..212].copy_from_slice(b"MAP ");
        out[212..216].copy_from_slice(if big_endian {
            &BE_MACHINE_STAMP
        } else {
            &LE_MACHINE_STAMP
        });
        out
    }

    fn append_f32_values(out: &mut Vec<u8>, values: &[f32], big_endian: bool) {
        for &v in values {
            out.extend_from_slice(&if big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            });
        }
    }

    fn orthogonal_header(
        nx: i32,
        ny: i32,
        nz: i32,
        mode: i32,
        mapc: i32,
        mapr: i32,
        maps: i32,
    ) -> Vec<u8> {
        build_header(
            nx,
            ny,
            nz,
            mode,
            0,
            0,
            0,
            nx,
            ny,
            nz,
            UnitCell::new(10.0, 20.0, 30.0, 90.0, 90.0, 90.0),
            mapc,
            mapr,
            maps,
            [0.0, 0.0, 0.0],
            false,
        )
    }

    #[test]
    fn test_an_orthogonal_round_trip_with_a_non_cubic_cell() {
        let grid = VolumeGrid::new(
            [2, 3, 4],
            Point3::new(1.0, 2.0, 3.0),
            [
                Point3::new(5.0, 0.0, 0.0),
                Point3::new(0.0, 20.0 / 3.0, 0.0),
                Point3::new(0.0, 0.0, 30.0 / 4.0),
            ],
            (0..24).map(|i| i as f64).collect(),
            Some(UnitCell::new(10.0, 20.0, 30.0, 90.0, 90.0, 90.0)),
        )
        .expect("valid grid");

        let bytes = write_ccp4(&grid);
        let back = parse_ccp4(&bytes).expect("valid CCP4");

        assert_eq!(back.dims(), grid.dims());
        let back_cell = back.cell().expect("cell survives");
        assert!((back_cell.a - 10.0).abs() < 1e-3);
        assert!((back_cell.b - 20.0).abs() < 1e-3);
        assert!((back_cell.c - 30.0).abs() < 1e-3);
        for ix in 0..2 {
            for iy in 0..3 {
                for iz in 0..4 {
                    assert!(
                        (back.value(ix, iy, iz) - grid.value(ix, iy, iz)).abs() < 1e-2,
                        "{} vs {} at {ix},{iy},{iz}",
                        back.value(ix, iy, iz),
                        grid.value(ix, iy, iz)
                    );
                }
            }
        }
        let back_origin = back.origin();
        assert!((back_origin.x - 1.0).abs() < 1e-2, "{}", back_origin.x);
        assert!((back_origin.y - 2.0).abs() < 1e-2, "{}", back_origin.y);
        assert!((back_origin.z - 3.0).abs() < 1e-2, "{}", back_origin.z);
    }

    #[test]
    fn test_a_permuted_axis_non_cubic_cell_matches_independently_computed_values() {
        // Hand-derived, independent of this module's own code: MAPC=3
        // (columns represent crystallographic c), MAPR=1 (rows -> a),
        // MAPS=2 (sections -> b); NX=2,NY=1,NZ=3; cell (10,20,30,
        // orthogonal). File values, columns-fastest: 1,2,3,4,5,6 at
        // (col,sec) = (0,0),(1,0),(0,1),(1,1),(0,2),(1,2).
        //
        // Working through `VolumeGrid::from_source_order`'s own documented
        // index algebra by hand gives canonical dims [1,3,2] (X from rows,
        // Y from sections, Z from columns) and
        // grid.value(0,iy,iz) = [[1,3,5],[2,4,6]][iz][iy] -- i.e.
        // value(0,0,0)=1, value(0,1,0)=3, value(0,2,0)=5,
        // value(0,0,1)=2, value(0,1,1)=4, value(0,2,1)=6.
        let mut bytes = orthogonal_header(2, 1, 3, 2, 3, 1, 2);
        append_f32_values(&mut bytes, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], false);

        let grid = parse_ccp4(&bytes).expect("valid CCP4");
        assert_eq!(grid.dims(), [1, 3, 2]);
        assert_eq!(grid.value(0, 0, 0), 1.0);
        assert_eq!(grid.value(0, 1, 0), 3.0);
        assert_eq!(grid.value(0, 2, 0), 5.0);
        assert_eq!(grid.value(0, 0, 1), 2.0);
        assert_eq!(grid.value(0, 1, 1), 4.0);
        assert_eq!(grid.value(0, 2, 1), 6.0);

        // Columns (MAPC=3=c) have length 2 -> their step vector is
        // basis_c/MZ = (0,0,30)/3 = (0,0,10), landing on canonical Z.
        let axes = grid.axes();
        assert!((axes[2].z - 10.0).abs() < 1e-6, "{}", axes[2].z);
        // Rows (MAPR=1=a) have length 1 -> basis_a/MX = (10,0,0)/2 =
        // (5,0,0), landing on canonical X.
        assert!((axes[0].x - 5.0).abs() < 1e-6, "{}", axes[0].x);
        // Sections (MAPS=2=b) have length 3 -> basis_b/MY = (0,20,0)/1 =
        // (0,20,0), landing on canonical Y.
        assert!((axes[1].y - 20.0).abs() < 1e-6, "{}", axes[1].y);
    }

    #[test]
    fn test_a_big_endian_file_reads_correctly() {
        let mut bytes = build_header(
            2,
            1,
            1,
            2,
            0,
            0,
            0,
            2,
            1,
            1,
            UnitCell::new(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
            1,
            2,
            3,
            [0.0, 0.0, 0.0],
            true,
        );
        append_f32_values(&mut bytes, &[42.0, 43.0], true);

        let grid = parse_ccp4(&bytes).expect("valid CCP4");
        assert_eq!(grid.value(0, 0, 0), 42.0);
        assert_eq!(grid.value(1, 0, 0), 43.0);
    }

    #[test]
    fn test_origin_wins_over_nxstart_when_stated() {
        let mut bytes = orthogonal_header(1, 1, 1, 2, 1, 2, 3);
        // ORIGIN nonzero -- must win over NXSTART, which stays 0 here.
        bytes[word_offset(50)..word_offset(50) + 4].copy_from_slice(&99.0f32.to_le_bytes());
        append_f32_values(&mut bytes, &[7.0], false);

        let grid = parse_ccp4(&bytes).expect("valid CCP4");
        assert!((grid.origin().x - 99.0).abs() < 1e-3, "{}", grid.origin().x);
    }

    #[test]
    fn test_nxstart_is_used_only_when_origin_is_all_zero() {
        let mut bytes = orthogonal_header(1, 1, 1, 2, 1, 2, 3);
        // NXSTART=2, ORIGIN left at zero -- must fall back to NXSTART*step.
        bytes[word_offset(5)..word_offset(5) + 4].copy_from_slice(&2i32.to_le_bytes());
        append_f32_values(&mut bytes, &[7.0], false);

        let grid = parse_ccp4(&bytes).expect("valid CCP4");
        // step_a = 10.0/nx(=1) = 10.0; origin.x = 2 * 10.0 = 20.0.
        assert!((grid.origin().x - 20.0).abs() < 1e-3, "{}", grid.origin().x);
    }

    #[test]
    fn test_every_supported_mode_decodes_correctly() {
        // Mode 0: signed byte.
        let mut bytes = orthogonal_header(2, 1, 1, 0, 1, 2, 3);
        bytes.extend_from_slice(&[100u8, (-50i8) as u8]);
        let grid = parse_ccp4(&bytes).expect("valid mode 0");
        assert_eq!(grid.value(0, 0, 0), 100.0);
        assert_eq!(grid.value(1, 0, 0), -50.0);

        // Mode 1: signed 16-bit.
        let mut bytes = orthogonal_header(2, 1, 1, 1, 1, 2, 3);
        bytes.extend_from_slice(&1000i16.to_le_bytes());
        bytes.extend_from_slice(&(-2000i16).to_le_bytes());
        let grid = parse_ccp4(&bytes).expect("valid mode 1");
        assert_eq!(grid.value(0, 0, 0), 1000.0);
        assert_eq!(grid.value(1, 0, 0), -2000.0);

        // Mode 6: unsigned 16-bit.
        let mut bytes = orthogonal_header(2, 1, 1, 6, 1, 2, 3);
        bytes.extend_from_slice(&40000u16.to_le_bytes());
        bytes.extend_from_slice(&0u16.to_le_bytes());
        let grid = parse_ccp4(&bytes).expect("valid mode 6");
        assert_eq!(grid.value(0, 0, 0), 40000.0);
        assert_eq!(grid.value(1, 0, 0), 0.0);
    }

    #[test]
    fn test_complex_modes_are_refused_by_name() {
        for mode in [3, 4] {
            let bytes = orthogonal_header(1, 1, 1, mode, 1, 2, 3);
            let err = parse_ccp4(&bytes).unwrap_err();
            assert!(
                matches!(err, Ccp4Error::UnsupportedMode(m) if m == mode),
                "{err}"
            );
        }
    }

    #[test]
    fn test_nonzero_nsymbt_is_skipped_before_the_density_data() {
        let mut bytes = orthogonal_header(1, 1, 1, 2, 1, 2, 3);
        bytes[word_offset(24)..word_offset(24) + 4].copy_from_slice(&8i32.to_le_bytes());
        bytes.extend_from_slice(&[0u8; 8]); // the symmetry block itself, discarded
        append_f32_values(&mut bytes, &[77.0], false);

        let grid = parse_ccp4(&bytes).expect("valid CCP4");
        assert_eq!(grid.value(0, 0, 0), 77.0);
    }

    #[test]
    fn test_ccp4_rms_is_standard_deviation_around_the_mean() {
        // Confirmed against gemmi's own written RMS header word during
        // this story's research: values [100,200,300,400] (plus 20 zeros
        // in a 24-point grid) give RMS = 103.749..., not sqrt(mean(x^2)).
        let mut values = vec![0.0; 20];
        values.extend_from_slice(&[100.0, 200.0, 300.0, 400.0]);
        let rms = ccp4_rms(&values);
        assert!((rms - 103.749).abs() < 1e-2, "{rms}");
    }
}
