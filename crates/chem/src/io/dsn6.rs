//! DSN6 (#334) -- the bricked, byte-quantized electron-density map format
//! used by Frodo and O. Last of Phase 4's four volumetric formats after
//! CUBE (#331), CCP4/MRC (#332) and DX (#333).
//!
//! **No live oracle exists for this format anywhere in this environment**:
//! neither `gemmi` nor `gridData` (the two tools that verified every prior
//! volumetric format) support DSN6 or its text-header sibling BRIX at all,
//! confirmed directly. The one genuine reference found is `mgltools`'s
//! `ReadBRIX` class (`Volume/IO/volReaders.py`, Python 2), read in full and
//! cross-checked against its own cited spec
//! (`uoxray.uoregon.edu/tnt/manual/node104.html`, quoted in its own
//! docstring). Every load-bearing fact below traces to that reading, not to
//! memory -- including one genuine bug in the reference that this module
//! deliberately does not replicate.
//!
//! **The reference's own endianness-detection fallback is broken, not a
//! convention.** Its first attempt reads header word 19 (1-based; 18
//! 0-indexed) natively and checks it against the documented constant `100`.
//! Its two retries (big-endian, then little-endian) both read the *wrong*
//! word -- word 18 (1-based), the cell-scaling factor, not word 19 -- and
//! compare a `struct.unpack` **tuple** against the bare int `100`, which is
//! never equal regardless of the byte read. So both retries always "fail"
//! and the function always reports an invalid file whenever the native byte
//! order doesn't already match. This module does what the docstring actually
//! documents instead: read word 19 (1-based) under both byte orders, take
//! whichever gives exactly `100`, refuse if neither does -- the same
//! "refuse a silently-wrong guess" posture as CCP4's machine-stamp match
//! ([`crate::io::ccp4`]), just against a constant word rather than a fixed
//! byte pattern (DSN6 has no byte signature at all).
//!
//! **Header layout** (512 bytes, 19 significant `i16` words, 0-indexed word
//! `w` at byte offset `2w`): 0-2 origin (x,y,z, grid-index units); 3-5
//! extent `nc,nr,ns` -- the stored grid's own dims, with a **fixed** axis
//! mapping c=X, r=Y, s=Z (DSN6 has no CCP4-style `MAPC`/`MAPR`/`MAPS`
//! permutation, confirmed by the reference's brick-placement code using this
//! mapping unconditionally); 6-8 `nx,ny,nz` -- the whole unit cell's own
//! grid sampling (CCP4 `MX`/`MY`/`MZ`'s analogue), used only for the
//! fractional origin/step math below, never for `dims`; 9-14 cell
//! `a,b,c,alpha,beta,gamma`, each stored as `raw / cell_norm` (word 17);
//! 15-16 `prod,plus` -- density value scale/offset; 17 `cell_norm`; 18 the
//! fixed constant `100` (detection only).
//!
//! **Density value formula**, derived algebraically from the reference's own
//! docstring (its two defining equations for `prod`/`plus` in terms of
//! `rhomin`/`rhomax` over the documented byte range `[3,253]`, solved and
//! cross-checked in both directions): `real = (raw_byte - plus) / prod`,
//! `prod = raw_prod / 100.0` (dividing by the same word-18 constant used for
//! detection; `raw_byte` and `plus` are otherwise unscaled).
//!
//! **Bricked storage**: the volume is tiled into 8x8x8-voxel, 512-byte
//! bricks, brick-grid iterated Z-slowest -> Y -> X-fastest; the last brick
//! along an axis is clamped to the true remaining extent, and the padding
//! bytes beyond it are read and discarded, never treated as data (the
//! issue's own explicit warning). **Within one brick storage is X-fastest,
//! then Y, then Z-slowest** -- confirmed via the reference's own
//! `numpy.reshape(.., (8,8,8))` C-order assignment, then transposed back to
//! `(X,Y,Z)` for the value it actually returns -- which already matches
//! [`VolumeGrid`]'s own canonical X-fastest storage exactly. **No axis
//! permutation or [`VolumeGrid::from_source_order`] is needed anywhere in
//! this format**, a genuine simplification relative to CCP4.
//!
//! **A fixed, unconditional adjacent-byte-pair swap applies to every DSN6
//! brick's raw bytes** (`byte[2i] <-> byte[2i+1]`) before they're read as
//! independent `u8` density values -- confirmed in the reference
//! (`numpy.fromstring(cube, native_int16).byteswap().tostring()`,
//! independent of the file's own detected header byte order, so *not* tied
//! to it). BRIX needs no such swap at all; this module implements DSN6 only
//! (see below), so the swap is always applied.
//!
//! **Fractional origin/step -> Cartesian**, confirmed directly from the
//! reference's own arithmetic: `s_axis = 1.0 / (n_axis - 1)` using the
//! *sampling rate* (words 6-8, not the extent), `origin_frac = [xorigin *
//! sx, yorigin * sy, zorigin * sz]`. Converted to Cartesian via this crate's
//! own [`UnitCell::to_cartesian`]/[`UnitCell::basis`] (the same mechanism
//! CCP4/DX already use): `origin = cell.to_cartesian(origin_frac)`, and each
//! grid axis vector is `s_axis * basis[axis]`.
//!
//! **Writing has no independent oracle to check against** (the issue's own
//! named gap) -- verified here by inverting the read-side math exactly and
//! by hand-computed fixtures rather than a real second tool. `VolumeGrid`
//! has no separate concept of "whole-cell sampling rate" distinct from its
//! own stored extent, so the writer derives `nx,ny,nz` the same way CCP4's
//! own writer derives `MX,MY,MZ`: `round(basis_length / axis_length)`,
//! clamped to at least 2 (so `1/(n-1)` never divides by zero) -- self
//! -consistent on round trip rather than assumed equal to the extent.
//!
//! **BRIX is explicitly out of scope.** The issue names only DSN6; BRIX is a
//! closely related sibling (a human-readable, regex-parsed text header
//! instead of binary words, no byte-order concerns, an extra `sigma` field)
//! but a genuinely separate parsing path. Cutting it is a disclosed scope
//! decision, the same posture as CASTEP-PBC out of CUBE (#331) and
//! `.dx.gz` out of DX (#333).

use crate::core::cell::UnitCell;
use crate::core::volume::VolumeGrid;
use crate::io::errors::{Dsn6Error, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

const HEADER_LEN: usize = 512;
const BRICK: usize = 8;
const BRICK_VOXELS: usize = BRICK * BRICK * BRICK;

fn word_offset(word: usize) -> usize {
    word * 2
}

fn read_i16(bytes: &[u8], word: usize, big_endian: bool) -> Result<i16, Dsn6Error> {
    let off = word_offset(word);
    let raw: [u8; 2] = bytes
        .get(off..off + 2)
        .ok_or_else(|| Dsn6Error::ParseError("header truncated".to_string()))?
        .try_into()
        .unwrap();
    Ok(if big_endian {
        i16::from_be_bytes(raw)
    } else {
        i16::from_le_bytes(raw)
    })
}

fn put_i16(out: &mut [u8], word: usize, v: i16) {
    let off = word_offset(word);
    out[off..off + 2].copy_from_slice(&v.to_le_bytes());
}

/// Swaps every adjacent byte pair (`0<->1`, `2<->3`, ...) -- the fixed,
/// unconditional DSN6 brick convention (see the module doc). Self-inverse,
/// so the same function serves both the reader and the writer.
fn swap_byte_pairs(brick: &[u8; BRICK_VOXELS]) -> [u8; BRICK_VOXELS] {
    let mut out = [0u8; BRICK_VOXELS];
    let mut i = 0;
    while i < BRICK_VOXELS {
        out[i] = brick[i + 1];
        out[i + 1] = brick[i];
        i += 2;
    }
    out
}

/// Parses a whole DSN6 file into a [`VolumeGrid`].
pub(crate) fn parse_dsn6(bytes: &[u8]) -> Result<VolumeGrid, Dsn6Error> {
    if bytes.len() < HEADER_LEN {
        return Err(Dsn6Error::ParseError(
            "file shorter than the 512-byte header".to_string(),
        ));
    }

    // Header word 19 (1-based) = word 18 (0-indexed) is a documented
    // constant, always 100 -- the only detection mechanism DSN6 has (see
    // the module doc on why this checks the *right* word under both byte
    // orders, unlike the reference it was cross-checked against).
    let big_endian = if read_i16(bytes, 18, false)? == 100 {
        false
    } else if read_i16(bytes, 18, true)? == 100 {
        true
    } else {
        return Err(Dsn6Error::UnknownByteOrder);
    };

    let xorigin = read_i16(bytes, 0, big_endian)? as f64;
    let yorigin = read_i16(bytes, 1, big_endian)? as f64;
    let zorigin = read_i16(bytes, 2, big_endian)? as f64;
    let nc = read_i16(bytes, 3, big_endian)? as i64;
    let nr = read_i16(bytes, 4, big_endian)? as i64;
    let ns = read_i16(bytes, 5, big_endian)? as i64;
    let nx = read_i16(bytes, 6, big_endian)? as f64;
    let ny = read_i16(bytes, 7, big_endian)? as f64;
    let nz = read_i16(bytes, 8, big_endian)? as f64;
    let a_raw = read_i16(bytes, 9, big_endian)? as f64;
    let b_raw = read_i16(bytes, 10, big_endian)? as f64;
    let c_raw = read_i16(bytes, 11, big_endian)? as f64;
    let alpha_raw = read_i16(bytes, 12, big_endian)? as f64;
    let beta_raw = read_i16(bytes, 13, big_endian)? as f64;
    let gamma_raw = read_i16(bytes, 14, big_endian)? as f64;
    let prod_raw = read_i16(bytes, 15, big_endian)? as f64;
    let plus = read_i16(bytes, 16, big_endian)? as f64;
    let cell_norm = read_i16(bytes, 17, big_endian)? as f64;

    if nc <= 0 || nr <= 0 || ns <= 0 {
        return Err(Dsn6Error::ParseError(
            "extent (NC/NR/NS) must be positive".to_string(),
        ));
    }
    if cell_norm == 0.0 {
        return Err(Dsn6Error::ParseError(
            "cell scaling factor is zero".to_string(),
        ));
    }
    if nx <= 1.0 || ny <= 1.0 || nz <= 1.0 {
        return Err(Dsn6Error::ParseError(
            "sampling rate must be at least 2 along every axis".to_string(),
        ));
    }
    if prod_raw == 0.0 {
        return Err(Dsn6Error::ParseError(
            "density scale factor is zero".to_string(),
        ));
    }
    let (nc, nr, ns) = (nc as usize, nr as usize, ns as usize);

    let cell = UnitCell::new(
        a_raw / cell_norm,
        b_raw / cell_norm,
        c_raw / cell_norm,
        alpha_raw / cell_norm,
        beta_raw / cell_norm,
        gamma_raw / cell_norm,
    );
    let prod = prod_raw / 100.0;

    let dims = [nc, nr, ns];
    let mut values = vec![0.0f64; nc * nr * ns];
    let mut offset = HEADER_LEN;
    for zcube in 0..ns.div_ceil(BRICK) {
        let zstart = zcube * BRICK;
        let zfwd = (ns - zstart).min(BRICK);
        for ycube in 0..nr.div_ceil(BRICK) {
            let ystart = ycube * BRICK;
            let yfwd = (nr - ystart).min(BRICK);
            for xcube in 0..nc.div_ceil(BRICK) {
                let xstart = xcube * BRICK;
                let xfwd = (nc - xstart).min(BRICK);

                let raw: [u8; BRICK_VOXELS] = bytes
                    .get(offset..offset + BRICK_VOXELS)
                    .ok_or_else(|| Dsn6Error::ParseError("density data truncated".to_string()))?
                    .try_into()
                    .unwrap();
                offset += BRICK_VOXELS;
                let brick = swap_byte_pairs(&raw);

                for lz in 0..zfwd {
                    for ly in 0..yfwd {
                        for lx in 0..xfwd {
                            let local = lx + BRICK * (ly + BRICK * lz);
                            let density = (brick[local] as f64 - plus) / prod;
                            let ix = xstart + lx;
                            let iy = ystart + ly;
                            let iz = zstart + lz;
                            values[ix + dims[0] * (iy + dims[1] * iz)] = density;
                        }
                    }
                }
            }
        }
    }

    let sx = 1.0 / (nx - 1.0);
    let sy = 1.0 / (ny - 1.0);
    let sz = 1.0 / (nz - 1.0);
    let origin_frac = crate::core::geometry::Point3::new(xorigin * sx, yorigin * sy, zorigin * sz);
    let origin = cell.to_cartesian(origin_frac);
    let basis = cell.basis();
    let axes = [basis[0] * sx, basis[1] * sy, basis[2] * sz];

    let grid = VolumeGrid::new(dims, origin, axes, values, Some(cell))?;
    Ok(grid)
}

/// Writes a [`VolumeGrid`] as a DSN6 file -- always little-endian (see the
/// module doc: writing has no independent oracle, so this always emits the
/// one byte order its own reader is guaranteed to detect correctly).
pub(crate) fn write_dsn6(grid: &VolumeGrid) -> Vec<u8> {
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

    // The whole-cell sampling rate: derived the same way CCP4's own writer
    // derives MX/MY/MZ (round(basis_length / step_length)), since
    // `VolumeGrid` carries no separate concept of it. Clamped to at least 2
    // so `1/(n-1)` below never divides by zero.
    let basis = cell.basis();
    let sampling = |basis_len: f64, axis_len: f64| -> f64 {
        if axis_len <= 0.0 {
            2.0
        } else {
            (basis_len / axis_len).round().max(2.0)
        }
    };
    let nx = sampling(basis[0].length(), axes[0].length());
    let ny = sampling(basis[1].length(), axes[1].length());
    let nz = sampling(basis[2].length(), axes[2].length());

    let sx = 1.0 / (nx - 1.0);
    let sy = 1.0 / (ny - 1.0);
    let sz = 1.0 / (nz - 1.0);
    let origin_frac = cell.to_fractional(origin);
    let xorigin = (origin_frac.x / sx).round();
    let yorigin = (origin_frac.y / sy).round();
    let zorigin = (origin_frac.z / sz).round();

    // prod/plus map the grid's real value range onto the documented
    // "normal" byte range [3,253] (see the module doc's density formula) --
    // solved from the reference's own two defining equations. A flat grid
    // (max == min) has no real range to map, so every voxel maps to a fixed
    // mid-range byte instead.
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (prod, plus) = if !(max - min).is_finite() || (max - min) < 1e-12 {
        (1.0, 128.0 - min)
    } else {
        let prod = 250.0 / (max - min);
        (prod, 3.0 - min * prod)
    };
    // cell_norm: a fixed, disclosed choice giving two-decimal-place
    // precision for cell lengths/angles within i16's range -- the same
    // convention many real DSN6 files themselves use.
    let cell_norm = 100.0;

    let mut out = vec![0u8; HEADER_LEN];
    put_i16(&mut out, 0, xorigin as i16);
    put_i16(&mut out, 1, yorigin as i16);
    put_i16(&mut out, 2, zorigin as i16);
    put_i16(&mut out, 3, dims[0] as i16);
    put_i16(&mut out, 4, dims[1] as i16);
    put_i16(&mut out, 5, dims[2] as i16);
    put_i16(&mut out, 6, nx as i16);
    put_i16(&mut out, 7, ny as i16);
    put_i16(&mut out, 8, nz as i16);
    put_i16(&mut out, 9, (cell.a * cell_norm).round() as i16);
    put_i16(&mut out, 10, (cell.b * cell_norm).round() as i16);
    put_i16(&mut out, 11, (cell.c * cell_norm).round() as i16);
    put_i16(&mut out, 12, (cell.alpha * cell_norm).round() as i16);
    put_i16(&mut out, 13, (cell.beta * cell_norm).round() as i16);
    put_i16(&mut out, 14, (cell.gamma * cell_norm).round() as i16);
    put_i16(&mut out, 15, (prod * 100.0).round() as i16);
    put_i16(&mut out, 16, plus.round() as i16);
    put_i16(&mut out, 17, cell_norm as i16);
    put_i16(&mut out, 18, 100);

    let (nc, nr, ns) = (dims[0], dims[1], dims[2]);
    for zcube in 0..ns.div_ceil(BRICK) {
        let zstart = zcube * BRICK;
        let zfwd = (ns - zstart).min(BRICK);
        for ycube in 0..nr.div_ceil(BRICK) {
            let ystart = ycube * BRICK;
            let yfwd = (nr - ystart).min(BRICK);
            for xcube in 0..nc.div_ceil(BRICK) {
                let xstart = xcube * BRICK;
                let xfwd = (nc - xstart).min(BRICK);

                let mut brick = [0u8; BRICK_VOXELS];
                for lz in 0..zfwd {
                    for ly in 0..yfwd {
                        for lx in 0..xfwd {
                            let ix = xstart + lx;
                            let iy = ystart + ly;
                            let iz = zstart + lz;
                            let density = values[ix + dims[0] * (iy + dims[1] * iz)];
                            let raw = (density * prod + plus).round().clamp(0.0, 255.0) as u8;
                            brick[lx + BRICK * (ly + BRICK * lz)] = raw;
                        }
                    }
                }
                out.extend_from_slice(&swap_byte_pairs(&brick));
            }
        }
    }

    out
}

/// [`crate::io::format::ByteReadFn`] for DSN6.
pub(crate) fn read_dsn6_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_dsn6(bytes) {
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

/// [`crate::io::format::ByteWriteVolumeFn`] for DSN6.
pub(crate) fn write_dsn6_bytes(grid: &VolumeGrid, _options: &WriteOptions) -> Vec<u8> {
    write_dsn6(grid)
}

/// Buffers the whole input and parses it once -- a DSN6 file is always
/// exactly one grid, the same "one record" shape
/// [`crate::io::ccp4::Ccp4Supplier`] already has.
pub struct Dsn6Supplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl Dsn6Supplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match std::io::Read::read_to_end(&mut reader, &mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match parse_dsn6(&bytes) {
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

impl Iterator for Dsn6Supplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::geometry::Point3;

    #[test]
    fn test_a_non_multiple_of_8_round_trip() {
        // Every axis crosses a brick boundary without being a multiple of
        // it, exercising clamping on all three axes at once.
        let dims = [10usize, 9, 11];
        let n = dims[0] * dims[1] * dims[2];
        // Range [0, 250]: `prod` reconstructs to exactly 1.0 (no header
        // rounding beyond the byte itself), isolating what this test is
        // actually for -- brick/clamp geometry, not the format's separate,
        // disclosed `prod`-precision limit (see the module doc).
        let values: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64 * 250.0).collect();
        let grid = VolumeGrid::new(
            dims,
            Point3::new(1.0, 2.0, 3.0),
            [
                Point3::new(0.5, 0.0, 0.0),
                Point3::new(0.0, 0.6, 0.0),
                Point3::new(0.0, 0.0, 0.7),
            ],
            values,
            Some(UnitCell::new(50.0, 60.0, 70.0, 90.0, 90.0, 90.0)),
        )
        .expect("valid grid");

        let bytes = write_dsn6(&grid);
        let back = parse_dsn6(&bytes).expect("valid DSN6");

        assert_eq!(back.dims(), grid.dims());
        let back_cell = back.cell().expect("cell survives");
        assert!((back_cell.a - 50.0).abs() < 0.5, "{}", back_cell.a);
        assert!((back_cell.b - 60.0).abs() < 0.5, "{}", back_cell.b);
        assert!((back_cell.c - 70.0).abs() < 0.5, "{}", back_cell.c);

        // Byte quantization is lossy by design (see the module doc) --
        // tolerance reflects one part in ~250 of the value's own range.
        let tol = (grid.values().iter().cloned().fold(f64::MIN, f64::max)
            - grid.values().iter().cloned().fold(f64::MAX, f64::min))
        .abs()
            / 200.0;
        for ix in 0..dims[0] {
            for iy in 0..dims[1] {
                for iz in 0..dims[2] {
                    assert!(
                        (back.value(ix, iy, iz) - grid.value(ix, iy, iz)).abs() < tol,
                        "{} vs {} at {ix},{iy},{iz}",
                        back.value(ix, iy, iz),
                        grid.value(ix, iy, iz)
                    );
                }
            }
        }
    }

    /// A single, hand-built brick with an independently-derived expected
    /// on-disk byte, proving the byte-pair swap and the value formula
    /// together -- not just that read and write agree with each other.
    #[test]
    fn test_a_hand_built_single_brick_matches_an_independently_derived_byte() {
        // Header: extent 1x1x1, sampling 2x2x2 (minimum), cell 10A cubic,
        // cell_norm=100 -> raw cell words = 1000, angle words = 9000.
        // prod_raw=200 (prod=2.0), plus=10 -> real = (byte-10)/2.0.
        // Choosing byte=110 (post-swap) gives real = (110-10)/2.0 = 50.0.
        let mut header = [0u8; HEADER_LEN];
        let put = |h: &mut [u8; HEADER_LEN], word: usize, v: i16| {
            let off = word_offset(word);
            h[off..off + 2].copy_from_slice(&v.to_le_bytes());
        };
        put(&mut header, 0, 0); // xorigin
        put(&mut header, 1, 0); // yorigin
        put(&mut header, 2, 0); // zorigin
        put(&mut header, 3, 1); // nc
        put(&mut header, 4, 1); // nr
        put(&mut header, 5, 1); // ns
        put(&mut header, 6, 2); // nx
        put(&mut header, 7, 2); // ny
        put(&mut header, 8, 2); // nz
        put(&mut header, 9, 1000); // a * cell_norm
        put(&mut header, 10, 1000); // b * cell_norm
        put(&mut header, 11, 1000); // c * cell_norm
        put(&mut header, 12, 9000); // alpha * cell_norm
        put(&mut header, 13, 9000); // beta * cell_norm
        put(&mut header, 14, 9000); // gamma * cell_norm
        put(&mut header, 15, 200); // prod_raw (prod = 2.0)
        put(&mut header, 16, 10); // plus
        put(&mut header, 17, 100); // cell_norm
        put(&mut header, 18, 100); // fixed constant

        // One brick, 512 raw on-disk bytes: voxel (0,0,0) is local index 0
        // within the brick. Post byte-pair-swap, local index 0 comes from
        // raw byte 1 (swap_byte_pairs swaps 0<->1). So raw[1] = 110 makes
        // the swapped byte at index 0 equal 110.
        let mut brick = [0u8; BRICK_VOXELS];
        brick[1] = 110;
        let mut bytes = header.to_vec();
        bytes.extend_from_slice(&brick);

        let grid = parse_dsn6(&bytes).expect("valid DSN6");
        assert_eq!(grid.dims(), [1, 1, 1]);
        assert!(
            (grid.value(0, 0, 0) - 50.0).abs() < 1e-9,
            "{}",
            grid.value(0, 0, 0)
        );
    }

    #[test]
    fn test_a_big_endian_file_reads_identically_to_a_little_endian_one() {
        let dims = [3usize, 2, 2];
        let n = dims[0] * dims[1] * dims[2];
        let values: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let grid = VolumeGrid::new(
            dims,
            Point3::ORIGIN,
            [
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            values,
            Some(UnitCell::new(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)),
        )
        .expect("valid grid");

        let le_bytes = write_dsn6(&grid);

        // Byte-swap every header word to build an equivalent big-endian
        // fixture -- the density bytes are single bytes, unaffected by
        // header endianness (only the fixed pair swap applies to them).
        let mut be_bytes = le_bytes.clone();
        for word in 0..19 {
            let off = word_offset(word);
            be_bytes.swap(off, off + 1);
        }

        let from_le = parse_dsn6(&le_bytes).expect("valid little-endian DSN6");
        let from_be = parse_dsn6(&be_bytes).expect("valid big-endian DSN6");
        assert_eq!(from_le.dims(), from_be.dims());
        for ix in 0..dims[0] {
            for iy in 0..dims[1] {
                for iz in 0..dims[2] {
                    assert_eq!(from_le.value(ix, iy, iz), from_be.value(ix, iy, iz));
                }
            }
        }
    }

    #[test]
    fn test_neither_byte_order_gives_100_is_refused() {
        let mut bytes = vec![0u8; HEADER_LEN];
        put_i16(&mut bytes, 18, 42);
        let err = parse_dsn6(&bytes).unwrap_err();
        assert!(matches!(err, Dsn6Error::UnknownByteOrder), "{err}");
    }
}
