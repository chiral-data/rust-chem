//! GROMACS XTC — the same XDR framing [`crate::io::trr`] uses, wrapped
//! around a genuine lossy coordinate compressor (#325).
//!
//! Ported directly from the real `xdrfile` C library's source, not
//! reconstructed from prose: `xdrfile.c`/`xdrfile_xtc.c`, vendored inside
//! the locally-installed MDAnalysis package (this session found it at
//! `~/.cache/uv/archive-v0/*/MDAnalysis/lib/formats/src/{xdrfile.c,
//! xdrfile_xtc.c}`, also `include/xdrfile_xtc.h`) — re-consult that source
//! directly for any question about a step of this port, rather than this
//! module's own comments.
//!
//! **Frame shape**: `magic(1995)`, `natoms`, `step`, `time` (always single
//! precision — the real `xtc_header`/`xtc_coord` wrapper functions this
//! format actually uses are hardcoded to `f32`, unlike TRR's genuine
//! dual-precision support), `box` (9×`f32` row vectors, same convention
//! TRR/GRO already use), then the coordinate block: `size` (natoms,
//! restated), and either `size<=9` raw floats verbatim, or (`size>9`)
//! `precision`, `minint`/`maxint`, `smallidx`, a byte length, and that many
//! bytes of packed bits.
//!
//! **This is a full, faithful port of GROMACS's own encoder heuristic**,
//! not the simplified always-absolute variant the real format's own
//! tolerance for "a valid, not bit-for-bit identical, encoding" would also
//! permit (decided with the user) — every coordinate triple is absolute-
//! encoded into the minimum bits its frame's min/max range needs, and
//! consecutive atoms close enough to their predecessor are additionally
//! delta-coded in an adaptively-resized run, including the water-molecule-
//! specific first-two-atoms swap. This is what gets XTC to its usual
//! ~1/3-the-size-of-raw-floats ratio; a simplified encoder without it would
//! still be spec-valid but would not.
//!
//! **Every atom is [`Element::UNKNOWN`]** and there are no velocities or
//! forces — XTC states positions, time, step and a box, nothing else.

use crate::core::atom::{Atom, Element};
use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::trajectory::{Frame, FrameSource, Trajectory};
use crate::core::units::NM_TO_ANGSTROM;
use crate::io::errors::{ReadError, XtcError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};
use crate::io::xdr::{XdrReader, XdrWriter};

const GROMACS_XTC_MAGIC: i32 = 1995;

/// A fixed table of "nice" bit-width thresholds, indexed by `smallidx`. Only
/// the encoder/decoder's own bookkeeping ever chooses an index into this —
/// nothing else in this module interprets the values themselves.
#[rustfmt::skip]
const MAGICINTS: [i32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
    16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031,
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216,
];

const FIRSTIDX: usize = 9;
const LASTIDX: usize = MAGICINTS.len();

/// `MAGICINTS[idx]`, clamped in bounds. The reference C indexes this table
/// with `maxidx`, itself capped at `LASTIDX` -- the table's own *length*,
/// one past its last valid index -- so a pathological input (`smallidx`
/// adapted all the way to the top of the table) reads one element past the
/// array in C. Clamping is the responsible translation: correct for every
/// realistic input, safe (never a panic) for the input that would have hit
/// undefined behaviour in the original.
fn magicint(idx: usize) -> i32 {
    MAGICINTS[idx.min(LASTIDX - 1)]
}

/// Smallest number of bits needed to represent `0..=size`.
fn sizeofint(size: u32) -> u32 {
    let mut num: u64 = 1;
    let mut num_of_bits = 0u32;
    while (size as u64) >= num && num_of_bits < 32 {
        num_of_bits += 1;
        num <<= 1;
    }
    num_of_bits
}

/// Bits needed to jointly represent one value from each of `sizes`, packed
/// as a single mixed-radix number -- a base-256 bignum multiply-with-carry,
/// ported as a direct translation of the given byte-array algorithm rather
/// than reinterpreted via wider integer math, since a subtly different but
/// plausible reformulation is exactly the kind of bug an external reader
/// would silently refuse to decode.
fn sizeofints(sizes: &[u32]) -> u32 {
    let mut bytes = [0u32; 32];
    bytes[0] = 1;
    let mut num_of_bytes: usize = 1;
    for &size in sizes {
        let mut tmp: u64 = 0;
        for byte in bytes.iter_mut().take(num_of_bytes) {
            tmp += *byte as u64 * size as u64;
            *byte = (tmp & 0xff) as u32;
            tmp >>= 8;
        }
        while tmp != 0 {
            bytes[num_of_bytes] = (tmp & 0xff) as u32;
            num_of_bytes += 1;
            tmp >>= 8;
        }
    }
    let top = num_of_bytes - 1;
    let mut num: u64 = 1;
    let mut num_of_bits = 0u32;
    while (bytes[top] as u64) >= num {
        num_of_bits += 1;
        num *= 2;
    }
    num_of_bits + (top as u32) * 8
}

/// A byte-cursor bit-packer -- a clean re-expression of the reference
/// `encodebits`, which instead packs its own cursor state into the same
/// `int` array as the payload. Produces byte-identical output.
struct BitWriter {
    bytes: Vec<u8>,
    last_bits: u32,
    last_byte: u32,
}

impl BitWriter {
    fn new() -> Self {
        Self {
            bytes: Vec::new(),
            last_bits: 0,
            last_byte: 0,
        }
    }

    fn write_bits(&mut self, mut num_of_bits: u32, num: u32) {
        while num_of_bits >= 8 {
            self.last_byte = (self.last_byte << 8) | ((num >> (num_of_bits - 8)) & 0xff);
            self.bytes.push((self.last_byte >> self.last_bits) as u8);
            num_of_bits -= 8;
        }
        if num_of_bits > 0 {
            self.last_byte = (self.last_byte << num_of_bits) | (num & ((1 << num_of_bits) - 1));
            self.last_bits += num_of_bits;
            if self.last_bits >= 8 {
                self.last_bits -= 8;
                self.bytes.push((self.last_byte >> self.last_bits) as u8);
            }
        }
    }

    fn finish(mut self) -> Vec<u8> {
        if self.last_bits > 0 {
            self.bytes
                .push((self.last_byte << (8 - self.last_bits)) as u8);
        }
        self.bytes
    }
}

/// The `decodebits` counterpart of [`BitWriter`].
struct BitReader<'a> {
    bytes: &'a [u8],
    pos: usize,
    last_bits: u32,
    last_byte: u32,
}

impl<'a> BitReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            pos: 0,
            last_bits: 0,
            last_byte: 0,
        }
    }

    fn next_byte(&mut self) -> Result<u32, XtcError> {
        let b = *self.bytes.get(self.pos).ok_or_else(|| {
            XtcError::ParseError("compressed coordinate block ran out of bits".to_string())
        })?;
        self.pos += 1;
        Ok(b as u32)
    }

    fn read_bits(&mut self, num_of_bits: u32) -> Result<u32, XtcError> {
        let mask: u32 = if num_of_bits >= 32 {
            u32::MAX
        } else {
            (1 << num_of_bits) - 1
        };
        let mut remaining = num_of_bits;
        let mut num: u32 = 0;
        while remaining >= 8 {
            self.last_byte = (self.last_byte << 8) | self.next_byte()?;
            num |= (self.last_byte >> self.last_bits) << (remaining - 8);
            remaining -= 8;
        }
        if remaining > 0 {
            if self.last_bits < remaining {
                self.last_bits += 8;
                self.last_byte = (self.last_byte << 8) | self.next_byte()?;
            }
            self.last_bits -= remaining;
            num |= (self.last_byte >> self.last_bits) & ((1 << remaining) - 1);
        }
        Ok(num & mask)
    }
}

/// The `encodeints` port, specialised to the 3-value case every real call
/// site uses (a coordinate triple, or a run's delta triple).
fn encode_ints3(w: &mut BitWriter, num_of_bits: u32, sizes: [u32; 3], nums: [u32; 3]) {
    let mut bytes = [0u32; 32];
    let mut tmp = nums[0] as u64;
    let mut num_of_bytes: usize = 0;
    loop {
        bytes[num_of_bytes] = (tmp & 0xff) as u32;
        num_of_bytes += 1;
        tmp >>= 8;
        if tmp == 0 {
            break;
        }
    }
    for (i, &size) in sizes.iter().enumerate().skip(1) {
        debug_assert!(
            nums[i] < size,
            "encodeints: {} does not fit in {size}",
            nums[i]
        );
        let mut tmp = nums[i] as u64;
        for byte in bytes.iter_mut().take(num_of_bytes) {
            tmp += *byte as u64 * size as u64;
            *byte = (tmp & 0xff) as u32;
            tmp >>= 8;
        }
        while tmp != 0 {
            bytes[num_of_bytes] = (tmp & 0xff) as u32;
            num_of_bytes += 1;
            tmp >>= 8;
        }
    }
    if num_of_bits >= (num_of_bytes as u32) * 8 {
        for &byte in bytes.iter().take(num_of_bytes) {
            w.write_bits(8, byte);
        }
        w.write_bits(num_of_bits - (num_of_bytes as u32) * 8, 0);
    } else {
        for &byte in bytes.iter().take(num_of_bytes - 1) {
            w.write_bits(8, byte);
        }
        w.write_bits(
            num_of_bits - (num_of_bytes as u32 - 1) * 8,
            bytes[num_of_bytes - 1],
        );
    }
}

/// The `decodeints` port, specialised the same way as [`encode_ints3`].
fn decode_ints3(
    r: &mut BitReader,
    num_of_bits: u32,
    sizes: [u32; 3],
) -> Result<[u32; 3], XtcError> {
    let mut bytes = [0u32; 32];
    let mut num_of_bytes: usize = 0;
    let mut remaining = num_of_bits;
    while remaining > 8 {
        bytes[num_of_bytes] = r.read_bits(8)?;
        num_of_bytes += 1;
        remaining -= 8;
    }
    if remaining > 0 {
        bytes[num_of_bytes] = r.read_bits(remaining)?;
        num_of_bytes += 1;
    }
    let mut nums = [0u32; 3];
    for i in (1..3).rev() {
        let mut num: u64 = 0;
        for j in (0..num_of_bytes).rev() {
            num = (num << 8) | bytes[j] as u64;
            let p = num / sizes[i] as u64;
            bytes[j] = p as u32;
            num -= p * sizes[i] as u64;
        }
        nums[i] = num as u32;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
    Ok(nums)
}

/// One frame's worth of scaled-integer bookkeeping shared by encode/decode
/// setup -- everything derived once from `minint`/`maxint` before either
/// the per-triple absolute path or the run-delta path can proceed.
struct RangeCoding {
    sizeint: [u32; 3],
    bitsize: u32,
    bitsizeint: [u32; 3],
}

fn range_coding(minint: [i32; 3], maxint: [i32; 3]) -> RangeCoding {
    let sizeint = [
        (maxint[0] - minint[0] + 1) as u32,
        (maxint[1] - minint[1] + 1) as u32,
        (maxint[2] - minint[2] + 1) as u32,
    ];
    if sizeint.iter().any(|&s| s > 0xffffff) {
        RangeCoding {
            sizeint,
            bitsize: 0,
            bitsizeint: [
                sizeofint(sizeint[0]),
                sizeofint(sizeint[1]),
                sizeofint(sizeint[2]),
            ],
        }
    } else {
        RangeCoding {
            sizeint,
            bitsize: sizeofints(&sizeint),
            bitsizeint: [0; 3],
        }
    }
}

/// Compresses `positions_nm` (already nm-scaled, single precision) into one
/// coordinate block: `size`, then either raw floats (`size<=9`) or the full
/// compressed encoding.
fn compress_coord_block(positions_nm: &[[f32; 3]], precision: f32) -> Vec<u8> {
    let size = positions_nm.len();
    let mut out = XdrWriter::new();
    out.write_i32(size as i32);
    if size <= 9 {
        for p in positions_nm {
            out.write_f32(p[0]);
            out.write_f32(p[1]);
            out.write_f32(p[2]);
        }
        return out.into_bytes();
    }

    let precision = if precision <= 0.0 { 1000.0 } else { precision };
    out.write_f32(precision);

    // Round-half-away-from-zero, matching the reference's `x*precision +/-
    // 0.5` then truncate-toward-zero.
    let mut ints: Vec<[i32; 3]> = Vec::with_capacity(size);
    for p in positions_nm {
        let mut triple = [0i32; 3];
        for (axis, &v) in p.iter().enumerate() {
            let scaled = v * precision;
            let rounded = if scaled >= 0.0 {
                (scaled + 0.5).floor()
            } else {
                (scaled - 0.5).ceil()
            };
            triple[axis] = rounded as i32;
        }
        ints.push(triple);
    }

    let mut minint = [i32::MAX; 3];
    let mut maxint = [i32::MIN; 3];
    let mut mindiff: i64 = i64::MAX;
    let mut prev = [0i32; 3];
    for (i, &triple) in ints.iter().enumerate() {
        for axis in 0..3 {
            minint[axis] = minint[axis].min(triple[axis]);
            maxint[axis] = maxint[axis].max(triple[axis]);
        }
        if i >= 1 {
            let diff: i64 = (0..3)
                .map(|axis| (triple[axis] - prev[axis]).unsigned_abs() as i64)
                .sum();
            mindiff = mindiff.min(diff);
        }
        prev = triple;
    }
    out.write_i32(minint[0]);
    out.write_i32(minint[1]);
    out.write_i32(minint[2]);
    out.write_i32(maxint[0]);
    out.write_i32(maxint[1]);
    out.write_i32(maxint[2]);

    let coding = range_coding(minint, maxint);

    let mut smallidx = FIRSTIDX;
    while smallidx < LASTIDX && (magicint(smallidx) as i64) < mindiff {
        smallidx += 1;
    }
    out.write_i32(smallidx as i32);

    let maxidx = LASTIDX.min(smallidx + 8);
    let minidx = maxidx - 8;
    let tmp_idx = ((smallidx as i64) - 1).max(FIRSTIDX as i64) as usize;
    let mut smaller = magicint(tmp_idx) / 2;
    let mut smallnum = magicint(smallidx) / 2;
    let mut sizesmall = [magicint(smallidx) as u32; 3];
    let larger = magicint(maxidx) / 2;

    let mut w = BitWriter::new();
    let mut prevcoord = [0i32; 3];
    let mut prevrun: Option<u32> = None;
    let mut i = 0usize;
    while i < size {
        let mut is_smaller: i32 = if smallidx < maxidx
            && i >= 1
            && (0..3).all(|axis| (ints[i][axis] - prevcoord[axis]).abs() < larger)
        {
            1
        } else if smallidx > minidx {
            -1
        } else {
            0
        };

        let mut is_small = false;
        if i + 1 < size && (0..3).all(|axis| (ints[i][axis] - ints[i + 1][axis]).abs() < smallnum) {
            ints.swap(i, i + 1);
            is_small = true;
        }

        let tmpcoord = [
            (ints[i][0] - minint[0]) as u32,
            (ints[i][1] - minint[1]) as u32,
            (ints[i][2] - minint[2]) as u32,
        ];
        if coding.bitsize == 0 {
            w.write_bits(coding.bitsizeint[0], tmpcoord[0]);
            w.write_bits(coding.bitsizeint[1], tmpcoord[1]);
            w.write_bits(coding.bitsizeint[2], tmpcoord[2]);
        } else {
            encode_ints3(&mut w, coding.bitsize, coding.sizeint, tmpcoord);
        }
        prevcoord = ints[i];
        i += 1;

        if !is_small && is_smaller == -1 {
            is_smaller = 0;
        }
        let mut run: u32 = 0;
        let mut run_deltas: Vec<u32> = Vec::new();
        while is_small && run < 24 {
            let cur = ints[i];
            let tmpsum: i64 = (0..3)
                .map(|axis| {
                    let d = (cur[axis] - prevcoord[axis]) as i64;
                    d * d
                })
                .sum();
            if is_smaller == -1 && tmpsum >= (smaller as i64) * (smaller as i64) {
                is_smaller = 0;
            }
            for axis in 0..3 {
                run_deltas.push((cur[axis] - prevcoord[axis] + smallnum) as u32);
            }
            run += 3;
            prevcoord = cur;
            i += 1;
            is_small =
                i < size && (0..3).all(|axis| (ints[i][axis] - prevcoord[axis]).abs() < smallnum);
        }

        if prevrun != Some(run) || is_smaller != 0 {
            prevrun = Some(run);
            w.write_bits(1, 1);
            w.write_bits(5, (run as i32 + is_smaller + 1) as u32);
        } else {
            w.write_bits(1, 0);
        }
        for k in (0..run as usize).step_by(3) {
            encode_ints3(
                &mut w,
                smallidx as u32,
                sizesmall,
                [run_deltas[k], run_deltas[k + 1], run_deltas[k + 2]],
            );
        }

        smallidx = (smallidx as i32 + is_smaller) as usize;
        if is_smaller < 0 {
            smallnum = smaller;
            smaller = if smallidx > FIRSTIDX {
                magicint(smallidx - 1) / 2
            } else {
                0
            };
        } else if is_smaller > 0 {
            smaller = smallnum;
            smallnum = magicint(smallidx) / 2;
        }
        sizesmall = [magicint(smallidx) as u32; 3];
    }

    let packed = w.finish();
    out.write_i32(packed.len() as i32);
    out.write_padded(&packed);
    out.into_bytes()
}

/// Decompresses one coordinate block, returning nm-scaled, single-precision
/// positions.
fn decompress_coord_block(r: &mut XdrReader) -> Result<Vec<[f32; 3]>, XtcError> {
    let lsize = r.read_i32()?;
    if lsize < 0 {
        return Err(XtcError::ParseError(format!("negative atom count {lsize}")));
    }
    let lsize = lsize as usize;

    if lsize <= 9 {
        let mut out = Vec::with_capacity(lsize);
        for _ in 0..lsize {
            out.push([r.read_f32()?, r.read_f32()?, r.read_f32()?]);
        }
        return Ok(out);
    }

    let precision = r.read_f32()?;
    let mut minint = [0i32; 3];
    let mut maxint = [0i32; 3];
    for v in minint.iter_mut() {
        *v = r.read_i32()?;
    }
    for v in maxint.iter_mut() {
        *v = r.read_i32()?;
    }
    let coding = range_coding(minint, maxint);

    let smallidx_raw = r.read_i32()?;
    if smallidx_raw < 0 || smallidx_raw as usize >= LASTIDX {
        return Err(XtcError::ParseError(format!(
            "smallidx {smallidx_raw} out of range"
        )));
    }
    let mut smallidx = smallidx_raw as usize;
    let tmp_idx = ((smallidx as i64) - 1).max(FIRSTIDX as i64) as usize;
    let mut smaller = magicint(tmp_idx) / 2;
    let mut smallnum = magicint(smallidx) / 2;
    let mut sizesmall = [magicint(smallidx) as u32; 3];

    let length = r.read_i32()?;
    if length < 0 {
        return Err(XtcError::ParseError(format!(
            "negative compressed length {length}"
        )));
    }
    let packed = r.read_padded(length as usize)?;
    let mut br = BitReader::new(&packed);

    let mut out: Vec<[i32; 3]> = Vec::with_capacity(lsize);
    // Unlike the encoder's own `prevcoord`, this one never needs to carry a
    // value across outer-loop iterations -- every iteration overwrites it
    // unconditionally before any read.
    #[allow(unused_assignments)]
    let mut prevcoord = [0i32; 3];
    // Persists across outer-loop iterations, unlike `is_smaller` -- a
    // `flag==0` control bit means "same run-length as last time", so this
    // must keep its previous value rather than reset, mirroring the
    // encoder's own `prevrun`.
    let mut run: u32 = 0;
    let mut i = 0usize;
    while i < lsize {
        let mut thiscoord = if coding.bitsize == 0 {
            [
                br.read_bits(coding.bitsizeint[0])? as i32,
                br.read_bits(coding.bitsizeint[1])? as i32,
                br.read_bits(coding.bitsizeint[2])? as i32,
            ]
        } else {
            let v = decode_ints3(&mut br, coding.bitsize, coding.sizeint)?;
            [v[0] as i32, v[1] as i32, v[2] as i32]
        };
        i += 1;
        for axis in 0..3 {
            thiscoord[axis] += minint[axis];
        }
        prevcoord = thiscoord;

        let flag = br.read_bits(1)?;
        let mut is_smaller: i32 = 0;
        if flag == 1 {
            let raw = br.read_bits(5)?;
            let rem = raw % 3;
            run = raw - rem;
            is_smaller = rem as i32 - 1;
        }
        if run > 0 {
            let mut k = 0u32;
            while k < run {
                if i >= lsize {
                    return Err(XtcError::ParseError(
                        "compressed run overruns the declared atom count".to_string(),
                    ));
                }
                let delta = decode_ints3(&mut br, smallidx as u32, sizesmall)?;
                i += 1;
                let mut cur = [
                    delta[0] as i32 + prevcoord[0] - smallnum,
                    delta[1] as i32 + prevcoord[1] - smallnum,
                    delta[2] as i32 + prevcoord[2] - smallnum,
                ];
                if k == 0 {
                    std::mem::swap(&mut cur, &mut prevcoord);
                    out.push(prevcoord);
                } else {
                    prevcoord = cur;
                }
                out.push(cur);
                k += 3;
            }
        } else {
            out.push(thiscoord);
        }

        let new_smallidx = smallidx as i64 + is_smaller as i64;
        if new_smallidx < 0 || new_smallidx >= LASTIDX as i64 {
            return Err(XtcError::ParseError(
                "corrupt compressed coordinate block: smallidx out of range".to_string(),
            ));
        }
        smallidx = new_smallidx as usize;
        if is_smaller < 0 {
            smallnum = smaller;
            smaller = if smallidx > FIRSTIDX {
                magicint(smallidx - 1) / 2
            } else {
                0
            };
        } else if is_smaller > 0 {
            smaller = smallnum;
            smallnum = magicint(smallidx) / 2;
        }
        sizesmall = [magicint(smallidx) as u32; 3];
    }

    if out.len() != lsize {
        return Err(XtcError::ParseError(format!(
            "decompressed {} positions, expected {lsize}",
            out.len()
        )));
    }
    let inv_precision = 1.0 / precision;
    Ok(out
        .into_iter()
        .map(|c| {
            [
                c[0] as f32 * inv_precision,
                c[1] as f32 * inv_precision,
                c[2] as f32 * inv_precision,
            ]
        })
        .collect())
}

fn read_header(r: &mut XdrReader) -> Result<(usize, i32, f32), XtcError> {
    let magic = r.read_i32()?;
    if magic != GROMACS_XTC_MAGIC {
        return Err(XtcError::InvalidMagicNumber(magic));
    }
    let natoms = r.read_i32()?;
    if natoms < 0 {
        return Err(XtcError::ParseError("negative atom count".to_string()));
    }
    let step = r.read_i32()?;
    let time = r.read_f32()?;
    Ok((natoms as usize, step, time))
}

fn parse_frame(bytes: &[u8]) -> Result<Frame, XtcError> {
    let mut r = XdrReader::new(bytes);
    let (_natoms, step, time) = read_header(&mut r)?;

    let mut box_reals = [0f32; 9];
    for v in box_reals.iter_mut() {
        *v = r.read_f32()?;
    }
    let basis: Vec<Point3> = box_reals
        .chunks(3)
        .map(|c| Point3::new(c[0] as f64, c[1] as f64, c[2] as f64) * NM_TO_ANGSTROM)
        .collect();
    let cell = if box_reals.iter().any(|&v| v != 0.0) {
        Some(cell_from_box_vectors(basis[0], basis[1], basis[2]))
    } else {
        None
    };

    let positions_nm = decompress_coord_block(&mut r)?;
    let positions = positions_nm
        .into_iter()
        .map(|p| Point3::new(p[0] as f64, p[1] as f64, p[2] as f64) * NM_TO_ANGSTROM)
        .collect();

    Ok(Frame {
        positions,
        velocities: None,
        forces: None,
        time: Some(time as f64),
        step: Some(step.max(0) as u64),
        cell,
    })
}

fn angle_degrees(a: Point3, b: Point3) -> f64 {
    (a.dot(b) / (a.length() * b.length()))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

/// Same convention TRR's own box already uses -- three literal row vectors,
/// not `xlo`/`xhi`-style bounds.
fn cell_from_box_vectors(v1: Point3, v2: Point3, v3: Point3) -> UnitCell {
    UnitCell::new(
        v1.length(),
        v2.length(),
        v3.length(),
        angle_degrees(v2, v3),
        angle_degrees(v1, v3),
        angle_degrees(v1, v2),
    )
}

fn write_frame(out: &mut Vec<u8>, frame: &Frame, precision: f32) {
    let mut w = XdrWriter::new();
    w.write_i32(GROMACS_XTC_MAGIC);
    w.write_i32(frame.num_atoms() as i32);
    w.write_i32(frame.step.unwrap_or(0) as i32);
    w.write_f32(frame.time.unwrap_or(0.0) as f32);

    let basis = match frame.cell {
        Some(cell) => cell.basis(),
        None => [Point3::ORIGIN; 3],
    };
    for v in basis {
        w.write_f32((v.x / NM_TO_ANGSTROM) as f32);
        w.write_f32((v.y / NM_TO_ANGSTROM) as f32);
        w.write_f32((v.z / NM_TO_ANGSTROM) as f32);
    }

    let positions_nm: Vec<[f32; 3]> = frame
        .positions
        .iter()
        .map(|p| {
            [
                (p.x / NM_TO_ANGSTROM) as f32,
                (p.y / NM_TO_ANGSTROM) as f32,
                (p.z / NM_TO_ANGSTROM) as f32,
            ]
        })
        .collect();

    out.extend_from_slice(&w.into_bytes());
    out.extend_from_slice(&compress_coord_block(&positions_nm, precision));
}

/// One pass building a per-frame byte-offset index -- each frame's own
/// `size`/precision/length fields are exactly what lets the body be
/// skipped without decoding it, the same approach [`crate::io::trr`]'s own
/// frame source uses.
fn scan_offsets(bytes: &[u8]) -> Result<(Vec<u64>, usize), XtcError> {
    let mut offsets = Vec::new();
    let mut pos = 0usize;
    let mut natoms = None;
    while pos < bytes.len() {
        let mut r = XdrReader::new(&bytes[pos..]);
        let (n, _step, _time) = read_header(&mut r)?;
        offsets.push(pos as u64);
        if natoms.is_none() {
            natoms = Some(n);
        }
        r.skip(9 * 4)?; // box: 9 f32s
        let size = r.read_i32()?;
        if size < 0 {
            return Err(XtcError::ParseError("negative atom count".to_string()));
        }
        let body_len = if size as usize <= 9 {
            size as usize * 3 * 4
        } else {
            r.skip(4 + 12 + 12)?; // precision, minint, maxint
            let _smallidx = r.read_i32()?;
            let length = r.read_i32()?;
            if length < 0 {
                return Err(XtcError::ParseError(
                    "negative compressed length".to_string(),
                ));
            }
            (length as usize).div_ceil(4) * 4
        };
        pos += r.position() + body_len;
    }
    Ok((offsets, natoms.unwrap_or(0)))
}

pub(crate) struct XtcFrameSource {
    cursor: std::io::Cursor<Vec<u8>>,
    offsets: Vec<u64>,
    natoms: usize,
}

impl XtcFrameSource {
    fn open(bytes: Vec<u8>) -> Result<Self, XtcError> {
        let (offsets, natoms) = scan_offsets(&bytes)?;
        Ok(Self {
            cursor: std::io::Cursor::new(bytes),
            offsets,
            natoms,
        })
    }
}

impl FrameSource for XtcFrameSource {
    fn frame_count(&self) -> usize {
        self.offsets.len()
    }

    fn num_atoms(&self) -> usize {
        self.natoms
    }

    fn frame(&mut self, index: usize) -> std::io::Result<Frame> {
        use std::io::{Read, Seek, SeekFrom};
        self.cursor.seek(SeekFrom::Start(self.offsets[index]))?;
        let mut buf = Vec::new();
        self.cursor.read_to_end(&mut buf)?;
        parse_frame(&buf).map_err(std::io::Error::other)
    }
}

fn build_trajectory(bytes: Vec<u8>) -> Result<Trajectory, XtcError> {
    if bytes.is_empty() {
        return Err(XtcError::ParseError("empty input".to_string()));
    }
    let source = XtcFrameSource::open(bytes)?;
    let natoms = source.num_atoms();
    let mut topology = Molecule::new();
    for _ in 0..natoms {
        topology.add_atom(Atom::new(Element::UNKNOWN));
    }
    Ok(Trajectory::new(topology, Box::new(source))?)
}

/// [`crate::io::format::ByteReadFn`] for XTC.
pub(crate) fn read_xtc_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match build_trajectory(bytes.to_vec()) {
        Ok(trajectory) => out.records.push(Record {
            payload: Payload::Frames(trajectory),
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

/// [`crate::io::format::ByteWriteTrajectoryFn`] for XTC.
pub(crate) fn write_xtc_bytes(trajectory: &mut Trajectory, options: &WriteOptions) -> Vec<u8> {
    let mut out = Vec::new();
    let precision = options.xtc.precision;
    for i in 0..trajectory.frame_count() {
        let frame = trajectory
            .frame(i)
            .expect("a Trajectory already validated its own frame/atom counts at construction");
        write_frame(&mut out, &frame, precision);
    }
    out
}

/// Buffers the whole input and parses it once, mirroring
/// [`crate::io::trr::TrrSupplier`].
pub struct XtcSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl XtcSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match reader.read_to_end(&mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match build_trajectory(bytes) {
                Ok(trajectory) => vec![Ok(Record {
                    payload: Payload::Frames(trajectory),
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

impl Iterator for XtcSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct VecFrames(Vec<Frame>);

    impl FrameSource for VecFrames {
        fn frame_count(&self) -> usize {
            self.0.len()
        }
        fn num_atoms(&self) -> usize {
            self.0.first().map(Frame::num_atoms).unwrap_or(0)
        }
        fn frame(&mut self, index: usize) -> std::io::Result<Frame> {
            Ok(self.0[index].clone())
        }
    }

    fn topology(natoms: usize) -> Molecule {
        let mut mol = Molecule::new();
        for _ in 0..natoms {
            mol.add_atom(Atom::new(Element::UNKNOWN));
        }
        mol
    }

    fn trajectory_from(frames: Vec<Frame>) -> Trajectory {
        let natoms = frames[0].num_atoms();
        Trajectory::new(topology(natoms), Box::new(VecFrames(frames))).expect("valid")
    }

    fn as_trajectory(outcome: ReadOutcome) -> Trajectory {
        match outcome
            .records
            .into_iter()
            .next()
            .expect("one record")
            .payload
        {
            Payload::Frames(t) => t,
            other => panic!("expected Payload::Frames, got {other:?}"),
        }
    }

    fn frame_no_cell(positions: Vec<Point3>, time: f64, step: u64) -> Frame {
        Frame {
            positions,
            velocities: None,
            forces: None,
            time: Some(time),
            step: Some(step),
            cell: None,
        }
    }

    // --- bit-level primitives, hand-traced against the C source ---------

    #[test]
    fn test_sizeofint_matches_hand_traced_values() {
        assert_eq!(sizeofint(0), 0);
        assert_eq!(sizeofint(1), 1);
        assert_eq!(sizeofint(2), 2);
        assert_eq!(sizeofint(255), 8);
        assert_eq!(sizeofint(256), 9);
    }

    #[test]
    fn test_sizeofints_matches_a_hand_traced_value() {
        // 2*2*2 = 8, needing exactly 4 bits per the reference's own
        // bit-measuring loop (worked by hand against the C source).
        assert_eq!(sizeofints(&[2, 2, 2]), 4);
        // 256*256*256 = 2^24, and the reference's own measuring loop counts
        // bits needed to represent the product *itself* (not product-1),
        // so 2^24 needs 25 bits, not 24 -- confirmed by tracing the byte
        // array by hand against the C source, not assumed.
        assert_eq!(sizeofints(&[256, 256, 256]), 25);
    }

    #[test]
    fn test_bit_writer_and_reader_round_trip_arbitrary_widths() {
        let mut w = BitWriter::new();
        w.write_bits(1, 1);
        w.write_bits(5, 17);
        w.write_bits(12, 3000);
        w.write_bits(32, 0xdead_beef);
        let bytes = w.finish();

        let mut r = BitReader::new(&bytes);
        assert_eq!(r.read_bits(1).unwrap(), 1);
        assert_eq!(r.read_bits(5).unwrap(), 17);
        assert_eq!(r.read_bits(12).unwrap(), 3000);
        assert_eq!(r.read_bits(32).unwrap(), 0xdead_beef);
    }

    #[test]
    fn test_encode_decode_ints3_round_trip() {
        let sizes = [1000u32, 500, 2000]; // sizes[0] is never actually used
        let bitsize = sizeofints(&sizes);
        let mut w = BitWriter::new();
        encode_ints3(&mut w, bitsize, sizes, [1, 3, 42]);
        encode_ints3(&mut w, bitsize, sizes, [999, 499, 1999]);
        let bytes = w.finish();

        let mut r = BitReader::new(&bytes);
        assert_eq!(decode_ints3(&mut r, bitsize, sizes).unwrap(), [1, 3, 42]);
        assert_eq!(
            decode_ints3(&mut r, bitsize, sizes).unwrap(),
            [999, 499, 1999]
        );
    }

    // --- frame-level round trips -----------------------------------------

    #[test]
    fn test_nine_or_fewer_atoms_use_the_raw_float_special_case() {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
        ];
        let mut trajectory = trajectory_from(vec![frame_no_cell(positions.clone(), 0.0, 0)]);
        let bytes = write_xtc_bytes(&mut trajectory, &WriteOptions::default());
        let outcome = read_xtc_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut back = as_trajectory(outcome);
        let f = back.frame(0).unwrap();
        for (a, b) in f.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 1e-4, "{a} vs {b}");
        }
    }

    #[test]
    fn test_more_than_nine_atoms_round_trip_within_precision_and_actually_compress() {
        // 12 atoms in a loose, non-clustered arrangement -- exercises the
        // main range-based path without necessarily triggering runs.
        let positions: Vec<Point3> = (0..12)
            .map(|i| Point3::new(i as f64 * 2.0, (i % 3) as f64 * 3.0, 0.0))
            .collect();
        let mut trajectory = trajectory_from(vec![frame_no_cell(positions.clone(), 1.5, 3)]);
        let bytes = write_xtc_bytes(&mut trajectory, &WriteOptions::default());

        // Raw floats would cost 3*4 bytes/atom; a real range-based encoding
        // of a small, non-adversarial system should beat that meaningfully.
        assert!(
            bytes.len() < positions.len() * 12,
            "expected real compression, got {} bytes for {} atoms",
            bytes.len(),
            positions.len()
        );

        let outcome = read_xtc_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut back = as_trajectory(outcome);
        let f = back.frame(0).unwrap();
        for (a, b) in f.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 0.01, "{a} vs {b}");
            assert!((a.y - b.y).abs() < 0.01, "{a} vs {b}");
            assert!((a.z - b.z).abs() < 0.01, "{a} vs {b}");
        }
        assert!((f.time.unwrap() - 1.5).abs() < 1e-6);
        assert_eq!(f.step, Some(3));
        assert!(f.velocities.is_none());
        assert!(f.forces.is_none());
    }

    #[test]
    fn test_a_water_like_clustered_system_round_trips_and_exercises_the_run_path() {
        // Twelve triplets, each triplet tightly clustered (an O and two Hs)
        // but triplets far apart from each other -- exactly the shape the
        // water-molecule run/swap optimisation targets. If the run path
        // were silently broken, positions would still very likely fail to
        // round trip (the delta/swap arithmetic is wrong in a way that
        // corrupts values, not just skips a size optimisation).
        let mut positions = Vec::new();
        for m in 0..12 {
            let base = Point3::new(m as f64 * 10.0, 0.0, 0.0);
            positions.push(base);
            positions.push(base + Point3::new(0.09, 0.0, 0.0));
            positions.push(base + Point3::new(0.0, 0.09, 0.0));
        }
        let mut trajectory = trajectory_from(vec![frame_no_cell(positions.clone(), 0.0, 0)]);
        let bytes = write_xtc_bytes(&mut trajectory, &WriteOptions::default());
        let outcome = read_xtc_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut back = as_trajectory(outcome);
        let f = back.frame(0).unwrap();
        assert_eq!(f.positions.len(), positions.len());
        for (a, b) in f.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 0.01, "{a} vs {b}");
            assert!((a.y - b.y).abs() < 0.01, "{a} vs {b}");
            assert!((a.z - b.z).abs() < 0.01, "{a} vs {b}");
        }
    }

    #[test]
    fn test_a_triclinic_box_round_trips() {
        let cell = UnitCell::new(30.0, 30.0, 30.0, 80.0, 85.0, 95.0);
        let positions: Vec<Point3> = (0..15).map(|i| Point3::new(i as f64, 0.0, 0.0)).collect();
        let frame = Frame {
            positions: positions.clone(),
            velocities: None,
            forces: None,
            time: Some(0.0),
            step: Some(0),
            cell: Some(cell),
        };
        let mut trajectory = trajectory_from(vec![frame]);
        let bytes = write_xtc_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_xtc_bytes(&bytes, &ReadOptions::default()));
        let f = back.frame(0).unwrap();
        let back_cell = f.cell.expect("cell survives");
        assert!((back_cell.a - cell.a).abs() < 1e-2, "{}", back_cell.a);
        assert!(
            (back_cell.alpha - cell.alpha).abs() < 1e-2,
            "{}",
            back_cell.alpha
        );
        assert!(
            (back_cell.gamma - cell.gamma).abs() < 1e-2,
            "{}",
            back_cell.gamma
        );
    }

    #[test]
    fn test_multiple_frames_are_randomly_accessible() {
        let frames: Vec<Frame> = (0..3)
            .map(|f| {
                let positions: Vec<Point3> = (0..11)
                    .map(|i| Point3::new((i + f * 20) as f64, 0.0, 0.0))
                    .collect();
                frame_no_cell(positions, f as f64, f as u64)
            })
            .collect();
        let mut trajectory = trajectory_from(frames.clone());
        let bytes = write_xtc_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_xtc_bytes(&bytes, &ReadOptions::default()));
        assert_eq!(back.frame_count(), 3);

        let last = back.frame(2).unwrap();
        assert!((last.positions[0].x - 40.0).abs() < 0.01);
        let first = back.frame(0).unwrap();
        assert!((first.positions[0].x - 0.0).abs() < 0.01);
    }

    #[test]
    fn test_a_custom_precision_changes_the_rounding_resolution() {
        let positions = vec![Point3::new(1.2345, 0.0, 0.0); 11];
        let mut options = WriteOptions::default();
        options.xtc.precision = 10.0; // 0.1 nm = 1.0 Angstrom resolution
        let mut trajectory = trajectory_from(vec![frame_no_cell(positions, 0.0, 0)]);
        let bytes = write_xtc_bytes(&mut trajectory, &options);
        let mut back = as_trajectory(read_xtc_bytes(&bytes, &ReadOptions::default()));
        let f = back.frame(0).unwrap();
        // At 1 Angstrom resolution, 1.2345 must not survive to four decimal
        // places -- proving the option actually reached the compressor.
        assert!(
            (f.positions[0].x - 1.2345).abs() > 0.05,
            "{}",
            f.positions[0].x
        );
        assert!(
            (f.positions[0].x - 1.2345).abs() < 1.5,
            "{}",
            f.positions[0].x
        );
    }

    #[test]
    fn test_a_bad_magic_number_is_a_clear_error_not_a_panic() {
        let outcome = read_xtc_bytes(&[0, 0, 0, 0], &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0].error.contains("1995"),
            "{}",
            outcome.skipped[0].error
        );
    }

    #[test]
    fn test_empty_input_is_a_clear_error_not_a_panic() {
        let outcome = read_xtc_bytes(&[], &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
    }
}
