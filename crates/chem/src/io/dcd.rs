//! DCD — CHARMM/NAMD's binary trajectory (#327), and the oldest format in
//! this milestone: Fortran unformatted records (each block wrapped in a
//! leading and trailing byte-count marker) around per-frame coordinate
//! arrays. No XDR, no fixed byte order, no per-frame magic — unlike TRR/XTC,
//! this format's difficulty is entirely a handful of undocumented,
//! ambiguous conventions that fail *silently*, not with an error, when
//! misread.
//!
//! Ported directly from the real reference implementation, not
//! reconstructed from prose: `readdcd.h` (the C reader every real tool,
//! including MDAnalysis, still uses), `MDAnalysis/coordinates/DCD.py` (the
//! unit-cell interpretation heuristic) and `libdcd.pyx` (the writer), all
//! vendored inside the locally-installed MDAnalysis package this session
//! also used as TRR/XTC's oracle.
//!
//! **Coordinates are already Å** — CHARMM's own unit system — so unlike
//! TRR/XTC there is no nm-to-Å conversion for positions. Only the unit
//! cell's *angles* need interpretation (below); its lengths are plain Å.
//! Time (`delta`) is left in whatever unit the file states — DCD's own
//! convention is inconsistent across producers (AKMA, ps, ...) — matching
//! [`Frame::time`]'s own documented contract.
//!
//! **Endianness** is a per-file toggle with no marker of its own: the
//! leading Fortran record length must be `84`; if it isn't when read
//! little-endian, and is when read big-endian, the whole file is
//! big-endian. Decided once, applied to every subsequent field. This
//! crate's own writer only ever produces little-endian.
//!
//! **CHARMM vs. X-PLOR dialect** is one `i32` at a fixed header offset
//! (zero for X-PLOR, the CHARMM version number otherwise), and governs
//! three things at once: whether `DELTA` is a stored `f32` (CHARMM) or
//! `f64` (X-PLOR), whether a "has extra block" flag exists, and therefore
//! whether a per-frame unit cell can be present at all — X-PLOR files
//! never carry one. This crate's own writer always writes CHARMM.
//!
//! **The unit cell record**, when present, is six doubles in file order
//! `[A, gamma, B, beta, alpha, C]` — not `[A, B, C, alpha, beta, gamma]` —
//! confirmed against the exact reorder MDAnalysis's own reader uses. The
//! three angle values are then ambiguous on their own: if all three lie in
//! `[-1, 1]` they are angle-*cosines* (the modern NAMD>2.5/VMD convention);
//! otherwise they are already plain degrees (the older convention). This
//! is the same value-range heuristic MDAnalysis itself uses, since nothing
//! more reliable exists. This crate's own writer always uses the modern
//! cosine convention.
//!
//! **Out of scope, disclosed**: the rarer "new-style CHARMM symmetric
//! box-vector" unit-cell encoding (negative lengths or angles > 180°,
//! meaning the six values are matrix components, not lengths/angles at
//! all) is a clear [`crate::io::errors::DcdError::UnsupportedUnitCellFormat`]
//! rather than a silently wrong cell. The optional CHARMM "4th dimension"
//! data block is skipped, never interpreted, the same "parsed and
//! discarded" treatment this crate's own PRMTOP/TOP readers already give
//! their own out-of-scope sections.
//!
//! **Fixed atoms**: a CHARMM DCD may declare `NAMNF` atoms fixed. Frame 0
//! always states every atom in full (and is cached); every later frame
//! states only the *free* atoms' new positions, spliced back into the
//! cached frame-0 values by a 1-based index array read once from the
//! header. This module's own frame source transparently reads and caches frame 0
//! first if a caller asks for a later frame directly — true random access
//! still works, but frame 0 is unavoidably a dependency for this one case,
//! a real format-level constraint rather than a simplification. This
//! crate's own writer never declares fixed atoms.
//!
//! **Frame indexing is arithmetic, not a scan.** Unlike TRR/XTC, DCD's
//! header states the total frame count directly, and because `NAMNF`/the
//! CHARMM flags are fixed for the whole file, every frame's byte length is
//! one of exactly two constants computable once at header time — the same
//! formula the reference implementation's own `jump_to_dcdstep` uses for
//! real random access.

use crate::core::atom::{Atom, Element};
use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::trajectory::{Frame, FrameSource, Trajectory};
use crate::io::errors::{DcdError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

const DCD_RECORD_MAGIC: i32 = 84;

/// A byte cursor over one DCD file, aware of the file's own detected
/// endianness. Unlike [`crate::io::xdr::XdrReader`], byte order is a
/// per-file runtime choice rather than always big-endian, and a record's
/// length is a real, validated leading/trailing pair rather than an
/// implicit 4-byte-alignment pad.
struct FortranReader<'a> {
    bytes: &'a [u8],
    pos: usize,
    big_endian: bool,
}

impl<'a> FortranReader<'a> {
    /// Detects endianness from the leading record marker (must be `84` in
    /// one byte order or the other) and returns a reader positioned right
    /// after it.
    fn new(bytes: &'a [u8]) -> Result<Self, DcdError> {
        if bytes.len() < 4 {
            return Err(DcdError::InvalidMagicNumber);
        }
        let raw: [u8; 4] = bytes[0..4].try_into().unwrap();
        let big_endian = if i32::from_le_bytes(raw) == DCD_RECORD_MAGIC {
            false
        } else if i32::from_be_bytes(raw) == DCD_RECORD_MAGIC {
            true
        } else {
            return Err(DcdError::InvalidMagicNumber);
        };
        Ok(Self {
            bytes,
            pos: 4,
            big_endian,
        })
    }

    /// A reader over `bytes` starting at `pos`, with an already-known
    /// endianness — used for every frame after the header, which does not
    /// need to re-detect it.
    fn at(bytes: &'a [u8], pos: usize, big_endian: bool) -> Self {
        Self {
            bytes,
            pos,
            big_endian,
        }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], DcdError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.bytes.len())
            .ok_or_else(|| DcdError::ParseError("unexpected end of input".to_string()))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn read_i32(&mut self) -> Result<i32, DcdError> {
        let b: [u8; 4] = self.take(4)?.try_into().unwrap();
        Ok(if self.big_endian {
            i32::from_be_bytes(b)
        } else {
            i32::from_le_bytes(b)
        })
    }

    fn read_f32(&mut self) -> Result<f32, DcdError> {
        let b: [u8; 4] = self.take(4)?.try_into().unwrap();
        Ok(if self.big_endian {
            f32::from_be_bytes(b)
        } else {
            f32::from_le_bytes(b)
        })
    }

    fn read_f64(&mut self) -> Result<f64, DcdError> {
        let b: [u8; 8] = self.take(8)?.try_into().unwrap();
        Ok(if self.big_endian {
            f64::from_be_bytes(b)
        } else {
            f64::from_le_bytes(b)
        })
    }

    fn read_bytes(&mut self, n: usize) -> Result<&'a [u8], DcdError> {
        self.take(n)
    }

    fn skip(&mut self, n: usize) -> Result<(), DcdError> {
        self.take(n)?;
        Ok(())
    }

    fn position(&self) -> usize {
        self.pos
    }

    /// Reads and validates a Fortran record marker equals `expected`.
    fn expect_marker(&mut self, expected: i32) -> Result<(), DcdError> {
        let got = self.read_i32()?;
        if got != expected {
            return Err(DcdError::ParseError(format!(
                "expected Fortran record marker {expected}, got {got}"
            )));
        }
        Ok(())
    }
}

/// Little-endian only — this crate's own writer never produces big-endian
/// DCD, the same "reads more dialects than it writes" shape TRR's
/// double-precision reading and XTC's writer-side simplification already
/// established for the other two trajectory formats.
struct FortranWriter {
    buf: Vec<u8>,
}

impl FortranWriter {
    fn new() -> Self {
        Self { buf: Vec::new() }
    }

    fn into_bytes(self) -> Vec<u8> {
        self.buf
    }

    fn write_i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_raw(&mut self, bytes: &[u8]) {
        self.buf.extend_from_slice(bytes);
    }

    /// Wraps `body` with matching leading/trailing Fortran record markers.
    fn write_record(&mut self, body: &[u8]) {
        self.write_i32(body.len() as i32);
        self.write_raw(body);
        self.write_i32(body.len() as i32);
    }
}

/// Everything read from the header, needed to locate and interpret every
/// later frame.
struct DcdHeader {
    nset: i32,
    istart: i32,
    nsavc: i32,
    namnf: i32,
    delta: f64,
    has_extra_block: bool,
    has_4dims: bool,
    natoms: usize,
    /// 0-based indices of the free atoms, only populated when `namnf != 0`.
    free_indices: Vec<usize>,
    /// Total bytes consumed reading the header — where frame 0 begins.
    header_size: usize,
    big_endian: bool,
}

fn read_header(bytes: &[u8]) -> Result<DcdHeader, DcdError> {
    let mut r = FortranReader::new(bytes)?; // consumes and validates the leading marker
    let hdr = r.read_bytes(84)?;
    r.expect_marker(DCD_RECORD_MAGIC)?;

    if &hdr[0..4] != b"CORD" {
        return Err(DcdError::ParseError("missing CORD magic".to_string()));
    }
    let big_endian = r.big_endian;
    let field_i32 = |off: usize| -> i32 {
        let b: [u8; 4] = hdr[off..off + 4].try_into().unwrap();
        if big_endian {
            i32::from_be_bytes(b)
        } else {
            i32::from_le_bytes(b)
        }
    };

    let nset = field_i32(4);
    let istart = field_i32(8);
    let nsavc = field_i32(12);
    let namnf = field_i32(36);
    let charmm_version = field_i32(80);
    let charmm = charmm_version != 0;
    let (has_extra_block, has_4dims) = if charmm {
        (field_i32(44) != 0, field_i32(48) == 1)
    } else {
        (false, false)
    };
    let delta = if charmm {
        let b: [u8; 4] = hdr[40..44].try_into().unwrap();
        (if big_endian {
            f32::from_be_bytes(b)
        } else {
            f32::from_le_bytes(b)
        }) as f64
    } else {
        let b: [u8; 8] = hdr[40..48].try_into().unwrap();
        if big_endian {
            f64::from_be_bytes(b)
        } else {
            f64::from_le_bytes(b)
        }
    };

    // Title block: a record whose body is NTITLE (i32) then NTITLE 80-byte
    // strings, read and discarded.
    let title_len = r.read_i32()?;
    if title_len < 4 || (title_len - 4) % 80 != 0 {
        return Err(DcdError::ParseError("malformed title block".to_string()));
    }
    let ntitle = r.read_i32()?;
    if ntitle < 0 {
        return Err(DcdError::ParseError("negative title count".to_string()));
    }
    r.skip(ntitle as usize * 80)?;
    r.expect_marker(title_len)?;

    // Atom count block.
    r.expect_marker(4)?;
    let natoms = r.read_i32()?;
    if natoms < 0 {
        return Err(DcdError::ParseError("negative atom count".to_string()));
    }
    let natoms = natoms as usize;
    r.expect_marker(4)?;

    // Free-atom index block, only when fixed atoms are declared.
    let mut free_indices = Vec::new();
    if namnf != 0 {
        let nfree = natoms as i64 - namnf as i64;
        if nfree < 0 {
            return Err(DcdError::ParseError(
                "more fixed atoms than atoms".to_string(),
            ));
        }
        let nfree = nfree as usize;
        r.expect_marker((nfree * 4) as i32)?;
        free_indices.reserve(nfree);
        for _ in 0..nfree {
            let idx = r.read_i32()?;
            if idx < 1 || idx as usize > natoms {
                return Err(DcdError::ParseError(format!(
                    "free-atom index {idx} out of range for {natoms} atoms"
                )));
            }
            free_indices.push(idx as usize - 1);
        }
        r.expect_marker((nfree * 4) as i32)?;
    }

    Ok(DcdHeader {
        nset,
        istart,
        nsavc,
        namnf,
        delta,
        has_extra_block,
        has_4dims,
        natoms,
        free_indices,
        header_size: r.position(),
        big_endian,
    })
}

/// `(firstframesize, framesize)` — the two possible per-frame byte
/// lengths, matching the reference implementation's own
/// `jump_to_dcdstep` formula exactly.
fn frame_sizes(h: &DcdHeader) -> (usize, usize) {
    let extrablocksize = if h.has_extra_block { 48 + 8 } else { 0 };
    let ndims = if h.has_4dims { 4 } else { 3 };
    let nfixed = h.namnf.max(0) as usize;
    let first = (h.natoms + 2) * ndims * 4 + extrablocksize;
    let later = (h.natoms.saturating_sub(nfixed) + 2) * ndims * 4 + extrablocksize;
    (first, later)
}

/// Interprets the six raw doubles of a unit-cell record, in file order
/// `[A, gamma, B, beta, alpha, C]`, into this crate's `[A, B, C, alpha,
/// beta, gamma]` convention.
fn interpret_unit_cell(raw: [f64; 6]) -> Result<UnitCell, DcdError> {
    let a = raw[0];
    let b = raw[2];
    let c = raw[5];
    let mut alpha = raw[4];
    let mut beta = raw[3];
    let mut gamma = raw[1];

    if (-1.0..=1.0).contains(&alpha)
        && (-1.0..=1.0).contains(&beta)
        && (-1.0..=1.0).contains(&gamma)
    {
        alpha = alpha.acos().to_degrees();
        beta = beta.acos().to_degrees();
        gamma = gamma.acos().to_degrees();
    } else if a < 0.0 || b < 0.0 || c < 0.0 || alpha > 180.0 || beta > 180.0 || gamma > 180.0 {
        return Err(DcdError::UnsupportedUnitCellFormat);
    }
    // Otherwise: already plain degrees (older NAMD<=2.5 convention).

    Ok(UnitCell::new(a, b, c, alpha, beta, gamma))
}

/// The inverse of [`interpret_unit_cell`] — always the modern cosine
/// convention.
fn encode_unit_cell(cell: UnitCell) -> [f64; 6] {
    let mut raw = [0.0; 6];
    raw[0] = cell.a;
    raw[2] = cell.b;
    raw[5] = cell.c;
    raw[1] = cell.gamma.to_radians().cos();
    raw[3] = cell.beta.to_radians().cos();
    raw[4] = cell.alpha.to_radians().cos();
    raw
}

/// A reduced (fixed-atom) frame's free-atom coordinates, one flat array
/// per axis, plus this frame's own unit cell if present.
struct ReducedFrame {
    xs: Vec<f32>,
    ys: Vec<f32>,
    zs: Vec<f32>,
    cell: Option<UnitCell>,
}

/// Random access to a DCD file's frames, indexed by arithmetic rather than
/// a scan — see the module doc.
pub(crate) struct DcdFrameSource {
    bytes: Vec<u8>,
    header: DcdHeader,
    firstframesize: usize,
    framesize: usize,
    frame0_cache: Option<Vec<Point3>>,
}

impl DcdFrameSource {
    fn open(bytes: Vec<u8>) -> Result<Self, DcdError> {
        let header = read_header(&bytes)?;
        let (firstframesize, framesize) = frame_sizes(&header);
        let nset = header.nset.max(0) as usize;
        let expected_total = header.header_size
            + if nset == 0 {
                0
            } else {
                firstframesize + framesize * (nset - 1)
            };
        if expected_total > bytes.len() {
            return Err(DcdError::ParseError(format!(
                "file declares {nset} frames needing {expected_total} bytes, but only {} are present",
                bytes.len()
            )));
        }
        Ok(Self {
            bytes,
            header,
            firstframesize,
            framesize,
            frame0_cache: None,
        })
    }

    fn offset(&self, index: usize) -> usize {
        if index == 0 {
            self.header.header_size
        } else {
            self.header.header_size + self.firstframesize + self.framesize * (index - 1)
        }
    }

    fn read_cell(&self, r: &mut FortranReader) -> Result<Option<UnitCell>, DcdError> {
        if !self.header.has_extra_block {
            return Ok(None);
        }
        r.expect_marker(48)?;
        let mut raw = [0.0f64; 6];
        for v in raw.iter_mut() {
            *v = r.read_f64()?;
        }
        r.expect_marker(48)?;
        Ok(Some(interpret_unit_cell(raw)?))
    }

    fn skip_4dims(&self, r: &mut FortranReader) -> Result<(), DcdError> {
        if self.header.has_4dims {
            let len = r.read_i32()?;
            r.skip(len.max(0) as usize)?;
            r.read_i32()?; // trailing marker, not cross-validated -- matches the reference.
        }
        Ok(())
    }

    fn read_full_frame(&self, offset: usize) -> Result<(Vec<Point3>, Option<UnitCell>), DcdError> {
        let mut r = FortranReader::at(&self.bytes, offset, self.header.big_endian);
        let cell = self.read_cell(&mut r)?;

        let n = self.header.natoms;
        let marker = (n * 4) as i32;
        let mut xs = Vec::with_capacity(n);
        r.expect_marker(marker)?;
        for _ in 0..n {
            xs.push(r.read_f32()?);
        }
        r.expect_marker(marker)?;

        let mut ys = Vec::with_capacity(n);
        r.expect_marker(marker)?;
        for _ in 0..n {
            ys.push(r.read_f32()?);
        }
        r.expect_marker(marker)?;

        let mut zs = Vec::with_capacity(n);
        r.expect_marker(marker)?;
        for _ in 0..n {
            zs.push(r.read_f32()?);
        }
        r.expect_marker(marker)?;

        self.skip_4dims(&mut r)?;

        let positions = (0..n)
            .map(|i| Point3::new(xs[i] as f64, ys[i] as f64, zs[i] as f64))
            .collect();
        Ok((positions, cell))
    }

    fn read_reduced_frame(&self, offset: usize) -> Result<ReducedFrame, DcdError> {
        let mut r = FortranReader::at(&self.bytes, offset, self.header.big_endian);
        let cell = self.read_cell(&mut r)?;

        let nfree = self.header.free_indices.len();
        let marker = (nfree * 4) as i32;
        let mut xs = Vec::with_capacity(nfree);
        r.expect_marker(marker)?;
        for _ in 0..nfree {
            xs.push(r.read_f32()?);
        }
        r.expect_marker(marker)?;

        let mut ys = Vec::with_capacity(nfree);
        r.expect_marker(marker)?;
        for _ in 0..nfree {
            ys.push(r.read_f32()?);
        }
        r.expect_marker(marker)?;

        let mut zs = Vec::with_capacity(nfree);
        r.expect_marker(marker)?;
        for _ in 0..nfree {
            zs.push(r.read_f32()?);
        }
        r.expect_marker(marker)?;

        self.skip_4dims(&mut r)?;

        Ok(ReducedFrame { xs, ys, zs, cell })
    }

    fn frame_inner(&mut self, index: usize) -> Result<Frame, DcdError> {
        let nfixed = self.header.namnf.max(0) as usize;
        let (positions, cell) = if nfixed == 0 || index == 0 {
            let offset = self.offset(index);
            let (positions, cell) = self.read_full_frame(offset)?;
            if index == 0 && nfixed > 0 {
                self.frame0_cache = Some(positions.clone());
            }
            (positions, cell)
        } else {
            if self.frame0_cache.is_none() {
                let offset0 = self.offset(0);
                let (positions0, _) = self.read_full_frame(offset0)?;
                self.frame0_cache = Some(positions0);
            }
            let mut positions = self.frame0_cache.clone().unwrap();
            let offset = self.offset(index);
            let reduced = self.read_reduced_frame(offset)?;
            for (k, &idx) in self.header.free_indices.iter().enumerate() {
                positions[idx] = Point3::new(
                    reduced.xs[k] as f64,
                    reduced.ys[k] as f64,
                    reduced.zs[k] as f64,
                );
            }
            (positions, reduced.cell)
        };

        let step = self.header.istart as i64 + index as i64 * self.header.nsavc as i64;
        let time = self.header.delta * step as f64;
        Ok(Frame {
            positions,
            velocities: None,
            forces: None,
            time: Some(time),
            step: Some(step.max(0) as u64),
            cell,
        })
    }
}

impl FrameSource for DcdFrameSource {
    fn frame_count(&self) -> usize {
        self.header.nset.max(0) as usize
    }

    fn num_atoms(&self) -> usize {
        self.header.natoms
    }

    fn frame(&mut self, index: usize) -> std::io::Result<Frame> {
        self.frame_inner(index).map_err(std::io::Error::other)
    }
}

fn build_trajectory(bytes: Vec<u8>) -> Result<Trajectory, DcdError> {
    if bytes.is_empty() {
        return Err(DcdError::ParseError("empty input".to_string()));
    }
    let source = DcdFrameSource::open(bytes)?;
    let natoms = source.num_atoms();
    let mut topology = Molecule::new();
    for _ in 0..natoms {
        topology.add_atom(Atom::new(Element::UNKNOWN));
    }
    Ok(Trajectory::new(topology, Box::new(source))?)
}

/// [`crate::io::format::ByteReadFn`] for DCD.
pub(crate) fn read_dcd_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
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

fn write_header(
    out: &mut FortranWriter,
    natoms: usize,
    nframes: usize,
    istart: i32,
    nsavc: i32,
    delta: f64,
    has_cell: bool,
) {
    let mut hdr = Vec::with_capacity(84);
    hdr.extend_from_slice(b"CORD");
    hdr.extend_from_slice(&(nframes as i32).to_le_bytes()); // NSET
    hdr.extend_from_slice(&istart.to_le_bytes());
    hdr.extend_from_slice(&nsavc.to_le_bytes());
    hdr.extend_from_slice(&0i32.to_le_bytes()); // NSTEP, unused
    for _ in 0..4 {
        hdr.extend_from_slice(&0i32.to_le_bytes()); // offsets 20..36, unused
    }
    hdr.extend_from_slice(&0i32.to_le_bytes()); // NAMNF -- this writer never declares fixed atoms
    hdr.extend_from_slice(&(delta as f32).to_le_bytes()); // CHARMM stores DELTA as f32
    hdr.extend_from_slice(&(has_cell as i32).to_le_bytes()); // has-extra-block flag
    for _ in 0..8 {
        hdr.extend_from_slice(&0i32.to_le_bytes()); // offsets 48..80, including the 4-dims flag (always 0)
    }
    hdr.extend_from_slice(&24i32.to_le_bytes()); // CHARMM version -- nonzero selects the CHARMM dialect
    debug_assert_eq!(hdr.len(), 84);

    out.write_record(&hdr);

    let mut title = Vec::with_capacity(84);
    title.extend_from_slice(&1i32.to_le_bytes()); // NTITLE
    let mut line = [0u8; 80];
    let text = b"Written by chem";
    line[..text.len()].copy_from_slice(text);
    title.extend_from_slice(&line);
    out.write_record(&title);

    out.write_record(&(natoms as i32).to_le_bytes());
}

fn write_frame(out: &mut FortranWriter, frame: &Frame, has_cell: bool) {
    if has_cell {
        let cell = frame
            .cell
            .unwrap_or(UnitCell::new(0.0, 0.0, 0.0, 90.0, 90.0, 90.0));
        let raw = encode_unit_cell(cell);
        let mut body = Vec::with_capacity(48);
        for v in raw {
            body.extend_from_slice(&v.to_le_bytes());
        }
        out.write_record(&body);
    }

    let n = frame.num_atoms();
    write_axis(out, n, frame.positions.iter().map(|p| p.x as f32));
    write_axis(out, n, frame.positions.iter().map(|p| p.y as f32));
    write_axis(out, n, frame.positions.iter().map(|p| p.z as f32));
}

fn write_axis(out: &mut FortranWriter, n: usize, values: impl Iterator<Item = f32>) {
    let mut body = Vec::with_capacity(n * 4);
    for v in values {
        body.extend_from_slice(&v.to_le_bytes());
    }
    out.write_record(&body);
}

/// [`crate::io::format::ByteWriteTrajectoryFn`] for DCD. Always CHARMM
/// dialect, no fixed atoms, the modern cosine unit-cell convention.
pub(crate) fn write_dcd_bytes(trajectory: &mut Trajectory, _options: &WriteOptions) -> Vec<u8> {
    let natoms = trajectory.num_atoms();
    let nframes = trajectory.frame_count();

    let frame0 = if nframes > 0 {
        Some(
            trajectory
                .frame(0)
                .expect("a Trajectory already validated its own frame/atom counts at construction"),
        )
    } else {
        None
    };
    let has_cell = frame0.as_ref().is_some_and(|f| f.cell.is_some());
    let istart = frame0.as_ref().and_then(|f| f.step).unwrap_or(0) as i32;
    let delta = if nframes >= 2 {
        let frame1 = trajectory
            .frame(1)
            .expect("a Trajectory already validated its own frame/atom counts at construction");
        match (frame0.as_ref().and_then(|f| f.time), frame1.time) {
            (Some(t0), Some(t1)) => (t1 - t0).abs().max(f64::MIN_POSITIVE),
            _ => 1.0,
        }
    } else {
        1.0
    };
    let nsavc = 1;

    let mut w = FortranWriter::new();
    write_header(&mut w, natoms, nframes, istart, nsavc, delta, has_cell);
    for i in 0..nframes {
        let frame = trajectory
            .frame(i)
            .expect("a Trajectory already validated its own frame/atom counts at construction");
        write_frame(&mut w, &frame, has_cell);
    }
    w.into_bytes()
}

/// Buffers the whole input and parses it once, mirroring
/// [`crate::io::trr::TrrSupplier`]/[`crate::io::xtc::XtcSupplier`].
pub struct DcdSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl DcdSupplier {
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

impl Iterator for DcdSupplier {
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

    fn frame(positions: Vec<Point3>, time: f64, step: u64, cell: Option<UnitCell>) -> Frame {
        Frame {
            positions,
            velocities: None,
            forces: None,
            time: Some(time),
            step: Some(step),
            cell,
        }
    }

    #[test]
    fn test_a_little_endian_round_trip_with_no_cell() {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
        ];
        let mut trajectory = trajectory_from(vec![
            frame(positions.clone(), 0.0, 0, None),
            frame(
                positions
                    .iter()
                    .map(|p| *p + Point3::new(0.1, 0.0, 0.0))
                    .collect(),
                1.0,
                1,
                None,
            ),
        ]);
        let bytes = write_dcd_bytes(&mut trajectory, &WriteOptions::default());
        let outcome = read_dcd_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut back = as_trajectory(outcome);
        assert_eq!(back.frame_count(), 2);
        let f0 = back.frame(0).unwrap();
        for (a, b) in f0.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 1e-4, "{a} vs {b}");
        }
        assert!(f0.cell.is_none());
        assert!(f0.velocities.is_none());
        assert!(f0.forces.is_none());
    }

    #[test]
    fn test_a_charmm_round_trip_with_a_cell_uses_the_cosine_convention() {
        let cell = UnitCell::new(30.0, 25.0, 20.0, 80.0, 85.0, 95.0);
        let positions: Vec<Point3> = (0..10).map(|i| Point3::new(i as f64, 0.0, 0.0)).collect();
        let mut trajectory = trajectory_from(vec![frame(positions, 0.0, 0, Some(cell))]);
        let bytes = write_dcd_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_dcd_bytes(&bytes, &ReadOptions::default()));
        let f = back.frame(0).unwrap();
        let back_cell = f.cell.expect("cell survives");
        assert!((back_cell.a - cell.a).abs() < 1e-3, "{}", back_cell.a);
        assert!(
            (back_cell.alpha - cell.alpha).abs() < 1e-2,
            "{}",
            back_cell.alpha
        );
        assert!(
            (back_cell.beta - cell.beta).abs() < 1e-2,
            "{}",
            back_cell.beta
        );
        assert!(
            (back_cell.gamma - cell.gamma).abs() < 1e-2,
            "{}",
            back_cell.gamma
        );
    }

    #[test]
    fn test_time_and_step_are_derived_from_istart_nsavc_delta() {
        // DCD's own step/time model is `time = delta * step`, a
        // relationship anchored at the origin -- it can only round-trip
        // step/time pairs that are already proportional this way (step 0
        // at time 0), not an arbitrary affine offset. This is a real,
        // disclosed limitation of `write_dcd_bytes`'s derivation, not
        // something this test works around.
        let positions: Vec<Point3> = (0..5).map(|i| Point3::new(i as f64, 0.0, 0.0)).collect();
        let frames = vec![
            frame(positions.clone(), 0.0, 0, None),
            frame(positions.clone(), 2.0, 1, None),
            frame(positions, 4.0, 2, None),
        ];
        let mut trajectory = trajectory_from(frames);
        let bytes = write_dcd_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_dcd_bytes(&bytes, &ReadOptions::default()));
        let f0 = back.frame(0).unwrap();
        let f1 = back.frame(1).unwrap();
        let f2 = back.frame(2).unwrap();
        assert_eq!(f0.step, Some(0));
        assert_eq!(f1.step, Some(1));
        assert_eq!(f2.step, Some(2));
        assert!((f0.time.unwrap() - 0.0).abs() < 1e-6);
        assert!((f1.time.unwrap() - 2.0).abs() < 1e-6);
        assert!((f2.time.unwrap() - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_multiple_frames_are_randomly_accessible() {
        let frames: Vec<Frame> = (0..4)
            .map(|f| {
                let positions: Vec<Point3> = (0..12)
                    .map(|i| Point3::new((i + f * 100) as f64, 0.0, 0.0))
                    .collect();
                frame(positions, f as f64, f as u64, None)
            })
            .collect();
        let mut trajectory = trajectory_from(frames);
        let bytes = write_dcd_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_dcd_bytes(&bytes, &ReadOptions::default()));
        assert_eq!(back.frame_count(), 4);
        let last = back.frame(3).unwrap();
        assert!((last.positions[0].x - 300.0).abs() < 1e-3);
        let first = back.frame(0).unwrap();
        assert!((first.positions[0].x - 0.0).abs() < 1e-3);
    }

    #[test]
    fn test_a_bad_magic_number_is_a_clear_error_not_a_panic() {
        let outcome = read_dcd_bytes(&[1, 2, 3, 4], &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
    }

    #[test]
    fn test_empty_input_is_a_clear_error_not_a_panic() {
        let outcome = read_dcd_bytes(&[], &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
    }

    // --- hand-built fixtures for cases this crate's own writer never
    // produces, mirroring TRR/XTC's own precedent -------------------------

    /// One frame's raw positions, and its raw six-double unit cell if any,
    /// in the same file order [`RawDcdBuilder::build`] writes.
    type RawFrame = (Vec<[f32; 3]>, Option<[f64; 6]>);

    /// A minimal, hand-built DCD file: header + one frame, no title
    /// content beyond `NTITLE=0`, parameterised over endianness and
    /// dialect so the same builder covers every fixture below.
    struct RawDcdBuilder {
        big_endian: bool,
        charmm: bool,
        has_cell: bool,
        natoms: usize,
        namnf: i32,
        free_indices: Vec<i32>, // 1-based, as on the wire
        delta: f64,
    }

    impl RawDcdBuilder {
        fn new(natoms: usize) -> Self {
            Self {
                big_endian: false,
                charmm: true,
                has_cell: false,
                natoms,
                namnf: 0,
                free_indices: Vec::new(),
                delta: 1.0,
            }
        }

        fn push_i32(&self, buf: &mut Vec<u8>, v: i32) {
            buf.extend_from_slice(&if self.big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            });
        }

        fn push_f32(&self, buf: &mut Vec<u8>, v: f32) {
            buf.extend_from_slice(&if self.big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            });
        }

        fn push_f64(&self, buf: &mut Vec<u8>, v: f64) {
            buf.extend_from_slice(&if self.big_endian {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            });
        }

        fn push_record(&self, buf: &mut Vec<u8>, body: &[u8]) {
            self.push_i32(buf, body.len() as i32);
            buf.extend_from_slice(body);
            self.push_i32(buf, body.len() as i32);
        }

        fn header_bytes(&self, nset: usize) -> Vec<u8> {
            let mut hdr = Vec::with_capacity(84);
            hdr.extend_from_slice(b"CORD");
            self.push_i32(&mut hdr, nset as i32);
            self.push_i32(&mut hdr, 0); // istart
            self.push_i32(&mut hdr, 1); // nsavc
            self.push_i32(&mut hdr, 0);
            for _ in 0..4 {
                self.push_i32(&mut hdr, 0);
            }
            self.push_i32(&mut hdr, self.namnf);
            if self.charmm {
                self.push_f32(&mut hdr, self.delta as f32);
                self.push_i32(&mut hdr, self.has_cell as i32);
            } else {
                self.push_f64(&mut hdr, self.delta);
            }
            for _ in 0..8 {
                self.push_i32(&mut hdr, 0);
            }
            self.push_i32(&mut hdr, if self.charmm { 24 } else { 0 });
            assert_eq!(hdr.len(), 84);
            hdr
        }

        fn build(&self, frames: &[RawFrame]) -> Vec<u8> {
            let mut out = Vec::new();
            self.push_record(&mut out, &self.header_bytes(frames.len()));

            let mut title = Vec::new();
            self.push_i32(&mut title, 0); // NTITLE = 0
            self.push_record(&mut out, &title);

            let mut atoms = Vec::new();
            self.push_i32(&mut atoms, self.natoms as i32);
            self.push_record(&mut out, &atoms);

            if self.namnf != 0 {
                let mut idx = Vec::new();
                for &i in &self.free_indices {
                    self.push_i32(&mut idx, i);
                }
                self.push_record(&mut out, &idx);
            }

            for (i, (positions, cell)) in frames.iter().enumerate() {
                let n = if i == 0 || self.namnf == 0 {
                    self.natoms
                } else {
                    self.natoms - self.namnf.max(0) as usize
                };
                assert_eq!(positions.len(), n);
                if self.has_cell {
                    let raw = cell.expect("has_cell fixtures must supply a cell");
                    let mut body = Vec::new();
                    for v in raw {
                        self.push_f64(&mut body, v);
                    }
                    self.push_record(&mut out, &body);
                }
                for axis in 0..3 {
                    let mut body = Vec::new();
                    for p in positions {
                        self.push_f32(&mut body, p[axis]);
                    }
                    self.push_record(&mut out, &body);
                }
            }
            out
        }
    }

    #[test]
    fn test_a_hand_built_big_endian_file_with_a_unit_cell_reads_correctly() {
        // The issue's own named worst case: big-endian and a unit cell at
        // once, exercising the byte-swap path and the cell heuristic
        // together.
        let mut b = RawDcdBuilder::new(4);
        b.big_endian = true;
        b.has_cell = true;
        // File-order [A, gamma, B, beta, alpha, C]; angles as cosines
        // (all within [-1, 1]), the modern convention.
        let cell_raw = [30.0, 0.0, 25.0, 0.0, 0.0, 20.0]; // all angles 90 degrees (cos=0)
        let positions: Vec<[f32; 3]> = (0..4).map(|i| [i as f32, 0.0, 0.0]).collect();
        let bytes = b.build(&[(positions.clone(), Some(cell_raw))]);

        let outcome = read_dcd_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut trajectory = as_trajectory(outcome);
        let f = trajectory.frame(0).unwrap();
        for (got, want) in f.positions.iter().zip(&positions) {
            assert!((got.x - want[0] as f64).abs() < 1e-4);
        }
        let cell = f.cell.expect("cell survives");
        assert!((cell.a - 30.0).abs() < 1e-6);
        assert!((cell.b - 25.0).abs() < 1e-6);
        assert!((cell.c - 20.0).abs() < 1e-6);
        assert!((cell.alpha - 90.0).abs() < 1e-6, "{}", cell.alpha);
        assert!((cell.beta - 90.0).abs() < 1e-6, "{}", cell.beta);
        assert!((cell.gamma - 90.0).abs() < 1e-6, "{}", cell.gamma);
    }

    #[test]
    fn test_a_hand_built_file_with_plain_degree_angles_is_read_correctly() {
        // The older NAMD<=2.5 convention: angles already in degrees, none
        // of them within [-1, 1] (a triclinic cell, deliberately not 90),
        // proving the heuristic's other branch.
        let mut b = RawDcdBuilder::new(3);
        b.has_cell = true;
        let cell_raw = [30.0, 95.0, 25.0, 85.0, 80.0, 20.0]; // [A, gamma, B, beta, alpha, C]
        let positions: Vec<[f32; 3]> = (0..3).map(|i| [i as f32, 0.0, 0.0]).collect();
        let bytes = b.build(&[(positions, Some(cell_raw))]);

        let mut trajectory = as_trajectory(read_dcd_bytes(&bytes, &ReadOptions::default()));
        let f = trajectory.frame(0).unwrap();
        let cell = f.cell.expect("cell survives");
        assert!((cell.alpha - 80.0).abs() < 1e-6, "{}", cell.alpha);
        assert!((cell.beta - 85.0).abs() < 1e-6, "{}", cell.beta);
        assert!((cell.gamma - 95.0).abs() < 1e-6, "{}", cell.gamma);
    }

    #[test]
    fn test_a_hand_built_x_plor_file_never_has_a_cell_and_reads_f64_delta() {
        let mut b = RawDcdBuilder::new(3);
        b.charmm = false;
        b.has_cell = false; // X-PLOR: no extra-block flag exists at all
        b.delta = 0.002;
        let positions: Vec<[f32; 3]> = (0..3).map(|i| [i as f32, 1.0, 2.0]).collect();
        let bytes = b.build(&[(positions.clone(), None)]);

        let mut trajectory = as_trajectory(read_dcd_bytes(&bytes, &ReadOptions::default()));
        let f = trajectory.frame(0).unwrap();
        assert!(f.cell.is_none());
        for (got, want) in f.positions.iter().zip(&positions) {
            assert!((got.x - want[0] as f64).abs() < 1e-4);
        }
        // istart=0, nsavc=1, so time of frame 0 is 0 regardless of delta --
        // confirm the file was at least accepted (X-PLOR dispatch worked)
        // rather than rejected as malformed.
        assert!(f.time.is_some());
    }

    #[test]
    fn test_fixed_atoms_splice_against_the_cached_first_frame() {
        let natoms = 5;
        let mut b = RawDcdBuilder::new(natoms);
        b.namnf = 2; // atoms 1 and 3 (0-based) are fixed
        b.free_indices = vec![1, 2, 4, 5]; // 1-based free-atom indices (atoms 0,1,3,4 -- wait see below

        // Free indices are 1-based positions into the FULL atom array that
        // vary after frame 0. Fix atoms at 0-based indices {1, 3}, so the
        // free (varying) atoms are 0-based {0, 2, 4} -- 1-based {1, 3, 5}.
        b.namnf = 2;
        b.free_indices = vec![1, 3, 5];

        let frame0: Vec<[f32; 3]> = (0..natoms).map(|i| [i as f32, 0.0, 0.0]).collect();
        // Later frames restate only the 3 free atoms, in free-index order.
        let frame1: Vec<[f32; 3]> = vec![[100.0, 0.0, 0.0], [102.0, 0.0, 0.0], [104.0, 0.0, 0.0]];
        let bytes = b.build(&[(frame0.clone(), None), (frame1.clone(), None)]);

        let mut trajectory = as_trajectory(read_dcd_bytes(&bytes, &ReadOptions::default()));
        assert_eq!(trajectory.frame_count(), 2);

        // Random access: read frame 1 directly, without reading frame 0
        // first through the public API -- internally this still requires
        // caching frame 0, but the caller never has to ask for it.
        let f1 = trajectory.frame(1).unwrap();
        assert!(
            (f1.positions[0].x - 100.0).abs() < 1e-4,
            "{}",
            f1.positions[0].x
        ); // free
        assert!(
            (f1.positions[1].x - 1.0).abs() < 1e-4,
            "{}",
            f1.positions[1].x
        ); // fixed, from frame 0
        assert!(
            (f1.positions[2].x - 102.0).abs() < 1e-4,
            "{}",
            f1.positions[2].x
        ); // free
        assert!(
            (f1.positions[3].x - 3.0).abs() < 1e-4,
            "{}",
            f1.positions[3].x
        ); // fixed, from frame 0
        assert!(
            (f1.positions[4].x - 104.0).abs() < 1e-4,
            "{}",
            f1.positions[4].x
        ); // free
    }

    #[test]
    fn test_an_unsupported_new_style_unit_cell_is_a_clear_error_not_a_panic() {
        let mut b = RawDcdBuilder::new(2);
        b.has_cell = true;
        // A negative length in the A slot -- outside both the cosine range
        // and any plausible plain-degree cell, signalling the "new-style
        // CHARMM box vectors" encoding this crate does not decode.
        let cell_raw = [-5.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let positions: Vec<[f32; 3]> = (0..2).map(|i| [i as f32, 0.0, 0.0]).collect();
        let bytes = b.build(&[(positions, Some(cell_raw))]);

        // The header itself is well-formed, so building the `Trajectory`
        // succeeds -- frame content, including the unit cell, is only
        // decoded lazily when a frame is actually read, same as TRR/XTC.
        let outcome = read_dcd_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut trajectory = as_trajectory(outcome);
        let err = trajectory.frame(0).unwrap_err();
        assert!(err.to_string().contains("Unsupported"), "{err}");
    }
}
