//! GROMACS TRR — an XDR-encoded (#326) sequence of self-describing frames,
//! each stating which of position/velocity/force/box arrays it carries.
//!
//! The classic `xdrfile` header layout is the one implemented against (not
//! the newer GROMACS-internal one with a trailing `fep_state`), since it is
//! the one every third-party tool (MDAnalysis, mdtraj) still replicates:
//! magic number (must be `1993`), a version string, a run of legacy
//! "backward compatibility" size fields (`ir_size`/`e_size`/`top_size`/
//! `sym_size`/`nre` — read and discarded, no implementation old or new ever
//! stores data for them), `box_size`/`vir_size`/`pres_size` (the last two
//! also legacy, but still gate a real skip when non-zero), `x_size`/
//! `v_size`/`f_size`, `natoms`, `step`, then `t`/`lambda` in whichever
//! precision the frame turns out to use.
//!
//! **Precision is detected per frame**, exactly like the reference
//! implementation: `box_size` decides it when a box is present (`36` bytes
//! = single/`f32`, `72` = double/`f64`), falling back to `x_size`/`v_size`/
//! `f_size` in that order otherwise. This reader accepts both; the writer
//! here only ever produces single precision, GROMACS's own overwhelmingly
//! common default — the same "reads more dialects than it writes" shape
//! [`crate::io::options::SdfWriteOptions`]'s `MolfileVersion` already has
//! for V3000.
//!
//! **Every atom is [`Element::UNKNOWN`]**: TRR states no chemical identity
//! at all, not even a numeric type — only a bare atom count.
//!
//! **Units**: nm (length/box) and nm/ps (velocity) both convert to Å by
//! multiplying by this crate's shared `NM_TO_ANGSTROM` constant. **Force is
//! the one place this inverts**: force is energy-per-length (kJ/mol/nm), so
//! its numeric value converts by *dividing* by the same constant on read,
//! and multiplying on write — the exact inverse of positions/velocities,
//! and an easy direction to get backwards.
//!
//! **No native random-access index exists in a TRR file.** This module's
//! own frame source builds one on construction with a single header-only
//! pass: each frame's header states its own array sizes, so the body can be
//! skipped without decoding it, matching [`FrameSource`]'s own documented
//! contract and the same approach MDAnalysis's `XDR.py` uses.

use std::io::{self, Read, Seek, SeekFrom};

use crate::core::atom::{Atom, Element};
use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::trajectory::{Frame, FrameSource, Trajectory};
use crate::core::units::NM_TO_ANGSTROM;
use crate::io::errors::{ReadError, TrrError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};
use crate::io::xdr::{XdrReader, XdrWriter};

const GROMACS_MAGIC: i32 = 1993;
const TRR_VERSION_STRING: &str = "GMX_trn_file";

struct FrameHeader {
    box_size: i32,
    vir_size: i32,
    pres_size: i32,
    x_size: i32,
    v_size: i32,
    f_size: i32,
    natoms: usize,
    step: i32,
    t: f64,
    floatsize: usize,
}

fn read_header(r: &mut XdrReader) -> Result<FrameHeader, TrrError> {
    let magic = r.read_i32()?;
    if magic != GROMACS_MAGIC {
        return Err(TrrError::InvalidMagicNumber(magic));
    }
    // The real TRR header writes a redundant outer length -- `strlen(s) +
    // 1`, as a plain integer -- immediately before the version string
    // itself, which then carries its own, separate XDR string length
    // (`strlen(s)`, no NUL). Confirmed against the classic `xdrfile`
    // library's own `do_trnheader` (`xdrfile_trr.c`), which reads/writes
    // this exact extra field before ever calling `xdrfile_read_string`/
    // `xdrfile_write_string`.
    let _outer_slen = r.read_i32()?;
    let _version = r.read_string()?;
    let _ir_size = r.read_i32()?;
    let _e_size = r.read_i32()?;
    let box_size = r.read_i32()?;
    let vir_size = r.read_i32()?;
    let pres_size = r.read_i32()?;
    let _top_size = r.read_i32()?;
    let _sym_size = r.read_i32()?;
    let x_size = r.read_i32()?;
    let v_size = r.read_i32()?;
    let f_size = r.read_i32()?;
    let natoms = r.read_i32()?.max(0) as usize;
    let step = r.read_i32()?;
    let _nre = r.read_i32()?;

    let per_axis = (natoms * 3).max(1);
    let floatsize = if box_size != 0 {
        box_size as usize / 9
    } else if x_size != 0 {
        x_size as usize / per_axis
    } else if v_size != 0 {
        v_size as usize / per_axis
    } else if f_size != 0 {
        f_size as usize / per_axis
    } else {
        4
    };
    if floatsize != 4 && floatsize != 8 {
        return Err(TrrError::UnsupportedPrecision);
    }

    let t = if floatsize == 4 {
        let t = r.read_f32()? as f64;
        let _lambda = r.read_f32()?;
        t
    } else {
        let t = r.read_f64()?;
        let _lambda = r.read_f64()?;
        t
    };

    Ok(FrameHeader {
        box_size,
        vir_size,
        pres_size,
        x_size,
        v_size,
        f_size,
        natoms,
        step,
        t,
        floatsize,
    })
}

fn read_reals(r: &mut XdrReader, n: usize, floatsize: usize) -> Result<Vec<f64>, TrrError> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push(if floatsize == 4 {
            r.read_f32()? as f64
        } else {
            r.read_f64()?
        });
    }
    Ok(out)
}

fn read_point(r: &mut XdrReader, floatsize: usize) -> Result<Point3, TrrError> {
    let v = read_reals(r, 3, floatsize)?;
    Ok(Point3::new(v[0], v[1], v[2]))
}

fn read_points(
    r: &mut XdrReader,
    natoms: usize,
    floatsize: usize,
) -> Result<Vec<Point3>, TrrError> {
    let v = read_reals(r, natoms * 3, floatsize)?;
    Ok(v.chunks_exact(3)
        .map(|c| Point3::new(c[0], c[1], c[2]))
        .collect())
}

fn angle_degrees(a: Point3, b: Point3) -> f64 {
    (a.dot(b) / (a.length() * b.length()))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

/// GROMACS stores its box as three literal row vectors, unlike LAMMPS's
/// `xlo`/`xhi`-style bounds — recovering `a`/`b`/`c`/`alpha`/`beta`/`gamma`
/// is direct vector length/angle math, no closed-form conversion needed.
fn cell_from_box_vectors(v1: Point3, v2: Point3, v3: Point3) -> UnitCell {
    UnitCell::new(
        v1.length() * NM_TO_ANGSTROM,
        v2.length() * NM_TO_ANGSTROM,
        v3.length() * NM_TO_ANGSTROM,
        angle_degrees(v2, v3),
        angle_degrees(v1, v3),
        angle_degrees(v1, v2),
    )
}

/// Parses one whole frame (header and body) from `bytes`, which must start
/// exactly at the frame's own magic number.
fn parse_frame(bytes: &[u8]) -> Result<Frame, TrrError> {
    let mut r = XdrReader::new(bytes);
    let h = read_header(&mut r)?;

    let cell = if h.box_size != 0 {
        let v1 = read_point(&mut r, h.floatsize)?;
        let v2 = read_point(&mut r, h.floatsize)?;
        let v3 = read_point(&mut r, h.floatsize)?;
        Some(cell_from_box_vectors(v1, v2, v3))
    } else {
        None
    };
    if h.vir_size != 0 {
        r.skip(h.vir_size.max(0) as usize)?;
    }
    if h.pres_size != 0 {
        r.skip(h.pres_size.max(0) as usize)?;
    }

    let positions = if h.x_size != 0 {
        read_points(&mut r, h.natoms, h.floatsize)?
            .into_iter()
            .map(|p| p * NM_TO_ANGSTROM)
            .collect()
    } else {
        return Err(TrrError::ParseError("frame has no positions".to_string()));
    };
    let velocities = if h.v_size != 0 {
        Some(
            read_points(&mut r, h.natoms, h.floatsize)?
                .into_iter()
                .map(|p| p * NM_TO_ANGSTROM)
                .collect(),
        )
    } else {
        None
    };
    let forces = if h.f_size != 0 {
        Some(
            read_points(&mut r, h.natoms, h.floatsize)?
                .into_iter()
                .map(|p| p / NM_TO_ANGSTROM)
                .collect(),
        )
    } else {
        None
    };

    Ok(Frame {
        positions,
        velocities,
        forces,
        time: Some(h.t),
        step: Some(h.step.max(0) as u64),
        cell,
    })
}

/// One pass over the whole file building a byte-offset per frame — each
/// frame's header states its own array sizes, so the body is skipped rather
/// than decoded.
fn scan_offsets(bytes: &[u8]) -> Result<(Vec<u64>, usize), TrrError> {
    let mut offsets = Vec::new();
    let mut pos = 0usize;
    let mut natoms = None;
    while pos < bytes.len() {
        let mut r = XdrReader::new(&bytes[pos..]);
        let h = read_header(&mut r)?;
        offsets.push(pos as u64);
        if natoms.is_none() {
            natoms = Some(h.natoms);
        }
        let body_len = h.box_size.max(0) as usize
            + h.vir_size.max(0) as usize
            + h.pres_size.max(0) as usize
            + h.x_size.max(0) as usize
            + h.v_size.max(0) as usize
            + h.f_size.max(0) as usize;
        pos += r.position() + body_len;
    }
    Ok((offsets, natoms.unwrap_or(0)))
}

/// Random access to a TRR file's frames, indexed once on construction.
pub(crate) struct TrrFrameSource {
    cursor: io::Cursor<Vec<u8>>,
    offsets: Vec<u64>,
    natoms: usize,
}

impl TrrFrameSource {
    fn open(bytes: Vec<u8>) -> Result<Self, TrrError> {
        let (offsets, natoms) = scan_offsets(&bytes)?;
        Ok(Self {
            cursor: io::Cursor::new(bytes),
            offsets,
            natoms,
        })
    }
}

impl FrameSource for TrrFrameSource {
    fn frame_count(&self) -> usize {
        self.offsets.len()
    }

    fn num_atoms(&self) -> usize {
        self.natoms
    }

    fn frame(&mut self, index: usize) -> io::Result<Frame> {
        self.cursor.seek(SeekFrom::Start(self.offsets[index]))?;
        let mut buf = Vec::new();
        self.cursor.read_to_end(&mut buf)?;
        parse_frame(&buf).map_err(io::Error::other)
    }
}

fn build_trajectory(bytes: Vec<u8>) -> Result<Trajectory, TrrError> {
    if bytes.is_empty() {
        return Err(TrrError::ParseError("empty input".to_string()));
    }
    let source = TrrFrameSource::open(bytes)?;
    let natoms = source.num_atoms();
    let mut topology = Molecule::new();
    for _ in 0..natoms {
        topology.add_atom(Atom::new(Element::UNKNOWN));
    }
    Ok(Trajectory::new(topology, Box::new(source))?)
}

/// [`crate::io::format::ByteReadFn`] for TRR.
pub(crate) fn read_trr_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match build_trajectory(bytes.to_vec()) {
        Ok(trajectory) => out.records.push(Record {
            payload: Payload::Frames(trajectory),
            // TRR states no name of its own -- the same `Molecule_{position}`
            // fallback PSF/PRMTOP use, but there is only ever one trajectory
            // per file.
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

fn write_frame(out: &mut Vec<u8>, frame: &Frame) {
    let mut w = XdrWriter::new();
    w.write_i32(GROMACS_MAGIC);
    w.write_i32(TRR_VERSION_STRING.len() as i32 + 1); // the header's own redundant outer length
    w.write_string(TRR_VERSION_STRING);
    w.write_i32(0); // ir_size
    w.write_i32(0); // e_size
    w.write_i32(if frame.cell.is_some() { 36 } else { 0 }); // box_size
    w.write_i32(0); // vir_size
    w.write_i32(0); // pres_size
    w.write_i32(0); // top_size
    w.write_i32(0); // sym_size

    let natoms = frame.num_atoms();
    let coord_size = (natoms * 3 * 4) as i32;
    w.write_i32(coord_size); // x_size -- positions are always present
    w.write_i32(if frame.velocities.is_some() {
        coord_size
    } else {
        0
    });
    w.write_i32(if frame.forces.is_some() {
        coord_size
    } else {
        0
    });
    w.write_i32(natoms as i32);
    w.write_i32(frame.step.unwrap_or(0) as i32);
    w.write_i32(0); // nre
    w.write_f32(frame.time.unwrap_or(0.0) as f32);
    w.write_f32(0.0); // lambda

    if let Some(cell) = frame.cell {
        for v in cell.basis() {
            w.write_f32((v.x / NM_TO_ANGSTROM) as f32);
            w.write_f32((v.y / NM_TO_ANGSTROM) as f32);
            w.write_f32((v.z / NM_TO_ANGSTROM) as f32);
        }
    }
    for p in &frame.positions {
        w.write_f32((p.x / NM_TO_ANGSTROM) as f32);
        w.write_f32((p.y / NM_TO_ANGSTROM) as f32);
        w.write_f32((p.z / NM_TO_ANGSTROM) as f32);
    }
    if let Some(velocities) = &frame.velocities {
        for v in velocities {
            w.write_f32((v.x / NM_TO_ANGSTROM) as f32);
            w.write_f32((v.y / NM_TO_ANGSTROM) as f32);
            w.write_f32((v.z / NM_TO_ANGSTROM) as f32);
        }
    }
    if let Some(forces) = &frame.forces {
        for f in forces {
            w.write_f32((f.x * NM_TO_ANGSTROM) as f32);
            w.write_f32((f.y * NM_TO_ANGSTROM) as f32);
            w.write_f32((f.z * NM_TO_ANGSTROM) as f32);
        }
    }
    out.extend_from_slice(&w.into_bytes());
}

/// [`crate::io::format::ByteWriteTrajectoryFn`] for TRR. Single precision
/// only -- see the module doc.
pub(crate) fn write_trr_bytes(trajectory: &mut Trajectory, _options: &WriteOptions) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..trajectory.frame_count() {
        let frame = trajectory
            .frame(i)
            .expect("a Trajectory already validated its own frame/atom counts at construction");
        write_frame(&mut out, &frame);
    }
    out
}

/// Buffers the whole input and parses it once -- TRR's frame count is not
/// known until the whole file is scanned, the same situation
/// [`crate::io::bcif::BcifSupplier`] is in for a different reason.
pub struct TrrSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl TrrSupplier {
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

impl Iterator for TrrSupplier {
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
        fn frame(&mut self, index: usize) -> io::Result<Frame> {
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

    #[test]
    fn test_positions_only_frames_round_trip() {
        let frames = vec![
            Frame {
                positions: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.5, 0.0, 0.0)],
                velocities: None,
                forces: None,
                time: Some(0.0),
                step: Some(0),
                cell: None,
            },
            Frame {
                positions: vec![Point3::new(0.1, 0.0, 0.0), Point3::new(1.6, 0.0, 0.0)],
                velocities: None,
                forces: None,
                time: Some(0.002),
                step: Some(1),
                cell: None,
            },
        ];
        let mut trajectory = trajectory_from(frames.clone());
        let bytes = write_trr_bytes(&mut trajectory, &WriteOptions::default());
        let outcome = read_trr_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);

        let mut back = as_trajectory(outcome);
        assert_eq!(back.frame_count(), 2);
        for (i, expected) in frames.iter().enumerate() {
            let f = back.frame(i).unwrap();
            for (a, b) in f.positions.iter().zip(&expected.positions) {
                assert!((a.x - b.x).abs() < 1e-3, "{a} vs {b}");
                assert!((a.y - b.y).abs() < 1e-3, "{a} vs {b}");
                assert!((a.z - b.z).abs() < 1e-3, "{a} vs {b}");
            }
            assert!(f.velocities.is_none());
            assert!(f.forces.is_none());
        }
    }

    #[test]
    fn test_velocities_and_forces_round_trip_with_correct_unit_directions() {
        let frame = Frame {
            positions: vec![Point3::new(1.0, 2.0, 3.0)],
            velocities: Some(vec![Point3::new(0.5, -0.5, 0.0)]),
            forces: Some(vec![Point3::new(10.0, 0.0, 0.0)]),
            time: Some(1.0),
            step: Some(5),
            cell: None,
        };
        let mut trajectory = trajectory_from(vec![frame]);
        let bytes = write_trr_bytes(&mut trajectory, &WriteOptions::default());
        let outcome = read_trr_bytes(&bytes, &ReadOptions::default());
        let mut back = as_trajectory(outcome);
        let f = back.frame(0).unwrap();

        assert!((f.positions[0].x - 1.0).abs() < 1e-3);
        let v = f.velocities.expect("velocities survive");
        assert!((v[0].x - 0.5).abs() < 1e-3);
        assert!((v[0].y - (-0.5)).abs() < 1e-3);
        let force = f.forces.expect("forces survive");
        assert!((force[0].x - 10.0).abs() < 1e-2, "{}", force[0].x);
        assert_eq!(f.step, Some(5));
        assert!((f.time.unwrap() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_presence_of_velocities_and_forces_varies_independently_per_frame() {
        let frames = vec![
            Frame {
                positions: vec![Point3::ORIGIN],
                velocities: None,
                forces: None,
                time: Some(0.0),
                step: Some(0),
                cell: None,
            },
            Frame {
                positions: vec![Point3::new(1.0, 0.0, 0.0)],
                velocities: Some(vec![Point3::new(0.1, 0.0, 0.0)]),
                forces: None,
                time: Some(0.002),
                step: Some(1),
                cell: None,
            },
            Frame {
                positions: vec![Point3::new(2.0, 0.0, 0.0)],
                velocities: None,
                forces: Some(vec![Point3::new(5.0, 0.0, 0.0)]),
                time: Some(0.004),
                step: Some(2),
                cell: None,
            },
        ];
        let mut trajectory = trajectory_from(frames);
        let bytes = write_trr_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_trr_bytes(&bytes, &ReadOptions::default()));

        assert_eq!(back.frame_count(), 3);
        assert!(back.frame(0).unwrap().velocities.is_none());
        assert!(back.frame(0).unwrap().forces.is_none());
        assert!(back.frame(1).unwrap().velocities.is_some());
        assert!(back.frame(1).unwrap().forces.is_none());
        assert!(back.frame(2).unwrap().velocities.is_none());
        assert!(back.frame(2).unwrap().forces.is_some());
    }

    #[test]
    fn test_a_triclinic_box_round_trips() {
        let cell = UnitCell::new(20.0, 20.0, 20.0, 80.0, 85.0, 75.0);
        let frame = Frame {
            positions: vec![Point3::ORIGIN],
            velocities: None,
            forces: None,
            time: Some(0.0),
            step: Some(0),
            cell: Some(cell),
        };
        let mut trajectory = trajectory_from(vec![frame]);
        let bytes = write_trr_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_trr_bytes(&bytes, &ReadOptions::default()));
        let f = back.frame(0).unwrap();
        let back_cell = f.cell.expect("cell survives");

        assert!((back_cell.a - cell.a).abs() < 1e-2, "{}", back_cell.a);
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
    fn test_frames_are_randomly_accessible_not_only_in_order() {
        let frames = vec![
            Frame {
                positions: vec![Point3::new(0.0, 0.0, 0.0)],
                velocities: None,
                forces: None,
                time: Some(0.0),
                step: Some(0),
                cell: None,
            },
            Frame {
                positions: vec![Point3::new(1.0, 0.0, 0.0)],
                velocities: None,
                forces: None,
                time: Some(1.0),
                step: Some(1),
                cell: None,
            },
            Frame {
                positions: vec![Point3::new(2.0, 0.0, 0.0)],
                velocities: None,
                forces: None,
                time: Some(2.0),
                step: Some(2),
                cell: None,
            },
        ];
        let mut trajectory = trajectory_from(frames);
        let bytes = write_trr_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_trr_bytes(&bytes, &ReadOptions::default()));

        let last = back.frame(2).unwrap();
        assert!((last.positions[0].x - 2.0).abs() < 1e-3);
        let first = back.frame(0).unwrap();
        assert!((first.positions[0].x - 0.0).abs() < 1e-3);
    }

    #[test]
    fn test_a_bad_magic_number_is_a_clear_error_not_a_panic() {
        let outcome = read_trr_bytes(&[0, 0, 0, 0], &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0].error.contains("1993"),
            "{}",
            outcome.skipped[0].error
        );
    }

    #[test]
    fn test_empty_input_is_a_clear_error_not_a_panic() {
        let outcome = read_trr_bytes(&[], &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
    }

    fn push_i32(buf: &mut Vec<u8>, v: i32) {
        buf.extend_from_slice(&v.to_be_bytes());
    }
    fn push_f64(buf: &mut Vec<u8>, v: f64) {
        buf.extend_from_slice(&v.to_be_bytes());
    }
    /// Pushes the exact bytes `read_header` expects for the version field:
    /// the header's own redundant outer length, then a plain XDR string.
    fn push_version_string(buf: &mut Vec<u8>, s: &str) {
        push_i32(buf, s.len() as i32 + 1);
        push_i32(buf, s.len() as i32);
        buf.extend_from_slice(s.as_bytes());
        let padded = s.len().div_ceil(4) * 4;
        buf.resize(buf.len() + (padded - s.len()), 0);
    }

    #[test]
    fn test_double_precision_frames_are_read_correctly() {
        // This crate's own writer only ever produces single precision (see
        // `write_trr_bytes`'s doc comment) -- a real double-precision
        // GROMACS build does not, so the reader must handle both. Built by
        // hand since `write_trr_bytes` cannot produce this.
        let natoms = 2usize;
        let positions = [Point3::new(1.0, 2.0, 3.0), Point3::new(-1.0, 0.5, 0.25)];
        let mut buf = Vec::new();
        push_i32(&mut buf, GROMACS_MAGIC);
        push_version_string(&mut buf, TRR_VERSION_STRING);
        push_i32(&mut buf, 0); // ir_size
        push_i32(&mut buf, 0); // e_size
        push_i32(&mut buf, 0); // box_size -- no box in this frame
        push_i32(&mut buf, 0); // vir_size
        push_i32(&mut buf, 0); // pres_size
        push_i32(&mut buf, 0); // top_size
        push_i32(&mut buf, 0); // sym_size
        push_i32(&mut buf, (natoms * 3 * 8) as i32); // x_size, double precision
        push_i32(&mut buf, 0); // v_size
        push_i32(&mut buf, 0); // f_size
        push_i32(&mut buf, natoms as i32);
        push_i32(&mut buf, 0); // step
        push_i32(&mut buf, 0); // nre
        push_f64(&mut buf, 0.0); // t
        push_f64(&mut buf, 0.0); // lambda
        for p in &positions {
            push_f64(&mut buf, p.x / NM_TO_ANGSTROM);
            push_f64(&mut buf, p.y / NM_TO_ANGSTROM);
            push_f64(&mut buf, p.z / NM_TO_ANGSTROM);
        }

        let outcome = read_trr_bytes(&buf, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut trajectory = as_trajectory(outcome);
        let frame = trajectory.frame(0).unwrap();
        for (got, want) in frame.positions.iter().zip(&positions) {
            assert!((got.x - want.x).abs() < 1e-9, "{got} vs {want}");
            assert!((got.y - want.y).abs() < 1e-9, "{got} vs {want}");
            assert!((got.z - want.z).abs() < 1e-9, "{got} vs {want}");
        }
    }

    #[test]
    fn test_an_unsupported_precision_is_a_clear_error_not_a_panic() {
        let mut buf = Vec::new();
        push_i32(&mut buf, GROMACS_MAGIC);
        push_version_string(&mut buf, TRR_VERSION_STRING);
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0); // box_size
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 5); // x_size: 5 bytes for 1 atom's 3 reals -- neither 4 nor 8
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 1); // natoms
        push_i32(&mut buf, 0);
        push_i32(&mut buf, 0);

        let outcome = read_trr_bytes(&buf, &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0]
                .error
                .to_lowercase()
                .contains("precision"),
            "{}",
            outcome.skipped[0].error
        );
    }
}
