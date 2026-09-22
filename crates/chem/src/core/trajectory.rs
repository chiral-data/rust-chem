//! One topology, many frames — a trajectory (#311).
//!
//! [`Molecule`] holds exactly one conformer. Reading a trajectory as
//! `Vec<Molecule>` — which is exactly how PDB's multi-`MODEL` and XYZ's and
//! GRO's multi-frame readers already work — stores the atoms, bonds, graph,
//! chains and residues once per frame, and the topology is the part that
//! does not change. [`Trajectory`] holds it once and the frames beside it.

use std::io;

use thiserror::Error;

use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;

/// One conformer of a trajectory's shared topology: positions, and whatever
/// else the format that produced it happened to carry.
///
/// Velocities, forces, simulation time, step number and a per-frame box are
/// each independently `Option`, the same shape [`crate::core::site::AtomSite`]'s
/// fields take and for the same reason: TRR carries velocities and forces,
/// XTC never does, and forcing a value into a column a format never wrote is
/// how a plausible-looking trajectory full of invented zeroes gets produced.
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    /// One position per atom, in the topology's atom order. Always present
    /// — a frame with no positions is not a frame.
    pub positions: Vec<Point3>,
    /// One velocity per atom, same order as `positions`. `None` for a format
    /// that does not carry it — XTC never does.
    pub velocities: Option<Vec<Point3>>,
    /// One force per atom, same order as `positions`. `None` for a format
    /// that does not carry it.
    pub forces: Option<Vec<Point3>>,
    /// Simulation time, in whatever unit the format states.
    pub time: Option<f64>,
    /// Simulation step number.
    pub step: Option<u64>,
    /// This frame's box.
    ///
    /// Per-frame rather than one cell on the shared topology (see
    /// [`Trajectory::new`]'s validation): a format with a single fixed box
    /// for its whole run just repeats the same [`UnitCell`] on every frame,
    /// which costs nothing — six `f64`s, `Copy` — unlike the per-frame
    /// topology duplication this type exists to avoid.
    pub cell: Option<UnitCell>,
}

impl Frame {
    pub fn num_atoms(&self) -> usize {
        self.positions.len()
    }
}

/// Random access to one trajectory's frames over a seekable byte source.
///
/// Implemented once per binary trajectory format (XTC, TRR, DCD, ...) in
/// `crate::io`, in a later story (#325-#329) — this trait only names the
/// contract, so [`Trajectory`] can hold one without `core` depending on
/// `crate::io` and without knowing which format backs it.
///
/// `Read + Seek` rather than `mmap`: `chem-app` ships a
/// `wasm32-unknown-unknown` build with no filesystem, so every implementor
/// has to work over an in-memory `std::io::Cursor<Vec<u8>>` there, not just
/// over a `File`. Nothing in this trait's signature requires that directly —
/// an implementor holds its own `Read + Seek` source internally and builds a
/// frame-offset index from it once, on construction, which is what makes
/// [`Self::frame_count`]/[`Self::num_atoms`] cheap and "read frame 9,000 of
/// 10,000" not require reading the 8,999 before it.
pub trait FrameSource {
    /// How many frames this source has.
    fn frame_count(&self) -> usize;

    /// The number of atoms every frame carries — fixed for the life of the
    /// source, and what [`Trajectory`] cross-checks a fetched [`Frame`]
    /// against.
    fn num_atoms(&self) -> usize;

    /// Reads frame `index` fresh from the underlying source.
    ///
    /// # Errors
    /// Any I/O failure from the seek/read. A caller should also treat an
    /// out-of-range `index` as an error, though [`Trajectory::frame`]
    /// already bounds-checks before calling this, so an implementation
    /// reached only through `Trajectory` need not re-check it.
    fn frame(&mut self, index: usize) -> io::Result<Frame>;
}

/// Errors constructing or reading a [`Trajectory`].
///
/// `#[non_exhaustive]` for the same reason [`crate::core::molecule::MoleculeError`]
/// is: adding a variant to a public enum is a breaking change unless callers
/// are told not to match it exhaustively, and that has to be true from the
/// first published version.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum TrajectoryError {
    /// The topology already carries a conformer of its own — ambiguous,
    /// since a trajectory's frames are its only legitimate source of
    /// positions and there is no way to say which frame the topology's own
    /// conformer would be.
    #[error(
        "the topology already carries its own 3D conformer; a trajectory's \
         frames are the only source of positions, so which frame would it be?"
    )]
    TopologyHasConformer,

    /// The topology already carries a unit cell of its own — a trajectory
    /// carries a cell per frame ([`Frame::cell`]) instead.
    #[error(
        "the topology already carries its own unit cell; a trajectory \
         carries a cell per frame instead"
    )]
    TopologyHasCell,

    #[error("the frame source declares {got} atoms per frame but the topology has {expected}")]
    AtomCountMismatch { expected: usize, got: usize },

    #[error("frame {index} is out of range: this trajectory has {frame_count} frames")]
    FrameIndexOutOfRange { index: usize, frame_count: usize },

    #[error("frame {index} carries {got} positions but the topology has {expected} atoms")]
    FrameAtomCountMismatch {
        index: usize,
        expected: usize,
        got: usize,
    },

    #[error("I/O error reading frame {index}: {source}")]
    Io {
        index: usize,
        #[source]
        source: io::Error,
    },
}

/// One topology, many frames.
///
/// Holds the shared atoms, bonds, graph, chains and residues exactly once,
/// as a [`Molecule`] whose own conformer and cell stay unset (see
/// [`Trajectory::new`]), plus lazy, indexed access to the frames beside it.
///
/// Reuses `Molecule` for the topology rather than a parallel, lighter type:
/// `Molecule`'s per-atom tables are already `Option`, so they simply stay
/// unset for a trajectory's shared topology, and every existing consumer of
/// atoms/bonds/residues keeps working on it unmodified.
pub struct Trajectory {
    topology: Molecule,
    frames: Box<dyn FrameSource>,
}

impl Trajectory {
    /// # Errors
    /// [`TrajectoryError::TopologyHasConformer`]/[`TrajectoryError::TopologyHasCell`]
    /// if `topology` already carries one, and
    /// [`TrajectoryError::AtomCountMismatch`] if `frames.num_atoms()`
    /// disagrees with `topology.num_atoms()`.
    pub fn new(topology: Molecule, frames: Box<dyn FrameSource>) -> Result<Self, TrajectoryError> {
        if topology.has_coords3() {
            return Err(TrajectoryError::TopologyHasConformer);
        }
        if topology.has_cell() {
            return Err(TrajectoryError::TopologyHasCell);
        }
        if frames.num_atoms() != topology.num_atoms() {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: topology.num_atoms(),
                got: frames.num_atoms(),
            });
        }
        Ok(Self { topology, frames })
    }

    pub fn topology(&self) -> &Molecule {
        &self.topology
    }

    pub fn frame_count(&self) -> usize {
        self.frames.frame_count()
    }

    pub fn num_atoms(&self) -> usize {
        self.topology.num_atoms()
    }

    /// Reads frame `index`.
    ///
    /// # Errors
    /// [`TrajectoryError::FrameIndexOutOfRange`], [`TrajectoryError::Io`], or
    /// [`TrajectoryError::FrameAtomCountMismatch`] if the fetched frame's own
    /// position count disagrees with the topology — a [`FrameSource`] that
    /// declares a count it does not actually deliver is a bug this catches
    /// rather than trusts, the same reasoning
    /// [`Molecule::set_coords3`](crate::core::molecule::Molecule::set_coords3)
    /// already applies to its own count.
    pub fn frame(&mut self, index: usize) -> Result<Frame, TrajectoryError> {
        let frame_count = self.frame_count();
        if index >= frame_count {
            return Err(TrajectoryError::FrameIndexOutOfRange { index, frame_count });
        }
        let frame = self
            .frames
            .frame(index)
            .map_err(|source| TrajectoryError::Io { index, source })?;
        let expected = self.topology.num_atoms();
        if frame.positions.len() != expected {
            return Err(TrajectoryError::FrameAtomCountMismatch {
                index,
                expected,
                got: frame.positions.len(),
            });
        }
        Ok(frame)
    }
}

impl std::fmt::Debug for Trajectory {
    /// The topology and the frame/atom counts, not the trait object —
    /// `Box<dyn FrameSource>` has no `Debug` bound and none is worth adding
    /// to the trait for a diagnostic print.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Trajectory")
            .field("topology", &self.topology)
            .field("frame_count", &self.frames.frame_count())
            .field("num_atoms", &self.frames.num_atoms())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{Atom, Element};
    use std::io::{Cursor, Read, Seek, SeekFrom};

    /// A minimal, test-only [`FrameSource`]: each frame on the wire is
    /// `[u32 length][length bytes of little-endian f64 position triples]`,
    /// genuinely variable-sized — proving a real "index of frame offsets
    /// built on open" rather than the arithmetic a fixed-stride format would
    /// let a reader skip. A fixed-stride format (TRR, DCD) is the strictly
    /// simpler case and needs no separate fixture once this one works.
    struct IndexedFrames {
        source: Cursor<Vec<u8>>,
        offsets: Vec<u64>,
        num_atoms: usize,
    }

    impl IndexedFrames {
        fn open(mut source: Cursor<Vec<u8>>, num_atoms: usize) -> io::Result<Self> {
            let mut offsets = Vec::new();
            let mut pos = 0u64;
            let len = source.seek(SeekFrom::End(0))?;
            loop {
                if pos >= len {
                    break;
                }
                source.seek(SeekFrom::Start(pos))?;
                let mut len_buf = [0u8; 4];
                source.read_exact(&mut len_buf)?;
                let payload_len = u32::from_le_bytes(len_buf) as u64;
                offsets.push(pos + 4);
                pos += 4 + payload_len;
            }
            Ok(Self {
                source,
                offsets,
                num_atoms,
            })
        }

        /// Encodes `frames` (one `Vec<Point3>` per frame) into the wire
        /// format [`Self::open`] reads back.
        fn encode(frames: &[Vec<Point3>]) -> Vec<u8> {
            let mut buf = Vec::new();
            for positions in frames {
                let mut payload = Vec::with_capacity(positions.len() * 24);
                for p in positions {
                    payload.extend_from_slice(&p.x.to_le_bytes());
                    payload.extend_from_slice(&p.y.to_le_bytes());
                    payload.extend_from_slice(&p.z.to_le_bytes());
                }
                buf.extend_from_slice(&(payload.len() as u32).to_le_bytes());
                buf.extend_from_slice(&payload);
            }
            buf
        }
    }

    impl FrameSource for IndexedFrames {
        fn frame_count(&self) -> usize {
            self.offsets.len()
        }

        fn num_atoms(&self) -> usize {
            self.num_atoms
        }

        fn frame(&mut self, index: usize) -> io::Result<Frame> {
            self.source.seek(SeekFrom::Start(self.offsets[index]))?;
            let mut positions = Vec::with_capacity(self.num_atoms);
            for _ in 0..self.num_atoms {
                let mut xyz = [0u8; 24];
                self.source.read_exact(&mut xyz)?;
                positions.push(Point3::new(
                    f64::from_le_bytes(xyz[0..8].try_into().unwrap()),
                    f64::from_le_bytes(xyz[8..16].try_into().unwrap()),
                    f64::from_le_bytes(xyz[16..24].try_into().unwrap()),
                ));
            }
            Ok(Frame {
                positions,
                velocities: None,
                forces: None,
                time: None,
                step: None,
                cell: None,
            })
        }
    }

    /// A [`FrameSource`] whose declared atom count does not match what
    /// [`FrameSource::frame`] actually delivers -- proving
    /// [`Trajectory::frame`]'s cross-check fires rather than trusting it.
    struct LyingFrameSource;

    impl FrameSource for LyingFrameSource {
        fn frame_count(&self) -> usize {
            1
        }
        fn num_atoms(&self) -> usize {
            1
        }
        fn frame(&mut self, _index: usize) -> io::Result<Frame> {
            // Two positions, though `num_atoms` declared one.
            Ok(Frame {
                positions: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                velocities: None,
                forces: None,
                time: None,
                step: None,
                cell: None,
            })
        }
    }

    fn two_atom_topology() -> Molecule {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));
        mol
    }

    fn sample_frames() -> Vec<Vec<Point3>> {
        vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.1, 0.0, 0.0), Point3::new(1.1, 0.0, 0.0)],
            vec![Point3::new(0.2, 0.0, 0.0), Point3::new(1.2, 0.0, 0.0)],
        ]
    }

    fn indexed_frames() -> IndexedFrames {
        let bytes = IndexedFrames::encode(&sample_frames());
        IndexedFrames::open(Cursor::new(bytes), 2).expect("well-formed fixture")
    }

    #[test]
    fn test_a_trajectory_reads_back_first_middle_and_last_frames() {
        let mut trajectory =
            Trajectory::new(two_atom_topology(), Box::new(indexed_frames())).expect("valid");
        assert_eq!(trajectory.frame_count(), 3);

        let expected = sample_frames();
        for index in [0, 1, 2] {
            let frame = trajectory.frame(index).expect("in range");
            assert_eq!(frame.positions, expected[index]);
            assert!(frame.velocities.is_none());
            assert!(frame.cell.is_none());
        }
    }

    #[test]
    fn test_an_out_of_range_frame_is_refused() {
        let mut trajectory =
            Trajectory::new(two_atom_topology(), Box::new(indexed_frames())).expect("valid");
        match trajectory.frame(3) {
            Err(TrajectoryError::FrameIndexOutOfRange { index, frame_count }) => {
                assert_eq!(index, 3);
                assert_eq!(frame_count, 3);
            }
            other => panic!("expected FrameIndexOutOfRange, got {other:?}"),
        }
    }

    #[test]
    fn test_a_topology_with_its_own_conformer_is_refused() {
        let mut topology = two_atom_topology();
        topology
            .set_coords3(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)])
            .expect("valid coords");
        match Trajectory::new(topology, Box::new(indexed_frames())) {
            Err(TrajectoryError::TopologyHasConformer) => {}
            other => panic!("expected TopologyHasConformer, got {other:?}"),
        }
    }

    #[test]
    fn test_a_topology_with_its_own_cell_is_refused() {
        let mut topology = two_atom_topology();
        topology
            .set_cell(UnitCell::cubic(10.0))
            .expect("valid cell");
        match Trajectory::new(topology, Box::new(indexed_frames())) {
            Err(TrajectoryError::TopologyHasCell) => {}
            other => panic!("expected TopologyHasCell, got {other:?}"),
        }
    }

    #[test]
    fn test_a_frame_source_disagreeing_with_the_topology_is_refused() {
        // The frame source declares 2 atoms; the topology (built here with
        // three) disagrees before a single frame is ever read.
        let mut topology = two_atom_topology();
        topology.add_atom(Atom::new(Element::nitrogen()));
        match Trajectory::new(topology, Box::new(indexed_frames())) {
            Err(TrajectoryError::AtomCountMismatch { expected, got }) => {
                assert_eq!(expected, 3);
                assert_eq!(got, 2);
            }
            other => panic!("expected AtomCountMismatch, got {other:?}"),
        }
    }

    #[test]
    fn test_a_lying_frame_source_is_caught_on_read_not_trusted() {
        // `LyingFrameSource::num_atoms` agrees with a one-atom topology, so
        // `Trajectory::new` accepts it -- the lie only shows up once a frame
        // is actually read and its position count is checked for real.
        let mut topology = Molecule::new();
        topology.add_atom(Atom::new(Element::carbon()));
        let mut trajectory =
            Trajectory::new(topology, Box::new(LyingFrameSource)).expect("counts agree up front");
        match trajectory.frame(0) {
            Err(TrajectoryError::FrameAtomCountMismatch {
                index,
                expected,
                got,
            }) => {
                assert_eq!(index, 0);
                assert_eq!(expected, 1);
                assert_eq!(got, 2);
            }
            other => panic!("expected FrameAtomCountMismatch, got {other:?}"),
        }
    }
}
