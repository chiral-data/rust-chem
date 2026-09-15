//! A minimal `Read + Seek`-backed [`FrameSource`], and reading frames 0, mid
//! and last out of it.
//!
//! ```sh
//! cargo run --example trajectory_frame_source
//! ```
//!
//! No registered format produces a [`Trajectory`] yet (#311 is the container
//! story; a real trajectory format is #325-#329), so this is where the
//! mechanism can actually be seen working: an index of frame offsets built
//! once when the source opens, then random access to any frame by seeking
//! straight to it -- never reading the frames before it. `chem-app`'s
//! `wasm32-unknown-unknown` build has no filesystem, so a real backend there
//! wraps browser-handed bytes in exactly the `std::io::Cursor` this example
//! uses, rather than a `std::fs::File`.

use chem::core::prelude::*;
use std::io::{self, Cursor, Read, Seek, SeekFrom};

/// Each frame on the wire is `[u32 length][length bytes of little-endian f64
/// position triples]` -- genuinely variable-sized, unlike a fixed-stride
/// format (TRR, DCD), which needs no index at all beyond a stride multiply.
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
        while pos < len {
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Ethane: two atoms, no positions of its own -- `Trajectory::new` would
    // refuse a topology that already carried a conformer.
    let mut topology = Molecule::new();
    topology.add_atom(Atom::new(Element::carbon()));
    topology.add_atom(Atom::new(Element::carbon()));

    let frames: Vec<Vec<Point3>> = (0..10)
        .map(|i| {
            let x = i as f64 * 0.1;
            vec![Point3::new(x, 0.0, 0.0), Point3::new(x + 1.5, 0.0, 0.0)]
        })
        .collect();
    let bytes = IndexedFrames::encode(&frames);
    println!(
        "encoded {} frames into {} bytes (variable-length records)",
        frames.len(),
        bytes.len()
    );

    let source = IndexedFrames::open(Cursor::new(bytes), 2)?;
    let mut trajectory = Trajectory::new(topology, Box::new(source))?;
    println!(
        "trajectory: {} atoms, {} frames",
        trajectory.num_atoms(),
        trajectory.frame_count()
    );

    let last = trajectory.frame_count() - 1;
    for index in [0, trajectory.frame_count() / 2, last] {
        let frame = trajectory.frame(index)?;
        println!("frame {index}: {:?}", frame.positions);
    }

    Ok(())
}
