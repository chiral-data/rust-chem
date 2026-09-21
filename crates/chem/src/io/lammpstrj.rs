//! LAMMPS's text dump trajectory (#329): `ITEM:`-delimited sections per
//! frame -- `TIMESTEP`, `NUMBER OF ATOMS`, `BOX BOUNDS`, then `ATOMS` with a
//! column list stated fresh on every frame's own header line, not once for
//! the whole file. The first `Kind::Frames` format that's text rather than
//! binary.
//!
//! **The column list is read every frame, never cached** -- nothing stops a
//! later frame from naming a different set of columns (e.g. velocities
//! appearing partway through a run), and the issue's own text calls this
//! out directly.
//!
//! **Four coordinate conventions share the same three-letter shape**:
//! `x/y/z` (unscaled), `xs/ys/zs` (scaled to the box), `xu/yu/zu`
//! (unwrapped), `ix/iy/iz` (image flags -- not a position by themselves,
//! out of scope). Convention priority when more than one is present:
//! unscaled, then scaled, then unwrapped -- confirmed against MDAnalysis's
//! real `DumpReader` (`coordinates/LAMMPS.py`, vendored locally), which
//! resolves the same three in the same order (plus `xsu/ysu/zsu`, out of
//! scope here -- not in the issue's text, and rare). Scaled coordinates
//! convert via [`crate::core::cell::UnitCell::to_cartesian`] plus the
//! frame's own box origin, kept locally for this conversion only --
//! `UnitCell` itself has no origin field, the same disclosed loss
//! [`crate::io::lammps`]'s module doc already states for data files.
//!
//! **`id`-based reordering.** Real dump files are frequently not
//! atom-id-ordered (parallel-MPI output writes whatever order a process
//! happens to hold) -- confirmed as a real concern, not a theoretical one,
//! by MDAnalysis's own reader doing exactly this (`np.argsort` on the `id`
//! column) before trusting row order. When an `id` column is present, rows
//! are sorted by it before being placed into the position array; when it
//! is absent, file order is trusted as-is, the same fallback MDAnalysis
//! itself uses.
//!
//! **Triclinic `BOX BOUNDS` states bounding-box values, not box edges** --
//! LAMMPS shifts the reported `xlo_bound`/`xhi_bound` etc. by the tilt
//! factors. Confirmed against two independent real implementations,
//! algebraically proven equivalent to each other during this story's
//! research: MDAnalysis's `DumpReader` and ASE's
//! `ase.io.lammpsrun.construct_cell`. Once inverted back to true box edges,
//! the shape conversion itself is
//! `crate::io::lammps::unit_cell_from_lammps_box` -- the exact same
//! closed-form math LAMMPS *data* files already use, since both share one
//! box representation.
//!
//! **No `Frame::time`.** A dump's `TIMESTEP` is a step count, not a time --
//! LAMMPS's `units`/`timestep` command lives in the input script, never the
//! dump itself, the same "units unrecoverable from the file alone"
//! precedent already stated for data files. `Frame::step` is populated;
//! `Frame::time` stays `None`, and [`crate::io::format::Carries::FRAME_TIME`]
//! (which names `Frame::time` specifically) is not in this format's mask.
//!
//! **Writing** always uses unscaled coordinates -- sidesteps the exact
//! ambiguity this module's own reader has to resolve, the same
//! least-ambiguous-convention posture NCTRAJ took always writing
//! `"angstrom"`. `id` is written `1..=natoms` in topology order (a round
//! trip through this crate alone needs no sorting); `type` is written as a
//! constant `1` for every atom -- this format's shared topology carries no
//! real LAMMPS type numbers, the same bare-topology shape every other
//! trajectory format's shared topology already has.

use crate::core::atom::{Atom, Element};
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::trajectory::{Frame, FrameSource, Trajectory};
use crate::io::errors::{LammpstrjError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

fn next_line<'a>(
    lines: &mut std::iter::Peekable<std::str::Lines<'a>>,
) -> Result<&'a str, LammpstrjError> {
    lines
        .next()
        .ok_or_else(|| LammpstrjError::ParseError("unexpected end of file".to_string()))
}

fn expect_item<'a>(
    lines: &mut std::iter::Peekable<std::str::Lines<'a>>,
    expected: &str,
) -> Result<&'a str, LammpstrjError> {
    let line = lines
        .next()
        .ok_or_else(|| LammpstrjError::MissingItemHeader(expected.to_string()))?;
    let trimmed = line.trim();
    let prefix = format!("ITEM: {expected}");
    if trimmed == prefix || trimmed.starts_with(&format!("{prefix} ")) {
        Ok(trimmed)
    } else {
        Err(LammpstrjError::MissingItemHeader(expected.to_string()))
    }
}

fn next_f64<'a>(it: &mut impl Iterator<Item = &'a str>) -> Result<f64, LammpstrjError> {
    it.next()
        .ok_or_else(|| LammpstrjError::ParseError("missing box bound value".to_string()))?
        .parse()
        .map_err(|_| LammpstrjError::ParseError("invalid box bound value".to_string()))
}

fn parse_two(line: &str) -> Result<(f64, f64), LammpstrjError> {
    let mut it = line.split_whitespace();
    Ok((next_f64(&mut it)?, next_f64(&mut it)?))
}

fn parse_three(line: &str) -> Result<(f64, f64, f64), LammpstrjError> {
    let mut it = line.split_whitespace();
    Ok((next_f64(&mut it)?, next_f64(&mut it)?, next_f64(&mut it)?))
}

fn min4(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a.min(b).min(c).min(d)
}

fn max4(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a.max(b).max(c).max(d)
}

#[derive(Clone, Copy)]
enum Coord {
    Unscaled([usize; 3]),
    Scaled([usize; 3]),
    Unwrapped([usize; 3]),
}

fn parse_one_frame<'a>(
    lines: &mut std::iter::Peekable<std::str::Lines<'a>>,
) -> Result<Frame, LammpstrjError> {
    expect_item(lines, "TIMESTEP")?;
    let step: u64 = next_line(lines)?
        .trim()
        .parse()
        .map_err(|_| LammpstrjError::ParseError("invalid TIMESTEP value".to_string()))?;

    expect_item(lines, "NUMBER OF ATOMS")?;
    let natoms: usize = next_line(lines)?
        .trim()
        .parse()
        .map_err(|_| LammpstrjError::ParseError("invalid NUMBER OF ATOMS value".to_string()))?;

    let box_header = expect_item(lines, "BOX BOUNDS")?;
    let box_tokens: Vec<&str> = box_header.split_whitespace().collect();
    let triclinic =
        box_tokens.contains(&"xy") && box_tokens.contains(&"xz") && box_tokens.contains(&"yz");

    let (cell, xlo, ylo, zlo) = if triclinic {
        let (xlo_b, xhi_b, xy) = parse_three(next_line(lines)?)?;
        let (ylo_b, yhi_b, xz) = parse_three(next_line(lines)?)?;
        let (zlo, zhi, yz) = parse_three(next_line(lines)?)?;
        // Inverts the dump format's bounding-box shift back to the true
        // box edges -- see the module doc's two-independent-readers note.
        let xlo = xlo_b - min4(0.0, xy, xz, xy + xz);
        let xhi = xhi_b - max4(0.0, xy, xz, xy + xz);
        let ylo = ylo_b - yz.min(0.0);
        let yhi = yhi_b - yz.max(0.0);
        let cell = crate::io::lammps::unit_cell_from_lammps_box(
            xhi - xlo,
            yhi - ylo,
            zhi - zlo,
            xy,
            xz,
            yz,
        );
        (cell, xlo, ylo, zlo)
    } else {
        let (xlo, xhi) = parse_two(next_line(lines)?)?;
        let (ylo, yhi) = parse_two(next_line(lines)?)?;
        let (zlo, zhi) = parse_two(next_line(lines)?)?;
        let cell = crate::io::lammps::unit_cell_from_lammps_box(
            xhi - xlo,
            yhi - ylo,
            zhi - zlo,
            0.0,
            0.0,
            0.0,
        );
        (cell, xlo, ylo, zlo)
    };

    let atoms_header = expect_item(lines, "ATOMS")?;
    let columns: Vec<&str> = atoms_header
        .strip_prefix("ITEM: ATOMS")
        .unwrap_or("")
        .split_whitespace()
        .collect();
    let col = |name: &str| columns.iter().position(|&c| c == name);

    let coord = if let (Some(a), Some(b), Some(c)) = (col("x"), col("y"), col("z")) {
        Coord::Unscaled([a, b, c])
    } else if let (Some(a), Some(b), Some(c)) = (col("xs"), col("ys"), col("zs")) {
        Coord::Scaled([a, b, c])
    } else if let (Some(a), Some(b), Some(c)) = (col("xu"), col("yu"), col("zu")) {
        Coord::Unwrapped([a, b, c])
    } else {
        return Err(LammpstrjError::MissingCoordinateColumns);
    };
    let id_col = col("id");
    let vel_cols = match (col("vx"), col("vy"), col("vz")) {
        (Some(a), Some(b), Some(c)) => Some([a, b, c]),
        _ => None,
    };
    let force_cols = match (col("fx"), col("fy"), col("fz")) {
        (Some(a), Some(b), Some(c)) => Some([a, b, c]),
        _ => None,
    };

    let mut rows: Vec<(i64, Vec<f64>)> = Vec::with_capacity(natoms);
    for _ in 0..natoms {
        let line = next_line(lines)?;
        let fields: Vec<f64> = line
            .split_whitespace()
            .map(|tok| {
                tok.parse::<f64>()
                    .map_err(|_| LammpstrjError::ParseError(format!("invalid atom field {tok:?}")))
            })
            .collect::<Result<_, _>>()?;
        let id = id_col
            .and_then(|i| fields.get(i))
            .map(|v| v.round() as i64)
            .unwrap_or(0);
        rows.push((id, fields));
    }
    if id_col.is_some() {
        rows.sort_by_key(|(id, _)| *id);
    }

    let get = |fields: &[f64], i: usize| -> Result<f64, LammpstrjError> {
        fields
            .get(i)
            .copied()
            .ok_or_else(|| LammpstrjError::ParseError("atom line has too few columns".to_string()))
    };

    let mut positions = Vec::with_capacity(natoms);
    let mut velocities = vel_cols.map(|_| Vec::with_capacity(natoms));
    let mut forces = force_cols.map(|_| Vec::with_capacity(natoms));

    for (_, fields) in &rows {
        let position = match coord {
            Coord::Unscaled(c) | Coord::Unwrapped(c) => {
                Point3::new(get(fields, c[0])?, get(fields, c[1])?, get(fields, c[2])?)
            }
            Coord::Scaled(c) => {
                let frac = Point3::new(get(fields, c[0])?, get(fields, c[1])?, get(fields, c[2])?);
                cell.to_cartesian(frac) + Point3::new(xlo, ylo, zlo)
            }
        };
        positions.push(position);
        if let (Some(cols), Some(out)) = (vel_cols, velocities.as_mut()) {
            out.push(Point3::new(
                get(fields, cols[0])?,
                get(fields, cols[1])?,
                get(fields, cols[2])?,
            ));
        }
        if let (Some(cols), Some(out)) = (force_cols, forces.as_mut()) {
            out.push(Point3::new(
                get(fields, cols[0])?,
                get(fields, cols[1])?,
                get(fields, cols[2])?,
            ));
        }
    }

    Ok(Frame {
        positions,
        velocities,
        forces,
        time: None,
        step: Some(step),
        cell: Some(cell),
    })
}

fn parse_frames(text: &str) -> Result<Vec<Frame>, LammpstrjError> {
    let mut lines = text.lines().peekable();
    let mut frames = Vec::new();
    while lines.peek().is_some() {
        while matches!(lines.peek(), Some(l) if l.trim().is_empty()) {
            lines.next();
        }
        if lines.peek().is_none() {
            break;
        }
        frames.push(parse_one_frame(&mut lines)?);
    }
    Ok(frames)
}

struct LammpstrjFrameSource(Vec<Frame>);

impl FrameSource for LammpstrjFrameSource {
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

fn build_trajectory(text: &str) -> Result<Trajectory, LammpstrjError> {
    let frames = parse_frames(text)?;
    if frames.is_empty() {
        return Err(LammpstrjError::ParseError("no frames found".to_string()));
    }
    let natoms = frames[0].num_atoms();
    let mut topology = Molecule::new();
    for _ in 0..natoms {
        topology.add_atom(Atom::new(Element::UNKNOWN));
    }
    Ok(Trajectory::new(
        topology,
        Box::new(LammpstrjFrameSource(frames)),
    )?)
}

/// [`crate::io::format::ReadFn`] for LAMMPS Trajectory -- plain text, unlike
/// every other `Kind::Frames` format so far, so this uses the same `&str`
/// entry point every `Kind::Molecules` text format already does rather than
/// the byte-level one (that keeps `Format::encoding` and the presence of a
/// byte reader in agreement, an invariant [`crate::io::format`]'s own tests
/// enforce across every non-binary format).
pub(crate) fn read_lammpstrj(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match build_trajectory(text) {
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

fn write_frame(out: &mut String, frame: &Frame) {
    out.push_str("ITEM: TIMESTEP\n");
    out.push_str(&format!("{}\n", frame.step.unwrap_or(0)));
    out.push_str("ITEM: NUMBER OF ATOMS\n");
    out.push_str(&format!("{}\n", frame.num_atoms()));

    let (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz) =
        crate::io::lammps::lammps_box_bounds(frame.cell.as_ref(), &frame.positions);
    let triclinic = xy.abs() > 1e-9 || xz.abs() > 1e-9 || yz.abs() > 1e-9;
    if triclinic {
        out.push_str("ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
        let xlo_b = xlo + min4(0.0, xy, xz, xy + xz);
        let xhi_b = xhi + max4(0.0, xy, xz, xy + xz);
        let ylo_b = ylo + yz.min(0.0);
        let yhi_b = yhi + yz.max(0.0);
        out.push_str(&format!("{xlo_b} {xhi_b} {xy}\n"));
        out.push_str(&format!("{ylo_b} {yhi_b} {xz}\n"));
        out.push_str(&format!("{zlo} {zhi} {yz}\n"));
    } else {
        out.push_str("ITEM: BOX BOUNDS pp pp pp\n");
        out.push_str(&format!("{xlo} {xhi}\n"));
        out.push_str(&format!("{ylo} {yhi}\n"));
        out.push_str(&format!("{zlo} {zhi}\n"));
    }

    let mut header = String::from("ITEM: ATOMS id type x y z");
    if frame.velocities.is_some() {
        header.push_str(" vx vy vz");
    }
    if frame.forces.is_some() {
        header.push_str(" fx fy fz");
    }
    out.push_str(&header);
    out.push('\n');

    for i in 0..frame.num_atoms() {
        let p = frame.positions[i];
        out.push_str(&format!("{} 1 {:.6} {:.6} {:.6}", i + 1, p.x, p.y, p.z));
        if let Some(vs) = &frame.velocities {
            let v = vs[i];
            out.push_str(&format!(" {:.6} {:.6} {:.6}", v.x, v.y, v.z));
        }
        if let Some(fs) = &frame.forces {
            let f = fs[i];
            out.push_str(&format!(" {:.6} {:.6} {:.6}", f.x, f.y, f.z));
        }
        out.push('\n');
    }
}

/// [`crate::io::format::ByteWriteTrajectoryFn`] for LAMMPS Trajectory.
/// Always unscaled coordinates -- see the module doc.
pub(crate) fn write_lammpstrj_bytes(
    trajectory: &mut Trajectory,
    _options: &WriteOptions,
) -> Vec<u8> {
    let nframes = trajectory.frame_count();
    let mut out = String::new();
    for i in 0..nframes {
        let frame = trajectory
            .frame(i)
            .expect("a Trajectory already validated its own frame/atom counts at construction");
        write_frame(&mut out, &frame);
    }
    out.into_bytes()
}

/// Buffers the whole input and parses it once, mirroring
/// [`crate::io::dcd::DcdSupplier`]/[`crate::io::nctraj::NctrajSupplier`].
pub struct LammpstrjSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl LammpstrjSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match reader.read_to_end(&mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match std::str::from_utf8(&bytes) {
                Err(e) => vec![Err(ReadError::Parse {
                    position: 1,
                    message: e.to_string(),
                })],
                Ok(text) => match build_trajectory(text) {
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
            },
        };
        Self {
            records: records.into_iter(),
        }
    }
}

impl Iterator for LammpstrjSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::cell::UnitCell;
    use crate::io::options::ReadOptions;

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

    fn frame(
        positions: Vec<Point3>,
        velocities: Option<Vec<Point3>>,
        forces: Option<Vec<Point3>>,
        step: Option<u64>,
        cell: Option<UnitCell>,
    ) -> Frame {
        Frame {
            positions,
            velocities,
            forces,
            time: None,
            step,
            cell,
        }
    }

    #[test]
    fn test_an_orthogonal_round_trip_with_velocities_and_forces() {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let velocities = vec![
            Point3::new(0.1, 0.0, 0.0),
            Point3::new(0.0, 0.1, 0.0),
            Point3::new(0.0, 0.0, 0.1),
        ];
        let forces = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let cell = UnitCell::new(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
        let mut trajectory = trajectory_from(vec![
            frame(
                positions.clone(),
                Some(velocities.clone()),
                Some(forces.clone()),
                Some(0),
                Some(cell),
            ),
            frame(
                positions
                    .iter()
                    .map(|p| *p + Point3::new(0.5, 0.0, 0.0))
                    .collect(),
                Some(velocities.clone()),
                Some(forces.clone()),
                Some(1),
                Some(cell),
            ),
        ]);
        let bytes = write_lammpstrj_bytes(&mut trajectory, &WriteOptions::default());
        let text = std::str::from_utf8(&bytes).expect("valid UTF-8");
        let outcome = read_lammpstrj(text, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut back = as_trajectory(outcome);
        assert_eq!(back.frame_count(), 2);

        let f0 = back.frame(0).unwrap();
        for (a, b) in f0.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 1e-4);
            assert!((a.y - b.y).abs() < 1e-4);
            assert!((a.z - b.z).abs() < 1e-4);
        }
        assert_eq!(f0.step, Some(0));
        assert!(f0.time.is_none());
        let back_vel = f0.velocities.expect("velocities survive");
        for (a, b) in back_vel.iter().zip(&velocities) {
            assert!((a.x - b.x).abs() < 1e-4);
            assert!((a.y - b.y).abs() < 1e-4);
            assert!((a.z - b.z).abs() < 1e-4);
        }
        let back_force = f0.forces.expect("forces survive");
        for (a, b) in back_force.iter().zip(&forces) {
            assert!((a.x - b.x).abs() < 1e-4);
        }
        let back_cell = f0.cell.expect("cell survives");
        assert!((back_cell.a - 10.0).abs() < 1e-3);
        assert!((back_cell.alpha - 90.0).abs() < 1e-3);
    }

    #[test]
    fn test_a_triclinic_box_matches_independently_verified_values() {
        // Hand-computed, and cross-checked against two real readers
        // (MDAnalysis's `DumpReader`, ASE's `read_lammps_dump_text`) during
        // this story's research -- not derived from the code under test.
        // True box edges lx=ly=lz=10, tilts xy=2, xz=1, yz=0.5; the dump's
        // own bounding-box shift turns that into the bounds below.
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS xy xz yz pp pp pp
0 13 2
0 10.5 1
0 10 0.5
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
";
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        let f0 = trajectory.frame(0).unwrap();
        let cell = f0.cell.expect("has a cell");
        assert!((cell.a - 10.0).abs() < 1e-6, "{}", cell.a);
        assert!((cell.b - 10.198_039_027_185_569).abs() < 1e-6, "{}", cell.b);
        assert!((cell.c - 10.062_305_898_749_054).abs() < 1e-6, "{}", cell.c);
        assert!(
            (cell.alpha - 86.088_495_040_323_6).abs() < 1e-6,
            "{}",
            cell.alpha
        );
        assert!(
            (cell.beta - 84.296_484_743_617).abs() < 1e-6,
            "{}",
            cell.beta
        );
        assert!(
            (cell.gamma - 78.690_067_525_979_79).abs() < 1e-6,
            "{}",
            cell.gamma
        );
    }

    #[test]
    fn test_scaled_coordinates_convert_using_the_box() {
        // Hand-computed and cross-checked against MDAnalysis's `DumpReader`
        // during this story's research: box (10, 20, 30), scaled
        // (0.5, 0.25, 0.1) -> real (5.0, 5.0, 3.0).
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type xs ys zs
1 1 0.5 0.25 0.1
";
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        let f0 = trajectory.frame(0).unwrap();
        assert!(
            (f0.positions[0].x - 5.0).abs() < 1e-6,
            "{}",
            f0.positions[0].x
        );
        assert!(
            (f0.positions[0].y - 5.0).abs() < 1e-6,
            "{}",
            f0.positions[0].y
        );
        assert!(
            (f0.positions[0].z - 3.0).abs() < 1e-6,
            "{}",
            f0.positions[0].z
        );
    }

    #[test]
    fn test_scaled_coordinates_in_a_triclinic_box_use_its_own_edge_vectors() {
        // A committed fixture and its generator (#341), combining what
        // `test_a_triclinic_box_matches_independently_verified_values`
        // and `test_scaled_coordinates_convert_using_the_box` each hand-
        // verify in isolation: a scaled coordinate's conversion genuinely
        // depends on the box being triclinic, not just orthogonal --
        // `real = origin + xs*a + ys*b + zs*c` over the triclinic edge
        // vectors, not `real = origin + frac * (hi - lo)`. Same triclinic
        // box both source tests already pin (lx=ly=lz=10, xy=2, xz=1,
        // yz=0.5) and the same fractional position (0.5, 0.25, 0.1) the
        // orthogonal-box test used; hand-derived here to
        // `(5.6, 2.55, 1.0)` and cross-checked against MDAnalysis's
        // `DumpReader` during this story's research. See
        // `tests/corpus/lammpstrj/generate_triclinic_scaled.py`.
        let text = include_str!("../../tests/corpus/lammpstrj/triclinic_scaled.lammpstrj");
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        let f0 = trajectory.frame(0).unwrap();
        assert!(
            (f0.positions[0].x - 5.6).abs() < 1e-6,
            "{}",
            f0.positions[0].x
        );
        assert!(
            (f0.positions[0].y - 2.55).abs() < 1e-6,
            "{}",
            f0.positions[0].y
        );
        assert!(
            (f0.positions[0].z - 1.0).abs() < 1e-6,
            "{}",
            f0.positions[0].z
        );
    }

    #[test]
    fn test_unwrapped_coordinates_pass_through_directly() {
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type xu yu zu
1 1 15.5 -3.25 100.0
";
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        let f0 = trajectory.frame(0).unwrap();
        assert!((f0.positions[0].x - 15.5).abs() < 1e-6);
        assert!((f0.positions[0].y - (-3.25)).abs() < 1e-6);
        assert!((f0.positions[0].z - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_the_atoms_header_is_read_fresh_every_frame() {
        // Frame 0 has no velocities; frame 1 does, with columns in a
        // different order -- proving the column map is rebuilt per frame,
        // not cached from frame 0 (the issue's own explicit warning).
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
1 1 1.0 2.0 3.0
ITEM: TIMESTEP
1
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id vx vy vz type x y z
1 9.0 8.0 7.0 1 1.5 2.5 3.5
";
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        let f0 = trajectory.frame(0).unwrap();
        assert!(f0.velocities.is_none());
        let f1 = trajectory.frame(1).unwrap();
        let v1 = f1.velocities.expect("frame 1 has velocities");
        assert!((v1[0].x - 9.0).abs() < 1e-6);
        assert!((f1.positions[0].x - 1.5).abs() < 1e-6);
    }

    #[test]
    fn test_out_of_order_atom_ids_are_sorted_before_use() {
        // File order is id 2 then id 1; sorted by id, index 0 must hold
        // id 1's position and index 1 must hold id 2's.
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
2 1 20.0 20.0 20.0
1 1 10.0 10.0 10.0
";
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        let f0 = trajectory.frame(0).unwrap();
        assert!(
            (f0.positions[0].x - 10.0).abs() < 1e-6,
            "{}",
            f0.positions[0].x
        );
        assert!(
            (f0.positions[1].x - 20.0).abs() < 1e-6,
            "{}",
            f0.positions[1].x
        );
    }

    #[test]
    fn test_missing_coordinate_columns_is_refused() {
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type q
1 1 0.5
";
        let outcome = read_lammpstrj(text, &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0].error.contains("coordinate"),
            "{}",
            outcome.skipped[0].error
        );
    }

    #[test]
    fn test_a_changing_atom_count_is_refused() {
        // Construction only checks frame 0 against the shared topology;
        // `Trajectory::frame` itself is what catches a later frame that
        // declared a different count (see `core::trajectory`).
        let text = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
ITEM: TIMESTEP
1
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
3 1 2.0 0.0 0.0
";
        let mut trajectory = as_trajectory(read_lammpstrj(text, &ReadOptions::default()));
        assert!(trajectory.frame(0).is_ok());
        let err = trajectory.frame(1).unwrap_err();
        assert!(
            matches!(
                err,
                crate::core::trajectory::TrajectoryError::FrameAtomCountMismatch { .. }
            ),
            "{err}"
        );
    }
}
