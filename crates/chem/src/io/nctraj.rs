//! NCTRAJ — Amber's NetCDF trajectory convention (#328), layered on top of
//! [`crate::io::netcdf3`]'s generic, Amber-unaware NetCDF-3 container.
//!
//! Confirmed field-for-field against a real, validated implementation:
//! `chemfiles`'s `AmberNetCDF.cpp` (vendored in the local cargo registry
//! cache under `chemfiles-sys`), not reconstructed from the wider Amber
//! NetCDF specification's prose. Scope is exactly what this format's own
//! issue names — `coordinates`, optional `velocities`, optional
//! `cell_lengths`/`cell_angles` together, and the `Conventions == "AMBER"`
//! check — not the spec's optional `time`/`forces` variables, which
//! chemfiles' own real implementation doesn't handle either.
//!
//! **Layout**: dimensions `frame` (the record dimension), `atom`,
//! `spatial` (3), and — only when a cell is present — `cell_spatial` (3),
//! `cell_angular` (3), and a `label` (5) dimension used only by the
//! `cell_angular` label variable. Label variables `spatial`/`cell_spatial`/
//! `cell_angular` spell out `"xyz"`/`"abc"`/`"alpha"`+`"beta "`+`"gamma"`.
//! `coordinates(frame, atom, spatial)` (units `"angstrom"`), optional
//! `velocities(frame, atom, spatial)` (units `"angstrom/picosecond"`), and
//! `cell_lengths(frame, cell_spatial)`/`cell_angles(frame, cell_angular)`
//! (units `"angstrom"`/`"degree"`) — required together, never just one
//! (stricter than chemfiles' own slightly asymmetric rule, matching this
//! crate's preference for an explicit error over a silently-accepted
//! half-complete cell). Global attributes `Conventions` (must read
//! `"AMBER"` exactly — the "refused by name" check this format's own issue
//! asks for) and `ConventionVersion` (must read `"1.0"`).
//!
//! **Coordinates are already Å — no unit conversion**, unlike TRR/XTC's
//! nm. The reader still honours a `units` attribute stating something
//! else (a small distance/velocity/angle table, mirroring chemfiles' own
//! — corrected here for nanometres, which chemfiles' own table gives the
//! wrong scale factor for; this crate's own `crate::core::units`
//! constant is used instead of trusting that value uncritically).
//!
//! No `FRAME_TIME`/`FORCES` in this format's [`crate::io::format::Carries`]
//! mask — both are out of scope, so `Frame::time`/`Frame::step` are always
//! `None` and `Frame::forces` is always `None`.

use crate::core::atom::{Atom, Element};
use crate::core::cell::UnitCell;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::trajectory::{Frame, FrameSource, Trajectory};
use crate::io::errors::{NctrajError, ReadError};
use crate::io::netcdf3::{Netcdf3Dimension, Netcdf3File, Netcdf3Value, Netcdf3Variable};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

fn units_of(var: &Netcdf3Variable) -> &str {
    var.attributes
        .iter()
        .find(|(name, _)| name == "units")
        .and_then(|(_, v)| v.as_str())
        .unwrap_or("angstrom")
}

fn distance_scale(units: &str) -> f64 {
    match units.to_ascii_lowercase().as_str() {
        "angstrom" | "angstroms" | "a" | "" => 1.0,
        "nanometer" | "nanometers" | "nm" => crate::core::units::NM_TO_ANGSTROM,
        "bohr" | "bohrs" => 0.529_177_210_903,
        "meter" | "meters" | "m" => 1e10,
        // Unknown unit: best-effort, matching chemfiles' own "warn and
        // default to 1.0" rather than refusing the whole file over it.
        _ => 1.0,
    }
}

fn angle_scale(units: &str) -> f64 {
    match units.to_ascii_lowercase().as_str() {
        "degree" | "degrees" | "" => 1.0,
        "radian" | "radians" => 180.0 / std::f64::consts::PI,
        _ => 1.0,
    }
}

fn velocity_scale(units: &str) -> f64 {
    let lower = units.to_ascii_lowercase();
    let Some((dist, time)) = lower.split_once('/') else {
        return 1.0;
    };
    let time_scale = match time {
        "picosecond" | "picoseconds" | "ps" | "" => 1.0,
        "femtosecond" | "femtoseconds" | "fs" => 1e3,
        "nanosecond" | "nanoseconds" | "ns" => 1e-3,
        "microsecond" | "microseconds" | "us" => 1e-6,
        "second" | "seconds" | "s" => 1e-12,
        _ => 1.0,
    };
    distance_scale(dist) * time_scale
}

fn get3(data: &[f64], offset: usize) -> Result<[f64; 3], NctrajError> {
    let s = data
        .get(offset..offset + 3)
        .ok_or_else(|| NctrajError::ParseError("data ends unexpectedly".to_string()))?;
    Ok([s[0], s[1], s[2]])
}

struct CellData {
    lengths: Vec<f64>,
    length_scale: f64,
    angles: Vec<f64>,
    angle_scale: f64,
}

struct VelocityData {
    data: Vec<f64>,
    scale: f64,
}

pub(crate) struct NctrajFrameSource {
    natoms: usize,
    frame_count: usize,
    coordinates: Vec<f64>,
    coord_scale: f64,
    velocities: Option<VelocityData>,
    cell: Option<CellData>,
}

impl NctrajFrameSource {
    fn frame_inner(&self, index: usize) -> Result<Frame, NctrajError> {
        let mut positions = Vec::with_capacity(self.natoms);
        for atom in 0..self.natoms {
            let [x, y, z] = get3(&self.coordinates, (index * self.natoms + atom) * 3)?;
            positions.push(Point3::new(
                x * self.coord_scale,
                y * self.coord_scale,
                z * self.coord_scale,
            ));
        }

        let velocities = match &self.velocities {
            Some(v) => {
                let mut out = Vec::with_capacity(self.natoms);
                for atom in 0..self.natoms {
                    let [x, y, z] = get3(&v.data, (index * self.natoms + atom) * 3)?;
                    out.push(Point3::new(x * v.scale, y * v.scale, z * v.scale));
                }
                Some(out)
            }
            None => None,
        };

        let cell = match &self.cell {
            Some(c) => {
                let [a, b, cc] = get3(&c.lengths, index * 3)?;
                let [alpha, beta, gamma] = get3(&c.angles, index * 3)?;
                Some(UnitCell::new(
                    a * c.length_scale,
                    b * c.length_scale,
                    cc * c.length_scale,
                    alpha * c.angle_scale,
                    beta * c.angle_scale,
                    gamma * c.angle_scale,
                ))
            }
            None => None,
        };

        Ok(Frame {
            positions,
            velocities,
            forces: None,
            time: None,
            step: None,
            cell,
        })
    }
}

impl FrameSource for NctrajFrameSource {
    fn frame_count(&self) -> usize {
        self.frame_count
    }

    fn num_atoms(&self) -> usize {
        self.natoms
    }

    fn frame(&mut self, index: usize) -> std::io::Result<Frame> {
        self.frame_inner(index).map_err(std::io::Error::other)
    }
}

fn require_floats(var: &Netcdf3Variable) -> Result<Vec<f64>, NctrajError> {
    var.data.as_f64_slice().ok_or_else(|| {
        NctrajError::ParseError(format!("'{}' must contain floating point data", var.name))
    })
}

fn build_trajectory(bytes: Vec<u8>) -> Result<Trajectory, NctrajError> {
    if bytes.is_empty() {
        return Err(NctrajError::ParseError("empty input".to_string()));
    }
    let file: Netcdf3File = crate::io::netcdf3::read(&bytes)?;

    let conventions = file.attribute("Conventions").and_then(Netcdf3Value::as_str);
    if conventions != Some("AMBER") {
        return Err(NctrajError::NotAmberConvention(format!(
            "'Conventions' attribute is {conventions:?}, expected \"AMBER\""
        )));
    }
    let version = file
        .attribute("ConventionVersion")
        .and_then(Netcdf3Value::as_str);
    if version != Some("1.0") {
        return Err(NctrajError::NotAmberConvention(format!(
            "'ConventionVersion' attribute is {version:?}, expected \"1.0\""
        )));
    }

    let atom_dim = file
        .dimension("atom")
        .ok_or_else(|| NctrajError::MissingDimension("atom".to_string()))?;
    let natoms = atom_dim.length.ok_or_else(|| {
        NctrajError::ParseError("'atom' must not be the record dimension".to_string())
    })?;

    let spatial_dim = file
        .dimension("spatial")
        .ok_or_else(|| NctrajError::MissingDimension("spatial".to_string()))?;
    if spatial_dim.length != Some(3) {
        return Err(NctrajError::ParseError(
            "'spatial' dimension must have size 3".to_string(),
        ));
    }

    let frame_dim = file
        .dimension("frame")
        .ok_or_else(|| NctrajError::MissingDimension("frame".to_string()))?;
    if frame_dim.length.is_some() {
        return Err(NctrajError::ParseError(
            "'frame' must be the record dimension".to_string(),
        ));
    }

    let coordinates_var = file
        .variable("coordinates")
        .ok_or_else(|| NctrajError::ParseError("missing 'coordinates' variable".to_string()))?;
    let coord_scale = distance_scale(units_of(coordinates_var));
    let coordinates = require_floats(coordinates_var)?;

    let velocities = match file.variable("velocities") {
        Some(v) => Some(VelocityData {
            scale: velocity_scale(units_of(v)),
            data: require_floats(v)?,
        }),
        None => None,
    };

    let cell = match (file.variable("cell_lengths"), file.variable("cell_angles")) {
        (Some(cl), Some(ca)) => Some(CellData {
            length_scale: distance_scale(units_of(cl)),
            lengths: require_floats(cl)?,
            angle_scale: angle_scale(units_of(ca)),
            angles: require_floats(ca)?,
        }),
        (None, None) => None,
        _ => {
            return Err(NctrajError::ParseError(
                "'cell_lengths' and 'cell_angles' must be present together".to_string(),
            ));
        }
    };

    let frame_count = file.num_records();
    let source = NctrajFrameSource {
        natoms,
        frame_count,
        coordinates,
        coord_scale,
        velocities,
        cell,
    };

    let mut topology = Molecule::new();
    for _ in 0..natoms {
        topology.add_atom(Atom::new(Element::UNKNOWN));
    }
    Ok(Trajectory::new(topology, Box::new(source))?)
}

/// [`crate::io::format::ByteReadFn`] for NCTRAJ.
pub(crate) fn read_nctraj_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
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

/// [`crate::io::format::ByteWriteTrajectoryFn`] for NCTRAJ.
pub(crate) fn write_nctraj_bytes(trajectory: &mut Trajectory, _options: &WriteOptions) -> Vec<u8> {
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
    let has_velocities = frame0.as_ref().is_some_and(|f| f.velocities.is_some());
    let has_cell = frame0.as_ref().is_some_and(|f| f.cell.is_some());

    let mut coordinates = Vec::with_capacity(nframes * natoms * 3);
    let mut velocities = has_velocities.then(|| Vec::with_capacity(nframes * natoms * 3));
    let mut cell_lengths = has_cell.then(|| Vec::with_capacity(nframes * 3));
    let mut cell_angles = has_cell.then(|| Vec::with_capacity(nframes * 3));

    for i in 0..nframes {
        let frame = trajectory
            .frame(i)
            .expect("a Trajectory already validated its own frame/atom counts at construction");
        for p in &frame.positions {
            coordinates.push(p.x as f32);
            coordinates.push(p.y as f32);
            coordinates.push(p.z as f32);
        }
        if let Some(out) = &mut velocities {
            match &frame.velocities {
                Some(vs) => {
                    for v in vs {
                        out.push(v.x as f32);
                        out.push(v.y as f32);
                        out.push(v.z as f32);
                    }
                }
                None => out.extend(std::iter::repeat_n(0.0f32, natoms * 3)),
            }
        }
        if let (Some(lengths), Some(angles)) = (&mut cell_lengths, &mut cell_angles) {
            let cell = frame
                .cell
                .unwrap_or(UnitCell::new(0.0, 0.0, 0.0, 90.0, 90.0, 90.0));
            lengths.push(cell.a as f32);
            lengths.push(cell.b as f32);
            lengths.push(cell.c as f32);
            angles.push(cell.alpha as f32);
            angles.push(cell.beta as f32);
            angles.push(cell.gamma as f32);
        }
    }

    let mut dimensions = vec![
        Netcdf3Dimension {
            name: "frame".to_string(),
            length: None,
        },
        Netcdf3Dimension {
            name: "atom".to_string(),
            length: Some(natoms),
        },
        Netcdf3Dimension {
            name: "spatial".to_string(),
            length: Some(3),
        },
    ];
    let mut variables = vec![Netcdf3Variable {
        name: "spatial".to_string(),
        dimensions: vec![2],
        attributes: vec![],
        data: Netcdf3Value::Char("xyz".to_string()),
    }];

    if has_cell {
        let cell_spatial_dim = dimensions.len();
        dimensions.push(Netcdf3Dimension {
            name: "cell_spatial".to_string(),
            length: Some(3),
        });
        let cell_angular_dim = dimensions.len();
        dimensions.push(Netcdf3Dimension {
            name: "cell_angular".to_string(),
            length: Some(3),
        });
        let label_dim = dimensions.len();
        dimensions.push(Netcdf3Dimension {
            name: "label".to_string(),
            length: Some(5),
        });
        variables.push(Netcdf3Variable {
            name: "cell_spatial".to_string(),
            dimensions: vec![cell_spatial_dim],
            attributes: vec![],
            data: Netcdf3Value::Char("abc".to_string()),
        });
        variables.push(Netcdf3Variable {
            name: "cell_angular".to_string(),
            dimensions: vec![cell_angular_dim, label_dim],
            attributes: vec![],
            data: Netcdf3Value::Char("alphabeta gamma".to_string()),
        });
    }

    variables.push(Netcdf3Variable {
        name: "coordinates".to_string(),
        dimensions: vec![0, 1, 2],
        attributes: vec![(
            "units".to_string(),
            Netcdf3Value::Char("angstrom".to_string()),
        )],
        data: Netcdf3Value::Float(coordinates),
    });
    if let Some(data) = velocities {
        variables.push(Netcdf3Variable {
            name: "velocities".to_string(),
            dimensions: vec![0, 1, 2],
            attributes: vec![(
                "units".to_string(),
                Netcdf3Value::Char("angstrom/picosecond".to_string()),
            )],
            data: Netcdf3Value::Float(data),
        });
    }
    if has_cell {
        let cell_spatial_dim = dimensions
            .iter()
            .position(|d| d.name == "cell_spatial")
            .unwrap();
        let cell_angular_dim = dimensions
            .iter()
            .position(|d| d.name == "cell_angular")
            .unwrap();
        variables.push(Netcdf3Variable {
            name: "cell_lengths".to_string(),
            dimensions: vec![0, cell_spatial_dim],
            attributes: vec![(
                "units".to_string(),
                Netcdf3Value::Char("angstrom".to_string()),
            )],
            data: Netcdf3Value::Float(cell_lengths.unwrap()),
        });
        variables.push(Netcdf3Variable {
            name: "cell_angles".to_string(),
            dimensions: vec![0, cell_angular_dim],
            attributes: vec![(
                "units".to_string(),
                Netcdf3Value::Char("degree".to_string()),
            )],
            data: Netcdf3Value::Float(cell_angles.unwrap()),
        });
    }

    let file = Netcdf3File {
        dimensions,
        attributes: vec![
            (
                "Conventions".to_string(),
                Netcdf3Value::Char("AMBER".to_string()),
            ),
            (
                "ConventionVersion".to_string(),
                Netcdf3Value::Char("1.0".to_string()),
            ),
            (
                "program".to_string(),
                Netcdf3Value::Char("chem".to_string()),
            ),
            (
                "programVersion".to_string(),
                Netcdf3Value::Char(env!("CARGO_PKG_VERSION").to_string()),
            ),
        ],
        variables,
    };
    crate::io::netcdf3::write(&file)
}

/// Buffers the whole input and parses it once, mirroring
/// [`crate::io::dcd::DcdSupplier`].
pub struct NctrajSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl NctrajSupplier {
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

impl Iterator for NctrajSupplier {
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

    fn frame(
        positions: Vec<Point3>,
        velocities: Option<Vec<Point3>>,
        cell: Option<UnitCell>,
    ) -> Frame {
        Frame {
            positions,
            velocities,
            forces: None,
            time: None,
            step: None,
            cell,
        }
    }

    #[test]
    fn test_a_coordinates_only_round_trip() {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
        ];
        let mut trajectory = trajectory_from(vec![
            frame(positions.clone(), None, None),
            frame(
                positions
                    .iter()
                    .map(|p| *p + Point3::new(0.5, 0.0, 0.0))
                    .collect(),
                None,
                None,
            ),
        ]);
        let bytes = write_nctraj_bytes(&mut trajectory, &WriteOptions::default());
        let outcome = read_nctraj_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let mut back = as_trajectory(outcome);
        assert_eq!(back.frame_count(), 2);
        let f0 = back.frame(0).unwrap();
        for (a, b) in f0.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 1e-4, "{a:?} vs {b:?}");
            assert!((a.y - b.y).abs() < 1e-4, "{a:?} vs {b:?}");
            assert!((a.z - b.z).abs() < 1e-4, "{a:?} vs {b:?}");
        }
        assert!(f0.velocities.is_none());
        assert!(f0.cell.is_none());
    }

    #[test]
    fn test_a_velocities_and_cell_round_trip() {
        let cell = UnitCell::new(30.0, 25.0, 20.0, 80.0, 85.0, 95.0);
        let positions: Vec<Point3> = (0..4).map(|i| Point3::new(i as f64, 0.0, 0.0)).collect();
        let velocities: Vec<Point3> = (0..4)
            .map(|i| Point3::new(0.0, i as f64 * 0.5, 0.0))
            .collect();
        let mut trajectory = trajectory_from(vec![frame(
            positions.clone(),
            Some(velocities.clone()),
            Some(cell),
        )]);
        let bytes = write_nctraj_bytes(&mut trajectory, &WriteOptions::default());
        let mut back = as_trajectory(read_nctraj_bytes(&bytes, &ReadOptions::default()));
        let f = back.frame(0).unwrap();
        for (a, b) in f.positions.iter().zip(&positions) {
            assert!((a.x - b.x).abs() < 1e-4);
        }
        let back_v = f.velocities.expect("velocities survive");
        for (a, b) in back_v.iter().zip(&velocities) {
            assert!((a.y - b.y).abs() < 1e-4);
        }
        let back_cell = f.cell.expect("cell survives");
        assert!((back_cell.a - cell.a).abs() < 1e-3, "{}", back_cell.a);
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
    fn test_missing_conventions_attribute_is_refused_by_name() {
        let file = Netcdf3File {
            dimensions: vec![
                Netcdf3Dimension {
                    name: "frame".to_string(),
                    length: None,
                },
                Netcdf3Dimension {
                    name: "atom".to_string(),
                    length: Some(2),
                },
                Netcdf3Dimension {
                    name: "spatial".to_string(),
                    length: Some(3),
                },
            ],
            attributes: vec![],
            variables: vec![Netcdf3Variable {
                name: "coordinates".to_string(),
                dimensions: vec![0, 1, 2],
                attributes: vec![],
                data: Netcdf3Value::Float(vec![0.0; 6]),
            }],
        };
        let bytes = crate::io::netcdf3::write(&file);
        let outcome = read_nctraj_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0].error.contains("Conventions"),
            "{}",
            outcome.skipped[0].error
        );
    }

    #[test]
    fn test_only_cell_lengths_without_cell_angles_is_refused() {
        let file = Netcdf3File {
            dimensions: vec![
                Netcdf3Dimension {
                    name: "frame".to_string(),
                    length: None,
                },
                Netcdf3Dimension {
                    name: "atom".to_string(),
                    length: Some(2),
                },
                Netcdf3Dimension {
                    name: "spatial".to_string(),
                    length: Some(3),
                },
                Netcdf3Dimension {
                    name: "cell_spatial".to_string(),
                    length: Some(3),
                },
            ],
            attributes: vec![
                (
                    "Conventions".to_string(),
                    Netcdf3Value::Char("AMBER".to_string()),
                ),
                (
                    "ConventionVersion".to_string(),
                    Netcdf3Value::Char("1.0".to_string()),
                ),
            ],
            variables: vec![
                Netcdf3Variable {
                    name: "coordinates".to_string(),
                    dimensions: vec![0, 1, 2],
                    attributes: vec![],
                    data: Netcdf3Value::Float(vec![0.0; 6]),
                },
                Netcdf3Variable {
                    name: "cell_lengths".to_string(),
                    dimensions: vec![0, 3],
                    attributes: vec![],
                    data: Netcdf3Value::Float(vec![10.0, 20.0, 30.0]),
                },
            ],
        };
        let bytes = crate::io::netcdf3::write(&file);
        let outcome = read_nctraj_bytes(&bytes, &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0].error.contains("cell_lengths"),
            "{}",
            outcome.skipped[0].error
        );
    }

    #[test]
    fn test_a_nanometer_units_attribute_scales_coordinates_to_angstrom() {
        let file = Netcdf3File {
            dimensions: vec![
                Netcdf3Dimension {
                    name: "frame".to_string(),
                    length: None,
                },
                Netcdf3Dimension {
                    name: "atom".to_string(),
                    length: Some(1),
                },
                Netcdf3Dimension {
                    name: "spatial".to_string(),
                    length: Some(3),
                },
            ],
            attributes: vec![
                (
                    "Conventions".to_string(),
                    Netcdf3Value::Char("AMBER".to_string()),
                ),
                (
                    "ConventionVersion".to_string(),
                    Netcdf3Value::Char("1.0".to_string()),
                ),
            ],
            variables: vec![Netcdf3Variable {
                name: "coordinates".to_string(),
                dimensions: vec![0, 1, 2],
                attributes: vec![(
                    "units".to_string(),
                    Netcdf3Value::Char("nanometer".to_string()),
                )],
                data: Netcdf3Value::Float(vec![1.0, 2.0, 3.0]),
            }],
        };
        let bytes = crate::io::netcdf3::write(&file);
        let mut back = as_trajectory(read_nctraj_bytes(&bytes, &ReadOptions::default()));
        let f0 = back.frame(0).unwrap();
        assert!(
            (f0.positions[0].x - 10.0).abs() < 1e-3,
            "{}",
            f0.positions[0].x
        );
        assert!(
            (f0.positions[0].y - 20.0).abs() < 1e-3,
            "{}",
            f0.positions[0].y
        );
        assert!(
            (f0.positions[0].z - 30.0).abs() < 1e-3,
            "{}",
            f0.positions[0].z
        );
    }
}
