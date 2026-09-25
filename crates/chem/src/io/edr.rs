//! GROMACS energy files (`.edr`, #399): a names block declaring each energy
//! term and its unit once, then XDR frames of per-term values.
//!
//! **Version 5 only**, which every GROMACS since 4.6 (2013) writes. Versions
//! 1-4 lay the frame header out differently -- no `dt`, no 64-bit `nsteps`,
//! an old-style restraint count -- and a reader that guessed would produce
//! plausible numbers, so they are refused by name instead.
//!
//! **Precision is not stated anywhere**; it is inferred the way GROMACS
//! itself does, from frame 0's `e_size` (`nre * 4 * sizeof(real)`), and held
//! for the whole file.
//!
//! **A [`Table`] shaped like [`crate::io::xvg`]'s**: `x` is the time and
//! each term is a `Float` column under GROMACS's own name. Metadata holds
//! `xaxis_label` (`Time (ps)`), `edr_version` and one `unit:<term>` per
//! term. Not kept: the step (t = step * dt), the running sums a frame
//! carries when `nsum > 0`, and the blocks' free-energy or restraint data --
//! each read past by its declared size. A truncated final frame, which a
//! crashed run leaves, keeps every complete frame and reports one skip, as
//! `gmx energy` does.
//!
//! **Writing** produces version 5 in single precision, what a mixed-precision
//! `mdrun` writes: the first numeric column is the time and the rest are
//! terms, with no running sums or blocks. The table carries no step, so
//! each frame's step is its index.

use crate::core::table::{Column, ColumnType, Table, Value};
use crate::io::errors::{EdrError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};
use crate::io::xdr::{XdrReader, XdrWriter};
use crate::io::xvg::is_numeric;

const NAMES_MAGIC: i32 = -55555;
const FRAME_MAGIC: i32 = -7777777;
const VERSION: i32 = 5;
/// The real every frame opens with; anything above `-1e10` is a version-1
/// frame, whose first word is the time.
const FIRST_REAL: f32 = -2e10;

#[derive(Clone, Copy, PartialEq, Debug)]
enum Precision {
    Single,
    Double,
}

impl Precision {
    fn size(self) -> usize {
        match self {
            Precision::Single => 4,
            Precision::Double => 8,
        }
    }
}

fn read_real(r: &mut XdrReader, p: Precision) -> Result<f64, EdrError> {
    Ok(match p {
        Precision::Single => f64::from(r.read_f32()?),
        Precision::Double => r.read_f64()?,
    })
}

struct Header {
    t: f64,
    nsum: i32,
    nre: usize,
    /// `(type, count)` per sub-block, in order.
    subs: Vec<(i32, usize)>,
    e_size: i32,
}

fn count(v: i32, frame: usize) -> Result<usize, EdrError> {
    usize::try_from(v).map_err(|_| EdrError::NegativeCount { frame })
}

fn read_header(r: &mut XdrReader, p: Precision, frame: usize) -> Result<Header, EdrError> {
    if read_real(r, p)? > -1e10 {
        return Err(EdrError::UnsupportedVersion(1));
    }
    let magic = r.read_i32()?;
    if magic != FRAME_MAGIC {
        return Err(EdrError::BadFrameMagic {
            frame,
            found: magic,
        });
    }
    let version = r.read_i32()?;
    if version != VERSION {
        return Err(EdrError::UnsupportedVersion(version));
    }
    let t = r.read_f64()?;
    let _step = r.read_i64()?;
    let nsum = r.read_i32()?;
    let _nsteps = r.read_i64()?;
    let _dt = r.read_f64()?;
    let nre = count(r.read_i32()?, frame)?;
    let _reserved = r.read_i32()?;
    let nblock = count(r.read_i32()?, frame)?;
    let mut subs = Vec::new();
    for _ in 0..nblock {
        let _id = r.read_i32()?;
        let nsub = count(r.read_i32()?, frame)?;
        for _ in 0..nsub {
            let kind = r.read_i32()?;
            let nr = count(r.read_i32()?, frame)?;
            subs.push((kind, nr));
        }
    }
    let e_size = r.read_i32()?;
    r.read_i32()?;
    r.read_i32()?;
    Ok(Header {
        t,
        nsum,
        nre,
        subs,
        e_size,
    })
}

/// Frame 0's header under whichever precision its `e_size` agrees with.
fn detect_precision(bytes: &[u8], start: usize) -> Result<Precision, EdrError> {
    let mut first_err = None;
    for p in [Precision::Single, Precision::Double] {
        let mut r = XdrReader::new(bytes);
        r.skip(start)?;
        match read_header(&mut r, p, 0) {
            Ok(h) if h.e_size as i64 == (h.nre * 4 * p.size()) as i64 => return Ok(p),
            Ok(h) => {
                first_err.get_or_insert(EdrError::UnknownPrecision {
                    nre: h.nre,
                    e_size: h.e_size,
                });
            }
            Err(e) => {
                first_err.get_or_insert(e);
            }
        }
    }
    Err(first_err.expect("two attempts were made"))
}

/// Skips one frame's block data; sizes per GROMACS's `XdrDataType`.
fn skip_blocks(r: &mut XdrReader, subs: &[(i32, usize)], frame: usize) -> Result<(), EdrError> {
    for &(kind, nr) in subs {
        let width = match kind {
            0 | 1 => 4, // int, float
            2 | 3 => 8, // double, int64
            _ => return Err(EdrError::UnsupportedBlockType { frame, kind }),
        };
        r.skip(width * nr)?;
    }
    Ok(())
}

/// Everything read, plus a note for a truncated final frame.
struct Parsed {
    table: Table,
    truncated: Option<String>,
}

fn parse(bytes: &[u8]) -> Result<Parsed, EdrError> {
    let mut r = XdrReader::new(bytes);
    let magic = r.read_i32()?;
    // Version 1 has no magic: its first word is the term count. Any other
    // file starting with a positive integer looks the same (XTC's 1995), so
    // the refusal names both readings rather than claiming an old EDR.
    if magic > 0 {
        return Err(EdrError::VersionOneOrNotEdr(magic));
    }
    if magic != NAMES_MAGIC {
        return Err(EdrError::NotEdr(magic));
    }
    let version = r.read_i32()?;
    if version != VERSION {
        return Err(EdrError::UnsupportedVersion(version));
    }
    let nre = count(r.read_i32()?, 0)?;
    let mut names = Vec::with_capacity(nre);
    let mut units = Vec::with_capacity(nre);
    for _ in 0..nre {
        names.push(r.read_string()?);
        units.push(r.read_string()?);
    }

    let mut times = Vec::new();
    let mut values: Vec<Vec<f64>> = vec![Vec::new(); nre];
    let mut truncated = None;
    if r.remaining() > 0 {
        let precision = detect_precision(bytes, r.position())?;
        let mut frame = 0;
        while r.remaining() > 0 {
            let start = r.position();
            match read_frame(&mut r, precision, frame, nre) {
                Ok(Some((t, e))) => {
                    times.push(t);
                    for (column, v) in values.iter_mut().zip(e) {
                        column.push(v);
                    }
                }
                Ok(None) => {}
                // Running out of bytes mid-frame is a crashed run's tail;
                // anything else is a corrupt or unknown file.
                Err(EdrError::Xdr(_)) => {
                    truncated = Some(format!(
                        "frame {frame} is incomplete ({} of its bytes present); kept the {} before it",
                        bytes.len() - start,
                        times.len()
                    ));
                    break;
                }
                Err(e) => return Err(e),
            }
            frame += 1;
        }
    }

    let float_column = |name: String, cells: Vec<f64>| Column {
        name,
        kind: if cells.is_empty() {
            ColumnType::Empty
        } else {
            ColumnType::Float
        },
        values: cells.into_iter().map(|v| Some(Value::Float(v))).collect(),
    };
    let mut columns = vec![float_column("x".to_string(), times)];
    columns.extend(
        names
            .iter()
            .cloned()
            .zip(values)
            .map(|(name, cells)| float_column(name, cells)),
    );
    let mut metadata = vec![
        ("xaxis_label".to_string(), "Time (ps)".to_string()),
        ("edr_version".to_string(), VERSION.to_string()),
    ];
    metadata.extend(
        names
            .iter()
            .zip(units)
            .map(|(name, unit)| (format!("unit:{name}"), unit)),
    );

    Ok(Parsed {
        table: Table::new(columns)
            .expect("every column has one value per frame")
            .with_metadata(metadata),
        truncated,
    })
}

/// One frame's time and term values, or `None` for a frame with no terms.
fn read_frame(
    r: &mut XdrReader,
    p: Precision,
    frame: usize,
    nre: usize,
) -> Result<Option<(f64, Vec<f64>)>, EdrError> {
    let h = read_header(r, p, frame)?;
    if h.nre != 0 && h.nre != nre {
        return Err(EdrError::TermCountMismatch {
            frame,
            expected: nre,
            found: h.nre,
        });
    }
    let mut e = Vec::with_capacity(h.nre);
    for _ in 0..h.nre {
        e.push(read_real(r, p)?);
        if h.nsum > 0 {
            // The running sum of squared deviations and the running sum.
            read_real(r, p)?;
            read_real(r, p)?;
        }
    }
    skip_blocks(r, &h.subs, frame)?;
    Ok((h.nre != 0).then_some((h.t, e)))
}

/// Parses a whole energy file.
///
/// # Errors
/// Any [`EdrError`]. A truncated final frame is not an error; use
/// `read_edr_bytes` to see it reported.
pub fn parse_edr(bytes: &[u8]) -> Result<Table, EdrError> {
    parse(bytes).map(|p| p.table)
}

/// Serialises `table` as a version-5, single-precision energy file; see the
/// module doc for which columns become what.
pub fn write_edr(table: &Table) -> Vec<u8> {
    let columns: Vec<&Column> = table.columns().iter().filter(|c| is_numeric(c)).collect();
    let (time, terms) = match columns.split_first() {
        Some((time, terms)) => (Some(*time), terms),
        None => (None, &[][..]),
    };
    let as_f64 = |v: &Option<Value>| match v {
        Some(Value::Integer(i)) => *i as f64,
        Some(Value::Float(f)) => *f,
        _ => unreachable!("only numeric columns are written"),
    };

    let mut w = XdrWriter::new();
    w.write_i32(NAMES_MAGIC);
    w.write_i32(VERSION);
    w.write_i32(terms.len() as i32);
    for term in terms {
        w.write_string(&term.name);
        let unit = table
            .metadata_value(&format!("unit:{}", term.name))
            .unwrap_or("");
        w.write_string(unit);
    }

    let Some(time) = time else {
        return w.into_bytes();
    };
    for (row, t) in time.values.iter().enumerate() {
        w.write_f32(FIRST_REAL);
        w.write_i32(FRAME_MAGIC);
        w.write_i32(VERSION);
        w.write_f64(as_f64(t));
        w.write_i64(row as i64);
        w.write_i32(0); // nsum: no running sums
        w.write_i64(0); // nsteps
        w.write_f64(0.0); // dt
        w.write_i32(terms.len() as i32);
        w.write_i32(0); // reserved
        w.write_i32(0); // nblock
        w.write_i32((terms.len() * 16) as i32); // e_size, single precision
        w.write_i32(0);
        w.write_i32(0);
        for term in terms {
            w.write_f32(as_f64(&term.values[row]) as f32);
        }
    }
    w.into_bytes()
}

/// [`crate::io::format::ByteReadFn`] for EDR.
pub(crate) fn read_edr_bytes(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse(bytes) {
        Ok(parsed) => {
            out.records.push(Record {
                payload: Payload::Table(parsed.table),
                name: "Molecule_1".to_string(),
                smiles: None,
            });
            if let Some(note) = parsed.truncated {
                out.skipped.push(Skipped {
                    position: 1,
                    input: String::new(),
                    error: note,
                });
            }
        }
        Err(e) => out.skipped.push(Skipped {
            position: 1,
            input: String::new(),
            error: e.to_string(),
        }),
    }
    out
}

/// [`crate::io::format::ByteWriteTableFn`] for EDR.
pub(crate) fn write_edr_table_bytes(table: &Table, _options: &WriteOptions) -> Vec<u8> {
    write_edr(table)
}

/// Buffers the whole input and parses it once: one file is one record.
pub struct EdrSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl EdrSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match std::io::Read::read_to_end(&mut reader, &mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => {
                let outcome = read_edr_bytes(&bytes, options);
                let mut results: Vec<Result<Record, ReadError>> =
                    outcome.records.into_iter().map(Ok).collect();
                for skipped in outcome.skipped {
                    results.push(Err(ReadError::Parse {
                        position: skipped.position,
                        message: skipped.error,
                    }));
                }
                results
            }
        };
        Self {
            records: records.into_iter(),
        }
    }
}

impl Iterator for EdrSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds an energy file by hand, independently of [`write_edr`], so a
    /// test of the reader does not just check the writer's own assumptions.
    struct Fixture {
        precision: Precision,
        names_version: i32,
        frame_version: i32,
        terms: Vec<(&'static str, &'static str)>,
    }

    struct FixtureFrame {
        t: f64,
        values: Vec<f64>,
        nsum: i32,
        /// `(type, values)` for one block of one sub-block each.
        blocks: Vec<(i32, usize)>,
        nre_override: Option<i32>,
    }

    fn frame(t: f64, values: &[f64]) -> FixtureFrame {
        FixtureFrame {
            t,
            values: values.to_vec(),
            nsum: 0,
            blocks: vec![],
            nre_override: None,
        }
    }

    impl Fixture {
        fn v5(terms: &[(&'static str, &'static str)]) -> Self {
            Self {
                precision: Precision::Single,
                names_version: 5,
                frame_version: 5,
                terms: terms.to_vec(),
            }
        }

        fn bytes(&self, frames: &[FixtureFrame]) -> Vec<u8> {
            let be = |v: i32| v.to_be_bytes().to_vec();
            let string = |s: &str| {
                let mut out = (s.len() as u32).to_be_bytes().to_vec();
                out.extend_from_slice(s.as_bytes());
                out.resize(4 + s.len().div_ceil(4) * 4, 0);
                out
            };
            let real = |v: f64| match self.precision {
                Precision::Single => (v as f32).to_be_bytes().to_vec(),
                Precision::Double => v.to_be_bytes().to_vec(),
            };
            let mut out = Vec::new();
            out.extend(be(NAMES_MAGIC));
            out.extend(be(self.names_version));
            out.extend(be(self.terms.len() as i32));
            for (name, unit) in &self.terms {
                out.extend(string(name));
                out.extend(string(unit));
            }
            for f in frames {
                let nre = f.nre_override.unwrap_or(f.values.len() as i32);
                out.extend(real(-2e10));
                out.extend(be(FRAME_MAGIC));
                out.extend(be(self.frame_version));
                out.extend(f.t.to_be_bytes());
                out.extend(7i64.to_be_bytes());
                out.extend(be(f.nsum));
                out.extend(1i64.to_be_bytes());
                out.extend(0.002f64.to_be_bytes());
                out.extend(be(nre));
                out.extend(be(0));
                out.extend(be(f.blocks.len() as i32));
                for (kind, nr) in &f.blocks {
                    out.extend(be(99)); // block id
                    out.extend(be(1)); // one sub-block
                    out.extend(be(*kind));
                    out.extend(be(*nr as i32));
                }
                out.extend(be(nre * 4 * self.precision.size() as i32));
                out.extend(be(0));
                out.extend(be(0));
                for v in &f.values {
                    out.extend(real(*v));
                    if f.nsum > 0 {
                        out.extend(real(123.0));
                        out.extend(real(456.0));
                    }
                }
                for (kind, nr) in &f.blocks {
                    let width = if *kind >= 2 { 8 } else { 4 };
                    out.extend(vec![0xAB; width * nr]);
                }
            }
            out
        }
    }

    fn floats(table: &Table, name: &str) -> Vec<f64> {
        table
            .column(name)
            .unwrap_or_else(|| panic!("no column {name}"))
            .values
            .iter()
            .map(|v| match v {
                Some(Value::Float(f)) => *f,
                other => panic!("not a float: {other:?}"),
            })
            .collect()
    }

    const TERMS: &[(&str, &str)] = &[("Potential", "kJ/mol"), ("Temperature", "K")];

    #[test]
    fn test_reads_terms_units_and_time() {
        let bytes =
            Fixture::v5(TERMS).bytes(&[frame(0.0, &[-1.5, 300.0]), frame(10.0, &[-2.5, 301.0])]);
        let table = parse_edr(&bytes).unwrap();
        let names: Vec<_> = table.columns().iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["x", "Potential", "Temperature"]);
        assert_eq!(floats(&table, "x"), [0.0, 10.0]);
        assert_eq!(floats(&table, "Potential"), [-1.5, -2.5]);
        assert_eq!(table.metadata_value("unit:Temperature"), Some("K"));
        assert_eq!(table.metadata_value("xaxis_label"), Some("Time (ps)"));
        assert_eq!(table.metadata_value("edr_version"), Some("5"));
    }

    #[test]
    fn test_double_precision_is_detected_from_the_energy_size() {
        let mut fixture = Fixture::v5(TERMS);
        fixture.precision = Precision::Double;
        let bytes = fixture.bytes(&[frame(0.0, &[-1.000000001, 300.0])]);
        let table = parse_edr(&bytes).unwrap();
        // A value only a double holds: the float reading would have rounded it.
        assert_eq!(floats(&table, "Potential"), [-1.000000001]);
    }

    #[test]
    fn test_running_sums_and_blocks_are_read_past() {
        let mut f = frame(5.0, &[-1.5, 300.0]);
        f.nsum = 50;
        f.blocks = vec![(1, 3), (2, 2), (0, 1), (3, 1)];
        let bytes = Fixture::v5(TERMS).bytes(&[f, frame(6.0, &[-2.5, 301.0])]);
        let table = parse_edr(&bytes).unwrap();
        assert_eq!(floats(&table, "Potential"), [-1.5, -2.5]);
        assert_eq!(floats(&table, "x"), [5.0, 6.0]);
    }

    #[test]
    fn test_old_versions_and_foreign_files_are_refused() {
        let mut v4 = Fixture::v5(TERMS);
        v4.names_version = 4;
        assert!(matches!(
            parse_edr(&v4.bytes(&[])),
            Err(EdrError::UnsupportedVersion(4))
        ));

        let mut frame_v4 = Fixture::v5(TERMS);
        frame_v4.frame_version = 4;
        assert!(matches!(
            parse_edr(&frame_v4.bytes(&[frame(0.0, &[1.0, 2.0])])),
            Err(EdrError::UnsupportedVersion(4))
        ));

        // Version 1 has no names magic: its first word is the term count,
        // which an XTC's own magic (1995) is indistinguishable from.
        assert!(matches!(
            parse_edr(&2i32.to_be_bytes()),
            Err(EdrError::VersionOneOrNotEdr(2))
        ));
        let xtc = parse_edr(&1995i32.to_be_bytes()).unwrap_err().to_string();
        assert!(xtc.contains("not a GROMACS energy file"), "{xtc}");
        assert!(xtc.contains("version 1"), "{xtc}");
        assert!(matches!(
            parse_edr(&(-1i32).to_be_bytes()),
            Err(EdrError::NotEdr(-1))
        ));
    }

    #[test]
    fn test_unknown_block_types_and_term_count_mismatches_are_refused() {
        let mut f = frame(0.0, &[1.0, 2.0]);
        f.blocks = vec![(5, 1)];
        assert!(matches!(
            parse_edr(&Fixture::v5(TERMS).bytes(&[f])),
            Err(EdrError::UnsupportedBlockType { frame: 0, kind: 5 })
        ));

        let bytes = Fixture::v5(TERMS).bytes(&[frame(0.0, &[1.0, 2.0]), frame(1.0, &[1.0])]);
        assert!(matches!(
            parse_edr(&bytes),
            Err(EdrError::TermCountMismatch {
                frame: 1,
                expected: 2,
                found: 1
            })
        ));
    }

    #[test]
    fn test_a_truncated_final_frame_keeps_the_complete_ones() {
        let bytes = Fixture::v5(TERMS).bytes(&[frame(0.0, &[1.0, 2.0]), frame(1.0, &[3.0, 4.0])]);
        let cut = &bytes[..bytes.len() - 3];
        let outcome = read_edr_bytes(cut, &ReadOptions::default());
        let table = outcome.records[0].table().unwrap();
        assert_eq!(floats(table, "x"), [0.0]);
        assert_eq!(outcome.skipped.len(), 1);
        assert!(
            outcome.skipped[0].error.contains("frame 1 is incomplete"),
            "{}",
            outcome.skipped[0].error
        );
    }

    #[test]
    fn test_writer_round_trips_and_names_nothing_it_did_not_have() {
        let bytes =
            Fixture::v5(TERMS).bytes(&[frame(0.0, &[-1.5, 300.0]), frame(10.0, &[-2.5, 301.0])]);
        let table = parse_edr(&bytes).unwrap();
        let written = write_edr(&table);
        assert_eq!(parse_edr(&written).unwrap(), table);
        assert_eq!(write_edr(&parse_edr(&written).unwrap()), written);
    }

    #[test]
    fn test_writer_uses_the_first_numeric_column_as_time() {
        let table = Table::from_csv("t,phase,E\n0,solid,1.5\n2,liquid,2.5\n").unwrap();
        let back = parse_edr(&write_edr(&table)).unwrap();
        assert_eq!(floats(&back, "x"), [0.0, 2.0]);
        assert_eq!(floats(&back, "E"), [1.5, 2.5]);
        assert!(back.column("phase").is_none());
        assert_eq!(back.metadata_value("unit:E"), Some(""));
    }

    #[test]
    fn test_real_gromacs_output_reads_as_gmx_energy_prints_it() {
        // The first three frames of a real energy minimisation, cut at a
        // frame boundary; the expected values are what `gmx energy` 2022.2
        // prints for the same bytes.
        let table = parse_edr(include_bytes!(
            "../../tests/corpus/edr/em-first-3-frames.edr"
        ))
        .unwrap();
        assert_eq!(table.num_columns(), 33);
        assert_eq!(floats(&table, "x"), [0.0, 1.0, 2.0]);
        assert_eq!(
            floats(&table, "Potential"),
            [-3412099.75, -3556133.0, -3727251.75]
        );
        assert_eq!(table.metadata_value("unit:Potential"), Some("kJ/mol"));
    }
}
