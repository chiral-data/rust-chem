//! Grace data sets (`.xvg`, #398), as `gmx energy` and the other GROMACS
//! analysis tools write them: `#` comments, `@` directives, then
//! whitespace-separated numeric columns.
//!
//! **A [`Table`] with metadata.** Column 0 is `x`; each later column takes
//! its `@ sN legend`, or `yN` where the file gives none. Every value is a
//! `Float`, `nan` and `inf` included, and never inferred as an integer.
//! The per-file `title`, `subtitle`, `xaxis_label` and `yaxis_label` -- the
//! axis labels carry the units -- and `@TYPE`, as `type`, go in
//! [`Table::metadata`]. Legends and labels keep Grace's escape codes
//! (`Rg\sX\N`) as written.
//!
//! **Directives count wherever they appear**, after data rows included: a
//! reader that stopped at the first number would miss a late legend.
//!
//! **One data set per file.** A trailing `&` is accepted; data after one is
//! refused, since several sets are not one wide table.
//!
//! **Writing** keeps only the columns whose every cell is a number;
//! [`dropped_columns`] names the rest, which `chem convert` reports.

use crate::core::table::{Column, ColumnType, Table, Value};
use crate::io::csv::value_to_string;
use crate::io::errors::{ReadError, XvgError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

/// Metadata keys with the directive each is written back as.
const LABELS: &[(&str, &str)] = &[
    ("title", "@    title"),
    ("subtitle", "@    subtitle"),
    ("xaxis_label", "@    xaxis  label"),
    ("yaxis_label", "@    yaxis  label"),
];
const TYPE: &str = "type";

/// The names after `x` that an error-bar `@TYPE` fixes, or `None` for the
/// types whose later columns are all y series.
fn type_columns(kind: &str) -> Result<Option<&'static [&'static str]>, XvgError> {
    match kind {
        "xy" | "nxy" | "bar" => Ok(None),
        "xydy" => Ok(Some(&["y", "dy"])),
        "xydx" => Ok(Some(&["y", "dx"])),
        "xydxdy" => Ok(Some(&["y", "dx", "dy"])),
        "xydydy" => Ok(Some(&["y", "dy1", "dy2"])),
        other => Err(XvgError::UnsupportedType(other.to_string())),
    }
}

/// The text between the first and last `"`, or the whole trimmed rest.
fn quoted(rest: &str) -> String {
    match (rest.find('"'), rest.rfind('"')) {
        (Some(a), Some(b)) if b > a => rest[a + 1..b].to_string(),
        _ => rest.trim().to_string(),
    }
}

/// The columns [`write_xvg`] skips: every one with a cell that is not a
/// number.
pub fn dropped_columns(table: &Table) -> Vec<&str> {
    table
        .columns()
        .iter()
        .filter(|c| !is_numeric(c))
        .map(|c| c.name.as_str())
        .collect()
}

fn is_numeric(column: &Column) -> bool {
    column
        .values
        .iter()
        .all(|v| matches!(v, Some(Value::Integer(_) | Value::Float(_))))
}

/// Sets `key`, replacing an earlier value: a repeated directive wins, as in
/// Grace itself.
fn set(metadata: &mut Vec<(String, String)>, key: &str, value: String) {
    match metadata.iter_mut().find(|(k, _)| k == key) {
        Some(entry) => entry.1 = value,
        None => metadata.push((key.to_string(), value)),
    }
}

/// Parses a whole XVG file.
///
/// # Errors
/// Any [`XvgError`]; one bad line fails the file.
pub fn parse_xvg(text: &str) -> Result<Table, XvgError> {
    let mut metadata: Vec<(String, String)> = Vec::new();
    let mut legends: Vec<(usize, String)> = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut ended = false;

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with('&') {
            ended = true;
            continue;
        }
        if let Some(directive) = line.strip_prefix('@') {
            let directive = directive.trim_start();
            let mut words = directive.split_whitespace();
            match words.next() {
                Some(w) if w.eq_ignore_ascii_case("type") => {
                    let kind = words.next().unwrap_or("").to_lowercase();
                    set(&mut metadata, TYPE, kind);
                }
                Some("title") => set(&mut metadata, "title", quoted(directive)),
                Some("subtitle") => set(&mut metadata, "subtitle", quoted(directive)),
                Some(axis @ ("xaxis" | "yaxis")) if words.next() == Some("label") => {
                    set(&mut metadata, &format!("{axis}_label"), quoted(directive));
                }
                Some(series) if series.starts_with('s') && words.next() == Some("legend") => {
                    if let Ok(n) = series[1..].parse::<usize>() {
                        let name = quoted(directive);
                        match legends.iter_mut().find(|(k, _)| *k == n) {
                            Some(entry) => entry.1 = name,
                            None => legends.push((n, name)),
                        }
                    }
                }
                _ => {}
            }
            continue;
        }

        if ended {
            return Err(XvgError::SecondDataSet { line: line_no });
        }
        let row = line
            .split_whitespace()
            .map(|t| {
                t.parse::<f64>().map_err(|_| XvgError::InvalidNumber {
                    line: line_no,
                    token: t.to_string(),
                })
            })
            .collect::<Result<Vec<f64>, _>>()?;
        if let Some(first) = rows.first()
            && row.len() != first.len()
        {
            return Err(XvgError::RaggedRow {
                line: line_no,
                expected: first.len(),
                got: row.len(),
            });
        }
        rows.push(row);
    }

    let kind = metadata
        .iter()
        .find(|(k, _)| k == TYPE)
        .map_or("xy".to_string(), |(_, v)| v.clone());
    let fixed = type_columns(&kind)?;
    let legend = |n: usize| {
        legends
            .iter()
            .find(|(k, _)| *k == n)
            .map(|(_, v)| v.clone())
    };

    let width = match rows.first() {
        Some(row) => row.len(),
        None => match fixed {
            Some(names) => 1 + names.len(),
            None => legends.iter().map(|(n, _)| n + 2).max().unwrap_or(0),
        },
    };
    if let Some(names) = fixed
        && width != 1 + names.len()
    {
        return Err(XvgError::TypeWidthMismatch {
            kind,
            expected: 1 + names.len(),
            got: width,
        });
    }

    let columns = (0..width)
        .map(|c| {
            let name = match (c, fixed) {
                (0, _) => "x".to_string(),
                (1, Some(_)) => legend(0).unwrap_or_else(|| "y".to_string()),
                (c, Some(names)) => names[c - 1].to_string(),
                (c, None) => legend(c - 1).unwrap_or_else(|| format!("y{}", c - 1)),
            };
            Column {
                name,
                kind: if rows.is_empty() {
                    ColumnType::Empty
                } else {
                    ColumnType::Float
                },
                values: rows.iter().map(|r| Some(Value::Float(r[c]))).collect(),
            }
        })
        .collect();

    Ok(Table::new(columns)
        .expect("every column has one value per row")
        .with_metadata(metadata))
}

/// Serialises `table`'s numeric columns as an XVG file: its known metadata
/// in stored order, then one legend per series, then the rows.
///
/// A column named exactly what the reader would call it with no legend
/// (`y0`, `y1`, …) gets none, so no label the source never had is invented.
pub fn write_xvg(table: &Table) -> String {
    let mut out = String::new();
    for (key, value) in table.metadata() {
        if let Some((_, directive)) = LABELS.iter().find(|(k, _)| k == key) {
            out.push_str(&format!("{directive} \"{value}\"\n"));
        }
    }
    let kind = table.metadata_value(TYPE);
    if let Some(kind) = kind {
        out.push_str(&format!("@TYPE {kind}\n"));
    }

    let columns: Vec<&Column> = table.columns().iter().filter(|c| is_numeric(c)).collect();
    let error_bars = kind.is_some_and(|k| matches!(type_columns(k), Ok(Some(_))));
    for (c, column) in columns.iter().enumerate().skip(1) {
        if error_bars && c > 1 {
            break;
        }
        let fallback = if error_bars {
            "y".to_string()
        } else {
            format!("y{}", c - 1)
        };
        if column.name != fallback {
            out.push_str(&format!("@ s{} legend \"{}\"\n", c - 1, column.name));
        }
    }

    for row in 0..table.num_rows() {
        let cells: Vec<String> = columns
            .iter()
            .map(|c| {
                c.values[row]
                    .as_ref()
                    .map(value_to_string)
                    .unwrap_or_default()
            })
            .collect();
        out.push_str(&cells.join(" "));
        out.push('\n');
    }
    out
}

/// [`crate::io::format::ReadFn`] for XVG.
pub(crate) fn read_xvg_with_options(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_xvg(text) {
        Ok(table) => out.records.push(Record {
            payload: Payload::Table(table),
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

/// [`crate::io::format::ByteWriteTableFn`] for XVG.
pub(crate) fn write_xvg_table_bytes(table: &Table, _options: &WriteOptions) -> Vec<u8> {
    write_xvg(table).into_bytes()
}

/// Buffers the whole input and parses it once, the same posture
/// [`crate::io::csv::CsvSupplier`] takes: one file is one record.
pub struct XvgSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl XvgSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => {
                let outcome = read_xvg_with_options(&text, options);
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

impl Iterator for XvgSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GYRATE: &str = r#"# This file was created by gmx gyrate
@    title "Radius of gyration (total and around axes)"
@    xaxis  label "Time (ps)"
@    yaxis  label "Rg (nm)"
@TYPE xy
@ view 0.15, 0.15, 0.75, 0.85
@ legend on
@ s0 legend "Rg"
@ s1 legend "Rg\sX\N"
         0     4.30449     3.41273
        10     4.30394      3.4171
"#;

    fn floats(table: &Table, name: &str) -> Vec<f64> {
        table
            .column(name)
            .unwrap()
            .values
            .iter()
            .map(|v| match v {
                Some(Value::Float(f)) => *f,
                other => panic!("not a float: {other:?}"),
            })
            .collect()
    }

    fn names(table: &Table) -> Vec<&str> {
        table.columns().iter().map(|c| c.name.as_str()).collect()
    }

    #[test]
    fn test_gmx_output_reads_as_columns_and_metadata() {
        let table = parse_xvg(GYRATE).unwrap();
        assert_eq!(names(&table), ["x", "Rg", r"Rg\sX\N"]);
        assert_eq!(floats(&table, "x"), [0.0, 10.0]);
        assert_eq!(table.column("x").unwrap().kind, ColumnType::Float);
        assert_eq!(floats(&table, "Rg"), [4.30449, 4.30394]);
        assert_eq!(table.metadata_value("xaxis_label"), Some("Time (ps)"));
        assert_eq!(table.metadata_value("yaxis_label"), Some("Rg (nm)"));
        assert_eq!(table.metadata_value("type"), Some("xy"));
    }

    #[test]
    fn test_a_legend_after_the_data_still_names_its_column() {
        let table = parse_xvg("0 1 2\n@ s1 legend \"late\"\n1 2 3\n").unwrap();
        assert_eq!(names(&table), ["x", "y0", "late"]);
    }

    #[test]
    fn test_nan_and_inf_are_values() {
        let table = parse_xvg("0 nan\n1 inf\n2 -inf\n").unwrap();
        let y = floats(&table, "y0");
        assert!(y[0].is_nan());
        assert_eq!(y[1..], [f64::INFINITY, f64::NEG_INFINITY]);
    }

    #[test]
    fn test_error_bar_types_name_their_columns() {
        let table = parse_xvg("@TYPE xydy\n@ s0 legend \"RMSD\"\n0 0.1 0.01\n").unwrap();
        assert_eq!(names(&table), ["x", "RMSD", "dy"]);
        assert!(matches!(
            parse_xvg("@TYPE xydy\n0 1\n"),
            Err(XvgError::TypeWidthMismatch { got: 2, .. })
        ));
        assert!(matches!(
            parse_xvg("@TYPE xyz\n0 1\n"),
            Err(XvgError::UnsupportedType(t)) if t == "xyz"
        ));
    }

    #[test]
    fn test_malformed_data_rejects_the_file() {
        assert!(matches!(
            parse_xvg("0 1\n1 2 3\n"),
            Err(XvgError::RaggedRow {
                line: 2,
                expected: 2,
                got: 3
            })
        ));
        assert!(matches!(
            parse_xvg("0 1\n1 two\n"),
            Err(XvgError::InvalidNumber { line: 2, token }) if token == "two"
        ));
    }

    #[test]
    fn test_a_trailing_ampersand_ends_the_set_and_a_second_set_is_refused() {
        assert_eq!(parse_xvg("0 1\n&\n").unwrap().num_rows(), 1);
        assert!(matches!(
            parse_xvg("0 1\n&\n0 2\n"),
            Err(XvgError::SecondDataSet { line: 3 })
        ));
    }

    #[test]
    fn test_round_trip_keeps_metadata_legends_and_values() {
        let table = parse_xvg(GYRATE).unwrap();
        let written = write_xvg(&table);
        assert_eq!(parse_xvg(&written).unwrap(), table);
        assert_eq!(write_xvg(&parse_xvg(&written).unwrap()), written);
        assert!(written.starts_with("@    title \"Radius of gyration"));
    }

    #[test]
    fn test_writer_keeps_metadata_order_and_invents_no_legend() {
        let text = "@    title \"T\"\n@    xaxis  label \"X\"\n@    subtitle \"S\"\n0 1\n";
        assert_eq!(write_xvg(&parse_xvg(text).unwrap()), text);
    }

    #[test]
    fn test_writer_keeps_only_numeric_columns_and_names_the_rest() {
        let table = Table::from_csv("t,label,E\n0,a,1.5\n1,b,\n").unwrap();
        assert_eq!(dropped_columns(&table), ["label", "E"]);
        assert_eq!(write_xvg(&table), "0\n1\n");
    }
}
