//! GROMACS run parameters (`.mdp`, #395): `key = value` lines with `;`
//! comments and no sections.
//!
//! The first format here that describes a computation rather than a system.
//! It is taken as a [`Table`] rather than a new core type: three string
//! columns, `key`, `value` and `comment`, one row per non-blank line in file
//! order. A comment on its own line is a row with only `comment` set.
//!
//! **Nothing is interpreted.** Keys keep the file's spelling;
//! [`parameter_name`] gives the identity `grompp` uses, since `ref-t`,
//! `ref_t` and `Ref-T` are one parameter to it. Values keep their text:
//! `ref_t = 300 300` is one value per coupling group and `define = -DPOSRES`
//! is compiler flags, so splitting either is the caller's call. An empty
//! value (`define =`) is an empty string, not an absent cell. Duplicate keys
//! are kept.
//!
//! **Writing** reads the same three columns back by name, so any table can be
//! written. [`dropped_columns`] names the ones it ignores, which is what
//! `chem convert` reports when a CSV with other columns becomes an MDP.

use crate::core::table::{Column, ColumnType, Table, Value};
use crate::io::csv::value_to_string;
use crate::io::errors::{MdpError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

const KEY: &str = "key";
const VALUE: &str = "value";
const COMMENT: &str = "comment";

/// The parameter a key names, as `grompp` reads it: case-insensitive, with
/// `-` and `_` interchangeable.
pub fn parameter_name(key: &str) -> String {
    key.to_lowercase().replace('-', "_")
}

/// The columns [`write_mdp`] ignores: every one not named `key`, `value` or
/// `comment`.
pub fn dropped_columns(table: &Table) -> Vec<&str> {
    table
        .columns()
        .iter()
        .map(|c| c.name.as_str())
        .filter(|n| ![KEY, VALUE, COMMENT].contains(n))
        .collect()
}

/// Parses a whole run-parameter file.
///
/// # Errors
/// [`MdpError`] for a line with no `=` or no key; one bad line fails the file.
pub fn parse_mdp(text: &str) -> Result<Table, MdpError> {
    let mut keys = Vec::new();
    let mut values = Vec::new();
    let mut comments = Vec::new();

    for (i, raw) in text.lines().enumerate() {
        let (content, comment) = match raw.split_once(';') {
            Some((content, comment)) => (content.trim(), Some(comment.trim())),
            None => (raw.trim(), None),
        };
        let comment = comment.map(|c| Value::String(c.to_string()));

        if content.is_empty() {
            if comment.is_some() {
                keys.push(None);
                values.push(None);
                comments.push(comment);
            }
            continue;
        }

        let Some((key, value)) = content.split_once('=') else {
            return Err(MdpError::MissingEquals {
                line: i + 1,
                text: raw.to_string(),
            });
        };
        let key = key.trim();
        if key.is_empty() {
            return Err(MdpError::EmptyKey { line: i + 1 });
        }
        keys.push(Some(Value::String(key.to_string())));
        values.push(Some(Value::String(value.trim().to_string())));
        comments.push(comment);
    }

    Ok(Table::new(vec![
        string_column(KEY, keys),
        string_column(VALUE, values),
        string_column(COMMENT, comments),
    ])
    .expect("the three columns grow together"))
}

fn string_column(name: &str, values: Vec<Option<Value>>) -> Column {
    let kind = if values.iter().any(Option::is_some) {
        ColumnType::String
    } else {
        ColumnType::Empty
    };
    Column {
        name: name.to_string(),
        kind,
        values,
    }
}

/// Serialises `table` as a run-parameter file from its `key`, `value` and
/// `comment` columns. A row with no key writes only its comment, or nothing.
pub fn write_mdp(table: &Table) -> String {
    let cell = |name: &str, row: usize| {
        table
            .column(name)
            .and_then(|c| c.values[row].as_ref())
            .map(value_to_string)
    };

    let mut out = String::new();
    for row in 0..table.num_rows() {
        let comment = cell(COMMENT, row);
        let line = match (cell(KEY, row), comment) {
            (Some(key), comment) => {
                let value = cell(VALUE, row).unwrap_or_default();
                let mut line = format!("{key} = {value}");
                if let Some(c) = comment {
                    line.push_str(&format!(" ; {c}"));
                }
                line.trim_end().to_string()
            }
            (None, Some(c)) => format!("; {c}").trim_end().to_string(),
            (None, None) => continue,
        };
        out.push_str(&line);
        out.push('\n');
    }
    out
}

/// [`crate::io::format::ReadFn`] for MDP.
pub(crate) fn read_mdp_with_options(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_mdp(text) {
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

/// [`crate::io::format::ByteWriteTableFn`] for MDP.
pub(crate) fn write_mdp_table_bytes(table: &Table, _options: &WriteOptions) -> Vec<u8> {
    write_mdp(table).into_bytes()
}

/// Buffers the whole input and parses it once, the same posture
/// [`crate::io::csv::CsvSupplier`] takes: one file is one record.
pub struct MdpSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl MdpSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => {
                let outcome = read_mdp_with_options(&text, options);
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

impl Iterator for MdpSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text(v: &Option<Value>) -> Option<&str> {
        match v {
            Some(Value::String(s)) => Some(s.as_str()),
            None => None,
            other => panic!("not a string cell: {other:?}"),
        }
    }

    fn rows(table: &Table) -> Vec<(Option<&str>, Option<&str>, Option<&str>)> {
        let col = |n| &table.column(n).unwrap().values;
        (0..table.num_rows())
            .map(|r| {
                (
                    text(&col(KEY)[r]),
                    text(&col(VALUE)[r]),
                    text(&col(COMMENT)[r]),
                )
            })
            .collect()
    }

    #[test]
    fn test_keys_are_kept_as_written_and_normalise_to_one_parameter() {
        let table = parse_mdp("ref-t = 300\nref_t = 310\nRef-T = 320\n").unwrap();
        let keys: Vec<_> = rows(&table).into_iter().map(|r| r.0.unwrap()).collect();
        assert_eq!(keys, ["ref-t", "ref_t", "Ref-T"]);
        assert!(keys.iter().all(|k| parameter_name(k) == "ref_t"));
    }

    #[test]
    fn test_values_are_verbatim_and_never_type_inferred() {
        let text = "\
ref_t = 300   300
define = -DPOSRES -DFLEXIBLE
include = -I../top
nsteps = 50000
define =
";
        let table = parse_mdp(text).unwrap();
        assert_eq!(table.column(VALUE).unwrap().kind, ColumnType::String);
        let values: Vec<_> = rows(&table).into_iter().map(|r| r.1).collect();
        assert_eq!(
            values,
            [
                Some("300   300"),
                Some("-DPOSRES -DFLEXIBLE"),
                Some("-I../top"),
                Some("50000"),
                Some(""),
            ]
        );
    }

    #[test]
    fn test_both_kinds_of_comment_and_blank_lines() {
        let text = "; header\n\ntitle\t= Minimization\t; Title of run\r\n;\ndt = 0.002\n";
        let table = parse_mdp(text).unwrap();
        assert_eq!(
            rows(&table),
            [
                (None, None, Some("header")),
                (Some("title"), Some("Minimization"), Some("Title of run")),
                (None, None, Some("")),
                (Some("dt"), Some("0.002"), None),
            ]
        );
    }

    #[test]
    fn test_malformed_lines_reject_the_file() {
        assert!(matches!(
            parse_mdp("dt = 0.002\nintegrator md\n"),
            Err(MdpError::MissingEquals { line: 2, .. })
        ));
        assert!(matches!(
            parse_mdp(" = 5 ; no name\n"),
            Err(MdpError::EmptyKey { line: 1 })
        ));
    }

    #[test]
    fn test_an_empty_file_is_an_empty_table_not_an_error() {
        let table = parse_mdp("").unwrap();
        assert_eq!(table.num_rows(), 0);
        assert_eq!(table.num_columns(), 3);
    }

    #[test]
    fn test_round_trip_is_exact() {
        let text = "; header\ntitle = Minimization ; Title of run\ndefine =\nref_t = 300 300\n;\n";
        let table = parse_mdp(text).unwrap();
        assert_eq!(write_mdp(&table), text);
        assert_eq!(parse_mdp(&write_mdp(&table)).unwrap(), table);
    }

    #[test]
    fn test_any_table_writes_and_names_what_it_drops() {
        let table = Table::from_csv("key,value,units\ndt,0.002,ps\nnsteps,500,\n").unwrap();
        assert_eq!(write_mdp(&table), "dt = 0.002\nnsteps = 500\n");
        assert_eq!(dropped_columns(&table), ["units"]);

        let unrelated = Table::from_csv("a,b\n1,2\n").unwrap();
        assert_eq!(write_mdp(&unrelated), "");
        assert_eq!(dropped_columns(&unrelated), ["a", "b"]);
    }
}
