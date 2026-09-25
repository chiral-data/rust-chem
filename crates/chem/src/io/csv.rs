//! CSV (#337) -- last story of Phase 5, and the only format whose parse is
//! trivial while its semantics are not.
//!
//! **The real tension**: a [`crate::io::format::FormatDescriptor`]'s `kind`
//! is one fixed value per format, yet a CSV with a SMILES column is, to
//! every cheminformatics tool, a set of molecules, while a CSV of assay
//! measurements is a table -- same format, same parser, two different
//! record kinds.
//!
//! **Resolution**: this format is declared [`crate::io::format::Kind::Table`]
//! -- the honest, always-safe default, since a bare CSV genuinely has no
//! molecule information without a stated structure column -- with molecule
//! production as an explicit, disclosed read option
//! ([`crate::io::options::CsvReadOptions::structure_column`], naming the
//! column, default `None`). This mirrors the exact mechanism #330 already
//! established and
//! proved safe: XYZ/PDB/PDBQT/GRO are all still declared `Kind::Molecules`
//! yet produce a `Kind::Frames`-shaped [`crate::io::reader::Payload::Frames`]
//! under `MultiFrameMode::Frames`, an opt-in a caller states rather than a
//! heuristic the reader guesses. Confirmed directly against `format.rs`:
//! nothing in the registry's invariants compares a produced `Payload`
//! variant against a format's declared `Kind` -- the fidelity matrix
//! filters by declared kind alone, and every read entry point is fully
//! generic over `Payload`. So this format declaring `Kind::Table` while its
//! reader can, under an explicit option, produce
//! [`crate::io::reader::Payload::Molecule`] records is exactly as safe as
//! XYZ's own opt-in.
//!
//! This also directly answers the named hazard of column-name detection
//! (`SMILES`/`smiles`/`Smiles`/`structure`/`canonical_smiles` all being real
//! conventions) by making the question moot: the caller names the exact
//! column, so nothing is guessed.
//!
//! **Reuses [`Table::from_csv`] (#314) directly, not a second RFC4180
//! parser.** It already handles quoting, doubled quotes, embedded
//! newlines and CRLF, and already implements the type-inference discipline
//! this crate needs (`"01"` stays a string, an empty cell stays absent, not
//! zero). The only genuinely new parsing surface here is this module's own
//! first-row field counter, for the headerless case -- `Table::from_csv`
//! always treats row 0 as a header, so a headerless file gets a synthesized
//! `column_0,column_1,...` header prepended before the rest is handed to
//! `Table::from_csv` unchanged.
//!
//! A malformed or empty structure-column cell skips that one row rather
//! than failing the whole file -- directly answering the "row 400 fails
//! after 399 succeed" question the issue itself raises.

use crate::core::molecule::Molecule;
use crate::core::table::{Column, Table, Value};
use crate::io::errors::{CsvError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};
use crate::io::smiles::parse_smiles;
use crate::io::smiles_writer::write_smiles_for_molecule;

pub(crate) fn value_to_string(value: &Value) -> String {
    match value {
        Value::Integer(i) => i.to_string(),
        Value::Float(f) => f.to_string(),
        Value::String(s) => s.clone(),
        Value::Boolean(b) => b.to_string(),
    }
}

/// Quotes and doubles embedded quotes only when a field actually needs it
/// (contains a comma, quote or newline) -- the plain RFC4180 write-side rule.
fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Counts the first logical row's fields -- quote-aware (a quoted field may
/// contain a literal embedded newline before the row actually ends), used
/// only to size the synthesized header for a headerless file. Everything
/// past this one row is `Table::from_csv`'s own job.
///
/// Deliberately mirrors that function's own private tokenizer bit for bit,
/// including its one divergence from strict RFC4180: a bare `"` toggles
/// quote mode whenever it appears outside an already-open quoted region,
/// not only when it is the very first character of a fresh field (Python's
/// own `csv` module, checked directly during this story's research, only
/// treats a leading `"` as a delimiter). Matching `table.rs` exactly here
/// matters more than matching the stricter reference, since this function
/// exists only to feed `Table::from_csv` a correctly-sized synthetic
/// header -- disagreeing with it would misalign the two.
fn first_row_field_count(text: &str) -> Result<usize, CsvError> {
    let mut chars = text.chars().peekable();
    if chars.peek().is_none() {
        return Err(CsvError::ParseError("the CSV file is empty".to_string()));
    }

    let mut field_count = 1usize;
    let mut in_quotes = false;

    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                } else {
                    in_quotes = false;
                }
            }
            continue;
        }
        match c {
            '"' => in_quotes = true,
            ',' => field_count += 1,
            '\r' | '\n' => return Ok(field_count),
            _ => {}
        }
    }

    if in_quotes {
        return Err(CsvError::ParseError(
            "a quoted field in the first row is never closed".to_string(),
        ));
    }
    Ok(field_count)
}

fn build_table(
    text: &str,
    options: &crate::io::options::CsvReadOptions,
) -> Result<Table, CsvError> {
    if options.has_header {
        Ok(Table::from_csv(text)?)
    } else {
        let count = first_row_field_count(text)?;
        let header = (0..count)
            .map(|i| format!("column_{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let full_text = format!("{header}\n{text}");
        Ok(Table::from_csv(&full_text)?)
    }
}

/// [`crate::io::format::ReadFn`] for CSV.
pub(crate) fn read_csv_with_options(text: &str, options: &ReadOptions) -> ReadOutcome {
    let opts = &options.csv;
    let mut out = ReadOutcome::default();

    if opts.structure_column.is_some() && !opts.has_header {
        out.skipped.push(Skipped {
            position: 1,
            input: String::new(),
            error: "structure_column requires has_header".to_string(),
        });
        return out;
    }

    let table = match build_table(text, opts) {
        Ok(t) => t,
        Err(e) => {
            out.skipped.push(Skipped {
                position: 1,
                input: String::new(),
                error: e.to_string(),
            });
            return out;
        }
    };

    let Some(col_name) = &opts.structure_column else {
        out.records.push(Record {
            payload: Payload::Table(table),
            name: "Molecule_1".to_string(),
            smiles: None,
        });
        return out;
    };

    let Some(structure_col) = table.column(col_name) else {
        out.skipped.push(Skipped {
            position: 1,
            input: String::new(),
            error: format!("no column named {col_name:?}"),
        });
        return out;
    };
    let other_columns: Vec<&Column> = table
        .columns()
        .iter()
        .filter(|c| c.name != *col_name)
        .collect();

    for row in 0..table.num_rows() {
        let smiles_str = match &structure_col.values[row] {
            Some(v) => value_to_string(v),
            None => {
                out.skipped.push(Skipped {
                    position: row + 1,
                    input: String::new(),
                    error: "empty structure column cell".to_string(),
                });
                continue;
            }
        };
        match parse_smiles(&smiles_str) {
            Ok(mut mol) => {
                for col in &other_columns {
                    if let Some(v) = &col.values[row] {
                        mol.set_property(col.name.clone(), value_to_string(v));
                    }
                }
                out.records.push(Record {
                    payload: Payload::Molecule(mol),
                    name: format!("Molecule_{}", row + 1),
                    smiles: Some(smiles_str),
                });
            }
            Err(e) => {
                out.skipped.push(Skipped {
                    position: row + 1,
                    input: smiles_str,
                    error: e.to_string(),
                });
            }
        }
    }

    out
}

/// [`crate::io::format::ByteWriteTableFn`] for CSV -- the declared/primary
/// write shape.
pub(crate) fn write_csv_table_bytes(table: &Table, _options: &WriteOptions) -> Vec<u8> {
    let columns = table.columns();
    let mut out = columns
        .iter()
        .map(|c| csv_escape(&c.name))
        .collect::<Vec<_>>()
        .join(",");
    out.push('\n');

    for row in 0..table.num_rows() {
        let cells: Vec<String> = columns
            .iter()
            .map(|c| match &c.values[row] {
                Some(v) => csv_escape(&value_to_string(v)),
                None => String::new(),
            })
            .collect();
        out.push_str(&cells.join(","));
        out.push('\n');
    }

    out.into_bytes()
}

/// [`crate::io::format::WriteFn`] for CSV -- the opt-in molecule-list write
/// shape: a `name`/`smiles` column plus one column per distinct molecule
/// property, sorted for deterministic output (`HashMap` iteration order is
/// not this crate's to expose).
pub(crate) fn write_csv_records(records: &[(String, Molecule)], _options: &WriteOptions) -> String {
    let mut keys: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for (_, mol) in records {
        for key in mol.properties().keys() {
            keys.insert(key.clone());
        }
    }
    let keys: Vec<String> = keys.into_iter().collect();

    let mut header = vec!["name".to_string(), "smiles".to_string()];
    header.extend(keys.iter().cloned());
    let mut out = header
        .iter()
        .map(|h| csv_escape(h))
        .collect::<Vec<_>>()
        .join(",");
    out.push('\n');

    for (name, mol) in records {
        let smiles = write_smiles_for_molecule(mol);
        let mut row = vec![csv_escape(name), csv_escape(&smiles)];
        for key in &keys {
            row.push(csv_escape(mol.property(key).unwrap_or("")));
        }
        out.push_str(&row.join(","));
        out.push('\n');
    }

    out
}

/// Buffers the whole input and parses it once -- the same "small enough to
/// hold in memory" posture [`crate::io::dx::DxSupplier`] takes, generalized
/// to however many records this module's own reader actually produces (one
/// `Table`, or N `Molecule`s) rather than always exactly one.
pub struct CsvSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl CsvSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => {
                let outcome = read_csv_with_options(&text, options);
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

impl Iterator for CsvSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::options::CsvReadOptions;

    fn opts(csv: CsvReadOptions) -> ReadOptions {
        ReadOptions {
            csv,
            ..Default::default()
        }
    }

    #[test]
    fn test_a_plain_table_round_trip_with_a_header() {
        let text = "name,age\nAlice,30\nBob,25\n";
        let outcome = read_csv_with_options(text, &ReadOptions::default());
        assert!(outcome.skipped.is_empty());
        assert_eq!(outcome.records.len(), 1);
        let table = outcome.records[0].table().expect("a table payload");
        assert_eq!(table.num_rows(), 2);
        assert_eq!(table.num_columns(), 2);
    }

    #[test]
    fn test_a_headerless_file_gets_synthesized_column_names() {
        let text = "Alice,30\nBob,25\n";
        let outcome = read_csv_with_options(
            text,
            &opts(CsvReadOptions {
                has_header: false,
                structure_column: None,
            }),
        );
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        let table = outcome.records[0].table().expect("a table payload");
        assert_eq!(table.num_rows(), 2);
        assert!(table.column("column_0").is_some());
        assert!(table.column("column_1").is_some());
        // Type inference still applies past the synthesized header.
        assert_eq!(
            table.column("column_1").unwrap().kind,
            crate::core::table::ColumnType::Integer
        );
    }

    #[test]
    fn test_structure_column_produces_molecules_with_properties() {
        let text = "name,smiles,activity\nethanol,CCO,0.5\nbenzene,c1ccccc1,1.2\n";
        let outcome = read_csv_with_options(
            text,
            &opts(CsvReadOptions {
                has_header: true,
                structure_column: Some("smiles".to_string()),
            }),
        );
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        assert_eq!(outcome.records.len(), 2);
        let ethanol = outcome.records[0].molecule().expect("a molecule payload");
        assert_eq!(ethanol.num_atoms(), 3);
        assert_eq!(outcome.records[0].smiles.as_deref(), Some("CCO"));
        // The "name" column becomes a property too -- only "smiles" itself
        // is consumed as structure, every other column is a property.
        assert_eq!(ethanol.property("name"), Some("ethanol"));
        assert_eq!(ethanol.property("activity"), Some("0.5"));
    }

    #[test]
    fn test_a_malformed_smiles_row_is_skipped_not_fatal() {
        let text = "smiles\nCCO\nnot a smiles(((\nc1ccccc1\n";
        let outcome = read_csv_with_options(
            text,
            &opts(CsvReadOptions {
                has_header: true,
                structure_column: Some("smiles".to_string()),
            }),
        );
        assert_eq!(outcome.records.len(), 2, "the two good rows still parse");
        assert_eq!(outcome.skipped.len(), 1);
        assert_eq!(outcome.skipped[0].position, 2);
    }

    #[test]
    fn test_an_empty_structure_cell_is_skipped() {
        let text = "smiles\nCCO\n\nCCN\n";
        let outcome = read_csv_with_options(
            text,
            &opts(CsvReadOptions {
                has_header: true,
                structure_column: Some("smiles".to_string()),
            }),
        );
        assert_eq!(outcome.records.len(), 2);
        assert_eq!(outcome.skipped.len(), 1);
        assert_eq!(outcome.skipped[0].position, 2);
    }

    #[test]
    fn test_an_unknown_structure_column_is_a_whole_file_error() {
        let text = "name,smiles\nethanol,CCO\n";
        let outcome = read_csv_with_options(
            text,
            &opts(CsvReadOptions {
                has_header: true,
                structure_column: Some("nonexistent".to_string()),
            }),
        );
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
    }

    #[test]
    fn test_structure_column_without_a_header_is_refused() {
        let text = "CCO\n";
        let outcome = read_csv_with_options(
            text,
            &opts(CsvReadOptions {
                has_header: false,
                structure_column: Some("smiles".to_string()),
            }),
        );
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(outcome.skipped[0].error.contains("has_header"));
    }

    #[test]
    fn test_a_molecule_write_then_read_round_trip_through_the_real_entry_points() {
        use crate::io::format::Format;

        let mut ethanol = parse_smiles("CCO").unwrap();
        ethanol.set_property("activity".to_string(), "0.5".to_string());
        let mut benzene = parse_smiles("c1ccccc1").unwrap();
        benzene.set_property("weight".to_string(), "78.11".to_string());
        let records = vec![
            ("ethanol".to_string(), ethanol),
            ("benzene".to_string(), benzene),
        ];

        let text = Format::CSV
            .write(&records)
            .expect("CSV can write molecules");
        assert!(text.starts_with("name,smiles,activity,weight\n"));

        let back = crate::io::reader::read_with_options(
            &text,
            Format::CSV,
            &opts(CsvReadOptions {
                has_header: true,
                structure_column: Some("smiles".to_string()),
            }),
        );
        assert!(back.skipped.is_empty(), "{:?}", back.skipped);
        assert_eq!(back.records.len(), 2);
        let back_ethanol = back.records[0].molecule().expect("a molecule");
        assert_eq!(back_ethanol.property("activity"), Some("0.5"));
        // benzene has no "activity" property -- the union column is empty
        // for it, not defaulted to some other value.
        let back_benzene = back.records[1].molecule().expect("a molecule");
        assert_eq!(back_benzene.property("activity"), None);
        assert_eq!(back_benzene.property("weight"), Some("78.11"));
    }

    #[test]
    fn test_a_table_write_round_trips_through_table_from_csv_directly() {
        let table = Table::from_csv("name,note\n\"Smith, J.\",\"a \"\"quote\"\"\"\nJones,plain\n")
            .expect("valid CSV");

        let bytes = write_csv_table_bytes(&table, &WriteOptions::default());
        let text = String::from_utf8(bytes).unwrap();

        // Cross-checked against #314's own trusted parser, not this
        // story's own reader.
        let back = Table::from_csv(&text).expect("valid CSV written back");
        assert_eq!(back, table);
    }

    #[test]
    fn test_a_bare_mid_field_quote_opens_a_region_the_same_way_table_from_csv_does() {
        // `a,b"c,d` -- the `"` after `b` is not at the start of a field.
        // Python's own `csv` module (checked directly) treats this as a
        // literal character and would report 4 fields; `table.rs`'s own
        // tokenizer does not make that distinction and opens a quoted
        // region right there, swallowing the rest of the line (including
        // the row-terminating newline) as part of one still-open field.
        // This function must agree with `table.rs`, not with Python.
        let err = first_row_field_count("a,b\"c,d\n").unwrap_err();
        assert!(matches!(err, CsvError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_properly_opened_quoted_field_with_an_embedded_comma_counts_as_one_field() {
        let text = "\"a,b\",c\n";
        assert_eq!(first_row_field_count(text).unwrap(), 2);
    }

    #[test]
    fn test_first_row_field_count_respects_a_quoted_embedded_newline() {
        let text = "a,\"b\nstill the same field\",c\nnext row\n";
        assert_eq!(first_row_field_count(text).unwrap(), 3);
    }

    #[test]
    fn test_first_row_field_count_refuses_an_unterminated_quote() {
        let err = first_row_field_count("a,\"unterminated\n").unwrap_err();
        assert!(matches!(err, CsvError::ParseError(_)), "{err}");
    }
}
