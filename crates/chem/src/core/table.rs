//! A column store: named columns, a type per column inferred on read, and
//! rows (#314). CSV is on the format list (#337) and has nowhere else to
//! land — this is the container, and the parse it needs.

use thiserror::Error;

/// What a column's non-empty cells collectively are.
///
/// `Empty` only when every cell in the column is empty — there is nothing
/// to infer a real type from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnType {
    Integer,
    Float,
    String,
    Boolean,
    Empty,
}

/// One cell's value, once inferred.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Integer(i64),
    Float(f64),
    String(String),
    Boolean(bool),
}

/// One column: a name, the type its cells were inferred to share, and the
/// cells themselves.
///
/// Plain and fully public — unlike [`Table`], nothing about one column on
/// its own can be internally inconsistent; the only invariant is between
/// columns (equal row counts), which [`Table::new`] enforces.
#[derive(Debug, Clone, PartialEq)]
pub struct Column {
    pub name: String,
    pub kind: ColumnType,
    /// One entry per row. `None` is an empty cell — never coerced to a
    /// zero, an empty string, or `false`; a cell with nothing in it is a
    /// genuinely absent value, not a default.
    pub values: Vec<Option<Value>>,
}

/// Errors constructing a [`Table`].
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum TableError {
    #[error("column {name:?} has {got} values, but the table has {expected} rows")]
    ColumnLengthMismatch {
        name: String,
        expected: usize,
        got: usize,
    },

    #[error("row {row} has {got} fields, but the header has {expected}")]
    RowLengthMismatch {
        row: usize,
        expected: usize,
        got: usize,
    },

    #[error("a quoted field starting around row {row} is never closed")]
    UnterminatedQuotedField { row: usize },
}

/// A column store: named columns, a type per column, and rows.
///
/// # Type inference
///
/// A column's type is inferred from its non-empty cells only, most specific
/// first: [`ColumnType::Boolean`] (exact, case-insensitive `true`/`false`
/// only — deliberately not `1`/`0`, which are already unambiguous integers),
/// then [`ColumnType::Integer`] (every cell parses as `i64` *and*
/// re-formatting it reproduces the original text exactly — `"01"` parses as
/// `1`, but `1.to_string()` is `"1"`, not `"01"`, so a column of `01, 02`
/// becomes [`ColumnType::String`] instead: the leading zero is significant),
/// then [`ColumnType::Float`] (every cell parses as `f64` — no equivalent
/// round-trip check; `0.5` is ordinary), then [`ColumnType::String`], the
/// always-succeeding fallback. A column whose cells are every one of them
/// empty is [`ColumnType::Empty`].
///
/// # This is not the CSV format
///
/// This type only holds the rows. Whether a given CSV reads as a set of
/// molecules or as a table, and registering [`crate::io::format::Kind::Table`]
/// against a real format, is #337's job.
///
/// # Metadata
///
/// Ordered key/value pairs describing the table as a whole rather than any
/// one column -- an XVG's title and axis labels (#398), which is where its
/// units live. Empty for every table built by [`Self::new`] or
/// [`Self::from_csv`]; CSV has nowhere to put it.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    columns: Vec<Column>,
    metadata: Vec<(String, String)>,
}

impl Table {
    /// # Errors
    /// [`TableError::ColumnLengthMismatch`] if the columns don't all have
    /// the same number of values.
    pub fn new(columns: Vec<Column>) -> Result<Self, TableError> {
        let expected = columns.first().map_or(0, |c| c.values.len());
        for column in &columns {
            if column.values.len() != expected {
                return Err(TableError::ColumnLengthMismatch {
                    name: column.name.clone(),
                    expected,
                    got: column.values.len(),
                });
            }
        }
        Ok(Self {
            columns,
            metadata: Vec::new(),
        })
    }

    /// Parses `text` as CSV: the first row is column names, every row after
    /// it is data, and each column's type is inferred from its cells (see
    /// the type's own doc comment).
    ///
    /// ```
    /// use chem::core::prelude::*;
    ///
    /// let table = Table::from_csv("name,count\n\"Smith, J.\",3\nJones,07\n").unwrap();
    ///
    /// // The quoted field's embedded comma survives; `07`'s leading zero
    /// // makes the whole column a string, not an integer.
    /// assert_eq!(table.column("name").unwrap().kind, ColumnType::String);
    /// assert_eq!(table.column("count").unwrap().kind, ColumnType::String);
    /// ```
    ///
    /// # Errors
    /// [`TableError::UnterminatedQuotedField`] if a quoted field never
    /// closes, [`TableError::RowLengthMismatch`] if a data row has a
    /// different number of fields than the header, and anything
    /// [`Self::new`] would reject (unreachable here in practice, since the
    /// columns built are already equal-length by construction).
    pub fn from_csv(text: &str) -> Result<Self, TableError> {
        let rows = tokenize_csv(text)?;
        let Some((header, data_rows)) = rows.split_first() else {
            return Ok(Self {
                columns: Vec::new(),
                metadata: Vec::new(),
            });
        };

        let mut raw_columns: Vec<Vec<String>> = vec![Vec::new(); header.len()];
        for (row_index, row) in data_rows.iter().enumerate() {
            if row.len() != header.len() {
                return Err(TableError::RowLengthMismatch {
                    row: row_index + 1,
                    expected: header.len(),
                    got: row.len(),
                });
            }
            for (column, cell) in raw_columns.iter_mut().zip(row) {
                column.push(cell.clone());
            }
        }

        let columns = header
            .iter()
            .zip(raw_columns)
            .map(|(name, cells)| infer_column(name.clone(), cells))
            .collect();
        Self::new(columns)
    }

    pub fn columns(&self) -> &[Column] {
        &self.columns
    }

    pub fn column(&self, name: &str) -> Option<&Column> {
        self.columns.iter().find(|c| c.name == name)
    }

    pub fn num_rows(&self) -> usize {
        self.columns.first().map_or(0, |c| c.values.len())
    }

    pub fn num_columns(&self) -> usize {
        self.columns.len()
    }

    /// This table with `metadata` in place of whatever it had, in order.
    pub fn with_metadata(mut self, metadata: Vec<(String, String)>) -> Self {
        self.metadata = metadata;
        self
    }

    pub fn metadata(&self) -> &[(String, String)] {
        &self.metadata
    }

    /// The first value stored under `key`, if any.
    pub fn metadata_value(&self, key: &str) -> Option<&str> {
        self.metadata
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

/// Builds one [`Column`] from its raw cell text, inferring the type — see
/// [`Table`]'s own doc comment for the precedence.
fn infer_column(name: String, raw: Vec<String>) -> Column {
    let non_empty: Vec<&str> = raw
        .iter()
        .map(String::as_str)
        .filter(|s| !s.is_empty())
        .collect();

    let kind = if non_empty.is_empty() {
        ColumnType::Empty
    } else if non_empty.iter().all(|s| is_boolean(s)) {
        ColumnType::Boolean
    } else if non_empty.iter().all(|s| is_round_trip_integer(s)) {
        ColumnType::Integer
    } else if non_empty
        .iter()
        .all(|s| s.parse::<f64>().is_ok() && !has_significant_leading_zero(s))
    {
        ColumnType::Float
    } else {
        ColumnType::String
    };

    let values = raw
        .into_iter()
        .map(|cell| {
            if cell.is_empty() {
                return None;
            }
            Some(match kind {
                ColumnType::Boolean => Value::Boolean(cell.eq_ignore_ascii_case("true")),
                ColumnType::Integer => Value::Integer(cell.parse().expect("checked above")),
                ColumnType::Float => Value::Float(cell.parse().expect("checked above")),
                ColumnType::String | ColumnType::Empty => Value::String(cell),
            })
        })
        .collect();

    Column { name, kind, values }
}

fn is_boolean(s: &str) -> bool {
    s.eq_ignore_ascii_case("true") || s.eq_ignore_ascii_case("false")
}

/// Whether `s` parses as an `i64` *and* re-formatting that integer
/// reproduces `s` exactly — the check that catches a leading zero
/// (`"01"`), which parses fine but is not what `1.to_string()` writes.
fn is_round_trip_integer(s: &str) -> bool {
    s.parse::<i64>().is_ok_and(|n| n.to_string() == s)
}

/// Whether `s`'s integer part has a leading zero followed by another digit
/// (`"01"`, `"007"`) — significant, not a formatting quirk, so this excludes
/// `s` from [`ColumnType::Float`] too, the same way [`is_round_trip_integer`]
/// already excludes it from [`ColumnType::Integer`]. Rust's own `f64`
/// parser is happy to read `"01"` as `1.0`, which is exactly the silent
/// corruption the issue this type exists for warns about. `"0.5"` is
/// unaffected: the digit after the leading `0` is `.`, not another digit.
fn has_significant_leading_zero(s: &str) -> bool {
    let digits = s.strip_prefix(['+', '-']).unwrap_or(s);
    let mut chars = digits.chars();
    matches!((chars.next(), chars.next()), (Some('0'), Some(c)) if c.is_ascii_digit())
}

/// Splits `text` into rows of raw field strings, handling quoted fields
/// (doubled `""` as one literal quote, embedded commas and newlines inside
/// quotes) and CR/CRLF/LF line endings uniformly outside them.
fn tokenize_csv(text: &str) -> Result<Vec<Vec<String>>, TableError> {
    let mut rows = Vec::new();
    let mut row = Vec::new();
    let mut field = String::new();
    let mut in_quotes = false;
    let mut chars = text.chars().peekable();

    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    field.push('"');
                    chars.next();
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(c);
            }
            continue;
        }
        match c {
            '"' => in_quotes = true,
            ',' => row.push(std::mem::take(&mut field)),
            '\r' => {
                if chars.peek() == Some(&'\n') {
                    chars.next();
                }
                row.push(std::mem::take(&mut field));
                rows.push(std::mem::take(&mut row));
            }
            '\n' => {
                row.push(std::mem::take(&mut field));
                rows.push(std::mem::take(&mut row));
            }
            _ => field.push(c),
        }
    }

    if in_quotes {
        return Err(TableError::UnterminatedQuotedField { row: rows.len() });
    }
    if !field.is_empty() || !row.is_empty() {
        row.push(field);
        rows.push(row);
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_well_formed_csv_builds_the_right_columns() {
        let table = Table::from_csv("name,age\nAlice,30\nBob,25\n").unwrap();
        assert_eq!(table.num_rows(), 2);
        assert_eq!(table.num_columns(), 2);

        let name = table.column("name").unwrap();
        assert_eq!(name.kind, ColumnType::String);
        assert_eq!(
            name.values,
            vec![
                Some(Value::String("Alice".to_string())),
                Some(Value::String("Bob".to_string())),
            ]
        );

        let age = table.column("age").unwrap();
        assert_eq!(age.kind, ColumnType::Integer);
        assert_eq!(
            age.values,
            vec![Some(Value::Integer(30)), Some(Value::Integer(25))]
        );
    }

    #[test]
    fn test_leading_zeros_infer_as_string() {
        let table = Table::from_csv("code\n01\n02\n").unwrap();
        assert_eq!(table.column("code").unwrap().kind, ColumnType::String);
    }

    #[test]
    fn test_mixed_int_and_float_infers_as_float() {
        let table = Table::from_csv("value\n1\n2\n3.5\n").unwrap();
        let column = table.column("value").unwrap();
        assert_eq!(column.kind, ColumnType::Float);
        assert_eq!(
            column.values,
            vec![
                Some(Value::Float(1.0)),
                Some(Value::Float(2.0)),
                Some(Value::Float(3.5)),
            ]
        );
    }

    #[test]
    fn test_an_empty_cell_is_not_zero() {
        let table = Table::from_csv("value\n1.0\n\n3.0\n").unwrap();
        let column = table.column("value").unwrap();
        assert_eq!(column.kind, ColumnType::Float);
        assert_eq!(
            column.values,
            vec![Some(Value::Float(1.0)), None, Some(Value::Float(3.0))]
        );
    }

    #[test]
    fn test_boolean_inference_is_exact_and_does_not_claim_zero_one() {
        let bools = Table::from_csv("flag\ntrue\nfalse\nTRUE\n").unwrap();
        assert_eq!(bools.column("flag").unwrap().kind, ColumnType::Boolean);

        let ones_and_zeros = Table::from_csv("flag\n1\n0\n").unwrap();
        assert_eq!(
            ones_and_zeros.column("flag").unwrap().kind,
            ColumnType::Integer
        );
    }

    #[test]
    fn test_an_all_empty_column_infers_as_empty() {
        let table = Table::from_csv("value\n\n\n\n").unwrap();
        assert_eq!(table.column("value").unwrap().kind, ColumnType::Empty);
        assert_eq!(
            table.column("value").unwrap().values,
            vec![None, None, None]
        );
    }

    #[test]
    fn test_quoting_handles_commas_newlines_and_doubled_quotes() {
        let table = Table::from_csv("text\n\"a, b\nc\"\"d\"\n").unwrap();
        assert_eq!(
            table.column("text").unwrap().values,
            vec![Some(Value::String("a, b\nc\"d".to_string()))]
        );
    }

    #[test]
    fn test_crlf_and_lf_line_endings_split_rows_the_same_way() {
        let lf = Table::from_csv("a,b\n1,2\n3,4\n").unwrap();
        let crlf = Table::from_csv("a,b\r\n1,2\r\n3,4\r\n").unwrap();
        assert_eq!(lf, crlf);
    }

    #[test]
    fn test_row_length_mismatch_is_refused() {
        match Table::from_csv("a,b\n1,2\n3\n") {
            Err(TableError::RowLengthMismatch { row, expected, got }) => {
                assert_eq!(row, 2);
                assert_eq!(expected, 2);
                assert_eq!(got, 1);
            }
            other => panic!("expected RowLengthMismatch, got {other:?}"),
        }
    }

    #[test]
    fn test_unterminated_quoted_field_is_refused() {
        match Table::from_csv("a\n\"unterminated\n") {
            Err(TableError::UnterminatedQuotedField { .. }) => {}
            other => panic!("expected UnterminatedQuotedField, got {other:?}"),
        }
    }

    #[test]
    fn test_metadata_defaults_to_empty_and_keeps_order() {
        let table = Table::from_csv("a\n1\n").unwrap();
        assert!(table.metadata().is_empty());

        let table = table.with_metadata(vec![
            ("title".to_string(), "Energies".to_string()),
            ("xaxis_label".to_string(), "Time (ps)".to_string()),
        ]);
        assert_eq!(table.metadata()[0].0, "title");
        assert_eq!(table.metadata_value("xaxis_label"), Some("Time (ps)"));
        assert_eq!(table.metadata_value("yaxis_label"), None);
    }
}
