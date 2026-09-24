//! GROMACS index file (`.ndx`, #394): `[ name ]` headers over
//! whitespace-separated atom indices, wrapped at no fixed width.
//!
//! **1-based in the file, 0-based in [`IndexGroups`].** The conversion
//! happens here and nowhere else; an off-by-one would not fail, it would
//! silently select the neighbouring atoms, which is why an index of `0` is
//! rejected rather than wrapped or skipped.
//!
//! **Nothing is normalised.** Duplicate group names (which `make_ndx`
//! produces), an index repeated inside a group or across groups, and empty
//! groups are all kept, in file order. `gmx` was not available to check which
//! of two same-named groups it selects by name; keeping both is the reading
//! that loses nothing, and [`IndexGroups::named`] returns them all.
//!
//! **No structure is consulted.** An index file names no coordinate file, so
//! its indices are not bounds-checked; [`IndexGroups::max_atom`] is the only
//! bound it carries.
//!
//! `;` starts a comment, as GROMACS's own line reader treats it. Any
//! malformed line rejects the whole file rather than one group: a partial
//! index is a wrong selection, not a smaller one.

use crate::core::index_groups::{IndexGroup, IndexGroups};
use crate::io::errors::{NdxError, ReadError};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

/// Indices per line on write -- what `gmx`'s own writer uses.
const INDICES_PER_LINE: usize = 15;

/// Parses a whole index file.
///
/// # Errors
/// Any [`NdxError`]; see the module doc for why one bad line fails the file.
pub fn parse_ndx(text: &str) -> Result<IndexGroups, NdxError> {
    let mut groups: Vec<IndexGroup> = Vec::new();

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.split(';').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }

        if let Some(rest) = line.strip_prefix('[') {
            let Some((name, trailing)) = rest.split_once(']') else {
                return Err(NdxError::InvalidHeader {
                    line: line_no,
                    text: line.to_string(),
                });
            };
            if !trailing.trim().is_empty() {
                return Err(NdxError::InvalidHeader {
                    line: line_no,
                    text: line.to_string(),
                });
            }
            groups.push(IndexGroup {
                name: name.trim().to_string(),
                atoms: Vec::new(),
            });
            continue;
        }

        let Some(group) = groups.last_mut() else {
            return Err(NdxError::IndexBeforeHeader { line: line_no });
        };
        for token in line.split_whitespace() {
            let index = token
                .parse::<usize>()
                .ok()
                .and_then(|n| n.checked_sub(1))
                .ok_or_else(|| NdxError::InvalidIndex {
                    line: line_no,
                    token: token.to_string(),
                })?;
            group.atoms.push(index);
        }
    }

    if groups.is_empty() {
        return Err(NdxError::NoGroups);
    }
    Ok(IndexGroups::new(groups))
}

/// Serialises `groups` as an index file, 1-based, [`INDICES_PER_LINE`] to a
/// line.
pub fn write_ndx(groups: &IndexGroups) -> String {
    let mut out = String::new();
    for group in groups.groups() {
        out.push_str(&format!("[ {} ]\n", group.name));
        for chunk in group.atoms.chunks(INDICES_PER_LINE) {
            let line: Vec<String> = chunk.iter().map(|a| format!("{:>4}", a + 1)).collect();
            out.push_str(&line.join(" "));
            out.push('\n');
        }
    }
    out
}

/// [`crate::io::format::ReadFn`] for NDX.
pub(crate) fn read_ndx_with_options(text: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match parse_ndx(text) {
        Ok(groups) => out.records.push(Record {
            payload: Payload::IndexGroups(groups),
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

/// [`crate::io::format::ByteWriteIndexGroupsFn`] for NDX.
pub(crate) fn write_ndx_bytes(groups: &IndexGroups, _options: &WriteOptions) -> Vec<u8> {
    write_ndx(groups).into_bytes()
}

/// Buffers the whole input and parses it once, the same posture
/// [`crate::io::csv::CsvSupplier`] takes: one file is one record.
pub struct NdxSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl NdxSupplier {
    pub fn new<R: std::io::BufRead>(mut reader: R, options: &ReadOptions) -> Self {
        let mut text = String::new();
        let records = match std::io::Read::read_to_string(&mut reader, &mut text) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => {
                let outcome = read_ndx_with_options(&text, options);
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

impl Iterator for NdxSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms_of(groups: &IndexGroups) -> Vec<(&str, Vec<usize>)> {
        groups
            .groups()
            .iter()
            .map(|g| (g.name.as_str(), g.atoms.clone()))
            .collect()
    }

    #[test]
    fn test_indices_are_converted_from_one_based() {
        let groups = parse_ndx("[ X ]\n1 2 3\n").unwrap();
        assert_eq!(atoms_of(&groups), vec![("X", vec![0, 1, 2])]);
    }

    #[test]
    fn test_nothing_is_normalised() {
        let text = "\
[ System ]
   1    2    3
   4    5
[ Empty ]
[ System ]
   2    2    9
";
        let groups = parse_ndx(text).unwrap();
        assert_eq!(
            atoms_of(&groups),
            vec![
                ("System", vec![0, 1, 2, 3, 4]),
                ("Empty", vec![]),
                ("System", vec![1, 1, 8]),
            ]
        );
        assert_eq!(groups.named("System").count(), 2);
        assert_eq!(groups.max_atom(), Some(8));
    }

    #[test]
    fn test_comments_blank_lines_and_names_with_spaces() {
        let text =
            "; made by hand\n\n[  Protein_&_!Water  ]  ; trailing\n10 ; eleventh? no, tenth\n";
        let groups = parse_ndx(text).unwrap();
        assert_eq!(atoms_of(&groups), vec![("Protein_&_!Water", vec![9])]);
    }

    #[test]
    fn test_zero_negative_and_non_numeric_indices_are_rejected() {
        for bad in ["0", "-3", "1.5", "abc"] {
            let err = parse_ndx(&format!("[ X ]\n1 {bad}\n")).unwrap_err();
            assert!(
                matches!(&err, NdxError::InvalidIndex { line: 2, token } if token == bad),
                "{bad}: {err}"
            );
        }
    }

    #[test]
    fn test_structural_errors() {
        assert!(matches!(
            parse_ndx("1 2\n[ X ]\n"),
            Err(NdxError::IndexBeforeHeader { line: 1 })
        ));
        assert!(matches!(
            parse_ndx("[ X\n"),
            Err(NdxError::InvalidHeader { line: 1, .. })
        ));
        assert!(matches!(
            parse_ndx("[ X ] 1 2\n"),
            Err(NdxError::InvalidHeader { line: 1, .. })
        ));
        assert!(matches!(parse_ndx(""), Err(NdxError::NoGroups)));
        assert!(matches!(
            parse_ndx("; only a comment\n"),
            Err(NdxError::NoGroups)
        ));
    }

    #[test]
    fn test_write_is_one_based_and_wraps_at_fifteen() {
        let groups = IndexGroups::new(vec![
            IndexGroup {
                name: "A".to_string(),
                atoms: (0..16).collect(),
            },
            IndexGroup {
                name: "E".to_string(),
                atoms: vec![],
            },
        ]);
        let text = write_ndx(&groups);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "[ A ]");
        assert_eq!(lines[1].split_whitespace().count(), 15);
        assert!(lines[1].starts_with("   1    2"));
        assert_eq!(lines[2].trim(), "16");
        assert_eq!(lines[3], "[ E ]");
        assert_eq!(lines.len(), 4);
    }

    #[test]
    fn test_round_trip() {
        let text = "[ A ]\n5 1 5\n[ E ]\n[ A ]\n100000\n";
        let groups = parse_ndx(text).unwrap();
        assert_eq!(parse_ndx(&write_ndx(&groups)).unwrap(), groups);
    }

    #[test]
    fn test_read_reports_a_bad_file_as_one_skip() {
        let outcome = read_ndx_with_options("[ X ]\n0\n", &ReadOptions::default());
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
        assert!(outcome.skipped[0].error.contains("\"0\""));
    }
}
