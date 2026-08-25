//! Reading whole files into molecules.
//!
//! [`crate::smiles::parse_smiles`] takes one string and
//! [`crate::sdf::parse_sdf`] takes one `$$$$`-terminated record. A
//! file is neither: it is a batch, and turning one into molecules means
//! splitting it, naming the records that carry no name, and deciding what to do
//! about the ones that will not parse.
//!
//! That logic used to live in the GUI crate, which meant nothing else could
//! read a file — the same trap depiction was in before it moved out. It is here
//! now so the application and the command line agree, by construction, about
//! how many molecules a file contains.
//!
//! # Failures are returned, not logged
//!
//! A bad record does not abort the file: a thousand-molecule dataset with one
//! malformed line is still 999 usable molecules, and refusing all of them helps
//! nobody. But the caller has to *know*, which a log line does not achieve —
//! it cannot be counted, reported in a summary, or turned into an exit code.
//! So [`ReadOutcome`] carries the skipped records alongside the good ones and
//! lets each front end decide what to do with them.

use crate::sdf::parse_sdf;
use crate::smiles::parse_smiles;
use chem_core::molecule::Molecule;

/// Which format a file is in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Smiles,
    Sdf,
}

impl Format {
    /// Picks a format from a filename.
    ///
    /// Only `.sdf` is SDF. `.smi`, `.smiles`, `.txt` and no extension at all
    /// take the SMILES path, which is deliberate rather than a gap: SMILES is
    /// the universal default and SDF is the opt-in. Anything stricter would
    /// reject the extensionless files people actually have.
    pub fn from_filename(name: &str) -> Self {
        if name.to_lowercase().ends_with(".sdf") {
            Format::Sdf
        } else {
            Format::Smiles
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Format::Smiles => "SMILES",
            Format::Sdf => "SDF",
        }
    }
}

/// One molecule read from a file.
#[derive(Debug, Clone)]
pub struct Record {
    pub molecule: Molecule,
    /// The record's own name, or `Molecule_N` where the file gave none.
    pub name: String,
    /// The SMILES the molecule was read from, where there was one.
    ///
    /// `None` for SDF, which stores coordinates and connectivity rather than a
    /// SMILES string. A placeholder belongs in whatever displays this, not
    /// here — a library should not invent text for a value it does not have.
    pub smiles: Option<String>,
}

/// A record that could not be read, and why.
#[derive(Debug, Clone)]
pub struct Skipped {
    /// Line number for SMILES, record number for SDF. One-based, so it matches
    /// what an editor or an error message would say.
    pub position: usize,
    /// The offending input, for SMILES. Empty for SDF, where the record is a
    /// multi-line block and quoting it back is noise rather than help.
    pub input: String,
    pub error: String,
}

/// Everything a file yielded, including what it did not.
///
/// Deliberately not a `Result`: reading a file cannot fail as a whole. An
/// unreadable file is the caller's problem before this is called, and every
/// per-record failure is already carried in [`Self::skipped`]. The previous
/// signature returned `anyhow::Result` and never once returned `Err`.
#[derive(Debug, Clone, Default)]
pub struct ReadOutcome {
    pub records: Vec<Record>,
    pub skipped: Vec<Skipped>,
}

impl ReadOutcome {
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }
}

/// Reads a file's contents in the given format.
pub fn read(content: &str, format: Format) -> ReadOutcome {
    match format {
        Format::Smiles => read_smiles(content),
        Format::Sdf => read_sdf(content),
    }
}

/// One molecule per line: the SMILES, then optionally a name.
///
/// Blank lines and `#` comments are skipped silently — they are not failures,
/// and counting them as such would make every commented file look broken.
pub fn read_smiles(content: &str) -> ReadOutcome {
    let mut out = ReadOutcome::default();

    for (index, raw) in content.lines().enumerate() {
        let position = index + 1;
        let line = raw.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut parts = line.split_whitespace();
        let Some(smiles) = parts.next() else {
            continue;
        };
        let rest: Vec<&str> = parts.collect();
        let name = if rest.is_empty() {
            format!("Molecule_{position}")
        } else {
            rest.join(" ")
        };

        match parse_smiles(smiles) {
            Ok(molecule) => out.records.push(Record {
                molecule,
                name,
                smiles: Some(smiles.to_owned()),
            }),
            Err(e) => out.skipped.push(Skipped {
                position,
                input: smiles.to_owned(),
                error: e.to_string(),
            }),
        }
    }

    out
}

/// One molecule per `$$$$`-terminated record.
pub fn read_sdf(content: &str) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        lines.push(line);
        if line.trim() == "$$$$" {
            position += 1;
            push_record(&mut out, &lines, position);
            lines.clear();
        }
    }

    // A trailing record with no `$$$$` — which a single-molecule file often is.
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_record(&mut out, &lines, position);
    }

    out
}

fn push_record(out: &mut ReadOutcome, lines: &[&str], position: usize) {
    let record = lines.join("\n");
    match parse_sdf(&record) {
        Ok(molecule) => {
            let name = molecule
                .name()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Molecule_{position}"));
            out.records.push(Record {
                molecule,
                name,
                smiles: None,
            });
        }
        Err(e) => out.skipped.push(Skipped {
            position,
            input: String::new(),
            error: e.to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_only_sdf_extension_means_sdf() {
        assert_eq!(Format::from_filename("a.sdf"), Format::Sdf);
        assert_eq!(Format::from_filename("A.SDF"), Format::Sdf);
        // Everything else takes the SMILES path, extensionless files included.
        for name in ["a.smi", "a.smiles", "a.txt", "a", "a.sdf.gz"] {
            assert_eq!(Format::from_filename(name), Format::Smiles, "{name}");
        }
    }

    #[test]
    fn test_blank_lines_and_comments_are_not_failures() {
        let out = read_smiles("# header\n\nCCO\n\n# trailing\n");
        assert_eq!(out.len(), 1);
        assert!(
            out.skipped.is_empty(),
            "a commented file must not look broken"
        );
    }

    #[test]
    fn test_name_is_the_rest_of_the_line() {
        let out = read_smiles("CCO ethyl alcohol\n");
        assert_eq!(out.records[0].name, "ethyl alcohol");
    }

    #[test]
    fn test_unnamed_molecules_take_their_line_number() {
        // Not their index among the *kept* records: the number has to point at
        // the line the reader actually saw, or it is useless for finding it.
        let out = read_smiles("# comment\n\nCCO\nCC\n");
        assert_eq!(out.records[0].name, "Molecule_3");
        assert_eq!(out.records[1].name, "Molecule_4");
    }

    #[test]
    fn test_a_bad_line_is_skipped_and_reported_with_its_position() {
        let out = read_smiles("CCO good\nnot-a-smiles bad\nCC also_good\n");
        assert_eq!(out.len(), 2, "the good lines survive");
        assert_eq!(out.skipped.len(), 1);
        assert_eq!(out.skipped[0].position, 2);
        assert_eq!(out.skipped[0].input, "not-a-smiles");
        assert!(!out.skipped[0].error.is_empty());
    }

    #[test]
    fn test_smiles_records_keep_their_source_and_sdf_records_do_not() {
        let smi = read_smiles("CCO\n");
        assert_eq!(smi.records[0].smiles.as_deref(), Some("CCO"));

        let sdf = read_sdf(TWO_RECORDS);
        assert_eq!(
            sdf.records[0].smiles, None,
            "SDF has no SMILES to report, and inventing one is the caller's job"
        );
    }

    const TWO_RECORDS: &str = "\
ethanol
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 C   0  0
    2.2500    1.3000    0.0000 O   0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
water
  test

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0
M  END
$$$$
";

    #[test]
    fn test_sdf_splits_on_the_record_terminator() {
        let out = read_sdf(TWO_RECORDS);
        assert_eq!(out.len(), 2);
        assert_eq!(out.records[0].name, "ethanol");
        assert_eq!(out.records[1].name, "water");
        assert!(out.skipped.is_empty());
    }

    #[test]
    fn test_a_trailing_record_without_a_terminator_still_counts() {
        // Single-molecule files often have no `$$$$` at all, and dropping the
        // last record would silently lose one molecule from every such file.
        let unterminated = TWO_RECORDS.trim_end().trim_end_matches("$$$$");
        let out = read_sdf(unterminated);
        assert_eq!(out.len(), 2);
        assert_eq!(out.records[1].name, "water");
    }

    #[test]
    fn test_reading_an_empty_file_yields_nothing_rather_than_failing() {
        for content in ["", "\n\n", "# only a comment\n"] {
            let out = read(content, Format::Smiles);
            assert!(out.is_empty());
            assert!(out.skipped.is_empty());
        }
    }

    #[test]
    fn test_read_dispatches_on_format() {
        let sdf = read(TWO_RECORDS, Format::Sdf);
        assert_eq!(sdf.len(), 2);
        assert!(sdf.records.iter().all(|r| r.molecule.num_atoms() > 0));

        // The same bytes read as SMILES yield nothing at all. Before #151 the
        // `$$$$` terminators survived as atomless molecules, so this reported
        // "read 2 molecules" and exited successfully on a wrong-format file.
        let wrong = read(TWO_RECORDS, Format::Smiles);
        assert!(wrong.is_empty(), "kept {} records", wrong.records.len());
        assert!(!wrong.skipped.is_empty());
    }

    #[test]
    fn test_a_dollar_line_is_skipped_rather_than_read_as_a_molecule() {
        // The flip side of #151. `$` is a legal SMILES token — the
        // quadruple-bond character — so `$$$$` tokenizes cleanly and used to
        // build an atomless molecule that was returned as a success. It is
        // also the SDF record terminator, so an SDF read as SMILES reported
        // one molecule per record, and a caller asking "did anything parse?"
        // was told yes.
        let out = read_smiles("$$$$\n");
        assert!(out.is_empty());
        assert_eq!(out.skipped.len(), 1);
        assert_eq!(out.skipped[0].position, 1);
    }
}
