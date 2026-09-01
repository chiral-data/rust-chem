//! Writing molecules back out.
//!
//! `chem fp` and `chem search` produce reports. `chem aromatic` and `chem
//! coords` produce *molecules*, which raises a question a report never does:
//! what format, and does that format survive the operation.

use anyhow::{Result, bail};
use chem::core::molecule::Molecule;
use chem::io::reader::Format;
use clap::ValueEnum;
use std::path::Path;

/// Where a command's molecules should be written, and in what format.
///
/// Output format is not a preference for these two commands, it is a
/// correctness question. `chem coords` computes geometry, and **SMILES cannot
/// carry coordinates** — writing SMILES would silently discard the entire
/// result. So the format follows what the data needs, and the reason is
/// reported rather than left for the user to notice from an unexpected file.
/// `chem::io::reader::Format` is a registry handle and not ours to derive
/// `ValueEnum` on, so this is the clap-facing mirror — the same arrangement
/// `FormatArg` uses for input. Kept adjacent to its conversion so the two
/// cannot drift.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OutputFormat {
    Smiles,
    Sdf,
}

impl From<OutputFormat> for Format {
    fn from(value: OutputFormat) -> Self {
        match value {
            OutputFormat::Smiles => Format::SMILES,
            OutputFormat::Sdf => Format::SDF,
        }
    }
}

impl OutputFormat {
    /// Picks a format, given what the caller asked for and what the data needs.
    ///
    /// `needs_coords` commands default to SDF whatever the input was, because
    /// the alternative is producing a file that is missing the thing that was
    /// just computed. An explicit `--out-format smiles` is still honoured — the
    /// user may want connectivity only — but it warns, because silently
    /// dropping a result is the failure this whole type exists to prevent.
    pub fn resolve(requested: Option<Self>, needs_coords: bool, path: Option<&Path>) -> Self {
        if let Some(explicit) = requested {
            if needs_coords && explicit == OutputFormat::Smiles {
                eprintln!(
                    "warning: --out-format smiles discards the coordinates just computed; \
                     SMILES cannot store them"
                );
            }
            return explicit;
        }

        // A named file's extension is a clear request, so honour it — with the
        // same warning if it throws the result away.
        if let Some(p) = path.filter(|p| p.as_os_str() != "-") {
            // Through the registry rather than a second copy of the
            // `.sdf`-or-else rule, which is what this used to be.
            let from_name = if Format::from_filename(&p.to_string_lossy()) == Format::SDF {
                OutputFormat::Sdf
            } else {
                OutputFormat::Smiles
            };
            if needs_coords && from_name == OutputFormat::Smiles {
                eprintln!(
                    "warning: {} is not .sdf, so the coordinates just computed will be \
                     discarded; pass --out-format sdf or use a .sdf name",
                    p.display()
                );
            }
            return from_name;
        }

        if needs_coords {
            OutputFormat::Sdf
        } else {
            OutputFormat::Smiles
        }
    }

    pub fn label(&self) -> &'static str {
        Format::from(*self).label()
    }
}

/// Serialises molecules, carrying their names.
///
/// The serialisation itself lives on the format descriptor now, so this is a
/// lookup rather than a `match` that grows a arm per format. Every registered
/// format writes today; `expect` documents that rather than hiding a `None`
/// the caller cannot act on.
pub fn render(format: OutputFormat, records: &[(String, Molecule)]) -> String {
    let format = Format::from(format);
    format
        .write(records)
        .unwrap_or_else(|| panic!("{} has no writer", format.name()))
}

/// Refuses to write over the file being read.
///
/// `chem coords mols.smi -o mols.smi` is a plausible mistake, and a bad one:
/// the input is read fully before the write begins today, so it would happen to
/// work — which is worse than failing, because it would keep working until the
/// day the reading became streamed. Requiring `--force` makes the intent
/// explicit rather than accidental.
pub fn refuse_to_clobber_input(
    input: Option<&Path>,
    output: Option<&Path>,
    force: bool,
) -> Result<()> {
    if force {
        return Ok(());
    }
    let (Some(inp), Some(out)) = (input, output) else {
        return Ok(());
    };
    if inp.as_os_str() == "-" || out.as_os_str() == "-" {
        return Ok(());
    }
    // Compared after canonicalising, so `./mols.smi` and `mols.smi` are caught.
    // A path that does not exist yet cannot be the input, so a failure to
    // canonicalise the output is not a problem.
    let same = match (inp.canonicalize(), out.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => inp == out,
    };
    if same {
        bail!(
            "{} is the input file; pass --force to overwrite it",
            out.display()
        );
    }
    Ok(())
}
