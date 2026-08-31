//! Where input comes from and where output goes.
//!
//! Settled once, here, so every subcommand inherits it rather than deciding
//! again: read the named file or stdin, write the named file or stdout, and put
//! everything else on stderr.
//!
//! The reason for the last one is that it makes `chem a | chem b` possible at
//! all. A progress line on stdout is indistinguishable from data to the next
//! process in the pipe, so the rule is not a style preference — it is what
//! keeps the tool composable.

use anyhow::{Context, Result};
use chem::io::reader::{self, Format, ReadOutcome};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// What was read, and what it was called.
pub struct Input {
    pub outcome: ReadOutcome,
    pub format: Format,
    /// The filename, or `-` for stdin — for messages, not for reopening.
    pub label: String,
}

/// Reads the file at `path`, or stdin when it is `None` or `-`.
///
/// `format` overrides detection. Detection needs a filename, so stdin without
/// an override falls back to SMILES — which is the same default a file with an
/// unrecognised extension gets, rather than a special case for pipes.
pub fn read_input(path: Option<&Path>, format: Option<Format>) -> Result<Input> {
    let (content, label) = read_text(path)?;
    // No special case for stdin any more. `-` has no extension, and
    // `from_filename` already answers SMILES for a name nothing claims — the
    // branch that used to be here was restating the fallback, in a second
    // place where it could drift from it.
    let format = format.unwrap_or_else(|| Format::from_filename(&label));
    Ok(Input {
        outcome: reader::read(&content, format),
        format,
        label,
    })
}

/// Reads a named file, or stdin when the path is absent or `-`.
///
/// Every input in this tool goes through here, so `-` means the same thing
/// everywhere. `chem search` needing it is why this is a function rather than
/// inline: without it `chem fp mols.smi | chem search -` fails, and that
/// pipeline is the point of the two commands being separate.
pub fn read_text(path: Option<&Path>) -> Result<(String, String)> {
    match path {
        Some(p) if p.as_os_str() != "-" => {
            let content =
                std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            Ok((content, p.display().to_string()))
        }
        _ => {
            let mut buf = String::new();
            std::io::stdin()
                .read_to_string(&mut buf)
                .context("reading standard input")?;
            Ok((buf, "-".to_owned()))
        }
    }
}

/// Writes to the file at `path`, or stdout when it is `None` or `-`.
pub fn write_output(path: Option<&PathBuf>, contents: &str) -> Result<()> {
    match path {
        Some(p) if p.as_os_str() != "-" => {
            std::fs::write(p, contents).with_context(|| format!("writing {}", p.display()))
        }
        _ => {
            let stdout = std::io::stdout();
            let mut lock = stdout.lock();
            lock.write_all(contents.as_bytes())
                .context("writing to standard output")?;
            lock.flush().context("flushing standard output")
        }
    }
}

/// Reports what was read, on stderr.
///
/// Every subcommand says this, in the same words, so a failure count never
/// depends on which command happened to be running.
pub fn report(input: &Input) {
    eprintln!(
        "read {} molecules from {} ({})",
        input.outcome.len(),
        input.label,
        input.format.label()
    );
    for skipped in &input.outcome.skipped {
        if skipped.input.is_empty() {
            eprintln!("  skipped record {}: {}", skipped.position, skipped.error);
        } else {
            eprintln!(
                "  skipped line {} ({}): {}",
                skipped.position, skipped.input, skipped.error
            );
        }
    }
}
