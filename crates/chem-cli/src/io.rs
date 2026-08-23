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
use chemio::reader::{self, Format, ReadOutcome};
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
    let from_stdin = match path {
        None => true,
        Some(p) => p.as_os_str() == "-",
    };

    let (content, label, detected) = if from_stdin {
        let mut buf = String::new();
        std::io::stdin()
            .read_to_string(&mut buf)
            .context("reading standard input")?;
        (buf, "-".to_owned(), Format::Smiles)
    } else {
        let p = path.expect("checked above");
        let content =
            std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
        let label = p.display().to_string();
        (content, label.clone(), Format::from_filename(&label))
    };

    let format = format.unwrap_or(detected);
    Ok(Input {
        outcome: reader::read(&content, format),
        format,
        label,
    })
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
