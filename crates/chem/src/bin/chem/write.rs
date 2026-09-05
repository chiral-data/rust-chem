//! Writing molecules back out.
//!
//! `chem fp` and `chem search` produce reports. `chem aromatic` and `chem
//! coords` produce *molecules*, which raises a question a report never does:
//! what format, and does that format survive the operation.

use anyhow::{Result, bail};
use chem::core::molecule::Molecule;
use chem::io::format::{self, Carries, held};
use chem::io::reader::Format;
use std::path::Path;

/// Picks a format, given what the caller asked for and what the command
/// just produced.
///
/// Output format is not a preference for `chem aromatic`/`chem coords`, it
/// is a correctness question. `chem coords` computes geometry, and **SMILES
/// cannot carry coordinates** — writing SMILES would silently discard the
/// entire result. So the format follows what the data needs, and the reason
/// is reported rather than left for the user to notice from an unexpected
/// file.
///
/// `produced` is what the operation added to the molecules — `COORDS_2D`
/// for `chem coords`, nothing for `chem aromatic`. A command that produced
/// something defaults to a format that can hold it, because the
/// alternative is writing a file missing the thing that was just computed.
///
/// It used to be a `needs_coords: bool` that the call site passed in, with
/// `main.rs` hardcoding `true` for one command and `false` for the other —
/// a fact about the *format* living everywhere except the format table. The
/// call site now states what it produced, which is the thing it actually
/// knows, and the registry answers which formats can keep it.
///
/// No warning here any more. An explicit request is still honoured, and
/// what a write drops is reported once, by [`report_drops`], from the data
/// rather than from a guess about it.
///
/// Operates on [`Format`] directly (#215) rather than a hand-written
/// two-variant mirror, so a third registered format that can hold
/// `produced` becomes reachable here with no further edit — the old
/// version's fallback only ever recognised SDF specifically, because SDF
/// and SMILES were the only two variants that existed to recognise.
pub fn resolve_output_format(
    requested: Option<Format>,
    produced: Carries,
    path: Option<&Path>,
) -> Format {
    if let Some(explicit) = requested {
        return explicit;
    }

    // A named file's extension is a clear request, so honour it.
    if let Some(p) = path.filter(|p| p.as_os_str() != "-") {
        return Format::from_filename(&p.to_string_lossy());
    }

    // Nothing named a format, so pick the first registered one that can
    // hold what was just computed. `contains` is vacuously true for an
    // empty `produced`, so a command that produced nothing lands on the
    // first entry — SMILES — exactly as it did before.
    format::all()
        .find(|f| f.carries().contains(produced))
        .unwrap_or(Format::SMILES)
}

/// Says on stderr what a write is about to throw away.
///
/// One mechanism instead of a message per format pair. There used to be two
/// hand-written warnings — one for `--out-format smiles` after `chem coords`,
/// one for a 3D conformer flattened by the same command — and at 154 formats
/// the pairs are not enumerable. Every silent loss is a bug nobody notices,
/// which is how OpenBabel came to zero the B-factor column on every PDB write.
///
/// Accumulates one record at a time rather than over a pre-built slice, so a
/// streaming writer (`chem convert`, #214) can call [`Self::record`] as it
/// goes without materializing the whole file just to report what it drops.
/// [`report_drops`] is the non-streaming convenience over the same thing.
pub struct DropTracker {
    target: Carries,
    // Accumulated in a Vec in first-seen order rather than a map, so the
    // report reads the same way every run and two invocations can be diffed.
    losses: Vec<(&'static str, usize)>,
    per_molecule: Vec<(String, Vec<&'static str>)>,
}

impl DropTracker {
    pub fn new(target: Carries) -> Self {
        Self {
            target,
            losses: Vec::new(),
            per_molecule: Vec::new(),
        }
    }

    pub fn record(&mut self, name: &str, molecule: &Molecule) {
        let dropped = held(molecule).difference(self.target);
        if dropped.is_empty() {
            return;
        }
        let names: Vec<&'static str> = dropped.names().collect();
        for attribute in &names {
            match self.losses.iter_mut().find(|(a, _)| a == attribute) {
                Some((_, count)) => *count += 1,
                None => self.losses.push((attribute, 1)),
            }
        }
        self.per_molecule.push((name.to_string(), names));
    }

    /// Summarised by default: a per-molecule line is right for a handful and
    /// unusable for a hundred thousand, so the default names each lost
    /// attribute once with a count, and `explain` gives the breakdown.
    ///
    /// Silent when nothing is lost. The common case must not gain noise.
    pub fn report(&self, format_label: &str, explain: bool) {
        if self.losses.is_empty() {
            return;
        }

        let summary = self
            .losses
            .iter()
            .map(|(attribute, count)| format!("{attribute} ({count})"))
            .collect::<Vec<_>>()
            .join(", ");
        eprintln!(
            "{format_label} cannot carry: {summary}{}",
            if explain {
                ""
            } else {
                " \u{2014} pass --explain-drops for the molecules"
            }
        );

        if explain {
            for (name, names) in &self.per_molecule {
                eprintln!("  {name}: {}", names.join(", "));
            }
        }
    }
}

/// [`DropTracker`] over a whole slice already in memory — what every
/// non-streaming write command (`chem aromatic`, `chem coords`) uses.
pub fn report_drops(format: Format, records: &[(String, Molecule)], explain: bool) {
    let mut tracker = DropTracker::new(format.carries());
    for (name, molecule) in records {
        tracker.record(name, molecule);
    }
    tracker.report(format.label(), explain);
}

/// Serialises molecules, carrying their names.
///
/// The serialisation itself lives on the format descriptor now, so this is a
/// lookup rather than a `match` that grows a arm per format. Every registered
/// format writes today; `expect` documents that rather than hiding a `None`
/// the caller cannot act on.
pub fn render(format: Format, records: &[(String, Molecule)]) -> String {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_an_explicit_request_wins_over_everything_else() {
        let path = PathBuf::from("out.sdf");
        assert_eq!(
            resolve_output_format(Some(Format::SMILES), Carries::COORDS_2D, Some(&path)),
            Format::SMILES
        );
    }

    #[test]
    fn test_a_named_files_extension_is_honoured() {
        let path = PathBuf::from("out.sdf");
        assert_eq!(
            resolve_output_format(None, Carries::empty(), Some(&path)),
            Format::SDF
        );
    }

    #[test]
    fn test_nothing_named_falls_back_to_whatever_registered_format_holds_produced() {
        // Not a hardcoded "is it SDF" check (that only worked because SDF and
        // SMILES were the only two variants that could exist) -- whichever
        // registered format actually carries what was produced.
        assert_eq!(
            resolve_output_format(None, Carries::COORDS_2D, None),
            Format::SDF
        );
    }

    #[test]
    fn test_producing_nothing_falls_back_to_smiles() {
        assert_eq!(
            resolve_output_format(None, Carries::empty(), None),
            Format::SMILES
        );
    }
}
