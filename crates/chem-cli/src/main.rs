//! `chem` — cheminformatics operations on the command line.
//!
//! The operations themselves live in the libraries; this is the front end that
//! makes them usable from a shell, over files, and without a display.
//!
//! # What this crate is for
//!
//! Three things, which are one design if it is built the ordinary Unix way:
//! processing a file in batch, composing with other tools through a pipe, and
//! reaching the GPU kernels on a machine with no window server. They differ in
//! defaults, not in structure.
//!
//! Being a readable *example* of the libraries is deliberately not among them.
//! That is what `examples/` in each crate is for, where docs.rs renders it —
//! shaping this around readability would trade away the streaming, error
//! handling and exit codes a tool actually needs.

mod backend;
mod exit;
mod fpfile;
mod io;

use anyhow::{Context, Result};
use backend::Backend;
use chemio::reader::Format;
use clap::{Parser, Subcommand, ValueEnum};
use fpfile::FingerprintFile;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "chem",
    version,
    about = "Cheminformatics operations on the command line",
    long_about = None,
)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Compute backend. `gpu` fails rather than falling back.
    #[arg(long, short = 'b', global = true, value_enum, default_value_t = Backend::Auto)]
    backend: Backend,

    /// Treat skipped records as a failure, exiting 4.
    ///
    /// Off by default: one bad line in a thousand-molecule file should not throw
    /// away the other 999. On for a job that must not half-succeed.
    #[arg(long, global = true)]
    strict: bool,
}

#[derive(Subcommand)]
enum Command {
    /// Read the input and report what is in it.
    ///
    /// The operations proper are separate subcommands. This one exists because
    /// checking that a file parses, and how, is the first thing anyone does
    /// with a new dataset — and because it exercises every I/O convention the
    /// other subcommands rely on.
    Info {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        /// Input format. Detected from the filename otherwise; standard input
        /// has no name to read, so it defaults to SMILES.
        #[arg(long, value_enum)]
        format: Option<FormatArg>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },

    /// Generate Morgan fingerprints for every molecule in the input.
    Fp {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long, value_enum)]
        format: Option<FormatArg>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Morgan radius.
        #[arg(long, default_value_t = 2)]
        radius: u32,

        /// Fingerprint length in bits.
        #[arg(long, default_value_t = 2048)]
        size: u32,
    },

    /// Rank a file of fingerprints against a query molecule.
    Search {
        /// A file written by `chem fp`. Reads standard input when absent or `-`,
        /// so `chem fp mols.smi | chem search --query ...` works.
        fingerprints: Option<PathBuf>,

        /// The query, as SMILES.
        #[arg(long)]
        query: String,

        /// How many results to return. 0 means all of them.
        #[arg(long, default_value_t = 10)]
        top: usize,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },
}

/// `chemio::reader::Format` is not ours to derive `ValueEnum` on, so this is the
/// clap-facing mirror. Kept adjacent to its conversion so the two cannot drift.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum FormatArg {
    Smiles,
    Sdf,
}

impl From<FormatArg> for Format {
    fn from(value: FormatArg) -> Self {
        match value {
            FormatArg::Smiles => Format::Smiles,
            FormatArg::Sdf => Format::Sdf,
        }
    }
}

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(code) => std::process::ExitCode::from(code as u8),
        // Anyhow's chain, one line each, so a wrapped cause is not lost to a
        // single-line summary.
        Err(e) => {
            eprintln!("error: {e}");
            for cause in e.chain().skip(1) {
                eprintln!("  caused by: {cause}");
            }
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> Result<i32> {
    match &cli.command {
        Command::Info {
            input,
            format,
            output,
        } => {
            let read = io::read_input(input.as_deref(), format.map(Into::into))?;
            io::report(&read);
            // Echoed rather than resolved: `info` computes nothing, and
            // probing for a device costs a device creation. A mistyped flag is
            // still visible here rather than three commands later.
            eprintln!("backend: {} (not used by `info`)", cli.backend.label());

            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }

            io::write_output(output.as_ref(), &describe(&read))?;

            if cli.strict && !read.outcome.skipped.is_empty() {
                return Ok(exit::PARTIAL);
            }
            Ok(exit::OK)
        }

        Command::Fp {
            input,
            format,
            output,
            radius,
            size,
        } => {
            let read = io::read_input(input.as_deref(), format.map(Into::into))?;
            io::report(&read);
            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }

            let search = cli.backend.open()?;
            let molecules: Vec<_> = read
                .outcome
                .records
                .iter()
                .map(|r| r.molecule.clone())
                .collect();

            let started = Instant::now();
            let fingerprints = pollster::block_on(
                search.generate_fingerprints_batch_async(&molecules, *radius, *size),
            )
            .context("generating fingerprints")?;
            // On stderr, so a timing line never lands in the data a pipe
            // carries. The Operations window reports the same thing per run.
            eprintln!(
                "fingerprinted {} molecules in {:.0} ms (radius {radius}, {size} bits)",
                fingerprints.len(),
                started.elapsed().as_secs_f64() * 1000.0
            );

            let file = FingerprintFile {
                radius: *radius,
                size: *size,
                names: read
                    .outcome
                    .records
                    .iter()
                    .map(|r| r.name.clone())
                    .collect(),
                fingerprints,
            };
            io::write_output(output.as_ref(), &file.to_text())?;

            if cli.strict && !read.outcome.skipped.is_empty() {
                return Ok(exit::PARTIAL);
            }
            Ok(exit::OK)
        }

        Command::Search {
            fingerprints,
            query,
            top,
            output,
        } => {
            let (text, label) = io::read_text(fingerprints.as_deref())?;
            let targets =
                FingerprintFile::parse(&text).with_context(|| format!("reading {label}"))?;
            eprintln!(
                "read {} fingerprints from {label} (radius {}, {} bits)",
                targets.fingerprints.len(),
                targets.radius,
                targets.size
            );

            if targets.fingerprints.is_empty() {
                eprintln!("no fingerprints to search");
                return Ok(exit::NO_INPUT);
            }

            let molecule = chemio::smiles::parse_smiles(query)
                .with_context(|| format!("parsing the query {query:?}"))?;

            let mut search = cli.backend.open()?;
            // The query must use the target file's parameters, not this
            // command's defaults. Fingerprints generated with a different
            // radius or size are not comparable, and comparing them anyway
            // yields plausible similarities rather than an error — which is
            // the failure mode the file header exists to prevent.
            let query_fp = pollster::block_on(search.generate_fingerprint_async(
                &molecule,
                targets.radius,
                targets.size,
            ))
            .context("fingerprinting the query")?;

            let limit = if *top == 0 {
                targets.fingerprints.len()
            } else {
                *top
            };

            let started = Instant::now();
            search
                .set_target_dataset(&targets.fingerprints)
                .context("uploading the target fingerprints")?;
            let results =
                pollster::block_on(search.search_async(&query_fp, &targets.fingerprints, limit))
                    .context("searching")?;
            eprintln!(
                "searched {} fingerprints in {:.0} ms",
                targets.fingerprints.len(),
                started.elapsed().as_secs_f64() * 1000.0
            );

            if results.is_empty() {
                eprintln!("no results");
                return Ok(exit::EMPTY_RESULT);
            }

            let mut out = String::from("rank\tname\tsimilarity\n");
            for (rank, result) in results.iter().enumerate() {
                let name = targets
                    .names
                    .get(result.index)
                    .map(String::as_str)
                    .unwrap_or("?");
                out.push_str(&format!(
                    "{}\t{}\t{:.6}\n",
                    rank + 1,
                    name,
                    result.similarity
                ));
            }
            io::write_output(output.as_ref(), &out)?;
            Ok(exit::OK)
        }
    }
}

/// A tab-separated row per molecule: parseable by `cut` and `awk`, which is the
/// point of stdout being data.
fn describe(read: &io::Input) -> String {
    let mut out = String::from("name\tatoms\tbonds\tcoords\n");
    for record in &read.outcome.records {
        let molecule = &record.molecule;
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            record.name,
            molecule.num_atoms(),
            molecule.num_bonds(),
            if molecule.has_coords() { "yes" } else { "no" },
        ));
    }
    out
}
