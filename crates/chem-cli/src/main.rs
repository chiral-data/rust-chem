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
mod io;

use anyhow::Result;
use backend::Backend;
use chemio::reader::Format;
use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

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
            // Echoed so a mistyped flag is visible here rather than three
            // commands later, when a pipeline has already run on the CPU.
            eprintln!("backend: {}", cli.backend.label());

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
