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
mod export;
mod fpfile;
mod stream;
mod write;

use anyhow::{Context, Result, bail};
use backend::Backend;
use chem::draw::structure::{StructureOptions, StructureTheme};
use chem::draw::svg::structure_to_svg;
use chem::io::reader::Format;
use clap::{Parser, Subcommand, ValueEnum};
use emath::Vec2;
use fpfile::FingerprintFile;
use std::path::PathBuf;
use std::time::Instant;
use write::OutputFormat;

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

    /// Perceive aromatic rings and write the molecules back out.
    ///
    /// Aromaticity survives a SMILES round trip as lowercase atoms, so the
    /// default output format follows the input.
    Aromatic {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long, value_enum)]
        format: Option<FormatArg>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output format. Defaults to the output file's extension, or SMILES.
        #[arg(long, value_enum)]
        out_format: Option<OutputFormat>,

        /// Allow writing over the input file.
        #[arg(long)]
        force: bool,
    },

    /// Generate 2D coordinates and write the molecules back out.
    ///
    /// Defaults to SDF, because SMILES cannot store coordinates — writing
    /// SMILES would discard the result.
    Coords {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long, value_enum)]
        format: Option<FormatArg>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output format. Defaults to SDF, which is the only one that can hold
        /// coordinates.
        #[arg(long, value_enum)]
        out_format: Option<OutputFormat>,

        /// Recompute coordinates for molecules that already have them, such as
        /// those read from an SDF.
        #[arg(long)]
        relayout: bool,

        /// Allow writing over the input file.
        #[arg(long)]
        force: bool,
    },

    /// Draw each molecule as an SVG.
    ///
    /// Coordinates are a prerequisite and SMILES carries none, so any molecule
    /// without a layout gets one — reported on stderr, since it is work the
    /// command did that was not asked for. `chem coords` is the explicit form.
    Draw {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long, value_enum)]
        format: Option<FormatArg>,

        /// Write one SVG per molecule into this directory, named after each
        /// molecule. Without it, a single structure goes to standard output.
        #[arg(long)]
        outdir: Option<PathBuf>,

        /// Write the single structure to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        #[arg(long, default_value_t = 360.0)]
        width: f32,

        #[arg(long, default_value_t = 300.0)]
        height: f32,

        /// Palette. Light by default: an SVG is bound for a document or a
        /// slide, and should not carry a dark background's colours there.
        #[arg(long, value_enum, default_value_t = Theme::Light)]
        theme: Theme,
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

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Theme {
    Light,
    Dark,
}

impl Theme {
    fn palette(self) -> StructureTheme {
        match self {
            Theme::Light => StructureTheme::light(),
            Theme::Dark => StructureTheme::dark(),
        }
    }
}

/// `chem::io::reader::Format` is not ours to derive `ValueEnum` on, so this is the
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
            let read = stream::read_input(input.as_deref(), format.map(Into::into))?;
            stream::report(&read);
            // Echoed rather than resolved: `info` computes nothing, and
            // probing for a device costs a device creation. A mistyped flag is
            // still visible here rather than three commands later.
            eprintln!("backend: {} (not used by `info`)", cli.backend.label());

            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }

            stream::write_output(output.as_ref(), &describe(&read))?;

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
            let read = stream::read_input(input.as_deref(), format.map(Into::into))?;
            stream::report(&read);
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
            stream::write_output(output.as_ref(), &file.to_text())?;

            if cli.strict && !read.outcome.skipped.is_empty() {
                return Ok(exit::PARTIAL);
            }
            Ok(exit::OK)
        }

        Command::Aromatic {
            input,
            format,
            output,
            out_format,
            force,
        } => {
            write::refuse_to_clobber_input(input.as_deref(), output.as_deref(), *force)?;
            let read = stream::read_input(input.as_deref(), format.map(Into::into))?;
            stream::report(&read);
            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }
            note_cpu_only(cli);

            let mut changed = 0;
            let mut records = Vec::with_capacity(read.outcome.records.len());
            for record in &read.outcome.records {
                let mut molecule = record.molecule.clone();
                let before = aromatic_atoms(&molecule);
                chem::io::aromaticity::detect_aromaticity(&mut molecule);
                if aromatic_atoms(&molecule) != before {
                    changed += 1;
                }
                records.push((record.name.clone(), molecule));
            }
            // Reported because "it did nothing" and "nothing needed doing" look
            // identical in the output file, and only one of them is a problem.
            eprintln!(
                "perceived aromaticity: {changed} of {} molecules changed",
                records.len()
            );

            let format = OutputFormat::resolve(*out_format, false, output.as_deref());
            eprintln!("writing {}", format.label());
            stream::write_output(output.as_ref(), &write::render(format, &records))?;

            if cli.strict && !read.outcome.skipped.is_empty() {
                return Ok(exit::PARTIAL);
            }
            Ok(exit::OK)
        }

        Command::Coords {
            input,
            format,
            output,
            out_format,
            relayout,
            force,
        } => {
            write::refuse_to_clobber_input(input.as_deref(), output.as_deref(), *force)?;
            let read = stream::read_input(input.as_deref(), format.map(Into::into))?;
            stream::report(&read);
            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }
            note_cpu_only(cli);

            let mut laid_out = 0;
            let mut kept = 0;
            let mut failed = 0;
            let mut records = Vec::with_capacity(read.outcome.records.len());
            for record in &read.outcome.records {
                let mut molecule = record.molecule.clone();
                let had = molecule.has_coords();
                let ok = if *relayout {
                    chem::core::layout::layout(&mut molecule)
                } else {
                    chem::core::layout::ensure_coords(&mut molecule)
                };
                if !ok {
                    failed += 1;
                } else if had && !*relayout {
                    kept += 1;
                } else {
                    laid_out += 1;
                }
                records.push((record.name.clone(), molecule));
            }
            eprintln!("laid out {laid_out}, kept {kept} existing, {failed} without coordinates");

            let format = OutputFormat::resolve(*out_format, true, output.as_deref());
            eprintln!("writing {}", format.label());
            stream::write_output(output.as_ref(), &write::render(format, &records))?;

            if cli.strict && !read.outcome.skipped.is_empty() {
                return Ok(exit::PARTIAL);
            }
            Ok(exit::OK)
        }

        Command::Draw {
            input,
            format,
            outdir,
            output,
            width,
            height,
            theme,
        } => {
            let read = stream::read_input(input.as_deref(), format.map(Into::into))?;
            stream::report(&read);
            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }
            note_cpu_only(cli);

            // Concatenated SVG documents are not a valid SVG, so a single
            // stream can only ever hold one structure. Saying which flag fixes
            // it beats letting someone discover it from a broken file.
            if outdir.is_none() && read.outcome.len() > 1 {
                bail!(
                    "{} molecules but no --outdir: an SVG stream holds one structure, so pass --outdir to write a file each",
                    read.outcome.len()
                );
            }

            let options = StructureOptions::default();
            let palette = theme.palette();
            let size = Vec2::new(*width, *height);

            let mut generated = 0;
            let mut rendered = Vec::with_capacity(read.outcome.len());
            for record in &read.outcome.records {
                let mut molecule = record.molecule.clone();
                if !molecule.has_coords() {
                    chem::core::layout::layout(&mut molecule);
                    generated += 1;
                }
                rendered.push(structure_to_svg(&molecule, size, &options, &palette));
            }
            if generated > 0 {
                eprintln!(
                    "generated coordinates for {generated} of {} molecules (use `chem coords` to do this explicitly)",
                    rendered.len()
                );
            }

            match outdir {
                Some(dir) => {
                    let names: Vec<String> = read
                        .outcome
                        .records
                        .iter()
                        .map(|r| r.name.clone())
                        .collect();
                    let filenames = export::unique_filenames(&names);
                    let renamed = filenames
                        .iter()
                        .zip(&names)
                        .filter(|(f, n)| *f != &chem::draw::svg::suggested_filename(n))
                        .count();
                    let files: Vec<(String, String)> =
                        filenames.into_iter().zip(rendered).collect();
                    export::write_directory(dir, &files)?;
                    eprintln!("wrote {} files to {}", files.len(), dir.display());
                    if renamed > 0 {
                        // Silence here would mean the caller believes the
                        // filenames match the molecule names, and acts on it.
                        eprintln!("{renamed} had duplicate names and were suffixed");
                    }
                }
                None => {
                    stream::write_output(output.as_ref(), &rendered[0])?;
                }
            }

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
            let (text, label) = stream::read_text(fingerprints.as_deref())?;
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

            let molecule = chem::io::smiles::parse_smiles(query)
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
            stream::write_output(output.as_ref(), &out)?;
            Ok(exit::OK)
        }
    }
}

fn aromatic_atoms(molecule: &chem::core::molecule::Molecule) -> usize {
    (0..molecule.num_atoms())
        .filter(|&i| molecule.atom(i).is_aromatic())
        .count()
}

/// Says so rather than erroring when `--backend` is passed to a CPU-only
/// command, since a script that pins the backend globally should not have to
/// special-case which subcommands can use it.
fn note_cpu_only(cli: &Cli) {
    if cli.backend != Backend::Cpu {
        eprintln!(
            "backend: cpu (this operation has no GPU path; --backend {} ignored)",
            cli.backend.label()
        );
    } else {
        eprintln!("backend: cpu");
    }
}

/// A tab-separated row per molecule: parseable by `cut` and `awk`, which is the
/// point of stdout being data.
fn describe(read: &stream::Input) -> String {
    // Two coordinate columns rather than one, because a layout and a conformer
    // are different things and a file can carry either. Folding them into a
    // single `coords` column would report "no" for a 3D structure that has
    // geometry but no depiction, which is the wrong answer to the question
    // anyone is actually asking.
    let mut out = String::from("name\tatoms\tbonds\tcoords2d\tcoords3d\n");
    for record in &read.outcome.records {
        let molecule = &record.molecule;
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            record.name,
            molecule.num_atoms(),
            molecule.num_bonds(),
            if molecule.has_coords() { "yes" } else { "no" },
            if molecule.has_coords3() { "yes" } else { "no" },
        ));
    }
    out
}
