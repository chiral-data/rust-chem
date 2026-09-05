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
use chem::io::format::Carries;
use chem::io::open::{open_supplier_as, open_writer_as};
use chem::io::options::{ReadOptions, WriteOptions};
use chem::io::reader::Format;
use chem::io::supplier::{Supplier, Writer};
use clap::{Parser, Subcommand, ValueEnum};
use emath::Vec2;
use fpfile::FingerprintFile;
use std::io::Cursor;
use std::path::{Path, PathBuf};
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

    /// List the molecules behind a "cannot carry" summary, one line each.
    ///
    /// Off by default because the summary is what a hundred-thousand-molecule
    /// file needs; on when a conversion surprised you and the question is
    /// which records it touched.
    #[arg(long, global = true)]
    explain_drops: bool,
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

    /// Read format A, write format B — a general converter for the pairs
    /// this crate registers a reader and a writer for.
    ///
    /// Streams rather than materializing the whole file, and reads/writes
    /// gzip transparently when the `gzip` feature is compiled in.
    ///
    /// `-o`/`--output` keeps this CLI's own established meaning (the output
    /// path, same as every other subcommand) rather than following
    /// `obabel`'s convention, where `-o` is the output *format* and `-O` is
    /// the path — `--from`/`--to` name the format instead, so nothing here
    /// means something different depending on which subcommand you're
    /// reading.
    Convert {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        /// Input format code (e.g. `smi`, `sdf`). Inferred from the input's
        /// extension when omitted, or SMILES when there is no filename to
        /// infer from (standard input, `--literal`).
        #[arg(long, short = 'i')]
        from: Option<String>,

        /// Output format code. Inferred from `--output`'s extension when
        /// omitted. Ambiguous, and an error, when neither is given and
        /// output is standard output — unlike `chem aromatic`/`chem
        /// coords`, a general converter has no format its own output
        /// naturally needs, so guessing would just be a guess.
        #[arg(long, short = 't')]
        to: Option<String>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Treat this text as the input directly, in `--from`'s format,
        /// instead of reading a file.
        #[arg(long, conflicts_with = "input")]
        literal: Option<String>,

        /// Rejected: coordinate generation is out of scope for this crate.
        /// `chem coords` computes a 2D layout; there is no 3D equivalent.
        #[arg(long)]
        gen3d: bool,

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

/// `chem::io::reader::Format` is a registry handle and not ours to derive
/// `ValueEnum` on, so this is the clap-facing mirror. Kept adjacent to its
/// conversion so the two cannot drift.
///
/// The mirror stays hand-written on purpose while the CLI exposes exactly two
/// input formats. Generating the value list from the registry is what the
/// `-L formats` story does, and doing it here first would mean two half-built
/// mechanisms.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum FormatArg {
    Smiles,
    Sdf,
}

impl From<FormatArg> for Format {
    fn from(value: FormatArg) -> Self {
        match value {
            FormatArg::Smiles => Format::SMILES,
            FormatArg::Sdf => Format::SDF,
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

            // Perceiving aromaticity adds nothing a format has to make room
            // for, so no format is required and the default stands.
            let format = OutputFormat::resolve(*out_format, Carries::empty(), output.as_deref());
            write::report_drops(format, &records, cli.explain_drops);
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
            let mut flattened = 0;
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

                // This command produces a depiction, and an SDF atom block
                // holds one set of positions. A conformer outranks a layout on
                // write — it is the more valuable data — so leaving it in place
                // would emit the input unchanged after reporting that a layout
                // was computed. Drop it, and say so rather than let someone
                // find a 3D file where they asked for a drawing.
                //
                // Reported here rather than by `write::report_drops`, and not
                // for want of trying: `Carries` is a set of per-attribute
                // capabilities, and SDF genuinely carries both COORDS_2D and
                // COORDS_3D. What it cannot do is carry them *at the same
                // time*, which a flat mask has no way to say. This loss is
                // also this command's own decision rather than the format's
                // limit, so deriving it from the mask would report a true
                // thing for a false reason.
                if ok && molecule.has_coords3() {
                    molecule.clear_coords3();
                    flattened += 1;
                }

                records.push((record.name.clone(), molecule));
            }
            eprintln!("laid out {laid_out}, kept {kept} existing, {failed} without coordinates");
            if flattened > 0 {
                eprintln!(
                    "discarded the 3D conformer of {flattened} molecules: this writes a 2D depiction, and one atom block cannot carry both"
                );
            }

            // A layout is the thing this command produced, so the default
            // format is whichever registered one can hold it — SDF — rather
            // than SDF by name.
            let format = OutputFormat::resolve(*out_format, Carries::COORDS_2D, output.as_deref());
            write::report_drops(format, &records, cli.explain_drops);
            eprintln!("writing {}", format.label());
            stream::write_output(output.as_ref(), &write::render(format, &records))?;

            if cli.strict && !read.outcome.skipped.is_empty() {
                return Ok(exit::PARTIAL);
            }
            Ok(exit::OK)
        }

        Command::Convert {
            input,
            from,
            to,
            output,
            literal,
            gen3d,
            force,
        } => {
            if *gen3d {
                bail!(
                    "chem convert does not generate 3D coordinates (out of scope for this crate) \
                     -- use `chem coords` for a 2D layout, or generate 3D externally before converting"
                );
            }
            write::refuse_to_clobber_input(input.as_deref(), output.as_deref(), *force)?;
            note_cpu_only(cli);

            let from_format = from
                .as_deref()
                .map(|code| {
                    Format::from_code(code)
                        .ok_or_else(|| anyhow::anyhow!("unrecognized format code: {code:?}"))
                })
                .transpose()?;
            let (mut supplier, label): (Box<dyn Supplier>, String) =
                resolve_convert_input(input.as_deref(), literal.as_deref(), from_format)?;

            let to_format = to
                .as_deref()
                .map(|code| {
                    Format::from_code(code)
                        .ok_or_else(|| anyhow::anyhow!("unrecognized format code: {code:?}"))
                })
                .transpose()?;
            let format = resolve_convert_output_format(to_format, output.as_deref())?;
            let mut writer = resolve_convert_output(format, output.as_deref())?;

            let mut tracker = write::DropTracker::new(format.carries());
            let mut written = 0usize;
            let mut skipped = 0usize;
            for (i, item) in supplier.by_ref().enumerate() {
                match item {
                    Ok(record) => {
                        tracker.record(&record.name, &record.molecule);
                        writer.write_molecule(&record.name, &record.molecule)?;
                        written += 1;
                    }
                    Err(e) => {
                        eprintln!("skipping record {} of {label}: {e}", i + 1);
                        skipped += 1;
                    }
                }
            }
            writer.finish()?;
            tracker.report(format.label(), cli.explain_drops);
            eprintln!("converted {written}, skipped {skipped}");

            if written == 0 {
                eprintln!("nothing readable in {label}");
                return Ok(exit::NO_INPUT);
            }
            if cli.strict && skipped > 0 {
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

/// Resolves `chem convert`'s input into a streaming [`Supplier`], and a
/// label for its skip/error messages.
///
/// `--literal` bypasses the filesystem entirely; a real path streams
/// through [`open_supplier`]/[`open_supplier_as`] (gzip-aware, per #213);
/// standard input has no name to infer a format from, so `--from` or the
/// SMILES fallback decides it, the same rule `stream::read_input` already
/// uses for `-`.
fn resolve_convert_input(
    input: Option<&Path>,
    literal: Option<&str>,
    from: Option<Format>,
) -> Result<(Box<dyn Supplier>, String)> {
    if let Some(text) = literal {
        let format = from.unwrap_or(Format::SMILES);
        let supplier = format
            .supplier(Cursor::new(text.as_bytes().to_vec()), &ReadOptions)
            .ok_or_else(|| anyhow::anyhow!("{} cannot be read, only written", format.name()))?;
        return Ok((supplier, "<literal>".to_string()));
    }

    match input.filter(|p| p.as_os_str() != "-") {
        Some(path) => {
            let supplier = match from {
                Some(format) => open_supplier_as(path, format, &ReadOptions)?,
                None => chem::io::open::open_supplier(path, &ReadOptions)?,
            };
            Ok((supplier, path.display().to_string()))
        }
        None => {
            let format = from.unwrap_or(Format::SMILES);
            let supplier = format
                .supplier(std::io::stdin().lock(), &ReadOptions)
                .ok_or_else(|| anyhow::anyhow!("{} cannot be read, only written", format.name()))?;
            Ok((supplier, "-".to_string()))
        }
    }
}

/// Resolves `chem convert`'s output format. Unlike `chem aromatic`/`chem
/// coords` (which fall back to whatever format can hold what they just
/// computed), a general converter has no such fallback — guessing here
/// would be exactly the kind of silent choice this crate's drop-report
/// philosophy exists to avoid, so an unresolvable case is an error rather
/// than a default.
fn resolve_convert_output_format(to: Option<Format>, output: Option<&Path>) -> Result<Format> {
    if let Some(format) = to {
        return Ok(format);
    }
    if let Some(path) = output.filter(|p| p.as_os_str() != "-") {
        let name = path.to_string_lossy();
        let format_name = name.strip_suffix(".gz").unwrap_or(&name);
        return Ok(Format::from_filename(format_name));
    }
    bail!(
        "output format is ambiguous — pass --to <code>, or --output <path> with a recognized extension"
    );
}

/// Resolves `chem convert`'s output into a streaming [`Writer`].
fn resolve_convert_output(format: Format, output: Option<&Path>) -> Result<Box<dyn Writer>> {
    match output.filter(|p| p.as_os_str() != "-") {
        Some(path) => Ok(open_writer_as(path, format, &WriteOptions::default())?),
        None => format
            .writer_stream(std::io::stdout().lock(), &WriteOptions::default())
            .ok_or_else(|| anyhow::anyhow!("{} cannot be written, only read", format.name())),
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
