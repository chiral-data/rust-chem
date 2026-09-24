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
use chem::io::format::{self, Carries, Kind};
use chem::io::open::{open_supplier_as, open_writer_as};
use chem::io::options::{ReadOptions, WriteOptions};
use chem::io::reader::{Format, Payload, Record};
use chem::io::supplier::{Supplier, Writer};
use clap::{Parser, Subcommand, ValueEnum};
use emath::Vec2;
use fpfile::FingerprintFile;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::time::Instant;
use write::resolve_output_format;

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

/// How `chem fp` writes its fingerprints.
#[derive(Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
enum FpFormat {
    /// This crate's own format, carrying radius and size so `chem search`
    /// can refuse a mismatched query rather than guess.
    Chem,
    /// chemfp's FPS format (#243).
    Fps,
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

        /// Input format code (e.g. `smi`, `sdf`; `chem convert -L formats`
        /// lists them). Detected from the filename otherwise; standard input
        /// has no name to read, so it defaults to SMILES.
        #[arg(long)]
        format: Option<String>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,
    },

    /// Generate Morgan fingerprints for every molecule in the input.
    Fp {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long)]
        format: Option<String>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Morgan radius.
        #[arg(long, default_value_t = 2)]
        radius: u32,

        /// Fingerprint length in bits.
        #[arg(long, default_value_t = 2048)]
        size: u32,

        /// Output format. `chem` is this crate's own self-describing file,
        /// the only one `chem search` reads. `fps` is chemfp's interchange
        /// format, for handing the bits to another tool -- write-only here,
        /// as it is everywhere else.
        #[arg(long, value_enum, default_value_t = FpFormat::Chem)]
        out_format: FpFormat,
    },

    /// Perceive aromatic rings and write the molecules back out.
    ///
    /// Aromaticity survives a SMILES round trip as lowercase atoms, so the
    /// default output format follows the input.
    Aromatic {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long)]
        format: Option<String>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output format code. Defaults to the output file's extension, or
        /// SMILES.
        #[arg(long)]
        out_format: Option<String>,

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

        #[arg(long)]
        format: Option<String>,

        /// Write to a file instead of standard output.
        #[arg(long, short = 'o')]
        output: Option<PathBuf>,

        /// Output format code. Defaults to SDF, which is the only one that
        /// can hold coordinates.
        #[arg(long)]
        out_format: Option<String>,

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

        /// List every registered format (`chem convert -L formats`), or
        /// detail one by code (`chem convert -L sdf`). A pure query — no
        /// conversion runs.
        #[arg(short = 'L')]
        list: Option<String>,

        /// Detail one format's options by code. A pure query — no
        /// conversion runs.
        #[arg(short = 'H')]
        describe: Option<String>,
    },

    /// Draw each molecule as an SVG.
    ///
    /// Coordinates are a prerequisite and SMILES carries none, so any molecule
    /// without a layout gets one — reported on stderr, since it is work the
    /// command did that was not asked for. `chem coords` is the explicit form.
    Draw {
        /// Input file. Reads standard input when absent or `-`.
        input: Option<PathBuf>,

        #[arg(long)]
        format: Option<String>,

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

/// Resolves a `--format`/`--from`/`--to`/etc. code against the registry
/// (#215) — a typo'd or not-yet-registered code is a clear error rather
/// than a guess, and every format the registry grows gains a working
/// `--format` value automatically, with no enum to edit.
fn resolve_format_code(code: Option<&str>) -> Result<Option<Format>> {
    code.map(|c| {
        Format::from_code(c).ok_or_else(|| anyhow::anyhow!("unrecognized format code: {c:?}"))
    })
    .transpose()
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
            let read =
                stream::read_input(input.as_deref(), resolve_format_code(format.as_deref())?)?;
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
            out_format,
        } => {
            let read =
                stream::read_input(input.as_deref(), resolve_format_code(format.as_deref())?)?;
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
                .filter_map(|r| stream::molecule_or_report(r, &read.label).cloned())
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
            let text = match out_format {
                FpFormat::Chem => file.to_text(),
                FpFormat::Fps => file.to_fps(),
            };
            stream::write_output(output.as_ref(), &text)?;

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
            let read =
                stream::read_input(input.as_deref(), resolve_format_code(format.as_deref())?)?;
            stream::report(&read);
            if read.outcome.is_empty() {
                eprintln!("nothing readable in {}", read.label);
                return Ok(exit::NO_INPUT);
            }
            note_cpu_only(cli);

            let mut changed = 0;
            let mut records = Vec::with_capacity(read.outcome.records.len());
            for record in &read.outcome.records {
                let Some(molecule) = stream::molecule_or_report(record, &read.label) else {
                    continue;
                };
                let mut molecule = molecule.clone();
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
            let out_format = resolve_format_code(out_format.as_deref())?;
            let format = resolve_output_format(out_format, Carries::empty(), output.as_deref());
            write::report_drops(format, &records, cli.explain_drops);
            eprintln!("writing {}", format.label());
            stream::write_output(output.as_ref(), &write::render(format, &records)?)?;

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
            let read =
                stream::read_input(input.as_deref(), resolve_format_code(format.as_deref())?)?;
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
                let Some(molecule) = stream::molecule_or_report(record, &read.label) else {
                    continue;
                };
                let mut molecule = molecule.clone();
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
            let out_format = resolve_format_code(out_format.as_deref())?;
            let format = resolve_output_format(out_format, Carries::COORDS_2D, output.as_deref());
            write::report_drops(format, &records, cli.explain_drops);
            eprintln!("writing {}", format.label());
            stream::write_output(output.as_ref(), &write::render(format, &records)?)?;

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
            list,
            describe,
        } => {
            if let Some(query) = list {
                print_format_listing(query)?;
                return Ok(exit::OK);
            }
            if let Some(code) = describe {
                print_format_options(code)?;
                return Ok(exit::OK);
            }
            if *gen3d {
                bail!(
                    "chem convert does not generate 3D coordinates (out of scope for this crate) \
                     -- use `chem coords` for a 2D layout, or generate 3D externally before converting"
                );
            }
            write::refuse_to_clobber_input(input.as_deref(), output.as_deref(), *force)?;
            note_cpu_only(cli);

            let from_format = resolve_format_code(from.as_deref())?;
            let to_format = resolve_format_code(to.as_deref())?;
            let format = resolve_convert_output_format(to_format, output.as_deref())?;

            // The Kind gate (#338): both format identities resolved with no
            // file opened on either side yet, so a genuine mismatch is
            // refused before any bytes move -- never a silent empty result.
            let source_kind_format =
                resolve_convert_input_format(input.as_deref(), literal.as_deref(), from_format);
            format::kinds_compatible(source_kind_format, format).map_err(|e| anyhow::anyhow!(e))?;

            match (source_kind_format.kind(), format.kind()) {
                (Kind::Molecules, Kind::Molecules) => {
                    let (mut supplier, label, source_format): (Box<dyn Supplier>, String, Format) =
                        resolve_convert_input(input.as_deref(), literal.as_deref(), from_format)?;
                    let mut writer = resolve_convert_output(format, output.as_deref())?;

                    let mut tracker = write::DropTracker::for_conversion(source_format, format);
                    let mut written = 0usize;
                    let mut skipped = 0usize;
                    for (i, item) in supplier.by_ref().enumerate() {
                        match item {
                            Ok(record) => match record.molecule() {
                                Some(molecule) => {
                                    tracker.record(&record.name, molecule);
                                    writer.write_molecule(&record.name, molecule)?;
                                    written += 1;
                                }
                                None => {
                                    eprintln!(
                                        "skipping record {} of {label}: not a molecule",
                                        i + 1
                                    );
                                    skipped += 1;
                                }
                            },
                            Err(e) => {
                                eprintln!("skipping record {} of {label}: {e}", i + 1);
                                skipped += 1;
                            }
                        }
                    }
                    writer.finish()?;
                    tracker.report(format.label(), cli.explain_drops);
                    // Not part of the tracker: losing atoms is a fact about
                    // the conversion rather than an attribute of a molecule,
                    // and `Carries` has no flag for it -- `held` sets
                    // TOPOLOGY on the atom count alone, so one atom of six
                    // satisfies every mask (#259).
                    if let Some(reason) = format::pair_gap(source_format, format) {
                        eprintln!("{} also loses atoms: {reason}", format.label());
                    }
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

                (Kind::Volume, Kind::Molecules) => {
                    // CUBE's own dual nature (#338) -- the only cross-Kind
                    // pair the gate allows. Keep the atoms, drop the grid,
                    // disclose it, then reuse the ordinary molecule-writing
                    // machinery for the single resulting record.
                    let Some(record) =
                        read_one_record(source_kind_format, input.as_deref(), literal.as_deref())?
                    else {
                        return Ok(exit::NO_INPUT);
                    };
                    let Payload::Volume(grid) = record.payload else {
                        bail!(
                            "internal error: {} did not produce a volume",
                            source_kind_format.name()
                        );
                    };
                    let Some(molecule) = grid.atoms() else {
                        bail!(
                            "{} claims atoms (Carries::TOPOLOGY) but this file states none",
                            source_kind_format.name()
                        );
                    };

                    let mut writer = resolve_convert_output(format, output.as_deref())?;
                    let mut tracker =
                        write::DropTracker::for_conversion(source_kind_format, format);
                    tracker.record(&record.name, molecule);
                    writer.write_molecule(&record.name, molecule)?;
                    writer.finish()?;
                    tracker.report(format.label(), cli.explain_drops);
                    // Not something `Carries` can express at all -- SAMPLES
                    // belongs to the Volume family, disjoint from a
                    // molecule's own attribute flags, so this is always true
                    // for this one exceptional path and stated directly
                    // rather than forced through the bitmask machinery.
                    eprintln!(
                        "{} discards {}'s density grid (Carries::SAMPLES)",
                        format.label(),
                        source_kind_format.label()
                    );
                    eprintln!("converted 1, skipped 0");
                    Ok(exit::OK)
                }

                (Kind::Frames, Kind::Frames) => {
                    let Some(record) =
                        read_one_record(source_kind_format, input.as_deref(), literal.as_deref())?
                    else {
                        return Ok(exit::NO_INPUT);
                    };
                    let Payload::Frames(mut trajectory) = record.payload else {
                        bail!(
                            "internal error: {} did not produce a trajectory",
                            source_kind_format.name()
                        );
                    };
                    let held = format::held_from_trajectory(&mut trajectory)?;
                    let bytes =
                        format
                            .write_trajectory_bytes(&mut trajectory)
                            .ok_or_else(|| {
                                anyhow::anyhow!("{} cannot be written, only read", format.name())
                            })?;
                    write_all_bytes(&bytes, output.as_deref())?;
                    write::report_kind_drop(format, held);
                    eprintln!("converted 1, skipped 0");
                    Ok(exit::OK)
                }

                (Kind::Volume, Kind::Volume) => {
                    let Some(record) =
                        read_one_record(source_kind_format, input.as_deref(), literal.as_deref())?
                    else {
                        return Ok(exit::NO_INPUT);
                    };
                    let Payload::Volume(grid) = record.payload else {
                        bail!(
                            "internal error: {} did not produce a volume",
                            source_kind_format.name()
                        );
                    };
                    let held = format::held_from_volume(&grid);
                    let bytes = format.write_volume_bytes(&grid).ok_or_else(|| {
                        anyhow::anyhow!("{} cannot be written, only read", format.name())
                    })?;
                    write_all_bytes(&bytes, output.as_deref())?;
                    write::report_kind_drop(format, held);
                    eprintln!("converted 1, skipped 0");
                    Ok(exit::OK)
                }

                (Kind::Mesh, Kind::Mesh) => {
                    let Some(record) =
                        read_one_record(source_kind_format, input.as_deref(), literal.as_deref())?
                    else {
                        return Ok(exit::NO_INPUT);
                    };
                    let Payload::Mesh(mesh) = record.payload else {
                        bail!(
                            "internal error: {} did not produce a mesh",
                            source_kind_format.name()
                        );
                    };
                    let bytes = format.write_mesh_bytes(&mesh).ok_or_else(|| {
                        anyhow::anyhow!("{} cannot be written, only read", format.name())
                    })?;
                    write_all_bytes(&bytes, output.as_deref())?;
                    // Neither registered Mesh format has any optional
                    // Carries flag beyond VERTICES/FACES to lose (#338) --
                    // nothing to disclose is honest, not a gap.
                    eprintln!("converted 1, skipped 0");
                    Ok(exit::OK)
                }

                (Kind::Table, Kind::Table) => {
                    let Some(record) =
                        read_one_record(source_kind_format, input.as_deref(), literal.as_deref())?
                    else {
                        return Ok(exit::NO_INPUT);
                    };
                    let Payload::Table(table) = record.payload else {
                        bail!(
                            "internal error: {} did not produce a table",
                            source_kind_format.name()
                        );
                    };
                    let bytes = format.write_table_bytes(&table).ok_or_else(|| {
                        anyhow::anyhow!("{} cannot be written, only read", format.name())
                    })?;
                    write_all_bytes(&bytes, output.as_deref())?;
                    // TABLE_METADATA (#398) is the one optional flag, and it
                    // goes through the same report as every other kind's.
                    // MDP and XVG also write only some columns, so the rest
                    // are named (#395, #398).
                    write::report_kind_drop(format, format::held_from_table(&table));
                    if format == Format::XVG {
                        let dropped = chem::io::xvg::dropped_columns(&table);
                        if !dropped.is_empty() {
                            eprintln!(
                                "{} discards column(s) {} (it writes only numeric columns)",
                                format.label(),
                                dropped.join(", ")
                            );
                        }
                    }
                    if format == Format::MDP {
                        let dropped = chem::io::mdp::dropped_columns(&table);
                        if !dropped.is_empty() {
                            eprintln!(
                                "{} discards column(s) {} (it writes only key, value, comment)",
                                format.label(),
                                dropped.join(", ")
                            );
                        }
                        if table.num_rows() > 0 && table.column("key").is_none() {
                            eprintln!(
                                "{} has no key column, so no parameters were written",
                                format.label()
                            );
                        }
                    }
                    eprintln!("converted 1, skipped 0");
                    Ok(exit::OK)
                }

                (Kind::IndexGroups, Kind::IndexGroups) => {
                    let Some(record) =
                        read_one_record(source_kind_format, input.as_deref(), literal.as_deref())?
                    else {
                        return Ok(exit::NO_INPUT);
                    };
                    let Payload::IndexGroups(groups) = record.payload else {
                        bail!(
                            "internal error: {} did not produce index groups",
                            source_kind_format.name()
                        );
                    };
                    let bytes = format.write_index_groups_bytes(&groups).ok_or_else(|| {
                        anyhow::anyhow!("{} cannot be written, only read", format.name())
                    })?;
                    write_all_bytes(&bytes, output.as_deref())?;
                    // NDX has no optional Carries flag beyond GROUPS to lose.
                    eprintln!("converted 1, skipped 0");
                    Ok(exit::OK)
                }

                (source_kind, target_kind) => bail!(
                    "internal error: kinds_compatible allowed {source_kind:?} -> {target_kind:?}, \
                     which has no dispatch arm"
                ),
            }
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
            let read =
                stream::read_input(input.as_deref(), resolve_format_code(format.as_deref())?)?;
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
            // Named alongside its rendering, not derived from `read.outcome`
            // again afterward: a name list built separately from the records
            // actually drawn would silently misalign with `rendered` the
            // moment a record here is skipped (#310) -- unreachable today,
            // since every registered format is `Kind::Molecules`, but wrong
            // to leave for whichever format makes it reachable first.
            let mut drawn: Vec<(String, String)> = Vec::with_capacity(read.outcome.len());
            for record in &read.outcome.records {
                let Some(molecule) = stream::molecule_or_report(record, &read.label) else {
                    continue;
                };
                let mut molecule = molecule.clone();
                if !molecule.has_coords() {
                    chem::core::layout::layout(&mut molecule);
                    generated += 1;
                }
                drawn.push((
                    record.name.clone(),
                    structure_to_svg(&molecule, size, &options, &palette),
                ));
            }
            if generated > 0 {
                eprintln!(
                    "generated coordinates for {generated} of {} molecules (use `chem coords` to do this explicitly)",
                    drawn.len()
                );
            }

            match outdir {
                Some(dir) => {
                    let names: Vec<String> = drawn.iter().map(|(name, _)| name.clone()).collect();
                    let rendered: Vec<String> = drawn.into_iter().map(|(_, svg)| svg).collect();
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
                    stream::write_output(output.as_ref(), &drawn[0].1)?;
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
/// Also returns the format the input was read *as*.
///
/// It was resolved here all along and thrown away, which is why the drop report
/// could only ever ask about the target (#276). `open_supplier` picks the
/// format from the path, so the same rule is applied here rather than guessed
/// at -- `chem::io::open::format_for_path` strips a trailing `.gz` first.
fn resolve_convert_input(
    input: Option<&Path>,
    literal: Option<&str>,
    from: Option<Format>,
) -> Result<(Box<dyn Supplier>, String, Format)> {
    if let Some(text) = literal {
        let format = from.unwrap_or(Format::SMILES);
        let supplier = format
            .supplier(
                Cursor::new(text.as_bytes().to_vec()),
                &ReadOptions::default(),
            )
            .ok_or_else(|| anyhow::anyhow!("{} cannot be read, only written", format.name()))?;
        return Ok((supplier, "<literal>".to_string(), format));
    }

    match input.filter(|p| p.as_os_str() != "-") {
        Some(path) => {
            let format = from.unwrap_or_else(|| chem::io::open::format_for_path(path));
            let supplier = match from {
                Some(format) => open_supplier_as(path, format, &ReadOptions::default())?,
                None => chem::io::open::open_supplier(path, &ReadOptions::default())?,
            };
            Ok((supplier, path.display().to_string(), format))
        }
        None => {
            let format = from.unwrap_or(Format::SMILES);
            let supplier = format
                .supplier(std::io::stdin().lock(), &ReadOptions::default())
                .ok_or_else(|| anyhow::anyhow!("{} cannot be read, only written", format.name()))?;
            Ok((supplier, "-".to_string(), format))
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

/// Resolves `chem convert`'s effective source format identity with no I/O at
/// all -- the exact decision `resolve_convert_input` makes internally
/// (explicit `--from`, else the path's extension, else SMILES), pulled out
/// standalone so the `Kind` gate (#338) can run before any file is opened,
/// not after `resolve_convert_input` has already read one.
fn resolve_convert_input_format(
    input: Option<&Path>,
    literal: Option<&str>,
    from: Option<Format>,
) -> Format {
    if literal.is_some() {
        return from.unwrap_or(Format::SMILES);
    }
    match input.filter(|p| p.as_os_str() != "-") {
        Some(path) => from.unwrap_or_else(|| chem::io::open::format_for_path(path)),
        None => from.unwrap_or(Format::SMILES),
    }
}

/// Reads the whole of `chem convert`'s input into memory -- what every
/// `Kind::Frames`/`Volume`/`Mesh`/`Table` conversion path needs (#338),
/// since each of those payload types is read via [`Format::read_bytes`] over
/// a whole buffer, unlike `Kind::Molecules`'s per-record streaming
/// [`Supplier`].
fn read_all_bytes(input: Option<&Path>, literal: Option<&str>) -> Result<(Vec<u8>, String)> {
    if let Some(text) = literal {
        return Ok((text.as_bytes().to_vec(), "<literal>".to_string()));
    }
    match input.filter(|p| p.as_os_str() != "-") {
        Some(path) => {
            let bytes =
                std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
            Ok((bytes, path.display().to_string()))
        }
        None => {
            let mut bytes = Vec::new();
            std::io::Read::read_to_end(&mut std::io::stdin().lock(), &mut bytes)
                .context("reading standard input")?;
            Ok((bytes, "-".to_string()))
        }
    }
}

/// Writes `bytes` to `chem convert`'s output -- the non-streaming
/// counterpart to [`resolve_convert_output`]'s `Box<dyn Writer>` (#338).
fn write_all_bytes(bytes: &[u8], output: Option<&Path>) -> Result<()> {
    match output.filter(|p| p.as_os_str() != "-") {
        Some(path) => {
            std::fs::write(path, bytes).with_context(|| format!("writing {}", path.display()))
        }
        None => std::io::Write::write_all(&mut std::io::stdout().lock(), bytes)
            .context("writing standard output"),
    }
}

/// Reads exactly one record for `chem convert`'s non-streaming `Kind` paths
/// (#338) -- every one of `Kind::Frames`/`Volume`/`Mesh`/`Table` is always
/// exactly one record per file, never a stream of many. `Ok(None)` means
/// nothing readable (already reported on stderr); the caller returns
/// `exit::NO_INPUT`.
fn read_one_record(
    source_format: Format,
    input: Option<&Path>,
    literal: Option<&str>,
) -> Result<Option<Record>> {
    let (bytes, label) = read_all_bytes(input, literal)?;
    let outcome = source_format
        .read_bytes(&bytes)
        .ok_or_else(|| anyhow::anyhow!("{} cannot be read, only written", source_format.name()))?;
    for s in &outcome.skipped {
        eprintln!("skipping record in {label}: {}", s.error);
    }
    match outcome.records.into_iter().next() {
        Some(record) => Ok(Some(record)),
        None => {
            eprintln!("nothing readable in {label}");
            Ok(None)
        }
    }
}

/// `chem convert -L formats` / `-L <code>` (#215).
fn print_format_listing(query: &str) -> Result<()> {
    if query.eq_ignore_ascii_case("formats") {
        println!("code\tname\tcategory\tread\twrite");
        for f in format::all() {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                f.codes().first().unwrap_or(&""),
                f.name(),
                f.category().label(),
                if f.can_read() { "yes" } else { "no" },
                if f.can_write() { "yes" } else { "no" },
            );
        }
        return Ok(());
    }

    if query.eq_ignore_ascii_case("matrix") {
        print_fidelity_matrix();
        return Ok(());
    }

    let format = Format::from_code(query)
        .ok_or_else(|| anyhow::anyhow!("unrecognized format code: {query:?}"))?;
    print_format_detail(format);
    Ok(())
}

/// One letter per attribute, for the matrix grid.
///
/// Single letters because a cell holds up to ten of them and eleven columns of
/// spelled-out names would not fit a terminal. `C` is formal charge and `P`
/// partial; `F` is the temperature factor, which is what the PDB spec calls
/// the B-factor column.
///
/// The seven entries below `PROPERTIES` cover `Kind::Frames`/`Volume`/`Mesh`/
/// `Table` (#339), picked from the letters the block above left free: `V`
/// velocities, `N` forces (Newtons; `F` was already b_factor), `M` frame time
/// (a moment; `T` was already topology), `L` a volume's sample values
/// (levels), `E` mesh vertices (`V` was already velocities), `K` mesh faces,
/// `Y` table columns (`C` was already formal charge). `J` is an index file's
/// groups (#394); `G` was already stereo_group. `H` is a table's header
/// metadata (#398).
const MATRIX_LEGEND: &[(Carries, char, &str)] = &[
    (Carries::TOPOLOGY, 'T', "topology"),
    (Carries::BONDS, 'B', "bonds"),
    (Carries::COORDS_2D, '2', "coords_2d"),
    (Carries::COORDS_3D, '3', "coords_3d"),
    (Carries::FORMAL_CHARGE, 'C', "formal_charge"),
    (Carries::PARTIAL_CHARGE, 'P', "partial_charge"),
    (Carries::ISOTOPE, 'I', "isotope"),
    (Carries::STEREO_ATOM, 'S', "stereo_atom"),
    (Carries::STEREO_BOND, 'D', "stereo_bond"),
    (Carries::STEREO_GROUP, 'G', "stereo_group"),
    (Carries::AROMATICITY, 'A', "aromaticity"),
    (Carries::RESIDUES, 'R', "residues"),
    (Carries::B_FACTOR, 'F', "b_factor"),
    (Carries::OCCUPANCY, 'O', "occupancy"),
    (Carries::UNIT_CELL, 'U', "unit_cell"),
    (Carries::PROPERTIES, 'X', "properties"),
    (Carries::VELOCITIES, 'V', "velocities"),
    (Carries::FORCES, 'N', "forces"),
    (Carries::FRAME_TIME, 'M', "frame_time"),
    (Carries::SAMPLES, 'L', "samples"),
    (Carries::VERTICES, 'E', "vertices"),
    (Carries::FACES, 'K', "faces"),
    (Carries::COLUMNS, 'Y', "columns"),
    (Carries::GROUPS, 'J', "groups"),
    (Carries::TABLE_METADATA, 'H', "table_metadata"),
];

/// `chem convert -L matrix` (#257) — what survives every registered conversion.
///
/// A pure query, like the rest of `-L`: the cells come from the registry, not
/// from running conversions. What makes them true is
/// `test_every_format_pair_delivers_what_the_matrix_says`, which measures
/// every pair `format::fidelity_pairs()` yields — 352 today (#339), same-kind
/// pairs across all five kinds plus CUBE's 17 `Kind::Volume` -> `Kind::
/// Molecules` cross-kind pairs — never the naive 29x29 = 841 full product.
fn print_fidelity_matrix() {
    let formats: Vec<Format> = format::all().collect();
    // Rows and columns still span every registered format (a format's own
    // diagonal cell always applies), but a cell is only ever computed for a
    // pair `fidelity_pairs()` actually yields -- otherwise blank. Without
    // this, a cell like CSV x CUBE would render `format::fidelity`'s answer
    // for a pair `kinds_compatible` itself refuses, the exact "growing to
    // 812 pairs" #339 exists to close.
    let valid_pairs: std::collections::HashSet<(Format, Format)> =
        format::fidelity_pairs().collect();
    let cell = |source: Format, target: Format| -> String {
        if !valid_pairs.contains(&(source, target)) {
            return String::new();
        }
        let delivered = format::fidelity(source, target);
        // Lowercase for an attribute the target's writer manufactures rather
        // than carries across. Without the distinction `smi -> pdb` reads as
        // though a B-factor survived, when the column was filled with 0.00.
        let mut text: String = MATRIX_LEGEND
            .iter()
            .filter(|(flag, _, _)| delivered.contains(*flag))
            .map(|(flag, letter, _)| {
                if source.carries().contains(*flag) {
                    *letter
                } else {
                    letter.to_ascii_lowercase()
                }
            })
            .collect();
        if !format::pair_loss(source, target).is_empty() {
            text.push('*');
        }
        if format::pair_gap(source, target).is_some() {
            text.push('!');
        }
        text
    };

    // Columns sized to their own content: the widest possible cell is ten
    // letters, but most targets are far narrower and a fixed width would push
    // the table past a terminal for no gain.
    let label = formats
        .iter()
        .map(|f| f.codes()[0].len())
        .max()
        .unwrap_or(0);
    let widths: Vec<usize> = formats
        .iter()
        .map(|target| {
            formats
                .iter()
                .map(|source| cell(*source, *target).len())
                .chain(std::iter::once(target.codes()[0].len()))
                .max()
                .unwrap_or(0)
        })
        .collect();

    println!("what survives a conversion, source row to target column\n");
    print!("{:width$}", "", width = label + 2);
    for (target, width) in formats.iter().zip(&widths) {
        print!("{:<width$} ", target.codes()[0], width = width);
    }
    println!();
    for source in &formats {
        print!("{:<width$}  ", source.codes()[0], width = label);
        for (target, width) in formats.iter().zip(&widths) {
            print!("{:<width$} ", cell(*source, *target), width = width);
        }
        println!();
    }

    println!();
    for (index, (_, letter, name)) in MATRIX_LEGEND.iter().enumerate() {
        print!("{letter} {name:<15}");
        if index % 4 == 3 {
            println!();
        }
    }
    println!("\n\nlowercase  supplied by the target's writer, not carried from the source");

    println!("\n* an attribute both formats claim, lost anyway:");
    for source in &formats {
        for target in &formats {
            let lost = format::pair_loss(*source, *target);
            if lost.is_empty() {
                continue;
            }
            println!(
                "    {:<10} -> {:<10} {}: {}",
                source.codes()[0],
                target.codes()[0],
                lost.names().collect::<Vec<_>>().join(", "),
                format::pair_loss_reason(*source, *target).unwrap_or(""),
            );
        }
    }

    println!("\n! the conversion also loses atoms:");
    for source in &formats {
        for target in &formats {
            if let Some(why) = format::pair_gap(*source, *target) {
                println!(
                    "    {:<10} -> {:<10} {why}",
                    source.codes()[0],
                    target.codes()[0]
                );
            }
        }
    }
}

/// `chem convert -H <code>` (#215).
fn print_format_options(code: &str) -> Result<()> {
    let format = Format::from_code(code)
        .ok_or_else(|| anyhow::anyhow!("unrecognized format code: {code:?}"))?;
    print_format_detail(format);
    // #212 built the option-bag mechanism, but nothing has wired a
    // per-format option to a CLI flag yet -- SDF's one option
    // (MolfileVersion) isn't reachable from `chem convert` today, and
    // SMILES has no options at all. A real "describe your options"
    // mechanism earns its place once there's a second thing to describe;
    // until then this says so honestly rather than inventing detail.
    println!("no format-specific options are exposed on the command line yet");
    Ok(())
}

fn print_format_detail(format: Format) {
    println!(
        "{} ({})",
        format.name(),
        format.codes().first().unwrap_or(&"")
    );
    println!("  codes: {}", format.codes().join(", "));
    println!("  extensions: {}", format.extensions().join(", "));
    println!("  category: {}", format.category().label());
    let carries: Vec<&str> = format.carries().names().collect();
    println!(
        "  carries: {}",
        if carries.is_empty() {
            "(nothing)".to_string()
        } else {
            carries.join(", ")
        }
    );
    println!(
        "  read: {}, write: {}",
        if format.can_read() { "yes" } else { "no" },
        if format.can_write() { "yes" } else { "no" }
    );
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
        let Some(molecule) = stream::molecule_or_report(record, &read.label) else {
            continue;
        };
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
