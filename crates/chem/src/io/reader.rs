//! Reading whole files into molecules.
//!
//! [`crate::io::smiles::parse_smiles`] takes one string and
//! [`crate::io::sdf::parse_sdf`] takes one `$$$$`-terminated record. A
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

use crate::core::index_groups::IndexGroups;
use crate::core::mesh::Mesh;
use crate::core::molecule::Molecule;
use crate::core::table::Table;
use crate::core::trajectory::{Frame, FrameSource, Trajectory};
use crate::core::volume::VolumeGrid;
use crate::io::options::{MultiFrameMode, ReadOptions};
use crate::io::sdf::parse_sdf;
use crate::io::smiles::parse_smiles;

/// Re-exported so `chem::io::reader::Format` keeps resolving.
///
/// The type moved to [`crate::io::format`] when it stopped being a
/// two-variant enum and became a handle into the registry. Keeping the old
/// path working means every `Format`-typed signature in the CLI and the
/// workbench stayed as written.
pub use crate::io::format::Format;

/// What a record's payload holds.
///
/// Mirrors [`crate::io::format::Kind`] on the format that produced it. Every
/// `Kind` from the v0.9.0 milestone (#307) has a variant here as of #314.
/// `#[non_exhaustive]` anyway: a variant added past this milestone would
/// still force every exhaustive match on this crate to be revisited rather
/// than silently miscompiling.
///
/// No longer `Clone` (#311): [`Trajectory`] holds a `Box<dyn FrameSource>`,
/// which cannot derive it without every future trajectory-format backend
/// committing to being cheaply cloneable before any of them exist to say
/// whether that's always possible. Nothing in this crate or `chem-app` ever
/// cloned a whole `Record`/`Payload` — only its individual `String` fields or
/// its `Molecule` — so this costs nothing today.
#[non_exhaustive]
#[derive(Debug)]
pub enum Payload {
    Molecule(Molecule),
    /// One topology, many frames (#311). Every trajectory format
    /// (TRR/XTC/DCD/NCTRAJ/LAMMPS Trajectory, #325-#329) produces this
    /// unconditionally; XYZ, PDB, PDBQT and GRO — already multi-frame
    /// before `Trajectory` existed — produce it too when
    /// [`crate::io::options::MultiFrameMode::Frames`] is requested (#330),
    /// instead of their default `Molecule`-per-frame reading.
    Frames(Trajectory),
    /// A scalar sampled on a regular 3D grid (#312). No registered format
    /// produces this yet — it lands with whichever of #331-#334 needs it
    /// first.
    Volume(VolumeGrid),
    /// A triangulated surface (#313). No registered format produces this
    /// yet — it lands with whichever of #335-#336 needs it first.
    Mesh(Mesh),
    /// A column store (#314). No registered format produces this yet — it
    /// lands with #337 (CSV).
    Table(Table),
    /// Named groups of atom indices (#394) -- a GROMACS index file.
    IndexGroups(IndexGroups),
}

/// One record read from a file.
#[derive(Debug)]
pub struct Record {
    pub payload: Payload,
    /// The record's own name, or `Molecule_N` where the file gave none.
    pub name: String,
    /// The SMILES the molecule was read from, where there was one.
    ///
    /// `None` for SDF, which stores coordinates and connectivity rather than a
    /// SMILES string. A placeholder belongs in whatever displays this, not
    /// here — a library should not invent text for a value it does not have.
    pub smiles: Option<String>,
}

impl Record {
    /// The molecule this record holds, or `None` if its payload is not one.
    pub fn molecule(&self) -> Option<&Molecule> {
        match &self.payload {
            Payload::Molecule(m) => Some(m),
            Payload::Frames(_) => None,
            Payload::Volume(_) => None,
            Payload::Mesh(_) => None,
            Payload::Table(_) => None,
            Payload::IndexGroups(_) => None,
        }
    }

    /// The trajectory this record holds, or `None` if its payload is not one.
    pub fn trajectory(&self) -> Option<&Trajectory> {
        match &self.payload {
            Payload::Molecule(_) => None,
            Payload::Frames(t) => Some(t),
            Payload::Volume(_) => None,
            Payload::Mesh(_) => None,
            Payload::Table(_) => None,
            Payload::IndexGroups(_) => None,
        }
    }

    /// The volume grid this record holds, or `None` if its payload is not one.
    pub fn volume(&self) -> Option<&VolumeGrid> {
        match &self.payload {
            Payload::Molecule(_) => None,
            Payload::Frames(_) => None,
            Payload::Volume(v) => Some(v),
            Payload::Mesh(_) => None,
            Payload::Table(_) => None,
            Payload::IndexGroups(_) => None,
        }
    }

    /// The mesh this record holds, or `None` if its payload is not one.
    pub fn mesh(&self) -> Option<&Mesh> {
        match &self.payload {
            Payload::Molecule(_) => None,
            Payload::Frames(_) => None,
            Payload::Volume(_) => None,
            Payload::Mesh(m) => Some(m),
            Payload::Table(_) => None,
            Payload::IndexGroups(_) => None,
        }
    }

    /// The table this record holds, or `None` if its payload is not one.
    pub fn table(&self) -> Option<&Table> {
        match &self.payload {
            Payload::Molecule(_) => None,
            Payload::Frames(_) => None,
            Payload::Volume(_) => None,
            Payload::Mesh(_) => None,
            Payload::Table(t) => Some(t),
            Payload::IndexGroups(_) => None,
        }
    }

    /// The index groups this record holds, or `None` if its payload is not
    /// them.
    pub fn index_groups(&self) -> Option<&IndexGroups> {
        match &self.payload {
            Payload::Molecule(_) => None,
            Payload::Frames(_) => None,
            Payload::Volume(_) => None,
            Payload::Mesh(_) => None,
            Payload::Table(_) => None,
            Payload::IndexGroups(g) => Some(g),
        }
    }
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
///
/// No longer `Clone` (#311) — see [`Payload`]'s doc comment.
#[derive(Debug, Default)]
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
///
/// Dispatch goes through the registry rather than a `match`, so adding a
/// format is adding a table entry rather than editing every call site.
///
/// A format with no reader — none today, but the descriptor allows it, and
/// several planned formats are write-only — yields an outcome carrying one
/// skipped record saying so. That keeps the "reading a file cannot fail as a
/// whole" contract [`ReadOutcome`] documents, rather than introducing a
/// `Result` for a case the caller can already see in `skipped`.
pub fn read(content: &str, format: Format) -> ReadOutcome {
    read_with_options(content, format, &ReadOptions::default())
}

/// [`read`], with explicit per-format options (#212). No format has a read
/// option yet; this exists so a future one only widens this signature once.
///
/// Routed through [`Format::read_bytes_with_options`] (#309) so there is one
/// canonical read path, text or binary — decoding `content` back to bytes is
/// always lossless since it started as a `&str`.
pub fn read_with_options(content: &str, format: Format, options: &ReadOptions) -> ReadOutcome {
    match format.read_bytes_with_options(content.as_bytes(), options) {
        Some(outcome) => outcome,
        None => ReadOutcome {
            records: Vec::new(),
            skipped: vec![Skipped {
                position: 1,
                input: String::new(),
                error: format!("{} cannot be read, only written", format.name()),
            }],
        },
    }
}

/// [`read_smiles_with_options`] with default options.
pub fn read_smiles(content: &str) -> ReadOutcome {
    read_smiles_with_options(content, &ReadOptions::default())
}

/// One molecule per line: the SMILES, then optionally a name.
///
/// Blank lines and `#` comments are skipped silently — they are not failures,
/// and counting them as such would make every commented file look broken.
pub fn read_smiles_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
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
                payload: Payload::Molecule(molecule),
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

/// [`read_cxsmiles_with_options`] with default options.
pub fn read_cxsmiles(content: &str) -> ReadOutcome {
    read_cxsmiles_with_options(content, &ReadOptions::default())
}

/// One molecule per line: a SMILES, optionally followed by a `|...|`
/// enhanced-stereo-group block, then optionally a name (#221). A strict
/// superset of [`read_smiles_with_options`]'s own convention — a line with
/// no block reads exactly the same way.
pub fn read_cxsmiles_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();

    for (index, raw) in content.lines().enumerate() {
        let position = index + 1;
        let line = raw.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let (smiles, block, name_parts) = crate::io::cxsmiles::split_cxsmiles_line(line);
        if smiles.is_empty() {
            continue;
        }
        let name = if name_parts.is_empty() {
            format!("Molecule_{position}")
        } else {
            name_parts.join(" ")
        };

        match crate::io::cxsmiles::parse_cxsmiles(smiles, block) {
            Ok(molecule) => out.records.push(Record {
                payload: Payload::Molecule(molecule),
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

/// [`read_sdf_with_options`] with default options.
pub fn read_sdf(content: &str) -> ReadOutcome {
    read_sdf_with_options(content, &ReadOptions::default())
}

/// One molecule per `$$$$`-terminated record.
pub fn read_sdf_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
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
                payload: Payload::Molecule(molecule),
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

/// Shared by every already-multi-frame text format (XYZ, PDB, PDBQT, GRO —
/// #330): a [`FrameSource`] over frames already fully decoded into memory —
/// the same shape every binary trajectory format's own test-only
/// `VecFrames` helper takes, promoted to real code here since four formats
/// need it at once.
struct MoleculeFrames(Vec<Frame>);

impl FrameSource for MoleculeFrames {
    fn frame_count(&self) -> usize {
        self.0.len()
    }

    fn num_atoms(&self) -> usize {
        self.0.first().map(Frame::num_atoms).unwrap_or(0)
    }

    fn frame(&mut self, index: usize) -> std::io::Result<Frame> {
        Ok(self.0[index].clone())
    }
}

/// Builds a [`Trajectory`] from a file's already-independently-parsed
/// frames (#330): the first molecule, stripped of its own coordinates and
/// cell, becomes the shared topology (whatever it actually parsed — atoms,
/// bonds, residues); a later molecule's own bonds are never consulted,
/// since a `Trajectory` holds exactly one topology. A molecule missing
/// coordinates, or an empty list, is a clear error rather than a silently
/// shorter trajectory. An atom-count mismatch between frames needs no
/// check here — `Trajectory::frame` already catches it lazily
/// ([`crate::core::trajectory::TrajectoryError::FrameAtomCountMismatch`]),
/// the same reused, not reimplemented, precedent LAMMPS Trajectory (#329)
/// and NCTRAJ (#328) already established.
fn molecules_to_trajectory(molecules: Vec<Molecule>) -> Result<Trajectory, String> {
    let mut iter = molecules.into_iter();
    let first = iter.next().ok_or_else(|| "no frames found".to_string())?;

    let mut topology = first.clone();
    topology.clear_coords3();
    topology.clear_cell();

    let mut frames = Vec::new();
    for (index, molecule) in std::iter::once(first).chain(iter).enumerate() {
        let positions = molecule
            .coords3()
            .ok_or_else(|| format!("frame {} has no coordinates", index + 1))?
            .to_vec();
        frames.push(Frame {
            positions,
            velocities: None,
            forces: None,
            time: None,
            step: Some(index as u64),
            cell: molecule.cell(),
        });
    }

    Trajectory::new(topology, Box::new(MoleculeFrames(frames))).map_err(|e| e.to_string())
}

/// Assembles a [`ReadOutcome`] from an already-multi-frame text format's
/// successfully-parsed molecules and skipped positions/errors, honoring
/// [`MultiFrameMode`] (#330). `Molecules` (the default) rebuilds exactly
/// the per-molecule `Record`/name output every one of these four formats
/// already produced before this option existed.
fn assemble_multi_frame_outcome(
    mode: MultiFrameMode,
    molecules: Vec<(usize, Molecule)>,
    skipped: Vec<Skipped>,
) -> ReadOutcome {
    match mode {
        MultiFrameMode::Molecules => {
            let mut out = ReadOutcome {
                records: Vec::new(),
                skipped,
            };
            for (position, molecule) in molecules {
                let name = molecule
                    .name()
                    .map(str::to_owned)
                    .unwrap_or_else(|| format!("Molecule_{position}"));
                out.records.push(Record {
                    payload: Payload::Molecule(molecule),
                    name,
                    smiles: None,
                });
            }
            out
        }
        MultiFrameMode::Frames => {
            let name = molecules
                .first()
                .and_then(|(_, m)| m.name())
                .map(str::to_owned)
                .unwrap_or_else(|| "Molecule_1".to_string());
            let mut out = ReadOutcome {
                records: Vec::new(),
                skipped,
            };
            let frames: Vec<Molecule> = molecules.into_iter().map(|(_, m)| m).collect();
            match molecules_to_trajectory(frames) {
                Ok(trajectory) => out.records.push(Record {
                    payload: Payload::Frames(trajectory),
                    name,
                    smiles: None,
                }),
                Err(e) => out.skipped.push(Skipped {
                    position: 1,
                    input: String::new(),
                    error: e,
                }),
            }
            out
        }
    }
}

/// [`read_xyz_with_options`] with default options.
pub fn read_xyz(content: &str) -> ReadOutcome {
    read_xyz_with_options(content, &ReadOptions::default())
}

/// One molecule per frame: a count line, a comment line, then that many
/// atom lines (#222). A frame's own declared count is its record boundary,
/// the role SDF's `$$$$` terminator plays -- so a trajectory (many frames
/// back to back) reads as many records by default, or one
/// [`Payload::Frames`] record when [`MultiFrameMode::Frames`] is requested
/// (#330).
pub fn read_xyz_with_options(content: &str, options: &ReadOptions) -> ReadOutcome {
    let mut molecules: Vec<(usize, Molecule)> = Vec::new();
    let mut skipped: Vec<Skipped> = Vec::new();
    let all_lines: Vec<&str> = content.lines().collect();
    let mut position = 0;
    let mut i = 0;

    while i < all_lines.len() {
        if all_lines[i].trim().is_empty() {
            i += 1;
            continue;
        }

        let count_line = all_lines[i];
        let count: usize = match count_line.trim().parse() {
            Ok(n) => n,
            Err(_) => {
                position += 1;
                skipped.push(Skipped {
                    position,
                    input: count_line.to_string(),
                    error: format!("invalid atom count: {count_line:?}"),
                });
                i += 1;
                continue;
            }
        };

        let frame_len = 2 + count;
        if i + frame_len > all_lines.len() {
            position += 1;
            skipped.push(Skipped {
                position,
                input: count_line.to_string(),
                error: format!(
                    "declared {count} atoms but only {} lines remain",
                    all_lines.len().saturating_sub(i + 2)
                ),
            });
            break;
        }

        let frame = all_lines[i..i + frame_len].join("\n");
        position += 1;
        match crate::io::xyz::parse_xyz(&frame) {
            Ok(molecule) => molecules.push((position, molecule)),
            Err(e) => skipped.push(Skipped {
                position,
                input: frame,
                error: e.to_string(),
            }),
        }
        i += frame_len;
    }

    assemble_multi_frame_outcome(options.xyz.multi_frame, molecules, skipped)
}

/// [`read_pdb_with_options`] with default options.
pub fn read_pdb(content: &str) -> ReadOutcome {
    read_pdb_with_options(content, &ReadOptions::default())
}

/// One molecule per structure. A file may hold several back to back via
/// `MODEL`/`ENDMDL` (an NMR ensemble) -- `ENDMDL` is the record boundary,
/// the same role SDF's `$$$$` and XYZ's atom count play; a file with
/// neither `MODEL` nor `ENDMDL` at all (today's common case, a single
/// deposited structure) is one implicit record, the whole file. Reads as
/// many records by default, or one [`Payload::Frames`] record (the first
/// model's own topology — including its `CONECT` bonds — shared by every
/// frame; a later model's `CONECT` is never consulted) when
/// [`MultiFrameMode::Frames`] is requested (#330).
pub fn read_pdb_with_options(content: &str, options: &ReadOptions) -> ReadOutcome {
    let mut molecules: Vec<(usize, Molecule)> = Vec::new();
    let mut skipped: Vec<Skipped> = Vec::new();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        lines.push(line);
        if line.trim() == "ENDMDL" {
            position += 1;
            push_pdb_record(&mut molecules, &mut skipped, &lines, position);
            lines.clear();
        }
    }

    // A trailing structure with no `ENDMDL` -- which a single-model file
    // always is.
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_pdb_record(&mut molecules, &mut skipped, &lines, position);
    }

    assemble_multi_frame_outcome(options.pdb.multi_frame, molecules, skipped)
}

fn push_pdb_record(
    molecules: &mut Vec<(usize, Molecule)>,
    skipped: &mut Vec<Skipped>,
    lines: &[&str],
    position: usize,
) {
    let record = lines.join("\n");
    match crate::io::pdb::parse_pdb(&record) {
        Ok(molecule) => molecules.push((position, molecule)),
        Err(e) => skipped.push(Skipped {
            position,
            input: String::new(),
            error: e.to_string(),
        }),
    }
}

/// [`read_mmcif_with_options`] with default options.
pub fn read_mmcif(content: &str) -> ReadOutcome {
    read_mmcif_with_options(content, &ReadOptions::default())
}

/// One molecule per `data_` block (#224). A block boundary is a line
/// starting with `data_`; a file with only one is one implicit record,
/// today's common case for a single deposited structure. Within a block, a
/// varying `pdbx_PDB_model_num` is [`crate::io::mmcif::parse_mmcif`]'s own
/// job (its own doc explains why), not this splitter's.
pub fn read_mmcif_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        let starts_new_block = line.trim_start().to_ascii_lowercase().starts_with("data_");
        if starts_new_block && lines.iter().any(|l: &&str| !l.trim().is_empty()) {
            position += 1;
            push_mmcif_record(&mut out, &lines, position);
            lines.clear();
        }
        lines.push(line);
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_mmcif_record(&mut out, &lines, position);
    }

    out
}

fn push_mmcif_record(out: &mut ReadOutcome, lines: &[&str], position: usize) {
    let record = lines.join("\n");
    match crate::io::mmcif::parse_mmcif(&record) {
        Ok(molecule) => {
            let name = molecule
                .name()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Molecule_{position}"));
            out.records.push(Record {
                payload: Payload::Molecule(molecule),
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

/// [`read_cif_core_with_options`] with default options.
pub fn read_cif_core(content: &str) -> ReadOutcome {
    read_cif_core_with_options(content, &ReadOptions::default())
}

/// One molecule per `data_` block (#320) -- the CIF-core analogue of
/// [`read_mmcif_with_options`], splitting the same dictionary-agnostic way:
/// a block boundary is a line starting with `data_`.
pub fn read_cif_core_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        let starts_new_block = line.trim_start().to_ascii_lowercase().starts_with("data_");
        if starts_new_block && lines.iter().any(|l: &&str| !l.trim().is_empty()) {
            position += 1;
            push_cif_core_record(&mut out, &lines, position);
            lines.clear();
        }
        lines.push(line);
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_cif_core_record(&mut out, &lines, position);
    }

    out
}

fn push_cif_core_record(out: &mut ReadOutcome, lines: &[&str], position: usize) {
    let record = lines.join("\n");
    match crate::io::cif_core::parse_cif_core(&record) {
        Ok(molecule) => {
            let name = molecule
                .name()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Molecule_{position}"));
            out.records.push(Record {
                payload: Payload::Molecule(molecule),
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

/// [`read_psf_with_options`] with default options.
pub fn read_psf(content: &str) -> ReadOutcome {
    read_psf_with_options(content, &ReadOptions::default())
}

/// A real PSF is always exactly one topology -- but this crate's own
/// invariant (every registered format writes and reads back as many
/// records as it was given, #267) means a *file this crate wrote* may hold
/// several, back to back, each with its own `PSF` header line acting as a
/// start marker -- mirroring `read_mmcif_with_options`'s `data_` splitting
/// exactly, since PSF has no real multi-topology convention of its own to
/// borrow instead.
pub fn read_psf_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        let starts_new_block = line.trim() == "PSF";
        if starts_new_block && lines.iter().any(|l: &&str| !l.trim().is_empty()) {
            position += 1;
            push_psf_record(&mut out, &lines, position);
            lines.clear();
        }
        lines.push(line);
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_psf_record(&mut out, &lines, position);
    }

    out
}

fn push_psf_record(out: &mut ReadOutcome, lines: &[&str], position: usize) {
    let record = lines.join("\n");
    match crate::io::psf::parse_psf(&record) {
        Ok(molecule) => {
            let name = molecule
                .name()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Molecule_{position}"));
            out.records.push(Record {
                payload: Payload::Molecule(molecule),
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

/// [`read_prmtop_with_options`] with default options.
pub fn read_prmtop(content: &str) -> ReadOutcome {
    read_prmtop_with_options(content, &ReadOptions::default())
}

/// A real PRMTOP is always exactly one topology -- but this crate's own
/// invariant (every registered format writes and reads back as many
/// records as it was given, #267) means a *file this crate wrote* may hold
/// several, back to back, each with its own `%VERSION` header line acting
/// as a start marker -- mirroring [`read_psf_with_options`]'s splitting
/// exactly, since PRMTOP has no real multi-topology convention of its own
/// to borrow instead. Unlike PSF's bare `PSF` line, a real `%VERSION` line
/// carries trailing content (`VERSION_STAMP = ...`), so the marker check is
/// a prefix match, not an exact one.
pub fn read_prmtop_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        let starts_new_block = line.trim_start().starts_with("%VERSION");
        if starts_new_block && lines.iter().any(|l: &&str| !l.trim().is_empty()) {
            position += 1;
            push_prmtop_record(&mut out, &lines, position);
            lines.clear();
        }
        lines.push(line);
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_prmtop_record(&mut out, &lines, position);
    }

    out
}

fn push_prmtop_record(out: &mut ReadOutcome, lines: &[&str], position: usize) {
    let record = lines.join("\n");
    match crate::io::prmtop::parse_prmtop(&record) {
        Ok(molecule) => {
            let name = molecule
                .name()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Molecule_{position}"));
            out.records.push(Record {
                payload: Payload::Molecule(molecule),
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

/// [`read_lammps_data_with_options`] with default options.
pub fn read_lammps_data(content: &str) -> ReadOutcome {
    read_lammps_data_with_options(content, &ReadOptions::default())
}

/// A real LAMMPS data file is always exactly one topology -- but this
/// crate's own invariant (every registered format writes and reads back as
/// many records as it was given, #267) means a *file this crate wrote* may
/// hold several, back to back. LAMMPS has no natural repeating start
/// marker the way PSF's bare `PSF` line or PRMTOP's `%VERSION` line do (a
/// real file's title/comment line is arbitrary text) -- so this crate's own
/// writer emits a fixed, literal title line, `"Written by chem"`, and that
/// is the marker this splits on, mirroring PSF/PRMTOP's exact splitting
/// idiom. A real third-party file's arbitrary title line only ever appears
/// once, so it is simply read as the ordinary skipped first line and never
/// triggers a split.
pub fn read_lammps_data_with_options(content: &str, options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        let starts_new_block = line.trim() == "Written by chem";
        if starts_new_block && lines.iter().any(|l: &&str| !l.trim().is_empty()) {
            position += 1;
            push_lammps_record(&mut out, &lines, position, &options.lammps);
            lines.clear();
        }
        lines.push(line);
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_lammps_record(&mut out, &lines, position, &options.lammps);
    }

    out
}

fn push_lammps_record(
    out: &mut ReadOutcome,
    lines: &[&str],
    position: usize,
    options: &crate::io::options::LammpsReadOptions,
) {
    let record = lines.join("\n");
    match crate::io::lammps::parse_lammps_data(&record, options) {
        Ok(molecule) => {
            out.records.push(Record {
                payload: Payload::Molecule(molecule),
                name: format!("Molecule_{position}"),
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

/// [`read_top_with_options`] with default options.
pub fn read_top(content: &str) -> ReadOutcome {
    read_top_with_options(content, &ReadOptions::default())
}

/// One record per `[ moleculetype ]` block (#323) -- unlike PSF/PRMTOP
/// (exactly one topology per call, split only to satisfy this crate's own
/// round-trip invariant), a GROMACS `.top` file legitimately defines
/// several moleculetypes already, so `top::parse_top` returns all of them
/// from one pass over the whole buffer -- the same "parse once, get many"
/// shape as `read_commonchem_with_options`, not a line-scan split.
pub fn read_top_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();

    match crate::io::top::parse_top(content) {
        Ok(records) => {
            for (name, molecule) in records {
                out.records.push(Record {
                    payload: Payload::Molecule(molecule),
                    name,
                    smiles: None,
                });
            }
        }
        Err(e) => out.skipped.push(Skipped {
            position: 1,
            input: content.to_string(),
            error: e.to_string(),
        }),
    }

    out
}

/// [`read_mol2_with_options`] with default options.
pub fn read_mol2(content: &str) -> ReadOutcome {
    read_mol2_with_options(content, &ReadOptions::default())
}

/// One molecule per `@<TRIPOS>MOLECULE` block (#225) -- a ligand library's
/// normal batch shape. The marker is a *start*, not an end, so a block
/// boundary is only recognised once a previous one has already begun
/// collecting content, the same reasoning `read_mmcif_with_options` uses
/// for `data_`.
pub fn read_mol2_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        let starts_new_block = line.trim() == "@<TRIPOS>MOLECULE";
        if starts_new_block && lines.iter().any(|l: &&str| !l.trim().is_empty()) {
            position += 1;
            push_mol2_record(&mut out, &lines, position);
            lines.clear();
        }
        lines.push(line);
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_mol2_record(&mut out, &lines, position);
    }

    out
}

fn push_mol2_record(out: &mut ReadOutcome, lines: &[&str], position: usize) {
    let record = lines.join("\n");
    match crate::io::mol2::parse_mol2(&record) {
        Ok(molecule) => {
            let name = molecule
                .name()
                .map(str::to_owned)
                .unwrap_or_else(|| format!("Molecule_{position}"));
            out.records.push(Record {
                payload: Payload::Molecule(molecule),
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

/// [`read_pdbqt_with_options`] with default options.
pub fn read_pdbqt(content: &str) -> ReadOutcome {
    read_pdbqt_with_options(content, &ReadOptions::default())
}

/// One molecule per structure. `ENDMDL` is the record boundary for a
/// `MODEL`/`ENDMDL`-wrapped multi-pose file (AutoDock Vina's real docked-
/// results output shape) (#226) -- mirrors [`read_pdb_with_options`]
/// exactly, since PDBQT borrows this framing directly from PDB, including
/// the same [`MultiFrameMode::Frames`] option (#330).
pub fn read_pdbqt_with_options(content: &str, options: &ReadOptions) -> ReadOutcome {
    let mut molecules: Vec<(usize, Molecule)> = Vec::new();
    let mut skipped: Vec<Skipped> = Vec::new();
    let mut lines: Vec<&str> = Vec::new();
    let mut position = 0;

    for line in content.lines() {
        lines.push(line);
        if line.trim() == "ENDMDL" {
            position += 1;
            push_pdbqt_record(&mut molecules, &mut skipped, &lines, position);
            lines.clear();
        }
    }
    if lines.iter().any(|line| !line.trim().is_empty()) {
        position += 1;
        push_pdbqt_record(&mut molecules, &mut skipped, &lines, position);
    }

    assemble_multi_frame_outcome(options.pdbqt.multi_frame, molecules, skipped)
}

fn push_pdbqt_record(
    molecules: &mut Vec<(usize, Molecule)>,
    skipped: &mut Vec<Skipped>,
    lines: &[&str],
    position: usize,
) {
    let record = lines.join("\n");
    match crate::io::pdbqt::parse_pdbqt(&record) {
        Ok(molecule) => molecules.push((position, molecule)),
        Err(e) => skipped.push(Skipped {
            position,
            input: String::new(),
            error: e.to_string(),
        }),
    }
}

/// [`read_gro_with_options`] with default options.
pub fn read_gro(content: &str) -> ReadOutcome {
    read_gro_with_options(content, &ReadOptions::default())
}

/// One molecule per frame: a title line, a count line, that many atom
/// lines, then a box-vector line (#227). Unlike every prior format's
/// reader, the declared count is trusted -- GRO has no structural
/// terminator at all separating the atom block from the box-vector line
/// that follows it, so the count is the only signal available; see
/// `io/gro.rs`'s module doc. Reads as many records by default, or one
/// [`Payload::Frames`] record when [`MultiFrameMode::Frames`] is requested
/// (#330) -- some tools concatenate `.gro` as an ad-hoc trajectory.
pub fn read_gro_with_options(content: &str, options: &ReadOptions) -> ReadOutcome {
    let mut molecules: Vec<(usize, Molecule)> = Vec::new();
    let mut skipped: Vec<Skipped> = Vec::new();
    let all_lines: Vec<&str> = content.lines().collect();
    let mut position = 0;
    let mut i = 0;

    while i < all_lines.len() {
        if all_lines[i].trim().is_empty() {
            i += 1;
            continue;
        }

        let count_line_idx = i + 1;
        let Some(count_line) = all_lines.get(count_line_idx) else {
            break;
        };
        let count: usize = match count_line.trim().parse() {
            Ok(n) => n,
            Err(_) => {
                position += 1;
                skipped.push(Skipped {
                    position,
                    input: (*count_line).to_string(),
                    error: format!("invalid atom count: {count_line:?}"),
                });
                i += 1;
                continue;
            }
        };

        // Title line + count line + count atom lines + one box-vector line.
        let frame_len = 3 + count;
        if i + frame_len > all_lines.len() {
            position += 1;
            skipped.push(Skipped {
                position,
                input: (*count_line).to_string(),
                error: format!(
                    "declared {count} atoms but only {} lines remain",
                    all_lines.len().saturating_sub(i + 2)
                ),
            });
            break;
        }

        let frame = all_lines[i..i + frame_len].join("\n");
        position += 1;
        match crate::io::gro::parse_gro(&frame) {
            Ok(molecule) => molecules.push((position, molecule)),
            Err(e) => skipped.push(Skipped {
                position,
                input: frame,
                error: e.to_string(),
            }),
        }
        i += frame_len;
    }

    assemble_multi_frame_outcome(options.gro.multi_frame, molecules, skipped)
}

/// [`read_cml_with_options`] with default options.
pub fn read_cml(content: &str) -> ReadOutcome {
    read_cml_with_options(content, &ReadOptions::default())
}

/// One molecule per `<molecule>` element (#228). A byte-offset scan over
/// the raw text, not a per-line scan like every prior format's framing
/// (mmCIF's `data_`, Mol2's `@<TRIPOS>MOLECULE`, PDB/PDBQT's `ENDMDL`) --
/// XML is not line-oriented, and a `<molecule>` start or `</molecule>` end
/// can fall anywhere on a line, including both on the same line for
/// compact output. Nested `<molecule>` elements are out of scope: each
/// record is assumed to run from its own `<molecule` start to the very
/// next `</molecule>` close, not a depth-balanced match -- see
/// `io/cml.rs`'s module doc.
pub fn read_cml_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    let mut position = 0;
    let mut search_from = 0;

    while let Some(start) = find_molecule_start(content, search_from) {
        let close_tag = "</molecule>";
        let (record, next_search_from) = match content[start..].find(close_tag) {
            Some(rel_end) => {
                let end = start + rel_end + close_tag.len();
                (&content[start..end], end)
            }
            None => (&content[start..], content.len()),
        };

        position += 1;
        match crate::io::cml::parse_cml(record) {
            Ok(molecule) => {
                let name = molecule
                    .name()
                    .map(str::to_owned)
                    .unwrap_or_else(|| format!("Molecule_{position}"));
                out.records.push(Record {
                    payload: Payload::Molecule(molecule),
                    name,
                    smiles: None,
                });
            }
            Err(e) => out.skipped.push(Skipped {
                position,
                input: record.to_string(),
                error: e.to_string(),
            }),
        }

        if next_search_from >= content.len() {
            break;
        }
        search_from = next_search_from;
    }

    out
}

/// [`read_commonchem_with_options`] with default options.
pub fn read_commonchem(content: &str) -> ReadOutcome {
    read_commonchem_with_options(content, &ReadOptions::default())
}

/// One molecule per entry in the document's `molecules` array (#229).
///
/// The only format here with no per-record framing at all: JSON has no
/// boundary a scan could find, and a document is valid only whole. So the
/// granularity of failure differs from every prior format -- a malformed
/// *document* is one `Skipped` for the entire file, since nothing smaller can
/// be salvaged, while a malformed *molecule* costs just its own record. That
/// is still this crate's "bad input costs one record, not the run" contract,
/// applied at the only two granularities this format has.
pub fn read_commonchem_with_options(content: &str, _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();

    match crate::io::commonchem::parse_commonchem(content) {
        Ok(records) => {
            for (name, molecule) in records {
                out.records.push(Record {
                    payload: Payload::Molecule(molecule),
                    name,
                    smiles: None,
                });
            }
        }
        Err(e) => out.skipped.push(Skipped {
            position: 1,
            input: content.to_string(),
            error: e.to_string(),
        }),
    }

    out
}

/// Finds the byte offset of the next `<molecule` start tag at or after
/// `from`, requiring the character right after `<molecule` to be a real
/// tag boundary (whitespace, `>`, `/`, or end of input) so `<moleculeList>`
/// and similar tags are never mistaken for a record start.
fn find_molecule_start(content: &str, from: usize) -> Option<usize> {
    let mut search_from = from;
    loop {
        let rel = content[search_from..].find("<molecule")?;
        let idx = search_from + rel;
        let after = &content[idx + "<molecule".len()..];
        match after.chars().next() {
            Some(c) if c.is_whitespace() || c == '>' || c == '/' => return Some(idx),
            None => return Some(idx),
            _ => search_from = idx + "<molecule".len(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_a_genuinely_unterminated_single_record_still_parses() {
        // Not `TWO_RECORDS` with its last `$$$$` trimmed off (that's the test
        // above) -- a real `.mol` file, the shape #318 fixed extension
        // resolution for, never had `$$$$` anywhere in it to begin with.
        let content = include_str!("../../tests/corpus/sdf/ethanol.mol");
        let out = read_sdf(content);
        assert_eq!(out.len(), 1);
        assert!(out.skipped.is_empty());
    }

    #[test]
    fn test_reading_an_empty_file_yields_nothing_rather_than_failing() {
        for content in ["", "\n\n", "# only a comment\n"] {
            let out = read(content, Format::SMILES);
            assert!(out.is_empty());
            assert!(out.skipped.is_empty());
        }
    }

    #[test]
    fn test_read_dispatches_on_format() {
        let sdf = read(TWO_RECORDS, Format::SDF);
        assert_eq!(sdf.len(), 2);
        assert!(
            sdf.records
                .iter()
                .all(|r| r.molecule().unwrap().num_atoms() > 0)
        );

        // The same bytes read as SMILES yield nothing at all. Before #151 the
        // `$$$$` terminators survived as atomless molecules, so this reported
        // "read 2 molecules" and exited successfully on a wrong-format file.
        let wrong = read(TWO_RECORDS, Format::SMILES);
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

    #[test]
    fn test_garbage_pdb_text_is_skipped_rather_than_read_as_an_empty_molecule() {
        let out = read("not a pdb file at all\n", Format::PDB);
        assert!(out.is_empty());
        assert_eq!(out.skipped.len(), 1);
        assert_eq!(out.skipped[0].position, 1);
    }

    #[test]
    fn test_garbage_mmcif_text_is_skipped_rather_than_read_as_an_empty_molecule() {
        let out = read("not an mmcif file at all\n", Format::MMCIF);
        assert!(out.is_empty());
        assert_eq!(out.skipped.len(), 1);
        assert_eq!(out.skipped[0].position, 1);
    }

    #[test]
    fn test_garbage_mol2_text_is_skipped_rather_than_read_as_an_empty_molecule() {
        let out = read("not a mol2 file at all\n", Format::MOL2);
        assert!(out.is_empty());
        assert_eq!(out.skipped.len(), 1);
        assert_eq!(out.skipped[0].position, 1);
    }

    #[test]
    fn test_garbage_pdbqt_text_is_skipped_rather_than_read_as_an_empty_molecule() {
        let out = read("not a pdbqt file at all\n", Format::PDBQT);
        assert!(out.is_empty());
        assert_eq!(out.skipped.len(), 1);
        assert_eq!(out.skipped[0].position, 1);
    }

    /// Two `MODEL`/`ENDMDL`-framed copies of the same bonded molecule, with
    /// different coordinates -- the shape #330's `MultiFrameMode::Frames`
    /// exists for.
    fn two_frame_bonded_pdb() -> String {
        use crate::core::atom::{Atom, Element};
        use crate::core::bond::{Bond, BondOrder};
        use crate::core::geometry::Point3;

        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.add_bond(Bond::new(0, 1, BondOrder::Single))
            .expect("valid bond");
        mol.set_coords3(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.5, 0.0, 0.0)])
            .expect("valid coords");

        let mut mol2 = mol.clone();
        mol2.set_coords3(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)])
            .expect("valid coords");

        crate::io::pdb::frame_models(&[
            crate::io::pdb::write_pdb(&mol),
            crate::io::pdb::write_pdb(&mol2),
        ])
    }

    #[test]
    fn test_pdb_default_mode_is_unchanged_multiple_molecule_records() {
        let framed = two_frame_bonded_pdb();
        let out = read_pdb(&framed);
        assert!(out.skipped.is_empty(), "{:?}", out.skipped);
        assert_eq!(out.records.len(), 2);
        assert!(
            out.records
                .iter()
                .all(|r| matches!(r.payload, Payload::Molecule(_)))
        );
    }

    #[test]
    fn test_pdb_frames_mode_shares_the_first_models_topology_and_positions_differ_per_frame() {
        let framed = two_frame_bonded_pdb();
        let options = ReadOptions {
            pdb: crate::io::options::PdbReadOptions {
                multi_frame: MultiFrameMode::Frames,
            },
            ..ReadOptions::default()
        };
        let outcome = read_pdb_with_options(&framed, &options);
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        assert_eq!(outcome.records.len(), 1);

        let mut records = outcome.records;
        let record = records.remove(0);
        let mut trajectory = match record.payload {
            Payload::Frames(t) => t,
            other => panic!("expected Payload::Frames, got {other:?}"),
        };
        assert_eq!(trajectory.num_atoms(), 2);
        assert_eq!(
            trajectory.topology().num_bonds(),
            1,
            "the shared topology keeps the bond"
        );
        assert_eq!(trajectory.frame_count(), 2);

        let f0 = trajectory.frame(0).unwrap();
        assert!(
            (f0.positions[1].x - 1.5).abs() < 1e-6,
            "{}",
            f0.positions[1].x
        );
        let f1 = trajectory.frame(1).unwrap();
        assert!(
            (f1.positions[1].x - 2.0).abs() < 1e-6,
            "{}",
            f1.positions[1].x
        );
    }

    #[test]
    fn test_xyz_frames_mode_a_later_frames_different_atom_count_is_caught_lazily() {
        // Construction only checks frame 0 against the shared topology --
        // `Trajectory::frame` itself is what catches a later frame that
        // declared a different count, the same reused-not-reimplemented
        // precedent LAMMPS Trajectory (#329) and NCTRAJ (#328) established.
        let text = "\
2
frame 0
C 0.0 0.0 0.0
O 1.0 0.0 0.0
3
frame 1
C 0.0 0.0 0.0
O 1.0 0.0 0.0
H 2.0 0.0 0.0
";
        let options = ReadOptions {
            xyz: crate::io::options::XyzReadOptions {
                multi_frame: MultiFrameMode::Frames,
            },
            ..ReadOptions::default()
        };
        let outcome = read_xyz_with_options(text, &options);
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        assert_eq!(outcome.records.len(), 1);

        let mut records = outcome.records;
        let record = records.remove(0);
        let mut trajectory = match record.payload {
            Payload::Frames(t) => t,
            other => panic!("expected Payload::Frames, got {other:?}"),
        };
        assert!(trajectory.frame(0).is_ok());
        let err = trajectory.frame(1).unwrap_err();
        assert!(
            matches!(
                err,
                crate::core::trajectory::TrajectoryError::FrameAtomCountMismatch { .. }
            ),
            "{err}"
        );
    }

    #[test]
    fn test_frames_mode_on_an_empty_input_is_one_clear_skipped_entry() {
        let options = ReadOptions {
            xyz: crate::io::options::XyzReadOptions {
                multi_frame: MultiFrameMode::Frames,
            },
            ..ReadOptions::default()
        };
        let outcome = read_xyz_with_options("", &options);
        assert!(outcome.records.is_empty());
        assert_eq!(outcome.skipped.len(), 1);
    }
}
