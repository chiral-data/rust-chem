//! Path-based, gzip-transparent entry points for the streaming API (#213).
//!
//! [`crate::io::format::Format::supplier`]/[`Format::writer_stream`] work
//! over anything [`BufRead`]/[`Write`] — these two functions are the
//! convenience that resolves a format from a path, opens the file, and
//! wraps it in gzip decompression/compression when the `gzip` feature is
//! compiled in and the file is (or should be) compressed.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::io::format::Format;
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::supplier::{Supplier, Writer};

/// Opens `path` and streams molecules from it, one at a time.
///
/// The format is resolved by extension first (with a trailing `.gz` removed
/// first, if present, so `ligand.sdf.gz` resolves as SDF), content-sniffing
/// second when nothing claims the extension, and the SMILES default last
/// (#317) — see this module's private `resolve_format`.
pub fn open_supplier(path: &Path, options: &ReadOptions) -> io::Result<Box<dyn Supplier>> {
    let file = BufReader::new(File::open(path)?);
    let mut reader = maybe_decompress(file)?;
    let format = resolve_format(path, reader.as_mut())?;
    supply(reader, format, options)
}

/// [`open_supplier`], with the format given explicitly rather than resolved
/// from `path` — for a caller that already knows (or was told, e.g. `chem
/// convert --from`) which format the bytes are in regardless of the name
/// on disk. Nothing to resolve, so neither the extension check nor the
/// content sniff in `resolve_format` applies here.
pub fn open_supplier_as(
    path: &Path,
    format: Format,
    options: &ReadOptions,
) -> io::Result<Box<dyn Supplier>> {
    let file = BufReader::new(File::open(path)?);
    let reader = maybe_decompress(file)?;
    supply(reader, format, options)
}

fn supply(
    reader: Box<dyn BufRead>,
    format: Format,
    options: &ReadOptions,
) -> io::Result<Box<dyn Supplier>> {
    format.supplier(reader, options).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!("{} cannot be read, only written", format.name()),
        )
    })
}

/// Resolves `path`'s format for [`open_supplier`]: by extension first
/// (`Format::from_filename_checked`, the same rule [`format_for_path`]
/// answers), content-sniffing second (peeking `reader`'s buffer through
/// `fill_buf`, which does not consume it — the same non-destructive check
/// [`maybe_decompress`] already relies on for gzip), and the SMILES default
/// last, same as `Format::from_filename` always has (#317).
///
/// A correctly-named file never pays for the peek: extension resolution
/// short-circuits before `fill_buf` is ever called.
///
/// Two exceptions: a bare `.cif` extension resolves to mmCIF *or* CIF core
/// (#320) depending on content — the two dictionaries share that extension,
/// and nothing about the name says which one a given file is. `.mmcif`
/// stays unambiguous and skips this check entirely. A bare `.top`
/// extension resolves to PRMTOP *or* GROMACS TOP (#323) the same way —
/// AMBER and GROMACS both use it.
fn resolve_format(path: &Path, reader: &mut dyn BufRead) -> io::Result<Format> {
    let name = degzipped_name(path);
    if let Some(format) = Format::from_filename_checked(&name) {
        if format == Format::MMCIF && has_extension(&name, "cif") {
            let buf = reader.fill_buf()?;
            if let Ok(text) = std::str::from_utf8(buf)
                && crate::io::cif_core::is_small_molecule_cif(text)
            {
                return Ok(Format::CIF_CORE);
            }
        }
        if format == Format::PRMTOP && has_extension(&name, "top") {
            let buf = reader.fill_buf()?;
            if let Ok(text) = std::str::from_utf8(buf)
                && crate::io::top::is_gromacs_top(text)
            {
                return Ok(Format::TOP);
            }
        }
        return Ok(format);
    }
    let buf = reader.fill_buf()?;
    Ok(crate::io::format::sniff(buf).unwrap_or(Format::SMILES))
}

fn has_extension(name: &str, extension: &str) -> bool {
    name.rsplit_once('.')
        .is_some_and(|(_, ext)| ext.eq_ignore_ascii_case(extension))
}

/// Creates (or truncates) `path` and streams molecules to it, one at a
/// time. Compresses on write when `path` ends `.gz` and the `gzip` feature
/// is compiled in.
pub fn open_writer(path: &Path, options: &WriteOptions) -> io::Result<Box<dyn Writer>> {
    open_writer_as(path, format_for_path(path), options)
}

/// [`open_writer`], with the format given explicitly rather than resolved
/// from `path` — for a caller that already knows (or was told, e.g. `chem
/// convert --to`) which format to write regardless of the name on disk.
pub fn open_writer_as(
    path: &Path,
    format: Format,
    options: &WriteOptions,
) -> io::Result<Box<dyn Writer>> {
    let name = path.to_string_lossy();
    let is_gz = name.ends_with(".gz");

    let file = BufWriter::new(File::create(path)?);
    let writer = maybe_compress(file, is_gz);

    format.writer_stream(writer, options).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!("{} cannot be written, only read", format.name()),
        )
    })
}

/// `path`'s name with a trailing `.gz` removed, if present — shared by
/// [`format_for_path`] (a pure, extension-only answer) and [`resolve_format`]
/// (which needs the same de-gzipped name before deciding whether to sniff,
/// #317).
fn degzipped_name(path: &Path) -> String {
    let name = path.to_string_lossy();
    name.strip_suffix(".gz").unwrap_or(&name).to_string()
}

/// The format a bare path resolves to: a trailing `.gz` removed first, if
/// present, so `ligand.sdf.gz` resolves as SDF rather than falling through
/// to the unrecognised-extension default.
///
/// Public because a caller that opens a supplier still has to know *which*
/// format it got, and re-deriving it means re-implementing the `.gz` rule
/// somewhere else -- which `chem convert` had already done once before #276
/// needed it a third time.
pub fn format_for_path(path: &Path) -> Format {
    Format::from_filename(&degzipped_name(path))
}

/// Peeks the first two bytes for gzip's magic number (`\x1f\x8b`) — the
/// extension is a hint at open time, the magic bytes are the actual check,
/// since a caller may have renamed a file. `fill_buf` does not consume
/// anything, so the bytes are still there for the format parser (or the
/// decompressor) to read for real afterward.
fn maybe_decompress(mut reader: BufReader<File>) -> io::Result<Box<dyn BufRead>> {
    let is_gzip = {
        let buf = reader.fill_buf()?;
        buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
    };

    if is_gzip {
        #[cfg(feature = "gzip")]
        {
            return Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(
                reader,
            ))));
        }
        // Without the `gzip` feature, the compressed bytes are handed to
        // the format parser as-is. It will fail with an ordinary parse
        // error rather than silently misreading them — no worse than
        // `Format::from_filename` already falling through on `.gz` today.
    }

    Ok(Box::new(reader))
}

fn maybe_compress(writer: BufWriter<File>, is_gz: bool) -> Box<dyn Write> {
    if is_gz {
        #[cfg(feature = "gzip")]
        {
            return Box::new(flate2::write::GzEncoder::new(
                writer,
                flate2::Compression::default(),
            ));
        }
        // Without the `gzip` feature, a `.gz`-named output is written
        // uncompressed rather than silently claiming compression it did
        // not do.
    }
    Box::new(writer)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_open_supplier_reads_a_plain_smiles_file() {
        let dir = std::env::temp_dir().join(format!("chem-open-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("plain.smi");
        std::fs::write(&path, "CCO ethanol\n").unwrap();

        let mut supplier = open_supplier(&path, &ReadOptions::default()).unwrap();
        let record = supplier.next().unwrap().unwrap();
        assert_eq!(record.molecule().unwrap().formula(), "C2H6O");
        assert_eq!(record.name, "ethanol");
        assert!(supplier.next().is_none());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_a_dot_mol_file_reads_as_sdf_not_smiles() {
        // The literal, reported bug (#318): before `mol` was added to
        // Format::SDF's extensions, this wrote a real molfile to a `.mol`
        // path and it came back as SMILES, every line skipped.
        let dir = std::env::temp_dir().join(format!("chem-mol-ext-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("ethanol.mol");
        std::fs::write(&path, include_str!("../../tests/corpus/sdf/ethanol.mol")).unwrap();

        let mut supplier = open_supplier(&path, &ReadOptions::default()).unwrap();
        let record = supplier.next().unwrap().unwrap();
        assert_eq!(record.molecule().unwrap().num_atoms(), 3);
        assert!(supplier.next().is_none());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_extension_wins_over_content() {
        // A `.smi` file whose actual bytes don't parse as SMILES at all
        // still resolves as `Format::SMILES` by extension -- sniffing is
        // never even attempted when the extension already claims something
        // (#317). It fails to *parse* (a separate concern from format
        // *resolution*), but it is read as SMILES, not defaulted elsewhere
        // by content.
        //
        // Deliberately not gzip's own magic bytes here: those would trigger
        // `maybe_decompress`'s unrelated, pre-existing detection and fail
        // for a different reason (bad deflate data), muddying what this
        // test is actually proving.
        let dir = std::env::temp_dir().join(format!("chem-ext-wins-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("probe.smi");
        std::fs::write(&path, [0x00u8, 0x01, 0x02, 0x03]).unwrap();

        let mut supplier = open_supplier(&path, &ReadOptions::default()).unwrap();
        assert!(
            supplier.next().unwrap().is_err(),
            "resolved as SMILES (by extension) and failed to parse, rather \
             than being resolved as anything else"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_unrecognized_extension_falls_back_after_a_failed_sniff() {
        // Today's real, unchanged observable behaviour, now flowing through
        // the new three-step resolver: no format claims `.probe`, `sniff`
        // finds no registered signature (none exists yet), so the SMILES
        // default still applies.
        let dir = std::env::temp_dir().join(format!("chem-unrec-ext-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("plain.probe");
        std::fs::write(&path, "CCO ethanol\n").unwrap();

        let mut supplier = open_supplier(&path, &ReadOptions::default()).unwrap();
        let record = supplier.next().unwrap().unwrap();
        assert_eq!(record.molecule().unwrap().formula(), "C2H6O");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_an_empty_file_with_an_unrecognized_extension_does_not_panic() {
        // `fill_buf` on an empty file returns an empty slice; `sniff` (and
        // everything under it) must handle that without panicking.
        let dir = std::env::temp_dir().join(format!("chem-empty-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("empty.probe");
        std::fs::write(&path, []).unwrap();

        let supplier = open_supplier(&path, &ReadOptions::default());
        assert!(supplier.is_ok());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn test_open_supplier_decompresses_a_gzipped_sdf_file() {
        use crate::io::smiles::parse_smiles;
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let dir = std::env::temp_dir().join(format!("chem-gzip-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("compressed.sdf.gz");

        let mol = parse_smiles("CCO").unwrap();
        let sdf_text = crate::io::sdf::write_sdf(&mol);

        {
            let file = File::create(&path).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder.write_all(sdf_text.as_bytes()).unwrap();
            encoder.finish().unwrap();
        }

        let mut supplier = open_supplier(&path, &ReadOptions::default()).unwrap();
        let read_back = supplier.next().unwrap().unwrap();
        assert_eq!(read_back.molecule().unwrap().num_atoms(), mol.num_atoms());
        assert!(supplier.next().is_none());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(not(feature = "gzip"))]
    #[test]
    fn test_gzipped_file_fails_cleanly_without_the_gzip_feature() {
        // A real gzip stream, still magic-byte-tagged, but nothing to
        // decompress it -- the raw compressed bytes reach the SDF parser
        // and fail as an ordinary parse error, not a panic.
        let dir = std::env::temp_dir().join(format!("chem-nogzip-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("compressed.sdf.gz");
        // Real gzip magic bytes, `1f 8b`, followed by nonsense -- enough to
        // be detected as gzip and enough to fail to parse as SDF.
        std::fs::write(&path, [0x1fu8, 0x8b, 0x08, 0x00, 0x00, 0x00]).unwrap();

        let mut supplier = open_supplier(&path, &ReadOptions::default()).unwrap();
        assert!(supplier.next().unwrap().is_err());

        std::fs::remove_dir_all(&dir).ok();
    }
}
