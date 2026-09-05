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
/// The format is resolved from the filename with a trailing `.gz` removed
/// first, if present, so `ligand.sdf.gz` resolves as SDF rather than
/// falling through to the unrecognised-extension default.
pub fn open_supplier(path: &Path, options: &ReadOptions) -> io::Result<Box<dyn Supplier>> {
    open_supplier_as(path, format_for_path(path), options)
}

/// [`open_supplier`], with the format given explicitly rather than resolved
/// from `path` — for a caller that already knows (or was told, e.g. `chem
/// convert --from`) which format the bytes are in regardless of the name
/// on disk.
pub fn open_supplier_as(
    path: &Path,
    format: Format,
    options: &ReadOptions,
) -> io::Result<Box<dyn Supplier>> {
    let file = BufReader::new(File::open(path)?);
    let reader = maybe_decompress(file)?;

    format.supplier(reader, options).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!("{} cannot be read, only written", format.name()),
        )
    })
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

/// The format a bare path resolves to: a trailing `.gz` removed first, if
/// present, so `ligand.sdf.gz` resolves as SDF rather than falling through
/// to the unrecognised-extension default.
fn format_for_path(path: &Path) -> Format {
    let name = path.to_string_lossy();
    let format_name = name.strip_suffix(".gz").unwrap_or(&name);
    Format::from_filename(format_name)
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

        let mut supplier = open_supplier(&path, &ReadOptions).unwrap();
        let record = supplier.next().unwrap().unwrap();
        assert_eq!(record.molecule.formula(), "C2H6O");
        assert_eq!(record.name, "ethanol");
        assert!(supplier.next().is_none());

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

        let mut supplier = open_supplier(&path, &ReadOptions).unwrap();
        let read_back = supplier.next().unwrap().unwrap();
        assert_eq!(read_back.molecule.num_atoms(), mol.num_atoms());
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

        let mut supplier = open_supplier(&path, &ReadOptions).unwrap();
        assert!(supplier.next().unwrap().is_err());

        std::fs::remove_dir_all(&dir).ok();
    }
}
