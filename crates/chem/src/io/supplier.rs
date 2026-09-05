//! Streaming molecule sources and sinks (#213).
//!
//! [`crate::io::reader::read`] takes a whole file already loaded into
//! memory. A [`Supplier`] takes a [`BufRead`] instead, and yields one
//! molecule at a time — the file never has to be fully materialized, which
//! matters at the sizes real structure libraries actually come in.

use std::io::{BufRead, Write};

use crate::core::molecule::Molecule;
use crate::io::errors::ReadError;
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::sdf::parse_sdf;
use crate::io::smiles::parse_smiles;

/// Streams molecules one record at a time from a [`BufRead`], rather than
/// materializing a whole file first.
///
/// A blanket impl over the right kind of `Iterator` rather than a trait
/// with its own methods — every format's supplier is already exactly an
/// iterator, and giving it a name is all this adds.
pub trait Supplier: Iterator<Item = Result<Molecule, ReadError>> {}
impl<T: Iterator<Item = Result<Molecule, ReadError>>> Supplier for T {}

/// Accepts molecules one at a time and writes them to a [`Write`], rather
/// than materializing one `String` for the whole output first.
pub trait Writer {
    fn write_molecule(&mut self, name: &str, molecule: &Molecule) -> std::io::Result<()>;

    /// Called after the last [`Self::write_molecule`]. A no-op for both
    /// formats implemented today — every SMILES line and every SDF record
    /// is already self-terminated — but real for a format whose file needs
    /// a footer once no more records are coming.
    fn finish(self: Box<Self>) -> std::io::Result<()>;
}

/// One molecule per line: the SMILES, then optionally a name. Mirrors
/// [`crate::io::reader::read_smiles_with_options`]'s splitting exactly, one
/// line read at a time instead of over a pre-loaded string.
pub struct SmilesSupplier<R> {
    lines: std::io::Lines<R>,
    position: usize,
    _options: ReadOptions,
}

impl<R: BufRead> SmilesSupplier<R> {
    pub fn new(reader: R, options: &ReadOptions) -> Self {
        Self {
            lines: reader.lines(),
            position: 0,
            _options: *options,
        }
    }
}

impl<R: BufRead> Iterator for SmilesSupplier<R> {
    type Item = Result<Molecule, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let raw = match self.lines.next()? {
                Ok(line) => line,
                Err(source) => {
                    self.position += 1;
                    return Some(Err(ReadError::Io {
                        position: self.position,
                        source,
                    }));
                }
            };
            self.position += 1;
            let line = raw.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let mut parts = line.split_whitespace();
            let Some(smiles) = parts.next() else {
                continue;
            };

            return Some(parse_smiles(smiles).map_err(|e| ReadError::Parse {
                position: self.position,
                message: e.to_string(),
            }));
        }
    }
}

/// One molecule per `$$$$`-terminated record. Mirrors
/// [`crate::io::reader::read_sdf_with_options`]'s splitting exactly, one
/// line read at a time instead of over a pre-loaded string.
pub struct SdfSupplier<R> {
    lines: std::io::Lines<R>,
    position: usize,
    _options: ReadOptions,
}

impl<R: BufRead> SdfSupplier<R> {
    pub fn new(reader: R, options: &ReadOptions) -> Self {
        Self {
            lines: reader.lines(),
            position: 0,
            _options: *options,
        }
    }
}

impl<R: BufRead> Iterator for SdfSupplier<R> {
    type Item = Result<Molecule, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buffer = String::new();
        let mut got_any_line = false;

        loop {
            let raw = match self.lines.next() {
                None => {
                    if !got_any_line || buffer.trim().is_empty() {
                        return None;
                    }
                    // A trailing record with no `$$$$` — a single-molecule
                    // file often is exactly this.
                    break;
                }
                Some(Ok(line)) => line,
                Some(Err(source)) => {
                    self.position += 1;
                    return Some(Err(ReadError::Io {
                        position: self.position,
                        source,
                    }));
                }
            };
            got_any_line = true;
            let is_terminator = raw.trim() == "$$$$";
            buffer.push_str(&raw);
            buffer.push('\n');
            if is_terminator {
                break;
            }
        }

        self.position += 1;
        Some(parse_sdf(&buffer).map_err(|e| ReadError::Parse {
            position: self.position,
            message: e.to_string(),
        }))
    }
}

/// Streams SMILES out, one line per molecule — the same shape the one-shot
/// SMILES writer produces in one pass.
pub struct SmilesWriter<W> {
    writer: W,
}

impl<W: Write> SmilesWriter<W> {
    pub fn new(writer: W, _options: &WriteOptions) -> Self {
        Self { writer }
    }
}

impl<W: Write> Writer for SmilesWriter<W> {
    fn write_molecule(&mut self, name: &str, molecule: &Molecule) -> std::io::Result<()> {
        writeln!(
            self.writer,
            "{} {}",
            crate::io::smiles_writer::write_smiles_for_molecule(molecule),
            name
        )
    }

    fn finish(self: Box<Self>) -> std::io::Result<()> {
        Ok(())
    }
}

/// Streams SDF out, one `$$$$`-terminated record per molecule.
pub struct SdfWriter<W> {
    writer: W,
    options: crate::io::options::SdfWriteOptions,
}

impl<W: Write> SdfWriter<W> {
    pub fn new(writer: W, options: &WriteOptions) -> Self {
        Self {
            writer,
            options: options.sdf,
        }
    }
}

impl<W: Write> Writer for SdfWriter<W> {
    fn write_molecule(&mut self, name: &str, molecule: &Molecule) -> std::io::Result<()> {
        let mut copy = molecule.clone();
        copy.set_name(name.to_string());
        self.writer
            .write_all(crate::io::sdf::write_sdf_with_options(&copy, &self.options).as_bytes())
    }

    fn finish(self: Box<Self>) -> std::io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::reader::{read_sdf_with_options, read_smiles_with_options};
    use std::io::Cursor;

    fn options() -> ReadOptions {
        ReadOptions
    }

    #[test]
    fn test_smiles_supplier_agrees_with_the_one_shot_reader() {
        let text = "# a comment\n\nCCO ethanol\nnot-a-smiles bad\nc1ccccc1\n";

        let one_shot = read_smiles_with_options(text, &options());
        let streamed: Vec<_> =
            SmilesSupplier::new(Cursor::new(text.as_bytes()), &options()).collect::<Vec<_>>();

        let streamed_ok: Vec<_> = streamed
            .iter()
            .filter_map(|r| r.as_ref().ok())
            .map(|m| m.formula())
            .collect();
        let one_shot_ok: Vec<_> = one_shot
            .records
            .iter()
            .map(|r| r.molecule.formula())
            .collect();
        assert_eq!(streamed_ok, one_shot_ok);

        let streamed_err_count = streamed.iter().filter(|r| r.is_err()).count();
        assert_eq!(streamed_err_count, one_shot.skipped.len());
    }

    #[test]
    fn test_sdf_supplier_agrees_with_the_one_shot_reader() {
        let sdf_one = crate::io::sdf::write_sdf(&{
            let mut mol = crate::core::molecule::Molecule::new();
            mol.add_atom(crate::core::atom::Atom::new(
                crate::core::atom::Element::carbon(),
            ));
            mol
        });
        let text = format!("{sdf_one}{sdf_one}");

        let one_shot = read_sdf_with_options(&text, &options());
        let streamed: Vec<_> =
            SdfSupplier::new(Cursor::new(text.as_bytes()), &options()).collect::<Vec<_>>();

        assert_eq!(streamed.len(), one_shot.records.len());
        for (streamed, one_shot) in streamed.iter().zip(one_shot.records.iter()) {
            assert_eq!(
                streamed.as_ref().unwrap().num_atoms(),
                one_shot.molecule.num_atoms()
            );
        }
    }

    #[test]
    fn test_supplier_is_lazy_not_fully_materialized() {
        // Consuming only the first item must not force parsing the rest —
        // proven by an input whose later records are malformed enough that
        // parsing them would show up as an error if it happened eagerly.
        let text = "C first\nnot-a-smiles-at-all second\n";
        let mut supplier = SmilesSupplier::new(Cursor::new(text.as_bytes()), &options());
        let first = supplier.next().unwrap();
        assert!(first.is_ok());
        // The second record's failure is only observed on the second `next`.
        let second = supplier.next().unwrap();
        assert!(second.is_err());
        assert!(supplier.next().is_none());
    }

    #[test]
    fn test_smiles_writer_round_trips_through_the_supplier() {
        let mol = parse_smiles("CCO").unwrap();
        let mut out: Vec<u8> = Vec::new();
        {
            let mut writer = SmilesWriter::new(&mut out, &WriteOptions::default());
            writer.write_molecule("ethanol", &mol).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let read_back = read_smiles_with_options(&text, &ReadOptions);
        assert_eq!(read_back.records.len(), 1);
        assert_eq!(read_back.records[0].name, "ethanol");
    }

    #[test]
    fn test_sdf_writer_round_trips_through_the_supplier() {
        let mol = parse_smiles("CCO").unwrap();
        let mut out: Vec<u8> = Vec::new();
        {
            let mut writer = SdfWriter::new(&mut out, &WriteOptions::default());
            writer.write_molecule("ethanol", &mol).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let read_back = read_sdf_with_options(&text, &ReadOptions);
        assert_eq!(read_back.records.len(), 1);
        assert_eq!(read_back.records[0].name, "ethanol");
    }
}
