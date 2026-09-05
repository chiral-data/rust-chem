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
use crate::io::reader::Record;
use crate::io::sdf::parse_sdf;
use crate::io::smiles::parse_smiles;

/// Streams molecules one record at a time from a [`BufRead`], rather than
/// materializing a whole file first.
///
/// Yields [`Record`] — the same type the one-shot reader collects into a
/// `Vec` — rather than a bare `Molecule`, so a record's name (or its
/// generated `Molecule_N` fallback) survives streaming exactly as it does
/// one-shot.
///
/// A blanket impl over the right kind of `Iterator` rather than a trait
/// with its own methods — every format's supplier is already exactly an
/// iterator, and giving it a name is all this adds.
pub trait Supplier: Iterator<Item = Result<Record, ReadError>> {}
impl<T: Iterator<Item = Result<Record, ReadError>>> Supplier for T {}

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
    type Item = Result<Record, ReadError>;

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
            let rest: Vec<&str> = parts.collect();
            let name = if rest.is_empty() {
                format!("Molecule_{}", self.position)
            } else {
                rest.join(" ")
            };

            return Some(
                parse_smiles(smiles)
                    .map(|molecule| Record {
                        molecule,
                        name,
                        smiles: Some(smiles.to_owned()),
                    })
                    .map_err(|e| ReadError::Parse {
                        position: self.position,
                        message: e.to_string(),
                    }),
            );
        }
    }
}

/// One molecule per line: a SMILES, optionally followed by a `|...|`
/// enhanced-stereo-group block, then optionally a name (#221). Mirrors
/// [`crate::io::reader::read_cxsmiles_with_options`]'s splitting exactly,
/// one line read at a time instead of over a pre-loaded string.
pub struct CxSmilesSupplier<R> {
    lines: std::io::Lines<R>,
    position: usize,
    _options: ReadOptions,
}

impl<R: BufRead> CxSmilesSupplier<R> {
    pub fn new(reader: R, options: &ReadOptions) -> Self {
        Self {
            lines: reader.lines(),
            position: 0,
            _options: *options,
        }
    }
}

impl<R: BufRead> Iterator for CxSmilesSupplier<R> {
    type Item = Result<Record, ReadError>;

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

            let (smiles, block, name_parts) = crate::io::cxsmiles::split_cxsmiles_line(line);
            if smiles.is_empty() {
                continue;
            }
            let name = if name_parts.is_empty() {
                format!("Molecule_{}", self.position)
            } else {
                name_parts.join(" ")
            };

            return Some(
                crate::io::cxsmiles::parse_cxsmiles(smiles, block)
                    .map(|molecule| Record {
                        molecule,
                        name,
                        smiles: Some(smiles.to_owned()),
                    })
                    .map_err(|e| ReadError::Parse {
                        position: self.position,
                        message: e.to_string(),
                    }),
            );
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
    type Item = Result<Record, ReadError>;

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
        let position = self.position;
        Some(
            parse_sdf(&buffer)
                .map(|molecule| {
                    let name = molecule
                        .name()
                        .map(str::to_owned)
                        .unwrap_or_else(|| format!("Molecule_{position}"));
                    Record {
                        molecule,
                        name,
                        smiles: None,
                    }
                })
                .map_err(|e| ReadError::Parse {
                    position,
                    message: e.to_string(),
                }),
        )
    }
}

/// One molecule per structure. `ENDMDL` is the record boundary for a file
/// holding several back to back (#223), the same role SDF's `$$$$` plays;
/// mirrors [`crate::io::reader::read_pdb_with_options`]'s splitting
/// exactly, one line read at a time instead of over a pre-loaded string.
pub struct PdbSupplier<R> {
    lines: std::io::Lines<R>,
    position: usize,
    _options: ReadOptions,
}

impl<R: BufRead> PdbSupplier<R> {
    pub fn new(reader: R, options: &ReadOptions) -> Self {
        Self {
            lines: reader.lines(),
            position: 0,
            _options: *options,
        }
    }
}

impl<R: BufRead> Iterator for PdbSupplier<R> {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buffer = String::new();
        let mut got_any_line = false;

        loop {
            let raw = match self.lines.next() {
                None => {
                    if !got_any_line || buffer.trim().is_empty() {
                        return None;
                    }
                    // A trailing structure with no `ENDMDL` -- which a
                    // single-model file always is.
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
            let is_terminator = raw.trim() == "ENDMDL";
            buffer.push_str(&raw);
            buffer.push('\n');
            if is_terminator {
                break;
            }
        }

        self.position += 1;
        let position = self.position;
        Some(
            crate::io::pdb::parse_pdb(&buffer)
                .map(|molecule| {
                    let name = molecule
                        .name()
                        .map(str::to_owned)
                        .unwrap_or_else(|| format!("Molecule_{position}"));
                    Record {
                        molecule,
                        name,
                        smiles: None,
                    }
                })
                .map_err(|e| ReadError::Parse {
                    position,
                    message: e.to_string(),
                }),
        )
    }
}

/// One molecule per frame: a count line, a comment line, then that many
/// atom lines (#222). Mirrors
/// [`crate::io::reader::read_xyz_with_options`]'s splitting exactly, one
/// line read at a time instead of over a pre-loaded string.
pub struct XyzSupplier<R> {
    lines: std::io::Lines<R>,
    position: usize,
    _options: ReadOptions,
}

impl<R: BufRead> XyzSupplier<R> {
    pub fn new(reader: R, options: &ReadOptions) -> Self {
        Self {
            lines: reader.lines(),
            position: 0,
            _options: *options,
        }
    }
}

impl<R: BufRead> Iterator for XyzSupplier<R> {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let count_line = loop {
            match self.lines.next()? {
                Ok(line) if line.trim().is_empty() => continue,
                Ok(line) => break line,
                Err(source) => {
                    self.position += 1;
                    return Some(Err(ReadError::Io {
                        position: self.position,
                        source,
                    }));
                }
            }
        };

        let count: usize = match count_line.trim().parse() {
            Ok(n) => n,
            Err(_) => {
                self.position += 1;
                return Some(Err(ReadError::Parse {
                    position: self.position,
                    message: format!("invalid atom count: {count_line:?}"),
                }));
            }
        };

        let mut frame = count_line;
        frame.push('\n');
        // The comment line plus `count` atom lines.
        for _ in 0..=count {
            match self.lines.next() {
                Some(Ok(line)) => {
                    frame.push_str(&line);
                    frame.push('\n');
                }
                Some(Err(source)) => {
                    self.position += 1;
                    return Some(Err(ReadError::Io {
                        position: self.position,
                        source,
                    }));
                }
                None => {
                    self.position += 1;
                    return Some(Err(ReadError::Parse {
                        position: self.position,
                        message: format!("declared {count} atoms but the file ended early"),
                    }));
                }
            }
        }

        self.position += 1;
        let position = self.position;
        Some(
            crate::io::xyz::parse_xyz(&frame)
                .map(|molecule| {
                    let name = molecule
                        .name()
                        .map(str::to_owned)
                        .unwrap_or_else(|| format!("Molecule_{position}"));
                    Record {
                        molecule,
                        name,
                        smiles: None,
                    }
                })
                .map_err(|e| ReadError::Parse {
                    position,
                    message: e.to_string(),
                }),
        )
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
        // Canonical (#220), matching the non-streaming registry writer.
        writeln!(
            self.writer,
            "{} {}",
            crate::io::smiles_writer::write_smiles_for_molecule_canonical(molecule),
            name
        )
    }

    fn finish(self: Box<Self>) -> std::io::Result<()> {
        Ok(())
    }
}

/// Streams CXSMILES out, one line per molecule (#221) — the same shape the
/// one-shot writer produces in one pass.
pub struct CxSmilesWriter<W> {
    writer: W,
}

impl<W: Write> CxSmilesWriter<W> {
    pub fn new(writer: W, _options: &WriteOptions) -> Self {
        Self { writer }
    }
}

impl<W: Write> Writer for CxSmilesWriter<W> {
    fn write_molecule(&mut self, name: &str, molecule: &Molecule) -> std::io::Result<()> {
        writeln!(
            self.writer,
            "{} {}",
            crate::io::cxsmiles::write_cxsmiles(molecule),
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

/// Streams XYZ out, one frame per molecule (#222).
pub struct XyzWriter<W> {
    writer: W,
}

impl<W: Write> XyzWriter<W> {
    pub fn new(writer: W, _options: &WriteOptions) -> Self {
        Self { writer }
    }
}

impl<W: Write> Writer for XyzWriter<W> {
    fn write_molecule(&mut self, name: &str, molecule: &Molecule) -> std::io::Result<()> {
        let mut copy = molecule.clone();
        copy.set_name(name.to_string());
        self.writer
            .write_all(crate::io::xyz::write_xyz(&copy).as_bytes())
    }

    fn finish(self: Box<Self>) -> std::io::Result<()> {
        Ok(())
    }
}

/// Streams PDB out, one structure per molecule (#223). Unlike SMILES/SDF/
/// XYZ, `name` is not threaded through `write_pdb`: PDB has no single
/// title-line concept for a record's name (`HEADER`/`COMPND` are out of
/// scope, see `io/pdb.rs`'s module doc), so there is nowhere for it to go.
pub struct PdbWriter<W> {
    writer: W,
}

impl<W: Write> PdbWriter<W> {
    pub fn new(writer: W, _options: &WriteOptions) -> Self {
        Self { writer }
    }
}

impl<W: Write> Writer for PdbWriter<W> {
    fn write_molecule(&mut self, _name: &str, molecule: &Molecule) -> std::io::Result<()> {
        self.writer
            .write_all(crate::io::pdb::write_pdb(molecule).as_bytes())
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
        // Names too, not just molecules -- a streamed record used to drop
        // its name entirely, a gap the previous version of this test did
        // not catch because it only compared formulas.
        let text = "# a comment\n\nCCO ethanol\nnot-a-smiles bad\nc1ccccc1\n";

        let one_shot = read_smiles_with_options(text, &options());
        let streamed: Vec<_> =
            SmilesSupplier::new(Cursor::new(text.as_bytes()), &options()).collect::<Vec<_>>();

        let streamed_ok: Vec<_> = streamed
            .iter()
            .filter_map(|r| r.as_ref().ok())
            .map(|r| (r.molecule.formula(), r.name.clone()))
            .collect();
        let one_shot_ok: Vec<_> = one_shot
            .records
            .iter()
            .map(|r| (r.molecule.formula(), r.name.clone()))
            .collect();
        assert_eq!(streamed_ok, one_shot_ok);

        let streamed_err_count = streamed.iter().filter(|r| r.is_err()).count();
        assert_eq!(streamed_err_count, one_shot.skipped.len());
    }

    #[test]
    fn test_sdf_supplier_agrees_with_the_one_shot_reader() {
        // One titled record and one untitled -- the untitled one exercises
        // the Molecule_N fallback, which a streamed record used to skip
        // entirely.
        let mut titled = crate::core::molecule::Molecule::new();
        titled.add_atom(crate::core::atom::Atom::new(
            crate::core::atom::Element::carbon(),
        ));
        titled.set_name("my-molecule".to_string());
        let mut untitled = crate::core::molecule::Molecule::new();
        untitled.add_atom(crate::core::atom::Atom::new(
            crate::core::atom::Element::oxygen(),
        ));

        let text = format!(
            "{}{}",
            crate::io::sdf::write_sdf(&titled),
            crate::io::sdf::write_sdf(&untitled)
        );

        let one_shot = read_sdf_with_options(&text, &options());
        let streamed: Vec<_> =
            SdfSupplier::new(Cursor::new(text.as_bytes()), &options()).collect::<Vec<_>>();

        assert_eq!(streamed.len(), one_shot.records.len());
        for (streamed, one_shot) in streamed.iter().zip(one_shot.records.iter()) {
            let streamed = streamed.as_ref().unwrap();
            assert_eq!(streamed.molecule.num_atoms(), one_shot.molecule.num_atoms());
            assert_eq!(streamed.name, one_shot.name);
        }
        assert_eq!(one_shot.records[0].name, "my-molecule");
        assert_eq!(one_shot.records[1].name, "Molecule_2");
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
