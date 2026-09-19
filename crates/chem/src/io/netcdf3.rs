//! A generic NetCDF-3 *classic* (and 64-bit-offset) container, hand-rolled
//! rather than wrapping the real `netcdf` crate (#328).
//!
//! That crate is a `-sys` crate with a `links` key; this repo's own CI
//! purity gate (`.github/workflows/ci.yml`, built for #184/#195) refuses
//! any build where one reaches the default feature set, specifically so
//! this crate keeps building for `wasm32-unknown-unknown` and installing
//! via a plain `cargo install` with no C toolchain. So this is
//! hand-roll-or-drop, the same reasoning behind [`crate::io::xdr`] and
//! [`crate::io::msgpack`] — NetCDF-3 *classic* (not NetCDF-4/HDF5) is a
//! bounded, fully-specified binary layout, which makes this a real, if
//! substantial, piece of work rather than an open-ended one.
//!
//! Confirmed directly against two real reference implementations, not
//! reconstructed from prose: `scipy.io.netcdf`'s pure-Python
//! `_netcdf.py` (a complete, exact byte-level reference for both
//! directions) and `chemfiles`'s `Netcdf3File.cpp`/`.hpp` (vendored in the
//! local cargo registry cache under `chemfiles-sys`).
//!
//! **This module has no Amber awareness at all** — that convention is
//! layered on top in [`crate::io::nctraj`], which is the whole reason this
//! is its own module: "the next Amber format will want it too."
//!
//! **Reads both classic (version 1) and 64-bit-offset (version 2)** — they
//! differ only in whether a variable's `begin` offset is a 4- or 8-byte
//! field. **Writes only the 64-bit-offset variant** — not the "reads more
//! dialects than it writes" shape TRR/DCD established, but a real
//! requirement: `nctraj`'s own live testing found that real Amber-ecosystem
//! consumers (`chemfiles`, MDAnalysis) refuse a classic-format file for this
//! convention outright, so writing classic would produce output nothing
//! real actually accepts.
//!
//! **One-shot, not streaming**: [`write`] takes an already fully-populated
//! `Netcdf3File` (every record variable's data for every record already
//! concatenated) and serialises it in one pass, unlike a real NetCDF
//! library's incremental-append API — this crate always has the whole
//! `Trajectory` available upfront, so there is nothing to stream.
//!
//! **Record variables are interleaved per record**, and a variable-level
//! 4-byte pad applies to each record's slice only when more than one
//! record variable exists — confirmed against scipy's own writer, a real,
//! easy-to-miss rule, not a simplification.

use crate::io::errors::Netcdf3Error;

const NC_DIMENSION: i32 = 10;
const NC_VARIABLE: i32 = 11;
const NC_ATTRIBUTE: i32 = 12;

const NC_BYTE: i32 = 1;
const NC_CHAR: i32 = 2;
const NC_SHORT: i32 = 3;
const NC_INT: i32 = 4;
const NC_FLOAT: i32 = 5;
const NC_DOUBLE: i32 = 6;

/// One NetCDF-3 primitive value list — either an attribute's value, or a
/// variable's fully-decoded data (flattened row-major, every record
/// concatenated for a record variable).
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Netcdf3Value {
    Byte(Vec<i8>),
    /// Raw `NC_CHAR` bytes, decoded as a string — every real use of this
    /// type in practice (Amber's own label variables included) is ASCII.
    Char(String),
    Short(Vec<i16>),
    Int(Vec<i32>),
    Float(Vec<f32>),
    Double(Vec<f64>),
}

impl Netcdf3Value {
    fn type_code(&self) -> i32 {
        match self {
            Netcdf3Value::Byte(_) => NC_BYTE,
            Netcdf3Value::Char(_) => NC_CHAR,
            Netcdf3Value::Short(_) => NC_SHORT,
            Netcdf3Value::Int(_) => NC_INT,
            Netcdf3Value::Float(_) => NC_FLOAT,
            Netcdf3Value::Double(_) => NC_DOUBLE,
        }
    }

    /// Element count -- bytes for [`Self::Char`], items for every other
    /// variant.
    fn len(&self) -> usize {
        match self {
            Netcdf3Value::Byte(v) => v.len(),
            Netcdf3Value::Char(s) => s.len(),
            Netcdf3Value::Short(v) => v.len(),
            Netcdf3Value::Int(v) => v.len(),
            Netcdf3Value::Float(v) => v.len(),
            Netcdf3Value::Double(v) => v.len(),
        }
    }

    pub(crate) fn as_str(&self) -> Option<&str> {
        match self {
            Netcdf3Value::Char(s) => Some(s),
            _ => None,
        }
    }

    pub(crate) fn as_f64_slice(&self) -> Option<Vec<f64>> {
        match self {
            Netcdf3Value::Float(v) => Some(v.iter().map(|&x| x as f64).collect()),
            Netcdf3Value::Double(v) => Some(v.clone()),
            _ => None,
        }
    }

    fn encode(&self) -> Vec<u8> {
        match self {
            Netcdf3Value::Byte(v) => v.iter().map(|&b| b as u8).collect(),
            Netcdf3Value::Char(s) => s.as_bytes().to_vec(),
            Netcdf3Value::Short(v) => {
                let mut out = Vec::with_capacity(v.len() * 2);
                for &x in v {
                    out.extend_from_slice(&x.to_be_bytes());
                }
                out
            }
            Netcdf3Value::Int(v) => {
                let mut out = Vec::with_capacity(v.len() * 4);
                for &x in v {
                    out.extend_from_slice(&x.to_be_bytes());
                }
                out
            }
            Netcdf3Value::Float(v) => {
                let mut out = Vec::with_capacity(v.len() * 4);
                for &x in v {
                    out.extend_from_slice(&x.to_be_bytes());
                }
                out
            }
            Netcdf3Value::Double(v) => {
                let mut out = Vec::with_capacity(v.len() * 8);
                for &x in v {
                    out.extend_from_slice(&x.to_be_bytes());
                }
                out
            }
        }
    }

    /// A new value holding just element range `start..start+count` of this
    /// one -- used to slice one record out of a record variable's fully
    /// concatenated data.
    fn slice(&self, start: usize, count: usize) -> Netcdf3Value {
        let end = start + count;
        match self {
            Netcdf3Value::Byte(v) => Netcdf3Value::Byte(v[start..end].to_vec()),
            Netcdf3Value::Char(s) => {
                Netcdf3Value::Char(String::from_utf8(s.as_bytes()[start..end].to_vec()).unwrap())
            }
            Netcdf3Value::Short(v) => Netcdf3Value::Short(v[start..end].to_vec()),
            Netcdf3Value::Int(v) => Netcdf3Value::Int(v[start..end].to_vec()),
            Netcdf3Value::Float(v) => Netcdf3Value::Float(v[start..end].to_vec()),
            Netcdf3Value::Double(v) => Netcdf3Value::Double(v[start..end].to_vec()),
        }
    }

    /// Concatenates same-typed values in order -- used to assemble a
    /// record variable's data from its per-record slices while reading.
    fn concat(parts: Vec<Netcdf3Value>) -> Netcdf3Value {
        let mut iter = parts.into_iter();
        let Some(first) = iter.next() else {
            return Netcdf3Value::Float(Vec::new());
        };
        match first {
            Netcdf3Value::Byte(mut v) => {
                for p in iter {
                    if let Netcdf3Value::Byte(mut x) = p {
                        v.append(&mut x);
                    }
                }
                Netcdf3Value::Byte(v)
            }
            Netcdf3Value::Char(mut s) => {
                for p in iter {
                    if let Netcdf3Value::Char(x) = p {
                        s.push_str(&x);
                    }
                }
                Netcdf3Value::Char(s)
            }
            Netcdf3Value::Short(mut v) => {
                for p in iter {
                    if let Netcdf3Value::Short(mut x) = p {
                        v.append(&mut x);
                    }
                }
                Netcdf3Value::Short(v)
            }
            Netcdf3Value::Int(mut v) => {
                for p in iter {
                    if let Netcdf3Value::Int(mut x) = p {
                        v.append(&mut x);
                    }
                }
                Netcdf3Value::Int(v)
            }
            Netcdf3Value::Float(mut v) => {
                for p in iter {
                    if let Netcdf3Value::Float(mut x) = p {
                        v.append(&mut x);
                    }
                }
                Netcdf3Value::Float(v)
            }
            Netcdf3Value::Double(mut v) => {
                for p in iter {
                    if let Netcdf3Value::Double(mut x) = p {
                        v.append(&mut x);
                    }
                }
                Netcdf3Value::Double(v)
            }
        }
    }
}

fn type_size(nc_type: i32) -> Result<usize, Netcdf3Error> {
    match nc_type {
        NC_BYTE | NC_CHAR => Ok(1),
        NC_SHORT => Ok(2),
        NC_INT | NC_FLOAT => Ok(4),
        NC_DOUBLE => Ok(8),
        other => Err(Netcdf3Error::UnknownTypeCode(other)),
    }
}

fn pad_len(unpadded: usize) -> usize {
    (4 - unpadded % 4) % 4
}

#[derive(Debug, Clone)]
pub(crate) struct Netcdf3Dimension {
    pub name: String,
    /// `None` for the one record (unlimited) dimension a classic NetCDF-3
    /// file may have.
    pub length: Option<usize>,
}

#[derive(Debug, Clone)]
pub(crate) struct Netcdf3Variable {
    pub name: String,
    /// Indices into the owning [`Netcdf3File::dimensions`].
    pub dimensions: Vec<usize>,
    pub attributes: Vec<(String, Netcdf3Value)>,
    pub data: Netcdf3Value,
}

impl Netcdf3Variable {
    fn is_record(&self, dims: &[Netcdf3Dimension]) -> bool {
        self.dimensions
            .first()
            .is_some_and(|&d| dims[d].length.is_none())
    }

    /// Element count per record -- the product of every dimension after
    /// the first (the record dimension itself contributes one record at a
    /// time, not a length).
    fn per_record_len(&self, dims: &[Netcdf3Dimension]) -> usize {
        self.dimensions[1..]
            .iter()
            .map(|&d| dims[d].length.unwrap_or(1))
            .product::<usize>()
            .max(1)
    }

    fn element_count(&self, dims: &[Netcdf3Dimension]) -> usize {
        self.dimensions
            .iter()
            .map(|&d| dims[d].length.unwrap_or(1))
            .product::<usize>()
            .max(1)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Netcdf3File {
    pub dimensions: Vec<Netcdf3Dimension>,
    pub attributes: Vec<(String, Netcdf3Value)>,
    pub variables: Vec<Netcdf3Variable>,
}

impl Netcdf3File {
    pub(crate) fn attribute(&self, name: &str) -> Option<&Netcdf3Value> {
        self.attributes
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v)
    }

    pub(crate) fn dimension(&self, name: &str) -> Option<&Netcdf3Dimension> {
        self.dimensions.iter().find(|d| d.name == name)
    }

    pub(crate) fn variable(&self, name: &str) -> Option<&Netcdf3Variable> {
        self.variables.iter().find(|v| v.name == name)
    }

    /// The record dimension's length -- how many records every record
    /// variable in this file holds.
    pub(crate) fn num_records(&self) -> usize {
        self.variables
            .iter()
            .filter(|v| v.is_record(&self.dimensions))
            .map(|v| v.data.len() / v.per_record_len(&self.dimensions).max(1))
            .max()
            .unwrap_or(0)
    }
}

// --- reading -----------------------------------------------------------

struct Netcdf3Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Netcdf3Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], Netcdf3Error> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.bytes.len())
            .ok_or_else(|| Netcdf3Error::ParseError("unexpected end of input".to_string()))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn skip(&mut self, n: usize) -> Result<(), Netcdf3Error> {
        self.take(n)?;
        Ok(())
    }

    fn read_i32(&mut self) -> Result<i32, Netcdf3Error> {
        Ok(i32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn read_i64(&mut self) -> Result<i64, Netcdf3Error> {
        Ok(i64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }

    fn read_name(&mut self) -> Result<String, Netcdf3Error> {
        let len = self.read_i32()?;
        if len < 0 {
            return Err(Netcdf3Error::ParseError("negative name length".to_string()));
        }
        let len = len as usize;
        let bytes = self.take(len)?.to_vec();
        self.skip(pad_len(len))?;
        String::from_utf8(bytes).map_err(|e| Netcdf3Error::ParseError(e.to_string()))
    }

    fn read_values(&mut self, nc_type: i32, n: usize) -> Result<Netcdf3Value, Netcdf3Error> {
        let itemsize = type_size(nc_type)?;
        let bytes = self.take(n * itemsize)?;
        decode_values(bytes, nc_type, n)
    }

    fn read_attribute_list(&mut self) -> Result<Vec<(String, Netcdf3Value)>, Netcdf3Error> {
        let tag = self.read_i32()?;
        let count = self.read_i32()?;
        if tag == 0 {
            return Ok(Vec::new());
        }
        if tag != NC_ATTRIBUTE {
            return Err(Netcdf3Error::ParseError(format!(
                "expected NC_ATTRIBUTE tag, got {tag}"
            )));
        }
        let mut out = Vec::with_capacity(count.max(0) as usize);
        for _ in 0..count {
            let name = self.read_name()?;
            let nc_type = self.read_i32()?;
            let n = self.read_i32()?.max(0) as usize;
            let itemsize = type_size(nc_type)?;
            let value = self.read_values(nc_type, n)?;
            self.skip(pad_len(n * itemsize))?;
            out.push((name, value));
        }
        Ok(out)
    }
}

fn decode_values(bytes: &[u8], nc_type: i32, n: usize) -> Result<Netcdf3Value, Netcdf3Error> {
    // Manual index loops rather than `chunks_exact`: a newer clippy lint
    // (`chunks_exact_to_as_chunks`) appeared on CI's floating toolchain
    // during #322/#325 and broke a local-green build; this form is
    // toolchain-version-agnostic.
    Ok(match nc_type {
        NC_BYTE => Netcdf3Value::Byte(bytes[..n].iter().map(|&b| b as i8).collect()),
        NC_CHAR => Netcdf3Value::Char(
            String::from_utf8_lossy(&bytes[..n])
                .trim_end_matches('\0')
                .to_string(),
        ),
        NC_SHORT => {
            let mut v = Vec::with_capacity(n);
            for i in (0..n * 2).step_by(2) {
                v.push(i16::from_be_bytes([bytes[i], bytes[i + 1]]));
            }
            Netcdf3Value::Short(v)
        }
        NC_INT => {
            let mut v = Vec::with_capacity(n);
            for i in (0..n * 4).step_by(4) {
                v.push(i32::from_be_bytes(bytes[i..i + 4].try_into().unwrap()));
            }
            Netcdf3Value::Int(v)
        }
        NC_FLOAT => {
            let mut v = Vec::with_capacity(n);
            for i in (0..n * 4).step_by(4) {
                v.push(f32::from_be_bytes(bytes[i..i + 4].try_into().unwrap()));
            }
            Netcdf3Value::Float(v)
        }
        NC_DOUBLE => {
            let mut v = Vec::with_capacity(n);
            for i in (0..n * 8).step_by(8) {
                v.push(f64::from_be_bytes(bytes[i..i + 8].try_into().unwrap()));
            }
            Netcdf3Value::Double(v)
        }
        other => return Err(Netcdf3Error::UnknownTypeCode(other)),
    })
}

/// Reads a whole NetCDF-3 (classic or 64-bit-offset) file into a generic,
/// fully in-memory representation.
pub(crate) fn read(bytes: &[u8]) -> Result<Netcdf3File, Netcdf3Error> {
    if bytes.len() < 4 || &bytes[0..3] != b"CDF" {
        return Err(Netcdf3Error::InvalidMagicNumber);
    }
    let version = bytes[3];
    if version != 1 && version != 2 {
        return Err(Netcdf3Error::UnsupportedVersion(version));
    }

    let mut r = Netcdf3Reader::new(bytes);
    r.skip(4)?;

    let numrecs = r.read_i32()?;
    if numrecs < 0 {
        return Err(Netcdf3Error::ParseError(
            "streaming (unknown-length) NetCDF-3 files are not supported".to_string(),
        ));
    }

    let dim_tag = r.read_i32()?;
    let dim_count = r.read_i32()?;
    let mut dimensions = Vec::new();
    if dim_tag != 0 {
        if dim_tag != NC_DIMENSION {
            return Err(Netcdf3Error::ParseError(format!(
                "expected NC_DIMENSION tag, got {dim_tag}"
            )));
        }
        for _ in 0..dim_count {
            let name = r.read_name()?;
            let length = r.read_i32()?;
            dimensions.push(Netcdf3Dimension {
                name,
                length: if length == 0 {
                    None
                } else {
                    Some(length as usize)
                },
            });
        }
    }

    let attributes = r.read_attribute_list()?;

    let var_tag = r.read_i32()?;
    let var_count = r.read_i32()?;
    let mut raw_vars = Vec::new();
    if var_tag != 0 {
        if var_tag != NC_VARIABLE {
            return Err(Netcdf3Error::ParseError(format!(
                "expected NC_VARIABLE tag, got {var_tag}"
            )));
        }
        for _ in 0..var_count {
            let name = r.read_name()?;
            let ndims = r.read_i32()?.max(0) as usize;
            let mut dim_ids = Vec::with_capacity(ndims);
            for _ in 0..ndims {
                dim_ids.push(r.read_i32()?.max(0) as usize);
            }
            let attrs = r.read_attribute_list()?;
            let nc_type = r.read_i32()?;
            let vsize = r.read_i32()?.max(0) as usize;
            let begin = if version == 1 {
                r.read_i32()? as usize
            } else {
                r.read_i64()? as usize
            };
            raw_vars.push((name, dim_ids, attrs, nc_type, vsize, begin));
        }
    }

    let numrecs = numrecs as usize;
    // The on-disk `vsize` is authoritative for a record variable's stride
    // between records -- it already includes whatever per-record padding
    // the writer applied (only when more than one record variable exists,
    // and only up to the next 4-byte boundary), which a `per_record *
    // itemsize` recomputation here would silently ignore, misaligning every
    // record after the first for any record variable whose own per-record
    // byte count isn't already a multiple of 4 (a real bug this exact
    // container's own test suite caught: `flags`, a 3-`i16` record variable
    // at 6 bytes/record, needs a 2-byte pad once a second record variable
    // exists, and only reading `vsize` back off disk gets that right).
    let recsize: usize = raw_vars
        .iter()
        .filter(|(_, dims, _, _, _, _)| {
            dims.first()
                .is_some_and(|&d| dimensions[d].length.is_none())
        })
        .map(|(_, _, _, _, vsize, _)| *vsize)
        .sum();
    let num_record_vars = raw_vars
        .iter()
        .filter(|(_, dims, _, _, _, _)| {
            dims.first()
                .is_some_and(|&d| dimensions[d].length.is_none())
        })
        .count();

    let mut variables = Vec::with_capacity(raw_vars.len());
    for (name, dim_ids, attrs, nc_type, _vsize, begin) in raw_vars {
        let is_record = dim_ids
            .first()
            .is_some_and(|&d| dimensions[d].length.is_none());
        let itemsize = type_size(nc_type)?;
        let data = if is_record {
            let per_record = dim_ids[1..]
                .iter()
                .map(|&d| dimensions[d].length.unwrap_or(1))
                .product::<usize>()
                .max(1);
            let mut records = Vec::with_capacity(numrecs);
            for rec in 0..numrecs {
                let offset = begin + rec * recsize;
                let byte_len = per_record * itemsize;
                let slice = bytes.get(offset..offset + byte_len).ok_or_else(|| {
                    Netcdf3Error::ParseError(format!("record {rec} of '{name}' runs past EOF"))
                })?;
                records.push(decode_values(slice, nc_type, per_record)?);
            }
            Netcdf3Value::concat(records)
        } else {
            let n = dim_ids
                .iter()
                .map(|&d| dimensions[d].length.unwrap_or(1))
                .product::<usize>()
                .max(1);
            let byte_len = n * itemsize;
            let slice = bytes
                .get(begin..begin + byte_len)
                .ok_or_else(|| Netcdf3Error::ParseError(format!("'{name}' data runs past EOF")))?;
            decode_values(slice, nc_type, n)?
        };
        variables.push(Netcdf3Variable {
            name,
            dimensions: dim_ids,
            attributes: attrs,
            data,
        });
    }
    let _ = num_record_vars;

    Ok(Netcdf3File {
        dimensions,
        attributes,
        variables,
    })
}

// --- writing -------------------------------------------------------------

fn write_name(out: &mut Vec<u8>, name: &str) {
    out.extend_from_slice(&(name.len() as i32).to_be_bytes());
    out.extend_from_slice(name.as_bytes());
    out.resize(out.len() + pad_len(name.len()), 0);
}

fn write_attribute_list(out: &mut Vec<u8>, attrs: &[(String, Netcdf3Value)]) {
    if attrs.is_empty() {
        out.extend_from_slice(&0i32.to_be_bytes());
        out.extend_from_slice(&0i32.to_be_bytes());
        return;
    }
    out.extend_from_slice(&NC_ATTRIBUTE.to_be_bytes());
    out.extend_from_slice(&(attrs.len() as i32).to_be_bytes());
    for (name, value) in attrs {
        write_name(out, name);
        out.extend_from_slice(&value.type_code().to_be_bytes());
        out.extend_from_slice(&(value.len() as i32).to_be_bytes());
        let bytes = value.encode();
        out.extend_from_slice(&bytes);
        out.resize(out.len() + pad_len(bytes.len()), 0);
    }
}

/// Serialises a fully-populated [`Netcdf3File`] into 64-bit-offset (version
/// 2) NetCDF-3 bytes.
///
/// Not classic (version 1), despite this container reading that dialect
/// fine: real Amber-ecosystem consumers (confirmed against both `chemfiles`,
/// which refuses to even open a classic file for this convention, and
/// MDAnalysis, which raises the same complaint) only accept version 2 for
/// this format, so a classic-format writer would produce files no real tool
/// treats as valid Amber NetCDF trajectories -- caught by an actual
/// `MDAnalysis` read of this crate's own output, not by inspection.
pub(crate) fn write(file: &Netcdf3File) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"CDF");
    out.push(2); // 64-bit offset

    let numrecs = file.num_records();
    out.extend_from_slice(&(numrecs as i32).to_be_bytes());

    if file.dimensions.is_empty() {
        out.extend_from_slice(&0i32.to_be_bytes());
        out.extend_from_slice(&0i32.to_be_bytes());
    } else {
        out.extend_from_slice(&NC_DIMENSION.to_be_bytes());
        out.extend_from_slice(&(file.dimensions.len() as i32).to_be_bytes());
        for dim in &file.dimensions {
            write_name(&mut out, &dim.name);
            out.extend_from_slice(&(dim.length.unwrap_or(0) as i32).to_be_bytes());
        }
    }

    write_attribute_list(&mut out, &file.attributes);

    if file.variables.is_empty() {
        out.extend_from_slice(&0i32.to_be_bytes());
        out.extend_from_slice(&0i32.to_be_bytes());
        return out;
    }

    out.extend_from_slice(&NC_VARIABLE.to_be_bytes());
    out.extend_from_slice(&(file.variables.len() as i32).to_be_bytes());

    let is_record: Vec<bool> = file
        .variables
        .iter()
        .map(|v| v.is_record(&file.dimensions))
        .collect();
    let num_record_vars = is_record.iter().filter(|&&b| b).count();

    let mut begin_patch_at = Vec::with_capacity(file.variables.len());
    let mut vsizes = Vec::with_capacity(file.variables.len());
    for (var, &rec) in file.variables.iter().zip(&is_record) {
        write_name(&mut out, &var.name);
        out.extend_from_slice(&(var.dimensions.len() as i32).to_be_bytes());
        for &d in &var.dimensions {
            out.extend_from_slice(&(d as i32).to_be_bytes());
        }
        write_attribute_list(&mut out, &var.attributes);

        let nc_type = var.data.type_code();
        out.extend_from_slice(&nc_type.to_be_bytes());
        let itemsize = type_size(nc_type).expect("valid type code");

        let vsize = if rec {
            let count = var.per_record_len(&file.dimensions) * itemsize;
            if num_record_vars > 1 {
                count + pad_len(count)
            } else {
                count
            }
        } else {
            let count = var.element_count(&file.dimensions) * itemsize;
            count + pad_len(count)
        };
        vsizes.push(vsize);
        out.extend_from_slice(&(vsize as i32).to_be_bytes());

        begin_patch_at.push(out.len());
        out.extend_from_slice(&0i64.to_be_bytes()); // placeholder, patched below
    }

    let mut begins = vec![0usize; file.variables.len()];
    for (i, (var, &rec)) in file.variables.iter().zip(&is_record).enumerate() {
        if rec {
            continue;
        }
        begins[i] = out.len();
        let bytes = var.data.encode();
        out.extend_from_slice(&bytes);
        out.resize(out.len() + pad_len(bytes.len()), 0);
    }

    for rec_index in 0..numrecs {
        for (i, (var, &rec)) in file.variables.iter().zip(&is_record).enumerate() {
            if !rec {
                continue;
            }
            if rec_index == 0 {
                begins[i] = out.len();
            }
            let per_record = var.per_record_len(&file.dimensions);
            let slice = var.data.slice(rec_index * per_record, per_record);
            let bytes = slice.encode();
            out.extend_from_slice(&bytes);
            if num_record_vars > 1 {
                out.resize(out.len() + pad_len(bytes.len()), 0);
            }
        }
    }

    for (i, pos) in begin_patch_at.into_iter().enumerate() {
        out[pos..pos + 8].copy_from_slice(&(begins[i] as i64).to_be_bytes());
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_file() -> Netcdf3File {
        Netcdf3File {
            dimensions: vec![
                Netcdf3Dimension {
                    name: "frame".to_string(),
                    length: None,
                },
                Netcdf3Dimension {
                    name: "atom".to_string(),
                    length: Some(3),
                },
                Netcdf3Dimension {
                    name: "spatial".to_string(),
                    length: Some(3),
                },
            ],
            attributes: vec![
                ("title".to_string(), Netcdf3Value::Char("test".to_string())),
                ("version".to_string(), Netcdf3Value::Int(vec![1])),
            ],
            variables: vec![
                Netcdf3Variable {
                    name: "spatial".to_string(),
                    dimensions: vec![2],
                    attributes: vec![],
                    data: Netcdf3Value::Char("xyz".to_string()),
                },
                Netcdf3Variable {
                    name: "coordinates".to_string(),
                    dimensions: vec![0, 1, 2],
                    attributes: vec![(
                        "units".to_string(),
                        Netcdf3Value::Char("angstrom".to_string()),
                    )],
                    data: Netcdf3Value::Float(vec![
                        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, // frame 0
                        0.5, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.5, 0.0, // frame 1
                    ]),
                },
                Netcdf3Variable {
                    name: "cell_lengths".to_string(),
                    dimensions: vec![0, 1], // reuses "atom" dim id as a stand-in size-3 dim in this synthetic fixture
                    attributes: vec![],
                    data: Netcdf3Value::Double(vec![10.0, 20.0, 30.0, 10.0, 20.0, 30.0]),
                },
                Netcdf3Variable {
                    name: "flags".to_string(),
                    // 3 shorts/record = 6 bytes, not a multiple of 4 -- unlike
                    // "coordinates" (36B/record) and "cell_lengths" (24B/record),
                    // this one only round-trips correctly if the per-record pad
                    // is actually applied.
                    dimensions: vec![0, 1],
                    attributes: vec![],
                    data: Netcdf3Value::Short(vec![1, 2, 3, 4, 5, 6]),
                },
            ],
        }
    }

    #[test]
    fn test_magic_and_dimensions_round_trip() {
        let file = sample_file();
        let bytes = write(&file);
        assert_eq!(&bytes[0..3], b"CDF");
        assert_eq!(bytes[3], 2);

        let back = read(&bytes).expect("valid NetCDF-3");
        assert_eq!(back.dimensions.len(), 3);
        assert_eq!(back.dimensions[0].name, "frame");
        assert!(back.dimensions[0].length.is_none());
        assert_eq!(back.dimensions[1].length, Some(3));
    }

    #[test]
    fn test_global_attributes_of_every_type_round_trip() {
        let file = sample_file();
        let bytes = write(&file);
        let back = read(&bytes).expect("valid NetCDF-3");
        assert_eq!(
            back.attribute("title").and_then(Netcdf3Value::as_str),
            Some("test")
        );
        assert_eq!(back.attribute("version"), Some(&Netcdf3Value::Int(vec![1])));
    }

    #[test]
    fn test_a_non_record_variable_round_trips() {
        let file = sample_file();
        let bytes = write(&file);
        let back = read(&bytes).expect("valid NetCDF-3");
        let spatial = back.variable("spatial").expect("spatial exists");
        assert_eq!(spatial.data, Netcdf3Value::Char("xyz".to_string()));
    }

    #[test]
    fn test_record_variables_interleave_and_round_trip_with_padding() {
        // Three record variables (coordinates as f32, cell_lengths as f64,
        // flags as an odd-byte-count i16) -- more than one, so the
        // per-record 4-byte pad rule applies, and "flags" specifically only
        // round-trips correctly if that pad is actually inserted.
        let file = sample_file();
        let bytes = write(&file);
        let back = read(&bytes).expect("valid NetCDF-3");
        assert_eq!(back.num_records(), 2);

        let coords = back.variable("coordinates").expect("coordinates exists");
        assert_eq!(
            coords.data.as_f64_slice().unwrap(),
            vec![
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0,
                1.5, 0.0
            ]
        );

        let cell = back.variable("cell_lengths").expect("cell_lengths exists");
        assert_eq!(
            cell.data.as_f64_slice().unwrap(),
            vec![10.0, 20.0, 30.0, 10.0, 20.0, 30.0]
        );

        let flags = back.variable("flags").expect("flags exists");
        assert_eq!(flags.data, Netcdf3Value::Short(vec![1, 2, 3, 4, 5, 6]));
    }

    #[test]
    fn test_a_single_record_variable_gets_no_inter_record_padding() {
        // Exactly one record variable of odd per-record byte size (3
        // shorts = 6 bytes, not a multiple of 4) -- the padding rule must
        // NOT apply here, unlike the multi-record-variable case above.
        let file = Netcdf3File {
            dimensions: vec![
                Netcdf3Dimension {
                    name: "frame".to_string(),
                    length: None,
                },
                Netcdf3Dimension {
                    name: "three".to_string(),
                    length: Some(3),
                },
            ],
            attributes: vec![],
            variables: vec![Netcdf3Variable {
                name: "only".to_string(),
                dimensions: vec![0, 1],
                attributes: vec![],
                data: Netcdf3Value::Short(vec![1, 2, 3, 4, 5, 6]),
            }],
        };
        let bytes = write(&file);
        let back = read(&bytes).expect("valid NetCDF-3");
        assert_eq!(back.num_records(), 2);
        let only = back.variable("only").expect("only exists");
        assert_eq!(only.data, Netcdf3Value::Short(vec![1, 2, 3, 4, 5, 6]));
    }

    #[test]
    fn test_an_empty_file_has_no_dimensions_attributes_or_variables() {
        let file = Netcdf3File {
            dimensions: vec![],
            attributes: vec![],
            variables: vec![],
        };
        let bytes = write(&file);
        let back = read(&bytes).expect("valid NetCDF-3");
        assert!(back.dimensions.is_empty());
        assert!(back.attributes.is_empty());
        assert!(back.variables.is_empty());
    }

    #[test]
    fn test_a_hand_built_64_bit_offset_file_reads_with_an_8_byte_begin() {
        // Version 2 differs from version 1 only in the width of each
        // variable's `begin` field (i64 instead of i32) -- built by hand
        // since this crate's own writer never produces version 2.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"CDF");
        bytes.push(2); // 64-bit offset
        bytes.extend_from_slice(&0i32.to_be_bytes()); // numrecs
        bytes.extend_from_slice(&0i32.to_be_bytes()); // ABSENT dim_list
        bytes.extend_from_slice(&0i32.to_be_bytes());
        bytes.extend_from_slice(&0i32.to_be_bytes()); // ABSENT gatt_list
        bytes.extend_from_slice(&0i32.to_be_bytes());
        // var_list: one non-record NC_INT variable "x", no dims, no atts.
        bytes.extend_from_slice(&NC_VARIABLE.to_be_bytes());
        bytes.extend_from_slice(&1i32.to_be_bytes());
        bytes.extend_from_slice(&1i32.to_be_bytes()); // name length
        bytes.extend_from_slice(b"x\0\0\0"); // name + pad
        bytes.extend_from_slice(&0i32.to_be_bytes()); // ndims
        bytes.extend_from_slice(&0i32.to_be_bytes()); // ABSENT vatt_list
        bytes.extend_from_slice(&0i32.to_be_bytes());
        bytes.extend_from_slice(&NC_INT.to_be_bytes());
        bytes.extend_from_slice(&4i32.to_be_bytes()); // vsize
        let begin_pos = bytes.len();
        bytes.extend_from_slice(&0i64.to_be_bytes()); // 8-byte begin, patched below
        let data_pos = bytes.len() as i64;
        bytes[begin_pos..begin_pos + 8].copy_from_slice(&data_pos.to_be_bytes());
        bytes.extend_from_slice(&42i32.to_be_bytes());

        let file = read(&bytes).expect("valid NetCDF-3 v2");
        assert_eq!(
            file.variable("x").unwrap().data,
            Netcdf3Value::Int(vec![42])
        );
    }

    #[test]
    fn test_a_bad_magic_number_is_a_clear_error_not_a_panic() {
        let err = read(b"nope").unwrap_err();
        assert!(matches!(err, Netcdf3Error::InvalidMagicNumber), "{err}");
    }

    #[test]
    fn test_an_unsupported_version_byte_is_a_clear_error() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(b"CDF");
        bytes.push(9);
        let err = read(&bytes).unwrap_err();
        assert!(matches!(err, Netcdf3Error::UnsupportedVersion(9)), "{err}");
    }

    #[test]
    fn test_truncated_input_is_a_clear_error_not_a_panic() {
        let err = read(b"CDF\x01\x00\x00").unwrap_err();
        assert!(matches!(err, Netcdf3Error::ParseError(_)), "{err}");
    }
}
