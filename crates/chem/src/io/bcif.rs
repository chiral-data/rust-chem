//! BinaryCIF — the MessagePack encoding of the same CIF data model
//! [`crate::io::mmcif`] reads and writes as text (#319).
//!
//! The container is a MessagePack map: `{version, encoder, dataBlocks:
//! [{header, categories: [{name, rowCount, columns: [{name, data: {encoding,
//! data}, mask?}]}]}]}`. `data`/`mask`/`offsets` are raw byte blobs; the
//! `encoding` array names a chain of transforms to undo, in reverse array
//! order, before the bytes mean anything — `ByteArray` (a typed
//! reinterpretation of raw bytes, always last in the array, so always
//! decoded first), `FixedPoint`, `IntervalQuantization`, `RunLength`,
//! `Delta`, `IntegerPacking`, `StringArray`. Get one step wrong and the
//! result is not a crash, it is a plausible, wrong coordinate — the whole
//! reason this module's algorithms are taken from BinaryCIF's own reference
//! implementation rather than re-derived from the prose spec.
//!
//! Reuses [`crate::io::cif_model`]'s `build_molecule`/`build_rows` for
//! everything past "what are the `_atom_site` tags/rows and the `_cell`/
//! `_symmetry` singles" — this module's own job is only the container and
//! the seven encodings, not atom-site/chain/cell interpretation again.
//!
//! **Writer's own encoding choice** (the spec allows any legal pipeline per
//! column; this picks one, always lossless, not the smallest possible
//! file): coordinates/cell/occupancy/B-factor as `ByteArray(Float64)` only;
//! integer columns (`id`, `auth_seq_id`, `pdbx_PDB_model_num`) as `Delta` →
//! `IntegerPacking(Int8)` → `ByteArray(Int8)`; string columns as
//! `StringArray` over `Delta`-packed indices and offsets. `FixedPoint`/
//! `IntervalQuantization` are decode-only here — real producers (RCSB,
//! PDBe) commonly use `FixedPoint` for coordinates, so this reader must
//! decode both correctly, but this writer never produces either: nothing it
//! writes needs to be lossy.
//!
//! **No mask is ever written.** Every value in `cif_model::build_rows`'s
//! output is already a concrete string — missing values are the literal
//! CIF tokens `.`/`?`/`UNK`, the same way mmCIF text represents them — so
//! there is never a masked-absent slot to encode. A real foreign file may
//! still use one, so the *reader* honors it (mask `1`/`2` both map to
//! mmCIF's own `?`-means-absent convention).
//!
//! **Streams the whole file, and says so.** Unlike every line-oriented
//! supplier in [`crate::io::supplier`], MessagePack has no record boundary
//! to scan for — a document is valid only whole, the same situation
//! commonchem JSON is in (see [`crate::io::supplier::CommonchemSupplier`]'s
//! own doc comment). `BcifSupplier`/`BcifWriter` buffer accordingly.

use std::collections::HashMap;
use std::io::BufRead;

use crate::core::molecule::Molecule;
use crate::io::cif_model::{self, AtomSiteRows, RowPrecision};
use crate::io::errors::{BcifError, ReadError};
use crate::io::msgpack::{self, Value};
use crate::io::options::{ReadOptions, WriteOptions};
use crate::io::reader::{Payload, ReadOutcome, Record, Skipped};

/// One decoded column, at whatever stage the encoding chain has reached.
///
/// The chain always ends on `Ints`, `Floats` or `Strings` — `Bytes` only
/// ever appears as the *input* to a `ByteArray` step.
#[derive(Debug)]
enum Column {
    Bytes(Vec<u8>),
    Ints(Vec<i64>),
    Floats(Vec<f64>),
    Strings(Vec<String>),
}

impl Column {
    fn into_ints(self) -> Result<Vec<i64>, BcifError> {
        match self {
            Column::Ints(v) => Ok(v),
            _ => Err(BcifError::ParseError(
                "expected an integer column at this point in the encoding chain".to_string(),
            )),
        }
    }

    fn into_bytes(self) -> Result<Vec<u8>, BcifError> {
        match self {
            Column::Bytes(b) => Ok(b),
            _ => Err(BcifError::ParseError(
                "expected raw bytes at this point in the encoding chain".to_string(),
            )),
        }
    }

    /// Stringifies the fully-decoded column, full precision -- the string
    /// bridge into [`cif_model`] must never be mmCIF text's truncated
    /// `{:.3}`/`{:.2}` (#319).
    fn to_strings(&self) -> Vec<String> {
        match self {
            Column::Ints(v) => v.iter().map(i64::to_string).collect(),
            Column::Floats(v) => v.iter().map(f64::to_string).collect(),
            Column::Strings(v) => v.clone(),
            Column::Bytes(_) => {
                unreachable!("a column's encoding chain must not end on raw bytes")
            }
        }
    }
}

fn le_chunks<T>(
    bytes: &[u8],
    width: usize,
    from_le: impl Fn(&[u8]) -> T,
) -> Result<Vec<T>, BcifError> {
    if !bytes.len().is_multiple_of(width) {
        return Err(BcifError::ParseError(format!(
            "{} bytes is not a multiple of the element width {width}",
            bytes.len()
        )));
    }
    Ok(bytes.chunks_exact(width).map(from_le).collect())
}

fn decode_byte_array(bytes: Vec<u8>, type_code: i64) -> Result<Column, BcifError> {
    match type_code {
        1 => Ok(Column::Ints(
            bytes.iter().map(|&b| b as i8 as i64).collect(),
        )),
        2 => Ok(Column::Ints(le_chunks(&bytes, 2, |c| {
            i16::from_le_bytes(c.try_into().unwrap()) as i64
        })?)),
        3 => Ok(Column::Ints(le_chunks(&bytes, 4, |c| {
            i32::from_le_bytes(c.try_into().unwrap()) as i64
        })?)),
        4 => Ok(Column::Ints(bytes.iter().map(|&b| b as i64).collect())),
        5 => Ok(Column::Ints(le_chunks(&bytes, 2, |c| {
            u16::from_le_bytes(c.try_into().unwrap()) as i64
        })?)),
        6 => Ok(Column::Ints(le_chunks(&bytes, 4, |c| {
            u32::from_le_bytes(c.try_into().unwrap()) as i64
        })?)),
        32 => Ok(Column::Floats(le_chunks(&bytes, 4, |c| {
            f32::from_le_bytes(c.try_into().unwrap()) as f64
        })?)),
        33 => Ok(Column::Floats(le_chunks(&bytes, 8, |c| {
            f64::from_le_bytes(c.try_into().unwrap())
        })?)),
        other => Err(BcifError::ParseError(format!(
            "unknown ByteArray type code {other}"
        ))),
    }
}

fn require_i64(step: &Value, key: &str) -> Result<i64, BcifError> {
    step.map_get(key)
        .and_then(Value::as_i64)
        .ok_or_else(|| BcifError::ParseError(format!("encoding step has no integer '{key}'")))
}

fn require_f64(step: &Value, key: &str) -> Result<f64, BcifError> {
    step.map_get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| BcifError::ParseError(format!("encoding step has no number '{key}'")))
}

fn decode_step(current: Column, step: &Value) -> Result<Column, BcifError> {
    let kind = step
        .map_get("kind")
        .and_then(Value::as_str)
        .ok_or_else(|| BcifError::ParseError("encoding step has no 'kind'".to_string()))?;

    match kind {
        "ByteArray" => decode_byte_array(current.into_bytes()?, require_i64(step, "type")?),

        "FixedPoint" => {
            let ints = current.into_ints()?;
            let factor = require_f64(step, "factor")?;
            Ok(Column::Floats(
                ints.iter().map(|&i| i as f64 / factor).collect(),
            ))
        }

        "IntervalQuantization" => {
            let ints = current.into_ints()?;
            let min = require_f64(step, "min")?;
            let max = require_f64(step, "max")?;
            let num_steps = require_i64(step, "numSteps")?;
            let delta = (max - min) / (num_steps as f64 - 1.0);
            Ok(Column::Floats(
                ints.iter().map(|&i| min + delta * i as f64).collect(),
            ))
        }

        "RunLength" => {
            let ints = current.into_ints()?;
            let src_size = require_i64(step, "srcSize")? as usize;
            let mut out = Vec::with_capacity(src_size);
            for pair in ints.chunks(2) {
                let [value, count] = pair else {
                    return Err(BcifError::ParseError(
                        "RunLength data has an odd number of elements".to_string(),
                    ));
                };
                for _ in 0..*count {
                    out.push(*value);
                }
            }
            Ok(Column::Ints(out))
        }

        "Delta" => {
            let ints = current.into_ints()?;
            let origin = step.map_get("origin").and_then(Value::as_i64).unwrap_or(0);
            let mut out = Vec::with_capacity(ints.len());
            if let Some(&first) = ints.first() {
                out.push(first + origin);
                for &d in &ints[1..] {
                    out.push(d + out[out.len() - 1]);
                }
            }
            Ok(Column::Ints(out))
        }

        "IntegerPacking" => {
            let packed = current.into_ints()?;
            let byte_count = require_i64(step, "byteCount")?;
            let is_unsigned = matches!(step.map_get("isUnsigned"), Some(Value::Bool(true)));
            let src_size = require_i64(step, "srcSize")? as usize;

            if packed.len() == src_size {
                return Ok(Column::Ints(packed));
            }

            let upper = if byte_count == 1 {
                if is_unsigned { 0xFF } else { 0x7F }
            } else if is_unsigned {
                0xFFFF
            } else {
                0x7FFF
            };
            let lower = if is_unsigned { i64::MIN } else { -upper - 1 };

            let mut out = Vec::with_capacity(src_size);
            let mut i = 0;
            while i < packed.len() {
                let mut value = 0i64;
                let mut t = packed[i];
                while t == upper || (!is_unsigned && t == lower) {
                    value += t;
                    i += 1;
                    t = *packed.get(i).ok_or_else(|| {
                        BcifError::ParseError("IntegerPacking ran past its input".to_string())
                    })?;
                }
                value += t;
                out.push(value);
                i += 1;
            }
            Ok(Column::Ints(out))
        }

        "StringArray" => {
            let raw = current.into_bytes()?;
            let data_encoding = step
                .map_get("dataEncoding")
                .and_then(Value::as_array)
                .ok_or_else(|| {
                    BcifError::ParseError("StringArray has no dataEncoding".to_string())
                })?;
            let mut indices_col = Column::Bytes(raw);
            for s in data_encoding.iter().rev() {
                indices_col = decode_step(indices_col, s)?;
            }
            let indices = indices_col.into_ints()?;

            let offset_bytes = step
                .map_get("offsets")
                .and_then(Value::as_bytes)
                .ok_or_else(|| BcifError::ParseError("StringArray has no offsets".to_string()))?;
            let offset_encoding = step
                .map_get("offsetEncoding")
                .and_then(Value::as_array)
                .ok_or_else(|| {
                    BcifError::ParseError("StringArray has no offsetEncoding".to_string())
                })?;
            let mut offsets_col = Column::Bytes(offset_bytes.to_vec());
            for s in offset_encoding.iter().rev() {
                offsets_col = decode_step(offsets_col, s)?;
            }
            let offsets = offsets_col.into_ints()?;

            let string_data = step
                .map_get("stringData")
                .and_then(Value::as_str)
                .ok_or_else(|| {
                    BcifError::ParseError("StringArray has no stringData".to_string())
                })?;

            let mut table = vec![String::new()];
            for w in 1..offsets.len() {
                let start = offsets[w - 1].max(0) as usize;
                let end = offsets[w].max(0) as usize;
                let slice = string_data.get(start..end).ok_or_else(|| {
                    BcifError::ParseError("StringArray offsets run past stringData".to_string())
                })?;
                table.push(slice.to_string());
            }

            let mut out = Vec::with_capacity(indices.len());
            for idx in indices {
                let position = usize::try_from(idx + 1).map_err(|_| {
                    BcifError::ParseError(format!("StringArray index {idx} is out of range"))
                })?;
                let s = table.get(position).ok_or_else(|| {
                    BcifError::ParseError(format!("StringArray index {idx} is out of range"))
                })?;
                out.push(s.clone());
            }
            Ok(Column::Strings(out))
        }

        other => Err(BcifError::UnknownEncoding(other.to_string())),
    }
}

fn decode_encoded_data(data: &Value) -> Result<Column, BcifError> {
    let encoding = data
        .map_get("encoding")
        .and_then(Value::as_array)
        .ok_or_else(|| BcifError::ParseError("column data has no encoding array".to_string()))?;
    let raw = data
        .map_get("data")
        .and_then(Value::as_bytes)
        .ok_or_else(|| BcifError::ParseError("column data has no data bytes".to_string()))?;

    let mut current = Column::Bytes(raw.to_vec());
    for step in encoding.iter().rev() {
        current = decode_step(current, step)?;
    }
    Ok(current)
}

/// Decodes one column into one string per row, `None` where the mask (if
/// any) says the value is not specified (`1`) or unknown (`2`) — both
/// collapse into the one "absent" state [`cif_model`] already knows, the
/// same way mmCIF text's `.`/`?` do.
fn decode_column(col: &Value) -> Result<Vec<Option<String>>, BcifError> {
    let data = col
        .map_get("data")
        .ok_or_else(|| BcifError::ParseError("column has no data".to_string()))?;
    let values = decode_encoded_data(data)?.to_strings();

    let mask = match col.map_get("mask") {
        Some(m) if !matches!(m, Value::Nil) => Some(decode_encoded_data(m)?.into_ints()?),
        _ => None,
    };

    Ok(values
        .into_iter()
        .enumerate()
        .map(|(i, v)| {
            let absent = mask
                .as_ref()
                .and_then(|m| m.get(i))
                .is_some_and(|&c| c != 0);
            if absent { None } else { Some(v) }
        })
        .collect())
}

/// Decodes a whole BinaryCIF document into one `(name, Molecule)` per data
/// block.
pub(crate) fn decode_bcif(bytes: &[u8]) -> Result<Vec<(String, Molecule)>, BcifError> {
    let top = msgpack::decode(bytes)?;
    let blocks = top
        .map_get("dataBlocks")
        .and_then(Value::as_array)
        .ok_or_else(|| BcifError::ParseError("BinaryCIF file has no dataBlocks".to_string()))?;

    let mut out = Vec::with_capacity(blocks.len());
    for (i, block) in blocks.iter().enumerate() {
        let header = block
            .map_get("header")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_string();
        let categories = block
            .map_get("categories")
            .and_then(Value::as_array)
            .unwrap_or(&[]);

        let mut singles: HashMap<String, String> = HashMap::new();
        let mut atom_site: Option<AtomSiteRows> = None;

        for category in categories {
            let name = category
                .map_get("name")
                .and_then(Value::as_str)
                .ok_or_else(|| BcifError::ParseError("category has no name".to_string()))?;
            let row_count = category
                .map_get("rowCount")
                .and_then(Value::as_i64)
                .unwrap_or(0) as usize;
            let columns = category
                .map_get("columns")
                .and_then(Value::as_array)
                .unwrap_or(&[]);

            if name == "_atom_site" {
                let mut tags = Vec::with_capacity(columns.len());
                let mut per_column: Vec<Vec<Option<String>>> = Vec::with_capacity(columns.len());
                for col in columns {
                    let col_name = col
                        .map_get("name")
                        .and_then(Value::as_str)
                        .ok_or_else(|| BcifError::ParseError("column has no name".to_string()))?;
                    tags.push(format!("_atom_site.{col_name}"));
                    per_column.push(decode_column(col)?);
                }
                let mut rows = Vec::with_capacity(row_count);
                for r in 0..row_count {
                    let row: Vec<String> = per_column
                        .iter()
                        .map(|values| {
                            values
                                .get(r)
                                .cloned()
                                .flatten()
                                .unwrap_or_else(|| "?".to_string())
                        })
                        .collect();
                    rows.push(row);
                }
                atom_site = Some((tags, rows));
            } else {
                for col in columns {
                    let col_name = col
                        .map_get("name")
                        .and_then(Value::as_str)
                        .ok_or_else(|| BcifError::ParseError("column has no name".to_string()))?;
                    let values = decode_column(col)?;
                    if let Some(Some(first)) = values.first() {
                        singles.insert(format!("{name}.{col_name}"), first.clone());
                    }
                }
            }
        }

        let molecule = cif_model::build_molecule(
            &singles,
            atom_site
                .as_ref()
                .map(|(t, r)| (t.as_slice(), r.as_slice())),
        )?;
        let name = if header.is_empty() {
            format!("Molecule_{}", i + 1)
        } else {
            header
        };
        out.push((name, molecule));
    }
    Ok(out)
}

// --- encoding (writer direction) -------------------------------------------

fn encode_float64_steps(values: &[f64]) -> (Vec<Value>, Vec<u8>) {
    let mut bytes = Vec::with_capacity(values.len() * 8);
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    (
        vec![Value::map(vec![
            ("kind", Value::str("ByteArray")),
            ("type", Value::Int(33)),
        ])],
        bytes,
    )
}

/// Packs `values` into signed `Int8` bytes, expanding any value outside
/// `[-128, 127]` into a run of `127`/`-128` sentinels plus a final,
/// unambiguous remainder — the exact inverse of `decode_step`'s
/// `IntegerPacking` branch. A value equal to a sentinel is *also* expanded
/// (never emitted as a lone byte): the decoder cannot otherwise tell "this
/// is the value" from "this starts an overflow run".
fn pack_i8(values: &[i64]) -> Vec<i8> {
    const UPPER: i64 = 0x7F;
    const LOWER: i64 = -UPPER - 1;

    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        if v > LOWER && v < UPPER {
            out.push(v as i8);
            continue;
        }
        let mut remaining = v;
        if remaining >= UPPER {
            while remaining >= UPPER {
                out.push(UPPER as i8);
                remaining -= UPPER;
            }
        } else {
            while remaining <= LOWER {
                out.push(LOWER as i8);
                remaining -= LOWER;
            }
        }
        out.push(remaining as i8);
    }
    out
}

fn encode_packed_int_steps(values: &[i64]) -> (Vec<Value>, Vec<u8>) {
    let n = values.len();
    let origin = values.first().copied().unwrap_or(0);
    let mut deltas = Vec::with_capacity(n);
    if let Some(&first) = values.first() {
        deltas.push(first - origin);
        for w in 1..n {
            deltas.push(values[w] - values[w - 1]);
        }
    }
    let packed = pack_i8(&deltas);
    let bytes: Vec<u8> = packed.iter().map(|&b| b as u8).collect();
    let steps = vec![
        Value::map(vec![
            ("kind", Value::str("Delta")),
            ("origin", Value::Int(origin)),
            ("srcType", Value::Int(3)),
        ]),
        Value::map(vec![
            ("kind", Value::str("IntegerPacking")),
            ("byteCount", Value::Int(1)),
            ("isUnsigned", Value::Bool(false)),
            ("srcSize", Value::Int(n as i64)),
        ]),
        Value::map(vec![
            ("kind", Value::str("ByteArray")),
            ("type", Value::Int(1)),
        ]),
    ];
    (steps, bytes)
}

/// `offsets` are monotonically non-decreasing and there is one per unique
/// string rather than one per atom, so plain `Delta` → `ByteArray(Int32)`
/// (no packing) is simple and small enough.
fn encode_offsets_steps(offsets: &[i64]) -> (Vec<Value>, Vec<u8>) {
    let n = offsets.len();
    let origin = offsets.first().copied().unwrap_or(0);
    let mut deltas = Vec::with_capacity(n);
    if let Some(&first) = offsets.first() {
        deltas.push(first - origin);
        for w in 1..n {
            deltas.push(offsets[w] - offsets[w - 1]);
        }
    }
    let mut bytes = Vec::with_capacity(deltas.len() * 4);
    for d in &deltas {
        bytes.extend_from_slice(&(*d as i32).to_le_bytes());
    }
    let steps = vec![
        Value::map(vec![
            ("kind", Value::str("Delta")),
            ("origin", Value::Int(origin)),
            ("srcType", Value::Int(3)),
        ]),
        Value::map(vec![
            ("kind", Value::str("ByteArray")),
            ("type", Value::Int(3)),
        ]),
    ];
    (steps, bytes)
}

fn encode_string_steps(values: &[String]) -> (Vec<Value>, Vec<u8>) {
    let mut table: Vec<&str> = Vec::new();
    let mut lookup: HashMap<&str, i64> = HashMap::new();
    let mut indices = Vec::with_capacity(values.len());
    for v in values {
        if v.is_empty() {
            indices.push(-1);
            continue;
        }
        let idx = *lookup.entry(v.as_str()).or_insert_with(|| {
            table.push(v.as_str());
            (table.len() - 1) as i64
        });
        indices.push(idx);
    }

    let mut string_data = String::new();
    let mut offsets: Vec<i64> = vec![0];
    for s in &table {
        string_data.push_str(s);
        offsets.push(string_data.len() as i64);
    }

    let (data_steps, data_bytes) = encode_packed_int_steps(&indices);
    let (offset_steps, offset_bytes) = encode_offsets_steps(&offsets);

    let step = Value::map(vec![
        ("kind", Value::str("StringArray")),
        ("dataEncoding", Value::Array(data_steps)),
        ("stringData", Value::str(string_data)),
        ("offsetEncoding", Value::Array(offset_steps)),
        ("offsets", Value::Bin(offset_bytes)),
    ]);
    (vec![step], data_bytes)
}

enum ColumnEncoding {
    Float,
    Int,
    Str,
}

fn atom_site_column_encoding(name: &str) -> ColumnEncoding {
    match name {
        "Cartn_x" | "Cartn_y" | "Cartn_z" | "occupancy" | "B_iso_or_equiv" => ColumnEncoding::Float,
        "id" | "auth_seq_id" | "pdbx_PDB_model_num" => ColumnEncoding::Int,
        _ => ColumnEncoding::Str,
    }
}

fn encode_column(name: &str, values: &[String], encoding: ColumnEncoding) -> Value {
    let (steps, bytes) = match encoding {
        ColumnEncoding::Float => {
            let floats: Vec<f64> = values.iter().map(|s| s.parse().unwrap_or(0.0)).collect();
            encode_float64_steps(&floats)
        }
        ColumnEncoding::Int => {
            let ints: Vec<i64> = values.iter().map(|s| s.parse().unwrap_or(0)).collect();
            encode_packed_int_steps(&ints)
        }
        ColumnEncoding::Str => encode_string_steps(values),
    };
    Value::map(vec![
        ("name", Value::str(name)),
        (
            "data",
            Value::map(vec![
                ("encoding", Value::Array(steps)),
                ("data", Value::Bin(bytes)),
            ]),
        ),
    ])
}

/// Encodes named molecules into one BinaryCIF document, one data block per
/// record.
pub(crate) fn encode_bcif(records: &[(String, Molecule)]) -> Vec<u8> {
    let mut blocks = Vec::with_capacity(records.len());
    for (name, mol) in records {
        let (singles, (tags, rows)) = cif_model::build_rows(mol, RowPrecision::Full);
        let mut categories = Vec::new();

        if singles.contains_key("_cell.length_a") {
            let cell_cols: Vec<Value> = [
                ("length_a", "_cell.length_a"),
                ("length_b", "_cell.length_b"),
                ("length_c", "_cell.length_c"),
                ("angle_alpha", "_cell.angle_alpha"),
                ("angle_beta", "_cell.angle_beta"),
                ("angle_gamma", "_cell.angle_gamma"),
            ]
            .iter()
            .map(|(col, key)| {
                encode_column(
                    col,
                    std::slice::from_ref(&singles[*key]),
                    ColumnEncoding::Float,
                )
            })
            .collect();
            categories.push(Value::map(vec![
                ("name", Value::str("_cell")),
                ("rowCount", Value::Int(1)),
                ("columns", Value::Array(cell_cols)),
            ]));
        }
        if let Some(sg) = singles.get("_symmetry.space_group_name_H-M") {
            categories.push(Value::map(vec![
                ("name", Value::str("_symmetry")),
                ("rowCount", Value::Int(1)),
                (
                    "columns",
                    Value::Array(vec![encode_column(
                        "space_group_name_H-M",
                        std::slice::from_ref(sg),
                        ColumnEncoding::Str,
                    )]),
                ),
            ]));
        }

        let atom_site_columns: Vec<Value> = tags
            .iter()
            .enumerate()
            .map(|(ci, tag)| {
                let col_name = tag.strip_prefix("_atom_site.").unwrap_or(tag);
                let values: Vec<String> = rows.iter().map(|r| r[ci].clone()).collect();
                encode_column(col_name, &values, atom_site_column_encoding(col_name))
            })
            .collect();
        categories.push(Value::map(vec![
            ("name", Value::str("_atom_site")),
            ("rowCount", Value::Int(rows.len() as i64)),
            ("columns", Value::Array(atom_site_columns)),
        ]));

        blocks.push(Value::map(vec![
            ("header", Value::str(name.as_str())),
            ("categories", Value::Array(categories)),
        ]));
    }

    let top = Value::map(vec![
        ("version", Value::str("0.3.0")),
        ("encoder", Value::str("chem")),
        ("dataBlocks", Value::Array(blocks)),
    ]);
    msgpack::encode(&top)
}

// --- registry entry points --------------------------------------------------

/// [`crate::io::format::ByteReadFn`] for BinaryCIF.
pub(crate) fn read_bcif_with_options(bytes: &[u8], _options: &ReadOptions) -> ReadOutcome {
    let mut out = ReadOutcome::default();
    match decode_bcif(bytes) {
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
            input: String::new(),
            error: e.to_string(),
        }),
    }
    out
}

/// [`crate::io::format::ByteWriteFn`] for BinaryCIF.
pub(crate) fn write_bcif_records(
    records: &[(String, Molecule)],
    _options: &WriteOptions,
) -> Vec<u8> {
    encode_bcif(records)
}

/// Buffers every molecule and decodes once the stream ends -- MessagePack
/// has no record boundary to scan for one at a time, the same situation
/// [`crate::io::supplier::CommonchemSupplier`] is in.
pub struct BcifSupplier {
    records: std::vec::IntoIter<Result<Record, ReadError>>,
}

impl BcifSupplier {
    pub fn new<R: BufRead>(mut reader: R, _options: &ReadOptions) -> Self {
        let mut bytes = Vec::new();
        let records = match reader.read_to_end(&mut bytes) {
            Err(source) => vec![Err(ReadError::Io {
                position: 1,
                source,
            })],
            Ok(_) => match decode_bcif(&bytes) {
                Ok(molecules) => molecules
                    .into_iter()
                    .map(|(name, molecule)| {
                        Ok(Record {
                            payload: Payload::Molecule(molecule),
                            name,
                            smiles: None,
                        })
                    })
                    .collect(),
                Err(e) => vec![Err(ReadError::Parse {
                    position: 1,
                    message: e.to_string(),
                })],
            },
        };
        Self {
            records: records.into_iter(),
        }
    }
}

impl Iterator for BcifSupplier {
    type Item = Result<Record, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next()
    }
}

/// Buffers every molecule and emits one document in
/// [`crate::io::supplier::Writer::finish`] -- the column encodings need a
/// whole column before any of them can be written, so there is no
/// incremental path even in principle.
pub struct BcifWriter<W> {
    writer: W,
    records: Vec<(String, Molecule)>,
}

impl<W: std::io::Write> BcifWriter<W> {
    pub fn new(writer: W, _options: &WriteOptions) -> Self {
        Self {
            writer,
            records: Vec::new(),
        }
    }
}

impl<W: std::io::Write> crate::io::supplier::Writer for BcifWriter<W> {
    fn write_molecule(&mut self, name: &str, molecule: &Molecule) -> std::io::Result<()> {
        self.records.push((name.to_string(), molecule.clone()));
        Ok(())
    }

    fn finish(mut self: Box<Self>) -> std::io::Result<()> {
        let bytes = encode_bcif(&self.records);
        self.writer.write_all(&bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{Atom, Element};
    use crate::core::cell::UnitCell;
    use crate::core::geometry::Point3;

    fn ethanol_like() -> Molecule {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.set_coords3(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(2.25, 1.3, 0.0),
        ])
        .unwrap();
        mol
    }

    #[test]
    fn test_a_molecule_round_trips_through_encode_and_decode() {
        let mol = ethanol_like();
        let bytes = encode_bcif(&[("ethanol".to_string(), mol.clone())]);
        let back = decode_bcif(&bytes).expect("valid BinaryCIF");
        assert_eq!(back.len(), 1);
        assert_eq!(back[0].0, "ethanol");
        assert_eq!(back[0].1.num_atoms(), 3);
        for i in 0..3 {
            assert_eq!(back[0].1.coord3(i), mol.coord3(i));
        }
    }

    #[test]
    fn test_a_coordinate_survives_full_precision_not_mmcif_texts_three_decimals() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.set_coords3(vec![Point3::new(1.234_567_89, 0.0, 0.0)])
            .unwrap();

        let bytes = encode_bcif(&[("x".to_string(), mol)]);
        let back = decode_bcif(&bytes).unwrap();
        assert_eq!(back[0].1.coord3(0).unwrap().x, 1.234_567_89);
    }

    #[test]
    fn test_cell_and_space_group_round_trip() {
        let mut mol = ethanol_like();
        mol.set_cell(UnitCell::new(10.0, 20.0, 30.0, 90.0, 90.0, 90.0))
            .unwrap();
        mol.set_space_group(crate::core::cell::SpaceGroup::from_symbol("P 21 21 21"));

        let bytes = encode_bcif(&[("x".to_string(), mol.clone())]);
        let back = decode_bcif(&bytes).unwrap();
        assert_eq!(back[0].1.cell(), mol.cell());
        assert_eq!(
            back[0].1.space_group().and_then(|g| g.symbol.as_deref()),
            Some("P 21 21 21")
        );
    }

    #[test]
    fn test_a_large_chain_boundary_delta_round_trips() {
        // auth_seq_id resetting from a large number back to 1 at a chain
        // boundary is a delta well outside Int8's range -- exactly the case
        // `pack_i8`'s sentinel-expansion exists for.
        let mut mol = Molecule::new();
        for _ in 0..5 {
            mol.add_atom(Atom::new(Element::carbon()));
        }
        mol.set_coords3(vec![Point3::new(0.0, 0.0, 0.0); 5])
            .unwrap();
        mol.set_topology(
            vec![
                crate::core::residue::Chain {
                    id: "A".to_string(),
                    label_id: None,
                    residues: 0..1,
                },
                crate::core::residue::Chain {
                    id: "B".to_string(),
                    label_id: None,
                    residues: 1..2,
                },
            ],
            vec![
                crate::core::residue::Residue {
                    name: "UNK".to_string(),
                    sequence: 300,
                    insertion_code: None,
                    label_seq: None,
                    chain_ix: 0,
                    is_hetero: false,
                    atoms: 0..3,
                },
                crate::core::residue::Residue {
                    name: "UNK".to_string(),
                    sequence: 1,
                    insertion_code: None,
                    label_seq: None,
                    chain_ix: 1,
                    is_hetero: false,
                    atoms: 3..5,
                },
            ],
        )
        .unwrap();

        let bytes = encode_bcif(&[("x".to_string(), mol.clone())]);
        let back = decode_bcif(&bytes).unwrap();
        assert_eq!(back[0].1.residues()[0].sequence, 300);
    }

    #[test]
    fn test_pack_i8_matches_the_specs_own_worked_example() {
        // From the BinaryCIF spec: decoding byteCount=1 packed bytes
        // [1, 2, -3, 127, 1] yields [1, 2, -3, 128]. `pack_i8` is the
        // encode direction, so it must produce exactly that packed form
        // for that input -- not merely *a* valid encoding.
        assert_eq!(pack_i8(&[1, 2, -3, 128]), vec![1, 2, -3, 127, 1]);
    }

    #[test]
    fn test_delta_decode_matches_the_specs_own_worked_example() {
        let step = Value::map(vec![
            ("kind", Value::str("Delta")),
            ("origin", Value::Int(1000)),
            ("srcType", Value::Int(3)),
        ]);
        let out = decode_step(Column::Ints(vec![0, 3, 2, 1]), &step).unwrap();
        assert_eq!(out.into_ints().unwrap(), vec![1000, 1003, 1005, 1006]);
    }

    #[test]
    fn test_fixed_point_decode_matches_the_specs_own_worked_example() {
        // The spec's own encode-side example rounds 0.123 * 100 to the
        // integer 12 -- FixedPoint is lossy by construction, so decoding
        // that same 12 back gives 0.12, not 0.123. This test is the decode
        // half only: 12 / 100 = 0.12, exactly.
        let step = Value::map(vec![
            ("kind", Value::str("FixedPoint")),
            ("factor", Value::Int(100)),
            ("srcType", Value::Int(33)),
        ]);
        let out = decode_step(Column::Ints(vec![120, 123, 12]), &step).unwrap();
        match out {
            Column::Floats(v) => {
                assert!((v[0] - 1.2).abs() < 1e-9);
                assert!((v[1] - 1.23).abs() < 1e-9);
                assert!((v[2] - 0.12).abs() < 1e-9);
            }
            _ => panic!("expected Floats"),
        }
    }

    #[test]
    fn test_an_unknown_encoding_kind_is_a_clear_error_not_a_panic() {
        let step = Value::map(vec![("kind", Value::str("NotARealEncoding"))]);
        let err = decode_step(Column::Ints(vec![1]), &step).unwrap_err();
        assert!(matches!(err, BcifError::UnknownEncoding(_)), "{err}");
    }

    #[test]
    fn test_a_document_with_no_data_blocks_key_is_a_clear_error() {
        let bytes = msgpack::encode(&Value::map(vec![("version", Value::str("0.3.0"))]));
        let err = decode_bcif(&bytes).unwrap_err();
        assert!(matches!(err, BcifError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_the_corpus_fixtures_this_crate_wrote_decode_correctly() {
        // Generated by this crate's own CLI (`chem convert *.cif --to bcif`)
        // from the existing, already-ours mmCIF fixtures -- not vendored,
        // and not hand-crafted bytes, so this proves the encoder and decoder
        // agree with each other on real structure-sized input, not just the
        // small hand-built molecules the other tests in this module use.
        let dipeptide = decode_bcif(include_bytes!(
            "../../tests/corpus/bcif/dipeptide-with-ligand.bcif"
        ))
        .expect("valid BinaryCIF");
        assert_eq!(dipeptide.len(), 1);
        assert_eq!(dipeptide[0].1.num_atoms(), 8);
        assert!(dipeptide[0].1.has_coords3());

        let water = decode_bcif(include_bytes!("../../tests/corpus/bcif/water-no-cell.bcif"))
            .expect("valid BinaryCIF");
        assert_eq!(water.len(), 1);
        assert_eq!(water[0].1.num_atoms(), 3);
        assert!(water[0].1.cell().is_none());
    }
}
