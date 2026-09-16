//! MessagePack, hand-rolled rather than taken as a dependency (#319).
//!
//! Only the subset BinaryCIF actually emits: nil, bool, fixint/uint8-64/
//! int8-64, float32/64, fixstr/str8-32, bin8-32, fixarray/array16-32,
//! fixmap/map16-32. Deliberately **not** implemented: ext types, timestamps,
//! str64/array64/map64 — no real BinaryCIF producer emits any of these (a
//! CIF file's largest arrays are per-structure atom counts, nowhere near
//! `u32::MAX`), so an unsupported tag is a hard [`BcifError::UnsupportedTag`]
//! rather than a silent misread. Consistent with this crate's existing
//! style of hand-rolling a small, bounded surface rather than depending on
//! a general-purpose crate for it — see [`crate::io::format::Carries`]'s own
//! doc comment.
//!
//! All multi-byte integers/floats are big-endian, per the MessagePack spec
//! (independent of BinaryCIF's own little-endian [`crate::io::bcif`] byte
//! arrays, which are a different layer entirely).

use crate::io::errors::BcifError;

/// One decoded MessagePack value.
///
/// A tree rather than direct field extraction: BinaryCIF's descriptors are
/// heterogeneous maps whose shape depends on which encoding `"kind"` a
/// column uses, so decoding into a value first and then matching on it is
/// far simpler than threading a byte cursor through CIF-specific code.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Value {
    Nil,
    Bool(bool),
    Int(i64),
    UInt(u64),
    Float32(f32),
    Float64(f64),
    Str(String),
    Bin(Vec<u8>),
    Array(Vec<Value>),
    /// Insertion order preserved. BinaryCIF's own maps always have string
    /// keys in practice, but nothing here assumes that.
    Map(Vec<(Value, Value)>),
}

impl Value {
    pub(crate) fn str(s: impl Into<String>) -> Value {
        Value::Str(s.into())
    }

    /// Builds a [`Value::Map`] from `&str` keys, the shape every call site
    /// in [`crate::io::bcif`] wants.
    pub(crate) fn map(entries: Vec<(&str, Value)>) -> Value {
        Value::Map(
            entries
                .into_iter()
                .map(|(k, v)| (Value::Str(k.to_string()), v))
                .collect(),
        )
    }

    pub(crate) fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    pub(crate) fn as_array(&self) -> Option<&[Value]> {
        match self {
            Value::Array(items) => Some(items),
            _ => None,
        }
    }

    pub(crate) fn as_map(&self) -> Option<&[(Value, Value)]> {
        match self {
            Value::Map(entries) => Some(entries),
            _ => None,
        }
    }

    /// Looks up a string key in a [`Value::Map`]. `None` for any other
    /// variant, or a map without that key.
    pub(crate) fn map_get(&self, key: &str) -> Option<&Value> {
        self.as_map()?
            .iter()
            .find(|(k, _)| k.as_str() == Some(key))
            .map(|(_, v)| v)
    }

    pub(crate) fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Int(i) => Some(*i),
            Value::UInt(u) => i64::try_from(*u).ok(),
            _ => None,
        }
    }

    pub(crate) fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Float32(f) => Some(*f as f64),
            Value::Float64(f) => Some(*f),
            Value::Int(i) => Some(*i as f64),
            Value::UInt(u) => Some(*u as f64),
            _ => None,
        }
    }

    pub(crate) fn as_bytes(&self) -> Option<&[u8]> {
        match self {
            Value::Bin(b) => Some(b),
            _ => None,
        }
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], BcifError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.bytes.len())
            .ok_or_else(|| BcifError::InvalidMessagePack("unexpected end of input".to_string()))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn u8(&mut self) -> Result<u8, BcifError> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16, BcifError> {
        Ok(u16::from_be_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32, BcifError> {
        Ok(u32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64, BcifError> {
        Ok(u64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn i8(&mut self) -> Result<i8, BcifError> {
        Ok(self.take(1)?[0] as i8)
    }
    fn i16(&mut self) -> Result<i16, BcifError> {
        Ok(i16::from_be_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn i32(&mut self) -> Result<i32, BcifError> {
        Ok(i32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn i64(&mut self) -> Result<i64, BcifError> {
        Ok(i64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn f32(&mut self) -> Result<f32, BcifError> {
        Ok(f32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn f64(&mut self) -> Result<f64, BcifError> {
        Ok(f64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }
}

fn decode_str(r: &mut Reader, len: usize) -> Result<Value, BcifError> {
    let bytes = r.take(len)?;
    String::from_utf8(bytes.to_vec())
        .map(Value::Str)
        .map_err(|e| BcifError::InvalidMessagePack(e.to_string()))
}

fn decode_bin(r: &mut Reader, len: usize) -> Result<Value, BcifError> {
    Ok(Value::Bin(r.take(len)?.to_vec()))
}

fn decode_array(r: &mut Reader, len: usize) -> Result<Value, BcifError> {
    let mut out = Vec::with_capacity(len.min(1 << 20));
    for _ in 0..len {
        out.push(decode_value(r)?);
    }
    Ok(Value::Array(out))
}

fn decode_map(r: &mut Reader, len: usize) -> Result<Value, BcifError> {
    let mut out = Vec::with_capacity(len.min(1 << 20));
    for _ in 0..len {
        let k = decode_value(r)?;
        let v = decode_value(r)?;
        out.push((k, v));
    }
    Ok(Value::Map(out))
}

fn decode_value(r: &mut Reader) -> Result<Value, BcifError> {
    let tag = r.u8()?;
    match tag {
        0x00..=0x7f => Ok(Value::Int(tag as i64)),
        0xe0..=0xff => Ok(Value::Int(tag as i8 as i64)),
        0x80..=0x8f => decode_map(r, (tag & 0x0f) as usize),
        0x90..=0x9f => decode_array(r, (tag & 0x0f) as usize),
        0xa0..=0xbf => decode_str(r, (tag & 0x1f) as usize),
        0xc0 => Ok(Value::Nil),
        0xc2 => Ok(Value::Bool(false)),
        0xc3 => Ok(Value::Bool(true)),
        0xc4 => {
            let n = r.u8()? as usize;
            decode_bin(r, n)
        }
        0xc5 => {
            let n = r.u16()? as usize;
            decode_bin(r, n)
        }
        0xc6 => {
            let n = r.u32()? as usize;
            decode_bin(r, n)
        }
        0xca => Ok(Value::Float32(r.f32()?)),
        0xcb => Ok(Value::Float64(r.f64()?)),
        0xcc => Ok(Value::UInt(r.u8()? as u64)),
        0xcd => Ok(Value::UInt(r.u16()? as u64)),
        0xce => Ok(Value::UInt(r.u32()? as u64)),
        0xcf => Ok(Value::UInt(r.u64()?)),
        0xd0 => Ok(Value::Int(r.i8()? as i64)),
        0xd1 => Ok(Value::Int(r.i16()? as i64)),
        0xd2 => Ok(Value::Int(r.i32()? as i64)),
        0xd3 => Ok(Value::Int(r.i64()?)),
        0xd9 => {
            let n = r.u8()? as usize;
            decode_str(r, n)
        }
        0xda => {
            let n = r.u16()? as usize;
            decode_str(r, n)
        }
        0xdb => {
            let n = r.u32()? as usize;
            decode_str(r, n)
        }
        0xdc => {
            let n = r.u16()? as usize;
            decode_array(r, n)
        }
        0xdd => {
            let n = r.u32()? as usize;
            decode_array(r, n)
        }
        0xde => {
            let n = r.u16()? as usize;
            decode_map(r, n)
        }
        0xdf => {
            let n = r.u32()? as usize;
            decode_map(r, n)
        }
        other => Err(BcifError::UnsupportedTag(other)),
    }
}

/// Decodes one MessagePack value from the start of `bytes`. Trailing bytes
/// after the value are ignored — BinaryCIF's container is exactly one
/// top-level value, so there is never anything after it to check.
pub(crate) fn decode(bytes: &[u8]) -> Result<Value, BcifError> {
    let mut r = Reader::new(bytes);
    decode_value(&mut r)
}

fn encode_uint(u: u64, out: &mut Vec<u8>) {
    if u <= 0x7f {
        out.push(u as u8);
    } else if u <= u8::MAX as u64 {
        out.push(0xcc);
        out.push(u as u8);
    } else if u <= u16::MAX as u64 {
        out.push(0xcd);
        out.extend_from_slice(&(u as u16).to_be_bytes());
    } else if u <= u32::MAX as u64 {
        out.push(0xce);
        out.extend_from_slice(&(u as u32).to_be_bytes());
    } else {
        out.push(0xcf);
        out.extend_from_slice(&u.to_be_bytes());
    }
}

fn encode_int(i: i64, out: &mut Vec<u8>) {
    if i >= 0 {
        return encode_uint(i as u64, out);
    }
    if i >= -32 {
        out.push(i as u8);
    } else if i >= i8::MIN as i64 {
        out.push(0xd0);
        out.push(i as i8 as u8);
    } else if i >= i16::MIN as i64 {
        out.push(0xd1);
        out.extend_from_slice(&(i as i16).to_be_bytes());
    } else if i >= i32::MIN as i64 {
        out.push(0xd2);
        out.extend_from_slice(&(i as i32).to_be_bytes());
    } else {
        out.push(0xd3);
        out.extend_from_slice(&i.to_be_bytes());
    }
}

fn encode_str(s: &str, out: &mut Vec<u8>) {
    let bytes = s.as_bytes();
    let len = bytes.len();
    if len <= 31 {
        out.push(0xa0 | len as u8);
    } else if len <= u8::MAX as usize {
        out.push(0xd9);
        out.push(len as u8);
    } else if len <= u16::MAX as usize {
        out.push(0xda);
        out.extend_from_slice(&(len as u16).to_be_bytes());
    } else {
        out.push(0xdb);
        out.extend_from_slice(&(len as u32).to_be_bytes());
    }
    out.extend_from_slice(bytes);
}

fn encode_bin(b: &[u8], out: &mut Vec<u8>) {
    let len = b.len();
    if len <= u8::MAX as usize {
        out.push(0xc4);
        out.push(len as u8);
    } else if len <= u16::MAX as usize {
        out.push(0xc5);
        out.extend_from_slice(&(len as u16).to_be_bytes());
    } else {
        out.push(0xc6);
        out.extend_from_slice(&(len as u32).to_be_bytes());
    }
    out.extend_from_slice(b);
}

fn encode_array(items: &[Value], out: &mut Vec<u8>) {
    let len = items.len();
    if len <= 15 {
        out.push(0x90 | len as u8);
    } else if len <= u16::MAX as usize {
        out.push(0xdc);
        out.extend_from_slice(&(len as u16).to_be_bytes());
    } else {
        out.push(0xdd);
        out.extend_from_slice(&(len as u32).to_be_bytes());
    }
    for item in items {
        encode_into(item, out);
    }
}

fn encode_map(entries: &[(Value, Value)], out: &mut Vec<u8>) {
    let len = entries.len();
    if len <= 15 {
        out.push(0x80 | len as u8);
    } else if len <= u16::MAX as usize {
        out.push(0xde);
        out.extend_from_slice(&(len as u16).to_be_bytes());
    } else {
        out.push(0xdf);
        out.extend_from_slice(&(len as u32).to_be_bytes());
    }
    for (k, v) in entries {
        encode_into(k, out);
        encode_into(v, out);
    }
}

fn encode_into(value: &Value, out: &mut Vec<u8>) {
    match value {
        Value::Nil => out.push(0xc0),
        Value::Bool(false) => out.push(0xc2),
        Value::Bool(true) => out.push(0xc3),
        Value::Int(i) => encode_int(*i, out),
        Value::UInt(u) => encode_uint(*u, out),
        Value::Float32(f) => {
            out.push(0xca);
            out.extend_from_slice(&f.to_be_bytes());
        }
        Value::Float64(f) => {
            out.push(0xcb);
            out.extend_from_slice(&f.to_be_bytes());
        }
        Value::Str(s) => encode_str(s, out),
        Value::Bin(b) => encode_bin(b, out),
        Value::Array(items) => encode_array(items, out),
        Value::Map(entries) => encode_map(entries, out),
    }
}

pub(crate) fn encode(value: &Value) -> Vec<u8> {
    let mut out = Vec::new();
    encode_into(value, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trips(v: Value) {
        let bytes = encode(&v);
        let back = decode(&bytes).unwrap_or_else(|e| panic!("{v:?} failed to decode: {e}"));
        assert_eq!(back, v, "{bytes:x?}");
    }

    #[test]
    fn test_scalars_round_trip() {
        round_trips(Value::Nil);
        round_trips(Value::Bool(true));
        round_trips(Value::Bool(false));
        round_trips(Value::Float32(1.5));
        round_trips(Value::Float64(std::f64::consts::PI));
    }

    #[test]
    fn test_every_integer_width_round_trips() {
        // Exact `Value` identity only within the range that is unambiguous:
        // negative-fixint and the signed int8/16/32/64 tags always decode
        // back as `Int`.
        for i in [
            0i64,
            1,
            100,
            127,
            -1,
            -32,
            -33,
            i8::MIN as i64,
            i16::MIN as i64,
            i32::MIN as i64,
            i64::MIN,
        ] {
            round_trips(Value::Int(i));
        }

        // Beyond positive-fixint, the encoder picks the smallest tag
        // regardless of which variant a value started as -- a non-negative
        // `Int` above 127 comes back as `UInt`, since MessagePack's
        // uint8/16/32/64 tags carry no memory of "was originally signed".
        // Every call site in this crate reads an integer back through
        // `as_i64`, which treats both variants identically, so this checks
        // the *number* survives, not which variant names it.
        for i in [
            128i64,
            255,
            256,
            65535,
            65536,
            i16::MAX as i64,
            i32::MAX as i64,
            i64::MAX,
        ] {
            let bytes = encode(&Value::Int(i));
            let back = decode(&bytes).unwrap();
            assert_eq!(back.as_i64(), Some(i), "{bytes:x?}");
        }
        for u in [
            0u64,
            127,
            128,
            255,
            256,
            65535,
            65536,
            u32::MAX as u64,
            u64::MAX,
        ] {
            let bytes = encode(&Value::UInt(u));
            let back = decode(&bytes).unwrap();
            assert_eq!(back.as_i64(), i64::try_from(u).ok(), "{bytes:x?}");
        }
    }

    #[test]
    fn test_strings_of_every_length_class_round_trip() {
        for len in [0, 1, 31, 32, 255, 256, 70000] {
            round_trips(Value::str("x".repeat(len)));
        }
    }

    #[test]
    fn test_bin_of_every_length_class_round_trips() {
        for len in [0, 1, 255, 256, 70000] {
            round_trips(Value::Bin(vec![0xab; len]));
        }
    }

    #[test]
    fn test_arrays_and_maps_round_trip() {
        round_trips(Value::Array(vec![
            Value::Int(1),
            Value::str("two"),
            Value::Nil,
        ]));
        round_trips(Value::map(vec![
            ("a", Value::Int(1)),
            (
                "b",
                Value::Array(vec![Value::Bool(true), Value::Float64(2.5)]),
            ),
        ]));

        // Beyond `array16`'s 65535-element boundary; some of these values
        // exceed positive-fixint too, so this checks each element's number
        // rather than its exact `Value` variant (see
        // `test_every_integer_width_round_trips`).
        let big_array = Value::Array((0..70000).map(Value::Int).collect());
        let bytes = encode(&big_array);
        let back = decode(&bytes).unwrap();
        let back = back.as_array().unwrap();
        assert_eq!(back.len(), 70000);
        for (i, v) in back.iter().enumerate() {
            assert_eq!(v.as_i64(), Some(i as i64));
        }
    }

    #[test]
    fn test_map_get_finds_a_key_by_name() {
        let v = Value::map(vec![
            ("kind", Value::str("Delta")),
            ("origin", Value::Int(5)),
        ]);
        assert_eq!(v.map_get("kind").and_then(Value::as_str), Some("Delta"));
        assert_eq!(v.map_get("origin").and_then(Value::as_i64), Some(5));
        assert!(v.map_get("missing").is_none());
    }

    #[test]
    fn test_an_ext_type_tag_is_a_clear_error_not_a_panic() {
        // 0xd4 is fixext1 -- not implemented, and must fail loudly rather
        // than being silently misread as something else.
        let err = decode(&[0xd4, 0x01, 0x02]).unwrap_err();
        assert!(matches!(err, BcifError::UnsupportedTag(0xd4)), "{err}");
    }

    #[test]
    fn test_truncated_input_is_a_clear_error_not_a_panic() {
        // str8 claims a 10-byte string but only provides 2.
        let err = decode(&[0xd9, 10, b'h', b'i']).unwrap_err();
        assert!(matches!(err, BcifError::InvalidMessagePack(_)), "{err}");
    }
}
