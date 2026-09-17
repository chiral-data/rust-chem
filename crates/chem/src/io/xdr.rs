//! XDR (RFC 1014), hand-rolled rather than taken as a dependency.
//!
//! Only the subset TRR (and later XTC) actually need: signed/unsigned 32-bit
//! integers, IEEE single/double floats, and length-prefixed, zero-padded
//! opaque strings. Big-endian throughout, per the spec — independent of
//! [`crate::io::msgpack`]'s own big-endian encoding, which is a different
//! format entirely. Consistent with this crate's existing style of
//! hand-rolling a small, bounded surface rather than depending on a
//! general-purpose crate for it — see [`crate::io::format::Carries`]'s own
//! doc comment.
//!
//! Every primitive here is already a multiple of 4 bytes — an `i32`/`f32` is
//! exactly 4, an `f64` is exactly 8 (already a multiple of 4, so no padding
//! quirk exists for it) — so only this module's own `read_string`/
//! `write_string` need the length-prefix-plus-padding logic RFC 1014
//! describes for opaque data.

use crate::io::errors::TrrError;

pub(crate) struct XdrReader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> XdrReader<'a> {
    pub(crate) fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    pub(crate) fn position(&self) -> usize {
        self.pos
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], TrrError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.bytes.len())
            .ok_or_else(|| TrrError::ParseError("unexpected end of input".to_string()))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    pub(crate) fn read_i32(&mut self) -> Result<i32, TrrError> {
        Ok(i32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub(crate) fn read_u32(&mut self) -> Result<u32, TrrError> {
        Ok(u32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub(crate) fn read_f32(&mut self) -> Result<f32, TrrError> {
        Ok(f32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub(crate) fn read_f64(&mut self) -> Result<f64, TrrError> {
        Ok(f64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }

    /// RFC 1014's `xdr_string`: a length prefix (the string's own length,
    /// no NUL terminator), that many raw bytes, zero-padded to a multiple
    /// of 4. Confirmed against the real `xdr_string`/`xdr_opaque` in the
    /// classic `xdrfile` library's own C source (`xdrfile.c`) — not
    /// `xdrfile_write_string`'s own outer `slen` (`strlen(s) + 1`), which
    /// TRR's header writes as a *separate*, preceding plain integer before
    /// ever calling this (see [`crate::io::trr`]'s own header parsing).
    pub(crate) fn read_string(&mut self) -> Result<String, TrrError> {
        let len = self.read_u32()? as usize;
        let padded = len.div_ceil(4) * 4;
        let bytes = self.take(padded)?;
        String::from_utf8(bytes[..len].to_vec()).map_err(|e| TrrError::ParseError(e.to_string()))
    }

    /// Skips `n` bytes without decoding them — used for legacy header fields
    /// and body blocks (virial, pressure) this crate never stores.
    pub(crate) fn skip(&mut self, n: usize) -> Result<(), TrrError> {
        self.take(n)?;
        Ok(())
    }
}

pub(crate) struct XdrWriter {
    buf: Vec<u8>,
}

impl XdrWriter {
    pub(crate) fn new() -> Self {
        Self { buf: Vec::new() }
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.buf
    }

    pub(crate) fn write_i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_be_bytes());
    }

    pub(crate) fn write_f32(&mut self, v: f32) {
        self.buf.extend_from_slice(&v.to_be_bytes());
    }

    /// See [`XdrReader::read_string`]'s doc comment.
    pub(crate) fn write_string(&mut self, s: &str) {
        let bytes = s.as_bytes();
        self.write_i32(bytes.len() as i32);
        self.buf.extend_from_slice(bytes);
        let padded = bytes.len().div_ceil(4) * 4;
        self.buf.resize(self.buf.len() + (padded - bytes.len()), 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integers_and_floats_round_trip() {
        let mut w = XdrWriter::new();
        w.write_i32(1993);
        w.write_f32(1.5);
        w.write_i32(-7);
        let bytes = w.into_bytes();

        let mut r = XdrReader::new(&bytes);
        assert_eq!(r.read_i32().unwrap(), 1993);
        assert_eq!(r.read_f32().unwrap(), 1.5);
        assert_eq!(r.read_i32().unwrap(), -7);
    }

    #[test]
    fn test_strings_of_every_length_class_round_trip_and_pad_to_a_multiple_of_four() {
        for s in ["", "a", "ab", "abc", "abcd", "GMX_trn_file"] {
            let mut w = XdrWriter::new();
            w.write_string(s);
            let bytes = w.into_bytes();
            assert_eq!(
                bytes.len() % 4,
                0,
                "{s:?} did not pad to a multiple of 4: {} bytes",
                bytes.len()
            );

            let mut r = XdrReader::new(&bytes);
            assert_eq!(r.read_string().unwrap(), s);
        }
    }

    #[test]
    fn test_truncated_input_is_a_clear_error_not_a_panic() {
        let mut r = XdrReader::new(&[0, 0, 0]);
        let err = r.read_i32().unwrap_err();
        assert!(matches!(err, TrrError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_a_truncated_string_length_prefix_is_a_clear_error() {
        // Claims a 100-byte string but the buffer has nothing after the
        // length prefix.
        let mut w = XdrWriter::new();
        w.write_i32(100);
        let bytes = w.into_bytes();
        let mut r = XdrReader::new(&bytes);
        let err = r.read_string().unwrap_err();
        assert!(matches!(err, TrrError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_skip_advances_without_decoding() {
        let mut w = XdrWriter::new();
        w.write_i32(1);
        w.write_i32(2);
        w.write_i32(3);
        let bytes = w.into_bytes();
        let mut r = XdrReader::new(&bytes);
        r.skip(4).unwrap();
        assert_eq!(r.read_i32().unwrap(), 2);
        assert_eq!(r.position(), 8);
    }
}
