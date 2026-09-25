//! XDR (RFC 1014), hand-rolled rather than taken as a dependency.
//!
//! Only the subset TRR and XTC actually need: signed/unsigned 32-bit
//! integers, IEEE single/double floats, and length-prefixed, zero-padded
//! opaque data (both a UTF-8 string and an arbitrary byte blob — XTC's
//! packed-bit coordinate block is the latter). Big-endian throughout, per
//! the spec — independent of [`crate::io::msgpack`]'s own big-endian
//! encoding, which is a different format entirely. Consistent with this
//! crate's existing style of hand-rolling a small, bounded surface rather
//! than depending on a general-purpose crate for it — see
//! [`crate::io::format::Carries`]'s own doc comment.
//!
//! Every primitive here is already a multiple of 4 bytes — an `i32`/`f32` is
//! exactly 4, an `f64` is exactly 8 (already a multiple of 4, so no padding
//! quirk exists for it) — so only this module's own `read_opaque`/
//! `write_opaque` (and `read_string`/`write_string`, built on top of them)
//! need the length-prefix-plus-padding logic RFC 1014 describes.

use crate::io::errors::XdrError;

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

    fn take(&mut self, n: usize) -> Result<&'a [u8], XdrError> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&end| end <= self.bytes.len())
            .ok_or_else(|| XdrError::ParseError("unexpected end of input".to_string()))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    pub(crate) fn read_i32(&mut self) -> Result<i32, XdrError> {
        Ok(i32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub(crate) fn read_u32(&mut self) -> Result<u32, XdrError> {
        Ok(u32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub(crate) fn read_f32(&mut self) -> Result<f32, XdrError> {
        Ok(f32::from_be_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub(crate) fn read_f64(&mut self) -> Result<f64, XdrError> {
        Ok(f64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }

    /// GROMACS's `xdr_int64`: the high then the low 32 bits, which is an
    /// 8-byte big-endian integer (#399).
    pub(crate) fn read_i64(&mut self) -> Result<i64, XdrError> {
        Ok(i64::from_be_bytes(self.take(8)?.try_into().unwrap()))
    }

    pub(crate) fn remaining(&self) -> usize {
        self.bytes.len() - self.pos
    }

    /// RFC 1014's `xdr_opaque`: `len` raw bytes, zero-padded to a multiple
    /// of 4 — **no length of its own** on the wire. XTC's coordinate block
    /// calls exactly this (via `xdrfile_read_opaque`/`xdrfile_write_opaque`)
    /// with a length it already read/wrote as its own separate plain `int`
    /// field; [`Self::read_opaque`]/[`XdrWriter::write_opaque`] below are
    /// the different, length-*prefixed* shape `xdr_string` builds on top of
    /// this, which is what TRR's version string and every other string
    /// field actually need.
    pub(crate) fn read_padded(&mut self, len: usize) -> Result<Vec<u8>, XdrError> {
        let padded = len.div_ceil(4) * 4;
        let bytes = self.take(padded)?;
        Ok(bytes[..len].to_vec())
    }

    /// RFC 1014's `xdr_string`-style framing: a length prefix, then
    /// [`Self::read_padded`]. Confirmed against the real `xdr_string` in
    /// the classic `xdrfile` library's own C source (`xdrfile.c`).
    pub(crate) fn read_opaque(&mut self) -> Result<Vec<u8>, XdrError> {
        let len = self.read_u32()? as usize;
        self.read_padded(len)
    }

    /// [`Self::read_opaque`], then decoded as UTF-8 — RFC 1014's
    /// `xdr_string`. Not `xdrfile_write_string`'s own outer `slen`
    /// (`strlen(s) + 1`), which TRR's header writes as a *separate*,
    /// preceding plain integer before ever calling this (see
    /// [`crate::io::trr`]'s own header parsing).
    pub(crate) fn read_string(&mut self) -> Result<String, XdrError> {
        String::from_utf8(self.read_opaque()?).map_err(|e| XdrError::ParseError(e.to_string()))
    }

    /// Skips `n` bytes without decoding them — used for legacy header fields
    /// and body blocks (virial, pressure) this crate never stores.
    pub(crate) fn skip(&mut self, n: usize) -> Result<(), XdrError> {
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

    pub(crate) fn write_f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_be_bytes());
    }

    /// See [`XdrReader::read_i64`].
    pub(crate) fn write_i64(&mut self, v: i64) {
        self.buf.extend_from_slice(&v.to_be_bytes());
    }

    /// See [`XdrReader::read_padded`]'s doc comment.
    pub(crate) fn write_padded(&mut self, bytes: &[u8]) {
        self.buf.extend_from_slice(bytes);
        let padded = bytes.len().div_ceil(4) * 4;
        self.buf.resize(self.buf.len() + (padded - bytes.len()), 0);
    }

    /// See [`XdrReader::read_opaque`]'s doc comment.
    pub(crate) fn write_opaque(&mut self, bytes: &[u8]) {
        self.write_i32(bytes.len() as i32);
        self.write_padded(bytes);
    }

    /// See [`XdrReader::read_string`]'s doc comment.
    pub(crate) fn write_string(&mut self, s: &str) {
        self.write_opaque(s.as_bytes());
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
    fn test_opaque_bytes_of_every_length_class_round_trip_and_pad_to_a_multiple_of_four() {
        for len in [0, 1, 3, 4, 5, 17] {
            let data = vec![0xab; len];
            let mut w = XdrWriter::new();
            w.write_opaque(&data);
            let bytes = w.into_bytes();
            assert_eq!(
                bytes.len() % 4,
                0,
                "{len} did not pad: {} bytes",
                bytes.len()
            );

            let mut r = XdrReader::new(&bytes);
            assert_eq!(r.read_opaque().unwrap(), data);
        }
    }

    #[test]
    fn test_padded_bytes_carry_no_length_of_their_own_unlike_opaque() {
        for len in [0, 1, 3, 4, 5, 17] {
            let data = vec![0xcd; len];
            let mut w = XdrWriter::new();
            w.write_padded(&data);
            let bytes = w.into_bytes();
            assert_eq!(
                bytes.len(),
                len.div_ceil(4) * 4,
                "{len} should pad with no length prefix"
            );

            let mut r = XdrReader::new(&bytes);
            assert_eq!(r.read_padded(len).unwrap(), data);
        }
    }

    #[test]
    fn test_truncated_input_is_a_clear_error_not_a_panic() {
        let mut r = XdrReader::new(&[0, 0, 0]);
        let err = r.read_i32().unwrap_err();
        assert!(matches!(err, XdrError::ParseError(_)), "{err}");
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
        assert!(matches!(err, XdrError::ParseError(_)), "{err}");
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
