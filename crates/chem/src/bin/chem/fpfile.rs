//! The fingerprint file format.
//!
//! Fingerprints are binary, and the CLI has two audiences: a pipeline wants
//! text it can grep and diff, and `chem search` wants to read a file back and
//! fail loudly if the parameters disagree. One self-describing text format
//! serves both, so this story adds no `--binary` flag.
//!
//! ```text
//! # chem-fingerprints 1
//! # radius 2
//! # size   2048
//! name    fingerprint
//! ethanol 0000a0...
//! ```
//!
//! Hex rather than base64 deliberately: each character is exactly four bits, so
//! a bit position can be located by hand and two fingerprints diff at the
//! character that actually differs. Base64 is 25% smaller and neither of those.
//!
//! The `#` metadata is what makes `chem search` able to refuse rather than
//! guess. A query has to be fingerprinted with the same radius and size as the
//! targets or the comparison is meaningless — and meaningless in a way that
//! produces plausible-looking similarities rather than an error.
//!
//! A binary form would be about half the size. That is worth having when
//! someone has a million molecules and not before, and adding it later is
//! additive: this header already carries the version that would distinguish it.

use anyhow::{Context, Result, bail};
use bitvec::prelude::*;

const MAGIC: &str = "# chem-fingerprints 1";

/// Fingerprints plus the parameters they were generated with.
#[derive(Debug)]
pub struct FingerprintFile {
    pub radius: u32,
    pub size: u32,
    pub names: Vec<String>,
    pub fingerprints: Vec<BitVec>,
}

impl FingerprintFile {
    pub fn to_text(&self) -> String {
        let mut out =
            String::with_capacity(self.fingerprints.len() * (self.size as usize / 4 + 32));
        out.push_str(MAGIC);
        out.push('\n');
        out.push_str(&format!("# radius {}\n", self.radius));
        out.push_str(&format!("# size {}\n", self.size));
        out.push_str("name\tfingerprint\n");
        for (name, fp) in self.names.iter().zip(&self.fingerprints) {
            out.push_str(name);
            out.push('\t');
            out.push_str(&to_hex(fp));
            out.push('\n');
        }
        out
    }

    /// The same fingerprints as chemfp's FPS format (#243).
    ///
    /// A second serialiser rather than a replacement: `chem search` parses the
    /// native form above and nothing reads FPS back, here or in OpenBabel, which
    /// marks its own `fps` write-only. This exists so the bits can leave — FPS is
    /// what a similarity-search tool outside this crate expects.
    ///
    /// `#type` names *this crate's* Morgan bits deliberately. The hash differs
    /// from RDKit's on purpose (#192), so a search mixing the two would produce
    /// plausible similarities from unrelated bit positions. Namespacing the type
    /// string is what tells a reader not to.
    pub fn to_fps(&self) -> String {
        let mut out =
            String::with_capacity(self.fingerprints.len() * (self.size as usize / 4 + 32));
        // `#FPS1` is the only line the format requires; the rest are conventional
        // and carry what a reader needs to use these sensibly. FPS has no radius
        // field of its own, so it rides in `#type`, which is where chemfp's own
        // RDKit-Morgan strings put it.
        out.push_str("#FPS1\n");
        out.push_str(&format!("#num_bits={}\n", self.size));
        out.push_str(&format!(
            "#type=chem-Morgan/1 radius={} fpSize={}\n",
            self.radius, self.size
        ));
        out.push_str(&format!("#software=chem/{}\n", env!("CARGO_PKG_VERSION")));
        for (name, fp) in self.names.iter().zip(&self.fingerprints) {
            // Hex first, name second -- the reverse of the native format above.
            out.push_str(&to_fps_hex(fp));
            out.push('\t');
            out.push_str(name);
            out.push('\n');
        }
        out
    }

    pub fn parse(text: &str) -> Result<Self> {
        // Checked before anything else: pointing `chem search` at a SMILES
        // file is an ordinary mistake, and "line 2: expected name<TAB>
        // fingerprint" describes the symptom rather than the mistake.
        let first = text.lines().find(|l| !l.trim().is_empty());
        if first != Some(MAGIC) {
            bail!("not a chem fingerprint file: expected the first line to be `{MAGIC}`");
        }

        let mut radius = None;
        let mut size = None;
        let mut names = Vec::new();
        let mut fingerprints = Vec::new();
        let mut saw_magic = false;

        for (index, raw) in text.lines().enumerate() {
            let line = raw.trim_end();
            if line.is_empty() {
                continue;
            }
            if line == MAGIC {
                saw_magic = true;
                continue;
            }
            if let Some(rest) = line.strip_prefix("# radius ") {
                radius = Some(rest.trim().parse().context("radius in the header")?);
                continue;
            }
            if let Some(rest) = line.strip_prefix("# size ") {
                size = Some(rest.trim().parse().context("size in the header")?);
                continue;
            }
            if line.starts_with('#') || line == "name\tfingerprint" {
                continue;
            }

            let (name, hex) = line
                .split_once('\t')
                .with_context(|| format!("line {}: expected name<TAB>fingerprint", index + 1))?;
            let expected = size.context("the size header must come before the fingerprints")?;
            names.push(name.to_owned());
            fingerprints.push(
                from_hex(hex.trim(), expected).with_context(|| format!("line {}", index + 1))?,
            );
        }

        debug_assert!(saw_magic, "checked before the loop");

        Ok(Self {
            radius: radius.context("missing `# radius` header")?,
            size: size.context("missing `# size` header")?,
            names,
            fingerprints,
        })
    }
}

/// Least-significant-bit-first within each nibble, so bit *n* of the
/// fingerprint is at a position derivable from *n* alone.
/// Hex in FPS's bit order, which is not this file's own.
///
/// FPS packs bit *i* as bit *i % 8* of byte *i / 8* -- little-endian within
/// each byte -- while [`to_hex`] above is little-endian within each *nibble*.
/// The two differ by swapping the hex characters of every byte, and the
/// difference is invisible: both produce a plausible fingerprint that a reader
/// will happily compare against the wrong bits.
///
/// Verified against OpenBabel's own FPS output rather than derived. Decoding
/// one of its lines four ways against the bits it reports set:
///
/// ```text
/// msb          [516, 532, 665]   no
/// nibble_lsb   [519, 535, 666]   no      <- to_hex's convention
/// byte_le      [515, 531, 670]   yes, == OpenBabel's 1-indexed bits
/// byte_be      [516, 532, 665]   no
/// ```
fn to_fps_hex(fp: &BitVec) -> String {
    let mut out = String::with_capacity(fp.len() / 4 + 2);
    for byte in fp.chunks(8) {
        let mut value = 0u8;
        for (i, bit) in byte.iter().enumerate() {
            if *bit {
                value |= 1 << i;
            }
        }
        out.push_str(&format!("{value:02x}"));
    }
    out
}

fn to_hex(fp: &BitVec) -> String {
    let mut out = String::with_capacity(fp.len() / 4 + 1);
    for chunk in fp.chunks(4) {
        let mut nibble = 0u8;
        for (i, bit) in chunk.iter().enumerate() {
            if *bit {
                nibble |= 1 << i;
            }
        }
        out.push(char::from_digit(nibble as u32, 16).expect("nibble is 0..16"));
    }
    out
}

fn from_hex(hex: &str, size: u32) -> Result<BitVec> {
    let expected_chars = size.div_ceil(4) as usize;
    if hex.len() != expected_chars {
        bail!(
            "expected {expected_chars} hex characters for a {size}-bit fingerprint, got {}",
            hex.len()
        );
    }

    let mut fp = BitVec::repeat(false, size as usize);
    for (index, c) in hex.chars().enumerate() {
        let nibble = c
            .to_digit(16)
            .with_context(|| format!("{c:?} is not a hex digit"))? as u8;
        for i in 0..4 {
            let bit = index * 4 + i;
            if bit < fp.len() && nibble & (1 << i) != 0 {
                fp.set(bit, true);
            }
        }
    }
    Ok(fp)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> FingerprintFile {
        let mut a = BitVec::repeat(false, 16);
        a.set(0, true);
        a.set(5, true);
        a.set(15, true);
        FingerprintFile {
            radius: 2,
            size: 16,
            names: vec!["ethanol".into(), "empty".into()],
            fingerprints: vec![a, BitVec::repeat(false, 16)],
        }
    }

    /// Bit indices decoded from FPS hex: byte *k*, bit *i % 8* little-endian
    /// within it -- the convention `to_fps_hex` writes and `to_hex` does not.
    fn decode_fps_hex(hex: &str) -> Vec<usize> {
        let mut bits = Vec::new();
        for (k, pair) in hex.as_bytes().chunks(2).enumerate() {
            let byte = u8::from_str_radix(std::str::from_utf8(pair).expect("ascii"), 16)
                .expect("hex byte");
            for offset in 0..8 {
                if byte & (1 << offset) != 0 {
                    bits.push(k * 8 + offset);
                }
            }
        }
        bits
    }

    #[test]
    fn test_fps_header_names_this_crate_not_rdkit() {
        let text = sample().to_fps();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "#FPS1", "the one line the format requires");
        assert!(lines.contains(&"#num_bits=16"), "{text}");
        // The hash differs from RDKit's deliberately (#192), so a reader has
        // to be able to tell these are not RDKit-Morgan bits.
        assert!(
            lines
                .iter()
                .any(|l| l.starts_with("#type=chem-Morgan/1 radius=2")),
            "{text}"
        );
    }

    #[test]
    fn test_fps_uses_its_own_bit_order_not_ours() {
        // The failure this format can ship while looking entirely plausible:
        // reusing `to_hex` produces a well-formed file whose bits are in the
        // wrong places, and every reader would compare against them happily.
        let file = sample();
        let fps_hex = file.to_fps();
        let row = fps_hex
            .lines()
            .find(|l| l.starts_with(|c: char| c.is_ascii_hexdigit()));
        let hex = row
            .expect("a data row")
            .split('\t')
            .next()
            .expect("hex field");

        assert_eq!(
            decode_fps_hex(hex),
            vec![0, 5, 15],
            "decoded FPS bits must be the ones the fingerprint actually holds"
        );
        // And the two conventions genuinely differ here, or the test above
        // would pass for the wrong reason.
        assert_ne!(hex, to_hex(&file.fingerprints[0]), "to_hex was reused");
    }

    #[test]
    fn test_fps_puts_the_hex_first_and_the_name_second() {
        // Reversed from the native format, and easy to lose.
        let text = sample().to_fps();
        let row = text.lines().find(|l| l.contains("ethanol")).expect("a row");
        let (hex, name) = row.split_once('\t').expect("two fields");
        assert_eq!(name, "ethanol");
        assert!(hex.chars().all(|c| c.is_ascii_hexdigit()), "{row}");
    }

    #[test]
    fn test_fps_with_no_fingerprints_is_still_a_valid_file() {
        let empty = FingerprintFile {
            radius: 2,
            size: 16,
            names: vec![],
            fingerprints: vec![],
        };
        let text = empty.to_fps();
        assert!(text.starts_with("#FPS1\n"), "{text}");
        assert!(
            text.lines().all(|l| l.starts_with('#')),
            "a header and no rows, not a malformed file: {text}"
        );
    }

    #[test]
    fn test_a_file_survives_a_round_trip() {
        let original = sample();
        let parsed = FingerprintFile::parse(&original.to_text()).expect("round trip");
        assert_eq!(parsed.radius, 2);
        assert_eq!(parsed.size, 16);
        assert_eq!(parsed.names, original.names);
        assert_eq!(parsed.fingerprints, original.fingerprints);
    }

    #[test]
    fn test_the_parameters_are_in_the_file_not_assumed() {
        // The whole reason the header exists: `chem search` must be able to
        // refuse a query fingerprinted with different parameters, and it can
        // only do that if the file says what its own were.
        let text = sample().to_text();
        assert!(text.contains("# radius 2"));
        assert!(text.contains("# size 16"));
    }

    #[test]
    fn test_a_file_without_the_magic_line_is_refused() {
        // Otherwise pointing `chem search` at the wrong file gives a confusing
        // per-line parse error instead of "this is not a fingerprint file".
        let text = sample().to_text().replace(MAGIC, "# something else");
        let err = FingerprintFile::parse(&text).expect_err("should refuse");
        assert!(err.to_string().contains("not a chem fingerprint file"));
    }

    #[test]
    fn test_a_smiles_file_is_named_as_the_wrong_kind_of_file() {
        // The mistake this catches in practice: `chem search mols.smi` instead
        // of `chem search fps.tsv`. Identifying the file has to happen before
        // parsing its rows, or the error describes the symptom — a row that is
        // not a fingerprint — rather than the mistake.
        let err =
            FingerprintFile::parse("CCO ethanol\nc1ccccc1 benzene\n").expect_err("should refuse");
        assert!(
            err.to_string().contains("not a chem fingerprint file"),
            "got {err}"
        );
    }

    #[test]
    fn test_a_truncated_fingerprint_is_refused() {
        let text = sample().to_text().replace("\tf", "\t");
        let result = FingerprintFile::parse(&text);
        assert!(
            result.is_err(),
            "a short hex run must not silently zero-pad"
        );
    }

    #[test]
    fn test_hex_is_greppable_for_a_known_bit() {
        // Bit 0 set and nothing else in the first nibble means the first
        // character is `1`. This is the property base64 would cost.
        let mut fp = BitVec::repeat(false, 8);
        fp.set(0, true);
        assert_eq!(to_hex(&fp), "10");
        fp.set(3, true);
        assert_eq!(to_hex(&fp), "90");
    }

    #[test]
    fn test_a_size_that_is_not_a_multiple_of_four_still_round_trips() {
        let mut fp = BitVec::repeat(false, 10);
        fp.set(9, true);
        let hex = to_hex(&fp);
        assert_eq!(from_hex(&hex, 10).expect("round trip"), fp);
    }
}
