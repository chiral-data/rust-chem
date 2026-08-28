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
