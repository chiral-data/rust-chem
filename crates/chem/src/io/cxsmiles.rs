//! CXSMILES — RDKit's `|...|`-suffixed extension of SMILES (#221).
//!
//! This crate implements the enhanced-stereochemistry-group subset only:
//! `&n:`/`on:`/`a:` groups (see [`crate::core::stereo_group`]). Every other
//! CXSMILES feature — atom coordinates, atom labels/aliases, per-atom
//! properties (`atomProp`), aromatic-bond markers, radical/valence flags —
//! is not read or written; a block containing one is a parse error naming
//! the block, not a silent partial read. This targets RDKit's documented
//! syntax for the stereo-group subset specifically, not every variant
//! spelling that might appear in the wild — the same honesty already
//! applied to canonical SMILES not matching another toolkit's algorithm
//! byte for byte (#220).
//!
//! CXSMARTS (the same extension over SMARTS) is out of scope entirely:
//! this crate has no SMARTS/query-molecule support to extend — #173's own
//! roadmap places that at Wave 5 / v0.12.0, four releases past this
//! milestone.

use nom::{
    IResult, Parser,
    character::complete::{char, digit1, one_of},
    combinator::{all_consuming, map_res, opt},
    multi::separated_list1,
};

use crate::core::molecule::Molecule;
use crate::core::stereo_group::{StereoGroup, StereoGroupKind};
use crate::io::errors::SmilesError;
use crate::io::smiles::parse_smiles;
use crate::io::smiles_writer::write_smiles_for_molecule_canonical;

/// Parses a CXSMILES string: a plain SMILES, optionally followed by a
/// `|...|` enhanced-stereo-group block. `smiles` and `block` (the block's
/// content, with the surrounding `|`s already removed) are already split
/// apart by the caller — the file-line reader ([`crate::io::reader`]) or
/// the streaming supplier ([`crate::io::supplier`]) — so this function does
/// no whitespace splitting of its own, matching [`parse_smiles`]'s own
/// contract of taking exactly one token.
pub fn parse_cxsmiles(smiles: &str, block: Option<&str>) -> Result<Molecule, SmilesError> {
    let mut mol = parse_smiles(smiles)?;
    if let Some(block) = block {
        let groups = parse_stereo_groups(block)?;
        mol.set_stereo_groups(groups)
            .map_err(|e| SmilesError::ParseError(e.to_string()))?;
    }
    Ok(mol)
}

/// Writes `mol` as CXSMILES: #220's canonical SMILES, followed by a
/// `|...|` block naming its stereo groups — omitted entirely when there are
/// none, so a molecule with nothing to say here writes identically to plain
/// canonical SMILES.
pub fn write_cxsmiles(mol: &Molecule) -> String {
    let smiles = write_smiles_for_molecule_canonical(mol);
    if mol.stereo_groups().is_empty() {
        return smiles;
    }
    format!("{smiles} |{}|", write_stereo_groups(mol.stereo_groups()))
}

/// Splits a line into its SMILES token, its `|...|` block content (`|`s
/// stripped), if present, and whatever whitespace-separated tokens remain —
/// the name. The file-level convention CXSMILES adds on top of plain
/// SMILES's "first token is the SMILES, the rest is the name": the block,
/// when present, is always its own whitespace-delimited token immediately
/// after the SMILES (true for every group spec this story parses — none
/// contain a space), so this is plain token classification, not a
/// bracket-matching scan over the whole line.
pub fn split_cxsmiles_line(line: &str) -> (&str, Option<&str>, Vec<&str>) {
    let mut parts = line.split_whitespace();
    let Some(smiles) = parts.next() else {
        return ("", None, Vec::new());
    };
    let rest: Vec<&str> = parts.collect();
    match rest.split_first() {
        Some((first, tail))
            if first.len() >= 2 && first.starts_with('|') && first.ends_with('|') =>
        {
            (smiles, Some(&first[1..first.len() - 1]), tail.to_vec())
        }
        _ => (smiles, None, rest),
    }
}

fn stereo_group_kind(input: &str) -> IResult<&str, char> {
    one_of("&oa")(input)
}

fn atom_index(input: &str) -> IResult<&str, usize> {
    map_res(digit1, str::parse).parse(input)
}

/// One group spec: `&<n>:<atoms>`, `o<n>:<atoms>`, or `a:<atoms>`. `<n>` is
/// optional on `&`/`o` and defaults to 1 (RDKit permits omitting it when
/// there is only one group of that kind); `a` never carries a number.
fn stereo_group_spec(input: &str) -> IResult<&str, StereoGroup> {
    let (input, marker) = stereo_group_kind(input)?;
    let (input, number) = if marker == 'a' {
        (input, None)
    } else {
        let (input, n) = opt(map_res(digit1, str::parse::<u32>)).parse(input)?;
        (input, Some(n.unwrap_or(1)))
    };
    let (input, _) = char(':')(input)?;
    let (input, atoms) = separated_list1(char(','), atom_index).parse(input)?;

    let kind = match marker {
        '&' => StereoGroupKind::And,
        'o' => StereoGroupKind::Or,
        'a' => StereoGroupKind::Absolute,
        _ => unreachable!("one_of(\"&oa\") admits no other character"),
    };
    Ok((input, StereoGroup::new(kind, number, atoms)))
}

/// The comma between two groups and the comma between two atoms in one
/// group's atom list share the same character; nom's list combinators
/// disambiguate them correctly by backtracking — a comma followed by
/// something that isn't a valid next atom index simply isn't consumed as
/// part of the atom list, leaving it for the group-level separator to
/// claim instead.
fn parse_stereo_groups(block: &str) -> Result<Vec<StereoGroup>, SmilesError> {
    let (_, groups) = all_consuming(separated_list1(char(','), stereo_group_spec))
        .parse(block)
        .map_err(|e| SmilesError::ParseError(format!("invalid CXSMILES block {block:?}: {e}")))?;
    Ok(groups)
}

fn write_stereo_groups(groups: &[StereoGroup]) -> String {
    groups
        .iter()
        .map(|group| {
            let marker = match group.kind {
                StereoGroupKind::And => format!("&{}", group.number.unwrap_or(1)),
                StereoGroupKind::Or => format!("o{}", group.number.unwrap_or(1)),
                StereoGroupKind::Absolute => "a".to_string(),
            };
            let atoms = group
                .atoms
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            format!("{marker}:{atoms}")
        })
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_bare_smiles_with_no_block_parses_like_plain_smiles() {
        let mol = parse_cxsmiles("N[C@@H](C)C(=O)O", None).expect("valid SMILES");
        assert_eq!(mol.num_atoms(), 6);
        assert!(mol.stereo_groups().is_empty());
    }

    #[test]
    fn test_an_and_group_lands_on_the_named_atoms() {
        let mol = parse_cxsmiles("N[C@@H](C)C(=O)O", Some("&1:1")).expect("valid CXSMILES");
        assert_eq!(mol.stereo_groups().len(), 1);
        let group = &mol.stereo_groups()[0];
        assert_eq!(group.kind, StereoGroupKind::And);
        assert_eq!(group.number, Some(1));
        assert_eq!(group.atoms, vec![1]);
    }

    #[test]
    fn test_multiple_groups_of_different_kinds_all_parse() {
        // The comma before "o2" is the same character as the comma inside
        // "&1"'s own atom list -- this is the case that actually exercises
        // the disambiguation described on `parse_stereo_groups`.
        let mol = parse_cxsmiles("N[C@@H](C)C(=O)O", Some("&1:1,3,o2:5")).expect("valid CXSMILES");
        assert_eq!(mol.stereo_groups().len(), 2);
        assert_eq!(mol.stereo_groups()[0].atoms, vec![1, 3]);
        assert_eq!(mol.stereo_groups()[1].kind, StereoGroupKind::Or);
        assert_eq!(mol.stereo_groups()[1].atoms, vec![5]);
    }

    #[test]
    fn test_the_absolute_group_carries_no_number() {
        let mol = parse_cxsmiles("N[C@@H](C)C(=O)O", Some("a:1")).expect("valid CXSMILES");
        assert_eq!(mol.stereo_groups()[0].kind, StereoGroupKind::Absolute);
        assert_eq!(mol.stereo_groups()[0].number, None);
    }

    #[test]
    fn test_an_out_of_range_atom_index_is_a_clear_error_not_a_panic() {
        let err = parse_cxsmiles("CC", Some("&1:5")).unwrap_err();
        assert!(err.to_string().contains("atom 5"), "{err}");
    }

    #[test]
    fn test_an_unrecognised_block_is_a_clear_error() {
        let err = parse_cxsmiles("CC", Some("coords:(0,0,0)")).unwrap_err();
        assert!(err.to_string().contains("invalid CXSMILES block"), "{err}");
    }

    #[test]
    fn test_writing_omits_the_block_when_there_are_no_stereo_groups() {
        let mol = parse_smiles("CCO").expect("valid SMILES");
        let written = write_cxsmiles(&mol);
        assert!(!written.contains('|'), "{written}");
        assert_eq!(written, write_smiles_for_molecule_canonical(&mol));
    }

    #[test]
    fn test_stereo_groups_round_trip_through_writing_and_parsing() {
        let mut mol = parse_smiles("N[C@@H](C)C(=O)O").expect("valid SMILES");
        mol.set_stereo_groups(vec![StereoGroup::new(
            StereoGroupKind::And,
            Some(1),
            vec![1],
        )])
        .expect("valid stereo groups");

        let written = write_cxsmiles(&mol);
        let (smiles, block, _name) = split_cxsmiles_line(&written);
        let back = parse_cxsmiles(smiles, block).expect("round trips");
        assert_eq!(back.stereo_groups(), mol.stereo_groups());
    }

    #[test]
    fn test_split_cxsmiles_line_separates_smiles_block_and_name() {
        let (smiles, block, name) = split_cxsmiles_line("CCO |&1:0| ethanol note");
        assert_eq!(smiles, "CCO");
        assert_eq!(block, Some("&1:0"));
        assert_eq!(name, vec!["ethanol", "note"]);
    }

    #[test]
    fn test_split_cxsmiles_line_with_no_block_matches_plain_smiles_convention() {
        let (smiles, block, name) = split_cxsmiles_line("CCO ethanol");
        assert_eq!(smiles, "CCO");
        assert_eq!(block, None);
        assert_eq!(name, vec!["ethanol"]);
    }
}
