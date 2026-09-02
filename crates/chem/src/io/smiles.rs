use std::collections::HashMap;

use nom::{
    IResult, Parser,
    branch::alt,
    bytes::complete::tag,
    character::complete::{char, digit1, one_of},
    combinator::{map, map_res, opt},
    multi::many0,
    sequence::{delimited, preceded},
};

use crate::core::prelude::*;
use crate::io::errors::SmilesError;

#[derive(Debug, Clone)]
enum Token {
    Atom(AtomToken),
    Bond(BondToken),
    Branch,
    BranchEnd,
    Ring(u32),
    /// `.` — the next atom starts a new component rather than bonding to the
    /// previous one. Salts and co-crystals are ordinary catalogue entries, so
    /// this is not an exotic corner (#191).
    Dot,
}

#[derive(Debug, Clone)]
struct AtomToken {
    element: String,
    aromatic: bool,
    charge: i8,
    isotope: Option<u16>,
    h_count: Option<u8>,
    chirality: Chirality,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BondToken {
    Single,
    Double,
    Triple,
    Quadruble,
    Aromatic,
    /// `/` and `\` — single bonds that also say which way they lean.
    ///
    /// The direction is not a property of this bond at all: it describes the
    /// configuration of the *double* bond next to it. Carried through
    /// tokenisation and resolved once the whole molecule is built, in
    /// [`resolve_directional_bonds`].
    Up,
    Down,
}

impl BondToken {
    fn to_bond_order(self) -> BondOrder {
        match self {
            BondToken::Single => BondOrder::Single,
            BondToken::Double => BondOrder::Double,
            BondToken::Triple => BondOrder::Triple,
            BondToken::Quadruble => BondOrder::Quadruple,
            BondToken::Aromatic => BondOrder::Aromatic,
            // A directional bond is a single bond. The slash carries stereo,
            // not order.
            BondToken::Up | BondToken::Down => BondOrder::Single,
        }
    }

    /// Which way the marker leans, as written: `/` up, `\` down.
    fn direction(self) -> Option<i8> {
        match self {
            BondToken::Up => Some(1),
            BondToken::Down => Some(-1),
            _ => None,
        }
    }
}

pub fn parse_smiles(smiles: &str) -> Result<Molecule, SmilesError> {
    let tokens = tokenize_smiles(smiles)?;
    build_molecule(&tokens)
}

///////////////////
// TOKEN PARSING //
///////////////////

fn tokenize_smiles(input: &str) -> Result<Vec<Token>, SmilesError> {
    match parse_smiles_tokens(input) {
        Ok((remaining, tokens)) => {
            if !remaining.is_empty() {
                Err(SmilesError::ParseError(format!(
                    "Parsing {input} with error Unexpected characters: {remaining}",
                )))
            } else {
                Ok(tokens)
            }
        }
        Err(e) => Err(SmilesError::ParseError(e.to_string())),
    }
}

fn parse_smiles_tokens(input: &str) -> IResult<&str, Vec<Token>> {
    many0(alt((
        map(parse_bracket_atom, Token::Atom),
        map(parse_element_atom, Token::Atom),
        map(parse_bond, Token::Bond),
        map(char('('), |_| Token::Branch),
        map(char(')'), |_| Token::BranchEnd),
        map(parse_ring_number, Token::Ring),
        map(char('.'), |_| Token::Dot),
    )))
    .parse(input)
}

fn parse_bracket_atom(input: &str) -> IResult<&str, AtomToken> {
    delimited(
        char('['),
        |input| {
            let (input, isotope) =
                opt(map(digit1, |s: &str| s.parse::<u16>().unwrap())).parse(input)?;
            let (input, symbol) = parse_element_symbol(input)?;
            let (input, aromatic) = if symbol.chars().next().unwrap().is_lowercase() {
                (input, true)
            } else {
                (input, false)
            };

            let element = if aromatic {
                symbol.to_uppercase()
            } else {
                symbol.to_string()
            };

            // Order inside brackets is fixed by the spec: isotope, symbol,
            // chirality, hydrogen count, charge. Parsing chirality after the
            // H count would reject `[C@H]`, which is the common spelling.
            let (input, chirality) = opt(parse_chirality).parse(input)?;
            let (input, h_count) = opt(parse_h_count).parse(input)?;
            let (input, charge) = opt(parse_charge).parse(input)?;

            Ok((
                input,
                AtomToken {
                    element,
                    aromatic,
                    charge: charge.unwrap_or(0),
                    isotope,
                    h_count,
                    chirality: chirality.unwrap_or(Chirality::None),
                },
            ))
        },
        char(']'),
    )
    .parse(input)
}

fn parse_element_atom(input: &str) -> IResult<&str, AtomToken> {
    let (remaining_input, symbol) = parse_element_symbol(input)?;

    // Determine aromaticity: lowercase first character indicates aromatic atom
    let aromatic = symbol.chars().next().unwrap().is_lowercase();

    // Normalize element symbol to uppercase for consistency
    let element = if aromatic {
        symbol.to_uppercase()
    } else {
        symbol
    };

    Ok((
        remaining_input,
        AtomToken {
            element,
            aromatic,
            charge: 0,
            isotope: None,
            h_count: None,
            // The organic-subset shorthand has nowhere to write a chirality
            // marker; `[C@H]` is the only spelling.
            chirality: Chirality::None,
        },
    ))
}

/// Parses element symbols from SMILES strings.
///
/// This parser handles:
/// - All standard element symbols from the periodic table (case-sensitive)
/// - SMILES aromatic lowercase symbols (b, c, n, o, p, s, se, as)
///
/// The parser tries to match the longest possible symbol first (2 chars),
/// then falls back to single character symbols.
fn parse_element_symbol(input: &str) -> IResult<&str, String> {
    alt((
        // Try 2-character symbols first (to avoid matching "C" when input is "Cl")
        parse_two_char_element,
        // Then try 1-character symbols
        parse_one_char_element,
    ))
    .parse(input)
}

/// Parses 2-character element symbols.
///
/// This dynamically generates parsers for all 2-character elements from ELEMENT_SYMBOLS,
/// plus SMILES aromatic symbols (se, as).
fn parse_two_char_element(input: &str) -> IResult<&str, String> {
    // Need at least 2 characters
    if input.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let two_chars = &input[0..2];

    // Check if it matches any 2-character element from ELEMENT_SYMBOLS
    let is_valid_element = ELEMENT_SYMBOLS
        .iter()
        .any(|&symbol| symbol.len() == 2 && symbol == two_chars);

    // Also check for SMILES aromatic 2-character symbols
    let is_aromatic = matches!(two_chars, "se" | "as");

    if is_valid_element || is_aromatic {
        Ok((&input[2..], two_chars.to_string()))
    } else {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )))
    }
}

/// Parses 1-character element symbols.
///
/// This handles:
/// - Single-letter elements from ELEMENT_SYMBOLS (H, B, C, N, O, F, P, S, I, etc.)
/// - SMILES aromatic lowercase symbols (b, c, n, o, p, s)
fn parse_one_char_element(input: &str) -> IResult<&str, String> {
    if input.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let one_char = &input[0..1];

    // Check if it matches any 1-character element from ELEMENT_SYMBOLS
    let is_valid_element = ELEMENT_SYMBOLS
        .iter()
        .any(|&symbol| symbol.len() == 1 && symbol == one_char);

    // Also check for SMILES aromatic 1-character symbols (lowercase)
    let is_aromatic = matches!(one_char, "b" | "c" | "n" | "o" | "p" | "s");

    if is_valid_element || is_aromatic {
        Ok((&input[1..], one_char.to_string()))
    } else {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )))
    }
}

/// Tetrahedral chirality: `@` or `@@`.
///
/// Bracket atoms only — the organic-subset shorthand has nowhere to put it.
/// `@@` must be tried first or `alt` matches the single `@` and leaves the
/// second one to fail the closing bracket.
///
/// Only the tetrahedral case. `@AL1`, `@TB1` and `@OH1` describe allene,
/// trigonal-bipyramidal and octahedral centres, which the data model has no
/// room for and no format in this milestone needs.
fn parse_chirality(input: &str) -> IResult<&str, Chirality> {
    alt((
        map(tag("@@"), |_| Chirality::Clockwise),
        map(char('@'), |_| Chirality::CounterClockwise),
    ))
    .parse(input)
}

fn parse_h_count(input: &str) -> IResult<&str, u8> {
    preceded(
        char('H'),
        map(opt(digit1), |s: Option<&str>| {
            s.map(|d| d.parse::<u8>().unwrap()).unwrap_or(1)
        }),
    )
    .parse(input)
}

fn parse_charge(input: &str) -> IResult<&str, i8> {
    alt((
        map(tag("++"), |_| 2),
        map(tag("--"), |_| -2),
        map(preceded(char('+'), opt(digit1)), |d| {
            d.map(|s: &str| s.parse::<i8>().unwrap()).unwrap_or(1)
        }),
        map(preceded(char('-'), opt(digit1)), |d| {
            -(d.map(|s: &str| s.parse::<i8>().unwrap()).unwrap_or(1))
        }),
    ))
    .parse(input)
}

fn parse_bond(input: &str) -> IResult<&str, BondToken> {
    alt((
        map(char('-'), |_| BondToken::Single),
        map(char('='), |_| BondToken::Double),
        map(char('#'), |_| BondToken::Triple),
        map(char(':'), |_| BondToken::Aromatic),
        map(char('$'), |_| BondToken::Quadruble),
        map(char('/'), |_| BondToken::Up),
        map(char('\\'), |_| BondToken::Down),
    ))
    .parse(input)
}

/// Parses a ring-closure label.
///
/// Two forms, and the distinction between them is load-bearing: **a bare label
/// is exactly one digit.** `%NN` exists precisely so a label of 10 or more can
/// be written next to another digit without ambiguity.
///
/// Consuming more than one bare digit is what #69 was: `c12` means *close ring
/// 1, then close ring 2*, and reading it as a single ring 12 leaves both open.
/// An atom closing two labels is not exotic — it is the shared atom of any
/// fused ring system, so `c1ccnc2ccccc12` (quinoline) and
/// `c1ccc2c(c1)ccc1ccccc12` (anthracene) both failed while the same ring
/// systems spelled to avoid adjacency parsed fine.
fn parse_ring_number(input: &str) -> IResult<&str, u32> {
    alt((
        // `map_res` rather than `map` + `unwrap`: a long enough digit run
        // overflows u32, and a library that panics on a malformed line turns
        // one bad record into a dead process. `%` still takes any number of
        // digits, which is looser than the spec's exactly-two — deliberately
        // left alone, since no real molecule has 100 open rings and tightening
        // it is not what #69 is about.
        map_res(preceded(char('%'), digit1), |s: &str| s.parse()),
        // Exactly one digit. `one_of` cannot match a non-digit, so the
        // conversion is infallible.
        map(one_of("0123456789"), |c: char| {
            c.to_digit(10).expect("one_of matched a digit")
        }),
    ))
    .parse(input)
}

//////////////////////////
// MOLECULE BUILD LOGIC //
//////////////////////////

fn build_molecule(tokens: &[Token]) -> Result<Molecule, SmilesError> {
    let mut mol = Molecule::new();
    let mut stack: Vec<usize> = Vec::new();
    let mut current_atom: Option<usize> = None;
    let mut next_bond: Option<BondToken> = None;
    let mut dangling_bond = false;
    let mut rings: HashMap<u32, (usize, Option<BondToken>)> = HashMap::new();
    // Every `/` or `\` as written: (first atom, second atom, lean). Written
    // order matters and cannot be recovered from the finished graph, which is
    // why this is collected here rather than derived later.
    let mut directional: Vec<(usize, usize, i8)> = Vec::new();
    // A dot with nothing before it, or nothing after it by the time we finish.
    let mut dangling_dot = false;
    let mut after_dot = false;

    for token in tokens {
        match token {
            Token::Atom(atom_token) => {
                let atom_idx = add_atom_from_token(&mut mol, atom_token)?;

                if let Some(prev_atom) = current_atom {
                    // With no explicit bond symbol, a bond is aromatic only when
                    // **both** atoms are. Consulting only the atom being added
                    // makes the result depend on which end the author wrote
                    // first: `c1ccccc1O` gave the C-O bond Single and
                    // `Oc1ccccc1` gave it Aromatic, for the same molecule. That
                    // is #167, and it put a wrong fingerprint on 32 of the 122
                    // molecules in test.smi — every substituted aromatic written
                    // substituent-first.
                    //
                    // The ring-closure branch below has always had this right;
                    // its comment claims this path "already infers this", which
                    // is what kept the bug hidden.
                    let implicit_order = if atom_token.aromatic && mol.atom(prev_atom).is_aromatic()
                    {
                        BondOrder::Aromatic
                    } else {
                        BondOrder::Single
                    };
                    let bond_order = next_bond // the pending_bond, applied NOW
                        .map(|b| b.to_bond_order())
                        .unwrap_or(implicit_order);

                    let mut bond = Bond::new(prev_atom, atom_idx, bond_order);
                    if bond_order == BondOrder::Aromatic {
                        bond = bond.with_aromatic(true);
                    }
                    mol.add_bond(bond)
                        .map_err(|e| SmilesError::ParseError(e.to_string()))?;

                    // Kept in written order — `prev_atom` came first in the
                    // string. `resolve_directional_bonds` needs that, because
                    // the same two symbols mean opposite things depending on
                    // which side of the double bond they were written.
                    if let Some(lean) = next_bond.and_then(BondToken::direction) {
                        directional.push((prev_atom, atom_idx, lean));
                    }
                }

                current_atom = Some(atom_idx);
                after_dot = false;
                next_bond = None; // pending_bond: i.e. it will be applied once next atom arrives
            }

            Token::Bond(bond_token) => {
                if current_atom.is_none() {
                    // `=CC` — a bond before any atom. Recorded rather than
                    // returned here, so that a string of nothing but bond
                    // characters still reports NoAtoms below. `$` is the
                    // quadruple-bond token, so `$$$$` reaches this arm four
                    // times, and "No atoms in SMILES" names that mistake far
                    // better than "a bond has no atom to attach to" would.
                    dangling_bond = true;
                }
                next_bond = Some(*bond_token);
            }

            Token::Branch => {
                if let Some(atom) = current_atom {
                    stack.push(atom);
                }
            }

            Token::BranchEnd => {
                current_atom = Some(stack.pop().ok_or(SmilesError::MismatchedBranches)?);
                next_bond = None;
            }

            Token::Dot => {
                // Breaks the chain: the next atom starts a component rather
                // than bonding to this one.
                //
                // Ring labels deliberately survive. `C1.C1` is legal and bonds
                // the two components together — the dot ends the *chain*, not
                // the molecule, so `rings` is left alone.
                if current_atom.is_none() {
                    dangling_dot = true;
                }
                if next_bond.is_some() {
                    // `C=.C` — a bond symbol with a dot where its second atom
                    // should be.
                    return Err(SmilesError::DanglingBond);
                }
                current_atom = None;
                after_dot = true;
            }

            Token::Ring(ring_num) => {
                let current = current_atom.ok_or(SmilesError::InvalidRing(*ring_num))?;

                if let Some((ring_start, ring_bond)) = rings.remove(ring_num) {
                    // Same rule as the atom-to-atom path above: with no
                    // explicit bond symbol, a bond is aromatic only when both
                    // atoms are. Ring-closure bonds used to default to Single
                    // regardless, so `c1ccccc1` came out with five aromatic
                    // bonds and one single one (#94).
                    let implicit_order =
                        if mol.atom(ring_start).is_aromatic() && mol.atom(current).is_aromatic() {
                            BondOrder::Aromatic
                        } else {
                            BondOrder::Single
                        };
                    let bond_order = next_bond
                        .or(ring_bond)
                        .map(|b| b.to_bond_order())
                        .unwrap_or(implicit_order);

                    let mut bond = Bond::new(ring_start, current, bond_order);
                    if bond_order == BondOrder::Aromatic {
                        bond = bond.with_aromatic(true);
                    }
                    mol.add_bond(bond)
                        .map_err(|e| SmilesError::ParseError(e.to_string()))?;
                    next_bond = None;
                } else {
                    rings.insert(*ring_num, (current, next_bond));
                    next_bond = None;
                }
            }
        }
    }

    if !rings.is_empty() {
        // Sorted: `rings` is a HashMap, so taking an arbitrary key made the same
        // input report a different ring between runs (#153).
        let mut open: Vec<u32> = rings.keys().copied().collect();
        open.sort_unstable();
        return Err(SmilesError::UnclosedRings(open));
    }

    if !stack.is_empty() {
        return Err(SmilesError::MismatchedBranches);
    }

    // Last, so the more specific failures above keep their better messages.
    //
    // Bond characters (`-`, `=`, `#`, `:`, `$`) and `(` are all legal tokens
    // in isolation, so a string of nothing but those tokenized cleanly and
    // built an empty molecule, which was then returned as a success. `$` is
    // the one that bites: it is the quadruple-bond token *and* the character
    // SDF uses for its `$$$$` record terminator, so an SDF file read as SMILES
    // reported one atomless molecule per record rather than failing. A caller
    // checking "did anything parse?" was told yes.
    if mol.num_atoms() == 0 {
        return Err(SmilesError::NoAtoms);
    }

    // After the atom-count check, so an atomless string keeps the better
    // message. `next_bond` still set means a trailing bond: `CC=` (#155).
    if dangling_bond || next_bond.is_some() {
        return Err(SmilesError::DanglingBond);
    }

    // `after_dot` still set means the string ended on one: `CC.`
    if dangling_dot || after_dot {
        return Err(SmilesError::DanglingDot);
    }

    resolve_directional_bonds(&mut mol, &directional);

    mol.calculate_implicit_hydrogens();
    Ok(mol)
}

fn add_atom_from_token(mol: &mut Molecule, token: &AtomToken) -> Result<usize, SmilesError> {
    let atomic_number = ELEMENT_SYMBOLS
        .iter()
        .position(|&symbol| symbol == token.element.as_str())
        .ok_or_else(|| SmilesError::InvalidElement(token.element.clone()))?;

    let element = Element::new(atomic_number as u8)
        .ok_or_else(|| SmilesError::InvalidElement(token.element.clone()))?;

    let mut atom = Atom::new(element);

    if token.aromatic {
        atom = atom.with_aromatic(true);
    }

    if token.charge != 0 {
        atom = atom.with_charge(token.charge);
    }

    if let Some(isotope) = token.isotope {
        atom = atom.with_isotope(isotope);
    }

    if token.chirality != Chirality::None {
        atom = atom.with_chirality(token.chirality);
    }

    let idx = mol.add_atom(atom);

    if let Some(h_count) = token.h_count {
        mol.atom_mut(idx).set_explicit_hydrogens(h_count);
    }

    Ok(idx)
}

/// Turns `/` and `\` markers into [`BondStereo`] on the double bonds they
/// describe.
///
/// A directional marker says nothing about its own bond: it says which side of
/// an adjacent double bond a substituent sits on. So the markers are collected
/// during the token walk and resolved here, once every atom index is final.
///
/// # Written order is the whole difficulty
///
/// `F/C=C/F` is trans and `F/C=C\F` is cis — same symbols in the first
/// position, opposite meaning. But `C(/F)=C/F` is *cis* despite also carrying
/// two `/`, because the first marker is written pointing away from the double
/// bond rather than towards it.
///
/// The question a marker answers is *which side of the double-bond axis the
/// substituent sits on*, and that flips with the writing direction. In `F/C`
/// the slash rises from F into C, so F sits **below** the axis. In `C(/F)` it
/// rises from C out to F, so F sits **above** it — the same character, the
/// opposite side.
///
/// So: the substituent's side is the lean when the double-bond atom is written
/// first, and the negated lean when it is written second. Two substituents on
/// the same side is Z; opposite sides is E.
///
/// Anything the markers do not describe is left alone. A double bond with a
/// marker on only one end is under-specified, not wrong, and stays
/// [`BondStereo::None`] rather than being guessed at.
fn resolve_directional_bonds(mol: &mut Molecule, directional: &[(usize, usize, i8)]) {
    if directional.is_empty() {
        return;
    }

    /// Which side of the axis the substituent sits on, seen from `anchor` —
    /// one end of the double bond.
    fn side_at(anchor: usize, &(first, second, lean): &(usize, usize, i8)) -> Option<i8> {
        if first == anchor {
            // `C(/F)=...` — the slash rises away from the anchor, so the
            // substituent is above.
            Some(lean)
        } else if second == anchor {
            // `F/C=...` — the slash rises into the anchor, so the substituent
            // it came from is below.
            Some(-lean)
        } else {
            None
        }
    }

    let doubles: Vec<(usize, usize, usize)> = mol
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, bond)| bond.order() == BondOrder::Double)
        .map(|(idx, bond)| (idx, bond.atom1(), bond.atom2()))
        .collect();

    for (bond_idx, left, right) in doubles {
        let left_side = directional.iter().find_map(|marker| side_at(left, marker));
        let right_side = directional.iter().find_map(|marker| side_at(right, marker));

        let (Some(left_side), Some(right_side)) = (left_side, right_side) else {
            continue;
        };

        let stereo = if left_side == right_side {
            BondStereo::Z
        } else {
            BondStereo::E
        };
        let existing = mol.bond(bond_idx);
        let replacement = Bond::new(existing.atom1(), existing.atom2(), existing.order())
            .with_aromatic(existing.is_aromatic())
            .with_stereo(stereo);
        *mol.bond_mut(bond_idx) = replacement;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stereo_of_double(smiles: &str) -> BondStereo {
        let mol = parse_smiles(smiles).expect(smiles);
        mol.bonds()
            .iter()
            .find(|b| b.order() == BondOrder::Double)
            .expect("a double bond")
            .stereo()
    }

    #[test]
    fn test_dot_starts_a_new_component() {
        // Every salt and co-crystal in a catalogue was a skipped record before
        // this (#191).
        let mol = parse_smiles("[Na+].[Cl-]").expect("valid SMILES");
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 0);
        assert_eq!(mol.graph().connected_components().len(), 2);
        assert_eq!(mol.atom(0).formal_charge(), 1);
        assert_eq!(mol.atom(1).formal_charge(), -1);
    }

    #[test]
    fn test_a_ring_label_survives_a_dot() {
        // The dot ends the chain, not the molecule. `C1.C1` closes ring 1
        // across the break, so the two "components" are one bonded pair — a
        // parser that cleared its ring table on `.` would report an unclosed
        // ring instead.
        let mol = parse_smiles("C1.C1").expect("valid SMILES");
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.graph().connected_components().len(), 1);
    }

    #[test]
    fn test_a_dot_needs_something_on_both_sides() {
        for bad in [".CC", "CC."] {
            assert!(
                matches!(parse_smiles(bad), Err(SmilesError::DanglingDot)),
                "{bad} should be a dangling dot, got {:?}",
                parse_smiles(bad)
            );
        }
        // Not a DanglingDot: the bond symbol is the thing left unattached.
        assert!(matches!(
            parse_smiles("C=.C"),
            Err(SmilesError::DanglingBond)
        ));
    }

    #[test]
    fn test_chirality_is_read_and_the_two_markers_differ() {
        let anticlockwise = parse_smiles("N[C@H](C)C(=O)O").expect("valid SMILES");
        let clockwise = parse_smiles("N[C@@H](C)C(=O)O").expect("valid SMILES");

        assert_eq!(
            anticlockwise.atom(1).chirality(),
            Chirality::CounterClockwise
        );
        assert_eq!(clockwise.atom(1).chirality(), Chirality::Clockwise);
        // The whole point: the two spellings are different molecules. Dropping
        // the marker yields a valid molecule that is the wrong one. Compared
        // by the field rather than by the molecule — `Molecule` has no
        // `PartialEq`, and giving it one means deciding whether two graphs
        // with different atom orders are equal, which is a canonicalisation
        // question and not this story's.
        assert_ne!(
            anticlockwise.atom(1).chirality(),
            clockwise.atom(1).chirality()
        );
    }

    #[test]
    fn test_an_unmarked_atom_has_no_chirality() {
        let mol = parse_smiles("N[CH](C)C(=O)O").expect("valid SMILES");
        assert_eq!(mol.atom(1).chirality(), Chirality::None);
    }

    #[test]
    fn test_directional_bonds_resolve_to_e_and_z() {
        // The textbook pair.
        assert_eq!(stereo_of_double("F/C=C/F"), BondStereo::E, "trans");
        assert_eq!(stereo_of_double("F/C=C\\F"), BondStereo::Z, "cis");
    }

    #[test]
    fn test_the_same_markers_mean_the_opposite_when_written_reversed() {
        // The case that makes this more than a lookup, and the one a naive
        // implementation gets wrong: `C(/F)=C/F` carries two `/` just like
        // `F/C=C/F`, but the first marker is written pointing away from the
        // double bond rather than into it, so the substituents end up on the
        // same side. Same characters, opposite answer.
        assert_eq!(stereo_of_double("C(/F)=C/F"), BondStereo::Z);
        assert_eq!(stereo_of_double("C(\\F)=C/F"), BondStereo::E);
    }

    #[test]
    fn test_a_half_marked_double_bond_is_left_alone() {
        // Under-specified rather than wrong. Guessing a configuration from one
        // marker would invent stereochemistry the author did not write.
        assert_eq!(stereo_of_double("F/C=CF"), BondStereo::None);
        assert_eq!(stereo_of_double("FC=CF"), BondStereo::None);
    }

    #[test]
    fn test_methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol.num_atoms(), 1);
        assert_eq!(mol.formula(), "CH4");
    }

    #[test]
    fn test_ethane() {
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.formula(), "C2H6");
    }

    #[test]
    fn test_ethene() {
        let mol = parse_smiles("C=C").unwrap();
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.formula(), "C2H4");
    }

    #[test]
    fn test_water() {
        let mol = parse_smiles("O").unwrap();
        assert_eq!(mol.formula(), "H2O");
    }

    #[test]
    fn test_propane() {
        let mol = parse_smiles("CCC").unwrap();
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.formula(), "C3H8");
    }

    #[test]
    fn test_branch() {
        let mol = parse_smiles("CC(C)C").unwrap();
        assert_eq!(mol.num_atoms(), 4);
        assert_eq!(mol.formula(), "C4H10");
    }

    #[test]
    fn test_complex_mol() {
        let mol = parse_smiles("CC(=O)c1ccccc1").unwrap();
        assert_eq!(mol.num_atoms(), 9); // 8 C + 1 O == 9
        assert_eq!(mol.formula(), "C8H8O");
    }

    #[test]
    fn test_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.num_atoms(), 6);
        assert_eq!(mol.formula(), "C6H6");
        assert!(mol.atom(0).is_aromatic());
    }

    #[test]
    fn test_charged() {
        let mol = parse_smiles("[NH4+]").unwrap();
        assert_eq!(mol.atom(0).formal_charge(), 1);
        assert_eq!(mol.formula(), "H4N");
    }

    #[test]
    fn test_isotope() {
        let mol = parse_smiles("[13C]").unwrap();
        assert_eq!(mol.atom(0).isotope(), Some(13));
    }

    #[test]
    fn test_parse_organic_symbol() {
        let result = parse_element_atom("Cl123");
        assert!(result.is_ok());
        let (remaining, atom) = result.unwrap();
        assert_eq!(remaining, "123");
        assert_eq!(atom.element, "Cl");
        assert!(!atom.aromatic);
        assert_eq!(atom.charge, 0);
        assert_eq!(atom.isotope, None);
        assert_eq!(atom.h_count, None);

        let result = parse_element_atom("BrOH");
        assert!(result.is_ok());
        let (remaining, atom) = result.unwrap();
        assert_eq!(remaining, "OH");
        assert_eq!(atom.element, "Br");
        assert!(!atom.aromatic);
    }

    #[test]
    fn test_parse_h_count() {
        // Test 1: H with digit (H2)
        let result = parse_h_count("H2Random");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "Random");
        assert_eq!(h_count, 2);

        // Test 2: H with multiple digits (H34)
        let result = parse_h_count("H43Y");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "Y");
        assert_eq!(h_count, 43);

        // Test 3: H without digit (defaults to 1)
        let result = parse_h_count("HY");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "Y");
        assert_eq!(h_count, 1);

        // Test 4: H without digit followed by parenthesis
        let result = parse_h_count("H(OH)");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "(OH)");
        assert_eq!(h_count, 1);

        // Test 5: H with single digit
        let result = parse_h_count("H5");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(h_count, 5);

        // Test 6: H at end of string (no digit)
        let result = parse_h_count("H");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(h_count, 1);

        // Test 7: H with digit followed by special char
        let result = parse_h_count("H3+2");
        assert!(result.is_ok());
        let (remaining, h_count) = result.unwrap();
        assert_eq!(remaining, "+2");
        assert_eq!(h_count, 3);

        // Test 8: Should fail - no H at start
        let result = parse_h_count("2H");
        assert!(result.is_err());

        let result = parse_h_count("CH3");
        assert!(result.is_err());

        // Test 9: Should fail - empty input
        let result = parse_h_count("");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_charge() {
        // Test 1: Positive charge with digit (+2)
        let result = parse_charge("+2]");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "]");
        assert_eq!(charge, 2);

        // Test 2: Positive charge without digit (defaults to +1)
        let result = parse_charge("+]");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "]");
        assert_eq!(charge, 1);

        // Test 3: Negative charge with digit (-3)
        let result = parse_charge("-3]");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "]");
        assert_eq!(charge, -3);

        // Test 4: Negative charge without digit (defaults to -1)
        let result = parse_charge("-]");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "]");
        assert_eq!(charge, -1);

        // Test 5: Double positive (++)
        let result = parse_charge("++]");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "]");
        assert_eq!(charge, 2);

        // Test 6: Double negative (--)
        let result = parse_charge("--]");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "]");
        assert_eq!(charge, -2);

        // Test 7: Large positive charge
        let result = parse_charge("+5H");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "H");
        assert_eq!(charge, 5);

        // Test 8: Large negative charge
        let result = parse_charge("-4H");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "H");
        assert_eq!(charge, -4);

        // Test 9: Positive charge at end
        let result = parse_charge("+");
        assert!(result.is_ok());
        let (remaining, charge) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(charge, 1);

        // Test 10: Should fail - no charge symbol
        let result = parse_charge("2");
        assert!(result.is_err());

        let result = parse_charge("H");
        assert!(result.is_err());

        // Test 11: Should fail - empty input
        let result = parse_charge("");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_bond() {
        // Test 1: Single bond
        let result = parse_bond("-C");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "C");
        assert_eq!(bond, BondToken::Single);

        // Test 2: Double bond
        let result = parse_bond("=O");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "O");
        assert_eq!(bond, BondToken::Double);

        // Test 3: Triple bond
        let result = parse_bond("#N");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "N");
        assert_eq!(bond, BondToken::Triple);

        // Test 4: Aromatic bond
        let result = parse_bond(":c");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "c");
        assert_eq!(bond, BondToken::Aromatic);

        // Test 5: Single bond at end
        let result = parse_bond("-");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(bond, BondToken::Single);

        // Test 6: Double bond followed by parenthesis
        let result = parse_bond("=(C)");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "(C)");
        assert_eq!(bond, BondToken::Double);

        // Test 7: Triple bond followed by bracket
        let result = parse_bond("#[NH3+]");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "[NH3+]");
        assert_eq!(bond, BondToken::Triple);

        // Test 8: Aromatic bond followed by digit
        let result = parse_bond(":1");
        assert!(result.is_ok());
        let (remaining, bond) = result.unwrap();
        assert_eq!(remaining, "1");
        assert_eq!(bond, BondToken::Aromatic);

        // Test 9: Should fail - invalid bond character
        let result = parse_bond("~");
        assert!(result.is_err());

        let result = parse_bond("C");
        assert!(result.is_err());

        // Test 10: Should fail - empty input
        let result = parse_bond("");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_ring_number() {
        // Test 1: Single digit ring number
        let result = parse_ring_number("1ccccc");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "ccccc");
        assert_eq!(ring, 1);

        // Test 2: Another single digit
        let result = parse_ring_number("5C");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "C");
        assert_eq!(ring, 5);

        // Test 3: Two-digit ring number with % prefix
        let result = parse_ring_number("%10C");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "C");
        assert_eq!(ring, 10);

        // Test 4: Large ring number with %
        let result = parse_ring_number("%99cc");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "cc");
        assert_eq!(ring, 99);

        // Test 5: Ring number at end
        let result = parse_ring_number("3");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(ring, 3);

        // Test 6: % ring number at end
        let result = parse_ring_number("%15");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(ring, 15);

        // Test 7: a bare label is one digit, so "123" is three labels and this
        // parser takes only the first. The assertions here used to say `123`
        // and `""` while the comment above them said "only takes first digit"
        // — the test was written to match the bug rather than the spec it
        // describes, which is why #69 survived having a test at all.
        let result = parse_ring_number("123");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "23");
        assert_eq!(ring, 1);

        // Test 8: Ring number followed by parenthesis
        let result = parse_ring_number("2(C)");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "(C)");
        assert_eq!(ring, 2);

        // Test 9: % with three digits
        let result = parse_ring_number("%100X");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, "X");
        assert_eq!(ring, 100);

        // Test 10: Should fail - no digit
        let result = parse_ring_number("C");
        assert!(result.is_err());

        let result = parse_ring_number("(");
        assert!(result.is_err());

        // Test 11: Should fail - % without digit
        let result = parse_ring_number("%C");
        assert!(result.is_err());

        // Test 12: Should fail - empty input
        let result = parse_ring_number("");
        assert!(result.is_err());
    }

    /// A bare ring label is one digit, so an atom can close two rings by
    /// naming them back to back. That is not an exotic spelling — it is the
    /// shared atom of any fused ring system, which is why #69 cost real
    /// molecules.
    #[test]
    fn test_an_atom_can_close_two_rings_at_once() {
        // The minimal case: two fused six-rings sharing an edge, no aromatics
        // involved. 10 atoms, 11 bonds — one more bond than atoms is what
        // makes it two rings rather than one.
        let mol = parse_smiles("C12CCCCC1CCCC2").expect("bicyclic should parse");
        assert_eq!(mol.num_atoms(), 10);
        assert_eq!(mol.num_bonds(), 11);
    }

    #[test]
    fn test_the_molecules_69_lost_now_parse() {
        // Each pair is the same ring system spelled two ways. Before the fix
        // the first of each pair failed and the second worked, purely because
        // of where the closure digits happened to land.
        for (adjacent, spaced, atoms, bonds) in [
            ("c1ccnc2ccccc12", "c1ccc2ncccc2c1", 10, 11),
            ("c1ccc2c(c1)ccc1ccccc12", "c1ccc2c(c1)ccc1c2cccc1", 14, 16),
        ] {
            let a = parse_smiles(adjacent).unwrap_or_else(|e| panic!("{adjacent}: {e}"));
            let b = parse_smiles(spaced).unwrap_or_else(|e| panic!("{spaced}: {e}"));
            assert_eq!((a.num_atoms(), a.num_bonds()), (atoms, bonds), "{adjacent}");
            // The point of the pairing: two spellings of one ring system must
            // agree, not merely both succeed.
            assert_eq!(
                (a.num_atoms(), a.num_bonds()),
                (b.num_atoms(), b.num_bonds()),
                "{adjacent} and {spaced} describe the same ring system"
            );
        }
    }

    #[test]
    fn test_percent_still_reads_a_multi_digit_label() {
        // The bare branch got narrower; `%NN` must not have.
        let mol = parse_smiles("C%10CCCC%10").expect("%10 should still close");
        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 5);
    }

    #[test]
    fn test_an_oversized_label_errors_instead_of_panicking() {
        // Found while fixing #69, and worse than a wrong answer: `unwrap()` on
        // a u32 overflow aborted the process. Once a file reader exists, one
        // malformed line killing the whole run is the difference between a
        // skipped record and no output at all.
        let result = parse_smiles("C%99999999999999999999CC");
        assert!(
            result.is_err(),
            "must not panic, and must not succeed either"
        );
    }

    /// Bond and branch characters are legal tokens on their own, so a string
    /// of nothing but those used to tokenize cleanly and build a molecule with
    /// no atoms — returned as success. #151.
    #[test]
    fn test_a_string_with_no_atoms_is_an_error() {
        for input in ["$", "$$$$", "=", "#", "-", ":", "====", "(", ""] {
            let result = parse_smiles(input);
            assert!(
                matches!(result, Err(SmilesError::NoAtoms)),
                "{input:?} should be NoAtoms, got {:?}",
                result.map(|m| m.num_atoms())
            );
        }
    }

    #[test]
    fn test_the_sdf_terminator_is_the_case_that_mattered() {
        // `$$$$` ends every SDF record, so reading an SDF as SMILES used to
        // report one atomless molecule per record. That defeated the "did
        // anything parse?" check a caller relies on to notice a wrong format.
        assert!(parse_smiles("$$$$").is_err());
    }

    #[test]
    fn test_every_unclosed_ring_is_reported_and_sorted() {
        // It used to name one arbitrary HashMap key, so the same input reported
        // a different ring between runs (#153). Both are genuinely open, and
        // naming one tells the reader least.
        let Err(SmilesError::UnclosedRings(open)) = parse_smiles("C1CC2") else {
            panic!("expected UnclosedRings");
        };
        assert_eq!(open, vec![1, 2]);

        // The message is what a user actually sees, and sorting is only
        // observable through it — an unsorted Vec would still contain both.
        let message = parse_smiles("C1CC2").unwrap_err().to_string();
        assert_eq!(message, "Unclosed rings: 1, 2");
    }

    #[test]
    fn test_a_bond_with_nothing_to_attach_to_is_rejected() {
        // Both used to parse as ethane — a real molecule, quietly different
        // from the one written, which is worse than a rejection (#155).
        for input in ["=CC", "CC=", "#CC", "CCO="] {
            assert!(
                matches!(parse_smiles(input), Err(SmilesError::DanglingBond)),
                "{input:?} should be DanglingBond, got {:?}",
                parse_smiles(input).map(|m| (m.num_atoms(), m.num_bonds()))
            );
        }
    }

    #[test]
    fn test_real_bonds_are_untouched() {
        // The other direction: a check this aggressive could reject valid
        // input, and every one of these has a bond symbol in it.
        for input in ["C=C", "C#N", "CC(=O)O", "c1ccccc1-c1ccccc1", "C=1CCCCC1"] {
            assert!(parse_smiles(input).is_ok(), "{input:?} should parse");
        }
    }

    #[test]
    fn test_an_atomless_string_still_reports_no_atoms() {
        // Ordering, asserted rather than left to `test_a_string_with_no_atoms_
        // is_an_error` to notice. `$` is the quadruple-bond token, so `$$$$` —
        // the SDF record terminator — is four dangling bonds and no atoms.
        // Both complaints are true; "No atoms in SMILES" names the mistake
        // someone actually made, which is reading an SDF as SMILES.
        for input in ["$$$$", "=", "===="] {
            assert!(
                matches!(parse_smiles(input), Err(SmilesError::NoAtoms)),
                "{input:?} should prefer NoAtoms over DanglingBond"
            );
        }
    }

    #[test]
    fn test_the_more_specific_errors_still_win() {
        // The atom-count check is last on purpose: a molecule that also has an
        // unclosed ring or an unbalanced branch should say so, since that is
        // the more useful diagnosis.
        assert!(matches!(
            parse_smiles("1"),
            Err(SmilesError::InvalidRing(1))
        ));
        assert!(matches!(
            parse_smiles("()"),
            Err(SmilesError::MismatchedBranches)
        ));
    }

    #[test]
    fn test_real_molecules_are_unaffected() {
        // The narrowest possible molecule, to make sure the check is "no
        // atoms" and not "too few atoms".
        for input in ["C", "O", "[H]", "CCO", "c1ccccc1"] {
            let mol = parse_smiles(input).unwrap_or_else(|e| panic!("{input}: {e}"));
            assert!(mol.num_atoms() > 0, "{input}");
        }
    }

    /// An unspecified bond is aromatic only when **both** atoms are. Checking
    /// only the atom being added made the answer depend on which end the author
    /// wrote first, which is #167.
    #[test]
    fn test_a_substituent_bond_is_single_whichever_end_is_written_first() {
        for (smiles, aromatic_bonds) in [
            ("c1ccccc1O", 6),
            ("Oc1ccccc1", 6),
            ("c1ccc(O)cc1", 6),
            ("Cc1ccccc1", 6),
            ("c1ccccc1C", 6),
            ("CCc1ccccc1", 6),
        ] {
            let mol = parse_smiles(smiles).expect(smiles);
            let found = (0..mol.num_bonds())
                .filter(|&i| mol.bond(i).is_aromatic())
                .count();
            assert_eq!(found, aromatic_bonds, "{smiles}");
        }
    }

    #[test]
    fn test_no_aromatic_bond_has_a_non_aromatic_atom() {
        // The invariant, rather than a list of cases. It is what would have
        // caught both this and #94, and it keeps catching whatever comes next.
        for smiles in [
            "c1ccccc1O",
            "Oc1ccccc1",
            "Cc1ccccc1O",
            "CC(C)c1ccc(O)cc1",
            "c1ccc2ccccc2c1",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
            "CC(=O)Oc1ccccc1C(=O)O",
        ] {
            let mol = parse_smiles(smiles).expect(smiles);
            for i in 0..mol.num_bonds() {
                let bond = mol.bond(i);
                if bond.is_aromatic() {
                    assert!(
                        mol.atom(bond.atom1()).is_aromatic()
                            && mol.atom(bond.atom2()).is_aromatic(),
                        "{smiles}: bond {i} is aromatic but an endpoint is not"
                    );
                }
            }
        }
    }

    #[test]
    fn test_two_spellings_of_one_molecule_fingerprint_identically() {
        // The symptom that made #167 visible, asserted where it actually
        // matters. Equal bond counts would pass on a molecule whose bonds were
        // the wrong *kind*; equal fingerprints would not.
        use crate::fp::morgan::MorganFingerprint;

        let fp = |smiles: &str| {
            let mol = parse_smiles(smiles).expect(smiles);
            MorganFingerprint::get_fingerprint_as_bitvec(
                &mol, 2, 2048, None, None, false, true, false,
            )
            .expect("fingerprint")
        };

        let reference = fp("c1ccccc1O");
        for other in ["Oc1ccccc1", "c1ccc(O)cc1"] {
            assert_eq!(fp(other), reference, "{other} should match c1ccccc1O");
        }

        let toluene = fp("Cc1ccccc1");
        assert_eq!(fp("c1ccccc1C"), toluene);
        // ...and a different molecule must still differ, or the test above
        // would pass on a fingerprint function that returned a constant.
        assert_ne!(toluene, reference);
    }

    #[test]
    fn test_a_biaryl_bond_needs_its_explicit_single() {
        // Where the rule could look like it over-fires, and does not. Two
        // aromatic rings joined by a single bond must say so: SMILES treats an
        // unspecified bond between aromatic atoms as aromatic, which is exactly
        // why biphenyl is written with the dash.
        let explicit = parse_smiles("c1ccccc1-c1ccccc1").expect("biphenyl");
        assert_eq!(
            (0..explicit.num_bonds())
                .filter(|&i| explicit.bond(i).is_aromatic())
                .count(),
            12,
            "the linking bond must stay single"
        );

        let implicit = parse_smiles("c1ccccc1c1ccccc1").expect("dash-less");
        assert_eq!(
            (0..implicit.num_bonds())
                .filter(|&i| implicit.bond(i).is_aromatic())
                .count(),
            13,
            "without the dash the linking bond is aromatic, per the spec"
        );
    }

    #[test]
    fn test_tokenize_smiles() {
        // Test 1: Simple molecule - methane (C)
        let result = tokenize_smiles("C");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 1);
        match &tokens[0] {
            Token::Atom(atom) => {
                assert_eq!(atom.element, "C");
                assert!(!atom.aromatic);
            }
            _ => panic!("Expected Atom token"),
        }

        // Test 2: Ethane with bond (C-C)
        let result = tokenize_smiles("C-C");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 3);
        assert!(matches!(tokens[0], Token::Atom(_)));
        assert!(matches!(tokens[1], Token::Bond(BondToken::Single)));
        assert!(matches!(tokens[2], Token::Atom(_)));

        // Test 3: Benzene ring (c1ccccc1)
        let result = tokenize_smiles("c1ccccc1");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 8); // 6 atoms + 2 ring numbers
        assert!(matches!(tokens[0], Token::Atom(_)));
        assert!(matches!(tokens[1], Token::Ring(1)));

        // Test 4: Branched molecule (C(C)C)
        let result = tokenize_smiles("C(C)C");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 5);
        assert!(matches!(tokens[0], Token::Atom(_)));
        assert!(matches!(tokens[1], Token::Branch));
        assert!(matches!(tokens[2], Token::Atom(_)));
        assert!(matches!(tokens[3], Token::BranchEnd));
        assert!(matches!(tokens[4], Token::Atom(_)));

        // Test 5: Molecule with double bond (C=O)
        let result = tokenize_smiles("C=O");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 3);
        assert!(matches!(tokens[1], Token::Bond(BondToken::Double)));

        // Test 6: Bracket atom with charge ([NH4+])
        let result = tokenize_smiles("[NH4+]");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 1);
        match &tokens[0] {
            Token::Atom(atom) => {
                assert_eq!(atom.element, "N");
                assert_eq!(atom.h_count, Some(4));
                assert_eq!(atom.charge, 1);
            }
            _ => panic!("Expected Atom token"),
        }

        // Test 7: Complex molecule (CC(=O)O)
        let result = tokenize_smiles("CC(=O)O");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 7);

        // Test 8: Aromatic with bond (c:c)
        let result = tokenize_smiles("c:c");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 3);
        assert!(matches!(tokens[1], Token::Bond(BondToken::Aromatic)));

        // Test 9: Ring with % notation (C%10CCCCC%10)
        let result = tokenize_smiles("C%10C");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 3);
        assert!(matches!(tokens[1], Token::Ring(10)));

        // Test 10: Multiple bonds and branches (C#CC(C)C)
        let result = tokenize_smiles("C#CC(C)C");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert!(tokens.len() > 5);

        // Test 11: Empty string
        let result = tokenize_smiles("");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 0);

        // Test 12: Should fail - invalid character
        let result = tokenize_smiles("C@C");
        assert!(result.is_err());
        match result {
            Err(SmilesError::ParseError(msg)) => {
                assert!(msg.contains("Unexpected characters"));
            }
            _ => panic!("Expected ParseError"),
        }

        // Test 13: Should fail - unclosed bracket
        let result = tokenize_smiles("[NH4+");
        assert!(result.is_err());

        // Test 14: Should fail - invalid element
        let result = tokenize_smiles("X");
        assert!(result.is_err());

        // Test 15: Chlorine and Bromine
        let result = tokenize_smiles("ClBr");
        assert!(result.is_ok());
        let tokens = result.unwrap();
        assert_eq!(tokens.len(), 2);
    }

    #[test]
    fn test_parse_element_symbol() {
        // Test 1: Two-letter elements (uppercase) - should match first alternative
        let result = parse_element_symbol("Cl123");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "123");
        assert_eq!(element, "Cl");

        let result = parse_element_symbol("Br(OH)");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "(OH)");
        assert_eq!(element, "Br");

        // Test 3: Single letter elements - should match second alternative
        let result = parse_element_symbol("C");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(element, "C");

        let result = parse_element_symbol("N123");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "123");
        assert_eq!(element, "N");

        // Test 4: Two-letter element symbols (second alternative)
        let result = parse_element_symbol("Ca123");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "123");
        assert_eq!(element, "Ca");

        let result = parse_element_symbol("Fe+2");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "+2");
        assert_eq!(element, "Fe");

        // Test 5: Priority test - "Cl" should match as whole, not just "C"
        let result = parse_element_symbol("Cl");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(element, "Cl");

        // Test 6: Single aromatic elements
        let result = parse_element_symbol("c1");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "1");
        assert_eq!(element, "c");

        let result = parse_element_symbol("n");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(element, "n");

        // Test 7: Two-letter combinations with lowercase second char
        let result = parse_element_symbol("Si(");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "(");
        assert_eq!(element, "Si");

        let result = parse_element_symbol("Br");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(element, "Br");

        // Test 8: Stops at non-matching characters
        let result = parse_element_symbol("C(OH)");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "(OH)");
        assert_eq!(element, "C");

        let result = parse_element_symbol("O=C");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "=C");
        assert_eq!(element, "O");

        // Test 9: Invalid input - should fail
        let result = parse_element_symbol("[NH3]");
        assert!(
            result.is_err(),
            "Should fail when starting with non-element character"
        );

        let result = parse_element_symbol("(OH)");
        assert!(
            result.is_err(),
            "Should fail when starting with parenthesis"
        );

        let result = parse_element_symbol("");
        assert!(result.is_err(), "Should fail on empty input");

        let result = parse_element_symbol("X");
        assert!(result.is_err(), "Should fail on invalid element 'X'");

        let result = parse_element_symbol("123");
        assert!(result.is_err(), "Should fail when starting with digit");

        // Test 10: Case sensitivity verification
        let result_upper = parse_element_symbol("Br");
        assert!(result_upper.is_ok());
        assert_eq!(result_upper.unwrap().1, "Br");

        // Test 11: Elements that could be ambiguous
        let result = parse_element_symbol("C3");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "3");
        assert_eq!(element, "C");
    }

    #[test]
    fn test_aromatic_ring_closure_bond_is_aromatic() {
        // The ring-closure bond used to default to Single regardless of
        // aromaticity, so benzene came out with five aromatic bonds and one
        // single one -- visible as a missing inner line when drawn, and, since
        // Morgan invariants weight bond order, a wrong fingerprint.
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.num_bonds(), 6);
        for bond in mol.bonds() {
            assert_eq!(
                bond.order(),
                BondOrder::Aromatic,
                "bond {}-{} should be aromatic",
                bond.atom1(),
                bond.atom2()
            );
            assert!(bond.is_aromatic());
        }
    }

    #[test]
    fn test_fused_aromatic_ring_closures_are_aromatic() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        for bond in mol.bonds() {
            assert_eq!(bond.order(), BondOrder::Aromatic);
        }
    }

    #[test]
    fn test_non_aromatic_ring_closure_stays_single() {
        // Cyclohexane's atoms aren't aromatic, so its closure must not be
        // promoted -- the inference keys off the atoms, not off being a ring.
        let mol = parse_smiles("C1CCCCC1").unwrap();
        for bond in mol.bonds() {
            assert_eq!(bond.order(), BondOrder::Single);
        }
    }

    #[test]
    fn test_kekule_ring_closure_stays_single() {
        // Kekule benzene alternates explicitly; the closing bond is the single
        // one of the alternation and must stay that way.
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        let closure = mol.bonds().last().unwrap();
        assert_eq!(closure.order(), BondOrder::Single);
    }

    #[test]
    fn test_explicit_ring_closure_bond_symbol_wins() {
        // An explicit symbol on the closure must override the inferred order.
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        assert_eq!(mol.num_bonds(), 6);

        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.num_bonds(), 6);
    }

    #[test]
    fn test_equivalent_smiles_give_equal_bond_orders() {
        // Phenol two ways. The same molecule must parse to the same multiset of
        // bond orders however it was written -- when it didn't, the two forms
        // fingerprinted differently and scored 0.27 Tanimoto against each other.
        let sort_orders = |smi: &str| {
            let mol = parse_smiles(smi).unwrap();
            let mut v: Vec<String> = mol
                .bonds()
                .iter()
                .map(|b| format!("{:?}", b.order()))
                .collect();
            v.sort();
            v
        };
        assert_eq!(sort_orders("c1ccccc1O"), sort_orders("c1ccc(O)cc1"));
    }

    #[test]
    fn test_quadruble() {
        let mol = parse_smiles("[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl").unwrap();
        assert_eq!(mol.num_atoms(), 10);
        assert_eq!(mol.formula(), "Cl8Rh2");
    }
}
