use std::collections::HashMap;

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{char, digit1, one_of},
    combinator::{map, opt, recognize},
    multi::many0,
    sequence::{delimited, preceded},
    IResult, Parser,
};

use crate::errors::SmilesError;
use chemcore::{
    atom::{Atom, Element},
    bond::{Bond, BondOrder},
    molecule::Molecule,
};

#[derive(Debug, Clone)]
enum Token {
    Atom(AtomToken),
    Bond(BondToken),
    Branch,
    BranchEnd,
    Ring(u32),
}

#[derive(Debug, Clone)]
struct AtomToken {
    element: String,
    aromatic: bool,
    charge: i8,
    isotope: Option<u16>,
    h_count: Option<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BondToken {
    Single,
    Double,
    Triple,
    Quadruble,
    Aromatic,
}

impl BondToken {
    fn to_bond_order(&self) -> BondOrder {
        match self {
            BondToken::Single => BondOrder::Single,
            BondToken::Double => BondOrder::Double,
            BondToken::Triple => BondOrder::Triple,
            BondToken::Quadruble => BondOrder::Quadruple,
            BondToken::Aromatic => BondOrder::Aromatic,
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
                    "Unexpected characters: {}",
                    remaining
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
        map(parse_organic_atom, Token::Atom),
        map(parse_bond, Token::Bond),
        map(char('('), |_| Token::Branch),
        map(char(')'), |_| Token::BranchEnd),
        map(parse_ring_number, Token::Ring),
    )))
    .parse(input)
}

fn parse_organic_atom(input: &str) -> IResult<&str, AtomToken> {
    let (remaining_input, symbol) = alt((
        recognize(tag("Cl")),
        recognize(tag("Br")),
        recognize(one_of("BCNOPSFIbcnops")),
    ))
    .parse(input)?;

    let aromatic = symbol.chars().next().unwrap().is_lowercase();
    let element = if aromatic {
        symbol.to_uppercase()
    } else {
        symbol.to_string()
    };

    Ok((
        remaining_input,
        AtomToken {
            element,
            aromatic,
            charge: 0,
            isotope: None,
            h_count: None,
        },
    ))
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
                },
            ))
        },
        char(']'),
    )
    .parse(input)
}

fn parse_element_symbol(input: &str) -> IResult<&str, String> {
    alt((
        // First Alternative
        map(
            alt((tag("Cl"), tag("Br"), tag("cl"), tag("br"))), // element must match exactly one of the 4
            |s: &str| s.to_string(),
        ),
        // Second Alternative
        map(
            recognize(nom::sequence::pair(
                // Example: Fe ==> (F): Matched by 1st Parser || (e): Matched by 2nd Parser
                one_of("BCNOFPSIbcnofpsi"), // Matched 1st Char
                opt(one_of("abcdefghijklmnopqrstuvwxyz")), // Matches 2nd Char
            )),
            |s: &str| s.to_string(),
        ),
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
    ))
    .parse(input)
}

fn parse_ring_number(input: &str) -> IResult<&str, u32> {
    alt((
        map(preceded(char('%'), digit1), |s: &str| s.parse().unwrap()),
        map(digit1, |s: &str| s.parse().unwrap()),
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
    let mut rings: HashMap<u32, (usize, Option<BondToken>)> = HashMap::new();

    for token in tokens {
        match token {
            Token::Atom(atom_token) => {
                let atom_idx = add_atom_from_token(&mut mol, atom_token)?;

                if let Some(prev_atom) = current_atom {
                    let bond_order = next_bond // This is the previous pending_bond that will be applied NOW.
                        .map(|b| b.to_bond_order())
                        .unwrap_or(if atom_token.aromatic {
                            BondOrder::Aromatic
                        } else {
                            BondOrder::Single
                        });

                    let mut bond = Bond::new(prev_atom, atom_idx, bond_order);
                    if bond_order == BondOrder::Aromatic {
                        bond = bond.with_aromatic(true);
                    }
                    mol.add_bond(bond)
                        .map_err(|e| SmilesError::ParseError(e.to_string()))?;
                }

                current_atom = Some(atom_idx);
                next_bond = None; // pending_bond: i.e. it will be applied once next atom arrives
            }

            Token::Bond(bond_token) => {
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

            Token::Ring(ring_num) => {
                let current = current_atom.ok_or(SmilesError::InvalidRing(*ring_num))?;

                if let Some((ring_start, ring_bond)) = rings.remove(ring_num) {
                    let bond_order = next_bond
                        .or(ring_bond)
                        .map(|b| b.to_bond_order())
                        .unwrap_or(BondOrder::Single);

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
        return Err(SmilesError::UnclosedRing(*rings.keys().next().unwrap()));
    }

    if !stack.is_empty() {
        return Err(SmilesError::MismatchedBranches);
    }

    mol.calculate_implicit_hydrogens();
    Ok(mol)
}

fn add_atom_from_token(mol: &mut Molecule, token: &AtomToken) -> Result<usize, SmilesError> {
    // FIX:Hard coded... solve this
    let atomic_number = match token.element.as_str() {
        "H" => 1,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        "B" => 5,
        _ => return Err(SmilesError::InvalidElement(token.element.clone())),
    };

    let element = Element::new(atomic_number)
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

    let idx = mol.add_atom(atom);

    if let Some(h_count) = token.h_count {
        mol.atom_mut(idx).set_explicit_hydrogens(h_count);
    }

    Ok(idx)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let result = parse_organic_atom("Cl123");
        assert!(result.is_ok());
        let (remaining, atom) = result.unwrap();
        assert_eq!(remaining, "123");
        assert_eq!(atom.element, "Cl");
        assert!(!atom.aromatic);
        assert_eq!(atom.charge, 0);
        assert_eq!(atom.isotope, None);
        assert_eq!(atom.h_count, None);

        let result = parse_organic_atom("BrOH");
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

        // Test 7: Multiple digits without % (only takes first digit)
        let result = parse_ring_number("123");
        assert!(result.is_ok());
        let (remaining, ring) = result.unwrap();
        assert_eq!(remaining, ""); // Only consumes "1"
        assert_eq!(ring, 123);

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

        // Test 2: Two-letter elements (lowercase) - should match first alternative
        let result = parse_element_symbol("cl2");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "2");
        assert_eq!(element, "cl");

        let result = parse_element_symbol("br-");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "-");
        assert_eq!(element, "br");

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
        let result_lower = parse_element_symbol("br");
        assert!(result_upper.is_ok());
        assert!(result_lower.is_ok());
        assert_eq!(result_upper.unwrap().1, "Br");
        assert_eq!(result_lower.unwrap().1, "br");

        // Test 11: Elements that could be ambiguous
        let result = parse_element_symbol("C3");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "3");
        assert_eq!(element, "C");

        let result = parse_element_symbol("cl");
        assert!(result.is_ok());
        let (remaining, element) = result.unwrap();
        assert_eq!(remaining, "");
        assert_eq!(element, "cl");
    }
}
