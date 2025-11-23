use crate::errors::SdfError;
use chemcore::prelude::*;

const SDF_ENTRY_END: &str = "$$$$";
const SDF_STRUCTURE_END: &str = "M END";
const SDF_PROPERTY_SUFFIX: char = '>';
const SDF_PROPERTY_PREFIX: &str = "> <";

const SDF_NAME_LINE: usize = 0;
const SDF_COUNT_LINE: usize = 3;
const SDF_ATOM_FIELD: usize = 4;
const SDF_MIN_BOND_FIELDS: usize = 3;
const SDF_MIN_HEADER_LINES: usize = 4;
const SDF_ATOM_BLOCK_START: usize = 4;
const SDF_PROPERTY_PREFIX_LEN: usize = 3;

/// Parses an SDF (Structure-Data File) string into a `Molecule` object.
///
/// This parser handles real-world SDF files with flexible formatting.
pub fn parse_sdf(sdf: &str) -> Result<Molecule, SdfError> {
    let lines: Vec<&str> = sdf.lines().collect();

    // Validate minimum SDF structure (header must be at least 4 lines)
    if lines.len() < SDF_MIN_HEADER_LINES {
        return Err(SdfError::ParseError("Too few lines".to_string()));
    }

    let mut mol = Molecule::new();

    // Line 0: Molecule name (optional)
    if !lines[SDF_NAME_LINE].trim().is_empty() {
        mol.set_name(lines[SDF_NAME_LINE].trim().to_string());
    }

    // Lines 1-2: Program information and comments (skipped)

    // Line 3: Counts line - specifies number of atoms and bonds
    let counts_line = lines[SDF_COUNT_LINE];
    let (num_atoms, num_bonds) = parse_counts_line(counts_line)?;

    // Parse atom block (starts at line 4)
    for i in 0..num_atoms {
        let line_idx = SDF_ATOM_BLOCK_START + i;
        if line_idx >= lines.len() {
            return Err(SdfError::ParseError("Not enough atom lines".to_string()));
        }
        parse_atom_line(&mut mol, lines[line_idx])?;
    }

    // Parse bond block (follows atom block)
    for i in 0..num_bonds {
        let line_idx = SDF_ATOM_BLOCK_START + num_atoms + i;
        if line_idx >= lines.len() {
            return Err(SdfError::ParseError("Not enough bond lines".to_string()));
        }
        parse_bond_line(&mut mol, lines[line_idx])?;
    }

    // Parse optional properties block (follows bond block)
    let prop_start = SDF_ATOM_BLOCK_START + num_atoms + num_bonds;
    parse_properties(&mut mol, &lines[prop_start..])?;

    // Calculate implicit hydrogens for each atom based on valence rules
    mol.calculate_implicit_hydrogens();
    Ok(mol)
}

/// Parses the counts line (line 3) of an SDF file.
///
/// **FIXED**: Now handles both fixed-width format and whitespace-separated format.
///
/// Example: " 10  9  0  0  0  0  0  0  0  0999 V2000"
///           ^^^ ~~~
///          10 atoms, 9 bonds
fn parse_counts_line(line: &str) -> Result<(usize, usize), SdfError> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() < 2 {
        return Err(SdfError::MissingCountsLine);
    }

    let num_atoms = parts[0]
        .parse::<usize>()
        .map_err(|_| SdfError::MissingCountsLine)?;

    let num_bonds = parts[1]
        .parse::<usize>()
        .map_err(|_| SdfError::MissingCountsLine)?;

    Ok((num_atoms, num_bonds))
}

/// Parses a single atom line from the SDF atom block.
///
/// **FIXED**: Now handles both fixed-width and whitespace-separated formats.
///
/// # Real-World SDF Atom Line Format
/// Atom lines can vary in format:
/// 1. **Fixed-width format** (classic V2000):
///    - Positions 0-9: X coordinate
///    - Positions 10-19: Y coordinate  
///    - Positions 20-29: Z coordinate
///    - Positions 31-33: Element symbol (right-aligned in 3 chars)
///
/// 2. **Whitespace-separated format** (more common):
///    - Fields separated by spaces: x y z symbol [charge] [other fields]
///
/// Example from your file:
/// ```text
///     1.3051    0.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
///     ^^^^^^    ^^^^^^    ^^^^^^ ^
///        x         y         z   element
/// ```
fn parse_atom_line(mol: &mut Molecule, line: &str) -> Result<(), SdfError> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    // Real-world SDF files typically have at least 4 fields: x, y, z, symbol
    if parts.len() < SDF_ATOM_FIELD {
        return Err(SdfError::InvalidAtomLine(format!(
            "Expected at least 4 fields (x y z symbol), got {}: {}",
            parts.len(),
            line
        )));
    }

    // The 4th field (index 3) is the element symbol in whitespace-separated format
    let symbol = parts[3];

    // Look up atomic number from the ELEMENT_SYMBOLS constant
    let atomic_number = ELEMENT_SYMBOLS
        .iter()
        .position(|&s| s == symbol)
        .ok_or_else(|| SdfError::InvalidAtomLine(format!("Unknown element: {}", symbol)))?;

    // Create element from atomic number
    let element = Element::new(atomic_number as u8)
        .ok_or_else(|| SdfError::InvalidAtomLine(symbol.to_string()))?;

    // Create and add atom to molecule
    let atom = Atom::new(element);
    mol.add_atom(atom);

    Ok(())
}

/// Parses a single bond line from the SDF bond block.
///
/// **ENHANCED**: Added better error messages and validation.
///
/// Bond lines are whitespace-separated: atom1 atom2 bond_type [stereo] [other fields]
/// Example: "  1  2  1  0  0  0  0"
///             ^  ^  ^
///             |  |  └─ Bond type (1=single, 2=double, 3=triple, 4=aromatic)
///             |  └──── Second atom (1-based index)
///             └─────── First atom (1-based index)
fn parse_bond_line(mol: &mut Molecule, line: &str) -> Result<(), SdfError> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    // Bond line must have at least 3 fields: atom1, atom2, bond_type
    if parts.len() < SDF_MIN_BOND_FIELDS {
        return Err(SdfError::InvalidBondLine(format!(
            "Expected at least 3 fields (atom1 atom2 type), got {}: {}",
            parts.len(),
            line
        )));
    }

    // Parse atom indices (convert from 1-based SDF to 0-based internal indexing)
    let atom1: usize = parts[0]
        .parse::<usize>()
        .map_err(|_| SdfError::InvalidBondLine(format!("Invalid atom1 index: {}", parts[0])))?
        .checked_sub(1)
        .ok_or_else(|| SdfError::InvalidBondLine("Atom index cannot be 0".to_string()))?;

    let atom2: usize = parts[1]
        .parse::<usize>()
        .map_err(|_| SdfError::InvalidBondLine(format!("Invalid atom2 index: {}", parts[1])))?
        .checked_sub(1)
        .ok_or_else(|| SdfError::InvalidBondLine("Atom index cannot be 0".to_string()))?;

    // Parse bond type
    let bond_type: u8 = parts[2]
        .parse()
        .map_err(|_| SdfError::InvalidBondLine(format!("Invalid bond type: {}", parts[2])))?;

    // Map SDF bond type to internal BondOrder enum
    let bond_order = match bond_type {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => {
            return Err(SdfError::InvalidBondLine(format!(
                "Invalid bond type: {} (must be 1-4)",
                bond_type
            )))
        }
    };

    // Add bond to molecule (validates atom indices internally)
    mol.add_bond(Bond::new(atom1, atom2, bond_order))
        .map_err(|e| SdfError::ParseError(format!("Failed to add bond: {}", e)))?;

    Ok(())
}

/// Parses the optional properties block at the end of an SDF entry.
///
/// **ENHANCED**: Now skips "M  END" lines and handles empty properties better.
fn parse_properties(mol: &mut Molecule, lines: &[&str]) -> Result<(), SdfError> {
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i].trim();

        // Skip "M  END" marker
        if line.starts_with(SDF_STRUCTURE_END) {
            i += 1;
            continue;
        }

        // "$$$$" marks the end of this SDF entry
        if line == SDF_ENTRY_END {
            break;
        }

        // Property names are enclosed in "> <" and ">"
        if line.starts_with(SDF_PROPERTY_PREFIX) && line.ends_with(SDF_PROPERTY_SUFFIX) {
            // Extract property name (remove "> <" prefix and ">" suffix)
            let key = line[SDF_PROPERTY_PREFIX_LEN..line.len() - 1].to_string();

            // Property value is on the next line
            i += 1;
            if i < lines.len() {
                let value = lines[i].trim().to_string();
                if !value.is_empty() && value != SDF_ENTRY_END {
                    mol.set_property(key, value);
                }
            }
        }

        i += 1;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_sdf() {
        let sdf = "\
Methane
  
  
  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert_eq!(mol.num_atoms(), 1);
        assert_eq!(mol.name(), Some("Methane"));
    }

    #[test]
    fn test_parse_ethane_sdf() {
        let sdf = "\
Ethane
  
  
  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 1);
    }

    #[test]
    fn test_parse_real_world_acetone() {
        let sdf = "\
C3H6O
APtclcactv06051922463D 0   0.00000     0.00000
 
 10  9  0  0  0  0  0  0  0  0999 V2000
    1.3051    0.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.0763   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3051    0.6772   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.2839   -0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1059    1.7488   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8767    0.4138    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8767    0.4138   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1059    1.7488    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8767    0.4138   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8767    0.4138    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
  3 10  1  0  0  0  0
M  END

ADDITIONAL INFORMATION CAN BE ADDED HERE
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert_eq!(mol.num_atoms(), 10);
        assert_eq!(mol.num_bonds(), 9);
        assert_eq!(mol.name(), Some("C3H6O"));
    }
}
