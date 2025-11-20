// FIX: Cannot work with Realworld SDF files...

use chemcore::prelude::*;

use crate::errors::SdfError;

/// Parses an SDF (Structure-Data File) string into a `Molecule` object.
///
/// # SDF Format Overview
/// SDF is a chemical file format used to represent molecular structures and associated data.
/// The format consists of:
/// - Header block (lines 0-3): Name, program info, comments, and counts line
/// - Atom block: One line per atom with coordinates and element symbol
/// - Bond block: One line per bond with atom indices and bond order
/// - Properties block: Optional key-value pairs for molecular properties
///
/// # Parameters
/// - `sdf`: A string slice containing the SDF-formatted molecular data
///
/// # Returns
/// - `Ok(Molecule)`: Successfully parsed molecule with atoms, bonds, and properties
/// - `Err(SdfError)`: Parsing error with details about what went wrong
///
/// # Errors
/// Returns `SdfError` if:
/// - The SDF has fewer than 4 lines (minimum required for header)
/// - The counts line cannot be parsed
/// - There are insufficient atom or bond lines for the declared counts
/// - Atom or bond lines are malformed
/// - Unknown element symbols are encountered
///
/// # SDF Structure Example
/// ```text
/// Methane                      <- Line 0: Molecule name
///                               <- Line 1: Program/user info (optional)
///                               <- Line 2: Comment (optional)
///   1  0  0  0  0...           <- Line 3: Counts line (atoms, bonds, ...)
///     0.0000    0.0000... C    <- Line 4+: Atom lines (coordinates + symbol)
/// M  END                        <- Marks end of structure
/// > <PropertyName>              <- Optional properties
/// PropertyValue
/// $$$$                          <- Marks end of entry
/// ```
pub fn parse_sdf(sdf: &str) -> Result<Molecule, SdfError> {
    let lines: Vec<&str> = sdf.lines().collect();

    // Validate minimum SDF structure (header must be at least 4 lines)
    if lines.len() < 4 {
        return Err(SdfError::ParseError("Too few lines".to_string()));
    }

    let mut mol = Molecule::new();

    // Line 0: Molecule name (optional)
    if !lines[0].trim().is_empty() {
        mol.set_name(lines[0].trim().to_string());
    }

    // Lines 1-2: Program information and comments (skipped)

    // Line 3: Counts line - specifies number of atoms and bonds
    let counts_line = lines[3];
    let (num_atoms, num_bonds) = parse_counts_line(counts_line)?;

    // Parse atom block (starts at line 4)
    // Each atom line contains: coordinates (x, y, z) and element symbol
    for i in 0..num_atoms {
        let line_idx = 4 + i;
        if line_idx >= lines.len() {
            return Err(SdfError::ParseError("Not enough atom lines".to_string()));
        }
        parse_atom_line(&mut mol, lines[line_idx])?;
    }

    // Parse bond block (follows atom block)
    // Each bond line contains: atom1_idx, atom2_idx, bond_order
    for i in 0..num_bonds {
        let line_idx = 4 + num_atoms + i;
        if line_idx >= lines.len() {
            return Err(SdfError::ParseError("Not enough bond lines".to_string()));
        }
        parse_bond_line(&mut mol, lines[line_idx])?;
    }

    // Parse optional properties block (follows bond block)
    // Properties are key-value pairs in format: "> <KEY>" followed by value
    let prop_start = 4 + num_atoms + num_bonds;
    parse_properties(&mut mol, &lines[prop_start..])?;

    // Calculate implicit hydrogens for each atom based on valence rules
    mol.calculate_implicit_hydrogens();
    Ok(mol)
}

/// Parses the counts line (line 3) of an SDF file.
///
/// # SDF Counts Line Format
/// The counts line contains space-separated integers:
/// ```text
/// aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
/// ```
/// Where:
/// - `aaa`: Number of atoms (positions 0-2)
/// - `bbb`: Number of bonds (positions 3-5)
/// - Other fields are mostly obsolete or optional
///
/// # Parameters
/// - `line`: The counts line string (line 3 of SDF)
///
/// # Returns
/// - `Ok((num_atoms, num_bonds))`: Tuple containing atom and bond counts
/// - `Err(SdfError::MissingCountsLine)`: If parsing fails
///
/// # Example
/// ```text
/// Input:  "  2  1  0  0  0  0  0  0  0  0999 V2000"
/// Output: Ok((2, 1))  // 2 atoms, 1 bond
/// ```
fn parse_counts_line(line: &str) -> Result<(usize, usize), SdfError> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    // Counts line must have at least 2 fields (atom count and bond count)
    if parts.len() < 2 {
        return Err(SdfError::MissingCountsLine);
    }

    let num_atoms = parts[0].parse().map_err(|_| SdfError::MissingCountsLine)?;
    let num_bonds = parts[1].parse().map_err(|_| SdfError::MissingCountsLine)?;

    Ok((num_atoms, num_bonds))
}

/// Parses a single atom line from the SDF atom block.
///
/// # SDF Atom Line Format (V2000)
/// Fixed-width format (minimum 34 characters):
/// ```text
/// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
/// |--------10---------|----10---| |3|
///     x coord         y    z     symbol
/// ```
/// Positions:
/// - 0-9: X coordinate (10 characters, format: xxxxx.xxxx)
/// - 10-19: Y coordinate (10 characters)
/// - 20-29: Z coordinate (10 characters)
/// - 31-33: Atom symbol (3 characters, right-aligned, e.g., " C", " N", "Cl")
///
/// # Parameters
/// - `mol`: Mutable reference to the molecule being constructed
/// - `line`: The atom line string to parse
///
/// # Returns
/// - `Ok(())`: Atom successfully parsed and added to molecule
/// - `Err(SdfError::InvalidAtomLine)`: If line is malformed or contains unknown element
///
/// # Notes
/// - Coordinates are currently parsed but not stored (TODO: add coordinate support)
/// - Element symbol mapping is hardcoded (FIXME: use periodic table lookup)
///
/// # Example
/// ```text
/// Input: "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
///                                        ^^^
///                                      Carbon atom
/// ```
fn parse_atom_line(mol: &mut Molecule, line: &str) -> Result<(), SdfError> {
    // Validate minimum line length for V2000 format
    if line.len() < 34 {
        return Err(SdfError::InvalidAtomLine(line.to_string()));
    }

    // Extract element symbol from fixed positions 31-33 (right-aligned, trimmed)
    let symbol = line[31..34].trim();

    // FIXME: Replace hardcoded mapping with periodic table lookup
    // Map element symbol to atomic number
    let atomic_number = match symbol {
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
        _ => {
            return Err(SdfError::InvalidAtomLine(format!(
                "Unknown element: {}",
                symbol
            )))
        }
    };

    // Create element from atomic number
    let element =
        Element::new(atomic_number).ok_or_else(|| SdfError::InvalidAtomLine(symbol.to_string()))?;

    // Create and add atom to molecule
    let atom = Atom::new(element);
    mol.add_atom(atom);

    Ok(())
}

/// Parses a single bond line from the SDF bond block.
///
/// # SDF Bond Line Format (V2000)
/// Space-separated integers:
/// ```text
/// 111222tttsssxxxrrrccc
/// |3||3||3|
///  |  |  └─ Bond type (1=Single, 2=Double, 3=Triple, 4=Aromatic)
///  |  └──── Second atom index (1-based)
///  └─────── First atom index (1-based)
/// ```
///
/// # Parameters
/// - `mol`: Mutable reference to the molecule being constructed
/// - `line`: The bond line string to parse
///
/// # Returns
/// - `Ok(())`: Bond successfully parsed and added to molecule
/// - `Err(SdfError::InvalidBondLine)`: If line is malformed or indices are invalid
///
/// # Bond Type Mapping
/// - 1 → Single bond
/// - 2 → Double bond
/// - 3 → Triple bond
/// - 4 → Aromatic bond
/// - Other values → Error
///
/// # Important Notes
/// - SDF uses 1-based atom indexing; this function converts to 0-based
/// - Atom indices are validated when adding the bond to the molecule
///
/// # Example
/// ```text
/// Input: "  1  2  1  0  0  0  0"
///          ^  ^  ^
///          |  |  └─ Bond type 1 (single bond)
///          |  └──── Atom 2
///          └─────── Atom 1
/// Result: Single bond between atoms 0 and 1 (converted to 0-based)
/// ```
fn parse_bond_line(mol: &mut Molecule, line: &str) -> Result<(), SdfError> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    // Bond line must have at least 3 fields: atom1, atom2, bond_type
    if parts.len() < 3 {
        return Err(SdfError::InvalidBondLine(line.to_string()));
    }

    // Parse atom indices (convert from 1-based SDF to 0-based internal indexing)
    let atom1: usize = parts[0]
        .parse::<usize>()
        .map_err(|_| SdfError::InvalidBondLine(line.to_string()))?
        - 1; // SDF uses 1-based indexing

    let atom2: usize = parts[1]
        .parse::<usize>()
        .map_err(|_| SdfError::InvalidBondLine(line.to_string()))?
        - 1;

    // Parse bond type
    let bond_type: u8 = parts[2]
        .parse()
        .map_err(|_| SdfError::InvalidBondLine(line.to_string()))?;

    // Map SDF bond type to internal BondOrder enum
    let bond_order = match bond_type {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => {
            return Err(SdfError::InvalidBondLine(format!(
                "Invalid bond type: {}",
                bond_type
            )))
        }
    };

    // Add bond to molecule (validates atom indices internally)
    mol.add_bond(Bond::new(atom1, atom2, bond_order))
        .map_err(|e| SdfError::ParseError(e.to_string()))?;

    Ok(())
}

/// Parses the optional properties block at the end of an SDF entry.
///
/// # SDF Properties Format
/// Properties appear after the atom and bond blocks, using the format:
/// ```text
/// > <PropertyName>
/// PropertyValue
/// > <AnotherProperty>
/// AnotherValue
/// M  END
/// $$$$
/// ```
///
/// # Parameters
/// - `mol`: Mutable reference to the molecule to store properties
/// - `lines`: Slice of remaining lines after the bond block
///
/// # Returns
/// - `Ok(())`: Properties successfully parsed and stored
/// - `Err(SdfError)`: If property parsing encounters an error (currently unused)
///
/// # Property Format
/// - Property names: Enclosed in `> < >`, e.g., `> <MOLECULAR_WEIGHT>`
/// - Property values: On the next line after the property name
/// - Termination: `$$$$` marks the end of the SDF entry
///
/// # Behavior
/// - Stops parsing when encountering `$$$$` (entry terminator)
/// - Stores each property as a key-value pair in the molecule
/// - Empty or malformed properties are silently skipped
///
/// # Example
/// ```text
/// Input lines:
/// > <MOLECULAR_WEIGHT>
/// 16.04
/// > <FORMULA>
/// CH4
/// $$$$
///
/// Result: Molecule has properties:
/// - "MOLECULAR_WEIGHT" = "16.04"
/// - "FORMULA" = "CH4"
/// ```
fn parse_properties(mol: &mut Molecule, lines: &[&str]) -> Result<(), SdfError> {
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i].trim();

        // Property names are enclosed in "> <" and ">"
        if line.starts_with("> <") && line.ends_with('>') {
            // Extract property name (remove "> <" prefix and ">" suffix)
            let key = line[3..line.len() - 1].to_string();

            // Property value is on the next line
            i += 1;
            if i < lines.len() {
                let value = lines[i].trim().to_string();
                mol.set_property(key, value);
            }
        }

        // "$$$$" marks the end of this SDF entry
        if line == "$$$$" {
            break;
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
}
