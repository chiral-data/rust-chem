use crate::core::prelude::*;
use crate::io::errors::SdfError;

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
///
/// Atom coordinates are preserved: SDF carries per-atom x/y/z, and the x/y are
/// stored on the molecule for depiction. Coordinates are set only if every atom
/// supplied a parseable pair — otherwise the molecule parses normally but
/// without a layout.
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
    let mut coords: Vec<Point2> = Vec::with_capacity(num_atoms);
    let mut all_coords_parsed = true;
    for i in 0..num_atoms {
        let line_idx = SDF_ATOM_BLOCK_START + i;
        if line_idx >= lines.len() {
            return Err(SdfError::ParseError("Not enough atom lines".to_string()));
        }
        match parse_atom_line(&mut mol, lines[line_idx])? {
            Some(point) => coords.push(point),
            // An unparseable coordinate isn't fatal — the atom itself is
            // still valid, so the molecule parses as it always did, just
            // without a layout.
            None => all_coords_parsed = false,
        }
    }

    // Set coordinates before the bond block: `add_atom` discards them (the
    // new atom has no position), while `add_bond` preserves them, so this has
    // to come after every atom is added. Only set them if every atom supplied
    // one, matching Molecule's all-or-nothing coordinate model.
    if all_coords_parsed && coords.len() == num_atoms {
        mol.set_coords(coords)
            .map_err(|e| SdfError::ParseError(e.to_string()))?;
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
fn parse_atom_line(mol: &mut Molecule, line: &str) -> Result<Option<Point2>, SdfError> {
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

    // Fields 0-2 are x, y, z. Only x/y are kept — coordinate storage is 2D,
    // and for a 3D SDF this is a straight projection down the z axis, which is
    // a serviceable starting depiction but can superimpose atoms that are only
    // separated in z. A non-numeric coordinate yields None (no layout) rather
    // than an error, so files that parsed before this change still parse.
    let point = match (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
        (Ok(x), Ok(y)) => Some(Point2::new(x, y)),
        _ => None,
    };

    Ok(point)
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
            )));
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

// ============================================================================
// WRITING
// ============================================================================

/// Writes one molecule as a molfile V2000 record, terminated by `$$$$`.
///
/// This is the only format here that can carry coordinates, which is why it
/// exists: 2D layout is computed by `crate::core::layout`, and SMILES has nowhere
/// to put the result. Writing a laid-out molecule as SMILES silently discards
/// the geometry.
///
/// # Fixed columns
///
/// The molfile spec fixes the column of every field. [`parse_sdf`] is
/// deliberately lenient and splits on whitespace instead, so it would accept
/// looser output than this produces — but other software will not, and a file
/// only this crate can read is not an interchange format. So the widths here
/// follow the spec (`%10.4f` coordinates, a left-justified 3-column symbol,
/// `%3d` bond fields) rather than what the local parser happens to tolerate.
///
/// # What is not written
///
/// Charges, isotopes, chirality and bond stereo. The atom block emits the
/// zero-valued fields the spec expects in their place, so the record is
/// structurally complete and a reader that wants them finds defaults rather
/// than absent columns. Aromatic bonds are written as type 4, which round-trips
/// through [`parse_sdf`] but is a query bond type in strict molfile — a Kekulé
/// form would be more portable and needs a Kekulisation pass that does not
/// exist yet.
///
/// A molecule with no coordinates is written with zeros, and
/// [`molecule_has_coords_for_sdf`] lets a caller check first rather than
/// discovering it in the file.
pub fn write_sdf(mol: &Molecule) -> String {
    let mut out = String::with_capacity(128 + mol.num_atoms() * 70 + mol.num_bonds() * 22);

    // Line 0: the name. Lines 1 and 2 are the program line and a comment; both
    // may be blank, and blank is more honest than inventing content.
    out.push_str(mol.name().unwrap_or(""));
    out.push('\n');
    out.push('\n');
    out.push('\n');

    // Line 3: counts. The trailing fields are the spec's defaults —
    // `0999 V2000` is the version marker every V2000 file carries.
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        mol.num_atoms(),
        mol.num_bonds()
    ));

    let coords = mol.coords();
    for index in 0..mol.num_atoms() {
        // z is always 0: coordinate storage is 2D. A molecule read from a 3D
        // SDF already lost its z on the way in, so writing 0 is reporting what
        // is held rather than flattening something.
        let (x, y) = match coords {
            Some(points) => (points[index].x, points[index].y),
            None => (0.0, 0.0),
        };
        let symbol = ELEMENT_SYMBOLS
            .get(mol.atom(index).atomic_number() as usize)
            .copied()
            .unwrap_or("*");
        out.push_str(&format!(
            "{x:>10.4}{y:>10.4}{:>10.4} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            0.0
        ));
    }

    for index in 0..mol.num_bonds() {
        let bond = mol.bond(index);
        let order = match bond.order() {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            // Quadruple has no molfile bond type. Writing 1 would silently
            // change the chemistry; 8 is the spec's "any", which at least does
            // not assert something false.
            _ => 8,
        };
        // Molfile atom indices are 1-based.
        out.push_str(&format!(
            "{:>3}{:>3}{order:>3}  0  0  0  0\n",
            bond.atom1() + 1,
            bond.atom2() + 1
        ));
    }

    out.push_str("M  END\n");
    out.push_str(SDF_ENTRY_END);
    out.push('\n');
    out
}

/// Writes several molecules as one SDF file.
pub fn write_sdf_all<'a>(molecules: impl IntoIterator<Item = &'a Molecule>) -> String {
    let mut out = String::new();
    for mol in molecules {
        out.push_str(&write_sdf(mol));
    }
    out
}

/// Whether writing this molecule would record real coordinates.
///
/// Exposed so a caller can warn before writing rather than after: a file full
/// of zeros is a plausible-looking result that is entirely useless, and the
/// difference is invisible in the output.
pub fn molecule_has_coords_for_sdf(mol: &Molecule) -> bool {
    mol.has_coords()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- writing ----------------------------------------------------------
    //
    // Every test here goes through `parse_sdf`, because the writer's job is to
    // produce something a reader accepts. Asserting on expected text instead
    // would pass for a file nothing can read, and would need rewriting every
    // time a column moved.

    use crate::core::layout::ensure_coords;
    use crate::io::smiles::parse_smiles;

    fn laid_out(smiles: &str) -> Molecule {
        let mut mol = parse_smiles(smiles).expect("valid SMILES");
        assert!(ensure_coords(&mut mol), "layout should succeed");
        mol
    }

    #[test]
    fn test_a_molecule_survives_a_round_trip() {
        let original = laid_out("CCO");
        let back = parse_sdf(&write_sdf(&original)).expect("our own output should parse");

        assert_eq!(back.num_atoms(), original.num_atoms());
        assert_eq!(back.num_bonds(), original.num_bonds());
        for i in 0..original.num_atoms() {
            assert_eq!(
                back.atom(i).atomic_number(),
                original.atom(i).atomic_number(),
                "atom {i}"
            );
        }
    }

    #[test]
    fn test_coordinates_survive_a_round_trip() {
        // The whole reason this writer exists. Tolerance is the four decimal
        // places the format stores, not floating-point epsilon.
        let original = laid_out("c1ccccc1O");
        let back = parse_sdf(&write_sdf(&original)).expect("parse");

        let before = original.coords().expect("laid out");
        let after = back.coords().expect("coordinates should survive");
        assert_eq!(before.len(), after.len());
        for (i, (a, b)) in before.iter().zip(after).enumerate() {
            assert!(
                (a.x - b.x).abs() < 1e-4 && (a.y - b.y).abs() < 1e-4,
                "atom {i}: ({}, {}) became ({}, {})",
                a.x,
                a.y,
                b.x,
                b.y
            );
        }
    }

    #[test]
    fn test_bond_orders_survive_a_round_trip() {
        // Single, double, triple and aromatic, in one molecule each, since a
        // writer that mapped every order to 1 would still round-trip atoms
        // and coordinates perfectly.
        for smiles in ["CC", "C=C", "C#C", "c1ccccc1", "CC(=O)O"] {
            let original = laid_out(smiles);
            let back = parse_sdf(&write_sdf(&original)).expect(smiles);
            let mut before: Vec<_> = (0..original.num_bonds())
                .map(|i| original.bond(i).order())
                .collect();
            let mut after: Vec<_> = (0..back.num_bonds())
                .map(|i| back.bond(i).order())
                .collect();
            before.sort_by_key(|o| format!("{o:?}"));
            after.sort_by_key(|o| format!("{o:?}"));
            assert_eq!(before, after, "{smiles}");
        }
    }

    #[test]
    fn test_bond_topology_survives_not_just_the_count() {
        // Off-by-one on the 1-based index conversion would keep the bond count
        // and rewire the molecule.
        let original = laid_out("CC(C)C");
        let back = parse_sdf(&write_sdf(&original)).expect("parse");

        let pairs = |m: &Molecule| {
            let mut v: Vec<(usize, usize)> = (0..m.num_bonds())
                .map(|i| {
                    let b = m.bond(i);
                    let (x, y) = (b.atom1(), b.atom2());
                    if x <= y { (x, y) } else { (y, x) }
                })
                .collect();
            v.sort();
            v
        };
        assert_eq!(pairs(&back), pairs(&original));
    }

    #[test]
    fn test_the_name_survives() {
        let mut mol = laid_out("CCO");
        mol.set_name("ethanol".to_string());
        let back = parse_sdf(&write_sdf(&mol)).expect("parse");
        assert_eq!(back.name(), Some("ethanol"));
    }

    #[test]
    fn test_a_molecule_with_no_layout_writes_zeros_and_says_so() {
        // Not an error: a caller may legitimately want connectivity only. But
        // a file of zeros looks like a real result, so the caller needs to be
        // able to ask first.
        let mol = parse_smiles("CCO").expect("valid SMILES");
        assert!(!molecule_has_coords_for_sdf(&mol));
        let text = write_sdf(&mol);
        assert!(text.contains("0.0000    0.0000"));
        assert_eq!(parse_sdf(&text).expect("still parses").num_atoms(), 3);
    }

    #[test]
    fn test_a_single_atom_round_trips() {
        let mol = laid_out("C");
        let back = parse_sdf(&write_sdf(&mol)).expect("parse");
        assert_eq!(back.num_atoms(), 1);
        assert_eq!(back.num_bonds(), 0);
    }

    #[test]
    fn test_writing_many_molecules_produces_records_the_reader_splits() {
        // `write_sdf_all` has to agree with how `crate::io::reader` splits records,
        // or a multi-molecule file reads back as one molecule or none.
        let molecules: Vec<Molecule> = ["CCO", "c1ccccc1", "CC(=O)O"]
            .iter()
            .map(|s| laid_out(s))
            .collect();
        let text = write_sdf_all(&molecules);
        let outcome = crate::io::reader::read_sdf(&text);
        assert!(outcome.skipped.is_empty(), "{:?}", outcome.skipped);
        assert_eq!(outcome.len(), 3);
        for (read, original) in outcome.records.iter().zip(&molecules) {
            assert_eq!(read.molecule.num_atoms(), original.num_atoms());
            assert!(read.molecule.has_coords());
        }
    }

    #[test]
    fn test_the_counts_line_matches_the_blocks() {
        // The counts line is what the reader trusts to find the bond block, so
        // a mismatch shifts everything after it.
        let mol = laid_out("c1ccc2ccccc2c1");
        let text = write_sdf(&mol);
        let counts: Vec<usize> = text
            .lines()
            .nth(3)
            .expect("counts line")
            .split_whitespace()
            .take(2)
            .map(|f| f.parse().expect("numeric"))
            .collect();
        assert_eq!(counts, vec![mol.num_atoms(), mol.num_bonds()]);
    }

    #[test]
    fn test_the_record_is_terminated() {
        let text = write_sdf(&laid_out("CCO"));
        assert!(text.contains("M  END"));
        assert!(text.trim_end().ends_with("$$$$"));
    }

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

    #[test]
    fn test_coordinates_are_preserved() {
        let sdf = "\
Ethane


  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert!(mol.has_coords());
        assert_eq!(mol.coord(0), Some(Point2::new(0.0, 0.0)));
        assert_eq!(mol.coord(1), Some(Point2::new(1.5, 0.5)));
    }

    #[test]
    fn test_negative_and_fractional_coordinates() {
        // Real files routinely carry negative and multi-decimal coordinates.
        let sdf = "\
Acetone
APtclcactv06051922463D 0   0.00000     0.00000

  3  2  0  0  0  0  0  0  0  0999 V2000
    1.3051    0.6772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.0763   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3051    0.6772   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        let coords = mol.coords().expect("coordinates");
        assert_eq!(coords.len(), 3);
        assert_eq!(coords[0], Point2::new(1.3051, 0.6772));
        assert_eq!(coords[1], Point2::new(0.0, -0.0763));
        assert_eq!(coords[2], Point2::new(-1.3051, 0.6772));
    }

    #[test]
    fn test_coordinates_are_one_per_atom() {
        let sdf = "\
Acetone
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
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert_eq!(mol.coords().expect("coordinates").len(), mol.num_atoms());
        // A 3D file: z differs between atoms 6 and 7, but only x/y is stored,
        // so they project onto the same 2D point.
        assert_eq!(mol.coord(5), mol.coord(6));
    }

    #[test]
    fn test_unparseable_coordinate_yields_no_layout() {
        // "abc" in the x field: the atom is still valid (the symbol parses),
        // so the molecule parses as before -- just without coordinates,
        // rather than failing outright.
        let sdf = "\
Broken


  2  1  0  0  0  0  0  0  0  0999 V2000
      abc    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 1);
        assert!(!mol.has_coords());
    }
}
