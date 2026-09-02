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
/// Atom coordinates are preserved. SDF carries per-atom x/y/z, and what that
/// means depends on z: an all-zero z is a flat drawing and is stored as the
/// molecule's 2D layout, while any non-zero z is a conformer and is stored as
/// its 3D coordinates. Coordinates are set only if every atom supplied a
/// parseable triple — otherwise the molecule parses normally but carries
/// neither.
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
    let mut coords: Vec<Point3> = Vec::with_capacity(num_atoms);
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
            // without geometry.
            None => all_coords_parsed = false,
        }
    }

    // Set coordinates before the bond block: `add_atom` discards them (the
    // new atom has no position), while `add_bond` preserves them, so this has
    // to come after every atom is added. Only set them if every atom supplied
    // one, matching Molecule's all-or-nothing coordinate model.
    if all_coords_parsed && coords.len() == num_atoms {
        // A molfile's atom block is x/y/z whether the record is a drawing or a
        // conformer, and z is what distinguishes them. All-zero z is a 2D
        // depiction and belongs in the layout; any non-zero z is geometry and
        // belongs in the conformer.
        //
        // Filling both would claim a conformer for every flat drawing, and
        // storing a 3D record's x/y as a layout is the projection that
        // superimposes atoms differing only in depth — which is what this
        // parser did until now, silently.
        if coords.iter().all(|point| point.z == 0.0) {
            let flat: Vec<Point2> = coords.iter().map(|point| point.to_2d()).collect();
            mol.set_coords(flat)
                .map_err(|e| SdfError::ParseError(e.to_string()))?;
        } else {
            mol.set_coords3(coords)
                .map_err(|e| SdfError::ParseError(e.to_string()))?;
        }
    }

    // Parse bond block (follows atom block)
    for i in 0..num_bonds {
        let line_idx = SDF_ATOM_BLOCK_START + num_atoms + i;
        if line_idx >= lines.len() {
            return Err(SdfError::ParseError("Not enough bond lines".to_string()));
        }
        parse_bond_line(&mut mol, lines[line_idx])?;
    }

    // The `M  CHG` / `M  ISO` lines, before the data block.
    //
    // Written since #197 and, until #197, not read: the file was correct and
    // RDKit could load it, while this crate's own round trip still lost the
    // charge. `test_declared_masks_match_what_actually_survives` is what said
    // so — the mask cannot claim `FORMAL_CHARGE` while only half the trip
    // works.
    let prop_start = SDF_ATOM_BLOCK_START + num_atoms + num_bonds;
    parse_atom_property_lines(&mut mol, &lines[prop_start..]);

    // Parse optional properties block (follows bond block)
    parse_properties(&mut mol, &lines[prop_start..])?;

    // Calculate implicit hydrogens for each atom based on valence rules
    mol.calculate_implicit_hydrogens();

    // Perceive aromaticity, rather than relying on the file to assert it.
    //
    // A molfile normally states a Kekulé form — alternating single and double
    // bonds — because bond type 4 is a query type that belongs in a
    // substructure search, not a structure record. So benzene arrives as three
    // double bonds and three single ones, and a reader that takes the file
    // literally hands back a molecule that is chemically benzene but claims no
    // aromatic ring.
    //
    // That only became visible once `write_sdf` started emitting Kekulé forms
    // (#197): the round trip lost the aromatic flag, and
    // `test_declared_masks_match_what_actually_survives` refused the mask that
    // still claimed it. Perceiving on read is what every toolkit does and what
    // makes the claim true again.
    crate::io::aromaticity::detect_aromaticity(&mut mol);

    // Double-bond geometry, for the same reason and from the same place: a
    // molfile states where the atoms are, not which isomer they make, so a
    // reader that takes it literally returns a molecule that is drawn cis and
    // claims nothing. Left alone where the drawing does not say — a substituent
    // on one end only has no configuration to read (#198).
    crate::core::stereo::perceive_bond_stereo(&mut mol);
    Ok(mol)
}

/// Reads `M  CHG` and `M  ISO` lines into the atoms they name.
///
/// Both carry a count followed by that many `(atom, value)` pairs, atom indices
/// 1-based. A malformed line is skipped rather than fatal: these are optional
/// enrichments of an already-complete record, and a molfile with a mangled
/// charge line still describes a molecule.
///
/// An `M  CHG` line supersedes the atom block's old `ccc` charge column, which
/// this crate never wrote and does not read.
fn parse_atom_property_lines(mol: &mut Molecule, lines: &[&str]) {
    for line in lines {
        let tag = if line.starts_with("M  CHG") {
            true
        } else if line.starts_with("M  ISO") {
            false
        } else {
            if line.starts_with("M  END") {
                break;
            }
            continue;
        };

        let fields: Vec<i32> = line[6..]
            .split_whitespace()
            .filter_map(|f| f.parse().ok())
            .collect();
        // First field is the count; the rest are pairs.
        for pair in fields.iter().skip(1).collect::<Vec<_>>().chunks(2) {
            let [index, value] = pair else { continue };
            let Some(atom_idx) = usize::try_from(**index).ok().and_then(|i| i.checked_sub(1))
            else {
                continue;
            };
            if atom_idx >= mol.num_atoms() {
                continue;
            }
            let atom = mol.atom(atom_idx).clone();
            *mol.atom_mut(atom_idx) = if tag {
                atom.with_charge(**value as i8)
            } else {
                atom.with_isotope(**value as u16)
            };
        }
    }
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
fn parse_atom_line(mol: &mut Molecule, line: &str) -> Result<Option<Point3>, SdfError> {
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

    // Fields after the symbol are `dd ccc sss`: mass difference, old-style
    // charge, atom stereo parity. Only the parity is read — `M  CHG` and
    // `M  ISO` carry the other two and supersede them.
    //
    // The parity's convention and this crate's `Chirality` are anchored to the
    // same thing since #191's `normalise_chirality`: neighbours in increasing
    // index order, implicit hydrogen last. So it maps straight across, and
    // `write_sdf`'s `atom_parity` is the same mapping read backwards.
    let parity = parts
        .get(SDF_ATOM_FIELD + 2)
        .and_then(|field| field.parse::<u8>().ok());
    let atom = match parity {
        Some(1) => Atom::new(element).with_chirality(Chirality::Clockwise),
        Some(2) => Atom::new(element).with_chirality(Chirality::CounterClockwise),
        _ => Atom::new(element),
    };
    mol.add_atom(atom);

    // Fields 0-2 are x, y, z, all three kept. The caller decides from z
    // whether they are a layout or a conformer. A non-numeric coordinate
    // yields None (no geometry) rather than an error, so files that parsed
    // before this change still parse.
    let point = match (
        parts[0].parse::<f64>(),
        parts[1].parse::<f64>(),
        parts[2].parse::<f64>(),
    ) {
        (Ok(x), Ok(y), Ok(z)) => Some(Point3::new(x, y, z)),
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

    // The field after the type is bond stereo. Value 3 on a double bond means
    // "cis or trans, either" — the file declining to say, which is different
    // from the file not mentioning it. Recorded as `Unspecified` so
    // `perceive_bond_stereo` knows to leave the bond alone rather than read a
    // configuration out of coordinates the writer never intended as a claim.
    let declines_to_say = bond_order == BondOrder::Double
        && parts.get(3).and_then(|f| f.parse::<u8>().ok()) == Some(3);

    let bond = Bond::new(atom1, atom2, bond_order);
    let bond = if declines_to_say {
        bond.with_stereo(BondStereo::Unspecified)
    } else {
        bond
    };

    // Add bond to molecule (validates atom indices internally)
    mol.add_bond(bond)
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
/// Geometry is written from whichever coordinate set the molecule carries: a
/// conformer goes out as x/y/z with the header's dimensional code set to `3D`,
/// a layout goes out as x/y with a zero z and `2D`. A molecule with neither is
/// written with zeros, and [`molecule_has_coords_for_sdf`] lets a caller check
/// first rather than discovering it in the file.
pub fn write_sdf(mol: &Molecule) -> String {
    // A molecule that states double-bond geometry but carries no drawing is
    // laid out first, because geometry is the *only* channel V2000 has for it:
    // there is no field to write, so a record with no coordinates cannot
    // express cis at all.
    //
    // Without this the attribute would survive for a molecule that happens to
    // have been laid out and vanish for one that has not, and `Carries` is a
    // per-format claim with no way to say "sometimes". Generating the layout
    // makes the claim true for every molecule rather than most of them.
    //
    // Narrow on purpose: only when there is a configuration to lose and no
    // drawing to hold it. An ordinary molecule with no coordinates is still
    // written with zeros, as the caller expects.
    let laid_out;
    let mol = if mol.coords().is_none() && needs_geometry_for_stereo(mol) {
        laid_out = crate::core::layout::with_coords(mol);
        &laid_out
    } else {
        mol
    };

    let mut out = String::with_capacity(128 + mol.num_atoms() * 70 + mol.num_bonds() * 22);

    // Line 0: the name. Lines 1 and 2 are the program line and a comment; both
    // may be blank, and blank is more honest than inventing content.
    out.push_str(mol.name().unwrap_or(""));
    out.push('\n');

    // Line 1 is the program line, whose columns 21-22 carry the dimensional
    // code. A reader that trusts the header rather than sniffing z needs it to
    // say 3D, and that is the difference between a conformer surviving a trip
    // through another toolkit and being treated as a drawing. Blank when there
    // is nothing to declare — claiming 2D for a molecule with no coordinates
    // at all would be inventing content.
    out.push_str(match (mol.has_coords3(), mol.has_coords()) {
        (true, _) => "                    3D\n",
        (false, true) => "                    2D\n",
        (false, false) => "\n",
    });

    // Line 2 is a free-text comment; blank is more honest than inventing one.
    out.push('\n');

    // Line 3: counts. The trailing fields are the spec's defaults —
    // `0999 V2000` is the version marker every V2000 file carries.
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        mol.num_atoms(),
        mol.num_bonds()
    ));

    let coords3 = mol.coords3();
    let coords = mol.coords();
    for index in 0..mol.num_atoms() {
        // The conformer wins when there is one: it is the physical geometry,
        // and a layout alongside it is a drawing of the same molecule rather
        // than a competing set of positions. A zero z on the layout path is
        // reporting what is held, not flattening something.
        let (x, y, z) = match (coords3, coords) {
            (Some(points), _) => (points[index].x, points[index].y, points[index].z),
            (None, Some(points)) => (points[index].x, points[index].y, 0.0),
            (None, None) => (0.0, 0.0, 0.0),
        };
        let symbol = ELEMENT_SYMBOLS
            .get(mol.atom(index).atomic_number() as usize)
            .copied()
            .unwrap_or("*");
        // Columns after the symbol are `dd ccc sss ...`: mass difference,
        // old-style charge, atom stereo parity. `dd` and `ccc` stay zero — the
        // `M  CHG` and `M  ISO` lines above supersede them, and a reader seeing
        // both forms is entitled to believe either.
        let parity = atom_parity(mol, index);
        out.push_str(&format!(
            "{x:>10.4}{y:>10.4}{z:>10.4} {symbol:<3} 0  0{parity:>3}  0  0  0  0  0  0  0  0  0\n"
        ));
    }

    // A Kekulé form when one exists, so the record carries structure bonds
    // rather than the type-4 query bond. `None` means no valid assignment, and
    // then type 4 is still the honest answer — see `kekulize`.
    let kekulised = crate::io::aromaticity::kekulize(mol);
    let wedge = wedge_bonds(mol);

    for index in 0..mol.num_bonds() {
        let bond = mol.bond(index);
        let effective_order = kekulised
            .as_ref()
            .map_or_else(|| bond.order(), |orders| orders[index]);
        let order = match effective_order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            // Quadruple has no molfile bond type. Writing 1 would silently
            // change the chemistry; 8 is the spec's "any", which at least does
            // not assert something false.
            _ => 8,
        };
        // Molfile atom indices are 1-based. The field after the order is bond
        // stereo: 1 is a wedge rising towards the reader, 6 a hash falling away.
        //
        // A wedge is directional in a way the other fields are not: its narrow
        // end must sit on the stereocentre, which means the *first* atom of the
        // line. Written in storage order instead, `N[C@@H](C)C(=O)O` put the
        // narrow end on the nitrogen, and RDKit read it as a claim about the
        // nitrogen — parsed the record happily, and reported no stereocentre at
        // all.
        // A double bond that states no configuration, but whose ends both
        // carry a substituent, needs saying so explicitly. Any 2D drawing puts
        // those substituents on *some* side, so a reader perceiving geometry
        // from coordinates — ours included — would read a configuration the
        // molecule never claimed. Field 3 is "either", which is a statement of
        // ignorance rather than of geometry, and is what RDKit writes here.
        if bond.order() == BondOrder::Double
            && bond.stereo() == BondStereo::None
            && geometry_is_inferrable(mol, bond.atom1(), bond.atom2())
        {
            out.push_str(&format!(
                "{:>3}{:>3}{order:>3}  3  0  0  0\n",
                bond.atom1() + 1,
                bond.atom2() + 1
            ));
            continue;
        }

        let (stereo, (first, second)) = match wedge.get(&index) {
            Some(&(centre, direction)) if centre == bond.atom2() => {
                (direction, (bond.atom2(), bond.atom1()))
            }
            Some(&(_, direction)) => (direction, (bond.atom1(), bond.atom2())),
            None => (0, (bond.atom1(), bond.atom2())),
        };
        out.push_str(&format!(
            "{:>3}{:>3}{order:>3}{stereo:>3}  0  0  0\n",
            first + 1,
            second + 1
        ));
    }

    out.push_str(&property_lines(mol));
    out.push_str("M  END\n");
    out.push_str(&data_block(mol));
    out.push_str(SDF_ENTRY_END);
    out.push('\n');
    out
}

/// Whether a reader would derive a configuration for this double bond from a
/// drawing of it — that is, whether both ends carry something to be on a side.
fn geometry_is_inferrable(mol: &Molecule, left: usize, right: usize) -> bool {
    let has_other =
        |atom: usize, partner: usize| mol.neighbors(atom).iter().any(|n| n.atom_idx != partner);
    has_other(left, right) && has_other(right, left)
}

/// Whether this molecule states a double-bond configuration that only a
/// drawing can carry.
fn needs_geometry_for_stereo(mol: &Molecule) -> bool {
    (0..mol.num_bonds()).any(|i| {
        let bond = mol.bond(i);
        bond.order() == BondOrder::Double && matches!(bond.stereo(), BondStereo::E | BondStereo::Z)
    })
}

/// `M  CHG` and `M  ISO` lines, between the bond block and `M  END`.
///
/// These are what make a charged record *readable*, not merely lossless.
/// Without `M  CHG` a quaternary nitrogen has four bonds and no charge, which
/// is a valence error, so a strict reader refuses the whole record rather than
/// quietly returning a neutral amine — `C[N+](C)(C)C` and `[O-][N+](=O)c1ccccc1`
/// both did exactly that (#194).
///
/// Both lines carry a count in three columns followed by that many
/// `(atom, value)` pairs in four columns each, atom indices 1-based as in the
/// bond block. **Eight pairs per line at most**, so a molecule with more needs
/// several — a limit that is easy to miss because nothing in a small corpus
/// reaches it.
fn property_lines(mol: &Molecule) -> String {
    fn emit(tag: &str, entries: &[(usize, i32)]) -> String {
        let mut out = String::new();
        for chunk in entries.chunks(8) {
            out.push_str(tag);
            out.push_str(&format!("{:>3}", chunk.len()));
            for (index, value) in chunk {
                out.push_str(&format!("{:>4}{value:>4}", index + 1));
            }
            out.push('\n');
        }
        out
    }

    let charges: Vec<(usize, i32)> = (0..mol.num_atoms())
        .filter_map(|i| match mol.atom(i).formal_charge() {
            0 => None,
            c => Some((i, i32::from(c))),
        })
        .collect();

    // The absolute mass number, not the difference from natural abundance. The
    // atom line's `dd` column is the difference form and is superseded; it stays
    // zero, because a reader that sees both is entitled to believe either.
    let isotopes: Vec<(usize, i32)> = (0..mol.num_atoms())
        .filter_map(|i| mol.atom(i).isotope().map(|iso| (i, i32::from(iso))))
        .collect();

    // Nothing to say means no line. `M  CHG  0` is legal and is noise.
    format!("{}{}", emit("M  CHG", &charges), emit("M  ISO", &isotopes))
}

/// The data block: everything after `M  END` and before `$$$$`.
///
/// Data fields are why people choose SDF over SMILES — assay values, catalogue
/// IDs, docking scores. `parse_sdf` has always read them into
/// [`Molecule::properties`]; the writer never wrote them back, so an SDF to SDF
/// round trip through this crate silently lost every one (#188).
///
/// Sorted by key, because `properties` is a `HashMap` and iterating it writes
/// the same molecule differently between runs. #153 was that bug for ring
/// closures, and the fix is the same one.
fn data_block(mol: &Molecule) -> String {
    let mut keys: Vec<&String> = mol.properties().keys().collect();
    keys.sort();

    let mut out = String::new();
    for key in keys {
        // The blank line is the field terminator, not decoration: a reader
        // takes every line up to it as the value.
        out.push_str(&format!(
            "> <{key}>\n{}\n\n",
            mol.properties()[key.as_str()]
        ));
    }
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

/// Whether writing this molecule would record real coordinates, of either
/// kind.
///
/// Exposed so a caller can warn before writing rather than after: a file full
/// of zeros is a plausible-looking result that is entirely useless, and the
/// difference is invisible in the output.
pub fn molecule_has_coords_for_sdf(mol: &Molecule) -> bool {
    mol.has_coords() || mol.has_coords3()
}

/// The atom stereo parity for the atom line's `sss` column, or 0 for none.
///
/// Both conventions are anchored to the same thing since `normalise_chirality`:
/// neighbours in increasing index order, with an implicit hydrogen last. So
/// this is a direct mapping rather than a permutation problem — the permutation
/// was done once, at parse time, where the written order was still known.
///
/// Which way round the mapping goes is fixed by RDKit rather than by reading
/// the specification: both alanine enantiomers must keep their distinct `/m0`
/// and `/m1` InChI layers through an SDF trip, and a centre carrying a
/// ring-closure digit must survive too.
fn atom_parity(mol: &Molecule, atom_idx: usize) -> u8 {
    let base = match mol.atom(atom_idx).chirality() {
        Chirality::Clockwise => 1,
        Chirality::CounterClockwise => 2,
        Chirality::None | Chirality::Unspecified => return 0,
    };

    // A parity describes four things around a centre. Three neighbours and one
    // implicit hydrogen is the common shape; anything else is not a tetrahedral
    // centre this can describe.
    let neighbours = mol.neighbors(atom_idx).len();
    let hydrogens = usize::from(mol.atom(atom_idx).total_hydrogens());
    if neighbours + hydrogens != 4 {
        return 0;
    }
    base
}

/// Which bonds carry a wedge, in which direction, and around which centre.
///
/// A wedge is not another way of writing the parity — it is a claim about the
/// *drawing*, so it only means anything when there is one. `core/layout.rs`
/// places atoms with no knowledge of chirality, so deriving a wedge straight
/// from [`Chirality`] would produce a marker that contradicts the picture it
/// sits on, and a reader is entitled to believe either.
///
/// So the direction is not derived by formula. It is found by **simulating the
/// reader**: build the 3D arrangement each candidate direction implies, work
/// out the parity a reader would compute from it, and keep the one that agrees
/// with [`atom_parity`]. A formula was tried first and got alanine right and
/// *trans*-cyclohexanediol wrong — one centre of two — which is what a guess
/// that happens to fit one arrangement looks like.
///
/// When neither direction reproduces the parity, or there is no layout, the
/// parity stands alone: under-specified and readable, rather than confidently
/// contradictory.
fn wedge_bonds(mol: &Molecule) -> std::collections::HashMap<usize, (usize, u8)> {
    let mut wedges = std::collections::HashMap::new();
    let Some(points) = mol.coords() else {
        return wedges; // no drawing, so nothing a wedge could mean
    };

    for atom_idx in 0..mol.num_atoms() {
        let parity = atom_parity(mol, atom_idx);
        if parity == 0 {
            continue;
        }

        let neighbours = mol.neighbors(atom_idx);
        // Prefer a bond to a low-degree neighbour, so the wedge reads as coming
        // out of the centre rather than out of a ring.
        let Some(chosen) = neighbours
            .iter()
            .min_by_key(|n| mol.graph().degree(n.atom_idx))
        else {
            continue;
        };

        for direction in [1u8, 6u8] {
            if implied_parity(mol, points, atom_idx, chosen.atom_idx, direction) == Some(parity) {
                wedges.insert(chosen.bond_idx, (atom_idx, direction));
                break;
            }
        }
    }
    wedges
}

/// The parity a reader would compute from the drawing, if `wedged` were lifted
/// out of the page (`1`) or pushed behind it (`6`).
///
/// Everything not wedged lies in the plane. An implicit hydrogen has no drawn
/// position, so it is placed where a real one would go: opposite the sum of the
/// other bonds.
fn implied_parity(
    mol: &Molecule,
    points: &[Point2],
    centre: usize,
    wedged: usize,
    direction: u8,
) -> Option<u8> {
    let origin = points[centre];
    let mut vectors: Vec<(usize, [f64; 3])> = Vec::new();
    for neighbour in mol.neighbors(centre) {
        let p = points[neighbour.atom_idx];
        let z = if neighbour.atom_idx == wedged {
            if direction == 1 { 1.0 } else { -1.0 }
        } else {
            0.0
        };
        vectors.push((neighbour.atom_idx, [p.x - origin.x, p.y - origin.y, z]));
    }

    if mol.atom(centre).total_hydrogens() == 1 {
        // Opposite everything else, and numbered highest — molfile treats an
        // implicit hydrogen as the last neighbour.
        let sum = vectors.iter().fold([0.0; 3], |acc, (_, v)| {
            [acc[0] + v[0], acc[1] + v[1], acc[2] + v[2]]
        });
        vectors.push((usize::MAX, [-sum[0], -sum[1], -sum[2]]));
    }
    if vectors.len() != 4 {
        return None;
    }

    // Molfile numbering order.
    vectors.sort_by_key(|(idx, _)| *idx);

    // With the highest-numbered neighbour pointing away, do the first three
    // read clockwise? That is the sign of the determinant of the three, taken
    // relative to the fourth.
    let d = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
    };
    let det = d(vectors[0].1, vectors[1].1, vectors[2].1);
    if det.abs() < 1e-9 {
        return None; // the drawing is degenerate here
    }
    Some(if det < 0.0 { 1 } else { 2 })
}

#[cfg(test)]
mod tests {

    /// The `M  CHG` / `M  ISO` lines of a written record.
    fn property_lines_of(smiles: &str) -> Vec<String> {
        let mol = crate::io::smiles::parse_smiles(smiles).expect(smiles);
        write_sdf(&mol)
            .lines()
            .filter(|l| l.starts_with("M  CHG") || l.starts_with("M  ISO"))
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn test_charges_are_written_as_a_property_line() {
        assert_eq!(property_lines_of("[NH4+]"), ["M  CHG  1   1   1"]);
        assert_eq!(property_lines_of("[Cl-]"), ["M  CHG  1   1  -1"]);
        assert_eq!(property_lines_of("[Mg+2]"), ["M  CHG  1   1   2"]);
    }

    #[test]
    fn test_a_neutral_molecule_writes_no_property_line() {
        // `M  CHG  0` is legal and is noise.
        assert!(property_lines_of("CCO").is_empty());
    }

    #[test]
    fn test_more_than_eight_charges_wrap_onto_a_second_line() {
        // The limit is 8 pairs per line. Nothing in the corpus reaches it, so
        // without this the wrap would ship untested and break on the first real
        // polyelectrolyte.
        let nine = "[NH4+].".repeat(9);
        let lines = property_lines_of(nine.trim_end_matches('.'));
        assert_eq!(lines.len(), 2, "{lines:?}");
        assert!(lines[0].starts_with("M  CHG  8"), "{}", lines[0]);
        assert!(lines[1].starts_with("M  CHG  1"), "{}", lines[1]);
    }

    #[test]
    fn test_isotopes_are_written_as_the_absolute_mass_number() {
        // 13, not 1: `M  ISO` supersedes the atom block's difference-from-
        // natural-abundance column, which stays zero.
        assert_eq!(property_lines_of("[13CH4]"), ["M  ISO  1   1  13"]);
    }

    #[test]
    fn test_charge_and_isotope_survive_a_round_trip() {
        for smiles in ["[NH4+]", "[Cl-]", "[Mg+2]", "[13CH4]"] {
            let before = crate::io::smiles::parse_smiles(smiles).expect(smiles);
            let after = parse_sdf(&write_sdf(&before)).expect(smiles);
            assert_eq!(
                after.atom(0).formal_charge(),
                before.atom(0).formal_charge(),
                "{smiles} lost its charge"
            );
            assert_eq!(
                after.atom(0).isotope(),
                before.atom(0).isotope(),
                "{smiles} lost its isotope"
            );
        }
    }

    #[test]
    fn test_data_fields_are_written_sorted_and_terminated() {
        let mut mol = crate::io::smiles::parse_smiles("CCO").expect("valid SMILES");
        mol.set_property("zeta".into(), "last".into());
        mol.set_property("alpha".into(), "first".into());

        let text = write_sdf(&mol);
        let block = text.split("M  END\n").nth(1).expect("a data block");
        // Sorted, because `properties` is a HashMap and iteration order would
        // otherwise write the same molecule differently between runs (#153).
        assert_eq!(block, "> <alpha>\nfirst\n\n> <zeta>\nlast\n\n$$$$\n");

        let back = parse_sdf(&text).expect("round trips");
        assert_eq!(back.property("alpha"), Some("first"));
        assert_eq!(back.property("zeta"), Some("last"));
    }

    #[test]
    fn test_aromatic_bonds_are_written_as_a_kekule_form() {
        // Type 4 is a *query* bond type. Benzene survives it because a reader
        // can kekulise unambiguously; pyrrole cannot, and RDKit rejects the
        // record (#194).
        for smiles in [
            "c1ccccc1",
            "c1cc[nH]c1",
            "c1ccncc1",
            "c1ccoc1",
            "c1cnc[nH]1",
        ] {
            let mol = crate::io::smiles::parse_smiles(smiles).expect(smiles);
            let text = write_sdf(&mol);
            let type_four = text
                .lines()
                .skip(4 + mol.num_atoms())
                .take(mol.num_bonds())
                .filter(|l| l.split_whitespace().nth(2) == Some("4"))
                .count();
            assert_eq!(type_four, 0, "{smiles} still writes a query bond type");
        }
    }

    #[test]
    fn test_a_chiral_centre_writes_a_parity_and_reads_it_back() {
        let before = crate::io::smiles::parse_smiles("N[C@@H](C)C(=O)O").expect("valid SMILES");
        let after = parse_sdf(&write_sdf(&before)).expect("round trips");
        assert_ne!(after.atom(1).chirality(), Chirality::None);
        assert_eq!(after.atom(1).chirality(), before.atom(1).chirality());
    }

    #[test]
    fn test_the_two_enantiomers_do_not_collapse_into_each_other() {
        // The failure a self-consistent writer hides: emit one parity
        // regardless of input and every round trip still looks perfect.
        let l = crate::io::smiles::parse_smiles("N[C@@H](C)C(=O)O").expect("valid");
        let d = crate::io::smiles::parse_smiles("N[C@H](C)C(=O)O").expect("valid");
        let l_back = parse_sdf(&write_sdf(&l)).expect("round trips");
        let d_back = parse_sdf(&write_sdf(&d)).expect("round trips");
        assert_ne!(l_back.atom(1).chirality(), d_back.atom(1).chirality());
    }

    #[test]
    fn test_no_layout_means_a_parity_and_no_wedge() {
        // A wedge is a claim about a drawing. With no drawing it would be a
        // claim about nothing, and a reader that trusts it over the parity gets
        // a different molecule.
        let mol = crate::io::smiles::parse_smiles("N[C@@H](C)C(=O)O").expect("valid SMILES");
        assert!(mol.coords().is_none());
        let text = write_sdf(&mol);
        let wedges = text
            .lines()
            .skip(4 + mol.num_atoms())
            .take(mol.num_bonds())
            .filter(|l| matches!(l.split_whitespace().nth(3), Some("1") | Some("6")))
            .count();
        assert_eq!(wedges, 0, "wedge written with no layout to justify it");
    }
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

        // A 3D file, so the coordinates are a conformer and not a layout.
        // Until #174 this stored x/y as a 2D layout, which projected atoms 6
        // and 7 — identical but for z — onto the same point. The molecule is
        // left with no layout deliberately: `ensure_coords` can compute a
        // readable one, and a flattened conformer is not that.
        assert!(!mol.has_coords(), "a conformer is not a depiction");
        let conformer = mol.coords3().expect("conformer");
        assert_eq!(conformer.len(), mol.num_atoms());

        assert_ne!(mol.coord3(5), mol.coord3(6));
        assert_eq!(mol.coord3(5).unwrap().z, 0.8900);
        assert_eq!(mol.coord3(6).unwrap().z, -0.8900);
    }

    #[test]
    fn test_flat_sdf_is_a_layout_not_a_conformer() {
        // Every z is zero, so this is a drawing. Storing it as a conformer
        // would claim geometry the file does not assert, and would put a
        // planar "structure" into any 3D format written from it.
        let sdf = "\
Flat


  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$";

        let mol = parse_sdf(sdf).unwrap();
        assert!(mol.has_coords());
        assert!(!mol.has_coords3());
        assert_eq!(mol.coord(1), Some(Point2::new(1.5, 0.5)));
    }

    #[test]
    fn test_conformer_survives_a_write_read_round_trip() {
        // The point of the whole change: z has to come back. Before #174 the
        // writer emitted a hardcoded 0.0 in the z column, so this molecule
        // came back flat.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::oxygen()));
        mol.add_bond(Bond::new(0, 1, BondOrder::Single)).unwrap();
        mol.set_coords3(vec![
            Point3::new(0.1234, -0.5678, 1.2345),
            Point3::new(1.5000, 0.5000, -0.7500),
        ])
        .unwrap();

        let text = write_sdf(&mol);
        let back = parse_sdf(&text).unwrap();

        let conformer = back.coords3().expect("conformer should survive");
        assert_eq!(conformer[0], Point3::new(0.1234, -0.5678, 1.2345));
        assert_eq!(conformer[1], Point3::new(1.5, 0.5, -0.75));
        assert!(!back.has_coords(), "a conformer must not become a layout");
    }

    #[test]
    fn test_writer_declares_the_dimensional_code() {
        // Columns 21-22 of the program line. A reader that trusts the header
        // instead of sniffing z gets the right answer only if we write it.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));

        mol.set_coords3(vec![Point3::new(0.0, 0.0, 1.0)]).unwrap();
        assert_eq!(
            write_sdf(&mol).lines().nth(1).unwrap().trim_end(),
            "                    3D"
        );

        let mut flat = Molecule::new();
        flat.add_atom(Atom::new(Element::carbon()));
        flat.set_coords(vec![Point2::new(0.0, 0.0)]).unwrap();
        assert_eq!(
            write_sdf(&flat).lines().nth(1).unwrap().trim_end(),
            "                    2D"
        );

        let mut bare = Molecule::new();
        bare.add_atom(Atom::new(Element::carbon()));
        assert_eq!(write_sdf(&bare).lines().nth(1).unwrap(), "");
    }

    #[test]
    fn test_has_coords_for_sdf_accepts_either_kind() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        assert!(!molecule_has_coords_for_sdf(&mol));

        mol.set_coords3(vec![Point3::new(0.0, 0.0, 1.0)]).unwrap();
        assert!(molecule_has_coords_for_sdf(&mol));
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
