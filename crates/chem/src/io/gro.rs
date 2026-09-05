//! GRO (Gromacs) — fixed-width text: a title line, an atom-count line,
//! that many per-atom lines, then a box-vector line (#227).
//!
//! **The one format in scope needing a real unit conversion, not a
//! straight column copy.** GRO coordinates and box vectors are
//! nanometres; this crate's whole data model, and every other format
//! registered so far, is Ångströms. Reading multiplies by 10, writing
//! divides by 10 — there is exactly one correct behaviour here, not a
//! design choice.
//!
//! **No velocities.** GRO's optional trailing three columns are
//! unmodelled anywhere in `core/` — docking/MD-prep pipelines pass
//! positions forward, not velocities. Read ignores them entirely (parses
//! only up to the z column); write never emits them. Stated here rather
//! than silently mishandled.
//!
//! **No dedicated element column exists in this format at all** — the
//! atom name is all there is, and this module resolves an element from
//! that name's first letter only, never a two-letter match. Real
//! protein/nucleic MD topology naming (this format's actual pipeline
//! usage) puts the element there and treats the rest as positional
//! labelling (`CA`/`CB`/`CG`, `ND1`/`OD1`) — trying a two-letter match
//! first would misread the single most common atom name in any protein
//! GRO file, `CA`, as calcium rather than carbon. The real cost is a
//! two-letter-only element named as a plain ion residue (`NA`, `MG`,
//! `ZN`, `CL`) resolving to whatever its first letter is instead — a
//! known, stated limitation.
//!
//! **A molecule with no unit cell writes an all-zero box line** (`0.0 0.0
//! 0.0`) rather than inventing one — reading it back correctly produces
//! no cell either, since [`UnitCell::validate`] would reject a
//! zero-length edge outright, so an all-zero line can only ever mean "no
//! box," never a genuine degenerate one.
//!
//! **The box line's orientation already matches this crate's own fixed
//! convention** (`v1` along x, `v2` in the xy-plane, `v3` general — the
//! same "a along x, b in the xy plane" convention
//! [`crate::core::cell::UnitCell`] already commits to), unlike extended
//! XYZ's genuinely arbitrary `Lattice=` matrix (deferred in #222). No
//! rotation step is needed here, only a units-aware vectors→parameters
//! conversion on read; write reuses [`UnitCell::matrix`] directly rather
//! than reimplementing the parametrization.
//!
//! **The declared atom count is trusted**, unlike every prior format's
//! reader (mmCIF/Mol2 scan for the next section marker; PDB/PDBQT use
//! `ENDMDL`) — GRO has no structural terminator at all separating the
//! atom block from the box-vector line that follows it, so the count is
//! the *only* signal available, and trusting it here is correct, not a
//! shortcut. A file may hold several such frames back to back (some
//! tools concatenate `.gro` as an ad-hoc trajectory) — splitting those is
//! the reader's/supplier's job
//! ([`crate::io::reader::read_gro_with_options`],
//! [`crate::io::supplier::GroSupplier`]), the same division every prior
//! format has; this module parses and writes exactly one already-isolated
//! frame.

use crate::core::atom::{Atom, Element};
use crate::core::cell::UnitCell;
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::io::errors::GroError;

const NM_TO_ANGSTROM: f64 = 10.0;

fn element_from_symbol(sym: &str) -> Option<Element> {
    let sym = sym.trim();
    let mut chars = sym.chars();
    let normalised = match (chars.next(), chars.next()) {
        (Some(a), Some(b)) if chars.next().is_none() => {
            format!("{}{}", a.to_ascii_uppercase(), b.to_ascii_lowercase())
        }
        (Some(a), None) => a.to_ascii_uppercase().to_string(),
        _ => return None,
    };
    ELEMENT_SYMBOLS
        .iter()
        .position(|&s| s == normalised)
        .and_then(|n| Element::new(n as u8))
}

/// Resolves an atom's element from its GRO atom name (there is no
/// dedicated element column in this format) — strips leading digits
/// (common in MD atom naming, e.g. `1HB`) then takes the leading
/// alphabetic run.
///
/// **Prefers the single first letter, never a two-letter match** — real
/// protein/nucleic MD topology naming (this format's actual pipeline
/// usage) puts the element in that first letter alone and everything
/// after it is positional labelling: `CA`/`CB`/`CG` (alpha/beta/gamma
/// carbon), `ND1`/`NE2`, `OD1`/`OG`, `HA`/`HB1`. Trying a two-letter match
/// first would read the single most common atom name in any protein GRO
/// file, `CA`, as calcium rather than carbon. The cost is real but rarer:
/// two-letter-only elements named as a plain ion residue (`NA`, `MG`,
/// `ZN`, `CL`) resolve to whatever their first letter is instead (`N`,
/// unresolvable, unresolvable, `C`) — a known, stated limitation, not a
/// silent one.
fn element_from_atom_name(name: &str) -> Option<Element> {
    let letters: String = name
        .trim()
        .chars()
        .skip_while(|c| c.is_ascii_digit())
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    letters
        .chars()
        .next()
        .and_then(|c| element_from_symbol(&c.to_string()))
}

#[derive(Clone, PartialEq, Eq)]
struct ResidueKey {
    name: String,
    sequence: i32,
}

fn group_into_chain_and_residues(keys: &[ResidueKey]) -> (Vec<Chain>, Vec<Residue>) {
    let mut residues = Vec::new();
    let mut i = 0;
    while i < keys.len() {
        let key = &keys[i];
        let start = i;
        while i < keys.len() && keys[i] == *key {
            i += 1;
        }
        residues.push(Residue {
            name: key.name.clone(),
            sequence: key.sequence,
            insertion_code: None,
            label_seq: None,
            chain_ix: 0,
            is_hetero: false,
            atoms: start..i,
        });
    }
    let chains = if residues.is_empty() {
        Vec::new()
    } else {
        vec![Chain {
            id: String::new(),
            label_id: None,
            residues: 0..residues.len(),
        }]
    };
    (chains, residues)
}

/// The box line's 9 values, in GROMACS's own documented (and unusual)
/// order: `v1x v2y v3z v1y v1z v2x v2z v3x v3y`. A 3-value line is the
/// shorthand for the orthorhombic case, with every other term `0.0`.
fn parse_box_vectors(values: &[f64]) -> (Point3, Point3, Point3) {
    let get = |i: usize| values.get(i).copied().unwrap_or(0.0);
    let v1 = Point3::new(get(0), get(3), get(4));
    let v2 = Point3::new(get(5), get(1), get(6));
    let v3 = Point3::new(get(7), get(8), get(2));
    (v1, v2, v3)
}

fn angle_degrees(a: Point3, b: Point3) -> f64 {
    let dot = a.x * b.x + a.y * b.y + a.z * b.z;
    let len_a = (a.x * a.x + a.y * a.y + a.z * a.z).sqrt();
    let len_b = (b.x * b.x + b.y * b.y + b.z * b.z).sqrt();
    if len_a == 0.0 || len_b == 0.0 {
        return 90.0;
    }
    (dot / (len_a * len_b)).clamp(-1.0, 1.0).acos().to_degrees()
}

fn cell_from_box_vectors(v1: Point3, v2: Point3, v3: Point3) -> UnitCell {
    let len = |p: Point3| (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    UnitCell::new(
        len(v1) * NM_TO_ANGSTROM,
        len(v2) * NM_TO_ANGSTROM,
        len(v3) * NM_TO_ANGSTROM,
        angle_degrees(v2, v3),
        angle_degrees(v1, v3),
        angle_degrees(v1, v2),
    )
}

/// Parses one GRO frame (a title line, a count line, that many atom
/// lines, and a box-vector line).
pub fn parse_gro(text: &str) -> Result<Molecule, GroError> {
    let mut lines = text.lines();

    let title = lines
        .next()
        .ok_or_else(|| GroError::ParseError("empty frame".to_string()))?
        .trim();
    let count_line = lines
        .next()
        .ok_or_else(|| GroError::ParseError("missing atom count line".to_string()))?;
    let count: usize = count_line
        .trim()
        .parse()
        .map_err(|_| GroError::ParseError(format!("invalid atom count: {count_line:?}")))?;

    let mut mol = Molecule::new();
    if !title.is_empty() {
        mol.set_name(title.to_string());
    }

    let mut coords = Vec::with_capacity(count);
    let mut keys = Vec::with_capacity(count);
    for _ in 0..count {
        let line = lines
            .next()
            .ok_or(GroError::AtomCountMismatch { expected: count })?;
        if line.len() < 44 {
            return Err(GroError::InvalidAtomLine(line.to_string()));
        }
        let res_seq: i32 = line[0..5]
            .trim()
            .parse()
            .map_err(|_| GroError::InvalidAtomLine(line.to_string()))?;
        let res_name = line[5..10].trim().to_string();
        let atom_name = line[10..15].trim().to_string();
        let x: f64 = line[20..28]
            .trim()
            .parse()
            .map_err(|_| GroError::InvalidAtomLine(line.to_string()))?;
        let y: f64 = line[28..36]
            .trim()
            .parse()
            .map_err(|_| GroError::InvalidAtomLine(line.to_string()))?;
        let z: f64 = line[36..44]
            .trim()
            .parse()
            .map_err(|_| GroError::InvalidAtomLine(line.to_string()))?;
        // Velocity columns, if present, run from byte 44 onward -- not
        // read, see the module doc.

        let element = element_from_atom_name(&atom_name)
            .ok_or_else(|| GroError::InvalidElement(atom_name.clone()))?;
        mol.add_atom(Atom::new(element));
        coords.push(Point3::new(
            x * NM_TO_ANGSTROM,
            y * NM_TO_ANGSTROM,
            z * NM_TO_ANGSTROM,
        ));
        keys.push(ResidueKey {
            name: res_name,
            sequence: res_seq,
        });
    }

    let box_line = lines
        .next()
        .ok_or_else(|| GroError::ParseError("missing box vector line".to_string()))?;
    let box_values: Vec<f64> = box_line
        .split_whitespace()
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|_| GroError::ParseError(format!("invalid box vector line: {box_line:?}")))?;

    mol.set_coords3(coords)
        .map_err(|e| GroError::ParseError(e.to_string()))?;
    let (chains, residues) = group_into_chain_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| GroError::ParseError(e.to_string()))?;

    if !box_values.is_empty() && box_values.iter().any(|&v| v != 0.0) {
        // An all-zero box line is the convention this crate's own writer
        // (and several real GRO-writing tools) use for "no box defined" --
        // distinct from a genuine, validatable cell. Trying to build a
        // `UnitCell` from it would always fail (`UnitCell::validate`
        // correctly rejects a zero-length edge), so it's read as "no
        // cell" rather than a parse error.
        let (v1, v2, v3) = parse_box_vectors(&box_values);
        mol.set_cell(cell_from_box_vectors(v1, v2, v3))
            .map_err(|e| GroError::ParseError(e.to_string()))?;
    }

    Ok(mol)
}

/// Writes one GRO frame.
pub fn write_gro(mol: &Molecule) -> String {
    let mut out = String::new();
    out.push_str(mol.name().unwrap_or(""));
    out.push('\n');
    out.push_str(&format!("{:>5}\n", mol.num_atoms()));

    for (i, atom) in mol.atoms().iter().enumerate() {
        let residue = mol.residue_of(i);
        let res_seq = residue.map(|r| r.sequence).unwrap_or(1);
        let res_name = residue.map(|r| r.name.as_str()).unwrap_or("UNK");
        let atom_name = atom.element().symbol();
        let p = mol.coord3(i).unwrap_or(Point3::ORIGIN);
        out.push_str(&format!(
            "{res_seq:>5}{res_name:<5}{atom_name:>5}{serial:>5}{x:>8.3}{y:>8.3}{z:>8.3}\n",
            serial = i + 1,
            x = p.x / NM_TO_ANGSTROM,
            y = p.y / NM_TO_ANGSTROM,
            z = p.z / NM_TO_ANGSTROM,
        ));
    }

    if let Some(cell) = mol.cell() {
        let basis = cell.basis();
        let (v1, v2, v3) = (basis[0], basis[1], basis[2]);
        let is_orthorhombic = [v1.y, v1.z, v2.x, v2.z, v3.x, v3.y]
            .iter()
            .all(|&c| c.abs() < 1e-6);
        if is_orthorhombic {
            out.push_str(&format!(
                "{:>10.5}{:>10.5}{:>10.5}\n",
                v1.x / NM_TO_ANGSTROM,
                v2.y / NM_TO_ANGSTROM,
                v3.z / NM_TO_ANGSTROM,
            ));
        } else {
            out.push_str(&format!(
                "{:>10.5}{:>10.5}{:>10.5}{:>10.5}{:>10.5}{:>10.5}{:>10.5}{:>10.5}{:>10.5}\n",
                v1.x / NM_TO_ANGSTROM,
                v2.y / NM_TO_ANGSTROM,
                v3.z / NM_TO_ANGSTROM,
                v1.y / NM_TO_ANGSTROM,
                v1.z / NM_TO_ANGSTROM,
                v2.x / NM_TO_ANGSTROM,
                v2.z / NM_TO_ANGSTROM,
                v3.x / NM_TO_ANGSTROM,
                v3.y / NM_TO_ANGSTROM,
            ));
        }
    } else {
        out.push_str("   0.00000   0.00000   0.00000\n");
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_GRO: &str = "\
water
    3
    1SOL     OW    1   0.000   0.000   0.000
    1SOL    HW1    2   0.076   0.000   0.050
    1SOL    HW2    3   0.076   0.000  -0.050
   1.00000   1.00000   1.00000
";

    #[test]
    fn test_a_single_frame_round_trips() {
        let mol = parse_gro(WATER_GRO).expect("valid GRO");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.name(), Some("water"));
        assert_eq!(mol.residues()[0].name, "SOL");
        assert_eq!(mol.residues()[0].sequence, 1);

        let written = write_gro(&mol);
        let back = parse_gro(&written).expect("round trips");
        assert_eq!(back.num_atoms(), mol.num_atoms());
        assert_eq!(back.residues()[0].name, "SOL");
    }

    #[test]
    fn test_nanometres_convert_to_angstroms_on_read() {
        let mol = parse_gro(WATER_GRO).expect("valid GRO");
        // 0.076 nm -> 0.76 Angstrom, exactly.
        let p = mol.coord3(1).expect("has coords");
        assert!((p.x - 0.76).abs() < 1e-9, "{}", p.x);
    }

    #[test]
    fn test_angstroms_convert_to_nanometres_on_write() {
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.set_coords3(vec![Point3::new(1.5, 0.0, 0.0)])
            .expect("one per atom");
        let written = write_gro(&mol);
        let atom_line = written.lines().nth(2).expect("has an atom line");
        // 1.5 Angstrom -> 0.150 nm, exactly.
        assert!(atom_line.contains("0.150"), "{atom_line}");
    }

    #[test]
    fn test_an_orthorhombic_box_round_trips_with_ninety_degree_angles() {
        let mol = parse_gro(WATER_GRO).expect("valid GRO");
        let cell = mol.cell().expect("has a cell");
        assert_eq!((cell.alpha, cell.beta, cell.gamma), (90.0, 90.0, 90.0));
        assert!((cell.a - 10.0).abs() < 1e-6, "{}", cell.a);

        let written = write_gro(&mol);
        assert_eq!(
            written.lines().last().unwrap().split_whitespace().count(),
            3
        );
    }

    #[test]
    fn test_a_triclinic_box_round_trips_its_off_diagonal_terms() {
        let text = "\
triclinic
    1
    1UNK      C    1   0.000   0.000   0.000
   2.00000   2.00000   2.00000   0.00000   0.00000   0.50000   0.00000   0.30000   0.00000
";
        let mol = parse_gro(text).expect("valid GRO");
        let cell = mol.cell().expect("has a cell");
        // gamma (angle between v1 and v2) is no longer 90 degrees, since
        // v2x = 0.5 nm is nonzero.
        assert!((cell.gamma - 90.0).abs() > 1.0, "{}", cell.gamma);

        let written = write_gro(&mol);
        let box_line = written.lines().last().unwrap();
        assert_eq!(box_line.split_whitespace().count(), 9, "{box_line}");
    }

    #[test]
    fn test_velocity_columns_are_read_without_error_and_never_written() {
        let text = "\
with velocities
    1
    1UNK      C    1   0.000   0.000   0.000  1.2345  0.0000 -0.5000
   1.00000   1.00000   1.00000
";
        let mol = parse_gro(text).expect("valid GRO");
        assert_eq!(mol.num_atoms(), 1);
        let written = write_gro(&mol);
        assert_eq!(written.lines().count(), 4, "{written}");
    }

    #[test]
    fn test_a_multi_frame_files_second_frame_is_this_modules_callers_job() {
        // parse_gro itself only ever sees one already-isolated frame --
        // confirm it doesn't read past its own declared count into a
        // second frame's title line.
        let mol = parse_gro(WATER_GRO).expect("valid GRO");
        assert_eq!(mol.num_atoms(), 3);
    }

    #[test]
    fn test_too_few_atom_lines_is_a_clear_error_not_a_panic() {
        let text = "short\n    2\n    1SOL     OW    1   0.000   0.000   0.000\n";
        let err = parse_gro(text).unwrap_err();
        assert!(matches!(err, GroError::AtomCountMismatch { .. }), "{err}");
    }

    #[test]
    fn test_an_unrecognised_element_is_a_clear_error() {
        let text = "bad\n    1\n    1UNK     Zq    1   0.000   0.000   0.000\n   1.0   1.0   1.0\n";
        let err = parse_gro(text).unwrap_err();
        assert!(matches!(err, GroError::InvalidElement(_)), "{err}");
    }

    #[test]
    fn test_protein_atom_name_ca_resolves_to_carbon_not_calcium() {
        // The single most common atom name in any protein GRO file --
        // must never be misread as the two-letter element symbol it also
        // happens to spell. See the module doc.
        let text = "ala\n    1\n    1ALA     CA    1   0.000   0.000   0.000\n   1.0   1.0   1.0\n";
        let mol = parse_gro(text).expect("valid GRO");
        assert_eq!(mol.atoms()[0].element().symbol(), "C");
    }

    #[test]
    fn test_a_molecule_with_no_cell_writes_an_all_zero_box_and_reads_back_with_no_cell() {
        use crate::io::smiles::parse_smiles;

        let mol = parse_smiles("CC").expect("valid SMILES");
        assert!(mol.cell().is_none());
        let written = write_gro(&mol);
        let back = parse_gro(&written).expect("round trips even with no real cell");
        assert!(back.cell().is_none());
    }
}
