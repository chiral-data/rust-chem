//! Bonds for standard residues, from a built-in table of templates (#403).
//!
//! PDB and mmCIF leave polymer connectivity unstated: a reader is expected to
//! know that `CA` bonds to `CB` in every alanine. This table covers the 20
//! amino acids, HIS's protonation names (`HID`, `HIE`, `HIP`), `MSE`, the
//! terminal `OXT`, and the RNA/DNA nucleotides (`A C G U`, `DA DC DG DT`),
//! with heavy-atom bonds and the standard PDB v3 hydrogen names.
//!
//! Orders are Kekulé, never [`BondOrder::Aromatic`], so aromaticity stays the
//! writers' own perception. Formal charges are not set: every residue reads as
//! its neutral form.
//!
//! The function only adds bonds between atoms not yet bonded and never
//! changes an existing order, so it composes with whatever the file stated.

use crate::core::bond::{Bond, BondOrder};
use crate::core::molecule::Molecule;
use crate::core::residue::Residue;

/// Peptide-bond and phosphodiester partners are bonded only when closer than
/// this, so a gap in a model stays a gap.
pub const MAX_LINK_DISTANCE: f64 = 2.0;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Polymer {
    Peptide,
    NucleicAcid,
}

const PEPTIDE: &str = "N-CA CA-C C=O C-OXT OXT-HXT N-H N-H1 N-H2 N-H3 CA-HA";

const PHOSPHATE: &str = "OP3-P P=OP1 P-OP2 P-O5' OP2-HOP2 OP3-HOP3 O5'-HO5'";

const RIBOSE: &str = "O5'-C5' C5'-C4' C4'-O4' C4'-C3' C3'-O3' C3'-C2' C2'-C1' C1'-O4' \
    C2'-O2' C5'-H5' C5'-H5'' C4'-H4' C3'-H3' C2'-H2' O2'-HO2' C1'-H1' O3'-HO3'";

const DEOXYRIBOSE: &str = "O5'-C5' C5'-C4' C4'-O4' C4'-C3' C3'-O3' C3'-C2' C2'-C1' C1'-O4' \
    C5'-H5' C5'-H5'' C4'-H4' C3'-H3' C2'-H2' C2'-H2'' C1'-H1' O3'-HO3'";

const ADENINE: &str = "C1'-N9 N9-C8 N9-C4 C8=N7 N7-C5 C5-C6 C5=C4 C6-N6 C6=N1 N1-C2 C2=N3 \
    N3-C4 C8-H8 N6-H61 N6-H62 C2-H2";

const GUANINE: &str = "C1'-N9 N9-C8 N9-C4 C8=N7 N7-C5 C5-C6 C5=C4 C6=O6 C6-N1 N1-C2 C2-N2 \
    C2=N3 N3-C4 C8-H8 N1-H1 N2-H21 N2-H22";

const CYTOSINE: &str =
    "C1'-N1 N1-C2 N1-C6 C2=O2 C2-N3 N3=C4 C4-N4 C4-C5 C5=C6 N4-H41 N4-H42 C5-H5 C6-H6";

const URACIL: &str = "C1'-N1 N1-C2 N1-C6 C2=O2 C2-N3 N3-C4 C4=O4 C4-C5 C5=C6 N3-H3 C5-H5 C6-H6";

const THYMINE: &str = "C1'-N1 N1-C2 N1-C6 C2=O2 C2-N3 N3-C4 C4=O4 C4-C5 C5=C6 C5-C7 N3-H3 \
    C6-H6 C7-H71 C7-H72 C7-H73";

/// HIS with the proton on ND1.
const HIS_DELTA: &str = "CA-CB CB-CG CG-ND1 ND1-CE1 CE1=NE2 NE2-CD2 CD2=CG CB-HB2 CB-HB3 \
    ND1-HD1 CD2-HD2 CE1-HE1 NE2-HE2";

/// HIS with the proton on NE2.
const HIS_EPSILON: &str = "CA-CB CB-CG CG-ND1 ND1=CE1 CE1-NE2 NE2-CD2 CD2=CG CB-HB2 CB-HB3 \
    ND1-HD1 CD2-HD2 CE1-HE1 NE2-HE2";

const PHE_RING: &str = "CA-CB CB-CG CG=CD1 CG-CD2 CD1-CE1 CD2=CE2 CE1=CZ CE2-CZ CB-HB2 CB-HB3 \
    CD1-HD1 CD2-HD2 CE1-HE1 CE2-HE2";

/// The side chain (or base and sugar) of a residue, beyond the shared part.
fn side_chain(name: &str) -> Option<(Polymer, [&'static str; 3])> {
    let aa = |side: &'static str| Some((Polymer::Peptide, [PEPTIDE, side, ""]));
    let rna = |base: &'static str| Some((Polymer::NucleicAcid, [PHOSPHATE, RIBOSE, base]));
    let dna = |base: &'static str| Some((Polymer::NucleicAcid, [PHOSPHATE, DEOXYRIBOSE, base]));
    match name {
        "ALA" => aa("CA-CB CB-HB1 CB-HB2 CB-HB3"),
        "ARG" => aa(
            "CA-CB CB-CG CG-CD CD-NE NE-CZ CZ-NH1 CZ=NH2 CB-HB2 CB-HB3 CG-HG2 \
            CG-HG3 CD-HD2 CD-HD3 NE-HE NH1-HH11 NH1-HH12 NH2-HH21 NH2-HH22",
        ),
        "ASN" => aa("CA-CB CB-CG CG=OD1 CG-ND2 CB-HB2 CB-HB3 ND2-HD21 ND2-HD22"),
        "ASP" => aa("CA-CB CB-CG CG=OD1 CG-OD2 CB-HB2 CB-HB3 OD2-HD2"),
        "CYS" => aa("CA-CB CB-SG CB-HB2 CB-HB3 SG-HG"),
        "GLN" => aa(
            "CA-CB CB-CG CG-CD CD=OE1 CD-NE2 CB-HB2 CB-HB3 CG-HG2 CG-HG3 NE2-HE21 \
            NE2-HE22",
        ),
        "GLU" => aa("CA-CB CB-CG CG-CD CD=OE1 CD-OE2 CB-HB2 CB-HB3 CG-HG2 CG-HG3 OE2-HE2"),
        "GLY" => aa("CA-HA2 CA-HA3"),
        "HID" => aa(HIS_DELTA),
        "HIS" | "HIE" | "HIP" => aa(HIS_EPSILON),
        "ILE" => aa(
            "CA-CB CB-CG1 CB-CG2 CG1-CD1 CB-HB CG1-HG12 CG1-HG13 CG2-HG21 CG2-HG22 \
            CG2-HG23 CD1-HD11 CD1-HD12 CD1-HD13",
        ),
        "LEU" => aa(
            "CA-CB CB-CG CG-CD1 CG-CD2 CB-HB2 CB-HB3 CG-HG CD1-HD11 CD1-HD12 CD1-HD13 \
            CD2-HD21 CD2-HD22 CD2-HD23",
        ),
        "LYS" => aa(
            "CA-CB CB-CG CG-CD CD-CE CE-NZ CB-HB2 CB-HB3 CG-HG2 CG-HG3 CD-HD2 CD-HD3 \
            CE-HE2 CE-HE3 NZ-HZ1 NZ-HZ2 NZ-HZ3",
        ),
        "MET" => aa("CA-CB CB-CG CG-SD SD-CE CB-HB2 CB-HB3 CG-HG2 CG-HG3 CE-HE1 CE-HE2 CE-HE3"),
        "MSE" => aa("CA-CB CB-CG CG-SE SE-CE CB-HB2 CB-HB3 CG-HG2 CG-HG3 CE-HE1 CE-HE2 CE-HE3"),
        "PHE" => Some((Polymer::Peptide, [PEPTIDE, PHE_RING, "CZ-HZ"])),
        "PRO" => aa("CA-CB CB-CG CG-CD CD-N CB-HB2 CB-HB3 CG-HG2 CG-HG3 CD-HD2 CD-HD3"),
        "SER" => aa("CA-CB CB-OG CB-HB2 CB-HB3 OG-HG"),
        "THR" => aa("CA-CB CB-OG1 CB-CG2 CB-HB OG1-HG1 CG2-HG21 CG2-HG22 CG2-HG23"),
        "TRP" => aa(
            "CA-CB CB-CG CG=CD1 CG-CD2 CD1-NE1 NE1-CE2 CD2=CE2 CD2-CE3 CE2-CZ2 \
            CE3=CZ3 CZ2=CH2 CZ3-CH2 CB-HB2 CB-HB3 CD1-HD1 NE1-HE1 CE3-HE3 CZ2-HZ2 CZ3-HZ3 \
            CH2-HH2",
        ),
        "TYR" => Some((Polymer::Peptide, [PEPTIDE, PHE_RING, "CZ-OH OH-HH"])),
        "VAL" => aa(
            "CA-CB CB-CG1 CB-CG2 CB-HB CG1-HG11 CG1-HG12 CG1-HG13 CG2-HG21 CG2-HG22 \
            CG2-HG23",
        ),
        "A" => rna(ADENINE),
        "G" => rna(GUANINE),
        "C" => rna(CYTOSINE),
        "U" => rna(URACIL),
        "DA" => dna(ADENINE),
        "DG" => dna(GUANINE),
        "DC" => dna(CYTOSINE),
        "DT" => dna(THYMINE),
        _ => None,
    }
}

/// Whether `name` has a template.
pub fn has_residue_template(name: &str) -> bool {
    side_chain(name).is_some()
}

/// The atoms of `residue` named `name`, one per alternate location.
fn atoms_named(mol: &Molecule, residue: &Residue, name: &str) -> Vec<usize> {
    let Some(sites) = mol.sites() else {
        return Vec::new();
    };
    residue
        .atoms
        .clone()
        .filter(|&i| sites[i].name.as_deref() == Some(name))
        .collect()
}

/// Two atoms belong to one conformer when their `alt_loc`s match or either
/// is blank; otherwise an A/B split would be cross-bonded.
fn same_conformer(mol: &Molecule, a: usize, b: usize) -> bool {
    let alt = |i: usize| mol.site(i).and_then(|s| s.alt_loc);
    match (alt(a), alt(b)) {
        (Some(x), Some(y)) => x == y,
        _ => true,
    }
}

/// Bonds every same-conformer pair across `left` × `right` not yet bonded,
/// leaving an existing bond's order alone. Shared with the readers' links
/// (`SSBOND`/`LINK`, `_struct_conn`).
pub(crate) fn bond_conformer_pairs(
    mol: &mut Molecule,
    left: &[usize],
    right: &[usize],
    order: BondOrder,
) {
    for &a in left {
        for &b in right {
            if a != b && same_conformer(mol, a, b) && !mol.graph().has_edge(a, b) {
                // Both indices come from `mol`, so this cannot fail.
                let _ = mol.add_bond(Bond::new(a, b, order));
            }
        }
    }
}

/// The parts of a template that apply to `residue`, resolving HIS's
/// tautomer from which ring hydrogen the file carries.
fn template_for(mol: &Molecule, residue: &Residue) -> Option<(Polymer, [&'static str; 3])> {
    let (polymer, mut parts) = side_chain(&residue.name)?;
    if residue.name == "HIS"
        && !atoms_named(mol, residue, "HD1").is_empty()
        && atoms_named(mol, residue, "HE2").is_empty()
    {
        parts[1] = HIS_DELTA;
    }
    Some((polymer, parts))
}

fn parse_bond(spec: &str) -> Option<(&str, &str, BondOrder)> {
    if let Some((a, b)) = spec.split_once('=') {
        Some((a, b, BondOrder::Double))
    } else {
        spec.split_once('-').map(|(a, b)| (a, b, BondOrder::Single))
    }
}

/// Adds the bonds of every standard residue's template, then the peptide
/// C–N and phosphodiester O3'–P bonds between consecutive residues of one
/// chain that are within [`MAX_LINK_DISTANCE`].
///
/// Atoms are matched by the names in [`Molecule::sites`]; a residue not in
/// the table, or an atom with an unexpected name, gets no template bond.
pub fn add_residue_template_bonds(mol: &mut Molecule) {
    if !mol.has_sites() {
        return;
    }

    let residues = mol.residues().to_vec();
    let templates: Vec<_> = residues.iter().map(|r| template_for(mol, r)).collect();

    for (residue, template) in residues.iter().zip(&templates) {
        let Some((_, parts)) = template else {
            continue;
        };
        for spec in parts.iter().flat_map(|part| part.split_whitespace()) {
            let Some((a, b, order)) = parse_bond(spec) else {
                continue;
            };
            let left = atoms_named(mol, residue, a);
            let right = atoms_named(mol, residue, b);
            bond_conformer_pairs(mol, &left, &right, order);
        }
    }

    // Without coordinates there is no distance to gate on, so no link.
    if !mol.has_coords3() {
        return;
    }
    for chain in mol.chains().to_vec() {
        for ix in chain.residues.start..chain.residues.end.saturating_sub(1) {
            let (Some((first, _)), Some((second, _))) = (templates[ix], templates[ix + 1]) else {
                continue;
            };
            if first != second {
                continue;
            }
            let (from, to) = match first {
                Polymer::Peptide => ("C", "N"),
                Polymer::NucleicAcid => ("O3'", "P"),
            };
            let left = atoms_named(mol, &residues[ix], from);
            let right = atoms_named(mol, &residues[ix + 1], to);
            for &a in &left {
                for &b in &right {
                    if within_link_distance(mol, a, b) {
                        bond_conformer_pairs(mol, &[a], &[b], BondOrder::Single);
                    }
                }
            }
        }
    }
}

fn within_link_distance(mol: &Molecule, a: usize, b: usize) -> bool {
    match (mol.coord3(a), mol.coord3(b)) {
        (Some(p), Some(q)) => p.distance(q) < MAX_LINK_DISTANCE,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{Atom, Element};
    use crate::core::geometry::Point3;
    use crate::core::residue::Chain;
    use crate::core::site::AtomSite;

    /// One atom: name, element, position, alt-loc.
    type Spec = (&'static str, u8, [f64; 3], Option<char>);

    /// A one-chain molecule from `(residue name, atoms)` pairs.
    fn build(residues: &[(&str, &[Spec])]) -> Molecule {
        let mut mol = Molecule::new();
        let (mut sites, mut coords, mut ranges) = (Vec::new(), Vec::new(), Vec::new());
        for (ix, (name, atoms)) in residues.iter().enumerate() {
            let start = mol.num_atoms();
            for &(atom_name, z, [x, y, zc], alt_loc) in atoms.iter() {
                mol.add_atom(Atom::new(Element::new(z).unwrap()));
                coords.push(Point3::new(x, y, zc));
                sites.push(AtomSite {
                    name: Some(atom_name.to_string()),
                    alt_loc,
                    ..AtomSite::default()
                });
            }
            ranges.push(Residue {
                name: name.to_string(),
                sequence: ix as i32 + 1,
                insertion_code: None,
                label_seq: None,
                chain_ix: 0,
                is_hetero: false,
                atoms: start..mol.num_atoms(),
            });
        }
        mol.set_coords3(coords).unwrap();
        mol.set_sites(sites).unwrap();
        let chain = Chain {
            id: "A".to_string(),
            label_id: None,
            residues: 0..ranges.len(),
        };
        mol.set_topology(vec![chain], ranges).unwrap();
        mol
    }

    fn order_between(mol: &Molecule, a: usize, b: usize) -> Option<BondOrder> {
        mol.graph().get_bond(a, b).map(|ix| mol.bond(ix).order())
    }

    const ALA: &[Spec] = &[
        ("N", 7, [0.0, 0.0, 0.0], None),
        ("CA", 6, [1.46, 0.0, 0.0], None),
        ("C", 6, [2.0, 1.4, 0.0], None),
        ("O", 8, [1.3, 2.4, 0.0], None),
        ("CB", 6, [2.0, -0.8, 1.2], None),
    ];

    #[test]
    fn test_alanine_gets_its_four_heavy_atom_bonds_with_a_double_carbonyl() {
        let mut mol = build(&[("ALA", ALA)]);
        add_residue_template_bonds(&mut mol);
        assert_eq!(mol.num_bonds(), 4);
        assert_eq!(order_between(&mol, 0, 1), Some(BondOrder::Single));
        assert_eq!(order_between(&mol, 1, 2), Some(BondOrder::Single));
        assert_eq!(order_between(&mol, 2, 3), Some(BondOrder::Double));
        assert_eq!(order_between(&mol, 1, 4), Some(BondOrder::Single));
    }

    #[test]
    fn test_phenylalanine_ring_is_kekule_not_aromatic() {
        let ring = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"];
        let atoms: Vec<Spec> = ["N", "CA", "C", "O", "CB"]
            .iter()
            .chain(ring.iter())
            .map(|&n| (n, 6, [0.0; 3], None))
            .collect();
        let mut mol = build(&[("PHE", &atoms)]);
        add_residue_template_bonds(&mut mol);

        let orders: Vec<BondOrder> = (0..6)
            .map(|i| order_between(&mol, 5 + i, 5 + (i + 1) % 6).expect("ring bond"))
            .collect();
        let doubles = orders.iter().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 3);
        assert!(orders.windows(2).all(|w| w[0] != w[1]), "{orders:?}");
        assert!(mol.bonds().iter().all(|b| !b.is_aromatic()));
    }

    fn two_glycines(gap: f64) -> Molecule {
        let first: &[Spec] = &[
            ("N", 7, [0.0, 0.0, 0.0], None),
            ("CA", 6, [1.46, 0.0, 0.0], None),
            ("C", 6, [2.0, 1.4, 0.0], None),
        ];
        let second = [
            ("N", 7, [2.0 + gap, 1.4, 0.0], None),
            ("CA", 6, [3.5 + gap, 1.4, 0.0], None),
        ];
        build(&[("GLY", first), ("GLY", &second)])
    }

    #[test]
    fn test_consecutive_residues_bond_only_within_peptide_distance() {
        let mut close = two_glycines(1.33);
        add_residue_template_bonds(&mut close);
        assert_eq!(order_between(&close, 2, 3), Some(BondOrder::Single));

        let mut gapped = two_glycines(4.0);
        add_residue_template_bonds(&mut gapped);
        assert_eq!(order_between(&gapped, 2, 3), None);
    }

    #[test]
    fn test_alternate_locations_bond_only_within_their_own_conformer() {
        let atoms: &[Spec] = &[
            ("CA", 6, [0.0; 3], None),
            ("CB", 6, [1.5, 0.0, 0.0], Some('A')),
            ("CB", 6, [1.5, 0.3, 0.0], Some('B')),
            ("OG", 8, [2.5, 0.0, 0.0], Some('A')),
            ("OG", 8, [2.5, 0.3, 0.0], Some('B')),
        ];
        let mut mol = build(&[("SER", atoms)]);
        add_residue_template_bonds(&mut mol);
        assert!(order_between(&mol, 0, 1).is_some());
        assert!(order_between(&mol, 0, 2).is_some());
        assert!(order_between(&mol, 1, 3).is_some());
        assert!(order_between(&mol, 2, 4).is_some());
        assert_eq!(order_between(&mol, 1, 4), None);
        assert_eq!(order_between(&mol, 2, 3), None);
        assert_eq!(mol.num_bonds(), 4);
    }

    #[test]
    fn test_an_unknown_residue_gets_no_bonds() {
        let mut mol = build(&[("LIG", ALA)]);
        add_residue_template_bonds(&mut mol);
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_an_existing_bond_keeps_its_order() {
        let mut mol = build(&[("ALA", ALA)]);
        mol.add_bond(Bond::new(2, 3, BondOrder::Single)).unwrap();
        add_residue_template_bonds(&mut mol);
        assert_eq!(mol.num_bonds(), 4);
        assert_eq!(order_between(&mol, 2, 3), Some(BondOrder::Single));
    }

    #[test]
    fn test_his_takes_its_tautomer_from_the_ring_hydrogen_present() {
        let ring: Vec<Spec> = ["CB", "CG", "ND1", "CE1", "NE2", "CD2", "HD1"]
            .iter()
            .map(|&n| (n, if n.starts_with('H') { 1 } else { 6 }, [0.0; 3], None))
            .collect();
        let mut delta = build(&[("HIS", &ring)]);
        add_residue_template_bonds(&mut delta);
        assert_eq!(order_between(&delta, 3, 4), Some(BondOrder::Double)); // CE1=NE2

        let mut epsilon = build(&[("HIS", &ring[..6])]);
        add_residue_template_bonds(&mut epsilon);
        assert_eq!(order_between(&epsilon, 2, 3), Some(BondOrder::Double)); // ND1=CE1
    }

    #[test]
    fn test_nucleotides_link_through_o3_prime_and_p() {
        let first: &[Spec] = &[
            ("C3'", 6, [0.0; 3], None),
            ("O3'", 8, [1.4, 0.0, 0.0], None),
        ];
        let second: &[Spec] = &[
            ("P", 15, [3.0, 0.0, 0.0], None),
            ("OP1", 8, [3.5, 1.4, 0.0], None),
        ];
        let mut mol = build(&[("DA", first), ("DT", second)]);
        add_residue_template_bonds(&mut mol);
        assert_eq!(order_between(&mol, 1, 2), Some(BondOrder::Single));
        assert_eq!(order_between(&mol, 2, 3), Some(BondOrder::Double));
    }

    #[test]
    fn test_every_template_atom_pair_is_well_formed() {
        for name in [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
            "ILE", "LEU", "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "A", "C", "G", "U", "DA", "DC", "DG", "DT",
        ] {
            let (_, parts) = side_chain(name).expect(name);
            for spec in parts.iter().flat_map(|p| p.split_whitespace()) {
                let (a, b, _) = parse_bond(spec).expect(spec);
                assert!(!a.is_empty() && !b.is_empty() && a != b, "{name}: {spec}");
            }
        }
    }
}
