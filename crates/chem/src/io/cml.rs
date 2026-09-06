//! CML (Chemical Markup Language) — XML. The first format in this crate
//! needing an XML parser at all (#228); see [`roxmltree`] usage below and
//! the dependency's own doc comment in `Cargo.toml` for why it is always
//! compiled in rather than gated behind a feature.
//!
//! **Scope: the `<molecule>` element only**, matching what `core` already
//! expresses — atoms, bonds, 2D-or-3D coordinates, formal charge, isotope.
//! CML's full spec also covers reactions, crystal structures, sequences and
//! spectra; none of that is attempted here, the same curation stance as
//! Mol2's SYBYL subset (#225) and PDBQT's AutoDock-type subset (#226): a
//! curated slice, stated plainly, not the whole specification.
//!
//! - `<atomArray><atom id="a1" elementType="C" .../></atomArray>` — `id` is
//!   a string key resolved only through `<bond>` references, never stored
//!   on the atom itself. `x2`/`y2` and `x3`/`y3`/`z3` populate whichever
//!   coordinate channel(s) are present — both if both are, since
//!   [`crate::core::molecule::Molecule`] already carries them as
//!   independent optional channels. **Coordinates are only ever set when
//!   every atom in the record has them** — a file where some atoms carry
//!   `x2`/`y2` and others don't is unusual/malformed, and this reader does
//!   not attempt a partial recovery that could silently place some atoms
//!   at a false origin.
//! - `<bondArray><bond atomRefs2="a1 a2" order="1"/></bondArray>` — `order`
//!   accepts both conventions real CML files use: numeric (`1`/`2`/`3`/`4`)
//!   and single-letter (`S`/`D`/`T`/`A`), case-insensitively. `A`/aromatic
//!   is what lets [`crate::io::format::Carries::AROMATICITY`] survive at
//!   all here — through the same `bond.order() ==
//!   BondOrder::Aromatic` path [`crate::io::format::held`] already uses for
//!   SDF's type-4 bonds, not a separate atom-level aromaticity flag. Write
//!   always uses the numeric convention.
//!
//! **Namespace-agnostic.** Real CML files declare `xmlns` inconsistently —
//! the modern schema URI, an old DTD reference, or nothing at all — so
//! elements and attributes are matched by local name regardless of
//! namespace, not validated against any particular schema. Consistent with
//! this crate's established "pragmatic over spec-perfect" stance (GRO's
//! declared-count trust, Mol2's curated type table).
//!
//! **Out of scope, stated rather than silently dropped:** `<crystal>`
//! (CML's own periodic-structure extension — the same deferral reasoning
//! as extended XYZ's `Lattice=`, #222), reactions, sequences, spectra, atom
//! and bond stereo elements, and any markush/query construct. No
//! `RESIDUES`, `UNIT_CELL`, `PARTIAL_CHARGE`, `OCCUPANCY`, `B_FACTOR`,
//! `PROPERTIES`, `STEREO_ATOM` or `STEREO_BOND` is claimed.
//!
//! A file may hold several `<molecule>` elements — bare, or wrapped in a
//! container such as `<cml>`/`<list>` — splitting those is the reader's/
//! supplier's job ([`crate::io::reader::read_cml_with_options`],
//! [`crate::io::supplier::CmlSupplier`]), the same division every prior
//! format has. Nested `<molecule>` elements (CML technically allows this
//! for coordination complexes) are out of scope there — each record is
//! assumed flat. This module parses and writes exactly one already-
//! isolated `<molecule>` element.
//!
//! **A multi-molecule *write* wraps every record in one `<cml>` root**
//! (see [`crate::io::format`]'s `write_cml_records` and
//! [`crate::io::supplier::CmlWriter`]) — several sibling `<molecule>`
//! elements with no enclosing root is not valid XML at all (exactly one
//! root element is required), and a strict parser stops after the first.
//! Found during this story's own verification: OpenBabel's CML reader
//! silently returned only one molecule from a first, unwrapped attempt at
//! multi-record output, correctly rejecting the rest of the file as
//! "extra content at the end of the document."

use std::collections::HashMap;

use roxmltree::Document;

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::{Point2, Point3};
use crate::core::molecule::Molecule;
use crate::io::errors::CmlError;

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

fn bond_order_from_str(order: &str) -> Option<BondOrder> {
    match order.trim().to_ascii_uppercase().as_str() {
        "1" | "S" => Some(BondOrder::Single),
        "2" | "D" => Some(BondOrder::Double),
        "3" | "T" => Some(BondOrder::Triple),
        "4" => Some(BondOrder::Quadruple),
        "A" | "AROMATIC" => Some(BondOrder::Aromatic),
        _ => None,
    }
}

fn bond_order_to_str(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "1",
        BondOrder::Double => "2",
        BondOrder::Triple => "3",
        BondOrder::Quadruple => "4",
        BondOrder::Aromatic => "A",
    }
}

fn escape_xml_attr(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn find_child<'a, 'input>(
    node: roxmltree::Node<'a, 'input>,
    name: &str,
) -> Option<roxmltree::Node<'a, 'input>> {
    node.children()
        .find(|n| n.is_element() && n.has_tag_name(name))
}

/// Parses one already-isolated `<molecule>` element.
pub fn parse_cml(text: &str) -> Result<Molecule, CmlError> {
    let doc = Document::parse(text).map_err(|e| CmlError::ParseError(e.to_string()))?;
    let molecule_node = doc.root_element();
    if !molecule_node.has_tag_name("molecule") {
        return Err(CmlError::ParseError(format!(
            "expected a <molecule> root element, found <{}>",
            molecule_node.tag_name().name()
        )));
    }

    let mut mol = Molecule::new();
    if let Some(title) = molecule_node
        .attribute("title")
        .or_else(|| molecule_node.attribute("id"))
        && !title.is_empty()
    {
        mol.set_name(title.to_string());
    }

    let mut coords2: Vec<Option<Point2>> = Vec::new();
    let mut coords3: Vec<Option<Point3>> = Vec::new();
    let mut id_to_index: HashMap<String, usize> = HashMap::new();

    if let Some(atom_array) = find_child(molecule_node, "atomArray") {
        for atom_node in atom_array
            .children()
            .filter(|n| n.is_element() && n.has_tag_name("atom"))
        {
            let id = atom_node
                .attribute("id")
                .ok_or_else(|| CmlError::InvalidAtomElement("atom has no id".to_string()))?;
            let element_type = atom_node.attribute("elementType").ok_or_else(|| {
                CmlError::InvalidAtomElement(format!("atom {id:?} has no elementType"))
            })?;
            let element = element_from_symbol(element_type)
                .ok_or_else(|| CmlError::InvalidElement(element_type.to_string()))?;

            let mut atom = Atom::new(element);
            if let Some(charge) = atom_node
                .attribute("formalCharge")
                .and_then(|s| s.parse::<i8>().ok())
            {
                atom = atom.with_charge(charge);
            }
            if let Some(isotope) = atom_node
                .attribute("isotopeNumber")
                .and_then(|s| s.parse::<u16>().ok())
            {
                atom = atom.with_isotope(isotope);
            }
            if let Some(h) = atom_node
                .attribute("hydrogenCount")
                .and_then(|s| s.parse::<u8>().ok())
            {
                atom.set_implicit_hydrogens(h);
            }

            let atom_idx = mol.add_atom(atom);
            id_to_index.insert(id.to_string(), atom_idx);

            let xy2 = (
                atom_node
                    .attribute("x2")
                    .and_then(|s| s.parse::<f64>().ok()),
                atom_node
                    .attribute("y2")
                    .and_then(|s| s.parse::<f64>().ok()),
            );
            coords2.push(match xy2 {
                (Some(x), Some(y)) => Some(Point2::new(x, y)),
                _ => None,
            });

            let xyz3 = (
                atom_node
                    .attribute("x3")
                    .and_then(|s| s.parse::<f64>().ok()),
                atom_node
                    .attribute("y3")
                    .and_then(|s| s.parse::<f64>().ok()),
                atom_node
                    .attribute("z3")
                    .and_then(|s| s.parse::<f64>().ok()),
            );
            coords3.push(match xyz3 {
                (Some(x), Some(y), Some(z)) => Some(Point3::new(x, y, z)),
                _ => None,
            });
        }
    }

    if !coords2.is_empty() && coords2.iter().all(Option::is_some) {
        let flat: Vec<Point2> = coords2.into_iter().map(|p| p.unwrap()).collect();
        mol.set_coords(flat)
            .map_err(|e| CmlError::ParseError(e.to_string()))?;
    }
    if !coords3.is_empty() && coords3.iter().all(Option::is_some) {
        let flat: Vec<Point3> = coords3.into_iter().map(|p| p.unwrap()).collect();
        mol.set_coords3(flat)
            .map_err(|e| CmlError::ParseError(e.to_string()))?;
    }

    if let Some(bond_array) = find_child(molecule_node, "bondArray") {
        for bond_node in bond_array
            .children()
            .filter(|n| n.is_element() && n.has_tag_name("bond"))
        {
            let refs = bond_node
                .attribute("atomRefs2")
                .ok_or_else(|| CmlError::InvalidBondElement("bond has no atomRefs2".to_string()))?;
            let mut ids = refs.split_whitespace();
            let (Some(a), Some(b)) = (ids.next(), ids.next()) else {
                return Err(CmlError::InvalidBondElement(refs.to_string()));
            };
            let a_idx = *id_to_index
                .get(a)
                .ok_or_else(|| CmlError::UnknownAtomReference(a.to_string()))?;
            let b_idx = *id_to_index
                .get(b)
                .ok_or_else(|| CmlError::UnknownAtomReference(b.to_string()))?;
            let order_str = bond_node.attribute("order").unwrap_or("1");
            let order = bond_order_from_str(order_str).ok_or_else(|| {
                CmlError::InvalidBondElement(format!("unrecognised bond order {order_str:?}"))
            })?;
            mol.add_bond(Bond::new(a_idx, b_idx, order))
                .map_err(|e| CmlError::ParseError(e.to_string()))?;
        }
    }

    Ok(mol)
}

/// Writes one `<molecule>` element.
pub fn write_cml(mol: &Molecule) -> String {
    let mut out = String::from("<molecule");
    out.push_str(" xmlns=\"http://www.xml-cml.org/schema\"");
    if let Some(name) = mol.name() {
        out.push_str(&format!(" title=\"{}\"", escape_xml_attr(name)));
    }
    out.push_str(">\n");

    out.push_str("  <atomArray>\n");
    for (i, atom) in mol.atoms().iter().enumerate() {
        out.push_str(&format!(
            "    <atom id=\"a{}\" elementType=\"{}\"",
            i + 1,
            atom.element().symbol()
        ));
        if let Some(p) = mol.coord(i) {
            out.push_str(&format!(" x2=\"{:.4}\" y2=\"{:.4}\"", p.x, p.y));
        }
        if let Some(p) = mol.coord3(i) {
            out.push_str(&format!(
                " x3=\"{:.4}\" y3=\"{:.4}\" z3=\"{:.4}\"",
                p.x, p.y, p.z
            ));
        }
        if atom.formal_charge() != 0 {
            out.push_str(&format!(" formalCharge=\"{}\"", atom.formal_charge()));
        }
        if let Some(isotope) = atom.isotope() {
            out.push_str(&format!(" isotopeNumber=\"{isotope}\""));
        }
        if atom.implicit_hydrogens() > 0 {
            out.push_str(&format!(" hydrogenCount=\"{}\"", atom.implicit_hydrogens()));
        }
        out.push_str("/>\n");
    }
    out.push_str("  </atomArray>\n");

    out.push_str("  <bondArray>\n");
    for bond in mol.bonds() {
        out.push_str(&format!(
            "    <bond atomRefs2=\"a{} a{}\" order=\"{}\"/>\n",
            bond.atom1() + 1,
            bond.atom2() + 1,
            bond_order_to_str(bond.order()),
        ));
    }
    out.push_str("  </bondArray>\n");

    out.push_str("</molecule>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const ETHANE_CML: &str = "\
<molecule id=\"m1\" title=\"ethane\" xmlns=\"http://www.xml-cml.org/schema\">
  <atomArray>
    <atom id=\"a1\" elementType=\"C\" x3=\"0.0\" y3=\"0.0\" z3=\"0.0\"/>
    <atom id=\"a2\" elementType=\"C\" x3=\"1.5\" y3=\"0.0\" z3=\"0.0\"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2=\"a1 a2\" order=\"1\"/>
  </bondArray>
</molecule>
";

    #[test]
    fn test_a_single_molecule_round_trips() {
        let mol = parse_cml(ETHANE_CML).expect("valid CML");
        assert_eq!(mol.name(), Some("ethane"));
        assert_eq!(mol.num_atoms(), 2);
        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.coord3(0), Some(Point3::new(0.0, 0.0, 0.0)));
        assert_eq!(mol.bonds()[0].order(), BondOrder::Single);

        let written = write_cml(&mol);
        let back = parse_cml(&written).expect("round trips");
        assert_eq!(back.num_atoms(), 2);
        assert_eq!(back.num_bonds(), 1);
        assert_eq!(back.coord3(0), Some(Point3::new(0.0, 0.0, 0.0)));
    }

    #[test]
    fn test_2d_coordinates_round_trip() {
        let text = "\
<molecule>
  <atomArray>
    <atom id=\"a1\" elementType=\"C\" x2=\"0.0\" y2=\"0.0\"/>
    <atom id=\"a2\" elementType=\"O\" x2=\"1.2\" y2=\"0.5\"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2=\"a1 a2\" order=\"2\"/>
  </bondArray>
</molecule>
";
        let mol = parse_cml(text).expect("valid CML");
        assert_eq!(mol.coord(1), Some(Point2::new(1.2, 0.5)));
        assert_eq!(mol.bonds()[0].order(), BondOrder::Double);
    }

    #[test]
    fn test_charge_isotope_and_aromatic_bond_round_trip() {
        let text = "\
<molecule>
  <atomArray>
    <atom id=\"a1\" elementType=\"N\" formalCharge=\"1\" isotopeNumber=\"15\"/>
    <atom id=\"a2\" elementType=\"C\"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2=\"a1 a2\" order=\"A\"/>
  </bondArray>
</molecule>
";
        let mol = parse_cml(text).expect("valid CML");
        assert_eq!(mol.atoms()[0].formal_charge(), 1);
        assert_eq!(mol.atoms()[0].isotope(), Some(15));
        assert_eq!(mol.bonds()[0].order(), BondOrder::Aromatic);

        let written = write_cml(&mol);
        assert!(written.contains("formalCharge=\"1\""), "{written}");
        assert!(written.contains("isotopeNumber=\"15\""), "{written}");
        assert!(written.contains("order=\"A\""), "{written}");
        let back = parse_cml(&written).expect("round trips");
        assert_eq!(back.atoms()[0].isotope(), Some(15));
        assert_eq!(back.bonds()[0].order(), BondOrder::Aromatic);
    }

    #[test]
    fn test_lowercase_letter_bond_orders_are_accepted() {
        let text = "\
<molecule>
  <atomArray>
    <atom id=\"a1\" elementType=\"C\"/>
    <atom id=\"a2\" elementType=\"C\"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2=\"a1 a2\" order=\"d\"/>
  </bondArray>
</molecule>
";
        let mol = parse_cml(text).expect("valid CML");
        assert_eq!(mol.bonds()[0].order(), BondOrder::Double);
    }

    #[test]
    fn test_partial_2d_coordinates_are_not_claimed() {
        // One atom has x2/y2, the other doesn't -- an inconsistent, unusual
        // file this reader does not attempt to patch with a false origin.
        let text = "\
<molecule>
  <atomArray>
    <atom id=\"a1\" elementType=\"C\" x2=\"0.0\" y2=\"0.0\"/>
    <atom id=\"a2\" elementType=\"C\"/>
  </atomArray>
</molecule>
";
        let mol = parse_cml(text).expect("valid CML");
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_a_bond_referencing_an_unknown_atom_id_is_a_clear_error() {
        let text = "\
<molecule>
  <atomArray>
    <atom id=\"a1\" elementType=\"C\"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2=\"a1 a99\" order=\"1\"/>
  </bondArray>
</molecule>
";
        let err = parse_cml(text).unwrap_err();
        assert!(matches!(err, CmlError::UnknownAtomReference(_)), "{err}");
    }

    #[test]
    fn test_malformed_xml_is_a_clear_error_not_a_panic() {
        let err = parse_cml("<molecule><atomArray><atom id=\"a1\"").unwrap_err();
        assert!(matches!(err, CmlError::ParseError(_)), "{err}");
    }

    #[test]
    fn test_an_unrecognised_element_symbol_is_a_clear_error() {
        let text = "\
<molecule>
  <atomArray>
    <atom id=\"a1\" elementType=\"Xx\"/>
  </atomArray>
</molecule>
";
        let err = parse_cml(text).unwrap_err();
        assert!(matches!(err, CmlError::InvalidElement(_)), "{err}");
    }
}
