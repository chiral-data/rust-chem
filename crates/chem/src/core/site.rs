//! Per-atom data that comes from a file rather than from chemistry.

use std::fmt;

/// The columns a structural file carries about an atom that are not properties
/// of the atom itself.
///
/// An element, a formal charge and an isotope are facts about a species. A
/// temperature factor, an occupancy and a name like `CA` are facts about one
/// *observation* of it — they come from a PDB, an mmCIF, a Mol2 or a PQR, and
/// two files describing the same molecule will disagree about them.
///
/// Held on [`Molecule`] as a side table indexed in parallel with the atoms,
/// for the same reason coordinates are: [`Atom`] derives `Eq`, and these are
/// floats.
///
/// # Every field is independently optional
///
/// Formats supply different subsets. A PDB gives occupancy and B-factor but no
/// partial charge; a Mol2 gives a charge and neither of the others; a PQR gives
/// a charge and a radius. One `Option` around the whole struct would force a
/// writer to invent values for the columns its source never supplied, which is
/// how a plausible-looking file full of zeroes gets written.
///
/// [`Molecule`]: crate::core::molecule::Molecule
/// [`Atom`]: crate::core::atom::Atom
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AtomSite {
    /// The name as written — `CA`, `OD1`, `HB2`. Not an element symbol, and
    /// unique only within its residue rather than within the molecule.
    pub name: Option<String>,

    /// Alternate-location indicator, for atoms modelled in more than one
    /// position. `None` is the common case.
    pub alt_loc: Option<char>,

    /// Partial (fractional) charge, as distinct from the integer formal charge
    /// on [`Atom`](crate::core::atom::Atom).
    ///
    /// The value alone is close to meaningless without the model that produced
    /// it — Mol2 names the method in the file, and a docking format's charges
    /// come from a different calculation entirely. Recording that provenance is
    /// deliberately deferred until a format has to round-trip it.
    pub partial_charge: Option<f64>,

    /// Fraction of the structure in which this atom occupies this position.
    /// Usually 1.0; less where alternate locations split it.
    pub occupancy: Option<f64>,

    /// Temperature factor.
    ///
    /// Predicted structures reuse this column for a per-atom confidence score,
    /// which is why zeroing it on write is a data-loss bug rather than a
    /// rounding choice. It must survive a round trip.
    pub b_factor: Option<f64>,

    /// Atomic radius, in Ångström. Carried by PQR, and by little else.
    pub radius: Option<f64>,
}

impl AtomSite {
    /// A site with nothing recorded — the same as [`Default`], named for the
    /// places where `AtomSite::empty()` reads better than `default()`.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Whether this site records anything at all.
    ///
    /// A table of empty sites is indistinguishable from no table, and a writer
    /// can use this to decide whether a column is worth emitting.
    pub fn is_empty(&self) -> bool {
        self.name.is_none()
            && self.alt_loc.is_none()
            && self.partial_charge.is_none()
            && self.occupancy.is_none()
            && self.b_factor.is_none()
            && self.radius.is_none()
    }
}

impl fmt::Display for AtomSite {
    /// Only the fields that are set, so a debug print of a PDB atom does not
    /// scroll past four `None`s to reach the two values that matter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        let mut field = |f: &mut fmt::Formatter<'_>, text: String| -> fmt::Result {
            if !first {
                write!(f, " ")?;
            }
            first = false;
            write!(f, "{text}")
        };

        if let Some(name) = &self.name {
            field(f, format!("name={name}"))?;
        }
        if let Some(alt) = self.alt_loc {
            field(f, format!("alt={alt}"))?;
        }
        if let Some(q) = self.partial_charge {
            field(f, format!("q={q:.4}"))?;
        }
        if let Some(occ) = self.occupancy {
            field(f, format!("occ={occ:.2}"))?;
        }
        if let Some(b) = self.b_factor {
            field(f, format!("b={b:.2}"))?;
        }
        if let Some(r) = self.radius {
            field(f, format!("r={r:.3}"))?;
        }
        if first {
            write!(f, "(empty)")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_new_site_records_nothing() {
        let site = AtomSite::default();
        assert!(site.is_empty());
        assert_eq!(site, AtomSite::empty());
        assert_eq!(site.name, None);
        assert_eq!(site.b_factor, None);
    }

    #[test]
    fn test_fields_are_independently_optional() {
        // The point of six `Option`s rather than one around the struct: a PDB
        // supplies occupancy and B-factor, a Mol2 supplies a charge, and
        // neither should force the other to invent a value.
        let from_pdb = AtomSite {
            name: Some("CA".to_string()),
            occupancy: Some(1.0),
            b_factor: Some(23.45),
            ..AtomSite::default()
        };
        assert!(!from_pdb.is_empty());
        assert_eq!(from_pdb.partial_charge, None);

        let from_mol2 = AtomSite {
            partial_charge: Some(-0.4157),
            ..AtomSite::default()
        };
        assert!(!from_mol2.is_empty());
        assert_eq!(from_mol2.occupancy, None);
        assert_eq!(from_mol2.b_factor, None);
    }

    #[test]
    fn test_display_shows_only_what_is_set() {
        let site = AtomSite {
            name: Some("CA".to_string()),
            b_factor: Some(23.45),
            ..AtomSite::default()
        };
        assert_eq!(site.to_string(), "name=CA b=23.45");
        assert_eq!(AtomSite::empty().to_string(), "(empty)");
    }
}
