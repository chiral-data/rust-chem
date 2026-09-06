//! commonchem JSON (#229) — RDKit's structured interchange format: one JSON
//! document holding a `molecules` array, each molecule an atom list, a bond
//! list, optional conformers and optional properties.
//!
//! **Written as `commonchem` version 10; read as `commonchem` 10 or
//! `rdkitjson` 10–12.** These are the same schema family under two header
//! keys, and the split is real rather than cosmetic: RDKit 2025.03.3's
//! writer emits `{"rdkitjson":{"version":12}}` and nothing else, while its
//! reader accepts `commonchem` only at version 10. Writing the vendor-neutral
//! key is the point of the format; accepting the vendor one is the only way to
//! read what the reference implementation actually produces.
//!
//! **Aromaticity lives in a vendor extension, not in the bond orders.** This
//! is the format's single most surprising property and the reason this module
//! emits a `rdkitRepresentation` extension at all. Atoms and bonds are always
//! written in a Kekulé form; aromaticity is carried by `aromaticAtoms` and
//! `aromaticBonds` inside `extensions`. Drop the extension and RDKit reads
//! benzene back as `OC1=CC=CC=C1` — it does not re-perceive. And `bo: 4` is
//! **quadruple**, not aromatic: an all-`bo:4` ring reads back as
//! `[CH]1$[CH]$[CH]$[CH]$[CH]$[CH]$1`, so the molfile trick of a type-4 bond
//! (#194) has no analogue here and must not be improvised.
//!
//! **A ring this crate cannot kekulise still round-trips.** When
//! [`crate::io::aromaticity::kekulize`] returns `None` the aromatic bonds are
//! written as plain single bonds and the extension carries the truth anyway.
//! Verified against RDKit: an all-`bo:1` five-ring plus
//! `aromaticAtoms`/`aromaticBonds` reads back as `c1cc[nH]c1`, correctly
//! aromatic. That is what makes this the first format registered here where
//! aromaticity survives in both directions with no portability asterisk —
//! `write_sdf`'s equivalent fallback is a type-4 bond, which RDKit rejects
//! (#194).
//!
//! The fallback is rarely reached: `kekulize` solves every aromatic system
//! this crate has been pointed at, pyrrole included, so the neighbouring doc
//! comments in `io/sdf.rs` and `io/aromaticity.rs` naming pyrrole as the
//! standing failure are out of date. It is kept because the writer must have
//! *some* answer, and the point of this format is that its answer costs
//! nothing.
//!
//! **`defaults` is load-bearing on read, so it is always written.** A document
//! missing the block — or missing just `defaults.atom.stereo` inside it —
//! fails RDKit's reader with `Bad Format: bad stereo value for atom`. Every
//! other default is optional there; all of them are written here anyway, since
//! the atom and bond objects are encoded relative to them. A field absent from
//! an atom object means *the default*, never zero.
//!
//! **Bond stereo maps exactly, because neither side assumes CIP.** See
//! `bond_stereo_out` and `bond_stereo_in` below: this crate defines
//! [`BondStereo::E`]/[`BondStereo::Z`] relative to the lowest-indexed
//! substituent at each end (`smiles_writer::reference_substituent`), and
//! commonchem names its two reference atoms explicitly in `stereoAtoms`.
//! Those three are private, so they are named here rather than linked —
//! rustdoc rejects a public doc linking a private item, and CI builds docs
//! with `-D warnings`.
//! The two are the same convention, so `cis`/`trans` is a rename rather than a
//! reinterpretation — as long as the reference atoms travel with it, which is
//! what makes this the one place in the module that can be silently wrong.
//!
//! **What is not modelled.** `nRad` (radical count) is read and discarded —
//! `Atom` has no radical field. `cipRanks`, `cipCodes` and `atomRings` in the
//! extension are derived data this crate recomputes on demand, so they are
//! ignored on read and never written. `stereoGroups` and query features are
//! out of scope (they belong with #221's `StereoGroup` and the v0.12.0 query
//! wave respectively). `Molecule` holds one 3D conformer, so the first `dim:3`
//! conformer is kept and any others dropped; likewise the first `dim:2`.
//! An atom's implicit/explicit hydrogen split is a SMILES-syntax distinction
//! with no field here — both are written into `impHs` and read back as
//! implicit, which preserves the chemistry and not the spelling.
//!
//! **`impHs` reports what the molecule holds, including when that is wrong.**
//! This is the first format registered here with a per-atom hydrogen count, so
//! it is the first one that puts that number where another toolkit can see it
//! — and doing so immediately exposed #240, in which
//! [`Molecule::calculate_implicit_hydrogens`] inverts the charge sign and
//! hands every anion two extra hydrogens. `C(=O)[O-]` therefore writes
//! `impHs: 2` and RDKit reads it back as `[OH2-]`. That is deliberate: a
//! writer that quietly corrected the count on the way out would disagree with
//! this crate's own `formula()` and would have hidden the defect from the
//! oracle that found it. Fix #240, and this output becomes right with no
//! change here.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::core::atom::{Atom, Chirality, Element};
use crate::core::bond::{Bond, BondOrder, BondStereo};
use crate::core::geometry::{Point2, Point3};
use crate::core::molecule::Molecule;
use crate::io::errors::CommonchemError;
use crate::io::smiles_writer::reference_substituent;

/// The schema version this module writes, under the `commonchem` key.
const WRITTEN_VERSION: u32 = 10;

/// The extension RDKit uses to carry aromaticity, and the only one written
/// here. `formatVersion` is mandatory whenever an extension is present —
/// omitting it fails RDKit's reader with `Bad Format: missing format_version`.
const RDKIT_EXTENSION: &str = "rdkitRepresentation";
const RDKIT_EXTENSION_FORMAT_VERSION: u32 = 2;

// --- the document ------------------------------------------------------------
//
// These types describe the file, not the molecule. They are private on
// purpose: commonchem's shape (defaults-relative atoms, an index-pair bond
// list, a conformer array, a vendor extension) has nothing in common with
// `Molecule`'s, and deriving `Serialize` on the core types to avoid them would
// either produce the wrong document or push serde attributes onto types whose
// `Eq` derive and side-table design were settled in #174/#176.

#[derive(Debug, Serialize, Deserialize)]
struct Version {
    version: u32,
}

#[derive(Debug, Serialize, Deserialize)]
struct Document {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    commonchem: Option<Version>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    rdkitjson: Option<Version>,
    #[serde(default)]
    defaults: Defaults,
    #[serde(default)]
    molecules: Vec<JMol>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Defaults {
    atom: AtomDefaults,
    bond: BondDefaults,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct AtomDefaults {
    z: u8,
    imp_hs: u8,
    chg: i8,
    n_rad: u8,
    isotope: u16,
    stereo: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BondDefaults {
    bo: u32,
    stereo: String,
}

impl Default for Defaults {
    /// The schema's own defaults, which are also exactly what RDKit writes.
    /// Not `#[derive(Default)]`: a carbon's `z` is 6 and a single bond's `bo`
    /// is 1, so a derived all-zeroes default would silently turn every
    /// unwritten atom into a neutron.
    fn default() -> Self {
        Defaults {
            atom: AtomDefaults {
                z: 6,
                imp_hs: 0,
                chg: 0,
                n_rad: 0,
                isotope: 0,
                stereo: "unspecified".to_string(),
            },
            bond: BondDefaults {
                bo: 1,
                stereo: "unspecified".to_string(),
            },
        }
    }
}

#[derive(Debug, Default, Serialize, Deserialize)]
struct JMol {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    name: Option<String>,
    #[serde(default)]
    atoms: Vec<JAtom>,
    #[serde(default)]
    bonds: Vec<JBond>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    conformers: Vec<JConformer>,
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    properties: HashMap<String, serde_json::Value>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    extensions: Vec<JExtension>,
}

/// Every field is `Option` because absence means *the document's default*,
/// not zero. Resolution against [`Defaults`] happens in [`atom_in`].
#[derive(Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct JAtom {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    z: Option<u8>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    imp_hs: Option<u8>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    chg: Option<i8>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    n_rad: Option<u8>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    isotope: Option<u16>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    stereo: Option<String>,
}

#[derive(Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct JBond {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    bo: Option<u32>,
    atoms: Vec<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    stereo: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    stereo_atoms: Option<Vec<usize>>,
}

#[derive(Debug, Serialize, Deserialize)]
struct JConformer {
    dim: u8,
    coords: Vec<Vec<f64>>,
}

#[derive(Debug, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct JExtension {
    name: String,
    /// Mandatory when written; defaulted on read so an unfamiliar extension
    /// that omits it costs nothing — this module only ever reads the two
    /// aromatic arrays out of it.
    #[serde(default)]
    format_version: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    toolkit_version: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    aromatic_atoms: Vec<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    aromatic_bonds: Vec<usize>,
}

// --- reading -----------------------------------------------------------------

/// Parses a whole commonchem document into its molecules and their names.
///
/// Unlike `parse_gro`/`parse_xyz`/`parse_cml`, this takes the entire file
/// rather than one already-isolated record: `molecules` is an array inside a
/// single JSON value, and a JSON document is only valid whole. Splitting into
/// records is therefore this function's job, not the reader's or the
/// supplier's — see [`crate::io::reader::read_commonchem_with_options`].
pub fn parse_commonchem(text: &str) -> Result<Vec<(String, Molecule)>, CommonchemError> {
    let doc: Document = serde_json::from_str(text)?;
    check_version(&doc)?;

    let mut out = Vec::with_capacity(doc.molecules.len());
    for (ix, jmol) in doc.molecules.iter().enumerate() {
        let molecule = molecule_in(jmol, &doc.defaults)?;
        let name = jmol
            .name
            .clone()
            .unwrap_or_else(|| format!("Molecule_{}", ix + 1));
        out.push((name, molecule));
    }
    Ok(out)
}

/// `commonchem` is accepted at version 10 only and `rdkitjson` at 10 to 12 —
/// the same windows RDKit's own reader enforces, verified against 2025.03.3
/// rather than read off the schema page. Refusing `commonchem` 12 is agreement
/// with the reference implementation, not extra strictness.
fn check_version(doc: &Document) -> Result<(), CommonchemError> {
    match (&doc.commonchem, &doc.rdkitjson) {
        (Some(v), _) if v.version == WRITTEN_VERSION => Ok(()),
        (Some(v), _) => Err(CommonchemError::UnsupportedVersion {
            key: "commonchem".to_string(),
            version: v.version,
        }),
        (None, Some(v)) if (10..=12).contains(&v.version) => Ok(()),
        (None, Some(v)) => Err(CommonchemError::UnsupportedVersion {
            key: "rdkitjson".to_string(),
            version: v.version,
        }),
        (None, None) => Err(CommonchemError::MissingHeader),
    }
}

fn molecule_in(jmol: &JMol, defaults: &Defaults) -> Result<Molecule, CommonchemError> {
    let mut mol = Molecule::new();
    if let Some(name) = &jmol.name {
        mol.set_name(name.clone());
    }

    for jatom in &jmol.atoms {
        mol.add_atom(atom_in(jatom, &defaults.atom)?);
    }

    // Bonds first with no configuration, then stereo in a second pass: the
    // E/Z convention is defined by an atom's neighbours, so it cannot be
    // resolved until every bond exists.
    for jbond in &jmol.bonds {
        let (a, b) = bond_ends(jbond, mol.num_atoms())?;
        let order = bond_order_in(jbond.bo.unwrap_or(defaults.bond.bo))?;
        mol.add_bond(Bond::new(a, b, order))
            .map_err(|e| CommonchemError::ParseError(e.to_string()))?;
    }
    for (ix, jbond) in jmol.bonds.iter().enumerate() {
        let stated = jbond.stereo.as_deref().unwrap_or(&defaults.bond.stereo);
        let stereo = bond_stereo_in(&mol, ix, stated, jbond.stereo_atoms.as_deref())?;
        if stereo != BondStereo::None {
            let updated = mol.bond(ix).clone().with_stereo(stereo);
            *mol.bond_mut(ix) = updated;
        }
    }

    // The extension is applied last: it overrides the bond orders the Kekulé
    // form just installed, which is the whole reason it exists.
    apply_aromatic_extension(&mut mol, &jmol.extensions)?;

    let (two_d, three_d) = conformers_in(jmol, mol.num_atoms())?;
    if let Some(coords) = two_d {
        mol.set_coords(coords)
            .map_err(|e| CommonchemError::ParseError(e.to_string()))?;
    }
    if let Some(coords) = three_d {
        mol.set_coords3(coords)
            .map_err(|e| CommonchemError::ParseError(e.to_string()))?;
    }

    for (key, value) in &jmol.properties {
        mol.set_property(key.clone(), property_in(value));
    }

    Ok(mol)
}

fn atom_in(jatom: &JAtom, defaults: &AtomDefaults) -> Result<Atom, CommonchemError> {
    let z = jatom.z.unwrap_or(defaults.z);
    let element = Element::new(z).ok_or(CommonchemError::InvalidElement(z))?;

    let mut atom = Atom::new(element).with_charge(jatom.chg.unwrap_or(defaults.chg));

    // Isotope 0 is the schema's "no isotope", not mass number zero.
    let isotope = jatom.isotope.unwrap_or(defaults.isotope);
    if isotope != 0 {
        atom = atom.with_isotope(isotope);
    }

    let stereo = jatom.stereo.as_deref().unwrap_or(&defaults.stereo);
    atom = atom.with_chirality(match stereo {
        "unspecified" => Chirality::None,
        "cw" => Chirality::Clockwise,
        "ccw" => Chirality::CounterClockwise,
        other => return Err(CommonchemError::InvalidAtomStereo(other.to_string())),
    });

    atom.set_implicit_hydrogens(jatom.imp_hs.unwrap_or(defaults.imp_hs));
    Ok(atom)
}

fn bond_ends(jbond: &JBond, num_atoms: usize) -> Result<(usize, usize), CommonchemError> {
    let [a, b] = jbond.atoms[..] else {
        return Err(CommonchemError::ParseError(format!(
            "a bond needs exactly two atom indices, got {}",
            jbond.atoms.len()
        )));
    };
    for atom in [a, b] {
        if atom >= num_atoms {
            return Err(CommonchemError::BondIndexOutOfRange { atom, num_atoms });
        }
    }
    Ok((a, b))
}

fn bond_order_in(bo: u32) -> Result<BondOrder, CommonchemError> {
    // No aromatic case: `bo: 4` is a quadruple bond in this schema, and
    // aromaticity arrives through the extension instead.
    match bo {
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Quadruple),
        other => Err(CommonchemError::InvalidBondOrder(other)),
    }
}

/// `cis`/`trans` → [`BondStereo::Z`]/[`BondStereo::E`], reoriented onto this
/// crate's own reference substituents.
///
/// Both conventions say "these two named atoms are on the same side / opposite
/// sides"; they differ only in which atoms get named. Swapping one end's
/// reference atom for the other substituent inverts the relationship, so an
/// odd number of disagreeing ends flips the answer and an even number does not.
///
/// A stated `stereoAtoms` that this crate would not have chosen is normal, not
/// an error — RDKit picks its own. Its absence is treated as agreement, which
/// is the only assumption available and is what RDKit's own reader does.
fn bond_stereo_in(
    mol: &Molecule,
    bond_ix: usize,
    stated: &str,
    stereo_atoms: Option<&[usize]>,
) -> Result<BondStereo, CommonchemError> {
    let base = match stated {
        "unspecified" => return Ok(BondStereo::None),
        "either" => return Ok(BondStereo::Unspecified),
        "cis" => BondStereo::Z,
        "trans" => BondStereo::E,
        other => return Err(CommonchemError::InvalidBondStereo(other.to_string())),
    };

    let bond = mol.bond(bond_ix);
    let ours = [
        reference_substituent(mol, bond.atom1(), bond.atom2()),
        reference_substituent(mol, bond.atom2(), bond.atom1()),
    ];
    let flips = match stereo_atoms {
        Some([a, b]) => usize::from(ours[0] != Some(*a)) + usize::from(ours[1] != Some(*b)),
        _ => 0,
    };

    Ok(if flips % 2 == 1 {
        match base {
            BondStereo::Z => BondStereo::E,
            _ => BondStereo::Z,
        }
    } else {
        base
    })
}

/// Sets aromaticity from `rdkitRepresentation`, the only channel this format
/// has for it. Bonds get [`BondOrder::Aromatic`] rather than keeping the
/// Kekulé order they were read with, so the result matches what
/// `parse_smiles("c1ccccc1")` builds and a round trip is exact rather than
/// merely equivalent.
fn apply_aromatic_extension(
    mol: &mut Molecule,
    extensions: &[JExtension],
) -> Result<(), CommonchemError> {
    for ext in extensions.iter().filter(|e| e.name == RDKIT_EXTENSION) {
        for &ix in &ext.aromatic_atoms {
            if ix >= mol.num_atoms() {
                return Err(CommonchemError::ExtensionIndexOutOfRange {
                    what: "aromaticAtoms",
                    index: ix,
                });
            }
            mol.atom_mut(ix).set_aromatic(true);
        }
        for &ix in &ext.aromatic_bonds {
            if ix >= mol.num_bonds() {
                return Err(CommonchemError::ExtensionIndexOutOfRange {
                    what: "aromaticBonds",
                    index: ix,
                });
            }
            let bond = mol.bond_mut(ix);
            bond.set_order(BondOrder::Aromatic);
            bond.set_aromatic(true);
        }
    }
    Ok(())
}

type Conformers = (Option<Vec<Point2>>, Option<Vec<Point3>>);

/// The first `dim:2` and the first `dim:3` conformer. Later ones are dropped —
/// `Molecule` holds one of each.
fn conformers_in(jmol: &JMol, num_atoms: usize) -> Result<Conformers, CommonchemError> {
    let mut two_d = None;
    let mut three_d = None;

    for conf in &jmol.conformers {
        if conf.coords.len() != num_atoms {
            return Err(CommonchemError::ConformerLengthMismatch {
                dim: conf.dim,
                expected: num_atoms,
                got: conf.coords.len(),
            });
        }
        match conf.dim {
            2 if two_d.is_none() => {
                let mut points = Vec::with_capacity(num_atoms);
                for c in &conf.coords {
                    points.push(Point2::new(axis(c, 0, 2)?, axis(c, 1, 2)?));
                }
                two_d = Some(points);
            }
            3 if three_d.is_none() => {
                let mut points = Vec::with_capacity(num_atoms);
                for c in &conf.coords {
                    points.push(Point3::new(axis(c, 0, 3)?, axis(c, 1, 3)?, axis(c, 2, 3)?));
                }
                three_d = Some(points);
            }
            2 | 3 => {}
            other => return Err(CommonchemError::UnsupportedConformerDim(other)),
        }
    }
    Ok((two_d, three_d))
}

fn axis(coord: &[f64], ix: usize, dim: u8) -> Result<f64, CommonchemError> {
    coord
        .get(ix)
        .copied()
        .ok_or(CommonchemError::ConformerLengthMismatch {
            dim,
            expected: dim as usize,
            got: coord.len(),
        })
}

/// commonchem properties are typed; this crate's are `String`-valued, so a
/// number or a bool arrives here as its rendering. Lossy in this direction
/// only — see the module doc.
fn property_in(value: &serde_json::Value) -> String {
    match value {
        serde_json::Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

// --- writing -----------------------------------------------------------------

/// Serialises named molecules into one commonchem document.
///
/// One document for all records, not one per record — unlike every other
/// format registered here, the output has no concatenable per-record framing.
pub fn write_commonchem(records: &[(String, Molecule)]) -> String {
    let doc = Document {
        commonchem: Some(Version {
            version: WRITTEN_VERSION,
        }),
        rdkitjson: None,
        defaults: Defaults::default(),
        molecules: records
            .iter()
            .map(|(name, mol)| molecule_out(name, mol))
            .collect(),
    };
    // Compact, the way RDKit writes it. `jq` is the pretty-printer.
    let mut text =
        serde_json::to_string(&doc).expect("a document built from a Molecule serialises");
    text.push('\n');
    text
}

fn molecule_out(name: &str, mol: &Molecule) -> JMol {
    let defaults = Defaults::default();

    // A Kekulé form when one exists; plain single bonds for the aromatic ones
    // when it does not. Either way the extension below carries the truth, so
    // unlike `write_sdf` this fallback loses nothing.
    let kekulised = crate::io::aromaticity::kekulize(mol);

    let atoms = mol
        .atoms()
        .iter()
        .map(|a| atom_out(a, &defaults.atom))
        .collect();

    let bonds = (0..mol.num_bonds())
        .map(|ix| {
            let bond = mol.bond(ix);
            let order = match &kekulised {
                Some(orders) => orders[ix],
                None if bond.order() == BondOrder::Aromatic => BondOrder::Single,
                None => bond.order(),
            };
            let (stereo, stereo_atoms) = bond_stereo_out(mol, bond);
            JBond {
                bo: Some(bond_order_out(order)).filter(|bo| *bo != defaults.bond.bo),
                atoms: vec![bond.atom1(), bond.atom2()],
                stereo,
                stereo_atoms,
            }
        })
        .collect();

    let mut conformers = Vec::new();
    if let Some(coords) = mol.coords3() {
        conformers.push(JConformer {
            dim: 3,
            coords: coords.iter().map(|p| vec![p.x, p.y, p.z]).collect(),
        });
    }
    if let Some(coords) = mol.coords() {
        conformers.push(JConformer {
            dim: 2,
            coords: coords.iter().map(|p| vec![p.x, p.y]).collect(),
        });
    }

    JMol {
        name: (!name.is_empty()).then(|| name.to_string()),
        atoms,
        bonds,
        conformers,
        properties: mol
            .properties()
            .iter()
            .map(|(k, v)| (k.clone(), serde_json::Value::String(v.clone())))
            .collect(),
        extensions: aromatic_extension(mol).into_iter().collect(),
    }
}

/// Every field is omitted when it equals the document's default, which is what
/// makes the output the same shape RDKit's is.
fn atom_out(atom: &Atom, defaults: &AtomDefaults) -> JAtom {
    JAtom {
        z: Some(atom.atomic_number()).filter(|z| *z != defaults.z),
        // Implicit and explicit both land here: commonchem has one field for
        // "hydrogens not in the graph", and that is what they both are.
        imp_hs: Some(atom.total_hydrogens()).filter(|h| *h != defaults.imp_hs),
        chg: Some(atom.formal_charge()).filter(|c| *c != defaults.chg),
        n_rad: None,
        isotope: atom.isotope().filter(|i| *i != defaults.isotope),
        stereo: match atom.chirality() {
            Chirality::Clockwise => Some("cw".to_string()),
            Chirality::CounterClockwise => Some("ccw".to_string()),
            // `None` and `Unspecified` collapse: the schema has exactly three
            // atom stereo values (`unspecified`/`cw`/`ccw`), with nothing for
            // "a stereocentre whose configuration is unknown".
            Chirality::None | Chirality::Unspecified => None,
        },
    }
}

fn bond_order_out(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        // Unreachable via `molecule_out`, which kekulises first and falls back
        // to single. Mapping it to 1 rather than 4 keeps that true if a future
        // caller skips the fallback: `bo: 4` would claim a quadruple bond.
        BondOrder::Aromatic => 1,
    }
}

/// The `stereo`/`stereoAtoms` pair for one bond. See the module doc: naming
/// the reference atoms is what makes `cis`/`trans` mean the same thing as this
/// crate's `Z`/`E` rather than merely resemble it.
fn bond_stereo_out(mol: &Molecule, bond: &Bond) -> (Option<String>, Option<Vec<usize>>) {
    match bond.stereo() {
        BondStereo::None => (None, None),
        BondStereo::Unspecified => (Some("either".to_string()), None),
        stereo => {
            let ends = (
                reference_substituent(mol, bond.atom1(), bond.atom2()),
                reference_substituent(mol, bond.atom2(), bond.atom1()),
            );
            match ends {
                (Some(a), Some(b)) => {
                    let word = if stereo == BondStereo::Z {
                        "cis"
                    } else {
                        "trans"
                    };
                    (Some(word.to_string()), Some(vec![a, b]))
                }
                // A double bond with no substituent to reference at one end
                // cannot state a configuration at all. Writing `cis` with no
                // `stereoAtoms` would be an unanchored claim, so write nothing.
                _ => (None, None),
            }
        }
    }
}

/// The `rdkitRepresentation` extension, when there is aromaticity to record.
///
/// Only the two aromatic arrays are written. `cipRanks`/`cipCodes`/`atomRings`
/// are derived data this crate recomputes, and emitting a stale copy of
/// something a reader can recalculate is how a file starts lying.
fn aromatic_extension(mol: &Molecule) -> Option<JExtension> {
    let aromatic_atoms: Vec<usize> = (0..mol.num_atoms())
        .filter(|&i| mol.atom(i).is_aromatic())
        .collect();
    let aromatic_bonds: Vec<usize> = (0..mol.num_bonds())
        .filter(|&i| {
            let bond = mol.bond(i);
            bond.is_aromatic() || bond.order() == BondOrder::Aromatic
        })
        .collect();

    if aromatic_atoms.is_empty() && aromatic_bonds.is_empty() {
        return None;
    }
    Some(JExtension {
        name: RDKIT_EXTENSION.to_string(),
        format_version: RDKIT_EXTENSION_FORMAT_VERSION,
        toolkit_version: Some(format!("chem {}", env!("CARGO_PKG_VERSION"))),
        aromatic_atoms,
        aromatic_bonds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::smiles::parse_smiles;
    use crate::io::smiles_writer::write_smiles_for_molecule_canonical;

    /// Write one molecule, read it back, and hand back the single record.
    fn round_trip(mol: &Molecule) -> Molecule {
        let text = write_commonchem(&[("probe".to_string(), mol.clone())]);
        let mut back = parse_commonchem(&text).expect("our own output parses");
        assert_eq!(back.len(), 1, "one record in, one record out");
        back.remove(0).1
    }

    fn from_smiles(smiles: &str) -> Molecule {
        parse_smiles(smiles).expect("valid SMILES")
    }

    #[test]
    fn test_round_trip_preserves_the_molecule() {
        let mol = from_smiles("CCO");
        let back = round_trip(&mol);
        assert_eq!(back.num_atoms(), 3);
        assert_eq!(back.num_bonds(), 2);
        assert_eq!(
            write_smiles_for_molecule_canonical(&back),
            write_smiles_for_molecule_canonical(&mol)
        );
    }

    #[test]
    fn test_the_name_travels_with_the_record() {
        let text = write_commonchem(&[("ethanol".to_string(), from_smiles("CCO"))]);
        assert!(text.contains(r#""name":"ethanol""#), "{text}");
        assert_eq!(parse_commonchem(&text).unwrap()[0].0, "ethanol");
    }

    #[test]
    fn test_atoms_are_written_relative_to_the_defaults() {
        // A neutral, non-isotopic, unstereo carbon with no hydrogens is the
        // document's default in every field, so it must serialise as `{}`.
        // Getting this wrong is invisible in a round trip -- both sides agree
        // -- and only shows up as output RDKit's reader disagrees with.
        let mut mol = Molecule::new();
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_atom(Atom::new(Element::carbon()));
        mol.add_bond(Bond::new(0, 1, BondOrder::Single)).unwrap();
        let text = write_commonchem(&[(String::new(), mol)]);

        assert!(text.contains(r#""atoms":[{},{}]"#), "{text}");
        // A single bond is `bo: 1`, the default, so only `atoms` survives.
        assert!(text.contains(r#""bonds":[{"atoms":[0,1]}]"#), "{text}");
    }

    #[test]
    fn test_the_defaults_block_is_always_written() {
        // RDKit's reader fails with `bad stereo value for atom` on a document
        // with no `defaults`, or with `defaults.atom.stereo` missing from one.
        // Neither is optional, however sparse the rest of the output is.
        let text = write_commonchem(&[(String::new(), from_smiles("C"))]);
        assert!(text.contains(r#""defaults":{"atom":{"z":6"#), "{text}");
        assert!(text.contains(r#""stereo":"unspecified""#), "{text}");
    }

    #[test]
    fn test_formal_charge_survives() {
        let mol = from_smiles("CC(=O)[O-]");
        let back = round_trip(&mol);
        let charges: Vec<i8> = back.atoms().iter().map(|a| a.formal_charge()).collect();
        assert!(charges.contains(&-1), "{charges:?}");
    }

    #[test]
    fn test_isotope_survives_and_zero_means_absent() {
        let mol = from_smiles("[13CH4]");
        let back = round_trip(&mol);
        assert_eq!(back.atom(0).isotope(), Some(13));

        // Isotope 0 is the schema's "no isotope", not mass number zero, so an
        // ordinary atom must not come back claiming to be isotope 0.
        let plain = round_trip(&from_smiles("C"));
        assert_eq!(plain.atom(0).isotope(), None);
    }

    #[test]
    fn test_atom_chirality_survives_both_ways() {
        for smiles in ["N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"] {
            let mol = from_smiles(smiles);
            let back = round_trip(&mol);
            let before: Vec<_> = mol.atoms().iter().map(|a| a.chirality()).collect();
            let after: Vec<_> = back.atoms().iter().map(|a| a.chirality()).collect();
            assert_eq!(before, after, "{smiles}");
        }
    }

    #[test]
    fn test_bond_stereo_survives_both_ways() {
        for (smiles, expected) in [("F/C=C/F", BondStereo::E), ("F/C=C\\F", BondStereo::Z)] {
            let back = round_trip(&from_smiles(smiles));
            let stereo = back
                .bonds()
                .iter()
                .find(|b| b.order() == BondOrder::Double)
                .expect("a double bond")
                .stereo();
            assert_eq!(stereo, expected, "{smiles}");
        }
    }

    #[test]
    fn test_bond_stereo_is_anchored_to_its_reference_atoms() {
        // The case where lowest-index and CIP priority disagree, and the one
        // that would pass anyway if `stereoAtoms` were ignored on both sides.
        // Atom 0 is F and atom 4 is Cl on the first carbon; a reader that took
        // `cis`/`trans` as CIP-relative would flip this.
        let mol = from_smiles("F/C(Cl)=C(Br)/I");
        let text = write_commonchem(&[("probe".to_string(), mol.clone())]);
        assert!(text.contains(r#""stereoAtoms""#), "{text}");
        assert_eq!(
            round_trip(&mol)
                .bonds()
                .iter()
                .find(|b| b.order() == BondOrder::Double)
                .unwrap()
                .stereo(),
            mol.bonds()
                .iter()
                .find(|b| b.order() == BondOrder::Double)
                .unwrap()
                .stereo()
        );
    }

    #[test]
    fn test_a_disagreeing_stereo_atom_flips_the_configuration() {
        // One end's reference atom swapped for the other substituent inverts
        // the relationship; this is the arithmetic in `bond_stereo_in`, and
        // it is the whole reason `stereoAtoms` is read rather than assumed.
        let doc = |stereo_atoms: &str| {
            format!(
                r#"{{"commonchem":{{"version":10}},"molecules":[{{"atoms":[{{"z":9}},{{}},{{"z":17}},{{}},{{"z":35}},{{"z":53}}],
                "bonds":[{{"atoms":[0,1]}},{{"atoms":[1,2]}},{{"bo":2,"atoms":[1,3],"stereo":"cis","stereoAtoms":{stereo_atoms}}},{{"atoms":[3,4]}},{{"atoms":[3,5]}}]}}]}}"#
            )
        };
        let stereo_of = |text: String| {
            parse_commonchem(&text).expect("valid document")[0]
                .1
                .bonds()
                .iter()
                .find(|b| b.order() == BondOrder::Double)
                .unwrap()
                .stereo()
        };
        // [0, 4] are this crate's own reference substituents -- no flip.
        assert_eq!(stereo_of(doc("[0,4]")), BondStereo::Z);
        // [2, 4] disagrees at one end -- one flip.
        assert_eq!(stereo_of(doc("[2,4]")), BondStereo::E);
        // [2, 5] disagrees at both -- two flips, back to where it started.
        assert_eq!(stereo_of(doc("[2,5]")), BondStereo::Z);
    }

    #[test]
    fn test_aromaticity_rides_on_the_extension_not_the_bond_orders() {
        let text = write_commonchem(&[("benzene".to_string(), from_smiles("c1ccccc1"))]);
        // Kekulé in the bond list: `bo: 4` would claim a quadruple bond.
        assert!(!text.contains(r#""bo":4"#), "{text}");
        assert!(text.contains(r#""aromaticAtoms":[0,1,2,3,4,5]"#), "{text}");
        assert!(text.contains(r#""formatVersion":2"#), "{text}");

        let back = round_trip(&from_smiles("c1ccccc1"));
        assert!(back.atoms().iter().all(|a| a.is_aromatic()));
        assert!(
            back.bonds()
                .iter()
                .all(|b| b.order() == BondOrder::Aromatic)
        );
    }

    #[test]
    fn test_a_molecule_with_no_aromaticity_writes_no_extension() {
        let text = write_commonchem(&[("ethanol".to_string(), from_smiles("CCO"))]);
        assert!(!text.contains("extensions"), "{text}");
    }

    #[test]
    fn test_an_unkekulisable_ring_still_round_trips() {
        // Built rather than parsed, because `kekulize` solves every real
        // aromatic molecule this crate can read -- pyrrole included, despite
        // what `io/sdf.rs`'s doc comment still says. Three aromatic carbons
        // each needing one more bond is an odd number of atoms to pair up, so
        // no perfect matching exists and the writer must take its fallback.
        let mut mol = Molecule::new();
        for _ in 0..3 {
            let mut atom = Atom::new(Element::carbon());
            atom.set_implicit_hydrogens(1);
            mol.add_atom(atom);
        }
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            mol.add_bond(Bond::new(a, b, BondOrder::Aromatic).with_aromatic(true))
                .unwrap();
        }
        for ix in 0..3 {
            mol.atom_mut(ix).set_aromatic(true);
        }
        assert!(
            crate::io::aromaticity::kekulize(&mol).is_none(),
            "this fixture must be the case kekulize cannot solve, or the test \
             is checking the happy path instead of the fallback"
        );

        let text = write_commonchem(&[("probe".to_string(), mol.clone())]);
        // Every bond object carries no `bo` at all, which is the default of
        // 1. Matched against the bond list rather than the whole document,
        // since `defaults` legitimately contains `"bo":1` itself.
        assert!(
            text.contains(r#""bonds":[{"atoms":[0,1]},{"atoms":[1,2]},{"atoms":[2,0]}]"#),
            "aromatic bonds fall back to single, never to the quadruple `bo: 4`: {text}"
        );
        assert!(text.contains(r#""aromaticBonds":[0,1,2]"#), "{text}");

        let back = round_trip(&mol);
        assert_eq!(back.num_atoms(), 3);
        assert!(
            back.bonds()
                .iter()
                .all(|b| b.order() == BondOrder::Aromatic)
        );
        assert!(back.atoms().iter().all(|a| a.is_aromatic()));
    }

    #[test]
    fn test_2d_and_3d_conformers_coexist_on_one_molecule() {
        let mut mol = from_smiles("CO");
        mol.set_coords(vec![Point2::new(0.0, 0.0), Point2::new(1.4, 0.0)])
            .unwrap();
        mol.set_coords3(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.4, 0.0, 0.5)])
            .unwrap();

        let text = write_commonchem(&[("probe".to_string(), mol)]);
        assert!(text.contains(r#""dim":3"#), "{text}");
        assert!(text.contains(r#""dim":2"#), "{text}");

        let back = parse_commonchem(&text).unwrap().remove(0).1;
        assert_eq!(back.coord(1), Some(Point2::new(1.4, 0.0)));
        assert_eq!(back.coord3(1), Some(Point3::new(1.4, 0.0, 0.5)));
    }

    #[test]
    fn test_properties_survive_as_strings() {
        let mut mol = from_smiles("C");
        mol.set_property("comment".to_string(), "kept".to_string());
        let back = round_trip(&mol);
        assert_eq!(back.property("comment"), Some("kept"));
    }

    #[test]
    fn test_a_typed_property_arrives_as_its_rendering() {
        // commonchem's properties are typed and this crate's are strings, so
        // this direction is lossy by construction -- stated rather than
        // silently stringified.
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{}],"bonds":[],
            "properties":{"n":3,"ok":true,"s":"text"}}]}"#;
        let mol = &parse_commonchem(text).unwrap()[0].1;
        assert_eq!(mol.property("n"), Some("3"));
        assert_eq!(mol.property("ok"), Some("true"));
        assert_eq!(mol.property("s"), Some("text"));
    }

    #[test]
    fn test_one_document_holds_every_record() {
        let records = vec![
            ("a".to_string(), from_smiles("C")),
            ("b".to_string(), from_smiles("O")),
            ("c".to_string(), from_smiles("N")),
        ];
        let text = write_commonchem(&records);
        // One document, not three concatenated -- two JSON values back to
        // back are not a JSON document.
        assert_eq!(text.matches(r#""commonchem""#).count(), 1, "{text}");
        let back = parse_commonchem(&text).unwrap();
        assert_eq!(
            back.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>(),
            ["a", "b", "c"]
        );
    }

    #[test]
    fn test_an_unnamed_molecule_gets_a_positional_name() {
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{}],"bonds":[]},{"atoms":[{}],"bonds":[]}]}"#;
        let back = parse_commonchem(text).unwrap();
        assert_eq!(back[0].0, "Molecule_1");
        assert_eq!(back[1].0, "Molecule_2");
    }

    #[test]
    fn test_the_rdkitjson_header_is_accepted_on_read() {
        // What RDKit 2025.03.3 actually emits. Refusing it would mean this
        // reader could not read the reference implementation's own output.
        for version in [10, 11, 12] {
            let text = format!(
                r#"{{"rdkitjson":{{"version":{version}}},"molecules":[{{"atoms":[{{}}],"bonds":[]}}]}}"#
            );
            assert!(parse_commonchem(&text).is_ok(), "rdkitjson {version}");
        }
    }

    #[test]
    fn test_an_out_of_window_version_is_refused() {
        // RDKit rejects `commonchem` 12 too, so refusing it is agreement with
        // the reference implementation rather than extra strictness.
        let text = r#"{"commonchem":{"version":12},"molecules":[]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::UnsupportedVersion { .. })
        ));

        let text = r#"{"rdkitjson":{"version":13},"molecules":[]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::UnsupportedVersion { .. })
        ));
    }

    #[test]
    fn test_a_document_with_no_header_is_refused() {
        let text = r#"{"molecules":[{"atoms":[{}],"bonds":[]}]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::MissingHeader)
        ));
    }

    #[test]
    fn test_malformed_json_is_an_error_not_a_panic() {
        for text in [
            r#"{"commonchem":{"version":10},"#,
            "not json at all",
            "",
            "[]",
        ] {
            assert!(parse_commonchem(text).is_err(), "{text:?}");
        }
    }

    #[test]
    fn test_a_bond_index_past_the_end_is_refused() {
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{},{}],"bonds":[{"atoms":[0,7]}]}]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::BondIndexOutOfRange { atom: 7, .. })
        ));
    }

    #[test]
    fn test_an_extension_index_past_the_end_is_refused() {
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{}],"bonds":[],
            "extensions":[{"name":"rdkitRepresentation","formatVersion":2,"aromaticAtoms":[9]}]}]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::ExtensionIndexOutOfRange { .. })
        ));
    }

    #[test]
    fn test_a_short_conformer_is_refused() {
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{},{}],"bonds":[],
            "conformers":[{"dim":3,"coords":[[0,0,0]]}]}]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::ConformerLengthMismatch { .. })
        ));
    }

    #[test]
    fn test_bo_four_is_a_quadruple_bond_not_an_aromatic_one() {
        // The trap this format sets for anyone carrying molfile habits over:
        // RDKit reads an all-`bo:4` ring as `[CH]1$[CH]$...$1`, six quadruple
        // bonds, not benzene.
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{},{}],"bonds":[{"bo":4,"atoms":[0,1]}]}]}"#;
        let mol = &parse_commonchem(text).unwrap()[0].1;
        assert_eq!(mol.bond(0).order(), BondOrder::Quadruple);
        assert!(!mol.bond(0).is_aromatic());
    }

    #[test]
    fn test_an_unknown_stereo_word_is_refused() {
        let text = r#"{"commonchem":{"version":10},"molecules":[{"atoms":[{"stereo":"either"}],"bonds":[]}]}"#;
        assert!(matches!(
            parse_commonchem(text),
            Err(CommonchemError::InvalidAtomStereo(_))
        ));
    }
}
