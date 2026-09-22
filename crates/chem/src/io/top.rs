//! GROMACS topology (`.top`) — a plain-text, `[ section ]`-delimited
//! format: `[ moleculetype ]` blocks (each with its own local `[ atoms ]`/
//! `[ bonds ]`/`[ angles ]`/`[ dihedrals ]`/`[ exclusions ]`), preceded by
//! parameter-level directives (`[ defaults ]`, `[ atomtypes ]`, ...) and
//! followed by system-level `[ system ]`/`[ molecules ]` (#323).
//!
//! **Topology only**, the same milestone-wide scope PSF/PRMTOP already
//! settled on: no `[ atomtypes ]`/`[ bondtypes ]`/`[ pairtypes ]`/
//! `[ angletypes ]`/`[ dihedraltypes ]`/`[ nonbond_params ]`/`[ cmaptypes ]`
//! parameter tables, no `[ pairs ]` (1-4 nonbonded pairs — not the same
//! thing as `exclusions`), no `[ constraints ]`/`[ settles ]`/
//! `[ position_restraints ]`/`[ dihedral_restraints ]`/`[ virtual_sites* ]`/
//! `[ cmap ]`. Every one of these is recognized and skipped to the next
//! `[ ` header, never stored, the same treatment PRMTOP gives its own
//! force-constant sections. `[ atoms ]` must state mass inline as a
//! consequence: GROMACS lets a line omit mass/charge and fall back to
//! `[ atomtypes ]`'s own default, a table this reader never parses, so an
//! omitted mass is a clear [`TopError::InvalidElement`]-adjacent error
//! rather than a guess. Element inference reuses PSF's own
//! `element_from_mass` verbatim — GROMACS states no element any more
//! directly than PSF does (its `type` column is a force-field-specific
//! string like `opls_135`).
//!
//! **A `.top` file legitimately holds several `[ moleculetype ]` blocks
//! already** — real practice, not an artifact of this crate's own
//! round-trip invariant. [`parse_top`] returns one [`Molecule`] per
//! `[ moleculetype ]`, named after its own stated name (unlike PSF/PRMTOP,
//! which have no per-topology name at all). Parameter-level sections
//! before the first `[ moleculetype ]` and system-level sections after the
//! last are recognized and skipped inline.
//!
//! **`#include` is recorded as unresolved, never followed.** A `.top` that
//! `#include`s a force-field tree or a shared molecule fragment (e.g.
//! water from a `tip3p.itp`) simply won't produce whatever that include
//! would have defined — a real, honest limitation, not a bug, and the
//! option that keeps this reader pure (no filesystem access, so it works
//! the same in a wasm build).
//!
//! **`#ifdef`/`#ifndef`/`#else`/`#endif` are evaluated against only this
//! same buffer's own `#define`s/`#undef`s** — no externally-supplied `-D`
//! flags, since this crate never runs `grompp`; an unresolved macro is
//! simply "not defined", ordinary cpp semantics. See this module's own
//! `preprocess`. Reading every branch unconditionally was rejected: a
//! real conditional block can state different structural content (not
//! just parameters) per branch, which unconditional reading would
//! double-count.
//!
//! **`[ exclusions ]` maps directly onto
//! [`crate::core::force_field::ForceFieldTopology::exclusions`]** — one
//! atom followed by 1+ partners, expanded into pairs. This is only the
//! file's *explicit* extra exclusions: GROMACS also auto-generates
//! exclusions out to `nrexcl` bonds from the bond graph at `grompp` time,
//! which this reader does not compute (no data-model slot for `nrexcl`,
//! and deriving the implied set would be inventing data the file never
//! literally stated) — the same spirit as CIF core never expanding
//! symmetry.
//!
//! **Dihedral `funct` codes decide proper vs. improper**, unlike PRMTOP's
//! sign trick: funct 2 and 4 are improper, every other defined code is
//! proper -- see this module's own `is_improper_funct`.
//!
//! No chain concept — same placeholder-chain-id pattern as PRMTOP, via
//! `cif_model`'s own `group_into_chains_and_residues`. Bonds go through
//! [`Molecule`]'s ordinary bond list, not `ForceFieldTopology` — same
//! precedent as PSF/PRMTOP/PDB's `CONECT` — and
//! `calculate_implicit_hydrogens_where_bonded` runs once a moleculetype's
//! `[ bonds ]` (its complete real bond list) is loaded.
//!
//! **`.top` is shared with PRMTOP**, sometimes saved under the same
//! extension. [`is_gromacs_top`] disambiguates by content, mirroring
//! `cif_core::is_small_molecule_cif` exactly.

use std::collections::HashSet;

use crate::core::atom::Atom;
use crate::core::bond::{Bond, BondOrder};
use crate::core::elements::ATOMIC_MASSES;
use crate::core::force_field::{ForceFieldAtom, ForceFieldTopology};
use crate::core::molecule::Molecule;
use crate::core::site::AtomSite;
use crate::io::cif_model::{ResidueKey, group_into_chains_and_residues};
use crate::io::errors::TopError;

/// GROMACS states no chain/segment concept at all -- every atom in a
/// moleculetype gets this same placeholder when grouped into
/// chains/residues.
const PLACEHOLDER_CHAIN_ID: &str = "";

/// Whether `text` looks like a GROMACS topology rather than an AMBER-style
/// PRMTOP saved under the same `.top` extension. GROMACS content never
/// uses `%` anywhere in its grammar; a real (or this crate's own) PRMTOP
/// always starts with a mandatory `%VERSION` line. `true` (GROMACS) unless
/// a `%`-prefixed marker is seen first -- `.top` already resolves to
/// `Format::PRMTOP` by default (registered first), so this only overrides
/// that default with real evidence, the same "override, don't guess"
/// shape as `cif_core::is_small_molecule_cif`.
///
/// `pub`, not `pub(crate)`: called from both `io::open::resolve_format`
/// (inside this crate) and `bin/chem/stream.rs::resolve_format_from_content`
/// (the CLI binary, a separate crate target that only sees this library's
/// public API).
pub fn is_gromacs_top(text: &str) -> bool {
    for line in text.lines() {
        let trimmed = line.trim_start();
        if trimmed.is_empty() {
            continue;
        }
        return !trimmed.starts_with('%');
    }
    true
}

/// GROMACS's dihedral function-type table: 2 (harmonic improper) and 4
/// (periodic improper) are impropers; every other defined code (1 proper
/// periodic, 3 Ryckaert-Bellemans, 5 Fourier, 8 tabulated, 9 proper
/// "multiple" -- several lines naming the same four atoms, stacking
/// several Fourier terms, needing no special handling since
/// `ForceFieldTopology::dihedrals` already tolerates duplicate
/// quadruplets --, 10 restricted, 11 combined bending-torsion) is a
/// proper torsion.
fn is_improper_funct(funct: u32) -> bool {
    matches!(funct, 2 | 4)
}

/// Runs GROMACS's own preprocessor layer over `text` once, spanning the
/// whole buffer (preprocessing state can cross section boundaries): strips
/// `;` comments, joins `\`-continued lines, drops `#include` lines
/// (recorded as unresolved -- never opened, see the module doc), and
/// evaluates `#ifdef`/`#ifndef`/`#else`/`#endif` against only this same
/// buffer's own `#define`s/`#undef`s (a name never locally defined is
/// simply "not defined"). Output is plain text, one cleaned line per row,
/// ready for section-walking in [`parse_top`].
fn preprocess(text: &str) -> String {
    struct Frame {
        condition: bool,
        took_else: bool,
    }

    let mut defines: HashSet<String> = HashSet::new();
    let mut stack: Vec<Frame> = Vec::new();
    let mut pending: Option<String> = None;
    let mut out = String::new();

    for raw_line in text.lines() {
        let line = match raw_line.split_once(';') {
            Some((before, _)) => before,
            None => raw_line,
        };
        let line = line.trim_end();
        let (line, continues) = match line.strip_suffix('\\') {
            Some(stripped) => (stripped, true),
            None => (line, false),
        };

        let joined = match pending.take() {
            Some(mut acc) => {
                acc.push(' ');
                acc.push_str(line.trim_start());
                acc
            }
            None => line.to_string(),
        };
        if continues {
            pending = Some(joined);
            continue;
        }

        let trimmed = joined.trim();

        if let Some(name) = trimmed.strip_prefix("#ifdef") {
            let name = name.trim();
            stack.push(Frame {
                condition: defines.contains(name),
                took_else: false,
            });
            continue;
        }
        if let Some(name) = trimmed.strip_prefix("#ifndef") {
            let name = name.trim();
            stack.push(Frame {
                condition: !defines.contains(name),
                took_else: false,
            });
            continue;
        }
        if trimmed == "#else" {
            if let Some(frame) = stack.last_mut() {
                frame.took_else = true;
            }
            continue;
        }
        if trimmed == "#endif" {
            stack.pop();
            continue;
        }

        let active = stack.iter().all(|f| {
            if f.took_else {
                !f.condition
            } else {
                f.condition
            }
        });

        if let Some(rest) = trimmed.strip_prefix("#define") {
            if active && let Some(name) = rest.split_whitespace().next() {
                defines.insert(name.to_string());
            }
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("#undef") {
            if active && let Some(name) = rest.split_whitespace().next() {
                defines.remove(name);
            }
            continue;
        }
        if trimmed.starts_with("#include") {
            continue;
        }
        if active && !trimmed.is_empty() {
            out.push_str(trimmed);
            out.push('\n');
        }
    }

    out
}

/// Recognizes a `[ name ]` header line, returning its normalized
/// (lowercased, trimmed) name.
fn section_header(line: &str) -> Option<String> {
    let inner = line.strip_prefix('[')?.strip_suffix(']')?;
    Some(inner.trim().to_lowercase())
}

fn parse_i64(token: &str) -> Result<i64, TopError> {
    token
        .parse()
        .map_err(|_| TopError::ParseError(format!("expected an integer, got {token:?}")))
}

/// Converts a moleculetype-local 1-based atom serial into a 0-based index,
/// checked against `num_atoms` rather than trusted.
fn atom_index(serial: i64, num_atoms: usize) -> Result<usize, TopError> {
    if serial < 1 || serial as usize > num_atoms {
        return Err(TopError::ParseError(format!(
            "atom index {serial} is out of range for {num_atoms} atoms"
        )));
    }
    Ok(serial as usize - 1)
}

/// One `[ moleculetype ]` block under construction.
struct Building {
    name: String,
    mol: Molecule,
    sites: Vec<AtomSite>,
    ff_atoms: Vec<ForceFieldAtom>,
    keys: Vec<ResidueKey>,
    angles: Vec<[usize; 3]>,
    dihedrals: Vec<[usize; 4]>,
    impropers: Vec<[usize; 4]>,
    exclusions: Vec<[usize; 2]>,
}

impl Building {
    fn new() -> Self {
        Self {
            name: String::new(),
            mol: Molecule::new(),
            sites: Vec::new(),
            ff_atoms: Vec::new(),
            keys: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            impropers: Vec::new(),
            exclusions: Vec::new(),
        }
    }

    fn set_name_and_nrexcl(&mut self, line: &str) -> Result<(), TopError> {
        // `nrexcl` (the rest of the line) has no data-model slot -- this
        // crate does not compute nrexcl-implied exclusions, see the module
        // doc -- so only the name is kept.
        self.name = line
            .split_whitespace()
            .next()
            .ok_or_else(|| TopError::ParseError(format!("empty moleculetype header: {line:?}")))?
            .to_string();
        Ok(())
    }

    fn push_atom_line(&mut self, line: &str) -> Result<(), TopError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 8 {
            return Err(TopError::InvalidAtomLine(line.to_string()));
        }
        let stated_nr: usize = fields[0]
            .parse()
            .map_err(|_| TopError::InvalidAtomLine(line.to_string()))?;
        if stated_nr != self.sites.len() + 1 {
            return Err(TopError::InvalidAtomLine(line.to_string()));
        }
        let atom_type = fields[1];
        let resnr: i32 = fields[2]
            .parse()
            .map_err(|_| TopError::InvalidAtomLine(line.to_string()))?;
        let residue = fields[3];
        let atom_name = fields[4];
        // fields[5] (cgnr, charge-group number) has no data-model slot.
        let charge: f64 = fields[6]
            .parse()
            .map_err(|_| TopError::InvalidAtomLine(line.to_string()))?;
        let mass: f64 = fields[7]
            .parse()
            .map_err(|_| TopError::InvalidAtomLine(line.to_string()))?;

        let element =
            crate::io::psf::element_from_mass(mass).ok_or(TopError::InvalidElement(mass))?;
        self.mol.add_atom(Atom::new(element));
        self.sites.push(AtomSite {
            name: Some(atom_name.to_string()),
            ..AtomSite::empty()
        });
        self.ff_atoms.push(ForceFieldAtom {
            atom_type: Some(atom_type.to_string()),
            mass: Some(mass),
            partial_charge: Some(charge),
        });
        self.keys.push(ResidueKey {
            chain_id: PLACEHOLDER_CHAIN_ID.to_string(),
            name: residue.to_string(),
            sequence: resnr,
            insertion_code: None,
            is_hetero: false,
        });
        Ok(())
    }

    fn push_bond_line(&mut self, line: &str) -> Result<(), TopError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(TopError::ParseError(format!("invalid bond line: {line:?}")));
        }
        let num_atoms = self.mol.num_atoms();
        let a = atom_index(parse_i64(fields[0])?, num_atoms)?;
        let b = atom_index(parse_i64(fields[1])?, num_atoms)?;
        self.mol
            .add_bond(Bond::new(a, b, BondOrder::Single))
            .map_err(|e| TopError::ParseError(e.to_string()))?;
        Ok(())
    }

    fn push_angle_line(&mut self, line: &str) -> Result<(), TopError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            return Err(TopError::ParseError(format!(
                "invalid angle line: {line:?}"
            )));
        }
        let num_atoms = self.mol.num_atoms();
        self.angles.push([
            atom_index(parse_i64(fields[0])?, num_atoms)?,
            atom_index(parse_i64(fields[1])?, num_atoms)?,
            atom_index(parse_i64(fields[2])?, num_atoms)?,
        ]);
        Ok(())
    }

    fn push_dihedral_line(&mut self, line: &str) -> Result<(), TopError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 5 {
            return Err(TopError::ParseError(format!(
                "dihedral line missing its function-type column: {line:?}"
            )));
        }
        let num_atoms = self.mol.num_atoms();
        let atoms = [
            atom_index(parse_i64(fields[0])?, num_atoms)?,
            atom_index(parse_i64(fields[1])?, num_atoms)?,
            atom_index(parse_i64(fields[2])?, num_atoms)?,
            atom_index(parse_i64(fields[3])?, num_atoms)?,
        ];
        let funct: u32 = fields[4].parse().map_err(|_| {
            TopError::ParseError(format!("invalid dihedral function type: {line:?}"))
        })?;
        if is_improper_funct(funct) {
            self.impropers.push(atoms);
        } else {
            self.dihedrals.push(atoms);
        }
        Ok(())
    }

    fn push_exclusion_line(&mut self, line: &str) -> Result<(), TopError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(TopError::ParseError(format!(
                "invalid exclusions line: {line:?}"
            )));
        }
        let num_atoms = self.mol.num_atoms();
        let atom = atom_index(parse_i64(fields[0])?, num_atoms)?;
        for &partner_field in &fields[1..] {
            let partner = atom_index(parse_i64(partner_field)?, num_atoms)?;
            self.exclusions.push([atom, partner]);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<(String, Molecule), TopError> {
        if self.mol.num_atoms() == 0 {
            return Err(TopError::NoAtoms);
        }
        self.mol
            .set_sites(self.sites)
            .map_err(|e| TopError::ParseError(e.to_string()))?;
        let (chains, residues) = group_into_chains_and_residues(&self.keys);
        self.mol
            .set_topology(chains, residues)
            .map_err(|e| TopError::ParseError(e.to_string()))?;
        // `[ bonds ]` is a moleculetype's complete real bond list, not a
        // partial one -- the same precedent PSF/PRMTOP's readers already
        // established for implying hydrogens from a complete bond graph.
        self.mol.calculate_implicit_hydrogens_where_bonded();
        self.mol
            .set_force_field(ForceFieldTopology {
                atoms: Some(self.ff_atoms),
                angles: self.angles,
                dihedrals: self.dihedrals,
                impropers: self.impropers,
                exclusions: self.exclusions,
                donors: Vec::new(),
                acceptors: Vec::new(),
            })
            .map_err(|e| TopError::ParseError(e.to_string()))?;
        self.mol.set_name(self.name.clone());
        Ok((self.name, self.mol))
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Section {
    Skip,
    MoleculeTypeHeader,
    Atoms,
    Bonds,
    Angles,
    Dihedrals,
    Exclusions,
}

/// Parses a whole GROMACS topology buffer into one [`Molecule`] per
/// `[ moleculetype ]` block, named after each block's own stated name.
pub fn parse_top(text: &str) -> Result<Vec<(String, Molecule)>, TopError> {
    let cleaned = preprocess(text);

    let mut results = Vec::new();
    let mut current: Option<Building> = None;
    let mut section = Section::Skip;

    for line in cleaned.lines() {
        if let Some(name) = section_header(line) {
            section = match name.as_str() {
                "moleculetype" => {
                    if let Some(building) = current.take() {
                        results.push(building.finish()?);
                    }
                    current = Some(Building::new());
                    Section::MoleculeTypeHeader
                }
                "atoms" => Section::Atoms,
                "bonds" => Section::Bonds,
                "angles" => Section::Angles,
                "dihedrals" => Section::Dihedrals,
                "exclusions" => Section::Exclusions,
                _ => Section::Skip,
            };
            continue;
        }

        match section {
            Section::MoleculeTypeHeader => {
                let building = current
                    .as_mut()
                    .expect("just set to Some above when this section was entered");
                building.set_name_and_nrexcl(line)?;
                section = Section::Skip;
            }
            Section::Atoms => current
                .as_mut()
                .ok_or_else(|| TopError::ParseError("[ atoms ] outside a moleculetype".into()))?
                .push_atom_line(line)?,
            Section::Bonds => current
                .as_mut()
                .ok_or_else(|| TopError::ParseError("[ bonds ] outside a moleculetype".into()))?
                .push_bond_line(line)?,
            Section::Angles => current
                .as_mut()
                .ok_or_else(|| TopError::ParseError("[ angles ] outside a moleculetype".into()))?
                .push_angle_line(line)?,
            Section::Dihedrals => current
                .as_mut()
                .ok_or_else(|| TopError::ParseError("[ dihedrals ] outside a moleculetype".into()))?
                .push_dihedral_line(line)?,
            Section::Exclusions => current
                .as_mut()
                .ok_or_else(|| {
                    TopError::ParseError("[ exclusions ] outside a moleculetype".into())
                })?
                .push_exclusion_line(line)?,
            Section::Skip => {}
        }
    }

    if let Some(building) = current.take() {
        results.push(building.finish()?);
    }
    if results.is_empty() {
        return Err(TopError::NoAtoms);
    }
    Ok(results)
}

fn write_moleculetype(out: &mut String, name: &str, mol: &Molecule) {
    let empty_ff = ForceFieldTopology::default();
    let force_field = mol.force_field().unwrap_or(&empty_ff);

    out.push_str("[ moleculetype ]\n; name  nrexcl\n");
    out.push_str(&format!("{name}    3\n\n"));

    out.push_str("[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n");
    for (i, atom) in mol.atoms().iter().enumerate() {
        let site = mol.site(i);
        let residue = mol.residue_of(i);
        let ff_atom = force_field.atoms.as_ref().and_then(|atoms| atoms.get(i));

        let resnr = residue.map(|r| r.sequence).unwrap_or(1);
        let resname = residue
            .map(|r| r.name.as_str())
            .filter(|s| !s.is_empty())
            .unwrap_or("UNK");
        let atom_name = site
            .and_then(|s| s.name.as_deref())
            .filter(|s| !s.is_empty())
            .unwrap_or(atom.element().symbol());
        let atom_type = ff_atom
            .and_then(|a| a.atom_type.as_deref())
            .filter(|s| !s.is_empty())
            .unwrap_or(atom.element().symbol());
        let charge = ff_atom.and_then(|a| a.partial_charge).unwrap_or(0.0);
        let mass = ff_atom
            .and_then(|a| a.mass)
            .unwrap_or(ATOMIC_MASSES[atom.element().atomic_number as usize]);

        out.push_str(&format!(
            "{:>6} {atom_type:<10} {resnr:>6} {resname:<8} {atom_name:<8} {:>6} {charge:>10.6} {mass:>10.4}\n",
            i + 1,
            i + 1,
        ));
    }
    out.push('\n');

    if mol.num_bonds() > 0 {
        out.push_str("[ bonds ]\n");
        for bond in mol.bonds() {
            out.push_str(&format!(
                "{:>6} {:>6}\n",
                bond.atom1() + 1,
                bond.atom2() + 1
            ));
        }
        out.push('\n');
    }

    if !force_field.angles.is_empty() {
        out.push_str("[ angles ]\n");
        for &[a, b, c] in &force_field.angles {
            out.push_str(&format!("{:>6} {:>6} {:>6}\n", a + 1, b + 1, c + 1));
        }
        out.push('\n');
    }

    if !force_field.dihedrals.is_empty() {
        out.push_str("[ dihedrals ]\n; proper (funct 9)\n");
        for &[a, b, c, d] in &force_field.dihedrals {
            out.push_str(&format!(
                "{:>6} {:>6} {:>6} {:>6} 9\n",
                a + 1,
                b + 1,
                c + 1,
                d + 1
            ));
        }
        out.push('\n');
    }

    if !force_field.impropers.is_empty() {
        out.push_str("[ dihedrals ]\n; improper (funct 4)\n");
        for &[a, b, c, d] in &force_field.impropers {
            out.push_str(&format!(
                "{:>6} {:>6} {:>6} {:>6} 4\n",
                a + 1,
                b + 1,
                c + 1,
                d + 1
            ));
        }
        out.push('\n');
    }

    if !force_field.exclusions.is_empty() {
        out.push_str("[ exclusions ]\n");
        let mut by_atom: Vec<Vec<usize>> = vec![Vec::new(); mol.num_atoms()];
        for &[a, b] in &force_field.exclusions {
            by_atom[a].push(b);
        }
        for (i, partners) in by_atom.iter_mut().enumerate() {
            if partners.is_empty() {
                continue;
            }
            partners.sort_unstable();
            let cols: String = partners.iter().map(|p| format!(" {:>6}", p + 1)).collect();
            out.push_str(&format!("{:>6}{cols}\n", i + 1));
        }
        out.push('\n');
    }
}

/// Writes `records` as a GROMACS topology: one `[ moleculetype ]` block per
/// record, then a truthful `[ system ]`/`[ molecules ]` footer -- one line
/// per record, each stating count `1`, since this crate tracks no real
/// system multiplicity (see the module doc). Not simulation-ready: no
/// parameters, no defaults, no real composition -- the same acknowledgment
/// PRMTOP's own writer already makes.
pub fn write_top(records: &[(String, Molecule)]) -> String {
    let mut out = String::new();
    for (name, mol) in records {
        write_moleculetype(&mut out, name, mol);
    }
    out.push_str("[ system ]\nWritten by chem\n\n");
    out.push_str("[ molecules ]\n; moleculetype  count\n");
    for (name, _) in records {
        out.push_str(&format!("{name:<15} 1\n"));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const UREA_LIKE: &str = include_str!("../../tests/corpus/top/urea_like.top");
    const MULTI_MOLECULETYPE: &str = include_str!("../../tests/corpus/top/multi_moleculetype.top");
    const CONDITIONAL: &str = include_str!("../../tests/corpus/top/conditional.top");

    #[test]
    fn test_a_single_moleculetype_round_trips() {
        let records = parse_top(UREA_LIKE).expect("valid TOP");
        assert_eq!(records.len(), 1);
        let (name, mol) = &records[0];
        assert_eq!(name, "Urea");
        assert!(!mol.has_coords3(), "TOP states no coordinates");
        assert_eq!(mol.num_bonds(), 3);

        let force_field = mol.force_field().expect("has a force field");
        assert_eq!(force_field.dihedrals.len(), 1, "funct 9 is proper");
        assert_eq!(force_field.impropers.len(), 1, "funct 4 is improper");
        assert!(!force_field.exclusions.is_empty());

        let written = write_top(&records);
        let back = parse_top(&written).expect("round trips");
        assert_eq!(back.len(), 1);
        assert_eq!(back[0].1.num_atoms(), mol.num_atoms());
        assert_eq!(back[0].1.num_bonds(), mol.num_bonds());
        assert_eq!(
            back[0].1.force_field().unwrap().dihedrals,
            force_field.dihedrals
        );
        assert_eq!(
            back[0].1.force_field().unwrap().impropers,
            force_field.impropers
        );

        let mut original: Vec<[usize; 2]> = force_field.exclusions.clone();
        let mut round_tripped: Vec<[usize; 2]> =
            back[0].1.force_field().unwrap().exclusions.clone();
        for pair in original.iter_mut().chain(round_tripped.iter_mut()) {
            pair.sort_unstable();
        }
        original.sort_unstable();
        round_tripped.sort_unstable();
        assert_eq!(round_tripped, original);
    }

    #[test]
    fn test_multiple_moleculetypes_become_separate_records() {
        let records = parse_top(MULTI_MOLECULETYPE).expect("valid TOP");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].0, "Ion");
        assert_eq!(records[1].0, "Water");
        assert_eq!(records[0].1.num_atoms(), 1);
        assert_eq!(records[1].1.num_atoms(), 3);
    }

    #[test]
    fn test_an_include_is_recorded_as_unresolved_not_followed() {
        // The fixture #includes a nonexistent file; a reader that tried to
        // open it would fail outright. This one simply never sees whatever
        // that file would have defined.
        let records = parse_top(UREA_LIKE).expect("valid TOP despite the #include");
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn test_ifdef_is_evaluated_against_in_file_defines_only() {
        let records = parse_top(CONDITIONAL).expect("valid TOP");
        let (_, mol) = &records[0];
        let angles = &mol.force_field().unwrap().angles;
        // FLEXIBLE is #define'd in the fixture, so the #ifdef branch's
        // angle (1 2 3, atoms 0/1/2) survives. If both branches were read
        // (conditional evaluation broken, or reading unconditionally),
        // this would be 2 -- and if the wrong branch won, it would be
        // [2, 1, 0] instead.
        assert_eq!(angles.len(), 1);
        assert_eq!(angles[0], [0, 1, 2]);
    }

    #[test]
    fn test_no_coordinates_are_ever_invented() {
        let (_, mol) = &parse_top(UREA_LIKE).expect("valid TOP")[0];
        assert!(!mol.has_coords3());
        assert!(!mol.has_coords());
    }

    #[test]
    fn test_a_mass_matching_no_element_is_a_clear_error() {
        let text = "[ moleculetype ]\nM 3\n[ atoms ]\n1 X 1 RES A 1 0.0 999.0\n";
        let err = parse_top(text).unwrap_err();
        assert!(matches!(err, TopError::InvalidElement(_)), "{err}");
    }

    #[test]
    fn test_an_atom_line_missing_mass_is_a_clear_error_not_a_guess() {
        let text = "[ moleculetype ]\nM 3\n[ atoms ]\n1 X 1 RES A 1 0.0\n";
        let err = parse_top(text).unwrap_err();
        assert!(matches!(err, TopError::InvalidAtomLine(_)), "{err}");
    }

    #[test]
    fn test_empty_input_is_a_clear_error_not_a_panic() {
        let err = parse_top("").unwrap_err();
        assert!(matches!(err, TopError::NoAtoms), "{err}");
    }

    #[test]
    fn test_is_gromacs_top_distinguishes_from_prmtop() {
        assert!(is_gromacs_top(UREA_LIKE));
        assert!(!is_gromacs_top(
            "%VERSION  VERSION_STAMP = V0001.000  DATE = 00/00/00\n%FLAG TITLE\n"
        ));
    }
}
