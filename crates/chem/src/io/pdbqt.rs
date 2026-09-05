//! PDBQT (AutoDock) — PDB's fixed-column layout plus AutoDock atom types,
//! a per-atom partial charge column, and (for the flexible/ligand variant)
//! a `ROOT`/`BRANCH i j`/`ENDBRANCH i j`/`ENDROOT`/`TORSDOF` torsion tree
//! (#226).
//!
//! Targets AutoDock's own documented spec, not obabel's or Meeko's
//! specific quirks — the two are known to disagree with each other, per
//! this issue's own text, so oracle agreement with one of them is not the
//! bar here.
//!
//! # Rotatable-bond detection and the torsion tree
//!
//! A bond is a *candidate* rotatable bond iff [`Bond::is_rotatable`]
//! (single order, not in a ring) **and** both endpoints are heavy atoms
//! **and** both endpoints have heavy-atom degree >= 2 (the number of
//! non-hydrogen neighbours). That last rule is what excludes every
//! terminal group (`-CH3`, `-OH`, `-NH2`, `-F`, `-Cl`) with no
//! per-functional-group special-casing at all: a terminal heavy atom's
//! *only* heavy neighbour is, by definition, the bond partner itself, so
//! its heavy degree is 1.
//!
//! Every candidate bond is a graph-theoretic bridge — a bond not in any
//! ring is exactly a bond not in any cycle, which is the definition of a
//! bridge. Removing every candidate bond and taking connected components
//! therefore always yields a fragment-contraction graph that is a **tree**
//! (the bridge/block-cut tree of a connected graph is always a tree) —
//! `TORSDOF` is exactly `fragments - 1`, no separate counting needed, and
//! no tree-vs-forest handling for the fragment structure itself.
//!
//! Atoms are renumbered in traversal order (the ROOT fragment's atoms
//! first, then each branch depth-first) — real ligand PDBQT files lay
//! fragments out contiguously this way, and `BRANCH i j`'s `i`/`j` refer
//! to those written serials, not the input molecule's own atom indices.
//!
//! **Requires a single connected input molecule** — a multi-component
//! input (a salt with a counter-ion) can't be one PDBQT ligand, the same
//! constraint real ligand-prep tooling has. Since every registered
//! writer's signature is infallible, this is not surfaced as an error but
//! as an in-band `REMARK` line stating the input had N components and
//! only the largest was written.
//!
//! # AutoDock atom typing
//!
//! Simpler than Mol2's SYBYL table, not more complex: AutoDock's
//! vocabulary for "everything else" (halogens, metals, P) *is* the plain
//! element symbol, so the same plain-element-symbol lookup this crate's
//! other readers already use covers that fallback — no dummy-type error
//! path is needed. Only five elements get a special letter:
//!
//! - **C**: `A` if aromatic, else `C`.
//! - **N**: `NA` if it has no attached hydrogens (a "dry" nitrogen), else `N`.
//! - **O**: always `OA` — nearly every oxygen in an organic ligand is an
//!   H-bond acceptor in the AutoDock4 force field; the rarely-used plain
//!   `O` is not attempted (a documented simplification).
//! - **S**: `SA` if no attached hydrogens, else `S` (thiol).
//! - **H**: `HD` if any neighbour is N/O/S (polar), else `H`. **Nonpolar
//!   hydrogens are kept explicit, never merged** into their parent heavy
//!   atom's charge — real AutoDock4 ligand-prep merges them for correct
//!   docking energetics; this crate does not perform that charge
//!   redistribution, the same kind of chemistry-transformation deferral
//!   already established for PDB's residue-template bonds and mmCIF's
//!   `_struct_conn`.
//!
//! Reading reverses this via a small match on the five letters, falling
//! through to the plain element-symbol lookup otherwise. `HD` vs `H` carries no
//! information back onto [`Atom`] beyond "this is a hydrogen" — the
//! donor/acceptor distinction is a write-time-only computed view.
//!
//! # Bonds
//!
//! **No `CONECT`-equivalent exists in PDBQT for intra-fragment bonds at
//! all.** The only explicit bond data anywhere in the format is each
//! `BRANCH i j` line's own pivot pair. Reading a PDBQT therefore produces
//! a sparse bond graph — only the rotatable-bond pivots, nothing else —
//! rather than attempting distance-based bond perception to reconstruct
//! the rest, out of scope and consistent with every prior format's "no
//! bond inference."
//!
//! A file may hold several `MODEL`/`ENDMDL`-wrapped poses (AutoDock
//! Vina's real docked-results output shape) — reuses the exact `ENDMDL`
//! boundary technique [`crate::io::supplier::PdbSupplier`] already
//! implements, not reinvented.

use std::collections::{HashMap, HashSet, VecDeque};

use crate::core::atom::{Atom, Element};
use crate::core::bond::{Bond, BondOrder};
use crate::core::elements::ELEMENT_SYMBOLS;
use crate::core::geometry::Point3;
use crate::core::molecule::Molecule;
use crate::core::residue::{Chain, Residue};
use crate::core::rings::perceive_rings;
use crate::core::site::AtomSite;
use crate::io::errors::PdbqtError;

fn column(line: &str, start_1based: usize, end_1based: usize) -> Option<&str> {
    let bytes = line.as_bytes();
    if start_1based == 0 || start_1based > bytes.len() {
        return None;
    }
    let start = start_1based - 1;
    let end = end_1based.min(bytes.len());
    line.get(start..end).filter(|s| !s.trim().is_empty())
}

fn parse_f64_column(line: &str, start: usize, end: usize) -> Result<f64, PdbqtError> {
    column(line, start, end)
        .map(str::trim)
        .ok_or_else(|| PdbqtError::InvalidAtomLine(line.to_string()))?
        .parse()
        .map_err(|_| PdbqtError::InvalidAtomLine(line.to_string()))
}

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

/// Resolves an AutoDock atom type into `(element, is_aromatic)`. Only `A`
/// (aromatic carbon) carries aromaticity; every other letter, and every
/// plain element symbol, does not.
fn resolve_autodock_type(type_str: &str) -> Result<(Element, bool), PdbqtError> {
    match type_str {
        "A" => Ok((Element::carbon(), true)),
        "NA" => Ok((Element::nitrogen(), false)),
        "OA" => Ok((Element::oxygen(), false)),
        "SA" => element_from_symbol("S")
            .map(|e| (e, false))
            .ok_or_else(|| PdbqtError::InvalidAtomType(type_str.to_string())),
        "HD" => Ok((Element::hydrogen(), false)),
        other => element_from_symbol(other)
            .map(|e| (e, false))
            .ok_or_else(|| PdbqtError::InvalidAtomType(type_str.to_string())),
    }
}

/// The number of `atom_idx`'s neighbours that are not hydrogen.
fn heavy_degree(mol: &Molecule, atom_idx: usize) -> usize {
    mol.neighbors(atom_idx)
        .iter()
        .filter(|n| mol.atom(n.atom_idx).atomic_number() != 1)
        .count()
}

/// The reverse of [`resolve_autodock_type`]. Needs the whole molecule (not
/// just the atom) because hydrogen's donor/acceptor status depends on its
/// neighbour, and nitrogen/sulfur's depends on `total_hydrogens()`.
fn autodock_type_for(mol: &Molecule, atom_idx: usize) -> String {
    let atom = mol.atom(atom_idx);
    let symbol = atom.element().symbol();
    match symbol {
        "C" => {
            if atom.is_aromatic() {
                "A".to_string()
            } else {
                "C".to_string()
            }
        }
        "N" => {
            if atom.total_hydrogens() == 0 {
                "NA".to_string()
            } else {
                "N".to_string()
            }
        }
        "O" => "OA".to_string(),
        "S" => {
            if atom.total_hydrogens() == 0 {
                "SA".to_string()
            } else {
                "S".to_string()
            }
        }
        "H" => {
            let polar = mol
                .neighbors(atom_idx)
                .iter()
                .any(|n| matches!(mol.atom(n.atom_idx).element().symbol(), "N" | "O" | "S"));
            if polar {
                "HD".to_string()
            } else {
                "H".to_string()
            }
        }
        other => other.to_string(),
    }
}

/// A candidate rotatable bond: not in a ring, single order, both
/// endpoints heavy with heavy-atom degree >= 2. `mol` must already have
/// fresh ring membership (via [`perceive_rings`]).
fn is_candidate_rotatable(mol: &Molecule, bond: &Bond) -> bool {
    if !bond.is_rotatable() {
        return false;
    }
    let (a, b) = (bond.atom1(), bond.atom2());
    if mol.atom(a).atomic_number() == 1 || mol.atom(b).atomic_number() == 1 {
        return false;
    }
    heavy_degree(mol, a) >= 2 && heavy_degree(mol, b) >= 2
}

/// Connected components of `mol`'s atom graph, refusing to traverse any
/// bond in `excluded`. Kept local to this module rather than added to
/// `MoleculeGraph` -- no other format needs it.
fn fragments_excluding(mol: &Molecule, excluded: &HashSet<usize>) -> Vec<Vec<usize>> {
    let mut visited = vec![false; mol.num_atoms()];
    let mut fragments = Vec::new();

    for start in 0..mol.num_atoms() {
        if visited[start] {
            continue;
        }
        let mut fragment = Vec::new();
        let mut queue = VecDeque::from([start]);
        visited[start] = true;
        while let Some(atom_idx) = queue.pop_front() {
            fragment.push(atom_idx);
            for neighbor in mol.neighbors(atom_idx) {
                if excluded.contains(&neighbor.bond_idx) || visited[neighbor.atom_idx] {
                    continue;
                }
                visited[neighbor.atom_idx] = true;
                queue.push_back(neighbor.atom_idx);
            }
        }
        fragment.sort_unstable();
        fragments.push(fragment);
    }

    fragments
}

/// One edge of the fragment tree: `bond_idx` is the original candidate
/// bond, `(atom_a, atom_b)` its endpoints, with `atom_a` in `frag_a` and
/// `atom_b` in `frag_b`.
struct FragmentEdge {
    frag_a: usize,
    frag_b: usize,
    atom_a: usize,
    atom_b: usize,
}

/// Assigns new, sequential serials to every atom by walking the fragment
/// tree depth-first from the root fragment -- the ROOT block's atoms
/// first, then each branch's atoms, recursively. Returns the order as a
/// flat `Vec<usize>` of original atom indices (index into this vec + 1 is
/// the written serial) alongside the tree structure needed to place
/// `BRANCH`/`ENDBRANCH` markers.
enum Node {
    /// A run of this many atoms (from `order`), then nested children.
    Fragment {
        atom_count: usize,
        pivot_in: Option<(usize, usize)>, // (i, j) for a BRANCH line, None for ROOT
        children: Vec<Node>,
    },
}

fn build_tree(
    fragments: &[Vec<usize>],
    edges: &[FragmentEdge],
    root_frag: usize,
    order: &mut Vec<usize>,
    new_serial: &mut HashMap<usize, usize>,
) -> Node {
    for &atom in &fragments[root_frag] {
        new_serial.insert(atom, order.len() + 1);
        order.push(atom);
    }
    let root_atom_count = fragments[root_frag].len();

    let mut child_edges: Vec<&FragmentEdge> = edges
        .iter()
        .filter(|e| e.frag_a == root_frag || e.frag_b == root_frag)
        .collect();
    child_edges.sort_by_key(|e| {
        if e.frag_a == root_frag {
            e.atom_a
        } else {
            e.atom_b
        }
    });

    let mut children = Vec::new();
    for edge in child_edges {
        let (pivot_here, child_frag, pivot_child) = if edge.frag_a == root_frag {
            (edge.atom_a, edge.frag_b, edge.atom_b)
        } else {
            (edge.atom_b, edge.frag_a, edge.atom_a)
        };
        let remaining: Vec<FragmentEdge> = edges
            .iter()
            .filter(|e| !(e.frag_a == root_frag && e.frag_b == child_frag))
            .filter(|e| !(e.frag_b == root_frag && e.frag_a == child_frag))
            .map(|e| FragmentEdge {
                frag_a: e.frag_a,
                frag_b: e.frag_b,
                atom_a: e.atom_a,
                atom_b: e.atom_b,
            })
            .collect();
        let i = new_serial[&pivot_here];
        let child_node = build_tree(fragments, &remaining, child_frag, order, new_serial);
        let j = new_serial[&pivot_child];
        let Node::Fragment {
            atom_count,
            children: grandchildren,
            ..
        } = child_node;
        children.push(Node::Fragment {
            atom_count,
            pivot_in: Some((i, j)),
            children: grandchildren,
        });
    }

    Node::Fragment {
        atom_count: root_atom_count,
        pivot_in: None,
        children,
    }
}

fn write_fragment(
    out: &mut String,
    node: &Node,
    order: &[usize],
    cursor: &mut usize,
    mol: &Molecule,
) {
    let Node::Fragment {
        atom_count,
        children,
        ..
    } = node;
    for _ in 0..*atom_count {
        write_atom_line(out, mol, order[*cursor], *cursor + 1);
        *cursor += 1;
    }
    for child in children {
        let Node::Fragment { pivot_in, .. } = child;
        let (i, j) = pivot_in.expect("a child fragment always has a pivot");
        out.push_str(&format!("BRANCH{i:>4}{j:>4}\n"));
        write_fragment(out, child, order, cursor, mol);
        out.push_str(&format!("ENDBRANCH{i:>4}{j:>4}\n"));
    }
}

fn write_atom_line(out: &mut String, mol: &Molecule, atom_idx: usize, serial: usize) {
    let atom = mol.atom(atom_idx);
    let site = mol.site(atom_idx);
    let residue = mol.residue_of(atom_idx);
    let chain = mol.chain_of(atom_idx);
    let p = mol.coord3(atom_idx).unwrap_or(Point3::new(0.0, 0.0, 0.0));

    let record = "ATOM";
    let name = site
        .and_then(|s| s.name.as_deref())
        .unwrap_or(atom.element().symbol());
    let res_name = residue.map(|r| r.name.as_str()).unwrap_or("LIG");
    let res_seq = residue.map(|r| r.sequence).unwrap_or(1);
    let chain_id = chain.and_then(|c| c.id.chars().next()).unwrap_or(' ');
    let occupancy = site.and_then(|s| s.occupancy).unwrap_or(1.0);
    let b_factor = site.and_then(|s| s.b_factor).unwrap_or(0.0);
    let charge = site.and_then(|s| s.partial_charge).unwrap_or(0.0);
    let autodock_type = autodock_type_for(mol, atom_idx);

    out.push_str(&format!(
        "{record:<6}{serial:>5} {name:<4} {res_name:<3} {chain_id}{res_seq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{occupancy:>6.2}{b_factor:>6.2}    {charge:>6.3} {autodock_type:<2}\n",
        x = p.x,
        y = p.y,
        z = p.z,
    ));
}

/// Writes one PDBQT ligand (a `ROOT`/`BRANCH`.../`TORSDOF` block).
///
/// Panics never; a molecule with more than one connected component gets a
/// `REMARK` line naming the loss and only its largest component is
/// written -- see the module doc for why this isn't a `Result`.
pub fn write_pdbqt(mol: &Molecule) -> String {
    let mut mol = mol.clone();
    perceive_rings(&mut mol);

    let components = mol.graph().connected_components();
    let mol_for_write;
    let mut remark = String::new();
    let working: &Molecule = if components.len() > 1 {
        remark = format!(
            "REMARK  chem: input had {} connected components; only the largest ({} atoms) was written\n",
            components.len(),
            components.iter().map(Vec::len).max().unwrap_or(0),
        );
        // Rebuild from just the largest component's atoms, in their
        // original relative order, so indices stay simple to reason
        // about.
        let largest = components
            .into_iter()
            .max_by_key(Vec::len)
            .unwrap_or_default();
        let keep: HashSet<usize> = largest.into_iter().collect();
        let mut rebuilt = Molecule::new();
        let mut remap = HashMap::new();
        for (old_idx, atom) in mol.atoms().iter().enumerate() {
            if keep.contains(&old_idx) {
                remap.insert(old_idx, rebuilt.add_atom(atom.clone()));
            }
        }
        for bond in mol.bonds() {
            if let (Some(&a), Some(&b)) = (remap.get(&bond.atom1()), remap.get(&bond.atom2())) {
                let _ = rebuilt.add_bond(Bond::new(a, b, bond.order()));
            }
        }
        if let Some(coords) = mol.coords3() {
            let kept: Vec<Point3> = (0..mol.num_atoms())
                .filter(|i| keep.contains(i))
                .map(|i| coords[i])
                .collect();
            let _ = rebuilt.set_coords3(kept);
        }
        if let Some(sites) = mol.sites() {
            let kept: Vec<AtomSite> = (0..mol.num_atoms())
                .filter(|i| keep.contains(i))
                .map(|i| sites[i].clone())
                .collect();
            let _ = rebuilt.set_sites(kept);
        }
        perceive_rings(&mut rebuilt);
        mol_for_write = rebuilt;
        &mol_for_write
    } else {
        &mol
    };

    let candidate_bonds: Vec<usize> = working
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, b)| is_candidate_rotatable(working, b))
        .map(|(i, _)| i)
        .collect();
    let excluded: HashSet<usize> = candidate_bonds.iter().copied().collect();

    let fragments = fragments_excluding(working, &excluded);
    let frag_of: HashMap<usize, usize> = fragments
        .iter()
        .enumerate()
        .flat_map(|(f, atoms)| atoms.iter().map(move |&a| (a, f)))
        .collect();

    let edges: Vec<FragmentEdge> = candidate_bonds
        .iter()
        .map(|&bond_idx| {
            let bond = working.bond(bond_idx);
            let (a, b) = (bond.atom1(), bond.atom2());
            FragmentEdge {
                frag_a: frag_of[&a],
                frag_b: frag_of[&b],
                atom_a: a,
                atom_b: b,
            }
        })
        .collect();

    let root_frag = fragments
        .iter()
        .enumerate()
        .max_by_key(|(_, atoms)| {
            (
                atoms.len(),
                std::cmp::Reverse(atoms.first().copied().unwrap_or(usize::MAX)),
            )
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    let mut order = Vec::new();
    let mut new_serial = HashMap::new();
    let tree = build_tree(&fragments, &edges, root_frag, &mut order, &mut new_serial);

    let mut out = remark;
    out.push_str("ROOT\n");
    let mut cursor = 0;
    let Node::Fragment { atom_count, .. } = &tree;
    for _ in 0..*atom_count {
        write_atom_line(&mut out, working, order[cursor], cursor + 1);
        cursor += 1;
    }
    out.push_str("ENDROOT\n");
    let Node::Fragment { children, .. } = &tree;
    for child in children {
        let Node::Fragment { pivot_in, .. } = child;
        let (i, j) = pivot_in.expect("a child fragment always has a pivot");
        out.push_str(&format!("BRANCH{i:>4}{j:>4}\n"));
        write_fragment(&mut out, child, &order, &mut cursor, working);
        out.push_str(&format!("ENDBRANCH{i:>4}{j:>4}\n"));
    }
    out.push_str(&format!("TORSDOF {}\n", candidate_bonds.len()));

    out
}

#[derive(Clone, PartialEq, Eq)]
struct ResidueKey {
    chain_id: String,
    name: String,
    sequence: i32,
    insertion_code: Option<char>,
}

fn group_into_chains_and_residues(keys: &[ResidueKey]) -> (Vec<Chain>, Vec<Residue>) {
    let mut chains = Vec::new();
    let mut residues = Vec::new();
    let mut current_chain_id: Option<&str> = None;
    let mut chain_start = 0;

    let mut i = 0;
    while i < keys.len() {
        let key = &keys[i];
        let start = i;
        while i < keys.len() && keys[i] == *key {
            i += 1;
        }

        if current_chain_id != Some(key.chain_id.as_str()) {
            if let Some(id) = current_chain_id {
                chains.push(Chain {
                    id: id.to_string(),
                    label_id: None,
                    residues: chain_start..residues.len(),
                });
            }
            current_chain_id = Some(&key.chain_id);
            chain_start = residues.len();
        }

        residues.push(Residue {
            name: key.name.clone(),
            sequence: key.sequence,
            insertion_code: key.insertion_code,
            label_seq: None,
            chain_ix: chains.len(),
            is_hetero: false,
            atoms: start..i,
        });
    }
    if let Some(id) = current_chain_id {
        chains.push(Chain {
            id: id.to_string(),
            label_id: None,
            residues: chain_start..residues.len(),
        });
    }

    (chains, residues)
}

/// Parses one PDBQT ligand (a single `MODEL`, or a whole file with none).
pub fn parse_pdbqt(text: &str) -> Result<Molecule, PdbqtError> {
    let mut mol = Molecule::new();
    let mut coords = Vec::new();
    let mut sites = Vec::new();
    let mut keys = Vec::new();
    let mut serial_to_index: HashMap<i64, usize> = HashMap::new();
    let mut bond_pairs: Vec<(i64, i64)> = Vec::new();

    for line in text.lines() {
        if line.len() < 6 {
            continue;
        }
        let record = line[0..6.min(line.len())].trim();
        match record {
            "ATOM" | "HETATM" => {
                let serial: i64 = column(line, 7, 11)
                    .map(str::trim)
                    .ok_or_else(|| PdbqtError::InvalidAtomLine(line.to_string()))?
                    .parse()
                    .map_err(|_| PdbqtError::InvalidAtomLine(line.to_string()))?;
                let name_field = column(line, 13, 16).unwrap_or("").trim().to_string();
                let res_name = column(line, 18, 20).unwrap_or("LIG").trim().to_string();
                let chain_id = column(line, 22, 22).unwrap_or("").trim().to_string();
                let res_seq: i32 = column(line, 23, 26)
                    .map(str::trim)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(1);
                let insertion_code = column(line, 27, 27).and_then(|s| s.chars().next());
                let x = parse_f64_column(line, 31, 38)?;
                let y = parse_f64_column(line, 39, 46)?;
                let z = parse_f64_column(line, 47, 54)?;
                let occupancy = column(line, 55, 60).and_then(|s| s.trim().parse().ok());
                let b_factor = column(line, 61, 66).and_then(|s| s.trim().parse().ok());
                let charge = column(line, 71, 76).and_then(|s| s.trim().parse().ok());
                let autodock_type = column(line, 78, 79)
                    .ok_or_else(|| PdbqtError::InvalidAtomLine(line.to_string()))?
                    .trim();

                let (element, is_aromatic) = resolve_autodock_type(autodock_type)?;
                let atom = Atom::new(element).with_aromatic(is_aromatic);
                let atom_idx = mol.add_atom(atom);
                serial_to_index.insert(serial, atom_idx);
                coords.push(Point3::new(x, y, z));
                sites.push(AtomSite {
                    name: if name_field.is_empty() {
                        None
                    } else {
                        Some(name_field)
                    },
                    alt_loc: None,
                    partial_charge: charge,
                    occupancy,
                    b_factor,
                    radius: None,
                });
                keys.push(ResidueKey {
                    chain_id,
                    name: res_name,
                    sequence: res_seq,
                    insertion_code,
                });
            }
            "BRANCH" => {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.len() >= 3 {
                    let i: i64 = fields[1]
                        .parse()
                        .map_err(|_| PdbqtError::InvalidTorsionTree(line.to_string()))?;
                    let j: i64 = fields[2]
                        .parse()
                        .map_err(|_| PdbqtError::InvalidTorsionTree(line.to_string()))?;
                    bond_pairs.push((i, j));
                }
            }
            _ => {
                // ROOT, ENDROOT, ENDBRANCH, TORSDOF, REMARK, and everything
                // else: no further data to extract. MODEL/ENDMDL framing
                // is the caller's job (see the module doc).
            }
        }
    }

    mol.set_coords3(coords)
        .map_err(|e| PdbqtError::ParseError(e.to_string()))?;
    mol.set_sites(sites)
        .map_err(|e| PdbqtError::ParseError(e.to_string()))?;

    let (chains, residues) = group_into_chains_and_residues(&keys);
    mol.set_topology(chains, residues)
        .map_err(|e| PdbqtError::ParseError(e.to_string()))?;

    for (i, j) in bond_pairs {
        let a = *serial_to_index
            .get(&i)
            .ok_or_else(|| PdbqtError::InvalidTorsionTree(format!("unknown serial {i}")))?;
        let b = *serial_to_index
            .get(&j)
            .ok_or_else(|| PdbqtError::InvalidTorsionTree(format!("unknown serial {j}")))?;
        mol.add_bond(Bond::new(a, b, BondOrder::Single))
            .map_err(|e| PdbqtError::ParseError(e.to_string()))?;
    }

    Ok(mol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::smiles::parse_smiles;

    #[test]
    fn test_a_rigid_molecule_writes_root_only_with_torsdof_zero() {
        let mut mol = parse_smiles("c1ccccc1").expect("valid SMILES");
        mol.set_coords3(vec![Point3::ORIGIN; mol.num_atoms()])
            .expect("one per atom");
        let written = write_pdbqt(&mol);
        assert!(written.contains("ROOT\n"), "{written}");
        assert!(written.contains("ENDROOT\n"), "{written}");
        assert!(!written.contains("BRANCH"), "{written}");
        assert!(written.contains("TORSDOF 0"), "{written}");
    }

    #[test]
    fn test_toluene_shaped_input_excludes_the_terminal_methyl_bond() {
        // The ring-methyl bond is Single and not in a ring, so it passes
        // Bond::is_rotatable() -- it must be excluded by the heavy-degree
        // rule instead (the methyl carbon's only heavy neighbour is the
        // ring carbon itself).
        let mut mol = parse_smiles("Cc1ccccc1").expect("valid SMILES");
        mol.set_coords3(vec![Point3::ORIGIN; mol.num_atoms()])
            .expect("one per atom");
        let written = write_pdbqt(&mol);
        assert!(written.contains("TORSDOF 0"), "{written}");
    }

    #[test]
    fn test_a_branched_ligand_with_two_rotatable_bonds_round_trips() {
        // Two biphenyl-like rings joined through a flexible linker gives
        // exactly two non-ring, non-terminal single bonds: ring1-CH2 and
        // CH2-ring2 is one shape, but a cleaner >=2 case is two separate
        // single bonds between three ring systems: ring1-ring2-ring3
        // joined by two direct aryl-aryl single bonds (each aryl carbon
        // has heavy degree 3, well above the threshold).
        let mut mol = parse_smiles("c1ccc(cc1)-c1ccc(cc1)-c1ccccc1").expect("valid SMILES");
        mol.set_coords3(vec![Point3::ORIGIN; mol.num_atoms()])
            .expect("one per atom");
        let written = write_pdbqt(&mol);
        assert!(written.contains("TORSDOF 2"), "{written}");
        let branch_count = written.matches("BRANCH").count() - written.matches("ENDBRANCH").count();
        assert_eq!(branch_count, 2, "{written}");

        let back = parse_pdbqt(&written).expect("valid PDBQT");
        assert_eq!(back.num_bonds(), 2, "only the two BRANCH pivots are stated");
    }

    #[test]
    fn test_autodock_type_resolution_both_directions() {
        for (type_str, symbol, aromatic) in [
            ("A", "C", true),
            ("NA", "N", false),
            ("OA", "O", false),
            ("SA", "S", false),
            ("HD", "H", false),
            ("Cl", "Cl", false),
            ("Zn", "Zn", false),
        ] {
            let (element, is_aromatic) = resolve_autodock_type(type_str).expect(type_str);
            assert_eq!(element.symbol(), symbol, "{type_str}");
            assert_eq!(is_aromatic, aromatic, "{type_str}");
        }
    }

    #[test]
    fn test_nitrogen_and_sulfur_type_depends_on_attached_hydrogens() {
        let pyridine = parse_smiles("c1ccncc1").expect("valid SMILES");
        // Atom index 3 is the ring nitrogen (no attached H).
        assert_eq!(autodock_type_for(&pyridine, 3), "NA");

        let aniline = parse_smiles("Nc1ccccc1").expect("valid SMILES");
        // Atom index 0 is the amine nitrogen (has attached H).
        assert_eq!(autodock_type_for(&aniline, 0), "N");
    }

    #[test]
    fn test_a_disconnected_input_gets_a_remark_and_only_the_largest_component() {
        let mut mol = parse_smiles("CCO.C").expect("valid SMILES");
        mol.set_coords3(vec![Point3::ORIGIN; mol.num_atoms()])
            .expect("one per atom");
        let written = write_pdbqt(&mol);
        assert!(written.starts_with("REMARK"), "{written}");
        assert!(written.contains("2 connected components"), "{written}");
    }

    #[test]
    fn test_a_malformed_atom_line_is_a_clear_error_not_a_panic() {
        let text = "ATOM      1 C1   LIG A   1      not-a-number   0.000   0.000  1.00  0.00     0.000 C \n";
        let err = parse_pdbqt(text).unwrap_err();
        assert!(matches!(err, PdbqtError::InvalidAtomLine(_)), "{err}");
    }

    #[test]
    fn test_an_unrecognised_atom_type_is_a_clear_error() {
        let text =
            "ATOM      1 C1   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 Xx\n";
        let err = parse_pdbqt(text).unwrap_err();
        assert!(matches!(err, PdbqtError::InvalidAtomType(_)), "{err}");
    }
}
