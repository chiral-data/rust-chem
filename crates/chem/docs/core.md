# ChemCore

Core molecular data structures for cheminformatics applications.

## Overview

ChemCore provides the fundamental building blocks for representing molecules:

- **Atoms**: Elements with properties (charge, isotope, aromaticity)
- **Bonds**: Connections with order and stereochemistry
- **Molecules**: Complete molecular structures
- **Graph**: Efficient connectivity queries
- **Side tables**: Per-atom data a file supplied — 2D layout, 3D conformer, atom sites
- **Topology**: Chains and residues, for the formats organised by them

## Quick Start

```rust
use chem::core::prelude::*;

fn main() {
    // Create a molecule
    let mut mol = Molecule::new();

    // Add atoms
    let c1 = mol.add_atom(Atom::new(Element::carbon()));
    let c2 = mol.add_atom(Atom::new(Element::carbon()));

    // Add bond
    mol.add_bond(Bond::new(c1, c2, BondOrder::Single)).unwrap();

    // Calculate hydrogens
    mol.calculate_implicit_hydrogens();

    // Get properties
    println!("Formula: {}", mol.formula());  // C2H6
    println!("Weight: {:.2} amu", mol.molecular_weight());
}
```

## Molecule Representation Examples

### Example 1: Water (H₂O)

```rust
use chem::core::prelude::*;

let mut water = Molecule::new();

// Add oxygen atom
let o = water.add_atom(Atom::new(Element::oxygen()));

// Set 2 implicit hydrogens (not shown as explicit atoms)
water.atom_mut(o).set_implicit_hydrogens(2);

println!("{}", water.formula());  // H2O
println!("{:.2} amu", water.molecular_weight());  // 18.02 amu
```

**Internal Representation:**

```
Atoms: [O]
Bonds: []
Implicit H: 2 on oxygen
Graph: Single node (no edges)
```

### Example 2: Methane (CH₄)

```rust
let mut methane = Molecule::new();

// Add carbon
let c = methane.add_atom(Atom::new(Element::carbon()));

// Calculate implicit hydrogens automatically
methane.calculate_implicit_hydrogens();

println!("{}", methane.formula());  // CH4
```

**Internal Representation:**

```
Atoms: [C]
Bonds: []
Implicit H: 4 (calculated from carbon's valence)
Graph: Single node
```

### Example 3: Ethane (C₂H₆)

```rust
let mut ethane = Molecule::new();

// Add two carbons
let c1 = ethane.add_atom(Atom::new(Element::carbon()));
let c2 = ethane.add_atom(Atom::new(Element::carbon()));

// Connect with single bond
ethane.add_bond(Bond::new(c1, c2, BondOrder::Single)).unwrap();

// Calculate hydrogens
ethane.calculate_implicit_hydrogens();

println!("{}", ethane.formula());  // C2H6

// Query connectivity
println!("Neighbors of C1: {}", ethane.degree(c1));  // 1 (just C2)
```

**Internal Representation:**

```
Atoms: [C, C]
Bonds: [C0-C1 (single)]
Implicit H: 3 on each carbon
Graph:
  0 -- 1
  (adjacency list: 0:[1], 1:[0])
```

### Example 4: Ethene (C₂H₄) - Double Bond

```rust
let mut ethene = Molecule::new();

let c1 = ethene.add_atom(Atom::new(Element::carbon()));
let c2 = ethene.add_atom(Atom::new(Element::carbon()));

// Double bond
ethene.add_bond(Bond::new(c1, c2, BondOrder::Double)).unwrap();
ethene.calculate_implicit_hydrogens();

println!("{}", ethene.formula());  // C2H4
```

**Internal Representation:**

```
Atoms: [C, C]
Bonds: [C0=C1 (double)]
Implicit H: 2 on each carbon (valence 4 - 2 from double bond)
Graph:
  0 == 1
```

### Example 5: Propane (C₃H₈) - Linear Chain

```rust
let mut propane = Molecule::new();

// Add three carbons
let c1 = propane.add_atom(Atom::new(Element::carbon()));
let c2 = propane.add_atom(Atom::new(Element::carbon()));
let c3 = propane.add_atom(Atom::new(Element::carbon()));

// Connect in chain: C-C-C
propane.add_bond(Bond::new(c1, c2, BondOrder::Single)).unwrap();
propane.add_bond(Bond::new(c2, c3, BondOrder::Single)).unwrap();

propane.calculate_implicit_hydrogens();

println!("{}", propane.formula());  // C3H8

// Query neighbors
for neighbor in propane.neighbors(c2) {
    println!("C2 connected to: {}", neighbor.atom_idx);
}
// Output: C2 connected to: 0
//         C2 connected to: 2
```

**Internal Representation:**

```
Atoms: [C, C, C]
Bonds: [C0-C1, C1-C2]
Implicit H: [3, 2, 3] (terminal carbons have 3, middle has 2)
Graph:
  0 -- 1 -- 2
  (adjacency list: 0:[1], 1:[0,2], 2:[1])
```

### Example 6: Benzene (C₆H₆) - Aromatic Ring

```rust
let mut benzene = Molecule::new();

// Add 6 aromatic carbons
let carbons: Vec<usize> = (0..6)
    .map(|_| benzene.add_atom(
        Atom::new(Element::carbon()).with_aromatic(true)
    ))
    .collect();

// Connect in a ring with aromatic bonds
for i in 0..6 {
    let next = (i + 1) % 6;
    benzene.add_bond(
        Bond::new(carbons[i], carbons[next], BondOrder::Aromatic)
            .with_aromatic(true)
    ).unwrap();
}

benzene.calculate_implicit_hydrogens();

println!("{}", benzene.formula());  // C6H6

// Check ring structure
let neighbors = benzene.graph().neighbors_within_radius(carbons[0], 3);
println!("Atoms within 3 bonds: {}", neighbors.len());  // 5 (all others)
```

**Internal Representation:**

```
Atoms: [C, C, C, C, C, C] (all aromatic)
Bonds: [C0:C1, C1:C2, C2:C3, C3:C4, C4:C5, C5:C0] (aromatic)
Implicit H: 1 on each carbon
Graph (ring):
    0 --- 1
   /       \
  5         2
   \       /
    4 --- 3
  (each atom has degree 2)
```

### Example 7: Ammonium Ion (NH₄⁺) - Charged Molecule

```rust
let mut ammonium = Molecule::new();

// Add nitrogen with +1 charge
let n = ammonium.add_atom(
    Atom::new(Element::nitrogen()).with_charge(1)
);

// Calculate hydrogens (considers charge)
ammonium.calculate_implicit_hydrogens();

println!("{}", ammonium.formula());  // H4N
println!("Charge: {}", ammonium.atom(n).formal_charge());  // +1
```

**Internal Representation:**

```
Atoms: [N+]
Bonds: []
Implicit H: 4 (nitrogen valence 3 + 1 for positive charge)
Charge: +1
```

### Example 8: Graph Traversal - Finding Connected Atoms

```rust
let mut mol = Molecule::new();

// Create a branched molecule
let c1 = mol.add_atom(Atom::new(Element::carbon()));
let c2 = mol.add_atom(Atom::new(Element::carbon()));
let c3 = mol.add_atom(Atom::new(Element::carbon()));
let c4 = mol.add_atom(Atom::new(Element::carbon()));

mol.add_bond(Bond::new(c1, c2, BondOrder::Single)).unwrap();
mol.add_bond(Bond::new(c1, c3, BondOrder::Single)).unwrap();
mol.add_bond(Bond::new(c1, c4, BondOrder::Single)).unwrap();

// DFS from c1
let visited = mol.graph().dfs(c1);
println!("DFS order: {:?}", visited);  // [0, 1, 2, 3]

// BFS from c1
let visited = mol.graph().bfs(c1);
println!("BFS order: {:?}", visited);  // [0, 1, 2, 3]

// Neighbors within radius
let nearby = mol.graph().neighbors_within_radius(c1, 1);
println!("Direct neighbors: {}", nearby.len());  // 3 (c2, c3, c4)
```

**Internal Representation:**

```
Atoms: [C, C, C, C]
Bonds: [C0-C1, C0-C2, C0-C3]
Graph (star):
    1
    |
  2-0-3
  (center has degree 3, others degree 1)
```

## Key Concepts

### 1. Atoms Store Properties

- Element (atomic number)
- Charge (formal charge)
- Isotope (optional)
- Aromaticity flag
- Implicit hydrogens (calculated or set)

Note what is **not** here: no coordinates, no partial charge, no occupancy, no
temperature factor, no atom name. Those live on the molecule as side tables, for
two reasons.

The mechanical one is that `Atom` derives `Eq`, and every one of them is a float.
The honest one is that they are a different kind of fact. An element is a
property of a species; a B-factor is a property of one *observation* of it. Two
files describing the same molecule will agree on the first and disagree on the
second.

### 2. Bonds Store Connectivity

- Two atom indices
- Bond order (single, double, triple, aromatic)
- Stereochemistry (E/Z)

### 3. Molecule Manages Everything

- Collection of atoms
- Collection of bonds
- Graph for fast queries
- Properties (metadata)
- Up to three per-atom side tables, each independently present or absent
- Chain and residue topology, for structural formats

### 3a. The Side Tables

Each is `Option<Vec<T>>`, indexed in parallel with `atoms`, and all-or-nothing —
there is no meaningful state where only some atoms have a position.

| Table | Type | Comes from | Accessors |
|---|---|---|---|
| 2D layout | `Point2` | a layout pass, or a flat SDF | `coords`, `coord`, `set_coords`, … |
| 3D conformer | `Point3` | a file that states geometry: XYZ, PDB, Mol2, a 3D SDF | `coords3`, `coord3`, `set_coords3`, … |
| Atom sites | `AtomSite` | file columns: name, alt-loc, partial charge, occupancy, B-factor, radius | `sites`, `site`, `set_sites`, … |

**The layout and the conformer are different artefacts, and neither is derivable
from the other.** A layout is computed for drawing; a conformer is physical. A
molecule read from a 3D file and then laid out for display legitimately has both,
and depiction wants the first while a conversion to XYZ wants the second.
Projecting a conformer to get a layout superimposes atoms that differ only in
depth, which is why `Point3::to_2d` is explicit rather than a `From` impl.

`set_*` rejects a length that is not exactly one per atom, because the tables are
positional: a mismatch would silently attribute a charge or a B-factor to the
wrong atom.

### 3b. Residue and Chain Topology

Separate from the side tables, because chains and residues are not indexed in
parallel with the atoms — they are collections that *own* ranges of them.

```
Molecule
 ├── chains:   Vec<Chain>     — each owns a Range into residues
 └── residues: Vec<Residue>   — each owns a Range into atoms
```

Empty means absent; there is no `Option` wrapper, unlike the side tables above.

**Ranges, not a residue index on every atom.** PDB and mmCIF already order their
records this way and it costs nothing per atom. The price is real and enforced:
a residue's atoms must be *contiguous*, and the ranges must *ascend* — the
latter because `residue_of` binary-searches them. `set_topology` validates both,
plus that ranges stay in bounds and that a residue's `chain_ix` agrees with the
chain claiming it. An interleaved file is rejected rather than represented.

Ranges need not cover every atom. A ligand appended with no residue information
is a legal molecule, and `residue_of` answers `None` for it.

**Every residue is numbered twice.** `sequence` is the depositor's number
(`auth_seq_id`) — what PDB files carry and what papers cite — and `label_seq` is
mmCIF's canonical one, which is `None` for waters and ligands. They frequently
disagree and neither is derivable from the other, so reading one and writing the
other silently renumbers a structure. Chains carry the same pair as `id` and
`label_id`.

Residue identity is the **(chain, sequence, insertion code)** triple, never the
number alone: `58`, `58A` and `58B` are three residues, and two chains may each
have their own `58`.

**`add_atom` drops all three tables and the topology.** It is the only method
that can change the atom count — there is no `remove_atom`, and `atoms_mut()`
returns a slice, so the length cannot change through it. That single choke point
is what makes the parallel-array model safe. Appending does not strictly
invalidate a residue range, but it leaves an atom belonging to no residue that a
PDB write would silently drop, so one rule covers everything.

### 4. Graph Enables Queries

- O(1) neighbor lookup
- DFS/BFS traversal
- Radius-based search (for fingerprints)
- Connected components

## Data Flow

```
User Input
    ↓
Create Atoms → Add to Molecule
    ↓
Create Bonds → Add to Molecule → Updates Graph
    ↓
Calculate Implicit H → Uses Graph + Valence Rules
    ↓
Query Properties → Uses Atoms + Bonds + Graph
```

## Common Operations

```rust
// Check connectivity
let degree = mol.degree(atom_idx);
let neighbors = mol.neighbors(atom_idx);

// Get bond between atoms
if let Some(bond_idx) = mol.graph().get_bond(atom1, atom2) {
    let bond = mol.bond(bond_idx);
}

// Traverse molecule
let all_atoms = mol.graph().dfs(start_atom);

// Find atoms within distance
let nearby = mol.graph().neighbors_within_radius(atom_idx, 2);

// Get molecular properties
let formula = mol.formula();
let weight = mol.molecular_weight();
```

## Testing

```bash
cargo test -p chem
```

## Documentation

```bash
cargo doc -p chem --open
```

## Why This Design?

- **Index-based**: Bonds store atom indices (not references) → GPU-ready
- **Adjacency list**: Sparse graph representation → Memory efficient
- **Separation**: Atoms, bonds, graph kept separate → Clean, parallelizable
- **Builder pattern**: Fluent API → Ergonomic
- **Type safety**: Enums for bond order, element → Compile-time checks
