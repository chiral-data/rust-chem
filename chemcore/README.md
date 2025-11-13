# ChemCore

Core molecular data structures for cheminformatics applications.

## Overview

ChemCore provides the fundamental building blocks for representing molecules:

- **Atoms**: Elements with properties (charge, isotope, aromaticity)
- **Bonds**: Connections with order and stereochemistry
- **Molecules**: Complete molecular structures
- **Graph**: Efficient connectivity queries

## Quick Start

```rust
use chemcore::prelude::*;

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
use chemcore::prelude::*;

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

### 2. Bonds Store Connectivity

- Two atom indices
- Bond order (single, double, triple, aromatic)
- Stereochemistry (E/Z)

### 3. Molecule Manages Everything

- Collection of atoms
- Collection of bonds
- Graph for fast queries
- Properties (metadata)

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
cargo test -p chemcore
```

## Documentation

```bash
cargo doc -p chemcore --open
```

## Why This Design?

- **Index-based**: Bonds store atom indices (not references) → GPU-ready
- **Adjacency list**: Sparse graph representation → Memory efficient
- **Separation**: Atoms, bonds, graph kept separate → Clean, parallelizable
- **Builder pattern**: Fluent API → Ergonomic
- **Type safety**: Enums for bond order, element → Compile-time checks

## Next Steps

Week 2 will add SMILES parsing to create molecules from strings like `"CCO"` (ethanol).
