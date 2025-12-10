# Morgan Fingerprint Implementation

## Overview

This Rust implementation generates Morgan (ECFP) fingerprints for molecular structures. Morgan fingerprints are circular fingerprints that encode the local chemical environment around each atom up to a specified radius.

## Algorithm Overview

### 1. Initialization (Radius 0)

Each atom is assigned an initial hash based on its local properties:

```rust
fn compute_atom_hash(atom: &Atom, degree: usize) -> u32 {
    let atomic_num = atom.atomic_number() as u32;
    let valence = (degree as u32) + (atom.total_hydrogens() as u32);
    let charge = (atom.formal_charge() + 5) as u32;  // Offset by 5

    let mut hash = atomic_num;
    hash = hash.wrapping_mul(37).wrapping_add(valence);
    hash = hash.wrapping_mul(37).wrapping_add(charge);

    if atom.is_aromatic() {
        hash = hash.wrapping_mul(37).wrapping_add(1);
    }

    hash
}
```

**Properties encoded:**

- Atomic number
- Degree (number of heavy atom neighbors)
- Total hydrogen count (implicit + explicit)
- Formal charge
- Aromaticity flag

### 2. Iterative Updates (Radius 1 to N)

For each iteration, atom identifiers are updated by incorporating neighbor information:

```rust
fn update_identifiers(mol: &Molecule, old_ids: &[u32]) -> Vec<u32> {
    for each atom:
        1. Collect neighbor hashes combined with bond orders
        2. Sort neighbor hashes (canonical ordering)
        3. Combine current hash with sorted neighbor hashes

    return new_identifiers
}
```

**Bond encoding:**

- Single bond → 1
- Double bond → 2
- Triple bond → 3
- Aromatic bond → 4

**Combination formula:**

```
neighbor_combined = neighbor_id × 100 + bond_order
new_hash = old_hash × 37 + neighbor_combined[0]
new_hash = new_hash × 37 + neighbor_combined[1]
...
```

### 3. Feature Collection

All unique hashes across all radii (0 to N) are collected in a HashSet to eliminate duplicates.

### 4. Bit Vector Generation

Each unique hash is mapped to a bit position in the fingerprint:

```
bit_position = hash % nbits
```

The bit at that position is set to 1.

## Example: Ethane (CC)

### Radius 0 (Initial Hashes)

```
Atom 0 (C): degree=1, H=3
  hash = hash_combine(6, 4, 5, 0) = 0x05030106

Atom 1 (C): degree=1, H=3
  hash = hash_combine(6, 4, 5, 0) = 0x05030106

Unique hashes: [0x05030106]  → 1 unique feature
```

### Radius 1 (First Neighborhood)

```
Atom 0 (C): neighbor C via single bond
  neighbor_hash = 0x05030106 × 100 + 1
  new_hash = 0x05030106 × 37 + neighbor_hash = 0x175584b1

Atom 1 (C): identical environment
  new_hash = 0x175584b1

Unique hashes: [0x05030106, 0x175584b1]  → 2 unique features
```

### Radius 2 (Extended Neighborhood)

```
Atom 0 (C): neighbor C (which is bonded to C)
  new_hash = 0xf590b2a1

Atom 1 (C): identical environment
  new_hash = 0xf590b2a1

Unique hashes: [0x05030106, 0x175584b1, 0xf590b2a1]  → 3 unique features
```

### Final Fingerprint

```
3 unique hashes → 3 bit positions in 2048-bit vector
bit_positions = [821, 1406, 1669]
```

## Key Features

### Symmetry Preservation

Identical atoms in symmetric molecules produce identical hashes:

- Both carbons in ethane generate the same hash at each radius
- Only unique chemical environments contribute to the fingerprint

### Canonicalization

Neighbor hashes are **sorted** before combination to ensure:

- Atom indexing order doesn't affect the fingerprint
- Isomorphic molecules produce identical fingerprints

### Collision Handling

Multiple features may map to the same bit position:

- This is expected behavior (folding)
- Reduces fingerprint size while maintaining similarity
- Example: 6 features → 4 bit positions (2 collisions)

## Configuration Parameters

```rust
pub struct MorganOptions {
    pub radius: u32,           // Number of bond hops (2 = ECFP4)
    pub nbits: usize,          // Fingerprint length (typically 2048)
    pub bits_per_feature: usize, // Bits set per feature (1 for standard)
}
```

**Common settings:**

- ECFP4: radius=2, nbits=2048
- ECFP6: radius=3, nbits=2048

## Implementation Details

### Data Structures

- **BitVec**: Compact bit vector using `Vec<u64>` for storage
- **HashSet<u32>**: Deduplication of features across radii
- **Vec<u32>**: Atom identifier storage

### Performance Characteristics

- Time complexity: O(radius × atoms × avg_degree)
- Space complexity: O(atoms + bits)
- Typical performance: <1ms for small molecules (<100 atoms)

## Current Status

### Known Issues

- Generated features (count) match RDKit (sometimes) ✓
- Bit positions differ from RDKit (hash function mismatch)
- Similarity: 0% due to different hash→bit mapping

**Root cause**: The hash combination formula differs slightly from RDKit's implementation, causing different bit positions despite correct feature counts.

## Testing

```bash
# Build CLI tool
cargo build --bin cli

# Generate fingerprint for ethane
./target/debug/cli "CC" 2 2048

# Compare with RDKit
python scripts/validate_morgan.py
```

## References

1. Rogers, D. & Hahn, M. "Extended-Connectivity Fingerprints" (2010)
2. RDKit Morgan Fingerprint Documentation
3. Boost hash_combine algorithm (used in RDKit)
