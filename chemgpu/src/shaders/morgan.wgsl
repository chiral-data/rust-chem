// Morgan Fingerprint GPU Shader
// Implements parallel Morgan algorithm with iterative neighbor aggregation

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct MoleculeData {
    num_atoms: u32,
    num_bonds: u32,
    radius: u32,
    fp_size: u32,
    use_chirality: u32,
    use_bond_types: u32,
    only_nonzero_invariants: u32,
    _padding: u32,
}

struct IterationParams {
    current_layer: u32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

struct Atom {
    atomic_number: u32,
    degree: u32,
    is_aromatic: u32,
    chirality: u32,  // 0=None, 1=Unspecified, 2=CCW, 3=CW
    formal_charge: i32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

struct Bond {
    atom1: u32,
    atom2: u32,
    order: u32,      // 1=Single, 2=Double, 3=Triple, 12=Aromatic
    stereo: u32,     // 0=None, 1=E, 2=Z
    is_aromatic: u32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

struct Neighbor {
    atom_idx: u32,
    bond_idx: u32,
}

// ============================================================================
// BUFFER BINDINGS
// ============================================================================

@group(0) @binding(0)
var<storage, read> molecule_data: MoleculeData;

@group(0) @binding(1)
var<storage, read> atoms: array<Atom>;

@group(0) @binding(2)
var<storage, read> bonds: array<Bond>;

@group(0) @binding(3)
var<storage, read> adjacency: array<Neighbor>;  // Flat adjacency list

@group(0) @binding(4)
var<storage, read> adjacency_offsets: array<u32>;  // Start index for each atom

@group(0) @binding(5)
var<storage, read_write> current_invariants: array<atomic<u32>>;

@group(0) @binding(6)
var<storage, read_write> next_invariants: array<atomic<u32>>;

@group(0) @binding(7)
var<storage, read_write> fingerprint: array<atomic<u32>>;  // Packed bits

@group(0) @binding(8)
var<uniform> iter_params: IterationParams;

// ============================================================================
// HASH FUNCTIONS (matching Rust implementation)
// ============================================================================

fn hash_combine(seed: ptr<function, u32>, value: u32) {
    // let hash = value;
    // FIX: Why these are different
    // *seed = *seed ^ (hash + 0x9e3779b9u + (*seed << 6u) + (*seed >> 2u));

    // let hash = value;
    //let temp = hash + 0x9e3779b9u;
    //let temp2 = temp + (*seed << 6u);
    //let temp3 = temp2 + (*seed >> 2u);
    //*seed = *seed ^ temp3;
    *seed = (*seed ^ value);
}

fn hash_pair(seed: ptr<function, u32>, first: i32, second: u32) {
    // Cast i32 to u32 for hashing
    let first_u = bitcast<u32>(first);
    hash_combine(seed, first_u);
    hash_combine(seed, second);
}

fn hash_to_bit_position(hash: u32, fp_size: u32) -> u32 {
    return hash % fp_size;
}

// ============================================================================
// ATOMIC BIT OPERATIONS
// ============================================================================

fn set_bit_atomic(bit_pos: u32) {
    let word_idx = bit_pos / 32u;
    let bit_idx = bit_pos % 32u;
    let mask = 1u << bit_idx;
    
    atomicOr(&fingerprint[word_idx], mask);
}

// ============================================================================
// INITIALIZATION KERNEL (Layer 0)
// ============================================================================

@compute @workgroup_size(64)
fn init_invariants(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let atom_idx = global_id.x;
    
    if (atom_idx >= molecule_data.num_atoms) {
        return;
    }
    
    let atom = atoms[atom_idx];
    
    // Compute initial invariant (matching MorganAtomInvGenerator)
    var inv = atom.atomic_number;
    inv = inv * 100u + atom.degree;
    
    if (atom.is_aromatic != 0u) {
        inv = inv * 2u;
    }
    
    inv = inv * 31u;
    
    // Store initial invariant
    atomicStore(&current_invariants[atom_idx], inv);
    
    // Add to fingerprint if not filtering zero invariants
    if (molecule_data.only_nonzero_invariants == 0u || inv != 0u) {
        let bit_pos = hash_to_bit_position(inv, molecule_data.fp_size);
        set_bit_atomic(bit_pos);
    }

    // Debug: print first atom's invariant
    if (atom_idx == 0u) {
        // Can't printf in WGSL, but store in unused buffer slot
    }
}

// ============================================================================
// MORGAN ITERATION KERNEL (Layers 1..radius)
// ============================================================================

@compute @workgroup_size(64)
fn morgan_iteration(
    @builtin(global_invocation_id) global_id: vec3<u32>,
    @builtin(workgroup_id) workgroup_id: vec3<u32>
) {
    let atom_idx = global_id.x;
    
    if (atom_idx >= molecule_data.num_atoms) {
        return;
    }
    
    let atom = atoms[atom_idx];
    let degree = atom.degree;
    
    if (degree == 0u) {
        atomicStore(&next_invariants[atom_idx], 0u);
        return;
    }
    
    // Collect neighbor invariants
    let start_idx = adjacency_offsets[atom_idx];
    let end_idx = adjacency_offsets[atom_idx + 1u];
    
    // Use shared memory for neighbor sorting (workgroup local)
    var neighbor_data: array<vec2<u32>, 16>;  // (bond_inv, atom_inv), max 16 neighbors
    var neighbor_count = 0u;
    
    // Gather neighbors
    for (var i = start_idx; i < end_idx; i = i + 1u) {
        if (neighbor_count >= 16u) { break; }
        
        let neighbor = adjacency[i];
        let bond = bonds[neighbor.bond_idx];
        
        // Get bond invariant
        var bond_inv = 1u;
        if (molecule_data.use_bond_types != 0u) {
            bond_inv = bond.order;
        }
        
        let neighbor_inv = atomicLoad(&current_invariants[neighbor.atom_idx]);
        
        neighbor_data[neighbor_count] = vec2<u32>(bond_inv, neighbor_inv);
        neighbor_count = neighbor_count + 1u;
    }
    
    // Bubble sort neighbors (small N, simple is fine)
    for (var i = 0u; i < neighbor_count; i = i + 1u) {
        for (var j = i + 1u; j < neighbor_count; j = j + 1u) {
            let a = neighbor_data[i];
            let b = neighbor_data[j];
            
            // Sort by bond_inv first, then atom_inv
            if (a.x > b.x || (a.x == b.x && a.y > b.y)) {
                neighbor_data[i] = b;
                neighbor_data[j] = a;
            }
        }
    }
    
    // Hash the sorted neighborhood
    let current_inv = atomicLoad(&current_invariants[atom_idx]);
    var new_inv = iter_params.current_layer;// workgroup_id.y;  // This is the layer number passed in
    hash_combine(&new_inv, current_inv);
    
    for (var i = 0u; i < neighbor_count; i = i + 1u) {
        let bond_inv = i32(neighbor_data[i].x);
        let atom_inv = neighbor_data[i].y;
        hash_pair(&new_inv, bond_inv, atom_inv);
    }
    
    // Add chirality if needed
    if (molecule_data.use_chirality != 0u && atom.chirality > 1u) {
        hash_combine(&new_inv, atom.chirality);
    }
    
    // Store new invariant
    atomicStore(&next_invariants[atom_idx], new_inv);
    
    // Add to fingerprint
    if (molecule_data.only_nonzero_invariants == 0u || current_inv != 0u) {
        let bit_pos = hash_to_bit_position(new_inv, molecule_data.fp_size);
        set_bit_atomic(bit_pos);
    }
}
