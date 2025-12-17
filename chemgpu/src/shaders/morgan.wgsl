// Batch Morgan Fingerprint GPU Shader
// Processes multiple molecules in parallel for maximum GPU utilization

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct BatchParams {
    num_molecules: u32,
    max_atoms: u32,       // Maximum atoms across all molecules (for bounds)
    radius: u32,
    fp_size: u32,
    fp_words: u32,        // Number of u32 words per fingerprint
    use_chirality: u32,
    use_bond_types: u32,
    only_nonzero_invariants: u32,
}

struct MoleculeOffset {
    atom_start: u32,
    atom_count: u32,
    bond_start: u32,
    bond_count: u32,
    adj_start: u32,
    fp_start: u32,        // Start index in fingerprint output buffer
}

struct Atom {
    atomic_number: u32,
    degree: u32,
    is_aromatic: u32,
    chirality: u32,
    formal_charge: i32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

struct Bond {
    atom1: u32,           // Local atom index within molecule
    atom2: u32,
    order: u32,
    stereo: u32,
    is_aromatic: u32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

struct Neighbor {
    atom_idx: u32,        // Local atom index within molecule
    bond_idx: u32,        // Local bond index within molecule
}

struct IterationParams {
    current_layer: u32,
    _padding1: u32,
    _padding2: u32,
    _padding3: u32,
}

// ============================================================================
// BUFFER BINDINGS
// ============================================================================

@group(0) @binding(0)
var<uniform> params: BatchParams;

@group(0) @binding(1)
var<storage, read> molecule_offsets: array<MoleculeOffset>;

@group(0) @binding(2)
var<storage, read> atoms: array<Atom>;

@group(0) @binding(3)
var<storage, read> bonds: array<Bond>;

@group(0) @binding(4)
var<storage, read> adjacency: array<Neighbor>;

@group(0) @binding(5)
var<storage, read> adjacency_offsets: array<u32>;

@group(0) @binding(6)
var<storage, read_write> current_invariants: array<atomic<u32>>;

@group(0) @binding(7)
var<storage, read_write> next_invariants: array<atomic<u32>>;

@group(0) @binding(8)
var<storage, read_write> fingerprints: array<atomic<u32>>;

@group(0) @binding(9)
var<uniform> iter_params: IterationParams;

// ============================================================================
// HASH FUNCTIONS
// ============================================================================

fn hash_combine(seed: ptr<function, u32>, value: u32) {
    *seed = (*seed ^ value);
}

fn hash_pair(seed: ptr<function, u32>, first: i32, second: u32) {
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

fn set_bit_atomic(fp_start: u32, bit_pos: u32) {
    let word_idx = fp_start + (bit_pos / 32u);
    let bit_idx = bit_pos % 32u;
    let mask = 1u << bit_idx;
    atomicOr(&fingerprints[word_idx], mask);
}

// ============================================================================
// BATCH INITIALIZATION KERNEL
// Each thread processes one atom from one molecule
// ============================================================================

@compute @workgroup_size(256)
fn init_invariants_batch(@builtin(global_invocation_id) global_id: vec3<u32>) {
    // global_id.x is the global atom index across all molecules
    let global_atom_idx = global_id.x;

    // Find which molecule this atom belongs to (binary search would be faster for many molecules)
    var mol_idx = 0u;
    var local_atom_idx = global_atom_idx;

    for (var m = 0u; m < params.num_molecules; m = m + 1u) {
        let mol = molecule_offsets[m];
        if (local_atom_idx < mol.atom_count) {
            mol_idx = m;
            break;
        }
        local_atom_idx = local_atom_idx - mol.atom_count;
    }

    let mol = molecule_offsets[mol_idx];

    // Bounds check
    if (local_atom_idx >= mol.atom_count) {
        return;
    }

    let atom = atoms[mol.atom_start + local_atom_idx];

    // Compute initial invariant
    var inv = atom.atomic_number;
    inv = inv * 100u + atom.degree;

    if (atom.is_aromatic != 0u) {
        inv = inv * 2u;
    }

    inv = inv * 31u;

    // Store initial invariant (global position)
    atomicStore(&current_invariants[mol.atom_start + local_atom_idx], inv);

    // Add to fingerprint
    if (params.only_nonzero_invariants == 0u || inv != 0u) {
        let bit_pos = hash_to_bit_position(inv, params.fp_size);
        set_bit_atomic(mol.fp_start, bit_pos);
    }
}

// ============================================================================
// BATCH ITERATION KERNEL
// ============================================================================

@compute @workgroup_size(256)
fn morgan_iteration_batch(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let global_atom_idx = global_id.x;

    // Find which molecule this atom belongs to
    var mol_idx = 0u;
    var local_atom_idx = global_atom_idx;
    var found = false;

    for (var m = 0u; m < params.num_molecules; m = m + 1u) {
        let mol = molecule_offsets[m];
        if (local_atom_idx < mol.atom_count) {
            mol_idx = m;
            found = true;
            break;
        }
        local_atom_idx = local_atom_idx - mol.atom_count;
    }

    if (!found) {
        return;
    }

    let mol = molecule_offsets[mol_idx];

    if (local_atom_idx >= mol.atom_count) {
        return;
    }

    let atom = atoms[mol.atom_start + local_atom_idx];
    let degree = atom.degree;

    if (degree == 0u) {
        atomicStore(&next_invariants[mol.atom_start + local_atom_idx], 0u);
        return;
    }

    // Get adjacency range for this atom
    let adj_base = mol.adj_start;
    let atom_adj_start = adjacency_offsets[mol.atom_start + local_atom_idx];
    let atom_adj_end = adjacency_offsets[mol.atom_start + local_atom_idx + 1u];

    // Collect neighbor invariants
    var neighbor_data: array<vec2<u32>, 16>;
    var neighbor_count = 0u;

    for (var i = atom_adj_start; i < atom_adj_end; i = i + 1u) {
        if (neighbor_count >= 16u) { break; }

        let neighbor = adjacency[i];
        let bond = bonds[mol.bond_start + neighbor.bond_idx];

        var bond_inv = 1u;
        if (params.use_bond_types != 0u) {
            bond_inv = bond.order;
        }

        let neighbor_inv = atomicLoad(&current_invariants[mol.atom_start + neighbor.atom_idx]);

        neighbor_data[neighbor_count] = vec2<u32>(bond_inv, neighbor_inv);
        neighbor_count = neighbor_count + 1u;
    }

    // Sort neighbors
    for (var i = 0u; i < neighbor_count; i = i + 1u) {
        for (var j = i + 1u; j < neighbor_count; j = j + 1u) {
            let a = neighbor_data[i];
            let b = neighbor_data[j];
            if (a.x > b.x || (a.x == b.x && a.y > b.y)) {
                neighbor_data[i] = b;
                neighbor_data[j] = a;
            }
        }
    }

    // Hash the sorted neighborhood
    let current_inv = atomicLoad(&current_invariants[mol.atom_start + local_atom_idx]);
    var new_inv = iter_params.current_layer;
    hash_combine(&new_inv, current_inv);

    for (var i = 0u; i < neighbor_count; i = i + 1u) {
        let bond_inv = i32(neighbor_data[i].x);
        let atom_inv = neighbor_data[i].y;
        hash_pair(&new_inv, bond_inv, atom_inv);
    }

    if (params.use_chirality != 0u && atom.chirality > 1u) {
        hash_combine(&new_inv, atom.chirality);
    }

    // Store new invariant
    atomicStore(&next_invariants[mol.atom_start + local_atom_idx], new_inv);

    // Add to fingerprint
    if (params.only_nonzero_invariants == 0u || current_inv != 0u) {
        let bit_pos = hash_to_bit_position(new_inv, params.fp_size);
        set_bit_atomic(mol.fp_start, bit_pos);
    }
}

// ============================================================================
// COPY KERNEL (next -> current invariants)
// ============================================================================

@compute @workgroup_size(256)
fn copy_invariants(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    if (idx >= params.max_atoms) {
        return;
    }
    let val = atomicLoad(&next_invariants[idx]);
    atomicStore(&current_invariants[idx], val);
}
