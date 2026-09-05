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
    // Aggregate sizes needed to compute the read-write scratch regions'
    // offsets below — not derivable from the fields above alone, since they
    // depend on each molecule's own bond count (see the offset functions).
    total_neighborhood_words: u32,
    total_seen_words: u32,
}

struct MoleculeOffset {
    atom_start: u32,
    atom_count: u32,
    bond_start: u32,
    bond_count: u32,
    adj_start: u32,
    fp_start: u32,          // Start index in fingerprint output buffer
    neighborhood_start: u32, // Start word index in the neighborhood-bitset buffers
    seen_start: u32,         // Start word index in the seen_neighborhoods buffer
}

struct Atom {
    atomic_number: u32,
    degree: u32,
    is_aromatic: u32,
    chirality: u32,
    formal_charge: i32,
    total_hydrogens: u32,
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
//
// All read-write working state (current/next invariants, fingerprints,
// dedup bookkeeping, neighbor scratch — 9 logically separate arrays in the
// original design) is packed into one combined `scratch` storage buffer
// instead of one binding each. Some WebGPU implementations only guarantee 8
// storage-buffer bindings per compute stage (the wgpu/WebGPU baseline;
// browsers have been observed reporting as low as 8, see #50), well under
// what 9 separate read-write buffers plus 5 read-only ones would need.
// Every element uses atomic ops even where nothing actually contends for
// it (each index is only ever touched by one invocation per pass) — that's
// required for storage buffers, since a plain (non-atomic) read_write
// binding can't be built from parts written by different logical "buffers"
// without WGSL treating the whole declaration as one array of the same
// scalar type throughout.
//
// The offset functions below (`*_base()`) are the single source of truth
// for where each logical region starts; every kernel goes through them (or
// wrapper accessors built on them) rather than recomputing offsets inline,
// so all five kernels agree on the same layout.
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
var<storage, read_write> scratch: array<atomic<u32>>;

@group(0) @binding(7)
var<uniform> iter_params: IterationParams;

// ============================================================================
// SCRATCH ARENA LAYOUT (offsets, in words, into `scratch`)
// ============================================================================

fn current_invariants_base() -> u32 {
    return 0u;
}

fn next_invariants_base() -> u32 {
    return params.max_atoms;
}

fn fingerprints_base() -> u32 {
    return 2u * params.max_atoms;
}

fn dead_atoms_base() -> u32 {
    return fingerprints_base() + params.num_molecules * params.fp_words;
}

fn current_neighborhoods_base() -> u32 {
    return dead_atoms_base() + params.max_atoms;
}

fn next_neighborhoods_base() -> u32 {
    return current_neighborhoods_base() + params.total_neighborhood_words;
}

fn seen_neighborhoods_base() -> u32 {
    return next_neighborhoods_base() + params.total_neighborhood_words;
}

fn seen_count_base() -> u32 {
    return seen_neighborhoods_base() + params.total_seen_words;
}

fn neighbor_scratch_base() -> u32 {
    return seen_count_base() + params.num_molecules;
}

// ============================================================================
// SCRATCH ARENA ACCESSORS
// ============================================================================

fn get_current_invariant(idx: u32) -> u32 {
    return atomicLoad(&scratch[current_invariants_base() + idx]);
}

fn set_current_invariant(idx: u32, value: u32) {
    atomicStore(&scratch[current_invariants_base() + idx], value);
}

fn get_next_invariant(idx: u32) -> u32 {
    return atomicLoad(&scratch[next_invariants_base() + idx]);
}

fn set_next_invariant(idx: u32, value: u32) {
    atomicStore(&scratch[next_invariants_base() + idx], value);
}

fn set_fingerprint_bit(fp_start: u32, bit_pos: u32) {
    let word_idx = fingerprints_base() + fp_start + (bit_pos / 32u);
    let bit_idx = bit_pos % 32u;
    atomicOr(&scratch[word_idx], 1u << bit_idx);
}

fn get_dead_atom(idx: u32) -> u32 {
    return atomicLoad(&scratch[dead_atoms_base() + idx]);
}

fn set_dead_atom(idx: u32, value: u32) {
    atomicStore(&scratch[dead_atoms_base() + idx], value);
}

fn get_current_neighborhood(idx: u32) -> u32 {
    return atomicLoad(&scratch[current_neighborhoods_base() + idx]);
}

fn set_current_neighborhood(idx: u32, value: u32) {
    atomicStore(&scratch[current_neighborhoods_base() + idx], value);
}

fn get_next_neighborhood(idx: u32) -> u32 {
    return atomicLoad(&scratch[next_neighborhoods_base() + idx]);
}

fn set_next_neighborhood(idx: u32, value: u32) {
    atomicStore(&scratch[next_neighborhoods_base() + idx], value);
}

fn or_next_neighborhood(idx: u32, mask: u32) {
    atomicOr(&scratch[next_neighborhoods_base() + idx], mask);
}

fn get_seen_neighborhood(idx: u32) -> u32 {
    return atomicLoad(&scratch[seen_neighborhoods_base() + idx]);
}

fn set_seen_neighborhood(idx: u32, value: u32) {
    atomicStore(&scratch[seen_neighborhoods_base() + idx], value);
}

fn get_seen_count(idx: u32) -> u32 {
    return atomicLoad(&scratch[seen_count_base() + idx]);
}

fn set_seen_count(idx: u32, value: u32) {
    atomicStore(&scratch[seen_count_base() + idx], value);
}

fn get_neighbor_scratch_x(idx: u32) -> u32 {
    return atomicLoad(&scratch[neighbor_scratch_base() + 2u * idx]);
}

fn get_neighbor_scratch_y(idx: u32) -> u32 {
    return atomicLoad(&scratch[neighbor_scratch_base() + 2u * idx + 1u]);
}

fn set_neighbor_scratch(idx: u32, x: u32, y: u32) {
    atomicStore(&scratch[neighbor_scratch_base() + 2u * idx], x);
    atomicStore(&scratch[neighbor_scratch_base() + 2u * idx + 1u], y);
}

// ============================================================================
// HASH FUNCTIONS
// ============================================================================

// boost::hash_combine (what RDKit uses for std::pair/std::vector hashing) —
// NOT a bare XOR: XOR-only combine is self-cancelling when two combined
// values are equal (e.g. two identical neighbor invariants), silently
// dropping one of them from the resulting hash. Mirrors chemfp/src/hash.rs.
fn hash_combine(seed: ptr<function, u32>, value: u32) {
    *seed = (*seed ^ (value + 0x9e3779b9u + (*seed << 6u) + (*seed >> 2u)));
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

    // Compute initial invariant. Mirrors MorganAtomInvGenerator::
    // get_connectivity_invariants on the CPU side exactly, term for term and
    // in the same order, so CPU and GPU hash identically (#192).
    var inv = atom.atomic_number;
    inv = inv * 100u + atom.degree;
    inv = inv * 100u + atom.total_hydrogens;
    inv = inv * 20u + u32(atom.formal_charge + 5);

    if (atom.is_aromatic != 0u) {
        inv = inv * 2u;
    }

    inv = inv * 31u;

    // Store initial invariant (global position)
    set_current_invariant(mol.atom_start + local_atom_idx, inv);

    // Round 0 has no reachable-bond neighborhood yet (every atom starts with
    // an empty bitset), so — matching the CPU reference — it is never subject
    // to the redundant-environment check: every atom unconditionally
    // contributes its round-0 bit.
    if (params.only_nonzero_invariants == 0u || inv != 0u) {
        let bit_pos = hash_to_bit_position(inv, params.fp_size);
        set_fingerprint_bit(mol.fp_start, bit_pos);
    }
}

// ============================================================================
// BATCH ITERATION KERNEL
// Computes this round's invariant AND reachable-bond neighborhood per atom.
// Fingerprint-bit contribution is decided afterwards, in dedup_environments_batch,
// once the full set of this round's candidate environments is known.
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

    let global_idx = mol.atom_start + local_atom_idx;
    let bond_words = (mol.bond_count + 31u) / 32u;
    let neigh_base = mol.neighborhood_start + local_atom_idx * bond_words;

    // Atoms frozen by a previous round's redundant-environment check (or a
    // previous round's degree-0 check, below) simply carry their state
    // forward unchanged — they contribute no further environments.
    if (get_dead_atom(global_idx) != 0u) {
        let inv = get_current_invariant(global_idx);
        set_next_invariant(global_idx, inv);
        for (var w = 0u; w < bond_words; w = w + 1u) {
            set_next_neighborhood(neigh_base + w, get_current_neighborhood(neigh_base + w));
        }
        return;
    }

    let atom = atoms[mol.atom_start + local_atom_idx];
    let degree = atom.degree;

    if (degree == 0u) {
        set_next_invariant(global_idx, 0u);
        set_dead_atom(global_idx, 1u);
        return;
    }

    // Get adjacency range for this atom
    let atom_adj_start = adjacency_offsets[mol.atom_start + local_atom_idx];
    let atom_adj_end = adjacency_offsets[mol.atom_start + local_atom_idx + 1u];

    for (var w = 0u; w < bond_words; w = w + 1u) {
        set_next_neighborhood(neigh_base + w, 0u);
    }

    // Collect neighbor invariants directly into neighbor_scratch, indexed
    // identically to `adjacency` (so this atom's own slice is exactly
    // neighbor_scratch[atom_adj_start..atom_adj_end]) — unbounded, no cap on
    // degree, unlike a fixed-size function-local array.
    for (var i = atom_adj_start; i < atom_adj_end; i = i + 1u) {
        let neighbor = adjacency[i];
        let bond = bonds[mol.bond_start + neighbor.bond_idx];

        // Mark this incident bond as reached.
        let bit_word = neighbor.bond_idx / 32u;
        let bit_bit = neighbor.bond_idx % 32u;
        or_next_neighborhood(neigh_base + bit_word, 1u << bit_bit);

        // Fold in the neighbor's own (previously committed) reachable-bond set.
        let neighbor_neigh_base = mol.neighborhood_start + neighbor.atom_idx * bond_words;
        for (var w = 0u; w < bond_words; w = w + 1u) {
            or_next_neighborhood(neigh_base + w, get_current_neighborhood(neighbor_neigh_base + w));
        }

        var bond_inv = 1u;
        if (params.use_bond_types != 0u) {
            bond_inv = bond.order;
        }
        let neighbor_inv = get_current_invariant(mol.atom_start + neighbor.atom_idx);
        set_neighbor_scratch(i, bond_inv, neighbor_inv);
    }

    // Sort this atom's neighbor_scratch slice in place.
    for (var i = atom_adj_start; i < atom_adj_end; i = i + 1u) {
        for (var j = i + 1u; j < atom_adj_end; j = j + 1u) {
            let a_x = get_neighbor_scratch_x(i);
            let a_y = get_neighbor_scratch_y(i);
            let b_x = get_neighbor_scratch_x(j);
            let b_y = get_neighbor_scratch_y(j);
            if (a_x > b_x || (a_x == b_x && a_y > b_y)) {
                set_neighbor_scratch(i, b_x, b_y);
                set_neighbor_scratch(j, a_x, a_y);
            }
        }
    }

    // Hash the sorted neighborhood
    let current_inv = get_current_invariant(global_idx);
    var new_inv = iter_params.current_layer;
    hash_combine(&new_inv, current_inv);

    for (var i = atom_adj_start; i < atom_adj_end; i = i + 1u) {
        let bond_inv = i32(get_neighbor_scratch_x(i));
        let atom_inv = get_neighbor_scratch_y(i);
        hash_pair(&new_inv, bond_inv, atom_inv);
    }

    if (params.use_chirality != 0u && atom.chirality > 1u) {
        hash_combine(&new_inv, atom.chirality);
    }

    // Store new invariant; fingerprint-bit contribution happens later, in
    // dedup_environments_batch.
    set_next_invariant(global_idx, new_inv);
}

// ============================================================================
// DEDUP KERNEL
// One invocation per molecule (workgroup_size 1): decides, for this round's
// freshly computed environments, which are redundant (an identical reachable-
// bond neighborhood already accepted, this round or an earlier one) and
// therefore contribute no fingerprint bit and freeze that atom (dead_atoms)
// from all further rounds — mirroring the CPU reference implementation's
// `neighborhoods: FxHashSet<BitVec>` / `dead_atoms` logic exactly, just
// sequential-per-molecule instead of sequential-per-whole-batch, since
// different molecules' environments can never collide with each other.
// ============================================================================

const MAX_ATOMS_PER_MOL: u32 = 256u;

fn candidate_less(mol: MoleculeOffset, bond_words: u32, a_local: u32, b_local: u32) -> bool {
    let a_base = mol.neighborhood_start + a_local * bond_words;
    let b_base = mol.neighborhood_start + b_local * bond_words;

    for (var w = 0u; w < bond_words; w = w + 1u) {
        let av = get_next_neighborhood(a_base + w);
        let bv = get_next_neighborhood(b_base + w);
        if (av != bv) {
            return av < bv;
        }
    }

    let a_inv = get_next_invariant(mol.atom_start + a_local);
    let b_inv = get_next_invariant(mol.atom_start + b_local);
    if (a_inv != b_inv) {
        return a_inv < b_inv;
    }

    return a_local < b_local;
}

// workgroup_size(256) so the dispatch below chunks by molecule count instead
// of mapping num_molecules directly 1:1 to the WebGPU dispatch-group count,
// which WebGPU caps at 65,535 per dimension (see #27). The mol_idx bounds
// check below already tolerates the tail workgroup overshooting num_molecules.
@compute @workgroup_size(256)
fn dedup_environments_batch(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let mol_idx = global_id.x;
    if (mol_idx >= params.num_molecules) {
        return;
    }

    let mol = molecule_offsets[mol_idx];
    let bond_words = (mol.bond_count + 31u) / 32u;

    // Gather this round's non-frozen atoms as dedup candidates.
    // Molecules with more than MAX_ATOMS_PER_MOL atoms are truncated to the
    // first MAX_ATOMS_PER_MOL (by local index) — a documented limitation,
    // matching the codebase's existing 16-neighbor cap precedent.
    var cand_atom: array<u32, 256>;
    var cand_count = 0u;

    for (var i = 0u; i < mol.atom_count; i = i + 1u) {
        if (cand_count >= MAX_ATOMS_PER_MOL) {
            break;
        }
        let global_idx = mol.atom_start + i;
        if (get_dead_atom(global_idx) == 0u) {
            cand_atom[cand_count] = i;
            cand_count = cand_count + 1u;
        }
    }

    // Stable-equivalent selection sort by (neighborhood, invariant, atom_idx),
    // mirroring the CPU reference's `sort_by(|a,b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)))`.
    for (var i = 0u; i < cand_count; i = i + 1u) {
        var best = i;
        for (var j = i + 1u; j < cand_count; j = j + 1u) {
            if (candidate_less(mol, bond_words, cand_atom[j], cand_atom[best])) {
                best = j;
            }
        }
        if (best != i) {
            let tmp = cand_atom[i];
            cand_atom[i] = cand_atom[best];
            cand_atom[best] = tmp;
        }
    }

    // Walk in sorted order: the first atom to claim a given neighborhood
    // (this round, or a duplicate of one already accepted in an earlier
    // round) survives and contributes its bit; later duplicates are frozen.
    var seen_c = get_seen_count(mol_idx);

    for (var k = 0u; k < cand_count; k = k + 1u) {
        let local_idx = cand_atom[k];
        let global_idx = mol.atom_start + local_idx;
        let neigh_base = mol.neighborhood_start + local_idx * bond_words;

        var is_redundant = false;
        for (var s = 0u; s < seen_c; s = s + 1u) {
            let seen_base = mol.seen_start + s * bond_words;
            var equal = true;
            for (var w = 0u; w < bond_words; w = w + 1u) {
                if (get_next_neighborhood(neigh_base + w) != get_seen_neighborhood(seen_base + w)) {
                    equal = false;
                    break;
                }
            }
            if (equal) {
                is_redundant = true;
                break;
            }
        }

        if (is_redundant) {
            set_dead_atom(global_idx, 1u);
        } else {
            let seen_base = mol.seen_start + seen_c * bond_words;
            for (var w = 0u; w < bond_words; w = w + 1u) {
                set_seen_neighborhood(seen_base + w, get_next_neighborhood(neigh_base + w));
            }
            seen_c = seen_c + 1u;

            let inv = get_next_invariant(global_idx);
            if (params.only_nonzero_invariants == 0u || inv != 0u) {
                let bit_pos = hash_to_bit_position(inv, params.fp_size);
                set_fingerprint_bit(mol.fp_start, bit_pos);
            }
        }
    }

    set_seen_count(mol_idx, seen_c);
}

// ============================================================================
// COPY KERNELS (next -> current, run after dedup so dedup sees this round's
// freshly computed values before they become "current" for the next round)
// ============================================================================

@compute @workgroup_size(256)
fn copy_invariants(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    if (idx >= params.max_atoms) {
        return;
    }
    set_current_invariant(idx, get_next_invariant(idx));
}

@compute @workgroup_size(256)
fn copy_neighborhoods(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    if (idx >= params.total_neighborhood_words) {
        return;
    }
    set_current_neighborhood(idx, get_next_neighborhood(idx));
}
