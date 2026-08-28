// Tanimoto Similarity GPU Shader
// Computes all-pairs similarity between query and target fingerprints

// ============================================================================
// BUFFER BINDINGS
// ============================================================================

@group(0) @binding(0)
var<storage, read> queries: array<u32>;  // Query fingerprints (packed bits)

@group(0) @binding(1)
var<storage, read> targets: array<u32>;  // Target fingerprints (packed bits)

@group(0) @binding(2)
var<storage, read_write> results: array<f32>;  // Similarity matrix

@group(0) @binding(3)
var<uniform> params: TanimotoParams;

struct TanimotoParams {
    num_queries: u32,
    num_targets: u32,
    fp_words: u32,      // Number of u32 words per fingerprint
    _padding: u32,
}

// ============================================================================
// BIT OPERATIONS
// ============================================================================

fn popcount(x: u32) -> u32 {
    var count = x;
    count = (count & 0x55555555u) + ((count >> 1u) & 0x55555555u);
    count = (count & 0x33333333u) + ((count >> 2u) & 0x33333333u);
    count = (count & 0x0F0F0F0Fu) + ((count >> 4u) & 0x0F0F0F0Fu);
    count = (count & 0x00FF00FFu) + ((count >> 8u) & 0x00FF00FFu);
    count = (count & 0x0000FFFFu) + ((count >> 16u) & 0x0000FFFFu);
    return count;
}

// ============================================================================
// TANIMOTO KERNEL
// ============================================================================

@compute @workgroup_size(16, 16)
fn compute_tanimoto(
    @builtin(global_invocation_id) global_id: vec3<u32>
) {
    let query_idx = global_id.x;
    let target_idx = global_id.y;
    
    // Bounds check
    if (query_idx >= params.num_queries || target_idx >= params.num_targets) {
        return;
    }
    
    // Compute offsets into fingerprint arrays
    let query_offset = query_idx * params.fp_words;
    let target_offset = target_idx * params.fp_words;
    
    var intersection: u32 = 0u;
    var union_count: u32 = 0u;
    
    // Process all words in fingerprint
    for (var i = 0u; i < params.fp_words; i = i + 1u) {
        let q_word = queries[query_offset + i];
        let t_word = targets[target_offset + i];
        
        // Intersection: bits set in both
        let and_word = q_word & t_word;
        intersection = intersection + popcount(and_word);
        
        // Union: bits set in either
        let or_word = q_word | t_word;
        union_count = union_count + popcount(or_word);
    }
    
    // Compute Tanimoto coefficient
    var similarity: f32 = 0.0;
    if (union_count > 0u) {
        similarity = f32(intersection) / f32(union_count);
    }
    
    // Store result in row-major order
    let result_idx = query_idx * params.num_targets + target_idx;
    results[result_idx] = similarity;
}

// ============================================================================
// OPTIMIZED SINGLE-QUERY KERNEL (for 1-vs-N searches)
// ============================================================================

@compute @workgroup_size(256)
fn compute_tanimoto_single_query(
    @builtin(global_invocation_id) global_id: vec3<u32>
) {
    let target_idx = global_id.x;
    
    if (target_idx >= params.num_targets) {
        return;
    }
    
    // Query is always index 0
    let query_offset = 0u;
    let target_offset = target_idx * params.fp_words;
    
    var intersection: u32 = 0u;
    var union_count: u32 = 0u;
    
    for (var i = 0u; i < params.fp_words; i = i + 1u) {
        let q_word = queries[query_offset + i];
        let t_word = targets[target_offset + i];
        
        intersection = intersection + popcount(q_word & t_word);
        union_count = union_count + popcount(q_word | t_word);
    }
    
    var similarity: f32 = 0.0;
    if (union_count > 0u) {
        similarity = f32(intersection) / f32(union_count);
    }
    
    results[target_idx] = similarity;
}
