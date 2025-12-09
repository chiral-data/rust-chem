#[inline]
pub fn hash_combine<T: std::hash::Hash>(seed: &mut u32, value: &T) {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::Hasher;

    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    let hash = hasher.finish() as u32;

    *seed ^= hash
        .wrapping_add(0x9e3779b9)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
}

/// Hash a pair (for bond type + neighbor invariant)
#[inline]
pub fn hash_pair(seed: &mut u32, first: i32, second: u32) {
    // RDKit hashes std::pair by hashing each element
    hash_combine(seed, &first);
    hash_combine(seed, &second);
}
