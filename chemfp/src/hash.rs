#[inline]
pub fn hash_combine<T: std::hash::Hash>(seed: &mut u32, value: &T) {
    let hash = if std::mem::size_of::<T>() == 4 && std::mem::align_of::<T>() == 4 {
        unsafe { std::ptr::read(value as *const T as *const u32) }
    } else {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::Hasher;
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish() as u32
    };

    // boost::hash_combine (what RDKit uses for std::pair/std::vector hashing) —
    // NOT a bare XOR: XOR-only combine is self-cancelling when two combined
    // values are equal (e.g. two identical neighbor invariants), silently
    // dropping one of them from the resulting hash.
    *seed ^= hash
        .wrapping_add(0x9e3779b9)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
}

/// Hash a pair (for bond type + neighbor invariant)
#[inline]
pub fn hash_pair(seed: &mut u32, first: i32, second: u32) {
    // RDKit hashes std::pair by hashing each element
    let first_u32 = first as u32;
    hash_combine(seed, &first_u32);
    hash_combine(seed, &second);
}
