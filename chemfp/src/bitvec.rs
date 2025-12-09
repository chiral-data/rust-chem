//use std::fmt;
//
///// Number of bits in a 64-bit machine word.
//const WORD_BITS: usize = 64;
//
///// Mask used to round up word count when allocating bit storage.
/////
///// The formula `(size + WORD_ROUND_MASK) / WORD_BITS`
///// ensures integer division rounds up without using floats.
///// For 64-bit words, rounding mask = 63.
//const WORD_ROUND_MASK: usize = WORD_BITS - 1;
//
///// A compact bit vector optimized for chemical fingerprints.
/////
///// Internally stores bits packed into `u64` words to achieve:
///// - 8× smaller memory footprint than `Vec<bool>`
///// - Fast bitwise operations (`AND`, `OR`)
///// - Constant-time bit access and mutation
/////
///// This representation mirrors fingerprint storage in toolkits like RDKit and OpenBabel.
//#[derive(Debug, Clone, PartialEq, Eq)]
//pub struct BitVec {
//    bits: Vec<u64>,
//    size: usize,
//}
//
//impl BitVec {
//    /// Creates a new bit vector with space for `size` bits.
//    /// Bits are initialized to `0`.
//    /// ## Parameters
//    ///  - `size`: number of bits stored in the vector
//    ///
//    /// ## Notes
//    /// Bits are packed into 64-bit words for memory efficiency and fast bitwise math.
//    ///
//    pub fn new(size: usize) -> Self {
//        let num_words = (size + WORD_ROUND_MASK).div_ceil(WORD_BITS);
//        BitVec {
//            bits: vec![0u64; num_words],
//            size,
//        }
//    }
//
//    /// Returns the total number of addressable bits in this vector.
//    #[inline]
//    pub const fn size(&self) -> usize {
//        self.size
//    }
//
//    /// Sets the bit at position `pos` to `1`.
//    ///
//    /// ## Safety
//    /// If `pos` is out of bounds, the call is ignored (no panic).
//    /// ## Complexity
//    /// O(1)
//    ///
//    #[inline]
//    pub fn set(&mut self, pos: usize) {
//        if pos < self.size {
//            let word_idx = pos / WORD_BITS;
//            let bit_idx = pos % WORD_BITS;
//            self.bits[word_idx] |= 1u64 << bit_idx;
//        }
//    }
//
//    /// Returns the bit value at position `pos`.
//    /// ## Safety
//    /// If `pos` is out of bounds, returns `false`.
//    ///
//    /// ## Complexity
//    /// O(1)
//    ///
//    #[inline]
//    pub fn get(&self, pos: usize) -> bool {
//        if pos < self.size {
//            let word_idx = pos / WORD_BITS;
//            let bit_idx = pos % WORD_BITS;
//            (self.bits[word_idx] & (1u64 << bit_idx)) != 0
//        } else {
//            false
//        }
//    }
//
//    /// Counts how many bits are set to `1` (population count).
//    #[inline]
//    pub fn count_ones(&self) -> usize {
//        self.bits.iter().map(|w| w.count_ones() as usize).sum()
//    }
//
//    /// Performs bitwise AND between two bit vectors.
//    /// ## Preconditions
//    /// Both bit vectors must have the same length.
//    /// ## Complexity
//    /// O(n_words)
//    ///
//    pub fn and(&self, other: &BitVec) -> BitVec {
//        assert_eq!(self.size, other.size);
//        let mut result = BitVec::new(self.size);
//        for i in 0..self.bits.len() {
//            result.bits[i] = self.bits[i] & other.bits[i];
//        }
//        result
//    }
//
//    /// Performs bitwise AND between two bit vectors.
//    /// ## Preconditions
//    /// Both bit vectors must have the same length.
//    /// ## Complexity
//    /// O(n_words)
//    ///
//    pub fn or(&self, other: &BitVec) -> BitVec {
//        assert_eq!(self.size, other.size);
//        let mut result = BitVec::new(self.size);
//        for i in 0..self.bits.len() {
//            result.bits[i] = self.bits[i] | other.bits[i];
//        }
//        result
//    }
//
//    /// Performs bitwise OR between two bit vectors.
//    /// ## Preconditions
//    /// Both bit vectors must have the same length.
//    ///
//    /// ## Complexity
//    /// O(n_words)
//    ///
//    pub fn as_slice(&self) -> &[u64] {
//        &self.bits
//    }
//}
//
//impl fmt::Display for BitVec {
//    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//        for i in 0..self.size {
//            write!(f, "{}", if self.get(i) { '1' } else { '0' })?;
//        }
//        Ok(())
//    }
//}
//
///// Hash a 64-bit value to a bit position in the fingerprint.
/////
///// The result is always in the range `[0, nbits)`.
/////
///// * ## Intended use
///// Mapping atom-environment hashes to fingerprint positions.
/////
//pub fn hash_to_position(hash: u32, nbits: usize) -> usize {
//    (hash as u64 % nbits as u64) as usize
//}
//
///// Derives multiple fingerprint bit positions from a single 64-bit hash.
/////
///// Uses a deterministic multiplicative hash step to produce different
///// bit positions suitable for folding multiple radii of the Morgan algorithm.
/////
///// ## Parameters
///// - `hash`: base hash value
///// - `nbits`: fingerprint length
///// - `count`: number of positions to generate
/////
///// ## Returns
///// A list of bit positions in the range `[0, nbits)`.
/////
//pub fn multi_hash(hash: u32, nbits: usize, count: usize) -> Vec<usize> {
//    if count == 1 {
//        return vec![hash_to_position(hash, nbits)];
//    }
//
//    let mut positions = Vec::with_capacity(count);
//    let mut h = hash;
//
//    for i in 0..count {
//        positions.push(hash_to_position(h, nbits));
//        // Use 32-bit mixing
//        h = h.wrapping_mul(0x9e3779b9u32).wrapping_add((i + 1) as u32);
//    }
//
//    positions
//}
//
//#[cfg(test)]
//mod tests {
//    use super::*;
//
//    #[test]
//    fn test_bitvec_new() {
//        let bv = BitVec::new(100);
//        assert_eq!(bv.size(), 100);
//        assert_eq!(bv.count_ones(), 0);
//    }
//
//    #[test]
//    fn test_bitvec_set_get() {
//        let mut bv = BitVec::new(100);
//        bv.set(10);
//        bv.set(50);
//        bv.set(99);
//
//        assert!(bv.get(10));
//        assert!(bv.get(50));
//        assert!(bv.get(99));
//        assert!(!bv.get(0));
//        assert!(!bv.get(20));
//    }
//
//    #[test]
//    fn test_bitvec_count() {
//        let mut bv = BitVec::new(100);
//        bv.set(10);
//        bv.set(20);
//        bv.set(30);
//
//        assert_eq!(bv.count_ones(), 3);
//    }
//
//    #[test]
//    fn test_bitvec_and() {
//        let mut bv1 = BitVec::new(100);
//        bv1.set(10);
//        bv1.set(20);
//
//        let mut bv2 = BitVec::new(100);
//        bv2.set(10);
//        bv2.set(30);
//
//        let result = bv1.and(&bv2);
//        assert_eq!(result.count_ones(), 1);
//        assert!(result.get(10));
//    }
//
//    #[test]
//    fn test_bitvec_or() {
//        let mut bv1 = BitVec::new(100);
//        bv1.set(10);
//        bv1.set(20);
//
//        let mut bv2 = BitVec::new(100);
//        bv2.set(10);
//        bv2.set(30);
//
//        let result = bv1.or(&bv2);
//        assert_eq!(result.count_ones(), 3);
//        assert!(result.get(10));
//        assert!(result.get(20));
//        assert!(result.get(30));
//    }
//}

use std::fmt;

const WORD_BITS: usize = 64;
const WORD_ROUND_MASK: usize = WORD_BITS - 1;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitVec {
    bits: Vec<u64>,
    size: usize,
}

impl BitVec {
    pub fn new(size: usize) -> Self {
        let num_words = (size + WORD_ROUND_MASK) / WORD_BITS;
        BitVec {
            bits: vec![0u64; num_words],
            size,
        }
    }

    #[inline]
    pub const fn size(&self) -> usize {
        self.size
    }

    #[inline]
    pub fn set(&mut self, pos: usize) {
        if pos < self.size {
            let word_idx = pos / WORD_BITS;
            let bit_idx = pos % WORD_BITS;
            self.bits[word_idx] |= 1u64 << bit_idx;
        }
    }

    #[inline]
    pub fn get(&self, pos: usize) -> bool {
        if pos < self.size {
            let word_idx = pos / WORD_BITS;
            let bit_idx = pos % WORD_BITS;
            (self.bits[word_idx] & (1u64 << bit_idx)) != 0
        } else {
            false
        }
    }

    #[inline]
    pub fn count_ones(&self) -> usize {
        self.bits.iter().map(|w| w.count_ones() as usize).sum()
    }

    pub fn and(&self, other: &BitVec) -> BitVec {
        assert_eq!(self.size, other.size);
        let mut result = BitVec::new(self.size);
        for i in 0..self.bits.len() {
            result.bits[i] = self.bits[i] & other.bits[i];
        }
        result
    }

    pub fn or(&self, other: &BitVec) -> BitVec {
        assert_eq!(self.size, other.size);
        let mut result = BitVec::new(self.size);
        for i in 0..self.bits.len() {
            result.bits[i] = self.bits[i] | other.bits[i];
        }
        result
    }

    pub fn as_slice(&self) -> &[u64] {
        &self.bits
    }

    pub fn on_bits(&self) -> Vec<usize> {
        let mut positions = Vec::new();
        for i in 0..self.size {
            if self.get(i) {
                positions.push(i);
            }
        }
        positions
    }
}

impl fmt::Display for BitVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.size {
            write!(f, "{}", if self.get(i) { '1' } else { '0' })?;
        }
        Ok(())
    }
}
