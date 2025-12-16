use crate::errors::FingerprintError;
use bitvec::prelude::BitVec;

/// Computes the Tanimoto (Jaccard) similarity between two binary fingerprints.
///
/// # Definition
/// For bit vectors A and B:
///
/// ```text
/// Tanimoto(A, B) = |A ∩ B| / |A ∪ B|
/// ```
///
/// Where:
/// - `|A ∩ B|` is the number of bits set to 1 in both fingerprints
/// - `|A ∪ B|` is the number of bits set to 1 in either fingerprint
///
/// # Behavior
/// - Returns `0.0` if both fingerprints contain no set bits (undefined union case)
/// - Returns an error if fingerprint sizes differ
///
/// # Use Case
/// The most widely used similarity metric for chemical fingerprints,
/// used in virtual screening, clustering, and nearest neighbor search.
///
/// # Errors
/// - [`FingerprintError::SizeMismatch`] if fingerprints have different bit lengths.
///
/// # Example
/// ```ignore
/// let sim = tanimoto_similarity(&fp1, &fp2)?;
/// ```
pub fn tanimoto_similarity(fp1: &BitVec, fp2: &BitVec) -> Result<f64, FingerprintError> {
    if fp1.len() != fp2.len() {
        return Err(FingerprintError::SizeMismatch(fp1.len(), fp2.len()));
    }

    let intersection = fp1
        .iter()
        .zip(fp2.iter())
        .filter(|(a, b)| **a && **b)
        .count(); // (fp1 & fp2).count_ones(); // Bits active in both fingerprints
    let union = fp1
        .iter()
        .zip(fp2.iter())
        .filter(|(a, b)| **a || **b)
        .count(); // (fp1 | fp2).count_ones(); // Bits active in either fingerprint

    // Avoid division by zero when both fp1 and fp2 are all-zero (no features)
    if union == 0 {
        Ok(0.0)
    } else {
        Ok(intersection as f64 / union as f64)
    }
}

/// Computes the Dice (Sørensen–Dice) similarity between two binary fingerprints.
///
/// # Definition
/// For bit vectors A and B:
///
/// ```text
/// Dice(A, B) = 2 × |A ∩ B| / (|A| + |B|)
/// ```
///
/// Where:
/// - `|A ∩ B|` is the number of bits set in both fingerprints
/// - `|A| + |B|` is the sum of active bits in each fingerprint
///
/// # Comparison to Tanimoto
/// - Dice is more sensitive to small overlaps
/// - Sometimes used in scaffold hops or when fingerprints are sparse
///
/// # Behavior
/// - Returns `0.0` if both fingerprints contain no bits
/// - Returns an error if fingerprint sizes differ
///
/// # Errors
/// - [`FingerprintError::SizeMismatch`] if fingerprints differ in bit length.
pub fn dice_similarity(fp1: &BitVec, fp2: &BitVec) -> Result<f64, FingerprintError> {
    if fp1.len() != fp2.len() {
        return Err(FingerprintError::SizeMismatch(fp1.len(), fp2.len()));
    }

    let intersection = fp1
        .iter()
        .zip(fp2.iter())
        .filter(|(a, b)| **a && **b)
        .count(); // (fp1 & fp2).count_ones();
    let sum = fp1.count_ones() + fp2.count_ones();

    if sum == 0 {
        Ok(0.0)
    } else {
        Ok(2.0 * intersection as f64 / sum as f64)
    }
}

/// Computes Tanimoto similarity between a query fingerprint and multiple target fingerprints.
///
/// # Behavior
/// - Returns a similarity score for each target fingerprint
/// - Propagates any size mismatch errors
///
/// # Performance Notes
/// - `batch_tanimoto` is intentionally sequential; parallelism can be added using `rayon`
///   if the caller wishes:
/// ```ignore
/// use rayon::prelude::*;
/// let scores: Vec<f64> = targets.par_iter()
///     .map(|t| tanimoto_similarity(query, t).unwrap())
///     .collect();
/// ```
///
/// # Returns
/// A vector of similarity scores, one per target.
///
/// # Errors
/// - If any target fingerprint differs in size from the query fingerprint.
pub fn batch_tanimoto(query: &BitVec, targets: &[BitVec]) -> Result<Vec<f64>, FingerprintError> {
    targets
        .iter()
        .map(|target| tanimoto_similarity(query, target))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tanimoto_identical() {
        let mut fp1 = BitVec::repeat(false, 100);
        fp1.set(10, true);
        fp1.set(20, true);
        fp1.set(30, true);

        let similarity = tanimoto_similarity(&fp1, &fp1).unwrap();
        assert!((similarity - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_tanimoto_disjoint() {
        let mut fp1 = BitVec::repeat(false, 100);
        fp1.set(10, true);
        fp1.set(20, true);

        let mut fp2 = BitVec::repeat(false, 100);
        fp2.set(30, true);
        fp2.set(40, true);

        let similarity = tanimoto_similarity(&fp1, &fp2).unwrap();
        assert!((similarity - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_tanimoto_overlap() {
        let mut fp1 = BitVec::repeat(false, 100);
        fp1.set(10, true);
        fp1.set(20, true);
        fp1.set(30, true);

        let mut fp2 = BitVec::repeat(false, 100);
        fp2.set(20, true);
        fp2.set(30, true);
        fp2.set(40, true);

        // Intersection: 2 (bits 20, 30)
        // Union: 4 (bits 10, 20, 30, 40)
        // Tanimoto: 2/4 = 0.5
        let similarity = tanimoto_similarity(&fp1, &fp2).unwrap();
        assert!((similarity - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_dice_similarity() {
        let mut fp1 = BitVec::repeat(false, 100);
        fp1.set(10, true);
        fp1.set(20, true);

        let mut fp2 = BitVec::repeat(false, 100);
        fp2.set(20, true);
        fp2.set(30, true);

        // Intersection: 1, Sum: 4, Dice: 2*1/4 = 0.5
        let similarity = dice_similarity(&fp1, &fp2).unwrap();
        assert!((similarity - 0.5).abs() < 1e-10);
    }
}
