//! Enhanced stereochemistry groups — a set of stereocentres whose
//! configuration is asserted together, rather than each independently
//! absolute.
//!
//! [`crate::core::atom::Chirality`] on an [`crate::core::atom::Atom`] stores
//! a tetrahedral marker exactly as written (`@`/`@@`), and that has always
//! been read as an absolute configuration. CXSMILES's enhanced stereo
//! extension (`&1:`, `o1:`, `a:` groups) can instead say a *relative*
//! configuration is known while the absolute one is not — e.g. "these two
//! centres are enantiomeric, but which enantiomer this molecule actually is
//! wasn't determined." A [`StereoGroup`] is a layer on top of `Chirality`
//! for exactly that, not a replacement: the per-atom marker is unchanged,
//! this only says how a set of them should be read together.

/// What a [`StereoGroup`] asserts about the atoms it names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoGroupKind {
    /// These centres' configurations are as absolute as
    /// [`crate::core::atom::Chirality`] alone would suggest. CXSMILES writes
    /// at most one of these (`a:`), grouping every centre not otherwise
    /// claimed by an `And`/`Or` group.
    Absolute,
    /// These centres' configurations are known relative to one another, but
    /// which overall enantiomer this molecule is was not determined (`&n:`).
    And,
    /// These centres' configurations are known to be one of a set of
    /// possibilities, relative to one another (`on:`).
    Or,
}

/// A set of stereocentres, read together per [`StereoGroupKind`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    pub kind: StereoGroupKind,
    /// The group's number as written (`&1`, `o2`, ...). `None` for
    /// [`StereoGroupKind::Absolute`] — CXSMILES's `a:` group isn't numbered.
    pub number: Option<u32>,
    /// Atom indices in this group, in the order written.
    pub atoms: Vec<usize>,
}

impl StereoGroup {
    pub fn new(kind: StereoGroupKind, number: Option<u32>, atoms: Vec<usize>) -> Self {
        Self {
            kind,
            number,
            atoms,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_group_holds_the_atoms_it_was_built_with() {
        let group = StereoGroup::new(StereoGroupKind::And, Some(1), vec![1, 4]);
        assert_eq!(group.kind, StereoGroupKind::And);
        assert_eq!(group.number, Some(1));
        assert_eq!(group.atoms, vec![1, 4]);
    }

    #[test]
    fn test_the_absolute_group_carries_no_number() {
        let group = StereoGroup::new(StereoGroupKind::Absolute, None, vec![0]);
        assert_eq!(group.number, None);
    }
}
