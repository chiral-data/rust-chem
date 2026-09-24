//! Named groups of atom indices -- a GROMACS index file's contents (#394).
//!
//! Groups are selections, not a partition: an atom may sit in several groups
//! and more than once inside one, a group may be empty, and two groups may
//! share a name. None of that is normalised away here, because `make_ndx`
//! produces every one of those shapes and a reader that tidied them would
//! disagree with the file.

/// One named group. `atoms` are 0-based, like every other atom index in this
/// crate; the 1-based numbering is the file format's concern, not this type's.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexGroup {
    pub name: String,
    pub atoms: Vec<usize>,
}

/// An ordered list of [`IndexGroup`]s, in file order.
///
/// Holds no structure and validates no index against one: an index file
/// names no coordinate file, so [`Self::max_atom`] is the only bound it
/// carries.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IndexGroups {
    groups: Vec<IndexGroup>,
}

impl IndexGroups {
    pub fn new(groups: Vec<IndexGroup>) -> Self {
        Self { groups }
    }

    pub fn groups(&self) -> &[IndexGroup] {
        &self.groups
    }

    pub fn num_groups(&self) -> usize {
        self.groups.len()
    }

    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Every group called `name`, in file order -- more than one when the
    /// file repeats a name, which is legal.
    pub fn named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a IndexGroup> + 'a {
        self.groups.iter().filter(move |g| g.name == name)
    }

    /// The largest atom index in any group, or `None` if every group is empty.
    pub fn max_atom(&self) -> Option<usize> {
        self.groups
            .iter()
            .flat_map(|g| g.atoms.iter().copied())
            .max()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(name: &str, atoms: &[usize]) -> IndexGroup {
        IndexGroup {
            name: name.to_string(),
            atoms: atoms.to_vec(),
        }
    }

    #[test]
    fn test_named_returns_every_duplicate_in_order() {
        let groups = IndexGroups::new(vec![group("A", &[0]), group("B", &[1]), group("A", &[2])]);
        let atoms: Vec<_> = groups.named("A").map(|g| g.atoms.clone()).collect();
        assert_eq!(atoms, vec![vec![0], vec![2]]);
        assert_eq!(groups.named("C").count(), 0);
    }

    #[test]
    fn test_max_atom_spans_groups_and_ignores_empty_ones() {
        let groups = IndexGroups::new(vec![group("A", &[3, 1]), group("E", &[]), group("B", &[7])]);
        assert_eq!(groups.max_atom(), Some(7));
        assert_eq!(IndexGroups::new(vec![group("E", &[])]).max_atom(), None);
        assert_eq!(IndexGroups::default().max_atom(), None);
    }
}
