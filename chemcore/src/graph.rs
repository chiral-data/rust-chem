use std::collections::HashMap;

/// A neighbor in the molecular graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Neighbor {
    pub atom_idx: usize,
    pub bond_idx: usize,
}

impl Neighbor {
    pub const fn new(atom_idx: usize, bond_idx: usize) -> Self {
        Neighbor { atom_idx, bond_idx }
    }
}

/// Adjacency list representation of a molecular graph.
#[derive(Debug, Clone)]
pub struct MoleculeGraph {
    num_atoms: usize,
    adjacency: Vec<Vec<Neighbor>>,
    edge_map: HashMap<(usize, usize), usize>,
}

impl MoleculeGraph {
    pub fn new(num_atoms: usize) -> Self {
        MoleculeGraph {
            num_atoms,
            adjacency: vec![Vec::new(); num_atoms],
            edge_map: HashMap::new(),
        }
    }

    pub const fn num_atoms(&self) -> usize {
        self.num_atoms
    }

    pub fn num_bonds(&self) -> usize {
        self.edge_map.len() / 2
    }

    pub fn add_edge(&mut self, atom1: usize, atom2: usize, bond_idx: usize) {
        assert!(atom1 < self.num_atoms, "Atom1 index out of bounds");
        assert!(atom2 < self.num_atoms, "Atom2 index out of bounds");
        assert!(atom1 != atom2, "Self-loops not allowed");

        self.adjacency[atom1].push(Neighbor::new(atom2, bond_idx));
        self.adjacency[atom2].push(Neighbor::new(atom1, bond_idx));

        self.edge_map.insert((atom1, atom2), bond_idx);
        self.edge_map.insert((atom2, atom1), bond_idx);
    }

    pub fn neighbors(&self, atom_idx: usize) -> &[Neighbor] {
        assert!(atom_idx < self.num_atoms, "Atom index out of bounds");
        &self.adjacency[atom_idx]
    }

    pub fn degree(&self, atom_idx: usize) -> usize {
        assert!(atom_idx < self.num_atoms, "Atom index out of bounds");
        self.adjacency[atom_idx].len()
    }

    pub fn has_edge(&self, atom1: usize, atom2: usize) -> bool {
        self.edge_map.contains_key(&(atom1, atom2))
    }

    pub fn get_bond(&self, atom1: usize, atom2: usize) -> Option<usize> {
        self.edge_map.get(&(atom1, atom2)).copied()
    }

    pub fn edges(&self) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
        self.edge_map
            .iter()
            .filter(|((a1, a2), _)| a1 < a2)
            .map(|((a1, a2), bond_idx)| (*a1, *a2, *bond_idx))
    }

    pub fn dfs(&self, start: usize) -> Vec<usize> {
        let mut visited = vec![false; self.num_atoms];
        let mut result = Vec::new();
        self.dfs_recursive(start, &mut visited, &mut result);
        result
    }

    fn dfs_recursive(&self, node: usize, visited: &mut [bool], result: &mut Vec<usize>) {
        visited[node] = true;
        result.push(node);

        for neighbor in &self.adjacency[node] {
            if !visited[neighbor.atom_idx] {
                self.dfs_recursive(neighbor.atom_idx, visited, result);
            }
        }
    }

    pub fn bfs(&self, start: usize) -> Vec<usize> {
        let mut visited = vec![false; self.num_atoms];
        let mut queue = std::collections::VecDeque::new();
        let mut result = Vec::new();

        visited[start] = true;
        queue.push_back(start);

        while let Some(node) = queue.pop_front() {
            result.push(node);

            for neighbor in &self.adjacency[node] {
                if !visited[neighbor.atom_idx] {
                    visited[neighbor.atom_idx] = true;
                    queue.push_back(neighbor.atom_idx);
                }
            }
        }

        result
    }

    pub fn neighbors_within_radius(&self, atom_idx: usize, radius: u32) -> Vec<(usize, u32)> {
        let mut distances = vec![None; self.num_atoms];
        let mut queue = std::collections::VecDeque::new();
        let mut result = Vec::new();

        distances[atom_idx] = Some(0);
        queue.push_back((atom_idx, 0));

        while let Some((node, dist)) = queue.pop_front() {
            if dist > 0 {
                result.push((node, dist));
            }

            if dist < radius {
                for neighbor in &self.adjacency[node] {
                    if distances[neighbor.atom_idx].is_none() {
                        distances[neighbor.atom_idx] = Some(dist + 1);
                        queue.push_back((neighbor.atom_idx, dist + 1));
                    }
                }
            }
        }

        result
    }

    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.num_atoms];
        let mut components = Vec::new();

        for start in 0..self.num_atoms {
            if !visited[start] {
                let mut component = Vec::new();
                self.dfs_recursive(start, &mut visited, &mut component);
                components.push(component);
            }
        }

        components
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_creation() {
        let graph = MoleculeGraph::new(5);
        assert_eq!(graph.num_atoms(), 5);
        assert_eq!(graph.num_bonds(), 0);
    }

    #[test]
    fn test_add_edge() {
        let mut graph = MoleculeGraph::new(3);
        graph.add_edge(0, 1, 0);
        graph.add_edge(1, 2, 1);
        assert_eq!(graph.num_bonds(), 2);
        assert!(graph.has_edge(0, 1));
        assert!(graph.has_edge(1, 0));
        assert!(!graph.has_edge(0, 2));
    }

    #[test]
    fn test_neighbors() {
        let mut graph = MoleculeGraph::new(4);
        graph.add_edge(0, 1, 0);
        graph.add_edge(0, 2, 1);
        graph.add_edge(0, 3, 2);
        let neighbors = graph.neighbors(0);
        assert_eq!(neighbors.len(), 3);
        assert_eq!(graph.degree(0), 3);
        assert_eq!(graph.degree(1), 1);
    }

    #[test]
    fn test_dfs() {
        let mut graph = MoleculeGraph::new(4);
        graph.add_edge(0, 1, 0);
        graph.add_edge(1, 2, 1);
        graph.add_edge(2, 3, 2);
        let visited = graph.dfs(0);
        assert_eq!(visited.len(), 4);
        assert_eq!(visited[0], 0);
    }

    #[test]
    fn test_bfs() {
        let mut graph = MoleculeGraph::new(4);
        graph.add_edge(0, 1, 0);
        graph.add_edge(0, 2, 1);
        graph.add_edge(1, 3, 2);
        let visited = graph.bfs(0);
        assert_eq!(visited.len(), 4);
        assert_eq!(visited[0], 0);
        assert!(visited.iter().position(|&x| x == 3).unwrap() > 1);
    }

    #[test]
    fn test_neighbors_within_radius() {
        let mut graph = MoleculeGraph::new(5);
        graph.add_edge(0, 1, 0);
        graph.add_edge(1, 2, 1);
        graph.add_edge(2, 3, 2);
        graph.add_edge(3, 4, 3);
        let radius1 = graph.neighbors_within_radius(0, 1);
        assert_eq!(radius1.len(), 1);
        let radius2 = graph.neighbors_within_radius(0, 2);
        assert_eq!(radius2.len(), 2);
        let radius3 = graph.neighbors_within_radius(0, 3);
        assert_eq!(radius3.len(), 3);
    }

    #[test]
    fn test_connected_components() {
        let mut graph = MoleculeGraph::new(6);
        graph.add_edge(0, 1, 0);
        graph.add_edge(1, 2, 1);
        graph.add_edge(3, 4, 2);
        let components = graph.connected_components();
        assert_eq!(components.len(), 3);
    }

    #[test]
    #[should_panic(expected = "Self-loops not allowed")]
    fn test_no_self_loops() {
        let mut graph = MoleculeGraph::new(3);
        graph.add_edge(0, 0, 0);
    }
}
