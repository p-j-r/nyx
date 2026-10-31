//! Fill-reducing ordering by **graph nested dissection** (port item S2), for the sparse Cholesky of
//! [`crate::solvers::direct::ShiftedFactorization`].
//!
//! Recursive bisection with a level-set separator: BFS from a pseudo-peripheral node, the level where the
//! cumulative node count reaches half is the separator, both sides are ordered recursively, the separator
//! is numbered last; subgraphs below `leaf` nodes are numbered in BFS order. This is the simplest member of the
//! family (no multilevel refinement as in METIS). Measured against faer's AMD (`examples/ordering_study.rs`,
//! `docs/LEARNINGS.md`): on 3-D grids 30 % fewer factor non-zeros and a 2.6× faster numeric factorisation
//! (sequential); on 2-D grids *more* non-zeros than AMD, which is why
//! [`FillOrdering::Auto`] compares the two symbolic factorisations and keeps the smaller one.

use std::collections::VecDeque;

/// Which fill-reducing ordering the Cholesky path uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FillOrdering {
    /// faer's approximate minimum degree.
    Amd,
    NestedDissection,
    /// Run both symbolic analyses once and keep the one with fewer factor non-zeros (default).
    #[default]
    Auto,
}

struct Scratch {
    member: Vec<u32>,
    seen: Vec<u32>,
    level: Vec<u32>,
    stamp: u32,
    next_id: u32,
}

fn bfs_levels(adj: &[Vec<u32>], sc: &mut Scratch, id: u32, start: usize) -> Vec<Vec<usize>> {
    sc.stamp += 1;
    let stamp = sc.stamp;
    let mut levels: Vec<Vec<usize>> = vec![vec![start]];
    sc.seen[start] = stamp;
    let mut q = VecDeque::new();
    q.push_back(start);
    sc.level[start] = 0;
    while let Some(u) = q.pop_front() {
        let l = sc.level[u] as usize;
        for &v in &adj[u] {
            let v = v as usize;
            if sc.member[v] == id && sc.seen[v] != stamp {
                sc.seen[v] = stamp;
                sc.level[v] = (l + 1) as u32;
                if levels.len() <= l + 1 {
                    levels.push(Vec::new());
                }
                levels[l + 1].push(v);
                q.push_back(v);
            }
        }
    }
    levels
}

/// A node of the dissection (assembly) tree: the unknowns `order[start..end]` are eliminated together (a
/// separator, or a leaf subdomain), after all unknowns of the `children`.
#[derive(Debug, Clone)]
pub struct NdNode {
    pub start: usize,
    pub end: usize,
    pub children: Vec<usize>,
}

/// Elimination order together with the dissection tree. Nodes are numbered children-first (a post-order), so
/// processing them in increasing number is a valid elimination sequence; `roots` lists the independent
/// components (more than one if the graph is disconnected).
#[derive(Debug, Clone)]
pub struct NdTree {
    pub order: Vec<usize>,
    pub nodes: Vec<NdNode>,
    pub roots: Vec<usize>,
}

fn dissect(
    adj: &[Vec<u32>],
    nodes: Vec<usize>,
    leaf: usize,
    sc: &mut Scratch,
    out: &mut Vec<usize>,
    tree: &mut Vec<NdNode>,
) -> Vec<usize> {
    let leaf_node = |nodes: Vec<usize>, out: &mut Vec<usize>, tree: &mut Vec<NdNode>| -> Vec<usize> {
        let start = out.len();
        out.extend(nodes);
        tree.push(NdNode { start, end: out.len(), children: Vec::new() });
        vec![tree.len() - 1]
    };
    if nodes.len() <= leaf {
        // BFS order inside a leaf keeps neighbours close
        return leaf_node(nodes, out, tree);
    }
    sc.next_id += 1;
    let id = sc.next_id;
    for &v in &nodes {
        sc.member[v] = id;
    }
    let first = bfs_levels(adj, sc, id, nodes[0]);
    let far = first.last().and_then(|l| l.first()).copied().unwrap_or(nodes[0]);
    let levels = bfs_levels(adj, sc, id, far);
    let reached: usize = levels.iter().map(Vec::len).sum();
    if reached < nodes.len() {
        // disconnected: order the reached component and the rest separately
        let comp: Vec<usize> = levels.into_iter().flatten().collect();
        for &v in &comp {
            sc.member[v] = 0;
        }
        let rest: Vec<usize> = nodes.into_iter().filter(|&v| sc.member[v] == id).collect();
        let rest_ids: Vec<usize> = rest;
        for &v in &comp {
            sc.member[v] = id; // restored below by the recursive call's own marking
        }
        let mut roots = dissect(adj, comp, leaf, sc, out, tree);
        roots.extend(dissect(adj, rest_ids, leaf, sc, out, tree));
        return roots;
    }
    let mut cum = 0;
    let mut mid = levels.len() / 2;
    for (i, l) in levels.iter().enumerate() {
        cum += l.len();
        if cum * 2 >= nodes.len() {
            mid = i;
            break;
        }
    }
    if levels.len() < 3 {
        return leaf_node(nodes, out, tree);
    }
    let mid = mid.clamp(1, levels.len() - 2);
    let a: Vec<usize> = levels[..mid].iter().flatten().copied().collect();
    let b: Vec<usize> = levels[mid + 1..].iter().flatten().copied().collect();
    let sep = levels[mid].clone();
    let mut children = dissect(adj, a, leaf, sc, out, tree);
    children.extend(dissect(adj, b, leaf, sc, out, tree));
    let start = out.len();
    out.extend(sep);
    tree.push(NdNode { start, end: out.len(), children });
    vec![tree.len() - 1]
}

/// Elimination order (`order[k]` = the node eliminated `k`-th) of the graph given by symmetric adjacency lists.
pub fn nested_dissection(adj: &[Vec<u32>], leaf: usize) -> Vec<usize> {
    nested_dissection_tree(adj, leaf).order
}

/// [`nested_dissection`] together with its assembly tree (for multifrontal factorisations, see [`crate::solvers::blr`]).
pub fn nested_dissection_tree(adj: &[Vec<u32>], leaf: usize) -> NdTree {
    let n = adj.len();
    let mut sc = Scratch { member: vec![0; n], seen: vec![0; n], level: vec![0; n], stamp: 0, next_id: 0 };
    let mut order = Vec::with_capacity(n);
    let mut nodes = Vec::new();
    let roots = dissect(adj, (0..n).collect(), leaf.max(2), &mut sc, &mut order, &mut nodes);
    debug_assert_eq!(order.len(), n);
    NdTree { order, nodes, roots }
}

/// Symmetric adjacency lists from a lower-triangular CSC pattern (diagonal ignored).
pub fn adjacency_from_lower(n: usize, col_ptr: &[u32], row_idx: &[u32]) -> Vec<Vec<u32>> {
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for j in 0..n {
        for &ri in &row_idx[col_ptr[j] as usize..col_ptr[j + 1] as usize] {
            let i = ri as usize;
            if i != j {
                adj[j].push(i as u32);
                adj[i].push(j as u32);
            }
        }
    }
    for l in &mut adj {
        l.sort_unstable();
        l.dedup();
    }
    adj
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid_adj(nx: usize, ny: usize) -> Vec<Vec<u32>> {
        let mut adj = vec![Vec::new(); nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let p = j * nx + i;
                if i + 1 < nx {
                    adj[p].push((p + 1) as u32);
                    adj[p + 1].push(p as u32);
                }
                if j + 1 < ny {
                    adj[p].push((p + nx) as u32);
                    adj[p + nx].push(p as u32);
                }
            }
        }
        adj
    }

    #[test]
    fn order_is_a_permutation_also_for_disconnected_graphs() {
        let a = grid_adj(30, 30);
        let o = nested_dissection(&a, 20);
        let mut s = o.clone();
        s.sort_unstable();
        assert_eq!(s, (0..900).collect::<Vec<_>>());
        // two disconnected 15×15 grids
        let g = grid_adj(15, 15);
        let mut adj = g.clone();
        adj.extend(g.iter().map(|l| l.iter().map(|&v| v + 225).collect::<Vec<_>>()));
        let o = nested_dissection(&adj, 10);
        let mut s = o.clone();
        s.sort_unstable();
        assert_eq!(s, (0..450).collect::<Vec<_>>());
    }

    #[test]
    fn the_dissection_tree_partitions_the_order_children_first() {
        let a = grid_adj(30, 30);
        let t = nested_dissection_tree(&a, 20);
        assert_eq!(t.order, nested_dissection(&a, 20));
        // node ranges tile 0..n exactly once
        let mut covered = vec![0u8; 900];
        for node in &t.nodes {
            for c in &mut covered[node.start..node.end] {
                *c += 1;
            }
            // children are numbered before their parent and lie before it in the order
            for &ch in &node.children {
                assert!(t.nodes[ch].end <= node.start);
            }
        }
        assert!(covered.iter().all(|&c| c == 1));
        let child_count: usize = t.nodes.iter().map(|n| n.children.len()).sum();
        assert_eq!(child_count + t.roots.len(), t.nodes.len(), "every node is a child or a root");
        assert_eq!(t.roots.len(), 1);
    }

    #[test]
    fn separators_come_last_in_each_recursion() {
        // on a 1-D chain the separator of the first split is a single middle node and it is eliminated last
        let n = 101usize;
        let adj: Vec<Vec<u32>> = (0..n)
            .map(|i| [i.wrapping_sub(1), i + 1].into_iter().filter(|&v| v < n).map(|v| v as u32).collect())
            .collect();
        let o = nested_dissection(&adj, 4);
        let last = *o.last().unwrap();
        assert!((40..=60).contains(&last), "last eliminated node {last}");
    }
}
