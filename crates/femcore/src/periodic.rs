//! Floquet-periodic boundary conditions: the dofs on a `dest` boundary equal the dofs on a `source`
//! boundary (translated by the vector between the two) times `exp(-i phase)`. The two boundaries must
//! carry the same mesh up to the translation.
//!
//! Constraints are collected as relations `x_dest = f * x_source` between dofs, merged into classes
//! with a weighted union-find (so edges and corners shared by several periodic pairs resolve to one
//! free dof per class), and applied as a reduction `x = P x_red`: the reduced matrices are `P^H A P`.

use crate::csparse::{c, CCoo, CCsr, C};
use crate::mesh::SimplexMesh;
use crate::problem::Periodic;
use std::collections::{BTreeSet, HashMap};

/// The Floquet factor `exp(-i phase)`.
pub fn phase_factor(phase: f64) -> C {
    c(phase.cos(), -phase.sin())
}

/// Vertex correspondence of a periodic pair.
pub struct Match {
    /// source vertex -> destination vertex
    pub map: HashMap<u32, u32>,
    /// facets of the source boundary
    pub src_facets: Vec<usize>,
    /// translation from the source to the destination
    pub shift: Vec<f64>,
}

fn group_vertices(mesh: &SimplexMesh, labels: &[i32]) -> (Vec<usize>, Vec<u32>) {
    let mut facets = vec![];
    let mut verts = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        if labels.contains(&mesh.facet_labels[f]) {
            facets.push(f);
            for &v in mesh.facet(f) {
                verts.insert(v);
            }
        }
    }
    (facets, verts.into_iter().collect())
}

/// Matches the vertices of the source boundaries to those of the destination boundaries.
pub fn match_vertices(mesh: &SimplexMesh, p: &Periodic) -> Result<Match, String> {
    let d = mesh.dim;
    let (src_facets, sv) = group_vertices(mesh, &p.source);
    let (_, dv) = group_vertices(mesh, &p.dest);
    if sv.is_empty() || dv.is_empty() {
        return Err("a periodic condition needs boundary numbers on both sides that exist in the mesh".into());
    }
    if sv.len() != dv.len() {
        return Err(format!("the periodic boundaries have different numbers of vertices ({} and {}): their meshes must match", sv.len(), dv.len()));
    }
    let centroid = |vs: &[u32]| -> Vec<f64> {
        let mut c0 = vec![0.0; d];
        for &v in vs {
            for k in 0..d {
                c0[k] += mesh.point(v as usize)[k] / vs.len() as f64;
            }
        }
        c0
    };
    let (cs, cd) = (centroid(&sv), centroid(&dv));
    let shift: Vec<f64> = (0..d).map(|k| cd[k] - cs[k]).collect();
    let (lo, hi) = mesh.bbox();
    let extent = (0..d).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-300);
    let tol = 1e-6 * extent;
    let cell = 4.0 * tol;
    let key = |p: &[f64]| -> Vec<i64> { p.iter().map(|x| (x / cell).floor() as i64).collect() };
    let mut grid: HashMap<Vec<i64>, Vec<u32>> = HashMap::new();
    for &v in &dv {
        grid.entry(key(mesh.point(v as usize))).or_default().push(v);
    }
    let mut map = HashMap::new();
    for &s in &sv {
        let target: Vec<f64> = (0..d).map(|k| mesh.point(s as usize)[k] + shift[k]).collect();
        let kt = key(&target);
        let mut found = None;
        // neighbouring grid cells in every direction
        let mut offsets: Vec<Vec<i64>> = vec![vec![]];
        for _ in 0..d {
            offsets = offsets.into_iter().flat_map(|o| (-1..=1).map(move |dk| { let mut n = o.clone(); n.push(dk); n })).collect();
        }
        for off in offsets {
            let kk: Vec<i64> = kt.iter().zip(&off).map(|(a, b)| a + b).collect();
            if let Some(l) = grid.get(&kk) {
                for &t in l {
                    let pt = mesh.point(t as usize);
                    if (0..d).all(|k| (pt[k] - target[k]).abs() <= tol) {
                        found = Some(t);
                    }
                }
            }
        }
        match found {
            Some(t) => {
                map.insert(s, t);
            }
            None => return Err(format!("no vertex on the destination boundary matches source vertex {s} shifted by the translation: the meshes of the periodic boundaries must be identical")),
        }
    }
    Ok(Match { map, src_facets, shift })
}

fn try_edge(mesh: &SimplexMesh, a: u32, b: u32) -> Option<usize> {
    mesh.edges().binary_search(&(a.min(b), a.max(b))).ok()
}

/// Relations for nodal fields of `ncomp` component blocks (dof layout: component blocks of
/// `n_dofs(order)`): vertex dofs and, for order 2, the midside dofs of the boundary edges.
pub fn nodal_relations(mesh: &SimplexMesh, order: usize, ncomp: usize, m: &Match, phase: f64) -> Result<Vec<(usize, usize, C)>, String> {
    let f = phase_factor(phase);
    let n1 = mesh.n_dofs(order);
    let nv = mesh.n_points();
    let mut out = vec![];
    for (&s, &t) in &m.map {
        for k in 0..ncomp {
            out.push((t as usize + k * n1, s as usize + k * n1, f));
        }
    }
    if order == 2 {
        let mut done = BTreeSet::new();
        for &fc in &m.src_facets {
            let v = mesh.facet(fc);
            let nvf = v.len();
            for i in 0..nvf {
                for j in i + 1..nvf {
                    let (a, b) = (v[i], v[j]);
                    let key = (a.min(b), a.max(b));
                    if !done.insert(key) {
                        continue;
                    }
                    let es = mesh.edge_index(a, b);
                    let et = try_edge(mesh, m.map[&a], m.map[&b]).ok_or("an edge of the source boundary has no counterpart on the destination boundary")?;
                    for k in 0..ncomp {
                        out.push((nv + et + k * n1, nv + es + k * n1, f));
                    }
                }
            }
        }
    }
    Ok(out)
}

/// Weighted union-find: `x_i = ratio_i * x_parent_i`.
struct Dsu {
    parent: Vec<usize>,
    ratio: Vec<C>,
}

impl Dsu {
    fn find(&mut self, i: usize) -> (usize, C) {
        let mut path = vec![];
        let mut cur = i;
        while self.parent[cur] != cur {
            path.push(cur);
            cur = self.parent[cur];
        }
        let root = cur;
        let mut acc = c(1.0, 0.0);
        for &node in path.iter().rev() {
            acc = self.ratio[node] * acc;
            self.ratio[node] = acc;
            self.parent[node] = root;
        }
        (root, if path.is_empty() { c(1.0, 0.0) } else { self.ratio[i] })
    }
    /// Enforces `x_dest = f * x_src`.
    fn union(&mut self, dest: usize, src: usize, f: C) {
        let (rd, kd) = self.find(dest);
        let (rs, ks) = self.find(src);
        if rd == rs {
            return;
        }
        self.parent[rd] = rs;
        self.ratio[rd] = f * ks / kd;
    }
}

/// `x_full = P x_red`: the reduced dofs are the roots of the constraint classes that contain no
/// fixed dof; dofs of fixed classes are zero.
pub struct Reduction {
    pub n: usize,
    /// full dof of every reduced dof
    pub roots: Vec<usize>,
    /// reduced index of every full dof, or `usize::MAX` if it is zero
    pub index_of: Vec<usize>,
    /// `x_i = ratio_i x_root`
    pub ratio: Vec<C>,
}

impl Reduction {
    pub fn build(n: usize, relations: &[(usize, usize, C)], fixed: &[usize]) -> Reduction {
        let mut dsu = Dsu { parent: (0..n).collect(), ratio: vec![c(1.0, 0.0); n] };
        for &(d, s, f) in relations {
            dsu.union(d, s, f);
        }
        let mut dead_root = vec![false; n];
        for &i in fixed {
            let (r, _) = dsu.find(i);
            dead_root[r] = true;
        }
        let mut index_of = vec![usize::MAX; n];
        let mut roots = vec![];
        let mut ratio = vec![c(1.0, 0.0); n];
        let mut root_index = vec![usize::MAX; n];
        for i in 0..n {
            let (r, k) = dsu.find(i);
            ratio[i] = k;
            if dead_root[r] {
                continue;
            }
            if root_index[r] == usize::MAX {
                root_index[r] = roots.len();
                roots.push(r);
            }
            index_of[i] = root_index[r];
        }
        Reduction { n, roots, index_of, ratio }
    }

    pub fn n_reduced(&self) -> usize {
        self.roots.len()
    }

    /// The prolongation `P` (`n x n_reduced`).
    pub fn prolong(&self) -> CCsr {
        let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
        for i in 0..self.n {
            if self.index_of[i] != usize::MAX {
                indices.push(self.index_of[i]);
                data.push(self.ratio[i]);
            }
            indptr.push(indices.len());
        }
        CCsr { nrows: self.n, ncols: self.n_reduced(), indptr, indices, data }
    }

    pub fn expand(&self, xr: &[C]) -> Vec<C> {
        (0..self.n).map(|i| if self.index_of[i] == usize::MAX { c(0.0, 0.0) } else { self.ratio[i] * xr[self.index_of[i]] }).collect()
    }
}

/// Sparse complex product.
pub fn ccsr_mul(a: &CCsr, b: &CCsr) -> CCsr {
    assert_eq!(a.ncols, b.nrows);
    let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
    let mut mark = vec![usize::MAX; b.ncols];
    let mut acc = vec![c(0.0, 0.0); b.ncols];
    let mut touched: Vec<usize> = vec![];
    for i in 0..a.nrows {
        touched.clear();
        for ka in a.indptr[i]..a.indptr[i + 1] {
            let (k, av) = (a.indices[ka], a.data[ka]);
            for kb in b.indptr[k]..b.indptr[k + 1] {
                let j = b.indices[kb];
                if mark[j] != i {
                    mark[j] = i;
                    acc[j] = c(0.0, 0.0);
                    touched.push(j);
                }
                acc[j] += av * b.data[kb];
            }
        }
        touched.sort_unstable();
        for &j in &touched {
            indices.push(j);
            data.push(acc[j]);
        }
        indptr.push(indices.len());
    }
    CCsr { nrows: a.nrows, ncols: b.ncols, indptr, indices, data }
}

/// Conjugate transpose.
pub fn ccsr_adjoint(a: &CCsr) -> CCsr {
    let t = a.transpose();
    CCsr { data: t.data.iter().map(|z| z.conj()).collect(), ..t }
}

/// `P^H A P`.
pub fn congruence(p: &CCsr, a: &CCsr) -> CCsr {
    ccsr_mul(&ccsr_adjoint(p), &ccsr_mul(a, p))
}

#[allow(dead_code)]
fn _unused(_: &CCoo) {}
