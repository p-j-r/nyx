// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Floquet-periodic boundary conditions: the dofs on a `dest` boundary equal the dofs on a `source`
//! boundary (translated by the vector between the two) times `exp(-i phase)`. The two boundaries must
//! carry the same mesh up to the translation.
//!
//! Constraints are collected as relations `x_dest = f * x_source` between dofs, merged into classes
//! with a weighted union-find (so edges and corners shared by several periodic pairs resolve to one
//! free dof per class), and applied as a reduction `x = P x_red`: the reduced matrices are `P^H A P`.

use crate::linalg::csr::{c, CCoo, CCsr, Csr, RCsr, C};
use crate::linalg::par as lpar;
use crate::linalg::scalar::Scalar;
use crate::mesh::SimplexMesh;
use crate::props::Periodic;
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
    /// Enforces `x_dest = f * x_src` without checking cycles (legacy `Reduction::build`).
    fn union(&mut self, dest: usize, src: usize, f: C) {
        let (rd, kd) = self.find(dest);
        let (rs, ks) = self.find(src);
        if rd == rs {
            return;
        }
        self.parent[rd] = rs;
        self.ratio[rd] = f * ks / kd;
    }

    /// Enforces `x_dest = f * x_src`, returning an error on inconsistent cycles.
    fn try_union(&mut self, dest: usize, src: usize, f: C) -> Result<(), String> {
        let (rd, kd) = self.find(dest);
        let (rs, ks) = self.find(src);
        if rd == rs {
            let rhs = f * ks;
            let tol = 1e-9 * kd.norm().max(rhs.norm()).max(1.0);
            if (kd - rhs).norm() > tol {
                return Err(format!(
                    "inconsistent constraint cycle between dof {dest} and dof {src}: existing ratio {} vs {}",
                    kd / ks,
                    f
                ));
            }
            return Ok(());
        }
        if kd.norm() <= 1e-30 {
            return Err(format!("degenerate zero constraint ratio on dof {dest}"));
        }
        self.parent[rd] = rs;
        self.ratio[rd] = f * ks / kd;
        Ok(())
    }
}

/// `x_full = P x_red`: the reduced dofs are the roots of the constraint classes that contain no
/// fixed dof; dofs of fixed classes are zero (or lifted inhomogeneous Dirichlet values `g`).
#[derive(Clone, Debug)]
pub struct Reduction {
    pub n: usize,
    /// full dof of every reduced dof
    pub roots: Vec<usize>,
    /// reduced index of every full dof, or `usize::MAX` if it is zero
    pub index_of: Vec<usize>,
    /// `x_i = ratio_i x_root`
    pub ratio: Vec<C>,
}

const CONG_CHUNK: usize = 2048;

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
        Self::finish_from_dsu(n, &mut dsu, &dead_root)
    }

    /// Builds a `Reduction` while rejecting inconsistent periodic/equipotential constraint cycles.
    pub fn try_build(n: usize, relations: &[(usize, usize, C)], fixed: &[usize]) -> Result<Reduction, String> {
        let mut dsu = Dsu { parent: (0..n).collect(), ratio: vec![c(1.0, 0.0); n] };
        for &(d, s, f) in relations {
            dsu.try_union(d, s, f)?;
        }
        let mut dead_root = vec![false; n];
        for &i in fixed {
            let (r, _) = dsu.find(i);
            dead_root[r] = true;
        }
        Ok(Self::finish_from_dsu(n, &mut dsu, &dead_root))
    }

    /// Builds a `Reduction` and inhomogeneous Dirichlet lifting vector `g` (`x = P x_red + g`),
    /// detecting both constraint-cycle conflicts and conflicting fixed Dirichlet values within any class.
    pub fn try_build_with_values(
        n: usize,
        relations: &[(usize, usize, C)],
        fixed_vals: &[(usize, C)],
    ) -> Result<(Reduction, Vec<C>), String> {
        let mut dsu = Dsu { parent: (0..n).collect(), ratio: vec![c(1.0, 0.0); n] };
        for &(d, s, f) in relations {
            dsu.try_union(d, s, f)?;
        }
        let mut root_val: Vec<Option<C>> = vec![None; n];
        for &(i, val) in fixed_vals {
            let (r, k) = dsu.find(i);
            if k.norm() <= 1e-30 {
                return Err(format!("degenerate zero constraint ratio on fixed dof {i}"));
            }
            let v_root = val / k;
            match root_val[r] {
                None => root_val[r] = Some(v_root),
                Some(prev) => {
                    let implied = k * prev;
                    let tol = 1e-9 * val.norm().max(implied.norm()).max(1e-15);
                    if (implied - val).norm() > tol {
                        return Err(format!(
                            "conflicting Dirichlet values on constraint class containing dof {i}: {implied} vs {val}"
                        ));
                    }
                }
            }
        }
        let dead_root: Vec<bool> = root_val.iter().map(Option::is_some).collect();
        let red = Self::finish_from_dsu(n, &mut dsu, &dead_root);
        let mut g = vec![c(0.0, 0.0); n];
        for i in 0..n {
            let (r, _) = dsu.find(i);
            if let Some(vr) = root_val[r] {
                g[i] = red.ratio[i] * vr;
            }
        }
        Ok((red, g))
    }

    /// Real-valued convenience wrapper around [`Self::try_build_with_values`].
    pub fn try_build_real_with_values(
        n: usize,
        relations: &[(usize, usize, f64)],
        fixed_vals: &[(usize, f64)],
    ) -> Result<(Reduction, Vec<f64>), String> {
        let rel_c: Vec<(usize, usize, C)> = relations.iter().map(|&(d, s, f)| (d, s, c(f, 0.0))).collect();
        let fix_c: Vec<(usize, C)> = fixed_vals.iter().map(|&(i, v)| (i, c(v, 0.0))).collect();
        let (red, g_c) = Self::try_build_with_values(n, &rel_c, &fix_c)?;
        Ok((red, g_c.into_iter().map(|z| z.re).collect()))
    }

    fn finish_from_dsu(n: usize, dsu: &mut Dsu, dead_root: &[bool]) -> Reduction {
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
        CCsr::from_raw(self.n, self.n_reduced(), indptr, indices, data)
    }

    /// The real prolongation `P` (`n x n_reduced`); errors if any active class carries a non-real ratio.
    pub fn prolong_real(&self) -> Result<RCsr, String> {
        let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
        for i in 0..self.n {
            if self.index_of[i] != usize::MAX {
                if self.ratio[i].im.abs() > 1e-12 {
                    return Err(format!("complex periodic ratio {} on dof {i} cannot be represented as RCsr", self.ratio[i]));
                }
                indices.push(self.index_of[i]);
                data.push(self.ratio[i].re);
            }
            indptr.push(indices.len());
        }
        Ok(RCsr::from_raw(self.n, self.n_reduced(), indptr, indices, data))
    }

    fn ratio_scalars<T: Scalar>(&self) -> Result<Vec<T>, String> {
        let mut out = Vec::with_capacity(self.n);
        for (i, &r) in self.ratio.iter().enumerate() {
            if !T::IS_COMPLEX && self.index_of[i] != usize::MAX && r.im.abs() > 1e-12 {
                return Err(format!("complex constraint ratio {} on dof {i} requires complex arithmetic", r));
            }
            out.push(T::from_parts(r.re, r.im));
        }
        Ok(out)
    }

    fn build_groups(&self) -> (Vec<usize>, Vec<usize>) {
        let nr = self.n_reduced();
        let mut count = vec![0usize; nr + 1];
        for &idx in &self.index_of {
            if idx != usize::MAX {
                count[idx + 1] += 1;
            }
        }
        for p in 0..nr {
            count[p + 1] += count[p];
        }
        let mut next = count.clone();
        let mut members = vec![0usize; count[nr]];
        for (i, &idx) in self.index_of.iter().enumerate() {
            if idx != usize::MAX {
                members[next[idx]] = i;
                next[idx] += 1;
            }
        }
        (count, members)
    }

    fn congruence_fast_impl<T: Scalar>(&self, a: &Csr<T>, ratio_t: &[T]) -> Csr<T> {
        assert_eq!(a.nrows(), self.n);
        assert_eq!(a.ncols(), self.n);
        let nr = self.n_reduced();
        if nr == 0 {
            return Csr::zeros(0, 0);
        }
        let (grp_ptr, members) = self.build_groups();
        let nchunks = nr.div_ceil(CONG_CHUNK);
        let parallel = lpar::worth_parallel(a.nnz());
        let index_of = &self.index_of;

        type ChunkBlock<T> = (Vec<usize>, Vec<u32>, Vec<T>);
        let blocks: Vec<ChunkBlock<T>> = lpar::chunk_map(nchunks, parallel, |chunk_idx| {
            let p0 = chunk_idx * CONG_CHUNK;
            let p1 = ((chunk_idx + 1) * CONG_CHUNK).min(nr);
            let mut marker = vec![usize::MAX; nr];
            let mut acc = vec![T::ZERO; nr];
            let mut touched: Vec<u32> = Vec::new();
            let mut counts = Vec::with_capacity(p1 - p0);
            let mut out_idx = Vec::new();
            let mut out_val = Vec::new();
            for p in p0..p1 {
                touched.clear();
                for &i in &members[grp_ptr[p]..grp_ptr[p + 1]] {
                    let ri_conj = ratio_t[i].conj();
                    let (cols, vals) = a.row(i);
                    for (&j_u32, &a_ij) in cols.iter().zip(vals) {
                        let j = j_u32 as usize;
                        let q = index_of[j];
                        if q == usize::MAX {
                            continue;
                        }
                        let contrib = ri_conj * a_ij * ratio_t[j];
                        if marker[q] != p {
                            marker[q] = p;
                            acc[q] = contrib;
                            touched.push(q as u32);
                        } else {
                            acc[q] += contrib;
                        }
                    }
                }
                touched.sort_unstable();
                for &q in &touched {
                    out_idx.push(q);
                    out_val.push(acc[q as usize]);
                }
                counts.push(touched.len());
            }
            (counts, out_idx, out_val)
        });

        let total_nnz: usize = blocks.iter().map(|b| b.1.len()).sum();
        let mut indptr = Vec::with_capacity(nr + 1);
        indptr.push(0);
        let mut indices = Vec::with_capacity(total_nnz);
        let mut values = Vec::with_capacity(total_nnz);
        for (counts, idx, val) in blocks {
            for cnt in counts {
                indptr.push(*indptr.last().unwrap() + cnt);
            }
            indices.extend(idx);
            values.extend(val);
        }
        Csr::try_new(nr, nr, indptr, indices, values).expect("valid congruence CSR")
    }

    /// Single-pass `O(nnz)` Galerkin congruence `P^H A P` for complex matrices.
    pub fn congruence_fast(&self, a: &CCsr) -> CCsr {
        self.congruence_fast_impl(a, &self.ratio)
    }

    /// Single-pass `O(nnz)` Galerkin congruence `P^T A P` for real matrices; errors if any active ratio is complex.
    pub fn congruence_real(&self, a: &RCsr) -> Result<RCsr, String> {
        let rt: Vec<f64> = self.ratio_scalars()?;
        Ok(self.congruence_fast_impl(a, &rt))
    }

    /// Restricts a full RHS vector `f` via `f_red = P^H f`.
    pub fn restrict_rhs<T: Scalar>(&self, f: &[T]) -> Vec<T> {
        assert_eq!(f.len(), self.n);
        let nr = self.n_reduced();
        let mut out = vec![T::ZERO; nr];
        for i in 0..self.n {
            let p = self.index_of[i];
            if p != usize::MAX {
                let rc = T::from_parts(self.ratio[i].re, -self.ratio[i].im);
                out[p] += rc * f[i];
            }
        }
        out
    }

    /// Computes the reduced RHS with inhomogeneous Dirichlet lifting `f_red = P^H (f - A g)`.
    pub fn lift_rhs<T: Scalar>(&self, a: &Csr<T>, f: &[T], g: &[T]) -> Vec<T> {
        assert_eq!(a.nrows(), self.n);
        assert_eq!(a.ncols(), self.n);
        assert_eq!(f.len(), self.n);
        assert_eq!(g.len(), self.n);
        if g.iter().all(|&x| x == T::ZERO) {
            return self.restrict_rhs(f);
        }
        let ag = a.matvec(g);
        let eff: Vec<T> = f.iter().zip(&ag).map(|(&fi, &agi)| fi - agi).collect();
        self.restrict_rhs(&eff)
    }

    pub fn expand(&self, xr: &[C]) -> Vec<C> {
        (0..self.n).map(|i| if self.index_of[i] == usize::MAX { c(0.0, 0.0) } else { self.ratio[i] * xr[self.index_of[i]] }).collect()
    }

    /// Expands a reduced solution `xr` back to the full vector `x = P x_red + g`.
    pub fn expand_with<T: Scalar>(&self, xr: &[T], g: &[T]) -> Vec<T> {
        assert_eq!(xr.len(), self.n_reduced());
        assert_eq!(g.len(), self.n);
        (0..self.n)
            .map(|i| {
                let p = self.index_of[i];
                if p == usize::MAX {
                    g[i]
                } else {
                    T::from_parts(self.ratio[i].re, self.ratio[i].im) * xr[p]
                }
            })
            .collect()
    }
}

/// Single-pass `O(nnz)` congruence `P^H A P` from a `Reduction`.
pub fn congruence_fast(red: &Reduction, a: &CCsr) -> CCsr {
    red.congruence_fast(a)
}

/// Single-pass `O(nnz)` real congruence `P^T A P` from a `Reduction`.
pub fn congruence_real(red: &Reduction, a: &RCsr) -> Result<RCsr, String> {
    red.congruence_real(a)
}

/// Restricts a full RHS vector `f` via `P^H f`.
pub fn restrict_rhs<T: Scalar>(red: &Reduction, f: &[T]) -> Vec<T> {
    red.restrict_rhs(f)
}

/// Sparse complex product.
pub fn ccsr_mul(a: &CCsr, b: &CCsr) -> CCsr {
    assert_eq!(a.ncols(), b.nrows());
    let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
    let mut mark = vec![usize::MAX; b.ncols()];
    let mut acc = vec![c(0.0, 0.0); b.ncols()];
    let mut touched: Vec<usize> = vec![];
    for i in 0..a.nrows() {
        touched.clear();
        for ka in a.indptr()[i]..a.indptr()[i + 1] {
            let (k, av) = (a.indices()[ka] as usize, a.values()[ka]);
            for kb in b.indptr()[k]..b.indptr()[k + 1] {
                let j = b.indices()[kb] as usize;
                if mark[j] != i {
                    mark[j] = i;
                    acc[j] = c(0.0, 0.0);
                    touched.push(j);
                }
                acc[j] += av * b.values()[kb];
            }
        }
        touched.sort_unstable();
        for &j in &touched {
            indices.push(j);
            data.push(acc[j]);
        }
        indptr.push(indices.len());
    }
    CCsr::from_raw(a.nrows(), b.ncols(), indptr, indices, data)
}

/// Conjugate transpose.
pub fn ccsr_adjoint(a: &CCsr) -> CCsr {
    let t = a.transpose();
    t.conj()
}

/// `P^H A P`.
pub fn congruence(p: &CCsr, a: &CCsr) -> CCsr {
    ccsr_mul(&ccsr_adjoint(p), &ccsr_mul(a, p))
}

#[allow(dead_code)]
fn _unused(_: &CCoo) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::driver::DirectSolver;

    #[test]
    fn congruence_fast_matches_two_matmul_and_real_lifting_solves_exactly() {
        let a = RCsr::from_triplets(
            5,
            5,
            &[
                (0, 0, 4.0),
                (0, 1, -1.0),
                (1, 0, -1.0),
                (1, 1, 5.0),
                (1, 2, -2.0),
                (2, 1, -2.0),
                (2, 2, 6.0),
                (2, 3, -1.5),
                (3, 2, -1.5),
                (3, 3, 7.0),
                (3, 4, -0.5),
                (4, 3, -0.5),
                (4, 4, 3.0),
            ],
        )
        .unwrap();
        // Complex periodic relation + real equipotential relation
        let ph = phase_factor(0.7);
        let red_c = Reduction::try_build(5, &[(4, 1, ph), (3, 2, c(1.0, 0.0))], &[0]).unwrap();
        let ac = a.to_complex();
        let expected = congruence(&red_c.prolong(), &ac);
        let fast = congruence_fast(&red_c, &ac);
        assert_eq!(fast.nrows(), expected.nrows());
        for i in 0..fast.nrows() {
            for j in 0..fast.ncols() {
                assert!((fast.get(i, j) - expected.get(i, j)).norm() < 1e-13);
            }
        }

        // Real inhomogeneous lifting: x[0] = 2.5, x[3] = x[2] (equipotential), x[4] = 1.0 * x[1]
        let (red_r, g) = Reduction::try_build_real_with_values(5, &[(3, 2, 1.0), (4, 1, 1.0)], &[(0, 2.5)]).unwrap();
        let a_red = congruence_real(&red_r, &a).unwrap();
        let p_real = red_r.prolong_real().unwrap();
        assert_eq!(p_real.ncols(), red_r.n_reduced());
        let f = vec![1.0, -2.0, 3.0, 0.5, -1.5];
        let rhs_red = red_r.lift_rhs(&a, &f, &g);
        let x_red = DirectSolver::new(&a_red, false).unwrap().solve(&rhs_red);
        let x_full = red_r.expand_with(&x_red, &g);
        assert!((x_full[0] - 2.5).abs() < 1e-13);
        assert!((x_full[3] - x_full[2]).abs() < 1e-13);
        assert!((x_full[4] - x_full[1]).abs() < 1e-13);
        let ax = a.matvec(&x_full);
        // Check P^T (A x - f) == 0
        let residual_proj = red_r.restrict_rhs(&(0..5).map(|i| ax[i] - f[i]).collect::<Vec<_>>());
        assert!(residual_proj.iter().all(|v| v.abs() < 1e-12), "residual {residual_proj:?}");
    }

    #[test]
    fn dsu_detects_cycle_and_fixed_dirichlet_collisions() {
        // Conflicting cycle: x1 = x0, x2 = x1, x2 = -x0
        let err_cyc = Reduction::try_build(
            3,
            &[(1, 0, c(1.0, 0.0)), (2, 1, c(1.0, 0.0)), (2, 0, c(-1.0, 0.0))],
            &[],
        )
        .unwrap_err();
        assert!(err_cyc.contains("cycle"), "{err_cyc}");

        // Consistent cycle passes
        assert!(Reduction::try_build(
            3,
            &[(1, 0, c(1.0, 0.0)), (2, 1, c(1.0, 0.0)), (2, 0, c(1.0, 0.0))],
            &[],
        )
        .is_ok());

        // Short circuit: x1 = x0, x0 = 0 V, x1 = 1 V
        let err_fix = Reduction::try_build_real_with_values(2, &[(1, 0, 1.0)], &[(0, 0.0), (1, 1.0)]).unwrap_err();
        assert!(err_fix.contains("conflicting Dirichlet"), "{err_fix}");
    }
}
