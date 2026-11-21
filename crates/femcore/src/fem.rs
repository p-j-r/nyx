// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Lagrange P1/P2 finite elements on simplices: quadrature, shape functions, scalar matrices, boundary integrals and
//! Dirichlet constraints, and the choice between a sparse factorisation and multigrid for symmetric positive systems.
//!
//! Elements are affine triangles/tetrahedra; shape functions are written in barycentric coordinates
//! (P1: `l_i`; P2: `l_i (2 l_i - 1)` at vertices and `4 l_i l_j` at edge midpoints), the gradients of the barycentric
//! coordinates are constant per cell (`geometry`), and element integrals use symmetric quadrature rules (`quad_rule`,
//! exact to degree 2 or 4). `assemble_parallel` runs an element loop on the thread pool and adds the pieces up;
//! `assemble_laplace`/`assemble_mass` give `int c grad phi_i . grad phi_j` and `int c phi_i phi_j`. `solve_spd` solves
//! `K u = f` with prescribed values by a Cholesky factorisation or, for large 3D problems, by conjugate gradients
//! preconditioned with smoothed-aggregation multigrid (`AmgSystem`, `amg`), falling back to the factorisation if
//! multigrid does not converge. The thresholds that decide are `AMG_THRESHOLD*` and `LOBPCG_MIN_FREE`.

pub use crate::farfield::{
    assemble_boundary_mass_weighted, assemble_farfield_robin, farfield_centroid, farfield_positive_measure,
    FarFieldPermittivity,
};
use crate::linalg::driver::DirectSolver;
use crate::mesh::SimplexMesh;
use crate::linalg::csr::{Coo, RCsr as Csr};
use rayon::prelude::*;
use std::collections::BTreeMap;

// ------------------------------------------------------------------ quadrature

/// Quadrature on the reference simplex: barycentric points and weights summing to 1
/// (multiply by the cell volume).
#[derive(Clone, Debug)]
pub struct Quad {
    pub lam: Vec<Vec<f64>>,
    pub w: Vec<f64>,
}

fn bary(p: &[f64]) -> Vec<f64> {
    let mut v = vec![1.0 - p.iter().sum::<f64>()];
    v.extend_from_slice(p);
    v
}

pub fn quad_rule(dim: usize, degree: usize) -> Quad {
    match (dim, degree <= 2) {
        (2, true) => Quad {
            lam: [[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]].iter().map(|p| bary(p)).collect(),
            w: vec![1.0 / 3.0; 3],
        },
        (2, false) => {
            // degree 4, 6 points (Dunavant)
            let (a1, a2) = (0.445948490915965, 0.091576213509771);
            let (w1, w2) = (0.223381589678011, 0.109951743655322);
            let p = [[a1, a1], [1.0 - 2.0 * a1, a1], [a1, 1.0 - 2.0 * a1], [a2, a2], [1.0 - 2.0 * a2, a2], [a2, 1.0 - 2.0 * a2]];
            Quad { lam: p.iter().map(|q| bary(q)).collect(), w: vec![w1, w1, w1, w2, w2, w2] }
        }
        (3, true) => {
            let (a, b) = (0.1381966011250105, 0.5854101966249685);
            let p = [[a, a, a], [b, a, a], [a, b, a], [a, a, b]];
            Quad { lam: p.iter().map(|q| bary(q)).collect(), w: vec![0.25; 4] }
        }
        (3, false) => {
            // degree 4, 14 points (Walkington / Keast type symmetric rule)
            let (a1, a2) = (0.0927352503108912, 0.3108859192633006);
            let b1 = 0.4544962958743504;
            let (w1, w2, w3) = (0.0122488405193936 * 6.0, 0.0187813209530026 * 6.0, 0.0070910034628469 * 6.0);
            let mut pts: Vec<[f64; 3]> = vec![];
            let mut w = vec![];
            for (a, ww) in [(a1, w1), (a2, w2)] {
                let c = 1.0 - 3.0 * a;
                for p in [[a, a, a], [c, a, a], [a, c, a], [a, a, c]] {
                    pts.push(p);
                    w.push(ww);
                }
            }
            let h = 0.5 - b1;
            for p in [[b1, b1, h], [b1, h, b1], [h, b1, b1], [h, h, b1], [h, b1, h], [b1, h, h]] {
                pts.push(p);
                w.push(w3);
            }
            let s: f64 = w.iter().sum();
            Quad { lam: pts.iter().map(|q| bary(q)).collect(), w: w.iter().map(|x| x / s).collect() }
        }
        _ => panic!("unsupported dimension"),
    }
}

/// Rule on a boundary facet (Gauss on edges in 2D, degree-4 Dunavant on triangles in 3D).
pub fn facet_rule(facet_dim: usize) -> Quad {
    if facet_dim == 1 {
        let r = (3.0f64 / 5.0).sqrt();
        let x = [-r, 0.0, r];
        let w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        Quad {
            lam: x.iter().map(|&xi| vec![1.0 - (xi + 1.0) / 2.0, (xi + 1.0) / 2.0]).collect(),
            w: w.iter().map(|x| x / 2.0).collect(),
        }
    } else {
        quad_rule(2, 4)
    }
}

// ------------------------------------------------------------------ shape functions

/// Shape functions at barycentric point `lam`: values `n[a]` and derivatives with respect to
/// the barycentric coordinates `dn[a * (d+1) + b]`.
pub fn shape(lam: &[f64], order: usize) -> (Vec<f64>, Vec<f64>) {
    let n = lam.len();
    if order == 1 {
        let mut dn = vec![0.0; n * n];
        for i in 0..n {
            dn[i * n + i] = 1.0;
        }
        return (lam.to_vec(), dn);
    }
    let nloc = n + n * (n - 1) / 2;
    let mut v = vec![0.0; nloc];
    let mut dn = vec![0.0; nloc * n];
    for i in 0..n {
        v[i] = lam[i] * (2.0 * lam[i] - 1.0);
        dn[i * n + i] = 4.0 * lam[i] - 1.0;
    }
    let mut k = n;
    for i in 0..n {
        for j in i + 1..n {
            v[k] = 4.0 * lam[i] * lam[j];
            dn[k * n + i] = 4.0 * lam[j];
            dn[k * n + j] = 4.0 * lam[i];
            k += 1;
        }
    }
    (v, dn)
}

/// Quadrature rule with shape functions tabulated at its points.
pub struct RefElement {
    pub dim: usize,
    pub order: usize,
    pub nloc: usize,
    pub w: Vec<f64>,
    /// `n[q][a]`
    pub n: Vec<Vec<f64>>,
    /// `dn[q][a * (dim+1) + b]`
    pub dn: Vec<Vec<f64>>,
}

impl RefElement {
    pub fn new(dim: usize, order: usize, degree: usize) -> Self {
        let q = quad_rule(dim, degree);
        let mut n = vec![];
        let mut dn = vec![];
        for lam in &q.lam {
            let (a, b) = shape(lam, order);
            n.push(a);
            dn.push(b);
        }
        RefElement { dim, order, nloc: n[0].len(), w: q.w, n, dn }
    }
}

// ------------------------------------------------------------------ geometry

/// Per-cell volume and gradients of the barycentric coordinates: `g[c*(d+1)*d + b*d + k]`.
pub struct Geometry {
    pub vol: Vec<f64>,
    pub g: Vec<f64>,
}

pub fn geometry(mesh: &SimplexMesh) -> Geometry {
    let d = mesh.dim;
    let nc = mesh.n_cells();
    let mut vol = vec![0.0; nc];
    let mut g = vec![0.0; nc * (d + 1) * d];
    for c in 0..nc {
        let v = mesh.cell(c);
        let x0 = mesh.point(v[0] as usize);
        // rows of J are the edge vectors
        let mut j = [[0.0f64; 3]; 3];
        for i in 0..d {
            let xi = mesh.point(v[i + 1] as usize);
            for k in 0..d {
                j[i][k] = xi[k] - x0[k];
            }
        }
        let (det, inv) = if d == 2 { inv2(&j) } else { inv3(&j) };
        vol[c] = det.abs() / if d == 2 { 2.0 } else { 6.0 };
        let base = c * (d + 1) * d;
        // grad lambda_i (i >= 1) at component k is Jinv[k][i-1]
        for i in 1..=d {
            for k in 0..d {
                g[base + i * d + k] = inv[k][i - 1];
            }
        }
        for k in 0..d {
            let mut s = 0.0;
            for i in 1..=d {
                s += g[base + i * d + k];
            }
            g[base + k] = -s;
        }
    }
    Geometry { vol, g }
}

fn inv2(m: &[[f64; 3]; 3]) -> (f64, [[f64; 3]; 3]) {
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let mut r = [[0.0; 3]; 3];
    r[0][0] = m[1][1] / det;
    r[0][1] = -m[0][1] / det;
    r[1][0] = -m[1][0] / det;
    r[1][1] = m[0][0] / det;
    (det, r)
}

fn inv3(m: &[[f64; 3]; 3]) -> (f64, [[f64; 3]; 3]) {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2];
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    let mut r = [[0.0; 3]; 3];
    r[0][0] = c00 / det;
    r[1][0] = c01 / det;
    r[2][0] = c02 / det;
    r[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
    r[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
    r[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
    r[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
    r[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
    r[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
    (det, r)
}

// ------------------------------------------------------------------ curved tetrahedra (ADR 82)

/// Gauss-Legendre points and weights on [0, 1] (`n` points, exact to degree `2 n - 1`), by Newton iteration on the Legendre polynomial.
fn gauss_legendre_01(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        out.push((0.5 * (1.0 - x), 1.0 / ((1.0 - x * x) * dp * dp)));
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    out
}

/// A quadrature rule on the reference tetrahedron that is exact to at least `degree`, built by collapsing the unit cube (Duffy
/// transform) with Gauss-Legendre points: `n = (degree + 3) / 2` points per direction, `n^3` points in all, weights summing to 1 like
/// those of `quad_rule`. Used on curved cells, where the integrand is a rational function and a high degree is wanted.
pub fn quad_rule_collapsed(degree: usize) -> Quad {
    let n = (degree + 3).div_ceil(2).max(2);
    let g = gauss_legendre_01(n);
    let (mut lam, mut w) = (vec![], vec![]);
    for &(u, wu) in &g {
        for &(v, wv) in &g {
            for &(t, wt) in &g {
                // x = u, y = v (1 - u), z = t (1 - u)(1 - v); jacobian (1 - u)^2 (1 - v); volume of the reference tetrahedron 1/6
                let (x, y, z) = (u, v * (1.0 - u), t * (1.0 - u) * (1.0 - v));
                lam.push(bary(&[x, y, z]));
                w.push(6.0 * wu * wv * wt * (1.0 - u) * (1.0 - u) * (1.0 - v));
            }
        }
    }
    Quad { lam, w }
}

/// The local vertex pairs of a tetrahedron in the order of its edge numbering (and of the geometry nodes 4 to 9 of [`CurvedTet`]).
pub const TET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// A tetrahedron with quadratic (isoparametric) geometry: the map from the reference cell is the quadratic Lagrange interpolation of ten
/// nodes, the four vertices and the nodes of the six edges in the order of [`TET_EDGES`]. With every edge node at the midpoint of its
/// edge the map is affine and the cell is the straight tetrahedron.
///
/// `jacobian` returns, at a barycentric point, the physical gradients `J^-T grad_ref lambda_a` of the four barycentric functions and the
/// determinant of `J`. Used pointwise in place of the constant gradients of a straight cell, the Whitney functions
/// `lambda_a g_b - lambda_b g_a` are exactly the covariant Piola images of the reference ones, and `2 g_a x g_b` is their curl.
#[derive(Clone, Copy, Debug)]
pub struct CurvedTet {
    pub x: [[f64; 3]; 10],
}

impl CurvedTet {
    pub fn new(x: [[f64; 3]; 10]) -> Self {
        CurvedTet { x }
    }

    /// The straight cell with vertices `v` (edge nodes at the midpoints).
    pub fn straight(v: [[f64; 3]; 4]) -> Self {
        let mut x = [[0.0; 3]; 10];
        x[..4].copy_from_slice(&v);
        for (k, &(i, j)) in TET_EDGES.iter().enumerate() {
            for a in 0..3 {
                x[4 + k][a] = 0.5 * (v[i][a] + v[j][a]);
            }
        }
        CurvedTet { x }
    }

    /// Derivatives of the position with respect to each barycentric coordinate treated as independent: `d[b]` is the sum over nodes of
    /// the node position times the derivative of its shape function with respect to `lambda_b`.
    fn dx_dlam(&self, lam: &[f64; 4]) -> [[f64; 3]; 4] {
        let mut d = [[0.0; 3]; 4];
        for a in 0..4 {
            // vertex function l_a (2 l_a - 1): derivative (4 l_a - 1) with respect to l_a
            let s = 4.0 * lam[a] - 1.0;
            for k in 0..3 {
                d[a][k] += s * self.x[a][k];
            }
        }
        for (e, &(i, j)) in TET_EDGES.iter().enumerate() {
            // edge function 4 l_i l_j
            for k in 0..3 {
                d[i][k] += 4.0 * lam[j] * self.x[4 + e][k];
                d[j][k] += 4.0 * lam[i] * self.x[4 + e][k];
            }
        }
        d
    }

    /// Position at the barycentric point `lam`.
    pub fn position(&self, lam: &[f64; 4]) -> [f64; 3] {
        let mut p = [0.0; 3];
        for a in 0..4 {
            let n = lam[a] * (2.0 * lam[a] - 1.0);
            for k in 0..3 {
                p[k] += n * self.x[a][k];
            }
        }
        for (e, &(i, j)) in TET_EDGES.iter().enumerate() {
            let n = 4.0 * lam[i] * lam[j];
            for k in 0..3 {
                p[k] += n * self.x[4 + e][k];
            }
        }
        p
    }

    /// The Jacobian matrix `J[i][k] = dx_i / dxi_k` with respect to the reference coordinates `(lambda_1, lambda_2, lambda_3)`.
    pub fn jacobian_matrix(&self, lam: &[f64; 4]) -> [[f64; 3]; 3] {
        let d = self.dx_dlam(lam);
        let mut j = [[0.0; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                j[i][k] = d[k + 1][i] - d[0][i];
            }
        }
        j
    }

    /// Physical gradients of the four barycentric functions at `lam`, and the determinant of the Jacobian (negative for an inverted
    /// cell).
    pub fn jacobian(&self, lam: &[f64; 4]) -> ([[f64; 3]; 4], f64) {
        let j = self.jacobian_matrix(lam);
        // rows of `m` are the columns of J, so `inv3` of it is the inverse transposed... : grad lambda_k = row k of J^-1 (k = 1..3)
        let (det, inv) = inv3(&j);
        let mut g = [[0.0; 3]; 4];
        for k in 0..3 {
            for i in 0..3 {
                g[k + 1][i] = inv[k][i];
            }
        }
        for i in 0..3 {
            g[0][i] = -(g[1][i] + g[2][i] + g[3][i]);
        }
        (g, det)
    }

    /// The smallest determinant over the ten nodes of the cell (vertices and edge nodes), a cheap validity probe.
    pub fn min_det_at_nodes(&self) -> f64 {
        let mut lams: Vec<[f64; 4]> = vec![[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]];
        for &(i, j) in TET_EDGES.iter() {
            let mut l = [0.0; 4];
            l[i] = 0.5;
            l[j] = 0.5;
            lams.push(l);
        }
        lams.iter().map(|l| self.jacobian(l).1).fold(f64::INFINITY, f64::min)
    }

    /// Volume of the curved cell: the integral of `|det J|` over the reference tetrahedron, with a rule exact to degree 8.
    pub fn volume(&self) -> f64 {
        let q = quad_rule_collapsed(8);
        q.lam.iter().zip(&q.w).map(|(l, w)| w * self.jacobian(&[l[0], l[1], l[2], l[3]]).1.abs() / 6.0).sum()
    }
}

/// Physical gradients of the shape functions at quadrature point `q` of cell `c`:
/// `grad[a*d + k]`.
#[inline]
pub fn shape_gradients(re: &RefElement, geo: &Geometry, c: usize, q: usize, out: &mut [f64]) {
    let d = re.dim;
    let base = c * (d + 1) * d;
    let dn = &re.dn[q];
    for a in 0..re.nloc {
        for k in 0..d {
            let mut s = 0.0;
            for b in 0..=d {
                s += dn[a * (d + 1) + b] * geo.g[base + b * d + k];
            }
            out[a * d + k] = s;
        }
    }
}

// ------------------------------------------------------------------ parallel assembly helper

/// Runs `f(range, &mut coo)` over fixed cell chunks on all available threads and merges the results
/// in deterministic partition order so the assembled matrix is bit-identical across `NYX_THREADS`.
pub fn assemble_parallel<F>(nc: usize, nrows: usize, ncols: usize, f: F) -> Csr
where
    F: Fn(std::ops::Range<usize>, &mut Coo) + Sync,
{
    let nt = crate::par::threads_up_to(16);
    let nparts = if nc < 2000 { 1 } else { nc.div_ceil(2048).min(16) };
    let chunk = nc.div_ceil(nparts).max(1);
    let _stage = crate::progress::stage(format!("assembling element matrices ({nc} cells, {nt} threads)"));
    let mut tl = std::time::Instant::now();
    // Every partition works through its cells in about 8 sub-ranges; the triplets of a sub-range (about 16 bytes per
    // element-matrix entry, many of them duplicates) are turned into a CSR matrix and added to the partition's running
    // sum before the next sub-range starts. Keeping `nparts` independent of `nt` guarantees bit-exact assembly across
    // different thread counts while keeping memory bounded (~16 partial CSRs at most).
    let build_part = |t: usize| -> Csr {
        let (c0, c1) = ((t * chunk).min(nc), ((t + 1) * chunk).min(nc));
        if c0 >= c1 {
            return Csr::zeros(nrows, ncols);
        }
        let sub = (c1 - c0).div_ceil(8).max(1024);
        let mut acc: Option<Csr> = None;
        let mut s0 = c0;
        loop {
            let s1 = (s0 + sub).min(c1);
            let mut coo = Coo::new(nrows, ncols);
            f(s0..s1, &mut coo);
            let part = coo.to_csr();
            acc = Some(match acc {
                None => part,
                Some(a) => a.add_scaled(&part, 1.0),
            });
            if s1 >= c1 {
                break;
            }
            s0 = s1;
        }
        acc.unwrap()
    };
    let parts: Vec<Csr> = if nt <= 1 || nparts <= 1 {
        (0..nparts).map(build_part).collect()
    } else {
        (0..nparts).into_par_iter().map(build_part).collect()
    };
    crate::progress::lap_since(&mut tl, "    assembly: element loops and per-thread CSR");
    let c = Csr::sum_all(parts);
    crate::progress::lap_since(&mut tl, "    assembly: merge");
    c
}

// ------------------------------------------------------------------ scalar matrices

/// `K_ij = integral coeff grad(phi_i) . grad(phi_j)`; `coeff` holds one value per cell.
pub fn assemble_laplace(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Csr {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * d];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.iter_mut().for_each(|x| *x = 0.0);
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                for a in 0..nloc {
                    for b in 0..nloc {
                        let mut s = 0.0;
                        for k in 0..d {
                            s += grad[a * d + k] * grad[b * d + k];
                        }
                        ke[a * nloc + b] += re.w[q] * s;
                    }
                }
            }
            let sc = geo.vol[c] * coeff[c];
            ke.iter_mut().for_each(|x| *x *= sc);
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

/// `K_ij = integral grad(phi_i) . K grad(phi_j)` with a symmetric tensor `K` per cell (only
/// the leading `dim x dim` block is used).
pub fn assemble_laplace_tensor(mesh: &SimplexMesh, order: usize, k: &[[[f64; 3]; 3]]) -> Csr {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * d];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.iter_mut().for_each(|x| *x = 0.0);
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                for a in 0..nloc {
                    // K grad(phi_a)
                    let mut kg = [0.0; 3];
                    for i in 0..d {
                        for j in 0..d {
                            kg[i] += k[c][i][j] * grad[a * d + j];
                        }
                    }
                    for b in 0..nloc {
                        let mut s = 0.0;
                        for i in 0..d {
                            s += kg[i] * grad[b * d + i];
                        }
                        ke[a * nloc + b] += re.w[q] * s;
                    }
                }
            }
            let sc = geo.vol[c];
            ke.iter_mut().for_each(|x| *x *= sc);
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

/// `M_ij = integral coeff phi_i phi_j`.
pub fn assemble_mass(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Csr {
    let d = mesh.dim;
    let re = RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    let mut me = vec![0.0; nloc * nloc];
    for q in 0..re.w.len() {
        for a in 0..nloc {
            for b in 0..nloc {
                me[a * nloc + b] += re.w[q] * re.n[q][a] * re.n[q][b];
            }
        }
    }
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            let sc = geo.vol[c] * coeff[c];
            for (k, m) in ke.iter_mut().zip(&me) {
                *k = m * sc;
            }
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

/// Integral of `coeff phi_i` over every cell: the load of a per-cell source density.
pub fn assemble_source(mesh: &SimplexMesh, order: usize, coeff: &[f64]) -> Vec<f64> {
    let m = assemble_mass(mesh, order, coeff);
    m.matvec(&vec![1.0; mesh.n_dofs(order)])
}

// ------------------------------------------------------------------ boundaries

pub fn facet_measure(mesh: &SimplexMesh, f: usize) -> f64 {
    let v = mesh.facet(f);
    let (p, q) = (mesh.point(v[0] as usize), mesh.point(v[1] as usize));
    if v.len() == 2 {
        ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2)).sqrt()
    } else {
        let r = mesh.point(v[2] as usize);
        let a = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let b = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
        let cx = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        0.5 * (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt()
    }
}

pub fn facets_with<'a>(mesh: &'a SimplexMesh, labels: &'a [i32]) -> impl Iterator<Item = usize> + 'a {
    (0..mesh.n_facets()).filter(move |&f| labels.contains(&mesh.facet_labels[f]))
}

/// `integral coeff phi_i phi_j` over the facets carrying one of `labels` (convective term).
pub fn assemble_boundary_mass(mesh: &SimplexMesh, labels: &[i32], order: usize, coeff: f64) -> Csr {
    let n = mesh.n_dofs(order);
    let fd = mesh.dim - 1;
    let rule = facet_rule(fd);
    let nloc = if order == 1 { mesh.dim } else { mesh.dim + mesh.dim * (mesh.dim - 1) / 2 };
    let mut me = vec![0.0; nloc * nloc];
    for (q, lam) in rule.lam.iter().enumerate() {
        let (nv, _) = shape(lam, order);
        for a in 0..nloc {
            for b in 0..nloc {
                me[a * nloc + b] += rule.w[q] * nv[a] * nv[b];
            }
        }
    }
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    for f in facets_with(mesh, labels) {
        mesh.facet_dofs(f, order, &mut dofs);
        let s = coeff * facet_measure(mesh, f);
        let ke: Vec<f64> = me.iter().map(|m| m * s).collect();
        coo.push_block(&dofs, &dofs, &ke);
    }
    coo.to_csr()
}

/// Load `integral value phi_i` over the facets carrying one of `labels`.
pub fn assemble_boundary_load(mesh: &SimplexMesh, labels: &[i32], value: f64, order: usize) -> Vec<f64> {
    let mut f = vec![0.0; mesh.n_dofs(order)];
    add_boundary_load(mesh, labels, &|_| value, order, &mut f, 1, 0);
    f
}

/// Adds the boundary load of a per-facet value into component `comp` of a vector of `ncomp`
/// stacked dof blocks.
pub fn add_boundary_load(mesh: &SimplexMesh, labels: &[i32], value: &dyn Fn(usize) -> f64, order: usize, out: &mut [f64], ncomp: usize, comp: usize) {
    let n1 = mesh.n_dofs(order);
    debug_assert_eq!(out.len(), ncomp * n1);
    let rule = facet_rule(mesh.dim - 1);
    let nloc = if order == 1 { mesh.dim } else { mesh.dim + mesh.dim * (mesh.dim - 1) / 2 };
    let mut ri = vec![0.0; nloc];
    for (q, lam) in rule.lam.iter().enumerate() {
        let (nv, _) = shape(lam, order);
        for a in 0..nloc {
            ri[a] += rule.w[q] * nv[a];
        }
    }
    let mut dofs = Vec::new();
    for f in facets_with(mesh, labels) {
        mesh.facet_dofs(f, order, &mut dofs);
        let s = value(f) * facet_measure(mesh, f);
        for (a, &dof) in dofs.iter().enumerate() {
            out[comp * n1 + dof] += s * ri[a];
        }
    }
}

/// Dofs on the boundary facets with any of the given boundary numbers (sorted, unique).
pub fn boundary_dofs(mesh: &SimplexMesh, labels: &[i32], order: usize) -> Vec<usize> {
    let mut set = std::collections::BTreeSet::new();
    let mut dofs = Vec::new();
    for f in facets_with(mesh, labels) {
        mesh.facet_dofs(f, order, &mut dofs);
        set.extend(dofs.iter().copied());
    }
    set.into_iter().collect()
}

/// Fixed dofs and values for `{boundary number: value}`. A dof shared by boundaries with
/// different values gets the mean of the values weighted by the number of facets of each
/// boundary that touch it.
pub fn dirichlet_constraints(mesh: &SimplexMesh, values: &BTreeMap<i32, f64>, order: usize) -> (Vec<usize>, Vec<f64>) {
    let mut acc: BTreeMap<usize, (f64, usize)> = BTreeMap::new();
    let mut dofs = Vec::new();
    for f in 0..mesh.n_facets() {
        if let Some(&v) = values.get(&mesh.facet_labels[f]) {
            mesh.facet_dofs(f, order, &mut dofs);
            for &d in &dofs {
                let e = acc.entry(d).or_insert((0.0, 0));
                e.0 += v;
                e.1 += 1;
            }
        }
    }
    let fixed: Vec<usize> = acc.keys().copied().collect();
    let vals = acc.values().map(|(s, n)| s / *n as f64).collect();
    (fixed, vals)
}

/// Solves `K u = f` with `u[fixed] = values` by elimination.
pub fn solve_dirichlet(k: &Csr, f: &[f64], fixed: &[usize], values: &[f64]) -> Result<Vec<f64>, String> {
    let n = k.nrows();
    let mut is_fixed = vec![false; n];
    for &i in fixed {
        is_fixed[i] = true;
    }
    let free: Vec<usize> = (0..n).filter(|&i| !is_fixed[i]).collect();
    let mut u = vec![0.0; n];
    for (&i, &v) in fixed.iter().zip(values) {
        u[i] = v;
    }
    if free.is_empty() {
        return Ok(u);
    }
    let kff = k.submatrix(&free, &free);
    let kfx = k.submatrix(&free, fixed);
    let kx = kfx.matvec(values);
    let rhs: Vec<f64> = free.iter().enumerate().map(|(a, &i)| f[i] - kx[a]).collect();
    let sol = DirectSolver::new(&kff, true)?.solve(&rhs);
    for (a, &i) in free.iter().enumerate() {
        u[i] = sol[a];
    }
    Ok(u)
}

/// Dofs (of `ncomp` component blocks of `n_dofs(order)`) that no active cell touches: the cells
/// of the `inactive` domains do not take part in the problem, so nodes that only they use are held.
pub fn untouched_dofs(mesh: &SimplexMesh, order: usize, ncomp: usize, inactive: &[i32]) -> Vec<usize> {
    if inactive.is_empty() {
        return vec![];
    }
    let n1 = mesh.n_dofs(order);
    let mut used = vec![false; n1];
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        if !inactive.contains(&mesh.cell_labels[c]) {
            mesh.cell_dofs(c, order, &mut dofs);
            for &d in &dofs {
                used[d] = true;
            }
        }
    }
    let mut out = vec![];
    for k in 0..ncomp {
        out.extend((0..n1).filter(|&d| !used[d]).map(|d| d + k * n1));
    }
    out
}

/// Free dof indices (complement of `fixed`).
pub fn free_dofs(n: usize, fixed: &[usize]) -> Vec<usize> {
    let mut is_fixed = vec![false; n];
    for &i in fixed {
        is_fixed[i] = true;
    }
    (0..n).filter(|&i| !is_fixed[i]).collect()
}




/// Free-dof count above which `LinearSolver::Auto` switches to multigrid in 2D. In 3D multigrid
/// is faster than a sparse factorisation from a few tens of thousands of unknowns (a drilled
/// plate with 190,000 unknowns: 1.7 s against 28 s; 3D blocks of 680,000: 7 s against 340 s).
pub const AMG_THRESHOLD: usize = 400_000;
pub const AMG_THRESHOLD_3D: usize = 25_000;
/// Free dofs above which eigenfrequency studies (Auto) try LOBPCG instead of a factorisation. Measured with the shared
/// pool, 6 modes of a cube, P2, (factorisation / Auto): 47,000 unknowns 4.5 s / 4.9 s, 108,000: 25 s / 8.9 s, 207,000:
/// 93 s / 19.5 s; a slender beam where multigrid stalls, 70,000: 5.3 s / 28 s. Break-even is near 80,000.
pub const LOBPCG_MIN_FREE: usize = 80_000;
/// Above this many free dofs a factorisation is not attempted for eigenfrequencies.
pub const FACTORISE_MAX_FREE: usize = 1_200_000;
/// Multigrid is used for eigenfrequencies when its convergence factor is at most this.
pub const LOBPCG_MAX_FACTOR: f64 = 0.65;

/// Whether `LinearSolver::Auto` should use multigrid for a system with `n_free` unknowns.
pub fn amg_wanted(dim: usize, n_free: usize) -> bool {
    n_free > if dim == 3 { AMG_THRESHOLD_3D } else { AMG_THRESHOLD }
}

/// Near-nullspace of the operator in node-major dof order: constants for a scalar field, the rigid
/// body modes for displacements, with entries at constrained dofs set to zero.
fn near_nullspace(mesh: &SimplexMesh, order: usize, ncomp: usize, is_fixed: &[bool], perm: &[usize]) -> Vec<Vec<f64>> {
    let dim = mesh.dim;
    let n = is_fixed.len();
    let n1 = n / ncomp;
    let coords = mesh.dof_coordinates(order);
    let mut centre = vec![0.0; dim];
    for p in coords.chunks(dim) {
        for a in 0..dim {
            centre[a] += p[a] / n1 as f64;
        }
    }
    let mut ns: Vec<Vec<f64>> = vec![];
    if ncomp == 1 {
        ns.push(vec![1.0; n]);
    } else {
        for t in 0..dim {
            let mut v = vec![0.0; n];
            for node in 0..n1 {
                v[node * ncomp + t] = 1.0;
            }
            ns.push(v);
        }
        let nrot = if dim == 2 { 1 } else { 3 };
        for r in 0..nrot {
            let mut v = vec![0.0; n];
            for node in 0..n1 {
                let p: Vec<f64> = (0..dim).map(|a| coords[node * dim + a] - centre[a]).collect();
                let (x, y, z) = (p[0], p[1], if dim == 3 { p[2] } else { 0.0 });
                let rot: [f64; 3] = if dim == 2 {
                    [-y, x, 0.0]
                } else {
                    match r {
                        0 => [0.0, -z, y],
                        1 => [z, 0.0, -x],
                        _ => [-y, x, 0.0],
                    }
                };
                for c in 0..dim {
                    v[node * ncomp + c] = rot[c];
                }
            }
            ns.push(v);
        }
    }
    for v in ns.iter_mut() {
        for (r, x) in v.iter_mut().enumerate() {
            if is_fixed[perm[r]] {
                *x = 0.0;
            }
        }
    }
    ns
}

/// Auto-detects vertical prismatic columns on a 3D mesh by grouping P1 vertices with identical
/// quantized `(x, y)` coordinates (`1e-9 * extent`) into sorted-`z` columns of length `>= 3`.
/// When columns cover `> 50%` of P1 vertices, returns `(p1_columns, p2_columns)`.
pub fn detect_vertical_columns(mesh: &SimplexMesh, order: usize) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    if mesh.dim != 3 {
        return (Vec::new(), Vec::new());
    }
    let nv = mesh.n_points();
    if nv < 3 {
        return (Vec::new(), Vec::new());
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for v in 0..nv {
        let p = mesh.point(v);
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(hi[2] - lo[2]).max(1e-30);
    let q_step = 1e-9 * extent;

    let mut p1_groups: BTreeMap<(i64, i64), Vec<usize>> = BTreeMap::new();
    for v in 0..nv {
        let p = mesh.point(v);
        let qx = ((p[0] - lo[0]) / q_step).round() as i64;
        let qy = ((p[1] - lo[1]) / q_step).round() as i64;
        p1_groups.entry((qx, qy)).or_default().push(v);
    }
    let mut p1_cols = Vec::new();
    let mut covered_p1 = 0usize;
    for (_, mut col) in p1_groups {
        if col.len() >= 3 {
            col.sort_by(|&a, &b| {
                mesh.point(a)[2]
                    .partial_cmp(&mesh.point(b)[2])
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.cmp(&b))
            });
            covered_p1 += col.len();
            p1_cols.push(col);
        }
    }
    if covered_p1 * 2 <= nv {
        return (Vec::new(), Vec::new());
    }
    if order != 2 {
        return (p1_cols, Vec::new());
    }
    let coords = mesh.dof_coordinates(2);
    let n2 = mesh.n_dofs(2);
    let mut p2_groups: BTreeMap<(i64, i64), Vec<usize>> = BTreeMap::new();
    for d in 0..n2 {
        let p = &coords[d * 3..(d + 1) * 3];
        let qx = ((p[0] - lo[0]) / q_step).round() as i64;
        let qy = ((p[1] - lo[1]) / q_step).round() as i64;
        p2_groups.entry((qx, qy)).or_default().push(d);
    }
    let mut p2_cols = Vec::new();
    for (_, mut col) in p2_groups {
        if col.len() >= 3 {
            col.sort_by(|&a, &b| {
                coords[a * 3 + 2]
                    .partial_cmp(&coords[b * 3 + 2])
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.cmp(&b))
            });
            p2_cols.push(col);
        }
    }
    (p1_cols, p2_cols)
}

/// A symmetric positive definite system over all dofs with some dofs constrained, prepared for
/// multigrid-preconditioned conjugate gradients: constrained rows and columns are replaced by a
/// scaled identity (so the system keeps its scale), the dofs are ordered node by node, and the
/// multigrid hierarchy is built once and can serve many right-hand sides (eigenvalue iterations,
/// time steps) or serve as the preconditioner of nearby matrices.
pub struct AmgSystem {
    n: usize,
    perm: Vec<usize>,
    is_fixed: Vec<bool>,
    free: Vec<usize>,
    inv: Vec<usize>,
    dscale: f64,
    mp: Csr,
    amg: crate::linalg::amg_elastic::Amg,
    pub setup_seconds: f64,
}

impl AmgSystem {
    /// `a` is the matrix over all `n = ncomp * n_dofs(order)` dofs (entries in constrained rows and
    /// columns are ignored); `ncomp` is 1 for a scalar field and `dim` for displacements.
    pub fn new(mesh: &SimplexMesh, order: usize, ncomp: usize, a: &Csr, fixed: &[usize]) -> Result<AmgSystem, String> {
        Self::new_with(mesh, order, ncomp, a, fixed, crate::linalg::amg_elastic::AmgOptions::default())
    }

    /// Like [`AmgSystem::new`] with explicit [`crate::linalg::amg_elastic::AmgOptions`] (such as `line_smooth`).
    pub fn new_with(
        mesh: &SimplexMesh,
        order: usize,
        ncomp: usize,
        a: &Csr,
        fixed: &[usize],
        mut opts: crate::linalg::amg_elastic::AmgOptions,
    ) -> Result<AmgSystem, String> {
        let t_setup = std::time::Instant::now();
        let n = a.nrows();
        let n1 = n / ncomp;
        let mut is_fixed = vec![false; n];
        for &i in fixed {
            is_fixed[i] = true;
        }
        let diag = a.diagonal();
        let dscale = {
            let (mut s, mut c) = (0.0, 0usize);
            for i in 0..n {
                if !is_fixed[i] {
                    s += diag[i].abs();
                    c += 1;
                }
            }
            if c > 0 && s > 0.0 { s / c as f64 } else { 1.0 }
        };
        // node-major ordering (all components of a node together); the constrained copy is only an intermediate
        // (as large as the matrix itself), so it is freed as soon as the permuted one exists
        let perm: Vec<usize> = (0..n).map(|r| (r % ncomp) * n1 + r / ncomp).collect();
        let mp = {
            let m = Self::constrained(a, &is_fixed, dscale);
            m.submatrix(&perm, &perm)
        };
        let ns = near_nullspace(mesh, order, ncomp, &is_fixed, &perm);
        if opts.line_smooth && ncomp == 1 && mesh.dim == 3 && opts.p1_columns.is_empty() {
            let (p1_cols, p2_cols) = detect_vertical_columns(mesh, order);
            opts.p1_columns = p1_cols;
            opts.p2_columns = p2_cols;
        }
        let amg = if order == 2 && !opts.no_p1_coarse {
            // linear interpolation of P1 into P2: vertex nodes are the identity, an edge node the mean of its two ends
            let nv = mesh.n_points();
            let mut indptr = vec![0usize];
            let (mut indices, mut data) = (vec![], vec![]);
            for node in 0..n1 {
                if node < nv {
                    indices.push(node);
                    data.push(1.0);
                } else {
                    let (p, q) = mesh.edges()[node - nv];
                    indices.extend_from_slice(&[p as usize, q as usize]);
                    data.extend_from_slice(&[0.5, 0.5]);
                }
                indptr.push(indices.len());
            }
            let lin = Csr::from_raw(n1, nv, indptr, indices, data);
            crate::linalg::amg_elastic::Amg::new_p2(&mp, ns, ncomp, &opts, &lin)?
        } else {
            crate::linalg::amg_elastic::Amg::new(&mp, ns, ncomp, &opts)?
        };
        let free = free_dofs(n, fixed);
        let mut inv = vec![0usize; n];
        for (pos, &full) in perm.iter().enumerate() {
            inv[full] = pos;
        }
        Ok(AmgSystem { n, perm, is_fixed, free, inv, dscale, mp, amg, setup_seconds: t_setup.elapsed().as_secs_f64() })
    }

    fn constrained(a: &Csr, is_fixed: &[bool], dscale: f64) -> Csr {
        let n = a.nrows();
        let (mut indptr, mut indices, mut data) = (vec![0usize], Vec::<usize>::new(), Vec::<f64>::new());
        for r in 0..n {
            if is_fixed[r] {
                indices.push(r);
                data.push(dscale);
            } else {
                for p in a.indptr()[r]..a.indptr()[r + 1] {
                    if !is_fixed[a.indices()[p] as usize] {
                        indices.push(a.indices()[p] as usize);
                        data.push(a.values()[p]);
                    }
                }
            }
            indptr.push(indices.len());
        }
        Csr::from_raw(n, n, indptr, indices, data)
    }

    pub fn n_levels(&self) -> usize {
        self.amg.n_levels()
    }

    pub fn sizes(&self) -> Vec<usize> {
        self.amg.sizes()
    }

    /// Whether Thomas tridiagonal vertical-line smoothing is active in the hierarchy.
    pub fn line_smooth_active(&self) -> bool {
        self.amg.line_smooth_active()
    }

    /// Number of vertical columns detected and smoothed on the primary P1/P2 levels.
    pub fn n_line_columns(&self) -> usize {
        self.amg.n_line_columns()
    }

    /// Solves `A x = b` for the right-hand side over all dofs (entries at constrained dofs are taken
    /// as the prescribed values times the scale of the matrix); `matrix` is the system to solve,
    /// which may differ from the one the hierarchy was built from (it is then only the preconditioner).
    pub fn solve_with(&self, matrix: Option<&Csr>, b: &[f64], fixed_values: &[(usize, f64)], tol: f64, max_iter: usize) -> (Vec<f64>, usize, f64) {
        let mut rhs = b.to_vec();
        for i in 0..self.n {
            if self.is_fixed[i] {
                rhs[i] = 0.0;
            }
        }
        for &(i, v) in fixed_values {
            rhs[i] = self.dscale * v;
        }
        let rp: Vec<f64> = self.perm.iter().map(|&i| rhs[i]).collect();
        let mp_other;
        let a = match matrix {
            Some(m) => {
                mp_other = Self::constrained(m, &self.is_fixed, self.dscale).submatrix(&self.perm, &self.perm);
                &mp_other
            }
            None => &self.mp,
        };
        let res = crate::linalg::amg_elastic::pcg(a, &rp, &self.amg, tol, max_iter);
        let mut x = vec![0.0; self.n];
        for r in 0..self.n {
            x[self.perm[r]] = res.x[r];
        }
        (x, res.iterations, res.relative_residual)
    }

    /// Solves `A X = B` simultaneously for `K = b_list.len()` right-hand sides using multi-RHS
    /// Block-CG (`block_pcg`), traversing `A` and the AMG V-cycle once per iteration across all `K`
    /// vectors. Returns `(solutions, iterations, max_relative_residual)`.
    pub fn solve_block_with(
        &self,
        matrix: Option<&Csr>,
        b_list: &[Vec<f64>],
        fixed_values_list: &[Vec<(usize, f64)>],
        tol: f64,
        max_iter: usize,
    ) -> (Vec<Vec<f64>>, usize, f64) {
        let k = b_list.len();
        assert_eq!(k, fixed_values_list.len());
        if k == 0 {
            return (Vec::new(), 0, 0.0);
        }
        let mut rp = vec![0.0; self.n * k];
        let mut rhs = vec![0.0; self.n];
        for j in 0..k {
            rhs.copy_from_slice(&b_list[j]);
            for i in 0..self.n {
                if self.is_fixed[i] {
                    rhs[i] = 0.0;
                }
            }
            for &(i, v) in &fixed_values_list[j] {
                rhs[i] = self.dscale * v;
            }
            for (r, &full) in self.perm.iter().enumerate() {
                rp[r * k + j] = rhs[full];
            }
        }
        let mp_other;
        let a = match matrix {
            Some(m) => {
                mp_other = Self::constrained(m, &self.is_fixed, self.dscale).submatrix(&self.perm, &self.perm);
                &mp_other
            }
            None => &self.mp,
        };
        let res = crate::linalg::amg_elastic::block_pcg(a, &rp, k, &self.amg, tol, max_iter);
        let mut out = vec![vec![0.0; self.n]; k];
        for (r, &full) in self.perm.iter().enumerate() {
            for j in 0..k {
                out[j][full] = res.x[r * k + j];
            }
        }
        (out, res.iterations, res.relative_residual)
    }

    /// Average residual reduction per iteration of multigrid-preconditioned CG on a random
    /// right-hand side over a few iterations (smaller is better; 0.4 means about 15 iterations
    /// per six digits, 0.9 means about 130). Used to choose between multigrid and a factorisation.
    pub fn convergence_factor(&self) -> f64 {
        let mut seed = 0x2545F4914F6CDD1Du64;
        let b: Vec<f64> = (0..self.n)
            .map(|i| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                if self.is_fixed[i] { 0.0 } else { ((seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5 }
            })
            .collect();
        let rp: Vec<f64> = self.perm.iter().map(|&i| b[i]).collect();
        let its = 8;
        let res = crate::linalg::amg_elastic::pcg(&self.mp, &rp, &self.amg, 1e-300, its);
        res.relative_residual.max(1e-300).powf(1.0 / its as f64)
    }

    /// One multigrid cycle applied to a vector over the free dofs: an approximation of the inverse
    /// of the matrix, for use as a preconditioner.
    pub fn precondition_free(&self, r_free: &[f64]) -> Vec<f64> {
        let mut b = vec![0.0; self.n];
        for (a, &i) in self.free.iter().enumerate() {
            b[self.perm_inv_pos(i)] = r_free[a];
        }
        let z = self.amg.apply(&b);
        self.free.iter().map(|&i| z[self.perm_inv_pos(i)]).collect()
    }

    /// Position of full dof `i` in the node-major ordering.
    fn perm_inv_pos(&self, i: usize) -> usize {
        self.inv[i]
    }

    /// Solves for the unconstrained dofs only: `b` and the result are vectors over the free dofs.
    pub fn solve_free(&self, matrix: Option<&Csr>, b_free: &[f64], tol: f64) -> Result<Vec<f64>, String> {
        let mut b = vec![0.0; self.n];
        for (a, &i) in self.free.iter().enumerate() {
            b[i] = b_free[a];
        }
        let (x, its, rel) = self.solve_with(matrix, &b, &[], tol, 800);
        crate::progress::check()?;
        if rel > 1e-7 {
            return Err(format!("the iterative solver did not converge ({its} iterations, residual {rel:.2e}); choose the direct solver"));
        }
        Ok(self.free.iter().map(|&i| x[i]).collect())
    }
}

/// Solves `K u = f` with `u[fixed] = values` choosing between the sparse direct solver and
/// AMG-preconditioned conjugate gradients. `ncomp` is 1 for a scalar field and `dim` for
/// displacements (dof layout: component blocks of `n_dofs(order)` entries). Returns the solution
/// and a note about the method used when the iterative solver ran. With `LinearSolver::Auto` a
/// multigrid run that does not converge falls back to the direct solver.
#[allow(clippy::too_many_arguments)]
pub fn solve_spd(mesh: &SimplexMesh, order: usize, ncomp: usize, solver: crate::props::LinearSolver, k: &Csr, f: &[f64], fixed: &[usize], values: &[f64]) -> Result<(Vec<f64>, Option<String>), String> {
    solve_spd_with(mesh, order, ncomp, solver, k, f, fixed, values, 1e-10, true)
}

/// Like [`solve_spd`] with explicit iterative tolerance `tol` and vertical-line smoother toggle `line_smooth`.
#[allow(clippy::too_many_arguments)]
pub fn solve_spd_with(
    mesh: &SimplexMesh,
    order: usize,
    ncomp: usize,
    solver: crate::props::LinearSolver,
    k: &Csr,
    f: &[f64],
    fixed: &[usize],
    values: &[f64],
    tol: f64,
    line_smooth: bool,
) -> Result<(Vec<f64>, Option<String>), String> {
    use crate::props::LinearSolver;
    let n_free = k.nrows() - fixed.len();
    let iterative = match solver {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => true,
        LinearSolver::Auto => amg_wanted(mesh.dim, n_free),
    };
    if !iterative || fixed.is_empty() {
        return Ok((solve_dirichlet(k, f, fixed, values)?, None));
    }
    let mut tl = std::time::Instant::now();
    let n = k.nrows();
    // right-hand side with the constraints moved over
    let mut uvec = vec![0.0; n];
    for (&i, &v) in fixed.iter().zip(values) {
        uvec[i] = v;
    }
    let ku = k.matvec(&uvec);
    let rhs: Vec<f64> = (0..n).map(|i| f[i] - ku[i]).collect();
    let opts = crate::linalg::amg_elastic::AmgOptions {
        line_smooth,
        ..Default::default()
    };
    let sys = AmgSystem::new_with(mesh, order, ncomp, k, fixed, opts)?;
    crate::progress::lap_since(&mut tl, "    spd: system and hierarchy");
    let t_solve = std::time::Instant::now();
    let fv: Vec<(usize, f64)> = fixed.iter().copied().zip(values.iter().copied()).collect();
    let (mut u, its, rel) = sys.solve_with(None, &rhs, &fv, tol, 600);
    let solve_s = t_solve.elapsed().as_secs_f64();
    crate::progress::check()?;
    if rel > (tol * 1e3).max(1e-7) {
        if solver == LinearSolver::Auto {
            let u = solve_dirichlet(k, f, fixed, values)?;
            return Ok((u, Some(format!("multigrid did not converge ({its} iterations, residual {rel:.1e}); used the direct solver"))));
        }
        return Err(format!("the iterative solver did not converge ({its} iterations, residual {rel:.2e}); choose the direct solver"));
    }
    for (&i, &v) in fixed.iter().zip(values) {
        u[i] = v;
    }
    let note = format!("AMG-CG: {} levels {:?}, {} iterations, residual {:.1e} (setup {:.1} s, solve {:.1} s)", sys.n_levels(), sys.sizes(), its, rel, sys.setup_seconds, solve_s);
    Ok((u, Some(note)))
}

/// Refuses boundary numbers that the mesh does not have, naming the ones it has. Every physics checks its conditions with it.
pub fn check_boundary_labels(mesh: &SimplexMesh, labels: impl IntoIterator<Item = i32>) -> Result<(), String> {
    let known = mesh.boundary_labels();
    let mut bad: Vec<i32> = labels.into_iter().filter(|l| !known.contains(l)).collect();
    bad.sort();
    bad.dedup();
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!("boundary number(s) {bad:?} do not exist in this mesh (it has {known:?})"))
    }
}

/// The values of a dof field at the local dofs of every cell: `out[c*nloc + a]`.
pub fn cell_values_at_nodes(mesh: &SimplexMesh, u: &[f64], order: usize) -> Vec<f64> {
    let nloc = mesh.nloc(order);
    let mut out = Vec::with_capacity(mesh.n_cells() * nloc);
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        out.extend(dofs.iter().take(nloc).map(|&d| u[d]));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadrature_integrates_monomials() {
        // exact integral of l0^a l1^b l2^c over a triangle (weights sum to 1): 2! a! b! c! / (a+b+c+2)!
        fn fact(n: usize) -> f64 {
            (1..=n).map(|x| x as f64).product()
        }
        for (dim, deg) in [(2, 2), (2, 4), (3, 2), (3, 4)] {
            let q = quad_rule(dim, deg);
            for a in 0..=deg {
                for b in 0..=deg - a {
                    let c = deg - a - b;
                    let val: f64 = q.lam.iter().zip(&q.w).map(|(l, w)| w * l[0].powi(a as i32) * l[1].powi(b as i32) * l[2].powi(c as i32)).sum();
                    let exact = if dim == 2 {
                        2.0 * fact(a) * fact(b) * fact(c) / fact(a + b + c + 2)
                    } else {
                        6.0 * fact(a) * fact(b) * fact(c) / fact(a + b + c + 3)
                    };
                    assert!((val - exact).abs() < 1e-12, "dim {dim} deg {deg} ({a},{b},{c}): {val} vs {exact}");
                }
            }
        }
    }

    #[test]
    fn p2_laplace_reproduces_quadratic() {
        // -laplace(u) = -4 with u = x^2 + y^2 on the unit square, Dirichlet from the exact solution
        let m = crate::meshing::geomdrive::plate(1.0, 1.0, 0.2).unwrap().mesh;
        let k = assemble_laplace(&m, 2, &vec![1.0; m.n_cells()]);
        let f = assemble_source(&m, 2, &vec![-4.0; m.n_cells()]);
        let x = m.dof_coordinates(2);
        let exact: Vec<f64> = (0..m.n_dofs(2)).map(|i| x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1]).collect();
        let fixed = boundary_dofs(&m, &m.boundary_labels(), 2);
        let vals: Vec<f64> = fixed.iter().map(|&i| exact[i]).collect();
        let u = solve_dirichlet(&k, &f, &fixed, &vals).unwrap();
        let err = u.iter().zip(&exact).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(err < 1e-12, "{err}");
    }

    #[test]
    fn p2_laplace_3d_reproduces_quadratic() {
        let m = crate::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.34, 3).unwrap().mesh;
        let k = assemble_laplace(&m, 2, &vec![1.0; m.n_cells()]);
        let f = assemble_source(&m, 2, &vec![-6.0; m.n_cells()]);
        let x = m.dof_coordinates(2);
        let exact: Vec<f64> = (0..m.n_dofs(2)).map(|i| x[3 * i].powi(2) + x[3 * i + 1].powi(2) + x[3 * i + 2].powi(2)).collect();
        let fixed = boundary_dofs(&m, &m.boundary_labels(), 2);
        let vals: Vec<f64> = fixed.iter().map(|&i| exact[i]).collect();
        let u = solve_dirichlet(&k, &f, &fixed, &vals).unwrap();
        let err = u.iter().zip(&exact).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(err < 1e-11, "{err}");
    }
}
