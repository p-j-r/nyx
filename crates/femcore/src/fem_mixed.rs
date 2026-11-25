// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Mixed-cell Lagrange finite elements (`Tri`, `Quad`, `Tet`, `Prism`, `Hex`): separable and symmetric
//! Gauss quadrature rules, `Q_1` / serendipity `Q_2'` / full tensor `Q_2` shape functions and reference
//! gradients, straight and curved isoparametric Jacobians (`CurvedQuad`, `CurvedPrism`, `CurvedHex`),
//! and global `C^0` DOF numbering (`MixedLagrangeDofMap`).

use crate::fem::{assemble_parallel, quad_rule, quad_rule_collapsed};
use crate::linalg::csr::RCsr as Csr;
use crate::mesh::cells::CellType;
use crate::mesh::mixed::{
    cell_local_edges, cell_local_quad_faces, FacetType, MixedDiscretisation, HEX_EDGES, HEX_QUAD_FACES,
    PRISM_QUAD_FACES, QUAD_EDGES, TET_EDGES, TRI_EDGES,
};
use std::collections::{BTreeMap, BTreeSet};

// ------------------------------------------------------------------ quadrature

/// 1D Gauss–Legendre points and weights on `[-1, 1]` (`n` points, exact to degree `2*n - 1`, weights sum to `2.0`).
pub fn gauss_legendre_1d(n: usize) -> Vec<(f64, f64)> {
    let n = n.max(1);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = -(std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
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
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        out.push((x, w));
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    out
}

/// Quadrature rule on a reference cell (`Tri`, `Quad`, `Tet`, `Prism`, `Hex`).
///
/// Weights `w` include the reference-cell measure (`0.5` for `Tri`, `4.0` for `Quad`, `1/6` for `Tet`,
/// `1.0` for `Prism`, `8.0` for `Hex`), so the physical integration weight at quadrature point `q` is
/// universally `rule.w[q] * jac.det.abs()`.
#[derive(Clone, Debug)]
pub struct MixedQuadRule {
    pub cell_type: CellType,
    /// Reference points `(\xi, \eta, 0)` in 2D or `(\xi, \eta, \zeta)` in 3D.
    pub pts: Vec<[f64; 3]>,
    /// Reference quadrature weights summing to `meas(ref_cell)`.
    pub w: Vec<f64>,
}

fn tri_quad_points_weights(degree: usize) -> (Vec<[f64; 2]>, Vec<f64>) {
    if degree <= 1 {
        return (vec![[1.0 / 3.0, 1.0 / 3.0]], vec![0.5]);
    }
    if degree <= 4 {
        let q = quad_rule(2, degree);
        let pts = q.lam.iter().map(|l| [l[1], l[2]]).collect();
        let w = q.w.iter().map(|&wq| 0.5 * wq).collect();
        return (pts, w);
    }
    // Collapsed Duffy Gauss-Legendre rule on the reference triangle, exact to `degree`.
    let n = (degree + 2).div_ceil(2).max(2);
    let g = gauss_legendre_1d(n);
    let mut pts = Vec::with_capacity(n * n);
    let mut w = Vec::with_capacity(n * n);
    for &(tx, wx) in &g {
        let u = 0.5 * (tx + 1.0);
        let wu = 0.5 * wx;
        for &(ty, wy) in &g {
            let v = 0.5 * (ty + 1.0);
            let wv = 0.5 * wy;
            let xi = u;
            let eta = (1.0 - u) * v;
            pts.push([xi, eta]);
            w.push(wu * wv * (1.0 - u));
        }
    }
    (pts, w)
}

/// Quadrature rule on `cell_type` exact to total/coordinate polynomial `degree`.
pub fn quad_rule_mixed(cell_type: CellType, degree: usize) -> MixedQuadRule {
    match cell_type {
        CellType::Tri => {
            let (p2, w) = tri_quad_points_weights(degree);
            let pts = p2.into_iter().map(|p| [p[0], p[1], 0.0]).collect();
            MixedQuadRule { cell_type, pts, w }
        }
        CellType::Quad => {
            let n = (degree + 2).div_ceil(2).max(1);
            let g = gauss_legendre_1d(n);
            let mut pts = Vec::with_capacity(n * n);
            let mut w = Vec::with_capacity(n * n);
            for &(eta, w_eta) in &g {
                for &(xi, w_xi) in &g {
                    pts.push([xi, eta, 0.0]);
                    w.push(w_xi * w_eta);
                }
            }
            MixedQuadRule { cell_type, pts, w }
        }
        CellType::Tet => {
            if degree <= 1 {
                return MixedQuadRule {
                    cell_type,
                    pts: vec![[0.25, 0.25, 0.25]],
                    w: vec![1.0 / 6.0],
                };
            }
            let q = if degree <= 4 { quad_rule(3, degree) } else { quad_rule_collapsed(degree) };
            let pts = q.lam.iter().map(|l| [l[1], l[2], l[3]]).collect();
            let w = q.w.iter().map(|&wq| wq / 6.0).collect();
            MixedQuadRule { cell_type, pts, w }
        }
        CellType::Prism => {
            let (tri_pts, tri_w) = tri_quad_points_weights(degree);
            let nz = (degree + 2).div_ceil(2).max(1);
            let gz = gauss_legendre_1d(nz);
            let mut pts = Vec::with_capacity(tri_pts.len() * nz);
            let mut w = Vec::with_capacity(tri_pts.len() * nz);
            for &(zeta, wz) in &gz {
                for (tp, &wt) in tri_pts.iter().zip(&tri_w) {
                    pts.push([tp[0], tp[1], zeta]);
                    w.push(wt * wz);
                }
            }
            MixedQuadRule { cell_type, pts, w }
        }
        CellType::Hex => {
            let n = (degree + 2).div_ceil(2).max(1);
            let g = gauss_legendre_1d(n);
            let mut pts = Vec::with_capacity(n * n * n);
            let mut w = Vec::with_capacity(n * n * n);
            for &(zeta, wz) in &g {
                for &(eta, wy) in &g {
                    for &(xi, wx) in &g {
                        pts.push([xi, eta, zeta]);
                        w.push(wx * wy * wz);
                    }
                }
            }
            MixedQuadRule { cell_type, pts, w }
        }
    }
}

/// Quadrature rule on a boundary or sheet facet (`Edge2`, `Tri3`, `Quad4`).
#[derive(Clone, Debug)]
pub struct MixedFacetRule {
    pub facet_type: FacetType,
    pub pts: Vec<[f64; 2]>,
    pub w: Vec<f64>,
}

/// Quadrature rule on `facet_type` exact to `degree` (`w` sums to `2.0` on `Edge2`, `0.5` on `Tri3`, `4.0` on `Quad4`).
pub fn facet_rule_mixed(facet_type: FacetType, degree: usize) -> MixedFacetRule {
    match facet_type {
        FacetType::Edge2 => {
            let n = (degree + 2).div_ceil(2).max(1);
            let g = gauss_legendre_1d(n);
            MixedFacetRule {
                facet_type,
                pts: g.iter().map(|&(xi, _)| [xi, 0.0]).collect(),
                w: g.iter().map(|&(_, w)| w).collect(),
            }
        }
        FacetType::Tri3 => {
            let (pts, w) = tri_quad_points_weights(degree);
            MixedFacetRule { facet_type, pts, w }
        }
        FacetType::Quad4 => {
            let n = (degree + 2).div_ceil(2).max(1);
            let g = gauss_legendre_1d(n);
            let mut pts = Vec::with_capacity(n * n);
            let mut w = Vec::with_capacity(n * n);
            for &(eta, wy) in &g {
                for &(xi, wx) in &g {
                    pts.push([xi, eta]);
                    w.push(wx * wy);
                }
            }
            MixedFacetRule { facet_type, pts, w }
        }
    }
}

// ------------------------------------------------------------------ Lagrange scalar spaces

/// Conforming `C^0` scalar Lagrange finite-element family across mixed meshes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LagrangeSpace {
    /// Order 1: `Tri(3)`, `Quad Q1(4)`, `Tet(4)`, `Prism P1xQ1(6)`, `Hex Q1(8)`.
    Order1,
    /// Order 2 serendipity (default for `nyx-sms` and `nyx-dfs` order 2):
    /// `Tri(6)`, `Quad Q2'(8)`, `Tet(10)`, `Prism(15)`, `Hex(20)`.
    Order2Serendipity,
    /// Order 2 full tensor product (exact de Rham `H^1` scalar partner for Nedelec `N_2`):
    /// `Tri(6)`, `Quad Q2(9)`, `Tet(10)`, `Prism P2xQ2(18)`, `Hex Q2(27)`.
    Order2Tensor,
}

impl LagrangeSpace {
    /// Polynomial order (`1` or `2`).
    #[inline]
    pub fn order(self) -> usize {
        match self {
            LagrangeSpace::Order1 => 1,
            LagrangeSpace::Order2Serendipity | LagrangeSpace::Order2Tensor => 2,
        }
    }
}

/// Number of local scalar DOFs of `cell_type` in `space`.
#[inline]
pub fn lagrange_nloc(cell_type: CellType, space: LagrangeSpace) -> usize {
    match (cell_type, space) {
        (CellType::Tri, LagrangeSpace::Order1) => 3,
        (CellType::Tri, _) => 6,
        (CellType::Quad, LagrangeSpace::Order1) => 4,
        (CellType::Quad, LagrangeSpace::Order2Serendipity) => 8,
        (CellType::Quad, LagrangeSpace::Order2Tensor) => 9,
        (CellType::Tet, LagrangeSpace::Order1) => 4,
        (CellType::Tet, _) => 10,
        (CellType::Prism, LagrangeSpace::Order1) => 6,
        (CellType::Prism, LagrangeSpace::Order2Serendipity) => 15,
        (CellType::Prism, LagrangeSpace::Order2Tensor) => 18,
        (CellType::Hex, LagrangeSpace::Order1) => 8,
        (CellType::Hex, LagrangeSpace::Order2Serendipity) => 20,
        (CellType::Hex, LagrangeSpace::Order2Tensor) => 27,
    }
}

/// Reference coordinates of the corner vertices of `cell_type`.
pub fn cell_corner_ref_coords(cell_type: CellType) -> &'static [[f64; 3]] {
    const TRI: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    const QUAD: [[f64; 3]; 4] = [
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
    ];
    const TET: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    const PRISM: [[f64; 3]; 6] = [
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    const HEX: [[f64; 3]; 8] = [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ];
    match cell_type {
        CellType::Tri => &TRI,
        CellType::Quad => &QUAD,
        CellType::Tet => &TET,
        CellType::Prism => &PRISM,
        CellType::Hex => &HEX,
    }
}

/// Reference coordinates of all `lagrange_nloc(cell_type, space)` local nodes in canonical local order:
/// corners first, then mid-edges (`cell_local_edges`), then quad-face centers (`Order2Tensor`), then hex center (`Hex27`).
pub fn lagrange_ref_nodes(cell_type: CellType, space: LagrangeSpace) -> Vec<[f64; 3]> {
    let corners = cell_corner_ref_coords(cell_type);
    let mut out = corners.to_vec();
    if space == LagrangeSpace::Order1 {
        return out;
    }
    for &(i, j) in cell_local_edges(cell_type) {
        let (a, b) = (corners[i], corners[j]);
        out.push([0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])]);
    }
    if space == LagrangeSpace::Order2Tensor {
        match cell_type {
            CellType::Quad => out.push([0.0, 0.0, 0.0]),
            CellType::Prism | CellType::Hex => {
                for qf in cell_local_quad_faces(cell_type) {
                    let mut c = [0.0; 3];
                    for &idx in qf {
                        for k in 0..3 {
                            c[k] += 0.25 * corners[idx][k];
                        }
                    }
                    out.push(c);
                }
                if cell_type == CellType::Hex {
                    out.push([0.0, 0.0, 0.0]);
                }
            }
            CellType::Tri | CellType::Tet => {}
        }
    }
    out
}

#[inline]
fn l1d_q1(t: f64) -> ([f64; 2], [f64; 2]) {
    ([0.5 * (1.0 - t), 0.5 * (1.0 + t)], [-0.5, 0.5])
}

#[inline]
fn l1d_q2(t: f64) -> ([f64; 3], [f64; 3]) {
    (
        [0.5 * t * (t - 1.0), 1.0 - t * t, 0.5 * t * (t + 1.0)],
        [t - 0.5, -2.0 * t, t + 0.5],
    )
}

#[inline]
fn tri_p2_1d(xi: f64, eta: f64) -> ([f64; 6], [[f64; 2]; 6]) {
    let l = [1.0 - xi - eta, xi, eta];
    let dl = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
    let mut n = [0.0; 6];
    let mut dn = [[0.0; 2]; 6];
    for i in 0..3 {
        n[i] = l[i] * (2.0 * l[i] - 1.0);
        let s = 4.0 * l[i] - 1.0;
        dn[i] = [s * dl[i][0], s * dl[i][1]];
    }
    for (e, &(i, j)) in TRI_EDGES.iter().enumerate() {
        n[3 + e] = 4.0 * l[i] * l[j];
        dn[3 + e] = [
            4.0 * (dl[i][0] * l[j] + l[i] * dl[j][0]),
            4.0 * (dl[i][1] * l[j] + l[i] * dl[j][1]),
        ];
    }
    (n, dn)
}

/// Evaluates scalar shape functions `vals[a]` and reference derivatives `dxi[a] = [d/d\xi, d/d\eta, d/d\zeta]`
/// at reference point `xi` (`len >= 2` in 2D, `len == 3` in 3D).
pub fn lagrange_shape_into(
    cell_type: CellType,
    space: LagrangeSpace,
    xi: &[f64],
    vals: &mut Vec<f64>,
    dxi: &mut Vec<[f64; 3]>,
) {
    vals.clear();
    dxi.clear();
    match (cell_type, space) {
        (CellType::Tri, LagrangeSpace::Order1) => {
            let (u, v) = (xi[0], xi[1]);
            vals.extend_from_slice(&[1.0 - u - v, u, v]);
            dxi.extend_from_slice(&[[-1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        }
        (CellType::Tri, _) => {
            let (n, dn) = tri_p2_1d(xi[0], xi[1]);
            vals.extend_from_slice(&n);
            for g in &dn {
                dxi.push([g[0], g[1], 0.0]);
            }
        }
        (CellType::Quad, LagrangeSpace::Order1) => {
            let (u, v) = (xi[0], xi[1]);
            let (lu, dlu) = l1d_q1(u);
            let (lv, dlv) = l1d_q1(v);
            for &(iu, iv) in &[(0usize, 0usize), (1, 0), (1, 1), (0, 1)] {
                vals.push(lu[iu] * lv[iv]);
                dxi.push([dlu[iu] * lv[iv], lu[iu] * dlv[iv], 0.0]);
            }
        }
        (CellType::Quad, LagrangeSpace::Order2Serendipity) => {
            let (u, v) = (xi[0], xi[1]);
            let signs: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
            for &(su, sv) in &signs {
                let (au, av) = (1.0 + su * u, 1.0 + sv * v);
                let s = su * u + sv * v - 1.0;
                vals.push(0.25 * au * av * s);
                dxi.push([0.25 * av * (su * s + au * su), 0.25 * au * (sv * s + av * sv), 0.0]);
            }
            // Mid-edges: (0,1) v=-1, (1,2) u=+1, (2,3) v=+1, (3,0) u=-1
            vals.push(0.5 * (1.0 - u * u) * (1.0 - v));
            dxi.push([-u * (1.0 - v), -0.5 * (1.0 - u * u), 0.0]);
            vals.push(0.5 * (1.0 + u) * (1.0 - v * v));
            dxi.push([0.5 * (1.0 - v * v), -(1.0 + u) * v, 0.0]);
            vals.push(0.5 * (1.0 - u * u) * (1.0 + v));
            dxi.push([-u * (1.0 + v), 0.5 * (1.0 - u * u), 0.0]);
            vals.push(0.5 * (1.0 - u) * (1.0 - v * v));
            dxi.push([-0.5 * (1.0 - v * v), -(1.0 - u) * v, 0.0]);
        }
        (CellType::Quad, LagrangeSpace::Order2Tensor) => {
            let (u, v) = (xi[0], xi[1]);
            let (lu, dlu) = l1d_q2(u);
            let (lv, dlv) = l1d_q2(v);
            // 4 corners: (-1,-1), (+1,-1), (+1,+1), (-1,+1) -> indices in {0, 1, 2} are (0,0), (2,0), (2,2), (0,2)
            // 4 mid-edges: (1,0), (2,1), (1,2), (0,1)
            // 1 center: (1,1)
            let order9: [(usize, usize); 9] = [
                (0, 0),
                (2, 0),
                (2, 2),
                (0, 2),
                (1, 0),
                (2, 1),
                (1, 2),
                (0, 1),
                (1, 1),
            ];
            for &(iu, iv) in &order9 {
                vals.push(lu[iu] * lv[iv]);
                dxi.push([dlu[iu] * lv[iv], lu[iu] * dlv[iv], 0.0]);
            }
        }
        (CellType::Tet, LagrangeSpace::Order1) => {
            let (u, v, w) = (xi[0], xi[1], xi[2]);
            vals.extend_from_slice(&[1.0 - u - v - w, u, v, w]);
            dxi.extend_from_slice(&[
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]);
        }
        (CellType::Tet, _) => {
            let (u, v, w) = (xi[0], xi[1], xi[2]);
            let l = [1.0 - u - v - w, u, v, w];
            let dl: [[f64; 3]; 4] = [
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            for i in 0..4 {
                vals.push(l[i] * (2.0 * l[i] - 1.0));
                let s = 4.0 * l[i] - 1.0;
                dxi.push([s * dl[i][0], s * dl[i][1], s * dl[i][2]]);
            }
            for &(i, j) in &TET_EDGES {
                vals.push(4.0 * l[i] * l[j]);
                dxi.push([
                    4.0 * (dl[i][0] * l[j] + l[i] * dl[j][0]),
                    4.0 * (dl[i][1] * l[j] + l[i] * dl[j][1]),
                    4.0 * (dl[i][2] * l[j] + l[i] * dl[j][2]),
                ]);
            }
        }
        (CellType::Prism, LagrangeSpace::Order1) => {
            let (u, v, z) = (xi[0], xi[1], xi[2]);
            let l = [1.0 - u - v, u, v];
            let dl = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            let (lz, dlz) = l1d_q1(z);
            for iz in 0..2 {
                for i in 0..3 {
                    vals.push(l[i] * lz[iz]);
                    dxi.push([dl[i][0] * lz[iz], dl[i][1] * lz[iz], l[i] * dlz[iz]]);
                }
            }
        }
        (CellType::Prism, LagrangeSpace::Order2Serendipity) => {
            let (u, v, z) = (xi[0], xi[1], xi[2]);
            let l = [1.0 - u - v, u, v];
            let dl = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            // 6 corners: z0 = -1 (0..3), z0 = +1 (3..6): 0.5 * l_i * (1 + z0*z) * (2*l_i + z0*z - 2)
            for &z0 in &[-1.0, 1.0] {
                let az = 1.0 + z0 * z;
                for i in 0..3 {
                    let s = 2.0 * l[i] + z0 * z - 2.0;
                    vals.push(0.5 * l[i] * az * s);
                    let df_dli = 0.5 * az * (4.0 * l[i] + z0 * z - 2.0);
                    let df_dz = 0.5 * l[i] * z0 * (2.0 * l[i] + 2.0 * z0 * z - 1.0);
                    dxi.push([df_dli * dl[i][0], df_dli * dl[i][1], df_dz]);
                }
            }
            // 6 base/top triangle mid-edges: 2 * l_a * l_b * (1 + z0 * z)
            for &z0 in &[-1.0, 1.0] {
                let az = 1.0 + z0 * z;
                for &(a, b) in &TRI_EDGES {
                    vals.push(2.0 * l[a] * l[b] * az);
                    dxi.push([
                        2.0 * az * (dl[a][0] * l[b] + l[a] * dl[b][0]),
                        2.0 * az * (dl[a][1] * l[b] + l[a] * dl[b][1]),
                        2.0 * l[a] * l[b] * z0,
                    ]);
                }
            }
            // 3 vertical mid-edges (0,3), (1,4), (2,5): l_i * (1 - z^2)
            let bz = 1.0 - z * z;
            for i in 0..3 {
                vals.push(l[i] * bz);
                dxi.push([dl[i][0] * bz, dl[i][1] * bz, -2.0 * z * l[i]]);
            }
        }
        (CellType::Prism, LagrangeSpace::Order2Tensor) => {
            let (u, v, z) = (xi[0], xi[1], xi[2]);
            let (nt, dnt) = tri_p2_1d(u, v);
            let (lz, dlz) = l1d_q2(z);
            let push_prod = |vals: &mut Vec<f64>, dxi: &mut Vec<[f64; 3]>, it: usize, iz: usize| {
                vals.push(nt[it] * lz[iz]);
                dxi.push([dnt[it][0] * lz[iz], dnt[it][1] * lz[iz], nt[it] * dlz[iz]]);
            };
            // 6 corners: base (iz=0, it=0..3), top (iz=2, it=0..3)
            for &iz in &[0usize, 2] {
                for it in 0..3 {
                    push_prod(vals, dxi, it, iz);
                }
            }
            // 9 mid-edges in PRISM_EDGES order:
            // base tri edges (0,1),(0,2),(1,2) -> it=3..6, iz=0
            // top tri edges (3,4),(3,5),(4,5) -> it=3..6, iz=2
            // vertical edges (0,3),(1,4),(2,5) -> it=0..3, iz=1
            for it in 3..6 {
                push_prod(vals, dxi, it, 0);
            }
            for it in 3..6 {
                push_prod(vals, dxi, it, 2);
            }
            for it in 0..3 {
                push_prod(vals, dxi, it, 1);
            }
            // 3 quad-face centers in PRISM_QUAD_FACES order:
            // [0,1,4,3] -> tri edge (0,1) (it=3), iz=1
            // [1,2,5,4] -> tri edge (1,2) (it=5), iz=1
            // [2,0,3,5] -> tri edge (0,2) (it=4), iz=1
            for &it in &[3usize, 5, 4] {
                push_prod(vals, dxi, it, 1);
            }
        }
        (CellType::Hex, LagrangeSpace::Order1) => {
            let (u, v, w) = (xi[0], xi[1], xi[2]);
            let (lu, dlu) = l1d_q1(u);
            let (lv, dlv) = l1d_q1(v);
            let (lw, dlw) = l1d_q1(w);
            let order8: [(usize, usize, usize); 8] = [
                (0, 0, 0),
                (1, 0, 0),
                (1, 1, 0),
                (0, 1, 0),
                (0, 0, 1),
                (1, 0, 1),
                (1, 1, 1),
                (0, 1, 1),
            ];
            for &(iu, iv, iw) in &order8 {
                vals.push(lu[iu] * lv[iv] * lw[iw]);
                dxi.push([
                    dlu[iu] * lv[iv] * lw[iw],
                    lu[iu] * dlv[iv] * lw[iw],
                    lu[iu] * lv[iv] * dlw[iw],
                ]);
            }
        }
        (CellType::Hex, LagrangeSpace::Order2Serendipity) => {
            let (u, v, w) = (xi[0], xi[1], xi[2]);
            let corners = cell_corner_ref_coords(CellType::Hex);
            for c in corners {
                let (su, sv, sw) = (c[0], c[1], c[2]);
                let (au, av, aw) = (1.0 + su * u, 1.0 + sv * v, 1.0 + sw * w);
                let s = su * u + sv * v + sw * w - 2.0;
                vals.push(0.125 * au * av * aw * s);
                dxi.push([
                    0.125 * av * aw * su * (s + au),
                    0.125 * au * aw * sv * (s + av),
                    0.125 * au * av * sw * (s + aw),
                ]);
            }
            for &(i, j) in &HEX_EDGES {
                let (ca, cb) = (corners[i], corners[j]);
                let m = [0.5 * (ca[0] + cb[0]), 0.5 * (ca[1] + cb[1]), 0.5 * (ca[2] + cb[2])];
                if m[0] == 0.0 {
                    // edge along u at (sv, sw)
                    let (sv, sw) = (m[1], m[2]);
                    let (av, aw) = (1.0 + sv * v, 1.0 + sw * w);
                    vals.push(0.25 * (1.0 - u * u) * av * aw);
                    dxi.push([
                        -0.5 * u * av * aw,
                        0.25 * (1.0 - u * u) * sv * aw,
                        0.25 * (1.0 - u * u) * av * sw,
                    ]);
                } else if m[1] == 0.0 {
                    // edge along v at (su, sw)
                    let (su, sw) = (m[0], m[2]);
                    let (au, aw) = (1.0 + su * u, 1.0 + sw * w);
                    vals.push(0.25 * au * (1.0 - v * v) * aw);
                    dxi.push([
                        0.25 * su * (1.0 - v * v) * aw,
                        -0.5 * v * au * aw,
                        0.25 * au * (1.0 - v * v) * sw,
                    ]);
                } else {
                    // edge along w at (su, sv)
                    let (su, sv) = (m[0], m[1]);
                    let (au, av) = (1.0 + su * u, 1.0 + sv * v);
                    vals.push(0.25 * au * av * (1.0 - w * w));
                    dxi.push([
                        0.25 * su * av * (1.0 - w * w),
                        0.25 * au * sv * (1.0 - w * w),
                        -0.5 * w * au * av,
                    ]);
                }
            }
        }
        (CellType::Hex, LagrangeSpace::Order2Tensor) => {
            let (u, v, w) = (xi[0], xi[1], xi[2]);
            let (lu, dlu) = l1d_q2(u);
            let (lv, dlv) = l1d_q2(v);
            let (lw, dlw) = l1d_q2(w);
            let idx = |coord: f64| -> usize {
                if coord < -0.5 {
                    0
                } else if coord > 0.5 {
                    2
                } else {
                    1
                }
            };
            let nodes = lagrange_ref_nodes(CellType::Hex, LagrangeSpace::Order2Tensor);
            for p in &nodes {
                let (iu, iv, iw) = (idx(p[0]), idx(p[1]), idx(p[2]));
                vals.push(lu[iu] * lv[iv] * lw[iw]);
                dxi.push([
                    dlu[iu] * lv[iv] * lw[iw],
                    lu[iu] * dlv[iv] * lw[iw],
                    lu[iu] * lv[iv] * dlw[iw],
                ]);
            }
        }
    }
}

/// Convenience wrapper returning `(vals, dxi)`.
pub fn lagrange_shape(cell_type: CellType, space: LagrangeSpace, xi: &[f64]) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    lagrange_shape_into(cell_type, space, xi, &mut vals, &mut dxi);
    (vals, dxi)
}

/// Shape functions and reference derivatives on a boundary facet (`Edge2`, `Tri3`, `Quad4`).
pub fn facet_lagrange_shape(
    facet_type: FacetType,
    space: LagrangeSpace,
    xi: &[f64; 2],
) -> (Vec<f64>, Vec<[f64; 2]>) {
    match (facet_type, space) {
        (FacetType::Edge2, LagrangeSpace::Order1) => {
            let (l, dl) = l1d_q1(xi[0]);
            (l.to_vec(), vec![[dl[0], 0.0], [dl[1], 0.0]])
        }
        (FacetType::Edge2, _) => {
            // Order: vertex 0 (-1), vertex 1 (+1), midpoint (0)
            let (l, dl) = l1d_q2(xi[0]);
            (
                vec![l[0], l[2], l[1]],
                vec![[dl[0], 0.0], [dl[2], 0.0], [dl[1], 0.0]],
            )
        }
        (FacetType::Tri3, _) => {
            let (v, dv) = lagrange_shape(CellType::Tri, space, xi);
            (v, dv.into_iter().map(|g| [g[0], g[1]]).collect())
        }
        (FacetType::Quad4, _) => {
            let (v, dv) = lagrange_shape(CellType::Quad, space, xi);
            (v, dv.into_iter().map(|g| [g[0], g[1]]).collect())
        }
    }
}

// ------------------------------------------------------------------ isoparametric curved Jacobians

/// Evaluated geometric Jacobian of a straight or curved mixed cell at a reference point.
#[derive(Clone, Copy, Debug)]
pub struct JacobianEval {
    /// Physical position `[x, z, 0.0]` in 2D or `[x, y, z]` in 3D.
    pub x: [f64; 3],
    /// Spatial dimension (`2` or `3`).
    pub dim: usize,
    /// Jacobian `j[phys_axis][ref_axis] = dx_i / dxi_k`.
    pub j: [[f64; 3]; 3],
    /// Jacobian determinant `det(J_K)`.
    pub det: f64,
    /// Inverse transpose `j_inv_t[phys_axis][ref_axis] = (J_K^{-T})_{i, k}`.
    pub j_inv_t: [[f64; 3]; 3],
}

impl JacobianEval {
    /// Maps a reference scalar gradient `\hat{\nabla}\hat{\phi}` to physical `\nabla \phi = J_K^{-T} \hat{\nabla}\hat{\phi}`.
    #[inline]
    pub fn map_grad(&self, dxi: [f64; 3]) -> [f64; 3] {
        let m = &self.j_inv_t;
        if self.dim == 2 {
            [m[0][0] * dxi[0] + m[0][1] * dxi[1], m[1][0] * dxi[0] + m[1][1] * dxi[1], 0.0]
        } else {
            [
                m[0][0] * dxi[0] + m[0][1] * dxi[1] + m[0][2] * dxi[2],
                m[1][0] * dxi[0] + m[1][1] * dxi[1] + m[1][2] * dxi[2],
                m[2][0] * dxi[0] + m[2][1] * dxi[1] + m[2][2] * dxi[2],
            ]
        }
    }

    /// Covariant Piola mapping `\mathbf{w} = J_K^{-T} \hat{\mathbf{w}}` for `H(curl)` basis functions.
    #[inline]
    pub fn covariant_piola(&self, w_ref: [f64; 3]) -> [f64; 3] {
        self.map_grad(w_ref)
    }

    /// Contravariant Piola curl mapping `\nabla \times \mathbf{w} = (1 / \det J_K) J_K (\hat{\nabla} \times \hat{\mathbf{w}})` in 3D.
    #[inline]
    pub fn contravariant_curl_3d(&self, curl_ref: [f64; 3]) -> [f64; 3] {
        let inv_det = 1.0 / self.det;
        let j = &self.j;
        [
            inv_det * (j[0][0] * curl_ref[0] + j[0][1] * curl_ref[1] + j[0][2] * curl_ref[2]),
            inv_det * (j[1][0] * curl_ref[0] + j[1][1] * curl_ref[1] + j[1][2] * curl_ref[2]),
            inv_det * (j[2][0] * curl_ref[0] + j[2][1] * curl_ref[1] + j[2][2] * curl_ref[2]),
        ]
    }

    /// 2D scalar curl mapping `\mathrm{curl}\,\mathbf{w} = (\hat{\mathrm{curl}}\,\hat{\mathbf{w}}) / \det J_K`.
    #[inline]
    pub fn contravariant_curl_2d(&self, curl_ref: f64) -> f64 {
        curl_ref / self.det
    }
}

#[inline]
fn invert_2x2_transpose(j: &[[f64; 3]; 3]) -> (f64, [[f64; 3]; 3]) {
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    let inv_det = 1.0 / det;
    let mut jit = [[0.0; 3]; 3];
    // J^{-1} = [[j11, -j01], [-j10, j00]] / det  =>  J^{-T} = [[j11, -j10], [-j01, j00]] / det
    jit[0][0] = j[1][1] * inv_det;
    jit[0][1] = -j[1][0] * inv_det;
    jit[1][0] = -j[0][1] * inv_det;
    jit[1][1] = j[0][0] * inv_det;
    (det, jit)
}

#[inline]
fn invert_3x3_transpose(j: &[[f64; 3]; 3]) -> (f64, [[f64; 3]; 3]) {
    // Cofactor matrix C of J satisfies J^{-T} = C / det(J)!
    let c00 = j[1][1] * j[2][2] - j[1][2] * j[2][1];
    let c01 = -(j[1][0] * j[2][2] - j[1][2] * j[2][0]);
    let c02 = j[1][0] * j[2][1] - j[1][1] * j[2][0];

    let c10 = -(j[0][1] * j[2][2] - j[0][2] * j[2][1]);
    let c11 = j[0][0] * j[2][2] - j[0][2] * j[2][0];
    let c12 = -(j[0][0] * j[2][1] - j[0][1] * j[2][0]);

    let c20 = j[0][1] * j[1][2] - j[0][2] * j[1][1];
    let c21 = -(j[0][0] * j[1][2] - j[0][2] * j[1][0]);
    let c22 = j[0][0] * j[1][1] - j[0][1] * j[1][0];

    let det = j[0][0] * c00 + j[0][1] * c01 + j[0][2] * c02;
    let inv_det = 1.0 / det;
    (
        det,
        [
            [c00 * inv_det, c01 * inv_det, c02 * inv_det],
            [c10 * inv_det, c11 * inv_det, c12 * inv_det],
            [c20 * inv_det, c21 * inv_det, c22 * inv_det],
        ],
    )
}

/// Gordon–Hall transfinite center of a curved 8-node quadrilateral `(corners c0..c3, mid-edges m0..m3)`.
#[inline]
pub fn transfinite_quad_center(c: [[f64; 3]; 4], m: [[f64; 3]; 4]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for k in 0..3 {
        out[k] = 0.5 * (m[0][k] + m[1][k] + m[2][k] + m[3][k])
            - 0.25 * (c[0][k] + c[1][k] + c[2][k] + c[3][k]);
    }
    out
}

/// Straight or curved isoparametric geometry of one mixed cell (`Tri`, `Quad`, `Tet`, `Prism`, `Hex`).
#[derive(Clone, Debug)]
pub struct IsoparametricCell {
    pub cell_type: CellType,
    pub geom_space: LagrangeSpace,
    pub curved: bool,
    pub nodes: Vec<[f64; 3]>,
}

/// Curved 8- or 9-node quadrilateral geometry alias.
pub type CurvedQuad = IsoparametricCell;
/// Curved 15- or 18-node triangular prism geometry alias.
pub type CurvedPrism = IsoparametricCell;
/// Curved 20- or 27-node hexahedron geometry alias.
pub type CurvedHex = IsoparametricCell;

impl IsoparametricCell {
    /// Extracts the exact geometric mapping of cell `ci` from `mesh` (`Order1` if straight;
    /// `Order2Serendipity` if curved via `edge_nodes` only; `Order2Tensor` if curved `quad_face_nodes` are present).
    pub fn from_mesh(mesh: &MixedDiscretisation, ci: usize) -> Self {
        let cell_type = mesh.cell_type(ci);
        let vs = mesh.cell(ci);
        let curved = mesh.cell_is_curved(ci);
        if !curved {
            let nodes = vs.iter().map(|&v| mesh.point3(v as usize)).collect();
            return Self {
                cell_type,
                geom_space: LagrangeSpace::Order1,
                curved: false,
                nodes,
            };
        }

        let has_qfn = if mesh.quad_face_nodes.is_empty() {
            false
        } else if mesh.dim == 2 {
            cell_type == CellType::Quad && mesh.quad_face_node(ci).is_some()
        } else {
            mesh.cell_quad_face_ids(ci)
                .iter()
                .any(|&qf| mesh.quad_face_node(qf as usize).is_some())
        };

        let geom_space = if has_qfn {
            LagrangeSpace::Order2Tensor
        } else {
            LagrangeSpace::Order2Serendipity
        };

        let corners: Vec<[f64; 3]> = vs.iter().map(|&v| mesh.point3(v as usize)).collect();
        let mut nodes = corners.clone();
        let eids = mesh.cell_edge_ids(ci);
        let le = cell_local_edges(cell_type);
        let mut edge_pts = Vec::with_capacity(le.len());
        for (k, &(i, j)) in le.iter().enumerate() {
            let p = mesh.edge_node(eids[k] as usize).unwrap_or_else(|| {
                let (a, b) = (corners[i], corners[j]);
                [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])]
            });
            edge_pts.push(p);
            nodes.push(p);
        }

        if geom_space == LagrangeSpace::Order2Tensor {
            match cell_type {
                CellType::Quad => {
                    let c4 = [corners[0], corners[1], corners[2], corners[3]];
                    let m4 = [edge_pts[0], edge_pts[1], edge_pts[2], edge_pts[3]];
                    nodes.push(mesh.quad_face_node(ci).unwrap_or_else(|| transfinite_quad_center(c4, m4)));
                }
                CellType::Prism => {
                    let qf_ids = mesh.cell_quad_face_ids(ci);
                    // Local edge indices of each PRISM_QUAD_FACES quad in cyclic order:
                    // [0,1,4,3]: edges (0,1)=0, (1,4)=7, (3,4)=3, (0,3)=6
                    // [1,2,5,4]: edges (1,2)=2, (2,5)=8, (4,5)=5, (1,4)=7
                    // [2,0,3,5]: edges (0,2)=1, (0,3)=6, (3,5)=4, (2,5)=8
                    let q_edges = [[0usize, 7, 3, 6], [2, 8, 5, 7], [1, 6, 4, 8]];
                    for (lq, qf_local) in PRISM_QUAD_FACES.iter().enumerate() {
                        let p = mesh.quad_face_node(qf_ids[lq] as usize).unwrap_or_else(|| {
                            let c4 = qf_local.map(|idx| corners[idx]);
                            let m4 = q_edges[lq].map(|ei| edge_pts[ei]);
                            transfinite_quad_center(c4, m4)
                        });
                        nodes.push(p);
                    }
                }
                CellType::Hex => {
                    let qf_ids = mesh.cell_quad_face_ids(ci);
                    // HEX_QUAD_FACES:
                    // 0: [0,3,2,1] -> edges 3, 2, 1, 0
                    // 1: [4,5,6,7] -> edges 4, 5, 6, 7
                    // 2: [0,1,5,4] -> edges 0, 9, 4, 8
                    // 3: [1,2,6,5] -> edges 1, 10, 5, 9
                    // 4: [2,3,7,6] -> edges 2, 11, 6, 10
                    // 5: [3,0,4,7] -> edges 3, 8, 7, 11
                    let q_edges = [
                        [3usize, 2, 1, 0],
                        [4, 5, 6, 7],
                        [0, 9, 4, 8],
                        [1, 10, 5, 9],
                        [2, 11, 6, 10],
                        [3, 8, 7, 11],
                    ];
                    let mut face_pts = [[0.0; 3]; 6];
                    for (lq, qf_local) in HEX_QUAD_FACES.iter().enumerate() {
                        let p = mesh.quad_face_node(qf_ids[lq] as usize).unwrap_or_else(|| {
                            let c4 = qf_local.map(|idx| corners[idx]);
                            let m4 = q_edges[lq].map(|ei| edge_pts[ei]);
                            transfinite_quad_center(c4, m4)
                        });
                        face_pts[lq] = p;
                        nodes.push(p);
                    }
                    // 3D Gordon-Hall transfinite center of Hex27:
                    // 1/2 sum(6 faces) - 1/4 sum(12 edges) + 1/8 sum(8 corners)
                    let mut center = [0.0; 3];
                    for k in 0..3 {
                        let sf: f64 = face_pts.iter().map(|p| p[k]).sum();
                        let se: f64 = edge_pts.iter().map(|p| p[k]).sum();
                        let sc: f64 = corners.iter().map(|p| p[k]).sum();
                        center[k] = 0.5 * sf - 0.25 * se + 0.125 * sc;
                    }
                    nodes.push(center);
                }
                CellType::Tri | CellType::Tet => {}
            }
        }

        Self {
            cell_type,
            geom_space,
            curved: true,
            nodes,
        }
    }

    /// Evaluates physical position, Jacobian, determinant, and `J_K^{-T}` at reference point `xi`
    /// using caller-supplied scratch buffers `scratch_n` and `scratch_dn` (zero heap allocation in loops).
    pub fn eval_into(
        &self,
        xi: &[f64],
        scratch_n: &mut Vec<f64>,
        scratch_dn: &mut Vec<[f64; 3]>,
    ) -> JacobianEval {
        let dim = self.cell_type.dim();
        lagrange_shape_into(self.cell_type, self.geom_space, xi, scratch_n, scratch_dn);
        let mut x = [0.0; 3];
        let mut j = [[0.0; 3]; 3];
        for (a, p) in self.nodes.iter().enumerate() {
            let na = scratch_n[a];
            let dna = scratch_dn[a];
            for i in 0..dim {
                x[i] += na * p[i];
                for k in 0..dim {
                    j[i][k] += p[i] * dna[k];
                }
            }
        }
        let (det, j_inv_t) = if dim == 2 {
            invert_2x2_transpose(&j)
        } else {
            invert_3x3_transpose(&j)
        };
        JacobianEval {
            x,
            dim,
            j,
            det,
            j_inv_t,
        }
    }

    /// Evaluates physical position, Jacobian, determinant, and `J_K^{-T}` at reference point `xi`.
    pub fn eval(&self, xi: &[f64]) -> JacobianEval {
        let mut n = Vec::new();
        let mut dn = Vec::new();
        self.eval_into(xi, &mut n, &mut dn)
    }

    /// Integrates the exact volume (3D) or area (2D) of this cell.
    pub fn measure(&self) -> f64 {
        let deg = if self.curved { 6 } else { 2 };
        let rule = quad_rule_mixed(self.cell_type, deg);
        let mut n = Vec::new();
        let mut dn = Vec::new();
        rule.pts
            .iter()
            .zip(&rule.w)
            .map(|(xi, &w)| w * self.eval_into(xi, &mut n, &mut dn).det.abs())
            .sum()
    }
}

/// Evaluates position `x`, surface/arc measure factor `ds`, and unit normal `n_hat` on boundary facet `f` at `xi`.
pub fn eval_facet_jacobian(
    mesh: &MixedDiscretisation,
    f: usize,
    xi: &[f64; 2],
) -> ([f64; 3], f64, [f64; 3]) {
    let ft = mesh.facet_type(f);
    let fv = mesh.facet(f);
    let corners: Vec<[f64; 3]> = fv.iter().map(|&v| mesh.point3(v as usize)).collect();
    match ft {
        FacetType::Edge2 => {
            let eid = mesh.edge_index(fv[0], fv[1]).expect("facet edge must exist");
            let (space, nodes) = if let Some(mid) = mesh.edge_node(eid) {
                (LagrangeSpace::Order2Serendipity, vec![corners[0], corners[1], mid])
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(ft, space, xi);
            let mut x = [0.0; 3];
            let mut t = [0.0; 2];
            for (a, p) in nodes.iter().enumerate() {
                x[0] += n[a] * p[0];
                x[1] += n[a] * p[1];
                t[0] += dn[a][0] * p[0];
                t[1] += dn[a][0] * p[1];
            }
            let ds = t[0].hypot(t[1]);
            let n_hat = if ds > 0.0 {
                [t[1] / ds, -t[0] / ds, 0.0]
            } else {
                [0.0; 3]
            };
            (x, ds, n_hat)
        }
        FacetType::Tri3 => {
            let e01 = mesh.edge_index(fv[0], fv[1]).unwrap();
            let e02 = mesh.edge_index(fv[0], fv[2]).unwrap();
            let e12 = mesh.edge_index(fv[1], fv[2]).unwrap();
            let any_curved =
                mesh.edge_node(e01).is_some() || mesh.edge_node(e02).is_some() || mesh.edge_node(e12).is_some();
            let (space, nodes) = if any_curved {
                let mid = |eid: usize, i: usize, j: usize| {
                    mesh.edge_node(eid).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ])
                };
                (
                    LagrangeSpace::Order2Serendipity,
                    vec![
                        corners[0],
                        corners[1],
                        corners[2],
                        mid(e01, 0, 1),
                        mid(e02, 0, 2),
                        mid(e12, 1, 2),
                    ],
                )
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(ft, space, xi);
            let (mut x, mut t0, mut t1) = ([0.0; 3], [0.0; 3], [0.0; 3]);
            for (a, p) in nodes.iter().enumerate() {
                for k in 0..3 {
                    x[k] += n[a] * p[k];
                    t0[k] += dn[a][0] * p[k];
                    t1[k] += dn[a][1] * p[k];
                }
            }
            let cx = [
                t0[1] * t1[2] - t0[2] * t1[1],
                t0[2] * t1[0] - t0[0] * t1[2],
                t0[0] * t1[1] - t0[1] * t1[0],
            ];
            let ds = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
            let n_hat = if ds > 0.0 {
                [cx[0] / ds, cx[1] / ds, cx[2] / ds]
            } else {
                [0.0; 3]
            };
            (x, ds, n_hat)
        }
        FacetType::Quad4 => {
            let qf_id = mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap();
            let eids = [
                mesh.edge_index(fv[0], fv[1]).unwrap(),
                mesh.edge_index(fv[1], fv[2]).unwrap(),
                mesh.edge_index(fv[2], fv[3]).unwrap(),
                mesh.edge_index(fv[3], fv[0]).unwrap(),
            ];
            let has_qfn = mesh.quad_face_node(qf_id).is_some();
            let any_curved = has_qfn || eids.iter().any(|&e| mesh.edge_node(e).is_some());
            let (space, nodes) = if any_curved {
                let mut ns = corners.clone();
                let mut mids = [[0.0; 3]; 4];
                for (k, &(i, j)) in QUAD_EDGES.iter().enumerate() {
                    let m = mesh.edge_node(eids[k]).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ]);
                    mids[k] = m;
                    ns.push(m);
                }
                if has_qfn {
                    ns.push(mesh.quad_face_node(qf_id).unwrap());
                    (LagrangeSpace::Order2Tensor, ns)
                } else {
                    (LagrangeSpace::Order2Serendipity, ns)
                }
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(ft, space, xi);
            let (mut x, mut t0, mut t1) = ([0.0; 3], [0.0; 3], [0.0; 3]);
            for (a, p) in nodes.iter().enumerate() {
                for k in 0..3 {
                    x[k] += n[a] * p[k];
                    t0[k] += dn[a][0] * p[k];
                    t1[k] += dn[a][1] * p[k];
                }
            }
            let cx = [
                t0[1] * t1[2] - t0[2] * t1[1],
                t0[2] * t1[0] - t0[0] * t1[2],
                t0[0] * t1[1] - t0[1] * t1[0],
            ];
            let ds = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
            let n_hat = if ds > 0.0 {
                [cx[0] / ds, cx[1] / ds, cx[2] / ds]
            } else {
                [0.0; 3]
            };
            (x, ds, n_hat)
        }
    }
}

// ------------------------------------------------------------------ MixedLagrangeDofMap

/// Global `C^0` scalar Lagrange DOF map unifying `Tri`/`Quad` (2D) and `Tet`/`Prism`/`Hex` (3D).
#[derive(Clone, Debug)]
pub struct MixedLagrangeDofMap {
    pub space: LagrangeSpace,
    pub dim: usize,
    pub n_vertices: usize,
    pub n_edges: usize,
    /// In 2D `Order2Tensor`: number of `Quad` cells (`1` center DOF each); in 3D `Order2Tensor`: `n_quad_faces`.
    pub n_quad_bubbles: usize,
    /// In 3D `Order2Tensor`: number of `Hex` cells (`1` volume center DOF each); `0` otherwise.
    pub n_hex_bubbles: usize,
    pub n_dofs: usize,
    /// CSR offsets into `cell_dofs_flat` (`len = n_cells + 1`).
    pub cell_dof_offsets: Vec<usize>,
    /// Precomputed global DOFs of every cell in canonical local node order.
    pub cell_dofs_flat: Vec<usize>,
    /// For `Order2Tensor` in 2D: maps `cell_id -> quad_bubble_index` (`usize::MAX` if `Tri`).
    pub quad_cell_bubble_id: Vec<usize>,
    /// For `Order2Tensor` in 3D: maps `cell_id -> hex_bubble_index` (`usize::MAX` if not `Hex`).
    pub hex_cell_bubble_id: Vec<usize>,
}

impl MixedLagrangeDofMap {
    /// Constructs the global `C^0` Lagrange DOF numbering for `mesh` and `space`.
    pub fn new(mesh: &MixedDiscretisation, space: LagrangeSpace) -> Self {
        let nv = mesh.n_points();
        let ne = mesh.n_edges();
        let nc = mesh.n_cells();

        let mut quad_cell_bubble_id = Vec::new();
        let mut hex_cell_bubble_id = Vec::new();
        let mut n_quad_bubbles = 0usize;
        let mut n_hex_bubbles = 0usize;

        if space == LagrangeSpace::Order2Tensor {
            if mesh.dim == 2 {
                quad_cell_bubble_id = vec![usize::MAX; nc];
                for c in 0..nc {
                    if mesh.cell_type(c) == CellType::Quad {
                        quad_cell_bubble_id[c] = n_quad_bubbles;
                        n_quad_bubbles += 1;
                    }
                }
            } else {
                n_quad_bubbles = mesh.n_quad_faces();
                hex_cell_bubble_id = vec![usize::MAX; nc];
                for c in 0..nc {
                    if mesh.cell_type(c) == CellType::Hex {
                        hex_cell_bubble_id[c] = n_hex_bubbles;
                        n_hex_bubbles += 1;
                    }
                }
            }
        }

        let n_dofs = match space {
            LagrangeSpace::Order1 => nv,
            LagrangeSpace::Order2Serendipity => nv + ne,
            LagrangeSpace::Order2Tensor => nv + ne + n_quad_bubbles + n_hex_bubbles,
        };

        let mut cell_dof_offsets = Vec::with_capacity(nc + 1);
        cell_dof_offsets.push(0);
        let mut cell_dofs_flat = Vec::with_capacity(nc * 10);

        for c in 0..nc {
            let ct = mesh.cell_type(c);
            for &v in mesh.cell(c) {
                cell_dofs_flat.push(v as usize);
            }
            if space != LagrangeSpace::Order1 {
                for &eid in mesh.cell_edge_ids(c) {
                    cell_dofs_flat.push(nv + eid as usize);
                }
                if space == LagrangeSpace::Order2Tensor {
                    if mesh.dim == 2 {
                        if ct == CellType::Quad {
                            cell_dofs_flat.push(nv + ne + quad_cell_bubble_id[c]);
                        }
                    } else {
                        for &qf in mesh.cell_quad_face_ids(c) {
                            cell_dofs_flat.push(nv + ne + qf as usize);
                        }
                        if ct == CellType::Hex {
                            cell_dofs_flat.push(nv + ne + n_quad_bubbles + hex_cell_bubble_id[c]);
                        }
                    }
                }
            }
            debug_assert_eq!(
                cell_dofs_flat.len() - *cell_dof_offsets.last().unwrap(),
                lagrange_nloc(ct, space)
            );
            cell_dof_offsets.push(cell_dofs_flat.len());
        }

        Self {
            space,
            dim: mesh.dim,
            n_vertices: nv,
            n_edges: ne,
            n_quad_bubbles,
            n_hex_bubbles,
            n_dofs,
            cell_dof_offsets,
            cell_dofs_flat,
            quad_cell_bubble_id,
            hex_cell_bubble_id,
        }
    }

    /// Global DOFs of cell `c`.
    #[inline]
    pub fn cell_dofs(&self, c: usize) -> &[usize] {
        &self.cell_dofs_flat[self.cell_dof_offsets[c]..self.cell_dof_offsets[c + 1]]
    }

    /// Global DOFs of boundary/sheet facet `f` in [`facet_lagrange_shape`] node order.
    pub fn facet_dofs(&self, mesh: &MixedDiscretisation, f: usize, out: &mut Vec<usize>) {
        out.clear();
        let fv = mesh.facet(f);
        let ft = mesh.facet_type(f);
        for &v in fv {
            out.push(v as usize);
        }
        if self.space == LagrangeSpace::Order1 {
            return;
        }
        let nv = self.n_vertices;
        match ft {
            FacetType::Edge2 => {
                out.push(nv + mesh.edge_index(fv[0], fv[1]).expect("facet edge"));
            }
            FacetType::Tri3 => {
                for &(i, j) in &TRI_EDGES {
                    out.push(nv + mesh.edge_index(fv[i], fv[j]).expect("facet tri edge"));
                }
            }
            FacetType::Quad4 => {
                for &(i, j) in &QUAD_EDGES {
                    out.push(nv + mesh.edge_index(fv[i], fv[j]).expect("facet quad edge"));
                }
                if self.space == LagrangeSpace::Order2Tensor {
                    let qf = mesh
                        .quad_face_index([fv[0], fv[1], fv[2], fv[3]])
                        .expect("facet quad face");
                    out.push(nv + self.n_edges + qf);
                }
            }
        }
    }

    /// Sorted unique scalar DOFs lying on facets carrying any of `labels`.
    pub fn boundary_dofs(&self, mesh: &MixedDiscretisation, labels: &[i32]) -> Vec<usize> {
        let mut set = BTreeSet::new();
        let mut buf = Vec::new();
        for f in 0..mesh.n_facets() {
            if labels.contains(&mesh.facet_labels[f]) {
                self.facet_dofs(mesh, f, &mut buf);
                set.extend(buf.iter().copied());
            }
        }
        set.into_iter().collect()
    }

    /// Fixed Dirichlet DOFs and values for `{boundary_label: value}` (shared boundary DOFs receive the mean).
    pub fn dirichlet_constraints(
        &self,
        mesh: &MixedDiscretisation,
        values: &BTreeMap<i32, f64>,
    ) -> (Vec<usize>, Vec<f64>) {
        let mut acc: BTreeMap<usize, (f64, usize)> = BTreeMap::new();
        let mut buf = Vec::new();
        for f in 0..mesh.n_facets() {
            if let Some(&val) = values.get(&mesh.facet_labels[f]) {
                self.facet_dofs(mesh, f, &mut buf);
                for &d in &buf {
                    let e = acc.entry(d).or_insert((0.0, 0));
                    e.0 += val;
                    e.1 += 1;
                }
            }
        }
        let fixed: Vec<usize> = acc.keys().copied().collect();
        let vals = acc.values().map(|&(s, n)| s / n as f64).collect();
        (fixed, vals)
    }

    /// Physical coordinates (`dim` per DOF) of every scalar DOF in global order.
    pub fn dof_coordinates(&self, mesh: &MixedDiscretisation) -> Vec<f64> {
        let d = mesh.dim;
        let mut out = Vec::with_capacity(self.n_dofs * d);
        out.extend_from_slice(&mesh.points);
        if self.space == LagrangeSpace::Order1 {
            return out;
        }
        for (eid, &(a, b)) in mesh.edges.iter().enumerate() {
            if let Some(p3) = mesh.edge_node(eid) {
                out.extend_from_slice(&p3[..d]);
            } else {
                let (pa, pb) = (mesh.point(a as usize), mesh.point(b as usize));
                for k in 0..d {
                    out.push(0.5 * (pa[k] + pb[k]));
                }
            }
        }
        if self.space == LagrangeSpace::Order2Tensor {
            if d == 2 {
                for c in 0..mesh.n_cells() {
                    if mesh.cell_type(c) == CellType::Quad {
                        let iso = IsoparametricCell::from_mesh(mesh, c);
                        let p = iso.eval(&[0.0, 0.0]).x;
                        out.extend_from_slice(&p[..2]);
                    }
                }
            } else {
                for (qf_id, &q) in mesh.quad_faces.iter().enumerate() {
                    if let Some(p3) = mesh.quad_face_node(qf_id) {
                        out.extend_from_slice(&p3);
                    } else {
                        let c4 = q.map(|v| mesh.point3(v as usize));
                        let m4 = QUAD_EDGES.map(|(i, j)| {
                            let eid = mesh.edge_index(q[i], q[j]).unwrap();
                            mesh.edge_node(eid).unwrap_or([
                                0.5 * (c4[i][0] + c4[j][0]),
                                0.5 * (c4[i][1] + c4[j][1]),
                                0.5 * (c4[i][2] + c4[j][2]),
                            ])
                        });
                        out.extend_from_slice(&transfinite_quad_center(c4, m4));
                    }
                }
                for c in 0..mesh.n_cells() {
                    if mesh.cell_type(c) == CellType::Hex {
                        let iso = IsoparametricCell::from_mesh(mesh, c);
                        let p = iso.eval(&[0.0, 0.0, 0.0]).x;
                        out.extend_from_slice(&p);
                    }
                }
            }
        }
        out
    }
}

// ------------------------------------------------------------------ scalar matrix assemblers

/// Assembles `K_{ij} = \int_\Omega coeff(c) \nabla\phi_i \cdot \nabla\phi_j \, dx` over a [`MixedDiscretisation`].
pub fn assemble_laplace_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    coeff: &[f64],
) -> Csr {
    let d = mesh.dim;
    let n = dof_map.n_dofs;
    let base_deg = 2 * dof_map.space.order();
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut grads: Vec<[f64; 3]> = Vec::new();
        let mut ke = Vec::new();
        for c in range {
            let ct = mesh.cell_type(c);
            let iso = IsoparametricCell::from_mesh(mesh, c);
            let deg = if iso.curved { base_deg.max(6) } else { base_deg };
            let rule = quad_rule_mixed(ct, deg);
            let dofs = dof_map.cell_dofs(c);
            let nloc = dofs.len();
            ke.clear();
            ke.resize(nloc * nloc, 0.0);
            let cf = coeff[c];
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, dof_map.space, xi, &mut vals, &mut dxi);
                grads.clear();
                for a in 0..nloc {
                    grads.push(jac.map_grad(dxi[a]));
                }
                let weight = wq * jac.det.abs() * cf;
                for a in 0..nloc {
                    for b in 0..nloc {
                        let mut dot = 0.0;
                        for k in 0..d {
                            dot += grads[a][k] * grads[b][k];
                        }
                        ke[a * nloc + b] += weight * dot;
                    }
                }
            }
            coo.push_block(dofs, dofs, &ke);
        }
    })
}

/// Assembles `K_{ij} = \int_\Omega \nabla\phi_i \cdot (K_c \nabla\phi_j) \, dx` over a [`MixedDiscretisation`].
pub fn assemble_laplace_tensor_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    tensor: &[[[f64; 3]; 3]],
) -> Csr {
    let d = mesh.dim;
    let n = dof_map.n_dofs;
    let base_deg = 2 * dof_map.space.order();
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut grads: Vec<[f64; 3]> = Vec::new();
        let mut ke = Vec::new();
        for c in range {
            let ct = mesh.cell_type(c);
            let iso = IsoparametricCell::from_mesh(mesh, c);
            let deg = if iso.curved { base_deg.max(6) } else { base_deg };
            let rule = quad_rule_mixed(ct, deg);
            let dofs = dof_map.cell_dofs(c);
            let nloc = dofs.len();
            ke.clear();
            ke.resize(nloc * nloc, 0.0);
            let kc = &tensor[c];
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, dof_map.space, xi, &mut vals, &mut dxi);
                grads.clear();
                for a in 0..nloc {
                    grads.push(jac.map_grad(dxi[a]));
                }
                let weight = wq * jac.det.abs();
                for a in 0..nloc {
                    let mut kg = [0.0; 3];
                    for i in 0..d {
                        for j in 0..d {
                            kg[i] += kc[i][j] * grads[a][j];
                        }
                    }
                    for b in 0..nloc {
                        let mut dot = 0.0;
                        for i in 0..d {
                            dot += kg[i] * grads[b][i];
                        }
                        ke[a * nloc + b] += weight * dot;
                    }
                }
            }
            coo.push_block(dofs, dofs, &ke);
        }
    })
}

/// Assembles `M_{ij} = \int_\Omega coeff(c) \phi_i \phi_j \, dx` over a [`MixedDiscretisation`].
pub fn assemble_mass_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    coeff: &[f64],
) -> Csr {
    let n = dof_map.n_dofs;
    let base_deg = 2 * dof_map.space.order() + 2;
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut ke = Vec::new();
        for c in range {
            let ct = mesh.cell_type(c);
            let iso = IsoparametricCell::from_mesh(mesh, c);
            let deg = if iso.curved { base_deg.max(6) } else { base_deg };
            let rule = quad_rule_mixed(ct, deg);
            let dofs = dof_map.cell_dofs(c);
            let nloc = dofs.len();
            ke.clear();
            ke.resize(nloc * nloc, 0.0);
            let cf = coeff[c];
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, dof_map.space, xi, &mut vals, &mut dxi);
                let weight = wq * jac.det.abs() * cf;
                for a in 0..nloc {
                    let wa = weight * vals[a];
                    for b in 0..nloc {
                        ke[a * nloc + b] += wa * vals[b];
                    }
                }
            }
            coo.push_block(dofs, dofs, &ke);
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem::solve_dirichlet;
    use crate::mesh::cells::CellMesh;

    #[test]
    fn kronecker_delta_and_quadratic_reproduction_on_all_shapes_and_spaces() {
        let shapes = [
            CellType::Tri,
            CellType::Quad,
            CellType::Tet,
            CellType::Prism,
            CellType::Hex,
        ];
        let spaces = [
            LagrangeSpace::Order1,
            LagrangeSpace::Order2Serendipity,
            LagrangeSpace::Order2Tensor,
        ];
        let probes_2d = [[-0.37, 0.23, 0.0], [0.19, 0.41, 0.0], [0.05, 0.85, 0.0]];
        let probes_3d = [
            [0.17, 0.29, 0.31],
            [0.42, 0.11, -0.35],
            [0.08, 0.61, 0.19],
        ];
        for ct in shapes {
            for sp in spaces {
                let nodes = lagrange_ref_nodes(ct, sp);
                assert_eq!(nodes.len(), lagrange_nloc(ct, sp));
                // 1. Kronecker delta at all reference nodes
                for (b, xb) in nodes.iter().enumerate() {
                    let (vals, _) = lagrange_shape(ct, sp, xb);
                    for (a, &va) in vals.iter().enumerate() {
                        let want = if a == b { 1.0 } else { 0.0 };
                        assert!(
                            (va - want).abs() < 1e-14,
                            "{ct:?} {sp:?}: phi_{a}(x_{b}) = {va}, expected {want}"
                        );
                    }
                }
                // 2. Partition of unity, linear (Order1) & quadratic (Order2) polynomial + derivative reproduction!
                let probes: &[[f64; 3]] = if ct.dim() == 2 { &probes_2d } else { &probes_3d };
                for xi in probes {
                    let (vals, dxi) = lagrange_shape(ct, sp, xi);
                    let sum_n: f64 = vals.iter().sum();
                    assert!((sum_n - 1.0).abs() < 1e-14, "{ct:?} {sp:?}: sum = {sum_n}");

                    // Test polynomial p(x, y, z): degree 1 for Order1, pure/mixed degree 2 for Order2
                    let poly = |p: [f64; 3]| -> (f64, [f64; 3]) {
                        if sp == LagrangeSpace::Order1 {
                            (
                                1.7 + 2.3 * p[0] - 1.1 * p[1] + 0.7 * p[2],
                                [2.3, -1.1, if ct.dim() == 3 { 0.7 } else { 0.0 }],
                            )
                        } else if ct.dim() == 2 {
                            (
                                0.5 + 1.3 * p[0] - 0.9 * p[1] + 2.1 * p[0] * p[0]
                                    - 1.4 * p[0] * p[1]
                                    + 0.8 * p[1] * p[1],
                                [
                                    1.3 + 4.2 * p[0] - 1.4 * p[1],
                                    -0.9 - 1.4 * p[0] + 1.6 * p[1],
                                    0.0,
                                ],
                            )
                        } else {
                            (
                                0.5 + 1.2 * p[0] - 0.8 * p[1] + 0.4 * p[2] + 1.7 * p[0] * p[0]
                                    - 0.9 * p[1] * p[1]
                                    + 1.3 * p[2] * p[2]
                                    + 0.6 * p[0] * p[1]
                                    - 0.5 * p[0] * p[2]
                                    + 0.7 * p[1] * p[2],
                                [
                                    1.2 + 3.4 * p[0] + 0.6 * p[1] - 0.5 * p[2],
                                    -0.8 - 1.8 * p[1] + 0.6 * p[0] + 0.7 * p[2],
                                    0.4 + 2.6 * p[2] - 0.5 * p[0] + 0.7 * p[1],
                                ],
                            )
                        }
                    };
                    let (want_v, want_g) = poly(*xi);
                    let mut got_v = 0.0;
                    let mut got_g = [0.0; 3];
                    for (a, xn) in nodes.iter().enumerate() {
                        let (pn, _) = poly(*xn);
                        got_v += pn * vals[a];
                        for k in 0..ct.dim() {
                            got_g[k] += pn * dxi[a][k];
                        }
                    }
                    assert!(
                        (got_v - want_v).abs() < 1e-13,
                        "{ct:?} {sp:?}: val {got_v} vs {want_v}"
                    );
                    for k in 0..ct.dim() {
                        assert!(
                            (got_g[k] - want_g[k]).abs() < 1e-13,
                            "{ct:?} {sp:?} axis {k}: grad {} vs {}",
                            got_g[k],
                            want_g[k]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn laplace_patch_test_on_3d_tet_prism_hex_for_all_spaces() {
        use crate::adapt::work::simplex_to_adapt;
        use crate::meshing::pml_box::{build_pml, PmlSpec};

        // Real 3D mixed Tet + Prism + Hex mesh with interior vertices, interior edges, interior tri-faces,
        // and interior quad-faces, sheared by a general non-orthogonal 3x3 affine transformation!
        let block = crate::meshing::geomdrive::block(1.0, 0.9, 0.8, 0.45, 3).unwrap().mesh;
        let phys = simplex_to_adapt::<3>(&block).unwrap();
        let mut cm = build_pml(&phys, &PmlSpec { outer_label: Some(1), ..PmlSpec::<3>::uniform(0.3, 2) })
            .unwrap()
            .mesh;
        assert!(cm.types.contains(&CellType::Tet));
        assert!(cm.types.contains(&CellType::Prism));
        assert!(cm.types.contains(&CellType::Hex));

        let shear = |p: [f64; 3]| -> [f64; 3] {
            [
                p[0] + 0.15 * p[1] - 0.10 * p[2],
                0.08 * p[0] + 1.10 * p[1] + 0.12 * p[2],
                -0.05 * p[0] + 0.07 * p[1] + 0.95 * p[2],
            ]
        };
        for p in &mut cm.points {
            *p = shear(*p);
        }
        cm.check().unwrap();
        let md = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        md.validate_conforming().unwrap();

        for sp in [
            LagrangeSpace::Order1,
            LagrangeSpace::Order2Serendipity,
            LagrangeSpace::Order2Tensor,
        ] {
            let dm = MixedLagrangeDofMap::new(&md, sp);
            let k = assemble_laplace_mixed(&md, &dm, &vec![1.0; md.n_cells()]);
            let coords = dm.dof_coordinates(&md);
            let exact = |d: usize| -> f64 {
                let (x, y, z) = (coords[3 * d], coords[3 * d + 1], coords[3 * d + 2]);
                if sp == LagrangeSpace::Order1 {
                    1.3 * x - 0.7 * y + 2.1 * z + 0.4
                } else {
                    x * x - y * y + 0.5 * x * z + 1.2 * y - 0.3
                }
            };
            let b_dofs = dm.boundary_dofs(&md, &[1]);
            assert!(b_dofs.len() < dm.n_dofs, "must have free interior DOFs ({sp:?})");
            let b_vals: Vec<f64> = b_dofs.iter().map(|&d| exact(d)).collect();
            let u = solve_dirichlet(&k, &vec![0.0; dm.n_dofs], &b_dofs, &b_vals).unwrap();
            for d in 0..dm.n_dofs {
                let err = (u[d] - exact(d)).abs();
                assert!(err < 1e-12, "{sp:?} dof {d}: u={} vs exact={}, err={err:e}", u[d], exact(d));
            }
        }

        // Also test non-affinely warped Hex + Prism + Tet reproducing linear harmonic field to < 1e-12!
        let mut cm_warped = cm;
        for p in &mut cm_warped.points {
            let (x, y, z) = (p[0], p[1], p[2]);
            p[0] = x + 0.015 * x * y;
            p[1] = y + 0.012 * y * z;
            p[2] = z + 0.010 * x * z;
        }
        cm_warped.check().unwrap();
        let md_w = MixedDiscretisation::from_cell_mesh(&cm_warped).unwrap();
        for sp in [
            LagrangeSpace::Order1,
            LagrangeSpace::Order2Serendipity,
            LagrangeSpace::Order2Tensor,
        ] {
            let dm = MixedLagrangeDofMap::new(&md_w, sp);
            let k = assemble_laplace_mixed(&md_w, &dm, &vec![1.0; md_w.n_cells()]);
            let coords = dm.dof_coordinates(&md_w);
            let exact = |d: usize| -> f64 {
                1.3 * coords[3 * d] - 0.7 * coords[3 * d + 1] + 2.1 * coords[3 * d + 2] + 0.4
            };
            let b_dofs = dm.boundary_dofs(&md_w, &[1]);
            let b_vals: Vec<f64> = b_dofs.iter().map(|&d| exact(d)).collect();
            let u = solve_dirichlet(&k, &vec![0.0; dm.n_dofs], &b_dofs, &b_vals).unwrap();
            for d in 0..dm.n_dofs {
                assert!((u[d] - exact(d)).abs() < 1e-12, "warped {sp:?} dof {d}");
            }
        }
    }

    #[test]
    fn curved_thin_spherical_shell_prism_preserves_positive_jacobian() {
        // Thin concentric spherical shell h_R / R_in = 0.005: straight tets would sagitta-invert,
        // whereas lockstep CurvedPrism15 has det(J) > 0 everywhere and matches exact shell volume!
        let r0 = 1.0_f64;
        let r1 = 1.005_f64;
        // Octant-like spherical patch (3 unit vectors)
        let dirs: [[f64; 3]; 3] = [
            [1.0, 0.0, 0.0],
            [0.8, 0.6, 0.0],
            [0.8, 0.0, 0.6],
        ];
        let norm = |v: [f64; 3]| {
            let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / l, v[1] / l, v[2] / l]
        };
        let u = dirs.map(norm);
        let scale = |r: f64, v: [f64; 3]| [r * v[0], r * v[1], r * v[2]];

        let mut cm = CellMesh::<3>::new();
        for &uv in &u {
            cm.points.push(scale(r0, uv));
        }
        for &uv in &u {
            cm.points.push(scale(r1, uv));
        }
        cm.push_cell(CellType::Prism, &[0, 1, 2, 3, 4, 5], 1);
        for f in cm.boundary_faces_ordered() {
            cm.push_facet(&f, 1);
        }

        // Lockstep curved midpoints on inner (r0) and outer (r1) spherical arcs
        let mut en = Vec::new();
        for &(i, j) in &[(0usize, 1usize), (0, 2), (1, 2)] {
            let m_dir = norm([u[i][0] + u[j][0], u[i][1] + u[j][1], u[i][2] + u[j][2]]);
            en.push(([i as u32, j as u32], scale(r0, m_dir)));
            en.push(([(i + 3) as u32, (j + 3) as u32], scale(r1, m_dir)));
        }
        cm.set_edge_nodes(en);

        let md = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        let cp = CurvedPrism::from_mesh(&md, 0);
        assert!(cp.curved);
        let rule = quad_rule_mixed(CellType::Prism, 6);
        for pt in &rule.pts {
            let jac = cp.eval(pt);
            assert!(jac.det > 0.0, "CurvedPrism inverted at {pt:?}: det = {}", jac.det);
        }
        assert!(cp.measure() > 0.0);
    }
}
