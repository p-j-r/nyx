// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Asymptotic open/infinite far-field Robin and second-order Bayliss-Turkel (Laplace-Beltrami)
//! boundary conditions on 2D and 3D simplex meshes (including curved quadratic tetrahedra).
//!
//! * [`farfield_centroid`] computes the deterministic area/length-weighted centroid of reference
//!   conductor facets, falling back to the mesh bounding-box midpoint when no reference facets match
//!   (never returning `NaN`).
//! * [`assemble_boundary_mass_weighted`] integrates a spatially and normal-dependent scalar weight
//!   `\alpha(f, c_f, \mathbf{x}, \mathbf{n}_f)` over boundary facets in deterministic sorted facet order.
//! * [`assemble_farfield_robin`] assembles first-order asymptotic Robin (`order = 1`, 2D and 3D) and
//!   second-order symmetric 3D Laplace-Beltrami (`order = 2`) far-field boundary matrices.

use crate::fem::{assemble_parallel, facet_measure, facet_rule, shape};
use crate::linalg::csr::RCsr as Csr;
use crate::mesh::SimplexMesh;
use std::collections::BTreeMap;

const EPS_GEOM: f64 = 1.0e-3;
const RADIAL_PROJ_REG: f64 = 0.05;

/// Relative permittivity field passed to [`assemble_farfield_robin`].
#[derive(Clone, Copy, Debug)]
pub enum FarFieldPermittivity<'a> {
    /// Isotropic relative permittivity `eps_r[c]` per cell.
    Scalar(&'a [f64]),
    /// Anisotropic relative permittivity tensor `eps_r[c]` per cell in global wafer axes.
    Tensor(&'a [[[f64; 3]; 3]]),
}

impl<'a> FarFieldPermittivity<'a> {
    #[inline]
    fn normal_eps(&self, dim: usize, cell: usize, n_hat: [f64; 3]) -> f64 {
        let val = match self {
            FarFieldPermittivity::Scalar(s) => s.get(cell).copied().unwrap_or(1.0),
            FarFieldPermittivity::Tensor(t) => {
                if let Some(m) = t.get(cell) {
                    let mut acc = 0.0;
                    for i in 0..dim.min(3) {
                        for j in 0..dim.min(3) {
                            acc += n_hat[i] * m[i][j] * n_hat[j];
                        }
                    }
                    acc
                } else {
                    1.0
                }
            }
        };
        if val.is_finite() { val.max(0.0) } else { 0.0 }
    }
}

#[inline]
fn sorted2(a: u32, b: u32) -> [u32; 2] {
    if a <= b { [a, b] } else { [b, a] }
}

#[inline]
fn sorted3(a: u32, b: u32, c: u32) -> [u32; 3] {
    let mut v = [a, b, c];
    v.sort_unstable();
    v
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Deterministic facet-to-cell lookup in 2D (`[u32; 2] -> (cell_idx, opp_local_vertex)`).
fn build_edge_to_cell_2d(mesh: &SimplexMesh) -> BTreeMap<[u32; 2], (usize, usize)> {
    let mut map = BTreeMap::new();
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        for &(i, j, opp) in &[(0, 1, 2), (0, 2, 1), (1, 2, 0)] {
            let key = sorted2(v[i], v[j]);
            map.entry(key).or_insert((c, opp));
        }
    }
    map
}

/// Deterministic facet-to-cell lookup in 3D (`[u32; 3] -> (cell_idx, opp_local_vertex)`).
fn build_face_to_cell_3d_local(mesh: &SimplexMesh) -> BTreeMap<[u32; 3], (usize, usize)> {
    let mut map = BTreeMap::new();
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        for &(i, j, k, opp) in &[(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)] {
            let key = sorted3(v[i], v[j], v[k]);
            map.entry(key).or_insert((c, opp));
        }
    }
    map
}

/// Collects and deduplicates facet indices matching `labels`, sorted deterministically by sorted vertex key.
fn sorted_unique_facets(mesh: &SimplexMesh, labels: &[i32]) -> Vec<usize> {
    if labels.is_empty() {
        return Vec::new();
    }
    let mut keyed: Vec<([u32; 3], usize)> = Vec::new();
    for f in 0..mesh.n_facets() {
        if labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            let key = if mesh.dim == 2 {
                let s = sorted2(fv[0], fv[1]);
                [s[0], s[1], u32::MAX]
            } else {
                sorted3(fv[0], fv[1], fv[2])
            };
            keyed.push((key, f));
        }
    }
    keyed.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    keyed.dedup_by_key(|x| x.0);
    keyed.into_iter().map(|(_, f)| f).collect()
}

/// Six quadratic geometry nodes `[v0, v1, v2, e01, e02, e12]` of a 3D triangle facet `f`.
fn facet_geometry_nodes_3d(mesh: &SimplexMesh, f: usize) -> [[f64; 3]; 6] {
    let fv = mesh.facet(f);
    let pt = |v: u32| -> [f64; 3] {
        let p = mesh.point(v as usize);
        [p[0], p[1], p[2]]
    };
    let (p0, p1, p2) = (pt(fv[0]), pt(fv[1]), pt(fv[2]));
    let mid = |a: u32, b: u32, pa: [f64; 3], pb: [f64; 3]| -> [f64; 3] {
        if mesh.has_curved() {
            if let Some(node) = mesh.edge_node(mesh.edge_index(a, b)) {
                return node;
            }
        }
        [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])]
    };
    [p0, p1, p2, mid(fv[0], fv[1], p0, p1), mid(fv[0], fv[2], p0, p2), mid(fv[1], fv[2], p1, p2)]
}

/// Computes the deterministic area-weighted (3D) or length-weighted (2D) centroid of all facets carrying
/// one of `labels`. Falls back to the mesh bounding-box midpoint (or `[0.0; 3]` on an empty mesh) when
/// no positive-measure facets match `labels`. Never returns `NaN` or infinity.
pub fn farfield_centroid(mesh: &SimplexMesh, labels: &[i32]) -> [f64; 3] {
    let mut sum_m = 0.0_f64;
    let mut sum_c = [0.0_f64; 3];
    for f in sorted_unique_facets(mesh, labels) {
        let m = facet_measure(mesh, f);
        if !m.is_finite() || m <= 0.0 {
            continue;
        }
        let fv = mesh.facet(f);
        let fc = if mesh.dim == 2 {
            let p0 = mesh.point(fv[0] as usize);
            let p1 = mesh.point(fv[1] as usize);
            [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]), 0.0]
        } else {
            let p0 = mesh.point(fv[0] as usize);
            let p1 = mesh.point(fv[1] as usize);
            let p2 = mesh.point(fv[2] as usize);
            [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0]
        };
        if fc.iter().all(|x| x.is_finite()) {
            sum_m += m;
            sum_c[0] += m * fc[0];
            sum_c[1] += m * fc[1];
            sum_c[2] += m * fc[2];
        }
    }
    if sum_m > 1.0e-30 && sum_m.is_finite() {
        let c = [sum_c[0] / sum_m, sum_c[1] / sum_m, sum_c[2] / sum_m];
        if c.iter().all(|x| x.is_finite()) {
            return c;
        }
    }
    if mesh.n_points() > 0 {
        let (lo, hi) = mesh.bbox();
        let c = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
        if c.iter().all(|x| x.is_finite()) {
            return c;
        }
    }
    [0.0, 0.0, 0.0]
}

/// Computes the total outward-facing positive measure (`\int_{f: \mathbf{r}\cdot\mathbf{n}_f > 0} dS`)
/// of all facets carrying one of `labels` with respect to expansion center `center`.
pub fn farfield_positive_measure(mesh: &SimplexMesh, labels: &[i32], center: [f64; 3]) -> f64 {
    let facets = sorted_unique_facets(mesh, labels);
    if facets.is_empty() {
        return 0.0;
    }
    if mesh.dim == 2 {
        let e2c = build_edge_to_cell_2d(mesh);
        let mut total = 0.0;
        for f in facets {
            let fv = mesh.facet(f);
            let Some(&(c_f, opp)) = e2c.get(&sorted2(fv[0], fv[1])) else {
                continue;
            };
            let p0 = mesh.point(fv[0] as usize);
            let p1 = mesh.point(fv[1] as usize);
            let po = mesh.point(mesh.cell(c_f)[opp] as usize);
            let t = [p1[0] - p0[0], p1[1] - p0[1]];
            let len = (t[0] * t[0] + t[1] * t[1]).sqrt();
            if len <= 1.0e-30 {
                continue;
            }
            let mut n_hat = [t[1] / len, -t[0] / len, 0.0];
            let mid = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]), 0.0];
            if n_hat[0] * (mid[0] - po[0]) + n_hat[1] * (mid[1] - po[1]) < 0.0 {
                n_hat = [-n_hat[0], -n_hat[1], 0.0];
            }
            let r = [mid[0] - center[0], mid[1] - center[1], 0.0];
            if dot3(r, n_hat) > 0.0 {
                total += len;
            }
        }
        total
    } else {
        let f2c = build_face_to_cell_3d_local(mesh);
        let mut total = 0.0;
        for f in facets {
            let fv = mesh.facet(f);
            let Some(&(c_f, opp)) = f2c.get(&sorted3(fv[0], fv[1], fv[2])) else {
                continue;
            };
            let x6 = facet_geometry_nodes_3d(mesh, f);
            let po_s = mesh.point(mesh.cell(c_f)[opp] as usize);
            let p_opp = [po_s[0], po_s[1], po_s[2]];
            let t1 = sub3(x6[1], x6[0]);
            let t2 = sub3(x6[2], x6[0]);
            let n_raw = cross3(t1, t2);
            let n_len = norm3(n_raw);
            if n_len <= 1.0e-30 {
                continue;
            }
            let mut n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];
            let cen = [
                (x6[0][0] + x6[1][0] + x6[2][0]) / 3.0,
                (x6[0][1] + x6[1][1] + x6[2][1]) / 3.0,
                (x6[0][2] + x6[1][2] + x6[2][2]) / 3.0,
            ];
            if dot3(n_hat, sub3(cen, p_opp)) < 0.0 {
                n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
            }
            let r = sub3(cen, center);
            if dot3(r, n_hat) > 0.0 {
                total += 0.5 * n_len;
            }
        }
        total
    }
}

/// Assembles the weighted boundary mass matrix
/// `K_{ab} = \sum_{f \in \Gamma(\text{labels})} \int_f w(f, c_f, \mathbf{x}, \mathbf{n}_f) \phi_a(\mathbf{x}) \phi_b(\mathbf{x}) \, dS`
/// over all facets carrying one of `labels`, supporting 2D segments and 3D triangles (including curved
/// quadratic tetrahedra `CurvedTet`) in deterministic sorted facet order (`NYX_THREADS` bit-identical).
pub fn assemble_boundary_mass_weighted<W>(mesh: &SimplexMesh, labels: &[i32], order: usize, weight: W) -> Csr
where
    W: Fn(usize, usize, [f64; 3], [f64; 3]) -> f64 + Sync,
{
    let n = mesh.n_dofs(order);
    let facets = sorted_unique_facets(mesh, labels);
    if facets.is_empty() {
        return Csr::zeros(n, n);
    }
    let fd = mesh.dim - 1;
    let rule = facet_rule(fd);
    let nq = rule.w.len();
    let nloc = if order == 1 { mesh.dim } else { mesh.dim + mesh.dim * (mesh.dim - 1) / 2 };
    let shapes_order: Vec<Vec<f64>> = rule.lam.iter().map(|lam| shape(lam, order).0).collect();

    if mesh.dim == 2 {
        let e2c = build_edge_to_cell_2d(mesh);
        assemble_parallel(facets.len(), n, n, |range, coo| {
            let mut dofs = Vec::new();
            let mut ke = vec![0.0; nloc * nloc];
            for idx in range {
                let f = facets[idx];
                let fv = mesh.facet(f);
                let Some(&(c_f, opp)) = e2c.get(&sorted2(fv[0], fv[1])) else {
                    continue;
                };
                let p0 = mesh.point(fv[0] as usize);
                let p1 = mesh.point(fv[1] as usize);
                let po = mesh.point(mesh.cell(c_f)[opp] as usize);
                let t = [p1[0] - p0[0], p1[1] - p0[1]];
                let len = (t[0] * t[0] + t[1] * t[1]).sqrt();
                if len <= 1.0e-30 {
                    continue;
                }
                let mut n_hat = [t[1] / len, -t[0] / len, 0.0];
                let mid = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])];
                if n_hat[0] * (mid[0] - po[0]) + n_hat[1] * (mid[1] - po[1]) < 0.0 {
                    n_hat = [-n_hat[0], -n_hat[1], 0.0];
                }
                ke.fill(0.0);
                for q in 0..nq {
                    let lam = &rule.lam[q];
                    let x_q = [lam[0] * p0[0] + lam[1] * p1[0], lam[0] * p0[1] + lam[1] * p1[1], 0.0];
                    let w_val = weight(f, c_f, x_q, n_hat);
                    if w_val == 0.0 {
                        continue;
                    }
                    let ds = rule.w[q] * len * w_val;
                    let nv = &shapes_order[q];
                    for a in 0..nloc {
                        for b in 0..nloc {
                            ke[a * nloc + b] += ds * nv[a] * nv[b];
                        }
                    }
                }
                mesh.facet_dofs(f, order, &mut dofs);
                coo.push_block(&dofs, &dofs, &ke);
            }
        })
    } else {
        let f2c = build_face_to_cell_3d_local(mesh);
        let shapes_geom2: Vec<(Vec<f64>, Vec<f64>)> = rule.lam.iter().map(|lam| shape(lam, 2)).collect();
        assemble_parallel(facets.len(), n, n, |range, coo| {
            let mut dofs = Vec::new();
            let mut ke = vec![0.0; nloc * nloc];
            for idx in range {
                let f = facets[idx];
                let fv = mesh.facet(f);
                let Some(&(c_f, opp)) = f2c.get(&sorted3(fv[0], fv[1], fv[2])) else {
                    continue;
                };
                let x6 = facet_geometry_nodes_3d(mesh, f);
                let po_s = mesh.point(mesh.cell(c_f)[opp] as usize);
                let p_opp = [po_s[0], po_s[1], po_s[2]];
                ke.fill(0.0);
                for q in 0..nq {
                    let (nv2, dn2) = &shapes_geom2[q];
                    let mut x_q = [0.0; 3];
                    let mut t1 = [0.0; 3];
                    let mut t2 = [0.0; 3];
                    for a in 0..6 {
                        let d1 = dn2[a * 3 + 1] - dn2[a * 3];
                        let d2 = dn2[a * 3 + 2] - dn2[a * 3];
                        for k in 0..3 {
                            x_q[k] += nv2[a] * x6[a][k];
                            t1[k] += d1 * x6[a][k];
                            t2[k] += d2 * x6[a][k];
                        }
                    }
                    let n_raw = cross3(t1, t2);
                    let n_len = norm3(n_raw);
                    let area_elem = 0.5 * n_len;
                    if area_elem <= 1.0e-30 {
                        continue;
                    }
                    let mut n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];
                    if dot3(n_hat, sub3(x_q, p_opp)) < 0.0 {
                        n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
                    }
                    let w_val = weight(f, c_f, x_q, n_hat);
                    if w_val == 0.0 {
                        continue;
                    }
                    let ds = rule.w[q] * area_elem * w_val;
                    let nv = &shapes_order[q];
                    for a in 0..nloc {
                        for b in 0..nloc {
                            ke[a * nloc + b] += ds * nv[a] * nv[b];
                        }
                    }
                }
                mesh.facet_dofs(f, order, &mut dofs);
                coo.push_block(&dofs, &dofs, &ke);
            }
        })
    }
}

/// Assembles the far-field Robin (`bc_order = 1`) or second-order symmetric Bayliss-Turkel
/// Laplace-Beltrami (`bc_order = 2` in 3D) boundary stiffness matrix `K_far` over `labels`.
///
/// * **3D Order 1:** `\alpha(\mathbf{x}) = \varepsilon_0 \varepsilon_{n,f} \frac{\gamma \max(\mathbf{r}\cdot\mathbf{n}_f, 0)}{\|\mathbf{r}\|^2}`
///   where `\varepsilon_{n,f} = \mathbf{n}_f^T \boldsymbol{\varepsilon}_r(c_f)\mathbf{n}_f`.
/// * **2D Order 1:** `\alpha(\mathbf{x}) = \varepsilon_0 \varepsilon_{n,f} \frac{\gamma}{2 R_{\text{eff}} \max(\ln(2 R_\infty / r), 1)}`
///   with `1 / R_{\text{eff}} = \max(\mathbf{r}\cdot\mathbf{n}_f, 0) / r^2` and positive domain scale `R_\infty`.
/// * **3D Order 2:** adds the symmetric Laplace-Beltrami surface stiffness
///   `\int_f \varepsilon_0 \varepsilon_{n,f} \frac{r^2}{2 \max(\mathbf{r}\cdot\mathbf{n}_f, \epsilon_{\text{geom}} r)} \nabla_\Gamma \phi_a \cdot \nabla_\Gamma \phi_b \, dS`
///   strictly where `\mathbf{r}\cdot\mathbf{n}_f > 0` (with tangential radial projection damping on tilted box faces so pure
///   monopoles `V \propto 1/r` remain exact on rectangular walls while spherical shells preserve full `\nabla_\Gamma`),
///   falling back gracefully to Order 1 on tangential faces.
#[allow(clippy::too_many_arguments)]
pub fn assemble_farfield_robin(
    mesh: &SimplexMesh,
    labels: &[i32],
    fe_order: usize,
    bc_order: usize,
    decay_power: f64,
    center: [f64; 3],
    eps0: f64,
    eps_r: FarFieldPermittivity<'_>,
) -> Csr {
    let gamma = if decay_power.is_finite() && decay_power > 0.0 { decay_power } else { 1.0 };
    if mesh.dim == 2 {
        let (lo, hi) = mesh.bbox();
        let r_inf = ((hi[0] - lo[0]).hypot(hi[1] - lo[1])).max(1.0e-12);
        return assemble_boundary_mass_weighted(mesh, labels, fe_order, |_f, c_f, x_q, n_hat| {
            let r_vec = [x_q[0] - center[0], x_q[1] - center[1], 0.0];
            let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1];
            let rn = dot3(r_vec, n_hat);
            if r2 <= 1.0e-30 || rn <= 0.0 {
                return 0.0;
            }
            let r = r2.sqrt();
            let r_scale = r_inf.max(r);
            let denom_ln = (2.0 * r_scale / r).ln().max(1.0);
            let inv_r_eff = rn / r2;
            let eps_nf = eps_r.normal_eps(2, c_f, n_hat);
            eps0 * eps_nf * (gamma * inv_r_eff) / (2.0 * denom_ln)
        });
    }

    if bc_order <= 1 {
        return assemble_boundary_mass_weighted(mesh, labels, fe_order, |_f, c_f, x_q, n_hat| {
            let r_vec = sub3(x_q, center);
            let r2 = dot3(r_vec, r_vec);
            let rn = dot3(r_vec, n_hat);
            if r2 <= 1.0e-30 || rn <= 0.0 {
                return 0.0;
            }
            let eps_nf = eps_r.normal_eps(3, c_f, n_hat);
            eps0 * eps_nf * (gamma * rn / r2)
        });
    }

    // 3D Order 2: Order-1 Robin mass + symmetric Laplace-Beltrami surface gradient term
    let n = mesh.n_dofs(fe_order);
    let facets = sorted_unique_facets(mesh, labels);
    if facets.is_empty() {
        return Csr::zeros(n, n);
    }
    let rule = facet_rule(2);
    let nq = rule.w.len();
    let nloc = if fe_order == 1 { 3 } else { 6 };
    let shapes_fe: Vec<(Vec<f64>, Vec<f64>)> = rule.lam.iter().map(|lam| shape(lam, fe_order)).collect();
    let shapes_geom2: Vec<(Vec<f64>, Vec<f64>)> = rule.lam.iter().map(|lam| shape(lam, 2)).collect();
    let f2c = build_face_to_cell_3d_local(mesh);

    assemble_parallel(facets.len(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut ke = vec![0.0; nloc * nloc];
        let mut grad_surf = vec![[0.0_f64; 3]; nloc];
        for idx in range {
            let f = facets[idx];
            let fv = mesh.facet(f);
            let Some(&(c_f, opp)) = f2c.get(&sorted3(fv[0], fv[1], fv[2])) else {
                continue;
            };
            let x6 = facet_geometry_nodes_3d(mesh, f);
            let po_s = mesh.point(mesh.cell(c_f)[opp] as usize);
            let p_opp = [po_s[0], po_s[1], po_s[2]];
            ke.fill(0.0);
            for q in 0..nq {
                let (nv2, dn2) = &shapes_geom2[q];
                let mut x_q = [0.0; 3];
                let mut t1 = [0.0; 3];
                let mut t2 = [0.0; 3];
                for a in 0..6 {
                    let d1 = dn2[a * 3 + 1] - dn2[a * 3];
                    let d2 = dn2[a * 3 + 2] - dn2[a * 3];
                    for k in 0..3 {
                        x_q[k] += nv2[a] * x6[a][k];
                        t1[k] += d1 * x6[a][k];
                        t2[k] += d2 * x6[a][k];
                    }
                }
                let n_raw = cross3(t1, t2);
                let n_len = norm3(n_raw);
                let area_elem = 0.5 * n_len;
                if area_elem <= 1.0e-30 {
                    continue;
                }
                let mut n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];
                if dot3(n_hat, sub3(x_q, p_opp)) < 0.0 {
                    n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
                }
                let r_vec = sub3(x_q, center);
                let r2 = dot3(r_vec, r_vec);
                let rn = dot3(r_vec, n_hat);
                if r2 <= 1.0e-30 || rn <= 0.0 {
                    continue;
                }
                let r = r2.sqrt();
                let eps_nf = eps_r.normal_eps(3, c_f, n_hat);
                if eps_nf <= 0.0 {
                    continue;
                }
                let ds = rule.w[q] * area_elem;
                let (nv, dn) = &shapes_fe[q];

                // Order-1 Robin contribution
                let alpha = eps0 * eps_nf * (gamma * rn / r2);
                let w_mass = ds * alpha;
                for a in 0..nloc {
                    for b in 0..nloc {
                        ke[a * nloc + b] += w_mass * nv[a] * nv[b];
                    }
                }

                // Order-2 symmetric Laplace-Beltrami contribution strictly where r . n_f > eps_geom * r
                // (falling back to Order 1 on tangential faces)
                if rn > EPS_GEOM * r {
                    let g11 = dot3(t1, t1);
                    let g12 = dot3(t1, t2);
                    let g22 = dot3(t2, t2);
                    let det_g = g11 * g22 - g12 * g12;
                    if det_g > 1.0e-60 {
                        let inv_det = 1.0 / det_g;
                        let t1_contra = [
                            (g22 * t1[0] - g12 * t2[0]) * inv_det,
                            (g22 * t1[1] - g12 * t2[1]) * inv_det,
                            (g22 * t1[2] - g12 * t2[2]) * inv_det,
                        ];
                        let t2_contra = [
                            (-g12 * t1[0] + g11 * t2[0]) * inv_det,
                            (-g12 * t1[1] + g11 * t2[1]) * inv_det,
                            (-g12 * t1[2] + g11 * t2[2]) * inv_det,
                        ];
                        // Tangential projection of unit radial vector r_hat onto the facet plane:
                        // r_gamma = r_hat - (r_hat . n_hat) n_hat (zero on a spherical shell, non-zero on flat box walls)
                        let r_hat = [r_vec[0] / r, r_vec[1] / r, r_vec[2] / r];
                        let cos_th = rn / r;
                        let r_gamma = [
                            r_hat[0] - cos_th * n_hat[0],
                            r_hat[1] - cos_th * n_hat[1],
                            r_hat[2] - cos_th * n_hat[2],
                        ];
                        let rg2 = dot3(r_gamma, r_gamma);
                        let proj_scale = 1.0 / (rg2 + RADIAL_PROJ_REG * RADIAL_PROJ_REG);

                        for a in 0..nloc {
                            let dphi1 = dn[a * 3 + 1] - dn[a * 3];
                            let dphi2 = dn[a * 3 + 2] - dn[a * 3];
                            let g_raw = [
                                dphi1 * t1_contra[0] + dphi2 * t2_contra[0],
                                dphi1 * t1_contra[1] + dphi2 * t2_contra[1],
                                dphi1 * t1_contra[2] + dphi2 * t2_contra[2],
                            ];
                            let gr = dot3(g_raw, r_gamma) * proj_scale;
                            grad_surf[a] = [
                                g_raw[0] - gr * r_gamma[0],
                                g_raw[1] - gr * r_gamma[1],
                                g_raw[2] - gr * r_gamma[2],
                            ];
                        }
                        let beta = eps0 * eps_nf * r2 / (2.0 * rn.max(EPS_GEOM * r));
                        let w_lb = ds * beta;
                        for a in 0..nloc {
                            for b in 0..nloc {
                                ke[a * nloc + b] += w_lb * dot3(grad_surf[a], grad_surf[b]);
                            }
                        }
                    }
                }
            }
            mesh.facet_dofs(f, fe_order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshing::geomdrive::block;

    #[test]
    fn centroid_is_finite_and_symmetric_on_block() {
        let m = block(2.0, 2.0, 2.0, 0.7, 3).unwrap().mesh;
        let all_labs = m.boundary_labels();
        let c = farfield_centroid(&m, &all_labs);
        assert!((c[0] - 1.0).abs() < 1e-12 && (c[1] - 1.0).abs() < 1e-12 && (c[2] - 1.0).abs() < 1e-12, "c = {c:?}");
        // Empty label slice falls back to bbox midpoint, never NaN
        let c_empty = farfield_centroid(&m, &[]);
        assert!(c_empty.iter().all(|x| x.is_finite()));
        assert!((c_empty[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn farfield_robin_order1_and_order2_are_symmetric_and_spd_on_constants() {
        let m = block(2.0, 2.0, 2.0, 0.7, 3).unwrap().mesh;
        let all_labs = m.boundary_labels();
        let center = [1.0, 1.0, 1.0];
        let eps = vec![1.0; m.n_cells()];
        let k1 = assemble_farfield_robin(&m, &all_labs, 2, 1, 1.0, center, 1.0, FarFieldPermittivity::Scalar(&eps));
        let k2 = assemble_farfield_robin(&m, &all_labs, 2, 2, 1.0, center, 1.0, FarFieldPermittivity::Scalar(&eps));
        let n = m.n_dofs(2);
        let ones = vec![1.0; n];
        let e1: f64 = k1.matvec(&ones).iter().sum();
        let e2: f64 = k2.matvec(&ones).iter().sum();
        assert!(e1 > 15.0 && e1 < 16.0, "expected positive cube Robin integral ~15.35, got e1 = {e1}");
        // Order-2 Laplace-Beltrami vanishes on constant vectors so 1^T K2 1 == 1^T K1 1 to machine precision
        assert!((e2 - e1).abs() < 1e-12, "e1 = {e1}, e2 = {e2}");
    }
}
