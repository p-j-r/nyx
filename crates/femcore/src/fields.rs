// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Derived fields and what is drawn: element gradients at the element's own nodes, the node averaging of nodal data
//! and `Display`, the surface mesh used for plotting.
//!
//! `gradients_at_nodes` evaluates `grad u` of a P1/P2 field at the nodes of every cell (so stresses can be formed per
//! cell and then averaged over the cells that share a node). `Display` keeps only the boundary facets
//! (3D) or the cells (2D); for P2 results each triangle is split in four so the midpoint values show.

use crate::fem::{geometry, shape, CurvedTet};
use crate::mesh::SimplexMesh;
use rayon::prelude::*;

/// Barycentric coordinates of the local dofs (vertices, then edge midpoints).
pub fn local_bary(dim: usize, order: usize) -> Vec<Vec<f64>> {
    let n = dim + 1;
    let mut pts: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| (i == j) as u8 as f64).collect()).collect();
    if order == 2 {
        for i in 0..n {
            for j in i + 1..n {
                pts.push((0..n).map(|k| if k == i || k == j { 0.5 } else { 0.0 }).collect());
            }
        }
    }
    pts
}

/// Gradients of `ncomp` stacked scalar dof fields (each `n_dofs` long) at the local dofs of
/// every cell: `out[((c*nloc + a)*ncomp + i)*d + k] = d u_i / d x_k`.
///
/// On 3D curved/quarter-point cells (`mesh.dim == 3 && mesh.cell_is_curved(c)`), evaluates [`CurvedTet::jacobian`],
/// clamping barycentric coordinates `lam[k] = lam[k].max(1e-3)` and renormalizing `sum(lam) == 1` whenever
/// `|det J| < 1e-12 * |det_straight|` so nodal gradients at singular quarter-point tips (`xi_QP = 0.25`, `det J = 0`)
/// stay strictly finite (`0` `NaN`/`Inf`).
pub fn gradients_at_nodes(mesh: &SimplexMesh, u: &[f64], ncomp: usize, order: usize) -> Vec<f64> {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    debug_assert_eq!(u.len(), ncomp * n1);
    let geo = geometry(mesh);
    let lam = local_bary(d, order);
    let nloc = lam.len();
    let dn: Vec<Vec<f64>> = lam.iter().map(|l| shape(l, order).1).collect(); // dn[point][a*(d+1)+b]
    let per_cell = nloc * ncomp * d;
    let mut out = vec![0.0; mesh.n_cells() * per_cell];
    let nt = if mesh.n_cells() < 2000 { 1 } else { crate::par::threads() };
    let chunk = mesh.n_cells().div_ceil(nt).max(1);
    let (geo, dn, lam) = (&geo, &dn, &lam);
    out.par_chunks_mut((chunk * per_cell).max(1)).enumerate().for_each(|(ci, block)| {
        let mut dofs = Vec::new();
        let mut s = vec![0.0; ncomp * d];
        for (lc, cell_out) in block.chunks_mut(per_cell).enumerate() {
            let c = ci * chunk + lc;
            mesh.cell_dofs(c, order, &mut dofs);
            let base = c * (d + 1) * d;
            if d == 3 && mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                let v = [ct.x[0], ct.x[1], ct.x[2], ct.x[3]];
                let det_straight = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
                let thresh = 1e-12 * det_straight.abs();
                for p in 0..nloc {
                    s.iter_mut().for_each(|x| *x = 0.0);
                    let lp = [lam[p][0], lam[p][1], lam[p][2], lam[p][3]];
                    let (mut g_bary, mut det) = ct.jacobian(&lp);
                    let clamped_dn;
                    let dn_p: &[f64] = if !det.is_finite()
                        || det.abs() < thresh
                        || g_bary.iter().flatten().any(|x| !x.is_finite())
                    {
                        let mut lc = [lp[0].max(1e-3), lp[1].max(1e-3), lp[2].max(1e-3), lp[3].max(1e-3)];
                        let sum_l: f64 = lc.iter().sum();
                        for x in &mut lc {
                            *x /= sum_l;
                        }
                        let (gc, det_c) = ct.jacobian(&lc);
                        g_bary = gc;
                        det = det_c;
                        clamped_dn = shape(&lc, order).1;
                        &clamped_dn
                    } else {
                        &dn[p]
                    };
                    if !det.is_finite() || det.abs() <= 1e-30 || g_bary.iter().flatten().any(|x| !x.is_finite()) {
                        for b in 0..4 {
                            for k in 0..3 {
                                g_bary[b][k] = geo.g[base + b * 3 + k];
                            }
                        }
                    }
                    for a in 0..nloc {
                        let mut ga = [0.0; 3];
                        for k in 0..3 {
                            for b in 0..4 {
                                ga[k] += dn_p[a * 4 + b] * g_bary[b][k];
                            }
                        }
                        for i in 0..ncomp {
                            let ua = u[i * n1 + dofs[a]];
                            for k in 0..3 {
                                s[i * 3 + k] += ua * ga[k];
                            }
                        }
                    }
                    cell_out[p * ncomp * d..(p + 1) * ncomp * d].copy_from_slice(&s);
                }
            } else {
                for p in 0..nloc {
                    s.iter_mut().for_each(|x| *x = 0.0);
                    for a in 0..nloc {
                        let mut ga = [0.0; 3];
                        for k in 0..d {
                            for b in 0..=d {
                                ga[k] += dn[p][a * (d + 1) + b] * geo.g[base + b * d + k];
                            }
                        }
                        for i in 0..ncomp {
                            let ua = u[i * n1 + dofs[a]];
                            for k in 0..d {
                                s[i * d + k] += ua * ga[k];
                            }
                        }
                    }
                    cell_out[p * ncomp * d..(p + 1) * ncomp * d].copy_from_slice(&s);
                }
            }
        }
    });
    out
}

/// Mean over the cells touching each dof of per-cell nodal values `vals[(c*nloc + a)*m + j]`.
pub fn node_average(mesh: &SimplexMesh, vals: &[f64], m: usize, order: usize) -> Vec<f64> {
    let n = mesh.n_dofs(order);
    let nloc = mesh.nloc(order);
    let mut acc = vec![0.0; n * m];
    let mut cnt = vec![0.0f64; n];
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        for a in 0..nloc {
            let d = dofs[a];
            cnt[d] += 1.0;
            for j in 0..m {
                acc[d * m + j] += vals[(c * nloc + a) * m + j];
            }
        }
    }
    for d in 0..n {
        if cnt[d] > 0.0 {
            for j in 0..m {
                acc[d * m + j] /= cnt[d];
            }
        }
    }
    acc
}

/// Surface mesh for plotting: the cells in 2D, the boundary facets in 3D, with every P2
/// triangle split in four so quadratic fields show their midpoint values. Only dofs used by
/// the surface are kept (`dofs` maps display points to dof indices).
#[derive(Clone, Debug)]
pub struct Display {
    pub dim: usize,
    pub dofs: Vec<usize>,
    /// 3 coordinates per display point
    pub points: Vec<[f32; 3]>,
    pub tris: Vec<[u32; 3]>,
    /// boundary label (3D) or domain label (2D) of each display triangle
    pub tri_label: Vec<i32>,
}

impl Display {
    pub fn new(mesh: &SimplexMesh, order: usize) -> Self {
        let d = mesh.dim;
        let (base, labels): (Vec<[u32; 3]>, &Vec<i32>) = if d == 2 {
            (mesh.cells.chunks(3).map(|c| [c[0], c[1], c[2]]).collect(), &mesh.cell_labels)
        } else {
            (mesh.facets.chunks(3).map(|c| [c[0], c[1], c[2]]).collect(), &mesh.facet_labels)
        };
        let nv = mesh.n_points();
        let mut tris_dof: Vec<[usize; 3]> = vec![];
        let mut tri_label = vec![];
        for (t, tri) in base.iter().enumerate() {
            if order == 2 {
                let e = |a: u32, b: u32| nv + mesh.edge_index(a, b);
                let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
                let (e01, e02, e12) = (e(tri[0], tri[1]), e(tri[0], tri[2]), e(tri[1], tri[2]));
                for s in [[a, e01, e02], [e01, b, e12], [e02, e12, c], [e01, e12, e02]] {
                    tris_dof.push(s);
                    tri_label.push(labels[t]);
                }
            } else {
                tris_dof.push([tri[0] as usize, tri[1] as usize, tri[2] as usize]);
                tri_label.push(labels[t]);
            }
        }
        let mut used: Vec<usize> = tris_dof.iter().flatten().copied().collect();
        used.sort_unstable();
        used.dedup();
        let mut map = std::collections::HashMap::with_capacity(used.len());
        for (i, &d) in used.iter().enumerate() {
            map.insert(d, i as u32);
        }
        let co = mesh.dof_coordinates(order);
        let points = used
            .iter()
            .map(|&dof| {
                let mut p = [0f32; 3];
                for k in 0..d {
                    p[k] = co[dof * d + k] as f32;
                }
                p
            })
            .collect();
        let tris = tris_dof.iter().map(|t| [map[&t[0]], map[&t[1]], map[&t[2]]]).collect();
        Display { dim: d, dofs: used, points, tris, tri_label }
    }

    pub fn pick(&self, values: &[f64]) -> Vec<f32> {
        self.dofs.iter().map(|&d| values[d] as f32).collect()
    }
}

// ------------------------------------------------------------------ per-cell nodal fields (the "fields" stage)

/// One quantity at the local dofs of every cell, not averaged over the cells that share a node:
/// `re[c * nloc + a]` (and `im` when the quantity is complex, else empty). Names are wafer names (`Ex`, `uz`, `Sxz`, ...).
#[derive(Clone, Debug, Default)]
pub struct CellField {
    pub name: String,
    pub unit: String,
    pub re: Vec<f64>,
    pub im: Vec<f64>,
}

/// The un-averaged per-cell fields of one solution. This is what an error estimator needs (the jump of a flux
/// between neighbouring cells is exactly what averaging destroys); `Solution::fields` are the node-averaged plot fields.
#[derive(Clone, Debug, Default)]
pub struct CellFields {
    pub order: usize,
    /// local dofs per cell (vertices, then edge midpoints for P2)
    pub nloc: usize,
    pub fields: Vec<CellField>,
}

impl CellFields {
    pub fn new(mesh: &SimplexMesh, order: usize) -> Self {
        CellFields { order, nloc: mesh.nloc(order), fields: vec![] }
    }

    pub fn get(&self, name: &str) -> Option<&CellField> {
        self.fields.iter().find(|f| f.name == name)
    }

    pub fn push_real(&mut self, name: impl Into<String>, unit: impl Into<String>, re: Vec<f64>) {
        self.fields.push(CellField { name: name.into(), unit: unit.into(), re, im: vec![] });
    }

    pub fn push_complex(&mut self, name: impl Into<String>, unit: impl Into<String>, re: Vec<f64>, im: Vec<f64>) {
        self.fields.push(CellField { name: name.into(), unit: unit.into(), re, im });
    }

    /// The `ncomp` interleaved per-cell values `g[(c*nloc + a)*ncomp + j]` as one real field per component, named
    /// `{prefix}{names[j]}`.
    pub fn push_components(&mut self, prefix: &str, names: &[&str], unit: &str, g: &[f64]) {
        let m = names.len();
        for (j, n) in names.iter().enumerate() {
            let re = g.chunks(m).map(|p| p[j]).collect();
            self.push_real(format!("{prefix}{n}"), unit, re);
        }
    }

    /// The node average of a real component (the plot field it must reproduce).
    pub fn node_average(&self, mesh: &SimplexMesh, name: &str) -> Option<Vec<f64>> {
        self.get(name).map(|f| node_average(mesh, &f.re, 1, self.order))
    }
}
