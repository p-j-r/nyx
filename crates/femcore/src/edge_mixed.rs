// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Conforming mixed-cell Nedelec `H(curl)` finite-element spaces (`Order 1` and `Order 2`) on
//! 2D (`Tri` + `Quad`) and 3D (`Tet` + `Prism` + `Hex`) meshes (`MixedDiscretisation`), covariant
//! Piola mappings on straight and curved isoparametric cells, and exact algebraic discrete-gradient
//! (`gradient_p1`, `gradient_p2`) and `P_1 -> P_2` inclusion (`p1_to_p2`) operators for `Ams` and `AmsP2`.
//!
//! Element local dimensions:
//! - 2D: `Tri` (`3` / `8` DOFs), `Quad` (`4` / `12` DOFs, `Q_{0,1} x Q_{1,0}` and `Q_{1,2} x Q_{2,1}`).
//! - 3D: `Tet` (`6` / `20` DOFs, identical to `nyx_ems::cavity::Ned3`),
//!   `Prism` (`9` / `36` DOFs = `18` edge + `4` tri-face + `12` quad-face + `2` interior,
//!   exact cochain space `(N_2(T) \otimes Q_2(I)) \oplus (P_2(T) \otimes Q_1(I) d\zeta)`),
//!   `Hex` (`12` / `54` DOFs = `24` edge + `24` quad-face + `6` interior, `Q_{1,2,2} x Q_{2,1,2} x Q_{2,2,1}`).

use crate::fem::assemble_parallel;
use crate::fem_mixed::{cell_corner_ref_coords, quad_rule_mixed, IsoparametricCell, JacobianEval};
use crate::linalg::ams::AmsP2Topology;
use crate::linalg::csr::{c, CCoo, CCsr, Coo, RCsr as Csr, C};
use crate::mesh::cells::CellType;
use crate::mesh::mixed::{
    d3_apply, d3_perm_for_keys, d4_apply, d4_dihedral_for_keys, FacetType, MixedDiscretisation, HEX_EDGES,
    HEX_QUAD_FACES, PRISM_QUAD_FACES, PRISM_TRI_FACES, QUAD_EDGES, TET_EDGES, TET_TRI_FACES, TRI_EDGES,
};
use rayon::prelude::*;
use std::collections::BTreeSet;

/// Number of local `H(curl)` basis functions on `cell_type` at polynomial `order` (`1` or `2`).
#[inline]
pub fn nedelec_nloc(cell_type: CellType, order: usize) -> usize {
    match (cell_type, order == 1) {
        (CellType::Tri, true) => 3,
        (CellType::Tri, false) => 8,
        (CellType::Quad, true) => 4,
        (CellType::Quad, false) => 12,
        (CellType::Tet, true) => 6,
        (CellType::Tet, false) => 20,
        (CellType::Prism, true) => 9,
        (CellType::Prism, false) => 36,
        (CellType::Hex, true) => 12,
        (CellType::Hex, false) => 54,
    }
}

/// Global conforming Nedelec `H(curl)` DOF numbering over a [`MixedDiscretisation`] (2D or 3D, order 1 or 2).
#[derive(Clone, Debug)]
pub struct NedMixed {
    pub dim: usize,
    pub order: usize,
    pub n_edges: usize,
    pub n_tri_faces: usize,
    pub n_quad_faces: usize,
    pub n_dofs: usize,
    /// Ordering key per vertex (`vkey[v] == v` when non-periodic; periodic vertex pairs share the class root).
    pub vkey: Vec<u32>,
    /// CSR offsets into `cell_dofs_flat` (`len = n_cells + 1`).
    pub cell_dof_offsets: Vec<usize>,
    /// Flat global DOFs of every cell in canonical local basis order.
    pub cell_dofs_flat: Vec<usize>,
    /// Starting global DOF of cell `c`'s interior bubble block (`usize::MAX` if cell `c` has 0 interior bubbles).
    pub cell_interior_start: Vec<usize>,
}

impl NedMixed {
    /// Builds a `NedMixed` numbering with default vertex keys (`0..n_points`).
    pub fn new(mesh: &MixedDiscretisation, order: usize) -> Self {
        let vkey: Vec<u32> = (0..mesh.n_points() as u32).collect();
        Self::with_vkey(mesh, order, vkey)
    }

    /// Builds a `NedMixed` numbering with explicit vertex ordering keys `vkey` (supporting Floquet periodicity).
    pub fn with_vkey(mesh: &MixedDiscretisation, order: usize, vkey: Vec<u32>) -> Self {
        assert!(order == 1 || order == 2, "NedMixed supports order 1 and order 2");
        assert_eq!(vkey.len(), mesh.n_points(), "vkey length must match n_points");
        let nc = mesh.n_cells();
        let ne = mesh.n_edges();
        let ntf = mesh.n_tri_faces();
        let nqf = mesh.n_quad_faces();

        let mut cell_interior_start = vec![usize::MAX; nc];
        let n_dofs = if order == 1 {
            ne
        } else if mesh.dim == 2 {
            let mut cur = 2 * ne;
            for c in 0..nc {
                cell_interior_start[c] = cur;
                cur += match mesh.cell_type(c) {
                    CellType::Tri => 2,
                    CellType::Quad => 4,
                    other => unreachable!("unexpected 2D cell {other:?}"),
                };
            }
            cur
        } else {
            let mut cur = 2 * ne + 2 * ntf + 4 * nqf;
            for c in 0..nc {
                match mesh.cell_type(c) {
                    CellType::Tet => {}
                    CellType::Prism => {
                        cell_interior_start[c] = cur;
                        cur += 2;
                    }
                    CellType::Hex => {
                        cell_interior_start[c] = cur;
                        cur += 6;
                    }
                    other => unreachable!("unexpected 3D cell {other:?}"),
                }
            }
            cur
        };

        let mut cell_dof_offsets = Vec::with_capacity(nc + 1);
        cell_dof_offsets.push(0);
        let mut cell_dofs_flat = Vec::with_capacity(nc * if order == 1 { 8 } else { 32 });

        let tf_base = 2 * ne;
        let qf_base = 2 * ne + 2 * ntf;

        for c in 0..nc {
            let ct = mesh.cell_type(c);
            let eids = mesh.cell_edge_ids(c);
            if order == 1 {
                for &e in eids {
                    cell_dofs_flat.push(e as usize);
                }
            } else {
                for &e in eids {
                    let eu = e as usize;
                    cell_dofs_flat.push(2 * eu);
                    cell_dofs_flat.push(2 * eu + 1);
                }
                if mesh.dim == 2 {
                    let base = cell_interior_start[c];
                    let n_int = if ct == CellType::Tri { 2 } else { 4 };
                    for k in 0..n_int {
                        cell_dofs_flat.push(base + k);
                    }
                } else {
                    for &tf in mesh.cell_tri_face_ids(c) {
                        let tfu = tf as usize;
                        cell_dofs_flat.push(tf_base + 2 * tfu);
                        cell_dofs_flat.push(tf_base + 2 * tfu + 1);
                    }
                    for &qf in mesh.cell_quad_face_ids(c) {
                        let qfu = qf as usize;
                        for k in 0..4 {
                            cell_dofs_flat.push(qf_base + 4 * qfu + k);
                        }
                    }
                    match ct {
                        CellType::Tet => {}
                        CellType::Prism => {
                            let base = cell_interior_start[c];
                            cell_dofs_flat.push(base);
                            cell_dofs_flat.push(base + 1);
                        }
                        CellType::Hex => {
                            let base = cell_interior_start[c];
                            for k in 0..6 {
                                cell_dofs_flat.push(base + k);
                            }
                        }
                        _ => unreachable!(),
                    }
                }
            }
            debug_assert_eq!(
                cell_dofs_flat.len() - *cell_dof_offsets.last().unwrap(),
                nedelec_nloc(ct, order)
            );
            cell_dof_offsets.push(cell_dofs_flat.len());
        }

        Self {
            dim: mesh.dim,
            order,
            n_edges: ne,
            n_tri_faces: ntf,
            n_quad_faces: nqf,
            n_dofs,
            vkey,
            cell_dof_offsets,
            cell_dofs_flat,
            cell_interior_start,
        }
    }

    /// Global DOFs of cell `c`.
    #[inline]
    pub fn cell_dofs(&self, c: usize) -> &[usize] {
        &self.cell_dofs_flat[self.cell_dof_offsets[c]..self.cell_dof_offsets[c + 1]]
    }

    /// Composite vertex key `(vkey[v], v)` (uses `vkey` primarily and breaks ties deterministically).
    #[inline]
    pub fn key_of(&self, v: u32) -> (u32, u32) {
        (self.vkey[v as usize], v)
    }

    /// Tangential Nedelec DOFs owned by boundary/sheet facet `f` (edges + face bubbles).
    pub fn facet_dofs(&self, mesh: &MixedDiscretisation, f: usize) -> Vec<usize> {
        let fv = mesh.facet(f);
        let ft = mesh.facet_type(f);
        let mut out = Vec::new();
        let push_edge = |out: &mut Vec<usize>, e: usize| {
            if self.order == 1 {
                out.push(e);
            } else {
                out.push(2 * e);
                out.push(2 * e + 1);
            }
        };
        match ft {
            FacetType::Edge2 => {
                push_edge(&mut out, mesh.edge_index(fv[0], fv[1]).expect("facet edge"));
            }
            FacetType::Tri3 => {
                for &(i, j) in &TRI_EDGES {
                    push_edge(&mut out, mesh.edge_index(fv[i], fv[j]).expect("facet tri edge"));
                }
                if self.order == 2 {
                    let tf = mesh.tri_face_index([fv[0], fv[1], fv[2]]).expect("facet tri face");
                    let base = 2 * self.n_edges + 2 * tf;
                    out.push(base);
                    out.push(base + 1);
                }
            }
            FacetType::Quad4 => {
                for &(i, j) in &QUAD_EDGES {
                    push_edge(&mut out, mesh.edge_index(fv[i], fv[j]).expect("facet quad edge"));
                }
                if self.order == 2 {
                    let qf = mesh
                        .quad_face_index([fv[0], fv[1], fv[2], fv[3]])
                        .expect("facet quad face");
                    let base = 2 * self.n_edges + 2 * self.n_tri_faces + 4 * qf;
                    for k in 0..4 {
                        out.push(base + k);
                    }
                }
            }
        }
        out
    }

    /// Sorted unique tangential Nedelec DOFs on all facets carrying any of `labels`.
    pub fn boundary_dofs(&self, mesh: &MixedDiscretisation, labels: &[i32]) -> Vec<usize> {
        let mut set = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if labels.contains(&mesh.facet_labels[f]) {
                set.extend(self.facet_dofs(mesh, f));
            }
        }
        set.into_iter().collect()
    }

    /// Sorted unique tangential Nedelec DOFs on all exterior facets (`mesh.exterior_facets[f] == true`).
    pub fn exterior_pec_dofs(&self, mesh: &MixedDiscretisation) -> Vec<usize> {
        let mut set = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if mesh.exterior_facets[f] {
                set.extend(self.facet_dofs(mesh, f));
            }
        }
        set.into_iter().collect()
    }
}

// ------------------------------------------------------------------ reference Nedelec basis + Piola mapping

#[inline]
fn cross2(u: [f64; 2], v: [f64; 2]) -> f64 {
    u[0] * v[1] - u[1] * v[0]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn l_pm(s: f64, t: f64) -> (f64, f64) {
    (0.5 * (1.0 + s * t), 0.5 * s)
}

#[inline]
fn bubble_1d(t: f64) -> (f64, f64) {
    (0.25 * (1.0 - t * t), -0.5 * t)
}

/// Evaluates 2D reference Nedelec basis `w_ref[a] = [w_\xi, w_\eta]` and scalar curl `curl_ref[a]`
/// on `Tri` or `Quad` at reference point `xi = [\xi, \eta]`.
pub fn eval_ref_basis_2d(
    cell_type: CellType,
    order: usize,
    xi: &[f64],
    vkeys: &[(u32, u32)],
    w_ref: &mut Vec<[f64; 2]>,
    curl_ref: &mut Vec<f64>,
) {
    w_ref.clear();
    curl_ref.clear();
    match cell_type {
        CellType::Tri => {
            let (u, v) = (xi[0], xi[1]);
            let l = [1.0 - u - v, u, v];
            let g: [[f64; 2]; 3] = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            let whitney = |a: usize, b: usize| -> ([f64; 2], f64) {
                (
                    [l[a] * g[b][0] - l[b] * g[a][0], l[a] * g[b][1] - l[b] * g[a][1]],
                    2.0 * cross2(g[a], g[b]),
                )
            };
            for &(i, j) in &TRI_EDGES {
                let (a, b) = if vkeys[i] < vkeys[j] { (i, j) } else { (j, i) };
                let (w0, c0) = whitney(a, b);
                w_ref.push(w0);
                curl_ref.push(c0);
                if order == 2 {
                    w_ref.push([l[a] * g[b][0] + l[b] * g[a][0], l[a] * g[b][1] + l[b] * g[a][1]]);
                    curl_ref.push(0.0);
                }
            }
            if order == 2 {
                let perm = d3_perm_for_keys([vkeys[0], vkeys[1], vkeys[2]]);
                let [a, b, cc] = d3_apply(perm, [0usize, 1, 2]);
                let (nab, cab) = whitney(a, b);
                let (nbc, cbc) = whitney(b, cc);
                w_ref.push([l[cc] * nab[0], l[cc] * nab[1]]);
                curl_ref.push(cross2(g[cc], nab) + l[cc] * cab);
                w_ref.push([l[a] * nbc[0], l[a] * nbc[1]]);
                curl_ref.push(cross2(g[a], nbc) + l[a] * cbc);
            }
        }
        CellType::Quad => {
            let (u, v) = (xi[0], xi[1]);
            let (lm_u, _) = l_pm(-1.0, u);
            let (lp_u, _) = l_pm(1.0, u);
            let (lm_v, _) = l_pm(-1.0, v);
            let (lp_v, _) = l_pm(1.0, v);
            let (bu, dbu) = bubble_1d(u);
            let (bv, dbv) = bubble_1d(v);

            // Edge 0: (0 -> 1), v = -1, u: -1 -> +1
            let s0 = if vkeys[0] < vkeys[1] { 1.0 } else { -1.0 };
            w_ref.push([0.5 * s0 * lm_v, 0.0]);
            curl_ref.push(0.25 * s0);
            if order == 2 {
                w_ref.push([dbu * lm_v, -0.5 * bu]);
                curl_ref.push(0.0);
            }

            // Edge 1: (1 -> 2), u = +1, v: -1 -> +1
            let s1 = if vkeys[1] < vkeys[2] { 1.0 } else { -1.0 };
            w_ref.push([0.0, 0.5 * s1 * lp_u]);
            curl_ref.push(0.25 * s1);
            if order == 2 {
                w_ref.push([0.5 * bv, lp_u * dbv]);
                curl_ref.push(0.0);
            }

            // Edge 2: (2 -> 3), v = +1, u: +1 -> -1
            let s2 = if vkeys[2] < vkeys[3] { 1.0 } else { -1.0 };
            w_ref.push([-0.5 * s2 * lp_v, 0.0]);
            curl_ref.push(0.25 * s2);
            if order == 2 {
                w_ref.push([dbu * lp_v, 0.5 * bu]);
                curl_ref.push(0.0);
            }

            // Edge 3: (3 -> 0), u = -1, v: +1 -> -1
            let s3 = if vkeys[3] < vkeys[0] { 1.0 } else { -1.0 };
            w_ref.push([0.0, -0.5 * s3 * lm_u]);
            curl_ref.push(0.25 * s3);
            if order == 2 {
                w_ref.push([-0.5 * bv, lm_u * dbv]);
                curl_ref.push(0.0);
            }

            if order == 2 {
                // 4 interior functions: f0=(bv, 0), f1=(u*bv, 0), f2=(0, bu), f3=(0, v*bu)
                w_ref.push([bv, 0.0]);
                curl_ref.push(-dbv);
                w_ref.push([u * bv, 0.0]);
                curl_ref.push(-u * dbv);
                w_ref.push([0.0, bu]);
                curl_ref.push(dbu);
                w_ref.push([0.0, v * bu]);
                curl_ref.push(v * dbu);
            }
        }
        other => panic!("eval_ref_basis_2d called on 3D cell {other:?}"),
    }
}

/// Evaluates 3D reference Nedelec basis `w_ref[a]` and reference curl `curl_ref[a]`
/// on `Tet`, `Prism`, or `Hex` at reference point `xi = [\xi, \eta, \zeta]`.
pub fn eval_ref_basis_3d(
    cell_type: CellType,
    order: usize,
    xi: &[f64],
    vkeys: &[(u32, u32)],
    w_ref: &mut Vec<[f64; 3]>,
    curl_ref: &mut Vec<[f64; 3]>,
) {
    w_ref.clear();
    curl_ref.clear();
    match cell_type {
        CellType::Tet => {
            let (u, v, z) = (xi[0], xi[1], xi[2]);
            let lam = [1.0 - u - v - z, u, v, z];
            let g: [[f64; 3]; 4] = [
                [-1.0, -1.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            let whitney = |a: usize, b: usize| -> ([f64; 3], [f64; 3]) {
                let w = [
                    lam[a] * g[b][0] - lam[b] * g[a][0],
                    lam[a] * g[b][1] - lam[b] * g[a][1],
                    lam[a] * g[b][2] - lam[b] * g[a][2],
                ];
                let cx = cross3(g[a], g[b]);
                (w, [2.0 * cx[0], 2.0 * cx[1], 2.0 * cx[2]])
            };
            for &(i, j) in &TET_EDGES {
                let (a, b) = if vkeys[i] < vkeys[j] { (i, j) } else { (j, i) };
                let (n, cn) = whitney(a, b);
                w_ref.push(n);
                curl_ref.push(cn);
                if order == 2 {
                    w_ref.push([
                        lam[a] * g[b][0] + lam[b] * g[a][0],
                        lam[a] * g[b][1] + lam[b] * g[a][1],
                        lam[a] * g[b][2] + lam[b] * g[a][2],
                    ]);
                    curl_ref.push([0.0; 3]);
                }
            }
            if order == 2 {
                for f in &TET_TRI_FACES {
                    let perm = d3_perm_for_keys([vkeys[f[0]], vkeys[f[1]], vkeys[f[2]]]);
                    let [a, b, cc] = d3_apply(perm, *f);
                    let (nab, cab) = whitney(a, b);
                    let (nbc, cbc) = whitney(b, cc);
                    w_ref.push([lam[cc] * nab[0], lam[cc] * nab[1], lam[cc] * nab[2]]);
                    let x0 = cross3(g[cc], nab);
                    curl_ref.push([
                        lam[cc] * cab[0] + x0[0],
                        lam[cc] * cab[1] + x0[1],
                        lam[cc] * cab[2] + x0[2],
                    ]);
                    w_ref.push([lam[a] * nbc[0], lam[a] * nbc[1], lam[a] * nbc[2]]);
                    let x1 = cross3(g[a], nbc);
                    curl_ref.push([
                        lam[a] * cbc[0] + x1[0],
                        lam[a] * cbc[1] + x1[1],
                        lam[a] * cbc[2] + x1[2],
                    ]);
                }
            }
        }
        CellType::Prism => {
            let (u, v, z) = (xi[0], xi[1], xi[2]);
            let lam = [1.0 - u - v, u, v];
            let g2: [[f64; 2]; 3] = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
            let (lm_z, dlm_z) = l_pm(-1.0, z);
            let (lp_z, dlp_z) = l_pm(1.0, z);
            let (bz, dbz) = bubble_1d(z);

            let whitney2 = |a: usize, b: usize| -> ([f64; 2], f64) {
                (
                    [
                        lam[a] * g2[b][0] - lam[b] * g2[a][0],
                        lam[a] * g2[b][1] - lam[b] * g2[a][1],
                    ],
                    2.0 * cross2(g2[a], g2[b]),
                )
            };
            let grad_prod2 = |a: usize, b: usize| -> [f64; 2] {
                [
                    lam[a] * g2[b][0] + lam[b] * g2[a][0],
                    lam[a] * g2[b][1] + lam[b] * g2[a][1],
                ]
            };

            // 1. Base tri edges (0..3) & top tri edges (3..6)
            for (end_idx, &(lz, dlz)) in [(lm_z, dlm_z), (lp_z, dlp_z)].iter().enumerate() {
                let v_off = 3 * end_idx;
                for &(i, j) in &TRI_EDGES {
                    let (a, b) = if vkeys[v_off + i] < vkeys[v_off + j] {
                        (i, j)
                    } else {
                        (j, i)
                    };
                    let (nab, cab) = whitney2(a, b);
                    w_ref.push([lz * nab[0], lz * nab[1], 0.0]);
                    curl_ref.push([-dlz * nab[1], dlz * nab[0], lz * cab]);
                    if order == 2 {
                        let gp = grad_prod2(a, b);
                        w_ref.push([lz * gp[0], lz * gp[1], dlz * lam[a] * lam[b]]);
                        curl_ref.push([0.0; 3]);
                    }
                }
            }
            // 2. Vertical edges (6..9): (0,3), (1,4), (2,5)
            for i in 0..3 {
                let sgn = if vkeys[i] < vkeys[i + 3] { 1.0 } else { -1.0 };
                w_ref.push([0.0, 0.0, 0.5 * sgn * lam[i]]);
                curl_ref.push([0.5 * sgn * g2[i][1], -0.5 * sgn * g2[i][0], 0.0]);
                if order == 2 {
                    w_ref.push([bz * g2[i][0], bz * g2[i][1], lam[i] * dbz]);
                    curl_ref.push([0.0; 3]);
                }
            }

            if order == 2 {
                // 3. Base and top Tri faces (`PRISM_TRI_FACES = [[0,1,2], [3,4,5]]`)
                for (lf, &(lz, dlz)) in [(lm_z, dlm_z), (lp_z, dlp_z)].iter().enumerate() {
                    let f = PRISM_TRI_FACES[lf];
                    let perm = d3_perm_for_keys([vkeys[f[0]], vkeys[f[1]], vkeys[f[2]]]);
                    let [a, b, cc] = d3_apply(perm, [0usize, 1, 2]);
                    let (nab, cab) = whitney2(a, b);
                    let (nbc, cbc) = whitney2(b, cc);

                    let v0 = [lam[cc] * nab[0], lam[cc] * nab[1]];
                    let c0 = cross2(g2[cc], nab) + lam[cc] * cab;
                    w_ref.push([lz * v0[0], lz * v0[1], 0.0]);
                    curl_ref.push([-dlz * v0[1], dlz * v0[0], lz * c0]);

                    let v1 = [lam[a] * nbc[0], lam[a] * nbc[1]];
                    let c1 = cross2(g2[a], nbc) + lam[a] * cbc;
                    w_ref.push([lz * v1[0], lz * v1[1], 0.0]);
                    curl_ref.push([-dlz * v1[1], dlz * v1[0], lz * c1]);
                }

                // 4. Side Quad faces (`PRISM_QUAD_FACES = [[0,1,4,3], [1,2,5,4], [2,0,3,5]]`)
                const LOCAL_SZ: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
                for qf in &PRISM_QUAD_FACES {
                    let (i, j) = (qf[0], qf[1]);
                    let (nij, cij) = whitney2(i, j);
                    let gij = grad_prod2(i, j);
                    let lij = lam[i] * lam[j];

                    // Local lifts in (s, \zeta):
                    let wh0 = [2.0 * bz * nij[0], 2.0 * bz * nij[1], 0.0];
                    let ch0 = [-2.0 * dbz * nij[1], 2.0 * dbz * nij[0], 2.0 * bz * cij];

                    let wh1 = [-2.0 * bz * gij[0], -2.0 * bz * gij[1], 0.0];
                    let ch1 = [2.0 * dbz * gij[1], -2.0 * dbz * gij[0], 0.0];

                    let wv0 = [0.0, 0.0, lij];
                    let cv0 = [gij[1], -gij[0], 0.0];

                    let wv1 = [0.0, 0.0, lij * z];
                    let cv1 = [z * gij[1], -z * gij[0], 0.0];

                    let dih = d4_dihedral_for_keys([vkeys[qf[0]], vkeys[qf[1]], vkeys[qf[2]], vkeys[qf[3]]]);
                    let [ca, cb, _, cd] = d4_apply(dih, LOCAL_SZ);
                    let du = [0.5 * (cb.0 - ca.0), 0.5 * (cb.1 - ca.1)];
                    let dv = [0.5 * (cd.0 - ca.0), 0.5 * (cd.1 - ca.1)];

                    let scale = |s: f64, v: [f64; 3]| [s * v[0], s * v[1], s * v[2]];
                    if du[0] != 0.0 {
                        // u = sigma_s * s, v = sigma_z * z
                        let (ss, sz) = (du[0], dv[1]);
                        w_ref.push(scale(ss, wh0));
                        curl_ref.push(scale(ss, ch0));
                        w_ref.push(wh1);
                        curl_ref.push(ch1);
                        w_ref.push(scale(sz, wv0));
                        curl_ref.push(scale(sz, cv0));
                        w_ref.push(wv1);
                        curl_ref.push(cv1);
                    } else {
                        // u = sigma_z * z, v = sigma_s * s
                        let (sz, ss) = (du[1], dv[0]);
                        w_ref.push(scale(sz, wv0));
                        curl_ref.push(scale(sz, cv0));
                        w_ref.push(wv1);
                        curl_ref.push(cv1);
                        w_ref.push(scale(ss, wh0));
                        curl_ref.push(scale(ss, ch0));
                        w_ref.push(wh1);
                        curl_ref.push(ch1);
                    }
                }

                // 5. 2 Prism interior transverse bubbles: bz * l_2 * N_01 and bz * l_0 * N_12
                let (n01, c01) = whitney2(0, 1);
                let (n12, c12) = whitney2(1, 2);
                let v0 = [lam[2] * n01[0], lam[2] * n01[1]];
                let ci0 = cross2(g2[2], n01) + lam[2] * c01;
                w_ref.push([bz * v0[0], bz * v0[1], 0.0]);
                curl_ref.push([-dbz * v0[1], dbz * v0[0], bz * ci0]);

                let v1 = [lam[0] * n12[0], lam[0] * n12[1]];
                let ci1 = cross2(g2[0], n12) + lam[0] * c12;
                w_ref.push([bz * v1[0], bz * v1[1], 0.0]);
                curl_ref.push([-dbz * v1[1], dbz * v1[0], bz * ci1]);
            }
        }
        CellType::Hex => {
            #[inline]
            fn push_axis(
                w_ref: &mut Vec<[f64; 3]>,
                curl_ref: &mut Vec<[f64; 3]>,
                ax: usize,
                f: [(f64, f64); 3],
            ) {
                let val = f[0].0 * f[1].0 * f[2].0;
                let grad_f = [
                    f[0].1 * f[1].0 * f[2].0,
                    f[0].0 * f[1].1 * f[2].0,
                    f[0].0 * f[1].0 * f[2].1,
                ];
                let mut e_ax = [0.0; 3];
                e_ax[ax] = 1.0;
                let mut w = [0.0; 3];
                w[ax] = val;
                w_ref.push(w);
                curl_ref.push(cross3(grad_f, e_ax));
            }

            let corners = cell_corner_ref_coords(CellType::Hex);
            let c3 = [xi[0], xi[1], xi[2]];
            let b3 = [bubble_1d(c3[0]), bubble_1d(c3[1]), bubble_1d(c3[2])];

            // 1. 12 Edges (`HEX_EDGES`)
            for &(i, j) in &HEX_EDGES {
                let (ca, cb) = (corners[i], corners[j]);
                let ax = if ca[0] != cb[0] {
                    0
                } else if ca[1] != cb[1] {
                    1
                } else {
                    2
                };
                let (lo, hi) = if vkeys[i] < vkeys[j] { (ca, cb) } else { (cb, ca) };
                let dir = 0.5 * (hi[ax] - lo[ax]);
                let mut f_n1 = [(0.0, 0.0); 3];
                for a in 0..3 {
                    f_n1[a] = if a == ax {
                        (0.5 * dir, 0.0)
                    } else {
                        l_pm(ca[a], c3[a])
                    };
                }
                push_axis(w_ref, curl_ref, ax, f_n1);

                if order == 2 {
                    // Exact gradient of b(xi_ax) * l_{a1} * l_{a2}
                    let mut f_b = f_n1;
                    f_b[ax] = b3[ax];
                    w_ref.push([
                        f_b[0].1 * f_b[1].0 * f_b[2].0,
                        f_b[0].0 * f_b[1].1 * f_b[2].0,
                        f_b[0].0 * f_b[1].0 * f_b[2].1,
                    ]);
                    curl_ref.push([0.0; 3]);
                }
            }

            if order == 2 {
                // 2. 6 Quad faces (`HEX_QUAD_FACES`, 4 functions each)
                for qf in &HEX_QUAD_FACES {
                    let dih = d4_dihedral_for_keys([vkeys[qf[0]], vkeys[qf[1]], vkeys[qf[2]], vkeys[qf[3]]]);
                    let [ca, cb, _, cd] =
                        d4_apply(dih, [corners[qf[0]], corners[qf[1]], corners[qf[2]], corners[qf[3]]]);
                    let ax_u = (0..3).find(|&a| cb[a] != ca[a]).unwrap();
                    let ax_v = (0..3).find(|&a| cd[a] != ca[a]).unwrap();
                    let ax_n = 3 - ax_u - ax_v;
                    let su = 0.5 * (cb[ax_u] - ca[ax_u]);
                    let sv = 0.5 * (cd[ax_v] - ca[ax_v]);
                    let sn = ca[ax_n];
                    let ln = l_pm(sn, c3[ax_n]);

                    // Omega_0 = l_Q * b(v) * su * e_{ax_u}
                    let mut f0 = [(0.0, 0.0); 3];
                    f0[ax_u] = (su, 0.0);
                    f0[ax_v] = b3[ax_v];
                    f0[ax_n] = ln;
                    push_axis(w_ref, curl_ref, ax_u, f0);

                    // Omega_1 = l_Q * b(v) * xi_{ax_u} * e_{ax_u}
                    let mut f1 = [(0.0, 0.0); 3];
                    f1[ax_u] = (c3[ax_u], 1.0);
                    f1[ax_v] = b3[ax_v];
                    f1[ax_n] = ln;
                    push_axis(w_ref, curl_ref, ax_u, f1);

                    // Omega_2 = l_Q * b(u) * sv * e_{ax_v}
                    let mut f2 = [(0.0, 0.0); 3];
                    f2[ax_u] = b3[ax_u];
                    f2[ax_v] = (sv, 0.0);
                    f2[ax_n] = ln;
                    push_axis(w_ref, curl_ref, ax_v, f2);

                    // Omega_3 = l_Q * b(u) * xi_{ax_v} * e_{ax_v}
                    let mut f3 = [(0.0, 0.0); 3];
                    f3[ax_u] = b3[ax_u];
                    f3[ax_v] = (c3[ax_v], 1.0);
                    f3[ax_n] = ln;
                    push_axis(w_ref, curl_ref, ax_v, f3);
                }

                // 3. 6 Hex interior bubbles: for ax = 0, 1, 2 -> (1, b, b) and (xi_ax, b, b)
                for ax in 0..3 {
                    let mut fi0 = b3;
                    fi0[ax] = (1.0, 0.0);
                    push_axis(w_ref, curl_ref, ax, fi0);

                    let mut fi1 = b3;
                    fi1[ax] = (c3[ax], 1.0);
                    push_axis(w_ref, curl_ref, ax, fi1);
                }
            }
        }
        other => panic!("eval_ref_basis_3d called on 2D cell {other:?}"),
    }
}

/// Scratch workspace for allocation-free Nedelec basis + Jacobian evaluation inside quadrature loops.
#[derive(Default)]
pub struct NedEvalScratch {
    pub vkeys: Vec<(u32, u32)>,
    pub w_ref2: Vec<[f64; 2]>,
    pub curl_ref2: Vec<f64>,
    pub w_ref3: Vec<[f64; 3]>,
    pub curl_ref3: Vec<[f64; 3]>,
    pub shape_n: Vec<f64>,
    pub shape_dn: Vec<[f64; 3]>,
}

impl NedEvalScratch {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Evaluates physical 3D Nedelec basis `w[a]` and physical `curl[a]` on cell `ci` at `xi` using covariant Piola mapping.
#[allow(clippy::too_many_arguments)]
pub fn eval_cell_basis_3d(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    ci: usize,
    iso: &IsoparametricCell,
    xi: &[f64],
    scratch: &mut NedEvalScratch,
    w_out: &mut Vec<[f64; 3]>,
    curl_out: &mut Vec<[f64; 3]>,
) -> JacobianEval {
    scratch.vkeys.clear();
    for &v in mesh.cell(ci) {
        scratch.vkeys.push(ned.key_of(v));
    }
    let jac = iso.eval_into(xi, &mut scratch.shape_n, &mut scratch.shape_dn);
    eval_ref_basis_3d(
        iso.cell_type,
        ned.order,
        xi,
        &scratch.vkeys,
        &mut scratch.w_ref3,
        &mut scratch.curl_ref3,
    );
    w_out.clear();
    curl_out.clear();
    for (wr, cr) in scratch.w_ref3.iter().zip(&scratch.curl_ref3) {
        w_out.push(jac.covariant_piola(*wr));
        curl_out.push(jac.contravariant_curl_3d(*cr));
    }
    jac
}

/// Evaluates physical 2D Nedelec basis `w[a] = [w_x, w_z]` and scalar `curl[a]` on cell `ci` at `xi` using covariant Piola mapping.
#[allow(clippy::too_many_arguments)]
pub fn eval_cell_basis_2d(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    ci: usize,
    iso: &IsoparametricCell,
    xi: &[f64],
    scratch: &mut NedEvalScratch,
    w_out: &mut Vec<[f64; 2]>,
    curl_out: &mut Vec<f64>,
) -> JacobianEval {
    scratch.vkeys.clear();
    for &v in mesh.cell(ci) {
        scratch.vkeys.push(ned.key_of(v));
    }
    let jac = iso.eval_into(xi, &mut scratch.shape_n, &mut scratch.shape_dn);
    eval_ref_basis_2d(
        iso.cell_type,
        ned.order,
        xi,
        &scratch.vkeys,
        &mut scratch.w_ref2,
        &mut scratch.curl_ref2,
    );
    w_out.clear();
    curl_out.clear();
    for (wr, &cr) in scratch.w_ref2.iter().zip(&scratch.curl_ref2) {
        let wp = jac.covariant_piola([wr[0], wr[1], 0.0]);
        w_out.push([wp[0], wp[1]]);
        curl_out.push(jac.contravariant_curl_2d(cr));
    }
    jac
}

// ------------------------------------------------------------------ Periodic Floquet matching on MixedDiscretisation

/// Floquet-periodic Nedelec DOF relations `(dest_dof, src_dof, exp(-i phase))` across `Tri3`, `Quad4`, or `Edge2` facets.
pub fn ned_mixed_periodic_relations(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    src_facets: &[usize],
    vmap: &std::collections::HashMap<u32, u32>,
    phase: f64,
) -> Result<Vec<(usize, usize, C)>, String> {
    let f = crate::periodic::phase_factor(phase);
    let mut out = Vec::new();
    let mut done_edges = BTreeSet::new();
    let mut done_tris = BTreeSet::new();
    let mut done_quads = BTreeSet::new();

    for &fc in src_facets {
        let fv = mesh.facet(fc);
        let ft = mesh.facet_type(fc);
        let le: &[(usize, usize)] = match ft {
            FacetType::Edge2 => &[(0, 1)],
            FacetType::Tri3 => &TRI_EDGES,
            FacetType::Quad4 => &QUAD_EDGES,
        };
        for &(i, j) in le {
            let (asrc, bsrc) = (fv[i], fv[j]);
            let key = (asrc.min(bsrc), asrc.max(bsrc));
            if !done_edges.insert(key) {
                continue;
            }
            let es = mesh.edge_index(asrc, bsrc).unwrap();
            let (adst, bdst) = (vmap[&asrc], vmap[&bsrc]);
            let et = mesh.edge_index(adst, bdst).ok_or_else(|| {
                "an edge of the source periodic boundary has no counterpart on the destination boundary".to_string()
            })?;
            if ned.order == 1 {
                out.push((et, es, f));
            } else {
                out.push((2 * et, 2 * es, f));
                out.push((2 * et + 1, 2 * es + 1, f));
            }
        }
        if ned.order == 2 {
            match ft {
                FacetType::Edge2 => {}
                FacetType::Tri3 => {
                    let tfs = mesh.tri_face_index([fv[0], fv[1], fv[2]]).unwrap();
                    if done_tris.insert(tfs) {
                        let tft = mesh
                            .tri_face_index([vmap[&fv[0]], vmap[&fv[1]], vmap[&fv[2]]])
                            .ok_or_else(|| "periodic source Tri3 face has no destination counterpart".to_string())?;
                        let base = 2 * ned.n_edges;
                        out.push((base + 2 * tft, base + 2 * tfs, f));
                        out.push((base + 2 * tft + 1, base + 2 * tfs + 1, f));
                    }
                }
                FacetType::Quad4 => {
                    let qfs = mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap();
                    if done_quads.insert(qfs) {
                        let qft = mesh
                            .quad_face_index([vmap[&fv[0]], vmap[&fv[1]], vmap[&fv[2]], vmap[&fv[3]]])
                            .ok_or_else(|| "periodic source Quad4 face has no destination counterpart".to_string())?;
                        let base = 2 * ned.n_edges + 2 * ned.n_tri_faces;
                        for k in 0..4 {
                            out.push((base + 4 * qft + k, base + 4 * qfs + k, f));
                        }
                    }
                }
            }
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------ Exact algebraic de Rham & AMS/AmsP2 operators

/// Full discrete gradient `G_1` (`n_edges x n_points`): row `e` has `-1` at `lo` and `+1` at `hi` ordered by `(vkey[v], v)`.
pub fn gradient_p1_full(mesh: &MixedDiscretisation, vkey: &[u32]) -> Csr {
    let ne = mesh.n_edges();
    let nv = mesh.n_points();
    let mut coo = Coo::new(ne, nv);
    for (e, &(a, b)) in mesh.edges.iter().enumerate() {
        let (lo, hi) = if (vkey[a as usize], a) < (vkey[b as usize], b) {
            (a, b)
        } else {
            (b, a)
        };
        coo.push(e, hi as usize, 1.0);
        coo.push(e, lo as usize, -1.0);
    }
    coo.to_csr()
}

/// Edge tangent vectors `x_hi - x_lo` (`len = n_edges`) matching the row orientation of [`gradient_p1_full`].
pub fn edge_vectors_full(mesh: &MixedDiscretisation, vkey: &[u32]) -> Vec<[f64; 3]> {
    mesh.edges
        .iter()
        .map(|&(a, b)| {
            let (lo, hi) = if (vkey[a as usize], a) < (vkey[b as usize], b) {
                (a, b)
            } else {
                (b, a)
            };
            let (pa, pb) = (mesh.point3(lo as usize), mesh.point3(hi as usize));
            [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]]
        })
        .collect()
}

/// Exact `P_1 -> P_2` inclusion matrix (`ned2.n_dofs x n_edges`): embeds lowest-order `N_1` edge `e` at row `2 * e`.
pub fn p1_to_p2_full(mesh: &MixedDiscretisation, ned2: &NedMixed) -> Csr {
    assert_eq!(ned2.order, 2, "p1_to_p2_full requires order == 2");
    let ne = mesh.n_edges();
    let mut coo = Coo::new(ned2.n_dofs, ne);
    for e in 0..ne {
        coo.push(2 * e, e, 1.0);
    }
    coo.to_csr()
}

/// Number of scalar `H^1` DOFs in the exact order-2 de Rham partner space `S_2` (`Order2Tensor` hierarchical space).
pub fn scalar_p2_tensor_count(mesh: &MixedDiscretisation) -> usize {
    let nv = mesh.n_points();
    let ne = mesh.n_edges();
    if mesh.dim == 2 {
        let n_quads = mesh.cell_types.iter().filter(|&&t| t == CellType::Quad).count();
        nv + ne + n_quads
    } else {
        let n_hexes = mesh.cell_types.iter().filter(|&&t| t == CellType::Hex).count();
        nv + ne + mesh.n_quad_faces() + n_hexes
    }
}

/// Full exact discrete gradient `G_2` (`ned2.n_dofs x n_scalar_p2_tensor`) mapping the full tensor-product
/// `P_2` / `Q_2` hierarchical scalar `H^1` space (vertices + edge bubbles + quad-face bubbles + hex volume bubbles)
/// into `ker(curl) \subset V_{N_2}`, satisfying `K * G_2 == 0` identically on straight and curved mixed cells.
pub fn gradient_p2_full(mesh: &MixedDiscretisation, ned2: &NedMixed) -> Csr {
    assert_eq!(ned2.order, 2, "gradient_p2_full requires order == 2");
    let nv = mesh.n_points();
    let ne = mesh.n_edges();
    let ncols = scalar_p2_tensor_count(mesh);
    let mut coo = Coo::new(ned2.n_dofs, ncols);

    // 1. Vertex columns `0 .. nv` and Edge-bubble columns `nv .. nv + ne`
    for (e, &(a, b)) in mesh.edges.iter().enumerate() {
        let (lo, hi) = if ned2.key_of(a) < ned2.key_of(b) {
            (a, b)
        } else {
            (b, a)
        };
        coo.push(2 * e, hi as usize, 1.0);
        coo.push(2 * e, lo as usize, -1.0);
        coo.push(2 * e + 1, nv + e, 1.0);
    }

    if mesh.dim == 2 {
        // 2D Quad-cell center bubbles: \nabla(b(u)*b(v)) = -0.5 * f_1 - 0.5 * f_3
        let mut q_idx = 0usize;
        for c in 0..mesh.n_cells() {
            if mesh.cell_type(c) == CellType::Quad {
                let col = nv + ne + q_idx;
                let base = ned2.cell_interior_start[c];
                coo.push(base + 1, col, -0.5);
                coo.push(base + 3, col, -0.5);
                q_idx += 1;
            }
        }
    } else {
        // 3D Quad-face center bubbles `nv + ne .. nv + ne + n_quad_faces`
        let qf_dof_base = 2 * ne + 2 * ned2.n_tri_faces;
        for qf in 0..ned2.n_quad_faces {
            let col = nv + ne + qf;
            coo.push(qf_dof_base + 4 * qf + 1, col, -0.5);
            coo.push(qf_dof_base + 4 * qf + 3, col, -0.5);
        }
        // Normal derivative of \Phi_Q^{Hex} inside each adjacent Hex cell, plus Hex volume bubbles!
        let hex_corners = cell_corner_ref_coords(CellType::Hex);
        let mut hex_idx = 0usize;
        for c in 0..mesh.n_cells() {
            if mesh.cell_type(c) == CellType::Hex {
                let int_base = ned2.cell_interior_start[c];
                let qf_ids = mesh.cell_quad_face_ids(c);
                for (lq, &qf_u32) in qf_ids.iter().enumerate() {
                    let qf_col = nv + ne + qf_u32 as usize;
                    // Normal axis of HEX_QUAD_FACES[lq]:
                    // 0: [0,3,2,1] z=-1, 1: [4,5,6,7] z=+1, 2: [0,1,5,4] y=-1,
                    // 3: [1,2,6,5] x=+1, 4: [2,3,7,6] y=+1, 5: [3,0,4,7] x=-1
                    let lqf = HEX_QUAD_FACES[lq];
                    let ax_n = (0..3)
                        .find(|&a| {
                            let c0 = hex_corners[lqf[0]][a];
                            lqf.iter().all(|&vi| hex_corners[vi][a] == c0)
                        })
                        .unwrap();
                    let sn = hex_corners[lqf[0]][ax_n];
                    coo.push(int_base + 2 * ax_n, qf_col, 0.5 * sn);
                }

                // Hex volume bubble: \nabla(b(\xi)*b(\eta)*b(\zeta)) = -0.5 * (I_{\xi,1} + I_{\eta,1} + I_{\zeta,1})
                let hex_col = nv + ne + ned2.n_quad_faces + hex_idx;
                coo.push(int_base + 1, hex_col, -0.5);
                coo.push(int_base + 3, hex_col, -0.5);
                coo.push(int_base + 5, hex_col, -0.5);
                hex_idx += 1;
            }
        }
    }
    coo.to_csr()
}

/// Disjoint local DOF blocks in full `NedMixed` order-2 space (for block Gauss–Seidel smoothing in `AmsP2`).
pub fn ams_p2_blocks_full(mesh: &MixedDiscretisation, ned2: &NedMixed) -> Vec<Vec<usize>> {
    assert_eq!(ned2.order, 2);
    let mut blocks = Vec::new();
    for e in 0..ned2.n_edges {
        blocks.push(vec![2 * e, 2 * e + 1]);
    }
    if mesh.dim == 2 {
        for c in 0..mesh.n_cells() {
            let base = ned2.cell_interior_start[c];
            let m = if mesh.cell_type(c) == CellType::Tri { 2 } else { 4 };
            blocks.push((base..base + m).collect());
        }
    } else {
        let tf_base = 2 * ned2.n_edges;
        for tf in 0..ned2.n_tri_faces {
            blocks.push(vec![tf_base + 2 * tf, tf_base + 2 * tf + 1]);
        }
        let qf_base = tf_base + 2 * ned2.n_tri_faces;
        for qf in 0..ned2.n_quad_faces {
            blocks.push((0..4).map(|k| qf_base + 4 * qf + k).collect());
        }
        for c in 0..mesh.n_cells() {
            match mesh.cell_type(c) {
                CellType::Tet => {}
                CellType::Prism => {
                    let base = ned2.cell_interior_start[c];
                    blocks.push(vec![base, base + 1]);
                }
                CellType::Hex => {
                    let base = ned2.cell_interior_start[c];
                    blocks.push((base..base + 6).collect());
                }
                _ => {}
            }
        }
    }
    blocks
}

/// Constrained lowest-order (`N_1`) auxiliary-space inputs ready for [`crate::linalg::ams::Ams::new`].
#[derive(Clone, Debug)]
pub struct AmsP1Data {
    /// Free `N_1` edge indices in `0..n_edges`.
    pub free_edges: Vec<usize>,
    /// Free scalar vertex indices in `0..n_points`.
    pub free_vertices: Vec<usize>,
    /// Interior discrete gradient `G_1` (`n_free_edges x n_free_vertices`).
    pub gradient_p1: Csr,
    /// Edge tangent vectors `x_hi - x_lo` (`len = n_free_edges`).
    pub edge_vectors: Vec<[f64; 3]>,
}

/// Builds [`AmsP1Data`] after eliminating constrained PEC facets (`pec_facets` are indices in `0..mesh.n_facets()`).
pub fn build_ams_p1(mesh: &MixedDiscretisation, vkey: &[u32], pec_facets: &[usize]) -> AmsP1Data {
    let mut fixed_edges = vec![false; mesh.n_edges()];
    let mut fixed_verts = vec![false; mesh.n_points()];
    for &f in pec_facets {
        let fv = mesh.facet(f);
        for &v in fv {
            fixed_verts[v as usize] = true;
        }
        let le: &[(usize, usize)] = match mesh.facet_type(f) {
            FacetType::Edge2 => &[(0, 1)],
            FacetType::Tri3 => &TRI_EDGES,
            FacetType::Quad4 => &QUAD_EDGES,
        };
        for &(i, j) in le {
            let eid = mesh.edge_index(fv[i], fv[j]).expect("facet edge");
            fixed_edges[eid] = true;
        }
    }
    build_ams_p1_from_masks(mesh, vkey, &fixed_edges, &fixed_verts)
}

/// Builds [`AmsP1Data`] from explicit boolean masks of constrained edges and constrained vertices.
pub fn build_ams_p1_from_masks(
    mesh: &MixedDiscretisation,
    vkey: &[u32],
    fixed_edges: &[bool],
    fixed_verts: &[bool],
) -> AmsP1Data {
    let free_edges: Vec<usize> = (0..mesh.n_edges()).filter(|&e| !fixed_edges[e]).collect();
    let free_vertices: Vec<usize> = (0..mesh.n_points()).filter(|&v| !fixed_verts[v]).collect();
    let g_full = gradient_p1_full(mesh, vkey);
    let gradient_p1 = g_full.submatrix(&free_edges, &free_vertices);
    let ev_full = edge_vectors_full(mesh, vkey);
    let edge_vectors = free_edges.iter().map(|&e| ev_full[e]).collect();
    AmsP1Data {
        free_edges,
        free_vertices,
        gradient_p1,
        edge_vectors,
    }
}

/// Constrained second-order (`N_2`) auxiliary-space topology ready for [`crate::linalg::ams::AmsP2::new`].
#[derive(Clone, Debug)]
pub struct AmsP2Data {
    /// Free `N_2` DOF indices in `0..ned2.n_dofs`.
    pub free_dofs_p2: Vec<usize>,
    /// Free `S_2` scalar DOF indices in `0..scalar_p2_tensor_count(mesh)`.
    pub free_scalar_p2: Vec<usize>,
    /// Algebraic operators ready to pass directly into [`crate::linalg::ams::AmsP2::new`].
    pub topology: AmsP2Topology,
}

/// Builds [`AmsP2Data`] after eliminating constrained PEC facets (`pec_facets` are indices in `0..mesh.n_facets()`).
pub fn build_ams_p2(mesh: &MixedDiscretisation, ned2: &NedMixed, pec_facets: &[usize]) -> AmsP2Data {
    assert_eq!(ned2.order, 2);
    let mut fixed_p2 = vec![false; ned2.n_dofs];
    let mut fixed_edges = vec![false; mesh.n_edges()];
    let mut fixed_verts = vec![false; mesh.n_points()];
    let mut fixed_quad_faces = vec![false; mesh.n_quad_faces()];

    for &f in pec_facets {
        for d in ned2.facet_dofs(mesh, f) {
            fixed_p2[d] = true;
        }
        let fv = mesh.facet(f);
        for &v in fv {
            fixed_verts[v as usize] = true;
        }
        let le: &[(usize, usize)] = match mesh.facet_type(f) {
            FacetType::Edge2 => &[(0, 1)],
            FacetType::Tri3 => &TRI_EDGES,
            FacetType::Quad4 => &QUAD_EDGES,
        };
        for &(i, j) in le {
            fixed_edges[mesh.edge_index(fv[i], fv[j]).unwrap()] = true;
        }
        if mesh.facet_type(f) == FacetType::Quad4 {
            fixed_quad_faces[mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap()] = true;
        }
    }
    build_ams_p2_from_masks(mesh, ned2, &fixed_p2, &fixed_edges, &fixed_verts, &fixed_quad_faces)
}

/// Builds [`AmsP2Data`] from explicit boolean masks of constrained `N_2` DOFs, edges, vertices, and 3D quad faces.
pub fn build_ams_p2_from_masks(
    mesh: &MixedDiscretisation,
    ned2: &NedMixed,
    fixed_p2: &[bool],
    fixed_edges: &[bool],
    fixed_verts: &[bool],
    fixed_quad_faces: &[bool],
) -> AmsP2Data {
    assert_eq!(ned2.order, 2);
    let p1 = build_ams_p1_from_masks(mesh, &ned2.vkey, fixed_edges, fixed_verts);
    let free_dofs_p2: Vec<usize> = (0..ned2.n_dofs).filter(|&d| !fixed_p2[d]).collect();

    let n_scalar = scalar_p2_tensor_count(mesh);
    let nv = mesh.n_points();
    let ne = mesh.n_edges();
    let mut fixed_scalar = vec![false; n_scalar];
    for v in 0..nv {
        if fixed_verts[v] {
            fixed_scalar[v] = true;
        }
    }
    for e in 0..ne {
        if fixed_edges[e] {
            fixed_scalar[nv + e] = true;
        }
    }
    if mesh.dim == 3 {
        for qf in 0..mesh.n_quad_faces() {
            if fixed_quad_faces[qf] {
                fixed_scalar[nv + ne + qf] = true;
            }
        }
    }
    let free_scalar_p2: Vec<usize> = (0..n_scalar).filter(|&s| !fixed_scalar[s]).collect();

    let p12_full = p1_to_p2_full(mesh, ned2);
    let p1_to_p2 = p12_full.submatrix(&free_dofs_p2, &p1.free_edges);

    let g2_full = gradient_p2_full(mesh, ned2);
    let gradient_p2 = g2_full.submatrix(&free_dofs_p2, &free_scalar_p2);

    let mut global_to_free = vec![usize::MAX; ned2.n_dofs];
    for (fi, &gd) in free_dofs_p2.iter().enumerate() {
        global_to_free[gd] = fi;
    }
    let blocks_full = ams_p2_blocks_full(mesh, ned2);
    let mut blocks = Vec::with_capacity(blocks_full.len());
    for blk in blocks_full {
        let fb: Vec<usize> = blk
            .into_iter()
            .filter_map(|gd| {
                let fi = global_to_free[gd];
                (fi != usize::MAX).then_some(fi)
            })
            .collect();
        if !fb.is_empty() {
            blocks.push(fb);
        }
    }

    AmsP2Data {
        free_dofs_p2,
        free_scalar_p2,
        topology: AmsP2Topology {
            p1_to_p2,
            gradient_p1: p1.gradient_p1,
            gradient_p2: Some(gradient_p2),
            edge_vectors: p1.edge_vectors,
            blocks,
        },
    }
}

// ------------------------------------------------------------------ Mixed Maxwell curl-curl & mass assemblers

/// Assembles real 3D curl-curl `K_{ij} = \int (\nabla\times w_i) \cdot (\mu^{-1}(x_q) \nabla\times w_j)`
/// and mass `M_{ij} = \int w_i \cdot (\varepsilon(x_q) w_j)` matrices over a 3D [`MixedDiscretisation`].
pub fn assemble_maxwell_3d_real<FMu, FEps>(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    quad_degree: usize,
    mu_inv_at: FMu,
    eps_at: FEps,
) -> (Csr, Csr)
where
    FMu: Fn(usize, [f64; 3]) -> [[f64; 3]; 3] + Sync,
    FEps: Fn(usize, [f64; 3]) -> [[f64; 3]; 3] + Sync,
{
    assert_eq!(mesh.dim, 3);
    let nc = mesh.n_cells();
    let n = ned.n_dofs;
    let run = |which: usize| -> Csr {
        assemble_parallel(nc, n, n, |range, coo| {
            let mut scratch = NedEvalScratch::new();
            let mut w = Vec::new();
            let mut cu = Vec::new();
            let mut ke = Vec::new();
            for ci in range {
                let iso = IsoparametricCell::from_mesh(mesh, ci);
                let deg = if iso.curved { quad_degree.max(6) } else { quad_degree };
                let rule = quad_rule_mixed(iso.cell_type, deg);
                let dofs = ned.cell_dofs(ci);
                let nb = dofs.len();
                ke.clear();
                ke.resize(nb * nb, 0.0);
                for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                    let jac = eval_cell_basis_3d(mesh, ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                    let weight = wq * jac.det.abs();
                    let tens = if which == 0 {
                        mu_inv_at(ci, jac.x)
                    } else {
                        eps_at(ci, jac.x)
                    };
                    let vecs = if which == 0 { &cu } else { &w };
                    for i in 0..nb {
                        let vi = vecs[i];
                        let tvi = [
                            tens[0][0] * vi[0] + tens[0][1] * vi[1] + tens[0][2] * vi[2],
                            tens[1][0] * vi[0] + tens[1][1] * vi[1] + tens[1][2] * vi[2],
                            tens[2][0] * vi[0] + tens[2][1] * vi[1] + tens[2][2] * vi[2],
                        ];
                        for j in 0..nb {
                            let vj = vecs[j];
                            ke[i * nb + j] += weight * (tvi[0] * vj[0] + tvi[1] * vj[1] + tvi[2] * vj[2]);
                        }
                    }
                }
                coo.push_block(dofs, dofs, &ke);
            }
        })
    };
    (run(0), run(1))
}

/// Assembles complex symmetric 3D curl-curl `K` and mass `M` matrices over a 3D [`MixedDiscretisation`],
/// evaluating `mu_inv_at(ci, x_q)` and `eps_at(ci, x_q)` pointwise at every Gauss quadrature node.
pub fn assemble_maxwell_3d_complex<FMu, FEps>(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    quad_degree: usize,
    mu_inv_at: FMu,
    eps_at: FEps,
) -> (CCsr, CCsr)
where
    FMu: Fn(usize, [f64; 3]) -> [[C; 3]; 3] + Sync,
    FEps: Fn(usize, [f64; 3]) -> [[C; 3]; 3] + Sync,
{
    assert_eq!(mesh.dim, 3);
    let nc = mesh.n_cells();
    let n = ned.n_dofs;
    let nt = crate::par::threads_up_to(16);
    let nparts = if nc < 1000 { 1 } else { nc.div_ceil(1024).min(16) };
    let chunk = nc.div_ceil(nparts).max(1);

    let run = |which: usize| -> CCsr {
        let build_part = |t: usize| -> CCsr {
            let (c0, c1) = ((t * chunk).min(nc), ((t + 1) * chunk).min(nc));
            if c0 >= c1 {
                return CCsr::zeros(n, n);
            }
            let mut coo = CCoo::new(n, n);
            let mut scratch = NedEvalScratch::new();
            let mut w = Vec::new();
            let mut cu = Vec::new();
            let mut ke = Vec::new();
            for ci in c0..c1 {
                let iso = IsoparametricCell::from_mesh(mesh, ci);
                let deg = if iso.curved { quad_degree.max(6) } else { quad_degree };
                let rule = quad_rule_mixed(iso.cell_type, deg);
                let dofs = ned.cell_dofs(ci);
                let nb = dofs.len();
                ke.clear();
                ke.resize(nb * nb, c(0.0, 0.0));
                for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                    let jac = eval_cell_basis_3d(mesh, ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                    let weight = wq * jac.det.abs();
                    let tens = if which == 0 {
                        mu_inv_at(ci, jac.x)
                    } else {
                        eps_at(ci, jac.x)
                    };
                    let vecs = if which == 0 { &cu } else { &w };
                    for i in 0..nb {
                        let vi = vecs[i];
                        let tvi = [
                            tens[0][0] * vi[0] + tens[0][1] * vi[1] + tens[0][2] * vi[2],
                            tens[1][0] * vi[0] + tens[1][1] * vi[1] + tens[1][2] * vi[2],
                            tens[2][0] * vi[0] + tens[2][1] * vi[1] + tens[2][2] * vi[2],
                        ];
                        for j in 0..nb {
                            let vj = vecs[j];
                            ke[i * nb + j] += (tvi[0] * vj[0] + tvi[1] * vj[1] + tvi[2] * vj[2]) * weight;
                        }
                    }
                }
                coo.push_block(dofs, dofs, &ke);
            }
            coo.to_csr()
        };
        let parts: Vec<CCsr> = if nt <= 1 || nparts <= 1 {
            (0..nparts).map(build_part).collect()
        } else {
            (0..nparts).into_par_iter().map(build_part).collect()
        };
        CCsr::sum_all(parts)
    };
    (run(0), run(1))
}

/// Assembles real 2D curl-curl `K_{ij} = \int \mu^{-1}(x_q) (\mathrm{curl}\,w_i)(\mathrm{curl}\,w_j)`
/// and mass `M_{ij} = \int w_i \cdot (\varepsilon(x_q) w_j)` matrices over a 2D [`MixedDiscretisation`].
pub fn assemble_maxwell_2d_real<FMu, FEps>(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    quad_degree: usize,
    mu_inv_at: FMu,
    eps_at: FEps,
) -> (Csr, Csr)
where
    FMu: Fn(usize, [f64; 2]) -> f64 + Sync,
    FEps: Fn(usize, [f64; 2]) -> [[f64; 2]; 2] + Sync,
{
    assert_eq!(mesh.dim, 2);
    let nc = mesh.n_cells();
    let n = ned.n_dofs;
    let run = |which: usize| -> Csr {
        assemble_parallel(nc, n, n, |range, coo| {
            let mut scratch = NedEvalScratch::new();
            let mut w = Vec::new();
            let mut cu = Vec::new();
            let mut ke = Vec::new();
            for ci in range {
                let iso = IsoparametricCell::from_mesh(mesh, ci);
                let deg = if iso.curved { quad_degree.max(6) } else { quad_degree };
                let rule = quad_rule_mixed(iso.cell_type, deg);
                let dofs = ned.cell_dofs(ci);
                let nb = dofs.len();
                ke.clear();
                ke.resize(nb * nb, 0.0);
                for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                    let jac = eval_cell_basis_2d(mesh, ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                    let weight = wq * jac.det.abs();
                    let x2 = [jac.x[0], jac.x[1]];
                    if which == 0 {
                        let mui = mu_inv_at(ci, x2) * weight;
                        for i in 0..nb {
                            let ci_w = mui * cu[i];
                            for j in 0..nb {
                                ke[i * nb + j] += ci_w * cu[j];
                            }
                        }
                    } else {
                        let eps = eps_at(ci, x2);
                        for i in 0..nb {
                            let vi = w[i];
                            let evi = [
                                eps[0][0] * vi[0] + eps[0][1] * vi[1],
                                eps[1][0] * vi[0] + eps[1][1] * vi[1],
                            ];
                            for j in 0..nb {
                                let vj = w[j];
                                ke[i * nb + j] += weight * (evi[0] * vj[0] + evi[1] * vj[1]);
                            }
                        }
                    }
                }
                coo.push_block(dofs, dofs, &ke);
            }
        })
    };
    (run(0), run(1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapt::work::simplex_to_adapt;
    use crate::linalg::ams::{Ams, AmsOptions, AmsP2, AmsP2Options};
    use crate::linalg::ksp::{Cg, CgKind, CgOptions, Tolerances};
    use crate::meshing::pml_box::{build_pml, PmlSpec};

    fn csr_inf_norm(m: &Csr) -> f64 {
        (0..m.nrows())
            .map(|r| m.row(r).1.iter().map(|v| v.abs()).sum::<f64>())
            .fold(0.0, f64::max)
    }

    /// Build a rich 3D mixed mesh (`Tet` physical core + `Prism` face layers + `Hex` edge & corner layers)
    /// via `pml_box::build_pml` with non-trivial anisotropic distortion AND curved interior + PML edges!
    fn build_skewed_3d_tet_prism_hex_mesh(curved: bool) -> MixedDiscretisation {
        let block = crate::meshing::geomdrive::block(1.0, 0.9, 0.8, 0.45, 3).unwrap().mesh;
        let phys = simplex_to_adapt::<3>(&block).unwrap();
        let spec = PmlSpec {
            outer_label: Some(1),
            grading: 1.8,
            ..PmlSpec::<3>::uniform(0.35, 2)
        };
        let mut cm = build_pml(&phys, &spec).unwrap().mesh;
        assert!(cm.types.contains(&CellType::Tet));
        assert!(cm.types.contains(&CellType::Prism));
        assert!(cm.types.contains(&CellType::Hex));

        // Apply a smooth non-affine warp to every point so all Prisms and Hexes have non-constant, non-orthogonal Jacobians!
        for p in &mut cm.points {
            let (x, y, z) = (p[0], p[1], p[2]);
            p[0] = x + 0.04 * y - 0.03 * z + 0.015 * x * y;
            p[1] = y - 0.03 * x + 0.05 * z + 0.012 * y * z;
            p[2] = z + 0.025 * x - 0.02 * y + 0.01 * x * z;
        }
        cm.check().unwrap();

        if curved {
            // Curve several interior & interface edges and a shared Prism-Hex quad-face center!
            let md_tmp = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
            let mut en = Vec::new();
            for (eid, &(a, b)) in md_tmp.edges.iter().enumerate() {
                if eid % 11 == 0 {
                    let (pa, pb) = (cm.points[a as usize], cm.points[b as usize]);
                    let mid = [
                        0.5 * (pa[0] + pb[0]) + 0.004 * ((eid as f64).sin()),
                        0.5 * (pa[1] + pb[1]) + 0.004 * ((eid as f64).cos()),
                        0.5 * (pa[2] + pb[2]) - 0.003 * ((eid as f64 * 0.7).sin()),
                    ];
                    en.push(([a, b], mid));
                }
            }
            cm.set_edge_nodes(en);
            let mut qfn = Vec::new();
            for (qfid, &q) in md_tmp.quad_faces.iter().enumerate() {
                if qfid % 7 == 0 {
                    let c4 = q.map(|v| cm.points[v as usize]);
                    let avg = [
                        0.25 * (c4[0][0] + c4[1][0] + c4[2][0] + c4[3][0]) + 0.003,
                        0.25 * (c4[0][1] + c4[1][1] + c4[2][1] + c4[3][1]) - 0.002,
                        0.25 * (c4[0][2] + c4[1][2] + c4[2][2] + c4[3][2]) + 0.002,
                    ];
                    qfn.push((q, avg));
                }
            }
            cm.set_quad_face_nodes(qfn);
        }

        let md = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        md.validate_conforming().unwrap();
        md
    }

    #[test]
    fn de_rham_exactness_3d_tet_prism_hex_straight_and_curved() {
        // Anisotropic full 3x3 mu^{-1} and eps tensors varying across domains to stress all cross-derivatives!
        let mu_inv = |ci: usize, x: [f64; 3]| -> [[f64; 3]; 3] {
            let s = 1.0 + 0.2 * ((ci % 5) as f64) + 0.1 * x[0];
            [
                [1.1 * s, 0.15 * s, -0.08 * s],
                [0.15 * s, 0.95 * s, 0.12 * s],
                [-0.08 * s, 0.12 * s, 1.25 * s],
            ]
        };
        let eps = |_ci: usize, _x: [f64; 3]| -> [[f64; 3]; 3] {
            [[2.1, 0.1, 0.0], [0.1, 2.4, -0.05], [0.0, -0.05, 1.9]]
        };

        let id3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for curved in [false, true] {
            let md = build_skewed_3d_tet_prism_hex_mesh(curved);

            // 1. Unit curl-curl operator (\mu^{-1} = I): verify unnormalized ||K_1 G_1||_\infty < 1e-13 and ||K_2 G_2||_\infty < 1e-13!
            let ned1 = NedMixed::new(&md, 1);
            let (k1_id, _) = assemble_maxwell_3d_real(&md, &ned1, 4, |_, _| id3, |_, _| id3);
            let g1 = gradient_p1_full(&md, &ned1.vkey);
            let raw_id_kg1 = csr_inf_norm(&k1_id.matmul(&g1));
            assert!(
                raw_id_kg1 < 1e-13,
                "3D curved={curved} unit ||K1 G1||_inf = {raw_id_kg1:e}"
            );

            let ned2 = NedMixed::new(&md, 2);
            let (k2_id, _) = assemble_maxwell_3d_real(&md, &ned2, 6, |_, _| id3, |_, _| id3);
            let g2 = gradient_p2_full(&md, &ned2);
            let raw_id_kg2 = csr_inf_norm(&k2_id.matmul(&g2));
            assert!(
                raw_id_kg2 < 1e-13,
                "3D curved={curved} unit ||K2 G2||_inf = {raw_id_kg2:e}"
            );

            // 2. Spatially varying full 3x3 anisotropic \mu^{-1}(x_q): verify relative ||K G||_\infty / ||K||_\infty < 1e-14 (~1e-16)!
            let (k1, m1) = assemble_maxwell_3d_real(&md, &ned1, 4, mu_inv, eps);
            let k1_norm = csr_inf_norm(&k1).max(1.0);
            let kg1 = k1.matmul(&g1);
            let raw_kg1 = csr_inf_norm(&kg1);
            let err_g1 = raw_kg1 / k1_norm;
            assert!(
                err_g1 < 1e-14,
                "3D curved={curved} Order 1 de Rham violation: ||K1 G1||_inf = {raw_kg1:e} (rel {err_g1:e})"
            );

            let (k2, m2) = assemble_maxwell_3d_real(&md, &ned2, 6, mu_inv, eps);
            let k2_norm = csr_inf_norm(&k2).max(1.0);
            let kg2 = k2.matmul(&g2);
            let raw_kg2 = csr_inf_norm(&kg2);
            let err_g2 = raw_kg2 / k2_norm;
            assert!(
                err_g2 < 1e-14,
                "3D curved={curved} Order 2 de Rham violation: ||K2 G2||_inf = {raw_kg2:e} (rel {err_g2:e})"
            );

            // Verify exact hierarchical subspace identity: P_{1->2}^T K_2 P_{1->2} == K_1 and P_{1->2}^T M_2 P_{1->2} == M_1!
            let p12 = p1_to_p2_full(&md, &ned2);
            let p21 = p12.transpose();
            let k1_galerkin = p21.matmul(&k2.matmul(&p12));
            let m1_galerkin = p21.matmul(&m2.matmul(&p12));
            let diff_k = csr_inf_norm(&k1.add_scaled(&k1_galerkin, -1.0));
            let diff_m = csr_inf_norm(&m1.add_scaled(&m1_galerkin, -1.0));
            let (k1_d6, m1_d6) = assemble_maxwell_3d_real(&md, &ned1, 6, mu_inv, eps);
            let diff_k6 = csr_inf_norm(&k1_d6.add_scaled(&k1_galerkin, -1.0));
            let diff_m6 = csr_inf_norm(&m1_d6.add_scaled(&m1_galerkin, -1.0));
            let _ = (diff_k, diff_m);
            println!(
                "[de_rham_3d curved={curved}] cells={} n_p1={} n_p2={} unit(||K1 G1||={raw_id_kg1:.3e}, ||K2 G2||={raw_id_kg2:.3e}) aniso_rel(g1={err_g1:.3e}, g2={err_g2:.3e}) dk6={diff_k6:.3e}",
                md.n_cells(),
                ned1.n_dofs,
                ned2.n_dofs
            );
            assert!(
                diff_k6 < 1e-13 && diff_m6 < 1e-13,
                "P1->P2 Galerkin mismatch curved={curved}: dk={diff_k6:e}, dm={diff_m6:e}"
            );
        }
    }

    #[test]
    fn de_rham_exactness_2d_tri_quad_straight_and_curved() {
        let plate = crate::meshing::geomdrive::plate(1.0, 0.8, 0.25).unwrap().mesh;
        let phys = simplex_to_adapt::<2>(&plate).unwrap();
        let spec = PmlSpec {
            outer_label: Some(1),
            grading: 1.5,
            ..PmlSpec::<2>::uniform(0.3, 2)
        };
        let mut cm = build_pml(&phys, &spec).unwrap().mesh;
        assert!(cm.types.contains(&CellType::Tri) && cm.types.contains(&CellType::Quad));

        for p in &mut cm.points {
            let (x, z) = (p[0], p[1]);
            p[0] = x + 0.05 * z + 0.02 * x * z;
            p[1] = z - 0.04 * x + 0.015 * x * x;
        }
        cm.check().unwrap();
        let md_s = MixedDiscretisation::from_cell_mesh(&cm).unwrap();
        md_s.validate_conforming().unwrap();

        // Add curved edges too
        let mut en = Vec::new();
        for (eid, &(a, b)) in md_s.edges.iter().enumerate() {
            if eid % 5 == 0 {
                let (pa, pb) = (cm.points[a as usize], cm.points[b as usize]);
                en.push(([a, b], [0.5 * (pa[0] + pb[0]) + 0.003, 0.5 * (pa[1] + pb[1]) - 0.002]));
            }
        }
        cm.set_edge_nodes(en);
        let md_c = MixedDiscretisation::from_cell_mesh(&cm).unwrap();

        for md in [md_s, md_c] {
            let ned1 = NedMixed::new(&md, 1);
            let (k1, _) = assemble_maxwell_2d_real(&md, &ned1, 4, |_, _| 1.3, |_, _| [[2.0, 0.2], [0.2, 1.8]]);
            let g1 = gradient_p1_full(&md, &ned1.vkey);
            let raw1 = csr_inf_norm(&k1.matmul(&g1));
            assert!(raw1 < 1e-13, "2D ||K1 G1||_inf = {raw1:e}");

            let ned2 = NedMixed::new(&md, 2);
            let (k2, _) = assemble_maxwell_2d_real(&md, &ned2, 6, |_, _| 1.3, |_, _| [[2.0, 0.2], [0.2, 1.8]]);
            let g2 = gradient_p2_full(&md, &ned2);
            let raw2 = csr_inf_norm(&k2.matmul(&g2));
            assert!(raw2 < 1e-13, "2D ||K2 G2||_inf = {raw2:e}");
        }
    }

    #[test]
    fn ams_and_amsp2_converge_rapidly_on_3d_mixed_tet_prism_hex() {
        let md = build_skewed_3d_tet_prism_hex_mesh(false);
        let pec_facets: Vec<usize> = (0..md.n_facets()).filter(|&f| md.exterior_facets[f]).collect();
        assert!(!pec_facets.is_empty());

        let id3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        // 1. Order 1 SPD curl-curl + beta * M solved with Ams + PCG
        let ned1 = NedMixed::new(&md, 1);
        let (k1, m1) = assemble_maxwell_3d_real(&md, &ned1, 4, |_, _| id3, |_, _| id3);
        let a1_full = k1.add_scaled(&m1, 1.0);
        let ams1_in = build_ams_p1(&md, &ned1.vkey, &pec_facets);
        let a1 = a1_full.submatrix(&ams1_in.free_edges, &ams1_in.free_edges);
        let k1_free = k1.submatrix(&ams1_in.free_edges, &ams1_in.free_edges);
        let kg1_free = csr_inf_norm(&k1_free.matmul(&ams1_in.gradient_p1));
        assert!(kg1_free < 1e-13, "constrained ||K1_free G1_free||_inf = {kg1_free:e}");

        let ams1 = Ams::new(&a1, &ams1_in.gradient_p1, &ams1_in.edge_vectors, &AmsOptions::default()).unwrap();
        let rhs1: Vec<f64> = (0..a1.nrows()).map(|i| ((i as f64 + 1.0) * 0.37).sin()).collect();
        let mut x1 = vec![0.0; a1.nrows()];
        let tol = Tolerances {
            rtol: 1e-10,
            max_it: 100,
            ..Default::default()
        };
        let st1 = Cg::new(CgOptions {
            tol,
            kind: CgKind::Hermitian,
            ..Default::default()
        })
        .solve(&a1, &ams1, &rhs1, &mut x1);
        assert!(
            st1.reason.is_converged() && st1.iterations <= 30,
            "Ams P1 failed or slow on mixed mesh: iters={}, res={:e}, reason={:?}",
            st1.iterations,
            st1.residual_norm,
            st1.reason
        );

        // 2. Order 2 SPD & complex-symmetric (PML-like) curl-curl + M solved with AmsP2 + PCG / COCG!
        let ned2 = NedMixed::new(&md, 2);
        let (k2, m2) = assemble_maxwell_3d_real(&md, &ned2, 6, |_, _| id3, |_, _| id3);
        let ams2_in = build_ams_p2(&md, &ned2, &pec_facets);
        let k2_free = k2.submatrix(&ams2_in.free_dofs_p2, &ams2_in.free_dofs_p2);
        let g2_free = ams2_in.topology.gradient_p2.as_ref().unwrap();
        let kg2_free = csr_inf_norm(&k2_free.matmul(g2_free));
        assert!(kg2_free < 1e-13, "constrained ||K2_free G2_free||_inf = {kg2_free:e}");

        let a2_full = k2.add_scaled(&m2, 1.0);
        let a2 = a2_full.submatrix(&ams2_in.free_dofs_p2, &ams2_in.free_dofs_p2);
        let ams2 = AmsP2::new(&a2, &ams2_in.topology, &AmsP2Options::default()).unwrap();
        let rhs2: Vec<f64> = (0..a2.nrows()).map(|i| ((i as f64 + 1.0) * 0.23).cos()).collect();
        let mut x2 = vec![0.0; a2.nrows()];
        let st2 = Cg::new(CgOptions {
            tol,
            kind: CgKind::Hermitian,
            ..Default::default()
        })
        .solve(&a2, &ams2, &rhs2, &mut x2);
        assert!(
            st2.reason.is_converged() && st2.iterations <= 35,
            "AmsP2 real PCG failed or slow on mixed mesh: iters={}, res={:e}, reason={:?}",
            st2.iterations,
            st2.residual_norm,
            st2.reason
        );

        // Complex-symmetric damped Maxwell operator: K_2 + (1 - 0.5 i) M_2 with AmsP2<C> + COCG!
        let a2_c = k2_free.to_complex().add_scaled(
            &m2.submatrix(&ams2_in.free_dofs_p2, &ams2_in.free_dofs_p2).to_complex(),
            c(1.0, -0.5),
        );
        let ams2_c = AmsP2::new(&a2_c, &ams2_in.topology, &AmsP2Options::default()).unwrap();
        let rhs2_c: Vec<C> = rhs2.iter().map(|&v| c(v, 0.25 * v)).collect();
        let mut x2_c = vec![c(0.0, 0.0); a2_c.nrows()];
        let st2_c = Cg::new(CgOptions {
            tol,
            kind: CgKind::Symmetric,
            ..Default::default()
        })
        .solve(&a2_c, &ams2_c, &rhs2_c, &mut x2_c);
        println!(
            "[ams_mixed_3d] Ams(P1, {} free): {} iters; AmsP2(real, {} free): {} iters; AmsP2(complex PML, {} free): {} iters; ||K2_free G2_free||_inf={kg2_free:.3e}",
            a1.nrows(),
            st1.iterations,
            a2.nrows(),
            st2.iterations,
            a2_c.nrows(),
            st2_c.iterations
        );
        assert!(
            st2_c.reason.is_converged() && st2_c.iterations <= 40,
            "AmsP2 complex COCG failed or slow on mixed mesh: iters={}, res={:e}, reason={:?}",
            st2_c.iterations,
            st2_c.residual_norm,
            st2_c.reason
        );
    }
}
