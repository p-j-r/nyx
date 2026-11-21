// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Kinetic-inductance surface-impedance boundary conditions (`surface_impedance`, KICQ IBC)
//! and single-facet deduplication for 3D `cavity` and 2D `modes` / `axisymmetric_modes`.

use crate::cavity::{basis3, Ned3};
use crate::modes::C0;
use crate::problem::SurfaceImpedanceBc;
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{facet_rule, geometry, quad_rule, shape, CurvedTet};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use nyx_physics::Discretisation;
use std::collections::{BTreeSet, HashMap};
use std::f64::consts::PI;

/// Vacuum magnetic permeability `mu_0` in H/m.
pub const MU0: f64 = 4.0 * PI * 1e-7;

const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Series-to-parallel conversion of `Z_s(omega_ref) = R_s + i omega_ref L_k_sq`:
/// returns `(1 / L_eff, G_eff)` where `L_eff = L_k_sq + R_s^2 / (omega_ref^2 L_k_sq)` and
/// `G_eff = R_s / (R_s^2 + omega_ref^2 L_k_sq^2)`.
#[inline]
pub fn ibc_parallel_params(r_s: f64, l_k_sq: f64, omega_ref: f64) -> (f64, f64) {
    if l_k_sq > 0.0 {
        let wl = omega_ref * l_k_sq;
        let denom = r_s * r_s + wl * wl;
        let inv_l_eff = if denom > 0.0 { (omega_ref * wl) / denom } else { 1.0 / l_k_sq };
        let g_eff = if denom > 0.0 { r_s / denom } else { 0.0 };
        (inv_l_eff, g_eff)
    } else if r_s > 0.0 {
        (0.0, 1.0 / r_s)
    } else {
        (0.0, 0.0)
    }
}

#[inline]
pub fn sorted3(a: u32, b: u32, c: u32) -> [u32; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn pt3(mesh: &Discretisation, v: u32) -> [f64; 3] {
    let p = mesh.point(v as usize);
    [p[0], p[1], p[2]]
}

/// 3D assembly of `surface_impedance` boundary conditions.
pub struct IbcAssembly3D {
    /// Unique sorted vertex triples of `surface_impedance` facets.
    pub ibc_tris: Vec<[u32; 3]>,
    /// Full-DOF complex IBC contribution `K_IBC + i k_{0,ref} C_IBC` (`n x n`).
    pub kc_ibc_full: CCsr,
}

/// Builds a map from each sorted triangle face `[v0, v1, v2]` in a 3D tetrahedral mesh to one adjacent
/// cell index `ci` and its 3 local vertex indices `[la, lb, lc]` in `0..4`.
pub fn build_face_to_cell_3d(mesh: &Discretisation) -> HashMap<[u32; 3], (usize, [usize; 3])> {
    let mut map = HashMap::with_capacity(mesh.n_cells() * 2);
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for f in FACES {
            let key = sorted3(v[f[0]], v[f[1]], v[f[2]]);
            map.entry(key).or_insert((ci, f));
        }
    }
    map
}

fn eval_triangle_tangential_basis(
    mesh: &Discretisation,
    ned: &Ned3,
    geo: &femcore::fem::Geometry,
    face_to_cell: &HashMap<[u32; 3], (usize, [usize; 3])>,
    tri_sorted: [u32; 3],
) -> Result<(f64, [f64; 3], Vec<usize>, Vec<f64>, Vec<Vec<[f64; 3]>>), String> {
    let &(ci, lf) = face_to_cell
        .get(&tri_sorted)
        .ok_or_else(|| format!("facet {:?} is not a face of any tetrahedron in the mesh", tri_sorted))?;
    let v = mesh.cell(ci);
    let (p0, p1, p2) = (pt3(mesh, v[lf[0]]), pt3(mesh, v[lf[1]]), pt3(mesh, v[lf[2]]));
    let n_raw = cross3(sub3(p1, p0), sub3(p2, p0));
    let n_len = norm3(n_raw);
    let area = 0.5 * n_len;
    if area <= 0.0 {
        return Err("degenerate zero-area boundary triangle".into());
    }
    let n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];

    let bnd_dofs = ned.boundary_dofs(mesh, tri_sorted);
    let cell_dofs = &ned.cell_dofs[ci];
    let local_slots: Vec<usize> = bnd_dofs
        .iter()
        .map(|&gd| {
            cell_dofs
                .iter()
                .position(|&cd| cd == gd)
                .expect("every boundary dof of a face belongs to its adjacent tetrahedron")
        })
        .collect();

    let gid = ned.keys(mesh, ci);
    let mut gs = [[0.0; 3]; 4];
    for b in 0..4 {
        for k in 0..3 {
            gs[b][k] = geo.g[ci * 12 + b * 3 + k];
        }
    }
    let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
    let quad2 = quad_rule(2, 4);
    let nq = quad2.w.len();
    let nb_face = bnd_dofs.len();
    let mut weights = Vec::with_capacity(nq);
    let mut psi_q = Vec::with_capacity(nq);
    let (mut w, mut cu) = (Vec::new(), Vec::new());

    for q in 0..nq {
        let l2 = &quad2.lam[q];
        let mut lam4 = [0.0; 4];
        lam4[lf[0]] = l2[0];
        lam4[lf[1]] = l2[1];
        lam4[lf[2]] = l2[2];
        match &curved {
            Some(ct) => basis3(ned.order, &lam4, &ct.jacobian(&lam4).0, &gid, &mut w, &mut cu),
            None => basis3(ned.order, &lam4, &gs, &gid, &mut w, &mut cu),
        }
        let mut row = Vec::with_capacity(nb_face);
        for &s in &local_slots {
            let ws = w[s];
            let wn = dot3(ws, n_hat);
            row.push([ws[0] - wn * n_hat[0], ws[1] - wn * n_hat[1], ws[2] - wn * n_hat[2]]);
        }
        weights.push(quad2.w[q] * area);
        psi_q.push(row);
    }

    Ok((area, n_hat, bnd_dofs, weights, psi_q))
}

/// Deduplicates and assembles 3D `surface_impedance` boundary conditions.
pub fn assemble_ibc_3d(
    mesh: &Discretisation,
    ned: &Ned3,
    ibc_bcs: &[SurfaceImpedanceBc],
    f0: f64,
) -> Result<IbcAssembly3D, String> {
    let n = ned.n_dofs();
    let geo = geometry(mesh);
    let face_to_cell = build_face_to_cell_3d(mesh);
    let mut seen_ibc = BTreeSet::new();
    let mut ibc_tris = Vec::new();
    let mut kc_ibc_coo = CCoo::new(n, n);

    for bc in ibc_bcs {
        let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { f0 };
        let omega_ref = 2.0 * PI * f_ref;
        let k0_ref = omega_ref / C0;
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
        let alpha_k = MU0 * inv_l_eff;
        let alpha_c = MU0 * C0 * g_eff;

        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = sorted3(fv[0], fv[1], fv[2]);
            if !seen_ibc.insert(key) {
                continue;
            }
            ibc_tris.push(key);
            let (_, _, bnd_dofs, weights, psi_q) = eval_triangle_tangential_basis(mesh, ned, &geo, &face_to_cell, key)?;

            let nb = bnd_dofs.len();
            let mut loc_g = vec![0.0; nb * nb];
            for (q, &wq) in weights.iter().enumerate() {
                let pq = &psi_q[q];
                for a in 0..nb {
                    for b in 0..nb {
                        loc_g[a * nb + b] += wq * dot3(pq[a], pq[b]);
                    }
                }
            }
            if alpha_k > 0.0 || alpha_c > 0.0 {
                let kce: Vec<C> = loc_g.iter().map(|&g| c(alpha_k * g, k0_ref * alpha_c * g)).collect();
                kc_ibc_coo.push_block(&bnd_dofs, &bnd_dofs, &kce);
            }
        }
    }

    Ok(IbcAssembly3D { ibc_tris, kc_ibc_full: kc_ibc_coo.to_csr() })
}

fn build_edge_to_cell_2d(mesh: &Discretisation) -> HashMap<(u32, u32), (usize, [usize; 2])> {
    let mut map = HashMap::with_capacity(mesh.n_edges());
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for [la, lb] in [[0, 1], [1, 2], [0, 2]] {
            let key = (v[la].min(v[lb]), v[la].max(v[lb]));
            map.entry(key).or_insert((ci, [la, lb]));
        }
    }
    map
}

/// Collects PEC boundary edge indices `be` and vertex indices `bn` for a 2D mesh,
/// excluding `surface_impedance` and `pmc` facets from exterior default PEC walls.
pub fn pec_edges_and_nodes_2d(
    mesh: &Discretisation,
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
) -> (Vec<usize>, Vec<usize>) {
    let exterior = mesh.exterior_facets();
    let ibc_labels: Vec<i32> = ibc_bcs.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let is_pec_facet = |f: usize| -> bool {
        let l = mesh.facet_labels[f];
        exterior[f] && !ibc_labels.contains(&l) && !pmc_labels.contains(&l)
    };
    let mut be: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
        .collect();
    be.sort_unstable();
    be.dedup();
    let mut bn: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
        .collect();
    bn.sort_unstable();
    bn.dedup();
    (be, bn)
}

/// Assembles 2D `surface_impedance` contributions into `(a_coo, b_coo)` for `modes_complex`.
pub fn add_ibc_2d_waveguide(
    mesh: &Discretisation,
    order: usize,
    freq: f64,
    ibc_bcs: &[SurfaceImpedanceBc],
    a_coo: &mut CCoo,
    b_coo: &mut CCoo,
) -> Result<(), String> {
    if ibc_bcs.is_empty() {
        return Ok(());
    }
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let rule = facet_rule(1);
    let edge_to_cell = build_edge_to_cell_2d(mesh);

    let mut seen_ibc = BTreeSet::new();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();

    for bc in ibc_bcs {
        let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { freq };
        let omega_ref = 2.0 * PI * f_ref;
        let k0_ref = omega_ref / C0;
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
        let alpha = c(MU0 * inv_l_eff, k0_ref * MU0 * C0 * g_eff);
        if alpha.norm() == 0.0 {
            continue;
        }

        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = (fv[0].min(fv[1]), fv[0].max(fv[1]));
            if !seen_ibc.insert(key) {
                continue;
            }
            let &(ci, [la, lb]) = edge_to_cell
                .get(&key)
                .ok_or_else(|| format!("1D facet {:?} is not an edge of any 2D triangle", key))?;
            let (p0, p1) = (mesh.point(fv[0] as usize), mesh.point(fv[1] as usize));
            let dx = p1[0] - p0[0];
            let dy = p1[1] - p0[1];
            let len = dx.hypot(dy);
            if len <= 0.0 {
                continue;
            }
            let tau = [dx / len, dy / len];
            let g = tri_g(&geo, ci);
            mesh.cell_dofs(ci, order, &mut pd);
            let rows_e = &cd.gdofs[ci];
            let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();

            let mut ke_e = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut ke_p = vec![c(0.0, 0.0); nl * nl];
            for (q, l1) in rule.lam.iter().enumerate() {
                let mut l2 = [0.0; 3];
                l2[la] = l1[0];
                l2[lb] = l1[1];
                edge_basis(&l2, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv, _) = shape(&l2, order);
                let wq = rule.w[q] * len;
                for s in 0..nloc_e {
                    let ws_t = w[s][0] * tau[0] + w[s][1] * tau[1];
                    for t in 0..nloc_e {
                        let wt_t = w[t][0] * tau[0] + w[t][1] * tau[1];
                        ke_e[s * nloc_e + t] += alpha * (wq * ws_t * wt_t);
                    }
                }
                for a in 0..nl {
                    for b in 0..nl {
                        ke_p[a * nl + b] += alpha * (wq * nv[a] * nv[b]);
                    }
                }
            }
            a_coo.push_block(rows_e, rows_e, &ke_e);
            b_coo.push_block(&cols_p, &cols_p, &ke_p);
        }
    }
    Ok(())
}

/// Assembles 2D `surface_impedance` contributions into `a_coo` for `axisymmetric_modes` (`(r, z)` half-plane).
pub fn add_ibc_2d_axisym(
    mesh: &Discretisation,
    order: usize,
    f_search: f64,
    ibc_bcs: &[SurfaceImpedanceBc],
    a_coo: &mut CCoo,
) -> Result<(), String> {
    if ibc_bcs.is_empty() {
        return Ok(());
    }
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let rule = facet_rule(1);
    let edge_to_cell = build_edge_to_cell_2d(mesh);

    let mut seen_ibc = BTreeSet::new();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();

    for bc in ibc_bcs {
        let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { f_search };
        let omega_ref = 2.0 * PI * f_ref;
        let k0_ref = omega_ref / C0;
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
        let alpha = c(MU0 * inv_l_eff, k0_ref * MU0 * C0 * g_eff);
        if alpha.norm() == 0.0 {
            continue;
        }

        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = (fv[0].min(fv[1]), fv[0].max(fv[1]));
            if !seen_ibc.insert(key) {
                continue;
            }
            let &(ci, [la, lb]) = edge_to_cell
                .get(&key)
                .ok_or_else(|| format!("1D facet {:?} is not an edge of any 2D triangle", key))?;
            let (p0, p1) = (mesh.point(fv[0] as usize), mesh.point(fv[1] as usize));
            let dx = p1[0] - p0[0];
            let dy = p1[1] - p0[1];
            let len = dx.hypot(dy);
            if len <= 0.0 {
                continue;
            }
            let tau = [dx / len, dy / len];
            let g = tri_g(&geo, ci);
            mesh.cell_dofs(ci, order, &mut pd);
            let rows_e = &cd.gdofs[ci];
            let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
            let v_cell = mesh.cell(ci);

            let mut ke_e = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut ke_p = vec![c(0.0, 0.0); nl * nl];
            for (q, l1) in rule.lam.iter().enumerate() {
                let mut l2 = [0.0; 3];
                l2[la] = l1[0];
                l2[lb] = l1[1];
                let r: f64 = (0..3).map(|k| l2[k] * mesh.point(v_cell[k] as usize)[0]).sum();
                edge_basis(&l2, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv, _) = shape(&l2, order);
                let wr = rule.w[q] * len * r;
                for s in 0..nloc_e {
                    let ws_t = w[s][0] * tau[0] + w[s][1] * tau[1];
                    for t in 0..nloc_e {
                        let wt_t = w[t][0] * tau[0] + w[t][1] * tau[1];
                        ke_e[s * nloc_e + t] += alpha * (wr * ws_t * wt_t);
                    }
                }
                for a in 0..nl {
                    for b in 0..nl {
                        ke_p[a * nl + b] += alpha * (wr * nv[a] * nv[b]);
                    }
                }
            }
            a_coo.push_block(rows_e, rows_e, &ke_e);
            a_coo.push_block(&cols_p, &cols_p, &ke_p);
        }
    }
    Ok(())
}
