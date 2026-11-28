// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Volumetric electromagnetic loss and boundary surface-impedance dissipation exported as the
//! conservative nodal `"Q_opt"` (`"W/m^3"`, role `"heat_source"`) field across 3D `driven`,
//! 3D `cavity`, and 2D `modes` (`modes_complex`) studies (ADR 139, FR-38).
//!
//! Under Nyx's $\exp(+i\omega t)$ time-harmonic convention with passive relative permittivity
//! $\tilde{\boldsymbol{\varepsilon}}_r = \boldsymbol{\varepsilon}_{r,\text{re}} - i\boldsymbol{\varepsilon}_{r,\text{im}}$
//! ($\boldsymbol{\varepsilon}_{r,\text{im}} \succeq 0$), the cycle-averaged volumetric electromagnetic
//! heating density in non-PML physical domains is:
//! $$Q_{\text{vol}}(\mathbf{x}) = \frac{1}{2}\omega\varepsilon_0 \operatorname{Im}\!\big(-\mathbf{E}(\mathbf{x})^H \tilde{\boldsymbol{\varepsilon}}_r(\mathbf{x}) \mathbf{E}(\mathbf{x})\big) = \frac{1}{2}\omega\varepsilon_0 \sum_{a,b=0}^{2} \varepsilon_{r,\text{im},ab}\,\operatorname{Re}\!\big(E_a(\mathbf{x})^*\,E_b(\mathbf{x})\big),$$
//! while artificial absorbing-layer (`pml`) cells contribute strictly zero physical heat.
//! Lossy `surface_impedance` sheets ($R_s > 0$) add cycle-averaged surface power density
//! $q_{\text{surf}}(\mathbf{x}) = \frac{1}{2} G_{\text{eff}} |\mathbf{E}_\parallel(\mathbf{x})|^2$ (or
//! $\frac{1}{2} R_s |\mathbf{K}_s(\mathbf{x})|^2$ in perturbative cavity mode), assembled directly into the
//! adjacent boundary vertex load functional and projected via lumped $P_1$ mass (with exact linear
//! prolongation to $P_2$ edge midpoints when `target_order == 2`) so nodal `Q_opt` is pointwise
//! non-negative and $\int_\Omega Q_{\text{opt}}\,dV = P_{\text{vol}} + P_{\text{surf}}$ to machine precision.

use crate::cavity::{basis3_keyed, safe_curved_jacobian, Ned3, CURVED_QUAD_DEGREE};
use crate::ibc::{
    build_face_to_cells_3d, collect_sheet2d_edges_2d_mixed, collect_sheet2d_edges_2d_simplex, eval_potential_edge_2d,
    eval_potential_tri_3d, ibc_parallel_params, sheet2d_tangent_frame_3d, sorted3, PerturbativeIbcTri, EPS0, ETA0, MU0,
};
use crate::modes::C0;
use crate::ports::eval_triangle_tangential_basis_ordered;
use crate::problem::{Problem, Sheet2dBc, SurfaceImpedanceBc};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::edge_mixed::{eval_cell_basis_2d, eval_cell_basis_3d, nedelec_nloc, NedEvalScratch, NedMixed};
use femcore::fem::{build_face_to_cell_3d, geometry, quad_rule, quad_rule_collapsed, shape, CurvedTet};
use femcore::fem_mixed::{
    cell_corner_ref_coords, eval_facet_jacobian, facet_rule_mixed, lagrange_shape_into, quad_rule_mixed,
    IsoparametricCell, MixedLagrangeDofMap,
};
use femcore::linalg::csr::{c, C};
use femcore::material::M3;
use femcore::mesh::mixed::{FacetType, MixedDiscretisation};
use femcore::props::Pml;
use nyx_physics::{Discretisation, MixedToSimplexProjector};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;

const HEAT_CELL_CHUNK: usize = 512;
const FACES_TET: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

#[inline]
fn eval_quad_loss_3d(omega: f64, eps_im: &M3, eq: &[C; 3]) -> f64 {
    let mut quad_im = 0.0_f64;
    for a in 0..3 {
        for b in 0..3 {
            quad_im += eps_im[a][b] * (eq[a].conj() * eq[b]).re;
        }
    }
    (0.5 * omega * EPS0 * quad_im).max(0.0)
}

#[inline]
fn eval_sheet2d_loss_uv(sig: &[[C; 2]; 2], eu: C, ev: C) -> f64 {
    let ju = sig[0][0] * eu + sig[0][1] * ev;
    let jv = sig[1][0] * eu + sig[1][1] * ev;
    (0.5 * (eu.conj() * ju + ev.conj() * jv).re).max(0.0)
}

#[inline]
fn prolongate_p1_to_order(mesh: &Discretisation, q_p1: Vec<f64>, target_order: usize) -> Vec<f64> {
    if target_order <= 1 {
        return q_p1;
    }
    let nv = mesh.n_points();
    let edges = mesh.edges();
    let mut out = Vec::with_capacity(nv + edges.len());
    out.extend_from_slice(&q_p1);
    for &(a, b) in edges {
        out.push(0.5 * (q_p1[a as usize] + q_p1[b as usize]));
    }
    out
}

/// Computes conservative nodal `Q_opt` (`"W/m^3"`, length `mesh.n_dofs(target_order)`) and total physical
/// dissipated power `P_loss` (`"W"`) on a 3D `SimplexMesh` (`Discretisation`) for physical electric-field
/// DOFs `x_phys` (`"V/m"`).
#[allow(clippy::too_many_arguments)]
pub fn compute_3d_heat_source_simplex(
    mesh: &Discretisation,
    ned: &Ned3,
    x_phys: &[C],
    freq_hz: f64,
    pml_cells: &[bool],
    eps_phys_im: &[M3],
    mur_phys: &[f64],
    pec_labels: &[i32],
    surface_impedance: &[SurfaceImpedanceBc],
    f_ibc_default: f64,
    pert_ibc_tris: &[PerturbativeIbcTri],
    pert_inv_scale: f64,
    target_order: usize,
) -> Result<(Vec<f64>, f64), String> {
    compute_3d_heat_source_simplex_ext(
        mesh,
        ned,
        x_phys,
        freq_hz,
        pml_cells,
        eps_phys_im,
        mur_phys,
        pec_labels,
        surface_impedance,
        f_ibc_default,
        pert_ibc_tris,
        pert_inv_scale,
        &[],
        None,
        None,
        None,
        target_order,
    )
}

/// Extended 3D `SimplexMesh` `Q_opt` evaluator adding penetrable `sheet_2d` surface dissipation
/// $\frac{1}{2}\operatorname{Re}(\mathbf{E}_t^H \boldsymbol{\sigma}_{2\text{D}} \mathbf{E}_t)$ and optional
/// incident electric field `inc_e` (so scattered-field excitations evaluate total field $\mathbf{E} = \mathbf{E}_{\text{scat}} + \mathbf{E}_{\text{inc}}$).
#[allow(clippy::too_many_arguments)]
pub fn compute_3d_heat_source_simplex_ext(
    mesh: &Discretisation,
    ned: &Ned3,
    x_phys: &[C],
    freq_hz: f64,
    pml_cells: &[bool],
    eps_phys_im: &[M3],
    mur_phys: &[f64],
    pec_labels: &[i32],
    surface_impedance: &[SurfaceImpedanceBc],
    f_ibc_default: f64,
    pert_ibc_tris: &[PerturbativeIbcTri],
    pert_inv_scale: f64,
    sheet_2d: &[Sheet2dBc],
    imported_potential: Option<&[f64]>,
    temperature_field: Option<&[f64]>,
    inc_e: Option<&(dyn Fn([f64; 3]) -> [C; 3] + Sync)>,
    target_order: usize,
) -> Result<(Vec<f64>, f64), String> {
    let nc = mesh.n_cells();
    let nv = mesh.n_points();
    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    let cquad = quad_rule_collapsed(if mesh.has_qp() { 8 } else { CURVED_QUAD_DEGREE });
    let nb = ned.nb();

    let cell_ids: Vec<usize> = (0..nc).collect();
    let chunks = femcore::par::map_chunks(&cell_ids, HEAT_CELL_CHUNK, |slice| {
        let mut m_loc = vec![0.0_f64; nv];
        let mut b_loc = vec![0.0_f64; nv];
        let mut p_vol_loc = 0.0_f64;
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];

        for &ci in slice {
            let cv = mesh.cell(ci);
            let is_pml = pml_cells.get(ci).copied().unwrap_or(false);
            let e_im = &eps_phys_im[ci];
            let cell_has_loss = !is_pml && e_im.iter().flatten().any(|&v| v.abs() > 0.0);
            let gid = ned.key_pairs(mesh, ci);
            for b in 0..4 {
                for k in 0..3 {
                    gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                }
            }

            if mesh.cell_is_curved(ci) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    for a in 0..4 {
                        m_loc[cv[a] as usize] += wq * lam[a];
                    }
                    if cell_has_loss {
                        basis3_keyed(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                        let mut eq = [c(0.0, 0.0); 3];
                        for s in 0..nb {
                            let cf = x_phys[ned.cell_dofs[ci][s]];
                            for k in 0..3 {
                                eq[k] += cf * w[s][k];
                            }
                        }
                        if let Some(f_inc) = inc_e {
                            let xq = ct.position(&lam);
                            let ei = f_inc(xq);
                            for k in 0..3 {
                                eq[k] += ei[k];
                            }
                        }
                        let q_vol = eval_quad_loss_3d(omega, e_im, &eq);
                        p_vol_loc += wq * q_vol;
                        for a in 0..4 {
                            b_loc[cv[a] as usize] += wq * q_vol * lam[a];
                        }
                    }
                }
            } else {
                let vol = geo.vol[ci];
                for a in 0..4 {
                    m_loc[cv[a] as usize] += 0.25 * vol;
                }
                if cell_has_loss {
                    let pts = [
                        mesh.point(cv[0] as usize),
                        mesh.point(cv[1] as usize),
                        mesh.point(cv[2] as usize),
                        mesh.point(cv[3] as usize),
                    ];
                    for q in 0..quad.w.len() {
                        let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                        basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                        let mut eq = [c(0.0, 0.0); 3];
                        for s in 0..nb {
                            let cf = x_phys[ned.cell_dofs[ci][s]];
                            for k in 0..3 {
                                eq[k] += cf * w[s][k];
                            }
                        }
                        if let Some(f_inc) = inc_e {
                            let xq = [
                                (0..4).map(|b| lam[b] * pts[b][0]).sum(),
                                (0..4).map(|b| lam[b] * pts[b][1]).sum(),
                                (0..4).map(|b| lam[b] * pts[b][2]).sum(),
                            ];
                            let ei = f_inc(xq);
                            for k in 0..3 {
                                eq[k] += ei[k];
                            }
                        }
                        let wq = quad.w[q] * vol;
                        let q_vol = eval_quad_loss_3d(omega, e_im, &eq);
                        p_vol_loc += wq * q_vol;
                        for a in 0..4 {
                            b_loc[cv[a] as usize] += wq * q_vol * lam[a];
                        }
                    }
                }
            }
        }
        (m_loc, b_loc, p_vol_loc)
    });

    let mut m_v = vec![0.0_f64; nv];
    let mut b_v = vec![0.0_f64; nv];
    let mut p_loss = 0.0_f64;
    for (ml, bl, pl) in chunks {
        for v in 0..nv {
            m_v[v] += ml[v];
            b_v[v] += bl[v];
        }
        p_loss += pl;
    }

    // Direct IBC surface dissipation: q_surf = 0.5 * g_eff * |E_|||^2
    let any_ibc_loss = surface_impedance.iter().any(|b| b.r_s > 0.0);
    let any_sheet2d = !sheet_2d.is_empty();
    if any_ibc_loss || any_sheet2d {
        let face_to_cell = build_face_to_cell_3d(mesh);
        let quad2 = quad_rule(2, 4);
        let pert_set: BTreeSet<[u32; 3]> = pert_ibc_tris.iter().map(|pt| pt.tri).collect();
        let mut seen_pec = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if pec_labels.contains(&mesh.facet_labels[f]) {
                let fv = mesh.facet(f);
                seen_pec.insert(sorted3(fv[0], fv[1], fv[2]));
            }
        }
        let mut seen_ibc = BTreeSet::new();
        if any_ibc_loss {
            for bc in surface_impedance {
                let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { f_ibc_default };
                let omega_ref = 2.0 * PI * f_ref;
                let (_, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_ref);
                for f in 0..mesh.n_facets() {
                    if !bc.labels.contains(&mesh.facet_labels[f]) {
                        continue;
                    }
                    let fv = mesh.facet(f);
                    let key = sorted3(fv[0], fv[1], fv[2]);
                    if seen_pec.contains(&key) || !seen_ibc.insert(key) || pert_set.contains(&key) {
                        continue;
                    }
                    if g_eff <= 0.0 {
                        continue;
                    }
                    let (_, _, n_hat, bnd_dofs, weights, psi_q) =
                        eval_triangle_tangential_basis_ordered(mesh, ned, &geo, &face_to_cell, key)?;
                    let p0 = mesh.point(key[0] as usize);
                    let p1 = mesh.point(key[1] as usize);
                    let p2 = mesh.point(key[2] as usize);
                    for (q, &wq) in weights.iter().enumerate() {
                        let pq = &psi_q[q];
                        let l2 = &quad2.lam[q];
                        let mut e_par = [c(0.0, 0.0); 3];
                        for (k_slot, &gd) in bnd_dofs.iter().enumerate() {
                            let cf = x_phys[gd];
                            e_par[0] += cf * pq[k_slot][0];
                            e_par[1] += cf * pq[k_slot][1];
                            e_par[2] += cf * pq[k_slot][2];
                        }
                        if let Some(f_inc) = inc_e {
                            let xq = [
                                l2[0] * p0[0] + l2[1] * p1[0] + l2[2] * p2[0],
                                l2[0] * p0[1] + l2[1] * p1[1] + l2[2] * p2[1],
                                l2[0] * p0[2] + l2[1] * p1[2] + l2[2] * p2[2],
                            ];
                            let ei = f_inc(xq);
                            let en = ei[0] * n_hat[0] + ei[1] * n_hat[1] + ei[2] * n_hat[2];
                            e_par[0] += ei[0] - en * n_hat[0];
                            e_par[1] += ei[1] - en * n_hat[1];
                            e_par[2] += ei[2] - en * n_hat[2];
                        }
                        let e_par2 = e_par[0].norm_sqr() + e_par[1].norm_sqr() + e_par[2].norm_sqr();
                        let q_surf = 0.5 * g_eff * e_par2;
                        p_loss += wq * q_surf;
                        b_v[key[0] as usize] += wq * q_surf * l2[0];
                        b_v[key[1] as usize] += wq * q_surf * l2[1];
                        b_v[key[2] as usize] += wq * q_surf * l2[2];
                    }
                }
            }
        }

        // Penetrable `sheet_2d` surface dissipation: q_surf = 0.5 * Re(E_t^H * sigma_2D * E_t)
        if any_sheet2d {
            let mut seen_sh = seen_pec;
            seen_sh.extend(seen_ibc);
            seen_sh.extend(pert_set);
            for bc in sheet_2d {
                let f_eval = if bc.f_ref > 0.0 { bc.f_ref } else { freq_hz.max(f_ibc_default) };
                for f in 0..mesh.n_facets() {
                    if !bc.labels.contains(&mesh.facet_labels[f]) {
                        continue;
                    }
                    let fv = mesh.facet(f);
                    let key = sorted3(fv[0], fv[1], fv[2]);
                    if !seen_sh.insert(key) {
                        continue;
                    }
                    let (_, _, n_hat, bnd_dofs, weights, psi_q) =
                        eval_triangle_tangential_basis_ordered(mesh, ned, &geo, &face_to_cell, key)?;
                    let (u_hat, v_hat, _) = sheet2d_tangent_frame_3d(n_hat);
                    let p0 = mesh.point(key[0] as usize);
                    let p1 = mesh.point(key[1] as usize);
                    let p2 = mesh.point(key[2] as usize);
                    for (q, &wq) in weights.iter().enumerate() {
                        let pq = &psi_q[q];
                        let l2 = &quad2.lam[q];
                        let local_v = eval_potential_tri_3d(mesh, key, l2, imported_potential);
                        let local_t = eval_potential_tri_3d(mesh, key, l2, temperature_field);
                        let sig = bc.eval_tensor_at(f_eval, local_t, local_v)?;
                        let mut e_par = [c(0.0, 0.0); 3];
                        for (k_slot, &gd) in bnd_dofs.iter().enumerate() {
                            let cf = x_phys[gd];
                            e_par[0] += cf * pq[k_slot][0];
                            e_par[1] += cf * pq[k_slot][1];
                            e_par[2] += cf * pq[k_slot][2];
                        }
                        if let Some(f_inc) = inc_e {
                            let xq = [
                                l2[0] * p0[0] + l2[1] * p1[0] + l2[2] * p2[0],
                                l2[0] * p0[1] + l2[1] * p1[1] + l2[2] * p2[1],
                                l2[0] * p0[2] + l2[1] * p1[2] + l2[2] * p2[2],
                            ];
                            let ei = f_inc(xq);
                            let en = ei[0] * n_hat[0] + ei[1] * n_hat[1] + ei[2] * n_hat[2];
                            e_par[0] += ei[0] - en * n_hat[0];
                            e_par[1] += ei[1] - en * n_hat[1];
                            e_par[2] += ei[2] - en * n_hat[2];
                        }
                        let eu = e_par[0] * u_hat[0] + e_par[1] * u_hat[1] + e_par[2] * u_hat[2];
                        let ev = e_par[0] * v_hat[0] + e_par[1] * v_hat[1] + e_par[2] * v_hat[2];
                        let q_surf = eval_sheet2d_loss_uv(&sig, eu, ev);
                        p_loss += wq * q_surf;
                        b_v[key[0] as usize] += wq * q_surf * l2[0];
                        b_v[key[1] as usize] += wq * q_surf * l2[1];
                        b_v[key[2] as usize] += wq * q_surf * l2[2];
                    }
                }
            }
        }
    }

    // Perturbative IBC surface dissipation (3D cavity): q_surf = 0.5 * inv_scale * R_s * |K_s|^2
    if !pert_ibc_tris.is_empty() && pert_inv_scale > 0.0 && (omega * MU0) >= 1e-100 {
        let face_to_cells = build_face_to_cells_3d(mesh);
        let quad2 = quad_rule(2, 4);
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let pref = 0.5 * pert_inv_scale / (omega * omega * MU0 * MU0);
        for pt in pert_ibc_tris {
            if pt.r_s <= 0.0 {
                continue;
            }
            let Some(adj) = face_to_cells.get(&pt.tri) else { continue };
            let p0 = mesh.point(pt.tri[0] as usize);
            let p1 = mesh.point(pt.tri[1] as usize);
            let p2 = mesh.point(pt.tri[2] as usize);
            let sub = |a: &[f64], b: &[f64]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let cross = |a: [f64; 3], b: [f64; 3]| {
                [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
            };
            let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            let n_raw = cross(sub(p1, p0), sub(p2, p0));
            let n_len = dot(n_raw, n_raw).sqrt();
            if n_len <= 0.0 {
                continue;
            }
            let area = 0.5 * n_len;
            let n_hat = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];
            for q in 0..quad2.w.len() {
                let l2 = &quad2.lam[q];
                let mut area_q = area;
                let mut n_q = n_hat;
                let mut k_tilde = [c(0.0, 0.0); 3];
                for (ac_idx, &(ci, lf)) in adj.iter().enumerate() {
                    let cv = mesh.cell(ci);
                    let opp = 6 - lf[0] - lf[1] - lf[2];
                    let p_opp = mesh.point(cv[opp] as usize);
                    let sign = if dot(sub(p_opp, p0), n_hat) >= 0.0 { 1.0 } else { -1.0 };
                    let sign_inv_mur = sign / mur_phys[ci];
                    let gid = ned.key_pairs(mesh, ci);
                    let mut gs = [[0.0; 3]; 4];
                    for b in 0..4 {
                        for k in 0..3 {
                            gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                        }
                    }
                    let mut lam4 = [0.0; 4];
                    lam4[lf[0]] = l2[0];
                    lam4[lf[1]] = l2[1];
                    lam4[lf[2]] = l2[2];
                    if mesh.cell_is_curved(ci) {
                        let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
                        let gq = safe_curved_jacobian(&ct, &lam4, &gs);
                        basis3_keyed(ned.order, &lam4, &gq, &gid, &mut w, &mut cu);
                        if ac_idx == 0 {
                            let d = ct.dx_dlam(&lam4);
                            let nr = cross(
                                [d[lf[1]][0] - d[lf[0]][0], d[lf[1]][1] - d[lf[0]][1], d[lf[1]][2] - d[lf[0]][2]],
                                [d[lf[2]][0] - d[lf[0]][0], d[lf[2]][1] - d[lf[0]][1], d[lf[2]][2] - d[lf[0]][2]],
                            );
                            let nl = dot(nr, nr).sqrt();
                            if nl > 1e-30 {
                                area_q = 0.5 * nl;
                                n_q = [nr[0] / nl, nr[1] / nl, nr[2] / nl];
                                if dot(n_q, n_hat) < 0.0 {
                                    n_q = [-n_q[0], -n_q[1], -n_q[2]];
                                }
                            }
                        }
                    } else {
                        basis3_keyed(ned.order, &lam4, &gs, &gid, &mut w, &mut cu);
                    }
                    let mut curl_e = [c(0.0, 0.0); 3];
                    for s in 0..nb {
                        let cf = x_phys[ned.cell_dofs[ci][s]];
                        curl_e[0] += cf * cu[s][0];
                        curl_e[1] += cf * cu[s][1];
                        curl_e[2] += cf * cu[s][2];
                    }
                    let n_cross_c = [
                        c(n_q[1], 0.0) * curl_e[2] - c(n_q[2], 0.0) * curl_e[1],
                        c(n_q[2], 0.0) * curl_e[0] - c(n_q[0], 0.0) * curl_e[2],
                        c(n_q[0], 0.0) * curl_e[1] - c(n_q[1], 0.0) * curl_e[0],
                    ];
                    k_tilde[0] += n_cross_c[0] * sign_inv_mur;
                    k_tilde[1] += n_cross_c[1] * sign_inv_mur;
                    k_tilde[2] += n_cross_c[2] * sign_inv_mur;
                }
                let kt2 = k_tilde[0].norm_sqr() + k_tilde[1].norm_sqr() + k_tilde[2].norm_sqr();
                let wq = quad2.w[q] * area_q;
                let q_surf = pref * pt.r_s * kt2;
                p_loss += wq * q_surf;
                b_v[pt.tri[0] as usize] += wq * q_surf * l2[0];
                b_v[pt.tri[1] as usize] += wq * q_surf * l2[1];
                b_v[pt.tri[2] as usize] += wq * q_surf * l2[2];
            }
        }
    }

    let q_p1: Vec<f64> = (0..nv).map(|v| if m_v[v] > 1e-300 { (b_v[v] / m_v[v]).max(0.0) } else { 0.0 }).collect();
    let q_nodal = prolongate_p1_to_order(mesh, q_p1, target_order);
    let _ = ETA0;
    let _ = FACES_TET;
    Ok((q_nodal, p_loss))
}

/// Computes conservative nodal `Q_opt` (`"W/m^3"`, length `simplex.n_dofs(target_order)`) and total physical
/// dissipated power `P_loss` (`"W"`) on a 3D `MixedDiscretisation` projected onto `simplex`.
#[allow(clippy::too_many_arguments)]
pub fn compute_3d_heat_source_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    x_phys: &[C],
    freq_hz: f64,
    pml_cells: &[bool],
    eps_phys_im: &[M3],
    pec_labels: &[i32],
    surface_impedance: &[SurfaceImpedanceBc],
    f_ibc_default: f64,
    pert_ibc_loss_w: f64,
    pert_ibc_tris: &[PerturbativeIbcTri],
    target_order: usize,
) -> Result<(Vec<f64>, f64), String> {
    compute_3d_heat_source_mixed_ext(
        mixed,
        ned,
        simplex,
        proj,
        x_phys,
        freq_hz,
        pml_cells,
        eps_phys_im,
        pec_labels,
        surface_impedance,
        f_ibc_default,
        pert_ibc_loss_w,
        pert_ibc_tris,
        &[],
        None,
        None,
        target_order,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn compute_3d_heat_source_mixed_ext(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    x_phys: &[C],
    freq_hz: f64,
    pml_cells: &[bool],
    eps_phys_im: &[M3],
    pec_labels: &[i32],
    surface_impedance: &[SurfaceImpedanceBc],
    f_ibc_default: f64,
    pert_ibc_loss_w: f64,
    pert_ibc_tris: &[PerturbativeIbcTri],
    sheet_2d: &[Sheet2dBc],
    imported_potential: Option<&[f64]>,
    temperature_field: Option<&[f64]>,
    target_order: usize,
) -> Result<(Vec<f64>, f64), String> {
    let nv_s = simplex.n_points();
    let nsc = simplex.n_cells();
    let geo_s = geometry(simplex);
    let mut m_v = vec![0.0_f64; nv_s];
    for sc in 0..nsc {
        let cv = simplex.cell(sc);
        let vol = geo_s.vol[sc];
        for &v in cv {
            m_v[v as usize] += 0.25 * vol;
        }
    }

    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let mut b_v = vec![0.0_f64; nv_s];
    let mut p_loss = 0.0_f64;

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut cu = Vec::new();

    // Map each mixed cell to its child simplex vertices for conservative load distribution
    let mut mixed_cell_subtets: Vec<Vec<usize>> = vec![Vec::new(); mixed.n_cells()];
    for sc in 0..nsc {
        mixed_cell_subtets[proj.parent_cell(sc)].push(sc);
    }

    for mc in 0..mixed.n_cells() {
        let is_pml = pml_cells.get(mc).copied().unwrap_or(false);
        let e_im = &eps_phys_im[mc];
        if is_pml || !e_im.iter().flatten().any(|&v| v.abs() > 0.0) {
            continue;
        }
        let iso = IsoparametricCell::from_mesh(mixed, mc);
        let deg = if iso.curved { CURVED_QUAD_DEGREE } else { 4 };
        let rule = quad_rule_mixed(iso.cell_type, deg);
        let nb = nedelec_nloc(iso.cell_type, ned.order);
        let cdofs = ned.cell_dofs(mc);
        let mut p_mc = 0.0_f64;
        for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
            let jac = eval_cell_basis_3d(mixed, ned, mc, &iso, xi, &mut scratch, &mut w, &mut cu);
            let mut eq = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let cf = x_phys[cdofs[s]];
                for k in 0..3 {
                    eq[k] += cf * w[s][k];
                }
            }
            let wq = w_ref * jac.det.abs();
            p_mc += wq * eval_quad_loss_3d(omega, e_im, &eq);
        }
        p_loss += p_mc;
        let subtets = &mixed_cell_subtets[mc];
        let vol_tot: f64 = subtets.iter().map(|&sc| geo_s.vol[sc]).sum::<f64>().max(1e-300);
        for &sc in subtets {
            let share = p_mc * (geo_s.vol[sc] / vol_tot);
            for &v in simplex.cell(sc) {
                b_v[v as usize] += 0.25 * share;
            }
        }
    }

    // Direct IBC and `sheet_2d` facets on MixedDiscretisation
    if surface_impedance.iter().any(|b| b.r_s > 0.0) || !sheet_2d.is_empty() {
        let mut facet_to_cell: HashMap<usize, usize> = HashMap::new();
        let mut tri_map: HashMap<[u32; 3], usize> = HashMap::new();
        let mut quad_map: HashMap<[u32; 4], usize> = HashMap::new();
        for f in 0..mixed.n_facets() {
            let fv = mixed.facet(f);
            match mixed.facet_type(f) {
                FacetType::Tri3 => {
                    tri_map.insert(sorted3(fv[0], fv[1], fv[2]), f);
                }
                FacetType::Quad4 => {
                    let mut q = [fv[0], fv[1], fv[2], fv[3]];
                    q.sort_unstable();
                    quad_map.insert(q, f);
                }
                FacetType::Edge2 => {}
            }
        }
        for mc in 0..mixed.n_cells() {
            let cv = mixed.cell(mc);
            match mixed.cell_type(mc) {
                femcore::mesh::cells::CellType::Tet => {
                    for fl in &femcore::mesh::mixed::TET_TRI_FACES {
                        if let Some(&fid) = tri_map.get(&sorted3(cv[fl[0]], cv[fl[1]], cv[fl[2]])) {
                            facet_to_cell.entry(fid).or_insert(mc);
                        }
                    }
                }
                femcore::mesh::cells::CellType::Prism => {
                    for fl in &femcore::mesh::mixed::PRISM_TRI_FACES {
                        if let Some(&fid) = tri_map.get(&sorted3(cv[fl[0]], cv[fl[1]], cv[fl[2]])) {
                            facet_to_cell.entry(fid).or_insert(mc);
                        }
                    }
                    for fl in &femcore::mesh::mixed::PRISM_QUAD_FACES {
                        let mut q = [cv[fl[0]], cv[fl[1]], cv[fl[2]], cv[fl[3]]];
                        q.sort_unstable();
                        if let Some(&fid) = quad_map.get(&q) {
                            facet_to_cell.entry(fid).or_insert(mc);
                        }
                    }
                }
                femcore::mesh::cells::CellType::Hex => {
                    for fl in &femcore::mesh::mixed::HEX_QUAD_FACES {
                        let mut q = [cv[fl[0]], cv[fl[1]], cv[fl[2]], cv[fl[3]]];
                        q.sort_unstable();
                        if let Some(&fid) = quad_map.get(&q) {
                            facet_to_cell.entry(fid).or_insert(mc);
                        }
                    }
                }
                _ => {}
            }
        }
        let pert_set: BTreeSet<[u32; 3]> = pert_ibc_tris.iter().map(|pt| pt.tri).collect();
        let mut seen_facets: BTreeSet<Vec<u32>> = BTreeSet::new();
        for f in 0..mixed.n_facets() {
            if pec_labels.contains(&mixed.facet_labels[f]) {
                let mut key = mixed.facet(f).to_vec();
                key.sort_unstable();
                seen_facets.insert(key);
            }
        }
        for bc in surface_impedance {
            let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { f_ibc_default };
            let (_, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, 2.0 * PI * f_ref);
            for f in 0..mixed.n_facets() {
                if !bc.labels.contains(&mixed.facet_labels[f]) {
                    continue;
                }
                let fv = mixed.facet(f);
                let mut key = fv.to_vec();
                key.sort_unstable();
                if !seen_facets.insert(key) {
                    continue;
                }
                if fv.len() == 3 && pert_set.contains(&sorted3(fv[0], fv[1], fv[2])) {
                    continue;
                }
                if g_eff <= 0.0 {
                    continue;
                }
                let Some(&mc) = facet_to_cell.get(&f) else { continue };
                let iso = IsoparametricCell::from_mesh(mixed, mc);
                let corners = cell_corner_ref_coords(iso.cell_type);
                let cv = mixed.cell(mc);
                let Some(lf) = fv.iter().map(|&v| cv.iter().position(|&u| u == v)).collect::<Option<Vec<usize>>>()
                else {
                    continue;
                };
                let ftype = mixed.facet_type(f);
                let frule = facet_rule_mixed(ftype, 4);
                let cdofs = ned.cell_dofs(mc);
                let nb = cdofs.len();
                let mut p_f = 0.0_f64;
                for (uv, &w_ref) in frule.pts.iter().zip(&frule.w) {
                    let (fvals, _) =
                        femcore::fem_mixed::facet_lagrange_shape(ftype, femcore::fem_mixed::LagrangeSpace::Order1, uv);
                    let mut xi = [0.0; 3];
                    for (a, &la) in lf.iter().enumerate() {
                        for k in 0..3 {
                            xi[k] += fvals[a] * corners[la][k];
                        }
                    }
                    let (_, ds, n_hat) = eval_facet_jacobian(mixed, f, uv);
                    eval_cell_basis_3d(mixed, ned, mc, &iso, &xi, &mut scratch, &mut w, &mut cu);
                    let mut eq = [c(0.0, 0.0); 3];
                    for s in 0..nb {
                        let cf = x_phys[cdofs[s]];
                        for k in 0..3 {
                            eq[k] += cf * w[s][k];
                        }
                    }
                    let en = eq[0] * n_hat[0] + eq[1] * n_hat[1] + eq[2] * n_hat[2];
                    let e_par = [eq[0] - en * n_hat[0], eq[1] - en * n_hat[1], eq[2] - en * n_hat[2]];
                    let e_par2 = e_par[0].norm_sqr() + e_par[1].norm_sqr() + e_par[2].norm_sqr();
                    let wq = w_ref * ds;
                    let q_surf = 0.5 * g_eff * e_par2;
                    p_f += wq * q_surf;
                    for (a, &v_id) in fv.iter().enumerate() {
                        if (v_id as usize) < nv_s {
                            b_v[v_id as usize] += wq * q_surf * fvals[a];
                        }
                    }
                }
                p_loss += p_f;
            }
        }
        for bc in sheet_2d {
            let f_eval = if bc.f_ref > 0.0 { bc.f_ref } else { freq_hz.max(f_ibc_default) };
            for f in 0..mixed.n_facets() {
                if !bc.labels.contains(&mixed.facet_labels[f]) {
                    continue;
                }
                let fv = mixed.facet(f);
                let mut key = fv.to_vec();
                key.sort_unstable();
                if !seen_facets.insert(key) {
                    continue;
                }
                let Some(&mc) = facet_to_cell.get(&f) else { continue };
                let iso = IsoparametricCell::from_mesh(mixed, mc);
                let corners = cell_corner_ref_coords(iso.cell_type);
                let cv = mixed.cell(mc);
                let Some(lf) = fv.iter().map(|&v| cv.iter().position(|&u| u == v)).collect::<Option<Vec<usize>>>()
                else {
                    continue;
                };
                let ftype = mixed.facet_type(f);
                let frule = facet_rule_mixed(ftype, 4);
                let cdofs = ned.cell_dofs(mc);
                let nb = cdofs.len();
                let mut p_f = 0.0_f64;
                for (uv, &w_ref) in frule.pts.iter().zip(&frule.w) {
                    let (fvals, _) =
                        femcore::fem_mixed::facet_lagrange_shape(ftype, femcore::fem_mixed::LagrangeSpace::Order1, uv);
                    let mut xi = [0.0; 3];
                    for (a, &la) in lf.iter().enumerate() {
                        for k in 0..3 {
                            xi[k] += fvals[a] * corners[la][k];
                        }
                    }
                    let (_, ds, n_hat) = eval_facet_jacobian(mixed, f, uv);
                    let (u_hat, v_hat, _) = sheet2d_tangent_frame_3d(n_hat);
                    let local_v = imported_potential.and_then(|pot| {
                        (pot.len() >= mixed.n_points())
                            .then(|| fvals.iter().zip(fv).map(|(&sh, &vid)| sh * pot[vid as usize]).sum::<f64>())
                    });
                    let local_t = temperature_field.and_then(|tf| {
                        (tf.len() >= mixed.n_points())
                            .then(|| fvals.iter().zip(fv).map(|(&sh, &vid)| sh * tf[vid as usize]).sum::<f64>())
                    });
                    let sig = bc.eval_tensor_at(f_eval, local_t, local_v)?;
                    eval_cell_basis_3d(mixed, ned, mc, &iso, &xi, &mut scratch, &mut w, &mut cu);
                    let mut eq = [c(0.0, 0.0); 3];
                    for s in 0..nb {
                        let cf = x_phys[cdofs[s]];
                        for k in 0..3 {
                            eq[k] += cf * w[s][k];
                        }
                    }
                    let en = eq[0] * n_hat[0] + eq[1] * n_hat[1] + eq[2] * n_hat[2];
                    let e_par = [eq[0] - en * n_hat[0], eq[1] - en * n_hat[1], eq[2] - en * n_hat[2]];
                    let eu = e_par[0] * u_hat[0] + e_par[1] * u_hat[1] + e_par[2] * u_hat[2];
                    let ev = e_par[0] * v_hat[0] + e_par[1] * v_hat[1] + e_par[2] * v_hat[2];
                    let wq = w_ref * ds;
                    let q_surf = eval_sheet2d_loss_uv(&sig, eu, ev);
                    p_f += wq * q_surf;
                    for (a, &v_id) in fv.iter().enumerate() {
                        if (v_id as usize) < nv_s {
                            b_v[v_id as usize] += wq * q_surf * fvals[a];
                        }
                    }
                }
                p_loss += p_f;
            }
        }
    }

    if pert_ibc_loss_w > 0.0 && !pert_ibc_tris.is_empty() {
        p_loss += pert_ibc_loss_w;
        let share_tri = pert_ibc_loss_w / (pert_ibc_tris.len() as f64);
        for pt in pert_ibc_tris {
            for &v in &pt.tri {
                if (v as usize) < nv_s {
                    b_v[v as usize] += share_tri / 3.0;
                }
            }
        }
    }

    let mut q_p1: Vec<f64> =
        (0..nv_s).map(|v| if m_v[v] > 1e-300 { (b_v[v] / m_v[v]).max(0.0) } else { 0.0 }).collect();
    let int_s: f64 = (0..nv_s).map(|v| q_p1[v] * m_v[v]).sum();
    if int_s > 1e-300 && p_loss > 0.0 {
        let scale = p_loss / int_s;
        for qv in &mut q_p1 {
            *qv *= scale;
        }
    }
    Ok((prolongate_p1_to_order(simplex, q_p1, target_order), p_loss))
}

/// Computes conservative nodal `Q_opt` (`"W/m^3"`, normalized to 1 W axial guided Poynting power,
/// or 1 J/m electric energy when evanescent) for every mode of a 2D `SimplexMesh` complex waveguide solve.
pub fn compute_2d_waveguide_heat_source_simplex(
    mesh: &Discretisation,
    problem: &Problem,
    freq_hz: f64,
    neff: &[C],
    et: &[Vec<C>],
    phi: &[Vec<C>],
) -> Result<Vec<(Vec<f64>, f64)>, String> {
    let n_modes = neff.len();
    if n_modes == 0 {
        return Ok(Vec::new());
    }
    let nv = mesh.n_points();
    let nc = mesh.n_cells();
    let order = problem.order;
    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let k0 = omega / C0;
    let tens = problem.cell_complex_permittivity(mesh)?;
    let mur_phys = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let pml_cells: Vec<bool> =
        (0..nc).map(|ci| problem.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(Pml::is_active)).collect();

    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };

    let mut m_v = vec![0.0_f64; nv];
    for ci in 0..nc {
        let area = geo.vol[ci];
        for &v in mesh.cell(ci) {
            m_v[v as usize] += area / 3.0;
        }
    }

    let mut p_y = vec![0.0_f64; n_modes];
    let mut u_e = vec![0.0_f64; n_modes];
    let mut b_modes = vec![vec![0.0_f64; nv]; n_modes];
    let mut p_loss_raw = vec![0.0_f64; n_modes];

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];

    for ci in 0..nc {
        if pml_cells[ci] {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let (er, e_im) = (&tens[ci].0, &tens[ci].1);
        let has_im = e_im.iter().flatten().any(|&v| v.abs() > 0.0);
        let mur = mur_phys[ci];
        let inv_om_mu = 1.0 / (omega * MU0 * mur);
        let cv = mesh.cell(ci);
        mesh.cell_dofs(ci, order, &mut pd);
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (nv_s, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let wq = quad.w[q] * area;
            for j in 0..n_modes {
                let beta = neff[j] * k0;
                let j_beta = c(0.0, 1.0) * beta;
                let (mut ex, mut ez) = (c(0.0, 0.0), c(0.0, 0.0));
                for s in 0..nloc_e {
                    let cf = et[j][cd.gdofs[ci][s]];
                    ex += cf * w[s][0];
                    ez += cf * w[s][1];
                }
                let (mut phi_q, mut dphi_dx, mut dphi_dz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                for a in 0..nl {
                    let pf = phi[j][pd[a]];
                    phi_q += pf * nv_s[a];
                    dphi_dx += pf * gn[a][0];
                    dphi_dz += pf * gn[a][1];
                }
                let ey = j_beta * phi_q;
                let eq = [ex, ey, ez];
                let hx = beta * (ez + dphi_dz) * inv_om_mu;
                let hz = -beta * (ex + dphi_dx) * inv_om_mu;
                p_y[j] += 0.5 * (ez.conj() * hx - ex.conj() * hz).re * wq;
                let mut e_eps_e = 0.0_f64;
                for a in 0..3 {
                    for b in 0..3 {
                        e_eps_e += er[a][b] * (eq[a].conj() * eq[b]).re;
                    }
                }
                u_e[j] += 0.25 * EPS0 * e_eps_e * wq;
                if has_im {
                    let q_vol = eval_quad_loss_3d(omega, e_im, &eq);
                    p_loss_raw[j] += wq * q_vol;
                    for a in 0..3 {
                        b_modes[j][cv[a] as usize] += wq * q_vol * l[a];
                    }
                }
            }
        }
    }

    // 1D boundary/sheet `surface_impedance` dissipation in 2D waveguide cross-sections
    if problem.surface_impedance.iter().any(|b| b.r_s > 0.0) {
        let mut edge_to_cell: BTreeMap<(u32, u32), (usize, usize, usize)> = BTreeMap::new();
        for ci in 0..nc {
            let cv = mesh.cell(ci);
            for (la, lb) in [(0, 1), (1, 2), (0, 2)] {
                let key = (cv[la].min(cv[lb]), cv[la].max(cv[lb]));
                edge_to_cell.entry(key).or_insert((ci, la, lb));
            }
        }
        let pec_set: BTreeSet<(u32, u32)> = (0..mesh.n_facets())
            .filter(|&f| problem.pec.contains(&mesh.facet_labels[f]))
            .map(|f| {
                let fv = mesh.facet(f);
                (fv[0].min(fv[1]), fv[0].max(fv[1]))
            })
            .collect();
        let mut seen_edges = BTreeSet::new();
        const GP: [f64; 3] = [0.112_701_665_379_258_3, 0.5, 0.887_298_334_620_741_7];
        const GW: [f64; 3] = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0];
        for bc in &problem.surface_impedance {
            if bc.r_s <= 0.0 {
                continue;
            }
            let f_ref = if bc.f_ref > 0.0 { bc.f_ref } else { freq_hz };
            let (_, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, 2.0 * PI * f_ref);
            if g_eff <= 0.0 {
                continue;
            }
            for f in 0..mesh.n_facets() {
                if !bc.labels.contains(&mesh.facet_labels[f]) {
                    continue;
                }
                let fv = mesh.facet(f);
                let key = (fv[0].min(fv[1]), fv[0].max(fv[1]));
                if pec_set.contains(&key) || !seen_edges.insert(key) {
                    continue;
                }
                let Some(&(ci, la, lb)) = edge_to_cell.get(&key) else { continue };
                let cv = mesh.cell(ci);
                let p0 = mesh.point(cv[la] as usize);
                let p1 = mesh.point(cv[lb] as usize);
                let dx = p1[0] - p0[0];
                let dz = p1[1] - p0[1];
                let len = dx.hypot(dz);
                if len <= 1e-30 {
                    continue;
                }
                let t_hat = [dx / len, dz / len];
                let g = tri_g(&geo, ci);
                mesh.cell_dofs(ci, order, &mut pd);
                for (&t, &gw) in GP.iter().zip(&GW) {
                    let mut lam = [0.0; 3];
                    lam[la] = 1.0 - t;
                    lam[lb] = t;
                    edge_basis(&lam, &g, &cd.sign[ci], order, &mut w, &mut curl);
                    let (nv_s, _) = shape(&lam, order);
                    let ds = gw * len;
                    for j in 0..n_modes {
                        let beta = neff[j] * k0;
                        let j_beta = c(0.0, 1.0) * beta;
                        let (mut ex, mut ez, mut phi_q) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                        for s in 0..nloc_e {
                            let cf = et[j][cd.gdofs[ci][s]];
                            ex += cf * w[s][0];
                            ez += cf * w[s][1];
                        }
                        for a in 0..nl {
                            phi_q += phi[j][pd[a]] * nv_s[a];
                        }
                        let e_tan = ex * t_hat[0] + ez * t_hat[1];
                        let ey = j_beta * phi_q;
                        let q_line = 0.5 * g_eff * (e_tan.norm_sqr() + ey.norm_sqr());
                        p_loss_raw[j] += ds * q_line;
                        b_modes[j][cv[la] as usize] += ds * q_line * (1.0 - t);
                        b_modes[j][cv[lb] as usize] += ds * q_line * t;
                    }
                }
            }
        }
    }

    // 1D boundary/interface `sheet_2d` dissipation in 2D waveguide cross-sections
    if !problem.sheet_2d.is_empty() {
        let mut edge_to_cell: BTreeMap<(u32, u32), (usize, usize, usize)> = BTreeMap::new();
        for ci in 0..nc {
            let cv = mesh.cell(ci);
            for (la, lb) in [(0, 1), (1, 2), (0, 2)] {
                let key = (cv[la].min(cv[lb]), cv[la].max(cv[lb]));
                edge_to_cell.entry(key).or_insert((ci, la, lb));
            }
        }
        let mut seen_edges: BTreeSet<(u32, u32)> = (0..mesh.n_facets())
            .filter(|&f| problem.pec.contains(&mesh.facet_labels[f]))
            .map(|f| {
                let fv = mesh.facet(f);
                (fv[0].min(fv[1]), fv[0].max(fv[1]))
            })
            .collect();
        const GP: [f64; 3] = [0.112_701_665_379_258_3, 0.5, 0.887_298_334_620_741_7];
        const GW: [f64; 3] = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0];
        for bc in &problem.sheet_2d {
            let f_eval = if bc.f_ref > 0.0 { bc.f_ref } else { freq_hz };
            for key in collect_sheet2d_edges_2d_simplex(mesh, bc, &mut seen_edges) {
                let Some(&(ci, la, lb)) = edge_to_cell.get(&key) else { continue };
                let cv = mesh.cell(ci);
                let p0 = mesh.point(cv[la] as usize);
                let p1 = mesh.point(cv[lb] as usize);
                let dx = p1[0] - p0[0];
                let dz = p1[1] - p0[1];
                let len = dx.hypot(dz);
                if len <= 1e-30 {
                    continue;
                }
                let sign = if dx > 1e-12 || (dx.abs() <= 1e-12 && dz >= 0.0) { 1.0 } else { -1.0 };
                let tau = [sign * dx / len, sign * dz / len];
                let g = tri_g(&geo, ci);
                mesh.cell_dofs(ci, order, &mut pd);
                for (&t, &gw) in GP.iter().zip(&GW) {
                    let l0 = 1.0 - t;
                    let l1 = t;
                    let mut lam = [0.0; 3];
                    lam[la] = l0;
                    lam[lb] = l1;
                    let local_v =
                        eval_potential_edge_2d(mesh, cv[la], cv[lb], l0, l1, problem.imported_potential.as_deref());
                    let local_t =
                        eval_potential_edge_2d(mesh, cv[la], cv[lb], l0, l1, problem.temperature_field.as_deref());
                    let sig = bc.eval_tensor_at(f_eval, local_t, local_v)?;
                    edge_basis(&lam, &g, &cd.sign[ci], order, &mut w, &mut curl);
                    let (nv_s, _) = shape(&lam, order);
                    let ds = gw * len;
                    for j in 0..n_modes {
                        let beta = neff[j] * k0;
                        let j_beta = c(0.0, 1.0) * beta;
                        let (mut ex, mut ez, mut phi_q) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                        for s in 0..nloc_e {
                            let cf = et[j][cd.gdofs[ci][s]];
                            ex += cf * w[s][0];
                            ez += cf * w[s][1];
                        }
                        for a in 0..nl {
                            phi_q += phi[j][pd[a]] * nv_s[a];
                        }
                        let eu = ex * tau[0] + ez * tau[1];
                        let ev = j_beta * phi_q;
                        let q_line = eval_sheet2d_loss_uv(&sig, eu, ev);
                        p_loss_raw[j] += ds * q_line;
                        b_modes[j][cv[la] as usize] += ds * q_line * l0;
                        b_modes[j][cv[lb] as usize] += ds * q_line * l1;
                    }
                }
            }
        }
    }

    let mut out = Vec::with_capacity(n_modes);
    for j in 0..n_modes {
        let norm_power = if p_y[j].abs() > 1e-30 { p_y[j].abs() } else { u_e[j].max(1e-300) };
        let scale2 = 1.0 / norm_power;
        let q_p1: Vec<f64> =
            (0..nv).map(|v| if m_v[v] > 1e-300 { (b_modes[j][v] * scale2 / m_v[v]).max(0.0) } else { 0.0 }).collect();
        out.push((prolongate_p1_to_order(mesh, q_p1, order), p_loss_raw[j] * scale2));
    }
    Ok(out)
}

/// Mixed-cell (`Tri + Quad`) 2D waveguide `Q_opt` (`"W/m^3"`, 1 W guided power normalization) projected onto `simplex`.
#[allow(clippy::too_many_arguments)]
pub fn compute_2d_waveguide_heat_source_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    phi_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    problem: &Problem,
    freq_hz: f64,
    neff: &[C],
    et: &[Vec<C>],
    phi: &[Vec<C>],
    mats: &crate::mixed_ops::MixedCellMaterials,
) -> Result<Vec<(Vec<f64>, f64)>, String> {
    let n_modes = neff.len();
    if n_modes == 0 {
        return Ok(Vec::new());
    }
    let nv_s = simplex.n_points();
    let nsc = simplex.n_cells();
    let geo_s = geometry(simplex);
    let mut m_v = vec![0.0_f64; nv_s];
    let mut mixed_subtris: Vec<Vec<usize>> = vec![Vec::new(); mixed.n_cells()];
    for sc in 0..nsc {
        mixed_subtris[proj.parent_cell(sc)].push(sc);
        for &v in simplex.cell(sc) {
            m_v[v as usize] += geo_s.vol[sc] / 3.0;
        }
    }

    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let k0 = omega / C0;
    let pml_cells: Vec<bool> = (0..mixed.n_cells())
        .map(|ci| problem.materials.pml.get(&mixed.cell_labels[ci]).is_some_and(Pml::is_active))
        .collect();

    let mut p_y = vec![0.0_f64; n_modes];
    let mut u_e = vec![0.0_f64; n_modes];
    let mut b_modes = vec![vec![0.0_f64; nv_s]; n_modes];
    let mut p_loss_raw = vec![0.0_f64; n_modes];

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut curl = Vec::new();
    let mut nv_sh = Vec::new();
    let mut dxi_sh = Vec::new();
    let mut gn = Vec::new();

    for mc in 0..mixed.n_cells() {
        if pml_cells[mc] {
            continue;
        }
        let iso = IsoparametricCell::from_mesh(mixed, mc);
        let rule = quad_rule_mixed(iso.cell_type, if iso.curved { 6 } else { 4 });
        let cdofs = ned.cell_dofs(mc);
        let pd = phi_map.cell_dofs(mc);
        let nloc_e = cdofs.len();
        let nl = pd.len();
        gn.resize(nl, [0.0; 2]);
        let er = &mats.eps_phys_re[mc];
        let e_im = &mats.eps_phys_im[mc];
        let has_im = e_im.iter().flatten().any(|&v| v.abs() > 0.0);
        let inv_om_mu = 1.0 / (omega * MU0 * mats.mur_phys[mc]);
        let mut cell_loss = vec![0.0_f64; n_modes];

        for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
            let jac = eval_cell_basis_2d(mixed, ned, mc, &iso, xi, &mut scratch, &mut w, &mut curl);
            lagrange_shape_into(iso.cell_type, phi_map.space, xi, &mut nv_sh, &mut dxi_sh);
            for a in 0..nl {
                let g3 = jac.map_grad(dxi_sh[a]);
                gn[a] = [g3[0], g3[1]];
            }
            let wq = w_ref * jac.det.abs();
            for j in 0..n_modes {
                let beta = neff[j] * k0;
                let j_beta = c(0.0, 1.0) * beta;
                let (mut ex, mut ez) = (c(0.0, 0.0), c(0.0, 0.0));
                for s in 0..nloc_e {
                    let cf = et[j][cdofs[s]];
                    ex += cf * w[s][0];
                    ez += cf * w[s][1];
                }
                let (mut phi_q, mut dphi_dx, mut dphi_dz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                for a in 0..nl {
                    let pf = phi[j][pd[a]];
                    phi_q += pf * nv_sh[a];
                    dphi_dx += pf * gn[a][0];
                    dphi_dz += pf * gn[a][1];
                }
                let ey = j_beta * phi_q;
                let eq = [ex, ey, ez];
                let hx = beta * (ez + dphi_dz) * inv_om_mu;
                let hz = -beta * (ex + dphi_dx) * inv_om_mu;
                p_y[j] += 0.5 * (ez.conj() * hx - ex.conj() * hz).re * wq;
                let mut e_eps_e = 0.0_f64;
                for a in 0..3 {
                    for b in 0..3 {
                        e_eps_e += er[a][b] * (eq[a].conj() * eq[b]).re;
                    }
                }
                u_e[j] += 0.25 * EPS0 * e_eps_e * wq;
                if has_im {
                    cell_loss[j] += wq * eval_quad_loss_3d(omega, e_im, &eq);
                }
            }
        }
        if has_im {
            let subtris = &mixed_subtris[mc];
            let area_tot: f64 = subtris.iter().map(|&sc| geo_s.vol[sc]).sum::<f64>().max(1e-300);
            for j in 0..n_modes {
                p_loss_raw[j] += cell_loss[j];
                for &sc in subtris {
                    let share = cell_loss[j] * (geo_s.vol[sc] / area_tot);
                    for &v in simplex.cell(sc) {
                        b_modes[j][v as usize] += share / 3.0;
                    }
                }
            }
        }
    }

    if !problem.sheet_2d.is_empty() {
        let mut edge_to_cell: HashMap<(u32, u32), usize> = HashMap::with_capacity(mixed.n_edges());
        for ci in 0..mixed.n_cells() {
            let cv = mixed.cell(ci);
            let le: &[(usize, usize)] = match mixed.cell_type(ci) {
                femcore::mesh::cells::CellType::Tri => &femcore::mesh::mixed::TRI_EDGES,
                femcore::mesh::cells::CellType::Quad => &femcore::mesh::mixed::QUAD_EDGES,
                _ => &[],
            };
            for &(i, j) in le {
                edge_to_cell.entry((cv[i].min(cv[j]), cv[i].max(cv[j]))).or_insert(ci);
            }
        }
        let mut seen_edges: BTreeSet<(u32, u32)> = (0..mixed.n_facets())
            .filter(|&f| problem.pec.contains(&mixed.facet_labels[f]))
            .map(|f| {
                let fv = mixed.facet(f);
                (fv[0].min(fv[1]), fv[0].max(fv[1]))
            })
            .collect();
        for bc in &problem.sheet_2d {
            let f_eval = if bc.f_ref > 0.0 { bc.f_ref } else { freq_hz };
            for key in collect_sheet2d_edges_2d_mixed(mixed, bc, &mut seen_edges) {
                let Some(&mc) = edge_to_cell.get(&key) else { continue };
                let fv = [key.0, key.1];
                let iso = IsoparametricCell::from_mesh(mixed, mc);
                let cv = mixed.cell(mc);
                let corners_ref = cell_corner_ref_coords(iso.cell_type);
                let Some(la) = cv.iter().position(|&v| v == fv[0]) else { continue };
                let Some(lb) = cv.iter().position(|&v| v == fv[1]) else { continue };
                let frule = facet_rule_mixed(FacetType::Edge2, if iso.curved { 6 } else { 4 });
                let cdofs = ned.cell_dofs(mc);
                let pd = phi_map.cell_dofs(mc);
                let p0 = mixed.point(fv[0] as usize);
                let p1 = mixed.point(fv[1] as usize);
                let dx = p1[0] - p0[0];
                let dz = p1[1] - p0[1];
                let len = dx.hypot(dz);
                if len <= 1e-30 {
                    continue;
                }
                let sign = if dx > 1e-12 || (dx.abs() <= 1e-12 && dz >= 0.0) { 1.0 } else { -1.0 };
                let tau = [sign * dx / len, sign * dz / len];
                for (xi_f, &wf) in frule.pts.iter().zip(&frule.w) {
                    let (fn1, _) = femcore::fem_mixed::facet_lagrange_shape(
                        FacetType::Edge2,
                        femcore::fem_mixed::LagrangeSpace::Order1,
                        xi_f,
                    );
                    let xi_cell = [
                        fn1[0] * corners_ref[la][0] + fn1[1] * corners_ref[lb][0],
                        fn1[0] * corners_ref[la][1] + fn1[1] * corners_ref[lb][1],
                    ];
                    let local_v = problem.imported_potential.as_deref().and_then(|pot| {
                        (pot.len() >= mixed.n_points())
                            .then_some(fn1[0] * pot[fv[0] as usize] + fn1[1] * pot[fv[1] as usize])
                    });
                    let local_t = problem.temperature_field.as_deref().and_then(|tf| {
                        (tf.len() >= mixed.n_points())
                            .then_some(fn1[0] * tf[fv[0] as usize] + fn1[1] * tf[fv[1] as usize])
                    });
                    let sig = bc.eval_tensor_at(f_eval, local_t, local_v)?;
                    eval_cell_basis_2d(mixed, ned, mc, &iso, &xi_cell, &mut scratch, &mut w, &mut curl);
                    lagrange_shape_into(iso.cell_type, phi_map.space, &xi_cell, &mut nv_sh, &mut dxi_sh);
                    let ds = wf * len;
                    for j in 0..n_modes {
                        let beta = neff[j] * k0;
                        let j_beta = c(0.0, 1.0) * beta;
                        let (mut ex, mut ez, mut phi_q) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                        for (s, &gd) in cdofs.iter().enumerate() {
                            let cf = et[j][gd];
                            ex += cf * w[s][0];
                            ez += cf * w[s][1];
                        }
                        for (a, &gp) in pd.iter().enumerate() {
                            phi_q += phi[j][gp] * nv_sh[a];
                        }
                        let eu = ex * tau[0] + ez * tau[1];
                        let ev = j_beta * phi_q;
                        let q_line = eval_sheet2d_loss_uv(&sig, eu, ev);
                        p_loss_raw[j] += ds * q_line;
                        if (fv[0] as usize) < nv_s {
                            b_modes[j][fv[0] as usize] += ds * q_line * fn1[0];
                        }
                        if (fv[1] as usize) < nv_s {
                            b_modes[j][fv[1] as usize] += ds * q_line * fn1[1];
                        }
                    }
                }
            }
        }
    }

    let mut out = Vec::with_capacity(n_modes);
    for j in 0..n_modes {
        let norm_power = if p_y[j].abs() > 1e-30 { p_y[j].abs() } else { u_e[j].max(1e-300) };
        let scale2 = 1.0 / norm_power;
        let q_p1: Vec<f64> =
            (0..nv_s).map(|v| if m_v[v] > 1e-300 { (b_modes[j][v] * scale2 / m_v[v]).max(0.0) } else { 0.0 }).collect();
        out.push((prolongate_p1_to_order(simplex, q_p1, problem.order), p_loss_raw[j] * scale2));
    }
    Ok(out)
}
