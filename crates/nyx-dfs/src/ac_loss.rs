// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! AC quasi-static dielectric loss (`tan_delta`, `eps_r_im`, `sigma`, `freq`) for `dfs` `electrostatic`:
//! assembles the real symmetric conductance loss matrix `K_cond_loss` (scalar, anisotropic tensor, and curved 3D),
//! computes the terminal conductance matrix `G_ij(omega)`, the positive 2nd-order capacitance shift
//! `Delta C_ij^(2) = +(K_loss v_i)_f^T K_ff^-1 (K_loss v_j)_f > 0` (and converged complex Schur complement
//! when `max_c tan_delta_eff(c) > 1e-2`), and the safe dielectric quality factor matrix `Q_ij = omega |C_ij| / |G_ij|`.

use crate::common::plane_tensors;
use crate::electro::{assemble_electrostatics, enrich_quarter_point_mesh, solve_electrostatics, EPS0};
use crate::problem::{BcKind, Problem};
use femcore::fem::*;
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::DirectSolver;
use femcore::material::{self, M3};
use nyx_physics::{CapacitanceMatrixData, Discretisation, SolveResult};
use std::collections::BTreeSet;
use std::f64::consts::PI;

/// Output of AC quasi-static loss and complex-permittivity Schur complement evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct AcLossResult {
    /// Operating frequency `f` in Hz (`None` or `Some(0.0)` at DC).
    pub freq_hz: Option<f64>,
    /// Maximum effective loss tangent `max_c tan_delta_eff(c)` across all mesh cells.
    pub max_tan_delta_eff: f64,
    /// AC terminal conductance matrix `G_ij(omega)` in Siemens (or S/m in 2D); `None` when total loss is zero.
    pub conductance_matrix_s: Option<Vec<Vec<f64>>>,
    /// Per-entry dielectric quality factor `Q_ij = omega |C_ij| / |G_ij|`; `None` when `freq == 0` or lossless,
    /// and `None` at entry `(i, j)` when `|G_ij| <= 1e-300` (never `f64::INFINITY` or `NaN`).
    pub q_dielectric_matrix: Option<Vec<Vec<Option<f64>>>>,
    /// Pure 2nd-order positive capacitance perturbation `Delta C_ij^(2) = +(K_loss v_i)_f^T K_ff^-1 (K_loss v_j)_f`
    /// evaluated on free DOFs `f` when `max_tan_delta_eff > 1e-2` and `freq > 0`.
    pub delta_c_second_order_f: Option<Vec<Vec<f64>>>,
    /// Complex-Schur-corrected Maxwell capacitance matrix `C_ij(omega)` in Farads (or F/m in 2D).
    /// Equals the unperturbed `base_maxwell_matrix_f` when `max_tan_delta_eff <= 1e-2` or `freq == 0`.
    pub corrected_maxwell_matrix_f: Vec<Vec<f64>>,
    /// Assembled real symmetric loss-conductance matrix `K_cond_loss` (`None` when lossless).
    pub k_cond_loss: Option<Csr>,
}

fn assemble_laplace_curved_3d(mesh: &Discretisation, order: usize, coeff: &[f64]) -> Csr {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(6);
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * 3];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.iter_mut().for_each(|x| *x = 0.0);
            if mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    let dn = &cshapes[q];
                    for a in 0..nloc {
                        for k in 0..3 {
                            let mut s = 0.0;
                            for b in 0..4 {
                                s += dn[a * 4 + b] * gq[b][k];
                            }
                            grad[a * 3 + k] = s;
                        }
                    }
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let s = grad[a * 3] * grad[b * 3]
                                + grad[a * 3 + 1] * grad[b * 3 + 1]
                                + grad[a * 3 + 2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += wq * s;
                        }
                    }
                }
                let sc = coeff[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            } else {
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, c, q, &mut grad);
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let s = grad[a * 3] * grad[b * 3]
                                + grad[a * 3 + 1] * grad[b * 3 + 1]
                                + grad[a * 3 + 2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += re.w[q] * s;
                        }
                    }
                }
                let sc = geo.vol[c] * coeff[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            }
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

fn assemble_laplace_tensor_curved_3d(mesh: &Discretisation, order: usize, kt: &[M3]) -> Csr {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(6);
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * 3];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.iter_mut().for_each(|x| *x = 0.0);
            if mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    let dn = &cshapes[q];
                    for a in 0..nloc {
                        for k in 0..3 {
                            let mut s = 0.0;
                            for b in 0..4 {
                                s += dn[a * 4 + b] * gq[b][k];
                            }
                            grad[a * 3 + k] = s;
                        }
                    }
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    for a in 0..nloc {
                        let mut kg = [0.0; 3];
                        for i in 0..3 {
                            for j in 0..3 {
                                kg[i] += kt[c][i][j] * grad[a * 3 + j];
                            }
                        }
                        for b in 0..nloc {
                            let s = kg[0] * grad[b * 3] + kg[1] * grad[b * 3 + 1] + kg[2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += wq * s;
                        }
                    }
                }
            } else {
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, c, q, &mut grad);
                    for a in 0..nloc {
                        let mut kg = [0.0; 3];
                        for i in 0..3 {
                            for j in 0..3 {
                                kg[i] += kt[c][i][j] * grad[a * 3 + j];
                            }
                        }
                        for b in 0..nloc {
                            let s = kg[0] * grad[b * 3] + kg[1] * grad[b * 3 + 1] + kg[2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += re.w[q] * s;
                        }
                    }
                }
                let sc = geo.vol[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            }
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

/// Evaluates the per-cell real permittivity tensor `eps_r(c)`, effective loss conductivity tensor
/// `sigma_eff(c) = omega * eps_0 * eps_r_im(c) + sigma(c)`, and maximum effective loss tangent `max_c tan_delta_eff(c)`.
///
/// Returns `(Option<Csr>, max_tan_delta_eff)`, where the `Csr` is `None` when total loss is zero across all cells.
pub fn assemble_cond_loss_matrix(mesh: &Discretisation, spec: &Problem, freq_hz: Option<f64>) -> (Option<Csr>, f64) {
    let order = spec.order.max(1);
    let nc = mesh.n_cells();
    let omega = freq_hz.filter(|&f| f > 0.0 && f.is_finite()).map_or(0.0, |f| 2.0 * PI * f);
    let has_curved = mesh.dim == 3 && mesh.has_curved();

    let mut has_tensor_loss = false;
    for an in spec.materials.aniso.values() {
        if an.permittivity_im.is_some() || an.conductivity.is_some() || an.permittivity.is_some() {
            has_tensor_loss = true;
            break;
        }
    }

    let mut max_tan_delta_eff = 0.0_f64;
    let mut max_cond_loss = 0.0_f64;

    if !has_tensor_loss {
        let mut cond_coeff = vec![0.0; nc];
        for c in 0..nc {
            let lab = mesh.cell_labels[c];
            let eps_r = spec.materials.get(lab, "eps_r").unwrap_or(1.0);
            let td_scalar = spec.materials.get(lab, "tan_delta").unwrap_or(0.0);
            let eps_im = spec.materials.get(lab, "eps_r_im").unwrap_or(eps_r * td_scalar);
            let sigma = spec.materials.get(lab, "sigma").unwrap_or(0.0);
            let c_loss = omega * EPS0 * eps_im + sigma;
            cond_coeff[c] = c_loss;
            max_cond_loss = max_cond_loss.max(c_loss.abs());
            if eps_r.abs() > 0.0 {
                let td_eff = if omega > 0.0 { c_loss / (omega * EPS0 * eps_r) } else { eps_im / eps_r };
                max_tan_delta_eff = max_tan_delta_eff.max(td_eff.abs());
            }
        }
        if max_cond_loss <= 0.0 {
            return (None, 0.0);
        }
        let k_cond = if has_curved {
            assemble_laplace_curved_3d(mesh, order, &cond_coeff)
        } else {
            assemble_laplace(mesh, order, &cond_coeff)
        };
        return (Some(k_cond), max_tan_delta_eff);
    }

    let mut loss_tensors = vec![[[0.0; 3]; 3]; nc];
    for c in 0..nc {
        let lab = mesh.cell_labels[c];
        let eps_r_scalar = spec.materials.get(lab, "eps_r").unwrap_or(1.0);
        let td_scalar = spec.materials.get(lab, "tan_delta").unwrap_or(0.0);
        let eps_im_scalar = spec.materials.get(lab, "eps_r_im").unwrap_or(0.0);
        let sigma_scalar = spec.materials.get(lab, "sigma").unwrap_or(0.0);

        let an = spec.materials.aniso.get(&lab);
        let euler = an.map_or([0.0; 3], |a| a.euler_deg);

        let eps_re_m = an.and_then(|a| a.permittivity.as_ref()).map(|t| material::tensor_global(t, euler)).unwrap_or([
            [eps_r_scalar, 0.0, 0.0],
            [0.0, eps_r_scalar, 0.0],
            [0.0, 0.0, eps_r_scalar],
        ]);

        let eps_im_m = if let Some(t_im) = an.and_then(|a| a.permittivity_im.as_ref()) {
            material::tensor_global(t_im, euler)
        } else if eps_im_scalar != 0.0 {
            [[eps_im_scalar, 0.0, 0.0], [0.0, eps_im_scalar, 0.0], [0.0, 0.0, eps_im_scalar]]
        } else {
            let mut m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] = eps_re_m[i][j] * td_scalar;
                }
            }
            m
        };

        let sigma_m = if let Some(t_sig) = an.and_then(|a| a.conductivity.as_ref()) {
            material::tensor_global(t_sig, euler)
        } else {
            [[sigma_scalar, 0.0, 0.0], [0.0, sigma_scalar, 0.0], [0.0, 0.0, sigma_scalar]]
        };

        let mut ct = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let val = omega * EPS0 * eps_im_m[i][j] + sigma_m[i][j];
                ct[i][j] = val;
                max_cond_loss = max_cond_loss.max(val.abs());
            }
        }
        for d in 0..3 {
            let er_d = eps_re_m[d][d].abs();
            if er_d > 0.0 {
                let td_d =
                    if omega > 0.0 { ct[d][d].abs() / (omega * EPS0 * er_d) } else { eps_im_m[d][d].abs() / er_d };
                max_tan_delta_eff = max_tan_delta_eff.max(td_d);
            }
        }
        loss_tensors[c] = ct;
    }

    if max_cond_loss <= 0.0 {
        return (None, 0.0);
    }

    let t_plane = plane_tensors(mesh, loss_tensors);
    let k_cond = if has_curved {
        assemble_laplace_tensor_curved_3d(mesh, order, &t_plane)
    } else {
        assemble_laplace_tensor(mesh, order, &t_plane)
    };
    (Some(k_cond), max_tan_delta_eff)
}

/// Computes AC terminal conductance `G_ij(omega)`, positive 2nd-order capacitance shift `Delta C_ij^(2)`,
/// converged complex Schur-complement capacitance `C_ij(omega)` (when `max_tan_delta_eff > 1e-2`),
/// and safe `q_dielectric_matrix` from the solved electrostatic terminal potentials `v_all`.
pub fn compute_ac_loss(
    mesh: &Discretisation,
    spec: &Problem,
    freq_hz: Option<f64>,
    k: &Csr,
    free: &[usize],
    v_all: &[Vec<f64>],
    base_maxwell_matrix_f: &[Vec<f64>],
) -> Result<AcLossResult, String> {
    let n_term = v_all.len();
    let omega = freq_hz.filter(|&f| f > 0.0 && f.is_finite()).map_or(0.0, |f| 2.0 * PI * f);
    let (k_cond_opt, max_tan_delta_eff) = assemble_cond_loss_matrix(mesh, spec, freq_hz);

    let Some(k_cond) = k_cond_opt else {
        return Ok(AcLossResult {
            freq_hz,
            max_tan_delta_eff: 0.0,
            conductance_matrix_s: None,
            q_dielectric_matrix: None,
            delta_c_second_order_f: None,
            corrected_maxwell_matrix_f: base_maxwell_matrix_f.to_vec(),
            k_cond_loss: None,
        });
    };

    let kcond_v_all: Vec<Vec<f64>> = v_all.iter().map(|v| k_cond.matvec(v)).collect();
    let mut g_mat = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in i..n_term {
            let gij = femcore::linalg::vecops::dotu(&v_all[i], &kcond_v_all[j]);
            let gji = if i == j { gij } else { femcore::linalg::vecops::dotu(&v_all[j], &kcond_v_all[i]) };
            let sym = 0.5 * (gij + gji);
            g_mat[i][j] = sym;
            g_mat[j][i] = sym;
        }
    }

    let mut c_mat = base_maxwell_matrix_f.to_vec();
    let mut delta_c_second_order_f = None;

    // When inhomogeneous loss is strong (max_c tan_delta_eff(c) > 1e-2) and omega > 0, evaluate:
    // 1. The positive 2nd-order capacitance shift:
    //      Delta C_ij^(2) = + (K_loss v_i)_f^T K_ff^-1 (K_loss v_j)_f > 0,   K_loss = K_cond_loss / omega
    // 2. The converged coupled real/imaginary Schur-complement correction on free DOFs f:
    //      K_ff x_{j,f} + K_{loss,ff} y_{j,f} = 0
    //      K_ff y_{j,f} - K_{loss,ff} x_{j,f} = (K_loss v_j)_f
    if max_tan_delta_eff > 1.0e-2 && omega > 0.0 && !free.is_empty() && n_term > 0 {
        let inv_omega = 1.0 / omega;
        let n_free = free.len();
        let r_free_list: Vec<Vec<f64>> =
            (0..n_term).map(|j| free.iter().map(|&fi| kcond_v_all[j][fi] * inv_omega).collect()).collect();
        let any_rhs_nonzero = r_free_list.iter().any(|r| r.iter().any(|&v| v.abs() > 1e-300));
        if any_rhs_nonzero {
            let kff = k.submatrix(free, free);
            let kloss_ff = k_cond.submatrix(free, free);
            let solver = DirectSolver::new(&kff, true)?;

            // Step 1: Pure 2nd-order perturbation y^(1) = K_ff^-1 r_free
            let y1_list = solver.solve_many(&r_free_list);
            let mut dc2 = vec![vec![0.0; n_term]; n_term];
            for i in 0..n_term {
                for j in i..n_term {
                    let d_ij = femcore::linalg::vecops::dotu(&r_free_list[i], &y1_list[j]);
                    let d_ji = if i == j { d_ij } else { femcore::linalg::vecops::dotu(&r_free_list[j], &y1_list[i]) };
                    let sym = 0.5 * (d_ij + d_ji);
                    dc2[i][j] = sym;
                    dc2[j][i] = sym;
                }
            }
            delta_c_second_order_f = Some(dc2);

            // Step 2: Converge coupled real/imaginary Schur system for exact high-tan_delta admittance Y(omega) = G(omega) + i omega C(omega)
            let mut x_list = vec![vec![0.0; n_free]; n_term];
            let mut y_list = y1_list;
            for _iter in 0..16 {
                // x_new = -K_ff^-1 (K_{loss,ff} y)
                let kloss_y_list: Vec<Vec<f64>> = y_list
                    .iter()
                    .map(|y| {
                        let ky = kloss_ff.matvec(y);
                        ky.into_iter().map(|v| -v * inv_omega).collect()
                    })
                    .collect();
                let x_new = solver.solve_many(&kloss_y_list);

                // y_new = K_ff^-1 (r_free + K_{loss,ff} x_new)
                let rhs_y_list: Vec<Vec<f64>> = (0..n_term)
                    .map(|j| {
                        let kx = kloss_ff.matvec(&x_new[j]);
                        (0..n_free).map(|a| r_free_list[j][a] + kx[a] * inv_omega).collect()
                    })
                    .collect();
                let y_new = solver.solve_many(&rhs_y_list);

                let mut max_diff = 0.0_f64;
                let mut max_norm = 1e-300_f64;
                for j in 0..n_term {
                    for a in 0..n_free {
                        max_diff =
                            max_diff.max((x_new[j][a] - x_list[j][a]).abs()).max((y_new[j][a] - y_list[j][a]).abs());
                        max_norm = max_norm.max(x_new[j][a].abs()).max(y_new[j][a].abs());
                    }
                }
                x_list = x_new;
                y_list = y_new;
                if max_diff <= 1e-13 * max_norm {
                    break;
                }
            }

            for i in 0..n_term {
                for j in i..n_term {
                    let dc_ij = femcore::linalg::vecops::dotu(&r_free_list[i], &y_list[j]);
                    let dc_ji = if i == j { dc_ij } else { femcore::linalg::vecops::dotu(&r_free_list[j], &y_list[i]) };
                    let dc_sym = 0.5 * (dc_ij + dc_ji);
                    c_mat[i][j] = base_maxwell_matrix_f[i][j] + dc_sym;
                    c_mat[j][i] = c_mat[i][j];

                    let dg_ij = omega * femcore::linalg::vecops::dotu(&r_free_list[i], &x_list[j]);
                    let dg_ji =
                        if i == j { dg_ij } else { omega * femcore::linalg::vecops::dotu(&r_free_list[j], &x_list[i]) };
                    let dg_sym = 0.5 * (dg_ij + dg_ji);
                    g_mat[i][j] += dg_sym;
                    g_mat[j][i] = g_mat[i][j];
                }
            }
        }
    }

    let q_dielectric_matrix = if omega > 0.0 {
        let mut q_mat = vec![vec![None; n_term]; n_term];
        for i in 0..n_term {
            for j in 0..n_term {
                let g_abs = g_mat[i][j].abs();
                if g_abs > 1.0e-300 {
                    let q_val = omega * c_mat[i][j].abs() / g_abs;
                    if q_val.is_finite() {
                        q_mat[i][j] = Some(q_val);
                    }
                }
            }
        }
        Some(q_mat)
    } else {
        None
    };

    Ok(AcLossResult {
        freq_hz,
        max_tan_delta_eff,
        conductance_matrix_s: Some(g_mat),
        q_dielectric_matrix,
        delta_c_second_order_f,
        corrected_maxwell_matrix_f: c_mat,
        k_cond_loss: Some(k_cond),
    })
}

/// Runs `solve_electrostatics` and evaluates the AC quasi-static conductance matrix `G_ij(omega)`,
/// positive 2nd-order capacitance shift `Delta C_ij^(2)`, and dielectric quality factor matrix `Q_ij`.
pub fn solve_electrostatics_with_ac_loss(
    mesh: &Discretisation,
    spec: &Problem,
    freq_hz: Option<f64>,
) -> Result<(SolveResult, AcLossResult), String> {
    let res = solve_electrostatics(mesh, spec)?;
    let qp_info = enrich_quarter_point_mesh(mesh, spec);
    let work_mesh = qp_info.as_ref().map_or(mesh, |(m, _, _)| m);
    let order = spec.order.max(1);
    let asm = assemble_electrostatics(work_mesh, spec)?;

    let mut terminals: Vec<(String, Vec<i32>)> = Vec::new();
    for bc in &spec.bcs {
        if let BcKind::Terminal(name) = &bc.kind {
            let name = name.trim();
            if let Some((_, labs)) = terminals.iter_mut().find(|(n, _)| n == name) {
                for &l in &bc.labels {
                    if !labs.contains(&l) {
                        labs.push(l);
                    }
                }
            } else {
                terminals.push((name.to_string(), bc.labels.clone()));
            }
        }
    }

    let n = asm.k.nrows();
    if terminals.is_empty() {
        let free = free_dofs(n, &asm.fixed);
        let v0 = res
            .solutions
            .first()
            .and_then(|s| s.fields.first())
            .map(|f| f.values.clone())
            .unwrap_or_else(|| vec![0.0; n]);
        let c0 = res.scalar("capacitance C").unwrap_or_else(|| {
            let kv = asm.k.matvec(&v0);
            femcore::linalg::vecops::dotu(&v0, &kv)
        });
        let ac = compute_ac_loss(work_mesh, spec, freq_hz, &asm.k, &free, &[v0], &[vec![c0]])?;
        return Ok((res, ac));
    }

    let term_dofs: Vec<BTreeSet<usize>> =
        terminals.iter().map(|(_, labs)| boundary_dofs(work_mesh, labs, order).into_iter().collect()).collect();
    let mut bg_labels = Vec::new();
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Ground | BcKind::Potential(_)) {
            for &l in &bc.labels {
                if !bg_labels.contains(&l) {
                    bg_labels.push(l);
                }
            }
        }
    }
    let bg_dofs: BTreeSet<usize> = boundary_dofs(work_mesh, &bg_labels, order).into_iter().collect();
    let mut all_fixed_set = bg_dofs;
    for td in &term_dofs {
        all_fixed_set.extend(td.iter().copied());
    }
    let fixed: Vec<usize> = all_fixed_set.into_iter().collect();
    let free = free_dofs(n, &fixed);

    let v_all: Vec<Vec<f64>> = res
        .solutions
        .iter()
        .map(|s| s.fields.first().map(|f| f.values.clone()).unwrap_or_else(|| vec![0.0; n]))
        .collect();
    let base_c = CapacitanceMatrixData::of(&res)
        .map(|cd| cd.maxwell_matrix_f)
        .unwrap_or_else(|| vec![vec![0.0; v_all.len()]; v_all.len()]);

    let ac = compute_ac_loss(work_mesh, spec, freq_hz, &asm.k, &free, &v_all, &base_c)?;
    Ok((res, ac))
}
