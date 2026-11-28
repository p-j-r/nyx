// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Complex-valued mode analysis of a 2D cross-section: lossy or absorbing media and perfectly
//! matched layers (PML) by complex coordinate stretching.
//!
//! The formulation is that of `ems`, with complex tensors: permittivity `eps = [[xx, xy], [xy, yy]]`,
//! `zz` and the couplings `xz`, `yz` to the propagation axis (a crystal with its axes tilted against the
//! guide, such as sapphire), permeability `diag(mu_x, mu_y, mu_z)`. The curl-curl
//! equation becomes (with `mu_t^-1` rotated by 90 degrees to `diag(1/mu_y, 1/mu_x)` in the terms
//! that couple to the transverse curl)
//!
//!   S = int (1/mu_z) curl w curl w',  Tm = int w . diag(1/mu_y, 1/mu_x) w',
//!   G = int w . diag(1/mu_y, 1/mu_x) grad(phi),  Kp = int grad(phi) . diag(1/mu_y, 1/mu_x) grad(phi'),
//!   Te = int w . eps_t w',  Tp = int eps_zz phi phi',
//!
//! and the pencil  [S - k0^2 Te, 0; 0, 0] x = -b^2 [Tm, G; G^T, Kp - k0^2 Tp] x  is complex
//! symmetric. A PML along x or y with stretch factors `sx`, `sy` (and `sz = 1`) is the material
//! `eps_ij = det(S) eps_ij / (s_i s_j)`, `mu = mu0 diag(sy/sx, sx/sy, sx sy)`.
//!
//! With `xz`, `yz` different from zero the longitudinal field couples to the transverse one in the
//! electric term and the pencil becomes quadratic in the propagation constant,
//! `(A0 + j b A1 + b^2 B) x = 0` with `A1 = k0^2 [[0, -Mg], [Mg^T, 0]]`, `Mg = int w . (xz, yz) N`.
//! It is solved by fixed-point iteration on `b` in `A1`: each step is the linear pencil
//! `(A0 + j b_prev A1) x = -b^2 B x`, and the iteration stops when `b` changes by less than
//! `COUPLED_BETA_TOL` (relative) between two steps.

use crate::ibc::{
    add_ibc_2d_waveguide, add_ibc_2d_waveguide_mixed, pec_edges_and_nodes_2d, pec_edges_and_nodes_2d_mixed, EPS0, MU0,
};
use crate::mixed_ops::{EmPmlTable, MixedCellMaterials, C3};
use crate::modes::C0;
use crate::problem::{pec_edges_and_nodes_2d_with_symmetry, Problem, Sheet2dBc, Study, SurfaceImpedanceBc, Symmetry};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::edge_mixed::{eval_cell_basis_2d, nedelec_nloc, NedEvalScratch, NedMixed};
use femcore::fem::*;
use femcore::fem_mixed::{lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::eigs_complex;
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::optical::{cell_materials, EpsC, MuC};
use femcore::props::{Pml, TensorKind};
use femcore::section::local_component;
use nyx_physics::{
    cell_fields_of, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use std::time::Instant;

/// Relative change of the propagation constant at which the coupled (tilted-crystal) fixed-point iteration stops.
///
/// Each step solves the linear pencil with a tight eigensolver tolerance (1e-11), but the Krylov solution of a
/// 160,000-unknown complex system carries a noise floor of about 1e-9 in `b`, so two successive values never agree
/// better than that: a tolerance of 1e-10 was never met and the loop always ran its full 10 steps (10 factorisations
/// of ~13 s each on the 300TE chip). The iteration is a strong contraction (the coupling `A1` is a small
/// correction), so the value after the first step is already within the noise floor of the fixed point; 1e-8 stops
/// after one or two steps without changing the result beyond ~1e-9.
const COUPLED_BETA_TOL: f64 = 1e-8;

pub struct ModesC {
    pub neff: Vec<C>,
    pub et: Vec<Vec<C>>,
    pub phi: Vec<Vec<C>>,
    /// True left-eigenvector transverse edge DOFs (`A^H y = conj(lambda) B^H y` when `is_non_symmetric`,
    /// or `y = conj(x)` complex-symmetric shortcut when reciprocal).
    pub et_left: Vec<Vec<C>>,
    /// True left-eigenvector longitudinal scalar potential DOFs.
    pub phi_left: Vec<Vec<C>>,
    /// Biorthogonal inner product `y^H B x` on reduced free DOFs (`x^T B x` when complex-symmetric).
    pub biorth_norms: Vec<C>,
    /// Whether the mode pencil is non-symmetric (`eps_{ij} != eps_{ji}` or `mu_{ij} != mu_{ji}`).
    pub is_non_symmetric: bool,
    pub residuals: Vec<f64>,
    pub restarts: usize,
    pub matvecs: usize,
    pub converged: bool,
}

/// Modes with complex effective index nearest `n_target` (complex).
pub fn solve_modes_complex(
    mesh: &Discretisation,
    freq: f64,
    eps: &[EpsC],
    mu: &[MuC],
    n_target: C,
    k: usize,
    order: usize,
) -> Result<ModesC, String> {
    let asm = assemble_modes_complex(mesh, freq, eps, mu, order)?;
    solve_mode_pairs_complex(&asm, n_target, k)
}

/// Stage 1 of the complex mode solver: `A x = lambda B x` (`lambda = -beta^2`) on the free dofs, plus the part of `A`
/// that is linear in `beta` when the permittivity couples the transverse and the longitudinal field (tilted crystal or Voigt gyrotropy).
pub struct ModeAssemblyC {
    pub a: CCsr,
    pub b: CCsr,
    /// `A1` of `(A + j beta A1) x = -beta^2 B x`, present when longitudinal-transverse coupling exists
    pub a1: Option<CCsr>,
    pub free: Vec<usize>,
    pub n_tot: usize,
    pub net: usize,
    pub k0: f64,
    pub is_non_symmetric: bool,
}

pub fn assemble_modes_complex(
    mesh: &Discretisation,
    freq: f64,
    eps: &[EpsC],
    mu: &[MuC],
    order: usize,
) -> Result<ModeAssemblyC, String> {
    assemble_modes_complex_with_bcs(mesh, freq, eps, mu, order, &[], &[], &[])
}

pub fn assemble_modes_complex_with_bcs(
    mesh: &Discretisation,
    freq: f64,
    eps: &[EpsC],
    mu: &[MuC],
    order: usize,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
) -> Result<ModeAssemblyC, String> {
    assemble_modes_complex_with_bcs_and_symmetry(
        mesh,
        freq,
        eps,
        mu,
        order,
        pec_labels,
        ibc_bcs,
        pmc_labels,
        Symmetry::None,
    )
}

pub(crate) fn assemble_modes_complex_with_bcs_and_symmetry(
    mesh: &Discretisation,
    freq: f64,
    eps: &[EpsC],
    mu: &[MuC],
    order: usize,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
    symmetry: Symmetry,
) -> Result<ModeAssemblyC, String> {
    assemble_modes_complex_with_bcs_and_symmetry_and_sheets(
        mesh,
        freq,
        eps,
        mu,
        order,
        pec_labels,
        ibc_bcs,
        &[],
        None,
        pmc_labels,
        symmetry,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_modes_complex_with_bcs_and_symmetry_and_sheets(
    mesh: &Discretisation,
    freq: f64,
    eps: &[EpsC],
    mu: &[MuC],
    order: usize,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    sheet_bcs: &[crate::problem::Sheet2dBc],
    imported_potential: Option<&[f64]>,
    pmc_labels: &[i32],
    symmetry: Symmetry,
) -> Result<ModeAssemblyC, String> {
    if mesh.has_curved() {
        return Err(
            "curved meshes are not supported by simplex 2D complex mode kernels; use mixed/isoparametric path or a straight mesh"
                .into(),
        );
    }
    if mesh.dim != 2 {
        return Err("mode analysis needs a 2D cross-section mesh".into());
    }
    let (nc, nv, ne) = (mesh.n_cells(), mesh.n_points(), mesh.n_edges());
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nphi = mesh.n_dofs(order);
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let n_tot = net + nphi;
    let is_non_symmetric =
        femcore::optical::is_non_symmetric_cells(eps, mu) || sheet_bcs.iter().any(Sheet2dBc::is_non_symmetric);
    let has_off_diag_mu = mu.iter().any(MuC::has_off_diagonal);
    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut a1_coo = CCoo::new(n_tot, n_tot);
    let mut has_coupling = false;
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    for ci in 0..nc {
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let (e, m) = (eps[ci], mu[ci]);
        if !is_non_symmetric && !has_off_diag_mu {
            let (m1, m2) = (c(1.0, 0.0) / m.y, c(1.0, 0.0) / m.x); // rotated mu_t^-1
            let imz = c(1.0, 0.0) / m.z;
            let n2 = nloc_e * nloc_e;
            let (mut wxx, mut wyy, mut wxy, mut wyx, mut se) =
                (vec![0.0; n2], vec![0.0; n2], vec![0.0; n2], vec![0.0; n2], vec![0.0; n2]);
            let (mut gx, mut gy) = (vec![0.0; nloc_e * nl], vec![0.0; nloc_e * nl]);
            let (mut kxx, mut kyy, mut tp) = (vec![0.0; nl * nl], vec![0.0; nl * nl], vec![0.0; nl * nl]);
            let (mut nxm, mut nym) = (vec![0.0; nloc_e * nl], vec![0.0; nloc_e * nl]);
            for (q, l) in quad.lam.iter().enumerate() {
                edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv_, dn) = &tab[q];
                let mut gn = vec![[0.0; 2]; nl];
                for a in 0..nl {
                    for kk in 0..2 {
                        gn[a][kk] = (0..3).map(|b| dn[a * 3 + b] * g[b][kk]).sum();
                    }
                }
                let wq = quad.w[q];
                for s in 0..nloc_e {
                    for t in 0..nloc_e {
                        wxx[s * nloc_e + t] += wq * w[s][0] * w[t][0];
                        wyy[s * nloc_e + t] += wq * w[s][1] * w[t][1];
                        wxy[s * nloc_e + t] += wq * w[s][0] * w[t][1];
                        wyx[s * nloc_e + t] += wq * w[s][1] * w[t][0];
                        se[s * nloc_e + t] += wq * curl[s] * curl[t];
                    }
                    for a in 0..nl {
                        gx[s * nl + a] += wq * w[s][0] * gn[a][0];
                        gy[s * nl + a] += wq * w[s][1] * gn[a][1];
                        nxm[s * nl + a] += wq * w[s][0] * nv_[a];
                        nym[s * nl + a] += wq * w[s][1] * nv_[a];
                    }
                }
                for a in 0..nl {
                    for b in 0..nl {
                        kxx[a * nl + b] += wq * gn[a][0] * gn[b][0];
                        kyy[a * nl + b] += wq * gn[a][1] * gn[b][1];
                        tp[a * nl + b] += wq * nv_[a] * nv_[b];
                    }
                }
            }
            mesh.cell_dofs(ci, order, &mut pd);
            let rows_e = &cd.gdofs[ci];
            let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
            let ar = c(area, 0.0);
            let k2 = c(k0 * k0, 0.0);
            let amat: Vec<C> = (0..n2)
                .map(|i| (imz * se[i] - k2 * (e.xx * wxx[i] + e.yy * wyy[i] + e.xy * (wxy[i] + wyx[i]))) * ar)
                .collect();
            a_coo.push_block(rows_e, rows_e, &amat);
            let tm: Vec<C> = (0..n2).map(|i| (m1 * wxx[i] + m2 * wyy[i]) * ar).collect();
            b_coo.push_block(rows_e, rows_e, &tm);
            let gm: Vec<C> = (0..nloc_e * nl).map(|i| (m1 * gx[i] + m2 * gy[i]) * ar).collect();
            b_coo.push_block(rows_e, &cols_p, &gm);
            let gt: Vec<C> = (0..nl * nloc_e).map(|i| gm[(i % nloc_e) * nl + i / nloc_e]).collect();
            b_coo.push_block(&cols_p, rows_e, &gt);
            let pp: Vec<C> = (0..nl * nl).map(|i| (m1 * kxx[i] + m2 * kyy[i] - k2 * e.zz * tp[i]) * ar).collect();
            b_coo.push_block(&cols_p, &cols_p, &pp);
            if e.xz.norm() > 0.0 || e.yz.norm() > 0.0 {
                has_coupling = true;
                // A1 = k0^2 [[0, -Mg], [Mg^T, 0]]
                let mg: Vec<C> = (0..nloc_e * nl).map(|i| (e.xz * nxm[i] + e.yz * nym[i]) * ar * k2).collect();
                let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
                a1_coo.push_block(rows_e, &cols_p, &neg);
                let mgt: Vec<C> = (0..nl * nloc_e).map(|i| mg[(i % nloc_e) * nl + i / nloc_e]).collect();
                a1_coo.push_block(&cols_p, rows_e, &mgt);
            }
        } else {
            let ep = e.to_matrix();
            let mui = m.inv_matrix();
            let k2 = c(k0 * k0, 0.0);
            let mut amat = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut tm = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut gm = vec![c(0.0, 0.0); nloc_e * nl];
            let mut gt = vec![c(0.0, 0.0); nl * nloc_e];
            let mut pp = vec![c(0.0, 0.0); nl * nl];
            let mut a1_ee = vec![c(0.0, 0.0); nloc_e * nloc_e];
            let mut mg = vec![c(0.0, 0.0); nloc_e * nl];
            let mut mgt = vec![c(0.0, 0.0); nl * nloc_e];
            let mut cell_cpl = false;
            let mut cell_mu_cpl = false;
            let mut gn = vec![[0.0; 2]; nl];
            for (q, l) in quad.lam.iter().enumerate() {
                edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv_, dn) = &tab[q];
                for a in 0..nl {
                    for kk in 0..2 {
                        gn[a][kk] = (0..3).map(|b| dn[a * 3 + b] * g[b][kk]).sum();
                    }
                }
                add_waveguide_complex_point_full(
                    nloc_e,
                    nl,
                    &w,
                    &curl,
                    nv_,
                    &gn,
                    quad.w[q] * area,
                    k2,
                    &ep,
                    &mui,
                    &mut amat,
                    &mut tm,
                    &mut gm,
                    &mut gt,
                    &mut pp,
                    &mut a1_ee,
                    &mut mg,
                    &mut mgt,
                    &mut cell_cpl,
                    &mut cell_mu_cpl,
                );
            }
            mesh.cell_dofs(ci, order, &mut pd);
            let rows_e = &cd.gdofs[ci];
            let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
            a_coo.push_block(rows_e, rows_e, &amat);
            b_coo.push_block(rows_e, rows_e, &tm);
            b_coo.push_block(rows_e, &cols_p, &gm);
            b_coo.push_block(&cols_p, rows_e, &gt);
            b_coo.push_block(&cols_p, &cols_p, &pp);
            if cell_cpl || cell_mu_cpl {
                has_coupling = true;
                if cell_mu_cpl {
                    a1_coo.push_block(rows_e, rows_e, &a1_ee);
                }
                let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
                a1_coo.push_block(rows_e, &cols_p, &neg);
                a1_coo.push_block(&cols_p, rows_e, &mgt);
            }
        }
    }
    add_ibc_2d_waveguide(mesh, order, freq, pec_labels, ibc_bcs, &mut a_coo, &mut b_coo)?;
    let sheet_cpl = crate::ibc::add_sheet2d_2d_waveguide(
        mesh,
        order,
        freq,
        pec_labels,
        sheet_bcs,
        imported_potential,
        &mut a_coo,
        &mut b_coo,
        &mut a1_coo,
    )?;
    has_coupling |= sheet_cpl;
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();
    let mut non_pec_labels: Vec<i32> = pmc_labels.to_vec();
    for sb in sheet_bcs {
        non_pec_labels.extend(sb.labels.iter().copied());
    }
    // PEC walls
    let (be, bn) =
        if pec_labels.is_empty() && ibc_bcs.is_empty() && non_pec_labels.is_empty() && symmetry == Symmetry::None {
            // walls are the exterior facets only: facets on interfaces between domains are not conductors
            let exterior = mesh.exterior_facets();
            let mut be: Vec<usize> = (0..mesh.n_facets())
                .filter(|&f| exterior[f])
                .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
                .collect();
            be.sort_unstable();
            be.dedup();
            let mut bn: Vec<usize> = (0..mesh.n_facets())
                .filter(|&f| exterior[f])
                .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
                .collect();
            bn.sort_unstable();
            bn.dedup();
            (be, bn)
        } else if symmetry == Symmetry::None {
            pec_edges_and_nodes_2d(mesh, pec_labels, ibc_bcs, &non_pec_labels)
        } else {
            pec_edges_and_nodes_2d_with_symmetry(mesh, pec_labels, ibc_bcs, &non_pec_labels, symmetry)?
        };
    let mut fixed: Vec<usize> = be.clone();
    if order == 2 {
        fixed.extend(be.iter().map(|&e| ne + e));
    }
    fixed.extend(bn.iter().map(|&v| net + v));
    if order == 2 {
        fixed.extend(be.iter().map(|&e| net + nv + e));
    }
    let free = free_dofs(n_tot, &fixed);
    // scale the phi block by 1/k0
    let scale: Vec<f64> = (0..n_tot).map(|i| if i < net { 1.0 } else { 1.0 / k0 }).collect();
    for r in 0..b_full.nrows() {
        for p in b_full.indptr()[r]..b_full.indptr()[r + 1] {
            let s = scale[r] * scale[b_full.indices()[p] as usize];
            b_full.values_mut()[p] *= s;
        }
    }
    let af = a_full.submatrix(&free, &free);
    let bf = b_full.submatrix(&free, &free);
    let a1 = if has_coupling {
        let mut a1 = a1_coo.to_csr();
        for r in 0..a1.nrows() {
            for p in a1.indptr()[r]..a1.indptr()[r + 1] {
                let s = scale[r] * scale[a1.indices()[p] as usize];
                a1.values_mut()[p] *= s;
            }
        }
        Some(a1.submatrix(&free, &free))
    } else {
        None
    };
    Ok(ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0, is_non_symmetric })
}

/// Stage 2: the `k` eigenpairs nearest the complex effective index `n_target`, as effective indices, edge coefficients
/// and `phi = Ez / (j beta)`, sorted by closeness to it. Pairs the Krylov process could not converge are dropped.
pub fn solve_mode_pairs_complex(asm: &ModeAssemblyC, n_target: C, k: usize) -> Result<ModesC, String> {
    let (af, bf, free, n_tot, net, k0) = (&asm.a, &asm.b, &asm.free, asm.n_tot, asm.net, asm.k0);
    let sigma = -(n_target * n_target) * (k0 * k0);
    let sign_re = if n_target.re < 0.0 { -1.0 } else { 1.0 };
    let mut pairs = if asm.is_non_symmetric {
        if let Some(a1f) = &asm.a1 {
            let beta_init = n_target * k0;
            let cmat0 = af.add_scaled(a1f, c(0.0, 1.0) * beta_init);
            eigs_complex(&cmat0, bf, sigma, k, 1e-10, 400)?
        } else {
            eigs_complex(af, bf, sigma, k, 1e-10, 400)?
        }
    } else {
        eigs_complex(af, bf, sigma, k, 1e-10, 400)?
    };
    let mut extra_restarts = 0usize;
    let mut extra_matvecs = 0usize;
    if let Some(a1f) = &asm.a1 {
        let beta_tol = if asm.is_non_symmetric { 1e-11 } else { COUPLED_BETA_TOL };
        // refine every mode by the fixed-point iteration on b in the linear term
        for i in 0..pairs.values.len() {
            let pick = |z: C| -> C {
                // forward-propagating (matching sign of n_target.re) or, for an evanescent mode, decaying root
                if z.re.abs() < 1e-9 * z.norm() {
                    if z.im * sign_re > 0.0 {
                        -z
                    } else {
                        z
                    }
                } else if z.re * sign_re < 0.0 {
                    -z
                } else {
                    z
                }
            };
            let mut beta = pick((-pairs.values[i]).sqrt());
            if beta.re.abs() < 0.1 * beta.im.abs() {
                // an evanescent mode of the uncoupled problem: left as it is (no field is carried along the guide)
                continue;
            }
            let mut last = None;
            // fixed-point iteration beta -> beta(A1 at beta) with Aitken acceleration every third step
            let mut hist: Vec<C> = vec![];
            for it in 0..10 {
                femcore::progress::check()?;
                let cmat = af.add_scaled(a1f, c(0.0, 1.0) * beta);
                let one = eigs_complex(&cmat, bf, -(beta * beta), 1, 1e-11, 400)?;
                extra_restarts += one.restarts;
                extra_matvecs += one.matvecs;
                let nb = pick((-one.values[0]).sqrt());
                let change = (nb - beta).norm() / beta.norm();
                femcore::progress::note(format!("tilted-crystal coupling, optical mode {}: refinement step {}, relative change of the propagation constant {change:.1e} (stop below {beta_tol:.0e})", i + 1, it + 1));
                hist.push(beta);
                beta = nb;
                last = Some(one);
                if change < beta_tol {
                    break;
                }
                if hist.len() >= 2 && hist.len() % 3 == 2 {
                    let (b0, b1, b2) = (hist[hist.len() - 2], hist[hist.len() - 1], beta);
                    let den = b2 - b1 * 2.0 + b0;
                    if den.norm() > 1e-300 {
                        let acc = b2 - (b2 - b1) * (b2 - b1) / den;
                        if (acc - b2).norm() / b2.norm() < 1e-2 {
                            beta = acc;
                        }
                    }
                }
            }
            let one = last.unwrap();
            pairs.values[i] = -(beta * beta);
            pairs.vectors[i] = one.vectors[0].clone();
            if let Some(r) = pairs.residuals.get_mut(i) {
                *r = one.residuals.first().copied().unwrap_or(0.0);
            }
        }
    }
    let mut neff = vec![];
    let (mut et_all, mut phi_all): (Vec<Vec<C>>, Vec<Vec<C>>) = (vec![], vec![]);
    let (mut et_left_all, mut phi_left_all): (Vec<Vec<C>>, Vec<Vec<C>>) = (vec![], vec![]);
    let mut biorth_all: Vec<C> = Vec::new();
    let mut res_all: Vec<f64> = Vec::new();
    let n_wanted = pairs.values.len();
    let frob = |m: &CCsr| m.values().iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    let norm_a0 = frob(af);
    let norm_a1 = asm.a1.as_ref().map(frob).unwrap_or(0.0);
    let norm_b = frob(bf);
    for (i, &lam) in pairs.values.iter().enumerate() {
        let beta2 = -lam;
        let mut nf = beta2.sqrt() / k0;
        if nf.re * sign_re < 0.0 {
            nf = -nf; // choose the root matching sign_re (forward by default, backward if n_target.re < 0)
        }
        // Check true generalized residual of original pencil (af + j beta a1f) u = lam B u
        let u = &pairs.vectors[i];
        let beta = nf * k0;
        let ax = if let Some(a1f) = &asm.a1 {
            let a0u = af.matvec(u);
            let a1u = a1f.matvec(u);
            let jb = c(0.0, 1.0) * beta;
            a0u.into_iter().zip(a1u).map(|(x, y)| x + jb * y).collect::<Vec<C>>()
        } else {
            af.matvec(u)
        };
        let bx = bf.matvec(u);
        let res_num: f64 = ax.iter().zip(&bx).map(|(&a, &b)| (a - lam * b).norm_sqr()).sum::<f64>().sqrt();
        let norm_u: f64 = u.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt().max(1e-300);
        // Normwise backward error (Rigal-Gaches) for `(A0 + j beta A1) u = lam B u`: immune to
        // nanophotonic vs RF block-scaling skew while still rejecting unconverged spurious vectors.
        let res_true = res_num / ((norm_a0 + beta.norm() * norm_a1 + lam.norm() * norm_b) * norm_u + 1e-300);
        if !res_true.is_finite() || res_true > 1e-6 {
            continue;
        }
        let (u_left, biorth_ip) = if asm.is_non_symmetric {
            let cmat_opt = asm.a1.as_ref().map(|a1f| af.add_scaled(a1f, c(0.0, 1.0) * beta));
            let a_op = cmat_opt.as_ref().unwrap_or(af);
            let mut y = femcore::linalg::sens::left_eigenvector(a_op, Some(bf), lam, u, 3)
                .unwrap_or_else(|_| u.iter().map(|z| z.conj()).collect());
            // Biorthogonal inner product <y, B u> = y^H B u; normalize y so y^H B u = 1 when non-zero.
            let ip: C = y.iter().zip(&bx).map(|(&yi, &bxi)| yi.conj() * bxi).sum();
            if ip.norm() > 1e-30 && ip.re.is_finite() && ip.im.is_finite() {
                let scale_y = c(1.0, 0.0) / ip.conj();
                for yi in &mut y {
                    *yi *= scale_y;
                }
                (y, c(1.0, 0.0))
            } else {
                (y, ip)
            }
        } else {
            // Zero-overhead complex-symmetric shortcut y = conj(u), y^H B u = u^T B u
            let y: Vec<C> = u.iter().map(|z| z.conj()).collect();
            let ip: C = u.iter().zip(&bx).map(|(&ui, &bxi)| ui * bxi).sum();
            (y, ip)
        };
        neff.push(nf);
        res_all.push(res_true);
        biorth_all.push(biorth_ip);
        let mut x = vec![c(0.0, 0.0); n_tot];
        let mut x_left = vec![c(0.0, 0.0); n_tot];
        for (a, &dof) in free.iter().enumerate() {
            x[dof] = u[a];
            x_left[dof] = u_left[a];
        }
        et_all.push(x[..net].to_vec());
        phi_all.push(x[net..].iter().map(|v| *v / k0).collect::<Vec<C>>());
        et_left_all.push(x_left[..net].to_vec());
        phi_left_all.push(x_left[net..].iter().map(|v| *v / k0).collect::<Vec<C>>());
    }
    let _ = n_wanted;
    // sort deterministically by closeness to the target (tie-break by -Re(neff))
    let mut idx: Vec<usize> = (0..neff.len()).collect();
    idx.sort_by(|&a, &b| {
        let da = (neff[a] - n_target).norm();
        let db = (neff[b] - n_target).norm();
        da.total_cmp(&db).then_with(|| neff[b].re.total_cmp(&neff[a].re))
    });
    Ok(ModesC {
        neff: idx.iter().map(|&i| neff[i]).collect(),
        et: idx.iter().map(|&i| et_all[i].clone()).collect(),
        phi: idx.iter().map(|&i| phi_all[i].clone()).collect(),
        et_left: idx.iter().map(|&i| et_left_all[i].clone()).collect(),
        phi_left: idx.iter().map(|&i| phi_left_all[i].clone()).collect(),
        biorth_norms: idx.iter().map(|&i| biorth_all[i]).collect(),
        is_non_symmetric: asm.is_non_symmetric,
        residuals: idx.iter().map(|&i| res_all[i]).collect(),
        restarts: pairs.restarts + extra_restarts,
        matvecs: pairs.matvecs + extra_matvecs,
        converged: pairs.converged,
    })
}

/// Nodal magnitude fields of a complex mode, scaled to max |E| = 1, under wafer names: `|E|, |Exz|, |Ey|, Re Ex,
/// Re Ez, Re Ey` (real parts after rotating the phase so that the largest component is real).
pub fn mode_fields_complex(mesh: &Discretisation, modes: &ModesC, j: usize, freq: f64, order: usize) -> Vec<Field> {
    mode_fields_complex_cells(mesh, modes, j, freq, order).0
}

/// Stage 3: `mode_fields_complex` and, next to them, the complex components `Ex`, `Ez`, `Ey` at the local dofs of every
/// cell, in the phase of the plot fields and with their scale (the node average of the real parts is `Re E..`).
pub fn mode_fields_complex_cells(
    mesh: &Discretisation,
    modes: &ModesC,
    j: usize,
    freq: f64,
    order: usize,
) -> (Vec<Field>, CellFields) {
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    let (nc, nloc) = (mesh.n_cells(), mesh.nloc(order));
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let lam = local_bary(2, order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut ex = vec![c(0.0, 0.0); nc * nloc];
    let mut ey = vec![c(0.0, 0.0); nc * nloc];
    let mut ez = vec![c(0.0, 0.0); nc * nloc];
    let beta = modes.neff[j] * k0;
    let jb = c(0.0, 1.0) * beta;
    let mut pd = Vec::new();
    for ci in 0..nc {
        let g = tri_g(&geo, ci);
        mesh.cell_dofs(ci, order, &mut pd);
        for (a, l) in lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (mut x, mut y) = (c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let cf = modes.et[j][cd.gdofs[ci][s]];
                x += cf * w[s][0];
                y += cf * w[s][1];
            }
            ex[ci * nloc + a] = x;
            ey[ci * nloc + a] = y;
            ez[ci * nloc + a] = jb * modes.phi[j][pd[a]];
        }
    }
    let big =
        ex.iter()
            .chain(ey.iter())
            .chain(ez.iter())
            .cloned()
            .fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let n = mesh.n_dofs(order);
    let avg = |v: &[f64]| node_average(mesh, v, 1, order);
    let mag = |a: &[C]| -> Vec<f64> { a.iter().map(|z| z.norm()).collect() };
    let re = |a: &[C]| -> Vec<f64> { a.iter().map(|z| (*z * rot).re).collect() };
    let (axm, aym, azm) = (avg(&mag(&ex)), avg(&mag(&ey)), avg(&mag(&ez)));
    let etn: Vec<f64> = (0..nc * nloc).map(|i| (ex[i].norm_sqr() + ey[i].norm_sqr()).sqrt()).collect();
    let etn = avg(&etn);
    let _ = (axm, aym);
    let ezn = azm;
    let tot: Vec<f64> = (0..n).map(|i| etn[i].hypot(ezn[i])).collect();
    let s = 1.0 / tot.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let f = |name: &str, v: Vec<f64>| Field::scalar(name, "a.u.", v.iter().map(|x| x * s).collect());
    let (lat, nor, prop) = (local_component(0), local_component(1), local_component(2));
    let signed = |c: &[C], sign: f64| -> Vec<f64> { avg(&re(c)).iter().map(|x| x * sign).collect() };
    let mut cells = CellFields::new(mesh, order);
    for (comp, v, sign) in [(lat, &ex, lat.1), (nor, &ey, nor.1), (prop, &ez, prop.1)] {
        let rotated: Vec<C> = v.iter().map(|z| *z * rot * s * sign).collect();
        cells.push_complex(
            format!("E{}", comp.0),
            "a.u.",
            rotated.iter().map(|z| z.re).collect(),
            rotated.iter().map(|z| z.im).collect(),
        );
    }
    let fields = vec![
        f("|E|", tot.clone()),
        f("|Exz|", etn),
        f(&format!("|E{}|", prop.0), ezn),
        f(&format!("Re E{}", lat.0), signed(&ex, lat.1)),
        f(&format!("Re E{}", nor.0), signed(&ey, nor.1)),
        f(&format!("Re E{}", prop.0), signed(&ez, prop.1)),
    ];
    (fields, cells)
}

/// The modes of a `Modes` study in the complex solver: absorbing layers, lossy or complex materials, a crystal tilted against the guide.
pub fn solve_waveguide_complex(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    if mesh.has_curved() {
        return Err(
            "curved meshes are not supported by simplex 2D complex mode kernels; use mixed/isoparametric path or a straight mesh"
                .into(),
        );
    }
    let t0 = Instant::now();
    let Study::Modes { freq, n_target, evanescent, n_modes } = problem.study.clone() else {
        return Err("mode analysis needs a Modes study".into());
    };
    if freq <= 0.0 || n_modes < 1 {
        return Err("frequency must be positive and at least one mode is needed".into());
    }
    let (eps, mu) = cell_materials(mesh, &problem.materials)?;
    let pml_table = EmPmlTable::from_simplex(mesh, &problem.materials, freq, true, false);
    let needs_quad_pml = pml_table.has_any()
        && problem.materials.pml.values().any(|p| {
            p.is_active()
                && (p.coord != femcore::props::PmlCoord::Cartesian
                    || p.bounds.is_some()
                    || p.auto_scale
                    || matches!(p.stabilization, femcore::props::PmlStabilization::Mpml { ratio } if ratio != 0.0)
                    || p.profile != femcore::props::PmlProfile::KappaPolynomial { kappa_max: 1.0 })
        });
    let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
    let has_pml = pml_table.has_any();
    let asm = if needs_quad_pml || problem.has_cell_material_perturbation() {
        assemble_modes_complex_quad_pml(mesh, problem, freq, &pml_table)?
    } else if problem.pec.is_empty()
        && problem.surface_impedance.is_empty()
        && problem.sheet_2d.is_empty()
        && problem.pmc.is_empty()
        && problem.symmetry == Symmetry::None
    {
        let (ep_loc, mu_loc) = femcore::optical::cell_materials(mesh, &problem.materials)?;
        assemble_modes_complex(mesh, freq, &ep_loc, &mu_loc, problem.order)?
    } else {
        assemble_modes_complex_with_bcs_and_symmetry_and_sheets(
            mesh,
            freq,
            &eps,
            &mu,
            problem.order,
            &problem.pec,
            &problem.surface_impedance,
            &problem.sheet_2d,
            problem.imported_potential.as_deref(),
            &problem.pmc,
            problem.symmetry,
        )?
    };
    let n_ask = if has_pml && n_modes <= 4 && asm.a.nrows() <= 4_000 {
        (n_modes + 4).min(asm.a.nrows().saturating_sub(1)).max(n_modes)
    } else {
        n_modes
    };
    let mut modes = solve_mode_pairs_complex(&asm, target, n_ask)?;
    let mut solutions = vec![];
    let mut scalars = vec![];
    let mut mode_metrics = compute_waveguide_2d_metrics(mesh, problem, freq, &modes.neff, &modes.et, &modes.phi)?;
    if has_pml && !modes.neff.is_empty() {
        let mut keep: Vec<usize> =
            (0..modes.neff.len()).filter(|&i| mode_metrics[i].layer_fraction < 0.25).take(n_modes).collect();
        if keep.is_empty() && n_ask < n_modes + 4 && asm.a.nrows() > n_modes + 1 {
            let n_retry = (n_modes + 4).min(asm.a.nrows().saturating_sub(1)).max(n_modes);
            modes = solve_mode_pairs_complex(&asm, target, n_retry)?;
            mode_metrics = compute_waveguide_2d_metrics(mesh, problem, freq, &modes.neff, &modes.et, &modes.phi)?;
            keep = (0..modes.neff.len()).filter(|&i| mode_metrics[i].layer_fraction < 0.25).take(n_modes).collect();
        }
        if !keep.is_empty() {
            modes.neff = keep.iter().map(|&i| modes.neff[i]).collect();
            modes.et = keep.iter().map(|&i| modes.et[i].clone()).collect();
            modes.phi = keep.iter().map(|&i| modes.phi[i].clone()).collect();
            modes.residuals = keep.iter().map(|&i| modes.residuals[i]).collect();
            mode_metrics = keep
                .into_iter()
                .map(|i| {
                    std::mem::replace(
                        &mut mode_metrics[i],
                        WaveguideModeMetrics {
                            vg: 0.0,
                            comp_shares: [0.0; 3],
                            domain_shares: Vec::new(),
                            layer_fraction: 1.0,
                        },
                    )
                })
                .collect();
        } else {
            modes.neff.truncate(n_modes);
            modes.et.truncate(n_modes);
            modes.phi.truncate(n_modes);
            modes.residuals.truncate(n_modes);
            mode_metrics.truncate(n_modes);
        }
    }
    let q_opt_modes = if problem.has_loss() {
        Some(crate::heat_source::compute_2d_waveguide_heat_source_simplex(
            mesh,
            problem,
            freq,
            &modes.neff,
            &modes.et,
            &modes.phi,
        )?)
    } else {
        None
    };
    for j in 0..modes.neff.len() {
        let z = modes.neff[j];
        let txt = format!("{:.6} {:+.3e}i", z.re, z.im);
        let (mut fields, cell_fields) = mode_fields_complex_cells(mesh, &modes, j, freq, problem.order);
        if let Some(ref qs) = q_opt_modes {
            fields.push(Field::scalar("Q_opt", "W/m^3", qs[j].0.clone()));
        }
        solutions.push(Solution {
            label: format!("mode {}: neff = {txt}", j + 1),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            eigenvalue: Some([z.re, z.im]),
            ..Default::default()
        });
        scalars.push((format!("neff {}", j + 1), Scalar::Text(txt)));
        let q = if z.im.abs() > 0.0 && z.re > 0.0 { z.re / (2.0 * z.im.abs()) } else { f64::INFINITY };
        if q.is_finite() {
            scalars.push((format!("Q {} = Re/(2|Im|) of neff", j + 1), Scalar::Num(q)));
        }
        mode_metrics[j].push_scalars(j + 1, &mut scalars);
    }
    let mut log = vec!["complex mode analysis (PML or complex materials)".to_string()];
    if !modes.converged {
        log.push("warning: the eigensolver did not converge every mode to the requested tolerance".into());
    }
    scalars.push(crate::optical::mode_frame(problem).scalar());
    let mut diag = crate::diag::EmsDiag::new("shift_invert_complex");
    diag.record_mode_residuals("shift_invert_complex", modes.residuals.clone(), modes.converged, Some(modes.restarts));
    let (diagnostics, error_budget) = diag.finish();
    let result = SolveResult {
        physics: "ems".into(),
        study: "modes".into(),
        order: problem.order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: mesh.n_dofs(problem.order),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    };
    Ok(result)
}

/// Poynting group velocity / group index and non-PML electric energy shares for a 2D waveguide mode.
pub(crate) struct WaveguideModeMetrics {
    pub vg: f64,
    pub comp_shares: [f64; 3],
    pub domain_shares: Vec<(i32, f64)>,
    pub layer_fraction: f64,
}

impl WaveguideModeMetrics {
    pub fn push_scalars(&self, mode_num: usize, scalars: &mut Vec<(String, Scalar)>) {
        scalars.push((format!("v_g {mode_num} (m/s)"), Scalar::Num(self.vg)));
        if self.vg.abs() > 1e-9 * C0 {
            scalars.push((format!("n_g {mode_num}"), Scalar::Num(C0 / self.vg)));
        }
        scalars.push((format!("Ex share {mode_num}"), Scalar::Num(self.comp_shares[0])));
        scalars.push((format!("Ey share {mode_num}"), Scalar::Num(self.comp_shares[1])));
        scalars.push((format!("Ez share {mode_num}"), Scalar::Num(self.comp_shares[2])));
        for &(label, share) in &self.domain_shares {
            scalars.push((format!("share in domain {label} {mode_num}"), Scalar::Num(share)));
        }
        scalars.push((format!("layer fraction {mode_num}"), Scalar::Num(self.layer_fraction)));
    }
}

/// Evaluates `P_y`, dispersive `U_EM`, `v_g`, `n_g`, component shares (`Ex`, `Ey`, `Ez`), domain shares,
/// and PML `layer fraction` at `quad_rule(2, 4)` using unstretched physical `eps_r`, `d(omega eps_r)/domega`, and `mu_r`.
pub(crate) fn compute_waveguide_2d_metrics(
    mesh: &Discretisation,
    problem: &Problem,
    freq: f64,
    neff: &[C],
    et: &[Vec<C>],
    phi: &[Vec<C>],
) -> Result<Vec<WaveguideModeMetrics>, String> {
    let n_modes = neff.len();
    if n_modes == 0 {
        return Ok(Vec::new());
    }
    let omega = (2.0 * std::f64::consts::PI * freq).max(1e-300);
    let k0 = omega / C0;
    let eps_phys_re = problem.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?;
    let d_omega_eps = problem.cell_d_omega_eps(mesh)?;
    let mur_phys = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let pml_cells: Vec<bool> = (0..mesh.n_cells())
        .map(|ci| problem.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(Pml::is_active))
        .collect();
    let any_pml = pml_cells.iter().any(|&p| p);
    let mut domains: Vec<i32> =
        (0..mesh.n_cells()).filter(|&ci| !pml_cells[ci]).map(|ci| mesh.cell_labels[ci]).collect();
    domains.sort_unstable();
    domains.dedup();
    let cell_dom_idx: Vec<usize> = (0..mesh.n_cells())
        .map(|ci| {
            if pml_cells[ci] {
                usize::MAX
            } else {
                domains.binary_search(&mesh.cell_labels[ci]).unwrap_or(usize::MAX)
            }
        })
        .collect();

    let geo = geometry(mesh);
    let cd = cell_data(mesh, problem.order);
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, problem.order)).collect();
    let nl = mesh.nloc(problem.order);
    let nloc_e = if problem.order == 1 { 3 } else { 8 };

    let mut p_y = vec![0.0_f64; n_modes];
    let mut u_em = vec![0.0_f64; n_modes];
    let mut comp_e = vec![[0.0_f64; 3]; n_modes];
    let mut dom_e = vec![vec![0.0_f64; domains.len()]; n_modes];
    let mut u_e_phys = vec![0.0_f64; n_modes];
    let mut u_e_pml = vec![0.0_f64; n_modes];

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];

    for ci in 0..mesh.n_cells() {
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let er = &eps_phys_re[ci];
        let de = &d_omega_eps[ci];
        let mur = mur_phys[ci];
        let inv_om_mu = 1.0 / (omega * MU0 * mur);
        let d_idx = cell_dom_idx[ci];
        mesh.cell_dofs(ci, problem.order, &mut pd);
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], problem.order, &mut w, &mut curl);
            let (nv_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let wq_area = quad.w[q] * area;
            for j in 0..n_modes {
                let beta = neff[j] * k0;
                let j_beta = c(0.0, 1.0) * beta;
                let (mut ex, mut ez, mut curl_2d) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                for s in 0..nloc_e {
                    let cf = et[j][cd.gdofs[ci][s]];
                    ex += cf * w[s][0];
                    ez += cf * w[s][1];
                    curl_2d += cf * curl[s];
                }
                let (mut phi_q, mut dphi_dx, mut dphi_dz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                for a in 0..nl {
                    let pf = phi[j][pd[a]];
                    phi_q += pf * nv_[a];
                    dphi_dx += pf * gn[a][0];
                    dphi_dz += pf * gn[a][1];
                }
                let ey = j_beta * phi_q;
                let e_w = [ex, ey, ez];
                let mut e_eps_e = 0.0_f64;
                for a in 0..3 {
                    for b in 0..3 {
                        e_eps_e += er[a][b] * (e_w[a].conj() * e_w[b]).re;
                    }
                }
                if pml_cells[ci] {
                    u_e_pml[j] += e_eps_e * wq_area;
                } else {
                    u_e_phys[j] += e_eps_e * wq_area;
                    for k in 0..3 {
                        comp_e[j][k] += er[k][k].abs() * e_w[k].norm_sqr() * wq_area;
                    }
                    if d_idx != usize::MAX {
                        dom_e[j][d_idx] += e_eps_e * wq_area;
                    }
                    // Wafer right-handed (x, y, z) magnetic field H = (j / (omega mu0 mu_r)) curl E:
                    // (curl E)_x = -j beta (e_z + dphi/dz)  => H_x = +(beta / (omega mu)) (e_z + dphi/dz)
                    // (curl E)_y = -(de_z/dx - de_x/dz)     => H_y = -(j / (omega mu)) curl_2d
                    // (curl E)_z = +j beta (e_x + dphi/dx)  => H_z = -(beta / (omega mu)) (e_x + dphi/dx)
                    let hx = beta * (ez + dphi_dz) * inv_om_mu;
                    let hy = -c(0.0, 1.0) * curl_2d * inv_om_mu;
                    let hz = -beta * (ex + dphi_dx) * inv_om_mu;
                    p_y[j] += 0.5 * (ez.conj() * hx - ex.conj() * hz).re * wq_area;
                    let mut e_de_e = 0.0_f64;
                    for a in 0..3 {
                        for b in 0..3 {
                            e_de_e += de[a][b] * (e_w[a].conj() * e_w[b]).re;
                        }
                    }
                    let h_mu_h = mur * (hx.norm_sqr() + hy.norm_sqr() + hz.norm_sqr());
                    u_em[j] += 0.25 * (EPS0 * e_de_e + MU0 * h_mu_h) * wq_area;
                }
            }
        }
    }

    let mut out = Vec::with_capacity(n_modes);
    for j in 0..n_modes {
        let vg_raw = if u_em[j] > 0.0 { p_y[j] / u_em[j] } else { 0.0 };
        let vg = if vg_raw.abs() <= 1e-9 * C0 { 0.0 } else { vg_raw };
        let comp_tot = comp_e[j][0] + comp_e[j][1] + comp_e[j][2];
        let comp_shares = if comp_tot > 0.0 {
            [comp_e[j][0] / comp_tot, comp_e[j][1] / comp_tot, comp_e[j][2] / comp_tot]
        } else {
            [0.0; 3]
        };
        let dom_tot: f64 = dom_e[j].iter().sum();
        let domain_shares = domains
            .iter()
            .enumerate()
            .map(|(d, &label)| (label, if dom_tot > 0.0 { dom_e[j][d] / dom_tot } else { 0.0 }))
            .collect();
        let u_all = u_e_phys[j] + u_e_pml[j];
        let layer_fraction = if any_pml && u_all > 0.0 { (u_e_pml[j] / u_all).clamp(0.0, 1.0) } else { 0.0 };
        out.push(WaveguideModeMetrics { vg, comp_shares, domain_shares, layer_fraction });
    }
    Ok(out)
}

// ------------------------------------------------------------------ Quadrature-point PML & Native 2D MixedDiscretisation (`Tri + Quad`)

#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_waveguide_complex_point_full(
    nloc_e: usize,
    nl: usize,
    w: &[[f64; 2]],
    curl: &[f64],
    nv: &[f64],
    gn: &[[f64; 2]],
    wq: f64,
    k2: C,
    ep: &C3,
    mui: &C3,
    amat: &mut [C],
    tm: &mut [C],
    gm: &mut [C],
    gt: &mut [C],
    pp: &mut [C],
    a1_ee: &mut [C],
    mg: &mut [C],
    mgt: &mut [C],
    has_coupling: &mut bool,
    has_mu_coupling: &mut bool,
) {
    // Rotated transverse inverse permeability M_mu = R_90^T [[nu00, nu01], [nu10, nu11]] R_90 = [[nu11, -nu10], [-nu01, nu00]]
    let m00 = mui[1][1];
    let m01 = -mui[1][0];
    let m10 = -mui[0][1];
    let m11 = mui[0][0];
    let imz = mui[2][2];
    let cpl = ep[0][2].norm_sqr() + ep[1][2].norm_sqr() + ep[2][0].norm_sqr() + ep[2][1].norm_sqr() > 0.0;
    if cpl {
        *has_coupling = true;
    }
    let mu_cpl = mui[0][2].norm_sqr() + mui[1][2].norm_sqr() + mui[2][0].norm_sqr() + mui[2][1].norm_sqr() > 0.0;
    if mu_cpl {
        *has_mu_coupling = true;
    }

    for s in 0..nloc_e {
        let ws = w[s];
        // Row s (test function ws): ws^T M_mu = [ws0*m00 + ws1*m10, ws0*m01 + ws1*m11]
        let mw_s = [m00 * ws[0] + m10 * ws[1], m01 * ws[0] + m11 * ws[1]];
        // Row s (test function ws): ws^T eps_t = [ws0*ep00 + ws1*ep10, ws0*ep01 + ws1*ep11]
        let ew_s = [ep[0][0] * ws[0] + ep[1][0] * ws[1], ep[0][1] * ws[0] + ep[1][1] * ws[1]];
        for t in 0..nloc_e {
            let wt = w[t];
            let idx = s * nloc_e + t;
            amat[idx] += (imz * (curl[s] * curl[t]) - k2 * (ew_s[0] * wt[0] + ew_s[1] * wt[1])) * wq;
            tm[idx] += (mw_s[0] * wt[0] + mw_s[1] * wt[1]) * wq;
            if mu_cpl {
                let t1 = curl[s] * (mui[2][0] * (-wt[1]) + mui[2][1] * wt[0]);
                let t2 = ((-ws[1]) * mui[0][2] + ws[0] * mui[1][2]) * curl[t];
                a1_ee[idx] += (t1 - t2) * wq;
            }
        }
        for a in 0..nl {
            let ga = gn[a];
            gm[s * nl + a] += (mw_s[0] * ga[0] + mw_s[1] * ga[1]) * wq;
            if cpl {
                mg[s * nl + a] += k2 * (ep[0][2] * ws[0] + ep[1][2] * ws[1]) * (nv[a] * wq);
            }
            if mu_cpl {
                // Note: mg is pushed into A1 with a minus sign (-mg), so subtract the +j beta curl[s]*(nu_z . R90 grad phi_a) term
                mg[s * nl + a] -= curl[s] * (mui[2][0] * (-ga[1]) + mui[2][1] * ga[0]) * wq;
            }
        }
    }
    for a in 0..nl {
        let ga = gn[a];
        let mga = [m00 * ga[0] + m10 * ga[1], m01 * ga[0] + m11 * ga[1]];
        for t in 0..nloc_e {
            let wt = w[t];
            gt[a * nloc_e + t] += (mga[0] * wt[0] + mga[1] * wt[1]) * wq;
            if cpl {
                mgt[a * nloc_e + t] += k2 * (ep[2][0] * wt[0] + ep[2][1] * wt[1]) * (nv[a] * wq);
            }
            if mu_cpl {
                mgt[a * nloc_e + t] -= ((-ga[1]) * mui[0][2] + ga[0] * mui[1][2]) * (curl[t] * wq);
            }
        }
        for b in 0..nl {
            let gb = gn[b];
            pp[a * nl + b] += (mga[0] * gb[0] + mga[1] * gb[1] - k2 * ep[2][2] * (nv[a] * nv[b])) * wq;
        }
    }
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn add_waveguide_complex_point(
    nloc_e: usize,
    nl: usize,
    w: &[[f64; 2]],
    curl: &[f64],
    nv: &[f64],
    gn: &[[f64; 2]],
    wq: f64,
    k2: C,
    ep: &C3,
    mui: &C3,
    amat: &mut [C],
    tm: &mut [C],
    gm: &mut [C],
    gt: &mut [C],
    pp: &mut [C],
    a1_ee: &mut [C],
    mg: &mut [C],
    mgt: &mut [C],
    has_coupling: &mut bool,
    has_mu_coupling: &mut bool,
) {
    add_waveguide_complex_point_full(
        nloc_e,
        nl,
        w,
        curl,
        nv,
        gn,
        wq,
        k2,
        ep,
        mui,
        amat,
        tm,
        gm,
        gt,
        pp,
        a1_ee,
        mg,
        mgt,
        has_coupling,
        has_mu_coupling,
    );
}

/// Quadrature-point curvilinear PML assembly on a 2D `SimplexMesh` (`Discretisation`).
pub(crate) fn assemble_modes_complex_quad_pml(
    mesh: &Discretisation,
    problem: &Problem,
    freq: f64,
    pml: &EmPmlTable,
) -> Result<ModeAssemblyC, String> {
    if mesh.has_curved() {
        return Err(
            "curved meshes are not supported by simplex 2D complex mode kernels; use mixed/isoparametric path or a straight mesh"
                .into(),
        );
    }
    if mesh.dim != 2 {
        return Err("mode analysis needs a 2D cross-section mesh".into());
    }
    let order = problem.order;
    let (nc, nv, ne) = (mesh.n_cells(), mesh.n_points(), mesh.n_edges());
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    let k2 = c(k0 * k0, 0.0);
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let net = cd.net;
    let nphi = mesh.n_dofs(order);
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let quad4 = quad_rule(2, 4);
    let quad6 = quad_rule(2, 6);
    let tab4: Vec<(Vec<f64>, Vec<f64>)> = quad4.lam.iter().map(|l| shape(l, order)).collect();
    let tab6: Vec<(Vec<f64>, Vec<f64>)> = quad6.lam.iter().map(|l| shape(l, order)).collect();

    let mut eps_w = problem.materials.cell_full_complex_permittivity(mesh, Some(1.0))?;
    if problem.has_cell_material_perturbation() {
        let tens = problem.cell_complex_permittivity(mesh)?;
        let eo_pe = problem.domain_eo_pe_global();
        for ci in 0..nc {
            let dom = mesh.cell_labels[ci];
            if problem.thermo_optic_domains.contains_key(&dom) || eo_pe.contains_key(&dom) {
                eps_w[ci] = crate::mixed_ops::apply_thermo_optic_to_c3(eps_w[ci], &tens[ci].0, &tens[ci].1);
            }
        }
    }
    let mui_w = problem.materials.cell_full_complex_permeability_inv(mesh, Some(1.0))?;
    let mut eps_loc = Vec::with_capacity(nc);
    let mut mui_loc = Vec::with_capacity(nc);
    let perm_loc = |t: &C3| -> C3 {
        let axes = [0usize, 2, 1];
        let mut r = [[c(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = t[axes[i]][axes[j]];
            }
        }
        r
    };
    for ci in 0..nc {
        eps_loc.push(perm_loc(&eps_w[ci]));
        mui_loc.push(perm_loc(&mui_w[ci]));
    }
    let is_non_symmetric = problem.materials.is_non_symmetric()
        || crate::mixed_ops::is_c3_slice_non_symmetric(&eps_loc, &mui_loc)
        || problem.sheet_2d.iter().any(Sheet2dBc::is_non_symmetric);

    let n_tot = net + nphi;
    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut a1_coo = CCoo::new(n_tot, n_tot);
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut has_coupling = false;

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];

    for ci in 0..nc {
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let is_pml = pml.is_pml(ci);
        let vs = mesh.cell(ci);
        let pts = [mesh.point(vs[0] as usize), mesh.point(vs[1] as usize), mesh.point(vs[2] as usize)];
        let (qrule, tab) = if is_pml { (&quad6, &tab6) } else { (&quad4, &tab4) };

        let mut amat = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut tm = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut gm = vec![c(0.0, 0.0); nloc_e * nl];
        let mut gt = vec![c(0.0, 0.0); nl * nloc_e];
        let mut pp = vec![c(0.0, 0.0); nl * nl];
        let mut a1_ee = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut mg = vec![c(0.0, 0.0); nloc_e * nl];
        let mut mgt = vec![c(0.0, 0.0); nl * nloc_e];
        let mut cell_cpl = false;
        let mut cell_mu_cpl = false;

        for (q, l) in qrule.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (nv_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let xq = [
                l[0] * pts[0][0] + l[1] * pts[1][0] + l[2] * pts[2][0],
                l[0] * pts[0][1] + l[1] * pts[1][1] + l[2] * pts[2][1],
            ];
            let (ep, mui, _) = pml.eval_2d(ci, eps_loc[ci], mui_loc[ci], xq, false);
            add_waveguide_complex_point(
                nloc_e,
                nl,
                &w,
                &curl,
                nv_,
                &gn,
                qrule.w[q] * area,
                k2,
                &ep,
                &mui,
                &mut amat,
                &mut tm,
                &mut gm,
                &mut gt,
                &mut pp,
                &mut a1_ee,
                &mut mg,
                &mut mgt,
                &mut cell_cpl,
                &mut cell_mu_cpl,
            );
        }

        mesh.cell_dofs(ci, order, &mut pd);
        let rows_e = &cd.gdofs[ci];
        let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
        a_coo.push_block(rows_e, rows_e, &amat);
        b_coo.push_block(rows_e, rows_e, &tm);
        b_coo.push_block(rows_e, &cols_p, &gm);
        b_coo.push_block(&cols_p, rows_e, &gt);
        b_coo.push_block(&cols_p, &cols_p, &pp);
        if cell_cpl || cell_mu_cpl {
            has_coupling = true;
            if cell_mu_cpl {
                a1_coo.push_block(rows_e, rows_e, &a1_ee);
            }
            let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
            a1_coo.push_block(rows_e, &cols_p, &neg);
            a1_coo.push_block(&cols_p, rows_e, &mgt);
        }
    }

    add_ibc_2d_waveguide(mesh, order, freq, &problem.pec, &problem.surface_impedance, &mut a_coo, &mut b_coo)?;
    let sheet_cpl = crate::ibc::add_sheet2d_2d_waveguide(
        mesh,
        order,
        freq,
        &problem.pec,
        &problem.sheet_2d,
        problem.imported_potential.as_deref(),
        &mut a_coo,
        &mut b_coo,
        &mut a1_coo,
    )?;
    has_coupling |= sheet_cpl;
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();
    let mut non_pec_labels: Vec<i32> = problem.pmc.clone();
    for sb in &problem.sheet_2d {
        non_pec_labels.extend(sb.labels.iter().copied());
    }
    let (be, bn) = if problem.pec.is_empty()
        && problem.surface_impedance.is_empty()
        && non_pec_labels.is_empty()
        && problem.symmetry == Symmetry::None
    {
        let exterior = mesh.exterior_facets();
        let mut be: Vec<usize> = (0..mesh.n_facets())
            .filter(|&f| exterior[f])
            .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
            .collect();
        be.sort_unstable();
        be.dedup();
        let mut bn: Vec<usize> = (0..mesh.n_facets())
            .filter(|&f| exterior[f])
            .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
            .collect();
        bn.sort_unstable();
        bn.dedup();
        (be, bn)
    } else if problem.symmetry == Symmetry::None {
        pec_edges_and_nodes_2d(mesh, &problem.pec, &problem.surface_impedance, &non_pec_labels)
    } else {
        pec_edges_and_nodes_2d_with_symmetry(
            mesh,
            &problem.pec,
            &problem.surface_impedance,
            &non_pec_labels,
            problem.symmetry,
        )?
    };

    let mut fixed: Vec<usize> = be.clone();
    if order == 2 {
        fixed.extend(be.iter().map(|&e| ne + e));
    }
    fixed.extend(bn.iter().map(|&v| net + v));
    if order == 2 {
        fixed.extend(be.iter().map(|&e| net + nv + e));
    }
    let free = free_dofs(n_tot, &fixed);
    let scale: Vec<f64> = (0..n_tot).map(|i| if i < net { 1.0 } else { 1.0 / k0 }).collect();
    for r in 0..b_full.nrows() {
        for p in b_full.indptr()[r]..b_full.indptr()[r + 1] {
            let s = scale[r] * scale[b_full.indices()[p] as usize];
            b_full.values_mut()[p] *= s;
        }
    }
    let af = a_full.submatrix(&free, &free);
    let bf = b_full.submatrix(&free, &free);
    let a1 = if has_coupling {
        let mut a1 = a1_coo.to_csr();
        for r in 0..a1.nrows() {
            for p in a1.indptr()[r]..a1.indptr()[r + 1] {
                let s = scale[r] * scale[a1.indices()[p] as usize];
                a1.values_mut()[p] *= s;
            }
        }
        Some(a1.submatrix(&free, &free))
    } else {
        None
    };
    Ok(ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0, is_non_symmetric })
}

/// Solves a 2D complex `Modes` study natively on `Tri + Quad` (`MixedDiscretisation`).
pub fn solve_waveguide_complex_mixed(
    mesh: &Mesh,
    mixed: &MixedDiscretisation,
    problem: &Problem,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let Study::Modes { freq, n_target, evanescent, n_modes } = problem.study.clone() else {
        return Err("mode analysis needs a Modes study".into());
    };
    if freq <= 0.0 || n_modes < 1 {
        return Err("frequency must be positive and at least one mode is needed".into());
    }
    let (asm, ned, phi_map) = assemble_modes_complex_mixed(mixed, problem, freq)?;
    let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
    let has_pml = problem.materials.has_pml();
    let n_ask = if has_pml && n_modes <= 4 && asm.a.nrows() <= 4_000 {
        (n_modes + 4).min(asm.a.nrows().saturating_sub(1)).max(n_modes)
    } else {
        n_modes
    };
    let mut modes = solve_mode_pairs_complex(&asm, target, n_ask)?;

    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, mixed).map_err(|e| e.to_string())?;
    let mut mode_metrics =
        compute_waveguide_2d_metrics_mixed(mixed, &ned, &phi_map, problem, freq, &modes.neff, &modes.et, &modes.phi)?;
    if has_pml && !modes.neff.is_empty() {
        let mut keep: Vec<usize> =
            (0..modes.neff.len()).filter(|&i| mode_metrics[i].layer_fraction < 0.25).take(n_modes).collect();
        if keep.is_empty() && n_ask < n_modes + 4 && asm.a.nrows() > n_modes + 1 {
            let n_retry = (n_modes + 4).min(asm.a.nrows().saturating_sub(1)).max(n_modes);
            modes = solve_mode_pairs_complex(&asm, target, n_retry)?;
            mode_metrics = compute_waveguide_2d_metrics_mixed(
                mixed,
                &ned,
                &phi_map,
                problem,
                freq,
                &modes.neff,
                &modes.et,
                &modes.phi,
            )?;
            keep = (0..modes.neff.len()).filter(|&i| mode_metrics[i].layer_fraction < 0.25).take(n_modes).collect();
        }
        if !keep.is_empty() {
            modes.neff = keep.iter().map(|&i| modes.neff[i]).collect();
            modes.et = keep.iter().map(|&i| modes.et[i].clone()).collect();
            modes.phi = keep.iter().map(|&i| modes.phi[i].clone()).collect();
            modes.residuals = keep.iter().map(|&i| modes.residuals[i]).collect();
            mode_metrics = keep
                .into_iter()
                .map(|i| {
                    std::mem::replace(
                        &mut mode_metrics[i],
                        WaveguideModeMetrics {
                            vg: 0.0,
                            comp_shares: [0.0; 3],
                            domain_shares: Vec::new(),
                            layer_fraction: 1.0,
                        },
                    )
                })
                .collect();
        } else {
            modes.neff.truncate(n_modes);
            modes.et.truncate(n_modes);
            modes.phi.truncate(n_modes);
            modes.residuals.truncate(n_modes);
            mode_metrics.truncate(n_modes);
        }
    }

    let q_opt_modes = if problem.has_loss() {
        let mats_q = MixedCellMaterials::from_problem(mixed, problem)?;
        Some(crate::heat_source::compute_2d_waveguide_heat_source_mixed(
            mixed,
            &ned,
            &phi_map,
            &simplex,
            &proj,
            problem,
            freq,
            &modes.neff,
            &modes.et,
            &modes.phi,
            &mats_q,
        )?)
    } else {
        None
    };
    let mut solutions = Vec::with_capacity(modes.neff.len());
    let mut scalars = Vec::new();
    for j in 0..modes.neff.len() {
        let z = modes.neff[j];
        let txt = format!("{:.6} {:+.3e}i", z.re, z.im);
        let (mut fields, cell_fields) =
            mode_fields_complex_cells_mixed(mixed, &ned, &phi_map, &simplex, &proj, &modes, j, freq, problem.order);
        if let Some(ref qs) = q_opt_modes {
            fields.push(Field::scalar("Q_opt", "W/m^3", qs[j].0.clone()));
        }
        solutions.push(Solution {
            label: format!("mode {}: neff = {txt}", j + 1),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            eigenvalue: Some([z.re, z.im]),
            ..Default::default()
        });
        scalars.push((format!("neff {}", j + 1), Scalar::Text(txt)));
        let q = if z.im.abs() > 0.0 && z.re > 0.0 { z.re / (2.0 * z.im.abs()) } else { f64::INFINITY };
        if q.is_finite() {
            scalars.push((format!("Q {} = Re/(2|Im|) of neff", j + 1), Scalar::Num(q)));
        }
        mode_metrics[j].push_scalars(j + 1, &mut scalars);
    }

    let mut log = vec!["complex mixed mode analysis (Tri+Quad, PML or complex materials)".to_string()];
    if !modes.converged {
        log.push("warning: the eigensolver did not converge every mode to the requested tolerance".into());
    }
    scalars.push(crate::optical::mode_frame(problem).scalar());
    let mut diag = crate::diag::EmsDiag::new("shift_invert_complex");
    diag.record_mode_residuals("shift_invert_complex", modes.residuals.clone(), modes.converged, Some(modes.restarts));
    let (diagnostics, error_budget) = diag.finish();
    Ok(SolveResult {
        physics: "ems".into(),
        study: "modes".into(),
        order: problem.order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.n_tot,
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

pub fn assemble_modes_complex_mixed(
    mesh: &MixedDiscretisation,
    problem: &Problem,
    freq: f64,
) -> Result<(ModeAssemblyC, NedMixed, MixedLagrangeDofMap), String> {
    if mesh.dim != 2 {
        return Err("mode analysis needs a 2D cross-section mesh".into());
    }
    let order = problem.order;
    let ned = NedMixed::new(mesh, order);
    let phi_space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Tensor };
    let phi_map = MixedLagrangeDofMap::new(mesh, phi_space);
    let net = ned.n_dofs;
    let nphi = phi_map.n_dofs;
    let n_tot = net + nphi;
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    let k2 = c(k0 * k0, 0.0);

    let mats = MixedCellMaterials::from_problem(mesh, problem)?;
    let pml_table = EmPmlTable::from_mixed(mesh, &problem.materials, freq, true, false);
    let nc = mesh.n_cells();
    let mut eps_loc = Vec::with_capacity(nc);
    let mut mui_loc = Vec::with_capacity(nc);
    let perm_loc = |t: &C3| -> C3 {
        let axes = [0usize, 2, 1];
        let mut r = [[c(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = t[axes[i]][axes[j]];
            }
        }
        r
    };
    for ci in 0..nc {
        eps_loc.push(perm_loc(&mats.eps_0[ci]));
        mui_loc.push(perm_loc(&mats.mu_inv_0[ci]));
    }
    let is_non_symmetric = mats.is_non_symmetric || problem.sheet_2d.iter().any(Sheet2dBc::is_non_symmetric);

    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut a1_coo = CCoo::new(n_tot, n_tot);
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut has_coupling = false;

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut curl = Vec::new();
    let mut nv_s = Vec::new();
    let mut dxi_s = Vec::new();
    let mut gn = Vec::new();

    for ci in 0..nc {
        let iso = IsoparametricCell::from_mesh(mesh, ci);
        let is_pml = pml_table.is_pml(ci);
        let deg = if iso.curved || is_pml { 6 } else { 4 };
        let rule = quad_rule_mixed(iso.cell_type, deg);
        let nloc_e = nedelec_nloc(iso.cell_type, order);
        let pd = phi_map.cell_dofs(ci);
        let nl = pd.len();
        gn.resize(nl, [0.0; 2]);

        let mut amat = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut tm = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut gm = vec![c(0.0, 0.0); nloc_e * nl];
        let mut gt = vec![c(0.0, 0.0); nl * nloc_e];
        let mut pp = vec![c(0.0, 0.0); nl * nl];
        let mut a1_ee = vec![c(0.0, 0.0); nloc_e * nloc_e];
        let mut mg = vec![c(0.0, 0.0); nloc_e * nl];
        let mut mgt = vec![c(0.0, 0.0); nl * nloc_e];
        let mut cell_cpl = false;
        let mut cell_mu_cpl = false;

        for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
            let jac = eval_cell_basis_2d(mesh, &ned, ci, &iso, xi, &mut scratch, &mut w, &mut curl);
            lagrange_shape_into(iso.cell_type, phi_space, xi, &mut nv_s, &mut dxi_s);
            for a in 0..nl {
                let g3 = jac.map_grad(dxi_s[a]);
                gn[a] = [g3[0], g3[1]];
            }
            let (ep, mui, _) = pml_table.eval_2d(ci, eps_loc[ci], mui_loc[ci], [jac.x[0], jac.x[1]], false);
            add_waveguide_complex_point(
                nloc_e,
                nl,
                &w,
                &curl,
                &nv_s,
                &gn,
                w_ref * jac.det.abs(),
                k2,
                &ep,
                &mui,
                &mut amat,
                &mut tm,
                &mut gm,
                &mut gt,
                &mut pp,
                &mut a1_ee,
                &mut mg,
                &mut mgt,
                &mut cell_cpl,
                &mut cell_mu_cpl,
            );
        }

        let rows_e = ned.cell_dofs(ci);
        let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
        a_coo.push_block(rows_e, rows_e, &amat);
        b_coo.push_block(rows_e, rows_e, &tm);
        b_coo.push_block(rows_e, &cols_p, &gm);
        b_coo.push_block(&cols_p, rows_e, &gt);
        b_coo.push_block(&cols_p, &cols_p, &pp);
        if cell_cpl || cell_mu_cpl {
            has_coupling = true;
            if cell_mu_cpl {
                a1_coo.push_block(rows_e, rows_e, &a1_ee);
            }
            let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
            a1_coo.push_block(rows_e, &cols_p, &neg);
            a1_coo.push_block(&cols_p, rows_e, &mgt);
        }
    }

    add_ibc_2d_waveguide_mixed(
        mesh,
        &ned,
        &phi_map,
        freq,
        &problem.pec,
        &problem.surface_impedance,
        &mut a_coo,
        &mut b_coo,
    )?;
    let sheet_cpl = crate::ibc::add_sheet2d_2d_waveguide_mixed(
        mesh,
        &ned,
        &phi_map,
        freq,
        &problem.pec,
        &problem.sheet_2d,
        problem.imported_potential.as_deref(),
        &mut a_coo,
        &mut b_coo,
        &mut a1_coo,
    )?;
    has_coupling |= sheet_cpl;
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();

    let mut non_pec_labels: Vec<i32> = problem.pmc.clone();
    for sb in &problem.sheet_2d {
        non_pec_labels.extend(sb.labels.iter().copied());
    }
    let (be, bn) = pec_edges_and_nodes_2d_mixed(
        mesh,
        &problem.pec,
        &problem.surface_impedance,
        &non_pec_labels,
        problem.symmetry,
    )?;
    let nv = mesh.n_points();
    let mut fixed = Vec::new();
    if order == 1 {
        fixed.extend(be.iter().copied());
    } else {
        for &e in &be {
            fixed.push(2 * e);
            fixed.push(2 * e + 1);
        }
    }
    fixed.extend(bn.iter().map(|&v| net + v));
    if order == 2 {
        fixed.extend(be.iter().map(|&e| net + nv + e));
    }
    let free = free_dofs(n_tot, &fixed);
    let scale: Vec<f64> = (0..n_tot).map(|i| if i < net { 1.0 } else { 1.0 / k0 }).collect();
    for r in 0..b_full.nrows() {
        for p in b_full.indptr()[r]..b_full.indptr()[r + 1] {
            let s = scale[r] * scale[b_full.indices()[p] as usize];
            b_full.values_mut()[p] *= s;
        }
    }
    let af = a_full.submatrix(&free, &free);
    let bf = b_full.submatrix(&free, &free);
    let a1 = if has_coupling {
        let mut a1 = a1_coo.to_csr();
        for r in 0..a1.nrows() {
            for p in a1.indptr()[r]..a1.indptr()[r + 1] {
                let s = scale[r] * scale[a1.indices()[p] as usize];
                a1.values_mut()[p] *= s;
            }
        }
        Some(a1.submatrix(&free, &free))
    } else {
        None
    };
    Ok((ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0, is_non_symmetric }, ned, phi_map))
}

#[allow(clippy::too_many_arguments)]
fn mode_fields_complex_cells_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    phi_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    modes: &ModesC,
    j: usize,
    freq: f64,
    order: usize,
) -> (Vec<Field>, CellFields) {
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    let jb = c(0.0, 1.0) * (modes.neff[j] * k0);
    let isos: Vec<IsoparametricCell> = (0..mixed.n_cells()).map(|ci| IsoparametricCell::from_mesh(mixed, ci)).collect();
    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut curl = Vec::new();
    let mut nv_s = Vec::new();
    let mut dxi_s = Vec::new();

    let sampled = proj.sample_simplex_cell_nodes(mixed, order, |mc, xi| {
        eval_cell_basis_2d(mixed, ned, mc, &isos[mc], xi, &mut scratch, &mut w, &mut curl);
        lagrange_shape_into(isos[mc].cell_type, phi_map.space, xi, &mut nv_s, &mut dxi_s);
        let cdofs = ned.cell_dofs(mc);
        let pd = phi_map.cell_dofs(mc);
        let (mut x, mut y, mut pval) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        for (s, &gd) in cdofs.iter().enumerate() {
            let cf = modes.et[j][gd];
            x += cf * w[s][0];
            y += cf * w[s][1];
        }
        for (a, &gp) in pd.iter().enumerate() {
            pval += modes.phi[j][gp] * nv_s[a];
        }
        [x, y, jb * pval]
    });

    let ex: Vec<C> = sampled.iter().map(|v| v[0]).collect();
    let ey: Vec<C> = sampled.iter().map(|v| v[1]).collect();
    let ez: Vec<C> = sampled.iter().map(|v| v[2]).collect();

    let big =
        ex.iter()
            .chain(ey.iter())
            .chain(ez.iter())
            .cloned()
            .fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let n = simplex.n_dofs(order);
    let avg = |v: &[f64]| node_average(simplex, v, 1, order);
    let mag = |a: &[C]| -> Vec<f64> { a.iter().map(|z| z.norm()).collect() };
    let re = |a: &[C]| -> Vec<f64> { a.iter().map(|z| (*z * rot).re).collect() };
    let ezn = avg(&mag(&ez));
    let etn_c: Vec<f64> = (0..ex.len()).map(|i| (ex[i].norm_sqr() + ey[i].norm_sqr()).sqrt()).collect();
    let etn = avg(&etn_c);
    let tot: Vec<f64> = (0..n).map(|i| etn[i].hypot(ezn[i])).collect();
    let s = 1.0 / tot.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let f = |name: &str, v: Vec<f64>| Field::scalar(name, "a.u.", v.iter().map(|x| x * s).collect());
    let (lat, nor, prop) = (local_component(0), local_component(1), local_component(2));
    let signed = |c_arr: &[C], sign: f64| -> Vec<f64> { avg(&re(c_arr)).iter().map(|x| x * sign).collect() };

    let mut cells = CellFields::new(simplex, order);
    for (comp, v, sign) in [(lat, &ex, lat.1), (nor, &ey, nor.1), (prop, &ez, prop.1)] {
        let rotated: Vec<C> = v.iter().map(|z| *z * rot * s * sign).collect();
        cells.push_complex(
            format!("E{}", comp.0),
            "a.u.",
            rotated.iter().map(|z| z.re).collect(),
            rotated.iter().map(|z| z.im).collect(),
        );
    }
    let fields = vec![
        f("|E|", tot),
        f("|Exz|", etn),
        f(&format!("|E{}|", prop.0), ezn),
        f(&format!("Re E{}", lat.0), signed(&ex, lat.1)),
        f(&format!("Re E{}", nor.0), signed(&ey, nor.1)),
        f(&format!("Re E{}", prop.0), signed(&ez, prop.1)),
    ];
    (fields, cells)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_waveguide_2d_metrics_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    phi_map: &MixedLagrangeDofMap,
    problem: &Problem,
    freq: f64,
    neff: &[C],
    et: &[Vec<C>],
    phi: &[Vec<C>],
) -> Result<Vec<WaveguideModeMetrics>, String> {
    let n_modes = neff.len();
    if n_modes == 0 {
        return Ok(Vec::new());
    }
    let omega = (2.0 * std::f64::consts::PI * freq).max(1e-300);
    let k0 = omega / C0;
    let mats = MixedCellMaterials::from_problem(mesh, problem)?;
    let pml_cells: Vec<bool> = (0..mesh.n_cells())
        .map(|ci| problem.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(Pml::is_active))
        .collect();
    let any_pml = pml_cells.iter().any(|&p| p);
    let mut domains: Vec<i32> =
        (0..mesh.n_cells()).filter(|&ci| !pml_cells[ci]).map(|ci| mesh.cell_labels[ci]).collect();
    domains.sort_unstable();
    domains.dedup();
    let cell_dom_idx: Vec<usize> = (0..mesh.n_cells())
        .map(|ci| {
            if pml_cells[ci] {
                usize::MAX
            } else {
                domains.binary_search(&mesh.cell_labels[ci]).unwrap_or(usize::MAX)
            }
        })
        .collect();

    let mut p_y = vec![0.0_f64; n_modes];
    let mut u_em = vec![0.0_f64; n_modes];
    let mut comp_e = vec![[0.0_f64; 3]; n_modes];
    let mut dom_e = vec![vec![0.0_f64; domains.len()]; n_modes];
    let mut u_e_phys = vec![0.0_f64; n_modes];
    let mut u_e_pml = vec![0.0_f64; n_modes];

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut curl = Vec::new();
    let mut nv_s = Vec::new();
    let mut dxi_s = Vec::new();
    let mut gn = Vec::new();

    for ci in 0..mesh.n_cells() {
        let iso = IsoparametricCell::from_mesh(mesh, ci);
        let deg = if iso.curved { 6 } else { 4 };
        let rule = quad_rule_mixed(iso.cell_type, deg);
        let cdofs = ned.cell_dofs(ci);
        let pd = phi_map.cell_dofs(ci);
        let nloc_e = cdofs.len();
        let nl = pd.len();
        gn.resize(nl, [0.0; 2]);

        let er = &mats.eps_phys_re[ci];
        let de = &mats.d_omega_eps[ci];
        let mur = mats.mur_phys[ci];
        let inv_om_mu = 1.0 / (omega * MU0 * mur);
        let d_idx = cell_dom_idx[ci];

        for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
            let jac = eval_cell_basis_2d(mesh, ned, ci, &iso, xi, &mut scratch, &mut w, &mut curl);
            lagrange_shape_into(iso.cell_type, phi_map.space, xi, &mut nv_s, &mut dxi_s);
            for a in 0..nl {
                let g3 = jac.map_grad(dxi_s[a]);
                gn[a] = [g3[0], g3[1]];
            }
            let wq_area = w_ref * jac.det.abs();
            for j in 0..n_modes {
                let beta = neff[j] * k0;
                let j_beta = c(0.0, 1.0) * beta;
                let (mut ex, mut ez, mut curl_2d) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                for s in 0..nloc_e {
                    let cf = et[j][cdofs[s]];
                    ex += cf * w[s][0];
                    ez += cf * w[s][1];
                    curl_2d += cf * curl[s];
                }
                let (mut phi_q, mut dphi_dx, mut dphi_dz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
                for a in 0..nl {
                    let pf = phi[j][pd[a]];
                    phi_q += pf * nv_s[a];
                    dphi_dx += pf * gn[a][0];
                    dphi_dz += pf * gn[a][1];
                }
                let ey = j_beta * phi_q;
                let e_w = [ex, ey, ez];
                let mut e_eps_e = 0.0_f64;
                for a in 0..3 {
                    for b in 0..3 {
                        e_eps_e += er[a][b] * (e_w[a].conj() * e_w[b]).re;
                    }
                }
                if pml_cells[ci] {
                    u_e_pml[j] += e_eps_e * wq_area;
                } else {
                    u_e_phys[j] += e_eps_e * wq_area;
                    for k in 0..3 {
                        comp_e[j][k] += er[k][k].abs() * e_w[k].norm_sqr() * wq_area;
                    }
                    if d_idx != usize::MAX {
                        dom_e[j][d_idx] += e_eps_e * wq_area;
                    }
                    let hx = beta * (ez + dphi_dz) * inv_om_mu;
                    let hy = -c(0.0, 1.0) * curl_2d * inv_om_mu;
                    let hz = -beta * (ex + dphi_dx) * inv_om_mu;
                    p_y[j] += 0.5 * (ez.conj() * hx - ex.conj() * hz).re * wq_area;
                    let mut e_de_e = 0.0_f64;
                    for a in 0..3 {
                        for b in 0..3 {
                            e_de_e += de[a][b] * (e_w[a].conj() * e_w[b]).re;
                        }
                    }
                    let h_mu_h = mur * (hx.norm_sqr() + hy.norm_sqr() + hz.norm_sqr());
                    u_em[j] += 0.25 * (EPS0 * e_de_e + MU0 * h_mu_h) * wq_area;
                }
            }
        }
    }

    let mut out = Vec::with_capacity(n_modes);
    for j in 0..n_modes {
        let vg_raw = if u_em[j] > 0.0 { p_y[j] / u_em[j] } else { 0.0 };
        let vg = if vg_raw.abs() <= 1e-9 * C0 { 0.0 } else { vg_raw };
        let comp_tot = comp_e[j][0] + comp_e[j][1] + comp_e[j][2];
        let comp_shares = if comp_tot > 0.0 {
            [comp_e[j][0] / comp_tot, comp_e[j][1] / comp_tot, comp_e[j][2] / comp_tot]
        } else {
            [0.0; 3]
        };
        let dom_tot: f64 = dom_e[j].iter().sum();
        let domain_shares = domains
            .iter()
            .enumerate()
            .map(|(d, &label)| (label, if dom_tot > 0.0 { dom_e[j][d] / dom_tot } else { 0.0 }))
            .collect();
        let u_all = u_e_phys[j] + u_e_pml[j];
        let layer_fraction = if any_pml && u_all > 0.0 { (u_e_pml[j] / u_all).clamp(0.0, 1.0) } else { 0.0 };
        out.push(WaveguideModeMetrics { vg, comp_shares, domain_shares, layer_fraction });
    }
    Ok(out)
}

/// Computes the biorthogonal cross-section Poynting inner product
/// `1/2 int_Omega (E_{t,left}^* x H_{t,right}) . y_hat dOmega` between left mode `left_idx`
/// (`modes.et_left[left_idx]`) and right mode `right_idx` (`modes.et[right_idx]`, `modes.phi[right_idx]`)
/// with full 3x3 inverse permeability tensors `mu[ci].inv_matrix()`.
pub fn biorthogonal_poynting_overlap(
    mesh: &Discretisation,
    freq: f64,
    _eps: &[EpsC],
    mu: &[MuC],
    order: usize,
    modes: &ModesC,
    left_idx: usize,
    right_idx: usize,
) -> C {
    if left_idx >= modes.et_left.len() || right_idx >= modes.et.len() {
        return c(0.0, 0.0);
    }
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let omega = (2.0 * std::f64::consts::PI * freq).max(1e-300);
    let k0 = omega / C0;
    let beta_r = modes.neff[right_idx] * k0;
    let j_beta = c(0.0, 1.0) * beta_r;
    let h_pref = c(0.0, 1.0 / (omega * MU0));
    let ps = femcore::section::PROP_SIGN;
    let et_l = &modes.et_left[left_idx];
    let et_r = &modes.et[right_idx];
    let phi_r = &modes.phi[right_idx];

    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];

    let mut total = c(0.0, 0.0);
    for ci in 0..mesh.n_cells() {
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        mesh.cell_dofs(ci, order, &mut pd);
        let mui = mu[ci].inv_matrix();
        let mui_local: C3 = [
            [mui[0][0], mui[0][1], mui[0][2] * ps],
            [mui[1][0], mui[1][1], mui[1][2] * ps],
            [mui[2][0] * ps, mui[2][1] * ps, mui[2][2]],
        ];
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
            let (_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let wq = quad.w[q] * area;
            let (mut el0, mut el1) = (c(0.0, 0.0), c(0.0, 0.0));
            let (mut er0, mut er1, mut curl_2d) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let dof = cd.gdofs[ci][s];
                let cl = et_l[dof];
                let cr = et_r[dof];
                el0 += cl * w[s][0];
                el1 += cl * w[s][1];
                er0 += cr * w[s][0];
                er1 += cr * w[s][1];
                curl_2d += cr * curl[s];
            }
            let (mut dp0, mut dp1) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nl {
                let pr = phi_r[pd[a]];
                dp0 += pr * gn[a][0];
                dp1 += pr * gn[a][1];
            }
            let curl_e = [j_beta * (dp1 + er1), -j_beta * (dp0 + er0), curl_2d];
            let mce = crate::mixed_ops::c3_mul_vec(&mui_local, curl_e);
            let h0 = mce[0] * h_pref;
            let h1 = mce[1] * h_pref;
            total += c(0.5 * wq, 0.0) * (el1.conj() * h0 - el0.conj() * h1);
        }
    }
    total
}
