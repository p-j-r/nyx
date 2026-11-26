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
use crate::mixed_ops::{complex_eps_c3, iso_mu_inv_c3, EmPmlTable, MixedCellMaterials, C3};
use crate::modes::C0;
use crate::problem::{pec_edges_and_nodes_2d_with_symmetry, Problem, Study, SurfaceImpedanceBc, Symmetry};
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
use femcore::section::{local_component, local_tensor};
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
/// that is linear in `beta` when the permittivity couples the transverse and the longitudinal field (tilted crystal).
pub struct ModeAssemblyC {
    pub a: CCsr,
    pub b: CCsr,
    /// `A1` of `(A + j beta A1) x = -beta^2 B x`, present only for a tilted crystal
    pub a1: Option<CCsr>,
    pub free: Vec<usize>,
    pub n_tot: usize,
    pub net: usize,
    pub k0: f64,
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
    }
    add_ibc_2d_waveguide(mesh, order, freq, pec_labels, ibc_bcs, &mut a_coo, &mut b_coo)?;
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();
    // PEC walls
    let (be, bn) = if pec_labels.is_empty() && ibc_bcs.is_empty() && pmc_labels.is_empty() && symmetry == Symmetry::None
    {
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
        pec_edges_and_nodes_2d(mesh, pec_labels, ibc_bcs, pmc_labels)
    } else {
        pec_edges_and_nodes_2d_with_symmetry(mesh, pec_labels, ibc_bcs, pmc_labels, symmetry)?
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
    Ok(ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0 })
}

/// Stage 2: the `k` eigenpairs nearest the complex effective index `n_target`, as effective indices, edge coefficients
/// and `phi = Ez / (j beta)`, sorted by closeness to it. Pairs the Krylov process could not converge are dropped.
pub fn solve_mode_pairs_complex(asm: &ModeAssemblyC, n_target: C, k: usize) -> Result<ModesC, String> {
    let (af, bf, free, n_tot, net, k0) = (&asm.a, &asm.b, &asm.free, asm.n_tot, asm.net, asm.k0);
    let sigma = -(n_target * n_target) * (k0 * k0);
    let mut pairs = eigs_complex(af, bf, sigma, k, 1e-10, 400)?;
    if let Some(a1f) = &asm.a1 {
        // refine every mode by the fixed-point iteration on b in the linear term
        for i in 0..pairs.values.len() {
            let pick = |z: C| -> C {
                // forward-propagating (positive real part) or, for an evanescent mode, decaying (negative imaginary part) root
                if z.re.abs() < 1e-9 * z.norm() {
                    if z.im > 0.0 {
                        -z
                    } else {
                        z
                    }
                } else if z.re < 0.0 {
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
                let nb = pick((-one.values[0]).sqrt());
                let change = (nb - beta).norm() / beta.norm();
                femcore::progress::note(format!("tilted-crystal coupling, optical mode {}: refinement step {}, relative change of the propagation constant {change:.1e} (stop below {COUPLED_BETA_TOL:.0e})", i + 1, it + 1));
                hist.push(beta);
                beta = nb;
                last = Some(one);
                if change < COUPLED_BETA_TOL {
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
    let n_wanted = pairs.values.len();
    let frob = |m: &CCsr| m.values().iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
    let norm_a0 = frob(af);
    let norm_a1 = asm.a1.as_ref().map(frob).unwrap_or(0.0);
    let norm_b = frob(bf);
    for (i, &lam) in pairs.values.iter().enumerate() {
        // pairs the Krylov process could not converge (clustered radiation modes in the absorbing layers) are dropped
        if pairs.residuals.get(i).is_some_and(|r| *r > 1e-6) {
            continue;
        }
        let beta2 = -lam;
        let mut nf = beta2.sqrt() / k0;
        if nf.re < 0.0 {
            nf = -nf; // choose the forward-propagating (or decaying) root
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
        if res_true > 1e-6 {
            continue;
        }
        neff.push(nf);
        let mut x = vec![c(0.0, 0.0); n_tot];
        for (a, &dof) in free.iter().enumerate() {
            x[dof] = u[a];
        }
        et_all.push(x[..net].to_vec());
        phi_all.push(x[net..].iter().map(|v| *v / k0).collect::<Vec<C>>());
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
    let asm = if needs_quad_pml {
        assemble_modes_complex_quad_pml(mesh, problem, freq, &pml_table)?
    } else if problem.pec.is_empty()
        && problem.surface_impedance.is_empty()
        && problem.pmc.is_empty()
        && problem.symmetry == Symmetry::None
    {
        let (ep_loc, mu_loc) = femcore::optical::cell_materials(mesh, &problem.materials)?;
        assemble_modes_complex(mesh, freq, &ep_loc, &mu_loc, problem.order)?
    } else {
        assemble_modes_complex_with_bcs_and_symmetry(
            mesh,
            freq,
            &eps,
            &mu,
            problem.order,
            &problem.pec,
            &problem.surface_impedance,
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
            mode_metrics.truncate(n_modes);
        }
    }
    for j in 0..modes.neff.len() {
        let z = modes.neff[j];
        let txt = format!("{:.6} {:+.3e}i", z.re, z.im);
        let (fields, cell_fields) = mode_fields_complex_cells(mesh, &modes, j, freq, problem.order);
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
        ..Default::default()
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
    pp: &mut [C],
    mg: &mut [C],
    has_coupling: &mut bool,
) {
    // Rotated transverse inverse permeability R^T [[nu00, nu01], [nu01, nu11]] R = [[nu11, -nu01], [-nu01, nu00]]
    let m00 = mui[1][1];
    let m01 = -mui[0][1];
    let m11 = mui[0][0];
    let imz = mui[2][2];
    let cpl = ep[0][2].norm_sqr() + ep[1][2].norm_sqr() > 0.0;
    if cpl {
        *has_coupling = true;
    }

    for s in 0..nloc_e {
        let ws = w[s];
        let mw_s = [m00 * ws[0] + m01 * ws[1], m01 * ws[0] + m11 * ws[1]];
        let ew_s = [ep[0][0] * ws[0] + ep[0][1] * ws[1], ep[1][0] * ws[0] + ep[1][1] * ws[1]];
        for t in 0..nloc_e {
            let wt = w[t];
            let idx = s * nloc_e + t;
            amat[idx] += (imz * (curl[s] * curl[t]) - k2 * (ew_s[0] * wt[0] + ew_s[1] * wt[1])) * wq;
            tm[idx] += (mw_s[0] * wt[0] + mw_s[1] * wt[1]) * wq;
        }
        for a in 0..nl {
            let ga = gn[a];
            gm[s * nl + a] += (mw_s[0] * ga[0] + mw_s[1] * ga[1]) * wq;
            if cpl {
                mg[s * nl + a] += k2 * (ep[0][2] * ws[0] + ep[1][2] * ws[1]) * (nv[a] * wq);
            }
        }
    }
    for a in 0..nl {
        let ga = gn[a];
        let mga = [m00 * ga[0] + m01 * ga[1], m01 * ga[0] + m11 * ga[1]];
        for b in 0..nl {
            let gb = gn[b];
            pp[a * nl + b] += (mga[0] * gb[0] + mga[1] * gb[1] - k2 * ep[2][2] * (nv[a] * nv[b])) * wq;
        }
    }
}

/// Quadrature-point curvilinear PML assembly on a 2D `SimplexMesh` (`Discretisation`).
fn assemble_modes_complex_quad_pml(
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

    let tens = problem.materials.cell_complex_permittivity(mesh, Some(1.0))?;
    let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let mut eps_loc = Vec::with_capacity(nc);
    let mut mui_loc = Vec::with_capacity(nc);
    for ci in 0..nc {
        let re_l = local_tensor(&tens[ci].0);
        let im_l = local_tensor(&tens[ci].1);
        eps_loc.push(complex_eps_c3(&re_l, &im_l));
        mui_loc.push(iso_mu_inv_c3(mur[ci]));
    }

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
        let mut pp = vec![c(0.0, 0.0); nl * nl];
        let mut mg = vec![c(0.0, 0.0); nloc_e * nl];
        let mut cell_cpl = false;

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
                &mut pp,
                &mut mg,
                &mut cell_cpl,
            );
        }

        mesh.cell_dofs(ci, order, &mut pd);
        let rows_e = &cd.gdofs[ci];
        let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
        a_coo.push_block(rows_e, rows_e, &amat);
        b_coo.push_block(rows_e, rows_e, &tm);
        b_coo.push_block(rows_e, &cols_p, &gm);
        let gt: Vec<C> = (0..nl * nloc_e).map(|i| gm[(i % nloc_e) * nl + i / nloc_e]).collect();
        b_coo.push_block(&cols_p, rows_e, &gt);
        b_coo.push_block(&cols_p, &cols_p, &pp);
        if cell_cpl {
            has_coupling = true;
            let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
            a1_coo.push_block(rows_e, &cols_p, &neg);
            let mgt: Vec<C> = (0..nl * nloc_e).map(|i| mg[(i % nloc_e) * nl + i / nloc_e]).collect();
            a1_coo.push_block(&cols_p, rows_e, &mgt);
        }
    }

    add_ibc_2d_waveguide(mesh, order, freq, &problem.pec, &problem.surface_impedance, &mut a_coo, &mut b_coo)?;
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();
    let (be, bn) = if problem.pec.is_empty()
        && problem.surface_impedance.is_empty()
        && problem.pmc.is_empty()
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
        pec_edges_and_nodes_2d(mesh, &problem.pec, &problem.surface_impedance, &problem.pmc)
    } else {
        pec_edges_and_nodes_2d_with_symmetry(
            mesh,
            &problem.pec,
            &problem.surface_impedance,
            &problem.pmc,
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
    Ok(ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0 })
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
            mode_metrics.truncate(n_modes);
        }
    }

    let mut solutions = Vec::with_capacity(modes.neff.len());
    let mut scalars = Vec::new();
    for j in 0..modes.neff.len() {
        let z = modes.neff[j];
        let txt = format!("{:.6} {:+.3e}i", z.re, z.im);
        let (fields, cell_fields) =
            mode_fields_complex_cells_mixed(mixed, &ned, &phi_map, &simplex, &proj, &modes, j, freq, problem.order);
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
        ..Default::default()
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
    for ci in 0..nc {
        let re_l = local_tensor(&mats.eps_phys_re[ci]);
        let im_l = local_tensor(&mats.eps_phys_im[ci]);
        eps_loc.push(complex_eps_c3(&re_l, &im_l));
        mui_loc.push(iso_mu_inv_c3(mats.mur_phys[ci]));
    }

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
        let mut pp = vec![c(0.0, 0.0); nl * nl];
        let mut mg = vec![c(0.0, 0.0); nloc_e * nl];
        let mut cell_cpl = false;

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
                &mut pp,
                &mut mg,
                &mut cell_cpl,
            );
        }

        let rows_e = ned.cell_dofs(ci);
        let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
        a_coo.push_block(rows_e, rows_e, &amat);
        b_coo.push_block(rows_e, rows_e, &tm);
        b_coo.push_block(rows_e, &cols_p, &gm);
        let gt: Vec<C> = (0..nl * nloc_e).map(|i| gm[(i % nloc_e) * nl + i / nloc_e]).collect();
        b_coo.push_block(&cols_p, rows_e, &gt);
        b_coo.push_block(&cols_p, &cols_p, &pp);
        if cell_cpl {
            has_coupling = true;
            let neg: Vec<C> = mg.iter().map(|v| -*v).collect();
            a1_coo.push_block(rows_e, &cols_p, &neg);
            let mgt: Vec<C> = (0..nl * nloc_e).map(|i| mg[(i % nloc_e) * nl + i / nloc_e]).collect();
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
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();

    let (be, bn) =
        pec_edges_and_nodes_2d_mixed(mesh, &problem.pec, &problem.surface_impedance, &problem.pmc, problem.symmetry)?;
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
    Ok((ModeAssemblyC { a: af, b: bf, a1, free, n_tot, net, k0 }, ned, phi_map))
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
