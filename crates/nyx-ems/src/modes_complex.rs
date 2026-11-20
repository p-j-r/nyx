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

use crate::modes::C0;
use crate::problem::{Problem, Study};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::*;
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::eigs_complex;
use femcore::optical::{cell_materials, EpsC, MuC};
use femcore::section::local_component;
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
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
    let a_full = a_coo.to_csr();
    let mut b_full = b_coo.to_csr();
    // PEC walls
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
        neff.push(nf);
        let mut x = vec![c(0.0, 0.0); n_tot];
        for (a, &dof) in free.iter().enumerate() {
            x[dof] = pairs.vectors[i][a];
        }
        et_all.push(x[..net].to_vec());
        phi_all.push(x[net..].iter().map(|v| *v / k0).collect::<Vec<C>>());
    }
    let _ = n_wanted;
    // sort by closeness to the target
    let mut idx: Vec<usize> = (0..neff.len()).collect();
    idx.sort_by(|&a, &b| (neff[a] - n_target).norm().partial_cmp(&(neff[b] - n_target).norm()).unwrap());
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
    let t0 = Instant::now();
    let Study::Modes { freq, n_target, evanescent, n_modes } = problem.study.clone() else {
        return Err("mode analysis needs a Modes study".into());
    };
    if freq <= 0.0 || n_modes < 1 {
        return Err("frequency must be positive and at least one mode is needed".into());
    }
    let (eps, mu) = cell_materials(mesh, &problem.materials)?;
    let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
    let modes = solve_modes_complex(mesh, freq, &eps, &mu, target, n_modes, problem.order)?;
    let mut solutions = vec![];
    let mut scalars = vec![];
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
    };
    Ok(result)
}
