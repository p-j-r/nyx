// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Full-vector mode analysis of a 2D waveguide cross-section (inhomogeneous media allowed),
//! perfectly conducting walls.
//!
//! The solver works in the local axes of the guide, (lateral, normal, propagation) = (x, z, `PROP_SIGN` y) of the wafer (`section`); its
//! internal names for them are x, y, z below. The problem comes in wafer axes (tensors are turned into the local axes by
//! `section::local_tensor`) and the result fields leave under wafer names (`Ex`, `Ez` in the plane, `Ey` along the guide).
//!
//! Fields E = (Et + z Ez) exp(-j beta z). Inserting this in the curl-curl equation documented
//! for the electromagnetic-waves problem, curl(mu_r^-1 curl E) - k0^2 eps_r E = 0, gives
//! the textbook mixed formulation with Nedelec edge elements for Et and Lagrange elements of
//! one degree higher for phi = Ez / (j beta):
//!
//!   [ S - k0^2 Te      0 ] [et ]          [ Tm    G            ] [et ]
//!   [      0           0 ] [phi] = -b^2   [ G^T   Kp - k0^2 Tp ] [phi]
//!
//! (S curl-curl with 1/mu_r, Te and Tm edge mass matrices with eps_r and 1/mu_r, G the coupling
//! integral of (1/mu_r) w . grad(lambda), Kp and Tp the Lagrange stiffness and mass.) The
//! spurious gradient modes sit at b^2 = 0, so a shift-and-invert search near the wanted
//! effective index finds physical modes only. `order` 1 is Whitney edge elements with P1 phi;
//! order 2 (default) is hierarchical degree-2 Nedelec elements (8 per triangle: 3 Whitney, 3
//! gradients of edge bubbles, 2 interior) with P2 phi, the accuracy of quadratic
//! elements on the same mesh.

use crate::problem::{Problem, Study};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::*;
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{Coo, RCsr as Csr};
use femcore::linalg::driver::eigs_general;
use femcore::linalg::driver::EigOpts;
use femcore::props::TensorKind;
use femcore::section::{local_component, local_tensor};
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use std::time::Instant;

pub use femcore::optical::C0;
/// Relative permittivity of one cell in the local axes (lateral x, normal y, propagation z): the transverse 2x2 block (xx, xy, yy) and zz. The
/// cross-section must be oriented so that the propagation axis is a principal axis of the
/// permittivity (no xz or yz terms).
#[derive(Clone, Copy, Debug)]
pub struct EpsCell {
    pub xx: f64,
    pub xy: f64,
    pub yy: f64,
    pub zz: f64,
}

pub struct Modes {
    /// effective indices (re, im), sorted by closeness to the target
    pub neff: Vec<(f64, f64)>,
    /// edge-element coefficients per mode (length = number of edge functions)
    pub et: Vec<Vec<f64>>,
    /// phi = Ez / (j beta) per mode at the P1/P2 dofs
    pub phi: Vec<Vec<f64>>,
    pub converged: bool,
    /// a converged complex pair of the real pencil (a leaky mode or one with a complex propagation constant) lies among the wanted
    /// eigenvalues: this real solve cannot represent it, and `solve_emw` repeats the study in the complex solver of `ems_pml`.
    pub complex_wanted: bool,
}

/// Modes nearest the effective index `n_target` (real, or imaginary if `evanescent`) at
/// frequency `freq` (Hz). `eps_r` and `mu_r` hold one value per cell.
#[allow(clippy::too_many_arguments)]
pub fn solve_modes(
    mesh: &Discretisation,
    freq: f64,
    eps_r: &[f64],
    mu_r: &[f64],
    n_target: f64,
    evanescent: bool,
    k: usize,
    order: usize,
) -> Result<Modes, String> {
    let eps: Vec<EpsCell> = eps_r.iter().map(|&e| EpsCell { xx: e, xy: 0.0, yy: e, zz: e }).collect();
    solve_modes_tensor(mesh, freq, &eps, mu_r, n_target, evanescent, k, order)
}

/// As `solve_modes` with a permittivity tensor per cell.
#[allow(clippy::too_many_arguments)]
pub fn solve_modes_tensor(
    mesh: &Discretisation,
    freq: f64,
    eps_r: &[EpsCell],
    mu_r: &[f64],
    n_target: f64,
    evanescent: bool,
    k: usize,
    order: usize,
) -> Result<Modes, String> {
    let asm = assemble_modes(mesh, freq, eps_r, mu_r, order)?;
    solve_mode_pairs(&asm, n_target, evanescent, k)
}

/// Stage 1 of the real mode solver: the generalized eigenproblem `A x = lambda B x` (`lambda = -beta^2`) of the
/// edge-element + P2-phi formulation, with the PEC walls removed and the phi block scaled by `1/k0`.
pub struct ModeAssembly {
    /// free-dof system matrices
    pub a: Csr,
    pub b: Csr,
    /// the free dofs of the full vector (edge functions first, then phi)
    pub free: Vec<usize>,
    pub n_tot: usize,
    /// edge-function dofs (the first block of the full vector)
    pub net: usize,
    pub k0: f64,
}

pub fn assemble_modes(
    mesh: &Discretisation,
    freq: f64,
    eps_r: &[EpsCell],
    mu_r: &[f64],
    order: usize,
) -> Result<ModeAssembly, String> {
    if mesh.dim != 2 {
        return Err("mode analysis needs a 2D cross-section mesh".into());
    }
    if order != 1 && order != 2 {
        return Err("element order must be 1 or 2".into());
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
    let mut a_coo = Coo::new(n_tot, n_tot);
    let mut b_coo = Coo::new(n_tot, n_tot);
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    for c in 0..nc {
        let g = tri_g(&geo, c);
        let area = geo.vol[c];
        let (ec, imu) = (eps_r[c], 1.0 / mu_r[c]);
        let eps = ec.zz;
        let mut wm = vec![0.0; nloc_e * nloc_e];
        let mut wme = vec![0.0; nloc_e * nloc_e];
        let mut se = vec![0.0; nloc_e * nloc_e];
        let mut ge = vec![0.0; nloc_e * nl];
        let mut kp = vec![0.0; nl * nl];
        let mut tp = vec![0.0; nl * nl];
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[c], order, &mut w, &mut curl);
            let (nv_, dn) = &tab[q];
            // gradients of the phi shape functions
            let mut gn = vec![[0.0; 2]; nl];
            for a in 0..nl {
                for kk in 0..2 {
                    gn[a][kk] = (0..3).map(|b| dn[a * 3 + b] * g[b][kk]).sum();
                }
            }
            let wq = quad.w[q];
            for s in 0..nloc_e {
                for t in 0..nloc_e {
                    wm[s * nloc_e + t] += wq * (w[s][0] * w[t][0] + w[s][1] * w[t][1]);
                    wme[s * nloc_e + t] += wq
                        * (w[s][0] * (ec.xx * w[t][0] + ec.xy * w[t][1])
                            + w[s][1] * (ec.xy * w[t][0] + ec.yy * w[t][1]));
                    se[s * nloc_e + t] += wq * curl[s] * curl[t];
                }
                for a in 0..nl {
                    ge[s * nl + a] += wq * (w[s][0] * gn[a][0] + w[s][1] * gn[a][1]);
                }
            }
            for a in 0..nl {
                for b in 0..nl {
                    kp[a * nl + b] += wq * (gn[a][0] * gn[b][0] + gn[a][1] * gn[b][1]);
                    tp[a * nl + b] += wq * nv_[a] * nv_[b];
                }
            }
        }
        mesh.cell_dofs(c, order, &mut pd);
        let rows_e = &cd.gdofs[c];
        let cols_p: Vec<usize> = pd.iter().map(|&p| net + p).collect();
        let sc = |m: &[f64], f: f64| -> Vec<f64> { m.iter().map(|x| x * f * area).collect() };
        // A = [S - k0^2 Te, 0; 0, 0]
        let s_minus: Vec<f64> = (0..nloc_e * nloc_e).map(|i| (se[i] * imu - k0 * k0 * wme[i]) * area).collect();
        a_coo.push_block(rows_e, rows_e, &s_minus);
        // B = [Tm, G; G^T, Kphi - k0^2 Tphi]
        b_coo.push_block(rows_e, rows_e, &sc(&wm, imu));
        let gm = sc(&ge, imu);
        b_coo.push_block(rows_e, &cols_p, &gm);
        let gt: Vec<f64> = (0..nl * nloc_e).map(|i| gm[(i % nloc_e) * nl + i / nloc_e]).collect();
        b_coo.push_block(&cols_p, rows_e, &gt);
        let pp: Vec<f64> = (0..nl * nl).map(|i| (kp[i] * imu - k0 * k0 * tp[i] * eps) * area).collect();
        b_coo.push_block(&cols_p, &cols_p, &pp);
    }
    let a_full = a_coo.to_csr();
    let b_full = b_coo.to_csr();

    // PEC: tangential Et = 0 on boundary edges (both edge functions of order 2), phi = 0 on
    // boundary nodes (and boundary edge midpoints for P2)
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
    // scale the phi block by 1/k0 so both blocks have comparable magnitude
    let scale: Vec<f64> = (0..n_tot).map(|i| if i < net { 1.0 } else { 1.0 / k0 }).collect();
    let mut bs = b_full.clone();
    for r in 0..bs.nrows() {
        for p in bs.row_range(r) {
            let s = scale[r] * scale[bs.indices()[p] as usize];
            bs.values_mut()[p] *= s;
        }
    }
    Ok(ModeAssembly { a: a_full.submatrix(&free, &free), b: bs.submatrix(&free, &free), free, n_tot, net, k0 })
}

/// Stage 2: the `k` eigenpairs nearest the effective index `n_target` (real, or imaginary if `evanescent`), sorted by
/// closeness to it, as effective indices, edge coefficients and `phi = Ez / (j beta)`.
pub fn solve_mode_pairs(asm: &ModeAssembly, n_target: f64, evanescent: bool, k: usize) -> Result<Modes, String> {
    let (k0, net, n_tot, free) = (asm.k0, asm.net, asm.n_tot, &asm.free);
    let sigma = if evanescent { (n_target * k0).powi(2) } else { -(n_target * k0).powi(2) };
    let pairs = eigs_general(&asm.a, &asm.b, sigma, k, &EigOpts::default())?;
    let mut out = vec![];
    for (i, &lam) in pairs.values.iter().enumerate() {
        let beta2 = -lam;
        let neff = if beta2 >= 0.0 { (beta2.sqrt() / k0, 0.0) } else { (0.0, (-beta2).sqrt() / k0) };
        let mut x = vec![0.0; n_tot];
        for (a, &dof) in free.iter().enumerate() {
            x[dof] = pairs.vectors[i][a];
        }
        let et = x[..net].to_vec();
        let phi: Vec<f64> = x[net..].iter().map(|v| v / k0).collect();
        out.push((neff, et, phi));
    }
    // sort by closeness to the (real or imaginary) target
    let target = if evanescent { (0.0, n_target) } else { (n_target, 0.0) };
    out.sort_by(|a, b| {
        let da = ((a.0 .0 - target.0).powi(2) + (a.0 .1 - target.1).powi(2)).sqrt();
        let db = ((b.0 .0 - target.0).powi(2) + (b.0 .1 - target.1).powi(2)).sqrt();
        da.partial_cmp(&db).unwrap()
    });
    Ok(Modes {
        neff: out.iter().map(|o| o.0).collect(),
        et: out.iter().map(|o| o.1.clone()).collect(),
        phi: out.iter().map(|o| o.2.clone()).collect(),
        converged: pairs.converged,
        complex_wanted: pairs.complex_found,
    })
}

/// Nodal fields of one mode: transverse and longitudinal field at every dof (node average of
/// the element values), scaled to max |E| = 1, under wafer names: `|E|`, `|Exz|` (in the plane), `|Ey|` and `Ey` (along the guide), `Ex`, `Ez`.
pub fn mode_fields(mesh: &Discretisation, modes: &Modes, j: usize, freq: f64, order: usize) -> Vec<Field> {
    mode_fields_cells(mesh, modes, j, freq, order).0
}

/// Stage 3: `mode_fields` and, next to them, the components `Ex`, `Ez` (in the plane) and `Ey` (along the guide) at the
/// local dofs of every cell, with the same scale and sign as the plot fields (whose node average they reproduce).
pub fn mode_fields_cells(
    mesh: &Discretisation,
    modes: &Modes,
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
    let mut etx = vec![0.0; nc * nloc];
    let mut ety = vec![0.0; nc * nloc];
    for c in 0..nc {
        let g = tri_g(&geo, c);
        for (a, l) in lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[c], order, &mut w, &mut curl);
            let (mut x, mut y) = (0.0, 0.0);
            for s in 0..nloc_e {
                let coef = modes.et[j][cd.gdofs[c][s]];
                x += coef * w[s][0];
                y += coef * w[s][1];
            }
            etx[c * nloc + a] = x;
            ety[c * nloc + a] = y;
        }
    }
    // Ez = j beta phi: magnitude |beta| |phi|; beta is real for propagating modes, imaginary
    // (Ez real) for evanescent ones
    let beta = (modes.neff[j].0 * k0, modes.neff[j].1 * k0);
    let bmag = beta.0.hypot(beta.1);
    let ezs: Vec<f64> = modes.phi[j].iter().map(|p| bmag * p).collect(); // signed profile of Ez up to a phase
    let n = mesh.n_dofs(order);
    let ez_nodal_signed = ezs;
    let mut pd = Vec::new();
    let mut etn_c = vec![0.0; nc * nloc];
    let mut ezn_c = vec![0.0; nc * nloc];
    for c in 0..nc {
        mesh.cell_dofs(c, order, &mut pd);
        for a in 0..nloc {
            etn_c[c * nloc + a] = etx[c * nloc + a].hypot(ety[c * nloc + a]);
            ezn_c[c * nloc + a] = ez_nodal_signed[pd[a]].abs();
        }
    }
    let etn = node_average(mesh, &etn_c, 1, order);
    let ezn = node_average(mesh, &ezn_c, 1, order);
    let etxn = node_average(mesh, &etx, 1, order);
    let etyn = node_average(mesh, &ety, 1, order);
    let tot: Vec<f64> = (0..n).map(|i| etn[i].hypot(ezn[i])).collect();
    let s = 1.0 / tot.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let sc = |v: &[f64]| v.iter().map(|x| x * s).collect::<Vec<_>>();
    // sign convention: the largest transverse component positive
    let big = etxn.iter().chain(etyn.iter()).cloned().fold(0.0f64, |m, x| if x.abs() > m.abs() { x } else { m });
    let sg = if big < 0.0 { -1.0 } else { 1.0 };
    let f = |name: &str, v: Vec<f64>| Field::scalar(name, "a.u.", v);
    let (lat, nor, prop) = (local_component(0), local_component(1), local_component(2));
    let mut cells = CellFields::new(mesh, order);
    let ez_cell: Vec<f64> = (0..nc)
        .flat_map(|c| {
            let mut pd = Vec::new();
            mesh.cell_dofs(c, order, &mut pd);
            pd.into_iter().take(nloc).map(|d| ez_nodal_signed[d] * s * sg * prop.1).collect::<Vec<_>>()
        })
        .collect();
    cells.push_real(format!("E{}", lat.0), "a.u.", etx.iter().map(|x| x * s * sg * lat.1).collect());
    cells.push_real(format!("E{}", nor.0), "a.u.", ety.iter().map(|x| x * s * sg * nor.1).collect());
    cells.push_real(format!("E{}", prop.0), "a.u.", ez_cell);
    let fields = vec![
        f("|E|", sc(&tot)),
        f("|Exz|", sc(&etn)),
        f(&format!("|E{}|", prop.0), sc(&ezn)),
        f(&format!("E{}", prop.0), ez_nodal_signed.iter().map(|x| x * s * sg * prop.1).collect()),
        f(&format!("E{}", lat.0), etxn.iter().map(|x| x * s * sg * lat.1).collect()),
        f(&format!("E{}", nor.0), etyn.iter().map(|x| x * s * sg * nor.1).collect()),
    ];
    (fields, cells)
}

/// The guided modes of a 2D cross-section for a `Modes` study: stages 1 to 3 of the real solver, or the complex solver (`modes_complex`)
/// when the problem has absorbing layers, a tilted crystal, or a complex pair among the wanted eigenvalues. The result carries the
/// polarization of the modes (scalars and log, as before the move; the analysis itself is in `Solved::analysis`).
pub fn solve_waveguide(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    if problem.materials.has_pml()
        || problem.materials.has_complex_permittivity()
        || !problem.surface_impedance.is_empty()
    {
        return crate::modes_complex::solve_waveguide_complex(mesh, problem);
    }
    let t0 = Instant::now();
    let Study::Modes { freq, n_target, evanescent, n_modes } = problem.study.clone() else {
        return Err("mode analysis needs a Modes study".into());
    };
    if freq <= 0.0 {
        return Err("frequency must be positive".into());
    }
    if n_modes < 1 {
        return Err("at least one mode is needed".into());
    }
    let tens: Vec<_> =
        problem.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?.iter().map(local_tensor).collect();
    for t in &tens {
        let scale = t[0][0].abs().max(t[1][1].abs()).max(t[2][2].abs()).max(1e-300);
        if t[0][2].abs() > 1e-9 * scale || t[1][2].abs() > 1e-9 * scale {
            // a crystal whose axes are tilted against the propagation axis: the complex solver handles the
            // coupling of the longitudinal and the transverse field exactly
            return crate::modes_complex::solve_waveguide_complex(mesh, problem);
        }
    }
    let eps: Vec<EpsCell> =
        tens.iter().map(|t| EpsCell { xx: t[0][0], xy: t[0][1], yy: t[1][1], zz: t[2][2] }).collect();
    let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let modes = solve_modes_tensor(mesh, freq, &eps, &mur, n_target, evanescent, n_modes, problem.order)?;
    if modes.complex_wanted {
        // a leaky or complex-propagation-constant mode is among the nearest eigenvalues: the complex solver represents it (complex vectors)
        femcore::progress::note(
            "a complex eigenvalue pair of the real mode pencil is among the wanted modes: solving again in complex arithmetic",
        );
        return crate::modes_complex::solve_waveguide_complex(mesh, problem);
    }
    let mut solutions = vec![];
    let mut scalars = vec![];
    for j in 0..modes.neff.len() {
        let (re, im) = modes.neff[j];
        let txt = if im.abs() > 1e-9 * re.abs().max(1.0) { format!("{re:.6} {im:+.4}i") } else { format!("{re:.6}") };
        let (fields, cell_fields) = mode_fields_cells(mesh, &modes, j, freq, problem.order);
        solutions.push(Solution {
            label: format!("mode {}: neff = {txt}", j + 1),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            eigenvalue: Some([re, im]),
            ..Default::default()
        });
        scalars.push((format!("neff {}", j + 1), Scalar::Text(txt)));
    }
    let mut log = vec![];
    if !modes.converged {
        log.push("warning: the eigensolver did not converge every mode to the requested tolerance".into());
    }
    if modes.neff.len() < n_modes {
        log.push(format!("only {} of {} real-valued modes were found near the target", modes.neff.len(), n_modes));
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

/// Number of edge-element dofs plus phi dofs (for reporting).
pub fn n_unknowns(mesh: &Discretisation, order: usize) -> usize {
    cell_data(mesh, order).net + mesh.n_dofs(order)
}
