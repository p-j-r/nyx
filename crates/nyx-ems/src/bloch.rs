// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Fixed-frequency 3D Floquet-Bloch waveguide and periodic grating eigenmodes (`Study::Bloch`) on tetrahedral Nedelec
//! elements (`Ned3`), with optional transverse perfectly matched layers (`x`/`z` PML), transverse Floquet periodicity,
//! midplane `symmetry`, magnetic walls (`pmc`), and dispersive single-solve Poynting group index `n_g`.
//!
//! Writing the physical electric field as `E(x, y, z) = e_tilde(x, y, z) exp(-j k_y y)` with `e_tilde` periodic along `+y`
//! over the unit-cell pitch `Lambda` (`phase = 0` on the `y`-periodic pair) and testing against `N_a exp(+j k_y y)` yields
//! the dimensionless complex quadratic eigenproblem `(A0 + nu A1 + nu^2 A2) u = 0` in the Bloch index `nu = k_y / k_0 = n_B`,
//! solved directly near `sigma = n_target + (m - 1) * (lambda_0 / Lambda)` via `eigs_quadratic_complex`.

use crate::cavity::{
    assemble_complex, basis3, compute_3d_mode_metrics, ned_relations, vertex_keys, Ned3, C3, CURVED_QUAD_DEGREE,
};
use crate::ibc::{assemble_ibc_and_jj_3d, MU0};
use crate::modes::C0;
use crate::problem::{symmetry_facet_flags, Problem, Study};
use femcore::fem::{geometry, quad_rule, quad_rule_collapsed, CurvedTet};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCsr, C};
use femcore::linalg::driver::eigs_quadratic_complex;
use femcore::pml::cell_stretches;
use femcore::props::PmlDir;
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use std::f64::consts::PI;
use std::time::Instant;

fn lap(t: &mut Instant, what: &str) {
    femcore::progress::lap(what, t.elapsed().as_secs_f64());
    *t = Instant::now();
}

/// Conversion factor `20 / (100 * ln 10)` turning `-Im(k_y)` (`1/m` amplitude decay) into power attenuation in `dB/cm`.
const ALPHA_DB_PER_CM_FROM_NEG_IM_KY: f64 = 0.086_858_896_380_650_36;

/// Stage 1 of the 3D fixed-frequency Bloch kernel: Nedelec numbering, reduced quadratic pencil `(A0, A1, A2)` in `nu = k_y / k_0`,
/// unit-cell pitch `period_y`, reference origin `y_min`, and unstretched material data for Poynting `v_g`/`n_g` and energy shares.
pub struct BlochAssembly {
    pub order: usize,
    pub n: usize,
    pub nred: usize,
    pub ned: Ned3,
    pub red: femcore::periodic::Reduction,
    pub a0_red: CCsr,
    pub a1_red: CCsr,
    pub a2_red: CCsr,
    pub freq: f64,
    pub k0: f64,
    pub lambda0: f64,
    pub period_y: f64,
    pub y_min: f64,
    pub n_target: f64,
    pub brillouin_zone: i64,
    pub sigma: f64,
    pub n_modes: usize,
    pub mu_inv: Vec<[C; 3]>,
    pub pml_cells: Vec<bool>,
    pub eps_phys_re: Vec<[[f64; 3]; 3]>,
    pub d_omega_eps: Vec<[[f64; 3]; 3]>,
    pub mur_phys: Vec<f64>,
}

/// Stage 2 output: converged complex Bloch indices `(n_B, vec_index)` sorted by distance to `sigma`.
pub struct BlochSolve {
    pub modes: Vec<(C, usize)>,
    pub vectors: Vec<Vec<C>>,
    pub residuals: Vec<f64>,
    pub dim: usize,
    pub converged: bool,
}

/// Solves a 3D `Study::Bloch` problem through its three stages (`assemble_bloch`, `solve_bloch_eigen`, `bloch_fields_cells`).
pub fn solve_bloch(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_bloch(mesh, problem)?;
    let solved = solve_bloch_eigen(&asm)?;
    let mut r = bloch_fields_cells(mesh, problem, &asm, &solved)?;
    r.seconds = t0.elapsed().as_secs_f64();
    Ok(r)
}

#[allow(clippy::too_many_arguments)]
fn add_bloch_point(
    which: usize,
    nb: usize,
    ke: &mut [C],
    w: &[[f64; 3]],
    cu: &[[f64; 3]],
    wq: f64,
    inv_k0: f64,
    inv_k0_2: f64,
    eps: &C3,
    mu_inv: &[C; 3],
) {
    match which {
        0 => {
            for i in 0..nb {
                let ci_ = cu[i];
                let wi = w[i];
                for j in 0..nb {
                    let cj = cu[j];
                    let wj = w[j];
                    let k_curl =
                        mu_inv[0] * (ci_[0] * cj[0]) + mu_inv[1] * (ci_[1] * cj[1]) + mu_inv[2] * (ci_[2] * cj[2]);
                    let mut m_eps = c(0.0, 0.0);
                    for a in 0..3 {
                        for b in 0..3 {
                            m_eps += eps[a][b] * (wi[a] * wj[b]);
                        }
                    }
                    ke[i * nb + j] += (k_curl * inv_k0_2 - m_eps) * wq;
                }
            }
        }
        1 => {
            // A1 = (j / k0) int [ (y_hat x N_i) . mu_r^-1 (curl N_j) - (curl N_i) . mu_r^-1 (y_hat x N_j) ] dOmega
            // with y_hat x N = (N_z, 0, -N_x)
            let pref = c(0.0, inv_k0 * wq);
            for i in 0..nb {
                let wi = w[i];
                let ci_ = cu[i];
                for j in 0..nb {
                    let wj = w[j];
                    let cj = cu[j];
                    let term =
                        mu_inv[0] * (wi[2] * cj[0] - ci_[0] * wj[2]) + mu_inv[2] * (ci_[2] * wj[0] - wi[0] * cj[2]);
                    ke[i * nb + j] += pref * term;
                }
            }
        }
        _ => {
            // A2 = int (y_hat x N_i) . mu_r^-1 (y_hat x N_j) dOmega
            for i in 0..nb {
                let wi = w[i];
                for j in 0..nb {
                    let wj = w[j];
                    let term = mu_inv[0] * (wi[2] * wj[2]) + mu_inv[2] * (wi[0] * wj[0]);
                    ke[i * nb + j] += term * wq;
                }
            }
        }
    }
}

/// Assembles `(A0, A1, A2)` over all full Nedelec DOFs on straight and curved tetrahedra.
fn assemble_bloch_matrices(
    mesh: &Discretisation,
    ned: &Ned3,
    k0: f64,
    eps: &[C3],
    mu_inv: &[[C; 3]],
) -> (CCsr, CCsr, CCsr) {
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    let curved_cell: Vec<bool> = (0..mesh.n_cells()).map(|ci| mesh.cell_is_curved(ci)).collect();
    let cquad = quad_rule_collapsed(CURVED_QUAD_DEGREE);
    let nb = ned.nb();
    let n = ned.n_dofs();
    let nc = mesh.n_cells();
    let inv_k0 = 1.0 / k0;
    let inv_k0_2 = inv_k0 * inv_k0;
    let run = |which: usize| -> CCsr {
        assemble_complex(nc, n, |range, coo| {
            let (mut w, mut cu) = (Vec::new(), Vec::new());
            let mut ke = vec![c(0.0, 0.0); nb * nb];
            let mut gs = [[0.0; 3]; 4];
            for ci in range {
                let gid = ned.keys(mesh, ci);
                for b in 0..4 {
                    for k in 0..3 {
                        gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                    }
                }
                ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
                if curved_cell[ci] {
                    let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
                    for q in 0..cquad.w.len() {
                        let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                        let (gq, det) = ct.jacobian(&lam);
                        basis3(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                        let wq = cquad.w[q] * det.abs() / 6.0;
                        add_bloch_point(which, nb, &mut ke, &w, &cu, wq, inv_k0, inv_k0_2, &eps[ci], &mu_inv[ci]);
                    }
                } else {
                    for q in 0..quad.w.len() {
                        let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                        basis3(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                        let wq = quad.w[q] * geo.vol[ci];
                        add_bloch_point(which, nb, &mut ke, &w, &cu, wq, inv_k0, inv_k0_2, &eps[ci], &mu_inv[ci]);
                    }
                }
                coo.push_block(&ned.cell_dofs[ci], &ned.cell_dofs[ci], &ke);
            }
        })
    };
    (run(0), run(1), run(2))
}

/// Stage 1: validates `y`-periodicity and transverse PMLs, assembles `(A0, A1, A2)` on `Ned3`, and reduces them via `P^H (.) P`
/// with `phase = 0` forced on the `y`-periodic pair (periodic envelope `e_tilde`) and user `phase` preserved on transverse pairs.
pub fn assemble_bloch(mesh: &Discretisation, problem: &Problem) -> Result<BlochAssembly, String> {
    if mesh.dim != 3 {
        return Err("the bloch study of EM waves needs a 3D mesh".into());
    }
    let Study::Bloch { freq, n_target, n_modes, brillouin_zone } = problem.study.clone() else {
        return Err("assemble_bloch requires a Bloch study".into());
    };
    if !(freq > 0.0 && freq.is_finite()) {
        return Err("bloch frequency must be positive and finite".into());
    }
    if n_modes < 1 {
        return Err("at least one mode is needed".into());
    }
    if let Some((dom, _)) = problem.materials.pml.iter().find(|(_, p)| p.y != PmlDir::Off) {
        return Err(format!(
            "the bloch study is periodic along y and cannot use a y-absorbing PML (domain {dom}); use x or z for transverse radiation"
        ));
    }

    let mut matches = Vec::with_capacity(problem.periodic.len());
    let mut is_y_pair = Vec::with_capacity(problem.periodic.len());
    let mut period_y: Option<f64> = None;
    for p in &problem.periodic {
        let mt = femcore::periodic::match_vertices(mesh, p)?;
        let sn = (mt.shift[0] * mt.shift[0] + mt.shift[1] * mt.shift[1] + mt.shift[2] * mt.shift[2]).sqrt();
        if sn <= 0.0 {
            return Err("periodic pair has zero translation vector".into());
        }
        if mt.shift[1].abs() > 0.99 * sn {
            let ly = mt.shift[1].abs();
            if let Some(prev) = period_y {
                if (ly - prev).abs() > 1e-6 * prev.max(1e-30) {
                    return Err(format!("conflicting y-periodic pitches in bloch study ({prev:.6e} m vs {ly:.6e} m)"));
                }
            } else {
                period_y = Some(ly);
            }
            is_y_pair.push(true);
        } else {
            if mt.shift[1].abs() > 1e-6 * sn {
                return Err(
                    "transverse periodic pairs in the bloch study must have zero translation along the y propagation axis"
                        .into(),
                );
            }
            is_y_pair.push(false);
        }
        matches.push(mt);
    }
    let period_y = period_y.ok_or_else(|| {
        "the bloch study requires a periodic boundary condition along the y propagation axis".to_string()
    })?;

    let order = if problem.order == 1 { 1 } else { 2 };
    let tens = problem.materials.cell_complex_permittivity(mesh, Some(1.0))?;
    let eps_phys_re: Vec<[[f64; 3]; 3]> = tens.iter().map(|(re, _)| *re).collect();
    let d_omega_eps = problem.cell_d_omega_eps(mesh)?;
    let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let has_pml = problem.materials.has_pml();
    let pml_cells: Vec<bool> = (0..mesh.n_cells())
        .map(|ci| {
            problem
                .materials
                .pml
                .get(&mesh.cell_labels[ci])
                .is_some_and(|p| p.x != PmlDir::Off || p.y != PmlDir::Off || p.z != PmlDir::Off)
        })
        .collect();
    let stretch = if has_pml { Some(cell_stretches(mesh, &problem.materials.pml)) } else { None };
    let mut eps: Vec<C3> = Vec::with_capacity(mesh.n_cells());
    let mut mu_inv: Vec<[C; 3]> = Vec::with_capacity(mesh.n_cells());
    for ci in 0..mesh.n_cells() {
        let (sx, sy, sz) =
            stretch.as_ref().map_or((c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1], s[ci][2]));
        let lam = [sy * sz / sx, sx * sz / sy, sx * sy / sz];
        let pml_here = stretch.is_some()
            && (sx.re != 1.0 || sx.im != 0.0 || sy.re != 1.0 || sy.im != 0.0 || sz.re != 1.0 || sz.im != 0.0);
        let (t_re, t_im) = tens[ci];
        let (sv, det) = ([sx, sy, sz], sx * sy * sz);
        let mut e = [[c(0.0, 0.0); 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let ec = c(t_re[a][b], -t_im[a][b]);
                e[a][b] = if pml_here { det * ec / (sv[a] * sv[b]) } else { ec };
            }
        }
        eps.push(e);
        let one = c(1.0, 0.0);
        let mu0 = c(1.0 / mur[ci], 0.0);
        mu_inv.push([
            mu0 / if pml_here { lam[0] } else { one },
            mu0 / if pml_here { lam[1] } else { one },
            mu0 / if pml_here { lam[2] } else { one },
        ]);
    }

    let mut tl = Instant::now();
    let vkey = vertex_keys(mesh, &matches);
    let ned = Ned3::new(mesh, order, vkey);
    let n = ned.n_dofs();
    let k0 = 2.0 * PI * freq / C0;
    let lambda0 = C0 / freq;
    let (mut a0_full, a1_full, a2_full) = assemble_bloch_matrices(mesh, &ned, k0, &eps, &mu_inv);
    lap(&mut tl, "bloch3d: numbering and assembly");

    let has_sheets = !problem.pec.is_empty() || !problem.surface_impedance.is_empty();
    let ibc_jj = if has_sheets {
        Some(assemble_ibc_and_jj_3d(mesh, &ned, &problem.pec, &problem.surface_impedance, &[], freq)?)
    } else {
        None
    };
    if let Some(ref ij) = ibc_jj {
        if ij.kc_ibc_full.nnz() > 0 {
            a0_full = a0_full.add_scaled(&ij.kc_ibc_full, c(1.0 / (k0 * k0), 0.0));
        }
    }

    let exterior = mesh.exterior_facets();
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let (sym_pec_facets, sym_pmc_facets) = symmetry_facet_flags(mesh, problem.symmetry)?;

    let mut fixed_set = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        let is_ext_pec = exterior[f]
            && !periodic_labels.contains(&label)
            && (sym_pec_facets[f]
                || (!problem.pmc.contains(&label) && !sym_pmc_facets[f] && !ibc_labels.contains(&label)));
        if is_ext_pec {
            let fv = mesh.facet(f);
            for d in ned.boundary_dofs(mesh, [fv[0], fv[1], fv[2]]) {
                fixed_set.insert(d);
            }
        }
    }
    if let Some(ref ij) = ibc_jj {
        for &fv in &ij.pec_tris {
            for d in ned.boundary_dofs(mesh, fv) {
                fixed_set.insert(d);
            }
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();

    let mut relations = vec![];
    for ((p, mt), &y_axis) in problem.periodic.iter().zip(&matches).zip(&is_y_pair) {
        let phase_eff = if y_axis { 0.0 } else { p.phase };
        relations.extend(ned_relations(mesh, &ned, mt, phase_eff)?);
    }
    let red = femcore::periodic::Reduction::build(n, &relations, &fixed);
    let pmat = red.prolong();
    let nred = red.n_reduced();
    if nred == 0 || n_modes >= 2 * nred {
        return Err("number of modes exceeds the reduced Bloch degree-of-freedom count".into());
    }

    let a0_red = femcore::periodic::congruence(&pmat, &a0_full);
    let a1_red = femcore::periodic::congruence(&pmat, &a1_full);
    let a2_red = femcore::periodic::congruence(&pmat, &a2_full);
    lap(&mut tl, "bloch3d: periodic congruence reduction");

    let (lo, _hi) = mesh.bbox();
    let y_min = lo[1];
    let sigma = n_target + (brillouin_zone - 1) as f64 * (lambda0 / period_y);

    Ok(BlochAssembly {
        order,
        n,
        nred,
        ned,
        red,
        a0_red,
        a1_red,
        a2_red,
        freq,
        k0,
        lambda0,
        period_y,
        y_min,
        n_target,
        brillouin_zone,
        sigma,
        n_modes,
        mu_inv,
        pml_cells,
        eps_phys_re,
        d_omega_eps,
        mur_phys: mur,
    })
}

/// Stage 2: solves `(A0 + nu A1 + nu^2 A2) u = 0` near `sigma`, filters true QEP relative residuals `< 1e-5`,
/// and sorts converged Bloch modes nearest `sigma`.
pub fn solve_bloch_eigen(asm: &BlochAssembly) -> Result<BlochSolve, String> {
    let mut tl = Instant::now();
    let sigma_c = c(asm.sigma, 0.0);
    let pairs = eigs_quadratic_complex(&asm.a0_red, &asm.a1_red, &asm.a2_red, sigma_c, asm.n_modes, 1e-9)?;
    lap(&mut tl, "bloch3d: quadratic eigenvalue iteration");

    let mut modes: Vec<(C, usize)> = Vec::new();
    for (i, &nu) in pairs.values.iter().enumerate() {
        if pairs.residuals.get(i).is_some_and(|&r| r >= 1e-5) {
            continue;
        }
        modes.push((nu, i));
    }
    modes.sort_by(|a, b| (a.0 - sigma_c).norm().total_cmp(&(b.0 - sigma_c).norm()));
    Ok(BlochSolve {
        modes,
        vectors: pairs.vectors,
        residuals: pairs.residuals,
        dim: pairs.dim,
        converged: pairs.converged,
    })
}

/// Stage 3: reconstructs physical `E = e_tilde exp(-j k_y (y - y_min))` and
/// `H = (j / (omega mu_0)) mu_r^-1 (curl e_tilde - j k_y y_hat x e_tilde) exp(-j k_y (y - y_min))`
/// at all 10 P2 nodes per tetrahedron, builds plot & cell fields, and computes `v_g`, `n_g`, attenuation `alpha`,
/// energy shares, and `FieldFrame`.
pub fn bloch_fields_cells(
    mesh: &Discretisation,
    problem: &Problem,
    asm: &BlochAssembly,
    solved: &BlochSolve,
) -> Result<SolveResult, String> {
    if solved.modes.is_empty() {
        return Err("no Bloch modes converged near the target index".into());
    }
    let geo = geometry(mesh);
    let lamn = local_bary(3, 2);
    let nloc = lamn.len();
    let nc = mesh.n_cells();
    let nb = asm.ned.nb();
    let omega = (2.0 * PI * asm.freq).max(1e-300);

    let mut solutions = Vec::with_capacity(solved.modes.len());
    let mut scalars = Vec::new();

    for (j, &(n_b, vec_idx)) in solved.modes.iter().enumerate() {
        let x = asm.red.expand(&solved.vectors[vec_idx]);
        let k_y = n_b * asm.k0;
        let j_ky = c(0.0, 1.0) * k_y;

        let mut e_phys = vec![[c(0.0, 0.0); 3]; nc * nloc];
        let mut curl_k_phys = vec![[c(0.0, 0.0); 3]; nc * nloc];
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];

        for ci in 0..nc {
            let gid = asm.ned.keys(mesh, ci);
            let v_cell = mesh.cell(ci);
            for b in 0..4 {
                for kk in 0..3 {
                    gs[b][kk] = geo.g[ci * 12 + b * 3 + kk];
                }
            }
            let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
            for (a, l) in lamn.iter().enumerate() {
                let lam = [l[0], l[1], l[2], l[3]];
                let y_node = match &curved {
                    Some(ct) => {
                        let (gq, _) = ct.jacobian(&lam);
                        basis3(asm.order, &lam, &gq, &gid, &mut w, &mut cu);
                        ct.position(&lam)[1]
                    }
                    None => {
                        basis3(asm.order, &lam, &gs, &gid, &mut w, &mut cu);
                        (0..4).map(|b| lam[b] * mesh.point(v_cell[b] as usize)[1]).sum()
                    }
                };
                let mut val = [c(0.0, 0.0); 3];
                let mut cval = [c(0.0, 0.0); 3];
                for s in 0..nb {
                    let coef = x[asm.ned.cell_dofs[ci][s]];
                    for kk in 0..3 {
                        val[kk] += coef * w[s][kk];
                        cval[kk] += coef * cu[s][kk];
                    }
                }
                let dy = y_node - asm.y_min;
                let phase_y = (-j_ky * dy).exp();
                let idx = ci * nloc + a;
                e_phys[idx] = [val[0] * phase_y, val[1] * phase_y, val[2] * phase_y];
                curl_k_phys[idx] =
                    [(cval[0] - j_ky * val[2]) * phase_y, cval[1] * phase_y, (cval[2] + j_ky * val[0]) * phase_y];
            }
        }

        let big =
            e_phys
                .iter()
                .flat_map(|v| v.iter())
                .cloned()
                .fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let mag: Vec<f64> =
            e_phys.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
        let mag_xz: Vec<f64> = e_phys.iter().map(|v| (v[0].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
        let mag_y: Vec<f64> = e_phys.iter().map(|v| v[1].norm()).collect();
        let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
        let avg = |vals: Vec<f64>| node_average(mesh, &vals.iter().map(|v| v * s).collect::<Vec<f64>>(), 1, 2);

        let fields = vec![
            Field::scalar("|E|", "a.u.", avg(mag)),
            Field::scalar("|Exz|", "a.u.", avg(mag_xz)),
            Field::scalar("|Ey|", "a.u.", avg(mag_y)),
            Field::scalar("Re Ex", "a.u.", avg(e_phys.iter().map(|v| (v[0] * rot).re).collect())),
            Field::scalar("Re Ey", "a.u.", avg(e_phys.iter().map(|v| (v[1] * rot).re).collect())),
            Field::scalar("Re Ez", "a.u.", avg(e_phys.iter().map(|v| (v[2] * rot).re).collect())),
        ];

        let e_scale = rot * s;
        let mut cells = CellFields::new(mesh, 2);
        for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
            let z: Vec<C> = e_phys.iter().map(|v| v[k] * e_scale).collect();
            cells.push_complex(*name, "a.u.", z.iter().map(|v| v.re).collect(), z.iter().map(|v| v.im).collect());
        }

        let h_prefactor = c(0.0, 1.0 / (omega * MU0)) * e_scale;
        let mut h_phys = vec![[c(0.0, 0.0); 3]; nc * nloc];
        for ci in 0..nc {
            let mui = c(1.0 / asm.mur_phys[ci], 0.0);
            for a in 0..nloc {
                let idx = ci * nloc + a;
                let ck = curl_k_phys[idx];
                h_phys[idx] = [ck[0] * mui * h_prefactor, ck[1] * mui * h_prefactor, ck[2] * mui * h_prefactor];
            }
        }
        for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
            let z: Vec<C> = h_phys.iter().map(|v| v[k]).collect();
            cells.push_complex(*name, "A/m", z.iter().map(|v| v.re).collect(), z.iter().map(|v| v.im).collect());
        }

        let m_metrics = compute_3d_mode_metrics(
            mesh,
            &asm.ned,
            &geo,
            &x,
            k_y,
            omega,
            &asm.pml_cells,
            &asm.eps_phys_re,
            &asm.d_omega_eps,
            &asm.mur_phys,
            true,
        );

        let mode_num = j + 1;
        let lossy = n_b.im.abs() > 1e-12 * n_b.re.abs().max(1.0);
        let label = if lossy {
            format!("mode {mode_num}: n_B = {:.6} {:+.3e}i", n_b.re, n_b.im)
        } else {
            format!("mode {mode_num}: n_B = {:.6}", n_b.re)
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(n_b.re),
            eigenvalue: Some([n_b.re, n_b.im]),
        });

        scalars.push((format!("n_B re {mode_num}"), Scalar::Num(n_b.re)));
        scalars.push((format!("n_B im {mode_num}"), Scalar::Num(n_b.im)));
        if lossy {
            scalars.push((format!("n_B {mode_num}"), Scalar::Text(format!("{:.8} {:+.4e}i", n_b.re, n_b.im))));
        } else {
            scalars.push((format!("n_B {mode_num}"), Scalar::Num(n_b.re)));
        }
        let alpha_db_cm = ALPHA_DB_PER_CM_FROM_NEG_IM_KY * (-k_y.im);
        scalars.push((format!("alpha {mode_num} (dB/cm)"), Scalar::Num(alpha_db_cm)));
        scalars.push((format!("k_y re {mode_num} (rad/m)"), Scalar::Num(k_y.re)));
        scalars.push((format!("k_y im {mode_num} (rad/m)"), Scalar::Num(k_y.im)));
        scalars.push((format!("phase {mode_num} (rad)"), Scalar::Num(k_y.re * asm.period_y)));
        scalars.push((format!("v_g {mode_num} (m/s)"), Scalar::Num(m_metrics.vg)));
        if m_metrics.vg.abs() > 1e-9 * C0 {
            scalars.push((format!("n_g {mode_num}"), Scalar::Num(C0 / m_metrics.vg)));
        }
        scalars.push((format!("Ex share {mode_num}"), Scalar::Num(m_metrics.comp_shares[0])));
        scalars.push((format!("Ey share {mode_num}"), Scalar::Num(m_metrics.comp_shares[1])));
        scalars.push((format!("Ez share {mode_num}"), Scalar::Num(m_metrics.comp_shares[2])));
        for &(dom_label, share) in &m_metrics.domain_shares {
            scalars.push((format!("share in domain {dom_label} {mode_num}"), Scalar::Num(share)));
        }
        scalars.push((format!("layer fraction {mode_num}"), Scalar::Num(m_metrics.layer_fraction)));
    }

    scalars.push(
        nyx_physics::FieldFrame {
            lateral: 0,
            normal: 2,
            propagation: 1,
            left_handed: true,
            outside: crate::optical::absorbing_domains(problem),
        }
        .scalar(),
    );

    let log = vec![
        format!(
            "3D Bloch Nedelec order {}, {} dofs ({} after constraints), f = {:.4e} Hz, pitch = {:.4e} m, shift n_B = {:.6} (BZ {}), Krylov dim {}{}",
            asm.order,
            asm.n,
            asm.nred,
            asm.freq,
            asm.period_y,
            asm.sigma,
            asm.brillouin_zone,
            solved.dim,
            if solved.converged { "" } else { " (not every pair converged; those are omitted)" }
        ),
        "polarization axes: lateral x, normal z, propagation y".into(),
    ];

    Ok(SolveResult {
        physics: "ems".into(),
        study: "bloch".into(),
        order: 2,
        dim: 3,
        solutions,
        scalars,
        log,
        dofs: asm.n,
        seconds: 0.0,
    })
}
