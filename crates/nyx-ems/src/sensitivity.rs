// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Discrete adjoint and modal eigenvalue sensitivities for `nyx-ems` (ADR 142):
//!
//! * `modes` / `modes_complex`: 2D waveguide effective-index `n_eff` and modal `Q` sensitivities via
//!   `femcore::linalg::sens::{eigenvalue_derivative, nonlinear_eigenvalue_derivative, left_eigenvector,
//!   match_modes}` across lossless, PML, tilted-crystal (`a1` quadratic coupling), gyrotropic (`FR-34`),
//!   and 2D Hall sheet (`FR-35`) cross-sections.
//! * `cavity`: 3D resonant frequency `f_0` and `Q` sensitivities on straight and curved (`ADR 82`)
//!   Nedelec tetrahedra, including periodic/PEC reductions and non-symmetric media.
//! * `bloch`: 3D Floquet-Bloch quadratic pencil `(A0 + nu A1 + nu^2 A2) u = 0` derivatives via
//!   `NonlinearEigenProblem::polynomial` and `nonlinear_eigenvalue_derivative`.
//! * `driven`: multi-port scattering (`|S_{pq}|^2`, `SplitterBalance`) and domain `FieldPower`
//!   discrete adjoint gradients, reusing the cached Stage-2 sparse factorization in full-order HDM
//!   mode or solving the `r x r` reduced transpose system in Hermite Tangential Galerkin PROM mode.

use crate::bloch::{assemble_bloch, solve_bloch_eigen};
use crate::cavity::{
    assemble_emw_eigen, assemble_maxwell_with_pml, basis3_keyed, safe_curved_jacobian, solve_emw_modes, Ned3,
};
use crate::driven::{assemble_emw_driven, DrivenAssembly};
use crate::mixed_ops::{simplex_cell_materials_3d, EmPmlTable};
use crate::modes::C0;
use crate::modes_complex::{
    assemble_modes_complex, assemble_modes_complex_quad_pml, assemble_modes_complex_with_bcs_and_symmetry_and_sheets,
    solve_mode_pairs_complex, ModeAssemblyC,
};
use crate::problem::{decode, Problem, Study, Symmetry};
use crate::Ems;
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{geometry, quad_rule, quad_rule_collapsed, shape, CurvedTet};
use femcore::linalg::csr::{c, CCsr, C};
use femcore::linalg::equilibrate::scale_symmetric;
use femcore::linalg::nep::NonlinearEigenProblem;
use femcore::linalg::prom::{
    run_greedy_prom_sweep, solve_dense_complex_multi, ImageProm, PromLinearSolver, PromOperator, SparseDrivenSolver,
};
use femcore::linalg::sens::{eigenvalue_derivative, left_eigenvector, match_modes, nonlinear_eigenvalue_derivative};
use femcore::optical::{cell_materials, EpsC};
use femcore::topology::{displace_simplex_mesh, optical_eps_interpolate, OpticalInterpMode};
use nyx_physics::{
    check_mesh, discretise, validate_spec, Ctx, DesignParam, Discretisation, Error, EvaluatedTopologyParam, Mesh,
    Physics, ProblemSpec, SensitivityObjective, SensitivityReport, SensitivityRequest, TrackedModeBasis,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

type C3 = [[C; 3]; 3];

/// Evaluates `Ems::sensitivity` on `mesh`.
pub fn compute_sensitivity(
    physics: &Ems,
    mesh: &Mesh,
    spec: &ProblemSpec,
    ctx: &mut Ctx,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    req.validate_finite()?;
    validate_spec(physics, spec)?;
    physics.validate(spec, mesh)?;
    check_mesh(mesh)?;
    ctx.check()?;

    let mut base_problem = decode(spec, ctx)?;
    base_problem.apply_transferred_inputs(mesh, ctx)?;

    match base_problem.study {
        Study::Modes { .. } => waveguide_sensitivity(mesh, ctx, &base_problem, req),
        Study::Cavity { .. } => cavity_sensitivity(mesh, ctx, &base_problem, req),
        Study::Bloch { .. } => bloch_sensitivity(mesh, ctx, &base_problem, req),
        Study::Driven { .. } => driven_sensitivity(mesh, ctx, &base_problem, req),
        Study::Eme { .. } => eme_sensitivity(mesh, ctx, &base_problem, req),
        Study::AxisymmetricModes { .. }
        | Study::NonlinearDriven(_)
        | Study::NonlinearModes(_)
        | Study::Dgtd(_)
        | Study::Fdtd(..) => Err(Error::Unsupported(
            "sensitivity is supported on linear modes, cavity, bloch, driven, and eme studies".into(),
        )),
    }
}

fn build_topology_states(
    simplex: &Discretisation,
    ctx: &Ctx,
    req: &SensitivityRequest,
) -> Result<BTreeMap<usize, EvaluatedTopologyParam>, Error> {
    let mut out = BTreeMap::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity { domain, filter_radius, beta, eta, initial_density, densities, .. } = param
        {
            let dom_labels = ctx.resolver.domains(domain)?;
            if dom_labels.is_empty() {
                return Err(Error::invalid(format!("TopologyDensity domain '{domain}' resolved to zero mesh domains")));
            }
            let cached = req.topology_filters.get(&idx);
            let st = EvaluatedTopologyParam::build(
                simplex,
                &dom_labels,
                *filter_radius,
                *beta,
                *eta,
                *initial_density,
                densities,
                cached,
            )?;
            out.insert(idx, st);
        }
    }
    Ok(out)
}

#[inline]
fn optical_interp_of_param(
    rho_bar: f64,
    simp_penalty: f64,
    ramp_q: Option<f64>,
    val_min: f64,
    val_max: f64,
    val_min_im: f64,
    val_max_im: f64,
    plasmon_damping: f64,
) -> (C, C) {
    let eps_void = c(val_min, val_min_im);
    let eps_solid = c(val_max, val_max_im);
    let mode = match ramp_q {
        Some(q) => OpticalInterpMode::DielectricRamp { q },
        None => OpticalInterpMode::Auto { penalty: simp_penalty, damping_floor: plasmon_damping },
    };
    optical_eps_interpolate(rho_bar, eps_void, eps_solid, mode)
}

fn target_mode_index(obj: &SensitivityObjective) -> usize {
    match obj {
        SensitivityObjective::Eigenfrequency { mode, .. }
        | SensitivityObjective::EffectiveIndex { mode, .. }
        | SensitivityObjective::QualityFactor { mode, .. } => *mode,
        _ => 0,
    }
}

fn eval_neff_or_q_objective(obj: &SensitivityObjective, neff: C) -> Result<(f64, C, BTreeMap<String, f64>), Error> {
    let mut scalars = BTreeMap::new();
    scalars.insert("neff_re".into(), neff.re);
    scalars.insert("neff_im".into(), neff.im);
    match obj {
        SensitivityObjective::EffectiveIndex { target_neff, .. } => match target_neff {
            None => Ok((neff.re, c(1.0, 0.0), scalars)),
            Some(nt) => {
                let diff = neff.re - nt;
                Ok((0.5 * diff * diff, c(diff, 0.0), scalars))
            }
        },
        SensitivityObjective::QualityFactor { inverse, .. } => {
            let re = neff.re;
            let im = neff.im.abs();
            let sgn_im = if neff.im >= 0.0 { 1.0 } else { -1.0 };
            let sgn_re = if re >= 0.0 { 1.0 } else { -1.0 };
            if *inverse {
                let inv_q = 2.0 * im / re.abs().max(1e-30);
                scalars.insert("inv_q".into(), inv_q);
                let d_re = -2.0 * im * sgn_re / (re * re).max(1e-30);
                let d_im = 2.0 * sgn_im / re.abs().max(1e-30);
                Ok((inv_q, c(d_re, d_im), scalars))
            } else {
                if im <= 1e-30 * re.abs() {
                    return Err(Error::invalid(
                        "QualityFactor sensitivity requires a lossy/PML mode with non-zero Im(neff) (or set inverse = true)",
                    ));
                }
                let q = re.abs() / (2.0 * im);
                scalars.insert("Q".into(), q);
                let d_re = sgn_re / (2.0 * im);
                let d_im = -re.abs() * sgn_im / (2.0 * im * im);
                Ok((q, c(d_re, d_im), scalars))
            }
        }
        other => Err(Error::Unsupported(format!(
            "ems waveguide/bloch sensitivity supports EffectiveIndex and QualityFactor, got '{other:?}'"
        ))),
    }
}

fn eval_freq_or_q_objective(obj: &SensitivityObjective, f_hz: C) -> Result<(f64, C, BTreeMap<String, f64>), Error> {
    let mut scalars = BTreeMap::new();
    scalars.insert("f_re_hz".into(), f_hz.re);
    scalars.insert("f_im_hz".into(), f_hz.im);
    match obj {
        SensitivityObjective::Eigenfrequency { target_hz, scale_hz, .. } => match target_hz {
            None => Ok((f_hz.re, c(1.0, 0.0), scalars)),
            Some(ft) => {
                let sc = scale_hz.unwrap_or_else(|| ft.abs().max(1.0));
                let diff = (f_hz.re - ft) / sc;
                Ok((0.5 * diff * diff, c(diff / sc, 0.0), scalars))
            }
        },
        SensitivityObjective::QualityFactor { inverse, .. } => {
            let re = f_hz.re;
            let im = f_hz.im.abs();
            let sgn_im = if f_hz.im >= 0.0 { 1.0 } else { -1.0 };
            let sgn_re = if re >= 0.0 { 1.0 } else { -1.0 };
            if *inverse {
                let inv_q = 2.0 * im / re.abs().max(1e-30);
                scalars.insert("inv_q".into(), inv_q);
                let d_re = -2.0 * im * sgn_re / (re * re).max(1e-30);
                let d_im = 2.0 * sgn_im / re.abs().max(1e-30);
                Ok((inv_q, c(d_re, d_im), scalars))
            } else {
                if im <= 1e-30 * re.abs() {
                    return Err(Error::invalid(
                        "QualityFactor sensitivity requires a lossy/PML cavity mode with non-zero Im(f) (or set inverse = true)",
                    ));
                }
                let q = re.abs() / (2.0 * im);
                scalars.insert("Q".into(), q);
                let d_re = sgn_re / (2.0 * im);
                let d_im = -re.abs() * sgn_im / (2.0 * im * im);
                Ok((q, c(d_re, d_im), scalars))
            }
        }
        other => Err(Error::Unsupported(format!(
            "ems cavity sensitivity supports Eigenfrequency and QualityFactor, got '{other:?}'"
        ))),
    }
}

// ---------------------------------------------------------------------------
// 1. 2D Waveguide Modes (`modes` / `modes_complex`)
// ---------------------------------------------------------------------------

fn assemble_2d_waveguide_for_sens(
    mesh: &Discretisation,
    problem: &Problem,
    freq: f64,
    eps_overrides: &[(usize, C)],
) -> Result<ModeAssemblyC, String> {
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

    if (needs_quad_pml || problem.has_cell_material_perturbation()) && eps_overrides.is_empty() {
        return assemble_modes_complex_quad_pml(mesh, problem, freq, &pml_table);
    }

    let (mut eps, mu) = cell_materials(mesh, &problem.materials)?;
    for &(ci, ec) in eps_overrides {
        eps[ci] = EpsC {
            xx: ec,
            xy: c(0.0, 0.0),
            yx: c(0.0, 0.0),
            yy: ec,
            zz: ec,
            xz: c(0.0, 0.0),
            zx: c(0.0, 0.0),
            yz: c(0.0, 0.0),
            zy: c(0.0, 0.0),
        };
    }

    if problem.pec.is_empty()
        && problem.surface_impedance.is_empty()
        && problem.sheet_2d.is_empty()
        && problem.pmc.is_empty()
        && problem.symmetry == Symmetry::None
    {
        assemble_modes_complex(mesh, freq, &eps, &mu, problem.order)
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
        )
    }
}

/// Elementwise 2D waveguide `y_c^H (dA_c - lambda dB_c) x_c / d(eps_r)` for isotropic permittivity:
/// `A_c` contributes `-k0^2 int (w_y^* . w_x) dA`; `B_c` (after `1/k0` scaling on `phi` block!)
/// contributes `-int (phi_y_scaled^* phi_x_scaled) dA`, so `dA_c - lambda dB_c` per unit `eps_r` is
/// `-k0^2 int (Et_left^* . Et_right) dA + lambda * k0^2 int (phi_left_unscaled^* phi_right_unscaled) dA`!
/// Notice in `solve_mode_pairs_complex`: `phi[j]` = `x[net..] / k0`, so `x_free[phi_dof] = k0 * phi[dof]`,
/// meaning `(x_scaled)^* (-1) (x_scaled) = -k0^2 phi_left^* phi_right`, and `-lambda * (-k0^2 phi_left^* phi_right) = +lambda k0^2 phi_left^* phi_right`!
fn cell_2d_waveguide_eps_dots(
    mesh: &Discretisation,
    order: usize,
    k0: f64,
    lambda: C,
    et_left: &[C],
    phi_left: &[C],
    et_right: &[C],
    phi_right: &[C],
    cells: &[usize],
) -> Vec<C> {
    let geo = geometry(mesh);
    let cd = cell_data(mesh, order);
    let nl = mesh.nloc(order);
    let nloc_e = if order == 1 { 3 } else { 8 };
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, order)).collect();
    let k2 = k0 * k0;
    const CHUNK: usize = 256;
    let parts = femcore::par::map_chunks(cells, CHUNK, |slice| {
        let mut out = Vec::with_capacity(slice.len());
        let mut w = vec![[0.0; 2]; nloc_e];
        let mut curl = vec![0.0; nloc_e];
        let mut pd = Vec::with_capacity(nl);
        for &ci in slice {
            let g = tri_g(&geo, ci);
            let area = geo.vol[ci];
            let rows_e = &cd.gdofs[ci];
            mesh.cell_dofs(ci, order, &mut pd);
            let mut int_et = c(0.0, 0.0);
            let mut int_phi = c(0.0, 0.0);
            for (q, l) in quad.lam.iter().enumerate() {
                edge_basis(l, &g, &cd.sign[ci], order, &mut w, &mut curl);
                let (nv_, _) = &tab[q];
                let mut ey_q = [c(0.0, 0.0); 2];
                let mut ex_q = [c(0.0, 0.0); 2];
                for s in 0..nloc_e {
                    let gd = rows_e[s];
                    let yl = et_left[gd];
                    let xr = et_right[gd];
                    ey_q[0] += yl * w[s][0];
                    ey_q[1] += yl * w[s][1];
                    ex_q[0] += xr * w[s][0];
                    ex_q[1] += xr * w[s][1];
                }
                let mut py_q = c(0.0, 0.0);
                let mut px_q = c(0.0, 0.0);
                for a in 0..nl {
                    let gd = pd[a];
                    py_q += phi_left[gd] * nv_[a];
                    px_q += phi_right[gd] * nv_[a];
                }
                let wq = quad.w[q] * area;
                int_et += (ey_q[0].conj() * ex_q[0] + ey_q[1].conj() * ex_q[1]) * wq;
                int_phi += (py_q.conj() * px_q) * wq;
            }
            out.push(int_et * (-k2) + int_phi * (lambda * k2));
        }
        out
    });
    parts.into_iter().flatten().collect()
}

fn waveguide_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let Study::Modes { freq, n_target, evanescent, n_modes } = base_problem.study else {
        return Err(Error::invalid("expected Modes study"));
    };
    let simplex = discretise(mesh)?;
    if simplex.dim != 2 {
        return Err(Error::invalid("2D waveguide sensitivity requires a 2D mesh"));
    }
    let topo_states = build_topology_states(&simplex, ctx, req)?;

    let mut eps_overrides: Vec<(usize, C)> = Vec::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity {
            simp_penalty,
            ramp_q,
            val_min,
            val_max,
            val_min_im,
            val_max_im,
            plasmon_damping,
            ..
        } = param
        {
            let st = &topo_states[&idx];
            for (&ci, &rb) in st.filter.design_cells().iter().zip(&st.rho_bar) {
                let (ec, _) = optical_interp_of_param(
                    rb,
                    *simp_penalty,
                    *ramp_q,
                    *val_min,
                    *val_max,
                    *val_min_im,
                    *val_max_im,
                    *plasmon_damping,
                );
                eps_overrides.push((ci, ec));
            }
        }
    }

    let asm =
        assemble_2d_waveguide_for_sens(&simplex, base_problem, freq, &eps_overrides).map_err(Error::from_kernel)?;
    let mode_idx = target_mode_index(&req.objective);
    let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
    let ask = n_modes.max(mode_idx + 1);
    let modes = solve_mode_pairs_complex(&asm, target, ask).map_err(Error::from_kernel)?;
    if modes.neff.is_empty() {
        return Err(Error::Solver("no converged 2D waveguide modes found".into()));
    }

    // Reconstruct reduced vectors `x_free` and `y_free` (`phi` block multiplied back by `k0`)
    let k0 = asm.k0;
    let extract_free = |et: &[C], phi: &[C]| -> Vec<C> {
        asm.free.iter().map(|&dof| if dof < asm.net { et[dof] } else { phi[dof - asm.net] * k0 }).collect()
    };
    let mut neffs = modes.neff.clone();
    let mut x_vecs: Vec<Vec<C>> = (0..neffs.len()).map(|i| extract_free(&modes.et[i], &modes.phi[i])).collect();
    let mut y_vecs: Vec<Vec<C>> =
        (0..neffs.len()).map(|i| extract_free(&modes.et_left[i], &modes.phi_left[i])).collect();
    let mut et_r = modes.et.clone();
    let mut phi_r = modes.phi.clone();
    let mut et_l = modes.et_left.clone();
    let mut phi_l = modes.phi_left.clone();

    if let Some(ref ref_basis) = req.reference_modes {
        if ref_basis.eigenvectors.len() == x_vecs.len()
            && ref_basis.eigenvectors.first().is_some_and(|v| v.len() == asm.free.len())
        {
            let matched = match_modes(&ref_basis.eigenvectors, &x_vecs, Some(&asm.b), 0.1, 0.05);
            let (old_n, old_x, old_y) = (neffs.clone(), x_vecs.clone(), y_vecs.clone());
            let (old_etr, old_phir, old_etl, old_phil) = (et_r.clone(), phi_r.clone(), et_l.clone(), phi_l.clone());
            for (ok, &nj) in matched.assignment.iter().enumerate() {
                neffs[ok] = old_n[nj];
                x_vecs[ok] = old_x[nj].clone();
                y_vecs[ok] = old_y[nj].clone();
                et_r[ok] = old_etr[nj].clone();
                phi_r[ok] = old_phir[nj].clone();
                et_l[ok] = old_etl[nj].clone();
                phi_l[ok] = old_phil[nj].clone();
            }
        }
    }

    if mode_idx >= neffs.len() {
        return Err(Error::invalid(format!(
            "requested mode {mode_idx} out of range (converged {} modes)",
            neffs.len()
        )));
    }
    let sign_re = if target.re < 0.0 { -1.0 } else { 1.0 };
    let mut neff = neffs[mode_idx];
    let mut beta = neff * k0;
    let mut lambda = -(beta * beta);

    // For complex-symmetric (`!asm.is_non_symmetric && asm.a1.is_none()`) pencils, `modes_complex`
    // skips `left_eigenvector` (`y = conj(x)`), leaving Krylov eigenvector residual `~1e-6` (`~1e-12`
    // in Rayleigh quotient). Two inverse-iteration steps polish `(x, y)` to machine precision (`~1e-15`).
    if asm.a1.is_none() && !asm.is_non_symmetric {
        if let Ok(y_pol) = left_eigenvector(&asm.a, Some(&asm.b), lambda, &y_vecs[mode_idx], 2) {
            let x_pol: Vec<C> = y_pol.iter().map(|z| z.conj()).collect();
            let mut xf_full = vec![c(0.0, 0.0); asm.n_tot];
            let mut yf_full = vec![c(0.0, 0.0); asm.n_tot];
            for (a, &dof) in asm.free.iter().enumerate() {
                xf_full[dof] = x_pol[a];
                yf_full[dof] = y_pol[a];
            }
            x_vecs[mode_idx] = x_pol;
            y_vecs[mode_idx] = y_pol;
            et_r[mode_idx] = xf_full[..asm.net].to_vec();
            phi_r[mode_idx] = xf_full[asm.net..].iter().map(|&v| v / k0).collect();
            et_l[mode_idx] = yf_full[..asm.net].to_vec();
            phi_l[mode_idx] = yf_full[asm.net..].iter().map(|&v| v / k0).collect();
        }
    }

    let x_free = &x_vecs[mode_idx];
    let y_free = &y_vecs[mode_idx];
    let bx_f = asm.b.matvec(x_free);
    let denom_ybx: C = y_free.iter().zip(&bx_f).map(|(&yi, &bxi)| yi.conj() * bxi).sum();

    // Precompute quadratic NEP `(A + beta * (j A1) + beta^2 * B) x = 0` when `a1` coupling is active
    let nep_coupled = if let Some(ref a1) = asm.a1 {
        let ja1 = a1.scaled(c(0.0, 1.0));
        let nep = NonlinearEigenProblem::polynomial(vec![asm.a.clone(), ja1, asm.b.clone()])
            .map_err(|e| Error::Solver(e.to_string()))?;
        let t_beta = nep.assemble(beta);
        let y_nep =
            left_eigenvector(&t_beta, None, c(0.0, 0.0), x_free, 4).map_err(|e| Error::Solver(e.to_string()))?;
        let mut r_vec = vec![c(0.0, 0.0); x_free.len()];
        nep.apply(beta, x_free, &mut r_vec);
        let a1x = a1.matvec(x_free);
        let num_r: C = y_nep.iter().zip(&r_vec).map(|(&yi, &ri)| yi.conj() * ri).sum();
        let den_tp: C =
            y_nep.iter().enumerate().map(|(i, &yi)| yi.conj() * (c(0.0, 1.0) * a1x[i] + beta * 2.0 * bx_f[i])).sum();
        if den_tp.norm() > 1e-30 && num_r.re.is_finite() && num_r.im.is_finite() {
            beta -= num_r / den_tp;
            neff = beta / k0;
            lambda = -(beta * beta);
            neffs[mode_idx] = neff;
        }
        Some((nep, y_nep))
    } else {
        if denom_ybx.norm() > 1e-30 {
            let ax_f = asm.a.matvec(x_free);
            let num_yax: C = y_free.iter().zip(&ax_f).map(|(&yi, &axi)| yi.conj() * axi).sum();
            let lam_rq = num_yax / denom_ybx;
            if lam_rq.re.is_finite() && lam_rq.im.is_finite() {
                lambda = lam_rq;
                let mut nf = (-lambda).sqrt() / k0;
                if nf.re * sign_re < 0.0 {
                    nf = -nf;
                }
                neff = nf;
                beta = neff * k0;
                neffs[mode_idx] = neff;
            }
        }
        None
    };

    let (objective, dneff_weight, mut aux_scalars) = eval_neff_or_q_objective(&req.objective, neff)?;

    // lambda = -(k0 * neff)^2 => dneff/dp = -dlambda/dp / (2 * k0^2 * neff)
    let dlam_to_dneff = c(-1.0, 0.0) / (neff * (2.0 * k0 * k0));

    let eval_dneff_from_asm = |asm_p: &ModeAssemblyC, asm_m: &ModeAssemblyC, inv_2h: C| -> Result<C, Error> {
        let da = CCsr::axpby(inv_2h, &asm_p.a, -inv_2h, &asm_m.a);
        let db = CCsr::axpby(inv_2h, &asm_p.b, -inv_2h, &asm_m.b);
        if let Some((ref nep, ref y_nep)) = nep_coupled {
            let z = CCsr::zeros(asm.a.nrows(), asm.a.ncols());
            let a1p = asm_p.a1.as_ref().unwrap_or(&z);
            let a1m = asm_m.a1.as_ref().unwrap_or(&z);
            let dja1 = CCsr::axpby(inv_2h * c(0.0, 1.0), a1p, -inv_2h * c(0.0, 1.0), a1m);
            let dt = NonlinearEigenProblem::polynomial(vec![da, dja1, db]).map_err(|e| Error::Solver(e.to_string()))?;
            let dbeta = nonlinear_eigenvalue_derivative(nep, &dt, beta, x_free, y_nep);
            Ok(dbeta / k0)
        } else {
            let dlam = eigenvalue_derivative(lambda, x_free, y_free, Some(&asm.b), &da, Some(&db));
            Ok(dlam * dlam_to_dneff)
        }
    };

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..2).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());
    let mut projected_densities = BTreeMap::new();
    let mut volume_fractions = BTreeMap::new();
    let mut volume_fraction_gradients = BTreeMap::new();

    for (idx, param) in req.params.iter().enumerate() {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = base_problem
                    .materials
                    .get(doms[0], property)
                    .unwrap_or(if property == "eps_r" || property == "mu_r" { 1.0 } else { 0.0 });
                let h = 1e-6 * p0.abs().max(1.0);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = assemble_2d_waveguide_for_sens(&simplex, &prob_p, freq, &eps_overrides)
                    .map_err(Error::from_kernel)?;
                let asm_m = assemble_2d_waveguide_for_sens(&simplex, &prob_m, freq, &eps_overrides)
                    .map_err(Error::from_kernel)?;
                let dneff = eval_dneff_from_asm(&asm_p, &asm_m, c(1.0 / (2.0 * h), 0.0))?;
                let g = dneff_weight.re * dneff.re + dneff_weight.im * dneff.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::ParamVar { name, .. } => {
                let Some(vel) = req.shape_velocities.get(name) else {
                    return Err(Error::invalid(format!("missing precomputed shape velocity for parameter '{name}'")));
                };
                let h = vel.optimal_step(char_len, 1e-6);
                let mp = displace_simplex_mesh(&simplex, h, vel).map_err(Error::from_kernel)?;
                let mm = displace_simplex_mesh(&simplex, -h, vel).map_err(Error::from_kernel)?;
                let asm_p = assemble_2d_waveguide_for_sens(&mp, base_problem, freq, &eps_overrides)
                    .map_err(Error::from_kernel)?;
                let asm_m = assemble_2d_waveguide_for_sens(&mm, base_problem, freq, &eps_overrides)
                    .map_err(Error::from_kernel)?;
                let dneff = eval_dneff_from_asm(&asm_p, &asm_m, c(1.0 / (2.0 * h), 0.0))?;
                let g = dneff_weight.re * dneff.re + dneff_weight.im * dneff.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity {
                simp_penalty,
                ramp_q,
                val_min,
                val_max,
                val_min_im,
                val_max_im,
                plasmon_damping,
                ..
            } => {
                let st = &topo_states[&idx];
                let dcells = st.filter.design_cells();
                let edots = cell_2d_waveguide_eps_dots(
                    &simplex,
                    base_problem.order,
                    k0,
                    lambda,
                    &et_l[mode_idx],
                    &phi_l[mode_idx],
                    &et_r[mode_idx],
                    &phi_r[mode_idx],
                    dcells,
                );
                let dj_dbar: Vec<f64> = st
                    .rho_bar
                    .iter()
                    .zip(&edots)
                    .map(|(&rb, &edot)| {
                        let (_, deps) = optical_interp_of_param(
                            rb,
                            *simp_penalty,
                            *ramp_q,
                            *val_min,
                            *val_max,
                            *val_min_im,
                            *val_max_im,
                            *plasmon_damping,
                        );
                        let dlam = (edot * deps) / denom_ybx;
                        let dneff = dlam * dlam_to_dneff;
                        dneff_weight.re * dneff.re + dneff_weight.im * dneff.im
                    })
                    .collect();
                let g_rho = st.pullback(&dj_dbar)?;
                let (vf, vf_grad) = st.volume_fraction_and_grad()?;
                aux_scalars.insert(format!("volume_fraction_{idx}"), vf);
                projected_densities.insert(idx, st.rho_bar.clone());
                volume_fractions.insert(idx, vf);
                volume_fraction_gradients.insert(idx, vf_grad);
                flat_gradient.extend_from_slice(&g_rho);
                gradients.push(g_rho);
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to study 'eme'"));
            }
        }
    }

    let evals: Vec<C> = neffs.iter().map(|&n| -((n * k0) * (n * k0))).collect();
    let rep = SensitivityReport {
        objective,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities,
        volume_fractions,
        volume_fraction_gradients,
        mode_basis: Some(TrackedModeBasis { eigenvalues: evals, eigenvectors: x_vecs, mass_matrix: Some(asm.b) }),
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}

// ---------------------------------------------------------------------------
// Shared 3D Nedelec elementwise permittivity contraction helpers
// ---------------------------------------------------------------------------

fn build_3d_delta_eps_mass_full(
    mesh: &Discretisation,
    ned: &Ned3,
    problem: &Problem,
    eps_overrides: &[(usize, C)],
) -> Result<CCsr, String> {
    let (base_eps, _) = simplex_cell_materials_3d(mesh, problem)?;
    let zero_c3: C3 = [[c(0.0, 0.0); 3]; 3];
    let mut delta_eps = vec![zero_c3; mesh.n_cells()];
    let mut zero_mui = vec![zero_c3; mesh.n_cells()];
    for &(ci, ec) in eps_overrides {
        let mut d = zero_c3;
        for a in 0..3 {
            for b in 0..3 {
                let target = if a == b { ec } else { c(0.0, 0.0) };
                d[a][b] = target - base_eps[ci][a][b];
            }
        }
        delta_eps[ci] = d;
    }
    let _ = &mut zero_mui;
    let (_, dm_full) = assemble_maxwell_with_pml(mesh, ned, &delta_eps, &zero_mui, None);
    Ok(dm_full)
}

/// Evaluates `integral_{Omega_c} (E_a(x) . E_b(x)) dV` (or `conj(E_a) . E_b` when `conj_first == true`)
/// for full-space Nedelec vectors `vec_a`, `vec_b` on each cell in `cells` (supporting straight and `CurvedTet` cells).
fn cell_3d_nedelec_mass_dots(
    mesh: &Discretisation,
    ned: &Ned3,
    vec_a: &[C],
    vec_b: &[C],
    conj_first: bool,
    cells: &[usize],
) -> Vec<C> {
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    let cquad = quad_rule_collapsed(6);
    let nb = ned.nb();
    const CHUNK: usize = 256;
    let parts = femcore::par::map_chunks(cells, CHUNK, |slice| {
        let mut out = Vec::with_capacity(slice.len());
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];
        for &ci in slice {
            let gid = ned.key_pairs(mesh, ci);
            for b in 0..4 {
                for k in 0..3 {
                    gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                }
            }
            let cdofs = &ned.cell_dofs[ci];
            let mut acc = c(0.0, 0.0);
            if mesh.cell_is_curved(ci) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (_, det) = ct.jacobian(&lam);
                    let gq = safe_curved_jacobian(&ct, &lam, &gs);
                    basis3_keyed(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    let mut ea = [c(0.0, 0.0); 3];
                    let mut eb = [c(0.0, 0.0); 3];
                    for s in 0..nb {
                        let ca = vec_a[cdofs[s]];
                        let cb = vec_b[cdofs[s]];
                        for k in 0..3 {
                            ea[k] += ca * w[s][k];
                            eb[k] += cb * w[s][k];
                        }
                    }
                    let dot = if conj_first {
                        ea[0].conj() * eb[0] + ea[1].conj() * eb[1] + ea[2].conj() * eb[2]
                    } else {
                        ea[0] * eb[0] + ea[1] * eb[1] + ea[2] * eb[2]
                    };
                    acc += dot * wq;
                }
            } else {
                let vol = geo.vol[ci];
                for q in 0..quad.w.len() {
                    let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                    basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                    let wq = quad.w[q] * vol;
                    let mut ea = [c(0.0, 0.0); 3];
                    let mut eb = [c(0.0, 0.0); 3];
                    for s in 0..nb {
                        let ca = vec_a[cdofs[s]];
                        let cb = vec_b[cdofs[s]];
                        for k in 0..3 {
                            ea[k] += ca * w[s][k];
                            eb[k] += cb * w[s][k];
                        }
                    }
                    let dot = if conj_first {
                        ea[0].conj() * eb[0] + ea[1].conj() * eb[1] + ea[2].conj() * eb[2]
                    } else {
                        ea[0] * eb[0] + ea[1] * eb[1] + ea[2] * eb[2]
                    };
                    acc += dot * wq;
                }
            }
            out.push(acc);
        }
        out
    });
    parts.into_iter().flatten().collect()
}

// ---------------------------------------------------------------------------
// 2. 3D Cavity Eigenmodes (`cavity`)
// ---------------------------------------------------------------------------

fn assemble_cavity_for_sens(
    mesh: &Discretisation,
    problem: &Problem,
    eps_overrides: &[(usize, C)],
) -> Result<crate::cavity::EmwAssembly, String> {
    let mut asm = assemble_emw_eigen(mesh, problem)?;
    if !eps_overrides.is_empty() {
        let dm_full = build_3d_delta_eps_mass_full(mesh, &asm.ned, problem, eps_overrides)?;
        let dm_red = femcore::periodic::congruence(&asm.red.prolong(), &dm_full);
        asm.mc = asm.mc.add_scaled(&dm_red, c(1.0, 0.0));
    }
    Ok(asm)
}

fn cavity_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let simplex = discretise(mesh)?;
    let topo_states = build_topology_states(&simplex, ctx, req)?;
    let mut eps_overrides: Vec<(usize, C)> = Vec::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity {
            simp_penalty,
            ramp_q,
            val_min,
            val_max,
            val_min_im,
            val_max_im,
            plasmon_damping,
            ..
        } = param
        {
            let st = &topo_states[&idx];
            for (&ci, &rb) in st.filter.design_cells().iter().zip(&st.rho_bar) {
                let (ec, _) = optical_interp_of_param(
                    rb,
                    *simp_penalty,
                    *ramp_q,
                    *val_min,
                    *val_max,
                    *val_min_im,
                    *val_max_im,
                    *plasmon_damping,
                );
                eps_overrides.push((ci, ec));
            }
        }
    }

    let mode_idx = target_mode_index(&req.objective);
    let mut work_prob = base_problem.clone();
    if let Study::Cavity { ref mut n_modes, .. } = work_prob.study {
        *n_modes = (*n_modes).max(mode_idx + 1);
    }

    let asm = assemble_cavity_for_sens(&simplex, &work_prob, &eps_overrides).map_err(Error::from_kernel)?;
    let solved = solve_emw_modes(&asm).map_err(Error::from_kernel)?;
    if solved.freqs.is_empty() {
        return Err(Error::Solver("no converged 3D cavity modes found".into()));
    }

    let mut freqs: Vec<C> = solved.freqs.iter().map(|&(f, _)| f).collect();
    let mut x_vecs: Vec<Vec<C>> = solved.freqs.iter().map(|&(_, vi)| solved.vectors[vi].clone()).collect();

    if let Some(ref ref_basis) = req.reference_modes {
        if ref_basis.eigenvectors.len() == x_vecs.len()
            && ref_basis.eigenvectors.first().is_some_and(|v| v.len() == asm.nred)
        {
            let matched = match_modes(&ref_basis.eigenvectors, &x_vecs, Some(&asm.mc), 0.1, 0.05);
            let (old_f, old_x) = (freqs.clone(), x_vecs.clone());
            for (ok, &nj) in matched.assignment.iter().enumerate() {
                freqs[ok] = old_f[nj];
                x_vecs[ok] = old_x[nj].clone();
            }
        }
    }

    if mode_idx >= freqs.len() {
        return Err(Error::invalid(format!(
            "requested cavity mode {mode_idx} out of range (converged {} modes)",
            freqs.len()
        )));
    }
    let f_hz = freqs[mode_idx];
    let k0_eig = f_hz * (2.0 * PI / C0);
    let lambda = k0_eig * k0_eig;
    let x_red = &x_vecs[mode_idx];
    let y_red = if !asm.is_non_symmetric && asm.kc.is_hermitian(1e-10) && asm.mc.is_hermitian(1e-10) {
        x_red.clone()
    } else if !asm.is_non_symmetric && asm.kc.is_symmetric(1e-10) && asm.mc.is_symmetric(1e-10) {
        x_red.iter().map(|z| z.conj()).collect()
    } else {
        left_eigenvector(&asm.kc, Some(&asm.mc), lambda, x_red, 4).map_err(|e| Error::Solver(e.to_string()))?
    };

    let (objective, df_weight, mut aux_scalars) = eval_freq_or_q_objective(&req.objective, f_hz)?;
    // lambda = (2 pi f / c0)^2 => df/dp = (c0^2 / (8 pi^2 f)) * dlambda/dp
    let dlam_to_df = c(C0 * C0 / (8.0 * PI * PI), 0.0) / f_hz;

    let mx_red = asm.mc.matvec(x_red);
    let denom_ymx: C = y_red.iter().zip(&mx_red).map(|(&yi, &mxi)| yi.conj() * mxi).sum();
    let x_full = asm.red.expand(x_red);
    let y_full = asm.red.expand(&y_red); // For `y_red^H (P^H dM P) x_red == (P y_red)^H dM (P x_red)`!

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..3).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());
    let mut projected_densities = BTreeMap::new();
    let mut volume_fractions = BTreeMap::new();
    let mut volume_fraction_gradients = BTreeMap::new();

    for (idx, param) in req.params.iter().enumerate() {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = work_prob
                    .materials
                    .get(doms[0], property)
                    .unwrap_or(if property == "eps_r" || property == "mu_r" { 1.0 } else { 0.0 });
                let h = 1e-6 * p0.abs().max(1.0);
                let mut prob_p = work_prob.clone();
                let mut prob_m = work_prob.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = assemble_cavity_for_sens(&simplex, &prob_p, &eps_overrides).map_err(Error::from_kernel)?;
                let asm_m = assemble_cavity_for_sens(&simplex, &prob_m, &eps_overrides).map_err(Error::from_kernel)?;
                let inv_2h = c(1.0 / (2.0 * h), 0.0);
                let dk = CCsr::axpby(inv_2h, &asm_p.kc, -inv_2h, &asm_m.kc);
                let dm = CCsr::axpby(inv_2h, &asm_p.mc, -inv_2h, &asm_m.mc);
                let dlam = eigenvalue_derivative(lambda, x_red, &y_red, Some(&asm.mc), &dk, Some(&dm));
                let df = dlam * dlam_to_df;
                let g = df_weight.re * df.re + df_weight.im * df.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::ParamVar { name, .. } => {
                let Some(vel) = req.shape_velocities.get(name) else {
                    return Err(Error::invalid(format!("missing precomputed shape velocity for parameter '{name}'")));
                };
                let h = vel.optimal_step(char_len, 1e-6);
                let mp = displace_simplex_mesh(&simplex, h, vel).map_err(Error::from_kernel)?;
                let mm = displace_simplex_mesh(&simplex, -h, vel).map_err(Error::from_kernel)?;
                let asm_p = assemble_cavity_for_sens(&mp, &work_prob, &eps_overrides).map_err(Error::from_kernel)?;
                let asm_m = assemble_cavity_for_sens(&mm, &work_prob, &eps_overrides).map_err(Error::from_kernel)?;
                let inv_2h = c(1.0 / (2.0 * h), 0.0);
                let dk = CCsr::axpby(inv_2h, &asm_p.kc, -inv_2h, &asm_m.kc);
                let dm = CCsr::axpby(inv_2h, &asm_p.mc, -inv_2h, &asm_m.mc);
                let dlam = eigenvalue_derivative(lambda, x_red, &y_red, Some(&asm.mc), &dk, Some(&dm));
                let df = dlam * dlam_to_df;
                let g = df_weight.re * df.re + df_weight.im * df.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity {
                simp_penalty,
                ramp_q,
                val_min,
                val_max,
                val_min_im,
                val_max_im,
                plasmon_damping,
                ..
            } => {
                let st = &topo_states[&idx];
                let dcells = st.filter.design_cells();
                let mdots = cell_3d_nedelec_mass_dots(&simplex, &asm.ned, &y_full, &x_full, true, dcells);
                let dj_dbar: Vec<f64> = st
                    .rho_bar
                    .iter()
                    .zip(&mdots)
                    .map(|(&rb, &mdot)| {
                        let (_, deps) = optical_interp_of_param(
                            rb,
                            *simp_penalty,
                            *ramp_q,
                            *val_min,
                            *val_max,
                            *val_min_im,
                            *val_max_im,
                            *plasmon_damping,
                        );
                        // dK/deps = 0, dM/deps = mdot => dlambda = -lambda * deps * mdot / denom_ymx
                        let dlam = -(lambda * deps * mdot) / denom_ymx;
                        let df = dlam * dlam_to_df;
                        df_weight.re * df.re + df_weight.im * df.im
                    })
                    .collect();
                let g_rho = st.pullback(&dj_dbar)?;
                let (vf, vf_grad) = st.volume_fraction_and_grad()?;
                aux_scalars.insert(format!("volume_fraction_{idx}"), vf);
                projected_densities.insert(idx, st.rho_bar.clone());
                volume_fractions.insert(idx, vf);
                volume_fraction_gradients.insert(idx, vf_grad);
                flat_gradient.extend_from_slice(&g_rho);
                gradients.push(g_rho);
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to study 'eme'"));
            }
        }
    }

    let evals: Vec<C> = freqs
        .iter()
        .map(|&f| {
            let k = f * (2.0 * PI / C0);
            k * k
        })
        .collect();
    let rep = SensitivityReport {
        objective,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities,
        volume_fractions,
        volume_fraction_gradients,
        mode_basis: Some(TrackedModeBasis { eigenvalues: evals, eigenvectors: x_vecs, mass_matrix: Some(asm.mc) }),
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}

// ---------------------------------------------------------------------------
// 3. 3D Floquet-Bloch Modes (`bloch`)
// ---------------------------------------------------------------------------

fn assemble_bloch_for_sens(
    mesh: &Discretisation,
    problem: &Problem,
    eps_overrides: &[(usize, C)],
) -> Result<crate::bloch::BlochAssembly, String> {
    let mut asm = assemble_bloch(mesh, problem)?;
    if !eps_overrides.is_empty() {
        let dm_full = build_3d_delta_eps_mass_full(mesh, &asm.ned, problem, eps_overrides)?;
        let dm_red = femcore::periodic::congruence(&asm.red.prolong(), &dm_full);
        asm.a0_red = asm.a0_red.add_scaled(&dm_red, c(-1.0, 0.0));
    }
    Ok(asm)
}

fn bloch_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let simplex = discretise(mesh)?;
    let topo_states = build_topology_states(&simplex, ctx, req)?;
    let mut eps_overrides: Vec<(usize, C)> = Vec::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity {
            simp_penalty,
            ramp_q,
            val_min,
            val_max,
            val_min_im,
            val_max_im,
            plasmon_damping,
            ..
        } = param
        {
            let st = &topo_states[&idx];
            for (&ci, &rb) in st.filter.design_cells().iter().zip(&st.rho_bar) {
                let (ec, _) = optical_interp_of_param(
                    rb,
                    *simp_penalty,
                    *ramp_q,
                    *val_min,
                    *val_max,
                    *val_min_im,
                    *val_max_im,
                    *plasmon_damping,
                );
                eps_overrides.push((ci, ec));
            }
        }
    }

    let mode_idx = target_mode_index(&req.objective);
    let mut work_prob = base_problem.clone();
    if let Study::Bloch { ref mut n_modes, .. } = work_prob.study {
        *n_modes = (*n_modes).max(mode_idx + 1);
    }

    let asm = assemble_bloch_for_sens(&simplex, &work_prob, &eps_overrides).map_err(Error::from_kernel)?;
    let solved = solve_bloch_eigen(&asm).map_err(Error::from_kernel)?;
    if solved.modes.is_empty() {
        return Err(Error::Solver("no converged 3D Bloch modes found".into()));
    }

    let mut nus: Vec<C> = solved.modes.iter().map(|&(nb, _)| nb).collect();
    let mut u_vecs: Vec<Vec<C>> = solved.modes.iter().map(|&(_, vi)| solved.vectors[vi].clone()).collect();

    if let Some(ref ref_basis) = req.reference_modes {
        if ref_basis.eigenvectors.len() == u_vecs.len()
            && ref_basis.eigenvectors.first().is_some_and(|v| v.len() == asm.nred)
        {
            let matched = match_modes(&ref_basis.eigenvectors, &u_vecs, Some(&asm.a2_red), 0.1, 0.05);
            let (old_n, old_u) = (nus.clone(), u_vecs.clone());
            for (ok, &nj) in matched.assignment.iter().enumerate() {
                nus[ok] = old_n[nj];
                u_vecs[ok] = old_u[nj].clone();
            }
        }
    }

    if mode_idx >= nus.len() {
        return Err(Error::invalid(format!(
            "requested Bloch mode {mode_idx} out of range (converged {} modes)",
            nus.len()
        )));
    }
    let nu = nus[mode_idx];
    let u_red = &u_vecs[mode_idx];
    let nep = NonlinearEigenProblem::polynomial(vec![asm.a0_red.clone(), asm.a1_red.clone(), asm.a2_red.clone()])
        .map_err(|e| Error::Solver(e.to_string()))?;
    let t_nu = nep.assemble(nu);
    let y_red = left_eigenvector(&t_nu, None, c(0.0, 0.0), u_red, 4).map_err(|e| Error::Solver(e.to_string()))?;

    // Denominator `y^H (A1 + 2 nu A2) u` for topology fast path
    let a1u = asm.a1_red.matvec(u_red);
    let a2u = asm.a2_red.matvec(u_red);
    let denom_tprime: C = y_red.iter().enumerate().map(|(i, &yi)| yi.conj() * (a1u[i] + nu * 2.0 * a2u[i])).sum();
    let u_full = asm.red.expand(u_red);
    let y_full = asm.red.expand(&y_red);

    let (objective, dnu_weight, mut aux_scalars) = eval_neff_or_q_objective(&req.objective, nu)?;

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..3).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());
    let mut projected_densities = BTreeMap::new();
    let mut volume_fractions = BTreeMap::new();
    let mut volume_fraction_gradients = BTreeMap::new();

    for (idx, param) in req.params.iter().enumerate() {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = work_prob
                    .materials
                    .get(doms[0], property)
                    .unwrap_or(if property == "eps_r" || property == "mu_r" { 1.0 } else { 0.0 });
                let h = 1e-6 * p0.abs().max(1.0);
                let mut prob_p = work_prob.clone();
                let mut prob_m = work_prob.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = assemble_bloch_for_sens(&simplex, &prob_p, &eps_overrides).map_err(Error::from_kernel)?;
                let asm_m = assemble_bloch_for_sens(&simplex, &prob_m, &eps_overrides).map_err(Error::from_kernel)?;
                let inv_2h = c(1.0 / (2.0 * h), 0.0);
                let da0 = CCsr::axpby(inv_2h, &asm_p.a0_red, -inv_2h, &asm_m.a0_red);
                let da1 = CCsr::axpby(inv_2h, &asm_p.a1_red, -inv_2h, &asm_m.a1_red);
                let da2 = CCsr::axpby(inv_2h, &asm_p.a2_red, -inv_2h, &asm_m.a2_red);
                let dt =
                    NonlinearEigenProblem::polynomial(vec![da0, da1, da2]).map_err(|e| Error::Solver(e.to_string()))?;
                let dnu = nonlinear_eigenvalue_derivative(&nep, &dt, nu, u_red, &y_red);
                let g = dnu_weight.re * dnu.re + dnu_weight.im * dnu.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::ParamVar { name, .. } => {
                let Some(vel) = req.shape_velocities.get(name) else {
                    return Err(Error::invalid(format!("missing precomputed shape velocity for parameter '{name}'")));
                };
                let h = vel.optimal_step(char_len, 1e-6);
                let mp = displace_simplex_mesh(&simplex, h, vel).map_err(Error::from_kernel)?;
                let mm = displace_simplex_mesh(&simplex, -h, vel).map_err(Error::from_kernel)?;
                let asm_p = assemble_bloch_for_sens(&mp, &work_prob, &eps_overrides).map_err(Error::from_kernel)?;
                let asm_m = assemble_bloch_for_sens(&mm, &work_prob, &eps_overrides).map_err(Error::from_kernel)?;
                let inv_2h = c(1.0 / (2.0 * h), 0.0);
                let da0 = CCsr::axpby(inv_2h, &asm_p.a0_red, -inv_2h, &asm_m.a0_red);
                let da1 = CCsr::axpby(inv_2h, &asm_p.a1_red, -inv_2h, &asm_m.a1_red);
                let da2 = CCsr::axpby(inv_2h, &asm_p.a2_red, -inv_2h, &asm_m.a2_red);
                let dt =
                    NonlinearEigenProblem::polynomial(vec![da0, da1, da2]).map_err(|e| Error::Solver(e.to_string()))?;
                let dnu = nonlinear_eigenvalue_derivative(&nep, &dt, nu, u_red, &y_red);
                let g = dnu_weight.re * dnu.re + dnu_weight.im * dnu.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity {
                simp_penalty,
                ramp_q,
                val_min,
                val_max,
                val_min_im,
                val_max_im,
                plasmon_damping,
                ..
            } => {
                let st = &topo_states[&idx];
                let dcells = st.filter.design_cells();
                let mdots = cell_3d_nedelec_mass_dots(&simplex, &asm.ned, &y_full, &u_full, true, dcells);
                let dj_dbar: Vec<f64> = st
                    .rho_bar
                    .iter()
                    .zip(&mdots)
                    .map(|(&rb, &mdot)| {
                        let (_, deps) = optical_interp_of_param(
                            rb,
                            *simp_penalty,
                            *ramp_q,
                            *val_min,
                            *val_max,
                            *val_min_im,
                            *val_max_im,
                            *plasmon_damping,
                        );
                        // dT/d(rho_bar) = dA0/d(rho_bar) = -deps * M_c => dnu = -y^H (-deps M_c) u / y^H T'(nu) u = +deps * mdot / denom_tprime
                        let dnu = (deps * mdot) / denom_tprime;
                        dnu_weight.re * dnu.re + dnu_weight.im * dnu.im
                    })
                    .collect();
                let g_rho = st.pullback(&dj_dbar)?;
                let (vf, vf_grad) = st.volume_fraction_and_grad()?;
                aux_scalars.insert(format!("volume_fraction_{idx}"), vf);
                projected_densities.insert(idx, st.rho_bar.clone());
                volume_fractions.insert(idx, vf);
                volume_fraction_gradients.insert(idx, vf_grad);
                flat_gradient.extend_from_slice(&g_rho);
                gradients.push(g_rho);
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to study 'eme'"));
            }
        }
    }

    let rep = SensitivityReport {
        objective,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities,
        volume_fractions,
        volume_fraction_gradients,
        mode_basis: Some(TrackedModeBasis { eigenvalues: nus, eigenvectors: u_vecs, mass_matrix: Some(asm.a2_red) }),
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}

// ---------------------------------------------------------------------------
// 4. 3D Driven Frequency-Domain Scattering & Field Adjoints (`driven`)
// ---------------------------------------------------------------------------

fn assemble_driven_for_sens(
    mesh: &Discretisation,
    problem: &Problem,
    eps_overrides: &[(usize, C)],
) -> Result<DrivenAssembly, String> {
    let mut asm = assemble_emw_driven(mesh, problem)?;
    if !eps_overrides.is_empty() {
        let dm_full = build_3d_delta_eps_mass_full(mesh, &asm.ned, problem, eps_overrides)?;
        let dm_red = femcore::periodic::congruence(&asm.red.prolong(), &dm_full);
        let dm_hat = scale_symmetric(&dm_red, &asm.d_ruiz);
        asm.terms_hat[2] = asm.terms_hat[2].add_scaled(&dm_hat, c(1.0, 0.0));
    }
    Ok(asm)
}

fn driven_freq_grid(asm: &DrivenAssembly) -> Vec<f64> {
    if let Some(ref fs) = asm.freqs_hz {
        if !fs.is_empty() {
            return fs.clone();
        }
    }
    if asm.n_points <= 1 {
        vec![asm.f_min]
    } else {
        (0..asm.n_points)
            .map(|k| asm.f_min + (asm.f_max - asm.f_min) * (k as f64) / ((asm.n_points - 1) as f64))
            .collect()
    }
}

fn unscale_sym_matrix(mat: &CCsr, d_ruiz: &[f64]) -> CCsr {
    let mut out = mat.clone();
    for i in 0..out.nrows() {
        let inv_di = 1.0 / d_ruiz[i];
        let rg = out.indptr()[i]..out.indptr()[i + 1];
        let cols: Vec<u32> = out.indices()[rg.clone()].to_vec();
        for (k, &col) in cols.iter().enumerate() {
            out.values_mut()[rg.start + k] *= inv_di / d_ruiz[col as usize];
        }
    }
    out
}

fn expand_conj_prolong(red: &femcore::periodic::Reduction, v_red: &[C]) -> Vec<C> {
    let v_conj: Vec<C> = v_red.iter().map(|z| z.conj()).collect();
    red.expand(&v_conj).into_iter().map(|z| z.conj()).collect()
}

struct DrivenFreqAdjointPair {
    freq_idx: usize,
    omega: f64,
    in_port: usize,
    /// Unscaled reduced forward solution `x_red` (`nred`) for `in_port`.
    x_red: Vec<C>,
    /// Unscaled reduced adjoint vector `lam_red` (`nred`) solving `A_red^T lam_red = g_red`.
    lam_red: Vec<C>,
    /// Per active `out_port`, `(mu_pq, z_red_p)` where `A_red^T z_red_p = z_wave_p^(-1/2) t_red_p`
    /// so explicit port-functional derivatives under shape perturbations can be added cleanly.
    port_weights: Vec<(usize, C)>,
}

struct PromSensSolver {
    inner: SparseDrivenSolver,
}

impl PromLinearSolver for PromSensSolver {
    fn update_and_solve(&mut self, _omega: f64, a_mat: &CCsr, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        self.inner.update_matrix(a_mat)?;
        self.inner.solve_many(rhs_cols)
    }

    fn solve_same_factor(&self, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        self.inner.solve_many(rhs_cols)
    }
}

fn driven_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let simplex = discretise(mesh)?;
    let topo_states = build_topology_states(&simplex, ctx, req)?;
    let mut eps_overrides: Vec<(usize, C)> = Vec::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity {
            simp_penalty,
            ramp_q,
            val_min,
            val_max,
            val_min_im,
            val_max_im,
            plasmon_damping,
            ..
        } = param
        {
            let st = &topo_states[&idx];
            for (&ci, &rb) in st.filter.design_cells().iter().zip(&st.rho_bar) {
                let (ec, _) = optical_interp_of_param(
                    rb,
                    *simp_penalty,
                    *ramp_q,
                    *val_min,
                    *val_max,
                    *val_min_im,
                    *val_max_im,
                    *plasmon_damping,
                );
                eps_overrides.push((ci, ec));
            }
        }
    }

    let asm = assemble_driven_for_sens(&simplex, base_problem, &eps_overrides).map_err(Error::from_kernel)?;
    let p_ports = asm.ports.len();
    if p_ports == 0 {
        return Err(Error::invalid("driven sensitivity requires at least one port"));
    }
    let all_freqs = driven_freq_grid(&asm);

    let (active_freq_indices, in_port, needed_out_ports, field_power_domain): (
        Vec<usize>,
        usize,
        Vec<usize>,
        Option<Vec<usize>>,
    ) = match &req.objective {
        SensitivityObjective::Transmission { out_port, in_port, freq_index, .. } => {
            if *out_port >= p_ports || *in_port >= p_ports {
                return Err(Error::invalid(format!(
                    "port index ({out_port}, {in_port}) out of range for {p_ports} ports"
                )));
            }
            let f_idxs = match freq_index {
                Some(fi) => {
                    if *fi >= all_freqs.len() {
                        return Err(Error::invalid(format!(
                            "freq_index {fi} out of range for {} frequencies",
                            all_freqs.len()
                        )));
                    }
                    vec![*fi]
                }
                None => (0..all_freqs.len()).collect(),
            };
            (f_idxs, *in_port, vec![*out_port], None)
        }
        SensitivityObjective::SplitterBalance { in_port, out_ports, freq_index, .. } => {
            if *in_port >= p_ports || out_ports.iter().any(|&p| p >= p_ports) {
                return Err(Error::invalid("splitter_balance port index out of range"));
            }
            let f_idxs = match freq_index {
                Some(fi) => {
                    if *fi >= all_freqs.len() {
                        return Err(Error::invalid("splitter_balance freq_index out of range"));
                    }
                    vec![*fi]
                }
                None => (0..all_freqs.len()).collect(),
            };
            let mut uniq = out_ports.clone();
            uniq.sort_unstable();
            uniq.dedup();
            (f_idxs, *in_port, uniq, None)
        }
        SensitivityObjective::FieldPower { domain, in_port, freq_index } => {
            if *in_port >= p_ports {
                return Err(Error::invalid("field_power in_port out of range"));
            }
            let doms = ctx.resolver.domains(domain)?;
            if doms.is_empty() {
                return Err(Error::invalid(format!("field_power domain '{domain}' resolved to zero domains")));
            }
            let cells: Vec<usize> =
                (0..simplex.n_cells()).filter(|&c| doms.contains(&simplex.cell_labels[c])).collect();
            let f_idxs = match freq_index {
                Some(fi) => {
                    if *fi >= all_freqs.len() {
                        return Err(Error::invalid("field_power freq_index out of range"));
                    }
                    vec![*fi]
                }
                None => (0..all_freqs.len()).collect(),
            };
            (f_idxs, *in_port, Vec::new(), Some(cells))
        }
        other => {
            return Err(Error::Unsupported(format!(
                "ems driven sensitivity supports Transmission, SplitterBalance, and FieldPower, got '{other:?}'"
            )))
        }
    };

    let active_set: BTreeSet<usize> = active_freq_indices.iter().copied().collect();

    // Per active frequency, store `(omega, x_red_in, s_row_for_in, unit_adj_z_red_by_out_port)`
    type FreqSolveRecord = (usize, f64, Vec<C>, Vec<C>, BTreeMap<usize, Vec<C>>);
    let mut freq_records: Vec<FreqSolveRecord> = Vec::with_capacity(active_freq_indices.len());

    let use_prom = asm.prom_tol.is_some() && field_power_domain.is_none();
    if use_prom {
        let prom_tol = asm.prom_tol.unwrap_or(1e-4);
        let mut solver = PromSensSolver {
            inner: SparseDrivenSolver::new_with_symmetry(&asm.solver, asm.blr_tol, !PromOperator::is_symmetric(&asm)),
        };
        run_greedy_prom_sweep(
            &asm,
            &mut solver,
            &all_freqs,
            prom_tol,
            asm.max_hdm_solves,
            asm.strict_prom_tol,
            &asm.seed_freqs_hz,
            |band, b_omegas, prom: &ImageProm<'_, DrivenAssembly>, xr_band| {
                let r = prom.dim();
                let m_terms = asm.terms_hat.len();
                let b_recv = prom.b_r_recv();
                for (k_loc, &w) in b_omegas.iter().enumerate() {
                    let g_idx = band.start + k_loc;
                    if !active_set.contains(&g_idx) {
                        continue;
                    }
                    let xr_k = &xr_band[k_loc];
                    let x_red = prom.reconstruct_unscaled(&xr_k[in_port]);
                    let mut s_col = vec![c(0.0, 0.0); p_ports];
                    let mut inv_sqrt_zw = vec![0.0; p_ports];
                    for p in 0..p_ports {
                        let (zw, _, _, _) = asm.ports[p].excitation_coeffs(w, &asm.waveports, &asm.floquet_ports);
                        inv_sqrt_zw[p] = 1.0 / zw.max(1e-30).sqrt();
                        let mut v_pq = c(0.0, 0.0);
                        for i in 0..r {
                            v_pq += xr_k[in_port][i] * b_recv[i][p];
                        }
                        let delta = if p == in_port { c(1.0, 0.0) } else { c(0.0, 0.0) };
                        s_col[p] = v_pq * inv_sqrt_zw[p] - delta;
                    }

                    // Solve reduced transpose system `A_r(w)^T z_{r,p} = inv_sqrt_zw[p] * b_{r,recv}[., p]`
                    let coeffs = asm.term_coeffs(w);
                    let mut tr_trans = vec![c(0.0, 0.0); r * r];
                    for m in 0..m_terms {
                        let alpha = coeffs[m].0;
                        if alpha.re == 0.0 && alpha.im == 0.0 {
                            continue;
                        }
                        let am = &prom.a_r[m];
                        for i in 0..r {
                            for j in 0..r {
                                tr_trans[i * r + j] += alpha * am[j * r + i];
                            }
                        }
                    }
                    let rhs_adj_r: Vec<Vec<C>> = needed_out_ports
                        .iter()
                        .map(|&p| {
                            let sc = inv_sqrt_zw[p];
                            (0..r).map(|i| b_recv[i][p] * sc).collect()
                        })
                        .collect();
                    let zr_cols = solve_dense_complex_multi(&tr_trans, &rhs_adj_r, r)
                        .ok_or_else(|| "singular reduced PROM transpose matrix in sensitivity".to_string())?;
                    let mut z_map = BTreeMap::new();
                    for (idx_p, &p) in needed_out_ports.iter().enumerate() {
                        z_map.insert(p, prom.reconstruct_unscaled(&zr_cols[idx_p]));
                    }
                    freq_records.push((g_idx, w, x_red, s_col, z_map));
                }
                Ok(())
            },
        )
        .map_err(Error::from_kernel)?;
    } else {
        let is_sym = PromOperator::is_symmetric(&asm);
        let mut fwd_solver = SparseDrivenSolver::new_with_symmetry(&asm.solver, asm.blr_tol, !is_sym);
        let mut trans_solver = (!is_sym).then(|| SparseDrivenSolver::new_with_symmetry(&asm.solver, asm.blr_tol, true));

        // Preassemble domain mass matrix if `FieldPower` is used
        let m_dom_full_opt = if let Some(ref dom_cells) = field_power_domain {
            let zero_c3: C3 = [[c(0.0, 0.0); 3]; 3];
            let one_c3: C3 = [
                [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)],
                [c(0.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)],
                [c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)],
            ];
            let mut eps_mask = vec![zero_c3; simplex.n_cells()];
            let zero_mui = vec![zero_c3; simplex.n_cells()];
            for &ci in dom_cells {
                eps_mask[ci] = one_c3;
            }
            let (_, m_dom) = assemble_maxwell_with_pml(&simplex, &asm.ned, &eps_mask, &zero_mui, None);
            Some(femcore::periodic::congruence(&asm.red.prolong(), &m_dom))
        } else {
            None
        };

        for &fi in &active_freq_indices {
            let w = 2.0 * PI * all_freqs[fi];
            let a_hat = asm.assemble_matrix_at(w);
            fwd_solver.update_matrix(&a_hat).map_err(Error::from_kernel)?;

            let (_, g_in, _, _) = asm.ports[in_port].excitation_coeffs(w, &asm.waveports, &asm.floquet_ports);
            let rhs_in: Vec<C> = asm.ports[in_port].drive_hat().into_iter().map(|b| g_in * b).collect();
            let y_in_cols = fwd_solver.solve_many(&[rhs_in]).map_err(Error::from_kernel)?;
            let y_in = &y_in_cols[0];
            let x_red: Vec<C> = y_in.iter().zip(&asm.d_ruiz).map(|(&y, &d)| y * d).collect();

            let mut s_col = vec![c(0.0, 0.0); p_ports];
            let mut inv_sqrt_zw = vec![0.0; p_ports];
            for p in 0..p_ports {
                let (zw, _, _, _) = asm.ports[p].excitation_coeffs(w, &asm.waveports, &asm.floquet_ports);
                inv_sqrt_zw[p] = 1.0 / zw.max(1e-30).sqrt();
                let v_pq = asm.ports[p].recv_dot(y_in);
                let delta = if p == in_port { c(1.0, 0.0) } else { c(0.0, 0.0) };
                s_col[p] = v_pq * inv_sqrt_zw[p] - delta;
            }

            let mut adj_rhs_cols: Vec<Vec<C>> = Vec::new();
            if let Some(ref m_dom_red) = m_dom_full_opt {
                // J = x_red^H M_dom x_red => dJ = 2 Re( (M_dom^T conj(x_red))^T dx_red )
                // Since M_dom is real symmetric (or Hermitian `M_dom^T = conj(M_dom)`):
                // `(x_red^H M_dom) dx_red = sum_j (sum_i conj(x_i) M_{ij}) dx_j = (M_dom^T conj(x_red))^T dx_red`!
                let x_conj: Vec<C> = x_red.iter().map(|z| z.conj()).collect();
                let g_red = m_dom_red.transpose().matvec(&x_conj);
                let g_hat: Vec<C> = g_red.iter().zip(&asm.d_ruiz).map(|(&g, &d)| g * d).collect();
                adj_rhs_cols.push(g_hat);
            } else {
                for &p in &needed_out_ports {
                    let sc = inv_sqrt_zw[p];
                    let t_hat: Vec<C> = match &asm.ports[p].c_hat {
                        Some((_, r)) => r.iter().map(|&ri| ri * sc).collect(),
                        None => asm.ports[p].b_hat.iter().map(|&bi| c(bi * sc, 0.0)).collect(),
                    };
                    adj_rhs_cols.push(t_hat);
                }
            }

            let z_hat_cols = if let Some(ref mut ts) = trans_solver {
                let a_hat_t = a_hat.transpose();
                ts.update_matrix(&a_hat_t).map_err(Error::from_kernel)?;
                ts.solve_many(&adj_rhs_cols).map_err(Error::from_kernel)?
            } else {
                // Symmetric operator `A_hat^T == A_hat`: reuse the exact same cached Stage-2 factorization!
                fwd_solver.solve_many(&adj_rhs_cols).map_err(Error::from_kernel)?
            };

            let mut z_map = BTreeMap::new();
            if m_dom_full_opt.is_some() {
                let z_unscaled: Vec<C> = z_hat_cols[0].iter().zip(&asm.d_ruiz).map(|(&zh, &d)| zh * d).collect();
                z_map.insert(0, z_unscaled);
            } else {
                for (idx_p, &p) in needed_out_ports.iter().enumerate() {
                    let z_unscaled: Vec<C> =
                        z_hat_cols[idx_p].iter().zip(&asm.d_ruiz).map(|(&zh, &d)| zh * d).collect();
                    z_map.insert(p, z_unscaled);
                }
            }
            freq_records.push((fi, w, x_red, s_col, z_map));
        }
    }

    // Combine frequency records into scalar `objective` and per-frequency `DrivenFreqAdjointPair`
    let nf = freq_records.len().max(1) as f64;
    let mut aux_scalars = BTreeMap::new();
    let mut adj_pairs: Vec<DrivenFreqAdjointPair> = Vec::with_capacity(freq_records.len());

    let objective = match &req.objective {
        SensitivityObjective::Transmission { out_port, target, .. } => {
            let p = *out_port;
            let mean_mag_sq: f64 = freq_records.iter().map(|(_, _, _, s_col, _)| s_col[p].norm_sqr()).sum::<f64>() / nf;
            aux_scalars.insert("S_mag_sq".into(), mean_mag_sq);
            let (obj, chain) = match target {
                None => (mean_mag_sq, 1.0),
                Some(t) => {
                    let diff = mean_mag_sq - t;
                    (0.5 * diff * diff, diff)
                }
            };
            for (fi, w, x_red, s_col, z_map) in freq_records {
                let mu = s_col[p].conj() * (chain / nf);
                let z_p = &z_map[&p];
                let lam_red: Vec<C> = z_p.iter().map(|&zi| mu * zi).collect();
                adj_pairs.push(DrivenFreqAdjointPair {
                    freq_idx: fi,
                    omega: w,
                    in_port,
                    x_red,
                    lam_red,
                    port_weights: vec![(p, mu)],
                });
            }
            obj
        }
        SensitivityObjective::SplitterBalance { out_ports, targets, .. } => {
            let m_out = out_ports.len();
            let eff_targets: Vec<f64> =
                if targets.is_empty() { vec![1.0 / m_out as f64; m_out] } else { targets.clone() };
            let mut mean_mags = vec![0.0; m_out];
            for (_, _, _, s_col, _) in &freq_records {
                for (mi, &p) in out_ports.iter().enumerate() {
                    mean_mags[mi] += s_col[p].norm_sqr() / nf;
                }
            }
            let mut obj = 0.0;
            for (mi, &p) in out_ports.iter().enumerate() {
                aux_scalars.insert(format!("S_{p}_{in_port}_mag_sq"), mean_mags[mi]);
                obj += 0.5 * (mean_mags[mi] - eff_targets[mi]).powi(2);
            }
            for (fi, w, x_red, s_col, z_map) in freq_records {
                let mut lam_red = vec![c(0.0, 0.0); asm.nred];
                let mut pw = Vec::with_capacity(m_out);
                for (mi, &p) in out_ports.iter().enumerate() {
                    let mu = s_col[p].conj() * ((mean_mags[mi] - eff_targets[mi]) / nf);
                    let z_p = &z_map[&p];
                    for i in 0..asm.nred {
                        lam_red[i] += mu * z_p[i];
                    }
                    pw.push((p, mu));
                }
                adj_pairs.push(DrivenFreqAdjointPair {
                    freq_idx: fi,
                    omega: w,
                    in_port,
                    x_red,
                    lam_red,
                    port_weights: pw,
                });
            }
            obj
        }
        SensitivityObjective::FieldPower { domain, .. } => {
            let doms = ctx.resolver.domains(domain)?;
            let dom_cells: Vec<usize> =
                (0..simplex.n_cells()).filter(|&c| doms.contains(&simplex.cell_labels[c])).collect();
            let mut total_pow = 0.0;
            for (fi, w, x_red, _, z_map) in freq_records {
                let x_full = asm.red.expand(&x_red);
                let cell_pows = cell_3d_nedelec_mass_dots(&simplex, &asm.ned, &x_full, &x_full, true, &dom_cells);
                let pow_k: f64 = cell_pows.iter().map(|z| z.re).sum();
                total_pow += pow_k / nf;
                let lam_red: Vec<C> = z_map[&0].iter().map(|&zi| zi / nf).collect();
                adj_pairs.push(DrivenFreqAdjointPair {
                    freq_idx: fi,
                    omega: w,
                    in_port,
                    x_red,
                    lam_red,
                    port_weights: Vec::new(),
                });
            }
            aux_scalars.insert("field_power".into(), total_pow);
            total_pow
        }
        other => {
            return Err(Error::Unsupported(format!(
                "ems driven sensitivity supports Transmission, SplitterBalance, and FieldPower, got '{other:?}'"
            )))
        }
    };

    // Unscaled reduced operator contraction helper for scalar `MaterialScalar` / `ParamVar`
    let eval_driven_asm_contraction = |mp_opt: Option<(&Discretisation, &Discretisation)>,
                                       asm_p: &DrivenAssembly,
                                       asm_m: &DrivenAssembly,
                                       inv_2h: f64|
     -> f64 {
        let mut total = 0.0;
        for pair in &adj_pairs {
            let w = pair.omega;
            let ap = unscale_sym_matrix(&asm_p.assemble_matrix_at(w), &asm_p.d_ruiz);
            let am = unscale_sym_matrix(&asm_m.assemble_matrix_at(w), &asm_m.d_ruiz);
            let da = CCsr::axpby(c(inv_2h, 0.0), &ap, c(-inv_2h, 0.0), &am);
            let dax = da.matvec(&pair.x_red);

            // Drive RHS in unscaled reduced coordinates: `r_red = gamma_q * D^-1 drive_hat`
            let (_, gp, _, _) = asm_p.ports[pair.in_port].excitation_coeffs(w, &asm_p.waveports, &asm_p.floquet_ports);
            let (_, gm, _, _) = asm_m.ports[pair.in_port].excitation_coeffs(w, &asm_m.waveports, &asm_m.floquet_ports);
            let dp_hat = asm_p.ports[pair.in_port].drive_hat();
            let dm_hat = asm_m.ports[pair.in_port].drive_hat();

            let mut inner = c(0.0, 0.0);
            for i in 0..asm.nred {
                let r_p = gp * (dp_hat[i] / asm_p.d_ruiz[i]);
                let r_m = gm * (dm_hat[i] / asm_m.d_ruiz[i]);
                let dr_i = (r_p - r_m) * inv_2h;
                inner += pair.lam_red[i] * (dr_i - dax[i]);
            }
            // Add explicit port receive-functional derivative `d(z_wave_p^(-1/2) t_red_p^T)/dp x_red`
            for &(p, mu_pq) in &pair.port_weights {
                let (zwp, _, _, _) = asm_p.ports[p].excitation_coeffs(w, &asm_p.waveports, &asm_p.floquet_ports);
                let (zwm, _, _, _) = asm_m.ports[p].excitation_coeffs(w, &asm_m.waveports, &asm_m.floquet_ports);
                let y_on_p: Vec<C> = pair.x_red.iter().zip(&asm_p.d_ruiz).map(|(&x, &d)| x / d).collect();
                let y_on_m: Vec<C> = pair.x_red.iter().zip(&asm_m.d_ruiz).map(|(&x, &d)| x / d).collect();
                let sp_exp = asm_p.ports[p].recv_dot(&y_on_p) / zwp.max(1e-30).sqrt();
                let sm_exp = asm_m.ports[p].recv_dot(&y_on_m) / zwm.max(1e-30).sqrt();
                inner += mu_pq * (sp_exp - sm_exp) * inv_2h;
            }
            total += 2.0 * inner.re;

            // Explicit geometric variation of `FieldPower` domain volume if `ParamVar` deforms `domain`
            if let (Some((mp, mm)), Some(ref dom_cells)) = (mp_opt, &field_power_domain) {
                let x_full = asm.red.expand(&pair.x_red);
                let pp: f64 = cell_3d_nedelec_mass_dots(mp, &asm_p.ned, &x_full, &x_full, true, dom_cells)
                    .iter()
                    .map(|z| z.re)
                    .sum();
                let pm: f64 = cell_3d_nedelec_mass_dots(mm, &asm_m.ned, &x_full, &x_full, true, dom_cells)
                    .iter()
                    .map(|z| z.re)
                    .sum();
                total += ((pp - pm) * inv_2h) / nf;
            }
        }
        let _ = pair_unused_guard(&adj_pairs);
        total
    };

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..3).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());
    let mut projected_densities = BTreeMap::new();
    let mut volume_fractions = BTreeMap::new();
    let mut volume_fraction_gradients = BTreeMap::new();

    for (idx, param) in req.params.iter().enumerate() {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = base_problem
                    .materials
                    .get(doms[0], property)
                    .unwrap_or(if property == "eps_r" || property == "mu_r" { 1.0 } else { 0.0 });
                let h = 1e-6 * p0.abs().max(1.0);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = assemble_driven_for_sens(&simplex, &prob_p, &eps_overrides).map_err(Error::from_kernel)?;
                let asm_m = assemble_driven_for_sens(&simplex, &prob_m, &eps_overrides).map_err(Error::from_kernel)?;
                let g = eval_driven_asm_contraction(None, &asm_p, &asm_m, 1.0 / (2.0 * h));
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::ParamVar { name, .. } => {
                let Some(vel) = req.shape_velocities.get(name) else {
                    return Err(Error::invalid(format!("missing precomputed shape velocity for parameter '{name}'")));
                };
                let h = vel.optimal_step(char_len, 1e-6);
                let mp = displace_simplex_mesh(&simplex, h, vel).map_err(Error::from_kernel)?;
                let mm = displace_simplex_mesh(&simplex, -h, vel).map_err(Error::from_kernel)?;
                let asm_p = assemble_driven_for_sens(&mp, base_problem, &eps_overrides).map_err(Error::from_kernel)?;
                let asm_m = assemble_driven_for_sens(&mm, base_problem, &eps_overrides).map_err(Error::from_kernel)?;
                let g = eval_driven_asm_contraction(Some((&mp, &mm)), &asm_p, &asm_m, 1.0 / (2.0 * h));
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity {
                simp_penalty,
                ramp_q,
                val_min,
                val_max,
                val_min_im,
                val_max_im,
                plasmon_damping,
                ..
            } => {
                let st = &topo_states[&idx];
                let dcells = st.filter.design_cells();
                let nc_des = dcells.len();
                let mut dj_dbar = vec![0.0; nc_des];

                for pair in &adj_pairs {
                    let k0 = pair.omega / C0;
                    let k0_2 = k0 * k0;
                    let x_full = asm.red.expand(&pair.x_red);
                    let lam_full = expand_conj_prolong(&asm.red, &pair.lam_red);
                    let edots = cell_3d_nedelec_mass_dots(&simplex, &asm.ned, &lam_full, &x_full, false, dcells);
                    for (k, (&rb, &edot)) in st.rho_bar.iter().zip(&edots).enumerate() {
                        let (_, deps) = optical_interp_of_param(
                            rb,
                            *simp_penalty,
                            *ramp_q,
                            *val_min,
                            *val_max,
                            *val_min_im,
                            *val_max_im,
                            *plasmon_damping,
                        );
                        // A = K - k0^2 M => -dA/d(rho_bar) = +k0^2 * deps * M_c^(1)
                        dj_dbar[k] += 2.0 * (deps * edot * k0_2).re;
                    }
                }

                let g_rho = st.pullback(&dj_dbar)?;
                let (vf, vf_grad) = st.volume_fraction_and_grad()?;
                aux_scalars.insert(format!("volume_fraction_{idx}"), vf);
                projected_densities.insert(idx, st.rho_bar.clone());
                volume_fractions.insert(idx, vf);
                volume_fraction_gradients.insert(idx, vf_grad);
                flat_gradient.extend_from_slice(&g_rho);
                gradients.push(g_rho);
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to study 'eme'"));
            }
        }
    }

    let rep = SensitivityReport {
        objective,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities,
        volume_fractions,
        volume_fraction_gradients,
        mode_basis: None,
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}

#[inline]
fn pair_unused_guard(pairs: &[DrivenFreqAdjointPair]) -> usize {
    pairs.iter().map(|p| p.freq_idx).sum()
}

// ---------------------------------------------------------------------------
// 5. 3D Eigenmode Expansion (`eme`, ADR 140 + ADR 142)
// ---------------------------------------------------------------------------

fn eme_sensitivity(
    mesh: &Mesh,
    ctx: &Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let asm = crate::eme::assemble_eme(mesh, base_problem).map_err(Error::from_kernel)?;
    let cascade = crate::eme::solve_eme_cascade(&asm).map_err(Error::from_kernel)?;
    let s_term = &cascade.s_terminal;
    let p_ports = s_term.len();
    if p_ports == 0 {
        return Err(Error::invalid("eme sensitivity requires at least one modal or waveport channel"));
    }

    let mut aux_scalars = BTreeMap::new();
    let objective = match &req.objective {
        SensitivityObjective::Transmission { out_port, in_port, target, .. } => {
            if *in_port >= p_ports || *out_port >= p_ports {
                return Err(Error::invalid(format!(
                    "transmission port index out of range (in_port={in_port}, out_port={out_port}, ports={p_ports})"
                )));
            }
            let mag_sq = s_term[*out_port][*in_port].norm_sqr();
            aux_scalars.insert("S_mag_sq".into(), mag_sq);
            match target {
                None => mag_sq,
                Some(t) => {
                    let diff = mag_sq - t;
                    0.5 * diff * diff
                }
            }
        }
        SensitivityObjective::SplitterBalance { in_port, out_ports, targets, .. } => {
            if *in_port >= p_ports || out_ports.iter().any(|&p| p >= p_ports) {
                return Err(Error::invalid("splitter_balance port index out of range"));
            }
            let m_out = out_ports.len();
            let eff_targets: Vec<f64> =
                if targets.is_empty() { vec![1.0 / m_out as f64; m_out] } else { targets.clone() };
            let mut obj = 0.0;
            for (mi, &p) in out_ports.iter().enumerate() {
                let mag_sq = s_term[p][*in_port].norm_sqr();
                aux_scalars.insert(format!("S_{p}_{in_port}_mag_sq"), mag_sq);
                obj += 0.5 * (mag_sq - eff_targets[mi]).powi(2);
            }
            obj
        }
        SensitivityObjective::EffectiveIndex { mode, .. } => {
            let sm0 = cascade.slice_modes.first().ok_or_else(|| Error::invalid("EME cascade has zero slices"))?;
            let &neff = sm0.neff_fwd.get(*mode).ok_or_else(|| {
                Error::invalid(format!("requested mode {mode}, EME slice solved {} mode(s)", sm0.neff_fwd.len()))
            })?;
            let (obj, _, sc) = eval_neff_or_q_objective(&req.objective, neff)?;
            aux_scalars.extend(sc);
            obj
        }
        other => {
            return Err(Error::Unsupported(format!(
                "ems eme sensitivity supports Transmission, SplitterBalance, and EffectiveIndex, got '{other:?}'"
            )));
        }
    };

    let contract_ds = |ds_term: &[Vec<C>]| -> f64 {
        match &req.objective {
            SensitivityObjective::Transmission { out_port, in_port, target, .. } => {
                let s_pq = s_term[*out_port][*in_port];
                let d_mag_sq = 2.0 * (s_pq.conj() * ds_term[*out_port][*in_port]).re;
                match target {
                    None => d_mag_sq,
                    Some(t) => (s_pq.norm_sqr() - t) * d_mag_sq,
                }
            }
            SensitivityObjective::SplitterBalance { in_port, out_ports, targets, .. } => {
                let m_out = out_ports.len();
                let eff_targets: Vec<f64> =
                    if targets.is_empty() { vec![1.0 / m_out as f64; m_out] } else { targets.clone() };
                let mut g = 0.0;
                for (mi, &p) in out_ports.iter().enumerate() {
                    let s_pq = s_term[p][*in_port];
                    let d_mag_sq = 2.0 * (s_pq.conj() * ds_term[p][*in_port]).re;
                    g += (s_pq.norm_sqr() - eff_targets[mi]) * d_mag_sq;
                }
                g
            }
            _ => 0.0,
        }
    };

    let char_len = {
        let (lo, hi) = asm.mesh_3d.bbox();
        (0..3).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());

    for param in &req.params {
        param_names.push(param.label());
        match param {
            DesignParam::EmeSegmentLength { segment, .. } => {
                if *segment >= cascade.data.total_segments {
                    return Err(Error::invalid(format!(
                        "eme_segment_length segment {segment} out of range (cascade has {} segment(s))",
                        cascade.data.total_segments
                    )));
                }
                let g = if matches!(req.objective, SensitivityObjective::EffectiveIndex { .. }) {
                    0.0
                } else {
                    let mut d_lengths = vec![0.0; cascade.data.total_segments];
                    d_lengths[*segment] = 1.0;
                    let (_, ds_term) =
                        cascade.data.evaluate_terminal_s_and_grad(None, &d_lengths).map_err(Error::from_kernel)?;
                    contract_ds(&ds_term)
                };
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::ParamVar { name, .. } => {
                let Some(vel) = req.shape_velocities.get(name) else {
                    return Err(Error::invalid(format!("missing precomputed shape velocity for parameter '{name}'")));
                };
                // Detect pure longitudinal y-stretch on y-invariant cross-sections without y-varying imported fields:
                // evaluate exact analytical Redheffer derivative via dP_c/dL_c = -j beta_c P_c with zero extra eigensolves!
                let mut max_v = 0.0_f64;
                let mut max_xz_v = 0.0_f64;
                for v in &vel.vertex_vel {
                    let (vx, vy, vz) = (v[0].abs(), v[1].abs(), v[2].abs());
                    max_v = max_v.max(vx).max(vy).max(vz);
                    max_xz_v = max_xz_v.max(vx).max(vz);
                }
                for (_, ev) in &vel.edge_mid_vel {
                    max_v = max_v.max(ev[0].abs()).max(ev[1].abs()).max(ev[2].abs());
                    max_xz_v = max_xz_v.max(ev[0].abs()).max(ev[2].abs());
                }
                let pure_y_stretch = max_xz_v <= 1e-13 * max_v.max(1.0)
                    && cascade.data.unique_slices.len() == 1
                    && !base_problem.has_cell_material_perturbation()
                    && base_problem.imported_potential.is_none()
                    && base_problem.transferred_potential.is_none();

                let h = vel.optimal_step(char_len, 1e-6);
                let mp = femcore::topology::displace_mesh(mesh, h, vel).map_err(Error::from_kernel)?;
                let mm = femcore::topology::displace_mesh(mesh, -h, vel).map_err(Error::from_kernel)?;
                let inv_2h = 1.0 / (2.0 * h);

                let g = if pure_y_stretch {
                    if matches!(req.objective, SensitivityObjective::EffectiveIndex { .. }) {
                        0.0
                    } else {
                        let asm_p = crate::eme::assemble_eme(&mp, base_problem).map_err(Error::from_kernel)?;
                        let asm_m = crate::eme::assemble_eme(&mm, base_problem).map_err(Error::from_kernel)?;
                        let d_lengths: Vec<f64> = asm_p
                            .mesh_segments
                            .iter()
                            .zip(&asm_m.mesh_segments)
                            .map(|(sp, sm)| (sp.length_m - sm.length_m) * inv_2h)
                            .collect();
                        let (_, ds_term) =
                            cascade.data.evaluate_terminal_s_and_grad(None, &d_lengths).map_err(Error::from_kernel)?;
                        contract_ds(&ds_term)
                    }
                } else {
                    let asm_p = crate::eme::assemble_eme(&mp, base_problem).map_err(Error::from_kernel)?;
                    let cas_p = crate::eme::solve_eme_cascade(&asm_p).map_err(Error::from_kernel)?;
                    let asm_m = crate::eme::assemble_eme(&mm, base_problem).map_err(Error::from_kernel)?;
                    let cas_m = crate::eme::solve_eme_cascade(&asm_m).map_err(Error::from_kernel)?;
                    if let SensitivityObjective::EffectiveIndex { mode, .. } = &req.objective {
                        let neff0 = cascade.slice_modes[0].neff_fwd[*mode];
                        let (_, w_chain, _) = eval_neff_or_q_objective(&req.objective, neff0)?;
                        let dneff =
                            (cas_p.slice_modes[0].neff_fwd[*mode] - cas_m.slice_modes[0].neff_fwd[*mode]) * inv_2h;
                        w_chain.re * dneff.re + w_chain.im * dneff.im
                    } else {
                        let ds_term: Vec<Vec<C>> = (0..p_ports)
                            .map(|r| {
                                (0..p_ports)
                                    .map(|col| (cas_p.s_terminal[r][col] - cas_m.s_terminal[r][col]) * inv_2h)
                                    .collect()
                            })
                            .collect();
                        contract_ds(&ds_term)
                    }
                };
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = base_problem
                    .materials
                    .get(doms[0], property)
                    .unwrap_or(if property == "eps_r" || property == "mu_r" { 1.0 } else { 0.0 });
                let h = 1e-6 * p0.abs().max(1.0);
                let inv_2h = 1.0 / (2.0 * h);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = crate::eme::assemble_eme(mesh, &prob_p).map_err(Error::from_kernel)?;
                let cas_p = crate::eme::solve_eme_cascade(&asm_p).map_err(Error::from_kernel)?;
                let asm_m = crate::eme::assemble_eme(mesh, &prob_m).map_err(Error::from_kernel)?;
                let cas_m = crate::eme::solve_eme_cascade(&asm_m).map_err(Error::from_kernel)?;
                let g = if let SensitivityObjective::EffectiveIndex { mode, .. } = &req.objective {
                    let neff0 = cascade.slice_modes[0].neff_fwd[*mode];
                    let (_, w_chain, _) = eval_neff_or_q_objective(&req.objective, neff0)?;
                    let dneff = (cas_p.slice_modes[0].neff_fwd[*mode] - cas_m.slice_modes[0].neff_fwd[*mode]) * inv_2h;
                    w_chain.re * dneff.re + w_chain.im * dneff.im
                } else {
                    let ds_term: Vec<Vec<C>> = (0..p_ports)
                        .map(|r| {
                            (0..p_ports)
                                .map(|col| (cas_p.s_terminal[r][col] - cas_m.s_terminal[r][col]) * inv_2h)
                                .collect()
                        })
                        .collect();
                    contract_ds(&ds_term)
                };
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity { .. } => {
                return Err(Error::Unsupported(
                    "topology_density on Study::Eme is not applicable (EME segments are y-invariant); use driven or modes for cellwise topology optimization, or eme_segment_length / param_var / material_scalar on eme".into(),
                ));
            }
        }
    }

    let rep = SensitivityReport {
        objective,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities: BTreeMap::new(),
        volume_fractions: BTreeMap::new(),
        volume_fraction_gradients: BTreeMap::new(),
        mode_basis: None,
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}
