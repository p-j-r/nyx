// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Discrete adjoint and eigenvalue sensitivities for `nyx-sms` (ADR 142):
//!
//! * `static`: exact discrete compliance sensitivity (`f^T u`, self-adjoint `-u^T dK/dp u + 2 df/dp^T u`
//!   via `femcore::linalg::opt::adjoint_gradient` and cached Cholesky solves) with respect to CAD shape
//!   parameters (`ParamVar`), isotropic elastic moduli (`MaterialScalar`), and Helmholtz-filtered SIMP/RAMP
//!   stiffness topology densities (`TopologyDensity`).
//! * `eigenfrequency`, `modes` (2D elastic waveguide), and `axisymmetric_eigenfrequency`: exact modal
//!   eigenfrequency `Re(f)` and quality-factor `Q` sensitivities via `femcore::linalg::sens::{eigenvalue_derivative,
//!   left_eigenvector, match_modes}` across real symmetric, PML complex-symmetric, and Floquet-reduced pencils.

use crate::cross_section::{assemble_cross_section, solve_cross_section_modes, CrossSection};
use crate::problem::{decode, Problem, Study};
use crate::solid::{assemble_solid, assemble_static, b_matrix, cell_rho, d_matrix, dof_layout, setup};
use crate::Sms;
use femcore::fem::{assemble_mass, free_dofs, geometry, shape_gradients, RefElement};
use femcore::linalg::csr::{c, CCsr, RCsr as Csr, C};
use femcore::linalg::direct::Definiteness;
use femcore::linalg::driver::{eigs_sym, DirectSolver, EigOpts};
use femcore::linalg::opt::adjoint_gradient;
use femcore::linalg::sens::{eigenvalue_derivative, left_eigenvector, match_modes};
use femcore::material::M6;
use femcore::topology::{displace_simplex_mesh, ramp_interpolate, simp_interpolate};
use nyx_physics::{
    check_mesh, discretise, validate_spec, Ctx, DesignParam, Discretisation, Error, EvaluatedTopologyParam, Mesh,
    Physics, ProblemSpec, SensitivityObjective, SensitivityReport, SensitivityRequest, TrackedModeBasis,
};
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Evaluates `Sms::sensitivity` on `mesh`.
pub fn compute_sensitivity(
    physics: &Sms,
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

    let base_problem = decode(spec, ctx)?;

    match base_problem.study {
        Study::Static => static_sensitivity(mesh, ctx, &base_problem, req),
        Study::Eigenfrequency { .. } => eigen_solid_sensitivity(mesh, ctx, &base_problem, req),
        Study::Modes { wavevector, .. } => {
            cross_section_sensitivity(mesh, ctx, &base_problem, req, CrossSection::Waveguide(wavevector))
        }
        Study::AxisymmetricEigenfrequency { azimuthal_order, .. } => {
            cross_section_sensitivity(mesh, ctx, &base_problem, req, CrossSection::Axisymmetric(azimuthal_order as f64))
        }
        ref other => {
            Err(Error::Unsupported(format!("sms study '{}' does not support adjoint sensitivities", other.kind())))
        }
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
fn eval_scalar_interp(rho_bar: f64, val_min: f64, val_max: f64, simp_penalty: f64, ramp_q: Option<f64>) -> (f64, f64) {
    match ramp_q {
        Some(q) => ramp_interpolate(rho_bar, val_min, val_max, q),
        None => simp_interpolate(rho_bar, val_min, val_max, simp_penalty),
    }
}

fn scale_m6(m: &M6, s: f64) -> M6 {
    let mut out = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            out[i][j] = s * m[i][j];
        }
    }
    out
}

/// Elementwise unscaled elastic strain energy bilinear form `u_c^T K_c(base_c6[c]) v_c`.
fn cell_elastic_strain_dots(
    mesh: &Discretisation,
    order: usize,
    plane: crate::problem::Plane,
    thickness: f64,
    base_c6: &[M6],
    u: &[f64],
    v: &[f64],
    cells: &[usize],
) -> Vec<f64> {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let thk = if d == 2 { thickness } else { 1.0 };
    const CHUNK: usize = 256;
    let parts = femcore::par::map_chunks(cells, CHUNK, |slice| {
        let mut out = Vec::with_capacity(slice.len());
        let (mut scratch, mut full) = (Vec::new(), Vec::new());
        let mut grad = vec![0.0; nloc * d];
        let mut b = vec![0.0; ns * nd];
        let mut eu = [0.0; 6];
        let mut ev = [0.0; 6];
        for &c in slice {
            let dm = d_matrix(&base_c6[c], d, plane);
            dof_layout(mesh, order, c, &mut scratch, &mut full);
            let mut acc = 0.0;
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                b_matrix(&grad, nloc, d, &mut b);
                for s in 0..ns {
                    let mut su = 0.0;
                    let mut sv = 0.0;
                    for j in 0..nd {
                        let b_sj = b[s * nd + j];
                        su += b_sj * u[full[j]];
                        sv += b_sj * v[full[j]];
                    }
                    eu[s] = su;
                    ev[s] = sv;
                }
                let mut quad_q = 0.0;
                for s in 0..ns {
                    for t in 0..ns {
                        quad_q += eu[s] * dm[s * ns + t] * ev[t];
                    }
                }
                acc += re.w[q] * quad_q;
            }
            out.push(geo.vol[c] * thk * acc);
        }
        out
    });
    parts.into_iter().flatten().collect()
}

/// Elementwise complex elastic bilinear form `y_c^H K_c(base_c6[c]) x_c` and mass `y_c^H M_c(base_rho[c]) x_c`.
fn cell_complex_modal_dots(
    mesh: &Discretisation,
    order: usize,
    plane: crate::problem::Plane,
    thickness: f64,
    base_c6: &[M6],
    base_rho: &[f64],
    y_full: &[C],
    x_full: &[C],
    cells: &[usize],
) -> Vec<(C, C)> {
    let d = mesh.dim;
    let re_k = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let re_m = RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re_k.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let thk = if d == 2 { thickness } else { 1.0 };
    const CHUNK: usize = 256;
    let parts = femcore::par::map_chunks(cells, CHUNK, |slice| {
        let mut out = Vec::with_capacity(slice.len());
        let (mut scratch, mut full) = (Vec::new(), Vec::new());
        let mut grad = vec![0.0; nloc * d];
        let mut b = vec![0.0; ns * nd];
        let mut ey = [c(0.0, 0.0); 6];
        let mut ex = [c(0.0, 0.0); 6];
        for &ci in slice {
            let dm = d_matrix(&base_c6[ci], d, plane);
            dof_layout(mesh, order, ci, &mut scratch, &mut full);
            let mut k_acc = c(0.0, 0.0);
            for q in 0..re_k.w.len() {
                shape_gradients(&re_k, &geo, ci, q, &mut grad);
                b_matrix(&grad, nloc, d, &mut b);
                for s in 0..ns {
                    let mut sy = c(0.0, 0.0);
                    let mut sx = c(0.0, 0.0);
                    for j in 0..nd {
                        let b_sj = b[s * nd + j];
                        sy += y_full[full[j]] * b_sj;
                        sx += x_full[full[j]] * b_sj;
                    }
                    ey[s] = sy;
                    ex[s] = sx;
                }
                let mut quad_q = c(0.0, 0.0);
                for s in 0..ns {
                    for t in 0..ns {
                        quad_q += ey[s].conj() * ex[t] * dm[s * ns + t];
                    }
                }
                k_acc += quad_q * re_k.w[q];
            }
            let mut m_acc = c(0.0, 0.0);
            let n1 = mesh.n_dofs(order);
            for q in 0..re_m.w.len() {
                let mut dot_q = c(0.0, 0.0);
                for comp in 0..d {
                    let mut uq_y = c(0.0, 0.0);
                    let mut uq_x = c(0.0, 0.0);
                    for a in 0..nloc {
                        let na = re_m.n[q][a];
                        let gd = comp * n1 + scratch[a];
                        uq_y += y_full[gd] * na;
                        uq_x += x_full[gd] * na;
                    }
                    dot_q += uq_y.conj() * uq_x;
                }
                m_acc += dot_q * re_m.w[q];
            }
            let vol_t = geo.vol[ci] * thk;
            out.push((k_acc * vol_t, m_acc * (vol_t * base_rho[ci])));
        }
        out
    });
    parts.into_iter().flatten().collect()
}

fn static_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let SensitivityObjective::Compliance { target } = &req.objective else {
        return Err(Error::Unsupported(format!("sms static supports Compliance objective, got '{:?}'", req.objective)));
    };
    let simplex = discretise(mesh)?;
    let order = base_problem.order;
    let topo_states = build_topology_states(&simplex, ctx, req)?;

    let assemble_with_topo =
        |m_cur: &Discretisation, prob: &Problem| -> Result<(crate::solid::Setup, Csr, Vec<f64>, Vec<M6>), String> {
            let mut st = setup(m_cur, prob)?;
            let base_c6 = st.c6.clone();
            for (idx, param) in req.params.iter().enumerate() {
                if let DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, .. } = param {
                    let ts = &topo_states[&idx];
                    for (&c_idx, &rb) in ts.filter.design_cells().iter().zip(&ts.rho_bar) {
                        let (scale, _) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                        st.c6[c_idx] = scale_m6(&base_c6[c_idx], scale);
                    }
                }
            }
            let k = assemble_solid(m_cur, &st.c6, order, prob.plane, prob.thickness);
            let asm = assemble_static(m_cur, prob, &st, k)?;
            Ok((st, asm.k, asm.f, base_c6))
        };

    let (st, k_mat, f_vec, base_c6) = assemble_with_topo(&simplex, base_problem).map_err(Error::from_kernel)?;
    let n = k_mat.nrows();
    let free = free_dofs(n, &st.fixed);
    let mut u_full = vec![0.0; n];
    for (&d, &v) in st.fixed.iter().zip(&st.vals) {
        u_full[d] = v;
    }
    let mut lam_full = vec![0.0; n];

    let comp = if st.vals.iter().all(|&v| v == 0.0) && !free.is_empty() {
        // Homogeneous Dirichlet static compliance `f^T u`: wire `femcore::linalg::opt::adjoint_gradient`
        let kff = k_mat.submatrix(&free, &free);
        let ff: Vec<f64> = free.iter().map(|&gi| f_vec[gi]).collect();
        let (c_val, _, uf) =
            adjoint_gradient(&kff, &ff, &ff, 0, |_, _| Vec::new(), |_| None, Definiteness::PositiveDefinite)
                .map_err(|e| Error::Solver(e.to_string()))?;
        for (a, &gi) in free.iter().enumerate() {
            u_full[gi] = uf[a];
            lam_full[gi] = uf[a];
        }
        c_val
    } else {
        if !free.is_empty() {
            let kff = k_mat.submatrix(&free, &free);
            let kfx = k_mat.submatrix(&free, &st.fixed);
            let kx = kfx.matvec(&st.vals);
            let rhs: Vec<f64> = free.iter().enumerate().map(|(a, &gi)| f_vec[gi] - kx[a]).collect();
            let solver = DirectSolver::new(&kff, true).map_err(Error::from_kernel)?;
            let uf = solver.solve(&rhs);
            let gf: Vec<f64> = free.iter().map(|&gi| f_vec[gi]).collect();
            let lf = solver.solve(&gf);
            for (a, &gi) in free.iter().enumerate() {
                u_full[gi] = uf[a];
                lam_full[gi] = lf[a];
            }
        }
        f_vec.iter().zip(&u_full).map(|(fi, ui)| fi * ui).sum()
    };

    let mut aux_scalars = BTreeMap::new();
    aux_scalars.insert("compliance".into(), comp);
    let (objective, chain_mult) = match target {
        Some(t) => (0.5 * (comp - t).powi(2), comp - t),
        None => (comp, 1.0),
    };

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..simplex.dim).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());
    let mut projected_densities = BTreeMap::new();
    let mut volume_fractions = BTreeMap::new();
    let mut volume_fraction_gradients = BTreeMap::new();

    // For J = f^T u with K_{ff} u_f = f_f - K_{fd} u_d (u_d fixed):
    // d(f^T u)/dp = (df/dp)^T u + lam_f^T (df_f/dp - (dK_{ff}/dp u_f + dK_{fd}/dp u_d))
    //             = (df/dp)^T u + lam_full^T (df/dp - (dK/dp) u_full).
    let contract_static = |kp: &Csr, km: &Csr, fp: &[f64], fm: &[f64], inv_2h: f64| -> f64 {
        let dk = Csr::axpby(inv_2h, kp, -inv_2h, km);
        let dku = dk.matvec(&u_full);
        let mut acc = 0.0;
        for i in 0..n {
            let df_i = (fp[i] - fm[i]) * inv_2h;
            acc += df_i * u_full[i] + lam_full[i] * (df_i - dku[i]);
        }
        chain_mult * acc
    };

    for (idx, param) in req.params.iter().enumerate() {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = base_problem.materials.get(doms[0], property).ok_or_else(|| {
                    Error::invalid(format!("material property '{property}' not set on region '{region}'"))
                })?;
                let h = 3e-5 * p0.abs().max(1.0);
                let inv_2h = 1.0 / (2.0 * h);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let (st_p, kp, fp, _) = assemble_with_topo(&simplex, &prob_p).map_err(Error::from_kernel)?;
                let (st_m, km, fm, _) = assemble_with_topo(&simplex, &prob_m).map_err(Error::from_kernel)?;
                let load_unchanged =
                    fp.iter().zip(&fm).all(|(a, b)| (a - b).abs() == 0.0) && base_problem.infinite_elements.is_empty();
                let g = if load_unchanged {
                    // Evaluate directly in element strain space `(B u)^T dC (B lam)` so rigid-body
                    // displacements cancel identically before constitutive differentiation.
                    let all_cells: Vec<usize> = (0..simplex.n_cells()).collect();
                    let dc6: Vec<M6> = st_p
                        .c6
                        .iter()
                        .zip(&st_m.c6)
                        .map(|(cp, cm)| {
                            let mut d6 = [[0.0; 6]; 6];
                            for r in 0..6 {
                                for c_i in 0..6 {
                                    d6[r][c_i] = (cp[r][c_i] - cm[r][c_i]) * inv_2h;
                                }
                            }
                            d6
                        })
                        .collect();
                    let edots = cell_elastic_strain_dots(
                        &simplex,
                        order,
                        base_problem.plane,
                        base_problem.thickness,
                        &dc6,
                        &lam_full,
                        &u_full,
                        &all_cells,
                    );
                    -chain_mult * edots.iter().sum::<f64>()
                } else {
                    contract_static(&kp, &km, &fp, &fm, inv_2h)
                };
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
                let (_, kp, fp, _) = assemble_with_topo(&mp, base_problem).map_err(Error::from_kernel)?;
                let (_, km, fm, _) = assemble_with_topo(&mm, base_problem).map_err(Error::from_kernel)?;
                let g = contract_static(&kp, &km, &fp, &fm, 1.0 / (2.0 * h));
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, .. } => {
                let ts = &topo_states[&idx];
                let dcells = ts.filter.design_cells();
                let edots = cell_elastic_strain_dots(
                    &simplex,
                    order,
                    base_problem.plane,
                    base_problem.thickness,
                    &base_c6,
                    &lam_full,
                    &u_full,
                    dcells,
                );
                let dj_dbar: Vec<f64> = ts
                    .rho_bar
                    .iter()
                    .zip(&edots)
                    .map(|(&rb, &edot)| {
                        let (_, dw) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                        chain_mult * (-dw * edot)
                    })
                    .collect();
                let g_rho = ts.pullback(&dj_dbar)?;
                let (vf, vf_grad) = ts.volume_fraction_and_grad()?;
                aux_scalars.insert(format!("volume_fraction_{idx}"), vf);
                projected_densities.insert(idx, ts.rho_bar.clone());
                volume_fractions.insert(idx, vf);
                volume_fraction_gradients.insert(idx, vf_grad);
                flat_gradient.extend_from_slice(&g_rho);
                gradients.push(g_rho);
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to ems 'eme' studies"));
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

/// Converts `(f, df/dp)` in complex Hz into objective value and scalar derivative `dJ/dp`.
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
            if *inverse {
                let inv_q = 2.0 * im / re.abs().max(1e-30);
                scalars.insert("inv_q".into(), inv_q);
                let sgn_re = if re >= 0.0 { 1.0 } else { -1.0 };
                let d_re = -2.0 * im * sgn_re / (re * re).max(1e-30);
                let d_im = 2.0 * sgn_im / re.abs().max(1e-30);
                Ok((inv_q, c(d_re, d_im), scalars))
            } else {
                if im <= 1e-30 * re.abs() {
                    return Err(Error::invalid(
                        "QualityFactor sensitivity requires a lossy/PML mode with non-zero Im(f) (or set inverse = true)",
                    ));
                }
                let q = re.abs() / (2.0 * im);
                scalars.insert("Q".into(), q);
                let sgn_re = if re >= 0.0 { 1.0 } else { -1.0 };
                let d_re = sgn_re / (2.0 * im);
                let d_im = -re.abs() * sgn_im / (2.0 * im * im);
                Ok((q, c(d_re, d_im), scalars))
            }
        }
        other => Err(Error::Unsupported(format!(
            "sms modal sensitivity supports Eigenfrequency and QualityFactor, got '{other:?}'"
        ))),
    }
}

fn target_mode_index(obj: &SensitivityObjective) -> usize {
    match obj {
        SensitivityObjective::Eigenfrequency { mode, .. } | SensitivityObjective::QualityFactor { mode, .. } => *mode,
        _ => 0,
    }
}

fn eigen_solid_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let Study::Eigenfrequency { n_modes, f_search } = base_problem.study else {
        return Err(Error::invalid("expected Eigenfrequency study"));
    };
    if base_problem.materials.has_pml() || !base_problem.periodic.is_empty() {
        return Err(Error::Unsupported(
            "3D/2D solid eigenfrequency sensitivity with PML or periodic BCs is supported on cross_section ('modes' / 'axisymmetric_eigenfrequency') studies"
                .into(),
        ));
    }

    let simplex = discretise(mesh)?;
    let order = base_problem.order;
    let d = simplex.dim;
    let thk = if d == 2 { base_problem.thickness } else { 1.0 };
    let topo_states = build_topology_states(&simplex, ctx, req)?;

    let assemble_km =
        |m_cur: &Discretisation, prob: &Problem| -> Result<(Vec<usize>, Csr, Csr, Vec<M6>, Vec<f64>), String> {
            let mut st = setup(m_cur, prob)?;
            let base_c6 = st.c6.clone();
            let base_rho = cell_rho(&prob.materials, m_cur, &prob.inactive)?;
            let mut rho = base_rho.clone();
            for (idx, param) in req.params.iter().enumerate() {
                if let DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, .. } = param {
                    let ts = &topo_states[&idx];
                    for (&c_idx, &rb) in ts.filter.design_cells().iter().zip(&ts.rho_bar) {
                        let (k_scale, _) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                        let (m_scale, _) = simp_interpolate(rb, (*val_min).max(1e-6), *val_max, 1.0);
                        st.c6[c_idx] = scale_m6(&base_c6[c_idx], k_scale);
                        rho[c_idx] = base_rho[c_idx] * m_scale;
                    }
                }
            }
            let k = assemble_solid(m_cur, &st.c6, order, prob.plane, prob.thickness);
            let rho_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
            let m1 = assemble_mass(m_cur, order, &rho_t);
            let blocks: Vec<&Csr> = (0..d).map(|_| &m1).collect();
            let m = Csr::block_diag(&blocks);
            let free = free_dofs(k.nrows(), &st.fixed);
            Ok((free, k, m, base_c6, base_rho))
        };

    let (free, k_full, m_full, base_c6, base_rho) = assemble_km(&simplex, base_problem).map_err(Error::from_kernel)?;
    let mode_idx = target_mode_index(&req.objective);
    let ask_modes = n_modes.max(mode_idx + 1);
    if ask_modes + 1 >= free.len() {
        return Err(Error::invalid("requested mode count exceeds free degrees of freedom"));
    }
    let kff = k_full.submatrix(&free, &free);
    let mff = m_full.submatrix(&free, &free);
    let sigma = if f_search > 0.0 { (2.0 * PI * f_search).powi(2) } else { -1.0 };
    let pairs = eigs_sym(&kff, &mff, sigma, ask_modes, &EigOpts::default(), None).map_err(Error::from_kernel)?;

    let mut order_idx: Vec<usize> = (0..pairs.values.len()).collect();
    order_idx.sort_by(|&a, &b| pairs.values[a].total_cmp(&pairs.values[b]));
    let mut evals: Vec<C> = order_idx.iter().map(|&i| c(pairs.values[i], 0.0)).collect();
    let mut evecs: Vec<Vec<C>> =
        order_idx.iter().map(|&i| pairs.vectors[i].iter().map(|&v| c(v, 0.0)).collect()).collect();

    let mf_c = mff.to_complex();
    if let Some(ref ref_basis) = req.reference_modes {
        if ref_basis.eigenvectors.len() == evecs.len()
            && ref_basis.eigenvectors.first().is_some_and(|v| v.len() == free.len())
        {
            let matched = match_modes(&ref_basis.eigenvectors, &evecs, Some(&mf_c), 0.1, 0.05);
            let old_evals = evals.clone();
            let old_evecs = evecs.clone();
            for (old_k, &new_j) in matched.assignment.iter().enumerate() {
                evals[old_k] = old_evals[new_j];
                evecs[old_k] = old_evecs[new_j].clone();
            }
        }
    }

    if mode_idx >= evals.len() {
        return Err(Error::invalid(format!("requested mode {mode_idx} not found (converged {} modes)", evals.len())));
    }
    let lambda = evals[mode_idx];
    let x_free = &evecs[mode_idx];
    let y_free = x_free.clone(); // Real symmetric `(K, M)` -> `y = x`
    let f_hz = c(lambda.re.max(0.0).sqrt() / (2.0 * PI), 0.0);
    if f_hz.re <= 1e-12 {
        return Err(Error::invalid(
            "target mode is a zero-frequency rigid-body mode; add fixed boundary constraints or select a higher mode",
        ));
    }

    let (objective, df_weight, mut aux_scalars) = eval_freq_or_q_objective(&req.objective, f_hz)?;
    // lambda = (2 pi f)^2 => df/dp = dlambda/dp / (8 pi^2 f)
    let dlam_to_df = c(1.0 / (8.0 * PI * PI * f_hz.re), 0.0);

    let mut x_full = vec![c(0.0, 0.0); k_full.nrows()];
    for (a, &gi) in free.iter().enumerate() {
        x_full[gi] = x_free[a];
    }
    let mx_f = mf_c.matvec(x_free);
    let denom_ymx: C = y_free.iter().zip(&mx_f).map(|(&yi, &mxi)| yi.conj() * mxi).sum();

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..simplex.dim).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
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
                let p0 = base_problem.materials.get(doms[0], property).ok_or_else(|| {
                    Error::invalid(format!("material property '{property}' not set on region '{region}'"))
                })?;
                let h = 1e-6 * p0.abs().max(1e-6);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &dom in &doms {
                    prob_p.materials.domains.entry(dom).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(dom).or_default().insert(property.clone(), p0 - h);
                }
                let (_, kp, mp_mat, _, _) = assemble_km(&simplex, &prob_p).map_err(Error::from_kernel)?;
                let (_, km, mm_mat, _, _) = assemble_km(&simplex, &prob_m).map_err(Error::from_kernel)?;
                let inv_2h = 1.0 / (2.0 * h);
                let dk_f = kp
                    .submatrix(&free, &free)
                    .add_scaled(&km.submatrix(&free, &free), -1.0)
                    .scaled(inv_2h)
                    .to_complex();
                let dm_f = mp_mat
                    .submatrix(&free, &free)
                    .add_scaled(&mm_mat.submatrix(&free, &free), -1.0)
                    .scaled(inv_2h)
                    .to_complex();
                let dlam = eigenvalue_derivative(lambda, x_free, &y_free, Some(&mf_c), &dk_f, Some(&dm_f));
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
                let (_, kp, mp_mat, _, _) = assemble_km(&mp, base_problem).map_err(Error::from_kernel)?;
                let (_, km, mm_mat, _, _) = assemble_km(&mm, base_problem).map_err(Error::from_kernel)?;
                let inv_2h = 1.0 / (2.0 * h);
                let dk_f = kp
                    .submatrix(&free, &free)
                    .add_scaled(&km.submatrix(&free, &free), -1.0)
                    .scaled(inv_2h)
                    .to_complex();
                let dm_f = mp_mat
                    .submatrix(&free, &free)
                    .add_scaled(&mm_mat.submatrix(&free, &free), -1.0)
                    .scaled(inv_2h)
                    .to_complex();
                let dlam = eigenvalue_derivative(lambda, x_free, &y_free, Some(&mf_c), &dk_f, Some(&dm_f));
                let df = dlam * dlam_to_df;
                let g = df_weight.re * df.re + df_weight.im * df.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, .. } => {
                let ts = &topo_states[&idx];
                let dcells = ts.filter.design_cells();
                let km_dots = cell_complex_modal_dots(
                    &simplex,
                    order,
                    base_problem.plane,
                    base_problem.thickness,
                    &base_c6,
                    &base_rho,
                    &x_full,
                    &x_full,
                    dcells,
                );
                let dj_dbar: Vec<f64> = ts
                    .rho_bar
                    .iter()
                    .zip(&km_dots)
                    .map(|(&rb, &(k_dot, m_dot))| {
                        let (_, dk_w) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                        let (_, dm_w) = simp_interpolate(rb, (*val_min).max(1e-6), *val_max, 1.0);
                        let num = k_dot * dk_w - lambda * m_dot * dm_w;
                        let dlam = num / denom_ymx;
                        let df = dlam * dlam_to_df;
                        df_weight.re * df.re + df_weight.im * df.im
                    })
                    .collect();
                let g_rho = ts.pullback(&dj_dbar)?;
                let (vf, vf_grad) = ts.volume_fraction_and_grad()?;
                aux_scalars.insert(format!("volume_fraction_{idx}"), vf);
                projected_densities.insert(idx, ts.rho_bar.clone());
                volume_fractions.insert(idx, vf);
                volume_fraction_gradients.insert(idx, vf_grad);
                flat_gradient.extend_from_slice(&g_rho);
                gradients.push(g_rho);
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to ems 'eme' studies"));
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
        mode_basis: Some(TrackedModeBasis { eigenvalues: evals, eigenvectors: evecs, mass_matrix: Some(mf_c) }),
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}

fn cross_section_sensitivity(
    mesh: &Mesh,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
    kind: CrossSection,
) -> Result<SensitivityReport, Error> {
    if req.params.iter().any(|p| matches!(p, DesignParam::TopologyDensity { .. })) {
        return Err(Error::Unsupported(
            "TopologyDensity on sms cross-section waveguide/axisymmetric modes is not yet enabled; use MaterialScalar or ParamVar"
                .into(),
        ));
    }
    let simplex = discretise(mesh)?;
    let asm = assemble_cross_section(&simplex, base_problem, kind).map_err(Error::from_kernel)?;
    if asm.cf.is_some() || asm.cf_complex.is_some() {
        return Err(Error::Unsupported(
            "cross-section sensitivity with quadratic low-reflecting boundary damping is not supported yet; use PML absorbing layers"
                .into(),
        ));
    }
    let solved = solve_cross_section_modes(&asm).map_err(Error::from_kernel)?;
    let mode_idx = target_mode_index(&req.objective);

    // Reconstruct reduced eigenvectors from `solved.modes.vectors`
    let n_conv = solved.modes.freqs.len();
    if mode_idx >= n_conv {
        return Err(Error::invalid(format!("requested mode {mode_idx} out of range (converged {n_conv} modes)")));
    }
    let mut freqs = solved.modes.freqs.clone();
    let mut red_vecs: Vec<Vec<C>> = solved
        .modes
        .vectors
        .iter()
        .map(|u_full| {
            if let Some(ref red) = asm.red {
                red.roots.iter().map(|&r| u_full[r]).collect()
            } else {
                asm.free.iter().map(|&gi| u_full[gi]).collect()
            }
        })
        .collect();

    if let Some(ref ref_basis) = req.reference_modes {
        if ref_basis.eigenvectors.len() == red_vecs.len()
            && ref_basis.eigenvectors.first().is_some_and(|v| v.len() == asm.kf.nrows())
        {
            let matched = match_modes(&ref_basis.eigenvectors, &red_vecs, Some(&asm.mf), 0.1, 0.05);
            let old_f = freqs.clone();
            let old_v = red_vecs.clone();
            for (old_k, &new_j) in matched.assignment.iter().enumerate() {
                freqs[old_k] = old_f[new_j];
                red_vecs[old_k] = old_v[new_j].clone();
            }
        }
    }

    let f_hz = freqs[mode_idx];
    let omega = f_hz * (2.0 * PI);
    let lambda = omega * omega;
    let x_red = &red_vecs[mode_idx];
    let y_red = if asm.kf.is_hermitian(1e-10) && asm.mf.is_hermitian(1e-10) {
        x_red.clone()
    } else if asm.kf.is_symmetric(1e-10) && asm.mf.is_symmetric(1e-10) {
        x_red.iter().map(|z| z.conj()).collect()
    } else {
        left_eigenvector(&asm.kf, Some(&asm.mf), lambda, x_red, 4).map_err(|e| Error::Solver(e.to_string()))?
    };

    let (objective, df_weight, aux_scalars) = eval_freq_or_q_objective(&req.objective, f_hz)?;
    let dlam_to_df = c(1.0, 0.0) / (f_hz * (8.0 * PI * PI));

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..simplex.dim).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());

    for param in &req.params {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = base_problem.materials.get(doms[0], property).ok_or_else(|| {
                    Error::invalid(format!("material property '{property}' not set on region '{region}'"))
                })?;
                let h = 1e-6 * p0.abs().max(1e-6);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &dom in &doms {
                    prob_p.materials.domains.entry(dom).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(dom).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = assemble_cross_section(&simplex, &prob_p, kind).map_err(Error::from_kernel)?;
                let asm_m = assemble_cross_section(&simplex, &prob_m, kind).map_err(Error::from_kernel)?;
                let inv_2h = c(1.0 / (2.0 * h), 0.0);
                let dk = CCsr::axpby(inv_2h, &asm_p.kf, -inv_2h, &asm_m.kf);
                let dm = CCsr::axpby(inv_2h, &asm_p.mf, -inv_2h, &asm_m.mf);
                let dlam = eigenvalue_derivative(lambda, x_red, &y_red, Some(&asm.mf), &dk, Some(&dm));
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
                let asm_p = assemble_cross_section(&mp, base_problem, kind).map_err(Error::from_kernel)?;
                let asm_m = assemble_cross_section(&mm, base_problem, kind).map_err(Error::from_kernel)?;
                let inv_2h = c(1.0 / (2.0 * h), 0.0);
                let dk = CCsr::axpby(inv_2h, &asm_p.kf, -inv_2h, &asm_m.kf);
                let dm = CCsr::axpby(inv_2h, &asm_p.mf, -inv_2h, &asm_m.mf);
                let dlam = eigenvalue_derivative(lambda, x_red, &y_red, Some(&asm.mf), &dk, Some(&dm));
                let df = dlam * dlam_to_df;
                let g = df_weight.re * df.re + df_weight.im * df.im;
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity { .. } => {
                return Err(Error::Unsupported(
                    "topology_density is not supported on 2D sms cross-section modes".into(),
                ));
            }
            DesignParam::EmeSegmentLength { .. } => {
                return Err(Error::invalid("eme_segment_length applies only to ems 'eme' studies"));
            }
        }
    }

    let evals: Vec<C> = freqs.iter().map(|&f| (f * (2.0 * PI)) * (f * (2.0 * PI))).collect();
    let rep = SensitivityReport {
        objective,
        gradients,
        flat_gradient,
        param_names,
        aux_scalars,
        projected_densities: BTreeMap::new(),
        volume_fractions: BTreeMap::new(),
        volume_fraction_gradients: BTreeMap::new(),
        mode_basis: Some(TrackedModeBasis { eigenvalues: evals, eigenvectors: red_vecs, mass_matrix: Some(asm.mf) }),
        log: Vec::new(),
    };
    rep.check_finite()?;
    Ok(rep)
}
