// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Discrete adjoint sensitivities for `nyx-dfs` (`electrostatic` and `heat_stationary`, ADR 142):
//!
//! * `electrostatic`: exact discrete derivatives of multi-terminal Maxwell capacitance `C_{ij}`,
//!   two-conductor capacitance `C = 2 W_e / (dV)^2`, and stored dielectric energy / compliance
//!   with respect to CAD shape parameters (`ParamVar`), scalar permittivities/charges (`MaterialScalar`),
//!   and Helmholtz-filtered SIMP/RAMP topology densities (`TopologyDensity`).
//! * `heat_stationary`: discrete adjoint gradients (via `femcore::linalg::opt::adjoint_gradient` and
//!   cached Cholesky back-substitutions) for thermal compliance `f^T T`, domain-averaged temperature,
//!   squared temperature tracking, and smooth `L_p` peak temperature.

use crate::electro::{assemble_electrostatics, prepare_electro_mesh, EPS0};
use crate::heat::assemble_heat;
use crate::problem::{decode, BcKind, Problem, Study};
use crate::Dfs;
use femcore::fem::{free_dofs, geometry, shape_gradients, RefElement};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::driver::DirectSolver;
use femcore::linalg::opt::adjoint_gradient;
use femcore::topology::{displace_simplex_mesh, ramp_interpolate, simp_interpolate};
use nyx_physics::{
    check_mesh, discretise, validate_spec, Ctx, DesignParam, Discretisation, Error, EvaluatedTopologyParam, Mesh,
    Physics, ProblemSpec, SensitivityObjective, SensitivityReport, SensitivityRequest,
};
use std::collections::BTreeMap;

/// Evaluates `Dfs::sensitivity` on `mesh`.
pub fn compute_sensitivity(
    physics: &Dfs,
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
        Study::Electrostatic => electro_sensitivity(mesh, spec, ctx, &base_problem, req),
        Study::HeatStationary => heat_stationary_sensitivity(mesh, spec, ctx, &base_problem, req),
        ref other => Err(Error::Unsupported(format!(
            "dfs study '{}' does not support adjoint sensitivities (supported: electrostatic, heat_stationary)",
            other.kind()
        ))),
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

fn collect_terminals(problem: &Problem) -> Vec<(String, Vec<i32>)> {
    let mut terminals: Vec<(String, Vec<i32>)> = Vec::new();
    for bc in &problem.bcs {
        if let BcKind::Terminal(name) = &bc.kind {
            let trimmed = name.trim();
            if let Some((_, labs)) = terminals.iter_mut().find(|(n, _)| n == trimmed) {
                for &l in &bc.labels {
                    if !labs.contains(&l) {
                        labs.push(l);
                    }
                }
            } else {
                terminals.push((trimmed.to_string(), bc.labels.clone()));
            }
        }
    }
    terminals
}

fn assemble_electro_with_cell_overrides(
    m_work: &Discretisation,
    problem: &Problem,
    eps_override: Option<&[(usize, f64)]>,
    rho_override: Option<&[(usize, f64)]>,
) -> Result<crate::electro::ElectroAssembly, String> {
    let mut asm = assemble_electrostatics(m_work, problem)?;
    if eps_override.is_none() && rho_override.is_none() {
        return Ok(asm);
    }
    // Apply per-cell topology density overrides `(eps_topo - eps_base)` and `(rho_topo - rho_base)`
    let order = problem.order;
    if let Some(ov) = eps_override {
        let base_eps = problem.materials.cell_values(m_work, "eps_r", Some(1.0))?;
        let mut delta_coeff = vec![0.0; m_work.n_cells()];
        for &(c, eps_val) in ov {
            delta_coeff[c] = EPS0 * (eps_val - base_eps[c]);
        }
        let dk = if m_work.dim == 3 && m_work.has_curved() {
            crate::electro::assemble_laplace_curved_3d(m_work, order, &delta_coeff)?
        } else {
            femcore::fem::assemble_laplace(m_work, order, &delta_coeff)
        };
        asm.k = asm.k.add_scaled(&dk, 1.0);
    }
    if let Some(ov) = rho_override {
        let base_rho = problem.materials.cell_values(m_work, "rho_v", Some(0.0))?;
        let mut delta_rho = vec![0.0; m_work.n_cells()];
        for &(c, rho_val) in ov {
            delta_rho[c] = rho_val - base_rho[c];
        }
        let df = femcore::fem::assemble_source(m_work, order, &delta_rho);
        for (fi, &dfi) in asm.f.iter_mut().zip(&df) {
            *fi += dfi;
        }
    }
    Ok(asm)
}

/// Elementwise unit-permittivity bilinear form `EPS0 * integral_{Omega_c} grad(u) . grad(v) dV`.
fn cell_grad_dots(mesh: &Discretisation, order: usize, u: &[f64], v: &[f64], cells: &[usize]) -> Vec<f64> {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    const CHUNK: usize = 512;
    let parts = femcore::par::map_chunks(cells, CHUNK, |slice| {
        let mut out = Vec::with_capacity(slice.len());
        let mut dofs = Vec::with_capacity(nloc);
        let mut grad = vec![0.0; nloc * d];
        for &c in slice {
            mesh.cell_dofs(c, order, &mut dofs);
            let mut acc = 0.0;
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                let mut gu = [0.0; 3];
                let mut gv = [0.0; 3];
                for a in 0..nloc {
                    let ua = u[dofs[a]];
                    let va = v[dofs[a]];
                    for k in 0..d {
                        let gak = grad[a * d + k];
                        gu[k] += ua * gak;
                        gv[k] += va * gak;
                    }
                }
                let mut dot = 0.0;
                for k in 0..d {
                    dot += gu[k] * gv[k];
                }
                acc += re.w[q] * dot;
            }
            out.push(geo.vol[c] * acc);
        }
        out
    });
    parts.into_iter().flatten().collect()
}

fn electro_sensitivity(
    mesh: &Mesh,
    _spec: &ProblemSpec,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    let raw_simplex = discretise(mesh)?;
    let (qp_mesh, _) = prepare_electro_mesh(&raw_simplex, base_problem);
    let m_work = qp_mesh.as_ref().map_or(&raw_simplex, |(m, _, _)| m);
    let topo_states = build_topology_states(m_work, ctx, req)?;

    let mut eps_overrides: Vec<(usize, f64)> = Vec::new();
    let mut rho_overrides: Vec<(usize, f64)> = Vec::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, property, .. } = param {
            let st = &topo_states[&idx];
            let is_rho = property.as_deref() == Some("rho_v");
            for (&c, &rb) in st.filter.design_cells().iter().zip(&st.rho_bar) {
                let (val, _) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                if is_rho {
                    rho_overrides.push((c, val));
                } else {
                    eps_overrides.push((c, val));
                }
            }
        }
    }

    let eps_ov_ref = (!eps_overrides.is_empty()).then_some(eps_overrides.as_slice());
    let rho_ov_ref = (!rho_overrides.is_empty()).then_some(rho_overrides.as_slice());
    let asm = assemble_electro_with_cell_overrides(m_work, base_problem, eps_ov_ref, rho_ov_ref)
        .map_err(Error::from_kernel)?;

    let n = asm.k.nrows();
    let free = free_dofs(n, &asm.fixed);
    let terminals = collect_terminals(base_problem);
    let mut aux_scalars = BTreeMap::new();

    // Determine forward state(s) `(u_state, v_state)`, adjoint `w_state`, scalar multiplier `chain_mult`, and `objective`.
    // Sensitivity of raw metric `M` is: `dM/dp = quad_scale * u^T (dK/dp) v + w^T (df/dp - (dK/dp) u)`.
    let (objective, chain_mult, quad_scale, u_state, v_state, w_state): (f64, f64, f64, Vec<f64>, Vec<f64>, Vec<f64>) =
        match &req.objective {
            SensitivityObjective::Capacitance { term_i, term_j, target_f } => {
                if !terminals.is_empty() {
                    if *term_i >= terminals.len() || *term_j >= terminals.len() {
                        return Err(Error::invalid(format!(
                            "capacitance terminal index ({term_i}, {term_j}) out of range for {} terminals",
                            terminals.len()
                        )));
                    }
                    let kff = asm.k.submatrix(&free, &free);
                    let solver = if free.is_empty() {
                        None
                    } else {
                        Some(DirectSolver::new(&kff, true).map_err(Error::from_kernel)?)
                    };
                    let solve_term = |ti: usize| -> Vec<f64> {
                        let mut vmap: BTreeMap<i32, f64> = asm.volt.keys().map(|&l| (l, 0.0)).collect();
                        for &l in &terminals[ti].1 {
                            vmap.insert(l, 1.0);
                        }
                        let (t_fixed, t_vals) = femcore::fem::dirichlet_constraints(m_work, &vmap, base_problem.order);
                        let mut v = vec![0.0; n];
                        for (&d, &val) in t_fixed.iter().zip(&t_vals) {
                            v[d] = val;
                        }
                        if let Some(ref s) = solver {
                            let kfx = asm.k.submatrix(&free, &t_fixed);
                            let kx = kfx.matvec(&t_vals);
                            let rhs: Vec<f64> = kx.into_iter().map(|x| -x).collect();
                            let vf = s.solve(&rhs);
                            for (a, &gi) in free.iter().enumerate() {
                                v[gi] = vf[a];
                            }
                        }
                        v
                    };
                    let vi = solve_term(*term_i);
                    let vj = if term_i == term_j { vi.clone() } else { solve_term(*term_j) };
                    let kvj = asm.k.matvec(&vj);
                    let cap_ij: f64 = vi.iter().zip(&kvj).map(|(a, b)| a * b).sum();
                    aux_scalars.insert("capacitance_f".into(), cap_ij);
                    let (obj, mult) = match target_f {
                        Some(ct) => {
                            let sc = ct.abs().max(1e-18);
                            let diff = (cap_ij - ct) / sc;
                            (0.5 * diff * diff, diff / sc)
                        }
                        None => (cap_ij, 1.0),
                    };
                    (obj, mult, 1.0, vi, vj, vec![0.0; n])
                } else {
                    // Two-conductor Potential / Ground study: C = (V^T K V) / (dV)^2
                    let mut v = vec![0.0; n];
                    for (&d, &val) in asm.fixed.iter().zip(&asm.vals) {
                        v[d] = val;
                    }
                    let mut w = vec![0.0; n];
                    let v_min = asm.volt.values().copied().fold(f64::INFINITY, f64::min);
                    let v_max = asm.volt.values().copied().fold(f64::NEG_INFINITY, f64::max);
                    let dv =
                        if (v_max - v_min).is_finite() && (v_max - v_min).abs() > 1e-15 { v_max - v_min } else { 1.0 };
                    let q_scale = 1.0 / (dv * dv);
                    if !free.is_empty() {
                        let kff = asm.k.submatrix(&free, &free);
                        let kfx = asm.k.submatrix(&free, &asm.fixed);
                        let kx = kfx.matvec(&asm.vals);
                        let rhs: Vec<f64> = free.iter().enumerate().map(|(a, &gi)| asm.f[gi] - kx[a]).collect();
                        let solver = DirectSolver::new(&kff, true).map_err(Error::from_kernel)?;
                        let vf = solver.solve(&rhs);
                        for (a, &gi) in free.iter().enumerate() {
                            v[gi] = vf[a];
                        }
                        if free.iter().any(|&gi| asm.f[gi] != 0.0) {
                            let adj_rhs: Vec<f64> = free.iter().map(|&gi| 2.0 * q_scale * asm.f[gi]).collect();
                            let wf = solver.solve(&adj_rhs);
                            for (a, &gi) in free.iter().enumerate() {
                                w[gi] = wf[a];
                            }
                        }
                    }
                    let kv = asm.k.matvec(&v);
                    let cap: f64 = q_scale * v.iter().zip(&kv).map(|(a, b)| a * b).sum::<f64>();
                    aux_scalars.insert("capacitance_f".into(), cap);
                    let (obj, mult) = match target_f {
                        Some(ct) => {
                            let sc = ct.abs().max(1e-18);
                            let diff = (cap - ct) / sc;
                            (0.5 * diff * diff, diff / sc)
                        }
                        None => (cap, 1.0),
                    };
                    (obj, mult, q_scale, v.clone(), v, w)
                }
            }
            SensitivityObjective::Compliance { target } => {
                let mut v = vec![0.0; n];
                for (&d, &val) in asm.fixed.iter().zip(&asm.vals) {
                    v[d] = val;
                }
                let has_load = free.iter().any(|&gi| asm.f[gi] != 0.0);
                if has_load && asm.vals.iter().all(|&x| x == 0.0) && !free.is_empty() {
                    // Pure charge-loaded compliance `f^T V` via `femcore::linalg::opt::adjoint_gradient`
                    let kff = asm.k.submatrix(&free, &free);
                    let ff: Vec<f64> = free.iter().map(|&gi| asm.f[gi]).collect();
                    let (comp, _, vf) = adjoint_gradient(
                        &kff,
                        &ff,
                        &ff,
                        0,
                        |_, _| Vec::new(),
                        |_| None,
                        Definiteness::PositiveDefinite,
                    )
                    .map_err(|e| Error::Solver(e.to_string()))?;
                    for (a, &gi) in free.iter().enumerate() {
                        v[gi] = vf[a];
                    }
                    aux_scalars.insert("compliance".into(), comp);
                    let (obj, mult) = match target {
                        Some(t) => (0.5 * (comp - t).powi(2), comp - t),
                        None => (comp, 1.0),
                    };
                    // d(f^T v)/dp = 2 df^T v - v^T dK v -> quad_scale = -1, w = v (Wait: explicit + df^T v + w^T df = 2 df^T v; handle below via quad_scale = 0, w = v when f_fixed == 0? Note: 2 df^T v - v^T dK v == (-1) v^T dK v + (2 v)^T (df) ? Wait: w^T (df - dK v) with w = v gives v^T df - v^T dK v, plus explicit df^T v!)
                    // Easier: for `Compliance` with `vals == 0`, `f^T v = v^T K v`, so `d(f^T v)/dp = 2 df^T v - v^T dK v = (-1.0) * v^T dK v + 2 v^T df`!
                    // Using `quad_scale = 1.0` and `w = 2 * v`: `quad_scale * v^T dK v + w^T (df - dK v) = v^T dK v + 2 v^T df - 2 v^T dK v = 2 v^T df - v^T dK v`! Exact!
                    let w: Vec<f64> = v.iter().map(|&vi| 2.0 * vi).collect();
                    (obj, mult, 1.0, v.clone(), v, w)
                } else {
                    // Stored dielectric energy `0.5 * V^T K V`
                    let mut w = vec![0.0; n];
                    if !free.is_empty() {
                        let kff = asm.k.submatrix(&free, &free);
                        let kfx = asm.k.submatrix(&free, &asm.fixed);
                        let kx = kfx.matvec(&asm.vals);
                        let rhs: Vec<f64> = free.iter().enumerate().map(|(a, &gi)| asm.f[gi] - kx[a]).collect();
                        let solver = DirectSolver::new(&kff, true).map_err(Error::from_kernel)?;
                        let vf = solver.solve(&rhs);
                        for (a, &gi) in free.iter().enumerate() {
                            v[gi] = vf[a];
                        }
                        if has_load {
                            let adj_rhs: Vec<f64> = free.iter().map(|&gi| asm.f[gi]).collect();
                            let wf = solver.solve(&adj_rhs);
                            for (a, &gi) in free.iter().enumerate() {
                                w[gi] = wf[a];
                            }
                        }
                    }
                    let kv = asm.k.matvec(&v);
                    let energy = 0.5 * v.iter().zip(&kv).map(|(a, b)| a * b).sum::<f64>();
                    aux_scalars.insert("stored_energy_j".into(), energy);
                    let (obj, mult) = match target {
                        Some(t) => (0.5 * (energy - t).powi(2), energy - t),
                        None => (energy, 1.0),
                    };
                    (obj, mult, 0.5, v.clone(), v, w)
                }
            }
            other => {
                return Err(Error::Unsupported(format!("dfs electrostatic does not support objective '{other:?}'")))
            }
        };

    // Evaluate parameter gradients
    let mut gradients = Vec::with_capacity(req.params.len());
    let mut flat_gradient = Vec::new();
    let mut param_names = Vec::with_capacity(req.params.len());
    let mut projected_densities = BTreeMap::new();
    let mut volume_fractions = BTreeMap::new();
    let mut volume_fraction_gradients = BTreeMap::new();

    let eval_asm_contraction =
        |asm_p: &crate::electro::ElectroAssembly, asm_m: &crate::electro::ElectroAssembly, inv_2h: f64| -> f64 {
            let dk = Csr::axpby(inv_2h, &asm_p.k, -inv_2h, &asm_m.k);
            let dkv = dk.matvec(&v_state);
            let dkv_u: f64 = u_state.iter().zip(&dkv).map(|(&ui, &dvi)| ui * dvi).sum();
            let mut val = quad_scale * dkv_u;
            if w_state.iter().any(|&wi| wi != 0.0) {
                let dku = dk.matvec(&u_state);
                let adj_term: f64 = w_state
                    .iter()
                    .enumerate()
                    .map(|(i, &wi)| {
                        let df_i = (asm_p.f[i] - asm_m.f[i]) * inv_2h;
                        wi * (df_i - dku[i])
                    })
                    .sum();
                val += adj_term;
            }
            chain_mult * val
        };

    let char_len = {
        let (lo, hi) = m_work.bbox();
        (0..m_work.dim).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    for (idx, param) in req.params.iter().enumerate() {
        param_names.push(param.label());
        match param {
            DesignParam::MaterialScalar { region, property, .. } => {
                let doms = ctx.resolver.domains(region)?;
                if doms.is_empty() {
                    return Err(Error::invalid(format!("region '{region}' resolved to no mesh domains")));
                }
                let p0 = base_problem.materials.get(doms[0], property).unwrap_or(if property == "eps_r" {
                    1.0
                } else {
                    0.0
                });
                let h = 1e-6 * p0.abs().max(1.0);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let asm_p = assemble_electro_with_cell_overrides(m_work, &prob_p, eps_ov_ref, rho_ov_ref)
                    .map_err(Error::from_kernel)?;
                let asm_m = assemble_electro_with_cell_overrides(m_work, &prob_m, eps_ov_ref, rho_ov_ref)
                    .map_err(Error::from_kernel)?;
                let g = eval_asm_contraction(&asm_p, &asm_m, 1.0 / (2.0 * h));
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::ParamVar { name, .. } => {
                let Some(vel) = req.shape_velocities.get(name) else {
                    return Err(Error::invalid(format!("missing precomputed shape velocity for parameter '{name}'")));
                };
                let h = vel.optimal_step(char_len, 1e-6);
                let mp = displace_simplex_mesh(m_work, h, vel).map_err(Error::from_kernel)?;
                let mm = displace_simplex_mesh(m_work, -h, vel).map_err(Error::from_kernel)?;
                let asm_p = assemble_electro_with_cell_overrides(&mp, base_problem, eps_ov_ref, rho_ov_ref)
                    .map_err(Error::from_kernel)?;
                let asm_m = assemble_electro_with_cell_overrides(&mm, base_problem, eps_ov_ref, rho_ov_ref)
                    .map_err(Error::from_kernel)?;
                let g = eval_asm_contraction(&asm_p, &asm_m, 1.0 / (2.0 * h));
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, property, .. } => {
                let st = &topo_states[&idx];
                let dcells = st.filter.design_cells();
                let is_rho = property.as_deref() == Some("rho_v");
                let dj_dbar: Vec<f64> = if is_rho {
                    let re = RefElement::new(m_work.dim, base_problem.order, 2 * base_problem.order);
                    let geo = geometry(m_work);
                    let mut dofs = Vec::new();
                    dcells
                        .iter()
                        .zip(&st.rho_bar)
                        .map(|(&c, &rb)| {
                            let (_, dval) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                            m_work.cell_dofs(c, base_problem.order, &mut dofs);
                            let mut int_w = 0.0;
                            for q in 0..re.w.len() {
                                let wq: f64 = (0..re.nloc).map(|a| re.n[q][a] * w_state[dofs[a]]).sum();
                                int_w += re.w[q] * wq;
                            }
                            chain_mult * geo.vol[c] * int_w * dval
                        })
                        .collect()
                } else {
                    let gdots_uv = cell_grad_dots(m_work, base_problem.order, &u_state, &v_state, dcells);
                    let has_w = w_state.iter().any(|&wi| wi != 0.0);
                    let gdots_wu = if has_w {
                        cell_grad_dots(m_work, base_problem.order, &w_state, &u_state, dcells)
                    } else {
                        vec![0.0; dcells.len()]
                    };
                    st.rho_bar
                        .iter()
                        .enumerate()
                        .map(|(k, &rb)| {
                            let (_, deps) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                            let d_ke = EPS0 * deps * (quad_scale * gdots_uv[k] - gdots_wu[k]);
                            chain_mult * d_ke
                        })
                        .collect()
                };
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

fn assemble_heat_with_cell_overrides(
    mesh: &Discretisation,
    problem: &Problem,
    k_override: Option<&[(usize, f64)]>,
    q_override: Option<&[(usize, f64)]>,
) -> Result<crate::heat::HeatSystem, String> {
    let mut sys = assemble_heat(mesh, problem)?;
    let order = problem.order;
    if let Some(ov) = k_override {
        let base_k: Vec<f64> = sys.ktensor.iter().map(|t| t[0][0]).collect();
        let mut dk_cell = vec![0.0; mesh.n_cells()];
        for &(c, kval) in ov {
            dk_cell[c] = kval - base_k[c];
        }
        let dk = femcore::fem::assemble_laplace(mesh, order, &dk_cell);
        sys.k = sys.k.add_scaled(&dk, 1.0);
    }
    if let Some(ov) = q_override {
        let base_q = problem.materials.cell_values(mesh, "Q", Some(0.0))?;
        let mut dq_cell = vec![0.0; mesh.n_cells()];
        for &(c, qval) in ov {
            dq_cell[c] = qval - base_q[c];
        }
        let df = femcore::fem::assemble_source(mesh, order, &dq_cell);
        for (fi, &dfi) in sys.f.iter_mut().zip(&df) {
            *fi += dfi;
        }
    }
    Ok(sys)
}

/// Evaluates `TargetTemperature` objective `(J, dJ/dT)` at fixed temperature vector `t_vec` on `mesh`.
fn eval_temperature_objective(
    mesh: &Discretisation,
    order: usize,
    t_vec: &[f64],
    cells: &[usize],
    target_k: Option<f64>,
    p_norm: Option<f64>,
) -> (f64, Vec<f64>) {
    let n = mesh.n_dofs(order);
    let re = RefElement::new(mesh.dim, order, (2 * order).max(4));
    let geo = geometry(mesh);
    let mut v_dom = 0.0;
    for &c in cells {
        v_dom += geo.vol[c];
    }
    let v_dom = v_dom.max(1e-30);
    let mut grad_t = vec![0.0; n];
    let mut dofs = Vec::new();

    if let Some(t_ref) = target_k {
        let mut int_sq = 0.0;
        for &c in cells {
            mesh.cell_dofs(c, order, &mut dofs);
            let vol = geo.vol[c];
            for q in 0..re.w.len() {
                let tq: f64 = (0..re.nloc).map(|a| re.n[q][a] * t_vec[dofs[a]]).sum();
                let diff = tq - t_ref;
                let wq = re.w[q] * vol;
                int_sq += wq * 0.5 * diff * diff;
                for a in 0..re.nloc {
                    grad_t[dofs[a]] += (wq / v_dom) * diff * re.n[q][a];
                }
            }
        }
        (int_sq / v_dom, grad_t)
    } else if let Some(p) = p_norm {
        // Shift/scale by max |T| so T^p cannot overflow for large p
        let t_scale = t_vec.iter().copied().map(f64::abs).fold(0.0f64, f64::max).max(1.0);
        let mut int_tp = 0.0;
        let mut raw_grad = vec![0.0; n];
        for &c in cells {
            mesh.cell_dofs(c, order, &mut dofs);
            let vol = geo.vol[c];
            for q in 0..re.w.len() {
                let tq: f64 = (0..re.nloc).map(|a| re.n[q][a] * t_vec[dofs[a]]).sum();
                let ts = (tq / t_scale).max(1e-30);
                let wq = re.w[q] * vol / v_dom;
                int_tp += wq * ts.powf(p);
                let d_ts = p * ts.powf(p - 1.0) / t_scale;
                for a in 0..re.nloc {
                    raw_grad[dofs[a]] += wq * d_ts * re.n[q][a];
                }
            }
        }
        let mean_tp = int_tp.max(1e-300);
        let val = t_scale * mean_tp.powf(1.0 / p);
        let factor = t_scale * (1.0 / p) * mean_tp.powf(1.0 / p - 1.0);
        for (g, &rg) in grad_t.iter_mut().zip(&raw_grad) {
            *g = factor * rg;
        }
        (val, grad_t)
    } else {
        let mut int_t = 0.0;
        for &c in cells {
            mesh.cell_dofs(c, order, &mut dofs);
            let vol = geo.vol[c];
            for q in 0..re.w.len() {
                let tq: f64 = (0..re.nloc).map(|a| re.n[q][a] * t_vec[dofs[a]]).sum();
                let wq = re.w[q] * vol / v_dom;
                int_t += wq * tq;
                for a in 0..re.nloc {
                    grad_t[dofs[a]] += wq * re.n[q][a];
                }
            }
        }
        (int_t, grad_t)
    }
}

fn heat_stationary_sensitivity(
    mesh: &Mesh,
    _spec: &ProblemSpec,
    ctx: &mut Ctx,
    base_problem: &Problem,
    req: &SensitivityRequest,
) -> Result<SensitivityReport, Error> {
    if !base_problem.k_curves.is_empty()
        || base_problem.bcs.iter().any(
            |b| matches!(b.kind, BcKind::Radiation { emissivity, view_factor, .. } if emissivity * view_factor > 0.0),
        )
    {
        return Err(Error::Unsupported(
            "stationary heat sensitivity currently requires linear conductivity and linear boundary conditions".into(),
        ));
    }

    let simplex = discretise(mesh)?;
    let order = base_problem.order;
    let topo_states = build_topology_states(&simplex, ctx, req)?;

    let mut k_overrides: Vec<(usize, f64)> = Vec::new();
    let mut q_overrides: Vec<(usize, f64)> = Vec::new();
    for (idx, param) in req.params.iter().enumerate() {
        if let DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, property, .. } = param {
            let st = &topo_states[&idx];
            let is_q = property.as_deref() == Some("Q");
            for (&c, &rb) in st.filter.design_cells().iter().zip(&st.rho_bar) {
                let (val, _) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                if is_q {
                    q_overrides.push((c, val));
                } else {
                    k_overrides.push((c, val));
                }
            }
        }
    }
    let k_ov_ref = (!k_overrides.is_empty()).then_some(k_overrides.as_slice());
    let q_ov_ref = (!q_overrides.is_empty()).then_some(q_overrides.as_slice());

    let sys =
        assemble_heat_with_cell_overrides(&simplex, base_problem, k_ov_ref, q_ov_ref).map_err(Error::from_kernel)?;
    let n = sys.k.nrows();
    let free = free_dofs(n, &sys.fixed);
    let mut t_full = vec![0.0; n];
    for (&d, &val) in sys.fixed.iter().zip(&sys.vals) {
        t_full[d] = val;
    }

    let mut aux_scalars = BTreeMap::new();
    let mut w_full = vec![0.0; n];
    let mut is_compliance_obj = false;
    let mut comp_chain_mult = 1.0;

    // Solve forward state and adjoint state `K_{ff}^T w_f = (dJ / dT)_f`
    type TempObjCells = (Vec<usize>, Option<f64>, Option<f64>);
    let (objective, obj_cells_opt): (f64, Option<TempObjCells>) = match &req.objective {
        SensitivityObjective::Compliance { target } => {
            is_compliance_obj = true;
            if sys.vals.iter().all(|&v| v == 0.0) && !free.is_empty() {
                // Direct call to `femcore::linalg::opt::adjoint_gradient`
                let kff = sys.k.submatrix(&free, &free);
                let ff: Vec<f64> = free.iter().map(|&gi| sys.f[gi]).collect();
                let (comp, _, tf) =
                    adjoint_gradient(&kff, &ff, &ff, 0, |_, _| Vec::new(), |_| None, Definiteness::PositiveDefinite)
                        .map_err(|e| Error::Solver(e.to_string()))?;
                for (a, &gi) in free.iter().enumerate() {
                    t_full[gi] = tf[a];
                    w_full[gi] = tf[a];
                }
                aux_scalars.insert("thermal_compliance".into(), comp);
                let (obj, mult) = match target {
                    Some(t) => (0.5 * (comp - t).powi(2), comp - t),
                    None => (comp, 1.0),
                };
                comp_chain_mult = mult;
                for wi in &mut w_full {
                    *wi *= mult;
                }
                (obj, None)
            } else {
                let solver_opt = if free.is_empty() {
                    None
                } else {
                    let kff = sys.k.submatrix(&free, &free);
                    let kfx = sys.k.submatrix(&free, &sys.fixed);
                    let kx = kfx.matvec(&sys.vals);
                    let rhs: Vec<f64> = free.iter().enumerate().map(|(a, &gi)| sys.f[gi] - kx[a]).collect();
                    let s = DirectSolver::new(&kff, true).map_err(Error::from_kernel)?;
                    let tf = s.solve(&rhs);
                    for (a, &gi) in free.iter().enumerate() {
                        t_full[gi] = tf[a];
                    }
                    Some(s)
                };
                let comp: f64 = sys.f.iter().zip(&t_full).map(|(fi, ti)| fi * ti).sum();
                aux_scalars.insert("thermal_compliance".into(), comp);
                let (obj, mult) = match target {
                    Some(t) => (0.5 * (comp - t).powi(2), comp - t),
                    None => (comp, 1.0),
                };
                comp_chain_mult = mult;
                if let Some(s) = solver_opt {
                    let gf: Vec<f64> = free.iter().map(|&gi| mult * sys.f[gi]).collect();
                    let wf = s.solve(&gf);
                    for (a, &gi) in free.iter().enumerate() {
                        w_full[gi] = wf[a];
                    }
                }
                (obj, None)
            }
        }
        SensitivityObjective::TargetTemperature { domain, target_k, p_norm } => {
            let cells: Vec<usize> = match domain {
                Some(dom_name) => {
                    let doms = ctx.resolver.domains(dom_name)?;
                    if doms.is_empty() {
                        return Err(Error::invalid(format!(
                            "TargetTemperature domain '{dom_name}' resolved to zero mesh domains"
                        )));
                    }
                    (0..simplex.n_cells()).filter(|&c| doms.contains(&simplex.cell_labels[c])).collect()
                }
                None => (0..simplex.n_cells()).collect(),
            };
            if cells.is_empty() {
                return Err(Error::invalid("TargetTemperature domain contains zero mesh cells"));
            }
            let solver_opt = if free.is_empty() {
                None
            } else {
                let kff = sys.k.submatrix(&free, &free);
                let kfx = sys.k.submatrix(&free, &sys.fixed);
                let kx = kfx.matvec(&sys.vals);
                let rhs: Vec<f64> = free.iter().enumerate().map(|(a, &gi)| sys.f[gi] - kx[a]).collect();
                let s = DirectSolver::new(&kff, true).map_err(Error::from_kernel)?;
                let tf = s.solve(&rhs);
                for (a, &gi) in free.iter().enumerate() {
                    t_full[gi] = tf[a];
                }
                Some(s)
            };
            let (obj, grad_t) = eval_temperature_objective(&simplex, order, &t_full, &cells, *target_k, *p_norm);
            aux_scalars.insert("temperature_objective".into(), obj);
            if let Some(s) = solver_opt {
                let gf: Vec<f64> = free.iter().map(|&gi| grad_t[gi]).collect();
                let wf = s.solve(&gf);
                for (a, &gi) in free.iter().enumerate() {
                    w_full[gi] = wf[a];
                }
            }
            (obj, Some((cells, *target_k, *p_norm)))
        }
        other => return Err(Error::Unsupported(format!("dfs heat_stationary does not support objective '{other:?}'"))),
    };

    let char_len = {
        let (lo, hi) = simplex.bbox();
        (0..simplex.dim).map(|a| hi[a] - lo[a]).fold(0.0f64, f64::max).max(1e-9)
    };

    let eval_sys_contraction = |mp_opt: Option<(&Discretisation, &Discretisation)>,
                                sys_p: &crate::heat::HeatSystem,
                                sys_m: &crate::heat::HeatSystem,
                                inv_2h: f64|
     -> f64 {
        let mut explicit_dj = 0.0;
        if is_compliance_obj {
            let df_t: f64 = (0..n).map(|i| (sys_p.f[i] - sys_m.f[i]) * inv_2h * t_full[i]).sum();
            explicit_dj += comp_chain_mult * df_t;
        } else if let (Some((mp, mm)), Some((ref cells, tk, pn))) = (mp_opt, &obj_cells_opt) {
            let (jp, _) = eval_temperature_objective(mp, order, &t_full, cells, *tk, *pn);
            let (jm, _) = eval_temperature_objective(mm, order, &t_full, cells, *tk, *pn);
            explicit_dj += (jp - jm) * inv_2h;
        }

        let dk = Csr::axpby(inv_2h, &sys_p.k, -inv_2h, &sys_m.k);
        let dkt = dk.matvec(&t_full);
        let adj_term: f64 = (0..n)
            .map(|i| {
                let wi = w_full[i];
                if wi == 0.0 {
                    0.0
                } else {
                    let df_i = (sys_p.f[i] - sys_m.f[i]) * inv_2h;
                    wi * (df_i - dkt[i])
                }
            })
            .sum();
        explicit_dj + adj_term
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
                let p0 = base_problem.materials.get(doms[0], property).unwrap_or(0.0);
                let h = 1e-6 * p0.abs().max(1.0);
                let mut prob_p = base_problem.clone();
                let mut prob_m = base_problem.clone();
                for &d in &doms {
                    prob_p.materials.domains.entry(d).or_default().insert(property.clone(), p0 + h);
                    prob_m.materials.domains.entry(d).or_default().insert(property.clone(), p0 - h);
                }
                let sys_p = assemble_heat_with_cell_overrides(&simplex, &prob_p, k_ov_ref, q_ov_ref)
                    .map_err(Error::from_kernel)?;
                let sys_m = assemble_heat_with_cell_overrides(&simplex, &prob_m, k_ov_ref, q_ov_ref)
                    .map_err(Error::from_kernel)?;
                let g = eval_sys_contraction(None, &sys_p, &sys_m, 1.0 / (2.0 * h));
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
                let sys_p = assemble_heat_with_cell_overrides(&mp, base_problem, k_ov_ref, q_ov_ref)
                    .map_err(Error::from_kernel)?;
                let sys_m = assemble_heat_with_cell_overrides(&mm, base_problem, k_ov_ref, q_ov_ref)
                    .map_err(Error::from_kernel)?;
                let g = eval_sys_contraction(Some((&mp, &mm)), &sys_p, &sys_m, 1.0 / (2.0 * h));
                gradients.push(vec![g]);
                flat_gradient.push(g);
            }
            DesignParam::TopologyDensity { simp_penalty, ramp_q, val_min, val_max, property, .. } => {
                let st = &topo_states[&idx];
                let dcells = st.filter.design_cells();
                let is_q = property.as_deref() == Some("Q");
                let dj_dbar: Vec<f64> = if is_q {
                    let re = RefElement::new(simplex.dim, order, 2 * order);
                    let geo = geometry(&simplex);
                    let mut dofs = Vec::new();
                    dcells
                        .iter()
                        .zip(&st.rho_bar)
                        .map(|(&c, &rb)| {
                            let (_, dq) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                            simplex.cell_dofs(c, order, &mut dofs);
                            let mut int_wt = 0.0;
                            for q in 0..re.w.len() {
                                let wq: f64 = (0..re.nloc).map(|a| re.n[q][a] * w_full[dofs[a]]).sum();
                                let tq: f64 = if is_compliance_obj {
                                    (0..re.nloc).map(|a| re.n[q][a] * t_full[dofs[a]]).sum()
                                } else {
                                    0.0
                                };
                                int_wt += re.w[q] * (wq + comp_chain_mult * tq);
                            }
                            geo.vol[c] * int_wt * dq
                        })
                        .collect()
                } else {
                    let gdots = cell_grad_dots(&simplex, order, &w_full, &t_full, dcells);
                    st.rho_bar
                        .iter()
                        .enumerate()
                        .map(|(k, &rb)| {
                            let (_, dk) = eval_scalar_interp(rb, *val_min, *val_max, *simp_penalty, *ramp_q);
                            -dk * gdots[k]
                        })
                        .collect()
                };
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
