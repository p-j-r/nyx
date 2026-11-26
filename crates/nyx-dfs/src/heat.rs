// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Heat transfer on P1/P2 simplices and native mixed cells (`Quad`, `Prism`, `Hex`), stationary and transient:
//! `rho Cp dT/dt - div(k grad T) = Q`, `T = T0`, `-n.q = q0`, `-n.q = h (Text - T)`, plus stationary
//! `far_field` Robin / Bayliss-Turkel and domain-applied static `infinite_element` conditions.
//! Three stages (ADR 9): `assemble_heat`, `solve_heat_states`, `heat_fields`; `solve_heat` runs them.

use crate::common::{
    assemble_boundary_load_field, assemble_boundary_load_mixed, assemble_boundary_mass_mixed,
    assemble_farfield_robin_mixed, assemble_ie_outer_robin_mixed, assemble_laplace_ie_mixed, assemble_source_field,
    assemble_source_mixed, cell_tensor_gradient, cell_tensor_mixed, cell_values_mixed, check_boundary_labels_mixed,
    classify_ie_outer_facets_mixed, components, dirichlet_constraints_field, farfield_centroid_mixed,
    mesh_has_curved_edges, minmax, plane_tensors, plane_tensors_mixed, resolve_ie_domains_mixed,
    sample_mixed_gradients_on_simplex,
};
use crate::diag::{amg_solve_free_checked, solve_spd_checked, solve_spd_mixed_checked, DfsDiag};
use crate::problem::{BcKind, Problem, Study};
use femcore::expr::eval_at;
use femcore::fem::*;
use femcore::fem_mixed::{assemble_mass_mixed, LagrangeSpace, MixedLagrangeDofMap};
use femcore::fields::{node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::{backward_residual_many, dense_gen_eigh as gen_eigh, DirectSolver};
use femcore::material::M3;
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::props::{LinearSolver, TensorKind};
use femcore::vvuq::{SourceField, UniformSource};
use nyx_physics::{
    cell_fields_of, discretise, discretise_mixed, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar,
    Solution, SolveResult,
};
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

fn check_transient_heat_guards(spec: &Problem) -> Result<(), String> {
    if matches!(spec.study, Study::HeatTransient { .. }) {
        if spec.has_infinite_element() {
            return Err(
                "infinite_element cannot be used with heat_transient: frequency-independent spatial coordinate stretching distorts volumetric heat capacity (rho * Cp * det(J_IE)) across transient diffusion timescales; use unstretched (s_u = 1) geometrically graded layers of diffusion-length thickness d = sqrt(4 * D_max * t_end * ln(1/R_0)) instead"
                    .into(),
            );
        }
        if spec.bcs.iter().any(|bc| matches!(bc.kind, BcKind::FarField { .. })) {
            return Err(
                "far_field asymptotic Laplace Robin boundary condition is only valid for stationary heat transfer (heat_stationary), not heat_transient; use unstretched geometrically graded diffusion-length layers d = sqrt(4 * D_max * t_end * ln(1/R_0))"
                    .into(),
            );
        }
    }
    Ok(())
}

fn reject_curved_simplex(mesh: &Discretisation) -> Result<(), String> {
    if mesh.has_curved() {
        return Err("unsupported: heat transfer does not support curved simplex edge nodes (R32)".into());
    }
    Ok(())
}

/// Options controlling transient BDF2 time integration in [`solve_heat_states_with_source`] and
/// [`solve_heat_bdf2_with_source`].
#[derive(Clone, Copy, Debug)]
pub struct TransientHeatOptions {
    /// Substeps per output interval (`8` by default).
    pub bdf2_substeps: usize,
    /// Force BDF2 time integration even when `free <= 2500` (bypassing dense `gen_eigh`).
    pub force_bdf2: bool,
}

impl Default for TransientHeatOptions {
    fn default() -> Self {
        Self { bdf2_substeps: 8, force_bdf2: false }
    }
}

/// Stage 1 of heat transfer: conductivity matrix, load, Dirichlet dofs and the conductivity tensors (for the flux).
pub struct HeatSystem {
    pub k: Csr,
    pub f: Vec<f64>,
    pub ktensor: Vec<M3>,
    pub fixed: Vec<usize>,
    pub vals: Vec<f64>,
    pub has_convection: bool,
}

pub fn assemble_heat(mesh: &Discretisation, spec: &Problem) -> Result<HeatSystem, String> {
    match mesh.dim {
        2 => assemble_heat_with_source::<2, UniformSource>(mesh, spec, &UniformSource),
        3 => assemble_heat_with_source::<3, UniformSource>(mesh, spec, &UniformSource),
        d => Err(format!("unsupported mesh dimension {d} in heat transfer")),
    }
}

pub fn assemble_heat_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<HeatSystem, String> {
    reject_curved_simplex(mesh)?;
    check_transient_heat_guards(spec)?;
    if !source.is_uniform() && mesh.dim != D {
        return Err(format!("SourceField<{D}> called on {}D mesh", mesh.dim));
    }
    let order = spec.order;
    check_boundary_labels(
        mesh,
        spec.bcs
            .iter()
            .filter(|b| !matches!(b.kind, BcKind::InfiniteElement(_)))
            .flat_map(|b| b.labels.iter().copied()),
    )?;
    let ktensor = plane_tensors(mesh, spec.materials.cell_tensor(mesh, TensorKind::Conductivity, None)?);
    let k_scalar: Vec<f64> = ktensor.iter().map(|t| t[0][0]).collect();

    let mut ie_rat_facets = Vec::new();
    let mut ie_exp_facets = Vec::new();
    let mut ie_map_opt = None;
    let mut mixed_opt = None;
    let mut dof_map_opt = None;

    let mut k = if spec.has_infinite_element() {
        let mixed = MixedDiscretisation::from_simplex_mesh(mesh)?;
        let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
        let dof_map = MixedLagrangeDofMap::new(&mixed, space);
        let ie_map = resolve_ie_domains_mixed(&mixed, spec)?;
        let (rf, ef) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
        let km = assemble_laplace_ie_mixed(&mixed, &dof_map, &k_scalar, Some(&ktensor), &ie_map);
        ie_rat_facets = rf;
        ie_exp_facets = ef;
        ie_map_opt = Some(ie_map);
        mixed_opt = Some(mixed);
        dof_map_opt = Some(dof_map);
        km
    } else if spec.materials.has_tensor(TensorKind::Conductivity) {
        assemble_laplace_tensor(mesh, order, &ktensor)
    } else {
        assemble_laplace(mesh, order, &k_scalar)
    };

    let n = k.nrows();
    let q = spec.materials.cell_values(mesh, "Q", Some(0.0))?;
    let mut f = if source.is_uniform() {
        if q.iter().any(|&x| x != 0.0) {
            assemble_source(mesh, order, &q)
        } else {
            vec![0.0; n]
        }
    } else {
        let extra = q.iter().any(|&x| x != 0.0).then_some(q.as_slice());
        assemble_source_field::<D, S>(mesh, order, source, 0.0, extra)?
    };
    let mut tfix: BTreeMap<i32, f64> = BTreeMap::new();
    let mut dirichlet_labels: BTreeSet<i32> = BTreeSet::new();
    let mut has_convection = false;
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            if !(1..=2).contains(ff_order) {
                return Err(format!("far_field condition 'order' must be 1 or 2, got {ff_order}"));
            }
            if !decay_power.is_finite() || *decay_power <= 0.0 {
                return Err(format!("far_field condition 'decay_power' must be positive, got {decay_power}"));
            }
            if let Some(c) = center {
                if c.iter().any(|x| !x.is_finite()) {
                    return Err("far_field condition 'center' has a NaN or infinite coordinate".into());
                }
            }
        }
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Temperature(t) => {
                    tfix.insert(lab, *t);
                    dirichlet_labels.insert(lab);
                }
                BcKind::Flux(q0) => {
                    if source.is_uniform() {
                        let l = assemble_boundary_load(mesh, &[lab], *q0, order);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    } else {
                        let l = assemble_boundary_load_field::<D, S>(mesh, &[lab], order, source, 0.0)?;
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                BcKind::Convection { h, text } => {
                    has_convection |= *h != 0.0;
                    k = k.add_scaled(&assemble_boundary_mass(mesh, &[lab], order, *h), 1.0);
                    if source.is_uniform() {
                        let l = assemble_boundary_load(mesh, &[lab], h * text, order);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    } else {
                        let l = assemble_boundary_load_field::<D, S>(mesh, &[lab], order, source, 0.0)?;
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                BcKind::FarField { .. } | BcKind::InfiniteElement(_) => {}
                other => return Err(format!("'{}' is not a heat-transfer boundary condition", other.name())),
            }
        }
    }

    let ref_labels: Vec<i32> = tfix.keys().copied().collect();
    let k_perm = FarFieldPermittivity::Tensor(&ktensor);
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !tfix.contains_key(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid(mesh, &ref_labels));
                let k_far = assemble_farfield_robin(mesh, &ff_labs, order, *ff_order, *decay_power, c, 1.0, k_perm);
                let ones = vec![1.0; n];
                let row_sums = k_far.matvec(&ones);
                if row_sums.iter().any(|&v| v > 0.0) {
                    has_convection = true;
                }
                k = k.add_scaled(&k_far, 1.0);
                let t_inf = spec.farfield_t_inf;
                if t_inf != 0.0 {
                    for (fi, &rs) in f.iter_mut().zip(&row_sums) {
                        *fi += rs * t_inf;
                    }
                }
            }
        }
    }

    if !ie_exp_facets.is_empty() {
        if let (Some(mixed), Some(dof_map), Some(ie_map)) = (&mixed_opt, &dof_map_opt, &ie_map_opt) {
            let k_ie_rob = assemble_ie_outer_robin_mixed(mixed, dof_map, &ie_exp_facets, ie_map, 1.0, k_perm);
            let ones = vec![1.0; n];
            let row_sums = k_ie_rob.matvec(&ones);
            if row_sums.iter().any(|&v| v > 0.0) {
                has_convection = true;
            }
            k = k.add_scaled(&k_ie_rob, 1.0);
            let t_inf = spec.farfield_t_inf;
            if t_inf != 0.0 {
                for (fi, &rs) in f.iter_mut().zip(&row_sums) {
                    *fi += rs * t_inf;
                }
            }
        }
    }

    let (mut fixed, mut vals) = if tfix.is_empty() {
        (vec![], vec![])
    } else if source.is_uniform() {
        dirichlet_constraints(mesh, &tfix, order)
    } else {
        dirichlet_constraints_field::<D, S>(mesh, &dirichlet_labels, order, source, 0.0)
    };

    if !ie_rat_facets.is_empty() {
        if let (Some(mixed), Some(dof_map)) = (&mixed_opt, &dof_map_opt) {
            let fixed_set: BTreeSet<usize> = fixed.iter().copied().collect();
            let mut ie_dofs = BTreeSet::new();
            let mut fbuf = Vec::new();
            for &f_idx in &ie_rat_facets {
                dof_map.facet_dofs(mixed, f_idx, &mut fbuf);
                for &d in &fbuf {
                    if !fixed_set.contains(&d) {
                        ie_dofs.insert(d);
                    }
                }
            }
            if !ie_dofs.is_empty() {
                let mut pairs: Vec<(usize, f64)> = fixed.into_iter().zip(vals).collect();
                for d in ie_dofs {
                    pairs.push((d, spec.farfield_t_inf));
                }
                pairs.sort_by_key(|p| p.0);
                fixed = pairs.iter().map(|p| p.0).collect();
                vals = pairs.iter().map(|p| p.1).collect();
            }
        }
    }

    if fixed.is_empty() && !has_convection {
        return Err("heat transfer needs a fixed-temperature or convection boundary (an all-insulated body has no unique solution)".into());
    }
    Ok(HeatSystem { k, f, ktensor, fixed, vals, has_convection })
}

/// Stage 3 of heat transfer: `T` and `q` (flat names `|q|`, `qx`, `qz` come from `Solution::flat_fields`) and the same per cell without
/// averaging.
pub fn heat_fields(mesh: &Discretisation, t: &[f64], ktensor: &[M3], order: usize) -> (Vec<Field>, CellFields) {
    let qcell = cell_tensor_gradient(mesh, t, order, ktensor, -1.0);
    let q = node_average(mesh, &qcell, mesh.dim, order);
    let qc = components(&q, mesh.dim);
    let names = femcore::section::component_names(mesh.dim);
    let fields = vec![Field::scalar("T", "K", t.to_vec()), Field::vector("q", "W/m^2", names, &qc)];
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("T", "K", cell_values_at_nodes(mesh, t, order));
    cf.push_components("q", femcore::section::component_names(mesh.dim), "W/m^2", &qcell);
    (fields, cf)
}

/// Stage 2 of heat transfer: the temperature at every dof, for the stationary study (one state) or at each output
/// time of the transient one.
pub struct HeatStates {
    /// output times, empty for the stationary study
    pub times: Vec<f64>,
    pub temps: Vec<Vec<f64>>,
    pub log: Vec<String>,
}

fn assemble_heat_forcing_and_dirichlet_at<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
    t: f64,
    dirichlet_labels: &BTreeSet<i32>,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let order = spec.order;
    let q = spec.materials.cell_values(mesh, "Q", Some(0.0))?;
    let extra = q.iter().any(|&x| x != 0.0).then_some(q.as_slice());
    let mut f = assemble_source_field::<D, S>(mesh, order, source, t, extra)?;
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Flux(_) | BcKind::Convection { .. }) {
            let l = assemble_boundary_load_field::<D, S>(mesh, &bc.labels, order, source, t)?;
            f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
        }
    }
    let (_, u_dir) = if dirichlet_labels.is_empty() {
        (vec![], vec![])
    } else {
        dirichlet_constraints_field::<D, S>(mesh, dirichlet_labels, order, source, t)
    };
    Ok((f, u_dir))
}

pub fn solve_heat_states(mesh: &Discretisation, spec: &Problem, sys: &HeatSystem) -> Result<HeatStates, String> {
    Ok(solve_heat_states_checked(mesh, spec, sys)?.0)
}

pub(crate) fn solve_heat_states_checked(
    mesh: &Discretisation,
    spec: &Problem,
    sys: &HeatSystem,
) -> Result<(HeatStates, DfsDiag), String> {
    match mesh.dim {
        2 => solve_heat_states_with_source_checked::<2, UniformSource>(mesh, spec, sys, &UniformSource, None),
        3 => solve_heat_states_with_source_checked::<3, UniformSource>(mesh, spec, sys, &UniformSource, None),
        d => Err(format!("unsupported mesh dimension {d} in heat transfer")),
    }
}

pub fn solve_heat_states_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    sys: &HeatSystem,
    source: &S,
    opts: Option<TransientHeatOptions>,
) -> Result<HeatStates, String> {
    Ok(solve_heat_states_with_source_checked::<D, S>(mesh, spec, sys, source, opts)?.0)
}

pub(crate) fn solve_heat_states_with_source_checked<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    sys: &HeatSystem,
    source: &S,
    opts: Option<TransientHeatOptions>,
) -> Result<(HeatStates, DfsDiag), String> {
    reject_curved_simplex(mesh)?;
    if !source.is_uniform() && mesh.dim != D {
        return Err(format!("SourceField<{D}> called on {}D mesh", mesh.dim));
    }
    let order = spec.order;
    let n = sys.k.nrows();
    let mut log = vec![];
    match &spec.study {
        Study::HeatStationary => {
            let tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
            let (t, note, diag) =
                solve_spd_checked(mesh, order, 1, spec.solver, &sys.k, &sys.f, &sys.fixed, &sys.vals, tol, true)?;
            log.extend(note);
            Ok((HeatStates { times: vec![], temps: vec![t], log }, diag))
        }
        Study::HeatTransient { t_end, n_out, t0: t0_expr } => {
            check_transient_heat_guards(spec)?;
            if *t_end <= 0.0 || *n_out < 1 {
                return Err("end time must be positive and at least one output time is needed".into());
            }
            let rho = spec.materials.cell_values(mesh, "rho", Some(1.0))?;
            let cp = spec.materials.cell_values(mesh, "Cp", Some(1.0))?;
            let rcp: Vec<f64> = rho.iter().zip(&cp).map(|(a, b)| a * b).collect();
            let m = assemble_mass(mesh, order, &rcp);
            let times: Vec<f64> = (0..=*n_out).map(|j| t_end * j as f64 / *n_out as f64).collect();
            let free = free_dofs(n, &sys.fixed);
            let kff = sys.k.submatrix(&free, &free);
            let mff = m.submatrix(&free, &free);
            let kfx = sys.k.submatrix(&free, &sys.fixed);
            let force_bdf2 = opts.is_some_and(|o| o.force_bdf2) || !source.is_uniform();
            let dense_limit = 2500;
            let use_amg = match spec.solver {
                LinearSolver::Direct => false,
                LinearSolver::Iterative => true,
                LinearSolver::Auto => femcore::fem::amg_wanted(mesh.dim, free.len()) && free.len() > 2500,
            } && !sys.fixed.is_empty();

            if free.len() <= dense_limit && !force_bdf2 {
                let kx = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&sys.vals) };
                let rhs0: Vec<f64> = free.iter().enumerate().map(|(a, &i)| sys.f[i] - kx[a]).collect();
                let tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
                let mut diag = DfsDiag::new("dense_eigh");
                let ts: Vec<f64> = if use_amg {
                    let (u, _note, d0) = solve_spd_checked(
                        mesh,
                        order,
                        1,
                        spec.solver,
                        &sys.k,
                        &sys.f,
                        &sys.fixed,
                        &sys.vals,
                        tol,
                        true,
                    )?;
                    diag.max_algebraic_residual = diag.max_algebraic_residual.max(d0.max_algebraic_residual);
                    diag.iterations = d0.iterations;
                    free.iter().map(|&i| u[i]).collect()
                } else {
                    let s0 = DirectSolver::new(&kff, true)?;
                    let (u, chk) = s0.solve_checked(&kff, &rhs0);
                    diag.max_algebraic_residual = diag.max_algebraic_residual.max(chk.algebraic_residual);
                    diag.rcond_estimate = chk.rcond_estimate;
                    diag.cholesky_fallbacks += chk.cholesky_fallbacks;
                    u
                };
                let coords = mesh.dof_coordinates(order);
                let mut y0 = eval_at(t0_expr, &coords, mesh.dim).map_err(|e| format!("initial temperature: {e}"))?;
                for (&i, &v) in sys.fixed.iter().zip(&sys.vals) {
                    y0[i] = v;
                }
                let y0f: Vec<f64> = free.iter().enumerate().map(|(a, &i)| y0[i] - ts[a]).collect();
                let mut out = vec![vec![0.0; n]; times.len()];
                for row in out.iter_mut() {
                    for (&i, &v) in sys.fixed.iter().zip(&sys.vals) {
                        row[i] = v;
                    }
                }
                let (lam, vecs) = gen_eigh(&kff.to_dense(), &mff.to_dense(), free.len())?;
                for (k_i, (vk, &lk)) in vecs.iter().zip(&lam).enumerate().take(16) {
                    let kv = kff.matvec(vk);
                    let mv = mff.matvec(vk);
                    let kv_n = kv.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let mv_n = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let r2: f64 = kv.iter().zip(&mv).map(|(&a, &b)| (a - lk * b).powi(2)).sum();
                    let res = r2.sqrt() / (kv_n + lk.abs() * mv_n + 1e-30);
                    if k_i == 0 {
                        diag.per_mode_residuals.push(res);
                    }
                    diag.max_algebraic_residual = diag.max_algebraic_residual.max(res);
                }
                let myf = mff.matvec(&y0f);
                let c: Vec<f64> = vecs.iter().map(|v| v.iter().zip(&myf).map(|(a, b)| a * b).sum()).collect();
                for (j, &t) in times.iter().enumerate() {
                    for (a, &i) in free.iter().enumerate() {
                        out[j][i] = ts[a];
                    }
                    for (k, v) in vecs.iter().enumerate() {
                        let w = (-lam[k] * t).exp() * c[k];
                        for (a, &i) in free.iter().enumerate() {
                            out[j][i] += w * v[a];
                        }
                    }
                }
                log.push(format!("exact semi-discrete solution (dense eigendecomposition, {} free dofs)", free.len()));
                Ok((HeatStates { times, temps: out, log }, diag))
            } else {
                let sub = opts.map_or(8, |o| o.bdf2_substeps).max(1);
                let dt = times[1] / sub as f64;
                let a1 = mff.scaled(1.0 / dt).add_scaled(&kff, 1.0);
                let a2 = mff.scaled(1.5 / dt).add_scaled(&kff, 1.0);
                let amg_sys = if use_amg {
                    let a2_full = m.scaled(1.5 / dt).add_scaled(&sys.k, 1.0);
                    femcore::fem::AmgSystem::new(mesh, order, 1, &a2_full, &sys.fixed).ok()
                } else {
                    None
                };
                let amg_sys = match amg_sys {
                    Some(a) if spec.solver == LinearSolver::Iterative || a.convergence_factor() < 0.7 => Some(a),
                    _ => None,
                };
                let a1_full = if amg_sys.is_some() { Some(m.scaled(1.0 / dt).add_scaled(&sys.k, 1.0)) } else { None };
                let (s1, s2) = if amg_sys.is_some() {
                    (None, None)
                } else {
                    (Some(DirectSolver::new(&a1, true)?), Some(DirectSolver::new(&a2, true)?))
                };
                let amg_tol = if spec.tol > 0.0 { spec.tol.min(1e-12) } else { 1e-12 };
                let mut diag = DfsDiag::new(if amg_sys.is_some() { "amg_pcg" } else { "direct_chol" });

                let run_step = |use_s1: bool,
                                rhs: &[f64],
                                mat1: Option<&Csr>,
                                mat2: Option<&Csr>,
                                s1_d: Option<&DirectSolver<f64>>,
                                s2_d: Option<&DirectSolver<f64>>,
                                d_track: Option<&mut DfsDiag>|
                 -> Result<Vec<f64>, String> {
                    if let Some(a) = &amg_sys {
                        let mat = if use_s1 { mat1 } else { mat2 };
                        let (sol, its, rel) = amg_solve_free_checked(a, mat, n, &free, rhs, amg_tol, 800)?;
                        let a_step = if use_s1 { &a1 } else { &a2 };
                        let bw_rel = if d_track.is_some() {
                            backward_residual_many(
                                a_step,
                                std::slice::from_ref(&rhs.to_vec()),
                                std::slice::from_ref(&sol),
                            )
                        } else {
                            rel
                        };
                        if let Some(dtk) = d_track {
                            dtk.record_amg(its, rel.max(bw_rel), rel <= 1e-7);
                        }
                        if rel > 1e-7 {
                            return Err(format!("the iterative solver did not converge ({its} iterations, residual {rel:.2e}); choose the direct solver"));
                        }
                        Ok(sol)
                    } else {
                        let s = if use_s1 { s1_d.unwrap() } else { s2_d.unwrap() };
                        let a_step = if use_s1 { &a1 } else { &a2 };
                        let (sol, chk) = if d_track.is_some() {
                            s.solve_checked(a_step, rhs)
                        } else {
                            (
                                s.solve(rhs),
                                femcore::linalg::driver::DirectSolveCheck {
                                    algebraic_residual: 0.0,
                                    rcond_estimate: None,
                                    stats: s.stats(),
                                    diagnostics: s.diagnostics(),
                                    cholesky_fallbacks: 0,
                                    pivot_growth: 1.0,
                                    solver_kind: "direct_chol",
                                },
                            )
                        };
                        if let Some(dtk) = d_track {
                            dtk.absorb_direct(&chk);
                        }
                        Ok(sol)
                    }
                };
                if let Some(a) = &amg_sys {
                    log.push(format!(
                        "time steps by AMG-CG: {} levels {:?}, setup {:.1} s",
                        a.n_levels(),
                        a.sizes(),
                        a.setup_seconds
                    ));
                }

                let mut out = vec![vec![0.0; n]; times.len()];
                let mut temporal_rel_est = None;
                let can_step_halve = free.len() <= 25_000 || amg_sys.is_some();

                if source.is_uniform() {
                    let kx = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&sys.vals) };
                    let rhs0: Vec<f64> = free.iter().enumerate().map(|(a, &i)| sys.f[i] - kx[a]).collect();
                    let ts: Vec<f64> = if use_amg {
                        match solve_spd_checked(
                            mesh,
                            order,
                            1,
                            spec.solver,
                            &sys.k,
                            &sys.f,
                            &sys.fixed,
                            &sys.vals,
                            amg_tol,
                            true,
                        ) {
                            Ok((u, _, _)) => free.iter().map(|&i| u[i]).collect(),
                            Err(_) if spec.solver == LinearSolver::Auto => DirectSolver::new(&kff, true)?.solve(&rhs0),
                            Err(e) => return Err(e),
                        }
                    } else {
                        DirectSolver::new(&kff, true)?.solve(&rhs0)
                    };
                    let coords = mesh.dof_coordinates(order);
                    let mut y0 =
                        eval_at(t0_expr, &coords, mesh.dim).map_err(|e| format!("initial temperature: {e}"))?;
                    for (&i, &v) in sys.fixed.iter().zip(&sys.vals) {
                        y0[i] = v;
                    }
                    let y0f: Vec<f64> = free.iter().enumerate().map(|(a, &i)| y0[i] - ts[a]).collect();
                    for row in out.iter_mut() {
                        for (&i, &v) in sys.fixed.iter().zip(&sys.vals) {
                            row[i] = v;
                        }
                    }
                    let mut yprev = y0f.clone();
                    let my0 = mff.matvec(&y0f);
                    let mut y = run_step(
                        true,
                        &my0.iter().map(|v| v / dt).collect::<Vec<_>>(),
                        a1_full.as_ref(),
                        None,
                        s1.as_ref(),
                        s2.as_ref(),
                        Some(&mut diag),
                    )?;
                    for (a, &i) in free.iter().enumerate() {
                        out[0][i] = ts[a] + y0f[a];
                    }
                    for s in 1..=sub * n_out {
                        if s > 1 {
                            let my = mff.matvec(&y);
                            let mp = mff.matvec(&yprev);
                            let rhs: Vec<f64> = my.iter().zip(&mp).map(|(a, b)| (2.0 * a - 0.5 * b) / dt).collect();
                            let ynew = run_step(
                                false,
                                &rhs,
                                a1_full.as_ref(),
                                None,
                                s1.as_ref(),
                                s2.as_ref(),
                                Some(&mut diag),
                            )?;
                            yprev = std::mem::replace(&mut y, ynew);
                        }
                        if s % sub == 0 {
                            for (a, &i) in free.iter().enumerate() {
                                out[s / sub][i] = ts[a] + y[a];
                            }
                        }
                    }
                    if can_step_halve {
                        let sub_h = 2 * sub;
                        let dt_h = dt * 0.5;
                        let a1_h_full = amg_sys.as_ref().map(|_| m.scaled(1.0 / dt_h).add_scaled(&sys.k, 1.0));
                        let a2_h_full = amg_sys.as_ref().map(|_| m.scaled(1.5 / dt_h).add_scaled(&sys.k, 1.0));
                        let (s1_h, s2_h) = if amg_sys.is_some() {
                            (None, None)
                        } else {
                            let a1h = mff.scaled(1.0 / dt_h).add_scaled(&kff, 1.0);
                            let a2h = mff.scaled(1.5 / dt_h).add_scaled(&kff, 1.0);
                            (Some(DirectSolver::new(&a1h, true)?), Some(DirectSolver::new(&a2h, true)?))
                        };
                        let mut yp_h = y0f.clone();
                        let mut y_h = if let Some(a) = &amg_sys {
                            amg_solve_free_checked(
                                a,
                                a1_h_full.as_ref(),
                                n,
                                &free,
                                &my0.iter().map(|v| v / dt_h).collect::<Vec<_>>(),
                                amg_tol,
                                800,
                            )?
                            .0
                        } else {
                            s1_h.as_ref().unwrap().solve(&my0.iter().map(|v| v / dt_h).collect::<Vec<_>>())
                        };
                        for _step in 2..=sub_h * n_out {
                            let my = mff.matvec(&y_h);
                            let mp = mff.matvec(&yp_h);
                            let rhs: Vec<f64> = my.iter().zip(&mp).map(|(a, b)| (2.0 * a - 0.5 * b) / dt_h).collect();
                            let ynew = if let Some(a) = &amg_sys {
                                amg_solve_free_checked(a, a2_h_full.as_ref(), n, &free, &rhs, amg_tol, 800)?.0
                            } else {
                                s2_h.as_ref().unwrap().solve(&rhs)
                            };
                            yp_h = std::mem::replace(&mut y_h, ynew);
                        }
                        let diff_sq: f64 = y_h.iter().zip(&y).map(|(a, b)| (a - b).powi(2)).sum();
                        let norm_h = (0..free.len()).map(|a| (ts[a] + y_h[a]).powi(2)).sum::<f64>().sqrt();
                        let _ = out.last();
                        temporal_rel_est = Some(diff_sq.sqrt() / (3.0 * norm_h.max(1e-30)));
                    }
                } else {
                    let mut dirichlet_labels = BTreeSet::new();
                    for bc in &spec.bcs {
                        if matches!(bc.kind, BcKind::Temperature(_)) {
                            dirichlet_labels.extend(bc.labels.iter().copied());
                        }
                    }
                    let mfx = m.submatrix(&free, &sys.fixed);
                    let coords = mesh.dof_coordinates(order);
                    let mut t0_vec =
                        eval_at(t0_expr, &coords, mesh.dim).map_err(|e| format!("initial temperature: {e}"))?;
                    let (_, mut ud_prev) = if dirichlet_labels.is_empty() {
                        (vec![], vec![])
                    } else {
                        dirichlet_constraints_field::<D, S>(mesh, &dirichlet_labels, order, source, 0.0)
                    };
                    for (&i, &v) in sys.fixed.iter().zip(&ud_prev) {
                        t0_vec[i] = v;
                    }
                    out[0].copy_from_slice(&t0_vec);
                    let mut t_prev: Vec<f64> = free.iter().map(|&i| t0_vec[i]).collect();

                    let (f1, mut ud_cur) =
                        assemble_heat_forcing_and_dirichlet_at::<D, S>(mesh, spec, source, dt, &dirichlet_labels)?;
                    let kx1 = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&ud_cur) };
                    let ud_diff: Vec<f64> = ud_cur.iter().zip(&ud_prev).map(|(a, b)| (a - b) / dt).collect();
                    let mxd1 = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { mfx.matvec(&ud_diff) };
                    let mt0 = mff.matvec(&t_prev);
                    let rhs1: Vec<f64> =
                        free.iter().enumerate().map(|(a, &i)| mt0[a] / dt + f1[i] - kx1[a] - mxd1[a]).collect();
                    let mut t_cur =
                        run_step(true, &rhs1, a1_full.as_ref(), None, s1.as_ref(), s2.as_ref(), Some(&mut diag))?;
                    if 1 % sub == 0 {
                        for (&i, &v) in sys.fixed.iter().zip(&ud_cur) {
                            out[1 / sub][i] = v;
                        }
                        for (a, &i) in free.iter().enumerate() {
                            out[1 / sub][i] = t_cur[a];
                        }
                    }

                    for s in 2..=sub * n_out {
                        let ts_eval = s as f64 * dt;
                        let (fs, ud_s) = assemble_heat_forcing_and_dirichlet_at::<D, S>(
                            mesh,
                            spec,
                            source,
                            ts_eval,
                            &dirichlet_labels,
                        )?;
                        let kxs = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&ud_s) };
                        let ud_bdf: Vec<f64> = (0..sys.fixed.len())
                            .map(|k| (1.5 * ud_s[k] - 2.0 * ud_cur[k] + 0.5 * ud_prev[k]) / dt)
                            .collect();
                        let mxs = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { mfx.matvec(&ud_bdf) };
                        let mtc = mff.matvec(&t_cur);
                        let mtp = mff.matvec(&t_prev);
                        let rhs_s: Vec<f64> = free
                            .iter()
                            .enumerate()
                            .map(|(a, &i)| (2.0 * mtc[a] - 0.5 * mtp[a]) / dt + fs[i] - kxs[a] - mxs[a])
                            .collect();
                        let t_new =
                            run_step(false, &rhs_s, a1_full.as_ref(), None, s1.as_ref(), s2.as_ref(), Some(&mut diag))?;
                        t_prev = std::mem::replace(&mut t_cur, t_new);
                        ud_prev = std::mem::replace(&mut ud_cur, ud_s);
                        if s % sub == 0 {
                            let idx = s / sub;
                            for (&i, &v) in sys.fixed.iter().zip(&ud_cur) {
                                out[idx][i] = v;
                            }
                            for (a, &i) in free.iter().enumerate() {
                                out[idx][i] = t_cur[a];
                            }
                        }
                    }

                    if can_step_halve {
                        let sub_h = 2 * sub;
                        let dt_h = dt * 0.5;
                        let a1_h_full = amg_sys.as_ref().map(|_| m.scaled(1.0 / dt_h).add_scaled(&sys.k, 1.0));
                        let a2_h_full = amg_sys.as_ref().map(|_| m.scaled(1.5 / dt_h).add_scaled(&sys.k, 1.0));
                        let (s1_h, s2_h) = if amg_sys.is_some() {
                            (None, None)
                        } else {
                            let a1h = mff.scaled(1.0 / dt_h).add_scaled(&kff, 1.0);
                            let a2h = mff.scaled(1.5 / dt_h).add_scaled(&kff, 1.0);
                            (Some(DirectSolver::new(&a1h, true)?), Some(DirectSolver::new(&a2h, true)?))
                        };
                        let mut ud_p_h = if dirichlet_labels.is_empty() {
                            vec![]
                        } else {
                            dirichlet_constraints_field::<D, S>(mesh, &dirichlet_labels, order, source, 0.0).1
                        };
                        let mut tp_h: Vec<f64> = free.iter().map(|&i| t0_vec[i]).collect();
                        let (f1h, mut ud_c_h) = assemble_heat_forcing_and_dirichlet_at::<D, S>(
                            mesh,
                            spec,
                            source,
                            dt_h,
                            &dirichlet_labels,
                        )?;
                        let kx1h = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&ud_c_h) };
                        let ud_dh: Vec<f64> = ud_c_h.iter().zip(&ud_p_h).map(|(a, b)| (a - b) / dt_h).collect();
                        let mxd1h = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { mfx.matvec(&ud_dh) };
                        let rhs1h: Vec<f64> = free
                            .iter()
                            .enumerate()
                            .map(|(a, &i)| mt0[a] / dt_h + f1h[i] - kx1h[a] - mxd1h[a])
                            .collect();
                        let mut tc_h = if let Some(a) = &amg_sys {
                            amg_solve_free_checked(a, a1_h_full.as_ref(), n, &free, &rhs1h, amg_tol, 800)?.0
                        } else {
                            s1_h.as_ref().unwrap().solve(&rhs1h)
                        };
                        for s in 2..=sub_h * n_out {
                            let ts_eval = s as f64 * dt_h;
                            let (fs, ud_s) = assemble_heat_forcing_and_dirichlet_at::<D, S>(
                                mesh,
                                spec,
                                source,
                                ts_eval,
                                &dirichlet_labels,
                            )?;
                            let kxs = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&ud_s) };
                            let ud_bdf: Vec<f64> = (0..sys.fixed.len())
                                .map(|k| (1.5 * ud_s[k] - 2.0 * ud_c_h[k] + 0.5 * ud_p_h[k]) / dt_h)
                                .collect();
                            let mxs = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { mfx.matvec(&ud_bdf) };
                            let mtc = mff.matvec(&tc_h);
                            let mtp = mff.matvec(&tp_h);
                            let rhs_s: Vec<f64> = free
                                .iter()
                                .enumerate()
                                .map(|(a, &i)| (2.0 * mtc[a] - 0.5 * mtp[a]) / dt_h + fs[i] - kxs[a] - mxs[a])
                                .collect();
                            let t_new = if let Some(a) = &amg_sys {
                                amg_solve_free_checked(a, a2_h_full.as_ref(), n, &free, &rhs_s, amg_tol, 800)?.0
                            } else {
                                s2_h.as_ref().unwrap().solve(&rhs_s)
                            };
                            tp_h = std::mem::replace(&mut tc_h, t_new);
                            ud_p_h = std::mem::replace(&mut ud_c_h, ud_s);
                        }
                        let diff_sq: f64 = tc_h.iter().zip(&t_cur).map(|(a, b)| (a - b).powi(2)).sum();
                        let norm_h = tc_h.iter().map(|x| x * x).sum::<f64>().sqrt();
                        temporal_rel_est = Some(diff_sq.sqrt() / (3.0 * norm_h.max(1e-30)));
                    }
                }
                if let Some(t_rel) = temporal_rel_est {
                    diag.set_temporal_rel(t_rel);
                }
                log.push(format!(
                    "BDF2, {} steps ({} free dofs; {})",
                    sub * n_out,
                    free.len(),
                    if force_bdf2 && free.len() <= dense_limit {
                        "opt-in BDF2 mode"
                    } else {
                        "too large for the exact dense method"
                    }
                ));
                if amg_sys.is_some() {
                    if let Some(its) = diag.iterations {
                        log.push(format!(
                            "AMG-CG BDF2 total iterations = {its}, max residual = {:.2e}",
                            diag.max_algebraic_residual
                        ));
                    }
                }
                Ok((HeatStates { times, temps: out, log }, diag))
            }
        }
        Study::Electrostatic | Study::Magnetostatic { .. } => {
            Err("heat transfer supports the stationary and transient studies".into())
        }
    }
}

pub fn solve_heat_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<SolveResult, String> {
    solve_heat_bdf2_with_source::<D, S>(mesh, spec, source, TransientHeatOptions::default())
}

pub fn solve_heat_bdf2_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
    opts: TransientHeatOptions,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let sys = assemble_heat_with_source::<D, S>(mesh, spec, source)?;
    let use_opts = if opts.force_bdf2 || !source.is_uniform() { Some(opts) } else { None };
    let (st, diag) = solve_heat_states_with_source_checked::<D, S>(mesh, spec, &sys, source, use_opts)?;
    let (diagnostics, error_budget) = diag.finish();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut solutions = vec![];
    for (j, row) in st.temps.iter().enumerate() {
        let (a, b) = minmax(row);
        lo = lo.min(a);
        hi = hi.max(b);
        let (fields, cell_fields) = heat_fields(mesh, row, &sys.ktensor, order);
        let (label, value) = match st.times.get(j) {
            Some(&t) => (format!("t = {:.6} s", t).replace(".000000", ""), Some(t)),
            None => ("stationary".to_string(), None),
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value,
            ..Default::default()
        });
    }
    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: mesh.dim as u8,
        solutions,
        scalars: vec![("T min (K)".into(), Scalar::Num(lo)), ("T max (K)".into(), Scalar::Num(hi))],
        log: st.log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: sys.k.nrows(),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

pub fn solve_heat(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let sys = assemble_heat(mesh, spec)?;
    let (st, diag) = solve_heat_states_checked(mesh, spec, &sys)?;
    let (diagnostics, error_budget) = diag.finish();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut solutions = vec![];
    for (j, row) in st.temps.iter().enumerate() {
        let (a, b) = minmax(row);
        lo = lo.min(a);
        hi = hi.max(b);
        let (fields, cell_fields) = heat_fields(mesh, row, &sys.ktensor, order);
        let (label, value) = match st.times.get(j) {
            Some(&t) => (format!("t = {:.6} s", t).replace(".000000", ""), Some(t)),
            None => ("stationary".to_string(), None),
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value,
            ..Default::default()
        });
    }
    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: mesh.dim as u8,
        solutions,
        scalars: vec![("T min (K)".into(), Scalar::Num(lo)), ("T max (K)".into(), Scalar::Num(hi))],
        log: st.log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: sys.k.nrows(),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

fn heat_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    t_mixed: &[f64],
    ktensor: &[M3],
    order: usize,
) -> (Vec<Field>, CellFields, Vec<f64>) {
    let t_simplex = proj.project_lagrange_blocks(mixed, dof_map, t_mixed, 1, order);
    let scale = vec![-1.0; mixed.n_cells()];
    let qcell = sample_mixed_gradients_on_simplex(mixed, dof_map, proj, t_mixed, Some(ktensor), &scale, order);
    let q = node_average(simplex, &qcell, simplex.dim, order);
    let qc = components(&q, simplex.dim);
    let names = femcore::section::component_names(simplex.dim);
    let fields = vec![Field::scalar("T", "K", t_simplex.clone()), Field::vector("q", "W/m^2", names, &qc)];
    let mut cf = CellFields::new(simplex, order);
    cf.push_real("T", "K", cell_values_at_nodes(simplex, &t_simplex, order));
    cf.push_components("q", names, "W/m^2", &qcell);
    (fields, cf, t_simplex)
}

/// Solves stationary or transient heat transfer natively on a mixed mesh (`Quad`, `Prism`, `Hex`)
/// and/or with domain-applied static Infinite Elements (`BcKind::InfiniteElement`), preserving the
/// pure-simplex fast path when `mesh.is_pure_simplex() && !spec.has_infinite_element()`.
pub fn solve_heat_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String> {
    check_transient_heat_guards(spec)?;
    if mesh.is_pure_simplex() {
        if mesh_has_curved_edges(mesh) {
            return Err("unsupported: heat transfer does not support curved simplex edge nodes (R32)".into());
        }
        if !spec.has_infinite_element() {
            let d = discretise(mesh).map_err(|e| e.to_string())?;
            return solve_heat(&d, spec);
        }
    }

    let t0 = Instant::now();
    let order = spec.order;
    let mixed = discretise_mixed(mesh).map_err(|e| e.to_string())?;
    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, &mixed).map_err(|e| e.to_string())?;
    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(&mixed, space);

    check_boundary_labels_mixed(
        &mixed,
        spec.bcs
            .iter()
            .filter(|b| !matches!(b.kind, BcKind::InfiniteElement(_)))
            .flat_map(|b| b.labels.iter().copied()),
    )?;
    let ie_map = resolve_ie_domains_mixed(&mixed, spec)?;
    let ktensor =
        plane_tensors_mixed(&mixed, cell_tensor_mixed(&spec.materials, &mixed, TensorKind::Conductivity, None)?);
    let k_scalar: Vec<f64> = ktensor.iter().map(|t| t[0][0]).collect();
    let has_aniso_k = spec.materials.has_tensor(TensorKind::Conductivity);

    let mut k =
        assemble_laplace_ie_mixed(&mixed, &dof_map, &k_scalar, has_aniso_k.then_some(ktensor.as_slice()), &ie_map);
    let n = dof_map.n_dofs;
    let q = cell_values_mixed(&spec.materials, &mixed, "Q", Some(0.0))?;
    let mut f = if q.iter().any(|&x| x != 0.0) { assemble_source_mixed(&mixed, &dof_map, &q) } else { vec![0.0; n] };

    let mut tfix: BTreeMap<i32, f64> = BTreeMap::new();
    let mut has_convection = false;
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            if !(1..=2).contains(ff_order) {
                return Err(format!("far_field condition 'order' must be 1 or 2, got {ff_order}"));
            }
            if !decay_power.is_finite() || *decay_power <= 0.0 {
                return Err(format!("far_field condition 'decay_power' must be positive, got {decay_power}"));
            }
            if let Some(c) = center {
                if c.iter().any(|x| !x.is_finite()) {
                    return Err("far_field condition 'center' has a NaN or infinite coordinate".into());
                }
            }
        }
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Temperature(t) => {
                    tfix.insert(lab, *t);
                }
                BcKind::Flux(q0) => {
                    let l = assemble_boundary_load_mixed(&mixed, &dof_map, &[lab], *q0);
                    f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                }
                BcKind::Convection { h, text } => {
                    has_convection |= *h != 0.0;
                    k = k.add_scaled(&assemble_boundary_mass_mixed(&mixed, &dof_map, &[lab], *h), 1.0);
                    let l = assemble_boundary_load_mixed(&mixed, &dof_map, &[lab], h * text);
                    f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                }
                BcKind::FarField { .. } | BcKind::InfiniteElement(_) => {}
                other => return Err(format!("'{}' is not a heat-transfer boundary condition", other.name())),
            }
        }
    }

    let ref_labels: Vec<i32> = tfix.keys().copied().collect();
    let k_perm = FarFieldPermittivity::Tensor(&ktensor);
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !tfix.contains_key(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid_mixed(&mixed, &ref_labels));
                let k_far =
                    assemble_farfield_robin_mixed(&mixed, &dof_map, &ff_labs, *ff_order, *decay_power, c, 1.0, k_perm);
                let ones = vec![1.0; n];
                let row_sums = k_far.matvec(&ones);
                if row_sums.iter().any(|&v| v > 0.0) {
                    has_convection = true;
                }
                k = k.add_scaled(&k_far, 1.0);
                let t_inf = spec.farfield_t_inf;
                if t_inf != 0.0 {
                    for (fi, &rs) in f.iter_mut().zip(&row_sums) {
                        *fi += rs * t_inf;
                    }
                }
            }
        }
    }

    let (ie_rat_facets, ie_exp_facets) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
    if !ie_exp_facets.is_empty() {
        let k_ie_rob = assemble_ie_outer_robin_mixed(&mixed, &dof_map, &ie_exp_facets, &ie_map, 1.0, k_perm);
        let ones = vec![1.0; n];
        let row_sums = k_ie_rob.matvec(&ones);
        if row_sums.iter().any(|&v| v > 0.0) {
            has_convection = true;
        }
        k = k.add_scaled(&k_ie_rob, 1.0);
        let t_inf = spec.farfield_t_inf;
        if t_inf != 0.0 {
            for (fi, &rs) in f.iter_mut().zip(&row_sums) {
                *fi += rs * t_inf;
            }
        }
    }

    let (mut fixed, mut vals) = dof_map.dirichlet_constraints(&mixed, &tfix);
    if !ie_rat_facets.is_empty() {
        let fixed_set: BTreeSet<usize> = fixed.iter().copied().collect();
        let mut ie_dofs = BTreeSet::new();
        let mut fbuf = Vec::new();
        for &f_idx in &ie_rat_facets {
            dof_map.facet_dofs(&mixed, f_idx, &mut fbuf);
            for &d in &fbuf {
                if !fixed_set.contains(&d) {
                    ie_dofs.insert(d);
                }
            }
        }
        if !ie_dofs.is_empty() {
            let mut pairs: Vec<(usize, f64)> = fixed.into_iter().zip(vals).collect();
            for d in ie_dofs {
                pairs.push((d, spec.farfield_t_inf));
            }
            pairs.sort_by_key(|p| p.0);
            fixed = pairs.iter().map(|p| p.0).collect();
            vals = pairs.iter().map(|p| p.1).collect();
        }
    }

    if fixed.is_empty() && !has_convection {
        return Err("heat transfer needs a fixed-temperature or convection boundary (an all-insulated body has no unique solution)".into());
    }

    let mut log = Vec::new();
    let mut diag = DfsDiag::new("direct_chol");
    let (times, temps_mixed) = match &spec.study {
        Study::HeatStationary => {
            let (t_sol, note, d0) = solve_spd_mixed_checked(spec.solver, &k, &f, &fixed, &vals)?;
            diag = d0;
            log.extend(note);
            (vec![], vec![t_sol])
        }
        Study::HeatTransient { t_end, n_out, t0: t0_expr } => {
            if *t_end <= 0.0 || *n_out < 1 {
                return Err("end time must be positive and at least one output time is needed".into());
            }
            let rho = cell_values_mixed(&spec.materials, &mixed, "rho", Some(1.0))?;
            let cp = cell_values_mixed(&spec.materials, &mixed, "Cp", Some(1.0))?;
            let rcp: Vec<f64> = rho.iter().zip(&cp).map(|(a, b)| a * b).collect();
            let m = assemble_mass_mixed(&mixed, &dof_map, &rcp);
            let times: Vec<f64> = (0..=*n_out).map(|j| t_end * j as f64 / *n_out as f64).collect();
            let free = free_dofs(n, &fixed);
            let kff = k.submatrix(&free, &free);
            let mff = m.submatrix(&free, &free);
            let kfx = k.submatrix(&free, &fixed);
            let kx = if fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&vals) };
            let rhs0: Vec<f64> = free.iter().enumerate().map(|(a, &i)| f[i] - kx[a]).collect();
            let s0 = DirectSolver::new(&kff, true)?;
            let (ts, chk0) = s0.solve_checked(&kff, &rhs0);
            diag.absorb_direct(&chk0);
            let coords = dof_map.dof_coordinates(&mixed);
            let mut y0 = eval_at(t0_expr, &coords, mixed.dim).map_err(|e| format!("initial temperature: {e}"))?;
            for (&i, &v) in fixed.iter().zip(&vals) {
                y0[i] = v;
            }
            let y0f: Vec<f64> = free.iter().enumerate().map(|(a, &i)| y0[i] - ts[a]).collect();
            let mut out = vec![vec![0.0; n]; times.len()];
            for row in out.iter_mut() {
                for (&i, &v) in fixed.iter().zip(&vals) {
                    row[i] = v;
                }
            }
            let dense_limit = 2500;
            if free.len() <= dense_limit {
                diag.solver_kind = "dense_eigh".into();
                let (lam, vecs) = gen_eigh(&kff.to_dense(), &mff.to_dense(), free.len())?;
                for (k_i, (vk, &lk)) in vecs.iter().zip(&lam).enumerate().take(16) {
                    let kv = kff.matvec(vk);
                    let mv = mff.matvec(vk);
                    let kv_n = kv.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let mv_n = mv.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let r2: f64 = kv.iter().zip(&mv).map(|(&a, &b)| (a - lk * b).powi(2)).sum();
                    let res = r2.sqrt() / (kv_n + lk.abs() * mv_n + 1e-30);
                    if k_i == 0 {
                        diag.per_mode_residuals.push(res);
                    }
                    diag.max_algebraic_residual = diag.max_algebraic_residual.max(res);
                }
                let myf = mff.matvec(&y0f);
                let c: Vec<f64> = vecs.iter().map(|v| v.iter().zip(&myf).map(|(a, b)| a * b).sum()).collect();
                for (j, &t) in times.iter().enumerate() {
                    for (a, &i) in free.iter().enumerate() {
                        out[j][i] = ts[a];
                    }
                    for (idx, v) in vecs.iter().enumerate() {
                        let w = (-lam[idx] * t).exp() * c[idx];
                        for (a, &i) in free.iter().enumerate() {
                            out[j][i] += w * v[a];
                        }
                    }
                }
                log.push(format!(
                    "exact semi-discrete solution (dense eigendecomposition, {} free mixed-cell dofs)",
                    free.len()
                ));
            } else {
                let sub = 8usize;
                let dt = times[1] / sub as f64;
                let a1 = mff.scaled(1.0 / dt).add_scaled(&kff, 1.0);
                let a2 = mff.scaled(1.5 / dt).add_scaled(&kff, 1.0);
                let s1 = DirectSolver::new(&a1, true)?;
                let s2 = DirectSolver::new(&a2, true)?;
                let mut yprev = y0f.clone();
                let my0 = mff.matvec(&y0f);
                let rhs1: Vec<f64> = my0.iter().map(|v| v / dt).collect();
                let (mut y, chk1) = s1.solve_checked(&a1, &rhs1);
                diag.absorb_direct(&chk1);
                for (a, &i) in free.iter().enumerate() {
                    out[0][i] = ts[a] + y0f[a];
                }
                for s in 1..=sub * n_out {
                    if s > 1 {
                        let my = mff.matvec(&y);
                        let mp = mff.matvec(&yprev);
                        let rhs: Vec<f64> = my.iter().zip(&mp).map(|(a, b)| (2.0 * a - 0.5 * b) / dt).collect();
                        let (ynew, chks) = s2.solve_checked(&a2, &rhs);
                        diag.absorb_direct(&chks);
                        yprev = std::mem::replace(&mut y, ynew);
                    }
                    if s % sub == 0 {
                        for (a, &i) in free.iter().enumerate() {
                            out[s / sub][i] = ts[a] + y[a];
                        }
                    }
                }
                if free.len() <= 25_000 {
                    let sub_h = 2 * sub;
                    let dt_h = dt * 0.5;
                    let a1h = mff.scaled(1.0 / dt_h).add_scaled(&kff, 1.0);
                    let a2h = mff.scaled(1.5 / dt_h).add_scaled(&kff, 1.0);
                    let s1h = DirectSolver::new(&a1h, true)?;
                    let s2h = DirectSolver::new(&a2h, true)?;
                    let mut yp_h = y0f.clone();
                    let mut y_h = s1h.solve(&my0.iter().map(|v| v / dt_h).collect::<Vec<_>>());
                    for _ in 2..=sub_h * n_out {
                        let my = mff.matvec(&y_h);
                        let mp = mff.matvec(&yp_h);
                        let rhs: Vec<f64> = my.iter().zip(&mp).map(|(a, b)| (2.0 * a - 0.5 * b) / dt_h).collect();
                        let ynew = s2h.solve(&rhs);
                        yp_h = std::mem::replace(&mut y_h, ynew);
                    }
                    let diff_sq: f64 = y_h.iter().zip(&y).map(|(a, b)| (a - b).powi(2)).sum();
                    let norm_h = (0..free.len()).map(|a| (ts[a] + y_h[a]).powi(2)).sum::<f64>().sqrt();
                    diag.set_temporal_rel(diff_sq.sqrt() / (3.0 * norm_h.max(1e-30)));
                }
                log.push(format!("BDF2, {} steps ({} free mixed-cell dofs)", sub * n_out, free.len()));
            }
            (times, out)
        }
        Study::Electrostatic | Study::Magnetostatic { .. } => {
            return Err("heat transfer supports the stationary and transient studies".into());
        }
    };
    let (diagnostics, error_budget) = diag.finish();

    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut solutions = Vec::with_capacity(temps_mixed.len());
    for (j, row) in temps_mixed.iter().enumerate() {
        let (fields, cell_fields, t_simplex) =
            heat_fields_mixed(&mixed, &dof_map, &simplex, &proj, row, &ktensor, order);
        let (a, b) = minmax(&t_simplex);
        lo = lo.min(a);
        hi = hi.max(b);
        let (label, value) = match times.get(j) {
            Some(&t) => (format!("t = {:.6} s", t).replace(".000000", ""), Some(t)),
            None => ("stationary".to_string(), None),
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value,
            ..Default::default()
        });
    }

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: mixed.dim as u8,
        solutions,
        scalars: vec![("T min (K)".into(), Scalar::Num(lo)), ("T max (K)".into(), Scalar::Num(hi))],
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}
