// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Heat transfer on P1/P2 simplices, stationary and transient: `rho Cp dT/dt - div(k grad T) = Q`, `T = T0`, `-n.q = q0`,
//! `-n.q = h (Text - T)`. Three stages (ADR 9): `assemble_heat`, `solve_heat_states`, `heat_fields`; `solve_heat` runs them.

use crate::common::{cell_tensor_gradient, components, minmax, plane_tensors};
use crate::problem::{BcKind, Problem, Study};
use femcore::expr::eval_at;
use femcore::fem::*;
use femcore::fields::{node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::{dense_gen_eigh as gen_eigh, DirectSolver};
use femcore::material::M3;
use femcore::props::{LinearSolver, TensorKind};
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use std::collections::BTreeMap;
use std::time::Instant;

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
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
    let ktensor = plane_tensors(mesh, spec.materials.cell_tensor(mesh, TensorKind::Conductivity, None)?);
    let mut k = if spec.materials.has_tensor(TensorKind::Conductivity) {
        assemble_laplace_tensor(mesh, order, &ktensor)
    } else {
        assemble_laplace(mesh, order, &ktensor.iter().map(|t| t[0][0]).collect::<Vec<_>>())
    };
    let n = k.nrows();
    let q = spec.materials.cell_values(mesh, "Q", Some(0.0))?;
    let mut f = if q.iter().any(|&x| x != 0.0) { assemble_source(mesh, order, &q) } else { vec![0.0; n] };
    let mut tfix: BTreeMap<i32, f64> = BTreeMap::new();
    let mut has_convection = false;
    for bc in &spec.bcs {
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Temperature(t) => {
                    tfix.insert(lab, *t);
                }
                BcKind::Flux(q0) => {
                    let l = assemble_boundary_load(mesh, &[lab], *q0, order);
                    f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                }
                BcKind::Convection { h, text } => {
                    has_convection |= *h != 0.0;
                    k = k.add_scaled(&assemble_boundary_mass(mesh, &[lab], order, *h), 1.0);
                    let l = assemble_boundary_load(mesh, &[lab], h * text, order);
                    f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                }
                other => return Err(format!("'{}' is not a heat-transfer boundary condition", other.name())),
            }
        }
    }
    let (fixed, vals) = if tfix.is_empty() { (vec![], vec![]) } else { dirichlet_constraints(mesh, &tfix, order) };
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

pub fn solve_heat_states(mesh: &Discretisation, spec: &Problem, sys: &HeatSystem) -> Result<HeatStates, String> {
    let order = spec.order;
    let n = sys.k.nrows();
    let mut log = vec![];
    match &spec.study {
        Study::HeatStationary => {
            let (t, note) = solve_spd(mesh, order, 1, spec.solver, &sys.k, &sys.f, &sys.fixed, &sys.vals)?;
            log.extend(note);
            Ok(HeatStates { times: vec![], temps: vec![t], log })
        }
        Study::HeatTransient { t_end, n_out, t0: t0_expr } => {
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
            let kx = if sys.fixed.is_empty() { vec![0.0; free.len()] } else { kfx.matvec(&sys.vals) };
            let rhs0: Vec<f64> = free.iter().enumerate().map(|(a, &i)| sys.f[i] - kx[a]).collect();
            let use_amg = match spec.solver {
                LinearSolver::Direct => false,
                LinearSolver::Iterative => true,
                LinearSolver::Auto => femcore::fem::amg_wanted(mesh.dim, free.len()) && free.len() > 2500,
            } && !sys.fixed.is_empty();
            let ts: Vec<f64> = if use_amg {
                // steady part by multigrid-preconditioned CG (falls back to a factorisation if that fails)
                match solve_spd(mesh, order, 1, spec.solver, &sys.k, &sys.f, &sys.fixed, &sys.vals) {
                    Ok((u, _)) => free.iter().map(|&i| u[i]).collect(),
                    Err(_) if spec.solver == LinearSolver::Auto => DirectSolver::new(&kff, true)?.solve(&rhs0),
                    Err(e) => return Err(e),
                }
            } else {
                DirectSolver::new(&kff, true)?.solve(&rhs0)
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
            let dense_limit = 2500;
            if free.len() <= dense_limit {
                let (lam, vecs) = gen_eigh(&kff.to_dense(), &mff.to_dense(), free.len())?;
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
            } else {
                // BDF2 with 8 substeps per output interval
                let sub = 8usize;
                let dt = times[1] / sub as f64;
                let a1 = mff.scaled(1.0 / dt).add_scaled(&kff, 1.0);
                let a2 = mff.scaled(1.5 / dt).add_scaled(&kff, 1.0);
                // one multigrid hierarchy (built for the second matrix) serves every step of both systems
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
                let solve1 = |b: &[f64]| -> Result<Vec<f64>, String> {
                    match (&amg_sys, &s1) {
                        (Some(a), _) => a.solve_free(a1_full.as_ref(), b, 1e-10),
                        (None, Some(s)) => Ok(s.solve(b)),
                        _ => unreachable!(),
                    }
                };
                let solve2 = |b: &[f64]| -> Result<Vec<f64>, String> {
                    match (&amg_sys, &s2) {
                        (Some(a), _) => a.solve_free(None, b, 1e-10),
                        (None, Some(s)) => Ok(s.solve(b)),
                        _ => unreachable!(),
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
                let mut yprev = y0f.clone();
                let my0 = mff.matvec(&y0f);
                let mut y = solve1(&my0.iter().map(|v| v / dt).collect::<Vec<_>>())?;
                for (a, &i) in free.iter().enumerate() {
                    out[0][i] = ts[a] + y0f[a];
                }
                for s in 1..=sub * n_out {
                    if s > 1 {
                        let my = mff.matvec(&y);
                        let mp = mff.matvec(&yprev);
                        let rhs: Vec<f64> = my.iter().zip(&mp).map(|(a, b)| (2.0 * a - 0.5 * b) / dt).collect();
                        let ynew = solve2(&rhs)?;
                        yprev = std::mem::replace(&mut y, ynew);
                    }
                    if s % sub == 0 {
                        for (a, &i) in free.iter().enumerate() {
                            out[s / sub][i] = ts[a] + y[a];
                        }
                    }
                }
                log.push(format!(
                    "BDF2, {} steps ({} free dofs; too large for the exact dense method)",
                    sub * n_out,
                    free.len()
                ));
            }
            Ok(HeatStates { times, temps: out, log })
        }
        Study::Electrostatic => Err("heat transfer supports the stationary and transient studies".into()),
    }
}

pub fn solve_heat(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let sys = assemble_heat(mesh, spec)?;
    let st = solve_heat_states(mesh, spec, &sys)?;
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
    })
}
