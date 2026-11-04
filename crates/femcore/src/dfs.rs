//! Scalar diffusion physics on P1/P2 simplices: electrostatics and heat transfer
//! (stationary and transient):
//!
//!   electrostatics   -div(eps0 eps_r grad V) = rho_v,       V = V0,  n.D = -rho_s
//!   heat             rho Cp dT/dt - div(k grad T) = Q,      T = T0,  -n.q = q0,  -n.q = h (Text - T)

use crate::expr::eval_at;
use crate::fem::*;
use crate::fields::{gradients_at_nodes, node_average, CellFields};
use crate::solvers::driver::{dense_gen_eigh as gen_eigh, DirectSolver};
use crate::mesh::SimplexMesh;
use crate::problem::*;
use crate::material::M3;
use crate::solvers::csr::RCsr as Csr;
use std::collections::BTreeMap;
use std::time::Instant;

pub const EPS0: f64 = 8.8541878128e-12;

/// The wafer tensors of the cells in the axes of the mesh (on a cross-section the x-z block comes first, see `section::plane_tensor`).
fn plane_tensors(mesh: &SimplexMesh, t: Vec<M3>) -> Vec<M3> {
    t.iter().map(|m| crate::section::plane_tensor(m, mesh.dim)).collect()
}

pub(crate) fn check_labels(mesh: &SimplexMesh, bcs: &[Bc]) -> Result<(), String> {
    let known = mesh.boundary_labels();
    let mut bad: Vec<i32> = bcs.iter().flat_map(|b| b.labels.iter().copied()).filter(|l| !known.contains(l)).collect();
    bad.sort();
    bad.dedup();
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!("boundary number(s) {bad:?} do not exist in this mesh (it has {known:?})"))
    }
}

/// `scale[c] * grad u` at the local dofs of every cell, not averaged: `g[(c*nloc + a)*dim + k]`.
pub(crate) fn cell_gradient(mesh: &SimplexMesh, u: &[f64], order: usize, scale: &[f64]) -> Vec<f64> {
    let d = mesh.dim;
    let nloc = mesh.nloc(order);
    let mut g = gradients_at_nodes(mesh, u, 1, order);
    for c in 0..mesh.n_cells() {
        for a in 0..nloc {
            for k in 0..d {
                g[(c * nloc + a) * d + k] *= scale[c];
            }
        }
    }
    g
}

/// `sign * K grad(u)` at the local dofs of every cell, not averaged, `dim` values per local dof.
pub(crate) fn cell_tensor_gradient(mesh: &SimplexMesh, u: &[f64], order: usize, k: &[M3], sign: f64) -> Vec<f64> {
    let d = mesh.dim;
    let nloc = mesh.nloc(order);
    let mut g = gradients_at_nodes(mesh, u, 1, order);
    for c in 0..mesh.n_cells() {
        for a in 0..nloc {
            let base = (c * nloc + a) * d;
            let gv: Vec<f64> = g[base..base + d].to_vec();
            for i in 0..d {
                g[base + i] = sign * (0..d).map(|j| k[c][i][j] * gv[j]).sum::<f64>();
            }
        }
    }
    g
}

/// The values of a dof field at the local dofs of every cell: `out[c*nloc + a]`.
pub(crate) fn cell_values_at_nodes(mesh: &SimplexMesh, u: &[f64], order: usize) -> Vec<f64> {
    let nloc = mesh.nloc(order);
    let mut out = Vec::with_capacity(mesh.n_cells() * nloc);
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        out.extend(dofs.iter().take(nloc).map(|&d| u[d]));
    }
    out
}

fn components(g: &[f64], d: usize) -> Vec<Vec<f64>> {
    (0..d).map(|k| g.chunks(d).map(|p| p[k]).collect()).collect()
}

fn norm_of(comps: &[Vec<f64>]) -> Vec<f64> {
    (0..comps[0].len()).map(|i| comps.iter().map(|c| c[i] * c[i]).sum::<f64>().sqrt()).collect()
}

fn minmax(v: &[f64]) -> (f64, f64) {
    v.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| (a.min(x), b.max(x)))
}

// ------------------------------------------------------------------ electrostatics

/// Stage 1 of electrostatics: the stiffness matrix, the load and the Dirichlet dofs.
pub struct ElectroAssembly {
    pub k: Csr,
    pub f: Vec<f64>,
    pub fixed: Vec<usize>,
    pub vals: Vec<f64>,
    /// fixed potential per boundary label
    pub volt: BTreeMap<i32, f64>,
}

pub fn assemble_electrostatics(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<ElectroAssembly, String> {
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let eps = spec.materials.cell_values(mesh, "eps_r", Some(1.0))?;
    let rho = spec.materials.cell_values(mesh, "rho_v", Some(0.0))?;
    let k = if spec.materials.has_tensor(TensorKind::Permittivity) {
        let mut t = plane_tensors(mesh, spec.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?);
        t.iter_mut().for_each(|m| m.iter_mut().for_each(|r| r.iter_mut().for_each(|x| *x *= EPS0)));
        assemble_laplace_tensor(mesh, order, &t)
    } else {
        assemble_laplace(mesh, order, &eps.iter().map(|e| EPS0 * e).collect::<Vec<_>>())
    };
    let n = k.nrows();
    let mut f = if rho.iter().any(|&r| r != 0.0) { assemble_source(mesh, order, &rho) } else { vec![0.0; n] };
    let mut volt: BTreeMap<i32, f64> = BTreeMap::new();
    for bc in &spec.bcs {
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Potential(v) => {
                    volt.insert(lab, *v);
                }
                BcKind::Ground => {
                    volt.insert(lab, 0.0);
                }
                BcKind::SurfaceCharge(s) => {
                    if *s != 0.0 {
                        let l = assemble_boundary_load(mesh, &[lab], *s, order);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                other => return Err(format!("'{}' is not an electrostatics boundary condition", other.name())),
            }
        }
    }
    if volt.is_empty() {
        return Err("electrostatics needs at least one potential or ground boundary".into());
    }
    let (fixed, vals) = dirichlet_constraints(mesh, &volt, order);
    Ok(ElectroAssembly { k, f, fixed, vals, volt })
}

/// Stage 2: the potential at every dof (and the note of the linear solver, if it has one). The backend is the one
/// `spec.solver` names.
pub fn solve_potential(mesh: &SimplexMesh, spec: &ProblemSpec, asm: &ElectroAssembly) -> Result<(Vec<f64>, Option<String>), String> {
    solve_spd(mesh, spec.order, 1, spec.solver, &asm.k, &asm.f, &asm.fixed, &asm.vals)
}

/// Stage 3: the plot fields `V`, `|E|`, `E{x,z}` and the same without averaging over the cells at a node.
pub fn electrostatics_fields(mesh: &SimplexMesh, order: usize, v: &[f64]) -> (Vec<Field>, CellFields) {
    let scale = vec![-1.0; mesh.n_cells()];
    let ecell = cell_gradient(mesh, v, order, &scale);
    let e = node_average(mesh, &ecell, mesh.dim, order);
    let ec = components(&e, mesh.dim);
    let mut fields = vec![Field { name: "V".into(), unit: "V".into(), values: v.to_vec() }, Field { name: "|E|".into(), unit: "V/m".into(), values: norm_of(&ec) }];
    for (name, c) in crate::section::component_names(mesh.dim).iter().zip(ec) {
        fields.push(Field { name: format!("E{name}"), unit: "V/m".into(), values: c });
    }
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("V", "V", cell_values_at_nodes(mesh, v, order));
    cf.push_components("E", crate::section::component_names(mesh.dim), "V/m", &ecell);
    (fields, cf)
}

pub fn solve_electrostatics(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let asm = assemble_electrostatics(mesh, spec)?;
    let (v, solver_note) = solve_potential(mesh, spec, &asm)?;
    let (fields, cell_fields) = electrostatics_fields(mesh, order, &v);
    let en = &fields[1].values;
    let kv = asm.k.matvec(&v);
    let w = 0.5 * v.iter().zip(&kv).map(|(a, b)| a * b).sum::<f64>();
    let (vmin, vmax) = minmax(&v);
    let mut scalars = vec![
        ("stored energy W_e (J; per metre in 2D)".to_string(), Scalar::Num(w)),
        ("max |E| (V/m)".into(), Scalar::Num(minmax(en).1)),
        ("V min (V)".into(), Scalar::Num(vmin)),
        ("V max (V)".into(), Scalar::Num(vmax)),
    ];
    let mut dv: Vec<f64> = asm.volt.values().map(|x| (x * 1e12).round() / 1e12).collect();
    dv.sort_by(|a, b| a.partial_cmp(b).unwrap());
    dv.dedup();
    if dv.len() == 2 {
        scalars.push(("capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(), Scalar::Num(2.0 * w / (dv[1] - dv[0]).powi(2))));
    }
    Ok(SolveResult {
        order,
        dim: mesh.dim,
        solutions: vec![Solution { label: "stationary".into(), fields, cell_fields: Some(cell_fields), ..Default::default() }],
        scalars,
        log: solver_note.into_iter().collect(),
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.k.nrows(),
        ..Default::default()
    })
}

// ------------------------------------------------------------------ heat transfer

/// Stage 1 of heat transfer: conductivity matrix, load, Dirichlet dofs and the conductivity tensors (for the flux).
pub struct HeatSystem {
    pub k: Csr,
    pub f: Vec<f64>,
    pub ktensor: Vec<M3>,
    pub fixed: Vec<usize>,
    pub vals: Vec<f64>,
    pub has_convection: bool,
}

pub fn assemble_heat(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<HeatSystem, String> {
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
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

/// Stage 3 of heat transfer: `T`, `|q|`, `q{x,z}` and the same per cell without averaging.
pub fn heat_fields(mesh: &SimplexMesh, t: &[f64], ktensor: &[M3], order: usize) -> (Vec<Field>, CellFields) {
    let qcell = cell_tensor_gradient(mesh, t, order, ktensor, -1.0);
    let q = node_average(mesh, &qcell, mesh.dim, order);
    let qc = components(&q, mesh.dim);
    let mut fields = vec![Field { name: "T".into(), unit: "K".into(), values: t.to_vec() }, Field { name: "|q|".into(), unit: "W/m^2".into(), values: norm_of(&qc) }];
    for (name, c) in crate::section::component_names(mesh.dim).iter().zip(qc) {
        fields.push(Field { name: format!("q{name}"), unit: "W/m^2".into(), values: c });
    }
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("T", "K", cell_values_at_nodes(mesh, t, order));
    cf.push_components("q", crate::section::component_names(mesh.dim), "W/m^2", &qcell);
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

pub fn solve_heat_states(mesh: &SimplexMesh, spec: &ProblemSpec, sys: &HeatSystem) -> Result<HeatStates, String> {
    let order = spec.order;
    let n = sys.k.nrows();
    let mut log = vec![];
    match &spec.study {
        Study::Stationary => {
            let (t, note) = solve_spd(mesh, order, 1, spec.solver, &sys.k, &sys.f, &sys.fixed, &sys.vals)?;
            log.extend(note);
            Ok(HeatStates { times: vec![], temps: vec![t], log })
        }
        Study::Transient { t_end, n_out, t0: t0_expr } => {
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
                LinearSolver::Auto => crate::fem::amg_wanted(mesh.dim, free.len()) && free.len() > 2500,
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
                    crate::fem::AmgSystem::new(mesh, order, 1, &a2_full, &sys.fixed).ok()
                } else {
                    None
                };
                let amg_sys = match amg_sys {
                    Some(a) if spec.solver == LinearSolver::Iterative || a.convergence_factor() < 0.7 => Some(a),
                    _ => None,
                };
                let a1_full = if amg_sys.is_some() { Some(m.scaled(1.0 / dt).add_scaled(&sys.k, 1.0)) } else { None };
                let (s1, s2) = if amg_sys.is_some() { (None, None) } else { (Some(DirectSolver::new(&a1, true)?), Some(DirectSolver::new(&a2, true)?)) };
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
                    log.push(format!("time steps by AMG-CG: {} levels {:?}, setup {:.1} s", a.n_levels(), a.sizes(), a.setup_seconds));
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
                log.push(format!("BDF2, {} steps ({} free dofs; too large for the exact dense method)", sub * n_out, free.len()));
            }
            Ok(HeatStates { times, temps: out, log })
        }
        _ => Err("heat transfer supports Stationary and Transient studies".into()),
    }
}

pub fn solve_heat(mesh: &SimplexMesh, spec: &ProblemSpec, state: &mut SolveState) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let sys = assemble_heat(mesh, spec)?;
    let st = solve_heat_states(mesh, spec, &sys)?;
    state.heat_t = st.temps.last().cloned();
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
        solutions.push(Solution { label, fields, cell_fields: Some(cell_fields), value, ..Default::default() });
    }
    Ok(SolveResult {
        order,
        dim: mesh.dim,
        solutions,
        scalars: vec![("T min (K)".into(), Scalar::Num(lo)), ("T max (K)".into(), Scalar::Num(hi))],
        log: st.log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: sys.k.nrows(),
        ..Default::default()
    })
}
