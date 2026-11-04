//! Electrostatics on P1/P2 simplices: `-div(eps0 eps_r grad V) = rho_v`, `V = V0` on fixed-potential boundaries, `n.D = -rho_s` on
//! charged ones. Three stages (ADR 9): `assemble_electrostatics`, `solve_potential`, `electrostatics_fields`; `solve_electrostatics` runs them.

use crate::common::{cell_gradient, components, minmax, norm_of, plane_tensors};
use crate::problem::{BcKind, Problem};
use femcore::fem::*;
use femcore::fields::{node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::props::TensorKind;
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use std::collections::BTreeMap;
use std::time::Instant;

pub const EPS0: f64 = 8.8541878128e-12;

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

pub fn assemble_electrostatics(mesh: &Discretisation, spec: &Problem) -> Result<ElectroAssembly, String> {
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
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
pub fn solve_potential(
    mesh: &Discretisation,
    spec: &Problem,
    asm: &ElectroAssembly,
) -> Result<(Vec<f64>, Option<String>), String> {
    solve_spd(mesh, spec.order, 1, spec.solver, &asm.k, &asm.f, &asm.fixed, &asm.vals)
}

/// Stage 3: the fields `V` and `E` (the flat names `|E|`, `Ex`, `Ez` come from `Solution::flat_fields`) and the same without averaging over
/// the cells at a node.
pub fn electrostatics_fields(mesh: &Discretisation, order: usize, v: &[f64]) -> (Vec<Field>, CellFields) {
    let scale = vec![-1.0; mesh.n_cells()];
    let ecell = cell_gradient(mesh, v, order, &scale);
    let e = node_average(mesh, &ecell, mesh.dim, order);
    let ec = components(&e, mesh.dim);
    let names = femcore::section::component_names(mesh.dim);
    let fields = vec![Field::scalar("V", "V", v.to_vec()), Field::vector("E", "V/m", names, &ec)];
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("V", "V", cell_values_at_nodes(mesh, v, order));
    cf.push_components("E", femcore::section::component_names(mesh.dim), "V/m", &ecell);
    (fields, cf)
}

pub fn solve_electrostatics(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let asm = assemble_electrostatics(mesh, spec)?;
    let (v, solver_note) = solve_potential(mesh, spec, &asm)?;
    let (fields, cell_fields) = electrostatics_fields(mesh, order, &v);
    let en = norm_of(&components(&fields[1].values, mesh.dim));
    let kv = asm.k.matvec(&v);
    let w = 0.5 * v.iter().zip(&kv).map(|(a, b)| a * b).sum::<f64>();
    let (vmin, vmax) = minmax(&v);
    let mut scalars = vec![
        ("stored energy W_e (J; per metre in 2D)".to_string(), Scalar::Num(w)),
        ("max |E| (V/m)".into(), Scalar::Num(minmax(&en).1)),
        ("V min (V)".into(), Scalar::Num(vmin)),
        ("V max (V)".into(), Scalar::Num(vmax)),
    ];
    let mut dv: Vec<f64> = asm.volt.values().map(|x| (x * 1e12).round() / 1e12).collect();
    dv.sort_by(|a, b| a.partial_cmp(b).unwrap());
    dv.dedup();
    if dv.len() == 2 {
        scalars.push((
            "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(),
            Scalar::Num(2.0 * w / (dv[1] - dv[0]).powi(2)),
        ));
    }
    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: mesh.dim as u8,
        solutions: vec![Solution {
            label: "stationary".into(),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            ..Default::default()
        }],
        scalars,
        log: solver_note.into_iter().collect(),
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.k.nrows(),
    })
}
