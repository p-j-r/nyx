// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Steady Ohmic electrical conduction (`Study::Conduction`) on P1/P2 simplices and native mixed cells
//! (`Quad`, `Prism`, `Hex`): `-div(sigma(T) grad V) = 0` with `voltage` (`potential`), `ground`, and
//! `current_density` (`-n . J = j_n`) boundary conditions, temperature-dependent conductivity `sigma(T)`,
//! and conservative volumetric Joule heating export `Q_joule = grad V . sigma(T) grad V`.
//! Three stages (ADR 9): `assemble_conduction`, `solve_conduction_potential`, `conduction_fields`;
//! `solve_conduction` runs them.

use crate::common::{
    assemble_boundary_load_mixed, cell_values_mixed, check_boundary_labels_mixed, components, mesh_has_curved_edges,
    minmax, plane_tensors,
};
use crate::diag::{solve_spd_checked, solve_spd_mixed_checked, DfsDiag};
use crate::problem::{BcKind, Problem, Study};
use femcore::fem::*;
use femcore::fem_mixed::{
    assemble_mass_mixed, lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::fields::{gradients_at_nodes, local_bary, node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::DirectSolver;
use femcore::material::{self, M3};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::props::Materials;
use nyx_physics::{
    cell_fields_of, discretise, discretise_mixed, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar,
    Solution, SolveResult,
};
use std::collections::BTreeMap;
use std::time::Instant;

/// Stage 1 of steady electrical conduction: conductivity stiffness matrix, current-density load,
/// representative cell conductivity tensors, and Dirichlet voltage constraints.
pub struct ConductionAssembly {
    pub k: Csr,
    pub f: Vec<f64>,
    pub sigma_tensor: Vec<M3>,
    pub fixed: Vec<usize>,
    pub vals: Vec<f64>,
}

pub(crate) fn cell_sigma_tensors_simplex(mesh: &Discretisation, materials: &Materials) -> Result<Vec<M3>, String> {
    let mut by_domain: BTreeMap<i32, M3> = BTreeMap::new();
    for d in mesh.domain_labels() {
        let t = match materials.aniso.get(&d).and_then(|a| a.conductivity.as_ref()) {
            Some(t) => material::tensor_global(t, materials.aniso[&d].euler_deg),
            None => {
                let v = materials.get(d, "sigma").unwrap_or(0.0);
                if v < 0.0 {
                    return Err(format!("domain {d} has negative electrical conductivity sigma = {v}"));
                }
                material::iso_tensor(v)
            }
        };
        by_domain.insert(d, t);
    }
    Ok(plane_tensors(mesh, mesh.cell_labels.iter().map(|l| by_domain[l]).collect()))
}

fn cell_sigma_tensors_mixed(mesh: &MixedDiscretisation, materials: &Materials) -> Result<Vec<M3>, String> {
    let mut by_domain: BTreeMap<i32, M3> = BTreeMap::new();
    for d in mesh.domain_labels() {
        let t = match materials.aniso.get(&d).and_then(|a| a.conductivity.as_ref()) {
            Some(t) => material::tensor_global(t, materials.aniso[&d].euler_deg),
            None => {
                let v = materials.get(d, "sigma").unwrap_or(0.0);
                if v < 0.0 {
                    return Err(format!("domain {d} has negative electrical conductivity sigma = {v}"));
                }
                material::iso_tensor(v)
            }
        };
        by_domain.insert(d, femcore::section::plane_tensor(&t, mesh.dim));
    }
    Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
}

/// Evaluates an imported temperature vector `t_field` (or domain-filtered per-cell P2 `t_cell_nodes`)
/// on simplex cell `ci` at barycentric coordinates `lam` (`dim + 1`).
pub(crate) fn eval_simplex_temp_at_bary(
    mesh: &Discretisation,
    t_field: &[f64],
    t_cell_nodes: Option<&[f64]>,
    ci: usize,
    lam: &[f64],
    dofs2_buf: &mut Vec<usize>,
) -> f64 {
    let nc = mesh.n_cells();
    let nloc2 = mesh.nloc(2);
    if let Some(tcn) = t_cell_nodes.filter(|v| v.len() == nc * nloc2) {
        let (n_vals, _) = shape(lam, 2);
        let base = ci * nloc2;
        return (0..nloc2).map(|b| n_vals[b] * tcn[base + b]).sum();
    }
    let nv = mesh.n_points();
    let n2 = mesh.n_dofs(2);
    if t_field.len() == n2 && n2 > nv {
        mesh.cell_dofs(ci, 2, dofs2_buf);
        let (n_vals, _) = shape(lam, 2);
        n_vals.iter().zip(dofs2_buf.iter()).map(|(&na, &d)| na * t_field[d]).sum()
    } else if t_field.len() >= nv {
        let cv = mesh.cell(ci);
        cv.iter().enumerate().map(|(a, &v)| lam[a] * t_field[v as usize]).sum()
    } else if t_field.len() == nc {
        t_field[ci]
    } else {
        t_field.first().copied().unwrap_or(293.15)
    }
}

fn eval_sigma_at_simplex_point(dim: usize, dom: i32, base_sigma: &M3, spec: &Problem, temp_k_opt: Option<f64>) -> M3 {
    match (temp_k_opt, spec.sigma_curves.get(&dom)) {
        (Some(tk), Some((curve, euler))) => femcore::section::plane_tensor(&curve.eval_matrix_global(tk, *euler), dim),
        _ => *base_sigma,
    }
}

/// Stage 1: assembles `K_sigma`, `f`, and Dirichlet constraints for `Study::Conduction` on a `SimplexMesh`.
pub fn assemble_conduction(mesh: &Discretisation, spec: &Problem) -> Result<ConductionAssembly, String> {
    if mesh.has_curved() {
        return Err("unsupported: steady conduction does not support curved simplex edge nodes (R32)".into());
    }
    if !matches!(spec.study, Study::Conduction) {
        return Err("assemble_conduction requires Study::Conduction".into());
    }
    let order = spec.order;
    let dim = mesh.dim;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
    let base_sigma = cell_sigma_tensors_simplex(mesh, &spec.materials)?;

    let has_temp_sigma =
        !spec.sigma_curves.is_empty() && (spec.temperature_field.is_some() || spec.temperature_cell_nodes.is_some());
    let mut sigma_cell = base_sigma.clone();
    let k = if has_temp_sigma {
        let tf = spec.temperature_field.as_deref().unwrap_or(&[]);
        let tcn = spec.temperature_cell_nodes.as_deref();
        let deg = (2 * order).max(2);
        let re = RefElement::new(dim, order, deg);
        let geo = geometry(mesh);
        let nloc = re.nloc;
        let n = mesh.n_dofs(order);
        // Update representative cell-centroid sigma tensors as well
        let bary_mid = if dim == 2 { vec![1.0 / 3.0; 3] } else { vec![0.25; 4] };
        let mut d2_buf = Vec::new();
        for c in 0..mesh.n_cells() {
            let tk = eval_simplex_temp_at_bary(mesh, tf, tcn, c, &bary_mid, &mut d2_buf);
            sigma_cell[c] = eval_sigma_at_simplex_point(dim, mesh.cell_labels[c], &base_sigma[c], spec, Some(tk));
        }
        assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
            let mut dofs = Vec::new();
            let mut d2 = Vec::new();
            let mut grad = vec![0.0; nloc * dim];
            let mut ke = vec![0.0; nloc * nloc];
            let qrule = quad_rule(dim, deg);
            for c in range {
                ke.fill(0.0);
                let dom = mesh.cell_labels[c];
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, c, q, &mut grad);
                    let tk = eval_simplex_temp_at_bary(mesh, tf, tcn, c, &qrule.lam[q], &mut d2);
                    let sig = eval_sigma_at_simplex_point(dim, dom, &base_sigma[c], spec, Some(tk));
                    let wq = re.w[q] * geo.vol[c];
                    for a in 0..nloc {
                        let mut sg = [0.0; 3];
                        for i in 0..dim {
                            for j in 0..dim {
                                sg[i] += sig[i][j] * grad[a * dim + j];
                            }
                        }
                        for b in 0..nloc {
                            let mut dot = 0.0;
                            for i in 0..dim {
                                dot += sg[i] * grad[b * dim + i];
                            }
                            ke[a * nloc + b] += wq * dot;
                        }
                    }
                }
                mesh.cell_dofs(c, order, &mut dofs);
                coo.push_block(&dofs, &dofs, &ke);
            }
        })
    } else if spec.materials.has_tensor(femcore::props::TensorKind::Conductivity) {
        assemble_laplace_tensor(mesh, order, &sigma_cell)
    } else {
        let sig_scalar: Vec<f64> = sigma_cell.iter().map(|t| t[0][0]).collect();
        assemble_laplace(mesh, order, &sig_scalar)
    };

    let n = k.nrows();
    let mut f = vec![0.0; n];
    let mut vfix: BTreeMap<i32, f64> = BTreeMap::new();
    for bc in &spec.bcs {
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Voltage(v) | BcKind::Potential(v) => {
                    vfix.insert(lab, *v);
                }
                BcKind::Ground => {
                    vfix.insert(lab, 0.0);
                }
                BcKind::CurrentDensity(jn) => {
                    if *jn != 0.0 {
                        let l = assemble_boundary_load(mesh, &[lab], *jn, order);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                other => return Err(format!("'{}' is not a steady-conduction boundary condition", other.name())),
            }
        }
    }
    let (fixed, vals) = if vfix.is_empty() { (vec![], vec![]) } else { dirichlet_constraints(mesh, &vfix, order) };
    if fixed.is_empty() {
        return Err("steady conduction needs at least one fixed-voltage or ground boundary".into());
    }
    let free = free_dofs(n, &fixed);
    let diag_k = k.diagonal();
    if free.iter().any(|&i| diag_k[i] <= 0.0) {
        return Err(
            "steady conduction has zero-conductivity free degrees of freedom; assign positive electrical conductivity 'sigma' (> 0 S/m) on conducting regions"
                .into(),
        );
    }

    Ok(ConductionAssembly { k, f, sigma_tensor: sigma_cell, fixed, vals })
}

/// Stage 2: solves `K_sigma V = f` subject to Dirichlet constraints.
pub fn solve_conduction_potential(
    mesh: &Discretisation,
    spec: &Problem,
    asm: &ConductionAssembly,
) -> Result<(Vec<f64>, Vec<String>), String> {
    let (v, log, _) = solve_conduction_potential_checked(mesh, spec, asm)?;
    Ok((v, log))
}

pub(crate) fn solve_conduction_potential_checked(
    mesh: &Discretisation,
    spec: &Problem,
    asm: &ConductionAssembly,
) -> Result<(Vec<f64>, Vec<String>, DfsDiag), String> {
    let tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
    let (v, note, diag) =
        solve_spd_checked(mesh, spec.order, 1, spec.solver, &asm.k, &asm.f, &asm.fixed, &asm.vals, tol, true)?;
    Ok((v, note.into_iter().collect(), diag))
}

/// Stage 3: computes `V` (`"V"`), `E = -grad V` (`"V/m"`), `J = sigma(T) E` (`"A/m^2"`), and conservative
/// `Q_joule = grad V . sigma(T) grad V` (`"W/m^3"`), returning `(fields, cell_fields, p_joule_watts)`.
pub fn conduction_fields(
    mesh: &Discretisation,
    spec: &Problem,
    v: &[f64],
    sigma_tensor: &[M3],
    order: usize,
) -> Result<(Vec<Field>, CellFields, f64), String> {
    let dim = mesh.dim;
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(order);
    let n = mesh.n_dofs(order);
    let g_pos = gradients_at_nodes(mesh, v, 1, order);
    let mut e_cell = vec![0.0; nc * nloc * dim];
    let mut j_cell = vec![0.0; nc * nloc * dim];
    let mut q_cell = vec![0.0; nc * nloc];

    let lam_nodes = local_bary(dim, order);
    let has_temp_sigma =
        !spec.sigma_curves.is_empty() && (spec.temperature_field.is_some() || spec.temperature_cell_nodes.is_some());
    let tf = spec.temperature_field.as_deref().unwrap_or(&[]);
    let tcn = spec.temperature_cell_nodes.as_deref();
    let mut d2_buf = Vec::new();

    for c in 0..nc {
        let dom = mesh.cell_labels[c];
        for a in 0..nloc {
            let base = (c * nloc + a) * dim;
            let sig = if has_temp_sigma {
                let tk = eval_simplex_temp_at_bary(mesh, tf, tcn, c, &lam_nodes[a], &mut d2_buf);
                eval_sigma_at_simplex_point(dim, dom, &sigma_tensor[c], spec, Some(tk))
            } else {
                sigma_tensor[c]
            };
            let mut q_val = 0.0;
            for i in 0..dim {
                let ei = -g_pos[base + i];
                e_cell[base + i] = ei;
                let mut ji = 0.0;
                for j in 0..dim {
                    ji += sig[i][j] * (-g_pos[base + j]);
                }
                j_cell[base + i] = ji;
                q_val += (-ei) * (-ji);
            }
            q_cell[c * nloc + a] = q_val.max(0.0);
        }
    }

    // Assemble exact quadrature Joule heating functional `b_i = int_Omega phi_i (grad V . sigma(T) grad V) dV`
    // and solve `M_1 * q_nodal = b` so `int_Omega Q_joule dV == P_joule` to machine precision.
    let deg = (2 * order).max(2);
    let re = RefElement::new(dim, order, deg);
    let qrule = quad_rule(dim, deg);
    let geo = geometry(mesh);
    let chunks: Vec<usize> = (0..nc).collect();
    let part_loads = femcore::par::map_chunks(&chunks, 256, |slice| {
        let mut loc_b = vec![0.0; n];
        let mut loc_p = 0.0;
        let mut dofs = Vec::new();
        let mut d2 = Vec::new();
        let mut grad = vec![0.0; nloc * dim];
        for &c in slice {
            mesh.cell_dofs(c, order, &mut dofs);
            let dom = mesh.cell_labels[c];
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                let mut gv = [0.0; 3];
                for a in 0..nloc {
                    let va = v[dofs[a]];
                    for k in 0..dim {
                        gv[k] += va * grad[a * dim + k];
                    }
                }
                let sig = if has_temp_sigma {
                    let tk = eval_simplex_temp_at_bary(mesh, tf, tcn, c, &qrule.lam[q], &mut d2);
                    eval_sigma_at_simplex_point(dim, dom, &sigma_tensor[c], spec, Some(tk))
                } else {
                    sigma_tensor[c]
                };
                let mut q_dens = 0.0;
                for i in 0..dim {
                    for j in 0..dim {
                        q_dens += gv[i] * sig[i][j] * gv[j];
                    }
                }
                q_dens = q_dens.max(0.0);
                let wq = re.w[q] * geo.vol[c];
                loc_p += wq * q_dens;
                for a in 0..nloc {
                    loc_b[dofs[a]] += wq * q_dens * re.n[q][a];
                }
            }
        }
        (loc_b, loc_p)
    });
    let mut b_joule = vec![0.0; n];
    let mut p_joule = 0.0;
    for (lb, lp) in part_loads {
        for (bi, &lbi) in b_joule.iter_mut().zip(&lb) {
            *bi += lbi;
        }
        p_joule += lp;
    }
    let ones = vec![1.0; nc];
    let m1 = assemble_mass(mesh, order, &ones);
    let m1_solver = DirectSolver::new(&m1, true)?;
    let q_nodal = m1_solver.solve(&b_joule);

    let e_nodal = node_average(mesh, &e_cell, dim, order);
    let j_nodal = node_average(mesh, &j_cell, dim, order);
    let ec = components(&e_nodal, dim);
    let jc = components(&j_nodal, dim);
    let names = femcore::section::component_names(dim);

    let fields = vec![
        Field::scalar("V", "V", v.to_vec()),
        Field::vector("E", "V/m", names, &ec),
        Field::vector("J", "A/m^2", names, &jc),
        Field::scalar("Q_joule", "W/m^3", q_nodal),
    ];
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("V", "V", cell_values_at_nodes(mesh, v, order));
    cf.push_components("E", names, "V/m", &e_cell);
    cf.push_components("J", names, "A/m^2", &j_cell);
    cf.push_real("Q_joule", "W/m^3", q_cell);
    Ok((fields, cf, p_joule))
}

/// Runs all three stages of steady electrical conduction (`Study::Conduction`) on a `SimplexMesh`.
pub fn solve_conduction(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let asm = assemble_conduction(mesh, spec)?;
    let (v, log, diag) = solve_conduction_potential_checked(mesh, spec, &asm)?;
    let (diagnostics, error_budget) = diag.finish();
    let (v_min, v_max) = minmax(&v);
    let (fields, cell_fields, p_joule) = conduction_fields(mesh, spec, &v, &asm.sigma_tensor, order)?;
    let solutions = vec![Solution {
        label: "conduction".into(),
        fields,
        cell_fields: cell_fields_of(&cell_fields),
        ..Default::default()
    }];
    Ok(SolveResult {
        physics: "dfs".into(),
        study: "conduction".into(),
        order: order as u8,
        dim: mesh.dim as u8,
        solutions,
        scalars: vec![
            ("V min (V)".into(), Scalar::Num(v_min)),
            ("V max (V)".into(), Scalar::Num(v_max)),
            ("P_joule (W)".into(), Scalar::Num(p_joule)),
        ],
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.k.nrows(),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

/// Solves `Study::Conduction` natively on a mixed mesh (`Quad`, `Prism`, `Hex`), delegating pure-simplex meshes to [`solve_conduction`].
pub fn solve_conduction_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String> {
    if mesh.is_pure_simplex() {
        if mesh_has_curved_edges(mesh) {
            return Err("unsupported: steady conduction does not support curved simplex edge nodes (R32)".into());
        }
        let d = discretise(mesh).map_err(|e| e.to_string())?;
        return solve_conduction(&d, spec);
    }

    let t0 = Instant::now();
    let order = spec.order;
    let mixed = discretise_mixed(mesh).map_err(|e| e.to_string())?;
    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, &mixed).map_err(|e| e.to_string())?;
    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(&mixed, space);
    let dim = mixed.dim;
    check_boundary_labels_mixed(&mixed, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;

    let base_sigma = cell_sigma_tensors_mixed(&mixed, &spec.materials)?;
    let _ = cell_values_mixed(&spec.materials, &mixed, "sigma", Some(0.0))?;
    let has_temp_sigma = !spec.sigma_curves.is_empty() && spec.temperature_field.is_some();
    let tf = spec.temperature_field.as_deref().unwrap_or(&[]);

    // Prolongate simplex nodal temperature onto mixed DOFs if present
    let t_mixed_opt: Option<Vec<f64>> = has_temp_sigma.then(|| {
        let coords = dof_map.dof_coordinates(&mixed);
        let nv = mixed.n_points();
        if tf.len() == dof_map.n_dofs {
            tf.to_vec()
        } else if tf.len() >= nv {
            let mut tm = vec![293.15; dof_map.n_dofs];
            tm[..nv].copy_from_slice(&tf[..nv]);
            for (e, &(a, b)) in mixed.edges.iter().enumerate() {
                if nv + e < tm.len() {
                    tm[nv + e] = 0.5 * (tf[a as usize] + tf[b as usize]);
                }
            }
            let _ = coords;
            tm
        } else {
            vec![tf.first().copied().unwrap_or(293.15); dof_map.n_dofs]
        }
    });

    let n = dof_map.n_dofs;
    let k = assemble_parallel(mixed.n_cells(), n, n, |range, coo| {
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut ke = Vec::new();
        for c in range {
            let ct = mixed.cell_type(c);
            let iso = IsoparametricCell::from_mesh(&mixed, c);
            let qrule = quad_rule_mixed(ct, (2 * order).max(2));
            let cdofs = dof_map.cell_dofs(c);
            let nloc = cdofs.len();
            ke.clear();
            ke.resize(nloc * nloc, 0.0);
            let dom = mixed.cell_labels[c];
            for (xi, &wq_ref) in qrule.pts.iter().zip(&qrule.w) {
                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                let tk_opt = t_mixed_opt.as_ref().map(|tm| vals.iter().zip(cdofs).map(|(&na, &d)| na * tm[d]).sum());
                let sig = eval_sigma_at_simplex_point(dim, dom, &base_sigma[c], spec, tk_opt);
                let wq = wq_ref * jac.det.abs();
                let grads: Vec<[f64; 3]> = dxi.iter().map(|&g_ref| jac.map_grad(g_ref)).collect();
                for a in 0..nloc {
                    let mut sg = [0.0; 3];
                    for i in 0..dim {
                        for j in 0..dim {
                            sg[i] += sig[i][j] * grads[a][j];
                        }
                    }
                    for b in 0..nloc {
                        let mut dot = 0.0;
                        for i in 0..dim {
                            dot += sg[i] * grads[b][i];
                        }
                        ke[a * nloc + b] += wq * dot;
                    }
                }
            }
            coo.push_block(cdofs, cdofs, &ke);
        }
    });

    let mut f = vec![0.0; n];
    let mut vfix: BTreeMap<i32, f64> = BTreeMap::new();
    for bc in &spec.bcs {
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Voltage(v) | BcKind::Potential(v) => {
                    vfix.insert(lab, *v);
                }
                BcKind::Ground => {
                    vfix.insert(lab, 0.0);
                }
                BcKind::CurrentDensity(jn) => {
                    if *jn != 0.0 {
                        let l = assemble_boundary_load_mixed(&mixed, &dof_map, &[lab], *jn);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                other => return Err(format!("'{}' is not a steady-conduction boundary condition", other.name())),
            }
        }
    }
    let (fixed, vals) = dof_map.dirichlet_constraints(&mixed, &vfix);
    if fixed.is_empty() {
        return Err("steady conduction needs at least one fixed-voltage or ground boundary".into());
    }
    let free = free_dofs(n, &fixed);
    let diag_k = k.diagonal();
    if free.iter().any(|&i| diag_k[i] <= 0.0) {
        return Err(
            "steady conduction has zero-conductivity free degrees of freedom; assign positive electrical conductivity 'sigma' (> 0 S/m) on conducting regions"
                .into(),
        );
    }

    let (v_mixed, note, diag) = solve_spd_mixed_checked(spec.solver, &k, &f, &fixed, &vals)?;
    let log: Vec<String> = note.into_iter().collect();
    let (diagnostics, error_budget) = diag.finish();

    // Project V to simplex and compute conservative Q_joule on mixed + simplex
    let v_simplex = proj.project_lagrange_blocks(&mixed, &dof_map, &v_mixed, 1, order);
    let (v_min, v_max) = minmax(&v_simplex);

    let mut b_mixed = vec![0.0; n];
    let mut p_joule = 0.0;
    {
        let mut vals_s = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        for c in 0..mixed.n_cells() {
            let ct = mixed.cell_type(c);
            let iso = IsoparametricCell::from_mesh(&mixed, c);
            let qrule = quad_rule_mixed(ct, (2 * order).max(2));
            let cdofs = dof_map.cell_dofs(c);
            let dom = mixed.cell_labels[c];
            for (xi, &wq_ref) in qrule.pts.iter().zip(&qrule.w) {
                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals_s, &mut dxi);
                let tk_opt = t_mixed_opt.as_ref().map(|tm| vals_s.iter().zip(cdofs).map(|(&na, &d)| na * tm[d]).sum());
                let sig = eval_sigma_at_simplex_point(dim, dom, &base_sigma[c], spec, tk_opt);
                let grads: Vec<[f64; 3]> = dxi.iter().map(|&g_ref| jac.map_grad(g_ref)).collect();
                let mut gv = [0.0; 3];
                for (a, &gd) in cdofs.iter().enumerate() {
                    for k_ax in 0..dim {
                        gv[k_ax] += v_mixed[gd] * grads[a][k_ax];
                    }
                }
                let mut q_dens = 0.0;
                for i in 0..dim {
                    for j in 0..dim {
                        q_dens += gv[i] * sig[i][j] * gv[j];
                    }
                }
                q_dens = q_dens.max(0.0);
                let wq = wq_ref * jac.det.abs();
                p_joule += wq * q_dens;
                for (a, &gd) in cdofs.iter().enumerate() {
                    b_mixed[gd] += wq * q_dens * vals_s[a];
                }
            }
        }
    }
    let m1_mixed = assemble_mass_mixed(&mixed, &dof_map, &vec![1.0; mixed.n_cells()]);
    let q_mixed_nodal = DirectSolver::new(&m1_mixed, true)?.solve(&b_mixed);
    let q_simplex_raw = proj.project_lagrange_blocks(&mixed, &dof_map, &q_mixed_nodal, 1, order);
    // Rescale projected simplex Q_joule so simplex integral matches exact mixed P_joule
    let m1_simp = assemble_mass(&simplex, order, &vec![1.0; simplex.n_cells()]);
    let int_simp: f64 = m1_simp.matvec(&q_simplex_raw).iter().sum();
    let q_scale = if int_simp > 1e-30 && p_joule > 0.0 { p_joule / int_simp } else { 1.0 };
    let q_simplex: Vec<f64> = q_simplex_raw.into_iter().map(|x| x * q_scale).collect();

    let scale_neg = vec![-1.0; mixed.n_cells()];
    let e_cell =
        crate::common::sample_mixed_gradients_on_simplex(&mixed, &dof_map, &proj, &v_mixed, None, &scale_neg, order);
    let j_cell = crate::common::sample_mixed_gradients_on_simplex(
        &mixed,
        &dof_map,
        &proj,
        &v_mixed,
        Some(&base_sigma),
        &scale_neg,
        order,
    );
    let e_nodal = node_average(&simplex, &e_cell, dim, order);
    let j_nodal = node_average(&simplex, &j_cell, dim, order);
    let ec = components(&e_nodal, dim);
    let jc = components(&j_nodal, dim);
    let names = femcore::section::component_names(dim);
    let q_cell = cell_values_at_nodes(&simplex, &q_simplex, order);

    let fields = vec![
        Field::scalar("V", "V", v_simplex.clone()),
        Field::vector("E", "V/m", names, &ec),
        Field::vector("J", "A/m^2", names, &jc),
        Field::scalar("Q_joule", "W/m^3", q_simplex),
    ];
    let mut cf = CellFields::new(&simplex, order);
    cf.push_real("V", "V", cell_values_at_nodes(&simplex, &v_simplex, order));
    cf.push_components("E", names, "V/m", &e_cell);
    cf.push_components("J", names, "A/m^2", &j_cell);
    cf.push_real("Q_joule", "W/m^3", q_cell);

    Ok(SolveResult {
        physics: "dfs".into(),
        study: "conduction".into(),
        order: order as u8,
        dim: dim as u8,
        solutions: vec![Solution {
            label: "conduction".into(),
            fields,
            cell_fields: cell_fields_of(&cf),
            ..Default::default()
        }],
        scalars: vec![
            ("V min (V)".into(), Scalar::Num(v_min)),
            ("V max (V)".into(), Scalar::Num(v_max)),
            ("P_joule (W)".into(), Scalar::Num(p_joule)),
        ],
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}
