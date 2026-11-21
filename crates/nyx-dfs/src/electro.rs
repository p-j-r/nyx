// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Electrostatics on P1/P2 simplices: `-div(eps0 eps_r grad V) = rho_v`, `V = V0` on fixed-potential boundaries, `n.D = -rho_s` on
//! charged ones. Three stages (ADR 9): `assemble_electrostatics`, `solve_potential`, `electrostatics_fields`; `solve_electrostatics` runs them.

use crate::common::{cell_gradient, components, minmax, norm_of, plane_tensors};
use crate::problem::{BcKind, Problem};
use femcore::fem::*;
use femcore::fields::{node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::{dense_eigh, DirectSolver};
use femcore::material::{self, M3};
use femcore::props::{LinearSolver, TensorKind};
use nyx_physics::{
    cell_fields_of, CapacitanceMatrixData, Discretisation, Field, Scalar, Solution, SolveResult, TlsMeta,
};
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

pub const EPS0: f64 = 8.8541878128e-12;
const E_CHARGE: f64 = 1.602_176_634e-19;
const PLANCK_H: f64 = 6.626_070_15e-34;

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

fn assemble_laplace_curved_3d(mesh: &Discretisation, order: usize, coeff: &[f64]) -> Csr {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(6);
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * 3];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.iter_mut().for_each(|x| *x = 0.0);
            if mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    let dn = &cshapes[q];
                    for a in 0..nloc {
                        for k in 0..3 {
                            let mut s = 0.0;
                            for b in 0..4 {
                                s += dn[a * 4 + b] * gq[b][k];
                            }
                            grad[a * 3 + k] = s;
                        }
                    }
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let s = grad[a * 3] * grad[b * 3]
                                + grad[a * 3 + 1] * grad[b * 3 + 1]
                                + grad[a * 3 + 2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += wq * s;
                        }
                    }
                }
                let sc = coeff[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            } else {
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, c, q, &mut grad);
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let s = grad[a * 3] * grad[b * 3]
                                + grad[a * 3 + 1] * grad[b * 3 + 1]
                                + grad[a * 3 + 2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += re.w[q] * s;
                        }
                    }
                }
                let sc = geo.vol[c] * coeff[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            }
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

fn assemble_laplace_tensor_curved_3d(mesh: &Discretisation, order: usize, kt: &[M3]) -> Csr {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(6);
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut dofs = Vec::new();
        let mut grad = vec![0.0; nloc * 3];
        let mut ke = vec![0.0; nloc * nloc];
        for c in range {
            ke.iter_mut().for_each(|x| *x = 0.0);
            if mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    let dn = &cshapes[q];
                    for a in 0..nloc {
                        for k in 0..3 {
                            let mut s = 0.0;
                            for b in 0..4 {
                                s += dn[a * 4 + b] * gq[b][k];
                            }
                            grad[a * 3 + k] = s;
                        }
                    }
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    for a in 0..nloc {
                        let mut kg = [0.0; 3];
                        for i in 0..3 {
                            for j in 0..3 {
                                kg[i] += kt[c][i][j] * grad[a * 3 + j];
                            }
                        }
                        for b in 0..nloc {
                            let s = kg[0] * grad[b * 3] + kg[1] * grad[b * 3 + 1] + kg[2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += wq * s;
                        }
                    }
                }
            } else {
                for q in 0..re.w.len() {
                    shape_gradients(&re, &geo, c, q, &mut grad);
                    for a in 0..nloc {
                        let mut kg = [0.0; 3];
                        for i in 0..3 {
                            for j in 0..3 {
                                kg[i] += kt[c][i][j] * grad[a * 3 + j];
                            }
                        }
                        for b in 0..nloc {
                            let s = kg[0] * grad[b * 3] + kg[1] * grad[b * 3 + 1] + kg[2] * grad[b * 3 + 2];
                            ke[a * nloc + b] += re.w[q] * s;
                        }
                    }
                }
                let sc = geo.vol[c];
                ke.iter_mut().for_each(|x| *x *= sc);
            }
            mesh.cell_dofs(c, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    })
}

/// Enriches a 3D mesh for P2 Barsoum quarter-point singularity capture (`singular_edges = "quarter_point"`):
/// finds all open perimeter vertices of interior conductor sheet facets (edges shared by exactly 1
/// interior conductor facet carrying `Terminal`, `Ground`, or `Potential`), shifts mid-edge geometry
/// nodes on tet edges `(a, b)` with `is_rim[a] ^ is_rim[b]` to `0.75 * p_rim + 0.25 * p_other`, and
/// validates positive `CurvedTet` Jacobians across all collapsed-quadrature points on every affected cell.
pub fn enrich_quarter_point_mesh(mesh: &Discretisation, spec: &Problem) -> Option<(Discretisation, usize, usize)> {
    if spec.order != 2 || mesh.dim != 3 || spec.singular_edges != "quarter_point" {
        return None;
    }
    let mut cond_labels = BTreeSet::new();
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Potential(_) | BcKind::Ground | BcKind::Terminal(_)) {
            cond_labels.extend(bc.labels.iter().copied());
        }
    }
    if cond_labels.is_empty() {
        return None;
    }
    let exterior = mesh.exterior_facets();
    let mut edge_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for f in 0..mesh.n_facets() {
        if !exterior[f] && cond_labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            for &(i, j) in &[(0, 1), (1, 2), (2, 0)] {
                let (a, b) = (fv[i].min(fv[j]), fv[i].max(fv[j]));
                *edge_count.entry((a, b)).or_insert(0) += 1;
            }
        }
    }
    if edge_count.is_empty() {
        return None;
    }
    let nv = mesh.n_points();
    let mut is_rim = vec![false; nv];
    for (&(a, b), &cnt) in &edge_count {
        if cnt == 1 {
            is_rim[a as usize] = true;
            is_rim[b as usize] = true;
        }
    }
    let n_rim_verts = is_rim.iter().filter(|&&r| r).count();
    if n_rim_verts == 0 {
        return None;
    }

    let mut edge_map: BTreeMap<[u32; 2], [f64; 3]> = mesh.edge_node_list().into_iter().collect();
    let mut qp_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
    for &(a, b) in mesh.edges() {
        let (ua, ub) = (a as usize, b as usize);
        if (is_rim[ua] ^ is_rim[ub]) && !edge_map.contains_key(&[a, b]) {
            let (p_rim, p_other) =
                if is_rim[ua] { (mesh.point(ua), mesh.point(ub)) } else { (mesh.point(ub), mesh.point(ua)) };
            let qpt = [
                0.75 * p_rim[0] + 0.25 * p_other[0],
                0.75 * p_rim[1] + 0.25 * p_other[1],
                0.75 * p_rim[2] + 0.25 * p_other[2],
            ];
            edge_map.insert([a, b], qpt);
            qp_edges.insert([a, b]);
        }
    }
    if qp_edges.is_empty() {
        return None;
    }

    let mut enriched = mesh.clone();
    let nodes_vec: Vec<([u32; 2], [f64; 3])> = edge_map.iter().map(|(&k, &v)| (k, v)).collect();
    enriched.set_edge_nodes(&nodes_vec);

    // Validate positive CurvedTet Jacobian across all quadrature points on every affected cell
    let cquad = quad_rule_collapsed(6);
    let mut bad_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
    for c in 0..enriched.n_cells() {
        if !enriched.cell_is_curved(c) {
            continue;
        }
        let cv = enriched.cell(c);
        let mut v4 = [[0.0; 3]; 4];
        for i in 0..4 {
            let p = enriched.point(cv[i] as usize);
            v4[i] = [p[0], p[1], p[2]];
        }
        let det0 = CurvedTet::straight(v4).jacobian(&[0.25, 0.25, 0.25, 0.25]).1;
        let ct = CurvedTet::new(enriched.cell_geometry_nodes(c));
        let mut valid = det0.abs() > 1e-30;
        if valid {
            for lam in &cquad.lam {
                let (_, det) = ct.jacobian(&[lam[0], lam[1], lam[2], lam[3]]);
                if det * det0 <= 1e-12 * det0 * det0 {
                    valid = false;
                    break;
                }
            }
        }
        if !valid {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    let key = [cv[i].min(cv[j]), cv[i].max(cv[j])];
                    if qp_edges.contains(&key) {
                        bad_edges.insert(key);
                    }
                }
            }
        }
    }
    if !bad_edges.is_empty() {
        for k in &bad_edges {
            edge_map.remove(k);
            qp_edges.remove(k);
        }
        if qp_edges.is_empty() {
            return None;
        }
        let nodes_vec: Vec<([u32; 2], [f64; 3])> = edge_map.iter().map(|(&k, &v)| (k, v)).collect();
        enriched.set_edge_nodes(&nodes_vec);
    }
    Some((enriched, n_rim_verts, qp_edges.len()))
}

pub fn assemble_electrostatics(mesh: &Discretisation, spec: &Problem) -> Result<ElectroAssembly, String> {
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
    let qp_mesh = enrich_quarter_point_mesh(mesh, spec);
    let m_work = qp_mesh.as_ref().map_or(mesh, |(m, _, _)| m);
    let eps = spec.materials.cell_values(m_work, "eps_r", Some(1.0))?;
    let rho = spec.materials.cell_values(m_work, "rho_v", Some(0.0))?;
    let k = if spec.materials.has_tensor(TensorKind::Permittivity) {
        let mut t = plane_tensors(m_work, spec.materials.cell_tensor(m_work, TensorKind::Permittivity, Some(1.0))?);
        t.iter_mut().for_each(|m| m.iter_mut().for_each(|r| r.iter_mut().for_each(|x| *x *= EPS0)));
        if m_work.dim == 3 && m_work.has_curved() {
            assemble_laplace_tensor_curved_3d(m_work, order, &t)
        } else {
            assemble_laplace_tensor(m_work, order, &t)
        }
    } else {
        let coeff: Vec<f64> = eps.iter().map(|e| EPS0 * e).collect();
        if m_work.dim == 3 && m_work.has_curved() {
            assemble_laplace_curved_3d(m_work, order, &coeff)
        } else {
            assemble_laplace(m_work, order, &coeff)
        }
    };
    let n = k.nrows();
    let mut f = if rho.iter().any(|&r| r != 0.0) { assemble_source(m_work, order, &rho) } else { vec![0.0; n] };
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
                BcKind::Terminal(name) => {
                    if name.trim().is_empty() {
                        return Err("terminal condition requires a non-empty 'name'".into());
                    }
                    volt.entry(lab).or_insert(0.0);
                }
                BcKind::SurfaceCharge(s) => {
                    if *s != 0.0 {
                        let l = assemble_boundary_load(m_work, &[lab], *s, order);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                other => return Err(format!("'{}' is not an electrostatics boundary condition", other.name())),
            }
        }
    }
    if volt.is_empty() {
        return Err("electrostatics needs at least one potential, ground or terminal boundary".into());
    }
    let (fixed, vals) = dirichlet_constraints(m_work, &volt, order);
    Ok(ElectroAssembly { k, f, fixed, vals, volt })
}

/// Stage 2: the potential at every dof (and the note of the linear solver, if it has one). The backend is the one
/// `spec.solver` names.
pub fn solve_potential(
    mesh: &Discretisation,
    spec: &Problem,
    asm: &ElectroAssembly,
) -> Result<(Vec<f64>, Option<String>), String> {
    let tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
    solve_spd_with(mesh, spec.order, 1, spec.solver, &asm.k, &asm.f, &asm.fixed, &asm.vals, tol, spec.amg_line_smooth)
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

fn build_tls_meta(mesh: &Discretisation, spec: &Problem) -> TlsMeta {
    let mut domain_eps_r = BTreeMap::new();
    for d in mesh.domain_labels() {
        let eps = if let Some(t) = spec.materials.aniso.get(&d).and_then(|a| a.permittivity.as_ref()) {
            let m = material::tensor_global(t, spec.materials.aniso[&d].euler_deg);
            (m[0][0] + m[1][1] + m[2][2]) / 3.0
        } else {
            spec.materials.get(d, "eps_r").unwrap_or(1.0)
        };
        domain_eps_r.insert(d, eps);
    }
    let mut cond_set = BTreeSet::new();
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Potential(_) | BcKind::Ground | BcKind::Terminal(_)) {
            cond_set.extend(bc.labels.iter().copied());
        }
    }
    let conductor_labels: Vec<i32> = cond_set.iter().copied().collect();
    let exterior = mesh.exterior_facets();
    let mut sheet_set = BTreeSet::new();
    let mut pec_set = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if cond_set.contains(&lab) {
            if exterior[f] {
                pec_set.insert(lab);
            } else {
                sheet_set.insert(lab);
            }
        }
    }
    TlsMeta {
        domain_eps_r,
        domain_names: spec.domain_names.clone(),
        conductor_labels,
        sheet_labels: sheet_set.into_iter().collect(),
        sheet_thickness_m: BTreeMap::new(),
        pec_labels: pec_set.into_iter().collect(),
    }
}

fn invert_spd_ec_matrix(maxwell: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = maxwell.len();
    if n == 0 {
        return Ok(vec![]);
    }
    let mut flat = Vec::with_capacity(n * n);
    for row in maxwell {
        flat.extend_from_slice(row);
    }
    let (eigs, vecs) = dense_eigh(&flat, n)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; n]; n];
    if lam_max <= 0.0 {
        return Ok(ec);
    }
    for k in 0..n {
        if eigs[k] > 1e-13 * lam_max {
            let inv_l = scale / eigs[k];
            for i in 0..n {
                for j in 0..n {
                    ec[i][j] += inv_l * vecs[k][i] * vecs[k][j];
                }
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let s = 0.5 * (ec[i][j] + ec[j][i]);
            ec[i][j] = s;
            ec[j][i] = s;
        }
    }
    Ok(ec)
}

pub fn solve_electrostatics(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    match spec.formulation.as_str() {
        "" | "volumetric_fem" => {}
        "planar_bem" => return crate::planar_bem::solve_planar_bem(mesh, spec),
        "auto" => {
            if crate::planar_bem::can_use_planar_bem(mesh, spec) {
                return crate::planar_bem::solve_planar_bem(mesh, spec);
            }
        }
        other => return Err(format!("unknown electrostatics formulation '{other}'")),
    }
    let t0 = Instant::now();
    let order = spec.order;
    let qp_info = enrich_quarter_point_mesh(mesh, spec);
    let asm = assemble_electrostatics(mesh, spec)?;

    // Group terminal facets by unique terminal name in first-seen order
    let mut terminals: Vec<(String, Vec<i32>)> = Vec::new();
    for bc in &spec.bcs {
        if let BcKind::Terminal(name) = &bc.kind {
            let name = name.trim();
            if name.is_empty() {
                return Err("terminal condition requires a non-empty 'name'".into());
            }
            if let Some((_, labs)) = terminals.iter_mut().find(|(n, _)| n == name) {
                for &l in &bc.labels {
                    if !labs.contains(&l) {
                        labs.push(l);
                    }
                }
            } else {
                terminals.push((name.to_string(), bc.labels.clone()));
            }
        }
    }

    let tls = build_tls_meta(mesh, spec);

    if terminals.is_empty() {
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
        scalars.push(tls.scalar());
        return Ok(SolveResult {
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
        });
    }

    // Multi-RHS terminal capacitance matrix extraction (N >= 1 terminals):
    let n_term = terminals.len();
    let n = asm.k.nrows();
    let term_dofs: Vec<BTreeSet<usize>> =
        terminals.iter().map(|(_, labs)| boundary_dofs(mesh, labs, order).into_iter().collect()).collect();
    let mut bg_labels = Vec::new();
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Ground | BcKind::Potential(_)) {
            for &l in &bc.labels {
                if !bg_labels.contains(&l) {
                    bg_labels.push(l);
                }
            }
        }
    }
    let bg_dofs: BTreeSet<usize> = boundary_dofs(mesh, &bg_labels, order).into_iter().collect();
    if n_term == 1 && bg_dofs.is_empty() {
        return Err(
            "electrostatics with a single terminal needs at least one ground or potential reference boundary".into()
        );
    }

    let mut all_fixed_set = bg_dofs.clone();
    for td in &term_dofs {
        all_fixed_set.extend(td.iter().copied());
    }
    let fixed: Vec<usize> = all_fixed_set.into_iter().collect();
    let free = free_dofs(n, &fixed);
    let n_free = free.len();

    let vals_list: Vec<Vec<f64>> = (0..n_term)
        .map(|j| {
            fixed.iter().map(|&d| if term_dofs[j].contains(&d) && !bg_dofs.contains(&d) { 1.0 } else { 0.0 }).collect()
        })
        .collect();

    let use_amg = match spec.solver {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => true,
        LinearSolver::Auto => amg_wanted(mesh.dim, n_free),
    } && !fixed.is_empty()
        && n_free > 0;

    let mut log = Vec::new();
    if let Some((_, n_rim, n_edges)) = &qp_info {
        log.push(format!(
            "singular_edges=quarter_point: shifted {n_edges} tet edges meeting {n_rim} conductor sheet perimeter vertices to Barsoum 1/4-points (all affected cells verified det(J) > 0)"
        ));
    }
    let mut v_all = vec![vec![0.0; n]; n_term];
    for j in 0..n_term {
        for (&d, &val) in fixed.iter().zip(&vals_list[j]) {
            v_all[j][d] = val;
        }
    }

    if n_free > 0 {
        let solve_direct = |log: &mut Vec<String>| -> Result<Vec<Vec<f64>>, String> {
            let kff = asm.k.submatrix(&free, &free);
            let kfx = asm.k.submatrix(&free, &fixed);
            let rhs_list: Vec<Vec<f64>> = vals_list
                .iter()
                .map(|vals| {
                    let kx = kfx.matvec(vals);
                    free.iter().enumerate().map(|(a, &i)| asm.f[i] - kx[a]).collect()
                })
                .collect();
            let solver = DirectSolver::new(&kff, true)?;
            let sols = solver.solve_many(&rhs_list);
            log.push(format!(
                "direct Cholesky: 1 factorization ({n_free} free dofs) served {n_term} terminal right-hand side(s)"
            ));
            Ok(sols)
        };

        if use_amg {
            let opts =
                femcore::linalg::amg_elastic::AmgOptions { line_smooth: spec.amg_line_smooth, ..Default::default() };
            let sys = AmgSystem::new_with(mesh, order, 1, &asm.k, &fixed, opts)?;
            let tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
            let mut rhs_list = Vec::with_capacity(n_term);
            let mut fv_list = Vec::with_capacity(n_term);
            for j in 0..n_term {
                let ku = asm.k.matvec(&v_all[j]);
                let rhs: Vec<f64> = (0..n).map(|i| asm.f[i] - ku[i]).collect();
                let fv: Vec<(usize, f64)> = fixed.iter().copied().zip(vals_list[j].iter().copied()).collect();
                rhs_list.push(rhs);
                fv_list.push(fv);
            }
            let (u_all, its, worst_rel) = sys.solve_block_with(None, &rhs_list, &fv_list, tol, 600);
            let amg_ok = worst_rel <= (tol * 1e3).max(1e-7);
            if amg_ok {
                for j in 0..n_term {
                    for &i in &free {
                        v_all[j][i] = u_all[j][i];
                    }
                }
                log.push(format!(
                    "AMG-CG: 1 hierarchy setup ({} levels {:?}{}, {:.1} s) served {n_term} terminal right-hand side(s) ({its} block iterations, max residual {worst_rel:.1e})",
                    sys.n_levels(),
                    sys.sizes(),
                    if sys.line_smooth_active() {
                        format!(", {} vertical-line columns", sys.n_line_columns())
                    } else {
                        String::new()
                    },
                    sys.setup_seconds
                ));
            } else if spec.solver == LinearSolver::Auto {
                let sols = solve_direct(&mut log)?;
                for j in 0..n_term {
                    for (a, &i) in free.iter().enumerate() {
                        v_all[j][i] = sols[j][a];
                    }
                }
            } else {
                return Err(format!(
                    "the iterative solver did not converge (residual {worst_rel:.2e}); choose the direct solver"
                ));
            }
        } else {
            let sols = solve_direct(&mut log)?;
            for j in 0..n_term {
                for (a, &i) in free.iter().enumerate() {
                    v_all[j][i] = sols[j][a];
                }
            }
        }
    } else {
        log.push(format!("0 free dofs: {n_term} terminal right-hand side(s) determined by boundary constraints"));
    }

    let kv_all: Vec<Vec<f64>> = v_all.iter().map(|v| asm.k.matvec(v)).collect();
    let mut maxwell_matrix_f = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in i..n_term {
            let cij = femcore::linalg::vecops::dotu(&v_all[i], &kv_all[j]);
            let cji = if i == j { cij } else { femcore::linalg::vecops::dotu(&v_all[j], &kv_all[i]) };
            let sym = 0.5 * (cij + cji);
            maxwell_matrix_f[i][j] = sym;
            maxwell_matrix_f[j][i] = sym;
        }
    }

    let mut mutual_matrix_f = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        let mut sum_row = 0.0;
        for k in 0..n_term {
            sum_row += maxwell_matrix_f[i][k];
            if i != k {
                mutual_matrix_f[i][k] = -maxwell_matrix_f[i][k];
            }
        }
        mutual_matrix_f[i][i] = sum_row;
    }

    let ec_matrix_hz = invert_spd_ec_matrix(&maxwell_matrix_f)?;

    let mut solutions = Vec::with_capacity(n_term);
    let mut scalars = Vec::new();
    for j in 0..n_term {
        let (fields, cell_fields) = electrostatics_fields(mesh, order, &v_all[j]);
        if j == 0 {
            let en = norm_of(&components(&fields[1].values, mesh.dim));
            let w = 0.5 * maxwell_matrix_f[0][0];
            let (vmin, vmax) = minmax(&v_all[0]);
            scalars.push(("stored energy W_e (J; per metre in 2D)".to_string(), Scalar::Num(w)));
            scalars.push(("max |E| (V/m)".into(), Scalar::Num(minmax(&en).1)));
            scalars.push(("V min (V)".into(), Scalar::Num(vmin)));
            scalars.push(("V max (V)".into(), Scalar::Num(vmax)));
            if n_term == 1 {
                scalars.push((
                    "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(),
                    Scalar::Num(maxwell_matrix_f[0][0]),
                ));
            }
        }
        solutions.push(Solution {
            label: format!("terminal {}", terminals[j].0),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(maxwell_matrix_f[j][j]),
            ..Default::default()
        });
    }

    for i in 0..n_term {
        for j in 0..n_term {
            let name_key = format!("C_{}_{} (F)", terminals[i].0, terminals[j].0);
            let idx_key = format!("C_{i}_{j} (F)");
            scalars.push((name_key.clone(), Scalar::Num(maxwell_matrix_f[i][j])));
            if idx_key != name_key {
                scalars.push((idx_key, Scalar::Num(maxwell_matrix_f[i][j])));
            }
        }
    }
    for i in 0..n_term {
        for j in 0..n_term {
            let name_key = format!("E_C_{}_{} (Hz)", terminals[i].0, terminals[j].0);
            let idx_key = format!("E_C_{i}_{j} (Hz)");
            scalars.push((name_key.clone(), Scalar::Num(ec_matrix_hz[i][j])));
            if idx_key != name_key {
                scalars.push((idx_key, Scalar::Num(ec_matrix_hz[i][j])));
            }
        }
    }

    let cap_data = CapacitanceMatrixData {
        terminals: terminals.into_iter().map(|(n, _)| n).collect(),
        maxwell_matrix_f,
        mutual_matrix_f,
        ec_matrix_hz,
        ..Default::default()
    };
    scalars.push(cap_data.scalar());
    scalars.push(tls.scalar());

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: mesh.dim as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
    })
}
