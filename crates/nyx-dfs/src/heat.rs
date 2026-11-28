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
use crate::problem::{BcKind, Problem, Study, SIGMA_SB};
use femcore::expr::eval_at;
use femcore::fem::*;
use femcore::fem_mixed::{
    assemble_mass_mixed, eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, lagrange_shape_into,
    quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::fields::{gradients_at_nodes, node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::driver::{backward_residual_many, dense_gen_eigh as gen_eigh, DirectSolver};
use femcore::linalg::snes::{newton, LineSearch, NewtonOptions};
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
    if let Some(ref cv) = spec.heat_source_cell_nodes {
        let f_imp = assemble_imported_cell_heat_source_simplex(mesh, order, cv)?;
        for (fi, &bi) in f.iter_mut().zip(&f_imp) {
            *fi += bi;
        }
    } else if let Some(ref hs) = spec.heat_source {
        let f_imp = assemble_imported_heat_source_simplex(mesh, order, hs)?;
        for (fi, &bi) in f.iter_mut().zip(&f_imp) {
            *fi += bi;
        }
    }
    let mut tfix: BTreeMap<i32, f64> = BTreeMap::new();
    let mut dirichlet_labels: BTreeSet<i32> = BTreeSet::new();
    let mut has_convection = false;
    let mut has_radiation = false;
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
                BcKind::Radiation { emissivity, view_factor, .. } => {
                    has_radiation |= *emissivity * *view_factor > 0.0;
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

    if fixed.is_empty() && !has_convection && !has_radiation {
        return Err("heat transfer needs a fixed-temperature, convection, or radiation boundary (an all-insulated body has no unique solution)".into());
    }
    Ok(HeatSystem { k, f, ktensor, fixed, vals, has_convection })
}

/// Assembles a domain-wise conservatively projected per-cell heat source `cv` (`nc * mesh.nloc(order)`)
/// into a consistent Lagrange load vector `f_i = sum_c int_{K_c} phi_i Q_c(x) dV` of length `mesh.n_dofs(order)`.
pub(crate) fn assemble_imported_cell_heat_source_simplex(
    mesh: &Discretisation,
    order: usize,
    cv: &[f64],
) -> Result<Vec<f64>, String> {
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(order);
    let n = mesh.n_dofs(order);
    if cv.len() != nc * nloc {
        return Err(format!(
            "projected per-cell 'heat_source' length {} does not match cells * nloc ({} * {} = {})",
            cv.len(),
            nc,
            nloc,
            nc * nloc
        ));
    }
    let re = RefElement::new(mesh.dim, order, (2 * order).max(2));
    let geo = geometry(mesh);
    let mut out = vec![0.0; n];
    let mut dofs = Vec::new();
    for c in 0..nc {
        mesh.cell_dofs(c, order, &mut dofs);
        let base = c * nloc;
        let vol = geo.vol[c];
        for q in 0..re.w.len() {
            let mut q_val = 0.0;
            for b in 0..nloc {
                q_val += re.n[q][b] * cv[base + b];
            }
            let wq = re.w[q] * vol * q_val;
            for a in 0..nloc {
                out[dofs[a]] += wq * re.n[q][a];
            }
        }
    }
    Ok(out)
}

/// Assembles an imported volumetric `heat_source` vector (`W/m^3`, cellwise or P1/P2 nodal) into a consistent
/// Lagrange load vector `f_i = int_Omega phi_i Q_imp(x) dV` of length `mesh.n_dofs(order)`.
pub(crate) fn assemble_imported_heat_source_simplex(
    mesh: &Discretisation,
    order: usize,
    hs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = mesh.n_dofs(order);
    let nv = mesh.n_points();
    let nc = mesh.n_cells();
    let n2 = mesh.n_dofs(2);
    if hs.len() == n {
        let m1 = assemble_mass(mesh, order, &vec![1.0; nc]);
        Ok(m1.matvec(hs))
    } else if hs.len() == nc {
        Ok(assemble_source(mesh, order, hs))
    } else if order == 2 && hs.len() == nv {
        let mut q2 = vec![0.0; n];
        q2[..nv].copy_from_slice(hs);
        for (e, &(a, b)) in mesh.edges().iter().enumerate() {
            q2[nv + e] = 0.5 * (hs[a as usize] + hs[b as usize]);
        }
        let m1 = assemble_mass(mesh, 2, &vec![1.0; nc]);
        Ok(m1.matvec(&q2))
    } else if order == 1 && hs.len() == n2 {
        let re1 = RefElement::new(mesh.dim, 1, 4);
        let qrule = quad_rule(mesh.dim, 4);
        let geo = geometry(mesh);
        let shapes2: Vec<Vec<f64>> = qrule.lam.iter().map(|lam| shape(lam, 2).0).collect();
        let mut out = vec![0.0; n];
        let mut d1 = Vec::new();
        let mut d2 = Vec::new();
        for c in 0..nc {
            mesh.cell_dofs(c, 1, &mut d1);
            mesh.cell_dofs(c, 2, &mut d2);
            for q in 0..re1.w.len() {
                let q_val: f64 = shapes2[q].iter().zip(&d2).map(|(&nb, &gd)| nb * hs[gd]).sum();
                let wq = re1.w[q] * geo.vol[c] * q_val;
                for a in 0..re1.nloc {
                    out[d1[a]] += wq * re1.n[q][a];
                }
            }
        }
        Ok(out)
    } else {
        Err(format!(
            "imported 'heat_source' length {} does not match mesh DOFs ({n}), vertices ({nv}), P2 DOFs ({n2}), or cells ({nc})",
            hs.len()
        ))
    }
}

fn assemble_imported_heat_source_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    hs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = dof_map.n_dofs;
    let nv = mixed.n_points();
    let nc = mixed.n_cells();
    let n_simp2 = simplex.n_dofs(2);
    if hs.len() == n {
        let m1 = assemble_mass_mixed(mixed, dof_map, &vec![1.0; nc]);
        Ok(m1.matvec(hs))
    } else if hs.len() == nc {
        Ok(assemble_source_mixed(mixed, dof_map, hs))
    } else if hs.len() == nv || hs.len() == n_simp2 {
        // Integrate simplex P1/P2 nodal heat_source on `simplex` to get total & vertex-prolongated load on `mixed`
        let order_s = if hs.len() == n_simp2 && n_simp2 > nv { 2 } else { 1 };
        let m1_s = assemble_mass(simplex, order_s, &vec![1.0; simplex.n_cells()]);
        let total_p: f64 = m1_s.matvec(hs).iter().sum();
        let mut qm = vec![0.0; n];
        qm[..nv].copy_from_slice(&hs[..nv]);
        for (e, &(a, b)) in mixed.edges.iter().enumerate() {
            if nv + e < n {
                if order_s == 2 {
                    if let Ok(se) = simplex.edges().binary_search(&(a.min(b), a.max(b))) {
                        qm[nv + e] = hs[nv + se];
                        continue;
                    }
                }
                qm[nv + e] = 0.5 * (hs[a as usize] + hs[b as usize]);
            }
        }
        let m1_m = assemble_mass_mixed(mixed, dof_map, &vec![1.0; nc]);
        let mut fm = m1_m.matvec(&qm);
        let sum_m: f64 = fm.iter().sum();
        if sum_m.abs() > 1e-30 && total_p.abs() > 0.0 {
            let sc = total_p / sum_m;
            fm.iter_mut().for_each(|v| *v *= sc);
        }
        Ok(fm)
    } else {
        Err(format!(
            "imported 'heat_source' length {} does not match mixed DOFs ({n}), vertices ({nv}), P2 simplex DOFs ({n_simp2}), or cells ({nc})",
            hs.len()
        ))
    }
}

fn heat_fields_with_spec(
    mesh: &Discretisation,
    spec: &Problem,
    t: &[f64],
    ktensor: &[M3],
    order: usize,
) -> (Vec<Field>, CellFields) {
    if spec.k_curves.is_empty() {
        return heat_fields(mesh, t, ktensor, order);
    }
    let d = mesh.dim;
    let nc = mesh.n_cells();
    let nloc = mesh.nloc(order);
    let mut g = gradients_at_nodes(mesh, t, 1, order);
    let mut pd = Vec::new();
    for c in 0..nc {
        mesh.cell_dofs(c, order, &mut pd);
        let dom = mesh.cell_labels[c];
        for a in 0..nloc {
            let base = (c * nloc + a) * d;
            let gv: Vec<f64> = g[base..base + d].to_vec();
            let kc = match spec.k_curves.get(&dom) {
                Some((curve, euler)) => femcore::section::plane_tensor(&curve.eval_matrix_global(t[pd[a]], *euler), d),
                None => ktensor[c],
            };
            for i in 0..d {
                g[base + i] = -(0..d).map(|j| kc[i][j] * gv[j]).sum::<f64>();
            }
        }
    }
    let q = node_average(mesh, &g, d, order);
    let qc = components(&q, d);
    let names = femcore::section::component_names(d);
    let fields = vec![Field::scalar("T", "K", t.to_vec()), Field::vector("q", "W/m^2", names, &qc)];
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("T", "K", cell_values_at_nodes(mesh, t, order));
    cf.push_components("q", names, "W/m^2", &g);
    (fields, cf)
}

fn radiation_facet_rule(facet_dim: usize) -> Quad {
    if facet_dim == 1 {
        // 5-point Gauss-Legendre rule on [0, 1] (exact to degree 9 for T^4 * phi_a)
        let s = 2.0 * (10.0_f64 / 7.0).sqrt();
        let x1 = (5.0 - s).sqrt() / 3.0;
        let x2 = (5.0 + s).sqrt() / 3.0;
        let w0 = 128.0 / 225.0;
        let w1 = (322.0 + 13.0 * 70.0_f64.sqrt()) / 900.0;
        let w2 = (322.0 - 13.0 * 70.0_f64.sqrt()) / 900.0;
        let xs = [-x2, -x1, 0.0, x1, x2];
        let ws = [w2, w1, w0, w1, w2];
        Quad {
            lam: xs.iter().map(|&xi| vec![0.5 * (1.0 - xi), 0.5 * (1.0 + xi)]).collect(),
            w: ws.iter().map(|&w| 0.5 * w).collect(),
        }
    } else {
        quad_rule(2, 4)
    }
}

#[derive(Clone)]
struct RadFacetSimplex {
    dofs: Vec<usize>,
    measure: f64,
    alpha_rad: f64,
    t_amb4: f64,
    t_ambient: f64,
}

fn collect_radiation_facets_simplex(mesh: &Discretisation, spec: &Problem) -> Vec<RadFacetSimplex> {
    let mut out = Vec::new();
    let mut dofs = Vec::new();
    for bc in &spec.bcs {
        if let BcKind::Radiation { emissivity, t_ambient, view_factor } = bc.kind {
            let alpha = emissivity * view_factor * SIGMA_SB;
            if alpha <= 0.0 {
                continue;
            }
            let t_amb4 = t_ambient.powi(4);
            for f in facets_with(mesh, &bc.labels) {
                mesh.facet_dofs(f, spec.order, &mut dofs);
                let m = facet_measure(mesh, f);
                if m > 0.0 {
                    out.push(RadFacetSimplex { dofs: dofs.clone(), measure: m, alpha_rad: alpha, t_amb4, t_ambient });
                }
            }
        }
    }
    out
}

fn assemble_linear_boundary_robin_simplex(mesh: &Discretisation, spec: &Problem, ktensor: &[M3]) -> Csr {
    let order = spec.order;
    let n = mesh.n_dofs(order);
    let mut k_bnd = Csr::zeros(n, n);
    let mut tfix: BTreeSet<i32> = BTreeSet::new();
    for bc in &spec.bcs {
        for &lab in &bc.labels {
            match &bc.kind {
                BcKind::Temperature(_) => {
                    tfix.insert(lab);
                }
                BcKind::Convection { h, .. } if *h != 0.0 => {
                    k_bnd = k_bnd.add_scaled(&assemble_boundary_mass(mesh, &[lab], order, *h), 1.0);
                }
                _ => {}
            }
        }
    }
    let ref_labels: Vec<i32> = tfix.iter().copied().collect();
    let k_perm = FarFieldPermittivity::Tensor(ktensor);
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !tfix.contains(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid(mesh, &ref_labels));
                let k_far = assemble_farfield_robin(mesh, &ff_labs, order, *ff_order, *decay_power, c, 1.0, k_perm);
                k_bnd = k_bnd.add_scaled(&k_far, 1.0);
            }
        }
    }
    if spec.has_infinite_element() {
        if let Ok(mixed) = MixedDiscretisation::from_simplex_mesh(mesh) {
            let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
            let dof_map = MixedLagrangeDofMap::new(&mixed, space);
            if let Ok(ie_map) = resolve_ie_domains_mixed(&mixed, spec) {
                let (_, ie_exp_facets) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
                if !ie_exp_facets.is_empty() {
                    let k_ie_rob =
                        assemble_ie_outer_robin_mixed(&mixed, &dof_map, &ie_exp_facets, &ie_map, 1.0, k_perm);
                    k_bnd = k_bnd.add_scaled(&k_ie_rob, 1.0);
                }
            }
        }
    }
    k_bnd
}

fn build_simplex_free_csr_template(
    mesh: &Discretisation,
    order: usize,
    n_free: usize,
    inv_free: &[usize],
    rad_facets: &[RadFacetSimplex],
    extra_ff: &Csr,
) -> (Vec<usize>, Vec<u32>) {
    let mut rows: Vec<Vec<u32>> = vec![Vec::new(); n_free];
    for r in 0..n_free {
        rows[r].push(r as u32);
    }
    for r in 0..extra_ff.nrows() {
        rows[r].extend_from_slice(extra_ff.row(r).0);
    }
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        for &ga in &dofs {
            let ra = inv_free[ga];
            if ra == usize::MAX {
                continue;
            }
            for &gb in &dofs {
                let rb = inv_free[gb];
                if rb != usize::MAX {
                    rows[ra].push(rb as u32);
                }
            }
        }
    }
    for rf in rad_facets {
        for &ga in &rf.dofs {
            let ra = inv_free[ga];
            if ra == usize::MAX {
                continue;
            }
            for &gb in &rf.dofs {
                let rb = inv_free[gb];
                if rb != usize::MAX {
                    rows[ra].push(rb as u32);
                }
            }
        }
    }
    let mut indptr = Vec::with_capacity(n_free + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    for row in &mut rows {
        row.sort_unstable();
        row.dedup();
        indices.extend_from_slice(row);
        indptr.push(indices.len());
    }
    (indptr, indices)
}

#[inline]
fn add_to_csr_template(indptr: &[usize], indices: &[u32], values: &mut [f64], row: usize, col: usize, val: f64) {
    if val == 0.0 {
        return;
    }
    let range = indptr[row]..indptr[row + 1];
    let col_u32 = col as u32;
    if let Ok(pos) = indices[range.clone()].binary_search(&col_u32) {
        values[range.start + pos] += val;
    }
}

fn copy_csr_into_template(src: &Csr, indptr: &[usize], indices: &[u32], values: &mut [f64]) {
    for r in 0..src.nrows() {
        let (cols, vals) = src.row(r);
        let range = indptr[r]..indptr[r + 1];
        for (&c, &v) in cols.iter().zip(vals) {
            if v != 0.0 {
                if let Ok(pos) = indices[range.clone()].binary_search(&c) {
                    values[range.start + pos] += v;
                }
            }
        }
    }
}

/// Solves nonlinear stationary and transient (BDF2) heat transfer via `femcore::linalg::snes::newton`
/// with exact analytical radiation (`+4 eps F sigma_SB T^3`) and material (`k'(T)`, `Cp'(T)`) tangents
/// and fixed symbolic factorization reuse across Newton steps.
fn solve_nonlinear_heat_states_simplex<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    sys: &HeatSystem,
    source: &S,
    opts: Option<TransientHeatOptions>,
) -> Result<(HeatStates, DfsDiag), String> {
    let dim = mesh.dim;
    let order = spec.order;
    let n = sys.k.nrows();
    let free = free_dofs(n, &sys.fixed);
    let n_free = free.len();
    let mut inv_free = vec![usize::MAX; n];
    for (r, &gd) in free.iter().enumerate() {
        inv_free[gd] = r;
    }

    let rad_facets = collect_radiation_facets_simplex(mesh, spec);
    let has_rad = !rad_facets.is_empty();
    let has_k_curve = !spec.k_curves.is_empty();
    let frule = radiation_facet_rule(dim - 1);
    let fshapes: Vec<Vec<f64>> = frule.lam.iter().map(|lam| shape(lam, order).0).collect();

    let deg = (2 * order).max(2);
    let re = RefElement::new(dim, order, deg);
    let geo = geometry(mesh);

    // Base linear matrix on free DOFs when k is constant, or linear Robin boundary matrix when k(T) varies
    let k_lin_full =
        if has_k_curve { assemble_linear_boundary_robin_simplex(mesh, spec, &sys.ktensor) } else { sys.k.clone() };
    let k_lin_ff = k_lin_full.submatrix(&free, &free);
    let k_lin_fx = k_lin_full.submatrix(&free, &sys.fixed);

    let (indptr, indices) = build_simplex_free_csr_template(mesh, order, n_free, &inv_free, &rad_facets, &k_lin_ff);

    // Helper evaluating stationary internal + radiation forces on `u_free` (given fixed `u_dir`)
    let eval_stat_res = |u_free: &[f64], u_dir: &[f64], f_load: &[f64], out: &mut [f64]| {
        if n_free == 0 {
            return;
        }
        let mut t_full = vec![0.0; n];
        for (&gd, &v) in sys.fixed.iter().zip(u_dir) {
            t_full[gd] = v;
        }
        for (r, &gd) in free.iter().enumerate() {
            t_full[gd] = u_free[r];
        }

        // Guard against non-physical negative absolute temperature on radiation facets during line search
        if has_rad {
            for rf in &rad_facets {
                for fshape in &fshapes {
                    let tq: f64 = fshape.iter().zip(&rf.dofs).map(|(&na, &gd)| na * t_full[gd]).sum();
                    if !tq.is_finite() || tq < 0.0 {
                        out.fill(f64::NAN);
                        return;
                    }
                }
            }
        }

        let kx = if sys.fixed.is_empty() { vec![0.0; n_free] } else { k_lin_fx.matvec(u_dir) };
        let ku = k_lin_ff.matvec(u_free);
        for r in 0..n_free {
            out[r] = ku[r] + kx[r] - f_load[free[r]];
        }

        if has_k_curve {
            let chunks: Vec<usize> = (0..mesh.n_cells()).collect();
            let parts = femcore::par::map_chunks(&chunks, 256, |slice| {
                let mut loc = vec![0.0; n_free];
                let mut dofs = Vec::new();
                let mut grad = vec![0.0; re.nloc * dim];
                for &c in slice {
                    mesh.cell_dofs(c, order, &mut dofs);
                    let dom = mesh.cell_labels[c];
                    for q in 0..re.w.len() {
                        shape_gradients(&re, &geo, c, q, &mut grad);
                        let mut tq = 0.0;
                        let mut gt = [0.0; 3];
                        for a in 0..re.nloc {
                            let ta = t_full[dofs[a]];
                            tq += ta * re.n[q][a];
                            for i in 0..dim {
                                gt[i] += ta * grad[a * dim + i];
                            }
                        }
                        let kc = match spec.k_curves.get(&dom) {
                            Some((curve, euler)) => {
                                femcore::section::plane_tensor(&curve.eval_matrix_global(tq, *euler), dim)
                            }
                            None => sys.ktensor[c],
                        };
                        let mut flux = [0.0; 3];
                        for i in 0..dim {
                            for j in 0..dim {
                                flux[i] += kc[i][j] * gt[j];
                            }
                        }
                        let wq = re.w[q] * geo.vol[c];
                        for a in 0..re.nloc {
                            let ra = inv_free[dofs[a]];
                            if ra != usize::MAX {
                                let mut dot = 0.0;
                                for i in 0..dim {
                                    dot += grad[a * dim + i] * flux[i];
                                }
                                loc[ra] += wq * dot;
                            }
                        }
                    }
                }
                loc
            });
            for part in parts {
                for r in 0..n_free {
                    out[r] += part[r];
                }
            }
        }

        for rf in &rad_facets {
            for (q, fshape) in fshapes.iter().enumerate() {
                let tq: f64 = fshape.iter().zip(&rf.dofs).map(|(&na, &gd)| na * t_full[gd]).sum();
                let q_rad = rf.alpha_rad * (tq.powi(4) - rf.t_amb4);
                let wq = frule.w[q] * rf.measure * q_rad;
                for (a, &gd) in rf.dofs.iter().enumerate() {
                    let ra = inv_free[gd];
                    if ra != usize::MAX {
                        out[ra] += wq * fshape[a];
                    }
                }
            }
        }
    };

    let eval_stat_jac_values = |u_free: &[f64], u_dir: &[f64], values: &mut [f64]| {
        values.fill(0.0);
        copy_csr_into_template(&k_lin_ff, &indptr, &indices, values);

        let mut t_full = vec![0.0; n];
        for (&gd, &v) in sys.fixed.iter().zip(u_dir) {
            t_full[gd] = v;
        }
        for (r, &gd) in free.iter().enumerate() {
            t_full[gd] = u_free[r];
        }

        if has_k_curve {
            let nloc = re.nloc;
            let chunks: Vec<usize> = (0..mesh.n_cells()).collect();
            let parts = femcore::par::map_chunks(&chunks, 256, |slice| {
                let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(slice.len() * nloc * nloc);
                let mut dofs = Vec::new();
                let mut grad = vec![0.0; nloc * dim];
                let mut ke = vec![0.0; nloc * nloc];
                for &c in slice {
                    ke.fill(0.0);
                    mesh.cell_dofs(c, order, &mut dofs);
                    let dom = mesh.cell_labels[c];
                    for q in 0..re.w.len() {
                        shape_gradients(&re, &geo, c, q, &mut grad);
                        let mut tq = 0.0;
                        let mut gt = [0.0; 3];
                        for a in 0..nloc {
                            let ta = t_full[dofs[a]];
                            tq += ta * re.n[q][a];
                            for i in 0..dim {
                                gt[i] += ta * grad[a * dim + i];
                            }
                        }
                        let (kc, dkc) = match spec.k_curves.get(&dom) {
                            Some((curve, euler)) => (
                                femcore::section::plane_tensor(&curve.eval_matrix_global(tq, *euler), dim),
                                femcore::section::plane_tensor(&curve.eval_matrix_deriv_global(tq, *euler), dim),
                            ),
                            None => (sys.ktensor[c], [[0.0; 3]; 3]),
                        };
                        let mut dflux = [0.0; 3];
                        for i in 0..dim {
                            for j in 0..dim {
                                dflux[i] += dkc[i][j] * gt[j];
                            }
                        }
                        let wq = re.w[q] * geo.vol[c];
                        for a in 0..nloc {
                            let mut kg_a = [0.0; 3];
                            let mut ga_dot_dflux = 0.0;
                            for i in 0..dim {
                                ga_dot_dflux += grad[a * dim + i] * dflux[i];
                                for j in 0..dim {
                                    kg_a[i] += kc[j][i] * grad[a * dim + j];
                                }
                            }
                            for b in 0..nloc {
                                let mut diff_ab = 0.0;
                                for i in 0..dim {
                                    diff_ab += kg_a[i] * grad[b * dim + i];
                                }
                                ke[a * nloc + b] += wq * (diff_ab + ga_dot_dflux * re.n[q][b]);
                            }
                        }
                    }
                    for a in 0..nloc {
                        let ra = inv_free[dofs[a]];
                        if ra == usize::MAX {
                            continue;
                        }
                        for b in 0..nloc {
                            let rb = inv_free[dofs[b]];
                            if rb != usize::MAX {
                                triplets.push((ra, rb, ke[a * nloc + b]));
                            }
                        }
                    }
                }
                triplets
            });
            for part in parts {
                for (ra, rb, v) in part {
                    add_to_csr_template(&indptr, &indices, values, ra, rb, v);
                }
            }
        }

        for rf in &rad_facets {
            let nfloc = rf.dofs.len();
            for (q, fshape) in fshapes.iter().enumerate() {
                let tq: f64 = fshape.iter().zip(&rf.dofs).map(|(&na, &gd)| na * t_full[gd]).sum();
                let h_tan = 4.0 * rf.alpha_rad * tq.max(0.0).powi(3);
                let wq = frule.w[q] * rf.measure * h_tan;
                for a in 0..nfloc {
                    let ra = inv_free[rf.dofs[a]];
                    if ra == usize::MAX {
                        continue;
                    }
                    for b in 0..nfloc {
                        let rb = inv_free[rf.dofs[b]];
                        if rb != usize::MAX {
                            add_to_csr_template(&indptr, &indices, values, ra, rb, wq * fshape[a] * fshape[b]);
                        }
                    }
                }
            }
        }
    };

    let mut log = Vec::new();
    match &spec.study {
        Study::HeatStationary => {
            if n_free == 0 {
                let mut t = vec![0.0; n];
                for (&gd, &v) in sys.fixed.iter().zip(&sys.vals) {
                    t[gd] = v;
                }
                return Ok((HeatStates { times: vec![], temps: vec![t], log }, DfsDiag::new("snes_newton")));
            }
            // Linearized positive-definite warm start `u0`
            let mut ref_sum = 0.0;
            let mut ref_cnt = 0usize;
            for &v in &sys.vals {
                if v > 0.0 {
                    ref_sum += v;
                    ref_cnt += 1;
                }
            }
            for rf in &rad_facets {
                if rf.t_ambient > 0.0 {
                    ref_sum += rf.t_ambient;
                    ref_cnt += 1;
                }
            }
            let t_ref_init = if ref_cnt > 0 { (ref_sum / ref_cnt as f64).max(4.0) } else { 293.15 };

            let mut k_init = sys.k.clone();
            let mut f_init = sys.f.clone();
            for bc in &spec.bcs {
                if let BcKind::Radiation { emissivity, t_ambient, view_factor } = bc.kind {
                    let alpha = emissivity * view_factor * SIGMA_SB;
                    if alpha > 0.0 {
                        let h0 = 4.0 * alpha * t_ref_init.powi(3);
                        let q0 = alpha * (3.0 * t_ref_init.powi(4) + t_ambient.powi(4));
                        k_init = k_init.add_scaled(&assemble_boundary_mass(mesh, &bc.labels, order, h0), 1.0);
                        let l0 = assemble_boundary_load(mesh, &bc.labels, q0, order);
                        for (fi, &li) in f_init.iter_mut().zip(&l0) {
                            *fi += li;
                        }
                    }
                }
            }
            let k_init_ff = k_init.submatrix(&free, &free);
            let k_init_fx = k_init.submatrix(&free, &sys.fixed);
            let kx0 = if sys.fixed.is_empty() { vec![0.0; n_free] } else { k_init_fx.matvec(&sys.vals) };
            let rhs0: Vec<f64> = free.iter().enumerate().map(|(r, &gd)| f_init[gd] - kx0[r]).collect();
            let mut u0 = DirectSolver::new(&k_init_ff, true)
                .map(|s| s.solve(&rhs0))
                .unwrap_or_else(|_| vec![t_ref_init; n_free]);
            if has_rad {
                for v in &mut u0 {
                    if !v.is_finite() || *v <= 0.0 {
                        *v = t_ref_init.max(1.0);
                    }
                }
            }

            let rhs_scale = rhs0.iter().map(|x| x * x).sum::<f64>().sqrt().max(1.0);
            let diag_ku_scale: f64 =
                k_init_ff.diagonal().iter().zip(&u0).map(|(d, u)| (d * u).powi(2)).sum::<f64>().sqrt();
            let res_scale = rhs_scale.max(diag_ku_scale);
            let user_tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
            let n_opts = NewtonOptions {
                rtol: user_tol.min(1e-10),
                atol: (user_tol * 1e-2).clamp(1e-13, 1e-10) * res_scale,
                stol: 1e-14,
                max_it: 50,
                line_search: LineSearch::Backtracking { c: 1e-4, alpha_min: 1e-8 },
                definiteness: if has_k_curve { Definiteness::General } else { Definiteness::PositiveDefinite },
                jacobian_lag: 1,
            };
            let mut f_fn = |u: &[f64], fx: &mut [f64]| eval_stat_res(u, &sys.vals, &sys.f, fx);
            let mut jac_fn = |u: &[f64]| -> Csr {
                let mut values = vec![0.0; indptr[n_free]];
                eval_stat_jac_values(u, &sys.vals, &mut values);
                Csr::try_new(n_free, n_free, indptr.clone(), indices.clone(), values)
                    .expect("fixed-pattern nonlinear thermal Jacobian CSR")
            };
            let nl = newton(&mut f_fn, &mut jac_fn, &u0, &n_opts)
                .map_err(|e| format!("nonlinear thermal Newton linear solve failed: {e}"))?;
            let last_norm = nl.residual_norms.last().copied().unwrap_or(f64::NAN);
            if !nl.reason.is_converged() {
                return Err(format!(
                    "nonlinear thermal Newton solver failed to converge ({:?} after {} iterations, residual {:.3e})",
                    nl.reason, nl.iterations, last_norm
                ));
            }

            let mut t_sol = vec![0.0; n];
            for (&gd, &v) in sys.fixed.iter().zip(&sys.vals) {
                t_sol[gd] = v;
            }
            for (r, &gd) in free.iter().enumerate() {
                t_sol[gd] = nl.x[r];
            }

            let mut diag = DfsDiag::new("snes_newton");
            diag.iterations = Some(nl.iterations);
            diag.converged = true;
            diag.max_algebraic_residual = last_norm / res_scale;
            log.push(format!(
                "nonlinear heat Newton: {} iteration(s), {} symbolic factorisation(s), residual {:.2e}",
                nl.iterations, nl.symbolic_analyses, diag.max_algebraic_residual
            ));
            Ok((HeatStates { times: vec![], temps: vec![t_sol], log }, diag))
        }
        Study::HeatTransient { t_end, n_out, t0: t0_expr } => {
            check_transient_heat_guards(spec)?;
            if *t_end <= 0.0 || *n_out < 1 {
                return Err("end time must be positive and at least one output time is needed".into());
            }
            let rho = spec.materials.cell_values(mesh, "rho", Some(1.0))?;
            let cp_const = spec.materials.cell_values(mesh, "Cp", Some(1.0))?;
            let has_cp_curve = !spec.cp_curves.is_empty();
            let rcp_const: Vec<f64> = rho.iter().zip(&cp_const).map(|(a, b)| a * b).collect();
            let m_const = (!has_cp_curve).then(|| assemble_mass(mesh, order, &rcp_const));
            let m_const_ff = m_const.as_ref().map(|m| m.submatrix(&free, &free));
            let m_const_fx = m_const.as_ref().map(|m| m.submatrix(&free, &sys.fixed));

            let times: Vec<f64> = (0..=*n_out).map(|j| t_end * j as f64 / *n_out as f64).collect();
            let sub = opts.map_or(8, |o| o.bdf2_substeps).max(1);

            let mut dirichlet_labels = BTreeSet::new();
            for bc in &spec.bcs {
                if matches!(bc.kind, BcKind::Temperature(_)) {
                    dirichlet_labels.extend(bc.labels.iter().copied());
                }
            }
            let coords = mesh.dof_coordinates(order);
            let mut t0_vec = eval_at(t0_expr, &coords, dim).map_err(|e| format!("initial temperature: {e}"))?;
            let ud_prev = if source.is_uniform() {
                sys.vals.clone()
            } else if dirichlet_labels.is_empty() {
                vec![]
            } else {
                dirichlet_constraints_field::<D, S>(mesh, &dirichlet_labels, order, source, 0.0).1
            };
            for (&gd, &v) in sys.fixed.iter().zip(&ud_prev) {
                t0_vec[gd] = v;
            }

            let mut out = vec![vec![0.0; n]; times.len()];
            out[0].copy_from_slice(&t0_vec);
            if n_free == 0 {
                for j in 1..times.len() {
                    let ud_j = if source.is_uniform() {
                        sys.vals.clone()
                    } else {
                        dirichlet_constraints_field::<D, S>(mesh, &dirichlet_labels, order, source, times[j]).1
                    };
                    for (&gd, &v) in sys.fixed.iter().zip(&ud_j) {
                        out[j][gd] = v;
                    }
                }
                return Ok((HeatStates { times, temps: out, log }, DfsDiag::new("snes_newton_bdf2")));
            }

            let def = if has_k_curve || has_cp_curve { Definiteness::General } else { Definiteness::PositiveDefinite };

            let run_bdf_chain = |sub_steps: usize,
                                 mut record_out: Option<&mut Vec<Vec<f64>>>|
             -> Result<(Vec<f64>, usize, f64), String> {
                let dt_s = times[1] / sub_steps as f64;
                let mut u_prev: Vec<f64> = free.iter().map(|&gd| t0_vec[gd]).collect();
                let mut u_cur = u_prev.clone();
                let mut d_prev = ud_prev.clone();
                let mut d_cur = d_prev.clone();
                let mut chain_iters = 0usize;
                let mut chain_max_rel = 0.0_f64;

                for step in 1..=sub_steps * n_out {
                    femcore::progress::check()?;
                    let t_time = step as f64 * dt_s;
                    let (f_step, d_next) = if source.is_uniform() {
                        (sys.f.clone(), sys.vals.clone())
                    } else {
                        let (mut fs, uds) = assemble_heat_forcing_and_dirichlet_at::<D, S>(
                            mesh,
                            spec,
                            source,
                            t_time,
                            &dirichlet_labels,
                        )?;
                        if let Some(ref cv) = spec.heat_source_cell_nodes {
                            let f_imp = assemble_imported_cell_heat_source_simplex(mesh, order, cv)?;
                            for (fi, &bi) in fs.iter_mut().zip(&f_imp) {
                                *fi += bi;
                            }
                        } else if let Some(ref hs) = spec.heat_source {
                            let f_imp = assemble_imported_heat_source_simplex(mesh, order, hs)?;
                            for (fi, &bi) in fs.iter_mut().zip(&f_imp) {
                                *fi += bi;
                            }
                        }
                        (fs, uds)
                    };
                    // BDF1 on step 1: alpha_0 = 1.0, u_hist = u_cur; BDF2 on step >= 2: alpha_0 = 1.5, u_hist = 2 u_cur - 0.5 u_prev
                    let alpha_0 = if step == 1 { 1.0 } else { 1.5 };
                    let mut t_hist_full = vec![0.0; n];
                    if step == 1 {
                        for (&gd, &v) in sys.fixed.iter().zip(&d_cur) {
                            t_hist_full[gd] = v;
                        }
                        for (r, &gd) in free.iter().enumerate() {
                            t_hist_full[gd] = u_cur[r];
                        }
                    } else {
                        for (k_i, &gd) in sys.fixed.iter().enumerate() {
                            t_hist_full[gd] = 2.0 * d_cur[k_i] - 0.5 * d_prev[k_i];
                        }
                        for (r, &gd) in free.iter().enumerate() {
                            t_hist_full[gd] = 2.0 * u_cur[r] - 0.5 * u_prev[r];
                        }
                    }
                    let u_hist_free: Vec<f64> = free.iter().map(|&gd| t_hist_full[gd]).collect();
                    let d_bdf: Vec<f64> = (0..sys.fixed.len())
                        .map(|k_i| (alpha_0 * d_next[k_i] - t_hist_full[sys.fixed[k_i]]) / dt_s)
                        .collect();

                    let mut step_res = |u: &[f64], fx: &mut [f64]| {
                        eval_stat_res(u, &d_next, &f_step, fx);
                        if fx.first().is_some_and(|v| !v.is_finite()) {
                            return;
                        }
                        if let (Some(mff), Some(mfx)) = (&m_const_ff, &m_const_fx) {
                            let du: Vec<f64> =
                                u.iter().zip(&u_hist_free).map(|(&ui, &uhi)| (alpha_0 * ui - uhi) / dt_s).collect();
                            let m_du = mff.matvec(&du);
                            let m_dd = if sys.fixed.is_empty() { vec![0.0; n_free] } else { mfx.matvec(&d_bdf) };
                            for r in 0..n_free {
                                fx[r] += m_du[r] + m_dd[r];
                            }
                        } else {
                            let mut t_full = vec![0.0; n];
                            for (&gd, &v) in sys.fixed.iter().zip(&d_next) {
                                t_full[gd] = v;
                            }
                            for (r, &gd) in free.iter().enumerate() {
                                t_full[gd] = u[r];
                            }
                            let chunks: Vec<usize> = (0..mesh.n_cells()).collect();
                            let parts = femcore::par::map_chunks(&chunks, 256, |slice| {
                                let mut loc = vec![0.0; n_free];
                                let mut dofs = Vec::new();
                                for &c in slice {
                                    mesh.cell_dofs(c, order, &mut dofs);
                                    let dom = mesh.cell_labels[c];
                                    let rho_c = rho[c];
                                    for q in 0..re.w.len() {
                                        let mut tq = 0.0;
                                        let mut thq = 0.0;
                                        for a in 0..re.nloc {
                                            let na = re.n[q][a];
                                            tq += t_full[dofs[a]] * na;
                                            thq += t_hist_full[dofs[a]] * na;
                                        }
                                        let cp_q = match spec.cp_curves.get(&dom) {
                                            Some(curve) => curve.eval(tq),
                                            None => cp_const[c],
                                        };
                                        let rate = (alpha_0 * tq - thq) / dt_s;
                                        let wq = re.w[q] * geo.vol[c] * rho_c * cp_q * rate;
                                        for a in 0..re.nloc {
                                            let ra = inv_free[dofs[a]];
                                            if ra != usize::MAX {
                                                loc[ra] += wq * re.n[q][a];
                                            }
                                        }
                                    }
                                }
                                loc
                            });
                            for part in parts {
                                for r in 0..n_free {
                                    fx[r] += part[r];
                                }
                            }
                        }
                    };

                    let mut step_jac = |u: &[f64]| -> Csr {
                        let mut values = vec![0.0; indptr[n_free]];
                        eval_stat_jac_values(u, &d_next, &mut values);
                        if let Some(mff) = &m_const_ff {
                            let scale_m = alpha_0 / dt_s;
                            for r in 0..mff.nrows() {
                                let (cols, vals) = mff.row(r);
                                let range = indptr[r]..indptr[r + 1];
                                for (&c_col, &v) in cols.iter().zip(vals) {
                                    if v != 0.0 {
                                        if let Ok(pos) = indices[range.clone()].binary_search(&c_col) {
                                            values[range.start + pos] += scale_m * v;
                                        }
                                    }
                                }
                            }
                        } else {
                            let mut t_full = vec![0.0; n];
                            for (&gd, &v) in sys.fixed.iter().zip(&d_next) {
                                t_full[gd] = v;
                            }
                            for (r, &gd) in free.iter().enumerate() {
                                t_full[gd] = u[r];
                            }
                            let nloc = re.nloc;
                            let chunks: Vec<usize> = (0..mesh.n_cells()).collect();
                            let parts = femcore::par::map_chunks(&chunks, 256, |slice| {
                                let mut triplets = Vec::with_capacity(slice.len() * nloc * nloc);
                                let mut dofs = Vec::new();
                                let mut me = vec![0.0; nloc * nloc];
                                for &c in slice {
                                    me.fill(0.0);
                                    mesh.cell_dofs(c, order, &mut dofs);
                                    let dom = mesh.cell_labels[c];
                                    let rho_c = rho[c];
                                    for q in 0..re.w.len() {
                                        let mut tq = 0.0;
                                        let mut thq = 0.0;
                                        for a in 0..nloc {
                                            let na = re.n[q][a];
                                            tq += t_full[dofs[a]] * na;
                                            thq += t_hist_full[dofs[a]] * na;
                                        }
                                        let (cp_q, dcp_q) = match spec.cp_curves.get(&dom) {
                                            Some(curve) => (curve.eval(tq), curve.eval_deriv(tq)),
                                            None => (cp_const[c], 0.0),
                                        };
                                        let coeff = (rho_c / dt_s) * (alpha_0 * cp_q + dcp_q * (alpha_0 * tq - thq));
                                        let wq = re.w[q] * geo.vol[c] * coeff;
                                        for a in 0..nloc {
                                            let na = re.n[q][a];
                                            for b in 0..nloc {
                                                me[a * nloc + b] += wq * na * re.n[q][b];
                                            }
                                        }
                                    }
                                    for a in 0..nloc {
                                        let ra = inv_free[dofs[a]];
                                        if ra == usize::MAX {
                                            continue;
                                        }
                                        for b in 0..nloc {
                                            let rb = inv_free[dofs[b]];
                                            if rb != usize::MAX {
                                                triplets.push((ra, rb, me[a * nloc + b]));
                                            }
                                        }
                                    }
                                }
                                triplets
                            });
                            for part in parts {
                                for (ra, rb, v) in part {
                                    add_to_csr_template(&indptr, &indices, &mut values, ra, rb, v);
                                }
                            }
                        }
                        Csr::try_new(n_free, n_free, indptr.clone(), indices.clone(), values)
                            .expect("fixed-pattern nonlinear transient Jacobian CSR")
                    };

                    let ref_norm = u_cur.iter().map(|x| x * x).sum::<f64>().sqrt().max(1.0);
                    let diag_ku_scale: f64 =
                        k_lin_ff.diagonal().iter().zip(&u_cur).map(|(d, u)| (d * u).powi(2)).sum::<f64>().sqrt();
                    let res_scale = (ref_norm / dt_s.max(1e-12)).max(diag_ku_scale);
                    let n_opts = NewtonOptions {
                        rtol: 1e-11,
                        atol: 1e-12 * res_scale,
                        stol: 1e-14,
                        max_it: 30,
                        line_search: LineSearch::Backtracking { c: 1e-4, alpha_min: 1e-8 },
                        definiteness: def,
                        jacobian_lag: 1,
                    };
                    let nl = newton(&mut step_res, &mut step_jac, &u_cur, &n_opts)
                        .map_err(|e| format!("nonlinear transient Newton step {step} failed: {e}"))?;
                    let last_norm = nl.residual_norms.last().copied().unwrap_or(f64::NAN);
                    if !nl.reason.is_converged() {
                        return Err(format!(
                            "nonlinear transient Newton step {step} failed to converge ({:?} after {} iterations, residual {:.3e})",
                            nl.reason, nl.iterations, last_norm
                        ));
                    }
                    chain_iters += nl.iterations;
                    let f0 = nl.residual_norms.first().copied().unwrap_or(1.0).max(res_scale);
                    chain_max_rel = chain_max_rel.max(last_norm / f0);
                    u_prev = std::mem::replace(&mut u_cur, nl.x);
                    d_prev = std::mem::replace(&mut d_cur, d_next);
                    if let Some(out_rows) = record_out.as_deref_mut() {
                        if step % sub_steps == 0 {
                            let idx = step / sub_steps;
                            for (&gd, &v) in sys.fixed.iter().zip(&d_cur) {
                                out_rows[idx][gd] = v;
                            }
                            for (r, &gd) in free.iter().enumerate() {
                                out_rows[idx][gd] = u_cur[r];
                            }
                        }
                    }
                }
                Ok((u_cur, chain_iters, chain_max_rel))
            };

            let (u_end, total_newton_iters, max_rel_res) = run_bdf_chain(sub, Some(&mut out))?;
            let mut diag = DfsDiag::new("snes_newton_bdf2");
            diag.iterations = Some(total_newton_iters);
            diag.converged = true;
            diag.max_algebraic_residual = max_rel_res;
            if n_free <= 500 {
                if let Ok((u_half, _, _)) = run_bdf_chain(2 * sub, None) {
                    let diff_sq: f64 = u_half.iter().zip(&u_end).map(|(a, b)| (a - b).powi(2)).sum();
                    let norm_h = u_half.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
                    diag.set_temporal_rel(diff_sq.sqrt() / (3.0 * norm_h));
                }
            }
            log.push(format!(
                "nonlinear BDF2, {} steps ({} free dofs, {} total Newton iterations)",
                sub * n_out,
                n_free,
                total_newton_iters
            ));
            Ok((HeatStates { times, temps: out, log }, diag))
        }
        Study::Electrostatic | Study::Conduction | Study::Magnetostatic { .. } => {
            Err("heat transfer supports the stationary and transient studies".into())
        }
    }
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
    if spec.is_nonlinear_heat() {
        return solve_nonlinear_heat_states_simplex::<D, S>(mesh, spec, sys, source, opts);
    }
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
        Study::Electrostatic | Study::Conduction | Study::Magnetostatic { .. } => {
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
        let (fields, cell_fields) = heat_fields_with_spec(mesh, spec, row, &sys.ktensor, order);
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
        let (fields, cell_fields) = heat_fields_with_spec(mesh, spec, row, &sys.ktensor, order);
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
    if let Some(ref cv) = spec.heat_source_cell_nodes {
        let f_s = assemble_imported_cell_heat_source_simplex(&simplex, order, cv)?;
        let f_imp = assemble_imported_heat_source_mixed(&mixed, &dof_map, &simplex, &{
            let m1_s = assemble_mass(&simplex, order, &vec![1.0; simplex.n_cells()]);
            let solver_s = DirectSolver::new(&m1_s, true)?;
            solver_s.solve(&f_s)
        })?;
        for (fi, &bi) in f.iter_mut().zip(&f_imp) {
            *fi += bi;
        }
    } else if let Some(ref hs) = spec.heat_source {
        let f_imp = assemble_imported_heat_source_mixed(&mixed, &dof_map, &simplex, hs)?;
        for (fi, &bi) in f.iter_mut().zip(&f_imp) {
            *fi += bi;
        }
    }

    let mut tfix: BTreeMap<i32, f64> = BTreeMap::new();
    let mut has_convection = false;
    let mut has_radiation = false;
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
                BcKind::Radiation { emissivity, view_factor, .. } => {
                    has_radiation |= *emissivity * *view_factor > 0.0;
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

    if fixed.is_empty() && !has_convection && !has_radiation {
        return Err("heat transfer needs a fixed-temperature, convection, or radiation boundary (an all-insulated body has no unique solution)".into());
    }

    if spec.is_nonlinear_heat() {
        return solve_nonlinear_heat_mixed_impl(
            t0, &mixed, &dof_map, &simplex, &proj, &ie_map, spec, &k, &f, &ktensor, &fixed, &vals,
        );
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
        Study::Electrostatic | Study::Conduction | Study::Magnetostatic { .. } => {
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

#[allow(clippy::too_many_arguments)]
fn solve_nonlinear_heat_mixed_impl(
    t0: Instant,
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    ie_map: &BTreeMap<i32, crate::common::ResolvedIeDomain>,
    spec: &Problem,
    k_base: &Csr,
    f: &[f64],
    ktensor: &[M3],
    fixed: &[usize],
    vals: &[f64],
) -> Result<SolveResult, String> {
    let dim = mixed.dim;
    let order = spec.order;
    let space = dof_map.space;
    let n = dof_map.n_dofs;
    let free = free_dofs(n, fixed);
    let n_free = free.len();
    let mut inv_free = vec![usize::MAX; n];
    for (r, &gd) in free.iter().enumerate() {
        inv_free[gd] = r;
    }

    // Collect mixed radiation facets
    let mut rad_facets: Vec<(usize, f64, f64, f64)> = Vec::new(); // (facet_idx, alpha_rad, t_amb4, t_amb)
    for bc in &spec.bcs {
        if let BcKind::Radiation { emissivity, t_ambient, view_factor } = bc.kind {
            let alpha = emissivity * view_factor * SIGMA_SB;
            if alpha <= 0.0 {
                continue;
            }
            let t_amb4 = t_ambient.powi(4);
            for fi in 0..mixed.n_facets() {
                if bc.labels.contains(&mixed.facet_labels[fi]) {
                    rad_facets.push((fi, alpha, t_amb4, t_ambient));
                }
            }
        }
    }
    let has_rad = !rad_facets.is_empty();
    let has_k_curve = !spec.k_curves.is_empty();

    // Linear boundary matrix when k(T) varies, else full k_base
    let k_lin_full = if has_k_curve {
        let zeros = vec![0.0; mixed.n_cells()];
        let k0_vol = assemble_laplace_ie_mixed(mixed, dof_map, &zeros, None, ie_map);
        let mut kb = k0_vol;
        let mut tfix_set = BTreeSet::new();
        for bc in &spec.bcs {
            for &lab in &bc.labels {
                match &bc.kind {
                    BcKind::Temperature(_) => {
                        tfix_set.insert(lab);
                    }
                    BcKind::Convection { h, .. } if *h != 0.0 => {
                        kb = kb.add_scaled(&assemble_boundary_mass_mixed(mixed, dof_map, &[lab], *h), 1.0);
                    }
                    _ => {}
                }
            }
        }
        let ref_labels: Vec<i32> = tfix_set.iter().copied().collect();
        let k_perm = FarFieldPermittivity::Tensor(ktensor);
        for bc in &spec.bcs {
            if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
                let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !tfix_set.contains(l)).collect();
                if !ff_labs.is_empty() {
                    let c = center.unwrap_or_else(|| farfield_centroid_mixed(mixed, &ref_labels));
                    let k_far = assemble_farfield_robin_mixed(
                        mixed,
                        dof_map,
                        &ff_labs,
                        *ff_order,
                        *decay_power,
                        c,
                        1.0,
                        k_perm,
                    );
                    kb = kb.add_scaled(&k_far, 1.0);
                }
            }
        }
        let (_, ie_exp_facets) = classify_ie_outer_facets_mixed(mixed, ie_map);
        if !ie_exp_facets.is_empty() {
            let k_ie_rob = assemble_ie_outer_robin_mixed(mixed, dof_map, &ie_exp_facets, ie_map, 1.0, k_perm);
            kb = kb.add_scaled(&k_ie_rob, 1.0);
        }
        kb
    } else {
        k_base.clone()
    };
    let k_lin_ff = k_lin_full.submatrix(&free, &free);
    let k_lin_fx = k_lin_full.submatrix(&free, fixed);

    // Build fixed CSR sparsity pattern on free DOFs
    let mut rows: Vec<Vec<u32>> = vec![Vec::new(); n_free];
    for r in 0..n_free {
        rows[r].push(r as u32);
        rows[r].extend_from_slice(k_lin_ff.row(r).0);
    }
    for c in 0..mixed.n_cells() {
        let cdofs = dof_map.cell_dofs(c);
        for &ga in cdofs {
            let ra = inv_free[ga];
            if ra == usize::MAX {
                continue;
            }
            for &gb in cdofs {
                let rb = inv_free[gb];
                if rb != usize::MAX {
                    rows[ra].push(rb as u32);
                }
            }
        }
    }
    let mut fbuf = Vec::new();
    for &(fi, _, _, _) in &rad_facets {
        dof_map.facet_dofs(mixed, fi, &mut fbuf);
        for &ga in &fbuf {
            let ra = inv_free[ga];
            if ra == usize::MAX {
                continue;
            }
            for &gb in &fbuf {
                let rb = inv_free[gb];
                if rb != usize::MAX {
                    rows[ra].push(rb as u32);
                }
            }
        }
    }
    let mut indptr = Vec::with_capacity(n_free + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    for row in &mut rows {
        row.sort_unstable();
        row.dedup();
        indices.extend_from_slice(row);
        indptr.push(indices.len());
    }

    let eval_mixed_stat_res = |u_free: &[f64], out: &mut [f64]| {
        if n_free == 0 {
            return;
        }
        let mut t_full = vec![0.0; n];
        for (&gd, &v) in fixed.iter().zip(vals) {
            t_full[gd] = v;
        }
        for (r, &gd) in free.iter().enumerate() {
            t_full[gd] = u_free[r];
        }
        let kx = if fixed.is_empty() { vec![0.0; n_free] } else { k_lin_fx.matvec(vals) };
        let ku = k_lin_ff.matvec(u_free);
        for r in 0..n_free {
            out[r] = ku[r] + kx[r] - f[free[r]];
        }
        if has_k_curve {
            let mut vals_s = Vec::new();
            let mut dxi = Vec::new();
            let mut sn = Vec::new();
            let mut sdn = Vec::new();
            for c in 0..mixed.n_cells() {
                let ct = mixed.cell_type(c);
                let iso = IsoparametricCell::from_mesh(mixed, c);
                let rule = quad_rule_mixed(ct, (2 * order).max(2));
                let cdofs = dof_map.cell_dofs(c);
                let dom = mixed.cell_labels[c];
                for (xi, &wq_ref) in rule.pts.iter().zip(&rule.w) {
                    let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                    lagrange_shape_into(ct, space, xi, &mut vals_s, &mut dxi);
                    let grads: Vec<[f64; 3]> = dxi.iter().map(|&g_ref| jac.map_grad(g_ref)).collect();
                    let mut tq = 0.0;
                    let mut gt = [0.0; 3];
                    for (a, &gd) in cdofs.iter().enumerate() {
                        let ta = t_full[gd];
                        tq += ta * vals_s[a];
                        for i in 0..dim {
                            gt[i] += ta * grads[a][i];
                        }
                    }
                    let kc = match spec.k_curves.get(&dom) {
                        Some((curve, euler)) => {
                            femcore::section::plane_tensor(&curve.eval_matrix_global(tq, *euler), dim)
                        }
                        None => ktensor[c],
                    };
                    let mut flux = [0.0; 3];
                    for i in 0..dim {
                        for j in 0..dim {
                            flux[i] += kc[i][j] * gt[j];
                        }
                    }
                    let wq = wq_ref * jac.det.abs();
                    for (a, &gd) in cdofs.iter().enumerate() {
                        let ra = inv_free[gd];
                        if ra != usize::MAX {
                            let mut dot = 0.0;
                            for i in 0..dim {
                                dot += grads[a][i] * flux[i];
                            }
                            out[ra] += wq * dot;
                        }
                    }
                }
            }
        }
        let fdeg = 2 * order + 4;
        let mut fdofs = Vec::new();
        for &(fi, alpha, t_amb4, _) in &rad_facets {
            let ft = mixed.facet_type(fi);
            let frule = facet_rule_mixed(ft, fdeg);
            dof_map.facet_dofs(mixed, fi, &mut fdofs);
            for (xi, &wq_ref) in frule.pts.iter().zip(&frule.w) {
                let (_, ds, _) = eval_facet_jacobian(mixed, fi, xi);
                let (nv, _) = facet_lagrange_shape(ft, space, xi);
                let tq: f64 = nv.iter().zip(&fdofs).map(|(&na, &gd)| na * t_full[gd]).sum();
                if !tq.is_finite() || tq < 0.0 {
                    out.fill(f64::NAN);
                    return;
                }
                let wq = wq_ref * ds * alpha * (tq.powi(4) - t_amb4);
                for (a, &gd) in fdofs.iter().enumerate() {
                    let ra = inv_free[gd];
                    if ra != usize::MAX {
                        out[ra] += wq * nv[a];
                    }
                }
            }
        }
    };

    let eval_mixed_stat_jac = |u_free: &[f64], values: &mut [f64]| {
        values.fill(0.0);
        copy_csr_into_template(&k_lin_ff, &indptr, &indices, values);
        let mut t_full = vec![0.0; n];
        for (&gd, &v) in fixed.iter().zip(vals) {
            t_full[gd] = v;
        }
        for (r, &gd) in free.iter().enumerate() {
            t_full[gd] = u_free[r];
        }
        if has_k_curve {
            let mut vals_s = Vec::new();
            let mut dxi = Vec::new();
            let mut sn = Vec::new();
            let mut sdn = Vec::new();
            for c in 0..mixed.n_cells() {
                let ct = mixed.cell_type(c);
                let iso = IsoparametricCell::from_mesh(mixed, c);
                let rule = quad_rule_mixed(ct, (2 * order).max(2));
                let cdofs = dof_map.cell_dofs(c);
                let nloc = cdofs.len();
                let dom = mixed.cell_labels[c];
                for (xi, &wq_ref) in rule.pts.iter().zip(&rule.w) {
                    let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                    lagrange_shape_into(ct, space, xi, &mut vals_s, &mut dxi);
                    let grads: Vec<[f64; 3]> = dxi.iter().map(|&g_ref| jac.map_grad(g_ref)).collect();
                    let mut tq = 0.0;
                    let mut gt = [0.0; 3];
                    for (a, &gd) in cdofs.iter().enumerate() {
                        let ta = t_full[gd];
                        tq += ta * vals_s[a];
                        for i in 0..dim {
                            gt[i] += ta * grads[a][i];
                        }
                    }
                    let (kc, dkc) = match spec.k_curves.get(&dom) {
                        Some((curve, euler)) => (
                            femcore::section::plane_tensor(&curve.eval_matrix_global(tq, *euler), dim),
                            femcore::section::plane_tensor(&curve.eval_matrix_deriv_global(tq, *euler), dim),
                        ),
                        None => (ktensor[c], [[0.0; 3]; 3]),
                    };
                    let mut dflux = [0.0; 3];
                    for i in 0..dim {
                        for j in 0..dim {
                            dflux[i] += dkc[i][j] * gt[j];
                        }
                    }
                    let wq = wq_ref * jac.det.abs();
                    for a in 0..nloc {
                        let ra = inv_free[cdofs[a]];
                        if ra == usize::MAX {
                            continue;
                        }
                        let mut ga_dot_df = 0.0;
                        let mut kg_a = [0.0; 3];
                        for i in 0..dim {
                            ga_dot_df += grads[a][i] * dflux[i];
                            for j in 0..dim {
                                kg_a[i] += kc[j][i] * grads[a][j];
                            }
                        }
                        for b in 0..nloc {
                            let rb = inv_free[cdofs[b]];
                            if rb != usize::MAX {
                                let mut diff_ab = 0.0;
                                for i in 0..dim {
                                    diff_ab += kg_a[i] * grads[b][i];
                                }
                                add_to_csr_template(
                                    &indptr,
                                    &indices,
                                    values,
                                    ra,
                                    rb,
                                    wq * (diff_ab + ga_dot_df * vals_s[b]),
                                );
                            }
                        }
                    }
                }
            }
        }
        let fdeg = 2 * order + 4;
        let mut fdofs = Vec::new();
        for &(fi, alpha, _, _) in &rad_facets {
            let ft = mixed.facet_type(fi);
            let frule = facet_rule_mixed(ft, fdeg);
            dof_map.facet_dofs(mixed, fi, &mut fdofs);
            for (xi, &wq_ref) in frule.pts.iter().zip(&frule.w) {
                let (_, ds, _) = eval_facet_jacobian(mixed, fi, xi);
                let (nv, _) = facet_lagrange_shape(ft, space, xi);
                let tq: f64 = nv.iter().zip(&fdofs).map(|(&na, &gd)| na * t_full[gd]).sum();
                let wq = wq_ref * ds * 4.0 * alpha * tq.max(0.0).powi(3);
                for (a, &ga) in fdofs.iter().enumerate() {
                    let ra = inv_free[ga];
                    if ra == usize::MAX {
                        continue;
                    }
                    for (b, &gb) in fdofs.iter().enumerate() {
                        let rb = inv_free[gb];
                        if rb != usize::MAX {
                            add_to_csr_template(&indptr, &indices, values, ra, rb, wq * nv[a] * nv[b]);
                        }
                    }
                }
            }
        }
    };

    let mut log = Vec::new();
    let mut diag = DfsDiag::new("snes_newton");
    let (times, temps_mixed) = match &spec.study {
        Study::HeatStationary => {
            let mut ref_sum = 0.0;
            let mut ref_cnt = 0usize;
            for &v in vals {
                if v > 0.0 {
                    ref_sum += v;
                    ref_cnt += 1;
                }
            }
            for &(_, _, _, ta) in &rad_facets {
                if ta > 0.0 {
                    ref_sum += ta;
                    ref_cnt += 1;
                }
            }
            let t_ref_init = if ref_cnt > 0 { (ref_sum / ref_cnt as f64).max(4.0) } else { 293.15 };
            let mut k_init = k_base.clone();
            let mut f_init = f.to_vec();
            for bc in &spec.bcs {
                if let BcKind::Radiation { emissivity, t_ambient, view_factor } = bc.kind {
                    let alpha = emissivity * view_factor * SIGMA_SB;
                    if alpha > 0.0 {
                        let h0 = 4.0 * alpha * t_ref_init.powi(3);
                        let q0 = alpha * (3.0 * t_ref_init.powi(4) + t_ambient.powi(4));
                        k_init = k_init.add_scaled(&assemble_boundary_mass_mixed(mixed, dof_map, &bc.labels, h0), 1.0);
                        let l0 = assemble_boundary_load_mixed(mixed, dof_map, &bc.labels, q0);
                        for (fi, &li) in f_init.iter_mut().zip(&l0) {
                            *fi += li;
                        }
                    }
                }
            }
            let k_init_ff = k_init.submatrix(&free, &free);
            let k_init_fx = k_init.submatrix(&free, fixed);
            let kx0 = if fixed.is_empty() { vec![0.0; n_free] } else { k_init_fx.matvec(vals) };
            let rhs0: Vec<f64> = free.iter().enumerate().map(|(r, &gd)| f_init[gd] - kx0[r]).collect();
            let mut u0 = DirectSolver::new(&k_init_ff, true)
                .map(|s| s.solve(&rhs0))
                .unwrap_or_else(|_| vec![t_ref_init; n_free]);
            if has_rad {
                for v in &mut u0 {
                    if !v.is_finite() || *v <= 0.0 {
                        *v = t_ref_init.max(1.0);
                    }
                }
            }
            let rhs_scale = rhs0.iter().map(|x| x * x).sum::<f64>().sqrt().max(1.0);
            let diag_ku_scale: f64 =
                k_init_ff.diagonal().iter().zip(&u0).map(|(d, u)| (d * u).powi(2)).sum::<f64>().sqrt();
            let res_scale = rhs_scale.max(diag_ku_scale);
            let user_tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
            let n_opts = NewtonOptions {
                rtol: user_tol.min(1e-10),
                atol: (user_tol * 1e-2).clamp(1e-13, 1e-10) * res_scale,
                stol: 1e-14,
                max_it: 50,
                line_search: LineSearch::Backtracking { c: 1e-4, alpha_min: 1e-8 },
                definiteness: if has_k_curve { Definiteness::General } else { Definiteness::PositiveDefinite },
                jacobian_lag: 1,
            };
            let mut f_fn = |u: &[f64], fx: &mut [f64]| eval_mixed_stat_res(u, fx);
            let mut jac_fn = |u: &[f64]| -> Csr {
                let mut values = vec![0.0; indptr[n_free]];
                eval_mixed_stat_jac(u, &mut values);
                Csr::try_new(n_free, n_free, indptr.clone(), indices.clone(), values)
                    .expect("fixed-pattern mixed thermal Jacobian CSR")
            };
            let nl = newton(&mut f_fn, &mut jac_fn, &u0, &n_opts)
                .map_err(|e| format!("nonlinear mixed thermal Newton failed: {e}"))?;
            let last_norm = nl.residual_norms.last().copied().unwrap_or(f64::NAN);
            if !nl.reason.is_converged() {
                return Err(format!(
                    "nonlinear mixed thermal Newton failed to converge ({:?} after {} iterations, residual {:.3e})",
                    nl.reason, nl.iterations, last_norm
                ));
            }
            let mut t_sol = vec![0.0; n];
            for (&gd, &v) in fixed.iter().zip(vals) {
                t_sol[gd] = v;
            }
            for (r, &gd) in free.iter().enumerate() {
                t_sol[gd] = nl.x[r];
            }
            diag.iterations = Some(nl.iterations);
            diag.converged = true;
            diag.max_algebraic_residual = last_norm / res_scale;
            log.push(format!(
                "nonlinear mixed heat Newton: {} iteration(s), residual {:.2e}",
                nl.iterations, diag.max_algebraic_residual
            ));
            (vec![], vec![t_sol])
        }
        Study::HeatTransient { t_end, n_out, t0: t0_expr } => {
            if *t_end <= 0.0 || *n_out < 1 {
                return Err("end time must be positive and at least one output time is needed".into());
            }
            let rho = cell_values_mixed(&spec.materials, mixed, "rho", Some(1.0))?;
            let cp_const = cell_values_mixed(&spec.materials, mixed, "Cp", Some(1.0))?;
            let has_cp_curve = !spec.cp_curves.is_empty();
            let rcp: Vec<f64> = rho.iter().zip(&cp_const).map(|(a, b)| a * b).collect();
            let m_const = assemble_mass_mixed(mixed, dof_map, &rcp);
            let mff = m_const.submatrix(&free, &free);

            let times: Vec<f64> = (0..=*n_out).map(|j| t_end * j as f64 / *n_out as f64).collect();
            let sub = 8usize;
            let dt = times[1] / sub as f64;
            let coords = dof_map.dof_coordinates(mixed);
            let mut y0 = eval_at(t0_expr, &coords, dim).map_err(|e| format!("initial temperature: {e}"))?;
            for (&gd, &v) in fixed.iter().zip(vals) {
                y0[gd] = v;
            }
            let mut out = vec![vec![0.0; n]; times.len()];
            out[0].copy_from_slice(&y0);
            let mut u_prev: Vec<f64> = free.iter().map(|&gd| y0[gd]).collect();
            let mut u_cur = u_prev.clone();
            let def = if has_k_curve || has_cp_curve { Definiteness::General } else { Definiteness::PositiveDefinite };
            let mut total_iters = 0usize;
            let mut max_res = 0.0_f64;
            for step in 1..=sub * n_out {
                femcore::progress::check()?;
                let alpha_0 = if step == 1 { 1.0 } else { 1.5 };
                let u_hist: Vec<f64> = if step == 1 {
                    u_cur.clone()
                } else {
                    (0..n_free).map(|r| 2.0 * u_cur[r] - 0.5 * u_prev[r]).collect()
                };
                let mut t_hist_full = vec![0.0; n];
                for (&gd, &v) in fixed.iter().zip(vals) {
                    t_hist_full[gd] = (alpha_0 - 1.0) * v + (2.0 - alpha_0) * v;
                }
                if step == 1 {
                    for (&gd, &v) in fixed.iter().zip(vals) {
                        t_hist_full[gd] = v;
                    }
                } else {
                    for (&gd, &v) in fixed.iter().zip(vals) {
                        t_hist_full[gd] = 1.5 * v;
                    }
                }
                for (r, &gd) in free.iter().enumerate() {
                    t_hist_full[gd] = u_hist[r];
                }
                let mut step_res = |u: &[f64], fx: &mut [f64]| {
                    eval_mixed_stat_res(u, fx);
                    if fx.first().is_some_and(|v| !v.is_finite()) {
                        return;
                    }
                    if !has_cp_curve {
                        let du: Vec<f64> = u.iter().zip(&u_hist).map(|(&ui, &uhi)| (alpha_0 * ui - uhi) / dt).collect();
                        let m_du = mff.matvec(&du);
                        for r in 0..n_free {
                            fx[r] += m_du[r];
                        }
                    } else {
                        let mut t_full = vec![0.0; n];
                        for (&gd, &v) in fixed.iter().zip(vals) {
                            t_full[gd] = v;
                        }
                        for (r, &gd) in free.iter().enumerate() {
                            t_full[gd] = u[r];
                        }
                        let mut vals_s = Vec::new();
                        let mut dxi = Vec::new();
                        let mut sn = Vec::new();
                        let mut sdn = Vec::new();
                        for c in 0..mixed.n_cells() {
                            let ct = mixed.cell_type(c);
                            let iso = IsoparametricCell::from_mesh(mixed, c);
                            let rule = quad_rule_mixed(ct, (2 * order).max(2));
                            let cdofs = dof_map.cell_dofs(c);
                            let dom = mixed.cell_labels[c];
                            let rho_c = rho[c];
                            for (xi, &wq_ref) in rule.pts.iter().zip(&rule.w) {
                                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                                lagrange_shape_into(ct, space, xi, &mut vals_s, &mut dxi);
                                let mut tq = 0.0;
                                let mut thq = 0.0;
                                for (a, &gd) in cdofs.iter().enumerate() {
                                    tq += t_full[gd] * vals_s[a];
                                    thq += t_hist_full[gd] * vals_s[a];
                                }
                                let cp_q = spec.cp_curves.get(&dom).map_or(cp_const[c], |cv| cv.eval(tq));
                                let wq = wq_ref * jac.det.abs() * rho_c * cp_q * (alpha_0 * tq - thq) / dt;
                                for (a, &gd) in cdofs.iter().enumerate() {
                                    let ra = inv_free[gd];
                                    if ra != usize::MAX {
                                        fx[ra] += wq * vals_s[a];
                                    }
                                }
                            }
                        }
                    }
                };
                let mut step_jac = |u: &[f64]| -> Csr {
                    let mut values = vec![0.0; indptr[n_free]];
                    eval_mixed_stat_jac(u, &mut values);
                    if !has_cp_curve {
                        let sm = alpha_0 / dt;
                        for r in 0..mff.nrows() {
                            let (cols, vs) = mff.row(r);
                            let range = indptr[r]..indptr[r + 1];
                            for (&col, &v) in cols.iter().zip(vs) {
                                if v != 0.0 {
                                    if let Ok(pos) = indices[range.clone()].binary_search(&col) {
                                        values[range.start + pos] += sm * v;
                                    }
                                }
                            }
                        }
                    } else {
                        let mut t_full = vec![0.0; n];
                        for (&gd, &v) in fixed.iter().zip(vals) {
                            t_full[gd] = v;
                        }
                        for (r, &gd) in free.iter().enumerate() {
                            t_full[gd] = u[r];
                        }
                        let mut vals_s = Vec::new();
                        let mut dxi = Vec::new();
                        let mut sn = Vec::new();
                        let mut sdn = Vec::new();
                        for c in 0..mixed.n_cells() {
                            let ct = mixed.cell_type(c);
                            let iso = IsoparametricCell::from_mesh(mixed, c);
                            let rule = quad_rule_mixed(ct, (2 * order).max(2));
                            let cdofs = dof_map.cell_dofs(c);
                            let dom = mixed.cell_labels[c];
                            let rho_c = rho[c];
                            for (xi, &wq_ref) in rule.pts.iter().zip(&rule.w) {
                                let jac = iso.eval_into(xi, &mut sn, &mut sdn);
                                lagrange_shape_into(ct, space, xi, &mut vals_s, &mut dxi);
                                let mut tq = 0.0;
                                let mut thq = 0.0;
                                for (a, &gd) in cdofs.iter().enumerate() {
                                    tq += t_full[gd] * vals_s[a];
                                    thq += t_hist_full[gd] * vals_s[a];
                                }
                                let (cp_q, dcp_q) = spec
                                    .cp_curves
                                    .get(&dom)
                                    .map_or((cp_const[c], 0.0), |cv| (cv.eval(tq), cv.eval_deriv(tq)));
                                let wq = wq_ref
                                    * jac.det.abs()
                                    * (rho_c / dt)
                                    * (alpha_0 * cp_q + dcp_q * (alpha_0 * tq - thq));
                                for (a, &ga) in cdofs.iter().enumerate() {
                                    let ra = inv_free[ga];
                                    if ra == usize::MAX {
                                        continue;
                                    }
                                    for (b, &gb) in cdofs.iter().enumerate() {
                                        let rb = inv_free[gb];
                                        if rb != usize::MAX {
                                            add_to_csr_template(
                                                &indptr,
                                                &indices,
                                                &mut values,
                                                ra,
                                                rb,
                                                wq * vals_s[a] * vals_s[b],
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    Csr::try_new(n_free, n_free, indptr.clone(), indices.clone(), values)
                        .expect("fixed-pattern mixed transient Jacobian CSR")
                };
                let ref_norm = u_cur.iter().map(|x| x * x).sum::<f64>().sqrt().max(1.0);
                let diag_ku_scale: f64 =
                    k_lin_ff.diagonal().iter().zip(&u_cur).map(|(d, u)| (d * u).powi(2)).sum::<f64>().sqrt();
                let res_scale = (ref_norm / dt.max(1e-12)).max(diag_ku_scale);
                let n_opts = NewtonOptions {
                    rtol: 1e-11,
                    atol: 1e-12 * res_scale,
                    stol: 1e-14,
                    max_it: 30,
                    line_search: LineSearch::Backtracking { c: 1e-4, alpha_min: 1e-8 },
                    definiteness: def,
                    jacobian_lag: 1,
                };
                let nl = newton(&mut step_res, &mut step_jac, &u_cur, &n_opts)
                    .map_err(|e| format!("nonlinear mixed transient Newton step {step} failed: {e}"))?;
                let last_norm = nl.residual_norms.last().copied().unwrap_or(f64::NAN);
                if !nl.reason.is_converged() {
                    return Err(format!(
                        "nonlinear mixed transient Newton step {step} failed to converge ({:?})",
                        nl.reason
                    ));
                }
                total_iters += nl.iterations;
                max_res = max_res.max(last_norm / nl.residual_norms.first().copied().unwrap_or(1.0).max(res_scale));
                u_prev = std::mem::replace(&mut u_cur, nl.x);
                if step % sub == 0 {
                    let idx = step / sub;
                    for (&gd, &v) in fixed.iter().zip(vals) {
                        out[idx][gd] = v;
                    }
                    for (r, &gd) in free.iter().enumerate() {
                        out[idx][gd] = u_cur[r];
                    }
                }
            }
            diag.solver_kind = "snes_newton_bdf2".into();
            diag.iterations = Some(total_iters);
            diag.converged = true;
            diag.max_algebraic_residual = max_res;
            log.push(format!("nonlinear mixed BDF2, {} steps ({} total Newton iterations)", sub * n_out, total_iters));
            (times, out)
        }
        Study::Electrostatic | Study::Conduction | Study::Magnetostatic { .. } => {
            return Err("heat transfer supports the stationary and transient studies".into());
        }
    };
    let (diagnostics, error_budget) = diag.finish();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut solutions = Vec::with_capacity(temps_mixed.len());
    for (j, row) in temps_mixed.iter().enumerate() {
        let (fields, cell_fields, t_simplex) = heat_fields_mixed(mixed, dof_map, simplex, proj, row, ktensor, order);
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
        dim: dim as u8,
        solutions,
        scalars: vec![("T min (K)".into(), Scalar::Num(lo)), ("T max (K)".into(), Scalar::Num(hi))],
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}
