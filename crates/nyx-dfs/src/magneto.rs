// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Gauge-stabilized London magnetostatics and multi-terminal inductance matrix extraction (`study = "magnetostatic"`).
//!
//! Three formulations are provided:
//! * `cross_section` (2D): out-of-plane vector potential `A_z in H^1(Omega_2D)` with `(mu_0 mu_r)^(-1) grad A_z . grad v`
//!   and London sheet kinetic admittance `int_{Gamma_SC} L_k_sq^(-1) A_z v ds`, yielding per-unit-length inductance `L'_{ij}` (H/m).
//! * `planar_stream` (2D): Brandt-Clem stream function `g in H_0^1(Gamma_SC)` (`K_s = curl(g z_hat)`) combining kinetic sheet
//!   inductance `int L_k_sq grad g . grad psi dA` and the symmetric 2D Biot-Savart geometric kernel with exact analytical polar
//!   self-term integration on coincident triangles, yielding `L_{ij}` (H).
//! * `mixed_3d` (3D): coupled lowest-order Nedelec `A in H(curl, Omega)` and 2D surface scalar potential `phi_s in H^1(Gamma_SC)`
//!   on single-facet deduplicated `london_sheet` triangles with gauge regularization `eps_reg = 1e-6 / L_domain^2`, pinning one
//!   reference `phi_s` vertex per connected sheet component and supporting `current_terminal` and `flux_loop` cohomology excitations,
//!   yielding `L_{ij}` (H).

use crate::common::{
    assemble_boundary_load_mixed, assemble_boundary_mass_mixed, assemble_farfield_robin_mixed,
    assemble_ie_outer_robin_mixed, assemble_laplace_ie_mixed, boundary_measure_mixed, build_facet_to_cell_mixed,
    cell_gradient, cell_values_mixed, check_boundary_labels_mixed, classify_ie_outer_facets_mixed, components,
    farfield_centroid_mixed, mesh_has_curved_edges, minmax, norm_of, quad_seg_high, resolve_ie_domains_mixed,
    sample_mixed_gradients_on_simplex, ResolvedIeDomain,
};
use crate::problem::{BcKind, Problem, Study};
use femcore::edge_mixed::{eval_cell_basis_3d, NedEvalScratch, NedMixed};
use femcore::fem::*;
use femcore::fem_mixed::{
    cell_corner_ref_coords, eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, quad_rule_mixed,
    IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{Coo, RCsr as Csr};
use femcore::linalg::driver::{dense_eigh, DirectSolver};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::props::LinearSolver;
use femcore::vvuq::{SourceField, UniformSource};
use nyx_physics::{
    cell_fields_of, discretise, discretise_mixed, Discretisation, Field, InductanceMatrixData, Mesh,
    MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::time::Instant;

/// Vacuum magnetic permeability `mu_0` in H/m.
pub const MU0: f64 = 4.0 * PI * 1e-7;

const TET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

#[derive(Clone, Debug)]
struct Excitation {
    name: String,
    labels: Vec<i32>,
    amplitude: f64,
    direction: [f64; 3],
    is_flux_loop: bool,
}

fn collect_excitations(spec: &Problem) -> Result<Vec<Excitation>, String> {
    let mut out: Vec<Excitation> = Vec::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::CurrentTerminal { name, current, direction } => {
                let name = name.trim();
                if name.is_empty() {
                    return Err("current_terminal condition requires a non-empty 'name'".into());
                }
                let amp = if current.abs() > 1e-30 { *current } else { 1.0 };
                if let Some(ex) = out.iter_mut().find(|e| e.name == name) {
                    for &l in &bc.labels {
                        if !ex.labels.contains(&l) {
                            ex.labels.push(l);
                        }
                    }
                } else {
                    out.push(Excitation {
                        name: name.to_string(),
                        labels: bc.labels.clone(),
                        amplitude: amp,
                        direction: *direction,
                        is_flux_loop: false,
                    });
                }
            }
            BcKind::FluxLoop { name, flux_wb } => {
                let name = name.trim();
                if name.is_empty() {
                    return Err("flux_loop condition requires a non-empty 'name'".into());
                }
                let amp = if flux_wb.abs() > 1e-30 { *flux_wb } else { 1.0 };
                if let Some(ex) = out.iter_mut().find(|e| e.name == name) {
                    for &l in &bc.labels {
                        if !ex.labels.contains(&l) {
                            ex.labels.push(l);
                        }
                    }
                } else {
                    out.push(Excitation {
                        name: name.to_string(),
                        labels: bc.labels.clone(),
                        amplitude: amp,
                        direction: [0.0, 0.0, 1.0],
                        is_flux_loop: true,
                    });
                }
            }
            BcKind::FarField { order: ff_order, decay_power, center } => {
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
            BcKind::LondonSheet { .. } | BcKind::Ground | BcKind::InfiniteElement(_) => {}
            other => return Err(format!("'{}' is not a magnetostatics boundary condition", other.name())),
        }
    }
    if out.is_empty() {
        return Err("magnetostatics needs at least one current_terminal or flux_loop condition".into());
    }
    Ok(out)
}

fn push_inductance_scalars(
    scalars: &mut Vec<(String, Scalar)>,
    excitations: &[Excitation],
    l_mat: &[Vec<f64>],
    formulation: &str,
) {
    let n = excitations.len();
    for i in 0..n {
        for j in 0..n {
            let name_key = format!("L_{}_{} (H)", excitations[i].name, excitations[j].name);
            let idx_key = format!("L_{i}_{j} (H)");
            scalars.push((name_key.clone(), Scalar::Num(l_mat[i][j])));
            if idx_key != name_key {
                scalars.push((idx_key, Scalar::Num(l_mat[i][j])));
            }
        }
    }
    let data = InductanceMatrixData {
        terminals: excitations.iter().map(|e| e.name.clone()).collect(),
        inductance_matrix_h: l_mat.to_vec(),
        formulation: formulation.to_string(),
    };
    scalars.push(data.scalar());
}

// ------------------------------------------------------------------ 1) 2D cross_section

fn cross_section_fields(mesh: &Discretisation, order: usize, az: &[f64]) -> (Vec<Field>, CellFields) {
    let ones = vec![1.0; mesh.n_cells()];
    let grad_az = cell_gradient(mesh, az, order, &ones);
    // In 2D (x, z), A = A_z y_hat => B = curl(A_z y_hat) = (-dA_z/dz, dA_z/dx)
    let mut bcell = vec![0.0; grad_az.len()];
    for p in 0..(grad_az.len() / 2) {
        let dax = grad_az[2 * p];
        let daz = grad_az[2 * p + 1];
        bcell[2 * p] = -daz;
        bcell[2 * p + 1] = dax;
    }
    let b_avg = node_average(mesh, &bcell, 2, order);
    let bc = components(&b_avg, 2);
    let names = femcore::section::component_names(2);
    let fields = vec![Field::scalar("A", "T m", az.to_vec()), Field::vector("B", "T", names, &bc)];
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("A", "T m", cell_values_at_nodes(mesh, az, order));
    cf.push_components("B", names, "T", &bcell);
    (fields, cf)
}

fn reject_curved_simplex(mesh: &Discretisation) -> Result<(), String> {
    if mesh.has_curved() {
        return Err("unsupported: magnetostatics does not support curved edges (R32)".into());
    }
    Ok(())
}

fn solve_cross_section(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    reject_curved_simplex(mesh)?;
    if spec.has_infinite_element() {
        let cm = Mesh::from_solver_mesh(mesh).map_err(|e| e.to_string())?;
        return solve_cross_section_mixed(&cm, spec);
    }
    let t0 = Instant::now();
    if mesh.dim != 2 {
        return Err("the 'cross_section' magnetostatic formulation requires a 2D mesh".into());
    }
    let order = spec.order;
    check_boundary_labels(
        mesh,
        spec.bcs
            .iter()
            .filter(|b| !matches!(b.kind, BcKind::InfiniteElement(_)))
            .flat_map(|b| b.labels.iter().copied()),
    )?;
    let excitations = collect_excitations(spec)?;
    let mu_r = spec.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    if mu_r.iter().any(|&m| m <= 0.0) {
        return Err("material property 'mu_r' must be positive".into());
    }
    let inv_mur: Vec<f64> = mu_r.iter().map(|&m| 1.0 / m).collect();
    let coeff: Vec<f64> = inv_mur.iter().map(|&im| im / MU0).collect();
    let mut k = assemble_laplace(mesh, order, &coeff);

    let mut sheet_lk: BTreeMap<i32, f64> = BTreeMap::new();
    let mut ground_labels: Vec<i32> = Vec::new();
    let mut ff_label_set: BTreeSet<i32> = BTreeSet::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::LondonSheet { l_k_sq } => {
                for &lab in &bc.labels {
                    sheet_lk.insert(lab, *l_k_sq);
                }
            }
            BcKind::Ground => {
                for &lab in &bc.labels {
                    if !ground_labels.contains(&lab) {
                        ground_labels.push(lab);
                    }
                }
            }
            BcKind::FarField { .. } => {
                for &lab in &bc.labels {
                    ff_label_set.insert(lab);
                }
            }
            _ => {}
        }
    }

    let term_label_set: BTreeSet<i32> = excitations.iter().flat_map(|e| e.labels.iter().copied()).collect();
    let ref_labels: Vec<i32> = term_label_set.iter().copied().collect();
    let mut has_far_field = false;
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !ground_labels.contains(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid(mesh, &ref_labels));
                let k_far = assemble_farfield_robin(
                    mesh,
                    &ff_labs,
                    order,
                    *ff_order,
                    *decay_power,
                    c,
                    1.0 / MU0,
                    FarFieldPermittivity::Scalar(&inv_mur),
                );
                k = k.add_scaled(&k_far, 1.0);
                has_far_field = true;
            }
        }
    }

    // Passive London sheets (not driven as current terminals): add int_{Gamma_SC} L_k_sq^(-1) A_z v ds
    for (&lab, &l_k_sq) in &sheet_lk {
        if !term_label_set.contains(&lab) {
            if l_k_sq > 0.0 {
                k = k.add_scaled(&assemble_boundary_mass(mesh, &[lab], order, 1.0 / l_k_sq), 1.0);
            } else if !ground_labels.contains(&lab) {
                ground_labels.push(lab);
            }
        }
    }

    if ground_labels.is_empty() && !has_far_field {
        let exterior = mesh.exterior_facets();
        for f in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f];
            if exterior[f]
                && !term_label_set.contains(&lab)
                && !sheet_lk.contains_key(&lab)
                && !ff_label_set.contains(&lab)
                && !ground_labels.contains(&lab)
            {
                ground_labels.push(lab);
            }
        }
    }
    let fixed = boundary_dofs(mesh, &ground_labels, order);
    if fixed.is_empty() && !has_far_field {
        return Err(
            "2D cross_section magnetostatics needs an exterior, ground, or far_field reference boundary (A_z = 0)"
                .into(),
        );
    }

    let n = k.nrows();
    let free = free_dofs(n, &fixed);
    let n_free = free.len();
    let n_exc = excitations.len();

    let mut rhs_full_list = Vec::with_capacity(n_exc);
    let mut term_kin_self = vec![0.0; n_exc];
    for (p, ex) in excitations.iter().enumerate() {
        let width: f64 = facets_with(mesh, &ex.labels).map(|f| facet_measure(mesh, f)).sum();
        if width <= 0.0 {
            return Err(format!("current_terminal '{}' has zero boundary measure", ex.name));
        }
        let j_surf = ex.amplitude / width;
        rhs_full_list.push(assemble_boundary_load(mesh, &ex.labels, j_surf, order));

        // When a driven current terminal also carries LondonSheet kinetic inductance L_k,sq,
        // the longitudinal kinetic energy 1/2 int_{Gamma_p} L_k,sq (I_p / w_p)^2 ds adds series L'_k = L_k,sq / w_p.
        let mut lk_int = 0.0;
        for f in facets_with(mesh, &ex.labels) {
            if let Some(&lk) = sheet_lk.get(&mesh.facet_labels[f]) {
                if lk > 0.0 {
                    lk_int += lk * facet_measure(mesh, f);
                }
            }
        }
        term_kin_self[p] = lk_int / (width * width);
    }

    let kff = k.submatrix(&free, &free);
    let rhs_free_list: Vec<Vec<f64>> = rhs_full_list.iter().map(|r| free.iter().map(|&i| r[i]).collect()).collect();
    let solver = DirectSolver::new(&kff, true)?;
    let sols_free = solver.solve_many(&rhs_free_list);
    let log = vec![format!(
        "2D cross_section magnetostatics: 1 factorization ({n_free} free dofs) served {n_exc} excitation(s)"
    )];

    let mut az_all = vec![vec![0.0; n]; n_exc];
    for p in 0..n_exc {
        for (a, &i) in free.iter().enumerate() {
            az_all[p][i] = sols_free[p][a];
        }
    }

    let kaz_all: Vec<Vec<f64>> = az_all.iter().map(|a| k.matvec(a)).collect();
    let mut l_mat = vec![vec![0.0; n_exc]; n_exc];
    for i in 0..n_exc {
        for j in i..n_exc {
            let e_ij = az_all[i].iter().zip(&kaz_all[j]).map(|(a, b)| a * b).sum::<f64>();
            let e_ji = if i == j { e_ij } else { az_all[j].iter().zip(&kaz_all[i]).map(|(a, b)| a * b).sum::<f64>() };
            let mut lij = 0.5 * (e_ij + e_ji) / (excitations[i].amplitude * excitations[j].amplitude);
            if i == j {
                lij += term_kin_self[i];
            }
            l_mat[i][j] = lij;
            l_mat[j][i] = lij;
        }
    }

    let mut solutions = Vec::with_capacity(n_exc);
    let mut scalars = Vec::new();
    for p in 0..n_exc {
        let (fields, cell_fields) = cross_section_fields(mesh, order, &az_all[p]);
        if p == 0 {
            let bn = norm_of(&components(&fields[1].values, 2));
            let w_m = 0.5 * l_mat[0][0] * excitations[0].amplitude.powi(2);
            scalars.push(("stored magnetic energy W_m (J/m)".into(), Scalar::Num(w_m)));
            scalars.push(("max |B| (T)".into(), Scalar::Num(minmax(&bn).1)));
        }
        solutions.push(Solution {
            label: format!("terminal {}", excitations[p].name),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(l_mat[p][p]),
            ..Default::default()
        });
    }

    push_inductance_scalars(&mut scalars, &excitations, &l_mat, "cross_section");

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
    })
}

// ------------------------------------------------------------------ 2) 2D planar_stream

/// Exact analytical integral `int_T dA' / |p - r'|` for an interior point `p` inside triangle `tri = [v0, v1, v2]`
/// in 2D polar coordinates `(r, theta)` around `p`, where `int_0^{R(theta)} (1/r) r dr = R(theta)` cancels the `1/r`
/// singularity identically.
fn triangle_potential_at(p: [f64; 2], tri: [[f64; 2]; 3]) -> f64 {
    let mut acc = 0.0;
    for i in 0..3 {
        let a = [tri[i][0] - p[0], tri[i][1] - p[1]];
        let b = [tri[(i + 1) % 3][0] - p[0], tri[(i + 1) % 3][1] - p[1]];
        let e = [b[0] - a[0], b[1] - a[1]];
        let len = (e[0] * e[0] + e[1] * e[1]).sqrt();
        if len <= 1e-30 {
            continue;
        }
        let t = [e[0] / len, e[1] / len];
        let h_perp = (a[0] * e[1] - a[1] * e[0]).abs() / len;
        if h_perp <= 1e-30 {
            continue;
        }
        let s0 = a[0] * t[0] + a[1] * t[1];
        let s1 = b[0] * t[0] + b[1] * t[1];
        acc += h_perp * ((s1 / h_perp).asinh() - (s0 / h_perp).asinh()).abs();
    }
    acc
}

fn solve_planar_stream(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    reject_curved_simplex(mesh)?;
    let t0 = Instant::now();
    if mesh.dim != 2 {
        return Err("the 'planar_stream' magnetostatic formulation requires a 2D mesh".into());
    }
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
    let excitations = collect_excitations(spec)?;

    let mut l_k_sq = 0.0_f64;
    let mut ground_labels: Vec<i32> = Vec::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::LondonSheet { l_k_sq: lk } => {
                if *lk > l_k_sq {
                    l_k_sq = *lk;
                }
            }
            BcKind::Ground => {
                for &lab in &bc.labels {
                    if !ground_labels.contains(&lab) {
                        ground_labels.push(lab);
                    }
                }
            }
            _ => {}
        }
    }

    let term_label_set: BTreeSet<i32> = excitations.iter().flat_map(|e| e.labels.iter().copied()).collect();
    if ground_labels.is_empty() {
        // Choose the longest non-excited boundary label as reference g = 0
        let mut measure_by_lab: BTreeMap<i32, f64> = BTreeMap::new();
        for f in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f];
            if !term_label_set.contains(&lab) {
                *measure_by_lab.entry(lab).or_insert(0.0) += facet_measure(mesh, f);
            }
        }
        if let Some((&best_lab, _)) =
            measure_by_lab.iter().max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        {
            ground_labels.push(best_lab);
        }
    }
    if ground_labels.is_empty() && excitations.len() == 1 {
        return Err("2D planar_stream magnetostatics needs a ground or outer reference boundary (g = 0)".into());
    }

    let n = mesh.n_dofs(order);
    let nc = mesh.n_cells();
    let re = RefElement::new(2, order, if order > 1 { 4 } else { 2 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nq = re.w.len();

    // Tabulate physical quadrature points, weights, and shape gradients per cell
    let mut q_pts = vec![[0.0_f64; 2]; nc * nq];
    let mut q_w = vec![0.0_f64; nc * nq];
    let mut q_grad = vec![[0.0_f64; 2]; nc * nq * nloc];
    let mut cell_tris = Vec::with_capacity(nc);
    let mut all_dofs = vec![0usize; nc * nloc];
    let mut tmp_dofs = Vec::new();
    let mut grad_buf = vec![0.0; nloc * 2];

    for c in 0..nc {
        let v = mesh.cell(c);
        let p0 = [mesh.point(v[0] as usize)[0], mesh.point(v[0] as usize)[1]];
        let p1 = [mesh.point(v[1] as usize)[0], mesh.point(v[1] as usize)[1]];
        let p2 = [mesh.point(v[2] as usize)[0], mesh.point(v[2] as usize)[1]];
        cell_tris.push([p0, p1, p2]);
        mesh.cell_dofs(c, order, &mut tmp_dofs);
        all_dofs[c * nloc..(c + 1) * nloc].copy_from_slice(&tmp_dofs);

        for q in 0..nq {
            let lam = &re.n[0]; // Note: for order=1, re.n[q] is barycentric; get barycentric from quad_rule
            let _ = lam;
            shape_gradients(&re, &geo, c, q, &mut grad_buf);
            q_w[c * nq + q] = re.w[q] * geo.vol[c];
            for a in 0..nloc {
                q_grad[(c * nq + q) * nloc + a] = [grad_buf[2 * a], grad_buf[2 * a + 1]];
            }
        }
    }
    let qrule = quad_rule(2, if order > 1 { 4 } else { 2 });
    for c in 0..nc {
        let [p0, p1, p2] = cell_tris[c];
        for q in 0..nq {
            let l = &qrule.lam[q];
            q_pts[c * nq + q] =
                [l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0], l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1]];
        }
    }

    // Start with kinetic inductance matrix K_kin = int_{Gamma_SC} L_k,sq grad g . grad psi dA
    let mut k_dense = vec![0.0_f64; n * n];
    if l_k_sq > 0.0 {
        let k_kin = assemble_laplace(mesh, order, &vec![l_k_sq; nc]);
        k_dense = k_kin.to_dense();
    }

    // Add symmetric 2D Biot-Savart geometric matrix:
    // K_BS(i, j) = (mu_0 / 4 pi) int int (curl(phi_i z) . curl(phi_j z)) / |r - r'| dA dA'
    //            = (mu_0 / 4 pi) int int (grad phi_i(r) . grad phi_j(r')) / |r - r'| dA dA'
    let bs_pref = MU0 / (4.0 * PI);
    let h_soft2 = {
        let min_vol = geo.vol.iter().copied().fold(f64::INFINITY, f64::min);
        1e-6 * min_vol.max(1e-30)
    };

    for c in 0..nc {
        let dofs_c = &all_dofs[c * nloc..(c + 1) * nloc];
        // Self-cell term (c == c') using exact analytical polar inner integral over cell_tris[c]
        let mut pot_self = 0.0;
        for q in 0..nq {
            pot_self += q_w[c * nq + q] * triangle_potential_at(q_pts[c * nq + q], cell_tris[c]);
        }
        let area_c = geo.vol[c].max(1e-30);
        let kernel_avg = bs_pref * (pot_self / (area_c * area_c));
        for q in 0..nq {
            let wq = q_w[c * nq + q];
            let gc = &q_grad[(c * nq + q) * nloc..(c * nq + q + 1) * nloc];
            for qp in 0..nq {
                let wqp = q_w[c * nq + qp];
                let gcp = &q_grad[(c * nq + qp) * nloc..(c * nq + qp + 1) * nloc];
                let w_pair = wq * wqp * kernel_avg;
                for a in 0..nloc {
                    let da = dofs_c[a];
                    for b in 0..nloc {
                        let db = dofs_c[b];
                        let dot = gc[a][0] * gcp[b][0] + gc[a][1] * gcp[b][1];
                        k_dense[da * n + db] += w_pair * dot;
                    }
                }
            }
        }

        // Distinct cell pairs (cp > c)
        for cp in (c + 1)..nc {
            let dofs_cp = &all_dofs[cp * nloc..(cp + 1) * nloc];
            for q in 0..nq {
                let rq = q_pts[c * nq + q];
                let wq = q_w[c * nq + q];
                let gc = &q_grad[(c * nq + q) * nloc..(c * nq + q + 1) * nloc];
                for qp in 0..nq {
                    let rqp = q_pts[cp * nq + qp];
                    let wqp = q_w[cp * nq + qp];
                    let dx = rq[0] - rqp[0];
                    let dy = rq[1] - rqp[1];
                    let dist = (dx * dx + dy * dy + h_soft2).sqrt();
                    let w_pair = bs_pref * wq * wqp / dist;
                    let gcp = &q_grad[(cp * nq + qp) * nloc..(cp * nq + qp + 1) * nloc];
                    for a in 0..nloc {
                        let da = dofs_c[a];
                        for b in 0..nloc {
                            let db = dofs_cp[b];
                            let val = w_pair * (gc[a][0] * gcp[b][0] + gc[a][1] * gcp[b][1]);
                            k_dense[da * n + db] += val;
                            k_dense[db * n + da] += val;
                        }
                    }
                }
            }
        }
    }

    let n_exc = excitations.len();
    let bg_dofs: BTreeSet<usize> = boundary_dofs(mesh, &ground_labels, order).into_iter().collect();
    let term_dofs: Vec<BTreeSet<usize>> =
        excitations.iter().map(|ex| boundary_dofs(mesh, &ex.labels, order).into_iter().collect()).collect();

    let mut fixed_set = bg_dofs.clone();
    for td in &term_dofs {
        fixed_set.extend(td.iter().copied());
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let free = free_dofs(n, &fixed);
    let n_free = free.len();

    let mut g_all = vec![vec![0.0; n]; n_exc];
    for p in 0..n_exc {
        let amp = excitations[p].amplitude;
        for &d in &fixed {
            if term_dofs[p].contains(&d) && !bg_dofs.contains(&d) {
                g_all[p][d] = amp;
            }
        }
    }

    if n_free > 0 {
        let mut kff = vec![0.0; n_free * n_free];
        for (a, &i) in free.iter().enumerate() {
            for (b, &j) in free.iter().enumerate() {
                kff[a * n_free + b] = 0.5 * (k_dense[i * n + j] + k_dense[j * n + i]);
            }
        }
        let (eigs, vecs) = dense_eigh(&kff, n_free)?;
        let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max).max(1e-30);
        for p in 0..n_exc {
            let mut rhs = vec![0.0; n_free];
            for (a, &i) in free.iter().enumerate() {
                let mut s = 0.0;
                for &d in &fixed {
                    s += k_dense[i * n + d] * g_all[p][d];
                }
                rhs[a] = -s;
            }
            let mut sol = vec![0.0; n_free];
            for k in 0..n_free {
                if eigs[k] > 1e-13 * lam_max {
                    let dot: f64 = vecs[k].iter().zip(&rhs).map(|(u, r)| u * r).sum();
                    let c = dot / eigs[k];
                    for a in 0..n_free {
                        sol[a] += c * vecs[k][a];
                    }
                }
            }
            for (a, &i) in free.iter().enumerate() {
                g_all[p][i] = sol[a];
            }
        }
    }

    let mut kg_all = vec![vec![0.0; n]; n_exc];
    for p in 0..n_exc {
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..n {
                s += k_dense[i * n + j] * g_all[p][j];
            }
            kg_all[p][i] = s;
        }
    }

    let mut l_mat = vec![vec![0.0; n_exc]; n_exc];
    for i in 0..n_exc {
        for j in i..n_exc {
            let e_ij: f64 = g_all[i].iter().zip(&kg_all[j]).map(|(a, b)| a * b).sum();
            let e_ji: f64 = if i == j { e_ij } else { g_all[j].iter().zip(&kg_all[i]).map(|(a, b)| a * b).sum() };
            let lij = 0.5 * (e_ij + e_ji) / (excitations[i].amplitude * excitations[j].amplitude);
            l_mat[i][j] = lij;
            l_mat[j][i] = lij;
        }
    }

    let mut solutions = Vec::with_capacity(n_exc);
    let mut scalars = Vec::new();
    for p in 0..n_exc {
        // If driven by a FluxLoop with prescribed flux_wb != 0, scale g so I_loop = flux_wb / L_pp
        let scale_g = if excitations[p].is_flux_loop && l_mat[p][p] > 0.0 { MU0 / l_mat[p][p] } else { MU0 };
        let a_scaled: Vec<f64> = g_all[p].iter().map(|&x| scale_g * x).collect();
        let (fields, cell_fields) = cross_section_fields(mesh, order, &a_scaled);
        if p == 0 {
            let w_m = 0.5 * l_mat[0][0] * excitations[0].amplitude.powi(2);
            scalars.push(("stored magnetic energy W_m (J)".into(), Scalar::Num(w_m)));
        }
        solutions.push(Solution {
            label: format!("terminal {}", excitations[p].name),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(l_mat[p][p]),
            ..Default::default()
        });
    }

    push_inductance_scalars(&mut scalars, &excitations, &l_mat, "planar_stream");
    let log = vec![format!(
        "2D planar_stream Brandt-Clem magnetostatics: 1 eigendecomposition ({n_free} free dofs, {nc} sheet triangles) served {n_exc} excitation(s)"
    )];

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
    })
}

// ------------------------------------------------------------------ 3) 3D mixed_3d

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn scale3(s: f64, a: [f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

/// Rigid-Euclidean and scaling-invariant 3D domain diameter (`2 * max_v |x_v - centroid|`).
fn domain_diameter_3d<F: Fn(usize) -> [f64; 3]>(n_pts: usize, pt: F) -> f64 {
    if n_pts == 0 {
        return 1.0;
    }
    let mut c = [0.0; 3];
    for i in 0..n_pts {
        let p = pt(i);
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    let inv = 1.0 / (n_pts as f64);
    c[0] *= inv;
    c[1] *= inv;
    c[2] *= inv;
    let mut r2_max = 0.0_f64;
    for i in 0..n_pts {
        let d = sub3(pt(i), c);
        r2_max = r2_max.max(dot3(d, d));
    }
    (2.0 * r2_max.sqrt()).max(1e-12)
}

/// Evaluates the 6 lowest-order Nedelec edge basis functions `w[0..6]` and their curls `curl[0..6]`
/// at barycentric point `lam` of tetrahedron `v`, oriented by global vertex order `v[i] < v[j]`.
fn whitney_tet(lam: &[f64; 4], g: &[[f64; 3]; 4], v: &[u32]) -> ([[f64; 3]; 6], [[f64; 3]; 6]) {
    let mut w = [[0.0; 3]; 6];
    let mut curl = [[0.0; 3]; 6];
    for (k, &(i, j)) in TET_EDGES.iter().enumerate() {
        let (a, b) = if v[i] < v[j] { (i, j) } else { (j, i) };
        w[k] = sub3(scale3(lam[a], g[b]), scale3(lam[b], g[a]));
        curl[k] = scale3(2.0, cross3(g[a], g[b]));
    }
    (w, curl)
}

/// Computes the area, unit normal, and tangential barycentric gradients `[g_0, g_1, g_2]` of a 3D triangle `[p0, p1, p2]`.
fn triangle_surface_gradients(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let e1 = sub3(p1, p0);
    let e2 = sub3(p2, p0);
    let n_raw = cross3(e1, e2);
    let n_norm2 = dot3(n_raw, n_raw);
    let n_norm = n_norm2.sqrt();
    if n_norm2 <= 1e-30 {
        return (0.0, [0.0, 0.0, 1.0], [[0.0; 3]; 3]);
    }
    let area = 0.5 * n_norm;
    let normal = scale3(1.0 / n_norm, n_raw);
    let g1 = scale3(1.0 / n_norm2, cross3(n_raw, sub3(p0, p2)));
    let g2 = scale3(1.0 / n_norm2, cross3(n_raw, sub3(p1, p0)));
    let g0 = [-g1[0] - g2[0], -g1[1] - g2[1], -g1[2] - g2[2]];
    (area, normal, [g0, g1, g2])
}

struct ExcitationFluxoidMeta {
    tris: BTreeSet<[u32; 3]>,
    is_flux_loop: bool,
    cx: [f64; 3],
    u_ax: [f64; 3],
    v_ax: [f64; 3],
    d_unit: [f64; 3],
    s_min: f64,
    span: f64,
    l_dom: f64,
}

fn canonical_tri_normal(normal: [f64; 3]) -> [f64; 3] {
    if normal[2] < -1e-12
        || (normal[2].abs() <= 1e-12 && (normal[1] < -1e-12 || (normal[1].abs() <= 1e-12 && normal[0] < 0.0)))
    {
        scale3(-1.0, normal)
    } else {
        normal
    }
}

impl ExcitationFluxoidMeta {
    /// Exact Whitney edge cochains `[c_01, c_02, c_12]` (ordered `v0 < v1 < v2`) representing the
    /// unit-fluxoid cohomology 1-form on sheet triangles without quadrature projection locking.
    fn edge_cochains(&self, v0: u32, v1: u32, v2: u32, pt3: &impl Fn(u32) -> [f64; 3]) -> [f64; 3] {
        if self.is_flux_loop {
            let theta = |v: u32| -> f64 {
                let rel = sub3(pt3(v), self.cx);
                dot3(rel, self.v_ax).atan2(dot3(rel, self.u_ax))
            };
            let wrap_frac = |a: f64, b: f64| -> f64 {
                let mut d = b - a;
                while d >= PI {
                    d -= 2.0 * PI;
                }
                while d < -PI {
                    d += 2.0 * PI;
                }
                d / (2.0 * PI)
            };
            let (t0, t1, t2) = (theta(v0), theta(v1), theta(v2));
            [wrap_frac(t0, t1), wrap_frac(t0, t2), wrap_frac(t1, t2)]
        } else {
            let alpha = |v: u32| -> f64 {
                if self.span > 1e-12 {
                    (dot3(pt3(v), self.d_unit) - self.s_min) / self.span
                } else {
                    dot3(pt3(v), self.d_unit) / self.l_dom
                }
            };
            let (a0, a1, a2) = (alpha(v0), alpha(v1), alpha(v2));
            [a1 - a0, a2 - a0, a2 - a1]
        }
    }
}

fn invert_spd_matrix(a: &[f64], n: usize) -> Result<Vec<Vec<f64>>, String> {
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut l = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.5 * (a[i * n + j] + a[j * n + i]);
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if !sum.is_finite() || sum <= 0.0 {
                    return Err(format!("inductance reluctance matrix is not positive definite at index {i}: {sum:e}"));
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    let mut inv = vec![vec![0.0_f64; n]; n];
    for col in 0..n {
        let mut y = vec![0.0_f64; n];
        for i in 0..n {
            let rhs = if i == col { 1.0 } else { 0.0 };
            let mut sum = rhs;
            for k in 0..i {
                sum -= l[i][k] * y[k];
            }
            y[i] = sum / l[i][i];
        }
        for i in (0..n).rev() {
            let mut sum = y[i];
            for k in (i + 1)..n {
                sum -= l[k][i] * inv[k][col];
            }
            inv[i][col] = sum / l[i][i];
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let sym = 0.5 * (inv[i][j] + inv[j][i]);
            inv[i][j] = sym;
            inv[j][i] = sym;
        }
    }
    Ok(inv)
}

type Mixed3dFluxoidRhs = (Vec<ExcitationFluxoidMeta>, Vec<Vec<f64>>, Vec<Vec<f64>>);

#[allow(clippy::too_many_arguments)]
fn build_mixed_3d_fluxoid_rhs<F, P, E>(
    excitations: &[Excitation],
    dedup_sheets: &BTreeMap<[u32; 3], f64>,
    vert_to_sidx: &BTreeMap<u32, usize>,
    facets_for_labels: F,
    pt3: P,
    edge_of: E,
    n_edges: usize,
    n_total: usize,
    l_dom: f64,
    fixed_set: &mut BTreeSet<usize>,
) -> Result<Mixed3dFluxoidRhs, String>
where
    F: Fn(&[i32]) -> Vec<[u32; 3]>,
    P: Fn(u32) -> [f64; 3],
    E: Fn(u32, u32) -> Option<usize>,
{
    let n_exc = excitations.len();
    let qtri = quad_rule(2, 2);
    let mut metas = Vec::with_capacity(n_exc);

    for ex in excitations {
        let mut ex_tris: BTreeSet<[u32; 3]> = BTreeSet::new();
        for fv in facets_for_labels(&ex.labels) {
            let mut key = fv;
            key.sort_unstable();
            ex_tris.insert(key);
        }
        if ex_tris.is_empty() {
            return Err(format!("excitation '{}' has no boundary facets", ex.name));
        }

        let mut cx = [0.0; 3];
        let mut norm_acc = [0.0; 3];
        let mut total_area = 0.0;
        for &[v0, v1, v2] in &ex_tris {
            let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
            let (area, normal, _) = triangle_surface_gradients(p0, p1, p2);
            let cen = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
            let n_up = canonical_tri_normal(normal);
            for k in 0..3 {
                cx[k] += area * cen[k];
                norm_acc[k] += area * n_up[k];
            }
            total_area += area;
        }
        if total_area <= 1e-30 {
            return Err(format!("excitation '{}' has zero facet area", ex.name));
        }
        for k in 0..3 {
            cx[k] /= total_area;
        }
        let n_acc_len = dot3(norm_acc, norm_acc).sqrt();
        let loop_normal = if n_acc_len > 1e-30 { scale3(1.0 / n_acc_len, norm_acc) } else { [0.0, 0.0, 1.0] };
        let seed = if loop_normal[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u_raw = cross3(seed, loop_normal);
        let u_ax = scale3(1.0 / dot3(u_raw, u_raw).sqrt().max(1e-30), u_raw);
        let v_ax = cross3(loop_normal, u_ax);

        let d_norm = dot3(ex.direction, ex.direction).sqrt().max(1e-30);
        let d_unit = scale3(1.0 / d_norm, ex.direction);
        let (mut s_min, mut s_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &[v0, v1, v2] in &ex_tris {
            for v in [v0, v1, v2] {
                let s = dot3(pt3(v), d_unit);
                s_min = s_min.min(s);
                s_max = s_max.max(s);
            }
        }
        let span = (s_max - s_min).max(0.0);

        if !ex.is_flux_loop && span > 1e-12 {
            let end_tol = 0.05 * span;
            for &[v0, v1, v2] in &ex_tris {
                for v in [v0, v1, v2] {
                    if let Some(&sidx) = vert_to_sidx.get(&v) {
                        let s = dot3(pt3(v), d_unit);
                        if (s - s_min).abs() <= end_tol || (s_max - s).abs() <= end_tol {
                            fixed_set.insert(n_edges + sidx);
                        }
                    }
                }
            }
        }

        let meta = ExcitationFluxoidMeta {
            tris: ex_tris,
            is_flux_loop: ex.is_flux_loop,
            cx,
            u_ax,
            v_ax,
            d_unit,
            s_min,
            span,
            l_dom,
        };

        metas.push(meta);
    }

    let mut rhs_u_list = vec![vec![0.0; n_total]; n_exc];
    let mut c_flux_list = vec![vec![0.0; n_total]; n_exc];

    for (p, meta) in metas.iter().enumerate() {
        let rhs = &mut rhs_u_list[p];
        let c_vec = &mut c_flux_list[p];
        for &key in &meta.tris {
            let [v0, v1, v2] = key;
            let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
            let (area, _normal, g_surf) = triangle_surface_gradients(p0, p1, p2);
            if area <= 1e-30 {
                continue;
            }
            let lk = dedup_sheets.get(&key).copied().unwrap_or(1e-6 * MU0 * l_dom);
            let (Some(e01), Some(e02), Some(e12)) = (edge_of(v0, v1), edge_of(v0, v2), edge_of(v1, v2)) else {
                continue;
            };
            let c_e = meta.edge_cochains(v0, v1, v2, &pt3);
            c_vec[e01] = c_e[0];
            c_vec[e02] = c_e[1];
            c_vec[e12] = c_e[2];
            let w_coeff = (MU0 / lk) * area;

            for q in 0..qtri.w.len() {
                let l = &qtri.lam[q];
                let w01 = sub3(scale3(l[0], g_surf[1]), scale3(l[1], g_surf[0]));
                let w02 = sub3(scale3(l[0], g_surf[2]), scale3(l[2], g_surf[0]));
                let w12 = sub3(scale3(l[1], g_surf[2]), scale3(l[2], g_surf[1]));
                let mut h_p = [0.0; 3];
                for k in 0..3 {
                    h_p[k] = c_e[0] * w01[k] + c_e[1] * w02[k] + c_e[2] * w12[k];
                }
                let wq = w_coeff * qtri.w[q];

                rhs[e01] += wq * dot3(h_p, w01);
                rhs[e02] += wq * dot3(h_p, w02);
                rhs[e12] += wq * dot3(h_p, w12);

                if let (Some(&i0), Some(&i1), Some(&i2)) =
                    (vert_to_sidx.get(&v0), vert_to_sidx.get(&v1), vert_to_sidx.get(&v2))
                {
                    rhs[n_edges + i0] -= wq * dot3(h_p, g_surf[0]);
                    rhs[n_edges + i1] -= wq * dot3(h_p, g_surf[1]);
                    rhs[n_edges + i2] -= wq * dot3(h_p, g_surf[2]);
                }
            }
        }
    }

    Ok((metas, rhs_u_list, c_flux_list))
}

type Mixed3dFluxoidResult = (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<String>);

#[allow(clippy::too_many_arguments)]
fn solve_mixed_3d_fluxoids_and_inductance<P, E>(
    excitations: &[Excitation],
    metas: &[ExcitationFluxoidMeta],
    dedup_sheets: &BTreeMap<[u32; 3], f64>,
    vert_to_sidx: &BTreeMap<u32, usize>,
    k_vr: &Csr,
    k_mag: &Csr,
    rhs_u_list: &[Vec<f64>],
    c_flux_list: &[Vec<f64>],
    fixed: &[usize],
    pt3: P,
    edge_of: E,
    n_edges: usize,
    n_sverts: usize,
) -> Result<Mixed3dFluxoidResult, String>
where
    P: Fn(u32) -> [f64; 3],
    E: Fn(u32, u32) -> Option<usize>,
{
    let n_exc = excitations.len();
    let n_total = k_mag.nrows();
    let free = free_dofs(n_total, fixed);
    let n_free = free.len();

    // Shift edge DOFs by `c_p` (zeroed on any fixed DOF) so we solve directly for the
    // kinetic residual `w_p = u_p - c_p` without catastrophic cancellation when `L_k` is small.
    let mut c_masked_list = Vec::with_capacity(n_exc);
    let mut rhs_w_free_list = Vec::with_capacity(n_exc);
    for p in 0..n_exc {
        let mut c_m = c_flux_list[p].clone();
        for &d in fixed {
            c_m[d] = 0.0;
        }
        let kc = k_mag.matvec(&c_m);
        let rhs_w_free: Vec<f64> = free.iter().map(|&i| rhs_u_list[p][i] - kc[i]).collect();
        c_masked_list.push(c_m);
        rhs_w_free_list.push(rhs_w_free);
    }

    let kff = k_mag.submatrix(&free, &free);
    let solver = DirectSolver::new(&kff, true)?;
    let w_sols_free = solver.solve_many(&rhs_w_free_list);
    let log = vec![format!(
        "3D mixed_3d London magnetostatics: 1 factorization ({n_free} free dofs: {n_edges} Nedelec edges + {n_sverts} sheet phi_s nodes on {} deduplicated triangles) served {n_exc} excitation(s)",
        dedup_sheets.len()
    )];

    let mut w_flux = vec![vec![0.0; n_total]; n_exc];
    let mut u_flux = Vec::with_capacity(n_exc);
    for p in 0..n_exc {
        for (a, &i) in free.iter().enumerate() {
            w_flux[p][i] = w_sols_free[p][a];
        }
        let mut u_p = c_masked_list[p].clone();
        for d in 0..n_total {
            u_p[d] += w_flux[p][d];
        }
        u_flux.push(u_p);
    }

    let k_vr_u: Vec<Vec<f64>> = u_flux.iter().map(|u| k_vr.matvec(u)).collect();
    let mut l_inv = vec![0.0_f64; n_exc * n_exc];

    for i in 0..n_exc {
        for j in i..n_exc {
            let vij: f64 = u_flux[i].iter().zip(&k_vr_u[j]).map(|(a, b)| a * b).sum();
            let vji: f64 = if i == j { vij } else { u_flux[j].iter().zip(&k_vr_u[i]).map(|(a, b)| a * b).sum() };
            let val = 0.5 * (vij + vji) / MU0;
            l_inv[i * n_exc + j] += val;
            if j != i {
                l_inv[j * n_exc + i] += val;
            }
        }
    }

    let qtri = quad_rule(2, 2);
    let mut r_q = vec![[0.0_f64; 3]; n_exc];
    for (&key, &lk) in dedup_sheets {
        let [v0, v1, v2] = key;
        let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
        let (area, _normal, g_surf) = triangle_surface_gradients(p0, p1, p2);
        if area <= 1e-30 {
            continue;
        }
        let (Some(e01), Some(e02), Some(e12)) = (edge_of(v0, v1), edge_of(v0, v2), edge_of(v1, v2)) else {
            continue;
        };
        let s0 = n_edges + vert_to_sidx[&v0];
        let s1 = n_edges + vert_to_sidx[&v1];
        let s2 = n_edges + vert_to_sidx[&v2];
        let w_base = area / lk;
        let edge_res_coefs: Vec<[f64; 3]> = (0..n_exc)
            .map(|i| {
                let src = if metas[i].tris.contains(&key) { &w_flux[i] } else { &u_flux[i] };
                [src[e01], src[e02], src[e12]]
            })
            .collect();

        for q in 0..qtri.w.len() {
            let l = &qtri.lam[q];
            let w01 = sub3(scale3(l[0], g_surf[1]), scale3(l[1], g_surf[0]));
            let w02 = sub3(scale3(l[0], g_surf[2]), scale3(l[2], g_surf[0]));
            let w12 = sub3(scale3(l[1], g_surf[2]), scale3(l[2], g_surf[1]));

            for i in 0..n_exc {
                let w = &w_flux[i];
                let ec = edge_res_coefs[i];
                let mut res = [0.0; 3];
                for k in 0..3 {
                    res[k] = ec[0] * w01[k] + ec[1] * w02[k] + ec[2] * w12[k]
                        - w[s0] * g_surf[0][k]
                        - w[s1] * g_surf[1][k]
                        - w[s2] * g_surf[2][k];
                }
                r_q[i] = res;
            }

            let wq = w_base * qtri.w[q];
            for i in 0..n_exc {
                for j in i..n_exc {
                    let contrib = wq * dot3(r_q[i], r_q[j]);
                    l_inv[i * n_exc + j] += contrib;
                    if j != i {
                        l_inv[j * n_exc + i] += contrib;
                    }
                }
            }
        }
    }

    let l_mat = invert_spd_matrix(&l_inv, n_exc)?;
    let mut u_all = vec![vec![0.0; n_total]; n_exc];
    for p in 0..n_exc {
        if excitations[p].is_flux_loop {
            let amp = excitations[p].amplitude;
            for d in 0..n_total {
                u_all[p][d] = amp * u_flux[p][d];
            }
        } else {
            let amp = excitations[p].amplitude;
            for q in 0..n_exc {
                let weight = l_mat[q][p] * amp;
                for d in 0..n_total {
                    u_all[p][d] += weight * u_flux[q][d];
                }
            }
        }
    }

    Ok((l_mat, u_all, log))
}

fn mixed_3d_fields(mesh: &Discretisation, order: usize, a_edges: &[f64]) -> (Vec<Field>, CellFields) {
    let geo = geometry(mesh);
    let lamn = local_bary(3, order);
    let nloc = lamn.len();
    let nc = mesh.n_cells();
    let mut a_mag_cell = vec![0.0; nc * nloc];
    let mut b_cell = vec![0.0; nc * nloc * 3];

    let mut gs = [[0.0; 3]; 4];
    for c in 0..nc {
        let v = mesh.cell(c);
        let eids = mesh.cell_edge_ids(c);
        for b in 0..4 {
            for k in 0..3 {
                gs[b][k] = geo.g[c * 12 + b * 3 + k];
            }
        }
        for (a, l) in lamn.iter().enumerate() {
            let lam = [l[0], l[1], l[2], l[3]];
            let (w, cu) = whitney_tet(&lam, &gs, v);
            let mut aval = [0.0; 3];
            let mut bval = [0.0; 3];
            for s in 0..6 {
                let coef = a_edges[eids[s] as usize];
                for k in 0..3 {
                    aval[k] += coef * w[s][k];
                    bval[k] += coef * cu[s][k];
                }
            }
            let idx = c * nloc + a;
            a_mag_cell[idx] = dot3(aval, aval).sqrt();
            b_cell[3 * idx] = bval[0];
            b_cell[3 * idx + 1] = bval[1];
            b_cell[3 * idx + 2] = bval[2];
        }
    }

    let a_nodal = node_average(mesh, &a_mag_cell, 1, order);
    let b_nodal = node_average(mesh, &b_cell, 3, order);
    let bc = components(&b_nodal, 3);
    let names = femcore::section::component_names(3);
    let fields = vec![Field::scalar("A", "T m", a_nodal), Field::vector("B", "T", names, &bc)];
    let mut cf = CellFields::new(mesh, order);
    cf.push_real("A", "T m", a_mag_cell);
    cf.push_components("B", names, "T", &b_cell);
    (fields, cf)
}

fn solve_mixed_3d(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    reject_curved_simplex(mesh)?;
    if spec.has_infinite_element() || spec.bcs.iter().any(|b| matches!(b.kind, BcKind::FarField { .. })) {
        let cm = Mesh::from_solver_mesh(mesh).map_err(|e| e.to_string())?;
        return solve_mixed_3d_mixed(&cm, spec);
    }
    solve_mixed_3d_impl::<UniformSource>(mesh, spec, &UniformSource).map(|(r, _)| r)
}

/// Solves 3D Nedelec `mixed_3d` magnetostatics with a generic monomorphized [`SourceField<3>`],
/// returning both the [`SolveResult`] and raw Nedelec edge coefficients `a_edges` (`0..n_edges`).
pub fn solve_mixed_3d_edges_with_source<S: SourceField<3>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<(SolveResult, Vec<f64>), String> {
    reject_curved_simplex(mesh)?;
    if !source.is_uniform()
        && (spec.has_infinite_element() || spec.bcs.iter().any(|b| matches!(b.kind, BcKind::FarField { .. })))
    {
        return Err("mixed_3d SourceField MMS does not support infinite_element or far_field".into());
    }
    solve_mixed_3d_impl::<S>(mesh, spec, source)
}

/// Solves 3D Nedelec `mixed_3d` magnetostatics with a generic monomorphized [`SourceField<3>`].
pub fn solve_mixed_3d_with_source<S: SourceField<3>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<SolveResult, String> {
    solve_mixed_3d_edges_with_source::<S>(mesh, spec, source).map(|(r, _)| r)
}

fn solve_mixed_3d_impl<S: SourceField<3>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<(SolveResult, Vec<f64>), String> {
    let t0 = Instant::now();
    if mesh.dim != 3 {
        return Err("the 'mixed_3d' magnetostatic formulation requires a 3D mesh".into());
    }
    let order = spec.order;
    check_boundary_labels(
        mesh,
        spec.bcs
            .iter()
            .filter(|b| !matches!(b.kind, BcKind::InfiniteElement(_)))
            .flat_map(|b| b.labels.iter().copied()),
    )?;
    let excitations = if source.is_uniform() { collect_excitations(spec)? } else { vec![] };
    let mu_r = spec.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    if mu_r.iter().any(|&m| m <= 0.0) {
        return Err("material property 'mu_r' must be positive".into());
    }

    let l_dom = domain_diameter_3d(mesh.n_points(), |i| {
        let p = mesh.point(i);
        [p[0], p[1], p[2]]
    });
    let eps_reg = 1e-6 / (l_dom * l_dom);

    // Collect LondonSheet and Ground BCs
    let mut sheet_lk_by_label: BTreeMap<i32, f64> = BTreeMap::new();
    let mut ground_labels: Vec<i32> = Vec::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::LondonSheet { l_k_sq } => {
                for &lab in &bc.labels {
                    sheet_lk_by_label.insert(lab, *l_k_sq);
                }
            }
            BcKind::Ground => {
                for &lab in &bc.labels {
                    if !ground_labels.contains(&lab) {
                        ground_labels.push(lab);
                    }
                }
            }
            _ => {}
        }
    }

    // Single-facet deduplication of london_sheet triangles by sorted [v0, v1, v2]
    let mut dedup_sheets: BTreeMap<[u32; 3], f64> = BTreeMap::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if let Some(&lk) = sheet_lk_by_label.get(&lab) {
            let fv = mesh.facet(f);
            let mut key = [fv[0], fv[1], fv[2]];
            key.sort_unstable();
            let lk_eff = if lk > 0.0 { lk } else { 1e-6 * MU0 * l_dom };
            dedup_sheets.entry(key).or_insert(lk_eff);
        }
    }

    // Map unique sheet vertices to surface scalar phi_s DOFs [n_edges .. n_edges + n_sverts]
    let n_edges = mesh.n_edges();
    let mut sheet_verts_set = BTreeSet::new();
    let mut sheet_edges_set = BTreeSet::new();
    for &[v0, v1, v2] in dedup_sheets.keys() {
        sheet_verts_set.insert(v0);
        sheet_verts_set.insert(v1);
        sheet_verts_set.insert(v2);
        sheet_edges_set.insert(mesh.edge_index(v0, v1));
        sheet_edges_set.insert(mesh.edge_index(v0, v2));
        sheet_edges_set.insert(mesh.edge_index(v1, v2));
    }
    let sheet_verts: Vec<u32> = sheet_verts_set.into_iter().collect();
    let vert_to_sidx: BTreeMap<u32, usize> = sheet_verts.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let n_sverts = sheet_verts.len();
    let n_total = n_edges + n_sverts;

    // 1. Assemble 3D volume curl-curl + gauge regularization matrix over Nedelec edges 0..n_edges
    let geo = geometry(mesh);
    let qtet = quad_rule(3, 2);
    let nc = mesh.n_cells();
    let k_vol = assemble_parallel(nc, n_total, n_total, |range, coo| {
        let mut gs = [[0.0; 3]; 4];
        let mut ke = [0.0; 36];
        for c in range {
            let v = mesh.cell(c);
            let eids = mesh.cell_edge_ids(c);
            let dofs: [usize; 6] = [
                eids[0] as usize,
                eids[1] as usize,
                eids[2] as usize,
                eids[3] as usize,
                eids[4] as usize,
                eids[5] as usize,
            ];
            for b in 0..4 {
                for k in 0..3 {
                    gs[b][k] = geo.g[c * 12 + b * 3 + k];
                }
            }
            ke.fill(0.0);
            let inv_mur = 1.0 / mu_r[c];
            let vol = geo.vol[c];
            for q in 0..qtet.w.len() {
                let lam = [qtet.lam[q][0], qtet.lam[q][1], qtet.lam[q][2], qtet.lam[q][3]];
                let (w, cu) = whitney_tet(&lam, &gs, v);
                let wq = qtet.w[q] * vol;
                for i in 0..6 {
                    for j in 0..6 {
                        let curl_dot = dot3(cu[i], cu[j]);
                        let mass_dot = dot3(w[i], w[j]);
                        ke[i * 6 + j] += wq * (inv_mur * curl_dot + eps_reg * mass_dot);
                    }
                }
            }
            coo.push_block(&dofs, &dofs, &ke);
        }
    });

    // 2. Assemble 2D surface kinetic coupling term on deduplicated london_sheet triangles:
    //    int_{Gamma_SC} (mu_0 / L_k,sq) (A_|| - grad_|| phi_s) . (v_|| - grad_|| psi_s) dS
    let qtri = quad_rule(2, 2);
    let mut coo_surf = Coo::new(n_total, n_total);
    let pt3 = |v: u32| -> [f64; 3] {
        let p = mesh.point(v as usize);
        [p[0], p[1], p[2]]
    };

    for (&[v0, v1, v2], &lk) in &dedup_sheets {
        let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
        let (area, _normal, g_surf) = triangle_surface_gradients(p0, p1, p2);
        if area <= 1e-30 {
            continue;
        }
        let e01 = mesh.edge_index(v0, v1);
        let e02 = mesh.edge_index(v0, v2);
        let e12 = mesh.edge_index(v1, v2);
        let s0 = n_edges + vert_to_sidx[&v0];
        let s1 = n_edges + vert_to_sidx[&v1];
        let s2 = n_edges + vert_to_sidx[&v2];
        let dofs = [e01, e02, e12, s0, s1, s2];
        let coeff = (MU0 / lk) * area;

        let mut ke = [0.0; 36];
        for q in 0..qtri.w.len() {
            let l = &qtri.lam[q];
            let w01 = sub3(scale3(l[0], g_surf[1]), scale3(l[1], g_surf[0]));
            let w02 = sub3(scale3(l[0], g_surf[2]), scale3(l[2], g_surf[0]));
            let w12 = sub3(scale3(l[1], g_surf[2]), scale3(l[2], g_surf[1]));
            let basis = [w01, w02, w12, scale3(-1.0, g_surf[0]), scale3(-1.0, g_surf[1]), scale3(-1.0, g_surf[2])];
            let wq = coeff * qtri.w[q];
            for i in 0..6 {
                for j in 0..6 {
                    ke[i * 6 + j] += wq * dot3(basis[i], basis[j]);
                }
            }
        }
        coo_surf.push_block(&dofs, &dofs, &ke);
    }

    let k_mag = k_vol.add_scaled(&coo_surf.to_csr(), 1.0);

    // 3. Pin one reference phi_s vertex per connected sheet component (plus terminal end-pins)
    let mut fixed_set = BTreeSet::new();
    if n_sverts > 0 {
        let mut parent: Vec<usize> = (0..n_sverts).collect();
        fn find(p: &mut [usize], mut i: usize) -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        }
        for &[v0, v1, v2] in dedup_sheets.keys() {
            let (i0, i1, i2) = (vert_to_sidx[&v0], vert_to_sidx[&v1], vert_to_sidx[&v2]);
            for (a, b) in [(i0, i1), (i1, i2), (i0, i2)] {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        let mut seen_roots = BTreeSet::new();
        for i in 0..n_sverts {
            let r = find(&mut parent, i);
            if seen_roots.insert(r) {
                fixed_set.insert(n_edges + r);
            }
        }
    }

    // 3b & 4. Assemble fixed DOFs & RHS depending on uniform excitations vs monomorphized SourceField<3>
    if !source.is_uniform() {
        let mut dirichlet_edge_acc: BTreeMap<usize, (f64, usize)> = BTreeMap::new();
        let qedge = quad_seg_high();
        for f in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f];
            if ground_labels.contains(&lab) {
                let fv = mesh.facet(f);
                for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                    let (v_lo, v_hi) = (a.min(b), a.max(b));
                    let e = mesh.edge_index(v_lo, v_hi);
                    if sheet_edges_set.contains(&e) {
                        continue;
                    }
                    let pa = pt3(v_lo);
                    let pb = pt3(v_hi);
                    let dl = sub3(pb, pa);
                    let mut val = 0.0;
                    for (q, lam) in qedge.lam.iter().enumerate() {
                        let xq = [
                            lam[0] * pa[0] + lam[1] * pb[0],
                            lam[0] * pa[1] + lam[1] * pb[1],
                            lam[0] * pa[2] + lam[1] * pb[2],
                        ];
                        let a_vec = source.dirichlet_re(xq, 0.0, lab as u32);
                        val += qedge.w[q] * dot3(a_vec, dl);
                    }
                    let entry = dirichlet_edge_acc.entry(e).or_insert((0.0, 0));
                    entry.0 += val;
                    entry.1 += 1;
                }
            }
        }
        let mut fixed_map: BTreeMap<usize, f64> = BTreeMap::new();
        for d in fixed_set {
            fixed_map.insert(d, 0.0);
        }
        for (e, (sum, cnt)) in dirichlet_edge_acc {
            fixed_map.insert(e, sum / cnt as f64);
        }
        let fixed: Vec<usize> = fixed_map.keys().copied().collect();
        let fixed_vals: Vec<f64> = fixed_map.values().copied().collect();
        let free = free_dofs(n_total, &fixed);
        let n_free = free.len();

        let mut rhs = vec![0.0; n_total];
        let q64 = quad_rule_collapsed(4);
        let mut gs = [[0.0; 3]; 4];
        for c in 0..nc {
            let v = mesh.cell(c);
            let eids = mesh.cell_edge_ids(c);
            let dom = mesh.cell_labels[c] as u32;
            for b in 0..4 {
                for k in 0..3 {
                    gs[b][k] = geo.g[c * 12 + b * 3 + k];
                }
            }
            let pts = [pt3(v[0]), pt3(v[1]), pt3(v[2]), pt3(v[3])];
            let vol = geo.vol[c];
            for q in 0..q64.w.len() {
                let lam = [q64.lam[q][0], q64.lam[q][1], q64.lam[q][2], q64.lam[q][3]];
                let mut xq = [0.0; 3];
                for b in 0..4 {
                    for k in 0..3 {
                        xq[k] += lam[b] * pts[b][k];
                    }
                }
                let (w, _) = whitney_tet(&lam, &gs, v);
                let j_vec = source.volume_rhs_re(xq, 0.0, dom);
                let wq = q64.w[q] * vol * MU0;
                for i in 0..6 {
                    rhs[eids[i] as usize] += wq * dot3(j_vec, w[i]);
                }
            }
        }

        let mut u = vec![0.0; n_total];
        for (&d, &val) in fixed.iter().zip(&fixed_vals) {
            u[d] = val;
        }
        if n_free > 0 {
            let kff = k_mag.submatrix(&free, &free);
            let kfx = k_mag.submatrix(&free, &fixed);
            let kx = if fixed.is_empty() { vec![0.0; n_free] } else { kfx.matvec(&fixed_vals) };
            let rhs_free: Vec<f64> = free.iter().enumerate().map(|(a, &i)| rhs[i] - kx[a]).collect();
            let solver = DirectSolver::new(&kff, true)?;
            let sol = solver.solve(&rhs_free);
            for (a, &i) in free.iter().enumerate() {
                u[i] = sol[a];
            }
        }
        let a_edges = u[..n_edges].to_vec();
        let (fields, cell_fields) = mixed_3d_fields(mesh, order, &a_edges);
        let bn = norm_of(&components(&fields[1].values, 3));
        let res = SolveResult {
            physics: "dfs".into(),
            study: spec.study.kind().into(),
            order: order as u8,
            dim: 3,
            solutions: vec![Solution {
                label: "source_field".into(),
                fields,
                cell_fields: cell_fields_of(&cell_fields),
                ..Default::default()
            }],
            scalars: vec![("max |B| (T)".into(), Scalar::Num(minmax(&bn).1))],
            log: vec![format!("3D mixed_3d SourceField MMS ({n_free} free Nedelec DOFs)")],
            seconds: t0.elapsed().as_secs_f64(),
            dofs: n_total,
        };
        return Ok((res, a_edges));
    }

    // Exterior PEC / Ground Nedelec edge constraints (never fixing edges on london_sheet facets)
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if ground_labels.contains(&lab) {
            let fv = mesh.facet(f);
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                let e = mesh.edge_index(a, b);
                if !sheet_edges_set.contains(&e) {
                    fixed_set.insert(e);
                }
            }
        }
    }

    // 4. Build unit-fluxoid cohomology right-hand sides & solve non-cancelling residual system
    let (metas, rhs_u_list, c_flux_list) = build_mixed_3d_fluxoid_rhs(
        &excitations,
        &dedup_sheets,
        &vert_to_sidx,
        |labs| {
            facets_with(mesh, labs)
                .map(|f| {
                    let fv = mesh.facet(f);
                    [fv[0], fv[1], fv[2]]
                })
                .collect()
        },
        pt3,
        |a, b| Some(mesh.edge_index(a, b)),
        n_edges,
        n_total,
        l_dom,
        &mut fixed_set,
    )?;

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let (l_mat, u_all, log) = solve_mixed_3d_fluxoids_and_inductance(
        &excitations,
        &metas,
        &dedup_sheets,
        &vert_to_sidx,
        &k_vol,
        &k_mag,
        &rhs_u_list,
        &c_flux_list,
        &fixed,
        pt3,
        |a, b| Some(mesh.edge_index(a, b)),
        n_edges,
        n_sverts,
    )?;

    let n_exc = excitations.len();
    let mut solutions = Vec::with_capacity(n_exc);
    let mut scalars = Vec::new();
    for p in 0..n_exc {
        let a_edges = &u_all[p][..n_edges];
        let (fields, cell_fields) = mixed_3d_fields(mesh, order, a_edges);
        if p == 0 {
            let bn = norm_of(&components(&fields[1].values, 3));
            let w_m = 0.5 * l_mat[0][0] * excitations[0].amplitude.powi(2);
            scalars.push(("stored magnetic energy W_m (J)".into(), Scalar::Num(w_m)));
            scalars.push(("max |B| (T)".into(), Scalar::Num(minmax(&bn).1)));
        }
        solutions.push(Solution {
            label: format!("terminal {}", excitations[p].name),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(l_mat[p][p]),
            ..Default::default()
        });
    }

    push_inductance_scalars(&mut scalars, &excitations, &l_mat, "mixed_3d");
    let _ = LinearSolver::Auto;

    let first_edges = u_all[0][..n_edges].to_vec();
    Ok((
        SolveResult {
            physics: "dfs".into(),
            study: spec.study.kind().into(),
            order: order as u8,
            dim: 3,
            solutions,
            scalars,
            log,
            seconds: t0.elapsed().as_secs_f64(),
            dofs: n_total,
        },
        first_edges,
    ))
}

/// Solves a magnetostatic problem using the requested or auto-selected formulation.
pub fn solve_magnetostatics(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    reject_curved_simplex(mesh)?;
    let Study::Magnetostatic { formulation } = &spec.study else {
        return Err("solve_magnetostatics requires a magnetostatic study".into());
    };
    let lower = formulation.trim().to_ascii_lowercase();
    let form = match lower.as_str() {
        "" | "auto" => {
            if mesh.dim == 2 {
                "cross_section"
            } else {
                "mixed_3d"
            }
        }
        other => other,
    };
    match form {
        "cross_section" => solve_cross_section(mesh, spec),
        "planar_stream" => solve_planar_stream(mesh, spec),
        "mixed_3d" => solve_mixed_3d(mesh, spec),
        other => Err(format!(
            "unknown magnetostatic formulation '{other}' (expected 'auto', 'cross_section', 'planar_stream', or 'mixed_3d')"
        )),
    }
}

// ------------------------------------------------------------------ Native Mixed-Cell & Open-Boundary Magnetostatics

fn cross_section_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    order: usize,
    az_mixed: &[f64],
) -> (Vec<Field>, CellFields) {
    let az_simplex = proj.project_lagrange_blocks(mixed, dof_map, az_mixed, 1, order);
    let ones = vec![1.0; mixed.n_cells()];
    let grad_az = sample_mixed_gradients_on_simplex(mixed, dof_map, proj, az_mixed, None, &ones, order);
    let mut bcell = vec![0.0; grad_az.len()];
    for p in 0..(grad_az.len() / 2) {
        let dax = grad_az[2 * p];
        let daz = grad_az[2 * p + 1];
        bcell[2 * p] = -daz;
        bcell[2 * p + 1] = dax;
    }
    let b_avg = node_average(simplex, &bcell, 2, order);
    let bc = components(&b_avg, 2);
    let names = femcore::section::component_names(2);
    let fields = vec![Field::scalar("A", "T m", az_simplex.clone()), Field::vector("B", "T", names, &bc)];
    let mut cf = CellFields::new(simplex, order);
    cf.push_real("A", "T m", cell_values_at_nodes(simplex, &az_simplex, order));
    cf.push_components("B", names, "T", &bcell);
    (fields, cf)
}

fn solve_cross_section_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String> {
    if mesh.is_pure_simplex() && !spec.has_infinite_element() {
        let d = discretise(mesh).map_err(|e| e.to_string())?;
        return solve_cross_section(&d, spec);
    }
    let t0 = Instant::now();
    if mesh.dim() != 2 {
        return Err("the 'cross_section' magnetostatic formulation requires a 2D mesh".into());
    }
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
    let excitations = collect_excitations(spec)?;
    let ie_map = resolve_ie_domains_mixed(&mixed, spec)?;
    let mu_r = cell_values_mixed(&spec.materials, &mixed, "mu_r", Some(1.0))?;
    if mu_r.iter().any(|&m| m <= 0.0) {
        return Err("material property 'mu_r' must be positive".into());
    }
    let inv_mur: Vec<f64> = mu_r.iter().map(|&m| 1.0 / m).collect();
    let coeff: Vec<f64> = inv_mur.iter().map(|&im| im / MU0).collect();
    let mut k = assemble_laplace_ie_mixed(&mixed, &dof_map, &coeff, None, &ie_map);

    let mut sheet_lk: BTreeMap<i32, f64> = BTreeMap::new();
    let mut ground_labels: Vec<i32> = Vec::new();
    let mut ff_label_set: BTreeSet<i32> = BTreeSet::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::LondonSheet { l_k_sq } => {
                for &lab in &bc.labels {
                    sheet_lk.insert(lab, *l_k_sq);
                }
            }
            BcKind::Ground => {
                for &lab in &bc.labels {
                    if !ground_labels.contains(&lab) {
                        ground_labels.push(lab);
                    }
                }
            }
            BcKind::FarField { .. } => {
                for &lab in &bc.labels {
                    ff_label_set.insert(lab);
                }
            }
            _ => {}
        }
    }

    let term_label_set: BTreeSet<i32> = excitations.iter().flat_map(|e| e.labels.iter().copied()).collect();
    let ref_labels: Vec<i32> = term_label_set.iter().copied().collect();
    let mut has_far_field = false;
    let mur_perm = FarFieldPermittivity::Scalar(&inv_mur);
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !ground_labels.contains(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid_mixed(&mixed, &ref_labels));
                let k_far = assemble_farfield_robin_mixed(
                    &mixed,
                    &dof_map,
                    &ff_labs,
                    *ff_order,
                    *decay_power,
                    c,
                    1.0 / MU0,
                    mur_perm,
                );
                k = k.add_scaled(&k_far, 1.0);
                has_far_field = true;
            }
        }
    }

    let (ie_rat_facets, ie_exp_facets) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
    let mut ie_exp_facet_set = BTreeSet::new();
    if !ie_exp_facets.is_empty() {
        for &(f_idx, _, _) in &ie_exp_facets {
            ie_exp_facet_set.insert(f_idx);
        }
        let k_ie_rob = assemble_ie_outer_robin_mixed(&mixed, &dof_map, &ie_exp_facets, &ie_map, 1.0 / MU0, mur_perm);
        k = k.add_scaled(&k_ie_rob, 1.0);
        has_far_field = true;
    }

    for (&lab, &l_k_sq) in &sheet_lk {
        if !term_label_set.contains(&lab) {
            if l_k_sq > 0.0 {
                k = k.add_scaled(&assemble_boundary_mass_mixed(&mixed, &dof_map, &[lab], 1.0 / l_k_sq), 1.0);
            } else if !ground_labels.contains(&lab) {
                ground_labels.push(lab);
            }
        }
    }

    let mut fixed_set: BTreeSet<usize> = dof_map.boundary_dofs(&mixed, &ground_labels).into_iter().collect();
    if !ie_rat_facets.is_empty() {
        let mut fbuf = Vec::new();
        for &f_idx in &ie_rat_facets {
            dof_map.facet_dofs(&mixed, f_idx, &mut fbuf);
            fixed_set.extend(fbuf.iter().copied());
        }
    }

    if fixed_set.is_empty() && !has_far_field {
        let mut fbuf = Vec::new();
        for f in 0..mixed.n_facets() {
            let lab = mixed.facet_labels[f];
            if mixed.exterior_facets[f]
                && !term_label_set.contains(&lab)
                && !sheet_lk.contains_key(&lab)
                && !ff_label_set.contains(&lab)
                && !ie_exp_facet_set.contains(&f)
            {
                dof_map.facet_dofs(&mixed, f, &mut fbuf);
                fixed_set.extend(fbuf.iter().copied());
            }
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    if fixed.is_empty() && !has_far_field {
        return Err(
            "2D cross_section magnetostatics needs an exterior, ground, or far_field reference boundary (A_z = 0)"
                .into(),
        );
    }

    let n = dof_map.n_dofs;
    let free = free_dofs(n, &fixed);
    let n_free = free.len();
    let n_exc = excitations.len();

    let mut rhs_full_list = Vec::with_capacity(n_exc);
    let mut term_kin_self = vec![0.0; n_exc];
    for (p, ex) in excitations.iter().enumerate() {
        let width = boundary_measure_mixed(&mixed, &ex.labels);
        if width <= 0.0 {
            return Err(format!("current_terminal '{}' has zero boundary measure", ex.name));
        }
        let j_surf = ex.amplitude / width;
        rhs_full_list.push(assemble_boundary_load_mixed(&mixed, &dof_map, &ex.labels, j_surf));

        let mut lk_int = 0.0;
        for (&lab, &lk) in &sheet_lk {
            if lk > 0.0 && ex.labels.contains(&lab) {
                lk_int += lk * boundary_measure_mixed(&mixed, &[lab]);
            }
        }
        term_kin_self[p] = lk_int / (width * width);
    }

    let kff = k.submatrix(&free, &free);
    let rhs_free_list: Vec<Vec<f64>> = rhs_full_list.iter().map(|r| free.iter().map(|&i| r[i]).collect()).collect();
    let solver = DirectSolver::new(&kff, true)?;
    let sols_free = solver.solve_many(&rhs_free_list);
    let log = vec![format!(
        "2D cross_section magnetostatics: 1 factorization ({n_free} free mixed-cell dofs) served {n_exc} excitation(s)"
    )];

    let mut az_all = vec![vec![0.0; n]; n_exc];
    for p in 0..n_exc {
        for (a, &i) in free.iter().enumerate() {
            az_all[p][i] = sols_free[p][a];
        }
    }

    let kaz_all: Vec<Vec<f64>> = az_all.iter().map(|a| k.matvec(a)).collect();
    let mut l_mat = vec![vec![0.0; n_exc]; n_exc];
    for i in 0..n_exc {
        for j in i..n_exc {
            let e_ij = az_all[i].iter().zip(&kaz_all[j]).map(|(a, b)| a * b).sum::<f64>();
            let e_ji = if i == j { e_ij } else { az_all[j].iter().zip(&kaz_all[i]).map(|(a, b)| a * b).sum::<f64>() };
            let mut lij = 0.5 * (e_ij + e_ji) / (excitations[i].amplitude * excitations[j].amplitude);
            if i == j {
                lij += term_kin_self[i];
            }
            l_mat[i][j] = lij;
            l_mat[j][i] = lij;
        }
    }

    let mut solutions = Vec::with_capacity(n_exc);
    let mut scalars = Vec::new();
    for p in 0..n_exc {
        let (fields, cell_fields) = cross_section_fields_mixed(&mixed, &dof_map, &simplex, &proj, order, &az_all[p]);
        if p == 0 {
            let bn = norm_of(&components(&fields[1].values, 2));
            let w_m = 0.5 * l_mat[0][0] * excitations[0].amplitude.powi(2);
            scalars.push(("stored magnetic energy W_m (J/m)".into(), Scalar::Num(w_m)));
            scalars.push(("max |B| (T)".into(), Scalar::Num(minmax(&bn).1)));
        }
        solutions.push(Solution {
            label: format!("terminal {}", excitations[p].name),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(l_mat[p][p]),
            ..Default::default()
        });
    }

    push_inductance_scalars(&mut scalars, &excitations, &l_mat, "cross_section");

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
    })
}

/// Assembles the 3D 1st-order asymptotic curl-Robin boundary stiffness
/// `K_{far,ab} = \int_{\Gamma} (\gamma s_\Omega / (\mu_r r)) (\hat{\mathbf{n}} \times \mathbf{w}_a) \cdot (\hat{\mathbf{n}} \times \mathbf{w}_b) dS`
/// (`\mu_0`-scaled to match `k_vol`) over 3D `Tri3` and `Quad4` exterior facets.
#[allow(clippy::too_many_arguments)]
fn assemble_curl_robin_3d_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    n_total: usize,
    mu_r: &[f64],
    facets_with_meta: &[(usize, usize, f64, [f64; 3], f64)], // (facet_idx, cell_idx, gamma, center, s_omega)
) -> Csr {
    if facets_with_meta.is_empty() {
        return Csr::zeros(n_total, n_total);
    }
    let mut coo = Coo::new(n_total, n_total);
    let mut scratch = NedEvalScratch::default();
    let mut w_out = Vec::new();
    let mut curl_out = Vec::new();
    let mut nw: Vec<[f64; 3]> = Vec::new();
    let mut ke = Vec::new();

    for &(f, c_f, gamma, center, s_omega) in facets_with_meta {
        let ct = mixed.cell_type(c_f);
        let vs = mixed.cell(c_f);
        let fv = mixed.facet(f);
        let ft = mixed.facet_type(f);
        let corner_refs = cell_corner_ref_coords(ct);

        let mut facet_corner_refs = Vec::with_capacity(fv.len());
        let mut all_found = true;
        for &v_id in fv {
            if let Some(pos) = vs.iter().position(|&x| x == v_id) {
                facet_corner_refs.push(corner_refs[pos]);
            } else {
                all_found = false;
                break;
            }
        }
        if !all_found {
            continue;
        }

        let iso = IsoparametricCell::from_mesh(mixed, c_f);
        let dofs = ned.cell_dofs(c_f);
        let nloc = dofs.len();
        ke.clear();
        ke.resize(nloc * nloc, 0.0);

        let rule = facet_rule_mixed(ft, if iso.curved { 4 } else { 2 });
        let inv_mur = 1.0 / mu_r[c_f];

        for (eta, &wq) in rule.pts.iter().zip(&rule.w) {
            let (x_q, ds_factor, n_hat) = eval_facet_jacobian(mixed, f, eta);
            if ds_factor <= 1e-30 {
                continue;
            }
            let r_vec = sub3(x_q, center);
            let r = dot3(r_vec, r_vec).sqrt();
            if r <= 1e-30 {
                continue;
            }

            let (n_f, _) = facet_lagrange_shape(ft, LagrangeSpace::Order1, eta);
            let mut xi_cell = [0.0; 3];
            for (k, &nk) in n_f.iter().enumerate() {
                xi_cell[0] += nk * facet_corner_refs[k][0];
                xi_cell[1] += nk * facet_corner_refs[k][1];
                xi_cell[2] += nk * facet_corner_refs[k][2];
            }

            eval_cell_basis_3d(mixed, ned, c_f, &iso, &xi_cell, &mut scratch, &mut w_out, &mut curl_out);
            nw.clear();
            for w_a in &w_out {
                nw.push(cross3(n_hat, *w_a));
            }

            let w = wq * ds_factor * (gamma * s_omega * inv_mur / r);
            for a in 0..nloc {
                for b in 0..nloc {
                    ke[a * nloc + b] += w * dot3(nw[a], nw[b]);
                }
            }
        }
        coo.push_block(dofs, dofs, &ke);
    }
    coo.to_csr()
}

fn mixed_3d_fields_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    order: usize,
    a_edges: &[f64],
    ie_map: &BTreeMap<i32, ResolvedIeDomain>,
) -> (Vec<Field>, CellFields) {
    let iso_cells: Vec<IsoparametricCell> =
        (0..mixed.n_cells()).map(|c| IsoparametricCell::from_mesh(mixed, c)).collect();
    let mut scratch = NedEvalScratch::default();
    let mut w_out = Vec::new();
    let mut curl_out = Vec::new();

    let samples: Vec<(f64, [f64; 3])> = proj.sample_simplex_cell_nodes(mixed, order, |mc, xi| {
        let jac_k = eval_cell_basis_3d(mixed, ned, mc, &iso_cells[mc], xi, &mut scratch, &mut w_out, &mut curl_out);
        let dofs = ned.cell_dofs(mc);
        let mut aval = [0.0; 3];
        let mut bval = [0.0; 3];
        for (s, &gd) in dofs.iter().enumerate() {
            let coef = a_edges[gd];
            for k in 0..3 {
                aval[k] += coef * w_out[s][k];
                bval[k] += coef * curl_out[s][k];
            }
        }
        if let Some(ie) = ie_map.get(&mixed.cell_labels[mc]) {
            let jac_ie = ie.eval_jacobian(jac_k.x, 3);
            let inv_det = 1.0 / jac_ie.det.re.abs().max(1e-30);
            let mut b_ie = [0.0; 3];
            for i in 0..3 {
                for m in 0..3 {
                    b_ie[i] += jac_ie.j[i][m].re * bval[m] * inv_det;
                }
            }
            bval = b_ie;
        }
        (dot3(aval, aval).sqrt(), bval)
    });

    let mut a_mag_cell = Vec::with_capacity(samples.len());
    let mut b_cell = Vec::with_capacity(samples.len() * 3);
    for (amag, b) in samples {
        a_mag_cell.push(amag);
        b_cell.extend_from_slice(&b);
    }

    let a_nodal = node_average(simplex, &a_mag_cell, 1, order);
    let b_nodal = node_average(simplex, &b_cell, 3, order);
    let bc = components(&b_nodal, 3);
    let names = femcore::section::component_names(3);
    let fields = vec![Field::scalar("A", "T m", a_nodal), Field::vector("B", "T", names, &bc)];
    let mut cf = CellFields::new(simplex, order);
    cf.push_real("A", "T m", a_mag_cell);
    cf.push_components("B", names, "T", &b_cell);
    (fields, cf)
}

fn solve_mixed_3d_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String> {
    if mesh.is_pure_simplex()
        && !spec.has_infinite_element()
        && !spec.bcs.iter().any(|b| matches!(b.kind, BcKind::FarField { .. }))
    {
        let d = discretise(mesh).map_err(|e| e.to_string())?;
        return solve_mixed_3d(&d, spec);
    }
    let t0 = Instant::now();
    if mesh.dim() != 3 {
        return Err("the 'mixed_3d' magnetostatic formulation requires a 3D mesh".into());
    }
    let order = spec.order;
    let mixed = discretise_mixed(mesh).map_err(|e| e.to_string())?;
    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, &mixed).map_err(|e| e.to_string())?;
    let ned = NedMixed::new(&mixed, 1);

    check_boundary_labels_mixed(
        &mixed,
        spec.bcs
            .iter()
            .filter(|b| !matches!(b.kind, BcKind::InfiniteElement(_)))
            .flat_map(|b| b.labels.iter().copied()),
    )?;
    let excitations = collect_excitations(spec)?;
    let ie_map = resolve_ie_domains_mixed(&mixed, spec)?;
    let mu_r = cell_values_mixed(&spec.materials, &mixed, "mu_r", Some(1.0))?;
    if mu_r.iter().any(|&m| m <= 0.0) {
        return Err("material property 'mu_r' must be positive".into());
    }

    let l_dom = domain_diameter_3d(mixed.n_points(), |i| mixed.point3(i));
    let eps_reg = 1e-6 / (l_dom * l_dom);

    let mut sheet_lk_by_label: BTreeMap<i32, f64> = BTreeMap::new();
    let mut ground_labels: Vec<i32> = Vec::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::LondonSheet { l_k_sq } => {
                for &lab in &bc.labels {
                    sheet_lk_by_label.insert(lab, *l_k_sq);
                }
            }
            BcKind::Ground => {
                for &lab in &bc.labels {
                    if !ground_labels.contains(&lab) {
                        ground_labels.push(lab);
                    }
                }
            }
            _ => {}
        }
    }

    // Single-facet deduplication of london_sheet triangles by sorted [v0, v1, v2]
    let mut dedup_sheets: BTreeMap<[u32; 3], f64> = BTreeMap::new();
    for f in 0..mixed.n_facets() {
        let lab = mixed.facet_labels[f];
        if let Some(&lk) = sheet_lk_by_label.get(&lab) {
            let fv = mixed.facet(f);
            if fv.len() == 3 {
                let mut key = [fv[0], fv[1], fv[2]];
                key.sort_unstable();
                let lk_eff = if lk > 0.0 { lk } else { 1e-6 * MU0 * l_dom };
                dedup_sheets.entry(key).or_insert(lk_eff);
            }
        }
    }

    let n_edges = ned.n_dofs;
    let mut sheet_verts_set = BTreeSet::new();
    let mut sheet_edges_set = BTreeSet::new();
    for &[v0, v1, v2] in dedup_sheets.keys() {
        sheet_verts_set.insert(v0);
        sheet_verts_set.insert(v1);
        sheet_verts_set.insert(v2);
        if let Some(e) = mixed.edge_index(v0, v1) {
            sheet_edges_set.insert(e);
        }
        if let Some(e) = mixed.edge_index(v0, v2) {
            sheet_edges_set.insert(e);
        }
        if let Some(e) = mixed.edge_index(v1, v2) {
            sheet_edges_set.insert(e);
        }
    }
    let sheet_verts: Vec<u32> = sheet_verts_set.into_iter().collect();
    let vert_to_sidx: BTreeMap<u32, usize> = sheet_verts.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let n_sverts = sheet_verts.len();
    let n_total = n_edges + n_sverts;

    // 1. Assemble 3D volume curl-curl + gauge regularization matrix over Nedelec edges 0..n_edges
    //    (transforming curl-curl inside IE domains via (1 / det J_IE) (J_IE curl w_a)^T mu_r^-1 (J_IE curl w_b))
    let nc = mixed.n_cells();
    let k_vol = assemble_parallel(nc, n_total, n_total, |range, coo| {
        let mut scratch = NedEvalScratch::default();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut w_out = Vec::new();
        let mut curl_out = Vec::new();
        let mut tilde_c: Vec<[f64; 3]> = Vec::new();
        let mut ke = Vec::new();

        for c in range {
            let ct = mixed.cell_type(c);
            let iso = IsoparametricCell::from_mesh(&mixed, c);
            let ie_opt = ie_map.get(&mixed.cell_labels[c]);
            let deg = if ie_opt.is_some() || iso.curved { 4 } else { 2 };
            let rule = quad_rule_mixed(ct, deg);
            let dofs = ned.cell_dofs(c);
            let nloc = dofs.len();
            ke.clear();
            ke.resize(nloc * nloc, 0.0);
            let inv_mur = 1.0 / mu_r[c];

            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jac_k = eval_cell_basis_3d(&mixed, &ned, c, &iso, xi, &mut scratch, &mut w_out, &mut curl_out);
                let w_phys = wq * jac_k.det.abs();

                if let Some(ie) = ie_opt {
                    let _ = iso.eval_into(xi, &mut sn, &mut sdn);
                    let jac_ie = ie.eval_cell_jacobian(&iso, &jac_k, &sn, &sdn, 3);
                    let det_ie = jac_ie.det.re.abs().max(1e-30);
                    let w_curl = w_phys * inv_mur / det_ie;
                    tilde_c.clear();
                    for cu in &curl_out {
                        let mut tc = [0.0; 3];
                        for i in 0..3 {
                            for m in 0..3 {
                                tc[i] += jac_ie.j[i][m].re * cu[m];
                            }
                        }
                        tilde_c.push(tc);
                    }
                    for i in 0..nloc {
                        for j in 0..nloc {
                            let curl_dot = dot3(tilde_c[i], tilde_c[j]);
                            let mass_dot = dot3(w_out[i], w_out[j]);
                            ke[i * nloc + j] += w_curl * curl_dot + w_phys * eps_reg * mass_dot;
                        }
                    }
                } else {
                    for i in 0..nloc {
                        for j in 0..nloc {
                            let curl_dot = dot3(curl_out[i], curl_out[j]);
                            let mass_dot = dot3(w_out[i], w_out[j]);
                            ke[i * nloc + j] += w_phys * (inv_mur * curl_dot + eps_reg * mass_dot);
                        }
                    }
                }
            }
            coo.push_block(dofs, dofs, &ke);
        }
    });

    // 2. Assemble 2D surface kinetic coupling term on deduplicated london_sheet triangles
    let qtri = quad_rule(2, 2);
    let mut coo_surf = Coo::new(n_total, n_total);
    let pt3 = |v: u32| -> [f64; 3] { mixed.point3(v as usize) };

    for (&[v0, v1, v2], &lk) in &dedup_sheets {
        let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
        let (area, _normal, g_surf) = triangle_surface_gradients(p0, p1, p2);
        if area <= 1e-30 {
            continue;
        }
        let e01 = mixed.edge_index(v0, v1).unwrap();
        let e02 = mixed.edge_index(v0, v2).unwrap();
        let e12 = mixed.edge_index(v1, v2).unwrap();
        let s0 = n_edges + vert_to_sidx[&v0];
        let s1 = n_edges + vert_to_sidx[&v1];
        let s2 = n_edges + vert_to_sidx[&v2];
        let dofs = [e01, e02, e12, s0, s1, s2];
        let coeff = (MU0 / lk) * area;

        let mut ke = [0.0; 36];
        for q in 0..qtri.w.len() {
            let l = &qtri.lam[q];
            let w01 = sub3(scale3(l[0], g_surf[1]), scale3(l[1], g_surf[0]));
            let w02 = sub3(scale3(l[0], g_surf[2]), scale3(l[2], g_surf[0]));
            let w12 = sub3(scale3(l[1], g_surf[2]), scale3(l[2], g_surf[1]));
            let basis = [w01, w02, w12, scale3(-1.0, g_surf[0]), scale3(-1.0, g_surf[1]), scale3(-1.0, g_surf[2])];
            let wq = coeff * qtri.w[q];
            for i in 0..6 {
                for j in 0..6 {
                    ke[i * 6 + j] += wq * dot3(basis[i], basis[j]);
                }
            }
        }
        coo_surf.push_block(&dofs, &dofs, &ke);
    }

    let k_surf = coo_surf.to_csr();

    // 3. Assemble 1st-order asymptotic curl-Robin boundary term on `far_field` and `Exponential` IE exterior facets
    let ref_labels: Vec<i32> = excitations.iter().flat_map(|e| e.labels.iter().copied()).collect();
    let f2c = build_facet_to_cell_mixed(&mixed);
    let mut robin_facets: Vec<(usize, usize, f64, [f64; 3], f64)> = Vec::new();
    for bc in &spec.bcs {
        if let BcKind::FarField { decay_power, center, .. } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !ground_labels.contains(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid_mixed(&mixed, &ref_labels));
                let gamma = if decay_power.is_finite() && *decay_power > 0.0 { *decay_power } else { 1.0 };
                for f in 0..mixed.n_facets() {
                    if ff_labs.contains(&mixed.facet_labels[f]) {
                        if let Some(c_f) = f2c[f] {
                            robin_facets.push((f, c_f, gamma, c, 1.0));
                        }
                    }
                }
            }
        }
    }
    let (ie_rat_facets, ie_exp_facets) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
    for &(f_idx, c_f, s_omega) in &ie_exp_facets {
        let dom = mixed.cell_labels[c_f];
        let c = ie_map.get(&dom).map_or_else(|| farfield_centroid_mixed(&mixed, &ref_labels), |ie| ie.center);
        robin_facets.push((f_idx, c_f, 1.0, c, s_omega));
    }
    let k_vr = if !robin_facets.is_empty() {
        let k_rob = assemble_curl_robin_3d_mixed(&mixed, &ned, n_total, &mu_r, &robin_facets);
        k_vol.add_scaled(&k_rob, 1.0)
    } else {
        k_vol
    };
    let k_mag = k_vr.add_scaled(&k_surf, 1.0);

    // 4. Pin one reference phi_s vertex per connected sheet component (plus terminal end-pins)
    let mut fixed_set = BTreeSet::new();
    if n_sverts > 0 {
        let mut parent: Vec<usize> = (0..n_sverts).collect();
        fn find(p: &mut [usize], mut i: usize) -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        }
        for &[v0, v1, v2] in dedup_sheets.keys() {
            let (i0, i1, i2) = (vert_to_sidx[&v0], vert_to_sidx[&v1], vert_to_sidx[&v2]);
            for (a, b) in [(i0, i1), (i1, i2), (i0, i2)] {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        let mut seen_roots = BTreeSet::new();
        for i in 0..n_sverts {
            let r = find(&mut parent, i);
            if seen_roots.insert(r) {
                fixed_set.insert(n_edges + r);
            }
        }
    }

    // Exterior PEC / Ground Nedelec edge constraints + outer Rational IE facets (`\tilde{R} = \infty`)
    for f in 0..mixed.n_facets() {
        let lab = mixed.facet_labels[f];
        if ground_labels.contains(&lab) {
            for e in ned.facet_dofs(&mixed, f) {
                if !sheet_edges_set.contains(&e) {
                    fixed_set.insert(e);
                }
            }
        }
    }
    for &f_idx in &ie_rat_facets {
        for e in ned.facet_dofs(&mixed, f_idx) {
            if !sheet_edges_set.contains(&e) {
                fixed_set.insert(e);
            }
        }
    }

    // 5. Build unit-fluxoid cohomology right-hand sides & solve non-cancelling residual system
    let (metas, rhs_u_list, c_flux_list) = build_mixed_3d_fluxoid_rhs(
        &excitations,
        &dedup_sheets,
        &vert_to_sidx,
        |labs| {
            (0..mixed.n_facets())
                .filter(|&f| labs.contains(&mixed.facet_labels[f]))
                .filter_map(|f| {
                    let fv = mixed.facet(f);
                    if fv.len() == 3 {
                        Some([fv[0], fv[1], fv[2]])
                    } else {
                        None
                    }
                })
                .collect()
        },
        pt3,
        |a, b| mixed.edge_index(a, b),
        n_edges,
        n_total,
        l_dom,
        &mut fixed_set,
    )?;

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let (l_mat, u_all, log) = solve_mixed_3d_fluxoids_and_inductance(
        &excitations,
        &metas,
        &dedup_sheets,
        &vert_to_sidx,
        &k_vr,
        &k_mag,
        &rhs_u_list,
        &c_flux_list,
        &fixed,
        pt3,
        |a, b| mixed.edge_index(a, b),
        n_edges,
        n_sverts,
    )?;

    let n_exc = excitations.len();
    let mut solutions = Vec::with_capacity(n_exc);
    let mut scalars = Vec::new();
    for p in 0..n_exc {
        let a_edges = &u_all[p][..n_edges];
        let (fields, cell_fields) = mixed_3d_fields_mixed(&mixed, &ned, &simplex, &proj, order, a_edges, &ie_map);
        if p == 0 {
            let bn = norm_of(&components(&fields[1].values, 3));
            let w_m = 0.5 * l_mat[0][0] * excitations[0].amplitude.powi(2);
            scalars.push(("stored magnetic energy W_m (J)".into(), Scalar::Num(w_m)));
            scalars.push(("max |B| (T)".into(), Scalar::Num(minmax(&bn).1)));
        }
        solutions.push(Solution {
            label: format!("terminal {}", excitations[p].name),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(l_mat[p][p]),
            ..Default::default()
        });
    }

    push_inductance_scalars(&mut scalars, &excitations, &l_mat, "mixed_3d");

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: 3,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n_total,
    })
}

/// Solves a magnetostatic problem natively on a mixed mesh (`Quad`, `Prism`, `Hex`) and/or with
/// domain-applied static Infinite Elements (`BcKind::InfiniteElement`), preserving the pure-simplex
/// fast path when `mesh.is_pure_simplex() && !spec.has_infinite_element()`.
pub fn solve_magnetostatics_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String> {
    if mesh_has_curved_edges(mesh) {
        return Err("unsupported: magnetostatics does not support curved edges (R32)".into());
    }
    if mesh.is_pure_simplex()
        && !spec.has_infinite_element()
        && !spec.bcs.iter().any(|b| matches!(b.kind, BcKind::FarField { .. }))
    {
        let d = discretise(mesh).map_err(|e| e.to_string())?;
        return solve_magnetostatics(&d, spec);
    }
    let Study::Magnetostatic { formulation } = &spec.study else {
        return Err("solve_magnetostatics requires a magnetostatic study".into());
    };
    let lower = formulation.trim().to_ascii_lowercase();
    let form = match lower.as_str() {
        "" | "auto" => {
            if mesh.dim() == 2 {
                "cross_section"
            } else {
                "mixed_3d"
            }
        }
        other => other,
    };
    match form {
        "cross_section" => solve_cross_section_mixed(mesh, spec),
        "planar_stream" => {
            let d = discretise(mesh).map_err(|e| e.to_string())?;
            solve_planar_stream(&d, spec)
        }
        "mixed_3d" => solve_mixed_3d_mixed(mesh, spec),
        other => Err(format!(
            "unknown magnetostatic formulation '{other}' (expected 'auto', 'cross_section', 'planar_stream', or 'mixed_3d')"
        )),
    }
}
