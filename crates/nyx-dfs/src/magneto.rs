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

use crate::common::{cell_gradient, components, minmax, norm_of};
use crate::problem::{BcKind, Problem, Study};
use femcore::fem::*;
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::Coo;
use femcore::linalg::driver::{dense_eigh, DirectSolver};
use femcore::props::LinearSolver;
use nyx_physics::{cell_fields_of, Discretisation, Field, InductanceMatrixData, Scalar, Solution, SolveResult};
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
            BcKind::LondonSheet { .. } | BcKind::Ground => {}
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

fn solve_cross_section(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    if mesh.dim != 2 {
        return Err("the 'cross_section' magnetostatic formulation requires a 2D mesh".into());
    }
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
    let excitations = collect_excitations(spec)?;
    let mu_r = spec.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    if mu_r.iter().any(|&m| m <= 0.0) {
        return Err("material property 'mu_r' must be positive".into());
    }
    let coeff: Vec<f64> = mu_r.iter().map(|&m| 1.0 / (MU0 * m)).collect();
    let mut k = assemble_laplace(mesh, order, &coeff);

    let mut sheet_lk: BTreeMap<i32, f64> = BTreeMap::new();
    let mut ground_labels: Vec<i32> = Vec::new();
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
            _ => {}
        }
    }

    let term_label_set: BTreeSet<i32> = excitations.iter().flat_map(|e| e.labels.iter().copied()).collect();

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

    if ground_labels.is_empty() {
        let exterior = mesh.exterior_facets();
        for f in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f];
            if exterior[f]
                && !term_label_set.contains(&lab)
                && !sheet_lk.contains_key(&lab)
                && !ground_labels.contains(&lab)
            {
                ground_labels.push(lab);
            }
        }
    }
    let fixed = boundary_dofs(mesh, &ground_labels, order);
    if fixed.is_empty() {
        return Err("2D cross_section magnetostatics needs an exterior or ground reference boundary (A_z = 0)".into());
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
    let t0 = Instant::now();
    if mesh.dim != 3 {
        return Err("the 'mixed_3d' magnetostatic formulation requires a 3D mesh".into());
    }
    let order = spec.order;
    check_boundary_labels(mesh, spec.bcs.iter().flat_map(|b| b.labels.iter().copied()))?;
    let excitations = collect_excitations(spec)?;
    let mu_r = spec.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    if mu_r.iter().any(|&m| m <= 0.0) {
        return Err("material property 'mu_r' must be positive".into());
    }

    let (lo, hi) = mesh.bbox();
    let l_dom = ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt().max(1e-12);
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

    // 4. Build right-hand sides mu_0 * I_p * b_p for all current_terminal and flux_loop excitations
    let n_exc = excitations.len();
    let mut rhs_full_list = vec![vec![0.0; n_total]; n_exc];

    for (p, ex) in excitations.iter().enumerate() {
        // Deduplicate facets belonging to excitation `ex`
        let mut ex_tris: BTreeSet<[u32; 3]> = BTreeSet::new();
        for f in facets_with(mesh, &ex.labels) {
            let fv = mesh.facet(f);
            let mut key = [fv[0], fv[1], fv[2]];
            key.sort_unstable();
            ex_tris.insert(key);
        }
        if ex_tris.is_empty() {
            return Err(format!("excitation '{}' has no boundary facets", ex.name));
        }

        // Compute centroid and span of excitation facets
        let mut cx = [0.0; 3];
        let mut total_area = 0.0;
        for &[v0, v1, v2] in &ex_tris {
            let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
            let (area, _, _) = triangle_surface_gradients(p0, p1, p2);
            let cen = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
            for k in 0..3 {
                cx[k] += area * cen[k];
            }
            total_area += area;
        }
        if total_area <= 1e-30 {
            return Err(format!("excitation '{}' has zero facet area", ex.name));
        }
        for k in 0..3 {
            cx[k] /= total_area;
        }

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
        let density_scale = if span > 1e-12 { span / total_area } else { 1.0 / total_area };

        // Pin phi_s at the two end-caps of a linear current_terminal sheet so grad_|| phi_s
        // cannot open-circuit the longitudinal current injected along d_unit
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

        let rhs = &mut rhs_full_list[p];
        for &[v0, v1, v2] in &ex_tris {
            let (p0, p1, p2) = (pt3(v0), pt3(v1), pt3(v2));
            let (area, normal, g_surf) = triangle_surface_gradients(p0, p1, p2);
            if area <= 1e-30 {
                continue;
            }
            let e01 = mesh.edge_index(v0, v1);
            let e02 = mesh.edge_index(v0, v2);
            let e12 = mesh.edge_index(v1, v2);

            for q in 0..qtri.w.len() {
                let l = &qtri.lam[q];
                let rq = [
                    l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0],
                    l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1],
                    l[0] * p0[2] + l[1] * p1[2] + l[2] * p2[2],
                ];
                // Direction of current at rq: azimuthal around loop center cx for FluxLoop, or projected d_unit for CurrentTerminal
                let dir_vec = if ex.is_flux_loop {
                    let rel = sub3(rq, cx);
                    let tan = cross3(normal, rel);
                    let r_len = dot3(tan, tan).sqrt();
                    if r_len > 1e-30 {
                        // Harmonic cohomology 1-form scaled by density_scale
                        scale3(density_scale / r_len, tan)
                    } else {
                        [0.0; 3]
                    }
                } else {
                    let d_proj = sub3(d_unit, scale3(dot3(d_unit, normal), normal));
                    let dp_len = dot3(d_proj, d_proj).sqrt();
                    if dp_len > 1e-12 {
                        scale3(density_scale / dp_len, d_proj)
                    } else {
                        scale3(density_scale, d_unit)
                    }
                };

                let w01 = sub3(scale3(l[0], g_surf[1]), scale3(l[1], g_surf[0]));
                let w02 = sub3(scale3(l[0], g_surf[2]), scale3(l[2], g_surf[0]));
                let w12 = sub3(scale3(l[1], g_surf[2]), scale3(l[2], g_surf[1]));
                let wq = MU0 * ex.amplitude * area * qtri.w[q];

                rhs[e01] += wq * dot3(dir_vec, w01);
                rhs[e02] += wq * dot3(dir_vec, w02);
                rhs[e12] += wq * dot3(dir_vec, w12);

                if ex.is_flux_loop {
                    if let (Some(&i0), Some(&i1), Some(&i2)) =
                        (vert_to_sidx.get(&v0), vert_to_sidx.get(&v1), vert_to_sidx.get(&v2))
                    {
                        rhs[n_edges + i0] -= wq * dot3(dir_vec, g_surf[0]);
                        rhs[n_edges + i1] -= wq * dot3(dir_vec, g_surf[1]);
                        rhs[n_edges + i2] -= wq * dot3(dir_vec, g_surf[2]);
                    }
                }
            }
        }
    }

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let free = free_dofs(n_total, &fixed);
    let n_free = free.len();

    let kff = k_mag.submatrix(&free, &free);
    let rhs_free_list: Vec<Vec<f64>> = rhs_full_list.iter().map(|r| free.iter().map(|&i| r[i]).collect()).collect();

    let solver = DirectSolver::new(&kff, true)?;
    let sols_free = solver.solve_many(&rhs_free_list);
    let log = vec![format!(
        "3D mixed_3d London magnetostatics: 1 factorization ({n_free} free dofs: {n_edges} Nedelec edges + {n_sverts} sheet phi_s nodes on {} deduplicated triangles) served {n_exc} excitation(s)",
        dedup_sheets.len()
    )];

    let mut u_all = vec![vec![0.0; n_total]; n_exc];
    for p in 0..n_exc {
        for (a, &i) in free.iter().enumerate() {
            u_all[p][i] = sols_free[p][a];
        }
    }

    let ku_all: Vec<Vec<f64>> = u_all.iter().map(|u| k_mag.matvec(u)).collect();
    let mut l_mat = vec![vec![0.0; n_exc]; n_exc];
    for i in 0..n_exc {
        for j in i..n_exc {
            let e_ij: f64 = u_all[i].iter().zip(&ku_all[j]).map(|(a, b)| a * b).sum();
            let e_ji: f64 = if i == j { e_ij } else { u_all[j].iter().zip(&ku_all[i]).map(|(a, b)| a * b).sum() };
            let lij = 0.5 * (e_ij + e_ji) / (MU0 * excitations[i].amplitude * excitations[j].amplitude);
            l_mat[i][j] = lij;
            l_mat[j][i] = lij;
        }
    }

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

/// Solves a magnetostatic problem using the requested or auto-selected formulation.
pub fn solve_magnetostatics(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
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
