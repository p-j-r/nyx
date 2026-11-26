// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Electrostatics on P1/P2 simplices: `-div(eps0 eps_r grad V) = rho_v`, `V = V0` on fixed-potential boundaries, `n.D = -rho_s` on
//! charged ones. Three stages (ADR 9): `assemble_electrostatics`, `solve_potential`, `electrostatics_fields`; `solve_electrostatics` runs them.

use crate::common::{
    assemble_boundary_load_field, assemble_boundary_load_mixed, assemble_boundary_mass_curved_aware,
    assemble_farfield_robin_mixed, assemble_ie_outer_robin_mixed, assemble_laplace_ie_mixed, assemble_source_field,
    assemble_source_mixed, cell_gradient, cell_tensor_mixed, cell_values_mixed, check_boundary_labels_mixed,
    classify_ie_outer_facets_mixed, components, dirichlet_constraints_field, farfield_centroid_mixed,
    farfield_positive_measure_mixed, minmax, norm_of, plane_tensors, plane_tensors_mixed, resolve_ie_domains_mixed,
    sample_mixed_gradients_on_simplex, solve_spd_mixed,
};
use crate::problem::{BcKind, Problem};
use femcore::fem::*;
use femcore::fem_mixed::{LagrangeSpace, MixedLagrangeDofMap};
use femcore::fields::{node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::driver::{dense_eigh, DirectSolver};
use femcore::material::{self, M3};
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::props::{LinearSolver, TensorKind};
use femcore::vvuq::{SourceField, UniformSource};
use nyx_physics::{
    cell_fields_of, discretise, discretise_mixed, CapacitanceMatrixData, Discretisation, ElectroAmrData, Field, Mesh,
    MixedToSimplexProjector, Scalar, Solution, SolveResult, TlsMeta,
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
    /// whether at least one spherical 1/r asymptotic Robin open boundary condition is active
    pub has_open: bool,
}

fn validate_curved_tets_q8(mesh: &Discretisation, cquad: &Quad) -> Result<(), String> {
    for c in 0..mesh.n_cells() {
        if mesh.cell_is_curved(c) {
            let cv = mesh.cell(c);
            let mut v4 = [[0.0; 3]; 4];
            for i in 0..4 {
                let p = mesh.point(cv[i] as usize);
                v4[i] = [p[0], p[1], p[2]];
            }
            let det_straight = CurvedTet::straight(v4).jacobian(&[0.25, 0.25, 0.25, 0.25]).1;
            if !(det_straight.is_finite() && det_straight.abs() > 1e-300) {
                return Err(format!(
                    "curved tetrahedron {c} has degenerate straight base determinant ({det_straight:.3e})"
                ));
            }
            let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
            for q in 0..cquad.w.len() {
                let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                let (_, det) = ct.jacobian(&lam);
                if !(det.is_finite() && det * det_straight > 0.0) {
                    return Err(format!(
                        "curved tetrahedron {c} has non-positive or inverted Jacobian determinant at Gauss point {q} (det = {det:.3e}, det_straight = {det_straight:.3e})"
                    ));
                }
            }
        }
    }
    Ok(())
}

pub(crate) fn assemble_laplace_curved_3d(mesh: &Discretisation, order: usize, coeff: &[f64]) -> Result<Csr, String> {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(8);
    validate_curved_tets_q8(mesh, &cquad)?;
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    Ok(assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
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
    }))
}

pub(crate) fn assemble_laplace_tensor_curved_3d(mesh: &Discretisation, order: usize, kt: &[M3]) -> Result<Csr, String> {
    let re = RefElement::new(3, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let cquad = quad_rule_collapsed(8);
    validate_curved_tets_q8(mesh, &cquad)?;
    let cshapes: Vec<Vec<f64>> = cquad.lam.iter().map(|lam| shape(lam, order).1).collect();
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    Ok(assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
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
    }))
}

/// Enriches a 3D mesh for P2 Barsoum quarter-point singularity capture (`singular_edges = "quarter_point"`)
/// when `mesh` has not already been enriched by `sheet3d`:
/// finds all open perimeter vertices of interior zero-thickness conductor sheet facets (`Terminal`, `Ground`,
/// or `Potential`), shifts mid-edge geometry nodes via [`femcore::meshing::singularity::apply_quarter_point_nodes`],
/// clamps `xi_QP` upward on curved rims, and validates positive `CurvedTet` Jacobians across all 216 Gauss points
/// of `quad_rule_collapsed(8)`.
pub fn enrich_quarter_point_mesh(mesh: &Discretisation, spec: &Problem) -> Option<(Discretisation, usize, usize)> {
    if spec.order != 2 || mesh.dim != 3 || spec.singular_edges != "quarter_point" || mesh.has_qp() {
        return None;
    }
    let mut cond_labels = Vec::new();
    for bc in &spec.bcs {
        if matches!(bc.kind, BcKind::Potential(_) | BcKind::Ground | BcKind::Terminal(_)) {
            for &l in &bc.labels {
                if spec.sheet_thickness_m.get(&l).copied().unwrap_or(0.0) <= 0.0 && !cond_labels.contains(&l) {
                    cond_labels.push(l);
                }
            }
        }
    }
    if cond_labels.is_empty() {
        return None;
    }
    let rim_verts = femcore::meshing::singularity::sheet_rim_vertices(mesh, &cond_labels);
    if !rim_verts.iter().any(|&r| r) {
        return None;
    }
    let mut curved_rim_edges = BTreeSet::new();
    for ([a, b], _) in mesh.edge_node_list() {
        if rim_verts.get(a as usize).copied().unwrap_or(false) && rim_verts.get(b as usize).copied().unwrap_or(false) {
            curved_rim_edges.insert([a.min(b), a.max(b)]);
        }
    }
    let mut enriched = mesh.clone();
    let rep = femcore::meshing::singularity::apply_quarter_point_nodes(
        &mut enriched,
        &rim_verts,
        femcore::meshing::singularity::DEFAULT_QUARTER_POINT_XI,
        &curved_rim_edges,
    );
    if rep.shifted_edges == 0 {
        return None;
    }
    Some((enriched, rep.rim_vertices, rep.shifted_edges))
}

/// Prepares the working discretisation for electrostatics:
/// - Strips quarter-point edges if `spec.order < 2` (retaining smooth `AF-5` `curved` edges);
/// - Deduplicates quarter-point enrichment if `sheet3d` (`set_mesh`) already enriched `mesh`;
/// - Otherwise enriches `mesh` when `spec.singular_edges == "quarter_point"` and `spec.order == 2`.
pub fn prepare_electro_mesh(mesh: &Discretisation, spec: &Problem) -> (Option<(Discretisation, usize, usize)>, bool) {
    if mesh.dim == 3 && mesh.has_qp() {
        if spec.order < 2 {
            let qp_set: BTreeSet<[u32; 2]> = mesh.qp_edge_list().into_iter().collect();
            let smooth_nodes: Vec<([u32; 2], [f64; 3])> =
                mesh.edge_node_list().into_iter().filter(|(k, _)| !qp_set.contains(k)).collect();
            let mut stripped = mesh.clone();
            stripped.set_edge_nodes(&smooth_nodes);
            stripped.set_qp_edges(&[]);
            return (Some((stripped, 0, 0)), false);
        }
        return (None, true);
    }
    (enrich_quarter_point_mesh(mesh, spec), false)
}

/// Assembles the 3D spherical `1/r` asymptotic Robin surface stiffness matrix on boundary triangles
/// matching `open_origins` (`label -> origin r0`):
/// `K_open,ab = int_{Gamma_open} eps0 (n^T eps_r n) * max((x - r0) . n, 0) / |x - r0|^2 * phi_a(x) phi_b(x) dS`.
fn assemble_open_robin(
    mesh: &Discretisation,
    order: usize,
    open_origins: &BTreeMap<i32, [f64; 3]>,
    tensor_perm: Option<&[M3]>,
    scalar_perm: &[f64],
) -> Result<Csr, String> {
    if mesh.dim != 3 {
        return Err("open boundary condition (spherical 1/r asymptotic Robin) requires a 3D mesh".into());
    }
    let nv = mesh.n_points();
    let mut on_open_v = vec![false; nv];
    for f in 0..mesh.n_facets() {
        if open_origins.contains_key(&mesh.facet_labels[f]) {
            for &v in mesh.facet(f) {
                on_open_v[v as usize] = true;
            }
        }
    }
    let mut face_to_cell: std::collections::HashMap<[u32; 3], (usize, u32)> = std::collections::HashMap::new();
    for c in 0..mesh.n_cells() {
        let cv = mesh.cell(c);
        let cnt = cv.iter().filter(|&&v| on_open_v[v as usize]).count();
        if cnt < 3 {
            continue;
        }
        for skip in 0..4 {
            let [i0, i1, i2] = match skip {
                0 => [1, 2, 3],
                1 => [0, 2, 3],
                2 => [0, 1, 3],
                _ => [0, 1, 2],
            };
            let (u0, u1, u2) = (cv[i0], cv[i1], cv[i2]);
            if on_open_v[u0 as usize] && on_open_v[u1 as usize] && on_open_v[u2 as usize] {
                let mut key = [u0, u1, u2];
                key.sort_unstable();
                face_to_cell.entry(key).or_insert((c, cv[skip]));
            }
        }
    }

    let rule = facet_rule(2);
    let nloc = if order == 1 { 3 } else { 6 };
    let shapes: Vec<Vec<f64>> = rule.lam.iter().map(|lam| shape(lam, order).0).collect();
    let n = mesh.n_dofs(order);
    let mut coo = femcore::linalg::csr::Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut ke = vec![0.0; nloc * nloc];
    let mut positive_facets = 0usize;

    for f in 0..mesh.n_facets() {
        let Some(&r0) = open_origins.get(&mesh.facet_labels[f]) else {
            continue;
        };
        let fv = mesh.facet(f);
        let mut key = [fv[0], fv[1], fv[2]];
        key.sort_unstable();
        let &(c, opp) =
            face_to_cell.get(&key).ok_or_else(|| format!("open boundary facet {f} is not attached to any 3D cell"))?;
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        let p2 = mesh.point(fv[2] as usize);
        let p_opp = mesh.point(opp as usize);
        let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let cx = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        let norm_cx = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
        if norm_cx <= 1e-30 {
            continue;
        }
        let area = 0.5 * norm_cx;
        let cx_hat = [cx[0] / norm_cx, cx[1] / norm_cx, cx[2] / norm_cx];
        let mut n_hat = cx_hat;
        let out_dot = (p0[0] - p_opp[0]) * n_hat[0] + (p0[1] - p_opp[1]) * n_hat[1] + (p0[2] - p_opp[2]) * n_hat[2];
        if out_dot < 0.0 {
            n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
        }

        // Reject origin lying directly on the triangle facet
        let h_plane = (r0[0] - p0[0]) * n_hat[0] + (r0[1] - p0[1]) * n_hat[1] + (r0[2] - p0[2]) * n_hat[2];
        let h_tol = 1e-12 * area.sqrt().max(1e-15);
        if h_plane.abs() <= h_tol {
            let w = [
                r0[0] - h_plane * n_hat[0] - p0[0],
                r0[1] - h_plane * n_hat[1] - p0[1],
                r0[2] - h_plane * n_hat[2] - p0[2],
            ];
            let wxb = [w[1] * b[2] - w[2] * b[1], w[2] * b[0] - w[0] * b[2], w[0] * b[1] - w[1] * b[0]];
            let axw = [a[1] * w[2] - a[2] * w[1], a[2] * w[0] - a[0] * w[2], a[0] * w[1] - a[1] * w[0]];
            let lam1 = (wxb[0] * cx_hat[0] + wxb[1] * cx_hat[1] + wxb[2] * cx_hat[2]) / norm_cx;
            let lam2 = (axw[0] * cx_hat[0] + axw[1] * cx_hat[1] + axw[2] * cx_hat[2]) / norm_cx;
            let lam0 = 1.0 - lam1 - lam2;
            if lam0 >= -1e-9 && lam1 >= -1e-9 && lam2 >= -1e-9 {
                return Err("open boundary origin coincides with an open boundary facet".into());
            }
        }

        let eps_nn = match tensor_perm {
            Some(t) => {
                let mut s = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        s += n_hat[i] * t[c][i][j] * n_hat[j];
                    }
                }
                s
            }
            None => scalar_perm[c],
        };

        ke.fill(0.0);
        let r2_min = 1e-24 * area.max(1e-30);
        let mut facet_pos = false;
        for (q, lam) in rule.lam.iter().enumerate() {
            let xq = [
                lam[0] * p0[0] + lam[1] * p1[0] + lam[2] * p2[0],
                lam[0] * p0[1] + lam[1] * p1[1] + lam[2] * p2[1],
                lam[0] * p0[2] + lam[1] * p1[2] + lam[2] * p2[2],
            ];
            let dx = [xq[0] - r0[0], xq[1] - r0[1], xq[2] - r0[2]];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 <= r2_min {
                return Err("open boundary origin coincides with an open boundary facet".into());
            }
            let r_dot_n = (dx[0] * n_hat[0] + dx[1] * n_hat[1] + dx[2] * n_hat[2]).max(0.0);
            let w_alpha = area * rule.w[q] * eps_nn * (r_dot_n / r2);
            if w_alpha > 0.0 {
                facet_pos = true;
                let nv = &shapes[q];
                for ia in 0..nloc {
                    for ib in 0..nloc {
                        ke[ia * nloc + ib] += w_alpha * nv[ia] * nv[ib];
                    }
                }
            }
        }
        if facet_pos {
            positive_facets += 1;
            mesh.facet_dofs(f, order, &mut dofs);
            coo.push_block(&dofs, &dofs, &ke);
        }
    }
    if !open_origins.is_empty() && positive_facets == 0 {
        return Err("open boundary origin produced zero outward Robin weight on all open boundary facets".into());
    }
    Ok(coo.to_csr())
}

fn assemble_electrostatics_on_work_mesh_with_source<const D: usize, S: SourceField<D>>(
    m_work: &Discretisation,
    spec: &Problem,
    source: &S,
    robin_bcs: &[(Vec<i32>, f64)],
) -> Result<ElectroAssembly, String> {
    if !source.is_uniform() && m_work.dim != D {
        return Err(format!("SourceField<{D}> called on {}D mesh", m_work.dim));
    }
    let order = spec.order;
    check_boundary_labels(
        m_work,
        spec.bcs
            .iter()
            .filter(|b| !matches!(b.kind, BcKind::InfiniteElement(_)))
            .flat_map(|b| b.labels.iter().copied())
            .chain(robin_bcs.iter().flat_map(|(labs, _)| labs.iter().copied())),
    )?;
    let eps = spec.materials.cell_values(m_work, "eps_r", Some(1.0))?;
    let rho = spec.materials.cell_values(m_work, "rho_v", Some(0.0))?;
    let coeff: Vec<f64> = eps.iter().map(|e| EPS0 * e).collect();
    let eps_tensors_rel: Option<Vec<M3>> = if spec.materials.has_tensor(TensorKind::Permittivity) {
        Some(plane_tensors(m_work, spec.materials.cell_tensor(m_work, TensorKind::Permittivity, Some(1.0))?))
    } else {
        None
    };
    let tensor_perm: Option<Vec<M3>> = eps_tensors_rel.as_ref().map(|tr| {
        let mut t = tr.clone();
        t.iter_mut().for_each(|m| m.iter_mut().for_each(|r| r.iter_mut().for_each(|x| *x *= EPS0)));
        t
    });
    let mut ie_rat_facets = Vec::new();
    let mut ie_exp_facets = Vec::new();
    let mut ie_map_opt = None;
    let mut mixed_opt = None;
    let mut dof_map_opt = None;
    let mut k = if spec.has_infinite_element() {
        let mixed = MixedDiscretisation::from_simplex_mesh(m_work)?;
        let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
        let dof_map = MixedLagrangeDofMap::new(&mixed, space);
        let ie_map = resolve_ie_domains_mixed(&mixed, spec)?;
        let (rf, ef) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
        let km = assemble_laplace_ie_mixed(&mixed, &dof_map, &coeff, tensor_perm.as_deref(), &ie_map);
        ie_rat_facets = rf;
        ie_exp_facets = ef;
        ie_map_opt = Some(ie_map);
        mixed_opt = Some(mixed);
        dof_map_opt = Some(dof_map);
        km
    } else if let Some(ref t) = tensor_perm {
        if m_work.dim == 3 && m_work.has_curved() {
            assemble_laplace_tensor_curved_3d(m_work, order, t)?
        } else {
            assemble_laplace_tensor(m_work, order, t)
        }
    } else if m_work.dim == 3 && m_work.has_curved() {
        assemble_laplace_curved_3d(m_work, order, &coeff)?
    } else {
        assemble_laplace(m_work, order, &coeff)
    };
    let n = k.nrows();
    let mut f = if source.is_uniform() {
        if rho.iter().any(|&r| r != 0.0) {
            assemble_source(m_work, order, &rho)
        } else {
            vec![0.0; n]
        }
    } else {
        let extra = rho.iter().any(|&r| r != 0.0).then_some(rho.as_slice());
        assemble_source_field::<D, S>(m_work, order, source, 0.0, extra)?
    };
    let mut volt: BTreeMap<i32, f64> = BTreeMap::new();
    let mut dirichlet_labels: BTreeSet<i32> = BTreeSet::new();
    let (lo, hi) = m_work.bbox();
    let default_origin = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
    let mut open_origins: BTreeMap<i32, [f64; 3]> = BTreeMap::new();
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
                BcKind::Potential(v) => {
                    volt.insert(lab, *v);
                    dirichlet_labels.insert(lab);
                }
                BcKind::Ground => {
                    volt.insert(lab, 0.0);
                    dirichlet_labels.insert(lab);
                }
                BcKind::Terminal(name) => {
                    if name.trim().is_empty() {
                        return Err("terminal condition requires a non-empty 'name'".into());
                    }
                    volt.entry(lab).or_insert(0.0);
                    dirichlet_labels.insert(lab);
                }
                BcKind::Open { origin } => {
                    if m_work.dim != 3 {
                        return Err(
                            "open boundary condition (spherical 1/r asymptotic Robin) requires a 3D mesh".into()
                        );
                    }
                    let r0 = origin.unwrap_or(default_origin);
                    if r0.iter().any(|x| !x.is_finite()) {
                        return Err("open boundary 'origin' must be finite".into());
                    }
                    open_origins.insert(lab, r0);
                }
                BcKind::SurfaceCharge(s) => {
                    if source.is_uniform() {
                        if *s != 0.0 {
                            let l = assemble_boundary_load(m_work, &[lab], *s, order);
                            f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                        }
                    } else {
                        let l = assemble_boundary_load_field::<D, S>(m_work, &[lab], order, source, 0.0)?;
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                BcKind::FarField { .. } | BcKind::InfiniteElement(_) => {}
                other => return Err(format!("'{}' is not an electrostatics boundary condition", other.name())),
            }
        }
    }
    let mut has_open = !open_origins.is_empty() || !ie_rat_facets.is_empty() || !ie_exp_facets.is_empty();
    if !open_origins.is_empty() {
        let k_open = assemble_open_robin(m_work, order, &open_origins, tensor_perm.as_deref(), &coeff)?;
        k = k.add_scaled(&k_open, 1.0);
    }
    for (labs, alpha) in robin_bcs {
        if *alpha != 0.0 {
            has_open = true;
            let k_rob = assemble_boundary_mass_curved_aware(m_work, labs, order, *alpha);
            k = k.add_scaled(&k_rob, 1.0);
        }
        if !source.is_uniform() {
            let l = assemble_boundary_load_field::<D, S>(m_work, labs, order, source, 0.0)?;
            f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
        }
    }
    if volt.is_empty() && !has_open {
        return Err("electrostatics needs at least one potential, ground or terminal boundary".into());
    }
    let ref_labels: Vec<i32> = volt.keys().copied().collect();
    let eps_ff = match &eps_tensors_rel {
        Some(tr) => FarFieldPermittivity::Tensor(tr),
        None => FarFieldPermittivity::Scalar(&eps),
    };
    let need_rt = spec.flux_bounds != "none" || (spec.adapt_passes > 1 && spec.adapt_estimator == "prager_synge");
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !volt.contains_key(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid(m_work, &ref_labels));
                let eff_ff_order = if need_rt { 1 } else { *ff_order };
                let k_far =
                    assemble_farfield_robin(m_work, &ff_labs, order, eff_ff_order, *decay_power, c, EPS0, eps_ff);
                k = k.add_scaled(&k_far, 1.0);
                if !source.is_uniform() {
                    let l = assemble_boundary_load_field::<D, S>(m_work, &ff_labs, order, source, 0.0)?;
                    f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                }
            }
        }
    }
    if !ie_exp_facets.is_empty() {
        if let (Some(mixed), Some(dof_map), Some(ie_map)) = (&mixed_opt, &dof_map_opt, &ie_map_opt) {
            let k_ie_rob = assemble_ie_outer_robin_mixed(mixed, dof_map, &ie_exp_facets, ie_map, EPS0, eps_ff);
            k = k.add_scaled(&k_ie_rob, 1.0);
            has_open = true;
        }
    }
    let (mut fixed, mut vals) = if source.is_uniform() {
        dirichlet_constraints(m_work, &volt, order)
    } else {
        dirichlet_constraints_field::<D, S>(m_work, &dirichlet_labels, order, source, 0.0)
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
                    pairs.push((d, 0.0));
                }
                pairs.sort_by_key(|p| p.0);
                fixed = pairs.iter().map(|p| p.0).collect();
                vals = pairs.iter().map(|p| p.1).collect();
            }
        }
    }
    Ok(ElectroAssembly { k, f, fixed, vals, volt, has_open })
}

pub fn assemble_electrostatics_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<ElectroAssembly, String> {
    let (qp_mesh, _) = prepare_electro_mesh(mesh, spec);
    let m_work = qp_mesh.as_ref().map_or(mesh, |(m, _, _)| m);
    assemble_electrostatics_on_work_mesh_with_source::<D, S>(m_work, spec, source, &[])
}

pub fn assemble_electrostatics(mesh: &Discretisation, spec: &Problem) -> Result<ElectroAssembly, String> {
    let (qp_mesh, _) = prepare_electro_mesh(mesh, spec);
    let m_work = qp_mesh.as_ref().map_or(mesh, |(m, _, _)| m);
    match m_work.dim {
        2 => assemble_electrostatics_on_work_mesh_with_source::<2, UniformSource>(m_work, spec, &UniformSource, &[]),
        3 => assemble_electrostatics_on_work_mesh_with_source::<3, UniformSource>(m_work, spec, &UniformSource, &[]),
        d => Err(format!("unsupported mesh dimension {d} in electrostatics")),
    }
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
    let sheet_thickness_m: BTreeMap<i32, f64> =
        spec.sheet_thickness_m.iter().filter(|&(_, &t)| t > 0.0 && t.is_finite()).map(|(&k, &v)| (k, v)).collect();
    TlsMeta {
        domain_eps_r,
        domain_names: spec.domain_names.clone(),
        conductor_labels,
        sheet_labels: sheet_set.into_iter().collect(),
        sheet_thickness_m,
        pec_labels: pec_set.into_iter().collect(),
        interfaces: spec.tls_interfaces.clone(),
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

#[derive(Clone, Copy)]
struct FacetFarFieldSpec {
    cell: usize,
    gamma: f64,
    center: [f64; 3],
}

fn build_facet_farfield_specs(
    mesh: &Discretisation,
    spec: &Problem,
    volt: &BTreeMap<i32, f64>,
) -> Vec<Option<FacetFarFieldSpec>> {
    let nf = mesh.n_facets();
    let has_ff = spec.bcs.iter().any(|bc| matches!(bc.kind, BcKind::FarField { .. } | BcKind::Open { .. }));
    if !has_ff || nf == 0 {
        return Vec::new();
    }
    let (lo, hi) = mesh.bbox();
    let default_origin = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
    let ref_labels: Vec<i32> = volt.keys().copied().collect();
    let mut label_ff: BTreeMap<i32, (f64, [f64; 3])> = BTreeMap::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::FarField { decay_power, center, .. } => {
                let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !volt.contains_key(l)).collect();
                if !ff_labs.is_empty() {
                    let c = center.unwrap_or_else(|| farfield_centroid(mesh, &ref_labels));
                    let gamma = if decay_power.is_finite() && *decay_power > 0.0 { *decay_power } else { 1.0 };
                    for l in ff_labs {
                        label_ff.insert(l, (gamma, c));
                    }
                }
            }
            BcKind::Open { origin } if mesh.dim == 3 => {
                let c = origin.unwrap_or(default_origin);
                for &l in &bc.labels {
                    if !volt.contains_key(&l) {
                        label_ff.insert(l, (1.0, c));
                    }
                }
            }
            _ => {}
        }
    }
    if label_ff.is_empty() {
        return Vec::new();
    }

    let sorted_key = |vs: &[u32]| -> [u32; 3] {
        let mut k = [u32::MAX; 3];
        for (i, &v) in vs.iter().take(3).enumerate() {
            k[i] = v;
        }
        k[..vs.len().min(3)].sort_unstable();
        k
    };
    let mut face_to_cell: BTreeMap<[u32; 3], usize> = BTreeMap::new();
    for c in 0..mesh.n_cells() {
        let cv = mesh.cell(c);
        if mesh.dim == 2 {
            for &(a, b) in &[(0, 1), (0, 2), (1, 2)] {
                face_to_cell.entry(sorted_key(&[cv[a], cv[b]])).or_insert(c);
            }
        } else {
            for &(a, b, d) in &[(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)] {
                face_to_cell.entry(sorted_key(&[cv[a], cv[b], cv[d]])).or_insert(c);
            }
        }
    }

    let mut out = vec![None; nf];
    for f in 0..nf {
        if let Some(&(gamma, center)) = label_ff.get(&mesh.facet_labels[f]) {
            if let Some(&cell) = face_to_cell.get(&sorted_key(mesh.facet(f))) {
                out[f] = Some(FacetFarFieldSpec { cell, gamma, center });
            }
        }
    }
    out
}

fn eval_farfield_alpha_pointwise(
    dim: usize,
    r_inf_2d: f64,
    ff: FacetFarFieldSpec,
    eps_rel: &[f64],
    eps_tensors_rel: Option<&[M3]>,
    x_q: [f64; 3],
    n_hat: [f64; 3],
) -> f64 {
    let eps_nf = if let Some(tr) = eps_tensors_rel {
        let m = &tr[ff.cell];
        let mut s = 0.0;
        for i in 0..dim.min(3) {
            for j in 0..dim.min(3) {
                s += n_hat[i] * m[i][j] * n_hat[j];
            }
        }
        if s.is_finite() {
            s.max(0.0)
        } else {
            0.0
        }
    } else {
        let e = eps_rel.get(ff.cell).copied().unwrap_or(1.0);
        if e.is_finite() {
            e.max(0.0)
        } else {
            0.0
        }
    };
    if eps_nf <= 0.0 {
        return 0.0;
    }
    if dim == 2 {
        let rx = x_q[0] - ff.center[0];
        let ry = x_q[1] - ff.center[1];
        let r2 = rx * rx + ry * ry;
        let rn = rx * n_hat[0] + ry * n_hat[1];
        if r2 <= 1e-30 || rn <= 0.0 {
            return 0.0;
        }
        let r = r2.sqrt();
        let r_scale = r_inf_2d.max(r);
        let denom_ln = (2.0 * r_scale / r).ln().max(1.0);
        EPS0 * eps_nf * (ff.gamma * (rn / r2)) / (2.0 * denom_ln)
    } else {
        let rx = x_q[0] - ff.center[0];
        let ry = x_q[1] - ff.center[1];
        let rz = x_q[2] - ff.center[2];
        let r2 = rx * rx + ry * ry + rz * rz;
        let rn = rx * n_hat[0] + ry * n_hat[1] + rz * n_hat[2];
        if r2 <= 1e-30 || rn <= 0.0 {
            return 0.0;
        }
        EPS0 * eps_nf * (ff.gamma * rn / r2)
    }
}

fn build_facet_midpoint_robin_alpha(
    mesh: &Discretisation,
    ff_specs: &[Option<FacetFarFieldSpec>],
    eps_rel: &[f64],
    eps_tensors_rel: Option<&[M3]>,
) -> Option<Vec<f64>> {
    if ff_specs.is_empty() {
        return None;
    }
    let nf = mesh.n_facets();
    let d = mesh.dim;
    let (lo, hi) = mesh.bbox();
    let r_inf_2d = if d == 2 { (hi[0] - lo[0]).hypot(hi[1] - lo[1]).max(1e-12) } else { 1.0 };
    let mut alphas = vec![0.0; nf];
    for (f, spec_opt) in ff_specs.iter().copied().enumerate() {
        let Some(ff) = spec_opt else { continue };
        let fv = mesh.facet(f);
        let cv = mesh.cell(ff.cell);
        let opp_v = cv.iter().copied().find(|v| !fv.contains(v)).unwrap_or(cv[0]);
        let p_opp_s = mesh.point(opp_v as usize);
        let mut p_opp = [0.0; 3];
        p_opp[..d].copy_from_slice(&p_opp_s[..d]);

        if d == 2 {
            let p0 = mesh.point(fv[0] as usize);
            let p1 = mesh.point(fv[1] as usize);
            let xm = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]), 0.0];
            let dx = p1[0] - p0[0];
            let dy = p1[1] - p0[1];
            let len = dx.hypot(dy);
            if len <= 1e-30 {
                continue;
            }
            let mut n_hat = [dy / len, -dx / len, 0.0];
            if (xm[0] - p_opp[0]) * n_hat[0] + (xm[1] - p_opp[1]) * n_hat[1] < 0.0 {
                n_hat[0] = -n_hat[0];
                n_hat[1] = -n_hat[1];
            }
            alphas[f] = eval_farfield_alpha_pointwise(2, r_inf_2d, ff, eps_rel, eps_tensors_rel, xm, n_hat);
        } else {
            let p0 = mesh.point(fv[0] as usize);
            let p1 = mesh.point(fv[1] as usize);
            let p2 = mesh.point(fv[2] as usize);
            let xm = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
            let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            let len = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
            if len <= 1e-30 {
                continue;
            }
            let mut n_hat = [cr[0] / len, cr[1] / len, cr[2] / len];
            if (xm[0] - p_opp[0]) * n_hat[0] + (xm[1] - p_opp[1]) * n_hat[1] + (xm[2] - p_opp[2]) * n_hat[2] < 0.0 {
                n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
            }
            alphas[f] = eval_farfield_alpha_pointwise(3, 1.0, ff, eps_rel, eps_tensors_rel, xm, n_hat);
        }
    }
    Some(alphas)
}

fn compute_hessian_cell_eta2(mesh: &Discretisation, potentials: &[Vec<f64>]) -> Result<(Vec<f64>, f64), String> {
    let nc = mesh.n_cells();
    if nc == 0 || potentials.is_empty() {
        return Ok((vec![0.0; nc], 0.0));
    }
    let geo = geometry(mesh);
    let mut eta2 = vec![0.0; nc];
    if mesh.dim == 2 {
        let hs_all: Vec<Vec<[[f64; 2]; 2]>> = potentials
            .iter()
            .map(|v| femcore::adapt::estimate::recovered_hessians(mesh, v))
            .collect::<Result<Vec<_>, _>>()?;
        for c in 0..nc {
            let area = geo.vol[c].max(0.0);
            let h_c = (4.0 * area / 3.0_f64.sqrt()).sqrt();
            let cv = mesh.cell(c);
            let mut frob_sq = 0.0;
            for hs in &hs_all {
                let mut m_avg = [[0.0; 2]; 2];
                for &vi in cv {
                    let mv = &hs[vi as usize];
                    for r in 0..2 {
                        for col in 0..2 {
                            m_avg[r][col] += mv[r][col] / 3.0;
                        }
                    }
                }
                frob_sq += m_avg[0][0] * m_avg[0][0] + 2.0 * m_avg[0][1] * m_avg[0][1] + m_avg[1][1] * m_avg[1][1];
            }
            eta2[c] = area * h_c.powi(4) * frob_sq;
        }
    } else {
        let hs_all: Vec<Vec<[[f64; 3]; 3]>> = potentials
            .iter()
            .map(|v| femcore::adapt::estimate3::recovered_hessians(mesh, v))
            .collect::<Result<Vec<_>, _>>()?;
        for c in 0..nc {
            let vol = geo.vol[c].max(0.0);
            let h_c = (6.0 * 2.0_f64.sqrt() * vol).cbrt();
            let cv = mesh.cell(c);
            let mut frob_sq = 0.0;
            for hs in &hs_all {
                let mut m_avg = [[0.0; 3]; 3];
                for &vi in cv {
                    let mv = &hs[vi as usize];
                    for r in 0..3 {
                        for col in 0..3 {
                            m_avg[r][col] += 0.25 * mv[r][col];
                        }
                    }
                }
                for r in 0..3 {
                    for col in 0..3 {
                        frob_sq += m_avg[r][col] * m_avg[r][col];
                    }
                }
            }
            eta2[c] = vol * h_c.powi(4) * frob_sq;
        }
    }
    let total = eta2.iter().sum::<f64>().sqrt();
    Ok((eta2, total))
}

#[allow(clippy::too_many_arguments)]
fn compute_amr_indicators_on_work_mesh(
    m_work: &Discretisation,
    spec: &Problem,
    volt: &BTreeMap<i32, f64>,
    eps_rel: &[f64],
    eps_tensors_rel: Option<&[M3]>,
    ff_specs: &[Option<FacetFarFieldSpec>],
    potentials: &[Vec<f64>],
    cap_matrix_opt: Option<&[Vec<f64>]>,
    eq_opt: Option<&femcore::rt::EquilibratedFluxResult>,
) -> Result<ElectroAmrData, String> {
    match spec.adapt_estimator.as_str() {
        "hessian" => {
            let (eta2, total_eta) = compute_hessian_cell_eta2(m_work, potentials)?;
            Ok(ElectroAmrData { eta2, total_eta })
        }
        "prager_synge" => {
            if let Some(eq) = eq_opt {
                let total_eta = eq.cell_eta2.iter().sum::<f64>().sqrt();
                return Ok(ElectroAmrData { eta2: eq.cell_eta2.clone(), total_eta });
            }
            let dir_labels: Vec<i32> = volt.keys().copied().collect();
            let eps_phys: Vec<f64> = eps_rel.iter().map(|e| EPS0 * e).collect();
            let eps_t_phys: Option<Vec<M3>> = eps_tensors_rel.map(|tr| {
                let mut t = tr.to_vec();
                t.iter_mut().for_each(|m| m.iter_mut().for_each(|r| r.iter_mut().for_each(|x| *x *= EPS0)));
                t
            });
            let cell_perm = match &eps_t_phys {
                Some(t) => femcore::rt::CellPermittivity::Tensor(t),
                None => femcore::rt::CellPermittivity::Scalar(&eps_phys),
            };
            let (lo, hi) = m_work.bbox();
            let r_inf_2d = if m_work.dim == 2 { (hi[0] - lo[0]).hypot(hi[1] - lo[1]).max(1e-12) } else { 1.0 };
            let pw_fn = |f: usize, xp: [f64; 3], n_out: [f64; 3]| -> f64 {
                ff_specs.get(f).and_then(|x| *x).map_or(0.0, |ff| {
                    eval_farfield_alpha_pointwise(m_work.dim, r_inf_2d, ff, eps_rel, eps_tensors_rel, xp, n_out)
                })
            };
            let robin = if ff_specs.is_empty() {
                femcore::rt::RobinBoundary::None
            } else {
                femcore::rt::RobinBoundary::Pointwise(&pw_fn)
            };
            let eq = femcore::rt::equilibrate_electro_flux_bounds(
                m_work,
                spec.order,
                potentials,
                cell_perm,
                &dir_labels,
                robin,
                cap_matrix_opt,
            )?;
            let total_eta = eq.cell_eta2.iter().sum::<f64>().sqrt();
            Ok(ElectroAmrData { eta2: eq.cell_eta2, total_eta })
        }
        _ => {
            let nf = m_work.n_facets();
            let dirichlet_facets: Vec<bool> = (0..nf).map(|f| volt.contains_key(&m_work.facet_labels[f])).collect();
            let perm = match eps_tensors_rel {
                Some(t) => femcore::adapt::electro::Permittivity::Tensor(t),
                None => femcore::adapt::electro::Permittivity::Scalar(eps_rel),
            };
            let robin_alpha = build_facet_midpoint_robin_alpha(m_work, ff_specs, eps_rel, eps_tensors_rel);
            let mut opt = femcore::adapt::electro::ElectroEstimatorOptions::new(spec.order, perm, &dirichlet_facets);
            opt.robin_alpha = robin_alpha.as_deref();
            if let Some(cap) = cap_matrix_opt {
                let goal = femcore::adapt::electro::electro_goal_indicators(
                    m_work,
                    potentials,
                    cap,
                    femcore::adapt::electro::DEFAULT_EPS_FLOOR,
                    &opt,
                )?;
                Ok(ElectroAmrData { eta2: goal.eta2, total_eta: goal.total })
            } else {
                let terms = femcore::adapt::electro::electro_term_indicators(m_work, potentials, &opt)?;
                let eta2 = terms.eta2_by_term.into_iter().next().unwrap_or_else(|| vec![0.0; m_work.n_cells()]);
                let total_eta = eta2.iter().sum::<f64>().sqrt();
                Ok(ElectroAmrData { eta2, total_eta })
            }
        }
    }
}

pub fn solve_electrostatics(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    match mesh.dim {
        2 => solve_electrostatics_with_source_and_robin::<2, UniformSource>(mesh, spec, &UniformSource, &[]),
        3 => solve_electrostatics_with_source_and_robin::<3, UniformSource>(mesh, spec, &UniformSource, &[]),
        d => Err(format!("unsupported mesh dimension {d} in electrostatics")),
    }
}

pub fn solve_electrostatics_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<SolveResult, String> {
    solve_electrostatics_with_source_and_robin::<D, S>(mesh, spec, source, &[])
}

pub fn solve_electrostatics_with_source_and_robin<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
    robin_bcs: &[(Vec<i32>, f64)],
) -> Result<SolveResult, String> {
    match spec.formulation.as_str() {
        "" | "volumetric_fem" => {}
        "planar_bem" => {
            if !source.is_uniform() || !robin_bcs.is_empty() {
                return Err("planar_bem does not support volumetric SourceField or explicit Robin BCs".into());
            }
            return crate::planar_bem::solve_planar_bem(mesh, spec);
        }
        "auto" => {
            if source.is_uniform() && robin_bcs.is_empty() && crate::planar_bem::can_use_planar_bem(mesh, spec) {
                return crate::planar_bem::solve_planar_bem(mesh, spec);
            }
        }
        other => return Err(format!("unknown electrostatics formulation '{other}'")),
    }
    let t0 = Instant::now();
    let order = spec.order;
    let (qp_info, already_qp) = prepare_electro_mesh(mesh, spec);
    let m_work = qp_info.as_ref().map_or(mesh, |(m, _, _)| m);

    let need_rt = spec.flux_bounds != "none" || (spec.adapt_passes > 1 && spec.adapt_estimator == "prager_synge");
    if need_rt {
        if spec.singular_edges == "quarter_point" || (order >= 2 && m_work.has_qp()) {
            return Err(
                "flux_bounds / prager_synge cannot be combined with singular_edges='quarter_point': contravariant Piola-mapped Raviart-Thomas fluxes have divergent L^2 energy on singular quarter-point (det J -> 0) CurvedTets".into(),
            );
        }
        let rho = spec.materials.cell_values(m_work, "rho_v", Some(0.0))?;
        let has_sigma_s = spec.bcs.iter().any(|bc| matches!(bc.kind, BcKind::SurfaceCharge(s) if s.abs() > 0.0));
        if rho.iter().any(|&x| x.abs() > 0.0) || has_sigma_s || !source.is_uniform() {
            return Err(
                "flux_bounds / prager_synge requires source-free electrostatics (rho_v == 0 and surface_charge == 0)"
                    .into(),
            );
        }
    }

    let asm = assemble_electrostatics_on_work_mesh_with_source::<D, S>(m_work, spec, source, robin_bcs)?;
    let eps_rel = spec.materials.cell_values(m_work, "eps_r", Some(1.0))?;
    let eps_tensors_rel: Option<Vec<M3>> = if spec.materials.has_tensor(TensorKind::Permittivity) {
        Some(plane_tensors(m_work, spec.materials.cell_tensor(m_work, TensorKind::Permittivity, Some(1.0))?))
    } else {
        None
    };
    let ff_specs = build_facet_farfield_specs(m_work, spec, &asm.volt);

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

    let tls = build_tls_meta(m_work, spec);
    let ref_labels: Vec<i32> = asm.volt.keys().copied().collect();
    let has_positive_farfield = spec.bcs.iter().any(|bc| {
        if let BcKind::FarField { center, .. } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !asm.volt.contains_key(l)).collect();
            if ff_labs.is_empty() {
                return false;
            }
            let c = center.unwrap_or_else(|| farfield_centroid(m_work, &ref_labels));
            farfield_positive_measure(m_work, &ff_labs, c) > 0.0
        } else {
            false
        }
    });

    if terminals.is_empty() {
        let (v, solver_note) = solve_potential(m_work, spec, &asm)?;
        let (fields, cell_fields) = electrostatics_fields(m_work, order, &v);
        let en = norm_of(&components(&fields[1].values, m_work.dim));
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
        } else if dv.len() == 1 && (asm.has_open || has_positive_farfield) && dv[0].abs() > 1e-12 {
            scalars.push((
                "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(),
                Scalar::Num(2.0 * w / dv[0].powi(2)),
            ));
        }
        if spec.adapt_passes > 1 {
            let c_single = vec![vec![(2.0 * w).max(1e-30)]];
            let amr = compute_amr_indicators_on_work_mesh(
                m_work,
                spec,
                &asm.volt,
                &eps_rel,
                eps_tensors_rel.as_deref(),
                &ff_specs,
                std::slice::from_ref(&v),
                Some(&c_single),
                None,
            )?;
            scalars.push(amr.scalar());
        }
        scalars.push(tls.scalar());
        return Ok(SolveResult {
            physics: "dfs".into(),
            study: spec.study.kind().into(),
            order: order as u8,
            dim: m_work.dim as u8,
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
            ..Default::default()
        });
    }

    // Multi-RHS terminal capacitance matrix extraction (N >= 1 terminals):
    let n_term = terminals.len();
    let n = asm.k.nrows();
    let term_dofs: Vec<BTreeSet<usize>> =
        terminals.iter().map(|(_, labs)| boundary_dofs(m_work, labs, order).into_iter().collect()).collect();
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
    let mut bg_dofs: BTreeSet<usize> = boundary_dofs(m_work, &bg_labels, order).into_iter().collect();
    for &d in &asm.fixed {
        if !term_dofs.iter().any(|td| td.contains(&d)) {
            bg_dofs.insert(d);
        }
    }
    if n_term == 1 && bg_dofs.is_empty() && !asm.has_open && !has_positive_farfield {
        return Err(
            "electrostatics with a single terminal needs at least one ground, potential, open, or far_field reference boundary".into()
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
        LinearSolver::Auto => amg_wanted(m_work.dim, n_free),
    } && !fixed.is_empty()
        && n_free > 0;

    let mut log = Vec::new();
    if let Some((_, n_rim, n_edges)) = &qp_info {
        if *n_edges > 0 {
            log.push(format!(
                "singular_edges=quarter_point: shifted {n_edges} tet edges meeting {n_rim} conductor sheet perimeter vertices to Barsoum 1/4-points (all affected cells verified det(J) > 0)"
            ));
        }
    } else if already_qp {
        log.push(format!(
            "singular_edges=quarter_point: using {} pre-enriched quarter-point tet edges from mesher",
            m_work.qp_edge_list().len()
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
            let sys = AmgSystem::new_with(m_work, order, 1, &asm.k, &fixed, opts)?;
            let base_tol = if spec.tol > 0.0 { spec.tol } else { 1e-10 };
            let tol = if need_rt { base_tol.min(1e-11) } else { base_tol };
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
    let mut c_raw = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in 0..n_term {
            c_raw[i][j] = femcore::linalg::vecops::dotu(&v_all[i], &kv_all[j]);
        }
    }
    let raw_asym = if n_term <= 1 {
        0.0
    } else {
        let mut diff_sq = 0.0;
        let mut norm_sq = 0.0;
        for i in 0..n_term {
            for j in 0..n_term {
                let d = c_raw[i][j] - c_raw[j][i];
                diff_sq += d * d;
                norm_sq += c_raw[i][j] * c_raw[i][j];
            }
        }
        if norm_sq > 0.0 {
            (diff_sq / norm_sq).sqrt()
        } else {
            0.0
        }
    };
    let mut maxwell_matrix_f = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in i..n_term {
            let sym = 0.5 * (c_raw[i][j] + c_raw[j][i]);
            maxwell_matrix_f[i][j] = sym;
            maxwell_matrix_f[j][i] = sym;
        }
    }

    // AF-4: Prager-Synge H(div) Raviart-Thomas complementary flux equilibration & Loewner bounds
    let mut maxwell_lower_matrix_f = None;
    let mut maxwell_upper_matrix_f = None;
    let mut maxwell_mid_matrix_f = None;
    let mut maxwell_error_bound_f = None;
    let eq_opt = if need_rt {
        let dir_labels: Vec<i32> = asm.volt.keys().copied().collect();
        let eps_phys: Vec<f64> = eps_rel.iter().map(|e| EPS0 * e).collect();
        let eps_t_phys: Option<Vec<M3>> = eps_tensors_rel.as_ref().map(|tr| {
            let mut t = tr.clone();
            t.iter_mut().for_each(|m| m.iter_mut().for_each(|r| r.iter_mut().for_each(|x| *x *= EPS0)));
            t
        });
        let cell_perm = match &eps_t_phys {
            Some(t) => femcore::rt::CellPermittivity::Tensor(t),
            None => femcore::rt::CellPermittivity::Scalar(&eps_phys),
        };
        let (lo, hi) = m_work.bbox();
        let r_inf_2d = if m_work.dim == 2 { (hi[0] - lo[0]).hypot(hi[1] - lo[1]).max(1e-12) } else { 1.0 };
        let pw_fn = |f: usize, xp: [f64; 3], n_out: [f64; 3]| -> f64 {
            ff_specs.get(f).and_then(|x| *x).map_or(0.0, |ff| {
                eval_farfield_alpha_pointwise(m_work.dim, r_inf_2d, ff, &eps_rel, eps_tensors_rel.as_deref(), xp, n_out)
            })
        };
        let robin = if ff_specs.is_empty() {
            femcore::rt::RobinBoundary::None
        } else {
            femcore::rt::RobinBoundary::Pointwise(&pw_fn)
        };
        let eq = femcore::rt::equilibrate_electro_flux_bounds(
            m_work,
            order,
            &v_all,
            cell_perm,
            &dir_labels,
            robin,
            Some(&maxwell_matrix_f),
        )?;
        if spec.flux_bounds != "none" {
            log.push(format!(
                "flux_bounds={}: RT_{} equilibration solved {} vertex half-star patches ({} sheet splits, max ||div D_tilde||_L2 = {:.2e})",
                spec.flux_bounds, eq.rt_order, eq.n_patches, eq.n_sheet_splits, eq.max_div_l2
            ));
            if spec.flux_bounds == "midpoint" {
                maxwell_matrix_f = eq.c_mid_matrix.clone();
            }
            maxwell_lower_matrix_f = Some(eq.c_lower_matrix.clone());
            maxwell_upper_matrix_f = Some(eq.c_upper_matrix.clone());
            maxwell_mid_matrix_f = Some(eq.c_mid_matrix.clone());
            maxwell_error_bound_f = Some(eq.c_error_bound_matrix.clone());
        }
        Some(eq)
    } else {
        None
    };

    // AF-6: Quasi-electrostatic AC dielectric conductance, Q matrix, and strong-loss Delta C^(2) >= 0
    let ac = crate::ac_loss::compute_ac_loss(
        m_work,
        spec,
        (spec.freq > 0.0).then_some(spec.freq),
        &asm.k,
        &free,
        &v_all,
        &maxwell_matrix_f,
    )?;
    maxwell_matrix_f = ac.corrected_maxwell_matrix_f;
    if ac.conductance_matrix_s.is_some() || ac.q_dielectric_matrix.is_some() {
        log.push(format!(
            "ac_loss: evaluated dielectric conductance/Q matrices (freq = {:.3e} Hz, max tan_delta_eff = {:.3e}{})",
            spec.freq,
            ac.max_tan_delta_eff,
            if ac.delta_c_second_order_f.is_some() { ", applied complex Schur Delta C" } else { "" }
        ));
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
        let (fields, cell_fields) = electrostatics_fields(m_work, order, &v_all[j]);
        if j == 0 {
            let en = norm_of(&components(&fields[1].values, m_work.dim));
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

    scalars.push(("maxwell_raw_asymmetry_rel".into(), Scalar::Num(raw_asym)));

    if spec.adapt_passes > 1 {
        let amr = compute_amr_indicators_on_work_mesh(
            m_work,
            spec,
            &asm.volt,
            &eps_rel,
            eps_tensors_rel.as_deref(),
            &ff_specs,
            &v_all,
            Some(&maxwell_matrix_f),
            eq_opt.as_ref(),
        )?;
        scalars.push(amr.scalar());
    }
    let cap_data = CapacitanceMatrixData {
        terminals: terminals.into_iter().map(|(n, _)| n).collect(),
        maxwell_matrix_f,
        mutual_matrix_f,
        ec_matrix_hz,
        maxwell_lower_matrix_f,
        maxwell_upper_matrix_f,
        maxwell_mid_matrix_f,
        maxwell_error_bound_f,
        conductance_matrix_s: ac.conductance_matrix_s,
        q_dielectric_matrix: ac.q_dielectric_matrix,
        ..Default::default()
    };
    scalars.push(cap_data.scalar());
    scalars.push(tls.scalar());

    Ok(SolveResult {
        physics: "dfs".into(),
        study: spec.study.kind().into(),
        order: order as u8,
        dim: m_work.dim as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        ..Default::default()
    })
}

fn electrostatics_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    order: usize,
    v_mixed: &[f64],
) -> (Vec<Field>, CellFields, Vec<f64>) {
    let v_simplex = proj.project_lagrange_blocks(mixed, dof_map, v_mixed, 1, order);
    let scale = vec![-1.0; mixed.n_cells()];
    let ecell = sample_mixed_gradients_on_simplex(mixed, dof_map, proj, v_mixed, None, &scale, order);
    let e = node_average(simplex, &ecell, simplex.dim, order);
    let ec = components(&e, simplex.dim);
    let names = femcore::section::component_names(simplex.dim);
    let fields = vec![Field::scalar("V", "V", v_simplex.clone()), Field::vector("E", "V/m", names, &ec)];
    let mut cf = CellFields::new(simplex, order);
    cf.push_real("V", "V", cell_values_at_nodes(simplex, &v_simplex, order));
    cf.push_components("E", names, "V/m", &ecell);
    (fields, cf, v_simplex)
}

/// Solves electrostatics natively on a mixed mesh (`Tri`/`Quad` in 2D, `Tet`/`Prism`/`Hex` in 3D)
/// and/or with domain-applied static Infinite Elements (`BcKind::InfiniteElement`), preserving the
/// pure-simplex fast path when `mesh.is_pure_simplex() && !spec.has_infinite_element()`.
pub fn solve_electrostatics_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String> {
    if mesh.is_pure_simplex() && !spec.has_infinite_element() {
        let d = discretise(mesh).map_err(|e| e.to_string())?;
        return solve_electrostatics(&d, spec);
    }
    match spec.formulation.as_str() {
        "" | "volumetric_fem" | "auto" => {}
        "planar_bem" => {
            let d = discretise(mesh).map_err(|e| e.to_string())?;
            return crate::planar_bem::solve_planar_bem(&d, spec);
        }
        other => return Err(format!("unknown electrostatics formulation '{other}'")),
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
    let eps = cell_values_mixed(&spec.materials, &mixed, "eps_r", Some(1.0))?;
    let rho = cell_values_mixed(&spec.materials, &mixed, "rho_v", Some(0.0))?;
    let coeff: Vec<f64> = eps.iter().map(|e| EPS0 * e).collect();
    let eps_tensors_rel: Option<Vec<M3>> = if spec.materials.has_tensor(TensorKind::Permittivity) {
        Some(plane_tensors_mixed(
            &mixed,
            cell_tensor_mixed(&spec.materials, &mixed, TensorKind::Permittivity, Some(1.0))?,
        ))
    } else {
        None
    };
    let tensor_perm: Option<Vec<M3>> = eps_tensors_rel.as_ref().map(|tr| {
        let mut t = tr.clone();
        t.iter_mut().for_each(|m| m.iter_mut().for_each(|r| r.iter_mut().for_each(|x| *x *= EPS0)));
        t
    });

    let mut k = assemble_laplace_ie_mixed(&mixed, &dof_map, &coeff, tensor_perm.as_deref(), &ie_map);
    let n = dof_map.n_dofs;
    let mut f =
        if rho.iter().any(|&r| r != 0.0) { assemble_source_mixed(&mixed, &dof_map, &rho) } else { vec![0.0; n] };

    let mut volt: BTreeMap<i32, f64> = BTreeMap::new();
    let (lo, hi) = mixed.bbox();
    let default_origin = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
    let mut open_origins: BTreeMap<i32, [f64; 3]> = BTreeMap::new();
    let mut terminals: Vec<(String, Vec<i32>)> = Vec::new();

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
                BcKind::Open { origin } => {
                    if mixed.dim != 3 {
                        return Err(
                            "open boundary condition (spherical 1/r asymptotic Robin) requires a 3D mesh".into()
                        );
                    }
                    let r0 = origin.unwrap_or(default_origin);
                    if r0.iter().any(|x| !x.is_finite()) {
                        return Err("open boundary 'origin' must be finite".into());
                    }
                    open_origins.insert(lab, r0);
                }
                BcKind::SurfaceCharge(s) => {
                    if *s != 0.0 {
                        let l = assemble_boundary_load_mixed(&mixed, &dof_map, &[lab], *s);
                        f.iter_mut().zip(&l).for_each(|(a, b)| *a += b);
                    }
                }
                BcKind::FarField { .. } | BcKind::InfiniteElement(_) => {}
                other => return Err(format!("'{}' is not an electrostatics boundary condition", other.name())),
            }
        }
    }

    let eps_ff = match &eps_tensors_rel {
        Some(tr) => FarFieldPermittivity::Tensor(tr),
        None => FarFieldPermittivity::Scalar(&eps),
    };
    let mut has_open = !open_origins.is_empty();
    for (&lab, &r0) in &open_origins {
        let k_open = assemble_farfield_robin_mixed(&mixed, &dof_map, &[lab], 1, 1.0, r0, EPS0, eps_ff);
        k = k.add_scaled(&k_open, 1.0);
    }

    let ref_labels: Vec<i32> = volt.keys().copied().collect();
    let mut has_positive_farfield = false;
    for bc in &spec.bcs {
        if let BcKind::FarField { order: ff_order, decay_power, center } = &bc.kind {
            let ff_labs: Vec<i32> = bc.labels.iter().copied().filter(|l| !volt.contains_key(l)).collect();
            if !ff_labs.is_empty() {
                let c = center.unwrap_or_else(|| farfield_centroid_mixed(&mixed, &ref_labels));
                if farfield_positive_measure_mixed(&mixed, &ff_labs, c) > 0.0 {
                    has_positive_farfield = true;
                }
                let k_far =
                    assemble_farfield_robin_mixed(&mixed, &dof_map, &ff_labs, *ff_order, *decay_power, c, EPS0, eps_ff);
                k = k.add_scaled(&k_far, 1.0);
            }
        }
    }

    let (ie_rat_facets, ie_exp_facets) = classify_ie_outer_facets_mixed(&mixed, &ie_map);
    if !ie_exp_facets.is_empty() {
        let k_ie_rob = assemble_ie_outer_robin_mixed(&mixed, &dof_map, &ie_exp_facets, &ie_map, EPS0, eps_ff);
        k = k.add_scaled(&k_ie_rob, 1.0);
        has_open = true;
    }
    if !ie_rat_facets.is_empty() {
        has_open = true;
    }

    if volt.is_empty() && !has_open && !has_positive_farfield {
        return Err("electrostatics needs at least one potential, ground or terminal boundary".into());
    }

    let mut ie_zero_dofs = BTreeSet::new();
    if !ie_rat_facets.is_empty() {
        let mut fbuf = Vec::new();
        for &f_idx in &ie_rat_facets {
            dof_map.facet_dofs(&mixed, f_idx, &mut fbuf);
            ie_zero_dofs.extend(fbuf.iter().copied());
        }
    }

    let tls = build_tls_meta(&simplex, spec);

    if terminals.is_empty() {
        let (mut fixed, mut vals) = dof_map.dirichlet_constraints(&mixed, &volt);
        if !ie_zero_dofs.is_empty() {
            let fixed_set: BTreeSet<usize> = fixed.iter().copied().collect();
            let mut pairs: Vec<(usize, f64)> = fixed.into_iter().zip(vals).collect();
            for &d in &ie_zero_dofs {
                if !fixed_set.contains(&d) {
                    pairs.push((d, 0.0));
                }
            }
            pairs.sort_by_key(|p| p.0);
            fixed = pairs.iter().map(|p| p.0).collect();
            vals = pairs.iter().map(|p| p.1).collect();
        }
        let (v_mixed, solver_note) = solve_spd_mixed(spec.solver, &k, &f, &fixed, &vals)?;
        let (fields, cell_fields, v_simplex) =
            electrostatics_fields_mixed(&mixed, &dof_map, &simplex, &proj, order, &v_mixed);
        let en = norm_of(&components(&fields[1].values, mixed.dim));
        let kv = k.matvec(&v_mixed);
        let w = 0.5 * v_mixed.iter().zip(&kv).map(|(a, b)| a * b).sum::<f64>();
        let (vmin, vmax) = minmax(&v_simplex);
        let mut scalars = vec![
            ("stored energy W_e (J; per metre in 2D)".to_string(), Scalar::Num(w)),
            ("max |E| (V/m)".into(), Scalar::Num(minmax(&en).1)),
            ("V min (V)".into(), Scalar::Num(vmin)),
            ("V max (V)".into(), Scalar::Num(vmax)),
        ];
        let mut dv: Vec<f64> = volt.values().map(|x| (x * 1e12).round() / 1e12).collect();
        dv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        dv.dedup();
        if dv.len() == 2 {
            scalars.push((
                "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(),
                Scalar::Num(2.0 * w / (dv[1] - dv[0]).powi(2)),
            ));
        } else if dv.len() == 1 && (has_open || has_positive_farfield) && dv[0].abs() > 1e-12 {
            scalars.push((
                "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)".into(),
                Scalar::Num(2.0 * w / dv[0].powi(2)),
            ));
        }
        scalars.push(tls.scalar());
        return Ok(SolveResult {
            physics: "dfs".into(),
            study: spec.study.kind().into(),
            order: order as u8,
            dim: mixed.dim as u8,
            solutions: vec![Solution {
                label: "stationary".into(),
                fields,
                cell_fields: cell_fields_of(&cell_fields),
                ..Default::default()
            }],
            scalars,
            log: solver_note.into_iter().collect(),
            seconds: t0.elapsed().as_secs_f64(),
            dofs: n,
            ..Default::default()
        });
    }

    // Multi-RHS terminal capacitance matrix extraction on MixedDiscretisation:
    let n_term = terminals.len();
    let term_dofs: Vec<BTreeSet<usize>> =
        terminals.iter().map(|(_, labs)| dof_map.boundary_dofs(&mixed, labs).into_iter().collect()).collect();
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
    let mut bg_dofs: BTreeSet<usize> = dof_map.boundary_dofs(&mixed, &bg_labels).into_iter().collect();
    for &d in &ie_zero_dofs {
        if !term_dofs.iter().any(|td| td.contains(&d)) {
            bg_dofs.insert(d);
        }
    }
    if n_term == 1 && bg_dofs.is_empty() && !has_open && !has_positive_farfield {
        return Err(
            "electrostatics with a single terminal needs at least one ground, potential, open, or far_field reference boundary".into()
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

    let mut v_all = vec![vec![0.0; n]; n_term];
    for j in 0..n_term {
        for (&d, &val) in fixed.iter().zip(&vals_list[j]) {
            v_all[j][d] = val;
        }
    }

    let mut log = Vec::new();
    if n_free > 0 {
        let kff = k.submatrix(&free, &free);
        let kfx = k.submatrix(&free, &fixed);
        let rhs_list: Vec<Vec<f64>> = vals_list
            .iter()
            .map(|vals| {
                let kx = kfx.matvec(vals);
                free.iter().enumerate().map(|(a, &i)| f[i] - kx[a]).collect()
            })
            .collect();
        let solver = DirectSolver::new(&kff, true)?;
        let sols = solver.solve_many(&rhs_list);
        for j in 0..n_term {
            for (a, &i) in free.iter().enumerate() {
                v_all[j][i] = sols[j][a];
            }
        }
        log.push(format!(
            "direct Cholesky: 1 factorization ({n_free} free mixed-cell dofs) served {n_term} terminal right-hand side(s)"
        ));
    } else {
        log.push(format!("0 free dofs: {n_term} terminal right-hand side(s) determined by boundary constraints"));
    }

    let kv_all: Vec<Vec<f64>> = v_all.iter().map(|v| k.matvec(v)).collect();
    let mut c_raw = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in 0..n_term {
            c_raw[i][j] = femcore::linalg::vecops::dotu(&v_all[i], &kv_all[j]);
        }
    }
    let raw_asym = if n_term <= 1 {
        0.0
    } else {
        let mut diff_sq = 0.0;
        let mut norm_sq = 0.0;
        for i in 0..n_term {
            for j in 0..n_term {
                let d = c_raw[i][j] - c_raw[j][i];
                diff_sq += d * d;
                norm_sq += c_raw[i][j] * c_raw[i][j];
            }
        }
        if norm_sq > 0.0 {
            (diff_sq / norm_sq).sqrt()
        } else {
            0.0
        }
    };
    let mut maxwell_matrix_f = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        for j in i..n_term {
            let sym = 0.5 * (c_raw[i][j] + c_raw[j][i]);
            maxwell_matrix_f[i][j] = sym;
            maxwell_matrix_f[j][i] = sym;
        }
    }

    let mut mutual_matrix_f = vec![vec![0.0; n_term]; n_term];
    for i in 0..n_term {
        let mut sum_row = 0.0;
        for col in 0..n_term {
            sum_row += maxwell_matrix_f[i][col];
            if i != col {
                mutual_matrix_f[i][col] = -maxwell_matrix_f[i][col];
            }
        }
        mutual_matrix_f[i][i] = sum_row;
    }

    let ec_matrix_hz = invert_spd_ec_matrix(&maxwell_matrix_f)?;

    let mut solutions = Vec::with_capacity(n_term);
    let mut scalars = Vec::new();
    for j in 0..n_term {
        let (fields, cell_fields, v_simplex) =
            electrostatics_fields_mixed(&mixed, &dof_map, &simplex, &proj, order, &v_all[j]);
        if j == 0 {
            let en = norm_of(&components(&fields[1].values, mixed.dim));
            let w = 0.5 * maxwell_matrix_f[0][0];
            let (vmin, vmax) = minmax(&v_simplex);
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

    scalars.push(("maxwell_raw_asymmetry_rel".into(), Scalar::Num(raw_asym)));
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
        dim: mixed.dim as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        ..Default::default()
    })
}
