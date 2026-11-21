// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification tests for AF-2: Barsoum quarter-point (`xi_QP ~ 0.25`) `CurvedTet` sheet-edge singularity
//! enrichment, 216-point `quad_rule_collapsed(8)` positivity, curved-rim sagitta guard, singularity-safe
//! `gradients_at_nodes` (`0` `NaN`/`Inf` at `det J = 0` tips), and Morrison-Lewis / Babuska corner grading.

use femcore::fem::{quad_rule_collapsed, CurvedTet, TET_EDGES};
use femcore::fields::{gradients_at_nodes, node_average};
use femcore::mesh::SimplexMesh;
use femcore::meshing::curve3::{relax_edge_nodes, scaled_jacobian_min, scaled_jacobian_min_interior_q8};
use femcore::meshing::geomdrive::evaluate;
use femcore::meshing::request::{mesh_request, request_schema, CurvedSummary, MeshRequest, MesherRequest};
use femcore::meshing::singularity::{
    analytical_qp_jacobian_ratio, apply_quarter_point_nodes, collect_overlay_sheet_corners, corner_graded_h_at,
    corner_grading_h, coupled_corner_dz, measure_overlay_sheet_corners, measure_sheet_corner_angles_m2,
    morrison_lewis_exponent, sheet_rim_vertices, CORNER_H_MIN_FACTOR, CORNER_RADIUS_FACTOR,
};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::Model;
use std::collections::BTreeSet;
use std::f64::consts::PI;

fn pt2(x: f64, y: f64) -> [Expr; 2] {
    [Expr::Num(x), Expr::Num(y)]
}

/// Builds a single-tet `CurvedTet` with vertices `v`, marking vertices in `is_rim` as sheet-perimeter vertices
/// and shifting only edges `(i, j)` with `is_rim[i] ^ is_rim[j]` to `(1 - xi_qp) * v_rim + xi_qp * v_other`.
fn make_qp_tet(v: [[f64; 3]; 4], is_rim: [bool; 4], xi_qp: f64) -> CurvedTet {
    let mut x = CurvedTet::straight(v).x;
    for (e, &(i, j)) in TET_EDGES.iter().enumerate() {
        if is_rim[i] ^ is_rim[j] {
            let (p_rim, p_other) = if is_rim[i] { (v[i], v[j]) } else { (v[j], v[i]) };
            for k in 0..3 {
                x[4 + e][k] = (1.0 - xi_qp) * p_rim[k] + xi_qp * p_other[k];
            }
        }
    }
    CurvedTet::new(x)
}

#[test]
fn case_1_2_3_analytical_jacobian_ratios_match_curved_tet_exactly() {
    // Non-trivial skewed 3D tetrahedron
    let v = [
        [0.12, -0.18, 0.04],
        [1.35, 0.14, -0.08],
        [0.22, 1.15, 0.27],
        [-0.11, 0.24, 0.93],
    ];
    let det_s = CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
    assert!(det_s > 0.1);
    let q8 = quad_rule_collapsed(8);
    assert_eq!(q8.lam.len(), 216, "degree-8 collapsed rule must have 6^3 = 216 Gauss points");

    let cases: &[(usize, [bool; 4])] = &[
        (0, [false, false, false, false]),
        (1, [true, false, false, false]), // Case 2: |S_rim| = 1 (vertex-adjacent tet)
        (1, [false, false, true, false]),
        (2, [true, true, false, false]),  // Case 1: |S_rim| = 2 (edge-adjacent tet)
        (2, [false, true, false, true]),
        (3, [true, true, true, false]),   // Case 3: |S_rim| = 3 (corner/gap tet)
        (4, [true, true, true, true]),    // |S_rim| = 4 (all vertices on rim -> 0 shifted edges)
    ];

    for &xi_qp in &[0.25, 0.30, 0.375, 0.50] {
        for &(n_rim, mask) in cases {
            let ct = make_qp_tet(v, mask, xi_qp);
            for l in &q8.lam {
                let lam = [l[0], l[1], l[2], l[3]];
                let lam_out: f64 = (0..4).filter(|&i| !mask[i]).map(|i| lam[i]).sum();
                let (_, det) = ct.jacobian(&lam);
                let ratio_num = det / det_s;
                let ratio_ana = analytical_qp_jacobian_ratio(n_rim, lam_out, xi_qp);
                assert!(
                    (ratio_num - ratio_ana).abs() < 1e-12,
                    "n_rim={n_rim}, xi_qp={xi_qp}, lam_out={lam_out}: numerical={ratio_num} vs analytical={ratio_ana}"
                );
                assert!(ratio_num > 0.0, "interior Gauss point must have strictly positive Jacobian ratio");

                if (xi_qp - 0.25).abs() < 1e-15 {
                    let closed_form_025 = match n_rim {
                        1 => 2.0 * lam_out.powi(3),
                        2 => 2.0 * lam_out.powi(2) * (1.0 + lam_out),
                        3 => 2.0 * lam_out * (1.0 + lam_out).powi(2),
                        _ => 1.0,
                    };
                    assert!(
                        (ratio_num - closed_form_025).abs() < 1e-12,
                        "xi_qp=0.25 closed form mismatch for n_rim={n_rim}: {ratio_num} vs {closed_form_025}"
                    );
                }
            }
        }
    }

    // At xi_qp = 0.25, verify det(J) == 0 on the rim (lam_out == 0) for Cases 1, 2, 3
    let ct_case1 = make_qp_tet(v, [true, true, false, false], 0.25);
    assert!(ct_case1.jacobian(&[0.5, 0.5, 0.0, 0.0]).1.abs() < 1e-14 * det_s);
    let ct_case2 = make_qp_tet(v, [true, false, false, false], 0.25);
    assert!(ct_case2.jacobian(&[1.0, 0.0, 0.0, 0.0]).1.abs() < 1e-14 * det_s);
    let ct_case3 = make_qp_tet(v, [true, true, true, false], 0.25);
    assert!(ct_case3.jacobian(&[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0]).1.abs() < 1e-14 * det_s);
}

fn build_sheet_pad_model() -> Model {
    let mut m = Model::new_3d();
    let si = m.box_("si", (-120e-6, -120e-6, -80e-6), (120e-6, 120e-6, 0.0));
    let vac = m.box_("vac", (-120e-6, -120e-6, 0.0), (120e-6, 120e-6, 80e-6));

    // Ground plane with [-50..50] x [-50..50] um pocket cutout (270 deg re-entrant corners)
    let gnd_outer = vec![pt2(-120e-6, -120e-6), pt2(120e-6, -120e-6), pt2(120e-6, 120e-6), pt2(-120e-6, 120e-6)];
    let gnd_hole = vec![pt2(-50e-6, -50e-6), pt2(50e-6, -50e-6), pt2(50e-6, 50e-6), pt2(-50e-6, 50e-6)];
    let p_gnd = m.polygons("p_gnd", vec![gnd_outer], vec![gnd_hole]);

    // Inner pad [-30..30] x [-30..30] um (90 deg convex corners)
    let pad_outer = vec![pt2(-30e-6, -30e-6), pt2(30e-6, -30e-6), pt2(30e-6, 30e-6), pt2(-30e-6, 30e-6)];
    let p_pad = m.polygons("p_pad", vec![pad_outer], vec![]);

    m.sheet("gnd", p_gnd, 0.0);
    m.sheet("pad", p_pad, 0.0);
    m.overlay(&[("Si", si), ("VacTop", vac)]);
    m
}

#[test]
fn apply_quarter_point_nodes_enforces_216_pt_q8_positivity_and_curved_summary() {
    let m = build_sheet_pad_model();
    let mut g = mesh_model(&m.list, &MeshPlan::simple(35e-6)).unwrap();

    let mut sheet_labels = g.boundary("pad");
    sheet_labels.extend(g.boundary("gnd"));
    let rim_verts = sheet_rim_vertices(&g.mesh, &sheet_labels);
    let n_rim = rim_verts.iter().filter(|&&r| r).count();
    assert!(n_rim >= 8, "expected rim vertices on pad and ground pocket, got {n_rim}");

    let rep = apply_quarter_point_nodes(&mut g.mesh, &rim_verts, 0.25, &BTreeSet::new());
    assert_eq!(rep.rim_vertices, n_rim);
    assert!(rep.shifted_edges > 0, "expected quarter-point shifted edges");
    assert!(rep.qp_cells > 0, "expected quarter-point tetrahedra");
    assert_eq!(rep.unshifted_edges, 0, "straight-rim prismatic sheet mesh should have 0 unshifted QP edges");
    assert!(
        rep.min_interior_jacobian_ratio > 0.0,
        "min_interior_jacobian_ratio must be strictly positive, got {:e}",
        rep.min_interior_jacobian_ratio
    );
    assert!(g.mesh.has_curved() && g.mesh.has_qp());

    // Verify no rim-to-rim edge was shifted
    for (ei, &(a, b)) in g.mesh.edges().iter().enumerate() {
        let both_rim = rim_verts[a as usize] && rim_verts[b as usize];
        let neither_rim = !rim_verts[a as usize] && !rim_verts[b as usize];
        if both_rim || neither_rim {
            assert!(
                g.mesh.edge_node(ei).is_none(),
                "edge ({a}, {b}) with rim_verts=({}, {}) must remain straight",
                rim_verts[a as usize],
                rim_verts[b as usize]
            );
        } else {
            assert!(
                g.mesh.edge_is_qp(ei) && g.mesh.edge_node(ei).is_some(),
                "edge ({a}, {b}) with exactly one rim endpoint must be QP-shifted"
            );
        }
    }

    // Verify every curved/QP cell passes 216-pt q8 interior positivity, while vertex-based scaled_jacobian_min hits 0
    let mut min_vertex_sj = f64::INFINITY;
    for c in 0..g.mesh.n_cells() {
        if g.mesh.cell_is_curved(c) {
            assert!(g.mesh.cell_is_qp(c));
            let q8_min = scaled_jacobian_min_interior_q8(&g.mesh, c);
            assert!(q8_min > 0.0, "cell {c} failed 216-pt q8 positivity: {q8_min:e}");
            min_vertex_sj = min_vertex_sj.min(scaled_jacobian_min(&g.mesh, c));
        }
    }
    assert!(
        min_vertex_sj.abs() < 1e-12,
        "exact xi_QP = 0.25 has det(J) = 0 at rim vertices (got {min_vertex_sj:e})"
    );

    // CurvedSummary::from_mesh must validate QP cells via q8 interior points (> 0), not vertex det(J) = 0
    let summary = CurvedSummary::from_mesh(&g.mesh).expect("curved summary present");
    assert_eq!(summary.edges, rep.shifted_edges);
    assert_eq!(summary.cells, rep.qp_cells);
    assert!(
        (summary.min_scaled_jacobian - rep.min_interior_jacobian_ratio).abs() < 1e-15
            && summary.min_scaled_jacobian > 0.0
    );
}

#[test]
fn curved_rim_sagitta_guard_clamps_xi_qp_upward_and_prevents_inversion() {
    // Two tetrahedra sharing face (0, 1, 2) with rim edge (0, 1) along x in [0, 1]
    let points = vec![
        0.0, 0.0, 0.0, // v0 (rim)
        1.0, 0.0, 0.0, // v1 (rim)
        0.5, 1.0, 0.0, // v2 (interior sheet/dielectric)
        0.5, 0.4, 1.0, // v3 (upper tet apex)
        0.5, 0.4, -1.0, // v4 (lower tet apex)
    ];
    let cells = vec![0, 1, 2, 3, 0, 2, 1, 4];
    let cell_labels = vec![1, 2];
    let facets = vec![0, 1, 2];
    let facet_labels = vec![1];
    let mut mesh = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);

    // Give rim edge [0, 1] a positive inward sagitta delta = 0.08 in +y (convex toward v2):
    // without sagitta clamping, xi_QP = 0.25 gives (4*xi - 1)*((4*xi - 1) - 2*delta/h_y) <= 0 at the rim!
    let delta = 0.08;
    mesh.set_edge_nodes(&[([0, 1], [0.5, delta, 0.0])]);

    let rim_verts = vec![true, true, false, false, false];
    let mut curved_rim_edges = BTreeSet::new();
    curved_rim_edges.insert([0u32, 1u32]);

    let rep = apply_quarter_point_nodes(&mut mesh, &rim_verts, 0.25, &curved_rim_edges);
    assert_eq!(rep.rim_vertices, 2);
    assert!(rep.shifted_edges > 0, "transverse edges must be shifted with clamped xi_QP");
    assert_eq!(rep.clamped_curved_edges, rep.shifted_edges, "all edges touching curved rim [0, 1] must clamp xi_QP");
    assert!(
        rep.min_interior_jacobian_ratio > 0.0,
        "clamped curved-rim QP cells must stay strictly positive on all 216 q8 points (got {:e})",
        rep.min_interior_jacobian_ratio
    );

    // Check that the shifted position along edge (0, 2) satisfies 4 * xi_eff - 1 > 2 * delta / h_y + 0.05
    let e02 = mesh.edge_index(0, 2);
    let q02 = mesh.edge_node(e02).expect("edge (0, 2) must carry clamped QP node");
    let p0 = [0.0_f64, 0.0, 0.0];
    let p2 = [0.5_f64, 1.0, 0.0];
    let h02 = ((p2[0] - p0[0]).powi(2) + (p2[1] - p0[1]).powi(2)).sqrt();
    let dist_from_rim = ((q02[0] - p0[0]).powi(2) + (q02[1] - p0[1]).powi(2)).sqrt();
    let xi_eff = dist_from_rim / h02;
    assert!(
        4.0 * xi_eff - 1.0 > 2.0 * (delta / h02) + 0.05,
        "sagitta guard violated: xi_eff={xi_eff}, 4*xi_eff-1={}, rhs={}",
        4.0 * xi_eff - 1.0,
        2.0 * (delta / h02) + 0.05
    );

    // Also verify relax_edge_nodes extracted in curve3.rs relaxes overly-curved nodes cleanly
    let mut raw_nodes = vec![([0u32, 1u32], [0.5, 0.45, 0.0])];
    let tet_cells = [[0u32, 1, 2, 3], [0, 2, 1, 4]];
    let pos_fn = |v: u32| -> [f64; 3] {
        let p = mesh.point(v as usize);
        [p[0], p[1], p[2]]
    };
    let c_rep = relax_edge_nodes(&pos_fn, &tet_cells, &mut raw_nodes);
    assert!(c_rep.relaxed > 0 || c_rep.straightened > 0);
    assert!(c_rep.min_scaled_jacobian > 0.1 || raw_nodes.is_empty());
}

#[test]
fn gradients_at_nodes_at_exact_quarter_point_xi_025_produces_zero_nan_or_inf() {
    let m = build_sheet_pad_model();
    let mut g = mesh_model(&m.list, &MeshPlan::simple(35e-6)).unwrap();
    let mut sheet_labels = g.boundary("pad");
    sheet_labels.extend(g.boundary("gnd"));
    let rim_verts = sheet_rim_vertices(&g.mesh, &sheet_labels);
    let rep = apply_quarter_point_nodes(&mut g.mesh, &rim_verts, 0.25, &BTreeSet::new());
    assert!(rep.shifted_edges > 0);

    // Populate a P2 field with r^(1/2) cylindrical edge singularity + linear gradient around x = 30 um, z = 0
    let order = 2;
    let coords = g.mesh.dof_coordinates(order);
    let ndofs = g.mesh.n_dofs(order);
    let mut u = Vec::with_capacity(ndofs);
    for p in coords.chunks(3) {
        let dx = p[0] - 30e-6;
        let dz = p[2];
        let r = dx.hypot(dz);
        let theta = dz.atan2(dx);
        let v_sing = r.sqrt() * (0.5 * theta).sin() + 1.2e3 * p[0] - 0.8e3 * p[1] + 2.5e3 * p[2];
        u.push(v_sing);
    }

    let g_cell = gradients_at_nodes(&g.mesh, &u, 1, order);
    assert_eq!(g_cell.len(), g.mesh.n_cells() * 10 * 3);
    assert!(
        g_cell.iter().all(|x| x.is_finite()),
        "gradients_at_nodes produced non-finite values on xi_QP = 0.25 mesh!"
    );

    let g_nodal = node_average(&g.mesh, &g_cell, 3, order);
    assert_eq!(g_nodal.len(), ndofs * 3);
    assert!(
        g_nodal.iter().all(|x| x.is_finite()),
        "node_average of gradients_at_nodes produced non-finite values!"
    );
    let max_grad = g_nodal.iter().map(|x| x.abs()).fold(0.0f64, f64::max);
    assert!(max_grad > 1e3, "expected non-trivial electric field magnitude, got {max_grad:e}");
}

#[test]
fn morrison_lewis_corner_exponents_and_babuska_grading_floor() {
    // 1. Check Morrison-Lewis corner exponent table & interpolation
    let a_90 = morrison_lewis_exponent(0.5 * PI, 11.45);
    let a_180 = morrison_lewis_exponent(PI, 11.45);
    let a_270 = morrison_lewis_exponent(1.5 * PI, 11.45);
    assert!((a_90 - 0.2966).abs() < 1e-3, "90 deg convex corner alpha = {a_90}");
    assert!((a_180 - 0.5000).abs() < 1e-12, "180 deg straight edge alpha = {a_180}");
    assert!((a_270 - 0.8147).abs() < 1e-3, "270 deg re-entrant corner alpha = {a_270}");

    // Strict monotonicity across (0, 2*pi)
    let mut prev = 0.0;
    for k in 1..36 {
        let omega = (k as f64 / 36.0) * (2.0 * PI);
        let alpha = morrison_lewis_exponent(omega, 1.0);
        assert!(alpha > prev, "alpha({omega}) = {alpha} not > {prev}");
        prev = alpha;
    }

    // 2. Check Babuska radical corner grading law & mandatory 0.05 * h_edge floor
    let h_edge = 10e-6;
    let r_c = CORNER_RADIUS_FACTOR * h_edge;
    let h_min = CORNER_H_MIN_FACTOR * h_edge;
    assert_eq!(corner_grading_h(0.0, h_edge, a_90), h_min);
    assert_eq!(corner_grading_h(1e-12 * h_edge, h_edge, a_90), h_min);
    assert!((corner_grading_h(r_c, h_edge, a_90) - h_edge).abs() < 1e-15 * h_edge);
    assert!((corner_grading_h(3.0 * h_edge, h_edge, a_90) - h_edge).abs() < 1e-15 * h_edge);

    // At d = 0.5 * R_c = h_edge, h(d) = h_edge * 0.5^(1 - alpha / 2.5)
    let d_mid = h_edge;
    let expected_mid = h_edge * 0.5f64.powf(1.0 - a_90 / 2.5);
    assert!((corner_grading_h(d_mid, h_edge, a_90) - expected_mid).abs() < 1e-15 * h_edge);

    // 3. Measure overlay corner angles from Complex AND from 2D triangulation m2
    let m = build_sheet_pad_model();
    let ev = evaluate(&m.list).unwrap();
    let cx = ev.complex();
    let st = cx.stack.as_ref().unwrap();
    let ov = &st.overlay;
    let ov_edge_h = vec![h_edge; ov.edges.len()];
    let corners = collect_overlay_sheet_corners(cx, 0.0, &ov_edge_h, 11.45);
    // Should detect 4 convex 90 deg corners of `pad`, 4 re-entrant 270 deg pocket corners of `gnd`, and 4 outer 90 deg corners of `gnd`
    let n_convex_90 = corners.iter().filter(|c| (c.omega_rad - 0.5 * PI).abs() < 1e-6).count();
    let n_reentrant_270 = corners.iter().filter(|c| (c.omega_rad - 1.5 * PI).abs() < 1e-6).count();
    assert_eq!(n_convex_90, 8, "expected 8 convex 90 deg sheet corners (4 pad + 4 outer gnd)");
    assert_eq!(n_reentrant_270, 4, "expected 4 re-entrant 270 deg pocket corners");

    // Verify coupled_corner_dz and corner_graded_h_at
    let dz = coupled_corner_dz(h_edge, &corners);
    assert!((dz - h_min).abs() < 1e-15 * h_edge, "coupled_corner_dz must clamp to h_min = 0.05 * h_edge");
    let h_at_pad_corner = corner_graded_h_at([30e-6, 30e-6], &corners, h_edge);
    assert!((h_at_pad_corner - h_min).abs() < 1e-15 * h_edge);

    // Also verify 2D triangle-angle sum measurement on a 2D mesh of an L-shaped sheet (five 90 deg corners + one 270 deg corner)
    let pts2 = vec![
        0.0, 0.0, // 0: 90 deg
        2.0, 0.0, // 1: 90 deg
        2.0, 1.0, // 2: 90 deg
        1.0, 1.0, // 3: 270 deg re-entrant
        1.0, 2.0, // 4: 90 deg
        0.0, 2.0, // 5: 90 deg
        2.0, 2.0, // 6: non-sheet corner
    ];
    let tris2 = vec![
        0, 1, 3, // sheet cell 1
        1, 2, 3, // sheet cell 1
        0, 3, 5, // sheet cell 1
        3, 4, 5, // sheet cell 1
        3, 2, 6, // non-sheet cell 2
        3, 6, 4, // non-sheet cell 2
    ];
    let labels2 = vec![1, 1, 1, 1, 2, 2];
    let m2 = SimplexMesh::new(2, pts2, tris2, labels2, vec![0, 1], vec![1]);
    let m2_angles = measure_sheet_corner_angles_m2(&m2, &[true, false]);
    for &v90 in &[0usize, 1, 2, 4, 5] {
        assert!((m2_angles[v90].unwrap() - 0.5 * PI).abs() < 1e-12, "vertex {v90} angle {:?}", m2_angles[v90]);
    }
    assert!((m2_angles[3].unwrap() - 1.5 * PI).abs() < 1e-12, "vertex 3 re-entrant angle {:?}", m2_angles[3]);
    assert!(m2_angles[6].is_none(), "non-sheet vertex 6 must be None");
    assert_eq!(measure_overlay_sheet_corners(cx, 0.0).iter().flatten().count(), 12);

    // 4. Verify MeshRequest & request_schema() round-trip with AF-2 fields
    let req = MeshRequest {
        h: 50e-6,
        max_dev: None,
        mesher: MesherRequest::Simple,
        structured_pml: false,
        waves: None,
        curved: false,
        gap_cells: Some(2.0),
        z_coarsen: true,
        singular_edges: Some("quarter_point".into()),
        quarter_point_xi: Some(0.25),
        corner_grading: true,
    };
    let val = serde_json::to_value(&req).unwrap();
    assert_eq!(val["singular_edges"], "quarter_point");
    assert_eq!(val["quarter_point_xi"], 0.25);
    assert_eq!(val["corner_grading"], true);
    let back: MeshRequest = serde_json::from_value(val).unwrap();
    assert_eq!(back, req);
    let plan = req.plan();
    assert_eq!(plan.opts.singular_edges.as_deref(), Some("quarter_point"));
    assert_eq!(plan.opts.quarter_point_xi, Some(0.25));
    assert!(plan.opts.corner_grading);
    let schema = request_schema();
    assert!(schema["properties"]["singular_edges"].is_object());
    assert!(schema["properties"]["quarter_point_xi"].is_object());
    assert!(schema["properties"]["corner_grading"].is_object());
    let (_, report) = mesh_request(&m.list, &req).unwrap();
    assert_eq!(report.dim, 3);
}
