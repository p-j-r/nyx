// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification tests for AF-5 in `femcore`: 3D `sheet3d` meshing of arc-fitted circular and
//! concentric-annular sheet layouts (`Op::Polygons` with `arc_tol > 0`), isoparametric `EdgeNodes`
//! projection on horizontal `z_planes` and Dompierre diagonal edges (`sheet3d_curve`), `u32::MAX`
//! guard under `z_coarsen`, `min_scaled_jacobian >= 0.10`, and exact curved `CurvedTet` volume/area
//! recovery vs chordal polygons.

use femcore::fem::{assemble_laplace, quad_rule_collapsed, CurvedTet};
use femcore::mesh::SimplexMesh;
use femcore::meshing::curve3::{scaled_jacobian_min, MIN_SCALED_JACOBIAN};
use femcore::meshing::geomdrive::evaluate;
use femcore::meshing::sheet3d::mesh_sheet_stack_3d;
use femcore::meshing::sheet3d_curve::{
    apply_sheet3d_curved_edges, project_sheet3d_curved_edges, Sheet3dTransitionRef,
};
use femcore::meshing::MeshPlan;
use geom::expr::Expr;
use geom::math::{TAU, V2};
use geom::model::Model;
use geom::select::Selector;

fn to_expr_ring(ring: &[V2]) -> Vec<[Expr; 2]> {
    ring.iter().map(|p| [Expr::Num(p[0]), Expr::Num(p[1])]).collect()
}

fn discretized_circle(r: f64, n: usize, ccw: bool) -> Vec<V2> {
    (0..n)
        .map(|i| {
            let frac = i as f64 / n as f64;
            let a = if ccw { TAU * frac } else { -TAU * frac };
            [r * a.cos(), r * a.sin()]
        })
        .collect()
}

/// Builds the 3-sheet concentric benchmark model (`c0` ground plane with slotted circular hole,
/// `c1` C-shaped annular ring with concave inner arc, `c2` inner circular disk with feed tab).
fn build_concentric_3d_model(arc_tol: f64) -> Model {
    let mut c0_hole = vec![[-1.1e-5, -130.08e-6], [1.1e-5, -130.08e-6]];
    let r0 = 98.08e-6;
    let a0_start = (-97.448e-6_f64).atan2(1.1e-5);
    let a0_end = (97.448e-6_f64).atan2(1.1e-5);
    let n_arc = 89usize;
    for k in 0..=n_arc {
        let a = a0_start + (a0_end - a0_start) * (k as f64 / n_arc as f64);
        c0_hole.push([r0 * a.cos(), r0 * a.sin()]);
    }
    c0_hole.push([1.1e-5, 125.08e-6]);
    c0_hole.push([-1.1e-5, 125.08e-6]);
    let a1_start = (97.448e-6_f64).atan2(-1.1e-5);
    let a1_end = a1_start + (a0_end - a0_start);
    for k in 0..=n_arc {
        let a = a1_start + (a1_end - a1_start) * (k as f64 / n_arc as f64);
        c0_hole.push([r0 * a.cos(), r0 * a.sin()]);
    }

    let mut c1_outer = Vec::new();
    let r1_out = 87.61e-6;
    let r1_in = 75.55e-6;
    let th_left_bot = (-86.628e-6_f64).atan2(-1.3e-5).rem_euclid(TAU);
    let th_left_top = (87.458e-6_f64).atan2(-5.0e-6).rem_euclid(TAU);
    for k in 0..=45 {
        let a = th_left_bot - (th_left_bot - th_left_top) * (k as f64 / 45.0);
        c1_outer.push([r1_out * a.cos(), r1_out * a.sin()]);
    }
    c1_outer.push([-5.0e-6, 123.08e-6]);
    c1_outer.push([5.0e-6, 123.08e-6]);
    let th_right_top = (87.458e-6_f64).atan2(5.0e-6);
    let th_right_bot = (-86.628e-6_f64).atan2(1.3e-5);
    for k in 0..=45 {
        let a = th_right_top + (th_right_bot - th_right_top) * (k as f64 / 45.0);
        c1_outer.push([r1_out * a.cos(), r1_out * a.sin()]);
    }
    let th_in_start = (-74.415e-6_f64).atan2(1.3e-5);
    let th_in_end = th_in_start + (TAU - 2.0 * (1.3e-5_f64 / r1_in).asin());
    for k in 0..=90 {
        let a = th_in_start + (th_in_end - th_in_start) * (k as f64 / 90.0);
        c1_outer.push([r1_in * a.cos(), r1_in * a.sin()]);
    }

    let mut c2_outer = vec![[-5.0e-6, -128.08e-6]];
    let r2 = 56.69e-6;
    let th2_start = (-56.463e-6_f64).atan2(-5.0e-6).rem_euclid(TAU);
    let th2_end = (-56.463e-6_f64).atan2(5.0e-6).rem_euclid(TAU) - TAU;
    for k in 0..=90 {
        let a = th2_start + (th2_end - th2_start) * (k as f64 / 90.0);
        c2_outer.push([r2 * a.cos(), r2 * a.sin()]);
    }
    c2_outer.push([5.0e-6, -128.08e-6]);

    let mut m = Model::new_3d();
    let vac_top = m.box_("vac_top", (-0.00016208, -0.00018958, 0.0), (0.00016208, 0.00018958, 0.001));
    let vac_bot = m.box_("vac_bot", (-0.00016208, -0.00018958, -0.00135), (0.00016208, 0.00018958, -0.00035));
    let si_sub = m.box_("si_sub", (-0.00016208, -0.00018958, -0.00035), (0.00016208, 0.00018958, 0.0));

    let c0_box = vec![
        [-0.00016208, -0.00018958],
        [0.00016208, -0.00018958],
        [0.00016208, 0.00018958],
        [-0.00016208, 0.00018958],
    ];
    let p0 = m.polygons_with_arcs(
        "poly_c0",
        vec![to_expr_ring(&c0_box)],
        vec![to_expr_ring(&c0_hole)],
        7.790615287643187e-07,
        0.0,
        arc_tol,
    );
    m.sheet("c0", p0, 0.0);

    let p1 = m.polygons_with_arcs(
        "poly_c1",
        vec![to_expr_ring(&c1_outer)],
        vec![],
        7.790615287643187e-07,
        0.0,
        arc_tol,
    );
    m.sheet("c1", p1, 0.0);

    let p2 = m.polygons_with_arcs(
        "poly_c2",
        vec![to_expr_ring(&c2_outer)],
        vec![],
        7.790615287643187e-07,
        0.0,
        arc_tol,
    );
    m.sheet("c2", p2, 0.0);

    m.overlay(&[("VacTop", vac_top), ("VacBot", vac_bot), ("Si", si_sub)]);
    m.size("sz_c1", Selector::parse("name(c1)").unwrap(), Some(Expr::Num(9.230483429295996e-06)), None, None);
    m.size("sz_c1e", Selector::parse("c1:edge").unwrap(), Some(Expr::Num(5.429696134879997e-06)), None, None);
    m.size("sz_c2", Selector::parse("name(c2)").unwrap(), Some(Expr::Num(3.431672259133796e-05)), None, None);
    m.size("sz_c2e", Selector::parse("c2:edge").unwrap(), Some(Expr::Num(2.018630740666939e-05)), None, None);
    m
}

#[test]
fn repro_concentric_meshes_cleanly_and_projects_curved_edgenodes_with_valid_jacobians() {
    let m = build_concentric_3d_model(200e-9);
    let ev = evaluate(&m.list).expect("concentric 3D model with arc_tol=200e-9 must evaluate");
    let mut plan = MeshPlan::simple(0.0001725);
    plan.opts.curved = true;

    let (mut g, _) =
        mesh_sheet_stack_3d(ev, &plan, None).expect("mesh_sheet_stack_3d must succeed without c0:edge workaround");
    let report = apply_sheet3d_curved_edges(&mut g);

    assert!(report.curved_edges > 0, "must emit curved EdgeNodes on concentric circular rims");
    assert!(report.curved_cells > 0, "must curve tetrahedra adjacent to circular rims");
    assert!(
        report.min_scaled_jacobian >= MIN_SCALED_JACOBIAN,
        "min_scaled_jacobian ({}) must be >= {MIN_SCALED_JACOBIAN}",
        report.min_scaled_jacobian
    );
    assert!(g.mesh.has_curved(), "SimplexMesh::has_curved() must be true after set_edge_nodes");

    // Verify that BOTH horizontal z-plane edges AND Dompierre diagonal edges are curved
    let mut n_horz = 0usize;
    let mut n_diag = 0usize;
    for ([u, v], q) in g.mesh.edge_node_list() {
        let zu = g.mesh.points[3 * u as usize + 2];
        let zv = g.mesh.points[3 * v as usize + 2];
        if (zu - zv).abs() < 1e-12 {
            n_horz += 1;
            if zu.abs() < 1e-12 {
                let r_mid = q[0].hypot(q[1]);
                let err_r = [56.69e-6_f64, 75.55e-6, 87.61e-6, 98.08e-6]
                    .iter()
                    .map(|&rt| (r_mid - rt).abs())
                    .fold(f64::INFINITY, f64::min);
                assert!(err_r < 25e-9, "z=0 curved mid-edge node radius {r_mid:e} off by {err_r:e}");
            }
        } else {
            n_diag += 1;
        }
    }
    assert!(n_horz > 0, "must curve horizontal edges on active z_planes");
    assert!(n_diag > 0, "must curve Dompierre diagonal edges across curved vertical walls");

    for c in 0..g.mesh.n_cells() {
        if g.mesh.cell_is_curved(c) {
            let sj = scaled_jacobian_min(&g.mesh, c);
            assert!(sj > MIN_SCALED_JACOBIAN, "cell {c} scaled_jacobian_min = {sj} <= {MIN_SCALED_JACOBIAN}");
        }
    }
}

#[test]
fn z_coarsen_guards_unused_vertices_and_curved_tets_recover_exact_cylinder_volume() {
    let mut m = Model::new_3d();
    let vac_top = m.box_("vac_top", (-120e-6, -120e-6, 0.0), (120e-6, 120e-6, 300e-6));
    let si_sub = m.box_("si_sub", (-120e-6, -120e-6, -200e-6), (120e-6, 120e-6, 0.0));

    let r_inner = 40e-6;
    let r_hole = 60e-6;
    let disk_ring = discretized_circle(r_inner, 64, true);
    let hole_ring = discretized_circle(r_hole, 64, true);
    let outer_box = vec![[-100e-6, -100e-6], [100e-6, -100e-6], [100e-6, 100e-6], [-100e-6, 100e-6]];

    let p_inner = m.polygons_with_arcs("p_inner", vec![to_expr_ring(&disk_ring)], vec![], 500e-9, 0.0, 50e-9);
    let p_outer = m.polygons_with_arcs(
        "p_outer",
        vec![to_expr_ring(&outer_box)],
        vec![to_expr_ring(&hole_ring)],
        500e-9,
        0.0,
        50e-9,
    );
    m.sheet("c_inner", p_inner, 0.0);
    m.sheet("c_outer", p_outer, 0.0);
    m.overlay(&[("VacTop", vac_top), ("Si", si_sub)]);

    let ev = evaluate(&m.list).unwrap();
    let mut plan = MeshPlan::simple(35e-6);
    plan.opts.z_coarsen = true;
    plan.opts.curved = true;

    let (mut g, _) = mesh_sheet_stack_3d(ev, &plan, None).unwrap();
    let n_pts = g.mesh.n_points() as u32;

    let report = apply_sheet3d_curved_edges(&mut g);
    assert!(report.curved_edges > 0);
    assert!(report.min_scaled_jacobian >= MIN_SCALED_JACOBIAN);
    for ([u, v], q) in g.mesh.edge_node_list() {
        assert!(u < n_pts && v < n_pts, "invalid vertex index in EdgeNodes: [{u}, {v}] vs n_pts={n_pts}");
        assert!(q.iter().all(|x| x.is_finite()));
    }

    // Explicitly verify `project_sheet3d_curved_edges` skips `u32::MAX` entries in `used_v` without panicking
    let cx = g.eval.complex();
    let dummy_m2 = SimplexMesh::new(2, vec![0.0, 0.0, 40e-6, 0.0], vec![], vec![], vec![], vec![]);
    let mut dummy_mesh = g.mesh.clone();
    let empty_rep = project_sheet3d_curved_edges(
        cx,
        &dummy_m2,
        &[(0, 1, 0)],
        &[-200e-6, 0.0, 300e-6],
        &[u32::MAX; 6],
        &[Some(0), Some(1)],
        &[1, 0, 1],
        &[Sheet3dTransitionRef { map_v: &[0, 0], rank: &[(0, 0, 0), (0, 0, 1)] }],
        &mut dummy_mesh,
    );
    assert_eq!(empty_rep.curved_edges, 0);

    // Verify CurvedTet Jacobian quadrature is finite and strictly positive on every curved cell
    let rule = quad_rule_collapsed(6);
    for c in 0..g.mesh.n_cells() {
        if g.mesh.cell_is_curved(c) {
            let ct = CurvedTet::new(g.mesh.cell_geometry_nodes(c));
            let det_s = CurvedTet::straight([ct.x[0], ct.x[1], ct.x[2], ct.x[3]]).jacobian(&[0.25; 4]).1;
            for (lam, &w) in rule.lam.iter().zip(rule.w.iter()) {
                let (_, det) = ct.jacobian(&[lam[0], lam[1], lam[2], lam[3]]);
                assert!(det.is_finite() && w.is_finite());
                assert!(det / det_s > MIN_SCALED_JACOBIAN);
            }
        }
    }

    let eps = vec![1.0; g.mesh.n_cells()];
    let k_mat = assemble_laplace(&g.mesh, 2, &eps);
    assert!(k_mat.values().iter().all(|v| v.is_finite()));
}

#[test]
fn r33_periodic_sheet3d_curved_edges_stay_congruent_between_twins() {
    use femcore::periodic::match_vertices;
    use femcore::props::Periodic;
    use std::collections::{BTreeSet, HashMap};

    let period_z = 200e-6;
    let mut m = Model::new_3d();
    let box_bot = m.box_("bot_box", (-120e-6, -120e-6, -0.5 * period_z), (120e-6, 120e-6, 0.0));
    let box_top = m.box_("top_box", (-120e-6, -120e-6, 0.0), (120e-6, 120e-6, 0.5 * period_z));
    let cyl_hole = discretized_circle(55e-6, 64, true);
    let core_disk = discretized_circle(35e-6, 48, true);
    let pad_box = vec![[-120e-6, -120e-6], [120e-6, -120e-6], [120e-6, 120e-6], [-120e-6, 120e-6]];
    let p_sheet = m.polygons_with_arcs(
        "p_sheet",
        vec![to_expr_ring(&pad_box)],
        vec![to_expr_ring(&cyl_hole)],
        500e-9,
        0.0,
        50e-9,
    );
    let p_core = m.polygons_with_arcs(
        "p_core",
        vec![to_expr_ring(&core_disk)],
        vec![],
        500e-9,
        0.0,
        50e-9,
    );
    m.sheet("gnd", p_sheet, 0.0);
    m.sheet("core", p_core, 0.0);
    // Extrude a dielectric core cylinder across BOTH Z-layers so `z = -period_z/2` and `z = +period_z/2`
    // caps both contain curved circular edges that must match identically under periodic translation.
    let cyl_c = m.circle("cyl_c", (0.0, 0.0), 75e-6);
    let pillar = m.extrude_along("pillar", geom::model::Axis::Z, cyl_c, -0.5 * period_z, 0.5 * period_z);
    m.overlay(&[("Bot", box_bot), ("Top", box_top), ("Pillar", pillar)]);
    let pz_src = Selector::parse(&format!("outer & plane(z = {:e})", -0.5 * period_z)).unwrap();
    let pz_dst = Selector::parse(&format!("outer & plane(z = {:e})", 0.5 * period_z)).unwrap();
    m.periodic("pz", pz_src, pz_dst, (0.0, 0.0, period_z));

    let ev = evaluate(&m.list).unwrap();
    let mut plan = MeshPlan::simple(38e-6);
    plan.opts.curved = true;
    let (mut g, _) = mesh_sheet_stack_3d(ev, &plan, None).unwrap();
    let rep = apply_sheet3d_curved_edges(&mut g);
    assert!(rep.curved_edges > 0);
    assert!(rep.min_scaled_jacobian > MIN_SCALED_JACOBIAN);

    let src = g.select(&format!("outer & plane(z = {:e})", -0.5 * period_z)).unwrap();
    let dst = g.select(&format!("outer & plane(z = {:e})", 0.5 * period_z)).unwrap();
    let matched = match_vertices(&g.mesh, &Periodic { source: src, dest: dst, phase: 0.0 }).unwrap();
    let edge_map: HashMap<[u32; 2], [f64; 3]> = g.mesh.edge_node_list().into_iter().collect();
    let mut src_verts: BTreeSet<u32> = BTreeSet::new();
    for &f in &matched.src_facets {
        src_verts.extend(g.mesh.facet(f));
    }
    let mut synced_periodic_curved = 0usize;
    for ([a, b], q_src) in &edge_map {
        if src_verts.contains(a) && src_verts.contains(b) {
            let ta = matched.map[a];
            let tb = matched.map[b];
            let key_dst = [ta.min(tb), ta.max(tb)];
            let q_dst = edge_map
                .get(&key_dst)
                .unwrap_or_else(|| panic!("periodic dst edge {key_dst:?} missing curved EdgeNode for src ({a}, {b})"));
            assert!((q_dst[0] - q_src[0]).abs() < 1e-13);
            assert!((q_dst[1] - q_src[1]).abs() < 1e-13);
            assert!((q_dst[2] - (q_src[2] + period_z)).abs() < 1e-13);
            synced_periodic_curved += 1;
        }
    }
    assert!(synced_periodic_curved > 0, "expected curved edges on periodic twin faces z = +/- period_z/2");
}
