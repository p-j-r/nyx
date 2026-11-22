// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification suite for cQED electrostatic solvers (QFR-2, QFR-3 part 2, QFR-4):
//! * QFR-2: Two-level P1-iso-P2 + Thomas vertical-line AMG smoother + simultaneous multi-RHS Block-CG,
//!   matching Direct Cholesky to `< 1e-8` rel and bit-identical (`0.0` drift) across `NYX_THREADS=1` vs `4`.
//! * QFR-3: Barsoum quarter-point P2 mid-edge node enrichment (`singular_edges = "quarter_point"`) on
//!   3D tetrahedra meeting open conductor sheet perimeter edges, with positive `CurvedTet` Jacobian validation.
//! * QFR-4: 2D layered-dielectric boundary-element solver (`formulation = "planar_bem"`, `richardson = true`),
//!   exact layered Green's function limits (`t -> inf`, `t -> 0`, grounded backplane `k23 = -1`), analytical
//!   isolated disk `C = 8 eps0 eps_eff R` (`< 0.5%`), and `< 2%` agreement with 3D `volumetric_fem`.

use femcore::fem::{quad_rule_collapsed, CurvedTet};
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::props::LinearSolver;
use geom::expr::Expr;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_dfs::electro::enrich_quarter_point_mesh;
use nyx_dfs::planar_bem::{layered_greens_function, solve_planar_bem_with_stack, LayeredStack};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Dfs, Problem, Study, EPS0};
use nyx_physics::{CapacitanceMatrixData, Ctx, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, TlsMeta};
use serde_json::json;
use std::collections::BTreeSet;
use std::f64::consts::PI;

fn pt(x: f64, y: f64) -> [Expr; 2] {
    [Expr::Num(x), Expr::Num(y)]
}

/// Builds a 3-terminal planar transmon + coupler 3D sheet-stack model on a layered Si/Vacuum box.
fn build_transmon_3port_model(global_h: f64, pad_edge_h: f64, z_box: f64) -> (MeshedGeom, Problem) {
    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-300e-6, -300e-6, -z_box), (300e-6, 300e-6, 0.0));
    let vac = m.box_("vacuum", (-300e-6, -300e-6, 0.0), (300e-6, 300e-6, z_box));

    let p1 = m.polygons(
        "p1_poly",
        vec![vec![pt(-160e-6, -50e-6), pt(-60e-6, -50e-6), pt(-60e-6, 50e-6), pt(-160e-6, 50e-6)]],
        vec![],
    );
    let p2 = m.polygons(
        "p2_poly",
        vec![vec![pt(-20e-6, -50e-6), pt(80e-6, -50e-6), pt(80e-6, 50e-6), pt(-20e-6, 50e-6)]],
        vec![],
    );
    let c = m.polygons(
        "c_poly",
        vec![vec![pt(120e-6, -50e-6), pt(180e-6, -50e-6), pt(180e-6, 50e-6), pt(120e-6, 50e-6)]],
        vec![],
    );

    m.sheet("pad_1", p1, 0.0);
    m.sheet("pad_2", p2, 0.0);
    m.sheet("coupler", c, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    for s in ["pad_1:edge", "pad_2:edge", "coupler:edge"] {
        m.size(
            &format!("sz_{}", s.replace(':', "_")),
            Selector::parse(s).unwrap(),
            Some(Expr::Num(pad_edge_h)),
            None,
            None,
        );
    }
    for s in ["pad_1", "pad_2", "coupler"] {
        m.size(&format!("sz_face_{s}"), Selector::Name(s.into()), Some(Expr::Num(pad_edge_h * 1.5)), None, None);
    }

    let g = mesh_model(&m.list, &MeshPlan::simple(global_h)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", 11.45);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());

    let p1_labs = g.boundary("pad_1");
    let p2_labs = g.boundary("pad_2");
    let c_labs = g.boundary("coupler");
    let sheet_set: BTreeSet<i32> = p1_labs.iter().chain(&p2_labs).chain(&c_labs).copied().collect();
    let ext = g.mesh.exterior_facets();
    let mut ground_labs = Vec::new();
    for (f, &is_ext) in ext.iter().enumerate() {
        let lab = g.mesh.facet_labels[f];
        if is_ext && !sheet_set.contains(&lab) && !ground_labs.contains(&lab) {
            ground_labs.push(lab);
        }
    }

    prob.bcs = vec![
        Bc { labels: p1_labs, kind: BcKind::Terminal("pad_1".into()) },
        Bc { labels: p2_labs, kind: BcKind::Terminal("pad_2".into()) },
        Bc { labels: c_labs, kind: BcKind::Terminal("coupler".into()) },
        Bc { labels: ground_labs, kind: BcKind::Ground },
    ];
    (g, prob)
}

#[test]
fn test_qfr2_block_cg_line_amg_matches_direct_and_is_bit_deterministic_across_threads() {
    // Subprocess worker mode for multi-process NYX_THREADS=1 vs 4 bit-identity verification:
    if let Ok(mode) = std::env::var("NYX_QFR2_SUBPROC_CHECK") {
        let (g, mut prob) = build_transmon_3port_model(85e-6, 20e-6, 120e-6);
        prob.solver = if mode == "direct" { LinearSolver::Direct } else { LinearSolver::Iterative };
        prob.tol = 1e-11;
        prob.amg_line_smooth = true;
        let res = solve_electrostatics(&g.mesh, &prob).unwrap();
        let cap = CapacitanceMatrixData::of(&res).unwrap();
        let mut bits: Vec<u64> = Vec::new();
        for row in &cap.maxwell_matrix_f {
            for &v in row {
                bits.push(v.to_bits());
            }
        }
        // Also include potential field checksums from all 3 terminal solutions
        for sol in &res.solutions {
            for &v in sol.fields[0].values.iter().step_by(7) {
                bits.push(v.to_bits());
            }
        }
        println!("NYX_QFR2_BITS:{}", serde_json::to_string(&bits).unwrap());
        return;
    }

    // 1. In-process test on both standard and thin-layer (high-z-aspect-ratio) 3D sheet meshes:
    for &z_box in &[120e-6, 25e-6] {
        let (g, mut prob_direct) = build_transmon_3port_model(85e-6, 20e-6, z_box);
        let (p1_cols, p2_cols) = femcore::fem::detect_vertical_columns(&g.mesh, 2);
        assert!(
            !p1_cols.is_empty() && !p2_cols.is_empty(),
            "expected vertical prismatic columns to be auto-detected on 3D sheet mesh"
        );

        prob_direct.solver = LinearSolver::Direct;
        let res_direct = solve_electrostatics(&g.mesh, &prob_direct).unwrap();
        let cap_direct = CapacitanceMatrixData::of(&res_direct).unwrap();

        let mut prob_amg = prob_direct.clone();
        prob_amg.solver = LinearSolver::Iterative;
        prob_amg.tol = 1e-12;
        prob_amg.amg_line_smooth = true;
        let res_amg = solve_electrostatics(&g.mesh, &prob_amg).unwrap();
        let cap_amg = CapacitanceMatrixData::of(&res_amg).unwrap();

        assert!(
            res_amg.log.iter().any(|l| l.contains("vertical-line columns") && l.contains("block iterations")),
            "expected vertical-line smoother and block iterations in AMG log, got {:?}",
            res_amg.log
        );

        let mut max_rel_err = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let cd = cap_direct.maxwell_matrix_f[i][j];
                let ca = cap_amg.maxwell_matrix_f[i][j];
                let rel = (ca - cd).abs() / cd.abs().max(1e-30);
                max_rel_err = max_rel_err.max(rel);
            }
        }
        eprintln!(
            "QFR-2 (z_box={z_box:.1e} m): Block-CG + line-AMG vs Direct Cholesky max_rel_err = {max_rel_err:.3e}"
        );
        assert!(
            max_rel_err < 1e-8,
            "Block-CG + line-AMG vs Direct Cholesky relative error {max_rel_err:e} exceeds 1e-8"
        );
    }

    // 2. Spawn child processes with NYX_THREADS=1 and NYX_THREADS=4 to verify 0.0 bit drift:
    let exe = std::env::current_exe().unwrap();
    let run_subproc = |threads: &str, mode: &str| -> Vec<u64> {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("test_qfr2_block_cg_line_amg_matches_direct_and_is_bit_deterministic_across_threads")
            .arg("--nocapture")
            .env("NYX_THREADS", threads)
            .env("NYX_QFR2_SUBPROC_CHECK", mode)
            .output()
            .expect("failed to spawn child test process");
        assert!(out.status.success(), "child process failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8_lossy(&out.stdout);
        let line = stdout
            .lines()
            .find_map(|l| l.strip_prefix("NYX_QFR2_BITS:"))
            .expect("missing NYX_QFR2_BITS line in child stdout");
        serde_json::from_str(line).unwrap()
    };

    let bits_amg_t1 = run_subproc("1", "amg");
    let bits_amg_t4 = run_subproc("4", "amg");
    assert_eq!(
        bits_amg_t1, bits_amg_t4,
        "Block-CG + line-AMG solutions must be 100% bit-identical (0.0 drift) between NYX_THREADS=1 and NYX_THREADS=4"
    );

    let bits_dir_t1 = run_subproc("1", "direct");
    let bits_dir_t4 = run_subproc("4", "direct");
    assert_eq!(
        bits_dir_t1, bits_dir_t4,
        "Direct Cholesky solutions must be 100% bit-identical (0.0 drift) between NYX_THREADS=1 and NYX_THREADS=4"
    );
}

/// Builds a 3D coplanar strip capacitor on a Si/Vacuum stack for quarter-point singularity testing.
fn build_coplanar_strip_mesh(global_h: f64, edge_h: f64) -> (MeshedGeom, Problem) {
    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-180e-6, -120e-6, -100e-6), (180e-6, 120e-6, 0.0));
    let vac = m.box_("vacuum", (-180e-6, -120e-6, 0.0), (180e-6, 120e-6, 100e-6));

    let s1 = m.polygons(
        "s1_poly",
        vec![vec![pt(-90e-6, -40e-6), pt(-15e-6, -40e-6), pt(-15e-6, 40e-6), pt(-90e-6, 40e-6)]],
        vec![],
    );
    let s2 = m.polygons(
        "s2_poly",
        vec![vec![pt(15e-6, -40e-6), pt(90e-6, -40e-6), pt(90e-6, 40e-6), pt(15e-6, 40e-6)]],
        vec![],
    );
    m.sheet("strip_1", s1, 0.0);
    m.sheet("strip_2", s2, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    for s in ["strip_1:edge", "strip_2:edge"] {
        m.size(
            &format!("sz_{}", s.replace(':', "_")),
            Selector::parse(s).unwrap(),
            Some(Expr::Num(edge_h)),
            None,
            None,
        );
    }
    let g = mesh_model(&m.list, &MeshPlan::simple(global_h)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", 11.45);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());

    let s1_labs = g.boundary("strip_1");
    let s2_labs = g.boundary("strip_2");
    let sheet_set: BTreeSet<i32> = s1_labs.iter().chain(&s2_labs).copied().collect();
    let ext = g.mesh.exterior_facets();
    let mut ground_labs = Vec::new();
    for (f, &is_ext) in ext.iter().enumerate() {
        let lab = g.mesh.facet_labels[f];
        if is_ext && !sheet_set.contains(&lab) && !ground_labs.contains(&lab) {
            ground_labs.push(lab);
        }
    }
    prob.bcs = vec![
        Bc { labels: s1_labs, kind: BcKind::Terminal("strip_1".into()) },
        Bc { labels: s2_labs, kind: BcKind::Ground },
        Bc { labels: ground_labs, kind: BcKind::Ground },
    ];
    (g, prob)
}

#[test]
fn test_qfr3_quarter_point_singular_edge_enrichment_improves_p2_accuracy() {
    let (g_coarse, mut prob_none) = build_coplanar_strip_mesh(55e-6, 18e-6);
    let mut prob_qp = prob_none.clone().with_singular_edges("quarter_point");

    // 1. Verify `enrich_quarter_point_mesh` shifts rim-incident edges and every curved tet has det(J) * det0 > 0:
    let (enriched, n_rim, n_shifted) =
        enrich_quarter_point_mesh(&g_coarse.mesh, &prob_qp).expect("quarter-point enrichment must activate");
    assert!(n_rim > 0 && n_shifted > 0);
    assert!(enriched.has_curved());
    let cquad = quad_rule_collapsed(6);
    let mut curved_cells = 0usize;
    for c in 0..enriched.n_cells() {
        if enriched.cell_is_curved(c) {
            curved_cells += 1;
            let cv = enriched.cell(c);
            let mut v4 = [[0.0; 3]; 4];
            for i in 0..4 {
                let p = enriched.point(cv[i] as usize);
                v4[i] = [p[0], p[1], p[2]];
            }
            let det0 = CurvedTet::straight(v4).jacobian(&[0.25, 0.25, 0.25, 0.25]).1;
            let ct = CurvedTet::new(enriched.cell_geometry_nodes(c));
            for lam in &cquad.lam {
                let (_, det) = ct.jacobian(&[lam[0], lam[1], lam[2], lam[3]]);
                assert!(
                    det * det0 > 0.0,
                    "inverted or singular CurvedTet Jacobian on cell {c}: det={det:e}, det0={det0:e}"
                );
            }
        }
    }
    assert!(curved_cells > 0);

    // 2. Compare coarse standard P2 vs coarse quarter-point P2 against a fine P2 reference mesh:
    prob_none.singular_edges = "none".into();
    let res_none = solve_electrostatics(&g_coarse.mesh, &prob_none).unwrap();
    let c_none = CapacitanceMatrixData::of(&res_none).unwrap().maxwell_matrix_f[0][0];

    prob_qp.singular_edges = "quarter_point".into();
    let res_qp = solve_electrostatics(&g_coarse.mesh, &prob_qp).unwrap();
    assert!(
        res_qp.log.iter().any(|l| l.contains("singular_edges=quarter_point")),
        "expected quarter-point log entry, got {:?}",
        res_qp.log
    );
    let c_qp = CapacitanceMatrixData::of(&res_qp).unwrap().maxwell_matrix_f[0][0];

    let (g_fine, prob_fine) = build_coplanar_strip_mesh(36e-6, 8.0e-6);
    let res_fine = solve_electrostatics(&g_fine.mesh, &prob_fine).unwrap();
    let c_fine = CapacitanceMatrixData::of(&res_fine).unwrap().maxwell_matrix_f[0][0];

    let err_none = (c_none - c_fine).abs() / c_fine;
    let err_qp = (c_qp - c_fine).abs() / c_fine;
    eprintln!(
        "QFR-3 quarter-point enrichment: C_none={c_none:.6e} F (err {:.3}%), C_qp={c_qp:.6e} F (err {:.3}%), C_fine={c_fine:.6e} F",
        err_none * 100.0,
        err_qp * 100.0
    );
    assert!(
        err_qp < err_none,
        "quarter-point enrichment must reduce single-grid error against fine reference (err_qp={err_qp:.4e} vs err_none={err_none:.4e})"
    );
}

#[test]
fn test_qfr4_layered_greens_function_limits_and_isolated_disk_analytic() {
    let rho = 50.0e-6;
    let eps1 = 1.0;
    let eps2 = 11.45;
    let eps3 = 2.1;

    // Limit 1: t_sub -> inf (thick substrate half-space limit)
    let stack_inf =
        LayeredStack { eps_top: eps1, eps_sub: eps2, eps_bot: eps3, t_sub: 1.0e6, grounded_backplane: false };
    let g_inf = layered_greens_function(rho, &stack_inf);
    let g_inf_exact = 1.0 / (4.0 * PI * EPS0 * 0.5 * (eps1 + eps2) * rho);
    let rel_inf = (g_inf - g_inf_exact).abs() / g_inf_exact;
    assert!(rel_inf < 1e-10, "t -> inf Green's function limit failed: rel_err = {rel_inf:e}");

    // Limit 2: t_sub -> 0 (vanishing substrate slab limit -> interface between eps1 and eps3)
    let stack_zero =
        LayeredStack { eps_top: eps1, eps_sub: eps2, eps_bot: eps3, t_sub: 1.0e-13, grounded_backplane: false };
    let g_zero = layered_greens_function(rho, &stack_zero);
    let g_zero_exact = 1.0 / (4.0 * PI * EPS0 * 0.5 * (eps1 + eps3) * rho);
    let rel_zero = (g_zero - g_zero_exact).abs() / g_zero_exact;
    assert!(rel_zero < 1e-10, "t -> 0 Green's function limit failed: rel_err = {rel_zero:e}");

    // Limit 3: Grounded backplane (k23 = -1, eps1 = eps2 = eps_r -> single image at 2*t_sub)
    let t_gp = 40.0e-6;
    let stack_gp = LayeredStack { eps_top: eps2, eps_sub: eps2, eps_bot: 1.0, t_sub: t_gp, grounded_backplane: true };
    let g_gp = layered_greens_function(rho, &stack_gp);
    let g_gp_exact = (1.0 / (4.0 * PI * EPS0 * eps2)) * (1.0 / rho - 1.0 / (rho * rho + (2.0 * t_gp).powi(2)).sqrt());
    let rel_gp = (g_gp - g_gp_exact).abs() / g_gp_exact;
    assert!(rel_gp < 1e-14, "grounded backplane Green's function limit failed: rel_err = {rel_gp:e}");

    // Isolated circular disk of radius R = 100 um on interface between eps1 = 1.0 and eps2 = 11.45:
    // Exact analytical capacitance: C_exact = 8 * eps0 * eps_eff * R
    let r_disk = 100.0e-6;
    let eps_eff = 0.5 * (eps1 + eps2);
    let c_disk_exact = 8.0 * EPS0 * eps_eff * r_disk;

    let mut m2 = Model::new_2d();
    let d_id = m2.circle("disk", (0.0, 0.0), r_disk);
    m2.overlay_full(&[("disk", d_id, RegionKind::Material, None)], &[]);
    let mut opts = GeomMeshOpts::new(14.0e-6);
    opts.max_dev = Some(14.0e-6 / 60.0);
    let g2 = mesh_geom_2d(evaluate(&m2.list).unwrap(), &opts).unwrap();
    let disk_bnds = g2.mesh.boundary_labels();

    let mut prob_disk = Problem::new(Study::Electrostatic, 1).with_formulation("planar_bem").with_richardson(true);
    prob_disk.materials.default.insert("eps_r".into(), eps_eff);
    prob_disk.bcs = vec![Bc { labels: disk_bnds, kind: BcKind::Terminal("disk".into()) }];

    let res_disk = solve_planar_bem_with_stack(
        &g2.mesh,
        &prob_disk,
        Some(LayeredStack {
            eps_top: eps1,
            eps_sub: eps2,
            eps_bot: eps2,
            t_sub: f64::INFINITY,
            grounded_backplane: false,
        }),
    )
    .unwrap();
    let cap_disk = CapacitanceMatrixData::of(&res_disk).unwrap();
    let c_bem = cap_disk.maxwell_matrix_f[0][0];
    let rel_err_disk = (c_bem - c_disk_exact).abs() / c_disk_exact;
    eprintln!(
        "QFR-4 isolated circular disk (R=100 um, eps_eff={eps_eff}): C_bem = {c_bem:.6e} F, C_exact = {c_disk_exact:.6e} F, rel_err = {:.4}%",
        rel_err_disk * 100.0
    );
    assert!(
        rel_err_disk < 0.005,
        "planar_bem isolated disk relative error {:.4}% exceeds 0.5% threshold",
        rel_err_disk * 100.0
    );
}

#[test]
fn test_qfr4_planar_bem_agrees_with_3d_volumetric_fem_on_coplanar_strip_on_si() {
    // 3D coplanar strip on a Silicon substrate (eps_2 = 11.45, thickness t_sub = 40 um) with Vacuum above (eps_1 = 1.0)
    // and a grounded backplane at z = -40 um (`substrate.zmin`).
    let t_sub = 40.0e-6;
    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-240e-6, -240e-6, -t_sub), (240e-6, 240e-6, 0.0));
    let vac = m.box_("vacuum", (-240e-6, -240e-6, 0.0), (240e-6, 240e-6, 160e-6));

    let strip_poly = m.polygons(
        "strip_poly",
        vec![vec![pt(-60e-6, -60e-6), pt(60e-6, -60e-6), pt(60e-6, 60e-6), pt(-60e-6, 60e-6)]],
        vec![],
    );
    m.sheet("strip", strip_poly, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);
    m.size("sz_strip_edge", Selector::parse("strip:edge").unwrap(), Some(Expr::Num(11.0e-6)), None, None);
    m.size("sz_strip_face", Selector::Name("strip".into()), Some(Expr::Num(15.0e-6)), None, None);

    let g = mesh_model(&m.list, &MeshPlan::simple(48.0e-6)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();
    let strip_labs = g.boundary("strip");

    // Identify bottom backplane facets at z = -t_sub
    let mut bot_labs = Vec::new();
    for f in 0..g.mesh.n_facets() {
        let fv = g.mesh.facet(f);
        if fv.iter().all(|&v| (g.mesh.point(v as usize)[2] + t_sub).abs() < 1e-9) {
            let lab = g.mesh.facet_labels[f];
            if !bot_labs.contains(&lab) {
                bot_labs.push(lab);
            }
        }
    }
    assert!(!bot_labs.is_empty(), "expected bottom backplane facets at z = -t_sub");

    let mut prob_fem =
        Problem::new(Study::Electrostatic, 2).with_formulation("volumetric_fem").with_singular_edges("quarter_point");
    prob_fem.materials.set(sub_dom, "eps_r", 11.45);
    prob_fem.materials.set(vac_dom, "eps_r", 1.0);
    prob_fem.domain_names.insert(sub_dom, "substrate".into());
    prob_fem.domain_names.insert(vac_dom, "vacuum".into());
    prob_fem.bcs = vec![
        Bc { labels: strip_labs.clone(), kind: BcKind::Terminal("strip".into()) },
        Bc { labels: bot_labs.clone(), kind: BcKind::Ground },
    ];

    let res_fem = solve_electrostatics(&g.mesh, &prob_fem).unwrap();
    let c_fem = CapacitanceMatrixData::of(&res_fem).unwrap().maxwell_matrix_f[0][0];

    let prob_bem = prob_fem.clone().with_formulation("planar_bem").with_richardson(true);
    let res_bem = solve_electrostatics(&g.mesh, &prob_bem).unwrap();
    let cap_bem = CapacitanceMatrixData::of(&res_bem).unwrap();
    let c_bem = cap_bem.maxwell_matrix_f[0][0];
    assert!(TlsMeta::of(&res_bem).is_some());

    let rel_diff = (c_bem - c_fem).abs() / c_fem;
    eprintln!(
        "QFR-4 strip-on-Si stack (t_sub=40 um, eps_sub=11.45): C_bem = {c_bem:.6e} F, C_fem = {c_fem:.6e} F, rel_diff = {:.3}%",
        rel_diff * 100.0
    );
    assert!(
        rel_diff < 0.02,
        "planar_bem vs 3D volumetric_fem relative difference {:.3}% exceeds 2.0% threshold",
        rel_diff * 100.0
    );

    // Also verify `ProblemSpec` envelope validation and decoding of the new study and physics parameters
    struct SimpleResolver {
        sub_dom: i32,
        vac_dom: i32,
        strip_labs: Vec<i32>,
        bot_labs: Vec<i32>,
    }
    impl LabelResolver for SimpleResolver {
        fn domains(&self, r: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            match r {
                "substrate" => Ok(vec![self.sub_dom]),
                "vacuum" => Ok(vec![self.vac_dom]),
                _ => Err(nyx_physics::Error::invalid(format!("unknown domain {r}"))),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            match s {
                "strip" => Ok(self.strip_labs.clone()),
                "backplane" => Ok(self.bot_labs.clone()),
                _ => Err(nyx_physics::Error::invalid(format!("unknown boundary {s}"))),
            }
        }
    }
    let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
    spec.study.params = json!({
        "formulation": "planar_bem",
        "singular_edges": "quarter_point",
        "richardson": true
    });
    spec.params = json!({
        "solver": "iterative",
        "tol": 1e-10,
        "amg_line_smooth": true
    });
    spec.materials.insert("substrate".into(), MaterialAssign::Props([("eps_r".into(), 11.45)].into()));
    spec.materials.insert("vacuum".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    spec.bcs.push(nyx_physics::BcSpec {
        on: "strip".into(),
        kind: "terminal".into(),
        params: json!({ "name": "strip" }),
    });
    spec.bcs.push(nyx_physics::BcSpec {
        on: "backplane".into(),
        kind: "ground".into(),
        params: serde_json::Value::Null,
    });

    let resolver = SimpleResolver { sub_dom, vac_dom, strip_labs, bot_labs };
    let mut ctx = Ctx::new(&resolver);
    let res_env = Dfs.solve(&Mesh::from_solver_mesh(&g.mesh).unwrap(), &spec, &mut ctx).unwrap();
    let c_env = CapacitanceMatrixData::of(&res_env).unwrap().maxwell_matrix_f[0][0];
    assert_eq!(c_env.to_bits(), c_bem.to_bits());
}

#[test]
fn test_fr19_3d_isolated_disk_open_robin_bc_vs_analytic() {
    // 3D isolated conducting circular disk of radius a = 100 um at z = 0 in free space (eps_r = 1.0)
    // inside a compact box [-R_box, R_box]^3 with half-width R_box = 3*a = 300 um.
    // Exact free-space capacitance: C_exact = 8 * eps_0 * a = 7.08335025 fF.
    let a = 100.0e-6;
    let r_box = 3.0 * a;
    let c_exact = 8.0 * EPS0 * a;

    let mut m = Model::new_3d();
    let lower = m.box_("lower", (-r_box, -r_box, -r_box), (r_box, r_box, 0.0));
    let upper = m.box_("upper", (-r_box, -r_box, 0.0), (r_box, r_box, r_box));
    let disk_c = m.circle("disk_c", (0.0, 0.0), a);
    m.sheet("disk", disk_c, 0.0);
    m.overlay(&[("lower", lower), ("upper", upper)]);
    m.size("sz_disk_edge", Selector::parse("disk:edge").unwrap(), Some(Expr::Num(9.0e-6)), None, None);
    m.size("sz_disk_face", Selector::Name("disk".into()), Some(Expr::Num(16.0e-6)), None, None);

    let g = mesh_model(&m.list, &MeshPlan::simple(55.0e-6)).unwrap();
    let disk_labs = g.boundary("disk");
    assert!(!disk_labs.is_empty());
    let disk_set: BTreeSet<i32> = disk_labs.iter().copied().collect();
    let ext = g.mesh.exterior_facets();
    let mut outer_labs = Vec::new();
    for (f, &is_ext) in ext.iter().enumerate() {
        let lab = g.mesh.facet_labels[f];
        if is_ext && !disk_set.contains(&lab) && !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }
    assert!(!outer_labs.is_empty());

    // 1. Solve with spherical 1/r asymptotic Robin open boundary (and NO ground reference BC!)
    let mut prob_open = Problem::new(Study::Electrostatic, 2).with_singular_edges("quarter_point");
    prob_open.materials.default.insert("eps_r".into(), 1.0);
    prob_open.bcs = vec![
        Bc { labels: disk_labs.clone(), kind: BcKind::Terminal("disk".into()) },
        Bc { labels: outer_labs.clone(), kind: BcKind::Open { origin: None } },
    ];
    let res_open = solve_electrostatics(&g.mesh, &prob_open).unwrap();
    let cap_open = CapacitanceMatrixData::of(&res_open).unwrap();
    let c_open = cap_open.maxwell_matrix_f[0][0];
    let rel_err_open = (c_open - c_exact).abs() / c_exact;
    eprintln!(
        "FR-19 3D isolated disk (a=100 um, R_box=3a) with open Robin BC: C_open = {c_open:.6e} F, C_exact = {c_exact:.6e} F, rel_err = {:.4}%",
        rel_err_open * 100.0
    );
    assert!(
        rel_err_open < 0.008,
        "FR-19 3D isolated disk open-boundary capacitance relative error {:.4}% exceeds 0.8% threshold (C_open={c_open:e}, C_exact={c_exact:e})",
        rel_err_open * 100.0
    );
    assert_eq!(res_open.scalar("maxwell_raw_asymmetry_rel"), Some(0.0));

    // Compare against grounded outer box walls at R_box = 3a (which suffers ~23% artificial wall capacitance):
    let mut prob_gnd = prob_open.clone();
    prob_gnd.bcs = vec![
        Bc { labels: disk_labs.clone(), kind: BcKind::Terminal("disk".into()) },
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
    ];
    let res_gnd = solve_electrostatics(&g.mesh, &prob_gnd).unwrap();
    let c_gnd = CapacitanceMatrixData::of(&res_gnd).unwrap().maxwell_matrix_f[0][0];
    let rel_err_gnd = (c_gnd - c_exact).abs() / c_exact;
    assert!(rel_err_gnd > 0.15, "grounded box at R_box=3a should overestimate C by >15%, got {rel_err_gnd}");

    // 2. Verify envelope decoding of `kind: "open"` with default origin and explicit origin [0, 0, 0]:
    struct DiskResolver {
        lower_dom: i32,
        upper_dom: i32,
        disk_labs: Vec<i32>,
        outer_labs: Vec<i32>,
    }
    impl LabelResolver for DiskResolver {
        fn domains(&self, r: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            match r {
                "lower" => Ok(vec![self.lower_dom]),
                "upper" => Ok(vec![self.upper_dom]),
                _ => Err(nyx_physics::Error::invalid(format!("unknown domain {r}"))),
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            match s {
                "disk" => Ok(self.disk_labs.clone()),
                "outer" => Ok(self.outer_labs.clone()),
                _ => Err(nyx_physics::Error::invalid(format!("unknown boundary {s}"))),
            }
        }
    }
    let resolver = DiskResolver {
        lower_dom: g.domain("lower").unwrap(),
        upper_dom: g.domain("upper").unwrap(),
        disk_labs: disk_labs.clone(),
        outer_labs: outer_labs.clone(),
    };
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
    spec.study.params = json!({ "singular_edges": "quarter_point" });
    spec.materials.insert("lower".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    spec.materials.insert("upper".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    spec.bcs.push(nyx_physics::BcSpec {
        on: "disk".into(),
        kind: "terminal".into(),
        params: json!({ "name": "disk" }),
    });
    spec.bcs.push(nyx_physics::BcSpec {
        on: "outer".into(),
        kind: "open".into(),
        params: json!({ "origin": [0.0, 0.0, 0.0] }),
    });
    let res_env = Dfs.solve(&ex_mesh, &spec, &mut Ctx::new(&resolver)).unwrap();
    let c_env = CapacitanceMatrixData::of(&res_env).unwrap().maxwell_matrix_f[0][0];
    assert!((c_env - c_open).abs() <= 1e-14 * c_open);

    // 3. Reject coincident origin on an open boundary face and reject 2D open boundary:
    let mut prob_bad_origin = prob_open.clone();
    prob_bad_origin.bcs[1].kind = BcKind::Open { origin: Some([r_box, 0.0, 0.0]) };
    let err_origin = solve_electrostatics(&g.mesh, &prob_bad_origin).unwrap_err();
    assert!(err_origin.contains("origin coincides"), "unexpected error for coincident origin: {err_origin}");

    let mesh2d = femcore::testkit::rect(1.0, 1.0, 2, 2);
    let mut prob_2d = Problem::new(Study::Electrostatic, 1);
    prob_2d.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Terminal("t".into()) },
        Bc { labels: vec![4], kind: BcKind::Open { origin: None } },
    ];
    let err_2d = solve_electrostatics(&mesh2d, &prob_2d).unwrap_err();
    assert!(err_2d.contains("3D"), "expected 2D rejection error mentioning 3D, got: {err_2d}");
}

#[test]
fn test_fr19_and_fr23_raw_asymmetry_scalar_and_silicon_4k_tls_interfaces() {
    // 1. Multi-terminal P2 AMG-PCG solve emits "maxwell_raw_asymmetry_rel" < 1e-7
    let (g, mut prob_amg) = build_transmon_3port_model(85e-6, 20e-6, 120e-6);
    prob_amg.solver = LinearSolver::Iterative;
    prob_amg.tol = 1e-10;
    prob_amg.amg_line_smooth = true;
    let res_amg = solve_electrostatics(&g.mesh, &prob_amg).unwrap();
    let raw_asym = res_amg
        .scalar("maxwell_raw_asymmetry_rel")
        .expect("multi-terminal electrostatic solve must emit maxwell_raw_asymmetry_rel");
    eprintln!("FR-19 3-terminal P2 AMG-PCG maxwell_raw_asymmetry_rel = {raw_asym:.3e}");
    assert!((0.0..1.0e-7).contains(&raw_asym), "expected maxwell_raw_asymmetry_rel in [0, 1e-7), got {raw_asym:e}");

    // 2. FR-23: Solving through Dfs::solve with named material "silicon_4k" populates TlsMeta::interfaces (ms and sa)
    struct TransmonResolver<'a> {
        g: &'a MeshedGeom,
        outer_labs: Vec<i32>,
    }
    impl LabelResolver for TransmonResolver<'_> {
        fn domains(&self, r: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            self.g.domain(r).map(|d| vec![d]).ok_or_else(|| nyx_physics::Error::invalid(r))
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            if s == "outer" {
                return Ok(self.outer_labs.clone());
            }
            let v = self.g.boundary(s);
            if v.is_empty() {
                Err(nyx_physics::Error::invalid(s))
            } else {
                Ok(v)
            }
        }
    }
    let outer_labs = prob_amg.bcs[3].labels.clone();

    let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
    spec.materials.insert("substrate".into(), MaterialAssign::Named("silicon_4k".into()));
    spec.materials.insert("vacuum".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
    for name in ["pad_1", "pad_2", "coupler"] {
        spec.bcs.push(nyx_physics::BcSpec {
            on: name.into(),
            kind: "terminal".into(),
            params: json!({ "name": name }),
        });
    }
    spec.bcs.push(nyx_physics::BcSpec { on: "outer".into(), kind: "open".into(), params: serde_json::Value::Null });
    let resolver = TransmonResolver { g: &g, outer_labs };
    let mut ctx = Ctx::new(&resolver);
    let mat_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    ctx.library = Some(std::sync::Arc::new(femcore::material::Library::load_dir(&mat_dir).unwrap()));

    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let res = Dfs.solve(&ex_mesh, &spec, &mut ctx).unwrap();
    let tls = TlsMeta::of(&res).expect("TlsMeta must be attached to SolveResult::scalars");
    let ifaces = tls.interfaces.expect("TlsMeta::interfaces must be auto-populated from silicon_4k.toml");
    let ms = ifaces.ms.expect("silicon_4k must populate [interfaces.ms]");
    let sa = ifaces.sa.expect("silicon_4k must populate [interfaces.sa]");
    assert!((ms.thickness_m - 2.0e-9).abs() < 1e-18);
    assert!((ms.eps_r - 9.8).abs() < 1e-12);
    assert!((ms.tan_delta - 1.5e-3).abs() < 1e-12);
    assert!((sa.thickness_m - 2.0e-9).abs() < 1e-18);
    assert!((sa.eps_r - 3.9).abs() < 1e-12);
    assert!((sa.tan_delta - 1.5e-3).abs() < 1e-12);

    // Also verify planar_bem populates TlsMeta::interfaces with BcKind::Open present:
    spec.study.params = json!({ "formulation": "planar_bem" });
    let res_bem = Dfs.solve(&ex_mesh, &spec, &mut ctx).unwrap();
    let tls_bem = TlsMeta::of(&res_bem).expect("planar_bem TlsMeta must be attached");
    assert_eq!(tls_bem.interfaces.as_ref().and_then(|i| i.ms.clone()), Some(ms));
    assert_eq!(tls_bem.interfaces.as_ref().and_then(|i| i.sa.clone()), Some(sa));
}
