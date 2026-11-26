// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! VVUQ Phase 3 integration tests verifying `SolveResult::{diagnostics, error_budget}` across `nyx-dfs`
//! studies (`electrostatic`, `planar_bem`, `heat_stationary`, `heat_transient` BDF2, and `magnetostatic`),
//! strict AMG convergence enforcement (`worst_rel <= tol * 2.0`), and `algebraic_floor_exceeded` warnings.

use femcore::props::LinearSolver;
use femcore::testkit::{block, rect};
use femcore::vvuq::FnSource;
use nyx_dfs::heat::{solve_heat_bdf2_with_source, TransientHeatOptions};
use nyx_dfs::{solve_electrostatics, solve_heat, solve_magnetostatics, Bc, BcKind, Problem, Study};
use nyx_physics::DiagnosticSeverity;

#[test]
fn electro_direct_iterative_and_prager_synge_populate_diagnostics_and_budget() {
    let mesh = block(1.0, 1.0, 1.0, 3, 3, 3);
    let mut p = Problem::new(Study::Electrostatic, 1);
    p.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(5.0) }, Bc { labels: vec![6], kind: BcKind::Ground }];
    p.solver = LinearSolver::Direct;

    let r_dir = solve_electrostatics(&mesh, &p).unwrap();
    let d_dir = r_dir.diagnostics.as_ref().expect("electro direct diagnostics");
    let b_dir = r_dir.error_budget.as_ref().expect("electro direct error_budget");
    assert!(d_dir.converged);
    assert!(d_dir.max_algebraic_residual < 1e-11, "res={:e}", d_dir.max_algebraic_residual);
    assert!((b_dir.algebraic_rel - d_dir.max_algebraic_residual).abs() <= 1e-30);

    p.solver = LinearSolver::Iterative;
    p.tol = 1e-10;
    let r_it = solve_electrostatics(&mesh, &p).unwrap();
    let d_it = r_it.diagnostics.as_ref().expect("electro iterative diagnostics");
    assert!(d_it.converged);
    assert!(d_it.iterations.is_some_and(|it| it > 0));
    assert!(d_it.max_algebraic_residual <= p.tol * 2.0, "worst={:e}", d_it.max_algebraic_residual);

    // Prager-Synge H(div) flux bounds populate `spatial_discretization` and combined budget.
    let mut p_ps = p.clone().with_flux_bounds("report");
    p_ps.solver = LinearSolver::Direct;
    p_ps.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Terminal("top".into()) },
        Bc { labels: vec![6], kind: BcKind::Ground },
    ];
    let r_ps = solve_electrostatics(&mesh, &p_ps).unwrap();
    let b_ps = r_ps.error_budget.as_ref().expect("prager-synge budget");
    let sp = b_ps.spatial_discretization.as_ref().expect("prager_synge_flux_bounds estimate");
    assert_eq!(sp.method, "prager_synge_flux_bounds");
    assert!(sp.is_epistemic_bound);
    assert!(b_ps.combined_num_uncertainty_rel.is_some());
    assert!(b_ps.alg_below_disc_ratio.is_some());

    // Forcing impossible tolerance `1e-22` triggers `worst_rel > tol * 2.0`: falls back cleanly to direct
    // Cholesky with explicit `amg_fallback_to_direct` warning instead of silent AMG acceptance.
    let mesh_coarse = block(1.0, 1.0, 1.0, 2, 2, 2);
    let mut p_strict = Problem::new(Study::Electrostatic, 1);
    p_strict.bcs =
        vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![6], kind: BcKind::Ground }];
    p_strict.solver = LinearSolver::Iterative;
    p_strict.tol = 1e-22;
    let r_fb = solve_electrostatics(&mesh_coarse, &p_strict).unwrap();
    let d_fb = r_fb.diagnostics.as_ref().expect("fallback diagnostics");
    assert!(
        d_fb.warnings.iter().any(|w| w.code == "amg_fallback_to_direct" || w.code == "amg_not_converged"),
        "expected AMG fallback or non-convergence warning when tol=1e-22, got {d_fb:?}"
    );
}

#[test]
fn planar_bem_and_algebraic_floor_exceeded_are_surfaced() {
    let mut m = geom::model::Model::new_3d();
    let box_vol = m.box_("vac", (-200e-6, -200e-6, -200e-6), (200e-6, 200e-6, 200e-6));
    let pad = m.rect("pad", (-40e-6, -40e-6), (40e-6, 40e-6));
    m.sheet("sig", pad, 0.0);
    m.overlay(&[("vac", box_vol)]);
    let g = femcore::meshing::mesh_model(&m.list, &femcore::meshing::MeshPlan::simple(100e-6)).unwrap();
    let sig_lab = g.boundary("sig")[0];

    let mut prob = Problem::new(Study::Electrostatic, 1).with_formulation("planar_bem").with_richardson(true);
    prob.bcs = vec![Bc { labels: vec![sig_lab], kind: BcKind::Terminal("sig".into()) }];
    let r = solve_electrostatics(&g.mesh, &prob).unwrap();
    let d = r.diagnostics.as_ref().expect("planar_bem diagnostics");
    let b = r.error_budget.as_ref().expect("planar_bem budget");
    assert!(d.converged);
    assert_eq!(d.solver_kind, "dense_eigh");
    let sp = b.spatial_discretization.as_ref().expect("richardson 1->4 estimate");
    assert_eq!(sp.method, "richardson_1_to_4");
    assert!(sp.rel_bound >= 0.0 && sp.rel_bound.is_finite());

    // Loose iterative tol (~1e-2) with flux_bounds (~1e-1 or 1e-2) vs `tol`: if eps_alg > 0.01 * eps_disc,
    // `algebraic_floor_exceeded` warning is emitted.
    let mesh = block(1.0, 1.0, 1.0, 4, 4, 4);
    let mut p_floor = Problem::new(Study::Electrostatic, 2).with_flux_bounds("report");
    p_floor.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Terminal("top".into()) },
        Bc { labels: vec![6], kind: BcKind::Ground },
    ];
    p_floor.solver = LinearSolver::Iterative;
    p_floor.tol = 1e-1;
    let rf = solve_electrostatics(&mesh, &p_floor).unwrap();
    let df = rf.diagnostics.as_ref().unwrap();
    let bf = rf.error_budget.as_ref().unwrap();
    if let Some(ratio) = bf.alg_below_disc_ratio {
        if ratio > 0.01 {
            assert!(
                df.warnings
                    .iter()
                    .any(|w| w.code == "algebraic_floor_exceeded" && w.severity == DiagnosticSeverity::Warning),
                "expected algebraic_floor_exceeded warning when ratio={ratio:e}, got {:?}",
                df.warnings
            );
        }
    }
}

#[test]
fn heat_stationary_and_bdf2_report_amg_telemetry_and_temporal_error() {
    let mesh = rect(1.0, 1.0, 8, 8);
    let mut p_st = Problem::new(Study::HeatStationary, 1);
    p_st.solver = LinearSolver::Iterative;
    p_st.materials.default.insert("k".into(), 0.5);
    p_st.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![4], kind: BcKind::Temperature(350.0) },
    ];
    let r_st = solve_heat(&mesh, &p_st).unwrap();
    let d_st = r_st.diagnostics.as_ref().expect("heat stationary diagnostics");
    let b_st = r_st.error_budget.as_ref().expect("heat stationary error_budget");
    assert!(d_st.converged);
    assert!(d_st.iterations.is_some_and(|it| it > 0));
    assert!(b_st.algebraic_rel < 1e-8);

    // BDF2 transient heat via AmgSystem must report AMG iterations and temporal step-halving error estimate.
    let mut p_tr = Problem::new(Study::HeatTransient { t_end: 0.2, n_out: 2, t0: "sin(pi*x)*sin(pi*z)".into() }, 1);
    p_tr.solver = LinearSolver::Iterative;
    p_tr.tol = 1e-10;
    p_tr.materials.default.insert("k".into(), 0.5);
    p_tr.materials.default.insert("rho".into(), 1.0);
    p_tr.materials.default.insert("Cp".into(), 1.0);
    p_tr.bcs = vec![Bc { labels: vec![1, 2, 3, 4], kind: BcKind::Temperature(0.0) }];
    let source = FnSource {
        volume_rhs: |x: [f64; 2], _t: f64, _d: u32| [10.0 * (std::f64::consts::PI * x[0]).sin(), 0.0, 0.0],
        dirichlet: |_, _, _| [0.0; 3],
        neumann_robin: |_, _, _, _| [0.0; 3],
    };
    let r_tr = solve_heat_bdf2_with_source::<2, _>(
        &mesh,
        &p_tr,
        &source,
        TransientHeatOptions { bdf2_substeps: 2, force_bdf2: true },
    )
    .unwrap();
    let d_tr = r_tr.diagnostics.as_ref().expect("heat bdf2 diagnostics");
    let b_tr = r_tr.error_budget.as_ref().expect("heat bdf2 error_budget");
    assert!(d_tr.converged);
    assert!(d_tr.iterations.is_some_and(|it| it > 0));
    assert!(b_tr.temporal_rel.is_some_and(|e| e.is_finite() && e >= 0.0));
}

#[test]
fn magneto_paths_populate_diagnostics_and_budget() {
    // 1. cross_section on thin 2D strip
    let mesh2 = rect(1.0e-3, 0.4e-3, 6, 4);
    let mut p_cs = Problem::new(Study::Magnetostatic { formulation: "cross_section".into() }, 2);
    p_cs.bcs = vec![
        Bc {
            labels: vec![1],
            kind: BcKind::CurrentTerminal { name: "in".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc { labels: vec![4], kind: BcKind::Ground },
    ];
    let r_cs = solve_magnetostatics(&mesh2, &p_cs).unwrap();
    let d_cs = r_cs.diagnostics.as_ref().expect("magneto cross_section diagnostics");
    assert!(d_cs.converged && d_cs.max_algebraic_residual < 1e-10);
    assert!(r_cs.error_budget.is_some());

    // 2. planar_stream on 2D rectangle (London penetration sheet)
    let mut p_ps = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_ps.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::FluxLoop { name: "left".into(), flux_wb: 1e-12 } },
        Bc { labels: vec![4], kind: BcKind::Ground },
    ];
    let r_ps = solve_magnetostatics(&mesh2, &p_ps).unwrap();
    let d_ps = r_ps.diagnostics.as_ref().expect("magneto planar_stream diagnostics");
    assert!(d_ps.converged && d_ps.max_algebraic_residual < 1e-10);
    assert!(r_ps.error_budget.is_some());

    // 3. mixed_3d on small 3D square superconducting washer loop at z = 0
    let pt = |x: f64, y: f64| [geom::expr::Expr::Num(x), geom::expr::Expr::Num(y)];
    let mut m3 = geom::model::Model::new_3d();
    let sub = m3.box_("substrate", (-60e-6, -60e-6, -30e-6), (60e-6, 60e-6, 0.0));
    let vac = m3.box_("vacuum", (-60e-6, -60e-6, 0.0), (60e-6, 60e-6, 30e-6));
    let outer_ring = vec![pt(-35e-6, -35e-6), pt(35e-6, -35e-6), pt(35e-6, 35e-6), pt(-35e-6, 35e-6)];
    let inner_hole = vec![pt(-15e-6, -15e-6), pt(15e-6, -15e-6), pt(15e-6, 15e-6), pt(-15e-6, 15e-6)];
    let washer_poly = m3.polygons("washer_poly", vec![outer_ring], vec![inner_hole]);
    m3.sheet("washer", washer_poly, 0.0);
    m3.overlay(&[("substrate", sub), ("vacuum", vac)]);
    let g3 = femcore::meshing::mesh_model(&m3.list, &femcore::meshing::MeshPlan::simple(40e-6)).unwrap();
    let washer_labs = g3.boundary("washer");
    let mut p_m3 = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
    p_m3.bcs = vec![
        Bc { labels: washer_labs.clone(), kind: BcKind::LondonSheet { l_k_sq: 1.5e-12 } },
        Bc { labels: washer_labs, kind: BcKind::FluxLoop { name: "washer_loop".into(), flux_wb: 1.0 } },
    ];
    let r_m3 = solve_magnetostatics(&g3.mesh, &p_m3).unwrap();
    let d_m3 = r_m3.diagnostics.as_ref().expect("magneto mixed_3d diagnostics");
    assert!(d_m3.converged && d_m3.max_algebraic_residual.is_finite());
    assert!(r_m3.error_budget.is_some());
}
