// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! VVUQ Phase 3 integration tests verifying `SolveResult::{diagnostics, error_budget}` across all
//! `nyx-sms` solvers (`solid` static/eigen, `axisymmetric`, `waveguide`, `floquet`, `pml`, `lrb`,
//! `piezo` static/eigen, and `piezo_driven` HDM & PROM).

use femcore::material::{Anisotropy, Stiffness, Tensor};
use femcore::props::{LinearSolver, Periodic, Pml, PmlDir};
use femcore::testkit::{block, rect};
use nyx_sms::{run, solve_solid, Bc, BcKind, ConductorRole, Plane, Problem, Study};

fn base_steel(study: Study) -> Problem {
    let mut p = Problem::new(study, 2);
    p.materials.default.insert("E".into(), 200e9);
    p.materials.default.insert("nu".into(), 0.3);
    p.materials.default.insert("rho".into(), 7850.0);
    p.bcs = vec![Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    p
}

fn assert_diag(r: &nyx_physics::SolveResult, label: &str) {
    let d = r.diagnostics.as_ref().unwrap_or_else(|| panic!("{label}: missing diagnostics"));
    let b = r.error_budget.as_ref().unwrap_or_else(|| panic!("{label}: missing error_budget"));
    assert!(d.converged, "{label}: expected converged diagnostics, got {d:?}");
    assert!(!d.solver_kind.is_empty(), "{label}: empty solver_kind");
    assert!(d.max_algebraic_residual.is_finite(), "{label}: non-finite residual");
    assert!(b.algebraic_rel.is_finite(), "{label}: non-finite budget residual");
}

#[test]
fn sms_solid_axisym_waveguide_pml_lrb_floquet_populate_diagnostics() {
    let m2 = rect(0.4, 0.1, 8, 2);

    // 1. Solid Static Direct + Iterative
    let mut p_dir = base_steel(Study::Static);
    p_dir.plane = Plane::Stress;
    p_dir.solver = LinearSolver::Direct;
    p_dir.bcs.push(Bc { labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1e4]) });
    assert_diag(&solve_solid(&m2, &p_dir).unwrap(), "solid_static_direct");

    let mut p_it = p_dir.clone();
    p_it.solver = LinearSolver::Iterative;
    let r_it = solve_solid(&m2, &p_it).unwrap();
    assert_diag(&r_it, "solid_static_iterative");
    assert!(r_it.diagnostics.as_ref().unwrap().iterations.is_some());

    // 2. Solid Eigenfrequency Direct + LOBPCG
    let mut p_eig = base_steel(Study::Eigenfrequency { n_modes: 3, f_search: 0.0 });
    p_eig.plane = Plane::Stress;
    p_eig.solver = LinearSolver::Direct;
    let r_ed = solve_solid(&m2, &p_eig).unwrap();
    assert_diag(&r_ed, "solid_eigen_direct");
    assert!(!r_ed.diagnostics.as_ref().unwrap().per_mode_residuals.is_empty());

    p_eig.solver = LinearSolver::Iterative;
    let r_ei = solve_solid(&m2, &p_eig).unwrap();
    assert_diag(&r_ei, "solid_eigen_lobpcg");
    assert!(r_ei.diagnostics.as_ref().unwrap().restarts.is_some());

    // 3. Axisymmetric Eigenfrequency
    let m_ax = rect(0.1, 0.05, 4, 2);
    let mut p_ax =
        Problem::new(Study::AxisymmetricEigenfrequency { n_modes: 3, f_search: 1000.0, azimuthal_order: 0 }, 2);
    p_ax.materials.default.insert("E".into(), 200e9);
    p_ax.materials.default.insert("nu".into(), 0.3);
    p_ax.materials.default.insert("rho".into(), 7850.0);
    p_ax.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
    ];
    assert_diag(&solve_solid(&m_ax, &p_ax).unwrap(), "sms_axisym");

    // 4. Waveguide Modes
    let mut p_wg = Problem::new(Study::Modes { f_search: 5.0e6, wavevector: 100.0, n_modes: 2 }, 2);
    p_wg.materials.default.insert("E".into(), 70e9);
    p_wg.materials.default.insert("nu".into(), 0.33);
    p_wg.materials.default.insert("rho".into(), 2700.0);
    assert_diag(&solve_solid(&m_ax, &p_wg).unwrap(), "sms_waveguide");

    // 5. Floquet Periodic
    let mut p_fl = Problem::new(Study::Eigenfrequency { n_modes: 2, f_search: 5000.0 }, 2);
    p_fl.plane = Plane::Stress;
    p_fl.materials.default.insert("E".into(), 200e9);
    p_fl.materials.default.insert("nu".into(), 0.3);
    p_fl.materials.default.insert("rho".into(), 7850.0);
    p_fl.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 0.3 }];
    assert_diag(&solve_solid(&m2, &p_fl).unwrap(), "sms_floquet");

    // 6. PML
    let mut p_pml = base_steel(Study::Eigenfrequency { n_modes: 2, f_search: 1000.0 });
    p_pml.plane = Plane::Stress;
    p_pml.materials.pml.insert(
        1,
        Pml {
            x: PmlDir::Plus,
            y: PmlDir::Off,
            z: PmlDir::Off,
            strength: 0.5,
            order: 2.0,
            re: 1.0,
            clamp_outer: false,
            ..Default::default()
        },
    );
    assert_diag(&solve_solid(&m2, &p_pml).unwrap(), "sms_pml");

    // 7. Low-Reflecting Boundary (LRB / QEP)
    let c_bar = (200e9_f64 / 7850.0_f64).sqrt();
    let mut p_lrb = base_steel(Study::Eigenfrequency { n_modes: 2, f_search: 1.05 * c_bar / (2.0 * 0.4) });
    p_lrb.plane = Plane::Strain;
    p_lrb.materials.default.insert("nu".into(), 0.0);
    p_lrb.bcs.push(Bc { labels: vec![2, 3], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } });
    p_lrb.bcs.push(Bc { labels: vec![4], kind: BcKind::LowReflecting });
    assert_diag(&solve_solid(&m2, &p_lrb).unwrap(), "sms_lrb");
}

#[test]
fn sms_piezo_static_eigen_and_driven_populate_diagnostics_and_passivity() {
    let m3 = block(1.0e-6, 1.0e-6, 1.0e-6, 2, 2, 3);
    let mut p = Problem::new(Study::Static, 1);
    p.materials.default.insert("rho".into(), 3260.0);
    p.materials.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Hexagonal { c11: 410e9, c12: 149e9, c13: 99e9, c33: 389e9, c44: 125e9 }),
            permittivity: Some(Tensor::Diagonal { x: 9.21, y: 9.21, z: 10.12 }),
            piezo_e: Some([
                [0.0, 0.0, 0.0, 0.0, -0.48, 0.0],
                [0.0, 0.0, 0.0, -0.48, 0.0, 0.0],
                [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0],
            ]),
            ..Default::default()
        },
    );
    p.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Potential(2.0) },
    ];
    let mesh = nyx_physics::Mesh::from_solver_mesh(&m3).unwrap();
    let prog = femcore::progress::Progress::new();
    assert_diag(&run(&mesh, &p, &prog).unwrap(), "piezo_static");

    p.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Ground },
    ];
    p.study = Study::Eigenfrequency { n_modes: 3, f_search: 2.5e9 };
    assert_diag(&run(&mesh, &p, &prog).unwrap(), "piezo_eigen");

    // Driven (both HDM and PROM)
    p.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Terminal { name: "top".into(), z0_ohm: 50.0 } },
    ];
    let mut p_hdm = p.clone();
    p_hdm.study = Study::Driven {
        f_min: 2.0e9,
        f_max: 2.4e9,
        n_points: 5,
        freqs_hz: None,
        prom_tol: None,
        max_hdm_solves: 5,
        strict_prom_tol: false,
        seed_freqs_hz: vec![],
        solver: "auto".into(),
        blr_tol: 0.0,
    };
    let r_hdm = run(&mesh, &p_hdm, &prog).unwrap();
    assert_diag(&r_hdm, "piezo_driven_hdm");
    assert_eq!(r_hdm.diagnostics.as_ref().unwrap().per_mode_residuals.len(), 5);

    let mut p_prom = p.clone();
    p_prom.study = Study::Driven {
        f_min: 2.0e9,
        f_max: 2.4e9,
        n_points: 9,
        freqs_hz: None,
        prom_tol: Some(1e-3),
        max_hdm_solves: 6,
        strict_prom_tol: false,
        seed_freqs_hz: vec![],
        solver: "auto".into(),
        blr_tol: 0.0,
    };
    let r_prom = run(&mesh, &p_prom, &prog).unwrap();
    assert_diag(&r_prom, "piezo_driven_prom");
    assert!(
        r_prom.diagnostics.as_ref().unwrap().warnings.iter().any(|w| w.code == "s_matrix_passivity_reciprocity"),
        "expected s_matrix_passivity_reciprocity diagnostic"
    );
    // Unused ConductorRole import verification check
    let _ = ConductorRole::Ground;
}
