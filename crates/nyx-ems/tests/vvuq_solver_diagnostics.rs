// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! VVUQ Phase 3 integration tests verifying `SolveResult::{diagnostics, error_budget}` across all
//! `nyx-ems` studies (`modes`, `modes_complex`, `axisym`, `cavity`, `bloch`, `driven` HDM/PROM/bisection),
//! passivity/reciprocity telemetry, waveport cutoff & eigenresiduals, and ROM diagnostic checks.

use femcore::mesh::SimplexMesh;
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::props::{Periodic, Pml, PmlDir};
use femcore::testkit::{block, rect, rect_split_x};
use nyx_ems::bloch::{assemble_bloch, solve_bloch_eigen};
use nyx_ems::cavity::{assemble_emw_eigen, solve_emw_modes};
use nyx_ems::driven::{assemble_emw_driven, solve_driven_with_options, DrivenOptions};
use nyx_ems::ports::AssembledLumpedPort;
use nyx_ems::rom::build_circuit_rom_checked;
use nyx_ems::{run, LumpedPortBc, Problem, Study, WavePortBc, C0};
use nyx_physics::{Mesh, ModeEnergyRow};

fn assert_diag(r: &nyx_physics::SolveResult, label: &str) {
    let d = r.diagnostics.as_ref().unwrap_or_else(|| panic!("{label}: missing diagnostics"));
    let b = r.error_budget.as_ref().unwrap_or_else(|| panic!("{label}: missing error_budget"));
    assert!(d.converged, "{label}: expected converged=true, got {d:?}");
    assert!(!d.solver_kind.is_empty(), "{label}: empty solver_kind");
    assert!(d.max_algebraic_residual.is_finite(), "{label}: non-finite residual");
    assert!(b.algebraic_rel.is_finite(), "{label}: non-finite budget residual");
}

#[test]
fn ems_modes_axisym_cavity_and_bloch_populate_diagnostics_and_budget() {
    let prog = femcore::progress::Progress::new();

    // 1. Real 2D waveguide modes
    let m2 = rect(0.02286, 0.01016, 6, 3);
    let p_m = Problem::new(Study::Modes { freq: 10.0e9, n_target: 1.0, evanescent: false, n_modes: 2 }, 2);
    let r_m = run(&Mesh::from_solver_mesh(&m2).unwrap(), &p_m, &prog).unwrap();
    assert_diag(&r_m, "ems_modes_real");
    assert!(!r_m.diagnostics.as_ref().unwrap().per_mode_residuals.is_empty());

    // 2. Complex 2D waveguide modes with PML
    let m2_pml = rect_split_x(0.03, 0.015, 0.02, 0.005);
    let mut p_mc = Problem::new(Study::Modes { freq: 12.0e9, n_target: 1.5, evanescent: false, n_modes: 2 }, 2);
    p_mc.materials.set(1, "eps_r", 4.0);
    p_mc.materials.set(2, "eps_r", 1.0);
    p_mc.materials.pml.insert(2, Pml { x: PmlDir::Plus, ..Default::default() });
    let r_mc = run(&Mesh::from_solver_mesh(&m2_pml).unwrap(), &p_mc, &prog).unwrap();
    assert_diag(&r_mc, "ems_modes_complex");
    assert!(r_mc.diagnostics.as_ref().unwrap().restarts.is_some());

    // 3. Axisymmetric cavity modes
    let m_ax = rect(0.02, 0.02, 4, 4);
    let p_ax = Problem::new(Study::AxisymmetricModes { f_search: 8.0e9, azimuthal_order: 0, n_modes: 2 }, 2);
    let r_ax = run(&Mesh::from_solver_mesh(&m_ax).unwrap(), &p_ax, &prog).unwrap();
    assert_diag(&r_ax, "ems_axisym");

    // 4. 3D Cavity (real) and checked CircuitRom
    let m3 = block(0.02286, 0.02, 0.01016, 4, 4, 2);
    let p_cav = Problem::new(Study::Cavity { n_modes: 2, f_search: 8.0e9 }, 2);
    let r_cav = run(&Mesh::from_solver_mesh(&m3).unwrap(), &p_cav, &prog).unwrap();
    assert_diag(&r_cav, "ems_cavity_real");

    let mut p_rom = p_cav.clone();
    p_rom.lumped_ports = vec![LumpedPortBc {
        labels: vec![2],
        name: "lp".into(),
        z0: 50.0,
        direction: [0.0, 0.0, 1.0],
        complex_damping: false,
        mode: "single".into(),
    }];
    let asm_c = assemble_emw_eigen(&m3, &p_rom).unwrap();
    let solved_c = solve_emw_modes(&asm_c).unwrap();
    let modes_row: Vec<ModeEnergyRow> = (0..solved_c.freqs.len())
        .map(|i| ModeEnergyRow {
            freq_hz: solved_c.freqs[i].0.re,
            q_factor: None,
            u_elec_j: 1.0,
            u_mag_vol_j: 0.5,
            u_kin_ibc_j: 0.0,
            u_ind_jj_j: vec![],
            v_mj: vec![],
            i_mj: vec![],
        })
        .collect();
    let (_rom, rom_diag) =
        build_circuit_rom_checked(&m3, &p_rom, &asm_c, &solved_c.vectors, &solved_c.freqs, modes_row, None).unwrap();
    assert!(rom_diag.max_algebraic_residual.is_finite());

    // 5. 3D Bloch QEP (`residuals`/`restarts`/`matvecs` retained)
    let nb = m2.facet_labels.iter().copied().max().unwrap_or(0);
    let m3_bloch: SimplexMesh =
        extrude_y(&m2, &[Layer { thickness: 0.005, n: 1, domain_offset: 0, fill: None }], 0.0).unwrap();
    let (y0, y1) = (nb + 1, nb + 2);
    let mut p_bl = Problem::new(Study::Bloch { freq: 10.0e9, n_target: 0.8, n_modes: 2, brillouin_zone: 1 }, 2);
    p_bl.periodic = vec![Periodic { source: vec![y0], dest: vec![y1], phase: 0.0 }];
    let asm_b = assemble_bloch(&m3_bloch, &p_bl).unwrap();
    let sol_b = solve_bloch_eigen(&asm_b).unwrap();
    assert_eq!(sol_b.residuals.len(), sol_b.vectors.len());
    assert!(sol_b.matvecs > 0);
    let r_bl = run(&Mesh::from_solver_mesh(&m3_bloch).unwrap(), &p_bl, &prog).unwrap();
    assert_diag(&r_bl, "ems_bloch");
}

#[test]
fn ems_driven_hdm_prom_bisection_and_waveport_cutoff_diagnostics() {
    let (a, l, b) = (22.86e-3, 12.0e-3, 10.16e-3);
    let mesh = block(a, l, b, 5, 3, 2);

    // 1. Driven HDM (prom_tol = None) reports non-zero backward residual & per-freq residuals
    let mut p_hdm = Problem::new(
        Study::Driven {
            f_min: 8.0e9,
            f_max: 11.0e9,
            n_points: 4,
            freqs_hz: None,
            prom_tol: None,
            max_hdm_solves: 10,
            strict_prom_tol: false,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 0.0,
        },
        2,
    );
    p_hdm.waveports = vec![
        WavePortBc { labels: vec![2], name: "wp1".into(), mode: 1, n_eff_guess: None, z0_ohm: None },
        WavePortBc { labels: vec![5], name: "wp2".into(), mode: 1, n_eff_guess: None, z0_ohm: None },
    ];
    let r_hdm = run(&Mesh::from_solver_mesh(&mesh).unwrap(), &p_hdm, &femcore::progress::Progress::new()).unwrap();
    assert_diag(&r_hdm, "ems_driven_hdm");
    let d_hdm = r_hdm.diagnostics.as_ref().unwrap();
    assert_eq!(d_hdm.per_mode_residuals.len(), 4);
    assert!(d_hdm.max_algebraic_residual > 0.0, "expected real LU backward residual, got 0");
    assert!(
        d_hdm.warnings.iter().any(|w| w.code == "s_matrix_passivity_reciprocity"),
        "expected passivity/reciprocity diagnostic on driven HDM, got {:?}",
        d_hdm.warnings
    );

    // 2. Above-cutoff waveport assembly reports high mode overlap + eigenresidual, while evaluating
    //    `beta_and_derivs_checked` below cutoff emits explicit `waveport_below_cutoff` warning, and
    //    opt-in `DrivenOptions { resonance_bisection: true, .. }` populates diagnostics cleanly.
    let mut p_prom = p_hdm.clone();
    p_prom.study = Study::Driven {
        f_min: 8.0e9,
        f_max: 11.0e9,
        n_points: 5,
        freqs_hz: None,
        prom_tol: Some(1e-3),
        max_hdm_solves: 6,
        strict_prom_tol: false,
        seed_freqs_hz: vec![],
        solver: "auto".into(),
        blr_tol: 0.0,
    };
    let asm_wp = assemble_emw_driven(&mesh, &p_prom).unwrap();
    assert!(
        asm_wp.waveports.iter().all(|wp| wp.max_eigen_residual.is_finite() && wp.min_mode_overlap > 0.90),
        "waveport assembly telemetry missing: {:?}",
        asm_wp.waveports.iter().map(|w| (w.min_mode_overlap, w.max_eigen_residual)).collect::<Vec<_>>()
    );
    let f_c = C0 / (2.0 * a);
    let (_vals, warn) = asm_wp.waveports[0].beta_and_derivs_checked(2.0 * std::f64::consts::PI * (0.5 * f_c));
    assert!(warn.as_ref().is_some_and(|w| w.code == "waveport_below_cutoff"));

    let r_bis = solve_driven_with_options(
        &mesh,
        &p_prom,
        &DrivenOptions {
            resonance_bisection: true,
            bisection_max_depth: 1,
            bisection_phase_jump_rad: 0.05,
            bisection_s_dip_tol: 0.05,
        },
    )
    .unwrap();
    assert_diag(&r_bis, "ems_driven_prom_bisection");

    // Below-cutoff endpoint sweep tests `waveport_mode_overlap_low` warning surfacing.
    let mut p_low_ov = p_hdm.clone();
    p_low_ov.study = Study::Driven {
        f_min: 5.0e9,
        f_max: 10.0e9,
        n_points: 3,
        freqs_hz: None,
        prom_tol: Some(1e-3),
        max_hdm_solves: 6,
        strict_prom_tol: false,
        seed_freqs_hz: vec![],
        solver: "auto".into(),
        blr_tol: 0.0,
    };
    let asm_low = assemble_emw_driven(&mesh, &p_low_ov).unwrap();
    assert!(
        asm_low.waveports.iter().any(|wp| wp.warnings.iter().any(|w| w.code == "waveport_mode_overlap_low")),
        "expected waveport_mode_overlap_low warning across cutoff transition"
    );
    let _ = std::mem::size_of::<AssembledLumpedPort>();
}
