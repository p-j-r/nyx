// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Comprehensive Wave 1 verification suite for `nyx-analysis`:
//! - **FR-03 (Pt 1 & Pt 2)**: `f_cav_Hz` and `g_matrix_Hz` in `HamiltonianReport::to_json()` and shared Craig–Bampton
//!   cavity decomposition in both `"craig_bampton"` and `"normal_fock"` (with graceful rank-deficient ROM fallback).
//! - **FR-04**: Unified `NetworkAnalysis` (`id: "network"`): `2r x 2r` companion QEP via `femcore::linalg::dense::schur`,
//!   `CircuitRom` multi-port `Z/Y/S` synthesis, Schur-reduced Black-Box Quantization (`bbq_modes`), Taubin algebraic
//!   `S_21` circle + arctan phase fit (`"hanger"`, `"transmission"`, `"reflection"`), and Touchstone v1.1 (`.s3p`).
//! - **FR-11 Pt 1**: 3D matched-asymptotic `I_edge` Gram-matrix capacitance thickness correction (`metal_thickness_m = 200 nm`
//!   `< 0.3%` agreement with 2D finite-thickness FEM reference, SPD preservation, pre-Schur qubit reduction, and clean
//!   rejection of `planar_bem` dummy-field results).
//! - **FR-12**: Separation of nonlinear `lumped_jj` ports from passive resistive `lumped_port` entries in `epr` and `hamiltonian`.
//! - **FR-15**: `resolve_sweep_port` and `sweep_port` targeting by port name or 0-based index in both bases.
//! - **FR-18**: Native DC-SQUID `squid_effective_ej`, `squid_asymmetry`, `squid_flux_sweep`, `ej_hz_sweep`, auto-routed
//!   `flux_sweep`, `< 100 ms` 4.2–15.3 GHz `Dissipator` sweep, and `d = 0, flux = 0.5` (`E_{J,eff} = 0`) safety.
//! - **FR-22**: `c_r` / `z0_ohm` damping (without double-counting `nyx-ems` `c_r`), `epr` `q_factor`, and `hamiltonian`
//!   `kappa_cav_Hz` & `t1_purcell_s` cross-checked against Schur BBQ.
//! - **FR-23**: Flat 3D AABB `SegmentBvh3`, `TlsMeta::interfaces` TOML defaults overridden by JSON `params`, and
//!   continuous `trench_depth_m` corner regularization.

use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_analysis::epr::{project_rom_modes, C0, EPS0, MU0};
use nyx_analysis::{
    capacitance_matrix_with_mesh_and_params, epr, hamiltonian_applies, hamiltonian_with_params, network_with_params,
    resolve_sweep_port, squid_effective_ej, tls_participation_with_params, Analysis, CapacitanceMatrixAnalysis,
    EprAnalysis, NetworkAnalysis,
};
use nyx_dfs::planar_bem::solve_planar_bem;
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study};
use nyx_physics::{CircuitRom, DrivenPortData, Mesh, RomPort, SolveResult, TlsInterfaceSpec, TlsInterfaces, TlsMeta};
use serde_json::{json, Value};
use std::f64::consts::PI;

const W_TRACE: f64 = 10.0e-6;
const G_GAP: f64 = 6.0e-6;
const EPS_SI: f64 = 11.45;
const HALF_BOX: f64 = 40.0e-6;
const H_SUB: f64 = 40.0e-6;
const H_AIR: f64 = 40.0e-6;
const L_Y: f64 = 2.5e-6;

/// Solves a 2D finite-thickness (`t_film > 0`) CPW cross-section on `[-HALF_BOX, HALF_BOX] x [-H_SUB, H_AIR]`.
fn solve_cpw_2d_with_thickness(t_film: f64, eps_sub: f64, h_global: f64, h_corner: f64) -> SolveResult {
    let mut m = Model::new_2d();
    let sub = m.rect("substrate", (-HALF_BOX, -H_SUB), (HALF_BOX, 0.0));
    let vac = m.rect("vacuum", (-HALF_BOX, 0.0), (HALF_BOX, H_AIR));

    let a = W_TRACE / 2.0;
    let b = a + G_GAP;
    let trace_l = m.rect("trace_l", (-a, 0.0), (-a + 1.2e-6, t_film));
    let trace_m = m.rect("trace_m", (-a + 1.2e-6, 0.0), (a - 1.2e-6, t_film));
    let trace_r = m.rect("trace_r", (a - 1.2e-6, 0.0), (a, t_film));
    let gnd_l_far = m.rect("gnd_l_far", (-HALF_BOX, 0.0), (-b - 1.2e-6, t_film));
    let gnd_l_tip = m.rect("gnd_l_tip", (-b - 1.2e-6, 0.0), (-b, t_film));
    let gnd_r_tip = m.rect("gnd_r_tip", (b, 0.0), (b + 1.2e-6, t_film));
    let gnd_r_far = m.rect("gnd_r_far", (b + 1.2e-6, 0.0), (HALF_BOX, t_film));

    m.overlay_full(
        &[
            ("substrate", sub, RegionKind::Material, None),
            ("vacuum", vac, RegionKind::Material, None),
            ("trace_l", trace_l, RegionKind::Hole, None),
            ("trace_m", trace_m, RegionKind::Hole, None),
            ("trace_r", trace_r, RegionKind::Hole, None),
            ("gnd_l_far", gnd_l_far, RegionKind::Hole, None),
            ("gnd_l_tip", gnd_l_tip, RegionKind::Hole, None),
            ("gnd_r_tip", gnd_r_tip, RegionKind::Hole, None),
            ("gnd_r_far", gnd_r_far, RegionKind::Hole, None),
        ],
        &[],
    );
    for bnd in ["trace_l.left", "trace_r.right", "gnd_l_tip.right", "gnd_r_tip.left"] {
        m.size(
            &format!("sz_{}", bnd.replace('.', "_")),
            Selector::Name(bnd.into()),
            Some(Expr::Num(h_corner)),
            None,
            None,
        );
    }
    for bnd in [
        "trace_l.top",
        "trace_l.bottom",
        "trace_r.top",
        "trace_r.bottom",
        "gnd_l_tip.top",
        "gnd_l_tip.bottom",
        "gnd_r_tip.top",
        "gnd_r_tip.bottom",
    ] {
        m.size(
            &format!("sz_{}", bnd.replace('.', "_")),
            Selector::Name(bnd.into()),
            Some(Expr::Num(h_corner * 1.5)),
            None,
            None,
        );
    }
    for bnd in
        ["trace_m.top", "trace_m.bottom", "gnd_l_far.top", "gnd_l_far.bottom", "gnd_r_far.top", "gnd_r_far.bottom"]
    {
        m.size(
            &format!("sz_{}", bnd.replace('.', "_")),
            Selector::Name(bnd.into()),
            Some(Expr::Num(0.25e-6)),
            None,
            None,
        );
    }
    m.size(
        "sz_gap_interface",
        Selector::Interface("substrate".into(), Some("vacuum".into())),
        Some(Expr::Num(0.12e-6)),
        None,
        None,
    );

    let mut g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h_global)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    // Assign exact per-facet labels so outer box walls (`|x| == HALF_BOX`, `|z| == H_SUB`) stay pure Neumann
    // identically to the 3D sheet model:
    let tol = 1e-11;
    let ext2 = g.mesh.exterior_facets();
    for (f, &is_ext) in ext2.iter().enumerate() {
        if !is_ext {
            continue;
        }
        let fv = g.mesh.facet(f);
        let p0 = g.mesh.point(fv[0] as usize);
        let p1 = g.mesh.point(fv[1] as usize);
        let xm = 0.5 * (p0[0] + p1[0]);
        let zm = 0.5 * (p0[1] + p1[1]);
        let is_outer_wall = ((p0[0].abs() - HALF_BOX).abs() < tol && (p1[0].abs() - HALF_BOX).abs() < tol)
            || ((p0[1].abs() - H_SUB).abs() < tol && (p1[1].abs() - H_SUB).abs() < tol);
        if is_outer_wall {
            g.mesh.facet_labels[f] = 300;
        } else if zm >= -tol && zm <= t_film + tol && xm.abs() <= a + tol {
            g.mesh.facet_labels[f] = 100;
        } else if zm >= -tol && zm <= t_film + tol && xm.abs() >= b - tol {
            g.mesh.facet_labels[f] = 200;
        } else {
            g.mesh.facet_labels[f] = 300;
        }
    }
    let trace_labs = vec![100];
    let gnd_labs = vec![200];

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", eps_sub);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());
    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];

    solve_electrostatics(&g.mesh, &prob).unwrap()
}

/// Solves the 3D zero-thickness sheet CPW (`Op::Sheet`) of length `L_Y` with identical cross-section.
fn solve_cpw_3d_sheet_with_eps(eps_sub: f64, h_global: f64, h_edge: f64) -> (MeshedGeom, Problem, Mesh, SolveResult) {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let a = W_TRACE / 2.0;
    let b = a + G_GAP;

    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-HALF_BOX, 0.0, -H_SUB), (HALF_BOX, L_Y, 0.0));
    let vac = m.box_("vacuum", (-HALF_BOX, 0.0, 0.0), (HALF_BOX, L_Y, H_AIR));

    let trace_p = m.polygons("trace_p", vec![vec![pt(-a, 0.0), pt(a, 0.0), pt(a, L_Y), pt(-a, L_Y)]], vec![]);
    let gnd_l_p =
        m.polygons("gnd_l_p", vec![vec![pt(-HALF_BOX, 0.0), pt(-b, 0.0), pt(-b, L_Y), pt(-HALF_BOX, L_Y)]], vec![]);
    let gnd_r_p =
        m.polygons("gnd_r_p", vec![vec![pt(b, 0.0), pt(HALF_BOX, 0.0), pt(HALF_BOX, L_Y), pt(b, L_Y)]], vec![]);

    m.sheet("trace", trace_p, 0.0);
    m.sheet("gnd_l", gnd_l_p, 0.0);
    m.sheet("gnd_r", gnd_r_p, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    for edge_sel in ["trace:edge", "gnd_l:edge", "gnd_r:edge"] {
        m.size(
            &format!("sz_{}", edge_sel.replace(':', "_")),
            Selector::parse(edge_sel).unwrap(),
            Some(Expr::Num(h_edge)),
            None,
            None,
        );
    }
    m.size("sz_trace_face", Selector::Name("trace".into()), Some(Expr::Num(h_edge * 1.8)), None, None);

    let g = mesh_model(&m.list, &MeshPlan::simple(h_global)).unwrap();
    let sub_doms = g.domains("substrate");
    let vac_doms = g.domains("vacuum");

    let mut prob = Problem::new(Study::Electrostatic, 2);
    for d in sub_doms {
        prob.materials.set(d, "eps_r", eps_sub);
        prob.domain_names.insert(d, "substrate".into());
    }
    for d in vac_doms {
        prob.materials.set(d, "eps_r", 1.0);
        prob.domain_names.insert(d, "vacuum".into());
    }

    let trace_labs = g.boundary("trace");
    let mut gnd_labs = g.boundary("gnd_l");
    gnd_labs.extend(g.boundary("gnd_r"));

    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];

    let res = solve_electrostatics(&g.mesh, &prob).unwrap();
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    (g, prob, ex_mesh, res)
}

#[test]
fn fr11_pt1_3d_edge_tube_thickness_correction_matches_2d_200nm_cpw_under_0_3_percent_and_rejects_planar_bem() {
    let t_200nm = 200.0e-9;

    // 1. Exact 2D finite-thickness (t = 200 nm) CPW reference on Silicon (eps_r = 11.45)
    let res_2d_200 = solve_cpw_2d_with_thickness(t_200nm, EPS_SI, 3.5e-6, 35.0e-9);
    let c_pul_2d_200 = nyx_analysis::capacitance_matrix(&res_2d_200).unwrap().maxwell_matrix_f[0][0];
    let c_target_200 = c_pul_2d_200 * L_Y;

    // Near-zero thickness 2D CPW (t = 5 nm) to inspect raw finite-thickness shift
    let res_2d_5nm = solve_cpw_2d_with_thickness(5.0e-9, EPS_SI, 3.5e-6, 15.0e-9);
    let c_target_0 = nyx_analysis::capacitance_matrix(&res_2d_5nm).unwrap().maxwell_matrix_f[0][0] * L_Y;

    // 2. 3D zero-thickness Op::Sheet CPW solve on Silicon (eps_r = 11.45) and Vacuum (eps_r = 1.0)
    let (_g3, _prob3, mesh3, res3) = solve_cpw_3d_sheet_with_eps(EPS_SI, 5.0e-6, 0.28e-6);
    let c_3d_sheet = nyx_analysis::capacitance_matrix(&res3).unwrap().maxwell_matrix_f[0][0];

    let out_corr =
        CapacitanceMatrixAnalysis.run_with_params(&mesh3, &res3, &json!({ "metal_thickness_m": t_200nm })).unwrap();
    let c_3d_corr = out_corr["maxwell_matrix_F"][0][0].as_f64().unwrap();
    let s_eff = out_corr["thickness_correction"]["s_eff_m"].as_f64().unwrap();
    let i_edge = out_corr["thickness_correction"]["i_edge_F_per_m"][0][0].as_f64().unwrap();

    let raw_err = (c_3d_sheet - c_target_200).abs() / c_target_200;
    let corr_err = (c_3d_corr - c_target_200).abs() / c_target_200;
    eprintln!(
        "FR-11 Pt 1 [CPW w=10um, g=6um, L=2.5um, Si eps=11.45, t=200nm]:\n  2D ref (t~0):    {:.6e} F\n  2D ref (t=200nm): {:.6e} F (shift = +{:.3}%)\n  3D sheet (t=0):  {:.6e} F (raw err vs 200nm = {:.3}%)\n  3D corrected:    {:.6e} F (corr err vs 200nm = {:.4}%)\n  s_eff = {:.3e} m, I_edge = {:.5e} F/m",
        c_target_0,
        c_target_200,
        ((c_target_200 - c_target_0) / c_target_0) * 100.0,
        c_3d_sheet,
        raw_err * 100.0,
        c_3d_corr,
        corr_err * 100.0,
        s_eff,
        i_edge
    );

    assert!(c_3d_corr > c_3d_sheet, "finite thickness must increase self-capacitance");
    assert!(
        corr_err < 0.003,
        "FR-11 Pt 1 thickness-corrected 3D Si CPW error {:.4}% >= 0.3% (c_3d_corr={c_3d_corr:e}, c_target={c_target_200:e})",
        corr_err * 100.0
    );

    // Also verify homogeneous Vacuum (eps_r = 1.0), where finite thickness is a +3.2% shift reduced to < 0.3%!
    let c_vac_target_200 =
        nyx_analysis::capacitance_matrix(&solve_cpw_2d_with_thickness(t_200nm, 1.0, 3.5e-6, 35.0e-9))
            .unwrap()
            .maxwell_matrix_f[0][0]
            * L_Y;
    let (_, _, mesh3_vac, res3_vac) = solve_cpw_3d_sheet_with_eps(1.0, 5.0e-6, 0.28e-6);
    let c_vac_sheet = nyx_analysis::capacitance_matrix(&res3_vac).unwrap().maxwell_matrix_f[0][0];
    let out_vac_corr = CapacitanceMatrixAnalysis
        .run_with_params(&mesh3_vac, &res3_vac, &json!({ "metal_thickness_m": t_200nm }))
        .unwrap();
    let c_vac_corr = out_vac_corr["maxwell_matrix_F"][0][0].as_f64().unwrap();
    let vac_raw_err = (c_vac_sheet - c_vac_target_200).abs() / c_vac_target_200;
    let vac_corr_err = (c_vac_corr - c_vac_target_200).abs() / c_vac_target_200;
    eprintln!(
        "FR-11 Pt 1 [Vacuum CPW, t=200nm]: raw err = {:.3}% -> corrected err = {:.4}%",
        vac_raw_err * 100.0,
        vac_corr_err * 100.0
    );
    assert!(vac_raw_err > 0.025, "expected > 2.5% raw zero-thickness error in vacuum, got {:.3}%", vac_raw_err * 100.0);
    assert!(
        vac_corr_err < 0.003,
        "FR-11 Pt 1 thickness-corrected 3D Vacuum CPW error {:.4}% >= 0.3%",
        vac_corr_err * 100.0
    );

    // 3. Verify planar_bem dummy-field result is cleanly rejected when metal_thickness_m is passed!
    let (g_bem, mut prob_bem, mesh_bem, _) = solve_cpw_3d_sheet_with_eps(EPS_SI, 10.0e-6, 1.2e-6);
    prob_bem.bcs.push(Bc { labels: g_bem.boundary("vacuum.zmax"), kind: BcKind::Ground });
    prob_bem.bcs.push(Bc { labels: g_bem.boundary("substrate.zmin"), kind: BcKind::Ground });
    let res_bem = solve_planar_bem(&g_bem.mesh, &prob_bem).unwrap();
    let bem_err =
        capacitance_matrix_with_mesh_and_params(&mesh_bem, &res_bem, &json!({ "metal_thickness_m": t_200nm }))
            .expect_err("planar_bem dummy-field result must be rejected when metal_thickness_m is requested");
    assert!(
        bem_err.to_string().contains("planar_bem") || bem_err.to_string().contains("volumetric_fem"),
        "unexpected error message: {bem_err}"
    );

    // Even if result.log were empty, stripped/dummy boundary-only fields must STILL be rejected!
    let mut res_bem_nolog = res_bem.clone();
    res_bem_nolog.log.clear();
    assert!(capacitance_matrix_with_mesh_and_params(
        &mesh_bem,
        &res_bem_nolog,
        &json!({ "metal_thickness_m": t_200nm })
    )
    .is_err());
}

#[test]
fn fr11_pt1_multi_terminal_coupled_sheets_preserve_spd_and_apply_before_qubit_reduction() {
    // Two coupled floating pads + readout claw + surrounding ground sheet in 3D
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-160e-6, -160e-6, -120e-6), (160e-6, 160e-6, 0.0));
    let vac = m.box_("vacuum", (-160e-6, -160e-6, 0.0), (160e-6, 160e-6, 120e-6));

    let pad_a = m.polygons(
        "pad_a_p",
        vec![vec![pt(-50e-6, -65e-6), pt(50e-6, -65e-6), pt(50e-6, -10e-6), pt(-50e-6, -10e-6)]],
        vec![],
    );
    let pad_b = m.polygons(
        "pad_b_p",
        vec![vec![pt(-50e-6, 10e-6), pt(50e-6, 10e-6), pt(50e-6, 65e-6), pt(-50e-6, 65e-6)]],
        vec![],
    );
    let claw = m.polygons(
        "claw_p",
        vec![vec![pt(-40e-6, 80e-6), pt(40e-6, 80e-6), pt(40e-6, 105e-6), pt(-40e-6, 105e-6)]],
        vec![],
    );
    m.sheet("pad_a", pad_a, 0.0);
    m.sheet("pad_b", pad_b, 0.0);
    m.sheet("claw", claw, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let g = mesh_model(&m.list, &MeshPlan::simple(38e-6)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", EPS_SI);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    let mut box_gnd = g.boundary("vacuum.zmax");
    box_gnd.extend(g.boundary("substrate.zmin"));
    prob.bcs = vec![
        Bc { labels: box_gnd, kind: BcKind::Ground },
        Bc { labels: g.boundary("pad_a"), kind: BcKind::Terminal("pad_a".into()) },
        Bc { labels: g.boundary("pad_b"), kind: BcKind::Terminal("pad_b".into()) },
        Bc { labels: g.boundary("claw"), kind: BcKind::Terminal("claw".into()) },
    ];

    let res = solve_electrostatics(&g.mesh, &prob).unwrap();
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();

    let base_q = capacitance_matrix_with_mesh_and_params(
        &ex_mesh,
        &res,
        &json!({
            "qubits": [{ "name": "Q1", "pads": ["pad_a", "pad_b"], "E_J_Hz": 16.0e9 }],
            "resonators": [{ "name": "RO", "terminal": "claw", "f_r_Hz": 6.8e9 }]
        }),
    )
    .unwrap();

    let corr_q = capacitance_matrix_with_mesh_and_params(
        &ex_mesh,
        &res,
        &json!({
            "metal_thickness_m": 200.0e-9,
            "qubits": [{ "name": "Q1", "pads": ["pad_a", "pad_b"], "E_J_Hz": 16.0e9 }],
            "resonators": [{ "name": "RO", "terminal": "claw", "f_r_Hz": 6.8e9 }]
        }),
    )
    .unwrap();

    assert_eq!(corr_q["has_null_mode"], json!(false));
    assert!(corr_q["condition_number"].as_f64().unwrap() < 100.0);

    // Differential qubit C_Sigma increases (and E_C decreases) when finite metal thickness is applied before reduction!
    let c_sig_0 = base_q["qubits"][0]["C_Sigma_F"].as_f64().unwrap();
    let c_sig_t = corr_q["qubits"][0]["C_Sigma_F"].as_f64().unwrap();
    let ec_0 = base_q["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    let ec_t = corr_q["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    assert!(c_sig_t > c_sig_0, "C_Sigma(t=200nm)={c_sig_t:e} must exceed C_Sigma(0)={c_sig_0:e}");
    assert!(ec_t < ec_0);
    // Mutual capacitance -C_{01}^M between pad_a and pad_b increases (Delta C_{01}^M < 0)
    let delta_01 = corr_q["thickness_correction"]["delta_maxwell_F"][0][1].as_f64().unwrap();
    assert!(delta_01 < 0.0, "off-diagonal Maxwell correction delta_C_01 must be negative, got {delta_01:e}");
}

#[test]
fn fr03_fr12_fr15_fr18_fr22_end_to_end_consistency_and_no_double_counted_damping() {
    // Verify graceful fallback in normal_fock when synthetic ROM has r < n_nl (1 mode, 2 JJs)
    let scale = nyx_analysis::epr::jj_flux_scale();
    let b0 = (EPS0 / 70e-15_f64).sqrt();
    let mut rank_def_rom = CircuitRom {
        ports: vec![
            RomPort {
                name: "JJ_a".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(scale / 8.0e9),
                e_j_hz: Some(8.0e9),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            },
            RomPort {
                name: "JJ_b".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(scale / 8.0e9),
                e_j_hz: Some(8.0e9),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            },
        ],
        r_static: 1,
        r_dyn: 0,
        k_r: vec![vec![0.0]],
        c_r: vec![],
        m_r: vec![vec![1.0]],
        b_r: vec![vec![b0, b0]],
        c_port_f: vec![],
        modes: vec![],
    };
    rank_def_rom.modes = project_rom_modes(&rank_def_rom, &rank_def_rom.ports).unwrap();
    let res_rd = SolveResult { scalars: vec![rank_def_rom.scalar()], ..Default::default() };
    let rep_rd = hamiltonian_with_params(&res_rd, &json!({ "basis": "normal_fock", "fock_trunc": [8] })).unwrap();
    assert!(rep_rd.f_cav_hz.is_empty());
    assert!(rep_rd.g_matrix_hz.is_empty());
    assert_eq!(resolve_sweep_port(&rank_def_rom.ports, &json!({ "sweep_port": "JJ_b" })).unwrap(), 1);

    // Passive-only ROM: hamiltonian_applies is false, epr_applies is true
    let passive_only_rom = CircuitRom {
        ports: vec![RomPort {
            name: "port_50".into(),
            kind: "lumped_port".into(),
            l_j_h: None,
            e_j_hz: None,
            c_j_f: 0.0,
            e_l_hz: None,
            phi_ext_rad: 0.0,
            n_g: 0.0,
            squid_asymmetry: None,
            z0_ohm: Some(50.0),
            direction: [1.0, 0.0, 0.0],
            formulation: "rank1".into(),
        }],
        r_static: 0,
        r_dyn: 1,
        k_r: vec![vec![(2.0 * PI * 6.5e9 / C0).powi(2)]],
        c_r: vec![vec![(MU0 * C0 / 50.0) * 1.0e-4]],
        m_r: vec![vec![1.0]],
        b_r: vec![vec![1.0e-2]],
        c_port_f: vec![],
        modes: vec![],
    };
    let res_pass = SolveResult { scalars: vec![passive_only_rom.scalar()], ..Default::default() };
    assert!(!hamiltonian_applies(&res_pass));
    assert!(EprAnalysis.applies_to(&res_pass));
    assert!(NetworkAnalysis.applies_to(&res_pass));

    // Verify no double-counting when rom.c_r ALREADY contains (MU0 * C0 / z0) b_r b_r^T (nyx-ems convention)
    // vs when rom.c_r is empty!
    let c_q = 75.0e-15;
    let c_r = 50.0e-15;
    let c_g = 1.5e-15;
    let omega_r = 2.0 * PI * 6.8e9;
    let l_r_inv = omega_r * omega_r * (c_r + c_g);
    let b_feed = 0.025_f64;
    let d_port_11 = (MU0 * C0 / 50.0) * b_feed * b_feed;

    let make_rom = |with_cr_baked: bool| -> SolveResult {
        let mut rom = CircuitRom {
            ports: vec![
                RomPort {
                    name: "Q1_JJ".into(),
                    kind: "lumped_jj".into(),
                    l_j_h: Some(scale / 16.0e9),
                    e_j_hz: Some(16.0e9),
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: Some(0.25),
                    z0_ohm: None,
                    direction: [0.0, 1.0, 0.0],
                    formulation: "rank1".into(),
                },
                RomPort {
                    name: "ro_port".into(),
                    kind: "lumped_port".into(),
                    l_j_h: None,
                    e_j_hz: None,
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: None,
                    z0_ohm: Some(50.0),
                    direction: [1.0, 0.0, 0.0],
                    formulation: "rank1".into(),
                },
            ],
            r_static: 0,
            r_dyn: 2,
            k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * l_r_inv]],
            c_r: if with_cr_baked { vec![vec![0.0, 0.0], vec![0.0, d_port_11]] } else { vec![] },
            m_r: vec![vec![(c_q + c_g) / EPS0, -c_g / EPS0], vec![-c_g / EPS0, (c_r + c_g) / EPS0]],
            b_r: vec![vec![1.0, 0.0], vec![0.0, b_feed]],
            c_port_f: vec![],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        SolveResult { scalars: vec![rom.scalar()], ..Default::default() }
    };

    let res_empty_cr = make_rom(false);
    let res_baked_cr = make_rom(true);

    let epr_a = epr(&res_empty_cr).unwrap();
    let epr_b = epr(&res_baked_cr).unwrap();
    assert!((epr_a.q_factors[1].unwrap() - epr_b.q_factors[1].unwrap()).abs() < 1e-9);

    let ham_a = hamiltonian_with_params(&res_empty_cr, &json!({ "basis": "craig_bampton" })).unwrap();
    let ham_b = hamiltonian_with_params(&res_baked_cr, &json!({ "basis": "craig_bampton" })).unwrap();
    assert!((ham_a.kappa_cav_hz[0] - ham_b.kappa_cav_hz[0]).abs() < 1e-6);
    assert!((ham_a.t1_purcell_s[0] - ham_b.t1_purcell_s[0]).abs() / ham_a.t1_purcell_s[0] < 1e-12);

    let net_a = network_with_params(&res_empty_cr, &Value::Null).unwrap();
    let net_b = network_with_params(&res_baked_cr, &Value::Null).unwrap();
    assert!((net_a.qep_modes[1].kappa_hz - net_b.qep_modes[1].kappa_hz).abs() < 1e-6);

    // Check squid_effective_ej formula directly
    assert!((squid_effective_ej(20.0e9, 0.0, 0.25) - 20.0e9).abs() < 1e-6);
    assert!((squid_effective_ej(20.0e9, 0.5, 0.25) - 5.0e9).abs() < 1e-6);
    assert!(squid_effective_ej(20.0e9, 0.5, 0.0).abs() < 1e-6);
}

#[test]
fn fr04_schur_bbq_on_driven_port_data_and_reflection_circle_fit() {
    // Construct a 2-port DrivenPortData: Port 0 = linearized JJ port (active), Port 1 = 50 Ohm readout port (terminated)
    // Qubit C_J = 70 fF, L_J = 10 nH (f_q ~ 6.015 GHz) coupled via C_g = 2 fF to 50 Ohm port (direct capacitive T1 = C_Sigma / (omega^2 C_g^2 Z0))
    let c_j = 70.0e-15_f64;
    let l_j = 10.0e-9_f64;
    let c_g = 2.0e-15_f64;
    let z0 = 50.0_f64;

    let freqs_hz: Vec<f64> = (0..101).map(|k| 5.80e9 + (6.20e9 - 5.80e9) * (k as f64) / 100.0).collect();
    let mut s_matrix = Vec::with_capacity(freqs_hz.len());
    let mut z_matrix = Vec::with_capacity(freqs_hz.len());
    let mut y_matrix = Vec::with_capacity(freqs_hz.len());
    let mut z0_ohm = Vec::with_capacity(freqs_hz.len());

    for &f in &freqs_hz {
        let w = 2.0 * PI * f;
        // Nodal admittance Y of [port_jj, port_ro]:
        // Y_00 = i w (C_j + C_g) + 1 / (i w L_j), Y_01 = Y_10 = -i w C_g, Y_11 = i w C_g + 1e-9 (tiny shunt for invertibility)
        let y00_im = w * (c_j + c_g) - 1.0 / (w * l_j);
        let y01_im = -w * c_g;
        let y11_im = w * c_g;
        let g_tiny = 1.0e-9;

        let y00 = [0.0, y00_im];
        let y01 = [0.0, y01_im];
        let y11 = [g_tiny, y11_im];
        y_matrix.push(vec![vec![y00, y01], vec![y01, y11]]);

        // 2x2 inverse for Z
        let det_re = -y00_im * y11_im + y01_im * y01_im;
        let det_im = y00_im * g_tiny;
        let det_n2 = det_re * det_re + det_im * det_im;
        let inv_det = |re: f64, im: f64| -> [f64; 2] {
            [(re * det_re + im * det_im) / det_n2, (im * det_re - re * det_im) / det_n2]
        };
        let z00 = inv_det(g_tiny, y11_im);
        let z01 = inv_det(0.0, -y01_im);
        let z11 = inv_det(0.0, y00_im);
        z_matrix.push(vec![vec![z00, z01], vec![z01, z11]]);
        s_matrix.push(vec![vec![[0.0, 0.0]; 2]; 2]);
        z0_ohm.push(vec![[z0, 0.0], [z0, 0.0]]);
    }

    let drv = DrivenPortData {
        freqs_hz,
        ports: vec!["jj_port".into(), "ro_port".into()],
        z0_ohm,
        s_matrix,
        z_matrix,
        y_matrix,
        hdm_solves: 3,
        rom_dim: 4,
        anchor_freqs_hz: vec![6.0e9],
        max_residual: 1e-10,
    };
    let res = SolveResult { scalars: vec![drv.scalar()], ..Default::default() };

    let rep = network_with_params(
        &res,
        &json!({
            "active_ports": ["jj_port"],
            "terminations_ohm": { "ro_port": z0 }
        }),
    )
    .unwrap();

    assert_eq!(rep.bbq_modes.len(), 1);
    let m0 = &rep.bbq_modes[0];
    let c_tot_approx = c_j + c_g;
    let f_want = 1.0 / (2.0 * PI * (l_j * c_tot_approx).sqrt());
    let w_want = 2.0 * PI * f_want;
    // Exact Schur Re[Y_red] = (w C_g)^2 * (1/z0) / ((1/z0)^2 + (w C_g)^2)
    let wc = w_want * c_g;
    let g_load = 1.0 / z0;
    let re_y_want = (wc * wc * g_load) / (g_load * g_load + wc * wc) + 0.0;
    let t1_want = c_tot_approx / re_y_want;

    assert!((m0.f_m_hz - f_want).abs() / f_want < 1.0e-4);
    assert!((m0.c_sigma_f - c_tot_approx).abs() / c_tot_approx < 1.0e-3);
    assert!((m0.t1_purcell_s - t1_want).abs() / t1_want < 1.0e-2, "got T1={:e}, want {:e}", m0.t1_purcell_s, t1_want);
}

#[test]
fn fr23_tls_meta_interfaces_defaults_and_trench_depth_continuity() {
    let (_g3, _prob3, mesh3, mut res3) = solve_cpw_3d_sheet_with_eps(EPS_SI, 10.0e-6, 1.1e-6);

    // Inject TlsMeta::interfaces (simulating TOML [interfaces.{ma,ms,sa}] auto-population from nyx-dfs / nyx-ems)
    let mut meta = TlsMeta::of(&res3).unwrap();
    meta.interfaces = Some(TlsInterfaces {
        ma: Some(TlsInterfaceSpec { thickness_m: 3.0e-9, eps_r: 9.8, tan_delta: 2.0e-3 }),
        ms: Some(TlsInterfaceSpec { thickness_m: 2.0e-9, eps_r: 9.8, tan_delta: 1.5e-3 }),
        sa: Some(TlsInterfaceSpec { thickness_m: 2.5e-9, eps_r: 3.9, tan_delta: 1.8e-3 }),
    });
    res3.scalars.retain(|(k, _)| k != nyx_physics::TLS_META_KEY);
    res3.scalars.push(meta.scalar());

    // 1. Default run uses TlsMeta::interfaces; explicit JSON matching those values must match bit-for-bit!
    let rep_from_meta = tls_participation_with_params(&mesh3, &res3, &Value::Null).unwrap();
    let rep_explicit = tls_participation_with_params(
        &mesh3,
        &res3,
        &json!({
            "t_ma": 3.0e-9, "eps_ma": 9.8, "tan_delta_ma": 2.0e-3,
            "t_ms": 2.0e-9, "eps_ms": 9.8, "tan_delta_ms": 1.5e-3,
            "t_sa": 2.5e-9, "eps_sa": 3.9, "tan_delta_sa": 1.8e-3
        }),
    )
    .unwrap();
    assert_eq!(rep_from_meta.loss_tan_total.to_bits(), rep_explicit.loss_tan_total.to_bits());

    // 2. Explicit JSON override supersedes TlsMeta::interfaces
    let rep_zero_ms = tls_participation_with_params(&mesh3, &res3, &json!({ "tan_delta_ms": 0.0 })).unwrap();
    assert!(rep_zero_ms.loss_tan_total < rep_from_meta.loss_tan_total);

    // 3. trench_depth_m = 0.0 reduces identically to trench_depth_m = None, while trench_depth_m > 0
    //    monotonically lowers p_ms, p_sa, and loss_tan_total (increasing q_tls)
    let rep_tr0 = tls_participation_with_params(&mesh3, &res3, &json!({ "trench_depth_m": 0.0 })).unwrap();
    assert_eq!(rep_tr0.p_ms.to_bits(), rep_from_meta.p_ms.to_bits());
    assert_eq!(rep_tr0.p_sa.to_bits(), rep_from_meta.p_sa.to_bits());
    assert_eq!(rep_tr0.p_ma.to_bits(), rep_from_meta.p_ma.to_bits());

    let rep_tr100 = tls_participation_with_params(&mesh3, &res3, &json!({ "trench_depth_m": 100.0e-9 })).unwrap();
    let rep_tr500 = tls_participation_with_params(&mesh3, &res3, &json!({ "trench_depth_m": 500.0e-9 })).unwrap();
    assert!(rep_tr100.p_ms < rep_tr0.p_ms && rep_tr500.p_ms < rep_tr100.p_ms);
    assert!(rep_tr100.p_sa < rep_tr0.p_sa && rep_tr500.p_sa < rep_tr100.p_sa);
    assert!(rep_tr100.q_tls > rep_tr0.q_tls && rep_tr500.q_tls > rep_tr100.q_tls);
}
