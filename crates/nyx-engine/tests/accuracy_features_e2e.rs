// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end verification of `nyx3` Phase A electrostatic accuracy features (`AF-1..AF-6`) and all five
//! Cross-Feature Compatibility Rules (`Rules 1..5`) through `Engine::call`.

use geom::gds::encode_ibm_f64;
use nyx_engine::Engine;
use nyx_physics::{CapacitanceMatrixData, Error, TlsMeta, ELECTRO_AMR_KEY};
use serde_json::json;
use std::f64::consts::PI;
use std::path::PathBuf;

const EPS0: f64 = 8.854_187_812_8e-12;

fn temp_dir(tag: &str) -> PathBuf {
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!(
        "nyx-af-e2e-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

struct GdsBuilder {
    buf: Vec<u8>,
}

impl GdsBuilder {
    fn new(user_unit: f64, dbu_m: f64) -> Self {
        let mut b = Self { buf: Vec::new() };
        b.rec(0x00, 0x02, &600i16.to_be_bytes());
        b.rec(0x01, 0x02, &[0u8; 24]);
        b.rec_ascii(0x02, "af_e2e_lib");
        let mut u = Vec::with_capacity(16);
        u.extend_from_slice(&encode_ibm_f64(user_unit));
        u.extend_from_slice(&encode_ibm_f64(dbu_m));
        b.rec(0x03, 0x05, &u);
        b
    }

    fn rec(&mut self, rec_type: u8, data_type: u8, payload: &[u8]) {
        let len = (4 + payload.len()) as u16;
        self.buf.extend_from_slice(&len.to_be_bytes());
        self.buf.push(rec_type);
        self.buf.push(data_type);
        self.buf.extend_from_slice(payload);
    }

    fn rec_ascii(&mut self, rec_type: u8, s: &str) {
        let mut p = s.as_bytes().to_vec();
        if !p.len().is_multiple_of(2) {
            p.push(0);
        }
        self.rec(rec_type, 0x06, &p);
    }

    fn rec_xy(&mut self, pts: &[[i32; 2]]) {
        let mut p = Vec::with_capacity(pts.len() * 8);
        for &[x, y] in pts {
            p.extend_from_slice(&x.to_be_bytes());
            p.extend_from_slice(&y.to_be_bytes());
        }
        self.rec(0x10, 0x03, &p);
    }

    fn begin_cell(&mut self, name: &str) {
        self.rec(0x05, 0x02, &[0u8; 24]);
        self.rec_ascii(0x06, name);
    }

    fn end_cell(&mut self) {
        self.rec(0x07, 0x00, &[]);
    }

    fn boundary(&mut self, layer: i16, datatype: i16, open_pts: &[[i32; 2]]) {
        self.rec(0x08, 0x00, &[]);
        self.rec(0x0D, 0x02, &layer.to_be_bytes());
        self.rec(0x0E, 0x02, &datatype.to_be_bytes());
        let mut closed = open_pts.to_vec();
        if closed.first() != closed.last() {
            closed.push(open_pts[0]);
        }
        self.rec_xy(&closed);
        self.rec(0x11, 0x00, &[]);
    }

    fn text(&mut self, layer: i16, texttype: i16, xy: [i32; 2], label: &str) {
        self.rec(0x0C, 0x00, &[]);
        self.rec(0x0D, 0x02, &layer.to_be_bytes());
        self.rec(0x16, 0x02, &texttype.to_be_bytes());
        self.rec_xy(&[xy]);
        self.rec_ascii(0x19, label);
        self.rec(0x11, 0x00, &[]);
    }

    fn finish(mut self) -> Vec<u8> {
        self.rec(0x04, 0x00, &[]);
        self.buf
    }
}

#[test]
fn af1_multi_pass_amr_converges_without_manual_add_size_and_persists_in_results_h5() {
    let dir = temp_dir("af1");
    let h5_path = dir.join("amr_adapted.h5");

    let mut e = Engine::new();
    e.call("project_new", &json!({ "title": "AF1 AMR", "dim": 3 })).unwrap();
    e.call(
        "add_polygons",
        &json!({
            "name": "p1",
            "polygons": [[[-55e-6, -30e-6], [-10e-6, -30e-6], [-10e-6, 30e-6], [-55e-6, 30e-6]]]
        }),
    )
    .unwrap();
    e.call(
        "add_polygons",
        &json!({
            "name": "p2",
            "polygons": [[[10e-6, -30e-6], [55e-6, -30e-6], [55e-6, 30e-6], [10e-6, 30e-6]]]
        }),
    )
    .unwrap();
    e.call("add_sheet", &json!({ "name": "pad1", "input": "p1", "z": 0.0 })).unwrap();
    e.call("add_sheet", &json!({ "name": "pad2", "input": "p2", "z": 0.0 })).unwrap();
    e.call(
        "add_chip_stack",
        &json!({
            "lateral_pad": 65e-6,
            "t_sub": 80e-6,
            "h_vac_top": 80e-6,
            "h_vac_bot": 0.0
        }),
    )
    .unwrap();
    // Deliberately NO `add_size` rules: AMR starts from coarse global `h` alone!
    e.call(
        "set_mesh",
        &json!({
            "h": 42e-6,
            "adapt_grade": 0.35,
            "adapt_max_cells": 90000
        }),
    )
    .unwrap();
    e.call(
        "add_physics",
        &json!({
            "physics": "dfs",
            "study": "electrostatic",
            "order": 2,
            "study_params": {
                "adapt_passes": 3,
                "adapt_rel_tol": 1e-4,
                "adapt_growth": 1.6,
                "adapt_mode": "sheet_2d_z",
                "adapt_estimator": "flux_jump_goal"
            }
        }),
    )
    .unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "ground" })).unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "pad1", "kind": "terminal", "params": { "name": "c1" } }))
        .unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "pad2", "kind": "terminal", "params": { "name": "c2" } }))
        .unwrap();

    e.call("solve_run", &json!({})).unwrap();
    let (res, _) = e.result("dfs").expect("solved dfs");
    assert!(
        !res.scalars.iter().any(|(k, _)| k == ELECTRO_AMR_KEY),
        "transient ElectroAmrData scalar must be stripped before storing result"
    );
    let cap = CapacitanceMatrixData::of(res).expect("capacitance matrix");
    assert_eq!(cap.adapt_history.len(), 3, "expected 3 AMR pass records, got {:?}", cap.adapt_history);
    assert_eq!(cap.adapt_history[0].max_rel_delta, None);
    assert!(cap.adapt_history[1].max_rel_delta.is_some());
    assert!(cap.adapt_history[2].max_rel_delta.is_some());
    assert!(
        cap.adapt_history[2].n_cells > cap.adapt_history[0].n_cells,
        "AMR must refine cells across passes: {:?}",
        cap.adapt_history
    );
    assert!(
        cap.adapt_history[2].total_eta < cap.adapt_history[0].total_eta,
        "goal-oriented estimator must decay across passes: {:?}",
        cap.adapt_history
    );
    assert_eq!(e.model().unwrap().mesh().n_cells(), cap.adapt_history[2].n_cells);

    // `analyse("capacitance_matrix")` also exposes `adapt_history`
    let ana = e.call("analyse", &json!({ "analysis": "capacitance_matrix", "physics": "dfs" })).unwrap();
    assert_eq!(ana["adapt_history"].as_array().unwrap().len(), 3);

    // Export to HDF5 and round-trip via `open_results` onto a fresh Engine
    let h5_str = h5_path.to_str().unwrap();
    e.call("export_results", &json!({ "path": h5_str })).unwrap();
    let mut reopened = Engine::new();
    reopened.call("open_results", &json!({ "path": h5_str })).unwrap();
    assert_eq!(reopened.model().unwrap().mesh().n_cells(), cap.adapt_history[2].n_cells);
    let probe = reopened
        .call("probe", &json!({ "physics": "dfs", "expr": "V", "points": [[-30e-6, 0.0, 0.0], [0.0, 0.0, 20e-6]] }))
        .unwrap();
    assert!(probe["values"][0].as_f64().unwrap().is_finite());

    // `convergence_check` forces `adapt_passes = 1` on its scratch copies without mutating caller project
    let conv = e.call("convergence_check", &json!({ "physics": "dfs", "factor": 0.8 })).unwrap();
    assert!(conv["worst_relative"].as_f64().unwrap().is_finite());
    assert_eq!(e.project().spec("dfs").unwrap().study.params["adapt_passes"], json!(3));

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn af2_quarter_point_and_corner_grading_improve_single_grid_accuracy_with_finite_fields() {
    let run_case = |h: f64, qp: bool, corner: bool| -> (f64, Engine) {
        let mut e = Engine::new();
        e.call("project_new", &json!({ "title": "AF2 QP", "dim": 3 })).unwrap();
        e.call(
            "add_polygons",
            &json!({
                "name": "pad_poly",
                "polygons": [[[-35e-6, -35e-6], [35e-6, -35e-6], [35e-6, 35e-6], [-35e-6, 35e-6]]]
            }),
        )
        .unwrap();
        e.call("add_sheet", &json!({ "name": "pad", "input": "pad_poly", "z": 0.0 })).unwrap();
        e.call(
            "add_chip_stack",
            &json!({
                "bbox": [-90e-6, -90e-6, 90e-6, 90e-6],
                "t_sub": 70e-6,
                "h_vac_top": 70e-6,
                "h_vac_bot": 0.0
            }),
        )
        .unwrap();
        e.call(
            "set_mesh",
            &json!({
                "h": h,
                "singular_edges": if qp { "quarter_point" } else { "none" },
                "quarter_point_xi": 0.25,
                "corner_grading": corner
            }),
        )
        .unwrap();
        e.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 2 })).unwrap();
        e.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
        e.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
        e.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "ground" })).unwrap();
        e.call(
            "add_condition",
            &json!({ "physics": "dfs", "on": "pad", "kind": "terminal", "params": { "name": "p" } }),
        )
        .unwrap();
        e.call("solve_run", &json!({})).unwrap();
        let c = CapacitanceMatrixData::of(e.result("dfs").unwrap().0).unwrap().maxwell_matrix_f[0][0];
        (c, e)
    };

    let (c_ref, _) = run_case(11e-6, true, true);
    let (c_straight, _) = run_case(28e-6, false, false);
    let (c_qp, mut e_qp) = run_case(28e-6, true, true);

    let err_straight = ((c_straight - c_ref) / c_ref).abs();
    let err_qp = ((c_qp - c_ref) / c_ref).abs();
    assert!(
        err_qp < err_straight,
        "quarter_point + corner_grading error ({err_qp:.4e}) must beat straight P2 ({err_straight:.4e}) vs ref {c_ref:.6e}"
    );
    assert!(e_qp.model().unwrap().solver_mesh().has_qp(), "sheet3d must tag quarter-point edges");

    // Verify all fields and reductions remain strictly finite (no 0/0 NaN at singular rims/corners)
    let vf = e_qp.call("result_field", &json!({ "physics": "dfs", "field": "V", "max_values": 5000 })).unwrap();
    assert!(vf["values"].as_array().unwrap().iter().all(|x| x.as_f64().unwrap().is_finite()));
    let ef = e_qp.call("result_field", &json!({ "physics": "dfs", "field": "E", "max_values": 5000 })).unwrap();
    assert!(ef["values"].as_array().unwrap().iter().all(|x| x.as_f64().unwrap().is_finite()));
    let u_int = e_qp.call("integrate", &json!({ "physics": "dfs", "expr": "eps_r * (Ex^2 + Ey^2 + Ez^2)" })).unwrap();
    assert!(u_int["value"].as_f64().unwrap().is_finite() && u_int["value"].as_f64().unwrap() > 0.0);
}

#[test]
fn af3_farfield_isolated_sphere_and_disk_single_terminal_without_ground() {
    // 1. Isolated conducting sphere of radius R = 50 um inside [-120 um, 120 um]^3 vacuum box with `far_field`
    let r = 50e-6;
    let c_sphere_exact = 4.0 * PI * EPS0 * r; // 5.56325028 fF
    let mut e_sph = Engine::new();
    e_sph.call("project_new", &json!({ "title": "AF3 Sphere", "dim": 3 })).unwrap();
    e_sph
        .call(
            "add_box",
            &json!({ "name": "box_dom", "min": [-125e-6, -125e-6, -125e-6], "max": [125e-6, 125e-6, 125e-6] }),
        )
        .unwrap();
    e_sph.call("add_sphere", &json!({ "name": "sph_body", "center": [0.0, 0.0, 0.0], "r": r })).unwrap();
    e_sph
        .call(
            "add_overlay",
            &json!({
                "regions": [
                    { "name": "Vac", "input": "box_dom", "kind": "Material" },
                    { "name": "BallHole", "input": "sph_body", "kind": "Hole" }
                ]
            }),
        )
        .unwrap();
    e_sph.call("add_size", &json!({ "sel": "name(sph_body.surface)", "h": 12e-6 })).unwrap();
    e_sph.call("set_mesh", &json!({ "h": 24e-6, "curved": true })).unwrap();
    e_sph.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 2 })).unwrap();
    e_sph.call("set_material", &json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } })).unwrap();
    e_sph
        .call(
            "add_condition",
            &json!({ "physics": "dfs", "on": "sph_body.surface", "kind": "terminal", "params": { "name": "sph" } }),
        )
        .unwrap();
    e_sph
        .call(
            "add_condition",
            &json!({ "physics": "dfs", "on": "outer & !name(sph_body.surface)", "kind": "far_field", "params": { "order": 2 } }),
        )
        .unwrap();
    e_sph.call("solve_run", &json!({})).unwrap();
    let c_sph = CapacitanceMatrixData::of(e_sph.result("dfs").unwrap().0).unwrap().maxwell_matrix_f[0][0];
    let rel_sph = ((c_sph - c_sphere_exact) / c_sphere_exact).abs();
    assert!(rel_sph < 0.015, "isolated sphere C = {c_sph:.6e} F vs exact {c_sphere_exact:.6e} F (rel = {rel_sph:.4e})");

    // 2. Isolated circular disk of radius R = 50 um in vacuum using `add_chip_stack(lateral_pad=...)` + `far_field`
    let c_disk_exact = 8.0 * EPS0 * r; // 3.54167513 fF
    let mut e_dsk = Engine::new();
    e_dsk.call("project_new", &json!({ "title": "AF3 Disk", "dim": 3 })).unwrap();
    e_dsk.call("add_circle", &json!({ "name": "disk2d", "center": [0.0, 0.0], "r": r })).unwrap();
    e_dsk.call("add_sheet", &json!({ "name": "disk", "input": "disk2d", "z": 0.0 })).unwrap();
    let stack_out = e_dsk
        .call(
            "add_chip_stack",
            &json!({
                "lateral_pad": 90e-6,
                "t_sub": 130e-6,
                "h_vac_top": 130e-6,
                "h_vac_bot": 0.0
            }),
        )
        .unwrap();
    let pb: Vec<f64> = stack_out["padded_bbox"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    for (got, want) in pb.iter().zip([-140e-6, -140e-6, 140e-6, 140e-6]) {
        assert!((got - want).abs() < 1e-18, "padded_bbox entry {got} != {want}");
    }
    e_dsk.call("add_size", &json!({ "sel": "disk:edge", "h": 8e-6 })).unwrap();
    e_dsk
        .call(
            "set_mesh",
            &json!({
                "h": 32e-6,
                "curved": true,
                "singular_edges": "quarter_point"
            }),
        )
        .unwrap();
    e_dsk.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 2 })).unwrap();
    e_dsk.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 1.0 } })).unwrap();
    e_dsk.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
    e_dsk
        .call(
            "add_condition",
            &json!({ "physics": "dfs", "on": "disk", "kind": "terminal", "params": { "name": "d" } }),
        )
        .unwrap();
    e_dsk
        .call(
            "add_condition",
            &json!({ "physics": "dfs", "on": "outer", "kind": "far_field", "params": { "order": 2 } }),
        )
        .unwrap();
    e_dsk.call("solve_run", &json!({})).unwrap();
    let c_dsk = CapacitanceMatrixData::of(e_dsk.result("dfs").unwrap().0).unwrap().maxwell_matrix_f[0][0];
    let rel_dsk = ((c_dsk - c_disk_exact) / c_disk_exact).abs();
    assert!(rel_dsk < 0.015, "isolated disk C = {c_dsk:.6e} F vs exact {c_disk_exact:.6e} F (rel = {rel_dsk:.4e})");
}

#[test]
fn af4_prager_synge_flux_bounds_report_and_midpoint_synchronize_all_outputs() {
    let build_pocket = |h: f64, edge_h: f64, flux_bounds: &str| -> Engine {
        let mut e = Engine::new();
        e.call("project_new", &json!({ "title": "AF4 Bounds", "dim": 3 })).unwrap();
        e.call(
            "add_polygons",
            &json!({
                "name": "g_poly",
                "polygons": [[[-150e-6, -150e-6], [150e-6, -150e-6], [150e-6, 150e-6], [-150e-6, 150e-6]]],
                "holes": [[[-110e-6, -80e-6], [110e-6, -80e-6], [110e-6, 80e-6], [-110e-6, 80e-6]]]
            }),
        )
        .unwrap();
        e.call(
            "add_polygons",
            &json!({
                "name": "s_poly",
                "polygons": [[[-90e-6, -55e-6], [-15e-6, -55e-6], [-15e-6, 55e-6], [-90e-6, 55e-6]]]
            }),
        )
        .unwrap();
        e.call(
            "add_polygons",
            &json!({
                "name": "b_poly",
                "polygons": [[[15e-6, -55e-6], [90e-6, -55e-6], [90e-6, 55e-6], [15e-6, 55e-6]]]
            }),
        )
        .unwrap();
        e.call("add_sheet", &json!({ "name": "c0", "input": "g_poly", "z": 0.0 })).unwrap();
        e.call("add_sheet", &json!({ "name": "c1", "input": "s_poly", "z": 0.0 })).unwrap();
        e.call("add_sheet", &json!({ "name": "c2", "input": "b_poly", "z": 0.0 })).unwrap();
        e.call(
            "add_chip_stack",
            &json!({
                "bbox": [-170e-6, -170e-6, 170e-6, 170e-6],
                "t_sub": 90e-6,
                "h_vac_top": 90e-6,
                "h_vac_bot": 0.0
            }),
        )
        .unwrap();
        for s in ["c0", "c1", "c2"] {
            e.call("add_size", &json!({ "sel": format!("{s}:edge"), "h": edge_h })).unwrap();
        }
        e.call("set_mesh", &json!({ "h": h })).unwrap();
        e.call(
            "add_physics",
            &json!({
                "physics": "dfs",
                "study": "electrostatic",
                "order": 2,
                "study_params": { "flux_bounds": flux_bounds }
            }),
        )
        .unwrap();
        e.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
        e.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
        e.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "ground" })).unwrap();
        for s in ["c0", "c1", "c2"] {
            e.call("add_condition", &json!({ "physics": "dfs", "on": s, "kind": "terminal", "params": { "name": s } }))
                .unwrap();
        }
        e.call("solve_run", &json!({})).unwrap();
        e
    };

    let e_rep = build_pocket(55e-6, 16e-6, "report");
    let mut e_mid = build_pocket(55e-6, 16e-6, "midpoint");
    let mut e_ref = build_pocket(24e-6, 5e-6, "midpoint");

    let cap_rep = CapacitanceMatrixData::of(e_rep.result("dfs").unwrap().0).unwrap();
    let cap_mid = CapacitanceMatrixData::of(e_mid.result("dfs").unwrap().0).unwrap();
    let cap_ref = CapacitanceMatrixData::of(e_ref.result("dfs").unwrap().0).unwrap();

    let c_low = cap_rep.maxwell_lower_matrix_f.as_ref().unwrap();
    let c_up = cap_rep.maxwell_upper_matrix_f.as_ref().unwrap();
    let c_mid = cap_rep.maxwell_mid_matrix_f.as_ref().unwrap();
    for i in 0..3 {
        let cr = cap_ref.maxwell_matrix_f[i][i];
        assert!(
            c_low[i][i] < cr && cr < c_up[i][i],
            "terminal {i} self-capacitance not bracketed: low={:e}, ref={cr:e}, up={:e}",
            c_low[i][i],
            c_up[i][i]
        );
        // `report` keeps `maxwell_matrix_f == c_up`; `midpoint` promotes `c_mid` everywhere
        assert!((cap_rep.maxwell_matrix_f[i][i] - c_up[i][i]).abs() <= 1e-14 * c_up[i][i]);
        assert!((cap_mid.maxwell_matrix_f[i][i] - c_mid[i][i]).abs() <= 1e-14 * c_mid[i][i]);
    }

    // Verify `Solution::value` and `C_*` / `E_C_*` scalars match `cap_mid` when `flux_bounds == "midpoint"`
    let r_mid = e_mid.result("dfs").unwrap().0;
    for (j, row) in c_mid.iter().enumerate().take(3) {
        assert!((r_mid.solutions[j].value.unwrap() - row[j]).abs() <= 1e-14 * row[j]);
    }
    assert!((r_mid.scalar("C_c1_c2 (F)").unwrap() - c_mid[1][2]).abs() <= 1e-14 * c_mid[1][2].abs());
    assert!(
        (r_mid.scalar("E_C_c1_c1 (Hz)").unwrap() - cap_mid.ec_matrix_hz[1][1]).abs()
            <= 1e-12 * cap_mid.ec_matrix_hz[1][1]
    );

    // Verify `analyse("capacitance_matrix")` brackets `E_C_Hz`
    let q_params = json!({
        "ground": "c0",
        "qubits": [{ "name": "Q", "pads": ["c1", "c2"], "E_J_Hz": 16e9 }]
    });
    let ana_mid = e_mid
        .call("analyse", &json!({ "analysis": "capacitance_matrix", "physics": "dfs", "params": q_params }))
        .unwrap();
    let ana_ref = e_ref
        .call("analyse", &json!({ "analysis": "capacitance_matrix", "physics": "dfs", "params": q_params }))
        .unwrap();
    let ec_low = ana_mid["qubits"][0]["E_C_low_Hz"].as_f64().unwrap();
    let ec_up = ana_mid["qubits"][0]["E_C_up_Hz"].as_f64().unwrap();
    let ec_ref = ana_ref["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    assert!(ec_low < ec_ref && ec_ref < ec_up, "qubit E_C not bracketed: [{ec_low:e}, {ec_up:e}] vs {ec_ref:e}");
}

#[test]
fn af5_gds_arc_fitting_curved_sheet3d_and_project_script_roundtrip() {
    let dir = temp_dir("af5");
    let gds_path = dir.join("annular.gds");

    // 64-gon circular pad of radius R = 40 um in nm integer GDS coordinates (dbu = 1 nm)
    let r_nm = 40_000.0_f64;
    let ring64: Vec<[i32; 2]> = (0..64)
        .map(|k| {
            let th = 2.0 * PI * (k as f64) / 64.0;
            [(r_nm * th.cos()).round() as i32, (r_nm * th.sin()).round() as i32]
        })
        .collect();
    let mut gb = GdsBuilder::new(1e-3, 1e-9);
    gb.begin_cell("TOP");
    gb.boundary(1, 0, &ring64);
    gb.text(10, 0, [0, 0], "disk_pad");
    gb.end_cell();
    std::fs::write(&gds_path, gb.finish()).unwrap();

    let mut e = Engine::new();
    e.call("project_new", &json!({ "title": "AF5 Curved Arc GDS", "dim": 3 })).unwrap();
    let imp = e
        .call(
            "import_gds",
            &json!({
                "path": gds_path.to_str().unwrap(),
                "layers": { "metal": [1, 0] },
                "pin_layer": [10, 0],
                "arc_tol": 200e-9,
                "lateral_pad": 60e-6
            }),
        )
        .unwrap();
    assert_eq!(imp["sheets"], json!(["disk_pad"]));
    assert!(imp["padded_bbox"][2].as_f64().unwrap() > 99e-6);

    e.call(
        "add_chip_stack",
        &json!({
            "bbox": imp["unpadded_bbox"],
            "lateral_pad": 60e-6,
            "t_sub": 80e-6,
            "h_vac_top": 80e-6,
            "h_vac_bot": 0.0
        }),
    )
    .unwrap();
    e.call(
        "set_mesh",
        &json!({
            "h": 24e-6,
            "curved": true,
            "gap_cells": 2.5,
            "adapt_grade": 0.38,
            "adapt_max_cells": 250000
        }),
    )
    .unwrap();
    e.call("mesh_build", &json!({})).unwrap();
    assert!(
        e.model().unwrap().solver_mesh().has_curved(),
        "fitted GDS arc sheet with curved=true must produce CurvedTet EdgeNodes"
    );

    // Verify `project_script` exports and cleanly verifies all new `Polygons.arc_tol` and `set_mesh` fields
    let script_out = e.call("project_script", &json!({ "form": "commands" })).unwrap();
    let script_txt = script_out["script"].as_str().unwrap();
    assert!(script_txt.contains("arc_tol="), "project_script must serialize Polygons.arc_tol: {script_txt}");
    assert!(script_txt.contains("curved=True"), "project_script must serialize curved=True: {script_txt}");
    assert!(script_txt.contains("adapt_grade="), "project_script must serialize adapt_grade: {script_txt}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn af6_conformal_thick_metal_trench_tls_and_13ghz_ac_dielectric_loss() {
    // 1. Thick metal (`150 nm`) + substrate over-etch trench (`80 nm`) on `silicon` at `13 GHz`
    let mut e = Engine::new();
    e.call("project_new", &json!({ "title": "AF6 Thick AC", "dim": 3 })).unwrap();
    e.call(
        "add_polygons",
        &json!({
            "name": "p1",
            "polygons": [[[-45e-6, -25e-6], [-8e-6, -25e-6], [-8e-6, 25e-6], [-45e-6, 25e-6]]]
        }),
    )
    .unwrap();
    e.call(
        "add_polygons",
        &json!({
            "name": "p2",
            "polygons": [[[8e-6, -25e-6], [45e-6, -25e-6], [45e-6, 25e-6], [8e-6, 25e-6]]]
        }),
    )
    .unwrap();
    e.call("add_sheet", &json!({ "name": "pad1", "input": "p1", "z": 0.0, "thickness": 150e-9, "etch_depth": 80e-9 }))
        .unwrap();
    e.call("add_sheet", &json!({ "name": "pad2", "input": "p2", "z": 0.0, "thickness": 150e-9, "etch_depth": 80e-9 }))
        .unwrap();
    e.call(
        "add_chip_stack",
        &json!({
            "bbox": [-80e-6, -60e-6, 80e-6, 60e-6],
            "t_sub": 60e-6,
            "h_vac_top": 60e-6,
            "h_vac_bot": 0.0
        }),
    )
    .unwrap();
    e.call("set_mesh", &json!({ "h": 22e-6 })).unwrap();
    e.call(
        "add_physics",
        &json!({
            "physics": "dfs",
            "study": "electrostatic",
            "order": 2,
            "study_params": { "freq": 13e9 }
        }),
    )
    .unwrap();
    // Shipped `silicon` TOML material must decode at `freq = 13 GHz` without optical range error!
    e.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": "silicon" })).unwrap();
    e.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "ground" })).unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "pad1", "kind": "terminal", "params": { "name": "c1" } }))
        .unwrap();
    e.call("add_condition", &json!({ "physics": "dfs", "on": "pad2", "kind": "terminal", "params": { "name": "c2" } }))
        .unwrap();
    e.call("solve_run", &json!({})).unwrap();

    let (res, _) = e.result("dfs").unwrap();
    let cap = CapacitanceMatrixData::of(res).unwrap();
    let g_mat = cap.conductance_matrix_s.as_ref().expect("13 GHz silicon must produce conductance_matrix_s");
    let q_mat = cap.q_dielectric_matrix.as_ref().expect("13 GHz silicon must produce q_dielectric_matrix");
    assert!(g_mat[0][0] > 0.0 && q_mat[0][0].unwrap() > 1.0e5);

    let tls_meta = TlsMeta::of(res).expect("TlsMeta present");
    assert!(!tls_meta.sheet_thickness_m.is_empty(), "TlsMeta must record thick sheet thickness");
    let tls_json = e.call("analyse", &json!({ "analysis": "tls_participation", "physics": "dfs" })).unwrap();
    assert!(tls_json["p_ma"].as_f64().unwrap() > 0.0);
    assert!(tls_json["p_ms"].as_f64().unwrap() > 0.0);
    assert_eq!(
        tls_json["p_ma_edge"].as_f64().unwrap(),
        0.0,
        "thick metal must skip knife-edge asymptotic double-counting"
    );
    assert_eq!(tls_json["p_ms_edge"].as_f64().unwrap(), 0.0);

    // 2. Strong-loss two-layer series capacitor (`tan_delta = 0.5`) through `Engine::call` proves positive `Delta C > 0`
    let solve_series = |freq_hz: f64| -> f64 {
        let mut s = Engine::new();
        s.call("project_new", &json!({ "title": "Series Loss", "dim": 2 })).unwrap();
        s.call("add_rect", &json!({ "name": "top_r", "min": [0.0, 0.0], "max": [10e-6, 5e-6] })).unwrap();
        s.call("add_rect", &json!({ "name": "bot_r", "min": [0.0, -5e-6], "max": [10e-6, 0.0] })).unwrap();
        s.call(
            "add_overlay",
            &json!({
                "regions": [
                    { "name": "Top", "input": "top_r", "kind": "Material" },
                    { "name": "Bot", "input": "bot_r", "kind": "Material" }
                ]
            }),
        )
        .unwrap();
        s.call("set_mesh", &json!({ "h": 2.5e-6 })).unwrap();
        s.call(
            "add_physics",
            &json!({
                "physics": "dfs",
                "study": "electrostatic",
                "order": 2,
                "study_params": { "freq": freq_hz }
            }),
        )
        .unwrap();
        s.call(
            "set_material",
            &json!({ "physics": "dfs", "region": "Top", "material": { "eps_r": 4.0, "tan_delta": 0.0 } }),
        )
        .unwrap();
        s.call(
            "set_material",
            &json!({ "physics": "dfs", "region": "Bot", "material": { "eps_r": 4.0, "tan_delta": 0.5 } }),
        )
        .unwrap();
        s.call(
            "add_condition",
            &json!({ "physics": "dfs", "on": "top_r.top", "kind": "terminal", "params": { "name": "t" } }),
        )
        .unwrap();
        s.call("add_condition", &json!({ "physics": "dfs", "on": "bot_r.bottom", "kind": "ground" })).unwrap();
        s.call("solve_run", &json!({})).unwrap();
        CapacitanceMatrixData::of(s.result("dfs").unwrap().0).unwrap().maxwell_matrix_f[0][0]
    };
    let c_dc = solve_series(0.0);
    let c_ac = solve_series(10e9);
    assert!(
        c_ac > c_dc,
        "strong dielectric loss must increase effective AC capacitance (c_ac={c_ac:e} > c_dc={c_dc:e})"
    );
}

#[test]
fn cross_feature_compatibility_rules_1_through_5() {
    // Rule 1: `set_mesh` rejects invalid parameter ranges with `Error::Invalid`
    let mut e1 = Engine::new();
    e1.call("project_new", &json!({ "title": "Rule1", "dim": 3 })).unwrap();
    assert!(matches!(e1.call("set_mesh", &json!({ "quarter_point_xi": 0.1 })), Err(Error::Invalid(_))));
    assert!(matches!(e1.call("set_mesh", &json!({ "singular_edges": "bogus" })), Err(Error::Invalid(_))));
    assert!(matches!(e1.call("set_mesh", &json!({ "adapt_grade": 0.0 })), Err(Error::Invalid(_))));
    assert!(matches!(e1.call("set_mesh", &json!({ "adapt_max_cells": 0 })), Err(Error::Invalid(_))));

    // Rule 2: `singular_edges="quarter_point"` auto-skips sheets with `thickness > 0`
    let mut e2 = Engine::new();
    e2.call("project_new", &json!({ "title": "Rule2", "dim": 3 })).unwrap();
    e2.call("add_circle", &json!({ "name": "c2d", "center": [0.0, 0.0], "r": 30e-6 })).unwrap();
    e2.call("add_sheet", &json!({ "name": "thick_pad", "input": "c2d", "z": 0.0, "thickness": 200e-9 })).unwrap();
    e2.call("add_chip_stack", &json!({ "lateral_pad": 50e-6, "t_sub": 60e-6, "h_vac_top": 60e-6, "h_vac_bot": 0.0 }))
        .unwrap();
    e2.call("set_mesh", &json!({ "h": 25e-6, "singular_edges": "quarter_point" })).unwrap();
    e2.call("mesh_build", &json!({})).unwrap();
    assert!(!e2.model().unwrap().solver_mesh().has_qp(), "Rule 2: thick slab must skip quarter-point shifting");

    // Rule 3: `quarter_point` + `flux_bounds != "none"` AND `rho_v != 0` + `flux_bounds != "none"` -> `Error::Invalid`
    e2.call(
        "add_physics",
        &json!({
            "physics": "dfs",
            "study": "electrostatic",
            "order": 2,
            "study_params": { "flux_bounds": "report" }
        }),
    )
    .unwrap();
    e2.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
    e2.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
    e2.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "ground" })).unwrap();
    e2.call(
        "add_condition",
        &json!({ "physics": "dfs", "on": "thick_pad", "kind": "terminal", "params": { "name": "p" } }),
    )
    .unwrap();
    let err_r3 = e2.call("solve_run", &json!({})).unwrap_err();
    assert_eq!(err_r3.code(), "invalid", "Rule 3 expected Error::Invalid, got {err_r3:?}");

    e2.call("set_mesh", &json!({ "singular_edges": null })).unwrap();
    e2.call(
        "set_material",
        &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45, "rho_v": 1e-3 } }),
    )
    .unwrap();
    let err_rho = e2.call("solve_run", &json!({})).unwrap_err();
    assert_eq!(err_rho.code(), "invalid", "source charge + flux_bounds expected Error::Invalid, got {err_rho:?}");

    // Rule 4: `curved=true` + `singular_edges="quarter_point"` composes in strict order without Jacobian inversion
    let mut e4 = Engine::new();
    e4.call("project_new", &json!({ "title": "Rule4", "dim": 3 })).unwrap();
    e4.call("add_circle", &json!({ "name": "c2d", "center": [0.0, 0.0], "r": 35e-6 })).unwrap();
    e4.call("add_sheet", &json!({ "name": "disk", "input": "c2d", "z": 0.0 })).unwrap();
    e4.call("add_chip_stack", &json!({ "lateral_pad": 55e-6, "t_sub": 60e-6, "h_vac_top": 60e-6, "h_vac_bot": 0.0 }))
        .unwrap();
    e4.call(
        "set_mesh",
        &json!({
            "h": 24e-6,
            "curved": true,
            "singular_edges": "quarter_point",
            "quarter_point_xi": 0.25
        }),
    )
    .unwrap();
    e4.call("add_physics", &json!({ "physics": "dfs", "study": "electrostatic", "order": 2 })).unwrap();
    e4.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
    e4.call("set_material", &json!({ "physics": "dfs", "region": "VacTop", "material": { "eps_r": 1.0 } })).unwrap();
    e4.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "ground" })).unwrap();
    e4.call("add_condition", &json!({ "physics": "dfs", "on": "disk", "kind": "terminal", "params": { "name": "d" } }))
        .unwrap();
    e4.call("solve_run", &json!({})).unwrap();
    assert!(e4.model().unwrap().solver_mesh().has_curved());
    assert!(e4.model().unwrap().solver_mesh().has_qp());

    // Rule 5: `planar_bem` supports `freq` + `tan_delta`, and rejects `thickness > 0`, `flux_bounds != "none"`, and `far_field` with `Error::Unsupported`
    e4.call(
        "set_study",
        &json!({
            "physics": "dfs",
            "kind": "electrostatic",
            "order": 1,
            "params": { "formulation": "planar_bem", "freq": 13e9 }
        }),
    )
    .unwrap();
    e4.call("set_mesh", &json!({ "curved": false, "singular_edges": null })).unwrap();
    e4.call(
        "set_material",
        &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45, "tan_delta": 1.0e-4 } }),
    )
    .unwrap();
    e4.call("solve_run", &json!({})).unwrap();
    let bem_cap = CapacitanceMatrixData::of(e4.result("dfs").unwrap().0).unwrap();
    assert!(bem_cap.conductance_matrix_s.is_some() && bem_cap.q_dielectric_matrix.is_some());

    // `planar_bem` + `flux_bounds` -> `Error::Unsupported`
    e4.call(
        "set_study",
        &json!({
            "physics": "dfs",
            "kind": "electrostatic",
            "order": 1,
            "params": { "formulation": "planar_bem", "flux_bounds": "report" }
        }),
    )
    .unwrap();
    assert_eq!(e4.call("solve_run", &json!({})).unwrap_err().code(), "unsupported");

    // `planar_bem` + `far_field` -> `Error::Unsupported`
    e4.call(
        "set_study",
        &json!({
            "physics": "dfs",
            "kind": "electrostatic",
            "order": 1,
            "params": { "formulation": "planar_bem" }
        }),
    )
    .unwrap();
    e4.call("add_condition", &json!({ "physics": "dfs", "on": "outer", "kind": "far_field" })).unwrap();
    assert_eq!(e4.call("solve_run", &json!({})).unwrap_err().code(), "unsupported");

    // `planar_bem` + `thickness > 0` -> `Error::Unsupported`
    e2.call(
        "set_study",
        &json!({
            "physics": "dfs",
            "kind": "electrostatic",
            "order": 1,
            "params": { "formulation": "planar_bem" }
        }),
    )
    .unwrap();
    e2.call("set_material", &json!({ "physics": "dfs", "region": "Si", "material": { "eps_r": 11.45 } })).unwrap();
    assert_eq!(e2.call("solve_run", &json!({})).unwrap_err().code(), "unsupported");
}
