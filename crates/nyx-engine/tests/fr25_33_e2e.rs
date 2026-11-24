// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end integration tests through `Engine::call` covering SQuADDS Wave 3 / Post-SQuADDS
//! bug fixes and solver/meshing features (`FR-25`..`FR-33`, ADRs 125–127):
//! - `FR-30`: `ibc_mode: "perturbative"` and `"auto"` on `ems::cavity` (`surface_impedance`)
//!   shifting modal frequency downward via Wheeler kinetic inductance and yielding BCS `Q`
//!   without complex non-Hermitian Krylov assembly.
//! - `FR-32`: Hybrid quasi-electrostatic modal TLS participation (`analyse("tls_participation")`
//!   with `electrostatic_physics: "dfs"` and explicit `modal_voltages`), combining a fine
//!   quarter-point `P2` electrostatic solve with a coarse `P1` full-chip `ems::cavity` ROM,
//!   alongside coarse-`P1` `edge_singularity_underresolved` diagnostics.
//! - `FR-33`: `ems::driven` memory-scalable linear solver backends (`solver: "blr"` with
//!   `blr_tol` and symbolic reuse, `solver: "ams_fgmres"` at `order: 1` and `order: 2`) and
//!   `(solver_backend, order)`-keyed engine memory calibration (`estimate::chain_key_for_project`).

use nyx_engine::estimate::{MemoryProbe, Row};
use nyx_engine::Engine;
use nyx_physics::DrivenPortData;
use serde_json::{json, Value};
use std::sync::Arc;

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}): {err}"))
}

fn build_3d_transmon_cavity_project(e: &mut Engine, ibc_mode: Option<&str>) {
    call(e, "project_new", json!({ "title": "FR-30 Perturbative IBC Cavity", "dim": 3, "h": 240.0e-6 }));
    let (w_box, l_box, h_sub, h_air) = (1.4e-3, 1.4e-3, 300.0e-6, 300.0e-6);
    call(
        e,
        "add_box",
        json!({ "name": "si_box", "min": [-w_box / 2.0, -l_box / 2.0, -h_sub], "max": [w_box / 2.0, l_box / 2.0, 0.0] }),
    );
    call(
        e,
        "add_box",
        json!({ "name": "air_box", "min": [-w_box / 2.0, -l_box / 2.0, 0.0], "max": [w_box / 2.0, l_box / 2.0, h_air] }),
    );

    call(e, "add_rect", json!({ "name": "pad1_2d", "min": [-200.0e-6, -500.0e-6], "max": [200.0e-6, -60.0e-6] }));
    call(e, "add_sheet", json!({ "name": "res_pad1", "input": "pad1_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "pad2_2d", "min": [-200.0e-6, 60.0e-6], "max": [200.0e-6, 500.0e-6] }));
    call(e, "add_sheet", json!({ "name": "res_pad2", "input": "pad2_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "jj_2d", "min": [-60.0e-6, -60.0e-6], "max": [60.0e-6, 60.0e-6] }));
    call(e, "add_sheet", json!({ "name": "jj_sheet", "input": "jj_2d", "z": 0.0 }));

    call(
        e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "si", "input": "si_box" },
                { "name": "air", "input": "air_box" }
            ]
        }),
    );
    call(e, "set_mesh", json!({ "h": 240.0e-6 }));
    call(
        e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 1,
            "study_params": { "n_modes": 1, "f_search": 5.0e9 }
        }),
    );
    call(e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    call(e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));

    for sheet in ["res_pad1", "res_pad2"] {
        if let Some(mode) = ibc_mode {
            call(
                e,
                "add_condition",
                json!({
                    "physics": "ems",
                    "on": sheet,
                    "kind": "surface_impedance",
                    "params": {
                        "material": "aluminium_sc_emft",
                        "temperature": 0.45,
                        "thickness": 80.0e-9,
                        "ibc_mode": mode
                    }
                }),
            );
        } else {
            call(e, "add_condition", json!({ "physics": "ems", "on": sheet, "kind": "pec", "params": {} }));
        }
    }
    call(
        e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "jj_sheet",
            "kind": "lumped_jj",
            "params": {
                "name": "Q1_JJ",
                "E_j_Hz": 20.0e9,
                "C_j": 35.0e-15,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }
        }),
    );
}

#[test]
fn fr30_perturbative_bcs_ibc_cavity_through_engine_call() {
    let mut e_pec = Engine::new();
    build_3d_transmon_cavity_project(&mut e_pec, None);
    let r_pec = call(&mut e_pec, "solve_run", json!({ "isolate": false }));
    assert_eq!(r_pec["state"], json!("done"));
    let (res_pec, _) = e_pec.result("ems").unwrap();
    let f_pec = res_pec.scalar("f re 1 (Hz)").expect("PEC f re 1");

    let mut e_pert = Engine::new();
    build_3d_transmon_cavity_project(&mut e_pert, Some("perturbative"));
    let r_pert = call(&mut e_pert, "solve_run", json!({ "isolate": false }));
    assert_eq!(r_pert["state"], json!("done"));
    let (res_pert, _) = e_pert.result("ems").unwrap();
    let f_pert = res_pert.scalar("f re 1 (Hz)").expect("perturbative f re 1");
    let f_im_pert = res_pert.scalar("f im 1 (Hz)").expect("perturbative f im 1");
    let q_pert = res_pert.scalar("Q1").expect("perturbative Q1");

    assert!(
        f_pert < f_pec,
        "perturbative BCS kinetic inductance must shift mode downward: f_pert = {f_pert:.6e} Hz vs f_pec = {f_pec:.6e} Hz"
    );
    assert!(
        f_im_pert > 0.0 && q_pert > 1.0e2 && q_pert.is_finite(),
        "expected positive finite BCS Q1, got Q1={q_pert}, f_im={f_im_pert}"
    );
    assert!(
        res_pert.log.iter().any(|l| l.contains("surface_impedance perturbative mode active")),
        "expected perturbative IBC log line, got {:?}",
        res_pert.log
    );

    // CircuitRom & downstream cQED analyses remain valid on perturbative IBC cavity results
    let ham = call(&mut e_pert, "analyse", json!({ "analysis": "hamiltonian", "physics": "ems" }));
    let f_q = ham["f01_Hz"][0].as_f64().unwrap_or(0.0);
    assert!(
        f_q > 1.0e9 && f_q.is_finite(),
        "expected physical transmon f01_Hz from perturbative-IBC CircuitRom: {ham}"
    );
}

#[test]
fn fr32_hybrid_quasi_electrostatic_modal_tls_through_engine_call() {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "FR-32 Hybrid Modal TLS", "dim": 3, "h": 130.0e-6 }));
    let (w_box, l_box, h_sub, h_air) = (1.1e-3, 1.1e-3, 250.0e-6, 250.0e-6);
    call(
        &mut e,
        "add_box",
        json!({ "name": "si_box", "min": [-w_box / 2.0, -l_box / 2.0, -h_sub], "max": [w_box / 2.0, l_box / 2.0, 0.0] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "air_box", "min": [-w_box / 2.0, -l_box / 2.0, 0.0], "max": [w_box / 2.0, l_box / 2.0, h_air] }),
    );
    call(&mut e, "add_rect", json!({ "name": "pad1_2d", "min": [-180.0e-6, -380.0e-6], "max": [180.0e-6, -45.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "res_pad1", "input": "pad1_2d", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "pad2_2d", "min": [-180.0e-6, 45.0e-6], "max": [180.0e-6, 380.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "res_pad2", "input": "pad2_2d", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "jj_2d", "min": [-40.0e-6, -45.0e-6], "max": [40.0e-6, 45.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "jj_sheet", "input": "jj_2d", "z": 0.0 }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "si", "input": "si_box" },
                { "name": "air", "input": "air_box" }
            ]
        }),
    );

    // Stage 1: Fine order-2 quarter-point singular-edge electrostatic solve on `dfs`
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": 130.0e-6,
            "gap_cells": 3.0,
            "singular_edges": "quarter_point"
        }),
    );
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "air", "material": { "eps_r": 1.0 } }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "res_pad1", "kind": "terminal", "params": { "name": "pad1" } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "res_pad2", "kind": "terminal", "params": { "name": "pad2" } }),
    );
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "si_box.zmin", "kind": "ground", "params": {} }));
    call(&mut e, "solve_run", json!({ "physics": "dfs", "isolate": false }));
    let dfs_cells = e.model().unwrap().mesh().n_cells();

    // Stage 2: Switch to a coarse straight P1 mesh for full-chip `ems::cavity`
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": 280.0e-6,
            "gap_cells": null,
            "singular_edges": "none"
        }),
    );
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "cavity",
            "order": 1,
            "study_params": { "n_modes": 1, "f_search": 5.0e9 }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
    for sheet in ["res_pad1", "res_pad2"] {
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": sheet, "kind": "pec", "params": {} }));
    }
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "jj_sheet",
            "kind": "lumped_jj",
            "params": {
                "name": "Q1_JJ",
                "E_j_Hz": 20.0e9,
                "C_j": 35.0e-15,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }
        }),
    );
    call(&mut e, "solve_run", json!({ "physics": "ems", "isolate": false }));
    let ems_cells = e.model().unwrap().mesh().n_cells();
    assert!(
        dfs_cells > ems_cells,
        "fine electrostatic mesh ({dfs_cells}) must exceed coarse cavity mesh ({ems_cells})"
    );

    // Direct coarse-P1 TLS flags underresolved edge singularity with an actionable warning
    let tls_coarse = call(&mut e, "analyse", json!({ "analysis": "tls_participation", "physics": "ems" }));
    assert_eq!(tls_coarse["edge_singularity_underresolved"], json!(true), "{tls_coarse}");
    let warn = tls_coarse["warning"].as_str().unwrap_or("");
    assert!(!warn.is_empty(), "coarse P1 3D sheet TLS must emit a warning: {tls_coarse}");

    // Hybrid quasi-electrostatic modal TLS evaluates on dfs's stored P2 quarter-point mesh + ems's CircuitRom
    let tls_hybrid = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "tls_participation",
            "physics": "ems",
            "params": {
                "electrostatic_physics": "dfs",
                "port_terminals": { "Q1_JJ": ["pad1", "pad2"] },
                "d": 200.0e-9
            }
        }),
    );
    assert_eq!(tls_hybrid["quasi_electrostatic"], json!(true), "{tls_hybrid}");
    assert_eq!(tls_hybrid["electrostatic_physics"], json!("dfs"), "{tls_hybrid}");
    assert_eq!(tls_hybrid["edge_singularity_underresolved"], json!(false), "{tls_hybrid}");
    let p_tot = tls_hybrid["p_total"].as_f64().unwrap();
    let q_tls = tls_hybrid["q_tls"].as_f64().unwrap();
    let t1_tls = tls_hybrid["T1_TLS_s"].as_f64().unwrap();
    assert!(p_tot > 1.0e-6 && p_tot < 1.0e-1, "expected physical surface participation p_total, got {p_tot:.3e}");
    assert!(q_tls > 1.0e3 && q_tls.is_finite(), "expected finite q_tls, got {q_tls:.3e}");
    assert!(t1_tls > 1.0e-8 && t1_tls.is_finite(), "expected finite T1_TLS_s, got {t1_tls:.3e}");
    assert!(!tls_hybrid["modes"].as_array().unwrap().is_empty());

    // Explicit modal_voltages + modal_freqs_hz directly on `dfs` also works cleanly:
    // symmetric [-0.5, 0.5] agrees within < 1e-4 (mesh triangulation asymmetry), and exact
    // `terminal_voltages_V` matches `tls_hybrid` to machine precision.
    let tls_explicit = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "tls_participation",
            "physics": "dfs",
            "params": {
                "modal_voltages": [[-0.5, 0.5], tls_hybrid["terminal_voltages_V"]],
                "modal_freqs_hz": [5.0e9, tls_hybrid["freq_Hz"]],
                "d": 200.0e-9
            }
        }),
    );
    let p_exp = tls_explicit["p_total"].as_f64().unwrap();
    assert!(
        (p_exp - p_tot).abs() / p_tot < 1.0e-4,
        "symmetric differential [-0.5, 0.5] must agree with charge-balanced 2-pad hybrid TLS: {p_exp:.6e} vs {p_tot:.6e}"
    );
    let p_exact = tls_explicit["modes"][1]["p_total"].as_f64().unwrap();
    assert!(
        (p_exact - p_tot).abs() <= 1.0e-14 * p_tot,
        "explicit terminal_voltages_V must match CircuitRom-mapped hybrid TLS to machine precision: {p_exact:.16e} vs {p_tot:.16e}"
    );
}

#[test]
fn fr33_driven_blr_ams_fgmres_and_backend_memory_calibration_through_engine_call() {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "FR-33 Driven Solvers & Calibration", "dim": 3, "h": 240.0e-6 }));
    let (w_box, l_box, h_sub, h_air) = (1.4e-3, 1.4e-3, 300.0e-6, 300.0e-6);
    call(
        &mut e,
        "add_box",
        json!({ "name": "si_box", "min": [-w_box / 2.0, -l_box / 2.0, -h_sub], "max": [w_box / 2.0, l_box / 2.0, 0.0] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "air_box", "min": [-w_box / 2.0, -l_box / 2.0, 0.0], "max": [w_box / 2.0, l_box / 2.0, h_air] }),
    );
    call(&mut e, "add_rect", json!({ "name": "p1_2d", "min": [-520.0e-6, -700.0e-6], "max": [-340.0e-6, -460.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "port1_sheet", "input": "p1_2d", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "feed_2d", "min": [-520.0e-6, -460.0e-6], "max": [-340.0e-6, 460.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "feed_trace", "input": "feed_2d", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "p2_2d", "min": [-520.0e-6, 460.0e-6], "max": [-340.0e-6, 700.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "port2_sheet", "input": "p2_2d", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "pad1_2d", "min": [-140.0e-6, -500.0e-6], "max": [420.0e-6, -60.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "res_pad1", "input": "pad1_2d", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "pad2_2d", "min": [-140.0e-6, 60.0e-6], "max": [420.0e-6, 500.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "res_pad2", "input": "pad2_2d", "z": 0.0 }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "si", "input": "si_box" },
                { "name": "air", "input": "air_box" }
            ]
        }),
    );
    call(&mut e, "set_mesh", json!({ "h": 240.0e-6 }));
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "driven",
            "order": 1,
            "study_params": {
                "freqs_hz": [5.0e9, 6.0e9],
                "prom_tol": 1.0e-3,
                "solver": "auto"
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));
    for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
        call(&mut e, "add_condition", json!({ "physics": "ems", "on": sheet, "kind": "pec", "params": {} }));
    }
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port1_sheet",
            "kind": "lumped_port",
            "params": { "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port2_sheet",
            "kind": "lumped_port",
            "params": { "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }
        }),
    );

    // Baseline solve with `solver: "auto"`
    call(&mut e, "solve_run", json!({ "isolate": false }));
    let d_auto = DrivenPortData::of(e.result("ems").unwrap().0).unwrap();
    let cells = e.model().unwrap().mesh().n_cells() as u64;

    // Seed calibration history with 32 GB total RAM and >80 GB rows on legacy `"ems"` and `"ems:driven:direct:p2"`
    let probe: MemoryProbe = Arc::new(|| Some((32_000_000_000, 32_000_000_000)));
    e.set_memory_probe(probe);
    for chain in ["ems", "ems:driven:direct:p2"] {
        for (c, gb) in [(cells / 2, 80.0), (cells * 2, 180.0)] {
            e.history().add(Row {
                chain: chain.into(),
                dim: 3,
                cells: c,
                dofs: c * 6,
                peak_working_set: (gb * 1.0e9) as u64,
                peak_commit: (gb * 1.2e9) as u64,
                seconds: 25.0,
            });
        }
    }

    // `solver: "direct"`, `order: 2` matches `"ems:driven:direct:p2"` and is refused by REFUSE_ABOVE
    call(
        &mut e,
        "set_study",
        json!({
            "physics": "ems",
            "kind": "driven",
            "order": 2,
            "params": { "freqs_hz": [5.0e9, 6.0e9], "prom_tol": null, "solver": "direct" }
        }),
    );
    let err_direct = e.call("solve_run", &json!({ "isolate": false })).unwrap_err().to_string();
    assert!(
        err_direct.contains("more than this machine can hold"),
        "expected REFUSE_ABOVE on ems:driven:direct:p2, got: {err_direct}"
    );

    // `solver: "blr"` (order 1) and `solver: "ams_fgmres"` (order 1 & order 2) are NOT falsely blocked!
    for (solver_name, order, prom_tol, expected_backend) in [
        ("blr", 1, json!(1.0e-3), "backend=BlrNd"),
        ("ams_fgmres", 1, json!(null), "backend=AmsFgmres"),
        ("ams_fgmres", 2, json!(null), "backend=AmsFgmres"),
    ] {
        call(
            &mut e,
            "set_study",
            json!({
                "physics": "ems",
                "kind": "driven",
                "order": order,
                "params": {
                    "freqs_hz": [5.0e9, 6.0e9],
                    "prom_tol": prom_tol,
                    "solver": solver_name,
                    "blr_tol": 1.0e-6
                }
            }),
        );
        let reply = call(&mut e, "solve_run", json!({ "isolate": false }));
        assert_eq!(reply["state"], json!("done"), "{solver_name} order={order} must be admitted cleanly");
        let (res, _) = e.result("ems").unwrap();
        assert_eq!(res.scalar("symbolic_count"), Some(1.0), "{solver_name} order={order} symbolic reuse");
        assert!(
            res.log.iter().any(|l| l.contains(&format!("solver={solver_name}")) && l.contains(expected_backend)),
            "{solver_name} order={order} log: {:?}",
            res.log
        );
        if order == 1 {
            let d = DrivenPortData::of(res).unwrap();
            let mut max_diff = 0.0_f64;
            for k in 0..d_auto.freqs_hz.len() {
                for p in 0..2 {
                    for q in 0..2 {
                        let dr = d.s_matrix[k][p][q][0] - d_auto.s_matrix[k][p][q][0];
                        let di = d.s_matrix[k][p][q][1] - d_auto.s_matrix[k][p][q][1];
                        max_diff = max_diff.max(dr.hypot(di));
                    }
                }
            }
            assert!(max_diff < 1.0e-3, "{solver_name} S-matrix diff vs auto = {max_diff:.3e}");
        }
    }

    // Touchstone network analysis works on the final order-2 ams_fgmres driven result
    let net = call(&mut e, "analyse", json!({ "analysis": "network", "physics": "ems" }));
    assert!(net["touchstone_snp"].as_str().unwrap_or("").contains("# Hz S RI R 50"), "{net}");
}
