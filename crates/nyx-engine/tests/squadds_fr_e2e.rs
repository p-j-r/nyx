// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end integration tests through `Engine::call` covering SQuADDS Wave 1 & Wave 2 functional requirements:
//! - FR-01: `save_fields: false`, in-worker `post_analyses`, `< 256 KiB` `results.h5`, and per-mode HDF5 budget.
//! - FR-03 + FR-12 + FR-04: 3D cavity with `lumped_port` + `lumped_jj` -> `CircuitRom` -> `analyse("hamiltonian")`
//!   (`f_cav_Hz`, `g_matrix_Hz`, `kappa_cav_Hz`, `t1_purcell_s`) and `analyse("network")` (QEP, BBQ, Taubin, `.sNp`).
//! - FR-15 + FR-18: DC-SQUID `squid_flux_sweep`, `ej_hz_sweep`, and multi-port `sweep_port`.
//! - FR-19 + FR-11 Pt 1 + FR-17: 3D `dfs::electrostatic` with `kind: "open"`, `capacitance_matrix` with
//!   `metal_thickness_m: 200e-9`, and `< 1 pF` `convergence_check` with `adapt_scale` and `quantities`.
//! - FR-05 Pt 1 + FR-07 + FR-10 Pt 1 + FR-21 + FR-23: `sheet3d` `max_cells`/`dz_min`/`max_z_slabs`/`z_growth_ratio`,
//!   Ruppert `max_steiner_multiplier`/`strict_mesh`, `memory_wait_timeout_s`, `lossless: true` BCS IBC, and TOML-driven
//!   `tls_participation` with `trench_depth_m`.

use hdf5_pure::File;
use nyx_engine::estimate::{MemoryProbe, Row};
use nyx_engine::h5io::ReadBudgetConfig;
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

fn worker() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"))
}

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}): {err}"))
}

fn temp_dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-squadds-fr-e2e-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Builds a 3D sheet-stack coupled feedline + transmon cavity project on `Engine`.
fn build_3d_cavity_chip(e: &mut Engine, second_jj: bool, lossless_ibc: bool) {
    call(e, "project_new", json!({ "title": "SQuADDS Cavity Chip E2E", "dim": 3, "h": 240.0e-6 }));
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

    // Feedline + two 50-Ohm lumped ports
    call(e, "add_rect", json!({ "name": "p1_2d", "min": [-520.0e-6, -700.0e-6], "max": [-340.0e-6, -460.0e-6] }));
    call(e, "add_sheet", json!({ "name": "port1_sheet", "input": "p1_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "feed_2d", "min": [-520.0e-6, -460.0e-6], "max": [-340.0e-6, 460.0e-6] }));
    call(e, "add_sheet", json!({ "name": "feed_trace", "input": "feed_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "p2_2d", "min": [-520.0e-6, 460.0e-6], "max": [-340.0e-6, 700.0e-6] }));
    call(e, "add_sheet", json!({ "name": "port2_sheet", "input": "p2_2d", "z": 0.0 }));

    // Primary qubit/resonator pads + JJ bridge
    call(e, "add_rect", json!({ "name": "pad1_2d", "min": [-140.0e-6, -500.0e-6], "max": [160.0e-6, -60.0e-6] }));
    call(e, "add_sheet", json!({ "name": "res_pad1", "input": "pad1_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "pad2_2d", "min": [-140.0e-6, 60.0e-6], "max": [160.0e-6, 500.0e-6] }));
    call(e, "add_sheet", json!({ "name": "res_pad2", "input": "pad2_2d", "z": 0.0 }));
    call(e, "add_rect", json!({ "name": "jj_2d", "min": [-40.0e-6, -60.0e-6], "max": [80.0e-6, 60.0e-6] }));
    call(e, "add_sheet", json!({ "name": "jj_sheet", "input": "jj_2d", "z": 0.0 }));

    if second_jj {
        call(e, "add_rect", json!({ "name": "pad3_2d", "min": [260.0e-6, -500.0e-6], "max": [560.0e-6, -60.0e-6] }));
        call(e, "add_sheet", json!({ "name": "sq_pad1", "input": "pad3_2d", "z": 0.0 }));
        call(e, "add_rect", json!({ "name": "pad4_2d", "min": [260.0e-6, 60.0e-6], "max": [560.0e-6, 500.0e-6] }));
        call(e, "add_sheet", json!({ "name": "sq_pad2", "input": "pad4_2d", "z": 0.0 }));
        call(e, "add_rect", json!({ "name": "jj2_2d", "min": [340.0e-6, -60.0e-6], "max": [480.0e-6, 60.0e-6] }));
        call(e, "add_sheet", json!({ "name": "squid_sheet", "input": "jj2_2d", "z": 0.0 }));
    }

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
            "study_params": { "n_modes": 2, "f_search": 5.0e9 }
        }),
    );
    call(e, "set_material", json!({ "physics": "ems", "region": "si", "material": { "eps_r": 11.45 } }));
    call(e, "set_material", json!({ "physics": "ems", "region": "air", "material": { "eps_r": 1.0 } }));

    for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
        if lossless_ibc {
            call(
                e,
                "add_condition",
                json!({
                    "physics": "ems",
                    "on": sheet,
                    "kind": "surface_impedance",
                    "params": {
                        "material": "aluminium_sc_emft",
                        "temperature": 0.02,
                        "lossless": true
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
    if second_jj {
        for sheet in ["sq_pad1", "sq_pad2"] {
            call(e, "add_condition", json!({ "physics": "ems", "on": sheet, "kind": "pec", "params": {} }));
        }
        call(
            e,
            "add_condition",
            json!({
                "physics": "ems",
                "on": "squid_sheet",
                "kind": "lumped_jj",
                "params": {
                    "name": "SQUID_JJ",
                    "E_j_Hz": 18.0e9,
                    "C_j": 30.0e-15,
                    "squid_asymmetry": 0.25,
                    "direction": [0.0, 1.0, 0.0],
                    "formulation": "rank1"
                }
            }),
        );
    }

    call(
        e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port1_sheet",
            "kind": "lumped_port",
            "params": { "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }
        }),
    );
    call(
        e,
        "add_condition",
        json!({
            "physics": "ems",
            "on": "port2_sheet",
            "kind": "lumped_port",
            "params": { "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }
        }),
    );
}

#[test]
fn fr01_save_fields_false_post_analyses_isolated_worker_and_per_mode_h5_budget() {
    let ws = temp_dir("fr01");
    let mut e = Engine::new();
    e.set_worker(worker(), true);
    e.set_workspace(ws.clone(), true);

    // Compact 3D electrostatic model so worker solve is fast and deterministic
    call(&mut e, "project_new", json!({ "title": "FR-01 Worker E2E", "dim": 3, "h": 140.0e-6 }));
    call(
        &mut e,
        "add_box",
        json!({ "name": "sub", "min": [-400.0e-6, -400.0e-6, -200.0e-6], "max": [400.0e-6, 400.0e-6, 0.0] }),
    );
    call(
        &mut e,
        "add_box",
        json!({ "name": "air", "min": [-400.0e-6, -400.0e-6, 0.0], "max": [400.0e-6, 400.0e-6, 200.0e-6] }),
    );
    call(&mut e, "add_rect", json!({ "name": "p1", "min": [-180.0e-6, -120.0e-6], "max": [-30.0e-6, 120.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "pad1", "input": "p1", "z": 0.0 }));
    call(&mut e, "add_rect", json!({ "name": "p2", "min": [30.0e-6, -120.0e-6], "max": [180.0e-6, 120.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "pad2", "input": "p2", "z": 0.0 }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "sub", "input": "sub" }, { "name": "air", "input": "air" }] }),
    );
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "sub", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "air", "material": { "eps_r": 1.0 } }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "pad1", "kind": "terminal", "params": { "name": "pad1" } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "pad2", "kind": "terminal", "params": { "name": "pad2" } }),
    );
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "sub.zmin", "kind": "ground", "params": {} }));

    let reply = call(
        &mut e,
        "solve_run",
        json!({
            "isolate": true,
            "save_fields": false,
            "post_analyses": [
                { "label": "cap", "analysis": "capacitance_matrix" }
            ]
        }),
    );
    assert_eq!(reply["fields_saved"], json!(false), "{reply}");
    assert!(reply["worker"].is_object(), "{reply}");
    let pa = reply["post_analyses"].as_array().expect("post_analyses present");
    assert_eq!(pa.len(), 1);
    assert_eq!(pa[0]["label"], json!("cap"));
    let c00 = pa[0]["maxwell_matrix_F"][0][0].as_f64().unwrap();
    assert!(c00 > 1.0e-16 && c00.is_finite(), "expected finite C00, got {c00}");

    // Stored result stripped fields while run results.h5 stays < 256 KiB and skips results.xdmf
    let (res, _) = e.result("dfs").unwrap();
    assert!(res.solutions.iter().all(|s| s.fields.is_empty() && s.cell_fields.is_empty()));
    let runs = call(&mut e, "store_runs", json!({}));
    let artifacts = runs[0]["artifacts"].as_array().unwrap();
    assert!(!artifacts.iter().any(|a| a["path"].as_str().unwrap().ends_with("results.xdmf")));
    let h5_bytes = artifacts
        .iter()
        .find(|a| a["path"].as_str().unwrap().ends_with("results.h5"))
        .expect("results.h5")
        .get("bytes")
        .and_then(|v| v.as_u64())
        .unwrap();
    assert!(h5_bytes > 0 && h5_bytes < 256 * 1024, "results.h5 must be < 256 KiB, got {h5_bytes}");

    // Per-mode HDF5 budget check on a full 2-solution export
    let full_h5 = ws.join("full_two_modes.h5");
    call(&mut e, "solve_run", json!({ "isolate": false, "save_fields": true }));
    call(&mut e, "export_results", json!({ "path": full_h5.to_str().unwrap() }));

    let raw = File::open_streaming(&full_h5).unwrap();
    let ds = raw.group("results/dfs").unwrap().dataset("cell_fields_packed").unwrap();
    let sh = ds.shape().unwrap();
    assert_eq!(sh.len(), 5);
    let m_modes = sh[0];
    assert_eq!(m_modes, 2);
    let per_mode_bytes = sh[1..].iter().product::<u64>() * (std::mem::size_of::<f64>() as u64);
    let total_bytes = m_modes * per_mode_bytes;

    // Budget between per_mode_bytes and total_bytes: per-mode budget MUST ALLOW reading!
    let allow_budget = ReadBudgetConfig {
        max_dataset_bytes: per_mode_bytes + 64,
        max_total_result_bytes: total_bytes * 4,
        ..ReadBudgetConfig::default()
    };
    let stored = nyx_engine::h5io::read_results_with_budget(&full_h5, allow_budget).expect("per-mode budget allows");
    assert!(stored.fields_saved);
    assert_eq!(stored.results[0].1.solutions.len(), 2);

    // Budget strictly below per_mode_bytes MUST reject with "exceeds read budget"
    let deny_budget = ReadBudgetConfig {
        max_dataset_bytes: per_mode_bytes.saturating_sub(8).max(8),
        max_total_result_bytes: total_bytes * 4,
        ..ReadBudgetConfig::default()
    };
    let err = nyx_engine::h5io::read_results_with_budget(&full_h5, deny_budget).err().expect("oversized mode rejected");
    assert!(err.contains("exceeds read budget"), "{err}");

    let _ = std::fs::remove_dir_all(&ws);
}

#[test]
fn fr03_fr12_fr04_cavity_lumped_port_jj_hamiltonian_and_network_e2e() {
    let mut e = Engine::new();
    build_3d_cavity_chip(&mut e, false, false);

    call(&mut e, "solve_run", json!({ "isolate": false }));

    // list_analyses must expose `"network"` and mark it applicable to `"ems"`
    let list = call(&mut e, "list_analyses", json!({}));
    let net_row = list
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row["analysis"] == json!("network"))
        .expect("network analysis registered in list_analyses");
    assert_eq!(net_row["applies_to"], json!(["ems"]), "{net_row}");

    // FR-03 + FR-12 + FR-22: hamiltonian on mixed lumped_jj + lumped_port CircuitRom
    let ham = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "hamiltonian",
            "physics": "ems",
            "params": { "basis": "craig_bampton" }
        }),
    );
    let f_cav = ham["f_cav_Hz"].as_array().expect("f_cav_Hz array");
    let g_mat = ham["g_matrix_Hz"].as_array().expect("g_matrix_Hz array");
    let kappa = ham["kappa_cav_Hz"].as_array().expect("kappa_cav_Hz array");
    let t1_p = ham["t1_purcell_s"].as_array().expect("t1_purcell_s array");
    assert!(!f_cav.is_empty(), "expected non-empty f_cav_Hz: {ham}");
    assert!(!g_mat.is_empty(), "expected non-empty g_matrix_Hz: {ham}");
    assert!(!kappa.is_empty(), "expected non-empty kappa_cav_Hz: {ham}");
    assert!(!t1_p.is_empty(), "expected non-empty t1_purcell_s: {ham}");
    assert!(f_cav[0].as_f64().unwrap() > 1.0e8, "{ham}");
    assert!(kappa[0].as_f64().unwrap() > 0.0, "{ham}");
    assert!(t1_p[0].as_f64().unwrap() > 0.0, "{ham}");

    // FR-04: unified network analysis (QEP + BBQ + Taubin S21 circle fit + Touchstone .sNp)
    let net = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "network",
            "physics": "ems",
            "params": {
                "s21_ports": ["port1", "port2"],
                "resonator_branch": "hanger",
                "touchstone_format": "RI"
            }
        }),
    );
    assert_eq!(net["analysis"], json!("network"));
    let qep = net["qep_modes"].as_array().expect("qep_modes");
    assert!(!qep.is_empty(), "expected non-empty qep_modes: {net}");
    assert!(qep.iter().all(|m| m["f_re_Hz"].as_f64().unwrap() > 0.0));
    let bbq = net["bbq_modes"].as_array().expect("bbq_modes");
    assert!(!bbq.is_empty(), "expected non-empty bbq_modes: {net}");
    assert!(bbq[0]["C_Sigma_F"].as_f64().unwrap() > 0.0, "{net}");
    assert!(bbq[0]["t1_purcell_s"].as_f64().unwrap() > 0.0, "{net}");
    assert!(net["s21_fit"].is_object(), "expected s21_fit object: {net}");
    assert!(net["s21_fit"]["f_r_Hz"].as_f64().unwrap() > 0.0, "{net}");
    let snp = net["touchstone_snp"].as_str().expect("touchstone_snp string");
    assert!(snp.contains("# Hz S RI R"), "unexpected Touchstone header: {snp}");
}

#[test]
fn fr15_fr18_dc_squid_flux_sweep_and_multi_port_sweep_port_e2e() {
    let mut e = Engine::new();
    build_3d_cavity_chip(&mut e, true, false);
    call(&mut e, "solve_run", json!({ "isolate": false }));

    // FR-15 + FR-18: sweep_port by name ("SQUID_JJ") and native squid_flux_sweep
    let by_name = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "hamiltonian",
            "physics": "ems",
            "params": {
                "basis": "craig_bampton",
                "sweep_port": "SQUID_JJ",
                "squid_flux_sweep": [0.0, 0.25, 0.5]
            }
        }),
    );
    let by_idx = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "hamiltonian",
            "physics": "ems",
            "params": {
                "basis": "craig_bampton",
                "sweep_port": 1,
                "squid_flux_sweep": [0.0, 0.25, 0.5]
            }
        }),
    );
    let pts = by_name["sweep"].as_array().expect("sweep array");
    assert_eq!(pts.len(), 3, "{by_name}");
    assert_eq!(by_name["sweep"], by_idx["sweep"], "sweep_port by name vs index must match");

    let ej0 = pts[0]["e_j_Hz"].as_f64().unwrap();
    let ej_mid = pts[1]["e_j_Hz"].as_f64().unwrap();
    let ej_half = pts[2]["e_j_Hz"].as_f64().unwrap();
    assert!((ej0 - 18.0e9).abs() < 1.0, "got ej0={ej0}");
    assert!(ej_mid < ej0 && ej_mid > ej_half, "expected monotonic SQUID tuning: {ej0} > {ej_mid} > {ej_half}");
    assert!((ej_half - 4.5e9).abs() < 1.0, "d=0.25 at 0.5 Phi0 must yield 0.25*18 GHz = 4.5 GHz, got {ej_half}");

    // Also verify ej_hz_sweep through Engine::call
    let ej_sw = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "hamiltonian",
            "physics": "ems",
            "params": {
                "basis": "craig_bampton",
                "sweep_port": "Q1_JJ",
                "ej_hz_sweep": [12.0e9, 20.0e9]
            }
        }),
    );
    let ej_pts = ej_sw["sweep"].as_array().expect("ej_hz_sweep output");
    assert_eq!(ej_pts.len(), 2);
    assert!((ej_pts[0]["e_j_Hz"].as_f64().unwrap() - 12.0e9).abs() < 1.0);
    assert!((ej_pts[1]["e_j_Hz"].as_f64().unwrap() - 20.0e9).abs() < 1.0);
}

#[test]
fn fr19_fr11_pt1_fr17_open_robin_bc_thickness_correction_and_sub_pf_convergence_e2e() {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "FR-19/11/17 Open Sheet E2E", "dim": 3, "h": 150.0e-6 }));
    let (half, h_sub, h_air) = (350.0e-6, 250.0e-6, 250.0e-6);
    call(&mut e, "add_box", json!({ "name": "sub", "min": [-half, -half, -h_sub], "max": [half, half, 0.0] }));
    call(&mut e, "add_box", json!({ "name": "air", "min": [-half, -half, 0.0], "max": [half, half, h_air] }));
    call(&mut e, "add_rect", json!({ "name": "disk_r", "min": [-90.0e-6, -90.0e-6], "max": [90.0e-6, 90.0e-6] }));
    call(&mut e, "add_sheet", json!({ "name": "pad", "input": "disk_r", "z": 0.0 }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [{ "name": "sub", "input": "sub" }, { "name": "air", "input": "air" }] }),
    );
    call(&mut e, "set_mesh", json!({ "h": 150.0e-6, "adapt_scale": 1.0 }));

    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "sub", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "air", "material": { "eps_r": 1.0 } }));
    // Isolated conductor with NO ground BC, terminated solely by FR-19 spherical 1/r Robin "open" BC!
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "pad", "kind": "terminal", "params": { "name": "pad" } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer", "kind": "open", "params": { "origin": [0.0, 0.0, 0.0] } }),
    );

    call(&mut e, "solve_run", json!({ "isolate": false }));
    let (res, _) = e.result("dfs").unwrap();
    assert!(res.scalar("maxwell_raw_asymmetry_rel").is_some(), "expected maxwell_raw_asymmetry_rel scalar");

    // FR-11 Pt 1: 3D matched-asymptotic finite metal thickness correction (200 nm)
    let cap_raw = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
    let cap_thick = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "capacitance_matrix",
            "physics": "dfs",
            "params": { "metal_thickness_m": 200.0e-9 }
        }),
    );
    let c_raw = cap_raw["maxwell_matrix_F"][0][0].as_f64().unwrap();
    let c_thick = cap_thick["maxwell_matrix_F"][0][0].as_f64().unwrap();
    assert!(c_raw > 1.0e-15 && c_raw < 1.0e-12, "expected sub-pF isolated pad capacitance, got {c_raw:e} F");
    assert!(c_thick > c_raw, "200 nm finite metal thickness must increase capacitance: {c_thick:e} > {c_raw:e}");
    assert!(cap_thick["thickness_correction"]["s_eff_m"].as_f64().unwrap() > 0.0, "{cap_thick}");

    // FR-17: sub-pF convergence_check with quantities filter and adapt_scale (factor^1, not factor^2)
    let conv = call(
        &mut e,
        "convergence_check",
        json!({
            "physics": "dfs",
            "factor": 0.85,
            "scale_sizes": true,
            "tolerance": 0.5,
            "quantities": ["scalar C_pad_pad (F)"]
        }),
    );
    let rows = conv["rows"].as_array().expect("rows");
    assert_eq!(rows.len(), 1, "quantities filter must isolate C[0,0]: {conv}");
    let rel = rows[0]["relative"].as_f64().unwrap();
    assert!(
        rel > 1.0e-6 && rel.is_finite(),
        "sub-pF capacitance must use true relative error, not 1e-12 floor: {conv}"
    );
    let cells_coarse = conv["cells"]["coarse"].as_u64().unwrap() as f64;
    let cells_fine = conv["cells"]["refined"].as_u64().unwrap() as f64;
    let ratio = cells_fine / cells_coarse;
    // At factor=0.85 in 3D, factor^-3 ~ 1.63; factor^-6 (squared adapt_scale bug) would be ~ 2.65+
    assert!(ratio < 2.5, "adapt_scale must not square refinement factor (cells ratio={ratio:.2}): {conv}");
}

#[test]
fn fr05_pt1_fr07_fr10_pt1_fr21_fr23_mesh_guards_admission_timeout_lossless_ibc_and_tls_e2e() {
    // 1. FR-05 Pt 1: sheet3d max_cells guard + decoupled Z-slab controls + project_script round-trip
    let mut e_mesh = Engine::new();
    build_3d_cavity_chip(&mut e_mesh, false, true);
    call(
        &mut e_mesh,
        "set_mesh",
        json!({
            "h": 240.0e-6,
            "dz_min": 40.0e-6,
            "max_z_slabs": 8,
            "z_growth_ratio": 1.4,
            "max_cells": 500_000,
            "max_steiner_multiplier": 20,
            "strict_mesh": false
        }),
    );
    let script = call(&mut e_mesh, "project_script", json!({}));
    assert_eq!(script["form"], json!("commands"));
    let script_txt = script["script"].as_str().unwrap();
    for key in ["dz_min", "max_z_slabs", "z_growth_ratio", "max_cells", "max_steiner_multiplier"] {
        assert!(script_txt.contains(key), "project_script missing {key}: {script_txt}");
    }
    let built = call(&mut e_mesh, "mesh_build", json!({ "wait": true }));
    assert!(built["mesh"]["cells"].as_u64().unwrap() > 0);

    // Pre-allocation max_cells abort on sheet3d
    call(&mut e_mesh, "set_mesh", json!({ "max_cells": 50 }));
    let err_cap = e_mesh.call("mesh_build", &json!({ "wait": true })).unwrap_err().to_string();
    assert!(err_cap.contains("max_cells"), "expected max_cells guard error, got: {err_cap}");

    // 2. FR-21: restore sane max_cells, solve lossless BCS IBC on lossless substrate -> Im(f) == 0
    call(&mut e_mesh, "set_mesh", json!({ "max_cells": 500_000 }));
    call(&mut e_mesh, "solve_run", json!({ "isolate": false }));
    let (ems_res, _) = e_mesh.result("ems").unwrap();
    for sol in &ems_res.solutions {
        let ev = sol.eigenvalue.expect("cavity eigenvalue");
        assert_eq!(ev[1], 0.0, "FR-21 lossless: true BCS IBC must stay on real eigensolve path (Im f == 0)");
    }
    // FR-23: switch substrate to silicon_4k (TOML [interfaces.{ms,sa}]) alongside aluminium_sc_emft ([interfaces.ma]),
    // and check TOML auto-loading + continuous trench_depth_m regularization
    call(&mut e_mesh, "set_material", json!({ "physics": "ems", "region": "si", "material": "silicon_4k" }));
    call(&mut e_mesh, "solve_run", json!({ "isolate": false }));
    let tls_0 = call(&mut e_mesh, "analyse", json!({ "analysis": "tls_participation", "physics": "ems" }));
    let tls_trench = call(
        &mut e_mesh,
        "analyse",
        json!({
            "analysis": "tls_participation",
            "physics": "ems",
            "params": { "trench_depth_m": 80.0e-9 }
        }),
    );
    assert!(tls_0["p_ma"].as_f64().unwrap() > 0.0 && tls_0["q_tls"].as_f64().unwrap().is_finite(), "{tls_0}");
    assert!(
        tls_trench["q_tls"].as_f64().unwrap().is_finite()
            && (tls_trench["p_sa"].as_f64().unwrap() - tls_0["p_sa"].as_f64().unwrap()).abs() > 0.0,
        "trench_depth_m should continuously adjust surface participation: {tls_0} vs {tls_trench}"
    );

    // 3. FR-07: strict_mesh fails when 2D Ruppert refinement hits max_steiner_multiplier cap
    let mut e_strict = Engine::new();
    call(&mut e_strict, "project_new", json!({ "title": "FR-07 Strict Mesh", "dim": 2, "h": 0.05 }));
    call(&mut e_strict, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    call(&mut e_strict, "add_overlay", json!({ "regions": [{ "name": "box", "input": "box" }] }));
    call(&mut e_strict, "set_mesh", json!({ "h": 0.05, "max_steiner_multiplier": 0, "strict_mesh": true }));
    let strict_err = e_strict.call("mesh_build", &json!({ "wait": true })).unwrap_err().to_string();
    assert!(
        strict_err.contains("Ruppert") || strict_err.contains("Steiner") || strict_err.contains("strict"),
        "expected strict_mesh Ruppert cap failure, got: {strict_err}"
    );
    call(&mut e_strict, "set_mesh", json!({ "strict_mesh": false }));
    let relaxed = call(&mut e_strict, "mesh_build", json!({ "wait": true }));
    assert!(relaxed["mesh"]["cells"].as_u64().unwrap() > 0);

    // 4. FR-10 Pt 1: memory_wait_timeout_s through Engine::call("solve_run", ...)
    let mut e_adm = Engine::new();
    call(&mut e_adm, "project_new", json!({ "title": "FR-10 Admission", "dim": 2, "h": 0.25 }));
    call(&mut e_adm, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    call(&mut e_adm, "add_overlay", json!({ "regions": [{ "name": "box", "input": "box" }] }));
    call(&mut e_adm, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary", "order": 1 }));
    call(&mut e_adm, "set_material", json!({ "physics": "dfs", "region": "box", "material": { "k": 1.0 } }));
    call(
        &mut e_adm,
        "add_condition",
        json!({ "physics": "dfs", "on": "box.left", "kind": "temperature", "params": { "value": 300.0 } }),
    );
    call(
        &mut e_adm,
        "add_condition",
        json!({ "physics": "dfs", "on": "box.right", "kind": "temperature", "params": { "value": 310.0 } }),
    );
    let mesh_info = call(&mut e_adm, "mesh_build", json!({ "wait": true }));
    let cells = mesh_info["mesh"]["cells"].as_u64().unwrap();
    let mib: u64 = 1 << 20;
    e_adm.history().add(Row {
        chain: "dfs".into(),
        dim: 2,
        cells: (cells / 2).max(1),
        dofs: cells,
        peak_working_set: 100 * mib,
        peak_commit: 100 * mib,
        seconds: 0.1,
    });
    e_adm.history().add(Row {
        chain: "dfs".into(),
        dim: 2,
        cells: cells.max(2),
        dofs: cells * 2,
        peak_working_set: 400 * mib,
        peak_commit: 400 * mib,
        seconds: 0.2,
    });
    let low_ram_probe: MemoryProbe = Arc::new(move || Some((2000 * mib, 4 * mib)));
    e_adm.set_memory_probe(low_ram_probe);
    let t0 = Instant::now();
    let adm_err =
        e_adm.call("solve_run", &json!({ "isolate": false, "memory_wait_timeout_s": 0.2 })).unwrap_err().to_string();
    let elapsed = t0.elapsed().as_secs_f64();
    assert!(elapsed < 1.2, "admission timeout took too long: {elapsed:.2}s ({adm_err})");
    assert!(
        adm_err.contains("need") && adm_err.contains("free") && adm_err.contains("total") && adm_err.contains("0.2"),
        "timeout error must name need/free/total/timeout: {adm_err}"
    );
}
