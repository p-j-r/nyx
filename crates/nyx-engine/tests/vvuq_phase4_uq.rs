// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 VVUQ integration tests in `nyx-engine`:
//! - End-to-end `Engine::call("analyse", { "analysis": "uq_cqed", ... })` on solved `dfs` CPW (Prager–Synge
//!   `flux_bounds` vs `none`, verifying strict ASME V&V 20 epistemic/aleatoric separation);
//! - Deterministic `SplitMix64` LHS and low-discrepancy Sobol `sweep_run` streaming via HDF5 `SweepCheckpoint`,
//!   eliminating `Vec<Engine>` upfront cloning;
//! - Clean `skipped_invalid_geometry` recording and continuation on invalid sampled geometries (`SimplexMesh::check` /
//!   meshing / parameter check failures) plus `path_guard` confinement on `checkpoint_path`.

use nyx_engine::Engine;
use nyx_storage::SweepCheckpoint;
use serde_json::{json, Value};
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn tmp_dir(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("nyx-vvuq-p4-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

fn solve_cpw_2d_with_flux_bounds(flux_bounds: &str) -> Engine {
    let mut e = Engine::new();
    let call =
        |eng: &mut Engine, cmd: &str, args: Value| eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"));
    call(&mut e, "project_new", json!({ "title": "CPW UQ Test", "dim": 2, "h": 3.0e-6 }));
    call(&mut e, "add_rect", json!({ "name": "substrate", "min": [-40.0e-6, -40.0e-6], "max": [40.0e-6, 0.0] }));
    call(&mut e, "add_rect", json!({ "name": "vacuum", "min": [-40.0e-6, 0.0], "max": [40.0e-6, 40.0e-6] }));
    call(
        &mut e,
        "add_rounded_rect",
        json!({ "name": "trace", "min": [-5.0e-6, -200.0e-9], "max": [5.0e-6, 200.0e-9], "r": 200.0e-9 }),
    );
    call(
        &mut e,
        "add_rounded_rect",
        json!({ "name": "gnd_l", "min": [-38.0e-6, -200.0e-9], "max": [-11.0e-6, 200.0e-9], "r": 200.0e-9 }),
    );
    call(
        &mut e,
        "add_rounded_rect",
        json!({ "name": "gnd_r", "min": [11.0e-6, -200.0e-9], "max": [38.0e-6, 200.0e-9], "r": 200.0e-9 }),
    );
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "substrate", "input": "substrate", "kind": "Material" },
                { "name": "vacuum", "input": "vacuum", "kind": "Material" },
                { "name": "trace", "input": "trace", "kind": "Hole" },
                { "name": "gnd_l", "input": "gnd_l", "kind": "Hole" },
                { "name": "gnd_r", "input": "gnd_r", "kind": "Hole" }
            ]
        }),
    );
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call(
        &mut e,
        "set_study",
        json!({ "physics": "dfs", "kind": "electrostatic", "params": { "flux_bounds": flux_bounds } }),
    );
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": 11.45 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "vacuum", "material": { "eps_r": 1.0 } }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer & in_box([-5.1e-6, -1.0e-6, -3.0e-7], [5.1e-6, 1.0e-6, 3.0e-7])",
            "kind": "terminal",
            "params": { "name": "trace" }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer & (in_box([-38.1e-6, -1.0e-6, -3.0e-7], [-10.9e-6, 1.0e-6, 3.0e-7]) | in_box([10.9e-6, -1.0e-6, -3.0e-7], [38.1e-6, 1.0e-6, 3.0e-7]))",
            "kind": "ground",
            "params": {}
        }),
    );
    call(&mut e, "add_size", json!({ "sel": "kind(circle)", "h": 55.0e-9 }));
    call(&mut e, "solve_run", json!({ "physics": "dfs" }));
    e
}

#[test]
fn engine_analyse_uq_cqed_strict_epistemic_aleatoric_separation() {
    let mut e_none = solve_cpw_2d_with_flux_bounds("none");
    let mut e_bounds = solve_cpw_2d_with_flux_bounds("report");

    let list = e_bounds.call("list_analyses", &json!({})).unwrap();
    let has_uq = list
        .as_array()
        .unwrap()
        .iter()
        .any(|a| a["analysis"] == json!("uq_cqed") && a["applies_to"].as_array().unwrap().contains(&json!("dfs")));
    assert!(has_uq, "list_analyses missing applicable uq_cqed: {list}");

    let q_cfg = json!({
        "qubits": [{ "name": "Q_trace", "pads": ["trace"], "E_J_Hz": 15.0e9 }],
        "seed": 4242,
        "n_samples": 400
    });

    // Epistemic-only run (`uncertainties` omitted): NO fake Gaussian aleatoric sigma may be generated
    let ep_only =
        e_bounds.call("analyse", &json!({ "analysis": "uq_cqed", "physics": "dfs", "params": q_cfg })).unwrap();
    assert_eq!(ep_only["asme_vv20_separated"], json!(true));
    assert_eq!(ep_only["zero_fem_solves"], json!(true));
    let cs_ep = &ep_only["quantities"]["C_Sigma_F"];
    assert!(cs_ep.get("aleatoric").is_none());
    let lo = cs_ep["epistemic"]["low"].as_f64().unwrap();
    let up = cs_ep["epistemic"]["up"].as_f64().unwrap();
    let nom = cs_ep["nominal"].as_f64().unwrap();
    assert!(lo < up && lo <= nom && nom <= up, "expected Prager-Synge bracket [{lo}, {up}] enclosing {nom}");

    // Full joint UQ on both `"none"` and `"report"`: aleatoric sigmas must match identically!
    let mut joint_params = q_cfg.clone();
    joint_params["uncertainties"] = json!({
        "E_J_Hz": { "rel_tol": 0.025 },
        "tan_delta": { "rel_tol": 0.12 },
        "t_ms": { "rel_tol": 0.08 }
    });
    let uq_none =
        e_none.call("analyse", &json!({ "analysis": "uq_cqed", "physics": "dfs", "params": joint_params })).unwrap();
    let uq_bounds =
        e_bounds.call("analyse", &json!({ "analysis": "uq_cqed", "physics": "dfs", "params": joint_params })).unwrap();

    for key in ["f01_Hz", "Q_TLS", "loss_tan_total"] {
        let s_none = uq_none["quantities"][key]["aleatoric"]["linear_sigma"].as_f64().unwrap();
        let s_bounds = uq_bounds["quantities"][key]["aleatoric"]["linear_sigma"].as_f64().unwrap();
        let lin_rel_diff = (s_none - s_bounds).abs() / s_none.abs().max(1e-300);
        assert!(
            lin_rel_diff < 1.0e-8,
            "{key} linear_sigma polluted by Prager-Synge bounds: none={s_none}, report={s_bounds}, rel_diff={lin_rel_diff}"
        );
        let mc_none = uq_none["quantities"][key]["aleatoric"]["mc_sigma"].as_f64().unwrap();
        let mc_bounds = uq_bounds["quantities"][key]["aleatoric"]["mc_sigma"].as_f64().unwrap();
        let mc_rel_diff = (mc_none - mc_bounds).abs() / mc_none.abs().max(1e-300);
        assert!(
            mc_rel_diff < 1.0e-8,
            "{key} mc_sigma polluted by Prager-Synge bounds: none={mc_none}, report={mc_bounds}, rel_diff={mc_rel_diff}"
        );
        let ep_hw = uq_bounds["quantities"][key]["epistemic"]["half_width"].as_f64().unwrap();
        let ep_rel = uq_bounds["quantities"][key]["epistemic"]["rel_half_width"].as_f64().unwrap();
        assert!(ep_hw > 0.0 && ep_rel > 0.01, "expected multi-percent Prager-Synge epistemic bracket on {key}");
        assert!(uq_bounds["quantities"][key]["outer_bracket_95"].is_array());
    }
}

#[test]
fn qmc_lhs_and_sobol_streaming_sweep_determinism_and_hdf5_checkpoint() {
    let dir = tmp_dir("qmc-det");
    let ckpt_lhs1 = dir.join("lhs1.h5");
    let ckpt_lhs2 = dir.join("lhs2.h5");
    let ckpt_sobol = dir.join("sobol2d.h5");

    let run_lhs = |ckpt_path: &PathBuf| {
        let mut e = Engine::new();
        e.call("project_open", &json!({ "path": example("15_soi_wire") })).unwrap();
        let r = e
            .call(
                "sweep_run",
                &json!({
                    "name": "lhs-width",
                    "checkpoint_path": ckpt_path.to_str().unwrap(),
                    "sampling": {
                        "method": "lhs",
                        "n_samples": 6,
                        "seed": 20261126u64,
                        "parameters": [
                            { "name": "width_m", "cmd": "set_param", "args": { "name": "w" }, "field": "value", "min": 4.2e-7, "max": 5.8e-7 }
                        ]
                    }
                }),
            )
            .unwrap();
        assert_eq!(r["sampling"]["method"], json!("lhs"));
        assert_eq!(r["sampling"]["n_samples"], json!(6));
        e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap()
    };

    let t1 = run_lhs(&ckpt_lhs1);
    let t2 = run_lhs(&ckpt_lhs2);
    assert_eq!(t1["complete"], json!(true), "{t1}");
    assert_eq!(t1["solved"], json!(6), "{t1}");
    assert_eq!(t1["skipped_invalid_geometry"], json!(0), "{t1}");

    let pts1 = t1["points"].as_array().unwrap();
    let pts2 = t2["points"].as_array().unwrap();
    assert_eq!(pts1.len(), 6);
    for i in 0..6 {
        let w1 = pts1[i]["params"][0].as_f64().unwrap();
        let w2 = pts2[i]["params"][0].as_f64().unwrap();
        assert_eq!(w1.to_bits(), w2.to_bits(), "LHS param {i} not deterministic");
        let n1 = pts1[i]["numbers"]["ems"]["eigenvalue[0].re"].as_f64().unwrap();
        let n2 = pts2[i]["numbers"]["ems"]["eigenvalue[0].re"].as_f64().unwrap();
        assert_eq!(n1.to_bits(), n2.to_bits(), "LHS eigenvalue {i} not deterministic");
    }

    // Verify HDF5 SweepCheckpoint matches sweep_table 0-ULP on disk
    let h5_lhs = SweepCheckpoint::open(&ckpt_lhs1).unwrap();
    assert_eq!(h5_lhs.param_names(), &["width_m"]);
    assert_eq!(h5_lhs.n_points().unwrap(), 6);
    assert_eq!(h5_lhs.n_skipped().unwrap(), 0);
    let stored_pts = h5_lhs.read_points().unwrap();
    assert_eq!(stored_pts.len(), 6);
    for (i, sp) in stored_pts.iter().enumerate() {
        assert_eq!(sp.params[0].to_bits(), pts1[i]["params"][0].as_f64().unwrap().to_bits());
        assert_eq!(sp.eig_re[0].to_bits(), pts1[i]["numbers"]["ems"]["eigenvalue[0].re"].as_f64().unwrap().to_bits());
    }

    // 2D Sobol QMC sweep (1 uniform geometry param + 1 normal mesh param)
    let mut e_sob = Engine::new();
    e_sob.call("project_open", &json!({ "path": example("15_soi_wire") })).unwrap();
    let r_sob = e_sob
        .call(
            "sweep_run",
            &json!({
                "name": "sobol-2d",
                "checkpoint_path": ckpt_sobol.to_str().unwrap(),
                "checkpoint_lzf": true,
                "sampling": {
                    "method": "sobol",
                    "n_samples": 6,
                    "parameters": [
                        { "name": "width_m", "cmd": "set_param", "args": { "name": "w" }, "field": "value", "min": 4.3e-7, "max": 5.7e-7 },
                        { "name": "ppw", "cmd": "set_mesh", "field": "points_per_wavelength", "distribution": "normal", "mean": 10.0, "sigma": 1.0, "min": 7.0, "max": 14.0 }
                    ]
                }
            }),
        )
        .unwrap();
    let t_sob = e_sob.call("sweep_table", &json!({ "sweep": r_sob["sweep"], "wait_s": 300.0 })).unwrap();
    assert_eq!(t_sob["complete"], json!(true), "{t_sob}");
    assert_eq!(t_sob["solved"], json!(6), "{t_sob}");

    let h5_sob = SweepCheckpoint::open(&ckpt_sobol).unwrap();
    assert_eq!(h5_sob.param_names(), &["width_m", "ppw"]);
    assert_eq!(h5_sob.n_points().unwrap(), 6);

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn qmc_skipped_invalid_geometry_records_and_continues_ensemble_and_guards_path() {
    let dir = tmp_dir("qmc-skip");
    let ckpt_rel = "sub/skip_ckpt.h5";
    std::fs::create_dir_all(dir.join("sub")).unwrap();

    let mut e = Engine::new();
    let call =
        |eng: &mut Engine, cmd: &str, args: Value| eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"));
    call(&mut e, "project_new", json!({ "title": "QMC Invalid Circle", "dim": 2, "h": 0.8e-6 }));
    call(&mut e, "set_param", json!({ "name": "r", "value": 2.5e-6 }));
    call(&mut e, "add_rect", json!({ "name": "dielectric", "min": [-12.0e-6, -12.0e-6], "max": [12.0e-6, 12.0e-6] }));
    call(&mut e, "add_circle", json!({ "name": "pin", "center": [0.0, 0.0], "r": "r" }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "dielectric", "input": "dielectric", "kind": "Material" },
                { "name": "pin", "input": "pin", "kind": "Hole" }
            ]
        }),
    );
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "dielectric", "material": { "eps_r": 11.45 } }));
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer & in_box([-6.0e-6, -6.0e-6, -6.0e-6], [6.0e-6, 6.0e-6, 6.0e-6])",
            "kind": "terminal",
            "params": { "name": "pin" }
        }),
    );
    call(
        &mut e,
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer & !in_box([-6.0e-6, -6.0e-6, -6.0e-6], [6.0e-6, 6.0e-6, 6.0e-6])",
            "kind": "ground",
            "params": {}
        }),
    );
    e.set_path_root(Some(&dir)).unwrap();

    // Path guard must block checkpoint_path escaping set_path_root
    let err_esc = e
        .call(
            "sweep_run",
            &json!({
                "checkpoint_path": "../escape.h5",
                "sampling": {
                    "method": "lhs",
                    "n_samples": 2,
                    "parameters": [{ "name": "r", "cmd": "set_param", "args": { "name": "r" }, "field": "value", "min": 1.0e-6, "max": 4.0e-6 }]
                }
            }),
        )
        .unwrap_err()
        .to_string();
    assert!(err_esc.contains("checkpoint_path"), "{err_esc}");

    // Sample across negative & positive radii: r <= 0 fails geometry checks (skipped_invalid_geometry),
    // while r in [1.0um, 4.75um] evaluates, passes SimplexMesh::check, and solves cleanly.
    let r = e
        .call(
            "sweep_run",
            &json!({
                "name": "qmc-mixed-geom",
                "checkpoint_path": ckpt_rel,
                "sampling": {
                    "method": "sobol",
                    "n_samples": 6,
                    "parameters": [
                        { "name": "radius_m", "cmd": "set_param", "args": { "name": "r" }, "field": "value", "min": -4.0e-6, "max": 6.0e-6 }
                    ]
                }
            }),
        )
        .unwrap();
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap();
    assert_eq!(t["complete"], json!(true), "{t}");
    let n_solved = t["solved"].as_u64().unwrap();
    let n_skipped = t["skipped_invalid_geometry"].as_u64().unwrap();
    assert!(n_skipped > 0, "expected >= 1 skipped invalid geometry: {t}");
    assert!(n_solved > 0, "expected >= 1 solved valid point: {t}");
    assert_eq!(n_solved + n_skipped, 6, "{t}");

    let rows = t["points"].as_array().unwrap();
    assert!(rows.iter().any(|r| r["state"] == json!("skipped_invalid_geometry") && r["reason"].is_string()));
    assert!(rows
        .iter()
        .any(|r| r["state"] == json!("done") && r["numbers"]["dfs"].as_object().is_some_and(|m| !m.is_empty())));

    // Inspect HDF5 checkpoint on disk: solved points + skipped records both persisted
    let h5 = SweepCheckpoint::open(&dir.join(ckpt_rel)).unwrap();
    assert_eq!(h5.n_points().unwrap(), n_solved);
    assert_eq!(h5.n_skipped().unwrap(), n_skipped);
    let sk = h5.read_skipped().unwrap();
    assert_eq!(sk.len(), n_skipped as usize);
    for rec in &sk {
        assert!(rec.params[0] <= 0.0, "expected non-positive invalid radius for skipped pt: {:?}", rec);
        assert!(!rec.reason.is_empty());
    }

    let _ = std::fs::remove_dir_all(&dir);
}
