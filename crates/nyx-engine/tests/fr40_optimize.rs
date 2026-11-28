// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end verification of FR-40 (`sensitivity_run` and `optimize_run` in `nyx-engine`, ADR 142):
//! * `lbfgs` shape optimization tuning a photonic waveguide width `w` to match a target `n_eff`
//!   (with `checkpoint_path` confinement and `nyx-store` run recording).
//! * 2D and 3D `mma` + unstructured Helmholtz PDE filter + Heaviside `beta`-continuation topology
//!   optimization (2D structural/thermal splitter AND 3D 3-port electromagnetic Y-splitter).
//! * Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` bit-exact determinism.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

fn temp_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-fr40-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn setup_soi_waveguide(e: &mut Engine, w_m: f64) {
    e.call("project_new", &json!({ "title": "SOI Waveguide Shape Opt", "dim": 2, "h": 1.6e-7 })).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": w_m })).unwrap();
    e.call("add_rect", &json!({ "name": "clad", "min": [-1.5e-6, -1.0e-6], "max": [1.5e-6, 1.0e-6] })).unwrap();
    e.call("add_rect", &json!({ "name": "core", "min": ["-w/2", -1.1e-7], "max": ["w/2", 1.1e-7] })).unwrap();
    e.call(
        "add_overlay",
        &json!({
            "regions": [
                { "name": "oxide", "input": "clad" },
                { "name": "silicon", "input": "core" }
            ]
        }),
    )
    .unwrap();
    e.call(
        "add_physics",
        &json!({
            "physics": "ems",
            "study": "modes",
            "study_params": { "freq": 1.934e14, "n_target": 3.4, "n_modes": 2 },
            "order": 1
        }),
    )
    .unwrap();
    e.call("set_material", &json!({ "physics": "ems", "region": "oxide", "material": { "eps_r": 2.085 } })).unwrap();
    e.call("set_material", &json!({ "physics": "ems", "region": "silicon", "material": { "eps_r": 12.08 } })).unwrap();
}

#[test]
fn lbfgs_waveguide_shape_optimization_converges_and_records_checkpoint() {
    let root = temp_dir("lbfgs-wg");

    // First compute the exact reference neff at w_target = 520 nm on a reference engine
    let mut ref_eng = Engine::new();
    setup_soi_waveguide(&mut ref_eng, 520e-9);
    let ref_sens = ref_eng
        .call(
            "sensitivity_run",
            &json!({
                "physics": "ems",
                "objective": { "kind": "effective_index", "mode": 0 },
                "params": [{ "kind": "param_var", "name": "w" }]
            }),
        )
        .unwrap();
    let target_neff = ref_sens["objective"].as_f64().unwrap();
    assert!(target_neff > 2.0 && target_neff < 3.5);

    // Start from detuned width w_init = 400 nm and optimize to match target_neff via L-BFGS
    let mut e = Engine::new();
    setup_soi_waveguide(&mut e, 400e-9);
    e.call("store_save", &json!({ "workspace": root.join("ws").to_str().unwrap() })).unwrap();
    e.set_path_root(Some(&root)).unwrap();

    // Path guard must reject escaping checkpoint_path
    assert!(e
        .call(
            "optimize_run",
            &json!({
                "objective": { "kind": "effective_index", "mode": 0, "target_neff": target_neff },
                "params": [{ "kind": "param_var", "name": "w", "min": 300e-9, "max": 650e-9 }],
                "checkpoint_path": "../escaped_opt.json"
            })
        )
        .is_err());

    let out = e
        .call(
            "optimize_run",
            &json!({
                "physics": "ems",
                "algorithm": "lbfgs",
                "max_iters": 20,
                "gtol": 1e-9,
                "objective": { "kind": "effective_index", "mode": 0, "target_neff": target_neff },
                "params": [{ "kind": "param_var", "name": "w", "min": 300e-9, "max": 650e-9 }],
                "checkpoint_path": "opt_checkpoint.json"
            }),
        )
        .unwrap();

    let init_obj = out["initial_objective"].as_f64().unwrap();
    let final_obj = out["final_objective"].as_f64().unwrap();
    let w_final = out["final_params"]["w"].as_f64().unwrap();

    assert!(init_obj > 1e-3, "initial mismatch must be non-trivial, got {init_obj:e}");
    assert!(
        final_obj < 1e-8 * init_obj,
        "L-BFGS must reduce squared neff mismatch by > 8 orders of magnitude: {init_obj:e} -> {final_obj:e}"
    );
    assert!(
        (w_final - 520e-9).abs() < 35e-9,
        "optimized waveguide width {w_final:e} must converge cleanly into the target regime"
    );
    assert!(root.join("opt_checkpoint.json").is_file(), "checkpoint JSON must be written inside path_root");
    assert!(out.get("run").is_some_and(|r| r.get("id").is_some()), "nyx-store run must be recorded");

    // Verify the project's `w` parameter was updated and stored result is present
    let info = e.call("result_info", &json!({ "physics": "ems" })).unwrap();
    assert_eq!(info.as_array().unwrap().len(), 1);

    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn mma_2d_and_3d_splitter_topology_optimization_with_helmholtz_and_beta_continuation() {
    // 1. 2D Y-splitter / dual-sink thermal conduction topology optimization with MMA + Helmholtz PDE filter + beta continuation
    let mut e2 = Engine::new();
    e2.call("project_new", &json!({ "title": "2D Splitter Topology", "dim": 2, "h": 0.14 })).unwrap();
    e2.call("add_rect", &json!({ "name": "design", "min": [0.0, -0.5], "max": [1.0, 0.5] })).unwrap();
    e2.call("add_overlay", &json!({ "regions": [{ "name": "design", "input": "design" }] })).unwrap();
    e2.call("add_physics", &json!({ "physics": "dfs", "study": "heat_stationary", "order": 1 })).unwrap();
    e2.call("set_material", &json!({ "physics": "dfs", "region": "design", "material": { "k": 1.0, "Q": 1.0e4 } }))
        .unwrap();
    // Two symmetric cold output ports on x = 1.0 (top and bottom branches of the splitter)
    e2.call(
        "add_condition",
        &json!({
            "physics": "dfs",
            "on": "plane(x = 1.0)",
            "kind": "temperature",
            "params": { "value": 300.0 }
        }),
    )
    .unwrap();

    let out2 = e2
        .call(
            "optimize_run",
            &json!({
                "physics": "dfs",
                "algorithm": "mma",
                "sense": "minimize",
                "move_limit": 0.25,
                "xtol": 1e-5,
                "objective": { "kind": "compliance" },
                "params": [{
                    "kind": "topology_density",
                    "domain": "design",
                    "filter_radius": 0.16,
                    "simp_penalty": 3.0,
                    "beta": 1.0,
                    "eta": 0.5,
                    "val_min": 0.05,
                    "val_max": 20.0,
                    "property": "k"
                }],
                "constraints": [{
                    "kind": "volume_fraction",
                    "param_index": 0,
                    "max_fraction": 0.45
                }],
                "continuation_schedule": [
                    { "beta": 1.0, "max_iters": 6 },
                    { "beta": 4.0, "max_iters": 6 },
                    { "beta": 8.0, "max_iters": 6 }
                ]
            }),
        )
        .unwrap();

    let j0_2d = out2["initial_objective"].as_f64().unwrap();
    let jf_2d = out2["final_objective"].as_f64().unwrap();
    let vf_2d = out2["volume_fractions"]["0"].as_f64().unwrap();
    assert!(jf_2d < j0_2d, "2D MMA + Helmholtz + beta-continuation must reduce compliance: {j0_2d:e} -> {jf_2d:e}");
    assert!(vf_2d <= 0.45 + 1e-2, "2D MMA must satisfy volume fraction constraint (got {vf_2d})");
    assert_eq!(out2["stages"].as_array().unwrap().len(), 3);

    // Stored SolveResult must carry the `"topology_density"` cell field
    let info2 = e2.call("result_info", &json!({ "physics": "dfs" })).unwrap();
    let cfields2 = info2[0]["solutions"][0]["cell_fields"].as_array().unwrap();
    assert!(
        cfields2.iter().any(|cf| cf["name"] == "topology_density"),
        "stored SolveResult must expose `topology_density` CellField"
    );

    // 2. 3D Electromagnetic Y-Splitter (`ems::driven`, 1 input waveport -> 2 output waveports, `SplitterBalance` / `Transmission`)
    let mut e3 = Engine::new();
    e3.call("project_new", &json!({ "title": "3D EM Waveguide Splitter", "dim": 3, "h": 0.0085 })).unwrap();
    e3.call("add_box", &json!({ "name": "in_wg", "min": [-0.010, 0.0, 0.0], "max": [0.010, 0.008, 0.010] })).unwrap();
    e3.call("add_box", &json!({ "name": "core_box", "min": [-0.022, 0.008, 0.0], "max": [0.022, 0.020, 0.010] }))
        .unwrap();
    e3.call("add_box", &json!({ "name": "out_l", "min": [-0.022, 0.020, 0.0], "max": [-0.002, 0.028, 0.010] }))
        .unwrap();
    e3.call("add_box", &json!({ "name": "out_r", "min": [0.002, 0.020, 0.0], "max": [0.022, 0.028, 0.010] })).unwrap();
    e3.call(
        "add_overlay",
        &json!({
            "regions": [
                { "name": "in_wg", "input": "in_wg" },
                { "name": "design", "input": "core_box" },
                { "name": "out_l", "input": "out_l" },
                { "name": "out_r", "input": "out_r" }
            ]
        }),
    )
    .unwrap();
    e3.call(
        "set_mesh",
        &json!({
            "h": 0.007,
            "sweep": [
                { "thickness": 0.008, "n": 2 },
                { "thickness": 0.012, "n": 3 },
                { "thickness": 0.008, "n": 2 }
            ]
        }),
    )
    .unwrap();
    e3.call(
        "add_physics",
        &json!({
            "physics": "ems",
            "study": "driven",
            "study_params": {
                "f_min": 10.5e9,
                "f_max": 10.5e9,
                "n_points": 1,
                "prom_tol": 0.0
            },
            "order": 1
        }),
    )
    .unwrap();
    for reg in ["in_wg", "design", "out_l", "out_r"] {
        e3.call("set_material", &json!({ "physics": "ems", "region": reg, "material": { "eps_r": 1.0 } })).unwrap();
    }
    e3.call(
        "add_condition",
        &json!({
            "physics": "ems",
            "on": "boundary(region(in_wg)) & plane(y = 0.0)",
            "kind": "waveport",
            "params": { "name": "p_in", "mode": 1 }
        }),
    )
    .unwrap();
    e3.call(
        "add_condition",
        &json!({
            "physics": "ems",
            "on": "boundary(region(out_l)) & plane(y = 0.028)",
            "kind": "waveport",
            "params": { "name": "p_out1", "mode": 1 }
        }),
    )
    .unwrap();
    e3.call(
        "add_condition",
        &json!({
            "physics": "ems",
            "on": "boundary(region(out_r)) & plane(y = 0.028)",
            "kind": "waveport",
            "params": { "name": "p_out2", "mode": 1 }
        }),
    )
    .unwrap();

    let out3 = e3
        .call(
            "optimize_run",
            &json!({
                "physics": "ems",
                "algorithm": "mma",
                "sense": "minimize",
                "move_limit": 0.25,
                "xtol": 1e-5,
                "objective": {
                    "kind": "splitter_balance",
                    "in_port": 0,
                    "out_ports": [1, 2],
                    "targets": [0.45, 0.45],
                    "freq_index": 0
                },
                "params": [{
                    "kind": "topology_density",
                    "domain": "design",
                    "filter_radius": 0.009,
                    "simp_penalty": 1.0,
                    "beta": 1.0,
                    "eta": 0.5,
                    "val_min": 1.0,
                    "val_max": 4.0
                }],
                "constraints": [{
                    "kind": "volume_fraction",
                    "param_index": 0,
                    "max_fraction": 0.55
                }],
                "continuation_schedule": [
                    { "beta": 1.0, "max_iters": 4 },
                    { "beta": 3.0, "max_iters": 4 }
                ]
            }),
        )
        .unwrap();

    let j0_3d = out3["initial_objective"].as_f64().unwrap();
    let jf_3d = out3["final_objective"].as_f64().unwrap();
    let vf_3d = out3["volume_fractions"]["0"].as_f64().unwrap();
    assert!(
        jf_3d < j0_3d,
        "3D EM Y-splitter MMA + Helmholtz + beta-continuation must reduce splitter mismatch: {j0_3d:e} -> {jf_3d:e}"
    );
    assert!(vf_3d <= 0.55 + 1e-2, "3D EM Y-splitter volume fraction must satisfy constraint (got {vf_3d})");
}

#[test]
fn sensitivity_and_optimize_are_bit_deterministic_across_thread_counts() {
    if std::env::var("NYX_FR40_DET_SUBPROC").is_ok() {
        let mut e = Engine::new();
        setup_soi_waveguide(&mut e, 440e-9);
        let out = e
            .call(
                "optimize_run",
                &json!({
                    "physics": "ems",
                    "algorithm": "mma",
                    "sense": "maximize",
                    "max_iters": 3,
                    "objective": { "kind": "effective_index", "mode": 0 },
                    "params": [{
                        "kind": "topology_density",
                        "domain": "silicon",
                        "filter_radius": 1.5e-7,
                        "simp_penalty": 1.0,
                        "beta": 2.0,
                        "eta": 0.5,
                        "val_min": 2.1,
                        "val_max": 12.1,
                        "initial_density": 0.6
                    }],
                    "constraints": [{
                        "kind": "volume_fraction",
                        "param_index": 0,
                        "max_fraction": 0.7
                    }]
                }),
            )
            .unwrap();
        let jf: f64 = out["final_objective"].as_f64().unwrap();
        let hist: Vec<String> = out["history"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| format!("{:016x}", v.as_f64().unwrap().to_bits()))
            .collect();
        println!("FR40_DET={:016x}|{}", jf.to_bits(), hist.join(","));
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_nt = |nt: &str| -> Value {
        let out = std::process::Command::new(&exe)
            .args(["--exact", "sensitivity_and_optimize_are_bit_deterministic_across_thread_counts", "--nocapture"])
            .env("NYX_FR40_DET_SUBPROC", "1")
            .env("NYX_THREADS", nt)
            .output()
            .unwrap();
        assert!(out.status.success(), "subprocess NYX_THREADS={nt} failed: {}", String::from_utf8_lossy(&out.stderr));
        let line =
            String::from_utf8(out.stdout).unwrap().lines().find(|l| l.starts_with("FR40_DET=")).unwrap().to_string();
        json!(line)
    };

    let d1 = run_nt("1");
    let d8 = run_nt("8");
    assert_eq!(d1, d8, "optimize_run trajectory must be bit-identical under NYX_THREADS=1 vs 8");
}
