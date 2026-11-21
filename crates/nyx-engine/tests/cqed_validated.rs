// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end superconducting CPW workflow through `Engine::call`: builds the 2D cross-section,
//! solves `dfs` `electrostatic` with `terminal` + `ground`, confirms `convergence_check` moves
//! `validated` from `exploratory` to `checked`, runs `capacitance_matrix` and `tls_participation`,
//! and exercises `bcs_surface_impedance`.

use nyx_engine::Engine;
use serde_json::json;

#[test]
fn cpw_cross_section_workflow_reaches_checked_and_runs_cqed_analyses() {
    let mut e = Engine::new();
    let call = |eng: &mut Engine, cmd: &str, args: serde_json::Value| {
        eng.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"))
    };

    call(&mut e, "project_new", json!({ "title": "Superconducting CPW 2D", "dim": 2, "h": 2.5e-6 }));
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

    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
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
    call(
        &mut e,
        "add_size",
        json!({
            "sel": "kind(circle)",
            "h": 40.0e-9
        }),
    );

    call(&mut e, "check_project", json!({}));
    call(&mut e, "solve_run", json!({ "physics": "dfs" }));

    let info0 = call(&mut e, "result_info", json!({ "physics": "dfs" }));
    assert_eq!(info0[0]["validated"]["state"], "exploratory", "{info0}");

    let conv = call(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.6, "tolerance": 0.02 }));
    assert_eq!(conv["converged"], json!(true), "{conv}");
    assert_eq!(conv["validated"]["state"], "checked", "{conv}");

    let cap = call(&mut e, "analyse", json!({ "analysis": "capacitance_matrix", "physics": "dfs" }));
    assert_eq!(cap["validated"]["state"], "checked", "{cap}");
    let c11 = cap["maxwell_matrix_F"][0][0].as_f64().unwrap();
    assert!(c11 > 1.0e-12 && c11 < 1.0e-9, "unexpected CPW C' = {c11}");

    let tls = call(
        &mut e,
        "analyse",
        json!({
            "analysis": "tls_participation",
            "physics": "dfs",
            "params": {
                "solution": 0,
                "t_ma": 2.0e-9,
                "t_ms": 2.0e-9,
                "t_sa": 2.0e-9,
                "eps_ma": 10.0,
                "eps_ms": 11.45,
                "eps_sa": 10.0,
                "tan_delta_ma": 1.0e-3,
                "tan_delta_ms": 1.0e-3,
                "tan_delta_sa": 1.0e-3
            }
        }),
    );
    assert_eq!(tls["validated"]["state"], "checked", "{tls}");
    assert!(tls["p_ma"].as_f64().unwrap() > 0.0, "{tls}");
    assert!(tls["p_ms"].as_f64().unwrap() > 0.0, "{tls}");
    assert!(tls["p_sa"].as_f64().unwrap() > 0.0, "{tls}");
    assert!(tls["q_tls"].as_f64().unwrap() > 1.0e3, "{tls}");

    let bcs = call(
        &mut e,
        "bcs_surface_impedance",
        json!({
            "material": "niobium",
            "freq": 7.0e9,
            "thickness": 100.0e-9,
            "temperature": 0.01,
            "film_model": "kerr"
        }),
    );
    assert!(bcs["L_k_sq"].as_f64().unwrap() > 0.0, "{bcs}");
}
