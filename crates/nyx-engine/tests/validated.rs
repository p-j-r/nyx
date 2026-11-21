// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Validated is a state (productivity plan P6, ADR 91): a result is exploratory until a mesh-convergence check with a tolerance has passed on
//! its problem key, then checked; a stale result is never checked; a result whose check would not fit the machine says so; the verdict is
//! written into the run that solved the problem; every answer that carries `stale` carries `validated`.

use nyx_engine::estimate::Row;
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_str()
        .unwrap()
        .to_string()
}

/// The cantilever, solved: a fast 2D problem whose scalars (max displacement, max von Mises) move with the mesh. (The thermal bar would not do: its scalars are its boundary values and pass any check.)
fn solved() -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("02_cantilever") })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    e
}

fn state(e: &mut Engine) -> Value {
    e.call("result_info", &json!({ "physics": "sms" })).unwrap()[0]["validated"].clone()
}

fn tmp(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("nyx-validated-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&p);
    p
}

#[test]
fn a_result_is_exploratory_until_a_check_with_a_tolerance_passes() {
    let mut e = solved();
    let v = state(&mut e);
    assert_eq!(v["state"], "exploratory", "{v}");
    assert!(v["reason"].as_str().unwrap().contains("no mesh-convergence check"), "{v}");
    assert!(v.get("check").is_none());

    // a check without a tolerance passes nothing and fails nothing: the state does not move and nothing is recorded
    let r = e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6 })).unwrap();
    assert_eq!(r["recorded"], json!(false), "{r}");
    assert_eq!(r["validated"]["state"], "exploratory", "{r}");

    // a generous tolerance passes: checked, with the numbers the verdict rests on
    let r = e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 0.5 })).unwrap();
    assert_eq!(r["converged"], json!(true), "{r}");
    assert_eq!(r["recorded"], json!(true), "{r}");
    assert_eq!(r["validated"]["state"], "checked", "{r}");
    let c = &r["validated"]["check"];
    assert_eq!(
        (c["factor"].as_f64(), c["tolerance"].as_f64(), c["passed"].as_bool()),
        (Some(0.6), Some(0.5), Some(true))
    );
    assert_eq!(c["matching"], "index");
    assert!(c["observable"].as_str().unwrap().contains("matched by index"), "{c}");
    assert!(c["worst_relative"].as_f64().unwrap() <= 0.5);
    assert_eq!(state(&mut e), r["validated"], "result_info says what the check answer said");

    // the last check of the problem decides: an impossible tolerance refutes it
    let r = e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 1e-15 })).unwrap();
    assert_eq!(r["converged"], json!(false), "{r}");
    assert_eq!(r["validated"]["state"], "exploratory", "{r}");
    assert!(r["validated"]["reason"].as_str().unwrap().contains("did not pass"), "{r}");
    assert_eq!(r["validated"]["check"]["passed"], json!(false));
}

#[test]
fn every_answer_that_carries_stale_carries_validated() {
    let mut e = solved();
    e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 0.5 })).unwrap();
    for (name, v) in [
        ("project_info", e.call("project_info", &json!({})).unwrap()["results"][0].clone()),
        ("result_info", e.call("result_info", &json!({})).unwrap()[0].clone()),
        ("result_field", e.call("result_field", &json!({ "physics": "sms", "field": "von Mises" })).unwrap()),
        ("integrate", e.call("integrate", &json!({ "physics": "sms", "expr": "uz" })).unwrap()),
        ("extremum", e.call("extremum", &json!({ "physics": "sms", "expr": "uz", "op": "max" })).unwrap()),
    ] {
        assert!(v.get("stale").is_some(), "{name} has no stale: {v}");
        assert_eq!(v["validated"]["state"], "checked", "{name}: {v}");
    }
}

#[test]
fn a_stale_result_is_never_checked_and_the_verdict_belongs_to_the_problem() {
    let mut e = solved();
    let h0 = e.project().mesh.h;
    e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 0.5 })).unwrap();
    assert_eq!(state(&mut e)["state"], "checked");

    // the mesh changes: the result is stale, and what was checked is another problem
    e.call("set_mesh", &json!({ "h": h0 * 0.8 })).unwrap();
    let v = state(&mut e);
    assert_eq!(v["state"], "exploratory", "{v}");
    assert!(v["reason"].as_str().unwrap().contains("stale"), "{v}");
    // solved again on the new mesh: a new problem, not checked
    e.call("solve_run", &json!({})).unwrap();
    let v = state(&mut e);
    assert_eq!(v["state"], "exploratory", "{v}");
    assert!(v["reason"].as_str().unwrap().contains("no mesh-convergence check"), "{v}");
    // back to the first problem: its check still holds, because the key is the content
    e.call("set_mesh", &json!({ "h": h0 })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    assert_eq!(state(&mut e)["state"], "checked");
}

#[test]
fn a_check_that_cannot_fit_the_machine_is_called_that_and_a_verdict_beats_it() {
    let mut e = solved();
    let cells = e.model().unwrap().mesh().n_cells() as u64;
    // a machine with 4 GB, and two earlier solves of this chain that cost 2 and 4 GB at half and all of these cells (memory about linear in
    // cells): a mesh finer by 0.6 has 2.8 times the cells in 2D, about 11 GB, far above the 90 percent the machine can hold
    e.set_memory_probe(Arc::new(|| Some((4_000_000_000, 4_000_000_000))));
    let row = |cells: u64, bytes: u64| Row {
        chain: "sms".into(),
        dim: 2,
        cells,
        dofs: cells,
        peak_working_set: bytes,
        peak_commit: bytes,
        seconds: 1.0,
    };
    e.history().add(row(cells / 2, 2_000_000_000));
    e.history().add(row(cells, 4_000_000_000));
    let v = state(&mut e);
    assert_eq!(v["state"], "check_not_affordable", "{v}");
    let why = v["reason"].as_str().unwrap();
    assert!(why.contains("GB") && why.contains("earlier run"), "{why}");

    // with no estimate (no earlier runs) nothing is claimed: it is exploratory, not unaffordable
    let mut f = solved();
    f.set_memory_probe(Arc::new(|| Some((4_000_000_000, 4_000_000_000))));
    assert_eq!(f.call("result_info", &json!({})).unwrap()[0]["validated"]["state"], "exploratory");

    // a check that was made (here on a machine that holds it) is evidence; it decides, not the estimate
    e.set_memory_probe(Arc::new(|| Some((u64::MAX / 2, u64::MAX / 2))));
    e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 0.5 })).unwrap();
    e.set_memory_probe(Arc::new(|| Some((4_000_000_000, 4_000_000_000))));
    assert_eq!(state(&mut e)["state"], "checked");
}

#[test]
fn the_verdict_is_written_into_the_run_that_solved_the_problem() {
    let ws = tmp("run");
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("02_cantilever") })).unwrap();
    e.call("store_save", &json!({ "workspace": ws.to_str().unwrap() })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    let before = e.call("store_runs", &json!({})).unwrap();
    assert_eq!(before.as_array().unwrap().len(), 1, "{before}");
    assert!(before[0].get("convergence").is_none(), "a run starts without a verdict: {before}");

    let r = e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 0.5 })).unwrap();
    assert_eq!(r["run"], before[0]["id"], "the answer names the run it was written to: {r}");
    let runs = e.call("store_runs", &json!({})).unwrap();
    assert_eq!(runs.as_array().unwrap().len(), 1, "a check is not a run: {runs}");
    let c = &runs[0]["convergence"];
    assert_eq!(
        (c["factor"].as_f64(), c["tolerance"].as_f64(), c["passed"].as_bool()),
        (Some(0.6), Some(0.5), Some(true)),
        "{runs}"
    );
    assert!(c["observable"].as_str().unwrap().contains("matched by index"));
    // a later, failed check of the same problem replaces it
    e.call("convergence_check", &json!({ "physics": "sms", "factor": 0.6, "tolerance": 1e-15 })).unwrap();
    assert_eq!(e.call("store_runs", &json!({})).unwrap()[0]["convergence"]["passed"], json!(false));
    let _ = std::fs::remove_dir_all(&ws);
}

#[test]
fn a_lossless_tilted_crystal_waveguide_passes_convergence_check_with_coupled_beta_rounding_marked_as_noise() {
    let mut e = Engine::new();
    for (cmd, args) in [
        ("project_new", json!({ "title": "si_on_r_sapphire", "h": 1.2e-7 })),
        ("add_rect", json!({ "name": "sub", "min": [-1.2e-6, -1.2e-6], "max": [1.2e-6, 0.0] })),
        ("add_rect", json!({ "name": "clad", "min": [-1.2e-6, 0.0], "max": [1.2e-6, 0.8e-6] })),
        ("add_rect", json!({ "name": "core", "min": [-2.25e-7, 0.0], "max": [2.25e-7, 2.3e-7] })),
        (
            "add_overlay",
            json!({ "name": "regions", "regions": [
                { "name": "sapphire", "input": "sub" },
                { "name": "air", "input": "clad" },
                { "name": "silicon", "input": "core" }
            ] }),
        ),
        (
            "add_physics",
            json!({
                "physics": "ems",
                "study": "modes",
                "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.4, "n_modes": 2 },
                "order": 2
            }),
        ),
        (
            "set_material",
            json!({
                "physics": "ems",
                "region": "sapphire",
                "material": { "material": "sapphire", "euler_deg": [0.0, 57.6, 0.0] }
            }),
        ),
        ("set_material", json!({ "physics": "ems", "region": "air", "material": "air" })),
        ("set_material", json!({ "physics": "ems", "region": "silicon", "material": "silicon" })),
    ] {
        e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"));
    }
    e.call("solve_run", &json!({})).unwrap();
    let r = e.call("convergence_check", &json!({ "physics": "ems", "factor": 0.6, "tolerance": 1e-3 })).unwrap();
    assert_eq!(r["converged"], json!(true), "{r}");
    assert!(!r["noise_modes"].as_array().unwrap().is_empty(), "lossless tilted-crystal modes are in noise_modes: {r}");
    for row in r["rows"].as_array().unwrap() {
        let q = row["quantity"].as_str().unwrap();
        if q.starts_with("scalar Q ") || q.ends_with(".im") {
            assert_eq!(row["noise"], json!(true), "lossless Q/im row must be marked noise: {row}");
        }
    }
    assert_eq!(r["validated"]["state"], "checked", "{r}");
}
