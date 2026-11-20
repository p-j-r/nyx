// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Warm starts in a sweep (productivity plan P5b, ADR 96): `sweep_run {reuse: "warm"}` starts each point's eigensolves from the modes of the
//! point before. The numbers agree with a cold sweep to the solver's tolerance (not bit for bit, which is why a warm point says so), the first
//! point is cold, the iterations of the later ones are fewer, a point whose predecessor cannot finish is cold and says why, and the three ways
//! a solve can run (in the process, in a worker process, in the queue worker of the workspace) give one answer. A sweep that does not ask is
//! untouched (`sweep.rs`, `worker.rs` and the goldens cover that).

use nyx_engine::shared_queue::SharedQueue;
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn worker() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"))
}

/// The 2D waveguide of example 06, its mesh fixed by `h` so that a change of frequency keeps the mesh (the mesh is not sized by the physics).
fn guide() -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("06_wr90_waveguide") })).unwrap();
    e.call("set_mesh", &json!({ "physics_sizing": false, "h": 0.0025 })).unwrap();
    e
}

fn freq_series(args: Value) -> Value {
    let mut a = json!({
        "series": { "cmd": "set_study", "args": { "physics": "ems", "kind": "modes", "params": {} }, "field": "params", "values": [] },
    });
    // the study's own parameters are replaced as a whole, so the series gives them whole
    let freqs = [19.0e9, 19.4e9, 19.8e9, 20.2e9];
    a["series"]["values"] =
        json!(freqs.iter().map(|f| json!({ "freq": f, "n_modes": 4, "n_target": 0.9 })).collect::<Vec<_>>());
    for (k, v) in args.as_object().unwrap() {
        a[k] = v.clone();
    }
    a
}

fn table(e: &mut Engine, args: Value) -> Value {
    let r = e.call("sweep_run", &freq_series(args)).unwrap_or_else(|err| panic!("{err}"));
    e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap()
}

fn numbers(t: &Value) -> Vec<Value> {
    t["points"].as_array().unwrap().iter().map(|p| p["numbers"].clone()).collect()
}

/// Every number of every point agrees to `tol` relative.
fn agree(a: &[Value], b: &[Value], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        let (x, y) = (x.as_object().unwrap(), y.as_object().unwrap());
        assert_eq!(x.keys().collect::<Vec<_>>(), y.keys().collect::<Vec<_>>());
        for (phys, nx) in x {
            for (name, vx) in nx.as_object().unwrap() {
                let (vx, vy) = (vx.as_f64().unwrap(), y[phys][name].as_f64().unwrap());
                assert!((vx - vy).abs() <= tol * vy.abs().max(1e-300), "point {i} {phys} {name}: {vx} against {vy}");
            }
        }
    }
}

fn iterations(row: &Value) -> usize {
    row["reuse"]["eigensolves"].as_array().unwrap().iter().map(|c| c["iterations"].as_u64().unwrap() as usize).sum()
}

#[test]
fn a_warm_sweep_agrees_with_a_cold_one_starts_later_points_warm_and_says_so() {
    let mut cold_engine = guide();
    let cold = table(&mut cold_engine, json!({}));
    assert!(
        cold["points"].as_array().unwrap().iter().all(|p| p.get("reuse").is_none()),
        "a sweep that does not ask carries no reuse record"
    );
    let mut warm_engine = guide();
    let warm = table(&mut warm_engine, json!({ "reuse": "warm" }));
    assert_eq!(warm["complete"], json!(true), "{warm}");
    agree(&numbers(&warm), &numbers(&cold), 1e-7);
    let rows = warm["points"].as_array().unwrap();
    assert_eq!(rows[0]["reuse"]["warm"], json!(false));
    assert!(rows[0]["reuse"]["why_cold"].as_str().unwrap().contains("first point"), "{}", rows[0]["reuse"]);
    for (i, r) in rows.iter().enumerate().skip(1) {
        assert_eq!(r["reuse"]["warm"], json!(true), "point {i}: {}", r["reuse"]);
        assert_eq!(r["reuse"]["from_point"], json!(i - 1));
        assert_eq!(r["reuse"]["exploratory"], json!(true));
        assert!(
            iterations(r) < iterations(&rows[0]),
            "point {i} took {} cycles, the cold first point {}",
            iterations(r),
            iterations(&rows[0])
        );
    }
}

#[test]
fn a_worker_process_gives_the_numbers_of_the_same_warm_sweep_in_the_process_to_the_bit() {
    let mut here = guide();
    let in_process = table(&mut here, json!({ "reuse": "warm", "isolate": false, "shared": false }));
    let mut there = guide();
    there.set_worker(worker(), false);
    let in_worker = table(&mut there, json!({ "reuse": "warm", "isolate": true, "shared": false }));
    assert_eq!(in_worker["complete"], json!(true), "{in_worker}");
    assert_eq!(numbers(&in_process), numbers(&in_worker), "the vectors went to the worker and back as bytes");
    let rows = in_worker["points"].as_array().unwrap();
    assert!(rows.iter().skip(1).all(|r| r["reuse"]["warm"] == json!(true)), "{in_worker}");
    assert!(rows[0]["reuse"]["warm"] == json!(false));
}

#[test]
fn the_queue_worker_of_the_workspace_gives_the_same_warm_sweep_too() {
    let ws = std::env::temp_dir().join(format!("nyx-warm-{}", std::process::id()));
    std::fs::create_dir_all(&ws).unwrap();
    std::env::set_var("NYX_QUEUE_IDLE_SECS", "3");
    let mut here = guide();
    let in_process = table(&mut here, json!({ "reuse": "warm", "isolate": false, "shared": false }));
    let mut q = guide();
    q.set_worker(worker(), false);
    q.set_shared_queue(&ws, false);
    let shared = table(&mut q, json!({ "reuse": "warm", "shared": true }));
    assert_eq!(shared["complete"], json!(true), "{shared}");
    assert_eq!(numbers(&in_process), numbers(&shared));
    assert!(shared["points"].as_array().unwrap().iter().skip(1).all(|r| r["reuse"]["warm"] == json!(true)), "{shared}");
    assert!(SharedQueue::new(&ws).items().is_empty(), "the points were collected");
}

#[test]
fn a_point_whose_predecessor_cannot_finish_is_cold_and_the_lane_is_not_held_up() {
    let mut e = guide();
    // all four held: only the last is released, so its predecessor will not finish while it runs
    let r = e
        .call("sweep_run", &freq_series(json!({ "reuse": "warm", "hold": true, "isolate": false, "shared": false })))
        .unwrap();
    let jobs: Vec<u64> = r["points"].as_array().unwrap().iter().map(|p| p["job"].as_u64().unwrap()).collect();
    assert!(e.call("queue_release", &json!({ "id": jobs[3] })).unwrap()["released"] == json!(true));
    e.call("queue_wait", &json!({ "id": jobs[3], "timeout_s": 60.0 })).unwrap();
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"] })).unwrap();
    let last = &t["points"][3];
    assert_eq!(last["state"], json!("done"), "{last}");
    assert_eq!(last["reuse"]["warm"], json!(false), "{last}");
    assert!(last["reuse"]["why_cold"].as_str().unwrap().contains("had not finished"), "{last}");
    // and the sweep goes on when the others are released: point 0 first, then the rest start from their predecessors
    for j in &jobs[..3] {
        e.call("queue_release", &json!({ "id": j })).unwrap();
    }
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 120.0 })).unwrap();
    assert_eq!(t["complete"], json!(true), "{t}");
    assert_eq!(t["points"][1]["reuse"]["warm"], json!(true), "{}", t["points"][1]);
}
