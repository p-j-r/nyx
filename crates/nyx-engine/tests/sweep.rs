// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Sweeps (productivity plan Q6, ADR 87): one run per point, recorded; a point is a list of edits, so the case decides what varies. A sweep's
//! numbers are those of a hand-made solve with the same edits, in any order the queue runs them; it never touches the engine's own project,
//! mesh or results; a bad point stops the request before anything runs; its points wait behind a person's request.

use nyx_engine::{Engine, MeshSettings, MeshedModel, Mesher};
use nyx_physics::Error;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// A mesher that counts its calls.
struct Counting {
    inner: nyx_engine::FemcoreMesher,
    calls: AtomicUsize,
}

impl Mesher for Counting {
    fn mesh(
        &self,
        geometry: &geom::model::FeatureList,
        settings: &MeshSettings,
        sizing: Option<&nyx_physics::Sizing>,
    ) -> Result<MeshedModel, Error> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        self.inner.mesh(geometry, settings, sizing)
    }
}

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn open(name: &str) -> (Engine, Arc<Counting>) {
    let mut e = Engine::new();
    let m = Arc::new(Counting { inner: nyx_engine::FemcoreMesher, calls: AtomicUsize::new(0) });
    e.set_mesher(m.clone());
    e.call("project_open", &json!({ "path": example(name) })).unwrap();
    (e, m)
}

fn tmp(name: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!("nyx-sweep-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&p);
    p
}

/// The sweep of the wire's width as a series, and its table once every point has finished.
fn width_sweep(e: &mut Engine, widths: &[f64]) -> Value {
    let r = e
        .call(
            "sweep_run",
            &json!({ "series": { "cmd": "set_param", "args": { "name": "w" }, "field": "value", "values": widths } }),
        )
        .unwrap();
    e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap()
}

/// The first eigenvalue (real part) of the ems result a hand-made solve of the wire gives at width `w`.
fn by_hand(w: f64) -> f64 {
    let (mut e, _) = open("15_soi_wire");
    e.call("set_param", &json!({ "name": "w", "value": w })).unwrap();
    e.call("solve_run", &json!({})).unwrap();
    e.result("ems").unwrap().0.solutions[0].eigenvalue.unwrap()[0]
}

#[test]
fn a_sweep_gives_one_row_per_point_with_the_numbers_of_a_solve_made_by_hand() {
    let widths = [4.5e-7, 5.0e-7, 5.5e-7];
    let (mut e, _) = open("15_soi_wire");
    let t = width_sweep(&mut e, &widths);
    assert_eq!(t["complete"], json!(true), "{t}");
    let rows = t["points"].as_array().unwrap();
    assert_eq!(rows.len(), 3);
    for (row, w) in rows.iter().zip(widths) {
        assert_eq!(row["state"], json!("done"), "{row}");
        assert_eq!(row["edits"][0]["args"]["value"], json!(w));
        let n = row["numbers"]["ems"]["eigenvalue[0].re"].as_f64().unwrap();
        assert_eq!(n, by_hand(w), "width {w}: the sweep and a hand-made solve differ");
        // the metrics are what a recorded run of the solve carries (the runs table keys), for a front-end that records the point itself
        assert_eq!(row["metrics"]["eigenvalue_re"].as_f64(), Some(n), "{row}");
        assert!(row["metrics"]["dofs"].as_f64().is_some_and(|d| d > 0.0), "{row}");
    }
    let (a, c) = (&rows[0]["numbers"]["ems"]["eigenvalue[0].re"], &rows[2]["numbers"]["ems"]["eigenvalue[0].re"]);
    assert_ne!(a, c, "the widths differ, so the test can see a mix-up");
}

#[test]
fn a_point_may_vary_anything_an_edit_can_change() {
    // two edits in one point (a geometry parameter and the element size), and a series over a mesh field
    let (mut e, _) = open("15_soi_wire");
    let r = e
        .call(
            "sweep_run",
            &json!({ "points": [
                { "label": "narrow, coarse", "edits": [
                    { "cmd": "set_param", "args": { "name": "w", "value": 4.5e-7 } },
                    { "cmd": "set_mesh", "args": { "points_per_wavelength": 8.0 } } ] },
                { "label": "as it is", "edits": [] } ] }),
        )
        .unwrap();
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap();
    assert_eq!(t["complete"], json!(true), "{t}");
    assert_eq!(t["points"][0]["label"], json!("narrow, coarse"));
    assert!(t["points"][0]["state"] == json!("done") && t["points"][1]["state"] == json!("done"), "{t}");
    // the point with no edit is the project as it is
    assert_eq!(t["points"][1]["numbers"]["ems"]["eigenvalue[0].re"].as_f64().unwrap(), by_hand(5.0e-7));
    assert_ne!(t["points"][0]["cells"], t["points"][1]["cells"], "the first point has its own mesh");
}

#[test]
fn a_sweep_never_touches_the_project_the_mesh_or_the_results() {
    let (mut e, _) = open("15_soi_wire");
    e.call("solve_run", &json!({})).unwrap();
    let (stamp, mesh, before) =
        (e.stamp(), e.mesh_key(), e.result("ems").unwrap().0.solutions[0].eigenvalue.unwrap()[0]);
    let t = width_sweep(&mut e, &[4.5e-7, 5.5e-7]);
    assert_eq!(t["complete"], json!(true), "{t}");
    assert_eq!(e.stamp(), stamp, "the project changed");
    assert_eq!(e.mesh_key(), mesh);
    let (r, stale) = e.result("ems").unwrap();
    assert!(!stale, "the result became stale");
    assert_eq!(r.solutions[0].eigenvalue.unwrap()[0], before, "a point replaced the engine's result");
}

#[test]
fn a_bad_point_fails_the_whole_request_and_nothing_runs() {
    let (mut e, _) = open("15_soi_wire");
    let ask = |edits: Value| {
        json!({ "points": [
            { "edits": [{ "cmd": "set_param", "args": { "name": "w", "value": 5e-7 } }] },
            { "edits": edits } ] })
    };
    // an edit that fails, one that is not an edit, one that replaces the project, and a missing cmd: each names its point
    for (edits, word) in [
        (json!([{ "cmd": "set_param", "args": {} }]), "point 1"),
        (json!([{ "cmd": "solve_run", "args": {} }]), "not an edit"),
        (json!([{ "cmd": "project_new", "args": {} }]), "not an edit"),
        (json!([{ "cmd": "open_results", "args": { "path": "x.h5" } }]), "not an edit"),
        (json!([{ "cmd": "store_open", "args": { "folder": "x" } }]), "not an edit"),
        (json!([{ "args": {} }]), "cmd"),
    ] {
        let err = e.call("sweep_run", &ask(edits)).unwrap_err().to_string();
        assert!(err.contains("point 1") && err.contains(word), "{err}");
    }
    assert!(e.jobs().list().is_empty(), "something was submitted: {:?}", e.jobs().list());
    assert_eq!(e.call("sweep_table", &json!({})).unwrap(), json!([]));
    // no points, too many, both forms, neither
    for args in [
        json!({ "points": [] }),
        json!({ "series": { "cmd": "set_param", "args": { "name": "w" }, "field": "value", "values": vec![5e-7; 501] } }),
        json!({ "points": [{ "edits": [] }], "series": { "cmd": "set_param", "field": "value", "values": [1] } }),
        json!({}),
    ] {
        assert!(e.call("sweep_run", &args).is_err(), "{args}");
    }
    assert!(e.jobs().list().is_empty());
}

#[test]
fn points_that_share_a_mesh_build_it_once() {
    let (mut e, mesher) = open("15_soi_wire");
    let r = e
        .call("sweep_run", &json!({ "points": [{ "edits": [] }, { "edits": [] }, { "edits": [{ "cmd": "set_param", "args": { "name": "w", "value": 4.5e-7 } }] }] }))
        .unwrap();
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap();
    assert_eq!(t["complete"], json!(true), "{t}");
    assert_eq!(mesher.calls.load(Ordering::SeqCst), 2, "two meshes: the shared one and the narrow wire's");
    assert_eq!(t["points"][0]["numbers"], t["points"][1]["numbers"]);
}

#[test]
fn the_points_of_a_sweep_wait_behind_a_request_made_after_them() {
    let (mut e, _) = open("15_soi_wire");
    e.jobs().pause();
    let r = e
        .call("sweep_run", &json!({ "series": { "cmd": "set_param", "args": { "name": "w" }, "field": "value", "values": [4.5e-7, 5.5e-7] } }))
        .unwrap();
    let solve = e.start_solve(None).unwrap(); // an agent's or a script's single request: priority 2, after the sweep was queued
    let position = |id: u64| e.jobs().info(id).unwrap().position.unwrap();
    for p in r["points"].as_array().unwrap() {
        assert!(position(solve) < position(p["job"].as_u64().unwrap()), "the sweep's point went before the request");
    }
    e.jobs().resume();
    e.wait(solve).unwrap();
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap();
    assert_eq!(t["complete"], json!(true), "{t}");
}

#[test]
fn cancelled_points_never_run_and_leave_no_run_record() {
    let ws = tmp("cancel");
    let (mut e, _) = open("15_soi_wire");
    e.call("store_save", &json!({ "workspace": ws.to_str().unwrap() })).unwrap();
    e.jobs().pause();
    let r = e
        .call("sweep_run", &json!({ "series": { "cmd": "set_param", "args": { "name": "w" }, "field": "value", "values": [4.5e-7, 5.0e-7, 5.5e-7] } }))
        .unwrap();
    assert_eq!(e.call("sweep_cancel", &json!({ "sweep": r["sweep"] })).unwrap()["cancelled"], json!(3));
    e.jobs().resume();
    let t = e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 60.0 })).unwrap();
    assert!(t["points"].as_array().unwrap().iter().all(|p| p["state"] == json!("cancelled")), "{t}");
    assert_eq!(e.call("store_runs", &json!({})).unwrap(), json!([]), "a point that never ran left a record");
}

#[test]
fn each_point_is_recorded_as_a_run_and_the_experiment_keeps_its_own_project() {
    let ws = tmp("runs");
    let (mut e, _) = open("15_soi_wire");
    e.call("store_save", &json!({ "workspace": ws.to_str().unwrap() })).unwrap();
    let before = e.call("store_save", &json!({})).unwrap();
    assert_eq!(before["changed"], json!(false), "{before}");
    let t = width_sweep(&mut e, &[4.5e-7, 5.5e-7]);
    assert_eq!(t["complete"], json!(true), "{t}");
    let runs = e.call("store_runs", &json!({})).unwrap();
    let runs = runs.as_array().unwrap();
    assert_eq!(runs.len(), 2, "{runs:?}");
    // each run carries the parameter it was solved at, and the row says which run it is
    let widths: Vec<f64> = runs.iter().map(|r| r["params"]["w"].as_f64().unwrap()).collect();
    assert_eq!(widths, vec![4.5e-7, 5.5e-7]);
    assert_eq!(t["points"][0]["run"]["id"], runs[0]["id"]);
    // the experiment's project is still the engine's: recording a point did not commit the point's project
    let after = e.call("store_save", &json!({})).unwrap();
    assert_eq!(after["changed"], json!(false), "{after}");
    assert_eq!(after["version"], before["version"]);
}

#[test]
fn the_numbers_of_a_sweep_do_not_depend_on_the_order_of_its_points() {
    let table = |widths: &[f64]| {
        let (mut e, _) = open("15_soi_wire");
        let t = width_sweep(&mut e, widths);
        let mut out: Vec<(String, String)> = t["points"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| (p["edits"][0]["args"]["value"].to_string(), p["numbers"].to_string()))
            .collect();
        out.sort();
        out
    };
    assert_eq!(table(&[4.5e-7, 5.0e-7, 5.5e-7]), table(&[5.5e-7, 4.5e-7, 5.0e-7]));
}

#[test]
fn the_sweep_commands_are_in_the_catalog_and_a_table_without_a_sweep_lists_them() {
    let (mut e, _) = open("15_soi_wire");
    let names: Vec<String> = e.catalog().into_iter().map(|c| c.name).collect();
    for n in ["sweep_run", "sweep_table", "sweep_cancel"] {
        assert!(names.iter().any(|x| x == n), "{n}");
    }
    let r = e.call("sweep_run", &json!({ "name": "one", "points": [{ "edits": [] }] })).unwrap();
    let list = e.call("sweep_table", &json!({})).unwrap();
    assert_eq!(list[0]["name"], json!("one"));
    assert_eq!(list[0]["sweep"], r["sweep"]);
    e.call("sweep_table", &json!({ "sweep": r["sweep"], "wait_s": 300.0 })).unwrap();
    assert!(e.call("sweep_table", &json!({ "sweep": 99 })).is_err());
}
