// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What a person sees while editing (productivity plan P3, ADR 92): the project checked with edits that are not applied yet, the cells and the
//! cost of those edits from a trial mesh, and the polarization as a function off the drawing thread. None of it changes the engine, takes a
//! job or replaces the mesh: no mesh or solve starts without a request.

use nyx_engine::estimate::Row;
use nyx_engine::Engine;
use serde_json::json;
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_str()
        .unwrap()
        .to_string()
}

fn open(name: &str) -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example(name) })).unwrap();
    e
}

/// What would change if the engine had done something: the mesh it holds, its key, the jobs it knows, its project.
fn snapshot(e: &Engine) -> (String, bool, usize, String) {
    (e.mesh_key(), e.model().is_some(), e.jobs().list().len(), e.stamp())
}

#[test]
fn the_edits_are_checked_without_being_applied() {
    let e = open("02_cantilever");
    let before = snapshot(&e);
    let (g, mut m) = (e.project().geometry.clone(), e.project().mesh.clone());
    assert!(e.check_edits(&g, &m).is_ok());
    m.h = -1.0;
    let err = e.check_edits(&g, &m).unwrap_err().to_string();
    assert!(err.contains("element size"), "{err}");
    m.h = f64::NAN;
    assert!(e.check_edits(&g, &m).is_err());
    assert_eq!(snapshot(&e), before, "checking an edit changes nothing");
    assert!(e.project().mesh.h > 0.0, "the project keeps its own element size");
}

#[test]
fn the_trial_mesh_counts_the_cells_the_real_mesh_has_and_changes_nothing() {
    let mut e = open("02_cantilever");
    let before = snapshot(&e);
    let g = e.project().geometry.clone();
    let mut m = e.project().mesh.clone();
    m.h *= 0.6;
    let got = e.size_preview(g, m.clone()).unwrap().run().unwrap();
    assert_eq!(snapshot(&e), before, "a trial mesh takes no job, stores no mesh and moves no key");
    assert_eq!(got.dim, 2);
    assert!(got.mesh_seconds < 1.0, "{got:?}");
    assert_eq!(got.estimate, None, "no solve has run on this machine: no number is invented");

    // the same edit applied for real and meshed: the trial counted exactly those cells
    e.call("set_mesh", &json!({ "h": m.h })).unwrap();
    e.call("mesh_build", &json!({})).unwrap();
    assert_eq!(got.cells, e.model().unwrap().mesh().n_cells() as u64);
}

#[test]
fn the_estimate_comes_from_the_runs_of_this_machine_and_says_so() {
    let e = open("02_cantilever");
    let (g, m) = (e.project().geometry.clone(), e.project().mesh.clone());
    let chain = e.size_preview(g.clone(), m.clone()).unwrap().run().unwrap().chain;
    assert_eq!(chain, "sms");
    // two runs of different sizes make a fit (ADR 86); fewer make none
    let row = |cells: u64, mb: u64, s: f64| Row {
        chain: "sms".into(),
        dim: 2,
        cells,
        dofs: cells * 2,
        peak_working_set: mb << 20,
        peak_commit: mb << 20,
        seconds: s,
    };
    e.history().add(row(1_000, 100, 0.5));
    assert!(e.size_preview(g.clone(), m.clone()).unwrap().run().unwrap().estimate.is_none(), "one run is no estimate");
    e.history().add(row(4_000, 400, 2.0));
    let est = e.size_preview(g, m).unwrap().run().unwrap().estimate.expect("two sizes make an estimate");
    assert_eq!(est.runs, 2);
    assert!(est.bytes > 0 && est.bytes_with_margin >= est.bytes);
}

#[test]
fn a_model_whose_mesh_is_not_known_to_be_cheap_is_not_previewed() {
    let mut e = Engine::new();
    e.call("project_new", &json!({ "title": "t", "dim": 3 })).unwrap();
    let (g, m) = (e.project().geometry.clone(), e.project().mesh.clone());
    // an empty feature list is refused by the check first; either reason keeps the mesher from being called
    assert!(e.size_preview(g, m).is_err());

    // a 2D section swept along y is cheap (milliseconds, plan P0) and previews as a 3D mesh
    let e = open("05_block_cantilever");
    let (g, m) = (e.project().geometry.clone(), e.project().mesh.clone());
    assert!(m.sweep.is_some());
    let got = e.size_preview(g, m).unwrap().run().unwrap();
    assert_eq!(got.dim, 3);
}

#[test]
fn a_bad_edit_is_refused_before_any_mesher_runs() {
    let e = open("02_cantilever");
    let (g, mut m) = (e.project().geometry.clone(), e.project().mesh.clone());
    m.h = 0.0;
    let err = e.size_preview(g, m).err().expect("a zero element size is refused").to_string();
    assert!(err.contains("element size"), "{err}");
}

#[test]
fn the_polarization_runs_off_the_engine_and_equals_the_engine_command() {
    let mut e = open("06_wr90_waveguide");
    assert!(e.polarization_task("ems").is_none(), "no result yet");
    e.call("solve_run", &json!({})).unwrap();
    let (result, task) = e.polarization_task("ems").expect("a mode study has a polarization");
    assert!(std::sync::Arc::ptr_eq(&result, e.result("ems").unwrap().0), "the task says which result it is for");
    let before = snapshot(&e);
    let off_thread = std::thread::spawn(task).join().unwrap().unwrap();
    assert_eq!(snapshot(&e), before);
    let on_engine = e.polarization("ems").unwrap();
    assert_eq!(off_thread.csv(), on_engine.csv(), "the same function of the same inputs, bit for bit");
}
