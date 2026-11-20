// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Results are stale per physics (productivity plan P1, ADR 84): the key of a result is the problem key of its physics (its own problem, the project
//! materials, the mesh key and the keys of the physics it takes fields from), so an edit stales what it feeds and nothing else, and a job that
//! finishes late cannot replace a result that is current.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::PathBuf;

fn open(name: &str) -> Engine {
    let mut e = Engine::new();
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples").join(format!("{name}.nyx.json"));
    e.call("project_open", &json!({ "path": path.to_str().unwrap() })).unwrap_or_else(|err| panic!("{name}: {err}"));
    e
}

fn stale(e: &Engine, p: &str) -> bool {
    e.result(p).unwrap_or_else(|| panic!("no result for {p}")).1
}

fn neff(e: &Engine) -> String {
    let r = e.result("ems").expect("an ems result").0.clone();
    format!("{:?}", r.scalars.iter().find(|(n, _)| n == "neff 1").map(|(_, v)| v.clone()))
}

/// The order of `physics` the other way round (1 or 2): a change of its problem.
fn toggled(e: &Engine, p: &str) -> u8 {
    if e.project().spec(p).unwrap().order == 1 {
        2
    } else {
        1
    }
}

#[test]
fn a_coupled_physics_follows_its_import_and_the_import_does_not_follow_it() {
    // 04: `dfs` heat, `sms` takes the temperature from it
    let mut e = open("04_thermal_bar");
    e.call("solve_run", &json!({ "physics": ["dfs", "sms"] })).unwrap();
    assert!(!stale(&e, "dfs") && !stale(&e, "sms"));
    // a change in the solid problem leaves the heat result alone
    let (sms_order, dfs_order) = (e.project().spec("sms").unwrap().order, e.project().spec("dfs").unwrap().order);
    let order = toggled(&e, "sms");
    e.call("set_study", &json!({ "physics": "sms", "kind": "static", "order": order })).unwrap();
    assert!(stale(&e, "sms") && !stale(&e, "dfs"), "sms edited: only sms is stale");
    // staleness is a property of the content: putting the edit back makes the result current again
    e.call("set_study", &json!({ "physics": "sms", "kind": "static", "order": sms_order })).unwrap();
    assert!(!stale(&e, "sms") && !stale(&e, "dfs"));
    // a change in the heat problem stales the solid result that was made from it
    e.call("set_study", &json!({ "physics": "dfs", "kind": "heat_stationary", "order": 3 - dfs_order })).unwrap();
    assert!(stale(&e, "dfs") && stale(&e, "sms"), "dfs edited: the solid result took its temperature");
}

#[test]
fn a_geometry_edit_stales_every_result() {
    let mut e = open("04_thermal_bar");
    e.call("solve_run", &json!({ "physics": ["dfs", "sms"] })).unwrap();
    let info = e.call("project_info", &json!({})).unwrap();
    assert!(info["results"].as_array().unwrap().iter().all(|r| r["stale"] == false));
    e.call("set_param", &json!({ "name": "an_unused_length", "value": 0.123 })).unwrap();
    assert!(stale(&e, "dfs") && stale(&e, "sms"));
}

#[test]
fn choosing_the_active_physics_changes_no_key() {
    let mut e = open("04_thermal_bar");
    e.call("solve_run", &json!({ "physics": ["dfs", "sms"] })).unwrap();
    let before = e.problem_keys();
    let other = if e.project().active.as_deref() == Some("sms") { "dfs" } else { "sms" };
    let whole = e.stamp();
    e.call("set_active", &json!({ "physics": other })).unwrap();
    assert_ne!(whole, e.stamp(), "the whole-project hash did change: it is what the old staleness used");
    assert_eq!(before, e.problem_keys());
    assert!(!stale(&e, "dfs") && !stale(&e, "sms"));
}

#[test]
fn a_wavelength_is_in_the_mesh_key_only_when_the_mesh_is_sized_from_the_physics() {
    let mut e = open("15_soi_wire");
    let freq =
        |f: f64| json!({ "physics": "ems", "kind": "modes", "params": { "freq": f, "n_modes": 2, "n_target": 2.75 } });
    e.call("set_mesh", &json!({ "physics_sizing": false })).unwrap();
    let key = e.mesh_key();
    e.call("set_study", &freq(150.0e12)).unwrap();
    assert_eq!(key, e.mesh_key(), "not sized from the physics: the mesh does not depend on the wavelength");
    e.call("set_mesh", &json!({ "physics_sizing": true })).unwrap();
    let key = e.mesh_key();
    e.call("set_study", &freq(250.0e12)).unwrap();
    assert_ne!(
        key,
        e.mesh_key(),
        "sized from the physics: a wavelength changes the mesh, so every result on it is stale"
    );
}

#[test]
fn a_late_job_cannot_replace_a_result_that_is_current() {
    let mut e = open("15_soi_wire");
    let a = e.start_solve(None).unwrap(); // the wire as it is
    e.call("set_param", &json!({ "name": "w", "value": 6.0e-7 })).unwrap();
    let b = e.start_solve(None).unwrap(); // the wire after the edit: the one that is current
    assert!(b > a);
    // the newer job finishes first, the older one afterwards
    let rb = e.wait(b).unwrap();
    assert_eq!(rb["results"][0]["current"], json!(true));
    let from_b = neff(&e);
    let ra = e.wait(a).unwrap();
    assert_eq!(ra["results"][0]["current"], json!(false), "{ra}");
    assert_eq!(ra["results"][0]["stored"], json!(false), "{ra}");
    assert_eq!(neff(&e), from_b, "the stored result is still the one of the edited wire");
    assert!(!stale(&e, "ems"));
}

#[test]
fn of_two_old_jobs_the_one_started_later_is_kept() {
    let mut e = open("15_soi_wire");
    let a = e.start_solve(None).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": 5.5e-7 })).unwrap();
    let b = e.start_solve(None).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": 6.0e-7 })).unwrap(); // both jobs are old now
    e.wait(b).unwrap();
    let from_b = neff(&e);
    assert!(stale(&e, "ems"));
    let ra = e.wait(a).unwrap();
    assert_eq!(ra["results"][0]["stored"], json!(false), "the older job does not replace the newer one");
    assert_eq!(neff(&e), from_b);
    // the other order: the older job first, then the newer one replaces it
    let mut e = open("15_soi_wire");
    let a = e.start_solve(None).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": 5.5e-7 })).unwrap();
    let b = e.start_solve(None).unwrap();
    e.call("set_param", &json!({ "name": "w", "value": 6.0e-7 })).unwrap();
    e.wait(a).unwrap();
    let from_a = neff(&e);
    let rb = e.wait(b).unwrap();
    assert_eq!(rb["results"][0]["stored"], json!(true));
    assert_ne!(neff(&e), from_a);
}

#[test]
fn the_replies_of_the_result_commands_say_whether_the_result_is_stale() {
    let mut e = open("15_soi_wire");
    e.call("solve_run", &json!({})).unwrap();
    let info = e.call("result_info", &json!({})).unwrap();
    assert_eq!(info[0]["stale"], json!(false));
    assert!(info[0]["age_s"].as_f64().unwrap() >= 0.0);
    let field = e.call("result_field", &json!({ "physics": "ems", "field": "|E|" })).unwrap();
    assert_eq!(field["stale"], json!(false));
    let an = e.call("analyse", &json!({ "analysis": "polarization", "physics": "ems" })).unwrap();
    assert_eq!(an["stale"], json!(false));
    assert_eq!(e.call("list_analyses", &json!({})).unwrap()[0]["stale"]["ems"], json!(false));
    // none of them shows a key
    let all = format!("{info}{field}{an}");
    assert!(!all.contains("b3:"), "a cache key must not leave the process (ADR 64)");
    e.call("set_param", &json!({ "name": "w", "value": 6.0e-7 })).unwrap();
    assert_eq!(e.call("result_info", &json!({})).unwrap()[0]["stale"], json!(true));
    assert_eq!(e.call("result_field", &json!({ "physics": "ems", "field": "|E|" })).unwrap()["stale"], json!(true));
    assert_eq!(
        e.call("analyse", &json!({ "analysis": "polarization", "physics": "ems" })).unwrap()["stale"],
        json!(true)
    );
    assert_eq!(e.call("list_analyses", &json!({})).unwrap()[0]["stale"]["ems"], json!(true));
    assert_eq!(e.call("project_info", &json!({})).unwrap()["results"][0]["stale"], Value::Bool(true));
}
