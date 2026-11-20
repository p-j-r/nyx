// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Solving in a worker process (productivity plan P4, ADR 86): the numbers are those of the same solve in this process, to the bit, over the
//! examples and a coupled chain; the stage log of the child is the job's log; a cancel kills the child at once; the library and the mesh the
//! child used are the parent's; what cannot be isolated says so. The crash of a worker is in `worker_crash.rs` (it sets an environment variable,
//! so it is a test binary of its own).

use nyx_engine::{Engine, JobState, MeshSettings, MeshedModel, Mesher};
use nyx_physics::{Error, SolveResult};
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Duration, Instant};

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

fn open(name: &str) -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example(name) })).unwrap_or_else(|err| panic!("{name}: {err}"));
    e.set_worker(worker(), false);
    e
}

/// Everything of a result that is a number or a name (not the seconds it took, not the log, which carries timings).
fn same_numbers(a: &SolveResult, b: &SolveResult, what: &str) {
    assert_eq!(a.physics, b.physics, "{what}");
    assert_eq!(a.study, b.study, "{what}");
    assert_eq!(a.dofs, b.dofs, "{what}: unknowns");
    assert_eq!(a.scalars, b.scalars, "{what}: scalars");
    assert_eq!(a.solutions, b.solutions, "{what}: solutions (fields, eigenvalues)");
}

/// Solves `name` in this process and in a worker and compares every physics of the chain.
fn compare(name: &str, physics: Option<Value>) -> usize {
    let mut args_here = json!({ "isolate": false });
    let mut args_there = json!({ "isolate": true });
    if let Some(p) = physics {
        args_here["physics"] = p.clone();
        args_there["physics"] = p;
    }
    let mut here = open(name);
    here.call("solve_run", &args_here).unwrap_or_else(|e| panic!("{name} in process: {e}"));
    let mut there = open(name);
    let reply = there.call("solve_run", &args_there).unwrap_or_else(|e| panic!("{name} in a worker: {e}"));
    assert!(reply["worker"].is_object(), "{name}: the reply says a worker ran it: {reply}");
    let mut n = 0;
    for spec in here.project().physics.clone() {
        let Some((a, _)) = here.result(&spec.physics) else { continue };
        let (b, stale) = there
            .result(&spec.physics)
            .unwrap_or_else(|| panic!("{name}: the worker's run has no result for {}", spec.physics));
        assert!(!stale, "{name}: {} is current", spec.physics);
        same_numbers(a, b, &format!("{name} {}", spec.physics));
        n += 1;
    }
    assert!(n > 0, "{name}: something was solved");
    // the mesh of the parent is the mesh the results are on
    assert_eq!(here.model().unwrap().solver_mesh().n_cells(), there.model().unwrap().solver_mesh().n_cells());
    n
}

#[test]
fn a_worker_gives_the_numbers_of_the_same_solve_in_this_process_to_the_bit() {
    compare("01_capacitor", None);
    compare("02_cantilever", None);
    compare("06_wr90_waveguide", None);
    compare("15_soi_wire", None);
}

#[test]
fn a_coupled_chain_in_a_worker_equals_the_chain_in_this_process() {
    assert_eq!(
        compare("04_thermal_bar", Some(json!(["dfs", "sms"]))),
        2,
        "both physics of the chain were solved and compared"
    );
}

#[test]
fn the_stage_log_of_the_worker_is_the_log_of_the_job_and_the_run_says_how_much_memory_it_used() {
    let mut e = open("14_box_cavity");
    let id = e.call("solve_run", &json!({ "isolate": true, "wait": false })).unwrap()["job"].as_u64().unwrap();
    // while it runs the job reports the worker's stage
    let start = Instant::now();
    let mut stage = String::new();
    while start.elapsed() < Duration::from_secs(60) {
        let info = e.jobs().info(id).unwrap();
        if !info.stage.is_empty() && info.state == JobState::Running {
            stage = info.stage;
            break;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    assert!(!stage.is_empty(), "the job shows the stage of the worker while it runs");
    let reply = e.wait(id).unwrap();
    let w = &reply["worker"];
    assert!(w["peak_working_set"].as_u64().unwrap_or(0) > 10_000_000, "the worker reports its peak memory: {reply}");
    let log = e.call("queue_get", &json!({ "id": id })).unwrap()["log"].to_string();
    assert!(log.contains("worker process: peak memory"), "{log}");
    let full = e.jobs().log(id);
    assert!(full.len() > 2, "the worker stages are in the log: {full:?}");
}

#[test]
fn a_cancel_kills_the_worker_at_once_whatever_stage_it_is_in() {
    // 39k unknowns: the shifted matrix takes seconds to factorise, and an in-process cancel waits for it
    let mut e = open("14_box_cavity");
    e.call("set_mesh", &json!({ "h": 0.0028 })).unwrap();
    let id = e.call("solve_run", &json!({ "isolate": true, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let start = Instant::now();
    while !e.jobs().info(id).unwrap().stage.contains("factorising the shifted")
        && start.elapsed() < Duration::from_secs(60)
    {
        std::thread::sleep(Duration::from_millis(5));
    }
    assert!(
        e.jobs().info(id).unwrap().stage.contains("factorising the shifted"),
        "the worker reached the factorisation"
    );
    let asked = Instant::now();
    assert!(e.jobs().cancel(id));
    let err = e.wait(id).unwrap_err();
    assert!(matches!(err, Error::Cancelled), "{err}");
    // an in-process cancel would wait for the factorisation to end; a kill does not
    assert!(asked.elapsed() < Duration::from_millis(500), "stopped after {:?}", asked.elapsed());
    assert_eq!(e.jobs().info(id).unwrap().state, JobState::Cancelled);
    // the engine goes on: the next solve runs
    e.call("solve_run", &json!({ "isolate": false, "physics": "ems" })).unwrap();
}

#[test]
fn what_cannot_be_isolated_says_so() {
    // no worker set
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("02_cantilever") })).unwrap();
    let err = e.call("solve_run", &json!({ "isolate": true })).unwrap_err().to_string();
    assert!(err.contains("no worker") || err.contains("has none"), "{err}");
    // an engine with a mesher of its own
    struct Plain(nyx_engine::FemcoreMesher);
    impl Mesher for Plain {
        fn mesh(
            &self,
            g: &geom::model::FeatureList,
            s: &MeshSettings,
            z: Option<&nyx_physics::Sizing>,
        ) -> Result<MeshedModel, Error> {
            self.0.mesh(g, s, z)
        }
    }
    let mut e = open("02_cantilever");
    e.set_mesher(Arc::new(Plain(nyx_engine::FemcoreMesher)));
    let err = e.call("solve_run", &json!({ "isolate": true })).unwrap_err().to_string();
    assert!(err.contains("registry or mesher"), "{err}");
    // a worker program that is not there: the solve fails with a reason and the engine is intact
    let mut e = open("02_cantilever");
    e.set_worker(PathBuf::from("no-such-worker-program"), false);
    let err = e.call("solve_run", &json!({ "isolate": true })).unwrap_err().to_string();
    assert!(err.contains("could not start the worker"), "{err}");
    e.call("solve_run", &json!({})).unwrap();
}

#[test]
fn a_worker_by_default_isolates_every_solve_unless_it_says_not() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("01_capacitor") })).unwrap();
    e.set_worker(worker(), true);
    assert!(e.call("solve_run", &json!({})).unwrap()["worker"].is_object());
    assert!(e.call("solve_run", &json!({ "isolate": false })).unwrap()["worker"].is_null());
}

#[test]
fn the_materials_a_worker_is_given_read_back_to_the_same_library() {
    // the library a worker rebuilds from text is the library of the parent, material by material (a loss here would be a different number)
    let e = Engine::new();
    let lib = e.library().unwrap();
    let texts = nyx_engine::worker::material_texts(&lib);
    let back = femcore::material::Library::default().with_texts(&texts).unwrap();
    assert_eq!(lib.names(), back.names());
    for n in lib.names() {
        assert_eq!(lib.get(&n).unwrap().to_toml(), back.get(&n).unwrap().to_toml(), "{n}");
        assert_eq!(format!("{:?}", lib.get(&n).unwrap()), format!("{:?}", back.get(&n).unwrap()), "{n}");
    }
    assert!(lib.len() > 3, "the library was found ({} materials)", lib.len());
}

#[test]
fn sweep_points_can_run_in_workers_with_the_numbers_of_in_process_points() {
    let run = |isolate: bool| -> Value {
        let mut e = open("15_soi_wire");
        let id = e
            .call(
                "sweep_run",
                &json!({ "series": { "cmd": "set_param", "args": { "name": "w" }, "field": "value", "values": [4.5e-7, 5.5e-7] }, "isolate": isolate }),
            );
        let _ = id.expect("sweep starts");
        let t = e.call("sweep_table", &json!({ "sweep": 1, "wait_s": 120.0 })).unwrap();
        t["points"].clone()
    };
    let (a, b) = (run(false), run(true));
    let strip = |v: &Value| -> Vec<Value> {
        v.as_array()
            .unwrap()
            .iter()
            .map(|r| json!({ "numbers": r["numbers"], "unknowns": r["unknowns"], "state": r["state"] }))
            .collect()
    };
    assert_eq!(strip(&a), strip(&b));
    assert!(strip(&b).iter().all(|r| r["state"] == "done"), "{b}");
}
