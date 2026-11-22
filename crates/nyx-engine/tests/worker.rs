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

fn build_100k_tet_dfs_engine() -> Engine {
    let mut e = Engine::new();
    let mut call = |cmd: &str, args: Value| e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}: {err}"));
    // Unit cube with lattice h = 0.042 produces > 100,000 tetrahedra cleanly and solves P1 electrostatics in < 0.5 s
    call("project_new", json!({ "title": "qcap_100k_tets", "dim": 3, "h": 0.042 }));
    call("add_box", json!({ "name": "diel", "min": [0.0, 0.0, 0.0], "max": [1.0, 1.0, 1.0] }));
    call("add_box", json!({ "name": "pad", "min": [0.35, 0.35, 0.45], "max": [0.65, 0.65, 0.55] }));
    call(
        "add_overlay",
        json!({
            "name": "ov",
            "regions": [
                { "name": "substrate", "input": "diel", "kind": "Material" },
                { "name": "pad", "input": "pad", "kind": "Hole" }
            ]
        }),
    );
    call("add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    call("set_material", json!({ "physics": "dfs", "region": "substrate", "material": { "eps_r": 11.45 } }));
    call(
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer & in_box([0.30, 0.30, 0.40], [0.70, 0.70, 0.60])",
            "kind": "terminal",
            "params": { "name": "pad" }
        }),
    );
    call(
        "add_condition",
        json!({
            "physics": "dfs",
            "on": "outer - in_box([0.30, 0.30, 0.40], [0.70, 0.70, 0.60])",
            "kind": "ground",
            "params": {}
        }),
    );
    call("set_physics_params", json!({ "physics": "dfs", "params": { "solver": "iterative", "tol": 1.0e-10 } }));
    e.set_worker(worker(), false);
    e
}

#[test]
fn isolated_worker_save_fields_false_runs_post_analyses_and_keeps_100k_tet_results_h5_under_256_kib() {
    // 1. Upfront validation rejects non-read-only or malformed post_analyses before meshing/solving
    let mut bad_eng = open("01_capacitor");
    for bad_pa in [
        json!([{ "cmd": "set_param", "args": { "name": "d", "value": 1e-3 } }]),
        json!([{ "cmd": "solve_run", "args": {} }]),
        json!([{ "analysis": "no_such_analysis" }]),
        json!([{ "analysis": "capacitance_matrix", "cmd": "integrate" }]),
        json!([{}]),
    ] {
        assert!(
            bad_eng
                .call("solve_run", &json!({ "isolate": true, "save_fields": false, "post_analyses": bad_pa }))
                .is_err(),
            "expected upfront rejection for invalid post_analyses {bad_pa}"
        );
    }

    // 2. >100k-tet 3D solve in isolated worker + workspace store with save_fields: false + both shorthand & cmd post_analyses
    let ws = std::env::temp_dir().join(format!("nyx-worker-sf-false-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&ws);

    let post_reqs = json!([
        { "analysis": "capacitance_matrix", "label": "cap_shorthand" },
        { "cmd": "analyse", "args": { "analysis": "capacitance_matrix", "physics": "dfs" }, "label": "cap_cmd" },
        { "cmd": "integrate", "args": { "expr": "eps_r * (Ex*Ex + Ey*Ey + Ez*Ez)" }, "label": "field_energy" }
    ]);

    let mut iso = build_100k_tet_dfs_engine();
    iso.set_workspace(ws.clone(), true);
    let reply_iso = iso
        .call(
            "solve_run",
            &json!({
                "isolate": true,
                "save_fields": false,
                "post_analyses": post_reqs,
            }),
        )
        .unwrap();

    let n_tets = iso.model().unwrap().mesh().n_cells();
    assert!(n_tets >= 100_000, "expected >= 100,000 tets for large-mesh size check, got {n_tets}");
    assert_eq!(reply_iso["fields_saved"], json!(false), "{reply_iso}");
    assert!(reply_iso["worker"].is_object(), "{reply_iso}");

    let pa_iso = reply_iso["post_analyses"].as_array().expect("post_analyses surfaced in reply");
    assert_eq!(pa_iso.len(), 3, "{reply_iso}");
    assert_eq!(pa_iso[0]["label"], json!("cap_shorthand"));
    assert_eq!(pa_iso[0]["analysis"], json!("capacitance_matrix"));
    assert_eq!(pa_iso[1]["label"], json!("cap_cmd"));
    assert_eq!(pa_iso[0]["maxwell_matrix_F"], pa_iso[1]["maxwell_matrix_F"]);
    assert_eq!(pa_iso[2]["label"], json!("field_energy"));
    assert!(pa_iso[2]["value"].as_f64().unwrap() > 0.0, "{pa_iso:?}");

    // Verify stored result stripped fields while keeping scalars and solution values
    let res_iso = {
        let (r, stale_iso) = iso.result("dfs").unwrap();
        assert!(!stale_iso);
        r.clone()
    };
    assert!(!res_iso.scalars.is_empty());
    assert!(res_iso.solutions.iter().all(|s| s.fields.is_empty() && s.cell_fields.is_empty()));

    // Verify recorded run's results.h5 is < 256 KiB and results.xdmf was skipped!
    let runs = iso.call("store_runs", &json!({})).unwrap();
    let artifacts = runs[0]["artifacts"].as_array().unwrap();
    assert!(
        !artifacts.iter().any(|a| a["path"].as_str().unwrap().ends_with("results.xdmf")),
        "xdmf::write must be skipped when !fields_saved: {artifacts:?}"
    );
    let h5_art =
        artifacts.iter().find(|a| a["path"].as_str().unwrap().ends_with("results.h5")).expect("results.h5 artifact");
    let h5_bytes = h5_art["bytes"].as_u64().unwrap();
    assert!(h5_bytes > 0 && h5_bytes < 256 * 1024, "100k+ tet results.h5 must be < 256 KiB, got {h5_bytes} B");

    // Also verify direct worker job folder produces < 256 KiB results.h5
    let job_dir = ws.join("direct_worker_job");
    std::fs::create_dir_all(&job_dir).unwrap();
    let lib = iso.library().unwrap();
    std::fs::write(job_dir.join("project.nyx.json"), iso.project().to_json()).unwrap();
    std::fs::write(
        job_dir.join("job.json"),
        json!({
            "targets": ["dfs"],
            "materials": nyx_engine::worker::material_texts(&lib),
            "warm": false,
            "save_fields": false,
            "post_analyses": post_reqs,
        })
        .to_string(),
    )
    .unwrap();
    let mut child = std::process::Command::new(worker())
        .arg("--job")
        .arg(&job_dir)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let _stdin = child.stdin.take();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success(), "direct nyx-worker failed: {}", String::from_utf8_lossy(&out.stderr));
    let worker_h5_bytes = nyx_engine::h5io::file_bytes(&job_dir.join("results.h5"));
    assert!(
        worker_h5_bytes > 0 && worker_h5_bytes < 256 * 1024,
        "worker results.h5 for {n_tets} tets must be < 256 KiB, got {worker_h5_bytes} B"
    );
    assert!(job_dir.join("post_analyses.json").exists());

    // 3. In-process path (`isolate: false`) with `save_fields: false` behaves identically!
    let mut inproc = build_100k_tet_dfs_engine();
    let reply_inproc = inproc
        .call(
            "solve_run",
            &json!({
                "isolate": false,
                "save_fields": false,
                "post_analyses": post_reqs,
            }),
        )
        .unwrap();
    assert_eq!(reply_inproc["fields_saved"], json!(false));
    let pa_inproc = reply_inproc["post_analyses"].as_array().unwrap();
    assert_eq!(pa_inproc[0]["maxwell_matrix_F"], pa_iso[0]["maxwell_matrix_F"]);
    assert_eq!(pa_inproc[2]["value"], pa_iso[2]["value"]);
    let (res_inproc, _) = inproc.result("dfs").unwrap();
    same_numbers(&res_iso, res_inproc, "isolated vs in-process save_fields=false");

    let _ = std::fs::remove_dir_all(&ws);
}
