//! Every long operation is a job (productivity plan P2, ADR 85): a mesh, a solve, a mesh convergence check, opening a results file and writing
//! one run on a job that can be watched and cancelled, and a job started for a project that was replaced does not put its outcome into the new one.
//! Each cancel test fails when the job reaches `Done`: a cancel that did not happen cannot pass for one.

use nyx_engine::{Engine, JobState};
use serde_json::{json, Value};
use std::path::PathBuf;
use std::time::{Duration, Instant};

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn open(name: &str) -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example(name) })).unwrap();
    e
}

fn code<T>(r: &Result<T, nyx_physics::Error>) -> &'static str {
    match r {
        Ok(_) => "ok",
        Err(e) => e.code(),
    }
}

/// Waits until the job has a stage (it is really running), then asks it to stop; returns the time from the request to the outcome.
fn cancel_when_running(e: &mut Engine, id: u64) -> (Result<Value, nyx_physics::Error>, f64) {
    let t = Instant::now();
    while e.jobs().info(id).is_some_and(|i| i.stage.is_empty() && i.state == JobState::Running) {
        assert!(t.elapsed() < Duration::from_secs(30), "the job never reported a stage");
        std::thread::sleep(Duration::from_millis(2));
    }
    std::thread::sleep(Duration::from_millis(80));
    assert_eq!(
        e.jobs().info(id).unwrap().state,
        JobState::Running,
        "the job was over before the cancel: the test needs a larger model"
    );
    assert!(e.jobs().cancel(id));
    let asked = Instant::now();
    assert!(e.jobs().info(id).unwrap().cancel_requested, "the job info says a cancel was requested");
    let r = e.wait(id);
    (r, asked.elapsed().as_secs_f64())
}

#[test]
fn an_adaptive_mesh_job_is_cancelled_between_its_passes() {
    let mut e = open("15_soi_wire");
    e.call("set_mesh", &json!({ "h": 2e-8, "adapt_scale": 1.0 })).unwrap(); // about 52,000 cells, 0.8 s
    let id = e.start_mesh().unwrap();
    let (r, took) = cancel_when_running(&mut e, id);
    assert_eq!(code(&r), "cancelled", "{r:?}");
    assert!(took < 2.0, "the cancel took {took} s");
    assert!(e.model().is_none(), "a cancelled mesh job leaves the engine as it was");
    assert_eq!(e.jobs().info(id).unwrap().state, JobState::Cancelled);
}

#[test]
fn a_3d_mesh_job_is_cancelled_between_its_phases() {
    let mut e = Engine::new();
    e.call("project_new", &json!({ "dim": 3, "h": 0.01 })).unwrap(); // about 12 million tetrahedra, about 2 s
    e.call("add_box", &json!({ "name": "b", "min": [0.0, 0.0, 0.0], "max": [1.0, 1.0, 1.0] })).unwrap();
    e.call("add_overlay", &json!({ "name": "regions", "regions": [{ "name": "b", "input": "b" }] })).unwrap();
    let id = e.start_mesh().unwrap();
    let (r, took) = cancel_when_running(&mut e, id);
    assert_eq!(code(&r), "cancelled", "{r:?}");
    assert!(took < 3.0, "the cancel took {took} s");
}

#[test]
fn a_convergence_check_is_a_job_with_two_stages_and_a_cancel() {
    let mut e = open("14_box_cavity");
    e.call("set_mesh", &json!({ "h": 0.003, "sweep": [{ "domain_offset": 0, "n": 8, "thickness": 0.03 }] })).unwrap(); // a solve of several seconds
    let id = e.start_convergence("ems", 0.6, Some(1e-2)).unwrap();
    assert!(
        e.jobs().list().iter().any(|j| j.id == id && j.label.contains("convergence")),
        "the check is in the job list"
    );
    let (r, took) = cancel_when_running(&mut e, id);
    assert_eq!(code(&r), "cancelled", "{r:?}");
    assert!(took < 5.0, "the cancel took {took} s (the stage it was in cannot be interrupted)");
    assert!(
        e.jobs().log(id).iter().any(|l| l.contains("convergence check")),
        "the stage is in the log: {:?}",
        e.jobs().log(id)
    );
}

#[test]
fn a_convergence_job_gives_the_answer_of_the_command_with_wait() {
    let mut e = open("15_soi_wire");
    let direct = e.call("convergence_check", &json!({ "physics": "ems", "tolerance": 1e-2 })).unwrap();
    let started = e.call("convergence_check", &json!({ "physics": "ems", "tolerance": 1e-2, "wait": false })).unwrap();
    let job = started["job"].as_u64().expect("a job id");
    let via_job = e.wait(job).unwrap();
    assert_eq!(direct["rows"], via_job["rows"]);
    assert_eq!(direct["converged"], via_job["converged"]);
    // two stages, in the log of the job
    let log = e.jobs().log(job).join("\n");
    assert!(log.contains("present mesh") && log.contains("finer mesh"), "{log}");
    // the project, its mesh and its results were not touched
    assert!(e.model().is_none() && e.result("ems").is_none());
}

#[test]
fn a_results_file_is_written_and_opened_on_jobs_and_the_file_is_the_snapshot() {
    let dir = std::env::temp_dir().join(format!("nyx-jobs-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("r.h5").to_string_lossy().into_owned();
    let mut e = open("15_soi_wire");
    e.call("solve_run", &json!({})).unwrap();
    let id = e.call("export_results", &json!({ "path": path, "wait": false })).unwrap()["job"].as_u64().unwrap();
    // the engine goes on: this edit is not in the file
    e.call("set_param", &json!({ "name": "w", "value": 7.0e-7 })).unwrap();
    let r = e.wait(id).unwrap();
    assert!(r["bytes"].as_u64().unwrap() > 0);
    assert!(
        std::path::Path::new(&path).exists() && !std::path::Path::new(&format!("{path}.partial")).exists(),
        "the file appears under its name only when complete"
    );
    // opening: the engine is untouched until the job is applied
    let mut f = open("04_thermal_bar");
    let title = f.project().title.clone();
    let open_id = f.call("open_results", &json!({ "path": path, "wait": false })).unwrap()["job"].as_u64().unwrap();
    assert_eq!(f.project().title, title, "nothing changed yet");
    let o = f.wait(open_id).unwrap();
    assert_eq!(o["results"][0]["physics"], json!("ems"));
    assert_eq!(f.project().spec("ems").unwrap().physics, "ems");
    let width = f.project().geometry.params.iter().find(|p| p.name == "w").map(|p| p.value.clone()).expect("the width");
    let shown = format!("{width:?}");
    assert!(
        !(shown.contains("7e-7") || shown.contains("7.0e-7")),
        "the file holds the project of the moment of the request: {shown}"
    );
    assert!(!f.result("ems").unwrap().1, "opened results are not stale");
    // a file that is not a results file fails on the job and changes nothing
    let bad = dir.join("bad.h5");
    std::fs::write(&bad, b"not hdf5").unwrap();
    let t = f.project().title.clone();
    let id = f.start_open_results(bad.to_str().unwrap()).unwrap();
    assert!(f.wait(id).is_err());
    assert_eq!(f.project().title, t);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_job_started_for_a_project_that_was_replaced_is_dropped() {
    let mut e = open("15_soi_wire");
    let id = e.start_solve(None).unwrap();
    e.call("project_new", &json!({ "dim": 2 })).unwrap(); // another project now
    let r = e.wait(id);
    assert_eq!(code(&r), "invalid", "{r:?}");
    assert!(e.result("ems").is_none(), "the old project's result did not land in the new one");
    let id = e.start_mesh();
    assert!(id.is_err() || e.model().is_none());
}

/// The window calls no command that waits for the engine on its own thread: its long operations are jobs it polls. The tests of the window and
/// the headless `batch.rs` may wait.
#[test]
fn the_window_never_blocks_on_the_engine() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../app/src");
    let forbidden = [
        ".mesh_now(",
        ".solve_now(",
        ".engine.wait(",
        ".open_results(",
        ".export_results(",
        "\"open_results\"",
        "\"convergence_check\"",
        "\"export_results\"",
    ];
    let mut seen = 0;
    for entry in std::fs::read_dir(&dir).expect("crates/app/src") {
        let path = entry.unwrap().path();
        if path.extension().is_none_or(|e| e != "rs") || path.file_name().is_some_and(|n| n == "batch.rs") {
            continue;
        }
        let text = std::fs::read_to_string(&path).unwrap();
        let code = text.split("#[cfg(test)]").next().unwrap_or(&text);
        for (n, line) in code.lines().enumerate() {
            let line = line.split("//").next().unwrap_or(line);
            for f in forbidden {
                assert!(
                    !line.contains(f),
                    "{}:{}: the window must not block on the engine ({f}); start a job and poll it",
                    path.display(),
                    n + 1
                );
            }
        }
        seen += 1;
    }
    assert!(seen > 5, "the window's sources were found");
}

#[test]
fn job_status_says_what_a_waiting_person_needs() {
    let mut e = open("14_box_cavity");
    e.call("set_mesh", &json!({ "h": 0.003, "sweep": [{ "domain_offset": 0, "n": 8, "thickness": 0.03 }] })).unwrap();
    let id = e.start_solve(None).unwrap();
    let t = Instant::now();
    // wait for the factorisation, which cannot be interrupted
    loop {
        let i = e.jobs().info(id).unwrap();
        if i.stage.contains("factorising the shifted") {
            break;
        }
        assert!(t.elapsed() < Duration::from_secs(30), "no factorisation seen: {i:?}");
        std::thread::sleep(Duration::from_millis(5));
    }
    assert!(e.jobs().cancel(id));
    std::thread::sleep(Duration::from_millis(300));
    let s = e.call("job_status", &json!({ "id": id })).unwrap();
    assert_eq!(s["cancel_requested"], json!(true), "{s}");
    assert!(s["stage"].as_str().unwrap().contains("factorising"), "still in the factorisation: {s}");
    assert!(s["stage_seconds"].as_f64().unwrap() >= 0.25, "the time in the stage is reported: {s}");
    assert_eq!(code(&e.wait(id)), "cancelled");
}
