//! A worker that dies (here on purpose, by `NYX_WORKER_DRILL=abort`, which makes `nyx-worker` abort a moment after it starts, as the allocator does
//! when memory runs out) fails its job and nothing else: the engine, its results and the next job are intact (productivity plan P4, ADR 86).
//! The test sets an environment variable, so it is the only test of this binary.

use nyx_engine::{Engine, JobState};
use serde_json::json;
use std::path::PathBuf;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

#[test]
fn a_worker_that_crashes_fails_its_job_and_the_next_job_runs() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("02_cantilever") })).unwrap();
    // a result the crash must not touch
    e.call("solve_run", &json!({})).unwrap();
    let before = e.result("sms").unwrap().0.clone();
    e.set_worker(PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker")), false);
    std::env::set_var("NYX_WORKER_DRILL", "abort");
    let a = e.call("solve_run", &json!({ "isolate": true, "wait": false })).unwrap()["job"].as_u64().unwrap();
    // queued behind it, in this process
    let b = e.call("queue_submit", &json!({ "isolate": false })).unwrap()["id"].as_u64().unwrap();
    let err = e.wait(a).unwrap_err().to_string();
    assert!(err.contains("worker process stopped without an answer"), "{err}");
    assert!(err.contains("Nothing else was affected"), "{err}");
    assert!(matches!(e.jobs().info(a).unwrap().state, JobState::Failed(_)));
    // the next job starts by itself and finishes; the earlier result is still there
    let ok = e.wait(b).unwrap();
    assert_eq!(ok["state"], "done");
    assert!(
        std::sync::Arc::ptr_eq(&before, e.result("sms").unwrap().0) || e.result("sms").unwrap().0.dofs == before.dofs
    );
    std::env::remove_var("NYX_WORKER_DRILL");
}
