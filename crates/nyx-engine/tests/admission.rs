//! Estimates and admission (productivity plan Q2, ADR 86): what solves cost is learned from the runs of workers; a solve whose estimate does not
//! fit the memory free now waits and says why, one the machine could never hold is refused unless forced, a cancel ends the wait, and a
//! solve with no estimate (fewer than two sizes seen) starts as before. The machine's memory is played by a probe, so nothing here depends on
//! what the test machine has free.

use nyx_engine::estimate::Row;
use nyx_engine::{Engine, JobState};
use serde_json::json;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
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

/// An engine that has solved `name` once (so it holds the mesh), with two history rows around the size of that mesh whose memory is `gb_small`
/// and `gb_large` times 1e9 bytes: the estimate for this project is then in between, whatever the machine.
fn engine_with_history(name: &str, gb_small: f64, gb_large: f64) -> (Engine, String, u8, u64) {
    let mut e = open(name);
    e.call("solve_run", &json!({ "isolate": false })).unwrap();
    let cells = e.model().unwrap().mesh().n_cells() as u64;
    let dim = e.model().unwrap().mesh().dim() as u8;
    let chain = e.project().default_targets().join(",");
    for (c, gb) in [(cells / 2, gb_small), (cells * 2, gb_large)] {
        e.history().add(Row {
            chain: chain.clone(),
            dim,
            cells: c,
            dofs: c * 3,
            peak_working_set: (gb * 1e9) as u64,
            peak_commit: (gb * 1.2e9) as u64,
            seconds: 10.0,
        });
    }
    (e, chain, dim, cells)
}

fn machine(total_gb: f64, free_gb: f64) -> (Arc<AtomicU64>, nyx_engine::estimate::MemoryProbe) {
    let free = Arc::new(AtomicU64::new((free_gb * 1e9) as u64));
    let f = free.clone();
    let total = (total_gb * 1e9) as u64;
    (free, Arc::new(move || Some((total, f.load(Ordering::SeqCst)))))
}

#[test]
fn a_solve_with_no_estimate_starts_at_once() {
    // no history, or one size only: nothing is invented
    let mut e = open("02_cantilever");
    let (_, probe) = machine(16.0, 0.001);
    e.set_memory_probe(probe);
    let r = e.call("solve_run", &json!({ "isolate": false })).unwrap();
    assert_eq!(r["state"], "done");
    let cells = e.model().unwrap().mesh().n_cells() as u64;
    let dim = e.model().unwrap().mesh().dim() as u8;
    let chain = e.project().default_targets().join(",");
    e.history().add(Row {
        chain,
        dim,
        cells,
        dofs: 1,
        peak_working_set: 900_000_000_000,
        peak_commit: 1,
        seconds: 1.0,
    });
    // still one size
    assert_eq!(e.call("solve_run", &json!({ "isolate": false })).unwrap()["state"], "done");
}

#[test]
fn a_solve_waits_for_memory_says_why_and_starts_when_it_is_free() {
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    let (free, probe) = machine(2000.0, 4.0);
    e.set_memory_probe(probe);
    let id = e.call("solve_run", &json!({ "isolate": false, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let start = Instant::now();
    let mut stage = String::new();
    while start.elapsed() < Duration::from_secs(20) {
        stage = e.jobs().info(id).unwrap().stage;
        if stage.starts_with("waiting for memory") {
            break;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    assert!(stage.starts_with("waiting for memory: needs about") && stage.contains("4.0 GB free"), "{stage}");
    // the entry shows the estimate and its basis
    let entry = e.call("queue_get", &json!({ "id": id })).unwrap();
    assert_eq!(entry["state"], "running");
    assert_eq!(entry["estimate"]["runs"], 2, "{entry}");
    assert!(entry["estimate"]["bytes"].as_u64().unwrap() > 50_000_000_000, "{entry}");
    // memory frees up: it starts and finishes
    free.store(1_900_000_000_000, Ordering::SeqCst);
    let r = e.wait(id).unwrap();
    assert_eq!(r["state"], "done");
    assert_eq!(e.jobs().info(id).unwrap().state, JobState::Done);
}

#[test]
fn a_cancel_ends_the_wait() {
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    let (_, probe) = machine(2000.0, 4.0);
    e.set_memory_probe(probe);
    let id = e.call("solve_run", &json!({ "isolate": false, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let start = Instant::now();
    while !e.jobs().info(id).unwrap().stage.starts_with("waiting for memory")
        && start.elapsed() < Duration::from_secs(20)
    {
        std::thread::sleep(Duration::from_millis(20));
    }
    let asked = Instant::now();
    assert!(e.jobs().cancel(id));
    assert!(matches!(e.wait(id).unwrap_err(), nyx_physics::Error::Cancelled));
    assert!(asked.elapsed() < Duration::from_secs(2));
}

#[test]
fn what_the_machine_cannot_hold_is_refused_unless_forced() {
    let (mut e, ..) = engine_with_history("02_cantilever", 100.0, 400.0);
    // the estimate is about 100 to 400 GB; this machine has 64 GB
    let (_, probe) = machine(64.0, 64.0);
    e.set_memory_probe(probe);
    let err = e.call("solve_run", &json!({ "isolate": false })).unwrap_err().to_string();
    assert!(err.contains("more than this machine can hold") && err.contains("force"), "{err}");
    assert!(err.contains("from 2 earlier runs"), "{err}");
    // forced, it runs (here the estimate is wrong by design; the test only shows the way through)
    let r = e.call("solve_run", &json!({ "isolate": false, "force": true })).unwrap();
    assert_eq!(r["state"], "done");
}

#[test]
fn a_worker_run_teaches_the_next_estimate_and_the_history_survives_in_a_file() {
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"));
    let file = std::env::temp_dir().join(format!("nyx-admission-{}", std::process::id())).join("solve_history.jsonl");
    let mut e = open("14_box_cavity");
    e.set_worker(exe.clone(), false);
    e.set_history_file(&file);
    // two sizes
    for h in [0.006, 0.0045] {
        e.call("set_mesh", &json!({ "h": h })).unwrap();
        e.call("solve_run", &json!({ "isolate": true })).unwrap();
    }
    assert_eq!(e.history().rows().len(), 2);
    // a third size gets an estimate from them, in the entry of its job
    e.call("set_mesh", &json!({ "h": 0.005 })).unwrap();
    let id = e.call("solve_run", &json!({ "isolate": true, "wait": false })).unwrap()["job"].as_u64().unwrap();
    let r = e.wait(id).unwrap();
    assert_eq!(r["state"], "done");
    let entry = e.queue_entry(id).unwrap();
    let est = &entry["estimate"];
    assert_eq!(est["runs"], 2, "{entry}");
    assert_eq!(est["extrapolated"], false, "between the two sizes");
    let actual = r["worker"]["peak_working_set"].as_f64().unwrap();
    let predicted = est["bytes"].as_f64().unwrap();
    assert!(predicted / actual > 0.5 && predicted / actual < 2.0, "predicted {predicted}, measured {actual}");
    // a new engine reads what was learned
    let again = Engine::new();
    again.history().set_file(&file);
    assert_eq!(again.history().rows().len(), 3);
    let _ = std::fs::remove_dir_all(file.parent().unwrap());
}
