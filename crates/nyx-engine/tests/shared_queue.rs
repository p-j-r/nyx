//! The queue in the workspace folder (productivity plan Q3, ADR 95). Real processes where the gate is about processes: the queue worker is a
//! program of its own (`nyx-worker --worker FOLDER`), a submitter that exits is a child process (`nyx-worker queue submit`), and a worker is
//! killed with the system's own kill. Two `Engine`s in this test process stand for two front-ends (the window and the MCP server): they share
//! nothing but the folder, which is all two programs would share.

use nyx_engine::shared_queue::SharedQueue;
use nyx_engine::Engine;
use nyx_physics::{Error, SolveResult};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

fn worker() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"))
}

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples").join(format!("{name}.nyx.json"))
}

/// A folder of its own for each test, and a queue worker that ends soon after the work is done.
fn workspace(tag: &str) -> PathBuf {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| std::env::set_var("NYX_QUEUE_IDLE_SECS", "3"));
    static N: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let d = std::env::temp_dir().join(format!(
        "nyx-sq-{tag}-{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn engine(ws: &Path, owner: &str, name: &str) -> Engine {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example(name).to_string_lossy() }))
        .unwrap_or_else(|err| panic!("{name}: {err}"));
    e.set_worker(worker(), false);
    e.set_shared_queue(ws, false);
    e.set_owner(owner);
    e
}

/// 14_box_cavity at a mesh that takes seconds to factorise, saved as a file the command line can submit.
fn slow_project(ws: &Path) -> PathBuf {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("14_box_cavity").to_string_lossy() })).unwrap();
    e.call("set_mesh", &json!({ "h": 0.0028 })).unwrap();
    let p = ws.join("slow.nyx.json");
    e.call("project_save", &json!({ "path": p.to_string_lossy() })).unwrap();
    p
}

fn wait_for(what: &str, secs: u64, mut f: impl FnMut() -> bool) {
    let t = Instant::now();
    while !f() {
        assert!(t.elapsed() < Duration::from_secs(secs), "timed out after {secs} s waiting for: {what}");
        std::thread::sleep(Duration::from_millis(50));
    }
}

fn kill(pid: u32) {
    #[cfg(windows)]
    let _ = Command::new("taskkill").args(["/F", "/PID", &pid.to_string()]).output();
    #[cfg(not(windows))]
    let _ = Command::new("kill").args(["-9", &pid.to_string()]).output();
}

fn state_of(q: &SharedQueue, id: &str) -> String {
    q.item(id).map(|i| i.state_name().to_string()).unwrap_or_else(|| "gone".into())
}

/// Submits `project` from a program that exits at once (the window closed, the script ended) and returns the item's id.
fn submit_from_a_process(ws: &Path, project: &Path, extra: &[&str]) -> String {
    let out = Command::new(worker())
        .args(["queue", "submit", &ws.to_string_lossy(), &project.to_string_lossy()])
        .args(extra)
        .output()
        .expect("the program");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8_lossy(&out.stdout).trim().strip_prefix("queued ").expect("queued <id>").to_string()
}

fn same_numbers(a: &SolveResult, b: &SolveResult, what: &str) {
    assert_eq!(a.physics, b.physics, "{what}");
    assert_eq!(a.dofs, b.dofs, "{what}: unknowns");
    assert_eq!(a.scalars, b.scalars, "{what}: scalars");
    assert_eq!(a.solutions, b.solutions, "{what}: solutions");
}

#[test]
fn two_front_ends_see_one_ordered_queue_and_the_numbers_are_those_of_an_in_process_solve() {
    let ws = workspace("two");
    let q = SharedQueue::new(&ws);
    // the queue worker is started by the first job that needs it; paused, it ranks the items and starts none
    let mut window = engine(&ws, "gui", "01_capacitor");
    let mut agent = engine(&ws, "ai", "02_cantilever");
    q.pause();
    let a = agent.call("queue_submit", &json!({ "shared": true, "priority": 3 })).unwrap()["id"].as_u64().unwrap();
    let w = window.call("queue_submit", &json!({ "shared": true })).unwrap()["id"].as_u64().unwrap();
    // both appear in the folder, ranked: the person's request (priority 1) before the agent's (3), whoever submitted first
    wait_for("both items ranked", 60, || q.items().len() == 2 && q.items().iter().all(|i| i.position.is_some()));
    let items = q.items();
    let first = items.iter().find(|i| i.position == Some(1)).unwrap();
    assert_eq!(first.request.owner, "gui", "{items:?}");
    assert!(q.runner_alive(), "the queue worker runs");
    // each front-end lists the other's item, marked as foreign, and its own
    let seen_by_window = window.call("queue_list", &json!({})).unwrap();
    let seen_by_window = seen_by_window.as_array().unwrap();
    assert!(seen_by_window.iter().any(|e| e["foreign"] == json!(true) && e["owner"] == "ai"), "{seen_by_window:?}");
    assert!(
        seen_by_window.iter().any(|e| e["id"].as_u64() == Some(w) && e["shared_id"].is_string()),
        "{seen_by_window:?}"
    );
    // while the queue is paused nothing started
    assert!(items.iter().all(|i| i.state.is_none()), "{items:?}");
    q.resume();
    let rw = window.call("queue_wait", &json!({ "id": w, "timeout_s": 120.0 })).unwrap();
    let ra = agent.call("queue_wait", &json!({ "id": a, "timeout_s": 120.0 })).unwrap();
    assert!(rw["worker"].is_object() && ra["worker"].is_object(), "{rw} {ra}");
    // the same solves in this process, to the bit
    for (name, e) in [("01_capacitor", &window), ("02_cantilever", &agent)] {
        let mut here = Engine::new();
        here.call("project_open", &json!({ "path": example(name).to_string_lossy() })).unwrap();
        here.call("solve_run", &json!({ "isolate": false })).unwrap();
        for spec in here.project().physics.clone() {
            let (x, _) = here.result(&spec.physics).unwrap();
            let (y, stale) =
                e.result(&spec.physics).unwrap_or_else(|| panic!("{name}: no result for {}", spec.physics));
            assert!(!stale);
            same_numbers(x, y, name);
        }
    }
    // the front-ends collected their items: nothing is left but the history, which the queue worker wrote (one row per run, not two)
    assert!(q.items().is_empty(), "{:?}", q.items());
    let rows = std::fs::read_to_string(ws.join("solve_history.jsonl")).unwrap_or_default();
    assert_eq!(rows.lines().count(), 2, "one history row for each run: {rows}");
}

#[test]
fn closing_the_front_end_leaves_the_queued_item_running_and_its_results_are_kept() {
    let ws = workspace("close");
    let q = SharedQueue::new(&ws);
    // a program submits and exits: nothing of it is alive
    let id = submit_from_a_process(&ws, &example("01_capacitor"), &[]);
    wait_for("the item to be done", 90, || state_of(&q, &id) == "done");
    let item = q.item(&id).unwrap();
    assert!(item.state.as_ref().unwrap().has_results);
    assert!(q.item_dir(&id).join("results.h5").exists(), "the results file is there for whoever comes back");
    // a new front-end finds it, lists it as another program's, and opens its results
    let mut e = Engine::new();
    e.set_shared_queue(&ws, false);
    let list = e.call("queue_list", &json!({ "state": "done" })).unwrap();
    assert!(
        list.as_array().unwrap().iter().any(|v| v["shared_id"] == json!(id) && v["foreign"] == json!(true)),
        "{list}"
    );
    let path = item.to_json(&q)["results"].as_str().unwrap().to_string();
    e.set_worker(worker(), false);
    e.call("open_results", &json!({ "path": path })).unwrap();
    assert!(e.result("dfs").is_some(), "the results of the item are open");
    // the command line lists it too
    let out = Command::new(worker()).args(["queue", "list", &ws.to_string_lossy()]).output().unwrap();
    assert!(String::from_utf8_lossy(&out.stdout).contains(&id), "{}", String::from_utf8_lossy(&out.stdout));
}

#[test]
fn killing_a_solver_fails_only_its_item() {
    let ws = workspace("killsolver");
    let q = SharedQueue::new(&ws);
    let slow = slow_project(&ws);
    let a = submit_from_a_process(&ws, &slow, &[]);
    let b = submit_from_a_process(&ws, &example("01_capacitor"), &[]);
    wait_for("the slow item to run in a solver child", 90, || {
        q.item(&a).and_then(|i| i.state).is_some_and(|s| s.state == "running" && s.child_pid.is_some())
    });
    kill(q.item(&a).unwrap().state.unwrap().child_pid.unwrap());
    wait_for("the killed item to fail and the next to be done", 120, || {
        state_of(&q, &a) == "failed" && state_of(&q, &b) == "done"
    });
    let err = q.item(&a).unwrap().state.unwrap().error.unwrap().to_string();
    assert!(err.contains("stopped without an answer") && err.contains("Nothing else was affected"), "{err}");
    assert!(q.runner_alive() || state_of(&q, &b) == "done", "the queue worker went on");
}

#[test]
fn a_queue_worker_that_was_killed_fails_the_item_it_ran_and_the_next_one_runs_the_rest() {
    let ws = workspace("killrunner");
    let q = SharedQueue::new(&ws);
    let slow = slow_project(&ws);
    let a = submit_from_a_process(&ws, &slow, &[]);
    wait_for("the slow item to run", 90, || {
        q.item(&a).and_then(|i| i.state).is_some_and(|s| s.state == "running" && s.child_pid.is_some())
    });
    let lease: Value = serde_json::from_str(&std::fs::read_to_string(q.root().join("runner.lock")).unwrap()).unwrap();
    kill(lease["pid"].as_u64().unwrap() as u32);
    let b = submit_from_a_process(&ws, &example("01_capacitor"), &[]);
    // the dead lease goes stale; whoever looks next starts a queue worker (as a waiting job does every two seconds)
    wait_for("the new queue worker to fail the old item and run the new one", 120, || {
        let _ = q.ensure_runner(&worker());
        state_of(&q, &a) == "failed" && state_of(&q, &b) == "done"
    });
    let err = q.item(&a).unwrap().state.unwrap().error.unwrap().to_string();
    assert!(err.contains("queue worker"), "{err}");
}

#[test]
fn a_cancel_from_another_program_stops_a_running_solve_at_once_and_a_queued_one_never_runs() {
    let ws = workspace("cancel");
    let q = SharedQueue::new(&ws);
    let slow = slow_project(&ws);
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": slow.to_string_lossy() })).unwrap();
    e.set_worker(worker(), false);
    e.set_shared_queue(&ws, false);
    let running = e.call("queue_submit", &json!({ "shared": true })).unwrap()["id"].as_u64().unwrap();
    wait_for("it to run in a solver child", 90, || {
        q.items().first().and_then(|i| i.state.clone()).is_some_and(|s| s.state == "running" && s.child_pid.is_some())
    });
    // a second item behind it, from a process that exits
    let queued = submit_from_a_process(&ws, &example("01_capacitor"), &[]);
    // another program cancels both by their ids
    let mut other = Engine::new();
    other.set_shared_queue(&ws, false);
    let ids: Vec<String> = q.items().iter().map(|i| i.request.id.clone()).collect();
    let t = Instant::now();
    for id in ids.iter().rev() {
        assert_eq!(other.call("queue_cancel", &json!({ "shared_id": id })).unwrap()["cancelled"], json!(true));
    }
    let err = e.call("queue_wait", &json!({ "id": running, "timeout_s": 60.0 })).unwrap_err();
    assert!(matches!(err, Error::Cancelled), "{err}");
    assert!(t.elapsed() < Duration::from_secs(5), "a kill does not wait for a factorisation: {:?}", t.elapsed());
    wait_for("the queued item to be cancelled", 30, || state_of(&q, &queued) == "cancelled");
}

#[test]
fn hold_release_and_priority_work_on_items_by_their_shared_id() {
    let ws = workspace("steer");
    let q = SharedQueue::new(&ws);
    let held = submit_from_a_process(&ws, &example("01_capacitor"), &["--hold"]);
    let other = submit_from_a_process(&ws, &example("01_capacitor"), &["--priority", "3"]);
    wait_for("the unheld item to be done", 90, || state_of(&q, &other) == "done");
    assert_eq!(state_of(&q, &held), "held", "a held item does not start");
    let mut e = Engine::new();
    e.set_shared_queue(&ws, false);
    assert_eq!(
        e.call("queue_set_priority", &json!({ "shared_id": held, "priority": 1 })).unwrap()["changed"],
        json!(true)
    );
    assert_eq!(e.call("queue_release", &json!({ "shared_id": held })).unwrap()["released"], json!(true));
    wait_for("the released item to be done", 90, || state_of(&q, &held) == "done");
    assert_eq!(q.item(&held).unwrap().priority, 1);
    let got = e.call("queue_get", &json!({ "shared_id": held })).unwrap();
    assert_eq!(got["state"], json!("done"), "{got}");
}

#[test]
fn a_job_waiting_in_the_workspace_queue_does_not_block_this_engines_own_heavy_lane() {
    let ws = workspace("lanes");
    let q = SharedQueue::new(&ws);
    let mut e = engine(&ws, "gui", "01_capacitor");
    q.pause();
    let id = e.call("queue_submit", &json!({ "shared": true })).unwrap()["id"].as_u64().unwrap();
    wait_for("the item to be written", 60, || !q.items().is_empty());
    // the shared job is queued remotely; an in-process solve of the same engine runs at once
    let t = Instant::now();
    e.call("solve_run", &json!({ "isolate": false, "shared": false })).unwrap();
    assert!(t.elapsed() < Duration::from_secs(30));
    wait_for("the job to say it is queued in the workspace queue", 30, || {
        let entry = e.call("queue_get", &json!({ "id": id })).unwrap();
        entry["stage"].as_str().unwrap_or("").contains("queued")
            && entry["state"] == json!("running")
            && entry["lane"] == json!("shared")
    });
    q.resume();
    e.call("queue_wait", &json!({ "id": id, "timeout_s": 120.0 })).unwrap();
}

#[test]
fn a_wait_with_no_timeout_for_a_job_in_a_paused_or_held_queue_is_refused_and_does_not_hang() {
    let ws = workspace("blocked");
    let q = SharedQueue::new(&ws);
    let mut e = engine(&ws, "ai", "01_capacitor");
    q.pause();
    let id = e.call("queue_submit", &json!({ "shared": true })).unwrap()["id"].as_u64().unwrap();
    wait_for("the item to be written", 60, || !q.items().is_empty());
    assert!(e.job_blocked(id).is_some_and(|w| w.contains("paused")), "{:?}", e.job_blocked(id));
    let err = e.call("job_wait", &json!({ "id": id })).unwrap_err().to_string();
    assert!(err.contains("cannot start") && err.contains("paused"), "{err}");
    q.resume();
    let sid = q.items()[0].request.id.clone();
    assert!(q.hold(&sid));
    wait_for("the item to be held and unstarted", 10, || e.job_blocked(id).is_some_and(|w| w.contains("held")));
    assert!(q.release(&sid));
    e.call("queue_wait", &json!({ "id": id, "timeout_s": 120.0 })).unwrap();
}
