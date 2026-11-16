//! The queue (productivity plan Q1, ADR 87): jobs run in lanes, the heavy lane one at a time in order of priority, then round robin between
//! owners, then submission; a queued job starts by itself when the lane frees; requests of one project share the mesh; a convergence check
//! after a solve skips its coarse solve; the order of the queue never changes a number. The runner tests use closures that wait on a gate, so
//! no test depends on how long anything takes.

use femcore::progress::Progress;
use nyx_engine::{Engine, JobOutput, JobRunner, JobState, Lane, Mesher, QueueOpts, Submit};
use nyx_physics::Error;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

// ------------------------------------------------------------------------------------------------------------------ runner helpers

/// A gate that holds closures until it is opened.
#[derive(Clone, Default)]
struct Gate(Arc<(Mutex<bool>, Condvar)>);

impl Gate {
    fn closed() -> Gate {
        Gate::default()
    }
    fn open(&self) {
        *self.0 .0.lock().unwrap() = true;
        self.0 .1.notify_all();
    }
    fn pass(&self) {
        let mut g = self.0 .0.lock().unwrap();
        while !*g {
            g = self.0 .1.wait(g).unwrap();
        }
    }
}

type Work = Box<dyn FnOnce(Arc<Progress>) -> Result<JobOutput, Error> + Send>;

/// A job that waits for `gate` and then records its name.
fn gated(name: &str, gate: &Gate, order: &Arc<Mutex<Vec<String>>>) -> Work {
    let (name, gate, order) = (name.to_string(), gate.clone(), order.clone());
    Box::new(move |_| {
        order.lock().unwrap().push(name.clone());
        gate.pass();
        Ok(JobOutput::Exported(json!(name)))
    })
}

fn submit(r: &JobRunner, name: &str, owner: &str, priority: u8, w: Work) -> u64 {
    let mut s = Submit::new(name, Lane::Heavy);
    s.owner = owner.into();
    s.priority = priority;
    r.submit(s, w)
}

fn wait_for(what: &str, mut f: impl FnMut() -> bool) {
    let t = Instant::now();
    while !f() {
        assert!(t.elapsed() < Duration::from_secs(20), "timed out waiting for {what}");
        std::thread::sleep(Duration::from_millis(2));
    }
}

fn state(r: &JobRunner, id: u64) -> JobState {
    r.info(id).unwrap().state
}

/// Holds the heavy lane with one job so that what follows has to queue; returns the gate that frees it.
fn block_heavy(r: &JobRunner, order: &Arc<Mutex<Vec<String>>>) -> Gate {
    let g = Gate::closed();
    let id = submit(r, "blocker", "x", 1, gated("blocker", &g, order));
    wait_for("the blocker to run", || state(r, id) == JobState::Running);
    g
}

// ------------------------------------------------------------------------------------------------------------------ the runner

#[test]
fn a_queued_job_starts_by_itself_when_the_lane_frees() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    let g = Gate::closed();
    let a = submit(&r, "a", "gui", 2, gated("a", &g, &order));
    let b = submit(&r, "b", "gui", 2, gated("b", &Gate::default_open(), &order));
    wait_for("a to run", || state(&r, a) == JobState::Running);
    let info = r.info(b).unwrap();
    assert_eq!(info.state, JobState::Queued);
    assert_eq!(info.position, Some(1));
    assert!(info.seconds == 0.0);
    g.open();
    // nobody polls: `wait` only waits, the runner itself starts b when a ends
    assert!(r.wait(b).is_ok());
    assert!(r.wait(a).is_ok());
    assert_eq!(*order.lock().unwrap(), vec!["a", "b"]);
}

impl Gate {
    fn default_open() -> Gate {
        let g = Gate::closed();
        g.open();
        g
    }
}

#[test]
fn jobs_start_in_order_of_priority_then_submission() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    let blocker = block_heavy(&r, &order);
    let open = Gate::default_open();
    let ids: Vec<u64> = [("sweep", 3u8), ("agent", 2), ("person", 1), ("agent2", 2)]
        .iter()
        .map(|(n, p)| submit(&r, n, "y", *p, gated(n, &open, &order)))
        .collect();
    blocker.open();
    for id in &ids {
        r.wait(*id).unwrap();
    }
    assert_eq!(*order.lock().unwrap(), vec!["blocker", "person", "agent", "agent2", "sweep"]);
}

#[test]
fn owners_are_served_in_turn_within_a_priority() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    let blocker = block_heavy(&r, &order);
    let open = Gate::default_open();
    // an agent's sweep of three points is submitted first, then one request of the person at the same priority
    let ids: Vec<u64> = [("ai1", "ai"), ("ai2", "ai"), ("ai3", "ai"), ("gui1", "gui")]
        .iter()
        .map(|(n, o)| submit(&r, n, o, 2, gated(n, &open, &order)))
        .collect();
    blocker.open();
    for id in &ids {
        r.wait(*id).unwrap();
    }
    let got = order.lock().unwrap().clone();
    assert_eq!(got[1], "ai1");
    assert_eq!(got[2], "gui1", "the person is served after one item of the agent, not after all of them: {got:?}");
}

#[test]
fn a_job_that_waits_long_enough_rises_a_priority_level() {
    let r = JobRunner::new();
    r.set_age(Duration::from_millis(60));
    let order = Arc::new(Mutex::new(vec![]));
    let blocker = block_heavy(&r, &order);
    let open = Gate::default_open();
    let old = submit(&r, "old sweep", "y", 3, gated("old sweep", &open, &order));
    std::thread::sleep(Duration::from_millis(200)); // two levels: 3 becomes 1
    let fresh = submit(&r, "fresh", "y", 2, gated("fresh", &open, &order));
    blocker.open();
    r.wait(old).unwrap();
    r.wait(fresh).unwrap();
    assert_eq!(*order.lock().unwrap(), vec!["blocker", "old sweep", "fresh"], "ageing keeps a sweep from starving");
}

#[test]
fn a_queued_job_can_be_cancelled_held_and_superseded() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    let blocker = block_heavy(&r, &order);
    let open = Gate::default_open();
    let a = submit(&r, "a", "y", 2, gated("a", &open, &order));
    let b = submit(&r, "b", "y", 2, gated("b", &open, &order));
    let c = submit(&r, "c", "y", 2, gated("c", &open, &order));
    // cancel: never runs, finishes at once
    assert!(r.cancel(a));
    assert_eq!(state(&r, a), JobState::Cancelled);
    assert!(matches!(r.wait(a), Err(Error::Cancelled)));
    // hold: waits although the lane is free, until released
    assert!(r.hold(b));
    // supersede: the newer request replaces the queued one with the same key
    let mut s = Submit::new("d1", Lane::Heavy);
    s.owner = "y".into();
    s.key = Some("solve ems".into());
    // d1 was made without the flag: it can still be replaced by a newer request that has it
    let d1 = r.submit(s.clone(), gated("d1", &open, &order));
    s.label = "d2".into();
    s.supersede = true;
    let d2 = r.submit(s, gated("d2", &open, &order));
    assert_eq!(state(&r, d1), JobState::Cancelled);
    blocker.open();
    r.wait(c).unwrap();
    r.wait(d2).unwrap();
    assert_eq!(state(&r, b), JobState::Queued, "a held job does not start");
    assert!(r.info(b).unwrap().held);
    assert!(r.release(b));
    r.wait(b).unwrap();
    let got = order.lock().unwrap().clone();
    assert!(!got.contains(&"a".to_string()) && !got.contains(&"d1".to_string()), "cancelled jobs never ran: {got:?}");
}

#[test]
fn a_paused_queue_starts_nothing_and_resumes() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    r.pause();
    let a = submit(&r, "a", "y", 2, gated("a", &Gate::default_open(), &order));
    std::thread::sleep(Duration::from_millis(100));
    assert_eq!(state(&r, a), JobState::Queued);
    assert!(r.paused());
    r.resume();
    r.wait(a).unwrap();
}

#[test]
fn the_light_lane_runs_beside_a_busy_heavy_lane() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    let blocker = block_heavy(&r, &order);
    let mut s = Submit::new("io", Lane::Light);
    s.owner = "y".into();
    let io = r.submit(s, gated("io", &Gate::default_open(), &order));
    assert!(r.wait(io).is_ok(), "a file write does not wait for a solve");
    blocker.open();
}

#[test]
fn a_bounded_wait_returns_at_its_timeout_and_the_job_goes_on() {
    let r = JobRunner::new();
    let order = Arc::new(Mutex::new(vec![]));
    let g = block_heavy(&r, &order);
    let t = Instant::now();
    let id = r.running()[0];
    assert!(r.wait_timeout(id, Duration::from_millis(150)).is_none());
    assert!(t.elapsed() >= Duration::from_millis(150) && t.elapsed() < Duration::from_secs(5));
    assert_eq!(state(&r, id), JobState::Running, "the timeout did not stop the job");
    g.open();
    assert!(r.wait_timeout(id, Duration::from_secs(20)).unwrap().is_ok());
}

// ------------------------------------------------------------------------------------------------------------------ the engine

/// A mesher that counts its calls and, once closed, holds each mesh it has built until it is opened.
struct Probe {
    inner: nyx_engine::FemcoreMesher,
    calls: AtomicUsize,
    gate: Mutex<Gate>,
}

impl Probe {
    fn new() -> Arc<Probe> {
        let p = Probe {
            inner: nyx_engine::FemcoreMesher,
            calls: AtomicUsize::new(0),
            gate: Mutex::new(Gate::default_open()),
        };
        Arc::new(p)
    }
    fn close(&self) -> Gate {
        let g = Gate::closed();
        *self.gate.lock().unwrap() = g.clone();
        g
    }
    fn calls(&self) -> usize {
        self.calls.load(Ordering::SeqCst)
    }
}

impl Mesher for Probe {
    fn mesh(
        &self,
        geometry: &geom::model::FeatureList,
        settings: &nyx_engine::MeshSettings,
        sizing: Option<&nyx_physics::Sizing>,
    ) -> Result<nyx_engine::MeshedModel, Error> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        let m = self.inner.mesh(geometry, settings, sizing)?;
        let g = self.gate.lock().unwrap().clone();
        g.pass();
        Ok(m)
    }
}

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn open(name: &str) -> (Engine, Arc<Probe>) {
    let mut e = Engine::new();
    let p = Probe::new();
    e.set_mesher(p.clone());
    e.call("project_open", &json!({ "path": example(name) })).unwrap();
    (e, p)
}

fn neff_of(summary: &Value) -> String {
    summary["results"][0]["scalars"].to_string()
}

#[test]
fn two_requests_for_one_project_build_one_mesh() {
    let (mut e, probe) = open("15_soi_wire");
    let a = e.start_solve(None).unwrap();
    let b = e.start_solve(None).unwrap();
    let ra = e.wait(a).unwrap();
    let rb = e.wait(b).unwrap();
    assert_eq!(probe.calls(), 1, "the second request found the mesh the first built");
    assert_eq!(neff_of(&ra), neff_of(&rb));
    assert!(e.jobs().log(b).iter().any(|l| l.contains("earlier request")), "{:?}", e.jobs().log(b));
}

#[test]
fn a_convergence_check_after_a_solve_does_not_solve_the_coarse_mesh_again() {
    let (mut e, probe) = open("15_soi_wire");
    e.call("solve_run", &json!({})).unwrap();
    assert_eq!(probe.calls(), 1);
    let id = e.start_convergence("ems", 0.6, Some(1e-2)).unwrap();
    let r = e.wait(id).unwrap();
    assert_eq!(r["coarse_reused"], json!(true), "{r}");
    assert_eq!(probe.calls(), 2, "the mesh of the solve and the finer mesh; the present mesh was not built again");
    assert!(e.jobs().log(id).iter().any(|l| l.contains("reused from an earlier solve")));
    // the numbers are the ones of a check that solved both meshes itself
    let (mut cold, _) = open("15_soi_wire");
    let id = cold.start_convergence("ems", 0.6, Some(1e-2)).unwrap();
    let c = cold.wait(id).unwrap();
    assert_eq!(c["coarse_reused"], json!(false));
    assert_eq!(c["rows"], r["rows"], "reusing the coarse solve gives the same numbers");
}

#[test]
fn an_analysis_finishes_while_a_heavy_job_runs() {
    let (mut e, probe) = open("15_soi_wire");
    e.call("solve_run", &json!({})).unwrap();
    let gate = probe.close();
    e.call("set_param", &json!({ "name": "w", "value": 6.0e-7 })).unwrap();
    let id = e.start_solve(None).unwrap();
    wait_for("the solve to hold the heavy lane", || probe.calls() == 2);
    assert_eq!(e.jobs().info(id).unwrap().state, JobState::Running);
    // the engine answers reads and analyses of the result it holds while the heavy lane is busy
    let t = Instant::now();
    let a = e.call("analyse", &json!({ "analysis": "polarization", "physics": "ems" })).unwrap();
    assert_eq!(a["stale"], json!(true), "the edit made the held result stale");
    assert!(t.elapsed() < Duration::from_secs(5));
    assert_eq!(
        e.jobs().info(id).unwrap().state,
        JobState::Running,
        "the analysis did not wait for the solve, and the solve was still running"
    );
    gate.open();
    e.wait(id).unwrap();
}

#[test]
fn cancelling_one_of_two_requests_that_share_a_mesh_keeps_the_mesh() {
    let (mut e, probe) = open("15_soi_wire");
    let gate = probe.close();
    let a = e.start_solve(None).unwrap();
    let b = e.start_solve(None).unwrap();
    wait_for("the mesh to be built", || probe.calls() == 1);
    assert!(e.jobs().cancel(a));
    gate.open();
    assert!(matches!(e.wait(a), Err(Error::Cancelled)));
    let rb = e.wait(b).unwrap();
    assert_eq!(rb["state"], json!("done"));
    assert_eq!(probe.calls(), 1, "the mesh the cancelled request built served the other one");
}

#[test]
fn a_request_cancelled_in_the_queue_never_runs_and_leaves_no_trace() {
    let (mut e, probe) = open("15_soi_wire");
    let gate = probe.close();
    let a = e.start_solve(None).unwrap();
    let b = e.start_solve_with(None, &QueueOpts::default()).unwrap();
    wait_for("a to hold the lane", || probe.calls() == 1);
    assert_eq!(e.jobs().info(b).unwrap().state, JobState::Queued);
    assert!(e.jobs().cancel(b));
    gate.open();
    e.wait(a).unwrap();
    assert!(matches!(e.wait(b), Err(Error::Cancelled)));
    assert!(e.jobs().log(b).is_empty(), "it never ran");
}

#[test]
fn the_result_does_not_depend_on_the_order_the_queue_runs_the_requests() {
    // three widths of the wire, as three frozen requests; the second engine runs them in the opposite order
    let run = |order: [usize; 3], priorities: [u8; 3]| -> Vec<String> {
        let (mut e, _) = open("15_soi_wire");
        let widths = [4.5e-7, 5.0e-7, 5.5e-7];
        let mut ids = [0u64; 3];
        for (k, &i) in order.iter().enumerate() {
            e.call("set_param", &json!({ "name": "w", "value": widths[i] })).unwrap();
            ids[i] =
                e.start_solve_with(None, &QueueOpts { priority: Some(priorities[k]), ..Default::default() }).unwrap();
        }
        // the answers by width, whatever order they ran in
        (0..3).map(|i| neff_of(&e.wait(ids[i]).unwrap())).collect()
    };
    let forward = run([0, 1, 2], [2, 2, 2]);
    let reversed = run([2, 1, 0], [2, 2, 2]);
    let reprioritised = run([0, 1, 2], [3, 2, 1]); // the first one starts at once, then the last submitted runs before the middle one
    assert_eq!(forward, reversed, "submission order changed a number");
    assert_eq!(forward, reprioritised, "priority changed a number");
    assert_ne!(forward[0], forward[2], "the widths differ, so the test can see a mix-up");
}

#[test]
fn the_queue_commands_submit_list_hold_and_wait_with_a_timeout() {
    let (mut e, probe) = open("15_soi_wire");
    let gate = probe.close();
    let first = e.call("queue_submit", &json!({})).unwrap();
    let id1 = first["id"].as_u64().unwrap();
    wait_for("the first job to hold the lane", || probe.calls() == 1);
    let second =
        e.call("queue_submit", &json!({ "kind": "solve", "physics": "ems", "priority": 3, "hold": true })).unwrap();
    let id2 = second["id"].as_u64().unwrap();
    assert_eq!(second["state"], json!("queued"));
    assert_eq!(second["held"], json!(true));
    assert_eq!(second["kind"], json!("solve"));
    assert_eq!(second["targets"], json!(["ems"]));
    assert_eq!(second["project_changed"], json!(false));
    let list = e.call("queue_list", &json!({ "state": "queued" })).unwrap();
    assert_eq!(list.as_array().unwrap().len(), 1);
    // a bounded wait on a job that is still running returns at the timeout
    let t = Instant::now();
    let w = e.call("queue_wait", &json!({ "id": id1, "timeout_s": 0.2 })).unwrap();
    assert_eq!(w["waiting"], json!(true), "{w}");
    assert_eq!(w["state"], json!("running"));
    assert!(t.elapsed() < Duration::from_secs(5));
    // the timeout is required
    assert!(e.call("queue_wait", &json!({ "id": id1 })).is_err());
    // the project moves on: the queued request says so and keeps its own frozen copy
    e.call("set_param", &json!({ "name": "w", "value": 6.0e-7 })).unwrap();
    assert_eq!(e.call("queue_get", &json!({ "id": id2 })).unwrap()["project_changed"], json!(true));
    gate.open();
    let done = e.call("queue_wait", &json!({ "id": id1, "timeout_s": 60 })).unwrap();
    assert_eq!(done["state"], json!("done"));
    assert!(e.call("queue_release", &json!({ "id": id2 })).unwrap()["released"].as_bool().unwrap());
    let done2 = e.call("queue_wait", &json!({ "id": id2, "timeout_s": 60 })).unwrap();
    assert_eq!(done2["state"], json!("done"));
    // the frozen copy was solved: the width it was made with, not the edited one
    assert_eq!(done2["results"][0]["current"], json!(false));
    assert!(e.call("queue_set_priority", &json!({ "id": id1, "priority": 9 })).is_err());
}

#[test]
fn a_failed_physics_keeps_the_ones_solved_before_it() {
    let mut e = Engine::new();
    e.call("project_open", &json!({ "path": example("04_thermal_bar") })).unwrap();
    // the solid problem names a boundary the mesh does not have: it fails after the heat problem was solved
    let mut sms = serde_json::to_value(e.project().spec("sms").expect("example 04 has a solid problem")).unwrap();
    sms["bcs"] = json!([{ "on": "no_such_boundary", "kind": "fixed", "params": {} }]);
    if e.call("set_problem", &json!({ "physics": "sms", "problem": sms })).is_err() {
        panic!("set_problem refused the broken problem: this test needs another way to fail the second physics");
    }
    let id = e.start_solve(Some(vec!["dfs".into(), "sms".into()])).unwrap();
    let r = e.wait(id);
    let why = r.expect_err("the solid problem must fail").to_string();
    assert!(why.contains("no_such_boundary"), "it failed for another reason than the one the test made: {why}");
    assert!(e.result("dfs").is_some(), "the heat result that was solved before the failure is kept");
    assert!(e.result("sms").is_none());
}

#[test]
fn a_wait_on_a_job_that_cannot_start_refuses_instead_of_hanging() {
    let (mut e, _) = open("15_soi_wire");
    // a paused queue starts nothing: a wait with no timeout would never return
    e.call("queue_pause", &json!({})).unwrap();
    let a = e.start_solve(None).unwrap();
    let err = e.wait(a).expect_err("a wait on a paused queue must refuse").to_string();
    assert!(err.contains("queue_resume"), "{err}");
    let w = e.call("queue_wait", &json!({ "id": a, "timeout_s": 0.1 })).unwrap();
    assert_eq!(w["state"], json!("queued"), "a bounded wait is fine, and says the job is still queued: {w}");
    e.call("queue_resume", &json!({})).unwrap();
    // a held job is the same
    let b = e.call("queue_submit", &json!({ "hold": true })).unwrap()["id"].as_u64().unwrap();
    e.wait(a).unwrap();
    let err = e.call("job_wait", &json!({ "id": b })).expect_err("a wait on a held job must refuse").to_string();
    assert!(err.contains("queue_release"), "{err}");
    e.call("queue_release", &json!({ "id": b })).unwrap();
    assert_eq!(e.call("job_wait", &json!({ "id": b })).unwrap()["state"], json!("done"));
}

#[test]
fn the_convergence_check_command_goes_through_the_queue_and_reuses_the_coarse_solve() {
    let (mut e, probe) = open("15_soi_wire");
    e.call("solve_run", &json!({})).unwrap();
    let _ = probe;
    // the plain command (wait true, the default) is a job of the heavy lane like a solve, so it can reuse the solve before it
    let jobs_before = e.jobs().list().len();
    let r = e.call("convergence_check", &json!({ "physics": "ems", "factor": 0.6, "tolerance": 1e-2 })).unwrap();
    assert_eq!(r["coarse_reused"], json!(true), "{r}");
    assert_eq!(e.jobs().list().len(), jobs_before + 1, "the check was a job of the queue");
    // and with the queue paused it cannot start, which a wait with no timeout refuses
    e.call("queue_pause", &json!({})).unwrap();
    let err = e.call("convergence_check", &json!({ "physics": "ems" })).expect_err("a paused queue").to_string();
    assert!(err.contains("queue_resume"), "{err}");
}

#[test]
fn a_newer_request_replaces_a_queued_one_even_when_the_first_was_made_without_the_flag() {
    let (mut e, _) = open("15_soi_wire");
    e.call("queue_pause", &json!({})).unwrap();
    let a = e.start_solve(None).unwrap(); // no supersede flag
    let b = e.start_solve_with(None, &QueueOpts { supersede: true, ..Default::default() }).unwrap();
    assert_eq!(e.jobs().info(a).unwrap().state, JobState::Cancelled);
    assert_eq!(e.jobs().info(b).unwrap().state, JobState::Queued);
    e.call("queue_resume", &json!({})).unwrap();
    assert_eq!(e.wait(b).unwrap()["state"], json!("done"));
    assert!(matches!(e.wait(a), Err(Error::Cancelled)));
}
