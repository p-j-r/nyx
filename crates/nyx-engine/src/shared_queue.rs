// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The queue in the workspace folder (productivity plan Q3, ADR 95): the one queue that the window, the MCP server and the command line share, and
//! the queue worker that runs it. A front-end that has a workspace writes a solve as an *item* in `<workspace>/queue/items/<id>/` and watches it;
//! a program started with `--worker <workspace>` (the *queue worker*, a process of its own that outlives the window) runs the items one at a time,
//! each in a solver child exactly as `worker.rs` runs one (ADR 86), and writes what happens into the item. Closing a front-end leaves its items
//! running; killing a solver child fails its item and nothing else; a queue worker that died is replaced by the next front-end that needs one, and
//! the item it was running is found `failed`.
//!
//! **Who writes what** (so no file has two writers and nothing needs a lock but the queue worker's lease):
//! * a front-end writes an item once, as a whole folder renamed into place (`request.json`, `project.nyx.json`, `job.json`), and afterwards only
//!   *markers*: an empty file `cancel`, an empty file `hold` (removed to release), a file `priority` holding 1, 2 or 3;
//! * the queue worker writes `state.json` (queued items have none), `order.json` (the ids in the order they will start), the lease and the
//!   results (`results.h5`, `stderr.txt`) of an item; it also appends to `solve_history.jsonl`, so the estimates are learned whoever was watching.
//!
//! **The lease** is `queue/runner.lock`: created with `create_new` by the queue worker that takes the queue, kept alive by touching its modification
//! time every second. A lease not touched for six seconds is dead (a process id alone is not trusted: Windows reuses them). Two programs that find
//! a dead lease both try to take it; only one `create_new` succeeds, and each re-reads the file to see whose token it holds.
//!
//! Order, ageing and owners are the job runner's rules (ADR 87): priority first, lowered one level for each minute waited, then owners in turn,
//! then submission order. The order never changes a number (ADR 11): an item is a frozen project and nothing carries from one to the next.

use crate::estimate::{Estimate, History, MemoryProbe};
use nyx_physics::Error;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
#[cfg(feature = "hdf5")]
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

/// The folder the window and the MCP server share by default: `NYX_WORKSPACE`, else the person's `Documents/nyx` (the window's workspace).
pub fn default_workspace() -> PathBuf {
    std::env::var_os("NYX_WORKSPACE").map(PathBuf::from).unwrap_or_else(nyx_store::workspace::default_root)
}

/// A lease not touched for this long belongs to a queue worker that is gone.
pub const LEASE_DEAD: Duration = Duration::from_secs(6);
/// A finished item nobody collected is removed after this long.
pub const KEEP_FINISHED: Duration = Duration::from_secs(24 * 3600);

fn now_ms() -> u64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_millis() as u64).unwrap_or(0)
}

fn default_true() -> bool {
    true
}

/// What was asked, written once by the front-end.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Request {
    pub id: String,
    /// `gui`, `ai` or `cli`
    pub owner: String,
    /// 1 a person's request now, 2 an agent's, 3 sweeps and batches
    pub priority: u8,
    pub label: String,
    pub targets: Vec<String>,
    /// the hash of the project the request froze
    pub stamp: String,
    pub submitted_ms: u64,
    /// start even when the memory estimate says the machine cannot hold it
    pub force: bool,
    /// the estimate the front-end made once its mesh was known (the queue worker admits by it)
    pub estimate: Option<Estimate>,
    /// the process that submitted it, for the list ("from window, pid 1234"); not trusted for anything
    pub from_pid: u32,
    #[serde(default = "default_true")]
    pub save_fields: bool,
    #[serde(default)]
    pub post_analyses: Vec<Value>,
    #[serde(default)]
    pub memory_wait_timeout_s: Option<f64>,
}

/// What a front-end puts in an item besides the frozen project.
#[derive(Clone, Debug)]
pub(crate) struct NewItem {
    pub owner: String,
    pub priority: u8,
    pub label: String,
    pub stamp: String,
    pub hold: bool,
    pub force: bool,
    pub estimate: Option<Estimate>,
    pub save_fields: bool,
    pub post_analyses: Vec<Value>,
    pub memory_wait_timeout_s: Option<f64>,
}

/// What the queue worker says of an item. An item with no `state.json` is queued.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct State {
    /// `running`, `done`, `failed` or `cancelled`
    pub state: String,
    pub stage: String,
    pub stage_seconds: f64,
    pub fraction: Option<f64>,
    pub log: Vec<String>,
    /// `{code, message}` of a failure
    pub error: Option<Value>,
    pub has_results: bool,
    pub peak_working_set: Option<u64>,
    pub peak_commit: Option<u64>,
    pub started_ms: Option<u64>,
    pub ended_ms: Option<u64>,
    /// the solver child, so a person can see which process it is
    pub child_pid: Option<u32>,
    pub runner_pid: u32,
}

impl State {
    pub fn is_finished(&self) -> bool {
        matches!(self.state.as_str(), "done" | "failed" | "cancelled")
    }
}

/// An item as a front-end sees it: the request, the worker's state and the markers.
#[derive(Clone, Debug)]
pub struct Item {
    pub request: Request,
    pub state: Option<State>,
    pub held: bool,
    pub cancel_requested: bool,
    /// the priority marker, when one was written
    pub priority: u8,
    /// the place in the order of starting (1 is next), for a queued item the queue worker has ranked
    pub position: Option<usize>,
}

impl Item {
    /// `queued`, `held`, `running`, `done`, `failed` or `cancelled`
    pub fn state_name(&self) -> &str {
        match &self.state {
            Some(s) => &s.state,
            None if self.cancel_requested => "cancelled",
            None if self.held => "held",
            None => "queued",
        }
    }

    pub fn is_open(&self) -> bool {
        !matches!(self.state_name(), "done" | "failed" | "cancelled")
    }

    /// The entry the queue commands show.
    pub fn to_json(&self, workspace_queue: &SharedQueue) -> Value {
        let r = &self.request;
        let mut v = json!({
            "shared_id": r.id, "shared": true, "label": r.label, "owner": r.owner, "kind": "solve", "targets": r.targets,
            "priority": self.priority, "state": self.state_name(), "held": self.held,
            "cancel_requested": self.cancel_requested && self.is_open(),
            "submitted_ms": r.submitted_ms, "from_pid": r.from_pid, "stamp": r.stamp,
        });
        if let Some(p) = self.position {
            v["position"] = json!(p);
        }
        if let Some(e) = &r.estimate {
            v["estimate"] = json!(e);
        }
        if let Some(s) = &self.state {
            v["stage"] = json!(s.stage);
            v["stage_seconds"] = json!(s.stage_seconds);
            v["fraction"] = json!(s.fraction);
            v["error"] = json!(s.error);
            v["has_results"] = json!(s.has_results);
            v["child_pid"] = json!(s.child_pid);
            v["peak_working_set"] = json!(s.peak_working_set);
        }
        v["results"] = json!(workspace_queue.item_dir(&r.id).join("results.h5").to_string_lossy());
        v
    }
}

// ------------------------------------------------------------------------------------------------------------------ files

/// Writes `bytes` so that a reader sees the old file or the new one, never half of one: a temporary file renamed over the name. On Windows the
/// rename fails while another process has the file open, so it is tried again for a moment.
fn write_atomic(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    static N: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let tmp = path.with_extension(format!(
        "tmp{}-{}",
        std::process::id(),
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    std::fs::write(&tmp, bytes)?;
    let mut last = None;
    for _ in 0..50 {
        match std::fs::rename(&tmp, path) {
            Ok(()) => return Ok(()),
            Err(e) => {
                last = Some(e);
                std::thread::sleep(Duration::from_millis(10));
            }
        }
    }
    let _ = std::fs::remove_file(&tmp);
    Err(last.unwrap_or_else(|| std::io::Error::other("rename failed")))
}

/// Reads a JSON file another process may be replacing; a read that fails or does not parse is tried again for a moment.
fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Option<T> {
    for _ in 0..20 {
        match std::fs::read_to_string(path) {
            Ok(text) => {
                if let Ok(v) = serde_json::from_str::<T>(&text) {
                    return Some(v);
                }
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return None,
            Err(_) => {}
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    None
}

fn age_of(path: &Path) -> Option<Duration> {
    std::fs::metadata(path).ok()?.modified().ok().map(|t| SystemTime::now().duration_since(t).unwrap_or_default())
}

/// The priority a queued item competes with: its own, lowered one level for each `age` waited, never below 1 (as the job runner's).
fn effective(priority: u8, waited: Duration, age: Duration) -> u8 {
    let levels = if age.is_zero() { 0 } else { (waited.as_secs_f64() / age.as_secs_f64()) as u64 };
    (priority as u64).saturating_sub(levels).max(1) as u8
}

// ------------------------------------------------------------------------------------------------------------------ the front-end side

/// The queue folder of one workspace.
#[derive(Clone, Debug)]
pub struct SharedQueue {
    root: PathBuf,
}

impl SharedQueue {
    /// The queue of `workspace` (the folder is made when the first item is written).
    pub fn new(workspace: &Path) -> SharedQueue {
        SharedQueue { root: workspace.join("queue") }
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    fn items_dir(&self) -> PathBuf {
        self.root.join("items")
    }

    pub fn item_dir(&self, id: &str) -> PathBuf {
        self.items_dir().join(id)
    }

    fn lease_path(&self) -> PathBuf {
        self.root.join("runner.lock")
    }

    /// A new id: sortable by time of submission, unique over processes.
    fn new_id() -> String {
        static N: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        format!(
            "{:013}-{:06x}-{:03x}",
            now_ms(),
            std::process::id() & 0xff_ffff,
            N.fetch_add(1, std::sync::atomic::Ordering::SeqCst) & 0xfff
        )
    }

    /// Writes a solve as an item and returns its id. The folder appears under its name only when it is complete.
    #[cfg(feature = "hdf5")]
    pub(crate) fn submit(&self, n: &NewItem, spec: &crate::worker::RemoteSpec) -> Result<String, Error> {
        let io = |what: &str, e: std::io::Error| {
            Error::Solver(format!("the queue folder {}: {what}: {e}", self.root.display()))
        };
        std::fs::create_dir_all(self.items_dir()).map_err(|e| io("making it", e))?;
        let id = Self::new_id();
        let tmp = self.items_dir().join(format!(".new-{id}"));
        std::fs::create_dir_all(&tmp).map_err(|e| io("making an item", e))?;
        let request = Request {
            id: id.clone(),
            owner: n.owner.clone(),
            priority: n.priority,
            label: n.label.clone(),
            targets: spec.targets.clone(),
            stamp: n.stamp.clone(),
            submitted_ms: now_ms(),
            force: n.force,
            estimate: n.estimate.clone(),
            from_pid: std::process::id(),
            save_fields: n.save_fields,
            post_analyses: n.post_analyses.clone(),
            memory_wait_timeout_s: n.memory_wait_timeout_s,
        };
        let r = (|| -> Result<(), Error> {
            crate::worker::write_job_files(&tmp, spec)?;
            std::fs::write(tmp.join("request.json"), json!(request).to_string())
                .map_err(|e| io("writing the request", e))?;
            if n.hold {
                std::fs::write(tmp.join("hold"), b"").map_err(|e| io("writing the hold marker", e))?;
            }
            let mut last = None;
            for _ in 0..50 {
                match std::fs::rename(&tmp, self.item_dir(&id)) {
                    Ok(()) => return Ok(()),
                    Err(e) => {
                        last = Some(e);
                        std::thread::sleep(Duration::from_millis(10));
                    }
                }
            }
            Err(io("putting the item in place", last.unwrap_or_else(|| std::io::Error::other("rename failed"))))
        })();
        if r.is_err() {
            let _ = std::fs::remove_dir_all(&tmp);
        }
        r.map(|()| id)
    }

    /// The item `id`, if it is there.
    pub fn item(&self, id: &str) -> Option<Item> {
        let dir = self.item_dir(id);
        let request: Request = read_json(&dir.join("request.json"))?;
        let state: Option<State> = read_json(&dir.join("state.json"));
        let priority = std::fs::read_to_string(dir.join("priority"))
            .ok()
            .and_then(|t| t.trim().parse::<u8>().ok())
            .filter(|p| (1..=3).contains(p))
            .unwrap_or(request.priority);
        Some(Item {
            request,
            state,
            held: dir.join("hold").exists(),
            cancel_requested: dir.join("cancel").exists(),
            priority,
            position: None,
        })
    }

    /// Every item, oldest first, with its place in the order of starting where the queue worker has ranked it.
    pub fn items(&self) -> Vec<Item> {
        let Ok(rd) = std::fs::read_dir(self.items_dir()) else { return vec![] };
        let mut ids: Vec<String> = rd
            .filter_map(Result::ok)
            .filter_map(|e| e.file_name().into_string().ok())
            .filter(|n| !n.starts_with('.'))
            .collect();
        ids.sort();
        let order: Vec<String> = read_json(&self.root.join("order.json")).unwrap_or_default();
        ids.iter()
            .filter_map(|id| self.item(id))
            .map(|mut i| {
                i.position = order.iter().position(|o| *o == i.request.id).map(|p| p + 1);
                i
            })
            .collect()
    }

    fn marker(&self, id: &str, name: &str, content: Option<&str>) -> bool {
        let dir = self.item_dir(id);
        if !dir.join("request.json").exists() {
            return false;
        }
        match content {
            Some(c) => std::fs::write(dir.join(name), c).is_ok(),
            None => {
                let _ = std::fs::remove_file(dir.join(name));
                true
            }
        }
    }

    /// Asks for `id` to be cancelled: a queued item is cancelled at the next look of the queue worker, a running one has its solver killed.
    pub fn cancel(&self, id: &str) -> bool {
        self.marker(id, "cancel", Some(""))
    }

    pub fn hold(&self, id: &str) -> bool {
        self.marker(id, "hold", Some(""))
    }

    pub fn release(&self, id: &str) -> bool {
        self.marker(id, "hold", None)
    }

    pub fn set_priority(&self, id: &str, priority: u8) -> bool {
        (1..=3).contains(&priority) && self.marker(id, "priority", Some(&priority.to_string()))
    }

    /// The queue is paused: the queue worker finishes what it runs and starts nothing.
    pub fn pause(&self) -> bool {
        std::fs::create_dir_all(&self.root).is_ok() && std::fs::write(self.root.join("pause"), b"").is_ok()
    }

    pub fn resume(&self) {
        let _ = std::fs::remove_file(self.root.join("pause"));
    }

    pub fn paused(&self) -> bool {
        self.root.join("pause").exists()
    }

    /// Removes an item (the results of a finished one included); false when it is still open.
    pub fn remove(&self, id: &str) -> bool {
        match self.item(id) {
            Some(i) if i.is_open() => false,
            Some(_) => std::fs::remove_dir_all(self.item_dir(id)).is_ok(),
            None => true,
        }
    }

    /// Removes the finished items; returns how many.
    pub fn clear_finished(&self) -> usize {
        self.items().iter().filter(|i| !i.is_open()).filter(|i| self.remove(&i.request.id)).count()
    }

    /// Whether a queue worker holds the lease now.
    pub fn runner_alive(&self) -> bool {
        age_of(&self.lease_path()).is_some_and(|a| a < LEASE_DEAD)
    }

    /// Starts a queue worker (`program --worker <workspace>`) unless one holds the lease. The program runs on its own (`process::spawn_detached`:
    /// no console, no standard streams, no inherited handles and, on Windows, outside the job object of the program that starts it), so closing
    /// the window or the terminal leaves it running. The queue worker takes the lease itself: when two are started at once one of them ends.
    pub fn ensure_runner(&self, program: &Path) -> Result<(), Error> {
        if self.runner_alive() {
            return Ok(());
        }
        let workspace = self.root.parent().unwrap_or(&self.root);
        crate::process::spawn_detached(program, &["--worker".into(), workspace.as_os_str().to_owned()])
            .map_err(|e| Error::Solver(format!("could not start the queue worker {}: {e}", program.display())))
    }
}

// ------------------------------------------------------------------------------------------------------------------ the queue worker

/// What the queue worker is told.
pub struct RunnerOptions {
    /// the program an item's solver child is (`program --job FOLDER`)
    pub worker: PathBuf,
    /// the queue worker ends after this long with nothing queued and nothing running
    pub idle_exit: Duration,
    /// how often it looks at the folder
    pub tick: Duration,
    /// a waiting item rises one priority level for each time this passes
    pub age: Duration,
    pub probe: MemoryProbe,
    pub history: std::sync::Arc<History>,
}

impl RunnerOptions {
    /// The settings of a real queue worker: the program itself as the solver, the machine's memory, an end after 30 s idle (`NYX_QUEUE_IDLE_SECS`).
    pub fn standard(workspace: &Path) -> RunnerOptions {
        let secs = std::env::var("NYX_QUEUE_IDLE_SECS").ok().and_then(|s| s.parse::<f64>().ok()).unwrap_or(30.0);
        let history = std::sync::Arc::new(History::default());
        history.set_file(&workspace.join("solve_history.jsonl"));
        RunnerOptions {
            worker: std::env::current_exe().unwrap_or_else(|_| PathBuf::from("nyx")),
            idle_exit: Duration::from_secs_f64(secs),
            tick: Duration::from_millis(250),
            age: Duration::from_secs(60),
            probe: std::sync::Arc::new(crate::process::system_memory),
            history,
        }
    }
}

struct Lease {
    path: PathBuf,
    token: String,
}

impl Lease {
    /// Takes the lease when nobody holds a live one.
    fn take(q: &SharedQueue) -> Option<Lease> {
        std::fs::create_dir_all(&q.root).ok()?;
        let path = q.lease_path();
        if path.exists() {
            if age_of(&path).is_some_and(|a| a < LEASE_DEAD) {
                return None;
            }
            let _ = std::fs::remove_file(&path);
        }
        let token = format!("{}-{}", std::process::id(), now_ms());
        let mut f = std::fs::OpenOptions::new().write(true).create_new(true).open(&path).ok()?;
        use std::io::Write;
        f.write_all(
            json!({ "pid": std::process::id(), "token": token, "started_ms": now_ms() }).to_string().as_bytes(),
        )
        .ok()?;
        drop(f);
        // two programs that saw the same dead lease may both have removed it: only one create_new won, and it is the one whose token is there
        let held: Option<Value> = read_json(&path);
        (held.and_then(|v| v["token"].as_str().map(|t| t == token)) == Some(true)).then_some(Lease { path, token })
    }

    /// Touches the lease; false when it is no longer ours.
    fn beat(&self) -> bool {
        let ours =
            read_json::<Value>(&self.path).and_then(|v| v["token"].as_str().map(|t| t == self.token)) == Some(true);
        if ours {
            if let Ok(f) = std::fs::OpenOptions::new().write(true).open(&self.path) {
                let _ = f.set_modified(SystemTime::now());
            }
        }
        ours
    }

    fn release(&self) {
        if read_json::<Value>(&self.path).and_then(|v| v["token"].as_str().map(|t| t == self.token)) == Some(true) {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

fn write_state(dir: &Path, s: &State) {
    let _ = write_atomic(&dir.join("state.json"), json!(s).to_string().as_bytes());
}

/// Runs the queue of `workspace` until it has been idle for `idle_exit`, or another queue worker holds the lease. Returns at once (and `Ok`) when
/// one does. This is the whole of `nyx --worker <workspace>`.
#[cfg(feature = "hdf5")]
pub fn run_queue(workspace: &Path, o: &RunnerOptions) -> Result<(), String> {
    let q = SharedQueue::new(workspace);
    let Some(mut lease) = Lease::take(&q) else { return Ok(()) };
    // what the last queue worker left running is not running: it died with the lease it held, and its solver child went with it
    for i in q.items() {
        if matches!(i.state.as_ref().map(|s| s.state.as_str()), Some("running")) {
            let mut s = i.state.clone().unwrap_or_default();
            s.state = "failed".into();
            s.ended_ms = Some(now_ms());
            s.error = Some(
                json!({ "code": "solver", "message": "the queue worker that was running this item stopped (it was ended, or the machine restarted); nothing else was affected. Submit it again" }),
            );
            write_state(&q.item_dir(&i.request.id), &s);
        }
    }
    let mut served: BTreeMap<String, u64> = BTreeMap::new();
    let mut tick: u64 = 0;
    let mut idle_since = std::time::Instant::now();
    loop {
        if !lease.beat() {
            return Ok(());
        }
        let items = q.items();
        // cancelled before they started
        for i in items.iter().filter(|i| i.state.is_none() && i.cancel_requested) {
            write_state(
                &q.item_dir(&i.request.id),
                &State {
                    state: "cancelled".into(),
                    ended_ms: Some(now_ms()),
                    runner_pid: std::process::id(),
                    ..State::default()
                },
            );
        }
        let mut waiting: Vec<&Item> =
            items.iter().filter(|i| i.state.is_none() && !i.cancel_requested && !i.held).collect();
        let now = now_ms();
        waiting.sort_by_key(|i| {
            let waited = Duration::from_millis(now.saturating_sub(i.request.submitted_ms));
            (
                effective(i.priority, waited, o.age),
                served.get(&i.request.owner).copied().unwrap_or(0),
                i.request.id.clone(),
            )
        });
        let ids: Vec<&str> = waiting.iter().map(|i| i.request.id.as_str()).collect();
        let _ = write_atomic(&q.root.join("order.json"), json!(ids).to_string().as_bytes());
        // forget what nobody collected
        for i in items.iter().filter(|i| i.state.as_ref().is_some_and(State::is_finished)) {
            let ended = i.state.as_ref().and_then(|s| s.ended_ms).unwrap_or(now);
            if now.saturating_sub(ended) > KEEP_FINISHED.as_millis() as u64 {
                let _ = std::fs::remove_dir_all(q.item_dir(&i.request.id));
            }
        }
        match waiting.first() {
            Some(head) if !q.paused() => {
                tick += 1;
                served.insert(head.request.owner.clone(), tick);
                run_item(&q, head, o, &lease);
                idle_since = std::time::Instant::now();
            }
            Some(_) => idle_since = std::time::Instant::now(),
            None => {
                if idle_since.elapsed() >= o.idle_exit {
                    lease.release();
                    // an item that was written while the lease was being given up is taken by a queue worker of its own
                    if q.items().iter().any(|i| i.state.is_none() && !i.cancel_requested && !i.held) {
                        if let Some(again) = Lease::take(&q) {
                            lease = again;
                            idle_since = std::time::Instant::now();
                            continue;
                        }
                    }
                    return Ok(());
                }
            }
        }
        std::thread::sleep(o.tick);
    }
}

/// Runs one item to its end and records the outcome in the item.
#[cfg(feature = "hdf5")]
fn run_item(q: &SharedQueue, item: &Item, o: &RunnerOptions, lease: &Lease) {
    use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
    use std::sync::Arc;
    let dir = q.item_dir(&item.request.id);
    let progress = femcore::progress::Progress::new();
    let pid = Arc::new(AtomicU32::new(0));
    let stop = Arc::new(AtomicBool::new(false));
    let started = now_ms();
    let mk = |state: &str, p: &femcore::progress::Progress, child: u32| {
        let st = p.status();
        State {
            state: state.into(),
            stage: st.stage,
            stage_seconds: st.stage_seconds,
            fraction: st.fraction,
            log: p.log(),
            started_ms: Some(started),
            child_pid: (child != 0).then_some(child),
            runner_pid: std::process::id(),
            ..State::default()
        }
    };
    write_state(&dir, &mk("running", &progress, 0));
    // the one writer of state.json while the item runs: it shows the stage, keeps the lease alive and turns the cancel marker into a cancel
    let watcher = {
        let (progress, pid, stop, dir) = (progress.clone(), pid.clone(), stop.clone(), dir.clone());
        let lease = Lease { path: lease.path.clone(), token: lease.token.clone() };
        std::thread::spawn(move || {
            let mut last = String::new();
            let mut written = std::time::Instant::now();
            while !stop.load(Ordering::SeqCst) {
                lease.beat();
                if dir.join("cancel").exists() {
                    progress.cancel();
                }
                let status = progress.status();
                let st = State {
                    state: "running".into(),
                    stage: status.stage,
                    stage_seconds: status.stage_seconds,
                    fraction: status.fraction,
                    log: progress.log(),
                    started_ms: Some(started),
                    child_pid: Some(pid.load(Ordering::SeqCst)).filter(|p| *p != 0),
                    runner_pid: std::process::id(),
                    ..State::default()
                };
                // the seconds in a stage change at every look: write when something else did, and at least once a second
                let key = format!("{}|{:?}|{}|{:?}", st.stage, st.fraction, st.log.len(), st.child_pid);
                if key != last || written.elapsed() >= Duration::from_secs(1) {
                    let _ = write_atomic(&dir.join("state.json"), json!(st).to_string().as_bytes());
                    last = key;
                    written = std::time::Instant::now();
                }
                std::thread::sleep(Duration::from_millis(100));
            }
        })
    };
    let outcome: Result<Value, Error> = (|| {
        if let Some(est) = &item.request.estimate {
            crate::estimate::admit(est, item.request.force, &o.probe, &progress, item.request.memory_wait_timeout_s)?;
        }
        let base = progress.log();
        let end = crate::worker::run_child(
            &o.worker,
            &dir,
            &mut |s| {
                let log = s.get("log").and_then(|l| serde_json::from_value::<Vec<String>>(l.clone()).ok()).map(|l| {
                    let mut all = base.clone();
                    all.extend(l);
                    all
                });
                progress.mirror(
                    s.get("stage").and_then(Value::as_str).unwrap_or(""),
                    s.get("stage_seconds").and_then(Value::as_f64).unwrap_or(0.0),
                    s.get("fraction").and_then(Value::as_f64),
                    log,
                );
            },
            &mut |p| pid.store(p, Ordering::SeqCst),
            &|| progress.is_cancelled(),
        )?;
        end.done.ok_or_else(|| Error::Solver(crate::worker::stopped_words(&end.status, &end.stderr)))
    })();
    stop.store(true, Ordering::SeqCst);
    let _ = watcher.join();
    let mut s = mk("failed", &progress, pid.load(Ordering::SeqCst));
    s.ended_ms = Some(now_ms());
    match outcome {
        Ok(d) => {
            s.has_results = d["has_results"].as_bool() == Some(true);
            s.peak_working_set = d["peak_working_set"].as_u64();
            s.peak_commit = d["peak_commit"].as_u64();
            match d.get("error").filter(|e| !e.is_null()) {
                None => s.state = "done".into(),
                Some(e) => s.error = Some(e.clone()),
            }
            // what the machine learned: the same row a front-end would have added, added here so it is there whoever was watching
            if let (Some(ws), true, Some(sz)) = (s.peak_working_set, s.state == "done", d.get("size")) {
                if let (Some(chain), Some(dim), Some(cells), Some(dofs)) =
                    (sz["chain"].as_str(), sz["dim"].as_u64(), sz["cells"].as_u64(), sz["dofs"].as_u64())
                {
                    o.history.add(crate::estimate::Row {
                        chain: chain.into(),
                        dim: dim as u8,
                        cells,
                        dofs,
                        peak_working_set: ws,
                        peak_commit: s.peak_commit.unwrap_or(0),
                        seconds: sz["seconds"].as_f64().unwrap_or(0.0),
                    });
                    o.history.flush();
                }
            }
        }
        Err(Error::Cancelled) => s.state = "cancelled".into(),
        Err(e) => s.error = Some(json!({ "code": e.code(), "message": e.to_string() })),
    }
    write_state(&dir, &s);
}

/// Hands a solve to the queue of the workspace and waits for it there, as one job of this engine: writes the item, makes sure a queue worker is
/// running, shows where the item is (its place in the queue, then the queue worker's stage) as the job's progress, turns a cancel of the job into
/// the cancel marker, and when the item ends reads its results like a worker's. `id_slot` and `own` tell the engine which item is its own.
#[cfg(feature = "hdf5")]
pub(crate) fn run_shared(
    q: &SharedQueue,
    program: &Path,
    n: &NewItem,
    spec: &crate::worker::RemoteSpec,
    progress: &std::sync::Arc<femcore::progress::Progress>,
    id_slot: &std::sync::Mutex<Option<String>>,
    own: &std::sync::Mutex<std::collections::BTreeSet<String>>,
) -> Result<crate::worker::RemoteDone, Error> {
    if progress.is_cancelled() {
        return Err(Error::Cancelled);
    }
    let id = q.submit(n, spec)?;
    *id_slot.lock().expect("shared id") = Some(id.clone());
    own.lock().expect("own items").insert(id.clone());
    let forget = |id: &str| {
        own.lock().expect("own items").remove(id);
    };
    if let Err(e) = q.ensure_runner(program) {
        q.cancel(&id);
        q.remove(&id);
        forget(&id);
        return Err(e);
    }
    let base = progress.log();
    let (since, mut cancel_sent, mut last_start) = (std::time::Instant::now(), false, std::time::Instant::now());
    let finished = loop {
        if progress.is_cancelled() && !cancel_sent {
            q.cancel(&id);
            cancel_sent = true;
        }
        // a queue worker that ended (it was killed, or it ended just as the item was written) is replaced by the next one who looks
        if !q.runner_alive() && last_start.elapsed() > Duration::from_secs(2) {
            let _ = q.ensure_runner(program);
            last_start = std::time::Instant::now();
        }
        let Some(item) = q.item(&id) else {
            forget(&id);
            return Err(Error::Solver("the item was removed from the queue by someone else before it finished".into()));
        };
        match &item.state {
            None => {
                let text = match (item.held, item.position) {
                    (true, _) => "held in the queue: release it to run".to_string(),
                    (_, Some(p)) => format!("queued: place {p} in the queue of the workspace"),
                    _ => "queued: waiting for the queue worker".to_string(),
                };
                progress.mirror(&text, since.elapsed().as_secs_f64(), None, None);
            }
            Some(s) if !s.is_finished() => {
                let mut log = base.clone();
                log.extend(s.log.iter().cloned());
                progress.mirror(&s.stage, s.stage_seconds, s.fraction, Some(log));
            }
            Some(s) => break s.clone(),
        }
        std::thread::sleep(Duration::from_millis(100));
    };
    let dir = q.item_dir(&id);
    let result = if finished.state == "cancelled" {
        Err(Error::Cancelled)
    } else {
        let d = json!({
            "error": finished.error, "has_results": finished.has_results,
            "peak_working_set": finished.peak_working_set, "peak_commit": finished.peak_commit,
        });
        crate::worker::read_done(&d, &dir, Some(progress)).map(|mut done| {
            done.run.shared = true;
            done
        })
    };
    q.remove(&id);
    forget(&id);
    result
}

#[cfg(not(feature = "hdf5"))]
pub(crate) fn run_shared(
    _q: &SharedQueue,
    _program: &Path,
    _n: &NewItem,
    _spec: &crate::worker::RemoteSpec,
    _progress: &std::sync::Arc<femcore::progress::Progress>,
    _id_slot: &std::sync::Mutex<Option<String>>,
    _own: &std::sync::Mutex<std::collections::BTreeSet<String>>,
) -> Result<crate::worker::RemoteDone, Error> {
    Err(Error::Unsupported(
        "this build of nyx has no HDF5 support, which the shared queue needs to hand results back; solve in this process"
            .into(),
    ))
}

// ------------------------------------------------------------------------------------------------------------------ the command line

fn ago(ms: u64) -> String {
    let s = now_ms().saturating_sub(ms) / 1000;
    if s < 90 {
        format!("{s} s")
    } else if s < 5400 {
        format!("{} min", s / 60)
    } else {
        format!("{} h", s / 3600)
    }
}

/// `nyx queue list|cancel|hold|release|priority|pause|resume|clear WORKSPACE [ID [N]]`: what a person at a terminal does with the queue of a
/// workspace, without a window. Returns the exit code.
pub fn cli(args: &[String]) -> i32 {
    let usage = "usage: nyx queue submit WORKSPACE PROJECT.nyx.json [PHYSICS...] [--priority 1|2|3] [--hold]
       nyx queue list WORKSPACE [--json]
       nyx queue cancel|hold|release WORKSPACE ID
       nyx queue priority WORKSPACE ID 1|2|3
       nyx queue pause|resume|clear WORKSPACE";
    let [cmd, workspace, rest @ ..] = args else {
        eprintln!("{usage}");
        return 2;
    };
    let q = SharedQueue::new(Path::new(workspace));
    let with_id = |f: &dyn Fn(&str) -> bool, what: &str| match rest.first() {
        Some(id) if f(id) => {
            println!("{what} {id}");
            0
        }
        Some(id) => {
            eprintln!("nyx: no item {id} in the queue of {workspace}");
            1
        }
        None => {
            eprintln!("{usage}");
            2
        }
    };
    match cmd.as_str() {
        "submit" => cli_submit(workspace, rest),
        "list" => {
            let items = q.items();
            if rest.first().map(String::as_str) == Some("--json") {
                println!("{}", json!(items.iter().map(|i| i.to_json(&q)).collect::<Vec<_>>()));
                return 0;
            }
            println!(
                "queue worker: {}{}",
                if q.runner_alive() { "running" } else { "not running" },
                if q.paused() { ", queue paused" } else { "" }
            );
            if items.is_empty() {
                println!("(no items)");
            }
            for i in &items {
                let stage = i.state.as_ref().map(|s| s.stage.clone()).unwrap_or_default();
                println!(
                    "{:<26} {:<9} {:<4} p{} {:<24} {} ago{}",
                    i.request.id,
                    i.state_name(),
                    i.request.owner,
                    i.priority,
                    i.request.label,
                    ago(i.request.submitted_ms),
                    if stage.is_empty() { String::new() } else { format!("  {stage}") }
                );
            }
            0
        }
        "cancel" => with_id(&|id| q.cancel(id), "cancel requested for"),
        "hold" => with_id(&|id| q.hold(id), "held"),
        "release" => with_id(&|id| q.release(id), "released"),
        "priority" => match rest.get(1).and_then(|p| p.parse::<u8>().ok()) {
            Some(p) => with_id(&|id| q.set_priority(id, p), "priority changed for"),
            None => {
                eprintln!("{usage}");
                2
            }
        },
        "pause" => {
            q.pause();
            println!("queue paused");
            0
        }
        "resume" => {
            q.resume();
            println!("queue resumed");
            0
        }
        "clear" => {
            println!("removed {} finished item(s)", q.clear_finished());
            0
        }
        _ => {
            eprintln!("{usage}");
            2
        }
    }
}

/// `nyx queue submit WORKSPACE PROJECT [PHYSICS...]`: writes a solve of a project file into the queue and returns; the queue worker runs it.
fn cli_submit(workspace: &str, rest: &[String]) -> i32 {
    let mut opts = crate::queue::QueueOpts::default();
    let (mut project, mut targets) = (None::<&String>, vec![]);
    let mut it = rest.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--hold" => opts.hold = true,
            "--priority" => match it.next().and_then(|p| p.parse::<u8>().ok()).filter(|p| (1..=3).contains(p)) {
                Some(p) => opts.priority = Some(p),
                None => {
                    eprintln!("nyx: --priority is 1, 2 or 3");
                    return 2;
                }
            },
            _ if project.is_none() => project = Some(a),
            _ => targets.push(a.clone()),
        }
    }
    let Some(project) = project else {
        eprintln!("nyx: name the project file to solve");
        return 2;
    };
    let mut e = crate::engine::Engine::new();
    e.set_owner("cli");
    if !e.set_self_as_worker(false) {
        eprintln!("nyx: this system does not say which program is running, so no worker can be started");
        return 1;
    }
    e.set_shared_queue(Path::new(workspace), false);
    let done = e
        .call("project_open", &json!({ "path": project }))
        .and_then(|_| e.queue_detached((!targets.is_empty()).then_some(targets), &opts));
    match done {
        Ok(v) => {
            println!("queued {}", v["shared_id"].as_str().unwrap_or(""));
            0
        }
        Err(err) => {
            eprintln!("nyx: {err}");
            1
        }
    }
}
