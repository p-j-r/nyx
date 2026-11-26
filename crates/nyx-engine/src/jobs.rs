// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The job runner and its queue (productivity plan P2 and Q1, ADR 85 and 87): a mesh build, a solve, a convergence check, opening or writing a
//! results file is a job with a `Progress` handle, so it can be watched (stage, elapsed) and cancelled by the window, a script or the MCP server
//! alike. Jobs run in *lanes*: the heavy lane (mesh, solve, convergence check, open: they hold every core and the memory) runs one job at a
//! time, the light lane (writing results) beside it. A job submitted while its lane is busy waits in the queue; the runner starts the next one
//! itself when a job ends, so no front-end has to call anything for the queue to move. The order is by priority (a person's request before an
//! agent's single request before a sweep, with ageing so nothing waits forever), then round robin between owners, then submission order. The
//! runner knows nothing of windows: a front-end asks `info`, `cancel`, `try_take` (non-blocking, for a repaint loop) or `wait` (blocking).

use crate::meshing::MeshedModel;
use crate::project::Project;
use femcore::progress::Progress;
use nyx_physics::{Error, SolveResult};
use serde::Serialize;
use serde_json::Value;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

pub type JobId = u64;

/// What a finished job hands back. The engine applies it to its state (`Engine::apply`).
pub enum JobOutput {
    Mesh(Arc<MeshedModel>),
    /// the mesh the solve used (built inside the job when none was current) and one result per physics, in solve order. When a physics of the
    /// chain failed, `failed` holds why and `results` the physics that were solved before it (they stay valid results)
    Solved {
        model: Arc<MeshedModel>,
        results: Vec<SolveResult>,
        /// the chain the job solved (its targets), and the error that stopped it
        failed: Option<Error>,
        /// what the worker process reported, when the solve ran in one (`worker.rs`)
        worker: Option<WorkerRun>,
        /// how the solve was started when it was a point of a warm sweep (`reuse.rs`)
        reuse: Option<Value>,
        /// whether mesh and field arrays were kept (`false` when `save_fields: false`)
        fields_saved: bool,
        /// outputs of any read-only `post_analyses` run on full in-memory fields before dropping
        post_analyses: Vec<Value>,
    },
    /// the answer of a mesh convergence check
    Convergence(Value),
    /// a results file read, its project meshed again and checked against the stored mesh (applied by replacing the project)
    Opened(Box<Opened>),
    /// a results file written
    Exported(Value),
}

/// What a worker process measured of its own solve: the most memory it used, in bytes (the in-process job cannot know: the process is shared).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct WorkerRun {
    pub peak_working_set: Option<u64>,
    pub peak_commit: Option<u64>,
    /// the run was done by the queue worker of the workspace, which has added its row to the history itself (Q3)
    #[serde(skip)]
    pub shared: bool,
}

/// What opening a results file gives the engine to apply.
pub struct Opened {
    pub project: Project,
    pub model: Arc<MeshedModel>,
    /// one result per physics of the file
    pub results: Vec<(String, SolveResult)>,
    pub path: String,
    /// convergence checks stored in `/provenance/validated_json`, keyed by problem key
    pub checks: BTreeMap<String, crate::validated::Check>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case", tag = "state", content = "message")]
pub enum JobState {
    /// submitted and waiting for its lane
    Queued,
    Running,
    Done,
    Cancelled,
    Failed(String),
}

impl JobState {
    /// Queued or running: the job has not finished.
    pub fn is_open(&self) -> bool {
        matches!(self, JobState::Queued | JobState::Running)
    }
}

/// The resource class of a job. The heavy lane runs one job at a time; the light lane (I/O) runs beside it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Lane {
    Heavy,
    Light,
    /// a job that hands its work to the queue of the workspace (Q3) and waits for it there: it holds no core and no memory of this process, so
    /// many of them wait side by side while the queue worker runs one solve at a time
    Shared,
}

/// How a job is submitted: its name, lane, owner and priority.
#[derive(Clone, Debug)]
pub struct Submit {
    pub label: String,
    pub lane: Lane,
    /// who asked: `gui`, `ai` or `cli`; the queue serves owners in turn within one priority
    pub owner: String,
    /// 1 a person's request now, 2 an agent's single request (the default), 3 sweeps and batches; lower runs first
    pub priority: u8,
    /// what the request is, by kind and targets (`solve ems`): the key a newer request can replace a queued one by
    pub key: Option<String>,
    /// cancel a job of the same owner and key that is still queued when this one is submitted (a new Compute of the same target replaces the one
    /// that has not started); every job keeps its key, so a request made without this flag can still be replaced by one made with it
    pub supersede: bool,
    /// start held: it runs only after `release`
    pub hold: bool,
}

impl Submit {
    pub fn new(label: &str, lane: Lane) -> Submit {
        Submit { label: label.into(), lane, owner: "cli".into(), priority: 2, key: None, supersede: false, hold: false }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct JobInfo {
    pub id: JobId,
    pub label: String,
    #[serde(flatten)]
    pub state: JobState,
    pub lane: Lane,
    pub owner: String,
    pub priority: u8,
    /// held jobs wait until released
    pub held: bool,
    /// place in the queue of its lane (1 runs next) while the job is queued
    pub position: Option<usize>,
    /// seconds since the job was submitted and not started (while queued, and what it waited when it ran)
    pub queued_seconds: f64,
    /// the stage the job is in now (plain words), empty before the first
    pub stage: String,
    /// seconds the job has been in that stage: what a person waits for when the stage cannot be interrupted (a factorisation)
    pub stage_seconds: f64,
    pub fraction: Option<f64>,
    /// a cancel was asked for and the job has not stopped yet; it stops at its next check (a factorisation finishes first)
    pub cancel_requested: bool,
    /// seconds the job has run (0 while queued)
    pub seconds: f64,
}

type Work = Box<dyn FnOnce(Arc<Progress>) -> Result<JobOutput, Error> + Send>;

struct Inner {
    state: JobState,
    outcome: Option<Result<JobOutput, Error>>,
    seconds: f64,
    priority: u8,
    held: bool,
    started: Option<Instant>,
    work: Option<Work>,
}

struct Slot {
    label: String,
    lane: Lane,
    owner: String,
    key: Option<String>,
    progress: Arc<Progress>,
    submitted: Instant,
    inner: Mutex<Inner>,
    done: Condvar,
}

/// The queue's own settings and bookkeeping.
struct Sched {
    paused: bool,
    /// jobs that may run at once, per lane
    limit_heavy: usize,
    limit_light: usize,
    limit_shared: usize,
    /// a queued job rises one priority level after this long
    age: Duration,
    /// when each owner was last served (a counter), for the round robin
    served: BTreeMap<String, u64>,
    tick: u64,
}

struct Shared {
    next: AtomicU64,
    slots: Mutex<BTreeMap<JobId, Arc<Slot>>>,
    sched: Mutex<Sched>,
}

#[derive(Clone)]
pub struct JobRunner {
    shared: Arc<Shared>,
}

impl Default for JobRunner {
    fn default() -> Self {
        JobRunner::new()
    }
}

/// The priority a queued job competes with: its own, lowered one level for each `age` it has waited, never below 1.
fn effective(priority: u8, waited: Duration, age: Duration) -> u8 {
    let levels = if age.is_zero() { 0 } else { (waited.as_secs_f64() / age.as_secs_f64()) as u64 };
    (priority as u64).saturating_sub(levels).max(1) as u8
}

impl Shared {
    fn slot(&self, id: JobId) -> Option<Arc<Slot>> {
        self.slots.lock().expect("job table").get(&id).cloned()
    }

    /// The queued, not held jobs of `lane` in the order they will start.
    fn waiting(&self, lane: Lane, sched: &Sched) -> Vec<(JobId, Arc<Slot>)> {
        let all: Vec<(JobId, Arc<Slot>)> =
            self.slots.lock().expect("job table").iter().map(|(i, s)| (*i, s.clone())).collect();
        let mut out: Vec<(JobId, Arc<Slot>)> = all
            .into_iter()
            .filter(|(_, s)| {
                let g = s.inner.lock().expect("job slot");
                s.lane == lane && g.state == JobState::Queued && !g.held
            })
            .collect();
        let key = |(id, s): &(JobId, Arc<Slot>)| {
            let p = s.inner.lock().expect("job slot").priority;
            (effective(p, s.submitted.elapsed(), sched.age), sched.served.get(&s.owner).copied().unwrap_or(0), *id)
        };
        out.sort_by_key(key);
        out
    }

    fn running(&self, lane: Lane) -> usize {
        let all: Vec<Arc<Slot>> = self.slots.lock().expect("job table").values().cloned().collect();
        all.iter().filter(|s| s.lane == lane && s.inner.lock().expect("job slot").state == JobState::Running).count()
    }
}

/// Starts queued jobs while a lane has room. Called whenever something that can free a lane or add work happens.
fn pump(shared: &Arc<Shared>) {
    let mut sched = shared.sched.lock().expect("job queue");
    if sched.paused {
        return;
    }
    loop {
        let mut started = false;
        for lane in [Lane::Heavy, Lane::Light, Lane::Shared] {
            let limit = match lane {
                Lane::Heavy => sched.limit_heavy,
                Lane::Light => sched.limit_light,
                Lane::Shared => sched.limit_shared,
            };
            if shared.running(lane) >= limit {
                continue;
            }
            let Some((_, slot)) = shared.waiting(lane, &sched).into_iter().next() else { continue };
            let work = {
                let mut g = slot.inner.lock().expect("job slot");
                if g.state != JobState::Queued {
                    continue;
                }
                g.state = JobState::Running;
                g.started = Some(Instant::now());
                g.work.take()
            };
            sched.tick += 1;
            let tick = sched.tick;
            sched.served.insert(slot.owner.clone(), tick);
            if let Some(work) = work {
                spawn(shared.clone(), slot, work);
            }
            started = true;
        }
        if !started {
            break;
        }
    }
}

fn spawn(shared: Arc<Shared>, slot: Arc<Slot>, work: Work) {
    let progress = slot.progress.clone();
    std::thread::spawn(move || {
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(progress))).unwrap_or_else(|p| {
            let text = p.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| p.downcast_ref::<String>().cloned());
            Err(Error::Solver(format!("the job stopped unexpectedly: {}", text.unwrap_or_else(|| "panic".into()))))
        });
        {
            let mut g = slot.inner.lock().expect("job slot");
            g.seconds = g.started.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0);
            g.state = match &r {
                Ok(_) => JobState::Done,
                Err(Error::Cancelled) => JobState::Cancelled,
                Err(e) => JobState::Failed(e.to_string()),
            };
            g.outcome = Some(r);
            slot.done.notify_all();
        }
        // the lane is free: the next queued job starts without anyone asking
        pump(&shared);
    });
}

impl JobRunner {
    pub fn new() -> JobRunner {
        JobRunner {
            shared: Arc::new(Shared {
                next: AtomicU64::new(0),
                slots: Mutex::new(BTreeMap::new()),
                sched: Mutex::new(Sched {
                    paused: false,
                    limit_heavy: 1,
                    limit_light: 2,
                    limit_shared: 256,
                    age: Duration::from_secs(60),
                    served: BTreeMap::new(),
                    tick: 0,
                }),
            }),
        }
    }

    /// Submits `work` to the heavy lane with the default owner and priority.
    pub fn start<F>(&self, label: &str, work: F) -> JobId
    where
        F: FnOnce(Arc<Progress>) -> Result<JobOutput, Error> + Send + 'static,
    {
        self.submit(Submit::new(label, Lane::Heavy), work)
    }

    /// Submits `work`. It starts at once when its lane is free, otherwise it waits in the queue. The closure gets the job's `Progress`; a panic in
    /// it becomes a failed job, not a dead process.
    pub fn submit<F>(&self, spec: Submit, work: F) -> JobId
    where
        F: FnOnce(Arc<Progress>) -> Result<JobOutput, Error> + Send + 'static,
    {
        let id = self.shared.next.fetch_add(1, Ordering::SeqCst) + 1;
        if let (true, Some(key)) = (spec.supersede, &spec.key) {
            // the request this one replaces has not started: it never runs
            let old: Vec<(JobId, Arc<Slot>)> = self
                .shared
                .slots
                .lock()
                .expect("job table")
                .iter()
                .filter(|(_, s)| s.owner == spec.owner && s.key.as_ref() == Some(key))
                .map(|(i, s)| (*i, s.clone()))
                .collect();
            for (i, _) in old {
                self.cancel_queued(i);
            }
        }
        let slot = Arc::new(Slot {
            label: spec.label,
            lane: spec.lane,
            owner: spec.owner,
            key: spec.key,
            progress: Progress::new(),
            submitted: Instant::now(),
            inner: Mutex::new(Inner {
                state: JobState::Queued,
                outcome: None,
                seconds: 0.0,
                priority: spec.priority.max(1),
                held: spec.hold,
                started: None,
                work: Some(Box::new(work)),
            }),
            done: Condvar::new(),
        });
        self.shared.slots.lock().expect("job table").insert(id, slot);
        pump(&self.shared);
        id
    }

    fn slot(&self, id: JobId) -> Option<Arc<Slot>> {
        self.shared.slot(id)
    }

    pub fn info(&self, id: JobId) -> Option<JobInfo> {
        let s = self.slot(id)?;
        let position = {
            let queued = s.inner.lock().expect("job slot").state == JobState::Queued;
            if queued {
                let sched = self.shared.sched.lock().expect("job queue");
                self.shared.waiting(s.lane, &sched).iter().position(|(i, _)| *i == id).map(|p| p + 1)
            } else {
                None
            }
        };
        let g = s.inner.lock().expect("job slot");
        let st = s.progress.status();
        let waited = match g.started {
            Some(t) => t.duration_since(s.submitted).as_secs_f64(),
            None if g.state == JobState::Queued => s.submitted.elapsed().as_secs_f64(),
            None => 0.0,
        };
        Some(JobInfo {
            id,
            label: s.label.clone(),
            state: g.state.clone(),
            lane: s.lane,
            owner: s.owner.clone(),
            priority: g.priority,
            held: g.held,
            position,
            queued_seconds: waited,
            stage: st.stage,
            stage_seconds: st.stage_seconds,
            fraction: st.fraction,
            cancel_requested: st.cancelled,
            seconds: match g.state {
                JobState::Running => g.started.map(|t| t.elapsed().as_secs_f64()).unwrap_or(0.0),
                _ => g.seconds,
            },
        })
    }

    /// The progress handle of a job: the stage it is in, the log so far, the cancel flag. A front-end that shows a bar reads it.
    pub fn progress(&self, id: JobId) -> Option<Arc<Progress>> {
        self.slot(id).map(|s| s.progress.clone())
    }

    /// All jobs, oldest first.
    pub fn list(&self) -> Vec<JobInfo> {
        let ids: Vec<JobId> = self.shared.slots.lock().expect("job table").keys().copied().collect();
        ids.into_iter().filter_map(|i| self.info(i)).collect()
    }

    /// The jobs that have not finished (queued or running).
    pub fn running(&self) -> Vec<JobId> {
        self.list().into_iter().filter(|j| j.state.is_open()).map(|j| j.id).collect()
    }

    /// Cancels a job that has not started: it never runs and finishes as `Cancelled` at once.
    fn cancel_queued(&self, id: JobId) -> bool {
        let Some(s) = self.slot(id) else { return false };
        let mut g = s.inner.lock().expect("job slot");
        if g.state != JobState::Queued {
            return false;
        }
        g.state = JobState::Cancelled;
        g.work = None;
        g.held = false;
        g.outcome = Some(Err(Error::Cancelled));
        s.progress.cancel();
        s.done.notify_all();
        true
    }

    /// Asks the job to stop. A queued job stops at once and never runs; a running one notices at its next check (a factorisation finishes first).
    /// False when no such job, or it has finished.
    pub fn cancel(&self, id: JobId) -> bool {
        let Some(s) = self.slot(id) else { return false };
        if self.cancel_queued(id) {
            return true;
        }
        if s.inner.lock().expect("job slot").state == JobState::Running {
            s.progress.cancel();
            return true;
        }
        false
    }

    /// Keeps a queued job from starting until `release`. False when it is not queued.
    pub fn hold(&self, id: JobId) -> bool {
        match self.slot(id) {
            Some(s) => {
                let mut g = s.inner.lock().expect("job slot");
                if g.state == JobState::Queued {
                    g.held = true;
                    true
                } else {
                    false
                }
            }
            None => false,
        }
    }

    pub fn release(&self, id: JobId) -> bool {
        let ok = match self.slot(id) {
            Some(s) => {
                let mut g = s.inner.lock().expect("job slot");
                let ok = g.state == JobState::Queued;
                if ok {
                    g.held = false;
                }
                ok
            }
            None => false,
        };
        pump(&self.shared);
        ok
    }

    /// Changes the priority of a queued job. False when it is not queued.
    pub fn set_priority(&self, id: JobId, priority: u8) -> bool {
        match self.slot(id) {
            Some(s) => {
                let mut g = s.inner.lock().expect("job slot");
                if g.state == JobState::Queued {
                    g.priority = priority.max(1);
                    true
                } else {
                    false
                }
            }
            None => false,
        }
    }

    /// Stops the queue from starting anything new; what runs finishes.
    pub fn pause(&self) {
        self.shared.sched.lock().expect("job queue").paused = true;
    }

    pub fn resume(&self) {
        self.shared.sched.lock().expect("job queue").paused = false;
        pump(&self.shared);
    }

    pub fn paused(&self) -> bool {
        self.shared.sched.lock().expect("job queue").paused
    }

    /// How many jobs a lane runs at once (heavy 1, light 2 by default).
    pub fn set_limit(&self, lane: Lane, n: usize) {
        {
            let mut s = self.shared.sched.lock().expect("job queue");
            match lane {
                Lane::Heavy => s.limit_heavy = n.max(1),
                Lane::Light => s.limit_light = n.max(1),
                Lane::Shared => s.limit_shared = n.max(1),
            }
        }
        pump(&self.shared);
    }

    /// After this long in the queue a job competes one priority level higher.
    pub fn set_age(&self, age: Duration) {
        self.shared.sched.lock().expect("job queue").age = age;
    }

    /// The stage log of the job so far.
    pub fn log(&self, id: JobId) -> Vec<String> {
        self.slot(id).map(|s| s.progress.log()).unwrap_or_default()
    }

    /// The outcome if the job has finished and it has not been taken; `None` while queued or running, or after it was taken.
    pub fn try_take(&self, id: JobId) -> Option<Result<JobOutput, Error>> {
        self.slot(id)?.inner.lock().expect("job slot").outcome.take()
    }

    /// Why a queued job cannot start however long it is waited for: it is held, or the queue is paused. `None` for a job that is running or will
    /// start (behind a running one), and for a finished one. A caller that would block with no timeout asks first.
    pub fn blocked(&self, id: JobId) -> Option<&'static str> {
        let s = self.slot(id)?;
        let g = s.inner.lock().expect("job slot");
        if g.state != JobState::Queued {
            return None;
        }
        if g.held {
            Some("it is held (release it with queue_release, or wait with a timeout_s)")
        } else if self.paused() {
            Some("the queue is paused (resume it with queue_resume, or wait with a timeout_s)")
        } else {
            None
        }
    }

    /// Blocks until the job has finished and takes its outcome.
    pub fn wait(&self, id: JobId) -> Result<JobOutput, Error> {
        let s = self.slot(id).ok_or_else(|| Error::invalid(format!("no job {id}")))?;
        let mut g = s.inner.lock().expect("job slot");
        while g.state.is_open() {
            g = s.done.wait(g).expect("job slot");
        }
        g.outcome.take().unwrap_or_else(|| Err(Error::invalid(format!("the outcome of job {id} was already taken"))))
    }

    /// Like `wait`, but gives up after `timeout`: `None` when the job is still queued or running then (it keeps going and can be waited for again).
    pub fn wait_timeout(&self, id: JobId, timeout: Duration) -> Option<Result<JobOutput, Error>> {
        let s = match self.slot(id) {
            Some(s) => s,
            None => return Some(Err(Error::invalid(format!("no job {id}")))),
        };
        let g = s.inner.lock().expect("job slot");
        let (mut g, _) = s.done.wait_timeout_while(g, timeout, |g| g.state.is_open()).expect("job slot");
        if g.state.is_open() {
            return None;
        }
        Some(
            g.outcome
                .take()
                .unwrap_or_else(|| Err(Error::invalid(format!("the outcome of job {id} was already taken")))),
        )
    }
}
