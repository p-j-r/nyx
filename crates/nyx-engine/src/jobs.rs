//! The job runner: a mesh build or a solve runs on its own thread with a `Progress` handle, so it can be watched (stage, fraction, elapsed)
//! and cancelled, by the window, a script or the MCP server alike. The runner knows nothing of windows: a front-end asks `info`, `cancel`,
//! `try_take` (non-blocking, for a repaint loop) or `wait` (blocking, for a tool call).

use crate::meshing::MeshedModel;
use femcore::progress::Progress;
use nyx_physics::{Error, SolveResult};
use serde::Serialize;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::Instant;

pub type JobId = u64;

/// What a finished job hands back. The engine applies it to its state (`Engine::apply`).
pub enum JobOutput {
    Mesh(Arc<MeshedModel>),
    /// the mesh the solve used (built inside the job when none was current) and one result per physics, in solve order
    Solved {
        model: Arc<MeshedModel>,
        results: Vec<SolveResult>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case", tag = "state", content = "message")]
pub enum JobState {
    Running,
    Done,
    Cancelled,
    Failed(String),
}

#[derive(Clone, Debug, Serialize)]
pub struct JobInfo {
    pub id: JobId,
    pub label: String,
    #[serde(flatten)]
    pub state: JobState,
    /// the stage the job is in now (plain words), empty before the first
    pub stage: String,
    pub fraction: Option<f64>,
    pub seconds: f64,
}

struct Inner {
    state: JobState,
    outcome: Option<Result<JobOutput, Error>>,
    seconds: f64,
}

struct Slot {
    label: String,
    progress: Arc<Progress>,
    started: Instant,
    inner: Mutex<Inner>,
    done: Condvar,
}

#[derive(Default)]
pub struct JobRunner {
    next: AtomicU64,
    slots: Mutex<BTreeMap<JobId, Arc<Slot>>>,
}

impl JobRunner {
    pub fn new() -> JobRunner {
        JobRunner::default()
    }

    /// Starts `work` on a new thread. The closure gets the job's `Progress`; a panic in it becomes a failed job, not a dead process.
    pub fn start<F>(&self, label: &str, work: F) -> JobId
    where
        F: FnOnce(Arc<Progress>) -> Result<JobOutput, Error> + Send + 'static,
    {
        let id = self.next.fetch_add(1, Ordering::SeqCst) + 1;
        let slot = Arc::new(Slot {
            label: label.into(),
            progress: Progress::new(),
            started: Instant::now(),
            inner: Mutex::new(Inner { state: JobState::Running, outcome: None, seconds: 0.0 }),
            done: Condvar::new(),
        });
        self.slots.lock().expect("job table").insert(id, slot.clone());
        let progress = slot.progress.clone();
        std::thread::spawn(move || {
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(progress))).unwrap_or_else(|p| {
                let text =
                    p.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| p.downcast_ref::<String>().cloned());
                Err(Error::Solver(format!("the job stopped unexpectedly: {}", text.unwrap_or_else(|| "panic".into()))))
            });
            let mut g = slot.inner.lock().expect("job slot");
            g.seconds = slot.started.elapsed().as_secs_f64();
            g.state = match &r {
                Ok(_) => JobState::Done,
                Err(Error::Cancelled) => JobState::Cancelled,
                Err(e) => JobState::Failed(e.to_string()),
            };
            g.outcome = Some(r);
            slot.done.notify_all();
        });
        id
    }

    fn slot(&self, id: JobId) -> Option<Arc<Slot>> {
        self.slots.lock().expect("job table").get(&id).cloned()
    }

    pub fn info(&self, id: JobId) -> Option<JobInfo> {
        let s = self.slot(id)?;
        let g = s.inner.lock().expect("job slot");
        let st = s.progress.status();
        Some(JobInfo {
            id,
            label: s.label.clone(),
            state: g.state.clone(),
            stage: st.stage,
            fraction: st.fraction,
            seconds: if g.state == JobState::Running { s.started.elapsed().as_secs_f64() } else { g.seconds },
        })
    }

    /// The progress handle of a job: the stage it is in, the log so far, the cancel flag. A front-end that shows a bar reads it.
    pub fn progress(&self, id: JobId) -> Option<Arc<Progress>> {
        self.slot(id).map(|s| s.progress.clone())
    }

    /// All jobs, oldest first.
    pub fn list(&self) -> Vec<JobInfo> {
        let ids: Vec<JobId> = self.slots.lock().expect("job table").keys().copied().collect();
        ids.into_iter().filter_map(|i| self.info(i)).collect()
    }

    /// The jobs that have not finished.
    pub fn running(&self) -> Vec<JobId> {
        self.list().into_iter().filter(|j| j.state == JobState::Running).map(|j| j.id).collect()
    }

    /// Asks the job to stop; the solver notices at its next check. False when no such running job.
    pub fn cancel(&self, id: JobId) -> bool {
        match self.slot(id) {
            Some(s) if s.inner.lock().expect("job slot").state == JobState::Running => {
                s.progress.cancel();
                true
            }
            _ => false,
        }
    }

    /// The stage log of the job so far.
    pub fn log(&self, id: JobId) -> Vec<String> {
        self.slot(id).map(|s| s.progress.log()).unwrap_or_default()
    }

    /// The outcome if the job has finished and it has not been taken; `None` while running or after it was taken.
    pub fn try_take(&self, id: JobId) -> Option<Result<JobOutput, Error>> {
        self.slot(id)?.inner.lock().expect("job slot").outcome.take()
    }

    /// Blocks until the job has finished and takes its outcome.
    pub fn wait(&self, id: JobId) -> Result<JobOutput, Error> {
        let s = self.slot(id).ok_or_else(|| Error::invalid(format!("no job {id}")))?;
        let mut g = s.inner.lock().expect("job slot");
        while g.state == JobState::Running {
            g = s.done.wait(g).expect("job slot");
        }
        g.outcome.take().unwrap_or_else(|| Err(Error::invalid(format!("the outcome of job {id} was already taken"))))
    }
}
