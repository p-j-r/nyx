//! The queue of the engine (productivity plan Q1, ADR 87): what a front-end sees of the job runner's lanes (`queue_*` commands), the options a
//! request is submitted with, and the stage caches that let queued requests share work. A request that comes after another of the same project
//! finds the mesh the first built; a convergence check right after a solve finds its coarse result done. Stages are looked up by content
//! (the mesh key, the problem key and the chain that was solved: the problem key alone is not enough, because a physics solved alone and one
//! solved with the physics it imports from share a key and not their numbers, ADR 84), so the queue's order can never change a number.

use crate::engine::Engine;
use crate::jobs::{JobId, Lane, Submit};
use crate::meshing::MeshedModel;
use nyx_physics::{Error, SolveResult};
use serde_json::{json, Value};
use std::sync::{Arc, Mutex};
use std::time::Duration;

/// How a request is queued. The default is the engine's owner at its usual priority, with nothing superseded.
#[derive(Clone, Debug, Default)]
pub struct QueueOpts {
    /// 1 a person's request now, 2 an agent's single request, 3 sweeps and batches; default by owner (`gui` 1, others 2)
    pub priority: Option<u8>,
    /// replace a request of the same owner and kind and targets that has not started yet
    pub supersede: bool,
    /// start held: it runs only after `queue_release`
    pub hold: bool,
    /// run a solve in a worker process (a crash or running out of memory fails the job, not this process; a cancel stops it at once); `None`
    /// is what the engine was told to do by default (`Engine::set_worker`)
    pub isolate: Option<bool>,
    /// start a solve even when its memory estimate says the machine cannot hold it (admission, `estimate.rs`)
    pub force: bool,
    /// hand a solve to the queue of the workspace (Q3: the queue worker runs it, whichever program asked; it implies `isolate`); `None` is what the
    /// engine was told to do by default (`Engine::set_shared_queue`)
    pub shared: Option<bool>,
}

/// What the engine remembers of a job beyond the runner's own record.
pub(crate) struct JobMeta {
    pub kind: &'static str,
    pub targets: Vec<String>,
    /// the hash of the whole project when the request was made: the request carries a frozen copy, so it says when the project has moved on
    pub stamp: String,
    /// what the job is expected to cost, filled in by the job once its mesh is known (`estimate.rs`); empty while there is no estimate
    pub estimate: Arc<Mutex<Option<crate::estimate::Estimate>>>,
    /// the id of the item in the workspace queue, once the job has written it (a solve that was handed to the queue worker)
    pub shared_id: Option<Arc<Mutex<Option<String>>>>,
}

/// What a mesh convergence check needs of the coarse solve: the numbers that identify its solutions, the unknowns and the cells. A result's
/// fields are gigabytes at the sizes that matter, so the stage cache keeps this, not the result.
pub(crate) struct SolveSummary {
    pub dofs: usize,
    pub cells: usize,
    pub numbers: Vec<(String, f64)>,
}

impl SolveSummary {
    pub fn of(r: &SolveResult, cells: usize) -> SolveSummary {
        SolveSummary { dofs: r.dofs, cells, numbers: crate::convergence::numbers(r) }
    }
}

/// Recent stages, found by content. Small on purpose: a mesh of a large model is gigabytes.
#[derive(Default)]
pub(crate) struct Stages {
    meshes: Mutex<Vec<(String, Arc<MeshedModel>)>>,
    solves: Mutex<Vec<(String, Arc<SolveSummary>)>>,
}

const KEEP_MESHES: usize = 2;
const KEEP_SOLVES: usize = 8;

impl Stages {
    pub fn mesh(&self, key: &str) -> Option<Arc<MeshedModel>> {
        self.meshes.lock().expect("stage cache").iter().find(|(k, _)| k == key).map(|(_, m)| m.clone())
    }

    pub fn put_mesh(&self, key: &str, m: &Arc<MeshedModel>) {
        let mut g = self.meshes.lock().expect("stage cache");
        g.retain(|(k, _)| k != key);
        g.push((key.to_string(), m.clone()));
        let n = g.len();
        if n > KEEP_MESHES {
            g.drain(..n - KEEP_MESHES);
        }
    }

    pub fn solve(&self, key: &str) -> Option<Arc<SolveSummary>> {
        self.solves.lock().expect("stage cache").iter().find(|(k, _)| k == key).map(|(_, r)| r.clone())
    }

    pub fn put_solve(&self, key: &str, r: &Arc<SolveSummary>) {
        let mut g = self.solves.lock().expect("stage cache");
        g.retain(|(k, _)| k != key);
        g.push((key.to_string(), r.clone()));
        let n = g.len();
        if n > KEEP_SOLVES {
            g.drain(..n - KEEP_SOLVES);
        }
    }

    pub fn clear(&self) {
        self.meshes.lock().expect("stage cache").clear();
        self.solves.lock().expect("stage cache").clear();
    }
}

/// The cache key of a solve: the problem key of the physics and the chain it was solved in (the sorted targets of its job).
pub(crate) fn solve_key(problem_key: &str, targets: &[String]) -> String {
    let mut t: Vec<&str> = targets.iter().map(String::as_str).collect();
    t.sort_unstable();
    t.dedup();
    format!("{problem_key}|{}", t.join(","))
}

fn shared_id(args: &Value) -> Option<String> {
    args.get("shared_id").and_then(Value::as_str).map(String::from)
}

fn priority_of(args: &Value) -> Result<Option<u8>, Error> {
    match args.get("priority") {
        None | Some(Value::Null) => Ok(None),
        Some(v) => match v.as_u64() {
            Some(p @ 1..=3) => Ok(Some(p as u8)),
            _ => Err(Error::invalid("`priority` is 1 (a person's request now), 2 (an agent's request, the default) or 3 (sweeps and batches)")),
        },
    }
}

impl Engine {
    /// Names who this engine's requests are made for (`gui`, `ai`, `cli`): the queue serves owners in turn, and a person's requests run first.
    pub fn set_owner(&mut self, owner: &str) {
        self.owner = owner.to_string();
    }

    pub fn owner(&self) -> &str {
        &self.owner
    }

    /// The queue of the workspace, when the engine has a workspace.
    pub(crate) fn shared(&self) -> Option<crate::shared_queue::SharedQueue> {
        self.shared_root.as_deref().or(self.workspace.as_deref()).map(crate::shared_queue::SharedQueue::new)
    }

    /// The item id of a job of this engine that was handed to the workspace queue.
    pub(crate) fn shared_of_job(&self, id: JobId) -> Option<String> {
        self.meta.get(&id)?.shared_id.as_ref()?.lock().expect("shared id").clone()
    }

    /// Runs `f` on the item named by `shared_id` in the arguments; false when there is no such item.
    fn shared_op(
        &self,
        args: &Value,
        f: impl Fn(&crate::shared_queue::SharedQueue, &str) -> bool,
    ) -> Result<bool, Error> {
        let sid = shared_id(args).unwrap_or_default();
        let q =
            self.shared().ok_or_else(|| Error::invalid("this engine has no workspace, so there is no shared queue"))?;
        Ok(f(&q, &sid))
    }

    /// The submission of a job of this engine: owner, priority by owner unless given, and the key a newer request of the same kind replaces.
    pub(crate) fn submission(&self, label: &str, lane: Lane, opts: &QueueOpts, key: &str) -> Submit {
        let mut s = Submit::new(label, lane);
        s.owner = self.owner.clone();
        s.priority = opts.priority.unwrap_or(if self.owner == "gui" { 1 } else { 2 });
        s.hold = opts.hold;
        s.key = Some(key.to_string());
        s.supersede = opts.supersede;
        s
    }

    /// A job as the queue commands show it: the runner's record plus what was asked and whether the project has changed since.
    pub fn queue_entry(&self, id: JobId) -> Option<Value> {
        let info = self.jobs.info(id)?;
        let mut v = serde_json::to_value(&info).ok()?;
        if let Some(m) = self.meta.get(&id) {
            v["kind"] = json!(m.kind);
            v["targets"] = json!(m.targets);
            v["project_changed"] = json!(m.stamp != self.stamp());
            if let Some(e) = m.estimate.lock().expect("estimate").as_ref() {
                v["estimate"] = json!(e);
            }
        }
        if let Some(sid) =
            self.meta.get(&id).and_then(|m| m.shared_id.as_ref()).and_then(|s| s.lock().expect("shared id").clone())
        {
            v["shared"] = json!(true);
            v["shared_id"] = json!(sid);
            if let Some(item) = self.shared().and_then(|q| q.item(&sid)) {
                v["place"] = json!(item.position);
                v["queue_state"] = json!(item.state_name());
                v["held"] = json!(item.held);
            }
        }
        v["applied"] = json!(!self.pending.contains_key(&id));
        Some(v)
    }

    /// The queue, oldest first, optionally of one owner or one state (`queued`, `running`, `done`, `cancelled`, `failed`).
    pub fn queue_list(&self, owner: Option<&str>, state: Option<&str>) -> Vec<Value> {
        let mut out = self.queue_list_own(owner, state);
        out.extend(self.shared_foreign(owner, state));
        out
    }

    /// The items of the workspace queue that another program (or an earlier session of this one) wrote: what makes the queue one queue.
    fn shared_foreign(&self, owner: Option<&str>, state: Option<&str>) -> Vec<Value> {
        let Some(q) = self.shared() else { return vec![] };
        let own = self.shared_own.lock().expect("own items").clone();
        q.items()
            .into_iter()
            .filter(|i| !own.contains(&i.request.id))
            .filter(|i| owner.is_none_or(|o| i.request.owner == o))
            .filter(|i| state.is_none_or(|s| i.state_name() == s))
            .map(|i| {
                let mut v = i.to_json(&q);
                v["foreign"] = json!(true);
                v
            })
            .collect()
    }

    fn queue_list_own(&self, owner: Option<&str>, state: Option<&str>) -> Vec<Value> {
        self.jobs
            .list()
            .into_iter()
            .filter(|i| owner.is_none_or(|o| i.owner == o))
            .filter(|i| {
                state.is_none_or(|s| {
                    serde_json::to_value(&i.state)
                        .ok()
                        .and_then(|v| v["state"].as_str().map(|x| x == s))
                        .unwrap_or(false)
                })
            })
            .filter_map(|i| self.queue_entry(i.id))
            .collect()
    }

    /// Waits for a job up to `timeout` and applies its outcome; at the timeout the job's entry comes back with `waiting: true` and the job goes on.
    pub fn wait_timeout(&mut self, id: JobId, timeout: Duration) -> Result<Value, Error> {
        if !self.pending.contains_key(&id) {
            return Err(Error::invalid(format!(
                "job {id} is not waiting to be applied (it was applied already, or there is no such job)"
            )));
        }
        match self.jobs.wait_timeout(id, timeout) {
            Some(outcome) => self.apply(id, outcome),
            None => {
                let mut v = self.queue_entry(id).unwrap_or_else(|| json!({ "id": id }));
                v["waiting"] = json!(true);
                Ok(v)
            }
        }
    }

    /// The `queue_*` commands (see the catalog).
    pub(crate) fn queue_command(&mut self, name: &str, args: &Value) -> Result<Value, Error> {
        let id = || {
            args.get("id")
                .and_then(Value::as_u64)
                .ok_or_else(|| Error::invalid("`id` is required and must be a whole number"))
        };
        match name {
            "queue_submit" => {
                let opts = QueueOpts {
                    priority: priority_of(args)?,
                    supersede: args.get("supersede").and_then(Value::as_bool).unwrap_or(false),
                    hold: args.get("hold").and_then(Value::as_bool).unwrap_or(false),
                    isolate: args.get("isolate").and_then(Value::as_bool),
                    force: args.get("force").and_then(Value::as_bool).unwrap_or(false),
                    shared: args.get("shared").and_then(Value::as_bool),
                };
                let physics = args.get("physics").and_then(Value::as_str);
                let job = match args.get("kind").and_then(Value::as_str).unwrap_or("solve") {
                    "solve" => {
                        let targets = match args.get("physics") {
                            None | Some(Value::Null) => None,
                            Some(Value::String(s)) => Some(vec![s.clone()]),
                            Some(Value::Array(a)) => {
                                Some(a.iter().filter_map(|x| x.as_str().map(String::from)).collect())
                            }
                            Some(_) => return Err(Error::invalid("`physics` is a code or a list of codes")),
                        };
                        if args.get("detached").and_then(Value::as_bool) == Some(true) {
                            return self.queue_detached(targets, &opts);
                        }
                        self.start_solve_with(targets, &opts)?
                    }
                    "mesh" => self.start_mesh_with(&opts)?,
                    "convergence" => {
                        let physics = physics.ok_or_else(|| Error::invalid("a convergence check needs `physics`"))?;
                        let factor = args.get("factor").and_then(Value::as_f64).unwrap_or(0.6);
                        self.start_convergence_with(
                            physics,
                            factor,
                            args.get("tolerance").and_then(Value::as_f64),
                            &opts,
                        )?
                    }
                    other => {
                        return Err(Error::invalid(format!("`kind` is solve, mesh or convergence, not '{other}'")))
                    }
                };
                Ok(self.queue_entry(job).unwrap_or_else(|| json!({ "id": job })))
            }
            "queue_list" => {
                Ok(json!(self
                    .queue_list(args.get("owner").and_then(Value::as_str), args.get("state").and_then(Value::as_str))))
            }
            "queue_get" if shared_id(args).is_some() => {
                let sid = shared_id(args).unwrap_or_default();
                let q = self
                    .shared()
                    .ok_or_else(|| Error::invalid("this engine has no workspace, so there is no shared queue"))?;
                let item = q
                    .item(&sid)
                    .ok_or_else(|| Error::invalid(format!("no item {sid} in the queue of the workspace")))?;
                let mut v = item.to_json(&q);
                let log = item.state.as_ref().map(|s| s.log.clone()).unwrap_or_default();
                v["log"] = json!(log[log.len().saturating_sub(20)..]);
                Ok(v)
            }
            "queue_cancel" if shared_id(args).is_some() => {
                Ok(json!({ "cancelled": self.shared_op(args, |q, s| q.cancel(s))? }))
            }
            "queue_hold" if shared_id(args).is_some() => Ok(json!({ "held": self.shared_op(args, |q, s| q.hold(s))? })),
            "queue_release" if shared_id(args).is_some() => {
                Ok(json!({ "released": self.shared_op(args, |q, s| q.release(s))? }))
            }
            "queue_set_priority" if shared_id(args).is_some() => {
                let p = priority_of(args)?.ok_or_else(|| Error::invalid("`priority` is required"))?;
                Ok(json!({ "changed": self.shared_op(args, |q, s| q.set_priority(s, p))? }))
            }
            "queue_hold" | "queue_release" | "queue_set_priority"
                if args.get("id").and_then(Value::as_u64).and_then(|i| self.shared_of_job(i)).is_some() =>
            {
                // a job of this engine that is waiting in the workspace queue: its item carries the marker
                let sid = self.shared_of_job(id()?).unwrap_or_default();
                let q = self.shared().ok_or_else(|| Error::invalid("this engine has no workspace"))?;
                Ok(match name {
                    "queue_hold" => json!({ "held": q.hold(&sid) }),
                    "queue_release" => json!({ "released": q.release(&sid) }),
                    _ => {
                        let p = priority_of(args)?.ok_or_else(|| Error::invalid("`priority` is required"))?;
                        json!({ "changed": q.set_priority(&sid, p) })
                    }
                })
            }
            "queue_get" => {
                let id = id()?;
                let mut v = self.queue_entry(id).ok_or_else(|| Error::invalid(format!("no job {id}")))?;
                let log = self.jobs.log(id);
                let from = log.len().saturating_sub(20);
                v["log"] = json!(log[from..]);
                Ok(v)
            }
            "queue_cancel" => Ok(json!({ "cancelled": self.jobs.cancel(id()?) })),
            "queue_hold" => Ok(json!({ "held": self.jobs.hold(id()?) })),
            "queue_release" => Ok(json!({ "released": self.jobs.release(id()?) })),
            "queue_set_priority" => {
                let p = priority_of(args)?.ok_or_else(|| Error::invalid("`priority` is required"))?;
                Ok(json!({ "changed": self.jobs.set_priority(id()?, p) }))
            }
            "queue_pause" => {
                self.jobs.pause();
                if let Some(q) = self.shared() {
                    q.pause();
                }
                Ok(json!({ "paused": true }))
            }
            "queue_resume" => {
                self.jobs.resume();
                if let Some(q) = self.shared() {
                    q.resume();
                }
                Ok(json!({ "paused": false }))
            }
            "queue_wait" => {
                let t =
                    args.get("timeout_s").and_then(Value::as_f64).filter(|t| *t >= 0.0 && t.is_finite()).ok_or_else(
                        || {
                            Error::invalid(
                            "`timeout_s` is required: seconds to wait at most (a tool call must not block for ever)",
                        )
                        },
                    )?;
                self.wait_timeout(id()?, Duration::from_secs_f64(t))
            }
            other => Err(Error::invalid(format!("no command '{other}'"))),
        }
    }

    /// Remembers that the job `id` was handed to the workspace queue, and under which item.
    pub(crate) fn note_shared(&mut self, id: JobId, item: Option<Arc<Mutex<Option<String>>>>) {
        if let Some(m) = self.meta.get_mut(&id) {
            m.shared_id = item;
        }
    }

    pub(crate) fn note_job(&mut self, id: JobId, kind: &'static str, targets: Vec<String>) {
        self.note_job_estimate(id, kind, targets, Default::default());
    }

    /// `note_job` for a job that fills in its own estimate (a solve, once its mesh is known).
    pub(crate) fn note_job_estimate(
        &mut self,
        id: JobId,
        kind: &'static str,
        targets: Vec<String>,
        estimate: Arc<Mutex<Option<crate::estimate::Estimate>>>,
    ) {
        let stamp = self.stamp();
        self.meta.insert(id, JobMeta { kind, targets, stamp, estimate, shared_id: None });
    }
}
