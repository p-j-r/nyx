//! Sweeps (productivity plan Q6, ADR 87): one run per point, recorded. A point is a list of edits (any command of access `edit`: a geometry
//! parameter, a study value such as the wavelength, a material, a mesh setting) applied to a copy of the project; the case decides what varies,
//! the helper does not. Every point is a job of the queue at priority 3 (so a person's request and an agent's single request go before a
//! sweep and a sweep never starves them), carries its own frozen project, is solved on the same stage caches as any other request (points
//! that share a mesh share it) and is recorded in the run table when the project is bound to a workflow. A sweep never touches the engine's own
//! project, mesh or results: the table is what it returns. The edits of every point are applied and checked before anything is submitted, so
//! a bad point fails the whole request with its number and nothing runs.

use crate::engine::{Engine, Purpose};
use crate::jobs::{JobId, JobOutput};
use crate::queue::QueueOpts;
use crate::store_link::RunOutcome;
use nyx_physics::Error;
use serde_json::{json, Map, Value};
use std::time::{Duration, Instant};

/// The most points one request may make: a sweep is a queue of jobs, not a way to start ten thousand.
pub const MAX_POINTS: usize = 500;

pub(crate) struct Point {
    pub label: String,
    pub edits: Vec<Value>,
    pub job: JobId,
    /// what the point came to, once its job was applied
    pub row: Option<Value>,
}

pub(crate) struct Sweep {
    pub name: String,
    pub targets: Vec<String>,
    pub points: Vec<Point>,
}

/// A checked sweep request: the points, the physics named (if any), the queue options and one engine per point.
pub(crate) struct Prepared {
    asked: Vec<Asked>,
    pub targets: Option<Vec<String>>,
    opts: QueueOpts,
    pub engines: Vec<Engine>,
}

/// A point as asked for: its label and its edits (`{ "cmd": ..., "args": {...} }`).
struct Asked {
    label: String,
    edits: Vec<Value>,
}

fn asked_points(args: &Value) -> Result<Vec<Asked>, Error> {
    let (points, series) = (args.get("points").filter(|v| !v.is_null()), args.get("series").filter(|v| !v.is_null()));
    match (points, series) {
        (Some(_), Some(_)) => Err(Error::invalid("give `points` or `series`, not both")),
        (None, None) => Err(Error::invalid(
            "give `points` (a list of { label, edits: [{ cmd, args }] }) or `series` ({ cmd, args, field, values })",
        )),
        (Some(p), None) => {
            let list = p.as_array().ok_or_else(|| Error::invalid("`points` is a list"))?;
            list.iter()
                .enumerate()
                .map(|(i, p)| {
                    let edits = p
                        .get("edits")
                        .and_then(Value::as_array)
                        .ok_or_else(|| Error::invalid(format!("point {i}: `edits` is a list of {{ cmd, args }}")))?
                        .clone();
                    let label = p
                        .get("label")
                        .and_then(Value::as_str)
                        .map(String::from)
                        .unwrap_or_else(|| format!("point {i}"));
                    Ok(Asked { label, edits })
                })
                .collect()
        }
        (None, Some(s)) => {
            let cmd = s.get("cmd").and_then(Value::as_str).ok_or_else(|| {
                Error::invalid("`series.cmd` is the edit command to repeat (set_param, set_study ...)")
            })?;
            let field = s
                .get("field")
                .and_then(Value::as_str)
                .ok_or_else(|| Error::invalid("`series.field` is the argument of the command that takes each value"))?;
            let values = s
                .get("values")
                .and_then(Value::as_array)
                .ok_or_else(|| Error::invalid("`series.values` is a list of the values the field takes"))?;
            let base: Map<String, Value> = match s.get("args") {
                None | Some(Value::Null) => Map::new(),
                Some(Value::Object(o)) => o.clone(),
                Some(_) => return Err(Error::invalid("`series.args` is an object")),
            };
            Ok(values
                .iter()
                .map(|v| {
                    let mut a = base.clone();
                    a.insert(field.to_string(), v.clone());
                    Asked { label: format!("{field}={v}"), edits: vec![json!({ "cmd": cmd, "args": a })] }
                })
                .collect())
        }
    }
}

impl Engine {
    /// A copy of this engine for editing: the project is a copy, the mesher, the registry and the material library are shared. Nothing it does
    /// reaches this engine.
    pub(crate) fn scratch(&self) -> Engine {
        let mut e = Engine::with(self.registry.clone(), self.mesher.clone());
        e.base_library = self.base_library.clone();
        e.project = self.project.clone();
        e.owner = self.owner.clone();
        e
    }

    /// Applies the edits of point `i` to a copy of the project and checks the result. Only commands that change the project are allowed (never
    /// one that replaces it, starts work or writes a file).
    fn point_engine(&self, i: usize, edits: &[Value], allowed: &[String]) -> Result<Engine, Error> {
        let mut e = self.scratch();
        for (j, edit) in edits.iter().enumerate() {
            let at = |what: String| Error::invalid(format!("point {i}, edit {j}: {what}"));
            let cmd = edit
                .get("cmd")
                .and_then(Value::as_str)
                .ok_or_else(|| at("an edit is { \"cmd\": ..., \"args\": {...} }".into()))?;
            if !allowed.iter().any(|a| a == cmd) {
                return Err(at(format!(
                    "'{cmd}' is not an edit a point may make (the commands that change the project, except those that replace it)"
                )));
            }
            let args = edit.get("args").cloned().unwrap_or_else(|| json!({}));
            e.call(cmd, &args).map_err(|err| at(format!("{cmd}: {err}")))?;
        }
        e.project.check(&e.registry).map_err(|err| Error::invalid(format!("point {i}: {err}")))?;
        Ok(e)
    }

    /// The edit commands a point may use.
    fn point_commands(&self) -> Vec<String> {
        self.catalog()
            .into_iter()
            .filter(|c| c.access == crate::commands::Access::Edit)
            .map(|c| c.name)
            .filter(|n| {
                !matches!(n.as_str(), "project_new" | "project_load" | "project_open" | "store_open" | "open_results")
            })
            .collect()
    }

    pub(crate) fn sweep_command(&mut self, name: &str, args: &Value) -> Result<Value, Error> {
        match name {
            "sweep_run" => self.sweep_run(args),
            "sweep_table" => self.sweep_table(args),
            "sweep_cancel" => {
                let id = sweep_id(args)?;
                let sweep = self.sweeps.get(&id).ok_or_else(|| Error::invalid(format!("no sweep {id}")))?;
                let n = sweep.points.iter().filter(|p| p.row.is_none()).filter(|p| self.jobs.cancel(p.job)).count();
                Ok(json!({ "sweep": id, "cancelled": n }))
            }
            other => Err(Error::invalid(format!("no command '{other}'"))),
        }
    }

    /// Reads and checks a sweep request without submitting anything: the points, what they solve, how they queue, and the engine of every point
    /// (a copy of this one with the point's edits applied). `sweep_run` submits these; a plan of a command list counts them.
    pub(crate) fn sweep_prepare(&self, args: &Value) -> Result<Prepared, Error> {
        let asked = asked_points(args)?;
        if asked.is_empty() {
            return Err(Error::invalid("a sweep needs at least one point"));
        }
        if asked.len() > MAX_POINTS {
            return Err(Error::invalid(format!("a sweep has at most {MAX_POINTS} points; split it")));
        }
        let targets = match args.get("physics") {
            None | Some(Value::Null) => None,
            Some(Value::String(s)) => Some(vec![s.clone()]),
            Some(Value::Array(a)) => Some(a.iter().filter_map(|x| x.as_str().map(String::from)).collect()),
            Some(_) => return Err(Error::invalid("`physics` is a code or a list of codes")),
        };
        let priority = match args.get("priority") {
            None | Some(Value::Null) => 3,
            Some(v) => match v.as_u64() {
                Some(p @ 1..=3) => p as u8,
                _ => return Err(Error::invalid("`priority` is 1, 2 or 3 (3, the default, is for sweeps)")),
            },
        };
        let opts = QueueOpts {
            priority: Some(priority),
            supersede: false,
            hold: args.get("hold").and_then(Value::as_bool).unwrap_or(false),
            isolate: args.get("isolate").and_then(Value::as_bool),
            force: args.get("force").and_then(Value::as_bool).unwrap_or(false),
        };
        // every point is edited and checked first: a bad one fails the request and nothing is submitted
        let allowed = self.point_commands();
        let mut engines = Vec::with_capacity(asked.len());
        for (i, a) in asked.iter().enumerate() {
            let e = self.point_engine(i, &a.edits, &allowed)?;
            e.library().map_err(|err| Error::invalid(format!("point {i}: {err}")))?;
            for t in targets.as_deref().unwrap_or(&[]) {
                if e.project.spec(t).is_none() {
                    return Err(Error::invalid(format!("point {i}: the project has no problem for '{t}'")));
                }
            }
            engines.push(e);
        }
        Ok(Prepared { asked, targets, opts, engines })
    }

    fn sweep_run(&mut self, args: &Value) -> Result<Value, Error> {
        let Prepared { asked, targets, opts, engines } = self.sweep_prepare(args)?;
        // a host that records (the MCP server) gets its workflow before the first run, as a solve does
        if self.auto_store && self.link.is_none() && self.workspace.is_some() {
            self.save_to_workspace(None)?;
        }
        self.next_sweep += 1;
        let id = self.next_sweep;
        let mut points = vec![];
        let mut resolved = vec![];
        for (i, (a, e)) in asked.into_iter().zip(engines.iter()).enumerate() {
            match self.submit_solve(Some(e), targets.clone(), &opts, Purpose::Point { sweep: id, index: i }) {
                Ok(job) => {
                    resolved = targets.clone().unwrap_or_else(|| e.default_targets());
                    points.push(Point { label: a.label, edits: a.edits, job, row: None });
                }
                Err(err) => {
                    // nothing of a failed request stays: what was submitted is cancelled and forgotten
                    for p in &points {
                        self.jobs.cancel(p.job);
                        self.pending.remove(&p.job);
                    }
                    return Err(Error::invalid(format!("point {i}: {err}")));
                }
            }
        }
        let name = args.get("name").and_then(Value::as_str).map(String::from).unwrap_or_else(|| format!("sweep {id}"));
        let listed: Vec<Value> =
            points.iter().enumerate().map(|(i, p)| json!({ "index": i, "label": p.label, "job": p.job })).collect();
        self.sweeps.insert(id, Sweep { name: name.clone(), targets: resolved, points });
        Ok(json!({
            "sweep": id, "name": name, "points": listed,
            "note": "the points are queued at their priority and run one at a time; sweep_table shows them as they finish, and the current project and results are not touched",
        }))
    }

    fn sweep_table(&mut self, args: &Value) -> Result<Value, Error> {
        let Some(id) = args.get("sweep").filter(|v| !v.is_null()) else {
            return Ok(json!(self
                .sweeps
                .iter()
                .map(|(id, s)| {
                    let done = s.points.iter().filter(|p| p.row.is_some()).count();
                    json!({ "sweep": id, "name": s.name, "points": s.points.len(), "applied": done })
                })
                .collect::<Vec<_>>()));
        };
        let id = id.as_u64().ok_or_else(|| Error::invalid("`sweep` is a whole number"))?;
        let wait = args.get("wait_s").and_then(Value::as_f64).filter(|t| *t >= 0.0 && t.is_finite()).unwrap_or(0.0);
        let jobs: Vec<JobId> = self
            .sweeps
            .get(&id)
            .ok_or_else(|| Error::invalid(format!("no sweep {id}")))?
            .points
            .iter()
            .filter(|p| p.row.is_none())
            .map(|p| p.job)
            .collect();
        let deadline = Instant::now() + Duration::from_secs_f64(wait);
        for job in jobs {
            if !self.pending.contains_key(&job) {
                continue;
            }
            let outcome = if wait > 0.0 {
                self.jobs.wait_timeout(job, deadline.saturating_duration_since(Instant::now()))
            } else {
                self.jobs.try_take(job)
            };
            if let Some(o) = outcome {
                let _ = self.apply(job, o);
            }
        }
        let sweep = &self.sweeps[&id];
        let mut rows = vec![];
        for (i, p) in sweep.points.iter().enumerate() {
            rows.push(match &p.row {
                Some(r) => r.clone(),
                None => {
                    let live = self.queue_entry(p.job).unwrap_or_else(|| json!({ "job": p.job }));
                    json!({ "index": i, "label": p.label, "edits": p.edits, "job": p.job, "state": live["state"], "stage": live["stage"] })
                }
            });
        }
        let finished = sweep.points.iter().filter(|p| p.row.is_some()).count();
        Ok(json!({
            "sweep": id, "name": sweep.name, "physics": sweep.targets,
            "points": rows, "finished": finished, "total": sweep.points.len(), "complete": finished == sweep.points.len(),
        }))
    }

    /// Applies a finished point: records its run, keeps its numbers in the table and drops its result. The engine's own results and mesh are
    /// never replaced by a point.
    pub(crate) fn apply_point(
        &mut self,
        id: JobId,
        sweep: u64,
        index: usize,
        project: crate::project::Project,
        outcome: Result<JobOutput, Error>,
    ) -> Result<Value, Error> {
        let seconds = self.jobs.info(id).map(|i| i.seconds).unwrap_or(0.0);
        let mut log = self.jobs.log(id);
        let (state, error, run, numbers, dofs, cells) = match outcome {
            Ok(JobOutput::Solved { model, results, failed, worker }) => {
                let mesh = model.solver_mesh();
                for r in &results {
                    log.extend(r.log.iter().map(|l| format!("{}: {l}", r.physics)));
                }
                let mut numbers = Map::new();
                let mut dofs = Map::new();
                for r in &results {
                    let nums: Map<String, Value> =
                        crate::convergence::numbers(r).into_iter().map(|(k, v)| (k, json!(v))).collect();
                    numbers.insert(r.physics.clone(), Value::Object(nums));
                    dofs.insert(r.physics.clone(), json!(r.dofs));
                }
                let cells = model.mesh().n_cells();
                match failed {
                    Some(e) => {
                        log.push(format!("stopped: {e}"));
                        let run = self.record_point(&project, RunOutcome::Failed(&e), seconds, &log, Some(mesh));
                        ("failed", Some(e.to_string()), run, numbers, dofs, Some(cells))
                    }
                    None => {
                        self.learn(&results, cells as u64, model.mesh().dim() as u8, &worker);
                        let run = self.record_point(&project, RunOutcome::Solved(&results), seconds, &log, Some(mesh));
                        ("done", None, run, numbers, dofs, Some(cells))
                    }
                }
            }
            Ok(_) => (
                "failed",
                Some("a job returned something its kind does not produce".to_string()),
                Value::Null,
                Map::new(),
                Map::new(),
                None,
            ),
            Err(e) => {
                // a point cancelled before it started never ran: it leaves no record
                let never_ran = matches!(e, Error::Cancelled) && log.is_empty() && seconds == 0.0;
                let run = if never_ran {
                    Value::Null
                } else {
                    self.record_point(&project, RunOutcome::Failed(&e), seconds, &log, None)
                };
                let state = if matches!(e, Error::Cancelled) { "cancelled" } else { "failed" };
                (state, Some(e.to_string()), run, Map::new(), Map::new(), None)
            }
        };
        let (label, edits) = match self.sweeps.get(&sweep).and_then(|s| s.points.get(index)) {
            Some(p) => (p.label.clone(), p.edits.clone()),
            None => (format!("point {index}"), vec![]),
        };
        let row = json!({
            "index": index, "label": label, "edits": edits, "job": id, "state": state, "error": error,
            "numbers": numbers, "unknowns": dofs, "cells": cells, "seconds": seconds, "run": run,
        });
        if let Some(p) = self.sweeps.get_mut(&sweep).and_then(|s| s.points.get_mut(index)) {
            p.row = Some(row.clone());
        }
        Ok(row)
    }

    fn record_point(
        &mut self,
        project: &crate::project::Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<std::sync::Arc<femcore::mesh::SimplexMesh>>,
    ) -> Value {
        match self.link.as_mut() {
            None => Value::Null,
            Some(l) => match l.record_point(project, outcome, seconds, log, mesh.as_deref()) {
                Ok(r) => json!({ "id": r.id, "version": r.version, "status": r.status }),
                Err(e) => json!({ "error": e.to_string() }),
            },
        }
    }
}

fn sweep_id(args: &Value) -> Result<u64, Error> {
    args.get("sweep").and_then(Value::as_u64).ok_or_else(|| Error::invalid("`sweep` is required and a whole number"))
}
