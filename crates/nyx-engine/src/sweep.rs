// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Sweeps (productivity plan Q6, ADR 87; VVUQ Phase 4 LHS/Sobol QMC streaming §6.4):
//!
//! - **Explicit sweeps (`points` / `series`)**: validates each point one-at-a-time via a temporary scratch
//!   [`Engine`] (eliminating upfront `Vec<Engine>` cloning), then queues one priority-3 job per point on shared
//!   stage caches (`warm` reuse supported).
//! - **Deterministic QMC sweeps (`sampling`: `"lhs"` / `"sobol"`)**: generates deterministic `SplitMix64`
//!   Latin Hypercube Sampling or low-discrepancy Sobol coordinates over parametric geometry/material/study inputs
//!   and streams execution sample-by-sample inside a single priority-3 `Lane::Heavy` job in `O(1)` memory,
//!   incrementally checkpointing solved points and `skipped_invalid_geometry` diagnostics to HDF5 `SweepCheckpoint`
//!   without aborting the ensemble on invalid parameter combinations.

use crate::engine::{Engine, Pending, Purpose};
use crate::jobs::{JobId, JobOutput, Lane};
use crate::qmc::{generate_lhs_unit, generate_sobol_unit, MarginDistribution, QmcMethod, MAX_QMC_SAMPLES};
use crate::queue::QueueOpts;
use crate::store_link::RunOutcome;
use nyx_physics::Error;
use serde_json::{json, Map, Value};
use std::collections::BTreeSet;
use std::sync::Arc;
use std::time::{Duration, Instant};

/// The most points one explicit `points` / `series` request may make: a sweep is a queue of jobs, not a way to start ten thousand.
pub const MAX_POINTS: usize = 500;
const DEFAULT_QMC_SEED: u64 = 0x4E59_585F_514D_4304;

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
    pub qmc_job: Option<JobId>,
    pub qmc_summary: Option<Value>,
}

#[derive(Clone, Debug)]
pub(crate) struct QmcParamDim {
    pub col_name: String,
    pub cmd: String,
    pub base_args: Map<String, Value>,
    pub field: String,
    pub distribution: MarginDistribution,
}

#[derive(Clone, Debug)]
pub(crate) struct QmcPlan {
    pub method: QmcMethod,
    pub seed: u64,
    pub scramble: bool,
    pub offset: usize,
    pub skip_invalid_geometry: bool,
    pub dims: Vec<QmcParamDim>,
    pub samples: Vec<Vec<f64>>,
}

/// A checked sweep request: explicit points or streaming QMC plan, targets, queue options, and content mesh keys.
pub(crate) struct Prepared {
    asked: Vec<Asked>,
    pub targets: Option<Vec<String>>,
    opts: QueueOpts,
    /// `reuse: "warm"`: each point starts its eigensolves from the modes of the one before (P5b, ADR 96)
    warm: bool,
    pub mesh_keys: Vec<String>,
    checkpoint_path: Option<String>,
    checkpoint_lzf: bool,
    pub sampling: Option<QmcPlan>,
}

impl Prepared {
    pub fn point_count(&self) -> usize {
        match &self.sampling {
            Some(q) => q.samples.len(),
            None => self.asked.len(),
        }
    }

    pub fn is_sampling(&self) -> bool {
        self.sampling.is_some()
    }
}

/// A point as asked for: its label and its edits (`{ "cmd": ..., "args": {...} }`).
struct Asked {
    label: String,
    edits: Vec<Value>,
}

fn asked_points(args: &Value) -> Result<Vec<Asked>, Error> {
    let (points, series) = (args.get("points").filter(|v| !v.is_null()), args.get("series").filter(|v| !v.is_null()));
    match (points, series) {
        (Some(_), Some(_)) => Err(Error::invalid("give `points`, `series`, or `sampling`, not multiple")),
        (None, None) => Err(Error::invalid(
            "give `points` (a list of { label, edits: [{ cmd, args }] }), `series` ({ cmd, args, field, values }), or `sampling`",
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

fn parse_qmc_plan(s_val: &Value, allowed: &[String]) -> Result<QmcPlan, Error> {
    let obj = s_val.as_object().ok_or_else(|| Error::invalid("`sampling` must be a JSON object"))?;
    let method = match obj.get("method").filter(|v| !v.is_null()).and_then(Value::as_str) {
        None | Some("lhs") => QmcMethod::Lhs,
        Some("sobol") => QmcMethod::Sobol,
        Some(other) => return Err(Error::invalid(format!("sampling.method must be 'lhs' or 'sobol', not '{other}'"))),
    };
    let n_samples = obj
        .get("n_samples")
        .and_then(Value::as_u64)
        .filter(|&n| (1..=MAX_QMC_SAMPLES as u64).contains(&n))
        .ok_or_else(|| Error::invalid(format!("sampling.n_samples must be an integer in 1..={MAX_QMC_SAMPLES}")))?
        as usize;
    let seed = match obj.get("seed").filter(|v| !v.is_null()) {
        None => DEFAULT_QMC_SEED,
        Some(v) => v.as_u64().ok_or_else(|| Error::invalid("sampling.seed must be a non-negative integer"))?,
    };
    let scramble = obj.get("scramble").and_then(Value::as_bool).unwrap_or(false);
    let offset = obj
        .get("offset")
        .filter(|v| !v.is_null())
        .map(|v| v.as_u64().map(|n| n as usize).ok_or_else(|| Error::invalid("sampling.offset must be an integer")))
        .transpose()?
        .unwrap_or(1);
    let skip_invalid_geometry = obj.get("skip_invalid_geometry").and_then(Value::as_bool).unwrap_or(true);

    let p_arr = obj
        .get("parameters")
        .and_then(Value::as_array)
        .filter(|a| !a.is_empty())
        .ok_or_else(|| Error::invalid("sampling.parameters must be a non-empty array"))?;

    let mut used_cols = BTreeSet::new();
    let mut dims = Vec::with_capacity(p_arr.len());
    for (i, item) in p_arr.iter().enumerate() {
        let po = item
            .as_object()
            .ok_or_else(|| Error::invalid(format!("sampling.parameters[{i}] must be a JSON object")))?;
        let cmd = po.get("cmd").and_then(Value::as_str).unwrap_or("set_param").to_string();
        if !allowed.iter().any(|a| a == &cmd) {
            return Err(Error::invalid(format!(
                "sampling.parameters[{i}]: '{cmd}' is not an edit command a sweep point may execute"
            )));
        }
        let base_args: Map<String, Value> = match po.get("args") {
            None | Some(Value::Null) => Map::new(),
            Some(Value::Object(o)) => o.clone(),
            Some(_) => return Err(Error::invalid(format!("sampling.parameters[{i}].args must be an object"))),
        };
        let field = po.get("field").and_then(Value::as_str).unwrap_or("value").to_string();
        let raw_col = po
            .get("name")
            .and_then(Value::as_str)
            .map(String::from)
            .or_else(|| base_args.get("name").and_then(Value::as_str).map(String::from))
            .unwrap_or_else(|| format!("{field}_{i}"));
        let mut col_name = raw_col.clone();
        if !used_cols.insert(col_name.clone()) {
            col_name = format!("{raw_col}_{i}");
            used_cols.insert(col_name.clone());
        }

        let bounds_arr = po.get("bounds").and_then(Value::as_array);
        let min_opt = po
            .get("min")
            .and_then(Value::as_f64)
            .or_else(|| bounds_arr.and_then(|b| b.first().and_then(Value::as_f64)));
        let max_opt =
            po.get("max").and_then(Value::as_f64).or_else(|| bounds_arr.and_then(|b| b.get(1).and_then(Value::as_f64)));

        let dist_kind = po.get("distribution").and_then(Value::as_str).unwrap_or_else(|| {
            if po.contains_key("mean") || po.contains_key("sigma") || po.contains_key("rel_tol") {
                "normal"
            } else {
                "uniform"
            }
        });

        let distribution = match dist_kind {
            "uniform" => {
                let min = min_opt.ok_or_else(|| {
                    Error::invalid(format!(
                        "sampling.parameters[{i}] uniform distribution requires 'min'/'max' or 'bounds'"
                    ))
                })?;
                let max = max_opt.ok_or_else(|| {
                    Error::invalid(format!(
                        "sampling.parameters[{i}] uniform distribution requires 'min'/'max' or 'bounds'"
                    ))
                })?;
                if !min.is_finite() || !max.is_finite() || min > max {
                    return Err(Error::invalid(format!(
                        "sampling.parameters[{i}] uniform bounds must be finite with min <= max, got [{min}, {max}]"
                    )));
                }
                MarginDistribution::Uniform { min, max }
            }
            "normal" => {
                let mean = po.get("mean").and_then(Value::as_f64).filter(|x| x.is_finite()).ok_or_else(|| {
                    Error::invalid(format!("sampling.parameters[{i}] normal distribution requires finite 'mean'"))
                })?;
                let sigma = match (po.get("sigma").and_then(Value::as_f64), po.get("rel_tol").and_then(Value::as_f64)) {
                    (Some(s), _) if s.is_finite() && s >= 0.0 => s,
                    (None, Some(r)) if r.is_finite() && r >= 0.0 => mean.abs() * r,
                    _ => {
                        return Err(Error::invalid(format!(
                            "sampling.parameters[{i}] normal distribution requires non-negative 'sigma' or 'rel_tol'"
                        )))
                    }
                };
                if let (Some(lo), Some(hi)) = (min_opt, max_opt) {
                    if lo > hi {
                        return Err(Error::invalid(format!(
                            "sampling.parameters[{i}] normal clamp min ({lo}) > max ({hi})"
                        )));
                    }
                }
                MarginDistribution::Normal { mean, sigma, min: min_opt, max: max_opt }
            }
            other => {
                return Err(Error::invalid(format!(
                    "sampling.parameters[{i}].distribution must be 'uniform' or 'normal', not '{other}'"
                )))
            }
        };

        dims.push(QmcParamDim { col_name, cmd, base_args, field, distribution });
    }

    let unit_table = match method {
        QmcMethod::Lhs => generate_lhs_unit(n_samples, dims.len(), seed)?,
        QmcMethod::Sobol => generate_sobol_unit(n_samples, dims.len(), seed, scramble, offset)?,
    };
    let samples = unit_table
        .into_iter()
        .map(|row| row.iter().zip(&dims).map(|(&u, d)| d.distribution.map_unit(u)).collect())
        .collect();

    Ok(QmcPlan { method, seed, scramble, offset, skip_invalid_geometry, dims, samples })
}

fn sample_edits_for(dims: &[QmcParamDim], coords: &[f64]) -> Vec<Value> {
    dims.iter()
        .zip(coords)
        .map(|(d, &v)| {
            let mut args = d.base_args.clone();
            args.insert(d.field.clone(), json!(v));
            json!({ "cmd": d.cmd, "args": args })
        })
        .collect()
}

fn sample_label_for(idx: usize, dims: &[QmcParamDim], coords: &[f64]) -> String {
    let summary: Vec<String> = dims.iter().zip(coords).take(3).map(|(d, &v)| format!("{}={v}", d.col_name)).collect();
    if dims.len() > 3 {
        format!("qmc-{idx}: {}, ...", summary.join(", "))
    } else {
        format!("qmc-{idx}: {}", summary.join(", "))
    }
}

impl Engine {
    /// A copy of this engine for editing: the project is a copy, the mesher, the registry and the material library are shared. Nothing it does
    /// reaches this engine.
    pub(crate) fn scratch(&self) -> Engine {
        let mut e = Engine::with(self.registry.clone(), self.mesher.clone());
        e.set_base_library_arc(self.base_library().clone());
        e.project = self.project.clone();
        e.owner = self.owner.clone();
        e.memory_wait_timeout = self.memory_wait_timeout;
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
                let n = if let Some(qjob) = sweep.qmc_job {
                    if sweep.qmc_summary.is_none() && self.jobs.cancel(qjob) {
                        sweep.points.len()
                    } else {
                        0
                    }
                } else {
                    sweep.points.iter().filter(|p| p.row.is_none()).filter(|p| self.jobs.cancel(p.job)).count()
                };
                Ok(json!({ "sweep": id, "cancelled": n }))
            }
            other => Err(Error::invalid(format!("no command '{other}'"))),
        }
    }

    /// Reads and checks a sweep request without cloning `Vec<Engine>`:
    /// - For explicit `points` / `series`, validates each point via one reusable scratch [`Engine`], records its
    ///   `mesh_key`, and drops the scratch engine before moving on.
    /// - For QMC `sampling` (`lhs` / `sobol`), validates distributions/commands and generates pure numeric coordinates.
    pub(crate) fn sweep_prepare(&self, args: &Value) -> Result<Prepared, Error> {
        let has_sampling = args.get("sampling").is_some_and(|v| !v.is_null());
        let has_points = args.get("points").is_some_and(|v| !v.is_null());
        let has_series = args.get("series").is_some_and(|v| !v.is_null());
        if has_sampling && (has_points || has_series) {
            return Err(Error::invalid("give `points`, `series`, or `sampling`, not multiple"));
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
        let memory_wait_timeout_s = crate::queue::memory_wait_timeout_of(args)?;
        let _ = crate::estimate::resolve_memory_wait_timeout(memory_wait_timeout_s, self.memory_wait_timeout)?;
        let opts = QueueOpts {
            priority: Some(priority),
            supersede: false,
            hold: args.get("hold").and_then(Value::as_bool).unwrap_or(false),
            isolate: args.get("isolate").and_then(Value::as_bool),
            force: args.get("force").and_then(Value::as_bool).unwrap_or(false),
            shared: args.get("shared").and_then(Value::as_bool),
            warm: None,
            save_fields: crate::queue::save_fields_of(args)?,
            post_analyses: crate::queue::post_analyses_of(args)?,
            memory_wait_timeout_s,
            retain_amr_data: false,
        };
        let warm = match args.get("reuse").and_then(Value::as_str) {
            None | Some("none") => false,
            Some("warm") => true,
            Some(other) => return Err(Error::invalid(format!("`reuse` is `none` or `warm`, not '{other}'"))),
        };
        let checkpoint_path = match args.get("checkpoint_path") {
            None | Some(Value::Null) => None,
            Some(Value::String(s)) if !s.trim().is_empty() => {
                if cfg!(not(feature = "hdf5")) {
                    return Err(Error::Unsupported("checkpoint_path requires a build with the `hdf5` feature".into()));
                }
                Some(s.clone())
            }
            Some(_) => return Err(Error::invalid("`checkpoint_path` must be a non-empty file path string")),
        };
        let checkpoint_lzf = args.get("checkpoint_lzf").and_then(Value::as_bool).unwrap_or(false);

        let allowed = self.point_commands();
        if has_sampling {
            let sampling_val = args.get("sampling").unwrap();
            let qmc = parse_qmc_plan(sampling_val, &allowed)?;
            if self.project.physics.is_empty() {
                return Err(Error::invalid("the project has no physics; add one first (add_physics)"));
            }
            self.library().map_err(|err| Error::invalid(format!("sampling: {err}")))?;
            for t in targets.as_deref().unwrap_or(&[]) {
                if self.project.spec(t).is_none() {
                    return Err(Error::invalid(format!("sampling: the project has no problem for '{t}'")));
                }
            }
            crate::queue::validate_post_analyses(&opts.post_analyses)?;
            return Ok(Prepared {
                asked: Vec::new(),
                targets,
                opts,
                warm: false,
                mesh_keys: Vec::new(),
                checkpoint_path,
                checkpoint_lzf,
                sampling: Some(qmc),
            });
        }

        let asked = asked_points(args)?;
        if asked.is_empty() {
            return Err(Error::invalid("a sweep needs at least one point"));
        }
        if asked.len() > MAX_POINTS {
            return Err(Error::invalid(format!("a sweep has at most {MAX_POINTS} points; split it")));
        }
        // Validate each point one-at-a-time via a temporary scratch engine without retaining Vec<Engine>.
        let mut mesh_keys = Vec::with_capacity(asked.len());
        for (i, a) in asked.iter().enumerate() {
            let e = self.point_engine(i, &a.edits, &allowed)?;
            e.library().map_err(|err| Error::invalid(format!("point {i}: {err}")))?;
            for t in targets.as_deref().unwrap_or(&[]) {
                if e.project.spec(t).is_none() {
                    return Err(Error::invalid(format!("point {i}: the project has no problem for '{t}'")));
                }
            }
            mesh_keys.push(e.mesh_key());
        }
        Ok(Prepared { asked, targets, opts, warm, mesh_keys, checkpoint_path, checkpoint_lzf, sampling: None })
    }

    fn sweep_run(&mut self, args: &Value) -> Result<Value, Error> {
        let Prepared { asked, targets, opts, warm, mesh_keys: _, checkpoint_path, checkpoint_lzf, sampling } =
            self.sweep_prepare(args)?;
        if let Some(qmc) = sampling {
            return self.sweep_run_qmc(args, qmc, targets, &opts, checkpoint_path, checkpoint_lzf);
        }

        let n_pts = asked.len();
        let cells: Vec<std::sync::Arc<crate::reuse::WarmCell>> =
            (0..n_pts).map(|_| crate::reuse::WarmCell::new()).collect();
        // a host that records (the MCP server) gets its workflow before the first run, as a solve does
        if self.auto_store && self.link.is_none() && self.workspace.is_some() {
            self.save_to_workspace(None)?;
        }
        self.next_sweep += 1;
        let id = self.next_sweep;
        let allowed = self.point_commands();
        let mut points: Vec<Point> = vec![];
        let mut resolved = vec![];
        for (i, a) in asked.into_iter().enumerate() {
            let e = self.point_engine(i, &a.edits, &allowed)?;
            // a warm sweep links each point to the one before: its modes are this point's start block (the first point only records)
            let mut point_opts = opts.clone();
            if warm {
                let source = match points.last() {
                    Some(prev) => {
                        crate::reuse::WarmSource::Cell { cell: cells[i - 1].clone(), job: prev.job, index: i - 1 }
                    }
                    None => crate::reuse::WarmSource::Cold,
                };
                point_opts.warm = Some(crate::reuse::WarmLink { source, mine: cells[i].clone(), index: i });
            }
            match self.submit_solve(Some(&e), targets.clone(), &point_opts, Purpose::Point { sweep: id, index: i }) {
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
        // what the sweep is expected to take, from the present mesh and the solves this machine has run (L3): the points are solved one at a time;
        // a sweep that changes the geometry or the mesh meshes every point anew, so this is the size of one point at the mesh that is current now
        let cost = self.model().map(|m| {
            let per_point = self.cost_json(&resolved, m.mesh().dim() as u8, m.mesh().n_cells() as u64);
            let seconds = per_point["solve"]["seconds"].as_f64().map(|s| s * listed.len() as f64);
            json!({ "per_point": per_point, "points": listed.len(), "seconds_estimate": seconds })
        });
        self.sweeps
            .insert(id, Sweep { name: name.clone(), targets: resolved, points, qmc_job: None, qmc_summary: None });
        Ok(json!({
            "sweep": id, "name": name, "points": listed, "cost": cost,
            "note": "the points are queued at their priority and run one at a time; sweep_table shows them as they finish, and the current project and results are not touched",
        }))
    }

    fn sweep_run_qmc(
        &mut self,
        args: &Value,
        qmc: QmcPlan,
        targets: Option<Vec<String>>,
        opts: &QueueOpts,
        checkpoint_path: Option<String>,
        checkpoint_lzf: bool,
    ) -> Result<Value, Error> {
        let resolved_targets = targets.unwrap_or_else(|| self.default_targets());
        if self.auto_store && self.link.is_none() && self.workspace.is_some() {
            self.save_to_workspace(None)?;
        }
        self.next_sweep += 1;
        let id = self.next_sweep;
        let name = args.get("name").and_then(Value::as_str).map(String::from).unwrap_or_else(|| format!("sweep {id}"));
        let allowed = self.point_commands();
        let chain = resolved_targets.clone();
        let mut submit =
            self.submission(&format!("qmc sweep {id}: {}", chain.join(", ")), Lane::Heavy, opts, &format!("qmc:{id}"));
        submit.priority = opts.priority.unwrap_or(3).clamp(1, 3);

        let base_project = self.project.clone();
        let registry = self.registry.clone();
        let mesher = self.mesher.clone();
        let base_library = self.base_library().clone();
        let stages = self.stages.clone();
        let retain_amr_data = opts.retain_amr_data || base_project.mesh.retain_amr_data;
        let post_analyses = opts.post_analyses.clone();
        let qmc_for_job = qmc.clone();
        let ckpt_for_job = checkpoint_path.clone();

        let job_id = self.jobs.submit(submit, move |progress| {
            run_streaming_qmc_job(
                progress,
                qmc_for_job,
                chain,
                base_project,
                registry,
                mesher,
                base_library,
                stages,
                allowed,
                post_analyses,
                retain_amr_data,
                ckpt_for_job,
                checkpoint_lzf,
            )
        });
        self.pending.insert(job_id, Pending::QmcSweep { sweep: id });

        let mut points = Vec::with_capacity(qmc.samples.len());
        let mut listed = Vec::with_capacity(qmc.samples.len());
        for (i, coords) in qmc.samples.iter().enumerate() {
            let label = sample_label_for(i, &qmc.dims, coords);
            let edits = sample_edits_for(&qmc.dims, coords);
            listed.push(json!({ "index": i, "label": label, "params": coords, "job": job_id }));
            points.push(Point { label, edits, job: job_id, row: None });
        }

        self.sweeps.insert(
            id,
            Sweep { name: name.clone(), targets: resolved_targets, points, qmc_job: Some(job_id), qmc_summary: None },
        );

        let dim_meta: Vec<Value> =
            qmc.dims.iter().map(|d| json!({ "name": d.col_name, "cmd": d.cmd, "field": d.field })).collect();
        Ok(json!({
            "sweep": id,
            "name": name,
            "points": listed,
            "checkpoint_path": checkpoint_path,
            "sampling": {
                "method": qmc.method.as_str(),
                "n_samples": qmc.samples.len(),
                "seed": qmc.seed,
                "scramble": qmc.scramble,
                "offset": qmc.offset,
                "skip_invalid_geometry": qmc.skip_invalid_geometry,
                "parameters": dim_meta,
            },
            "note": "QMC sweep streams one sample at a time on the heavy lane; use sweep_table with wait_s to collect results and skipped_invalid_geometry counts",
        }))
    }

    pub(crate) fn apply_qmc_sweep(&mut self, sweep: u64, outcome: Result<JobOutput, Error>) -> Result<Value, Error> {
        match outcome {
            Ok(JobOutput::QmcSweep(summary)) => {
                if let Some(sw) = self.sweeps.get_mut(&sweep) {
                    if let Some(rows) = summary.get("rows").and_then(Value::as_array) {
                        for (i, row) in rows.iter().enumerate() {
                            if let Some(pt) = sw.points.get_mut(i) {
                                pt.row = Some(row.clone());
                            }
                        }
                    }
                    sw.qmc_summary = Some(summary.clone());
                }
                Ok(summary)
            }
            Ok(_) => Err(Error::invalid("QMC sweep job returned an unexpected JobOutput variant")),
            Err(e) => Err(e),
        }
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
        let wait = match args.get("wait_s") {
            None | Some(Value::Null) => Duration::ZERO,
            Some(v) => {
                let t = v
                    .as_f64()
                    .filter(|t| *t >= 0.0 && t.is_finite())
                    .ok_or_else(|| Error::invalid("`wait_s` must be a finite non-negative number of seconds"))?;
                Duration::try_from_secs_f64(t)
                    .map(|d| d.min(Duration::from_secs(86_400)))
                    .map_err(|_| Error::invalid("`wait_s` overflows Duration"))?
            }
        };
        let (qmc_job_opt, jobs): (Option<JobId>, Vec<JobId>) = {
            let sw = self.sweeps.get(&id).ok_or_else(|| Error::invalid(format!("no sweep {id}")))?;
            (sw.qmc_job, sw.points.iter().filter(|p| p.row.is_none()).map(|p| p.job).collect())
        };
        let deadline = Instant::now() + wait;
        if let Some(qjob) = qmc_job_opt {
            if self.pending.contains_key(&qjob) {
                let outcome = if !wait.is_zero() {
                    self.jobs.wait_timeout(qjob, deadline.saturating_duration_since(Instant::now()))
                } else {
                    self.jobs.try_take(qjob)
                };
                if let Some(o) = outcome {
                    let _ = self.apply(qjob, o);
                }
            }
        } else {
            for job in jobs {
                if !self.pending.contains_key(&job) {
                    continue;
                }
                let outcome = if !wait.is_zero() {
                    self.jobs.wait_timeout(job, deadline.saturating_duration_since(Instant::now()))
                } else {
                    self.jobs.try_take(job)
                };
                if let Some(o) = outcome {
                    let _ = self.apply(job, o);
                }
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
        let mut resp = json!({
            "sweep": id, "name": sweep.name, "physics": sweep.targets,
            "points": rows, "finished": finished, "total": sweep.points.len(), "complete": finished == sweep.points.len(),
        });
        if let Some(ref summary) = sweep.qmc_summary {
            if let Some(obj) = resp.as_object_mut() {
                for key in
                    ["skipped_invalid_geometry", "skipped_points", "checkpoint_path", "solved", "failed", "cancelled"]
                {
                    if let Some(val) = summary.get(key) {
                        obj.insert(key.to_string(), val.clone());
                    }
                }
            }
        }
        Ok(resp)
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
        let mut metrics = Map::new();
        let mut reuse: Option<Value> = None;
        let mut pt_fields_saved = true;
        let mut pt_post_analyses: Vec<Value> = Vec::new();
        let (state, error, run, numbers, dofs, cells) = match outcome {
            Ok(JobOutput::Solved { model, results, failed, worker, reuse: reuse_of, fields_saved, post_analyses }) => {
                reuse = reuse_of;
                pt_fields_saved = fields_saved;
                pt_post_analyses = post_analyses;
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
                if let Some(r) = results.last() {
                    for (k, v) in crate::store_link::metrics_of(r) {
                        metrics.insert(k, json!(v));
                    }
                }
                match failed {
                    Some(e) => {
                        log.push(format!("stopped: {e}"));
                        let run = self.record_point(
                            &project,
                            RunOutcome::Failed(&e),
                            seconds,
                            &log,
                            Some(mesh),
                            reuse.as_ref(),
                            fields_saved,
                            &[],
                        );
                        ("failed", Some(e.to_string()), run, numbers, dofs, Some(cells))
                    }
                    None => {
                        self.learn(&project, &results, cells as u64, model.mesh().dim() as u8, &worker);
                        let run = self.record_point(
                            &project,
                            RunOutcome::Solved(&results),
                            seconds,
                            &log,
                            Some(mesh),
                            reuse.as_ref(),
                            fields_saved,
                            &pt_post_analyses,
                        );
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
                    self.record_point(&project, RunOutcome::Failed(&e), seconds, &log, None, None, true, &[])
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
            "numbers": numbers, "metrics": metrics, "unknowns": dofs, "cells": cells, "seconds": seconds, "run": run,
        });
        let mut row = row;
        if let Some(r) = reuse {
            row["reuse"] = r;
        }
        if !pt_fields_saved {
            row["fields_saved"] = json!(false);
        }
        if !pt_post_analyses.is_empty() {
            row["post_analyses"] = json!(pt_post_analyses);
        }
        if let Some(p) = self.sweeps.get_mut(&sweep).and_then(|s| s.points.get_mut(index)) {
            p.row = Some(row.clone());
        }
        Ok(row)
    }

    #[allow(clippy::too_many_arguments)]
    fn record_point(
        &mut self,
        project: &crate::project::Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<std::sync::Arc<femcore::mesh::SimplexMesh>>,
        reuse: Option<&Value>,
        fields_saved: bool,
        post_analyses: &[Value],
    ) -> Value {
        let script = self.script.clone();
        match self.link.as_mut() {
            None => Value::Null,
            Some(l) => {
                l.script = script;
                match l.record_point_with_options(
                    project,
                    outcome,
                    seconds,
                    log,
                    mesh.as_deref(),
                    reuse,
                    fields_saved,
                    post_analyses,
                ) {
                    Ok(r) => json!({ "id": r.id, "version": r.version, "status": r.status }),
                    Err(e) => json!({ "error": e.to_string() }),
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_streaming_qmc_job(
    progress: Arc<femcore::progress::Progress>,
    qmc: QmcPlan,
    chain: Vec<String>,
    base_project: crate::project::Project,
    registry: nyx_physics::Registry,
    mesher: Arc<dyn crate::meshing::Mesher>,
    base_library: Arc<femcore::material::Library>,
    stages: Arc<crate::queue::Stages>,
    allowed: Vec<String>,
    post_analyses: Vec<Value>,
    retain_amr_data: bool,
    checkpoint_path: Option<String>,
    checkpoint_lzf: bool,
) -> Result<JobOutput, Error> {
    let n_samples = qmc.samples.len();
    let mut rows = Vec::with_capacity(n_samples);
    let mut skipped_records = Vec::new();
    let mut n_skipped = 0usize;
    let mut n_solved = 0usize;
    let mut n_failed = 0usize;
    let mut n_cancelled = 0usize;

    #[cfg(feature = "hdf5")]
    let col_names: Vec<String> = qmc.dims.iter().map(|d| d.col_name.clone()).collect();
    #[cfg(feature = "hdf5")]
    let mut ckpt: Option<nyx_storage::SweepCheckpoint> = None;
    #[cfg(feature = "hdf5")]
    let mut pending_skipped: Vec<nyx_storage::SkippedSweepPoint> = Vec::new();
    #[cfg(not(feature = "hdf5"))]
    let _ = checkpoint_lzf;

    for (idx, coords) in qmc.samples.iter().enumerate() {
        let label = sample_label_for(idx, &qmc.dims, coords);
        let edits = sample_edits_for(&qmc.dims, coords);

        if progress.is_cancelled() {
            n_cancelled += 1;
            rows.push(json!({
                "index": idx, "label": label, "params": coords, "edits": edits,
                "state": "cancelled", "error": "cancelled"
            }));
            continue;
        }

        let _stage = femcore::progress::stage(format!("qmc sample {}/{n_samples}", idx + 1));
        let t_sample = Instant::now();

        // Build ONE scratch Engine per sample and drop it right after solving/skipping (`O(1)` memory).
        let mut sample_eng = Engine::with(registry.clone(), mesher.clone());
        sample_eng.set_base_library_arc(base_library.clone());
        sample_eng.project = base_project.clone();

        let pre_check: Result<_, Error> = (|| {
            for edit in &edits {
                let cmd =
                    edit.get("cmd").and_then(Value::as_str).ok_or_else(|| Error::invalid("missing QMC edit cmd"))?;
                if !allowed.iter().any(|a| a == cmd) {
                    return Err(Error::invalid(format!("'{cmd}' is not an edit command")));
                }
                let args = edit.get("args").cloned().unwrap_or_else(|| json!({}));
                sample_eng.call(cmd, &args)?;
            }
            sample_eng.project.check(&sample_eng.registry)?;
            let library = sample_eng.library()?;
            let mesh_key = sample_eng.mesh_key();
            let model = match stages.mesh(&mesh_key) {
                Some(m) => m,
                None => {
                    let sizing = sample_eng.sizing();
                    let m = Arc::new(mesher.mesh(
                        &sample_eng.project.geometry,
                        &sample_eng.project.mesh,
                        sizing.as_ref(),
                    )?);
                    // Validate Exchange Mesh AND SimplexMesh integrity before accepting geometry.
                    nyx_physics::check_mesh(m.mesh())?;
                    let simplex = m.mesh().discretise().map_err(Error::Invalid)?;
                    simplex.check().map_err(|e| Error::invalid(format!("SimplexMesh::check failed: {e}")))?;
                    stages.put_mesh(&mesh_key, &m);
                    m
                }
            };
            Ok((model, library))
        })();

        let (model, library) = match pre_check {
            Ok(pair) => pair,
            Err(err) => {
                if qmc.skip_invalid_geometry {
                    let reason = err.to_string();
                    n_skipped += 1;
                    let sk_json = json!({
                        "index": idx,
                        "label": label,
                        "params": coords,
                        "reason": reason,
                        "status": "skipped_invalid_geometry",
                    });
                    skipped_records.push(sk_json);
                    #[cfg(feature = "hdf5")]
                    {
                        let rec = nyx_storage::SkippedSweepPoint {
                            index: idx,
                            label: label.clone(),
                            params: coords.clone(),
                            reason: reason.clone(),
                        };
                        if let Some(ref c) = ckpt {
                            c.record_skipped(&rec)
                                .map_err(|e| Error::invalid(format!("checkpoint skipped write: {e}")))?;
                        } else {
                            pending_skipped.push(rec);
                        }
                    }
                    rows.push(json!({
                        "index": idx,
                        "label": label,
                        "params": coords,
                        "edits": edits,
                        "state": "skipped_invalid_geometry",
                        "reason": reason,
                        "seconds": t_sample.elapsed().as_secs_f64(),
                    }));
                    continue;
                } else {
                    return Err(err);
                }
            }
        };

        let (solved_model, mut results, failed) = crate::solve::run_chain_partial_amr_with(
            &registry,
            &sample_eng.project,
            model,
            &chain,
            &library,
            Arc::clone(&progress),
            retain_amr_data,
        );

        let cells = solved_model.mesh().n_cells();
        let mut numbers = Map::new();
        let mut dofs = Map::new();
        for r in &results {
            let nums: Map<String, Value> =
                crate::convergence::numbers(r).into_iter().map(|(k, v)| (k, json!(v))).collect();
            numbers.insert(r.physics.clone(), Value::Object(nums));
            dofs.insert(r.physics.clone(), json!(r.dofs));
        }
        let mut metrics = Map::new();
        if let Some(r) = results.last() {
            for (k, v) in crate::store_link::metrics_of(r) {
                metrics.insert(k, json!(v));
            }
        }

        let mut pt_post_analyses = Vec::new();
        if failed.is_none() && !post_analyses.is_empty() {
            let solve_keys = sample_eng.problem_keys();
            let mesh_key = sample_eng.mesh_key();
            let mut tmp = Engine::with(registry.clone(), mesher.clone());
            tmp.set_base_library_arc(base_library.clone());
            tmp.project = sample_eng.project.clone();
            tmp.model = Some((mesh_key, solved_model.clone()));
            let now = Instant::now();
            for res in std::mem::take(&mut results) {
                let p_code = res.physics.clone();
                let k_code = solve_keys.get(&p_code).cloned().unwrap_or_default();
                tmp.results.insert(
                    p_code,
                    crate::engine::Stored {
                        result: Arc::new(res),
                        key: k_code,
                        job: 0,
                        at: now,
                        model: solved_model.clone(),
                        run: None,
                    },
                );
            }
            pt_post_analyses = tmp.run_post_analyses(&post_analyses, &chain)?;
            let mut restored: Vec<nyx_physics::SolveResult> = tmp
                .results
                .into_values()
                .map(|st| Arc::try_unwrap(st.result).unwrap_or_else(|a| (*a).clone()))
                .collect();
            let order: Vec<String> = registry
                .solve_order(&chain.iter().map(String::as_str).collect::<Vec<_>>())
                .unwrap_or_else(|_| chain.clone());
            restored.sort_by_key(|res| order.iter().position(|o| o == &res.physics).unwrap_or(usize::MAX));
            results = restored;
        }
        #[cfg(not(feature = "hdf5"))]
        let _ = &results;

        // Stream solved point incrementally into HDF5 SweepCheckpoint when configured.
        #[cfg(feature = "hdf5")]
        if let Some(ref cpath) = checkpoint_path {
            let (eig_re_raw, eig_im_raw): (Vec<f64>, Vec<f64>) = results
                .last()
                .map(|r| {
                    let eigs: Vec<[f64; 2]> = r.solutions.iter().filter_map(|s| s.eigenvalue).collect();
                    if !eigs.is_empty() {
                        (eigs.iter().map(|e| e[0]).collect(), eigs.iter().map(|e| e[1]).collect())
                    } else {
                        let vals: Vec<f64> = r.solutions.iter().filter_map(|s| s.value).collect();
                        let zeros = vec![0.0; vals.len()];
                        (vals, zeros)
                    }
                })
                .unwrap_or_default();

            if ckpt.is_none() {
                let col_refs: Vec<&str> = col_names.iter().map(String::as_str).collect();
                let created = nyx_storage::SweepCheckpoint::create_with_options(
                    std::path::Path::new(cpath),
                    &col_refs,
                    eig_re_raw.len(),
                    checkpoint_lzf,
                )
                .map_err(|e| Error::invalid(format!("failed to create SweepCheckpoint '{cpath}': {e}")))?;
                if !pending_skipped.is_empty() {
                    created
                        .record_skipped_batch(&pending_skipped)
                        .map_err(|e| Error::invalid(format!("failed to record leading skipped points: {e}")))?;
                    pending_skipped.clear();
                }
                ckpt = Some(created);
            }
            if let Some(ref c) = ckpt {
                let n_m = c.n_modes();
                let mut eig_re = eig_re_raw;
                let mut eig_im = eig_im_raw;
                eig_re.resize(n_m, f64::NAN);
                eig_im.resize(n_m, f64::NAN);
                c.append_point(&nyx_storage::SweepPoint { params: coords.clone(), eig_re, eig_im })
                    .map_err(|e| Error::invalid(format!("failed to stream SweepPoint to '{cpath}': {e}")))?;
            }
        }

        let elapsed = t_sample.elapsed().as_secs_f64();
        let (state, err_str) = match failed {
            None => {
                n_solved += 1;
                ("done", None)
            }
            Some(e) => {
                n_failed += 1;
                ("failed", Some(e.to_string()))
            }
        };
        let mut row = json!({
            "index": idx,
            "label": label,
            "params": coords,
            "edits": edits,
            "state": state,
            "error": err_str,
            "numbers": numbers,
            "metrics": metrics,
            "unknowns": dofs,
            "cells": cells,
            "seconds": elapsed,
        });
        if !pt_post_analyses.is_empty() {
            row["post_analyses"] = json!(pt_post_analyses);
        }
        rows.push(row);
    }

    // If every QMC sample was skipped before any solve initialized the checkpoint, materialize an n_modes=0 checkpoint with all skipped entries.
    #[cfg(feature = "hdf5")]
    if let Some(ref cpath) = checkpoint_path {
        if ckpt.is_none() {
            let col_refs: Vec<&str> = col_names.iter().map(String::as_str).collect();
            let created = nyx_storage::SweepCheckpoint::create_with_options(
                std::path::Path::new(cpath),
                &col_refs,
                0,
                checkpoint_lzf,
            )
            .map_err(|e| Error::invalid(format!("failed to create SweepCheckpoint '{cpath}': {e}")))?;
            if !pending_skipped.is_empty() {
                created
                    .record_skipped_batch(&pending_skipped)
                    .map_err(|e| Error::invalid(format!("failed to flush skipped points: {e}")))?;
            }
        }
    }

    Ok(JobOutput::QmcSweep(json!({
        "rows": rows,
        "skipped_invalid_geometry": n_skipped,
        "skipped_points": skipped_records,
        "checkpoint_path": checkpoint_path,
        "solved": n_solved,
        "failed": n_failed,
        "cancelled": n_cancelled,
    })))
}

fn sweep_id(args: &Value) -> Result<u64, Error> {
    args.get("sweep").and_then(Value::as_u64).ok_or_else(|| Error::invalid("`sweep` is required and a whole number"))
}
