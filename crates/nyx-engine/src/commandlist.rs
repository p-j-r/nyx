//! Command lists (productivity plan S1): a file of JSON lines, `{"cmd": "<command of Engine::catalog()>", "args": {...}}`, that a person or an agent
//! writes once and the engine checks, plans and runs. A list is a third front-end over `Engine::call` (ADR 57), beside the window and the MCP
//! server: nothing here knows a command's meaning except through the catalog and the few names that start work.
//!
//! * `parse` reads the lines (blank lines and `#` lines are skipped) and names every line that is not a command.
//! * `Engine::check_list` is the **static** stage (every `cmd` is in the catalog, every `args` fits that command's JSON Schema) and then the **plan**
//!   (a dry run on a scratch copy of the engine: edits are applied for real, so a feature that does not exist or a problem that does not check
//!   is found at its line; what would be meshed and solved is counted, with the meshes the stage cache would reuse; nothing is meshed, solved or
//!   written). It starts no job.
//! * `run` sends the lines to an engine one after another and reports each to an observer. A solve inside it goes through the queue like any
//!   other (`solve_run` waits for its job).
//!
//! Besides `cmd` and `args` a line may carry keys that only a replay uses: `note` (a label), `job_arg` (the argument that receives the id of the last
//! job a reply named), `after_stage` and `delay_s` (wait until a job's stage contains a text, then a time more: the way to cancel "during the
//! factorisation"), `expect_error` (an error reply is the expected outcome). The workflows of `tools/workflows/` use all of them.

use crate::commands::Access;
use crate::engine::Engine;
use nyx_physics::Error;
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::time::{Duration, Instant};

const KEYS: [&str; 7] = ["cmd", "args", "note", "job_arg", "after_stage", "delay_s", "expect_error"];

/// One command of a list.
#[derive(Clone, Debug)]
pub struct Line {
    /// the line number in the file (from 1)
    pub number: usize,
    pub cmd: String,
    pub args: Value,
    pub note: Option<String>,
    pub job_arg: Option<String>,
    pub after_stage: Option<String>,
    pub delay_s: Option<f64>,
    pub expect_error: bool,
}

/// What is wrong at a line, in words.
#[derive(Clone, Debug, Serialize)]
pub struct Problem {
    pub line: usize,
    pub cmd: String,
    pub message: String,
}

impl std::fmt::Display for Problem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.cmd.is_empty() {
            write!(f, "line {}: {}", self.line, self.message)
        } else {
            write!(f, "line {}, {}: {}", self.line, self.cmd, self.message)
        }
    }
}

/// One line of the plan: what the line would do.
#[derive(Clone, Debug, Serialize)]
pub struct Step {
    pub line: usize,
    pub cmd: String,
    /// `edit` (applied on the copy), `mesh`, `solve`, `convergence`, `sweep`, `read`, `queue`, `files` or `not planned`
    pub kind: &'static str,
    pub summary: String,
}

/// The outcome of `check_list`.
#[derive(Clone, Debug, Serialize)]
pub struct Report {
    /// commands in the list
    pub commands: usize,
    /// the hash of the text, a cache key (ADR 64)
    pub hash: String,
    /// what is wrong; when it is not empty nothing was planned
    pub problems: Vec<Problem>,
    pub steps: Vec<Step>,
    /// meshes the list would build (a mesh the cache holds is not counted) and solve jobs and mesh builds of a convergence check
    pub meshes: usize,
    pub solves: usize,
}

impl Report {
    pub fn ok(&self) -> bool {
        self.problems.is_empty()
    }
}

/// Reads a list. Problems come back with the lines that parsed, so a front-end can show them all at once.
pub fn parse(text: &str) -> (Vec<Line>, Vec<Problem>) {
    let (mut lines, mut problems) = (vec![], vec![]);
    for (i, raw) in text.lines().enumerate() {
        let n = i + 1;
        let raw = raw.trim();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let bad = |m: String| Problem { line: n, cmd: String::new(), message: m };
        let v: Value = match serde_json::from_str(raw) {
            Ok(v) => v,
            Err(e) => {
                problems.push(bad(format!("not JSON: {e}")));
                continue;
            }
        };
        let Some(o) = v.as_object() else {
            problems.push(bad("a line is a JSON object: { \"cmd\": ..., \"args\": {...} }".into()));
            continue;
        };
        let Some(cmd) = o.get("cmd").and_then(Value::as_str) else {
            problems.push(bad("no `cmd`: the name of a command (list_commands lists them)".into()));
            continue;
        };
        let at = |m: String| Problem { line: n, cmd: cmd.to_string(), message: m };
        if let Some(k) = o.keys().find(|k| !KEYS.contains(&k.as_str())) {
            problems.push(at(format!("unknown key '{k}' (a line has {})", KEYS.join(", "))));
            continue;
        }
        let args = match o.get("args") {
            None | Some(Value::Null) => json!({}),
            Some(a) if a.is_object() => a.clone(),
            Some(_) => {
                problems.push(at("`args` is an object".into()));
                continue;
            }
        };
        let text_of = |k: &str| o.get(k).and_then(Value::as_str).map(String::from);
        lines.push(Line {
            number: n,
            cmd: cmd.to_string(),
            args,
            note: text_of("note"),
            job_arg: text_of("job_arg"),
            after_stage: text_of("after_stage"),
            delay_s: o.get("delay_s").and_then(Value::as_f64),
            expect_error: o.get("expect_error").and_then(Value::as_bool).unwrap_or(false),
        });
    }
    (lines, problems)
}

// ------------------------------------------------------------------------------------------------------------------ static stage

fn type_ok(t: &str, v: &Value) -> bool {
    match t {
        "string" => v.is_string(),
        "number" => v.is_number(),
        "integer" => v.is_i64() || v.is_u64() || v.as_f64().is_some_and(|x| x.fract() == 0.0),
        "boolean" => v.is_boolean(),
        "object" => v.is_object(),
        "array" => v.is_array(),
        "null" => v.is_null(),
        _ => true,
    }
}

fn kind_of(v: &Value) -> &'static str {
    match v {
        Value::Null => "null",
        Value::Bool(_) => "a boolean",
        Value::Number(_) => "a number",
        Value::String(_) => "text",
        Value::Array(_) => "a list",
        Value::Object(_) => "an object",
    }
}

/// Checks `v` against the keywords of the JSON Schemas the catalog uses (`type`, `enum`, `oneOf`, `properties`, `required`,
/// `additionalProperties`, `items`); a keyword it does not know is not checked, so a valid command is never refused.
fn validate(schema: &Value, v: &Value, at: &str, out: &mut Vec<String>) {
    if let Some(alts) = schema.get("oneOf").and_then(Value::as_array) {
        let fits = alts.iter().any(|a| {
            let mut o = vec![];
            validate(a, v, at, &mut o);
            o.is_empty()
        });
        if !fits {
            out.push(format!("{at}: is {}, which none of the allowed forms accepts", kind_of(v)));
        }
        return;
    }
    let types: Vec<&str> = match schema.get("type") {
        Some(Value::String(t)) => vec![t.as_str()],
        Some(Value::Array(a)) => a.iter().filter_map(Value::as_str).collect(),
        _ => vec![],
    };
    if !types.is_empty() && !types.iter().any(|t| type_ok(t, v)) {
        out.push(format!("{at}: must be {}, not {}", types.join(" or "), kind_of(v)));
        return;
    }
    if let Some(allowed) = schema.get("enum").and_then(Value::as_array) {
        if !allowed.contains(v) {
            let shown: Vec<String> = allowed
                .iter()
                .take(12)
                .map(|a| a.as_str().map(String::from).unwrap_or_else(|| a.to_string()))
                .collect();
            out.push(format!("{at}: {v} is not one of {}", shown.join(", ")));
        }
        return;
    }
    if let (Some(o), Some(props)) = (v.as_object(), schema.get("properties").and_then(Value::as_object)) {
        for r in schema.get("required").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str) {
            if !o.contains_key(r) {
                out.push(format!("{at}: missing '{r}'"));
            }
        }
        for (k, x) in o {
            match props.get(k) {
                Some(s) => {
                    let required =
                        schema.get("required").and_then(Value::as_array).is_some_and(|r| r.iter().any(|n| n == k));
                    if !(x.is_null() && !required) {
                        validate(s, x, &format!("{at}.{k}"), out);
                    }
                }
                None if schema.get("additionalProperties") == Some(&json!(false)) => {
                    let names: Vec<&str> = props.keys().map(String::as_str).collect();
                    out.push(format!("{at}: no argument '{k}' (the arguments are {})", names.join(", ")));
                }
                None => {}
            }
        }
    }
    if let (Some(a), Some(items)) = (v.as_array(), schema.get("items")) {
        for (i, x) in a.iter().enumerate() {
            validate(items, x, &format!("{at}[{i}]"), out);
        }
    }
}

// ------------------------------------------------------------------------------------------------------------------ plan stage

/// The meshes the stage cache holds (the last two), as a plan sees them.
struct Cache {
    meshes: Vec<String>,
}

impl Cache {
    /// True when `key` has to be built.
    fn need(&mut self, key: String) -> bool {
        if self.meshes.contains(&key) {
            return false;
        }
        self.meshes.push(key);
        if self.meshes.len() > 2 {
            self.meshes.remove(0);
        }
        true
    }
}

fn physics_arg(args: &Value) -> Option<Vec<String>> {
    match args.get("physics") {
        Some(Value::String(s)) => Some(vec![s.clone()]),
        Some(Value::Array(a)) => Some(a.iter().filter_map(|x| x.as_str().map(String::from)).collect()),
        _ => None,
    }
}

impl Engine {
    /// Checks a list: first every line against the catalog, then, when those pass, a dry run on a copy of this engine. Nothing is meshed, solved or
    /// written and no job starts; this engine is not changed.
    pub fn check_list(&self, text: &str) -> Report {
        let (lines, mut problems) = parse(text);
        let mut report = Report {
            commands: lines.len(),
            hash: nyx_store::util::hash_bytes(text.as_bytes()),
            problems: vec![],
            steps: vec![],
            meshes: 0,
            solves: 0,
        };
        problems.extend(self.check_static(&lines));
        problems.sort_by_key(|p| p.line);
        if !problems.is_empty() {
            report.problems = problems;
            return report;
        }
        self.plan_into(&lines, &mut report);
        report
    }

    /// The static stage: every command is in the catalog and its arguments fit its schema.
    pub fn check_static(&self, lines: &[Line]) -> Vec<Problem> {
        let catalog: BTreeMap<String, Value> = self.catalog().into_iter().map(|c| (c.name, c.schema)).collect();
        let mut out = vec![];
        for l in lines {
            let Some(schema) = catalog.get(&l.cmd) else {
                let names: Vec<String> = catalog.keys().cloned().collect();
                let near = crate::apiref::nearest(&names, &l.cmd);
                let hint = if near.is_empty() {
                    "list_commands lists them".to_string()
                } else {
                    format!("did you mean {}?", near.join(", "))
                };
                out.push(Problem {
                    line: l.number,
                    cmd: l.cmd.clone(),
                    message: format!("there is no command '{}': {hint}", l.cmd),
                });
                continue;
            };
            // the argument a `job_arg` fills in at run time is there when the line runs
            let mut args = l.args.clone();
            if let Some(name) = &l.job_arg {
                args[name] = json!(0);
            }
            let mut found = vec![];
            validate(schema, &args, "args", &mut found);
            out.extend(found.into_iter().map(|m| Problem { line: l.number, cmd: l.cmd.clone(), message: m }));
        }
        out
    }

    fn plan_into(&self, lines: &[Line], report: &mut Report) {
        let access: BTreeMap<String, Access> = self.catalog().into_iter().map(|c| (c.name, c.access)).collect();
        let mut e = self.scratch();
        let mut cache = Cache { meshes: vec![] };
        let mut solved: BTreeSet<String> = BTreeSet::new();
        // a physics solved alone on a mesh, by problem key: the coarse solve a convergence check then reuses
        let mut alone: BTreeSet<String> = BTreeSet::new();
        let mut results_unknown = false;
        for l in lines {
            let problem = |m: String| Problem { line: l.number, cmd: l.cmd.clone(), message: m };
            let step = |kind: &'static str, summary: String| Step { line: l.number, cmd: l.cmd.clone(), kind, summary };
            // `queue_submit` is the solve, mesh or convergence it names
            let what = match l.cmd.as_str() {
                "queue_submit" => match l.args.get("kind").and_then(Value::as_str).unwrap_or("solve") {
                    "mesh" => "mesh_build",
                    "convergence" => "convergence_check",
                    _ => "solve_run",
                },
                other => other,
            };
            match what {
                "solve_run" => {
                    let targets = physics_arg(&l.args).unwrap_or_else(|| e.default_targets());
                    if let Err(m) = self.plan_ready(&e, &targets) {
                        report.problems.push(problem(m));
                        continue;
                    }
                    let new_mesh = cache.need(e.mesh_key());
                    report.meshes += usize::from(new_mesh);
                    report.solves += 1;
                    if targets.len() == 1 {
                        if let Some(k) = e.problem_key(&targets[0]) {
                            alone.insert(k);
                        }
                    }
                    solved.extend(targets.iter().cloned());
                    report.steps.push(step(
                        "solve",
                        format!("solve {} (mesh {})", targets.join(" + "), if new_mesh { "built" } else { "reused" }),
                    ));
                }
                "mesh_build" => {
                    if let Err(m) = e.project.geometry.validate().map_err(|err| format!("geometry: {err}")) {
                        report.problems.push(problem(m));
                        continue;
                    }
                    let new_mesh = cache.need(e.mesh_key());
                    report.meshes += usize::from(new_mesh);
                    report.steps.push(step("mesh", format!("mesh ({})", if new_mesh { "built" } else { "reused" })));
                }
                "convergence_check" => {
                    if l.args.get("dry_run").and_then(serde_json::Value::as_bool) == Some(true) {
                        report.steps.push(step("convergence", "convergence check, dry run: both meshes are built and their cost is told; nothing is solved".to_string()));
                        continue;
                    }
                    let Some(p) = physics_arg(&l.args).and_then(|v| v.into_iter().next()) else {
                        report.problems.push(problem("a convergence check needs `physics`".into()));
                        continue;
                    };
                    if let Err(m) = self.plan_ready(&e, std::slice::from_ref(&p)) {
                        report.problems.push(problem(m));
                        continue;
                    }
                    let coarse_mesh = cache.need(e.mesh_key());
                    let reused = e.problem_key(&p).is_some_and(|k| alone.contains(&k));
                    report.meshes += usize::from(coarse_mesh) + 1;
                    report.solves += if reused { 1 } else { 2 };
                    solved.insert(p.clone());
                    report.steps.push(step(
                        "convergence",
                        format!(
                            "convergence check of {p}: the finer mesh and its solve{}",
                            if reused {
                                " (the solve on the present mesh is reused)"
                            } else {
                                ", and the solve on the present mesh"
                            }
                        ),
                    ));
                }
                "sweep_run" => match e.sweep_prepare(&l.args) {
                    Err(err) => report.problems.push(problem(err.to_string())),
                    Ok(p) => {
                        let mut built = 0;
                        for pe in &p.engines {
                            built += usize::from(cache.need(pe.mesh_key()));
                        }
                        report.meshes += built;
                        report.solves += p.engines.len();
                        solved.extend(p.targets.clone().unwrap_or_else(|| e.default_targets()));
                        report.steps.push(step(
                            "sweep",
                            format!("{} point(s), {} mesh(es) to build, one recorded run each", p.engines.len(), built),
                        ));
                    }
                },
                "open_results" => {
                    // it replaces the project and the results with a file's; the lines after it cannot be checked against that
                    results_unknown = true;
                    report.steps.push(step(
                        "not planned",
                        "opens a results file: the project and results after it are the file's".into(),
                    ));
                }
                name => {
                    // a command that reads a result has to have one to read
                    if matches!(
                        name,
                        "result_info"
                            | "result_field"
                            | "analyse"
                            | "export_results"
                            | "field"
                            | "field_op"
                            | "integrate"
                            | "extremum"
                            | "probe"
                            | "line_cut"
                            | "overlap"
                            | "confinement"
                    ) && !results_unknown
                    {
                        if let Some(ps) = physics_arg(&l.args) {
                            if let Some(p) = ps.iter().find(|p| !solved.contains(*p)) {
                                report.problems.push(problem(format!(
                                    "there is no result of '{p}' at this point: nothing before this line solves it"
                                )));
                                continue;
                            }
                        }
                    }
                    match access.get(name) {
                        Some(Access::Edit) => match e.call(name, &l.args) {
                            Ok(_) => {
                                // a command that replaces the project drops the mesh, the results and the stage caches, as it does for real
                                if matches!(name, "project_new" | "project_load" | "project_open" | "store_open") {
                                    solved.clear();
                                    alone.clear();
                                    cache.meshes.clear();
                                }
                                report.steps.push(step("edit", "applied to the copy".into()))
                            }
                            Err(err) => report.problems.push(problem(err.to_string())),
                        },
                        Some(Access::Read) => report.steps.push(step("read", "not run in the plan".into())),
                        _ => {
                            let kind =
                                if name.starts_with("job_") || name.starts_with("queue_") || name.starts_with("sweep_")
                                {
                                    "queue"
                                } else {
                                    "files"
                                };
                            report.steps.push(step(kind, "not run in the plan".into()));
                        }
                    }
                }
            }
        }
    }

    /// Whether the copy's project can be solved for `targets` now: it has the physics and `Project::check` passes.
    fn plan_ready(&self, e: &Engine, targets: &[String]) -> Result<(), String> {
        if e.project.physics.is_empty() {
            return Err("the project has no physics at this point (nothing before this line adds one)".into());
        }
        for t in targets {
            if e.project.spec(t).is_none() {
                return Err(format!("the project has no problem for '{t}' at this point"));
            }
        }
        e.project.check(&e.registry).map_err(|err| err.to_string())
    }
}

// ------------------------------------------------------------------------------------------------------------------ run

/// What the observer of `run` is told after each command.
pub struct Done<'a> {
    pub line: &'a Line,
    pub reply: &'a Result<Value, Error>,
    pub seconds: f64,
    /// the stage of the awaited job when `after_stage` released the command
    pub stage: &'a str,
    /// seconds from a `job_cancel` to the end of the `job_wait` that follows it
    pub since_cancel: Option<f64>,
}

/// Why a run stopped.
#[derive(Clone, Debug)]
pub struct RunError {
    pub line: usize,
    pub cmd: String,
    pub message: String,
}

impl std::fmt::Display for RunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "line {}, {}: {}", self.line, self.cmd, self.message)
    }
}

/// Sends the lines to `engine` in order. The run stops at the first command that fails (unless the line expects an error) and at a replay key that
/// cannot be honoured (a cancel that found no running job, an `expect_error` command that succeeded), so a step that did not happen cannot pass.
pub fn run(engine: &mut Engine, lines: &[Line], observe: &mut dyn FnMut(&Done)) -> Result<(), RunError> {
    let (mut last_job, mut cancelled_at): (Option<u64>, Option<Instant>) = (None, None);
    for l in lines {
        let stop = |m: String| RunError { line: l.number, cmd: l.cmd.clone(), message: m };
        let mut args = l.args.clone();
        if let Some(name) = &l.job_arg {
            args[name] = json!(last_job.ok_or_else(|| stop("`job_arg`: no job was started before this line".into()))?);
        }
        let mut stage = String::new();
        if let Some(want) = &l.after_stage {
            let job = last_job.ok_or_else(|| stop("`after_stage`: no job was started before this line".into()))?;
            let t = Instant::now();
            loop {
                let info = engine.jobs().info(job).ok_or_else(|| stop("`after_stage`: the job is gone".into()))?;
                stage = info.stage.clone();
                if !info.state.is_open() {
                    return Err(stop(format!(
                        "`after_stage`: the job ended before stage '{want}' (last stage seen: '{stage}')"
                    )));
                }
                if stage.contains(want.as_str()) {
                    break;
                }
                if t.elapsed() > Duration::from_secs(900) {
                    return Err(stop(format!("`after_stage`: stage '{want}' not reached in 15 minutes")));
                }
                std::thread::sleep(Duration::from_millis(2));
            }
            if let Some(d) = l.delay_s {
                std::thread::sleep(Duration::from_secs_f64(d));
                let info = engine.jobs().info(job).ok_or_else(|| stop("`after_stage`: the job is gone".into()))?;
                if !info.state.is_open() {
                    return Err(stop("`delay_s`: the job ended during the delay".into()));
                }
                stage = info.stage;
            }
        }
        let t0 = Instant::now();
        if l.cmd == "job_cancel" {
            cancelled_at = Some(t0);
        }
        let reply = engine.call(&l.cmd, &args);
        let seconds = t0.elapsed().as_secs_f64();
        let since_cancel =
            if l.cmd == "job_wait" { cancelled_at.take().map(|t| t.elapsed().as_secs_f64()) } else { None };
        if let Ok(v) = &reply {
            if let Some(j) = v.get("job").and_then(Value::as_u64) {
                last_job = Some(j);
            }
        }
        observe(&Done { line: l, reply: &reply, seconds, stage: &stage, since_cancel });
        if l.cmd == "job_cancel" && reply.as_ref().is_ok_and(|v| v["cancelled"] != true) {
            return Err(stop("`job_cancel`: the job was no longer running, nothing was cancelled".into()));
        }
        match (&reply, l.expect_error) {
            (Ok(_), true) => {
                return Err(stop("`expect_error`: the command succeeded (a cancelled job completed instead)".into()))
            }
            (Err(e), false) => return Err(stop(e.to_string())),
            _ => {}
        }
    }
    Ok(())
}
