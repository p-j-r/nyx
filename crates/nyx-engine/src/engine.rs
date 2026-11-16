//! The engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them and the binding to a
//! workflow folder. The window and the MCP server hold an `Engine` and nothing else; every action is a method here or a named command
//! (`commands.rs`, which calls these methods).
//!
//! The mesh and the results are caches keyed by the content they were made from (a hash of the geometry and mesh settings; for a result the
//! *problem key* of its physics: its own problem, the project materials, the mesh key and the keys of the physics it takes fields from), so an edit
//! never has to remember to invalidate them: a result whose key is not the one its physics has now is *stale* and says so. A key is a cache key
//! (ADR 64, ADR 84): it is never stored in a file or sent to a front-end.

use crate::jobs::{JobId, JobOutput, JobRunner, Lane, Opened};
use crate::meshing::{FemcoreMesher, MeshedModel, Mesher};
use crate::project::Project;
use crate::queue::{solve_key, JobMeta, QueueOpts, SolveSummary, Stages};
use crate::store_link::{RunOutcome, StoreLink};
use nyx_physics::{Error, Registry, Sizing, SolveResult};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

pub(crate) struct Stored {
    pub result: Arc<SolveResult>,
    /// the problem key of its physics when the solve was started
    pub key: String,
    /// the job that made it (jobs are numbered in the order they start; 0 for results read from a file)
    pub job: JobId,
    /// when it was stored
    pub at: Instant,
}

pub(crate) enum Pending {
    Mesh { key: String, generation: u64 },
    Solve { key: String, keys: BTreeMap<String, String>, project: Box<Project>, generation: u64 },
    Convergence,
    Open,
    Export,
}

pub struct Engine {
    pub(crate) registry: Registry,
    pub(crate) mesher: Arc<dyn Mesher>,
    pub(crate) project: Project,
    pub(crate) model: Option<(String, Arc<MeshedModel>)>,
    pub(crate) results: BTreeMap<String, Stored>,
    pub(crate) jobs: JobRunner,
    pub(crate) pending: BTreeMap<JobId, Pending>,
    /// what each job was asked for (kind, targets, the project it was made for), for the queue commands
    pub(crate) meta: BTreeMap<JobId, JobMeta>,
    /// recent meshes and results by content, shared by queued requests (ADR 87)
    pub(crate) stages: Arc<Stages>,
    /// who this engine's requests are made for: `gui`, `ai` or `cli`
    pub(crate) owner: String,
    /// counts the projects this engine has held (`set_project`): a job started for an earlier one must not put its outcome into the present one
    pub(crate) generation: u64,
    pub(crate) link: Option<StoreLink>,
    pub(crate) workspace: Option<PathBuf>,
    /// bind the project to a workflow of `workspace` at the first solve, so what a host does shows in the dashboard (the MCP server turns it on)
    pub(crate) auto_store: bool,
    /// the material files of the program (`materials/`), before the project's own
    pub(crate) base_library: Arc<femcore::material::Library>,
    /// why the library is empty or partial, in words
    pub library_notes: Vec<String>,
}

fn hash_of(text: &str) -> String {
    nyx_store::util::hash_bytes(text.as_bytes())
}

impl Engine {
    /// An engine with the physics compiled into this build and an empty 2D project.
    pub fn new() -> Engine {
        let mut e = Engine::with(crate::registry(), Arc::new(FemcoreMesher));
        let (lib, notes) = crate::materials::load_library(&crate::materials::default_dirs());
        e.base_library = Arc::new(lib);
        e.library_notes = notes;
        e
    }

    pub fn with(registry: Registry, mesher: Arc<dyn Mesher>) -> Engine {
        Engine {
            registry,
            mesher,
            project: Project::new("Untitled", 2),
            model: None,
            results: BTreeMap::new(),
            jobs: JobRunner::new(),
            pending: BTreeMap::new(),
            meta: BTreeMap::new(),
            stages: Arc::new(Stages::default()),
            owner: "cli".into(),
            generation: 0,
            link: None,
            workspace: None,
            auto_store: false,
            base_library: Arc::new(femcore::material::Library::default()),
            library_notes: vec![],
        }
    }

    pub fn registry(&self) -> &Registry {
        &self.registry
    }

    /// Replaces the mesher (a front-end with its own mesher, or a test that counts or holds its calls). The mesh the engine holds is kept; the
    /// next one is built with the new mesher.
    pub fn set_mesher(&mut self, mesher: Arc<dyn Mesher>) {
        self.mesher = mesher;
    }

    pub fn project(&self) -> &Project {
        &self.project
    }

    pub fn jobs(&self) -> &JobRunner {
        &self.jobs
    }

    /// Names the workspace folder workflows are stored in. With `auto_store` the project is bound to a new workflow there at the first solve and
    /// every solve is recorded in its run table (origin `ai`).
    pub fn set_workspace(&mut self, root: PathBuf, auto_store: bool) {
        self.workspace = Some(root);
        self.auto_store = auto_store;
    }

    /// What the physics of the project want the mesh to resolve, merged: the shortest wavelength of each region, the highest order and the
    /// largest headroom. `None` when the mesh is not sized from the physics (`physics_sizing` off) or no physics states a wave.
    pub fn sizing(&self) -> Option<Sizing> {
        if !self.project.mesh.physics_sizing {
            return None;
        }
        let lib = self.library().ok();
        let mut out = Sizing { order: 1, headroom: 1.0, wavelengths: BTreeMap::new() };
        for spec in &self.project.physics {
            let Some(p) = self.registry.get(&spec.physics) else { continue };
            let s = p.sizing(spec, lib.as_ref());
            out.order = out.order.max(s.order);
            out.headroom = out.headroom.max(s.headroom);
            for (region, lam) in s.wavelengths {
                let e = out.wavelengths.entry(region).or_insert(lam);
                *e = e.min(lam);
            }
        }
        (!out.wavelengths.is_empty()).then_some(out)
    }

    /// The hash of the geometry, the mesh settings and the sizing: what a mesh depends on.
    pub fn mesh_key(&self) -> String {
        let v = json!({ "geometry": self.project.geometry, "mesh": self.project.mesh, "sizing": self.sizing() });
        hash_of(&nyx_physics::canonical_json(&v).to_string())
    }

    /// The hash of the whole project. Staleness no longer uses it (`problem_key` does); it is kept for tools that want an identity of the whole
    /// project in this process (a cache key, ADR 64).
    pub fn stamp(&self) -> String {
        hash_of(&self.project.to_json())
    }

    /// What a result of `physics` depends on, as a cache key: the problem of that physics, the materials the project carries, the mesh key and the
    /// problem keys of the physics it takes fields from. Settings that do not change a result (the active physics, the title, plot settings) are
    /// not in it. `None` when the project has no problem for `physics`. Not an identity: the physics a solve was *chained* with is not part of it
    /// (`solve_run` with one physics runs without the fields of the others), see ADR 84.
    pub fn problem_key(&self, physics: &str) -> Option<String> {
        self.project.spec(physics)?;
        Some(self.key_of(physics, &self.mesh_key(), &mut BTreeMap::new(), &mut vec![]))
    }

    /// The problem key of every physics of the project, with the mesh key computed once.
    pub fn problem_keys(&self) -> BTreeMap<String, String> {
        let mesh = self.mesh_key();
        let mut memo = BTreeMap::new();
        for s in &self.project.physics {
            self.key_of(&s.physics, &mesh, &mut memo, &mut vec![]);
        }
        memo
    }

    fn key_of(&self, id: &str, mesh: &str, memo: &mut BTreeMap<String, String>, stack: &mut Vec<String>) -> String {
        if let Some(k) = memo.get(id) {
            return k.clone();
        }
        let (Some(spec), Some(me)) = (self.project.spec(id), self.registry.get(id)) else { return String::new() };
        if stack.iter().any(|s| s == id) {
            return String::new(); // a cycle: solving refuses it, the key only has to terminate
        }
        stack.push(id.to_string());
        let mut from = std::collections::BTreeSet::new();
        for imp in me.imports() {
            for other in &self.project.physics {
                if other.physics != id
                    && self.registry.get(&other.physics).is_some_and(|q| q.exports().iter().any(|e| e.role == imp.role))
                {
                    from.insert(other.physics.clone());
                }
            }
        }
        let imports: BTreeMap<String, String> =
            from.into_iter().map(|o| (o.clone(), self.key_of(&o, mesh, memo, stack))).collect();
        stack.pop();
        let v = json!({ "physics": id, "problem": spec, "materials": self.project.materials, "mesh": mesh, "imports": imports });
        let key = hash_of(&nyx_physics::canonical_json(&v).to_string());
        memo.insert(id.to_string(), key.clone());
        key
    }

    /// True when the result stored for `physics` was not made from the problem the project has now (or the project has no such problem any more).
    pub fn is_stale(&self, physics: &str) -> bool {
        match (self.results.get(physics), self.problem_key(physics)) {
            (Some(s), Some(k)) => s.key != k,
            _ => true,
        }
    }

    /// Seconds since the result of `physics` was stored.
    pub fn result_age(&self, physics: &str) -> Option<f64> {
        self.results.get(physics).map(|s| s.at.elapsed().as_secs_f64())
    }

    /// Replaces the project (a new one, or one just opened). The mesh and results are dropped and the store binding is released.
    pub fn set_project(&mut self, project: Project) {
        self.project = project;
        self.generation += 1;
        self.model = None;
        self.results.clear();
        self.link = None;
        self.stages.clear();
    }

    /// The mesh if one has been built for the current geometry and settings.
    pub fn model(&self) -> Option<&Arc<MeshedModel>> {
        self.model.as_ref().filter(|(k, _)| *k == self.mesh_key()).map(|(_, m)| m)
    }

    /// The latest result of a physics and whether what it depends on has changed since (its problem key, not the whole project).
    pub fn result(&self, physics: &str) -> Option<(&Arc<SolveResult>, bool)> {
        let s = self.results.get(physics)?;
        Some((&s.result, self.is_stale(physics)))
    }

    // ----------------------------------------------------------------------------------------------------------------- jobs

    /// Starts meshing the project on a job (queued when the heavy lane is busy).
    pub fn start_mesh(&mut self) -> Result<JobId, Error> {
        self.start_mesh_with(&QueueOpts::default())
    }

    /// `start_mesh` with queue options.
    pub fn start_mesh_with(&mut self, opts: &QueueOpts) -> Result<JobId, Error> {
        self.project.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        let key = self.mesh_key();
        let (mesher, geometry, settings, sizing) =
            (self.mesher.clone(), self.project.geometry.clone(), self.project.mesh.clone(), self.sizing());
        let stages = self.stages.clone();
        let mesh_key = key.clone();
        let submit = self.submission("mesh", Lane::Heavy, opts, "mesh");
        let id = self.jobs.submit(submit, move |progress| {
            if progress.is_cancelled() {
                return Err(Error::Cancelled);
            }
            // the mesher checks the cancel flag between its passes and logs its stages, as a solve does
            let _installed = femcore::progress::install(&progress);
            if let Some(m) = stages.mesh(&mesh_key) {
                let _stage = femcore::progress::stage("mesh: the same mesh was built by an earlier request");
                return Ok(JobOutput::Mesh(m));
            }
            let m = Arc::new(mesher.mesh(&geometry, &settings, sizing.as_ref())?);
            stages.put_mesh(&mesh_key, &m);
            Ok(JobOutput::Mesh(m))
        });
        self.pending.insert(id, Pending::Mesh { key, generation: self.generation });
        self.note_job(id, "mesh", vec![]);
        Ok(id)
    }

    /// Starts solving `targets` (default: the active physics, else all) on a job. The mesh is built inside the job when none is current.
    pub fn start_solve(&mut self, targets: Option<Vec<String>>) -> Result<JobId, Error> {
        self.start_solve_with(targets, &QueueOpts::default())
    }

    /// `start_solve` with queue options. The job waits in the queue while the heavy lane is busy; it carries a frozen copy of the project.
    pub fn start_solve_with(&mut self, targets: Option<Vec<String>>, opts: &QueueOpts) -> Result<JobId, Error> {
        if self.project.physics.is_empty() {
            return Err(Error::invalid("the project has no physics; add one first (add_physics)"));
        }
        self.project.check(&self.registry)?;
        if self.auto_store && self.link.is_none() && self.workspace.is_some() {
            self.save_to_workspace(None)?;
        }
        let targets = targets.unwrap_or_else(|| self.default_targets());
        for t in &targets {
            if self.project.spec(t).is_none() {
                return Err(Error::invalid(format!("the project has no problem for '{t}'")));
            }
        }
        let (key, keys) = (self.mesh_key(), self.problem_keys());
        let (registry, mesher, project, sizing) =
            (self.registry.clone(), self.mesher.clone(), self.project.clone(), self.sizing());
        let library = self.library()?;
        let cached = self.model().cloned();
        let stages = self.stages.clone();
        let (mesh_key, solve_keys) = (key.clone(), keys.clone());
        let submit = self.submission(
            &format!("solve {}", targets.join(", ")),
            Lane::Heavy,
            opts,
            &format!("solve {}", targets.join(",")),
        );
        let chain = targets.clone();
        let id = self.jobs.submit(submit, {
            let project = project.clone();
            move |progress| {
                // the mesh: the one the engine held when the request was made, else one an earlier request built, else built here
                let model = match cached {
                    Some(m) => m,
                    None => {
                        let _stage = femcore::progress::install(&progress);
                        match stages.mesh(&mesh_key) {
                            Some(m) => {
                                let _s =
                                    femcore::progress::stage("mesh: the same mesh was built by an earlier request");
                                m
                            }
                            None => {
                                let m = Arc::new(mesher.mesh(&project.geometry, &project.mesh, sizing.as_ref())?);
                                stages.put_mesh(&mesh_key, &m);
                                m
                            }
                        }
                    }
                };
                stages.put_mesh(&mesh_key, &model);
                let (results, failed) =
                    crate::solve::run_chain_partial(&registry, &project, &model, &chain, &library, progress);
                match failed {
                    Some(Error::Cancelled) => return Err(Error::Cancelled),
                    Some(e) if results.is_empty() => return Err(e),
                    _ => {}
                }
                let cells = model.mesh().n_cells();
                for r in &results {
                    if let Some(k) = solve_keys.get(&r.physics) {
                        stages.put_solve(&solve_key(k, &chain), &Arc::new(SolveSummary::of(r, cells)));
                    }
                }
                Ok(JobOutput::Solved { model, results, failed })
            }
        });
        self.pending.insert(id, Pending::Solve { key, keys, project: Box::new(project), generation: self.generation });
        self.note_job(id, "solve", targets);
        Ok(id)
    }

    /// Starts a mesh convergence check (two solves, the second on a finer mesh) on a job. Both solves run on the job's thread under its progress,
    /// as two stages, so the check shows in `job_list`, reports where it is and can be cancelled.
    pub fn start_convergence(&mut self, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<JobId, Error> {
        self.start_convergence_with(physics, factor, tolerance, &QueueOpts::default())
    }

    /// `start_convergence` with queue options. When the solve of `physics` alone on the present mesh was done already by an earlier request, the
    /// check uses its numbers and solves only the finer mesh.
    pub fn start_convergence_with(
        &mut self,
        physics: &str,
        factor: f64,
        tolerance: Option<f64>,
        opts: &QueueOpts,
    ) -> Result<JobId, Error> {
        if !(factor > 0.0 && factor < 1.0) {
            return Err(Error::invalid("`factor` is the ratio of the finer element size to the present one, between 0 and 1 (0.6 is a good step)"));
        }
        if self.project.spec(physics).is_none() {
            return Err(Error::invalid(format!("the project has no problem for '{physics}'")));
        }
        let mut job = self.convergence_job();
        // the coarse solve of the check is the solve of this physics alone on the present mesh: its key is the problem key and that chain
        job.reuse = self.problem_key(physics).map(|k| (self.stages.clone(), solve_key(&k, &[physics.to_string()])));
        let physics = physics.to_string();
        let submit =
            self.submission(&format!("convergence {physics}"), Lane::Heavy, opts, &format!("convergence {physics}"));
        let name = physics.clone();
        let id = self.jobs.submit(submit, move |progress| {
            job.run_with(&progress, &physics, factor, tolerance).map(JobOutput::Convergence)
        });
        self.pending.insert(id, Pending::Convergence);
        self.note_job(id, "convergence", vec![name]);
        Ok(id)
    }

    /// Starts reading a results file, meshing its project again and checking the mesh against the stored one, on a job. The engine is not touched
    /// until the job is applied; then the file's project replaces the project and its results become the latest results.
    pub fn start_open_results(&mut self, path: &str) -> Result<JobId, Error> {
        let (registry, mesher, library, path) =
            (self.registry.clone(), self.mesher.clone(), self.base_library.clone(), path.to_string());
        let submit = self.submission("open results", Lane::Heavy, &QueueOpts::default(), "open results");
        let id = self.jobs.submit(submit, move |progress| {
            crate::results_file::read_and_mesh(&registry, &mesher, &library, &path, &progress)
                .map(|o| JobOutput::Opened(Box::new(o)))
        });
        self.pending.insert(id, Pending::Open);
        self.note_job(id, "open", vec![]);
        Ok(id)
    }

    /// Starts writing the results the engine holds now (the project, the mesh and the latest result of every physics) into an HDF5 file on a job.
    /// What is written is the snapshot of this moment, whatever the engine does next. The file appears under its name only when it is complete.
    pub fn start_export_results(&mut self, path: &str, parallel: bool) -> Result<JobId, Error> {
        let snapshot = self.results_snapshot()?;
        let path = path.to_string();
        // writing a file is the light lane: it runs beside a solve and never waits for one
        let submit = self.submission("write results", Lane::Light, &QueueOpts::default(), "write results");
        let id = self.jobs.submit(submit, move |progress| {
            if progress.is_cancelled() {
                return Err(Error::Cancelled);
            }
            let _stage = femcore::progress::stage("writing the results file");
            crate::results_file::write_snapshot(&snapshot, &path, parallel, true).map(JobOutput::Exported)
        });
        self.pending.insert(id, Pending::Export);
        self.note_job(id, "export", vec![]);
        Ok(id)
    }

    /// Takes the finished jobs and applies them: the mesh and results become current, a solve is recorded in the run table when the project
    /// is bound to a workflow. Returns what each job came to. A front-end calls this on every repaint.
    pub fn poll(&mut self) -> Vec<(JobId, Result<Value, Error>)> {
        let ids: Vec<JobId> = self.pending.keys().copied().collect();
        let mut out = vec![];
        for id in ids {
            if let Some(outcome) = self.jobs.try_take(id) {
                out.push((id, self.apply(id, outcome)));
            }
        }
        out
    }

    /// Waits for a job and applies it.
    pub fn wait(&mut self, id: JobId) -> Result<Value, Error> {
        if !self.pending.contains_key(&id) {
            return Err(Error::invalid(format!(
                "job {id} is not waiting to be applied (it was applied already, or there is no such job)"
            )));
        }
        // a wait with no timeout on a job that cannot start would never return (and a host that answers one call at a time could never send the
        // call that frees it)
        if let Some(why) = self.jobs.blocked(id) {
            return Err(Error::invalid(format!("job {id} cannot start: {why}")));
        }
        let outcome = self.jobs.wait(id);
        self.apply(id, outcome)
    }

    pub(crate) fn apply(&mut self, id: JobId, outcome: Result<JobOutput, Error>) -> Result<Value, Error> {
        let pending = self.pending.remove(&id).expect("a pending job");
        let seconds = self.jobs.info(id).map(|i| i.seconds).unwrap_or(0.0);
        match (pending, outcome) {
            (Pending::Mesh { key, generation }, Ok(JobOutput::Mesh(m))) => {
                let info = m.info();
                if generation != self.generation {
                    return Err(Error::invalid(
                        "the project was replaced while the mesh was built; the mesh was dropped",
                    ));
                }
                let current = self.install_model(key, m);
                Ok(json!({ "job": id, "state": "done", "mesh": info, "seconds": seconds, "current": current }))
            }
            (Pending::Solve { generation, .. }, Ok(JobOutput::Solved { .. })) if generation != self.generation => {
                Err(Error::invalid("the project was replaced while this solve ran; its results were dropped"))
            }
            (Pending::Solve { key, keys, project, .. }, Ok(JobOutput::Solved { model, results, failed })) => {
                let mesh = model.solver_mesh();
                self.install_model(key, model);
                // a job that finishes after the project changed, or after a newer job of the same physics, must not replace a result that is
                // current: it replaces the stored one only when its key is the one the physics has now, or when the stored one is stale and
                // this job started later (jobs are numbered in start order)
                let now = self.problem_keys();
                let mut summary = vec![];
                for r in &results {
                    let current = keys.contains_key(&r.physics) && keys.get(&r.physics) == now.get(&r.physics);
                    let replace = current
                        || match self.results.get(&r.physics) {
                            None => true,
                            Some(old) => Some(&old.key) != now.get(&r.physics) && id > old.job,
                        };
                    if replace {
                        let key = keys.get(&r.physics).cloned().unwrap_or_default();
                        self.results.insert(
                            r.physics.clone(),
                            Stored { result: Arc::new(r.clone()), key, job: id, at: Instant::now() },
                        );
                    }
                    summary.push(json!({
                        "physics": r.physics, "study": r.study, "dofs": r.dofs, "seconds": r.seconds, "solutions": r.solutions.len(),
                        "scalars": r.scalars, "current": current, "stored": replace,
                    }));
                }
                let mut log = self.jobs.log(id);
                for r in &results {
                    log.extend(r.log.iter().map(|l| format!("{}: {l}", r.physics)));
                }
                if let Some(e) = failed {
                    // a physics of the chain failed after others had been solved: what was solved stays (it is a valid result of its own
                    // problem), the run is recorded as failed and the caller gets the error
                    log.push(format!("stopped: {e}"));
                    let _ = self.record(&project, RunOutcome::Failed(&e), seconds, &log, Some(mesh));
                    return Err(e);
                }
                let run = self.record(&project, RunOutcome::Solved(&results), seconds, &log, Some(mesh));
                Ok(json!({ "job": id, "state": "done", "results": summary, "run": run, "seconds": seconds }))
            }
            (Pending::Convergence, Ok(JobOutput::Convergence(v))) => Ok(v),
            (Pending::Export, Ok(JobOutput::Exported(v))) => Ok(v),
            (Pending::Open, Ok(JobOutput::Opened(o))) => {
                let Opened { project, model, results, path } = *o;
                // the file's project replaces the present one (that is what opening means); the results are stored under the keys of the new project
                self.set_project(project);
                let key = self.mesh_key();
                let cells = model.mesh().n_cells();
                self.model = Some((key, model));
                let keys = self.problem_keys();
                let mut physics = vec![];
                for (code, r) in results {
                    physics.push(json!({
                        "physics": code, "study": r.study, "dofs": r.dofs, "seconds": r.seconds, "solutions": r.solutions.len(),
                        "scalars": r.scalars, "current": true, "stored": true,
                    }));
                    let k = keys.get(&code).cloned().unwrap_or_default();
                    self.results.insert(code, Stored { result: Arc::new(r), key: k, job: id, at: Instant::now() });
                }
                Ok(json!({ "title": self.project.title, "path": path, "cells": cells, "results": physics }))
            }
            (Pending::Solve { project, .. }, Err(e)) => {
                let log = self.jobs.log(id);
                // a failed or cancelled run is a record too: the table shows what was tried. A request cancelled before it started (replaced by a
                // newer one, or cancelled in the queue) never ran, so it leaves no record
                let never_ran = matches!(e, Error::Cancelled) && log.is_empty() && seconds == 0.0;
                if !never_ran {
                    let mesh = self.model.as_ref().map(|(_, m)| m.solver_mesh());
                    let _ = self.record(&project, RunOutcome::Failed(&e), seconds, &log, mesh);
                }
                Err(e)
            }
            (_, Err(e)) => Err(e),
            _ => Err(Error::invalid("a job returned something its kind does not produce")),
        }
    }

    /// Makes `model` (built for mesh key `key`) the engine's mesh unless a mesh that is current for the project is already there. True when it
    /// is the mesh the project has now.
    fn install_model(&mut self, key: String, model: Arc<MeshedModel>) -> bool {
        let now = self.mesh_key();
        let have_current = self.model.as_ref().is_some_and(|(k, _)| *k == now);
        if key == now || !have_current {
            self.model = Some((key.clone(), model));
        }
        key == now
    }

    fn record(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<Arc<femcore::mesh::SimplexMesh>>,
    ) -> Value {
        match self.link.as_mut() {
            None => Value::Null,
            Some(l) => match l.record(project, outcome, seconds, log, mesh.as_deref()) {
                Ok(r) => json!({ "id": r.id, "version": r.version, "status": r.status }),
                Err(e) => json!({ "error": e.to_string() }),
            },
        }
    }

    /// Mesh now (start and wait).
    pub fn mesh_now(&mut self) -> Result<Value, Error> {
        let id = self.start_mesh()?;
        self.wait(id)
    }

    /// Solve now (start and wait).
    pub fn solve_now(&mut self, targets: Option<Vec<String>>) -> Result<Value, Error> {
        let id = self.start_solve(targets)?;
        self.wait(id)
    }

    // ----------------------------------------------------------------------------------------------------------------- store

    /// Binds the project to a new workflow folder in `workspace` (or commits it, when already bound).
    pub fn save_to_workspace(&mut self, workspace: Option<PathBuf>) -> Result<Value, Error> {
        if let Some(l) = self.link.as_mut() {
            let v = l.commit(&self.project)?;
            return Ok(
                json!({ "workflow": l.dir.wf.id, "experiment": l.experiment, "version": v.unwrap_or(l.dir.wf.version), "changed": v.is_some() }),
            );
        }
        let root = workspace
            .or_else(|| self.workspace.clone())
            .ok_or_else(|| Error::invalid("name a workspace folder (`workspace`) to save into"))?;
        std::fs::create_dir_all(&root).map_err(|e| Error::invalid(format!("{}: {e}", root.display())))?;
        let link = StoreLink::create(&root, &self.project, "ai")?;
        let v = json!({ "workflow": link.dir.wf.id, "folder": link.dir.root, "experiment": link.experiment, "version": link.dir.wf.version, "changed": true });
        self.workspace = Some(root);
        self.link = Some(link);
        Ok(v)
    }

    /// Opens a workflow folder and makes its experiment the project.
    pub fn open_workflow(&mut self, folder: &std::path::Path, experiment: Option<&str>) -> Result<Value, Error> {
        let (link, project) = StoreLink::open(folder, experiment)?;
        project.check(&self.registry)?;
        let v = json!({ "workflow": link.dir.wf.id, "title": link.dir.wf.title, "experiment": link.experiment, "version": link.dir.wf.version });
        self.set_project(project);
        self.link = Some(link);
        Ok(v)
    }
}

impl Default for Engine {
    fn default() -> Self {
        Engine::new()
    }
}
