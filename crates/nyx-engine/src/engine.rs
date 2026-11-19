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
    /// the mesh it was solved on: the reductions (ADR 90) read its fields on exactly these cells, whatever the project's mesh is now
    pub model: Arc<MeshedModel>,
    /// the run record that holds this result's solve, when the engine wrote one (the verdict of a convergence check goes into it, P6)
    pub run: Option<String>,
}

/// What a solve job is for: the engine's own project, or one point of a sweep.
#[derive(Clone, Copy)]
pub(crate) enum Purpose {
    Current,
    Point { sweep: u64, index: usize },
}

pub(crate) enum Pending {
    Mesh {
        key: String,
        generation: u64,
    },
    Solve {
        key: String,
        keys: BTreeMap<String, String>,
        project: Box<Project>,
        generation: u64,
    },
    /// a point of a sweep: its frozen project, solved for the table and recorded; it never replaces the engine's own results
    Point {
        sweep: u64,
        index: usize,
        project: Box<Project>,
    },
    /// the physics checked and the problem key it was checked for (a verdict belongs to a problem, not to a result, ADR 91)
    Convergence {
        physics: String,
        key: String,
    },
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
    /// the sweeps started on this engine, by number (`sweep.rs`)
    pub(crate) sweeps: BTreeMap<u64, crate::sweep::Sweep>,
    pub(crate) next_sweep: u64,
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
    /// the worker program a solve runs in when it asks to be isolated (`set_worker`), and whether it does when it does not say (P4, ADR 86)
    pub(crate) worker: Option<PathBuf>,
    pub(crate) isolate_default: bool,
    /// false for an engine whose registry or mesher is not the compiled-in one: a worker could not reproduce it
    pub(crate) isolable: bool,
    /// what solves of this machine have cost, learned from the workers (`estimate.rs`), and how the engine asks the machine about its memory
    pub(crate) history: Arc<crate::estimate::History>,
    pub(crate) memory_probe: crate::estimate::MemoryProbe,
    /// what the mesh-convergence checks of this session found, by the problem key they were made for (`validated.rs`)
    pub(crate) checks: BTreeMap<String, crate::validated::Check>,
    /// when set, every file or folder argument of a command must lie inside this folder (`path_guard.rs`); a script host turns it on
    pub(crate) path_root: Option<PathBuf>,
    /// the script (name, text) driving this engine, when one does: the runs it makes keep the text and its hash (`store_link.rs`)
    pub(crate) script: Option<(String, String)>,
}

/// What a worker reported, in words for the log of the run.
fn worker_words(w: &crate::jobs::WorkerRun) -> String {
    let gb = |b: u64| format!("{:.2} GB", b as f64 / 1e9);
    match (w.peak_working_set, w.peak_commit) {
        (Some(ws), Some(c)) => format!("peak memory {} in RAM, {} committed", gb(ws), gb(c)),
        (Some(ws), None) => format!("peak memory {} in RAM", gb(ws)),
        _ => "peak memory not reported by this system".into(),
    }
}

fn hash_of(text: &str) -> String {
    nyx_store::util::hash_bytes(text.as_bytes())
}

impl Engine {
    /// An engine with the physics compiled into this build and an empty 2D project.
    pub fn new() -> Engine {
        let mut e = Engine::with(crate::registry(), Arc::new(FemcoreMesher));
        e.isolable = true;
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
            sweeps: BTreeMap::new(),
            next_sweep: 0,
            generation: 0,
            link: None,
            workspace: None,
            auto_store: false,
            base_library: Arc::new(femcore::material::Library::default()),
            library_notes: vec![],
            worker: None,
            isolate_default: false,
            isolable: false,
            history: Arc::new(crate::estimate::History::default()),
            memory_probe: Arc::new(crate::process::system_memory),
            checks: BTreeMap::new(),
            path_root: None,
            script: None,
        }
    }

    /// Names the script that drives this engine (`name`, `text`): every run recorded while it is set keeps the text, with its hash in the file name,
    /// beside the results. A script host sets it for the run of the script and clears it after; `None` clears it.
    pub fn set_script(&mut self, script: Option<(String, String)>) {
        self.script = script;
    }

    /// Keeps every file and folder argument of the commands inside `root` (S4, ADR [89]): a relative path is taken from `root`, a path that
    /// leaves it (by `..`, an absolute path or a link) makes the command fail before it runs. `None` turns the guard off. A front-end that
    /// runs a script turns it on for the run; the window and the MCP server do not use it.
    pub fn set_path_root(&mut self, root: Option<&std::path::Path>) -> Result<(), String> {
        self.path_root = match root {
            Some(r) => Some(crate::path_guard::absolute_root(r)?),
            None => None,
        };
        Ok(())
    }

    /// The folder the path guard keeps commands in, when it is on.
    pub fn path_root(&self) -> Option<&std::path::Path> {
        self.path_root.as_deref()
    }

    /// Keeps what solves cost in `path` (JSON lines, appended to) and learns from what it holds: the memory and time estimates of later solves
    /// come from it. A front-end names a file in its workspace; without one the engine learns for the session only.
    pub fn set_history_file(&mut self, path: &std::path::Path) {
        self.history.set_file(path);
    }

    /// What the engine has learned about the cost of solves.
    pub fn history(&self) -> &Arc<crate::estimate::History> {
        &self.history
    }

    /// Replaces how the engine asks the machine for `(total, available)` memory in bytes (a test that plays a machine that is full).
    pub fn set_memory_probe(&mut self, probe: crate::estimate::MemoryProbe) {
        self.memory_probe = probe;
    }

    /// Names the worker program (`nyx-worker`, built from this crate) that solves run in when they are isolated: a solve in a child process
    /// cannot take this process down (out of memory) and can be stopped at once (killed), at the price of a results round trip. With
    /// `by_default` every solve of this engine runs there unless it says `isolate: false`; without it only a solve that says `isolate: true`
    /// does. The solves of an engine with another registry or mesher than the compiled-in ones cannot be isolated.
    pub fn set_worker(&mut self, worker: PathBuf, by_default: bool) {
        self.worker = Some(worker);
        self.isolate_default = by_default;
    }

    /// Makes the running program its own worker: a front-end (`nyx`, `nyx-mcp`) that answers `--job FOLDER` with `worker::serve` calls this, so
    /// there is nothing else to build or to ship beside it. Never call it from a test: the program would be the test runner. False when the
    /// system does not say what the running program is.
    pub fn set_self_as_worker(&mut self, by_default: bool) -> bool {
        match std::env::current_exe() {
            Ok(exe) => {
                self.set_worker(exe, by_default);
                true
            }
            Err(_) => false,
        }
    }

    pub fn registry(&self) -> &Registry {
        &self.registry
    }

    /// Replaces the mesher (a front-end with its own mesher, or a test that counts or holds its calls). The mesh the engine holds is kept; the
    /// next one is built with the new mesher.
    pub fn set_mesher(&mut self, mesher: Arc<dyn Mesher>) {
        self.isolable = false;
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
        self.submit_solve(None, targets, opts, Purpose::Current)
    }

    /// The one place a solve is submitted. `src` is the engine whose project is solved (a point of a sweep is a scratch engine with the point's
    /// edits applied; `None` is this engine). The job is always run and queued by this engine.
    pub(crate) fn submit_solve(
        &mut self,
        src: Option<&Engine>,
        targets: Option<Vec<String>>,
        opts: &QueueOpts,
        purpose: Purpose,
    ) -> Result<JobId, Error> {
        let lane_key = |t: &[String]| match purpose {
            Purpose::Current => format!("solve {}", t.join(",")),
            Purpose::Point { sweep, index } => format!("sweep {sweep} point {index}"),
        };
        {
            let s = src.unwrap_or(&*self);
            if s.project.physics.is_empty() {
                return Err(Error::invalid("the project has no physics; add one first (add_physics)"));
            }
            s.project.check(&s.registry)?;
        }
        if matches!(purpose, Purpose::Current) && self.auto_store && self.link.is_none() && self.workspace.is_some() {
            self.save_to_workspace(None)?;
        }
        let s = src.unwrap_or(&*self);
        let targets = targets.unwrap_or_else(|| s.default_targets());
        for t in &targets {
            if s.project.spec(t).is_none() {
                return Err(Error::invalid(format!("the project has no problem for '{t}'")));
            }
        }
        let (key, keys) = (s.mesh_key(), s.problem_keys());
        let (registry, mesher, project, sizing) = (s.registry.clone(), s.mesher.clone(), s.project.clone(), s.sizing());
        let library = s.library()?;
        // a point of a sweep has its own geometry and mesh settings: it never takes the mesh this engine holds
        let cached = if matches!(purpose, Purpose::Current) { s.model().cloned() } else { None };
        let stages = self.stages.clone();
        let (mesh_key, solve_keys) = (key.clone(), keys.clone());
        let label = match purpose {
            Purpose::Current => format!("solve {}", targets.join(", ")),
            Purpose::Point { sweep, index } => format!("sweep {sweep} point {index}: {}", targets.join(", ")),
        };
        let submit = self.submission(&label, Lane::Heavy, opts, &lane_key(&targets));
        let chain = targets.clone();
        // in a worker process (P4)? the library the worker is given is the one this job would have used
        // asked for in the request, it must be possible or the request fails with the reason; only by the engine's default, a solve that cannot
        // be isolated (no worker, another registry or mesher, no HDF5) runs in this process as before
        let explicit = opts.isolate == Some(true);
        let remote = if opts.isolate.unwrap_or(self.isolate_default) {
            let why_not = if self.worker.is_none() {
                Some("this engine has no worker program (`Engine::set_worker`: nyx and nyx-mcp use themselves)")
            } else if !self.isolable {
                Some("this engine has a registry or mesher of its own, which a worker process cannot rebuild")
            } else if cfg!(not(feature = "hdf5")) {
                Some("this build has no HDF5 support, which the worker needs to hand its results back")
            } else {
                None
            };
            match (why_not, self.worker.clone()) {
                (None, Some(worker)) => Some(crate::worker::RemoteSpec {
                    worker,
                    project_json: project.to_json(),
                    materials: crate::worker::material_texts(&library),
                    targets: targets.clone(),
                }),
                (Some(why), _) if explicit => return Err(Error::invalid(format!("`isolate` is not possible: {why}"))),
                _ => None,
            }
        } else {
            None
        };
        let estimate_slot: Arc<std::sync::Mutex<Option<crate::estimate::Estimate>>> = Default::default();
        let (history, probe, force, slot) =
            (self.history.clone(), self.memory_probe.clone(), opts.force, estimate_slot.clone());
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
                // what this solve is expected to cost, from the solves this machine has run; it starts only when the memory is there
                if let Some(est) = history.estimate(&crate::estimate::chain_key(&chain), model.mesh().dim() as u8, model.mesh().n_cells() as u64) {
                    *slot.lock().expect("estimate") = Some(est.clone());
                    crate::estimate::admit(&est, force, &probe, &progress)?;
                }
                let (results, failed, worker) = match &remote {
                    None => {
                        let (r, f) =
                            crate::solve::run_chain_partial(&registry, &project, &model, &chain, &library, progress);
                        (r, f, None)
                    }
                    Some(spec) => {
                        let done = crate::worker::run_remote(spec, &progress)?;
                        {
                            // the last word of the log: what the process cost (the log shown is the worker's, mirrored, until now)
                            let _installed = femcore::progress::install(&progress);
                            femcore::progress::note(format!("worker process: {}", worker_words(&done.run)));
                        }
                        if let Some(m) = &done.mesh {
                            if !crate::worker::same_mesh(m, &model.solver_mesh()) {
                                return Err(Error::Solver(
                                    "the worker meshed the project to something else than this process did, so its results cannot be shown on this mesh (a worker of another version of nyx?)".into(),
                                ));
                            }
                        }
                        // in the order they were solved, which the file does not keep
                        let order: Vec<String> = registry
                            .solve_order(&chain.iter().map(String::as_str).collect::<Vec<_>>())
                            .unwrap_or_else(|_| chain.clone());
                        let mut stored = done.results;
                        stored.sort_by_key(|(code, _)| order.iter().position(|o| o == code).unwrap_or(usize::MAX));
                        (stored.into_iter().map(|(_, r)| r).collect(), done.failed, Some(done.run))
                    }
                };
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
                Ok(JobOutput::Solved { model, results, failed, worker })
            }
        });
        match purpose {
            Purpose::Current => {
                self.pending
                    .insert(id, Pending::Solve { key, keys, project: Box::new(project), generation: self.generation });
                self.note_job_estimate(id, "solve", targets, estimate_slot);
            }
            Purpose::Point { sweep, index } => {
                self.pending.insert(id, Pending::Point { sweep, index, project: Box::new(project) });
                self.note_job_estimate(id, "sweep", targets, estimate_slot);
            }
        }
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
        let key = self.problem_key(&physics).unwrap_or_default();
        let submit =
            self.submission(&format!("convergence {physics}"), Lane::Heavy, opts, &format!("convergence {physics}"));
        let name = physics.clone();
        let id = self.jobs.submit(submit, move |progress| {
            job.run_with(&progress, &physics, factor, tolerance).map(JobOutput::Convergence)
        });
        self.pending.insert(id, Pending::Convergence { physics: name.clone(), key });
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

    /// Why a queued job cannot start however long it is waited for (it is held, or the queue is paused); `None` for a job that is running or
    /// will start, and for a finished one. A host that polls a job instead of blocking on it asks, so that it does not poll for ever.
    pub fn job_blocked(&self, id: JobId) -> Option<&'static str> {
        self.jobs.blocked(id)
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
            (Pending::Solve { key, keys, project, .. }, Ok(JobOutput::Solved { model, results, failed, worker })) => {
                let (cells, dim) = (model.mesh().n_cells() as u64, model.mesh().dim() as u8);
                let mesh = model.solver_mesh();
                let solved_on = model.clone();
                self.install_model(key, model);
                // a job that finishes after the project changed, or after a newer job of the same physics, must not replace a result that is
                // current: it replaces the stored one only when its key is the one the physics has now, or when the stored one is stale and
                // this job started later (jobs are numbered in start order)
                let now = self.problem_keys();
                let mut summary = vec![];
                // which results replace the stored ones (with the key they are stored under); they are moved into the store after the run is
                // recorded, not copied: a result holds every field of every solution, and a copy on the thread that draws cost 6 ms (plan P3)
                let mut keep: Vec<Option<String>> = vec![];
                for r in &results {
                    let current = keys.contains_key(&r.physics) && keys.get(&r.physics) == now.get(&r.physics);
                    let replace = current
                        || match self.results.get(&r.physics) {
                            None => true,
                            Some(old) => Some(&old.key) != now.get(&r.physics) && id > old.job,
                        };
                    keep.push(replace.then(|| keys.get(&r.physics).cloned().unwrap_or_default()));
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
                    self.store_results(results, keep, id, &solved_on, None);
                    return Err(e);
                }
                self.learn(&results, cells, dim, &worker);
                let run = self.record(&project, RunOutcome::Solved(&results), seconds, &log, Some(mesh));
                let run_id = run["id"].as_str().map(str::to_string);
                self.store_results(results, keep, id, &solved_on, run_id);
                let mut answer =
                    json!({ "job": id, "state": "done", "results": summary, "run": run, "seconds": seconds });
                if let Some(w) = worker {
                    answer["worker"] = json!(w);
                }
                Ok(answer)
            }
            (Pending::Point { sweep, index, project }, outcome) => {
                self.apply_point(id, sweep, index, *project, outcome)
            }
            (Pending::Convergence { physics, key }, Ok(JobOutput::Convergence(v))) => {
                self.checked_answer(&physics, &key, v)
            }
            (Pending::Export, Ok(JobOutput::Exported(v))) => Ok(v),
            (Pending::Open, Ok(JobOutput::Opened(o))) => {
                let Opened { project, model, results, path } = *o;
                // the file's project replaces the present one (that is what opening means); the results are stored under the keys of the new project
                self.set_project(project);
                let key = self.mesh_key();
                let cells = model.mesh().n_cells();
                self.model = Some((key, model.clone()));
                let keys = self.problem_keys();
                let mut physics = vec![];
                for (code, r) in results {
                    physics.push(json!({
                        "physics": code, "study": r.study, "dofs": r.dofs, "seconds": r.seconds, "solutions": r.solutions.len(),
                        "scalars": r.scalars, "current": true, "stored": true,
                    }));
                    let k = keys.get(&code).cloned().unwrap_or_default();
                    self.results.insert(
                        code,
                        Stored {
                            result: Arc::new(r),
                            key: k,
                            job: id,
                            at: Instant::now(),
                            model: model.clone(),
                            run: None,
                        },
                    );
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

    /// Moves the results of job `id` that replace the stored ones (`keep` holds, per result, the key it is stored under) into the store, each
    /// with the run that recorded it. Nothing is copied.
    fn store_results(
        &mut self,
        results: Vec<SolveResult>,
        keep: Vec<Option<String>>,
        id: JobId,
        solved_on: &Arc<MeshedModel>,
        run: Option<String>,
    ) {
        for (r, key) in results.into_iter().zip(keep) {
            let Some(key) = key else { continue };
            self.results.insert(
                r.physics.clone(),
                Stored {
                    result: Arc::new(r),
                    key,
                    job: id,
                    at: Instant::now(),
                    model: solved_on.clone(),
                    run: run.clone(),
                },
            );
        }
    }

    /// Learns what a finished worker run cost: its peak memory against the size of the problem, for the estimates of later solves. Only runs
    /// whose chain finished and whose worker measured its memory count (an in-process run cannot say how much of the process was its own).
    pub(crate) fn learn(&self, results: &[SolveResult], cells: u64, dim: u8, worker: &Option<crate::jobs::WorkerRun>) {
        let Some(ws) = worker.as_ref().and_then(|w| w.peak_working_set) else { return };
        if results.is_empty() {
            return;
        }
        let chain: Vec<String> = results.iter().map(|r| r.physics.clone()).collect();
        self.history.add(crate::estimate::Row {
            chain: crate::estimate::chain_key(&chain),
            dim,
            cells,
            dofs: results.iter().map(|r| r.dofs as u64).max().unwrap_or(0),
            peak_working_set: ws,
            peak_commit: worker.as_ref().and_then(|w| w.peak_commit).unwrap_or(0),
            seconds: results.iter().map(|r| r.seconds).sum(),
        });
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
        let script = self.script.clone();
        match self.link.as_mut() {
            None => Value::Null,
            Some(l) => {
                l.script = script;
                match l.record(project, outcome, seconds, log, mesh.as_deref()) {
                    Ok(r) => json!({ "id": r.id, "version": r.version, "status": r.status }),
                    Err(e) => json!({ "error": e.to_string() }),
                }
            }
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
