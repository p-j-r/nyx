//! The engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them and the binding to a
//! workflow folder. The window and the MCP server hold an `Engine` and nothing else; every action is a method here or a named command
//! (`commands.rs`, which calls these methods).
//!
//! The mesh and the results are caches keyed by the content they were made from (a hash of the geometry and mesh settings; for a result the
//! *problem key* of its physics: its own problem, the project materials, the mesh key and the keys of the physics it takes fields from), so an edit
//! never has to remember to invalidate them: a result whose key is not the one its physics has now is *stale* and says so. A key is a cache key
//! (ADR 64, ADR 84): it is never stored in a file or sent to a front-end.

use crate::jobs::{JobId, JobOutput, JobRunner};
use crate::meshing::{FemcoreMesher, MeshedModel, Mesher};
use crate::project::Project;
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

enum Pending {
    Mesh { key: String },
    Solve { key: String, keys: BTreeMap<String, String>, project: Box<Project> },
}

pub struct Engine {
    pub(crate) registry: Registry,
    pub(crate) mesher: Arc<dyn Mesher>,
    pub(crate) project: Project,
    pub(crate) model: Option<(String, Arc<MeshedModel>)>,
    pub(crate) results: BTreeMap<String, Stored>,
    pub(crate) jobs: JobRunner,
    pending: BTreeMap<JobId, Pending>,
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
        self.model = None;
        self.results.clear();
        self.link = None;
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

    /// Starts meshing the project on a job.
    pub fn start_mesh(&mut self) -> Result<JobId, Error> {
        self.project.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        let key = self.mesh_key();
        let (mesher, geometry, settings, sizing) =
            (self.mesher.clone(), self.project.geometry.clone(), self.project.mesh.clone(), self.sizing());
        let id = self.jobs.start("mesh", move |progress| {
            if progress.is_cancelled() {
                return Err(Error::Cancelled);
            }
            Ok(JobOutput::Mesh(Arc::new(mesher.mesh(&geometry, &settings, sizing.as_ref())?)))
        });
        self.pending.insert(id, Pending::Mesh { key });
        Ok(id)
    }

    /// Starts solving `targets` (default: the active physics, else all) on a job. The mesh is built inside the job when none is current.
    pub fn start_solve(&mut self, targets: Option<Vec<String>>) -> Result<JobId, Error> {
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
        let id = self.jobs.start(&format!("solve {}", targets.join(", ")), {
            let project = project.clone();
            move |progress| {
                let model = match cached {
                    Some(m) => m,
                    None => {
                        let _stage = femcore::progress::install(&progress);
                        Arc::new(mesher.mesh(&project.geometry, &project.mesh, sizing.as_ref())?)
                    }
                };
                let results = crate::solve::run_chain(&registry, &project, &model, &targets, &library, progress)?;
                Ok(JobOutput::Solved { model, results })
            }
        });
        self.pending.insert(id, Pending::Solve { key, keys, project: Box::new(project) });
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
        let outcome = self.jobs.wait(id);
        self.apply(id, outcome)
    }

    fn apply(&mut self, id: JobId, outcome: Result<JobOutput, Error>) -> Result<Value, Error> {
        let pending = self.pending.remove(&id).expect("a pending job");
        let seconds = self.jobs.info(id).map(|i| i.seconds).unwrap_or(0.0);
        match (pending, outcome) {
            (Pending::Mesh { key }, Ok(JobOutput::Mesh(m))) => {
                let info = m.info();
                let current = self.install_model(key, m);
                Ok(json!({ "job": id, "state": "done", "mesh": info, "seconds": seconds, "current": current }))
            }
            (Pending::Solve { key, keys, project }, Ok(JobOutput::Solved { model, results })) => {
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
                let run = self.record(&project, RunOutcome::Solved(&results), seconds, &log, Some(mesh));
                Ok(json!({ "job": id, "state": "done", "results": summary, "run": run, "seconds": seconds }))
            }
            (Pending::Solve { project, .. }, Err(e)) => {
                let log = self.jobs.log(id);
                // a failed or cancelled run is a record too: the table shows what was tried
                let mesh = self.model.as_ref().map(|(_, m)| m.solver_mesh());
                let _ = self.record(&project, RunOutcome::Failed(&e), seconds, &log, mesh);
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
