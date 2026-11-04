//! The engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them and the binding to a
//! workflow folder. The window and the MCP server hold an `Engine` and nothing else; every action is a method here or a named command
//! (`commands.rs`, which calls these methods).
//!
//! The mesh and the results are caches keyed by the content they were made from (a hash of the geometry and mesh settings, a hash of the
//! whole project), so an edit never has to remember to invalidate them: a result whose stamp is not the project's is *stale* and says so.

use crate::jobs::{JobId, JobOutput, JobRunner};
use crate::meshing::{FemcoreMesher, MeshedModel, Mesher};
use crate::project::Project;
use crate::store_link::{RunOutcome, StoreLink};
use nyx_physics::{Error, Registry, SolveResult};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

pub(crate) struct Stored {
    pub result: Arc<SolveResult>,
    /// the hash of the project it was solved from
    pub stamp: String,
}

enum Pending {
    Mesh { key: String },
    Solve { key: String, stamp: String, project: Box<Project> },
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
}

fn hash_of(text: &str) -> String {
    nyx_store::util::hash_bytes(text.as_bytes())
}

impl Engine {
    /// An engine with the physics compiled into this build and an empty 2D project.
    pub fn new() -> Engine {
        Engine::with(crate::registry(), Arc::new(FemcoreMesher))
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

    /// The hash of the geometry and the mesh settings: what a mesh depends on.
    pub fn mesh_key(&self) -> String {
        let v = json!({ "geometry": self.project.geometry, "mesh": self.project.mesh });
        hash_of(&nyx_physics::canonical_json(&v).to_string())
    }

    /// The hash of the whole project: what a result depends on.
    pub fn stamp(&self) -> String {
        hash_of(&self.project.to_json())
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

    /// The latest result of a physics and whether the project has changed since.
    pub fn result(&self, physics: &str) -> Option<(&Arc<SolveResult>, bool)> {
        let s = self.results.get(physics)?;
        Some((&s.result, s.stamp != self.stamp()))
    }

    // ----------------------------------------------------------------------------------------------------------------- jobs

    /// Starts meshing the project on a job.
    pub fn start_mesh(&mut self) -> Result<JobId, Error> {
        self.project.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        let key = self.mesh_key();
        let (mesher, geometry, settings) = (self.mesher.clone(), self.project.geometry.clone(), self.project.mesh.clone());
        let id = self.jobs.start("mesh", move |progress| {
            if progress.is_cancelled() {
                return Err(Error::Cancelled);
            }
            Ok(JobOutput::Mesh(Arc::new(mesher.mesh(&geometry, &settings)?)))
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
        let targets = targets.unwrap_or_else(|| self.project.default_targets());
        for t in &targets {
            if self.project.spec(t).is_none() {
                return Err(Error::invalid(format!("the project has no problem for '{t}'")));
            }
        }
        let (key, stamp) = (self.mesh_key(), self.stamp());
        let (registry, mesher, project) = (self.registry.clone(), self.mesher.clone(), self.project.clone());
        let cached = self.model().cloned();
        let id = self.jobs.start(&format!("solve {}", targets.join(", ")), {
            let project = project.clone();
            move |progress| {
                let model = match cached {
                    Some(m) => m,
                    None => {
                        let _stage = femcore::progress::install(&progress);
                        Arc::new(mesher.mesh(&project.geometry, &project.mesh)?)
                    }
                };
                let results = crate::solve::run_chain(&registry, &project, &model, &targets, progress)?;
                Ok(JobOutput::Solved { model, results })
            }
        });
        self.pending.insert(id, Pending::Solve { key, stamp, project: Box::new(project) });
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
            return Err(Error::invalid(format!("job {id} is not waiting to be applied (it was applied already, or there is no such job)")));
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
                self.model = Some((key, m));
                Ok(json!({ "job": id, "state": "done", "mesh": info, "seconds": seconds }))
            }
            (Pending::Solve { key, stamp, project }, Ok(JobOutput::Solved { model, results })) => {
                self.model = Some((key, model));
                let mut summary = vec![];
                for r in &results {
                    summary.push(json!({
                        "physics": r.physics, "study": r.study, "dofs": r.dofs, "seconds": r.seconds, "solutions": r.solutions.len(),
                        "scalars": r.scalars,
                    }));
                    self.results.insert(r.physics.clone(), Stored { result: Arc::new(r.clone()), stamp: stamp.clone() });
                }
                let mut log = self.jobs.log(id);
                for r in &results {
                    log.extend(r.log.iter().map(|l| format!("{}: {l}", r.physics)));
                }
                let run = self.record(&project, RunOutcome::Solved(&results), seconds, &log);
                Ok(json!({ "job": id, "state": "done", "results": summary, "run": run, "seconds": seconds }))
            }
            (Pending::Solve { project, .. }, Err(e)) => {
                let log = self.jobs.log(id);
                // a failed or cancelled run is a record too: the table shows what was tried
                let _ = self.record(&project, RunOutcome::Failed(&e), seconds, &log);
                Err(e)
            }
            (_, Err(e)) => Err(e),
            _ => Err(Error::invalid("a job returned something its kind does not produce")),
        }
    }

    fn record(&mut self, project: &Project, outcome: RunOutcome, seconds: f64, log: &[String]) -> Value {
        match self.link.as_mut() {
            None => Value::Null,
            Some(l) => match l.record(project, outcome, seconds, log) {
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
            return Ok(json!({ "workflow": l.dir.wf.id, "experiment": l.experiment, "version": v.unwrap_or(l.dir.wf.version), "changed": v.is_some() }));
        }
        let root = workspace.or_else(|| self.workspace.clone()).ok_or_else(|| Error::invalid("name a workspace folder (`workspace`) to save into"))?;
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
