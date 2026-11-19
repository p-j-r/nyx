//! What a person sees while editing, before anything is built for real (productivity plan P3, ADR 92). Three questions about an edit that has
//! not been applied yet, answered without changing the engine, its mesh, its results or its queue:
//!
//! * `check_edits`: does the project with these geometry and mesh edits pass `Project::check` (the feature list, the element size, every
//!   problem against its physics)? Microseconds.
//! * `size_preview`: how many cells would these edits mesh to, and what does a solve of that size cost on this machine? A trial mesh of a copy
//!   of the edits, run by the caller on any thread (`SizeTask::run`), plus the estimate of `estimate.rs` (memory and seconds from the solves this
//!   machine has run: no number when there are too few). Only for models whose mesh is cheap: 2D, and 3D built by sweeping a 2D section.
//! * `polarization_task`: the polarization analysis of a result as a function that can run off the thread that draws (`nyx-analysis` is a
//!   pure function of the mesh and the result, so a front-end can compute it beside the frame).
//!
//! Nothing here starts a job, takes the heavy lane, stores a mesh in the stage cache or replaces `Engine::model`: a trial mesh is a drawing
//! of numbers, never the mesh a solve uses. "No mesh or solve without a request" (plan P3) is kept: the caller asks, the answer is a value.

use crate::engine::Engine;
use crate::estimate::{chain_key, Estimate, History};
use crate::meshing::Mesher;
use crate::project::MeshSettings;
use nyx_analysis::ModeAnalysis;
use nyx_physics::{Error, Sizing, SolveResult};
use serde::Serialize;
use std::sync::Arc;
use std::time::Instant;

use geom::model::FeatureList;

/// The size and cost of a model that is being edited.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SizePreview {
    /// dimension of the mesh a solve would use (a 2D section swept along y is 3)
    pub dim: u8,
    pub cells: u64,
    /// seconds the trial mesh took on this machine: a front-end that sees this above its budget stops previewing until asked
    pub mesh_seconds: f64,
    /// the physics a solve of the project would solve together (`default_targets`), the key the estimate is looked up under
    pub chain: String,
    /// memory and time of such a solve from the runs this machine has done, with the number of runs and the largest error of the fit; `None`
    /// (and no number invented) when there are fewer than two runs of different sizes
    pub estimate: Option<Estimate>,
}

/// A trial mesh waiting to run. It owns copies of everything it needs and is `Send`, so the caller chooses the thread.
pub struct SizeTask {
    mesher: Arc<dyn Mesher>,
    geometry: FeatureList,
    settings: MeshSettings,
    sizing: Option<Sizing>,
    history: Arc<History>,
    chain: String,
    sweep_dim: u8,
}

impl SizeTask {
    /// Meshes the copy of the edits and looks up the estimate. Blocks for as long as the mesher takes (milliseconds for the models this is
    /// offered for).
    pub fn run(self) -> Result<SizePreview, Error> {
        let t = Instant::now();
        let model = self.mesher.mesh(&self.geometry, &self.settings, self.sizing.as_ref())?;
        let mesh_seconds = t.elapsed().as_secs_f64();
        let cells = model.mesh().n_cells() as u64;
        let dim = model.mesh().dim() as u8;
        debug_assert!(dim == self.sweep_dim, "the dimension the task expected is the one the mesher made");
        let estimate = self.history.estimate(&self.chain, dim, cells);
        Ok(SizePreview { dim, cells, mesh_seconds, chain: self.chain, estimate })
    }
}

impl Engine {
    /// The first problem of this project with `geometry` and `mesh` in place of its own, in plain words, or `Ok` when it passes
    /// `Project::check`. The engine's project is not changed.
    pub fn check_edits(&self, geometry: &FeatureList, mesh: &MeshSettings) -> Result<(), Error> {
        let mut p = self.project.clone();
        p.geometry = geometry.clone();
        p.mesh = mesh.clone();
        p.check(&self.registry)
    }

    /// A trial mesh of the project with `geometry` and `mesh` applied, as a task that can run on any thread. Refused (with the reason, in plain
    /// words) when the edits do not pass `check_edits`, or when the mesh is not known to be cheap: only 2D models and 3D models swept from a
    /// 2D section are previewed, because any other 3D mesh can take seconds.
    pub fn size_preview(&self, geometry: FeatureList, mesh: MeshSettings) -> Result<SizeTask, Error> {
        self.check_edits(&geometry, &mesh)?;
        if geometry.dim == 3 && mesh.sweep.is_none() {
            return Err(Error::invalid(
                "a 3D mesh that is not swept along y can take seconds: it is built by Build mesh, not previewed",
            ));
        }
        let mut scratch = self.scratch();
        scratch.project.geometry = geometry.clone();
        scratch.project.mesh = mesh.clone();
        let sizing = scratch.sizing();
        let chain = chain_key(&scratch.default_targets());
        let sweep_dim = if mesh.sweep.is_some() { geometry.dim + 1 } else { geometry.dim };
        Ok(SizeTask {
            mesher: self.mesher.clone(),
            geometry,
            settings: mesh,
            sizing,
            history: self.history.clone(),
            chain,
            sweep_dim,
        })
    }

    /// The polarization analysis of the latest result of `physics` as a function that does not borrow the engine, so it can run beside the
    /// frame (the analysis is a pure function of the mesh and the result, ADR 56). The result it was made for is returned with it, to compare
    /// with what is on show when the answer comes back. `None` when the analysis does not apply, or there is no result or no mesh.
    pub fn polarization_task(
        &self,
        physics: &str,
    ) -> Option<(Arc<SolveResult>, impl FnOnce() -> Result<ModeAnalysis, Error> + Send + 'static)> {
        if !self.polarization_applies(physics) {
            return None;
        }
        let (model, result) = self.analysis_inputs(physics).ok()?;
        let r = result.clone();
        Some((r, move || nyx_analysis::polarization(model.mesh(), &result)))
    }
}
