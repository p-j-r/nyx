// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
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

/// A mesh above this many cells (2D, 3D) takes tens of seconds to solve with the quadratic elements and the complex eigensolver, and a
/// convergence check on a mesh made finer multiplies it again: said in the answer of `mesh_build`.
const HEAVY_CELLS: (u64, u64) = (15_000, 100_000);

/// What a solve of `targets` on a mesh of `cells` cells costs, from the solves `history` has seen (the answer of `mesh_build` and `convergence_check`
/// carries it): `solve` is `null` while there are too few runs (no number is invented), `warning` is there for a mesh that is large for the way it is
/// used (L3).
pub(crate) fn cost_of(history: &History, targets: &[String], dim: u8, cells: u64) -> serde_json::Value {
    use serde_json::json;
    let chain = chain_key(targets);
    let solve = history.estimate(&chain, dim, cells).map(
        |e| json!({ "seconds": e.seconds, "memory_bytes": e.bytes, "runs": e.runs, "extrapolated": e.extrapolated }),
    );
    let mut v = json!({ "chain": chain, "dim": dim, "cells": cells, "solve": solve });
    if let Some(w) = heavy_warning(dim, cells) {
        v["warning"] = json!(w);
    }
    v
}

fn heavy_warning(dim: u8, cells: u64) -> Option<String> {
    let limit = if dim >= 3 { HEAVY_CELLS.1 } else { HEAVY_CELLS.0 };
    (cells > limit).then(|| {
        format!(
            "{cells} cells is a large mesh: a solve takes long and a convergence_check at factor 0.6 makes the mesh several times larger (an adaptive 2D mesh 6 to 8 times). Coarsen h or adapt_scale first, or ask convergence_check with dry_run true what the finer mesh would cost"
        )
    })
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
    /// What a solve of `targets` on a mesh of `cells` cells costs, as the answer of `mesh_build` and `convergence_check` shows it: the estimate from the
    /// solves this machine has run (`null` while there are too few: no number is invented) and, for a mesh that is large for the way it is used, a
    /// warning in plain words (L3).
    pub(crate) fn cost_json(&self, targets: &[String], dim: u8, cells: u64) -> serde_json::Value {
        cost_of(&self.history, targets, dim, cells)
    }

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

#[cfg(test)]
mod cost_tests {
    use super::*;

    #[test]
    fn a_mesh_is_large_for_its_dimension_and_the_cost_has_no_number_without_runs() {
        assert!(heavy_warning(2, 15_000).is_none() && heavy_warning(2, 15_001).is_some());
        assert!(heavy_warning(3, 100_000).is_none() && heavy_warning(3, 100_001).is_some());
        let h = History::default();
        let v = cost_of(&h, &["ems".to_string()], 2, 22_714);
        assert!(v["solve"].is_null(), "no runs, no estimate: {v}");
        assert!(v["warning"].as_str().unwrap().contains("22714 cells"), "{v}");
        assert_eq!(v["chain"], "ems");
    }
}
