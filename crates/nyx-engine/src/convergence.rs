//! Mesh convergence as an engine command: the same project is solved on a finer mesh (on a scratch engine, so the project, the mesh and the
//! results of the caller are untouched) and every number that identifies a solution is compared: the eigenvalues (effective indices,
//! frequencies), the values that order the solutions and the numeric scalars of the result. This replaces the `mesh_converge` check of the
//! old session tools for every physics (ADR 61).

use crate::engine::Engine;
use crate::project::Project;
use nyx_physics::{Error, Scalar, SolveResult};
use serde_json::{json, Value};
use std::sync::Arc;

/// The project with a finer mesh: the element size is multiplied by `factor` (below one is finer), so are the adaptive size factor and the
/// wavelength rules; the layers of a sweep get `1 / factor` times as many cells.
fn refined(p: &Project, factor: f64) -> Project {
    let mut q = p.clone();
    q.mesh.h *= factor;
    if let Some(s) = q.mesh.adapt_scale.as_mut() {
        *s *= factor;
    }
    if let Some(d) = q.mesh.max_dev.as_mut() {
        *d *= factor;
    }
    if let Some(ppw) = q.mesh.points_per_wavelength.as_mut() {
        *ppw /= factor;
    }
    if let Some(layers) = q.mesh.sweep.as_mut() {
        for l in layers {
            l.n = ((l.n as f64) / factor).ceil() as usize;
        }
    }
    q
}

/// Every named number of a result, in a fixed order: `eigenvalue[i].re`, `eigenvalue[i].im`, `value[i]` and the numeric scalars.
fn numbers(r: &SolveResult) -> Vec<(String, f64)> {
    let mut out = vec![];
    for (i, s) in r.solutions.iter().enumerate() {
        if let Some(e) = s.eigenvalue {
            out.push((format!("eigenvalue[{i}].re"), e[0]));
            if e[1] != 0.0 {
                out.push((format!("eigenvalue[{i}].im"), e[1]));
            }
        }
        if let Some(v) = s.value {
            out.push((format!("value[{i}]"), v));
        }
    }
    for (name, s) in &r.scalars {
        if let Scalar::Num(x) = s {
            if x.is_finite() {
                out.push((format!("scalar {name}"), *x));
            }
        }
    }
    out
}

/// Everything a convergence check needs, owned, so that a front-end can run it on its own thread while the engine stays usable.
#[derive(Clone)]
pub struct ConvergenceJob {
    registry: nyx_physics::Registry,
    mesher: Arc<dyn crate::meshing::Mesher>,
    library: Arc<femcore::material::Library>,
    project: Project,
}

impl ConvergenceJob {
    fn scratch(&self, p: Project) -> Engine {
        let mut e = Engine::with(self.registry.clone(), self.mesher.clone());
        e.base_library = self.library.clone();
        e.set_project(p);
        e
    }

    /// Runs the check (see `Engine::convergence_check`); blocks until both solves are done.
    pub fn run(&self, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<Value, Error> {
        check(self, physics, factor, tolerance)
    }
}

impl Engine {
    /// A snapshot of what `convergence_check` needs.
    pub fn convergence_job(&self) -> ConvergenceJob {
        ConvergenceJob { registry: self.registry.clone(), mesher: self.mesher.clone(), library: self.base_library.clone(), project: self.project.clone() }
    }

    /// Solves `physics` on the current mesh and on one made finer by `factor`, and compares the numbers of the two results. `tolerance`,
    /// when given, is the largest relative difference accepted (the difference of a quantity over the larger of the two values, or over 1 for
    /// a value near zero); every row then says whether it passed and `converged` is true when all did.
    pub fn convergence_check(&self, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<Value, Error> {
        self.convergence_job().run(physics, factor, tolerance)
    }
}

fn check(job: &ConvergenceJob, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<Value, Error> {
    {
        let this = job;
        if !(factor > 0.0 && factor < 1.0) {
            return Err(Error::invalid("`factor` is the ratio of the finer element size to the present one, between 0 and 1 (0.6 is a good step)"));
        }
        this.project.spec(physics).ok_or_else(|| Error::invalid(format!("the project has no problem for '{physics}'")))?;
        let targets = vec![physics.to_string()];
        let mut coarse_engine = this.scratch(this.project.clone());
        let t0 = std::time::Instant::now();
        let coarse_summary = coarse_engine.solve_now(Some(targets.clone()))?;
        let coarse: Arc<SolveResult> = coarse_engine.result(physics).map(|(r, _)| r.clone()).ok_or_else(|| Error::invalid("the solve gave no result"))?;
        let cells_coarse = coarse_engine.model().map(|m| m.mesh().n_cells());
        let mut fine_engine = this.scratch(refined(&this.project, factor));
        fine_engine.solve_now(Some(targets))?;
        let fine: Arc<SolveResult> = fine_engine.result(physics).map(|(r, _)| r.clone()).ok_or_else(|| Error::invalid("the refined solve gave no result"))?;
        let cells_fine = fine_engine.model().map(|m| m.mesh().n_cells());

        let (a, b) = (numbers(&coarse), numbers(&fine));
        let mut rows = vec![];
        let mut worst = 0.0f64;
        let mut all = true;
        for (name, x) in &a {
            let Some((_, y)) = b.iter().find(|(n, _)| n == name) else { continue };
            let scale = x.abs().max(y.abs()).max(1e-300);
            let rel = (x - y).abs() / if scale < 1e-12 { 1.0 } else { scale };
            worst = worst.max(rel);
            let pass = tolerance.map(|t| rel <= t);
            all &= pass.unwrap_or(true);
            rows.push(json!({ "quantity": name, "coarse": x, "refined": y, "difference": y - x, "relative": rel, "converged": pass }));
        }
        if rows.is_empty() {
            return Err(Error::invalid("the result has no eigenvalue, value or numeric scalar to compare"));
        }
        let _ = coarse_summary;
        Ok(json!({
            "physics": physics, "factor": factor, "tolerance": tolerance,
            "cells": { "coarse": cells_coarse, "refined": cells_fine },
            "unknowns": { "coarse": coarse.dofs, "refined": fine.dofs },
            "seconds": t0.elapsed().as_secs_f64(),
            "worst_relative": worst,
            "converged": tolerance.map(|_| all),
            "rows": rows,
            "note": "the project was not changed; compare with a smaller factor if the difference is still falling",
        }))
    }
}
