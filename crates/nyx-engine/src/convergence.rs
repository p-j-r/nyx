// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Mesh convergence as an engine command: the same project is solved on a finer mesh (on a scratch engine, so the project, the mesh and the
//! results of the caller are untouched) and every number that identifies a solution is compared: the eigenvalues (effective indices,
//! frequencies), the values that order the solutions and the numeric scalars of the result. This replaces the `mesh_converge` check of the
//! old session tools for every physics (ADR 61).
//!
//! By default solutions are compared by index. When the order of the modes changes between the two meshes (the complex modes of an open
//! resonator do), `Matching::Nearest` pairs each mode with the nearest eigenvalue of the other mesh in the complex plane and compares the
//! pairs, and `Matching::Track` follows the one mode nearest a given value (ADR 80).

use crate::engine::Engine;
use crate::meshing::MeshedModel;
use crate::project::Project;
use crate::queue::{SolveSummary, Stages};
use femcore::progress::Progress;
use nyx_physics::{Error, Scalar, SolveResult};
use serde_json::{json, Value};
use std::sync::Arc;

/// The project with a finer mesh: the element size is multiplied by `factor` (below one is finer), so are the adaptive size factor and the
/// wavelength rules; the layers of a sweep get `1 / factor` times as many cells. Sizes that belong to the geometry (`add_size` features) are
/// not changed (rough spot R11b): a model whose resolution is set by them is refined only where the global size rules.
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
pub(crate) fn numbers(r: &SolveResult) -> Vec<(String, f64)> {
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

/// How the solutions of the two meshes are paired.
#[derive(Clone, Debug, PartialEq, Default)]
pub enum Matching {
    /// solution `i` of the present mesh against solution `i` of the finer one (ADR 61)
    #[default]
    Index,
    /// each solution against the nearest eigenvalue of the other mesh in the complex plane; only mutual nearest neighbours are a pair, the
    /// others are listed as unpaired
    Nearest,
    /// only the solution whose eigenvalue has its real part nearest this value (a frequency in Hz for a vibration study), against the
    /// nearest eigenvalue of the finer mesh
    Track(f64),
}

impl Matching {
    /// From the arguments of the command: `match` is `index` or `nearest`, `track` is `{ "f": <value> }` and implies `nearest`.
    pub fn from_args(args: &Value) -> Result<Matching, Error> {
        let named = match args.get("match").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => match v.as_str() {
                Some("index") => Some(Matching::Index),
                Some("nearest") => Some(Matching::Nearest),
                _ => {
                    return Err(Error::invalid(
                        "`match` is \"index\" (compare solutions by position) or \"nearest\" (pair each mode with the nearest eigenvalue of the finer mesh)",
                    ))
                }
            },
        };
        match args.get("track").filter(|v| !v.is_null()) {
            None => Ok(named.unwrap_or_default()),
            Some(t) => {
                let f = t.get("f").and_then(Value::as_f64).filter(|f| f.is_finite()).ok_or_else(|| {
                    Error::invalid(
                        "`track` is { \"f\": <value> }: the mode whose eigenvalue is nearest this value (Hz for a vibration study) is followed",
                    )
                })?;
                if named == Some(Matching::Index) {
                    return Err(Error::invalid(
                        "`track` follows a mode by its eigenvalue; it cannot be used with `match`: \"index\"",
                    ));
                }
                Ok(Matching::Track(f))
            }
        }
    }
}

/// The complex number a solution is matched by: its eigenvalue, or its value when it has none.
fn key(r: &SolveResult, i: usize) -> Option<[f64; 2]> {
    let s = &r.solutions[i];
    s.eigenvalue.or(s.value.map(|v| [v, 0.0]))
}

/// Distance of two eigenvalues in the complex plane over the larger modulus.
fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
    let scale = (a[0].hypot(a[1])).max(b[0].hypot(b[1]));
    if scale < 1e-300 {
        d
    } else {
        d / scale
    }
}

/// A scalar that belongs to one mode carries the one-based number of the mode at its end, before a unit in brackets if it has one
/// (`Q3`, `layer fraction 3`, `f re 3 (Hz)`): the name before the number, the unit, and the zero-based mode. `None` for any other scalar.
fn mode_scalar(name: &str, n_solutions: usize) -> Option<(String, String, usize)> {
    let (body, unit) = match name.rfind(" (") {
        Some(k) if name.ends_with(')') => (&name[..k], &name[k..]),
        _ => (name, ""),
    };
    let digits = body.len() - body.trim_end_matches(|c: char| c.is_ascii_digit()).len();
    if digits == 0 || digits == body.len() {
        return None;
    }
    let n: usize = body[body.len() - digits..].parse().ok()?;
    (n >= 1 && n <= n_solutions).then(|| (body[..body.len() - digits].to_string(), unit.to_string(), n - 1))
}

/// The numbers of solution `i`, named for solution `as_index` (so that a pair of different positions compares under one name).
fn numbers_of(r: &SolveResult, i: usize, as_index: usize) -> Vec<(String, f64)> {
    let mut out = vec![];
    let s = &r.solutions[i];
    if let Some(e) = s.eigenvalue {
        out.push((format!("eigenvalue[{as_index}].re"), e[0]));
        if e[1] != 0.0 {
            out.push((format!("eigenvalue[{as_index}].im"), e[1]));
        }
    }
    if let Some(v) = s.value {
        out.push((format!("value[{as_index}]"), v));
    }
    for (name, sc) in &r.scalars {
        if let (Scalar::Num(x), Some((head, unit, mode))) = (sc, mode_scalar(name, r.solutions.len())) {
            if mode == i && x.is_finite() {
                out.push((format!("scalar {head}{}{unit}", as_index + 1), *x));
            }
        }
    }
    out
}

/// An imaginary part below this times the real part of its eigenvalue, on both meshes, is rounding (a quality factor above 5e7).
const IM_FLOOR: f64 = 1e-8;

/// Zero-based mode index of a quality-factor scalar name (after `"scalar "`), such as `"Q1"` or `"Q 1 = Re/(2|Im|) of neff"`.
fn q_mode(rest: &str) -> Option<usize> {
    // the waveguide study names it `Q 1 = Re/(2|Im|) of neff`: the quality factor of mode 1, with the formula after the number
    if let Some(n) = rest
        .strip_prefix("Q ")
        .map(|r| r.chars().take_while(char::is_ascii_digit).collect::<String>())
        .and_then(|d| d.parse::<usize>().ok())
    {
        return (n >= 1).then(|| n - 1);
    }
    mode_scalar(rest, usize::MAX).and_then(|(head, _, mode)| (head.trim() == "Q").then_some(mode))
}

/// The modes (by the index their rows carry) whose `eigenvalue[i].im` (or quality factor `Q`, when the solver reports no
/// imaginary eigenvalue) is at the rounding level on both meshes.
fn noise_modes(compared: &[(String, f64, f64)]) -> std::collections::BTreeSet<usize> {
    let mut out = std::collections::BTreeSet::new();
    let q_floor = 1.0 / (2.0 * IM_FLOOR);
    for (name, x, y) in compared {
        if let Some(i) =
            name.strip_prefix("eigenvalue[").and_then(|s| s.strip_suffix("].im")).and_then(|s| s.parse::<usize>().ok())
        {
            let re_name = format!("eigenvalue[{i}].re");
            let Some(re) = compared.iter().find(|(n, _, _)| *n == re_name).map(|(_, a, b)| a.abs().max(b.abs())) else {
                continue;
            };
            if re > 0.0 && x.abs() <= IM_FLOOR * re && y.abs() <= IM_FLOOR * re {
                out.insert(i);
            }
            continue;
        }
        if let Some(i) = name.strip_prefix("scalar ").and_then(q_mode) {
            let im_name = format!("eigenvalue[{i}].im");
            if !compared.iter().any(|(n, _, _)| *n == im_name) && *x >= q_floor && *y >= q_floor {
                out.insert(i);
            }
        }
    }
    out
}

/// Whether a row belongs to a noise mode: its `im`, or a scalar made of it (the quality factor `Q<n>`, `f im <n> (Hz)`).
fn is_noise(name: &str, noise: &std::collections::BTreeSet<usize>) -> bool {
    if noise.is_empty() {
        return false;
    }
    if let Some(i) =
        name.strip_prefix("eigenvalue[").and_then(|s| s.strip_suffix("].im")).and_then(|s| s.parse::<usize>().ok())
    {
        return noise.contains(&i);
    }
    let Some(rest) = name.strip_prefix("scalar ") else { return false };
    if let Some(i) = q_mode(rest) {
        return noise.contains(&i);
    }
    mode_scalar(rest, usize::MAX).is_some_and(|(head, _, mode)| head.trim() == "f im" && noise.contains(&mode))
}

/// A pair of solutions: index on the present mesh, index on the finer mesh, distance of the eigenvalues.
type Pair = (usize, usize, f64);

/// The pairs of mutual nearest neighbours, and the solutions of each mesh left without a partner.
fn pair_nearest(coarse: &SolveResult, fine: &SolveResult) -> (Vec<Pair>, Vec<usize>, Vec<usize>) {
    let nearest = |from: &SolveResult, i: usize, to: &SolveResult| -> Option<(usize, f64)> {
        let a = key(from, i)?;
        (0..to.solutions.len())
            .filter_map(|j| key(to, j).map(|b| (j, distance(a, b))))
            .min_by(|x, y| x.1.total_cmp(&y.1))
    };
    let mut pairs = vec![];
    for i in 0..coarse.solutions.len() {
        if let Some((j, d)) = nearest(coarse, i, fine) {
            if nearest(fine, j, coarse).is_some_and(|(back, _)| back == i) {
                pairs.push((i, j, d));
            }
        }
    }
    let lone_coarse = (0..coarse.solutions.len()).filter(|i| !pairs.iter().any(|p| p.0 == *i)).collect();
    let lone_fine = (0..fine.solutions.len()).filter(|j| !pairs.iter().any(|p| p.1 == *j)).collect();
    (pairs, lone_coarse, lone_fine)
}

/// Everything a convergence check needs, owned, so that a front-end can run it on its own thread while the engine stays usable.
#[derive(Clone)]
pub struct ConvergenceJob {
    registry: nyx_physics::Registry,
    mesher: Arc<dyn crate::meshing::Mesher>,
    library: Arc<femcore::material::Library>,
    project: Project,
    /// the stage cache and the key under which the coarse solve may already be there (set by `Engine::start_convergence_with`)
    pub(crate) reuse: Option<(Arc<Stages>, String)>,
    /// what solves have cost on this machine: the cost of the finer mesh is looked up before it is solved
    history: Arc<crate::estimate::History>,
    /// the largest finer mesh, and the longest estimated check, the caller accepts (L3)
    pub(crate) limits: Limits,
}

/// What a convergence check may cost: it is refused, before any solve, when the finer mesh has more cells than `max_cells` or when the solves it
/// needs are estimated (from the solves this machine has run) to take longer than `max_seconds`. No estimate, no refusal on time.
#[derive(Clone, Debug, Default)]
pub struct Limits {
    pub max_cells: Option<u64>,
    pub max_seconds: Option<f64>,
}

impl Limits {
    /// `max_cells` and `max_seconds` of the arguments of the command.
    pub fn from_args(args: &Value) -> Result<Limits, Error> {
        let max_cells = match args.get("max_cells").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => Some(
                v.as_u64()
                    .filter(|n| *n > 0)
                    .ok_or_else(|| Error::invalid("`max_cells` is a number of cells, above zero"))?,
            ),
        };
        let max_seconds = match args.get("max_seconds").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => Some(
                v.as_f64()
                    .filter(|s| *s > 0.0 && s.is_finite())
                    .ok_or_else(|| Error::invalid("`max_seconds` is a number of seconds, above zero"))?,
            ),
        };
        Ok(Limits { max_cells, max_seconds })
    }
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
        self.run_with(&Progress::new(), physics, factor, tolerance, &Matching::Index)
    }

    /// The check with the solutions paired as `matching` says.
    pub fn run_matched(
        &self,
        physics: &str,
        factor: f64,
        tolerance: Option<f64>,
        matching: &Matching,
    ) -> Result<Value, Error> {
        self.run_with(&Progress::new(), physics, factor, tolerance, matching)
    }

    /// Runs the check on the calling thread under `progress`: its two solves are two stages of it, and its cancel flag stops either.
    pub fn run_with(
        &self,
        progress: &Arc<Progress>,
        physics: &str,
        factor: f64,
        tolerance: Option<f64>,
        matching: &Matching,
    ) -> Result<Value, Error> {
        let _installed = femcore::progress::install(progress);
        check(self, progress, physics, factor, tolerance, matching)
    }

    /// Meshes `project` on the calling thread.
    fn mesh_inline(&self, project: &Project) -> Result<Arc<MeshedModel>, Error> {
        let scratch = self.scratch(project.clone());
        project.check(&self.registry)?;
        let sizing = scratch.sizing();
        Ok(Arc::new(self.mesher.mesh(&project.geometry, &project.mesh, sizing.as_ref())?))
    }

    /// Solves `targets` of `project` on `model`, on the calling thread, under `progress`; returns the result of `physics`.
    fn solve_on(
        &self,
        project: Project,
        model: &Arc<MeshedModel>,
        targets: &[String],
        physics: &str,
        progress: &Arc<Progress>,
    ) -> Result<Arc<SolveResult>, Error> {
        let scratch = self.scratch(project.clone());
        let library = scratch.library()?;
        let results = crate::solve::run_chain(&self.registry, &project, model, targets, &library, progress.clone())?;
        let result = results
            .into_iter()
            .find(|r| r.physics == physics)
            .ok_or_else(|| Error::invalid("the solve gave no result"))?;
        Ok(Arc::new(result))
    }

    /// Meshes `project` and solves `targets` on the calling thread, under `progress`; returns the mesh and the result of `physics`.
    fn solve_inline(
        &self,
        project: Project,
        targets: &[String],
        physics: &str,
        progress: &Arc<Progress>,
    ) -> Result<(Arc<MeshedModel>, Arc<SolveResult>), Error> {
        let model = self.mesh_inline(&project)?;
        let result = self.solve_on(project, &model, targets, physics, progress)?;
        Ok((model, result))
    }

    /// What the check would cost, without solving anything: the present mesh and the finer one are built (milliseconds against a solve), and the
    /// answer gives their cells, what a solve of each is expected to cost from the solves this machine has run (`null` while there are too few: no
    /// number is invented) and the factor between the cells. A script that is unsure of a mesh asks this before `convergence_check` (L3): an
    /// adaptive 2D mesh gets 6 to 8 times as many cells at factor 0.6, not 2.8 times.
    pub fn plan(&self, physics: &str, factor: f64) -> Result<Value, Error> {
        if !(factor > 0.0 && factor < 1.0) {
            return Err(Error::invalid("`factor` is the ratio of the finer element size to the present one, between 0 and 1 (0.6 is a good step)"));
        }
        self.project
            .spec(physics)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{physics}'")))?;
        let targets = vec![physics.to_string()];
        let coarse = self.mesh_inline(&self.project)?;
        let fine = self.mesh_inline(&refined(&self.project, factor))?;
        let (c, f) = (coarse.mesh().n_cells() as u64, fine.mesh().n_cells() as u64);
        let (dc, df) = (coarse.mesh().dim() as u8, fine.mesh().dim() as u8);
        let (cost_c, cost_f) = (
            crate::preview::cost_of(&self.history, &targets, dc, c),
            crate::preview::cost_of(&self.history, &targets, df, f),
        );
        let secs = |v: &Value| v["solve"]["seconds"].as_f64();
        let total = secs(&cost_c).zip(secs(&cost_f)).map(|(a, b)| a + b);
        Ok(json!({
            "physics": physics, "factor": factor, "dry_run": true,
            "cells": { "coarse": c, "refined": f },
            "refined_over_present": f as f64 / c.max(1) as f64,
            "cost": { "coarse": cost_c, "refined": cost_f },
            "seconds_estimate": total,
            "note": "nothing was solved; the check solves both meshes (the present one is skipped when it was solved alone already). Run it with max_cells or max_seconds to have it refused when it is larger than you want",
        }))
    }
}

impl Engine {
    /// A snapshot of what `convergence_check` needs.
    pub fn convergence_job(&self) -> ConvergenceJob {
        ConvergenceJob {
            registry: self.registry.clone(),
            mesher: self.mesher.clone(),
            library: self.base_library.clone(),
            project: self.project.clone(),
            reuse: None,
            history: self.history.clone(),
            limits: Limits::default(),
        }
    }

    /// Solves `physics` on the current mesh and on one made finer by `factor`, and compares the numbers of the two results, on the calling thread
    /// and outside the queue (the `convergence_check` command goes through `start_convergence`). `tolerance`,
    /// when given, is the largest relative difference accepted (the difference of a quantity over the larger of the two values, or over 1 for
    /// a value near zero); every row then says whether it passed and `converged` is true when all did.
    pub fn convergence_check(&self, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<Value, Error> {
        self.convergence_job().run(physics, factor, tolerance)
    }

    /// `convergence_check` with the modes of the two meshes paired by their eigenvalues (`Matching::Nearest`) or one mode followed
    /// (`Matching::Track`): `worst_relative` and `converged` are then over the paired modes only, and the answer lists the pairs and the
    /// solutions without a partner.
    pub fn convergence_check_matched(
        &self,
        physics: &str,
        factor: f64,
        tolerance: Option<f64>,
        matching: &Matching,
    ) -> Result<Value, Error> {
        self.convergence_job().run_matched(physics, factor, tolerance, matching)
    }
}

fn check(
    job: &ConvergenceJob,
    progress: &Arc<Progress>,
    physics: &str,
    factor: f64,
    tolerance: Option<f64>,
    matching: &Matching,
) -> Result<Value, Error> {
    {
        let this = job;
        if !(factor > 0.0 && factor < 1.0) {
            return Err(Error::invalid("`factor` is the ratio of the finer element size to the present one, between 0 and 1 (0.6 is a good step)"));
        }
        this.project
            .spec(physics)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{physics}'")))?;
        let targets = vec![physics.to_string()];
        let t0 = std::time::Instant::now();
        // L3: the finer mesh is built first (milliseconds against a solve), so a check that is larger than the caller accepts is refused before
        // any solve, and the answer to a refusal says what the check would have cost
        let fine_project = refined(&this.project, factor);
        let fine_model = {
            let _stage =
                femcore::progress::stage(format!("convergence check: the finer mesh (element size times {factor})"));
            this.mesh_inline(&fine_project)?
        };
        let (fine_cells, fine_dim) = (fine_model.mesh().n_cells() as u64, fine_model.mesh().dim() as u8);
        if let Some(max) = this.limits.max_cells.filter(|m| fine_cells > *m) {
            return Err(Error::invalid(format!(
                "the finer mesh would have {fine_cells} cells, more than max_cells {max}: the check was not started. Coarsen the present mesh first (a larger h or adapt_scale), use a factor nearer 1 (now {factor}), or raise max_cells"
            )));
        }
        if let Some(max) = this.limits.max_seconds {
            let cost = crate::preview::cost_of(&this.history, &targets, fine_dim, fine_cells);
            if let Some(s) = cost["solve"]["seconds"].as_f64().filter(|s| *s > max) {
                return Err(Error::invalid(format!(
                    "the solve on the finer mesh ({fine_cells} cells) is expected to take about {s:.0} s from the solves this machine has run, more than max_seconds {max}: the check was not started. Coarsen the mesh, use a factor nearer 1 (now {factor}), or raise max_seconds"
                )));
            }
        }
        // the two solves are stages of the one job, run here under its progress, so the check can be watched and cancelled; the first is skipped
        // when an earlier request solved this physics alone on this mesh (found by content, so the numbers are the ones this solve would give)
        // a check that pairs modes by their eigenvalues needs the solutions of the present mesh, which the cache does not keep (only numbers)
        let reusable = if matches!(matching, Matching::Index) {
            this.reuse.as_ref().and_then(|(stages, key)| stages.solve(key))
        } else {
            None
        };
        let (coarse, self_reused, coarse_result) = match reusable {
            Some(s) => {
                let _stage = femcore::progress::stage(
                    "convergence check: the solve on the present mesh (reused from an earlier solve)",
                );
                (s, true, None)
            }
            None => {
                let _stage = femcore::progress::stage("convergence check: the solve on the present mesh");
                let (model, r) = this.solve_inline(this.project.clone(), &targets, physics, progress)?;
                (Arc::new(SolveSummary::of(&r, model.mesh().n_cells())), false, Some(r))
            }
        };
        let cells_coarse = Some(coarse.cells);
        let fine = {
            let _stage = femcore::progress::stage(format!(
                "convergence check: the solve on the finer mesh (element size times {factor})"
            ));
            this.solve_on(fine_project, &fine_model, &targets, physics, progress)?
        };
        let cells_fine = Some(fine_model.mesh().n_cells());

        // the rows to compare: (name, present, finer) and, when the modes are paired, which solutions were paired
        let mut compared: Vec<(String, f64, f64)> = vec![];
        let mut matched = None;
        let mut common = |a: Vec<(String, f64)>, b: Vec<(String, f64)>| {
            for (name, x) in a {
                if let Some((_, y)) = b.iter().find(|(n, _)| *n == name) {
                    compared.push((name, x, *y));
                }
            }
        };
        match matching {
            Matching::Index => common(coarse.numbers.clone(), numbers(&fine)),
            Matching::Nearest | Matching::Track(_) => {
                let coarse_r = coarse_result.as_ref().expect("a check that pairs modes solves the present mesh itself");
                let (pairs, lone_coarse, lone_fine) = match matching {
                    Matching::Track(f) => {
                        let i = (0..coarse_r.solutions.len())
                            .filter_map(|i| key(coarse_r, i).map(|k| (i, (k[0] - f).abs())))
                            .min_by(|x, y| x.1.total_cmp(&y.1))
                            .map(|(i, _)| i)
                            .ok_or_else(|| Error::invalid("the result has no eigenvalue or value to track"))?;
                        let a = key(coarse_r, i).expect("chosen by its key");
                        let (j, d) = (0..fine.solutions.len())
                            .filter_map(|j| key(&fine, j).map(|b| (j, distance(a, b))))
                            .min_by(|x, y| x.1.total_cmp(&y.1))
                            .ok_or_else(|| Error::invalid("the refined result has no eigenvalue or value to track"))?;
                        (vec![(i, j, d)], vec![], vec![])
                    }
                    _ => pair_nearest(coarse_r, &fine),
                };
                for &(i, j, _) in &pairs {
                    common(numbers_of(coarse_r, i, i), numbers_of(&fine, j, i));
                }
                matched = Some((pairs, lone_coarse, lone_fine));
            }
        }
        // L2: a mode whose imaginary part is zero up to rounding (a lossless guide solved by the complex solver, because a crystal axis is tilted)
        // has `im` of 1e-10 or so that differs between two meshes by 100 percent, and a quality factor of 1e10 made of it. Such a mode is not a
        // failure of the mesh: its `im` and its quality factors are listed but not counted. The floor is what the eigensolver resolves: an `im`
        // below it times the real part on both meshes (a Q above 5e7) is not told from rounding.
        let noise = noise_modes(&compared);
        let mut rows = vec![];
        let mut worst = 0.0f64;
        let mut worst_name = String::new();
        let mut all = true;
        for (name, x, y) in &compared {
            let scale = x.abs().max(y.abs()).max(1e-300);
            let rel = (x - y).abs() / if scale < 1e-12 { 1.0 } else { scale };
            if is_noise(name, &noise) {
                let mut row = json!({ "quantity": name, "coarse": x, "refined": y, "difference": y - x, "relative": rel, "converged": tolerance.map(|_| true), "noise": true });
                row["why"] = json!("imaginary part at the rounding level of this mode: not counted");
                rows.push(row);
                continue;
            }
            if rel > worst || worst_name.is_empty() {
                worst_name = name.clone();
            }
            worst = worst.max(rel);
            let pass = tolerance.map(|t| rel <= t);
            all &= pass.unwrap_or(true);
            rows.push(json!({ "quantity": name, "coarse": x, "refined": y, "difference": y - x, "relative": rel, "converged": pass }));
        }
        if rows.is_empty() {
            return Err(Error::invalid(match matching {
                Matching::Index => "the result has no eigenvalue, value or numeric scalar to compare",
                _ => "no mode of the present mesh could be paired with a mode of the finer mesh",
            }));
        }
        let mut out = json!({
            "physics": physics, "factor": factor, "tolerance": tolerance,
            "cells": { "coarse": cells_coarse, "refined": cells_fine },
            "unknowns": { "coarse": coarse.dofs, "refined": fine.dofs },
            "coarse_reused": self_reused,
            "seconds": t0.elapsed().as_secs_f64(),
            "worst_relative": worst,
            "worst_quantity": worst_name,
            "quantities": rows.len(),
            "noise_modes": noise,
            "matching": "index",
            "converged": tolerance.map(|_| all),
            "rows": rows,
            "note": "the project was not changed; compare with a smaller factor if the difference is still falling",
        });
        if let Some((pairs, lone_coarse, lone_fine)) = matched {
            out["match"] = json!(if matches!(matching, Matching::Track(_)) { "track" } else { "nearest" });
            if let Matching::Track(f) = matching {
                out["track"] = json!({ "f": f });
            }
            out["pairs"] = pairs.iter().map(|(i, j, d)| json!({ "coarse": i, "refined": j, "distance": d })).collect();
            out["unpaired"] = json!({ "coarse": lone_coarse, "refined": lone_fine });
            out["note"] = json!("the project was not changed; modes are paired by their eigenvalues: a quantity carries the index of the mode on the present mesh and `pairs` gives its index on the finer one; rows, worst_relative and converged cover the paired modes only");
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(v: &[(&str, f64, f64)]) -> Vec<(String, f64, f64)> {
        v.iter().map(|(n, a, b)| (n.to_string(), *a, *b)).collect()
    }

    #[test]
    fn an_imaginary_part_at_the_rounding_level_is_noise_and_so_is_the_quality_factor_made_of_it() {
        // a lossless guide solved by the complex solver: im of 1e-10 against 3e-10 and a quality factor of 1e10
        let c = rows(&[
            ("eigenvalue[0].re", 2.95, 2.9501),
            ("eigenvalue[0].im", 1.0e-10, 3.0e-10),
            ("scalar Q 1 = Re/(2|Im|) of neff", 1.5e10, 5.0e9),
            ("eigenvalue[1].re", 2.4, 2.4001),
            ("eigenvalue[1].im", 1.0e-4, 1.1e-4),
            ("scalar Q 2 = Re/(2|Im|) of neff", 1.2e4, 1.1e4),
            ("scalar Q1", 1.5e10, 5.0e9),
            ("scalar Q2", 1.2e4, 1.1e4),
        ]);
        let n = noise_modes(&c);
        assert_eq!(n.iter().copied().collect::<Vec<_>>(), vec![0], "only the mode with a lossless im: {n:?}");
        assert!(
            is_noise("eigenvalue[0].im", &n)
                && is_noise("scalar Q 1 = Re/(2|Im|) of neff", &n)
                && is_noise("scalar Q1", &n)
        );
        // the real part of the same mode and everything of a lossy mode count
        assert!(!is_noise("eigenvalue[0].re", &n));
        assert!(
            !is_noise("eigenvalue[1].im", &n)
                && !is_noise("scalar Q 2 = Re/(2|Im|) of neff", &n)
                && !is_noise("scalar Q2", &n)
        );
        // `Q 10` is mode 10, not mode 1
        assert!(!is_noise("scalar Q 10 = Re/(2|Im|) of neff", &n));
    }

    #[test]
    fn a_mode_with_loss_is_never_noise_and_a_vibration_mode_is_named_by_its_frequency_scalars() {
        // a cavity with a quality factor of 1e6: im / re = 5e-7, far above the floor
        let c = rows(&[("eigenvalue[0].re", 1.0, 1.0), ("eigenvalue[0].im", 5.0e-7, 5.4e-7)]);
        assert!(noise_modes(&c).is_empty());
        // the solid study names its parts `f im 1 (Hz)`
        let s = rows(&[
            ("eigenvalue[0].re", 4.6e8, 4.6e8),
            ("eigenvalue[0].im", 1.0e-2, 4.0e-2),
            ("scalar f im 1 (Hz)", 1.0e-2, 4.0e-2),
            ("scalar f re 1 (Hz)", 4.6e8, 4.6e8),
        ]);
        let n = noise_modes(&s);
        assert_eq!(n.len(), 1);
        assert!(is_noise("scalar f im 1 (Hz)", &n) && !is_noise("scalar f re 1 (Hz)", &n));
    }

    #[test]
    fn coupled_beta_rounding_and_quality_factors_without_imaginary_eigenvalue_are_noise() {
        // |im / re| = 5e-9 (between 1e-9 and 1e-8, from coupled-beta fixed-point rounding on a tilted crystal) is noise,
        // while |im / re| = 5e-7 is a physical loss and never noise
        let c = rows(&[
            ("eigenvalue[0].re", 2.4, 2.4001),
            ("eigenvalue[0].im", 1.2e-8, -1.2e-8),
            ("scalar Q 1 = Re/(2|Im|) of neff", 1.0e8, 1.0e8),
            ("eigenvalue[1].re", 2.0, 2.0),
            ("eigenvalue[1].im", 1.0e-6, 1.0e-6),
            ("scalar Q 2 = Re/(2|Im|) of neff", 1.0e6, 1.0e6),
        ]);
        let n = noise_modes(&c);
        assert_eq!(n.iter().copied().collect::<Vec<_>>(), vec![0], "{n:?}");
        assert!(is_noise("eigenvalue[0].im", &n) && is_noise("scalar Q 1 = Re/(2|Im|) of neff", &n));
        assert!(!is_noise("eigenvalue[1].im", &n) && !is_noise("scalar Q 2 = Re/(2|Im|) of neff", &n));

        // a solver that reports `scalar Q1` without `eigenvalue[0].im` (such as solid PML) marks Q >= 5e7 as noise
        let pml = rows(&[
            ("value[0]", 9.47e9, 9.4701e9),
            ("scalar Q1", 6.0e7, 1.5e8),
            ("value[1]", 9.55e9, 9.5501e9),
            ("scalar Q2", 1.2e4, 1.1e4),
        ]);
        let np = noise_modes(&pml);
        assert_eq!(np.iter().copied().collect::<Vec<_>>(), vec![0], "{np:?}");
        assert!(is_noise("scalar Q1", &np) && !is_noise("value[0]", &np) && !is_noise("scalar Q2", &np));
    }
}
