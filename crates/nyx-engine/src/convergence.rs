//! Mesh convergence as an engine command: the same project is solved on a finer mesh (on a scratch engine, so the project, the mesh and the
//! results of the caller are untouched) and every number that identifies a solution is compared: the eigenvalues (effective indices,
//! frequencies), the values that order the solutions and the numeric scalars of the result. This replaces the `mesh_converge` check of the
//! old session tools for every physics (ADR 61).
//!
//! By default solutions are compared by index. When the order of the modes changes between the two meshes (the complex modes of an open
//! resonator do), `Matching::Nearest` pairs each mode with the nearest eigenvalue of the other mesh in the complex plane and compares the
//! pairs, and `Matching::Track` follows the one mode nearest a given value (ADR 80).

use crate::engine::Engine;
use crate::project::Project;
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
        check(self, physics, factor, tolerance, &Matching::Index)
    }

    /// The check with the solutions paired as `matching` says.
    pub fn run_matched(
        &self,
        physics: &str,
        factor: f64,
        tolerance: Option<f64>,
        matching: &Matching,
    ) -> Result<Value, Error> {
        check(self, physics, factor, tolerance, matching)
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
        }
    }

    /// Solves `physics` on the current mesh and on one made finer by `factor`, and compares the numbers of the two results. `tolerance`,
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
        let mut coarse_engine = this.scratch(this.project.clone());
        let t0 = std::time::Instant::now();
        let coarse_summary = coarse_engine.solve_now(Some(targets.clone()))?;
        let coarse: Arc<SolveResult> = coarse_engine
            .result(physics)
            .map(|(r, _)| r.clone())
            .ok_or_else(|| Error::invalid("the solve gave no result"))?;
        let cells_coarse = coarse_engine.model().map(|m| m.mesh().n_cells());
        let mut fine_engine = this.scratch(refined(&this.project, factor));
        fine_engine.solve_now(Some(targets))?;
        let fine: Arc<SolveResult> = fine_engine
            .result(physics)
            .map(|(r, _)| r.clone())
            .ok_or_else(|| Error::invalid("the refined solve gave no result"))?;
        let cells_fine = fine_engine.model().map(|m| m.mesh().n_cells());

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
            Matching::Index => common(numbers(&coarse), numbers(&fine)),
            Matching::Nearest | Matching::Track(_) => {
                let (pairs, lone_coarse, lone_fine) = match matching {
                    Matching::Track(f) => {
                        let i = (0..coarse.solutions.len())
                            .filter_map(|i| key(&coarse, i).map(|k| (i, (k[0] - f).abs())))
                            .min_by(|x, y| x.1.total_cmp(&y.1))
                            .map(|(i, _)| i)
                            .ok_or_else(|| Error::invalid("the result has no eigenvalue or value to track"))?;
                        let a = key(&coarse, i).expect("chosen by its key");
                        let (j, d) = (0..fine.solutions.len())
                            .filter_map(|j| key(&fine, j).map(|b| (j, distance(a, b))))
                            .min_by(|x, y| x.1.total_cmp(&y.1))
                            .ok_or_else(|| Error::invalid("the refined result has no eigenvalue or value to track"))?;
                        (vec![(i, j, d)], vec![], vec![])
                    }
                    _ => pair_nearest(&coarse, &fine),
                };
                for &(i, j, _) in &pairs {
                    common(numbers_of(&coarse, i, i), numbers_of(&fine, j, i));
                }
                matched = Some((pairs, lone_coarse, lone_fine));
            }
        }
        let mut rows = vec![];
        let mut worst = 0.0f64;
        let mut all = true;
        for (name, x, y) in &compared {
            let scale = x.abs().max(y.abs()).max(1e-300);
            let rel = (x - y).abs() / if scale < 1e-12 { 1.0 } else { scale };
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
        let _ = coarse_summary;
        let mut out = json!({
            "physics": physics, "factor": factor, "tolerance": tolerance,
            "cells": { "coarse": cells_coarse, "refined": cells_fine },
            "unknowns": { "coarse": coarse.dofs, "refined": fine.dofs },
            "seconds": t0.elapsed().as_secs_f64(),
            "worst_relative": worst,
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
