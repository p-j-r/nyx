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
use nyx_physics::{CapacitanceMatrixData, Error, Scalar, SolveResult, CAPACITANCE_KEY};
use serde_json::{json, Value};
use std::sync::Arc;

const E_CHARGE: f64 = 1.602_176_634e-19;
const PLANCK_H: f64 = 6.626_070_15e-34;

fn check_factor(factor: f64) -> Result<(), Error> {
    if !(factor.is_finite() && factor > 0.0 && factor != 1.0) {
        return Err(Error::invalid(
            "`factor` is the ratio of the comparison element size to the present one: positive and not 1 (0.6 refines, 1.5 coarsens)",
        ));
    }
    Ok(())
}

fn project_has_multi_pass_amr(p: &Project) -> bool {
    p.physics.iter().any(|s| {
        s.physics == "dfs"
            && s.study.kind == "electrostatic"
            && s.study.params.get("adapt_passes").and_then(Value::as_i64).unwrap_or(1) > 1
    })
}

fn force_single_adapt_pass(p: &mut Project) {
    for s in &mut p.physics {
        if s.physics == "dfs" && s.study.kind == "electrostatic" {
            if let Some(obj) = s.study.params.as_object_mut() {
                if obj.get("adapt_passes").and_then(Value::as_i64).unwrap_or(1) > 1 {
                    obj.insert("adapt_passes".into(), json!(1));
                }
            }
        }
    }
}

/// The project with its element size multiplied by `factor` (`< 1` refines, `> 1` coarsens), along with the chord limit,
/// wavelength rules, and sweep layer counts. Note `adapt_scale` is intentionally NOT multiplied by `factor` because the adaptive mesher
/// already applies `adapt_scale` onto `h`, `Op::Size`, and `wave_size` (multiplying `adapt_scale` too would refine by `factor^2`).
/// When `scale_sizes` is true (the default), non-suppressed `add_size` (`Op::Size`) features on the scratch copy also scale `h` and
/// `max_deviation` by `factor` and `min_layers` by `1 / factor`.
fn refined(p: &Project, factor: f64, scale_sizes: bool, sizing_order: Option<u8>) -> Project {
    use geom::expr::Expr;
    use geom::model::Op;

    let mut q = p.clone();
    force_single_adapt_pass(&mut q);
    q.mesh.h *= factor;
    if let Some(d) = q.mesh.max_dev.as_mut() {
        *d *= factor;
    }
    if q.mesh.physics_sizing {
        if let Some(ppw) = q.mesh.points_per_wavelength.as_mut() {
            *ppw /= factor;
        } else if let Some(ord) = sizing_order {
            let base_ppw = if ord >= 2 { 6.0 } else { 12.0 };
            q.mesh.points_per_wavelength = Some(base_ppw / factor);
        }
    } else if let Some(ppw) = q.mesh.points_per_wavelength.as_mut() {
        *ppw /= factor;
    }
    if let Some(gc) = q.mesh.gap_cells.as_mut() {
        *gc /= factor;
    }
    if let Some(layers) = q.mesh.sweep.as_mut() {
        for l in layers {
            l.n = ((l.n as f64) / factor).ceil().max(1.0) as usize;
        }
    }
    if scale_sizes {
        if let Some(dz) = q.mesh.dz_min.as_mut() {
            *dz *= factor;
        }
        for feat in &mut q.geometry.features {
            if feat.suppressed {
                continue;
            }
            if let Op::Size { h, max_deviation, min_layers, .. } = &mut feat.op {
                if let Some(h) = h.as_mut() {
                    *h = match h {
                        Expr::Num(v) => Expr::Num(*v * factor),
                        Expr::Text(s) => Expr::Text(format!("({s}) * {factor:?}")),
                    };
                }
                if let Some(d) = max_deviation.as_mut() {
                    *d = match d {
                        Expr::Num(v) => Expr::Num(*v * factor),
                        Expr::Text(s) => Expr::Text(format!("({s}) * {factor:?}")),
                    };
                }
                if let Some(m) = min_layers.as_mut() {
                    *m = match m {
                        Expr::Num(v) => Expr::Num((*v / factor).ceil().max(1.0)),
                        Expr::Text(s) if factor < 1.0 => Expr::Text(format!("({s}) / {factor:?}")),
                        Expr::Text(s) => Expr::Text(format!("max(1, ({s}) / {factor:?})")),
                    };
                }
            } else if let Op::Pml { layers, .. } = &mut feat.op {
                *layers = match &*layers {
                    Expr::Num(v) if *v > 0.0 => Expr::Num((*v / factor).ceil().max(5.0)),
                    Expr::Num(v) => Expr::Num(*v),
                    Expr::Text(s) if s.trim().eq_ignore_ascii_case("auto") || s.trim().starts_with('-') => {
                        Expr::Text(s.clone())
                    }
                    Expr::Text(s) if factor < 1.0 => Expr::Text(format!("({s}) / {factor:?}")),
                    Expr::Text(s) => Expr::Text(format!("max(5, ({s}) / {factor:?})")),
                };
            }
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

/// Signature and eigenvalue target for following modes across meshes (`Matching::Track`).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct TrackFilter {
    /// Target real part of the eigenvalue (`f` in Hz, `n_b`, or `n_eff`).
    pub target: Option<f64>,
    /// Region name whose energy share (`share in <region> N`) must be at least `min_share`.
    pub region: Option<String>,
    /// Minimum energy share in `region`, in `[0.0, 1.0]`.
    pub min_share: Option<f64>,
    /// Component name whose energy share (`<component> share N`) must be at least `min_component`.
    pub component: Option<String>,
    /// Minimum energy share in `component`, in `[0.0, 1.0]`.
    pub min_component: Option<f64>,
    /// Maximum absorbing-layer energy fraction (`layer fraction N`), in `[0.0, 1.0]`.
    pub max_layer_fraction: Option<f64>,
}

fn scalar_num(r: &SolveResult, name: &str) -> Option<f64> {
    r.scalars.iter().find_map(|(k, v)| match v {
        Scalar::Num(x) if k == name && x.is_finite() => Some(*x),
        _ => None,
    })
}

impl TrackFilter {
    /// Whether any per-mode scalar signature predicate is set.
    pub fn has_signature(&self) -> bool {
        self.region.is_some()
            || self.min_share.is_some()
            || self.component.is_some()
            || self.min_component.is_some()
            || self.max_layer_fraction.is_some()
    }

    /// Validates all target and threshold invariants.
    pub fn validate(&self) -> Result<(), Error> {
        if let Some(t) = self.target {
            if !t.is_finite() {
                return Err(Error::invalid("`track` target (`f`, `n_b` or `n_eff`) must be a finite number"));
            }
        }
        match (&self.region, self.min_share) {
            (Some(r), _) if r.trim().is_empty() => {
                return Err(Error::invalid("`track.region` must be a non-empty region name"));
            }
            (Some(_), None) => {
                return Err(Error::invalid(
                    "`track.region` requires `min_share` in [0, 1] (for example, `{ \"region\": \"Si\", \"min_share\": 0.5 }`)",
                ));
            }
            (None, Some(_)) => {
                return Err(Error::invalid("`track.min_share` requires `region`"));
            }
            (Some(_), Some(s)) if !s.is_finite() || !(0.0..=1.0).contains(&s) => {
                return Err(Error::invalid("`track.min_share` must be between 0 and 1"));
            }
            _ => {}
        }
        match (&self.component, self.min_component) {
            (Some(c), _) if c.trim().is_empty() => {
                return Err(Error::invalid("`track.component` must be a non-empty component name"));
            }
            (Some(_), None) => {
                return Err(Error::invalid(
                    "`track.component` requires `min_component` in [0, 1] (for example, `{ \"component\": \"Ex\", \"min_component\": 0.5 }`)",
                ));
            }
            (None, Some(_)) => {
                return Err(Error::invalid("`track.min_component` requires `component`"));
            }
            (Some(_), Some(s)) if !s.is_finite() || !(0.0..=1.0).contains(&s) => {
                return Err(Error::invalid("`track.min_component` must be between 0 and 1"));
            }
            _ => {}
        }
        if let Some(m) = self.max_layer_fraction {
            if !m.is_finite() || !(0.0..=1.0).contains(&m) {
                return Err(Error::invalid("`track.max_layer_fraction` must be between 0 and 1"));
            }
        }
        if self.target.is_none() && !self.has_signature() {
            return Err(Error::invalid(
                "`track` needs at least one target (`f`, `n_b` or `n_eff`) or mode filter (`region`/`min_share`, `component`/`min_component`, `max_layer_fraction`)",
            ));
        }
        Ok(())
    }

    /// Parses a `track` JSON object.
    pub fn from_value(t: &Value) -> Result<TrackFilter, Error> {
        let Some(map) = t.as_object() else {
            return Err(Error::invalid(
                "`track` must be an object (e.g. `{ \"f\": 194.4e12 }` or `{ \"f\": 194.4e12, \"region\": \"Si\", \"min_share\": 0.5 }`)",
            ));
        };
        const KNOWN: &[&str] =
            &["f", "n_b", "n_eff", "region", "min_share", "component", "min_component", "max_layer_fraction"];
        for k in map.keys() {
            if !KNOWN.contains(&k.as_str()) {
                return Err(Error::invalid(format!(
                    "unknown field '{k}' in `track` (expected `f`, `n_b`, `n_eff`, `region`, `min_share`, `component`, `min_component`, or `max_layer_fraction`)"
                )));
            }
        }
        let mut target = None;
        for k in ["f", "n_b", "n_eff"] {
            if let Some(v) = map.get(k).filter(|v| !v.is_null()) {
                if target.is_some() {
                    return Err(Error::invalid("`track` takes at most one target (`f`, `n_b` or `n_eff`)"));
                }
                let x = v
                    .as_f64()
                    .filter(|x| x.is_finite())
                    .ok_or_else(|| Error::invalid(format!("`track.{k}` must be a finite number")))?;
                target = Some(x);
            }
        }
        let opt_text = |k: &str| -> Result<Option<String>, Error> {
            match map.get(k).filter(|v| !v.is_null()) {
                None => Ok(None),
                Some(v) => v
                    .as_str()
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .map(Some)
                    .ok_or_else(|| Error::invalid(format!("`track.{k}` must be non-empty text"))),
            }
        };
        let opt_unit = |k: &str| -> Result<Option<f64>, Error> {
            match map.get(k).filter(|v| !v.is_null()) {
                None => Ok(None),
                Some(v) => v
                    .as_f64()
                    .filter(|x| x.is_finite() && (0.0..=1.0).contains(x))
                    .map(Some)
                    .ok_or_else(|| Error::invalid(format!("`track.{k}` must be a number between 0 and 1"))),
            }
        };
        let filter = TrackFilter {
            target,
            region: opt_text("region")?,
            min_share: opt_unit("min_share")?,
            component: opt_text("component")?,
            min_component: opt_unit("min_component")?,
            max_layer_fraction: opt_unit("max_layer_fraction")?,
        };
        filter.validate()?;
        Ok(filter)
    }

    pub fn to_json(&self) -> Value {
        let mut m = serde_json::Map::new();
        if let Some(t) = self.target {
            m.insert("f".into(), json!(t));
        }
        if let Some(r) = &self.region {
            m.insert("region".into(), json!(r));
        }
        if let Some(s) = self.min_share {
            m.insert("min_share".into(), json!(s));
        }
        if let Some(c) = &self.component {
            m.insert("component".into(), json!(c));
        }
        if let Some(mc) = self.min_component {
            m.insert("min_component".into(), json!(mc));
        }
        if let Some(mlf) = self.max_layer_fraction {
            m.insert("max_layer_fraction".into(), json!(mlf));
        }
        Value::Object(m)
    }

    fn describe(&self) -> String {
        let mut parts = vec![];
        if let Some(t) = self.target {
            parts.push(format!("target={t}"));
        }
        if let (Some(r), Some(s)) = (&self.region, self.min_share) {
            parts.push(format!("share in {r} >= {s}"));
        }
        if let (Some(c), Some(s)) = (&self.component, self.min_component) {
            let base = c.strip_suffix(" share").unwrap_or(c.as_str());
            parts.push(format!("{base} share >= {s}"));
        }
        if let Some(m) = self.max_layer_fraction {
            parts.push(format!("layer fraction <= {m}"));
        }
        parts.join(", ")
    }

    fn matches_mode(&self, r: &SolveResult, i: usize) -> bool {
        if key(r, i).is_none() {
            return false;
        }
        let n = i + 1;
        if let (Some(region), Some(min_share)) = (&self.region, self.min_share) {
            if !scalar_num(r, &format!("share in {region} {n}")).is_some_and(|v| v >= min_share) {
                return false;
            }
        }
        if let (Some(comp), Some(min_comp)) = (&self.component, self.min_component) {
            let base = comp.strip_suffix(" share").unwrap_or(comp.as_str());
            if !scalar_num(r, &format!("{base} share {n}")).is_some_and(|v| v >= min_comp) {
                return false;
            }
        }
        if let Some(max_lf) = self.max_layer_fraction {
            if !scalar_num(r, &format!("layer fraction {n}")).is_some_and(|v| v <= max_lf) {
                return false;
            }
        }
        true
    }

    fn candidates(&self, r: &SolveResult, mesh_label: &str) -> Result<Vec<usize>, Error> {
        if !self.has_signature() {
            return Ok((0..r.solutions.len()).collect());
        }
        let out: Vec<usize> = (0..r.solutions.len()).filter(|&i| self.matches_mode(r, i)).collect();
        if out.is_empty() {
            let available: Vec<&str> = r
                .scalars
                .iter()
                .filter_map(|(name, sc)| match sc {
                    Scalar::Num(_) if mode_scalar(name, r.solutions.len()).is_some_and(|(_, _, m)| m == 0) => {
                        Some(name.as_str())
                    }
                    _ => None,
                })
                .collect();
            let hint = if available.is_empty() {
                String::new()
            } else {
                format!(" (mode 1 scalars: {})", available.join(", "))
            };
            return Err(Error::invalid(format!(
                "no mode on the {mesh_label} mesh ({} solved) matched track filter [{}]{hint}",
                r.solutions.len(),
                self.describe()
            )));
        }
        Ok(out)
    }
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
    /// follow mode(s) matching a target eigenvalue and/or per-mode signature (`share in <region>`, `<component> share`, `layer fraction`)
    Track(TrackFilter),
}

impl Matching {
    /// From the arguments of the command: `match` is `index` or `nearest`, `track` filters by signature and/or target (`f`, `n_b`, `n_eff`) and implies `nearest`.
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
                let filter = TrackFilter::from_value(t)?;
                if named == Some(Matching::Index) {
                    return Err(Error::invalid(
                        "`track` follows a mode by its eigenvalue; it cannot be used with `match`: \"index\"",
                    ));
                }
                Ok(Matching::Track(filter))
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

/// Whether a row belongs to a noise mode: its `im`, or a scalar made of it (the quality factor `Q<n>`, `f im <n> (Hz)`,
/// `n_B im <n>`, `k_y im <n> (rad/m)`, `alpha <n> (dB/cm)`, `loss <n> (dB/cm)`).
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
    mode_scalar(rest, usize::MAX).is_some_and(|(head, _, mode)| {
        matches!(head.trim(), "f im" | "n_B im" | "k_y im" | "alpha" | "loss") && noise.contains(&mode)
    })
}

/// A pair of solutions: index on the present mesh, index on the finer mesh, distance of the eigenvalues.
type Pair = (usize, usize, f64);

/// The pairs of mutual nearest neighbours among `coarse_indices` and `fine_indices`, and the candidate indices left without a partner.
fn pair_nearest_in(
    coarse: &SolveResult,
    coarse_indices: &[usize],
    fine: &SolveResult,
    fine_indices: &[usize],
) -> (Vec<Pair>, Vec<usize>, Vec<usize>) {
    let nearest = |from: &SolveResult, i: usize, to: &SolveResult, to_indices: &[usize]| -> Option<(usize, f64)> {
        let a = key(from, i)?;
        to_indices
            .iter()
            .copied()
            .filter_map(|j| key(to, j).map(|b| (j, distance(a, b))))
            .min_by(|x, y| x.1.total_cmp(&y.1))
    };
    let mut pairs = vec![];
    for &i in coarse_indices {
        if let Some((j, d)) = nearest(coarse, i, fine, fine_indices) {
            if nearest(fine, j, coarse, coarse_indices).is_some_and(|(back, _)| back == i) {
                pairs.push((i, j, d));
            }
        }
    }
    let lone_coarse = coarse_indices.iter().copied().filter(|i| !pairs.iter().any(|p| p.0 == *i)).collect();
    let lone_fine = fine_indices.iter().copied().filter(|j| !pairs.iter().any(|p| p.1 == *j)).collect();
    (pairs, lone_coarse, lone_fine)
}

/// The pairs of mutual nearest neighbours, and the solutions of each mesh left without a partner.
fn pair_nearest(coarse: &SolveResult, fine: &SolveResult) -> (Vec<Pair>, Vec<usize>, Vec<usize>) {
    let ci: Vec<usize> = (0..coarse.solutions.len()).collect();
    let fi: Vec<usize> = (0..fine.solutions.len()).collect();
    pair_nearest_in(coarse, &ci, fine, &fi)
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
    /// cost caps, local size scaling, and Richardson extrapolation options (L3, FR-6, QFR-3)
    pub limits: Limits,
}

/// What a convergence check may cost and how it refines and extrapolates: it is refused, before any solve, when the comparison mesh has more cells
/// than `max_cells` or when its solve is estimated to take longer than `max_seconds`. `scale_sizes` (default `true`) scales non-suppressed `add_size`
/// features proportionally on the scratch copy; `extrapolate` optionally runs two-grid Richardson extrapolation (`"richardson_log_spd"` or
/// `"richardson_linear"`); `quantities` optionally restricts the verdict to matching row names.
#[derive(Clone, Debug, PartialEq)]
pub struct Limits {
    pub max_cells: Option<u64>,
    pub max_seconds: Option<f64>,
    pub scale_sizes: bool,
    pub extrapolate: Option<String>,
    pub extrapolate_order: Option<f64>,
    pub apply_extrapolated: bool,
    pub quantities: Option<Vec<String>>,
}

impl Default for Limits {
    fn default() -> Self {
        Limits {
            max_cells: None,
            max_seconds: None,
            scale_sizes: true,
            extrapolate: None,
            extrapolate_order: None,
            apply_extrapolated: false,
            quantities: None,
        }
    }
}

impl Limits {
    /// Validates `Limits` invariants whether parsed from JSON or constructed in Rust.
    pub fn validate(&self) -> Result<(), Error> {
        if let Some(m) = self.max_cells {
            if m == 0 {
                return Err(Error::invalid("`max_cells` is a number of cells, above zero"));
            }
        }
        if let Some(s) = self.max_seconds {
            if !(s > 0.0 && s.is_finite()) {
                return Err(Error::invalid("`max_seconds` is a number of seconds, above zero"));
            }
        }
        if let Some(ext) = &self.extrapolate {
            if !matches!(ext.as_str(), "richardson_log_spd" | "richardson_linear") {
                return Err(Error::invalid("`extrapolate` must be \"richardson_log_spd\" or \"richardson_linear\""));
            }
        }
        if let Some(p) = self.extrapolate_order {
            if !(p > 0.0 && p.is_finite()) {
                return Err(Error::invalid("`extrapolate_order` must be a finite number above zero"));
            }
        }
        if let Some(qs) = &self.quantities {
            if qs.is_empty() {
                return Err(Error::invalid("`quantities` must be a non-empty list of quantity names"));
            }
            if qs.iter().any(|q| q.trim().is_empty()) {
                return Err(Error::invalid("every entry in `quantities` must be non-empty text"));
            }
        }
        Ok(())
    }

    /// `max_cells`, `max_seconds`, `scale_sizes`, `extrapolate`, `extrapolate_order`, `apply_extrapolated` and `quantities` of the command arguments.
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
        let scale_sizes = match args.get("scale_sizes").filter(|v| !v.is_null()) {
            None => true,
            Some(v) => v.as_bool().ok_or_else(|| Error::invalid("`scale_sizes` is true or false"))?,
        };
        let extrapolate = match args.get("extrapolate").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    Error::invalid("`extrapolate` must be \"richardson_log_spd\" or \"richardson_linear\"")
                })?;
                if !matches!(s, "richardson_log_spd" | "richardson_linear") {
                    return Err(Error::invalid(
                        "`extrapolate` must be \"richardson_log_spd\" or \"richardson_linear\"",
                    ));
                }
                Some(s.to_string())
            }
        };
        let extrapolate_order = match args.get("extrapolate_order").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => Some(
                v.as_f64()
                    .filter(|p| *p > 0.0 && p.is_finite())
                    .ok_or_else(|| Error::invalid("`extrapolate_order` must be a finite number above zero"))?,
            ),
        };
        let apply_extrapolated = match args.get("apply_extrapolated").filter(|v| !v.is_null()) {
            None => false,
            Some(v) => v.as_bool().ok_or_else(|| Error::invalid("`apply_extrapolated` is true or false"))?,
        };
        let quantities = match args.get("quantities").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => {
                let arr = v
                    .as_array()
                    .ok_or_else(|| Error::invalid("`quantities` must be a non-empty list of quantity names"))?;
                if arr.is_empty() {
                    return Err(Error::invalid("`quantities` must be a non-empty list of quantity names"));
                }
                let mut out = Vec::with_capacity(arr.len());
                for item in arr {
                    let s = item
                        .as_str()
                        .map(str::trim)
                        .filter(|s| !s.is_empty())
                        .ok_or_else(|| Error::invalid("every entry in `quantities` must be non-empty text"))?;
                    out.push(s.to_string());
                }
                Some(out)
            }
        };
        let limits = Limits {
            max_cells,
            max_seconds,
            scale_sizes,
            extrapolate,
            extrapolate_order,
            apply_extrapolated,
            quantities,
        };
        limits.validate()?;
        Ok(limits)
    }
}

#[allow(clippy::needless_range_loop)]
fn invert_spd_ec_matrix(maxwell: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, Error> {
    let n = maxwell.len();
    if n == 0 {
        return Ok(vec![]);
    }
    let mut flat = Vec::with_capacity(n * n);
    for row in maxwell {
        flat.extend_from_slice(row);
    }
    let (eigs, vecs) = femcore::linalg::driver::dense_eigh(&flat, n).map_err(Error::from_kernel)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; n]; n];
    if lam_max <= 0.0 {
        return Ok(ec);
    }
    for k in 0..n {
        if eigs[k] > 1e-13 * lam_max {
            let inv_l = scale / eigs[k];
            for i in 0..n {
                for j in 0..n {
                    ec[i][j] += inv_l * vecs[k][i] * vecs[k][j];
                }
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let s = 0.5 * (ec[i][j] + ec[j][i]);
            ec[i][j] = s;
            ec[j][i] = s;
        }
    }
    Ok(ec)
}

#[allow(clippy::needless_range_loop)]
fn mutual_from_maxwell(maxwell: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = maxwell.len();
    let mut mutual = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut sum_row = 0.0;
        for k in 0..n {
            sum_row += maxwell[i][k];
            if i != k {
                mutual[i][k] = -maxwell[i][k];
            }
        }
        mutual[i][i] = sum_row;
    }
    mutual
}

/// Two-grid Richardson extrapolation of a multi-terminal Maxwell capacitance matrix (`r = h_coarse / h_fine > 1.0`, `order = p > 0.0`).
/// Returns `(extrapolated_data, per_entry_rel_error, max_entry_rel_error)`.
#[allow(clippy::needless_range_loop)]
pub fn extrapolate_capacitance(
    coarse: &CapacitanceMatrixData,
    fine: &CapacitanceMatrixData,
    method: &str,
    r: f64,
    order: f64,
) -> Result<(CapacitanceMatrixData, Vec<Vec<f64>>, f64), Error> {
    if !(r.is_finite() && r > 1.0) {
        return Err(Error::invalid("Richardson extrapolation requires element-size ratio `r > 1`"));
    }
    if !(order.is_finite() && order > 0.0) {
        return Err(Error::invalid("`extrapolate_order` must be a finite number above zero"));
    }
    if coarse.terminals != fine.terminals {
        return Err(Error::invalid("the coarse and fine capacitance matrices do not have the same terminal list"));
    }
    let n = fine.terminals.len();
    if coarse.maxwell_matrix_f.len() != n
        || fine.maxwell_matrix_f.len() != n
        || coarse.maxwell_matrix_f.iter().any(|row| row.len() != n)
        || fine.maxwell_matrix_f.iter().any(|row| row.len() != n)
    {
        return Err(Error::invalid("capacitance matrix dimensions do not match terminal count"));
    }
    let rp = r.powf(order);
    let w1 = rp / (rp - 1.0);
    let w0 = w1 - 1.0;
    let c_c = &coarse.maxwell_matrix_f;
    let c_f = &fine.maxwell_matrix_f;
    let mut c_r = vec![vec![0.0; n]; n];

    match method {
        "richardson_log_spd" => {
            for i in 0..n {
                for j in (i + 1)..n {
                    let af = 0.5 * (c_f[i][j].abs() + c_f[j][i].abs());
                    let ac = 0.5 * (c_c[i][j].abs() + c_c[j][i].abs());
                    let m_r = if af > 1e-30 && ac > 1e-30 { (w1 * af.ln() - w0 * ac.ln()).exp() } else { af };
                    c_r[i][j] = -m_r;
                    c_r[j][i] = -m_r;
                }
            }
            for i in 0..n {
                let sum_off_c: f64 = (0..n).filter(|&j| j != i).map(|j| c_c[i][j].abs()).sum();
                let self_inf_c = (c_c[i][i] - sum_off_c).max(0.0);
                let sum_off_f: f64 = (0..n).filter(|&j| j != i).map(|j| c_f[i][j].abs()).sum();
                let self_inf_f = (c_f[i][i] - sum_off_f).max(0.0);
                let self_inf_r = if self_inf_f > 1e-30 && self_inf_c > 1e-30 {
                    (w1 * self_inf_f.ln() - w0 * self_inf_c.ln()).exp()
                } else {
                    self_inf_f.max(0.0)
                };
                let sum_off_r: f64 = (0..n).filter(|&j| j != i).map(|j| c_r[i][j].abs()).sum();
                c_r[i][i] = sum_off_r + self_inf_r;
            }
        }
        "richardson_linear" => {
            for i in 0..n {
                for j in i..n {
                    let cf = 0.5 * (c_f[i][j] + c_f[j][i]);
                    let cc = 0.5 * (c_c[i][j] + c_c[j][i]);
                    let val = w1 * cf - w0 * cc;
                    c_r[i][j] = val;
                    c_r[j][i] = val;
                }
            }
        }
        other => {
            return Err(Error::invalid(format!(
                "`extrapolate` must be \"richardson_log_spd\" or \"richardson_linear\", not '{other}'"
            )));
        }
    }

    let mutual_matrix_f = mutual_from_maxwell(&c_r);
    let ec_matrix_hz = invert_spd_ec_matrix(&c_r)?;
    let mut per_entry_rel_error = vec![vec![0.0; n]; n];
    let mut max_entry_rel_error = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let scale = c_r[i][j].abs().max(c_f[i][j].abs());
            let rel = if scale < 1e-30 { 0.0 } else { (c_r[i][j] - c_f[i][j]).abs() / scale };
            per_entry_rel_error[i][j] = rel;
            max_entry_rel_error = max_entry_rel_error.max(rel);
        }
    }

    let data = CapacitanceMatrixData {
        terminals: fine.terminals.clone(),
        maxwell_matrix_f: c_r,
        mutual_matrix_f,
        ec_matrix_hz,
        ..Default::default()
    };
    Ok((data, per_entry_rel_error, max_entry_rel_error))
}

/// Replaces the `CapacitanceMatrixData` payload (and matching `C_*` / `E_C_*` scalars and solution values) in `result` with `cap`.
pub fn apply_extrapolated_capacitance(result: &mut SolveResult, cap: &CapacitanceMatrixData) {
    let encoded = cap.scalar();
    if let Some(slot) = result.scalars.iter_mut().find(|(k, _)| k == CAPACITANCE_KEY) {
        *slot = encoded;
    } else {
        result.scalars.push(encoded);
    }

    let n = cap.terminals.len();
    let set_num = |scalars: &mut Vec<(String, Scalar)>, key: &str, val: f64| {
        if let Some((_, s)) = scalars.iter_mut().find(|(k, _)| k == key) {
            *s = Scalar::Num(val);
        }
    };
    if n >= 1 && cap.maxwell_matrix_f.len() == n {
        set_num(&mut result.scalars, "stored energy W_e (J; per metre in 2D)", 0.5 * cap.maxwell_matrix_f[0][0]);
        if n == 1 {
            set_num(
                &mut result.scalars,
                "capacitance C = 2 W_e / dV^2 (F; per metre in 2D)",
                cap.maxwell_matrix_f[0][0],
            );
        }
        for i in 0..n {
            for j in 0..n {
                set_num(
                    &mut result.scalars,
                    &format!("C_{}_{} (F)", cap.terminals[i], cap.terminals[j]),
                    cap.maxwell_matrix_f[i][j],
                );
                set_num(&mut result.scalars, &format!("C_{i}_{j} (F)"), cap.maxwell_matrix_f[i][j]);
                if cap.ec_matrix_hz.len() == n && cap.ec_matrix_hz[i].len() == n {
                    set_num(
                        &mut result.scalars,
                        &format!("E_C_{}_{} (Hz)", cap.terminals[i], cap.terminals[j]),
                        cap.ec_matrix_hz[i][j],
                    );
                    set_num(&mut result.scalars, &format!("E_C_{i}_{j} (Hz)"), cap.ec_matrix_hz[i][j]);
                }
            }
        }
        if result.solutions.len() == n {
            for j in 0..n {
                result.solutions[j].value = Some(cap.maxwell_matrix_f[j][j]);
            }
        }
    }
}

/// Extracts the extrapolated `CapacitanceMatrixData` from a finished `convergence_check` JSON reply, if one was computed.
pub fn extrapolated_capacitance(answer: &Value) -> Option<CapacitanceMatrixData> {
    let v = answer.get("extrapolated")?.get("capacitance_matrix")?;
    if v.is_null() {
        return None;
    }
    serde_json::from_value(v.clone()).ok()
}

/// Whether a finished `convergence_check` JSON reply asked to promote its extrapolated values onto the stored result.
pub fn should_apply_extrapolated(answer: &Value) -> bool {
    answer.get("extrapolated").and_then(|e| e.get("apply_extrapolated")).and_then(Value::as_bool).unwrap_or(false)
}

fn default_extrapolate_order(fine_model: &MeshedModel, fine_result: &SolveResult) -> f64 {
    let has_cap = CapacitanceMatrixData::of(fine_result).is_some();
    let has_sheets = nyx_physics::TlsMeta::of(fine_result).is_some_and(|t| !t.sheet_labels.is_empty())
        || fine_model
            .evaluated()
            .list
            .features
            .iter()
            .any(|f| !f.suppressed && matches!(f.op, geom::model::Op::Sheet { .. }));
    if has_cap && (fine_model.mesh().dim() == 3 || has_sheets) {
        1.0
    } else {
        (2 * fine_result.order.max(1) as usize) as f64
    }
}

impl ConvergenceJob {
    /// Returns this job configured with `limits`.
    pub fn with_limits(mut self, limits: Limits) -> Self {
        self.limits = limits;
        self
    }

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
        self.plan_with(physics, factor, self.limits.scale_sizes)
    }

    /// `plan` with an explicit `scale_sizes` toggle controlling whether local `add_size` features scale proportionally with `factor`.
    pub fn plan_with(&self, physics: &str, factor: f64, scale_sizes: bool) -> Result<Value, Error> {
        check_factor(factor)?;
        self.project
            .spec(physics)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{physics}'")))?;
        let targets = vec![physics.to_string()];
        let sizing_order = self.scratch(self.project.clone()).sizing().map(|s| s.order);
        let mut base = self.project.clone();
        force_single_adapt_pass(&mut base);
        let other_proj = refined(&base, factor, scale_sizes, sizing_order);
        let present = self.mesh_inline(&base)?;
        let other = self.mesh_inline(&other_proj)?;
        let (coarse_m, fine_m) = if factor < 1.0 { (&present, &other) } else { (&other, &present) };
        let (c, f) = (coarse_m.mesh().n_cells() as u64, fine_m.mesh().n_cells() as u64);
        let (dc, df) = (coarse_m.mesh().dim() as u8, fine_m.mesh().dim() as u8);
        let (cost_c, cost_f) = (
            crate::preview::cost_of_project(&self.history, &base, &targets, dc, c),
            crate::preview::cost_of_project(&self.history, &other_proj, &targets, df, f),
        );
        let secs = |v: &Value| v["solve"]["seconds"].as_f64();
        let total = secs(&cost_c).zip(secs(&cost_f)).map(|(a, b)| a + b);
        let p_cells = present.mesh().n_cells() as u64;
        let o_cells = other.mesh().n_cells() as u64;
        Ok(json!({
            "physics": physics, "factor": factor, "scale_sizes": scale_sizes, "dry_run": true,
            "cells": { "coarse": c, "refined": f },
            "refined_over_present": o_cells as f64 / p_cells.max(1) as f64,
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
    /// when given, is the largest relative difference accepted (the difference of a quantity over the larger of the two values, or the absolute
    /// difference for unit-fraction mode shares); every row then says whether it passed and `converged` is true when all did.
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

fn wildcard_match(text: &str, pat: &str) -> bool {
    if !pat.contains('*') {
        return text == pat;
    }
    let parts: Vec<&str> = pat.split('*').collect();
    let mut rest = text;
    for (idx, part) in parts.iter().enumerate() {
        if idx == 0 {
            let Some(after) = rest.strip_prefix(part) else { return false };
            rest = after;
        } else if idx == parts.len() - 1 {
            return rest.ends_with(part);
        } else {
            let Some(pos) = rest.find(part) else { return false };
            rest = &rest[pos + part.len()..];
        }
    }
    true
}

fn quantity_matches(row_name: &str, filter: &str) -> bool {
    let f = filter.trim();
    wildcard_match(row_name, f) || row_name.strip_prefix("scalar ").is_some_and(|raw| wildcard_match(raw, f))
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
        check_factor(factor)?;
        this.limits.validate()?;
        this.project
            .spec(physics)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{physics}'")))?;
        let targets = vec![physics.to_string()];
        let t0 = std::time::Instant::now();
        // L3: the comparison mesh is built first (milliseconds against a solve), so a check larger than the caller accepts is refused before
        // any solve, and the answer to a refusal says what the check would have cost
        let scale_sizes = this.limits.scale_sizes;
        let sizing_order = this.scratch(this.project.clone()).sizing().map(|s| s.order);
        let had_amr = project_has_multi_pass_amr(&this.project);
        let mut base_project = this.project.clone();
        force_single_adapt_pass(&mut base_project);
        let other_label = if factor < 1.0 { "finer" } else { "coarser" };
        let other_project = refined(&base_project, factor, scale_sizes, sizing_order);
        let other_model = {
            let _stage = femcore::progress::stage(format!(
                "convergence check: the {other_label} mesh (element size times {factor})"
            ));
            this.mesh_inline(&other_project)?
        };
        let (other_cells, other_dim) = (other_model.mesh().n_cells() as u64, other_model.mesh().dim() as u8);
        if let Some(max) = this.limits.max_cells.filter(|m| other_cells > *m) {
            return Err(Error::invalid(format!(
                "the {other_label} mesh would have {other_cells} cells, more than max_cells {max}: the check was not started. Coarsen the present mesh first (a larger h or adapt_scale), use a factor nearer 1 (now {factor}), or raise max_cells"
            )));
        }
        if let Some(max) = this.limits.max_seconds {
            let cost = crate::preview::cost_of_project(&this.history, &other_project, &targets, other_dim, other_cells);
            if let Some(s) = cost["solve"]["seconds"].as_f64().filter(|s| *s > max) {
                return Err(Error::invalid(format!(
                    "the solve on the {other_label} mesh ({other_cells} cells) is expected to take about {s:.0} s from the solves this machine has run, more than max_seconds {max}: the check was not started. Coarsen the mesh, use a factor nearer 1 (now {factor}), or raise max_seconds"
                )));
            }
        }
        // the two solves are stages of the one job, run here under its progress, so the check can be watched and cancelled; the first is skipped
        // when an earlier request solved this physics alone on this mesh (found by content, so the numbers are the ones this solve would give).
        // Mode pairing, Richardson extrapolation, and multi-pass AMR projects need fresh single-pass `SolveResult` payloads.
        let need_full_present = had_amr || !matches!(matching, Matching::Index) || this.limits.extrapolate.is_some();
        let reusable =
            if !need_full_present { this.reuse.as_ref().and_then(|(stages, key)| stages.solve(key)) } else { None };
        let (present_summary, self_reused, present_solved) = match reusable {
            Some(s) => {
                let _stage = femcore::progress::stage(
                    "convergence check: the solve on the present mesh (reused from an earlier solve)",
                );
                (s, true, None)
            }
            None => {
                let _stage = femcore::progress::stage("convergence check: the solve on the present mesh");
                let (model, r) = this.solve_inline(base_project, &targets, physics, progress)?;
                (Arc::new(SolveSummary::of(&r, model.mesh().n_cells())), false, Some((model, r)))
            }
        };
        let other_result = {
            let _stage = femcore::progress::stage(format!(
                "convergence check: the solve on the {other_label} mesh (element size times {factor})"
            ));
            this.solve_on(other_project, &other_model, &targets, physics, progress)?
        };
        let other_summary = Arc::new(SolveSummary::of(&other_result, other_model.mesh().n_cells()));

        // Identify `(coarse, fine)` regardless of whether `factor < 1` (refining) or `factor > 1` (coarsening)
        let (coarse_sum, fine_sum, coarse_pair, fine_pair) = if factor < 1.0 {
            (present_summary, other_summary, present_solved, Some((other_model, other_result)))
        } else {
            (other_summary, present_summary, Some((other_model, other_result)), present_solved)
        };
        let cells_coarse = Some(coarse_sum.cells);
        let cells_fine = Some(fine_sum.cells);

        // the rows to compare: (name, coarse, refined) and, when the modes are paired, which solutions were paired
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
            Matching::Index => common(coarse_sum.numbers.clone(), fine_sum.numbers.clone()),
            Matching::Nearest | Matching::Track(_) => {
                let (_, coarse_r) =
                    coarse_pair.as_ref().expect("a check that pairs modes solves both meshes with full results");
                let (_, fine) =
                    fine_pair.as_ref().expect("a check that pairs modes solves both meshes with full results");
                let (pairs, lone_coarse, lone_fine) = match matching {
                    Matching::Track(filter) => {
                        filter.validate()?;
                        let coarse_cands = filter.candidates(coarse_r, "present")?;
                        let fine_cands = filter.candidates(fine, "refined")?;
                        if let Some(f) = filter.target {
                            let i = coarse_cands
                                .iter()
                                .copied()
                                .filter_map(|i| key(coarse_r, i).map(|k| (i, (k[0] - f).abs())))
                                .min_by(|x, y| x.1.total_cmp(&y.1))
                                .map(|(i, _)| i)
                                .ok_or_else(|| Error::invalid("the result has no eigenvalue or value to track"))?;
                            let a = key(coarse_r, i).expect("chosen by its key");
                            let (j, d) = fine_cands
                                .iter()
                                .copied()
                                .filter_map(|j| key(fine, j).map(|b| (j, distance(a, b))))
                                .min_by(|x, y| x.1.total_cmp(&y.1))
                                .ok_or_else(|| {
                                    Error::invalid("the refined result has no eigenvalue or value to track")
                                })?;
                            (vec![(i, j, d)], vec![], vec![])
                        } else {
                            pair_nearest_in(coarse_r, &coarse_cands, fine, &fine_cands)
                        }
                    }
                    _ => pair_nearest(coarse_r, fine),
                };
                for &(i, j, _) in &pairs {
                    common(numbers_of(coarse_r, i, i), numbers_of(fine, j, i));
                }
                matched = Some((pairs, lone_coarse, lone_fine));
            }
        }
        // L2: a mode whose imaginary part is zero up to rounding (a lossless guide solved by the complex solver, because a crystal axis is tilted)
        // has `im` of 1e-10 or so that differs between two meshes by 100 percent, and a quality factor of 1e10 made of it. Such a mode is not a
        // failure of the mesh: its `im` and its quality factors are listed but not counted. The floor is what the eigensolver resolves: an `im`
        // below it times the real part on both meshes (a Q above 5e7) is not told from rounding.
        let noise = noise_modes(&compared);
        let active_compared: Vec<&(String, f64, f64)> = match &this.limits.quantities {
            None => compared.iter().collect(),
            Some(qs) => {
                let matched_rows: Vec<&(String, f64, f64)> =
                    compared.iter().filter(|(name, _, _)| qs.iter().any(|q| quantity_matches(name, q))).collect();
                if matched_rows.is_empty() {
                    let available: Vec<&str> = compared.iter().map(|(n, _, _)| n.as_str()).collect();
                    return Err(Error::invalid(format!(
                        "no quantity matched `quantities` {:?} (available quantities: {})",
                        qs,
                        if available.is_empty() { "none".to_string() } else { available.join(", ") }
                    )));
                }
                matched_rows
            }
        };
        let mut rows = vec![];
        let mut worst = 0.0f64;
        let mut worst_name = String::new();
        let mut all = true;
        for &(ref name, x, y) in &active_compared {
            let is_unit_fraction = name
                .strip_prefix("scalar ")
                .and_then(|rest| mode_scalar(rest, usize::MAX))
                .is_some_and(|(head, _, _)| {
                    let h = head.trim();
                    h.ends_with("share") || h.starts_with("share in ") || h == "layer fraction"
                });
            let scale = x.abs().max(y.abs());
            let rel = if is_unit_fraction {
                (x - y).abs()
            } else if scale <= 1e-30 {
                0.0
            } else {
                (x - y).abs() / scale
            };
            if is_noise(name, &noise) {
                let mut row = json!({ "quantity": name, "coarse": x, "refined": y, "difference": y - x, "relative": rel, "converged": tolerance.map(|_| true), "noise": true });
                row["why"] = json!("imaginary part at the rounding level of this mode: not counted");
                rows.push(row);
                continue;
            }
            if matches!(name.as_str(), "scalar max |E| (V/m)" | "scalar max |B| (T)") {
                let mut row = json!({ "quantity": name, "coarse": x, "refined": y, "difference": y - x, "relative": rel, "converged": tolerance.map(|_| true), "noise": true });
                row["why"] =
                    json!("pointwise field peak at conductor edges/corners is singular under refinement: not counted");
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
            "physics": physics, "factor": factor, "scale_sizes": scale_sizes, "tolerance": tolerance,
            "cells": { "coarse": cells_coarse, "refined": cells_fine },
            "unknowns": { "coarse": coarse_sum.dofs, "refined": fine_sum.dofs },
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
        if let Some(qs) = &this.limits.quantities {
            out["quantities_filter"] = json!(qs);
        }
        if let Some((pairs, lone_coarse, lone_fine)) = matched {
            let mode = if matches!(matching, Matching::Track(_)) { "track" } else { "nearest" };
            out["matching"] = json!(mode);
            out["match"] = json!(mode);
            if let Matching::Track(filter) = matching {
                out["track"] = filter.to_json();
            }
            out["pairs"] = pairs.iter().map(|(i, j, d)| json!({ "coarse": i, "refined": j, "distance": d })).collect();
            out["unpaired"] = json!({ "coarse": lone_coarse, "refined": lone_fine });
            out["note"] = json!("the project was not changed; modes are paired by their eigenvalues: a quantity carries the index of the mode on the present mesh and `pairs` gives its index on the finer one; rows, worst_relative and converged cover the paired modes only");
        }
        if let Some(method) = &this.limits.extrapolate {
            let (_, coarse_r) =
                coarse_pair.as_ref().expect("Richardson extrapolation solves both meshes with full results");
            let (fine_m, fine_r) =
                fine_pair.as_ref().expect("Richardson extrapolation solves both meshes with full results");
            let r = if factor < 1.0 { 1.0 / factor } else { factor };
            let order = this.limits.extrapolate_order.unwrap_or_else(|| default_extrapolate_order(fine_m, fine_r));
            let rp = r.powf(order);
            let w1 = rp / (rp - 1.0);
            let w0 = w1 - 1.0;

            let (cap_json, per_entry_rel, max_entry_rel) = match (
                CapacitanceMatrixData::of(coarse_r),
                CapacitanceMatrixData::of(fine_r),
            ) {
                (Some(cc), Some(cf)) => {
                    let (cap_r, per_rel, max_rel) = extrapolate_capacitance(&cc, &cf, method, r, order)?;
                    (serde_json::to_value(&cap_r).unwrap_or(Value::Null), json!(per_rel), json!(max_rel))
                }
                _ if method == "richardson_log_spd" => {
                    return Err(Error::invalid(
                            "`extrapolate: \"richardson_log_spd\"` requires a multi-terminal capacitance matrix (`dfs` electrostatic study with `terminal` conditions)",
                        ));
                }
                _ => (Value::Null, Value::Null, Value::Null),
            };

            let mut ext_scalars = serde_json::Map::new();
            for (name, x_c, x_f) in &compared {
                if is_noise(name, &noise) || matches!(name.as_str(), "scalar max |E| (V/m)" | "scalar max |B| (T)") {
                    continue;
                }
                let val = if method == "richardson_log_spd" && *x_f > 1e-30 && *x_c > 1e-30 {
                    (w1 * x_f.ln() - w0 * x_c.ln()).exp()
                } else if method == "richardson_log_spd" && *x_f < -1e-30 && *x_c < -1e-30 {
                    -(w1 * x_f.abs().ln() - w0 * x_c.abs().ln()).exp()
                } else {
                    w1 * x_f - w0 * x_c
                };
                if val.is_finite() {
                    ext_scalars.insert(name.clone(), json!(val));
                }
            }

            out["extrapolated"] = json!({
                "method": method,
                "order": order,
                "r": r,
                "weight_fine": w1,
                "apply_extrapolated": this.limits.apply_extrapolated,
                "capacitance_matrix": cap_json,
                "per_entry_rel_error": per_entry_rel,
                "max_entry_rel_error": max_entry_rel,
                "scalars": ext_scalars,
            });
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
        assert!(
            is_noise("scalar n_B im 1", &n)
                && is_noise("scalar k_y im 1 (rad/m)", &n)
                && is_noise("scalar alpha 1 (dB/cm)", &n)
                && !is_noise("scalar n_B re 1", &n)
                && !is_noise("scalar n_g 1", &n)
        );
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

    #[test]
    fn refined_scales_pml_layers_when_scale_sizes_is_on_and_preserves_auto_sentinel() {
        use geom::expr::Expr;
        use geom::model::{Feature, Op, Origin, PmlCoordKind, PmlGrowMode};

        let mut p = Project::new("pml_conv", 2);
        p.geometry.features.push(Feature {
            id: 1,
            name: "pml_explicit".into(),
            suppressed: false,
            origin: Origin::Code,
            tags: vec![],
            op: Op::Pml {
                window: None,
                sides: vec![],
                thickness: Expr::Num(1e-6),
                layers: Expr::Num(6.0),
                grading: Expr::Num(2.0),
                regions: vec![],
                mode: PmlGrowMode::Outward,
                coord: PmlCoordKind::Auto,
                center: None,
                axis: None,
                layer_grading: Some(Expr::Num(2.0)),
            },
        });
        p.geometry.features.push(Feature {
            id: 2,
            name: "pml_auto".into(),
            suppressed: false,
            origin: Origin::Code,
            tags: vec![],
            op: Op::Pml {
                window: None,
                sides: vec![],
                thickness: Expr::Text("auto".into()),
                layers: Expr::Num(-1.0),
                grading: Expr::Num(2.0),
                regions: vec![],
                mode: PmlGrowMode::Outward,
                coord: PmlCoordKind::Auto,
                center: None,
                axis: None,
                layer_grading: Some(Expr::Num(2.0)),
            },
        });
        let r = refined(&p, 0.5, true, Some(2));
        match &r.geometry.features[0].op {
            Op::Pml { layers, .. } => assert_eq!(*layers, Expr::Num(12.0)),
            other => panic!("unexpected op: {other:?}"),
        }
        match &r.geometry.features[1].op {
            Op::Pml { layers, .. } => assert_eq!(*layers, Expr::Num(-1.0)),
            other => panic!("unexpected op: {other:?}"),
        }
        let r_coarse = refined(&p, 2.0, true, Some(2));
        match &r_coarse.geometry.features[0].op {
            Op::Pml { layers, .. } => assert_eq!(*layers, Expr::Num(5.0)),
            other => panic!("unexpected op: {other:?}"),
        }
    }
}
