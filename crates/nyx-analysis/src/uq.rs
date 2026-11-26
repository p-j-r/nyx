// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Zero-solve post-processing uncertainty quantification for cQED observables (`UqCqedAnalysis`, §6.4).
//!
//! Propagates fabrication and constitutive `1 sigma` tolerances (interface thickness/permittivity/`tan_delta`,
//! film thickness `d`, trench depth, `E_J`/`L_J`/`C_J`/`E_L`, SQUID asymmetry and flux bias, resonator `f_r`/`Z_0`,
//! and 3D finite-thickness `metal_thickness_m`) using:
//! - closed-form / central finite-difference linearized covariance (`linear_sigma` and `variance_contributions`), and
//! - deterministic `SplitMix64` Box–Muller Monte Carlo (`mc_mean`, `mc_sigma`, `ci_68`, `ci_95`),
//!
//! all with **zero additional FEM solves**.
//!
//! **Strict ASME V&V 20 Separation (§6.1):** Deterministic discretization brackets (Prager–Synge `[low, up]`
//! capacitance matrices and `SolveResult::error_budget.spatial_discretization` bounds) are stored strictly in
//! [`EpistemicInterval`] and never folded into aleatoric Gaussian `sigma` estimates. When both are present,
//! an outer [`PBoxEnvelope`] (`outer_bracket_95`, `outer_bracket_68`) is reported alongside the two separate components.
#![allow(clippy::needless_range_loop)]

use crate::capacitance::{
    apply_3d_thickness_gram, capacitance_matrix, capacitance_matrix_applies, prepare_3d_thickness_gram,
    run_capacitance_matrix_from_thickness, PreparedThicknessGram,
};
use crate::epr::{
    epr_applies, epr_with_params, is_nonlinear_port, jj_flux_scale, resolve_jj_l_and_ej, squid_effective_ej, EprReport,
    C0, EPS0, MU0,
};
use crate::hamiltonian::{hamiltonian_applies, hamiltonian_with_params, HamiltonianReport};
use crate::tls::{
    prepare_quasi_electrostatic_tls, prepare_tls, tls_participation_applies, PreparedTls, TlsParams, TlsReport,
};
use crate::Analysis;
use femcore::linalg::rng::SplitMix64;
use femcore::par::map_chunks;
use nyx_physics::{CapacitanceMatrixData, CircuitRom, Error, Mesh, SolveResult, TlsMeta};
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::time::Instant;

const DEFAULT_SEED: u64 = 0x4E59_585F_5551_0004;
const DEFAULT_N_SAMPLES: usize = 512;
const MAX_N_SAMPLES: usize = 100_000;

/// Deterministic Box–Muller standard-normal generator seeded by [`SplitMix64`].
struct NormalRng {
    rng: SplitMix64,
    spare: Option<f64>,
}

impl NormalRng {
    fn new(seed: u64) -> Self {
        Self { rng: SplitMix64::new(seed), spare: None }
    }

    #[inline]
    fn uniform_open01(&mut self) -> f64 {
        // Map 53 top bits into (0, 1] for strictly positive log argument.
        let bits = self.rng.next_u64() >> 11;
        ((bits as f64) + 0.5) * (1.0 / ((1u64 << 53) as f64))
    }

    fn next_normal(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        let u1 = self.uniform_open01().max(1.0e-300);
        let u2 = self.uniform_open01();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * PI * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
}

/// Epistemic interval bracket `[low, up]` (from Prager–Synge, Roache GCI, or solver residuals).
///
/// Never folded into Gaussian `sigma` under ASME V&V 20 strict epistemic/aleatoric separation.
#[derive(Clone, Debug, PartialEq)]
pub struct EpistemicInterval {
    pub low: f64,
    pub up: f64,
    pub half_width: f64,
    pub midpoint: f64,
    pub rel_half_width: f64,
    pub source: String,
    pub is_strict_bound: bool,
}

impl EpistemicInterval {
    pub fn from_bounds(low: f64, up: f64, nominal: f64, source: impl Into<String>, is_strict_bound: bool) -> Self {
        let (lo, hi) = if low <= up { (low, up) } else { (up, low) };
        let half_width = 0.5 * (hi - lo);
        let midpoint = 0.5 * (lo + hi);
        let denom = nominal.abs().max(midpoint.abs()).max(1e-300);
        Self {
            low: lo,
            up: hi,
            half_width,
            midpoint,
            rel_half_width: half_width / denom,
            source: source.into(),
            is_strict_bound,
        }
    }

    pub fn from_rel_bound(
        nominal: f64,
        rel_bound: f64,
        non_negative: bool,
        source: String,
        is_strict_bound: bool,
    ) -> Self {
        let delta = nominal.abs() * rel_bound.abs();
        let lo = if non_negative && nominal >= 0.0 { (nominal - delta).max(0.0) } else { nominal - delta };
        let hi = nominal + delta;
        Self::from_bounds(lo, hi, nominal, source, is_strict_bound)
    }

    pub fn to_json(&self) -> Value {
        json!({
            "low": self.low,
            "up": self.up,
            "half_width": self.half_width,
            "midpoint": self.midpoint,
            "rel_half_width": self.rel_half_width,
            "source": self.source,
            "is_strict_bound": self.is_strict_bound,
        })
    }
}

/// Aleatoric statistical uncertainty summary (`1 sigma` linearized + Monte Carlo quantiles).
#[derive(Clone, Debug, PartialEq)]
pub struct AleatoricEstimate {
    pub linear_sigma: f64,
    pub rel_linear_sigma: f64,
    pub mc_mean: f64,
    pub mc_sigma: f64,
    pub rel_mc_sigma: f64,
    pub ci_68: [f64; 2],
    pub ci_95: [f64; 2],
    pub variance_contributions: BTreeMap<String, f64>,
    pub sensitivities: BTreeMap<String, f64>,
}

impl AleatoricEstimate {
    pub fn to_json(&self) -> Value {
        json!({
            "linear_sigma": self.linear_sigma,
            "rel_linear_sigma": self.rel_linear_sigma,
            "mc_mean": self.mc_mean,
            "mc_sigma": self.mc_sigma,
            "rel_mc_sigma": self.rel_mc_sigma,
            "ci_68": self.ci_68,
            "ci_95": self.ci_95,
            "variance_contributions": self.variance_contributions,
            "sensitivities": self.sensitivities,
        })
    }
}

/// Outer p-box envelope `[lower, upper]` combining separate aleatoric MC intervals with epistemic brackets.
#[derive(Clone, Debug, PartialEq)]
pub struct PBoxEnvelope {
    pub bracket_95: [f64; 2],
    pub bracket_68: [f64; 2],
}

impl PBoxEnvelope {
    pub fn to_json(&self) -> Value {
        json!({
            "bracket_95": self.bracket_95,
            "bracket_68": self.bracket_68,
        })
    }
}

/// Metadata describing an active uncertain input parameter.
#[derive(Clone, Debug, PartialEq)]
pub struct UqParameterSpec {
    pub name: String,
    pub domain: String,
    pub nominal: f64,
    pub sigma: f64,
    pub rel_tol: f64,
    pub min_bound: f64,
    pub max_bound: f64,
}

impl UqParameterSpec {
    pub fn to_json(&self) -> Value {
        json!({
            "name": self.name,
            "domain": self.domain,
            "nominal": self.nominal,
            "sigma": self.sigma,
            "rel_tol": self.rel_tol,
            "min_bound": self.min_bound.is_finite().then_some(self.min_bound),
            "max_bound": self.max_bound.is_finite().then_some(self.max_bound),
        })
    }
}

/// Uncertainty-quantified Quantity of Interest (QoI) under strict ASME V&V 20 separation.
#[derive(Clone, Debug, PartialEq)]
pub struct UqQuantity {
    pub name: String,
    pub unit: String,
    pub nominal: f64,
    pub aleatoric: Option<AleatoricEstimate>,
    pub epistemic: Option<EpistemicInterval>,
    pub outer_envelope: Option<PBoxEnvelope>,
}

impl UqQuantity {
    pub fn to_json(&self) -> Value {
        let mut obj = Map::new();
        obj.insert("name".into(), json!(self.name));
        obj.insert("unit".into(), json!(self.unit));
        obj.insert("nominal".into(), json!(self.nominal));
        if let Some(ref a) = self.aleatoric {
            obj.insert("aleatoric".into(), a.to_json());
        }
        if let Some(ref e) = self.epistemic {
            obj.insert("epistemic".into(), e.to_json());
        }
        if let Some(ref p) = self.outer_envelope {
            obj.insert("outer_envelope".into(), p.to_json());
            obj.insert("outer_bracket_95".into(), json!(p.bracket_95));
            obj.insert("outer_bracket_68".into(), json!(p.bracket_68));
        }
        Value::Object(obj)
    }
}

/// Typed output of [`uq_cqed`] and [`uq_cqed_with_params`].
#[derive(Clone, Debug, PartialEq)]
pub struct UqCqedReport {
    pub seed: u64,
    pub n_samples: usize,
    pub zero_fem_solves: bool,
    pub epistemic_aleatoric_separated: bool,
    pub domains_evaluated: Vec<String>,
    pub uncertain_parameters: Vec<UqParameterSpec>,
    pub quantities: Vec<UqQuantity>,
    pub diagnostics: Vec<String>,
    pub elapsed_s: f64,
}

impl UqCqedReport {
    /// Looks up a quantified observable by exact name (e.g. `"f01_Hz"`, `"Q_TLS"`, `"C_Sigma_F"`).
    pub fn quantity(&self, name: &str) -> Option<&UqQuantity> {
        self.quantities.iter().find(|q| q.name == name)
    }

    pub fn to_json(&self) -> Value {
        let mut q_map = Map::new();
        for q in &self.quantities {
            q_map.insert(q.name.clone(), q.to_json());
        }
        json!({
            "asme_vv20_separated": self.epistemic_aleatoric_separated,
            "zero_fem_solves": self.zero_fem_solves,
            "seed": self.seed,
            "n_samples": self.n_samples,
            "domains_evaluated": self.domains_evaluated,
            "uncertain_parameters": self.uncertain_parameters.iter().map(UqParameterSpec::to_json).collect::<Vec<_>>(),
            "quantities": Value::Object(q_map),
            "diagnostics": self.diagnostics,
            "elapsed_s": self.elapsed_s,
        })
    }
}

/// True when `result` carries at least one cQED payload (`CircuitRom`, `CapacitanceMatrixData`, or `TlsMeta`/E-field).
pub fn uq_cqed_applies(result: &SolveResult) -> bool {
    epr_applies(result)
        || hamiltonian_applies(result)
        || capacitance_matrix_applies(result)
        || tls_participation_applies(result)
}

#[derive(Clone, Debug)]
struct RawUncSpec {
    sigma: Option<f64>,
    rel_tol: Option<f64>,
    min_override: Option<f64>,
    max_override: Option<f64>,
}

#[derive(Clone, Debug)]
struct ResolvedVar {
    spec: UqParameterSpec,
    step_h: f64,
}

#[derive(Clone, Debug)]
struct EvalScalarDef {
    name: String,
    unit: String,
    non_negative: bool,
}

fn parse_raw_uncertainties(params: &Value) -> Result<BTreeMap<String, RawUncSpec>, Error> {
    let mut out = BTreeMap::new();
    let Some(unc_val) = params.get("uncertainties").filter(|v| !v.is_null()) else {
        return Ok(out);
    };
    let obj = unc_val.as_object().ok_or_else(|| Error::invalid("parameter 'uncertainties' must be a JSON object"))?;
    for (key, val) in obj {
        let spec = match val {
            Value::Number(n) => {
                let s = n.as_f64().filter(|x| x.is_finite() && *x >= 0.0).ok_or_else(|| {
                    Error::invalid(format!("uncertainties['{key}'] must be a non-negative finite number"))
                })?;
                RawUncSpec { sigma: Some(s), rel_tol: None, min_override: None, max_override: None }
            }
            Value::Object(m) => {
                for fk in m.keys() {
                    if !matches!(fk.as_str(), "sigma" | "rel_tol" | "min" | "max") {
                        return Err(Error::invalid(format!(
                            "unknown field '{fk}' in uncertainties['{key}'] (expected: sigma, rel_tol, min, max)"
                        )));
                    }
                }
                let sigma = match m.get("sigma").filter(|v| !v.is_null()) {
                    Some(v) => Some(v.as_f64().filter(|x| x.is_finite() && *x >= 0.0).ok_or_else(|| {
                        Error::invalid(format!("uncertainties['{key}'].sigma must be a non-negative finite number"))
                    })?),
                    None => None,
                };
                let rel_tol = match m.get("rel_tol").filter(|v| !v.is_null()) {
                    Some(v) => Some(v.as_f64().filter(|x| x.is_finite() && *x >= 0.0).ok_or_else(|| {
                        Error::invalid(format!("uncertainties['{key}'].rel_tol must be a non-negative finite number"))
                    })?),
                    None => None,
                };
                if sigma.is_none() && rel_tol.is_none() {
                    return Err(Error::invalid(format!(
                        "uncertainties['{key}'] must specify at least 'sigma' or 'rel_tol'"
                    )));
                }
                let min_override = match m.get("min").filter(|v| !v.is_null()) {
                    Some(v) => Some(
                        v.as_f64()
                            .filter(|x| x.is_finite())
                            .ok_or_else(|| Error::invalid(format!("uncertainties['{key}'].min must be finite")))?,
                    ),
                    None => None,
                };
                let max_override = match m.get("max").filter(|v| !v.is_null()) {
                    Some(v) => Some(
                        v.as_f64()
                            .filter(|x| x.is_finite())
                            .ok_or_else(|| Error::invalid(format!("uncertainties['{key}'].max must be finite")))?,
                    ),
                    None => None,
                };
                if let (Some(lo), Some(hi)) = (min_override, max_override) {
                    if lo > hi {
                        return Err(Error::invalid(format!("uncertainties['{key}'] has min ({lo}) > max ({hi})")));
                    }
                }
                RawUncSpec { sigma, rel_tol, min_override, max_override }
            }
            _ => {
                return Err(Error::invalid(format!(
                    "uncertainties['{key}'] must be a non-negative number or {{sigma, rel_tol, min, max}} object"
                )))
            }
        };
        out.insert(key.clone(), spec);
    }
    Ok(out)
}

fn resolve_var(
    raw: &BTreeMap<String, RawUncSpec>,
    raw_key: &str,
    out_name: impl Into<String>,
    domain: &str,
    nominal: f64,
    default_min: f64,
    default_max: f64,
) -> Result<Option<ResolvedVar>, Error> {
    let Some(spec) = raw.get(raw_key) else {
        return Ok(None);
    };
    let sigma = match (spec.sigma, spec.rel_tol) {
        (Some(s), Some(r)) => s.max(nominal.abs() * r),
        (Some(s), None) => s,
        (None, Some(r)) => {
            if nominal.abs() == 0.0 && r > 0.0 {
                return Err(Error::invalid(format!(
                    "uncertainties['{raw_key}'].rel_tol cannot be applied when nominal value is 0; supply explicit 'sigma'"
                )));
            }
            nominal.abs() * r
        }
        (None, None) => 0.0,
    };
    if !sigma.is_finite() || sigma <= 0.0 {
        return Ok(None);
    }
    let min_bound = spec.min_override.map_or(default_min, |m| m.max(default_min));
    let max_bound = spec.max_override.map_or(default_max, |m| m.min(default_max));
    if min_bound > max_bound {
        return Err(Error::invalid(format!(
            "effective bounds for uncertain parameter '{raw_key}' are inverted ({min_bound} > {max_bound})"
        )));
    }
    let nominal_clamped = nominal.clamp(min_bound, max_bound);
    let step_h = (1.0e-5 * nominal_clamped.abs()).max(1.0e-5 * sigma).max(1.0e-12);
    Ok(Some(ResolvedVar {
        spec: UqParameterSpec {
            name: out_name.into(),
            domain: domain.to_string(),
            nominal: nominal_clamped,
            sigma,
            rel_tol: if nominal_clamped.abs() > 0.0 { sigma / nominal_clamped.abs() } else { 0.0 },
            min_bound,
            max_bound,
        },
        step_h,
    }))
}

fn quantile_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let pos = q.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    let t = pos - (lo as f64);
    sorted[lo] * (1.0 - t) + sorted[hi] * t
}

fn outer_pbox(
    nominal: f64,
    aleatoric: Option<&AleatoricEstimate>,
    epistemic: Option<&EpistemicInterval>,
) -> Option<PBoxEnvelope> {
    let (Some(a), Some(e)) = (aleatoric, epistemic) else {
        return None;
    };
    let delta_lo = (nominal - e.low).max(0.0);
    let delta_hi = (e.up - nominal).max(0.0);
    Some(PBoxEnvelope {
        bracket_95: [a.ci_95[0] - delta_lo, a.ci_95[1] + delta_hi],
        bracket_68: [a.ci_68[0] - delta_lo, a.ci_68[1] + delta_hi],
    })
}

/// Evaluates finite-difference linearized sensitivities + parallelized deterministic Monte Carlo statistics
/// for a batch of domain observables without sharing or conflating epistemic brackets.
fn propagate_domain<E>(
    vars: &[ResolvedVar],
    defs: &[EvalScalarDef],
    nom_vals: &[f64],
    seed: u64,
    n_samples: usize,
    eval_fn: E,
) -> Result<(Vec<Option<AleatoricEstimate>>, usize), Error>
where
    E: Fn(&[f64]) -> Result<Vec<f64>, Error> + Sync,
{
    let n_q = defs.len();
    if vars.is_empty() || n_q == 0 {
        return Ok((vec![None; n_q], 0));
    }
    let n_v = vars.len();

    // 1. Linearized sensitivities (central FD where inside bounds, forward/backward at hard bounds).
    let mut grad = vec![vec![0.0_f64; n_v]; n_q];
    for (vi, var) in vars.iter().enumerate() {
        let mut p_plus = vars.iter().map(|v| v.spec.nominal).collect::<Vec<_>>();
        let mut p_minus = p_plus.clone();
        let x0 = var.spec.nominal;
        let h_req = var.step_h.min(0.25 * (var.spec.max_bound - var.spec.min_bound)).max(1.0e-18);
        let can_minus = x0 - h_req >= var.spec.min_bound;
        let can_plus = x0 + h_req <= var.spec.max_bound;
        let (denom, y_lo, y_hi) = if can_plus && can_minus {
            p_plus[vi] = x0 + h_req;
            p_minus[vi] = x0 - h_req;
            (2.0 * h_req, eval_fn(&p_minus)?, eval_fn(&p_plus)?)
        } else if can_plus {
            p_plus[vi] = x0 + h_req;
            (h_req, nom_vals.to_vec(), eval_fn(&p_plus)?)
        } else if can_minus {
            p_minus[vi] = x0 - h_req;
            (h_req, eval_fn(&p_minus)?, nom_vals.to_vec())
        } else {
            (1.0, nom_vals.to_vec(), nom_vals.to_vec())
        };
        for qi in 0..n_q {
            let (a, b) = (y_lo[qi], y_hi[qi]);
            if a.is_finite() && b.is_finite() && denom > 0.0 {
                grad[qi][vi] = (b - a) / denom;
            }
        }
    }

    // 2. Pre-generate deterministic SplitMix64 Box-Muller draws sequentially for 0-ULP thread invariance.
    let mut draws = Vec::with_capacity(n_samples * n_v);
    let mut rng = NormalRng::new(seed);
    let mut clamped_count = 0usize;
    for _ in 0..n_samples {
        for var in vars {
            let raw_val = var.spec.nominal + var.spec.sigma * rng.next_normal();
            let clamped = raw_val.clamp(var.spec.min_bound, var.spec.max_bound);
            if clamped != raw_val {
                clamped_count += 1;
            }
            draws.push(clamped);
        }
    }

    // 3. Evaluate Monte Carlo ensemble in fixed chunks on the shared pool.
    let mc_rows: Vec<Result<Vec<f64>, Error>> = if n_samples > 0 {
        let sample_indices: Vec<usize> = (0..n_samples).collect();
        const MC_CHUNK: usize = 16;
        map_chunks(&sample_indices, MC_CHUNK, |idxs| {
            let mut local = Vec::with_capacity(idxs.len());
            for &si in idxs {
                let x = &draws[si * n_v..(si + 1) * n_v];
                local.push(eval_fn(x));
            }
            local
        })
        .into_iter()
        .flatten()
        .collect()
    } else {
        Vec::new()
    };

    let mut mc_cols = vec![Vec::with_capacity(n_samples); n_q];
    for row_res in mc_rows {
        let row = row_res?;
        for qi in 0..n_q {
            if row[qi].is_finite() {
                mc_cols[qi].push(row[qi]);
            }
        }
    }

    // 4. Assemble linear + MC estimates cleanly.
    let mut estimates = Vec::with_capacity(n_q);
    for qi in 0..n_q {
        let nom = nom_vals[qi];
        let mut var_lin = 0.0;
        let mut sens = BTreeMap::new();
        for (vi, var) in vars.iter().enumerate() {
            let g = grad[qi][vi];
            sens.insert(var.spec.name.clone(), g);
            let contrib = (g * var.spec.sigma) * (g * var.spec.sigma);
            var_lin += contrib;
        }
        let lin_sigma = var_lin.sqrt();
        let mut var_contrib = BTreeMap::new();
        if var_lin > 0.0 {
            for (vi, var) in vars.iter().enumerate() {
                let c = (grad[qi][vi] * var.spec.sigma).powi(2) / var_lin;
                var_contrib.insert(var.spec.name.clone(), c);
            }
        }

        let (mc_mean, mc_sigma, ci_68, ci_95) = if !mc_cols[qi].is_empty() {
            let vals = &mut mc_cols[qi];
            let m = vals.len() as f64;
            let mean = vals.iter().sum::<f64>() / m;
            let var = if vals.len() > 1 {
                vals.iter().map(|&y| (y - mean) * (y - mean)).sum::<f64>() / (m - 1.0)
            } else {
                0.0
            };
            vals.sort_by(|a, b| a.total_cmp(b));
            let c68 = [quantile_sorted(vals, 0.158_655_253_931_457), quantile_sorted(vals, 0.841_344_746_068_543)];
            let c95 = [quantile_sorted(vals, 0.025), quantile_sorted(vals, 0.975)];
            (mean, var.sqrt(), c68, c95)
        } else {
            (
                nom,
                lin_sigma,
                [nom - lin_sigma, nom + lin_sigma],
                [nom - 1.959_963_984_540_054 * lin_sigma, nom + 1.959_963_984_540_054 * lin_sigma],
            )
        };

        let denom = nom.abs().max(1e-300);
        estimates.push(Some(AleatoricEstimate {
            linear_sigma: lin_sigma,
            rel_linear_sigma: if nom.is_finite() && nom.abs() > 0.0 { lin_sigma / denom } else { 0.0 },
            mc_mean,
            mc_sigma,
            rel_mc_sigma: if nom.is_finite() && nom.abs() > 0.0 { mc_sigma / denom } else { 0.0 },
            ci_68,
            ci_95,
            variance_contributions: var_contrib,
            sensitivities: sens,
        }));
    }

    Ok((estimates, clamped_count))
}

/// Applies capacitance delta `delta_c_f` to a nonlinear `CircuitRom` port cleanly across `c_port_f` and `m_r`.
fn shift_rom_port_capacitance(rom: &mut CircuitRom, orig_port_idx: usize, nl_idx: usize, delta_c_f: f64) {
    if delta_c_f == 0.0 {
        return;
    }
    let r = rom.k_r.len();
    for i in 0..r {
        let bi = rom.b_r[i].get(orig_port_idx).copied().unwrap_or(0.0);
        if bi == 0.0 {
            continue;
        }
        for j in 0..r {
            let bj = rom.b_r[j].get(orig_port_idx).copied().unwrap_or(0.0);
            rom.m_r[i][j] += (delta_c_f / EPS0) * bi * bj;
        }
    }
    if rom.c_port_f.len() == rom.ports.len() && rom.c_port_f.iter().all(|row| row.len() == rom.ports.len()) {
        rom.c_port_f[orig_port_idx][orig_port_idx] += delta_c_f;
    } else if rom.c_port_f.len() == nl_idx.max(1) && rom.c_port_f.get(nl_idx).is_some_and(|row| row.len() > nl_idx) {
        rom.c_port_f[nl_idx][nl_idx] += delta_c_f;
    }
}

/// Extracts a clean reduced damping matrix (`c_internal + sum(alpha * b b^T)`) so perturbing `z0_ohm` on `rom_clone`
/// does not mismatch `rom.ports` vs `rom.c_r`.
fn decompose_rom_internal_damping(rom: &CircuitRom) -> Vec<Vec<f64>> {
    let (c_internal, _, _) = crate::epr::assemble_rom_damping(rom, &rom.ports, None);
    c_internal
}

struct RomEvalCtx<'a> {
    base_rom: &'a CircuitRom,
    c_internal_clean: &'a [Vec<f64>],
    orig_template: &'a SolveResult,
    base_params: &'a Value,
    use_hamiltonian: bool,
    vars_meta: &'a [RomVarBinding],
    tan_delta_bulk_extra: Option<(f64, usize)>,
}

struct UqRunState<'a> {
    raw_unc: &'a BTreeMap<String, RawUncSpec>,
    seed: u64,
    n_samples: usize,
    all_var_specs: &'a mut Vec<UqParameterSpec>,
    out_quantities: &'a mut Vec<UqQuantity>,
    domains_evaluated: &'a mut Vec<String>,
    diagnostics: &'a mut Vec<String>,
}

impl UqRunState<'_> {
    fn push_with_eb(
        &mut self,
        defs: &[EvalScalarDef],
        nom_vals: &[f64],
        aleatoric_vec: &[Option<AleatoricEstimate>],
        eb_result: &SolveResult,
        explicit_ebs: Option<&[Option<EpistemicInterval>]>,
    ) {
        let eb_rel = eb_result.error_budget.as_ref().and_then(|eb| eb.spatial_discretization.as_ref());
        for (qi, def) in defs.iter().enumerate() {
            let nom = nom_vals[qi];
            let epistemic = explicit_ebs.and_then(|s| s.get(qi).and_then(|o| o.clone())).or_else(|| {
                eb_rel.filter(|s| s.rel_bound > 0.0).map(|se| {
                    EpistemicInterval::from_rel_bound(
                        nom,
                        se.rel_bound,
                        def.non_negative,
                        format!("error_budget:{}", se.method),
                        se.is_epistemic_bound,
                    )
                })
            });
            let aleatoric = aleatoric_vec.get(qi).and_then(|o| o.clone());
            let outer_envelope = outer_pbox(nom, aleatoric.as_ref(), epistemic.as_ref());
            self.out_quantities.push(UqQuantity {
                name: def.name.clone(),
                unit: def.unit.clone(),
                nominal: nom,
                aleatoric,
                epistemic,
                outer_envelope,
            });
        }
    }
}

fn eval_rom_with_vars(ctx: &RomEvalCtx<'_>, vals: &[f64]) -> Result<RomEvalValues, Error> {
    let (base_rom, c_internal_clean, orig_template, base_params, use_hamiltonian, vars_meta, tan_delta_bulk_extra) = (
        ctx.base_rom,
        ctx.c_internal_clean,
        ctx.orig_template,
        ctx.base_params,
        ctx.use_hamiltonian,
        ctx.vars_meta,
        ctx.tan_delta_bulk_extra,
    );
    let mut rom = base_rom.clone();
    let r = rom.k_r.len();
    let nl_indices: Vec<usize> =
        rom.ports.iter().enumerate().filter(|(_, p)| is_nonlinear_port(p)).map(|(i, _)| i).collect();
    let scale = jj_flux_scale();
    let mut eval_params = base_params.clone();
    let mut extra_tan_delta_bulk = 0.0_f64;

    if let Some((nom_td, idx)) = tan_delta_bulk_extra {
        let _ = nom_td;
        extra_tan_delta_bulk = vals[idx] - nom_td;
    }

    for (vi, b) in vars_meta.iter().enumerate() {
        let val = vals[vi];
        match b {
            RomVarBinding::LjPort(p_idx) => {
                let p = &mut rom.ports[*p_idx];
                p.l_j_h = Some(val);
                p.e_j_hz = Some(scale / val);
                clear_port_override(&mut eval_params, &p.name, "lj_override");
                clear_port_override(&mut eval_params, &p.name, "ej_hz_override");
            }
            RomVarBinding::EjPort(p_idx) => {
                let p = &mut rom.ports[*p_idx];
                p.e_j_hz = Some(val);
                p.l_j_h = Some(scale / val);
                clear_port_override(&mut eval_params, &p.name, "lj_override");
                clear_port_override(&mut eval_params, &p.name, "ej_hz_override");
            }
            RomVarBinding::ElPort(p_idx) => {
                rom.ports[*p_idx].e_l_hz = Some(val);
            }
            RomVarBinding::CjPort(p_idx) => {
                rom.ports[*p_idx].c_j_f = val;
            }
            RomVarBinding::SquidAsymmetry { p_idx, ej_sum_hz, flux_phi0 } => {
                let p = &mut rom.ports[*p_idx];
                p.squid_asymmetry = Some(val);
                let ej_eff = squid_effective_ej(*ej_sum_hz, *flux_phi0, val);
                p.e_j_hz = (ej_eff > 0.0).then_some(ej_eff);
                p.l_j_h = (ej_eff > 0.0).then(|| scale / ej_eff);
                clear_port_override(&mut eval_params, &p.name, "squid_asymmetry");
                clear_port_override(&mut eval_params, &p.name, "squid_flux_phi0");
                clear_port_override(&mut eval_params, &p.name, "ej_hz_override");
            }
            RomVarBinding::SquidFlux { p_idx, ej_sum_hz, asymmetry } => {
                let p = &mut rom.ports[*p_idx];
                let ej_eff = squid_effective_ej(*ej_sum_hz, val, *asymmetry);
                p.e_j_hz = (ej_eff > 0.0).then_some(ej_eff);
                p.l_j_h = (ej_eff > 0.0).then(|| scale / ej_eff);
                p.phi_ext_rad = 2.0 * PI * val;
                clear_port_override(&mut eval_params, &p.name, "squid_flux_phi0");
                clear_port_override(&mut eval_params, &p.name, "flux_phi0");
                clear_port_override(&mut eval_params, &p.name, "phi_ext_rad");
                clear_port_override(&mut eval_params, &p.name, "ej_hz_override");
            }
            RomVarBinding::LoopFluxPhi0(p_idx) => {
                let p = &mut rom.ports[*p_idx];
                p.phi_ext_rad = 2.0 * PI * val;
                clear_param_key(&mut eval_params, "flux_phi0");
                clear_param_key(&mut eval_params, "phi_ext_rad");
            }
            RomVarBinding::GateCharge(p_idx) => {
                rom.ports[*p_idx].n_g = val;
                clear_param_key(&mut eval_params, "n_g");
            }
            RomVarBinding::PortZ0(p_idx) => {
                rom.ports[*p_idx].z0_ohm = Some(val);
            }
            RomVarBinding::TanDeltaBulk => {}
        }
    }

    // Reconstruct c_r cleanly whenever any z0_ohm port was varied
    if vars_meta.iter().any(|b| matches!(b, RomVarBinding::PortZ0(_))) && r > 0 {
        let mut c_new = c_internal_clean.to_vec();
        for (p_idx, port) in rom.ports.iter().enumerate() {
            if port.kind == "lumped_port" {
                if let Some(z0) = port.z0_ohm.filter(|&z| z.is_finite() && z > 0.0) {
                    let alpha = MU0 * C0 / z0;
                    for i in 0..r {
                        let bi = rom.b_r[i].get(p_idx).copied().unwrap_or(0.0);
                        if bi == 0.0 {
                            continue;
                        }
                        for j in 0..r {
                            let bj = rom.b_r[j].get(p_idx).copied().unwrap_or(0.0);
                            c_new[i][j] += alpha * bi * bj;
                        }
                    }
                }
            }
        }
        rom.c_r = c_new;
    }

    // Clear modes so overrides reproject cleanly via r x r ROM eigenproblem
    if !vars_meta.is_empty() {
        rom.modes.clear();
    }

    let synth = SolveResult {
        physics: orig_template.physics.clone(),
        study: orig_template.study.clone(),
        order: orig_template.order,
        dim: orig_template.dim,
        scalars: vec![rom.scalar()],
        ..Default::default()
    };
    let epr_rep = epr_with_params(&synth, &eval_params)?;
    let ham_rep = if use_hamiltonian && !nl_indices.is_empty() {
        Some(hamiltonian_with_params(&synth, &eval_params)?)
    } else {
        None
    };

    Ok(extract_rom_eval_values(&epr_rep, ham_rep.as_ref(), extra_tan_delta_bulk))
}

fn clear_port_override(params: &mut Value, port_name: &str, key: &str) {
    if let Some(obj) = params.as_object_mut() {
        if let Some(sub) = obj.get_mut(key).and_then(Value::as_object_mut) {
            sub.remove(port_name);
            if sub.is_empty() {
                obj.remove(key);
            }
        } else if obj.get(key).is_some_and(|v| !v.is_object()) {
            obj.remove(key);
        }
    }
}

fn clear_param_key(params: &mut Value, key: &str) {
    if let Some(obj) = params.as_object_mut() {
        obj.remove(key);
    }
}

#[derive(Clone, Debug)]
enum RomVarBinding {
    LjPort(usize),
    EjPort(usize),
    ElPort(usize),
    CjPort(usize),
    SquidAsymmetry { p_idx: usize, ej_sum_hz: f64, flux_phi0: f64 },
    SquidFlux { p_idx: usize, ej_sum_hz: f64, asymmetry: f64 },
    LoopFluxPhi0(usize),
    GateCharge(usize),
    PortZ0(usize),
    TanDeltaBulk,
}

#[derive(Clone, Debug)]
struct RomEvalValues {
    defs: Vec<EvalScalarDef>,
    vals: Vec<f64>,
}

fn extract_rom_eval_values(epr: &EprReport, ham: Option<&HamiltonianReport>, extra_tan_delta: f64) -> RomEvalValues {
    let mut defs = Vec::new();
    let mut vals = Vec::new();
    if let Some(h) = ham {
        let primary_q_idx = if epr.ports.len() == 1 && h.f01_hz.len() > 1 {
            (0..h.f01_hz.len())
                .max_by(|&a, &b| {
                    let pa = epr.p_matrix.get(a).and_then(|r| r.first()).copied().unwrap_or(0.0);
                    let pb = epr.p_matrix.get(b).and_then(|r| r.first()).copied().unwrap_or(0.0);
                    pa.total_cmp(&pb)
                })
                .unwrap_or(0)
        } else {
            0
        };
        for (j, name) in h.subsystems.iter().take(h.f01_hz.len()).enumerate() {
            let is_single_primary = h.f01_hz.len() == 1 || (epr.ports.len() == 1 && j == primary_q_idx);
            let suffix = if is_single_primary { String::new() } else { format!("_{name}") };
            defs.push(EvalScalarDef { name: format!("f01_Hz{suffix}"), unit: "Hz".into(), non_negative: true });
            vals.push(h.f01_hz[j]);
            defs.push(EvalScalarDef { name: format!("alpha_Hz{suffix}"), unit: "Hz".into(), non_negative: false });
            vals.push(h.alpha_hz[j]);
            if j < h.f12_hz.len() {
                defs.push(EvalScalarDef { name: format!("f12_Hz{suffix}"), unit: "Hz".into(), non_negative: true });
                vals.push(h.f12_hz[j]);
            }
        }
        for (m, &fc) in h.f_cav_hz.iter().enumerate() {
            let cs = if h.f_cav_hz.len() > 1 { format!("_c{m}") } else { String::new() };
            defs.push(EvalScalarDef { name: format!("f_cav_Hz{cs}"), unit: "Hz".into(), non_negative: true });
            vals.push(fc);
            let q_rom_inv = if m < h.kappa_cav_hz.len() && fc > 0.0 && h.kappa_cav_hz[m] > 0.0 {
                h.kappa_cav_hz[m] / fc
            } else {
                0.0
            };
            let loss_tot = (q_rom_inv + extra_tan_delta).max(0.0);
            let kappa_tot =
                if extra_tan_delta != 0.0 { loss_tot * fc } else { h.kappa_cav_hz.get(m).copied().unwrap_or(0.0) };
            defs.push(EvalScalarDef { name: format!("kappa_cav_Hz{cs}"), unit: "Hz".into(), non_negative: true });
            vals.push(kappa_tot);
            if loss_tot > 0.0 {
                defs.push(EvalScalarDef { name: format!("Q_cav{cs}"), unit: "1".into(), non_negative: true });
                vals.push(1.0 / loss_tot);
            }
            for (j, pname) in epr.ports.iter().take(h.g_matrix_hz.len()).enumerate() {
                if let Some(&g) = h.g_matrix_hz[j].get(m) {
                    let gs = if h.g_matrix_hz.len() > 1 || h.f_cav_hz.len() > 1 {
                        format!("_{pname}_c{m}")
                    } else {
                        String::new()
                    };
                    defs.push(EvalScalarDef { name: format!("g_Hz{gs}"), unit: "Hz".into(), non_negative: false });
                    vals.push(g);
                }
            }
        }
        let n_chi = h.chi_matrix_hz.len();
        for a in 0..n_chi {
            for b in (a + 1)..h.chi_matrix_hz[a].len() {
                let name =
                    if n_chi == 2 && a == 0 && b == 1 { "chi_Hz".to_string() } else { format!("chi_{a}_{b}_Hz") };
                defs.push(EvalScalarDef { name, unit: "Hz".into(), non_negative: false });
                vals.push(h.chi_matrix_hz[a][b]);
            }
        }
        for (j, &t1) in h.t1_purcell_s.iter().enumerate() {
            if t1.is_finite() {
                let suffix = if h.t1_purcell_s.len() > 1 { format!("_q{j}") } else { String::new() };
                defs.push(EvalScalarDef {
                    name: format!("t1_purcell_s{suffix}"),
                    unit: "s".into(),
                    non_negative: true,
                });
                vals.push(t1);
            }
        }
    } else {
        // EPR-only fallback when Hamiltonian basis = "epr" or no nonlinear ports exist
        for (m, &fm) in epr.modes_hz.iter().enumerate() {
            let ms = if epr.modes_hz.len() > 1 { format!("_m{m}") } else { String::new() };
            defs.push(EvalScalarDef { name: format!("f_dressed_pert_Hz{ms}"), unit: "Hz".into(), non_negative: true });
            vals.push(epr.f_dressed_pert_hz.get(m).copied().unwrap_or(fm));
            if let Some(&a) = epr.alpha_hz.get(m) {
                defs.push(EvalScalarDef { name: format!("alpha_Hz{ms}"), unit: "Hz".into(), non_negative: false });
                vals.push(a);
            }
            if let Some(Some(q0)) = epr.q_factors.get(m) {
                let q_inv = (1.0 / *q0 + extra_tan_delta).max(1e-30);
                defs.push(EvalScalarDef { name: format!("Q_m{ms}"), unit: "1".into(), non_negative: true });
                vals.push(1.0 / q_inv);
            }
        }
        for a in 0..epr.chi_matrix_hz.len() {
            for b in (a + 1)..epr.chi_matrix_hz[a].len() {
                let name = if epr.chi_matrix_hz.len() == 2 && a == 0 && b == 1 {
                    "chi_Hz".to_string()
                } else {
                    format!("chi_{a}_{b}_Hz")
                };
                defs.push(EvalScalarDef { name, unit: "Hz".into(), non_negative: false });
                vals.push(epr.chi_matrix_hz[a][b]);
            }
        }
    }
    RomEvalValues { defs, vals }
}

fn eval_rom_bounds_from_companion(
    base_rom: &CircuitRom,
    orig_template: &SolveResult,
    base_params: &Value,
    use_hamiltonian: bool,
    companion: &SolveResult,
    cap_params: &Map<String, Value>,
) -> Result<Option<(RomEvalValues, RomEvalValues, String)>, Error> {
    if !capacitance_matrix_applies(companion) {
        return Ok(None);
    }
    let cap_nom = run_capacitance_matrix_from_thickness(&capacitance_matrix(companion)?, None, cap_params)?;
    let Some(qubits_arr) = cap_nom.get("qubits").and_then(Value::as_array) else {
        return Ok(None);
    };
    let nl_indices: Vec<usize> =
        base_rom.ports.iter().enumerate().filter(|(_, p)| is_nonlinear_port(p)).map(|(i, _)| i).collect();
    if nl_indices.is_empty() {
        return Ok(None);
    }

    let mut rom_lo_cap = base_rom.clone();
    let mut rom_hi_cap = base_rom.clone();
    let mut any_bounds = false;

    for (nl_idx, &orig_idx) in nl_indices.iter().enumerate() {
        let port_name = &base_rom.ports[orig_idx].name;
        let q_match = qubits_arr
            .iter()
            .find(|q| q.get("name").and_then(Value::as_str) == Some(port_name.as_str()))
            .or_else(|| (qubits_arr.len() == 1 && nl_indices.len() == 1).then(|| &qubits_arr[0]));
        let Some(q) = q_match else { continue };
        let (Some(c_nom), Some(c_lo), Some(c_hi)) = (
            q.get("C_Sigma_F").and_then(Value::as_f64),
            q.get("C_Sigma_low_F").and_then(Value::as_f64),
            q.get("C_Sigma_up_F").and_then(Value::as_f64),
        ) else {
            continue;
        };
        any_bounds = true;
        shift_rom_port_capacitance(&mut rom_lo_cap, orig_idx, nl_idx, c_lo - c_nom);
        shift_rom_port_capacitance(&mut rom_hi_cap, orig_idx, nl_idx, c_hi - c_nom);
    }
    if !any_bounds {
        return Ok(None);
    }

    rom_lo_cap.modes.clear();
    rom_hi_cap.modes.clear();

    let eval_one = |rom: CircuitRom| -> Result<RomEvalValues, Error> {
        let synth = SolveResult {
            physics: orig_template.physics.clone(),
            study: orig_template.study.clone(),
            order: orig_template.order,
            dim: orig_template.dim,
            scalars: vec![rom.scalar()],
            ..Default::default()
        };
        let epr_rep = epr_with_params(&synth, base_params)?;
        let ham_rep = if use_hamiltonian { Some(hamiltonian_with_params(&synth, base_params)?) } else { None };
        Ok(extract_rom_eval_values(&epr_rep, ham_rep.as_ref(), 0.0))
    };

    Ok(Some((eval_one(rom_lo_cap)?, eval_one(rom_hi_cap)?, "prager_synge_capacitance".to_string())))
}

fn run_rom_uq(
    result: &SolveResult,
    companion_result: Option<&SolveResult>,
    params: &Value,
    state: &mut UqRunState<'_>,
) -> Result<(), Error> {
    let raw_unc = state.raw_unc;
    let base_rom = CircuitRom::of(result).ok_or_else(|| Error::invalid("SolveResult has no CircuitRom"))?;
    let basis = params.get("basis").and_then(Value::as_str).unwrap_or("normal_fock");
    if !matches!(basis, "normal_fock" | "craig_bampton" | "epr") {
        return Err(Error::invalid(format!(
            "unknown basis '{basis}' (expected 'normal_fock', 'craig_bampton', or 'epr')"
        )));
    }
    let use_hamiltonian = basis != "epr" && hamiltonian_applies(result);

    // Prepare ROM nominal with explicit parameter overrides (E_J/L_J, SQUID, flux, n_g) applied
    let (eff_ports, _) = crate::epr::resolve_ports_with_overrides(&base_rom.ports, params)?;
    let mut rom_nom = base_rom.clone();
    rom_nom.ports = eff_ports;
    let c_internal_clean = decompose_rom_internal_damping(&base_rom);
    let nl_indices: Vec<usize> =
        rom_nom.ports.iter().enumerate().filter(|(_, p)| is_nonlinear_port(p)).map(|(i, _)| i).collect();

    let mut vars = Vec::new();
    let mut bindings = Vec::new();
    let mut tan_delta_bulk_extra = None;

    for (p_idx, port) in rom_nom.ports.iter().enumerate() {
        let qtag = if nl_indices.len() > 1 { format!(":{}", port.name) } else { String::new() };
        if is_nonlinear_port(port) {
            let (l_nom, e_nom) = resolve_jj_l_and_ej(port.l_j_h.unwrap_or(0.0), port.e_j_hz.unwrap_or(0.0));
            let ej_key_q = format!("E_J_Hz:{}", port.name);
            let lj_key_q = format!("L_J_H:{}", port.name);
            let has_ej =
                raw_unc.contains_key(&ej_key_q) || raw_unc.contains_key("E_J_Hz") || raw_unc.contains_key("ej_hz");
            let has_lj =
                raw_unc.contains_key(&lj_key_q) || raw_unc.contains_key("L_J_H") || raw_unc.contains_key("lj_h");
            if has_ej && has_lj {
                return Err(Error::invalid(format!(
                    "cannot specify aleatoric uncertainties on both E_J_Hz and L_J_H for port '{}'",
                    port.name
                )));
            }
            if has_ej && e_nom > 0.0 {
                let k = if raw_unc.contains_key(&ej_key_q) {
                    ej_key_q.as_str()
                } else if raw_unc.contains_key("E_J_Hz") {
                    "E_J_Hz"
                } else {
                    "ej_hz"
                };
                if let Some(rv) = resolve_var(raw_unc, k, format!("E_J_Hz{qtag}"), "rom", e_nom, 1e-12, f64::INFINITY)?
                {
                    vars.push(rv);
                    bindings.push(RomVarBinding::EjPort(p_idx));
                }
            } else if has_lj && l_nom > 0.0 {
                let k = if raw_unc.contains_key(&lj_key_q) {
                    lj_key_q.as_str()
                } else if raw_unc.contains_key("L_J_H") {
                    "L_J_H"
                } else {
                    "lj_h"
                };
                if let Some(rv) = resolve_var(raw_unc, k, format!("L_J_H{qtag}"), "rom", l_nom, 1e-24, f64::INFINITY)? {
                    vars.push(rv);
                    bindings.push(RomVarBinding::LjPort(p_idx));
                }
            }

            if let Some(el_nom) = port.e_l_hz.filter(|&v| v > 0.0) {
                let el_key_q = format!("E_L_Hz:{}", port.name);
                let k = if raw_unc.contains_key(&el_key_q) {
                    el_key_q.as_str()
                } else if raw_unc.contains_key("E_L_Hz") {
                    "E_L_Hz"
                } else {
                    "el_hz"
                };
                if let Some(rv) = resolve_var(raw_unc, k, format!("E_L_Hz{qtag}"), "rom", el_nom, 0.0, f64::INFINITY)? {
                    vars.push(rv);
                    bindings.push(RomVarBinding::ElPort(p_idx));
                }
            }

            let cj_key_q = format!("C_J_F:{}", port.name);
            let ck = if raw_unc.contains_key(&cj_key_q) {
                cj_key_q.as_str()
            } else if raw_unc.contains_key("C_J_F") {
                "C_J_F"
            } else {
                "c_j_f"
            };
            if let Some(rv) = resolve_var(raw_unc, ck, format!("C_J_F{qtag}"), "rom", port.c_j_f, 0.0, f64::INFINITY)? {
                vars.push(rv);
                bindings.push(RomVarBinding::CjPort(p_idx));
            }

            if let Some(asym_nom) = port.squid_asymmetry {
                let base_ej_sum = resolve_jj_l_and_ej(
                    base_rom.ports[p_idx].l_j_h.unwrap_or(0.0),
                    base_rom.ports[p_idx].e_j_hz.unwrap_or(0.0),
                )
                .1;
                let sflux_nom = params
                    .get("squid_flux_phi0")
                    .and_then(|v| v.as_f64().or_else(|| v.get(&port.name).and_then(Value::as_f64)))
                    .or_else(|| params.get("flux_phi0").and_then(Value::as_f64))
                    .unwrap_or(port.phi_ext_rad / (2.0 * PI));

                if let Some(rv) = resolve_var(
                    raw_unc,
                    "squid_asymmetry",
                    format!("squid_asymmetry{qtag}"),
                    "rom",
                    asym_nom,
                    0.0,
                    1.0,
                )? {
                    vars.push(rv);
                    bindings.push(RomVarBinding::SquidAsymmetry {
                        p_idx,
                        ej_sum_hz: base_ej_sum,
                        flux_phi0: sflux_nom,
                    });
                }
                if let Some(rv) = resolve_var(
                    raw_unc,
                    "squid_flux_phi0",
                    format!("squid_flux_phi0{qtag}"),
                    "rom",
                    sflux_nom,
                    f64::NEG_INFINITY,
                    f64::INFINITY,
                )? {
                    vars.push(rv);
                    bindings.push(RomVarBinding::SquidFlux { p_idx, ej_sum_hz: base_ej_sum, asymmetry: asym_nom });
                }
            }

            let phi_nom_phi0 = params
                .get("flux_phi0")
                .and_then(Value::as_f64)
                .or_else(|| params.get("phi_ext_rad").and_then(Value::as_f64).map(|r| r / (2.0 * PI)))
                .unwrap_or(port.phi_ext_rad / (2.0 * PI));
            if let Some(rv) = resolve_var(
                raw_unc,
                "flux_phi0",
                format!("flux_phi0{qtag}"),
                "rom",
                phi_nom_phi0,
                f64::NEG_INFINITY,
                f64::INFINITY,
            )? {
                vars.push(rv);
                bindings.push(RomVarBinding::LoopFluxPhi0(p_idx));
            }
            let ng_nom = params.get("n_g").and_then(Value::as_f64).unwrap_or(port.n_g);
            if let Some(rv) =
                resolve_var(raw_unc, "n_g", format!("n_g{qtag}"), "rom", ng_nom, f64::NEG_INFINITY, f64::INFINITY)?
            {
                vars.push(rv);
                bindings.push(RomVarBinding::GateCharge(p_idx));
            }
        } else if port.kind == "lumped_port" {
            if let Some(z0_nom) = port.z0_ohm.filter(|&z| z > 0.0) {
                let zq = format!("Z0_ohm:{}", port.name);
                let zk = if raw_unc.contains_key(&zq) {
                    zq.as_str()
                } else if raw_unc.contains_key("Z0_ohm") {
                    "Z0_ohm"
                } else {
                    "z0_ohm"
                };
                if let Some(rv) =
                    resolve_var(raw_unc, zk, format!("Z0_ohm:{}", port.name), "rom", z0_nom, 1e-6, f64::INFINITY)?
                {
                    vars.push(rv);
                    bindings.push(RomVarBinding::PortZ0(p_idx));
                }
            }
        }
    }

    let td_bulk_nom = params.get("tan_delta_bulk").and_then(Value::as_f64).unwrap_or(0.0);
    if let Some(rv) = resolve_var(raw_unc, "tan_delta_bulk", "tan_delta_bulk", "rom", td_bulk_nom, 0.0, f64::INFINITY)?
    {
        let idx = vars.len();
        tan_delta_bulk_extra = Some((rv.spec.nominal, idx));
        vars.push(rv);
        bindings.push(RomVarBinding::TanDeltaBulk);
    }

    let nom_vals_vec = vars.iter().map(|v| v.spec.nominal).collect::<Vec<_>>();
    let rom_ctx = RomEvalCtx {
        base_rom: &rom_nom,
        c_internal_clean: &c_internal_clean,
        orig_template: result,
        base_params: params,
        use_hamiltonian,
        vars_meta: &bindings,
        tan_delta_bulk_extra,
    };
    let nom_eval = eval_rom_with_vars(&rom_ctx, &nom_vals_vec)?;

    let (aleatoric_vec, clamped_cnt) =
        propagate_domain(&vars, &nom_eval.defs, &nom_eval.vals, state.seed ^ 0xA1E4_0001, state.n_samples, |x| {
            eval_rom_with_vars(&rom_ctx, x).map(|r| r.vals)
        })?;
    if clamped_cnt > 0 {
        state.diagnostics.push(format!("rom: clamped {clamped_cnt} normal MC draw(s) at physical domain bounds"));
    }

    let cap_map = params.as_object().cloned().unwrap_or_default();
    let companion_bounds = companion_result
        .map(|comp| eval_rom_bounds_from_companion(&rom_nom, result, params, use_hamiltonian, comp, &cap_map))
        .transpose()?
        .flatten();
    let companion_ep: Option<Vec<Option<EpistemicInterval>>> = companion_bounds.as_ref().map(|(lo_ev, hi_ev, src)| {
        nom_eval
            .defs
            .iter()
            .enumerate()
            .map(|(qi, _)| {
                lo_ev.vals.get(qi).zip(hi_ev.vals.get(qi)).map(|(&v_lo, &v_hi)| {
                    EpistemicInterval::from_bounds(v_lo, v_hi, nom_eval.vals[qi], src.clone(), true)
                })
            })
            .collect()
    });
    state.push_with_eb(&nom_eval.defs, &nom_eval.vals, &aleatoric_vec, result, companion_ep.as_deref());

    for rv in vars {
        state.all_var_specs.push(rv.spec);
    }
    state.domains_evaluated.push(if use_hamiltonian { format!("rom:{basis}") } else { "rom:epr".to_string() });
    Ok(())
}

#[derive(Clone, Debug)]
enum CapVarBinding {
    MetalThickness,
    C0,
    C1,
    QubitEj(usize),
    ResFr(usize),
    ResZ0(usize),
    ResCr(usize),
}

#[derive(Clone, Debug)]
struct CapEvalValues {
    defs: Vec<EvalScalarDef>,
    vals: Vec<f64>,
    epistemic_bounds: Vec<Option<EpistemicInterval>>,
}

fn patch_cap_params(
    base_params: &Map<String, Value>,
    bindings: &[CapVarBinding],
    vals: &[f64],
) -> (Option<f64>, Option<f64>, Option<f64>, Map<String, Value>) {
    let mut t_metal = base_params.get("metal_thickness_m").and_then(Value::as_f64);
    let mut c0 = base_params.get("c0").and_then(Value::as_f64);
    let mut c1 = base_params.get("c1").and_then(Value::as_f64);
    let mut patched = base_params.clone();

    for (vi, b) in bindings.iter().enumerate() {
        let v = vals[vi];
        match b {
            CapVarBinding::MetalThickness => t_metal = Some(v),
            CapVarBinding::C0 => c0 = Some(v),
            CapVarBinding::C1 => c1 = Some(v),
            CapVarBinding::QubitEj(qi) => {
                if let Some(arr) = patched.get_mut("qubits").and_then(Value::as_array_mut) {
                    if let Some(obj) = arr.get_mut(*qi).and_then(Value::as_object_mut) {
                        obj.insert("E_J_Hz".into(), json!(v));
                    }
                }
            }
            CapVarBinding::ResFr(ri) => {
                if let Some(arr) = patched.get_mut("resonators").and_then(Value::as_array_mut) {
                    if let Some(obj) = arr.get_mut(*ri).and_then(Value::as_object_mut) {
                        obj.insert("f_r_Hz".into(), json!(v));
                    }
                }
            }
            CapVarBinding::ResZ0(ri) => {
                if let Some(arr) = patched.get_mut("resonators").and_then(Value::as_array_mut) {
                    if let Some(obj) = arr.get_mut(*ri).and_then(Value::as_object_mut) {
                        obj.insert("Z0_ohm".into(), json!(v));
                    }
                }
            }
            CapVarBinding::ResCr(ri) => {
                if let Some(arr) = patched.get_mut("resonators").and_then(Value::as_array_mut) {
                    if let Some(obj) = arr.get_mut(*ri).and_then(Value::as_object_mut) {
                        obj.insert("C_r_F".into(), json!(v));
                    }
                }
            }
        }
    }
    (t_metal, c0, c1, patched)
}

fn eval_cap_with_vars(
    cap: &CapacitanceMatrixData,
    gram_opt: Option<&PreparedThicknessGram>,
    base_params: &Map<String, Value>,
    bindings: &[CapVarBinding],
    vals: &[f64],
    include_epistemic_bounds: bool,
) -> Result<CapEvalValues, Error> {
    let (t_metal, c0, c1, patched_map) = patch_cap_params(base_params, bindings, vals);
    let tc = match (gram_opt, t_metal) {
        (Some(gram), Some(t)) => Some(apply_3d_thickness_gram(cap, gram, t, c0, c1)?),
        _ => None,
    };
    let out = run_capacitance_matrix_from_thickness(cap, tc, &patched_map)?;

    let mut defs = Vec::new();
    let mut qvals = Vec::new();
    let mut ebs = Vec::new();

    if let Some(qubits) = out.get("qubits").and_then(Value::as_array) {
        let multi_q = qubits.len() > 1;
        for q in qubits {
            let qname = q.get("name").and_then(Value::as_str).unwrap_or("q");
            let prefix = if multi_q { format!("{qname}.") } else { String::new() };
            let c_nom = q.get("C_Sigma_F").and_then(Value::as_f64).unwrap_or(0.0);
            let ec_nom = q.get("E_C_Hz").and_then(Value::as_f64).unwrap_or(0.0);
            let cs_bounds =
                q.get("C_Sigma_low_F").and_then(Value::as_f64).zip(q.get("C_Sigma_up_F").and_then(Value::as_f64));
            let ec_bounds = q.get("E_C_low_Hz").and_then(Value::as_f64).zip(q.get("E_C_up_Hz").and_then(Value::as_f64));

            defs.push(EvalScalarDef { name: format!("{prefix}C_Sigma_F"), unit: "F".into(), non_negative: true });
            qvals.push(c_nom);
            ebs.push(if include_epistemic_bounds {
                cs_bounds.map(|(lo, hi)| EpistemicInterval::from_bounds(lo, hi, c_nom, "prager_synge", true))
            } else {
                None
            });

            defs.push(EvalScalarDef { name: format!("{prefix}E_C_Hz"), unit: "Hz".into(), non_negative: true });
            qvals.push(ec_nom);
            ebs.push(if include_epistemic_bounds {
                ec_bounds.map(|(lo, hi)| EpistemicInterval::from_bounds(lo, hi, ec_nom, "prager_synge", true))
            } else {
                None
            });

            if let (Some(f01), Some(ej)) =
                (q.get("f01_Hz").and_then(Value::as_f64), q.get("E_J_Hz").and_then(Value::as_f64))
            {
                defs.push(EvalScalarDef { name: format!("{prefix}f01_Hz"), unit: "Hz".into(), non_negative: true });
                qvals.push(f01);
                ebs.push(if include_epistemic_bounds {
                    ec_bounds.map(|(ec_lo, ec_hi)| {
                        let f_a = (8.0 * ej * ec_lo).sqrt() - ec_lo;
                        let f_b = (8.0 * ej * ec_hi).sqrt() - ec_hi;
                        EpistemicInterval::from_bounds(f_a, f_b, f01, "prager_synge", true)
                    })
                } else {
                    None
                });

                let alpha = q.get("anharmonicity_Hz").and_then(Value::as_f64).unwrap_or(-ec_nom);
                defs.push(EvalScalarDef {
                    name: format!("{prefix}anharmonicity_Hz"),
                    unit: "Hz".into(),
                    non_negative: false,
                });
                qvals.push(alpha);
                ebs.push(if include_epistemic_bounds {
                    ec_bounds.map(|(ec_lo, ec_hi)| {
                        EpistemicInterval::from_bounds(-ec_hi, -ec_lo, alpha, "prager_synge", true)
                    })
                } else {
                    None
                });
            }
        }
    } else if let Some(mat) = out.get("maxwell_reduced_F").and_then(Value::as_array) {
        let terms = out.get("reduced_terminals").and_then(Value::as_array);
        let lo_mat = out.get("maxwell_reduced_lower_F").and_then(Value::as_array);
        let hi_mat = out.get("maxwell_reduced_upper_F").and_then(Value::as_array);
        for (i, row) in mat.iter().enumerate() {
            if let Some(c_ii) = row.get(i).and_then(Value::as_f64) {
                let tname = terms.and_then(|t| t.get(i)).and_then(Value::as_str).unwrap_or("t");
                defs.push(EvalScalarDef { name: format!("C_{tname}_{tname}_F"), unit: "F".into(), non_negative: true });
                qvals.push(c_ii);
                let eb = if include_epistemic_bounds {
                    lo_mat
                        .and_then(|m| m.get(i).and_then(|r| r.get(i)).and_then(Value::as_f64))
                        .zip(hi_mat.and_then(|m| m.get(i).and_then(|r| r.get(i)).and_then(Value::as_f64)))
                        .map(|(lo, hi)| EpistemicInterval::from_bounds(lo, hi, c_ii, "prager_synge", true))
                } else {
                    None
                };
                ebs.push(eb);
            }
        }
    }

    if let Some(couplings) = out.get("couplings").and_then(Value::as_array) {
        let multi_c = couplings.len() > 1;
        let qubits = out.get("qubits").and_then(Value::as_array);
        let resonators = out.get("resonators").and_then(Value::as_array);
        for c_obj in couplings {
            let qname = c_obj.get("qubit").and_then(Value::as_str).unwrap_or("q");
            let rname = c_obj.get("resonator").and_then(Value::as_str).unwrap_or("r");
            let prefix = if multi_c { format!("{qname}_{rname}.") } else { String::new() };
            if let Some(g) = c_obj.get("g_Hz").and_then(Value::as_f64) {
                defs.push(EvalScalarDef { name: format!("{prefix}g_Hz"), unit: "Hz".into(), non_negative: true });
                qvals.push(g);
                let eb = if include_epistemic_bounds {
                    compute_coupling_epistemic_g(c_obj, qubits, resonators, qname, rname, g)
                } else {
                    None
                };
                ebs.push(eb);
            }
            if let Some(chi) = c_obj.get("chi_Hz").and_then(Value::as_f64) {
                defs.push(EvalScalarDef { name: format!("{prefix}chi_Hz"), unit: "Hz".into(), non_negative: false });
                qvals.push(chi);
                let eb = if include_epistemic_bounds {
                    compute_coupling_epistemic_chi(c_obj, qubits, resonators, qname, rname, chi)
                } else {
                    None
                };
                ebs.push(eb);
            }
        }
    }

    Ok(CapEvalValues { defs, vals: qvals, epistemic_bounds: ebs })
}

fn compute_coupling_epistemic_g(
    c_obj: &Value,
    qubits: Option<&Vec<Value>>,
    resonators: Option<&Vec<Value>>,
    qname: &str,
    _rname: &str,
    g_nom: f64,
) -> Option<EpistemicInterval> {
    let q = qubits?.iter().find(|q| q.get("name").and_then(Value::as_str) == Some(qname))?;
    let (cs_lo, cs_hi) =
        q.get("C_Sigma_low_F").and_then(Value::as_f64).zip(q.get("C_Sigma_up_F").and_then(Value::as_f64))?;
    let (ec_lo, ec_hi) = q.get("E_C_low_Hz").and_then(Value::as_f64).zip(q.get("E_C_up_Hz").and_then(Value::as_f64))?;
    let ej = q.get("E_J_Hz").and_then(Value::as_f64)?;
    let c_eff = c_obj.get("C_eff_F").and_then(Value::as_f64)?;
    let v_rms = c_obj
        .get("V_rms_V")
        .and_then(Value::as_f64)
        .or_else(|| resonators?.first().and_then(|r| r.get("V_rms_V").and_then(Value::as_f64)))?;
    let eval_g = |cs: f64, ec: f64| -> f64 {
        let beta = c_eff / cs.max(1e-30);
        let nzpf = (ej / (32.0 * ec)).powf(0.25);
        2.0 * crate::epr::E_CHARGE * beta * nzpf * v_rms / crate::epr::PLANCK_H
    };
    let g_a = eval_g(cs_hi, ec_lo);
    let g_b = eval_g(cs_lo, ec_hi);
    Some(EpistemicInterval::from_bounds(g_a, g_b, g_nom, "prager_synge", true))
}

fn compute_coupling_epistemic_chi(
    c_obj: &Value,
    qubits: Option<&Vec<Value>>,
    resonators: Option<&Vec<Value>>,
    qname: &str,
    rname: &str,
    chi_nom: f64,
) -> Option<EpistemicInterval> {
    let q = qubits?.iter().find(|q| q.get("name").and_then(Value::as_str) == Some(qname))?;
    let r = resonators?.iter().find(|r| r.get("name").and_then(Value::as_str) == Some(rname))?;
    let (cs_lo, cs_hi) =
        q.get("C_Sigma_low_F").and_then(Value::as_f64).zip(q.get("C_Sigma_up_F").and_then(Value::as_f64))?;
    let (ec_lo, ec_hi) = q.get("E_C_low_Hz").and_then(Value::as_f64).zip(q.get("E_C_up_Hz").and_then(Value::as_f64))?;
    let ej = q.get("E_J_Hz").and_then(Value::as_f64)?;
    let c_eff = c_obj.get("C_eff_F").and_then(Value::as_f64)?;
    let v_rms = r.get("V_rms_V").and_then(Value::as_f64)?;
    let f_r = r.get("f_r_Hz").and_then(Value::as_f64)?;
    let eval_chi = |cs: f64, ec: f64| -> f64 {
        let beta = c_eff / cs.max(1e-30);
        let nzpf = (ej / (32.0 * ec)).powf(0.25);
        let g = 2.0 * crate::epr::E_CHARGE * beta * nzpf * v_rms / crate::epr::PLANCK_H;
        let f01 = (8.0 * ej * ec).sqrt() - ec;
        let delta = f01 - f_r;
        -(g * g) / delta * (ec / (delta - ec))
    };
    let c_a = eval_chi(cs_hi, ec_lo);
    let c_b = eval_chi(cs_lo, ec_hi);
    if !(c_a.is_finite() && c_b.is_finite()) {
        return None;
    }
    Some(EpistemicInterval::from_bounds(c_a, c_b, chi_nom, "prager_synge", true))
}

fn run_capacitance_uq(
    mesh: Option<&Mesh>,
    result: &SolveResult,
    params: &Value,
    state: &mut UqRunState<'_>,
) -> Result<(), Error> {
    let raw_unc = state.raw_unc;
    let cap = capacitance_matrix(result)?;
    let mut cap_params = Map::new();
    for k in ["ground", "eliminate", "qubits", "resonators", "metal_thickness_m", "edge_cutoff_m", "c0", "c1"] {
        if let Some(v) = params.get(k).filter(|v| !v.is_null()) {
            cap_params.insert(k.to_string(), v.clone());
        }
    }

    let gram_opt = if let Some(t_val) = cap_params.get("metal_thickness_m").and_then(Value::as_f64) {
        let m = mesh.ok_or_else(|| {
            Error::invalid("metal_thickness_m in uq_cqed requires Mesh; pass uq_cqed_with_mesh_and_params")
        })?;
        Some(prepare_3d_thickness_gram(m, result, &cap, &cap_params, t_val)?)
    } else {
        None
    };

    let mut vars = Vec::new();
    let mut bindings = Vec::new();
    if let Some(ref gram) = gram_opt {
        let t_nom = cap_params.get("metal_thickness_m").and_then(Value::as_f64).unwrap_or(100.0e-9);
        if let Some(rv) = resolve_var(
            raw_unc,
            "metal_thickness_m",
            "metal_thickness_m",
            "capacitance",
            t_nom,
            1e-12,
            gram.s_eff * 0.999,
        )? {
            vars.push(rv);
            bindings.push(CapVarBinding::MetalThickness);
        }
        let c0_nom = cap_params.get("c0").and_then(Value::as_f64).unwrap_or(gram.def_c0);
        if let Some(rv) = resolve_var(raw_unc, "c0", "c0", "capacitance", c0_nom, f64::NEG_INFINITY, f64::INFINITY)? {
            vars.push(rv);
            bindings.push(CapVarBinding::C0);
        }
        let c1_nom = cap_params.get("c1").and_then(Value::as_f64).unwrap_or(gram.def_c1);
        if let Some(rv) = resolve_var(raw_unc, "c1", "c1", "capacitance", c1_nom, f64::NEG_INFINITY, f64::INFINITY)? {
            vars.push(rv);
            bindings.push(CapVarBinding::C1);
        }
    }

    if let Some(qubits) = cap_params.get("qubits").and_then(Value::as_array) {
        let multi_q = qubits.len() > 1;
        for (qi, q) in qubits.iter().enumerate() {
            if let Some(ej_nom) = q.get("E_J_Hz").and_then(Value::as_f64) {
                let qname = q.get("name").and_then(Value::as_str).unwrap_or("q");
                let key_q = format!("E_J_Hz:{qname}");
                let k = if raw_unc.contains_key(&key_q) {
                    key_q.as_str()
                } else if raw_unc.contains_key("E_J_Hz") {
                    "E_J_Hz"
                } else {
                    "ej_hz"
                };
                let label = if multi_q { format!("E_J_Hz:{qname}") } else { "E_J_Hz".to_string() };
                if let Some(rv) = resolve_var(raw_unc, k, label, "capacitance", ej_nom, 1e-12, f64::INFINITY)? {
                    vars.push(rv);
                    bindings.push(CapVarBinding::QubitEj(qi));
                }
            }
        }
    }

    if let Some(res_arr) = cap_params.get("resonators").and_then(Value::as_array) {
        let multi_r = res_arr.len() > 1;
        for (ri, r) in res_arr.iter().enumerate() {
            let rname = r.get("name").and_then(Value::as_str).unwrap_or("r");
            if let Some(fr_nom) = r.get("f_r_Hz").and_then(Value::as_f64) {
                let key_q = format!("f_r_Hz:{rname}");
                let k = if raw_unc.contains_key(&key_q) { key_q.as_str() } else { "f_r_Hz" };
                let label = if multi_r { format!("f_r_Hz:{rname}") } else { "f_r_Hz".to_string() };
                if let Some(rv) = resolve_var(raw_unc, k, label, "capacitance", fr_nom, 1.0, f64::INFINITY)? {
                    vars.push(rv);
                    bindings.push(CapVarBinding::ResFr(ri));
                }
            }
            let z0_nom = r.get("Z0_ohm").and_then(Value::as_f64).unwrap_or(50.0);
            let zq = format!("Z0_ohm:{rname}");
            let zk = if raw_unc.contains_key(&zq) {
                zq.as_str()
            } else if raw_unc.contains_key("Z0_ohm") {
                "Z0_ohm"
            } else {
                "z0_ohm"
            };
            let z_label = if multi_r { format!("Z0_ohm:{rname}") } else { "Z0_ohm".to_string() };
            if let Some(rv) = resolve_var(raw_unc, zk, z_label, "capacitance", z0_nom, 1e-6, f64::INFINITY)? {
                vars.push(rv);
                bindings.push(CapVarBinding::ResZ0(ri));
            }
            if let Some(cr_nom) = r.get("C_r_F").and_then(Value::as_f64) {
                let cq = format!("C_r_F:{rname}");
                let ck = if raw_unc.contains_key(&cq) { cq.as_str() } else { "C_r_F" };
                let c_label = if multi_r { format!("C_r_F:{rname}") } else { "C_r_F".to_string() };
                if let Some(rv) = resolve_var(raw_unc, ck, c_label, "capacitance", cr_nom, 1e-24, f64::INFINITY)? {
                    vars.push(rv);
                    bindings.push(CapVarBinding::ResCr(ri));
                }
            }
        }
    }

    let nom_vals_vec = vars.iter().map(|v| v.spec.nominal).collect::<Vec<_>>();
    let nom_eval = eval_cap_with_vars(&cap, gram_opt.as_ref(), &cap_params, &bindings, &nom_vals_vec, true)?;

    let (aleatoric_vec, clamped_cnt) =
        propagate_domain(&vars, &nom_eval.defs, &nom_eval.vals, state.seed ^ 0xCA9A_0002, state.n_samples, |x| {
            eval_cap_with_vars(&cap, gram_opt.as_ref(), &cap_params, &bindings, x, false).map(|r| r.vals)
        })?;
    if clamped_cnt > 0 {
        state
            .diagnostics
            .push(format!("capacitance: clamped {clamped_cnt} normal MC draw(s) at physical domain bounds"));
    }

    state.push_with_eb(&nom_eval.defs, &nom_eval.vals, &aleatoric_vec, result, Some(&nom_eval.epistemic_bounds));

    for rv in vars {
        state.all_var_specs.push(rv.spec);
    }
    state.domains_evaluated.push("capacitance_matrix".to_string());
    Ok(())
}

#[derive(Clone, Copy, Debug)]
enum TlsVarBinding {
    TMa,
    TMs,
    TSa,
    EpsMa,
    EpsMs,
    EpsSa,
    TanDeltaMa,
    TanDeltaMs,
    TanDeltaSa,
    TanDeltaSharedScale { base_ma: f64, base_ms: f64, base_sa: f64 },
    DFilm,
    TrenchDepth,
}

fn eval_tls_with_vars(
    prepared: &PreparedTls,
    base_params: &TlsParams,
    freq_hz_opt: Option<f64>,
    bindings: &[TlsVarBinding],
    vals: &[f64],
) -> Vec<f64> {
    let mut p = base_params.clone();
    let mut d_override = base_params.d;
    for (vi, &b) in bindings.iter().enumerate() {
        let v = vals[vi];
        match b {
            TlsVarBinding::TMa => p.t_ma = v,
            TlsVarBinding::TMs => p.t_ms = v,
            TlsVarBinding::TSa => p.t_sa = v,
            TlsVarBinding::EpsMa => p.eps_ma = v.max(1.0),
            TlsVarBinding::EpsMs => p.eps_ms = v.max(1.0),
            TlsVarBinding::EpsSa => p.eps_sa = v.max(1.0),
            TlsVarBinding::TanDeltaMa => p.tan_delta_ma = v,
            TlsVarBinding::TanDeltaMs => p.tan_delta_ms = v,
            TlsVarBinding::TanDeltaSa => p.tan_delta_sa = v,
            TlsVarBinding::TanDeltaSharedScale { base_ma, base_ms, base_sa } => {
                p.tan_delta_ma = base_ma * v;
                p.tan_delta_ms = base_ms * v;
                p.tan_delta_sa = base_sa * v;
            }
            TlsVarBinding::DFilm => d_override = Some(v),
            TlsVarBinding::TrenchDepth => p.trench_depth_m = Some(v),
        }
    }
    let rep: TlsReport = prepared.eval(&p, d_override);
    let mut out = vec![rep.p_ma, rep.p_ms, rep.p_sa, rep.p_total, rep.loss_tan_total, rep.q_tls];
    if let Some(f) = freq_hz_opt {
        let t1 = if f > 0.0 && rep.loss_tan_total > 0.0 && rep.q_tls.is_finite() {
            rep.q_tls / (2.0 * PI * f)
        } else {
            f64::INFINITY
        };
        out.push(t1);
    }
    out
}

fn tls_defs(prefix: &str, with_t1: bool) -> Vec<EvalScalarDef> {
    let mut defs = vec![
        EvalScalarDef { name: format!("{prefix}p_ma"), unit: "1".into(), non_negative: true },
        EvalScalarDef { name: format!("{prefix}p_ms"), unit: "1".into(), non_negative: true },
        EvalScalarDef { name: format!("{prefix}p_sa"), unit: "1".into(), non_negative: true },
        EvalScalarDef { name: format!("{prefix}p_total"), unit: "1".into(), non_negative: true },
        EvalScalarDef { name: format!("{prefix}loss_tan_total"), unit: "1".into(), non_negative: true },
        EvalScalarDef { name: format!("{prefix}Q_TLS"), unit: "1".into(), non_negative: true },
    ];
    if with_t1 {
        defs.push(EvalScalarDef { name: format!("{prefix}T1_TLS_s"), unit: "s".into(), non_negative: true });
    }
    defs
}

fn build_tls_vars(
    raw_unc: &BTreeMap<String, RawUncSpec>,
    p_nom: &TlsParams,
    d_default_opt: Option<f64>,
) -> Result<(Vec<ResolvedVar>, Vec<TlsVarBinding>), Error> {
    let mut vars = Vec::new();
    let mut bindings = Vec::new();

    for (key, nom, b) in [
        ("t_ma", p_nom.t_ma, TlsVarBinding::TMa),
        ("t_ms", p_nom.t_ms, TlsVarBinding::TMs),
        ("t_sa", p_nom.t_sa, TlsVarBinding::TSa),
    ] {
        if let Some(rv) = resolve_var(raw_unc, key, key, "tls", nom, 1e-12, f64::INFINITY)? {
            vars.push(rv);
            bindings.push(b);
        }
    }
    for (key, nom, b) in [
        ("eps_ma", p_nom.eps_ma, TlsVarBinding::EpsMa),
        ("eps_ms", p_nom.eps_ms, TlsVarBinding::EpsMs),
        ("eps_sa", p_nom.eps_sa, TlsVarBinding::EpsSa),
    ] {
        if let Some(rv) = resolve_var(raw_unc, key, key, "tls", nom, 1.0, f64::INFINITY)? {
            vars.push(rv);
            bindings.push(b);
        }
    }

    let has_shared_td = raw_unc.contains_key("tan_delta");
    let has_ind_td = raw_unc.contains_key("tan_delta_ma")
        || raw_unc.contains_key("tan_delta_ms")
        || raw_unc.contains_key("tan_delta_sa");
    if has_shared_td && has_ind_td {
        return Err(Error::invalid(
            "cannot specify both shared 'tan_delta' and individual 'tan_delta_{ma,ms,sa}' in uncertainties",
        ));
    }
    if has_shared_td {
        if let Some(rv) = resolve_var(raw_unc, "tan_delta", "tan_delta_scale", "tls", 1.0, 0.0, f64::INFINITY)? {
            vars.push(rv);
            bindings.push(TlsVarBinding::TanDeltaSharedScale {
                base_ma: p_nom.tan_delta_ma,
                base_ms: p_nom.tan_delta_ms,
                base_sa: p_nom.tan_delta_sa,
            });
        }
    } else {
        for (key, nom, b) in [
            ("tan_delta_ma", p_nom.tan_delta_ma, TlsVarBinding::TanDeltaMa),
            ("tan_delta_ms", p_nom.tan_delta_ms, TlsVarBinding::TanDeltaMs),
            ("tan_delta_sa", p_nom.tan_delta_sa, TlsVarBinding::TanDeltaSa),
        ] {
            if let Some(rv) = resolve_var(raw_unc, key, key, "tls", nom, 0.0, f64::INFINITY)? {
                vars.push(rv);
                bindings.push(b);
            }
        }
    }

    if let Some(d_nom) = p_nom.d.or(d_default_opt) {
        let dk = if raw_unc.contains_key("d_film_m") { "d_film_m" } else { "d" };
        if let Some(rv) = resolve_var(raw_unc, dk, "d", "tls", d_nom, 1e-12, f64::INFINITY)? {
            vars.push(rv);
            bindings.push(TlsVarBinding::DFilm);
        }
    }
    let trench_nom = p_nom.trench_depth_m.unwrap_or(0.0);
    if d_default_opt.is_some() {
        if let Some(rv) =
            resolve_var(raw_unc, "trench_depth_m", "trench_depth_m", "tls", trench_nom, 0.0, f64::INFINITY)?
        {
            vars.push(rv);
            bindings.push(TlsVarBinding::TrenchDepth);
        }
    }

    Ok((vars, bindings))
}

fn push_tls_mode_uq(
    prepared: &PreparedTls,
    p_nom: &TlsParams,
    freq_hz_opt: Option<f64>,
    prefix: &str,
    eb_result: &SolveResult,
    mode_seed: u64,
    state: &mut UqRunState<'_>,
) -> Result<(), Error> {
    let (vars, bindings) = build_tls_vars(state.raw_unc, p_nom, prepared.default_d_film())?;
    let defs = tls_defs(prefix, freq_hz_opt.is_some());
    let nom_x = vars.iter().map(|v| v.spec.nominal).collect::<Vec<_>>();
    let nom_vals = eval_tls_with_vars(prepared, p_nom, freq_hz_opt, &bindings, &nom_x);

    let (aleatoric_vec, clamped_cnt) = propagate_domain(&vars, &defs, &nom_vals, mode_seed, state.n_samples, |x| {
        Ok(eval_tls_with_vars(prepared, p_nom, freq_hz_opt, &bindings, x))
    })?;
    if clamped_cnt > 0 {
        state
            .diagnostics
            .push(format!("tls({prefix}): clamped {clamped_cnt} normal MC draw(s) at physical domain bounds"));
    }

    state.push_with_eb(&defs, &nom_vals, &aleatoric_vec, eb_result, None);

    for rv in vars {
        if !state.all_var_specs.iter().any(|s| s.domain == rv.spec.domain && s.name == rv.spec.name) {
            state.all_var_specs.push(rv.spec);
        }
    }
    Ok(())
}

fn run_tls_uq(
    mesh: &Mesh,
    result: &SolveResult,
    companion_mesh_and_res: Option<(&Mesh, &SolveResult)>,
    params: &Value,
    state: &mut UqRunState<'_>,
) -> Result<(), Error> {
    let mut tls_json = Map::new();
    for k in [
        "solution",
        "t_ma",
        "t_ms",
        "t_sa",
        "eps_ma",
        "eps_ms",
        "eps_sa",
        "tan_delta_ma",
        "tan_delta_ms",
        "tan_delta_sa",
        "d",
        "edge_cutoff",
        "trench_depth_m",
    ] {
        if let Some(v) = params.get(k).filter(|v| !v.is_null()) {
            tls_json.insert(k.to_string(), v.clone());
        }
    }
    let tls_val = Value::Object(tls_json);

    if let (Some((e_mesh, e_res)), Some(mv_val)) =
        (companion_mesh_and_res, params.get("modal_voltages").filter(|v| !v.is_null()))
    {
        let modal_voltages: Vec<Vec<f64>> = serde_json::from_value(mv_val.clone())
            .map_err(|e| Error::invalid(format!("invalid modal_voltages: {e}")))?;
        let modal_freqs_hz: Vec<f64> = match params.get("modal_freqs_hz").filter(|v| !v.is_null()) {
            Some(fv) => serde_json::from_value(fv.clone())
                .map_err(|e| Error::invalid(format!("invalid modal_freqs_hz: {e}")))?,
            None => vec![0.0; modal_voltages.len()],
        };
        let modes = prepare_quasi_electrostatic_tls(e_mesh, e_res, &modal_voltages, &modal_freqs_hz, &tls_val)?;
        let multi = modes.len() > 1;
        for m in &modes {
            let prefix = if multi { format!("m{}.", m.mode) } else { String::new() };
            push_tls_mode_uq(
                &m.prepared,
                &m.params,
                (m.freq_hz > 0.0).then_some(m.freq_hz),
                &prefix,
                e_res,
                state.seed ^ (0x7150_0000 + m.mode as u64),
                state,
            )?;
        }
        state.domains_evaluated.push("tls_quasi_electrostatic".to_string());
        return Ok(());
    }

    // Standard single-field TLS on primary or companion electrostatic result
    let (target_mesh, target_res) = match companion_mesh_and_res {
        Some((em, er)) if tls_participation_applies(er) && !tls_participation_applies(result) => (em, er),
        _ => (mesh, result),
    };
    let (prepared, p_nom) = prepare_tls(target_mesh, target_res, &tls_val)?;
    let _ = TlsMeta::of(target_res);
    let freq_hz_opt = CircuitRom::of(target_res)
        .and_then(|rom| rom.modes.get(p_nom.solution).map(|m| m.freq_hz))
        .or_else(|| params.get("freq_Hz").and_then(Value::as_f64))
        .or_else(|| target_res.solutions.get(p_nom.solution).and_then(|s| s.value))
        .filter(|&f| f.is_finite() && f > 1.0e5);

    push_tls_mode_uq(&prepared, &p_nom, freq_hz_opt, "", target_res, state.seed ^ 0x7153_0003, state)?;
    state.domains_evaluated.push("tls_participation".to_string());
    Ok(())
}

/// Runs zero-solve cQED post-processing UQ with default parameters (nominal + any epistemic bounds).
pub fn uq_cqed(mesh: &Mesh, result: &SolveResult) -> Result<UqCqedReport, Error> {
    uq_cqed_with_mesh_and_params(mesh, result, &Value::Null)
}

/// Runs zero-solve cQED post-processing UQ without a mesh (ROM and CapacitanceMatrix without 3D thickness).
pub fn uq_cqed_with_params(result: &SolveResult, params: &Value) -> Result<UqCqedReport, Error> {
    uq_cqed_impl(None, result, None, params)
}

/// Runs zero-solve cQED post-processing UQ with a mesh (enabling TLS and 3D thickness UQ).
pub fn uq_cqed_with_mesh_and_params(mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<UqCqedReport, Error> {
    uq_cqed_impl(Some(mesh), result, None, params)
}

/// Runs joint zero-solve cQED UQ across a primary `SolveResult` (e.g., `ems` `CircuitRom`) paired with a
/// companion electrostatic `SolveResult` (e.g., `dfs` `CapacitanceMatrixData` with Prager–Synge bounds and/or
/// electric fields for quasi-electrostatic TLS).
pub fn uq_cqed_with_companion_and_params(
    primary_mesh: Option<&Mesh>,
    primary_result: &SolveResult,
    companion: Option<(&Mesh, &SolveResult)>,
    params: &Value,
) -> Result<UqCqedReport, Error> {
    uq_cqed_impl(primary_mesh, primary_result, companion, params)
}

fn uq_cqed_impl(
    primary_mesh: Option<&Mesh>,
    result: &SolveResult,
    companion: Option<(&Mesh, &SolveResult)>,
    params: &Value,
) -> Result<UqCqedReport, Error> {
    let t0 = Instant::now();
    if !params.is_null() && !params.is_object() {
        return Err(Error::invalid("uq_cqed params must be a JSON object or null"));
    }

    let seed = match params.get("seed").filter(|v| !v.is_null()) {
        None => DEFAULT_SEED,
        Some(v) => v.as_u64().ok_or_else(|| Error::invalid("parameter 'seed' must be a non-negative integer"))?,
    };
    let n_samples = match params.get("n_samples").filter(|v| !v.is_null()) {
        None => DEFAULT_N_SAMPLES,
        Some(v) => {
            let n = v.as_u64().ok_or_else(|| Error::invalid("parameter 'n_samples' must be a non-negative integer"))?
                as usize;
            if n > MAX_N_SAMPLES {
                return Err(Error::invalid(format!(
                    "parameter 'n_samples' ({n}) exceeds maximum allowed ({MAX_N_SAMPLES})"
                )));
            }
            n
        }
    };

    let raw_unc = parse_raw_uncertainties(params)?;
    let has_rom = epr_applies(result);
    let has_cap = capacitance_matrix_applies(result);
    let has_tls_primary = tls_participation_applies(result);
    let has_tls_companion = companion.is_some_and(|(_, er)| tls_participation_applies(er))
        && (params.get("modal_voltages").is_some_and(|v| !v.is_null())
            || params.get("include_companion_tls").and_then(Value::as_bool) == Some(true));

    if !has_rom && !has_cap && !has_tls_primary && !has_tls_companion {
        return Err(Error::invalid(
            "the result carries no CircuitRom, CapacitanceMatrixData, or electric field suitable for uq_cqed",
        ));
    }

    let mut all_var_specs = Vec::new();
    let mut quantities = Vec::new();
    let mut domains_evaluated = Vec::new();
    let mut diagnostics = Vec::new();
    {
        let mut state = UqRunState {
            raw_unc: &raw_unc,
            seed,
            n_samples,
            all_var_specs: &mut all_var_specs,
            out_quantities: &mut quantities,
            domains_evaluated: &mut domains_evaluated,
            diagnostics: &mut diagnostics,
        };
        if has_rom {
            run_rom_uq(result, companion.map(|(_, er)| er), params, &mut state)?;
        }
        if has_cap {
            run_capacitance_uq(primary_mesh, result, params, &mut state)?;
        }
        if has_tls_primary || has_tls_companion {
            let m_for_tls = primary_mesh.or(companion.map(|(m, _)| m)).ok_or_else(|| {
                Error::invalid("TLS uncertainty quantification requires a Mesh; use uq_cqed_with_mesh_and_params")
            })?;
            run_tls_uq(m_for_tls, result, companion, params, &mut state)?;
        }
    }

    Ok(UqCqedReport {
        seed,
        n_samples,
        zero_fem_solves: true,
        epistemic_aleatoric_separated: true,
        domains_evaluated,
        uncertain_parameters: all_var_specs,
        quantities,
        diagnostics,
        elapsed_s: t0.elapsed().as_secs_f64(),
    })
}

/// Zero-solve cQED post-processing uncertainty quantification behind the [`Analysis`] trait (`id: "uq_cqed"`).
pub struct UqCqedAnalysis;

impl Analysis for UqCqedAnalysis {
    fn id(&self) -> &'static str {
        "uq_cqed"
    }

    fn name(&self) -> &'static str {
        "Zero-solve cQED uncertainty quantification (ASME V&V 20 epistemic/aleatoric separated)"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        uq_cqed_applies(result)
    }

    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        Ok(uq_cqed(mesh, result)?.to_json())
    }

    fn run_with_params(&self, mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        Ok(uq_cqed_with_mesh_and_params(mesh, result, params)?.to_json())
    }
}
