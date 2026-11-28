// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Active-study material provenance uncertainty and temperature-range guardrails (VVUQ Phase 4).
//!
//! Runs lightweight preflight and post-solve material audits over the active problem and its solved
//! temperature fields without adding any duplicate evaluation logic to solver kernels:
//! - surfaces `Source::flags_uncertainty()` (`material_approximate`, `material_unverified`) and
//!   out-of-range clamping (`temperature_clamped`) in `SolveResult::diagnostics.warnings`;
//! - upgrades temperature-range clamping to a refusal error when `strict_material_range` is enabled.

use nyx_physics::{
    audit_active_materials, Ctx, DiagnosticSeverity, DiagnosticWarning, Error, MaterialAudit, Physics, ProblemSpec,
    SolveResult, SolverDiagnostics,
};
use serde_json::Value;
use std::collections::BTreeSet;

/// Collects explicit active temperatures declared on `ctx` or `spec` (study params, uniform thermal load, imported fields, or heat BCs).
pub(crate) fn collect_pre_solve_temperatures(spec: &ProblemSpec, ctx: &Ctx) -> Vec<f64> {
    let mut temps = Vec::new();
    let mut push_finite = |t: f64| {
        if t.is_finite() && !temps.iter().any(|x: &f64| x.to_bits() == t.to_bits()) {
            temps.push(t);
        }
    };
    if let Some(t) = ctx.temperature {
        push_finite(t);
    }
    if let Some(t) = spec.params.get("temperature").and_then(Value::as_f64) {
        push_finite(t);
    } else if spec.params.get("thermal").and_then(Value::as_str) == Some("uniform") {
        // Sms uniform thermal expansion defaults `temperature` to 293.15 K when omitted.
        push_finite(293.15);
    }
    if let Some(imported) = ctx.inputs.get("temperature") {
        push_field_extrema(imported, &mut push_finite);
    }
    if matches!(spec.study.kind.as_str(), "heat_stationary" | "heat_transient") {
        if let Some(t0) = spec.study.params.get("t0").and_then(parse_numeric_value) {
            push_finite(t0);
        }
        let mut has_farfield = false;
        for bc in &spec.bcs {
            match bc.kind.as_str() {
                "temperature" => {
                    if let Some(v) = bc.params.get("value").and_then(Value::as_f64) {
                        push_finite(v);
                    }
                }
                "convection" => {
                    if let Some(v) = bc.params.get("t_ext").and_then(Value::as_f64) {
                        push_finite(v);
                    }
                }
                "radiation" => {
                    if let Some(v) = bc.params.get("t_ambient").and_then(Value::as_f64) {
                        push_finite(v);
                    }
                }
                "far_field" => {
                    has_farfield = true;
                    if let Some(v) = bc.params.get("t_inf").and_then(Value::as_f64) {
                        push_finite(v);
                    }
                }
                "infinite_element" => {
                    has_farfield = true;
                }
                _ => {}
            }
        }
        if has_farfield {
            if let Some(t_inf) = spec.study.params.get("farfield_t_inf").and_then(Value::as_f64) {
                push_finite(t_inf);
            }
        }
    }
    temps
}

fn parse_numeric_value(v: &Value) -> Option<f64> {
    v.as_f64().or_else(|| v.as_str().and_then(|s| s.trim().parse::<f64>().ok()))
}

fn push_field_extrema(values: &[f64], push_finite: &mut impl FnMut(f64)) {
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;
    for &t in values {
        if t.is_finite() {
            t_min = t_min.min(t);
            t_max = t_max.max(t);
        }
    }
    if t_min.is_finite() && t_max.is_finite() {
        push_finite(t_min);
        push_finite(t_max);
    }
}

/// Checks pre-solve active materials and temperatures, refusing immediately if `strict_material_range` is enabled and any curve clamps.
pub(crate) fn preflight(
    p: &dyn Physics,
    spec: &ProblemSpec,
    ctx: &Ctx,
    strict_material_range: bool,
) -> Result<MaterialAudit, Error> {
    let props = p.properties();
    let temps = collect_pre_solve_temperatures(spec, ctx);
    let audit = audit_active_materials(spec, ctx, &props, &temps)?;
    if strict_material_range {
        if let Some(err) = audit.strict_error() {
            return Err(err);
        }
    }
    Ok(audit)
}

/// Checks solved thermal fields (if any), enforces `strict_material_range`, and appends deduplicated material warnings into `result.diagnostics.warnings`.
pub(crate) fn finalize(
    result: &mut SolveResult,
    p: &dyn Physics,
    spec: &ProblemSpec,
    ctx: &Ctx,
    strict_material_range: bool,
    pre: MaterialAudit,
) -> Result<(), Error> {
    let mut solved_temps = Vec::new();
    if matches!(spec.study.kind.as_str(), "heat_stationary" | "heat_transient") {
        let mut push_finite = |t: f64| {
            if t.is_finite() && !solved_temps.iter().any(|x: &f64| x.to_bits() == t.to_bits()) {
                solved_temps.push(t);
            }
        };
        for sol in &result.solutions {
            if let Some(tf) = sol.field("T") {
                push_field_extrema(&tf.values, &mut push_finite);
            }
        }
    }
    let post = if solved_temps.is_empty() {
        None
    } else {
        let props = p.properties();
        Some(audit_active_materials(spec, ctx, &props, &solved_temps)?)
    };
    if strict_material_range {
        if let Some(err) = post.as_ref().and_then(MaterialAudit::strict_error) {
            return Err(err);
        }
    }

    let warnings = collect_deduped_warnings(&pre, post.as_ref());
    if !warnings.is_empty() {
        let diag = result.diagnostics.get_or_insert_with(|| SolverDiagnostics {
            solver_kind: format!("{}.{}", spec.physics, spec.study.kind),
            converged: true,
            iterations: None,
            restarts: None,
            max_algebraic_residual: 0.0,
            per_mode_residuals: Vec::new(),
            rcond_estimate: None,
            cholesky_fallbacks: 0,
            warnings: Vec::new(),
        });
        let mut existing: BTreeSet<(String, String)> =
            diag.warnings.iter().map(|w| (w.code.clone(), w.message.clone())).collect();
        for w in warnings {
            if existing.insert((w.code.clone(), w.message.clone())) {
                diag.warnings.push(DiagnosticWarning {
                    code: w.code,
                    severity: DiagnosticSeverity::Warning,
                    message: w.message,
                });
            }
        }
    }
    Ok(())
}

fn collect_deduped_warnings(pre: &MaterialAudit, post: Option<&MaterialAudit>) -> Vec<DiagnosticWarning> {
    let mut seen: BTreeSet<(String, String)> = BTreeSet::new();
    let mut out = Vec::new();
    let extra = post.into_iter().flat_map(|a| a.all_warnings());
    for w in pre.all_warnings().into_iter().chain(extra) {
        if seen.insert((w.code.clone(), w.message.clone())) {
            out.push(w);
        }
    }
    out
}
