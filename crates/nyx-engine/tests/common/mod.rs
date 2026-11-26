// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Shared Phase 4 ASME V&V 20 validation harness, dataset loaders, numerical uncertainty
//! (`u_num`) calculators, and quantitative example runners.

#![allow(dead_code)]

use femcore::vvuq::{roache_gci, RoacheGciRegime, RoacheGciReport};
use nyx_engine::Engine;
use serde::de::DeserializeOwned;
use serde::Deserialize;
use serde_json::{json, Value};
use std::path::PathBuf;

/// Default ASME V&V 20 coverage factor (`k = 2.0` -> ~95% coverage for aleatoric budget).
fn default_coverage_k() -> f64 {
    2.0
}

/// Versioned reference validation target loaded from `benchmarks/validation/*.json`.
#[derive(Clone, Debug, Deserialize)]
pub struct ValidationTarget {
    pub id: String,
    pub description: String,
    pub quantity: String,
    pub unit: String,
    pub reference_value: f64,
    pub u_input: f64,
    pub u_d: f64,
    pub tol_rel: f64,
    pub tol_abs: f64,
    #[serde(default = "default_coverage_k")]
    pub coverage_k: f64,
}

/// Detailed quantitative outcome of evaluating the explicit ASME V&V 20 `AND` rule:
/// `|E| <= tol_rel * |D| + tol_abs  AND  |E| <= coverage_k * sqrt(u_num^2 + u_input^2 + u_D^2)`.
#[derive(Clone, Debug)]
pub struct AsmeVv20Assessment {
    pub id: String,
    pub s: f64,
    pub d: f64,
    pub e: f64,
    pub abs_e: f64,
    pub rel_e: f64,
    pub u_num: f64,
    pub u_input: f64,
    pub u_d: f64,
    pub u_aleatoric: f64,
    pub u_val: f64,
    pub tol: f64,
    pub coverage_k: f64,
    pub passed_accuracy: bool,
    pub passed_uncertainty: bool,
    pub passed: bool,
}

/// Strict separation outcome when numerical uncertainty is a deterministic Prager-Synge
/// epistemic interval `[s_low, s_up]` (`delta_ep` is NEVER folded into Gaussian `sqrt` variances).
#[derive(Clone, Debug)]
pub struct EpistemicVv20Assessment {
    pub id: String,
    pub s_low: f64,
    pub s_up: f64,
    pub s_mid: f64,
    pub d: f64,
    pub delta_ep: f64,
    pub u_aleatoric: f64,
    pub abs_e_mid: f64,
    pub tol: f64,
    pub passed_accuracy: bool,
    pub passed_interval: bool,
    pub passed_bracket: bool,
    pub passed: bool,
}

/// Returns the root folder holding versioned `benchmarks/validation/*.json` datasets.
pub fn validation_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../benchmarks/validation")
}

/// Loads and deserializes a JSON reference dataset under `benchmarks/validation/`.
pub fn load_validation_json<T: DeserializeOwned>(filename: &str) -> T {
    let path = validation_dir().join(filename);
    let raw = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read validation dataset {}: {e}", path.display()));
    serde_json::from_str(&raw).unwrap_or_else(|e| panic!("failed to parse validation dataset {}: {e}", path.display()))
}

/// Builds a `ValidationTarget` inline for analytical closed-form checks (e.g. upgraded examples).
#[allow(clippy::too_many_arguments)]
pub fn inline_target(
    id: impl Into<String>,
    description: impl Into<String>,
    quantity: impl Into<String>,
    unit: impl Into<String>,
    reference_value: f64,
    u_input: f64,
    u_d: f64,
    tol_rel: f64,
    tol_abs: f64,
) -> ValidationTarget {
    ValidationTarget {
        id: id.into(),
        description: description.into(),
        quantity: quantity.into(),
        unit: unit.into(),
        reference_value,
        u_input,
        u_d,
        tol_rel,
        tol_abs,
        coverage_k: 2.0,
    }
}

/// Verifies the explicit ASME V&V 20 `AND` validation rule with no uncertainty-inflation loopholes:
///
/// 1. `E = S - D`, `u_aleatoric = sqrt(u_input^2 + u_D^2)`, `u_val = sqrt(u_num^2 + u_input^2 + u_D^2)`
/// 2. `u_num`, `u_input`, `u_D` must be finite and non-negative, and `u_val > 0`.
/// 3. **Domain accuracy ceiling**: `|E| <= tol_rel * |D| + tol_abs`
/// 4. **Uncertainty consistency**: `|E| <= coverage_k * u_val` (`1.0 <= coverage_k <= 2.0`)
///
/// Both clauses must hold simultaneously (`AND`).
pub fn verify_asme_vv20(target: &ValidationTarget, s: f64, u_num: f64) -> AsmeVv20Assessment {
    let d = target.reference_value;
    assert!(s.is_finite(), "[VV20 FAIL {}] non-finite S = {s}", target.id);
    assert!(d.is_finite(), "[VV20 FAIL {}] non-finite D = {d}", target.id);
    assert!(u_num.is_finite() && u_num >= 0.0, "[VV20 FAIL {}] invalid u_num = {u_num}", target.id);
    assert!(
        target.u_input.is_finite() && target.u_input >= 0.0,
        "[VV20 FAIL {}] invalid u_input = {}",
        target.id,
        target.u_input
    );
    assert!(target.u_d.is_finite() && target.u_d >= 0.0, "[VV20 FAIL {}] invalid u_D = {}", target.id, target.u_d);
    assert!(
        target.tol_rel.is_finite() && target.tol_rel >= 0.0 && target.tol_abs.is_finite() && target.tol_abs >= 0.0,
        "[VV20 FAIL {}] invalid tolerances tol_rel={}, tol_abs={}",
        target.id,
        target.tol_rel,
        target.tol_abs
    );
    assert!(
        target.coverage_k >= 1.0 && target.coverage_k <= 2.0,
        "[VV20 FAIL {}] coverage_k={} must lie in [1.0, 2.0]",
        target.id,
        target.coverage_k
    );

    let e = s - d;
    let abs_e = e.abs();
    let rel_e = abs_e / d.abs().max(1e-30);
    let u_aleatoric = target.u_input.hypot(target.u_d);
    let u_val = u_num.hypot(u_aleatoric);
    assert!(
        u_val > 0.0,
        "[VV20 FAIL {}] u_val must be strictly positive (u_num={u_num:e}, u_input={:e}, u_D={:e})",
        target.id,
        target.u_input,
        target.u_d
    );

    let tol = target.tol_rel * d.abs() + target.tol_abs;
    let u_bound = target.coverage_k * u_val;
    let passed_accuracy = abs_e <= tol;
    let passed_uncertainty = abs_e <= u_bound;
    let passed = passed_accuracy && passed_uncertainty;

    let status = if passed { "PASS" } else { "FAIL" };
    println!(
        "[VV20 {status} {}] {}: S={s:.6e} {}, D={d:.6e} {}, E={e:+.3e} ({:.4}%), u_num={u_num:.3e}, u_inp={:.3e}, u_D={:.3e}, u_val={u_val:.3e}, k*u_val={u_bound:.3e}, tol={tol:.3e}",
        target.id,
        target.quantity,
        target.unit,
        target.unit,
        rel_e * 100.0,
        target.u_input,
        target.u_d,
    );

    assert!(
        passed,
        "[VV20 FAIL {}] AND rule violated (accuracy={passed_accuracy} [|E|={abs_e:.4e} vs tol={tol:.4e}], uncertainty={passed_uncertainty} [|E|={abs_e:.4e} vs {:.1}*u_val={u_bound:.4e}], u_num={u_num:.4e}, u_val={u_val:.4e})",
        target.id,
        target.coverage_k
    );

    AsmeVv20Assessment {
        id: target.id.clone(),
        s,
        d,
        e,
        abs_e,
        rel_e,
        u_num,
        u_input: target.u_input,
        u_d: target.u_d,
        u_aleatoric,
        u_val,
        tol,
        coverage_k: target.coverage_k,
        passed_accuracy,
        passed_uncertainty,
        passed,
    }
}

/// Verifies an epistemic Prager-Synge dual interval `[s_low, s_up]` against reference `D`
/// using the §6.1 strict epistemic/aleatoric separation rule (`AND`):
///
/// 1. `s_low <= s_up`, `s_mid = 0.5*(s_low + s_up)`, `delta_ep = 0.5*(s_up - s_low)`
/// 2. `u_aleatoric = sqrt(u_input^2 + u_D^2)` (`delta_ep` is NEVER folded into variance)
/// 3. `|s_mid - D| <= tol_rel * |D| + tol_abs` **AND** `|s_mid - D| <= delta_ep + k * u_aleatoric`
/// 4. When `require_strict_bracket == true`, verifies `s_low - k*u_aleatoric <= D <= s_up + k*u_aleatoric`.
pub fn verify_epistemic_vv20(
    target: &ValidationTarget,
    s_low: f64,
    s_up: f64,
    require_strict_bracket: bool,
) -> EpistemicVv20Assessment {
    let d = target.reference_value;
    assert!(
        s_low.is_finite() && s_up.is_finite() && s_low <= s_up,
        "[VV20-EP FAIL {}] invalid interval [{s_low:e}, {s_up:e}]",
        target.id
    );
    let s_mid = 0.5 * (s_low + s_up);
    let delta_ep = 0.5 * (s_up - s_low);
    let u_aleatoric = target.u_input.hypot(target.u_d);
    let tol = target.tol_rel * d.abs() + target.tol_abs;
    let abs_e_mid = (s_mid - d).abs();
    let margin = target.coverage_k * u_aleatoric;

    let passed_accuracy = abs_e_mid <= tol;
    let passed_interval = abs_e_mid <= delta_ep + margin;
    let passed_bracket = !require_strict_bracket || (d >= s_low - margin && d <= s_up + margin);
    let passed = passed_accuracy && passed_interval && passed_bracket;

    let status = if passed { "PASS" } else { "FAIL" };
    println!(
        "[VV20-EP {status} {}] {}: [S_low, S_up]=[{s_low:.6e}, {s_up:.6e}] {}, S_mid={s_mid:.6e}, D={d:.6e}, delta_ep={delta_ep:.3e}, u_ale={u_aleatoric:.3e}, |E_mid|={abs_e_mid:.3e} ({:.4}%), tol={tol:.3e}",
        target.id,
        target.quantity,
        target.unit,
        (abs_e_mid / d.abs().max(1e-30)) * 100.0,
    );

    assert!(
        passed,
        "[VV20-EP FAIL {}] epistemic AND rule failed (acc={passed_accuracy}, interval={passed_interval}, bracket={passed_bracket}): mid={s_mid:.4e} vs D={d:.4e}, [low,up]=[{s_low:.4e},{s_up:.4e}], delta_ep={delta_ep:.3e}, u_ale={u_aleatoric:.3e}, tol={tol:.3e}",
        target.id
    );

    EpistemicVv20Assessment {
        id: target.id.clone(),
        s_low,
        s_up,
        s_mid,
        d,
        delta_ep,
        u_aleatoric,
        abs_e_mid,
        tol,
        passed_accuracy,
        passed_interval,
        passed_bracket,
        passed,
    }
}

/// Computes `u_num` in physical units from 3 ordered meshes (`phi1` finest -> `phi3` coarsest)
/// using `femcore::vvuq::roache_gci`.
pub fn u_num_three_grid(
    phi1_fine: f64,
    phi2_mid: f64,
    phi3_coarse: f64,
    r21: f64,
    r32: f64,
    p_formal: f64,
) -> (f64, RoacheGciReport) {
    let report =
        roache_gci([phi1_fine, phi2_mid, phi3_coarse], r21, r32, p_formal).expect("valid three-grid GCI inputs");
    let scale = phi1_fine.abs().max(1e-30);
    let u_rel = match report.regime {
        RoacheGciRegime::Monotonic | RoacheGciRegime::Oscillatory => report.gci_21,
        RoacheGciRegime::NonConvergent => report.gci_21.max((phi2_mid - phi1_fine).abs() / scale),
    };
    let u_num = (u_rel * scale).max(1e-12 * scale);
    (u_num, report)
}

/// Computes conservative 2-grid Roache GCI uncertainty (`F_s = 3.0`) from `(s_coarse, s_fine)`.
pub fn u_num_two_grid(s_coarse: f64, s_fine: f64, r_eff: f64, p_order: f64) -> f64 {
    assert!(
        r_eff.is_finite() && r_eff > 1.0 && p_order.is_finite() && p_order > 0.0,
        "invalid refinement ratio r_eff={r_eff} or order p_order={p_order}"
    );
    let denom = (r_eff.powf(p_order) - 1.0).max(0.1);
    let gci = 3.0 * (s_fine - s_coarse).abs() / denom;
    gci.max(1e-12 * s_fine.abs().max(1.0e-30))
}

/// Helper to execute an `Engine::call` and panic cleanly with context on error.
pub fn cmd(e: &mut Engine, name: &str, args: Value) -> Value {
    e.call(name, &args).unwrap_or_else(|err| panic!("engine command '{name}' failed: {err}"))
}

/// Looks up the first numeric scalar whose name equals `exact_name`.
pub fn scalar_exact(scalars: &[(String, nyx_physics::Scalar)], exact_name: &str) -> Option<f64> {
    scalars.iter().find_map(|(k, v)| match (k == exact_name, v) {
        (true, nyx_physics::Scalar::Num(x)) => Some(*x),
        _ => None,
    })
}

/// Looks up the first numeric scalar whose name starts with `prefix`.
pub fn scalar_prefix(scalars: &[(String, nyx_physics::Scalar)], prefix: &str) -> Option<f64> {
    scalars.iter().find_map(|(k, v)| match (k.starts_with(prefix), v) {
        (true, nyx_physics::Scalar::Num(x)) => Some(*x),
        _ => None,
    })
}

/// Looks up a scalar inside a JSON `scalars` array (`[[name, value], ...]`) by exact match.
pub fn json_scalar(scalars: &Value, name: &str) -> Option<f64> {
    scalars.as_array()?.iter().find_map(|pair| {
        let arr = pair.as_array()?;
        (arr.first()?.as_str()? == name).then(|| arr.get(1)?.as_f64()).flatten()
    })
}

/// Runs `convergence_check` on a solved `Engine` physics, returning `(refined_value, u_num)`
/// from either 3-grid Roache GCI or 2-grid conservative delta for the specified `qty`.
pub fn engine_convergence_u_num(
    e: &mut Engine,
    physics: &str,
    factor: f64,
    grids: usize,
    qty: &str,
    track: bool,
) -> (f64, f64, Value) {
    let mut args = json!({
        "physics": physics,
        "factor": factor,
        "grids": grids,
        "tolerance": 0.5,
    });
    if track {
        args["match"] = json!("track");
    }
    let cc = cmd(e, "convergence_check", args);
    let rows = cc["rows"].as_array().expect("convergence_check rows");
    let row = rows
        .iter()
        .find(|r| r["quantity"].as_str() == Some(qty))
        .unwrap_or_else(|| panic!("quantity '{qty}' not found in convergence_check rows: {cc}"));
    let refined = row["refined"].as_f64().expect("refined value");
    let u_num = if grids == 3 {
        if let Some(gci_rel) = row.get("gci_21").and_then(Value::as_f64) {
            (gci_rel * refined.abs().max(1e-30)).max(1e-12 * refined.abs().max(1e-30))
        } else {
            let coarse = row["coarse"].as_f64().expect("coarse value");
            let r_eff = cc["r_eff"].as_f64().unwrap_or(1.0 / factor).max(1.05);
            u_num_two_grid(coarse, refined, r_eff, 1.0)
        }
    } else {
        let coarse = row["coarse"].as_f64().expect("coarse value");
        let r_eff = cc["r_eff"].as_f64().unwrap_or(1.0 / factor).max(1.05);
        u_num_two_grid(coarse, refined, r_eff, 1.5)
    };
    (refined, u_num, cc)
}

/// Path to `examples/<name>.nyx.json`.
pub fn example_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples").join(format!("{name}.nyx.json"))
}

/// Opens, validates, meshes, and solves an example project through `Engine::call`.
pub fn solved_example(name: &str) -> (Engine, Value) {
    let mut e = Engine::new();
    let path = example_file(name);
    cmd(&mut e, "project_open", json!({ "path": path.to_str().unwrap() }));
    cmd(&mut e, "check_project", json!({}));
    let reply = cmd(&mut e, "solve_run", json!({}));
    (e, reply)
}
