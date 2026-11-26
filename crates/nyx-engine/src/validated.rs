// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Validated is a state, not a feeling (productivity plan P6, ADR 91). A result is `exploratory` until a mesh-convergence check with a stated
//! tolerance has passed on **the same problem key** (`Engine::problem_key`, ADR 84), then it is `checked`. A third state says why a result
//! that is not checked cannot be: `check_not_affordable` when the check, which solves again on a finer mesh, is expected not to fit this
//! machine's memory (the estimate of ADR 86). A stale result is never `checked`: its problem is not the project's any more.
//!
//! The checks the engine has seen are kept by problem key, not by result, so a result solved again from the same problem (another job, a file
//! opened) is checked as well; a check that did not pass is kept too and its numbers are the reason a result stays exploratory. The verdict
//! also goes into the run record that solved the problem when the engine writes runs (`Convergence` of `nyx_store::runs`).
//!
//! Limits, stated: the check compares every eigenvalue, value and numeric scalar of the result, matching solutions by index (the order the
//! solver returns them; tracking by overlap, plan R11, is not built), and a physics that takes fields from another is checked the way
//! `convergence_check` solves it, on its own. Cheaper evidence (the residual estimator of ADR 53, a ladder of two meshes) is not stored
//! anywhere yet, so it cannot count as partial support; `check_not_affordable` therefore says only what is missing.

use crate::engine::Engine;
use nyx_physics::Error;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

/// The element-size ratio a check is planned with when none was made (the step the command's own help recommends).
pub const PLANNED_FACTOR: f64 = 0.6;

fn default_grids() -> usize {
    2
}

/// What a mesh-convergence check with a tolerance found, for one problem key.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Check {
    pub factor: f64,
    pub tolerance: f64,
    pub matching: String,
    pub passed: bool,
    pub worst_relative: f64,
    pub worst_quantity: String,
    pub quantities: usize,
    pub cells_coarse: Option<u64>,
    pub cells_refined: Option<u64>,
    pub seconds: f64,
    #[serde(default = "default_grids")]
    pub grids: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cells_medium: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r_eff: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub observed_order: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gci_21: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gci_regime: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_algebraic_residual: Option<f64>,
}

impl Check {
    /// What was compared and how, in words (the `observable` of the run record).
    pub fn observable(&self) -> String {
        let how = match self.matching.as_str() {
            "nearest" => "solutions paired by nearest eigenvalue",
            "track" => "solutions paired by tracked mode",
            _ => "solutions matched by index",
        };
        let grid_desc = if self.grids >= 3 { ", 3-grid GCI" } else { "" };
        format!(
            "{} quantit{} (eigenvalues, values, numeric scalars), {how}{grid_desc}",
            self.quantities,
            if self.quantities == 1 { "y" } else { "ies" }
        )
    }

    /// The part of a convergence answer that is a verdict; `None` for a check that had no tolerance (it says nothing is passed or failed).
    pub(crate) fn from_answer(v: &Value) -> Option<Check> {
        let tolerance = v["tolerance"].as_f64()?;
        let rows = v["rows"].as_array()?;
        let matching = v
            .get("match")
            .and_then(Value::as_str)
            .or_else(|| v.get("matching").and_then(Value::as_str))
            .unwrap_or("index")
            .to_string();
        let grids = v.get("grids").and_then(Value::as_u64).unwrap_or(2) as usize;
        let r_eff = v
            .get("r_eff")
            .and_then(Value::as_f64)
            .or_else(|| v.get("gci").and_then(|g| g.get("r_eff")).and_then(Value::as_f64));
        let gci = v.get("gci");
        let observed_order = gci.and_then(|g| g.get("observed_order")).and_then(Value::as_f64);
        let gci_21 = gci.and_then(|g| g.get("gci_21")).and_then(Value::as_f64);
        let gci_regime = gci.and_then(|g| g.get("regime")).and_then(Value::as_str).map(str::to_string);
        let max_algebraic_residual = v.get("max_algebraic_residual").and_then(Value::as_f64);
        let rows_passed = !rows.is_empty() && rows.iter().all(|r| r["converged"] == json!(true));
        let top_passed = v.get("converged").and_then(Value::as_bool).unwrap_or(true);
        Some(Check {
            factor: v["factor"].as_f64()?,
            tolerance,
            matching,
            // every row passed, and any 3-grid GCI pass/regime check agreed
            passed: rows_passed && top_passed,
            worst_relative: v["worst_relative"].as_f64().unwrap_or(f64::NAN),
            worst_quantity: v["worst_quantity"].as_str().unwrap_or("").to_string(),
            quantities: v["quantities"].as_u64().map(|n| n as usize).unwrap_or(rows.len()),
            cells_coarse: v["cells"]["coarse"].as_u64(),
            cells_refined: v["cells"]["refined"].as_u64(),
            seconds: v["seconds"].as_f64().unwrap_or(0.0),
            grids,
            cells_medium: v["cells"]["medium"].as_u64(),
            r_eff,
            observed_order,
            gci_21,
            gci_regime,
            max_algebraic_residual,
        })
    }

    pub fn to_json(&self) -> Value {
        let mut cells = json!({ "coarse": self.cells_coarse, "refined": self.cells_refined });
        if let Some(m) = self.cells_medium {
            cells["medium"] = json!(m);
        }
        let mut v = json!({
            "factor": self.factor, "tolerance": self.tolerance, "observable": self.observable(),
            "matching": self.matching, "passed": self.passed, "grids": self.grids,
            "worst_relative": self.worst_relative, "worst_quantity": self.worst_quantity,
            "cells": cells, "seconds": self.seconds,
        });
        if let Some(r) = self.r_eff {
            v["r_eff"] = json!(r);
        }
        if let Some(p) = self.observed_order {
            v["observed_order"] = json!(p);
        }
        if let Some(g) = self.gci_21 {
            v["gci_21"] = json!(g);
        }
        if let Some(reg) = &self.gci_regime {
            v["gci_regime"] = json!(reg);
        }
        v
    }

    pub fn to_record(&self) -> nyx_store::runs::Convergence {
        let tier = if self.gci_regime.as_deref() == Some("non_convergent") {
            ValidationTier::NonConvergent
        } else if self.passed
            && self.r_eff.is_none_or(|r| r >= 1.3)
            && self.max_algebraic_residual.is_none_or(|a| a <= 0.01 * self.tolerance)
        {
            if self.grids >= 3 && self.gci_21.is_some() {
                ValidationTier::GciVerified
            } else {
                ValidationTier::TwoGridChecked
            }
        } else {
            ValidationTier::Exploratory
        };
        nyx_store::runs::Convergence {
            factor: self.factor,
            tolerance: self.tolerance,
            observable: self.observable(),
            passed: self.passed,
            worst_relative: self.worst_relative,
            grids: Some(self.grids),
            r_eff: self.r_eff,
            gci_fine: self.gci_21,
            gci_p_obs: self.observed_order,
            gci_regime: self.gci_regime.clone(),
            tier: Some(tier.as_str().to_string()),
        }
    }
}

/// The three states, in the words every front-end shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Validation {
    /// not checked (yet): the numbers are what one solve gave
    Exploratory,
    /// a convergence check with a tolerance passed on this problem (or flux bounds certified it)
    Checked,
    /// not checked, and the check is expected not to fit this machine
    NotAffordable,
}

impl Validation {
    pub fn as_str(self) -> &'static str {
        match self {
            Validation::Exploratory => "exploratory",
            Validation::Checked => "checked",
            Validation::NotAffordable => "check_not_affordable",
        }
    }
}

/// Granular validation tier ordered strictly from strongest certainty to non-convergent.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ValidationTier {
    BoundCertified,
    GciVerified,
    TwoGridChecked,
    ResidualChecked,
    Exploratory,
    NotAffordable,
    NonConvergent,
}

impl ValidationTier {
    pub fn as_str(self) -> &'static str {
        match self {
            ValidationTier::BoundCertified => "bound_certified",
            ValidationTier::GciVerified => "gci_verified",
            ValidationTier::TwoGridChecked => "two_grid_checked",
            ValidationTier::ResidualChecked => "residual_checked",
            ValidationTier::Exploratory => "exploratory",
            ValidationTier::NotAffordable => "not_affordable",
            ValidationTier::NonConvergent => "non_convergent",
        }
    }

    /// Ordering of certainty: `"bound_certified"` / `"gci_verified"` < `"two_grid_checked"` <
    /// `"residual_checked"` < `"exploratory"` < `"not_affordable"` < `"non_convergent"`.
    pub fn weakness(self) -> u8 {
        match self {
            ValidationTier::BoundCertified | ValidationTier::GciVerified => 0,
            ValidationTier::TwoGridChecked => 1,
            ValidationTier::ResidualChecked => 2,
            ValidationTier::Exploratory => 3,
            ValidationTier::NotAffordable => 4,
            ValidationTier::NonConvergent => 5,
        }
    }
}

/// The state of one result with its validation tier, reason, and the check it rests on or the last one that did not pass.
#[derive(Clone, Debug, PartialEq)]
pub struct Validated {
    pub state: Validation,
    pub tier: ValidationTier,
    pub reason: String,
    pub check: Option<Check>,
}

impl Validated {
    pub fn to_json(&self) -> Value {
        let mut v = json!({
            "state": self.state.as_str(),
            "tier": self.tier.as_str(),
            "reason": self.reason,
        });
        if let Some(c) = &self.check {
            v["check"] = c.to_json();
        }
        v
    }
}

fn gb(b: u64) -> String {
    format!("{:.1} GB", b as f64 / 1e9)
}

/// Returns maximum relative Prager–Synge capacitance bound `max_{i,j} dC_ij / max(|C_ij|, 1e-30)` when present.
fn flux_bound_max_rel(r: &nyx_physics::SolveResult) -> Option<f64> {
    let cap = nyx_physics::CapacitanceMatrixData::of(r)?;
    let eb = cap.maxwell_error_bound_f.as_ref()?;
    let n = cap.maxwell_matrix_f.len();
    if n == 0 || eb.len() != n {
        return None;
    }
    let mut worst = 0.0_f64;
    for (c_row, eb_row) in cap.maxwell_matrix_f.iter().zip(eb) {
        if c_row.len() != n || eb_row.len() != n {
            return None;
        }
        for (&c_ij, &dc) in c_row.iter().zip(eb_row) {
            if !dc.is_finite() || !c_ij.is_finite() {
                return None;
            }
            let denom = c_ij.abs().max(1.0e-30);
            worst = worst.max(dc.abs() / denom);
        }
    }
    Some(worst)
}

impl Engine {
    fn target_tolerance(&self, physics: &str, check: Option<&Check>) -> f64 {
        check
            .map(|c| c.tolerance)
            .or(self.project.mesh.target_tol)
            .or_else(|| {
                self.project.spec(physics).and_then(|s| s.study.params.get("adapt_rel_tol").and_then(Value::as_f64))
            })
            .filter(|t| t.is_finite() && *t > 0.0)
            .unwrap_or(1.0e-2)
    }

    /// The validated state of the latest result of `physics`. A physics without a result is `exploratory` with that as the reason.
    pub fn validated(&self, physics: &str) -> Validated {
        let none = |reason: &str| Validated {
            state: Validation::Exploratory,
            tier: ValidationTier::Exploratory,
            reason: reason.into(),
            check: None,
        };
        let Some(stored) = self.results.get(physics) else { return none("there is no result") };
        if self.is_stale(physics) {
            // the check, if any, was of another problem: nothing is claimed for a result the project no longer describes
            return none("the result is stale: the project has changed since it was solved, so no check applies");
        }
        let check_opt = self.checks.get(&stored.key);
        let eps_alg = stored.result.diagnostics.as_ref().map(|d| d.max_algebraic_residual).unwrap_or(0.0);
        let target_tol = self.target_tolerance(physics, check_opt);

        // Prager–Synge flux-bound certification: elevates electrostatic solves even when 3-grid refinement is unaffordable.
        if physics == "dfs" && stored.result.study == "electrostatic" {
            if let Some(max_dc_rel) = flux_bound_max_rel(&stored.result) {
                if max_dc_rel <= target_tol && eps_alg <= 0.01 * target_tol {
                    return Validated {
                        state: Validation::Checked,
                        tier: ValidationTier::BoundCertified,
                        reason: format!(
                            "Prager-Synge flux bounds within target_tol {target_tol:.2e} (max dC/C = {max_dc_rel:.3e}); algebraic floor {eps_alg:.2e} <= {:.2e}",
                            0.01 * target_tol
                        ),
                        check: check_opt.cloned(),
                    };
                }
            }
        }

        if let Some(c) = check_opt {
            let check_alg = c.max_algebraic_residual.unwrap_or(eps_alg).max(eps_alg);
            let r_ok = c.r_eff.is_none_or(|r| r >= 1.3);
            let alg_ok = check_alg <= 0.01 * c.tolerance;
            if c.gci_regime.as_deref() == Some("non_convergent") {
                return Validated {
                    state: Validation::Exploratory,
                    tier: ValidationTier::NonConvergent,
                    reason: format!(
                        "3-grid Roache GCI detected non-convergent behaviour (|R| >= 1): element size times {}, worst relative difference {:.2e} ({})",
                        c.factor, c.worst_relative, c.worst_quantity
                    ),
                    check: Some(c.clone()),
                };
            }
            if c.passed {
                if !r_ok {
                    return Validated {
                        state: Validation::Exploratory,
                        tier: ValidationTier::Exploratory,
                        reason: format!(
                            "mesh-convergence check passed tolerance {:.1e} but effective refinement ratio r_eff = {:.3} < 1.3 is too close to 1 to certify",
                            c.tolerance,
                            c.r_eff.unwrap_or(0.0)
                        ),
                        check: Some(c.clone()),
                    };
                }
                if !alg_ok {
                    return Validated {
                        state: Validation::Exploratory,
                        tier: ValidationTier::Exploratory,
                        reason: format!(
                            "mesh-convergence check passed tolerance {:.1e} but algebraic residual {check_alg:.2e} exceeds 0.01 * target_tol ({:.2e})",
                            c.tolerance,
                            0.01 * c.tolerance
                        ),
                        check: Some(c.clone()),
                    };
                }
                let tier = if c.grids >= 3 && c.gci_21.is_some() {
                    ValidationTier::GciVerified
                } else {
                    ValidationTier::TwoGridChecked
                };
                return Validated {
                    state: Validation::Checked,
                    tier,
                    reason: if tier == ValidationTier::GciVerified {
                        format!(
                            "3-grid Roache GCI check passed: element size times {}, GCI_21 = {:.2e} (p_obs = {}), worst relative difference {:.2e} ({}) within {:.1e}",
                            c.factor,
                            c.gci_21.unwrap_or(c.worst_relative),
                            c.observed_order.map_or_else(|| "n/a".into(), |p| format!("{p:.2}")),
                            c.worst_relative,
                            c.worst_quantity,
                            c.tolerance
                        )
                    } else {
                        format!(
                            "a mesh-convergence check passed: element size times {}, worst relative difference {:.2e} ({}) within {:.1e}",
                            c.factor, c.worst_relative, c.worst_quantity, c.tolerance
                        )
                    },
                    check: Some(c.clone()),
                };
            } else {
                return Validated {
                    state: Validation::Exploratory,
                    tier: ValidationTier::Exploratory,
                    reason: format!(
                        "the last mesh-convergence check did not pass: element size times {}, worst relative difference {:.2e} ({}) above {:.1e}",
                        c.factor, c.worst_relative, c.worst_quantity, c.tolerance
                    ),
                    check: Some(c.clone()),
                };
            }
        }
        // not checked: can it be? The check solves on a mesh finer by the planned factor; its size follows from the cells and the dimension
        let (cells, dim) = (stored.model.mesh().n_cells() as u64, stored.model.mesh().dim() as u8);
        let refined = (cells as f64 / PLANNED_FACTOR.powi(dim as i32)).ceil() as u64;
        let chain = crate::estimate::chain_key_for_project(&self.project, &[physics.to_string()]);
        if let (Some(est), Some((total, _))) = (self.history.estimate(&chain, dim, refined), (self.memory_probe)()) {
            if est.bytes_with_margin as f64 > total as f64 * crate::estimate::REFUSE_ABOVE {
                return Validated {
                    state: Validation::NotAffordable,
                    tier: ValidationTier::NotAffordable,
                    reason: format!(
                        "no convergence check yet, and one on a mesh with the element size times {PLANNED_FACTOR} (about {refined} cells) is expected to need about {} (from {} earlier run{}{}); this machine holds {}",
                        gb(est.bytes_with_margin),
                        est.runs,
                        if est.runs == 1 { "" } else { "s" },
                        if est.extrapolated { ", extrapolated" } else { "" },
                        gb(total)
                    ),
                    check: None,
                };
            }
        }
        if let Some(budget) = &stored.result.error_budget {
            if let Some(sp) = &budget.spatial_discretization {
                if (sp.method.contains("indicator") || sp.method.contains("residual"))
                    && sp.rel_bound <= target_tol
                    && eps_alg <= 0.01 * target_tol
                {
                    return Validated {
                        state: Validation::Exploratory,
                        tier: ValidationTier::ResidualChecked,
                        reason: format!(
                            "no mesh-convergence check yet; a posteriori residual indicator {:.2e} within target_tol {:.2e}",
                            sp.rel_bound, target_tol
                        ),
                        check: None,
                    };
                }
            }
        }
        none("no mesh-convergence check with a tolerance has been made on this problem")
    }

    /// `validated` as the JSON every result answer carries.
    pub fn validated_json(&self, physics: &str) -> Value {
        self.validated(physics).to_json()
    }

    /// The weakest state of several results (a reduction over more than one physics), with one row per physics.
    pub(crate) fn validated_of(&self, physics: &[String]) -> Value {
        match physics {
            [one] => self.validated_json(one),
            many => {
                let all: Vec<(&String, Validated)> = many.iter().map(|p| (p, self.validated(p))).collect();
                let Some((p, w)) = all.iter().max_by_key(|(_, v)| v.tier.weakness()) else { return json!(null) };
                let by: serde_json::Map<String, Value> = all
                    .iter()
                    .map(|(p, v)| ((*p).clone(), json!({ "state": v.state.as_str(), "tier": v.tier.as_str() })))
                    .collect();
                json!({
                    "state": w.state.as_str(),
                    "tier": w.tier.as_str(),
                    "reason": format!("the weakest of the results used: {p}: {}", w.reason),
                    "by_physics": by,
                })
            }
        }
    }

    /// Keeps what a finished convergence check found, under the problem key it was started for, and writes it into the run that solved that
    /// problem when the engine records runs. Nothing is kept for a check made without a tolerance (it passed nothing and failed nothing).
    /// Returns the check and the run it was written to.
    pub(crate) fn note_check(&mut self, physics: &str, key: &str, answer: &Value) -> Option<(Check, Option<String>)> {
        let check = Check::from_answer(answer)?;
        self.checks.insert(key.to_string(), check.clone());
        let run = self.results.get(physics).filter(|s| s.key == key).and_then(|s| s.run.clone());
        let written = match (&mut self.link, run) {
            (Some(link), Some(rid)) => {
                let patch = nyx_store::runs::RunPatch { convergence: Some(check.to_record()), ..Default::default() };
                link.dir.patch_run(&link.experiment, &rid, patch).ok().map(|_| rid)
            }
            _ => None,
        };
        Some((check, written))
    }

    /// The answer of a finished check with the state it leaves the result in.
    pub(crate) fn checked_answer(&mut self, physics: &str, key: &str, mut v: Value) -> Result<Value, Error> {
        let noted = self.note_check(physics, key, &v);
        v["validated"] = self.validated_json(physics);
        v["recorded"] = json!(noted.is_some());
        if let Some((_, Some(rid))) = noted {
            v["run"] = json!(rid);
        }
        Ok(v)
    }
}
