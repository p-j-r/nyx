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
use serde_json::{json, Value};

/// The element-size ratio a check is planned with when none was made (the step the command's own help recommends).
pub const PLANNED_FACTOR: f64 = 0.6;

/// What a mesh-convergence check with a tolerance found, for one problem key.
#[derive(Clone, Debug, PartialEq)]
pub struct Check {
    pub factor: f64,
    pub tolerance: f64,
    pub passed: bool,
    pub worst_relative: f64,
    pub worst_quantity: String,
    pub quantities: usize,
    pub cells_coarse: Option<u64>,
    pub cells_refined: Option<u64>,
    pub seconds: f64,
}

impl Check {
    /// What was compared and how, in words (the `observable` of the run record).
    pub fn observable(&self) -> String {
        format!(
            "{} quantit{} (eigenvalues, values, numeric scalars), solutions matched by index",
            self.quantities,
            if self.quantities == 1 { "y" } else { "ies" }
        )
    }

    /// The part of a convergence answer that is a verdict; `None` for a check that had no tolerance (it says nothing is passed or failed).
    pub(crate) fn from_answer(v: &Value) -> Option<Check> {
        let tolerance = v["tolerance"].as_f64()?;
        let rows = v["rows"].as_array()?;
        Some(Check {
            factor: v["factor"].as_f64()?,
            tolerance,
            // every row passed, which is what `converged` says; recomputed so a hand-made answer cannot disagree with its rows
            passed: !rows.is_empty() && rows.iter().all(|r| r["converged"] == json!(true)),
            worst_relative: v["worst_relative"].as_f64().unwrap_or(f64::NAN),
            worst_quantity: v["worst_quantity"].as_str().unwrap_or("").to_string(),
            quantities: v["quantities"].as_u64().map(|n| n as usize).unwrap_or(rows.len()),
            cells_coarse: v["cells"]["coarse"].as_u64(),
            cells_refined: v["cells"]["refined"].as_u64(),
            seconds: v["seconds"].as_f64().unwrap_or(0.0),
        })
    }

    pub fn to_json(&self) -> Value {
        json!({
            "factor": self.factor, "tolerance": self.tolerance, "observable": self.observable(), "matching": "index", "passed": self.passed,
            "worst_relative": self.worst_relative, "worst_quantity": self.worst_quantity,
            "cells": { "coarse": self.cells_coarse, "refined": self.cells_refined }, "seconds": self.seconds,
        })
    }

    pub fn to_record(&self) -> nyx_store::runs::Convergence {
        nyx_store::runs::Convergence {
            factor: self.factor,
            tolerance: self.tolerance,
            observable: self.observable(),
            passed: self.passed,
            worst_relative: self.worst_relative,
        }
    }
}

/// The three states, in the words every front-end shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Validation {
    /// not checked (yet): the numbers are what one solve gave
    Exploratory,
    /// a convergence check with a tolerance passed on this problem
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

    /// How much a reader may lean on it, for the weakest of several results: checked, then exploratory, then not affordable.
    fn weakness(self) -> u8 {
        match self {
            Validation::Checked => 0,
            Validation::Exploratory => 1,
            Validation::NotAffordable => 2,
        }
    }
}

/// The state of one result with the reason, and the check it rests on or the last one that did not pass.
#[derive(Clone, Debug, PartialEq)]
pub struct Validated {
    pub state: Validation,
    pub reason: String,
    pub check: Option<Check>,
}

impl Validated {
    pub fn to_json(&self) -> Value {
        let mut v = json!({ "state": self.state.as_str(), "reason": self.reason });
        if let Some(c) = &self.check {
            v["check"] = c.to_json();
        }
        v
    }
}

fn gb(b: u64) -> String {
    format!("{:.1} GB", b as f64 / 1e9)
}

impl Engine {
    /// The validated state of the latest result of `physics`. A physics without a result is `exploratory` with that as the reason.
    pub fn validated(&self, physics: &str) -> Validated {
        let none = |reason: &str| Validated { state: Validation::Exploratory, reason: reason.into(), check: None };
        let Some(stored) = self.results.get(physics) else { return none("there is no result") };
        if self.is_stale(physics) {
            // the check, if any, was of another problem: nothing is claimed for a result the project no longer describes
            return none("the result is stale: the project has changed since it was solved, so no check applies");
        }
        if let Some(c) = self.checks.get(&stored.key) {
            return if c.passed {
                Validated {
                    state: Validation::Checked,
                    reason: format!(
                        "a mesh-convergence check passed: element size times {}, worst relative difference {:.2e} ({}) within {:.1e}",
                        c.factor, c.worst_relative, c.worst_quantity, c.tolerance
                    ),
                    check: Some(c.clone()),
                }
            } else {
                Validated {
                    state: Validation::Exploratory,
                    reason: format!(
                        "the last mesh-convergence check did not pass: element size times {}, worst relative difference {:.2e} ({}) above {:.1e}",
                        c.factor, c.worst_relative, c.worst_quantity, c.tolerance
                    ),
                    check: Some(c.clone()),
                }
            };
        }
        // not checked: can it be? The check solves on a mesh finer by the planned factor; its size follows from the cells and the dimension
        let (cells, dim) = (stored.model.mesh().n_cells() as u64, stored.model.mesh().dim() as u8);
        let refined = (cells as f64 / PLANNED_FACTOR.powi(dim as i32)).ceil() as u64;
        let chain = crate::estimate::chain_key(&[physics.to_string()]);
        if let (Some(est), Some((total, _))) = (self.history.estimate(&chain, dim, refined), (self.memory_probe)()) {
            if est.bytes_with_margin as f64 > total as f64 * crate::estimate::REFUSE_ABOVE {
                return Validated {
                    state: Validation::NotAffordable,
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
                let Some((p, w)) = all.iter().max_by_key(|(_, v)| v.state.weakness()) else { return json!(null) };
                let by: serde_json::Map<String, Value> =
                    all.iter().map(|(p, v)| ((*p).clone(), json!(v.state.as_str()))).collect();
                json!({ "state": w.state.as_str(), "reason": format!("the weakest of the results used: {p}: {}", w.reason), "by_physics": by })
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
