//! nyx-analysis: analyses of results. An analysis is a function of a `SolveResult` (and the mesh it was solved on); results do not carry an
//! `analysis` of their own (restructure phase 5, ADR 56).
//!
//! * `polarization`: Stokes parameters, Poincare sphere, mode classes and the drawings of them (was `femcore::polarization`)
//! * `modes`: the polarization of a set of modes, the table, the CSV and the degenerate groups (was `femcore::analysis`)
//! * `samples`: the modes of a `SolveResult` as weighted samples of the field (what a mode study of `ems` states in its result)
//!
//! `polarization(mesh, result)` is the typed entry point (the window draws from it); `PolarizationAnalysis` is the same through the
//! `Analysis` trait, which returns plain JSON, so that the engine, the window and the MCP server can list and run analyses uniformly. Group
//! index is a later analysis.

pub mod modes;
pub mod polarization;
pub mod samples;

pub use modes::{analyze_modes, AnalysisOptions, GroupReport, ModeAnalysis, ModeInput, ModeRow, ResolveAxis, Resolved};
pub use samples::{polarization_input, PolarizationInput};

use nyx_physics::{Error, Mesh, SolveResult};
use serde_json::{json, Value};

/// One analysis: a name, a test of whether it applies to a result, and the run, which returns plain JSON data.
pub trait Analysis: Send + Sync {
    /// Short identifier (`polarization`).
    fn id(&self) -> &'static str;
    fn name(&self) -> &'static str;
    /// True when the result has what the analysis needs (a mode study with an electric field, ...).
    fn applies_to(&self, result: &SolveResult) -> bool;
    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error>;
}

/// The polarization of the modes of `result`, solved on `mesh`.
pub fn polarization(mesh: &Mesh, result: &SolveResult) -> Result<ModeAnalysis, Error> {
    let input = polarization_input(mesh, result)?;
    Ok(analyze_modes(&input.modes, &input.options))
}

/// Whether `polarization` can analyse `result`.
pub fn polarization_applies(result: &SolveResult) -> bool {
    samples::applies(result)
}

/// The polarization analysis behind the `Analysis` trait.
pub struct PolarizationAnalysis;

impl Analysis for PolarizationAnalysis {
    fn id(&self) -> &'static str {
        "polarization"
    }

    fn name(&self) -> &'static str {
        "Polarization of the modes"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        samples::applies(result)
    }

    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        let a = polarization(mesh, result)?;
        let report = |m: &polarization::ModeReport| {
            json!({
                "class": m.class.label(), "s": m.sphere, "dop": m.dop, "chi_deg": m.chi_deg, "psi_deg": m.psi_deg,
                "te_fraction": m.te_fraction, "long_fraction": m.long_fraction
            })
        };
        let rows: Vec<Value> = a
            .rows
            .iter()
            .enumerate()
            .map(|(i, r)| {
                json!({
                    "mode": i + 1, "label": r.label, "value": [r.neff.0, r.neff.1],
                    "report": r.report.as_ref().map(report), "group": r.group.map(|g| g + 1)
                })
            })
            .collect();
        let groups: Vec<Value> = a
            .groups
            .iter()
            .map(|g| {
                json!({
                    "members": g.members.iter().map(|i| i + 1).collect::<Vec<_>>(),
                    "invariant": g.invariant.as_ref().map(report), "overlap": g.overlap, "note": g.note
                })
            })
            .collect();
        Ok(json!({ "value_label": a.value_label, "rows": rows, "groups": groups, "table": a.table() }))
    }
}
