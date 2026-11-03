//! nyx-analysis: analyses of results. An analysis is a function of a `SolveResult` (and the mesh it was solved on); results no longer
//! carry an `analysis` of their own. Planned content (restructure phase 5): polarization (Stokes parameters, Poincare sphere, mode
//! classes, degenerate groups) moved from `femcore::{analysis, polarization, optical}`, with axis names from the wafer frame; group
//! index later.
//!
//! SKELETON (restructure foundation): only the shape of an analysis is fixed here, so the engine, the window and the MCP server can
//! list and run analyses uniformly. The phase 5 owner may reshape the trait; the contract `analysis_api` is then versioned.

use nyx_physics::{Error, Mesh, SolveResult};
use serde_json::Value;

/// One analysis: a name, a test of whether it applies to a result, and the run, which returns plain JSON data.
pub trait Analysis: Send + Sync {
    /// Short identifier (`polarization`).
    fn id(&self) -> &'static str;
    fn name(&self) -> &'static str;
    /// True when the result has what the analysis needs (a mode study with an electric field, ...).
    fn applies_to(&self, result: &SolveResult) -> bool;
    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error>;
}
