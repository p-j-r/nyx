//! The bridge from the pre-restructure window problem (`femcore::problem::ProblemSpec`, `SolveResult`) to the physics crates, for the
//! window and the MCP session until restructure phase 7 rewires them to engine commands (ADR 41). One file per physics converts the old
//! problem to that physics' typed problem (same mesh numbers, no names), runs it, and converts the result back to the old shape.
//! All three physics (`dfs`, `sms`, `ems`) run in their crates; `femcore::solve` is gone.
//!
//! Delete this module and `femcore::problem` together in phase 7.

use femcore::mesh::SimplexMesh;
use femcore::problem::{self as old, Physics as OldPhysics, SolveState};
use femcore::progress::{self, Progress};
use nyx_physics::{Error, Scalar, Solution, SolveResult};
use std::sync::Arc;

pub mod dfs;
pub mod ems;
pub mod sms;

/// The one entry point of the window: the kernels of every physics, which all live in their crates now.
pub fn solve(mesh: &SimplexMesh, spec: &old::ProblemSpec, state: &mut SolveState) -> Result<old::SolveResult, String> {
    // inside another solve (a study that solves sub-problems) the outer job's progress and cancel flag are kept
    let p = progress::current().unwrap_or_else(Progress::new);
    solve_with(mesh, spec, state, &p)
}

/// `solve` that reports to (and can be cancelled through) `progress`; a cancelled solve returns `Err("cancelled")`.
pub fn solve_with(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    match spec.physics {
        OldPhysics::Electrostatics | OldPhysics::Heat => dfs::solve(mesh, spec, state, progress),
        OldPhysics::Solid => sms::solve(mesh, spec, state, progress),
        OldPhysics::Emw => ems::solve(mesh, spec, state, progress),
    }
}

/// `solve_with` that also hands out the polarization analysis of the modes of an EM-wave study (`None` for every other study). The result
/// does not carry the analysis (phase 5); the window and the MCP session draw it from here.
pub fn solve_with_analysis(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<(old::SolveResult, Option<nyx_analysis::ModeAnalysis>), String> {
    match spec.physics {
        OldPhysics::Emw => ems::solve_analysed(mesh, spec, state, progress),
        _ => solve_with(mesh, spec, state, progress).map(|r| (r, None)),
    }
}

/// An `Error` of a physics as the `String` error of the old entry points (`"cancelled"` is the cancel signal there).
pub(crate) fn error_text(e: Error) -> String {
    e.to_string()
}

fn scalar_back(s: &Scalar) -> old::Scalar {
    match s {
        Scalar::Num(x) => old::Scalar::Num(*x),
        Scalar::Text(t) => old::Scalar::Text(t.clone()),
    }
}

/// A result of a physics crate in the shape the window reads: flat field names (`|E|`, `Ex`, ...), per-cell fields as `CellFields`.
/// `disp` and `neff` are filled by the physics whose results have them.
pub(crate) fn result_back(r: &SolveResult) -> old::SolveResult {
    old::SolveResult {
        order: r.order as usize,
        dim: r.dim as usize,
        solutions: r.solutions.iter().map(|s| solution_back(s, r.order as usize)).collect(),
        scalars: r.scalars.iter().map(|(n, s)| (n.clone(), scalar_back(s))).collect(),
        log: r.log.clone(),
        seconds: r.seconds,
        dofs: r.dofs,
    }
}

fn solution_back(s: &Solution, order: usize) -> old::Solution {
    old::Solution {
        label: s.label.clone(),
        fields: s.flat_fields().into_iter().map(|(name, unit, values)| old::Field { name, unit, values }).collect(),
        disp: None,
        value: s.value,
        neff: s.eigenvalue.map(|e| (e[0], e[1])),
        te_fraction: None,
        cell_fields: if s.cell_fields.is_empty() {
            None
        } else {
            Some(nyx_physics::cell_fields_back(order, &s.cell_fields))
        },
    }
}
