//! Shared helpers of the sms tests: the mesh fixtures of `femcore::testkit` and a solve helper over the typed problem (the entry the
//! engine uses: `nyx_sms::run`, which checks the numbers first).
#![allow(dead_code, unused_imports)]
use femcore::progress::Progress;
pub use femcore::testkit::*;
use nyx_physics::{Discretisation, Mesh, SolveResult};
use nyx_sms::Problem;

/// Solves a typed problem on a fixture mesh (what `Physics::solve` does after `decode`).
pub fn solve(mesh: &Discretisation, p: &Problem) -> SolveResult {
    try_solve(mesh, p).unwrap()
}

pub fn try_solve(mesh: &Discretisation, p: &Problem) -> Result<SolveResult, String> {
    let exchange = Mesh::from_solver_mesh(mesh)?;
    nyx_sms::run(&exchange, p, &Progress::new()).map_err(|e| e.to_string())
}
