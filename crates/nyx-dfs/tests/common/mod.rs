// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Shared helpers of the dfs tests: the mesh fixtures of `femcore::testkit`, a solve helper over the typed problem, and the flat field
//! names the golden files and the older tests use (`|E|`, `Ex`, `Ez`, `|q|`, ...).
#![allow(dead_code, unused_imports)]
use femcore::progress::Progress;
pub use femcore::testkit::*;
use nyx_dfs::Problem;
use nyx_physics::{Discretisation, Mesh, SolveResult};
use std::collections::BTreeSet;

/// Solves a typed problem on a fixture mesh (what `Physics::solve` does after `decode`).
pub fn solve(mesh: &Discretisation, p: &Problem) -> SolveResult {
    try_solve(mesh, p).unwrap()
}

pub fn try_solve(mesh: &Discretisation, p: &Problem) -> Result<SolveResult, String> {
    let exchange = Mesh::from_solver_mesh(mesh)?;
    nyx_dfs::run(&exchange, p, &Progress::new()).map_err(|e| e.to_string())
}

/// The flat names of the fields of the first solution, sorted.
pub fn names(r: &SolveResult) -> BTreeSet<String> {
    r.solutions[0].flat_fields().into_iter().map(|(n, _, _)| n).collect()
}

/// The mean of a field of the first solution by its flat name.
pub fn mean(r: &SolveResult, field: &str) -> f64 {
    let v = r.solutions[0].flat(field).unwrap_or_else(|| panic!("no field {field} in {:?}", names(r)));
    v.iter().sum::<f64>() / v.len() as f64
}
