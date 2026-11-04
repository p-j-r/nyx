//! Shared helpers of the ems tests: the mesh fixtures of `femcore::testkit`, a solve helper over the typed problem, the effective index of a
//! solution, and the conversion of an old-shape `ProblemSpec` (what the use cases build) to the typed problem.
#![allow(dead_code, unused_imports)]
use femcore::progress::Progress;
pub use femcore::testkit::*;
use nyx_ems::{Problem, Solved, Study};
use nyx_physics::{Discretisation, Mesh, Solution, SolveResult};

/// An EM-wave problem with the defaults of the old `ProblemSpec::new(Physics::Emw)`: order 2, a `Modes` study at 10 GHz.
pub fn emw() -> Problem {
    Problem::new(Study::Modes { freq: 10e9, n_target: 1.0, evanescent: false, n_modes: 4 }, 2)
}

/// Solves a typed problem on a fixture mesh (what `Physics::solve` does after `decode`); the stand-in of the old `solve::solve`.
pub fn solve(mesh: &Discretisation, p: &Problem) -> Result<SolveResult, String> {
    Ok(try_solve_analysed(mesh, p)?.result)
}

/// `solve` with the polarization analysis next to the result.
pub fn try_solve_analysed(mesh: &Discretisation, p: &Problem) -> Result<Solved, String> {
    let exchange = Mesh::from_solver_mesh(mesh)?;
    nyx_ems::run_analysed(&exchange, p, &Progress::new()).map_err(|e| e.to_string())
}

/// The effective index (re, im) of a mode solution.
pub fn neff(s: &Solution) -> (f64, f64) {
    let e = s.eigenvalue.expect("a mode with an effective index");
    (e[0], e[1])
}

/// The typed problem of an old-shape EM-wave spec built by a use case (same mesh numbers, no names).
pub fn from_old(spec: &femcore::problem::ProblemSpec) -> Problem {
    use femcore::problem::Study as Old;
    let study = match &spec.study {
        Old::Modes { freq, n_target, evanescent, n_modes } => {
            Study::Modes { freq: *freq, n_target: *n_target, evanescent: *evanescent, n_modes: *n_modes }
        }
        Old::Eigenfrequency { n_modes, f_search } => Study::Cavity { n_modes: *n_modes, f_search: *f_search },
        other => panic!("not an EM-wave study: {other:?}"),
    };
    Problem { order: spec.order, study, materials: spec.materials.clone(), periodic: spec.periodic.clone() }
}

/// Solves an old-shape spec (what the mesher's `physics_size` reads) through the typed problem.
pub fn solve_old(mesh: &Discretisation, spec: &femcore::problem::ProblemSpec) -> Result<SolveResult, String> {
    solve(mesh, &from_old(spec))
}
