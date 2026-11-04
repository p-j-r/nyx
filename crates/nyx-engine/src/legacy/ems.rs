//! The window's EM-wave problems as `nyx_ems::Problem` (ADR 41, 47). The old problem names no study kind for EM waves: a 2D mesh is a
//! waveguide cross-section (`Study::Modes`), a 3D mesh is a cavity (`Study::Eigenfrequency`), as the old dispatch of `femcore::solve` did.
//! The polarization analysis the window shows (`SolveResult::analysis`) comes back from `nyx_ems::run_analysed`.

use super::{error_text, result_back};
use femcore::mesh::SimplexMesh;
use femcore::problem::{self as old, SolveState};
use femcore::progress::Progress;
use std::sync::Arc;

#[cfg(feature = "ems")]
pub fn problem_from_old(spec: &old::ProblemSpec, dim: usize) -> Result<nyx_ems::Problem, String> {
    use nyx_ems::{Problem, Study};
    let study = match (dim, &spec.study) {
        (3, old::Study::Eigenfrequency { n_modes, f_search }) => Study::Cavity { n_modes: *n_modes, f_search: *f_search },
        (3, _) => return Err("EM waves on a 3D mesh use the Eigenfrequency study".into()),
        (_, old::Study::Modes { freq, n_target, evanescent, n_modes }) => {
            Study::Modes { freq: *freq, n_target: *n_target, evanescent: *evanescent, n_modes: *n_modes }
        }
        _ => return Err("mode analysis needs a Modes study".into()),
    };
    Ok(Problem { order: spec.order, study, materials: spec.materials.clone(), periodic: spec.periodic.clone() })
}

#[cfg(feature = "ems")]
pub fn solve(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    _state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    let problem = problem_from_old(spec, mesh.dim)?;
    let exchange = nyx_physics::Mesh::from_solver_mesh(mesh)?;
    let solved = nyx_ems::run_analysed(&exchange, &problem, progress).map_err(error_text)?;
    let mut back = result_back(&solved.result);
    back.analysis = solved.analysis;
    Ok(back)
}

#[cfg(not(feature = "ems"))]
pub fn solve(
    _mesh: &SimplexMesh,
    _spec: &old::ProblemSpec,
    _state: &mut SolveState,
    _progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    Err("this build does not contain the ems physics".into())
}
