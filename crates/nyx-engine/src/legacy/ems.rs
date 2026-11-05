//! The window's EM-wave problems as `nyx_ems::Problem` (ADR 41, 47). The old problem names no study kind for EM waves: a 2D mesh is a
//! waveguide cross-section (`Study::Modes`), a 3D mesh is a cavity (`Study::Eigenfrequency`), as the old dispatch of `femcore::solve` did.
//! The polarization analysis the window shows is not part of the result (phase 5, ADR 56): `solve_analysed` runs `nyx_analysis::polarization`
//! on the result of `nyx_ems::run` and hands it out next to the old-shape result.

use super::{error_text, result_back};
use femcore::mesh::SimplexMesh;
use femcore::problem::{self as old, SolveState};
use femcore::progress::Progress;
use std::sync::Arc;

#[cfg(feature = "ems")]
pub fn problem_from_old(spec: &old::ProblemSpec, dim: usize) -> Result<nyx_ems::Problem, String> {
    use nyx_ems::{Problem, Study};
    let study = match (dim, &spec.study) {
        (3, old::Study::Eigenfrequency { n_modes, f_search }) => {
            Study::Cavity { n_modes: *n_modes, f_search: *f_search }
        }
        (3, _) => return Err("EM waves on a 3D mesh use the Eigenfrequency study".into()),
        (_, old::Study::Modes { freq, n_target, evanescent, n_modes }) => {
            Study::Modes { freq: *freq, n_target: *n_target, evanescent: *evanescent, n_modes: *n_modes }
        }
        _ => return Err("mode analysis needs a Modes study".into()),
    };
    Ok(Problem { order: spec.order, study, materials: spec.materials.clone(), periodic: spec.periodic.clone() })
}

pub fn solve(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    solve_analysed(mesh, spec, state, progress).map(|(r, _)| r)
}

/// `solve` with the polarization analysis of the modes (none for a study whose result states no frame). The analysis table is also appended to
/// the log of the old-shape result, and `Solution::te_fraction` is set from it, for the consumers of the old result.
#[cfg(feature = "ems")]
pub fn solve_analysed(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    _state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<(old::SolveResult, Option<nyx_analysis::ModeAnalysis>), String> {
    let problem = problem_from_old(spec, mesh.dim)?;
    let exchange = nyx_physics::Mesh::from_solver_mesh(mesh)?;
    let result = nyx_ems::run(&exchange, &problem, progress).map_err(error_text)?;
    let analysis = if nyx_analysis::polarization_applies(&result) {
        Some(nyx_analysis::polarization(&exchange, &result).map_err(error_text)?)
    } else {
        None
    };
    let mut back = result_back(&result);
    if let Some(a) = &analysis {
        for (sol, row) in back.solutions.iter_mut().zip(&a.rows) {
            sol.te_fraction = row.report.as_ref().map(|m| m.te_fraction);
        }
        back.log.push(
            "polarization of the modes (Stokes parameters; s1 > 0 is x, the lateral direction, s1 < 0 is z, the surface normal):".to_string(),
        );
        back.log.extend(a.table().lines().filter(|l| !l.is_empty()).map(|l| l.to_string()));
    }
    if result.study == "cavity" {
        // the window never had an effective index for a cavity: its eigenvalue is a frequency
        back.solutions.iter_mut().for_each(|s| s.neff = None);
    }
    Ok((back, analysis))
}

#[cfg(not(feature = "ems"))]
pub fn solve_analysed(
    _mesh: &SimplexMesh,
    _spec: &old::ProblemSpec,
    _state: &mut SolveState,
    _progress: &Arc<Progress>,
) -> Result<(old::SolveResult, Option<nyx_analysis::ModeAnalysis>), String> {
    Err("this build does not contain the ems physics".into())
}
