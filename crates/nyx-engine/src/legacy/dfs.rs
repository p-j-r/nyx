//! The window's electrostatics and heat problems as `nyx_dfs::Problem` (ADR 41).

use super::{error_text, result_back};
use femcore::mesh::SimplexMesh;
use femcore::problem::{self as old, SolveState};
use femcore::progress::Progress;
use std::sync::Arc;

#[cfg(feature = "dfs")]
pub fn problem_from_old(spec: &old::ProblemSpec) -> Result<nyx_dfs::Problem, String> {
    use nyx_dfs::{Bc, BcKind, Study};
    let electro = spec.physics == old::Physics::Electrostatics;
    let study = match (&spec.study, electro) {
        (_, true) => Study::Electrostatic,
        (old::Study::Stationary, false) => Study::HeatStationary,
        (old::Study::Transient { t_end, n_out, t0 }, false) => {
            Study::HeatTransient { t_end: *t_end, n_out: *n_out, t0: t0.clone() }
        }
        (_, false) => return Err("heat transfer supports Stationary and Transient studies".into()),
    };
    let mut bcs = vec![];
    for bc in &spec.bcs {
        let kind = match &bc.kind {
            old::BcKind::Potential(v) => BcKind::Potential(*v),
            old::BcKind::Ground => BcKind::Ground,
            old::BcKind::SurfaceCharge(s) => BcKind::SurfaceCharge(*s),
            old::BcKind::Temperature(t) => BcKind::Temperature(*t),
            old::BcKind::Flux(q) => BcKind::Flux(*q),
            old::BcKind::Convection { h, text } => BcKind::Convection { h: *h, text: *text },
            other => {
                return Err(format!(
                    "'{}' is not {} boundary condition",
                    other.name(),
                    if electro { "an electrostatics" } else { "a heat-transfer" }
                ))
            }
        };
        bcs.push(Bc { labels: bc.labels.clone(), kind });
    }
    Ok(nyx_dfs::Problem { order: spec.order, study, materials: spec.materials.clone(), bcs, solver: spec.solver })
}

#[cfg(feature = "dfs")]
pub fn solve(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    let problem = problem_from_old(spec)?;
    let exchange = nyx_physics::Mesh::from_solver_mesh(mesh)?;
    let r = nyx_dfs::run(&exchange, &problem, progress).map_err(error_text)?;
    if spec.physics == old::Physics::Heat {
        // a later solid study can take the last temperature as its thermal load
        state.heat_t = r.solutions.last().and_then(|s| s.field("T")).map(|f| f.values.clone());
    }
    Ok(result_back(&r))
}

#[cfg(not(feature = "dfs"))]
pub fn solve(
    _mesh: &SimplexMesh,
    _spec: &old::ProblemSpec,
    _state: &mut SolveState,
    _progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    Err("this build does not contain the dfs physics".into())
}
