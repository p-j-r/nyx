//! The one entry point for all physics: `solve(mesh, spec, state)`.
//!
//! It checks the element order and the numbers of the problem (`validate`), dispatches to the module of the chosen physics
//! (`ems`/`emw3d`; `dfs` and `sms` are in their crates), which read the wafer-frame problem as it is and name their fields in wafer axes
//! (`section`), refuses a result that contains NaN or infinity, and appends the plain-
//! language stage log of the run to `SolveResult::log`. `solve_with` does the same for a `Progress` that a window can
//! watch and cancel; a cancelled solve returns `Err("cancelled")`.

use crate::mesh::SimplexMesh;
use crate::problem::*;
use crate::progress::{self, Progress};
use std::sync::Arc;

/// Solves `spec` on `mesh`. Heat-transfer results are kept in `state` so a later solid
/// mechanics study can use them as the thermal load.
pub fn solve(mesh: &SimplexMesh, spec: &ProblemSpec, state: &mut SolveState) -> Result<SolveResult, String> {
    // inside another solve (a study that solves sub-problems) the outer job's progress and cancel flag are kept
    let p = progress::current().unwrap_or_else(Progress::new);
    solve_with(mesh, spec, state, &p)
}

/// `solve` that reports to (and can be cancelled through) `progress`. The plain-language stage log of the run is
/// appended to `SolveResult::log`; a cancelled solve returns `Err("cancelled")`.
pub fn solve_with(mesh: &SimplexMesh, spec: &ProblemSpec, state: &mut SolveState, progress: &Arc<Progress>) -> Result<SolveResult, String> {
    let _installed = progress::install(progress);
    let first_line = progress.log().len();
    let mut result = {
        let _s = progress::stage(format!("solving {} on a {}D mesh of {} cells", spec.physics.name(), mesh.dim, mesh.n_cells()));
        solve_checked(mesh, spec, state)?
    };
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}

fn solve_checked(mesh: &SimplexMesh, spec: &ProblemSpec, state: &mut SolveState) -> Result<SolveResult, String> {
    if spec.order != 1 && spec.order != 2 {
        return Err("element order must be 1 or 2".into());
    }
    progress::check()?;
    crate::validate::check_inputs(mesh, spec)?;
    // every kernel reads the problem in wafer axes and names its fields in wafer axes (`section`)
    let result = dispatch(mesh, spec, state)?;
    crate::validate::check_result(&result)?;
    Ok(result)
}

fn dispatch(mesh: &SimplexMesh, spec: &ProblemSpec, _state: &mut SolveState) -> Result<SolveResult, String> {
    match spec.physics {
        // dfs lives in its crate (restructure phase 4): the window reaches it through `nyx_engine::legacy`
        Physics::Electrostatics | Physics::Heat => Err("electrostatics and heat transfer moved to the nyx-dfs crate; solve through nyx_engine::legacy::solve".into()),
        // sms lives in its crate too (restructure phase 4, ADR 42)
        Physics::Solid => Err("solid mechanics moved to the nyx-sms crate; solve through nyx_engine::legacy::solve".into()),
        Physics::Emw if mesh.dim == 3 => crate::ems3d::solve_emw_eigen(mesh, spec),
        Physics::Emw => crate::ems::solve_emw(mesh, spec),
    }
}
