//! The one entry point for all physics: `solve(mesh, spec, state)`.
//!
//! It checks the element order and the numbers of the problem (`validate`), translates a wafer-frame problem into the
//! solvers' model frame, dispatches to the module of the chosen physics (`scalar`, `solid`, `emw`/`emw3d`,
//! `brillouin`), renames the result fields back, refuses a result that contains NaN or infinity, and appends the plain-
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
    let result = if spec.materials.frame == Frame::Wafer && mesh.dim == 2 {
        // solve in the axes of the 2D model and name the results in the wafer axes
        let inner = spec.to_model_frame(2);
        let mut result = solve_model_frame(mesh, &inner, state)?;
        rename_wafer_fields(&mut result);
        result
    } else {
        solve_model_frame(mesh, spec, state)?
    };
    crate::validate::check_result(&result)?;
    Ok(result)
}

/// Result fields of a 2D model (x, z, out of plane = -y) under the names of the wafer axes (x, y, z).
pub fn rename_wafer_fields(r: &mut SolveResult) {
    // (internal name, wafer name, sign)
    const MAP: [(&str, &str, f64); 22] = [
        ("uy", "uz", 1.0),
        ("uz", "uy", -1.0),
        ("|uz|", "|uy|", 1.0),
        ("syy", "szz", 1.0),
        ("szz", "syy", 1.0),
        ("sxy", "sxz", 1.0),
        ("Ey", "Ez", 1.0),
        ("qy", "qz", 1.0),
        ("Etx", "Ex", 1.0),
        ("Ety", "Ez", 1.0),
        ("Ez", "Ey", -1.0),
        ("|Ez|", "|Ey|", 1.0),
        ("Re Etx", "Re Ex", 1.0),
        ("Re Ety", "Re Ez", 1.0),
        ("Re Ez", "Re Ey", -1.0),
        ("|Et|", "|Exz|", 1.0),
        ("Re ux", "Re ux", 1.0),
        ("Re uy", "Re uz", 1.0),
        ("Re uz", "Re uy", -1.0),
        ("|ux|", "|ux|", 1.0),
        ("|uy|", "|uz|", 1.0),
        ("ux", "ux", 1.0),
    ];
    for sol in &mut r.solutions {
        let originals: Vec<(String, Vec<f64>)> = sol.fields.iter().map(|f| (f.name.clone(), f.values.clone())).collect();
        for (f, (orig, vals)) in sol.fields.iter_mut().zip(&originals) {
            if let Some((_, to, sg)) = MAP.iter().find(|(from, _, _)| from == orig) {
                f.name = (*to).to_string();
                if *sg < 0.0 {
                    f.values = vals.iter().map(|v| -v).collect();
                }
            }
        }
    }
}

fn solve_model_frame(mesh: &SimplexMesh, spec: &ProblemSpec, state: &mut SolveState) -> Result<SolveResult, String> {
    match spec.physics {
        Physics::Electrostatics => crate::scalar::solve_electrostatics(mesh, spec),
        Physics::Heat => crate::scalar::solve_heat(mesh, spec, state),
        Physics::Solid => crate::solid::solve_solid(mesh, spec, state),
        Physics::Emw if mesh.dim == 3 => crate::emw3d::solve_emw_eigen(mesh, spec),
        Physics::Emw => crate::emw::solve_emw(mesh, spec),
        Physics::Brillouin => crate::brillouin::solve_brillouin(mesh, spec, state),
    }
}
