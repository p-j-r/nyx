//! The window's solid mechanics problem as `nyx_sms::Problem` (ADR 41, 42): the old wafer-frame `ProblemSpec` with the same mesh numbers
//! becomes the typed problem, the kernel runs, and the result goes back in the old shape. `disp`, the displacement per dof that the window
//! draws the deformed shape from, is read from the result's fields: the vector field `u` of a real study, the `Re u{x,z,y}` components
//! of a complex-frequency mode.

use super::{error_text, result_back};
use femcore::mesh::SimplexMesh;
use femcore::problem::{self as old, SolveState};
use femcore::progress::Progress;
use std::sync::Arc;

#[cfg(feature = "sms")]
pub fn problem_from_old(spec: &old::ProblemSpec, state: &SolveState) -> Result<nyx_sms::Problem, String> {
    use nyx_sms::{Bc, BcKind, Plane, Study, ThermalExpansion, ThermalSource};
    let study = match &spec.study {
        old::Study::Static => Study::Static,
        old::Study::Eigenfrequency { n_modes, f_search } => {
            Study::Eigenfrequency { n_modes: *n_modes, f_search: *f_search }
        }
        _ => return Err("solid mechanics supports Static and Eigenfrequency studies".into()),
    };
    let mut bcs = vec![];
    for bc in &spec.bcs {
        let kind = match &bc.kind {
            old::BcKind::Fixed { fixed, values } => BcKind::Fixed { fixed: *fixed, values: *values },
            old::BcKind::Load(t) => BcKind::Load(*t),
            old::BcKind::LowReflecting => BcKind::LowReflecting,
            other => return Err(format!("'{}' is not a solid mechanics boundary condition", other.name())),
        };
        bcs.push(Bc { labels: bc.labels.clone(), kind });
    }
    let thermal = match &spec.thermal {
        None => None,
        Some(th) => Some(ThermalExpansion {
            source: match &th.source {
                old::ThermalSource::Uniform(t) => ThermalSource::Uniform(*t),
                old::ThermalSource::FromHeat => match &state.heat_t {
                    Some(t) => ThermalSource::Imported(t.clone()),
                    None => {
                        return Err(
                            "no heat-transfer result on this mesh and element order; solve Heat transfer first".into()
                        )
                    }
                },
            },
            t_ref: th.t_ref,
        }),
    };
    Ok(nyx_sms::Problem {
        order: spec.order,
        study,
        materials: spec.materials.clone(),
        bcs,
        plane: match spec.plane {
            old::Plane::Strain => Plane::Strain,
            old::Plane::Stress => Plane::Stress,
        },
        thickness: spec.thickness,
        thermal,
        solver: spec.solver,
        inactive: spec.inactive.clone(),
        periodic: spec.periodic.clone(),
    })
}

/// The displacement per dof (interleaved by dof, one value per mesh axis) of a solution: the vector field `u`, or the real parts of the
/// complex mode shape (`Re ux`, `Re uz`, and `Re uy` in 3D).
#[cfg(feature = "sms")]
fn disp_of(sol: &nyx_physics::Solution, dim: usize) -> Option<Vec<f64>> {
    if let Some(u) = sol.field("u") {
        return Some(u.values.clone());
    }
    let comps: Option<Vec<&Vec<f64>>> = femcore::section::component_names(dim)
        .iter()
        .map(|ax| sol.field(&format!("Re u{ax}")).map(|f| &f.values))
        .collect();
    let comps = comps?;
    let n = comps.first()?.len();
    let mut disp = vec![0.0; n * dim];
    for (k, c) in comps.iter().enumerate() {
        for i in 0..n {
            disp[i * dim + k] = c[i];
        }
    }
    Some(disp)
}

#[cfg(feature = "sms")]
pub fn solve(
    mesh: &SimplexMesh,
    spec: &old::ProblemSpec,
    state: &mut SolveState,
    progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    let problem = problem_from_old(spec, state)?;
    let exchange = nyx_physics::Mesh::from_solver_mesh(mesh)?;
    let r = nyx_sms::run(&exchange, &problem, progress).map_err(error_text)?;
    let mut back = result_back(&r);
    for (b, s) in back.solutions.iter_mut().zip(&r.solutions) {
        b.disp = disp_of(s, mesh.dim);
    }
    Ok(back)
}

#[cfg(not(feature = "sms"))]
pub fn solve(
    _mesh: &SimplexMesh,
    _spec: &old::ProblemSpec,
    _state: &mut SolveState,
    _progress: &Arc<Progress>,
) -> Result<old::SolveResult, String> {
    Err("this build does not contain the sms physics".into())
}
