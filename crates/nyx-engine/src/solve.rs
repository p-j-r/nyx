//! Solving a project: the physics in an order that satisfies the one-way couplings, each checked against its descriptors and the mesh, then
//! solved with the fields it imports from the physics before it. Everything goes through `Physics`, so nothing here names a physics.

use crate::meshing::MeshedModel;
use crate::physics_for;
use crate::project::Project;
use femcore::progress::Progress;
use nyx_physics::{validate_spec, Ctx, Error, Registry, SolveResult};
use std::sync::Arc;

/// The values of the field a physics exports under `role`, taken from its last solution.
fn exported(registry: &Registry, physics: &str, result: &SolveResult, role: &str) -> Option<Vec<f64>> {
    let p = registry.get(physics)?;
    let link = p.exports().into_iter().find(|l| l.role == role)?;
    let name = link.field.split_once('.').map(|(_, n)| n).unwrap_or(&link.field);
    result.solutions.last()?.field(name).map(|f| f.values.clone())
}

/// Solves `targets` (and nothing else) in coupling order. The result of each is returned in the order it was solved.
pub fn run_chain(
    registry: &Registry,
    project: &Project,
    model: &MeshedModel,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
) -> Result<Vec<SolveResult>, Error> {
    let ids: Vec<&str> = targets.iter().map(String::as_str).collect();
    let order = registry.solve_order(&ids)?;
    let mut done: Vec<(String, SolveResult)> = vec![];
    for id in order {
        if progress.is_cancelled() {
            return Err(Error::Cancelled);
        }
        let p = physics_for(registry, &id)?;
        let spec = project
            .spec(&id)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{id}'; add the physics first")))?;
        validate_spec(p, spec)?;
        p.validate(spec, model.mesh())?;
        let mut ctx = Ctx::new(model);
        ctx.progress = progress.clone();
        ctx.library = Some(library.clone());
        ctx.wavelength = crate::materials::study_wavelength(p, spec);
        for imp in p.imports() {
            if let Some(values) = done.iter().find_map(|(o, r)| exported(registry, o, r, &imp.role)) {
                ctx.inputs.insert(imp.role.clone(), values);
            }
        }
        let result = p.solve(model.mesh(), spec, &mut ctx)?;
        done.push((id, result));
    }
    Ok(done.into_iter().map(|(_, r)| r).collect())
}
