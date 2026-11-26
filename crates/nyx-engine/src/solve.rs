// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
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

/// A kernel knows the regions by number and names a scalar of one region `... domain <number> ...` (`share in domain 3 1`); the model knows
/// the names, so the scalar becomes `share in disk 1` (ADR 80). A number that is not a region of the model is left as it is.
pub(crate) fn name_domains(result: &mut SolveResult, model: &MeshedModel) {
    let regions = model.regions();
    for (name, _) in result.scalars.iter_mut() {
        let words: Vec<&str> = name.split(' ').collect();
        let Some(k) = words.iter().position(|w| *w == "domain") else { continue };
        let Some(region) =
            words.get(k + 1).and_then(|w| w.parse::<i32>().ok()).and_then(|id| regions.iter().find(|(i, _)| *i == id))
        else {
            continue;
        };
        let mut out: Vec<&str> = words[..k].to_vec();
        out.push(&region.1);
        out.extend_from_slice(&words[k + 2..]);
        *name = out.join(" ");
    }
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
    match run_chain_partial(registry, project, model, targets, library, progress) {
        (results, None) => Ok(results),
        (_, Some(e)) => Err(e),
    }
}

/// Like `run_chain`, but when a physics fails the physics solved before it are returned with the error: a failed `sms` keeps the `ems` result
/// it took its field from (a failure fails the physics after it in the chain, not the ones before).
pub fn run_chain_partial(
    registry: &Registry,
    project: &Project,
    model: &MeshedModel,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
) -> (Vec<SolveResult>, Option<Error>) {
    run_chain_partial_with(registry, project, model, targets, library, progress, project.mesh.retain_amr_data)
}

pub(crate) fn run_chain_partial_with(
    registry: &Registry,
    project: &Project,
    model: &MeshedModel,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
    retain_amr_data: bool,
) -> (Vec<SolveResult>, Option<Error>) {
    let ids: Vec<&str> = targets.iter().map(String::as_str).collect();
    let order = match registry.solve_order(&ids) {
        Ok(o) => o,
        Err(e) => return (vec![], Some(e)),
    };
    let mut done: Vec<(String, SolveResult)> = vec![];
    let mut failed = None;
    for id in order {
        match solve_one_single(registry, project, model, &id, library, &progress, &done) {
            Ok(mut result) => {
                if !retain_amr_data {
                    nyx_physics::ElectroAmrData::strip_from(&mut result);
                }
                done.push((id, result));
            }
            Err(e) => {
                failed = Some(e);
                break;
            }
        }
    }
    (done.into_iter().map(|(_, r)| r).collect(), failed)
}

/// Like [`run_chain_partial`], and runs multi-pass AMR (`dfs`, `ems::cavity`, or `sms`) when configured,
/// returning the final adapted [`MeshedModel`] alongside the solved results.
pub fn run_chain_partial_amr(
    registry: &Registry,
    project: &Project,
    model: Arc<MeshedModel>,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
) -> (Arc<MeshedModel>, Vec<SolveResult>, Option<Error>) {
    run_chain_partial_amr_with(registry, project, model, targets, library, progress, project.mesh.retain_amr_data)
}

pub(crate) fn run_chain_partial_amr_with(
    registry: &Registry,
    project: &Project,
    mut model: Arc<MeshedModel>,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
    retain_amr_data: bool,
) -> (Arc<MeshedModel>, Vec<SolveResult>, Option<Error>) {
    let ids: Vec<&str> = targets.iter().map(String::as_str).collect();
    let order = match registry.solve_order(&ids) {
        Ok(o) => o,
        Err(e) => return (model, vec![], Some(e)),
    };
    let use_electro_amr = crate::amr::project_uses_electro_amr(project, targets);
    let use_ems_amr = crate::amr::project_uses_ems_cavity_amr(project, targets);
    let use_sms_amr = crate::amr::project_uses_sms_amr(project, targets);
    let mut done: Vec<(String, SolveResult)> = vec![];
    let mut failed = None;
    for id in order {
        let step = if id == "dfs" && use_electro_amr {
            crate::amr::solve_dfs_with_amr_opts(
                registry,
                project,
                model.clone(),
                library,
                &progress,
                &done,
                retain_amr_data,
            )
        } else if id == "ems" && use_ems_amr {
            crate::amr::solve_ems_cavity_with_amr(registry, project, model.clone(), library, &progress, &done)
        } else if id == "sms" && use_sms_amr {
            crate::amr::solve_sms_with_amr(registry, project, model.clone(), library, &progress, &done)
        } else {
            solve_one_single(registry, project, &model, &id, library, &progress, &done).map(|mut r| {
                if !retain_amr_data {
                    nyx_physics::ElectroAmrData::strip_from(&mut r);
                }
                (model.clone(), r)
            })
        };
        match step {
            Ok((next_model, result)) => {
                model = next_model;
                done.push((id, result));
            }
            Err(e) => {
                failed = Some(e);
                break;
            }
        }
    }
    (model, done.into_iter().map(|(_, r)| r).collect(), failed)
}

/// Builds a [`Ctx`] with imported fields from previously solved physics in `done`.
pub(crate) fn build_solve_ctx<'a>(
    registry: &Registry,
    model: &'a MeshedModel,
    p: &dyn nyx_physics::Physics,
    spec: &nyx_physics::ProblemSpec,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Ctx<'a> {
    let mut ctx = Ctx::new(model);
    ctx.progress = progress.clone();
    ctx.library = Some(library.clone());
    ctx.wavelength = crate::materials::study_wavelength(p, spec);
    ctx.frequency = crate::materials::study_frequency(p, spec);
    for imp in p.imports() {
        if let Some(values) = done.iter().find_map(|(o, r)| exported(registry, o, r, &imp.role)) {
            ctx.inputs.insert(imp.role.clone(), values);
        }
    }
    ctx
}

/// One physics of a chain, with the fields of the physics solved before it.
pub(crate) fn solve_one_single(
    registry: &Registry,
    project: &Project,
    model: &MeshedModel,
    id: &str,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<SolveResult, Error> {
    if progress.is_cancelled() {
        return Err(Error::Cancelled);
    }
    let p = physics_for(registry, id)?;
    let spec = project
        .spec(id)
        .ok_or_else(|| Error::invalid(format!("the project has no problem for '{id}'; add the physics first")))?;
    if id == "dfs"
        && spec.study.kind == "electrostatic"
        && project.mesh.singular_edges.as_deref() == Some("quarter_point")
    {
        let fb = spec.study.params.get("flux_bounds").and_then(|v| v.as_str()).unwrap_or("none");
        let ps_amr = spec.study.params.get("adapt_passes").and_then(|v| v.as_i64()).unwrap_or(1) > 1
            && spec.study.params.get("adapt_estimator").and_then(|v| v.as_str()) == Some("prager_synge");
        if fb != "none" || ps_amr {
            return Err(Error::invalid(
                "flux_bounds / prager_synge cannot be combined with singular_edges='quarter_point': contravariant Piola-mapped Raviart-Thomas fluxes have divergent L^2 energy on singular quarter-point (det J -> 0) CurvedTets",
            ));
        }
    }
    validate_spec(p, spec)?;
    p.validate(spec, model.mesh())?;
    let mut ctx = build_solve_ctx(registry, model, p, spec, library, progress, done);
    let mut result = p.solve(model.mesh(), spec, &mut ctx)?;
    name_domains(&mut result, model);
    Ok(result)
}
