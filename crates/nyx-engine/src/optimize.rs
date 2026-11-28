// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Discrete adjoint sensitivity evaluation (`sensitivity_run`) and gradient-based shape, material,
//! and PDE-filtered topology optimization (`optimize_run`, wiring `femcore::linalg::opt::{lbfgs,
//! nonlinear_cg, mma}`, ADR 142).
//!
//! Nothing here names a specific physics crate: every forward solve and discrete adjoint evaluation
//! dispatches through `dyn nyx_physics::Physics` via `Engine::registry()`.

use crate::engine::{Engine, Stored};
use crate::geometry;
use crate::meshing::MeshedModel;
use crate::physics_for;
use crate::store_link::RunOutcome;
use femcore::fem::geometry as simplex_geometry;
use femcore::linalg::opt::{lbfgs, mma, nonlinear_cg, MinOptions, MmaEval, MmaOptions, StopReason};
use femcore::progress::Progress;
use femcore::topology::{compute_shape_velocity, displace_mesh, HelmholtzFilter, ShapeVelocityField};
use nyx_physics::{
    validate_spec, CellField, ConstraintSpec, ContinuationStage, DesignParam, Error, LabelResolver, MaterialAssign,
    OptimizeSense, OptimizeSpec, ProblemSpec, PropValue, Scalar, SensitivityObjective, SensitivityReport,
    SensitivityRequest, SolveResult, TrackedModeBasis,
};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

fn stop_reason_name(r: StopReason) -> &'static str {
    match r {
        StopReason::GradientTolerance => "gradient_tolerance",
        StopReason::FunctionTolerance => "function_tolerance",
        StopReason::StepTolerance => "step_tolerance",
        StopReason::MaxIterations => "max_iterations",
        StopReason::LineSearchFailed => "line_search_failed",
        StopReason::KktTolerance => "kkt_tolerance",
    }
}

/// Metadata describing how one entry in `SensitivityRequest::params` maps onto the flat optimizer vector `x`.
#[derive(Clone)]
enum ParamSlot {
    ParamVar { name: String, offset: usize, base_val: f64, min: Option<f64>, max: Option<f64> },
    MaterialScalar { region: String, property: String, offset: usize, min: Option<f64>, max: Option<f64> },
    EmeSegmentLength { segment: usize, offset: usize, min: Option<f64>, max: Option<f64> },
    Topology { param_idx: usize, domain: String, offset: usize, len: usize, filter: Arc<HelmholtzFilter> },
}

fn read_and_promote_material_scalar(
    p_phys: &dyn nyx_physics::Physics,
    spec: &mut ProblemSpec,
    model: &MeshedModel,
    library: &Arc<femcore::material::Library>,
    region: &str,
    property: &str,
) -> Result<f64, Error> {
    if let Some(assign) = spec.materials.get(region) {
        match assign {
            MaterialAssign::Props(map) => {
                if let Some(&v) = map.get(property) {
                    return Ok(v);
                }
            }
            MaterialAssign::PropsTensor(map) => {
                if let Some(PropValue::Scalar(v)) = map.get(property) {
                    return Ok(*v);
                }
            }
            MaterialAssign::Named(_) | MaterialAssign::Oriented { .. } => {}
        }
    }
    let mut ctx = nyx_physics::Ctx::new(model);
    ctx.library = Some(library.clone());
    let mats = nyx_physics::decode::materials_by_label(spec, &ctx, &p_phys.properties())?;
    let doms = nyx_physics::LabelResolver::domains(model, region)?;
    if let Some(&d0) = doms.first() {
        if let Some(v) = mats.get(d0, property) {
            if !matches!(spec.materials.get(region), Some(MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_))) {
                let mut inline = mats.domains.get(&d0).cloned().unwrap_or_default();
                inline.insert(property.to_string(), v);
                spec.materials.insert(region.to_string(), MaterialAssign::Props(inline));
            }
            return Ok(v);
        }
    }
    Err(Error::invalid(format!(
        "MaterialScalar '{region}.{property}' has no initial scalar value on region '{region}' (assign an inline scalar property first)"
    )))
}

fn set_spec_material_scalar(spec: &mut ProblemSpec, region: &str, property: &str, val: f64) {
    match spec.materials.get_mut(region) {
        Some(MaterialAssign::Props(map)) => {
            map.insert(property.to_string(), val);
        }
        Some(MaterialAssign::PropsTensor(map)) => {
            map.insert(property.to_string(), PropValue::Scalar(val));
        }
        Some(MaterialAssign::Named(_)) | Some(MaterialAssign::Oriented { .. }) | None => {
            let mut props = BTreeMap::new();
            props.insert(property.to_string(), val);
            spec.materials.insert(region.to_string(), MaterialAssign::Props(props));
        }
    }
}

fn parse_optimize_spec(args: &Value) -> Result<OptimizeSpec, Error> {
    let mut opt: OptimizeSpec = match args.get("optimize").filter(|v| !v.is_null()) {
        Some(v) => serde_json::from_value(v.clone()).map_err(|e| Error::invalid(format!("`optimize`: {e}")))?,
        None => OptimizeSpec::default(),
    };
    if let Some(alg) = args.get("algorithm").and_then(Value::as_str) {
        opt.algorithm = alg.to_string();
    }
    if let Some(s) = args.get("sense").filter(|v| !v.is_null()) {
        opt.sense = serde_json::from_value(s.clone()).map_err(|e| Error::invalid(format!("`sense`: {e}")))?;
    }
    if let Some(n) = args.get("max_iters").and_then(Value::as_u64) {
        opt.max_iters = n as usize;
    }
    if let Some(g) = args.get("gtol").and_then(Value::as_f64) {
        opt.gtol = g;
    }
    if let Some(f) = args.get("ftol").and_then(Value::as_f64) {
        opt.ftol = f;
    }
    if let Some(xt) = args.get("xtol").and_then(Value::as_f64) {
        opt.xtol = xt;
    }
    if let Some(ml) = args.get("move_limit").and_then(Value::as_f64) {
        opt.move_limit = ml;
    }
    if let Some(cs) = args.get("constraints").filter(|v| !v.is_null()) {
        opt.constraints =
            serde_json::from_value(cs.clone()).map_err(|e| Error::invalid(format!("`constraints`: {e}")))?;
    }
    if let Some(sched) = args.get("continuation_schedule").filter(|v| !v.is_null()) {
        opt.continuation_schedule = serde_json::from_value(sched.clone())
            .map_err(|e| Error::invalid(format!("`continuation_schedule`: {e}")))?;
    }
    opt.validate_finite()?;
    Ok(opt)
}

fn parse_objective_and_params(args: &Value) -> Result<(SensitivityObjective, Vec<DesignParam>), Error> {
    let obj_val =
        args.get("objective").filter(|v| !v.is_null()).ok_or_else(|| Error::invalid("`objective` is required"))?;
    let objective: SensitivityObjective =
        serde_json::from_value(obj_val.clone()).map_err(|e| Error::invalid(format!("`objective`: {e}")))?;
    let params_val = args
        .get("params")
        .filter(|v| !v.is_null())
        .ok_or_else(|| Error::invalid("`params` is required (a non-empty list of design parameters)"))?;
    let params: Vec<DesignParam> =
        serde_json::from_value(params_val.clone()).map_err(|e| Error::invalid(format!("`params`: {e}")))?;
    Ok((objective, params))
}

impl Engine {
    fn resolve_target_physics(&self, args: &Value) -> Result<String, Error> {
        if self.project.physics.is_empty() {
            return Err(Error::invalid("the project has no physics; add one first (add_physics)"));
        }
        let target = match args.get("physics").filter(|v| !v.is_null()) {
            Some(v) => {
                v.as_str().ok_or_else(|| Error::invalid("`physics` must be a three-letter physics code"))?.to_string()
            }
            None => self
                .project
                .active
                .clone()
                .or_else(|| self.project.physics.first().map(|s| s.physics.clone()))
                .ok_or_else(|| Error::invalid("the project has no active physics"))?,
        };
        if self.project.spec(&target).is_none() {
            return Err(Error::invalid(format!("the project has no problem for '{target}'")));
        }
        Ok(target)
    }

    fn ensure_base_model(&mut self) -> Result<Arc<MeshedModel>, Error> {
        self.project.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        if let Some(m) = self.model() {
            return Ok(m.clone());
        }
        let key = self.mesh_key();
        if let Some(m) = self.stages.mesh(&key) {
            self.model = Some((key, m.clone()));
            return Ok(m);
        }
        let sizing = self.sizing();
        let mut m = self.mesher.mesh(&self.project.geometry, &self.project.mesh, sizing.as_ref())?;
        if let Some(active) = &self.project.mesh.active_domains {
            m = m.restrict_to_active(active)?;
        }
        let arc = Arc::new(m);
        self.stages.put_mesh(&key, &arc);
        self.model = Some((key, arc.clone()));
        Ok(arc)
    }

    fn ensure_physics_model(&mut self, physics_id: &str) -> Result<Arc<MeshedModel>, Error> {
        self.project.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        if self.project.has_own_mesh(physics_id) {
            let settings = self.project.mesh_for(physics_id).clone();
            let lib = self.library().ok();
            let sizing = if settings.physics_sizing {
                self.registry
                    .get(physics_id)
                    .and_then(|p| self.project.spec(physics_id).map(|sp| p.sizing(sp, lib.as_ref())))
                    .filter(|s| !s.wavelengths.is_empty())
            } else {
                None
            };
            let mut m = self.mesher.mesh(&self.project.geometry, &settings, sizing.as_ref())?;
            if let Some(active) = &settings.active_domains {
                m = m.restrict_to_active(active)?;
            }
            return Ok(Arc::new(m));
        }
        self.ensure_base_model()
    }

    #[allow(clippy::type_complexity)]
    fn upstream_done_for(
        &mut self,
        physics_id: &str,
        library: &Arc<femcore::material::Library>,
        progress: &Arc<Progress>,
    ) -> Result<(Vec<(String, SolveResult)>, BTreeMap<String, Arc<MeshedModel>>), Error> {
        let target_p = physics_for(&self.registry, physics_id)?;
        let target_study = self.project.spec(physics_id).map(|s| s.study.kind.as_str()).unwrap_or("");
        let target_roles: Vec<String> = target_p.study_imports(target_study).into_iter().map(|l| l.role).collect();
        if target_roles.is_empty() {
            return Ok((Vec::new(), BTreeMap::new()));
        }

        // Collect stationary/static upstream physics in the project whose exported roles feed `physics_id`
        let mut upstream_ids: Vec<&str> = Vec::new();
        for sp in &self.project.physics {
            if sp.physics == physics_id {
                continue;
            }
            if !matches!(sp.study.kind.as_str(), "electrostatic" | "heat_stationary" | "conduction" | "static") {
                continue;
            }
            if let Some(src_p) = self.registry.get(&sp.physics) {
                let exports_needed =
                    src_p.study_exports(&sp.study.kind).iter().any(|ex| target_roles.contains(&ex.role));
                if exports_needed {
                    upstream_ids.push(sp.physics.as_str());
                }
            }
        }
        if upstream_ids.is_empty() {
            return Ok((Vec::new(), BTreeMap::new()));
        }

        let mut all_ids = upstream_ids.clone();
        all_ids.push(physics_id);
        let order = crate::solve::project_solve_order(&self.registry, &self.project, &all_ids)?;
        let upstream: Vec<String> = order.into_iter().filter(|id| id != physics_id).collect();
        if upstream.is_empty() {
            return Ok((Vec::new(), BTreeMap::new()));
        }

        // Reuse fresh stored results & per-physics models when already solved and non-stale
        if upstream.iter().all(|id| self.results.contains_key(id) && !self.is_stale(id)) {
            let mut done = Vec::with_capacity(upstream.len());
            let mut models = BTreeMap::new();
            for id in upstream {
                if let Some(st) = self.results.get(&id) {
                    models.insert(id.clone(), st.model.clone());
                    done.push((id, (*st.result).clone()));
                }
            }
            return Ok((done, models));
        }

        let proj_base = self.ensure_base_model()?;
        let out = crate::solve::run_chain_models(
            &self.registry,
            &self.project,
            proj_base,
            &upstream,
            library,
            progress.clone(),
            false,
            &*self.mesher,
        );
        if let Some(e) = out.failed {
            return Err(e);
        }
        let done = out.results.into_iter().map(|r| (r.physics.clone(), r)).collect();
        Ok((done, out.models))
    }

    /// Evaluates exact discrete adjoint sensitivities (`J(p)` and `dJ/dp`) on the current project.
    pub fn sensitivity_run(&mut self, args: &Value) -> Result<Value, Error> {
        let t0 = Instant::now();
        let target = self.resolve_target_physics(args)?;
        self.project.check(&self.registry)?;
        let (objective, params) = parse_objective_and_params(args)?;
        let mut req = SensitivityRequest::new(objective, params);
        req.validate_finite()?;

        let step_rel = args.get("shape_step_rel").and_then(Value::as_f64);
        let include_densities = args.get("include_densities").and_then(Value::as_bool).unwrap_or(true);

        let model = self.ensure_physics_model(&target)?;
        let library = self.library()?;
        let progress = Progress::new();
        let (done, chain_models) = self.upstream_done_for(&target, &library, &progress)?;

        for param in &req.params {
            if let DesignParam::ParamVar { name, .. } = param {
                if !req.shape_velocities.contains_key(name) {
                    let vel =
                        compute_shape_velocity(&self.project.geometry, model.evaluated(), model.mesh(), name, step_rel)
                            .map_err(Error::invalid)?;
                    req.shape_velocities.insert(name.clone(), vel);
                }
            }
        }

        let p = physics_for(&self.registry, &target)?;
        let spec = self
            .project
            .spec(&target)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{target}'")))?;
        validate_spec(p, spec)?;
        let report = {
            let _chain = crate::solve::install_chain_models(chain_models);
            let mut ctx = crate::solve::build_solve_ctx(&self.registry, &model, p, spec, &library, &progress, &done)?;
            p.sensitivity(model.mesh(), spec, &mut ctx, &req)?
        };
        report.check_finite()?;

        let mut out = json!({
            "physics": target,
            "objective": report.objective,
            "gradients": report.gradients,
            "flat_gradient": report.flat_gradient,
            "param_names": report.param_names,
            "n_design_vars": report.flat_gradient.len(),
            "aux_scalars": report.aux_scalars,
            "volume_fractions": report.volume_fractions,
            "seconds": t0.elapsed().as_secs_f64(),
            "log": report.log,
        });
        if include_densities && !report.projected_densities.is_empty() {
            out["projected_densities"] = json!(report.projected_densities);
        }
        Ok(out)
    }

    /// Runs gradient-based shape, material, or PDE-filtered topology optimization (`lbfgs`, `ncg`, `mma`),
    /// optionally across a multi-stage Heaviside `beta`-continuation schedule, updates the project (unless
    /// `apply: false`), stores the final solved fields (including `"topology_density"` when applicable),
    /// and optionally writes an optimization checkpoint JSON file (`checkpoint_path`).
    pub fn optimize_run(&mut self, args: &Value) -> Result<Value, Error> {
        let t0 = Instant::now();
        let target = self.resolve_target_physics(args)?;
        self.project.check(&self.registry)?;
        if self.auto_store && self.link.is_none() && self.workspace.is_some() {
            self.save_to_workspace(None)?;
        }

        let (objective, mut params) = parse_objective_and_params(args)?;
        let probe_req = SensitivityRequest::new(objective.clone(), params.clone());
        probe_req.validate_finite()?;
        let opt_spec = parse_optimize_spec(args)?;

        let remesh = args.get("remesh").and_then(Value::as_bool).unwrap_or(false);
        let apply_project = args.get("apply").and_then(Value::as_bool).unwrap_or(true);
        let step_rel = args.get("shape_step_rel").and_then(Value::as_f64);
        let checkpoint_path = args
            .get("checkpoint_path")
            .filter(|v| !v.is_null())
            .map(|v| {
                v.as_str()
                    .filter(|s| !s.trim().is_empty())
                    .map(String::from)
                    .ok_or_else(|| Error::invalid("`checkpoint_path` must be a non-empty file path string"))
            })
            .transpose()?;

        let has_topology = params.iter().any(|p| matches!(p, DesignParam::TopologyDensity { .. }));
        if has_topology && remesh {
            return Err(Error::invalid(
                "`remesh: true` cannot be combined with `topology_density` parameters (cellwise densities require a fixed mesh)",
            ));
        }

        let registry = self.registry.clone();
        let p_phys = physics_for(&registry, &target)?;
        let base_model = self.ensure_physics_model(&target)?;
        let library = self.library()?;
        let progress = Progress::new();
        let (done, chain_models) = self.upstream_done_for(&target, &library, &progress)?;
        let mut base_spec = self
            .project
            .spec(&target)
            .ok_or_else(|| Error::invalid(format!("the project has no problem for '{target}'")))?
            .clone();

        // Build parameter slots, initial vector `x0`, box bounds `[xmin, xmax]`, cached Helmholtz filters,
        // and base harmonic shape velocities.
        let simplex0 = base_model.solver_mesh();
        let mut slots: Vec<ParamSlot> = Vec::with_capacity(params.len());
        let mut x0: Vec<f64> = Vec::new();
        let mut xmin: Vec<f64> = Vec::new();
        let mut xmax: Vec<f64> = Vec::new();
        let mut base_velocities: BTreeMap<String, ShapeVelocityField> = BTreeMap::new();
        let mut topology_filters: BTreeMap<usize, Arc<HelmholtzFilter>> = BTreeMap::new();
        let mut base_eme_lengths: Vec<f64> = Vec::new();
        let mut cached_eme_data: Option<nyx_ems::EmeData> = None;
        let pure_eme_lengths =
            !params.is_empty() && params.iter().all(|p| matches!(p, DesignParam::EmeSegmentLength { .. }));
        if params.iter().any(|p| matches!(p, DesignParam::EmeSegmentLength { .. })) {
            if target != "ems" {
                return Err(Error::invalid(format!(
                    "eme_segment_length applies only to physics 'ems' (study 'eme'), not '{target}'"
                )));
            }
            let _chain = crate::solve::install_chain_models(chain_models.clone());
            let mut ctx0 =
                crate::solve::build_solve_ctx(&registry, &base_model, p_phys, &base_spec, &library, &progress, &done)?;
            base_eme_lengths = nyx_ems::eme_nominal_segment_lengths(base_model.mesh(), &base_spec, &ctx0)?;
            if pure_eme_lengths {
                // Solve the 2D cross-section eigenmodes & interface overlaps ONCE: subsequent optimizer steps evaluate
                // Redheffer star products & analytical dP/dL tangents directly from `EmeData` in microseconds!
                let solved0 = p_phys.solve(base_model.mesh(), &base_spec, &mut ctx0)?;
                if let Some(edata) = nyx_ems::EmeData::of(&solved0) {
                    base_eme_lengths = edata.segment_lengths();
                    cached_eme_data = Some(edata);
                }
            }
        }

        for (idx, param) in params.iter().enumerate() {
            match param {
                DesignParam::ParamVar { name, min, max } => {
                    let val0 = base_model
                        .evaluated()
                        .scene()
                        .env
                        .get(name)
                        .ok_or_else(|| Error::invalid(format!("geometry parameter '{name}' does not exist")))?;
                    let vel = compute_shape_velocity(
                        &self.project.geometry,
                        base_model.evaluated(),
                        base_model.mesh(),
                        name,
                        step_rel,
                    )
                    .map_err(Error::invalid)?;
                    base_velocities.insert(name.clone(), vel);
                    let span = val0.abs().max(base_model.evaluated().bounds().diagonal()).max(1e-6);
                    let lo = min.unwrap_or(val0 - 0.5 * span);
                    let hi = max.unwrap_or(val0 + 0.5 * span);
                    let start = val0.clamp(lo, hi);
                    let offset = x0.len();
                    x0.push(start);
                    xmin.push(lo);
                    xmax.push(hi);
                    slots.push(ParamSlot::ParamVar {
                        name: name.clone(),
                        offset,
                        base_val: val0,
                        min: *min,
                        max: *max,
                    });
                }
                DesignParam::MaterialScalar { region, property, min, max } => {
                    let val0 = read_and_promote_material_scalar(
                        p_phys,
                        &mut base_spec,
                        &base_model,
                        &library,
                        region,
                        property,
                    )?;
                    let span = val0.abs().max(1.0);
                    let lo = min.unwrap_or((val0 - 0.9 * span).max(-1e12));
                    let hi = max.unwrap_or(val0 + 10.0 * span);
                    let start = val0.clamp(lo, hi);
                    let offset = x0.len();
                    x0.push(start);
                    xmin.push(lo);
                    xmax.push(hi);
                    slots.push(ParamSlot::MaterialScalar {
                        region: region.clone(),
                        property: property.clone(),
                        offset,
                        min: *min,
                        max: *max,
                    });
                }
                DesignParam::EmeSegmentLength { segment, min, max } => {
                    let val0 = *base_eme_lengths.get(*segment).ok_or_else(|| {
                        Error::invalid(format!(
                            "eme_segment_length segment {segment} out of range (cascade has {} segment(s))",
                            base_eme_lengths.len()
                        ))
                    })?;
                    let lo = min.unwrap_or((0.01 * val0).max(1e-12)).max(1e-12);
                    let hi = max.unwrap_or((20.0 * val0).max(lo + 1e-6));
                    let start = val0.clamp(lo, hi);
                    let offset = x0.len();
                    x0.push(start);
                    xmin.push(lo);
                    xmax.push(hi);
                    slots.push(ParamSlot::EmeSegmentLength { segment: *segment, offset, min: *min, max: *max });
                }
                DesignParam::TopologyDensity { domain, filter_radius, initial_density, densities, .. } => {
                    let dom_labels = base_model.domains(domain)?;
                    if dom_labels.is_empty() {
                        return Err(Error::invalid(format!(
                            "TopologyDensity domain '{domain}' resolved to zero mesh domains"
                        )));
                    }
                    let filter = Arc::new(
                        HelmholtzFilter::new(&simplex0, *filter_radius, Some(&dom_labels))
                            .map_err(Error::from_kernel)?,
                    );
                    let nc = filter.n_design_cells();
                    let default_d0 = initial_density.unwrap_or_else(|| {
                        opt_spec
                            .constraints
                            .iter()
                            .find_map(|c| match c {
                                ConstraintSpec::VolumeFraction { param_index, max_fraction, .. }
                                    if *param_index == idx =>
                                {
                                    Some(max_fraction.clamp(0.01, 0.99))
                                }
                                _ => None,
                            })
                            .unwrap_or(0.5)
                    });
                    let init_vec = if densities.is_empty() {
                        vec![default_d0.clamp(0.0, 1.0); nc]
                    } else if densities.len() == nc {
                        densities.iter().map(|&d| d.clamp(0.0, 1.0)).collect()
                    } else {
                        return Err(Error::invalid(format!(
                            "TopologyDensity densities length ({}) does not match design cell count ({nc}) in '{domain}'",
                            densities.len()
                        )));
                    };
                    let offset = x0.len();
                    x0.extend_from_slice(&init_vec);
                    xmin.extend(std::iter::repeat_n(0.0, nc));
                    xmax.extend(std::iter::repeat_n(1.0, nc));
                    topology_filters.insert(idx, filter.clone());
                    slots.push(ParamSlot::Topology { param_idx: idx, domain: domain.clone(), offset, len: nc, filter });
                }
            }
        }

        let n_vars = x0.len();
        if n_vars == 0 {
            return Err(Error::invalid("optimize_run resolved zero scalar or cell design variables"));
        }

        let sense_sign = match opt_spec.sense {
            OptimizeSense::Minimize => 1.0,
            OptimizeSense::Maximize => -1.0,
        };

        let stages: Vec<ContinuationStage> = if opt_spec.continuation_schedule.is_empty() {
            let b0 = params
                .iter()
                .find_map(|p| match p {
                    DesignParam::TopologyDensity { beta, .. } => Some(*beta),
                    _ => None,
                })
                .unwrap_or(1.0);
            vec![ContinuationStage { beta: b0, max_iters: opt_spec.max_iters }]
        } else {
            opt_spec.continuation_schedule.clone()
        };

        let base_geometry = self.project.geometry.clone();
        let base_mesh_settings = self.project.mesh_for(&target).clone();
        let mesher = self.mesher.clone();

        // Shared mutable state across optimizer callback invocations.
        let mut tracked_modes: Option<TrackedModeBasis> = None;
        let mut last_good_report: Option<(Vec<f64>, Vec<DesignParam>, SensitivityReport)> = None;
        let mut eval_error: Option<Error> = None;
        let mut total_evals: usize = 0;

        // Evaluates `(physical_report)` at candidate `x_cand` for current `params_stage`.
        let mut evaluate_candidate = |x_cand: &[f64],
                                      params_stage: &[DesignParam]|
         -> Result<SensitivityReport, Error> {
            if let Some((ref x_prev, ref p_prev, ref rep_prev)) = last_good_report {
                if x_prev == x_cand && p_prev == params_stage {
                    return Ok(rep_prev.clone());
                }
            }
            total_evals += 1;
            if x_cand.iter().any(|v| !v.is_finite()) {
                return Err(Error::NonFinite("optimizer proposed non-finite design coordinates".into()));
            }
            // Check explicit box bounds for unconstrained line-search trial steps (`lbfgs` / `ncg`).
            for slot in &slots {
                match slot {
                    ParamSlot::ParamVar { offset, min, max, .. }
                    | ParamSlot::MaterialScalar { offset, min, max, .. }
                    | ParamSlot::EmeSegmentLength { offset, min, max, .. } => {
                        let v = x_cand[*offset];
                        if matches!(slot, ParamSlot::EmeSegmentLength { .. }) && v <= 0.0 {
                            return Err(Error::invalid("trial EME segment length must be positive"));
                        }
                        if let Some(lo) = min {
                            if v < *lo - 1e-15 {
                                return Err(Error::invalid("trial step below min bound"));
                            }
                        }
                        if let Some(hi) = max {
                            if v > *hi + 1e-15 {
                                return Err(Error::invalid("trial step above max bound"));
                            }
                        }
                    }
                    ParamSlot::Topology { offset, len, .. } => {
                        for &v in &x_cand[*offset..*offset + *len] {
                            if !(-1e-12..=1.0 + 1e-12).contains(&v) {
                                return Err(Error::invalid("trial density outside [0, 1]"));
                            }
                        }
                    }
                }
            }

            let mut spec_cand = base_spec.clone();
            let mut params_cand = params_stage.to_vec();
            let mut eme_overrides: Vec<(usize, f64)> = Vec::new();
            for slot in &slots {
                match slot {
                    ParamSlot::ParamVar { .. } => {}
                    ParamSlot::MaterialScalar { region, property, offset, .. } => {
                        set_spec_material_scalar(&mut spec_cand, region, property, x_cand[*offset]);
                    }
                    ParamSlot::EmeSegmentLength { segment, offset, .. } => {
                        eme_overrides.push((*segment, x_cand[*offset]));
                    }
                    ParamSlot::Topology { param_idx, offset, len, .. } => {
                        if let DesignParam::TopologyDensity { densities, .. } = &mut params_cand[*param_idx] {
                            *densities = x_cand[*offset..*offset + *len].iter().map(|&d| d.clamp(0.0, 1.0)).collect();
                        }
                    }
                }
            }

            if let Some(ref edata) = cached_eme_data {
                let mut cand_lengths = base_eme_lengths.clone();
                for &(seg, val) in &eme_overrides {
                    if let Some(slot_l) = cand_lengths.get_mut(seg) {
                        *slot_l = val;
                    }
                }
                let rep = nyx_ems::evaluate_eme_length_sensitivity_from_data(
                    edata,
                    &objective,
                    &params_cand,
                    Some(&cand_lengths),
                )?;
                rep.check_finite()?;
                last_good_report = Some((x_cand.to_vec(), params_stage.to_vec(), rep.clone()));
                return Ok(rep);
            }

            // Build candidate mesh + shape velocities (`remesh: false` smoothly morphs `base_model.mesh()`).
            type CandGeom = (nyx_physics::Mesh, Option<Arc<MeshedModel>>, BTreeMap<String, ShapeVelocityField>);
            let (cand_mesh_holder, cand_model_opt, shape_vels): CandGeom = if remesh {
                let mut geom_cand = base_geometry.clone();
                for slot in &slots {
                    if let ParamSlot::ParamVar { name, offset, .. } = slot {
                        geometry::set_param(&mut geom_cand, name, &json!(x_cand[*offset]))?;
                    }
                }
                let mut m = mesher.mesh(&geom_cand, &base_mesh_settings, None)?;
                if let Some(active) = &base_mesh_settings.active_domains {
                    m = m.restrict_to_active(active)?;
                }
                let mut vels = BTreeMap::new();
                for slot in &slots {
                    if let ParamSlot::ParamVar { name, .. } = slot {
                        let v = compute_shape_velocity(&geom_cand, m.evaluated(), m.mesh(), name, step_rel)
                            .map_err(Error::invalid)?;
                        vels.insert(name.clone(), v);
                    }
                }
                let mesh_clone = (**m.mesh()).clone();
                (mesh_clone, Some(Arc::new(m)), vels)
            } else {
                let mut cur_mesh = (**base_model.mesh()).clone();
                for slot in &slots {
                    if let ParamSlot::ParamVar { name, offset, base_val, .. } = slot {
                        let delta = x_cand[*offset] - *base_val;
                        if delta.abs() > 0.0 {
                            if let Some(vel) = base_velocities.get(name) {
                                cur_mesh = displace_mesh(&cur_mesh, delta, vel).map_err(Error::invalid)?;
                            }
                        }
                    }
                }
                // Guard against element inversion during aggressive line-search trial steps.
                let simp = cur_mesh.discretise().map_err(Error::Invalid)?;
                let g = simplex_geometry(&simp);
                if g.vol.iter().any(|&v| !v.is_finite() || v <= 1e-30) {
                    return Err(Error::invalid("trial step inverted or collapsed a mesh element"));
                }
                (cur_mesh, None, base_velocities.clone())
            };

            let resolver_model: &MeshedModel = match &cand_model_opt {
                Some(m) => m.as_ref(),
                None => base_model.as_ref(),
            };
            let rep = {
                let _chain = crate::solve::install_chain_models(chain_models.clone());
                let mut ctx = crate::solve::build_solve_ctx(
                    &registry,
                    resolver_model,
                    p_phys,
                    &spec_cand,
                    &library,
                    &progress,
                    &done,
                )?;
                if !eme_overrides.is_empty() {
                    nyx_ems::eme_set_spec_segment_lengths(&mut spec_cand, &ctx, &base_eme_lengths, &eme_overrides)?;
                }

                let mut req = SensitivityRequest::new(objective.clone(), params_cand);
                req.shape_velocities = shape_vels;
                req.topology_filters = topology_filters.clone();
                req.reference_modes = tracked_modes.clone();

                p_phys.sensitivity(&cand_mesh_holder, &spec_cand, &mut ctx, &req)?
            };
            rep.check_finite()?;
            if rep.flat_gradient.len() != n_vars {
                return Err(Error::Solver(format!(
                    "sensitivity flat_gradient length ({}) != optimizer vector length ({n_vars})",
                    rep.flat_gradient.len()
                )));
            }
            if rep.mode_basis.is_some() {
                tracked_modes = rep.mode_basis.clone();
            }
            last_good_report = Some((x_cand.to_vec(), params_stage.to_vec(), rep.clone()));
            Ok(rep)
        };

        let mut x_cur = x0.clone();
        let mut initial_objective: Option<f64> = None;
        let mut cumulative_history: Vec<f64> = Vec::new();
        let mut stage_summaries: Vec<Value> = Vec::new();
        let mut total_iters: usize = 0;
        let mut final_reason = StopReason::MaxIterations;

        for (stage_idx, stage) in stages.iter().enumerate() {
            for p in &mut params {
                if let DesignParam::TopologyDensity { beta, .. } = p {
                    *beta = stage.beta;
                }
            }
            let params_stage = params.clone();

            match opt_spec.algorithm.as_str() {
                "lbfgs" | "ncg" => {
                    let min_opts = MinOptions {
                        gtol: opt_spec.gtol,
                        ftol: opt_spec.ftol,
                        max_it: stage.max_iters,
                        memory: 10,
                        c1: 1e-4,
                        c2: if opt_spec.algorithm == "ncg" { 0.2 } else { 0.9 },
                    };
                    // Scale EmeSegmentLength variables by their initial magnitude so L-BFGS / NCG initial step d0 = -g
                    // is O(1) in dimensionless length units while keeping all other parameter kinds at unit scale (1.0).
                    let mut var_scales = vec![1.0_f64; n_vars];
                    for slot in &slots {
                        if let ParamSlot::EmeSegmentLength { offset, .. } = slot {
                            var_scales[*offset] = x0[*offset].abs().max(1e-9);
                        }
                    }
                    let u_cur: Vec<f64> = x_cur.iter().zip(&var_scales).map(|(&x, &s)| x / s).collect();
                    let mut x_buf = vec![0.0_f64; n_vars];
                    let mut first_eval_in_stage = true;
                    let mut fg = |u_try: &[f64], g_out: &mut [f64]| -> f64 {
                        for k in 0..n_vars {
                            x_buf[k] = u_try[k] * var_scales[k];
                        }
                        match evaluate_candidate(&x_buf, &params_stage) {
                            Ok(rep) => {
                                if initial_objective.is_none() {
                                    initial_objective = Some(rep.objective);
                                }
                                first_eval_in_stage = false;
                                for (k, (&gk, &xk)) in rep.flat_gradient.iter().zip(&x_buf).enumerate() {
                                    let mut sg = sense_sign * gk * var_scales[k];
                                    // Project active box bounds so projected gradient norm vanishes cleanly at bounds.
                                    if (xk <= xmin[k] + 1e-14 * var_scales[k] && sg > 0.0)
                                        || (xk >= xmax[k] - 1e-14 * var_scales[k] && sg < 0.0)
                                    {
                                        sg = 0.0;
                                    }
                                    g_out[k] = sg;
                                }
                                sense_sign * rep.objective
                            }
                            Err(e) => {
                                if first_eval_in_stage {
                                    eval_error = Some(e);
                                }
                                g_out.fill(0.0);
                                f64::INFINITY
                            }
                        }
                    };
                    let res = if opt_spec.algorithm == "lbfgs" {
                        lbfgs(&mut fg, &u_cur, &min_opts)
                    } else {
                        nonlinear_cg(&mut fg, &u_cur, &min_opts)
                    };
                    if let Some(err) = eval_error.take() {
                        return Err(err);
                    }
                    x_cur = res.x.iter().zip(&var_scales).map(|(&u, &s)| u * s).collect();
                    total_iters += res.iterations;
                    final_reason = res.reason;
                    let phys_hist: Vec<f64> = res.history.iter().map(|&v| sense_sign * v).collect();
                    if stage_idx == 0 {
                        cumulative_history.extend_from_slice(&phys_hist);
                    } else if phys_hist.len() > 1 {
                        cumulative_history.extend_from_slice(&phys_hist[1..]);
                    } else {
                        cumulative_history.extend_from_slice(&phys_hist);
                    }
                    stage_summaries.push(json!({
                        "stage": stage_idx,
                        "beta": stage.beta,
                        "iterations": res.iterations,
                        "evaluations": res.evaluations,
                        "objective": sense_sign * res.f,
                        "grad_norm": res.grad_norm,
                        "stop_reason": stop_reason_name(res.reason),
                    }));
                }
                "mma" => {
                    let n_mma_constraints: usize = opt_spec
                        .constraints
                        .iter()
                        .map(|c| match c {
                            ConstraintSpec::VolumeFraction { min_fraction: Some(_), .. } => 2,
                            ConstraintSpec::VolumeFraction { min_fraction: None, .. } => 1,
                        })
                        .sum();
                    let mut mma_opts = MmaOptions::new(xmin.clone(), xmax.clone());
                    mma_opts.move_limit = opt_spec.move_limit;
                    mma_opts.max_it = stage.max_iters;
                    mma_opts.xtol = opt_spec.xtol;
                    let mut first_eval_in_stage = true;
                    // Normalize `(f0, g0)` by initial L1 gradient norm so `sum |g0_scaled| == 1 == sum |gc|`
                    // across arbitrary SI units (aF..GHz..MJ), keeping dual multipliers well inside `[0, penalty]`.
                    let mut stage_obj_scale: Option<f64> = None;
                    let mut mma_cb = |x_try: &[f64]| -> MmaEval {
                        match evaluate_candidate(x_try, &params_stage) {
                            Ok(rep) => {
                                if initial_objective.is_none() {
                                    initial_objective = Some(rep.objective);
                                }
                                let scale = *stage_obj_scale.get_or_insert_with(|| {
                                    let g_l1: f64 = rep.flat_gradient.iter().map(|g| g.abs()).sum();
                                    if g_l1 > 1e-30 {
                                        g_l1
                                    } else {
                                        rep.objective.abs().max(1.0)
                                    }
                                });
                                let inv_scale = 1.0 / scale;
                                first_eval_in_stage = false;
                                let g0: Vec<f64> =
                                    rep.flat_gradient.iter().map(|&g| sense_sign * g * inv_scale).collect();
                                let mut c_vals = Vec::new();
                                let mut gc_vecs = Vec::new();
                                for con in &opt_spec.constraints {
                                    match con {
                                        ConstraintSpec::VolumeFraction { param_index, max_fraction, min_fraction } => {
                                            let vf = rep.volume_fractions.get(param_index).copied().unwrap_or(0.0);
                                            let mut full_g = vec![0.0; n_vars];
                                            if let (
                                                Some(ParamSlot::Topology { offset, len, .. }),
                                                Some(sub_g),
                                            ) = (
                                                slots.iter().find(|s| matches!(s, ParamSlot::Topology { param_idx, .. } if *param_idx == *param_index)),
                                                rep.volume_fraction_gradients.get(param_index),
                                            ) {
                                                if sub_g.len() == *len {
                                                    full_g[*offset..*offset + *len].copy_from_slice(sub_g);
                                                }
                                            }
                                            c_vals.push(vf - *max_fraction);
                                            gc_vecs.push(full_g.clone());
                                            if let Some(lo) = min_fraction {
                                                c_vals.push(*lo - vf);
                                                gc_vecs.push(full_g.into_iter().map(|v| -v).collect());
                                            }
                                        }
                                    }
                                }
                                MmaEval { f0: sense_sign * rep.objective * inv_scale, g0, c: c_vals, gc: gc_vecs }
                            }
                            Err(e) => {
                                if first_eval_in_stage {
                                    eval_error = Some(e);
                                }
                                MmaEval {
                                    f0: 1e30,
                                    g0: vec![0.0; n_vars],
                                    c: vec![0.0; n_mma_constraints],
                                    gc: vec![vec![0.0; n_vars]; n_mma_constraints],
                                }
                            }
                        }
                    };
                    let res = mma(&mut mma_cb, &x_cur, &mma_opts);
                    if let Some(err) = eval_error.take() {
                        return Err(err);
                    }
                    let s_out = stage_obj_scale.unwrap_or(1.0);
                    x_cur = res.x;
                    total_iters += res.iterations;
                    final_reason = res.reason;
                    let phys_hist: Vec<f64> = res.history.iter().map(|&v| sense_sign * s_out * v).collect();
                    if stage_idx == 0 {
                        cumulative_history.extend_from_slice(&phys_hist);
                    } else if phys_hist.len() > 1 {
                        cumulative_history.extend_from_slice(&phys_hist[1..]);
                    } else {
                        cumulative_history.extend_from_slice(&phys_hist);
                    }
                    stage_summaries.push(json!({
                        "stage": stage_idx,
                        "beta": stage.beta,
                        "iterations": res.iterations,
                        "objective": sense_sign * s_out * res.f0,
                        "constraints": res.constraints,
                        "stop_reason": stop_reason_name(res.reason),
                    }));
                }
                other => {
                    return Err(Error::invalid(format!("unsupported algorithm '{other}'")));
                }
            }
        }

        // Ensure we have the exact sensitivity report at the final `x_cur` (cached hit if `x_cur` was last evaluated).
        let final_report = evaluate_candidate(&x_cur, &params)?;
        drop(evaluate_candidate);
        let init_obj = initial_objective.unwrap_or(final_report.objective);

        // Build final parameter summary map and optionally apply optimal scalar/CAD parameters to `self.project`.
        let mut final_params_map = serde_json::Map::new();
        let mut updated_spec = base_spec.clone();
        let mut has_spec_slots = false;
        let mut final_eme_overrides: Vec<(usize, f64)> = Vec::new();
        for slot in &slots {
            match slot {
                ParamSlot::ParamVar { name, offset, .. } => {
                    let v = x_cur[*offset];
                    final_params_map.insert(name.clone(), json!(v));
                    if apply_project {
                        geometry::set_param(&mut self.project.geometry, name, &json!(v))?;
                    }
                }
                ParamSlot::MaterialScalar { region, property, offset, .. } => {
                    let v = x_cur[*offset];
                    final_params_map.insert(format!("{region}.{property}"), json!(v));
                    set_spec_material_scalar(&mut updated_spec, region, property, v);
                    has_spec_slots = true;
                }
                ParamSlot::EmeSegmentLength { segment, offset, .. } => {
                    let v = x_cur[*offset];
                    final_params_map.insert(format!("eme_length[{segment}]"), json!(v));
                    final_eme_overrides.push((*segment, v));
                    has_spec_slots = true;
                }
                ParamSlot::Topology { param_idx, domain, offset, len, .. } => {
                    let vf = final_report.volume_fractions.get(param_idx).copied().unwrap_or(0.0);
                    let grey: f64 = final_report
                        .projected_densities
                        .get(param_idx)
                        .map(|rb| {
                            if rb.is_empty() {
                                0.0
                            } else {
                                rb.iter().map(|&r| 4.0 * r * (1.0 - r)).sum::<f64>() / rb.len() as f64
                            }
                        })
                        .unwrap_or(0.0);
                    final_params_map.insert(
                        format!("topology:{domain}"),
                        json!({
                            "cells": len,
                            "volume_fraction": vf,
                            "greyness": grey,
                            "min_raw": x_cur[*offset..*offset + *len].iter().copied().fold(1.0f64, f64::min),
                            "max_raw": x_cur[*offset..*offset + *len].iter().copied().fold(0.0f64, f64::max),
                        }),
                    );
                }
            }
        }
        if !final_eme_overrides.is_empty() {
            let _chain = crate::solve::install_chain_models(chain_models.clone());
            let ctx_tmp = crate::solve::build_solve_ctx(
                &self.registry,
                &base_model,
                p_phys,
                &updated_spec,
                &library,
                &progress,
                &done,
            )?;
            nyx_ems::eme_set_spec_segment_lengths(
                &mut updated_spec,
                &ctx_tmp,
                &base_eme_lengths,
                &final_eme_overrides,
            )?;
        }
        if apply_project && has_spec_slots {
            if let Some(sp) = self.project.spec_mut(&target) {
                *sp = updated_spec.clone();
            }
        }

        // Run a final forward solve to populate `self.results` (attaching `"topology_density"` CellField when applicable)
        // and record in `nyx-store` when bound to a workflow.
        let mut run_json = Value::Null;
        if apply_project {
            let final_model = if !has_topology && slots.iter().any(|s| matches!(s, ParamSlot::ParamVar { .. })) {
                self.model = None;
                self.ensure_physics_model(&target)?
            } else {
                base_model.clone()
            };
            let final_spec = self.project.spec(&target).cloned().unwrap_or(updated_spec);
            let _chain = crate::solve::install_chain_models(chain_models);
            let mut ctx_solve = crate::solve::build_solve_ctx(
                &self.registry,
                &final_model,
                p_phys,
                &final_spec,
                &library,
                &progress,
                &done,
            )?;
            // If only scalar/shape parameters exist, run `p_phys.solve` so fields reflect the optimized design;
            // if topology densities exist, run `p_phys.solve` and attach the projected `"topology_density"` cell field.
            if let Ok(mut solved) = p_phys.solve(final_model.mesh(), &final_spec, &mut ctx_solve) {
                crate::solve::name_domains(&mut solved, &final_model);
                solved.scalars.push(("opt_initial_objective".into(), Scalar::Num(init_obj)));
                solved.scalars.push(("opt_final_objective".into(), Scalar::Num(final_report.objective)));
                solved.scalars.push(("opt_iterations".into(), Scalar::Num(total_iters as f64)));

                if has_topology {
                    let simp = final_model.solver_mesh();
                    let n_cells = simp.n_cells();
                    let nloc = simp.dim + 1;
                    let mut cell_rho = vec![0.0f64; n_cells];
                    for slot in &slots {
                        if let ParamSlot::Topology { param_idx, filter, .. } = slot {
                            if let Some(rb) = final_report.projected_densities.get(param_idx) {
                                for (di, &ci) in filter.design_cells().iter().enumerate() {
                                    if let Some(&val) = rb.get(di) {
                                        cell_rho[ci] = val;
                                    }
                                }
                            }
                        }
                    }
                    let mut expanded = Vec::with_capacity(n_cells * nloc);
                    for &r in &cell_rho {
                        expanded.extend(std::iter::repeat_n(r, nloc));
                    }
                    let cf = CellField {
                        name: "topology_density".into(),
                        unit: "1".into(),
                        nloc,
                        re: expanded,
                        im: Vec::new(),
                    };
                    if solved.solutions.is_empty() {
                        solved.solutions.push(nyx_physics::Solution {
                            label: "optimized".into(),
                            cell_fields: vec![cf],
                            ..Default::default()
                        });
                    } else {
                        for sol in &mut solved.solutions {
                            sol.cell_fields.retain(|f| f.name != "topology_density");
                            sol.cell_fields.push(cf.clone());
                        }
                    }
                }

                let elapsed = t0.elapsed().as_secs_f64();
                let key = self.problem_key(&target).unwrap_or_default();
                let results_vec = vec![solved];
                if let Some(link) = self.link.as_mut() {
                    link.script = self.script.clone();
                    let log_lines = vec![format!(
                        "optimize_run ({}) finished in {total_iters} iterations: J = {:.6e} -> {:.6e}",
                        opt_spec.algorithm, init_obj, final_report.objective
                    )];
                    if let Ok(rec) = link.record_with_options(
                        &self.project,
                        RunOutcome::Solved(&results_vec),
                        elapsed,
                        &log_lines,
                        Some(&final_model.solver_mesh()),
                        true,
                        &[],
                    ) {
                        run_json = json!({ "id": rec.id, "version": rec.version, "status": rec.status });
                    }
                }
                let run_id = run_json.get("id").and_then(Value::as_str).map(String::from);
                if let Some(res_item) = results_vec.into_iter().next() {
                    self.results.insert(
                        target.clone(),
                        Stored {
                            result: Arc::new(res_item),
                            key,
                            job: 0,
                            at: Instant::now(),
                            model: final_model,
                            run: run_id,
                        },
                    );
                }
            }
        }

        let converged = final_reason.is_converged()
            || (opt_spec.algorithm == "mma"
                && !cumulative_history.is_empty()
                && ((opt_spec.sense == OptimizeSense::Minimize && final_report.objective <= init_obj)
                    || (opt_spec.sense == OptimizeSense::Maximize && final_report.objective >= init_obj)));

        let mut reply = json!({
            "physics": target,
            "algorithm": opt_spec.algorithm,
            "sense": opt_spec.sense,
            "converged": converged,
            "stop_reason": stop_reason_name(final_reason),
            "iterations": total_iters,
            "evaluations": total_evals,
            "initial_objective": init_obj,
            "final_objective": final_report.objective,
            "history": cumulative_history,
            "stages": stage_summaries,
            "final_params": Value::Object(final_params_map),
            "aux_scalars": final_report.aux_scalars,
            "volume_fractions": final_report.volume_fractions,
            "seconds": t0.elapsed().as_secs_f64(),
        });
        if !run_json.is_null() {
            reply["run"] = run_json;
        }
        if !final_report.projected_densities.is_empty() {
            reply["projected_densities"] = json!(final_report.projected_densities);
        }

        if let Some(cp) = checkpoint_path {
            if let Some(parent) = std::path::Path::new(&cp).parent().filter(|p| !p.as_os_str().is_empty()) {
                std::fs::create_dir_all(parent).map_err(|e| Error::invalid(format!("{}: {e}", parent.display())))?;
            }
            let text = serde_json::to_string_pretty(&nyx_physics::canonical_json(&reply))
                .map_err(|e| Error::invalid(format!("checkpoint serialization failed: {e}")))?;
            std::fs::write(&cp, text).map_err(|e| Error::invalid(format!("{cp}: {e}")))?;
            reply["checkpoint_path"] = json!(cp);
        }

        Ok(reply)
    }
}
