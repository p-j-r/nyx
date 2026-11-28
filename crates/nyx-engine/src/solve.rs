// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Solving a project: the physics in an order that satisfies the one-way couplings (or an explicit
//! opt-in Anderson-accelerated two-way coupled loop via [`solve_coupled`], ADR 139), each checked
//! against its descriptors and the mesh, then solved with the fields it imports from the physics
//! before it. Everything goes through `Physics`, so nothing here names a physics crate directly.
//!
//! Per-physics meshing (ADR 135): a physics with its own entry in `Project::mesh_by_physics` is meshed with those settings (cached by the
//! hash of geometry, settings and sizing, in the chain and process-wide while a model is alive) and solved on that mesh. A field it imports
//! from a physics solved on another mesh is transferred by canonical region name (`femcore::transfer`): `Ctx::inputs` gets the field sampled
//! at the dofs of the importing mesh and `Ctx::transferred_inputs` the field on its own mesh. Both one-way chains and Anderson coupled loops
//! (`solve_coupled`) operate seamlessly across independent per-physics meshes, while keeping the zero-cost path when a single mesh is shared.

use crate::meshing::{FemcoreMesher, MeshedModel, Mesher};
use crate::physics_for;
use crate::project::{CoupledLoopConfig, Project};
use femcore::linalg::snes::{anderson, FixedPointOptions};
use femcore::progress::Progress;
use nyx_physics::{
    validate_spec, Ctx, Error, Field, ProblemSpec, Registry, Scalar, Sizing, SolveResult, TransferredField,
};
use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, OnceLock, Weak};

/// Output of [`solve_coupled`]: the converged electromagnetic/conduction source result, the stationary
/// thermal result, the number of outer coupled evaluations performed, and the relative $L_\infty$
/// temperature residual history.
#[derive(Clone, Debug)]
pub struct CoupledSolveOutput {
    pub source_result: SolveResult,
    pub heat_result: SolveResult,
    pub iterations: usize,
    pub residual_history: Vec<f64>,
    pub rel_residual: f64,
}

/// The field a physics exports under `role`, taken from its last solution.
fn exported<'r>(registry: &Registry, physics: &str, result: &'r SolveResult, role: &str) -> Option<&'r Field> {
    let p = registry.get(physics)?;
    let link = p.exports().into_iter().find(|l| l.role == role)?;
    let name = link.field.split_once('.').map(|(_, n)| n).unwrap_or(&link.field).to_string();
    result.solutions.last()?.field(&name)
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

// ------------------------------------------------------------------------------------------------------------ per-physics meshes

thread_local! {
    /// The mesh each physics solved so far in the running chain was solved on (set by the chain around each solve; read by
    /// `build_solve_ctx`, which also serves the AMR loops of `amr.rs`). Empty outside a chain: imports then take the same-mesh path.
    static CHAIN_MODELS: RefCell<BTreeMap<String, Arc<MeshedModel>>> = const { RefCell::new(BTreeMap::new()) };
}

/// Restores the previous chain models when dropped.
pub(crate) struct ChainModelsGuard(Option<BTreeMap<String, Arc<MeshedModel>>>);

impl Drop for ChainModelsGuard {
    fn drop(&mut self) {
        if let Some(prev) = self.0.take() {
            CHAIN_MODELS.with(|c| *c.borrow_mut() = prev);
        }
    }
}

pub(crate) fn install_chain_models(models: BTreeMap<String, Arc<MeshedModel>>) -> ChainModelsGuard {
    ChainModelsGuard(Some(CHAIN_MODELS.with(|c| std::mem::replace(&mut *c.borrow_mut(), models))))
}

fn chain_model(physics: &str) -> Option<Arc<MeshedModel>> {
    CHAIN_MODELS.with(|c| c.borrow().get(physics).cloned())
}

/// Meshes made for physics of their own, by key, while some result or chain still holds them.
fn mesh_cache() -> &'static Mutex<BTreeMap<String, Weak<MeshedModel>>> {
    static CACHE: OnceLock<Mutex<BTreeMap<String, Weak<MeshedModel>>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(BTreeMap::new()))
}

/// The key of a per-physics mesh: the hash of the canonical geometry, mesh settings and sizing (what a mesh depends on).
pub(crate) fn physics_mesh_key(project: &Project, physics: &str, sizing: Option<&Sizing>) -> String {
    let v = serde_json::json!({ "geometry": project.geometry, "mesh": project.mesh_for(physics), "sizing": sizing });
    nyx_store::util::hash_bytes(nyx_physics::canonical_json(&v).to_string().as_bytes())
}

/// The sizing a physics with its own mesh asks for (only its own waves), when its settings size from the physics.
fn own_sizing(
    registry: &Registry,
    project: &Project,
    physics: &str,
    library: &Arc<femcore::material::Library>,
) -> Option<Sizing> {
    if !project.mesh_for(physics).physics_sizing {
        return None;
    }
    let p = registry.get(physics)?;
    let s = p.sizing(project.spec(physics)?, Some(library));
    (!s.wavelengths.is_empty()).then_some(s)
}

/// The mesh `physics` is solved on: `base` unless it has settings of its own, else its own mesh (from the chain's cache, the process cache,
/// or made here and cached).
fn model_for_physics(
    registry: &Registry,
    project: &Project,
    physics: &str,
    base: &Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    mesher: &dyn Mesher,
    local: &mut BTreeMap<String, Arc<MeshedModel>>,
) -> Result<Arc<MeshedModel>, Error> {
    if !project.has_own_mesh(physics) {
        return Ok(base.clone());
    }
    let settings = project.mesh_for(physics);
    let sizing = own_sizing(registry, project, physics, library);
    let key = physics_mesh_key(project, physics, sizing.as_ref());
    if let Some(m) = local.get(&key) {
        return Ok(m.clone());
    }
    let cached = mesh_cache().lock().unwrap_or_else(|e| e.into_inner()).get(&key).and_then(Weak::upgrade);
    if let Some(m) = cached {
        let _s =
            femcore::progress::stage(format!("mesh of '{physics}': the same mesh was built by an earlier request"));
        local.insert(key, m.clone());
        return Ok(m);
    }
    let _stage = femcore::progress::stage(format!(
        "meshing the geometry for '{physics}' (its own element size {} m)",
        settings.h
    ));
    let mut model = mesher.mesh(&project.geometry, settings, sizing.as_ref())?;
    if let Some(active) = &settings.active_domains {
        model = model.restrict_to_active(active)?;
    }
    let m = Arc::new(model);
    {
        let mut cache = mesh_cache().lock().unwrap_or_else(|e| e.into_inner());
        cache.retain(|_, w| w.strong_count() > 0);
        cache.insert(key.clone(), Arc::downgrade(&m));
    }
    local.insert(key, m.clone());
    Ok(m)
}

/// The project as `physics` sees it: its own mesh settings in place of the project's (so the checks and the AMR loops read them).
fn project_view<'p>(project: &'p Project, physics: &str) -> Cow<'p, Project> {
    if project.has_own_mesh(physics) {
        let mut p = project.clone();
        p.mesh = project.mesh_for(physics).clone();
        Cow::Owned(p)
    } else {
        Cow::Borrowed(project)
    }
}

// ------------------------------------------------------------------------------------------------------------ chains

/// What a chain hands back: the project mesh after the chain (adapted when AMR ran on it), the results in solve order, the mesh each
/// physics was solved on, and the error that stopped the chain.
pub struct ChainOutcome {
    pub model: Arc<MeshedModel>,
    pub results: Vec<SolveResult>,
    /// by physics code; a physics without settings of its own maps to the project mesh it was solved on
    pub models: BTreeMap<String, Arc<MeshedModel>>,
    pub failed: Option<Error>,
}

struct ExtraCoupledInput<'a> {
    role: &'a str,
    owner: &'a str,
    field: Field,
    src_model: &'a Arc<MeshedModel>,
    order: usize,
}

#[allow(clippy::too_many_arguments)]
fn solve_spec_with_inputs(
    registry: &Registry,
    project: &Project,
    model: &MeshedModel,
    spec: &ProblemSpec,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
    extra_inputs: &[ExtraCoupledInput<'_>],
) -> Result<SolveResult, Error> {
    if progress.is_cancelled() {
        return Err(Error::Cancelled);
    }
    let p = physics_for(registry, &spec.physics)?;
    validate_spec(p, spec)?;
    p.validate(spec, model.mesh())?;
    let mut ctx = build_solve_ctx(registry, model, p, spec, library, progress, done)?;
    for ex in extra_inputs {
        if ex.role == "temperature" {
            let (mut sum, mut count) = (0.0_f64, 0usize);
            for &v in &ex.field.values {
                if v.is_finite() && v > 0.0 {
                    sum += v;
                    count += 1;
                }
            }
            if count > 0 {
                ctx.temperature = Some(sum / count as f64);
            }
        }
        let same = std::ptr::eq(Arc::as_ptr(ex.src_model), model as *const MeshedModel)
            || femcore::transfer::same_simplex_mesh(&ex.src_model.solver_mesh(), &model.solver_mesh());
        if same {
            ctx.inputs.insert(ex.role.to_string(), ex.field.values.clone());
        } else {
            let tf = TransferredField {
                physics: ex.owner.to_string(),
                field: ex.field.name.clone(),
                mesh: ex.src_model.mesh().clone(),
                discretisation: Some(ex.src_model.solver_mesh()),
                domain_names: ex.src_model.domain_name_map(),
                boundary_names: ex.src_model.boundary_name_map(),
                order: ex.order.clamp(1, 2),
                components: ex.field.components.len().max(1),
                values: ex.field.values.clone(),
                cell_fields: Vec::new(),
            };
            let target_order = (spec.order as usize).clamp(1, 2);
            let sampled = tf.sample_dofs(&model.solver_mesh(), &model.domain_name_map(), target_order)?;
            if sampled.missing > 0 && sampled.found.iter().all(|&ok| !ok) {
                return Err(Error::invalid(format!(
                    "the {} imported from '{}' ({}) cannot be transferred to the mesh of '{}': none of its nodes lie in regions '{}' was solved on (check `active_domains` of the two meshes)",
                    ex.role,
                    ex.owner,
                    ex.field.name,
                    p.id(),
                    ex.owner
                )));
            }
            ctx.inputs.insert(ex.role.to_string(), sampled.values);
            ctx.transferred_inputs.insert(ex.role.to_string(), tf);
        }
    }
    let pre_audit = crate::material_warnings::preflight(p, spec, &ctx, project.strict_material_range)?;
    let mut result = p.solve(model.mesh(), spec, &mut ctx)?;
    name_domains(&mut result, model);
    crate::material_warnings::finalize(&mut result, p, spec, &ctx, project.strict_material_range, pre_audit)?;
    Ok(result)
}

/// Runs a self-consistent two-way Anderson-accelerated coupled loop (`femcore::linalg::snes::anderson`,
/// default depth `m = 5`) between `source_spec` (exporting `"heat_source"`, importing `"temperature"`)
/// and `heat_spec` (stationary thermal problem exporting `"temperature"`, importing `"heat_source"`),
/// terminating when `||T^{k+1} - T^k||_inf / max(1.0, ||T^k||_inf) <= cfg.rtol`.
#[allow(clippy::too_many_arguments)]
pub fn solve_coupled(
    registry: &Registry,
    project: &Project,
    model: &MeshedModel,
    source_spec: &ProblemSpec,
    heat_spec: &ProblemSpec,
    cfg: &CoupledLoopConfig,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done_before: &[(String, SolveResult)],
) -> Result<CoupledSolveOutput, Error> {
    let shared = Arc::new(model.shallow_clone());
    solve_coupled_models(
        registry,
        project,
        &shared,
        &shared,
        source_spec,
        heat_spec,
        cfg,
        library,
        progress,
        done_before,
    )
}

/// Like [`solve_coupled`], allowing `source_spec` and `heat_spec` to execute on distinct per-physics
/// [`MeshedModel`] instances (`ADR 135` + `ADR 139`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_coupled_models(
    registry: &Registry,
    project: &Project,
    src_model: &Arc<MeshedModel>,
    heat_model: &Arc<MeshedModel>,
    source_spec: &ProblemSpec,
    heat_spec: &ProblemSpec,
    cfg: &CoupledLoopConfig,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done_before: &[(String, SolveResult)],
) -> Result<CoupledSolveOutput, Error> {
    if cfg.max_iters < 1 {
        return Err(Error::invalid("coupled_loop.max_iters must be >= 1"));
    }
    if !(cfg.rtol > 0.0 && cfg.rtol.is_finite()) {
        return Err(Error::invalid("coupled_loop.rtol must be a positive finite number"));
    }
    if !(cfg.source_scale >= 0.0 && cfg.source_scale.is_finite()) {
        return Err(Error::invalid("coupled_loop.source_scale must be a non-negative finite number"));
    }

    let src_view = project_view(project, &source_spec.physics);
    let heat_view = project_view(project, &heat_spec.physics);

    // 1. Baseline thermal solve (zero imported heat source) to establish physical initial state T^(0).
    let heat_init =
        solve_spec_with_inputs(registry, &heat_view, heat_model, heat_spec, library, progress, done_before, &[])?;
    let t0_field = exported(registry, &heat_spec.physics, &heat_init, "temperature").ok_or_else(|| {
        Error::invalid(format!(
            "coupled thermal problem '{}' ('{}') did not export a 'temperature' field",
            heat_spec.physics, heat_spec.study.kind
        ))
    })?;
    let t0 = &t0_field.values;
    if t0.is_empty() || t0.iter().any(|v| !v.is_finite()) {
        return Err(Error::NonFinite("initial thermal field T^(0) is empty or non-finite".into()));
    }
    let t_name = t0_field.name.clone();
    let t_unit = t0_field.unit.clone();
    let heat_order = heat_init.order as usize;

    let captured_err: RefCell<Option<Error>> = RefCell::new(None);
    let converged_pair: RefCell<Option<(SolveResult, SolveResult, usize, f64)>> = RefCell::new(None);
    let last_pair: RefCell<Option<(SolveResult, SolveResult, usize, f64)>> = RefCell::new(None);
    let residual_history: RefCell<Vec<f64>> = RefCell::new(Vec::new());
    let eval_counter: RefCell<usize> = RefCell::new(0);

    // Smooth C^1 positive-Kelvin diffeomorphism `psi: R -> (0, inf)` that is the exact identity `psi(u) == u`
    // for all `u >= t_floor` (preserving exact linear/multisecant Anderson convergence in Kelvin `T`) while
    // smoothly mapping `u < t_floor` to `t_floor * exp((u - t_floor) / t_floor) > 0` without flat-clamp kinks.
    let t0_min = t0.iter().copied().fold(f64::INFINITY, f64::min).max(1e-3);
    let t_floor = (0.5 * t0_min).clamp(1e-3, 50.0);
    let psi_to_kelvin = |u: f64| -> f64 {
        if u >= t_floor {
            u
        } else {
            t_floor * ((u - t_floor) / t_floor).max(-40.0).exp()
        }
    };
    let psi_from_kelvin = |tk: f64| -> f64 {
        let tc = tk.max(1e-12);
        if tc >= t_floor {
            tc
        } else {
            t_floor * (1.0 + (tc / t_floor).ln())
        }
    };
    let u0: Vec<f64> = t0.iter().map(|&tk| psi_from_kelvin(tk)).collect();

    let mut g_map = |u: &[f64], gu: &mut [f64]| {
        if captured_err.borrow().is_some() {
            gu.fill(f64::NAN);
            return;
        }
        if u.iter().any(|v| !v.is_finite()) {
            *captured_err.borrow_mut() =
                Some(Error::NonFinite("non-finite temperature vector encountered in Anderson step".into()));
            gu.fill(f64::NAN);
            return;
        }
        let step_num = {
            let mut c = eval_counter.borrow_mut();
            *c += 1;
            *c
        };
        let t_cur: Vec<f64> = u.iter().map(|&ui| psi_to_kelvin(ui)).collect();

        let temp_extra = ExtraCoupledInput {
            role: "temperature",
            owner: &heat_spec.physics,
            field: Field::scalar(&t_name, &t_unit, t_cur.clone()),
            src_model: heat_model,
            order: heat_order,
        };
        let src_res = match solve_spec_with_inputs(
            registry,
            &src_view,
            src_model,
            source_spec,
            library,
            progress,
            done_before,
            &[temp_extra],
        ) {
            Ok(r) => r,
            Err(e) => {
                *captured_err.borrow_mut() = Some(e);
                gu.fill(f64::NAN);
                return;
            }
        };

        let (q_name, q_order, mut q_imp) = match exported(registry, &source_spec.physics, &src_res, "heat_source") {
            Some(qf) => (qf.name.clone(), src_res.order as usize, qf.values.clone()),
            None => ("Q_opt".to_string(), 1usize, vec![0.0; src_model.solver_mesh().n_points()]),
        };
        if (cfg.source_scale - 1.0).abs() > 0.0 {
            for qv in &mut q_imp {
                *qv *= cfg.source_scale;
            }
        }
        let heat_extra = ExtraCoupledInput {
            role: "heat_source",
            owner: &source_spec.physics,
            field: Field::scalar(q_name, "W/m^3", q_imp),
            src_model,
            order: q_order,
        };

        let heat_res = match solve_spec_with_inputs(
            registry,
            &heat_view,
            heat_model,
            heat_spec,
            library,
            progress,
            done_before,
            &[heat_extra],
        ) {
            Ok(r) => r,
            Err(e) => {
                *captured_err.borrow_mut() = Some(e);
                gu.fill(f64::NAN);
                return;
            }
        };

        let Some(t_next_field) = exported(registry, &heat_spec.physics, &heat_res, "temperature") else {
            *captured_err.borrow_mut() =
                Some(Error::Solver("thermal step inside coupled loop produced no 'temperature' field".into()));
            gu.fill(f64::NAN);
            return;
        };
        let t_next = &t_next_field.values;
        if t_next.len() != u.len() || t_next.iter().any(|v| !v.is_finite()) {
            *captured_err.borrow_mut() =
                Some(Error::NonFinite("thermal step inside coupled loop produced non-finite temperature".into()));
            gu.fill(f64::NAN);
            return;
        }

        let diff_inf = t_next.iter().zip(&t_cur).map(|(&tn, &tk)| (tn - tk).abs()).fold(0.0_f64, f64::max);
        let norm_t_inf = t_cur.iter().map(|&tk| tk.abs()).fold(1.0_f64, f64::max);
        let rel_inf = diff_inf / norm_t_inf;
        residual_history.borrow_mut().push(rel_inf);

        if rel_inf <= cfg.rtol {
            *converged_pair.borrow_mut() = Some((src_res, heat_res, step_num, rel_inf));
            // Zero fixed-point residual `fk = gu - u = 0` so `snes::anderson` exits cleanly without an extra solve.
            gu.copy_from_slice(u);
        } else {
            for (gui, &tn) in gu.iter_mut().zip(t_next) {
                *gui = psi_from_kelvin(tn);
            }
            *last_pair.borrow_mut() = Some((src_res, heat_res, step_num, rel_inf));
        }
    };

    let fp_opts = FixedPointOptions {
        rtol: 0.0,
        atol: 0.0,
        max_it: cfg.max_iters.saturating_sub(1),
        beta: 1.0,
        depth: cfg.anderson_m,
    };
    let _anderson_out = anderson(&mut g_map, &u0, &fp_opts);

    if let Some(err) = captured_err.into_inner() {
        return Err(err);
    }
    let hist = residual_history.into_inner();
    let Some((mut src_res, mut heat_res, iterations, rel_residual)) = converged_pair.into_inner() else {
        let last_rel = last_pair.into_inner().map_or(f64::NAN, |(_, _, _, r)| r);
        return Err(Error::Solver(format!(
            "coupled Anderson loop did not converge within {} iterations (relative L_inf temperature residual {:.3e} > rtol {:.3e})",
            cfg.max_iters, last_rel, cfg.rtol
        )));
    };

    let summary_line = format!(
        "coupled Anderson loop (m = {}): converged in {} iteration(s), ||dT||_inf / max(1, ||T||_inf) = {:.3e} <= {:.3e}",
        cfg.anderson_m, iterations, rel_residual, cfg.rtol
    );
    src_res.scalars.push(("coupled_iterations".into(), Scalar::Num(iterations as f64)));
    src_res.scalars.push(("coupled_rel_residual".into(), Scalar::Num(rel_residual)));
    src_res.log.push(summary_line.clone());

    heat_res.scalars.push(("coupled_iterations".into(), Scalar::Num(iterations as f64)));
    heat_res.scalars.push(("coupled_rel_residual".into(), Scalar::Num(rel_residual)));
    heat_res.log.push(summary_line);

    // When running intra-physics electrothermal (`dfs::conduction <-> dfs::heat_stationary`), expose
    // `P_joule (W)` and matching-order `V` / `Q_joule` directly on `heat_res` so single-physics `Project`
    // callers can inspect both electrical and thermal quantities on `"dfs"`.
    if source_spec.physics == heat_spec.physics {
        if let Some((_, Scalar::Num(pj))) = src_res.scalars.iter().find(|(k, _)| k == "P_joule (W)") {
            heat_res.scalars.push(("P_joule (W)".into(), Scalar::Num(*pj * cfg.source_scale)));
        }
        if let (Some(src_sol), Some(heat_sol)) = (src_res.solutions.last(), heat_res.solutions.last_mut()) {
            let n_expect = heat_sol.field("T").map_or(0, |f| f.values.len());
            for fname in ["V", "Q_joule"] {
                if heat_sol.field(fname).is_none() {
                    if let Some(f) = src_sol.field(fname).filter(|f| f.values.len() == n_expect) {
                        let vals = if fname == "Q_joule" && (cfg.source_scale - 1.0).abs() > 0.0 {
                            f.values.iter().map(|&v| v * cfg.source_scale).collect()
                        } else {
                            f.values.clone()
                        };
                        heat_sol.fields.push(Field::scalar(fname, &f.unit, vals));
                    }
                }
            }
        }
    }

    Ok(CoupledSolveOutput {
        source_result: src_res,
        heat_result: heat_res,
        iterations,
        residual_history: hist,
        rel_residual,
    })
}

pub(crate) fn project_solve_order(registry: &Registry, project: &Project, ids: &[&str]) -> Result<Vec<String>, Error> {
    let items: Vec<(&str, Option<&str>)> =
        ids.iter().map(|&id| (id, project.spec(id).map(|s| s.study.kind.as_str()))).collect();
    registry.solve_order_for_studies(&items)
}

fn coupled_mode<'a>(
    project: &'a Project,
    targets: &[String],
) -> Option<(&'a CoupledLoopConfig, &'a ProblemSpec, &'a ProblemSpec, bool)> {
    let cfg = project.coupled_loop.as_ref()?;
    let has_dfs = targets.iter().any(|t| t == "dfs");
    if !has_dfs {
        return None;
    }
    let heat_spec = project.spec("dfs")?;
    if heat_spec.study.kind != "heat_stationary" {
        return None;
    }
    if let Some(ref aux) = cfg.auxiliary {
        return Some((cfg, aux, heat_spec, true));
    }
    if targets.iter().any(|t| t == "ems") {
        let ems_spec = project.spec("ems")?;
        return Some((cfg, ems_spec, heat_spec, false));
    }
    None
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
    let base = Arc::new(model.shallow_clone());
    let out =
        run_chain_core(registry, project, base, targets, library, progress, retain_amr_data, false, &FemcoreMesher);
    (out.results, out.failed)
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
    model: Arc<MeshedModel>,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
    retain_amr_data: bool,
) -> (Arc<MeshedModel>, Vec<SolveResult>, Option<Error>) {
    let out =
        run_chain_core(registry, project, model, targets, library, progress, retain_amr_data, true, &FemcoreMesher);
    (out.model, out.results, out.failed)
}

/// The chain with AMR where configured, per-physics meshes made with `mesher`, and the mesh of every result.
#[allow(clippy::too_many_arguments)]
pub fn run_chain_models(
    registry: &Registry,
    project: &Project,
    model: Arc<MeshedModel>,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
    retain_amr_data: bool,
    mesher: &dyn Mesher,
) -> ChainOutcome {
    run_chain_core(registry, project, model, targets, library, progress, retain_amr_data, true, mesher)
}

#[allow(clippy::too_many_arguments)]
fn run_chain_core(
    registry: &Registry,
    project: &Project,
    mut model: Arc<MeshedModel>,
    targets: &[String],
    library: &Arc<femcore::material::Library>,
    progress: Arc<Progress>,
    retain_amr_data: bool,
    amr: bool,
    mesher: &dyn Mesher,
) -> ChainOutcome {
    let mut models: BTreeMap<String, Arc<MeshedModel>> = BTreeMap::new();
    let mut local_meshes: BTreeMap<String, Arc<MeshedModel>> = BTreeMap::new();
    let mut done: Vec<(String, SolveResult)> = vec![];

    let order = if let Some((cfg, src_spec, heat_spec, is_aux)) = coupled_mode(project, targets) {
        let heat_model = if project.has_own_mesh(&heat_spec.physics) {
            let _installed = femcore::progress::install(&progress);
            match model_for_physics(registry, project, &heat_spec.physics, &model, library, mesher, &mut local_meshes) {
                Ok(m) => m,
                Err(e) => return ChainOutcome { model, results: vec![], models, failed: Some(e) },
            }
        } else {
            model.clone()
        };
        let src_model = if is_aux && src_spec.physics == heat_spec.physics {
            heat_model.clone()
        } else if project.has_own_mesh(&src_spec.physics) {
            let _installed = femcore::progress::install(&progress);
            match model_for_physics(registry, project, &src_spec.physics, &model, library, mesher, &mut local_meshes) {
                Ok(m) => m,
                Err(e) => return ChainOutcome { model, results: vec![], models, failed: Some(e) },
            }
        } else {
            model.clone()
        };

        let coupled_res = {
            let _chain = install_chain_models(models.clone());
            solve_coupled_models(
                registry,
                project,
                &src_model,
                &heat_model,
                src_spec,
                heat_spec,
                cfg,
                library,
                &progress,
                &done,
            )
        };
        match coupled_res {
            Ok(mut out) => {
                if !retain_amr_data {
                    nyx_physics::ElectroAmrData::strip_from(&mut out.source_result);
                    nyx_physics::ElectroAmrData::strip_from(&mut out.heat_result);
                }
                if !is_aux {
                    models.insert(src_spec.physics.clone(), src_model);
                    done.push((src_spec.physics.clone(), out.source_result));
                }
                models.insert(heat_spec.physics.clone(), heat_model);
                done.push((heat_spec.physics.clone(), out.heat_result));
            }
            Err(e) => return ChainOutcome { model, results: vec![], models, failed: Some(e) },
        }

        let rem: Vec<&str> =
            targets.iter().map(String::as_str).filter(|&id| !done.iter().any(|(d, _)| d == id)).collect();
        if rem.is_empty() {
            Vec::new()
        } else {
            match project_solve_order(registry, project, &rem) {
                Ok(o) => o,
                Err(e) => {
                    return ChainOutcome {
                        model,
                        results: done.into_iter().map(|(_, r)| r).collect(),
                        models,
                        failed: Some(e),
                    };
                }
            }
        }
    } else {
        let ids: Vec<&str> = targets.iter().map(String::as_str).collect();
        match project_solve_order(registry, project, &ids) {
            Ok(o) => o,
            Err(e) => return ChainOutcome { model, results: vec![], models, failed: Some(e) },
        }
    };

    let mut failed = None;
    for id in order {
        if progress.is_cancelled() {
            failed = Some(Error::Cancelled);
            break;
        }
        let own = project.has_own_mesh(&id);
        let target = if own {
            let _installed = femcore::progress::install(&progress);
            match model_for_physics(registry, project, &id, &model, library, mesher, &mut local_meshes) {
                Ok(m) => m,
                Err(e) => {
                    failed = Some(e);
                    break;
                }
            }
        } else {
            model.clone()
        };
        let view = project_view(project, &id);
        let view: &Project = &view;
        let step = {
            let _chain = install_chain_models(models.clone());
            let use_electro_amr = amr && crate::amr::project_uses_electro_amr(view, targets);
            let use_ems_amr = amr && crate::amr::project_uses_ems_cavity_amr(view, targets);
            let use_sms_amr = amr && crate::amr::project_uses_sms_amr(view, targets);
            if id == "dfs" && use_electro_amr {
                crate::amr::solve_dfs_with_amr_opts(
                    registry,
                    view,
                    target.clone(),
                    library,
                    &progress,
                    &done,
                    retain_amr_data,
                )
            } else if id == "ems" && use_ems_amr {
                crate::amr::solve_ems_cavity_with_amr(registry, view, target.clone(), library, &progress, &done)
            } else if id == "sms" && use_sms_amr {
                crate::amr::solve_sms_with_amr(registry, view, target.clone(), library, &progress, &done)
            } else {
                solve_one_single(registry, view, &target, &id, library, &progress, &done).map(|mut r| {
                    if !retain_amr_data {
                        nyx_physics::ElectroAmrData::strip_from(&mut r);
                    }
                    (target.clone(), r)
                })
            }
        };
        match step {
            Ok((next_model, result)) => {
                if own {
                    models.insert(id.clone(), next_model);
                } else {
                    model = next_model;
                    models.insert(id.clone(), model.clone());
                }
                done.push((id, result));
            }
            Err(e) => {
                failed = Some(e);
                break;
            }
        }
    }
    ChainOutcome { model, results: done.into_iter().map(|(_, r)| r).collect(), models, failed }
}

/// Builds a [`Ctx`] with imported fields from previously solved physics in `done`. An import solved on another mesh than `model` (per-physics
/// meshing) is transferred by canonical region name: `inputs` gets it sampled at the dofs of `model`, `transferred_inputs` the field itself.
pub(crate) fn build_solve_ctx<'a>(
    registry: &Registry,
    model: &'a MeshedModel,
    p: &dyn nyx_physics::Physics,
    spec: &nyx_physics::ProblemSpec,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<Ctx<'a>, Error> {
    let mut ctx = Ctx::new(model);
    ctx.progress = progress.clone();
    ctx.library = Some(library.clone());
    ctx.wavelength = crate::materials::study_wavelength(p, spec);
    ctx.frequency = crate::materials::study_frequency(p, spec);
    for imp in p.imports() {
        let Some((owner, result, field)) =
            done.iter().find_map(|(o, r)| exported(registry, o, r, &imp.role).map(|f| (o, r, f)))
        else {
            continue;
        };
        let values = field.values.clone();
        if imp.role == "temperature" && ctx.temperature.is_none() {
            let (mut sum, mut count) = (0.0_f64, 0usize);
            for &v in &values {
                if v.is_finite() && v > 0.0 {
                    sum += v;
                    count += 1;
                }
            }
            if count > 0 {
                ctx.temperature = Some(sum / count as f64);
            }
        }
        // the zero-cost path: the exporter was solved on this very mesh (or the chain does not say otherwise)
        let other = chain_model(owner).filter(|sm| {
            !std::ptr::eq(Arc::as_ptr(sm), model as *const MeshedModel)
                && !femcore::transfer::same_simplex_mesh(&sm.solver_mesh(), &model.solver_mesh())
        });
        let Some(sm) = other else {
            ctx.inputs.insert(imp.role.clone(), values);
            continue;
        };
        let tf = TransferredField {
            physics: owner.clone(),
            field: field.name.clone(),
            mesh: sm.mesh().clone(),
            discretisation: Some(sm.solver_mesh()),
            domain_names: sm.domain_name_map(),
            boundary_names: sm.boundary_name_map(),
            order: result.order as usize,
            components: field.components.len().max(1),
            values,
            cell_fields: result.solutions.last().map(|s| s.cell_fields.clone()).unwrap_or_default(),
        };
        let target_order = (spec.order as usize).clamp(1, 2);
        let sampled = tf.sample_dofs(&model.solver_mesh(), &model.domain_name_map(), target_order)?;
        let allow_partial = imp.role == "heat_source" || p.id() == "ems" || p.id() == "dfs";
        if sampled.missing > 0 && (!allow_partial || sampled.found.iter().all(|&ok| !ok)) {
            return Err(Error::invalid(format!(
                "the {} imported from '{owner}' ({}) cannot be transferred to the mesh of '{}': {} of its nodes lie in regions '{owner}' was not solved on (check `active_domains` of the two meshes)",
                imp.role,
                field.name,
                p.id(),
                sampled.missing
            )));
        }
        {
            let _installed = femcore::progress::install(progress);
            femcore::progress::note(format!(
                "{}: '{}' of '{owner}' transferred from its own mesh by region name ({} nodes snapped to their region)",
                imp.role, field.name, sampled.snapped
            ));
        }
        ctx.inputs.insert(imp.role.clone(), sampled.values);
        ctx.transferred_inputs.insert(imp.role.clone(), tf);
    }
    Ok(ctx)
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
    let _ = physics_for(registry, id)?;
    let spec = project
        .spec(id)
        .ok_or_else(|| Error::invalid(format!("the project has no problem for '{id}'; add the physics first")))?;
    let mesh_settings = project.mesh_for(id);
    if id == "dfs"
        && spec.study.kind == "electrostatic"
        && mesh_settings.singular_edges.as_deref() == Some("quarter_point")
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
    solve_spec_with_inputs(registry, project, model, spec, library, progress, done, &[])
}
