// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-pass goal-oriented adaptive mesh refinement (`AF-1`) for electrostatic studies in `nyx-engine`.
//!
//! Orchestrates the outer solve -> estimate -> adapt loop when `dfs` `electrostatic` has `adapt_passes > 1`,
//! supporting column-projected 2D-overlay + `z`-extrusion (`sheet_2d_z`, default on 3D `sheet3d` models),
//! unstructured sheet-preserving 3D adaptation (`full_3d`), and 2D cross-section metric adaptation, while
//! keeping physics crates decoupled from `nyx-engine`.

use crate::meshing::{build_mesh_request, MeshedModel};
use crate::project::Project;
use femcore::adapt::electro::{
    electro_metrics_3d, electro_overlay_metrics_2d, project_indicators_to_overlay_2d, ElectroMetric2dOptions,
    ElectroMetric3dOptions,
};
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::sheet3d::{build_sized_overlay_mesh_2d, extrude_overlay_2d_with_plan};
use femcore::meshing::sheet3d_amr::{
    adapt_overlay_2d, adapt_sheet_full_3d_with_plan, extract_overlay_sheet_features, sheet_rim_vertices_3d,
};
use femcore::meshing::{adapt_to_metrics, Adaptive};
use femcore::progress::Progress;
use nyx_physics::{
    AdaptPassRecord, CapacitanceMatrixData, ElectroAmrData, Error, ProblemSpec, Registry, Scalar, SolveResult, TlsMeta,
    CAPACITANCE_KEY,
};
use serde_json::Value;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ResolvedAmrMode {
    Sheet2dZ,
    Full3d,
    Mesh2d,
}

struct AmrParams {
    passes: usize,
    rel_tol: f64,
    min_converged_passes: usize,
    growth: f64,
    mode: String,
}

fn parse_electro_amr_params(spec: &ProblemSpec) -> Option<AmrParams> {
    if spec.physics != "dfs" || spec.study.kind != "electrostatic" {
        return None;
    }
    let obj = spec.study.params.as_object()?;
    let passes = obj.get("adapt_passes").and_then(Value::as_i64).unwrap_or(1);
    if passes <= 1 {
        return None;
    }
    if obj.get("formulation").and_then(Value::as_str) == Some("planar_bem") {
        return None;
    }
    let rel_tol = obj.get("adapt_rel_tol").and_then(Value::as_f64).unwrap_or(1.0e-3);
    let min_converged_passes =
        obj.get("adapt_min_converged_passes").and_then(Value::as_i64).unwrap_or(1).max(1) as usize;
    let growth = obj.get("adapt_growth").and_then(Value::as_f64).unwrap_or(1.35);
    let mode = obj.get("adapt_mode").and_then(Value::as_str).unwrap_or("auto").to_string();
    Some(AmrParams { passes: passes as usize, rel_tol, min_converged_passes, growth, mode })
}

/// True when `project` and `targets` request multi-pass electrostatic AMR (`adapt_passes > 1`) on an unswept mesh.
pub fn project_uses_electro_amr(project: &Project, targets: &[String]) -> bool {
    if project.mesh.sweep.is_some() || !targets.iter().any(|t| t == "dfs") {
        return false;
    }
    project.spec("dfs").and_then(parse_electro_amr_params).is_some()
}

/// True when `project` carries a `dfs` `electrostatic` study configured with `adapt_passes > 1`.
pub fn project_has_electro_amr(project: &Project) -> bool {
    project.mesh.sweep.is_none() && project.spec("dfs").and_then(parse_electro_amr_params).is_some()
}

fn max_rel_capacitance_delta(prev: &SolveResult, curr: &SolveResult) -> Option<f64> {
    if let (Some(cp), Some(cc)) = (CapacitanceMatrixData::of(prev), CapacitanceMatrixData::of(curr)) {
        let n = cc.maxwell_matrix_f.len();
        if cp.maxwell_matrix_f.len() == n && n > 0 {
            let mut c_max = 0.0_f64;
            for i in 0..n {
                for j in 0..n {
                    c_max = c_max.max(cc.maxwell_matrix_f[i][j].abs()).max(cp.maxwell_matrix_f[i][j].abs());
                }
            }
            if c_max <= 0.0 {
                return Some(0.0);
            }
            let thresh = 1.0e-3 * c_max;
            let mut worst = 0.0_f64;
            for i in 0..n {
                for j in 0..n {
                    let a = cp.maxwell_matrix_f[i][j];
                    let b = cc.maxwell_matrix_f[i][j];
                    let denom = a.abs().max(b.abs());
                    if denom >= thresh {
                        worst = worst.max((b - a).abs() / denom);
                    }
                }
            }
            return Some(worst);
        }
    }
    let scalar_num = |r: &SolveResult, prefix: &str| -> Option<f64> {
        r.scalars.iter().find_map(|(k, v)| match v {
            Scalar::Num(x) if k.starts_with(prefix) && x.is_finite() => Some(*x),
            _ => None,
        })
    };
    let (a, b) = scalar_num(prev, "capacitance C")
        .zip(scalar_num(curr, "capacitance C"))
        .or_else(|| scalar_num(prev, "stored energy W_e").zip(scalar_num(curr, "stored energy W_e")))?;
    let denom = a.abs().max(b.abs()).max(1.0e-300);
    Some((b - a).abs() / denom)
}

fn write_adapt_history(result: &mut SolveResult, history: &[AdaptPassRecord]) {
    if let Some(mut cap) = CapacitanceMatrixData::of(result) {
        cap.adapt_history = history.to_vec();
        let encoded = cap.scalar();
        if let Some(slot) = result.scalars.iter_mut().find(|(k, _)| k == CAPACITANCE_KEY) {
            *slot = encoded;
        } else {
            result.scalars.push(encoded);
        }
    }
}

/// Runs the multi-pass goal-oriented AMR loop for `dfs` `electrostatic` and returns `(final_model, final_result)`.
pub(crate) fn solve_dfs_with_amr(
    registry: &Registry,
    project: &Project,
    initial_model: Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<(Arc<MeshedModel>, SolveResult), Error> {
    let spec = project
        .spec("dfs")
        .ok_or_else(|| Error::invalid("the project has no problem for 'dfs'; add the physics first"))?;
    let Some(amr) = parse_electro_amr_params(spec) else {
        let mut r = crate::solve::solve_one_single(registry, project, &initial_model, "dfs", library, progress, done)?;
        ElectroAmrData::strip_from(&mut r);
        return Ok((initial_model, r));
    };

    let req = build_mesh_request(&project.mesh, None)?;
    let plan = req.plan();
    let grade = project.mesh.adapt_grade.unwrap_or(0.4);
    let max_cells_3d = match (project.mesh.max_cells, project.mesh.adapt_max_cells) {
        (Some(a), Some(b)) => a.min(b),
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => 2_000_000,
    };
    let h_max = project.mesh.h;
    let h_min_fraction = 1.0e-5;
    let order = (spec.order.max(1) as usize).min(2);

    let dim = initial_model.mesh().dim();
    let has_stack_overlay =
        dim == 3 && project.mesh.sweep.is_none() && initial_model.evaluated().complex().stack.is_some();
    let has_sheets = has_stack_overlay && initial_model.evaluated().complex().has_sheets();
    let resolved_mode = if dim == 2 {
        ResolvedAmrMode::Mesh2d
    } else {
        match amr.mode.as_str() {
            "full_3d" => ResolvedAmrMode::Full3d,
            "sheet_2d_z" => {
                if has_stack_overlay {
                    ResolvedAmrMode::Sheet2dZ
                } else {
                    ResolvedAmrMode::Full3d
                }
            }
            _ => {
                if has_sheets {
                    ResolvedAmrMode::Sheet2dZ
                } else {
                    ResolvedAmrMode::Full3d
                }
            }
        }
    };

    let mut m2_opt = if resolved_mode == ResolvedAmrMode::Sheet2dZ {
        Some(build_sized_overlay_mesh_2d(initial_model.evaluated(), &plan, None).map_err(Error::Invalid)?)
    } else {
        None
    };
    let sheet_feats = if resolved_mode == ResolvedAmrMode::Sheet2dZ {
        extract_overlay_sheet_features(initial_model.evaluated().complex())
    } else {
        None
    };

    let mut current_model = initial_model;
    let mut prev_result: Option<SolveResult> = None;
    let mut history: Vec<AdaptPassRecord> = Vec::with_capacity(amr.passes);
    let mut converged_streak = 0usize;

    for pass_idx in 1..=amr.passes {
        if progress.is_cancelled() {
            return Err(Error::Cancelled);
        }
        let mut result = {
            let _s = femcore::progress::stage(format!(
                "electrostatic AMR pass {pass_idx}/{}: solving on {} cells",
                amr.passes,
                current_model.mesh().n_cells()
            ));
            crate::solve::solve_one_single(registry, project, &current_model, "dfs", library, progress, done)?
        };
        let amr_ind = ElectroAmrData::of(&result);
        ElectroAmrData::strip_from(&mut result);

        // If `formulation == "auto"` selected `planar_bem` or returned no volumetric indicators, return single-pass result cleanly.
        let Some(ind) = amr_ind else {
            return Ok((current_model, result));
        };

        let n_cells = current_model.mesh().n_cells();
        let max_rel_delta = prev_result.as_ref().and_then(|prev| max_rel_capacitance_delta(prev, &result));
        let rec =
            AdaptPassRecord { pass: pass_idx, n_cells, n_dofs: result.dofs, max_rel_delta, total_eta: ind.total_eta };
        history.push(rec);
        write_adapt_history(&mut result, &history);

        let delta_str = match max_rel_delta {
            Some(d) => format!(", max_rel_delta = {d:.3e}"),
            None => String::new(),
        };
        result.log.push(format!(
            "AMR pass {pass_idx}/{}: cells = {n_cells}, dofs = {}, total_eta = {:.3e}{delta_str}",
            amr.passes, result.dofs, ind.total_eta
        ));

        if max_rel_delta.is_some_and(|d| d <= amr.rel_tol) {
            converged_streak += 1;
        } else {
            converged_streak = 0;
        }
        if converged_streak >= amr.min_converged_passes || pass_idx == amr.passes {
            return Ok((current_model, result));
        }

        // Adapt mesh for pass `pass_idx + 1`
        let sm_cur = current_model.solver_mesh();
        let _s = femcore::progress::stage(format!(
            "electrostatic AMR pass {pass_idx}/{}: adapting mesh ({:?})",
            amr.passes, resolved_mode
        ));

        match resolved_mode {
            ResolvedAmrMode::Sheet2dZ => {
                let m2 = m2_opt.as_ref().expect("Sheet2dZ initializes m2");
                let plane_axes = sheet_feats.as_ref().map_or([0, 1], |f| f.plane_axes);
                let eta2_2d =
                    project_indicators_to_overlay_2d(&sm_cur, m2, &ind.eta2, plane_axes).map_err(Error::Invalid)?;
                let target_2d = ((m2.n_cells() as f64) * amr.growth).max(m2.n_cells() as f64 + 4.0);
                let mut mopt = ElectroMetric2dOptions::new(order, h_max, h_min_fraction, target_2d);
                mopt.grade = grade;
                mopt.max_aspect = 20.0;
                mopt.rim_anisotropy = 4.0;
                if let Some(f) = &sheet_feats {
                    mopt.rim_segments = &f.rim_segments;
                    mopt.corners = &f.corners;
                }
                let met2 = electro_overlay_metrics_2d(m2, &eta2_2d, &mopt, None).map_err(Error::Invalid)?;
                let slabs_per_tri = (n_cells as f64 / (m2.n_cells().max(1) as f64)).max(3.0);
                let cap_2d = ((max_cells_3d as f64) / slabs_per_tri).ceil().max(16.0) as usize;
                let a2 = Adaptive { scale: 1.0, grade, max_cells: Some(cap_2d) };
                let st = current_model
                    .evaluated()
                    .complex()
                    .stack
                    .as_ref()
                    .ok_or_else(|| Error::invalid("Sheet2dZ requires a 3D stack geometry"))?;
                let (m2_next, _) = adapt_overlay_2d(&st.overlay, m2, &met2, &a2).map_err(Error::Invalid)?;
                let ev_next = femcore::meshing::geomdrive::evaluate(&project.geometry).map_err(Error::Invalid)?;
                let g_next =
                    extrude_overlay_2d_with_plan(ev_next, &m2_next, &plan, None, None).map_err(Error::Invalid)?;
                m2_opt = Some(m2_next);
                current_model =
                    Arc::new(MeshedModel::new(g_next.mesh, g_next.eval, g_next.labels, false)?.with_amr_adapted());
            }
            ResolvedAmrMode::Full3d => {
                let sheet_labs = TlsMeta::of(&result).map(|t| t.sheet_labels).unwrap_or_default();
                let rim_v = sheet_rim_vertices_3d(&sm_cur, &sheet_labs);
                let target_3d = ((n_cells as f64) * amr.growth).min(max_cells_3d as f64).max(8.0);
                let mut opt3 = ElectroMetric3dOptions::new(order, h_max, h_min_fraction, target_3d);
                if rim_v.iter().any(|&b| b) {
                    opt3.rim_vertices = Some(&rim_v);
                }
                let met3 = electro_metrics_3d(&sm_cur, &ind.eta2, &opt3, None).map_err(Error::Invalid)?;
                let a3 = Adaptive { scale: 1.0, grade, max_cells: Some(max_cells_3d) };
                let g_in = MeshedGeom {
                    mesh: (*sm_cur).clone(),
                    cells: None,
                    labels: current_model.labels().clone(),
                    eval: current_model.evaluated().clone(),
                    warnings: vec![],
                };
                let (g_next, _) =
                    adapt_sheet_full_3d_with_plan(g_in, &met3, &a3, &plan, None).map_err(Error::Invalid)?;
                current_model =
                    Arc::new(MeshedModel::new(g_next.mesh, g_next.eval, g_next.labels, false)?.with_amr_adapted());
            }
            ResolvedAmrMode::Mesh2d => {
                let target_2d = ((n_cells as f64) * amr.growth).min(max_cells_3d as f64).max(8.0);
                let mut mopt = ElectroMetric2dOptions::new(order, h_max, h_min_fraction, target_2d);
                mopt.grade = grade;
                let met2 = electro_overlay_metrics_2d(&sm_cur, &ind.eta2, &mopt, None).map_err(Error::Invalid)?;
                let a2 = Adaptive { scale: 1.0, grade, max_cells: Some(max_cells_3d) };
                let g_in = MeshedGeom {
                    mesh: (*sm_cur).clone(),
                    cells: None,
                    labels: current_model.labels().clone(),
                    eval: current_model.evaluated().clone(),
                    warnings: vec![],
                };
                let (g_next, _) = adapt_to_metrics(g_in, &met2, &a2).map_err(Error::Invalid)?;
                current_model =
                    Arc::new(MeshedModel::new(g_next.mesh, g_next.eval, g_next.labels, false)?.with_amr_adapted());
            }
        }

        prev_result = Some(result);
    }

    unreachable!("AMR loop returns on final pass")
}
