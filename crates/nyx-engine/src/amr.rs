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

/// True when `project` and `targets` request multi-pass Maxwell residual AMR on `ems::cavity`.
pub fn project_uses_ems_cavity_amr(project: &Project, targets: &[String]) -> bool {
    cfg!(feature = "ems")
        && project.mesh.sweep.is_none()
        && project.mesh.adapt_passes.unwrap_or(1) > 1
        && targets.iter().any(|t| t == "ems")
        && project.spec("ems").is_some_and(|s| s.study.kind == "cavity")
}

/// True when `project` and `targets` request multi-pass recovered-Hessian AMR on `sms`.
pub fn project_uses_sms_amr(project: &Project, targets: &[String]) -> bool {
    cfg!(feature = "sms")
        && project.mesh.sweep.is_none()
        && project.mesh.adapt_passes.unwrap_or(1) > 1
        && targets.iter().any(|t| t == "sms")
        && project.spec("sms").is_some()
}

/// True when `project` and `targets` activate multi-pass AMR for any compiled physics.
pub fn project_uses_amr(project: &Project, targets: &[String]) -> bool {
    project_uses_electro_amr(project, targets)
        || project_uses_ems_cavity_amr(project, targets)
        || project_uses_sms_amr(project, targets)
}

/// True when `project` carries any physics configured for multi-pass AMR on an unswept mesh.
pub fn project_has_amr(project: &Project) -> bool {
    let all: Vec<String> = project.physics.iter().map(|s| s.physics.clone()).collect();
    project_uses_amr(project, &all)
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
#[allow(dead_code)]
pub(crate) fn solve_dfs_with_amr(
    registry: &Registry,
    project: &Project,
    initial_model: Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<(Arc<MeshedModel>, SolveResult), Error> {
    solve_dfs_with_amr_opts(registry, project, initial_model, library, progress, done, project.mesh.retain_amr_data)
}

pub(crate) fn solve_dfs_with_amr_opts(
    registry: &Registry,
    project: &Project,
    initial_model: Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
    retain_amr_data: bool,
) -> Result<(Arc<MeshedModel>, SolveResult), Error> {
    let spec = project
        .spec("dfs")
        .ok_or_else(|| Error::invalid("the project has no problem for 'dfs'; add the physics first"))?;
    let Some(amr) = parse_electro_amr_params(spec) else {
        let mut r = crate::solve::solve_one_single(registry, project, &initial_model, "dfs", library, progress, done)?;
        if !retain_amr_data {
            ElectroAmrData::strip_from(&mut r);
        }
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
        if !retain_amr_data {
            ElectroAmrData::strip_from(&mut result);
        }

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

fn geom_mesh_opts(project: &Project) -> femcore::meshing::geomdrive::GeomMeshOpts {
    let mut opts = femcore::meshing::geomdrive::GeomMeshOpts::new(project.mesh.h);
    opts.max_dev = project.mesh.max_dev;
    opts.curved = project.mesh.curved;
    opts.gap_cells = project.mesh.gap_cells;
    opts.z_coarsen = project.mesh.z_coarsen;
    opts.max_cells = project.mesh.max_cells;
    opts.dz_min = project.mesh.dz_min;
    opts.max_z_slabs = project.mesh.max_z_slabs;
    opts.z_growth_ratio = project.mesh.z_growth_ratio;
    opts.max_steiner_multiplier = project.mesh.max_steiner_multiplier;
    opts.strict_mesh = project.mesh.strict_mesh;
    opts.singular_edges = project.mesh.singular_edges.clone();
    opts.quarter_point_xi = project.mesh.quarter_point_xi;
    opts.corner_grading = project.mesh.corner_grading;
    opts
}

fn max_cells_cap(project: &Project) -> usize {
    match (project.mesh.max_cells, project.mesh.adapt_max_cells) {
        (Some(a), Some(b)) => a.min(b),
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => 2_000_000,
    }
}

fn max_rel_result_delta(prev: &SolveResult, curr: &SolveResult) -> Option<f64> {
    let mut worst: Option<f64> = None;
    let consider = |worst: &mut Option<f64>, a: f64, b: f64| {
        if a.is_finite() && b.is_finite() {
            let scale = a.abs().max(b.abs());
            if scale > 1.0e-30 {
                let d = (b - a).abs() / scale;
                *worst = Some(worst.map_or(d, |w| w.max(d)));
            }
        }
    };
    for (sp, sc) in prev.solutions.iter().zip(&curr.solutions) {
        if let (Some(a), Some(b)) = (sp.value, sc.value) {
            consider(&mut worst, a, b);
        }
        if let (Some([ar, ai]), Some([br, bi])) = (sp.eigenvalue, sc.eigenvalue) {
            let mag = (ar * ar + ai * ai).sqrt().max((br * br + bi * bi).sqrt());
            if mag > 1.0e-30 {
                let diff = ((br - ar).powi(2) + (bi - ai).powi(2)).sqrt() / mag;
                worst = Some(worst.map_or(diff, |w| w.max(diff)));
            }
        }
    }
    if worst.is_none() {
        for (name, v) in &prev.scalars {
            if let Scalar::Num(a) = v {
                if let Some((_, Scalar::Num(b))) = curr.scalars.iter().find(|(n, _)| n == name) {
                    consider(&mut worst, *a, *b);
                }
            }
        }
    }
    worst
}

#[cfg(feature = "ems")]
struct CavityCurlMode<'a> {
    mesh: &'a femcore::mesh::SimplexMesh,
    ned: &'a nyx_ems::cavity::Ned3,
    x: Vec<femcore::linalg::csr::C>,
    order: usize,
}

#[cfg(feature = "ems")]
impl femcore::adapt::maxwell::CurlField for CavityCurlMode<'_> {
    fn sample(&self, cell: usize, lam: [f64; 4]) -> ([femcore::linalg::csr::C; 3], [femcore::linalg::csr::C; 3]) {
        use femcore::linalg::csr::C;
        use nyx_ems::cavity::basis3_keyed;
        let v = self.mesh.cell(cell);
        let p: [[f64; 3]; 4] = std::array::from_fn(|i| {
            let q = self.mesh.point(v[i] as usize);
            [q[0], q[1], q[2]]
        });
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        };
        let (e1, e2, e3) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
        let c23 = cross(e2, e3);
        let six_v = e1[0] * c23[0] + e1[1] * c23[1] + e1[2] * c23[2];
        let g1 = c23.map(|x| x / six_v);
        let g2 = cross(e3, e1).map(|x| x / six_v);
        let g3 = cross(e1, e2).map(|x| x / six_v);
        let g0 = [-(g1[0] + g2[0] + g3[0]), -(g1[1] + g2[1] + g3[1]), -(g1[2] + g2[2] + g3[2])];
        let gid = self.ned.key_pairs(self.mesh, cell);
        let (mut w, mut cu) = (vec![], vec![]);
        basis3_keyed(self.order, &lam, &[g0, g1, g2, g3], &gid, &mut w, &mut cu);
        let mut e = [C::new(0.0, 0.0); 3];
        let mut c = [C::new(0.0, 0.0); 3];
        for (s, &dof) in self.ned.cell_dofs[cell].iter().enumerate() {
            let coef = self.x[dof];
            for k in 0..3 {
                e[k] += coef * w[s][k];
                c[k] += coef * cu[s][k];
            }
        }
        (e, c)
    }
}

/// Runs multi-pass goal-oriented Maxwell residual AMR (`femcore::adapt::maxwell`) on 3D `ems::cavity` studies.
#[cfg(feature = "ems")]
pub(crate) fn solve_ems_cavity_with_amr(
    registry: &Registry,
    project: &Project,
    initial_model: Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<(Arc<MeshedModel>, SolveResult), Error> {
    use femcore::adapt::maxwell::{indicator_metrics, maxwell_indicators, EquidistributeOptions, MaxwellOptions};
    use femcore::meshing::adapt3d::{adapt_to_metrics3, geometry_floor3};
    use std::collections::HashSet;
    use std::f64::consts::PI;

    let passes = project.mesh.adapt_passes.unwrap_or(1);
    if passes <= 1
        || project.mesh.sweep.is_some()
        || !initial_model.mesh().is_pure_simplex()
        || initial_model.mesh().dim() != 3
    {
        let r = crate::solve::solve_one_single(registry, project, &initial_model, "ems", library, progress, done)?;
        return Ok((initial_model, r));
    }

    let p = crate::physics_for(registry, "ems")?;
    let spec = project
        .spec("ems")
        .ok_or_else(|| Error::invalid("the project has no problem for 'ems'; add the physics first"))?;
    nyx_physics::validate_spec(p, spec)?;

    let rel_tol = project.mesh.adapt_rel_tol.unwrap_or(1.0e-3);
    let growth = project.mesh.adapt_growth.unwrap_or(1.5).max(1.05);
    let grade = project.mesh.adapt_grade.unwrap_or(0.4);
    let max_cells_3d = max_cells_cap(project);
    let h_max = project.mesh.h;
    let geom_opts = geom_mesh_opts(project);

    let mut current_model = initial_model;
    let mut prev_result: Option<SolveResult> = None;

    for pass_idx in 1..=passes {
        if progress.is_cancelled() {
            return Err(Error::Cancelled);
        }
        p.validate(spec, current_model.mesh())?;
        let ctx = crate::solve::build_solve_ctx(registry, &current_model, p, spec, library, progress, done);
        let pre_audit = crate::material_warnings::preflight(p, spec, &ctx, project.strict_material_range)?;
        let problem = nyx_ems::decode(spec, &ctx)?;
        let sm_cur = current_model.solver_mesh();
        let n_cells = sm_cur.n_cells();

        let _installed = femcore::progress::install(progress);
        let first_line = progress.log().len();
        let _s = femcore::progress::stage(format!(
            "Maxwell cavity AMR pass {pass_idx}/{passes}: solving on {n_cells} tetrahedra"
        ));
        let asm = nyx_ems::assemble_emw_eigen(&sm_cur, &problem).map_err(Error::from_kernel)?;
        let solved = nyx_ems::solve_emw_modes(&asm).map_err(Error::from_kernel)?;
        let mut result = nyx_ems::emw_eigen_fields(&sm_cur, &problem, &asm, &solved).map_err(Error::from_kernel)?;
        nyx_physics::check_result_finite(&result)?;
        result.log.extend(progress.log().into_iter().skip(first_line));
        crate::solve::name_domains(&mut result, &current_model);
        crate::material_warnings::finalize(&mut result, p, spec, &ctx, project.strict_material_range, pre_audit)?;

        let max_rel_delta = prev_result.as_ref().and_then(|prev| max_rel_result_delta(prev, &result));
        let eps_r: Vec<f64> =
            asm.eps_phys_re.iter().map(|t| ((t[0][0] + t[1][1] + t[2][2]) / 3.0).max(1.0e-12)).collect();
        let mut pml_doms: HashSet<i32> = problem.materials.pml.keys().copied().collect();
        for (c, &is_pml) in asm.pml_cells.iter().enumerate() {
            if is_pml {
                pml_doms.insert(sm_cur.cell_labels[c]);
            }
        }
        let mut natural_labels: HashSet<i32> = problem.pmc.iter().copied().collect();
        if let Ok((_, sym_pmc)) = nyx_ems::symmetry_facet_flags(&sm_cur, problem.symmetry) {
            for (f, &is_pmc) in sym_pmc.iter().enumerate() {
                if is_pmc {
                    natural_labels.insert(sm_cur.facet_labels[f]);
                }
            }
        }

        let mut combined_eta2 = vec![0.0_f64; n_cells];
        let mut primary_total = 0.0_f64;
        for (m_pos, &(freq, vec_idx)) in solved.freqs.iter().take(4).enumerate() {
            let k0 = 2.0 * PI * freq.re.abs() / nyx_ems::C0;
            if !(k0 > 0.0 && k0.is_finite()) {
                continue;
            }
            let mode = CavityCurlMode {
                mesh: &sm_cur,
                ned: &asm.ned,
                x: asm.red.expand(&solved.vectors[vec_idx]),
                order: asm.order,
            };
            let skip_fn = |dom: i32| pml_doms.contains(&dom);
            let nat_fn = |lab: i32| natural_labels.contains(&lab);
            let opt =
                MaxwellOptions { k0_sq: k0 * k0, eps_r: &eps_r, mu_r: &asm.mur_phys, skip: &skip_fn, natural: &nat_fn };
            if let Ok(ind) = maxwell_indicators(&sm_cur, &mode, &opt) {
                if m_pos == 0 {
                    primary_total = ind.total;
                }
                let denom = (ind.total * ind.total).max(1.0e-300);
                for (acc, &e2) in combined_eta2.iter_mut().zip(&ind.eta2) {
                    *acc += e2 / denom;
                }
            }
        }

        let delta_str = max_rel_delta.map_or(String::new(), |d| format!(", max_rel_delta = {d:.3e}"));
        result.log.push(format!(
            "Maxwell AMR pass {pass_idx}/{passes}: cells = {n_cells}, dofs = {}, total_eta = {primary_total:.3e}{delta_str}",
            result.dofs
        ));

        let has_signal = combined_eta2.iter().any(|&v| v > 0.0 && v.is_finite());
        if !has_signal || max_rel_delta.is_some_and(|d| d <= rel_tol) || pass_idx == passes {
            return Ok((current_model, result));
        }

        let target_3d = ((n_cells as f64) * growth).min(max_cells_3d as f64).max((n_cells as f64) + 8.0);
        let a3 = Adaptive { scale: 1.0, grade, max_cells: Some(max_cells_3d) };
        let g_in = MeshedGeom {
            mesh: (*sm_cur).clone(),
            cells: None,
            labels: current_model.labels().clone(),
            eval: current_model.evaluated().clone(),
            warnings: vec![],
        };
        let floor = geometry_floor3(&g_in, &geom_opts, &a3).map_err(Error::Invalid)?;
        let eq_opt = EquidistributeOptions { order: asm.order, h_max, h_min_fraction: 1.0e-4, target_cells: target_3d };
        let met3 = indicator_metrics(&sm_cur, &combined_eta2, &eq_opt, Some(&floor)).map_err(Error::Invalid)?;
        let (g_next, _) = adapt_to_metrics3(g_in, &met3, &a3).map_err(Error::Invalid)?;
        current_model = Arc::new(MeshedModel::new(g_next.mesh, g_next.eval, g_next.labels, false)?.with_amr_adapted());
        prev_result = Some(result);
    }

    unreachable!("Maxwell AMR loop returns on final pass")
}

#[cfg(not(feature = "ems"))]
pub(crate) fn solve_ems_cavity_with_amr(
    registry: &Registry,
    project: &Project,
    initial_model: Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<(Arc<MeshedModel>, SolveResult), Error> {
    let r = crate::solve::solve_one_single(registry, project, &initial_model, "ems", library, progress, done)?;
    Ok((initial_model, r))
}

/// Runs multi-pass recovered-Hessian AMR (`femcore::adapt::{estimate, estimate3}`) over `PlotFields` for `sms` studies.
pub(crate) fn solve_sms_with_amr(
    registry: &Registry,
    project: &Project,
    initial_model: Arc<MeshedModel>,
    library: &Arc<femcore::material::Library>,
    progress: &Arc<Progress>,
    done: &[(String, SolveResult)],
) -> Result<(Arc<MeshedModel>, SolveResult), Error> {
    use femcore::adapt::estimate::{estimator_metrics as estimator_metrics_2d, Estimator};
    use femcore::adapt::estimate3::estimator_metrics as estimator_metrics_3d;
    use femcore::meshing::adapt3d::{adapt_to_metrics3, geometry_floor3};

    let passes = project.mesh.adapt_passes.unwrap_or(1);
    let sms_phys = crate::physics_for(registry, "sms")?;
    let has_upstream = !done.is_empty()
        && sms_phys.imports().iter().any(|imp| {
            done.iter().any(|(c, _)| registry.get(c).is_some_and(|q| q.exports().iter().any(|e| e.role == imp.role)))
        });
    if passes <= 1 || has_upstream || project.mesh.sweep.is_some() || !initial_model.mesh().is_pure_simplex() {
        let r = crate::solve::solve_one_single(registry, project, &initial_model, "sms", library, progress, done)?;
        return Ok((initial_model, r));
    }

    let rel_tol = project.mesh.adapt_rel_tol.unwrap_or(1.0e-3);
    let growth = project.mesh.adapt_growth.unwrap_or(1.5).max(1.05);
    let grade = project.mesh.adapt_grade.unwrap_or(0.4);
    let max_cells_cap = max_cells_cap(project);
    let h_max = project.mesh.h;
    let geom_opts = geom_mesh_opts(project);

    let mut current_model = initial_model;
    let mut prev_result: Option<SolveResult> = None;

    for pass_idx in 1..=passes {
        if progress.is_cancelled() {
            return Err(Error::Cancelled);
        }
        let n_cells = current_model.mesh().n_cells();
        let mut result = {
            let _s = femcore::progress::stage(format!(
                "structural mechanics AMR pass {pass_idx}/{passes}: solving on {n_cells} cells"
            ));
            crate::solve::solve_one_single(registry, project, &current_model, "sms", library, progress, done)?
        };
        let max_rel_delta = prev_result.as_ref().and_then(|prev| max_rel_result_delta(prev, &result));
        let delta_str = max_rel_delta.map_or(String::new(), |d| format!(", max_rel_delta = {d:.3e}"));
        result
            .log
            .push(format!("SMS AMR pass {pass_idx}/{passes}: cells = {n_cells}, dofs = {}{delta_str}", result.dofs));

        if max_rel_delta.is_some_and(|d| d <= rel_tol) || pass_idx == passes {
            return Ok((current_model, result));
        }

        let est = if let Some(sol) = result.solutions.last() {
            let has_field = |n: &str| sol.fields.iter().any(|f| f.name == n);
            if has_field("von Mises") {
                Estimator { fields: vec!["von Mises".into()], ..Estimator::default() }
            } else {
                let comps: Vec<String> =
                    ["ux", "uy", "uz"].into_iter().filter(|n| has_field(n)).map(String::from).collect();
                if comps.len() >= 2 {
                    Estimator { fields: comps, combine: true, ..Estimator::default() }
                } else {
                    Estimator::default()
                }
            }
        } else {
            Estimator::default()
        };

        let sm_cur = current_model.solver_mesh();
        let g_in = MeshedGeom {
            mesh: (*sm_cur).clone(),
            cells: None,
            labels: current_model.labels().clone(),
            eval: current_model.evaluated().clone(),
            warnings: vec![],
        };
        let target_cells = ((n_cells as f64) * growth).min(max_cells_cap as f64).max((n_cells as f64) + 8.0);
        let adapt_cfg = Adaptive { scale: 1.0, grade, max_cells: Some(max_cells_cap) };

        let next_geom = if sm_cur.dim == 2 {
            let met2 = match estimator_metrics_2d(&g_in, &result, &est, &geom_opts, target_cells) {
                Ok(m) => m,
                Err(_) => return Ok((current_model, result)),
            };
            adapt_to_metrics(g_in, &met2, &adapt_cfg).map_err(Error::Invalid)?.0
        } else {
            let floor = geometry_floor3(&g_in, &geom_opts, &adapt_cfg).map_err(Error::Invalid)?;
            let met3 = match estimator_metrics_3d(&sm_cur, &result, &est, h_max, target_cells, Some(&floor)) {
                Ok(m) => m,
                Err(_) => return Ok((current_model, result)),
            };
            adapt_to_metrics3(g_in, &met3, &adapt_cfg).map_err(Error::Invalid)?.0
        };
        current_model =
            Arc::new(MeshedModel::new(next_geom.mesh, next_geom.eval, next_geom.labels, false)?.with_amr_adapted());
        prev_result = Some(result);
    }

    unreachable!("SMS AMR loop returns on final pass")
}
