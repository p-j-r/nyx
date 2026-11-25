//! The one way to mesh a `geom` model: a `MeshPlan` says how fine and with which mesher, `mesh_model` returns a `MeshedGeom` with the
//! names, labels and exact geometry of the model, whichever mesher made the cells.
//!
//! * `Mesher::Simple`: chords of every geom edge, constrained Delaunay triangulation, Ruppert refinement, lattice points in finer regions
//!   (`geomdrive::mesh_geom_2d`). Fast, predictable, what the app has always used. Right for most 2D cross-sections.
//! * `Mesher::Adaptive`: the same start, then the metric adaptation of `adapt` towards a size field built from the model alone
//!   (region sizes, `Size` rules, absorbing layers, curvature of the edges), graded. Vertices that move or are created on a boundary or an
//!   interface land on the *exact* geom edge (`GeomOracle2`), not on the chord polygon, and the vertices of the geometry never move. Labels
//!   of cells and facets are carried through unchanged, so names, selectors and boundary conditions mean the same on either mesh.
//!
//! `adapt_to_metrics` is the hook for an error estimator: give it a metric per vertex of an existing meshed model (from a recovered
//! Hessian, say) and it adapts the same way. Nothing else of the interface changes.

pub mod adapt3d;
pub mod curve3;
pub mod geomdrive;
pub mod geomdrive3;
pub mod mesher;
pub mod mesher3d;
pub mod physics_size;
pub mod sizing;
pub mod pml_geom;
pub mod pml_geom3;
pub mod pml3;
pub mod pml_box;
pub mod pml_cyl;
pub mod pml_sphere;
pub mod refine3d;
pub mod request;
pub mod sheet3d;
pub mod sheet3d_amr;
pub mod sheet3d_curve;
pub mod singularity;
pub mod validate;

use crate::meshing::geomdrive::{
    evaluate, mesh_geom_2d_with_skip, size_map_with_skip, GeomMeshOpts, MeshedGeom,
};
use crate::mesh::{Mesh, SimplexMesh};
pub use crate::meshing::pml_geom::{
    structured_pml, structured_pml_native_with_curved, structured_pml_with_curved, PmlMode, PmlReport,
};
use crate::adapt::drive::{adapt2_with, AdaptOptions, AdaptReport};
use crate::adapt::field::{BackgroundMetric2, MetricField};
use crate::adapt::geom::Geometry2;
use crate::adapt::metric::{intersect, isotropic, Sym};
use geomdrive3::Sizes3;
use physics_size::PhysicsSizing;
pub use sizing::Sizing;
use crate::adapt::topo::{Sig, Work2};
use crate::adapt::work::{simplex_from_adapt, simplex_to_adapt};
use geom::diff::{diff, ChangeKind};
use geom::eval::{Evaluated, Evaluator, Tier};
use geom::labels::LabelMap;
use geom::model::{FeatureList, PmlGrowMode};
use geom::select::EntityKind;
use std::collections::HashMap;

/// Settings of the adaptive mesher.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Adaptive {
    /// Factor on every size the model asks for: 1 keeps them, 0.5 makes the mesh about four times as fine. The start mesh is made at the
    /// model's own sizes, so a factor below one is what adaptation adds.
    pub scale: f64,
    /// Largest growth of the element size per unit distance (see `metric::grade`). Adapting to an ungraded field is unreliable; keep it.
    pub grade: f64,
    /// Stop refining beyond this many cells (the mesh stays valid, `AdaptReport::capped` says so).
    pub max_cells: Option<usize>,
}

impl Default for Adaptive {
    fn default() -> Self {
        Adaptive { scale: 1.0, grade: 0.4, max_cells: Some(2_000_000) }
    }
}

/// Which mesher makes the cells.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mesher {
    Simple,
    Adaptive(Adaptive),
}

/// How a model is meshed: the sizes (as for `mesh_geom_2d`) and the mesher.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshPlan {
    pub opts: GeomMeshOpts,
    pub mesher: Mesher,
    /// how the absorbing layers of the model are meshed
    pub pml: PmlMode,
}

impl MeshPlan {
    /// The simple mesher at element size `h` (physics-controlled `PmlMode::Auto` default).
    pub fn simple(h: f64) -> Self {
        MeshPlan { opts: GeomMeshOpts::new(h), mesher: Mesher::Simple, pml: PmlMode::Auto }
    }

    /// The adaptive mesher: start at element size `h`, adapt to `scale * h` (and to the curvature of the edges at that scale).
    pub fn adaptive(h: f64, scale: f64) -> Self {
        MeshPlan { opts: GeomMeshOpts::new(h), mesher: Mesher::Adaptive(Adaptive { scale, ..Adaptive::default() }), pml: PmlMode::Auto }
    }
}

impl From<GeomMeshOpts> for MeshPlan {
    fn from(opts: GeomMeshOpts) -> Self {
        MeshPlan { opts, mesher: Mesher::Simple, pml: PmlMode::Auto }
    }
}

/// Return tuple for zero-eager-split native mixed-cell meshing entry points.
pub type NativeMeshTuple = (
    Mesh,
    LabelMap,
    std::sync::Arc<Evaluated>,
    Option<AdaptReport>,
    crate::meshing::request::MeshReport,
    Vec<String>,
);

/// Internal raw tuple before `MeshReport` summarization.
pub type NativeRawTuple = (
    Mesh,
    LabelMap,
    std::sync::Arc<Evaluated>,
    Option<AdaptReport>,
    Option<PmlReport>,
    Vec<String>,
);

type NativeSizingRawTuple = (
    Mesh,
    LabelMap,
    std::sync::Arc<Evaluated>,
    Option<AdaptReport>,
    Option<PmlReport>,
    Vec<String>,
    PhysicsSizing,
);

fn check_outward_inline_guard(ev: &Evaluated, pml: PmlMode) -> Result<(), String> {
    if pml == PmlMode::Inline && ev.pml_regions().iter().any(|p| p.mode == PmlGrowMode::Outward) {
        return Err(
            "outward absorbing layers have no inline CAD cells; use PmlMode::Auto or PmlMode::Structured"
                .into(),
        );
    }
    Ok(())
}

fn pml_skip_mask(ev: &Evaluated, use_structured: bool) -> Vec<bool> {
    if !use_structured {
        return Vec::new();
    }
    let mut s = vec![false; ev.regions().len()];
    for p in ev.pml_regions() {
        s[p.region as usize] = true;
    }
    s
}

/// Meshes a model.
pub fn mesh_model(list: &FeatureList, plan: &MeshPlan) -> Result<MeshedGeom, String> {
    mesh_evaluated(evaluate(list)?, plan)
}

/// Meshes an evaluated model.
pub fn mesh_evaluated(ev: Evaluated, plan: &MeshPlan) -> Result<MeshedGeom, String> {
    Ok(mesh_evaluated_report(ev, plan)?.0)
}

/// Like `mesh_evaluated`, with the report of the adaptation when there was one.
pub fn mesh_evaluated_report(ev: Evaluated, plan: &MeshPlan) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    let (g, adapt, _) = mesh_evaluated_with_pml_report(ev, plan)?;
    Ok((g, adapt))
}

pub(crate) fn mesh_evaluated_with_pml_report(
    ev: Evaluated,
    plan: &MeshPlan,
) -> Result<(MeshedGeom, Option<AdaptReport>, Option<PmlReport>), String> {
    check_outward_inline_guard(&ev, plan.pml)?;
    if ev.list.dim == 3 {
        return mesh_3d_plan(ev, plan, None);
    }
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    let skip = pml_skip_mask(&ev, use_structured);
    let start = mesh_geom_2d_with_skip(ev, &plan.opts, &skip)?;
    let (g, report) = match plan.mesher {
        Mesher::Simple => (start, None),
        Mesher::Adaptive(a) => {
            let metrics = a_priori_metrics(&start, &plan.opts, &a)?;
            let (g, r) = adapt_to_metrics(start, &metrics, &a)?;
            (g, Some(r))
        }
    };
    if use_structured {
        let (g, pml_rep) = structured_pml_with_curved(g, plan.opts.curved)?;
        Ok((g, report, Some(pml_rep)))
    } else {
        Ok((g, report, None))
    }
}

/// Meshes an evaluated model directly into an uncut [`Mesh`] (`Quad`/`Prism`/`Hex` preserved without
/// eager simplex splitting), returning `(Mesh, LabelMap, Arc<Evaluated>, Option<AdaptReport>, MeshReport, Vec<String>)`.
pub fn mesh_evaluated_native_report(ev: Evaluated, plan: &MeshPlan) -> Result<NativeMeshTuple, String> {
    let (mesh, labels, eval, adapt, pml_rep, warnings) = mesh_evaluated_native_inner(ev, plan)?;
    let report = crate::meshing::request::summarise_native_mesh(
        &mesh,
        &labels,
        &eval,
        plan,
        adapt.as_ref(),
        pml_rep.as_ref(),
        &[],
        &warnings,
    );
    Ok((mesh, labels, eval, adapt, report, warnings))
}

pub(crate) fn mesh_evaluated_native_inner(
    ev: Evaluated,
    plan: &MeshPlan,
) -> Result<NativeRawTuple, String> {
    check_outward_inline_guard(&ev, plan.pml)?;
    if ev.list.dim == 3 {
        return mesh_3d_plan_native(ev, plan, None);
    }
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    let skip = pml_skip_mask(&ev, use_structured);
    let start = mesh_geom_2d_with_skip(ev, &plan.opts, &skip)?;
    let (g, adapt) = match plan.mesher {
        Mesher::Simple => (start, None),
        Mesher::Adaptive(a) => {
            let metrics = a_priori_metrics(&start, &plan.opts, &a)?;
            let (g, r) = adapt_to_metrics(start, &metrics, &a)?;
            (g, Some(r))
        }
    };
    if use_structured {
        let (labels, eval, warnings) = (g.labels.clone(), g.eval.clone(), g.warnings.clone());
        match structured_pml_native_with_curved(g, plan.opts.curved)? {
            pml_geom::Native::Unchanged(g_unch, rep) => {
                let m = g_unch.into_exchange_mesh()?;
                Ok((m, g_unch.labels, g_unch.eval, adapt, Some(rep), g_unch.warnings))
            }
            pml_geom::Native::Cells(cells, rep) => {
                Ok((Mesh::D2(cells), labels, eval, adapt, Some(rep), warnings))
            }
        }
    } else {
        let m = g.into_exchange_mesh()?;
        Ok((m, g.labels, g.eval, adapt, None, g.warnings))
    }
}

/// Meshes a model of either dimension for the neutral `Sizing` a physics reports (ADR 50): the sizes the model asks for, and in
/// addition the sizes the waves need. Region names of the sizing are resolved with the model's labels; inside structured PMLs
/// `pml_relax` coarsening is disabled (`1.0`).
pub fn mesh_for_sizing(list: &FeatureList, plan: &MeshPlan, sizing: &Sizing, sz: &PhysicsSizing) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    let (g, adapt, _) = mesh_for_sizing_with_pml_report(list, plan, sizing, sz)?;
    Ok((g, adapt))
}

pub(crate) fn mesh_for_sizing_with_pml_report(
    list: &FeatureList,
    plan: &MeshPlan,
    sizing: &Sizing,
    sz: &PhysicsSizing,
) -> Result<(MeshedGeom, Option<AdaptReport>, Option<PmlReport>), String> {
    let ev = evaluate(list)?;
    check_outward_inline_guard(&ev, plan.pml)?;
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    let mut eff_sz = *sz;
    if use_structured {
        eff_sz.pml_relax = 1.0;
    }
    let expanded = {
        let labels = ev.label_map().map_err(|e| e.to_string())?;
        let names: Vec<&str> = labels.domains.iter().map(|d| d.name.as_str()).collect();
        sizing.expand_wildcards(&names)
    };
    let sizing = &expanded;
    if ev.list.dim == 3 {
        return mesh_3d_plan(ev, plan, Some((sizing, &eff_sz)));
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let resolve = |n: &str| labels.domain(n);
    let pml: std::collections::HashSet<i32> = ev.pml_regions().iter().map(|p| labels.domains[p.region as usize].id).collect();
    let is_pml = |d: i32| pml.contains(&d);
    mesh_waves_2d(
        list,
        plan,
        &|m: &SimplexMesh| physics_size::vertex_sizes_for(sizing, &resolve, m, &is_pml, &eff_sz),
        &|m: &SimplexMesh, ceiling: f64| physics_size::metrics_for::<2>(sizing, &resolve, m, &is_pml, &eff_sz, ceiling),
    )
}

/// Meshes a model for study `Sizing` directly into an uncut [`Mesh`] (`Quad`/`Prism`/`Hex` preserved
/// with zero eager simplex splitting), returning `(Mesh, LabelMap, Arc<Evaluated>, Option<AdaptReport>, MeshReport, Vec<String>)`.
pub fn mesh_for_sizing_native(
    list: &FeatureList,
    plan: &MeshPlan,
    sizing: &Sizing,
    sz: &PhysicsSizing,
) -> Result<NativeMeshTuple, String> {
    let (mesh, labels, eval, adapt, pml_rep, warnings, eff_sz) =
        mesh_for_sizing_native_inner(list, plan, sizing, sz)?;
    let wave_sizes: Vec<Option<f64>> = wave_sizes_by_region(&eval, &labels, sizing, &eff_sz)?
        .into_iter()
        .map(|h| h.is_finite().then_some(h))
        .collect();
    let report = crate::meshing::request::summarise_native_mesh(
        &mesh,
        &labels,
        &eval,
        plan,
        adapt.as_ref(),
        pml_rep.as_ref(),
        &wave_sizes,
        &warnings,
    );
    Ok((mesh, labels, eval, adapt, report, warnings))
}

pub(crate) fn mesh_for_sizing_native_inner(
    list: &FeatureList,
    plan: &MeshPlan,
    sizing: &Sizing,
    sz: &PhysicsSizing,
) -> Result<NativeSizingRawTuple, String> {
    let ev = evaluate(list)?;
    check_outward_inline_guard(&ev, plan.pml)?;
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    let mut eff_sz = *sz;
    if use_structured {
        eff_sz.pml_relax = 1.0;
    }
    let expanded = {
        let labels = ev.label_map().map_err(|e| e.to_string())?;
        let names: Vec<&str> = labels.domains.iter().map(|d| d.name.as_str()).collect();
        sizing.expand_wildcards(&names)
    };
    if ev.list.dim == 3 {
        let (m, labels, eval, adapt, pml_rep, warnings) =
            mesh_3d_plan_native(ev, plan, Some((&expanded, &eff_sz)))?;
        return Ok((m, labels, eval, adapt, pml_rep, warnings, eff_sz));
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let resolve = |n: &str| labels.domain(n);
    let pml: std::collections::HashSet<i32> =
        ev.pml_regions().iter().map(|p| labels.domains[p.region as usize].id).collect();
    let is_pml = |d: i32| pml.contains(&d);
    let (m, labels_out, eval_out, adapt, pml_rep, warnings) = mesh_waves_2d_native(
        list,
        plan,
        &|m: &SimplexMesh| physics_size::vertex_sizes_for(&expanded, &resolve, m, &is_pml, &eff_sz),
        &|m: &SimplexMesh, ceiling: f64| {
            physics_size::metrics_for::<2>(&expanded, &resolve, m, &is_pml, &eff_sz, ceiling)
        },
    )?;
    Ok((m, labels_out, eval_out, adapt, pml_rep, warnings, eff_sz))
}

type VertexSizesFn<'a> = &'a dyn Fn(&SimplexMesh) -> Result<Vec<f64>, String>;
type MetricsFn<'a> = &'a dyn Fn(&SimplexMesh, f64) -> Result<Vec<Sym<2>>, String>;

fn mesh_waves_2d_phys(
    list: &FeatureList,
    plan: &MeshPlan,
    vertex_sizes: VertexSizesFn,
    wave_metrics: MetricsFn,
) -> Result<(MeshedGeom, Option<AdaptReport>, bool), String> {
    let mut opts = plan.opts.clone();
    let ev0 = evaluate(list)?;
    let use_structured = plan.pml.is_structured(!ev0.pml_regions().is_empty());
    let skip = pml_skip_mask(&ev0, use_structured);
    let (g, report) = match plan.mesher {
        Mesher::Simple => {
            let probe = mesh_geom_2d_with_skip(ev0, &opts, &skip)?;
            let hmin = vertex_sizes(&probe.mesh)?.into_iter().fold(f64::INFINITY, f64::min);
            if hmin < opts.h {
                opts.h = hmin;
                (mesh_geom_2d_with_skip(evaluate(list)?, &opts, &skip)?, None)
            } else {
                (probe, None)
            }
        }
        Mesher::Adaptive(a) => {
            let start = mesh_geom_2d_with_skip(ev0, &opts, &skip)?;
            let geo = a_priori_metrics(&start, &opts, &a)?;
            let phys = wave_metrics(&start.mesh, opts.h)?;
            let metrics: Vec<Sym<2>> = geo
                .iter()
                .zip(&phys)
                .map(|(g, p)| {
                    // the physics size is scaled like the geometry sizes, so `Adaptive::scale` stays one knob for the whole mesh
                    let p = [[p[0][0] / (a.scale * a.scale), p[0][1] / (a.scale * a.scale)], [p[1][0] / (a.scale * a.scale), p[1][1] / (a.scale * a.scale)]];
                    intersect(g, &p)
                })
                .collect();
            let (g, r) = adapt_to_metrics(start, &metrics, &a)?;
            (g, Some(r))
        }
    };
    Ok((g, report, use_structured))
}

fn mesh_waves_2d(
    list: &FeatureList,
    plan: &MeshPlan,
    vertex_sizes: VertexSizesFn,
    wave_metrics: MetricsFn,
) -> Result<(MeshedGeom, Option<AdaptReport>, Option<PmlReport>), String> {
    let (g, report, use_structured) = mesh_waves_2d_phys(list, plan, vertex_sizes, wave_metrics)?;
    if use_structured {
        let (g, pml_rep) = structured_pml_with_curved(g, plan.opts.curved)?;
        Ok((g, report, Some(pml_rep)))
    } else {
        Ok((g, report, None))
    }
}

fn mesh_waves_2d_native(
    list: &FeatureList,
    plan: &MeshPlan,
    vertex_sizes: VertexSizesFn,
    wave_metrics: MetricsFn,
) -> Result<NativeRawTuple, String> {
    let (g, report, use_structured) = mesh_waves_2d_phys(list, plan, vertex_sizes, wave_metrics)?;
    if use_structured {
        let (labels, eval, warnings) = (g.labels.clone(), g.eval.clone(), g.warnings.clone());
        match structured_pml_native_with_curved(g, plan.opts.curved)? {
            pml_geom::Native::Unchanged(g_unch, rep) => {
                let m = g_unch.into_exchange_mesh()?;
                Ok((m, g_unch.labels, g_unch.eval, report, Some(rep), g_unch.warnings))
            }
            pml_geom::Native::Cells(cells, rep) => {
                Ok((Mesh::D2(cells), labels, eval, report, Some(rep), warnings))
            }
        }
    } else {
        let m = g.into_exchange_mesh()?;
        Ok((m, g.labels, g.eval, report, None, g.warnings))
    }
}

/// The element size the waves of `sizing` need in every region of a 3D model, by region index; `f64::INFINITY` where a region carries no
/// wave. Absorbing pieces are relaxed by `sz.pml_relax`.
pub fn wave_sizes_by_region(ev: &Evaluated, labels: &LabelMap, sizing: &Sizing, sz: &PhysicsSizing) -> Result<Vec<f64>, String> {
    let pml: std::collections::HashSet<i32> = ev.pml_regions().iter().map(|p| labels.domains[p.region as usize].id).collect();
    let ids: Vec<i32> = labels.domains.iter().map(|d| d.id).collect();
    let names: Vec<&str> = labels.domains.iter().map(|d| d.name.as_str()).collect();
    physics_size::cell_sizes_for(&sizing.expand_wildcards(&names), &|n| labels.domain(n), &ids, &|d| pml.contains(&d), sz)
}

fn mesh_3d_physical_interior(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<(MeshedGeom, Option<AdaptReport>, bool), String> {
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let nr = ev.regions().len();
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    let mut eff_sz_opt: Option<PhysicsSizing> = None;
    let want = match waves {
        Some((sizing, sz)) => {
            let mut s = *sz;
            if use_structured {
                s.pml_relax = 1.0;
            }
            let v = wave_sizes_by_region(&ev, &labels, sizing, &s)?;
            eff_sz_opt = Some(s);
            v
        }
        None => vec![f64::INFINITY; nr],
    };
    let _ = eff_sz_opt;
    let skip = pml_skip_mask(&ev, use_structured);
    let start_sizes = match plan.mesher {
        Mesher::Simple => Sizes3 { region_h: want.clone(), skip: skip.clone() },
        Mesher::Adaptive(_) => Sizes3 { region_h: vec![], skip: skip.clone() },
    };
    // When adaptive meshing runs, curve the adapted physical interior BEFORE building structured PML shells
    // so physical tets are curved on the physical geometry and PML shells propagate `edge_nodes` in lockstep!
    let curve_after_adapt = plan.opts.curved && matches!(plan.mesher, Mesher::Adaptive(_));
    let mut start_opts = plan.opts.clone();
    if curve_after_adapt {
        start_opts.curved = false;
    }
    let (g, _, _) = geomdrive3::mesh_geom_3d_with(ev, &start_opts, &start_sizes)?;
    let (mut g, report) = match plan.mesher {
        Mesher::Simple => (g, None),
        Mesher::Adaptive(a) => {
            let geo = adapt3d::vertex_sizes3(&g, &plan.opts, &a, false)?;
            let mut size = geo;
            if waves.is_some() {
                let region_of: HashMap<i32, usize> = g.labels.domains.iter().map(|d| (d.id, d.region as usize)).collect();
                let mut vw = vec![f64::INFINITY; g.mesh.n_points()];
                for c in 0..g.mesh.n_cells() {
                    let h = region_of.get(&g.mesh.cell_labels[c]).map_or(f64::INFINITY, |&r| want[r]);
                    for &v in g.mesh.cell(c) {
                        vw[v as usize] = vw[v as usize].min(h);
                    }
                }
                for (s, w) in size.iter_mut().zip(&vw) {
                    *s = s.min(w * a.scale);
                }
            }
            let metrics: Vec<Sym<3>> = size.iter().map(|&h| isotropic::<3>(h.max(1e-300))).collect();
            let (g, r) = adapt3d::adapt_to_metrics3(g, &metrics, &a)?;
            (g, Some(r))
        }
    };
    if curve_after_adapt {
        let dom = geomdrive3::GeomDomain3::new_with(
            &g.eval,
            &g.labels,
            plan.opts.h,
            &Sizes3 { region_h: vec![], skip },
        )?;
        curve3::curve_solver_mesh(&mut g.mesh, &dom);
        if let Some(Mesh::D3(ref mut cm)) = g.cells {
            cm.set_edge_nodes(g.mesh.edge_node_list());
        }
    }
    Ok((g, report, use_structured))
}

fn uses_prismatic_stack_mesher(ev: &Evaluated, use_structured: bool) -> bool {
    let cx = ev.complex();
    if cx.has_sheets() {
        return true;
    }
    use_structured
        && ev
            .pml_regions()
            .iter()
            .any(|p| p.coord == geom::model::PmlCoordKind::Cylindrical && p.r_in > 0.0)
        && cx.stack.is_some()
        && cx.faces.iter().all(|f| {
            matches!(
                f.geom,
                geom::complex::FaceGeom::Cap { .. } | geom::complex::FaceGeom::Lateral { .. }
            )
        })
}

/// Meshes a 3D model (ADR 40, 50, 51, 52).
fn mesh_3d_plan(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<(MeshedGeom, Option<AdaptReport>, Option<PmlReport>), String> {
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    if uses_prismatic_stack_mesher(&ev, use_structured) {
        return sheet3d::mesh_sheet_stack_3d_with_pml(ev, plan, waves);
    }
    let (g, report, use_structured) = mesh_3d_physical_interior(ev, plan, waves)?;
    if use_structured {
        let (g, pml_rep) = pml_geom3::structured_pml_3d_with_curved(g, plan.opts.curved)?;
        Ok((g, report, Some(pml_rep)))
    } else {
        Ok((g, report, None))
    }
}

fn mesh_3d_plan_native(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<NativeRawTuple, String> {
    let use_structured = plan.pml.is_structured(!ev.pml_regions().is_empty());
    if uses_prismatic_stack_mesher(&ev, use_structured) {
        return sheet3d::mesh_sheet_stack_3d_native(ev, plan, waves);
    }
    let (g, report, use_structured) = mesh_3d_physical_interior(ev, plan, waves)?;
    if use_structured {
        let (labels, eval, warnings) = (g.labels.clone(), g.eval.clone(), g.warnings.clone());
        match pml_geom3::structured_pml_3d_native_with_curved(g, plan.opts.curved)? {
            pml_geom3::Native3::Unchanged(g_unch, rep) => {
                let m = g_unch.into_exchange_mesh()?;
                Ok((m, g_unch.labels, g_unch.eval, report, Some(rep), g_unch.warnings))
            }
            pml_geom3::Native3::Cells(cells, rep) => {
                Ok((Mesh::D3(cells), labels, eval, report, Some(rep), warnings))
            }
        }
    } else {
        let m = g.into_exchange_mesh()?;
        Ok((m, g.labels, g.eval, report, None, g.warnings))
    }
}

/// The isotropic metric at every vertex of a meshed model that the model itself asks for: the size of the region(s) the vertex belongs
/// to, the `Size` rules of the edges it lies on, and on a curved edge the size at which a chord stays within `max_dev` of the curve. All
/// times `a.scale`.
pub fn a_priori_metrics(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<2>>, String> {
    Ok(vertex_sizes(g, opts, a, false)?.iter().map(|&h| isotropic::<2>(h.max(1e-300))).collect())
}

/// The sizes the geometry insists on whatever the solution does: `Size` rules on regions and edges, the layers of an absorbing region, and
/// the curvature limit of the edges. Everywhere else the size is the global `opts.h` (a ceiling). Used as the floor of an estimator metric.
pub fn geometry_floor(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<2>>, String> {
    Ok(vertex_sizes(g, opts, a, true)?.iter().map(|&h| isotropic::<2>(h.max(1e-300))).collect())
}

fn vertex_sizes(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive, geometry_only: bool) -> Result<Vec<f64>, String> {
    if !(a.scale > 0.0 && a.scale.is_finite()) {
        return Err("the adaptation scale must be positive".into());
    }
    let ev = &g.eval;
    let m = &g.mesh;
    let pml_doms: std::collections::HashSet<i32> =
        ev.pml_regions().iter().map(|p| g.labels.domains[p.region as usize].id).collect();
    let has_pml_cells = !pml_doms.is_empty() && m.cell_labels.iter().any(|l| pml_doms.contains(l));
    let skip = pml_skip_mask(ev, !pml_doms.is_empty() && !has_pml_cells);
    let sizes = size_map_with_skip(ev, opts.h, &skip);
    let max_dev = opts.max_dev.unwrap_or(opts.h / 50.0);
    // region index of a domain label
    let region_of: HashMap<i32, usize> = g.labels.domains.iter().map(|d| (d.id, d.region as usize)).collect();
    let mut vh = vec![f64::INFINITY; m.n_points()];
    for c in 0..m.n_cells() {
        let mut h = region_of.get(&m.cell_labels[c]).map_or(opts.h, |&r| sizes.region_h[r]);
        if geometry_only && h >= opts.h * 0.999 {
            h = opts.h;
        }
        for &v in m.cell(c) {
            vh[v as usize] = vh[v as usize].min(h);
        }
    }
    // vertices on a boundary or an interface: the only ones an edge rule or a curvature limit can concern
    let constrained = constrained_vertices(m);
    let (lo, hi) = m.bbox();
    let tol = 1e-6 * (hi[0] - lo[0]).max(hi[1] - lo[1]).max(opts.h * 1e-3);
    for (i, ed) in ev.edges().iter().enumerate() {
        let curved = ev.min_radius_of_curvature(i as u32).is_finite();
        let mut h_edge = sizes.edge_rule[i].unwrap_or(f64::INFINITY);
        if curved {
            let r = ev.min_radius_of_curvature(i as u32);
            h_edge = h_edge.min(2.0 * (2.0 * r * max_dev).sqrt());
        }
        if !h_edge.is_finite() {
            continue;
        }
        let (blo, bhi) = (ed.bbox.lo, ed.bbox.hi);
        for &v in &constrained {
            let p = [m.points[2 * v as usize], m.points[2 * v as usize + 1]];
            if p[0] < blo[0] - tol || p[0] > bhi[0] + tol || p[1] < blo[2] - tol || p[1] > bhi[2] + tol {
                continue;
            }
            let (q, _) = ev.edge_closest(i as u32, [p[0], 0.0, p[1]]);
            if (q[0] - p[0]).hypot(q[2] - p[1]) <= tol {
                vh[v as usize] = vh[v as usize].min(h_edge);
            }
        }
    }
    Ok(vh.iter().map(|&h| h * a.scale).collect())
}

/// Vertices of the mesh that lie on a boundary facet or between cells of different domains.
fn constrained_vertices(m: &SimplexMesh) -> Vec<u32> {
    let mut seen: HashMap<(u32, u32), (u32, i32)> = HashMap::new();
    let mut out = vec![false; m.n_points()];
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for k in 0..3 {
            let (a, b) = (v[k].min(v[(k + 1) % 3]), v[k].max(v[(k + 1) % 3]));
            match seen.get_mut(&(a, b)) {
                None => {
                    seen.insert((a, b), (1, m.cell_labels[c]));
                }
                Some(e) => {
                    if e.1 != m.cell_labels[c] {
                        e.0 = 1;
                        e.1 = i32::MIN;
                    } else {
                        e.0 += 1;
                    }
                }
            }
        }
    }
    for ((a, b), (n, l)) in seen {
        if n == 1 || l == i32::MIN {
            out[a as usize] = true;
            out[b as usize] = true;
        }
    }
    (0..m.n_points() as u32).filter(|&v| out[v as usize]).collect()
}

/// Adapts a meshed model to a metric given at every vertex of its mesh. The metrics are graded and interpolated on the old mesh; the
/// vertices of the geometry stay, boundary and interface vertices are kept on the exact edges, and the labels of cells and facets are
/// carried through.
pub fn adapt_to_metrics(g: MeshedGeom, metrics: &[Sym<2>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String> {
    if g.mesh.dim != 2 {
        return Err("adaptation works on 2D meshes".into());
    }
    if metrics.len() != g.mesh.n_points() {
        return Err(format!("{} metrics for {} vertices", metrics.len(), g.mesh.n_points()));
    }
    let old = simplex_to_adapt::<2>(&g.mesh)?;
    let field = BackgroundMetric2::new_graded(&old, metrics, a.grade);
    adapt_to_field(g, &field, a)
}

/// Adapts a meshed model to a metric field defined anywhere in the plane, for example the final metric of an earlier, slightly different
/// model (`BackgroundMetric2` on that model's mesh). Same guarantees as `adapt_to_metrics`.
pub fn adapt_to_field(g: MeshedGeom, field: &dyn MetricField<2>, a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String> {
    if g.mesh.dim != 2 {
        return Err("adaptation works on 2D meshes".into());
    }
    let old = simplex_to_adapt::<2>(&g.mesh)?;
    let mut opts = AdaptOptions { max_cells: a.max_cells, ..AdaptOptions::default() };
    opts.topo.fixed_vertices = geometry_vertices(&g);
    opts.topo.frozen_labels = periodic_labels(&g);
    let result = {
        let ev = &g.eval;
        let labels = &g.labels;
        adapt2_with(&old, field, &opts, |w| GeomOracle2::from_work(ev, labels, w))?
    };
    let mesh = simplex_from_adapt(&result.mesh);
    let cells = Mesh::from_solver_mesh(&mesh)?;
    Ok((
        MeshedGeom {
            mesh,
            cells: Some(cells),
            labels: g.labels,
            eval: g.eval,
            warnings: g.warnings,
        },
        result.report,
    ))
}

/// Boundary labels of the edges (2D) or faces (3D) that take part in a periodic pair of the model. Their chains or patches stay
/// exactly as the start mesh has them, so the two sides of a pair keep matching vertex for vertex.
pub(crate) fn periodic_labels(g: &MeshedGeom) -> Vec<i32> {
    let mut v: Vec<i32> = vec![];
    for p in g.eval.periodic_pairs() {
        for set in [&p.a, &p.b] {
            if matches!(set.kind, EntityKind::Edges | EntityKind::Faces) {
                v.extend(set.ids.iter().filter_map(|&e| g.labels.label_of_edge(e)));
            }
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

/// Indices of the mesh vertices that sit on a vertex of the geometry.
fn geometry_vertices(g: &MeshedGeom) -> Vec<u32> {
    let m = &g.mesh;
    let (lo, hi) = m.bbox();
    let scale = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
    let tol = 1e-7 * scale;
    let mut out = vec![];
    for gv in g.eval.vertices() {
        let q = [gv.p[0], gv.p[2]];
        let best = (0..m.n_points())
            .map(|i| (i, (m.points[2 * i] - q[0]).hypot(m.points[2 * i + 1] - q[1])))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        if let Some((i, d)) = best {
            if d <= tol {
                out.push(i as u32);
            }
        }
    }
    out
}

/// The exact geometry as the oracle of the adaptation: where a vertex of a boundary or interface chain may go is the closest point on
/// the geom edge(s) of that chain, so refinement converges to the true circle or ellipse instead of the chord polygon it started from.
pub struct GeomOracle2<'a> {
    ev: &'a Evaluated,
    /// geom edges of every chain of the working mesh
    cands: Vec<Vec<u32>>,
}

impl<'a> GeomOracle2<'a> {
    pub fn from_work(ev: &'a Evaluated, labels: &LabelMap, w: &Work2) -> Self {
        let domain_of: Vec<i32> = {
            let mut v = vec![i32::MIN; ev.regions().len()];
            for d in &labels.domains {
                v[d.region as usize] = d.id;
            }
            v
        };
        let side = |r: Option<u32>| r.map_or(i32::MIN, |r| domain_of[r as usize]);
        let cands = w
            .chains
            .iter()
            .map(|c| match c.sig {
                Sig::Bnd(l) | Sig::Sheet(l) => labels.edges_of_label(l),
                Sig::Iface(a, b) => ev
                    .edges()
                    .iter()
                    .enumerate()
                    .filter(|(_, e)| {
                        let (l, r) = (side(e.left), side(e.right));
                        (l == a && r == b) || (l == b && r == a)
                    })
                    .map(|(i, _)| i as u32)
                    .collect(),
            })
            .collect();
        GeomOracle2 { ev, cands }
    }
}

impl Geometry2 for GeomOracle2<'_> {
    fn project(&self, p: [f64; 2], chain: u32) -> [f64; 2] {
        let mut best = (f64::INFINITY, p);
        for &e in &self.cands[chain as usize] {
            let bb = &self.ev.edges()[e as usize].bbox;
            // distance from the point to the box of the edge is a lower bound of the distance to the edge
            let dx = (bb.lo[0] - p[0]).max(p[0] - bb.hi[0]).max(0.0);
            let dz = (bb.lo[2] - p[1]).max(p[1] - bb.hi[2]).max(0.0);
            if dx.hypot(dz) >= best.0 {
                continue;
            }
            let (q, _) = self.ev.edge_closest(e, [p[0], 0.0, p[1]]);
            let d = (q[0] - p[0]).hypot(q[2] - p[1]);
            if d < best.0 {
                best = (d, [q[0], q[2]]);
            }
        }
        best.1
    }
}

// ------------------------------------------------------------------ caching

/// How a cached request was answered.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Reuse {
    /// the same model and plan as an earlier request: nothing was recomputed
    Hit,
    /// meshed afresh; `change` says how the model differs from the one meshed just before (`None` for the first request), which is what
    /// a mesh morph could use to skip the remeshing
    Meshed { change: Option<ChangeKind> },
}

/// Counters of a `MeshCache`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CacheStats {
    pub requests: usize,
    pub hits: usize,
    pub eval_seconds: f64,
    pub mesh_seconds: f64,
}

/// Meshes a series of related models (a parameter sweep, a refinement study, the same model for several physics) without repeating work.
///
/// * One geom `Evaluator` lives as long as the cache, so the primitives and pair intersections of features that did not change are reused
///   between evaluations.
/// * The 2D meshed model of every recent request is kept under the model and plan that made it; asking again, or extruding the same section
///   through different layers (`mesh_extruded`), costs a copy of the arrays.
///
/// The cache is keyed by content (the serialised feature list and the plan), never by identity, so equal models made separately hit.
pub struct MeshCache {
    evaluator: Evaluator,
    entries: Vec<(String, MeshedGeom)>,
    capacity: usize,
    last_eval: Option<std::sync::Arc<Evaluated>>,
    pub stats: CacheStats,
}

impl Default for MeshCache {
    fn default() -> Self {
        MeshCache::new(8)
    }
}

impl MeshCache {
    /// A cache that keeps the `capacity` most recent distinct requests.
    pub fn new(capacity: usize) -> Self {
        MeshCache { evaluator: Evaluator::new(), entries: vec![], capacity: capacity.max(1), last_eval: None, stats: CacheStats::default() }
    }

    fn key(list: &FeatureList, plan: &MeshPlan) -> Result<String, String> {
        let l = serde_json::to_string(list).map_err(|e| e.to_string())?;
        Ok(format!("{l}\u{1}{plan:?}"))
    }

    /// The meshed model for `list` and `plan`, from the cache when it is there.
    pub fn mesh(&mut self, list: &FeatureList, plan: &MeshPlan) -> Result<(MeshedGeom, Reuse), String> {
        self.stats.requests += 1;
        let key = Self::key(list, plan)?;
        if let Some(i) = self.entries.iter().position(|e| e.0 == key) {
            // most recently used last
            let e = self.entries.remove(i);
            let g = e.1.clone();
            self.entries.push(e);
            self.stats.hits += 1;
            return Ok((g, Reuse::Hit));
        }
        let t = std::time::Instant::now();
        let ev = self.evaluator.eval(list, Tier::Complex).map_err(|e| e.to_string())?;
        self.stats.eval_seconds += t.elapsed().as_secs_f64();
        let change = self.last_eval.as_ref().map(|prev| diff(prev, &ev).kind);
        let t = std::time::Instant::now();
        let g = mesh_evaluated(ev, plan)?;
        self.stats.mesh_seconds += t.elapsed().as_secs_f64();
        self.last_eval = Some(g.eval.clone());
        self.entries.push((key, g.clone()));
        if self.entries.len() > self.capacity {
            self.entries.remove(0);
        }
        Ok((g, Reuse::Meshed { change }))
    }

    /// The section of `list` swept along y through `layers`. The section is the cached part: changing only the layers re-extrudes it.
    pub fn mesh_extruded(&mut self, list: &FeatureList, plan: &MeshPlan, layers: &[crate::meshing::mesher::Layer]) -> Result<(MeshedGeom, Reuse), String> {
        let (mut g, reuse) = self.mesh(list, plan)?;
        g.mesh = crate::meshing::mesher::extrude_y(&g.mesh, layers, 0.0)?;
        g.cells = Some(Mesh::from_solver_mesh(&g.mesh)?);
        Ok((g, reuse))
    }
}
