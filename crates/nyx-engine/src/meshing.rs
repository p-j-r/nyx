// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Meshing seam. `MeshedModel` is a mesh together with the names of the model (regions, boundaries) and the label resolver the physics
//! use; `Mesher` builds one from a feature list. `FemcoreMesher` is the **only** code of the engine that calls femcore's mesher
//! (`meshing::mesh_model`, `mesher::extrude_y`): restructure phase 2 replaces the mesh types and moves those modules, and rewrites this
//! one function. Names and labels come from `geom` (`LabelMap`), which numbers domains from 1 in overlay order.

use crate::project::{MeshSettings, PmlModeRequest};
use geom::eval::Evaluated;
use geom::labels::LabelMap;
use geom::model::{FeatureList, PmlGrowMode};
use geom::select::{EntityKind, Selector};
use nyx_physics::{Error, LabelResolver, Mesh, Sizing};
use serde_json::{json, Value};
use std::sync::Arc;

pub struct MeshedModel {
    mesh: Arc<Mesh>,
    eval: Arc<Evaluated>,
    labels: LabelMap,
    /// the boundary labels of the two end faces of a 2D model swept along y (`ymin`, `ymax`)
    ends: Option<(i32, i32)>,
    /// the boundary label sets of the two end planes (`ymin`, `ymax`) of a 3D model swept along y
    ends_3d: Option<(Vec<i32>, Vec<i32>)>,
    /// what the mesher reports (domains, cells per domain, the element size the waves asked for, the adaptation); absent for a swept mesh
    report: Option<femcore::meshing::request::MeshReport>,
    /// true when produced by multi-pass a posteriori AMR (`AF-1`) rather than the initial `Mesher::mesh` call
    amr_adapted: bool,
    /// the solver view of the mesh (edges, P2 numbering), made on first use: what the plots read
    solver: std::sync::OnceLock<Arc<femcore::mesh::SimplexMesh>>,
}

impl MeshedModel {
    /// From an uncut mixed-cell exchange [`Mesh`] (`Quad`, `Prism`, `Hex` preserved) without eager simplex splitting.
    /// The simplex view [`MeshedModel::solver_mesh`] is built lazily on first display/export access.
    pub fn from_exchange(
        mesh: Mesh,
        eval: Arc<Evaluated>,
        labels: LabelMap,
        swept: bool,
    ) -> Result<MeshedModel, Error> {
        let top = match &mesh {
            Mesh::D2(cm) => cm.facet_label.iter().copied().max().unwrap_or(0),
            Mesh::D3(cm) => cm.facet_label.iter().copied().max().unwrap_or(0),
        };
        let (ends, ends_3d) = if swept && eval.list.dim == 2 {
            (Some((top - 1, top)), None)
        } else if swept && eval.list.dim == 3 {
            if let Mesh::D3(cm) = &mesh {
                if !cm.points.is_empty() {
                    let mut lo_y = f64::INFINITY;
                    let mut hi_y = f64::NEG_INFINITY;
                    for p in &cm.points {
                        lo_y = lo_y.min(p[1]);
                        hi_y = hi_y.max(p[1]);
                    }
                    let tol = 1e-8 * (hi_y - lo_y).abs().max(1e-30);
                    let mut ymin_l = Vec::new();
                    let mut ymax_l = Vec::new();
                    for f in 0..cm.n_facets() {
                        let fv = cm.facet(f);
                        if fv.iter().all(|&v| (cm.points[v as usize][1] - lo_y).abs() <= tol) {
                            ymin_l.push(cm.facet_label[f]);
                        } else if fv.iter().all(|&v| (cm.points[v as usize][1] - hi_y).abs() <= tol) {
                            ymax_l.push(cm.facet_label[f]);
                        }
                    }
                    ymin_l.sort_unstable();
                    ymin_l.dedup();
                    ymax_l.sort_unstable();
                    ymax_l.dedup();
                    (None, Some((ymin_l, ymax_l)))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            }
        } else {
            (None, None)
        };
        Ok(MeshedModel {
            mesh: Arc::new(mesh),
            eval,
            labels,
            ends,
            ends_3d,
            report: None,
            amr_adapted: false,
            solver: std::sync::OnceLock::new(),
        })
    }

    /// From the mesh the mesher made (the solver mesh of femcore) to the exchange mesh the physics are given.
    pub fn new(
        mesh: femcore::mesh::SimplexMesh,
        eval: Arc<Evaluated>,
        labels: LabelMap,
        swept: bool,
    ) -> Result<MeshedModel, Error> {
        let exchange = Mesh::from_solver_mesh(&mesh).map_err(Error::Invalid)?;
        let model = Self::from_exchange(exchange, eval, labels, swept)?;
        let _ = model.solver.set(Arc::new(mesh));
        Ok(model)
    }

    /// The same, with the report the mesher made of it.
    pub(crate) fn with_report(mut self, report: femcore::meshing::request::MeshReport) -> MeshedModel {
        self.report = Some(report);
        self
    }

    /// Marks this model as produced by multi-pass AMR (`AF-1`).
    pub(crate) fn with_amr_adapted(mut self) -> MeshedModel {
        self.amr_adapted = true;
        self
    }

    /// True when this model was adapted during multi-pass AMR rather than built solely from `MeshSettings`.
    pub fn is_amr_adapted(&self) -> bool {
        self.amr_adapted
    }

    pub fn mesh(&self) -> &Arc<Mesh> {
        &self.mesh
    }

    /// The mesh in the numbering the plots read (vertices, edges, P2 dofs): the exchange mesh cut into simplices, as every solve sees it.
    /// What a front-end draws results on.
    pub fn solver_mesh(&self) -> Arc<femcore::mesh::SimplexMesh> {
        self.solver
            .get_or_init(|| Arc::new(self.mesh.discretise().expect("a mesh the mesher made can be discretised")))
            .clone()
    }

    /// The geometry the mesh was made from, evaluated: what a front-end needs to name an entity it points at (selector suggestions).
    pub fn evaluated(&self) -> &Arc<Evaluated> {
        &self.eval
    }

    /// The names of the model: region and boundary numbers.
    pub fn labels(&self) -> &LabelMap {
        &self.labels
    }

    /// The boundary labels of the two end faces of a model swept along y.
    pub fn swept_ends(&self) -> Option<(i32, i32)> {
        self.ends
    }

    pub fn report(&self) -> Option<&femcore::meshing::request::MeshReport> {
        self.report.as_ref()
    }

    /// Region labels with their names.
    pub fn regions(&self) -> Vec<(i32, String)> {
        self.labels.domains.iter().map(|d| (d.id, d.name.clone())).collect()
    }

    /// What a person or a host needs to choose names: the regions, the boundary names and the size of the mesh.
    pub fn info(&self) -> Value {
        let mut boundaries: Vec<&str> = self.labels.boundaries.iter().map(|b| b.name.as_str()).collect();
        boundaries.sort_unstable();
        boundaries.dedup();
        json!({
            "cells": self.mesh.n_cells(),
            "dim": self.mesh.dim(),
            "regions": self.labels.domains.iter().map(|d| json!({ "label": d.id, "name": d.name })).collect::<Vec<_>>(),
            "boundaries": boundaries,
            "groups": self.labels.groups.iter().map(|g| json!({ "name": g.name, "rule": g.rule })).collect::<Vec<_>>(),
            "swept_ends": (self.ends.is_some() || self.ends_3d.is_some()).then(|| json!(["ymin", "ymax"])),
            "report": self.report,
        })
    }
}

fn resolve_domains(labels: &LabelMap, pattern: &str) -> Result<Vec<i32>, Error> {
    let mut v: Vec<i32> = vec![];
    for item in pattern.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let before = v.len();
        match item.strip_suffix('*') {
            // `sapphire*` is `sapphire` and its absorbing pieces `sapphire.pml_*`
            Some(base) => v.extend(
                labels
                    .domains
                    .iter()
                    .filter(|d| d.name == base || d.name.strip_prefix(base).is_some_and(|r| r.starts_with('.')))
                    .map(|d| d.id),
            ),
            None => v.extend(labels.domain(item)),
        }
        if v.len() == before {
            return Err(Error::invalid(format!(
                "no region '{item}' (the regions are: {})",
                labels.domains.iter().map(|d| d.name.as_str()).collect::<Vec<_>>().join(", ")
            )));
        }
    }
    v.sort_unstable();
    v.dedup();
    Ok(v)
}

impl LabelResolver for MeshedModel {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        resolve_domains(&self.labels, region)
    }

    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        let s = selector.trim();
        if s == "ymin" {
            if let Some((a, _)) = self.ends {
                return Ok(vec![a]);
            }
            if let Some((a, _)) = &self.ends_3d {
                return Ok(a.clone());
            }
            return Err(Error::invalid("ymin and ymax are the ends of a model swept along y; this model is not swept"));
        }
        if s == "ymax" {
            if let Some((_, b)) = self.ends {
                return Ok(vec![b]);
            }
            if let Some((_, b)) = &self.ends_3d {
                return Ok(b.clone());
            }
            return Err(Error::invalid("ymin and ymax are the ends of a model swept along y; this model is not swept"));
        }
        // a boundary name (`rib.top`, `air|rib`) or a tag first; these are not always valid selector text
        let named = self.labels.boundary(s);
        if !named.is_empty() {
            return Ok(named);
        }
        if let Some((lhs, rhs)) = s.split_once('|') {
            if !lhs.contains('(') && !rhs.contains('(') && !lhs.contains('|') && !rhs.contains('|') {
                let swapped = format!("{}|{}", rhs.trim(), lhs.trim());
                let named_swapped = self.labels.boundary(&swapped);
                if !named_swapped.is_empty() {
                    return Ok(named_swapped);
                }
            }
        }
        let sel = Selector::parse(s).map_err(|e| Error::invalid(format!("selector '{s}': {e}")))?;
        let r = self.eval.resolve(&sel).map_err(|e| Error::invalid(format!("selector '{s}': {e}")))?;
        match r.set.kind {
            EntityKind::Regions => Ok(r.set.ids.iter().map(|&i| self.labels.domains[i as usize].id).collect()),
            EntityKind::Vertices => Err(Error::invalid(format!("'{s}' selects vertices, which carry no condition"))),
            _ => {
                let mut v: Vec<i32> = r.set.ids.iter().filter_map(|&e| self.labels.label_of_edge(e)).collect();
                v.sort_unstable();
                v.dedup();
                if v.is_empty() {
                    return Err(Error::invalid(format!("'{s}' selects no boundary of this model")));
                }
                Ok(v)
            }
        }
    }

    fn sheet_thickness_m(&self) -> std::collections::BTreeMap<i32, f64> {
        self.labels.sheet_thickness_m.clone()
    }

    fn sheet_etch_depth_m(&self) -> std::collections::BTreeMap<i32, f64> {
        self.labels.sheet_etch_depth_m.clone()
    }
}

/// Builds the mesh of a model. `sizing` is what the physics of the project want resolved (`MeshSettings::physics_sizing`); the engine
/// passes it only when that flag is on.
pub trait Mesher: Send + Sync {
    fn mesh(
        &self,
        geometry: &FeatureList,
        settings: &MeshSettings,
        sizing: Option<&Sizing>,
    ) -> Result<MeshedModel, Error>;
}

/// A message of the mesher as an error: the text of a cancelled job is `Cancelled`, anything else an invalid model.
fn mesh_error(message: String) -> Error {
    if message == femcore::progress::CANCELLED {
        Error::Cancelled
    } else {
        Error::Invalid(message)
    }
}

/// The mesher of `femcore`, driven by a `MeshRequest` (ADR 54): 2D and 3D models, simple or adaptive, inline or structured absorbing
/// layers, element sizes from the waves of the study. A 2D model with `sweep` layers is extruded along y into a 3D mesh afterwards.
pub struct FemcoreMesher;

fn resolve_sweep_layers(
    layers: &[crate::project::SweepLayer],
    labels: &LabelMap,
) -> Result<Vec<femcore::meshing::mesher::RemappedLayer>, Error> {
    let all_names = || labels.domains.iter().map(|d| d.name.as_str()).collect::<Vec<_>>().join(", ");
    let mut out = Vec::with_capacity(layers.len());
    for l in layers {
        let fill = match &l.fill {
            Some(name) => Some(labels.domain(name).ok_or_else(|| {
                Error::invalid(format!("sweep layer fill: no region '{name}' (the regions are: {})", all_names()))
            })?),
            None => None,
        };
        let mut remap = std::collections::BTreeMap::new();
        if let Some(map) = &l.regions {
            for (src, dst) in map {
                let sid = labels.domain(src).ok_or_else(|| {
                    Error::invalid(format!("sweep layer regions: no region '{src}' (the regions are: {})", all_names()))
                })?;
                let did = labels.domain(dst).ok_or_else(|| {
                    Error::invalid(format!("sweep layer regions: no region '{dst}' (the regions are: {})", all_names()))
                })?;
                remap.insert(sid, did);
            }
        }
        out.push((
            femcore::meshing::mesher::Layer { thickness: l.thickness, n: l.n, domain_offset: l.domain_offset, fill },
            remap,
        ));
    }
    Ok(out)
}

/// Validates `MeshSettings` and builds the corresponding [`femcore::meshing::request::MeshRequest`].
pub(crate) fn build_mesh_request(
    s: &MeshSettings,
    sizing: Option<&Sizing>,
) -> Result<femcore::meshing::request::MeshRequest, Error> {
    use femcore::meshing::request::{MeshRequest, MesherRequest, WavesRequest};
    if !(s.h > 0.0 && s.h.is_finite()) {
        return Err(Error::invalid("the element size must be a positive finite number"));
    }
    if let Some(gc) = s.gap_cells {
        if !(gc > 0.0 && gc.is_finite()) {
            return Err(Error::invalid("`gap_cells` must be a positive finite number"));
        }
    }
    if let Some(mc) = s.max_cells {
        if mc < 1 {
            return Err(Error::invalid("`max_cells` must be an integer >= 1"));
        }
    }
    if let Some(dz) = s.dz_min {
        if !(dz > 0.0 && dz.is_finite()) {
            return Err(Error::invalid("`dz_min` must be a positive finite number"));
        }
    }
    if let Some(ms) = s.max_z_slabs {
        if ms < 1 {
            return Err(Error::invalid("`max_z_slabs` must be an integer >= 1"));
        }
    }
    if let Some(zr) = s.z_growth_ratio {
        if !(zr > 1.0 && zr.is_finite()) {
            return Err(Error::invalid("`z_growth_ratio` must be a finite number strictly greater than 1.0"));
        }
    }
    if let Some(gr) = s.adapt_grade {
        if !(gr > 0.0 && gr.is_finite()) {
            return Err(Error::invalid("`adapt_grade` must be a positive finite number"));
        }
    }
    if let Some(mc) = s.adapt_max_cells {
        if mc == 0 {
            return Err(Error::invalid("`adapt_max_cells` must be at least 1"));
        }
    }
    if let Some(se) = &s.singular_edges {
        if !matches!(se.as_str(), "none" | "quarter_point") {
            return Err(Error::invalid(format!(
                "`singular_edges` must be \"none\" or \"quarter_point\", got \"{se}\""
            )));
        }
    }
    if let Some(xi) = s.quarter_point_xi {
        if !(xi.is_finite() && (0.25..=0.50).contains(&xi)) {
            return Err(Error::invalid("`quarter_point_xi` must be in [0.25, 0.50]"));
        }
    }
    let pml_mode = s.resolved_pml_mode();
    let pml_relax = if pml_mode == PmlModeRequest::Inline { 1.5 } else { 1.0 };
    let waves = match (s.physics_sizing, sizing) {
        (true, Some(sz)) => Some(WavesRequest {
            sizing: sz.clone(),
            points_per_wavelength: s.points_per_wavelength,
            pml_relax,
            eigen_headroom: sz.headroom.max(1.0),
        }),
        _ => None,
    };
    let grade = s.adapt_grade.unwrap_or(0.4);
    let max_cells = match (s.max_cells, s.adapt_max_cells) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => Some(2_000_000),
    };
    // wave sizing adapts the mesh to the waves; an explicit `adapt_scale` sets the factor on every size
    let mesher = match (s.adapt_scale, waves.is_some()) {
        (Some(scale), _) => MesherRequest::Adaptive { scale, grade, max_cells },
        (None, true) => MesherRequest::Adaptive { scale: 1.0, grade, max_cells },
        (None, false) => MesherRequest::Simple,
    };
    Ok(MeshRequest {
        h: s.h,
        max_dev: s.max_dev,
        mesher,
        pml_mode,
        structured_pml: s.structured_pml,
        waves,
        curved: s.curved,
        gap_cells: s.gap_cells,
        z_coarsen: s.z_coarsen,
        max_cells: s.max_cells,
        dz_min: s.dz_min,
        max_z_slabs: s.max_z_slabs,
        z_growth_ratio: s.z_growth_ratio,
        max_steiner_multiplier: s.max_steiner_multiplier,
        strict_mesh: s.strict_mesh,
        singular_edges: s.singular_edges.clone(),
        quarter_point_xi: s.quarter_point_xi,
        corner_grading: s.corner_grading,
    })
}

impl Mesher for FemcoreMesher {
    fn mesh(&self, geometry: &FeatureList, s: &MeshSettings, sizing: Option<&Sizing>) -> Result<MeshedModel, Error> {
        let model = self.mesh_all(geometry, s, sizing)?;
        match &s.active_domains {
            Some(active) => model.restrict_to_active(active),
            None => Ok(model),
        }
    }
}

impl FemcoreMesher {
    /// Every region of the model meshed (`active_domains` is applied by `mesh`).
    fn mesh_all(
        &self,
        geometry: &FeatureList,
        s: &MeshSettings,
        sizing: Option<&Sizing>,
    ) -> Result<MeshedModel, Error> {
        use femcore::meshing::request::mesh_request_native;
        let req = build_mesh_request(s, sizing)?;
        if let Some(layers) = &s.sweep {
            if geometry.dim == 3 {
                let _stage = femcore::progress::stage(format!("sweeping 3D geometry along y (element size {} m)", s.h));
                let scale_factor = s.adapt_scale.unwrap_or(1.0);
                if !(scale_factor > 0.0 && scale_factor.is_finite()) {
                    return Err(Error::invalid("the adaptation scale must be positive"));
                }
                let ev = femcore::meshing::geomdrive::evaluate(geometry).map_err(mesh_error)?;
                let pml_mode = s.resolved_pml_mode();
                if pml_mode == PmlModeRequest::Inline && ev.pml_regions().iter().any(|p| p.mode == PmlGrowMode::Outward)
                {
                    return Err(Error::invalid(
                        "outward absorbing layers have no inline CAD cells; use PmlMode::Auto or PmlMode::Structured",
                    ));
                }
                let use_structured_pml = pml_mode != PmlModeRequest::Inline && !ev.pml_regions().is_empty();
                let pml_relax = if pml_mode == PmlModeRequest::Inline { 1.5 } else { 1.0 };
                let labels = ev.label_map().map_err(|e| Error::invalid(e.to_string()))?;
                let ls = resolve_sweep_layers(layers, &labels)?;
                let nr = ev.regions().len();
                let region_h = match (s.physics_sizing, sizing) {
                    (true, Some(sz)) => {
                        let ps = femcore::meshing::physics_size::PhysicsSizing {
                            points_per_wavelength: s.points_per_wavelength,
                            pml_relax,
                            eigen_headroom: sz.headroom.max(1.0),
                        };
                        let names: Vec<&str> = labels.domains.iter().map(|d| d.name.as_str()).collect();
                        femcore::meshing::wave_sizes_by_region(&ev, &labels, &sz.expand_wildcards(&names), &ps)
                            .map_err(mesh_error)?
                            .into_iter()
                            .map(|w| w * scale_factor)
                            .collect()
                    }
                    _ => vec![f64::INFINITY; nr],
                };
                let mut skip = vec![false; nr];
                if use_structured_pml {
                    for p in ev.pml_regions() {
                        skip[p.region as usize] = true;
                    }
                }
                let sizes = femcore::meshing::geomdrive3::Sizes3 { region_h, skip };
                let mesh = femcore::meshing::mesher::sweep_model_3d(
                    &ev,
                    &labels,
                    s.h * scale_factor,
                    s.max_dev,
                    25.0,
                    &sizes,
                    &ls,
                )
                .map_err(mesh_error)?;
                let ev_arc = Arc::new(ev);
                if use_structured_pml {
                    let g = femcore::meshing::geomdrive::MeshedGeom {
                        mesh,
                        cells: None,
                        labels: labels.clone(),
                        eval: ev_arc.clone(),
                        warnings: vec![],
                    };
                    return match femcore::meshing::pml_geom3::structured_pml_3d_native_with_curved(g, s.curved)
                        .map_err(mesh_error)?
                    {
                        femcore::meshing::pml_geom3::Native3::Unchanged(g, _) => {
                            MeshedModel::new(g.mesh, g.eval, g.labels, true)
                        }
                        femcore::meshing::pml_geom3::Native3::Cells(cells, _) => {
                            MeshedModel::from_exchange(Mesh::D3(cells), ev_arc, labels, true)
                        }
                    };
                }
                return MeshedModel::new(mesh, ev_arc, labels, true);
            }
        }
        // the stage a job reports while it meshes (its log gets the mesher's laps below it)
        let _stage = femcore::progress::stage(format!("meshing the geometry (element size {} m)", s.h));
        let (mesh, labels, eval, _adapt, report, _warnings) =
            mesh_request_native(geometry, &req).map_err(mesh_error)?;
        if let Some(layers) = &s.sweep {
            let ls = resolve_sweep_layers(layers, &labels)?;
            let mut sm = mesh.discretise().map_err(Error::Invalid)?;
            // numbered in the order of their coordinates, so that opposite side walls of a periodic cell are cut the same way
            // and a Floquet condition can pair them (`sort_nodes`)
            sm = femcore::meshing::geomdrive::sort_nodes(&sm);
            sm = femcore::meshing::mesher::extrude_y_with_remap(&sm, &ls, 0.0).map_err(mesh_error)?;
            return MeshedModel::new(sm, eval, labels, true);
        }
        Ok(MeshedModel::from_exchange(mesh, eval, labels, false)?.with_report(report))
    }
}

impl MeshedModel {
    /// A second handle on the same mesh, geometry and names (the arrays are shared, not copied): what a chain given a borrowed model
    /// keeps in its table of the meshes each physics was solved on.
    pub(crate) fn shallow_clone(&self) -> MeshedModel {
        let solver = std::sync::OnceLock::new();
        if let Some(s) = self.solver.get() {
            let _ = solver.set(s.clone());
        }
        MeshedModel {
            mesh: self.mesh.clone(),
            eval: self.eval.clone(),
            labels: self.labels.clone(),
            ends: self.ends,
            ends_3d: self.ends_3d.clone(),
            report: self.report.clone(),
            amr_adapted: self.amr_adapted,
            solver,
        }
    }

    /// The model with only the cells of the regions `active` names (region names or `name*` patterns, ADR 135): the other cells are dropped,
    /// the labels and names of the model are kept (so a field transfer still pairs regions by name). The model itself when every cell is
    /// already active.
    pub fn restrict_to_active(self, active: &[String]) -> Result<MeshedModel, Error> {
        let mut keep: Vec<i32> = vec![];
        for pattern in active {
            keep.extend(resolve_domains(&self.labels, pattern)?);
        }
        keep.sort_unstable();
        keep.dedup();
        let all_active = match self.mesh.as_ref() {
            Mesh::D2(cm) => cm.label.iter().all(|l| keep.binary_search(l).is_ok()),
            Mesh::D3(cm) => cm.label.iter().all(|l| keep.binary_search(l).is_ok()),
        };
        if all_active {
            return Ok(self);
        }
        let mut mesh = femcore::transfer::restrict_to_domains(&self.mesh, &keep)
            .map_err(|e| Error::invalid(format!("active_domains {active:?}: {e}")))?;
        match &mut mesh {
            Mesh::D2(cm) => {
                let existing: std::collections::HashSet<[u32; 2]> = (0..cm.n_facets())
                    .map(|f| {
                        let fv = cm.facet(f);
                        [fv[0].min(fv[1]), fv[0].max(fv[1])]
                    })
                    .collect();
                let missing: Vec<[u32; 2]> = cm
                    .boundary_faces_ordered()
                    .into_iter()
                    .filter_map(|fv| {
                        if fv.len() == 2 {
                            let k = [fv[0].min(fv[1]), fv[0].max(fv[1])];
                            (!existing.contains(&k)).then_some([fv[0], fv[1]])
                        } else {
                            None
                        }
                    })
                    .collect();
                if !missing.is_empty() {
                    let b = self.eval.bounds();
                    let scale = (b.hi[0] - b.lo[0]).max(b.hi[2] - b.lo[2]).max(1e-30);
                    let max_dev = 1e-4 * scale;
                    let tol = 1e-6 * scale;
                    let mut chains: Vec<(usize, i32, Vec<[f64; 2]>)> = Vec::new();
                    for (i, ed) in self.eval.edges().iter().enumerate() {
                        let Some(label) = self.labels.label_of_edge(i as u32) else {
                            continue;
                        };
                        let mut pts: Vec<[f64; 2]> = self
                            .eval
                            .discretize_edge(i as u32, geom::query::DiscretizeOpts { max_dev, max_len: scale })
                            .into_iter()
                            .map(|(_, p)| [p[0], p[2]])
                            .collect();
                        if let Some(first) = pts.first_mut() {
                            let p0 = self.eval.vertices()[ed.v0 as usize].p;
                            *first = [p0[0], p0[2]];
                        }
                        if let Some(last) = pts.last_mut() {
                            let p1 = self.eval.vertices()[ed.v1 as usize].p;
                            *last = [p1[0], p1[2]];
                        }
                        chains.push((i, label, pts));
                    }
                    let dist_seg = |p: [f64; 2], a: [f64; 2], b: [f64; 2]| -> f64 {
                        let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
                        let l2 = ab[0] * ab[0] + ab[1] * ab[1];
                        let t = if l2 > 0.0 { ((ap[0] * ab[0] + ap[1] * ab[1]) / l2).clamp(0.0, 1.0) } else { 0.0 };
                        ((a[0] + t * ab[0] - p[0]).powi(2) + (a[1] + t * ab[1] - p[1]).powi(2)).sqrt()
                    };
                    for [a, b_v] in missing {
                        let pa = cm.points[a as usize];
                        let pb = cm.points[b_v as usize];
                        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
                        let mut best: Option<(f64, usize, i32)> = None;
                        for &(part, label, ref pts) in &chains {
                            for seg in pts.windows(2) {
                                let d = dist_seg(mid, seg[0], seg[1]);
                                if best.is_none_or(|(bd, bp, _)| d < bd - tol || (d <= bd + tol && part >= bp)) {
                                    best = Some((d, part, label));
                                }
                            }
                        }
                        if let Some((_, _, lab)) = best {
                            cm.push_facet(&[a, b_v], lab);
                        }
                    }
                }
            }
            Mesh::D3(cm) => {
                let fkey = |vs: &[u32]| -> [u32; 4] {
                    let mut k = [u32::MAX; 4];
                    k[..vs.len().min(4)].copy_from_slice(&vs[..vs.len().min(4)]);
                    k[..vs.len().min(4)].sort_unstable();
                    k
                };
                let existing: std::collections::HashSet<[u32; 4]> =
                    (0..cm.n_facets()).map(|f| fkey(cm.facet(f))).collect();
                let missing: Vec<Vec<u32>> =
                    cm.boundary_faces_ordered().into_iter().filter(|fv| !existing.contains(&fkey(fv))).collect();
                if !missing.is_empty() && self.eval.list.dim == 3 {
                    let cx = self.eval.complex();
                    let nf = cx.facet_count();
                    for fv in missing {
                        let mut mid = [0.0_f64; 3];
                        for &v in &fv {
                            let p = cm.points[v as usize];
                            mid[0] += p[0] / (fv.len() as f64);
                            mid[1] += p[1] / (fv.len() as f64);
                            mid[2] += p[2] / (fv.len() as f64);
                        }
                        let mut best: Option<(f64, i32)> = None;
                        for fi in 0..nf as u32 {
                            let Some(lab) = self.labels.label_of_edge(fi) else {
                                continue;
                            };
                            let q = cx.facet_project(fi, mid);
                            let d =
                                ((q[0] - mid[0]).powi(2) + (q[1] - mid[1]).powi(2) + (q[2] - mid[2]).powi(2)).sqrt();
                            if best.is_none_or(|(bd, _)| d < bd) {
                                best = Some((d, lab));
                            }
                        }
                        if let Some((_, lab)) = best {
                            cm.push_facet(&fv, lab);
                        }
                    }
                }
            }
        }
        Ok(MeshedModel {
            mesh: Arc::new(mesh),
            eval: self.eval,
            labels: self.labels,
            ends: self.ends,
            ends_3d: self.ends_3d,
            report: None,
            amr_adapted: self.amr_adapted,
            solver: std::sync::OnceLock::new(),
        })
    }

    /// The canonical names of the regions by label (what a cross-mesh transfer pairs cells by).
    pub fn domain_name_map(&self) -> std::collections::BTreeMap<i32, String> {
        self.regions().into_iter().collect()
    }

    /// The canonical names of the boundaries by label.
    pub fn boundary_name_map(&self) -> std::collections::BTreeMap<i32, String> {
        self.boundary_names().into_iter().collect()
    }

    /// Selectors, simplest first, that pick exactly the boundaries with these numbers; empty when no selector does. What the window offers
    /// after a click on the mesh.
    pub fn suggest(&self, labels: &[i32]) -> Vec<String> {
        if let Some((a, b)) = self.ends {
            if labels == [a] {
                return vec!["ymin".into()];
            }
            if labels == [b] {
                return vec!["ymax".into()];
            }
        }
        let ids: Vec<u32> = self
            .labels
            .boundaries
            .iter()
            .filter(|b| labels.contains(&b.id))
            .flat_map(|b| b.entities.iter().copied())
            .collect();
        if ids.is_empty() {
            return vec![];
        }
        let set = geom::select::EntitySet::new(EntityKind::Edges, ids);
        let mut out: Vec<String> = self.eval.suggest(&set).iter().map(|s| s.to_string()).collect();
        out.dedup();
        out.truncate(6);
        out
    }

    /// The names of the boundaries with their numbers (a boundary name may stand for several numbers).
    pub fn boundary_names(&self) -> Vec<(i32, String)> {
        self.labels.boundaries.iter().map(|b| (b.id, b.name.clone())).collect()
    }
}
