// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Meshing seam. `MeshedModel` is a mesh together with the names of the model (regions, boundaries) and the label resolver the physics
//! use; `Mesher` builds one from a feature list. `FemcoreMesher` is the **only** code of the engine that calls femcore's mesher
//! (`meshing::mesh_model`, `mesher::extrude_y`): restructure phase 2 replaces the mesh types and moves those modules, and rewrites this
//! one function. Names and labels come from `geom` (`LabelMap`), which numbers domains from 1 in overlay order.

use crate::project::MeshSettings;
use geom::eval::Evaluated;
use geom::labels::LabelMap;
use geom::model::FeatureList;
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
    /// From the mesh the mesher made (the solver mesh of femcore) to the exchange mesh the physics are given.
    pub fn new(
        mesh: femcore::mesh::SimplexMesh,
        eval: Arc<Evaluated>,
        labels: LabelMap,
        swept: bool,
    ) -> Result<MeshedModel, Error> {
        let top = mesh.facet_labels.iter().copied().max().unwrap_or(0);
        let (ends, ends_3d) = if swept && eval.list.dim == 2 {
            (Some((top - 1, top)), None)
        } else if swept && eval.list.dim == 3 && mesh.dim == 3 && mesh.n_points() > 0 {
            let (lo, hi) = mesh.bbox();
            let tol = 1e-8 * (hi[1] - lo[1]).abs().max(1e-30);
            let mut ymin_l = Vec::new();
            let mut ymax_l = Vec::new();
            for f in 0..mesh.n_facets() {
                let fv = mesh.facet(f);
                if fv.iter().all(|&v| (mesh.point(v as usize)[1] - lo[1]).abs() <= tol) {
                    ymin_l.push(mesh.facet_labels[f]);
                } else if fv.iter().all(|&v| (mesh.point(v as usize)[1] - hi[1]).abs() <= tol) {
                    ymax_l.push(mesh.facet_labels[f]);
                }
            }
            ymin_l.sort_unstable();
            ymin_l.dedup();
            ymax_l.sort_unstable();
            ymax_l.dedup();
            (None, Some((ymin_l, ymax_l)))
        } else {
            (None, None)
        };
        let exchange = Mesh::from_solver_mesh(&mesh).map_err(Error::Invalid)?;
        let solver = std::sync::OnceLock::new();
        let _ = solver.set(Arc::new(mesh));
        Ok(MeshedModel {
            mesh: Arc::new(exchange),
            eval,
            labels,
            ends,
            ends_3d,
            report: None,
            amr_adapted: false,
            solver,
        })
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
    let waves = match (s.physics_sizing, sizing) {
        (true, Some(sz)) => Some(WavesRequest {
            sizing: sz.clone(),
            points_per_wavelength: s.points_per_wavelength,
            pml_relax: 1.5,
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
        use femcore::meshing::request::mesh_request;
        let req = build_mesh_request(s, sizing)?;
        if let Some(layers) = &s.sweep {
            if geometry.dim == 3 {
                let _stage = femcore::progress::stage(format!("sweeping 3D geometry along y (element size {} m)", s.h));
                let scale_factor = s.adapt_scale.unwrap_or(1.0);
                if !(scale_factor > 0.0 && scale_factor.is_finite()) {
                    return Err(Error::invalid("the adaptation scale must be positive"));
                }
                let ev = femcore::meshing::geomdrive::evaluate(geometry).map_err(mesh_error)?;
                let labels = ev.label_map().map_err(|e| Error::invalid(e.to_string()))?;
                let ls = resolve_sweep_layers(layers, &labels)?;
                let nr = ev.regions().len();
                let region_h = match (s.physics_sizing, sizing) {
                    (true, Some(sz)) => {
                        let ps = femcore::meshing::physics_size::PhysicsSizing {
                            points_per_wavelength: s.points_per_wavelength,
                            pml_relax: 1.5,
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
                let sizes = femcore::meshing::geomdrive3::Sizes3 { region_h, skip: vec![] };
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
                return MeshedModel::new(mesh, Arc::new(ev), labels, true);
            }
        }
        // the stage a job reports while it meshes (its log gets the mesher's laps below it)
        let _stage = femcore::progress::stage(format!("meshing the geometry (element size {} m)", s.h));
        let (g, report) = mesh_request(geometry, &req).map_err(mesh_error)?;
        let (mut mesh, eval, labels) = (g.mesh, g.eval, g.labels);
        if let Some(layers) = &s.sweep {
            let ls = resolve_sweep_layers(layers, &labels)?;
            // numbered in the order of their coordinates, so that opposite side walls of a periodic cell are cut the same way
            // and a Floquet condition can pair them (`sort_nodes`)
            mesh = femcore::meshing::geomdrive::sort_nodes(&mesh);
            mesh = femcore::meshing::mesher::extrude_y_with_remap(&mesh, &ls, 0.0).map_err(mesh_error)?;
            return MeshedModel::new(mesh, eval, labels, true);
        }
        Ok(MeshedModel::new(mesh, eval, labels, false)?.with_report(report))
    }
}

impl MeshedModel {
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
