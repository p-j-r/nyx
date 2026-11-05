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
    /// the boundary labels of the two end faces of a model swept along y (`ymin`, `ymax`)
    ends: Option<(i32, i32)>,
    /// what the mesher reports (domains, cells per domain, the element size the waves asked for, the adaptation); absent for a swept mesh
    report: Option<femcore::meshing::request::MeshReport>,
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
        let ends = swept.then_some((top - 1, top));
        let exchange = Mesh::from_solver_mesh(&mesh).map_err(Error::Invalid)?;
        let solver = std::sync::OnceLock::new();
        let _ = solver.set(Arc::new(mesh));
        Ok(MeshedModel { mesh: Arc::new(exchange), eval, labels, ends, report: None, solver })
    }

    /// The same, with the report the mesher made of it.
    pub(crate) fn with_report(mut self, report: femcore::meshing::request::MeshReport) -> MeshedModel {
        self.report = Some(report);
        self
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
            "swept_ends": self.ends.map(|_| json!(["ymin", "ymax"])),
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
        match (s, self.ends) {
            ("ymin", Some((a, _))) => return Ok(vec![a]),
            ("ymax", Some((_, b))) => return Ok(vec![b]),
            ("ymin" | "ymax", None) => {
                return Err(Error::invalid(
                    "ymin and ymax are the ends of a model swept along y; this model is not swept",
                ))
            }
            _ => {}
        }
        // a boundary name (`rib.top`, `air|rib`) or a tag first; these are not always valid selector text
        let named = self.labels.boundary(s);
        if !named.is_empty() {
            return Ok(named);
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
}

/// Builds the mesh of a model. `sizing` is what the physics of the project want resolved (`MeshSettings::physics_sizing`); the engine
/// passes it only when that flag is on.
pub trait Mesher: Send + Sync {
    fn mesh(&self, geometry: &FeatureList, settings: &MeshSettings, sizing: Option<&Sizing>) -> Result<MeshedModel, Error>;
}

/// The mesher of `femcore`, driven by a `MeshRequest` (ADR 54): 2D and 3D models, simple or adaptive, inline or structured absorbing
/// layers, element sizes from the waves of the study. A 2D model with `sweep` layers is extruded along y into a 3D mesh afterwards.
pub struct FemcoreMesher;

impl Mesher for FemcoreMesher {
    fn mesh(&self, geometry: &FeatureList, s: &MeshSettings, sizing: Option<&Sizing>) -> Result<MeshedModel, Error> {
        use femcore::meshing::request::{mesh_request, MeshRequest, MesherRequest, WavesRequest};
        if !(s.h > 0.0 && s.h.is_finite()) {
            return Err(Error::invalid("the element size must be a positive finite number"));
        }
        if s.sweep.is_some() && geometry.dim != 2 {
            return Err(Error::invalid("`sweep` extrudes a 2D cross-section; this model is 3D, mesh it directly"));
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
        // wave sizing adapts the mesh to the waves; an explicit `adapt_scale` sets the factor on every size
        let mesher = match (s.adapt_scale, waves.is_some()) {
            (Some(scale), _) => MesherRequest::Adaptive { scale, grade: 0.4, max_cells: Some(2_000_000) },
            (None, true) => MesherRequest::Adaptive { scale: 1.0, grade: 0.4, max_cells: Some(2_000_000) },
            (None, false) => MesherRequest::Simple,
        };
        let req = MeshRequest { h: s.h, max_dev: s.max_dev, mesher, structured_pml: s.structured_pml, waves };
        let (g, report) = mesh_request(geometry, &req).map_err(Error::Invalid)?;
        let (mut mesh, eval, labels) = (g.mesh, g.eval, g.labels);
        if let Some(layers) = &s.sweep {
            let layers: Vec<femcore::meshing::mesher::Layer> = layers
                .iter()
                .map(|l| femcore::meshing::mesher::Layer { thickness: l.thickness, n: l.n, domain_offset: l.domain_offset })
                .collect();
            mesh = femcore::meshing::mesher::extrude_y(&mesh, &layers, 0.0).map_err(Error::Invalid)?;
            return MeshedModel::new(mesh, eval, labels, true);
        }
        Ok(MeshedModel::new(mesh, eval, labels, false)?.with_report(report))
    }
}
