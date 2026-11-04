//! Meshing seam. `MeshedModel` is a mesh together with the names of the model (regions, boundaries) and the label resolver the physics
//! use; `Mesher` builds one from a feature list. `FemcoreMesher` is the **only** code of the engine that calls femcore's mesher
//! (`meshing::mesh_model`, `mesher::extrude_y`): restructure phase 2 replaces the mesh types and moves those modules, and rewrites this
//! one function. Names and labels come from `geom` (`LabelMap`), which numbers domains from 1 in overlay order.

use crate::project::MeshSettings;
use geom::eval::Evaluated;
use geom::labels::LabelMap;
use geom::model::FeatureList;
use geom::select::{EntityKind, Selector};
use nyx_physics::{Error, LabelResolver, Mesh};
use serde_json::{json, Value};
use std::sync::Arc;

pub struct MeshedModel {
    mesh: Arc<Mesh>,
    eval: Arc<Evaluated>,
    labels: LabelMap,
    /// the boundary labels of the two end faces of a model swept along y (`ymin`, `ymax`)
    ends: Option<(i32, i32)>,
}

impl MeshedModel {
    pub fn new(mesh: Mesh, eval: Arc<Evaluated>, labels: LabelMap, swept: bool) -> MeshedModel {
        let top = mesh.facet_labels.iter().copied().max().unwrap_or(0);
        let ends = swept.then_some((top - 1, top));
        MeshedModel { mesh: Arc::new(mesh), eval, labels, ends }
    }

    pub fn mesh(&self) -> &Arc<Mesh> {
        &self.mesh
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
            "dim": self.mesh.dim,
            "regions": self.labels.domains.iter().map(|d| json!({ "label": d.id, "name": d.name })).collect::<Vec<_>>(),
            "boundaries": boundaries,
            "groups": self.labels.groups.iter().map(|g| json!({ "name": g.name, "rule": g.rule })).collect::<Vec<_>>(),
            "swept_ends": self.ends.map(|_| json!(["ymin", "ymax"])),
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

/// Builds the mesh of a model.
pub trait Mesher: Send + Sync {
    fn mesh(&self, geometry: &FeatureList, settings: &MeshSettings) -> Result<MeshedModel, Error>;
}

/// The mesher of `femcore` as it is before restructure phase 2.
pub struct FemcoreMesher;

impl Mesher for FemcoreMesher {
    fn mesh(&self, geometry: &FeatureList, s: &MeshSettings) -> Result<MeshedModel, Error> {
        if geometry.dim != 2 {
            return Err(Error::Unsupported(
                "a 3D model cannot be meshed yet: model the x-z cross-section in 2D and sweep it along y with `sweep` layers".into(),
            ));
        }
        if !(s.h > 0.0 && s.h.is_finite()) {
            return Err(Error::invalid("the element size must be a positive finite number"));
        }
        if s.physics_sizing {
            return Err(Error::Unsupported(
                "physics-sized meshes wait for restructure phase 2 (the neutral Sizing); switch physics_sizing off"
                    .into(),
            ));
        }
        let plan = match s.adapt_scale {
            Some(scale) => femcore::meshing::MeshPlan::adaptive(s.h, scale),
            None => femcore::meshing::MeshPlan::simple(s.h),
        };
        let g = femcore::meshing::mesh_model(geometry, &plan).map_err(Error::Invalid)?;
        let (mut mesh, eval, labels) = (g.mesh, g.eval, g.labels);
        if let Some(layers) = &s.sweep {
            let layers: Vec<femcore::meshing::mesher::Layer> = layers
                .iter()
                .map(|l| femcore::meshing::mesher::Layer { thickness: l.thickness, n: l.n, domain_offset: l.domain_offset })
                .collect();
            mesh = femcore::meshing::mesher::extrude_y(&mesh, &layers, 0.0).map_err(Error::Invalid)?;
        }
        Ok(MeshedModel::new(mesh, eval, labels, s.sweep.is_some()))
    }
}
