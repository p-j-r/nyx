//! A mesh asked for and reported in data (ADR 54): the serialisable request that the engine, the MCP tools and the window can build, run
//! with [`mesh_request`], and the report that comes back, with JSON Schemas for both (draft 2020-12). Nothing here is new meshing: a request
//! names a `MeshPlan` (simple or adaptive, inline or structured absorbing layers), optionally the waves of a study (`Sizing`, points per
//! wavelength), and runs `mesh_for_sizing` or `mesh_model`, in 2D and 3D alike.

use crate::adapt::drive::AdaptReport;
use crate::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
use crate::meshing::physics_size::PhysicsSizing;
use crate::meshing::pml_geom::PmlMode;
use crate::meshing::{mesh_for_sizing, Adaptive, MeshPlan, Mesher, Sizing};
use geom::model::FeatureList;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

/// Which mesher makes the cells.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum MesherRequest {
    /// the lattice (3D) or constrained-Delaunay (2D) mesher at the global size
    Simple,
    /// the metric adaptation towards the sizes of the model and of the waves
    Adaptive {
        /// factor on every size (1 keeps them, 0.5 makes the mesh about 4 times (2D) or 8 times (3D) as fine)
        #[serde(default = "one")]
        scale: f64,
        /// largest growth of the element size per unit distance
        #[serde(default = "grade")]
        grade: f64,
        #[serde(default)]
        max_cells: Option<usize>,
    },
}

fn one() -> f64 {
    1.0
}
fn grade() -> f64 {
    0.4
}
fn relax() -> f64 {
    1.5
}
fn headroom() -> f64 {
    1.5
}

/// The waves of a study, as the physics reports them (`Sizing`) and how finely they are resolved.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WavesRequest {
    pub sizing: Sizing,
    /// points per shortest wavelength; absent: 12 for order 1, 6 for order 2 and above
    #[serde(default)]
    pub points_per_wavelength: Option<f64>,
    /// factor on the element size inside absorbing layers (at least 1)
    #[serde(default = "relax")]
    pub pml_relax: f64,
    /// factor on the search frequency of an eigenfrequency study
    #[serde(default = "headroom")]
    pub eigen_headroom: f64,
}

/// A mesh asked for.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MeshRequest {
    /// global element size (m)
    pub h: f64,
    /// largest distance of a chord from a curved edge or face (m); absent: h / 50
    #[serde(default)]
    pub max_dev: Option<f64>,
    pub mesher: MesherRequest,
    /// `true`: the absorbing layers of the model are rebuilt as structured layers (quadrilaterals, prisms, hexahedra, cut into simplices);
    /// `false`: they are meshed inline like any region
    #[serde(default)]
    pub structured_pml: bool,
    /// the waves of the study the mesh is for
    #[serde(default)]
    pub waves: Option<WavesRequest>,
}

impl MeshRequest {
    /// The plan this request names.
    pub fn plan(&self) -> MeshPlan {
        let mut opts = GeomMeshOpts::new(self.h);
        opts.max_dev = self.max_dev;
        MeshPlan {
            opts,
            mesher: match self.mesher {
                MesherRequest::Simple => Mesher::Simple,
                MesherRequest::Adaptive { scale, grade, max_cells } => Mesher::Adaptive(Adaptive { scale, grade, max_cells }),
            },
            pml: if self.structured_pml { PmlMode::Structured } else { PmlMode::Inline },
        }
    }
}

/// One domain of the meshed model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DomainReport {
    pub id: i32,
    pub name: String,
    pub cells: usize,
    /// the element size the waves asked for in this domain (m); absent where no wave was given
    #[serde(default)]
    pub wave_size: Option<f64>,
}

/// What the adaptation did.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AdaptSummary {
    pub iterations: usize,
    pub converged: bool,
    pub capped: bool,
    /// worst element quality at the end (0 to 1)
    pub q_min: f64,
    pub cells_start: usize,
    pub cells_end: usize,
}

/// A mesh reported.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MeshReport {
    pub dim: usize,
    pub vertices: usize,
    pub cells: usize,
    /// labelled boundary facets (3D: also the faces between regions)
    pub facets: usize,
    pub domains: Vec<DomainReport>,
    #[serde(default)]
    pub adaptation: Option<AdaptSummary>,
    /// layers per side of the structured absorbing layers, `[axis][side]` (x, y, z); absent for inline layers
    #[serde(default)]
    pub structured_layers: Option<Vec<[usize; 2]>>,
}

fn summarise(g: &MeshedGeom, req: &MeshRequest, report: Option<&AdaptReport>, sizes: &[Option<f64>]) -> MeshReport {
    let m = &g.mesh;
    let domains = g
        .labels
        .domains
        .iter()
        .map(|d| DomainReport { id: d.id, name: d.name.clone(), cells: m.cell_labels.iter().filter(|&&l| l == d.id).count(), wave_size: sizes.get(d.region as usize).copied().flatten() })
        .filter(|d| d.cells > 0 || d.wave_size.is_some())
        .collect();
    let adaptation = report.and_then(|r| {
        let (first, last) = (r.iters.first()?, r.iters.last()?);
        Some(AdaptSummary { iterations: r.iters.len(), converged: r.converged, capped: r.capped, q_min: last.q_min, cells_start: first.n_triangles, cells_end: last.n_triangles })
    });
    let structured_layers = if req.structured_pml && !g.eval.pml_regions().is_empty() {
        let mut layers = vec![[0usize; 2]; m.dim];
        for p in g.eval.pml_regions() {
            let nz: Vec<usize> = (0..3).filter(|&a| p.sides[a] != 0).collect();
            if nz.len() == 1 {
                let axis = if m.dim == 2 { [0usize, 2].iter().position(|&a| a == nz[0]) } else { Some(nz[0]) };
                if let Some(a) = axis {
                    layers[a][(p.sides[nz[0]] > 0) as usize] = (p.layers.round() as usize).max(1);
                }
            }
        }
        Some(layers)
    } else {
        None
    };
    MeshReport { dim: m.dim, vertices: m.n_points(), cells: m.n_cells(), facets: m.n_facets(), domains, adaptation, structured_layers }
}

/// Meshes `list` as `req` asks (any dimension) and reports what came out.
pub fn mesh_request(list: &FeatureList, req: &MeshRequest) -> Result<(MeshedGeom, MeshReport), String> {
    let plan = req.plan();
    let (g, adapt) = match &req.waves {
        Some(w) => {
            let sz = PhysicsSizing { points_per_wavelength: w.points_per_wavelength, pml_relax: w.pml_relax, eigen_headroom: w.eigen_headroom };
            mesh_for_sizing(list, &plan, &w.sizing, &sz)?
        }
        None => {
            let ev = crate::meshing::geomdrive::evaluate(list)?;
            crate::meshing::mesh_evaluated_report(ev, &plan)?
        }
    };
    let sizes: Vec<Option<f64>> = match &req.waves {
        Some(w) => {
            let sz = PhysicsSizing { points_per_wavelength: w.points_per_wavelength, pml_relax: w.pml_relax, eigen_headroom: w.eigen_headroom };
            crate::meshing::wave_sizes_by_region(&g.eval, &g.labels, &w.sizing, &sz)?.into_iter().map(|h| h.is_finite().then_some(h)).collect()
        }
        None => vec![],
    };
    let report = summarise(&g, req, adapt.as_ref(), &sizes);
    Ok((g, report))
}

/// JSON Schema of [`MeshRequest`].
pub fn request_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "nyx:mesh_request.v1",
        "title": "MeshRequest",
        "type": "object",
        "required": ["h", "mesher"],
        "additionalProperties": false,
        "properties": {
            "h": { "type": "number", "exclusiveMinimum": 0, "description": "global element size (m)" },
            "max_dev": { "type": ["number", "null"], "exclusiveMinimum": 0, "description": "largest distance of a chord from a curved edge or face (m); default h / 50" },
            "mesher": {
                "oneOf": [
                    { "type": "object", "required": ["kind"], "additionalProperties": false, "properties": { "kind": { "const": "simple" } } },
                    { "type": "object", "required": ["kind"], "additionalProperties": false, "properties": {
                        "kind": { "const": "adaptive" },
                        "scale": { "type": "number", "exclusiveMinimum": 0, "default": 1.0 },
                        "grade": { "type": "number", "exclusiveMinimum": 0, "default": 0.4 },
                        "max_cells": { "type": ["integer", "null"], "minimum": 1 }
                    } }
                ]
            },
            "structured_pml": { "type": "boolean", "default": false },
            "waves": {
                "type": ["object", "null"],
                "required": ["sizing"],
                "additionalProperties": false,
                "properties": {
                    "sizing": {
                        "type": "object",
                        "required": ["order", "headroom", "wavelengths"],
                        "additionalProperties": false,
                        "properties": {
                            "order": { "type": "integer", "minimum": 1 },
                            "headroom": { "type": "number", "exclusiveMinimum": 0 },
                            "wavelengths": { "type": "object", "additionalProperties": { "type": "number", "exclusiveMinimum": 0 } }
                        }
                    },
                    "points_per_wavelength": { "type": ["number", "null"], "exclusiveMinimum": 0 },
                    "pml_relax": { "type": "number", "minimum": 1, "default": 1.5 },
                    "eigen_headroom": { "type": "number", "exclusiveMinimum": 0, "default": 1.5 }
                }
            }
        }
    })
}

/// JSON Schema of [`MeshReport`].
pub fn report_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "nyx:mesh_report.v1",
        "title": "MeshReport",
        "type": "object",
        "required": ["dim", "vertices", "cells", "facets", "domains"],
        "additionalProperties": false,
        "properties": {
            "dim": { "enum": [2, 3] },
            "vertices": { "type": "integer", "minimum": 0 },
            "cells": { "type": "integer", "minimum": 0 },
            "facets": { "type": "integer", "minimum": 0 },
            "domains": { "type": "array", "items": {
                "type": "object", "required": ["id", "name", "cells"], "additionalProperties": false,
                "properties": {
                    "id": { "type": "integer", "minimum": 1 },
                    "name": { "type": "string" },
                    "cells": { "type": "integer", "minimum": 0 },
                    "wave_size": { "type": ["number", "null"] }
                } } },
            "adaptation": { "type": ["object", "null"], "required": ["iterations", "converged", "capped", "q_min", "cells_start", "cells_end"], "additionalProperties": false,
                "properties": {
                    "iterations": { "type": "integer", "minimum": 0 },
                    "converged": { "type": "boolean" },
                    "capped": { "type": "boolean" },
                    "q_min": { "type": "number" },
                    "cells_start": { "type": "integer", "minimum": 0 },
                    "cells_end": { "type": "integer", "minimum": 0 }
                } },
            "structured_layers": { "type": ["array", "null"], "items": { "type": "array", "items": { "type": "integer", "minimum": 0 }, "minItems": 2, "maxItems": 2 } }
        }
    })
}
