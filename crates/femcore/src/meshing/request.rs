// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! A mesh asked for and reported in data (ADR 54): the serialisable request that the engine, the MCP tools and the window can build, run
//! with [`mesh_request`] or zero-eager-split [`mesh_request_native`], and the report that comes back, with JSON Schemas for both
//! (draft 2020-12).

use crate::adapt::drive::AdaptReport;
use crate::mesh::cells::CellType;
use crate::mesh::mixed::MixedDiscretisation;
use crate::mesh::Mesh;
use crate::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
use crate::meshing::physics_size::PhysicsSizing;
use crate::meshing::pml_geom::{mixed_cell_scaled_jacobian_min, PmlMode, PmlReport};
use crate::meshing::{Adaptive, MeshPlan, Mesher, NativeMeshTuple, Sizing};
use geom::eval::Evaluated;
use geom::labels::LabelMap;
use geom::model::FeatureList;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeMap;

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

/// How absorbing layers are meshed in a [`MeshRequest`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PmlModeRequest {
    /// Physics-controlled default: builds structured native `Quad`/`Prism`/`Hex` layers whenever
    /// `!ev.pml_regions().is_empty()`, and does nothing otherwise.
    #[default]
    Auto,
    /// Structured mapped/swept absorbing layers (`Quad` in 2D; `Prism` + `Hex` in 3D).
    Structured,
    /// Legacy unstructured inline simplices inside CAD-sliced absorbing domains.
    Inline,
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
    /// factor on the element size inside inline absorbing layers (at least 1; ignored when structured PML is active)
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
    /// how absorbing layers are meshed (`"auto"` [default], `"structured"`, or `"inline"`)
    #[serde(default)]
    pub pml_mode: PmlModeRequest,
    /// legacy boolean alias (`true` -> `Structured`, `false` -> `Inline` when `pml_mode == Auto`) kept for backward compatibility
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub structured_pml: Option<bool>,
    /// the waves of the study the mesh is for
    #[serde(default)]
    pub waves: Option<WavesRequest>,
    /// 3D / curvilinear shells: quadratic (isoparametric) geometry for the cells on a curved surface (ADR 82)
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub curved: bool,
    /// 3D sheet mesher: target number of cells across narrow coplanar gaps; absent: disabled
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap_cells: Option<f64>,
    /// 3D sheet mesher: conforming pure-tetrahedral vertical coarsening away from sheet planes
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub z_coarsen: bool,
    /// maximum number of cells allowed before aborting meshing; absent: disabled (or `Adaptive::max_cells`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_cells: Option<usize>,
    /// 3D sheet mesher: floor on initial vertical z-slab thickness near sheet planes (m)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dz_min: Option<f64>,
    /// 3D sheet mesher: maximum total number of horizontal z-slabs across active stack layers
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_z_slabs: Option<usize>,
    /// 3D sheet mesher: geometric growth ratio for z-slab thickness away from sheet planes (default 1.35)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub z_growth_ratio: Option<f64>,
    /// 2D Ruppert refinement: cap on additional Steiner vertices as a multiple of initial vertices
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_steiner_multiplier: Option<usize>,
    /// fail meshing if 2D Ruppert refinement terminates early or drops unhealed collinear sliver triangles
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub strict_mesh: bool,
    /// 3D sheet mesher: edge singularity enrichment (`"none"` or `"quarter_point"`); absent: disabled
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub singular_edges: Option<String>,
    /// 3D sheet mesher: Barsoum quarter-point mid-edge parameter `xi_QP` in `[0.25, 0.50]`; absent: `0.25`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub quarter_point_xi: Option<f64>,
    /// 3D sheet mesher: radical corner grading toward sharp 2D sheet corners
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub corner_grading: bool,
}

impl Default for MeshRequest {
    fn default() -> Self {
        MeshRequest {
            h: 0.1,
            max_dev: None,
            mesher: MesherRequest::Simple,
            pml_mode: PmlModeRequest::Auto,
            structured_pml: None,
            waves: None,
            curved: false,
            gap_cells: None,
            z_coarsen: false,
            max_cells: None,
            dz_min: None,
            max_z_slabs: None,
            z_growth_ratio: None,
            max_steiner_multiplier: None,
            strict_mesh: false,
            singular_edges: None,
            quarter_point_xi: None,
            corner_grading: false,
        }
    }
}

impl MeshRequest {
    /// Resolves effective [`PmlModeRequest`], honoring legacy `structured_pml: Some(true | false)` when
    /// `pml_mode` was left at its default `Auto`.
    pub fn resolved_pml_mode(&self) -> PmlModeRequest {
        match (self.pml_mode, self.structured_pml) {
            (PmlModeRequest::Auto, Some(true)) => PmlModeRequest::Structured,
            (PmlModeRequest::Auto, Some(false)) => PmlModeRequest::Inline,
            (mode, _) => mode,
        }
    }

    /// The plan this request names.
    pub fn plan(&self) -> MeshPlan {
        let mut opts = GeomMeshOpts::new(self.h);
        opts.max_dev = self.max_dev;
        opts.curved = self.curved;
        opts.gap_cells = self.gap_cells;
        opts.z_coarsen = self.z_coarsen;
        opts.max_cells = self.max_cells;
        opts.dz_min = self.dz_min;
        opts.max_z_slabs = self.max_z_slabs;
        opts.z_growth_ratio = self.z_growth_ratio;
        opts.max_steiner_multiplier = self.max_steiner_multiplier;
        opts.strict_mesh = self.strict_mesh;
        opts.singular_edges = self.singular_edges.clone();
        opts.quarter_point_xi = self.quarter_point_xi;
        opts.corner_grading = self.corner_grading;
        MeshPlan {
            opts,
            mesher: match self.mesher {
                MesherRequest::Simple => Mesher::Simple,
                MesherRequest::Adaptive { scale, grade, max_cells } => {
                    Mesher::Adaptive(Adaptive { scale, grade, max_cells })
                }
            },
            pml: match self.resolved_pml_mode() {
                PmlModeRequest::Auto => PmlMode::Auto,
                PmlModeRequest::Structured => PmlMode::Structured,
                PmlModeRequest::Inline => PmlMode::Inline,
            },
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
    /// uncut element counts by shape (`"tri"`, `"quad"`, `"tet"`, `"prism"`, `"hex"`)
    #[serde(default)]
    pub element_counts: BTreeMap<String, usize>,
    pub domains: Vec<DomainReport>,
    #[serde(default)]
    pub adaptation: Option<AdaptSummary>,
    /// layers per side of the structured absorbing layers, `[axis][side]` (x, z in 2D; x, y, z in 3D); absent for inline layers
    #[serde(default)]
    pub structured_layers: Option<Vec<[usize; 2]>>,
    /// minimum normalized scaled Jacobian across all generated PML cells (curved and straight); absent when no structured PML cells exist
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_pml_scaled_jacobian: Option<f64>,
    /// absorbing-layer quality and resolution warnings (e.g. `< 5` layers); absent when none
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub pml_warnings: Vec<String>,
    /// the curved edges of a 3D mesh made with `curved` (ADR 82); absent for a straight-sided mesh
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub curved: Option<CurvedSummary>,
    /// things the request asked for that the mesher could not do (a size rule on an edge of a 3D model, say); absent when none
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

/// The curved geometry of a 3D mesh (ADR 82).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurvedSummary {
    /// edges that carry a geometry node
    pub edges: usize,
    /// cells with at least one such edge
    pub cells: usize,
    /// smallest ratio of the Jacobian determinant to the straight cell's over those cells
    pub min_scaled_jacobian: f64,
}

impl CurvedSummary {
    /// Summarizes the curved and/or quarter-point geometry of a 3D mesh, evaluating quarter-point cells
    /// (`m.cell_is_qp(c)`) on the 216 interior Gauss points of `quad_rule_collapsed(8)` and regular curved
    /// cells on the ten Lagrange nodes plus `quad_rule_collapsed(6)`.
    pub fn from_mesh(m: &crate::mesh::SimplexMesh) -> Option<Self> {
        (m.dim == 3 && m.has_curved()).then(|| {
            let cells: Vec<usize> = (0..m.n_cells()).filter(|&c| m.cell_is_curved(c)).collect();
            let min_sj = cells
                .iter()
                .map(|&c| {
                    if m.cell_is_qp(c) {
                        crate::meshing::curve3::scaled_jacobian_min_interior_q8(m, c)
                    } else {
                        crate::meshing::curve3::scaled_jacobian_min(m, c)
                    }
                })
                .fold(f64::INFINITY, f64::min);
            CurvedSummary {
                edges: m.edge_node_list().len(),
                cells: cells.len(),
                min_scaled_jacobian: if min_sj.is_finite() { min_sj } else { 1.0 },
            }
        })
    }
}

fn cell_type_key(t: CellType) -> &'static str {
    match t {
        CellType::Tri => "tri",
        CellType::Quad => "quad",
        CellType::Tet => "tet",
        CellType::Prism => "prism",
        CellType::Hex => "hex",
    }
}

fn collect_pml_warnings(ev: &Evaluated) -> Vec<String> {
    let mut min_layers: Option<usize> = None;
    for p in ev.pml_regions() {
        let n = (p.layers.round() as usize).max(1);
        if n < 5 {
            min_layers = Some(min_layers.map_or(n, |cur| cur.min(n)));
        }
    }
    let mut out = Vec::new();
    if let Some(n) = min_layers {
        out.push(format!(
            "absorbing layer uses {n} (< 5) layer(s); at least 5-6 layers are recommended to suppress discretization reflection"
        ));
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn summarise_native_mesh(
    mesh: &Mesh,
    labels: &LabelMap,
    eval: &Evaluated,
    plan: &MeshPlan,
    adapt: Option<&AdaptReport>,
    pml_rep: Option<&PmlReport>,
    sizes: &[Option<f64>],
    base_warnings: &[String],
) -> MeshReport {
    let (dim, vertices, cells, facets, cell_labels, element_counts) = match mesh {
        Mesh::D2(cm) => {
            let mut counts = BTreeMap::new();
            for (ct, cnt) in cm.counts() {
                if cnt > 0 {
                    counts.insert(cell_type_key(ct).to_string(), cnt);
                }
            }
            (2, cm.n_points(), cm.n_cells(), cm.n_facets(), &cm.label, counts)
        }
        Mesh::D3(cm) => {
            let mut counts = BTreeMap::new();
            for (ct, cnt) in cm.counts() {
                if cnt > 0 {
                    counts.insert(cell_type_key(ct).to_string(), cnt);
                }
            }
            (3, cm.n_points(), cm.n_cells(), cm.n_facets(), &cm.label, counts)
        }
    };

    let domains = labels
        .domains
        .iter()
        .map(|d| DomainReport {
            id: d.id,
            name: d.name.clone(),
            cells: cell_labels.iter().filter(|&&l| l == d.id).count(),
            wave_size: sizes.get(d.region as usize).copied().flatten(),
        })
        .filter(|d| d.cells > 0 || d.wave_size.is_some())
        .collect();

    let adaptation = adapt.and_then(|r| {
        let (first, last) = (r.iters.first()?, r.iters.last()?);
        Some(AdaptSummary {
            iterations: r.iters.len(),
            converged: r.converged,
            capped: r.capped,
            q_min: last.q_min,
            cells_start: first.n_triangles,
            cells_end: last.n_triangles,
        })
    });

    let use_structured = plan.pml.is_structured(!eval.pml_regions().is_empty());
    let structured_layers = if use_structured {
        if let Some(pr) = pml_rep {
            if dim == 2 {
                Some(pr.layers.to_vec())
            } else {
                Some(pr.layers3.to_vec())
            }
        } else {
            let mut layers = vec![[0usize; 2]; dim];
            for p in eval.pml_regions() {
                let nz: Vec<usize> = (0..3).filter(|&a| p.sides[a] != 0).collect();
                if nz.len() == 1 {
                    let axis = if dim == 2 {
                        [0usize, 2].iter().position(|&a| a == nz[0])
                    } else {
                        Some(nz[0])
                    };
                    if let Some(a) = axis {
                        layers[a][(p.sides[nz[0]] > 0) as usize] = (p.layers.round() as usize).max(1);
                    }
                }
            }
            Some(layers)
        }
    } else {
        None
    };

    let min_pml_scaled_jacobian = pml_rep.and_then(|pr| pr.min_pml_scaled_jacobian);
    let pml_warnings = collect_pml_warnings(eval);

    let mut warnings = base_warnings.to_vec();
    if dim == 3 {
        let has_sheets = eval.complex().has_sheets();
        for rule in eval.size_rules() {
            let name = eval
                .list
                .features
                .iter()
                .find(|f| f.id == rule.feature)
                .map_or(String::new(), |f| f.name.clone());
            let what = match rule.entities.kind {
                geom::select::EntityKind::Edges if !has_sheets => Some("crease edges"),
                geom::select::EntityKind::Vertices => Some("vertices"),
                _ => None,
            };
            if let Some(w) = what {
                warnings.push(format!(
                    "size rule '{name}' selects {w}, which a 3D mesh has no size for: it is ignored (select regions or faces)"
                ));
            } else if rule.min_layers.is_some() {
                warnings.push(format!(
                    "size rule '{name}': min_layers has no effect in the 3D lattice mesher and is ignored"
                ));
            } else if rule.entities.kind == geom::select::EntityKind::Regions && rule.max_dev.is_some() {
                warnings.push(format!(
                    "size rule '{name}': max_deviation applies to curved faces, not to regions, and is ignored here (select the faces)"
                ));
            }
        }
    }

    let curved = match mesh {
        Mesh::D3(cm) if !cm.edge_nodes.is_empty() => {
            if let Ok(md) = MixedDiscretisation::from_cell_mesh(cm) {
                let want_qp = !cm.qp_edges.is_empty()
                    || plan.opts.singular_edges.as_deref() == Some("quarter_point");
                let curved_cells: Vec<usize> =
                    (0..cm.n_cells()).filter(|&c| md.cell_is_curved(c)).collect();
                let mut tmp_tet_mesh: Option<crate::mesh::SimplexMesh> = None;
                let min_sj = curved_cells
                    .iter()
                    .map(|&c| {
                        if cm.types[c] == CellType::Tet && (want_qp || md.cell_is_qp(c)) {
                            let sm = tmp_tet_mesh.get_or_insert_with(|| {
                                let pts: Vec<f64> = cm.points.iter().flatten().copied().collect();
                                let mut tet_cells = Vec::new();
                                let mut tet_labels = Vec::new();
                                for ci in 0..cm.n_cells() {
                                    if cm.types[ci] == CellType::Tet {
                                        tet_cells.extend_from_slice(cm.cell(ci));
                                        tet_labels.push(cm.label[ci]);
                                    }
                                }
                                let mut s = crate::mesh::SimplexMesh::new(
                                    3,
                                    pts,
                                    tet_cells,
                                    tet_labels,
                                    Vec::new(),
                                    Vec::new(),
                                );
                                s.set_edge_nodes(&cm.edge_nodes);
                                if !cm.qp_edges.is_empty() {
                                    s.set_qp_edges(&cm.qp_edges);
                                }
                                s
                            });
                            if c < sm.n_cells() && (sm.cell_is_qp(c) || !sm.has_qp()) {
                                crate::meshing::curve3::scaled_jacobian_min_interior_q8(sm, c)
                            } else {
                                mixed_cell_scaled_jacobian_min(cm, &md, c)
                            }
                        } else if cm.types[c] == CellType::Tet {
                            // Exact parity with `curve3::scaled_jacobian_min` (`det J / det_straight` on 10 nodes + q6)
                            let vs = cm.cell(c);
                            let v = [
                                cm.points[vs[0] as usize],
                                cm.points[vs[1] as usize],
                                cm.points[vs[2] as usize],
                                cm.points[vs[3] as usize],
                            ];
                            let mut x = crate::fem::CurvedTet::straight(v).x;
                            const TET_E: [(usize, usize); 6] =
                                [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
                            for (k, &(i, j)) in TET_E.iter().enumerate() {
                                if let Some(q) = cm.edge_node(vs[i], vs[j]) {
                                    x[4 + k] = q;
                                }
                            }
                            let ct = crate::fem::CurvedTet::new(x);
                            let det_s = crate::fem::CurvedTet::straight(v).jacobian(&[0.25; 4]).1;
                            if det_s == 0.0 || !det_s.is_finite() {
                                1.0
                            } else {
                                let rule = crate::fem::quad_rule_collapsed(6);
                                let mut m_min = f64::INFINITY;
                                for l in &rule.lam {
                                    m_min = m_min.min(ct.jacobian(&[l[0], l[1], l[2], l[3]]).1 / det_s);
                                }
                                for (i, j) in TET_E {
                                    let mut l = [0.0; 4];
                                    l[i] = 0.5;
                                    l[j] = 0.5;
                                    m_min = m_min.min(ct.jacobian(&l).1 / det_s);
                                }
                                for i in 0..4 {
                                    let mut l = [0.0; 4];
                                    l[i] = 1.0;
                                    m_min = m_min.min(ct.jacobian(&l).1 / det_s);
                                }
                                m_min
                            }
                        } else {
                            mixed_cell_scaled_jacobian_min(cm, &md, c)
                        }
                    })
                    .fold(f64::INFINITY, f64::min);
                Some(CurvedSummary {
                    edges: md.edge_nodes.len(),
                    cells: curved_cells.len(),
                    min_scaled_jacobian: if min_sj.is_finite() { min_sj } else { 1.0 },
                })
            } else {
                None
            }
        }
        _ => None,
    };

    MeshReport {
        dim,
        vertices,
        cells,
        facets,
        element_counts,
        domains,
        adaptation,
        structured_layers,
        min_pml_scaled_jacobian,
        pml_warnings,
        curved,
        warnings,
    }
}

/// Meshes `list` as `req` asks (any dimension) and reports what came out.
pub fn mesh_request(list: &FeatureList, req: &MeshRequest) -> Result<(MeshedGeom, MeshReport), String> {
    let plan = req.plan();
    let (g, adapt, pml_rep) = match &req.waves {
        Some(w) => {
            let pml_relax = if req.resolved_pml_mode() == PmlModeRequest::Inline {
                w.pml_relax
            } else {
                1.0
            };
            let sz = PhysicsSizing {
                points_per_wavelength: w.points_per_wavelength,
                pml_relax,
                eigen_headroom: w.eigen_headroom,
            };
            crate::meshing::mesh_for_sizing_with_pml_report(list, &plan, &w.sizing, &sz)?
        }
        None => {
            let ev = crate::meshing::geomdrive::evaluate(list)?;
            crate::meshing::mesh_evaluated_with_pml_report(ev, &plan)?
        }
    };
    let sizes: Vec<Option<f64>> = match &req.waves {
        Some(w) => {
            let use_structured = plan.pml.is_structured(!g.eval.pml_regions().is_empty());
            let sz = PhysicsSizing {
                points_per_wavelength: w.points_per_wavelength,
                pml_relax: if use_structured { 1.0 } else { w.pml_relax },
                eigen_headroom: w.eigen_headroom,
            };
            crate::meshing::wave_sizes_by_region(&g.eval, &g.labels, &w.sizing, &sz)?
                .into_iter()
                .map(|h| h.is_finite().then_some(h))
                .collect()
        }
        None => vec![],
    };
    let ex = g.into_exchange_mesh()?;
    let report = summarise_native_mesh(
        &ex,
        &g.labels,
        &g.eval,
        &plan,
        adapt.as_ref(),
        pml_rep.as_ref(),
        &sizes,
        &g.warnings,
    );
    Ok((g, report))
}

/// Meshes `list` as `req` asks directly into an uncut mixed-cell [`Mesh`] (`Quad`, `Prism`, and `Hex`
/// preserved with zero eager simplex splitting).
pub fn mesh_request_native(list: &FeatureList, req: &MeshRequest) -> Result<NativeMeshTuple, String> {
    let plan = req.plan();
    match &req.waves {
        Some(w) => {
            let pml_relax = if req.resolved_pml_mode() == PmlModeRequest::Inline {
                w.pml_relax
            } else {
                1.0
            };
            let sz = PhysicsSizing {
                points_per_wavelength: w.points_per_wavelength,
                pml_relax,
                eigen_headroom: w.eigen_headroom,
            };
            crate::meshing::mesh_for_sizing_native(list, &plan, &w.sizing, &sz)
        }
        None => {
            let ev = crate::meshing::geomdrive::evaluate(list)?;
            crate::meshing::mesh_evaluated_native_report(ev, &plan)
        }
    }
}

/// JSON Schema of [`MeshRequest`].
pub fn request_schema() -> Value {
    let mesher_schema = json!({
        "oneOf": [
            { "type": "object", "required": ["kind"], "additionalProperties": false, "properties": { "kind": { "const": "simple" } } },
            { "type": "object", "required": ["kind"], "additionalProperties": false, "properties": {
                "kind": { "const": "adaptive" },
                "scale": { "type": "number", "exclusiveMinimum": 0, "default": 1.0 },
                "grade": { "type": "number", "exclusiveMinimum": 0, "default": 0.4 },
                "max_cells": { "type": ["integer", "null"], "minimum": 1 }
            } }
        ]
    });
    let waves_schema = json!({
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
    });
    let mut props = serde_json::Map::new();
    props.insert("h".into(), json!({ "type": "number", "exclusiveMinimum": 0, "description": "global element size (m)" }));
    props.insert("max_dev".into(), json!({ "type": ["number", "null"], "exclusiveMinimum": 0, "description": "largest distance of a chord from a curved edge or face (m); default h / 50" }));
    props.insert("mesher".into(), mesher_schema);
    props.insert(
        "pml_mode".into(),
        json!({
            "type": "string",
            "enum": ["auto", "structured", "inline"],
            "default": "auto",
            "description": "how absorbing layers are meshed ('auto' [default], 'structured', or 'inline')"
        }),
    );
    props.insert(
        "structured_pml".into(),
        json!({
            "type": ["boolean", "null"],
            "description": "legacy alias for pml_mode (true -> 'structured', false -> 'inline')"
        }),
    );
    props.insert("curved".into(), json!({ "type": "boolean", "default": false, "description": "3D / curvilinear shells: quadratic (isoparametric) geometry for the cells on a curved surface" }));
    props.insert("gap_cells".into(), json!({ "type": ["number", "null"], "exclusiveMinimum": 0, "description": "3D sheet mesher: target number of cells across narrow coplanar gaps" }));
    props.insert("z_coarsen".into(), json!({ "type": "boolean", "default": false, "description": "3D sheet mesher: conforming pure-tetrahedral vertical coarsening away from sheet planes" }));
    props.insert("max_cells".into(), json!({ "type": ["integer", "null"], "minimum": 1, "description": "maximum number of cells allowed before aborting meshing" }));
    props.insert("dz_min".into(), json!({ "type": ["number", "null"], "exclusiveMinimum": 0, "description": "3D sheet mesher: floor on initial vertical z-slab thickness near sheet planes (m)" }));
    props.insert("max_z_slabs".into(), json!({ "type": ["integer", "null"], "minimum": 1, "description": "3D sheet mesher: maximum total number of horizontal z-slabs across active stack layers" }));
    props.insert("z_growth_ratio".into(), json!({ "type": ["number", "null"], "exclusiveMinimum": 1.0, "description": "3D sheet mesher: geometric growth ratio for z-slab thickness away from sheet planes (default 1.35)" }));
    props.insert("max_steiner_multiplier".into(), json!({ "type": ["integer", "null"], "minimum": 0, "description": "2D Ruppert refinement: cap on additional Steiner vertices as a multiple of initial vertices" }));
    props.insert("strict_mesh".into(), json!({ "type": "boolean", "default": false, "description": "fail meshing if 2D Ruppert refinement terminates early or drops unhealed collinear sliver triangles" }));
    props.insert("singular_edges".into(), json!({ "type": ["string", "null"], "enum": ["none", "quarter_point", null], "description": "3D sheet mesher: edge singularity enrichment ('none' or 'quarter_point')" }));
    props.insert("quarter_point_xi".into(), json!({ "type": ["number", "null"], "minimum": 0.25, "maximum": 0.5, "description": "3D sheet mesher: Barsoum quarter-point mid-edge parameter in [0.25, 0.50] (default 0.25)" }));
    props.insert("corner_grading".into(), json!({ "type": "boolean", "default": false, "description": "3D sheet mesher: radical corner grading toward sharp 2D sheet corners" }));
    props.insert("waves".into(), waves_schema);
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "nyx:mesh_request.v1",
        "title": "MeshRequest",
        "type": "object",
        "required": ["h", "mesher"],
        "additionalProperties": false,
        "properties": Value::Object(props)
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
            "element_counts": {
                "type": "object",
                "additionalProperties": { "type": "integer", "minimum": 0 }
            },
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
            "structured_layers": { "type": ["array", "null"], "items": { "type": "array", "items": { "type": "integer", "minimum": 0 }, "minItems": 2, "maxItems": 2 } },
            "min_pml_scaled_jacobian": { "type": ["number", "null"] },
            "pml_warnings": { "type": "array", "items": { "type": "string" } },
            "curved": { "type": ["object", "null"], "required": ["edges", "cells", "min_scaled_jacobian"], "additionalProperties": false,
                "description": "the curved edges of a 3D mesh made with `curved`; absent for a straight-sided mesh",
                "properties": { "edges": { "type": "integer", "minimum": 0 }, "cells": { "type": "integer", "minimum": 0 }, "min_scaled_jacobian": { "type": "number" } } },
            "warnings": { "type": "array", "items": { "type": "string" }, "description": "requests the mesher could not honour; absent when none" }
        }
    })
}
