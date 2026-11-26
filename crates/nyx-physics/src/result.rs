// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What a solve returns: plain, serialisable data (flat arrays), so a result can be stored, sent to another process or compared with a
//! golden file. It does not carry an analysis: analyses are functions of a result (`nyx-analysis`).

use serde::{Deserialize, Serialize};

/// A field with one or more components, one value per dof (P1: vertices; P2: vertices then edge midpoints).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Field {
    pub name: String,
    pub unit: String,
    /// wafer-axis names of the components (`x`, `y`, `z`); empty for a scalar
    pub components: Vec<String>,
    /// `n_dofs` values for a scalar, `n_dofs * components.len()` for a vector, interleaved by dof
    pub values: Vec<f64>,
}

/// A field before it is averaged over the cells at a node, one value per cell node (`nloc` per cell).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CellField {
    pub name: String,
    pub unit: String,
    pub nloc: usize,
    pub re: Vec<f64>,
    /// empty for a real field
    pub im: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Scalar {
    Num(f64),
    Text(String),
}

/// One solution of a study: a time step, an eigenmode, a frequency.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Solution {
    pub label: String,
    pub fields: Vec<Field>,
    /// the real value that orders the solutions: time (s), frequency (Hz) or mode number
    pub value: Option<f64>,
    /// a complex eigenvalue in the study's own unit (effective index, complex frequency), as (re, im)
    pub eigenvalue: Option<[f64; 2]>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cell_fields: Vec<CellField>,
}

impl Solution {
    pub fn field(&self, name: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.name == name)
    }

    /// The values of a field by its flat name (`V`, `|E|`, `Ex`; see `flat_fields`).
    pub fn flat(&self, name: &str) -> Option<Vec<f64>> {
        self.flat_fields().into_iter().find(|(n, _, _)| n == name).map(|(_, _, v)| v)
    }
}

/// Severity level for structured VVUQ solver, material, or discretization warnings.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiagnosticSeverity {
    #[default]
    Info,
    Warning,
    Error,
}

/// Structured machine-readable diagnostic warning attached to [`SolveResult`].
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct DiagnosticWarning {
    pub code: String,
    pub severity: DiagnosticSeverity,
    pub message: String,
}

/// Algebraic solve, condition, and convergence telemetry emitted by every solver (`solve_result` `2.1`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SolverDiagnostics {
    pub solver_kind: String,
    pub converged: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub iterations: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub restarts: Option<usize>,
    pub max_algebraic_residual: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub per_mode_residuals: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rcond_estimate: Option<f64>,
    #[serde(default)]
    pub cholesky_fallbacks: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<DiagnosticWarning>,
}

/// Spatial discretization uncertainty estimate (Prager–Synge bound, 3-grid Roache GCI, 2-grid delta, or indicator).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SpatialErrorEstimate {
    pub method: String,
    pub rel_bound: f64,
    pub is_epistemic_bound: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub observed_order: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub refinement_ratio: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub safety_factor: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub extrapolated_scalars_only: Option<bool>,
}

/// Numerical error budget pairing algebraic floor, discretization uncertainty, temporal error, and model truncation.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct NumericalErrorBudget {
    pub algebraic_rel: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spatial_discretization: Option<SpatialErrorEstimate>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temporal_rel: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub model_truncation_rel: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub combined_num_uncertainty_rel: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub alg_below_disc_ratio: Option<f64>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SolveResult {
    /// the code of the physics that produced it and the study kind
    pub physics: String,
    pub study: String,
    pub order: u8,
    /// 2 or 3
    pub dim: u8,
    pub dofs: usize,
    pub seconds: f64,
    pub solutions: Vec<Solution>,
    pub scalars: Vec<(String, Scalar)>,
    /// the plain-language stage log of the run
    pub log: Vec<String>,
    /// algebraic solver diagnostics and warnings (additive `solve_result` `2.1`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diagnostics: Option<SolverDiagnostics>,
    /// decomposed numerical error budget (additive `solve_result` `2.1`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error_budget: Option<NumericalErrorBudget>,
}

impl SolveResult {
    /// The first numeric scalar whose name starts with `name`.
    pub fn scalar(&self, name: &str) -> Option<f64> {
        self.scalars.iter().find(|(n, _)| n.starts_with(name)).and_then(|(_, v)| match v {
            Scalar::Num(x) => Some(*x),
            Scalar::Text(_) => None,
        })
    }
}

impl Field {
    /// A scalar field, one value per dof.
    pub fn scalar(name: impl Into<String>, unit: impl Into<String>, values: Vec<f64>) -> Field {
        Field { name: name.into(), unit: unit.into(), components: vec![], values }
    }

    /// A vector field from one array per component (`comps[j][dof]`), stored interleaved by dof.
    pub fn vector(name: impl Into<String>, unit: impl Into<String>, names: &[&str], comps: &[Vec<f64>]) -> Field {
        assert_eq!(names.len(), comps.len(), "one array per component name");
        let n = comps.first().map_or(0, Vec::len);
        let mut values = Vec::with_capacity(n * comps.len());
        for i in 0..n {
            for c in comps {
                values.push(c[i]);
            }
        }
        Field {
            name: name.into(),
            unit: unit.into(),
            components: names.iter().map(|s| s.to_string()).collect(),
            values,
        }
    }

    /// The array of component `j` of a vector field (the whole array for a scalar).
    pub fn component(&self, j: usize) -> Vec<f64> {
        let m = self.components.len().max(1);
        self.values.chunks(m).map(|p| p[j]).collect()
    }
}

impl Solution {
    /// The fields in the flat form the window and the golden files use: a scalar keeps its name; a vector `E` with components `x`, `z`
    /// becomes `|E|` (the norm), `Ex`, `Ez`, in that order. This is the naming of the pre-restructure results (ADR 41).
    pub fn flat_fields(&self) -> Vec<(String, String, Vec<f64>)> {
        let mut out = vec![];
        for f in &self.fields {
            if f.components.is_empty() {
                out.push((f.name.clone(), f.unit.clone(), f.values.clone()));
            } else {
                let comps: Vec<Vec<f64>> = (0..f.components.len()).map(|j| f.component(j)).collect();
                let n = comps[0].len();
                let norm: Vec<f64> = (0..n).map(|i| comps.iter().map(|c| c[i] * c[i]).sum::<f64>().sqrt()).collect();
                out.push((format!("|{}|", f.name), f.unit.clone(), norm));
                for (name, c) in f.components.iter().zip(comps) {
                    out.push((format!("{}{}", f.name, name), f.unit.clone(), c));
                }
            }
        }
        out
    }
}

/// The per-cell fields of a kernel (`femcore::fields::CellFields`) as plain serialisable data.
pub fn cell_fields_of(cf: &femcore::fields::CellFields) -> Vec<CellField> {
    cf.fields
        .iter()
        .map(|f| CellField {
            name: f.name.clone(),
            unit: f.unit.clone(),
            nloc: cf.nloc,
            re: f.re.clone(),
            im: f.im.clone(),
        })
        .collect()
}

/// The inverse of `cell_fields_of`, for code that still reads `CellFields` (the window until phase 7, the analysis until phase 5).
pub fn cell_fields_back(order: usize, fields: &[CellField]) -> femcore::fields::CellFields {
    let nloc = fields.first().map_or(0, |f| f.nloc);
    femcore::fields::CellFields {
        order,
        nloc,
        fields: fields
            .iter()
            .map(|f| femcore::fields::CellField {
                name: f.name.clone(),
                unit: f.unit.clone(),
                re: f.re.clone(),
                im: f.im.clone(),
            })
            .collect(),
    }
}

impl femcore::adapt::PlotFields for SolveResult {
    fn plot_fields(&self) -> Vec<Vec<(String, Vec<f64>)>> {
        self.solutions.iter().map(|s| s.flat_fields().into_iter().map(|(n, _, v)| (n, v)).collect()).collect()
    }
}

// ------------------------------------------------------------------ Stage-3 MixedDiscretisation -> SimplexMesh projector

use crate::Error;
use femcore::fem_mixed::{cell_corner_ref_coords, lagrange_shape_into, LagrangeSpace, MixedLagrangeDofMap};
use femcore::mesh::cells::CellType;
use femcore::mesh::mixed::{cell_local_quad_faces, MixedDiscretisation};
use femcore::mesh::SimplexMesh;

fn minimal_sub_simplices_local(ct: CellType, vs: &[u32]) -> Vec<Vec<u8>> {
    let quad_tris_local = |cyc: [u8; 4]| -> [[u8; 3]; 2] {
        let m = (0..4).min_by_key(|&i| vs[cyc[i] as usize]).unwrap();
        let c = |k: usize| cyc[(m + k) % 4];
        [[c(0), c(1), c(2)], [c(0), c(2), c(3)]]
    };
    match ct {
        CellType::Tri => vec![vec![0, 1, 2]],
        CellType::Tet => vec![vec![0, 1, 2, 3]],
        CellType::Quad => quad_tris_local([0, 1, 2, 3]).into_iter().map(|t| t.to_vec()).collect(),
        CellType::Prism => {
            let v0 = (0..6).min_by_key(|&k| vs[k]).unwrap() as u8;
            const PRISM_QUADS: [[u8; 4]; 3] = [[0, 1, 3, 4], [1, 2, 4, 5], [0, 2, 3, 5]];
            let mut out = Vec::with_capacity(3);
            for f in &PRISM_QUADS {
                if !f.contains(&v0) {
                    for tr in quad_tris_local([f[0], f[1], f[3], f[2]]) {
                        out.push(vec![v0, tr[0], tr[1], tr[2]]);
                    }
                }
            }
            for tri in [[0u8, 1, 2], [3, 4, 5]] {
                if !tri.contains(&v0) {
                    out.push(vec![v0, tri[0], tri[1], tri[2]]);
                }
            }
            out
        }
        CellType::Hex => {
            // `SimplexMesh::from_cells_split_with` indexes `r = [v0, v1, v3, v2, v4, v5, v7, v6]`
            let r_map: [u8; 8] = [0, 1, 3, 2, 4, 5, 7, 6];
            let v0_r = (0..8).min_by_key(|&k| vs[r_map[k] as usize]).unwrap();
            let v0_loc = r_map[v0_r];
            const HEX_FACES: [[usize; 4]; 6] =
                [[0, 2, 4, 6], [1, 3, 5, 7], [0, 1, 4, 5], [2, 3, 6, 7], [0, 1, 2, 3], [4, 5, 6, 7]];
            let mut out = Vec::with_capacity(6);
            for f in &HEX_FACES {
                if !f.contains(&v0_r) {
                    let cyc = [r_map[f[0]], r_map[f[1]], r_map[f[3]], r_map[f[2]]];
                    for tr in quad_tris_local(cyc) {
                        out.push(vec![v0_loc, tr[0], tr[1], tr[2]]);
                    }
                }
            }
            out
        }
    }
}

/// Stage-3 projector mapping native [`MixedDiscretisation`] vertex, edge, and per-cell fields onto
/// a `SplitMode::Minimal` [`SimplexMesh`] (`Mesh::discretise()`) for `Display`, `h5io`, `xdmf`, and `view3d`.
#[derive(Clone, Debug)]
pub struct MixedToSimplexProjector {
    dim: usize,
    n_points: usize,
    pure_simplex: bool,
    n_simplex_cells: usize,
    n_simplex_edges: usize,
    /// For each `SimplexMesh` cell `sc`, the parent `MixedDiscretisation` cell index `mc`.
    parent_cells: Vec<usize>,
    /// For each `SimplexMesh` cell `sc`, the local vertex index (`0..mc.nv()`) in `mc` of each corner of `sc`.
    sub_local_corners: Vec<[u8; 4]>,
    /// For each `SimplexMesh` edge `se`, either `Some(mixed_edge_id)` when the edge exists natively in `mixed`,
    /// or `None` when `se` is a split face-diagonal or body-diagonal edge.
    simplex_edge_to_mixed: Vec<Option<usize>>,
    /// For split diagonal `SimplexMesh` edges (`simplex_edge_to_mixed[se] == None`), one `(parent_mc, midpoint_xi_ref)`
    /// pair evaluating the continuous field at the exact geometric/reference midpoint of the diagonal.
    diagonal_edge_eval: Vec<Option<(usize, [f64; 3])>>,
}

impl MixedToSimplexProjector {
    /// Builds a Stage-3 projector from `mixed` (`discretise_mixed(mesh)`) onto `simplex` (`discretise(mesh)`).
    pub fn new(mixed: &MixedDiscretisation, simplex: &SimplexMesh) -> Result<Self, Error> {
        if mixed.dim != simplex.dim {
            return Err(Error::invalid(format!(
                "projector dimension mismatch: mixed dim {} vs simplex dim {}",
                mixed.dim, simplex.dim
            )));
        }
        if mixed.n_points() != simplex.n_points() {
            return Err(Error::invalid(format!(
                "projector expects SplitMode::Minimal with identical vertex count: mixed {} vs simplex {}",
                mixed.n_points(),
                simplex.n_points()
            )));
        }
        let d = mixed.dim;
        let np = mixed.n_points();
        let nsc = simplex.n_cells();
        let pure_simplex = mixed.is_pure_simplex() && mixed.n_cells() == nsc;

        // Build vertex -> mixed cells adjacency list.
        let mut v_cells: Vec<Vec<u32>> = vec![Vec::new(); np];
        for mc in 0..mixed.n_cells() {
            for &v in mixed.cell(mc) {
                v_cells[v as usize].push(mc as u32);
            }
        }

        let mut parent_cells = Vec::with_capacity(nsc);
        let mut sub_local_corners = Vec::with_capacity(nsc);

        for sc in 0..nsc {
            let sv = simplex.cell(sc);
            let slab = simplex.cell_labels[sc];
            let mut chosen: Option<(usize, [u8; 4])> = None;

            for &mc_u32 in &v_cells[sv[0] as usize] {
                let mc = mc_u32 as usize;
                if mixed.cell_labels[mc] != slab {
                    continue;
                }
                let mvs = mixed.cell(mc);
                let mut lc = [0u8; 4];
                let mut all_in = true;
                for (k, &v_id) in sv.iter().enumerate() {
                    if let Some(pos) = mvs.iter().position(|&x| x == v_id) {
                        lc[k] = pos as u8;
                    } else {
                        all_in = false;
                        break;
                    }
                }
                if !all_in {
                    continue;
                }
                // In 3D, reject a candidate cell if all 4 vertices of `sv` lie on one single local quad face (coplanar).
                if d == 3 {
                    let coplanar_on_quad = cell_local_quad_faces(mixed.cell_type(mc))
                        .iter()
                        .any(|qf| (0..4).all(|k| qf.contains(&(lc[k] as usize))));
                    if coplanar_on_quad {
                        continue;
                    }
                }
                // Verify `lc[..d+1]` (up to permutation) matches one of `mc`'s SplitMode::Minimal sub-simplices.
                let mut sorted_lc = lc[..d + 1].to_vec();
                sorted_lc.sort_unstable();
                let matches_sub = minimal_sub_simplices_local(mixed.cell_type(mc), mvs).into_iter().any(|mut sub| {
                    sub.sort_unstable();
                    sub == sorted_lc
                });
                if matches_sub {
                    chosen = Some((mc, lc));
                    break;
                }
                if chosen.is_none() {
                    chosen = Some((mc, lc));
                }
            }

            let (mc, lc) = chosen.ok_or_else(|| {
                Error::invalid(format!("simplex cell {sc} ({sv:?}) has no matching parent cell in MixedDiscretisation"))
            })?;
            parent_cells.push(mc);
            sub_local_corners.push(lc);
        }

        let nse = simplex.n_edges();
        let mut simplex_edge_to_mixed = Vec::with_capacity(nse);
        let mut diagonal_edge_eval = vec![None; nse];
        for &(a, b) in simplex.edges() {
            let me = mixed.edge_index(a, b);
            simplex_edge_to_mixed.push(me);
        }

        // Populate exact reference midpoints for any split diagonal edges (`simplex_edge_to_mixed[se] == None`).
        if simplex_edge_to_mixed.iter().any(|x| x.is_none()) {
            let npair = d * (d + 1) / 2;
            for sc in 0..nsc {
                let eids = simplex.cell_edge_ids(sc);
                debug_assert_eq!(eids.len(), npair);
                let mc = parent_cells[sc];
                let ct = mixed.cell_type(mc);
                let corners = cell_corner_ref_coords(ct);
                let lc = sub_local_corners[sc];
                let mut k = 0usize;
                for i in 0..=d {
                    for j in i + 1..=d {
                        let se = eids[k] as usize;
                        if simplex_edge_to_mixed[se].is_none() && diagonal_edge_eval[se].is_none() {
                            let (ci, cj) = (corners[lc[i] as usize], corners[lc[j] as usize]);
                            let mid = [0.5 * (ci[0] + cj[0]), 0.5 * (ci[1] + cj[1]), 0.5 * (ci[2] + cj[2])];
                            diagonal_edge_eval[se] = Some((mc, mid));
                        }
                        k += 1;
                    }
                }
            }
        }

        Ok(Self {
            dim: d,
            n_points: np,
            pure_simplex,
            n_simplex_cells: nsc,
            n_simplex_edges: nse,
            parent_cells,
            sub_local_corners,
            simplex_edge_to_mixed,
            diagonal_edge_eval,
        })
    }

    /// Convenience constructor building both `SimplexMesh` and `MixedToSimplexProjector` from `mesh`.
    pub fn from_mesh(mesh: &crate::Mesh, mixed: &MixedDiscretisation) -> Result<(SimplexMesh, Self), Error> {
        let simplex = crate::discretise(mesh)?;
        let proj = Self::new(mixed, &simplex)?;
        Ok((simplex, proj))
    }

    /// Whether the underlying mesh is pure-simplex with 1-to-1 identity cell correspondence.
    #[inline]
    pub fn is_pure_simplex(&self) -> bool {
        self.pure_simplex
    }

    /// Number of cells in the target `SimplexMesh`.
    #[inline]
    pub fn n_simplex_cells(&self) -> usize {
        self.n_simplex_cells
    }

    /// Parent cell index in `MixedDiscretisation` of `SimplexMesh` cell `sc`.
    #[inline]
    pub fn parent_cell(&self, sc: usize) -> usize {
        self.parent_cells[sc]
    }

    /// Reference coordinates inside parent cell `parent_cell(sc)` for all `nloc(target_order)` nodes of `sc`
    /// (`corners` first, then local simplex edge midpoints `(i < j)` when `target_order == 2`).
    pub fn sub_simplex_ref_points(&self, mixed: &MixedDiscretisation, sc: usize, target_order: usize) -> Vec<[f64; 3]> {
        let d = self.dim;
        let mc = self.parent_cells[sc];
        let corners = cell_corner_ref_coords(mixed.cell_type(mc));
        let lc = self.sub_local_corners[sc];
        let mut pts = Vec::with_capacity(if target_order == 1 { d + 1 } else { (d + 1) * (d + 2) / 2 });
        for k in 0..=d {
            pts.push(corners[lc[k] as usize]);
        }
        if target_order == 2 {
            for i in 0..=d {
                for j in i + 1..=d {
                    let (a, b) = (pts[i], pts[j]);
                    pts.push([0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])]);
                }
            }
        }
        pts
    }

    /// Projects `ncomp` stacked scalar DOF blocks (`u_blocks[comp * dof_map.n_dofs + dof]`) computed on
    /// `MixedDiscretisation` onto `SimplexMesh` DOF blocks (`out[comp * simplex.n_dofs(target_order) + sdof]`).
    pub fn project_lagrange_blocks(
        &self,
        mixed: &MixedDiscretisation,
        dof_map: &MixedLagrangeDofMap,
        u_blocks: &[f64],
        ncomp: usize,
        target_order: usize,
    ) -> Vec<f64> {
        assert!(target_order == 1 || target_order == 2);
        let n_in = dof_map.n_dofs;
        assert_eq!(u_blocks.len(), ncomp * n_in);
        let nv = self.n_points;
        let n_out = nv + if target_order == 2 { self.n_simplex_edges } else { 0 };

        if self.pure_simplex && target_order == 1 && dof_map.space == LagrangeSpace::Order1 {
            return u_blocks.to_vec();
        }
        if self.pure_simplex && target_order == 2 && dof_map.space == LagrangeSpace::Order2Serendipity {
            return u_blocks.to_vec();
        }

        let mut out = vec![0.0; ncomp * n_out];
        // 1. Copy vertex DOFs directly (`0..n_points` are identical!).
        for comp in 0..ncomp {
            out[comp * n_out..comp * n_out + nv].copy_from_slice(&u_blocks[comp * n_in..comp * n_in + nv]);
        }
        if target_order == 1 {
            return out;
        }

        // 2. Populate edge midpoint DOFs (`nv .. nv + n_simplex_edges`).
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        for se in 0..self.n_simplex_edges {
            if dof_map.space != LagrangeSpace::Order1 {
                if let Some(me) = self.simplex_edge_to_mixed[se] {
                    for comp in 0..ncomp {
                        out[comp * n_out + nv + se] = u_blocks[comp * n_in + nv + me];
                    }
                    continue;
                }
            }
            if let Some((mc, xi_mid)) = self.diagonal_edge_eval[se] {
                let ct = mixed.cell_type(mc);
                lagrange_shape_into(ct, dof_map.space, &xi_mid, &mut vals, &mut dxi);
                let cdofs = dof_map.cell_dofs(mc);
                for comp in 0..ncomp {
                    let mut sum = 0.0;
                    for (a, &gd) in cdofs.iter().enumerate() {
                        sum += vals[a] * u_blocks[comp * n_in + gd];
                    }
                    out[comp * n_out + nv + se] = sum;
                }
            } else if let Some(me) = self.simplex_edge_to_mixed[se] {
                // `Order1` mixed field promoted to `target_order == 2` simplex display: average endpoints!
                let (va, vb) = mixed.edges[me];
                for comp in 0..ncomp {
                    out[comp * n_out + nv + se] =
                        0.5 * (u_blocks[comp * n_in + va as usize] + u_blocks[comp * n_in + vb as usize]);
                }
            }
        }
        out
    }

    /// Projects a scalar `MixedLagrangeDofMap` vector onto a [`Field`] over `SimplexMesh` (`target_order`).
    pub fn project_scalar_field(
        &self,
        mixed: &MixedDiscretisation,
        dof_map: &MixedLagrangeDofMap,
        name: impl Into<String>,
        unit: impl Into<String>,
        u_mixed: &[f64],
        target_order: usize,
    ) -> Field {
        Field::scalar(name, unit, self.project_lagrange_blocks(mixed, dof_map, u_mixed, 1, target_order))
    }

    /// Projects `ncomp` stacked scalar component blocks onto an interleaved vector [`Field`] over `SimplexMesh`.
    #[allow(clippy::too_many_arguments)]
    pub fn project_vector_field(
        &self,
        mixed: &MixedDiscretisation,
        dof_map: &MixedLagrangeDofMap,
        name: impl Into<String>,
        unit: impl Into<String>,
        comp_names: &[&str],
        u_blocks: &[f64],
        target_order: usize,
    ) -> Field {
        let ncomp = comp_names.len();
        let proj = self.project_lagrange_blocks(mixed, dof_map, u_blocks, ncomp, target_order);
        let n_out = proj.len() / ncomp.max(1);
        let comps: Vec<Vec<f64>> = (0..ncomp).map(|c| proj[c * n_out..(c + 1) * n_out].to_vec()).collect();
        Field::vector(name, unit, comp_names, &comps)
    }

    /// Projects discontinuous per-cell nodal values (`mixed_vals[mixed_cell_offsets[mc] + a]`, evaluated at
    /// `lagrange_ref_nodes(ct, space)`) onto flat `SimplexMesh` per-cell values (`n_simplex_cells * nloc(target_order)`).
    pub fn project_cell_nodal_values(
        &self,
        mixed: &MixedDiscretisation,
        space: LagrangeSpace,
        mixed_vals: &[f64],
        mixed_cell_offsets: &[usize],
        target_order: usize,
    ) -> Vec<f64> {
        let d = self.dim;
        let nloc_out = if target_order == 1 { d + 1 } else { (d + 1) * (d + 2) / 2 };
        let mut out = Vec::with_capacity(self.n_simplex_cells * nloc_out);
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        for sc in 0..self.n_simplex_cells {
            let mc = self.parent_cells[sc];
            let ct = mixed.cell_type(mc);
            let slice = &mixed_vals[mixed_cell_offsets[mc]..mixed_cell_offsets[mc + 1]];
            for xi in self.sub_simplex_ref_points(mixed, sc, target_order) {
                lagrange_shape_into(ct, space, &xi, &mut vals, &mut dxi);
                let mut sum = 0.0;
                for (a, &v) in slice.iter().enumerate() {
                    sum += vals[a] * v;
                }
                out.push(sum);
            }
        }
        out
    }

    /// Projects real or complex per-cell nodal arrays from `MixedDiscretisation` into a `SimplexMesh`-indexed [`CellField`].
    #[allow(clippy::too_many_arguments)]
    pub fn project_cell_field(
        &self,
        mixed: &MixedDiscretisation,
        space: LagrangeSpace,
        name: impl Into<String>,
        unit: impl Into<String>,
        re_mixed: &[f64],
        im_mixed: &[f64],
        mixed_cell_offsets: &[usize],
        target_order: usize,
    ) -> CellField {
        let d = self.dim;
        let nloc = if target_order == 1 { d + 1 } else { (d + 1) * (d + 2) / 2 };
        let re = self.project_cell_nodal_values(mixed, space, re_mixed, mixed_cell_offsets, target_order);
        let im = if im_mixed.is_empty() {
            Vec::new()
        } else {
            self.project_cell_nodal_values(mixed, space, im_mixed, mixed_cell_offsets, target_order)
        };
        CellField { name: name.into(), unit: unit.into(), nloc, re, im }
    }

    /// Evaluates `eval_fn(parent_mc, &xi_ref)` directly at every local node of every sub-simplex `sc`
    /// in `SimplexMesh` (`n_simplex_cells * nloc(target_order)` evaluations in total), ideal for Stage-3
    /// Nedelec `H(curl)` field extraction onto `SimplexMesh` `CellFields` and `node_average` plot fields.
    pub fn sample_simplex_cell_nodes<T, F>(
        &self,
        mixed: &MixedDiscretisation,
        target_order: usize,
        mut eval_fn: F,
    ) -> Vec<T>
    where
        F: FnMut(usize, &[f64; 3]) -> T,
    {
        let d = self.dim;
        let nloc = if target_order == 1 { d + 1 } else { (d + 1) * (d + 2) / 2 };
        let mut out = Vec::with_capacity(self.n_simplex_cells * nloc);
        for sc in 0..self.n_simplex_cells {
            let mc = self.parent_cells[sc];
            for xi in self.sub_simplex_ref_points(mixed, sc, target_order) {
                out.push(eval_fn(mc, &xi));
            }
        }
        out
    }
}

#[cfg(test)]
mod projector_tests {
    use super::*;
    use femcore::adapt::work::simplex_to_adapt;
    use femcore::fields::{node_average, Display};
    use femcore::meshing::pml_box::{build_pml, PmlSpec};

    #[test]
    fn mixed_to_simplex_projector_reproduces_quadratic_fields_and_display_on_3d_tet_prism_hex() {
        let block = femcore::meshing::geomdrive::block(1.0, 0.8, 0.6, 0.45, 3).unwrap().mesh;
        let phys = simplex_to_adapt::<3>(&block).unwrap();
        let cm = build_pml(&phys, &PmlSpec { outer_label: Some(1), ..PmlSpec::<3>::uniform(0.3, 2) }).unwrap().mesh;
        let mesh = crate::Mesh::D3(cm);
        assert!(!mesh.is_pure_simplex());

        let mixed = crate::discretise_mixed(&mesh).unwrap();
        let (simplex, proj) = MixedToSimplexProjector::from_mesh(&mesh, &mixed).unwrap();
        assert!(!proj.is_pure_simplex());
        assert_eq!(proj.n_simplex_cells(), simplex.n_cells());

        // Exact quadratic polynomial p(x, y, z)
        let poly = |x: f64, y: f64, z: f64| -> f64 {
            0.7 + 1.4 * x - 0.9 * y + 0.5 * z + 1.1 * x * x - 0.6 * x * y + 0.8 * z * z
        };

        for sp in [LagrangeSpace::Order2Serendipity, LagrangeSpace::Order2Tensor] {
            let dm = MixedLagrangeDofMap::new(&mixed, sp);
            let mcoords = dm.dof_coordinates(&mixed);
            let u_mixed: Vec<f64> =
                (0..dm.n_dofs).map(|d| poly(mcoords[3 * d], mcoords[3 * d + 1], mcoords[3 * d + 2])).collect();

            let f_proj = proj.project_scalar_field(&mixed, &dm, "u", "V", &u_mixed, 2);
            assert_eq!(f_proj.values.len(), simplex.n_dofs(2));

            let scoords = simplex.dof_coordinates(2);
            for sd in 0..simplex.n_dofs(2) {
                let want = poly(scoords[3 * sd], scoords[3 * sd + 1], scoords[3 * sd + 2]);
                let err = (f_proj.values[sd] - want).abs();
                assert!(err < 1e-13, "{sp:?} sdof {sd}: got {} vs want {want}, err={err:e}", f_proj.values[sd]);
            }

            // Verify Display::new(&simplex, 2) picks finite values cleanly without out-of-bounds!
            let disp = Display::new(&simplex, 2);
            let picked = disp.pick(&f_proj.values);
            assert_eq!(picked.len(), disp.dofs.len());
            assert!(picked.iter().all(|v| v.is_finite()));

            // Also test CellField projection + node_average consistency!
            let mut m_offsets = vec![0usize];
            let mut m_cell_vals = Vec::new();
            for mc in 0..mixed.n_cells() {
                for &gd in dm.cell_dofs(mc) {
                    m_cell_vals.push(u_mixed[gd]);
                }
                m_offsets.push(m_cell_vals.len());
            }
            let cf = proj.project_cell_field(&mixed, sp, "u_cell", "V", &m_cell_vals, &[], &m_offsets, 2);
            assert_eq!(cf.nloc, 10);
            assert_eq!(cf.re.len(), simplex.n_cells() * 10);
            let avg = node_average(&simplex, &cf.re, 1, 2);
            for (&a, &f) in avg.iter().zip(&f_proj.values) {
                assert!((a - f).abs() < 1e-13);
            }
        }
    }
}
