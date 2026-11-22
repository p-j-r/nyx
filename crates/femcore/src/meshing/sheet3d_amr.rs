// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! 2D-overlay (`OverlayOracle2`) and full-3D sheet-preserving adaptive mesh refinement helpers
//! for prismatic `sheet3d` models (`AF-1`).
//!
//! Provides:
//! - [`OverlayOracle2`]: exact 2D geometry oracle over `LayerStack::overlay` in the internal `(x, y)`
//!   sweep plane with pinned `overlay.vertices`, used to adapt the 2D `m2` triangulation without losing
//!   exact lines, circles, or ellipses (`st.overlay.edges`).
//! - [`extract_overlay_sheet_features`]: extracts 2D sheet perimeter segments and sharp `90/270 deg`
//!   corners from `Complex` for [`crate::adapt::electro::electro_overlay_metrics_2d`].
//! - [`adapt_overlay_2d`], [`build_overlay_mesh_2d`], and [`extrude_overlay_2d`]: building blocks for
//!   column-projected `sheet_2d_z` multi-pass AMR.
//! - [`adapt_sheet_full_3d`] and [`sheet_rim_vertices_3d`]: unstructured 3D `adapt3` preserving
//!   two-sided interior `Op::Sheet` patches (`Sig::Sheet`) and `Complex::edges3` rim creases.

use crate::adapt::drive::{adapt2_with, AdaptOptions, AdaptReport};
use crate::adapt::field::{BackgroundMetric2, MetricField};
use crate::adapt::geom::{Geometry2, PolylineGeometry};
use crate::adapt::metric::Sym;
use crate::adapt::topo::{Sig, Work2};
use crate::adapt::work::{simplex_from_adapt, simplex_to_adapt};
use crate::mesh::SimplexMesh;
use crate::meshing::adapt3d::adapt_to_metrics3;
use crate::meshing::geomdrive::MeshedGeom;
use crate::meshing::physics_size::PhysicsSizing;
pub use crate::meshing::sheet3d::{build_sized_overlay_mesh_2d, extrude_overlay_2d, extrude_overlay_2d_with_plan};
use crate::meshing::{Adaptive, MeshPlan, Sizing};
use geom::complex::{Complex, EdgeGeom};
use geom::curve::Curve;
use geom::eval::Evaluated;
use geom::math::{dist2, V2};
use std::collections::BTreeMap;

/// Exact 2D geometry oracle over `LayerStack::overlay` in the `(x, y)` sweep plane.
///
/// Maps each constrained chain of [`Work2`] (`m2` cell labels = 1-based overlay region indices)
/// to its exact `overlay.edges` curves so split midpoints and smoothed chain vertices land on
/// exact lines, circles, and ellipses while `overlay.vertices` remain pinned.
pub struct OverlayOracle2<'a> {
    pub ov: &'a Complex,
    cands: Vec<Vec<u32>>,
    base: PolylineGeometry,
}

impl<'a> OverlayOracle2<'a> {
    /// Builds the overlay oracle from a classified 2D working mesh `w` over `ov`.
    pub fn from_work(ov: &'a Complex, w: &Work2) -> Self {
        let bb = crate::adapt::work::bbox_of(&w.pts);
        let diag = (bb[0][1] - bb[0][0]).hypot(bb[1][1] - bb[1][0]).max(1e-30);
        let end_tol = 1e-5 * diag;

        let snap = w.snapshot();
        let mut chain_segs: Vec<Vec<(V2, V2)>> = vec![Vec::new(); w.chains.len()];
        let mut chain_owners: Vec<Vec<u32>> = vec![Vec::new(); w.chains.len()];
        let mut es: Vec<_> = w.chain_of.iter().collect();
        es.sort_unstable();
        for (e, &cid) in es {
            chain_segs[cid as usize].push((w.pts[e[0] as usize], w.pts[e[1] as usize]));
            if let Some(idx) = snap.index(e) {
                for &ti in &snap.edges[idx].1 {
                    if ti >= 0 {
                        let tl = w.tlabel[ti as usize];
                        if tl >= 1 {
                            let r = (tl - 1) as u32;
                            if !chain_owners[cid as usize].contains(&r) {
                                chain_owners[cid as usize].push(r);
                            }
                        }
                    }
                }
            }
        }

        let seg_matches_edge = |segs: &[(V2, V2)], ei: u32| -> bool {
            segs.iter().any(|&(pa, pb)| {
                dist2(ov.project_edge(ei, pa), pa) <= end_tol && dist2(ov.project_edge(ei, pb), pb) <= end_tol
            })
        };

        let cands = (0..w.chains.len())
            .map(|ci| {
                let segs = &chain_segs[ci];
                match w.chains[ci].sig {
                    Sig::Iface(a, b) if a >= 1 && b >= 1 => {
                        let (ra, rb) = ((a - 1) as u32, (b - 1) as u32);
                        let mut v: Vec<u32> = (0..ov.edges.len() as u32)
                            .filter(|&ei| {
                                let ed = &ov.edges[ei as usize];
                                (ed.left == Some(ra) && ed.right == Some(rb))
                                    || (ed.left == Some(rb) && ed.right == Some(ra))
                            })
                            .collect();
                        let matched: Vec<u32> = v.iter().copied().filter(|&ei| seg_matches_edge(segs, ei)).collect();
                        if !matched.is_empty() {
                            v = matched;
                        }
                        v
                    }
                    Sig::Bnd(l) | Sig::Sheet(l) => {
                        let mut out = Vec::new();
                        if l >= 1 && (l as usize) <= ov.edges.len() {
                            let ei = (l - 1) as u32;
                            if seg_matches_edge(segs, ei) {
                                out.push(ei);
                            }
                        }
                        if out.is_empty() {
                            let owners = &chain_owners[ci];
                            for (ei, ed) in ov.edges.iter().enumerate() {
                                let touches_owner = owners.is_empty()
                                    || owners.iter().any(|&r| ed.left == Some(r) || ed.right == Some(r));
                                if touches_owner && seg_matches_edge(segs, ei as u32) {
                                    out.push(ei as u32);
                                }
                            }
                        }
                        out
                    }
                    _ => Vec::new(),
                }
            })
            .collect();

        OverlayOracle2 { ov, cands, base: PolylineGeometry::from_work(w) }
    }

    /// Mesh vertex indices of `m2` that coincide with `ov.vertices` in `(x, y)` (pinned as `Corner`).
    pub fn pinned_vertices(ov: &Complex, m2: &SimplexMesh) -> Vec<u32> {
        if m2.dim != 2 || m2.n_points() == 0 {
            return Vec::new();
        }
        let (lo, hi) = m2.bbox();
        let scale = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
        let tol = 1e-6 * scale;
        let mut out = Vec::new();
        for gv in &ov.vertices {
            let q = [gv.p[0], gv.p[1]];
            let best = (0..m2.n_points())
                .map(|i| (i, (m2.points[2 * i] - q[0]).hypot(m2.points[2 * i + 1] - q[1])))
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            if let Some((i, d)) = best {
                if d <= tol {
                    out.push(i as u32);
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }
}

impl Geometry2 for OverlayOracle2<'_> {
    fn project(&self, p: [f64; 2], chain: u32) -> [f64; 2] {
        let cands = &self.cands[chain as usize];
        if cands.is_empty() || std::env::var_os("NYX_NO_ORACLE").is_some() {
            return self.base.project(p, chain);
        }
        let mut best = (f64::INFINITY, p);
        for &e in cands {
            let bb = &self.ov.edges[e as usize].bbox;
            let dx = (bb.lo[0] - p[0]).max(p[0] - bb.hi[0]).max(0.0);
            let dy = (bb.lo[1] - p[1]).max(p[1] - bb.hi[1]).max(0.0);
            if dx.hypot(dy) >= best.0 {
                continue;
            }
            let q = self.ov.project_edge(e, p);
            let d = dist2(q, p);
            if d < best.0 {
                best = (d, q);
            }
        }
        best.1
    }
}

/// Extracted 2D sheet perimeter segments, sharp corners, and sweep-plane world axes from a 3D stack `Complex`.
#[derive(Clone, Debug)]
pub struct OverlaySheetFeatures {
    /// Discretized 2D segments `(a, b)` along all `sheet:edge` perimeters in `st.overlay` coordinates.
    pub rim_segments: Vec<([f64; 2], [f64; 2])>,
    /// Sharp 2D corner coordinates (`90/270 deg` turns and curve junctions) on `sheet:edge` perimeters.
    pub corners: Vec<[f64; 2]>,
    /// World coordinate indices `[ax, ay]` corresponding to `st.overlay`'s `(x, y)` plane.
    pub plane_axes: [usize; 2],
}

/// Extracts 2D sheet rim segments and sharp corners from a 3D prismatic stack `Complex`.
pub fn extract_overlay_sheet_features(cx: &Complex) -> Option<OverlaySheetFeatures> {
    let st = cx.stack.as_ref()?;
    let ov = &st.overlay;
    let perm = geom::scene::axis_perm(cx.axis);
    let plane_axes = [perm[0], perm[1]];

    let mut rim_ov_edges: Vec<u32> = Vec::new();
    for sname in cx.sheet_names() {
        for ei in cx.sheet_perimeter_edges(&sname) {
            if let EdgeGeom::Horizontal { edge, .. } = cx.edges3[ei as usize].geom {
                rim_ov_edges.push(edge);
            }
        }
    }
    rim_ov_edges.sort_unstable();
    rim_ov_edges.dedup();

    let mut rim_segments = Vec::new();
    // Outgoing unit tangents at each overlay vertex along incident rim edges
    let mut vert_dirs: BTreeMap<u32, Vec<[f64; 2]>> = BTreeMap::new();
    for &oe in &rim_ov_edges {
        let ed = &ov.edges[oe as usize];
        let curve = &ov.curves[ed.curve as usize];
        let nsub = if matches!(curve, Curve::Line { .. }) { 1 } else { 16 };
        for k in 0..nsub {
            let ta = ed.t0 + (ed.t1 - ed.t0) * (k as f64 / nsub as f64);
            let tb = ed.t0 + (ed.t1 - ed.t0) * ((k + 1) as f64 / nsub as f64);
            let pa = if k == 0 {
                let p = ov.vertices[ed.v0 as usize].p;
                [p[0], p[1]]
            } else {
                curve.eval(ta)
            };
            let pb = if k + 1 == nsub {
                let p = ov.vertices[ed.v1 as usize].p;
                [p[0], p[1]]
            } else {
                curve.eval(tb)
            };
            rim_segments.push((pa, pb));
        }
        // Outgoing tangent from v0 (at t0, +tangent) and outgoing tangent from v1 (at t1, -tangent)
        let tg0 = curve.tangent(ed.t0);
        let tg1 = curve.tangent(ed.t1);
        if ed.v0 != ed.v1 {
            vert_dirs.entry(ed.v0).or_default().push(tg0);
            vert_dirs.entry(ed.v1).or_default().push([-tg1[0], -tg1[1]]);
        }
    }

    let cos_smooth = (160.0_f64.to_radians()).cos(); // collinear outgoing tangents point 180 deg apart (cos = -1)
    let mut corners = Vec::new();
    for (vid, dirs) in vert_dirs {
        let is_corner = if dirs.len() != 2 {
            true
        } else {
            let d = dirs[0][0] * dirs[1][0] + dirs[0][1] * dirs[1][1];
            // If smooth continuation, outgoing directions are opposite (`d ~ -1.0 <= cos(160 deg)`)
            d > cos_smooth
        };
        if is_corner {
            let p = ov.vertices[vid as usize].p;
            corners.push([p[0], p[1]]);
        }
    }

    Some(OverlaySheetFeatures { rim_segments, corners, plane_axes })
}

/// Adapts a 2D overlay triangulation `m2` over `ov` (`st.overlay`) to vertex metrics `metrics`,
/// grading `metrics` with `a.grade`, pinning `ov.vertices`, and projecting constrained chains onto
/// exact `ov.edges` via [`OverlayOracle2`].
pub fn adapt_overlay_2d(
    ov: &Complex,
    m2: &SimplexMesh,
    metrics: &[Sym<2>],
    a: &Adaptive,
) -> Result<(SimplexMesh, AdaptReport), String> {
    if m2.dim != 2 {
        return Err("adapt_overlay_2d requires a 2D overlay mesh".into());
    }
    if metrics.len() != m2.n_points() {
        return Err(format!("{} metrics for {} overlay vertices", metrics.len(), m2.n_points()));
    }
    let old = simplex_to_adapt::<2>(m2)?;
    let field = BackgroundMetric2::new_graded(&old, metrics, a.grade);
    adapt_overlay_2d_with_adapt_mesh(ov, m2, &old, &field, a)
}

/// Adapts a 2D overlay triangulation `m2` over `ov` (`st.overlay`) to a general [`MetricField<2>`].
pub fn adapt_overlay_2d_field(
    ov: &Complex,
    m2: &SimplexMesh,
    field: &dyn MetricField<2>,
    a: &Adaptive,
) -> Result<(SimplexMesh, AdaptReport), String> {
    if m2.dim != 2 {
        return Err("adapt_overlay_2d_field requires a 2D overlay mesh".into());
    }
    let old = simplex_to_adapt::<2>(m2)?;
    adapt_overlay_2d_with_adapt_mesh(ov, m2, &old, field, a)
}

fn adapt_overlay_2d_with_adapt_mesh(
    ov: &Complex,
    m2: &SimplexMesh,
    old: &crate::adapt::work::SimplexMesh<2>,
    field: &dyn MetricField<2>,
    a: &Adaptive,
) -> Result<(SimplexMesh, AdaptReport), String> {
    let mut opts = AdaptOptions { max_cells: a.max_cells, ..AdaptOptions::default() };
    opts.topo.corner_angle_deg = 75.0;
    opts.topo.fixed_vertices = OverlayOracle2::pinned_vertices(ov, m2);
    let result = adapt2_with(old, field, &opts, |w| OverlayOracle2::from_work(ov, w))?;
    Ok((simplex_from_adapt(&result.mesh), result.report))
}

/// Builds the initial 2D conformal triangulation `m2` over `ev.complex().stack.overlay` (`cell_labels`
/// = `ov_cell + 1`, boundary `facet_labels` = `ov_edge + 1`) using [`build_sized_overlay_mesh_2d`].
pub fn build_overlay_mesh_2d(ev: &Evaluated, plan: &MeshPlan) -> Result<SimplexMesh, String> {
    build_sized_overlay_mesh_2d(ev, plan, None)
}

/// Identifies vertices of a 3D `SimplexMesh` lying on open 1D perimeter rims of `sheet_labels`
/// (edges incident to exactly one sheet facet of a given sheet label).
pub fn sheet_rim_vertices_3d(mesh: &SimplexMesh, sheet_labels: &[i32]) -> Vec<bool> {
    let nv = mesh.n_points();
    if mesh.dim != 3 || nv == 0 || sheet_labels.is_empty() {
        return vec![false; nv];
    }
    let mut per_label_edges: BTreeMap<(i32, u32, u32), usize> = BTreeMap::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if sheet_labels.contains(&lab) {
            let fv = mesh.facet(f);
            for (i, j) in [(0, 1), (1, 2), (2, 0)] {
                let (a, b) = (fv[i].min(fv[j]), fv[i].max(fv[j]));
                *per_label_edges.entry((lab, a, b)).or_insert(0) += 1;
            }
        }
    }
    let mut is_rim = vec![false; nv];
    for (&(_, a, b), &cnt) in &per_label_edges {
        if cnt == 1 {
            is_rim[a as usize] = true;
            is_rim[b as usize] = true;
        }
    }
    is_rim
}

/// Adapts a 3D [`MeshedGeom`] containing interior 2D `Op::Sheet` boundaries (`l0 == l1` or `l0 != l1`)
/// to per-vertex metrics `metrics`, preserving `Sig::Sheet` patches, projecting `z = const` rim creases
/// onto exact `Complex::edges3`, and emitting each sheet facet once without duplicating interface facets.
pub fn adapt_sheet_full_3d(g: MeshedGeom, metrics: &[Sym<3>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String> {
    adapt_to_metrics3(g, metrics, a)
}

/// Same as [`adapt_sheet_full_3d`], and re-applies `AF-5` (`curved`) + `AF-2` (`quarter_point`)
/// `EdgeNodes` enrichment according to `plan`.
pub fn adapt_sheet_full_3d_with_plan(
    g: MeshedGeom,
    metrics: &[Sym<3>],
    a: &Adaptive,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<(MeshedGeom, AdaptReport), String> {
    let (mut g_out, report) = adapt_to_metrics3(g, metrics, a)?;
    crate::meshing::sheet3d::apply_sheet_post_enrichment(&mut g_out, plan, waves);
    Ok((g_out, report))
}
