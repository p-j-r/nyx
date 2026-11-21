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
use crate::meshing::mesher::{mesh_chains, Chain, Fill, MeshParams};
use crate::meshing::sheet3d::subdivide_z_interval;
use crate::meshing::{Adaptive, MeshPlan, PmlMode};
use geom::complex::{Complex, EdgeGeom, FaceGeom};
use geom::curve::Curve;
use geom::eval::Evaluated;
use geom::math::{cross3, dist2, dot3, sub3, Aabb, V2, V3};
use geom::model::RegionKind;
use geom::query::DiscretizeOpts;
use std::collections::{BTreeMap, HashMap};

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
/// = `ov_cell + 1`, boundary `facet_labels` = `ov_edge + 1`).
pub fn build_overlay_mesh_2d(ev: &Evaluated, plan: &MeshPlan) -> Result<SimplexMesh, String> {
    if ev.list.dim != 3 {
        return Err("build_overlay_mesh_2d requires a 3D stack model".into());
    }
    if !(plan.opts.h > 0.0 && plan.opts.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    let scale = match plan.mesher {
        crate::meshing::Mesher::Simple => 1.0,
        crate::meshing::Mesher::Adaptive(a) => a.scale,
    };
    let h_global = plan.opts.h * scale;
    let max_dev = plan.opts.max_dev.unwrap_or(h_global / 50.0);
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let cx = ev.complex();
    let st = cx.stack.as_ref().ok_or("model has no prismatic layer stack")?;
    let ov = &st.overlay;

    let nr = ev.regions().len();
    let mut skip_region = vec![false; nr];
    if matches!(plan.pml, PmlMode::Structured) {
        for p in ev.pml_regions() {
            skip_region[p.region as usize] = true;
        }
    }
    let meshable: Vec<bool> = (0..nr)
        .map(|r| !skip_region[r] && !matches!(labels.domains[r].kind, RegionKind::Hole | RegionKind::Void))
        .collect();

    let ncell = ov.regions.len();
    let nedge = ov.edges.len();
    let mut chains = Vec::with_capacity(nedge);
    for (e, ed) in ov.edges.iter().enumerate() {
        let t_pts = ov.discretize_edge(e as u32, DiscretizeOpts { max_dev, max_len: h_global });
        let mut pts: Vec<V2> = t_pts.into_iter().map(|(_, p)| [p[0], p[1]]).collect();
        let last = pts.len() - 1;
        let v0p = ov.vertices[ed.v0 as usize].p;
        let v1p = ov.vertices[ed.v1 as usize].p;
        pts[0] = [v0p[0], v0p[1]];
        pts[last] = [v1p[0], v1p[1]];
        chains.push(Chain { pts, label: (e + 1) as i32, part: e });
    }

    let mut ov_bb = Aabb::EMPTY;
    for r in &ov.regions {
        ov_bb = ov_bb.union(&r.bbox);
    }
    let cell_active: Vec<bool> =
        (0..ncell).map(|c| st.layers.iter().any(|l| l.region_of[c].is_some_and(|r| meshable[r as usize]))).collect();
    let domain_at = |cen: V2| -> Option<i32> {
        let c = ov.locate(cen)?;
        cell_active[c as usize].then_some(c as i32 + 1)
    };
    let fills: Vec<Fill> = Vec::new();
    mesh_chains(
        &chains,
        &fills,
        &MeshParams {
            h: h_global,
            h_fine_min: h_global * 1e-6,
            h_reach: h_global,
            min_angle_deg: plan.opts.min_angle_deg,
            lo: [ov_bb.lo[0], ov_bb.lo[1]],
            hi: [ov_bb.hi[0], ov_bb.hi[1]],
        },
        &domain_at,
    )
}

/// Extrudes a conformal 2D overlay mesh `m2` across the 3D layers of `ev` into a conformal 3D
/// tetrahedral [`MeshedGeom`] with graded vertical steps matched to the local sheet resolution in `m2`.
pub fn extrude_overlay_2d(
    ev: Evaluated,
    m2: &SimplexMesh,
    h_global: f64,
    sheet_dz_override: Option<f64>,
) -> Result<MeshedGeom, String> {
    if ev.list.dim != 3 || m2.dim != 2 {
        return Err("extrude_overlay_2d requires a 3D Evaluated model and a 2D overlay mesh".into());
    }
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let cx = ev.complex();
    let st = cx.stack.as_ref().ok_or("extrude_overlay_2d requires a prismatic layer stack")?;
    let ov = &st.overlay;
    let tol = cx.tol;

    let nr = ev.regions().len();
    let meshable: Vec<bool> =
        (0..nr).map(|r| !matches!(labels.domains[r].kind, RegionKind::Hole | RegionKind::Void)).collect();

    let ntri2 = m2.n_cells();
    let n2 = m2.n_points();
    let tri_ov_cell: Vec<usize> = m2.cell_labels.iter().map(|&l| (l - 1) as usize).collect();
    let pt2 = |v: u32| -> V2 { [m2.points[2 * v as usize], m2.points[2 * v as usize + 1]] };

    // Map 2D mesh edges separating distinct overlay cells (or outer boundary) to overlay edge ids
    let mut edge2_sides: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for (ti, &c) in tri_ov_cell.iter().enumerate() {
        let v = m2.cell(ti);
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            edge2_sides.entry((a.min(b), a.max(b))).or_default().push(c);
        }
    }
    let mut bnd_edges2: Vec<(u32, u32, u32)> = Vec::new();
    let mut edge_list: Vec<((u32, u32), Vec<usize>)> = edge2_sides.into_iter().collect();
    edge_list.sort_by_key(|x| x.0);
    for &((a, b), ref sides) in &edge_list {
        let (c0, c1) = match sides.as_slice() {
            [c0] => (Some(*c0 as u32), None),
            [c0, c1] if c0 != c1 => (Some(*c0 as u32), Some(*c1 as u32)),
            _ => continue,
        };
        let (pa, pb) = (pt2(a), pt2(b));
        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
        let mut best: Option<(f64, u32)> = None;
        for pass in 0..2 {
            for (ei, ed) in ov.edges.iter().enumerate() {
                if pass == 0 && !((ed.left == c0 && ed.right == c1) || (ed.left == c1 && ed.right == c0)) {
                    continue;
                }
                let d = dist2(ov.project_edge(ei as u32, mid), mid);
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, ei as u32));
                }
            }
            if best.is_some() {
                break;
            }
        }
        if let Some((_, oe)) = best {
            bnd_edges2.push((a, b, oe));
        }
    }

    // Infer characteristic sheet vertical step `dz` at each sheet plane `z` from `m2` triangle edges on the sheet
    let sheet_dz_at = |z: f64, default_h: f64| -> Option<f64> {
        let mut sheet_cells = vec![false; ov.regions.len()];
        let mut any = false;
        for f in &cx.faces {
            if f.is_sheet() {
                if let FaceGeom::Cap { cell, z: fz } = f.geom {
                    if (fz - z).abs() <= tol {
                        sheet_cells[cell as usize] = true;
                        any = true;
                    }
                }
            }
        }
        if !any {
            return None;
        }
        if let Some(ovr) = sheet_dz_override {
            return Some(ovr.min(default_h));
        }
        let mut min_edge = default_h;
        for (ti, &ov_c) in tri_ov_cell.iter().enumerate() {
            if ov_c < sheet_cells.len() && sheet_cells[ov_c] {
                let v = m2.cell(ti);
                for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                    min_edge = min_edge.min(dist2(pt2(a), pt2(b)));
                }
            }
        }
        // Prevent extreme pancake slab ratios while resolving fine rim scales
        Some(min_edge.max(default_h * 1e-3).min(default_h))
    };

    let mut z_planes: Vec<f64> = Vec::new();
    let mut slab_layer: Vec<Option<usize>> = Vec::new();
    for (k, l) in st.layers.iter().enumerate() {
        let active_in_layer = l.region_of.iter().flatten().any(|&r| meshable[r as usize]);
        if !active_in_layer {
            continue;
        }
        let dz0 = sheet_dz_at(l.z0, h_global).unwrap_or(h_global);
        let dz1 = sheet_dz_at(l.z1, h_global).unwrap_or(h_global);
        let zs = subdivide_z_interval(l.z0, l.z1, dz0, dz1, h_global, 1.35);
        if let Some(&last_z) = z_planes.last() {
            if (l.z0 - last_z).abs() > tol {
                z_planes.push(l.z0);
                slab_layer.push(None);
            }
        } else {
            z_planes.push(zs[0]);
        }
        for &zk in &zs[1..] {
            z_planes.push(zk);
            slab_layer.push(Some(k));
        }
    }
    if z_planes.len() < 2 {
        return Err("the 3D model has no active layers to mesh".into());
    }
    let nz = z_planes.len();
    let nslab = nz - 1;

    let nf = cx.faces.len();
    let face_label: Vec<i32> = (0..nf as u32)
        .map(|fi| labels.label_of_edge(fi).ok_or_else(|| format!("face '{}' has no boundary label", cx.facet_name(fi))))
        .collect::<Result<Vec<i32>, String>>()?;
    let face_normal: Vec<V3> = (0..nf as u32)
        .map(|fi| cx.facet_normals(fi).into_iter().next().unwrap_or([0.0, 0.0, 1.0]))
        .collect();

    let mut cap_at: HashMap<(u32, usize), usize> = HashMap::new();
    let mut lat_at: HashMap<(u32, usize), usize> = HashMap::new();
    for (fi, f) in cx.faces.iter().enumerate() {
        let eff = |r: Option<u32>| r.filter(|&r| meshable[r as usize]);
        if eff(f.first).is_none() && eff(f.second).is_none() {
            continue;
        }
        match f.geom {
            FaceGeom::Cap { cell, z } => {
                if let Some(zp) = z_planes.iter().position(|&zp| (zp - z).abs() <= tol) {
                    cap_at.insert((cell, zp), fi);
                }
            }
            FaceGeom::Lateral { edge, z0, z1 } => {
                for (k, l) in st.layers.iter().enumerate() {
                    if l.z0 >= z0 - tol && l.z1 <= z1 + tol {
                        lat_at.insert((edge, k), fi);
                    }
                }
            }
            FaceGeom::Quadric(_) | FaceGeom::Planar { .. } => {}
        }
    }

    let mut used_v = vec![u32::MAX; nz * n2];
    for (s, &layer_opt) in slab_layer.iter().enumerate() {
        let Some(k) = layer_opt else { continue };
        let l = &st.layers[k];
        for (ti, &ov_c) in tri_ov_cell.iter().enumerate().take(ntri2) {
            let Some(r) = l.region_of[ov_c] else { continue };
            if !meshable[r as usize] {
                continue;
            }
            for &v in m2.cell(ti) {
                used_v[s * n2 + v as usize] = 0;
                used_v[(s + 1) * n2 + v as usize] = 0;
            }
        }
    }

    let mut points: Vec<f64> = Vec::new();
    let mut n_pts3 = 0u32;
    for zp in 0..nz {
        let z = z_planes[zp];
        for v2 in 0..n2 {
            let slot = &mut used_v[zp * n2 + v2];
            if *slot == 0 {
                *slot = n_pts3;
                n_pts3 += 1;
                let p2 = pt2(v2 as u32);
                let pw = cx.wv([p2[0], p2[1], z]);
                points.extend_from_slice(&pw);
            }
        }
    }
    let vid = |v2: u32, zp: usize| -> u32 { used_v[zp * n2 + v2 as usize] };
    let pt3 = |v3: u32| -> V3 {
        let i = 3 * v3 as usize;
        [points[i], points[i + 1], points[i + 2]]
    };
    let det3 = |a: V3, b: V3, c: V3| -> f64 {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
    };

    let mut cells: Vec<u32> = Vec::new();
    let mut cell_labels: Vec<i32> = Vec::new();
    for (s, &layer_opt) in slab_layer.iter().enumerate() {
        let Some(k) = layer_opt else { continue };
        let l = &st.layers[k];
        for (ti, &ov_c) in tri_ov_cell.iter().enumerate().take(ntri2) {
            let Some(r) = l.region_of[ov_c] else { continue };
            if !meshable[r as usize] {
                continue;
            }
            let dom_id = labels.domains[r as usize].id;
            let mut v = [m2.cell(ti)[0], m2.cell(ti)[1], m2.cell(ti)[2]];
            v.sort_unstable();
            let (b0, b1, b2) = (vid(v[0], s), vid(v[1], s), vid(v[2], s));
            let (t0, t1, t2) = (vid(v[0], s + 1), vid(v[1], s + 1), vid(v[2], s + 1));
            for mut tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                let (p0, p1, p2, p3) = (pt3(tet[0]), pt3(tet[1]), pt3(tet[2]), pt3(tet[3]));
                if det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0)) < 0.0 {
                    tet.swap(1, 2);
                }
                cells.extend_from_slice(&tet);
                cell_labels.push(dom_id);
            }
        }
    }

    let mut facets: Vec<u32> = Vec::new();
    let mut facet_labels: Vec<i32> = Vec::new();
    let mut push_oriented_tri = |mut tri: [u32; 3], fi: usize| {
        let (p0, p1, p2) = (pt3(tri[0]), pt3(tri[1]), pt3(tri[2]));
        let n = cross3(sub3(p1, p0), sub3(p2, p0));
        if dot3(n, face_normal[fi]) < 0.0 {
            tri.swap(1, 2);
        }
        facets.extend_from_slice(&tri);
        facet_labels.push(face_label[fi]);
    };

    for zp in 0..nz {
        for (ti, &ov_c) in tri_ov_cell.iter().enumerate().take(ntri2) {
            if let Some(&fi) = cap_at.get(&(ov_c as u32, zp)) {
                let v = m2.cell(ti);
                let tri = [vid(v[0], zp), vid(v[1], zp), vid(v[2], zp)];
                if tri.iter().all(|&w| w != u32::MAX) {
                    push_oriented_tri(tri, fi);
                }
            }
        }
    }

    for &(a, b, oe) in &bnd_edges2 {
        let (a, b) = (a.min(b), a.max(b));
        for (s, &layer_opt) in slab_layer.iter().enumerate().take(nslab) {
            let Some(k) = layer_opt else { continue };
            if let Some(&fi) = lat_at.get(&(oe, k)) {
                let (ba, bb, ta, tb) = (vid(a, s), vid(b, s), vid(a, s + 1), vid(b, s + 1));
                if [ba, bb, ta, tb].iter().all(|&w| w != u32::MAX) {
                    push_oriented_tri([ba, bb, tb], fi);
                    push_oriented_tri([ba, tb, ta], fi);
                }
            }
        }
    }

    let mesh = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);
    Ok(MeshedGeom { mesh, labels, eval: std::sync::Arc::new(ev) })
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
