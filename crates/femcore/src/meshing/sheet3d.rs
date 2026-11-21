// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Conformal 3D sheet-stack mesher for layer-stack models containing 2D zero-thickness sheets (`Op::Sheet`).
//!
//! When a 3D model contains embedded 2D conductive sheets (superconducting transmon pads, CPW resonators,
//! Josephson junction straps, ground planes with flux-trap holes), an unstructured BCC lattice mesher cannot
//! conform to 100 nm planar features inside a 10 mm die without billions of cells. Instead, this mesher
//! exploits the exact prismatic layer-stack structure produced by `geom::stack::build`:
//!
//! 1. **Conformal 2D PSLG Triangulation**: `stack.overlay` is the common 2D planar arrangement of all 3D layers
//!    and all 2D sheets (with shared contact corners pinned by `geom::poly_batch`). Every 1D edge of `stack.overlay`
//!    is discretized once with graded chord spacing honoring global `h`, region sizes, face/sheet `Size` rules, and
//!    1D `"<sheet>:edge"` perimeter `Size` rules, then triangulated via `spade` CDT + Ruppert refinement
//!    (`mesher::mesh_chains`). Adjacent sheets and gaps share exact boundary vertices with 0 T-junctions.
//! 2. **Graded Z-Slab Subdivision**: Each 3D layer `[z0, z1]` is subdivided into horizontal z-slabs that match
//!    the local sheet size `dz` at sheet planes `z_s` and grow smoothly by `1.35x` per step away from `z_s` up to
//!    the layer's target element size `h_layer`.
//! 3. **Conformal Dompierre Prism Splitting**: Active 2D triangles in each z-slab are extruded into triangular
//!    prisms and split into 3 tetrahedra using deterministic lowest-global-vertex-id (Dompierre) diagonals so
//!    all vertical and horizontal faces conform identically across neighboring elements.
//! 4. **Single-Facet Emission**: Every outer boundary, material interface (`"rb|ra"` on uncovered cells), and
//!    2D sheet cap triangle (`"<sheet>"` on covered cells) is emitted **exactly once** into `SimplexMesh::facets`
//!    with its 1-based `LabelMap` boundary id.

use crate::adapt::drive::AdaptReport;
use crate::mesh::SimplexMesh;
use crate::meshing::geomdrive::MeshedGeom;
use crate::meshing::mesher::{mesh_chains, Chain, Fill, MeshParams};
use crate::meshing::physics_size::PhysicsSizing;
use crate::meshing::{wave_sizes_by_region, MeshPlan, Mesher, PmlMode, Sizing};
use geom::complex::{EdgeGeom, FaceGeom};
use geom::eval::Evaluated;
use geom::math::{cross3, dist2, dot3, sub3, Aabb, V2, V3};
use geom::model::RegionKind;
use geom::query::DiscretizeOpts;
use geom::select::EntityKind;
use std::collections::HashMap;

fn det3(a: V3, b: V3, c: V3) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn chord_size(r: f64, max_dev: f64) -> f64 {
    2.0 * (2.0 * r * max_dev).sqrt()
}

fn box_dist2(b: &Aabb, p: V2) -> f64 {
    let dx = (b.lo[0] - p[0]).max(p[0] - b.hi[0]).max(0.0);
    let dy = (b.lo[1] - p[1]).max(p[1] - b.hi[1]).max(0.0);
    dx.hypot(dy)
}

/// Subdivides `[z0, z1]` into graded z-levels starting at step `dz0` near `z0` and `dz1` near `z1`,
/// growing by `ratio` (typically `1.35`) per slab up to `h_max`.
pub fn subdivide_z_interval(z0: f64, z1: f64, dz0: f64, dz1: f64, h_max: f64, ratio: f64) -> Vec<f64> {
    let span = z1 - z0;
    if span <= 0.0 {
        return vec![z0];
    }
    let h_max = h_max.min(span).max(1e-15 * span);
    let mut s0 = dz0.min(h_max).max(1e-4 * span);
    let mut s1 = dz1.min(h_max).max(1e-4 * span);
    if span <= 1.25 * s0.min(s1) {
        return vec![z0, z1];
    }

    let mut lo = vec![0.0];
    let mut hi = vec![span];
    while hi.last().unwrap() - lo.last().unwrap() > s0 + s1 {
        if s0 <= s1 {
            let next = lo.last().unwrap() + s0;
            lo.push(next);
            s0 = (s0 * ratio).min(h_max);
        } else {
            let next = hi.last().unwrap() - s1;
            hi.push(next);
            s1 = (s1 * ratio).min(h_max);
        }
    }

    let lo_end = *lo.last().unwrap();
    let hi_end = *hi.last().unwrap();
    let rem = hi_end - lo_end;
    let step_mid = s0.min(s1).min(h_max);
    let n_mid = (rem / step_mid).round().max(1.0) as usize;
    for k in 1..n_mid {
        lo.push(lo_end + rem * (k as f64 / n_mid as f64));
    }
    lo.extend(hi.into_iter().rev());

    let last_idx = lo.len() - 1;
    lo.into_iter()
        .enumerate()
        .map(|(i, d)| {
            if i == 0 {
                z0
            } else if i == last_idx {
                z1
            } else {
                z0 + d
            }
        })
        .collect()
}

/// Meshes a 3D layer-stack model containing 2D sheets (`Op::Sheet`) into a conformal tetrahedral `MeshedGeom`.
pub fn mesh_sheet_stack_3d(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    if ev.list.dim != 3 {
        return Err("mesh_sheet_stack_3d requires a 3D model".into());
    }
    if !(plan.opts.h > 0.0 && plan.opts.h.is_finite()) {
        return Err("the element size must be positive".into());
    }
    let scale_factor = match plan.mesher {
        Mesher::Simple => 1.0,
        Mesher::Adaptive(a) => {
            if !(a.scale > 0.0 && a.scale.is_finite()) {
                return Err("the adaptation scale must be positive".into());
            }
            a.scale
        }
    };
    let h_global = plan.opts.h * scale_factor;
    let max_dev = plan.opts.max_dev.unwrap_or(h_global / 50.0);
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let cx = ev.complex();
    let st = cx.stack.as_ref().ok_or("3D sheet mesher requires a prismatic layer stack")?;
    let ov = &st.overlay;
    let tol = cx.tol;

    let nr = ev.regions().len();
    let mut region_h = vec![h_global; nr];
    for p in ev.pml_regions() {
        region_h[p.region as usize] = region_h[p.region as usize].min((p.thickness / p.layers) * scale_factor);
    }
    if let Some((sizing, sz)) = waves {
        for (r, w) in wave_sizes_by_region(&ev, &labels, sizing, sz)?.into_iter().enumerate().take(nr) {
            if w.is_finite() && w > 0.0 {
                region_h[r] = region_h[r].min(w * scale_factor);
            }
        }
    }

    let nf = cx.faces.len();
    let mut face_h = vec![f64::INFINITY; nf];
    let ne3 = cx.edges3.len();
    let mut edge3_h = vec![f64::INFINITY; ne3];

    for rule in ev.size_rules() {
        match rule.entities.kind {
            EntityKind::Regions => {
                if let Some(rh) = rule.h {
                    for &r in &rule.entities.ids {
                        region_h[r as usize] = region_h[r as usize].min(rh * scale_factor);
                    }
                }
            }
            EntityKind::Faces => {
                for &fi in &rule.entities.ids {
                    let fi = fi as usize;
                    if let Some(rh) = rule.h {
                        face_h[fi] = face_h[fi].min(rh * scale_factor);
                    }
                    if let (Some(md), Some(r)) = (rule.max_dev, cx.facet_radius(fi as u32)) {
                        face_h[fi] = face_h[fi].min(chord_size(r, md) * scale_factor);
                    }
                }
            }
            EntityKind::Edges => {
                if let Some(rh) = rule.h {
                    for &ei in &rule.entities.ids {
                        let ei = ei as usize;
                        if ei < ne3 {
                            edge3_h[ei] = edge3_h[ei].min(rh * scale_factor);
                        }
                    }
                }
            }
            _ => {}
        }
    }

    let mut skip_region = vec![false; nr];
    if matches!(plan.pml, PmlMode::Structured) {
        for p in ev.pml_regions() {
            skip_region[p.region as usize] = true;
        }
    }
    let meshable: Vec<bool> = (0..nr)
        .map(|r| !skip_region[r] && !matches!(labels.domains[r].kind, RegionKind::Hole | RegionKind::Void))
        .collect();

    // Map 3D region/face/edge3 sizes onto the 2D overlay cells and edges
    let ncell = ov.regions.len();
    let nedge = ov.edges.len();
    let mut cell_h = vec![h_global; ncell];
    for (c, ch) in cell_h.iter_mut().enumerate() {
        for l in &st.layers {
            if let Some(r) = l.region_of[c] {
                if meshable[r as usize] {
                    *ch = ch.min(region_h[r as usize]);
                }
            }
        }
    }
    for (fi, f) in cx.faces.iter().enumerate() {
        if let FaceGeom::Cap { cell, .. } = f.geom {
            cell_h[cell as usize] = cell_h[cell as usize].min(face_h[fi]);
        }
    }

    let mut ov_edge_h = vec![h_global; nedge];
    for (e, ed) in ov.edges.iter().enumerate() {
        let mut eh = [ed.left, ed.right].into_iter().flatten().map(|c| cell_h[c as usize]).fold(h_global, f64::min);
        for (fi, f) in cx.faces.iter().enumerate() {
            if let FaceGeom::Lateral { edge, .. } = f.geom {
                if edge == e as u32 {
                    eh = eh.min(face_h[fi]);
                }
            }
        }
        for (ei, e3) in cx.edges3.iter().enumerate() {
            if let EdgeGeom::Horizontal { edge, .. } = e3.geom {
                if edge == e as u32 {
                    eh = eh.min(edge3_h[ei]);
                }
            }
        }
        ov_edge_h[e] = eh;
    }

    let fine_edges: Vec<(usize, f64)> =
        (0..nedge).filter(|&e| ov_edge_h[e] < h_global * 0.999).map(|e| (e, ov_edge_h[e])).collect();
    let fine_cells: Vec<(usize, f64)> =
        (0..ncell).filter(|&c| cell_h[c] < h_global * 0.999).map(|c| (c, cell_h[c])).collect();
    let grade = 0.35;
    let local_h_2d = |p: V2| -> f64 {
        let mut h = h_global;
        for &(e, eh) in &fine_edges {
            let db = box_dist2(&ov.edges[e].bbox, p);
            if eh + grade * db < h {
                let d = dist2(ov.project_edge(e as u32, p), p);
                h = h.min(eh + grade * d);
            }
        }
        for &(c, ch) in &fine_cells {
            let db = box_dist2(&ov.regions[c].bbox, p);
            if ch + grade * db < h {
                h = h.min(ch + grade * db);
            }
        }
        h
    };

    // Discretize every 1D PSLG edge of stack.overlay once (graded along the curve)
    let mut chains = Vec::with_capacity(nedge);
    for (e, ed) in ov.edges.iter().enumerate() {
        let eh = ov_edge_h[e];
        let mut t_pts = ov.discretize_edge(e as u32, DiscretizeOpts { max_dev, max_len: eh });
        let curve = &ov.curves[ed.curve as usize];
        let mut k = 0;
        while k + 1 < t_pts.len() && t_pts.len() < 20_000 {
            let (ta, pa) = (t_pts[k].0, [t_pts[k].1[0], t_pts[k].1[1]]);
            let (tb, pb) = (t_pts[k + 1].0, [t_pts[k + 1].1[0], t_pts[k + 1].1[1]]);
            let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
            let want = local_h_2d(mid).min(eh);
            if dist2(pa, pb) > 1.25 * want {
                let tm = 0.5 * (ta + tb);
                let pm = curve.eval(tm);
                t_pts.insert(k + 1, (tm, [pm[0], pm[1], 0.0]));
            } else {
                k += 1;
            }
        }
        let mut pts: Vec<V2> = t_pts.into_iter().map(|(_, p)| [p[0], p[1]]).collect();
        let last = pts.len() - 1;
        let v0p = ov.vertices[ed.v0 as usize].p;
        let v1p = ov.vertices[ed.v1 as usize].p;
        pts[0] = [v0p[0], v0p[1]];
        pts[last] = [v1p[0], v1p[1]];
        chains.push(Chain { pts, label: (e + 1) as i32, part: e });
    }

    // Interior lattice fills for fine cells and graded collars around fine 1D edges
    let mut fills: Vec<Fill> = Vec::new();
    for &(c, ch) in &fine_cells {
        let bb = ov.regions[c].bbox;
        fills.push(Fill {
            lo: [bb.lo[0], bb.lo[1]],
            hi: [bb.hi[0], bb.hi[1]],
            h: ch,
            owner: Box::new(move |p| ov.locate(p) == Some(c as u32)),
        });
    }
    for &(e, eh) in &fine_edges {
        for (mult, reach_k) in [(1.2, 2.2), (2.4, 5.0)] {
            let h_band = eh * mult;
            if h_band >= h_global * 0.999 {
                continue;
            }
            let reach = eh * reach_k;
            let bb = ov.edges[e].bbox.pad(reach);
            let cell_h_ref = &cell_h;
            fills.push(Fill {
                lo: [bb.lo[0], bb.lo[1]],
                hi: [bb.hi[0], bb.hi[1]],
                h: h_band,
                owner: Box::new(move |p| {
                    ov.locate(p).is_some_and(|c| cell_h_ref[c as usize] > h_band)
                        && dist2(ov.project_edge(e as u32, p), p) <= reach
                }),
            });
        }
    }

    let mut ov_bb = Aabb::EMPTY;
    for r in &ov.regions {
        ov_bb = ov_bb.union(&r.bbox);
    }
    let h_fine_min =
        cell_h.iter().copied().chain(ov_edge_h.iter().copied()).fold(h_global, f64::min).max(1e-12 * h_global);
    let cell_active: Vec<bool> =
        (0..ncell).map(|c| st.layers.iter().any(|l| l.region_of[c].is_some_and(|r| meshable[r as usize]))).collect();
    let domain_at = |cen: V2| -> Option<i32> {
        let c = ov.locate(cen)?;
        if cell_active[c as usize] {
            Some(c as i32 + 1)
        } else {
            None
        }
    };
    let m2 = mesh_chains(
        &chains,
        &fills,
        &MeshParams {
            h: h_global,
            h_fine_min,
            h_reach: h_global,
            min_angle_deg: plan.opts.min_angle_deg,
            lo: [ov_bb.lo[0], ov_bb.lo[1]],
            hi: [ov_bb.hi[0], ov_bb.hi[1]],
        },
        &domain_at,
    )?;
    drop(fills);

    // Map every 2D mesh edge separating distinct overlay cells (or the outer boundary) to its overlay edge id
    let ntri2 = m2.n_cells();
    let tri_ov_cell: Vec<usize> = m2.cell_labels.iter().map(|&l| (l - 1) as usize).collect();
    let mut edge2_sides: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for (ti, &c) in tri_ov_cell.iter().enumerate() {
        let v = m2.cell(ti);
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            edge2_sides.entry((a.min(b), a.max(b))).or_default().push(c);
        }
    }
    let pt2 = |v: u32| -> V2 { [m2.points[2 * v as usize], m2.points[2 * v as usize + 1]] };
    let mut bnd_edges2: Vec<(u32, u32, u32)> = Vec::new();
    let mut edge_list: Vec<((u32, u32), Vec<usize>)> = edge2_sides.into_iter().collect();
    edge_list.sort_by_key(|x| x.0);
    for ((a, b), sides) in edge_list {
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

    // Build graded z-slabs across every 3D layer
    let sheet_dz_at = |z: f64, default_h: f64| -> Option<f64> {
        let mut has_sheet = false;
        let mut dz = default_h;
        for (fi, f) in cx.faces.iter().enumerate() {
            if let FaceGeom::Cap { cell, z: fz } = f.geom {
                if (fz - z).abs() <= tol && f.is_sheet() {
                    has_sheet = true;
                    dz = dz.min(face_h[fi]).min(cell_h[cell as usize]);
                }
            }
        }
        for (ei, e3) in cx.edges3.iter().enumerate() {
            if let EdgeGeom::Horizontal { edge, z: ez } = e3.geom {
                if (ez - z).abs() <= tol && e3.faces.iter().any(|&fi| cx.faces[fi as usize].is_sheet()) {
                    has_sheet = true;
                    dz = dz.min(edge3_h[ei]).min(ov_edge_h[edge as usize]);
                }
            }
        }
        has_sheet.then_some(dz)
    };

    let mut z_planes: Vec<f64> = Vec::new();
    // For each slab interval [z_planes[s], z_planes[s+1]], the 3D layer index `Some(k)` (or `None` across a void gap)
    let mut slab_layer: Vec<Option<usize>> = Vec::new();
    for (k, l) in st.layers.iter().enumerate() {
        let active_in_layer = l.region_of.iter().flatten().any(|&r| meshable[r as usize]);
        if !active_in_layer {
            continue;
        }
        let h_layer = l
            .region_of
            .iter()
            .flatten()
            .filter(|&&r| meshable[r as usize])
            .map(|&r| region_h[r as usize])
            .fold(h_global, f64::min);
        let dz0 = sheet_dz_at(l.z0, h_layer).unwrap_or(h_layer);
        let dz1 = sheet_dz_at(l.z1, h_layer).unwrap_or(h_layer);
        let zs = subdivide_z_interval(l.z0, l.z1, dz0, dz1, h_layer, 1.35);
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
    let n2 = m2.n_points();

    // Mark which 3D grid vertices (v2, zp) belong to at least one active prism
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

    // Extrude active 2D triangles across each z-slab into 3 Dompierre tetrahedra
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

    // Build lookup tables for horizontal cap faces (including 2D sheets) and vertical lateral faces
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

    // 1. Horizontal cap facets (outer caps, uncovered "rb|ra" interfaces, and 2D sheet caps):
    //    emitted ONCE per triangle with its 1-based LabelMap boundary id.
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

    // 2. Vertical lateral facets: split along the exact Dompierre diagonal (vid(a, s), vid(b, s + 1))
    //    where a < b, matching the adjacent prism tetrahedra on both sides.
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
    let mut g = MeshedGeom { mesh, labels, eval: std::sync::Arc::new(ev) };
    if matches!(plan.pml, PmlMode::Structured) && !g.eval.pml_regions().is_empty() {
        g = crate::meshing::pml_geom3::structured_pml_3d(g)?.0;
    }
    Ok((g, None))
}
