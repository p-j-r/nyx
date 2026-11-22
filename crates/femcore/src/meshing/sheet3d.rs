// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Conformal 3D sheet-stack mesher for layer-stack models containing 2D zero-thickness sheets (`Op::Sheet`).
//!
//! When a 3D model contains embedded 2D conductive sheets (superconducting transmon pads, CPW resonators,
//! Josephson junction straps, ground planes with flux-trap holes), an unstructured BCC lattice mesher cannot
//! conform to 100 nm planar features inside a 10 mm die without billions of cells. Instead, this mesher
//! exploits the exact prismatic layer-stack structure produced by `geom::stack::build`:
//!
//! 1. **Conformal 2D PSLG Triangulation & Coplanar Gap Proximity Sizing**: `stack.overlay` is the common 2D planar
//!    arrangement of all 3D layers and all 2D sheets (with shared contact corners pinned by `geom::poly_batch`).
//!    Optional `gap_cells` auto-detects narrow coplanar gaps between facing sheet perimeters via an AABB-pruned
//!    segment-distance pass and refines both rims symmetrically. Every 1D edge of `stack.overlay` is discretized
//!    once with graded chord spacing, then triangulated via `spade` CDT + Ruppert refinement (`mesher::mesh_chains`).
//!    Adjacent sheets and gaps share exact boundary vertices with 0 T-junctions.
//! 2. **Graded Z-Slab Subdivision**: Each 3D layer `[z0, z1]` is subdivided into horizontal z-slabs that match
//!    the local sheet size `dz` at sheet planes `z_s` and grow smoothly by `1.35x` per step away from `z_s` up to
//!    the layer's target element size `h_layer`.
//! 3. **Conformal Dompierre Prism Splitting & Vertical Coarsening (`z_coarsen`)**: Active 2D triangles in each z-slab
//!    are extruded into triangular prisms and split into 3 tetrahedra using deterministic Dompierre diagonals. When
//!    `z_coarsen` is enabled, slabs within `|z - z_sheet| <= 4 * h_sheet` (at least 2 slabs per side) keep the fine
//!    2D triangulation `M_0`, while outer slabs apply maximal independent-set 2D half-edge collapses bridged by 100%
//!    conforming pure-tetrahedral 3-tet / 2-tet transition slabs.
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
use geom::curve::Curve;
use geom::eval::Evaluated;
use geom::math::{cross3, dist2, dot3, sub3, Aabb, V2, V3};
use geom::model::RegionKind;
use geom::query::DiscretizeOpts;
use geom::select::EntityKind;
use std::collections::{HashMap, HashSet};

fn det3(a: V3, b: V3, c: V3) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn cross2_pts(o: V2, a: V2, b: V2) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

fn tri_quality_ccw(a: V2, b: V2, c: V2) -> f64 {
    let area2 = cross2_pts(a, b, c);
    if area2 <= 0.0 {
        return 0.0;
    }
    let l2 = (b[0] - a[0]).powi(2)
        + (b[1] - a[1]).powi(2)
        + (c[0] - b[0]).powi(2)
        + (c[1] - b[1]).powi(2)
        + (a[0] - c[0]).powi(2)
        + (a[1] - c[1]).powi(2);
    if l2 <= 0.0 {
        0.0
    } else {
        3.0f64.sqrt() * area2 / l2
    }
}

fn chord_size(r: f64, max_dev: f64) -> f64 {
    2.0 * (2.0 * r * max_dev).sqrt()
}

fn box_dist2(b: &Aabb, p: V2) -> f64 {
    let dx = (b.lo[0] - p[0]).max(p[0] - b.hi[0]).max(0.0);
    let dy = (b.lo[1] - p[1]).max(p[1] - b.hi[1]).max(0.0);
    dx.hypot(dy)
}

fn aabb_dist2(a: &Aabb, b: &Aabb) -> f64 {
    let dx = (a.lo[0] - b.hi[0]).max(b.lo[0] - a.hi[0]).max(0.0);
    let dy = (a.lo[1] - b.hi[1]).max(b.lo[1] - a.hi[1]).max(0.0);
    dx.hypot(dy)
}

fn closest_point_on_seg(p: V2, a: V2, b: V2) -> V2 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let l2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if l2 > 0.0 {
        (((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    [a[0] + t * ab[0], a[1] + t * ab[1]]
}

fn seg_seg_closest(a0: V2, a1: V2, b0: V2, b1: V2) -> (f64, V2, V2) {
    let c1 = cross2_pts(a0, a1, b0);
    let c2 = cross2_pts(a0, a1, b1);
    let c3 = cross2_pts(b0, b1, a0);
    let c4 = cross2_pts(b0, b1, a1);
    if c1 * c2 < 0.0 && c3 * c4 < 0.0 {
        let denom = c1 - c2;
        let s = if denom.abs() > 0.0 { (c1 / denom).clamp(0.0, 1.0) } else { 0.5 };
        let p = [b0[0] + s * (b1[0] - b0[0]), b0[1] + s * (b1[1] - b0[1])];
        return (0.0, p, p);
    }
    let candidates = [
        (a0, closest_point_on_seg(a0, b0, b1)),
        (a1, closest_point_on_seg(a1, b0, b1)),
        (closest_point_on_seg(b0, a0, a1), b0),
        (closest_point_on_seg(b1, a0, a1), b1),
    ];
    let mut best = (f64::INFINITY, a0, b0);
    for (pa, pb) in candidates {
        let d = dist2(pa, pb);
        if d < best.0 {
            best = (d, pa, pb);
        }
    }
    best
}

#[derive(Clone, Copy)]
struct BvhItem2 {
    id: usize,
    bbox: Aabb,
    min_h: f64,
}

enum BvhKind2 {
    Leaf { start: usize, end: usize },
    Branch { left: usize, right: usize },
}

struct BvhNode2 {
    bbox: Aabb,
    min_h: f64,
    kind: BvhKind2,
}

/// Flat `Vec`-backed 2D AABB BVH with subtree `min_h` augmentation for $O(\log N)$ rim pairing,
/// graded sizing (`local_h_2d`), and nearest-overlay-edge queries (`bnd_edges2`).
struct Bvh2 {
    items: Vec<BvhItem2>,
    nodes: Vec<BvhNode2>,
}

impl Bvh2 {
    fn build(mut items: Vec<BvhItem2>) -> Self {
        let mut nodes = Vec::with_capacity(items.len().saturating_mul(2));
        if !items.is_empty() {
            let n = items.len();
            Self::build_rec(&mut items, 0, n, &mut nodes);
        }
        Self { items, nodes }
    }

    fn build_rec(items: &mut [BvhItem2], start: usize, end: usize, nodes: &mut Vec<BvhNode2>) -> usize {
        let mut bbox = Aabb::EMPTY;
        let mut min_h = f64::INFINITY;
        let mut c_lo = [f64::INFINITY; 2];
        let mut c_hi = [f64::NEG_INFINITY; 2];
        for it in &items[start..end] {
            bbox = bbox.union(&it.bbox);
            min_h = min_h.min(it.min_h);
            let cx = 0.5 * (it.bbox.lo[0] + it.bbox.hi[0]);
            let cy = 0.5 * (it.bbox.lo[1] + it.bbox.hi[1]);
            c_lo[0] = c_lo[0].min(cx);
            c_lo[1] = c_lo[1].min(cy);
            c_hi[0] = c_hi[0].max(cx);
            c_hi[1] = c_hi[1].max(cy);
        }
        let idx = nodes.len();
        nodes.push(BvhNode2 { bbox, min_h, kind: BvhKind2::Leaf { start, end } });
        if end - start <= 8 {
            return idx;
        }
        let axis = usize::from(c_hi[1] - c_lo[1] > c_hi[0] - c_lo[0]);
        let mid = start + (end - start) / 2;
        items[start..end].select_nth_unstable_by(mid - start, |a, b| {
            let ca = 0.5 * (a.bbox.lo[axis] + a.bbox.hi[axis]);
            let cb = 0.5 * (b.bbox.lo[axis] + b.bbox.hi[axis]);
            ca.total_cmp(&cb).then_with(|| a.id.cmp(&b.id))
        });
        let left = Self::build_rec(items, start, mid, nodes);
        let right = Self::build_rec(items, mid, end, nodes);
        nodes[idx].kind = BvhKind2::Branch { left, right };
        idx
    }

    fn for_each_within_box(&self, query: &Aabb, margin: f64, mut visit: impl FnMut(usize)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = [0usize; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp]];
            if aabb_dist2(&node.bbox, query) >= margin {
                continue;
            }
            match node.kind {
                BvhKind2::Leaf { start, end } => {
                    for it in &self.items[start..end] {
                        if aabb_dist2(&it.bbox, query) < margin {
                            visit(it.id);
                        }
                    }
                }
                BvhKind2::Branch { left, right } => {
                    stack[sp] = right;
                    stack[sp + 1] = left;
                    sp += 2;
                }
            }
        }
    }

    fn min_graded(&self, p: V2, grade: f64, mut best_h: f64, mut eval_leaf: impl FnMut(usize, f64, f64) -> f64) -> f64 {
        if self.nodes.is_empty() {
            return best_h;
        }
        let mut stack = [0usize; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp]];
            if node.min_h + grade * box_dist2(&node.bbox, p) >= best_h {
                continue;
            }
            match node.kind {
                BvhKind2::Leaf { start, end } => {
                    for it in &self.items[start..end] {
                        let db = box_dist2(&it.bbox, p);
                        if it.min_h + grade * db < best_h {
                            best_h = best_h.min(eval_leaf(it.id, it.min_h, db));
                        }
                    }
                }
                BvhKind2::Branch { left, right } => {
                    let lb_l = self.nodes[left].min_h + grade * box_dist2(&self.nodes[left].bbox, p);
                    let lb_r = self.nodes[right].min_h + grade * box_dist2(&self.nodes[right].bbox, p);
                    if lb_l <= lb_r {
                        if lb_r < best_h {
                            stack[sp] = right;
                            sp += 1;
                        }
                        if lb_l < best_h {
                            stack[sp] = left;
                            sp += 1;
                        }
                    } else {
                        if lb_l < best_h {
                            stack[sp] = left;
                            sp += 1;
                        }
                        if lb_r < best_h {
                            stack[sp] = right;
                            sp += 1;
                        }
                    }
                }
            }
        }
        best_h
    }

    fn nearest_id(&self, p: V2, mut eval_dist: impl FnMut(usize) -> Option<f64>) -> Option<(f64, u32)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<(f64, u32)> = None;
        let mut stack = [0usize; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp]];
            if best.is_some_and(|(bd, _)| box_dist2(&node.bbox, p) > bd) {
                continue;
            }
            match node.kind {
                BvhKind2::Leaf { start, end } => {
                    for it in &self.items[start..end] {
                        if best.is_some_and(|(bd, _)| box_dist2(&it.bbox, p) > bd) {
                            continue;
                        }
                        if let Some(d) = eval_dist(it.id) {
                            let id32 = it.id as u32;
                            if best.is_none_or(|(bd, bid)| d < bd || (d == bd && id32 < bid)) {
                                best = Some((d, id32));
                            }
                        }
                    }
                }
                BvhKind2::Branch { left, right } => {
                    let dl = box_dist2(&self.nodes[left].bbox, p);
                    let dr = box_dist2(&self.nodes[right].bbox, p);
                    let keep = |d: f64, best: Option<(f64, u32)>| best.is_none_or(|(bd, _)| d <= bd);
                    if dl <= dr {
                        if keep(dr, best) {
                            stack[sp] = right;
                            sp += 1;
                        }
                        if keep(dl, best) {
                            stack[sp] = left;
                            sp += 1;
                        }
                    } else {
                        if keep(dl, best) {
                            stack[sp] = left;
                            sp += 1;
                        }
                        if keep(dr, best) {
                            stack[sp] = right;
                            sp += 1;
                        }
                    }
                }
            }
        }
        best
    }
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
        .map(|(i, d)| if i == 0 { z0 } else if i == last_idx { z1 } else { z0 + d })
        .collect()
}

/// Subdivides `[z0, z1]` for `z_coarsen = true`: uses `ratio` (default `1.35x`) geometric progression inside
/// `|z - z_sheet| <= 16 * h_sheet` (preserving near- and mid-field vertical resolution identically to
/// `subdivide_z_interval`), then accelerates to `ratio.max(1.52)` in the far dielectric/vacuum bulk up to `h_max`.
fn subdivide_z_interval_far(z0: f64, z1: f64, dz0: f64, dz1: f64, h_max: f64, ratio: f64) -> Vec<f64> {
    let span = z1 - z0;
    if span <= 0.0 {
        return vec![z0];
    }
    let h_max = h_max.min(span).max(1e-15 * span);
    let d0_init = dz0.min(h_max).max(1e-4 * span);
    let d1_init = dz1.min(h_max).max(1e-4 * span);
    let mut s0 = d0_init;
    let mut s1 = d1_init;
    if span <= 1.25 * s0.min(s1) {
        return vec![z0, z1];
    }

    let r_far = 1.52_f64.max(ratio);
    let mut lo = vec![0.0];
    let mut hi = vec![span];
    while hi.last().unwrap() - lo.last().unwrap() > s0 + s1 {
        if s0 <= s1 {
            let cur = *lo.last().unwrap();
            lo.push(cur + s0);
            let r = if cur + s0 > 16.0 * d0_init { r_far } else { ratio };
            s0 = (s0 * r).min(h_max);
        } else {
            let cur = *hi.last().unwrap();
            hi.push(cur - s1);
            let r = if span - (cur - s1) > 16.0 * d1_init { r_far } else { ratio };
            s1 = (s1 * r).min(h_max);
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
        .map(|(i, d)| if i == 0 { z0 } else if i == last_idx { z1 } else { z0 + d })
        .collect()
}

struct Level2D {
    /// CCW-oriented 2D triangles `(v0, v1, v2)` using original `m2` vertex indices, paired with overlay cell index.
    tris: Vec<([u32; 3], usize)>,
}

struct SlabTransition {
    map_v: Vec<u32>,
    rank: Vec<(u8, usize, u8)>,
    /// For triangles in `Level2D::tris` whose edge `(u, v)` collapsed `v -> u`, stores `Some((u, v, c))`.
    tri_collapse: Vec<Option<(u32, u32, u32)>>,
}

/// Applies one maximal independent-set 2D half-edge collapse pass `M_m -> M_{m+1}` for edges up to `h_target`,
/// guaranteeing:
/// 1. All frozen boundary/interface vertices remain untouched;
/// 2. Collapsed edges are vertex-disjoint and every collapsed-edge triangle `(u, v, c)` has `c` uncollapsed;
/// 3. The 2D Link Condition holds at every collapse so `M_{m+1}` is strictly manifold without duplicate edges;
/// 4. Every surviving 2D triangle stays strictly CCW with quality `>= 0.25`, and all 3 Dompierre tetrahedra on the
///    unit transition slab `[0, 1]` have consistent strictly positive volume.
fn coarsen_2d_once(cur: &Level2D, pts2: &[V2], frozen_v: &[bool], h_target: f64, h_max_edge: f64) -> (Level2D, SlabTransition) {
    let n2 = pts2.len();
    let mut work_tris = cur.tris.clone();
    let ntri = cur.tris.len();

    // Build initial candidate half-edges `(v -> u)` sorted deterministically by length
    let mut seen_undir: HashSet<(u32, u32)> = HashSet::new();
    let mut candidates: Vec<(u64, u32, u32)> = Vec::new();
    for &(v, _) in &cur.tris {
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            let key = (a.min(b), a.max(b));
            if !seen_undir.insert(key) {
                continue;
            }
            let len = dist2(pts2[a as usize], pts2[b as usize]);
            if len > h_target {
                continue;
            }
            let lbits = len.to_bits();
            if !frozen_v[b as usize] {
                candidates.push((lbits, b, a));
            }
            if !frozen_v[a as usize] {
                candidates.push((lbits, a, b));
            }
        }
    }
    candidates.sort_unstable();

    // Map each vertex in `cur.tris` to its incident original triangle indices (for fast 3D unit-slab prism check)
    let mut orig_v_tris: Vec<Vec<usize>> = vec![Vec::new(); n2];
    for (ti, &(v, _)) in cur.tris.iter().enumerate() {
        for &w in &v {
            orig_v_tris[w as usize].push(ti);
        }
    }

    // Track active triangles in `work_tris` and per-vertex active triangle lists
    let mut tri_alive = vec![true; ntri];
    let mut work_v_tris: Vec<Vec<usize>> = orig_v_tris.clone();

    let mut map_v: Vec<u32> = (0..n2 as u32).collect();
    let mut collapse_pair_of: Vec<Option<(usize, u8)>> = vec![None; n2];
    let mut n_collapses = 0usize;

    let rank_of = |w: u32, cp: &[Option<(usize, u8)>]| -> (u8, usize, u8) {
        match cp[w as usize] {
            Some((cid, sub)) => (0, cid, sub),
            None => (1, w as usize, 0),
        }
    };

    for (_, v_rem, u_keep) in candidates {
        if frozen_v[v_rem as usize] {
            continue;
        }
        if collapse_pair_of[v_rem as usize].is_some() || collapse_pair_of[u_keep as usize].is_some() {
            continue;
        }

        // Find triangles in `work_tris` currently incident to `(u_keep, v_rem)` and check Link Condition
        let mut shared_tris: Vec<(usize, u32)> = Vec::new();
        let mut nbrs_u: HashSet<u32> = HashSet::new();
        for &ti in &work_v_tris[u_keep as usize] {
            if !tri_alive[ti] {
                continue;
            }
            let tv = work_tris[ti].0;
            for &w in &tv {
                if w != u_keep {
                    nbrs_u.insert(w);
                }
            }
        }
        if !nbrs_u.contains(&v_rem) {
            continue;
        }
        let mut nbrs_v: HashSet<u32> = HashSet::new();
        for &ti in &work_v_tris[v_rem as usize] {
            if !tri_alive[ti] {
                continue;
            }
            let tv = work_tris[ti].0;
            let mut has_u = false;
            let mut opp = u32::MAX;
            for &w in &tv {
                if w != v_rem {
                    nbrs_v.insert(w);
                    if w == u_keep {
                        has_u = true;
                    } else {
                        opp = w;
                    }
                }
            }
            if has_u && opp != u32::MAX {
                shared_tris.push((ti, opp));
            }
        }
        // Interior manifold edge must have exactly 2 incident triangles in both `work_tris` and `cur.tris`
        if shared_tris.len() != 2 {
            continue;
        }
        let (c1, c2) = (shared_tris[0].1, shared_tris[1].1);
        if c1 == c2 {
            continue;
        }
        let orig_shared = orig_v_tris[v_rem as usize]
            .iter()
            .filter(|&&ti| cur.tris[ti].0.contains(&u_keep))
            .count();
        if orig_shared != 2 {
            continue;
        }

        // Exact 2D Link Condition: N(u) \cap N(v) == {c1, c2}
        let common_cnt = nbrs_v.iter().filter(|w| nbrs_u.contains(w)).count();
        if common_cnt != 2 {
            continue;
        }

        // Check that all non-collapsed triangles in `star_work(v_rem)` remain strictly CCW, well-shaped, and bounded in size
        let pu = pts2[u_keep as usize];
        let max_len2 = (1.6 * h_max_edge) * (1.6 * h_max_edge);
        let mut geom_ok = true;
        for &ti in &work_v_tris[v_rem as usize] {
            if !tri_alive[ti] || ti == shared_tris[0].0 || ti == shared_tris[1].0 {
                continue;
            }
            let mut tv = work_tris[ti].0;
            for w in &mut tv {
                if *w == v_rem {
                    *w = u_keep;
                }
            }
            let (p0, p1, p2) = (pts2[tv[0] as usize], pts2[tv[1] as usize], pts2[tv[2] as usize]);
            if tri_quality_ccw(p0, p1, p2) < 0.12 {
                geom_ok = false;
                break;
            }
            for q in [p0, p1, p2] {
                if (q[0] - pu[0]).powi(2) + (q[1] - pu[1]).powi(2) > max_len2 {
                    geom_ok = false;
                    break;
                }
            }
            if !geom_ok {
                break;
            }
        }
        if !geom_ok {
            continue;
        }

        // Check 3D Dompierre prism/wedge positivity on unit slab `[0, 1]` for all original triangles in `star(u) \cup star(v)`
        let cid = n_collapses;
        collapse_pair_of[u_keep as usize] = Some((cid, 0));
        collapse_pair_of[v_rem as usize] = Some((cid, 1));
        map_v[v_rem as usize] = u_keep;

        let mut prism_ok = true;
        for &ti in orig_v_tris[u_keep as usize].iter().chain(orig_v_tris[v_rem as usize].iter()) {
            let tv = cur.tris[ti].0;
            let mut col_edge = None;
            let mut n_col_edges = 0usize;
            for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
                let (a, b, c) = (tv[i], tv[j], tv[k]);
                if map_v[b as usize] == a && a != b {
                    col_edge = Some((a, b, c));
                    n_col_edges += 1;
                } else if map_v[a as usize] == b && a != b {
                    col_edge = Some((b, a, c));
                    n_col_edges += 1;
                }
            }
            if n_col_edges > 1 {
                prism_ok = false;
                break;
            }
            if let Some((u_col, v_col, c_opp)) = col_edge {
                let rc = rank_of(c_opp, &collapse_pair_of);
                if !(rank_of(u_col, &collapse_pair_of) < rc && rank_of(v_col, &collapse_pair_of) < rc) {
                    prism_ok = false;
                    break;
                }
                let c_top = map_v[c_opp as usize];
                if c_top == u_col || c_top == v_col {
                    prism_ok = false;
                    break;
                }
                let a_old = cross2_pts(pts2[u_col as usize], pts2[v_col as usize], pts2[c_opp as usize]);
                let a_new = cross2_pts(pts2[u_col as usize], pts2[v_col as usize], pts2[c_top as usize]);
                if a_old.signum() * a_new <= 1e-6 * a_old.abs() {
                    prism_ok = false;
                    break;
                }
                continue;
            }
            let mut x = tv;
            x.sort_by_key(|&w| rank_of(w, &collapse_pair_of));
            let b = [
                [pts2[x[0] as usize][0], pts2[x[0] as usize][1], 0.0],
                [pts2[x[1] as usize][0], pts2[x[1] as usize][1], 0.0],
                [pts2[x[2] as usize][0], pts2[x[2] as usize][1], 0.0],
            ];
            let (m0, m1, m2) = (map_v[x[0] as usize], map_v[x[1] as usize], map_v[x[2] as usize]);
            let t = [
                [pts2[m0 as usize][0], pts2[m0 as usize][1], 1.0],
                [pts2[m1 as usize][0], pts2[m1 as usize][1], 1.0],
                [pts2[m2 as usize][0], pts2[m2 as usize][1], 1.0],
            ];
            let area_sign = cross2_pts([b[0][0], b[0][1]], [b[1][0], b[1][1]], [b[2][0], b[2][1]]);
            let area_mag = area_sign.abs();
            let sigma = area_sign.signum();
            for [p0, p1, p2, p3] in [[b[0], b[1], b[2], t[2]], [b[0], b[1], t[2], t[1]], [b[0], t[0], t[1], t[2]]] {
                let d3 = det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0));
                if sigma * d3 <= 1e-6 * area_mag {
                    prism_ok = false;
                    break;
                }
            }
            if !prism_ok {
                break;
            }
        }
        if !prism_ok {
            collapse_pair_of[u_keep as usize] = None;
            collapse_pair_of[v_rem as usize] = None;
            map_v[v_rem as usize] = v_rem;
            continue;
        }

        // Commit collapse `v_rem -> u_keep` in `work_tris`
        n_collapses += 1;
        tri_alive[shared_tris[0].0] = false;
        tri_alive[shared_tris[1].0] = false;
        let moved_tris: Vec<usize> = work_v_tris[v_rem as usize]
            .iter()
            .copied()
            .filter(|&ti| tri_alive[ti])
            .collect();
        for ti in moved_tris {
            for w in &mut work_tris[ti].0 {
                if *w == v_rem {
                    *w = u_keep;
                }
            }
            work_v_tris[u_keep as usize].push(ti);
        }
        work_v_tris[v_rem as usize].clear();
    }

    let rank: Vec<(u8, usize, u8)> = (0..n2 as u32).map(|w| rank_of(w, &collapse_pair_of)).collect();
    let mut tri_collapse = vec![None; ntri];
    for (ti, &(tv, _)) in cur.tris.iter().enumerate() {
        for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
            let (a, b, c) = (tv[i], tv[j], tv[k]);
            if map_v[b as usize] == a && a != b {
                tri_collapse[ti] = Some((a, b, c));
                break;
            } else if map_v[a as usize] == b && a != b {
                tri_collapse[ti] = Some((b, a, c));
                break;
            }
        }
    }
    let next_tris: Vec<([u32; 3], usize)> = work_tris
        .into_iter()
        .enumerate()
        .filter_map(|(ti, t)| tri_alive[ti].then_some(t))
        .collect();

    (Level2D { tris: next_tris }, SlabTransition { map_v, rank, tri_collapse })
}

fn curve_sheet_mesh_edges(
    mesh: &mut SimplexMesh,
    cx: &geom::complex::Complex,
    pt2: impl Fn(u32) -> V2,
    bnd_edges2: &[(u32, u32, u32)],
    z_planes: &[f64],
    used_v: &[u32],
    n2: usize,
) {
    let Some(st) = cx.stack.as_ref() else { return };
    let ov = &st.overlay;
    let mut cand: crate::meshing::curve3::EdgeNodes = Vec::new();
    for &(a, b, oe) in bnd_edges2 {
        let ed = &ov.edges[oe as usize];
        if matches!(ov.curves[ed.curve as usize], Curve::Line { .. }) {
            continue;
        }
        let (pa2, pb2) = (pt2(a), pt2(b));
        let len2 = dist2(pa2, pb2);
        if !(len2 > 0.0) {
            continue;
        }
        let mid2 = [0.5 * (pa2[0] + pb2[0]), 0.5 * (pa2[1] + pb2[1])];
        let q2 = ov.project_edge(oe, mid2);
        let off2 = dist2(q2, mid2);
        if !(off2 > 1e-9 * len2 && off2 < 0.25 * len2 && q2[0].is_finite() && q2[1].is_finite()) {
            continue;
        }
        for (zp, &z) in z_planes.iter().enumerate() {
            let va = used_v[zp * n2 + a as usize];
            let vb = used_v[zp * n2 + b as usize];
            if va != u32::MAX && vb != u32::MAX && va != vb {
                let pr = [va.min(vb), va.max(vb)];
                let qw = cx.wv([q2[0], q2[1], z]);
                cand.push((pr, qw));
            }
        }
    }
    if !cand.is_empty() {
        crate::meshing::curve3::apply_edge_nodes_validated(mesh, cand);
    }
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
    if let Some(gc) = plan.opts.gap_cells {
        if !(gc > 0.0 && gc.is_finite()) {
            return Err("gap_cells must be positive".into());
        }
    }
    if let Some(dz_min) = plan.opts.dz_min {
        if !(dz_min > 0.0 && dz_min.is_finite()) {
            return Err("dz_min must be finite and positive".into());
        }
    }
    if let Some(z_ratio) = plan.opts.z_growth_ratio {
        if !(z_ratio > 1.0 && z_ratio.is_finite()) {
            return Err("z_growth_ratio must be finite and strictly greater than 1.0".into());
        }
    }
    if let Some(max_slabs) = plan.opts.max_z_slabs {
        if max_slabs < 1 {
            return Err("max_z_slabs must be at least 1".into());
        }
    }
    if let Some(mc) = plan.opts.max_cells {
        if mc < 1 {
            return Err("max_cells must be at least 1".into());
        }
    }
    let (scale_factor, adapt_max_cells) = match plan.mesher {
        Mesher::Simple => (1.0, None),
        Mesher::Adaptive(a) => {
            if !(a.scale > 0.0 && a.scale.is_finite()) {
                return Err("the adaptation scale must be positive".into());
            }
            if let Some(mc) = a.max_cells {
                if mc < 1 {
                    return Err("max_cells must be at least 1".into());
                }
            }
            (a.scale, a.max_cells)
        }
    };
    let eff_max_cells = match (plan.opts.max_cells, adapt_max_cells) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
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
        ov_edge_h[e] = [ed.left, ed.right].into_iter().flatten().map(|c| cell_h[c as usize]).fold(h_global, f64::min);
    }
    for (fi, f) in cx.faces.iter().enumerate() {
        if let FaceGeom::Lateral { edge, .. } = f.geom {
            let e = edge as usize;
            if e < nedge {
                ov_edge_h[e] = ov_edge_h[e].min(face_h[fi]);
            }
        }
    }
    for (ei, e3) in cx.edges3.iter().enumerate() {
        if let EdgeGeom::Horizontal { edge, .. } = e3.geom {
            let e = edge as usize;
            if e < nedge {
                ov_edge_h[e] = ov_edge_h[e].min(edge3_h[ei]);
            }
        }
    }

    let grade = 0.35;
    let has_explicit_fine =
        ov_edge_h.iter().any(|&h| h < h_global * 0.999) || cell_h.iter().any(|&h| h < h_global * 0.999);

    // Coplanar gap proximity pass across all 2D sheet planes:
    // 1. When `gap_cells = Some(gc)`, clamps `ov_edge_h[e] <= d_gap / gc` symmetrically on both facing rims.
    // 2. Even when `gap_cells = None`, prevents unsplit cyclic-trapezoid Ruppert encroachment across narrow gaps
    //    facing explicitly sized conductors or curved-hole pockets (QFR-5).
    let mut sheet_zs: Vec<f64> = Vec::new();
    for f in &cx.faces {
        if f.is_sheet() {
            if let FaceGeom::Cap { z, .. } = f.geom {
                if !sheet_zs.iter().any(|&zs| (zs - z).abs() <= tol) {
                    sheet_zs.push(z);
                }
            }
        }
    }
    if plan.opts.gap_cells.is_some() || has_explicit_fine {
        struct RimEdge {
            e: usize,
            v0: u32,
            v1: u32,
            bbox: Aabb,
            segs: Vec<(V2, V2, V2)>,
        }
        for &z_s in &sheet_zs {
            let mut cell_sheet: Vec<Option<&str>> = vec![None; ncell];
            for f in &cx.faces {
                if f.is_sheet() {
                    if let FaceGeom::Cap { cell, z } = f.geom {
                        if (z - z_s).abs() <= tol {
                            cell_sheet[cell as usize] = Some(f.name.as_str());
                        }
                    }
                }
            }
            let mut rims: Vec<RimEdge> = Vec::new();
            for (e, ed) in ov.edges.iter().enumerate() {
                if ed.left.is_none() || ed.right.is_none() {
                    continue;
                }
                let ls = ed.left.and_then(|c| cell_sheet[c as usize]);
                let rs = ed.right.and_then(|c| cell_sheet[c as usize]);
                if ls.is_some() == rs.is_some() {
                    continue;
                }
                let n_sign = if ls.is_some() { 1.0 } else { -1.0 };
                let curve = &ov.curves[ed.curve as usize];
                let nsub = if matches!(curve, Curve::Line { .. }) { 1 } else { 8 };
                let mut segs = Vec::with_capacity(nsub);
                for k in 0..nsub {
                    let ta = ed.t0 + (ed.t1 - ed.t0) * (k as f64 / nsub as f64);
                    let tb = ed.t0 + (ed.t1 - ed.t0) * ((k + 1) as f64 / nsub as f64);
                    let tm = 0.5 * (ta + tb);
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
                    let n = ov.edge_normal(e as u32, tm);
                    segs.push((pa, pb, [n_sign * n[0], n_sign * n[1]]));
                }
                rims.push(RimEdge { e, v0: ed.v0, v1: ed.v1, bbox: ed.bbox, segs });
            }

            let rim_bvh = Bvh2::build(
                rims.iter().enumerate().map(|(i, r)| BvhItem2 { id: i, bbox: r.bbox, min_h: 0.0 }).collect(),
            );
            let mut gap_clamps: Vec<(usize, f64)> = Vec::new();
            for (i, ra) in rims.iter().enumerate() {
                rim_bvh.for_each_within_box(&ra.bbox, h_global, |j| {
                    if j <= i {
                        return;
                    }
                    let rb = &rims[j];
                    if ra.v0 == rb.v0 || ra.v0 == rb.v1 || ra.v1 == rb.v0 || ra.v1 == rb.v1 {
                        return;
                    }
                    let mut best: Option<(f64, V2, V2, V2, V2)> = None;
                    for &(a0, a1, na) in &ra.segs {
                        for &(b0, b1, nb) in &rb.segs {
                            let (d, pa, pb) = seg_seg_closest(a0, a1, b0, b1);
                            if best.is_none_or(|b| d < b.0) {
                                best = Some((d, pa, pb, na, nb));
                            }
                        }
                    }
                    let Some((d, pa, pb, na, nb)) = best else { return };
                    if !(d > 1e-15 * h_global && d < h_global) {
                        return;
                    }
                    let dp = [pb[0] - pa[0], pb[1] - pa[1]];
                    let dot_a = dp[0] * na[0] + dp[1] * na[1];
                    let dot_b = -dp[0] * nb[0] - dp[1] * nb[1];
                    let dot_nn = na[0] * nb[0] + na[1] * nb[1];
                    if !(dot_a > 0.15 * d && dot_b > 0.15 * d && dot_nn < 0.25) {
                        return;
                    }
                    if let Some(gc) = plan.opts.gap_cells {
                        let h_gap = d / gc;
                        gap_clamps.push((ra.e, h_gap));
                        gap_clamps.push((rb.e, h_gap));
                    }
                    let h_pair_min = ov_edge_h[ra.e].min(ov_edge_h[rb.e]);
                    if has_explicit_fine && h_pair_min < h_global * 0.999 && d <= 2.5 * h_pair_min {
                        let h_cap = d.min(h_pair_min + grade * d);
                        gap_clamps.push((ra.e, h_cap));
                        gap_clamps.push((rb.e, h_cap));
                    }
                });
            }
            for (ei, hc) in gap_clamps {
                ov_edge_h[ei] = ov_edge_h[ei].min(hc);
            }
        }
    }

    let mut fine_edges: Vec<(usize, f64)> =
        (0..nedge).filter(|&e| ov_edge_h[e] < h_global * 0.999).map(|e| (e, ov_edge_h[e])).collect();
    fine_edges.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let fine_cells: Vec<(usize, f64)> =
        (0..ncell).filter(|&c| cell_h[c] < h_global * 0.999).map(|c| (c, cell_h[c])).collect();
    let fine_edge_bvh = Bvh2::build(
        fine_edges
            .iter()
            .map(|&(e, eh)| BvhItem2 { id: e, bbox: ov.edges[e].bbox, min_h: eh })
            .collect(),
    );
    let fine_cell_bvh = Bvh2::build(
        fine_cells
            .iter()
            .map(|&(c, ch)| BvhItem2 { id: c, bbox: ov.regions[c].bbox, min_h: ch })
            .collect(),
    );
    let local_h_2d = |p: V2| -> f64 {
        let h = fine_edge_bvh.min_graded(p, grade, h_global, |e, eh, _db| {
            let d = dist2(ov.project_edge(e as u32, p), p);
            eh + grade * d
        });
        fine_cell_bvh.min_graded(p, grade, h, |_c, ch, db| ch + grade * db)
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
            owner: Box::new(move |p| ov.locate(p) == Some(c as u32) && local_h_2d(p) >= 0.4 * ch),
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
                        && local_h_2d(p) >= 0.4 * h_band
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
    let (m2, mc_rep) = mesh_chains(
        &chains,
        &fills,
        &MeshParams {
            h: h_global,
            h_fine_min,
            min_angle_deg: plan.opts.min_angle_deg,
            lo: [ov_bb.lo[0], ov_bb.lo[1]],
            hi: [ov_bb.hi[0], ov_bb.hi[1]],
            max_steiner_multiplier: plan.opts.max_steiner_multiplier,
            strict_mesh: plan.opts.strict_mesh,
            max_cells: None,
        },
        &domain_at,
    )?;
    drop(fills);
    let warnings = mc_rep.warnings;

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
    let canon_sides = |a: Option<u32>, b: Option<u32>| -> (Option<u32>, Option<u32>) {
        if a <= b { (a, b) } else { (b, a) }
    };
    let mut ov_edges_by_pair: HashMap<(Option<u32>, Option<u32>), Vec<usize>> = HashMap::new();
    for (ei, ed) in ov.edges.iter().enumerate() {
        ov_edges_by_pair.entry(canon_sides(ed.left, ed.right)).or_default().push(ei);
    }
    let mut ov_pair_bvh: HashMap<(Option<u32>, Option<u32>), Bvh2> = HashMap::new();
    for (&pair_k, list) in &ov_edges_by_pair {
        if list.len() > 16 {
            let bvh = Bvh2::build(
                list.iter().map(|&ei| BvhItem2 { id: ei, bbox: ov.edges[ei].bbox, min_h: 0.0 }).collect(),
            );
            ov_pair_bvh.insert(pair_k, bvh);
        }
    }
    let all_ov_edge_bvh = Bvh2::build(
        ov.edges.iter().enumerate().map(|(ei, ed)| BvhItem2 { id: ei, bbox: ed.bbox, min_h: 0.0 }).collect(),
    );

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
        let pair_k = canon_sides(c0, c1);
        let best = if let Some(bvh) = ov_pair_bvh.get(&pair_k) {
            bvh.nearest_id(mid, |ei| Some(dist2(ov.project_edge(ei as u32, mid), mid)))
        } else if let Some(list) = ov_edges_by_pair.get(&pair_k) {
            let mut b0: Option<(f64, u32)> = None;
            for &ei in list {
                let d = dist2(ov.project_edge(ei as u32, mid), mid);
                if b0.is_none_or(|(bd, _)| d < bd) {
                    b0 = Some((d, ei as u32));
                }
            }
            b0
        } else {
            None
        }
        .or_else(|| all_ov_edge_bvh.nearest_id(mid, |ei| Some(dist2(ov.project_edge(ei as u32, mid), mid))));
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
        if has_sheet {
            if let Some(dz_min) = plan.opts.dz_min {
                dz = dz.max(dz_min).min(default_h);
            }
            Some(dz)
        } else {
            None
        }
    };

    let z_ratio = plan.opts.z_growth_ratio.unwrap_or(1.35);
    let mut active_layers: Vec<(usize, f64, f64, f64, f64, f64)> = Vec::new();
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
        active_layers.push((k, l.z0, l.z1, dz0, dz1, h_layer));
    }
    if active_layers.is_empty() {
        return Err("the 3D model has no active layers to mesh".into());
    }
    if let Some(max_slabs) = plan.opts.max_z_slabs {
        if max_slabs < active_layers.len() {
            return Err(format!(
                "max_z_slabs ({max_slabs}) is less than the number of active stack layers ({}) required to preserve layer interfaces",
                active_layers.len()
            ));
        }
    }
    let build_layer_zs = |alpha: f64| -> Vec<(usize, Vec<f64>)> {
        active_layers
            .iter()
            .map(|&(k, z0, z1, dz0, dz1, h_layer)| {
                let (d0, d1, hl, rt) = if alpha <= 1.0 {
                    (dz0, dz1, h_layer, z_ratio)
                } else {
                    let span = (z1 - z0).max(0.0);
                    ((dz0 * alpha).min(span), (dz1 * alpha).min(span), (h_layer * alpha).min(span), z_ratio + 0.15 * (alpha - 1.0))
                };
                let zs = if plan.opts.z_coarsen {
                    subdivide_z_interval_far(z0, z1, d0, d1, hl, rt)
                } else {
                    subdivide_z_interval(z0, z1, d0, d1, hl, rt)
                };
                (k, zs)
            })
            .collect()
    };
    let mut per_layer_zs = build_layer_zs(1.0);
    if let Some(max_slabs) = plan.opts.max_z_slabs {
        let count_slabs = |v: &[(usize, Vec<f64>)]| -> usize { v.iter().map(|(_, zs)| zs.len().saturating_sub(1)).sum() };
        if count_slabs(&per_layer_zs) > max_slabs {
            let mut lo_a = 1.0_f64;
            let mut hi_a = 2.0_f64;
            while count_slabs(&build_layer_zs(hi_a)) > max_slabs && hi_a < 1e12 {
                lo_a = hi_a;
                hi_a *= 2.0;
            }
            let mut best_zs = build_layer_zs(hi_a);
            for _ in 0..48 {
                let mid_a = 0.5 * (lo_a + hi_a);
                let cand = build_layer_zs(mid_a);
                if count_slabs(&cand) <= max_slabs {
                    hi_a = mid_a;
                    best_zs = cand;
                } else {
                    lo_a = mid_a;
                }
            }
            per_layer_zs = best_zs;
            while count_slabs(&per_layer_zs) > max_slabs {
                let mut best_drop: Option<(f64, usize, usize)> = None;
                for (li, (_, zs)) in per_layer_zs.iter().enumerate() {
                    if zs.len() <= 2 {
                        continue;
                    }
                    for zi in 1..(zs.len() - 1) {
                        let merged_step = zs[zi + 1] - zs[zi - 1];
                        if best_drop.is_none_or(|b| merged_step < b.0) {
                            best_drop = Some((merged_step, li, zi));
                        }
                    }
                }
                let Some((_, li, zi)) = best_drop else { break };
                per_layer_zs[li].1.remove(zi);
            }
        }
    }

    let mut z_planes: Vec<f64> = Vec::new();
    // For each slab interval [z_planes[s], z_planes[s+1]], the 3D layer index `Some(k)` (or `None` across a void gap)
    let mut slab_layer: Vec<Option<usize>> = Vec::new();
    for (k, zs) in per_layer_zs {
        let z0 = st.layers[k].z0;
        if let Some(&last_z) = z_planes.last() {
            if (z0 - last_z).abs() > tol {
                z_planes.push(z0);
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
    let n_active_slabs = slab_layer.iter().filter(|x| x.is_some()).count();
    let n2 = m2.n_points();

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

    let active_tris_per_layer: Vec<usize> = (0..st.layers.len())
        .map(|k| {
            tri_ov_cell
                .iter()
                .filter(|&&ov_c| st.layers[k].region_of[ov_c].is_some_and(|r| meshable[r as usize]))
                .count()
        })
        .collect();

    if !plan.opts.z_coarsen {
        if let Some(max_cells) = eff_max_cells {
            let mut projected_tets = 0usize;
            for &layer_opt in &slab_layer {
                if let Some(k) = layer_opt {
                    projected_tets = projected_tets.saturating_add(3usize.saturating_mul(active_tris_per_layer[k]));
                }
            }
            if projected_tets > max_cells {
                return Err(format!(
                    "3D sheet mesh exceeds max_cells limit: n_tri2={ntri2}, n_slabs={n_active_slabs}, projected_tets={projected_tets} > max_cells={max_cells}; increase h / dz_min / z_growth_ratio, lower max_z_slabs, enable z_coarsen, or raise max_cells"
                ));
            }
        }
        // Default exact prismatic Dompierre path (bit-identical to baseline)
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
        let mut g = MeshedGeom { mesh, labels, eval: std::sync::Arc::new(ev), warnings };
        if matches!(plan.pml, PmlMode::Structured) && !g.eval.pml_regions().is_empty() {
            g = crate::meshing::pml_geom3::structured_pml_3d(g)?.0;
        }
        if plan.opts.curved {
            let gcx = g.eval.complex();
            curve_sheet_mesh_edges(&mut g.mesh, gcx, pt2, &bnd_edges2, &z_planes, &used_v, n2);
        }
        return Ok((g, None));
    }

    // `z_coarsen = true`: Conforming pure-tetrahedral vertical coarsening away from sheet planes
    let mut sheet_plane_info: Vec<(usize, f64, f64)> = Vec::new();
    for (zp, &z) in z_planes.iter().enumerate() {
        if let Some(hs) = sheet_dz_at(z, h_global) {
            sheet_plane_info.push((zp, z, hs));
        }
    }
    let pts2: Vec<V2> = (0..n2 as u32).map(pt2).collect();
    let mut frozen_v = vec![false; n2];
    for &((a, b), ref sides) in &edge_list {
        if sides.len() != 2 {
            frozen_v[a as usize] = true;
            frozen_v[b as usize] = true;
            continue;
        }
        let (c0, c1) = (sides[0], sides[1]);
        if c0 != c1 {
            let diff_layer = st.layers.iter().any(|l| l.region_of[c0] != l.region_of[c1]);
            let diff_non_sheet_cap = (0..nz).any(|zp| {
                if sheet_plane_info.iter().any(|x| x.0 == zp) {
                    return false;
                }
                let f0 = cap_at.get(&(c0 as u32, zp)).filter(|&&fi| !cx.faces[fi].is_sheet());
                let f1 = cap_at.get(&(c1 as u32, zp)).filter(|&&fi| !cx.faces[fi].is_sheet());
                match (f0, f1) {
                    (None, None) => false,
                    (Some(&i0), Some(&i1)) => cx.faces[i0].name != cx.faces[i1].name,
                    _ => true,
                }
            });
            if diff_layer || diff_non_sheet_cap {
                frozen_v[a as usize] = true;
                frozen_v[b as usize] = true;
            }
        }
    }
    for &(a, b, oe) in &bnd_edges2 {
        if (0..st.layers.len()).any(|k| lat_at.contains_key(&(oe, k))) {
            frozen_v[a as usize] = true;
            frozen_v[b as usize] = true;
        }
    }

    // Determine coarsening level `level_of_zp[zp]` (`0` inside `|z - z_sheet| <= 4 * h_sheet` & >= 2 slabs per side)
    let mut level_of_zp = if sheet_plane_info.is_empty() {
        vec![0usize; nz]
    } else {
        let mut lvl = vec![usize::MAX / 2; nz];
        for &(zp_s, z_s, h_s) in &sheet_plane_info {
            lvl[zp_s] = 0;
            for s in 0..nslab {
                if slab_layer[s].is_none() {
                    continue;
                }
                let d_min = (z_planes[s] - z_s).abs().min((z_planes[s + 1] - z_s).abs());
                let slab_dist = if s >= zp_s { s - zp_s } else { zp_s - 1 - s };
                if slab_dist < 2 || d_min < 4.0 * h_s - tol {
                    lvl[s] = 0;
                    lvl[s + 1] = 0;
                }
            }
        }
        for s in 0..nslab {
            if slab_layer[s].is_some() {
                lvl[s + 1] = lvl[s + 1].min(lvl[s].saturating_add(1));
            }
        }
        for s in (0..nslab).rev() {
            if slab_layer[s].is_some() {
                lvl[s] = lvl[s].min(lvl[s + 1].saturating_add(1));
            }
        }
        for v in &mut lvl {
            if *v >= usize::MAX / 4 {
                *v = 0;
            }
        }
        lvl
    };

    if let Some(max_cells) = eff_max_cells {
        let mut min_uncoarsened_tets = 0usize;
        for (s, &layer_opt) in slab_layer.iter().enumerate() {
            if let Some(k) = layer_opt {
                if level_of_zp[s] == 0 && level_of_zp[s + 1] == 0 {
                    min_uncoarsened_tets = min_uncoarsened_tets.saturating_add(3usize.saturating_mul(active_tris_per_layer[k]));
                }
            }
        }
        if min_uncoarsened_tets > max_cells {
            return Err(format!(
                "3D sheet mesh exceeds max_cells limit: n_tri2={ntri2}, n_slabs={n_active_slabs}, projected_tets={min_uncoarsened_tets} > max_cells={max_cells}; increase h / dz_min / z_growth_ratio, lower max_z_slabs, or raise max_cells"
            ));
        }
    }

    let max_lvl = *level_of_zp.iter().max().unwrap_or(&0);

    // Level 0: fine 2D triangulation `M_0` oriented strictly CCW
    let mut tris0 = Vec::with_capacity(ntri2);
    for (ti, &ov_c) in tri_ov_cell.iter().enumerate().take(ntri2) {
        let mut v = [m2.cell(ti)[0], m2.cell(ti)[1], m2.cell(ti)[2]];
        if cross2_pts(pts2[v[0] as usize], pts2[v[1] as usize], pts2[v[2] as usize]) < 0.0 {
            v.swap(1, 2);
        }
        tris0.push((v, ov_c));
    }
    let mut levels_2d = Vec::with_capacity(max_lvl + 1);
    levels_2d.push(Level2D { tris: tris0 });
    let mut transitions = Vec::with_capacity(max_lvl);

    for m in 0..max_lvl {
        let mut h_layer_max = 0.0f64;
        for &layer_opt in &slab_layer {
            let Some(k) = layer_opt else { continue };
            let hl = st.layers[k]
                .region_of
                .iter()
                .flatten()
                .filter(|&&r| meshable[r as usize])
                .map(|&r| region_h[r as usize])
                .fold(h_global, f64::min);
            h_layer_max = h_layer_max.max(hl);
        }
        if h_layer_max <= 0.0 {
            h_layer_max = h_global;
        }
        let mut z_dist_min = f64::INFINITY;
        for (zp, &lv) in level_of_zp.iter().enumerate() {
            if lv == m + 1 {
                for &(_, z_s, _) in &sheet_plane_info {
                    z_dist_min = z_dist_min.min((z_planes[zp] - z_s).abs());
                }
            }
        }
        let h_target = if z_dist_min.is_finite() {
            (1.15 * z_dist_min).min(0.95 * h_layer_max)
        } else {
            0.95 * h_layer_max
        };
        let (next_lvl, tr) = coarsen_2d_once(&levels_2d[m], &pts2, &frozen_v, h_target, h_layer_max);
        levels_2d.push(next_lvl);
        transitions.push(tr);
    }

    // Clamp `level_of_zp` if coarsening plateaued early (identical triangle count)
    let mut effective_max = max_lvl;
    while effective_max > 0 && levels_2d[effective_max].tris.len() == levels_2d[effective_max - 1].tris.len() {
        effective_max -= 1;
    }
    for v in &mut level_of_zp {
        if *v > effective_max {
            *v = effective_max;
        }
    }

    if let Some(max_cells) = eff_max_cells {
        let mut projected_tets = 0usize;
        for (s, &layer_opt) in slab_layer.iter().enumerate() {
            let Some(k) = layer_opt else { continue };
            let l = &st.layers[k];
            let (ls, lt) = (level_of_zp[s], level_of_zp[s + 1]);
            if ls == lt {
                for &(_, ov_c) in &levels_2d[ls].tris {
                    if l.region_of[ov_c].is_some_and(|r| meshable[r as usize]) {
                        projected_tets = projected_tets.saturating_add(3);
                    }
                }
            } else {
                let m = ls.min(lt);
                let tr = &transitions[m];
                for (ti, &(_, ov_c)) in levels_2d[m].tris.iter().enumerate() {
                    if l.region_of[ov_c].is_some_and(|r| meshable[r as usize]) {
                        projected_tets = projected_tets.saturating_add(if tr.tri_collapse[ti].is_some() { 2 } else { 3 });
                    }
                }
            }
        }
        if projected_tets > max_cells {
            return Err(format!(
                "3D sheet mesh exceeds max_cells limit: n_tri2={ntri2}, n_slabs={n_active_slabs}, projected_tets={projected_tets} > max_cells={max_cells}; increase h / dz_min / z_growth_ratio, lower max_z_slabs, or raise max_cells"
            ));
        }
    }

    // Mark active 3D vertices `(v2, zp)` used across constant or coarsening slabs
    let mut used_v = vec![u32::MAX; nz * n2];
    for (s, &layer_opt) in slab_layer.iter().enumerate() {
        let Some(k) = layer_opt else { continue };
        let l = &st.layers[k];
        let (ls, lt) = (level_of_zp[s], level_of_zp[s + 1]);
        if ls == lt {
            for &(v, ov_c) in &levels_2d[ls].tris {
                let Some(r) = l.region_of[ov_c] else { continue };
                if !meshable[r as usize] {
                    continue;
                }
                for &w in &v {
                    used_v[s * n2 + w as usize] = 0;
                    used_v[(s + 1) * n2 + w as usize] = 0;
                }
            }
        } else {
            let (m, f_plane, c_plane) = if ls < lt { (ls, s, s + 1) } else { (lt, s + 1, s) };
            let tr = &transitions[m];
            for &(v, ov_c) in &levels_2d[m].tris {
                let Some(r) = l.region_of[ov_c] else { continue };
                if !meshable[r as usize] {
                    continue;
                }
                for &w in &v {
                    used_v[f_plane * n2 + w as usize] = 0;
                    used_v[c_plane * n2 + tr.map_v[w as usize] as usize] = 0;
                }
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
                let p2 = pts2[v2];
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

    let mut cells: Vec<u32> = Vec::new();
    let mut cell_labels: Vec<i32> = Vec::new();
    let mut push_oriented_tet = |mut tet: [u32; 4], dom_id: i32| {
        let (p0, p1, p2, p3) = (pt3(tet[0]), pt3(tet[1]), pt3(tet[2]), pt3(tet[3]));
        if det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0)) < 0.0 {
            tet.swap(1, 2);
        }
        cells.extend_from_slice(&tet);
        cell_labels.push(dom_id);
    };

    for (s, &layer_opt) in slab_layer.iter().enumerate() {
        let Some(k) = layer_opt else { continue };
        let l = &st.layers[k];
        let (ls, lt) = (level_of_zp[s], level_of_zp[s + 1]);
        if ls == lt {
            for &(mut v, ov_c) in &levels_2d[ls].tris {
                let Some(r) = l.region_of[ov_c] else { continue };
                if !meshable[r as usize] {
                    continue;
                }
                let dom_id = labels.domains[r as usize].id;
                v.sort_unstable();
                let (b0, b1, b2) = (vid(v[0], s), vid(v[1], s), vid(v[2], s));
                let (t0, t1, t2) = (vid(v[0], s + 1), vid(v[1], s + 1), vid(v[2], s + 1));
                for tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                    push_oriented_tet(tet, dom_id);
                }
            }
        } else {
            let (m, f_plane, c_plane) = if ls < lt { (ls, s, s + 1) } else { (lt, s + 1, s) };
            let tr = &transitions[m];
            for (ti, &(v, ov_c)) in levels_2d[m].tris.iter().enumerate() {
                let Some(r) = l.region_of[ov_c] else { continue };
                if !meshable[r as usize] {
                    continue;
                }
                let dom_id = labels.domains[r as usize].id;
                if let Some((u, v_rem, c)) = tr.tri_collapse[ti] {
                    let c_top = tr.map_v[c as usize];
                    let (uf, vf, cf) = (vid(u, f_plane), vid(v_rem, f_plane), vid(c, f_plane));
                    let (uc, cc) = (vid(u, c_plane), vid(c_top, c_plane));
                    push_oriented_tet([uf, vf, cf, cc], dom_id);
                    push_oriented_tet([uf, vf, cc, uc], dom_id);
                } else {
                    let mut x = v;
                    x.sort_by_key(|&w| tr.rank[w as usize]);
                    let (b0, b1, b2) = (vid(x[0], f_plane), vid(x[1], f_plane), vid(x[2], f_plane));
                    let (t0, t1, t2) = (
                        vid(tr.map_v[x[0] as usize], c_plane),
                        vid(tr.map_v[x[1] as usize], c_plane),
                        vid(tr.map_v[x[2] as usize], c_plane),
                    );
                    for tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                        push_oriented_tet(tet, dom_id);
                    }
                }
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

    // 1. Horizontal cap facets on each `zp` using `levels_2d[level_of_zp[zp]]`
    for zp in 0..nz {
        let m = level_of_zp[zp];
        for &(v, ov_c) in &levels_2d[m].tris {
            if let Some(&fi) = cap_at.get(&(ov_c as u32, zp)) {
                let tri = [vid(v[0], zp), vid(v[1], zp), vid(v[2], zp)];
                if tri.iter().all(|&w| w != u32::MAX) {
                    push_oriented_tri(tri, fi);
                }
            }
        }
    }

    // 2. Vertical lateral facets matching exact slab `(p, q, f_plane, c_plane)` diagonal
    for &(a, b, oe) in &bnd_edges2 {
        for (s, &layer_opt) in slab_layer.iter().enumerate().take(nslab) {
            let Some(k) = layer_opt else { continue };
            if let Some(&fi) = lat_at.get(&(oe, k)) {
                let (ls, lt) = (level_of_zp[s], level_of_zp[s + 1]);
                let (p, q, f_plane, c_plane) = if ls == lt {
                    (a.min(b), a.max(b), s, s + 1)
                } else {
                    let (m, fp, cp) = if ls < lt { (ls, s, s + 1) } else { (lt, s + 1, s) };
                    let tr = &transitions[m];
                    let (p, q) = if tr.rank[a as usize] <= tr.rank[b as usize] { (a, b) } else { (b, a) };
                    (p, q, fp, cp)
                };
                let (bp, bq, tp, tq) = (vid(p, f_plane), vid(q, f_plane), vid(p, c_plane), vid(q, c_plane));
                if [bp, bq, tp, tq].iter().all(|&w| w != u32::MAX) {
                    push_oriented_tri([bp, bq, tq], fi);
                    push_oriented_tri([bp, tq, tp], fi);
                }
            }
        }
    }

    let mesh = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);
    let mut g = MeshedGeom { mesh, labels, eval: std::sync::Arc::new(ev), warnings };
    if matches!(plan.pml, PmlMode::Structured) && !g.eval.pml_regions().is_empty() {
        g = crate::meshing::pml_geom3::structured_pml_3d(g)?.0;
    }
    if plan.opts.curved {
        let gcx = g.eval.complex();
        curve_sheet_mesh_edges(&mut g.mesh, gcx, pt2, &bnd_edges2, &z_planes, &used_v, n2);
    }
    Ok((g, None))
}
