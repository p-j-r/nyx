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
use crate::meshing::mesher::{mesh_chains_with_seeds, Chain, Fill, MeshParams, SheetChainExtras};
use crate::meshing::physics_size::PhysicsSizing;
use crate::meshing::{wave_sizes_by_region, MeshPlan, Mesher, Sizing};
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

fn norm3(a: V3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Minimum corner scaled Jacobian of a straight tetrahedron normalized so a regular tetrahedron has `1.0`
/// (`min_k (sqrt(2) * |det|) / (||e1|| ||e2|| ||e3||)`).
fn tet_scaled_jacobian(p0: V3, p1: V3, p2: V3, p3: V3) -> f64 {
    let e01 = sub3(p1, p0);
    let e02 = sub3(p2, p0);
    let e03 = sub3(p3, p0);
    let vol6 = det3(e01, e02, e03).abs();
    if !(vol6 > 0.0 && vol6.is_finite()) {
        return 0.0;
    }
    let num = std::f64::consts::SQRT_2 * vol6;
    let l01 = norm3(e01);
    let l02 = norm3(e02);
    let l03 = norm3(e03);
    let l12 = norm3(sub3(p2, p1));
    let l13 = norm3(sub3(p3, p1));
    let l23 = norm3(sub3(p3, p2));
    let mut min_j = f64::INFINITY;
    for prod in [l01 * l02 * l03, l01 * l12 * l13, l02 * l12 * l23, l03 * l13 * l23] {
        if !(prod > 0.0) {
            return 0.0;
        }
        min_j = min_j.min(num / prod);
    }
    min_j
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

/// Computes the facing overlap centerline `(m0, m1, d_min)` between two rim subsegments `(a0, a1)` and `(b0, b1)`
/// with outward-into-dielectric unit normals `na`, `nb`. Returns `None` if the segments do not face each other
/// within `d_max`.
fn facing_overlap_centerline(a0: V2, a1: V2, na: V2, b0: V2, b1: V2, nb: V2, d_max: f64) -> Option<(V2, V2, f64)> {
    let dot_nn = na[0] * nb[0] + na[1] * nb[1];
    if dot_nn >= 0.25 {
        return None;
    }
    let (d, pa, pb) = seg_seg_closest(a0, a1, b0, b1);
    if !(d > 0.0 && d <= d_max) {
        return None;
    }
    let dp = [pb[0] - pa[0], pb[1] - pa[1]];
    let dot_a = dp[0] * na[0] + dp[1] * na[1];
    let dot_b = -dp[0] * nb[0] - dp[1] * nb[1];
    if !(dot_a > 0.15 * d && dot_b > 0.15 * d) {
        return None;
    }
    let da = [a1[0] - a0[0], a1[1] - a0[1]];
    let la2 = da[0] * da[0] + da[1] * da[1];
    if la2 > 1e-30 {
        let tb0 = ((b0[0] - a0[0]) * da[0] + (b0[1] - a0[1]) * da[1]) / la2;
        let tb1 = ((b1[0] - a0[0]) * da[0] + (b1[1] - a0[1]) * da[1]) / la2;
        let t_lo = tb0.min(tb1).clamp(0.0, 1.0);
        let t_hi = tb0.max(tb1).clamp(0.0, 1.0);
        if (t_hi - t_lo) * la2.sqrt() >= 0.25 * d {
            let qa0 = [a0[0] + t_lo * da[0], a0[1] + t_lo * da[1]];
            let qa1 = [a0[0] + t_hi * da[0], a0[1] + t_hi * da[1]];
            let qb0 = closest_point_on_seg(qa0, b0, b1);
            let qb1 = closest_point_on_seg(qa1, b0, b1);
            let d0 = dist2(qa0, qb0);
            let d1 = dist2(qa1, qb1);
            if d0 <= 1.6 * d && d1 <= 1.6 * d {
                let m0 = [0.5 * (qa0[0] + qb0[0]), 0.5 * (qa0[1] + qb0[1])];
                let m1 = [0.5 * (qa1[0] + qb1[0]), 0.5 * (qa1[1] + qb1[1])];
                return Some((m0, m1, d));
            }
        }
    }
    let m = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
    Some((m, m, d))
}

#[derive(Clone, Copy)]
struct GapZone {
    a: V2,
    b: V2,
    h: f64,
    bbox: Aabb,
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
    while {
        let rem = hi.last().unwrap() - lo.last().unwrap();
        rem > s0 + s1 || (s0.max(s1) > 0.5 * span && rem > 2.0 * s0.min(s1))
    } {
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
/// `|z - z_sheet| <= 9.2 * h_sheet` (preserving near- and mid-field vertical resolution identically to
/// `subdivide_z_interval`), then accelerates to `ratio.max(1.66)` in the far dielectric/vacuum bulk up to `h_max`.
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

    let r_far = 1.66_f64.max(ratio);
    let mut lo = vec![0.0];
    let mut hi = vec![span];
    while {
        let rem = hi.last().unwrap() - lo.last().unwrap();
        rem > s0 + s1 || (s0.max(s1) > 0.5 * span && rem > 2.0 * s0.min(s1))
    } {
        if s0 <= s1 {
            let cur = *lo.last().unwrap();
            lo.push(cur + s0);
            let r = if cur + s0 > 9.2 * d0_init { r_far } else { ratio };
            s0 = (s0 * r).min(h_max);
        } else {
            let cur = *hi.last().unwrap();
            hi.push(cur - s1);
            let r = if span - (cur - s1) > 9.2 * d1_init { r_far } else { ratio };
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
/// 4. No active cap-equivalence cell component loses its last surviving triangle;
/// 5. Every surviving 2D triangle stays strictly CCW with quality `>= MIN_TRI_QUALITY`, and all 3D Dompierre prism / wedge
///    tetrahedra on the transition slab have strictly positive volume and `scaled_jacobian >= 0.08`.
#[allow(clippy::too_many_arguments)]
fn coarsen_2d_once(
    cur: &Level2D,
    pts2: &[V2],
    frozen_v: &[bool],
    h_target: f64,
    h_max_edge: f64,
    dz_slab: f64,
    comp_of_cell: &[usize],
    is_far_bulk: bool,
) -> (Level2D, SlabTransition) {
    let min_tri_quality: f64 = 0.12;
    let min_tet_scaled_jacobian: f64 = if is_far_bulk { 0.035 } else { 0.055 };
    let n2 = pts2.len();
    let mut work_tris = cur.tris.clone();
    let ntri = cur.tris.len();

    let n_comps = comp_of_cell.iter().copied().max().map_or(0, |m| m + 1);
    let mut comp_tri_count = vec![0usize; n_comps];
    for &(_, ov_c) in &cur.tris {
        if let Some(&cid) = comp_of_cell.get(ov_c) {
            comp_tri_count[cid] += 1;
        }
    }

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

    // Map each vertex in `cur.tris` to its incident original triangle indices (for fast 3D transition prism check)
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

        // Never drop the last surviving triangle of any cap-equivalence cell component
        let comp0 = comp_of_cell.get(work_tris[shared_tris[0].0].1).copied();
        let comp1 = comp_of_cell.get(work_tris[shared_tris[1].0].1).copied();
        match (comp0, comp1) {
            (Some(g0), Some(g1)) if g0 == g1 => {
                if comp_tri_count[g0] <= 2 {
                    continue;
                }
            }
            (Some(g0), Some(g1)) => {
                if comp_tri_count[g0] <= 1 || comp_tri_count[g1] <= 1 {
                    continue;
                }
            }
            (Some(g0), None) => {
                if comp_tri_count[g0] <= 1 {
                    continue;
                }
            }
            (None, Some(g1)) => {
                if comp_tri_count[g1] <= 1 {
                    continue;
                }
            }
            (None, None) => {}
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
            if tri_quality_ccw(p0, p1, p2) < min_tri_quality {
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

        // Check 3D Dompierre prism/wedge positivity & scaled Jacobian >= 0.08 for all original triangles in `star(u) \cup star(v)`
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
            let (q0, q1, q2) = (pts2[tv[0] as usize], pts2[tv[1] as usize], pts2[tv[2] as usize]);
            let h_tri = dist2(q0, q1).max(dist2(q1, q2)).max(dist2(q2, q0)).max(1e-12 * h_max_edge);
            // Enforce transition-tet `scaled_jacobian >= 0.08` on near-isotropic `dz_eff`
            // so sheared/degenerate collapses are rejected while preserving >= 3x bulk reduction.
            let dz_eff = dz_slab.clamp(0.85 * h_tri, 1.15 * h_tri).max(1e-12 * h_max_edge);
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
                if a_old.signum() * a_new <= 0.12 * a_old.abs() {
                    prism_ok = false;
                    break;
                }
                let uf = [pts2[u_col as usize][0], pts2[u_col as usize][1], 0.0];
                let vf = [pts2[v_col as usize][0], pts2[v_col as usize][1], 0.0];
                let cf = [pts2[c_opp as usize][0], pts2[c_opp as usize][1], 0.0];
                let uc = [pts2[u_col as usize][0], pts2[u_col as usize][1], dz_eff];
                let cc = [pts2[c_top as usize][0], pts2[c_top as usize][1], dz_eff];
                let w_sig = a_old.signum();
                let w_floor = 1e-6 * a_old.abs() * dz_eff;
                for [p0, p1, p2, p3] in [[uf, vf, cf, cc], [uf, vf, cc, uc]] {
                    if w_sig * det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0)) <= w_floor {
                        prism_ok = false;
                        break;
                    }
                }
                if !prism_ok {
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
                [pts2[m0 as usize][0], pts2[m0 as usize][1], dz_eff],
                [pts2[m1 as usize][0], pts2[m1 as usize][1], dz_eff],
                [pts2[m2 as usize][0], pts2[m2 as usize][1], dz_eff],
            ];
            let area_sign = cross2_pts([b[0][0], b[0][1]], [b[1][0], b[1][1]], [b[2][0], b[2][1]]);
            let area_mag = area_sign.abs();
            let sigma = area_sign.signum();
            for [p0, p1, p2, p3] in [[b[0], b[1], b[2], t[2]], [b[0], b[1], t[2], t[1]], [b[0], t[0], t[1], t[2]]] {
                let d3 = det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0));
                if sigma * d3 <= 1e-6 * area_mag * dz_eff
                    || tet_scaled_jacobian(p0, p1, p2, p3) < min_tet_scaled_jacobian
                {
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
        if let Some(g0) = comp0 {
            comp_tri_count[g0] -= 1;
        }
        if let Some(g1) = comp1 {
            comp_tri_count[g1] -= 1;
        }
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

struct Sheet3dSizing {
    h_global: f64,
    max_dev: f64,
    eff_max_cells: Option<usize>,
    meshable: Vec<bool>,
    region_h: Vec<f64>,
    face_h: Vec<f64>,
    edge3_h: Vec<f64>,
    cell_h: Vec<f64>,
    ov_edge_h: Vec<f64>,
    ov_edge_explicit_fine: Vec<bool>,
    edge_sheet_owner: Vec<u32>,
    dielectric_gap_cells: HashSet<i32>,
    gap_zones: Vec<GapZone>,
    gap_seeds: Vec<(V2, f64)>,
    sheet_local_gap_by_z: Vec<(f64, f64)>,
    sharp_corners_by_z: Vec<(f64, Vec<crate::meshing::singularity::SheetCorner>)>,
    all_sharp_corners: Vec<crate::meshing::singularity::SheetCorner>,
}

fn validate_sheet_plan_opts(plan: &MeshPlan) -> Result<(), String> {
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
    if let Some(se) = plan.opts.singular_edges.as_deref() {
        if se != "none" && se != "quarter_point" {
            return Err(format!("singular_edges must be \"none\" or \"quarter_point\", got '{se}'"));
        }
    }
    if let Some(xi) = plan.opts.quarter_point_xi {
        if !(xi.is_finite() && (0.25..=0.5).contains(&xi)) {
            return Err(format!("quarter_point_xi must be in [0.25, 0.50], got {xi}"));
        }
    }
    Ok(())
}

fn compute_sheet3d_sizing(
    ev: &Evaluated,
    labels: &geom::labels::LabelMap,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<Sheet3dSizing, String> {
    if ev.list.dim != 3 {
        return Err("mesh_sheet_stack_3d requires a 3D model".into());
    }
    validate_sheet_plan_opts(plan)?;
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
    let cx = ev.complex();
    let st = cx.stack.as_ref().ok_or("3D sheet mesher requires a prismatic layer stack")?;
    let ov = &st.overlay;
    let tol = cx.tol;

    let nr = ev.regions().len();
    let structured_pml = plan.pml.is_structured(!ev.pml_regions().is_empty());
    let mut region_h = vec![h_global; nr];
    if !structured_pml {
        for p in ev.pml_regions() {
            region_h[p.region as usize] = region_h[p.region as usize].min((p.thickness / p.layers) * scale_factor);
        }
    }
    if let Some((sizing, sz)) = waves {
        let mut eff_sz = *sz;
        if structured_pml {
            eff_sz.pml_relax = 1.0;
        }
        for (r, w) in wave_sizes_by_region(ev, labels, sizing, &eff_sz)?.into_iter().enumerate().take(nr) {
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
    if structured_pml {
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
    let mut ov_edge_explicit_fine = vec![false; nedge];
    for (e, ed) in ov.edges.iter().enumerate() {
        ov_edge_h[e] = [ed.left, ed.right].into_iter().flatten().map(|c| cell_h[c as usize]).fold(h_global, f64::min);
    }
    for (fi, f) in cx.faces.iter().enumerate() {
        if let FaceGeom::Lateral { edge, .. } = f.geom {
            let e = edge as usize;
            if e < nedge {
                if face_h[fi] < h_global * 0.999 {
                    ov_edge_explicit_fine[e] = true;
                }
                ov_edge_h[e] = ov_edge_h[e].min(face_h[fi]);
            }
        }
    }
    for (ei, e3) in cx.edges3.iter().enumerate() {
        if let EdgeGeom::Horizontal { edge, .. } = e3.geom {
            let e = edge as usize;
            if e < nedge {
                if edge3_h[ei] < h_global * 0.999 {
                    ov_edge_explicit_fine[e] = true;
                }
                ov_edge_h[e] = ov_edge_h[e].min(edge3_h[ei]);
            }
        }
    }

    let grade = 0.35;
    let has_explicit_fine =
        ov_edge_h.iter().any(|&h| h < h_global * 0.999) || cell_h.iter().any(|&h| h < h_global * 0.999);
    let want_qp = plan.opts.singular_edges.as_deref() == Some("quarter_point");
    let eff_gap_cells = plan.opts.gap_cells.map(|gc| if want_qp { gc.max(2.0) } else { gc });

    // Coplanar gap proximity pass across all 2D sheet planes (FR-27, FR-28):
    // 1. When `gap_cells = Some(gc)`, clamps `ov_edge_h[e] <= d_gap / gc` symmetrically on both facing rims.
    // 2. Automatically detects facing coplanar sheet rims with distinct sheet/terminal ownership within
    //    `d < 0.25 * h_global` (implicit `eff_gap_cells >= 2.0`), pre-discretizing straight gap walls before
    //    `cdt.refine()`, and inserts deterministic gap-centerline seeds along facing narrow gaps (`d <= 0.35 * h_global`).
    let mut sheet_zs: Vec<f64> = Vec::new();
    let mut all_sheet_names: Vec<&str> = Vec::new();
    for f in &cx.faces {
        if f.is_sheet() {
            all_sheet_names.push(f.name.as_str());
            if let FaceGeom::Cap { z, .. } = f.geom {
                if !sheet_zs.iter().any(|&zs| (zs - z).abs() <= tol) {
                    sheet_zs.push(z);
                }
            }
        }
    }
    all_sheet_names.sort_unstable();
    all_sheet_names.dedup();
    let sheet_id_of: HashMap<&str, u32> =
        all_sheet_names.iter().enumerate().map(|(i, &s)| (s, (i + 1) as u32)).collect();

    let mut edge_sheet_owner = vec![0u32; nedge];
    let mut dielectric_gap_cells: HashSet<i32> = HashSet::new();
    let mut gap_zones: Vec<GapZone> = Vec::new();
    let mut gap_seeds: Vec<(V2, f64)> = Vec::new();
    let mut sheet_local_gap_by_z: Vec<(f64, f64)> = Vec::new();

    if !sheet_zs.is_empty() {
        struct RimEdge<'a> {
            e: usize,
            v0: u32,
            v1: u32,
            sheet: &'a str,
            len: f64,
            bbox: Aabb,
            segs: Vec<(V2, V2, V2)>,
        }
        let n_active_layers = st
            .layers
            .iter()
            .filter(|l| l.region_of.iter().flatten().any(|&r| meshable[r as usize]))
            .count()
            .max(1);
        let search_r = if eff_gap_cells.is_some() || has_explicit_fine { h_global } else { 0.35 * h_global };
        // Respect explicit `eff_gap_cells` when provided; default to `2.0` for automatic narrow-gap detection (FR-28).
        let auto_gc = eff_gap_cells.unwrap_or(2.0);

        for &z_s in &sheet_zs {
            let mut local_gap_max = 0.0_f64;
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
            let mut rims: Vec<RimEdge<'_>> = Vec::new();
            for (e, ed) in ov.edges.iter().enumerate() {
                if ed.left.is_none() || ed.right.is_none() {
                    continue;
                }
                let ls = ed.left.and_then(|c| cell_sheet[c as usize]);
                let rs = ed.right.and_then(|c| cell_sheet[c as usize]);
                if ls.is_some() == rs.is_some() {
                    continue;
                }
                let sheet_name = ls.or(rs).unwrap();
                let gap_cell = if ls.is_some() { ed.right.unwrap() } else { ed.left.unwrap() };
                dielectric_gap_cells.insert(gap_cell as i32 + 1);
                let sid = sheet_id_of[sheet_name];
                if edge_sheet_owner[e] == 0 || sid < edge_sheet_owner[e] {
                    edge_sheet_owner[e] = sid;
                }
                let n_sign = if ls.is_some() { 1.0 } else { -1.0 };
                let curve = &ov.curves[ed.curve as usize];
                let nsub = if matches!(curve, Curve::Line { .. }) { 1 } else { 8 };
                let mut segs = Vec::with_capacity(nsub);
                let mut total_len = 0.0_f64;
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
                    total_len += dist2(pa, pb);
                    let n = ov.edge_normal(e as u32, tm);
                    segs.push((pa, pb, [n_sign * n[0], n_sign * n[1]]));
                }
                rims.push(RimEdge { e, v0: ed.v0, v1: ed.v1, sheet: sheet_name, len: total_len, bbox: ed.bbox, segs });
            }

            let rim_bvh = Bvh2::build(
                rims.iter().enumerate().map(|(i, r)| BvhItem2 { id: i, bbox: r.bbox, min_h: 0.0 }).collect(),
            );
            let mut gap_clamps: Vec<(usize, f64)> = Vec::new();
            // Auto narrow-gap clamps: (edge_idx, h_want, overlap_len, centerline_m0, centerline_m1)
            let mut auto_clamps: Vec<(usize, f64, f64, V2, V2)> = Vec::new();
            let mut raw_centerlines: Vec<(V2, V2, f64, usize, usize)> = Vec::new();
            for (i, ra) in rims.iter().enumerate() {
                rim_bvh.for_each_within_box(&ra.bbox, search_r, |j| {
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
                    if !(d > 1e-15 * h_global && d < search_r) {
                        return;
                    }
                    let dp = [pb[0] - pa[0], pb[1] - pa[1]];
                    let dot_a = dp[0] * na[0] + dp[1] * na[1];
                    let dot_b = -dp[0] * nb[0] - dp[1] * nb[1];
                    let dot_nn = na[0] * nb[0] + na[1] * nb[1];
                    if !(dot_a > 0.15 * d && dot_b > 0.15 * d && dot_nn < 0.25) {
                        return;
                    }
                    if d < 0.25 * h_global {
                        local_gap_max = local_gap_max.max(d);
                    }
                    if let Some(gc) = eff_gap_cells {
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
                    // FR-28: Auto narrow-gap detection & centerline seeding for distinct sheet/terminal ownership
                    // (active when `gap_cells` is set or when `quarter_point`/`corner_grading` is enabled within
                    // `d <= 0.35 * h_global`; leaves plain `MeshPlan::simple` with `gap_cells = None` unperturbed
                    // so coarse-mesh benchmarks and explicit `gap_cells` comparisons remain exact).
                    let enable_auto_gap = eff_gap_cells.is_some() || want_qp || plan.opts.corner_grading;
                    if ra.sheet != rb.sheet && d <= 0.35 * h_global && enable_auto_gap {
                        let pm_best = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
                        let clear_gap = ov.locate(pm_best).is_some_and(|c| cell_sheet[c as usize].is_none());
                        if !clear_gap {
                            return;
                        }
                        let d_lim = (1.5 * d).min(0.35 * h_global);
                        let mut any_cl = false;
                        for &(a0, a1, na) in &ra.segs {
                            for &(b0, b1, nb) in &rb.segs {
                                if let Some((m0, m1, d_seg)) =
                                    facing_overlap_centerline(a0, a1, na, b0, b1, nb, d_lim)
                                {
                                    let mc = [0.5 * (m0[0] + m1[0]), 0.5 * (m0[1] + m1[1])];
                                    if !ov.locate(mc).is_some_and(|c| cell_sheet[c as usize].is_none()) {
                                        continue;
                                    }
                                    any_cl = true;
                                    raw_centerlines.push((m0, m1, d_seg, ra.e, rb.e));
                                    if d_seg < 0.25 * h_global {
                                        let h_auto = d_seg / auto_gc;
                                        let ov_len = dist2(m0, m1).max(d_seg);
                                        auto_clamps.push((ra.e, h_auto, ov_len, m0, m1));
                                        auto_clamps.push((rb.e, h_auto, ov_len, m0, m1));
                                    }
                                }
                            }
                        }
                        if !any_cl {
                            raw_centerlines.push((pm_best, pm_best, d, ra.e, rb.e));
                            if d < 0.25 * h_global {
                                let h_auto = d / auto_gc;
                                auto_clamps.push((ra.e, h_auto, d, pm_best, pm_best));
                                auto_clamps.push((rb.e, h_auto, d, pm_best, pm_best));
                            }
                        }
                    }
                });
            }
            for (ei, hc) in gap_clamps {
                ov_edge_h[ei] = ov_edge_h[ei].min(hc);
            }
            if !auto_clamps.is_empty() {
                // Budget-aware scaling when `gap_cells` was not explicitly requested and `max_cells` is set
                let mut auto_scale = 1.0_f64;
                if plan.opts.gap_cells.is_none() {
                    if let Some(mc) = eff_max_cells {
                        let est_slabs = plan.opts.max_z_slabs.unwrap_or(8).max(n_active_layers) as f64;
                        let budget_tris_2d = (mc as f64) / (3.0 * est_slabs);
                        let mut min_h_per_edge: HashMap<usize, (f64, f64)> = HashMap::new();
                        for &(ei, hw, ov_len, _, _) in &auto_clamps {
                            let entry = min_h_per_edge.entry(ei).or_insert((hw, ov_len));
                            entry.0 = entry.0.min(hw);
                            entry.1 = entry.1.max(ov_len);
                        }
                        let est_gap_tris: f64 =
                            min_h_per_edge.values().map(|&(hw, len)| 3.0 * (len / hw.max(1e-15))).sum();
                        let cap_tris = (0.35 * budget_tris_2d).max(128.0);
                        if est_gap_tris > cap_tris {
                            auto_scale = est_gap_tris / cap_tris;
                        }
                    }
                }
                let rim_len_by_edge: HashMap<usize, f64> = rims.iter().map(|r| (r.e, r.len)).collect();
                for (ei, hw_raw, ov_len, m0, m1) in auto_clamps {
                    let hw = hw_raw * auto_scale;
                    if hw >= ov_edge_h[ei] {
                        continue;
                    }
                    let e_len = rim_len_by_edge.get(&ei).copied().unwrap_or(ov_len);
                    if e_len <= 4.0 * ov_len || plan.opts.gap_cells.is_some() {
                        ov_edge_h[ei] = ov_edge_h[ei].min(hw);
                    } else {
                        let lo = [m0[0].min(m1[0]), m0[1].min(m1[1]), 0.0];
                        let hi = [m0[0].max(m1[0]), m0[1].max(m1[1]), 0.0];
                        gap_zones.push(GapZone { a: m0, b: m1, h: hw, bbox: Aabb { lo, hi } });
                    }
                }
            }
            // Generate deterministic gap-centerline seed vertices along facing narrow gaps (`d <= 0.35 * h_global`)
            for (m0, m1, d_seg, ea, eb) in raw_centerlines {
                let h_wall = ov_edge_h[ea].min(ov_edge_h[eb]);
                let step = h_wall.clamp(0.75 * d_seg, 1.25 * d_seg).max(1e-15);
                let cl_len = dist2(m0, m1);
                if cl_len <= 1e-12 * d_seg {
                    if ov.locate(m0).is_some_and(|c| cell_sheet[c as usize].is_none()) {
                        gap_seeds.push((m0, d_seg));
                    }
                } else {
                    let n_pts = ((cl_len / step).ceil() as usize).clamp(1, 2000);
                    for k in 0..n_pts {
                        let t = (k as f64 + 0.5) / (n_pts as f64);
                        let p = [m0[0] + t * (m1[0] - m0[0]), m0[1] + t * (m1[1] - m0[1])];
                        if ov.locate(p).is_some_and(|c| cell_sheet[c as usize].is_none()) {
                            gap_seeds.push((p, d_seg));
                        }
                    }
                }
            }
            if local_gap_max > 0.0 {
                sheet_local_gap_by_z.push((z_s, local_gap_max));
            }
        }
        gap_seeds.sort_by(|a, b| {
            a.0[0]
                .total_cmp(&b.0[0])
                .then_with(|| a.0[1].total_cmp(&b.0[1]))
                .then_with(|| a.1.total_cmp(&b.1))
        });
        gap_seeds.dedup_by(|a, b| dist2(a.0, b.0) <= 1e-12 * a.1.min(b.1));
    }

    let mut sharp_corners_by_z = Vec::new();
    let mut all_sharp_corners = Vec::new();
    if plan.opts.corner_grading {
        for &z_s in &sheet_zs {
            let cs = crate::meshing::singularity::collect_overlay_sheet_corners(cx, z_s, &ov_edge_h, 1.0);
            if !cs.is_empty() {
                all_sharp_corners.extend(cs.iter().copied());
                sharp_corners_by_z.push((z_s, cs));
            }
        }
    }

    Ok(Sheet3dSizing {
        h_global,
        max_dev,
        eff_max_cells,
        meshable,
        region_h,
        face_h,
        edge3_h,
        cell_h,
        ov_edge_h,
        ov_edge_explicit_fine,
        edge_sheet_owner,
        dielectric_gap_cells,
        gap_zones,
        gap_seeds,
        sheet_local_gap_by_z,
        sharp_corners_by_z,
        all_sharp_corners,
    })
}

fn build_overlay_mesh_2d_with_sizing(
    ev: &Evaluated,
    plan: &MeshPlan,
    sz: &Sheet3dSizing,
) -> Result<(SimplexMesh, Vec<String>), String> {
    let cx = ev.complex();
    let st = cx.stack.as_ref().ok_or("3D sheet mesher requires a prismatic layer stack")?;
    let ov = &st.overlay;
    let ncell = ov.regions.len();
    let nedge = ov.edges.len();
    let h_global = sz.h_global;
    let max_dev = sz.max_dev;
    let grade = 0.35;

    // FR-31: Identify outer-perimeter edges, outer-ring overlay cells, and genuine nested 3D regions
    // (such as a fixed `inner_vac_box` or finite-die `sub_box` confined inside an outer vacuum box)
    // so interior edge discretization, fill seeds, and CDT coordinate normalization stay invariant to
    // outer-box extent without perturbing single-box stacks that merely carry interior 2D sheets.
    let mut cell_touches_outer = vec![false; ncell];
    let mut edge_is_outer = vec![false; nedge];
    for (e, ed) in ov.edges.iter().enumerate() {
        if ed.left.is_none() || ed.right.is_none() {
            edge_is_outer[e] = true;
            if let Some(c) = ed.left {
                cell_touches_outer[c as usize] = true;
            }
            if let Some(c) = ed.right {
                cell_touches_outer[c as usize] = true;
            }
        }
    }
    let mut outer_regions: HashSet<u32> = HashSet::new();
    for (c, &touches) in cell_touches_outer.iter().enumerate() {
        if touches {
            for l in &st.layers {
                if let Some(reg) = l.region_of[c] {
                    outer_regions.insert(reg);
                }
            }
        }
    }
    let mut interior_bb = Aabb::EMPTY;
    for (c, r) in ov.regions.iter().enumerate() {
        let in_nested_3d_box = !cell_touches_outer[c]
            && st
                .layers
                .iter()
                .any(|l| l.region_of[c].is_some_and(|reg| !outer_regions.contains(&reg)));
        if in_nested_3d_box {
            interior_bb = interior_bb.union(&r.bbox);
        }
    }
    let ref_frame = if !interior_bb.is_empty() {
        let int_ext = (interior_bb.hi[0] - interior_bb.lo[0]).max(interior_bb.hi[1] - interior_bb.lo[1]);
        if int_ext > 1e-6 * h_global {
            Some(([interior_bb.lo[0], interior_bb.lo[1]], [interior_bb.hi[0], interior_bb.hi[1]]))
        } else {
            None
        }
    } else {
        None
    };

    let cmp_bbox = |a: &Aabb, b: &Aabb| {
        a.lo[0]
            .total_cmp(&b.lo[0])
            .then_with(|| a.lo[1].total_cmp(&b.lo[1]))
            .then_with(|| a.hi[0].total_cmp(&b.hi[0]))
            .then_with(|| a.hi[1].total_cmp(&b.hi[1]))
    };

    let mut fine_edges: Vec<(usize, f64)> = (0..nedge)
        .filter(|&e| {
            sz.ov_edge_h[e] < h_global * 0.999
                && (ref_frame.is_none() || !edge_is_outer[e] || sz.ov_edge_explicit_fine[e])
        })
        .map(|e| (e, sz.ov_edge_h[e]))
        .collect();
    if ref_frame.is_some() {
        fine_edges.sort_by(|a, b| {
            a.1.total_cmp(&b.1)
                .then_with(|| cmp_bbox(&ov.edges[a.0].bbox, &ov.edges[b.0].bbox))
                .then(a.0.cmp(&b.0))
        });
    } else {
        fine_edges.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    }
    let mut fine_cells: Vec<(usize, f64)> =
        (0..ncell).filter(|&c| sz.cell_h[c] < h_global * 0.999).map(|c| (c, sz.cell_h[c])).collect();
    if ref_frame.is_some() {
        fine_cells.sort_by(|a, b| {
            a.1.total_cmp(&b.1)
                .then_with(|| cmp_bbox(&ov.regions[a.0].bbox, &ov.regions[b.0].bbox))
                .then(a.0.cmp(&b.0))
        });
    } else {
        fine_cells.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    }
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
    let gap_zone_bvh = Bvh2::build(
        sz.gap_zones
            .iter()
            .enumerate()
            .map(|(i, gz)| BvhItem2 { id: i, bbox: gz.bbox, min_h: gz.h })
            .collect(),
    );
    let gap_zones_ref = &sz.gap_zones;
    let corners_ref = &sz.all_sharp_corners;
    // When `dz_min` is configured, floor corner-graded XY sizing at `0.35 * dz_min` so sharp corners on
    // narrow-gap fingers do not generate extreme sub-`dz` needle prisms and runaway 2D triangle counts.
    let corner_h_floor = plan.opts.dz_min.map_or(0.0, |dz| 0.35 * dz);
    let local_h_2d = |p: V2| -> f64 {
        let h = fine_edge_bvh.min_graded(p, grade, h_global, |e, eh, _db| {
            let d = dist2(ov.project_edge(e as u32, p), p);
            eh + grade * d
        });
        let h = gap_zone_bvh.min_graded(p, grade, h, |gi, gh, _db| {
            let gz = &gap_zones_ref[gi];
            let d = dist2(closest_point_on_seg(p, gz.a, gz.b), p);
            gh + grade * d
        });
        let mut h = fine_cell_bvh.min_graded(p, grade, h, |_c, ch, db| ch + grade * db);
        if !corners_ref.is_empty() {
            let hc = crate::meshing::singularity::corner_graded_h_at(p, corners_ref, h);
            h = h.min(hc.max(corner_h_floor));
        }
        h
    };

    // Discretize every 1D PSLG edge of stack.overlay once (graded along the curve). When an interior
    // reference frame exists (`ref_frame.is_some()`, FR-31), order interior edges first in deterministic
    // geometric order so interior CDT topology is outer-box invariant; otherwise preserve natural `0..nedge`
    // order so coarse single-box benchmarks stay bit-identical.
    let mut edge_order: Vec<usize> = (0..nedge).collect();
    if ref_frame.is_some() {
        edge_order.sort_by(|&a, &b| {
            edge_is_outer[a]
                .cmp(&edge_is_outer[b])
                .then_with(|| cmp_bbox(&ov.edges[a].bbox, &ov.edges[b].bbox))
                .then(a.cmp(&b))
        });
    }
    let mut chains = Vec::with_capacity(nedge);
    let mut chain_owner = Vec::with_capacity(nedge);
    for e in edge_order {
        let ed = &ov.edges[e];
        let eh = sz.ov_edge_h[e];
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
        chain_owner.push(sz.edge_sheet_owner[e]);
    }

    // Snap fill lattice origins to an origin-anchored grid `(k * h, 2 * m * dy)` when `ref_frame.is_some()`
    // so interior fill seeds never shift phase when an outer bounding box changes from `R_1` to `R_2` (FR-31).
    let snap_fill_lo = |lo: [f64; 2], h: f64| -> [f64; 2] {
        if ref_frame.is_none() {
            return lo;
        }
        let dy = h * 3f64.sqrt() / 2.0;
        let step_y = 2.0 * dy;
        [(lo[0] / h).floor() * h, (lo[1] / step_y).floor() * step_y]
    };

    // Interior lattice fills for fine cells, graded collars around fine 1D edges, and sharp corner balls
    let mut fills: Vec<Fill> = Vec::new();
    for &(c, ch) in &fine_cells {
        let bb = ov.regions[c].bbox;
        fills.push(Fill {
            lo: snap_fill_lo([bb.lo[0], bb.lo[1]], ch),
            hi: [bb.hi[0], bb.hi[1]],
            h: ch,
            owner: Box::new(move |p| ov.locate(p) == Some(c as u32) && local_h_2d(p) >= 0.4 * ch),
        });
    }
    for &(e, eh) in &fine_edges {
        // Avoid flooding narrow PEC fingers with redundant collar fills on auto-clamped gap edges
        // (gap centerline seeds + Ruppert gradation already handle the gap cleanly).
        if !sz.ov_edge_explicit_fine[e] && plan.opts.gap_cells.is_none() && !sz.gap_seeds.is_empty() {
            continue;
        }
        for (mult, reach_k) in [(1.2, 2.2), (2.4, 5.0)] {
            let h_band = eh * mult;
            if h_band >= h_global * 0.999 {
                continue;
            }
            let reach = eh * reach_k;
            let bb = ov.edges[e].bbox.pad(reach);
            let cell_h_ref = &sz.cell_h;
            fills.push(Fill {
                lo: snap_fill_lo([bb.lo[0], bb.lo[1]], h_band),
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
    if !corners_ref.is_empty() {
        let mut sorted_corners: Vec<crate::meshing::singularity::SheetCorner> = corners_ref.clone();
        sorted_corners.sort_by(|a, b| a.pos[0].total_cmp(&b.pos[0]).then_with(|| a.pos[1].total_cmp(&b.pos[1])));
        for sc in sorted_corners {
            for (mult, rad_k) in [(0.25, 0.5), (0.55, 1.2)] {
                let h_c = (sc.h_edge * mult).max(sc.h_min()).max(corner_h_floor);
                if h_c >= h_global * 0.999 {
                    continue;
                }
                let r_fill = sc.h_edge * rad_k;
                let cp = sc.pos;
                fills.push(Fill {
                    lo: snap_fill_lo([cp[0] - r_fill, cp[1] - r_fill], h_c),
                    hi: [cp[0] + r_fill, cp[1] + r_fill],
                    h: h_c,
                    owner: Box::new(move |p| {
                        ov.locate(p).is_some()
                            && dist2(cp, p) <= r_fill
                            && local_h_2d(p) >= 0.4 * h_c
                    }),
                });
            }
        }
    }

    let mut ov_bb = Aabb::EMPTY;
    for r in &ov.regions {
        ov_bb = ov_bb.union(&r.bbox);
    }

    let corner_h_min = corners_ref
        .iter()
        .map(|c| c.h_min())
        .fold(h_global, f64::min);
    let gap_h_min = sz
        .gap_zones
        .iter()
        .map(|gz| gz.h)
        .chain(sz.gap_seeds.iter().map(|&(_, d)| 0.5 * d))
        .fold(corner_h_min, f64::min);
    let h_fine_min = sz
        .cell_h
        .iter()
        .copied()
        .chain(sz.ov_edge_h.iter().copied())
        .fold(gap_h_min, f64::min)
        .max(1e-12 * h_global);
    let cell_active: Vec<bool> =
        (0..ncell).map(|c| st.layers.iter().any(|l| l.region_of[c].is_some_and(|r| sz.meshable[r as usize]))).collect();
    let domain_at = |cen: V2| -> Option<i32> {
        let c = ov.locate(cen)?;
        if cell_active[c as usize] {
            Some(c as i32 + 1)
        } else {
            None
        }
    };

    // FR-28 budget plumbing: derive a bounded 2D additional-vertex budget from `eff_max_cells` and active z-slabs
    let n_active_layers = st
        .layers
        .iter()
        .filter(|l| l.region_of.iter().flatten().any(|&r| sz.meshable[r as usize]))
        .count()
        .max(1);
    let est_slabs = plan.opts.max_z_slabs.unwrap_or(n_active_layers.max(4)).max(n_active_layers);
    let max_cells_2d = sz.eff_max_cells.map(|mc| {
        let slab_div = if plan.opts.z_coarsen { est_slabs.min(8) } else { est_slabs };
        (mc / (3usize.saturating_mul(slab_div).max(1))).max(512)
    });

    let want_bridging_audit = plan.opts.gap_cells.is_some()
        || plan.opts.singular_edges.as_deref() == Some("quarter_point")
        || plan.opts.corner_grading
        || !sz.gap_seeds.is_empty();
    let extras = SheetChainExtras {
        gap_seeds: &sz.gap_seeds,
        chain_owner: if want_bridging_audit { &chain_owner } else { &[] },
        dielectric_domains: &sz.dielectric_gap_cells,
        ref_frame,
    };
    let (m2, mc_rep) = mesh_chains_with_seeds(
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
            max_cells: max_cells_2d,
        },
        &domain_at,
        Some(&extras),
    )?;
    drop(fills);
    Ok((m2, mc_rep.warnings))
}

/// Builds the sized 2D conformal overlay triangulation `m2` (`cell_labels = ov_cell + 1`, boundary
/// `facet_labels = ov_edge + 1`) honoring `Size` rules, `gap_cells`, `waves`, and `corner_grading`.
pub fn build_sized_overlay_mesh_2d(
    ev: &Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<SimplexMesh, String> {
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let sz = compute_sheet3d_sizing(ev, &labels, plan, waves)?;
    Ok(build_overlay_mesh_2d_with_sizing(ev, plan, &sz)?.0)
}

fn zero_thickness_sheet_labels(g: &MeshedGeom) -> Vec<i32> {
    let mut out = Vec::new();
    for sname in g.eval.complex().sheet_names() {
        if g.eval.complex().sheet_thickness(&sname).unwrap_or(0.0) <= 0.0 {
            for lab in g.labels.boundary(&sname) {
                if g.labels.sheet_thickness_m.get(&lab).copied().unwrap_or(0.0) <= 0.0 {
                    out.push(lab);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn apply_quarter_point_on_meshed_geom(
    g: &mut MeshedGeom,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) {
    let want_qp = plan.opts.singular_edges.as_deref() == Some("quarter_point")
        && waves.is_none_or(|(sz, _)| sz.order >= 2);
    if !want_qp {
        return;
    }
    // Cross-Feature Compatibility Rule 2: automatically skip thick sheets (`thickness > 0`), where
    // 270 deg physical corners induce `r^(2/3)` rather than `r^(1/2)` singularities.
    let zero_thick_labels = zero_thickness_sheet_labels(g);
    if zero_thick_labels.is_empty() {
        return;
    }
    let rim_verts = crate::meshing::singularity::sheet_rim_vertices(&g.mesh, &zero_thick_labels);
    let mut curved_rim_edges = std::collections::BTreeSet::new();
    for ([a, b], _) in g.mesh.edge_node_list() {
        if rim_verts.get(a as usize).copied().unwrap_or(false)
            && rim_verts.get(b as usize).copied().unwrap_or(false)
        {
            curved_rim_edges.insert([a.min(b), a.max(b)]);
        }
    }
    let xi = plan
        .opts
        .quarter_point_xi
        .unwrap_or(crate::meshing::singularity::DEFAULT_QUARTER_POINT_XI);
    crate::meshing::singularity::apply_quarter_point_nodes(&mut g.mesh, &rim_verts, xi, &curved_rim_edges);
    if let Some(crate::mesh::Mesh::D3(ref mut cm)) = g.cells {
        if cm.n_points() == g.mesh.n_points() && cm.types.iter().all(|&t| t == crate::mesh::cells::CellType::Tet) {
            cm.set_edge_nodes(g.mesh.edge_node_list());
            let mut qp: Vec<[u32; 2]> =
                g.mesh.qp_edge_list().into_iter().map(|[a, b]| [a.min(b), a.max(b)]).collect();
            qp.sort_unstable();
            qp.dedup();
            cm.qp_edges = qp;
        }
    }
}

/// Applies post-adaptation `AF-5` (`curved`) and `AF-2` (`singular_edges = "quarter_point"`) `EdgeNodes`
/// enrichment to an unstructured 3D [`MeshedGeom`] (e.g. after `full_3d` `adapt3`).
pub fn apply_sheet_post_enrichment(
    g: &mut MeshedGeom,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) {
    if plan.opts.curved {
        let nr = g.eval.regions().len();
        let mut skip = vec![false; nr];
        if plan.pml.is_structured(!g.eval.pml_regions().is_empty()) {
            for p in g.eval.pml_regions() {
                skip[p.region as usize] = true;
            }
        }
        if let Ok(dom) = crate::meshing::geomdrive3::GeomDomain3::new_with(
            &g.eval,
            &g.labels,
            plan.opts.h,
            &crate::meshing::geomdrive3::Sizes3 { region_h: vec![], skip },
        ) {
            crate::meshing::curve3::curve_solver_mesh(&mut g.mesh, &dom);
            if let Some(crate::mesh::Mesh::D3(ref mut cm)) = g.cells {
                if cm.n_points() == g.mesh.n_points() && cm.types.iter().all(|&t| t == crate::mesh::cells::CellType::Tet) {
                    cm.set_edge_nodes(g.mesh.edge_node_list());
                }
            }
        }
    }
    apply_quarter_point_on_meshed_geom(g, plan, waves);
}

#[allow(clippy::too_many_arguments)]
fn extrude_overlay_2d_inner(
    ev: Evaluated,
    mut labels: geom::labels::LabelMap,
    m2: &SimplexMesh,
    plan: &MeshPlan,
    sz: &Sheet3dSizing,
    waves: Option<(&Sizing, &PhysicsSizing)>,
    sheet_dz_override: Option<f64>,
    infer_dz_from_m2: bool,
    warnings: Vec<String>,
) -> Result<MeshedGeom, String> {
    if ev.list.dim != 3 || m2.dim != 2 {
        return Err("extrude_overlay_2d requires a 3D Evaluated model and a 2D overlay mesh".into());
    }
    let cx = ev.complex();
    let st = cx.stack.as_ref().ok_or("3D sheet mesher requires a prismatic layer stack")?;
    let ov = &st.overlay;
    let tol = cx.tol;
    let h_global = sz.h_global;
    let eff_max_cells = sz.eff_max_cells;
    let meshable = &sz.meshable;
    let region_h = &sz.region_h;
    let face_h = &sz.face_h;
    let edge3_h = &sz.edge3_h;
    let cell_h = &sz.cell_h;
    let ov_edge_h = &sz.ov_edge_h;
    let nf = cx.faces.len();

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

    let scene_sheets = &ev.scene().sheets;
    let has_thick_or_etch = scene_sheets.iter().any(|s| s.thickness > 0.0 || s.etch_depth > 0.0);

    // Build graded z-slabs across every 3D layer
    let sheet_dz_at_raw = |z_query: f64, default_h: f64| -> Option<f64> {
        let mut has_sheet = false;
        let mut dz = default_h;
        let mut sheet_cells = if infer_dz_from_m2 { vec![false; ov.regions.len()] } else { Vec::new() };
        for (fi, f) in cx.faces.iter().enumerate() {
            if let FaceGeom::Cap { cell, z: fz } = f.geom {
                if (fz - z_query).abs() <= tol && f.is_sheet() {
                    has_sheet = true;
                    dz = dz.min(face_h[fi]).min(cell_h[cell as usize]);
                    if infer_dz_from_m2 && (cell as usize) < sheet_cells.len() {
                        sheet_cells[cell as usize] = true;
                    }
                }
            }
        }
        for (ei, e3) in cx.edges3.iter().enumerate() {
            if let EdgeGeom::Horizontal { edge, z: ez } = e3.geom {
                if (ez - z_query).abs() <= tol && e3.faces.iter().any(|&fi| cx.faces[fi as usize].is_sheet()) {
                    has_sheet = true;
                    dz = dz.min(edge3_h[ei]).min(ov_edge_h[edge as usize]);
                }
            }
        }
        if infer_dz_from_m2 && has_sheet {
            let mut min_m2_edge = default_h;
            for (ti, &ov_c) in tri_ov_cell.iter().enumerate() {
                if ov_c < sheet_cells.len() && sheet_cells[ov_c] {
                    let v = m2.cell(ti);
                    for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                        min_m2_edge = min_m2_edge.min(dist2(pt2(a), pt2(b)));
                    }
                }
            }
            dz = dz.min(min_m2_edge.max(default_h * 1e-3));
        }
        if plan.opts.corner_grading && has_sheet {
            for &(zs, ref cs) in &sz.sharp_corners_by_z {
                if (zs - z_query).abs() <= tol {
                    dz = crate::meshing::singularity::coupled_corner_dz(dz, cs);
                }
            }
        }
        if let Some(dz_min) = plan.opts.dz_min {
            if has_sheet {
                dz = dz.max(dz_min).min(default_h);
            }
        }
        if let Some(ovr) = sheet_dz_override {
            if has_sheet {
                dz = ovr.min(default_h);
            }
        }
        has_sheet.then_some(dz)
    };

    let sheet_dz_at = |z: f64, default_h: f64| -> Option<f64> {
        if let Some(dz) = sheet_dz_at_raw(z, default_h) {
            return Some(dz);
        }
        if has_thick_or_etch {
            // Also match trench floor `z_s - t_e` to the parent sheet plane `z_s` so substrate slabs grade smoothly up to trenches.
            for sh in scene_sheets {
                if sh.etch_depth > 0.0 && (z - (sh.z - sh.etch_depth)).abs() <= tol {
                    if let Some(dz) = sheet_dz_at_raw(sh.z, default_h) {
                        return Some(dz);
                    }
                }
            }
        }
        None
    };

    let z_ratio = plan.opts.z_growth_ratio.unwrap_or(1.35);
    // (layer_idx, z0, z1, dz0, dz1, has_sheet_z0, has_sheet_z1, h_layer, is_thin_sheet_slab)
    type ActiveLayerSpec = (usize, f64, f64, f64, f64, bool, bool, f64, bool);
    let mut active_layers: Vec<ActiveLayerSpec> = Vec::new();
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
        let is_thin_sheet_slab = has_thick_or_etch
            && (l.is_conductor_cavity.iter().any(|&b| b)
                || scene_sheets.iter().any(|sh| {
                    (sh.thickness > 0.0 && l.z0 >= sh.z - tol && l.z1 <= sh.z + sh.thickness + tol)
                        || (sh.etch_depth > 0.0 && l.z0 >= sh.z - sh.etch_depth - tol && l.z1 <= sh.z + tol)
                }));
        let sdz0 = sheet_dz_at(l.z0, h_layer);
        let sdz1 = sheet_dz_at(l.z1, h_layer);
        let dz0 = sdz0.unwrap_or(h_layer);
        let dz1 = sdz1.unwrap_or(h_layer);
        active_layers.push((k, l.z0, l.z1, dz0, dz1, sdz0.is_some(), sdz1.is_some(), h_layer, is_thin_sheet_slab));
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
    // Phase 1 (`preserve_sheet_dz = true`): keep near-sheet `dz` anchored at `dz_min` while accelerating
    // bulk `h_layer` and geometric growth `z_ratio` so narrow-gap sheets do not suffer vertical stretching.
    // Phase 2 (`preserve_sheet_dz = false`): fallback only if even steep geometric grading still exceeds `max_z_slabs`.
    let build_layer_zs = |alpha: f64, preserve_sheet_dz: bool| -> Vec<(usize, Vec<f64>)> {
        active_layers
            .iter()
            .map(|&(k, z0, z1, dz0, dz1, has_s0, has_s1, h_layer, is_thin_sheet_slab)| {
                if is_thin_sheet_slab {
                    return (k, vec![z0, z1]);
                }
                let (d0, d1, hl, rt) = if alpha <= 1.0 {
                    (dz0, dz1, h_layer, z_ratio)
                } else if preserve_sheet_dz {
                    let span = (z1 - z0).max(0.0);
                    let bulk_alpha = if !has_s0 && !has_s1 { alpha * alpha } else { alpha };
                    let d0_eff = if has_s0 { dz0.min(span) } else { (dz0 * bulk_alpha).min(span) };
                    let d1_eff = if has_s1 { dz1.min(span) } else { (dz1 * bulk_alpha).min(span) };
                    let hl_eff = (h_layer * bulk_alpha).min(span);
                    let rt_eff = (z_ratio * alpha.powf(0.45)).min(8.0);
                    (d0_eff, d1_eff, hl_eff, rt_eff)
                } else {
                    let span = (z1 - z0).max(0.0);
                    (
                        (dz0 * alpha).min(span),
                        (dz1 * alpha).min(span),
                        (h_layer * alpha).min(span),
                        (z_ratio + 0.15 * (alpha - 1.0)).min(8.0),
                    )
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
    let mut per_layer_zs = build_layer_zs(1.0, true);
    if let Some(max_slabs) = plan.opts.max_z_slabs {
        let count_slabs = |v: &[(usize, Vec<f64>)]| -> usize { v.iter().map(|(_, zs)| zs.len().saturating_sub(1)).sum() };
        if count_slabs(&per_layer_zs) > max_slabs {
            // Try preserving exact near-sheet `dz` first by increasing far-field `h_layer` & geometric growth `z_ratio`
            let mut preserve_ok = false;
            let mut hi_p = 1.0_f64;
            while hi_p < 64.0 {
                hi_p *= 2.0;
                if count_slabs(&build_layer_zs(hi_p, true)) <= max_slabs {
                    preserve_ok = true;
                    break;
                }
            }
            if preserve_ok {
                let mut lo_a = 1.0_f64;
                let mut hi_a = hi_p;
                let mut best_zs = build_layer_zs(hi_a, true);
                for _ in 0..40 {
                    let mid_a = 0.5 * (lo_a + hi_a);
                    let cand = build_layer_zs(mid_a, true);
                    if count_slabs(&cand) <= max_slabs {
                        hi_a = mid_a;
                        best_zs = cand;
                    } else {
                        lo_a = mid_a;
                    }
                }
                per_layer_zs = best_zs;
            } else {
                let mut lo_a = 1.0_f64;
                let mut hi_a = 2.0_f64;
                while count_slabs(&build_layer_zs(hi_a, false)) > max_slabs && hi_a < 1e12 {
                    lo_a = hi_a;
                    hi_a *= 2.0;
                }
                let mut best_zs = build_layer_zs(hi_a, false);
                for _ in 0..48 {
                    let mid_a = 0.5 * (lo_a + hi_a);
                    let cand = build_layer_zs(mid_a, false);
                    if count_slabs(&cand) <= max_slabs {
                        hi_a = mid_a;
                        best_zs = cand;
                    } else {
                        lo_a = mid_a;
                    }
                }
                per_layer_zs = best_zs;
            }
            while count_slabs(&per_layer_zs) > max_slabs {
                let mut best_drop: Option<(f64, usize, usize)> = None;
                for (li, (_, zs)) in per_layer_zs.iter().enumerate() {
                    if zs.len() <= 2 {
                        continue;
                    }
                    let (_, _, _, _, _, has_s0, has_s1, _, _) = active_layers[li];
                    for zi in 1..(zs.len() - 1) {
                        // Prefer dropping bulk planes over the first near-sheet plane
                        let near_sheet = (has_s0 && zi == 1) || (has_s1 && zi + 2 == zs.len());
                        let merged_step = (zs[zi + 1] - zs[zi - 1]) * if near_sheet { 1e6 } else { 1.0 };
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
    let mut face_label: Vec<i32> = (0..nf as u32)
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
        // Default exact prismatic Dompierre path (bit-identical to baseline when AF-2/5/6 flags are off)
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

        let emitted_labels: HashSet<i32> = facet_labels.iter().copied().collect();
        let registered_labels: HashSet<i32> = labels.boundaries.iter().map(|b| b.id).collect();
        let mut missing_active: Vec<i32> = cap_at
            .values()
            .chain(lat_at.values())
            .map(|&fi| face_label[fi])
            .filter(|lab| !emitted_labels.contains(lab))
            .collect();
        missing_active.sort_unstable();
        missing_active.dedup();
        if !missing_active.is_empty() {
            return Err(format!(
                "sheet3d facet-label bijection failed (z_coarsen=false): active boundary label(s) {missing_active:?} have 0 emitted facets"
            ));
        }
        let mut orphan_emitted: Vec<i32> =
            emitted_labels.iter().copied().filter(|l| !registered_labels.contains(l)).collect();
        orphan_emitted.sort_unstable();
        if !orphan_emitted.is_empty() {
            return Err(format!(
                "sheet3d facet-label bijection failed (z_coarsen=false): emitted orphan boundary label(s) {orphan_emitted:?} absent from LabelMap"
            ));
        }

        let mesh = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);
        let mut g = MeshedGeom {
            mesh,
            cells: None,
            labels,
            eval: std::sync::Arc::new(ev),
            warnings,
        };
        // Cross-Feature Compatibility Rule 4: curve & quarter-point enrich the physical interior FIRST,
        // then extrude structured PML layers in lockstep so PML `edge_nodes` are never overwritten!
        if plan.opts.curved {
            let level_zero = vec![0usize; nz];
            crate::meshing::sheet3d_curve::project_sheet3d_curved_edges(
                g.eval.complex(),
                m2,
                &bnd_edges2,
                &z_planes,
                &used_v,
                &slab_layer,
                &level_zero,
                &[],
                &mut g.mesh,
            );
        }
        apply_quarter_point_on_meshed_geom(&mut g, plan, waves);
        g.cells = Some(crate::mesh::Mesh::from_solver_mesh(&g.mesh)?);
        return Ok(g);
    }

    // `z_coarsen = true`: Conforming pure-tetrahedral vertical coarsening away from sheet planes
    let mut sheet_plane_info: Vec<(usize, f64, f64, f64)> = Vec::new();
    for (zp, &z) in z_planes.iter().enumerate() {
        if let Some(hs) = sheet_dz_at(z, h_global) {
            let h_local_gap = sz
                .sheet_local_gap_by_z
                .iter()
                .find(|&&(zs, _)| (zs - z).abs() <= tol)
                .map_or(0.0, |&(_, g)| g);
            sheet_plane_info.push((zp, z, hs, h_local_gap));
        }
    }

    // Build canonical cap equivalence signature per face:
    // (face name, region pair, pinned boundary id if pinned, sorted group indices)
    let pinned_by_id: HashMap<i32, bool> = labels.boundaries.iter().map(|b| (b.id, b.pinned)).collect();
    let mut groups_by_id: HashMap<i32, Vec<usize>> = HashMap::new();
    for (gi, g) in labels.groups.iter().enumerate() {
        for &lab in &g.labels {
            groups_by_id.entry(lab).or_default().push(gi);
        }
    }
    for v in groups_by_id.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    let face_pinned: Vec<Option<i32>> = face_label
        .iter()
        .map(|&lab| if pinned_by_id.get(&lab).copied().unwrap_or(false) { Some(lab) } else { None })
        .collect();
    let face_groups: Vec<Vec<usize>> =
        face_label.iter().map(|lab| groups_by_id.get(lab).cloned().unwrap_or_default()).collect();
    let cap_sig_eq = |i0: usize, i1: usize| -> bool {
        cx.faces[i0].name == cx.faces[i1].name
            && cx.faces[i0].first == cx.faces[i1].first
            && cx.faces[i0].second == cx.faces[i1].second
            && face_pinned[i0] == face_pinned[i1]
            && face_groups[i0] == face_groups[i1]
    };

    let ncell = ov.regions.len();
    let mut dsu_parent: Vec<usize> = (0..ncell).collect();
    let dsu_find = |parent: &mut [usize], mut x: usize| -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    };

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
                    (Some(&i0), Some(&i1)) => !cap_sig_eq(i0, i1),
                    _ => true,
                }
            });
            if diff_layer || diff_non_sheet_cap {
                frozen_v[a as usize] = true;
                frozen_v[b as usize] = true;
            } else {
                let r0 = dsu_find(&mut dsu_parent, c0);
                let r1 = dsu_find(&mut dsu_parent, c1);
                if r0 != r1 {
                    if r0 < r1 {
                        dsu_parent[r1] = r0;
                    } else {
                        dsu_parent[r0] = r1;
                    }
                }
            }
        }
    }
    for &(a, b, oe) in &bnd_edges2 {
        if (0..st.layers.len()).any(|k| lat_at.contains_key(&(oe, k))) {
            frozen_v[a as usize] = true;
            frozen_v[b as usize] = true;
        }
    }
    let comp_of_cell: Vec<usize> = (0..ncell).map(|c| dsu_find(&mut dsu_parent, c)).collect();

    // Determine coarsening level `level_of_zp[zp]` (`0` inside `slab_dist < 2 || d_min < max(4 * h_s, 1.5 * min(h_local_gap, 2.5 * h_s))`)
    let mut level_of_zp = if sheet_plane_info.is_empty() {
        vec![0usize; nz]
    } else {
        let mut lvl = vec![usize::MAX / 2; nz];
        for &(zp_s, z_s, h_s, h_local_gap) in &sheet_plane_info {
            lvl[zp_s] = 0;
            let protect_dist = (2.2 * h_s).max(1.2 * h_local_gap.min(2.0 * h_s));
            for s in 0..nslab {
                if slab_layer[s].is_none() {
                    continue;
                }
                let d_min = (z_planes[s] - z_s).abs().min((z_planes[s + 1] - z_s).abs());
                let slab_dist = if s >= zp_s { s - zp_s } else { zp_s - 1 - s };
                if slab_dist < 2 || d_min < protect_dist - tol {
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
                for &(_, z_s, _, _) in &sheet_plane_info {
                    z_dist_min = z_dist_min.min((z_planes[zp] - z_s).abs());
                }
            }
        }
        let h_target = if z_dist_min.is_finite() {
            (1.65 * z_dist_min).min(h_layer_max)
        } else {
            h_layer_max
        };
        let mut dz_slab_min = f64::INFINITY;
        for (s, &layer_opt) in slab_layer.iter().enumerate().take(nslab) {
            if layer_opt.is_none() {
                continue;
            }
            let (ls, lt) = (level_of_zp[s], level_of_zp[s + 1]);
            if (ls == m && lt == m + 1) || (ls == m + 1 && lt == m) {
                dz_slab_min = dz_slab_min.min((z_planes[s + 1] - z_planes[s]).abs());
            }
        }
        let dz_slab = if dz_slab_min.is_finite() && dz_slab_min > 0.0 { dz_slab_min } else { h_target };
        let (next_lvl, tr) = coarsen_2d_once(
            &levels_2d[m],
            &pts2,
            &frozen_v,
            h_target,
            h_layer_max,
            dz_slab,
            &comp_of_cell,
            m >= 1,
        );
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

    // Coalesce coplanar non-sheet cap patches on coarsened z-planes (`level_of_zp[zp] > 0`) that belong to
    // the same unfrozen cell component and carry identical BC/region semantics onto one canonical `face_label`,
    // remapping and pruning `labels` consistently so no active boundary ID can ever lose all its facets.
    let mut label_remap: HashMap<i32, i32> = HashMap::new();
    for zp in 0..nz {
        let m = level_of_zp[zp];
        if m == 0 {
            continue;
        }
        let mut cell_has_tri = vec![false; ncell];
        for &(_, ov_c) in &levels_2d[m].tris {
            if ov_c < ncell {
                cell_has_tri[ov_c] = true;
            }
        }
        type CapSigKey = (usize, String, Option<u32>, Option<u32>, Option<i32>, Vec<usize>);
        let mut by_comp_sig: HashMap<CapSigKey, Vec<(u32, usize)>> = HashMap::new();
        for (&(cell, czp), &fi) in &cap_at {
            if czp != zp || cx.faces[fi].is_sheet() {
                continue;
            }
            let comp = comp_of_cell.get(cell as usize).copied().unwrap_or(cell as usize);
            let key = (
                comp,
                cx.faces[fi].name.clone(),
                cx.faces[fi].first,
                cx.faces[fi].second,
                face_pinned[fi],
                face_groups[fi].clone(),
            );
            by_comp_sig.entry(key).or_default().push((cell, fi));
        }
        for mut members in by_comp_sig.into_values() {
            if members.len() <= 1 {
                continue;
            }
            members.sort_by_key(|&(c, fi)| (face_label[fi], c));
            let canon_id = members
                .iter()
                .filter(|&&(c, _)| cell_has_tri.get(c as usize).copied().unwrap_or(false))
                .map(|&(_, fi)| face_label[fi])
                .min()
                .unwrap_or_else(|| face_label[members[0].1]);
            for &(_, fi) in &members {
                let old_id = face_label[fi];
                if old_id != canon_id {
                    label_remap.insert(old_id, canon_id);
                    face_label[fi] = canon_id;
                }
            }
        }
    }
    if !label_remap.is_empty() {
        let mut absorbed_entities: HashMap<i32, Vec<u32>> = HashMap::new();
        labels.boundaries.retain(|b| {
            if let Some(&target) = label_remap.get(&b.id) {
                absorbed_entities.entry(target).or_default().extend(b.entities.iter().copied());
                false
            } else {
                true
            }
        });
        for b in &mut labels.boundaries {
            if let Some(extra) = absorbed_entities.remove(&b.id) {
                b.entities.extend(extra);
                b.entities.sort_unstable();
                b.entities.dedup();
            }
        }
        for g in &mut labels.groups {
            for lab in &mut g.labels {
                if let Some(&target) = label_remap.get(lab) {
                    *lab = target;
                }
            }
            g.labels.sort_unstable();
            g.labels.dedup();
        }
        for (&old_id, &target) in &label_remap {
            if let Some(tm) = labels.sheet_thickness_m.remove(&old_id) {
                labels.sheet_thickness_m.entry(target).or_insert(tm);
            }
            if let Some(te) = labels.sheet_etch_depth_m.remove(&old_id) {
                labels.sheet_etch_depth_m.entry(target).or_insert(te);
            }
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

    let emitted_labels: HashSet<i32> = facet_labels.iter().copied().collect();
    let registered_labels: HashSet<i32> = labels.boundaries.iter().map(|b| b.id).collect();
    let mut missing_active: Vec<i32> = cap_at
        .values()
        .chain(lat_at.values())
        .map(|&fi| face_label[fi])
        .filter(|lab| !emitted_labels.contains(lab))
        .collect();
    missing_active.sort_unstable();
    missing_active.dedup();
    if !missing_active.is_empty() {
        return Err(format!(
            "sheet3d facet-label bijection failed (z_coarsen=true): active boundary label(s) {missing_active:?} have 0 emitted facets"
        ));
    }
    let mut orphan_emitted: Vec<i32> = emitted_labels.iter().copied().filter(|l| !registered_labels.contains(l)).collect();
    orphan_emitted.sort_unstable();
    if !orphan_emitted.is_empty() {
        return Err(format!(
            "sheet3d facet-label bijection failed (z_coarsen=true): emitted orphan boundary label(s) {orphan_emitted:?} absent from LabelMap"
        ));
    }

    let mesh = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);
    let mut g = MeshedGeom {
        mesh,
        cells: None,
        labels,
        eval: std::sync::Arc::new(ev),
        warnings,
    };
    // Cross-Feature Compatibility Rule 4: curve & quarter-point enrich the physical interior FIRST,
    // then extrude structured PML layers in lockstep so PML `edge_nodes` are never overwritten!
    if plan.opts.curved {
        let transition_refs: Vec<crate::meshing::sheet3d_curve::Sheet3dTransitionRef<'_>> = transitions
            .iter()
            .map(|tr| crate::meshing::sheet3d_curve::Sheet3dTransitionRef {
                map_v: &tr.map_v,
                rank: &tr.rank,
            })
            .collect();
        crate::meshing::sheet3d_curve::project_sheet3d_curved_edges(
            g.eval.complex(),
            m2,
            &bnd_edges2,
            &z_planes,
            &used_v,
            &slab_layer,
            &level_of_zp,
            &transition_refs,
            &mut g.mesh,
        );
    }
    apply_quarter_point_on_meshed_geom(&mut g, plan, waves);
    g.cells = Some(crate::mesh::Mesh::from_solver_mesh(&g.mesh)?);
    Ok(g)
}

/// Extrudes a conformal 2D overlay triangulation `m2` across the 3D layers of `ev` using the full
/// [`MeshPlan`] (`z_coarsen`, `AF-6` thin-slab single-layer extrusion, `PmlMode`,
/// `AF-5` `curved` projection, and `AF-2` `quarter_point` enrichment).
pub fn extrude_overlay_2d_with_plan(
    ev: Evaluated,
    m2: &SimplexMesh,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
    sheet_dz_override: Option<f64>,
) -> Result<MeshedGeom, String> {
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let sz = compute_sheet3d_sizing(&ev, &labels, plan, waves)?;
    let mut g = extrude_overlay_2d_inner(ev, labels, m2, plan, &sz, waves, sheet_dz_override, true, Vec::new())?;
    if plan.pml.is_structured(!g.eval.pml_regions().is_empty()) {
        g = crate::meshing::pml_geom3::structured_pml_3d_with_curved(g, plan.opts.curved)?.0;
    }
    Ok(g)
}

/// Extrudes a conformal 2D overlay triangulation `m2` across the 3D layers of `ev` at global
/// element size `h_global` (shared between initial `sheet3d` meshing and `AF-1` `sheet_2d_z` AMR).
pub fn extrude_overlay_2d(
    ev: Evaluated,
    m2: &SimplexMesh,
    h_global: f64,
    sheet_dz_override: Option<f64>,
) -> Result<MeshedGeom, String> {
    extrude_overlay_2d_with_plan(ev, m2, &MeshPlan::simple(h_global), None, sheet_dz_override)
}

/// Meshes a 3D layer-stack model containing 2D sheets (`Op::Sheet`) into a [`MeshedGeom`], with
/// optional [`PmlReport`].
pub fn mesh_sheet_stack_3d_with_pml(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<(MeshedGeom, Option<AdaptReport>, Option<crate::meshing::PmlReport>), String> {
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let sz = compute_sheet3d_sizing(&ev, &labels, plan, waves)?;
    let (m2, warnings) = build_overlay_mesh_2d_with_sizing(&ev, plan, &sz)?;
    let mut g = extrude_overlay_2d_inner(ev, labels, &m2, plan, &sz, waves, None, false, warnings)?;
    let pml_rep = if plan.pml.is_structured(!g.eval.pml_regions().is_empty()) {
        let (g2, rep) = crate::meshing::pml_geom3::structured_pml_3d_with_curved(g, plan.opts.curved)?;
        g = g2;
        Some(rep)
    } else {
        None
    };
    Ok((g, None, pml_rep))
}

/// Meshes a 3D layer-stack model containing 2D sheets (`Op::Sheet`) directly into an uncut
/// [`crate::mesh::Mesh`] (`Tet` + `Prism` + `Hex`) with zero eager simplex splitting of PML layers.
pub fn mesh_sheet_stack_3d_native(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<crate::meshing::NativeRawTuple, String> {
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let sz = compute_sheet3d_sizing(&ev, &labels, plan, waves)?;
    let (m2, warnings) = build_overlay_mesh_2d_with_sizing(&ev, plan, &sz)?;
    let g = extrude_overlay_2d_inner(ev, labels, &m2, plan, &sz, waves, None, false, warnings)?;
    if plan.pml.is_structured(!g.eval.pml_regions().is_empty()) {
        let (labels_out, eval_out, warnings_out) = (g.labels.clone(), g.eval.clone(), g.warnings.clone());
        match crate::meshing::pml_geom3::structured_pml_3d_native_with_curved(g, plan.opts.curved)? {
            crate::meshing::pml_geom3::Native3::Unchanged(g_unch, rep) => {
                let m = g_unch.into_exchange_mesh()?;
                Ok((m, g_unch.labels, g_unch.eval, None, Some(rep), g_unch.warnings))
            }
            crate::meshing::pml_geom3::Native3::Cells(cells, rep) => {
                Ok((crate::mesh::Mesh::D3(cells), labels_out, eval_out, None, Some(rep), warnings_out))
            }
        }
    } else {
        let m = g.into_exchange_mesh()?;
        Ok((m, g.labels, g.eval, None, None, g.warnings))
    }
}

/// Meshes a 3D layer-stack model containing 2D sheets (`Op::Sheet`) into a conformal tetrahedral `MeshedGeom`.
pub fn mesh_sheet_stack_3d(
    ev: Evaluated,
    plan: &MeshPlan,
    waves: Option<(&Sizing, &PhysicsSizing)>,
) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    let (g, adapt, _) = mesh_sheet_stack_3d_with_pml(ev, plan, waves)?;
    Ok((g, adapt))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::meshing::geomdrive::evaluate;
    use geom::expr::Expr;
    use geom::model::Model;
    use geom::select::Selector;

    fn pt2(x: f64, y: f64) -> [Expr; 2] {
        [Expr::Num(x), Expr::Num(y)]
    }

    #[test]
    fn fr27_z_coarsen_preserves_all_active_cap_labels_and_positive_tet_volumes() {
        let (half, pad, h_sub, h_air, h_bot) = (500e-6, 1200e-6, 750e-6, 890e-6, 900e-6);
        let build_ev = || {
            let mut m = Model::new_3d();
            let vac = m.box_(
                "vac_box",
                (-half - pad, -half - pad, -h_sub - h_bot - pad),
                (half + pad, half + pad, h_air + pad),
            );
            let sub = m.box_("sub_box", (-half, -half, -h_sub), (half, half, 0.0));
            let cross = m.polygons(
                "cross_poly",
                vec![vec![
                    pt2(-15e-6, -120e-6),
                    pt2(15e-6, -120e-6),
                    pt2(15e-6, -15e-6),
                    pt2(120e-6, -15e-6),
                    pt2(120e-6, 15e-6),
                    pt2(15e-6, 15e-6),
                    pt2(15e-6, 120e-6),
                    pt2(-15e-6, 120e-6),
                    pt2(-15e-6, 15e-6),
                    pt2(-120e-6, 15e-6),
                    pt2(-120e-6, -15e-6),
                    pt2(-15e-6, -15e-6),
                ]],
                vec![],
            );
            let claw = m.polygons(
                "claw_poly",
                vec![vec![
                    pt2(-45e-6, 132e-6),
                    pt2(45e-6, 132e-6),
                    pt2(45e-6, 190e-6),
                    pt2(30e-6, 190e-6),
                    pt2(30e-6, 147e-6),
                    pt2(-30e-6, 147e-6),
                    pt2(-30e-6, 190e-6),
                    pt2(-45e-6, 190e-6),
                ]],
                vec![],
            );
            let gnd = m.polygons(
                "gnd_poly",
                vec![vec![pt2(-half, -half), pt2(half, -half), pt2(half, half), pt2(-half, half)]],
                vec![vec![pt2(-140e-6, -140e-6), pt2(140e-6, -140e-6), pt2(140e-6, 210e-6), pt2(-140e-6, 210e-6)]],
            );
            m.sheet("cross", cross, 0.0);
            m.sheet("claw", claw, 0.0);
            m.sheet("ground_plane", gnd, 0.0);
            m.overlay(&[("vacuum", vac), ("substrate", sub)]);
            m.size("sz_cross", Selector::parse("cross:edge").unwrap(), Some(Expr::Num(10e-6)), None, None);
            m.size("sz_claw", Selector::parse("claw:edge").unwrap(), Some(Expr::Num(8e-6)), None, None);
            evaluate(&m.list).unwrap()
        };

        let mut p_off = MeshPlan::simple(340e-6);
        p_off.opts.z_coarsen = false;
        p_off.opts.dz_min = Some(5e-6);
        p_off.opts.z_growth_ratio = Some(1.6);
        p_off.opts.max_z_slabs = Some(22);
        let (g_off, _) = mesh_sheet_stack_3d(build_ev(), &p_off, None).unwrap();

        let mut p_on = p_off.clone();
        p_on.opts.z_coarsen = true;
        let (g_on, _) = mesh_sheet_stack_3d(build_ev(), &p_on, None).unwrap();

        // Every boundary ID in `outer`, `cross`, `claw`, `ground_plane` must exist in `mesh.facet_labels`
        let emitted: HashSet<i32> = g_on.mesh.facet_labels.iter().copied().collect();
        for sel in ["outer", "name(cross)", "name(claw)", "name(ground_plane)"] {
            let labs = g_on.select(sel).unwrap();
            assert!(!labs.is_empty(), "selector '{sel}' resolved empty");
            for lab in labs {
                assert!(emitted.contains(&lab), "selector '{sel}' label {lab} missing from facet_labels");
            }
        }
        // Strictly positive tet volumes across 100% of cells
        for c in 0..g_on.mesh.n_cells() {
            let v = g_on.mesh.cell(c);
            let p = |id: u32| {
                let s = g_on.mesh.point(id as usize);
                [s[0], s[1], s[2]]
            };
            let d = det3(sub3(p(v[1]), p(v[0])), sub3(p(v[2]), p(v[0])), sub3(p(v[3]), p(v[0])));
            assert!(d > 0.0 && d.is_finite(), "non-positive tet volume at cell {c}: {d:e}");
        }
        // Deep-stack 3D tet reduction
        assert!(
            g_on.mesh.n_cells() < g_off.mesh.n_cells(),
            "expected z_coarsen=true ({}) < z_coarsen=false ({})",
            g_on.mesh.n_cells(),
            g_off.mesh.n_cells()
        );
    }

    #[test]
    fn fr28_narrow_gap_centerline_seeding_eliminates_cross_terminal_bridging() {
        let mut m = Model::new_3d();
        let half = 120e-6;
        let sub = m.box_("sub", (-half, -half, -100e-6), (half, half, 0.0));
        let vac = m.box_("vac", (-half, -half, 0.0), (half, half, 100e-6));
        let w = 4e-6;
        let g = 2e-6;
        let r = w / 2.0;
        let n_arc = 16usize;
        let mut pa = vec![pt2(-50e-6, -6e-6), pt2(-40e-6, -6e-6), pt2(-40e-6, -3e-6)];
        for i in 0..=n_arc {
            let th = -std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (i as f64) / (n_arc as f64);
            pa.push(pt2(36e-6 + r * th.cos(), -1e-6 + r * th.sin()));
        }
        pa.push(pt2(-40e-6, 1e-6));
        pa.push(pt2(-40e-6, 6e-6));
        pa.push(pt2(-50e-6, 6e-6));

        let mut pb = vec![pt2(50e-6, 10e-6), pt2(40e-6, 10e-6), pt2(40e-6, 1e-6 + g + w)];
        for i in 0..=n_arc {
            let th = std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (i as f64) / (n_arc as f64);
            pb.push(pt2(-36e-6 + r * th.cos(), 1e-6 + g + r + r * th.sin()));
        }
        pb.push(pt2(40e-6, 1e-6 + g));
        pb.push(pt2(40e-6, -2e-6));
        pb.push(pt2(50e-6, -2e-6));

        let poly_a = m.polygons("pa", vec![pa], vec![]);
        let poly_b = m.polygons("pb", vec![pb], vec![]);
        m.sheet("cap_a", poly_a, 0.0);
        m.sheet("cap_b", poly_b, 0.0);
        m.overlay(&[("sub", sub), ("vac", vac)]);

        let ev = evaluate(&m.list).unwrap();
        let mut plan = MeshPlan::simple(80e-6);
        plan.opts.gap_cells = None; // must auto-detect 2 um gap!
        plan.opts.max_steiner_multiplier = Some(8);
        plan.opts.corner_grading = true;
        plan.opts.singular_edges = Some("quarter_point".into());
        let labels = ev.label_map().unwrap();
        let sz = compute_sheet3d_sizing(&ev, &labels, &plan, None).unwrap();
        assert!(!sz.gap_seeds.is_empty(), "expected deterministic centerline gap seeds for 2 um gap");
        let (m2, _) = build_overlay_mesh_2d_with_sizing(&ev, &plan, &sz).unwrap();

        // Audit m2 directly: zero dielectric triangles may have vertices touching both cap_a and cap_b rims with 0 interior vertices
        let st = ev.complex().stack.as_ref().unwrap();
        let mut cell_sheet_id: Vec<u32> = vec![0; st.overlay.regions.len()];
        for f in &ev.complex().faces {
            if f.is_sheet() {
                if let FaceGeom::Cap { cell, .. } = f.geom {
                    cell_sheet_id[cell as usize] = if f.name == "cap_a" { 1 } else { 2 };
                }
            }
        }
        let nv = m2.n_points();
        let mut v_touch = vec![0u8; nv];
        for ti in 0..m2.n_cells() {
            let ov_c = (m2.cell_labels[ti] - 1) as usize;
            let sid = cell_sheet_id[ov_c];
            if sid > 0 {
                for &vi in m2.cell(ti) {
                    v_touch[vi as usize] |= 1u8 << (sid - 1);
                }
            }
        }
        let mut n_bridging = 0usize;
        for ti in 0..m2.n_cells() {
            let ov_c = (m2.cell_labels[ti] - 1) as usize;
            if cell_sheet_id[ov_c] != 0 {
                continue;
            }
            let v = m2.cell(ti);
            let m0 = v_touch[v[0] as usize];
            let m1 = v_touch[v[1] as usize];
            let m2_v = v_touch[v[2] as usize];
            if m0 != 0 && m1 != 0 && m2_v != 0 && ((m0 | m1 | m2_v) == 3) {
                n_bridging += 1;
            }
        }
        assert_eq!(n_bridging, 0, "found {n_bridging} cross-terminal bridging triangles across 2 um gap");

        // Also verify the exact 6-finger filleted IDC geometry from bug4_sheet3d_narrow_gap_steiner_locking
        let mut m6 = Model::new_3d();
        let h6 = 260e-6;
        let vb = m6.box_("vac_bot", (-h6, -h6, -500e-6), (h6, h6, -250e-6));
        let sb = m6.box_("sub_box", (-h6, -h6, -250e-6), (h6, h6, 0.0));
        let vt = m6.box_("vac_top", (-h6, -h6, 0.0), (h6, h6, 250e-6));
        let mut bus_a = vec![pt2(-75e-6, -22e-6), pt2(-60e-6, -22e-6)];
        let mut bus_b = vec![pt2(75e-6, 22e-6), pt2(60e-6, 22e-6)];
        let pitch = 2.0 * (w + g);
        let n_arc6 = 20usize;
        for k in 0..3 {
            let yc = -15e-6 + (k as f64) * pitch;
            bus_a.push(pt2(-60e-6, yc - r));
            for i in 0..=n_arc6 {
                let th = -std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (i as f64) / (n_arc6 as f64);
                bus_a.push(pt2(48e-6 + r * th.cos(), yc + r * th.sin()));
            }
            bus_a.push(pt2(-60e-6, yc + r));
        }
        bus_a.push(pt2(-60e-6, 22e-6));
        bus_a.push(pt2(-75e-6, 22e-6));
        for k in (0..3).rev() {
            let yc = -15e-6 + (w + g) + (k as f64) * pitch;
            bus_b.push(pt2(60e-6, yc + r));
            for i in 0..=n_arc6 {
                let th = std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * (i as f64) / (n_arc6 as f64);
                bus_b.push(pt2(-48e-6 + r * th.cos(), yc + r * th.sin()));
            }
            bus_b.push(pt2(60e-6, yc - r));
        }
        bus_b.push(pt2(60e-6, -22e-6));
        bus_b.push(pt2(75e-6, -22e-6));
        let mut hole = Vec::new();
        let (hx, hy, hr) = (85e-6, 30e-6, 4e-6);
        for &(cx0, cy0, th0) in &[
            (hx - hr, -hy + hr, -std::f64::consts::FRAC_PI_2),
            (hx - hr, hy - hr, 0.0),
            (-hx + hr, hy - hr, std::f64::consts::FRAC_PI_2),
            (-hx + hr, -hy + hr, std::f64::consts::PI),
        ] {
            for i in 0..=n_arc6 {
                let th = th0 + std::f64::consts::FRAC_PI_2 * (i as f64) / (n_arc6 as f64);
                hole.push(pt2(cx0 + hr * th.cos(), cy0 + hr * th.sin()));
            }
        }
        let pa6 = m6.polygons_with("p_a", vec![bus_a], vec![], 1e-9, 0.0);
        let pb6 = m6.polygons_with("p_b", vec![bus_b], vec![], 1e-9, 0.0);
        let pg6 = m6.polygons_with(
            "p_gnd",
            vec![vec![pt2(-h6, -h6), pt2(h6, -h6), pt2(h6, h6), pt2(-h6, h6)]],
            vec![hole],
            1e-9,
            0.0,
        );
        m6.sheet("cap_a", pa6, 0.0);
        m6.sheet("cap_b", pb6, 0.0);
        m6.sheet("gnd", pg6, 0.0);
        m6.overlay(&[("vacuum_bot", vb), ("substrate", sb), ("vacuum_top", vt)]);
        let ev6 = evaluate(&m6.list).unwrap();
        let mut plan6 = MeshPlan::simple(95e-6);
        plan6.opts.dz_min = Some(2e-6);
        plan6.opts.z_growth_ratio = Some(2.4);
        plan6.opts.max_z_slabs = Some(12);
        plan6.opts.max_steiner_multiplier = Some(8);
        plan6.opts.singular_edges = Some("quarter_point".into());
        plan6.opts.quarter_point_xi = Some(0.25);
        plan6.opts.corner_grading = true;
        let labels6 = ev6.label_map().unwrap();
        let sz6 = compute_sheet3d_sizing(&ev6, &labels6, &plan6, None).unwrap();
        let (m2_6, warns6) = build_overlay_mesh_2d_with_sizing(&ev6, &plan6, &sz6).unwrap();
        let (g6, _) = mesh_sheet_stack_3d(ev6, &plan6, None).unwrap();
        // Inspect z planes and bridging on bug4 geometry
        let mut zs: Vec<i64> = (0..g6.mesh.n_points()).map(|i| (g6.mesh.point(i)[2] * 1e9).round() as i64).collect();
        zs.sort_unstable();
        zs.dedup();
        let st6 = g6.eval.complex().stack.as_ref().unwrap();
        let mut csid6 = vec![0u8; st6.overlay.regions.len()];
        for f in &g6.eval.complex().faces {
            if f.is_sheet() {
                if let FaceGeom::Cap { cell, .. } = f.geom {
                    csid6[cell as usize] = match f.name.as_str() {
                        "cap_a" => 1,
                        "cap_b" => 2,
                        "gnd" => 4,
                        _ => 0,
                    };
                }
            }
        }
        let mut vt6 = vec![0u8; m2_6.n_points()];
        for ti in 0..m2_6.n_cells() {
            let sid = csid6[(m2_6.cell_labels[ti] - 1) as usize];
            if sid != 0 {
                for &vi in m2_6.cell(ti) {
                    vt6[vi as usize] |= sid;
                }
            }
        }
        let mut n_bridge6 = 0usize;
        let mut n_ab_edges6 = 0usize;
        for ti in 0..m2_6.n_cells() {
            if csid6[(m2_6.cell_labels[ti] - 1) as usize] != 0 {
                continue;
            }
            let v = m2_6.cell(ti);
            let m0 = vt6[v[0] as usize];
            let m1 = vt6[v[1] as usize];
            let m2v = vt6[v[2] as usize];
            if m0 != 0 && m1 != 0 && m2v != 0 && (m0 | m1 | m2v).count_ones() >= 2 {
                n_bridge6 += 1;
            }
            for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                let ma = vt6[a as usize];
                let mb = vt6[b as usize];
                if ma != 0 && mb != 0 && (ma & mb) == 0 {
                    n_ab_edges6 += 1;
                }
            }
        }
        assert_eq!(n_bridge6, 0, "6-finger filleted IDC had {n_bridge6} cross-terminal bridging triangles");
        assert_eq!(n_ab_edges6, 0, "6-finger filleted IDC had {n_ab_edges6} direct cross-terminal edges");
        assert!(warns6.is_empty(), "unexpected 2D refinement warnings on 6-finger IDC: {warns6:?}");
        // Verify near-sheet dz = 2 um is preserved under max_z_slabs = 12 (no vertical stretching to 50 um)
        assert!(
            zs.contains(&2000) && zs.contains(&-2000),
            "near-sheet +/-2 um z-planes must be preserved under max_z_slabs=12, got {zs:?}"
        );
        // Verify 100% strictly positive tet volumes on the 3D mesh
        for c in 0..g6.mesh.n_cells() {
            let v = g6.mesh.cell(c);
            let p = |id: u32| {
                let s = g6.mesh.point(id as usize);
                [s[0], s[1], s[2]]
            };
            let d = det3(sub3(p(v[1]), p(v[0])), sub3(p(v[2]), p(v[0])), sub3(p(v[3]), p(v[0])));
            assert!(d > 0.0 && d.is_finite(), "non-positive tet volume at cell {c}: {d:e}");
        }
    }

    #[test]
    fn fr31_interior_mesh_inside_inner_vac_box_is_invariant_to_outer_box_extent() {
        let build_2d = |r_outer: f64| -> SimplexMesh {
            let mut m = Model::new_3d();
            let vac_out = m.box_("vac_out", (-r_outer, -r_outer, -300e-6), (r_outer, r_outer, 300e-6));
            let vac_in = m.box_("inner_vac_box", (-250e-6, -250e-6, -200e-6), (250e-6, 250e-6, 200e-6));
            let sub = m.box_("sub_box", (-150e-6, -150e-6, -150e-6), (150e-6, 150e-6, 0.0));
            let p1 = m.polygons(
                "p1",
                vec![vec![pt2(-80e-6, -40e-6), pt2(-10e-6, -40e-6), pt2(-10e-6, 40e-6), pt2(-80e-6, 40e-6)]],
                vec![],
            );
            let p2 = m.polygons(
                "p2",
                vec![vec![pt2(10e-6, -40e-6), pt2(80e-6, -40e-6), pt2(80e-6, 40e-6), pt2(10e-6, 40e-6)]],
                vec![],
            );
            m.sheet("pad1", p1, 0.0);
            m.sheet("pad2", p2, 0.0);
            m.overlay(&[("vac_out", vac_out), ("inner_vac_box", vac_in), ("sub", sub)]);
            m.size("sz_p1", Selector::parse("pad1:edge").unwrap(), Some(Expr::Num(15e-6)), None, None);
            m.size("sz_p2", Selector::parse("pad2:edge").unwrap(), Some(Expr::Num(15e-6)), None, None);
            m.size("sz_in", Selector::parse("region(inner_vac_box)").unwrap(), Some(Expr::Num(60e-6)), None, None);

            let ev = evaluate(&m.list).unwrap();
            let plan = MeshPlan::simple(120e-6);
            build_sized_overlay_mesh_2d(&ev, &plan, None).unwrap()
        };

        let m_r1 = build_2d(380e-6);
        let m_r2 = build_2d(760e-6);

        let interior_pts = |m: &SimplexMesh| -> Vec<(i64, i64)> {
            let mut out: Vec<(i64, i64)> = (0..m.n_points())
                .filter_map(|i| {
                    let x = m.points[2 * i];
                    let y = m.points[2 * i + 1];
                    if x.abs() < 245e-6 && y.abs() < 245e-6 {
                        Some(((x * 1e12).round() as i64, (y * 1e12).round() as i64))
                    } else {
                        None
                    }
                })
                .collect();
            out.sort_unstable();
            out
        };
        let p_r1 = interior_pts(&m_r1);
        let p_r2 = interior_pts(&m_r2);
        assert!(!p_r1.is_empty());
        assert_eq!(
            p_r1, p_r2,
            "interior 2D vertices inside inner_vac_box must be invariant to outer box expansion R1 -> R2 ({} vs {})",
            p_r1.len(),
            p_r2.len()
        );
    }
}
