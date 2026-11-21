//! Batch polygon preprocessing for `Op::Polygons`: closed-ring Ramer-Douglas-Peucker (RDP)
//! simplification with cross-ring shared-vertex pinning and minimum-hole-area filtering.
//!
//! When multiple polygons (such as transmon capacitor pads and a Josephson junction strap, or
//! adjacent GDSII polygons exported from KLayout/gdsfactory) touch at shared vertices or along
//! contact corners, simplifying each ring independently could drop a slightly collinear contact
//! vertex from one ring while keeping it in the other, opening a sub-tolerance gap or sliver.
//! Pinning every vertex that lies within `pin_tol` of a vertex of another ring (or an external
//! pin list from sibling features) guarantees that shared contact corners survive RDP untouched.

use crate::error::{GeomError, GeomResult};
use crate::math::*;
use crate::prim::polygon_area;

/// Unsigned shoelace area of a closed 2D polygon ring.
pub fn shoelace_area(ring: &[V2]) -> f64 {
    polygon_area(ring).abs()
}

/// Removes consecutive duplicate vertices (and a closing duplicate of the first vertex).
pub fn clean_ring(ring: &[V2], eps: f64) -> Vec<V2> {
    let mut out: Vec<V2> = Vec::with_capacity(ring.len());
    for &p in ring {
        if out.last().is_none_or(|&q| dist2(p, q) > eps) {
            out.push(p);
        }
    }
    while out.len() > 1 && dist2(*out.last().unwrap(), out[0]) <= eps {
        out.pop();
    }
    out
}

/// Distance from point `p` to the line segment `a -> b`.
fn point_segment_dist(p: V2, a: V2, b: V2) -> f64 {
    let ab = sub2(b, a);
    let l2 = dot2(ab, ab);
    let t = if l2 > 0.0 { (dot2(sub2(p, a), ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    dist2(p, add2(a, mul2(ab, t)))
}

/// Recursive Ramer-Douglas-Peucker simplification on `pts[0..=m]` (endpoints `0` and `m` are
/// always kept). Marks surviving indices in `keep`.
fn rdp_chain(pts: &[V2], indices: &[usize], tol: f64, keep: &mut [bool]) {
    if indices.len() <= 2 {
        return;
    }
    let a = pts[indices[0]];
    let b = pts[*indices.last().unwrap()];
    let mut best_d = -1.0;
    let mut best_k = 0;
    for (k, &idx) in indices.iter().enumerate().take(indices.len() - 1).skip(1) {
        let d = point_segment_dist(pts[idx], a, b);
        if d > best_d {
            best_d = d;
            best_k = k;
        }
    }
    if best_d > tol {
        keep[indices[best_k]] = true;
        rdp_chain(pts, &indices[..=best_k], tol, keep);
        rdp_chain(pts, &indices[best_k..], tol, keep);
    }
}

/// Simplifies a single closed ring with mandatory `pinned` vertices using Ramer-Douglas-Peucker
/// at tolerance `tol`. Always preserves at least 3 non-degenerate vertices.
pub fn simplify_closed_ring(ring: &[V2], pinned: &[bool], tol: f64) -> Vec<V2> {
    let n = ring.len();
    if n <= 3 || tol <= 0.0 {
        return ring.to_vec();
    }
    let mut keep = pinned.to_vec();
    keep.resize(n, false);

    // Ensure at least two convex-hull extreme anchors so we can split the closed cycle into open chains.
    let n_pinned = keep.iter().filter(|&&b| b).count();
    if n_pinned == 0 {
        // Lexicographically smallest (x, y) is an extreme vertex of the convex hull.
        let i0 = (0..n)
            .min_by(|&a, &b| ring[a][0].total_cmp(&ring[b][0]).then(ring[a][1].total_cmp(&ring[b][1])))
            .unwrap_or(0);
        keep[i0] = true;
    }
    if keep.iter().filter(|&&b| b).count() == 1 {
        let i0 = keep.iter().position(|&b| b).unwrap();
        // Farthest vertex from i0 is also an extreme vertex of the convex hull.
        let i1 = (0..n)
            .filter(|&i| i != i0)
            .max_by(|&a, &b| dist2(ring[a], ring[i0]).total_cmp(&dist2(ring[b], ring[i0])).then(b.cmp(&a)))
            .unwrap_or((i0 + 1) % n);
        keep[i1] = true;
    }

    let anchors: Vec<usize> = (0..n).filter(|&i| keep[i]).collect();
    for k in 0..anchors.len() {
        let start = anchors[k];
        let end = anchors[(k + 1) % anchors.len()];
        let mut chain = vec![start];
        let mut cur = (start + 1) % n;
        while cur != end {
            chain.push(cur);
            cur = (cur + 1) % n;
        }
        chain.push(end);
        rdp_chain(ring, &chain, tol, &mut keep);
    }

    // A closed polygon must retain at least 3 vertices.
    if keep.iter().filter(|&&b| b).count() < 3 {
        let kept: Vec<usize> = (0..n).filter(|&i| keep[i]).collect();
        let (a, b) = (ring[kept[0]], ring[kept[1 % kept.len()]]);
        if let Some(i2) = (0..n)
            .filter(|&i| !keep[i])
            .max_by(|&x, &y| point_segment_dist(ring[x], a, b).total_cmp(&point_segment_dist(ring[y], a, b)))
        {
            keep[i2] = true;
        }
    }

    (0..n).filter(|&i| keep[i]).map(|i| ring[i]).collect()
}

/// Simplified outer polygon rings and surviving hole rings.
pub type RingBatch = (Vec<Vec<V2>>, Vec<Vec<V2>>);

/// Simplifies outer polygon rings and hole rings while pinning vertices shared across distinct
/// rings (or matching `external_pins`) within `simplify_tol.max(1e-12)`, and drops any hole ring
/// whose unsigned shoelace area is `< min_hole_area`.
pub fn simplify_rings_with_pins(
    polygons: &[Vec<V2>],
    holes: &[Vec<V2>],
    external_pins: &[V2],
    simplify_tol: f64,
    min_hole_area: f64,
) -> GeomResult<RingBatch> {
    if simplify_tol < 0.0 || !simplify_tol.is_finite() {
        return Err(GeomError(format!("simplify_tol ({simplify_tol}) must be non-negative and finite")));
    }
    if min_hole_area < 0.0 || !min_hole_area.is_finite() {
        return Err(GeomError(format!("min_hole_area ({min_hole_area}) must be non-negative and finite")));
    }
    let clean_eps = (simplify_tol * 1e-3).clamp(1e-15, 1e-11);
    let mut clean_outers: Vec<Vec<V2>> = Vec::with_capacity(polygons.len());
    for (i, r) in polygons.iter().enumerate() {
        let c = clean_ring(r, clean_eps);
        if c.len() < 3 {
            return Err(GeomError(format!("outer polygon ring {i} has fewer than 3 distinct points")));
        }
        if shoelace_area(&c) <= 1e-30 {
            return Err(GeomError(format!("outer polygon ring {i} has zero area")));
        }
        clean_outers.push(c);
    }

    let mut clean_holes: Vec<Vec<V2>> = Vec::new();
    for (i, r) in holes.iter().enumerate() {
        let c = clean_ring(r, clean_eps);
        if c.len() < 3 {
            if min_hole_area > 0.0 {
                continue;
            }
            return Err(GeomError(format!("hole ring {i} has fewer than 3 distinct points")));
        }
        let area = shoelace_area(&c);
        if (min_hole_area > 0.0 && area < min_hole_area) || area <= 1e-30 {
            continue;
        }
        clean_holes.push(c);
    }

    let pin_tol = simplify_tol.max(1e-12);
    let all_rings: Vec<&[V2]> =
        clean_outers.iter().map(Vec::as_slice).chain(clean_holes.iter().map(Vec::as_slice)).collect();

    let is_pinned = |ring_idx: usize, p: V2| -> bool {
        if external_pins.iter().any(|&q| dist2(p, q) <= pin_tol) {
            return true;
        }
        for (other_idx, other) in all_rings.iter().enumerate() {
            if other_idx == ring_idx {
                continue;
            }
            if other.iter().any(|&q| dist2(p, q) <= pin_tol) {
                return true;
            }
        }
        false
    };

    let mut out_polys = Vec::with_capacity(clean_outers.len());
    for (ri, r) in clean_outers.iter().enumerate() {
        let pinned: Vec<bool> = r.iter().map(|&p| is_pinned(ri, p)).collect();
        let simp = simplify_closed_ring(r, &pinned, simplify_tol);
        if simp.len() < 3 || shoelace_area(&simp) <= 1e-30 {
            return Err(GeomError(format!(
                "outer polygon ring {ri} collapsed under simplify_tol = {simplify_tol:.3e}"
            )));
        }
        out_polys.push(simp);
    }

    let n_out = clean_outers.len();
    let mut out_holes = Vec::with_capacity(clean_holes.len());
    for (hi, r) in clean_holes.iter().enumerate() {
        let ri = n_out + hi;
        let pinned: Vec<bool> = r.iter().map(|&p| is_pinned(ri, p)).collect();
        let simp = simplify_closed_ring(r, &pinned, simplify_tol);
        let area = shoelace_area(&simp);
        if simp.len() >= 3 && area > 1e-30 && (min_hole_area == 0.0 || area >= min_hole_area) {
            out_holes.push(simp);
        }
    }

    Ok((out_polys, out_holes))
}

/// Simplifies outer polygon rings and hole rings with cross-ring vertex pinning and hole area filtering.
pub fn simplify_rings(
    polygons: &[Vec<V2>],
    holes: &[Vec<V2>],
    simplify_tol: f64,
    min_hole_area: f64,
) -> GeomResult<RingBatch> {
    simplify_rings_with_pins(polygons, holes, &[], simplify_tol, min_hole_area)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rdp_preserves_shared_contact_corners_and_filters_small_holes() {
        // Pad ring has a slightly deflected contact vertex at (0.0, 20e-9) shared with a junction strap.
        // Without cross-ring pinning, RDP with simplify_tol = 50e-9 would drop (0.0, 20e-9) from the pad
        // because 20 nm < 50 nm from the (-10e-6, 0) -> (10e-6, 0) chord.
        let pad = vec![
            [-10e-6, 0.0],
            [-5e-6, 10e-9], // unshared jitter (10 nm < 50 nm): should be simplified away
            [0.0, 20e-9],   // shared with jj strap (20 nm < 50 nm): MUST be pinned and preserved!
            [10e-6, 0.0],
            [10e-6, 10e-6],
            [-10e-6, 10e-6],
        ];
        let jj = vec![[-1e-6, -5e-6], [1e-6, -5e-6], [0.0, 20e-9]];
        // Two holes: one tiny (1 um x 1 um = 1e-12 m^2), one large (4 um x 4 um = 16e-12 m^2)
        let small_hole = vec![[-2e-6, 2e-6], [-1e-6, 2e-6], [-1e-6, 3e-6], [-2e-6, 3e-6]];
        let big_hole = vec![[2e-6, 2e-6], [6e-6, 2e-6], [6e-6, 6e-6], [2e-6, 6e-6]];

        let (polys, holes) = simplify_rings(&[pad, jj], &[small_hole, big_hole], 50e-9, 4e-12).unwrap();
        assert_eq!(polys.len(), 2);
        // Pad kept the 4 main corners + the shared contact vertex [0.0, 20e-9], but dropped [-5e-6, 10e-9]
        assert_eq!(polys[0].len(), 5, "simplified pad: {:?}", polys[0]);
        assert!(polys[0].contains(&[0.0, 20e-9]), "shared contact corner must be pinned");
        assert!(!polys[0].contains(&[-5e-6, 10e-9]), "unshared collinear jitter must be removed");
        // Only the big hole (16e-12 >= 4e-12) survives
        assert_eq!(holes.len(), 1);
        assert!((shoelace_area(&holes[0]) - 16e-12).abs() < 1e-24);
    }
}
