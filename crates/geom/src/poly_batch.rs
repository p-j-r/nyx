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
use crate::prim::{polygon_area, DGeom, DSeg};

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

    let cell_of = |p: V2| -> (i64, i64) { ((p[0] / pin_tol).floor() as i64, (p[1] / pin_tol).floor() as i64) };
    let mut pin_grid: std::collections::HashMap<(i64, i64), Vec<(usize, V2)>> = std::collections::HashMap::new();
    for &q in external_pins {
        pin_grid.entry(cell_of(q)).or_default().push((usize::MAX, q));
    }
    for (ring_idx, ring) in all_rings.iter().enumerate() {
        for &q in *ring {
            pin_grid.entry(cell_of(q)).or_default().push((ring_idx, q));
        }
    }

    let is_pinned = |ring_idx: usize, p: V2| -> bool {
        let (cx, cy) = cell_of(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(bucket) = pin_grid.get(&(cx + dx, cy + dy)) {
                    if bucket.iter().any(|&(other_idx, q)| other_idx != ring_idx && dist2(p, q) <= pin_tol) {
                        return true;
                    }
                }
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

/// Fitted outer contour segment lists and surviving hole contour segment lists.
pub type ContourBatch = (Vec<Vec<DSeg>>, Vec<Vec<DSeg>>);

/// Corner turn-angle threshold (25 degrees) separating smooth circular arcs from polygon corners.
const ARC_CORNER_RAD: f64 = 0.436_332_312_998_582_4;

/// Evaluates the exact start and end coordinates of `DGeom::Arc { c, r, a0, sweep }` matching `geom::arrange`.
pub fn arc_endpoints(c: V2, r: f64, a0: f64, sweep: f64) -> (V2, V2) {
    let eval = |t: f64| [c[0] + r * t.cos(), c[1] + r * t.sin()];
    if sweep >= 0.0 {
        let a_lo = a0.rem_euclid(TAU);
        (eval(a_lo), eval(a_lo + sweep))
    } else {
        let a_lo = (a0 + sweep).rem_euclid(TAU);
        (eval(a_lo - sweep), eval(a_lo))
    }
}

/// Fits a closed full circle (`sweep = TAU`) to an unpinned closed CCW ring if all vertices lie within `arc_tol`.
fn fit_closed_full_circle(pts: &[V2], arc_tol: f64) -> Option<DSeg> {
    let n = pts.len();
    if n < 6 || arc_tol <= 0.0 || !arc_tol.is_finite() {
        return None;
    }
    let mut lens = Vec::with_capacity(n);
    for i in 0..n {
        lens.push(dist2(pts[i], pts[(i + 1) % n]));
    }
    let mut sorted_lens = lens.clone();
    sorted_lens.sort_by(|a, b| a.total_cmp(b));
    let med_len = sorted_lens[n / 2];
    if med_len <= 0.0 || lens.iter().any(|&l| l > 3.5 * med_len && l > 4.0 * arc_tol) {
        return None;
    }

    let mean = mul2(pts.iter().fold([0.0, 0.0], |acc, &p| add2(acc, p)), 1.0 / n as f64);
    let mut qs = Vec::with_capacity(n);
    let mut zs = Vec::with_capacity(n);
    let mut z_mean = 0.0;
    for &p in pts {
        let q = sub2(p, mean);
        let z = dot2(q, q);
        qs.push(q);
        zs.push(z);
        z_mean += z;
    }
    z_mean /= n as f64;

    let (mut sxx, mut sxy, mut syy, mut sxz, mut syz) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (q, &z) in qs.iter().zip(&zs) {
        let dz = z - z_mean;
        sxx += q[0] * q[0];
        sxy += q[0] * q[1];
        syy += q[1] * q[1];
        sxz += 0.5 * q[0] * dz;
        syz += 0.5 * q[1] * dz;
    }
    let det = sxx * syy - sxy * sxy;
    if det.abs() <= 1e-30 {
        return None;
    }
    let mut c = [mean[0] + (syy * sxz - sxy * syz) / det, mean[1] + (sxx * syz - sxy * sxz) / det];

    for _ in 0..8 {
        let r: f64 = pts.iter().map(|&p| dist2(p, c)).sum::<f64>() / n as f64;
        if r <= 1e-15 {
            return None;
        }
        let (mut j00, mut j01, mut j11, mut g0, mut g1) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for &p in pts {
            let d = sub2(c, p);
            let rho = norm2(d);
            if rho <= 1e-15 {
                continue;
            }
            let ux = d[0] / rho;
            let uy = d[1] / rho;
            let f = rho - r;
            j00 += ux * ux;
            j01 += ux * uy;
            j11 += uy * uy;
            g0 += f * ux;
            g1 += f * uy;
        }
        let d_j = j00 * j11 - j01 * j01;
        if d_j.abs() <= 1e-30 {
            break;
        }
        let dx = (j11 * g0 - j01 * g1) / d_j;
        let dy = (j00 * g1 - j01 * g0) / d_j;
        c[0] -= dx;
        c[1] -= dy;
        if dx.hypot(dy) <= 1e-14 * r {
            break;
        }
    }

    let r = pts.iter().map(|&p| dist2(p, c)).sum::<f64>() / n as f64;
    if !r.is_finite() || r <= arc_tol {
        return None;
    }
    let max_err = pts.iter().map(|&p| (dist2(p, c) - r).abs()).fold(0.0f64, f64::max);
    if max_err > arc_tol {
        return None;
    }

    let mut sum_dphi = 0.0f64;
    for i in 0..n {
        let d0 = sub2(pts[i], c);
        let d1 = sub2(pts[(i + 1) % n], c);
        let dphi = cross2(d0, d1).atan2(dot2(d0, d1));
        if !dphi.is_finite() || dphi <= 1e-9 || dphi >= std::f64::consts::FRAC_PI_2 {
            return None;
        }
        sum_dphi += dphi;
    }
    if (sum_dphi - TAU).abs() > 1e-2 {
        return None;
    }
    Some(DSeg { geom: DGeom::Arc { c, r, a0: 0.0, sweep: TAU }, role: "e0".to_string() })
}

/// Endpoint-anchored 1-DOF circular arc fitter on an open chain `chain[0..k]` (`k >= 4`).
/// Constrains the circle center strictly to the perpendicular bisector of `chain[0] -> chain[k-1]`
/// (`c(d) = m + d * n_perp`, `R(d) = ||p_0 - c(d)||`) and optimizes scalar `d` via 1D exact algebraic
/// (Taubin) seed + 1D Gauss-Newton iterations.
fn fit_open_arc_1dof(chain: &[V2], arc_tol: f64) -> Option<(V2, f64, f64, f64)> {
    let k = chain.len();
    if k < 4 || arc_tol <= 0.0 || !arc_tol.is_finite() {
        return None;
    }
    let (p0, p1) = (chain[0], chain[k - 1]);
    let l_chord = dist2(p0, p1);
    if l_chord <= 1e-15 {
        return None;
    }
    let mut max_line_dev = 0.0f64;
    let mut path_len = 0.0f64;
    let mut lens = Vec::with_capacity(k - 1);
    for i in 0..(k - 1) {
        let sl = dist2(chain[i], chain[i + 1]);
        path_len += sl;
        lens.push(sl);
        if i > 0 {
            max_line_dev = max_line_dev.max(point_segment_dist(chain[i], p0, p1));
        }
    }
    if max_line_dev <= arc_tol && path_len <= 1.5 * l_chord {
        return None;
    }

    let mut sorted_lens = lens.clone();
    sorted_lens.sort_by(|a, b| a.total_cmp(b));
    let med_len = sorted_lens[sorted_lens.len() / 2];
    if med_len <= 0.0 || lens.iter().any(|&l| l > 3.5 * med_len && l > 4.0 * arc_tol) {
        return None;
    }

    let m = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])];
    let u = sub2(p1, p0);
    let a = 0.5 * l_chord;
    let u_hat = [u[0] / l_chord, u[1] / l_chord];
    let n_perp = [-u_hat[1], u_hat[0]];

    // 1D algebraic (Taubin/Kasa) seed along the perpendicular bisector:
    // ||p_i - c(d)||^2 - R(d)^2 = (||q_i||^2 - a^2) - 2 * d * (q_i . n_perp)
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for &pi in &chain[1..(k - 1)] {
        let qi = sub2(pi, m);
        let vi = dot2(qi, n_perp);
        let ri2 = dot2(qi, qi);
        num += vi * (ri2 - a * a);
        den += 2.0 * vi * vi;
    }
    if den <= 1e-30 || !num.is_finite() || !den.is_finite() {
        return None;
    }
    let mut d = num / den;

    let eval_loss = |d_val: f64| -> f64 {
        let r_val = a.hypot(d_val);
        let cd = add2(m, mul2(n_perp, d_val));
        chain[1..(k - 1)]
            .iter()
            .map(|&pi| {
                let f = dist2(pi, cd) - r_val;
                f * f
            })
            .sum()
    };
    let mut cur_loss = eval_loss(d);
    for _ in 0..8 {
        let r_val = a.hypot(d);
        if r_val <= 1e-15 {
            break;
        }
        let mut g = 0.0f64;
        let mut h = 0.0f64;
        for &pi in &chain[1..(k - 1)] {
            let qi = sub2(pi, m);
            let vi = dot2(qi, n_perp);
            let rho = dist2(qi, mul2(n_perp, d));
            if rho <= 1e-15 {
                continue;
            }
            let f = rho - r_val;
            let df = (d - vi) / rho - d / r_val;
            g += f * df;
            h += df * df;
        }
        if h <= 1e-30 {
            break;
        }
        let step = g / h;
        if !step.is_finite() {
            break;
        }
        let mut accepted = false;
        for alpha in [1.0, 0.5, 0.25] {
            let cand_d = d - alpha * step;
            let cand_loss = eval_loss(cand_d);
            if cand_loss <= cur_loss {
                d = cand_d;
                cur_loss = cand_loss;
                accepted = true;
                break;
            }
        }
        if !accepted || step.abs() <= 1e-14 * r_val {
            break;
        }
    }

    let mut c = add2(m, mul2(n_perp, d));
    // Eliminate 0.5-ULP midpoint rounding along u_hat so ||p0 - c|| == ||p1 - c|| to the last bit
    let r0_sq = dot2(sub2(p0, c), sub2(p0, c));
    let r1_sq = dot2(sub2(p1, c), sub2(p1, c));
    let du = (r1_sq - r0_sq) / (2.0 * l_chord);
    if du.is_finite() {
        c = add2(c, mul2(u_hat, du));
    }
    let r = 0.5 * (dist2(p0, c) + dist2(p1, c));
    if !r.is_finite() || r <= arc_tol {
        return None;
    }

    let max_err = chain.iter().map(|&p| (dist2(p, c) - r).abs()).fold(0.0f64, f64::max);
    if max_err > arc_tol {
        return None;
    }

    let mut sum_dphi = 0.0f64;
    let mut dir_sign = 0.0f64;
    for j in 0..(k - 1) {
        let dj0 = sub2(chain[j], c);
        let dj1 = sub2(chain[j + 1], c);
        let dphi = cross2(dj0, dj1).atan2(dot2(dj0, dj1));
        if !dphi.is_finite() || dphi.abs() <= 1e-9 || dphi.abs() >= std::f64::consts::FRAC_PI_2 {
            return None;
        }
        if dir_sign == 0.0 {
            dir_sign = dphi.signum();
        } else if dphi * dir_sign <= 0.0 {
            return None;
        }
        sum_dphi += dphi;
    }
    let abs_sweep = sum_dphi.abs();
    if !(abs_sweep > 1e-3 && abs_sweep < TAU - 1e-6) {
        return None;
    }
    let sagitta = if abs_sweep >= std::f64::consts::PI { r + d.abs() } else { r * (1.0 - (0.5 * abs_sweep).cos()) };
    if sagitta <= arc_tol {
        return None;
    }

    let d0 = sub2(p0, c);
    let d1 = sub2(p1, c);
    let a0 = d0[1].atan2(d0[0]).rem_euclid(TAU);
    let a1 = d1[1].atan2(d1[0]).rem_euclid(TAU);
    let sweep = if sum_dphi > 0.0 { (a1 - a0).rem_euclid(TAU) } else { -((a0 - a1).rem_euclid(TAU)) };
    if (sweep - sum_dphi).abs() > 1e-2 {
        return None;
    }
    Some((c, r, a0, sweep))
}

/// Identifies mandatory break indices (pinned vertices, sharp corners > 25 deg, and tangential
/// long-edge transitions) on a closed CCW ring `pts`.
fn find_chain_breaks(pts: &[V2], pins: &[bool], arc_tol: f64) -> Vec<usize> {
    let n = pts.len();
    let mut is_break = vec![false; n];
    for i in 0..n {
        if pins.get(i).copied().unwrap_or(false) {
            is_break[i] = true;
            continue;
        }
        let prev = pts[(i + n - 1) % n];
        let cur = pts[i];
        let next = pts[(i + 1) % n];
        let d_in = sub2(cur, prev);
        let d_out = sub2(next, cur);
        let l_in = norm2(d_in);
        let l_out = norm2(d_out);
        if l_in <= 1e-15 || l_out <= 1e-15 {
            is_break[i] = true;
            continue;
        }
        let turn = cross2(d_in, d_out).atan2(dot2(d_in, d_out)).abs();
        if turn > ARC_CORNER_RAD {
            is_break[i] = true;
            continue;
        }
        // Tangential transition between a long straight edge and short arc chords
        let l_min = l_in.min(l_out);
        let l_max = l_in.max(l_out);
        if l_min > 4.0 * arc_tol && l_max > 3.5 * l_min {
            is_break[i] = true;
        }
    }

    let mut breaks: Vec<usize> = (0..n).filter(|&i| is_break[i]).collect();
    if breaks.len() == 1 {
        let i0 = breaks[0];
        let i1 = (0..n)
            .filter(|&i| i != i0)
            .max_by(|&a, &b| dist2(pts[a], pts[i0]).total_cmp(&dist2(pts[b], pts[i0])).then(b.cmp(&a)))
            .unwrap_or((i0 + n / 2) % n);
        is_break[i1] = true;
        breaks = (0..n).filter(|&i| is_break[i]).collect();
    } else if breaks.is_empty() {
        let i0 =
            (0..n).min_by(|&a, &b| pts[a][0].total_cmp(&pts[b][0]).then(pts[a][1].total_cmp(&pts[b][1]))).unwrap_or(0);
        let i1 = (0..n)
            .filter(|&i| i != i0)
            .max_by(|&a, &b| dist2(pts[a], pts[i0]).total_cmp(&dist2(pts[b], pts[i0])).then(b.cmp(&a)))
            .unwrap_or((i0 + n / 2) % n);
        is_break[i0] = true;
        is_break[i1] = true;
        breaks = (0..n).filter(|&i| is_break[i]).collect();
    }
    breaks
}

/// Emits simplified straight `DGeom::Line` segments for a non-arc open chain `chain[0..m]` at tolerance `line_tol`.
fn emit_simplified_lines(chain: &[V2], line_tol: f64, out: &mut Vec<DGeom>) {
    let m = chain.len();
    if m < 2 {
        return;
    }
    if m == 2 || line_tol <= 0.0 {
        for w in chain.windows(2) {
            if dist2(w[0], w[1]) > 0.0 {
                out.push(DGeom::Line { p0: w[0], p1: w[1] });
            }
        }
        return;
    }
    let indices: Vec<usize> = (0..m).collect();
    let mut keep = vec![false; m];
    keep[0] = true;
    keep[m - 1] = true;
    rdp_chain(chain, &indices, line_tol, &mut keep);
    let kept: Vec<V2> = (0..m).filter(|&i| keep[i]).map(|i| chain[i]).collect();
    for w in kept.windows(2) {
        if dist2(w[0], w[1]) > 0.0 {
            out.push(DGeom::Line { p0: w[0], p1: w[1] });
        }
    }
}

/// Segments an open chain `chain` between two mandatory break vertices into `DGeom::Arc` and `DGeom::Line` pieces.
fn segment_open_chain(chain: &[V2], arc_tol: f64, line_tol: f64, out: &mut Vec<DGeom>) {
    let m = chain.len();
    if m < 2 {
        return;
    }
    if m == 2 {
        out.push(DGeom::Line { p0: chain[0], p1: chain[1] });
        return;
    }
    if let Some((c, r, a0, sweep)) = fit_open_arc_1dof(chain, arc_tol) {
        out.push(DGeom::Arc { c, r, a0, sweep });
        return;
    }

    // Check if the entire chain is collinear within line_tol.
    let max_dev = (1..(m - 1)).map(|i| point_segment_dist(chain[i], chain[0], chain[m - 1])).fold(0.0f64, f64::max);
    if max_dev <= line_tol.max(arc_tol) {
        out.push(DGeom::Line { p0: chain[0], p1: chain[m - 1] });
        return;
    }

    // Compound smooth chain: greedily extract maximal valid arcs (>= 4 vertices) and simplify remaining linear runs.
    let mut s = 0usize;
    let mut line_start = 0usize;
    while s + 1 < m {
        let mut best_arc: Option<(usize, (V2, f64, f64, f64))> = None;
        if s + 3 < m {
            for e in ((s + 3)..m).rev() {
                if let Some(arc) = fit_open_arc_1dof(&chain[s..=e], arc_tol) {
                    best_arc = Some((e, arc));
                    break;
                }
            }
        }
        if let Some((e, (c, r, a0, sweep))) = best_arc {
            if s > line_start {
                emit_simplified_lines(&chain[line_start..=s], line_tol, out);
            }
            out.push(DGeom::Arc { c, r, a0, sweep });
            s = e;
            line_start = e;
        } else {
            s += 1;
        }
    }
    if line_start + 1 < m {
        emit_simplified_lines(&chain[line_start..m], line_tol, out);
    }
}

/// Fits circular arcs (`DGeom::Arc` within `arc_tol`) and simplified straight lines (`DGeom::Line` within `line_tol`)
/// to a closed polygon ring while respecting `pinned` vertices as mandatory chain delimiters across wrap-around.
pub fn fit_contour_arcs_with_line_tol(ring: &[V2], pinned: &[bool], arc_tol: f64, line_tol: f64) -> Vec<DSeg> {
    let n = ring.len();
    if n < 3 {
        return Vec::new();
    }
    let mut pts = ring.to_vec();
    let mut pins = pinned.to_vec();
    pins.resize(n, false);
    if polygon_area(&pts) < 0.0 {
        pts.reverse();
        pins.reverse();
    }

    if arc_tol <= 0.0 || !arc_tol.is_finite() {
        let simp = simplify_closed_ring(&pts, &pins, line_tol);
        return (0..simp.len())
            .map(|i| DSeg { geom: DGeom::Line { p0: simp[i], p1: simp[(i + 1) % simp.len()] }, role: format!("e{i}") })
            .collect();
    }

    if !pins.iter().any(|&b| b) {
        if let Some(full) = fit_closed_full_circle(&pts, arc_tol) {
            return vec![full];
        }
    }

    let breaks = find_chain_breaks(&pts, &pins, arc_tol);
    let mut geoms: Vec<DGeom> = Vec::new();
    for k in 0..breaks.len() {
        let start = breaks[k];
        let end = breaks[(k + 1) % breaks.len()];
        let mut chain = vec![pts[start]];
        let mut cur = (start + 1) % n;
        while cur != end {
            chain.push(pts[cur]);
            cur = (cur + 1) % n;
        }
        chain.push(pts[end]);
        segment_open_chain(&chain, arc_tol, line_tol, &mut geoms);
    }

    // Snap adjacent DGeom::Line endpoints to exact evaluated DGeom::Arc endpoints so ||arc_end - line_end|| == 0.0.
    let m = geoms.len();
    if m >= 2 {
        let arc_ends: Vec<Option<(V2, V2)>> = geoms
            .iter()
            .map(|g| match g {
                DGeom::Arc { c, r, a0, sweep } => Some(arc_endpoints(*c, *r, *a0, *sweep)),
                _ => None,
            })
            .collect();
        for (i, g_item) in geoms.iter_mut().enumerate().take(m) {
            if let DGeom::Line { p0, p1 } = g_item {
                let prev = (i + m - 1) % m;
                let next = (i + 1) % m;
                if let Some((_, prev_end)) = arc_ends[prev] {
                    *p0 = prev_end;
                }
                if let Some((next_start, _)) = arc_ends[next] {
                    *p1 = next_start;
                }
            }
        }
    }

    geoms.into_iter().enumerate().map(|(i, geom)| DSeg { geom, role: format!("e{i}") }).collect()
}

/// Endpoint-anchored 1-DOF perpendicular-bisector circular arc fitter for a closed ring with mandatory `pinned` vertices.
pub fn fit_contour_arcs(ring: &[V2], pinned: &[bool], arc_tol: f64) -> Vec<DSeg> {
    fit_contour_arcs_with_line_tol(ring, pinned, arc_tol, arc_tol)
}

/// Cleans outer polygon rings and hole rings, filters holes with area `< min_hole_area`, pins shared
/// vertices across rings and `external_pins`, and fits `LShape::Contour` segments (`DGeom::Arc` + `DGeom::Line`).
pub fn fit_rings_with_pins(
    polygons: &[Vec<V2>],
    holes: &[Vec<V2>],
    external_pins: &[V2],
    simplify_tol: f64,
    min_hole_area: f64,
    arc_tol: f64,
) -> GeomResult<ContourBatch> {
    if simplify_tol < 0.0 || !simplify_tol.is_finite() {
        return Err(GeomError(format!("simplify_tol ({simplify_tol}) must be non-negative and finite")));
    }
    if min_hole_area < 0.0 || !min_hole_area.is_finite() {
        return Err(GeomError(format!("min_hole_area ({min_hole_area}) must be non-negative and finite")));
    }
    if arc_tol < 0.0 || !arc_tol.is_finite() {
        return Err(GeomError(format!("arc_tol ({arc_tol}) must be non-negative and finite")));
    }
    let clean_eps = (simplify_tol.max(arc_tol) * 1e-3).clamp(1e-15, 1e-11);
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

    let pin_tol = simplify_tol.max(arc_tol).max(1e-12);
    let all_rings: Vec<&[V2]> =
        clean_outers.iter().map(Vec::as_slice).chain(clean_holes.iter().map(Vec::as_slice)).collect();

    let cell_of = |p: V2| -> (i64, i64) { ((p[0] / pin_tol).floor() as i64, (p[1] / pin_tol).floor() as i64) };
    let mut pin_grid: std::collections::HashMap<(i64, i64), Vec<(usize, V2)>> = std::collections::HashMap::new();
    for &q in external_pins {
        pin_grid.entry(cell_of(q)).or_default().push((usize::MAX, q));
    }
    for (ring_idx, ring) in all_rings.iter().enumerate() {
        for &q in *ring {
            pin_grid.entry(cell_of(q)).or_default().push((ring_idx, q));
        }
    }

    let is_pinned = |ring_idx: usize, p: V2| -> bool {
        let (cx, cy) = cell_of(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(bucket) = pin_grid.get(&(cx + dx, cy + dy)) {
                    if bucket.iter().any(|&(other_idx, q)| other_idx != ring_idx && dist2(p, q) <= pin_tol) {
                        return true;
                    }
                }
            }
        }
        false
    };

    let line_tol = if simplify_tol > 0.0 { simplify_tol } else { arc_tol };
    let mut out_contours = Vec::with_capacity(clean_outers.len());
    for (ri, r) in clean_outers.iter().enumerate() {
        let pinned: Vec<bool> = r.iter().map(|&p| is_pinned(ri, p)).collect();
        let segs = fit_contour_arcs_with_line_tol(r, &pinned, arc_tol, line_tol);
        if segs.is_empty() {
            return Err(GeomError(format!("outer polygon ring {ri} collapsed under arc_tol = {arc_tol:.3e}")));
        }
        out_contours.push(segs);
    }

    let n_out = clean_outers.len();
    let mut out_holes = Vec::with_capacity(clean_holes.len());
    for (hi, r) in clean_holes.iter().enumerate() {
        let ri = n_out + hi;
        let pinned: Vec<bool> = r.iter().map(|&p| is_pinned(ri, p)).collect();
        let segs = fit_contour_arcs_with_line_tol(r, &pinned, arc_tol, line_tol);
        if !segs.is_empty() {
            out_holes.push(segs);
        }
    }

    Ok((out_contours, out_holes))
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
