//! Local mesh operators for 2D metric adaptation: edge split, edge collapse, edge flip, and
//! length-equalising smoothing.
//!
//! Every pass takes a snapshot of the topology, orders its candidates deterministically (by
//! decreasing or increasing metric length, ties by edge id) and applies them greedily. The
//! triangles an operator touches are *locked* for the rest of the pass, so two operators in one
//! pass never overlap (the same scheme extends to independent sets for threading). All operators
//! are quality gated: the result must have positive area everywhere and must not create
//! triangles worse than the thresholds in [`AdaptOptions`], and constrained edges, corners and
//! frozen chains are never violated. The rules are the standard ones in the metric-based
//! adaptation literature (length thresholds `1/sqrt 2` and `sqrt 2`, link condition for collapses,
//! quality-improving flips).

use super::drive::AdaptOptions;
use super::field::MetricField;
use super::geom::Geometry2;
use super::metric::{Sym, edge_length};
use super::topo::{Edge, Snap, VClass, Work2, ekey, tri_area};
use crate::meshers::quality::simplex_quality;

/// Metric length of the edge `(a, b)`.
pub fn elen(w: &Work2, a: u32, b: u32) -> f64 {
    let (pa, pb) = (w.pts[a as usize], w.pts[b as usize]);
    edge_length(&w.met[a as usize], &w.met[b as usize], &[pb[0] - pa[0], pb[1] - pa[1]])
}

/// Mean-ratio quality in the metric of a triangle given by positions and vertex metrics.
pub fn tri_quality(p: [[f64; 2]; 3], m: [&Sym<2>; 3]) -> f64 {
    simplex_quality::<2>(&p, Some(&[*m[0], *m[1], *m[2]]))
}

fn tq(w: &Work2, t: [u32; 3]) -> f64 {
    let g = |i: usize| w.pts[t[i] as usize];
    tri_quality([g(0), g(1), g(2)], [&w.met[t[0] as usize], &w.met[t[1] as usize], &w.met[t[2] as usize]])
}

fn lock_get(lock: &[bool], t: i32) -> bool {
    t >= 0 && lock[t as usize]
}

/// Split every edge longer than `opt.lmax`. Returns the number of splits.
pub fn split_pass(
    w: &mut Work2,
    field: &dyn MetricField<2>,
    geom: &dyn Geometry2,
    opt: &AdaptOptions,
    sample_mid: bool,
) -> usize {
    let snap = w.snapshot();
    let mut cand: Vec<(f64, usize)> = Vec::new();
    for (i, (e, _)) in snap.edges.iter().enumerate() {
        let mut l = elen(w, e[0], e[1]);
        if sample_mid && l > 0.3 && l <= opt.lmax {
            // a feature narrower than the edge can sit between its endpoints: measure the edge as
            // two segments through the field value at its midpoint
            let (pa, pb) = (w.pts[e[0] as usize], w.pts[e[1] as usize]);
            let m = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
            let mm = field.at(&m);
            let h = [0.5 * (pb[0] - pa[0]), 0.5 * (pb[1] - pa[1])];
            l = edge_length(&w.met[e[0] as usize], &mm, &h) + edge_length(&mm, &w.met[e[1] as usize], &h);
        }
        if l > opt.lmax {
            if w.chain_of.get(e).is_some_and(|&c| w.chains[c as usize].frozen) {
                continue;
            }
            cand.push((l, i));
        }
    }
    cand.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
    let mut lock = vec![false; w.tris.len()];
    let mut n = 0;
    for (_, i) in cand {
        let (e, ts) = snap.edges[i];
        if lock_get(&lock, ts[0]) || lock_get(&lock, ts[1]) {
            continue;
        }
        // A triangle that contains a frozen edge longer than `lmax` can never get its other two
        // edges into band (triangle inequality), so splitting them only produces slivers forever.
        let blocked = ts.iter().filter(|&&t| t >= 0).any(|&t| {
            let tr = w.tris[t as usize];
            (0..3).any(|k| {
                let (a, b) = (tr[k], tr[(k + 1) % 3]);
                w.chain_of.get(&ekey(a, b)).is_some_and(|&c| w.chains[c as usize].frozen) && elen(w, a, b) > opt.lmax
            })
        });
        if blocked {
            continue;
        }
        let (pa, pb) = (w.pts[e[0] as usize], w.pts[e[1] as usize]);
        let chain = w.chain_of.get(&e).copied();
        let mut m = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
        if let Some(c) = chain {
            m = geom.project(m, c);
        }
        let met_m = field.at(&m);
        // validity of the sub-triangles
        let mut parts: Vec<(usize, [u32; 3])> = Vec::new(); // (tri, [apex, p, q])
        let mut ok = true;
        for &t in ts.iter().filter(|&&t| t >= 0) {
            let tr = w.tris[t as usize];
            let k = (0..3).find(|&k| ekey(tr[(k + 1) % 3], tr[(k + 2) % 3]) == e).unwrap();
            let (apex, p, q) = (tr[k], tr[(k + 1) % 3], tr[(k + 2) % 3]);
            let ar = |x: [f64; 2], y: [f64; 2], z: [f64; 2]| {
                0.5 * ((y[0] - x[0]) * (z[1] - x[1]) - (y[1] - x[1]) * (z[0] - x[0]))
            };
            let (pa_, pp, pq) = (w.pts[apex as usize], w.pts[p as usize], w.pts[q as usize]);
            let old = tri_area(&w.pts, &tr);
            if ar(pa_, pp, m) < 1e-6 * old || ar(pa_, m, pq) < 1e-6 * old {
                ok = false;
                break;
            }
            parts.push((t as usize, [apex, p, q]));
        }
        if !ok {
            continue;
        }
        let mid = w.pts.len() as u32;
        w.pts.push(m);
        w.met.push(met_m);
        w.vclass.push(match chain {
            Some(c) => VClass::Chain(c),
            None => VClass::Interior,
        });
        for (t, [apex, p, q]) in parts {
            let label = w.tlabel[t];
            w.tris[t] = [apex, p, mid];
            w.tris.push([apex, mid, q]);
            w.tlabel.push(label);
            w.talive.push(true);
            w.touched[t] = true;
            w.touched.push(true);
            lock[t] = true;
            lock.push(true);
        }
        if let Some(c) = chain {
            w.chain_of.remove(&e);
            w.chain_of.insert(ekey(e[0], mid), c);
            w.chain_of.insert(ekey(mid, e[1]), c);
        }
        n += 1;
    }
    n
}

fn neighbours(w: &Work2, snap: &Snap, v: u32) -> Vec<u32> {
    let mut n: Vec<u32> = Vec::new();
    for &t in snap.star(v) {
        for &x in &w.tris[t as usize] {
            if x != v {
                n.push(x);
            }
        }
    }
    n.sort_unstable();
    n.dedup();
    n
}

/// Can vertex `a` be merged into `b`? Returns the worst new quality if so.
fn collapse_eval(w: &Work2, snap: &Snap, a: u32, b: u32, opt: &AdaptOptions) -> Option<f64> {
    let e = ekey(a, b);
    match w.vclass[a as usize] {
        VClass::Corner => return None,
        VClass::Chain(c) => {
            if w.chains[c as usize].frozen || w.chain_of.get(&e) != Some(&c) {
                return None;
            }
            match w.vclass[b as usize] {
                VClass::Chain(cb) if cb == c => {}
                VClass::Corner => {}
                _ => return None,
            }
        }
        VClass::Interior => {}
    }
    let ts = snap.edges[snap.index(&e).unwrap()].1;
    let shared = if ts[1] >= 0 { 2 } else { 1 };
    let (na, nb) = (neighbours(w, snap, a), neighbours(w, snap, b));
    let common = na.iter().filter(|x| nb.binary_search(x).is_ok()).count();
    if common != shared {
        return None;
    }
    let mut worst = f64::INFINITY;
    for &t in snap.star(a) {
        let tr = w.tris[t as usize];
        if tr.contains(&b) {
            continue;
        }
        let old_q = tq(w, tr);
        let nt = [
            if tr[0] == a { b } else { tr[0] },
            if tr[1] == a { b } else { tr[1] },
            if tr[2] == a { b } else { tr[2] },
        ];
        if !(tri_area(&w.pts, &nt) > 0.0) {
            return None;
        }
        let q = tq(w, nt);
        if q < opt.collapse_qmin.min(old_q) {
            return None;
        }
        for &x in &nt {
            if x != b && elen(w, b, x) > opt.lmax {
                return None;
            }
        }
        worst = worst.min(q);
    }
    Some(worst)
}

/// Collapse every edge shorter than `opt.lmin` that can be collapsed. Returns the number done.
pub fn collapse_pass(w: &mut Work2, opt: &AdaptOptions) -> usize {
    let snap = w.snapshot();
    let mut cand: Vec<(f64, usize)> = Vec::new();
    for (i, (e, _)) in snap.edges.iter().enumerate() {
        let l = elen(w, e[0], e[1]);
        if l < opt.lmin {
            cand.push((l, i));
        }
    }
    cand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
    let mut lock = vec![false; w.tris.len()];
    let mut n = 0;
    for (_, i) in cand {
        let (e, _) = snap.edges[i];
        let (a, b) = (e[0], e[1]);
        if snap.star(a).iter().chain(snap.star(b)).any(|&t| lock[t as usize]) {
            continue;
        }
        let ra = collapse_eval(w, &snap, a, b, opt);
        let rb = collapse_eval(w, &snap, b, a, opt);
        let (from, to) = match (ra, rb) {
            (None, None) => continue,
            (Some(_), None) => (a, b),
            (None, Some(_)) => (b, a),
            (Some(qa), Some(qb)) => {
                if qa >= qb {
                    (a, b)
                } else {
                    (b, a)
                }
            }
        };
        for &x in &neighbours(w, &snap, from) {
            let k = ekey(from, x);
            // the chain entry is removed whether or not it moves on to the edge (to, x)
            if let Some(c) = w.chain_of.remove(&k) {
                if x != to {
                    w.chain_of.entry(ekey(to, x)).or_insert(c);
                }
            }
        }
        for &t in snap.star(from) {
            let tr = w.tris[t as usize];
            if tr.contains(&to) {
                w.talive[t as usize] = false;
            } else {
                w.tris[t as usize] = [
                    if tr[0] == from { to } else { tr[0] },
                    if tr[1] == from { to } else { tr[1] },
                    if tr[2] == from { to } else { tr[2] },
                ];
                w.touched[t as usize] = true;
            }
        }
        for &t in snap.star(from).iter().chain(snap.star(to)) {
            lock[t as usize] = true;
        }
        n += 1;
    }
    n
}

/// Flip interior edges when that improves the worst quality of the two triangles.
pub fn flip_pass(w: &mut Work2) -> usize {
    let snap = w.snapshot();
    let mut lock = vec![false; w.tris.len()];
    let active = std::mem::replace(&mut w.touched, vec![false; w.tris.len()]);
    let mut n = 0;
    for (e, ts) in &snap.edges {
        if ts[1] < 0 || !(active[ts[0] as usize] || active[ts[1] as usize]) || w.chain_of.contains_key(e) {
            continue;
        }
        let (t1, t2) = (ts[0] as usize, ts[1] as usize);
        if lock[t1] || lock[t2] {
            continue;
        }
        let tr1 = w.tris[t1];
        let k = (0..3).find(|&k| ekey(tr1[(k + 1) % 3], tr1[(k + 2) % 3]) == *e).unwrap();
        let (c, p, q) = (tr1[k], tr1[(k + 1) % 3], tr1[(k + 2) % 3]);
        let tr2 = w.tris[t2];
        let k2 = match (0..3).find(|&k| tr2[(k + 1) % 3] == q && tr2[(k + 2) % 3] == p) {
            Some(k) => k,
            None => continue,
        };
        let d = tr2[k2];
        if snap.index(&ekey(c, d)).is_some() {
            continue;
        }
        let (n1, n2) = ([c, p, d], [c, d, q]);
        if !(tri_area(&w.pts, &n1) > 0.0 && tri_area(&w.pts, &n2) > 0.0) {
            continue;
        }
        let old = tq(w, tr1).min(tq(w, tr2));
        let new = tq(w, n1).min(tq(w, n2));
        if new > old * (1.0 + 1e-6) + 1e-12 {
            w.tris[t1] = n1;
            w.tris[t2] = n2;
            lock[t1] = true;
            lock[t2] = true;
            w.touched[t1] = true;
            w.touched[t2] = true;
            n += 1;
        }
    }
    n
}

/// One Gauss-Seidel sweep of length-equalising smoothing over the movable vertices. A vertex moves
/// to the relaxed average of the points that would give each incident edge unit metric length;
/// chain vertices are projected back onto their chain. Moves that create an inverted triangle or
/// lower the worst incident quality are rejected. Returns the number of vertices moved.
pub fn smooth_pass(w: &mut Work2, field: &dyn MetricField<2>, geom: &dyn Geometry2, opt: &AdaptOptions) -> usize {
    let snap = w.snapshot();
    let mut moved = 0;
    for v in 0..w.pts.len() as u32 {
        let chain = match w.vclass[v as usize] {
            VClass::Corner => continue,
            VClass::Chain(c) => {
                if w.chains[c as usize].frozen {
                    continue;
                }
                Some(c)
            }
            VClass::Interior => None,
        };
        let star = snap.star(v);
        if star.is_empty() {
            continue;
        }
        let xv = w.pts[v as usize];
        let nb = neighbours(w, &snap, v);
        let (mut sx, mut sy, mut cnt) = (0.0, 0.0, 0.0);
        for &n in &nb {
            let xn = w.pts[n as usize];
            let l = elen(w, v, n);
            if l < 1e-12 {
                continue;
            }
            sx += xn[0] + (xv[0] - xn[0]) / l;
            sy += xn[1] + (xv[1] - xn[1]) / l;
            cnt += 1.0;
        }
        if cnt == 0.0 {
            continue;
        }
        let r = opt.smooth_relax;
        let mut cands: Vec<[f64; 2]> = vec![[xv[0] + r * (sx / cnt - xv[0]), xv[1] + r * (sy / cnt - xv[1])]];
        // evaluate a candidate: (min quality over the star after the move, metric there), None if invalid
        let old_min = star.iter().map(|&t| tq(w, w.tris[t as usize])).fold(f64::INFINITY, f64::min);
        if old_min < 0.4 {
            // poorly shaped star: also try the centroid of the ring (pulls caps and slivers inward)
            let (mut cx, mut cy) = (0.0, 0.0);
            for &n in &nb {
                cx += w.pts[n as usize][0] / nb.len() as f64;
                cy += w.pts[n as usize][1] / nb.len() as f64;
            }
            cands.push([xv[0] + r * (cx - xv[0]), xv[1] + r * (cy - xv[1])]);
        }
        let mut best: Option<(f64, [f64; 2], Sym<2>)> = None;
        for mut np in cands {
            if let Some(c) = chain {
                np = geom.project(np, c);
            }
            let new_met = field.at(&np);
            let mut new_min = f64::INFINITY;
            let mut ok = true;
            for &t in star {
                let tr = w.tris[t as usize];
                let g = |i: usize| if tr[i] == v { np } else { w.pts[tr[i] as usize] };
                let m = |i: usize| if tr[i] == v { &new_met } else { &w.met[tr[i] as usize] };
                let (p0, p1, p2) = (g(0), g(1), g(2));
                let area = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]));
                if !(area > 0.0) {
                    ok = false;
                    break;
                }
                new_min = new_min.min(tri_quality([p0, p1, p2], [m(0), m(1), m(2)]));
            }
            if ok && new_min >= old_min - 1e-12 && best.as_ref().is_none_or(|b| new_min > b.0) {
                best = Some((new_min, np, new_met));
            }
        }
        if let Some((_, np, new_met)) = best {
            w.pts[v as usize] = np;
            w.met[v as usize] = new_met;
            for &t in star {
                w.touched[t as usize] = true;
            }
            moved += 1;
        }
    }
    moved
}

/// Number of frozen edges longer than `lmax` (they cannot be split, so the input must resolve them).
pub fn unresolved_frozen(w: &Work2, lmax: f64) -> usize {
    w.snapshot()
        .edges
        .iter()
        .filter(|(e, ts)| {
            ts[0] >= 0 && w.chain_of.get(e).is_some_and(|&c| w.chains[c as usize].frozen) && elen(w, e[0], e[1]) > lmax
        })
        .count()
}

/// Sorted list of the edges (helper for statistics).
pub fn all_edges(w: &Work2) -> Vec<Edge> {
    w.snapshot().edges.iter().map(|x| x.0).collect()
}
