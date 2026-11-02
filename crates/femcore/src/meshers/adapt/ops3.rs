//! Local mesh operators for 3D metric adaptation: edge split, collapse of interior vertices, 3-2 and 2-3 flips, and length-equalising
//! smoothing (the 3D counterparts of `ops`).
//!
//! Every pass takes a snapshot of the incidence, orders its candidates deterministically (by decreasing or increasing metric length,
//! ties by edge) and applies them greedily. The tetrahedra an operator touches are *locked* for the rest of the pass, so two operators
//! in one pass never overlap. All operators are quality gated: the result must have positive volume everywhere and must not create
//! tetrahedra worse than the thresholds in `AdaptOptions3`; surface faces, creases and corners are never violated: a flip never removes
//! a surface face or a crease edge, a split puts the new vertex on the surface or crease it lies on (through the `Geometry3` oracle) and
//! a smoothing step slides a surface or crease vertex over its patch or curve.
//!
//! **Limits.** Collapse removes interior vertices only, so a surface or crease that is finer than the metric asks is not coarsened
//! (refinement towards a metric is complete; coarsening of the surface is not). Flips are 3-2 and 2-3; edge removal with more tets
//! around the edge, and face-edge flips on the surface, are not implemented. Serial.

use super::field::MetricField;
use super::geom3::Geometry3;
use super::metric::{edge_length, Sym};
use super::topo::{ekey, Edge, FxMap};
use super::topo3::{fkey, tet_faces, tet_volume, VClass3, Work3};
use crate::meshers::quality::simplex_quality;

/// Thresholds of the operators.
#[derive(Clone, Debug)]
pub struct AdaptOptions3 {
    /// Edges longer than this (in the metric) are split. Default `sqrt 2`.
    pub lmax: f64,
    /// Interior edges shorter than this are collapsed. Default `1 / sqrt 2`.
    pub lmin: f64,
    pub max_iter: usize,
    pub smooth_sweeps: usize,
    pub smooth_relax: f64,
    /// A collapse may not create a tetrahedron below `min(collapse_qmin, quality it replaces)`.
    pub collapse_qmin: f64,
    /// A split may not create a tetrahedron below this quality.
    pub split_qmin: f64,
    /// During the first `sample_iters` iterations edges are also measured through the field value at their midpoint.
    pub sample_iters: usize,
    /// Runaway guard: stop when an iteration ends with more tetrahedra than this (a split cascade is stopped at twenty times it).
    pub max_cells: Option<usize>,
    /// Stop when an iteration's splits plus collapses are below this fraction of the tetrahedra.
    pub converge_frac: f64,
}

impl Default for AdaptOptions3 {
    fn default() -> Self {
        AdaptOptions3 {
            lmax: std::f64::consts::SQRT_2,
            lmin: std::f64::consts::FRAC_1_SQRT_2,
            max_iter: 30,
            smooth_sweeps: 3,
            smooth_relax: 0.7,
            collapse_qmin: 0.1,
            split_qmin: 1e-3,
            sample_iters: 2,
            max_cells: None,
            converge_frac: 0.005,
        }
    }
}

/// Incidence of the alive tetrahedra.
pub struct Snap3 {
    pub edge_tets: FxMap<Edge, Vec<u32>>,
    pub vert_tets: Vec<Vec<u32>>,
}

pub fn snapshot(w: &Work3) -> Snap3 {
    let mut edge_tets: FxMap<Edge, Vec<u32>> = FxMap::default();
    let mut vert_tets = vec![vec![]; w.pts.len()];
    for (t, tt) in w.tets.iter().enumerate() {
        if !w.talive[t] {
            continue;
        }
        for i in 0..4 {
            vert_tets[tt[i] as usize].push(t as u32);
            for j in i + 1..4 {
                edge_tets.entry(ekey(tt[i], tt[j])).or_default().push(t as u32);
            }
        }
    }
    Snap3 { edge_tets, vert_tets }
}

/// Metric length of the edge `(a, b)`.
pub fn elen(w: &Work3, a: u32, b: u32) -> f64 {
    let (pa, pb) = (w.pts[a as usize], w.pts[b as usize]);
    edge_length(&w.met[a as usize], &w.met[b as usize], &[pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]])
}

/// Mean-ratio quality in the metric of the tetrahedron `t` with the points and metrics of `w`.
pub fn tet_q(w: &Work3, t: &[u32; 4]) -> f64 {
    let p: [[f64; 3]; 4] = std::array::from_fn(|i| w.pts[t[i] as usize]);
    let m: [Sym<3>; 4] = std::array::from_fn(|i| w.met[t[i] as usize]);
    simplex_quality::<3>(&p, Some(&m))
}

/// All edges of the alive tetrahedra, sorted.
pub fn all_edges_of(w: &Work3) -> Vec<Edge> {
    all_edges(&snapshot(w))
}

fn all_edges(snap: &Snap3) -> Vec<Edge> {
    let mut e: Vec<Edge> = snap.edge_tets.keys().copied().collect();
    e.sort_unstable();
    e
}

fn replace(t: &[u32; 4], from: u32, to: u32) -> [u32; 4] {
    std::array::from_fn(|i| if t[i] == from { to } else { t[i] })
}

fn locked(lock: &[bool], ts: &[u32]) -> bool {
    ts.iter().any(|&t| lock[t as usize])
}

/// The surface patch an edge lies on, if one of the faces on it is a surface face.
fn patch_of_edge(w: &Work3, snap: &Snap3, e: &Edge) -> Option<u32> {
    for &t in &snap.edge_tets[e] {
        for (f, _) in tet_faces(&w.tets[t as usize]) {
            if f.contains(&e[0]) && f.contains(&e[1]) {
                if let Some(&p) = w.patch_of.get(&f) {
                    return Some(p);
                }
            }
        }
    }
    None
}

/// Number of edges on a frozen curve or patch that are longer than `lmax`.
pub fn unresolved_frozen(w: &Work3, lmax: f64) -> usize {
    let snap = snapshot(w);
    all_edges(&snap).iter().filter(|e| edge_frozen(w, &snap, e) && elen(w, e[0], e[1]) > lmax).count()
}

fn edge_frozen(w: &Work3, snap: &Snap3, e: &Edge) -> bool {
    if let Some(&c) = w.curve_of.get(e) {
        return w.curves[c as usize].frozen;
    }
    patch_of_edge(w, snap, e).is_some_and(|p| w.patches[p as usize].frozen)
}

/// Split every edge longer than `opt.lmax`. Returns the number of splits.
pub fn split_pass(w: &mut Work3, field: &dyn MetricField<3>, geom: &dyn Geometry3, opt: &AdaptOptions3, sample_mid: bool) -> usize {
    let snap = snapshot(w);
    let mut cand: Vec<(f64, Edge)> = vec![];
    for e in all_edges(&snap) {
        let mut l = elen(w, e[0], e[1]);
        if sample_mid && l > 0.3 && l <= opt.lmax {
            let (pa, pb) = (w.pts[e[0] as usize], w.pts[e[1] as usize]);
            let m = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            let mm = field.at(&m);
            let h = [0.5 * (pb[0] - pa[0]), 0.5 * (pb[1] - pa[1]), 0.5 * (pb[2] - pa[2])];
            l = edge_length(&w.met[e[0] as usize], &mm, &h) + edge_length(&mm, &w.met[e[1] as usize], &h);
        }
        if l > opt.lmax && !edge_frozen(w, &snap, &e) {
            cand.push((l, e));
        }
    }
    cand.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
    let mut lock = vec![false; w.tets.len()];
    let mut n = 0;
    for (_, e) in cand {
        let ts = snap.edge_tets[&e].clone();
        if locked(&lock, &ts) {
            continue;
        }
        let [a, b] = e;
        let curve = w.curve_of.get(&e).copied();
        let patch = if curve.is_none() { patch_of_edge(w, &snap, &e) } else { None };
        let (pa, pb) = (w.pts[a as usize], w.pts[b as usize]);
        let mut m = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
        if let Some(c) = curve {
            m = geom.project_crease(m, c);
        } else if let Some(p) = patch {
            m = geom.project_surface(m, p);
        }
        let met_m = field.at(&m);
        // validity of the two halves of every tetrahedron on the edge
        let mut ok = true;
        'chk: for &t in &ts {
            let tt = w.tets[t as usize];
            let old = tet_volume(&w.pts, &tt);
            for from in [b, a] {
                let child = replace(&tt, from, w.pts.len() as u32);
                let p: [[f64; 3]; 4] = std::array::from_fn(|i| if child[i] == w.pts.len() as u32 { m } else { w.pts[child[i] as usize] });
                let vol = crate::meshers::mesh::simplex_signed_measure::<3>(&p);
                let mm: [Sym<3>; 4] = std::array::from_fn(|i| if child[i] == w.pts.len() as u32 { met_m } else { w.met[child[i] as usize] });
                if !(vol > 1e-6 * old) || simplex_quality::<3>(&p, Some(&mm)) < opt.split_qmin {
                    ok = false;
                    break 'chk;
                }
            }
        }
        if !ok {
            continue;
        }
        let mid = w.pts.len() as u32;
        // surface faces on the edge, before the tetrahedra change
        let mut sf: Vec<[u32; 3]> = vec![];
        for &t in &ts {
            for (f, _) in tet_faces(&w.tets[t as usize]) {
                if f.contains(&a) && f.contains(&b) && w.patch_of.contains_key(&f) && !sf.contains(&f) {
                    sf.push(f);
                }
            }
        }
        w.pts.push(m);
        w.met.push(met_m);
        w.vclass.push(match (curve, patch) {
            (Some(c), _) => VClass3::Crease(c),
            (None, Some(p)) => VClass3::Surface(p),
            _ => VClass3::Interior,
        });
        for &t in &ts {
            let tt = w.tets[t as usize];
            let label = w.tlabel[t as usize];
            w.tets[t as usize] = replace(&tt, b, mid);
            w.tets.push(replace(&tt, a, mid));
            w.tlabel.push(label);
            w.talive.push(true);
            lock[t as usize] = true;
            lock.push(true);
        }
        for f in sf {
            let c = *f.iter().find(|&&x| x != a && x != b).unwrap();
            let p = w.patch_of.remove(&f).unwrap();
            w.patch_of.insert(fkey(a, mid, c), p);
            w.patch_of.insert(fkey(mid, b, c), p);
        }
        if let Some(c) = w.curve_of.remove(&e) {
            w.curve_of.insert(ekey(a, mid), c);
            w.curve_of.insert(ekey(mid, b), c);
        }
        n += 1;
    }
    n
}

fn neighbours(w: &Work3, snap: &Snap3, v: u32) -> Vec<u32> {
    let mut n: Vec<u32> = vec![];
    for &t in &snap.vert_tets[v as usize] {
        n.extend(w.tets[t as usize].iter().copied().filter(|&x| x != v));
    }
    n.sort_unstable();
    n.dedup();
    n
}

/// Can interior vertex `a` be merged into `b`? Returns the worst new quality if so.
fn collapse_eval(w: &Work3, snap: &Snap3, a: u32, b: u32, opt: &AdaptOptions3) -> Option<f64> {
    if w.vclass[a as usize] != VClass3::Interior {
        return None;
    }
    let e = ekey(a, b);
    // link condition: the common neighbours of a and b are exactly the vertices opposite the edge in its tetrahedra
    let mut link: Vec<u32> = vec![];
    for &t in &snap.edge_tets[&e] {
        link.extend(w.tets[t as usize].iter().copied().filter(|&x| x != a && x != b));
    }
    link.sort_unstable();
    link.dedup();
    let (na, nb) = (neighbours(w, snap, a), neighbours(w, snap, b));
    let common: Vec<u32> = na.iter().copied().filter(|x| nb.binary_search(x).is_ok()).collect();
    if common != link {
        return None;
    }
    // The vertex test above is not enough in 3D (the links must also agree on edges and triangles): no face at b may end up with more
    // than two tetrahedra, counting the ones that stay and the ones that move from a to b.
    let mut faces_at_b: FxMap<[u32; 3], u8> = FxMap::default();
    for &t in &snap.vert_tets[b as usize] {
        let tt = w.tets[t as usize];
        if !tt.contains(&a) {
            for (f, _) in tet_faces(&tt) {
                if f.contains(&b) {
                    *faces_at_b.entry(f).or_default() += 1;
                }
            }
        }
    }
    let mut worst = f64::INFINITY;
    for &t in &snap.vert_tets[a as usize] {
        let tt = w.tets[t as usize];
        if tt.contains(&b) {
            continue;
        }
        let nt = replace(&tt, a, b);
        for (f, _) in tet_faces(&nt) {
            if f.contains(&b) {
                let c = faces_at_b.entry(f).or_default();
                *c += 1;
                if *c > 2 {
                    return None;
                }
            }
        }
        if !(tet_volume(&w.pts, &nt) > 0.0) {
            return None;
        }
        let q = tet_q(w, &nt);
        if q < opt.collapse_qmin.min(tet_q(w, &tt)) {
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

/// Collapse every interior edge shorter than `opt.lmin` that can be collapsed (the interior end moves onto the other). Returns the number.
pub fn collapse_pass(w: &mut Work3, opt: &AdaptOptions3) -> usize {
    let snap = snapshot(w);
    let mut cand: Vec<(f64, Edge)> = vec![];
    for e in all_edges(&snap) {
        let l = elen(w, e[0], e[1]);
        if l < opt.lmin && (w.vclass[e[0] as usize] == VClass3::Interior || w.vclass[e[1] as usize] == VClass3::Interior) {
            cand.push((l, e));
        }
    }
    cand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1)));
    let mut lock = vec![false; w.tets.len()];
    let mut n = 0;
    for (_, e) in cand {
        let mut best: Option<(f64, u32, u32)> = None;
        for (a, b) in [(e[0], e[1]), (e[1], e[0])] {
            if locked(&lock, &snap.vert_tets[a as usize]) || locked(&lock, &snap.vert_tets[b as usize]) {
                continue;
            }
            if let Some(q) = collapse_eval(w, &snap, a, b, opt) {
                if best.is_none_or(|x| q > x.0) {
                    best = Some((q, a, b));
                }
            }
        }
        let Some((_, a, b)) = best else { continue };
        for &t in &snap.vert_tets[a as usize] {
            let tt = w.tets[t as usize];
            if tt.contains(&b) {
                w.talive[t as usize] = false;
            } else {
                w.tets[t as usize] = replace(&tt, a, b);
            }
            lock[t as usize] = true;
        }
        for &t in &snap.vert_tets[b as usize] {
            lock[t as usize] = true;
        }
        n += 1;
    }
    n
}

/// Quality-improving 3-2 and 2-3 flips of interior edges and faces. Returns the number of flips.
pub fn flip_pass(w: &mut Work3) -> usize {
    let snap = snapshot(w);
    let mut lock = vec![false; w.tets.len()];
    let mut n = 0;
    // 3-2: an interior edge with three tetrahedra around it
    for e in all_edges(&snap) {
        let ts = &snap.edge_tets[&e];
        if ts.len() != 3 || locked(&lock, ts) || w.curve_of.contains_key(&e) {
            continue;
        }
        let label = w.tlabel[ts[0] as usize];
        if ts.iter().any(|&t| w.tlabel[t as usize] != label) {
            continue;
        }
        if ts.iter().any(|&t| tet_faces(&w.tets[t as usize]).iter().any(|(f, _)| f.contains(&e[0]) && f.contains(&e[1]) && w.patch_of.contains_key(f))) {
            continue;
        }
        let mut ring: Vec<u32> = ts.iter().flat_map(|&t| w.tets[t as usize].iter().copied().filter(|&x| x != e[0] && x != e[1])).collect();
        ring.sort_unstable();
        ring.dedup();
        if ring.len() != 3 {
            continue;
        }
        // the two new tetrahedra fill the bipyramid only if the ends of the edge lie on opposite sides of the ring's plane
        let sides = [e[0], e[1]].map(|apex| tet_volume(&w.pts, &[ring[0], ring[1], ring[2], apex]));
        if !(sides[0] * sides[1] < 0.0) {
            continue;
        }
        // and the ring triangle must not already be a face of the mesh
        if snap.vert_tets[ring[0] as usize].iter().any(|&t| ring.iter().all(|r| w.tets[t as usize].contains(r))) {
            continue;
        }
        let new: Vec<[u32; 4]> = [e[0], e[1]].iter().map(|&apex| orient(w, [ring[0], ring[1], ring[2], apex])).collect();
        let old_min = ts.iter().map(|&t| tet_q(w, &w.tets[t as usize])).fold(f64::INFINITY, f64::min);
        let new_min = new.iter().map(|t| tet_q(w, t)).fold(f64::INFINITY, f64::min);
        if new_min <= old_min {
            continue;
        }
        for &t in ts {
            w.talive[t as usize] = false;
            lock[t as usize] = true;
        }
        for t in new {
            w.tets.push(t);
            w.tlabel.push(label);
            w.talive.push(true);
            lock.push(true);
        }
        n += 1;
    }
    // 2-3: an interior face shared by two tetrahedra of one label
    let mut faces: Vec<([u32; 3], [u32; 2])> = vec![];
    {
        let mut by_face: FxMap<[u32; 3], Vec<u32>> = FxMap::default();
        for (t, tt) in w.tets.iter().enumerate() {
            if t < lock.len() && w.talive[t] && !lock[t] {
                for (f, _) in tet_faces(tt) {
                    by_face.entry(f).or_default().push(t as u32);
                }
            }
        }
        for (f, ts) in by_face {
            if ts.len() == 2 && !w.patch_of.contains_key(&f) && w.tlabel[ts[0] as usize] == w.tlabel[ts[1] as usize] {
                faces.push((f, [ts[0], ts[1]]));
            }
        }
        faces.sort_unstable();
    }
    for (f, [t1, t2]) in faces {
        if lock[t1 as usize] || lock[t2 as usize] || !w.talive[t1 as usize] || !w.talive[t2 as usize] {
            continue;
        }
        let opp = |t: u32| *w.tets[t as usize].iter().find(|x| !f.contains(x)).unwrap();
        let (p, q) = (opp(t1), opp(t2));
        if snap.edge_tets.contains_key(&ekey(p, q)) {
            continue;
        }
        let raw = [[p, q, f[0], f[1]], [p, q, f[1], f[2]], [p, q, f[2], f[0]]];
        let vols: Vec<f64> = raw.iter().map(|t| tet_volume(&w.pts, t)).collect();
        let sign = vols[0].signum();
        if vols.iter().any(|v| v.signum() != sign || v.abs() < 1e-14) {
            continue;
        }
        let new: Vec<[u32; 4]> = raw.iter().map(|&t| if sign > 0.0 { t } else { [t[1], t[0], t[2], t[3]] }).collect();
        let old_min = tet_q(w, &w.tets[t1 as usize]).min(tet_q(w, &w.tets[t2 as usize]));
        let new_min = new.iter().map(|t| tet_q(w, t)).fold(f64::INFINITY, f64::min);
        if new_min <= old_min {
            continue;
        }
        let label = w.tlabel[t1 as usize];
        for t in [t1, t2] {
            w.talive[t as usize] = false;
            lock[t as usize] = true;
        }
        for t in new {
            w.tets.push(t);
            w.tlabel.push(label);
            w.talive.push(true);
            lock.push(true);
        }
        n += 1;
    }
    n
}


/// The tetrahedron with the vertex order that gives it positive volume.
fn orient(w: &Work3, t: [u32; 4]) -> [u32; 4] {
    if tet_volume(&w.pts, &t) < 0.0 { [t[1], t[0], t[2], t[3]] } else { t }
}

/// One sweep of length-equalising smoothing: every free vertex moves towards the position where its edges would have unit metric length
/// (projected onto its patch or curve when it lies on one), if that does not worsen the worst tetrahedron of its star.
pub fn smooth_pass(w: &mut Work3, field: &dyn MetricField<3>, geom: &dyn Geometry3, opt: &AdaptOptions3) -> usize {
    let snap = snapshot(w);
    let mut moved = 0;
    for v in 0..w.pts.len() as u32 {
        let class = w.vclass[v as usize];
        match class {
            VClass3::Corner => continue,
            VClass3::Surface(p) if w.patches[p as usize].frozen => continue,
            VClass3::Crease(c) if w.curves[c as usize].frozen => continue,
            _ => {}
        }
        let star = &snap.vert_tets[v as usize];
        if star.is_empty() {
            continue;
        }
        // the neighbours that move it: along the curve, over the patch, or through the volume
        let mut nb: Vec<u32> = vec![];
        match class {
            VClass3::Crease(_) => {
                for e in w.curve_of.keys() {
                    if e[0] == v {
                        nb.push(e[1]);
                    } else if e[1] == v {
                        nb.push(e[0]);
                    }
                }
            }
            VClass3::Surface(_) => {
                for &t in star {
                    for (f, _) in tet_faces(&w.tets[t as usize]) {
                        if f.contains(&v) && w.patch_of.contains_key(&f) {
                            nb.extend(f.iter().copied().filter(|&x| x != v));
                        }
                    }
                }
                nb.sort_unstable();
                nb.dedup();
            }
            _ => nb = neighbours(w, &snap, v),
        }
        let xv = w.pts[v as usize];
        let (mut s, mut cnt) = ([0.0; 3], 0.0);
        for &n in &nb {
            let xn = w.pts[n as usize];
            let l = elen(w, v, n);
            if l < 1e-12 {
                continue;
            }
            for k in 0..3 {
                s[k] += xn[k] + (xv[k] - xn[k]) / l;
            }
            cnt += 1.0;
        }
        if cnt == 0.0 {
            continue;
        }
        let r = opt.smooth_relax;
        let mut np: [f64; 3] = std::array::from_fn(|k| xv[k] + r * (s[k] / cnt - xv[k]));
        match class {
            VClass3::Crease(c) => np = geom.project_crease(np, c),
            VClass3::Surface(p) => np = geom.project_surface(np, p),
            _ => {}
        }
        let new_met = field.at(&np);
        let old_min = star.iter().map(|&t| tet_q(w, &w.tets[t as usize])).fold(f64::INFINITY, f64::min);
        let mut new_min = f64::INFINITY;
        let mut ok = true;
        for &t in star {
            let tt = w.tets[t as usize];
            let p: [[f64; 3]; 4] = std::array::from_fn(|i| if tt[i] == v { np } else { w.pts[tt[i] as usize] });
            let m: [Sym<3>; 4] = std::array::from_fn(|i| if tt[i] == v { new_met } else { w.met[tt[i] as usize] });
            if !(crate::meshers::mesh::simplex_signed_measure::<3>(&p) > 0.0) {
                ok = false;
                break;
            }
            new_min = new_min.min(simplex_quality::<3>(&p, Some(&m)));
        }
        if ok && new_min >= old_min - 1e-12 {
            w.pts[v as usize] = np;
            w.met[v as usize] = new_met;
            moved += 1;
        }
    }
    moved
}
