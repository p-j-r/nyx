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
//! Edge removal replaces an interior edge with 4 to 7 tets around it by the best triangulation of its ring (`edge_removal_pass`), and a
//! surface flip changes the diagonal of two coplanar boundary faces (`surface_flip_pass`).
//!
//! **Limits.** Edge removal stops at 7 tets around the edge; a surface flip needs a flat boundary patch (not an interface) with two
//! tets on the edge, so curved surfaces and interfaces keep their triangulation. A surface collapse refuses to turn a surface face by more than about 25 degrees, so a curved surface
//! coarsens only where it is nearly flat at the scale of the edge. Serial.

use super::field::MetricField;
use super::geom3::Geometry3;
use super::metric::{edge_length, Sym};
use super::topo::{ekey, Edge, FxMap, Sig};
use super::topo3::{fkey, tet_faces, tet_volume, VClass3, Work3};
use crate::mesh::quality::simplex_quality;

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
    /// Also try edge removal and surface flips after the 3-2 / 2-3 flips.
    pub edge_removal: bool,
    /// The collapse and flip sweeps of an iteration repeat until one changes at most this fraction of the tetrahedra (0: until one changes
    /// nothing). Every sweep rebuilds the incidence of the whole mesh, so on a large mesh the last sweeps that change a handful of
    /// tetrahedra cost as much as the first.
    pub sweep_frac: f64,
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
            edge_removal: true,
            sweep_frac: 0.0,
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
                let vol = crate::adapt::work::simplex_signed_measure::<3>(&p);
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

/// Vertices joined to `v` by an edge of a surface face.
fn surface_nbrs(w: &Work3, snap: &Snap3, v: u32) -> Vec<u32> {
    let mut n: Vec<u32> = vec![];
    for &t in &snap.vert_tets[v as usize] {
        for (f, _) in tet_faces(&w.tets[t as usize]) {
            if f.contains(&v) && w.patch_of.contains_key(&f) {
                n.extend(f.iter().copied().filter(|&x| x != v));
            }
        }
    }
    n.sort_unstable();
    n.dedup();
    n
}

/// The surface faces on the edge `e`, as `(face, patch)`.
fn surface_faces_on(w: &Work3, snap: &Snap3, e: &Edge) -> Vec<([u32; 3], u32)> {
    let mut out: Vec<([u32; 3], u32)> = vec![];
    for &t in &snap.edge_tets[e] {
        for (f, _) in tet_faces(&w.tets[t as usize]) {
            if f.contains(&e[0]) && f.contains(&e[1]) {
                if let Some(&p) = w.patch_of.get(&f) {
                    if !out.iter().any(|x| x.0 == f) {
                        out.push((f, p));
                    }
                }
            }
        }
    }
    out
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Unit normal of the face `f` of tetrahedron `t` pointing out of `t`.
fn outward_normal(pts: &[[f64; 3]], t: &[u32; 4], f: &[u32; 3]) -> [f64; 3] {
    let opp = *t.iter().find(|x| !f.contains(x)).unwrap();
    let (a, b, c) = (pts[f[0] as usize], pts[f[1] as usize], pts[f[2] as usize]);
    let (u, v) = (sub3(b, a), sub3(c, a));
    let mut n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-300);
    let d = sub3(pts[opp as usize], a);
    if n[0] * d[0] + n[1] * d[1] + n[2] * d[2] > 0.0 {
        n = [-n[0], -n[1], -n[2]];
    }
    [n[0] / l, n[1] / l, n[2] / l]
}

fn has_sheet_patches(w: &Work3) -> bool {
    w.patches.iter().any(|p| matches!(p.sig, Sig::Sheet(_)))
}

/// Can vertex `a` be merged into `b`? An interior vertex moves onto anything; a surface vertex only along an edge of its own patch;
/// a crease vertex only along its own curve; a corner never. Returns the worst new quality if so.
fn collapse_eval(w: &Work3, snap: &Snap3, a: u32, b: u32, opt: &AdaptOptions3, has_sheets: bool) -> Option<f64> {
    let e = ekey(a, b);
    match w.vclass[a as usize] {
        VClass3::Corner => return None,
        VClass3::Surface(p) => {
            if w.patches[p as usize].frozen || !surface_faces_on(w, snap, &e).iter().any(|x| x.1 == p) {
                return None;
            }
        }
        VClass3::Crease(c) => {
            if w.curves[c as usize].frozen || w.curve_of.get(&e) != Some(&c) {
                return None;
            }
            match w.vclass[b as usize] {
                VClass3::Crease(cb) if cb == c => {}
                VClass3::Corner => {}
                _ => return None,
            }
        }
        VClass3::Interior => {}
    }
    if w.vclass[a as usize] != VClass3::Interior {
        // the same link condition on the surface: the common surface neighbours are the apexes of the surface faces on the edge
        let (sa, sb) = (surface_nbrs(w, snap, a), surface_nbrs(w, snap, b));
        let common: Vec<u32> = sa.iter().copied().filter(|x| sb.binary_search(x).is_ok()).collect();
        let mut opp: Vec<u32> = surface_faces_on(w, snap, &e).iter().map(|(f, _)| *f.iter().find(|&&x| x != a && x != b).unwrap()).collect();
        opp.sort_unstable();
        if has_sheets {
            opp.dedup();
        }
        if common != opp {
            return None;
        }
        if has_sheets {
            if let VClass3::Crease(_) = w.vclass[a as usize] {
                for k in w.curve_of.keys() {
                    if k.contains(&a) && !k.contains(&b) {
                        let o = if k[0] == a { k[1] } else { k[0] };
                        if w.curve_of.contains_key(&ekey(b, o)) {
                            return None;
                        }
                    }
                }
            }
        }
    }
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

    // Track post-collapse surface patches on faces containing `b` so sheet/interface/boundary topology is strictly preserved.
    let mut post_patch_at_b: FxMap<[u32; 3], u32> = FxMap::default();
    for &t in &snap.vert_tets[b as usize] {
        let tt = w.tets[t as usize];
        for (f, _) in tet_faces(&tt) {
            if f.contains(&b) && !f.contains(&a) {
                if let Some(&p) = w.patch_of.get(&f) {
                    post_patch_at_b.insert(f, p);
                }
            }
        }
    }
    if w.vclass[a as usize] != VClass3::Interior {
        let mut moved_sf: Vec<([u32; 3], u32)> = Vec::new();
        for &t in &snap.vert_tets[a as usize] {
            let tt = w.tets[t as usize];
            if !tt.contains(&b) {
                for (f, _) in tet_faces(&tt) {
                    if f.contains(&a) {
                        if let Some(&p) = w.patch_of.get(&f) {
                            if !moved_sf.iter().any(|x| x.0 == f) {
                                moved_sf.push((f, p));
                            }
                        }
                    }
                }
            }
        }
        for (f, p) in moved_sf {
            let g: Vec<u32> = f.iter().map(|&x| if x == a { b } else { x }).collect();
            let f_new = fkey(g[0], g[1], g[2]);
            if w.patch_of.contains_key(&f_new) || post_patch_at_b.insert(f_new, p).is_some() {
                return None;
            }
        }
    }

    // The vertex test above is not enough in 3D (the links must also agree on edges and triangles): no face at b may end up with more
    // than two tetrahedra, two tetrahedra sharing a face must lie on strictly opposite sides of it, and boundary/sheet/interface counts must hold.
    let mut faces_at_b: FxMap<[u32; 3], (u8, f64, i32)> = FxMap::default();
    let mut reg_face = |f: [u32; 3], opp: u32, lab: i32| -> bool {
        let side = tet_volume(&w.pts, &[f[0], f[1], f[2], opp]);
        if side == 0.0 {
            return false;
        }
        match faces_at_b.get_mut(&f) {
            None => {
                faces_at_b.insert(f, (1, side, lab));
                true
            }
            Some(entry) if entry.0 == 1 => {
                if entry.1 * side >= 0.0 {
                    return false;
                }
                entry.0 = 2;
                if entry.2 != lab {
                    post_patch_at_b
                        .get(&f)
                        .is_some_and(|&p| matches!(w.patches[p as usize].sig, Sig::Iface(_, _) | Sig::Sheet(_)))
                } else {
                    post_patch_at_b.get(&f).is_none_or(|&p| matches!(w.patches[p as usize].sig, Sig::Sheet(_)))
                }
            }
            Some(_) => false,
        }
    };
    for &t in &snap.vert_tets[b as usize] {
        let tt = w.tets[t as usize];
        if !tt.contains(&a) {
            let lab = w.tlabel[t as usize];
            for (f, opp) in tet_faces(&tt) {
                if f.contains(&b) && !reg_face(f, opp, lab) {
                    return None;
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
        if !(tet_volume(&w.pts, &nt) > 0.0) {
            return None;
        }
        let lab = w.tlabel[t as usize];
        for (f, opp) in tet_faces(&nt) {
            if f.contains(&b) && !reg_face(f, opp, lab) {
                return None;
            }
        }
        if w.vclass[a as usize] != VClass3::Interior {
            // a surface face that moves with a must not turn by more than about 25 degrees (flat patches do not turn at all)
            for (f, _) in tet_faces(&tt) {
                if f.contains(&a) && !f.contains(&b) && w.patch_of.contains_key(&f) {
                    let g: Vec<u32> = f.iter().map(|&x| if x == a { b } else { x }).collect();
                    let (n0, n1) = (outward_normal(&w.pts, &tt, &f), outward_normal(&w.pts, &nt, &fkey(g[0], g[1], g[2])));
                    if n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < 0.9 {
                        return None;
                    }
                }
            }
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
    for (f, &(cnt, _, _)) in &faces_at_b {
        let ok = match post_patch_at_b.get(f) {
            Some(&p) => match w.patches[p as usize].sig {
                Sig::Bnd(_) => cnt == 1,
                Sig::Iface(_, _) | Sig::Sheet(_) => cnt == 2,
            },
            None => cnt == 2,
        };
        if !ok {
            return None;
        }
    }
    Some(worst)
}

/// Collapse every edge shorter than `opt.lmin` that can be collapsed (one end moves onto the other, along its patch or curve when it is on one). Returns the number.
pub fn collapse_pass(w: &mut Work3, opt: &AdaptOptions3) -> usize {
    let has_sheets = has_sheet_patches(w);
    let snap = snapshot(w);
    let mut cand: Vec<(f64, Edge)> = vec![];
    for e in all_edges(&snap) {
        let l = elen(w, e[0], e[1]);
        if l < opt.lmin {
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
            if let Some(q) = collapse_eval(w, &snap, a, b, opt, has_sheets) {
                if best.is_none_or(|x| q > x.0) {
                    best = Some((q, a, b));
                }
            }
        }
        let Some((_, a, b)) = best else { continue };
        // surface bookkeeping: faces on the edge vanish, the other surface faces and crease edges at a move to b
        if w.vclass[a as usize] != VClass3::Interior {
            let mut moved: Vec<[u32; 3]> = vec![];
            for &t in &snap.vert_tets[a as usize] {
                for (f, _) in tet_faces(&w.tets[t as usize]) {
                    if f.contains(&a) && w.patch_of.contains_key(&f) && !moved.contains(&f) {
                        moved.push(f);
                    }
                }
            }
            for f in moved {
                let p = w.patch_of.remove(&f).unwrap();
                if !f.contains(&b) {
                    let g: Vec<u32> = f.iter().map(|&x| if x == a { b } else { x }).collect();
                    w.patch_of.insert(fkey(g[0], g[1], g[2]), p);
                }
            }
            if let VClass3::Crease(_) = w.vclass[a as usize] {
                let ces: Vec<(Edge, u32)> = w.curve_of.iter().filter(|(k, _)| k.contains(&a)).map(|(k, &c)| (*k, c)).collect();
                for (k, c) in ces {
                    w.curve_of.remove(&k);
                    if !k.contains(&b) {
                        let o = if k[0] == a { k[1] } else { k[0] };
                        w.curve_of.insert(ekey(b, o), c);
                    }
                }
            }
        }
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
    // edges and faces that tetrahedra made earlier in this pass have: the snapshot does not know them, and a flip that made one a second
    // time would leave a face with four tetrahedra
    let mut new_edges: std::collections::HashSet<Edge> = Default::default();
    let mut new_faces: std::collections::HashSet<[u32; 3]> = Default::default();
    let note = |t: &[u32; 4], edges: &mut std::collections::HashSet<Edge>, faces: &mut std::collections::HashSet<[u32; 3]>| {
        for i in 0..4 {
            for j in i + 1..4 {
                edges.insert(ekey(t[i], t[j]));
            }
        }
        for (f, _) in tet_faces(t) {
            faces.insert(f);
        }
    };    // 3-2: an interior edge with three tetrahedra around it
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
        if new_faces.contains(&fkey(ring[0], ring[1], ring[2])) {
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
            note(&t, &mut new_edges, &mut new_faces);
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
        if snap.edge_tets.contains_key(&ekey(p, q)) || new_edges.contains(&ekey(p, q)) {
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
            note(&t, &mut new_edges, &mut new_faces);
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
    let has_sheets = has_sheet_patches(w);
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
                if has_sheets {
                    nb.sort_unstable();
                }
            }
            VClass3::Surface(p) => {
                for &t in star {
                    for (f, _) in tet_faces(&w.tets[t as usize]) {
                        let on_patch = if has_sheets {
                            w.patch_of.get(&f) == Some(&p)
                        } else {
                            w.patch_of.contains_key(&f)
                        };
                        if f.contains(&v) && on_patch {
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
            if !(crate::adapt::work::simplex_signed_measure::<3>(&p) > 0.0) {
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

/// The vertices around the interior edge `e` in cyclic order (each tetrahedron on the edge contributes the pair of its other vertices);
/// `None` when the edge is on the boundary of the mesh.
fn ring_of(w: &Work3, snap: &Snap3, e: &Edge) -> Option<Vec<u32>> {
    let ts = &snap.edge_tets[e];
    let mut pairs: Vec<[u32; 2]> = vec![];
    for &t in ts {
        let o: Vec<u32> = w.tets[t as usize].iter().copied().filter(|&x| x != e[0] && x != e[1]).collect();
        pairs.push([o[0], o[1]]);
    }
    // walk the cycle
    let mut ring = vec![pairs[0][0], pairs[0][1]];
    let mut used = vec![false; pairs.len()];
    used[0] = true;
    for _ in 1..pairs.len() {
        let last = *ring.last().unwrap();
        let (i, next) = pairs.iter().enumerate().find_map(|(i, p)| {
            if used[i] {
                None
            } else if p[0] == last {
                Some((i, p[1]))
            } else if p[1] == last {
                Some((i, p[0]))
            } else {
                None
            }
        })?;
        used[i] = true;
        ring.push(next);
    }
    // closed: the walk ends where it began
    if ring.pop()? != ring[0] || ring.len() != pairs.len() {
        return None;
    }
    Some(ring)
}

/// Edge removal: an interior edge with `n` (4 to 7) tetrahedra around it is replaced by the `2(n - 2)` tetrahedra of the best
/// triangulation of its ring, if that raises the worst quality. (The 3-2 flip of `flip_pass` is the case `n = 3`.) Returns the number.
pub fn edge_removal_pass(w: &mut Work3) -> usize {
    let snap = snapshot(w);
    let mut lock = vec![false; w.tets.len()];
    let mut n_done = 0;
    let mut new_edges: std::collections::HashSet<Edge> = Default::default();
    for e in all_edges(&snap) {
        let ts = &snap.edge_tets[&e];
        let n = ts.len();
        if !(4..=7).contains(&n) || locked(&lock, ts) || w.curve_of.contains_key(&e) {
            continue;
        }
        let label = w.tlabel[ts[0] as usize];
        if ts.iter().any(|&t| w.tlabel[t as usize] != label) {
            continue;
        }
        if ts.iter().any(|&t| tet_faces(&w.tets[t as usize]).iter().any(|(f, _)| f.contains(&e[0]) && f.contains(&e[1]) && w.patch_of.contains_key(f))) {
            continue;
        }
        let Some(ring) = ring_of(w, &snap, &e) else { continue };
        let old_min = ts.iter().map(|&t| tet_q(w, &w.tets[t as usize])).fold(f64::INFINITY, f64::min);
        // quality of the two tetrahedra over a ring triangle (i, j, k): the ends of the edge must lie on opposite sides
        let tri_q = |i: usize, j: usize, k: usize| -> f64 {
            let (ci, cj, ck) = (ring[i], ring[j], ring[k]);
            let va = tet_volume(&w.pts, &[ci, cj, ck, e[0]]);
            let vb = tet_volume(&w.pts, &[ci, cj, ck, e[1]]);
            if !(va * vb < 0.0) {
                return f64::NEG_INFINITY;
            }
            let ta = orient(w, [ci, cj, ck, e[0]]);
            let tb = orient(w, [ci, cj, ck, e[1]]);
            tet_q(w, &ta).min(tet_q(w, &tb))
        };
        // best triangulation of the polygon by dynamic programming over its intervals: best[i][j] = (worst quality, split vertex)
        let mut best = vec![vec![(f64::INFINITY, 0usize); n]; n];
        for len in 2..n {
            for i in 0..n - len {
                let j = i + len;
                let mut b = (f64::NEG_INFINITY, 0usize);
                for k in i + 1..j {
                    let q = tri_q(i, k, j).min(best[i][k].0).min(best[k][j].0);
                    if q > b.0 {
                        b = (q, k);
                    }
                }
                best[i][j] = b;
            }
        }
        let (q_new, _) = best[0][n - 1];
        if !(q_new > old_min) {
            continue;
        }
        let mut tris: Vec<(usize, usize, usize)> = vec![];
        let mut stack = vec![(0usize, n - 1)];
        while let Some((i, j)) = stack.pop() {
            if j - i < 2 {
                continue;
            }
            let k = best[i][j].1;
            tris.push((i, k, j));
            stack.push((i, k));
            stack.push((k, j));
        }
        // the triangles of the ring must all face the same way as the ring does, or the new tetrahedra overlap (a triangle of a non-convex
        // ring can be turned inside out and still have the ends of the edge on its two sides)
        let fan: f64 = (1..n - 1).map(|m| tet_volume(&w.pts, &[ring[0], ring[m], ring[m + 1], e[0]])).sum();
        if tris.iter().any(|&(i, k, j)| tet_volume(&w.pts, &[ring[i], ring[k], ring[j], e[0]]) * fan <= 0.0) {
            continue;
        }
        // a diagonal of the ring that is already an edge of the mesh elsewhere would be doubled: the new tetrahedra around it would
        // share a face with the ones that are already there (a face with four tetrahedra)
        let adjacent = |a: usize, b: usize| b == a + 1 || (a == 0 && b == n - 1);
        if tris.iter().any(|&(i, k, j)| [(i, k), (k, j), (i, j)].iter().any(|&(a, b)| !adjacent(a, b) && (snap.edge_tets.contains_key(&ekey(ring[a], ring[b])) || new_edges.contains(&ekey(ring[a], ring[b]))))) {
            continue;
        }
        for &t in ts.iter() {
            w.talive[t as usize] = false;
            lock[t as usize] = true;
        }
        for (i, k, j) in tris {
            for apex in [e[0], e[1]] {
                let t = orient(w, [ring[i], ring[k], ring[j], apex]);
                for (a, b) in [(i, k), (k, j), (i, j)] {
                    new_edges.insert(ekey(ring[a], ring[b]));
                }
                w.tets.push(t);
                w.tlabel.push(label);
                w.talive.push(true);
                lock.push(true);
            }
        }
        n_done += 1;
    }
    n_done
}

/// Flip of a boundary edge between two coplanar faces of a boundary patch with exactly two tetrahedra on it: the diagonal of the
/// quadrilateral changes, the surface (flat there) does not. Returns the number of flips.
pub fn surface_flip_pass(w: &mut Work3) -> usize {
    let snap = snapshot(w);
    let mut lock = vec![false; w.tets.len()];
    let mut new_faces: std::collections::HashSet<[u32; 3]> = Default::default();
    let mut n_done = 0;
    for e in all_edges(&snap) {
        let ts = &snap.edge_tets[&e];
        if ts.len() != 2 || locked(&lock, ts) || w.curve_of.contains_key(&e) {
            continue;
        }
        let sf = surface_faces_on(w, &snap, &e);
        if sf.len() != 2 || sf[0].1 != sf[1].1 {
            continue;
        }
        let p = sf[0].1;
        if !matches!(w.patches[p as usize].sig, Sig::Bnd(_)) || w.patches[p as usize].frozen {
            continue;
        }
        if w.tlabel[ts[0] as usize] != w.tlabel[ts[1] as usize] {
            continue;
        }
        let (c, d) = (*sf[0].0.iter().find(|&&x| x != e[0] && x != e[1]).unwrap(), *sf[1].0.iter().find(|&&x| x != e[0] && x != e[1]).unwrap());
        if c == d {
            continue;
        }
        // the two faces must be coplanar: the four points of the quadrilateral a, c, b, d lie in one plane
        let (pa, pb, pc, pd) = (w.pts[e[0] as usize], w.pts[e[1] as usize], w.pts[c as usize], w.pts[d as usize]);
        let (u, v, x) = (sub3(pb, pa), sub3(pc, pa), sub3(pd, pa));
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let scale = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt() * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        if (n[0] * x[0] + n[1] * x[1] + n[2] * x[2]).abs() > 1e-9 * scale {
            continue;
        }
        // the third vertex x of the two tetrahedra (the one they share besides the edge)
        let others = |t: u32| -> Vec<u32> { w.tets[t as usize].iter().copied().filter(|&y| y != e[0] && y != e[1]).collect() };
        let (o0, o1) = (others(ts[0]), others(ts[1]));
        let Some(&xv) = o0.iter().find(|y| o1.contains(y)) else { continue };
        let cand_faces = [fkey(c, d, xv), fkey(c, d, e[0]), fkey(c, d, e[1])];
        if cand_faces.iter().any(|f| {
            w.patch_of.contains_key(f)
                || new_faces.contains(f)
                || snap.vert_tets[c as usize].iter().any(|&t| {
                    t != ts[0]
                        && t != ts[1]
                        && w.talive[t as usize]
                        && f.iter().all(|r| w.tets[t as usize].contains(r))
                })
        }) {
            continue;
        }
        // old tets: (a,b,c,x) and (a,b,x,d) in some order; new: (c,d,x,a) and (c,d,x,b)
        let new = [orient(w, [c, d, xv, e[0]]), orient(w, [c, d, xv, e[1]])];
        let (va, vb) = (tet_volume(&w.pts, &[c, d, xv, e[0]]), tet_volume(&w.pts, &[c, d, xv, e[1]]));
        if !(va * vb < 0.0) {
            continue;
        }
        let old_min = ts.iter().map(|&t| tet_q(w, &w.tets[t as usize])).fold(f64::INFINITY, f64::min);
        let new_min = new.iter().map(|t| tet_q(w, t)).fold(f64::INFINITY, f64::min);
        if new_min <= old_min {
            continue;
        }
        for &t in ts {
            w.talive[t as usize] = false;
            lock[t as usize] = true;
        }
        let label = w.tlabel[ts[0] as usize];
        for t in new {
            w.tets.push(t);
            w.tlabel.push(label);
            w.talive.push(true);
            lock.push(true);
        }
        new_faces.insert(fkey(c, d, xv));
        new_faces.insert(fkey(c, d, e[0]));
        new_faces.insert(fkey(c, d, e[1]));
        for (f, _) in &sf {
            w.patch_of.remove(f);
        }
        w.patch_of.insert(fkey(c, d, e[0]), p);
        w.patch_of.insert(fkey(c, d, e[1]), p);
        n_done += 1;
    }
    n_done
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapt::topo3::PatchInfo;

    #[test]
    fn collapse_eval_rejects_same_side_interior_face_overlap_without_sheets() {
        // Without sheets, collapsing 0 -> 1 leaves tets [1,2,3,4] and [1,2,3,5] on the SAME side
        // of shared face [1,2,3] while satisfying the vertex link condition and det > 0.
        let pts = vec![
            [0.0, 0.0, -0.5],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.2, 0.2, 0.8],
            [0.2, 0.2, 0.4],
        ];
        let mut tets = vec![[0, 1, 2, 3], [1, 2, 3, 4], [0, 2, 3, 5]];
        for t in &mut tets {
            if tet_volume(&pts, t) < 0.0 {
                t.swap(2, 3);
            }
            assert!(tet_volume(&pts, t) > 0.0);
        }
        let w = Work3 {
            pts,
            met: vec![super::super::metric::identity::<3>(); 6],
            vclass: vec![VClass3::Interior; 6],
            tets,
            tlabel: vec![1, 1, 1],
            talive: vec![true, true, true],
            patches: vec![],
            curves: vec![],
            patch_of: FxMap::default(),
            curve_of: FxMap::default(),
        };
        let snap = snapshot(&w);
        let opt = AdaptOptions3 { collapse_qmin: 0.0, ..AdaptOptions3::default() };
        assert!(collapse_eval(&w, &snap, 0, 1, &opt, false).is_none());
    }

    #[test]
    fn surface_flip_pass_blocks_colliding_face_without_sheets() {
        // Long diagonal (0,1) would strictly improve quality by flipping to short diagonal (2,3),
        // but candidate face [2,3,4] already exists on tet [2,3,4,5]. Even with zero sheet patches,
        // the existing-face collision guard must reject the flip.
        let pts = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, -0.3, 0.0],
            [1.0, 0.3, 0.0],
            [1.0, 0.0, 0.8],
            [2.0, 0.0, 1.2],
        ];
        let mut tets = vec![[0, 1, 2, 4], [0, 3, 1, 4], [2, 3, 4, 5]];
        for t in &mut tets {
            if tet_volume(&pts, t) < 0.0 {
                t.swap(2, 3);
            }
            assert!(tet_volume(&pts, t) > 0.0);
        }
        let mut patch_of = FxMap::default();
        patch_of.insert(fkey(0, 1, 2), 0);
        patch_of.insert(fkey(0, 1, 3), 0);
        let mut w = Work3 {
            pts,
            met: vec![super::super::metric::identity::<3>(); 6],
            vclass: vec![VClass3::Surface(0); 6],
            tets,
            tlabel: vec![1, 1, 1],
            talive: vec![true, true, true],
            patches: vec![PatchInfo { sig: Sig::Bnd(1), frozen: false, normal: None }],
            curves: vec![],
            patch_of,
            curve_of: FxMap::default(),
        };
        assert_eq!(surface_flip_pass(&mut w), 0);
    }
}
