//! The planar arrangement: all primitive boundary pieces are interned onto shared curves, split at every crossing and
//! touching point, identical pieces merge (keeping every origin), each piece is classified on both sides with the painter rule
//! (owners decide their own side from orientation; other primitives are sampled at the piece midpoint, which cannot lie on
//! their boundary after the split), pieces with the same region on both sides vanish, collinear leftovers are re-merged,
//! and the surviving pieces become the edges of the region complex with their loops.
//!
//! Decisions are tolerance based on interned geometry (`tol = 1e-9 * L`). Anything the tolerance cannot decide (tangent
//! contact, near-coincident curves, a point that nearly touches a curve, self-intersecting profiles) is refused with an
//! error that names the pieces involved, instead of guessing.

use crate::complex::*;
use crate::curve::*;
use crate::error::{GeomError, GeomResult};
use crate::math::*;
use crate::model::{NameScheme, RegionKind};
use crate::prim::*;
use crate::scene::{Cache, RegionX, Scene};
use std::collections::HashMap;

/// Width of the "undecidable" band around the tolerance, in units of `tol`.
pub const GRAZE_K: f64 = 64.0;

/// Raw curve-curve intersection of two pieces (cached by the content hash of both).
#[derive(Debug, Default, Clone)]
pub struct PairRes {
    pub hits: Vec<V2>,
    pub graze: Vec<V2>,
    pub near: bool,
}

struct Seg {
    cid: u32,
    a: f64,
    b: f64,
    in_left: bool,
    prim: u32,
    role: u32,
    closed: bool,
    bb: Aabb,
    va: u32,
    vb: u32,
    h: u64,
}

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    prim: u32,
    role: u32,
    in_left: bool,
    /// position of the piece along its originating segment, in `[0, 1)`: continuous in the parameters
    frac: f64,
}

struct Sub {
    cid: u32,
    t0: f64,
    t1: f64,
    v0: u32,
    v1: u32,
    owners: Vec<Owner>,
}

struct KEdge {
    cid: u32,
    t0: f64,
    t1: f64,
    v0: u32,
    v1: u32,
    left: Option<usize>,
    right: Option<usize>,
    owners: Vec<Owner>,
}

struct Curves {
    list: Vec<Curve>,
    tol: f64,
}

impl Curves {
    fn add(&mut self, c: Curve) -> u32 {
        let tol = self.tol;
        let found = self.list.iter().position(|o| match (o, &c) {
            (Curve::Circle { c: c1, r: r1 }, Curve::Circle { c: c2, r: r2 }) => {
                dist2(*c1, *c2) <= tol && (r1 - r2).abs() <= tol
            }
            (Curve::Ellipse { c: c1, a: a1, b: b1, phi: p1 }, Curve::Ellipse { c: c2, a: a2, b: b2, phi: p2 }) => {
                dist2(*c1, *c2) <= tol
                    && (a1 - a2).abs() <= tol
                    && (b1 - b2).abs() <= tol
                    && (a1 - b1) * (p1 - p2).sin().abs() <= tol
            }
            _ => false,
        });
        match found {
            Some(i) => i as u32,
            None => {
                self.list.push(c);
                (self.list.len() - 1) as u32
            }
        }
    }
    /// Interned line through `p`, `q`; the flag says that `p -> q` runs against the canonical direction.
    fn line(&mut self, p: V2, q: V2) -> (u32, bool) {
        for (i, c) in self.list.iter().enumerate() {
            if let Curve::Line { o, d } = c {
                if cross2(*d, sub2(p, *o)).abs() <= self.tol && cross2(*d, sub2(q, *o)).abs() <= self.tol {
                    return (i as u32, dot2(sub2(q, p), *d) < 0.0);
                }
            }
        }
        let (c, rev) = canonical_line(p, sub2(q, p));
        self.list.push(c);
        ((self.list.len() - 1) as u32, rev)
    }
}

struct Verts {
    p: Vec<V2>,
    who: Vec<String>,
    grid: HashMap<(i64, i64), Vec<u32>>,
    tol: f64,
    band: f64,
}

impl Verts {
    fn cell(&self, p: V2) -> (i64, i64) {
        ((p[0] / self.band).floor() as i64, (p[1] / self.band).floor() as i64)
    }
    fn intern(&mut self, p: V2, who: &dyn Fn() -> String) -> GeomResult<u32> {
        let (cx, cy) = self.cell(p);
        let mut best: Option<(f64, u32)> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(v) = self.grid.get(&(cx + dx, cy + dy)) {
                    for &k in v {
                        let d = dist2(self.p[k as usize], p);
                        if d <= self.band && best.is_none_or(|b| d < b.0) {
                            best = Some((d, k));
                        }
                    }
                }
            }
        }
        match best {
            Some((d, k)) if d <= self.tol => Ok(k),
            Some((d, k)) => Err(GeomError(format!(
                "{} and {} have corners {:.3e} m apart, which is close to but not within the tolerance {:.1e} m: make them coincide or move them apart",
                who(),
                self.who[k as usize],
                d,
                self.tol
            ))),
            None => {
                self.p.push(p);
                self.who.push(who());
                let id = (self.p.len() - 1) as u32;
                self.grid.entry((cx, cy)).or_default().push(id);
                Ok(id)
            }
        }
    }
}

fn speed(c: &Curve) -> f64 {
    match c {
        Curve::Line { .. } => 1.0,
        Curve::Circle { r, .. } => *r,
        Curve::Ellipse { b, .. } => *b,
    }
}

#[cfg(feature = "parallel")]
fn par_map<T: Send, R: Send>(items: Vec<T>, f: impl Fn(T) -> R + Sync + Send) -> Vec<R> {
    use rayon::prelude::*;
    items.into_par_iter().map(f).collect()
}
#[cfg(not(feature = "parallel"))]
fn par_map<T, R>(items: Vec<T>, f: impl Fn(T) -> R) -> Vec<R> {
    items.into_iter().map(f).collect()
}

type PairWork = (Curve, Curve, (f64, f64), (f64, f64));

/// Builds the region complex of a scene.
pub fn arrange(scene: &Scene, cache: &mut Cache) -> GeomResult<Complex> {
    arrange_with(scene, cache, &|bit| scene.region_of(bit), &|| scene.regions.clone())
}

/// Classifier of a side: the region index (any numbering) for a set of inside bits, `None` for the void.
pub type Classify<'a> = &'a dyn Fn(&dyn Fn(usize) -> bool) -> Option<usize>;

/// The arrangement with a caller-supplied classifier (the layer stack uses composite classes for the common refinement);
/// `regions` is called once the classes are known and returns the region table the class indices refer to.
pub fn arrange_with(
    scene: &Scene,
    cache: &mut Cache,
    classify: Classify,
    regions: &dyn Fn() -> Vec<RegionX>,
) -> GeomResult<Complex> {
    let diag = scene.bounds.diagonal();
    if diag.is_nan() || diag <= 0.0 || !diag.is_finite() {
        return Err(GeomError("the model has no extent: no region has any area".into()));
    }
    let (tol, band) = (1e-9 * diag, GRAZE_K * 1e-9 * diag);
    let prims = &scene.prims;
    let mut curves = Curves { list: vec![], tol };
    let mut verts = Verts { p: vec![], who: vec![], grid: HashMap::new(), tol, band };
    let pname =
        |pi: u32, role: u32| format!("{}.{}", prims[pi as usize].name, prims[pi as usize].segs[role as usize].role);

    // 1. interned boundary pieces
    let mut segs: Vec<Seg> = vec![];
    for (pi, prim) in prims.iter().enumerate() {
        for (ri, ds) in prim.segs.iter().enumerate() {
            let (pi, ri) = (pi as u32, ri as u32);
            let who = || pname(pi, ri);
            let (cid, a, b, in_left, closed, ends);
            match &ds.geom {
                DGeom::Line { p0, p1 } => {
                    if dist2(*p0, *p1) <= band {
                        return Err(GeomError(format!(
                            "{} is shorter than {:.1e} m: it has no usable length",
                            who(),
                            band
                        )));
                    }
                    let (id, rev) = curves.line(*p0, *p1);
                    let c = &curves.list[id as usize];
                    let (ta, tb) = (c.closest(*p0).0, c.closest(*p1).0);
                    let ends_ordered = if ta <= tb { (*p0, *p1) } else { (*p1, *p0) };
                    (cid, a, b, in_left, closed, ends) = (id, ta.min(tb), ta.max(tb), !rev, false, Some(ends_ordered));
                }
                DGeom::Arc { c, r, a0, sweep } => {
                    let id = curves.add(Curve::Circle { c: *c, r: *r });
                    let a0 = a0.rem_euclid(TAU);
                    closed = *sweep >= TAU - 1e-12;
                    let cv = &curves.list[id as usize];
                    ends = if closed { None } else { Some((cv.eval(a0), cv.eval(a0 + sweep))) };
                    (cid, a, b, in_left) = (id, a0, a0 + sweep, true);
                }
                DGeom::Ell { c, a: ea, b: eb, phi } => {
                    let id = curves.add(canonical_ellipse(*c, *ea, *eb, *phi, tol));
                    (cid, a, b, in_left, closed, ends) = (id, 0.0, TAU, true, true, None);
                }
            }
            let cv = &curves.list[cid as usize];
            let (va, vb) = match ends {
                Some((p, q)) => (verts.intern(p, &who)?, verts.intern(q, &who)?),
                None => (0, 0),
            };
            let mut h = Hasher64::new();
            cv.hash_into(&mut h);
            h.f64(a).f64(b);
            let bb = cv.bbox(a, b).pad(band);
            segs.push(Seg { cid, a, b, in_left, prim: pi, role: ri, closed, bb, va, vb, h: h.finish() });
        }
    }

    // a primitive thinner than the undecidable band cannot be classified: refuse it
    {
        let o = {
            let c = scene.bounds.center();
            [c[0], c[1]]
        };
        let mut acc: Vec<(f64, f64)> = vec![(0.0, 0.0); prims.len()];
        for sg in &segs {
            let c = &curves.list[sg.cid as usize];
            let g = c.green(sg.a, sg.b, o);
            acc[sg.prim as usize].0 += if sg.in_left { g } else { -g };
            acc[sg.prim as usize].1 += c.length(sg.a, sg.b);
        }
        for (pi, (a2, per)) in acc.iter().enumerate() {
            let thick = a2.abs() / per;
            if thick <= band {
                return Err(GeomError(format!(
                    "'{}' is only about {:.2e} m thick, which is below the smallest distinguishable feature ({:.1e} m = 64 tolerances of this model); make it thicker or remove it",
                    prims[pi].name, thick, band
                )));
            }
        }
    }

    // 2. candidate pairs by sort-and-sweep on x
    let mut order: Vec<usize> = (0..segs.len()).collect();
    order.sort_by(|&i, &j| segs[i].bb.lo[0].total_cmp(&segs[j].bb.lo[0]).then(i.cmp(&j)));
    let mut pairs: Vec<(usize, usize)> = vec![];
    for (k, &i) in order.iter().enumerate() {
        for &j in &order[k + 1..] {
            if segs[j].bb.lo[0] > segs[i].bb.hi[0] {
                break;
            }
            if segs[i].bb.overlaps(&segs[j].bb) {
                pairs.push((i.min(j), i.max(j)));
            }
        }
    }
    pairs.sort_unstable();

    // 3. raw intersections, cached by content hash; misses computed (in parallel) and stored in pair order
    let key = |i: usize, j: usize| (segs[i].h.min(segs[j].h), segs[i].h.max(segs[j].h));
    let mut missing: Vec<(usize, usize)> = vec![];
    for &(i, j) in &pairs {
        if segs[i].cid != segs[j].cid {
            if cache.pairs.contains_key(&key(i, j)) {
                cache.stats.pairs_cached += 1;
            } else if !missing.iter().any(|&(a, b)| key(a, b) == key(i, j)) {
                missing.push((i, j));
            }
        }
    }
    {
        let work: Vec<PairWork> = missing
            .iter()
            .map(|&(i, j)| {
                let (si, sj) = (&segs[i], &segs[j]);
                (curves.list[si.cid as usize].clone(), curves.list[sj.cid as usize].clone(), (si.a, si.b), (sj.a, sj.b))
            })
            .collect();
        let res = par_map(work, |(ca, cb, ra, rb)| {
            let r = intersect(&ca, &cb, ra, rb, tol, GRAZE_K);
            PairRes { hits: r.hits.iter().map(|h| h.p).collect(), graze: r.graze, near: r.near }
        });
        for (&(i, j), r) in missing.iter().zip(res) {
            cache.stats.pairs_computed += 1;
            cache.pairs.insert(key(i, j), std::sync::Arc::new(r));
        }
    }

    // 4. split points, in pair order
    let mut splits: Vec<Vec<u32>> = vec![vec![]; segs.len()];
    let eps_t = |s: &Seg| 4.0 * tol / speed(&curves.list[s.cid as usize]);
    let near_seg = |s: &Seg, p: V2| -> bool {
        // p lies on the piece (within the band, parameter within slack)
        let c = &curves.list[s.cid as usize];
        let (t, q) = c.closest(p);
        dist2(p, q) <= band && c.in_range(t, s.a, s.b, band / speed(c))
    };
    for &(i, j) in &pairs {
        let (si, sj) = (&segs[i], &segs[j]);
        let adjacent = si.prim == sj.prim
            && (j == i + 1
                || (prims[si.prim as usize].segs.len() > 2
                    && i + prims[si.prim as usize].segs.len() - 1 == j
                    && segs[i].role == 0));
        if adjacent {
            continue;
        }
        let who_i = || pname(si.prim, si.role);
        let who_j = || pname(sj.prim, sj.role);
        let (ci, cj) = (&curves.list[si.cid as usize], &curves.list[sj.cid as usize]);
        if si.cid == sj.cid {
            // the same curve: split each piece at the other's end points that fall strictly inside it
            let mut overlap = false;
            for (a, b, ia, ib) in [(si, sj, i, j), (sj, si, j, i)] {
                if b.closed {
                    continue;
                }
                let c = &curves.list[a.cid as usize];
                for (vb, tb) in [(b.va, b.a), (b.vb, b.b)] {
                    let t = c.param_in(tb, a.a);
                    let _ = ib;
                    if a.closed || (t > a.a + eps_t(a) && t < a.b - eps_t(a)) {
                        splits[ia].push(verts.intern(c.eval(tb), &|| pname(b.prim, b.role))?);
                        let _ = vb;
                        overlap = true;
                    } else if c.in_range(tb, a.a, a.b, eps_t(a)) {
                        overlap = true;
                    }
                }
            }
            if si.closed || sj.closed {
                overlap = true;
            }
            if overlap && si.prim == sj.prim {
                return Err(GeomError(format!(
                    "the profile of '{}' overlaps itself ({} and {})",
                    prims[si.prim as usize].name,
                    who_i(),
                    who_j()
                )));
            }
            continue;
        }
        let raw = cache.pairs.get(&key(i, j)).cloned().unwrap_or_default();
        // a tangent contact at the end point of a piece whose own profile continues smoothly along the other curve is not a grazing
        // of two different shapes (an arc that starts tangent to the edge it follows)
        let smooth = |a: &Seg, b: &Seg, g: V2| -> bool {
            segs.iter().any(|o| {
                o.prim == a.prim
                    && o.cid == b.cid
                    && !o.closed
                    && [curves.list[o.cid as usize].eval(o.a), curves.list[o.cid as usize].eval(o.b)]
                        .iter()
                        .any(|e| dist2(*e, g) <= band)
            })
        };
        for g in &raw.graze {
            if smooth(si, sj, *g) || smooth(sj, si, *g) {
                continue;
            }
            if near_seg(si, *g) && near_seg(sj, *g) {
                return Err(GeomError(format!(
                    "{} touches {} tangentially near ({:.6e}, {:.6e}): move one of them by more than {:.1e} m, or make the curves coincide",
                    who_i(), who_j(), g[0], g[1], band
                )));
            }
        }
        if raw.near {
            let ov = |s: &Seg, o: &Seg| {
                let c = &curves.list[s.cid as usize];
                [o.a, o.b].iter().any(|&t| c.in_range(c.closest(curves.list[o.cid as usize].eval(t)).0, s.a, s.b, 0.0))
                    || near_seg(s, curves.list[o.cid as usize].eval(o.a))
            };
            if ov(si, sj) || ov(sj, si) || si.closed || sj.closed {
                return Err(GeomError(format!(
                    "{} and {} are almost, but not exactly, on top of each other (within {:.1e} m): make them coincide or separate them",
                    who_i(), who_j(), band
                )));
            }
        }
        for &h in &raw.hits {
            let (ti, tj) = (ci.closest(h).0, cj.closest(h).0);
            if ci.in_range(ti, si.a, si.b, eps_t(si)) && cj.in_range(tj, sj.a, sj.b, eps_t(sj)) {
                if si.prim == sj.prim {
                    return Err(GeomError(format!(
                        "the profile of '{}' crosses itself ({} and {})",
                        prims[si.prim as usize].name,
                        who_i(),
                        who_j()
                    )));
                }
                let v = verts.intern(h, &|| format!("the crossing of {} and {}", who_i(), who_j()))?;
                splits[i].push(v);
                splits[j].push(v);
            }
        }
        // end points lying on the other curve (T junctions), and near misses
        for (a, b, ia, ib, cb) in [(si, sj, i, j, cj), (sj, si, j, i, ci)] {
            if a.closed {
                continue;
            }
            let ca = &curves.list[a.cid as usize];
            for (v, t) in [(a.va, a.a), (a.vb, a.b)] {
                let p = ca.eval(t);
                let (tb, q) = cb.closest(p);
                let d = dist2(p, q);
                let _ = ia;
                if d <= band && cb.in_range(tb, b.a, b.b, eps_t(b)) {
                    if d <= tol {
                        splits[ib].push(v);
                    } else {
                        return Err(GeomError(format!(
                            "the corner of {} lies {:.3e} m from {}: closer than a model feature but not on it (tolerance {:.1e} m); put it exactly on the curve or move it away",
                            pname(a.prim, a.role), d, pname(b.prim, b.role), tol
                        )));
                    }
                }
            }
        }
    }

    // 5. sub-pieces
    let mut subs: Vec<Sub> = vec![];
    for (i, s) in segs.iter().enumerate() {
        let c = &curves.list[s.cid as usize];
        let eps = eps_t(s);
        let mut inner: Vec<(f64, u32)> = splits[i]
            .iter()
            .filter_map(|&v| {
                let t = c.param_in(c.closest(verts.p[v as usize]).0, s.a);
                if s.closed || (t > s.a + eps && t < s.b - eps) {
                    Some((t, v))
                } else {
                    None
                }
            })
            .collect();
        inner.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        inner.dedup_by_key(|x| x.1);
        let mut push = |t0: f64, t1: f64, v0: u32, v1: u32| {
            let owner = Owner { prim: s.prim, role: s.role, in_left: s.in_left, frac: (t0 - s.a) / (s.b - s.a) };
            subs.push(Sub { cid: s.cid, t0, t1, v0, v1, owners: vec![owner] })
        };
        if s.closed {
            if inner.is_empty() {
                let v = verts.intern(c.eval(s.a), &|| pname(s.prim, s.role))?;
                push(s.a, s.a + TAU, v, v);
            } else {
                for k in 0..inner.len() {
                    let (t0, v0) = inner[k];
                    let (t1, v1) = if k + 1 < inner.len() { inner[k + 1] } else { (inner[0].0 + TAU, inner[0].1) };
                    push(t0, t1, v0, v1);
                }
            }
        } else {
            let mut pts = vec![(s.a, s.va)];
            pts.extend(inner);
            pts.push((s.b, s.vb));
            for w in pts.windows(2) {
                push(w[0].0, w[1].0, w[0].1, w[1].1);
            }
        }
    }

    // 6. merge identical pieces from different primitives
    let mut idx: Vec<usize> = (0..subs.len()).collect();
    idx.sort_by(|&x, &y| {
        let (a, b) = (&subs[x], &subs[y]);
        (a.cid, a.v0, a.v1).cmp(&(b.cid, b.v0, b.v1)).then((a.t1 - a.t0).total_cmp(&(b.t1 - b.t0)))
    });
    let mut merged: Vec<Sub> = vec![];
    let mut taken: Vec<Option<Sub>> = subs.into_iter().map(Some).collect();
    for k in idx {
        let s = taken[k].take().unwrap();
        if let Some(m) = merged.last_mut() {
            let pt = 4.0 * tol / speed(&curves.list[s.cid as usize]);
            if m.cid == s.cid && m.v0 == s.v0 && m.v1 == s.v1 && ((m.t1 - m.t0) - (s.t1 - s.t0)).abs() <= pt {
                m.owners.extend(s.owners);
                continue;
            }
        }
        merged.push(s);
    }

    // 7. classify both sides
    let np = prims.len();
    let (mut bl, mut br) = (vec![false; np], vec![false; np]);
    let mut kept: Vec<KEdge> = vec![];
    for s in merged {
        let c = &curves.list[s.cid as usize];
        let mid = c.eval(0.5 * (s.t0 + s.t1));
        for p in 0..np {
            let own = s.owners.iter().find(|o| o.prim as usize == p);
            if let Some(o) = own {
                bl[p] = o.in_left;
                br[p] = !o.in_left;
            } else {
                let d = prims[p].sdf_near(mid, band);
                if d.abs() <= band {
                    let o = &s.owners[0];
                    return Err(GeomError(format!(
                        "{} runs along the boundary of '{}' without sharing its curve (within {:.1e} m); make the two coincide exactly",
                        pname(o.prim, o.role), prims[p].name, band
                    )));
                }
                bl[p] = d < 0.0;
                br[p] = bl[p];
            }
        }
        let left = classify(&|i| bl[i]);
        let right = classify(&|i| br[i]);
        if left != right {
            kept.push(KEdge { cid: s.cid, t0: s.t0, t1: s.t1, v0: s.v0, v1: s.v1, left, right, owners: s.owners });
        }
    }

    // 8. re-merge collinear fragments that the dropped edges split (same curve, sides and origins)
    let owner_key = |o: &[Owner]| -> Vec<(u32, u32)> {
        let mut v: Vec<(u32, u32)> = o.iter().map(|x| (x.prim, x.role)).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    loop {
        let mut inc: HashMap<u32, Vec<usize>> = HashMap::new();
        for (e, k) in kept.iter().enumerate() {
            inc.entry(k.v0).or_default().push(e);
            inc.entry(k.v1).or_default().push(e);
        }
        let mut vs: Vec<u32> = inc.keys().copied().collect();
        vs.sort_unstable();
        let mut did = None;
        for v in vs {
            let l = &inc[&v];
            if l.len() != 2 || l[0] == l[1] {
                continue;
            }
            let (x, y) = if kept[l[0]].v1 == v && kept[l[1]].v0 == v {
                (l[0], l[1])
            } else if kept[l[1]].v1 == v && kept[l[0]].v0 == v {
                (l[1], l[0])
            } else {
                continue;
            };
            let (a, b) = (&kept[x], &kept[y]);
            let pt = 4.0 * tol / speed(&curves.list[a.cid as usize]);
            let c = &curves.list[a.cid as usize];
            let consecutive = if c.periodic() {
                (wrap_from(b.t0, a.t1) - a.t1).min(TAU - (wrap_from(b.t0, a.t1) - a.t1)) <= pt
            } else {
                (b.t0 - a.t1).abs() <= pt
            };
            if a.cid == b.cid
                && a.left == b.left
                && a.right == b.right
                && consecutive
                && owner_key(&a.owners) == owner_key(&b.owners)
            {
                did = Some((x, y));
                break;
            }
        }
        let Some((x, y)) = did else { break };
        let b = kept.remove(y);
        let a = &mut kept[if y < x { x - 1 } else { x }];
        let sweep = (a.t1 - a.t0) + (b.t1 - b.t0);
        a.t1 = a.t0 + sweep;
        a.v1 = b.v1;
    }
    finish(scene, &regions(), curves, verts, kept, tol, prims)
}

/// Orders edges and vertices, builds loops and regions, names everything.
fn finish(
    scene: &Scene,
    sreg: &[RegionX],
    curves: Curves,
    verts: Verts,
    mut kept: Vec<KEdge>,
    tol: f64,
    prims: &[std::sync::Arc<Prim>],
) -> GeomResult<Complex> {
    let curves = curves.list;
    let area_origin = {
        let c = scene.bounds.center();
        [c[0], c[1]]
    };
    // stable edge order, independent of parameter values: first origin (feature id, instance, role), then the position of the
    // piece along that origin's segment (a fraction, continuous in the parameters), then the curve id
    let key_of = |k: &KEdge| -> (u32, String, String, i64, u32) {
        k.owners
            .iter()
            .map(|o| {
                let p = &prims[o.prim as usize];
                (
                    p.feature,
                    p.name.clone(),
                    p.segs[o.role as usize].role.clone(),
                    (o.frac * 1.0e12).round() as i64,
                    k.cid,
                )
            })
            .min()
            .unwrap()
    };
    let keys: Vec<_> = kept.iter().map(key_of).collect();
    let mut ord: Vec<usize> = (0..kept.len()).collect();
    ord.sort_by(|&a, &b| keys[a].cmp(&keys[b]).then(a.cmp(&b)));
    let mut slots: Vec<Option<KEdge>> = kept.drain(..).map(Some).collect();
    let kept: Vec<KEdge> = ord.iter().map(|&i| slots[i].take().unwrap()).collect();

    // vertices in order of first appearance
    let mut vmap: HashMap<u32, u32> = HashMap::new();
    let mut vpos: Vec<V2> = vec![];
    for k in &kept {
        for v in [k.v0, k.v1] {
            vmap.entry(v).or_insert_with(|| {
                vpos.push(verts.p[v as usize]);
                (vpos.len() - 1) as u32
            });
        }
    }
    let nv = vpos.len();
    let ne = kept.len();

    // half-edge rotation system
    let mut out: Vec<Vec<(f64, usize)>> = vec![vec![]; nv];
    let origin = |h: usize| -> u32 {
        let k = &kept[h / 2];
        vmap[&if h % 2 == 0 { k.v0 } else { k.v1 }]
    };
    for (e, k) in kept.iter().enumerate() {
        let c = &curves[k.cid as usize];
        let (df, db) = (c.tangent(k.t0), mul2(c.tangent(k.t1), -1.0));
        out[vmap[&k.v0] as usize].push((df[1].atan2(df[0]), 2 * e));
        out[vmap[&k.v1] as usize].push((db[1].atan2(db[0]), 2 * e + 1));
    }
    let mut pos = vec![0usize; 2 * ne];
    for (v, list) in out.iter_mut().enumerate() {
        list.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        for i in 0..list.len() {
            let (a, b) = (list[i].0, list[(i + 1) % list.len()].0);
            let gap = if list.len() > 1 { (b - a).rem_euclid(TAU) } else { TAU };
            if list.len() > 1 && gap < 1e-9 {
                let name = |h: usize| {
                    let k = &kept[h / 2];
                    let o = &k.owners[0];
                    format!("{}.{}", prims[o.prim as usize].name, prims[o.prim as usize].segs[o.role as usize].role)
                };
                return Err(GeomError(format!(
                    "{} and {} leave the point ({:.6e}, {:.6e}) in the same direction (tangent contact); the arrangement cannot order them",
                    name(list[i].1), name(list[(i + 1) % list.len()].1), vpos[v][0], vpos[v][1]
                )));
            }
            pos[list[i].1] = i;
        }
    }
    let next = |h: usize| -> usize {
        let tw = h ^ 1;
        let v = origin(tw) as usize;
        let l = &out[v];
        l[(pos[tw] + l.len() - 1) % l.len()].1
    };
    let side = |h: usize| if h % 2 == 0 { kept[h / 2].left } else { kept[h / 2].right };

    // cycles
    let mut seen = vec![false; 2 * ne];
    struct Cyc {
        region: usize,
        edges: Vec<(u32, bool)>,
        area: f64,
    }
    let mut cycles: Vec<Cyc> = vec![];
    for h0 in 0..2 * ne {
        if seen[h0] {
            continue;
        }
        let (mut h, mut list, mut a2) = (h0, vec![], 0.0);
        loop {
            seen[h] = true;
            let k = &kept[h / 2];
            let g = curves[k.cid as usize].green(k.t0, k.t1, area_origin);
            a2 += if h % 2 == 0 { g } else { -g };
            list.push(((h / 2) as u32, h % 2 == 0));
            h = next(h);
            if h == h0 {
                break;
            }
            if seen[h] {
                return Err(GeomError("internal error: the boundary walk revisited an edge".into()));
            }
        }
        let Some(r) = side(h0) else { continue };
        if list.iter().any(|&(e, f)| side(2 * e as usize + usize::from(!f)) != Some(r)) {
            if std::env::var("GEOM_DEBUG").is_ok() {
                for (e, k) in kept.iter().enumerate() {
                    eprintln!(
                        "edge {e}: cid {} t {:.6}..{:.6} v{}->v{} L{:?} R{:?}",
                        k.cid, k.t0, k.t1, k.v0, k.v1, k.left, k.right
                    );
                }
                for (v, p) in vpos.iter().enumerate() {
                    eprintln!("vertex {v}: {:?} out {:?}", p, out[v]);
                }
                eprintln!("cycle {:?}", list);
            }
            return Err(GeomError("internal error: a boundary loop changes region".into()));
        }
        cycles.push(Cyc { region: r, edges: list, area: 0.5 * a2 });
    }

    // regions that exist: all declared ones, and generated (Pml) ones only when they have area
    let nreg = sreg.len();
    let mut has_loop = vec![false; nreg];
    for c in &cycles {
        has_loop[c.region] = true;
    }
    let mut warnings_region = vec![];
    let mut fin: Vec<Option<u32>> = vec![None; nreg];
    let mut regions: Vec<Region> = vec![];
    for (i, r) in sreg.iter().enumerate() {
        let generated = r.pml.is_some();
        if matches!(r.kind, RegionKind::Hole | RegionKind::Void) {
            continue;
        }
        if generated && !has_loop[i] {
            continue;
        }
        if !has_loop[i] {
            warnings_region.push(r.name.clone());
        }
        fin[i] = Some(regions.len() as u32);
        regions.push(Region {
            name: r.name.clone(),
            kind: r.kind,
            feature: r.feature,
            pinned_id: r.pinned_id,
            scene_index: i,
            measure: 0.0,
            bbox: Aabb::EMPTY,
            loops: vec![],
            parent: None,
        });
    }
    for i in 0..regions.len() {
        let si = regions[i].scene_index;
        if let Some(pn) = &sreg[si].parent {
            regions[i].parent = regions.iter().position(|r| &r.name == pn).map(|p| p as u32);
        }
    }

    // edges
    let leaves: Vec<Vec<usize>> = sreg
        .iter()
        .map(|r| {
            let mut v = vec![];
            r.expr.leaves(&mut v);
            v
        })
        .collect();
    let mut edges: Vec<Edge> = Vec::with_capacity(ne);
    let mut vedges: Vec<Vec<u32>> = vec![vec![]; nv];
    for (e, k) in kept.iter().enumerate() {
        let c = &curves[k.cid as usize];
        let (v0, v1) = (vmap[&k.v0], vmap[&k.v1]);
        let mut prov: Vec<Prov> = k
            .owners
            .iter()
            .map(|o| {
                let p = &prims[o.prim as usize];
                Prov { feature: p.feature, instance: p.name.clone(), role: p.segs[o.role as usize].role.clone() }
            })
            .collect();
        prov.sort();
        prov.dedup();
        let (left, right) = (k.left.and_then(|r| fin[r]), k.right.and_then(|r| fin[r]));
        // auto-name: interfaces between two regions are `A|B` (sorted); boundaries use the first origin on the region's side
        let name = match (left, right) {
            (Some(a), Some(b)) => {
                let (mut x, mut y) = (regions[a as usize].name.clone(), regions[b as usize].name.clone());
                if x > y {
                    std::mem::swap(&mut x, &mut y);
                }
                format!("{x}|{y}")
            }
            _ => {
                let inside_left = left.is_some();
                // the region on the inside owns the name: prefer pieces of its own shapes, then pieces whose interior faces it
                let own_leaves = &leaves[k.left.or(k.right).unwrap_or(0)];
                let mut best: Option<(u8, u32, String, String)> = None;
                for o in &k.owners {
                    let p = &prims[o.prim as usize];
                    let rank =
                        u8::from(!own_leaves.contains(&(o.prim as usize))) * 2 + u8::from(o.in_left != inside_left);
                    let cand = (rank, p.feature, p.name.clone(), p.segs[o.role as usize].role.clone());
                    if best.as_ref().is_none_or(|b| cand < *b) {
                        best = Some(cand);
                    }
                }
                let b = best.unwrap();
                match scene.naming {
                    NameScheme::Feature => format!("{}.{}", b.2, b.3),
                    _ => format!("{}.{}", regions[left.or(right).unwrap() as usize].name, b.3),
                }
            }
        };
        let length = c.length(k.t0, k.t1);
        edges.push(Edge {
            name,
            curve: k.cid,
            t0: k.t0,
            t1: k.t1,
            v0,
            v1,
            left,
            right,
            prov,
            owner_prims: {
                let mut o: Vec<u32> = k.owners.iter().map(|x| x.prim).collect();
                o.sort_unstable();
                o.dedup();
                o
            },
            length,
            bbox: c.bbox(k.t0, k.t1),
        });
        vedges[v0 as usize].push(e as u32);
        if v1 != v0 {
            vedges[v1 as usize].push(e as u32);
        }
    }

    if scene.naming == NameScheme::Numbered {
        let mut count = vec![0usize; regions.len()];
        for ed in edges.iter_mut().filter(|e| e.left.is_none() || e.right.is_none()) {
            let r = ed.left.or(ed.right).unwrap() as usize;
            count[r] += 1;
            ed.name = format!("{}.{}", regions[r].name, count[r]);
        }
    }

    // loops, holes, region measures
    let mut loops: Vec<Loop> = vec![];
    for c in &cycles {
        let Some(r) = fin[c.region] else {
            return Err(GeomError(format!("internal error: boundary of dropped region '{}'", sreg[c.region].name)));
        };
        loops.push(Loop { region: r, edges: c.edges.clone(), area: c.area, outer: None });
    }
    let lbox = |l: &Loop| -> Aabb { l.edges.iter().fold(Aabb::EMPTY, |b, &(e, _)| b.union(&edges[e as usize].bbox)) };
    let boxes: Vec<Aabb> = loops.iter().map(lbox).collect();
    for li in 0..loops.len() {
        regions[loops[li].region as usize].loops.push(li as u32);
        regions[loops[li].region as usize].measure += loops[li].area;
        if loops[li].area > 0.0 {
            let rb = regions[loops[li].region as usize].bbox;
            regions[loops[li].region as usize].bbox = rb.union(&boxes[li]);
        }
    }
    for li in 0..loops.len() {
        if loops[li].area >= 0.0 {
            continue;
        }
        let r = loops[li].region;
        let cands: Vec<usize> = regions[r as usize]
            .loops
            .iter()
            .map(|&x| x as usize)
            .filter(|&o| {
                loops[o].area > 0.0
                    && boxes[o].overlaps(&boxes[li])
                    && boxes[o].lo[0] <= boxes[li].lo[0] + tol
                    && boxes[o].lo[1] <= boxes[li].lo[1] + tol
                    && boxes[o].hi[0] >= boxes[li].hi[0] - tol
                    && boxes[o].hi[1] >= boxes[li].hi[1] - tol
            })
            .collect();
        let pick = if cands.len() == 1 {
            Some(cands[0])
        } else {
            let (e, f) = loops[li].edges[0];
            let ed = &edges[e as usize];
            let m = curves[ed.curve as usize].eval(0.5 * (ed.t0 + ed.t1));
            let _ = f;
            let mut inside: Vec<usize> =
                cands.into_iter().filter(|&o| point_in_loop(&loops[o], &edges, &curves, m)).collect();
            inside.sort_by(|&a, &b| loops[a].area.total_cmp(&loops[b].area));
            inside.first().copied()
        };
        loops[li].outer = pick.map(|p| p as u32);
    }
    let vertices: Vec<Vertex> = (0..nv)
        .map(|v| {
            let mut names: Vec<&str> = vedges[v].iter().map(|&e| edges[e as usize].name.as_str()).collect();
            names.sort_unstable();
            let name = names.join("/");
            Vertex { p: to3(vpos[v]), name, edges: vedges[v].clone() }
        })
        .collect();
    let _ = warnings_region;
    Ok(Complex {
        dim: 2,
        tol,
        origin: area_origin,
        curves,
        regions,
        edges,
        vertices,
        loops,
        scene_to_region: fin,
        faces: vec![],
        edges3: vec![],
        planar: vec![],
        axis: 2,
        stack: None,
    })
}

/// Even-odd test of `p` against a loop flattened to chords (used only to decide which outer loop a hole belongs to).
pub(crate) fn point_in_loop(l: &Loop, edges: &[Edge], curves: &[Curve], p: V2) -> bool {
    let mut inside = false;
    for &(e, fwd) in &l.edges {
        let ed = &edges[e as usize];
        let c = &curves[ed.curve as usize];
        let n = if matches!(c, Curve::Line { .. }) { 1 } else { 256 };
        let mut prev = c.eval(if fwd { ed.t0 } else { ed.t1 });
        for k in 1..=n {
            let u = k as f64 / n as f64;
            let t = if fwd { ed.t0 + (ed.t1 - ed.t0) * u } else { ed.t1 - (ed.t1 - ed.t0) * u };
            let q = c.eval(t);
            if (prev[1] > p[1]) != (q[1] > p[1])
                && p[0] < (q[0] - prev[0]) * (p[1] - prev[1]) / (q[1] - prev[1]) + prev[0]
            {
                inside = !inside;
            }
            prev = q;
        }
    }
    inside
}
