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

const LINE_THETA_STEP: f64 = 1e-4;

struct Curves {
    list: Vec<Curve>,
    tol: f64,
    center: V2,
    line_ids: Vec<u32>,
    line_grid: HashMap<(i64, i64), Vec<u32>>,
    circle_grid: HashMap<(i64, i64, i64), Vec<u32>>,
    ellipse_grid: HashMap<(i64, i64, i64), Vec<u32>>,
}

impl Curves {
    fn new(tol: f64, center: V2) -> Self {
        Self {
            list: Vec::new(),
            tol,
            center,
            line_ids: Vec::new(),
            line_grid: HashMap::new(),
            circle_grid: HashMap::new(),
            ellipse_grid: HashMap::new(),
        }
    }

    fn add(&mut self, c: Curve) -> u32 {
        let tol = self.tol;
        let step = (4.0 * tol).max(1e-30);
        match &c {
            Curve::Circle { c: c2, r: r2 } => {
                let (kx, ky, kr) =
                    ((c2[0] / step).floor() as i64, (c2[1] / step).floor() as i64, (*r2 / step).floor() as i64);
                let mut best: Option<u32> = None;
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dr in -1..=1 {
                            if let Some(bucket) = self.circle_grid.get(&(kx + dx, ky + dy, kr + dr)) {
                                for &id in bucket {
                                    if best.is_some_and(|b| id >= b) {
                                        continue;
                                    }
                                    if let Curve::Circle { c: c1, r: r1 } = &self.list[id as usize] {
                                        if dist2(*c1, *c2) <= tol && (r1 - r2).abs() <= tol {
                                            best = Some(id);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if let Some(id) = best {
                    return id;
                }
                self.list.push(c);
                let id = (self.list.len() - 1) as u32;
                self.circle_grid.entry((kx, ky, kr)).or_default().push(id);
                id
            }
            Curve::Ellipse { c: c2, a: a2, b: b2, phi: p2 } => {
                let (kx, ky, ka) =
                    ((c2[0] / step).floor() as i64, (c2[1] / step).floor() as i64, (*a2 / step).floor() as i64);
                let mut best: Option<u32> = None;
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for da in -1..=1 {
                            if let Some(bucket) = self.ellipse_grid.get(&(kx + dx, ky + dy, ka + da)) {
                                for &id in bucket {
                                    if best.is_some_and(|b| id >= b) {
                                        continue;
                                    }
                                    if let Curve::Ellipse { c: c1, a: a1, b: b1, phi: p1 } = &self.list[id as usize] {
                                        if dist2(*c1, *c2) <= tol
                                            && (a1 - a2).abs() <= tol
                                            && (b1 - b2).abs() <= tol
                                            && (a1 - b1) * (p1 - p2).sin().abs() <= tol
                                        {
                                            best = Some(id);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if let Some(id) = best {
                    return id;
                }
                self.list.push(c);
                let id = (self.list.len() - 1) as u32;
                self.ellipse_grid.entry((kx, ky, ka)).or_default().push(id);
                id
            }
            Curve::Line { .. } => {
                self.list.push(c);
                (self.list.len() - 1) as u32
            }
        }
    }

    /// Interned line through `p`, `q`; the flag says that `p -> q` runs against the canonical direction.
    fn line(&mut self, p: V2, q: V2) -> (u32, bool) {
        let (c_cand, rev) = canonical_line(p, sub2(q, p));
        let Curve::Line { o: o_cand, d: d_cand } = c_cand else { unreachable!() };
        let len = dist2(p, q).max(self.tol);
        let r_p = dist2(p, self.center).max(dist2(q, self.center));
        let delta_theta = (2.01 * self.tol / len) + 1e-11;
        let delta_s = delta_theta * r_p + 1.01 * self.tol;
        let s_step = (16.0 * self.tol).max(1e-30);
        let theta = d_cand[1].atan2(d_cand[0]);
        let s_cand = cross2(d_cand, sub2(o_cand, self.center));
        let kt0 = ((theta - delta_theta) / LINE_THETA_STEP).floor() as i64;
        let kt1 = ((theta + delta_theta) / LINE_THETA_STEP).floor() as i64;
        let ks0 = ((s_cand - delta_s) / s_step).floor() as i64;
        let ks1 = ((s_cand + delta_s) / s_step).floor() as i64;

        let mut best: Option<u32> = None;
        if (kt1 - kt0 + 1).saturating_mul(ks1 - ks0 + 1) <= 64 {
            for kt in kt0..=kt1 {
                for ks in ks0..=ks1 {
                    if let Some(bucket) = self.line_grid.get(&(kt, ks)) {
                        for &id in bucket {
                            if best.is_some_and(|b| id >= b) {
                                continue;
                            }
                            if let Curve::Line { o, d } = &self.list[id as usize] {
                                if cross2(*d, sub2(p, *o)).abs() <= self.tol
                                    && cross2(*d, sub2(q, *o)).abs() <= self.tol
                                {
                                    best = Some(id);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            for &id in &self.line_ids {
                if let Curve::Line { o, d } = &self.list[id as usize] {
                    if cross2(*d, sub2(p, *o)).abs() <= self.tol && cross2(*d, sub2(q, *o)).abs() <= self.tol {
                        best = Some(id);
                        break;
                    }
                }
            }
        }
        if let Some(id) = best {
            let Curve::Line { d, .. } = &self.list[id as usize] else { unreachable!() };
            return (id, dot2(sub2(q, p), *d) < 0.0);
        }
        self.list.push(c_cand);
        let id = (self.list.len() - 1) as u32;
        self.line_ids.push(id);
        let kt = (theta / LINE_THETA_STEP).floor() as i64;
        let ks = (s_cand / s_step).floor() as i64;
        self.line_grid.entry((kt, ks)).or_default().push(id);
        if d_cand[0].abs() <= 0.1 {
            let theta_opp = (-d_cand[1]).atan2(-d_cand[0]);
            let s_opp = -s_cand;
            let kt_opp = (theta_opp / LINE_THETA_STEP).floor() as i64;
            let ks_opp = (s_opp / s_step).floor() as i64;
            if (kt_opp, ks_opp) != (kt, ks) {
                self.line_grid.entry((kt_opp, ks_opp)).or_default().push(id);
            }
        }
        (id, rev)
    }
}

struct Dsu {
    parent: Vec<u32>,
    size: Vec<u32>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Self { parent: (0..n as u32).collect(), size: vec![1; n] }
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let p = self.parent[x as usize];
            self.parent[x as usize] = self.parent[p as usize];
            x = p;
        }
        x
    }
    fn union(&mut self, a: u32, b: u32) {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.size[ra as usize] < self.size[rb as usize] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb as usize] = ra;
        self.size[ra as usize] += self.size[rb as usize];
    }
}

struct Verts {
    p: Vec<V2>,
    who: Vec<String>,
    tol: f64,
    band: f64,
}

impl Verts {
    fn new(tol: f64, band: f64) -> Self {
        Self { p: vec![], who: vec![], tol, band }
    }
    fn push_raw(&mut self, p: V2, who: &dyn Fn() -> String) -> u32 {
        self.p.push(p);
        self.who.push(who());
        (self.p.len() - 1) as u32
    }
    fn cell_with(p: V2, width: f64) -> (i64, i64) {
        ((p[0] / width).floor() as i64, (p[1] / width).floor() as i64)
    }
    /// Clusters raw vertices by union-find over `d <= tol`, rejects clusters whose diameter exceeds
    /// `band` (`64 * tol`) and near-miss vertex/centroid pairs, snaps every cluster to its deterministic
    /// canonical centroid, and returns the canonicalized `Verts` plus the `raw_id -> cluster_id` map.
    fn cluster(self) -> GeomResult<(Self, Vec<u32>)> {
        let n = self.p.len();
        if n == 0 {
            return Ok((Self { p: vec![], who: vec![], tol: self.tol, band: self.band }, vec![]));
        }
        let cell_w = (self.band * (1.0 + 1e-9)).max(1e-300);
        let mut grid: HashMap<(i64, i64), Vec<u32>> = HashMap::with_capacity(n);
        for (i, &pt) in self.p.iter().enumerate() {
            grid.entry(Self::cell_with(pt, cell_w)).or_default().push(i as u32);
        }
        let mut dsu = Dsu::new(n);
        for (i, &pi) in self.p.iter().enumerate() {
            let (cx, cy) = Self::cell_with(pi, cell_w);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = grid.get(&(cx + dx, cy + dy)) {
                        for &j in bucket {
                            if (i as u32) < j && dist2(pi, self.p[j as usize]) <= self.tol {
                                dsu.union(i as u32, j);
                            }
                        }
                    }
                }
            }
        }
        // Reject raw points in distinct clusters that fall within the undecidable band `(tol, band]`.
        let mut comp_map: HashMap<u32, Vec<usize>> = HashMap::new();
        for i in 0..n {
            comp_map.entry(dsu.find(i as u32)).or_default().push(i);
        }
        struct Cluster {
            centroid: V2,
            who: String,
            members: Vec<usize>,
        }
        let mut clusters: Vec<Cluster> = Vec::with_capacity(comp_map.len());
        for (_root, members) in comp_map {
            // Check exact Euclidean diameter <= band (`64 * tol`).
            let mut bb = Aabb::EMPTY;
            for &m in &members {
                bb.grow(to3(self.p[m]));
            }
            let bb_diag = norm2(sub2([bb.hi[0], bb.hi[1]], [bb.lo[0], bb.lo[1]]));
            if bb_diag > self.band || members.len() > 1 {
                let mut worst: Option<(f64, usize, usize)> = None;
                for (a_idx, &u) in members.iter().enumerate() {
                    for &v in &members[a_idx + 1..] {
                        let d = dist2(self.p[u], self.p[v]);
                        if d > self.band && worst.is_none_or(|w| d > w.0) {
                            worst = Some((d, u, v));
                        }
                    }
                }
                if let Some((d, u, v)) = worst {
                    return Err(GeomError(format!(
                        "{} and {} span a near-coincident vertex cluster of diameter {:.3e} m (> {:.1e} m = 64 tolerances): make vertices coincide cleanly or separate them",
                        self.who[u], self.who[v], d, self.band
                    )));
                }
            }
            // Deterministic order-independent centroid over sorted unique cluster coordinates.
            let mut uniq: Vec<V2> = members
                .iter()
                .map(|&m| {
                    let p = self.p[m];
                    [if p[0] == 0.0 { 0.0 } else { p[0] }, if p[1] == 0.0 { 0.0 } else { p[1] }]
                })
                .collect();
            uniq.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
            uniq.dedup_by(|a, b| a[0].to_bits() == b[0].to_bits() && a[1].to_bits() == b[1].to_bits());
            let (mut sx, mut sy) = (0.0, 0.0);
            for p in &uniq {
                sx += p[0];
                sy += p[1];
            }
            let inv = 1.0 / (uniq.len() as f64);
            let centroid = [sx * inv, sy * inv];
            let rep_who = members.iter().map(|&m| &self.who[m]).min().unwrap().clone();
            clusters.push(Cluster { centroid, who: rep_who, members });
        }
        clusters.sort_by(|a, b| {
            a.centroid[0]
                .total_cmp(&b.centroid[0])
                .then(a.centroid[1].total_cmp(&b.centroid[1]))
                .then(a.who.cmp(&b.who))
        });

        // Assert minimum inter-cluster centroid separation > tol first, then reject pairs in `(tol, band]`.
        let mut cgrid: HashMap<(i64, i64), Vec<usize>> = HashMap::with_capacity(clusters.len());
        for (ci, c) in clusters.iter().enumerate() {
            cgrid.entry(Self::cell_with(c.centroid, cell_w)).or_default().push(ci);
        }
        for (ci, ca) in clusters.iter().enumerate() {
            let (cx, cy) = Self::cell_with(ca.centroid, cell_w);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = cgrid.get(&(cx + dx, cy + dy)) {
                        for &cj in bucket {
                            if ci < cj {
                                let d = dist2(ca.centroid, clusters[cj].centroid);
                                if d <= self.tol {
                                    return Err(GeomError(format!(
                                        "{} and {} have vertex cluster centroids only {:.3e} m apart (<= tolerance {:.1e} m), preventing a near-miss split",
                                        ca.who, clusters[cj].who, d, self.tol
                                    )));
                                }
                            }
                        }
                    }
                }
            }
        }
        // Reject raw points in distinct clusters that fall within `(tol, band]`.
        for (i, &pi) in self.p.iter().enumerate() {
            let (cx, cy) = Self::cell_with(pi, cell_w);
            let ri = dsu.find(i as u32);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = grid.get(&(cx + dx, cy + dy)) {
                        for &j in bucket {
                            if (i as u32) < j && ri != dsu.find(j) {
                                let d = dist2(pi, self.p[j as usize]);
                                if d <= self.band {
                                    return Err(GeomError(format!(
                                        "{} and {} have corners {:.3e} m apart, which is close to but not within the tolerance {:.1e} m: make them coincide or move them apart",
                                        self.who[i], self.who[j as usize], d, self.tol
                                    )));
                                }
                            }
                        }
                    }
                }
            }
        }
        for (ci, ca) in clusters.iter().enumerate() {
            let (cx, cy) = Self::cell_with(ca.centroid, cell_w);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = cgrid.get(&(cx + dx, cy + dy)) {
                        for &cj in bucket {
                            if ci < cj {
                                let d = dist2(ca.centroid, clusters[cj].centroid);
                                if d <= self.band {
                                    return Err(GeomError(format!(
                                        "{} and {} have corners {:.3e} m apart, which is close to but not within the tolerance {:.1e} m: make them coincide or move them apart",
                                        ca.who, clusters[cj].who, d, self.tol
                                    )));
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut remap = vec![0u32; n];
        let mut out_p = Vec::with_capacity(clusters.len());
        let mut out_who = Vec::with_capacity(clusters.len());
        for (cid, cl) in clusters.into_iter().enumerate() {
            for m in cl.members {
                remap[m] = cid as u32;
            }
            out_p.push(cl.centroid);
            out_who.push(cl.who);
        }
        Ok((Self { p: out_p, who: out_who, tol: self.tol, band: self.band }, remap))
    }
}

fn snapped_param(c: &Curve, p: V2, t_ref: f64) -> f64 {
    let t = c.closest(p).0;
    if c.periodic() {
        t_ref + (t - t_ref + std::f64::consts::PI).rem_euclid(TAU) - std::f64::consts::PI
    } else {
        t
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
    let scene_center = {
        let c = scene.bounds.center();
        [c[0], c[1]]
    };
    let mut curves = Curves::new(tol, scene_center);
    let mut verts = Verts::new(tol, band);
    let pname =
        |pi: u32, role: u32| format!("{}.{}", prims[pi as usize].name, prims[pi as usize].segs[role as usize].role);

    // 1. boundary pieces with raw endpoint observations
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
                    let (a_lo, span, left) = if *sweep >= 0.0 {
                        (a0.rem_euclid(TAU), *sweep, true)
                    } else {
                        ((a0 + sweep).rem_euclid(TAU), -*sweep, false)
                    };
                    closed = span >= TAU - 1e-12;
                    let cv = &curves.list[id as usize];
                    ends = if closed { None } else { Some((cv.eval(a_lo), cv.eval(a_lo + span))) };
                    (cid, a, b, in_left) = (id, a_lo, a_lo + span, left);
                }
                DGeom::Ell { c, a: ea, b: eb, phi } => {
                    let id = curves.add(canonical_ellipse(*c, *ea, *eb, *phi, tol));
                    (cid, a, b, in_left, closed, ends) = (id, 0.0, TAU, true, true, None);
                }
            }
            let cv = &curves.list[cid as usize];
            let (va, vb) = match ends {
                Some((p, q)) => (verts.push_raw(p, &who), verts.push_raw(q, &who)),
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
    let mut missing_keys: std::collections::HashSet<(u64, u64)> = std::collections::HashSet::new();
    for &(i, j) in &pairs {
        if segs[i].cid != segs[j].cid {
            let k = key(i, j);
            if cache.pairs.contains_key(&k) {
                cache.stats.pairs_cached += 1;
            } else if missing_keys.insert(k) {
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
        let who_i = || pname(si.prim, si.role);
        let who_j = || pname(sj.prim, sj.role);
        let (ci, cj) = (&curves.list[si.cid as usize], &curves.list[sj.cid as usize]);
        let adjacent = si.prim == sj.prim
            && (j == i + 1
                || (prims[si.prim as usize].segs.len() > 2
                    && i + prims[si.prim as usize].segs.len() - 1 == j
                    && segs[i].role == 0));
        if adjacent {
            if si.cid == sj.cid && si.in_left != sj.in_left {
                return Err(GeomError(format!(
                    "the profile of '{}' overlaps itself ({} and {})",
                    prims[si.prim as usize].name,
                    who_i(),
                    who_j()
                )));
            }
            if si.cid != sj.cid
                && !si.closed
                && !sj.closed
                && matches!(ci, Curve::Line { .. })
                && matches!(cj, Curve::Line { .. })
            {
                let (pa0, pa1) = (ci.eval(si.a), ci.eval(si.b));
                let (pb0, pb1) = (cj.eval(sj.a), cj.eval(sj.b));
                let d_i0 = dist2(pa0, pb0).min(dist2(pa0, pb1));
                let d_i1 = dist2(pa1, pb0).min(dist2(pa1, pb1));
                let p_i_free = if d_i0 <= d_i1 { pa1 } else { pa0 };
                let d_j0 = dist2(pb0, pa0).min(dist2(pb0, pa1));
                let d_j1 = dist2(pb1, pa0).min(dist2(pb1, pa1));
                let p_j_free = if d_j0 <= d_j1 { pb1 } else { pb0 };
                for (p, other, co) in [(p_i_free, sj, cj), (p_j_free, si, ci)] {
                    let (to, q) = co.closest(p);
                    if dist2(p, q) <= band && co.in_range(to, other.a, other.b, band / speed(co)) {
                        return Err(GeomError(format!(
                            "the profile of '{}' folds back on itself ({} and {})",
                            prims[si.prim as usize].name,
                            who_i(),
                            who_j()
                        )));
                    }
                }
            }
            continue;
        }
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
                        splits[ia].push(verts.push_raw(c.eval(tb), &|| pname(b.prim, b.role)));
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
                let v = verts.push_raw(h, &|| format!("the crossing of {} and {}", who_i(), who_j()));
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

    // 4b. canonical union-find vertex clustering over `d <= tol`, then snap endpoints and update `t0`/`t1`
    let mut closed_seam: Vec<Option<u32>> = vec![None; segs.len()];
    for (i, s) in segs.iter().enumerate() {
        if s.closed && splits[i].is_empty() {
            let p = curves.list[s.cid as usize].eval(s.a);
            closed_seam[i] = Some(verts.push_raw(p, &|| pname(s.prim, s.role)));
        }
    }
    let (verts, remap) = verts.cluster()?;
    for (i, s) in segs.iter_mut().enumerate() {
        for v in &mut splits[i] {
            *v = remap[*v as usize];
        }
        if let Some(seam) = &mut closed_seam[i] {
            *seam = remap[*seam as usize];
        }
        if !s.closed {
            s.va = remap[s.va as usize];
            s.vb = remap[s.vb as usize];
            let c = &curves.list[s.cid as usize];
            s.a = snapped_param(c, verts.p[s.va as usize], s.a);
            s.b = snapped_param(c, verts.p[s.vb as usize], s.b);
        }
    }

    // 5. sub-pieces
    let mut subs: Vec<Sub> = vec![];
    for (i, s) in segs.iter().enumerate() {
        let c = &curves.list[s.cid as usize];
        let eps = eps_t(s);
        let mut uniq_splits = splits[i].clone();
        uniq_splits.sort_unstable();
        uniq_splits.dedup();
        let ref_mid = 0.5 * (s.a + s.b);
        let mut inner: Vec<(f64, u32)> = uniq_splits
            .iter()
            .filter_map(|&v| {
                if !s.closed && (v == s.va || v == s.vb) {
                    return None;
                }
                let t = if s.closed {
                    c.param_in(c.closest(verts.p[v as usize]).0, s.a)
                } else {
                    snapped_param(c, verts.p[v as usize], ref_mid)
                };
                if s.closed || (t > s.a + eps && t < s.b - eps) {
                    Some((t, v))
                } else {
                    None
                }
            })
            .collect();
        inner.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        inner.dedup_by_key(|x| x.1);
        let span_s = (s.b - s.a).max(1e-300);
        let mut push = |t0: f64, t1: f64, v0: u32, v1: u32| {
            let owner = Owner { prim: s.prim, role: s.role, in_left: s.in_left, frac: (t0 - s.a) / span_s };
            subs.push(Sub { cid: s.cid, t0, t1, v0, v1, owners: vec![owner] })
        };
        if s.closed {
            if inner.is_empty() {
                let v = closed_seam[i].expect("closed seam vertex clustered");
                let t0 = snapped_param(c, verts.p[v as usize], s.a);
                push(t0, t0 + TAU, v, v);
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

    // 7. classify both sides (accelerated via 2D AABB BVH over planar primitives)
    let np = prims.len();
    let prim_boxes: Vec<crate::bvh::Bbox2> = prims
        .iter()
        .map(|p| {
            if p.is_curved_solid() || p.poly.is_some() {
                crate::bvh::Bbox2::EMPTY
            } else {
                crate::bvh::Bbox2::from_aabb_xy(&p.bbox)
            }
        })
        .collect();
    let prim_bvh = crate::bvh::AabbBvh2::build(prim_boxes);
    let (mut bl, mut br) = (vec![false; np], vec![false; np]);
    let mut cands: Vec<usize> = Vec::new();
    let mut kept: Vec<KEdge> = vec![];
    for s in merged {
        let c = &curves.list[s.cid as usize];
        let mid = c.eval(0.5 * (s.t0 + s.t1));
        cands.clear();
        prim_bvh.query_point(mid, band, |p| cands.push(p as usize));
        for o in &s.owners {
            cands.push(o.prim as usize);
        }
        cands.sort_unstable();
        cands.dedup();
        for &p in &cands {
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
        for &p in &cands {
            bl[p] = false;
            br[p] = false;
        }
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
    let okeys: Vec<Vec<(u32, u32)>> = kept.iter().map(|k| owner_key(&k.owners)).collect();
    let mut alive = vec![true; kept.len()];
    let mut inc: Vec<Vec<usize>> = vec![Vec::new(); verts.p.len()];
    for (e, k) in kept.iter().enumerate() {
        inc[k.v0 as usize].push(e);
        inc[k.v1 as usize].push(e);
    }
    let mut worklist: std::collections::BTreeSet<u32> = (0..verts.p.len() as u32)
        .filter(|&v| inc[v as usize].len() == 2 && inc[v as usize][0] != inc[v as usize][1])
        .collect();
    while let Some(v) = worklist.pop_first() {
        let l = &inc[v as usize];
        if l.len() != 2 || l[0] == l[1] {
            continue;
        }
        let (l0, l1) = (l[0], l[1]);
        if !alive[l0] || !alive[l1] {
            continue;
        }
        let (x, y) = if kept[l0].v1 == v && kept[l1].v0 == v {
            (l0, l1)
        } else if kept[l1].v1 == v && kept[l0].v0 == v {
            (l1, l0)
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
        if a.cid == b.cid && a.left == b.left && a.right == b.right && consecutive && okeys[x] == okeys[y] {
            alive[y] = false;
            let b_span = b.t1 - b.t0;
            let b_v1 = b.v1;
            let ax = &mut kept[x];
            let target_t1 = ax.t0 + (ax.t1 - ax.t0) + b_span;
            ax.t1 = if c.periodic() && ax.v0 == b_v1 {
                ax.t0 + TAU
            } else {
                snapped_param(c, verts.p[b_v1 as usize], target_t1)
            };
            ax.v1 = b_v1;
            inc[v as usize].clear();
            for slot in &mut inc[b_v1 as usize] {
                if *slot == y {
                    *slot = x;
                    break;
                }
            }
            if inc[b_v1 as usize].len() == 2 && inc[b_v1 as usize][0] != inc[b_v1 as usize][1] {
                worklist.insert(b_v1);
            }
        }
    }
    let kept: Vec<KEdge> = kept.into_iter().enumerate().filter_map(|(i, k)| alive[i].then_some(k)).collect();
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
        vmap[&if h.is_multiple_of(2) { k.v0 } else { k.v1 }]
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
    let side = |h: usize| if h.is_multiple_of(2) { kept[h / 2].left } else { kept[h / 2].right };

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
    let vertices: Vec<Vertex> = (0..nv)
        .map(|v| {
            let mut names: Vec<&str> = vedges[v].iter().map(|&e| edges[e as usize].name.as_str()).collect();
            names.sort_unstable();
            let name = names.join("/");
            Vertex { p: to3(vpos[v]), name, edges: vedges[v].clone() }
        })
        .collect();
    let _ = warnings_region;
    let mut cx = Complex {
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
        bvh: crate::bvh::ComplexBvh::default(),
    };
    cx.rebuild_bvh();
    for li in 0..cx.loops.len() {
        if cx.loops[li].area >= 0.0 {
            continue;
        }
        let r = cx.loops[li].region;
        let hole_box = cx.bvh.loops.item_boxes[li];
        let mut cands: Vec<usize> = Vec::new();
        if cx.regions[r as usize].loops.len() >= 32 {
            cx.bvh.loops.query_box(&hole_box, tol, |o| {
                let o = o as usize;
                if cx.loops[o].region == r
                    && cx.loops[o].area > 0.0
                    && boxes[o].overlaps(&boxes[li])
                    && boxes[o].lo[0] <= boxes[li].lo[0] + tol
                    && boxes[o].lo[1] <= boxes[li].lo[1] + tol
                    && boxes[o].hi[0] >= boxes[li].hi[0] - tol
                    && boxes[o].hi[1] >= boxes[li].hi[1] - tol
                {
                    cands.push(o);
                }
            });
            cands.sort_unstable();
        } else {
            for &x in &cx.regions[r as usize].loops {
                let o = x as usize;
                if cx.loops[o].area > 0.0
                    && boxes[o].overlaps(&boxes[li])
                    && boxes[o].lo[0] <= boxes[li].lo[0] + tol
                    && boxes[o].lo[1] <= boxes[li].lo[1] + tol
                    && boxes[o].hi[0] >= boxes[li].hi[0] - tol
                    && boxes[o].hi[1] >= boxes[li].hi[1] - tol
                {
                    cands.push(o);
                }
            }
        }
        let pick = if cands.len() == 1 {
            Some(cands[0])
        } else {
            let (e, f) = cx.loops[li].edges[0];
            let ed = &cx.edges[e as usize];
            let m = cx.curves[ed.curve as usize].eval(0.5 * (ed.t0 + ed.t1));
            let _ = f;
            let mut inside: Vec<usize> = cands.into_iter().filter(|&o| cx.point_in_loop_idx(o, m)).collect();
            inside.sort_by(|&a, &b| cx.loops[a].area.total_cmp(&cx.loops[b].area));
            inside.first().copied()
        };
        cx.loops[li].outer = pick.map(|p| p as u32);
    }
    Ok(cx)
}

/// Parity of crossings of one directed edge against the horizontal ray `[p[0], +inf) x {p[1]}`.
/// Uses exact linear/quadratic conic root solving over y-monotone sub-arcs for circles and ellipses.
pub(crate) fn edge_crosses_ray_x(ed: &Edge, curves: &[Curve], fwd: bool, p: V2) -> bool {
    let c = &curves[ed.curve as usize];
    match *c {
        Curve::Line { .. } => {
            let (a, b) = if fwd { (c.eval(ed.t0), c.eval(ed.t1)) } else { (c.eval(ed.t1), c.eval(ed.t0)) };
            (a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]
        }
        Curve::Circle { .. } | Curve::Ellipse { .. } => ellipse_arc_crosses_ray_x(ed, c, p),
    }
}

/// Exact quadratic horizontal ray-crossing test on an elliptical/circular arc `[t0, t1]` partitioned at vertical extrema.
fn ellipse_arc_crosses_ray_x(ed: &Edge, curve: &Curve, p: V2) -> bool {
    let (cc, rx, ry, phi) = match *curve {
        Curve::Circle { c, r } => (c, r, r, 0.0),
        Curve::Ellipse { c, a, b, phi } => (c, a, b, phi),
        Curve::Line { .. } => return false,
    };
    let (t0, t1) = (ed.t0, ed.t1);
    let (sphi, cphi) = phi.sin_cos();
    let ay = rx * sphi;
    let by = ry * cphi;
    let ry_ext_sq = ay * ay + by * by;
    if ry_ext_sq <= 0.0 {
        return false;
    }
    let ty = by.atan2(ay);
    let dy = p[1] - cc[1];
    // Collect vertical extrema ty + k*PI strictly inside (t0, t1) to split into at most 3 y-monotone pieces.
    let mut cuts = [t0, 0.0, 0.0, 0.0, t1];
    let mut n = 1usize;
    let mut k = ((t0 - ty) / std::f64::consts::PI).floor() as i64 - 1;
    loop {
        let tk = ty + (k as f64) * std::f64::consts::PI;
        if tk >= t1 {
            break;
        }
        if tk > t0 {
            cuts[n] = tk;
            n += 1;
            if n == 4 {
                break;
            }
        }
        k += 1;
    }
    cuts[n] = t1;
    n += 1;
    cuts[..n].sort_by(f64::total_cmp);

    let shift_x = sphi * cphi * (rx * rx - ry * ry) * dy / ry_ext_sq;
    let disc_scale = rx * ry / ry_ext_sq;
    let mut inside = false;
    let mut prev_y = curve.eval(cuts[0])[1];
    for w in cuts[..n].windows(2) {
        let (ta, tb) = (w[0], w[1]);
        let next_y = curve.eval(tb)[1];
        if ta < tb && (prev_y > p[1]) != (next_y > p[1]) {
            let mid_sin = (0.5 * (ta + tb) - ty).sin();
            let sgn = if mid_sin <= 0.0 { 1.0 } else { -1.0 };
            let dx = sgn * disc_scale * (ry_ext_sq - dy * dy).max(0.0).sqrt();
            let x_hit = cc[0] + shift_x + dx;
            if p[0] < x_hit {
                inside = !inside;
            }
        }
        prev_y = next_y;
    }
    inside
}

/// Even-odd test of `p` against a loop flattened to chords (used only to decide which outer loop a hole belongs to).
pub(crate) fn point_in_loop(l: &Loop, edges: &[Edge], curves: &[Curve], p: V2) -> bool {
    let mut inside = false;
    for &(e, fwd) in &l.edges {
        if edge_crosses_ray_x(&edges[e as usize], curves, fwd, p) {
            inside = !inside;
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verts_cluster_order_independent_canonical_centroid() {
        let tol = 1e-9;
        let band = GRAZE_K * tol;
        let raw = [
            ([0.0, 0.0], "a.p0"),
            ([0.4e-9, 0.0], "b.p0"),
            ([0.0, 0.3e-9], "c.p0"),
            ([0.5e-9, 0.5e-9], "d.p0"),
            ([1.0, 2.0], "far.p0"),
        ];
        let mut results = vec![];
        for perm in [[0, 1, 2, 3, 4], [3, 4, 1, 0, 2], [4, 2, 0, 3, 1]] {
            let mut v = Verts::new(tol, band);
            let mut ids = [0u32; 5];
            for idx in perm {
                let (pt, nm) = raw[idx];
                ids[idx] = v.push_raw(pt, &|| nm.to_string());
            }
            let (cv, remap) = v.cluster().unwrap();
            assert_eq!(cv.p.len(), 2);
            let c0 = remap[ids[0] as usize];
            for k in 1..4 {
                assert_eq!(remap[ids[k] as usize], c0);
            }
            assert_ne!(remap[ids[4] as usize], c0);
            results.push(cv.p[c0 as usize]);
        }
        assert_eq!(results[0][0].to_bits(), results[1][0].to_bits());
        assert_eq!(results[0][1].to_bits(), results[1][1].to_bits());
        assert_eq!(results[0][0].to_bits(), results[2][0].to_bits());
        assert_eq!(results[0][1].to_bits(), results[2][1].to_bits());
    }

    #[test]
    fn verts_cluster_rejects_chain_diameter_exceeding_64_tol() {
        let tol = 1e-9;
        let band = GRAZE_K * tol;
        let mut v = Verts::new(tol, band);
        for k in 0..80 {
            let _ = v.push_raw([(k as f64) * 0.9 * tol, 0.0], &|| format!("chain.{k}"));
        }
        let err = v.cluster().err().expect("chain > 64*tol diameter must fail").0;
        assert!(err.contains("diameter"), "{err}");
    }

    #[test]
    fn verts_cluster_rejects_inter_cluster_centroid_within_tol() {
        let tol = 1e-9;
        let band = GRAZE_K * tol;
        let mut v = Verts::new(tol, band);
        // Symmetric ring of 12 points at radius 1.05 * tol around the origin: adjacent ring points are ~0.54*tol
        // apart (connected within tol, diameter = 2.1*tol <= 64*tol), while every ring point is 1.05*tol > tol
        // from the center point (0, 0). The canonical centroid of the ring collides with (0, 0)!
        let r = 1.05 * tol;
        for k in 0..12 {
            let th = TAU * (k as f64) / 12.0;
            v.push_raw([r * th.cos(), r * th.sin()], &|| format!("ring.{k}"));
        }
        v.push_raw([0.0, 0.0], &|| "center".to_string());
        let err = v.cluster().err().expect("centroid collision within tol must fail").0;
        assert!(err.contains("centroids"), "{err}");
    }

    #[test]
    fn exact_circle_and_ellipse_ray_crossing_sub_chord_precision() {
        // Test sub-sagitta (1e-11) accuracy around both Circle and rotated Ellipse.
        // At 256 chords on R=1, sagitta error was ~7.5e-5; exact solver must classify points at 1e-11 offset.
        let curves =
            vec![Curve::Circle { c: [0.0, 0.0], r: 1.0 }, Curve::Ellipse { c: [2.0, 1.0], a: 1.8, b: 0.6, phi: 0.55 }];
        let mk_edge = |cid: u32, t0: f64, t1: f64| Edge {
            name: String::new(),
            curve: cid,
            t0,
            t1,
            v0: 0,
            v1: 0,
            left: None,
            right: None,
            prov: vec![],
            owner_prims: vec![],
            length: 0.0,
            bbox: Aabb::EMPTY,
        };
        let circle_edges = [mk_edge(0, 0.0, TAU)];
        let ell_edges = [mk_edge(1, 0.0, 2.3), mk_edge(1, 2.3, TAU)];
        for &ang in &[0.123, 1.4, 2.7, 4.1, 5.5] {
            let p_c = curves[0].eval(ang);
            let r_in = [p_c[0] * (1.0 - 1e-11), p_c[1] * (1.0 - 1e-11)];
            let r_out = [p_c[0] * (1.0 + 1e-11), p_c[1] * (1.0 + 1e-11)];
            assert!(edge_crosses_ray_x(&circle_edges[0], &curves, true, r_in), "circle inside at {ang}");
            assert!(!edge_crosses_ray_x(&circle_edges[0], &curves, true, r_out), "circle outside at {ang}");

            let p_e = curves[1].eval(ang);
            let d_in = [2.0 + (p_e[0] - 2.0) * (1.0 - 1e-11), 1.0 + (p_e[1] - 1.0) * (1.0 - 1e-11)];
            let d_out = [2.0 + (p_e[0] - 2.0) * (1.0 + 1e-11), 1.0 + (p_e[1] - 1.0) * (1.0 + 1e-11)];
            let parity =
                |pt: V2| ell_edges.iter().fold(false, |acc, ed| acc ^ edge_crosses_ray_x(ed, &curves, true, pt));
            assert!(parity(d_in), "ellipse inside at {ang}");
            assert!(!parity(d_out), "ellipse outside at {ang}");
        }
    }
}
