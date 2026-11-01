//! Exact queries on an evaluated 2D model: measures, curve evaluation with first and second derivatives, curvature radii,
//! chord / facet deviation, edge discretisation by deviation and size, midpoint projection (where a quadratic element puts its
//! mid-edge node), parametric fractions of boundary nodes, and the topology signature. Everything is analytic or tabulated from
//! the exact curves; nothing here meshes.
//!
//! Parametrisation conventions (they vary continuously with the feature parameters):
//! * line: `t` is arc length along the canonical direction (angle in (-90, 90] degrees), measured from the foot of the
//!   perpendicular from the origin; a segment keeps its end points as `t0 < t1`;
//! * circle: `t` is the counter-clockwise angle from +x; an arc is `[a, a + sweep]` with `a` in `[0, 2 pi)`; a full circle has
//!   `t1 - t0 = 2 pi` and its seam at `t0`, which is the angle where the circle's first boundary piece starts (0 unless a transform
//!   rotated it);
//! * ellipse: `t` is the eccentric angle from the major axis; the seam is at the major-axis end.
//!
//! The fraction of a point on an edge is `(t - t0) / (t1 - t0)` in `[0, 1]`.

use crate::curve::Curve;
use crate::eval::Evaluated;
use crate::math::*;
use crate::sim::EntityRef;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscretizeOpts {
    pub max_dev: f64,
    pub max_len: f64,
}

/// Parametric position on an entity: one fraction along an edge (`v` is 0), two for a face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fractions {
    pub u: f64,
    pub v: f64,
}

impl Evaluated {
    // ---- measures (analytic) -------------------------------------------------------------------------------

    /// Length of an edge.
    pub fn length(&self, edge: u32) -> f64 {
        self.complex().edges[edge as usize].length
    }
    /// Area of a region (volume in 3D).
    pub fn volume(&self, region: u32) -> f64 {
        self.complex().regions[region as usize].measure
    }
    /// Measure of the interface between two regions: length of the shared edges in 2D, area in 3D.
    pub fn interface_area(&self, a: u32, b: u32) -> f64 {
        self.complex().interface_length(a, b)
    }

    // ---- curves ---------------------------------------------------------------------------------------------

    // Points and directions are in world coordinates: a 2D model lies in the x-z plane, `(x, 0, z)`.

    pub fn edge_eval(&self, e: u32, t: f64) -> V3 {
        let cx = self.complex();
        cx.wv(cx.edge_point(e, t))
    }
    pub fn edge_tangent(&self, e: u32, t: f64) -> V3 {
        let cx = self.complex();
        cx.wv(to3(cx.curve_of(e).0.tangent(t)))
    }
    /// Signed curvature along the canonical direction (positive: turns left).
    pub fn edge_curvature(&self, e: u32, t: f64) -> f64 {
        self.complex().curve_of(e).0.curvature(t)
    }
    /// Closest point of the edge (between its end points) and its parameter.
    pub fn edge_closest(&self, e: u32, p: V3) -> (V3, f64) {
        let cx = self.complex();
        let p = cx.iv(p);
        let (t, q) = cx.edge_closest(e, [p[0], p[1]]);
        (cx.wv(to3(q)), t)
    }

    /// Smallest radius of curvature on the edge (infinite for a line); a mesher derives `n = pi / acos(1 - d / r)` chords from it.
    pub fn min_radius_of_curvature(&self, e: u32) -> f64 {
        let (c, t0, t1) = self.complex().curve_of(e);
        let mut k = c.curvature(t0).abs().max(c.curvature(t1).abs());
        if c.periodic() {
            for q in 0..4 {
                let t = q as f64 * std::f64::consts::FRAC_PI_2;
                if c.in_range(t, t0, t1, 0.0) {
                    k = k.max(c.curvature(t).abs());
                }
            }
        }
        if k == 0.0 {
            f64::INFINITY
        } else {
            1.0 / k
        }
    }

    /// Largest distance between the curve on `[ta, tb]` and its chord (exact for lines and circles, sampled for ellipses).
    pub fn chord_deviation(&self, e: u32, ta: f64, tb: f64) -> f64 {
        chord_dev(self.complex().curve_of(e).0, ta, tb)
    }

    /// Largest distance from the polyline through `pts` (straight facets) to the exact edge, sampled along each segment.
    pub fn facet_deviation(&self, e: u32, pts: &[V3]) -> f64 {
        let mut worst = 0.0f64;
        for w in pts.windows(2) {
            for k in 0..=16 {
                let p = self.complex().iv(add3(mul3(w[0], 1.0 - k as f64 / 16.0), mul3(w[1], k as f64 / 16.0)));
                let q = self.complex().edge_closest(e, [p[0], p[1]]).1;
                worst = worst.max(dist2([p[0], p[1]], q));
            }
        }
        worst
    }

    /// The point of the exact edge nearest the midpoint of the chord between parameters `ta` and `tb` (the mid-edge node of a
    /// quadratic element). Returns `(t, point)` with `t` between `ta` and `tb`; seam-safe (parameters are taken relative to `ta`).
    pub fn point_at_mid(&self, e: u32, ta: f64, tb: f64) -> (f64, V3) {
        let c = self.complex().curve_of(e).0;
        let m = mul2(add2(c.eval(ta), c.eval(tb)), 0.5);
        let (t, q) = c.closest(m);
        let u = c.param_in(t, ta.min(tb));
        let cx = self.complex();
        if u >= ta.min(tb) && u <= ta.max(tb) {
            (u, cx.wv(to3(q)))
        } else {
            let tm = 0.5 * (ta + tb);
            (tm, cx.wv(to3(c.eval(tm))))
        }
    }

    /// The point of the exact edge nearest the chord midpoint of two positions on (or near) the edge; valid for chords up to a half turn.
    pub fn project_midpoint(&self, e: u32, a: V3, b: V3) -> V3 {
        self.edge_closest(e, mul3(add3(a, b), 0.5)).0
    }

    /// Points `(t, p)` along the edge, end points included, such that the chord deviation is at most `max_dev` and each chord is at
    /// most `max_len` long. Depends only on the edge, so both neighbouring regions get identical points; deterministic.
    pub fn discretize_edge(&self, e: u32, o: DiscretizeOpts) -> Vec<(f64, V3)> {
        let cx = self.complex();
        cx.discretize_edge(e, o).into_iter().map(|(t, p)| (t, cx.wv(p))).collect()
    }

    // ---- parametric positions (parameter sweeps) -----------------------------------------------------------------------------

    /// Fractions of the point of the entity closest to `p`. Edges: `u` in `[0, 1]` along `[t0, t1]`. Returns `None` for entities
    /// without a parametrisation (2D regions) or for unknown ids.
    pub fn fractions_of(&self, ent: EntityRef, p: V3) -> Option<Fractions> {
        match ent {
            EntityRef::Edge(i) if self.complex().dim == 3 && (i as usize) < self.complex().edges3.len() => {
                Some(Fractions { u: self.complex().edge3_closest(i, p).1, v: 0.0 })
            }
            EntityRef::Edge(i) if (i as usize) < self.complex().edges.len() => {
                let ed = &self.complex().edges[i as usize];
                let p = self.complex().iv(p);
                let (t, _) = self.complex().edge_closest(i, [p[0], p[1]]);
                Some(Fractions { u: ((t - ed.t0) / (ed.t1 - ed.t0)).clamp(0.0, 1.0), v: 0.0 })
            }
            EntityRef::Face(i) if (i as usize) < self.complex().faces.len() => {
                let cx = self.complex();
                let q = cx.iv(cx.facet_project(i, p));
                let f = &cx.faces[i as usize];
                Some(match f.geom {
                    crate::complex::FaceGeom::Lateral { edge, z0, z1 } => {
                        let ov = &cx.stack.as_ref()?.overlay;
                        let ed = &ov.edges[edge as usize];
                        let (t, _) = ov.edge_closest(edge, [q[0], q[1]]);
                        Fractions {
                            u: ((t - ed.t0) / (ed.t1 - ed.t0)).clamp(0.0, 1.0),
                            v: ((q[2] - z0) / (z1 - z0)).clamp(0.0, 1.0),
                        }
                    }
                    crate::complex::FaceGeom::Quadric(b) => {
                        let (u, v) = b.fractions(q);
                        Fractions { u, v }
                    }
                    crate::complex::FaceGeom::Planar { plane, cell } => {
                        let pd = &cx.planar[plane as usize];
                        let (b, d) = (pd.cx.regions[cell as usize].bbox, sub3(q, pd.o));
                        let fr = |x: f64, k: usize| {
                            if b.hi[k] > b.lo[k] {
                                ((x - b.lo[k]) / (b.hi[k] - b.lo[k])).clamp(0.0, 1.0)
                            } else {
                                0.0
                            }
                        };
                        Fractions { u: fr(dot3(d, pd.u), 0), v: fr(dot3(d, pd.v), 1) }
                    }
                    crate::complex::FaceGeom::Cap { cell, .. } => {
                        let b = cx.stack.as_ref()?.overlay.regions[cell as usize].bbox;
                        let frac = |x: f64, k: usize| {
                            if b.hi[k] > b.lo[k] {
                                ((x - b.lo[k]) / (b.hi[k] - b.lo[k])).clamp(0.0, 1.0)
                            } else {
                                0.0
                            }
                        };
                        Fractions { u: frac(q[0], 0), v: frac(q[1], 1) }
                    }
                })
            }
            EntityRef::Vertex(i) if (i as usize) < self.complex().vertices.len() => Some(Fractions { u: 0.0, v: 0.0 }),
            _ => None,
        }
    }

    /// The exact point of the entity at the given fractions.
    pub fn position_at(&self, ent: EntityRef, fr: Fractions) -> Option<V3> {
        match ent {
            EntityRef::Edge(i) if self.complex().dim == 3 && (i as usize) < self.complex().edges3.len() => {
                Some(self.complex().edge3_point(i, fr.u))
            }
            EntityRef::Edge(i) if (i as usize) < self.complex().edges.len() => {
                let ed = &self.complex().edges[i as usize];
                let cx = self.complex();
                Some(cx.wv(cx.edge_point(i, ed.t0 + fr.u.clamp(0.0, 1.0) * (ed.t1 - ed.t0))))
            }
            EntityRef::Face(i) if (i as usize) < self.complex().faces.len() => {
                let cx = self.complex();
                let f = &cx.faces[i as usize];
                let p = match f.geom {
                    crate::complex::FaceGeom::Lateral { edge, z0, z1 } => {
                        let ed = &cx.stack.as_ref()?.overlay.edges[edge as usize];
                        let q = cx.stack.as_ref()?.overlay.edge_point(edge, ed.t0 + f.u_clamp(fr.u) * (ed.t1 - ed.t0));
                        [q[0], q[1], z0 + f.u_clamp(fr.v) * (z1 - z0)]
                    }
                    crate::complex::FaceGeom::Quadric(b) => b.point_at(fr.u.clamp(0.0, 1.0), fr.v.clamp(0.0, 1.0)),
                    crate::complex::FaceGeom::Planar { plane, cell } => {
                        let b = cx.planar[plane as usize].cx.regions[cell as usize].bbox;
                        cx.plane_lift(
                            plane,
                            [
                                b.lo[0] + fr.u.clamp(0.0, 1.0) * (b.hi[0] - b.lo[0]),
                                b.lo[1] + fr.v.clamp(0.0, 1.0) * (b.hi[1] - b.lo[1]),
                            ],
                        )
                    }
                    crate::complex::FaceGeom::Cap { cell, z } => {
                        let b = cx.stack.as_ref()?.overlay.regions[cell as usize].bbox;
                        [
                            b.lo[0] + f.u_clamp(fr.u) * (b.hi[0] - b.lo[0]),
                            b.lo[1] + f.u_clamp(fr.v) * (b.hi[1] - b.lo[1]),
                            z,
                        ]
                    }
                };
                // `p` is in the internal frame (identity for polyhedral models)
                Some(cx.facet_project(i, cx.wv(p)))
            }
            EntityRef::Vertex(i) => self.complex().vertices.get(i as usize).map(|v| v.p),
            _ => None,
        }
    }

    // ---- topology signature ------------------------------------------------------------------------------------------------

    /// Hash of the combinatorial structure only: regions (kinds), edges (curve kinds, end vertices, sides, origins and names),
    /// vertices, loops and their incidence. Independent of parameter values; equal signatures mean identical id sets, in identical
    /// order, with the same meaning (names and provenance included).
    pub fn topology_signature(&self) -> u64 {
        let cx = self.complex();
        let mut h = Hasher64::new();
        h.u64(cx.dim as u64)
            .u64(cx.regions.len() as u64)
            .u64(cx.edges.len() as u64)
            .u64(cx.vertices.len() as u64)
            .u64(cx.loops.len() as u64);
        for r in &cx.regions {
            h.u64(r.kind as u64)
                .u64(r.feature as u64)
                .u64(r.parent.map_or(u64::MAX, |p| p as u64))
                .u64(r.loops.len() as u64);
            h.str(&r.name);
            r.loops.iter().for_each(|&l| {
                h.u64(l as u64);
            });
        }
        for e in &cx.edges {
            h.str(crate::complex::curve_kind(&cx.curves[e.curve as usize])).u64(e.v0 as u64).u64(e.v1 as u64);
            h.u64(e.left.map_or(u64::MAX, |x| x as u64)).u64(e.right.map_or(u64::MAX, |x| x as u64)).str(&e.name);
            for p in &e.prov {
                h.u64(p.feature as u64).str(&p.instance).str(&p.role);
            }
        }
        for l in &cx.loops {
            h.u64(l.region as u64).u64((l.area < 0.0) as u64).u64(l.outer.map_or(u64::MAX, |x| x as u64));
            for &(e, f) in &l.edges {
                h.u64(e as u64 * 2 + f as u64);
            }
        }
        for r in &cx.regions {
            h.u64(r.kind as u64);
        }
        if let Some(st) = &cx.stack {
            h.u64(st.layers.len() as u64).u64(cx.faces.len() as u64).u64(cx.edges3.len() as u64);
            for e in &cx.edges3 {
                h.str(&e.name).u64(e.v0 as u64).u64(e.v1 as u64);
                e.faces.iter().for_each(|&f| {
                    h.u64(f as u64);
                });
            }
            for f in &cx.faces {
                h.str(&f.name).u64(f.first.map_or(u64::MAX, |x| x as u64)).u64(f.second.map_or(u64::MAX, |x| x as u64));
                match f.geom {
                    crate::complex::FaceGeom::Lateral { edge, .. } => h.u64(1).u64(edge as u64),
                    crate::complex::FaceGeom::Cap { cell, .. } => h.u64(2).u64(cell as u64),
                    crate::complex::FaceGeom::Quadric(b) => h.u64(3).str(b.kind()),
                    crate::complex::FaceGeom::Planar { plane, cell } => h.u64(4).u64(plane as u64).u64(cell as u64),
                };
                for p in &f.prov {
                    h.u64(p.feature as u64).str(&p.instance).str(&p.role);
                }
            }
            for l in &st.layers {
                h.u64(l.region_ids.len() as u64);
                l.region_of.iter().for_each(|r| {
                    h.u64(r.map_or(u64::MAX, |x| x as u64));
                });
                h.u64(Evaluated::partition_sig(&l.partition));
            }
            h.u64(Evaluated::partition_sig(&st.overlay));
        }
        h.finish()
    }

    /// Hash of the structure of a 2D partition (used for the layers and the overlay of a stack).
    fn partition_sig(cx: &crate::complex::Complex) -> u64 {
        let mut h = Hasher64::new();
        h.u64(cx.regions.len() as u64).u64(cx.edges.len() as u64).u64(cx.vertices.len() as u64);
        for e in &cx.edges {
            h.str(crate::complex::curve_kind(&cx.curves[e.curve as usize])).u64(e.v0 as u64).u64(e.v1 as u64);
            h.u64(e.left.map_or(u64::MAX, |x| x as u64)).u64(e.right.map_or(u64::MAX, |x| x as u64));
        }
        for l in &cx.loops {
            h.u64(l.region as u64).u64(l.edges.len() as u64);
        }
        h.finish()
    }

    /// Hash of the xy geometry of the overlay of a layer stack (curves and edge ranges): unchanged when only z extents move.
    pub fn overlay_geometry_hash(&self) -> Option<u64> {
        let ov = &self.complex.as_ref()?.stack.as_ref()?.overlay;
        let mut h = Hasher64::new();
        for e in &ov.edges {
            ov.curves[e.curve as usize].hash_into(&mut h);
            h.f64(e.t0).f64(e.t1);
        }
        // curved bodies are not part of the overlay: any change to one is more than a z change
        for p in self.scene.as_ref()?.prims.iter().filter(|p| p.is_curved_solid()) {
            h.u64(p.key);
        }
        Some(h.finish())
    }
}

impl crate::complex::Complex {
    pub fn discretize_edge(&self, e: u32, o: DiscretizeOpts) -> Vec<(f64, V3)> {
        let (c, t0, t1) = self.curve_of(e);
        let sweep = t1 - t0;
        let len = self.edges[e as usize].length;
        let mut n = (len / o.max_len).ceil().max(1.0) as usize;
        match c {
            Curve::Line { .. } => {}
            Curve::Circle { r, .. } => {
                let half = (1.0 - o.max_dev / r).clamp(-1.0, 1.0).acos();
                n = n
                    .max((sweep / (2.0 * half).max(1e-3)).ceil() as usize)
                    .max((sweep / std::f64::consts::FRAC_PI_2).ceil() as usize);
            }
            Curve::Ellipse { .. } => n = n.max((sweep / (std::f64::consts::PI / 8.0)).ceil() as usize),
        }
        let mut ts: Vec<f64> = (0..=n).map(|k| t0 + sweep * k as f64 / n as f64).collect();
        if let Curve::Ellipse { .. } = c {
            // refine where the sampled chord deviation is too large
            let mut k = 0;
            while k + 1 < ts.len() && ts.len() < 100_000 {
                if chord_dev(c, ts[k], ts[k + 1]) > o.max_dev || dist2(c.eval(ts[k]), c.eval(ts[k + 1])) > o.max_len {
                    ts.insert(k + 1, 0.5 * (ts[k] + ts[k + 1]));
                } else {
                    k += 1;
                }
            }
        }
        let last = ts.len() - 1;
        ts.iter().enumerate().map(|(i, &t)| (t, to3(c.eval(if i == last { t1 } else { t })))).collect()
    }
}

fn chord_dev(c: &Curve, ta: f64, tb: f64) -> f64 {
    match c {
        Curve::Line { .. } => 0.0,
        Curve::Circle { r, .. } => r * (1.0 - (0.5 * (tb - ta).abs().min(std::f64::consts::TAU)).cos()),
        Curve::Ellipse { .. } => {
            let (pa, pb) = (c.eval(ta), c.eval(tb));
            let d = sub2(pb, pa);
            let l = norm2(d);
            if l == 0.0 {
                // a closed chord: distance from its end point to the farthest sample
                return (1..32).map(|k| dist2(c.eval(ta + (tb - ta) * k as f64 / 32.0), pa)).fold(0.0, f64::max) * 0.5;
            }
            (1..32)
                .map(|k| {
                    let p = c.eval(ta + (tb - ta) * k as f64 / 32.0);
                    cross2(d, sub2(p, pa)).abs() / l
                })
                .fold(0.0, f64::max)
        }
    }
}
