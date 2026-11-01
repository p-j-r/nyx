//! 2D primitives in the world frame: a local shape (rectangle, rounded rectangle, circle, ellipse, polygon), the similarity
//! transform that placed it, its directed boundary pieces with role names (interior always on the left), an exact signed
//! distance in the local frame, and a bounding box. Curves are *not* interned here; the arrangement does that once the
//! model tolerance is known.

use crate::curve::Curve;
use crate::error::GeomResult;
use crate::math::*;

/// Local-frame shape. Polygons are counter-clockwise with `names[i]` the role of the edge `pts[i] -> pts[i+1]`.
#[derive(Clone, Debug, PartialEq)]
pub enum LShape {
    Rect { lo: V2, hi: V2 },
    Rounded { lo: V2, hi: V2, r: f64 },
    Circle { c: V2, r: f64 },
    Ellipse { c: V2, a: f64, b: f64 },
    Poly { pts: Vec<V2>, names: Vec<String> },
}

/// A directed boundary piece in world coordinates; the primitive's interior is on the left of the travel direction.
#[derive(Clone, Debug, PartialEq)]
pub enum DGeom {
    Line {
        p0: V2,
        p1: V2,
    },
    /// counter-clockwise arc from angle `a0` over `sweep` (`TAU`: the full circle)
    Arc {
        c: V2,
        r: f64,
        a0: f64,
        sweep: f64,
    },
    /// full ellipse
    Ell {
        c: V2,
        a: f64,
        b: f64,
        phi: f64,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct DSeg {
    pub geom: DGeom,
    pub role: String,
}

/// A convex polyhedron: planar faces whose rings are counter-clockwise seen from outside.
#[derive(Clone, Debug)]
pub struct Poly3 {
    pub faces: Vec<PFace>,
}

#[derive(Clone, Debug)]
pub struct PFace {
    pub ring: Vec<V3>,
    pub role: String,
    /// outward unit normal and plane offset (`n . p = d`)
    pub n: V3,
    pub d: f64,
}

impl Poly3 {
    /// Builds a polyhedron from face rings in any orientation; each ring is turned to point away from the centroid.
    pub fn new(rings: Vec<(Vec<V3>, String)>) -> Poly3 {
        let pts: Vec<V3> = rings.iter().flat_map(|r| r.0.iter().copied()).collect();
        let c = mul3(pts.iter().fold([0.0; 3], |a, p| add3(a, *p)), 1.0 / pts.len() as f64);
        let faces = rings
            .into_iter()
            .map(|(mut ring, role)| {
                let newell = |r: &[V3]| {
                    let mut n = [0.0; 3];
                    for i in 0..r.len() {
                        let (a, b) = (r[i], r[(i + 1) % r.len()]);
                        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
                        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
                        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
                    }
                    unit3(n)
                };
                let mut n = newell(&ring);
                let fc = mul3(ring.iter().fold([0.0; 3], |a, p| add3(a, *p)), 1.0 / ring.len() as f64);
                if dot3(n, sub3(fc, c)) < 0.0 {
                    ring.reverse();
                    n = mul3(n, -1.0);
                }
                let d = ring.iter().map(|p| dot3(n, *p)).sum::<f64>() / ring.len() as f64;
                PFace { ring, role, n, d }
            })
            .collect();
        Poly3 { faces }
    }
    pub fn bbox(&self) -> Aabb {
        let mut b = Aabb::EMPTY;
        self.faces.iter().flat_map(|f| f.ring.iter()).for_each(|p| b.grow(*p));
        b
    }
    /// Largest signed plane distance: exact inside, a lower bound of the distance outside.
    pub fn sdf(&self, p: V3) -> f64 {
        self.faces.iter().map(|f| dot3(f.n, p) - f.d).fold(f64::NEG_INFINITY, f64::max)
    }
    pub fn key(&self, h: &mut Hasher64) {
        for f in &self.faces {
            h.str(&f.role);
            f.ring.iter().for_each(|p| {
                h.f64(p[0]).f64(p[1]).f64(p[2]);
            });
        }
    }
}

/// A curved solid and the curved faces it is made of: sphere and torus (one closed face, no edges), a cone frustum along z, and the
/// faces of cones and revolved profiles (`ConeSide`, `CylSide`, `Disk`, `Annulus`). Exact volume, area, distance, projection and
/// parametrisation (`u`: azimuth, `v`: polar / minor angle / height / radius, as fractions).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Body {
    Sphere {
        c: V3,
        r: f64,
    },
    Torus {
        c: V3,
        big: f64,
        small: f64,
    },
    /// frustum from the base centre `c` up by `h`, radius `r0` at the bottom and `r1` at the top (either may be 0)
    Cone {
        c: V3,
        r0: f64,
        r1: f64,
        h: f64,
    },
    /// `flip`: the outward normal of the solid points towards the axis
    ConeSide {
        c: V3,
        r0: f64,
        r1: f64,
        h: f64,
        flip: bool,
    },
    Disk {
        c: V3,
        r: f64,
        up: bool,
    },
    /// a cylinder wall from `c` up by `h`; `inward`: the outward normal of the solid points towards the axis
    CylSide {
        c: V3,
        r: f64,
        h: f64,
        inward: bool,
    },
    /// a flat ring between two radii at the height of `c`
    Annulus {
        c: V3,
        r_in: f64,
        r_out: f64,
        up: bool,
    },
}

impl Body {
    pub fn kind(&self) -> &'static str {
        match self {
            Body::Sphere { .. } => "sphere",
            Body::Torus { .. } => "torus",
            Body::Cone { .. } | Body::ConeSide { .. } => "cone",
            Body::Disk { .. } | Body::Annulus { .. } => "plane",
            Body::CylSide { .. } => "cylinder",
        }
    }
    pub fn volume(&self) -> f64 {
        let pi = std::f64::consts::PI;
        match self {
            Body::Sphere { r, .. } => 4.0 / 3.0 * pi * r.powi(3),
            Body::Torus { big, small, .. } => 2.0 * pi.powi(2) * big * small * small,
            Body::Cone { r0, r1, h, .. } => pi * h * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0,
            _ => 0.0,
        }
    }
    pub fn area(&self) -> f64 {
        let pi = std::f64::consts::PI;
        match self {
            Body::Sphere { r, .. } => 4.0 * pi * r * r,
            Body::Torus { big, small, .. } => 4.0 * pi.powi(2) * big * small,
            Body::Cone { .. } => 0.0,
            Body::ConeSide { r0, r1, h, .. } => pi * (r0 + r1) * h.hypot(r0 - r1),
            Body::Disk { r, .. } => pi * r * r,
            Body::CylSide { r, h, .. } => TAU * r * h,
            Body::Annulus { r_in, r_out, .. } => pi * (r_out * r_out - r_in * r_in),
        }
    }
    /// Total length of the circle edges around a face.
    pub fn perimeter(&self) -> f64 {
        match self {
            Body::ConeSide { r0, r1, .. } => TAU * (r0 + r1),
            Body::Disk { r, .. } => TAU * r,
            Body::CylSide { r, .. } => 2.0 * TAU * r,
            Body::Annulus { r_in, r_out, .. } => TAU * (r_in + r_out),
            _ => 0.0,
        }
    }
    pub fn bbox(&self) -> Aabb {
        match self {
            Body::Sphere { c, r } => Aabb::new(sub3(*c, [*r; 3]), add3(*c, [*r; 3])),
            Body::Torus { c, big, small } => {
                let e = big + small;
                Aabb::new([c[0] - e, c[1] - e, c[2] - small], [c[0] + e, c[1] + e, c[2] + small])
            }
            Body::Cone { c, r0, r1, h } | Body::ConeSide { c, r0, r1, h, .. } => {
                let r = r0.max(*r1);
                Aabb::new([c[0] - r, c[1] - r, c[2]], [c[0] + r, c[1] + r, c[2] + h])
            }
            Body::Disk { c, r, .. } => Aabb::new([c[0] - r, c[1] - r, c[2]], [c[0] + r, c[1] + r, c[2]]),
            Body::CylSide { c, r, h, .. } => Aabb::new([c[0] - r, c[1] - r, c[2]], [c[0] + r, c[1] + r, c[2] + h]),
            Body::Annulus { c, r_out, .. } => {
                Aabb::new([c[0] - r_out, c[1] - r_out, c[2]], [c[0] + r_out, c[1] + r_out, c[2]])
            }
        }
    }
    /// Distance from `q = (rho, z)` to the segment `a -> b`, and the nearest point.
    pub(crate) fn seg(q: V2, a: V2, b: V2) -> (f64, V2) {
        let d = sub2(b, a);
        let l2 = dot2(d, d);
        let t = if l2 > 0.0 { (dot2(sub2(q, a), d) / l2).clamp(0.0, 1.0) } else { 0.0 };
        let n = add2(a, mul2(d, t));
        (dist2(q, n), n)
    }
    /// Signed distance of a solid (negative inside); infinite for face-only variants.
    pub fn sdf(&self, p: V3) -> f64 {
        match self {
            Body::Sphere { c, r } => norm3(sub3(p, *c)) - r,
            Body::Torus { c, big, small } => {
                let d = sub3(p, *c);
                (((d[0].hypot(d[1]) - big).powi(2) + d[2] * d[2]).sqrt()) - small
            }
            Body::Cone { c, r0, r1, h } => {
                let d = sub3(p, *c);
                let q = [d[0].hypot(d[1]), d[2]];
                let dist = [(0.0, 0.0, *r0, 0.0), (*r0, 0.0, *r1, *h), (*r1, *h, 0.0, *h)]
                    .iter()
                    .map(|&(x0, y0, x1, y1)| Body::seg(q, [x0, y0], [x1, y1]).0)
                    .fold(f64::INFINITY, f64::min);
                let inside = q[1] >= 0.0 && q[1] <= *h && q[0] <= r0 + (r1 - r0) * q[1] / h;
                if inside {
                    -dist
                } else {
                    dist
                }
            }
            _ => f64::INFINITY,
        }
    }
    /// Closest point of the surface (of the face, for the face-only variants).
    pub fn project(&self, p: V3) -> V3 {
        match self {
            Body::Sphere { c, r } => {
                let d = sub3(p, *c);
                let n = norm3(d);
                if n == 0.0 {
                    add3(*c, [*r, 0.0, 0.0])
                } else {
                    add3(*c, mul3(d, r / n))
                }
            }
            Body::Torus { c, big, small } => {
                let d = sub3(p, *c);
                let rho = d[0].hypot(d[1]);
                let (ux, uy) = if rho == 0.0 { (1.0, 0.0) } else { (d[0] / rho, d[1] / rho) };
                let ring = [c[0] + big * ux, c[1] + big * uy, c[2]];
                let e = sub3(p, ring);
                let n = norm3(e);
                if n == 0.0 {
                    add3(ring, [small * ux, small * uy, 0.0])
                } else {
                    add3(ring, mul3(e, small / n))
                }
            }
            Body::Cone { c, r0, r1, h } => {
                let parts = [
                    Body::ConeSide { c: *c, r0: *r0, r1: *r1, h: *h, flip: false },
                    Body::Disk { c: *c, r: *r0, up: false },
                    Body::Disk { c: add3(*c, [0.0, 0.0, *h]), r: *r1, up: true },
                ];
                parts
                    .iter()
                    .map(|b| b.project(p))
                    .min_by(|a, b| norm3(sub3(*a, p)).total_cmp(&norm3(sub3(*b, p))))
                    .unwrap()
            }
            Body::ConeSide { c, r0, r1, h, .. } => {
                let d = sub3(p, *c);
                let rho = d[0].hypot(d[1]);
                let (ux, uy) = if rho == 0.0 { (1.0, 0.0) } else { (d[0] / rho, d[1] / rho) };
                let (_, n) = Body::seg([rho, d[2]], [*r0, 0.0], [*r1, *h]);
                [c[0] + n[0] * ux, c[1] + n[0] * uy, c[2] + n[1]]
            }
            Body::Disk { c, r, .. } => {
                let d = sub3(p, *c);
                let rho = d[0].hypot(d[1]);
                let k = if rho > *r { r / rho } else { 1.0 };
                [c[0] + d[0] * k, c[1] + d[1] * k, c[2]]
            }
            Body::CylSide { c, r, h, .. } => {
                let d = sub3(p, *c);
                let rho = d[0].hypot(d[1]);
                let (ux, uy) = if rho == 0.0 { (1.0, 0.0) } else { (d[0] / rho, d[1] / rho) };
                [c[0] + r * ux, c[1] + r * uy, c[2] + d[2].clamp(0.0, *h)]
            }
            Body::Annulus { c, r_in, r_out, .. } => {
                let d = sub3(p, *c);
                let rho = d[0].hypot(d[1]);
                let (ux, uy) = if rho == 0.0 { (1.0, 0.0) } else { (d[0] / rho, d[1] / rho) };
                let k = rho.clamp(*r_in, *r_out);
                [c[0] + k * ux, c[1] + k * uy, c[2]]
            }
        }
    }
    /// Outward unit normal (outward from the solid) at the surface point nearest `p`.
    pub fn normal(&self, p: V3) -> V3 {
        match self {
            Body::Sphere { c, .. } => unit3(sub3(self.project(p), *c)),
            Body::Torus { c, big, .. } => {
                let q = self.project(p);
                let d = sub3(q, *c);
                let rho = d[0].hypot(d[1]).max(1e-300);
                unit3(sub3(q, [c[0] + big * d[0] / rho, c[1] + big * d[1] / rho, c[2]]))
            }
            Body::Cone { .. } => unit3(sub3(p, self.project(p))),
            Body::ConeSide { c, r0, r1, h, flip } => {
                let d = sub3(self.project(p), *c);
                let rho = d[0].hypot(d[1]);
                let (ux, uy) = if rho == 0.0 { (1.0, 0.0) } else { (d[0] / rho, d[1] / rho) };
                let s = h.hypot(r0 - r1);
                let n = [h / s * ux, h / s * uy, (r0 - r1) / s];
                if *flip {
                    mul3(n, -1.0)
                } else {
                    n
                }
            }
            Body::CylSide { c, inward, .. } => {
                let d = sub3(self.project(p), *c);
                let rho = d[0].hypot(d[1]).max(1e-300);
                let s = if *inward { -1.0 } else { 1.0 };
                [s * d[0] / rho, s * d[1] / rho, 0.0]
            }
            Body::Disk { up, .. } | Body::Annulus { up, .. } => [0.0, 0.0, if *up { 1.0 } else { -1.0 }],
        }
    }
    /// Point at the fractions `(u, v)`.
    pub fn point_at(&self, u: f64, v: f64) -> V3 {
        let (a, b) = (u * TAU, v * TAU);
        match self {
            Body::Sphere { c, r } => {
                // v runs from the south pole (0) to the north pole (1)
                let th = v * std::f64::consts::PI;
                [c[0] + r * th.sin() * a.cos(), c[1] + r * th.sin() * a.sin(), c[2] - r * th.cos()]
            }
            Body::Torus { c, big, small } => {
                let rho = big + small * b.cos();
                [c[0] + rho * a.cos(), c[1] + rho * a.sin(), c[2] + small * b.sin()]
            }
            Body::Cone { c, r0, r1, h } | Body::ConeSide { c, r0, r1, h, .. } => {
                let rho = r0 + (r1 - r0) * v;
                [c[0] + rho * a.cos(), c[1] + rho * a.sin(), c[2] + h * v]
            }
            Body::Disk { c, r, .. } => [c[0] + v * r * a.cos(), c[1] + v * r * a.sin(), c[2]],
            Body::CylSide { c, r, h, .. } => [c[0] + r * a.cos(), c[1] + r * a.sin(), c[2] + h * v],
            Body::Annulus { c, r_in, r_out, .. } => {
                let rho = r_in + (r_out - r_in) * v;
                [c[0] + rho * a.cos(), c[1] + rho * a.sin(), c[2]]
            }
        }
    }
    /// Fractions of the surface point nearest `p`.
    pub fn fractions(&self, p: V3) -> (f64, f64) {
        let q = self.project(p);
        let turn = |d: V3| (d[1].atan2(d[0]) / TAU).rem_euclid(1.0);
        match self {
            Body::Sphere { c, r } => {
                let d = sub3(q, *c);
                (turn(d), ((-d[2] / r).clamp(-1.0, 1.0).acos()) / std::f64::consts::PI)
            }
            Body::Torus { c, big, .. } => {
                let d = sub3(q, *c);
                let rho = d[0].hypot(d[1]);
                (turn(d), (d[2].atan2(rho - big) / TAU).rem_euclid(1.0))
            }
            Body::Cone { c, h, .. } | Body::ConeSide { c, h, .. } | Body::CylSide { c, h, .. } => {
                let d = sub3(q, *c);
                (turn(d), (d[2] / h).clamp(0.0, 1.0))
            }
            Body::Disk { c, r, .. } => {
                let d = sub3(q, *c);
                (turn(d), if *r > 0.0 { (d[0].hypot(d[1]) / r).clamp(0.0, 1.0) } else { 0.0 })
            }
            Body::Annulus { c, r_in, r_out, .. } => {
                let d = sub3(q, *c);
                let w = r_out - r_in;
                (turn(d), if w > 0.0 { ((d[0].hypot(d[1]) - r_in) / w).clamp(0.0, 1.0) } else { 0.0 })
            }
        }
    }
    pub fn radius_of_curvature(&self) -> f64 {
        match self {
            Body::Sphere { r, .. } => *r,
            Body::Torus { small, .. } => *small,
            Body::Cone { r0, r1, .. } | Body::ConeSide { r0, r1, .. } => r0.max(*r1),
            Body::Disk { r, .. } | Body::CylSide { r, .. } => *r,
            Body::Annulus { r_out, .. } => *r_out,
        }
    }
    pub fn key(&self, h: &mut Hasher64) {
        match self {
            Body::Sphere { c, r } => {
                h.u64(7).f64(c[0]).f64(c[1]).f64(c[2]).f64(*r);
            }
            Body::Torus { c, big, small } => {
                h.u64(8).f64(c[0]).f64(c[1]).f64(c[2]).f64(*big).f64(*small);
            }
            Body::Cone { c, r0, r1, h: hh } | Body::ConeSide { c, r0, r1, h: hh, .. } => {
                h.u64(9).f64(c[0]).f64(c[1]).f64(c[2]).f64(*r0).f64(*r1).f64(*hh);
            }
            Body::Disk { c, r, up } => {
                h.u64(10).f64(c[0]).f64(c[1]).f64(c[2]).f64(*r).u64(*up as u64);
            }
            Body::CylSide { c, r, h: hh, inward } => {
                h.u64(11).f64(c[0]).f64(c[1]).f64(c[2]).f64(*r).f64(*hh).u64(*inward as u64);
            }
            Body::Annulus { c, r_in, r_out, up } => {
                h.u64(12).f64(c[0]).f64(c[1]).f64(c[2]).f64(*r_in).f64(*r_out).u64(*up as u64);
            }
        }
    }
}

/// A solid of revolution: a closed counter-clockwise profile in the `(rho, z)` half plane (`rho >= 0`) turned once about the axis through
/// `c`. Exact volume (Pappus), distance and faces: one per profile segment off the axis (cone side, cylinder wall or flat ring), with a
/// circle edge at every profile corner off the axis.
#[derive(Clone, Debug, PartialEq)]
pub struct Revolved {
    pub c: V3,
    pub pts: Vec<V2>,
    /// role of the segment `pts[i] -> pts[i+1]`
    pub roles: Vec<String>,
}

/// One face of a revolved solid: its role, geometry, and the profile corner (index) at its end.
pub struct RevFace {
    pub role: String,
    pub body: Body,
    /// index of the profile segment
    pub seg: usize,
}

impl Revolved {
    pub fn volume(&self) -> f64 {
        let n = self.pts.len();
        let s: f64 = (0..n)
            .map(|i| {
                let (a, b) = (self.pts[i], self.pts[(i + 1) % n]);
                (a[0] + b[0]) * (a[0] * b[1] - b[0] * a[1])
            })
            .sum();
        TAU / 6.0 * s
    }
    pub fn bbox(&self) -> Aabb {
        let rm = self.pts.iter().map(|p| p[0]).fold(0.0, f64::max);
        let (z0, z1) = self.pts.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, p| (a.0.min(p[1]), a.1.max(p[1])));
        Aabb::new([self.c[0] - rm, self.c[1] - rm, self.c[2] + z0], [self.c[0] + rm, self.c[1] + rm, self.c[2] + z1])
    }
    fn on_axis(a: V2, b: V2, tol: f64) -> bool {
        a[0] <= tol && b[0] <= tol
    }
    /// Signed distance (negative inside) to the surface, ignoring the segments on the axis.
    pub fn sdf(&self, p: V3) -> f64 {
        let d = sub3(p, self.c);
        let q = [d[0].hypot(d[1]), d[2]];
        let n = self.pts.len();
        let (mut dist, mut inside) = (f64::INFINITY, false);
        for i in 0..n {
            let (a, b) = (self.pts[i], self.pts[(i + 1) % n]);
            if (a[1] > q[1]) != (b[1] > q[1]) && q[0] < (b[0] - a[0]) * (q[1] - a[1]) / (b[1] - a[1]) + a[0] {
                inside = !inside;
            }
            if !Revolved::on_axis(a, b, 0.0) {
                dist = dist.min(Body::seg(q, a, b).0);
            }
        }
        if inside {
            -dist
        } else {
            dist
        }
    }
    /// The faces (profile segments off the axis).
    pub fn faces(&self) -> Vec<RevFace> {
        let n = self.pts.len();
        let mut out = vec![];
        for i in 0..n {
            let (a, b) = (self.pts[i], self.pts[(i + 1) % n]);
            if Revolved::on_axis(a, b, 0.0) {
                continue;
            }
            let (dr, dz) = (b[0] - a[0], b[1] - a[1]);
            let role = format!("swept:{}", self.roles[i]);
            let at = |z: f64| [self.c[0], self.c[1], self.c[2] + z];
            let body = if dz == 0.0 {
                // flat ring; the solid is on the left of the travel direction, so the outward normal is +z when moving to smaller rho
                Body::Annulus { c: at(a[1]), r_in: a[0].min(b[0]), r_out: a[0].max(b[0]), up: dr < 0.0 }
            } else if dr == 0.0 {
                Body::CylSide { c: at(a[1].min(b[1])), r: a[0], h: dz.abs(), inward: dz < 0.0 }
            } else {
                let (lo, hi) = if a[1] < b[1] { (a, b) } else { (b, a) };
                Body::ConeSide { c: at(lo[1]), r0: lo[0], r1: hi[0], h: hi[1] - lo[1], flip: dz < 0.0 }
            };
            out.push(RevFace { role, body, seg: i });
        }
        out
    }
    /// A point on the surface away from the edges (the middle of the longest segment off the axis).
    pub fn surface_point(&self) -> V3 {
        let n = self.pts.len();
        let (a, b) = (0..n)
            .map(|i| (self.pts[i], self.pts[(i + 1) % n]))
            .filter(|(a, b)| !Revolved::on_axis(*a, *b, 0.0))
            .max_by(|x, y| dist2(x.0, x.1).total_cmp(&dist2(y.0, y.1)))
            .unwrap_or((self.pts[0], self.pts[1]));
        [self.c[0] + 0.5 * (a[0] + b[0]), self.c[1], self.c[2] + 0.5 * (a[1] + b[1])]
    }
    pub fn key(&self, h: &mut Hasher64) {
        h.u64(13).f64(self.c[0]).f64(self.c[1]).f64(self.c[2]);
        for (p, r) in self.pts.iter().zip(&self.roles) {
            h.f64(p[0]).f64(p[1]).str(r);
        }
    }
}

/// How the boundary pieces are named in 3D: a rectangle read as a box, a circle read as a cylinder, any profile extruded, or plain 2D.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Style {
    Plain,
    Box,
    Cylinder,
    Extrude,
}

#[derive(Clone, Debug)]
pub struct Prim {
    /// a solid of revolution (then `segs` is empty and the xy shape is only a placeholder)
    pub rev: Option<Revolved>,
    /// a polyhedral solid (then `segs` is empty and the xy shape is only a placeholder)
    pub poly: Option<Poly3>,
    /// a curved body (then `segs` is empty and the xy shape is only a placeholder)
    pub body: Option<Body>,
    /// z interval of the solid (a prism over this profile); `(0, 0)` for 2D models
    pub z: (f64, f64),
    /// roles of the bottom and top faces
    pub caps: [String; 2],
    pub feature: u32,
    /// instance name: feature name plus pattern indices (`hole[3]`)
    pub name: String,
    pub shape: LShape,
    pub xf: Xf2,
    pub segs: Vec<DSeg>,
    pub bbox: Aabb,
    /// helper primitives (cut half-planes) do not count towards the model size
    pub aux: bool,
    pub key: u64,
}

fn poly_area(p: &[V2]) -> f64 {
    0.5 * (0..p.len()).map(|i| cross2(p[i], p[(i + 1) % p.len()])).sum::<f64>()
}

impl LShape {
    /// Local segments (counter-clockwise) as `(geometry, role)`.
    fn local_segs(&self) -> Vec<(DGeom, String)> {
        let line = |p0: V2, p1: V2, r: &str| (DGeom::Line { p0, p1 }, r.to_string());
        let arc = |c: V2, r: f64, a0: f64, name: &str| {
            (DGeom::Arc { c, r, a0, sweep: std::f64::consts::FRAC_PI_2 }, name.to_string())
        };
        match self {
            LShape::Rect { lo, hi } => vec![
                line([lo[0], lo[1]], [hi[0], lo[1]], "bottom"),
                line([hi[0], lo[1]], [hi[0], hi[1]], "right"),
                line([hi[0], hi[1]], [lo[0], hi[1]], "top"),
                line([lo[0], hi[1]], [lo[0], lo[1]], "left"),
            ],
            LShape::Rounded { lo, hi, r } => {
                let r = *r;
                let (x0, y0, x1, y1) = (lo[0], lo[1], hi[0], hi[1]);
                let mut v = vec![];
                let h = std::f64::consts::FRAC_PI_2;
                let push_line = |v: &mut Vec<(DGeom, String)>, p0: V2, p1: V2, role: &str| {
                    if dist2(p0, p1) > 0.0 {
                        v.push(line(p0, p1, role));
                    }
                };
                push_line(&mut v, [x0 + r, y0], [x1 - r, y0], "bottom");
                v.push(arc([x1 - r, y0 + r], r, -h, "c_br"));
                push_line(&mut v, [x1, y0 + r], [x1, y1 - r], "right");
                v.push(arc([x1 - r, y1 - r], r, 0.0, "c_tr"));
                push_line(&mut v, [x1 - r, y1], [x0 + r, y1], "top");
                v.push(arc([x0 + r, y1 - r], r, h, "c_tl"));
                push_line(&mut v, [x0, y1 - r], [x0, y0 + r], "left");
                v.push(arc([x0 + r, y0 + r], r, 2.0 * h, "c_bl"));
                v
            }
            LShape::Circle { c, r } => vec![(DGeom::Arc { c: *c, r: *r, a0: 0.0, sweep: TAU }, "edge".into())],
            LShape::Ellipse { c, a, b } => vec![(DGeom::Ell { c: *c, a: *a, b: *b, phi: 0.0 }, "edge".into())],
            LShape::Poly { pts, names } => (0..pts.len())
                .map(|i| (DGeom::Line { p0: pts[i], p1: pts[(i + 1) % pts.len()] }, names[i].clone()))
                .collect(),
        }
    }

    /// Signed distance in the local frame (negative inside).
    fn sdf(&self, p: V2) -> f64 {
        let boxd = |lo: V2, hi: V2| {
            let (c, h) =
                ([0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])], [0.5 * (hi[0] - lo[0]), 0.5 * (hi[1] - lo[1])]);
            let q = [(p[0] - c[0]).abs() - h[0], (p[1] - c[1]).abs() - h[1]];
            q[0].max(0.0).hypot(q[1].max(0.0)) + q[0].max(q[1]).min(0.0)
        };
        match self {
            LShape::Rect { lo, hi } => boxd(*lo, *hi),
            LShape::Rounded { lo, hi, r } => boxd([lo[0] + r, lo[1] + r], [hi[0] - r, hi[1] - r]) - r,
            LShape::Circle { c, r } => dist2(p, *c) - r,
            LShape::Ellipse { c, a, b } => Curve::Ellipse { c: *c, a: *a, b: *b, phi: 0.0 }.sdf(p).unwrap_or(0.0),
            LShape::Poly { pts, .. } => {
                let n = pts.len();
                let (mut d2, mut inside) = (f64::INFINITY, false);
                for i in 0..n {
                    let (a, b) = (pts[i], pts[(i + 1) % n]);
                    let ab = sub2(b, a);
                    let t = (dot2(sub2(p, a), ab) / dot2(ab, ab)).clamp(0.0, 1.0);
                    let q = add2(a, mul2(ab, t));
                    d2 = d2.min(dot2(sub2(p, q), sub2(p, q)));
                    if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0] {
                        inside = !inside;
                    }
                }
                if inside {
                    -d2.sqrt()
                } else {
                    d2.sqrt()
                }
            }
        }
    }

    fn hash(&self, h: &mut Hasher64) {
        match self {
            LShape::Rect { lo, hi } => {
                h.u64(1).f64(lo[0]).f64(lo[1]).f64(hi[0]).f64(hi[1]);
            }
            LShape::Rounded { lo, hi, r } => {
                h.u64(2).f64(lo[0]).f64(lo[1]).f64(hi[0]).f64(hi[1]).f64(*r);
            }
            LShape::Circle { c, r } => {
                h.u64(3).f64(c[0]).f64(c[1]).f64(*r);
            }
            LShape::Ellipse { c, a, b } => {
                h.u64(4).f64(c[0]).f64(c[1]).f64(*a).f64(*b);
            }
            LShape::Poly { pts, names } => {
                h.u64(5);
                for (p, n) in pts.iter().zip(names) {
                    h.f64(p[0]).f64(p[1]).str(n);
                }
            }
        }
    }
}

impl DGeom {
    pub fn curve_range(&self) -> (Curve, f64, f64) {
        match self {
            DGeom::Line { p0, p1 } => {
                let (c, rev) = crate::curve::canonical_line(*p0, sub2(*p1, *p0));
                let (a, b) = (c.closest(*p0).0, c.closest(*p1).0);
                let _ = rev;
                (c, a.min(b), a.max(b))
            }
            DGeom::Arc { c, r, a0, sweep } => {
                let a = a0.rem_euclid(TAU);
                (Curve::Circle { c: *c, r: *r }, a, a + sweep)
            }
            DGeom::Ell { c, a, b, phi } => (crate::curve::canonical_ellipse(*c, *a, *b, *phi, 0.0), 0.0, TAU),
        }
    }
    fn map(&self, xf: &Xf2) -> DGeom {
        match self {
            DGeom::Line { p0, p1 } => DGeom::Line { p0: xf.apply(*p0), p1: xf.apply(*p1) },
            DGeom::Arc { c, r, a0, sweep } => {
                let c2 = xf.apply(*c);
                let at = |a: f64| xf.apply([c[0] + r * a.cos(), c[1] + r * a.sin()]);
                let start = if xf.mirrored() { at(a0 + sweep) } else { at(*a0) };
                let v = sub2(start, c2);
                DGeom::Arc { c: c2, r: r * xf.scale(), a0: v[1].atan2(v[0]), sweep: *sweep }
            }
            DGeom::Ell { c, a, b, .. } => {
                let e1 = xf.apply_dir([1.0, 0.0]);
                DGeom::Ell { c: xf.apply(*c), a: a * xf.scale(), b: b * xf.scale(), phi: e1[1].atan2(e1[0]) }
            }
        }
    }
}

impl Prim {
    pub fn new(feature: u32, name: String, shape: LShape, xf: Xf2, aux: bool) -> GeomResult<Prim> {
        Prim::new_styled(feature, name, shape, xf, aux, Style::Plain, (0.0, 0.0), 2)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_styled(
        feature: u32,
        name: String,
        shape: LShape,
        xf: Xf2,
        aux: bool,
        style: Style,
        z: (f64, f64),
        axis: usize,
    ) -> GeomResult<Prim> {
        let mut segs: Vec<DSeg> =
            shape.local_segs().into_iter().map(|(g, role)| DSeg { geom: g.map(&xf), role }).collect();
        if xf.mirrored() {
            // reflection reverses the travel direction of lines (arcs were re-anchored above); restore counter-clockwise order
            for s in &mut segs {
                if let DGeom::Line { p0, p1 } = &mut s.geom {
                    std::mem::swap(p0, p1);
                }
            }
            segs.reverse();
        }
        let names = ["x", "y", "z"];
        for s in &mut segs {
            s.role = match style {
                Style::Plain => continue,
                Style::Extrude => format!("side:{}", s.role),
                Style::Cylinder => "side".to_string(),
                Style::Box => {
                    // internal x and y are the world axes `axis_perm` assigns; roles carry the world names
                    let pm = crate::scene::axis_perm(axis);
                    let (nx, ny) = (names[pm[0]], names[pm[1]]);
                    match s.role.as_str() {
                        "left" => format!("{nx}min"),
                        "right" => format!("{nx}max"),
                        "bottom" => format!("{ny}min"),
                        _ => format!("{ny}max"),
                    }
                }
            };
        }
        let caps = if style == Style::Box {
            [format!("{}min", names[axis]), format!("{}max", names[axis])]
        } else {
            [crate::scene::cap_role(axis, false), crate::scene::cap_role(axis, true)]
        };
        let mut bbox = Aabb::EMPTY;
        for s in &segs {
            let (c, a, b) = s.geom.curve_range();
            bbox = bbox.union(&c.bbox(a, b));
        }
        let mut h = Hasher64::new();
        shape.hash(&mut h);
        h.u64(axis as u64)
            .u64(style as u64)
            .f64(z.0)
            .f64(z.1)
            .f64(xf.m[0][0])
            .f64(xf.m[0][1])
            .f64(xf.m[1][0])
            .f64(xf.m[1][1])
            .f64(xf.t[0])
            .f64(xf.t[1])
            .str(&name)
            .u64(aux as u64)
            .u64(feature as u64);
        Ok(Prim {
            rev: None,
            poly: None,
            body: None,
            z,
            caps,
            feature,
            name,
            shape,
            xf,
            segs,
            bbox,
            aux,
            key: h.finish(),
        })
    }

    pub fn new_body(feature: u32, name: String, body: Body) -> Prim {
        let b = body.bbox();
        let mut h = Hasher64::new();
        body.key(&mut h);
        h.str(&name).u64(feature as u64);
        Prim {
            rev: None,
            poly: None,
            body: Some(body),
            z: (b.lo[2], b.hi[2]),
            caps: ["surface".to_string(), "surface".to_string()],
            feature,
            name,
            shape: LShape::Circle { c: [b.lo[0], b.lo[1]], r: 0.0 },
            xf: Xf2::ID,
            segs: vec![],
            bbox: Aabb::new([b.lo[0], b.lo[1], 0.0], [b.hi[0], b.hi[1], 0.0]),
            aux: false,
            key: h.finish(),
        }
    }

    pub fn new_revolved(feature: u32, name: String, rev: Revolved) -> Prim {
        let b = rev.bbox();
        let mut h = Hasher64::new();
        rev.key(&mut h);
        h.str(&name).u64(feature as u64);
        Prim {
            rev: Some(rev),
            poly: None,
            body: None,
            z: (b.lo[2], b.hi[2]),
            caps: ["surface".to_string(), "surface".to_string()],
            feature,
            name,
            shape: LShape::Circle { c: [b.lo[0], b.lo[1]], r: 0.0 },
            xf: Xf2::ID,
            segs: vec![],
            bbox: Aabb::new([b.lo[0], b.lo[1], 0.0], [b.hi[0], b.hi[1], 0.0]),
            aux: false,
            key: h.finish(),
        }
    }

    /// True for solids that take no part in the planar / prismatic arrangement (bodies and solids of revolution).
    pub fn is_curved_solid(&self) -> bool {
        self.body.is_some() || self.rev.is_some()
    }

    pub fn new_poly(feature: u32, name: String, poly: Poly3, aux: bool) -> Prim {
        let b = poly.bbox();
        let mut h = Hasher64::new();
        poly.key(&mut h);
        h.str(&name).u64(feature as u64).u64(aux as u64);
        Prim {
            rev: None,
            poly: Some(poly),
            body: None,
            z: (b.lo[2], b.hi[2]),
            caps: ["bottom".to_string(), "top".to_string()],
            feature,
            name,
            shape: LShape::Circle { c: [b.lo[0], b.lo[1]], r: 0.0 },
            xf: Xf2::ID,
            segs: vec![],
            bbox: Aabb::new([b.lo[0], b.lo[1], 0.0], [b.hi[0], b.hi[1], 0.0]),
            aux,
            key: h.finish(),
        }
    }

    /// Bounding box including the z interval.
    pub fn bbox3(&self) -> Aabb {
        Aabb { lo: [self.bbox.lo[0], self.bbox.lo[1], self.z.0], hi: [self.bbox.hi[0], self.bbox.hi[1], self.z.1] }
    }

    /// Signed distance estimate of the prism in 3D (negative inside).
    pub fn sdf3(&self, p: V3, reach: f64) -> f64 {
        if let Some(b) = &self.body {
            return b.sdf(p);
        }
        if let Some(q) = &self.poly {
            return q.sdf(p);
        }
        if let Some(r) = &self.rev {
            return r.sdf(p);
        }
        let dxy = self.sdf_near([p[0], p[1]], reach);
        if self.z.0 == self.z.1 {
            return dxy;
        }
        let dz = (self.z.0 - p[2]).max(p[2] - self.z.1);
        if dxy > 0.0 || dz > 0.0 {
            dxy.max(0.0).hypot(dz.max(0.0))
        } else {
            dxy.max(dz)
        }
    }

    /// Signed distance estimate in world units (negative inside); exact for polygons, rectangles and circles.
    pub fn sdf(&self, p: V2) -> f64 {
        self.xf.scale() * self.shape.sdf(self.xf.inverse().apply(p))
    }

    /// Signed distance, skipping the evaluation when `p` is farther than `reach` from the bounding box.
    pub fn sdf_near(&self, p: V2, reach: f64) -> f64 {
        if self.is_curved_solid() || self.poly.is_some() {
            return f64::INFINITY; // curved bodies take no part in the planar arrangement
        }
        let b = &self.bbox;
        let d = (b.lo[0] - p[0]).max(p[0] - b.hi[0]).max(b.lo[1] - p[1]).max(p[1] - b.hi[1]);
        if d > reach {
            d
        } else {
            self.sdf(p)
        }
    }
}

/// Counter-clockwise copy of a polygon with the edge names reordered to match.
pub fn ccw_polygon(pts: &[V2], names: &[String]) -> (Vec<V2>, Vec<String>) {
    if poly_area(pts) >= 0.0 {
        return (pts.to_vec(), names.to_vec());
    }
    let n = pts.len();
    let q: Vec<V2> = pts.iter().rev().copied().collect();
    // edge j of the reversed polygon joins pts[n-1-j] and pts[n-2-j], which is the original edge n-2-j
    let nm = (0..n).map(|j| names[(2 * n - 2 - j) % n].clone()).collect();
    (q, nm)
}

pub fn polygon_area(p: &[V2]) -> f64 {
    poly_area(p)
}
