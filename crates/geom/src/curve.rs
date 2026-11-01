//! Planar curves of the 2D kernel: lines, circles and ellipses in one canonical parametrisation each (so that coincident
//! curves compare equal), with evaluation, tangent, curvature, closest point, bounding boxes, length, the Green-theorem
//! boundary integral used for exact areas, and curve-curve intersection with tangency / near-coincidence reporting.
//!
//! Parametrisation: `Line { o, d }` has `o` the point nearest the origin and `d` a sign-canonical unit direction, `t = (p - o) . d`.
//! `Circle` and `Ellipse` run counter-clockwise from the major axis, `p = c + a cos(t) e1 + b sin(t) e2`.

use crate::math::*;

#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Line { o: V2, d: V2 },
    Circle { c: V2, r: f64 },
    Ellipse { c: V2, a: f64, b: f64, phi: f64 },
}

/// One intersection point with its parameter on each curve.
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    pub p: V2,
    pub ta: f64,
    pub tb: f64,
}

/// Result of intersecting two curves. `graze` lists tangent contacts (grazing within the tolerance band), `near` is set when
/// the curves are nearly but not exactly coincident (parallel lines, concentric circles, almost equal ellipses).
#[derive(Clone, Debug, Default)]
pub struct Inter {
    pub hits: Vec<Hit>,
    pub graze: Vec<V2>,
    pub near: bool,
}

fn ell_axes(phi: f64) -> (V2, V2) {
    let (s, c) = phi.sin_cos();
    ([c, s], [-s, c])
}

impl Curve {
    pub fn periodic(&self) -> bool {
        !matches!(self, Curve::Line { .. })
    }

    pub fn eval(&self, t: f64) -> V2 {
        match self {
            Curve::Line { o, d } => [o[0] + t * d[0], o[1] + t * d[1]],
            Curve::Circle { c, r } => [c[0] + r * t.cos(), c[1] + r * t.sin()],
            Curve::Ellipse { c, a, b, phi } => {
                let (e1, e2) = ell_axes(*phi);
                let (ca, sb) = (a * t.cos(), b * t.sin());
                [c[0] + ca * e1[0] + sb * e2[0], c[1] + ca * e1[1] + sb * e2[1]]
            }
        }
    }

    /// Velocity `dp/dt`.
    pub fn deriv(&self, t: f64) -> V2 {
        match self {
            Curve::Line { d, .. } => *d,
            Curve::Circle { r, .. } => [-r * t.sin(), r * t.cos()],
            Curve::Ellipse { a, b, phi, .. } => {
                let (e1, e2) = ell_axes(*phi);
                let (ca, sb) = (-a * t.sin(), b * t.cos());
                [ca * e1[0] + sb * e2[0], ca * e1[1] + sb * e2[1]]
            }
        }
    }

    /// Unit tangent in the direction of increasing `t`.
    pub fn tangent(&self, t: f64) -> V2 {
        let v = self.deriv(t);
        mul2(v, 1.0 / norm2(v))
    }

    /// Signed curvature (positive when the curve turns left in the direction of increasing `t`).
    pub fn curvature(&self, t: f64) -> f64 {
        match self {
            Curve::Line { .. } => 0.0,
            Curve::Circle { r, .. } => 1.0 / r,
            Curve::Ellipse { a, b, .. } => {
                let q = a * a * t.sin().powi(2) + b * b * t.cos().powi(2);
                a * b / (q * q.sqrt())
            }
        }
    }

    /// Closest point on the whole curve (parameter unrestricted; periodic parameters lie in `[0, 2 pi)`).
    pub fn closest(&self, p: V2) -> (f64, V2) {
        let t = match self {
            Curve::Line { o, d } => dot2(sub2(p, *o), *d),
            Curve::Circle { c, .. } => {
                let v = sub2(p, *c);
                if v == [0.0, 0.0] {
                    0.0
                } else {
                    v[1].atan2(v[0]).rem_euclid(TAU)
                }
            }
            Curve::Ellipse { c, a, b, phi } => {
                let (e1, e2) = ell_axes(*phi);
                let v = sub2(p, *c);
                let q = [dot2(v, e1), dot2(v, e2)];
                let f = |t: f64| (b * b - a * a) * t.sin() * t.cos() + a * q[0] * t.sin() - b * q[1] * t.cos();
                let df = |t: f64| (b * b - a * a) * (2.0 * t).cos() + a * q[0] * t.cos() + b * q[1] * t.sin();
                let dist2 = |t: f64| (a * t.cos() - q[0]).powi(2) + (b * t.sin() - q[1]).powi(2);
                let mut best = (f64::INFINITY, 0.0);
                for k in 0..5 {
                    let mut t = match k {
                        0 => (q[1] / b).atan2(q[0] / a),
                        k => (k as f64 - 1.0) * std::f64::consts::FRAC_PI_2,
                    };
                    for _ in 0..40 {
                        let d = df(t);
                        if d == 0.0 {
                            break;
                        }
                        let step = (f(t) / d).clamp(-0.5, 0.5);
                        t -= step;
                        if step.abs() < 1e-15 {
                            break;
                        }
                    }
                    let d2 = dist2(t);
                    if d2 < best.0 {
                        best = (d2, t);
                    }
                }
                best.1.rem_euclid(TAU)
            }
        };
        (t, self.eval(t))
    }

    /// Signed distance estimate for closed curves (negative inside); `None` for lines.
    pub fn sdf(&self, p: V2) -> Option<f64> {
        match self {
            Curve::Line { .. } => None,
            Curve::Circle { c, r } => Some(dist2(p, *c) - r),
            Curve::Ellipse { c, a, b, phi } => {
                let (e1, e2) = ell_axes(*phi);
                let v = sub2(p, *c);
                let (qx, qy) = (dot2(v, e1) / a, dot2(v, e2) / b);
                let qn = qx.hypot(qy);
                if qn == 0.0 {
                    return Some(-a.min(*b));
                }
                let g = [qx / (a * qn), qy / (b * qn)];
                Some((qn - 1.0) / norm2(g))
            }
        }
    }

    /// True if `t` (any representative) lies in the segment `[a, b]`, with a parameter slack `eps`.
    pub fn in_range(&self, t: f64, a: f64, b: f64, eps: f64) -> bool {
        if self.periodic() {
            if b - a >= TAU - eps {
                return true;
            }
            let u = wrap_from(t, a);
            u <= b + eps || u >= a + TAU - eps
        } else {
            t >= a - eps && t <= b + eps
        }
    }

    /// Parameter of a point of the curve relative to a segment start (`a`), wrapped for periodic curves.
    pub fn param_in(&self, t: f64, a: f64) -> f64 {
        if self.periodic() {
            wrap_from(t, a)
        } else {
            t
        }
    }

    /// Bounding box of the segment `[a, b]`.
    pub fn bbox(&self, a: f64, b: f64) -> Aabb {
        let mut bb = Aabb::EMPTY;
        for t in [a, b] {
            let p = self.eval(t);
            bb.grow([p[0], p[1], 0.0]);
        }
        if let Curve::Circle { .. } | Curve::Ellipse { .. } = self {
            let (aa, bb_, phi) = match self {
                Curve::Circle { r, .. } => (*r, *r, 0.0),
                Curve::Ellipse { a, b, phi, .. } => (*a, *b, *phi),
                _ => unreachable!(),
            };
            let tx = (-bb_ * phi.sin()).atan2(aa * phi.cos());
            let ty = (bb_ * phi.cos()).atan2(aa * phi.sin());
            for t in [tx, tx + std::f64::consts::PI, ty, ty + std::f64::consts::PI] {
                if self.in_range(t, a, b, 0.0) {
                    let p = self.eval(t);
                    bb.grow([p[0], p[1], 0.0]);
                }
            }
        }
        bb
    }

    /// Arc length of the segment `[a, b]`.
    pub fn length(&self, a: f64, b: f64) -> f64 {
        match self {
            Curve::Line { .. } => (b - a).abs(),
            Curve::Circle { r, .. } => r * (b - a).abs(),
            Curve::Ellipse { .. } => {
                let panels = 8;
                let (xs, ws) = gauss_legendre(16);
                let h = (b - a) / panels as f64;
                let mut s = 0.0;
                for k in 0..panels {
                    let lo = a + h * k as f64;
                    for (x, w) in xs.iter().zip(&ws) {
                        s += w * 0.5 * h * norm2(self.deriv(lo + 0.5 * h * (1.0 + x)));
                    }
                }
                s.abs()
            }
        }
    }

    /// `integral of (x dy - y dx)` along `[a, b]` about the point `o` (twice the signed area swept); Green's theorem gives exact
    /// areas. `o` should be near the model so that small regions far from the global origin do not lose digits.
    pub fn green(&self, a: f64, b: f64, o: V2) -> f64 {
        match self {
            Curve::Line { .. } => cross2(sub2(self.eval(a), o), sub2(self.eval(b), o)),
            Curve::Circle { c, r } => {
                let c = sub2(*c, o);
                r * (c[0] * (b.sin() - a.sin()) - c[1] * (b.cos() - a.cos())) + r * r * (b - a)
            }
            Curve::Ellipse { c, a: sa, b: sb, phi } => {
                let c = &sub2(*c, o);
                let (e1, e2) = ell_axes(*phi);
                let (qa, qb) = ([a.cos(), a.sin()], [b.cos(), b.sin()]);
                let dq = sub2(qb, qa);
                let aq = [sa * e1[0] * dq[0] + sb * e2[0] * dq[1], sa * e1[1] * dq[0] + sb * e2[1] * dq[1]];
                cross2(*c, aq) + sa * sb * (b - a)
            }
        }
    }

    /// `integral of x dA`-type first moments are not needed; centroid of the boundary curve for naming and ordering.
    pub fn mid_point(&self, a: f64, b: f64) -> V2 {
        self.eval(0.5 * (a + b))
    }

    /// Stable content hash of the geometry (cache keys).
    pub fn hash_into(&self, h: &mut Hasher64) {
        match self {
            Curve::Line { o, d } => {
                h.u64(1).f64(o[0]).f64(o[1]).f64(d[0]).f64(d[1]);
            }
            Curve::Circle { c, r } => {
                h.u64(2).f64(c[0]).f64(c[1]).f64(*r);
            }
            Curve::Ellipse { c, a, b, phi } => {
                h.u64(3).f64(c[0]).f64(c[1]).f64(*a).f64(*b).f64(*phi);
            }
        }
    }
}

/// Gauss-Legendre nodes and weights on `[-1, 1]` by Newton iteration on the Legendre polynomial.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (mut xs, mut ws) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        xs[i] = x;
        ws[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (xs, ws)
}

/// Canonical line through `p` with direction `dir` and whether `dir` is reversed with respect to the canonical direction.
pub fn canonical_line(p: V2, dir: V2) -> (Curve, bool) {
    let d = mul2(dir, 1.0 / norm2(dir));
    let rev = d[0] < -1e-12 || (d[0].abs() <= 1e-12 && d[1] < 0.0);
    let d = if rev { mul2(d, -1.0) } else { d };
    let s = dot2(p, d);
    (Curve::Line { o: sub2(p, mul2(d, s)), d }, rev)
}

/// Canonical ellipse; collapses to a circle when the semi-axes agree to `tol`.
pub fn canonical_ellipse(c: V2, a: f64, b: f64, phi: f64, tol: f64) -> Curve {
    if (a - b).abs() <= tol {
        return Curve::Circle { c, r: 0.5 * (a + b) };
    }
    let (a, b, phi) = if a >= b { (a, b, phi) } else { (b, a, phi + std::f64::consts::FRAC_PI_2) };
    Curve::Ellipse { c, a, b, phi: phi.rem_euclid(std::f64::consts::PI) }
}

/// Intersection of two curves. `k` is the grazing band in units of `tol`: contacts closer than `k * tol` are tangent.
/// Segment ranges are used only to bound the numeric search for ellipses.
pub fn intersect(ca: &Curve, cb: &Curve, ra: (f64, f64), rb: (f64, f64), tol: f64, k: f64) -> Inter {
    let band = k * tol;
    let mut out = Inter::default();
    let push = |ca: &Curve, cb: &Curve, p: V2, out: &mut Inter| {
        out.hits.push(Hit { p, ta: ca.closest(p).0, tb: cb.closest(p).0 });
    };
    match (ca, cb) {
        (Curve::Line { o: o1, d: d1 }, Curve::Line { o: o2, d: d2 }) => {
            let den = cross2(*d1, *d2);
            if den.abs() < 1e-10 {
                // parallel: distinct lines (equal ones are interned together)
                if cross2(*d1, sub2(*o2, *o1)).abs() <= band {
                    out.near = true;
                }
            } else {
                let t = cross2(sub2(*o2, *o1), *d2) / den;
                push(ca, cb, add2(*o1, mul2(*d1, t)), &mut out);
            }
        }
        (Curve::Line { .. }, Curve::Circle { c, r }) | (Curve::Circle { c, r }, Curve::Line { .. }) => {
            let l = if matches!(ca, Curve::Line { .. }) { ca } else { cb };
            let Curve::Line { o, d } = l else { unreachable!() };
            let s = dot2(sub2(*c, *o), *d);
            let dist = cross2(*d, sub2(*c, *o)).abs();
            if dist > r + band {
            } else if (dist - r).abs() <= band {
                out.graze.push(add2(*o, mul2(*d, s)));
            } else {
                let h = (r * r - dist * dist).sqrt();
                for sg in [-1.0, 1.0] {
                    push(ca, cb, add2(*o, mul2(*d, s + sg * h)), &mut out);
                }
            }
        }
        (Curve::Circle { c: c1, r: r1 }, Curve::Circle { c: c2, r: r2 }) => {
            let d = dist2(*c1, *c2);
            if d <= band {
                if (r1 - r2).abs() <= band {
                    out.near = true;
                }
            } else if (d - (r1 + r2)).abs() <= band || (d - (r1 - r2).abs()).abs() <= band {
                let u = mul2(sub2(*c2, *c1), 1.0 / d);
                let toward = if (d - (r1 + r2)).abs() <= band || r1 > r2 { *r1 } else { -*r1 };
                out.graze.push(add2(*c1, mul2(u, toward)));
            } else if d < r1 + r2 && d > (r1 - r2).abs() {
                let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
                let h = (r1 * r1 - x * x).max(0.0).sqrt();
                let u = mul2(sub2(*c2, *c1), 1.0 / d);
                let m = add2(*c1, mul2(u, x));
                for sg in [-1.0, 1.0] {
                    push(ca, cb, add2(m, mul2(perp(u), sg * h)), &mut out);
                }
            }
        }
        _ => numeric(ca, cb, ra, rb, tol, k, &mut out),
    }
    out
}

/// Intersection when an ellipse is involved: sample one curve over its segment, bracket sign changes of the other's implicit
/// distance, bisect, and flag near-zero minima without a sign change as grazing.
fn numeric(ca: &Curve, cb: &Curve, ra: (f64, f64), rb: (f64, f64), tol: f64, k: f64, out: &mut Inter) {
    // sample the curve that is not an ellipse over its range when possible, else the first ellipse over its range
    let (sa, ra_, sb) = match (ca, cb) {
        (Curve::Ellipse { .. }, Curve::Ellipse { .. }) => (ca, ra, cb),
        (Curve::Ellipse { .. }, _) => (cb, rb, ca),
        _ => (ca, ra, cb),
    };
    let signed = |p: V2| -> f64 {
        match sb {
            Curve::Line { o, d } => cross2(*d, sub2(p, *o)),
            _ => sb.sdf(p).unwrap_or(0.0),
        }
    };
    let span = if sa.periodic() { TAU } else { ra_.1 - ra_.0 };
    let n = 720usize;
    let t0 = ra_.0;
    let at = |i: usize| t0 + span * i as f64 / n as f64;
    let f: Vec<f64> = (0..=n).map(|i| signed(sa.eval(at(i)))).collect();
    let size = sa.bbox(t0, t0 + span).diagonal();
    let record = |t: f64, out: &mut Inter| {
        let p = sa.eval(t);
        out.hits.push(Hit { p, ta: ca.closest(p).0, tb: cb.closest(p).0 });
    };
    let band = k * tol;
    for i in 0..n {
        let (f0, f1) = (f[i], f[i + 1]);
        if f0 == 0.0 || (f0 < 0.0) != (f1 < 0.0) {
            let (mut lo, mut hi, mut flo) = (at(i), at(i + 1), f0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let fm = signed(sa.eval(mid));
                if (fm < 0.0) == (flo < 0.0) {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
            }
            record(0.5 * (lo + hi), out);
        } else if i > 0
            && f0.abs() <= f[i - 1].abs()
            && f0.abs() <= f1.abs()
            && (f[i - 1] < 0.0) == (f0 < 0.0)
            && f0.abs() < 1e-3 * size
        {
            // local minimum of |f| without a sign change: refine by golden section and test the band
            let (mut lo, mut hi) = (at(i - 1), at(i + 1));
            let g = (5f64.sqrt() - 1.0) / 2.0;
            for _ in 0..80 {
                let (x1, x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
                if signed(sa.eval(x1)).abs() < signed(sa.eval(x2)).abs() {
                    hi = x2;
                } else {
                    lo = x1;
                }
            }
            let tm = 0.5 * (lo + hi);
            if signed(sa.eval(tm)).abs() <= band {
                out.graze.push(sa.eval(tm));
            }
        }
    }
    // an ellipse against an equal-ish ellipse never crosses: report near-coincidence from the sampled distances
    if out.hits.is_empty() && out.graze.is_empty() {
        let mx = f.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if mx <= band {
            out.near = true;
        }
    }
}
