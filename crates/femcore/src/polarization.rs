//! Polarization analysis of optical modes: Stokes parameters, Poincaré sphere and a classification of the mode.
//!
//! A mode propagating along z has a transverse field (Ex, Ey) that is complex at every point of the cross-section;
//! in time its tip draws an ellipse. The four Stokes parameters describe that ellipse:
//!
//! * `S0 = |Ex|^2 + |Ey|^2` intensity
//! * `S1 = |Ex|^2 - |Ey|^2` x versus y
//! * `S2 = 2 Re(conj(Ex) Ey)` +45 degrees versus -45 degrees
//! * `S3 = -+2 Im(conj(Ex) Ey)` circular part
//!
//! The sign of S3 depends on the time convention of the field, so it is defined by what can be seen: `S3 > 0` means
//! that the E-vector turns counter-clockwise in the (x, y) plane as time advances (x to the right, y up). Whether
//! that is called "left" or "right" circular differs between optics and IEEE usage; this module only says
//! counter-clockwise (CCW) or clockwise (CW). Divided by S0 the parameters give a point (s1, s2, s3) on the unit
//! Poincaré sphere: the equator is linear polarization, the poles are circular, latitude is twice the ellipticity
//! angle chi and longitude twice the orientation angle psi.
//!
//! One point of a mode is always fully polarized. The whole mode is not: `analyze` adds the Stokes parameters of all
//! samples (weighted by area), and the degree of polarization `P = |(S1, S2, S3)| / S0` of the sum is below 1 when
//! the polarization changes across the cross-section (a hybrid mode). The Stokes parameters do not change when the
//! eigensolver multiplies the whole mode by a complex number.
//!
//! The module works on plain samples (complex components and a weight) so that it does not depend on how a solver
//! stores its fields. `p1_vertex_weights` makes area weights for fields given at the vertices of a triangle mesh. Only
//! the transverse components enter the Stokes parameters; the longitudinal part is reported as a fraction.
//!
//! Axes: x is the lateral (in-plane) direction, y the normal direction and z the propagation direction, as in the
//! model frame of the solvers. "Quasi-TE" means x-dominated and "quasi-TM" y-dominated. For a ring (r, z) mesh pass
//! E_r as `ex` and E_z as `ey` and E_theta as `ez`, and decide how that maps to TE and TM.

use std::fmt::Write;

/// Complex number as [re, im].
pub type Cx = [f64; 2];

fn conj_mul(a: Cx, b: Cx) -> Cx {
    [a[0] * b[0] + a[1] * b[1], a[0] * b[1] - a[1] * b[0]]
}

fn abs2(a: Cx) -> f64 {
    a[0] * a[0] + a[1] * a[1]
}

/// Time dependence of the complex field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum TimeConvention {
    /// `E(t) = Re(E exp(+i w t))` (engineering). A wave along +z is `exp(-i kz z)`, which is what the solvers in nyx use for
    /// the Floquet and kz factors; this is the default.
    #[default]
    PlusJ,
    /// `E(t) = Re(E exp(-i w t))` (physics)
    MinusI,
}

impl TimeConvention {
    /// +1 when the Stokes S3 is `+2 Im(conj(Ex) Ey)` for counter-clockwise rotation.
    fn s3_sign(self) -> f64 {
        match self {
            TimeConvention::MinusI => 1.0,
            TimeConvention::PlusJ => -1.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Stokes {
    pub s0: f64,
    pub s1: f64,
    pub s2: f64,
    pub s3: f64,
}

impl Stokes {
    pub fn from_fields(ex: Cx, ey: Cx, conv: TimeConvention) -> Stokes {
        let cross = conj_mul(ex, ey);
        Stokes { s0: abs2(ex) + abs2(ey), s1: abs2(ex) - abs2(ey), s2: 2.0 * cross[0], s3: conv.s3_sign() * 2.0 * cross[1] }
    }

    fn add_scaled(&mut self, o: &Stokes, w: f64) {
        self.s0 += w * o.s0;
        self.s1 += w * o.s1;
        self.s2 += w * o.s2;
        self.s3 += w * o.s3;
    }

    /// (s1, s2, s3) divided by S0; zeros when there is no intensity. A point of a single sample has unit length.
    pub fn normalized(&self) -> [f64; 3] {
        if self.s0 > 0.0 {
            [self.s1 / self.s0, self.s2 / self.s0, self.s3 / self.s0]
        } else {
            [0.0; 3]
        }
    }

    /// Degree of polarization in [0, 1].
    pub fn dop(&self) -> f64 {
        let s = self.normalized();
        (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt().min(1.0)
    }

    /// Orientation of the major axis in degrees, in (-90, 90], measured from x towards y.
    pub fn psi_deg(&self) -> f64 {
        0.5 * self.s2.atan2(self.s1).to_degrees()
    }

    /// Ellipticity angle in degrees, in [-45, 45]: 0 linear, +45 CCW circular, -45 CW circular. Uses the direction of
    /// the polarized part, so it is meaningful for partially polarized sums too.
    pub fn chi_deg(&self) -> f64 {
        let p = (self.s1 * self.s1 + self.s2 * self.s2 + self.s3 * self.s3).sqrt();
        if p > 0.0 {
            0.5 * (self.s3 / p).clamp(-1.0, 1.0).asin().to_degrees()
        } else {
            0.0
        }
    }
}

/// Angle in radians between two points of the Poincaré sphere (pi for orthogonal polarizations).
pub fn sphere_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let na = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let nb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
    if na == 0.0 || nb == 0.0 {
        return 0.0;
    }
    ((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (na * nb)).clamp(-1.0, 1.0).acos()
}

/// The field of a mode at one point of the cross-section and the area (or volume share) that point stands for.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    pub ex: Cx,
    pub ey: Cx,
    /// longitudinal component (zero when not available)
    pub ez: Cx,
    pub w: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PolClass {
    /// linear, along x
    QuasiTe,
    /// linear, along y
    QuasiTm,
    /// linear, along +45 degrees
    LinearPlus45,
    /// linear, along -45 degrees
    LinearMinus45,
    Elliptical { ccw: bool },
    Circular { ccw: bool },
    /// polarization changes across the mode (degree of polarization below the limit)
    Hybrid,
}

impl PolClass {
    pub fn label(&self) -> &'static str {
        match self {
            PolClass::QuasiTe => "quasi-TE (x)",
            PolClass::QuasiTm => "quasi-TM (z)",
            PolClass::LinearPlus45 => "linear +45 deg",
            PolClass::LinearMinus45 => "linear -45 deg",
            PolClass::Elliptical { ccw: true } => "elliptical CCW",
            PolClass::Elliptical { ccw: false } => "elliptical CW",
            PolClass::Circular { ccw: true } => "circular CCW",
            PolClass::Circular { ccw: false } => "circular CW",
            PolClass::Hybrid => "hybrid",
        }
    }
}

/// Limits of the classification, all on the integrated Stokes vector of the mode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClassOptions {
    /// below this degree of polarization the mode is hybrid
    pub p_min: f64,
    /// `|s3| / P` at or above this is circular
    pub circ_min: f64,
    /// `|s3| / P` at or above this (and below `circ_min`) is elliptical; below it the mode is linear
    pub ell_min: f64,
    pub time: TimeConvention,
}

impl Default for ClassOptions {
    fn default() -> Self {
        ClassOptions { p_min: 0.8, circ_min: 0.9, ell_min: 0.3, time: TimeConvention::default() }
    }
}

#[derive(Clone, Debug)]
pub struct ModeReport {
    /// Stokes parameters summed over all samples with their weights
    pub integrated: Stokes,
    /// (s1, s2, s3) of the integrated vector; its length is the degree of polarization
    pub sphere: [f64; 3],
    pub dop: f64,
    pub chi_deg: f64,
    pub psi_deg: f64,
    /// `integral |Ex|^2 / integral (|Ex|^2 + |Ey|^2)`; equals `(1 + s1) / 2` of the integrated vector
    pub te_fraction: f64,
    /// share of the field along the propagation axis (y in the wafer frame): `integral |E_y|^2 / integral |E|^2`
    pub long_fraction: f64,
    pub class: PolClass,
}

impl ModeReport {
    pub fn summary_line(&self) -> String {
        format!(
            "s = ({:+.3}, {:+.3}, {:+.3})  P = {:.3}  chi = {:+.1} deg  psi = {:+.1} deg  TE = {:.3}  Ey = {:.3}  {}",
            self.sphere[0], self.sphere[1], self.sphere[2], self.dop, self.chi_deg, self.psi_deg, self.te_fraction, self.long_fraction, self.class.label()
        )
    }
}

/// Stokes parameters of every sample, in the order of `samples` (for plots on the sphere and maps of the cross-section).
pub fn local_stokes(samples: &[Sample], conv: TimeConvention) -> Vec<Stokes> {
    samples.iter().map(|s| Stokes::from_fields(s.ex, s.ey, conv)).collect()
}

pub fn classify(integrated: &Stokes, opts: &ClassOptions) -> PolClass {
    let s = integrated.normalized();
    let p = integrated.dop();
    if p < opts.p_min {
        return PolClass::Hybrid;
    }
    let circ = s[2] / p;
    if circ.abs() >= opts.circ_min {
        return PolClass::Circular { ccw: s[2] > 0.0 };
    }
    if circ.abs() >= opts.ell_min {
        return PolClass::Elliptical { ccw: s[2] > 0.0 };
    }
    if s[0].abs() >= s[1].abs() {
        if s[0] >= 0.0 { PolClass::QuasiTe } else { PolClass::QuasiTm }
    } else if s[1] >= 0.0 {
        PolClass::LinearPlus45
    } else {
        PolClass::LinearMinus45
    }
}

/// The report of one mode; `None` if the weighted intensity is zero or not finite.
pub fn analyze(samples: &[Sample], opts: &ClassOptions) -> Option<ModeReport> {
    let mut sum = Stokes::default();
    let mut ez2 = 0.0;
    let mut ex2 = 0.0;
    let mut ey2 = 0.0;
    for s in samples {
        sum.add_scaled(&Stokes::from_fields(s.ex, s.ey, opts.time), s.w);
        ex2 += s.w * abs2(s.ex);
        ey2 += s.w * abs2(s.ey);
        ez2 += s.w * abs2(s.ez);
    }
    if !sum.s0.is_finite() || sum.s0 <= 0.0 {
        return None;
    }
    Some(ModeReport {
        integrated: sum,
        sphere: sum.normalized(),
        dop: sum.dop(),
        chi_deg: sum.chi_deg(),
        psi_deg: sum.psi_deg(),
        te_fraction: ex2 / (ex2 + ey2),
        long_fraction: ez2 / (ex2 + ey2 + ez2),
        class: classify(&sum, opts),
    })
}

/// Groups of modes whose effective indices differ by at most `tol` (consecutive in sorted order of the real part). The
/// eigensolver may return any mixture of degenerate modes, so their polarization points are not meaningful on their own.
pub fn degenerate_groups(neff: &[(f64, f64)], tol: f64) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..neff.len()).collect();
    order.sort_by(|&a, &b| neff[a].0.total_cmp(&neff[b].0));
    let mut groups: Vec<Vec<usize>> = vec![];
    let mut cur: Vec<usize> = vec![];
    for &i in &order {
        let near = cur.last().is_some_and(|&j| ((neff[i].0 - neff[j].0).powi(2) + (neff[i].1 - neff[j].1).powi(2)).sqrt() <= tol);
        if near {
            cur.push(i);
        } else {
            if cur.len() > 1 {
                groups.push(std::mem::take(&mut cur));
            }
            cur = vec![i];
        }
    }
    if cur.len() > 1 {
        groups.push(cur);
    }
    for g in &mut groups {
        g.sort_unstable();
    }
    groups
}

/// Area weight of every vertex of a 2D triangle mesh (a third of each adjacent triangle), skipping triangles whose
/// label is in `skip` (for example PML domains). `points` holds 2 numbers per vertex and `cells` 3 vertex numbers per
/// triangle, as in `SimplexMesh`. For fields stored at vertices (P1) this is the lumped integral.
pub fn p1_vertex_weights(points: &[f64], cells: &[u32], cell_labels: &[i32], skip: &[i32]) -> Vec<f64> {
    let mut w = vec![0.0; points.len() / 2];
    for (t, c) in cells.chunks(3).enumerate() {
        if cell_labels.get(t).is_some_and(|l| skip.contains(l)) {
            continue;
        }
        let p = |i: u32| (points[2 * i as usize], points[2 * i as usize + 1]);
        let (a, b, d) = (p(c[0]), p(c[1]), p(c[2]));
        let area = 0.5 * ((b.0 - a.0) * (d.1 - a.1) - (d.0 - a.0) * (b.1 - a.1)).abs();
        for &v in c {
            w[v as usize] += area / 3.0;
        }
    }
    w
}

/// Projection of the sphere on the page: orthographic, S3 up, seen from 35 degrees azimuth and 20 degrees elevation.
struct View {
    cx: f64,
    cy: f64,
    r: f64,
    right: [f64; 3],
    up: [f64; 3],
    cam: [f64; 3],
}

impl View {
    fn new(size: f64) -> View {
        let (az, el) = (35f64.to_radians(), 20f64.to_radians());
        View {
            cx: 0.5 * size,
            cy: 0.5 * size,
            r: 0.37 * size,
            cam: [el.cos() * az.cos(), el.cos() * az.sin(), el.sin()],
            right: [-az.sin(), az.cos(), 0.0],
            up: [-el.sin() * az.cos(), -el.sin() * az.sin(), el.cos()],
        }
    }

    /// page x, page y and depth (positive towards the viewer)
    fn proj(&self, p: [f64; 3]) -> (f64, f64, f64) {
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        (self.cx + self.r * dot(p, self.right), self.cy - self.r * dot(p, self.up), dot(p, self.cam))
    }
}

/// Draws a curve on the sphere: the part facing the viewer solid, the part behind it dashed and pale.
fn push_curve(v: &View, pts: &[[f64; 3]], back: &mut String, front: &mut String, width: f64, op: f64) {
    let mut run: Vec<(f64, f64)> = vec![];
    let mut is_front = true;
    let flush = |run: &mut Vec<(f64, f64)>, is_front: bool, back: &mut String, front: &mut String| {
        if run.len() > 1 {
            let d: Vec<String> = run.iter().map(|(x, y)| format!("{x:.1},{y:.1}")).collect();
            let (out, extra, o) = if is_front { (front, "", op) } else { (back, " stroke-dasharray=\"2 3\"", op * 0.6) };
            let _ = write!(out, "<polyline points=\"{}\" fill=\"none\" stroke=\"#3b5b82\" stroke-width=\"{width}\" stroke-opacity=\"{o:.2}\" stroke-linecap=\"round\"{extra}/>", d.join(" "));
        }
        run.clear();
    };
    for p in pts {
        let (x, y, d) = v.proj(*p);
        let f = d >= 0.0;
        if f != is_front {
            run.push((x, y));
            flush(&mut run, is_front, back, front);
            is_front = f;
        }
        run.push((x, y));
    }
    flush(&mut run, is_front, back, front);
}

/// The curves of the empty sphere in drawing order, as `(points on the unit sphere, stroke width, opacity)`: a fine grid (meridians every 30
/// degrees, parallels at +-30 and +-60 degrees), the two vertical great circles and the bold equator of the linear polarizations.
fn frame_curves() -> Vec<(Vec<[f64; 3]>, f64, f64)> {
    let great = |k: usize, t: f64| -> [f64; 3] {
        let (c, s) = (t.cos(), t.sin());
        match k {
            0 => [c, s, 0.0],
            1 => [c, 0.0, s],
            _ => [0.0, c, s],
        }
    };
    let n = 120;
    let ring = |f: &dyn Fn(f64) -> [f64; 3]| -> Vec<[f64; 3]> { (0..=n).map(|i| f(std::f64::consts::TAU * i as f64 / n as f64)).collect() };
    let mut out = vec![];
    for m in 0..6 {
        if m % 3 == 0 {
            continue; // the planes s2 = 0 and s1 = 0 are drawn stronger below
        }
        let phi = std::f64::consts::PI * m as f64 / 6.0;
        out.push((ring(&|t| [t.cos() * phi.cos(), t.cos() * phi.sin(), t.sin()]), 0.5, 0.28));
    }
    for lat in [-60.0f64, -30.0, 30.0, 60.0] {
        let (z, rr) = (lat.to_radians().sin(), lat.to_radians().cos());
        out.push((ring(&|t| [rr * t.cos(), rr * t.sin(), z]), 0.5, 0.28));
    }
    out.push((ring(&|t| great(1, t)), 0.8, 0.5));
    out.push((ring(&|t| great(2, t)), 0.8, 0.5));
    out.push((ring(&|t| great(0, t)), 1.4, 0.85)); // equator: linear polarizations
    out
}

/// The empty sphere: returns (defs, back layer, front layer, labels). Layers are drawn back, then the translucent ball,
/// then front, then labels, so that everything behind the ball is seen through it.
fn sphere_frame(v: &View) -> (String, String, String, String) {
    let mut back = String::new();
    let mut front = String::new();
    for (curve, width, opacity) in frame_curves() {
        push_curve(v, &curve, &mut back, &mut front, width, opacity);
    }
    // axes through the centre
    for a in 0..3 {
        let mut dir = [0.0; 3];
        dir[a] = 1.0;
        let (x0, y0, _) = v.proj(dir.map(|c| -1.12 * c));
        let (x1, y1, _) = v.proj(dir.map(|c| 1.12 * c));
        let (xc, yc, _) = v.proj([0.0; 3]);
        // half towards the viewer solid, the other half dashed
        let (_, _, d1) = v.proj(dir);
        let ((fx, fy), (bx, by)) = if d1 >= 0.0 { ((x1, y1), (x0, y0)) } else { ((x0, y0), (x1, y1)) };
        let _ = write!(front, "<line x1=\"{xc:.1}\" y1=\"{yc:.1}\" x2=\"{fx:.1}\" y2=\"{fy:.1}\" stroke=\"#3b5b82\" stroke-width=\"0.9\" stroke-opacity=\"0.7\"/>");
        let _ = write!(back, "<line x1=\"{xc:.1}\" y1=\"{yc:.1}\" x2=\"{bx:.1}\" y2=\"{by:.1}\" stroke=\"#3b5b82\" stroke-width=\"0.9\" stroke-opacity=\"0.45\" stroke-dasharray=\"2 3\"/>");
    }
    // labels
    let mut labels = String::new();
    for (p, text) in [([1.0, 0.0, 0.0], "x"), ([-1.0, 0.0, 0.0], "y"), ([0.0, 1.0, 0.0], "+45\u{b0}"), ([0.0, -1.0, 0.0], "\u{2212}45\u{b0}"), ([0.0, 0.0, 1.0], "CCW"), ([0.0, 0.0, -1.0], "CW")] {
        let (x, y, d) = v.proj(p.map(|c| 1.27 * c));
        let op = if d >= 0.0 { 0.9 } else { 0.5 };
        let _ = write!(labels, "<text x=\"{x:.1}\" y=\"{y:.1}\" text-anchor=\"middle\" dominant-baseline=\"central\" fill=\"currentColor\" fill-opacity=\"{op}\">{text}</text>");
    }
    let defs = format!(
        "<defs><radialGradient id=\"pg-ball\" cx=\"36%\" cy=\"30%\" r=\"75%\"><stop offset=\"0\" stop-color=\"#ffffff\" stop-opacity=\"0.95\"/><stop offset=\"0.55\" stop-color=\"#cfe3f7\" stop-opacity=\"0.55\"/><stop offset=\"1\" stop-color=\"#8db3dd\" stop-opacity=\"0.5\"/></radialGradient><radialGradient id=\"pg-floor\"><stop offset=\"0\" stop-color=\"#000\" stop-opacity=\"0.14\"/><stop offset=\"1\" stop-color=\"#000\" stop-opacity=\"0\"/></radialGradient></defs><ellipse cx=\"{:.1}\" cy=\"{:.1}\" rx=\"{:.1}\" ry=\"{:.1}\" fill=\"url(#pg-floor)\"/>",
        v.cx,
        v.cy + 1.12 * v.r,
        0.8 * v.r,
        0.14 * v.r
    );
    (defs, back, front, labels)
}

fn svg_open(size: f64) -> String {
    format!("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {size} {size}\" width=\"{size}\" height=\"{size}\" font-family=\"-apple-system, Segoe UI, Helvetica, Arial, sans-serif\" font-size=\"{:.1}\">", size * 0.05)
}

fn ball(v: &View) -> String {
    format!("<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"url(#pg-ball)\" stroke=\"#3b5b82\" stroke-opacity=\"0.6\" stroke-width=\"1\"/>", v.cx, v.cy, v.r)
}

/// A Poincar\u{e9} sphere picture in the style of a Bloch sphere (translucent ball, light grid, bold equator, axes with
/// labels), with a cloud of points; each point is drawn with opacity in proportion to its weight and those behind
/// the ball are paler. `points` are (s1, s2, s3, weight). A size of 200 to 260 is enough to read it.
pub fn poincare_svg(points: &[(f64, f64, f64, f64)], size: f64) -> String {
    let v = View::new(size);
    let (defs, back, front, labels) = sphere_frame(&v);
    let wmax = points.iter().map(|p| p.3).filter(|w| w.is_finite()).fold(0.0, f64::max);
    let mut drawn: Vec<(f64, f64, f64, f64)> = points
        .iter()
        .filter(|p| wmax > 0.0 && p.3 >= 1e-4 * wmax)
        .map(|p| {
            let n = (p.0 * p.0 + p.1 * p.1 + p.2 * p.2).sqrt().max(1.0);
            let (x, y, d) = v.proj([p.0 / n, p.1 / n, p.2 / n]);
            (x, y, d, p.3 / wmax)
        })
        .collect();
    drawn.sort_by(|a, b| a.2.total_cmp(&b.2));
    let dot = |s: &mut String, x: f64, y: f64, w: f64, k: f64| {
        let _ = write!(s, "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"{:.1}\" fill=\"#e8590c\" fill-opacity=\"{:.2}\"/>", size * 0.01, (0.2 + 0.7 * w) * k);
    };
    let (mut back_pts, mut front_pts) = (String::new(), String::new());
    for (x, y, d, w) in drawn {
        if d >= 0.0 {
            dot(&mut front_pts, x, y, w, 1.0);
        } else {
            dot(&mut back_pts, x, y, w, 0.6);
        }
    }
    let mut s = svg_open(size);
    s.push_str(&defs);
    s.push_str(&back);
    s.push_str(&back_pts);
    s.push_str(&ball(&v));
    s.push_str(&front);
    s.push_str(&front_pts);
    s.push_str(&labels);
    s.push_str("</svg>");
    s
}

/// One mode on the sphere: a label, its (s1, s2, s3) and whether it belongs to a degenerate group (drawn dashed
/// because its polarization is then not unique).
#[derive(Clone, Debug)]
pub struct SphereState {
    pub label: String,
    pub s: [f64; 3],
    pub dashed: bool,
}

/// The same sphere with one numbered marker per mode on a line from the centre to its integrated state, like a Bloch vector: a fully
/// polarized mode ends on the surface, a hybrid mode ends inside it, and the length is the degree of polarization.
pub fn poincare_modes_svg(states: &[SphereState], size: f64) -> String {
    const PALETTE: [&str; 8] = ["#d9480f", "#1971c2", "#2f9e44", "#9c36b5", "#e67700", "#0c8599", "#c2255c", "#5f3dc4"];
    let v = View::new(size);
    let (defs, back, front, labels) = sphere_frame(&v);
    let (mut back_v, mut front_v) = (String::new(), String::new());
    let (cx, cy, _) = v.proj([0.0; 3]);
    for (i, st) in states.iter().enumerate() {
        if !st.s.iter().all(|c| c.is_finite()) {
            continue;
        }
        let n = (st.s[0] * st.s[0] + st.s[1] * st.s[1] + st.s[2] * st.s[2]).sqrt();
        let p = if n > 1.0 { st.s.map(|c| c / n) } else { st.s };
        let (x, y, d) = v.proj(p);
        let (dx, dy) = (x - cx, y - cy);
        let len = (dx * dx + dy * dy).sqrt();
        let col = PALETTE[i % PALETTE.len()];
        let out = if d >= 0.0 { &mut front_v } else { &mut back_v };
        let k = if d >= 0.0 { 1.0 } else { 0.6 };
        let dash = if st.dashed { " stroke-dasharray=\"4 3\"" } else { "" };
        if len >= 1.0 {
            let _ = write!(out, "<line x1=\"{cx:.1}\" y1=\"{cy:.1}\" x2=\"{x:.1}\" y2=\"{y:.1}\" stroke=\"{col}\" stroke-width=\"1.8\" stroke-opacity=\"{k}\" stroke-linecap=\"round\"{dash}/>");
        }
        // numbered badge at the tip
        let _ = write!(out, "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"{:.1}\" fill=\"{col}\" fill-opacity=\"{k}\" stroke=\"#ffffff\" stroke-width=\"0.8\"/>", size * 0.03);
        let _ = write!(out, "<text x=\"{x:.1}\" y=\"{y:.1}\" text-anchor=\"middle\" dominant-baseline=\"central\" font-weight=\"700\" font-size=\"{:.1}\" fill=\"#ffffff\">{}</text>", size * 0.036, st.label.replace('&', "&amp;").replace('<', "&lt;"));
    }
    let mut s = svg_open(size);
    s.push_str(&defs);
    s.push_str(&back);
    s.push_str(&back_v);
    s.push_str(&ball(&v));
    s.push_str(&front);
    s.push_str(&front_v);
    s.push_str(&labels);
    s.push_str("</svg>");
    s
}

// ------------------------------------------------------------------ the sphere as primitives

/// A polyline of the sphere picture in picture coordinates (x to the right, y down, the picture is a square of side `size`).
#[derive(Clone, Debug)]
pub struct WireLine {
    pub pts: Vec<(f64, f64)>,
    /// on the side facing the viewer (the far side is drawn dashed and paler)
    pub front: bool,
    pub width: f64,
    pub opacity: f64,
}

/// A label at a point of the picture.
#[derive(Clone, Debug)]
pub struct WireLabel {
    pub text: String,
    pub x: f64,
    pub y: f64,
    pub front: bool,
}

/// The empty Poincare sphere as lines and labels, for programs that draw it themselves (the SVG functions draw the same picture).
#[derive(Clone, Debug)]
pub struct SphereDrawing {
    pub size: f64,
    pub centre: (f64, f64),
    pub radius: f64,
    pub lines: Vec<WireLine>,
    pub labels: Vec<WireLabel>,
}

/// One mode on the picture: its tip, the colour (index into `MODE_COLOURS`) and whether it ends on the surface (`len` = 1) or inside.
#[derive(Clone, Debug)]
pub struct ModeMarker {
    pub label: String,
    pub centre: (f64, f64),
    pub tip: (f64, f64),
    pub front: bool,
    pub dashed: bool,
    pub colour: usize,
}

/// The colours of the mode markers, in the order of the numbers (as `#rrggbb`).
pub const MODE_COLOURS: [&str; 8] = ["#d9480f", "#1971c2", "#2f9e44", "#9c36b5", "#e67700", "#0c8599", "#c2255c", "#5f3dc4"];

/// Splits a curve into runs on the near and the far side of the sphere.
fn curve_runs(v: &View, pts: &[[f64; 3]]) -> Vec<(Vec<(f64, f64)>, bool)> {
    let mut out: Vec<(Vec<(f64, f64)>, bool)> = vec![];
    let mut run: Vec<(f64, f64)> = vec![];
    let mut is_front = true;
    for p in pts {
        let (x, y, d) = v.proj(*p);
        let f = d >= 0.0;
        if f != is_front {
            run.push((x, y));
            if run.len() > 1 {
                out.push((std::mem::take(&mut run), is_front));
            }
            run.clear();
            is_front = f;
        }
        run.push((x, y));
    }
    if run.len() > 1 {
        out.push((run, is_front));
    }
    out
}

/// The picture of the empty sphere in a square of side `size`.
pub fn sphere_drawing(size: f64) -> SphereDrawing {
    let v = View::new(size);
    let mut lines = vec![];
    for (curve, width, opacity) in frame_curves() {
        for (pts, front) in curve_runs(&v, &curve) {
            lines.push(WireLine { pts, front, width, opacity });
        }
    }
    let (xc, yc, _) = v.proj([0.0; 3]);
    for a in 0..3 {
        let mut dir = [0.0; 3];
        dir[a] = 1.0;
        for sign in [1.0, -1.0] {
            let (x, y, d) = v.proj(dir.map(|c| sign * 1.12 * c));
            lines.push(WireLine { pts: vec![(xc, yc), (x, y)], front: d >= 0.0, width: 0.9, opacity: if d >= 0.0 { 0.7 } else { 0.45 } });
        }
    }
    let mut labels = vec![];
    for (p, text) in [([1.0, 0.0, 0.0], "x"), ([-1.0, 0.0, 0.0], "y"), ([0.0, 1.0, 0.0], "+45\u{b0}"), ([0.0, -1.0, 0.0], "\u{2212}45\u{b0}"), ([0.0, 0.0, 1.0], "CCW"), ([0.0, 0.0, -1.0], "CW")] {
        let (x, y, d) = v.proj(p.map(|c| 1.27 * c));
        labels.push(WireLabel { text: text.to_string(), x, y, front: d >= 0.0 });
    }
    SphereDrawing { size, centre: (v.cx, v.cy), radius: v.r, lines, labels }
}

/// The markers of `states` on the picture of `sphere_drawing(size)`; states with a non-finite coordinate are skipped.
pub fn mode_markers(states: &[SphereState], size: f64) -> Vec<ModeMarker> {
    let v = View::new(size);
    let (cx, cy, _) = v.proj([0.0; 3]);
    states
        .iter()
        .enumerate()
        .filter(|(_, st)| st.s.iter().all(|c| c.is_finite()))
        .map(|(i, st)| {
            let n = (st.s[0] * st.s[0] + st.s[1] * st.s[1] + st.s[2] * st.s[2]).sqrt();
            let p = if n > 1.0 { st.s.map(|c| c / n) } else { st.s };
            let (x, y, d) = v.proj(p);
            ModeMarker { label: st.label.clone(), centre: (cx, cy), tip: (x, y), front: d >= 0.0, dashed: st.dashed, colour: i % MODE_COLOURS.len() }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const ONE: Cx = [1.0, 0.0];
    const ZERO: Cx = [0.0; 2];
    const I: Cx = [0.0, 1.0];

    fn smp(ex: Cx, ey: Cx, w: f64) -> Sample {
        Sample { ex, ey, ez: ZERO, w }
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn linear_x_and_y() {
        let s = Stokes::from_fields(ONE, ZERO, TimeConvention::PlusJ);
        assert_eq!((s.s0, s.s1, s.s2, s.s3), (1.0, 1.0, 0.0, 0.0));
        let r = analyze(&[smp(ONE, ZERO, 1.0)], &ClassOptions::default()).unwrap();
        assert_eq!(r.class, PolClass::QuasiTe);
        assert!(close(r.te_fraction, 1.0));
        let r = analyze(&[smp(ZERO, [0.0, 2.0], 1.0)], &ClassOptions::default()).unwrap();
        assert_eq!(r.class, PolClass::QuasiTm);
        assert!(close(r.te_fraction, 0.0));
    }

    #[test]
    fn diagonal_states() {
        let s = Stokes::from_fields(ONE, ONE, TimeConvention::PlusJ).normalized();
        assert!(close(s[1], 1.0) && close(s[0], 0.0) && close(s[2], 0.0));
        let r = analyze(&[smp(ONE, ONE, 1.0)], &ClassOptions::default()).unwrap();
        assert_eq!(r.class, PolClass::LinearPlus45);
        assert!(close(r.psi_deg, 45.0));
        let r = analyze(&[smp(ONE, [-1.0, 0.0], 1.0)], &ClassOptions::default()).unwrap();
        assert_eq!(r.class, PolClass::LinearMinus45);
        assert!(close(r.psi_deg, -45.0));
    }

    #[test]
    fn linear_at_any_angle() {
        for deg in [0.0f64, 10.0, 30.0, 60.0, 80.0] {
            let t = deg.to_radians();
            let s = Stokes::from_fields([t.cos(), 0.0], [t.sin(), 0.0], TimeConvention::MinusI);
            assert!(close(s.s1 / s.s0, (2.0 * t).cos()) && close(s.s2 / s.s0, (2.0 * t).sin()));
            assert!(close(s.psi_deg(), deg));
        }
    }

    #[test]
    fn circular_sense_depends_on_time_convention() {
        // Ey = i Ex: with exp(-iwt) the field is (cos wt, sin wt), counter-clockwise; with exp(+iwt) it is (cos wt, -sin wt)
        let m = Stokes::from_fields(ONE, I, TimeConvention::MinusI);
        let p = Stokes::from_fields(ONE, I, TimeConvention::PlusJ);
        assert!(close(m.s3 / m.s0, 1.0) && close(p.s3 / p.s0, -1.0));
        let r = analyze(&[smp(ONE, I, 1.0)], &ClassOptions { time: TimeConvention::MinusI, ..Default::default() }).unwrap();
        assert_eq!(r.class, PolClass::Circular { ccw: true });
        assert!(close(r.chi_deg, 45.0));
        let r = analyze(&[smp(ONE, I, 1.0)], &ClassOptions::default()).unwrap();
        assert_eq!(r.class, PolClass::Circular { ccw: false });
        assert!(close(r.chi_deg, -45.0));
    }

    /// S3 must agree with the rotation seen in the time domain, for random fields and both conventions
    #[test]
    fn s3_sign_matches_time_domain_rotation() {
        let mut seed = 12345u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        };
        for _ in 0..50 {
            let (ex, ey) = ([rnd(), rnd()], [rnd(), rnd()]);
            for (conv, sgn) in [(TimeConvention::PlusJ, 1.0), (TimeConvention::MinusI, -1.0)] {
                // field(t) = Re(E exp(+i sgn w t)), w t = phase
                let n = 2000;
                let pos = |k: usize| {
                    let ph = sgn * std::f64::consts::TAU * k as f64 / n as f64;
                    let f = |e: Cx| e[0] * ph.cos() - e[1] * ph.sin();
                    (f(ex), f(ey))
                };
                // twice the swept area; positive for counter-clockwise rotation in time
                let mut area = 0.0;
                for k in 0..n {
                    let (a, b) = (pos(k), pos(k + 1));
                    area += a.0 * b.1 - a.1 * b.0;
                }
                // time runs with k only when sgn*phase is increasing in time: PlusJ has phase = +w t, MinusI has phase = -w t,
                // and pos() already includes sgn, so k is time in both cases
                let s3 = Stokes::from_fields(ex, ey, conv).s3;
                if s3.abs() > 1e-6 {
                    assert_eq!(area > 0.0, s3 > 0.0, "ex {ex:?} ey {ey:?} {conv:?}");
                }
            }
        }
    }

    #[test]
    fn global_phase_does_not_change_the_result() {
        let ex = [0.7, -0.2];
        let ey = [-0.3, 0.5];
        let ph = [0.3f64.cos(), 0.3f64.sin()];
        let rot = |a: Cx| [a[0] * ph[0] - a[1] * ph[1], a[0] * ph[1] + a[1] * ph[0]];
        let a = Stokes::from_fields(ex, ey, TimeConvention::PlusJ);
        let b = Stokes::from_fields(rot(ex), rot(ey), TimeConvention::PlusJ);
        assert!(close(a.s0, b.s0) && close(a.s1, b.s1) && close(a.s2, b.s2) && close(a.s3, b.s3));
    }

    #[test]
    fn a_single_point_is_fully_polarized() {
        for &(ex, ey) in &[([0.3, 0.9], [-0.5, 0.1]), ([1.0, 0.0], [0.0, 0.0]), ([0.2, -0.4], [0.7, 0.7])] {
            let s = Stokes::from_fields(ex, ey, TimeConvention::MinusI);
            assert!(close(s.dop(), 1.0));
            let n = s.normalized();
            assert!(close((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt(), 1.0));
        }
    }

    #[test]
    fn mixture_of_x_and_y_is_hybrid_and_unpolarized() {
        let r = analyze(&[smp(ONE, ZERO, 1.0), smp(ZERO, ONE, 1.0)], &ClassOptions::default()).unwrap();
        assert!(close(r.dop, 0.0));
        assert_eq!(r.class, PolClass::Hybrid);
        assert!(close(r.te_fraction, 0.5));
    }

    #[test]
    fn weights_scale_the_contributions_and_te_fraction_matches_s1() {
        let r = analyze(&[smp(ONE, ZERO, 3.0), smp(ZERO, ONE, 1.0)], &ClassOptions { p_min: 0.4, ..Default::default() }).unwrap();
        assert!(close(r.sphere[0], 0.5));
        assert!(close(r.te_fraction, 0.75));
        assert!(close(r.te_fraction, 0.5 * (1.0 + r.sphere[0])));
        assert_eq!(r.class, PolClass::QuasiTe);
    }

    #[test]
    fn elliptical_state() {
        // Ey = 0.5 i Ex: s3 = 2 * 0.5 / 1.25 = 0.8
        let r = analyze(&[smp(ONE, [0.0, 0.5], 1.0)], &ClassOptions { time: TimeConvention::MinusI, ..Default::default() }).unwrap();
        assert!(close(r.sphere[2], 0.8));
        assert_eq!(r.class, PolClass::Elliptical { ccw: true });
        assert!(r.chi_deg > 0.0 && r.chi_deg < 45.0);
    }

    #[test]
    fn longitudinal_fraction() {
        let s = Sample { ex: ONE, ey: ZERO, ez: [0.0, 1.0], w: 1.0 };
        let r = analyze(&[s], &ClassOptions::default()).unwrap();
        assert!(close(r.long_fraction, 0.5));
    }

    #[test]
    fn no_intensity_gives_none() {
        assert!(analyze(&[smp(ZERO, ZERO, 1.0)], &ClassOptions::default()).is_none());
        assert!(analyze(&[], &ClassOptions::default()).is_none());
        assert!(analyze(&[smp(ONE, ZERO, 0.0)], &ClassOptions::default()).is_none());
    }

    #[test]
    fn orthogonal_states_are_pi_apart() {
        assert!(close(sphere_distance([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]), std::f64::consts::PI));
        assert!(close(sphere_distance([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]), std::f64::consts::FRAC_PI_2));
        assert!(close(sphere_distance([0.3, 0.0, 0.0], [2.0, 0.0, 0.0]), 0.0));
    }

    #[test]
    fn degenerate_modes_are_grouped() {
        let neff = [(2.40, 0.0), (1.80, 0.0), (2.4000001, 0.0), (1.2, 0.0), (1.8000002, 1e-7), (2.4000002, 0.0)];
        let g = degenerate_groups(&neff, 1e-5);
        assert_eq!(g, vec![vec![1, 4], vec![0, 2, 5]]);
        assert!(degenerate_groups(&neff, 0.0).is_empty());
    }

    #[test]
    fn vertex_weights_sum_to_the_area_and_skip_labels() {
        let pts = [0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 0.0, 1.0];
        let cells = [0u32, 1, 2, 0, 2, 3];
        let w = p1_vertex_weights(&pts, &cells, &[1, 2], &[]);
        assert!(close(w.iter().sum::<f64>(), 2.0));
        let w = p1_vertex_weights(&pts, &cells, &[1, 2], &[2]);
        assert!(close(w.iter().sum::<f64>(), 1.0));
        assert!(close(w[3], 0.0));
    }

    #[test]
    fn svg_is_well_formed_and_draws_the_points() {
        let pts = [(1.0, 0.0, 0.0, 1.0), (0.0, 0.0, -1.0, 0.5), (0.0, 1.0, 0.0, 0.0)];
        let svg = poincare_svg(&pts, 400.0);
        assert!(svg.starts_with("<svg") && svg.ends_with("</svg>"));
        // two points are drawn; the weightless one is left out
        assert_eq!(svg.matches("fill=\"#e8590c\"").count(), 2);
        assert!(!svg.contains("NaN"));
        assert!(poincare_svg(&[], 400.0).ends_with("</svg>"));
    }

    #[test]
    fn the_drawing_primitives_match_the_svg_picture() {
        let d = sphere_drawing(300.0);
        assert!((d.centre.0 - 150.0).abs() < 1e-9 && d.radius > 100.0);
        // every polyline point lies within the ball or on its outline, the equator is among the bold lines
        assert!(d.lines.iter().all(|l| l.pts.iter().all(|p| ((p.0 - d.centre.0).powi(2) + (p.1 - d.centre.1).powi(2)).sqrt() <= 1.13 * d.radius)));
        assert!(d.lines.iter().any(|l| l.width > 1.3 && l.front));
        assert_eq!(d.labels.len(), 6);
        let m = mode_markers(&[SphereState { label: "1".into(), s: [1.0, 0.0, 0.0], dashed: false }, SphereState { label: "2".into(), s: [f64::NAN, 0.0, 0.0], dashed: true }], 300.0);
        assert_eq!(m.len(), 1);
        // the x pole is drawn to the left of the centre (view from 35 degrees azimuth), facing the viewer, at radius * sin(35 degrees) from it
        assert!(m[0].front && (m[0].centre.0 - m[0].tip.0 - d.radius * 35f64.to_radians().sin()).abs() < 1e-9, "{:?}", m[0]);
        let svg = poincare_modes_svg(&[SphereState { label: "1".into(), s: [1.0, 0.0, 0.0], dashed: false }], 300.0);
        assert!(svg.contains(&format!("cx=\"{:.1}\" cy=\"{:.1}\"", m[0].tip.0, m[0].tip.1)));
    }
}
