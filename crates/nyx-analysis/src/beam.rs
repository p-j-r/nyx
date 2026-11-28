// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Analytic beams and the vectorial beam overlap of 3D frequency-domain electromagnetic results (FR-36, ADR 138).
//!
//! * beams ([`Beam`]): paraxial Gaussian (complex beam parameter, Gouy phase), thin-lens focus (scalar Debye integral
//!   over the aperture cone `sin a = NA / n`) and weakly guiding step-index fiber modes `LP_lm` ([`solve_lp_mode`]:
//!   the scalar dispersion relation `u J_{l+1}(u) K_l(w) = w K_{l+1}(w) J_l(u)` with `u^2 + w^2 = V^2`, solved with
//!   pure-Rust Bessel functions [`bessel_j`], [`bessel_k_scaled`]). Conventions as the far-field module: `exp(+j w t)`,
//!   propagation `exp(-j beta s)` along the beam direction, `H = d x E / eta` for every beam.
//! * overlap ([`beam_overlap`]): on a plane of located points, the four-field Poynting coupling efficiency
//!   `eta = |A B| / |Re P_f Re P_b|` with `A = int (E x Hb*) . n`, `B = int (Eb x H*) . n`, `P_f = int (E x H*) . n`,
//!   `P_b = int (Eb x Hb*) . n`, and the E-only scalar version `|int E . Eb*|^2 / (int |E|^2 int |Eb|^2)` (used when the
//!   result carries no H).

#![allow(clippy::needless_range_loop)]

use crate::farfield::{
    axis_param, cconj, ccross, cnorm2, cross3, dot3, f64_param, gauss_legendre, medium_from_params, opt_f64, par_map,
    rcross, scale3, sub3, unit3, usize_param, uv_axes, vec2_param, vec3_param, CVec3, CellSampler, FieldSelection,
    Locator, Medium, SurfaceSample, AXIS_NAMES, C0, C64, EPS0, MU0,
};
use nyx_physics::{Error, Mesh, SolveResult};
use serde_json::{json, Value};
use std::f64::consts::PI;

// ------------------------------------------------------------------------------------------------ Bessel functions

fn bessel_j_series(n: usize, x: f64) -> f64 {
    let h = 0.5 * x;
    let mut t = 1.0;
    for i in 1..=n {
        t *= h / i as f64;
    }
    let mut s = t;
    let q = -h * h;
    for k in 1..500 {
        t *= q / (k as f64 * (k + n) as f64);
        s += t;
        if t.abs() <= 1e-17 * s.abs() {
            break;
        }
    }
    s
}

fn bessel_j_trapezoid(n: usize, x: f64) -> f64 {
    // J_n(x) = 1/(2 pi) int_0^{2 pi} cos(n t - x sin t) dt: the trapezoid rule on a periodic analytic integrand converges
    // exponentially once the node count exceeds |x| + n by a margin
    let m = 2 * (x.ceil() as usize + n) + 64;
    let h = 2.0 * PI / m as f64;
    let mut s = 0.0;
    for i in 0..m {
        let t = i as f64 * h;
        s += (n as f64 * t - x * t.sin()).cos();
    }
    s / m as f64
}

/// Bessel function of the first kind `J_n(x)` for integer `n` and real `x` (power series where it does not cancel,
/// exponentially convergent trapezoid rule on the Bessel integral otherwise; absolute accuracy ~1e-15).
pub fn bessel_j(n: i32, x: f64) -> f64 {
    if !x.is_finite() {
        return f64::NAN;
    }
    let na = n.unsigned_abs() as usize;
    let mut sign = if n < 0 && na % 2 == 1 { -1.0 } else { 1.0 };
    if x < 0.0 && na % 2 == 1 {
        sign = -sign;
    }
    let ax = x.abs();
    if ax == 0.0 {
        return if na == 0 { 1.0 } else { 0.0 };
    }
    let v = if ax * ax / 4.0 < (na + 1) as f64 || ax < na as f64 {
        bessel_j_series(na, ax)
    } else {
        bessel_j_trapezoid(na, ax)
    };
    sign * v
}

/// Exponentially scaled modified Bessel function of the second kind `exp(x) K_n(x)`, `x > 0`, from
/// `K_n(x) = int_0^inf exp(-x cosh t) cosh(n t) dt` by the (exponentially convergent) trapezoid rule.
pub fn bessel_k_scaled(n: i32, x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    if x.is_infinite() {
        return 0.0;
    }
    let nn = n.unsigned_abs() as f64;
    let h = (0.5 / x.sqrt()).min(0.1);
    let f = |t: f64| {
        let s = (0.5 * t).sinh();
        (-2.0 * x * s * s + nn * t).exp() * 0.5 * (1.0 + (-2.0 * nn * t).exp())
    };
    let mut s = 0.5 * f(0.0);
    let mut k = 1usize;
    loop {
        let t = k as f64 * h;
        let v = f(t);
        s += v;
        if (v <= 1e-18 * s && t > 1.0) || k > 1_000_000 || !v.is_finite() {
            break;
        }
        k += 1;
    }
    s * h
}

/// Modified Bessel function of the second kind `K_n(x)`, `x > 0`.
pub fn bessel_k(n: i32, x: f64) -> f64 {
    (-x).exp() * bessel_k_scaled(n, x)
}

// ------------------------------------------------------------------------------------------------ LP fiber modes

/// A guided scalar mode `LP_lm` of a weakly guiding step-index fiber.
#[derive(Clone, Debug, PartialEq)]
pub struct LpMode {
    pub l: u32,
    pub m: u32,
    pub core_radius_m: f64,
    pub n_core: f64,
    pub n_clad: f64,
    pub wavelength_m: f64,
    /// normalized frequency `V = k0 a NA`
    pub v: f64,
    /// core and cladding parameters `u = a sqrt(k0^2 n_core^2 - beta^2)`, `w = a sqrt(beta^2 - k0^2 n_clad^2)`
    pub u: f64,
    pub w: f64,
    /// normalized propagation constant `b = 1 - u^2 / V^2`
    pub b: f64,
    pub n_eff: f64,
}

/// Solves the scalar LP dispersion relation for the `m`-th (1-based, by increasing `u`) root of azimuthal order `l`.
pub fn solve_lp_mode(
    l: u32,
    m: u32,
    core_radius_m: f64,
    n_core: f64,
    n_clad: f64,
    wavelength_m: f64,
) -> Result<LpMode, Error> {
    if !(core_radius_m > 0.0 && wavelength_m > 0.0 && n_clad > 0.0 && n_core > n_clad && m >= 1) {
        return Err(Error::invalid(format!(
            "an LP mode needs core_radius_m > 0, wavelength_m > 0, n_core > n_clad > 0 and m >= 1 (got a = {core_radius_m}, \
             lambda = {wavelength_m}, n_core = {n_core}, n_clad = {n_clad}, m = {m})"
        )));
    }
    let v = 2.0 * PI / wavelength_m * core_radius_m * (n_core * n_core - n_clad * n_clad).sqrt();
    let li = l as i32;
    let g = |u: f64| {
        let w = (v * v - u * u).max(0.0).sqrt();
        u * bessel_j(li + 1, u) * bessel_k_scaled(li, w) - w * bessel_k_scaled(li + 1, w) * bessel_j(li, u)
    };
    let (u0, u1) = (v * 1e-6, v * (1.0 - 1e-9));
    let nscan = ((400.0 * (v / PI + 1.0)).ceil() as usize).max(2000);
    let mut pu = u0;
    let mut pg = g(pu);
    let mut found = 0;
    for i in 1..=nscan {
        let u = u0 + (u1 - u0) * i as f64 / nscan as f64;
        let gv = g(u);
        if pg.is_finite() && gv.is_finite() && pg != 0.0 && (pg < 0.0) != (gv < 0.0) {
            found += 1;
            if found == m {
                let (mut a, mut b, mut ga) = (pu, u, pg);
                for _ in 0..200 {
                    let c = 0.5 * (a + b);
                    let gc = g(c);
                    if gc == 0.0 {
                        a = c;
                        b = c;
                        break;
                    }
                    if (gc < 0.0) == (ga < 0.0) {
                        a = c;
                        ga = gc;
                    } else {
                        b = c;
                    }
                    if b - a <= 1e-15 * v {
                        break;
                    }
                }
                let u = 0.5 * (a + b);
                let w = (v * v - u * u).max(0.0).sqrt();
                let bn = 1.0 - u * u / (v * v);
                let n_eff = (n_clad * n_clad + bn * (n_core * n_core - n_clad * n_clad)).sqrt();
                return Ok(LpMode { l, m, core_radius_m, n_core, n_clad, wavelength_m, v, u, w, b: bn, n_eff });
            }
        }
        pu = u;
        pg = gv;
    }
    Err(Error::invalid(format!("LP{l}{m} is not guided at V = {v:.6} (only {found} LP{l}x modes are)")))
}

impl LpMode {
    /// The scalar transverse profile at radius `r` and azimuth `phi` (1 at the core edge for `l = 0`, continuous with a
    /// continuous derivative at `r = a`).
    pub fn profile(&self, r: f64, phi: f64, sin_orientation: bool) -> f64 {
        let l = self.l as i32;
        let ang = if self.l == 0 {
            1.0
        } else if sin_orientation {
            (self.l as f64 * phi).sin()
        } else {
            (self.l as f64 * phi).cos()
        };
        let a = self.core_radius_m;
        let rad = if r < a {
            bessel_j(l, self.u * r / a) / bessel_j(l, self.u)
        } else {
            let x = self.w * r / a;
            (-(x - self.w)).exp() * bessel_k_scaled(l, x) / bessel_k_scaled(l, self.w)
        };
        rad * ang
    }

    /// Propagation constant `beta = k0 n_eff` (1/m).
    pub fn beta(&self) -> f64 {
        2.0 * PI / self.wavelength_m * self.n_eff
    }
}

// ------------------------------------------------------------------------------------------------ beams

/// The transverse shape of an analytic beam.
#[derive(Clone, Debug, PartialEq)]
pub enum BeamShape {
    /// paraxial Gaussian beam of waist radius (1/e field) `waist_m` at the origin
    Gaussian { waist_m: f64 },
    /// ideal thin-lens focus of numerical aperture `na` at the origin (uniform pupil, scalar Debye integral)
    ThinLens { na: f64 },
    /// LP fiber mode centred on the origin; `sin_orientation` picks `sin(l phi)` instead of `cos(l phi)`
    Lp { mode: LpMode, sin_orientation: bool },
}

/// An analytic beam: its shape, origin (focus or fiber centre), unit propagation direction and polarization (the
/// transverse unit vector; also the `phi = 0` reference of LP modes).
#[derive(Clone, Debug, PartialEq)]
pub struct Beam {
    pub shape: BeamShape,
    pub origin: [f64; 3],
    pub direction: [f64; 3],
    pub polarization: [f64; 3],
    e2: [f64; 3],
    gl: Vec<(f64, f64)>,
}

impl Beam {
    /// A beam; the polarization is made transverse to the direction.
    pub fn new(shape: BeamShape, origin: [f64; 3], direction: [f64; 3], polarization: [f64; 3]) -> Result<Beam, Error> {
        let d = unit3(direction).ok_or_else(|| Error::invalid("beam direction must be a non-zero vector"))?;
        let p = unit3(sub3(polarization, scale3(d, dot3(polarization, d))))
            .ok_or_else(|| Error::invalid("beam polarization must not be parallel to its direction"))?;
        match &shape {
            BeamShape::Gaussian { waist_m } if !waist_m.is_finite() || *waist_m <= 0.0 => {
                return Err(Error::invalid("gaussian beam needs waist_m > 0"));
            }
            BeamShape::ThinLens { na } if !na.is_finite() || *na <= 0.0 => {
                return Err(Error::invalid("thin_lens beam needs na > 0"));
            }
            _ => {}
        }
        let gl = match shape {
            BeamShape::ThinLens { .. } => {
                let (x, w) = gauss_legendre(64);
                x.iter().zip(&w).map(|(x, w)| (0.5 * (x + 1.0), 0.5 * w)).collect()
            }
            _ => Vec::new(),
        };
        Ok(Beam { shape, origin, direction: d, polarization: p, e2: cross3(d, p), gl })
    }

    /// `(E, H)` of the beam at `x` in `medium` (unit amplitude scale).
    pub fn fields(&self, medium: &Medium, x: [f64; 3]) -> (CVec3, CVec3) {
        let d = self.direction;
        let rel = sub3(x, self.origin);
        let s = dot3(rel, d);
        let t = sub3(rel, scale3(d, s));
        let rho2 = dot3(t, t);
        let (psi, eta) = match &self.shape {
            BeamShape::Gaussian { waist_m } => {
                let k = medium.k();
                let zr = 0.5 * k * waist_m * waist_m;
                let q = C64::new(s, zr);
                let amp = C64::new(0.0, zr) / q;
                let z = C64::new(0.0, -0.5 * k * rho2) / q;
                (amp * C64::cis(-k * s) * C64::cis(z.im) * z.re.exp(), medium.eta())
            }
            BeamShape::ThinLens { na } => {
                let k = medium.k();
                let alpha = (na / medium.index()).min(1.0).asin();
                let rho = rho2.sqrt();
                let mut acc = C64::ZERO;
                for &(xi, wi) in &self.gl {
                    let th = alpha * xi;
                    let (st, ct) = th.sin_cos();
                    acc += C64::cis(-k * s * ct) * (wi * alpha * st * bessel_j(0, k * rho * st));
                }
                (acc, medium.eta())
            }
            BeamShape::Lp { mode, sin_orientation } => {
                let r = rho2.sqrt();
                let phi = dot3(t, self.e2).atan2(dot3(t, self.polarization));
                let eta0 = (MU0 / EPS0).sqrt();
                (C64::cis(-mode.beta() * s) * mode.profile(r, phi, *sin_orientation), eta0 / mode.n_eff)
            }
        };
        let e = self.polarization.map(|p| psi * p);
        let h = rcross(d, &e).map(|z| z * (1.0 / eta));
        (e, h)
    }
}

/// The beam of a JSON object (`kind`: gaussian | thin_lens | lp); `default_origin`, `default_dir`, `default_pol` apply
/// when `focus`/`center`, `direction`, `polarization` are absent.
pub fn beam_from_params(
    v: &Value,
    medium: &Medium,
    default_origin: [f64; 3],
    default_dir: [f64; 3],
    default_pol: [f64; 3],
) -> Result<Beam, Error> {
    let kind = v.get("kind").and_then(Value::as_str).unwrap_or("gaussian");
    let origin = vec3_param(v, "focus")?.or(vec3_param(v, "center")?).unwrap_or(default_origin);
    let dir = vec3_param(v, "direction")?.unwrap_or(default_dir);
    let pol = vec3_param(v, "polarization")?.unwrap_or(default_pol);
    let shape = match kind {
        "gaussian" => BeamShape::Gaussian {
            waist_m: opt_f64(v, "waist_m")?.ok_or_else(|| Error::invalid("gaussian beam needs waist_m"))?,
        },
        "thin_lens" => {
            let na = opt_f64(v, "na")?.ok_or_else(|| Error::invalid("thin_lens beam needs na"))?;
            if na >= medium.index() {
                return Err(Error::invalid(format!(
                    "thin_lens na {na} must be below the medium index {}",
                    medium.index()
                )));
            }
            BeamShape::ThinLens { na }
        }
        "lp" | "fiber" => {
            let lambda = f64_param(v, "wavelength_m", C0 / medium.freq_hz)?;
            let mode = solve_lp_mode(
                usize_param(v, "l", 0)? as u32,
                usize_param(v, "m", 1)? as u32,
                opt_f64(v, "core_radius_m")?.ok_or_else(|| Error::invalid("lp beam needs core_radius_m"))?,
                opt_f64(v, "n_core")?.ok_or_else(|| Error::invalid("lp beam needs n_core"))?,
                opt_f64(v, "n_clad")?.ok_or_else(|| Error::invalid("lp beam needs n_clad"))?,
                lambda,
            )?;
            let sin_orientation = match v.get("orientation").and_then(Value::as_str) {
                None | Some("cos") => false,
                Some("sin") => true,
                Some(o) => return Err(Error::invalid(format!("lp orientation '{o}' is not cos or sin"))),
            };
            BeamShape::Lp { mode, sin_orientation }
        }
        _ => return Err(Error::invalid(format!("beam kind '{kind}' is not one of gaussian, thin_lens, lp"))),
    };
    Beam::new(shape, origin, dir, pol)
}

// ------------------------------------------------------------------------------------------------ overlap

/// The coupling of the sampled field into one beam.
#[derive(Clone, Debug, PartialEq)]
pub struct OverlapReport {
    /// four-field Poynting coupling efficiency (the E-only one when `uses_h` is false)
    pub eta: f64,
    /// E-only scalar overlap `|int E . Eb*|^2 / (int |E|^2 int |Eb|^2)`
    pub eta_e: f64,
    pub uses_h: bool,
    /// `int (E x Hb*) . n` and `int (Eb x H*) . n`
    pub a: C64,
    pub b: C64,
    /// `1/2 Re int (E x H*) . n` of the field and of the beam (W)
    pub p_field_w: f64,
    pub p_beam_w: f64,
}

/// The overlap of the sampled fields with `beam` (`use_h`: the four-field Poynting form; E-only otherwise).
pub fn beam_overlap(
    samples: &[SurfaceSample],
    beam: &Beam,
    medium: &Medium,
    use_h: bool,
    threads: usize,
) -> OverlapReport {
    let bf = par_map(samples.len(), threads, |i| beam.fields(medium, samples[i].x));
    let dotn = |v: CVec3, n: [f64; 3]| v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
    let (mut a, mut b, mut p1, mut p2, mut ee) = (C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO);
    let (mut e1, mut e2) = (0.0, 0.0);
    for (s, (eb, hb)) in samples.iter().zip(&bf) {
        a += dotn(ccross(&s.e, &cconj(hb)), s.n) * s.w;
        b += dotn(ccross(eb, &cconj(&s.h)), s.n) * s.w;
        p1 += dotn(ccross(&s.e, &cconj(&s.h)), s.n) * s.w;
        p2 += dotn(ccross(eb, &cconj(hb)), s.n) * s.w;
        for k in 0..3 {
            ee += s.e[k] * eb[k].conj() * s.w;
        }
        e1 += s.w * cnorm2(&s.e);
        e2 += s.w * cnorm2(eb);
    }
    let eta_e = if e1 > 0.0 && e2 > 0.0 { ee.norm_sqr() / (e1 * e2) } else { 0.0 };
    let den = (p1.re * p2.re).abs();
    let eta = if use_h {
        if den > 0.0 {
            (a * b).abs() / den
        } else {
            0.0
        }
    } else {
        eta_e
    };
    OverlapReport {
        eta,
        eta_e,
        uses_h: use_h,
        a,
        b,
        p_field_w: if use_h { 0.5 * p1.re } else { f64::NAN },
        p_beam_w: 0.5 * p2.re,
    }
}

/// Whether `result` carries 3D per-cell `E` fields (and so the beam overlap applies).
pub fn overlap_beam_applies(result: &SolveResult) -> bool {
    result.dim == 3 && result.solutions.iter().any(|s| s.cell_fields.iter().any(|f| f.name == "Ex" || f.name == "Esx"))
}

/// The beam-overlap analysis with JSON parameters (see [`BeamOverlapAnalysis`]).
pub fn overlap_beam_with_params(mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
    let beams_json: Vec<Value> =
        match (params.get("beams").filter(|v| !v.is_null()), params.get("beam").filter(|v| !v.is_null())) {
            (Some(Value::Array(a)), _) => a.clone(),
            (Some(_), _) => return Err(Error::invalid("params.beams must be an array of beam objects")),
            (None, Some(b)) => vec![b.clone()],
            (None, None) => {
                return Err(Error::invalid("params.beam (or beams) is needed: {kind: gaussian | thin_lens | lp, ...}"))
            }
        };
    let sel = FieldSelection::from_params(params)?;
    let sampler = CellSampler::new(mesh, result, &sel)?;
    let threads = usize_param(params, "threads", femcore::par::threads())?.max(1);
    let medium = medium_from_params(params, sampler.solution_value())?;
    let pl = params.get("plane").cloned().unwrap_or(Value::Null);
    let axis = axis_param(&pl, "axis", 2)?;
    let pos = f64_param(&pl, "position", 0.0)?;
    let sign = if f64_param(&pl, "normal_sign", 1.0)? < 0.0 { -1.0 } else { 1.0 };
    let (iu, iv) = uv_axes(axis);
    let (blo, bhi) = sampler.disc().bbox();
    let lo = vec2_param(&pl, "lo")?.unwrap_or([blo[iu], blo[iv]]);
    let hi = vec2_param(&pl, "hi")?.unwrap_or([bhi[iu], bhi[iv]]);
    if !(hi[0] > lo[0] && hi[1] > lo[1]) {
        return Err(Error::invalid("plane needs hi > lo on both in-plane axes"));
    }
    let nuv = match pl.get("n").and_then(Value::as_array) {
        Some(a) if a.len() == 2 => {
            [a[0].as_u64().unwrap_or(1).max(1) as usize, a[1].as_u64().unwrap_or(1).max(1) as usize]
        }
        _ => {
            let h = opt_f64(&pl, "h")?.unwrap_or_else(|| sampler.default_h()).max(1e-300);
            [((hi[0] - lo[0]) / h).ceil() as usize, ((hi[1] - lo[1]) / h).ceil() as usize]
        }
    };
    let loc = Locator::new(sampler.disc());
    let mut samples = Vec::new();
    let (area, ncell, outside) = sampler.plane_grid(&loc, axis, pos, lo, hi, nuv, sign, &mut samples);
    if samples.is_empty() {
        return Err(Error::invalid("the overlap plane has no point inside the mesh"));
    }
    let use_h = match params.get("mode").and_then(Value::as_str) {
        None | Some("auto") => sampler.has_h(),
        Some("poynting") => {
            if !sampler.has_h() {
                return Err(Error::invalid(format!("mode poynting needs the H fields {:?}", sel.h_names)));
            }
            true
        }
        Some("e_only") => false,
        Some(m) => return Err(Error::invalid(format!("mode '{m}' is not one of auto, poynting, e_only"))),
    };
    let mut centre = [0.0; 3];
    centre[axis] = pos;
    centre[iu] = 0.5 * (lo[0] + hi[0]);
    centre[iv] = 0.5 * (lo[1] + hi[1]);
    let mut ndir = [0.0; 3];
    ndir[axis] = sign;
    let mut pol = [0.0; 3];
    pol[iu] = 1.0;
    let mut reports = Vec::new();
    for bj in &beams_json {
        let beam = beam_from_params(bj, &medium, centre, ndir, pol)?;
        let r = beam_overlap(&samples, &beam, &medium, use_h, threads);
        let mut o = json!({
            "kind": bj.get("kind").and_then(Value::as_str).unwrap_or("gaussian"),
            "eta": r.eta, "eta_db": 10.0 * r.eta.log10(), "eta_e": r.eta_e, "uses_h": r.uses_h,
            "A_re": r.a.re, "A_im": r.a.im, "B_re": r.b.re, "B_im": r.b.im, "P_beam_w": r.p_beam_w,
            "origin_m": beam.origin, "direction": beam.direction, "polarization": beam.polarization
        });
        if let (BeamShape::Lp { mode, .. }, Some(m)) = (&beam.shape, o.as_object_mut()) {
            m.insert(
                "lp".into(),
                json!({ "l": mode.l, "m": mode.m, "V": mode.v, "u": mode.u, "w": mode.w, "b": mode.b, "n_eff": mode.n_eff,
                        "wavelength_m": mode.wavelength_m }),
            );
        }
        reports.push(o);
    }
    let p_field = if use_h { crate::farfield::surface_flux(&samples) } else { f64::NAN };
    let (en, hn) = sampler.field_names();
    Ok(json!({
        "analysis": "overlap_beam", "solution": sel.solution, "e_fields": en, "h_fields": hn,
        "freq_hz": medium.freq_hz, "eps_r": medium.eps_r, "mu_r": medium.mu_r, "wavelength_m": medium.wavelength(),
        "plane": { "axis": AXIS_NAMES[axis], "position_m": pos, "lo_m": lo, "hi_m": hi, "n": nuv,
                   "normal_sign": sign, "n_samples": samples.len(), "n_cells": ncell, "n_outside": outside, "area_m2": area },
        "P_field_w": p_field,
        "eta": reports.first().and_then(|r| r.get("eta")).cloned().unwrap_or(Value::Null),
        "beams": reports
    }))
}

/// The vectorial beam-overlap analysis behind the `Analysis` trait.
///
/// Parameters: `solution`, `field_prefix`, `e_fields`/`h_fields`, `freq_hz`, `eps_r`, `mu_r` (as `farfield`), `plane`
/// {axis (z), position (0), lo/hi (mesh extent), n: [nu, nv] or h, normal_sign (1)}, `beam` or `beams`
/// ({kind: gaussian, waist_m} | {kind: thin_lens, na} | {kind: lp, l, m, core_radius_m, n_core, n_clad, wavelength_m,
/// orientation}, each with optional focus/center, direction (plane normal), polarization (plane u axis)), `mode`
/// (auto | poynting | e_only), `threads`.
pub struct BeamOverlapAnalysis;

impl crate::Analysis for BeamOverlapAnalysis {
    fn id(&self) -> &'static str {
        "overlap_beam"
    }
    fn name(&self) -> &'static str {
        "Vectorial overlap with Gaussian, thin-lens and fiber-mode beams"
    }
    fn applies_to(&self, result: &SolveResult) -> bool {
        overlap_beam_applies(result)
    }
    fn run(&self, mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        overlap_beam_with_params(mesh, result, &Value::Null)
    }
    fn run_with_params(&self, mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        overlap_beam_with_params(mesh, result, params)
    }
}

/// Field samples of an analytic beam on a plane grid (for tests and for beam-to-beam overlaps).
pub fn beam_plane_samples(
    beam: &Beam,
    medium: &Medium,
    axis: usize,
    position: f64,
    lo: [f64; 2],
    hi: [f64; 2],
    nuv: [usize; 2],
) -> Vec<SurfaceSample> {
    let mut n = [0.0; 3];
    n[axis.min(2)] = 1.0;
    crate::farfield::plane_grid_points(axis, position, lo, hi, nuv)
        .into_iter()
        .map(|(x, w)| {
            let (e, h) = beam.fields(medium, x);
            SurfaceSample { x, w, n, e, h }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bessel_values_match_tables() {
        // Abramowitz & Stegun table 9.1 / 9.8 values
        assert!((bessel_j(0, 1.0) - 0.765_197_686_557_966_6).abs() < 1e-14);
        assert!((bessel_j(1, 1.0) - 0.440_050_585_744_933_5).abs() < 1e-14);
        assert!((bessel_j(0, 10.0) - (-0.245_935_764_451_348_3)).abs() < 1e-14);
        assert!((bessel_j(2, 5.0) - 0.046_565_116_277_752_2).abs() < 1e-14);
        assert!(bessel_j(1, crate::farfield::AIRY_ZERO).abs() < 1e-12);
        assert!(bessel_j(0, 2.404_825_557_695_773).abs() < 1e-14);
        assert!((bessel_k(0, 1.0) - 0.421_024_438_240_708_3).abs() < 1e-13);
        assert!((bessel_k(1, 1.0) - 0.601_907_230_197_234_6).abs() < 1e-13);
        assert!((bessel_k(1, 0.1) - 9.853_844_780_870_606).abs() < 1e-11);
    }

    #[test]
    fn lp11_cutoff_is_the_first_zero_of_j0() {
        let a = 4.0e-6;
        let (nco, ncl) = (1.45f64, 1.44f64);
        let na = (nco * nco - ncl * ncl).sqrt();
        let lam_cut = 2.0 * PI * a * na / 2.404_825_557_695_773;
        assert!(solve_lp_mode(1, 1, a, nco, ncl, lam_cut * 1.01).is_err());
        assert!(solve_lp_mode(1, 1, a, nco, ncl, lam_cut * 0.99).is_ok());
        let m = solve_lp_mode(0, 1, a, nco, ncl, lam_cut * 1.2).unwrap();
        assert!(m.n_eff > ncl && m.n_eff < nco && m.b > 0.0 && m.b < 1.0);
    }
}
